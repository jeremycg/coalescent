#!/usr/bin/env python3
"""Enforce the production-core boundary used by standalone DSP tests."""

from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]

# The standalone test must compile the same SDK-free implementation included by
# the Rack wrapper. Tests may provide independent oracles and explicit historical
# mutants, but must not include a module .cpp or maintain another production path.
CONTRACTS = (
    ("GENDYN", "src/GENDYN.cpp", "tools/stability/gendyn.cpp", "src/dsp/gendyn_core.hpp"),
    ("Haptik", "src/Haptik.cpp", "tools/stability/haptik.cpp", "src/dsp/haptik_core.hpp"),
    ("Axon", "src/neuron/Axon.cpp", "tools/stability/axon.cpp", "src/dsp/neuron_models.hpp"),
    ("Soma", "src/neuron/Soma.cpp", "tools/stability/soma.cpp", "src/dsp/neuron_models.hpp"),
    ("Operon", "src/Operon.cpp", "tools/stability/operon.cpp", "src/dsp/operon_core.hpp"),
    ("Bunnies", "src/Bunnies.cpp", "tools/stability/bunnies.cpp", "src/dsp/bunnies_core.hpp"),
    ("Foxes", "src/Foxes.cpp", "tools/stability/foxes.cpp", "src/dsp/foxes_core.hpp"),
    ("Finches", "src/Finches.cpp", "tools/stability/finches.cpp", "src/dsp/finches_field.hpp"),
    ("Islands", "src/Islands.cpp", "tools/stability/islands.cpp", "src/dsp/islands_model.hpp"),
    ("Archipelago", "src/Archipelago.cpp", "tools/stability/archipelago.cpp", "src/dsp/archipelago_field.hpp"),
    ("Lineages", "src/Lineages.cpp", "tools/stability/lineages.cpp", "src/lineages/state.hpp"),
)

# Direct includes alone are too weak: an unused include could coexist with a
# private copy. Require representative production entry points in both sides of
# every mapping. These are deliberately source-level guardrails, not a claim of
# semantic proof; behavioral contracts still provide the proof.
REQUIRED_USAGES = {
    "GENDYN": (
        ("completeCycleAndRetarget", "quantizeDuration"),
        ("completeCycleAndRetarget", "quantizeDuration"),
    ),
    "Haptik": (
        ("captureInterpolatedFrame", "shouldStep"),
        ("captureInterpolatedFrame", "shouldStep"),
    ),
    "Axon": (
        ("Core::advanceObservation", "Core::repair"),
        ("Core::advanceObservation", "Core::repair"),
    ),
    "Soma": (
        ("Core::advanceObservation", "Core::repair"),
        ("Core::advanceObservation", "Core::repair"),
    ),
    "Operon": (
        ("advanceAcceptedSubstep", "HillLut"),
        ("advanceAcceptedSubstep", "HillLut"),
    ),
    "Bunnies": (
        ("coalescent::bunnies::step", "coalescent::bunnies::repairState"),
        ("model::step", "model::repairState"),
    ),
    "Foxes": (
        ("coalescent::foxes::step", "coalescent::foxes::repairState"),
        ("model::step", "model::repairState"),
    ),
    "Finches": (
        ("field.advance", "candidate.restore"),
        ("field.advance", "resumed.restore"),
    ),
    "Islands": (
        ("model.advance", "candidate.restore"),
        ("neutral.advance", "resumed.restore"),
    ),
    "Archipelago": (
        ("field.advance", "candidate.restore"),
        ("field.advance", "resumed.restore"),
    ),
    "Lineages": (
        ("generator.generate", "playback.advance"),
        ("generator.generate", "playback.advance"),
    ),
}

# Auxiliary verification and diagnostic tools are subject to the same boundary.
# They may add independent measurements or analytical oracles, but model
# transitions must still come from the production core.
AUXILIARY_CONTRACTS = (
    ("SIMD neurons", "tools/simd_equiv.cpp", "src/dsp/neuron_models.hpp",
     ("AxonCore::advanceObservation", "SomaCore::advanceObservation")),
    ("SIMD Haptik", "tools/simd_equiv.cpp", "src/dsp/haptik_core.hpp",
     ("coalescent::haptik::advanceStateValue",)),
    ("MIN_SUB profiler", "tools/perf_minsub.cpp", "src/dsp/neuron_models.hpp",
     ("Core::advanceObservation", "scalarScheduleWithMinimum")),
    ("Axon renderer", "tools/render_wav_axon.cpp", "src/dsp/neuron_models.hpp",
     ("Core::advanceObservation", "Core::repair")),
    ("Soma renderer", "tools/render_wav_soma.cpp", "src/dsp/neuron_models.hpp",
     ("Core::advanceObservation", "Core::repair")),
    ("RK4 analytic contract", "tools/integrator_equiv.cpp", "src/dsp/rk4.hpp",
     ("coalescent::rk4<",)),
)

KNOWN_COPY_PATTERNS = (
    (re.compile(r"\bstepDirect\b"), "test-only model step"),
    (re.compile(r"\brk4(?:A|S)_old\b"), "copied historical RK4 stepper"),
    (re.compile(r"\b(?:fFHN|fHR)\s*\("), "copied neuron derivative"),
)

INCLUDE_RE = re.compile(r'^\s*#\s*include\s+"([^"]+)"', re.MULTILINE)


def resolved_local_includes(source: Path):
    text = source.read_text(encoding="utf-8")
    for include in INCLUDE_RE.findall(text):
        relative = (source.parent / include).resolve()
        include_root = (ROOT / "src" / include).resolve()
        yield include, relative
        if include_root != relative:
            yield include, include_root


def includes(source: Path, target: Path) -> bool:
    wanted = target.resolve()
    return any(candidate == wanted for _, candidate in resolved_local_includes(source))


def fail(message: str) -> None:
    print(f"shared-core architecture: {message}", file=sys.stderr)


def main() -> int:
    failures = 0
    checked_cores = set()

    for name, module_name, test_name, core_name in CONTRACTS:
        module = ROOT / module_name
        test = ROOT / test_name
        core = ROOT / core_name

        for path in (module, test, core):
            if not path.is_file():
                fail(f"{name}: missing {path.relative_to(ROOT)}")
                failures += 1

        if not all(path.is_file() for path in (module, test, core)):
            continue

        if not includes(module, core):
            fail(f"{name}: Rack wrapper does not include {core_name}")
            failures += 1
        if not includes(test, core):
            fail(f"{name}: stability test does not include {core_name}")
            failures += 1

        module_text = module.read_text(encoding="utf-8")
        test_text = test.read_text(encoding="utf-8")
        if re.search(r'#\s*include\s+"[^"]+\.cpp"', test_text):
            fail(f"{name}: stability tests must not include a module .cpp")
            failures += 1

        module_usages, test_usages = REQUIRED_USAGES[name]
        for token in module_usages:
            if token not in module_text:
                fail(f"{name}: Rack wrapper does not call shared entry point {token}")
                failures += 1
        for token in test_usages:
            if token not in test_text:
                fail(f"{name}: stability test does not call shared entry point {token}")
                failures += 1

        checked_cores.add(core)

    mapped_tool_texts = {}
    for name, tool_name, core_name, required_usages in AUXILIARY_CONTRACTS:
        tool = ROOT / tool_name
        core = ROOT / core_name

        for path in (tool, core):
            if not path.is_file():
                fail(f"{name}: missing {path.relative_to(ROOT)}")
                failures += 1

        if not tool.is_file() or not core.is_file():
            continue

        if not includes(tool, core):
            fail(f"{name}: {tool_name} does not include {core_name}")
            failures += 1

        tool_text = tool.read_text(encoding="utf-8")
        mapped_tool_texts[tool] = tool_text
        if re.search(r'#\s*include\s+"[^"]+\.cpp"', tool_text):
            fail(f"{name}: auxiliary tools must not include a module .cpp")
            failures += 1

        for token in required_usages:
            if token not in tool_text:
                fail(f"{name}: {tool_name} does not call shared entry point {token}")
                failures += 1

        checked_cores.add(core)

    copied_sources = {ROOT / contract[2] for contract in CONTRACTS}
    copied_sources.update(mapped_tool_texts)
    for source in sorted(copied_sources):
        text = mapped_tool_texts.get(source)
        if text is None:
            text = source.read_text(encoding="utf-8")
        for pattern, description in KNOWN_COPY_PATTERNS:
            if pattern.search(text):
                fail(f"{source.relative_to(ROOT)} retains a {description}")
                failures += 1

    for core in sorted(checked_cores):
        text = core.read_text(encoding="utf-8")
        forbidden = ("plugin.hpp", "rack.hpp", "<rack/", '"rack/')
        for token in forbidden:
            if token in text:
                fail(f"{core.relative_to(ROOT)} is not SDK-free ({token})")
                failures += 1

    if failures:
        return 1

    print(
        "Shared-core architecture: all 11 Rack wrappers and stability suites "
        "plus 6 auxiliary contracts use SDK-free production cores"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
