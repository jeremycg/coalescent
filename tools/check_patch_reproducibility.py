#!/usr/bin/env python3
"""Regenerate every scripted demo in isolation and compare it with the repo."""

import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parent.parent
EXPECTED_PATCH_COUNT = 49


class ReproducibilityError(RuntimeError):
    pass


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def decompress_archive(path):
    result = subprocess.run(
        ["zstd", "-q", "-d", "-c", str(path)],
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        message = result.stderr.decode("utf-8", errors="replace").strip()
        raise ReproducibilityError(f"could not decompress {path.name}: {message}")
    return result.stdout


def run_generators(isolated_root, generators):
    environment = os.environ.copy()
    environment["COALESCENT_SKIP_PATCH_INSTALL"] = "1"
    environment["PYTHONDONTWRITEBYTECODE"] = "1"

    for source in generators:
        script = isolated_root / source.relative_to(ROOT)
        result = subprocess.run(
            [sys.executable, str(script)],
            cwd=isolated_root,
            env=environment,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            details = "\n".join(part for part in (result.stdout, result.stderr) if part.strip())
            raise ReproducibilityError(
                f"{source.name} failed with exit status {result.returncode}\n{details}"
            )


def require_same_names(expected, generated):
    expected_names = set(expected)
    generated_names = set(generated)
    if expected_names == generated_names:
        return

    missing = sorted(expected_names - generated_names)
    unexpected = sorted(generated_names - expected_names)
    details = []
    if missing:
        details.append("missing: " + ", ".join(missing))
    if unexpected:
        details.append("unexpected: " + ", ".join(unexpected))
    raise ReproducibilityError("generated patch set differs (" + "; ".join(details) + ")")


def main():
    generators = sorted((ROOT / "tools").glob("make_patch*.py"))
    expected = {path.name: path for path in sorted((ROOT / "patches").glob("*.vcv"))}
    if len(expected) != EXPECTED_PATCH_COUNT:
        raise ReproducibilityError(
            f"expected {EXPECTED_PATCH_COUNT} checked-in generated patches, found {len(expected)}"
        )
    if not generators:
        raise ReproducibilityError("no patch generators found")

    with tempfile.TemporaryDirectory(prefix="coalescent-patches-") as temporary:
        isolated_root = Path(temporary) / "coalescent"
        shutil.copytree(
            ROOT / "tools",
            isolated_root / "tools",
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
        )
        shutil.copy2(ROOT / "plugin.json", isolated_root / "plugin.json")
        generated_dir = isolated_root / "patches"
        generated_dir.mkdir()

        run_generators(isolated_root, generators)
        generated = {path.name: path for path in sorted(generated_dir.glob("*.vcv"))}
        require_same_names(expected, generated)

        first_pass = {name: path.read_bytes() for name, path in generated.items()}
        mismatches = []
        for name, expected_path in expected.items():
            expected_tar = decompress_archive(expected_path)
            generated_tar = decompress_archive(generated[name])
            if expected_tar != generated_tar:
                mismatches.append(
                    f"{name}: checked-in tar {sha256(expected_tar)}, "
                    f"generated tar {sha256(generated_tar)}"
                )
        if mismatches:
            raise ReproducibilityError(
                "canonical archive content differs:\n  " + "\n  ".join(mismatches)
            )

        # A second pass proves byte stability with the current zstd toolchain.
        # Cross-version comparison uses uncompressed tar bytes above because zstd
        # is free to change its compressed representation between releases.
        for path in generated_dir.glob("*.vcv"):
            path.unlink()
        run_generators(isolated_root, generators)
        regenerated = {
            path.name: path for path in sorted(generated_dir.glob("*.vcv"))
        }
        require_same_names(expected, regenerated)
        second_pass = {name: path.read_bytes() for name, path in regenerated.items()}
        unstable = [name for name in sorted(first_pass) if first_pass[name] != second_pass[name]]
        if unstable:
            raise ReproducibilityError(
                "same-toolchain archive bytes changed on the second pass: " + ", ".join(unstable)
            )

    print(
        f"Patch reproducibility: {len(expected)} archives match canonical checked-in "
        "content; repeated generation is byte-identical"
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ReproducibilityError as error:
        print(f"Patch reproducibility FAILED: {error}", file=sys.stderr)
        sys.exit(1)
