#!/usr/bin/env python3
"""Fast SDK-free checks for manifests, documentation, and Rack demo patches."""

import io
import json
import re
import subprocess
import sys
import tarfile
from pathlib import Path
from urllib.parse import unquote, urlsplit

from patch_utils import patch_json_bytes


ROOT = Path(__file__).resolve().parent.parent

# Widths are Rack grid units (HP). Port counts are (inputs, outputs). Keeping the
# third-party data here makes this check independent of an installed Rack library.
MODULES = {
    ("Core", "AudioInterface"): (10, 8, 8),
    ("Core", "MIDIToCVInterface"): (8, 0, 6),
    ("Fundamental", "8vert"): (8, 8, 8),
    ("Fundamental", "ADSR"): (9, 5, 1),
    ("Fundamental", "LFO"): (9, 5, 4),
    ("Fundamental", "Merge"): (5, 16, 1),
    ("Fundamental", "Mixer"): (3, 6, 1),
    ("Fundamental", "Sum"): (3, 1, 1),
    ("Fundamental", "VCA-1"): (3, 2, 1),
    ("Fundamental", "VCMixer"): (9, 9, 5),
    ("Fundamental", "VCO"): (9, 6, 4),
    ("Coalescent", "GENDYN"): (12, 4, 3),
    ("Coalescent", "Haptik"): (18, 7, 2),
    ("Coalescent", "Axon"): (12, 5, 3),
    ("Coalescent", "Soma"): (12, 5, 3),
    ("Coalescent", "Operon"): (14, 5, 6),
    ("Coalescent", "Bunnies"): (12, 4, 4),
    ("Coalescent", "Foxes"): (12, 4, 6),
    ("Coalescent", "Finches"): (14, 6, 7),
    ("Coalescent", "Islands"): (16, 8, 10),
    ("Coalescent", "Archipelago"): (18, 8, 7),
    ("Coalescent", "Lineages"): (16, 4, 6),
}

LINK_RE = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
HEADING_RE = re.compile(r"^#{1,6}\s+(.+?)\s*#*\s*$", re.MULTILINE)
MANUAL_RE = re.compile(
    r"^https://github\.com/jeremycg/coalescent/blob/main/(docs/[^?#]+\.md)$"
)


class ValidationError(Exception):
    pass


def require(condition, message):
    if not condition:
        raise ValidationError(message)


def github_anchor(text):
    text = re.sub(r"<[^>]+>", "", text).strip().lower()
    text = re.sub(r"[`*_~]", "", text)
    text = re.sub(r"[^\w\- ]", "", text, flags=re.UNICODE)
    return re.sub(r"\s+", "-", text)


def markdown_anchors(text):
    counts = {}
    anchors = set()
    for heading in HEADING_RE.findall(text):
        base = github_anchor(heading)
        index = counts.get(base, 0)
        counts[base] = index + 1
        anchors.add(base if index == 0 else f"{base}-{index}")
    return anchors


def validate_manifest():
    manifest = json.loads((ROOT / "plugin.json").read_text(encoding="utf-8"))
    modules = manifest.get("modules")
    require(isinstance(modules, list) and modules, "plugin.json has no modules")
    slugs = [module.get("slug") for module in modules]
    require(len(slugs) == len(set(slugs)), "plugin.json contains duplicate module slugs")

    plugin_cpp = (ROOT / "src/plugin.cpp").read_text(encoding="utf-8")
    plugin_hpp = (ROOT / "src/plugin.hpp").read_text(encoding="utf-8")
    for module in modules:
        slug = module.get("slug")
        require(isinstance(slug, str) and slug, "plugin.json module has no slug")
        match = MANUAL_RE.match(module.get("manualUrl", ""))
        require(match is not None, f"{slug}: manualUrl is not a repository docs URL")
        manual = ROOT / match.group(1)
        require(manual.is_file(), f"{slug}: missing manual {manual.relative_to(ROOT)}")
        require(f"p->addModel(model{slug});" in plugin_cpp, f"{slug}: not registered in plugin.cpp")
        require(f"extern Model* model{slug};" in plugin_hpp, f"{slug}: not declared in plugin.hpp")
    return len(modules)


def validate_markdown():
    files = [ROOT / "README.md", *sorted((ROOT / "docs").glob("*.md"))]
    cache = {path.resolve(): path.read_text(encoding="utf-8") for path in files}
    checked = 0
    for source in files:
        text = cache[source.resolve()]
        for raw_target in LINK_RE.findall(text):
            target = raw_target.strip()
            if target.startswith("<") and target.endswith(">"):
                target = target[1:-1]
            parts = urlsplit(target)
            if parts.scheme or target.startswith("//"):
                continue
            path_part = unquote(parts.path)
            destination = (source.parent / path_part).resolve() if path_part else source.resolve()
            require(destination.exists(), f"{source.relative_to(ROOT)}: broken link {raw_target}")
            require(
                destination == ROOT.resolve() or ROOT.resolve() in destination.parents,
                f"{source.relative_to(ROOT)}: link escapes repository: {raw_target}",
            )
            if parts.fragment:
                require(destination.is_file(), f"{source.relative_to(ROOT)}: anchor target is not a file")
                destination_text = cache.get(destination)
                if destination_text is None:
                    destination_text = destination.read_text(encoding="utf-8")
                    cache[destination] = destination_text
                anchor = unquote(parts.fragment).lower()
                require(
                    anchor in markdown_anchors(destination_text),
                    f"{source.relative_to(ROOT)}: missing anchor #{parts.fragment} in {destination.relative_to(ROOT)}",
                )
            checked += 1
    return len(files), checked


def read_patch(path, generated):
    result = subprocess.run(
        ["zstd", "-q", "-d", "-c", str(path)], capture_output=True, check=False
    )
    require(result.returncode == 0, f"{path.relative_to(ROOT)}: invalid zstd archive")
    try:
        with tarfile.open(fileobj=io.BytesIO(result.stdout), mode="r:") as archive:
            members = archive.getmembers()
            require(members, f"{path.relative_to(ROOT)}: empty archive")
            require(
                all(not member.name.startswith("/") and ".." not in Path(member.name).parts for member in members),
                f"{path.relative_to(ROOT)}: unsafe archive member path",
            )
            patch_members = [
                member for member in members
                if member.isfile() and member.name.removeprefix("./") == "patch.json"
            ]
            require(len(patch_members) == 1, f"{path.relative_to(ROOT)}: expected one patch.json")
            if generated:
                require(len(members) == 1, f"{path.relative_to(ROOT)}: generated archive must contain one member")
                require(members[0].name == "patch.json", f"{path.relative_to(ROOT)}: non-canonical patch path")
            member = patch_members[0]
            payload_file = archive.extractfile(member)
            require(payload_file is not None, f"{path.relative_to(ROOT)}: unreadable patch.json")
            payload = payload_file.read()
            patch = json.loads(payload)
            return patch, member, payload
    except (tarfile.TarError, json.JSONDecodeError) as exc:
        raise ValidationError(f"{path.relative_to(ROOT)}: malformed archive: {exc}") from exc


def validate_patch(path, generated):
    patch, member, payload = read_patch(path, generated)
    if generated:
        require(member.mtime == 0, f"{path.name}: non-canonical tar mtime")
        require(member.mode == 0o644, f"{path.name}: non-canonical tar mode")
        require(member.uid == 0 and member.gid == 0, f"{path.name}: non-canonical tar owner")
        require(member.uname == "" and member.gname == "", f"{path.name}: non-canonical tar owner name")
        require(payload == patch_json_bytes(patch), f"{path.name}: non-canonical patch.json encoding")

    modules = patch.get("modules")
    cables = patch.get("cables")
    require(isinstance(modules, list), f"{path.name}: modules is not an array")
    require(isinstance(cables, list), f"{path.name}: cables is not an array")
    module_ids = [module.get("id") for module in modules]
    cable_ids = [cable.get("id") for cable in cables]
    require(len(module_ids) == len(set(module_ids)), f"{path.name}: duplicate module ID")
    require(len(cable_ids) == len(set(cable_ids)), f"{path.name}: duplicate cable ID")
    require(all(isinstance(value, int) for value in module_ids), f"{path.name}: invalid module ID")
    require(all(isinstance(value, int) for value in cable_ids), f"{path.name}: invalid cable ID")
    by_id = {module["id"]: module for module in modules}
    master_id = patch.get("masterModuleId", -1)
    require(
        isinstance(master_id, int) and (master_id == -1 or master_id in by_id),
        f"{path.name}: invalid masterModuleId",
    )

    for module in modules:
        pos = module.get("pos")
        require(
            isinstance(pos, list) and len(pos) == 2 and all(isinstance(v, (int, float)) for v in pos),
            f"{path.name}: module {module['id']} has invalid position",
        )
        if generated:
            require(pos[0] >= 0 and pos[1] >= 0, f"{path.name}: module {module['id']} has negative position")
            key = (module.get("plugin"), module.get("model"))
            require(key in MODULES, f"{path.name}: no validation data for {key[0]}/{key[1]}")

    for cable in cables:
        output_module = cable.get("outputModuleId")
        input_module = cable.get("inputModuleId")
        output_id = cable.get("outputId")
        input_id = cable.get("inputId")
        require(output_module in by_id, f"{path.name}: cable references missing output module")
        require(input_module in by_id, f"{path.name}: cable references missing input module")
        require(isinstance(output_id, int) and output_id >= 0, f"{path.name}: invalid output port")
        require(isinstance(input_id, int) and input_id >= 0, f"{path.name}: invalid input port")
        output_key = (by_id[output_module].get("plugin"), by_id[output_module].get("model"))
        input_key = (by_id[input_module].get("plugin"), by_id[input_module].get("model"))
        if output_key in MODULES:
            require(output_id < MODULES[output_key][2], f"{path.name}: output port out of range")
        if input_key in MODULES:
            require(input_id < MODULES[input_key][1], f"{path.name}: input port out of range")

    if generated:
        rows = {}
        for module in modules:
            rows.setdefault(module["pos"][1], []).append(module)
        for row, row_modules in rows.items():
            ordered = sorted(row_modules, key=lambda module: module["pos"][0])
            for left, right in zip(ordered, ordered[1:]):
                key = (left["plugin"], left["model"])
                left_end = left["pos"][0] + MODULES[key][0]
                require(
                    left_end <= right["pos"][0],
                    f"{path.name}: {key[0]}/{key[1]} at {left['pos']} overlaps "
                    f"{right['plugin']}/{right['model']} at {right['pos']} on row {row}",
                )
    return len(modules), len(cables)


def main():
    try:
        module_count = validate_manifest()
        markdown_files, link_count = validate_markdown()
        generated = sorted((ROOT / "patches").glob("*.vcv"))
        community = sorted((ROOT / "patches/community").glob("*.vcv"))
        require(len(generated) == 49, f"expected 49 generated patches, found {len(generated)}")
        total_modules = 0
        total_cables = 0
        for path in generated:
            modules, cables = validate_patch(path, generated=True)
            total_modules += modules
            total_cables += cables
        for path in community:
            modules, cables = validate_patch(path, generated=False)
            total_modules += modules
            total_cables += cables
    except (OSError, ValidationError, ValueError) as exc:
        print(f"asset validation failed: {exc}", file=sys.stderr)
        return 1

    print(
        f"asset validation: PASS ({module_count} manifest modules, "
        f"{markdown_files} Markdown files/{link_count} local links, "
        f"{len(generated)} generated + {len(community)} community patches, "
        f"{total_modules} patch modules/{total_cables} cables)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
