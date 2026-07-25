#!/usr/bin/env python3
"""Shared deterministic writer for generated Rack patch archives."""

import glob
import io
import json
import os
import subprocess
import tarfile


def patch_json_bytes(patch):
    """Return the canonical UTF-8 representation stored in generated patches."""
    return (json.dumps(patch, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


def patch_tar_bytes(patch):
    """Build a tar containing only patch.json with reproducible metadata."""
    payload = patch_json_bytes(patch)
    archive = io.BytesIO()
    with tarfile.open(fileobj=archive, mode="w:", format=tarfile.USTAR_FORMAT) as tar:
        info = tarfile.TarInfo("patch.json")
        info.size = len(payload)
        info.mtime = 0
        info.mode = 0o644
        info.uid = 0
        info.gid = 0
        info.uname = ""
        info.gname = ""
        tar.addfile(info, io.BytesIO(payload))
    return archive.getvalue()


def write_patch_archive(path, patch, compression_level=19):
    """Write a deterministic tar+zstd Rack patch archive."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    result = subprocess.run(
        ["zstd", f"-{compression_level}", "-q", "-f", "-o", path],
        input=patch_tar_bytes(patch),
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        message = result.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"zstd failed while writing {path}: {message}")


def windows_patch_directories():
    """Return optional WSL Rack install targets unless explicitly disabled."""
    if os.environ.get("COALESCENT_SKIP_PATCH_INSTALL"):
        return []
    return glob.glob("/mnt/c/Users/*/AppData/Local/Rack2/patches")
