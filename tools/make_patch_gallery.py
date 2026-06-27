#!/usr/bin/env python3
"""Gallery patch: all four Coalescent modules side by side for a quick visual
check (panels, labels, layout, live scopes). No cables — each module free-runs
on its configured default voicing.

Widths: GENDYN / Axon / Soma = 12 HP, Haptik = 18 HP.
"""

import json, os, io, glob, shutil, subprocess, sys, tarfile, tempfile, random

random.seed(40)
def uid():
    return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)

ROOT = os.path.dirname(os.path.abspath(__file__))
VER = json.load(open(os.path.join(ROOT, "..", "plugin.json")))["version"]

def mod(model, x):
    # params=[] → each module uses its configured defaults.
    return {"id": uid(), "plugin": "Coalescent", "model": model,
            "version": VER, "params": [], "pos": [x, 0]}

modules = [mod("GENDYN", 0), mod("Haptik", 12), mod("Axon", 30), mod("Soma", 42)]

patch = {"version": "2.6.6", "zoom": 0.5, "gridOffset": [0.0, 0.0],
         "modules": modules, "cables": [], "masterModuleId": -1}

name = "coalescent_gallery.vcv"
out_dir = os.path.join(ROOT, "..", "patches")
os.makedirs(out_dir, exist_ok=True)
out_file = os.path.join(out_dir, name)
with tempfile.TemporaryDirectory() as tmp:
    jp = os.path.join(tmp, "patch.json")
    with open(jp, "w") as f:
        json.dump(patch, f, indent=2)
    tar_buf = io.BytesIO()
    with tarfile.open(fileobj=tar_buf, mode="w:") as tf:
        tf.add(jp, arcname="patch.json")
    r = subprocess.run(["zstd", "-19", "-o", out_file, "-f"],
                       input=tar_buf.getvalue(), capture_output=True)
    if r.returncode != 0:
        print("zstd error:", r.stderr.decode(), file=sys.stderr); sys.exit(1)

print(f"  {name}: {len(modules)} modules (v{VER}), {os.path.getsize(out_file)} bytes")
for win in glob.glob("/mnt/c/Users/*/AppData/Local/Rack2/patches"):
    shutil.copy2(out_file, os.path.join(win, name))
    print(f"    installed -> {win}/{name}")
