#!/usr/bin/env python3
"""Generate the four §7 Haptik smoke-test patches as .vcv files.

Each patch is Haptik → Core AudioInterface (L+R). The driven-resonator patch
also wires a Fundamental VCO SAW output into Haptik's EXT IN.

Haptik positional ids (must match enum order in src/Haptik.cpp):
  ParamId : 0 N, 1 PITCH, 2 RATE(log2 Hz), 3 COUPLE, 4 DAMP, 5 INJECT,
            6 EXCITE, 7 FREEZE, 8 RATE_ATT, 9 COUPLE_ATT, 10 DAMP_ATT, 11 INJECT_ATT,
            12 DRIVER, 13 MODE
  InputId : 0 VOCT, 1 RATE, 2 COUPLE, 3 DAMP, 4 INJECT, 5 TRIG, 6 EXT
  OutputId: 0 OUT, 1 MOTION
Fundamental VCO (from source): FREQ_PARAM=2, PITCH_INPUT=0, SAW_OUTPUT=2.
"""

import json, math, os, io, glob, shutil, subprocess, sys, tarfile, tempfile, random
ROOT = os.path.dirname(os.path.abspath(__file__))
VER = json.load(open(os.path.join(ROOT, "..", "plugin.json")))["version"]  # Coalescent version from the manifest

random.seed(7)
def uid():
    return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)

LOG2_3 = math.log2(3.0)   # RATE = 3 Hz expressed as log2(Hz)

def haptik(params, pos):
    return {"id": uid(), "plugin": "Coalescent", "model": "Haptik",
            "version": VER, "params": params, "pos": pos}

def hp(N=64, pitch=0.0, rate=LOG2_3, couple=0.3, damp=0.35, inject=0.6,
       excite=1, freeze=0.0, driver=0.25, mode=0.0):
    return [
        {"id": 0,  "value": float(N)},
        {"id": 1,  "value": float(pitch)},
        {"id": 2,  "value": float(rate)},
        {"id": 3,  "value": float(couple)},
        {"id": 4,  "value": float(damp)},
        {"id": 5,  "value": float(inject)},
        {"id": 6,  "value": float(excite)},
        {"id": 7,  "value": float(freeze)},
        {"id": 8,  "value": 0.0}, {"id": 9, "value": 0.0},
        {"id": 10, "value": 0.0}, {"id": 11, "value": 0.0},
        {"id": 12, "value": float(driver)},   # DRIVER position (added after v1)
        {"id": 13, "value": float(mode)},      # MODE: 0 Fast / 1 Slow
    ]

def audio(pos):
    return {"id": uid(), "plugin": "Core", "model": "AudioInterface",
            "version": "2.6.6", "params": [],
            "data": {"audio": {"driver": -1, "deviceName": "", "sampleRate": 44100.0,
                               "blockSize": 256, "inputOffset": 0, "outputOffset": 0},
                     "dcFilter": True},
            "pos": pos}

def cable(om, oid, im, iid, ci):
    colors = ["#f3374b", "#ffb437", "#00b56e", "#3695ef"]
    return {"id": uid(), "outputModuleId": om, "outputId": oid,
            "inputModuleId": im, "inputId": iid, "color": colors[ci % len(colors)],
            "inputPlugOrder": ci, "outputPlugOrder": ci}

def write_patch(name, modules, cables, master_id):
    patch = {"version": "2.6.6", "zoom": 0.5, "gridOffset": [0.0, 0.0],
             "modules": modules, "cables": cables, "masterModuleId": master_id}
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "patches")
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
    print(f"  {name}: {len(modules)} modules, {len(cables)} cables, {os.path.getsize(out_file)} bytes")
    for win in glob.glob("/mnt/c/Users/*/AppData/Local/Rack2/patches"):
        shutil.copy2(out_file, os.path.join(win, name))
        print(f"    installed -> {win}/{name}")

# ── 1. Plucked, settling ──────────────────────────────────────────────────────
def patch_plucked():
    h = haptik(hp(excite=1, damp=0.35, couple=0.3), [0, 0])
    a = audio([18, 0])   # Haptik is 18 HP wide; place downstream at x >= 18
    cs = [cable(h["id"], 0, a["id"], 0, 0), cable(h["id"], 0, a["id"], 1, 1)]
    write_patch("haptik_1_plucked.vcv", [h, a], cs, a["id"])

# ── 2. Drone (DAMP=0 → seed bump rings forever) ───────────────────────────────
def patch_drone():
    h = haptik(hp(excite=2, damp=0.0, couple=0.6), [0, 0])
    a = audio([18, 0])
    cs = [cable(h["id"], 0, a["id"], 0, 0), cable(h["id"], 0, a["id"], 1, 1)]
    write_patch("haptik_2_drone.vcv", [h, a], cs, a["id"])

# ── 3. Frozen wavetable (FREEZE on; play V/OCT) ───────────────────────────────
def patch_frozen():
    h = haptik(hp(excite=1, damp=0.35, couple=0.3, freeze=1.0), [0, 0])
    a = audio([18, 0])
    cs = [cable(h["id"], 0, a["id"], 0, 0), cable(h["id"], 0, a["id"], 1, 1)]
    write_patch("haptik_3_frozen.vcv", [h, a], cs, a["id"])

# ── 4. Driven resonator (VCO SAW → EXT IN) ────────────────────────────────────
def patch_driven():
    vco = {"id": uid(), "plugin": "Fundamental", "model": "VCO", "version": "2.6.4",
           "params": [{"id": 2, "value": 0.0}], "pos": [0, 0]}
    h = haptik(hp(excite=1, damp=0.4, couple=0.3, inject=0.6), [10, 0])  # VCO is 10 HP
    a = audio([28, 0])   # Haptik spans 10..28
    cs = [
        cable(vco["id"], 2, h["id"], 6, 0),   # VCO SAW -> Haptik EXT IN
        cable(h["id"], 0, a["id"], 0, 1),     # Haptik OUT -> L
        cable(h["id"], 0, a["id"], 1, 2),     # Haptik OUT -> R
    ]
    write_patch("haptik_4_driven.vcv", [vco, h, a], cs, a["id"])

if __name__ == "__main__":
    print("Generating Haptik smoke-test patches:")
    patch_plucked(); patch_drone(); patch_frozen(); patch_driven()
    print("Done.")
