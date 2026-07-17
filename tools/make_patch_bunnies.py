#!/usr/bin/env python3
"""Bunnies demo patches. Minimal (Fundamental VCMixer + Core AudioInterface only).
Bunnies params: 0 RATE, 1 BALANCE, 2 WILD, 3 MODE(0=LV,1=RM), 4 BAL_ATT, 5 WILD_ATT
Bunnies outputs: 0 PREY, 1 PRED, 2 PREY_POP, 3 PRED_POP
"""
import json, os, shutil, random

from patch_utils import windows_patch_directories, write_patch_archive
random.seed(7)
def uid(): return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)
ROOT = os.path.dirname(os.path.abspath(__file__))
VER = json.load(open(os.path.join(ROOT, "..", "plugin.json")))["version"]

def bunnies(params, pos):
    return {"id": uid(), "plugin": "Coalescent", "model": "Bunnies", "version": VER,
            "params": [{"id": i, "value": float(v)} for i, v in params], "pos": pos}
def vcmixer(pos):
    return {"id": uid(), "plugin": "Fundamental", "model": "VCMixer", "version": "2.6.4",
            "params": [{"id": i, "value": 1.0} for i in range(5)], "pos": pos}
def audio(pos):
    return {"id": uid(), "plugin": "Core", "model": "AudioInterface", "version": "2.6.6", "params": [], "pos": pos,
            "data": {"audio": {"driver": -1, "deviceName": "", "sampleRate": 44100.0,
                               "blockSize": 256, "inputOffset": 0, "outputOffset": 0}, "dcFilter": True}}
def cable(om, oid, im, iid, ci):
    cols = ["#f3374b", "#ffb437", "#00b56e", "#3695ef"]
    return {"id": uid(), "outputModuleId": om, "outputId": oid, "inputModuleId": im,
            "inputId": iid, "color": cols[ci % len(cols)], "inputPlugOrder": ci, "outputPlugOrder": ci}
def write(name, mods, cables, master):
    patch = {"version": "2.6.6", "zoom": 0.5, "gridOffset": [0.0, 0.0], "modules": mods, "cables": cables, "masterModuleId": master}
    out = os.path.join(ROOT, "..", "patches", name)
    write_patch_archive(out, patch)
    print(f"  {name}: {len(mods)} modules, {len(cables)} cables")
    for w in windows_patch_directories():
        try: shutil.copy2(out, os.path.join(w, name))
        except OSError as error: print(f"    (skipped Windows copy: {error})")

# 1. LV predator-prey pair: PREY/PRED hard-panned — the quarter-cycle chase in stereo.
def p_lv():
    b = bunnies([(0, 0.0), (1, 0.5), (2, 0.4), (3, 0.0)], [0, 0]); a = audio([14, 0])
    write("bunnies_1_lv_pair.vcv", [b, a], [cable(b["id"], 0, a["id"], 0, 0), cable(b["id"], 1, a["id"], 1, 2)], a["id"])

# 2. RM enrichment cycle: raise WILD from low → silence blooms into a limit cycle.
def p_rm():
    b = bunnies([(0, 0.0), (1, 0.5), (2, 0.6), (3, 1.0)], [0, 0]); mx = vcmixer([14, 0]); a = audio([25, 0])
    cs = [cable(b["id"], 0, mx["id"], 1, 0), cable(b["id"], 1, mx["id"], 2, 1),
          cable(mx["id"], 0, a["id"], 0, 2), cable(mx["id"], 0, a["id"], 1, 2)]
    write("bunnies_2_rm_cycle.vcv", [b, mx, a], cs, a["id"])

# 3. Two-phase clock: low rate (LV), PREY_POP/PRED_POP → offset ticks.
def p_clock():
    b = bunnies([(0, -6.0), (1, 0.5), (2, 0.5), (3, 0.0)], [0, 0]); mx = vcmixer([14, 0]); a = audio([25, 0])
    cs = [cable(b["id"], 2, mx["id"], 1, 0), cable(b["id"], 3, mx["id"], 2, 1),
          cable(mx["id"], 0, a["id"], 0, 2), cable(mx["id"], 0, a["id"], 1, 2)]
    write("bunnies_3_clock.vcv", [b, mx, a], cs, a["id"])

if __name__ == "__main__":
    print("Generating Bunnies demo patches:")
    p_lv(); p_rm(); p_clock()
    print("Done.")
