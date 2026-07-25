#!/usr/bin/env python3
"""Foxes demo patches. Minimal (Fundamental VCMixer + Core AudioInterface only).
Foxes params:  0 RATE, 1 BALANCE, 2 WILD, 3 BAL_ATT, 4 WILD_ATT
Foxes outputs: 0 GRASS, 1 BUNNY, 2 FOX, 3 GRASS_PEAK, 4 BUNNY_PEAK, 5 FOX_PEAK
WILD → b1 = 1 + 5.2·wild²:  0.50→b1=2.3 (periodic default), 0.62→b1≈3.0 (canonical
chaos), 0.57→b1≈2.7 (period-doubling/transition). BALANCE 0.5 → b2=2.0 (canonical).
"""
import json, os, shutil, sys, random

from patch_utils import windows_patch_directories, write_patch_archive
random.seed(11)
def uid(): return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)
ROOT = os.path.dirname(os.path.abspath(__file__))
VER = json.load(open(os.path.join(ROOT, "..", "plugin.json")))["version"]

def foxes(params, pos):
    return {"id": uid(), "plugin": "Coalescent", "model": "Foxes", "version": VER,
            "params": [{"id": i, "value": float(v)} for i, v in params], "pos": pos}
def vcmixer(levels, pos):
    # param 0 = master, 1..4 = channel levels
    return {"id": uid(), "plugin": "Fundamental", "model": "VCMixer", "version": "2.6.4",
            "params": [{"id": i, "value": float(v)} for i, v in enumerate([1.0] + levels + [0.0] * (4 - len(levels)))],
            "pos": pos}
def audio(pos):
    return {"id": uid(), "plugin": "Core", "model": "AudioInterface", "version": "2.6.6", "params": [], "pos": pos,
            "data": {"audio": {"driver": -1, "deviceName": "", "sampleRate": 44100.0,
                               "blockSize": 256, "inputOffset": 0, "outputOffset": 0}, "dcFilter": True}}
def cable(om, oid, im, iid, ci):
    cols = ["#8ad06a", "#f3a0a8", "#f0a83c", "#3695ef"]   # grass, bunny, fox, blue
    return {"id": uid(), "outputModuleId": om, "outputId": oid, "inputModuleId": im,
            "inputId": iid, "color": cols[ci % len(cols)], "inputPlugOrder": ci, "outputPlugOrder": ci}
def write(name, mods, cables, master):
    patch = {"version": "2.6.6", "zoom": 0.5, "gridOffset": [0.0, 0.0], "modules": mods, "cables": cables, "masterModuleId": master}
    out = os.path.join(ROOT, "..", "patches", name)
    write_patch_archive(out, patch)
    print(f"  {name}: {len(mods)} modules, {len(cables)} cables")
    for w in windows_patch_directories():
        try: shutil.copy2(out, os.path.join(w, name))   # optional convenience copy — never fail generation on it
        except OSError as e: print(f"    (skipped Windows copy: {e})", file=sys.stderr)

# Three population outputs → a conservative mono mix → audio L+R. Shared by the
# tone demos (out ids differ only in the source: populations vs peak events).
def mix3(fx, out_ids, levels, name):
    mx = vcmixer(levels, [13, 0]); a = audio([24, 0])
    cs = [cable(fx["id"], out_ids[i], mx["id"], i + 1, i) for i in range(3)]
    cs += [cable(mx["id"], 0, a["id"], 0, 3), cable(mx["id"], 0, a["id"], 1, 3)]
    write(name, [fx, mx, a], cs, a["id"])

# 1. Food chain: the periodic default (WILD 0.5). GRASS/BUNNY/FOX mixed at
#    conservative levels — the regular three-population chase as a drifting drone.
def p_food_chain():
    fx = foxes([(0, 0.0), (1, 0.5), (2, 0.5)], [0, 0])
    mix3(fx, [0, 1, 2], [0.5, 0.5, 0.5], "foxes_1_food_chain.vcv")

# 2. Teacup: canonical chaos (WILD 0.62 → b1≈3). Slightly below C4 for a textural,
#    broadband voice — the deterministic strange attractor as irregular correlated tone.
def p_teacup():
    fx = foxes([(0, -2.0), (1, 0.5), (2, 0.62)], [0, 0])
    mix3(fx, [0, 1, 2], [0.45, 0.45, 0.45], "foxes_2_teacup.vcv")

# 3. Transition: a period-doubling/multipeak point (WILD 0.57 → b1≈2.7) between the
#    regular chase and full chaos.
def p_transition():
    fx = foxes([(0, -1.0), (1, 0.5), (2, 0.57)], [0, 0])
    mix3(fx, [0, 1, 2], [0.5, 0.5, 0.5], "foxes_3_transition.vcv")

# 4. Events: slow RATE, the three peak gates (G^/B^/F^) as offset ticks. The manual
#    explains that users normally route these to envelopes/sequencers, not straight
#    to audio — here they audition as a rhythmic pattern that turns irregular in chaos.
def p_events():
    fx = foxes([(0, -6.0), (1, 0.5), (2, 0.5)], [0, 0])
    mix3(fx, [3, 4, 5], [0.6, 0.6, 0.6], "foxes_4_events.vcv")

if __name__ == "__main__":
    print("Generating Foxes demo patches:")
    p_food_chain(); p_teacup(); p_transition(); p_events()
    print("Done.")
