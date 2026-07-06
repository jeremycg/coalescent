#!/usr/bin/env python3
"""Operon smoke-test / demo patches. Minimal, dependency-light (only Fundamental
VCMixer + Core AudioInterface, both universal), so they load anywhere and are
easy to verify. Port maps checked against Fundamental source:
  VCMixer  inputs: 0=MIX_CV, 1-4=CH audio, 5-8=level CV ; outputs: 0=MIX, 1-4=CH
  AudioInterface inputs: 0=L(ch1), 1=R(ch2)
Operon params: 0 PITCH, 1 ALPHA, 2 HILL, 3 BETA, 4 LEAK, 5-7 att
Operon outputs: 0-2 OUT1/2/3, 3-5 GATE1/2/3
"""
import json, os, io, glob, shutil, subprocess, sys, tarfile, tempfile, random

random.seed(12)
def uid(): return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)
ROOT = os.path.dirname(os.path.abspath(__file__))
VER = json.load(open(os.path.join(ROOT, "..", "plugin.json")))["version"]

def operon(params, pos):
    return {"id": uid(), "plugin": "Coalescent", "model": "Operon", "version": VER,
            "params": [{"id": i, "value": float(v)} for i, v in params], "pos": pos}
def vcmixer(pos):
    return {"id": uid(), "plugin": "Fundamental", "model": "VCMixer", "version": "2.6.4",
            "params": [{"id": i, "value": 1.0} for i in range(5)], "pos": pos}
def audio(pos):
    return {"id": uid(), "plugin": "Core", "model": "AudioInterface", "version": "2.6.6",
            "params": [], "pos": pos,
            "data": {"audio": {"driver": -1, "deviceName": "", "sampleRate": 44100.0,
                               "blockSize": 256, "inputOffset": 0, "outputOffset": 0},
                     "dcFilter": True}}
def cable(om, oid, im, iid, ci):
    cols = ["#f3374b", "#ffb437", "#00b56e", "#3695ef", "#e6e6e6"]
    return {"id": uid(), "outputModuleId": om, "outputId": oid, "inputModuleId": im,
            "inputId": iid, "color": cols[ci % len(cols)], "inputPlugOrder": ci, "outputPlugOrder": ci}

def write_patch(name, modules, cables, master):
    patch = {"version": "2.6.6", "zoom": 0.5, "gridOffset": [0.0, 0.0],
             "modules": modules, "cables": cables, "masterModuleId": master}
    out = os.path.join(ROOT, "..", "patches", name)
    with tempfile.TemporaryDirectory() as tmp:
        jp = os.path.join(tmp, "patch.json"); json.dump(patch, open(jp, "w"), indent=2)
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:") as tf: tf.add(jp, arcname="patch.json")
        r = subprocess.run(["zstd", "-19", "-o", out, "-f"], input=buf.getvalue(), capture_output=True)
        if r.returncode: print("zstd error:", r.stderr.decode(), file=sys.stderr); sys.exit(1)
    print(f"  {name}: {len(modules)} modules, {len(cables)} cables")
    for win in glob.glob("/mnt/c/Users/*/AppData/Local/Rack2/patches"):
        shutil.copy2(out, os.path.join(win, name))

DEF = [(0, 0.0), (1, 12.0), (2, 2.5), (3, 1.0), (4, 0.05)]   # default voicing

# 1. Three-phase tone: OUT1/2/3 summed through a mixer.
def p_threephase():
    op = operon(DEF, [0, 0]); mx = vcmixer([14, 0]); a = audio([25, 0])
    cs = [cable(op["id"], k, mx["id"], k + 1, k) for k in range(3)]
    cs += [cable(mx["id"], 0, a["id"], 0, 3), cable(mx["id"], 0, a["id"], 1, 3)]
    write_patch("operon_1_threephase.vcv", [op, mx, a], cs, a["id"])

# 2. Stereo: two phases (OUT1, OUT3) hard-panned L/R for a wide three-phase field.
def p_stereo():
    op = operon(DEF, [0, 0]); a = audio([14, 0])
    cs = [cable(op["id"], 0, a["id"], 0, 0), cable(op["id"], 2, a["id"], 1, 2)]
    write_patch("operon_2_stereo.vcv", [op, a], cs, a["id"])

# 3. Relaxation pulses: high HILL + high ALPHA → sharp three-phase (watch aliasing high).
def p_pulses():
    op = operon([(0, 0.0), (1, 40.0), (2, 6.0), (3, 0.6), (4, 0.05)], [0, 0])
    mx = vcmixer([14, 0]); a = audio([25, 0])
    cs = [cable(op["id"], k, mx["id"], k + 1, k) for k in range(3)]
    cs += [cable(mx["id"], 0, a["id"], 0, 3), cable(mx["id"], 0, a["id"], 1, 3)]
    write_patch("operon_3_pulses.vcv", [op, mx, a], cs, a["id"])

# 4. Three-phase clock: low pitch (~2 Hz), GATE1/2/3 summed → a three-phase tick.
#    Route the gates to your own envelopes/sequencers for a real three-phase clock.
def p_clock():
    op = operon([(0, -7.0), (1, 12.0), (2, 4.0), (3, 1.0), (4, 0.05)], [0, 0])
    mx = vcmixer([14, 0]); a = audio([25, 0])
    cs = [cable(op["id"], 3 + k, mx["id"], k + 1, k) for k in range(3)]   # GATE1/2/3
    cs += [cable(mx["id"], 0, a["id"], 0, 3), cable(mx["id"], 0, a["id"], 1, 3)]
    write_patch("operon_4_clock.vcv", [op, mx, a], cs, a["id"])

if __name__ == "__main__":
    print("Generating Operon demo patches:")
    p_threephase(); p_stereo(); p_pulses(); p_clock()
    print("Done.")
