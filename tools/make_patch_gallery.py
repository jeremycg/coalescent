#!/usr/bin/env python3
"""Gallery/showcase patch: all four Coalescent modules in one row, each doing
its best trick — for screenshots and a quick visual check.

- Axon: 4 poly voices (Cmaj7) with CURRENT spread via 8vert->Merge, so the scope
  draws four differently-sized coloured orbits (values proven in axon_6_poly).
- Soma: 4 poly voices (Cmin7) with currents walking tonic->chaos (soma_6_poly).
- GENDYN: sine seed with SCALE turned up so the polygon visibly morphs, N=24.
- Haptik: continuous drive, low damp, N=32 so the mass nodes render.

Top row (y=0) holds ONLY the Coalescent modules; all utility modules (8vert /
Merge poly sources) live on the second row (y=1) so a screenshot can crop to
the top row cleanly. No audio cables — the displays animate regardless.
"""

import json, math, os, io, glob, shutil, subprocess, sys, tarfile, tempfile, random

random.seed(40)
def uid():
    return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)

ROOT = os.path.dirname(os.path.abspath(__file__))
VER = json.load(open(os.path.join(ROOT, "..", "plugin.json")))["version"]

def mod(model, x, y=0, params=None, data=None):
    m = {"id": uid(), "plugin": "Coalescent", "model": model,
         "version": VER, "params": params or [], "pos": [x, y]}
    if data: m["data"] = data
    return m

def pv(i, v):
    return {"id": i, "value": float(v)}

def eightvert(gains, x, y=1):
    params = [pv(i, gains[i] if i < len(gains) else 0.0) for i in range(8)]
    return {"id": uid(), "plugin": "Fundamental", "model": "8vert", "version": "2.6.4",
            "params": params, "pos": [x, y]}

def merge(x, y=1):
    return {"id": uid(), "plugin": "Fundamental", "model": "Merge", "version": "2.6.4",
            "params": [], "pos": [x, y]}

def cable(om, oid, im, iid, ci):
    colors = ["#f3374b", "#ffb437", "#00b56e", "#3695ef"]
    return {"id": uid(), "outputModuleId": om, "outputId": oid,
            "inputModuleId": im, "inputId": iid, "color": colors[ci % len(colors)],
            "inputPlugOrder": ci, "outputPlugOrder": ci}

# ── Top row: the four Coalescent modules, dressed to impress ──────────────────
# GENDYN: N=24 breakpoints, SCALE up so the sine seed morphs into a living
# polygon within seconds (params: 0 N, 1 SCALE_AMP, 2 SCALE_DUR).
gendyn = mod("GENDYN", 0, 0, params=[pv(0, 24), pv(1, 0.006), pv(2, 0.006)])
# Haptik: drive excitation, low damp, N=32 (nodes render at N<=48).
# (params: 0 N, 4 DAMP, 5 INJECT, 6 EXCITE)
haptik = mod("Haptik", 12, 0, params=[pv(0, 32), pv(4, 0.10), pv(5, 0.7), pv(6, 3)])
# Axon: defaults + CURRENT attenuverter open for the poly spread.
# (params: 1 CURRENT, 2 EPS, 3 SHAPE, 4 CURRENT_ATT)
axon = mod("Axon", 30, 0, params=[pv(1, 0.6), pv(2, 0.08), pv(3, 0.7), pv(4, 1.0)])
# Soma: bursting BURST rate so the current spread walks tonic -> chaos.
# (params: 1 CURRENT, 2 BURST=log2 r, 3 ADAPT, 4 CURRENT_ATT)
soma = mod("Soma", 42, 0, params=[pv(1, 2.0), pv(2, math.log2(0.006)), pv(3, 4.0), pv(4, 1.0)])

# ── Second row: poly sources (proven values from axon_6_poly / soma_6_poly) ──
axon_semis, axon_cur = [0, 4, 7, 11], [-1.5, 1.0, 4.0, 7.5]     # I ≈ 0.45..1.35
soma_semis, soma_cur = [0, 3, 7, 10], [-2.5, 0.0, 4.0, 6.25]    # I ≈ 1.5..3.25
evPa = eightvert([s / 120.0 for s in axon_semis], 0)
evCa = eightvert([v / 10.0 for v in axon_cur], 8)
mgPa, mgCa = merge(16), merge(18)
evPs = eightvert([s / 120.0 for s in soma_semis], 24)
evCs = eightvert([v / 10.0 for v in soma_cur], 32)
mgPs, mgCs = merge(40), merge(42)

modules = [gendyn, haptik, axon, soma, evPa, evCa, mgPa, mgCa, evPs, evCs, mgPs, mgCs]

cables  = [cable(evPa["id"], i, mgPa["id"], i, i) for i in range(4)]
cables += [cable(evCa["id"], i, mgCa["id"], i, i) for i in range(4)]
cables += [cable(mgPa["id"], 0, axon["id"], 0, 0),   # -> Axon V/OCT (4ch)
           cable(mgCa["id"], 0, axon["id"], 1, 1)]   # -> Axon CURRENT CV
cables += [cable(evPs["id"], i, mgPs["id"], i, i) for i in range(4)]
cables += [cable(evCs["id"], i, mgCs["id"], i, i) for i in range(4)]
cables += [cable(mgPs["id"], 0, soma["id"], 0, 2),   # -> Soma V/OCT (4ch)
           cable(mgCs["id"], 0, soma["id"], 1, 3)]   # -> Soma CURRENT CV

patch = {"version": "2.6.6", "zoom": 0.5, "gridOffset": [0.0, 0.0],
         "modules": modules, "cables": cables, "masterModuleId": -1}

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

print(f"  {name}: {len(modules)} modules, {len(cables)} cables (v{VER}), {os.path.getsize(out_file)} bytes")
for win in glob.glob("/mnt/c/Users/*/AppData/Local/Rack2/patches"):
    shutil.copy2(out_file, os.path.join(win, name))
    print(f"    installed -> {win}/{name}")
