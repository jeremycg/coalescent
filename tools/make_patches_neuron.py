#!/usr/bin/env python3
"""Generate the §7 Axon smoke-test patches as .vcv files.

Axon positional ids (must match enum order in src/Axon.cpp):
  ParamId : 0 PITCH, 1 CURRENT, 2 EPS, 3 SHAPE, 4 CURRENT_ATT, 5 EPS_ATT
  InputId : 0 VOCT, 1 CURRENT, 2 EPS, 3 TRIG, 4 SYNC
  OutputId: 0 OUT, 1 SPIKE, 2 W
Axon is 12 HP wide → place downstream modules at x >= 12.

Third-party enum orders (verified from Fundamental v2 source):
  LFO : FREQ_PARAM=2; SQR_OUTPUT=3                (6 HP)
  VCO : FREQ_PARAM=2; PITCH_INPUT=0; SAW_OUTPUT=2; SQR_OUTPUT=3  (10 HP)
"""

import json, os, io, glob, shutil, subprocess, sys, tarfile, tempfile, random

random.seed(11)
def uid():
    return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)

PLUGIN_VERSION = "2.0.0"

def axon(params, pos):
    return {"id": uid(), "plugin": "Coalescent", "model": "Axon",
            "version": PLUGIN_VERSION, "params": params, "pos": pos}

def vco(freq=0.0, pos=(0, 0)):
    # Fundamental VCO as a hard-sync master: FREQ_PARAM=2, SQR_OUTPUT=3. 10 HP.
    return {"id": uid(), "plugin": "Fundamental", "model": "VCO", "version": "2.6.4",
            "params": [{"id": 2, "value": float(freq)}], "pos": list(pos)}

# ── Poly-source utilities (verified from Fundamental v2 source) ────────────────
#   8vert : GAIN_PARAMS 0..7; IN_INPUTS 0..7; OUT_OUTPUTS 0..7.  With an output
#           connected but its input unpatched, out = gain*10 V (a mono constant) — so
#           the gain knob doubles as a fixed V/oct source. 8 HP.
#   Merge : MONO_INPUTS 0..15; POLY_OUTPUT 0. Channels auto = last connected +1. 2 HP.
#   Sum   : LEVEL_PARAM 0 (0..1); POLY_INPUT 0; MONO_OUTPUT 0. out = sum(poly)*level. 2 HP.
def eightvert(gains, pos):
    params = [{"id": i, "value": float(gains[i] if i < len(gains) else 0.0)} for i in range(8)]
    return {"id": uid(), "plugin": "Fundamental", "model": "8vert", "version": "2.6.4",
            "params": params, "pos": list(pos)}

def merge(pos):
    return {"id": uid(), "plugin": "Fundamental", "model": "Merge", "version": "2.6.4",
            "params": [], "pos": list(pos)}

def summ(level, pos):
    return {"id": uid(), "plugin": "Fundamental", "model": "Sum", "version": "2.6.4",
            "params": [{"id": 0, "value": float(level)}], "pos": list(pos)}

# ── Playable-poly utilities (verified from Core / Fundamental v2 source) ───────
#   MIDIToCVInterface : PITCH_OUTPUT 0, GATE_OUTPUT 1. Poly channel count lives in
#       midiParser (right-click "Polyphony channels"); the user still selects their
#       MIDI device. "channels" in data is read by MidiParser.fromJson.
#   ADSR  : A/D/S/R params 0..3; GATE_INPUT 4; ENVELOPE_OUTPUT 0. Poly via GATE.
#   VCA-1 : LEVEL_PARAM 0, EXP_PARAM 1; CV_INPUT 0 (0..10 V → 0..1), IN_INPUT 1;
#       OUT_OUTPUT 0. Poly, follows IN channels.
def midicv(channels=4, pos=(0, 0)):
    return {"id": uid(), "plugin": "Core", "model": "MIDIToCVInterface", "version": "2.6.6",
            "params": [], "data": {"channels": int(channels)}, "pos": list(pos)}

def adsr(a=0.2, d=0.3, s=0.8, r=0.4, pos=(0, 0)):
    return {"id": uid(), "plugin": "Fundamental", "model": "ADSR", "version": "2.6.4",
            "params": [{"id": 0, "value": float(a)}, {"id": 1, "value": float(d)},
                       {"id": 2, "value": float(s)}, {"id": 3, "value": float(r)}],
            "pos": list(pos)}

def vca1(level=1.0, pos=(0, 0)):
    return {"id": uid(), "plugin": "Fundamental", "model": "VCA-1", "version": "2.6.4",
            "params": [{"id": 0, "value": float(level)}, {"id": 1, "value": 0.0}],  # LEVEL, EXP off
            "pos": list(pos)}

def ap(pitch=0.0, current=0.6, eps=0.08, shape=0.7, current_att=0.0, eps_att=0.0):
    return [
        {"id": 0, "value": float(pitch)},
        {"id": 1, "value": float(current)},
        {"id": 2, "value": float(eps)},
        {"id": 3, "value": float(shape)},
        {"id": 4, "value": float(current_att)},
        {"id": 5, "value": float(eps_att)},
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

# ── 1. Free-run tone: default voicing → audio ────────────────────────────────
def patch_freerun():
    x = axon(ap(current=0.6, eps=0.08, shape=0.7), [0, 0])
    a = audio([12, 0])   # Axon is 12 HP
    cs = [cable(x["id"], 0, a["id"], 0, 0), cable(x["id"], 0, a["id"], 1, 1)]
    write_patch("axon_1_freerun.vcv", [x, a], cs, a["id"])

# ── 2. Excitable blips: LFO square clocks TRIG, sub-threshold current ─────────
def patch_blips():
    lfo = {"id": uid(), "plugin": "Fundamental", "model": "LFO", "version": "2.6.4",
           "params": [{"id": 2, "value": 1.0}], "pos": [0, 0]}        # FREQ=1 → 2 Hz
    x = axon(ap(current=0.1, eps=0.08, shape=0.7), [6, 0])            # LFO is 6 HP
    a = audio([18, 0])                                                # Axon spans 6..18
    cs = [
        cable(lfo["id"], 3, x["id"], 3, 0),   # LFO SQR -> Axon TRIG
        cable(x["id"], 0, a["id"], 0, 1),     # Axon OUT -> L
        cable(x["id"], 0, a["id"], 1, 2),     # Axon OUT -> R
    ]
    write_patch("axon_2_blips.vcv", [lfo, x, a], cs, a["id"])

# ── 3. Self-evolving: W feeds CURRENT CV (rides its own recovery variable) ────
def patch_selfevolving():
    x = axon(ap(current=0.7, eps=0.08, shape=0.7, current_att=0.8), [0, 0])
    a = audio([12, 0])
    cs = [
        cable(x["id"], 2, x["id"], 1, 0),     # Axon W -> Axon CURRENT CV (self-patch)
        cable(x["id"], 0, a["id"], 0, 1),
        cable(x["id"], 0, a["id"], 1, 2),
    ]
    write_patch("axon_3_selfevolving.vcv", [x, a], cs, a["id"])

# ── 4. Cross-mod: VCO SAW into CURRENT CV for FM-like sidebands ───────────────
def patch_crossmod():
    vco = {"id": uid(), "plugin": "Fundamental", "model": "VCO", "version": "2.6.4",
           "params": [{"id": 2, "value": 0.0}], "pos": [0, 0]}        # VCO is 10 HP
    x = axon(ap(current=0.6, eps=0.08, shape=0.7, current_att=0.4), [10, 0])
    a = audio([22, 0])                                                # Axon spans 10..22
    cs = [
        cable(vco["id"], 2, x["id"], 1, 0),   # VCO SAW -> Axon CURRENT CV
        cable(x["id"], 0, a["id"], 0, 1),
        cable(x["id"], 0, a["id"], 1, 2),
    ]
    write_patch("axon_4_crossmod.vcv", [vco, x, a], cs, a["id"])

# ── 5. Hard sync: a master VCO resets Axon each cycle (classic sync sweep) ─────
def patch_sync():
    master = vco(freq=-1.0, pos=(0, 0))                              # VCO master, ~C3
    x = axon(ap(pitch=0.0, current=0.6, eps=0.08, shape=0.7), [10, 0])  # free-running
    a = audio([22, 0])                                              # Axon spans 10..22
    cs = [
        cable(master["id"], 3, x["id"], 4, 0),   # VCO SQR -> Axon SYNC
        cable(x["id"], 0, a["id"], 0, 1),         # Axon OUT -> L
        cable(x["id"], 0, a["id"], 1, 2),         # Axon OUT -> R
    ]
    write_patch("axon_5_sync.vcv", [master, x, a], cs, a["id"])

# ── 6. Polyphony: 4 voices spread across PITCH *and* CURRENT ───────────────────
# The phase portrait is pitch-invariant (pitch only rescales time), so voices at
# the same CURRENT trace the same orbit and pile up. Driving a second poly CV into
# CURRENT gives each voice a differently-sized limit cycle — that is what makes the
# coloured orbits visibly separate on the display.
def patch_poly():
    semis  = [0, 4, 7, 11]                                  # Cmaj7 pitches
    cur_cv = [-1.5, 1.0, 4.0, 7.5]                          # ×0.1 → I ≈ 0.45,0.70,1.00,1.35
    evP = eightvert([s / 120.0 for s in semis], (0, 0))    # pitch source (8 HP)
    evC = eightvert([v / 10.0 for v in cur_cv], (8, 0))    # current-CV source (8 HP)
    mgP = merge((16, 0)); mgC = merge((18, 0))             # 2 HP each
    x   = axon(ap(current=0.6, eps=0.08, shape=0.7, current_att=1.0), [20, 0])
    sm  = summ(0.25, (32, 0))
    a   = audio([34, 0])
    cs  = [cable(evP["id"], i, mgP["id"], i, i) for i in range(4)]   # pitch chord
    cs += [cable(evC["id"], i, mgC["id"], i, i) for i in range(4)]   # current spread
    cs += [
        cable(mgP["id"], 0, x["id"], 0, 0),   # Merge -> Axon V/OCT (4 ch)
        cable(mgC["id"], 0, x["id"], 1, 1),   # Merge -> Axon CURRENT CV (att = 1.0)
        cable(x["id"], 0, sm["id"], 0, 2),    # Axon OUT (poly) -> Sum
        cable(sm["id"], 0, a["id"], 0, 3),    # Sum -> L
        cable(sm["id"], 0, a["id"], 1, 3),    # Sum -> R
    ]
    write_patch("axon_6_poly.vcv", [evP, evC, mgP, mgC, x, sm, a], cs, a["id"])

# ── 7. Playable chords: MIDI->CV gates a poly VCA so notes start/stop ──────────
# Axon free-runs (a drone), so polyphony alone just sounds all voices forever.
# Gate the audio AFTER Axon: MIDI GATE -> ADSR -> VCA level, so each voice only
# sounds while its key is held. (Open in Rack, right-click the MIDI module to pick
# your device and set Polyphony channels to taste.)
def patch_midipoly():
    m  = midicv(4, (0, 0))                                  # 8 HP — pick device, poly=4
    x  = axon(ap(current=0.6, eps=0.08, shape=0.7), [8, 0])   # 12 HP → 4 voices
    en = adsr(0.2, 0.3, 0.8, 0.4, (20, 0))                 # 6 HP — gentle pad
    vc = vca1(1.0, (26, 0))                                 # 3 HP
    sm = summ(0.3, (29, 0))                                 # 2 HP
    a  = audio([31, 0])
    cs = [
        cable(m["id"],  0, x["id"],  0, 0),   # MIDI PITCH -> Axon V/OCT
        cable(m["id"],  1, en["id"], 4, 1),   # MIDI GATE  -> ADSR GATE
        cable(x["id"],  0, vc["id"], 1, 2),   # Axon OUT   -> VCA IN
        cable(en["id"], 0, vc["id"], 0, 3),   # ADSR ENV   -> VCA level CV
        cable(vc["id"], 0, sm["id"], 0, 4),   # VCA OUT    -> Sum (poly->mono)
        cable(sm["id"], 0, a["id"],  0, 5),   # Sum -> L
        cable(sm["id"], 0, a["id"],  1, 5),   # Sum -> R
    ]
    write_patch("axon_7_midipoly.vcv", [m, x, en, vc, sm, a], cs, a["id"])

# ── Soma (Hindmarsh-Rose) ─────────────────────────────────────────────────────
# Soma positional ids (match enum order in src/Soma.cpp):
#   ParamId : 0 PITCH, 1 CURRENT, 2 BURST(=log2 r), 3 ADAPT, 4 CURRENT_ATT, 5 BURST_ATT
#   InputId : 0 VOCT, 1 CURRENT, 2 BURST, 3 TRIG, 4 SYNC
#   OutputId: 0 OUT, 1 SPIKE, 2 Z
import math

def soma(params, pos):
    return {"id": uid(), "plugin": "Coalescent", "model": "Soma",
            "version": PLUGIN_VERSION, "params": params, "pos": pos}

def sp(pitch=0.0, current=2.0, r=0.006, adapt=4.0, current_att=0.0, burst_att=0.0):
    return [
        {"id": 0, "value": float(pitch)},
        {"id": 1, "value": float(current)},
        {"id": 2, "value": float(math.log2(r))},   # BURST stores log2(r)
        {"id": 3, "value": float(adapt)},
        {"id": 4, "value": float(current_att)},
        {"id": 5, "value": float(burst_att)},
    ]

def patch_soma_bursting():
    x = soma(sp(current=2.0, r=0.006, adapt=4.0), [0, 0])
    a = audio([12, 0])
    cs = [cable(x["id"], 0, a["id"], 0, 0), cable(x["id"], 0, a["id"], 1, 1)]
    write_patch("soma_1_bursting.vcv", [x, a], cs, a["id"])

def patch_soma_chaos():
    x = soma(sp(current=3.25, r=0.006, adapt=4.0), [0, 0])
    a = audio([12, 0])
    cs = [cable(x["id"], 0, a["id"], 0, 0), cable(x["id"], 0, a["id"], 1, 1)]
    write_patch("soma_2_chaos.vcv", [x, a], cs, a["id"])

def patch_soma_blips():
    lfo = {"id": uid(), "plugin": "Fundamental", "model": "LFO", "version": "2.6.4",
           "params": [{"id": 2, "value": -0.5}], "pos": [0, 0]}       # FREQ ~ 0.7 Hz
    x = soma(sp(current=0.6, r=0.006, adapt=4.0), [6, 0])             # sub-threshold; trig fires a burst
    a = audio([18, 0])
    cs = [
        cable(lfo["id"], 3, x["id"], 3, 0),   # LFO SQR -> Soma TRIG
        cable(x["id"], 0, a["id"], 0, 1),
        cable(x["id"], 0, a["id"], 1, 2),
    ]
    write_patch("soma_3_blips.vcv", [lfo, x, a], cs, a["id"])

def patch_soma_zmod():
    # Burst-envelope feedback: Z drives CURRENT CV, so the slow adaptation steers
    # the module between regimes — a self-evolving bursting texture.
    x = soma(sp(current=2.2, r=0.004, adapt=4.0, current_att=0.6), [0, 0])
    a = audio([12, 0])
    cs = [
        cable(x["id"], 2, x["id"], 1, 0),     # Soma Z -> Soma CURRENT CV (self-patch)
        cable(x["id"], 0, a["id"], 0, 1),
        cable(x["id"], 0, a["id"], 1, 2),
    ]
    write_patch("soma_4_zmod.vcv", [x, a], cs, a["id"])

# ── 5. Rhythmic sync: an LFO clocks SYNC, restarting each burst in time ────────
def patch_soma_sync():
    lfo = {"id": uid(), "plugin": "Fundamental", "model": "LFO", "version": "2.6.4",
           "params": [{"id": 2, "value": 1.0}], "pos": [0, 0]}       # FREQ=1 → 2 Hz clock
    x = soma(sp(current=2.0, r=0.006, adapt=4.0), [6, 0])            # bursting voicing
    a = audio([18, 0])
    cs = [
        cable(lfo["id"], 3, x["id"], 4, 0),   # LFO SQR -> Soma SYNC (restarts the burst)
        cable(x["id"], 0, a["id"], 0, 1),
        cable(x["id"], 0, a["id"], 1, 2),
    ]
    write_patch("soma_5_sync.vcv", [lfo, x, a], cs, a["id"])

# ── 6. Polyphony: 4 voices spread across PITCH *and* CURRENT (tonic→chaos) ──────
# As with Axon, CURRENT (not pitch) is what separates the (x,z) attractors. The
# four current values walk from tonic spiking up into the chaotic window, so each
# coloured voice draws a distinctly different attractor.
def patch_soma_poly():
    semis  = [0, 3, 7, 10]                                  # Cmin7 pitches
    cur_cv = [-2.5, 0.0, 4.0, 6.25]                         # ×0.2 → I ≈ 1.5,2.0,2.8,3.25
    evP = eightvert([s / 120.0 for s in semis], (0, 0))
    evC = eightvert([v / 10.0 for v in cur_cv], (8, 0))
    mgP = merge((16, 0)); mgC = merge((18, 0))
    x   = soma(sp(current=2.0, r=0.006, adapt=4.0, current_att=1.0), [20, 0])
    sm  = summ(0.3, (32, 0))
    a   = audio([34, 0])
    cs  = [cable(evP["id"], i, mgP["id"], i, i) for i in range(4)]
    cs += [cable(evC["id"], i, mgC["id"], i, i) for i in range(4)]
    cs += [
        cable(mgP["id"], 0, x["id"], 0, 0),   # Merge -> Soma V/OCT
        cable(mgC["id"], 0, x["id"], 1, 1),   # Merge -> Soma CURRENT CV (att = 1.0)
        cable(x["id"], 0, sm["id"], 0, 2),    # Soma OUT (poly) -> Sum
        cable(sm["id"], 0, a["id"], 0, 3),    # Sum -> L
        cable(sm["id"], 0, a["id"], 1, 3),    # Sum -> R
    ]
    write_patch("soma_6_poly.vcv", [evP, evC, mgP, mgC, x, sm, a], cs, a["id"])

# ── 8. Polyphonic triggers: proves per-voice poly TRIG. Three LFOs at different
# rates → Merge (a 3-channel poly gate) → Axon TRIG, with CURRENT at REST so each
# gate channel fires ONE spike on its own voice (not a free-running drone). The
# poly OUT is summed, so you hear three independent rhythms. ─────────────────────
def patch_polytrig():
    freqs = [-1.0, 0.5, 1.6]          # LFO FREQ (V) → ~1, 3, 6 Hz: three distinct rates
    lfos = [{"id": uid(), "plugin": "Fundamental", "model": "LFO", "version": "2.6.4",
             "params": [{"id": 2, "value": f}], "pos": [i * 6, 0]} for i, f in enumerate(freqs)]
    mg = merge((18, 0))
    x  = axon(ap(current=0.1, eps=0.08, shape=0.7), [20, 0])   # CURRENT=0.1 → rest (silent until triggered)
    sm = summ(0.3, (32, 0))
    a  = audio([34, 0])
    cs  = [cable(lfos[i]["id"], 3, mg["id"], i, i) for i in range(3)]  # LFO SQR → Merge ch i
    cs += [
        cable(mg["id"], 0, x["id"], 3, 3),    # Merge poly (3ch) → Axon TRIG
        cable(x["id"], 0, sm["id"], 0, 4),    # Axon OUT (poly) → Sum
        cable(sm["id"], 0, a["id"], 0, 5),    # Sum → L
        cable(sm["id"], 0, a["id"], 1, 5),    # Sum → R
    ]
    write_patch("axon_8_polytrig.vcv", lfos + [mg, x, sm, a], cs, a["id"])

# ── 9. Polyphonic gated voices: gate → a ringing NOTE per channel. Axon is the
# oscillator (a pitched voice per V/OCT channel, CURRENT in the oscillating band);
# the poly gates drive a poly ADSR→VCA so each note attacks, holds while held, and
# RINGS OUT on release. TRIG is left unpatched — it's a percussion ping, not a note
# gate; note shaping is external, exactly as after any VCO. ────────────────────────
def patch_polyvoice():
    gfreq = [-1.0, 0.4, 1.3]           # three gate rates → notes at different times
    lfos = [{"id": uid(), "plugin": "Fundamental", "model": "LFO", "version": "2.6.4",
             "params": [{"id": 2, "value": f}], "pos": [i * 6, 0]} for i, f in enumerate(gfreq)]
    mgG = merge((18, 0))                                     # 3 gates → poly gate
    ev  = eightvert([s / 120.0 for s in (0, 4, 7)], (20, 0))  # C-E-G: 3 poly V/OCT pitches
    mgP = merge((28, 0))
    x   = axon(ap(current=0.6, eps=0.08, shape=0.7), [30, 0])  # oscillating → a pitched voice
    env = adsr(a=0.01, d=0.3, s=0.6, r=0.9, pos=(42, 0))       # release 0.9 s → ring-out
    vca = vca1(1.0, (52, 0))
    sm  = summ(0.3, (56, 0))
    a   = audio([58, 0])
    cs  = [cable(lfos[i]["id"], 3, mgG["id"], i, i) for i in range(3)]   # gates → Merge
    cs += [cable(ev["id"], i, mgP["id"], i, i) for i in range(3)]        # pitches → Merge
    cs += [
        cable(mgP["id"], 0, x["id"], 0, 0),    # poly V/OCT → Axon (3 pitched voices)
        cable(mgG["id"], 0, env["id"], 4, 1),  # poly gate → ADSR GATE
        cable(x["id"], 0, vca["id"], 1, 2),    # Axon OUT (poly) → VCA IN
        cable(env["id"], 0, vca["id"], 0, 3),  # ADSR env (poly) → VCA CV → gated notes
        cable(vca["id"], 0, sm["id"], 0, 0),   # VCA out (poly) → Sum
        cable(sm["id"], 0, a["id"], 0, 1), cable(sm["id"], 0, a["id"], 1, 1),
    ]
    write_patch("axon_9_polyvoice.vcv", lfos + [mgG, ev, mgP, x, env, vca, sm, a], cs, a["id"])

if __name__ == "__main__":
    print("Generating Axon smoke-test patches:")
    patch_freerun(); patch_blips(); patch_selfevolving(); patch_crossmod()
    patch_sync(); patch_poly(); patch_midipoly(); patch_polytrig(); patch_polyvoice()
    print("Generating Soma smoke-test patches:")
    patch_soma_bursting(); patch_soma_chaos(); patch_soma_blips(); patch_soma_zmod()
    patch_soma_sync(); patch_soma_poly()
    print("Done.")
