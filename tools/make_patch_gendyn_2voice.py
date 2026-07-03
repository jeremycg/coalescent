#!/usr/bin/env python3
"""Generate GENDY3_2voice.vcv — two slow, deep GENDYN oscillators.

Two voices at 60 Hz and 80 Hz, N=7, Gaussian distribution.
Scale is chosen so each voice evolves over roughly 30–40 seconds.
With the second-order walk (persistent steps), drift is ~8x faster than the
old first-order walk at equal scale, so the old diffusion formula gains a
0.35 factor (measured by simulation, constant across the usable range):
  scale = 0.35/sqrt(evolution_secs * update_rate)
  60 Hz, 35 sec target → scale = 0.35/sqrt(35*60) ≈ 0.0076
  80 Hz, 35 sec target → scale = 0.35/sqrt(35*80) ≈ 0.0066
Wide barriers (0.5) so drift is audible when it happens.
"""

import json, os, sys, io, tempfile, tarfile, subprocess, random

random.seed(99)

def uid():
    return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)

modules = []
cables  = []

# ── 2 GENDYN voices ──────────────────────────────────────────────────────────
# (center_hz, scale_amp, scale_dur)
VOICES = [
    (60,  0.004, 0.004),
    (80,  0.0035, 0.0035),
]

gendyn_ids = []
for i, (center_hz, s_amp, s_dur) in enumerate(VOICES):
    mid = uid()
    gendyn_ids.append(mid)
    modules.append({
        "id": mid,
        "plugin": "Coalescent",
        "model": "GENDYN",
        "version": "2.0.0",
        "params": [
            {"id": 0,  "value": 7.0},           # N (breakpoints)
            {"id": 1,  "value": s_amp},          # SCALE_AMP
            {"id": 2,  "value": s_dur},          # SCALE_DUR
            {"id": 3,  "value": 0.9},            # B_AMP (wide amplitude range)
            {"id": 4,  "value": float(center_hz)},# B_DUR_CENTER
            {"id": 5,  "value": 0.5},            # B_DUR_WIDTH (±50% pitch range)
            {"id": 6,  "value": 1.0},            # DISTRIBUTION: Gaussian
            {"id": 7,  "value": 0.0},            # SCALE_AMP_ATT
            {"id": 8,  "value": 0.0},            # SCALE_DUR_ATT
            {"id": 9,  "value": 0.0},            # B_AMP_ATT
            {"id": 10, "value": 0.0},            # B_DUR_ATT
            {"id": 11, "value": 0.3},            # PERSIST (glide persistence)
            {"id": 12, "value": 0.0},            # LOCK off — pitch wanders
        ],
        "pos": [i * 12, 0],
    })

# ── Fundamental::Mixer (2-channel mix → audio out) ───────────────────────────
mixer_id = uid()
modules.append({
    "id": mixer_id,
    "plugin": "Fundamental",
    "model": "Mixer",
    "version": "2.6.4",
    # Mixer's only param is one level knob; it sums all connected inputs.
    "params": [
        {"id": 0, "value": 0.25},
    ],
    "pos": [26, 0],
})

# ── Core::AudioInterface ──────────────────────────────────────────────────────
audio_id = uid()
modules.append({
    "id": audio_id,
    "plugin": "Core",
    "model": "AudioInterface",
    "version": "2.6.6",
    "params": [],
    "data": {
        "audio": {
            "driver": -1,
            "deviceName": "",
            "sampleRate": 44100.0,
            "blockSize": 256,
            "inputOffset": 0,
            "outputOffset": 0
        },
        "dcFilter": True
    },
    "pos": [40, 0],
})

# ── Cables ────────────────────────────────────────────────────────────────────
COLORS = ["#f3374b", "#3695ef"]

for i, gid in enumerate(gendyn_ids):
    cables.append({
        "id": uid(),
        "outputModuleId": gid,
        "outputId": 0,
        "inputModuleId": mixer_id,
        "inputId": i,
        "color": COLORS[i],
        "inputPlugOrder": i,
        "outputPlugOrder": i,
    })

# Mixer out → AudioInterface L and R
for ch in range(2):
    cables.append({
        "id": uid(),
        "outputModuleId": mixer_id,
        "outputId": 0,
        "inputModuleId": audio_id,
        "inputId": ch,
        "color": "#ffffff",
        "inputPlugOrder": ch,
        "outputPlugOrder": ch,
    })

# ── Write patch ───────────────────────────────────────────────────────────────
patch = {
    "version": "2.6.6",
    "zoom": 1.0,
    "gridOffset": [0.0, 0.0],
    "modules": modules,
    "cables": cables,
    "masterModuleId": audio_id,
}

out_dir  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "patches")
out_file = os.path.join(out_dir, "GENDY3_2voice.vcv")

with tempfile.TemporaryDirectory() as tmp:
    json_path = os.path.join(tmp, "patch.json")
    with open(json_path, "w") as f:
        json.dump(patch, f, indent=2)

    tar_buf = io.BytesIO()
    with tarfile.open(fileobj=tar_buf, mode="w:") as tf:
        tf.add(json_path, arcname="patch.json")

    result = subprocess.run(
        ["zstd", "-19", "-o", out_file, "-f"],
        input=tar_buf.getvalue(),
        capture_output=True,
    )
    if result.returncode != 0:
        print("zstd error:", result.stderr.decode(), file=sys.stderr)
        sys.exit(1)

print(f"Written: {out_file}")
print(f"  {len(modules)} modules, {len(cables)} cables")

import glob
_win = glob.glob("/mnt/c/Users/*/AppData/Local/Rack2/patches")
if _win:
    import shutil
    dst = os.path.join(_win[0], "GENDY3_2voice.vcv")
    shutil.copy2(out_file, dst)
    print(f"  Installed to Windows Rack: {dst}")
