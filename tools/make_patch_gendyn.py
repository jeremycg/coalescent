#!/usr/bin/env python3
"""Generate GENDY3_16voice.vcv patch file for VCV Rack 2."""

import json, random, os
from patch_utils import windows_patch_directories, write_patch_archive
VER = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "plugin.json")))["version"]  # Coalescent version from the manifest

random.seed(42)

def uid():
    """Generate a positive 53-bit int as a Rack-style module/cable ID."""
    return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)

# ── Voice table ──────────────────────────────────────────────────────────────
# (N, center_hz)
VOICES = [
    ( 7,   80),   # 1
    (11,   95),   # 2
    (13,  110),   # 3
    (17,  130),   # 4
    (19,  155),   # 5
    (23,  185),   # 6
    (29,  220),   # 7
    (31,  260),   # 8
    ( 7,  310),   # 9
    (11,  370),   # 10
    (13,  440),   # 11
    (17,  520),   # 12
    (19,  620),   # 13
    (23,  740),   # 14
    (29,  880),   # 15
    (31, 1050),   # 16
]

CABLE_COLORS = [
    "#f3374b", "#ffb437", "#00b56e", "#3695ef",
    "#8b4ade", "#ffffff", "#c9b70e", "#e55728",
]

def cable_color(i):
    return CABLE_COLORS[i % len(CABLE_COLORS)]

modules = []
cables  = []

# ── GENDYN voices ─────────────────────────────────────────────────────────────
# Row 0: voices 1-8  (col 0,8,16,...,56)
# Row 1: voices 9-16 (col 0,8,16,...,56)
gendyn_ids = []
for i, (N, center_hz) in enumerate(VOICES):
    row = i // 8
    col = (i % 8) * 12
    mid = uid()
    gendyn_ids.append(mid)
    modules.append({
        "id": mid,
        "plugin": "Coalescent",
        "model": "GENDYN",
        "version": VER,
        "params": [
            {"id": 0,  "value": float(N)},          # N (breakpoints)
            {"id": 1,  "value": 0.004},            # SCALE_AMP (slower, stately morph)
            {"id": 2,  "value": 0.004},            # SCALE_DUR (slower, stately morph)
            {"id": 3,  "value": 0.8},              # B_AMP
            {"id": 4,  "value": float(center_hz)}, # B_DUR_CENTER
            {"id": 5,  "value": 0.4},              # B_DUR_WIDTH
            {"id": 6,  "value": 3.0},              # DISTRIBUTION (Logistic)
            {"id": 7,  "value": 0.0},              # SCALE_AMP_ATT
            {"id": 8,  "value": 0.0},              # SCALE_DUR_ATT
            {"id": 9,  "value": 0.0},              # B_AMP_ATT
            {"id": 10, "value": 0.0},              # B_DUR_ATT
            {"id": 11, "value": 0.3},              # PERSIST (glide persistence)
            {"id": 12, "value": 0.0},              # LOCK off — pitch wanders
        ],
        "pos": [col, row],
    })

# ── 4 × Fundamental::VCMixer  (groups of 4 voices, 9 HP each) ────────────────
# Inputs: 0=MIX_CV, 1-4=CH1-CH4 audio  Outputs: 0=MIX
# Params: 0=master, 1-4=channel levels
vcmixer_ids = []
for g in range(4):
    mid = uid()
    vcmixer_ids.append(mid)
    modules.append({
        "id": mid,
        "plugin": "Fundamental",
        "model": "VCMixer",
        "version": "2.6.4",
        "params": [
            {"id": 0, "value": 0.7},   # master
            {"id": 1, "value": 0.7},
            {"id": 2, "value": 0.7},
            {"id": 3, "value": 0.7},
            {"id": 4, "value": 0.7},
        ],
        "pos": [98 + g * 9, 0],
    })

# ── Fundamental::Mixer  (sums the 4 group outputs, 3 HP) ─────────────────────
# Inputs: 0-5=channels  Outputs: 0=MIX
# Params: 0=level (its only param — it sums all inputs through one level knob)
mixer_id = uid()
modules.append({
    "id": mixer_id,
    "plugin": "Fundamental",
    "model": "Mixer",
    "version": "2.6.4",
    "params": [
        {"id": 0, "value": 0.25},
    ],
    "pos": [135, 0],
})

# ── Core::AudioInterface ──────────────────────────────────────────────────────
# Inputs: 0=L, 1=R  (signals going OUT to audio hardware)
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
    "pos": [140, 0],
})

# ── Cables ────────────────────────────────────────────────────────────────────
ci = 0  # color index

# GENDYN AUDIO_OUT(0) → VCMixer channel inputs
for i, gid in enumerate(gendyn_ids):
    group  = i // 4       # which VCMixer (0-3)
    ch_in  = (i % 4) + 1  # ids 1-4 = CH_INPUTS audio (0=mix CV, 5-8=level CV) — verified from Fundamental src
    cables.append({
        "id": uid(),
        "outputModuleId": gid,
        "outputId": 0,   # AUDIO_OUTPUT
        "inputModuleId": vcmixer_ids[group],
        "inputId": ch_in,
        "color": cable_color(ci),
        "inputPlugOrder": ci,
        "outputPlugOrder": ci,
    })
    ci += 1

# VCMixer outputs(0) → Mixer channels 0-3
for g, vmid in enumerate(vcmixer_ids):
    cables.append({
        "id": uid(),
        "outputModuleId": vmid,
        "outputId": 0,
        "inputModuleId": mixer_id,
        "inputId": g,
        "color": cable_color(ci),
        "inputPlugOrder": ci,
        "outputPlugOrder": ci,
    })
    ci += 1

# Mixer out(0) → AudioInterface L and R
for ch in range(2):
    cables.append({
        "id": uid(),
        "outputModuleId": mixer_id,
        "outputId": 0,
        "inputModuleId": audio_id,
        "inputId": ch,
        "color": cable_color(ci),
        "inputPlugOrder": ci,
        "outputPlugOrder": ci,
    })
    ci += 1

# ── Assemble patch JSON ───────────────────────────────────────────────────────
patch = {
    "version": "2.6.6",
    "zoom": 0.5,
    "gridOffset": [0.0, 0.0],
    "modules": modules,
    "cables": cables,
    "masterModuleId": audio_id,
}

# ── Write as tar+zstd .vcv ────────────────────────────────────────────────────
out_dir  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "patches")
out_file = os.path.join(out_dir, "GENDY3_16voice.vcv")

write_patch_archive(out_file, patch)

print(f"Written: {out_file}")
print(f"  {len(modules)} modules, {len(cables)} cables")
print(f"  Size: {os.path.getsize(out_file)} bytes")

# Also copy to Windows Rack patches folder
_win = windows_patch_directories()
if _win:
    import shutil
    dst = os.path.join(_win[0], "GENDY3_16voice.vcv")
    try:
        shutil.copy2(out_file, dst)
        print(f"  Installed to Windows Rack: {dst}")
    except OSError as error:
        print(f"  (skipped Windows copy: {error})")
