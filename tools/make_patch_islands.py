#!/usr/bin/env python3
"""Generate four Core/Fundamental-only Islands demo patches.

Islands positional IDs (must match enum order in src/Islands.cpp):
  ParamId : 0 SIZE, 1 SELECT, 2 MUTATE, 3 MIGRATE, 4 GENERATIONS,
            5 FOUNDER, 6 RESET, 7 SELECT_ATT, 8 MUTATE_ATT, 9 MIGRATE_ATT
  InputId : 0 SIZE, 1 SELECT, 2 MUTATE, 3 MIGRATE, 4 GENERATIONS,
            5 STEP, 6 FOUNDER, 7 RESET
  OutputId: 0 I1, 1 I2, 2 I3, 3 I4, 4 MEAN, 5 HET,
            6 FIX_A, 7 FIX_B, 8 LOSS, 9 SWEEP

Fundamental v2 IDs used here:
  LFO     : FREQ_PARAM=2; TRI_OUTPUT=1; SQR_OUTPUT=3 (6 HP)
  8vert   : LEVEL_PARAM=channel; IN_INPUT=channel; OUT_OUTPUT=channel (8 HP)
  VCO     : FREQ_PARAM=2; PITCH_INPUT=0; SIN_OUTPUT=0 (10 HP)
  VCMixer : MASTER_PARAM=0, CH_LEVEL_PARAMS=1..4; CH_INPUTS=1..4; MIX_OUTPUT=0
"""

import glob
import io
import json
import os
import random
import shutil
import subprocess
import sys
import tarfile
import tempfile


random.seed(37)


def uid():
    return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)


ROOT = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(ROOT, "..", "plugin.json"), encoding="utf-8") as manifest:
    PLUGIN_VERSION = json.load(manifest)["version"]


def param(param_id, value):
    return {"id": param_id, "value": float(value)}


def islands(
    pos,
    size_log2=6.0,
    selection=0.0,
    mutation=0.35,
    migration=0.15,
    generations_log2=1.0,
    select_att=0.0,
    mutate_att=0.0,
    migrate_att=0.0,
):
    values = [
        size_log2,
        selection,
        mutation,
        migration,
        generations_log2,
        0.0,  # FOUNDER button
        0.0,  # RESET button
        select_att,
        mutate_att,
        migrate_att,
    ]
    return {
        "id": uid(),
        "plugin": "Coalescent",
        "model": "Islands",
        "version": PLUGIN_VERSION,
        "params": [param(i, value) for i, value in enumerate(values)],
        "pos": list(pos),
    }


def lfo(freq, pos):
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "LFO",
        "version": "2.6.4",
        "params": [param(2, freq)],
        "pos": list(pos),
    }


def eightvert(gains, pos):
    values = list(gains) + [0.0] * (8 - len(gains))
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "8vert",
        "version": "2.6.4",
        "params": [param(i, values[i]) for i in range(8)],
        "pos": list(pos),
    }


def vco(freq, pos):
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "VCO",
        "version": "2.6.4",
        "params": [param(2, freq)],
        "pos": list(pos),
    }


def vcmixer(levels, pos):
    values = [1.0] + list(levels) + [0.0] * (4 - len(levels))
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "VCMixer",
        "version": "2.6.4",
        "params": [param(i, values[i]) for i in range(5)],
        "pos": list(pos),
    }


def audio(pos):
    return {
        "id": uid(),
        "plugin": "Core",
        "model": "AudioInterface",
        "version": "2.6.6",
        "params": [],
        "pos": list(pos),
        "data": {
            "audio": {
                "driver": -1,
                "deviceName": "",
                "sampleRate": 44100.0,
                "blockSize": 256,
                "inputOffset": 0,
                "outputOffset": 0,
            },
            "dcFilter": True,
        },
    }


def cable(output_module, output_id, input_module, input_id, color_index):
    colors = ["#69c8ee", "#f0b45d", "#c77ae8", "#7bd477", "#e8e8e8"]
    return {
        "id": uid(),
        "outputModuleId": output_module,
        "outputId": output_id,
        "inputModuleId": input_module,
        "inputId": input_id,
        "color": colors[color_index % len(colors)],
        "inputPlugOrder": color_index,
        "outputPlugOrder": color_index,
    }


def write_patch(name, modules, cables, master_id):
    patch = {
        "version": "2.6.6",
        "zoom": 0.5,
        "gridOffset": [0.0, 0.0],
        "modules": modules,
        "cables": cables,
        "masterModuleId": master_id,
    }
    out_dir = os.path.join(ROOT, "..", "patches")
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, name)
    with tempfile.TemporaryDirectory() as tmp:
        json_path = os.path.join(tmp, "patch.json")
        with open(json_path, "w", encoding="utf-8") as handle:
            json.dump(patch, handle, indent=2)
        archive = io.BytesIO()
        with tarfile.open(fileobj=archive, mode="w:") as tar:
            tar.add(json_path, arcname="patch.json")
        result = subprocess.run(
            ["zstd", "-19", "-o", out_file, "-f"],
            input=archive.getvalue(),
            capture_output=True,
            check=False,
        )
        if result.returncode:
            print("zstd error:", result.stderr.decode(), file=sys.stderr)
            sys.exit(1)

    print(f"  {name}: {len(modules)} modules, {len(cables)} cables")
    for windows_dir in glob.glob("/mnt/c/Users/*/AppData/Local/Rack2/patches"):
        try:
            shutil.copy2(out_file, os.path.join(windows_dir, name))
        except OSError as error:
            print(f"    (skipped Windows copy: {error})", file=sys.stderr)


def add_four_voice_instrument(modules, cables, source, first_x):
    # 0..10 V island frequencies become a one-octave pitch span. Four base
    # pitches make the correlated movement audible as a slowly changing chord.
    attenuator = eightvert([0.1, 0.1, 0.1, 0.1], (first_x, 0))
    oscillators = [
        vco(-12.0, (first_x + 8, 0)),
        vco(-5.0, (first_x + 18, 0)),
        vco(0.0, (first_x + 28, 0)),
        vco(7.0, (first_x + 38, 0)),
    ]
    mixer = vcmixer([0.18, 0.18, 0.18, 0.18], (first_x + 48, 0))
    interface = audio((first_x + 55, 0))
    modules.extend([attenuator] + oscillators + [mixer, interface])

    for i, oscillator in enumerate(oscillators):
        cables.append(cable(source["id"], i, attenuator["id"], i, i))
        cables.append(cable(attenuator["id"], i, oscillator["id"], 0, i))
        cables.append(cable(oscillator["id"], 0, mixer["id"], i + 1, i))
    cables.append(cable(mixer["id"], 0, interface["id"], 0, 4))
    cables.append(cable(mixer["id"], 0, interface["id"], 1, 4))
    return interface["id"]


# Four uncoupled, unbiased populations make a chord whose voices wander within
# one octave. Symmetric mutation is high enough to keep a 64-copy demo moving
# after it touches a boundary, without choosing A or B.
def patch_neutral_lanes():
    source = islands(
        (0, 0), size_log2=6.0, selection=0.0, mutation=0.70,
        migration=0.0, generations_log2=2.0
    )
    modules = [source]
    cables = []
    master = add_four_voice_instrument(modules, cables, source, 16)
    write_patch("islands_1_neutral_lanes.vcv", modules, cables, master)


# A slow triangle moves migration between nearly independent and strongly
# coupled regimes. The four voices approach and separate without being copied.
def patch_migration():
    modulator = lfo(-5.0, (0, 0))
    source = islands(
        (6, 0), size_log2=6.0, selection=0.0, mutation=0.70,
        migration=0.5, generations_log2=3.0, migrate_att=0.8
    )
    modules = [modulator, source]
    cables = [cable(modulator["id"], 1, source["id"], 3, 4)]
    master = add_four_voice_instrument(modules, cables, source, 22)
    write_patch("islands_2_migration.vcv", modules, cables, master)


# Positive selection repeatedly carries A toward fixation. A slow square wave
# resets the experiment; SWEEP anticipates the extreme and LOSS marks exact global
# fixation. The event voltages are reduced before the stereo audio interface.
def patch_selection_sweep():
    reset_clock = lfo(-4.0, (0, 0))
    source = islands(
        (6, 0), size_log2=6.0, selection=0.12, mutation=0.0,
        migration=0.45, generations_log2=4.0
    )
    event_level = eightvert([0.25, 0.25], (22, 0))
    interface = audio((30, 0))
    modules = [reset_clock, source, event_level, interface]
    cables = [
        cable(reset_clock["id"], 3, source["id"], 7, 4),
        cable(source["id"], 9, event_level["id"], 0, 0),
        cable(source["id"], 8, event_level["id"], 1, 1),
        cable(event_level["id"], 0, interface["id"], 0, 0),
        cable(event_level["id"], 1, interface["id"], 1, 1),
    ]
    write_patch("islands_3_selection_sweep.vcv", modules, cables, interface["id"])


# At large SIZE the four lanes normally move in fine steps. Periodic FOUNDER
# triggers rotate an eight-copy bottleneck through them, making one voice take a
# conspicuous jump while the others continue smoothly.
def patch_founder():
    founder_clock = lfo(-3.0, (0, 0))
    source = islands(
        (6, 0), size_log2=10.0, selection=0.0, mutation=0.35,
        migration=0.02, generations_log2=2.0
    )
    modules = [founder_clock, source]
    cables = [cable(founder_clock["id"], 3, source["id"], 6, 4)]
    master = add_four_voice_instrument(modules, cables, source, 22)
    write_patch("islands_4_founder.vcv", modules, cables, master)


if __name__ == "__main__":
    print("Generating Islands demo patches:")
    patch_neutral_lanes()
    patch_migration()
    patch_selection_sweep()
    patch_founder()
    print("Done.")
