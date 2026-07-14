#!/usr/bin/env python3
"""Generate four Core/Fundamental-only Finches demo patches.

Finches positional IDs (must match enum order in src/Finches.cpp):
  ParamId : 0 RATE, 1 MUTATE, 2 COMPETE, 3 NICHE, 4 MUTANT, 5 SEED,
            6 MUTATE_ATT, 7 COMPETE_ATT
  InputId : 0 RATE, 1 MUTATE, 2 COMPETE, 3 ENV, 4 SEED, 5 RESET
  OutputId: 0 MASS_L, 1 MASS_R, 2 PITCH_L, 3 PITCH_R, 4 SPREAD,
            5 SPLIT, 6 MERGE

Fundamental v2 IDs used here:
  LFO   : FREQ_PARAM=2; TRI_OUTPUT=1; SQR_OUTPUT=3 (6 HP)
  8vert : LEVEL_PARAM=0; IN_INPUT=0; OUT_OUTPUT=0 (8 HP; first channel)
  VCO   : FREQ_PARAM=2; PITCH_INPUT=0; SIN_OUTPUT=0 (10 HP)
  VCA-1 : LEVEL_PARAM=0; EXP_PARAM=1; CV_INPUT=0; IN_INPUT=1; OUT_OUTPUT=0
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


random.seed(23)


def uid():
    return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)


ROOT = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(ROOT, "..", "plugin.json"), encoding="utf-8") as manifest:
    PLUGIN_VERSION = json.load(manifest)["version"]


def param(param_id, value):
    return {"id": param_id, "value": float(value)}


def finches(values, pos, mutate_att=0.0, compete_att=0.0):
    return {
        "id": uid(),
        "plugin": "Coalescent",
        "model": "Finches",
        "version": PLUGIN_VERSION,
        "params": [param(i, value) for i, value in enumerate(values)]
        + [param(6, mutate_att), param(7, compete_att)],
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


def vco(pos):
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "VCO",
        "version": "2.6.4",
        "params": [param(2, 0.0)],
        "pos": list(pos),
    }


def eightvert(level, pos):
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "8vert",
        "version": "2.6.4",
        "params": [param(0, level)],
        "pos": list(pos),
    }


def vca(pos):
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "VCA-1",
        "version": "2.6.4",
        "params": [param(0, 1.0), param(1, 0.0)],
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
    colors = ["#f3374b", "#ffb437", "#00b56e", "#3695ef"]
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


def add_two_voice_instrument(modules, cables, finches_module, first_x):
    left_vco = vco((first_x, 0))
    right_vco = vco((first_x + 10, 0))
    left_vca = vca((first_x + 20, 0))
    right_vca = vca((first_x + 25, 0))
    interface = audio((first_x + 30, 0))
    modules.extend([left_vco, right_vco, left_vca, right_vca, interface])
    cables.extend(
        [
            cable(finches_module["id"], 2, left_vco["id"], 0, 0),
            cable(finches_module["id"], 3, right_vco["id"], 0, 1),
            cable(finches_module["id"], 0, left_vca["id"], 0, 0),
            cable(finches_module["id"], 1, right_vca["id"], 0, 1),
            cable(left_vco["id"], 0, left_vca["id"], 1, 2),
            cable(right_vco["id"], 0, right_vca["id"], 1, 3),
            cable(left_vca["id"], 0, interface["id"], 0, 0),
            cable(right_vca["id"], 0, interface["id"], 1, 1),
        ]
    )
    return interface["id"]


# One central phenotype becomes two audible VCO voices. The cluster masses open
# their VCAs while the cluster positions set pitch.
def patch_branching():
    fn = finches([1.0, 0.28, 0.78, 0.60, 0.72], (0, 0))
    modules = [fn]
    cables = []
    master = add_two_voice_instrument(modules, cables, fn, 14)
    write_patch("finches_1_branching.vcv", modules, cables, master)


# A very slow, attenuated LFO moves the environmental optimum. Once branched,
# both clusters follow it while retaining their ecological spacing.
def patch_moving_niche():
    modulator = lfo(-5.0, (0, 0))
    environment_depth = eightvert(0.273, (6, 0))
    fn = finches([1.0, 0.28, 0.90, 0.62, 0.72], (14, 0))
    modules = [modulator, environment_depth, fn]
    cables = [
        cable(modulator["id"], 1, environment_depth["id"], 0, 0),
        cable(environment_depth["id"], 0, fn["id"], 3, 1),
    ]
    master = add_two_voice_instrument(modules, cables, fn, 28)
    write_patch("finches_2_moving_niche.vcv", modules, cables, master)


# A slow square wave repeatedly challenges a branching population with a small,
# nearby mutant cohort. The aggregate field does not retain lineage labels.
def patch_mutant_invasion():
    seeder = lfo(-4.0, (0, 0))
    fn = finches([1.0, 0.20, 0.78, 0.58, 0.20], (6, 0))
    modules = [seeder, fn]
    cables = [cable(seeder["id"], 3, fn["id"], 4, 0)]
    master = add_two_voice_instrument(modules, cables, fn, 20)
    write_patch("finches_3_mutant_invasion.vcv", modules, cables, master)


# Slow COMPETE modulation repeatedly crosses the branch/merge region. The two
# accepted state transitions are auditioned as separate left/right trigger clicks.
def patch_events():
    modulator = lfo(-5.0, (0, 0))
    fn = finches(
        [1.0, 0.38, 0.52, 0.60, 0.72], (6, 0), compete_att=1.0
    )
    interface = audio((20, 0))
    modules = [modulator, fn, interface]
    cables = [
        cable(modulator["id"], 1, fn["id"], 2, 0),
        cable(fn["id"], 5, interface["id"], 0, 1),
        cable(fn["id"], 6, interface["id"], 1, 2),
    ]
    write_patch("finches_4_split_merge.vcv", modules, cables, interface["id"])


if __name__ == "__main__":
    print("Generating Finches demo patches:")
    patch_branching()
    patch_moving_niche()
    patch_mutant_invasion()
    patch_events()
    print("Done.")
