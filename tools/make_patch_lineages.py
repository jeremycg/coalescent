#!/usr/bin/env python3
"""Generate three Core/Fundamental-only Lineages demo patches.

Lineages positional IDs (must match enum order in src/Lineages.cpp):
  ParamId : 0 RATE, 1 SAMPLES, 2 MUTATE, 3 STEP, 4 DIRECTION,
            5 LOOP, 6 NEW, 7 RESET
  InputId : 0 RATE, 1 TIME, 2 NEW, 3 RESET
  OutputId: 0 TRAITS (poly), 1 NODE, 2 MUTATION, 3 MRCA,
            4 LINEAGES, 5 DIVERSITY

Fundamental v2 IDs used here:
  LFO   : OFFSET_PARAM=0 (1=unipolar); FREQ_PARAM=2; TRI_OUTPUT=1
  8vert : LEVEL_PARAMS/INPUTS/OUTPUTS=0..7
  VCO   : FREQ_PARAM=2; PITCH_INPUT=0; SIN_OUTPUT=0
  ADSR  : A/D/S/R_PARAMS=0..3; GATE_INPUT=4; ENVELOPE_OUTPUT=0
  VCA-1 : LEVEL_PARAM=0; EXP_PARAM=1; CV_INPUT=0; IN_INPUT=1;
          OUT_OUTPUT=0
  Sum   : LEVEL_PARAM=0; POLY_INPUT=0; MONO_OUTPUT=0

SAMPLES and MUTATE remain at their constructor defaults in every patch. The
approved Lineages constructor therefore supplies each newly loaded demo with a
valid, clean tree; this generator deliberately does not duplicate or freeze the
genealogy mathematics in Python.
"""

import json
import os
import random
import shutil
import sys

from patch_utils import windows_patch_directories, write_patch_archive


random.seed(73)


def uid():
    return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)


ROOT = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(ROOT, "..", "plugin.json"), encoding="utf-8") as manifest:
    PLUGIN_VERSION = json.load(manifest)["version"]


def param(param_id, value):
    return {"id": param_id, "value": float(value)}


def lineages(pos, rate=-3.0, step=3.0, direction=1.0, loop=1.0):
    values = [
        rate,
        8.0,   # SAMPLES: match the generated constructor tree
        0.4,   # MUTATE: match the generated constructor tree
        step,
        direction,
        loop,
        0.0,   # NEW button
        0.0,   # RESET button
    ]
    return {
        "id": uid(),
        "plugin": "Coalescent",
        "model": "Lineages",
        "version": PLUGIN_VERSION,
        "params": [param(i, value) for i, value in enumerate(values)],
        "pos": list(pos),
    }


def lfo(frequency, pos, unipolar=False):
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "LFO",
        "version": "2.6.4",
        "params": [param(0, 1.0 if unipolar else 0.0), param(2, frequency)],
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


def vco(frequency, pos):
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "VCO",
        "version": "2.6.4",
        "params": [param(2, frequency)],
        "pos": list(pos),
    }


def adsr(attack, decay, sustain, release, pos):
    values = [attack, decay, sustain, release]
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "ADSR",
        "version": "2.6.4",
        "params": [param(i, value) for i, value in enumerate(values)],
        "pos": list(pos),
    }


def vca(level, pos):
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "VCA-1",
        "version": "2.6.4",
        "params": [param(0, level), param(1, 0.0)],
        "pos": list(pos),
    }


def summ(level, pos):
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "Sum",
        "version": "2.6.4",
        "params": [param(0, level)],
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
    write_patch_archive(out_file, patch)

    print(f"  {name}: {len(modules)} modules, {len(cables)} cables")
    for windows_dir in windows_patch_directories():
        try:
            shutil.copy2(out_file, os.path.join(windows_dir, name))
        except OSError as error:
            print(f"    (skipped Windows copy: {error})", file=sys.stderr)


def add_poly_voice(modules, cables, source, first_x):
    """Turn Lineages' fixed-polyphony trait vector into a quiet sine cluster."""
    oscillator = vco(-12.0, (first_x, 0))
    voice_sum = summ(0.12, (first_x + 9, 0))
    interface = audio((first_x + 12, 0))
    modules.extend([oscillator, voice_sum, interface])
    cables.extend(
        [
            cable(source["id"], 0, oscillator["id"], 0, 0),
            cable(oscillator["id"], 0, voice_sum["id"], 0, 1),
            cable(voice_sum["id"], 0, interface["id"], 0, 2),
            cable(voice_sum["id"], 0, interface["id"], 1, 2),
        ]
    )
    return interface["id"]


# Starting at the MRCA, eight equal sine voices descend slowly through the tree.
# Signed branch mutations turn the unison into a related pitch cluster; Sum keeps
# the worst-case in-phase level conservative without collapsing source polyphony.
def patch_polyphonic_descent():
    source = lineages((0, 0), rate=-3.0, step=3.0, direction=1.0, loop=1.0)
    modules = [source]
    cables = []
    master = add_poly_voice(modules, cables, source, 16)
    write_patch("lineages_1_polyphonic_descent.vcv", modules, cables, master)


# NODE and MUTATION crossings excite separate envelope/VCA voices in the left and
# right channels. LINEAGES bends the low node voice while DIVERSITY bends the
# higher mutation voice. A normalled 10 V lane gives RESET one startup edge after
# Rack applies ANCESTRY, moving the existing tree to its present-day source.
def patch_ancestral_rhythm():
    source = lineages((0, 0), rate=-2.0, step=2.0, direction=0.0, loop=1.0)
    modulation = eightvert([0.015, 0.02, 1.0], (16, 0))
    node_oscillator = vco(-24.0, (24, 0))
    node_envelope = adsr(0.0, 0.08, 0.0, 0.12, (33, 0))
    node_vca = vca(0.25, (42, 0))
    mutation_oscillator = vco(-5.0, (45, 0))
    mutation_envelope = adsr(0.0, 0.18, 0.0, 0.30, (54, 0))
    mutation_vca = vca(0.20, (63, 0))
    interface = audio((66, 0))
    modules = [
        source,
        modulation,
        node_oscillator,
        node_envelope,
        node_vca,
        mutation_oscillator,
        mutation_envelope,
        mutation_vca,
        interface,
    ]
    cables = [
        cable(source["id"], 1, node_envelope["id"], 4, 0),
        cable(source["id"], 2, mutation_envelope["id"], 4, 1),
        cable(source["id"], 4, modulation["id"], 0, 2),
        cable(source["id"], 5, modulation["id"], 1, 3),
        cable(modulation["id"], 2, source["id"], 3, 4),
        cable(modulation["id"], 0, node_oscillator["id"], 0, 2),
        cable(modulation["id"], 1, mutation_oscillator["id"], 0, 3),
        cable(node_oscillator["id"], 0, node_vca["id"], 1, 0),
        cable(node_envelope["id"], 0, node_vca["id"], 0, 0),
        cable(mutation_oscillator["id"], 0, mutation_vca["id"], 1, 1),
        cable(mutation_envelope["id"], 0, mutation_vca["id"], 0, 1),
        cable(node_vca["id"], 0, interface["id"], 0, 0),
        cable(mutation_vca["id"], 0, interface["id"], 1, 1),
    ]
    write_patch(
        "lineages_2_ancestral_rhythm.vcv", modules, cables, interface["id"]
    )


# A unipolar 0..10 V triangle owns TIME absolutely. Its rising and falling halves
# traverse the same stored tree in opposite directions, audibly separating and
# rejoining the fixed eight-channel voice without NEW or transport regeneration.
def patch_cv_scrub():
    scrubber = lfo(-5.0, (0, 0), unipolar=True)
    source = lineages((9, 0), rate=0.0, step=3.0, direction=1.0, loop=1.0)
    modules = [scrubber, source]
    cables = [cable(scrubber["id"], 1, source["id"], 1, 0)]
    master = add_poly_voice(modules, cables, source, 25)
    write_patch("lineages_3_cv_scrub.vcv", modules, cables, master)


if __name__ == "__main__":
    print("Generating Lineages demo patches:")
    patch_polyphonic_descent()
    patch_ancestral_rhythm()
    patch_cv_scrub()
    print("Done.")
