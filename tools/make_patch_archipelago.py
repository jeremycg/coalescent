#!/usr/bin/env python3
"""Generate four Core/Fundamental-only Archipelago demo patches.

Archipelago positional IDs (must match enum order in src/Archipelago.cpp):
  ParamId : 0 RATE, 1 SELECT, 2 MUTATE, 3 MIGRATE, 4 GRADIENT,
            5 BARRIER, 6 CLIMATE, 7 RESET, 8 SELECT_ATT, 9 MUTATE_ATT,
            10 MIGRATE_ATT, 11 GRADIENT_ATT, 12 BARRIER_ATT, 13 TOPOLOGY
  InputId : 0 RATE, 1 SELECT, 2 MUTATE, 3 MIGRATE, 4 GRADIENT,
            5 BARRIER, 6 CLIMATE, 7 RESET
  OutputId: 0 TRAIT (poly 8), 1 MASS (poly 8), 2 MEAN, 3 DIFF, 4 FLUX,
            5 COLONIZE, 6 EXTINCT

Fundamental v2 IDs used here:
  LFO     : OFFSET_PARAM=0 (0=bipolar); FREQ_PARAM=2; TRI_OUTPUT=1 (6 HP)
  VCO     : FREQ_PARAM=2; PITCH_INPUT=0; SIN_OUTPUT=0 (10 HP, poly)
  VCA-1   : LEVEL_PARAM=0; EXP_PARAM=1; CV_INPUT=0; IN_INPUT=1;
            OUT_OUTPUT=0 (5 HP, poly)
  Sum     : LEVEL_PARAM=0; POLY_INPUT=0; MONO_OUTPUT=0 (2 HP)
  8vert   : LEVEL_PARAMS/INPUTS/OUTPUTS=0..7 (8 HP)
  VCMixer : MASTER_PARAM=0; CH_LEVEL_PARAMS=1..4; CH_INPUTS=1..4;
            MIX_OUTPUT=0 (9 HP)

Archipelago parameter storage:
  RATE      : -8..+4 octaves (0 = 4 evolutionary tau/s)
  SELECT    : unit 0..1, physical strength 8*u^2 / tau
  MUTATE    : unit 0..1, true zero then 1e-5*300^u / tau
  MIGRATE   : unit 0..1, true zero then 0.002*1000^u / tau
  GRADIENT  : -0.85..+0.85 trait
  BARRIER   : 0..1; central-link permeability is (1-b)^2
  CLIMATE   : -0.85..+0.85 trait
  TOPOLOGY  : 0 row, 1 ring
"""

import glob
import io
import json
import math
import os
import random
import shutil
import subprocess
import sys
import tarfile


random.seed(61)


def uid():
    return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)


ROOT = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(ROOT, "..", "plugin.json"), encoding="utf-8") as manifest:
    PLUGIN_VERSION = json.load(manifest)["version"]


def param(param_id, value):
    return {"id": param_id, "value": float(value)}


def adapted_state(gradient, climate):
    """Return the same deterministic field that RESET builds for these controls."""
    bins = 32
    sigma = 0.065
    density = []
    reported_traits = []
    global_moment = 0.0
    total_mass = 0.0
    for habitat in range(8):
        position = -1.0 + 2.0 * habitat / 7.0
        optimum = max(-1.0, min(climate + gradient * position, 1.0))
        traits = [-1.0 + (trait_bin + 0.5) * 2.0 / bins
                  for trait_bin in range(bins)]
        population = [math.exp(-0.5 * ((trait - optimum) / sigma) ** 2)
                      for trait in traits]
        scale = 0.75 / sum(population)
        population = [value * scale for value in population]
        moment = sum(value * trait for value, trait in zip(population, traits))
        density.extend(population)
        reported_traits.append(moment / 0.75)
        global_moment += moment
        total_mass += 0.75
    return {
        "archipelagoVersion": 1,
        "density": density,
        "reportedTrait": reported_traits,
        "reportedGlobalMean": global_moment / total_mass,
        "occupiedMask": 0xff,
    }


def archipelago(
    pos,
    rate=0.0,
    selection=0.60,
    mutation=0.45,
    migration=0.40,
    gradient=0.55,
    barrier=0.0,
    climate=0.0,
    select_att=0.0,
    mutate_att=0.0,
    migrate_att=0.0,
    gradient_att=0.0,
    barrier_att=0.0,
    topology=0.0,
):
    values = [
        rate,
        selection,
        mutation,
        migration,
        gradient,
        barrier,
        climate,
        0.0,  # RESET button
        select_att,
        mutate_att,
        migrate_att,
        gradient_att,
        barrier_att,
        topology,
    ]
    return {
        "id": uid(),
        "plugin": "Coalescent",
        "model": "Archipelago",
        "version": PLUGIN_VERSION,
        "params": [param(i, value) for i, value in enumerate(values)],
        "pos": list(pos),
        "data": adapted_state(gradient, climate),
    }


def lfo(frequency, pos):
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "LFO",
        "version": "2.6.4",
        "params": [param(0, 0.0), param(2, frequency)],
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


def vca(pos):
    return {
        "id": uid(),
        "plugin": "Fundamental",
        "model": "VCA-1",
        "version": "2.6.4",
        "params": [param(0, 1.0), param(1, 0.0)],
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
    patch_json = json.dumps(patch, indent=2).encode("utf-8")
    archive = io.BytesIO()
    with tarfile.open(fileobj=archive, mode="w:") as tar:
        info = tarfile.TarInfo("patch.json")
        info.size = len(patch_json)
        info.mtime = 0
        info.mode = 0o644
        tar.addfile(info, io.BytesIO(patch_json))
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


def add_poly_local_adaptation_instrument(modules, cables, source, first_x):
    """Map the eight local populations literally to one polyphonic instrument."""
    oscillator = vco(-12.0, (first_x, 0))
    amplifier = vca((first_x + 10, 0))
    voice_sum = summ(0.12, (first_x + 15, 0))
    interface = audio((first_x + 17, 0))
    modules.extend([oscillator, amplifier, voice_sum, interface])
    cables.extend(
        [
            cable(source["id"], 0, oscillator["id"], 0, 0),
            cable(source["id"], 1, amplifier["id"], 0, 1),
            cable(oscillator["id"], 0, amplifier["id"], 1, 2),
            cable(amplifier["id"], 0, voice_sum["id"], 0, 3),
            cable(voice_sum["id"], 0, interface["id"], 0, 4),
            cable(voice_sum["id"], 0, interface["id"], 1, 4),
        ]
    )
    return interface["id"]


# Eight local traits become oscillator pitches while the corresponding local
# population masses open eight VCA lanes. Sum reduces the polyphonic result to a
# conservative mono level without erasing the local voices inside the patch.
def patch_local_adaptation():
    source = archipelago(
        (0, 0), rate=0.5, selection=0.68, mutation=0.42,
        migration=0.38, gradient=0.62, barrier=0.0, climate=0.0,
        topology=0.0,
    )
    modules = [source]
    cables = []
    master = add_poly_local_adaptation_instrument(modules, cables, source, 18)
    write_patch("archipelago_1_local_adaptation.vcv", modules, cables, master)


# A slow triangle moves migration around a moderate baseline. The audible lanes
# alternately retain local pitch identities and move as a more coherent group.
def patch_migration_modulation():
    modulator = lfo(-5.0, (0, 0))
    source = archipelago(
        (6, 0), rate=1.0, selection=0.62, mutation=0.42,
        migration=0.50, gradient=0.58, barrier=0.0, climate=0.0,
        migrate_att=0.75, topology=0.0,
    )
    modules = [modulator, source]
    cables = [cable(modulator["id"], 1, source["id"], 3, 0)]
    master = add_poly_local_adaptation_instrument(modules, cables, source, 24)
    write_patch("archipelago_2_migration.vcv", modules, cables, master)


# A slow climate sweep moves the habitable range through the eight local
# populations. MASS shapes the poly voice continuously; COLONIZE and EXTINCT are
# mixed quietly as distinct event clicks around the summed population sound.
def patch_climate_range_shift():
    climate = lfo(-5.0, (0, 0))
    source = archipelago(
        (6, 0), rate=1.0, selection=0.78, mutation=0.40,
        migration=0.42, gradient=0.82, barrier=0.0, climate=0.0,
        topology=0.0,
    )
    oscillator = vco(-12.0, (24, 0))
    amplifier = vca((34, 0))
    voice_sum = summ(0.10, (39, 0))
    event_level = eightvert([0.08, 0.08], (41, 0))
    mixer = vcmixer([0.75, 0.45, 0.45], (49, 0))
    interface = audio((58, 0))
    modules = [climate, source, oscillator, amplifier, voice_sum,
               event_level, mixer, interface]
    cables = [
        cable(climate["id"], 1, source["id"], 6, 0),
        cable(source["id"], 0, oscillator["id"], 0, 0),
        cable(source["id"], 1, amplifier["id"], 0, 1),
        cable(oscillator["id"], 0, amplifier["id"], 1, 2),
        cable(amplifier["id"], 0, voice_sum["id"], 0, 3),
        cable(source["id"], 5, event_level["id"], 0, 0),
        cable(source["id"], 6, event_level["id"], 1, 1),
        cable(voice_sum["id"], 0, mixer["id"], 1, 2),
        cable(event_level["id"], 0, mixer["id"], 2, 0),
        cable(event_level["id"], 1, mixer["id"], 3, 1),
        cable(mixer["id"], 0, interface["id"], 0, 4),
        cable(mixer["id"], 0, interface["id"], 1, 4),
    ]
    write_patch("archipelago_3_climate_shift.vcv", modules, cables, interface["id"])


# High migration would normally couple the range strongly; a nearly closed
# central barrier lets the two halves retain different locally adapted traits.
def patch_central_barrier():
    source = archipelago(
        (0, 0), rate=1.0, selection=0.72, mutation=0.38,
        migration=0.75, gradient=0.72, barrier=0.88, climate=0.0,
        topology=0.0,
    )
    modules = [source]
    cables = []
    master = add_poly_local_adaptation_instrument(modules, cables, source, 18)
    write_patch("archipelago_4_barrier.vcv", modules, cables, master)


if __name__ == "__main__":
    print("Generating Archipelago demo patches:")
    patch_local_adaptation()
    patch_migration_modulation()
    patch_climate_range_shift()
    patch_central_barrier()
    print("Done.")
