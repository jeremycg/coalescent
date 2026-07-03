# Changelog

## 2.0.3

Performance and discoverability pass, prompted by user feedback (issue #2).

**Performance** (Axon/Soma especially — no audible change):
- Lower the RK4 substep floor `MIN_SUB` 4 → 2: the adaptive stepper was pinned at
  the floor across most of the normal range (worst with oversampling), doing ~2×
  the needed integration. Profiled at 0-cent / −74 dB, ~50% less RK4 CPU.
- Replace `std::tanh` in the output soft-clip with a Padé[7/6] approximation
  (~14× faster, −64 dB error, inaudible).
- Use `dsp::approxExp2_taylor5` for the per-sample V/OCT pitch (and Soma BURST /
  Haptik evolution-rate) conversions; compute Haptik's damping coefficient with
  the fast exp too.

**Defaults & discoverability:**
- **GENDYN** now seeds a clean **sine** (not random noise), so it boots as a
  pitched tone that the stochastic walk evolves. Right-click **Initial waveform**
  (Sine/Triangle/Saw/Square/Random) + **Re-seed waveform**; centre default → C4;
  internal ~5 Hz DC blocker on OUT.
- **Soma** default is now **tonic spiking** (a clear C4 tone) rather than bursting;
  turn BURST down for bursting/chaos.
- **Haptik** default **sustains** (EXCITE = drive) instead of decaying after one
  pluck.

**Anti-aliasing:** the Axon/Soma menu adds a **×2** option (Off / ×2 / ×4 / ×8,
default ×4) for lower CPU. Saved patches store the actual factor now (with a
legacy remap), so pre-2.0.3 patches keep their setting despite the menu change.

**Fixes:**
- **GENDY3_cluster / GENDY3_16voice demo patches were silent** — their GENDYN
  audio was wired into the VCMixer's channel-*CV* inputs (a regression introduced
  by an old "wiring fix"); rewired to the channel audio inputs and all GENDYN
  demo morph rates retuned slower.
- **Haptik**: a NaN on V/OCT could permanently stick the scan phase and read out
  of bounds (potential crash); guarded.
- **GENDYN**: Initialize now always re-seeds (and resets the seed shape); the
  FREQ CV no longer over-reports pitch under LOCK at extreme centre-freq × N;
  re-seeding starts exactly on the chosen waveform (no first-cycle click).
- Refreshed panel screenshots; new `coalescent_gallery.vcv` showcase patch (all
  four modules, poly neurons, croppable layout).

## 2.0.2

More panel legibility (continuing the report in issue #1):

- Control labels now use Rack's **Nunito Bold** system font — bold reads far
  better at small sizes than the previous regular weight; nothing is bundled.
- Labels unified to a single larger size, near-white on the dark panels.
- Knob labels moved **above** their controls, consistent across all four modules
  (matching the CV-strip labels), so no label/label collisions.
- Centered the knob and CV rows on the 12 HP panels; trimmed the GENDYN/Haptik
  scopes ~2 mm to make room.
- Added `patches/coalescent_gallery.vcv` — all four modules side by side for a
  quick visual check.

No DSP or behaviour change.

## 2.0.1

Panel legibility (reported by a user): the control labels were small and
low-contrast (~4:1) against the dark panels. Brightened the label colour to
`#c8c8e4` (~11.5:1) and scaled all panel labels ~22% larger, across all four
modules. No DSP or behaviour change.

## 2.0.0

Initial release of **Coalescent** — the merge of four previously separate plugins
into one. The modules (GENDYN, Haptik, Axon, Soma) are unchanged in behaviour;
this release is the unification:

- Single plugin slug/brand `Coalescent`, one manifest, one Makefile.
- Four modules under one brand: **GENDYN**, **Haptik**, and the **Neuron** pair
  **Axon** (FitzHugh–Nagumo) and **Soma** (Hindmarsh–Rose), grouped by a shared
  name prefix and accent/panel language (and both tagged `Polyphonic`).
- Axon and Soma share an RK4 integrator (`src/neuron/integrator.hpp`).
- Demo patches regenerated against the `Coalescent` plugin slug.
- Per-module documentation under `docs/`; new top-level README.

Module behaviour carried forward from the prior standalone plugins:

- **Axon / Soma** — polyphony (up to 16 voices), hard SYNC input, and a right-click
  anti-aliasing (oversampling) option.
- **GENDYN** — Xenakis GENDY3 dynamic stochastic synthesis with pitch LOCK.
- **Haptik** — scanned synthesis with Fast/Slow lattice modes and FREEZE.
