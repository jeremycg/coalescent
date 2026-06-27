# Changelog

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
