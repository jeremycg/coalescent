# Changelog

## 2.2.1 (unreleased)

- **GENDYN — DUR WID = 0 tuning fix**: a zero-width frequency barrier is
  documented to hold pitch fixed, but a stray `+1`-sample floor widened the window
  so the duration walk kept moving — playing ~41 cents flat at the default. Removed
  it; the walk now pins to `durCenter` and the error-diffused playback lands on the
  exact centre pitch (verified 0 cents across N and sample rates).
- **Axon/Soma — extreme-pitch step guard**: the per-sample sim-time is now capped
  (like Operon/Bunnies) so the integration step `h` stays ≤ `HSUB_MAX` and the
  substep count can't overflow the float→int conversion at very large V/OCT
  (undefined behaviour). Trade-off: with oversampling **Off**, pitch goes flat
  above ~+3.3 oct instead of integrating with an oversized step; the ×4 default is
  unaffected.
- **Haptik — Slow-mode FREEZE click**: engaging FREEZE in Slow mode jumped the
  readout from the interpolated inter-frame shape to the raw frame endpoint (a
  click). FREEZE now captures the shape currently being heard, so it holds
  seamlessly (and un-freezes cleanly).
- **Tests / tooling**: added a GENDYN stability test (barrier degeneracy, LOCK
  pitch, `fs/N` max-frequency floor, `reflect()` bounds — it would have caught the
  DUR WID bug) and wired it + the SIMD-equivalence check into `make check`; the
  SIMD check now covers Soma's 3-state HR path (0.000 cents scalar-vs-SIMD, incl.
  mixed-lane groups); the Soma/Axon offline WAV renderers use current production
  constants (Soma `RATE_CAL` 14.9→55.36, `MIN_SUB` 4→2) and no longer overstate
  output-stage parity.
- **Docs**: GENDYN scope is sampled at ~45 Hz, not low-pass filtered (claim
  corrected); documented GENDYN's `fs/N` frequency ceiling; Axon manual lists all
  nine patches (adds polytrig/polyvoice) and describes poly CV normalling correctly
  (monophonic broadcasts; partial-poly channels read 0 V); softened Axon's
  "one spike per trigger" wording; softened Haptik's display "race-free" comment.
- **Bunnies — tuning fix**: the LV default voicing was ~37 cents sharp
  (`RATE_CAL` 7.49, which put PITCH=0 at ~267 Hz instead of C4's 261.6 Hz). The
  constant was left stale after `LV_V0_RANGE` was retuned. Corrected to 7.33
  (verified −0.3 cents), and the stability test now **asserts** the constant
  tracks its measured value so this can't silently drift again. **Patch note:**
  saved Bunnies patches will play back ~37 cents lower (now in tune).
- **Bunnies — long-idle overflow fix**: the display's cycle-phase sample counter
  was an unbounded `int`; with no zero-crossing (e.g. RM at rest) it overflowed
  after ~13.5 h at 44.1 kHz (signed overflow = UB). Now saturates well below the
  limit.
- **Operon — HILL performance (CV path)**: the Hill lookup table was disabled
  whenever a cable was patched into HILL, forcing per-sample `pow()` in the
  derivative even when the CV was static or slow (and even when only DRIVE was the
  thing being modulated). The gate now keys on whether `n` is *actually moving*
  per sample: a settled `n` — from a knob at rest **or** a static/slow HILL CV —
  uses the LUT; only genuine audio-rate HILL modulation falls back to direct
  `pow()` (where a lagging LUT would detune). Measured ~80× less transcendental
  work in the patched-but-static case (≈2.1M → 27k `pow`/s at default pitch); the
  rate-limiter still backstops the knob-smoothing tail. No behaviour change.
- **Manifest / docs housekeeping**: `keywords` changed from arrays to the
  spec-compliant space-separated strings (VCV Manifest); README now points at
  `src/dsp/rk4.hpp`/`coalescent::rk4`; Operon scope comments corrected to ~25 Hz;
  Bunnies documents its residual DC; the Bunnies patch generator now checks `zstd`.
- **Haptik — DAMP retaper**: the damping knob is now musical across its whole
  travel. `DAMP_MAX_HZ` 800 → 250 with a quadratic knob taper (default 0.35 →
  0.30, ~44 ms decay). Previously everything above ~15 % of the knob collapsed to
  a click. **Patch note:** the stored value is unchanged but the *mapping* is, so
  saved Haptik patches will play back noticeably more resonant (a saved DAMP=0.5
  is now ~31 Hz effective, was ~400 Hz).
- **Haptik — ring display**: the scope now removes the ring's DC wander and scales
  to deviation-from-mean, so it shows the wave *shape* instead of reading as a
  circle at normal damping; brightness follows activity (dim near-circle at rest).

## 2.2.0

- **New module — Operon**: a three-phase oscillator on the Elowitz–Leibler
  repressilator (three genes repressing each other in a ring). Three protein
  outputs ~120° apart, three phase gates, DRIVE/HILL/DECAY/LEAK + PERTURB, and a
  three-protein time scope; usable as an audio voice, a three-phase LFO, or a clock.
  Math verified against BioModels BIOMD0000000012.
- **New module — Bunnies**: a predator–prey oscillator with two modes,
  Lotka–Volterra (conservative, amplitude set by a conserved-quantity servo) and
  Rosenzweig–MacArthur (self-correcting limit cycle). PREY/PRED outputs a quarter
  cycle apart, two peak gates, KICK prey-force, and a phase-orbit display with a
  bunny. Usable as an audio voice, a two-phase LFO, or a clock.
- The shared RK4 integrator moved from `src/neuron/integrator.hpp` to
  `src/dsp/rk4.hpp` (`coalescent::rk4`), now that a non-neuron module uses it.


## 2.1.0

- **Axon/Soma polyphony now runs 4 voices per SIMD group** (`simd::float_4`) —
  the whole oversampled chain (RK4 → spike detect → DC block → tanh → decimate)
  is vectorised. Measured ~4× on the integration chain at 16 voices; mono
  unchanged. Verified 0-cent equivalent to the scalar path (incl. mixed-lane
  groups); inactive lanes are masked so silent voices keep their state.
- **Community patches**: `patches/community/` ships two self-playing generative
  systems contributed by **fractalgee (Georg Carlson)** (CC BY 4.0 — thanks!),
  with a README covering credits, required plugins, and the press-RUN /
  select-your-audio-device gotchas.

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
