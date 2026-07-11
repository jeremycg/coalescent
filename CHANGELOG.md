# Changelog

## 2.2.1 (unreleased)

**What you'll notice.** **Bunnies** now plays in tune (was ~37¢ sharp) and
**Haptik's DAMP** knob is musical across its whole travel. **GENDYN `DUR WID = 0`**
holds pitch exactly and **GENDYN LOCK** holds pitch even at high centre frequencies
/ large N (it used to drift flat near the sample floor). **Haptik FREEZE** no longer
clicks, and **GENDYN/Haptik now save their evolved or frozen sound** with the patch.
Everything else is internal robustness, performance, and concurrency work with no
audible effect.

**Heads-up for saved patches** (the *mapping* changed, not your stored settings):

- **Bunnies** plays ~37¢ lower than before — that's it now being in tune.
- **Haptik DAMP** is more resonant: a saved `DAMP = 0.5` is now ~31 Hz effective (was ~400 Hz).
- **GENDYN/Haptik** restore their evolved/frozen state; pre-2.2.1 or malformed patches fall back to a fresh seed.

### Audible fixes

- **Bunnies tuning**: `RATE_CAL` 7.49 → 7.33 (default was ~37¢ sharp); the stability test now asserts the constant tracks its measured value so it can't silently drift again.
- **Haptik DAMP retaper**: `DAMP_MAX_HZ` 800 → 250 with a quadratic taper, so the top ~85% of the knob no longer collapses to a click.
- **GENDYN `DUR WID = 0`**: removed a stray `+1`-sample floor that kept the walk moving (~41¢ flat); it now lands on the exact centre pitch.
- **GENDYN LOCK**: a per-cycle servo trims the residual rounding bias so the *measured* period tracks the target through the whole reachable range (best-effort only below the physical `fs/N` floor). It feeds back against the target each cycle was rendered against — not the live knob — and resets the trim when the target changes or drops below the floor, so a B DUR CTR / LOCK change lands within ~1 sample of pitch on the first cycle (exact within a few) instead of carrying a stale correction. The converged trim is saved so a near-floor sound recalls its pitch on load. FREQ reports the measured period.
- **Haptik Slow-mode FREEZE**: captures the interpolated shape currently being heard instead of jumping to the frame endpoint, so it holds (and un-freezes) without a click.
- **GENDYN & Haptik state**: breakpoints + walk velocities (GENDYN) and lattice + scan phase (Haptik) persist via `dataToJson`/`dataFromJson`.

### Correctness & robustness

- **Race-free saves**: the audio thread publishes a coherent state snapshot; `dataToJson()` reads that instead of racing `process()`.
- **Sanitize before fast `exp2`**: Axon/Soma/Haptik pitch CV is NaN/inf-guarded before `approxExp2` (its internal float→int shift was UB on a non-finite CV); Axon/Soma also cap sim-time so the substep count can't overflow.
- **Operon**: PERTURB sanitized before the integrator and the Hill argument flushed to 0 (a corrupt CV could reach the LUT index — an OOB read); LUT-rebuild counter saturates (was signed-overflow UB after hours).
- **Bunnies**: NaN backstop runs before the state clamp (the clamp was hiding NaN); the display's idle phase counter saturates (was UB after ~13.5 h).
- **Corrupt-patch bounds** on restored GENDYN/Haptik state (GENDYN amplitude ±1.5 — headroom over the dimensionless ±1, not a volt claim); GENDYN **Initialize** clears any in-flight cycle pulse.
- **GENDYN reseed** is atomic; **Haptik Fast/Slow** reconciles interpolation state; **Haptik EXT** sums all polyphonic channels (was channel 0 only).

### Performance

- **Operon**: the Hill LUT gates on whether `n` is *actually moving*, so a static/held HILL CV uses the LUT (~80× less `pow` vs the old cable-present gate); the `pStar` centre solve warm-starts Newton from the previous root (~8× fewer `pow` under DRIVE modulation, matches bisection to ~1e-6); and the LUT rebuild fills incrementally (256 entries/sample, spread over ~33 module frames) so no single `process()` call pays the whole ~60 µs table build, and a slow HILL LFO no longer thrashes the rebuild (it stays on direct `pow`).
- **Axon/Soma**: phase-trail rendering batched into 32 alpha bands — ~16× fewer NanoVG strokes (511 → 32 per voice) for a visually identical fade.
- **Haptik**: the lattice state/clamp pass is vectorized 4-wide, bit-identical to the scalar path (proven in `check-simd`).

### UI

- **Lock-free triple buffer** for every display snapshot (`coalescent::DisplaySnapshot<T>`), replacing an atomic-index double buffer that could tear mid-draw; each module defines only its payload.
- Fonts load per-frame (a cached `Font` can dangle across a Window recreation); panel captions moved to a dedicated child widget.
- **GENDYN scope** pairs each duration with the segment it plays and reflects the LOCK-normalized, error-diffused schedule (a clean-start approximation).
- **Haptik ring display** shows wave *shape* (DC removed, deviation-scaled) instead of reading as a circle.

### Tests, manifest & docs

- New **GENDYN** stability test (barriers + walk catches the DUR WID regression; LOCK servo near/below floor; target-change / LOCK-toggle / restore) and **Haptik FREEZE**-continuity test; `check` (SDK-free) and `check-simd` (Rack headers) are separate targets; offline WAV renderers use production constants and mirror the `subTau` cap.
- **Manifest**: `keywords` are spec-compliant space-separated strings; Operon/Bunnies tagged `Low-frequency oscillator` + `Clock generator`. Patch generators derive their Coalescent version from `plugin.json`.
- **Docs**: corrected CPU budgets (i5-9600K lower bounds; +6 oct via CV can saturate a core) incl. Operon's moving-HILL corner; manual install path (`plugins-<os>-<cpu>/`); GENDYN scope-sampling, `fs/N` ceiling, and state-recall wording; `Help → Open user folder`; poly-CV normalling; genetic/ecological taxonomy.

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
