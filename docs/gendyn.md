# GENDYN

A VCV Rack 2 implementation of Iannis Xenakis's dynamic stochastic synthesis
algorithm (GENDY3, 1991). Part of the **Coalescent** plugin's *Fluctuations*
series — see the [main README](../README.md).

![GENDYN module screenshot](img/gendyn.png)

> Xenakis's original: [GENDY3 on YouTube](https://www.youtube.com/watch?v=5qS5lqbx9H0)

---

## Algorithm

GENDYN implements a piecewise-linear oscillator whose waveform shape
evolves stochastically each cycle. Each of the N breakpoints has an
amplitude and a duration (in samples) governed by a **second-order
random walk** (Serra 1993, eqs. 2 & 4; Hoffmann 2023):

    step_{i,j+1} = MIR(step_{i,j} + f(z), −1, +1)
    y_{i,j+1}    = MIR(y_{i,j} + scale·barrier·step_{i,j+1}, barrier_min, barrier_max)

A fresh draw from the chosen distribution nudges a persistent,
normalized step variable (the primary walk), mirrored into [−1, 1];
the step then moves the breakpoint (the secondary walk) scaled by
SCALE × barrier, mirrored into the barrier range. The spread of the
draws f(z) is set by the PERSIST knob and determines how many cycles
the step keeps its direction. Serra's eq. 4 mirror structure — one
mirror pair on the step values, one on the positions — corresponds to
the primary and secondary barriers respectively.

This second-order structure is what distinguishes GENDY3 from S.709,
where draws move breakpoints directly (Hoffmann 2023). Because steps
are correlated cycle to cycle, duration changes accumulate in one
direction for many cycles, producing GENDY3's characteristic directed
glissandi; the step limits cap the maximum glide rate (SCALE DUR) and
the rate of timbral change (SCALE AMP).

**Continuity condition** (Serra eq. 1): every segment starts from the
value the output just reached, so the segment that crosses the cycle
boundary interpolates from the previous cycle's last breakpoint to the
freshly walked breakpoint 0 — the waveform is seamless by construction
and all N segments carry the stochastic walk.

Frequency is emergent — it equals `sampleRate / sum(dur[i])` across
all N breakpoints.

## Controls

| Control | Description |
|---------|-------------|
| N | Number of breakpoints (integer, 2–64) |
| SCALE AMP | Step scale for amplitude random walk (0–1) |
| SCALE DUR | Step scale for duration random walk (0–1) |
| B AMP | Amplitude barrier half-width (0=frozen, 1=full ±5V range) |
| B DUR CTR | Center frequency for duration barriers (20–5000 Hz) |
| B DUR WID | Duration barrier half-width around center (0=fixed pitch, 1=wide) |
| DIST | Distribution: 0=Cauchy, 1=Gaussian, 2=Uniform, 3=Logistic (default) |
| PERSIST | Glide persistence: how many cycles a step keeps its direction. 0% ≈ uncorrelated jitter (first-order / SC Gendy feel), 30% (default) ≈ 16 cycles, 100% = very long steady glides |
| LOCK | Pitch lock: normalizes durations each cycle so pitch holds exactly at B DUR CTR while the waveform keeps evolving (SC `Gendy3` behaviour) |

## CV Inputs

All CV inputs are ±5V with attenuverter knobs (±5V × attenuverter × 0.1 = ±0.5 modulation depth).

| Input | Target |
|-------|--------|
| SCALE AMP CV | Amplitude step scale |
| SCALE DUR CV | Duration step scale |
| B AMP CV | Amplitude barrier |
| B DUR CV | Duration barrier width |

## Outputs

| Output | Description |
|--------|-------------|
| OUT | Audio output, ±5V |
| TRIG | 10V / 1ms trigger on each complete waveform cycle |
| FREQ | Current frequency as 1V/oct CV (0V = C4 = 261.626 Hz) |

## Display

The screen up top draws the live waveform polygon: the N breakpoints as a
piecewise-linear curve over one cycle (x = cumulative duration, y = amplitude).
As the two random walks evolve, the vertices drift — vertically (the amplitude
walk) and horizontally (the duration walk) — so the polygon slowly morphs in
place. The faint horizontal band marks the ±B AMP barriers the amplitude walk
reflects within. The display low-passes the breakpoints (~0.4 s) so the picture
eases at a watchable pace — the walk actually updates every cycle (audio rate),
which is ~an order of magnitude too fast to follow, so the smoothing removes the
per-cycle shimmer while leaving the slow morph (and the audio) untouched. The
per-sample playhead and the individual random draws are deliberately not shown.

## Tuning notes

- **SCALE AMP / SCALE DUR** set the *maximum rate of change*: the
  per-cycle move is capped at SCALE × barrier. Steps are correlated
  between cycles, so even small values produce audible directed
  motion — glissandi for SCALE DUR, timbral drift for SCALE AMP.
  For a target evolution time, `scale ≈ 0.35 / sqrt(seconds × Hz)`.
- **PERSIST** trades jitter against glide at constant step size: low
  settings decorrelate the steps every cycle (rough, noisy, close to
  the classic SuperCollider Gendy texture), high settings hold a
  direction for hundreds of cycles (steady glissandi).
- **LOCK** rescales each cycle's durations to sum to exactly
  `sampleRate / B DUR CTR`: pitch is constant, but the *relative*
  durations keep walking, so the waveform — and therefore timbre —
  still evolves. This matches SuperCollider's `Gendy3.ar`. The wider
  B DUR WID, the more timbral duration variety under the locked pitch.
- **B DUR WID = 0** also fixes pitch (all durations equal) but kills
  duration-walk timbre with it; LOCK with a wide B DUR WID is how
  Xenakis got the "beautiful clear tones" — stable pitch, living
  waveform.
- Logistic distribution (DIST=3) is the closest match to Xenakis's original and is the default.
- **SuperCollider Gendy mode:** PERSIST 0%, SCALE ≈ 0.1, DIST Cauchy,
  B DUR CTR ≈ 520 Hz with B DUR WID ≈ 0.2 lands close to `Gendy1.ar()`
  at its defaults.
- **Aliasing is intended character, not a defect.** The output is a piecewise-
  linear waveform with randomised breakpoints — it is *not* a band-limited
  oscillator. At high B DUR CTR or extreme duration/amplitude settings the sharp
  corners produce audible aliasing. That raw digital edge is part of the GENDY3
  sound; tame it with a downstream low-pass if you want it cleaner.

## Compositional technique

Xenakis composed GENDY3 by varying barrier widths per voice per section. Wide barriers allow chaotic drift; narrow barriers freeze or constrain the walk. In a multi-voice patch, slow modulation of B AMP and B DUR CV inputs shapes large-scale form.

See the `patches/` folder for reference patches:
- `GENDY3_2voice.vcv` — two slow low-register voices (60 Hz, 80 Hz)
- `GENDY3_16voice.vcv` — 16-voice spread across the spectrum
- `GENDY3_cluster.vcv` — Xenakis's actual 14-voice pitch cluster from the score (Hoffmann 2022, Table 1), voiced like SuperCollider's `Gendy3` at its defaults: Cauchy, uncorrelated steps, pitch LOCKed per voice with the waveform free to evolve

Patch files can be regenerated from the scripts in `tools/`:

```bash
python3 tools/make_patch_gendyn_2voice.py
python3 tools/make_patch_gendyn.py
python3 tools/make_patch_gendyn_cluster.py
```

## References

- Serra, M.-H. (1993). Stochastic Composition and Stochastic Timbre:
  GENDY3 by Iannis Xenakis. *Perspectives of New Music*, 31(1), 236–257.
- Hoffmann, P. (2023). Stochastic Synthesis. iannis-xenakis.org/en/stochastic-synthesis/
- Hoffmann, P. (2022). The Genesis of GENDY3. *Computer Music Journal*, 46(1).
- Xenakis, I. (1992). *Formalized Music*. Pendragon Press.
