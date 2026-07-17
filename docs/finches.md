# Finches (evolutionary branching)

An evolutionary CV instrument for VCV Rack 2. Finches represents a population
as a density along one continuous trait — body size, beak size, feeding
preference, oscillator pitch, or any other quantity you choose to hear. Under
stabilizing selection it forms one cluster. When competition between similar
phenotypes becomes sufficiently local, that cluster broadens and branches into
two persistent morphs. Part of the **Coalescent** plugin's *Fluctuations* series
— see the [main README](../README.md).

![Finches](img/finches.png)

The density display makes the defining gesture literal: one hill becomes two.
The two accepted peak positions leave as paired pitch CVs, their population
fractions leave as paired VCA-ready CVs, and dedicated events mark an accepted
split or merge.

## How it works

Finches discretizes trait space `x` from −1 to +1 into 64 bins. Its state `p(x)`
is a probability distribution, so the bin masses remain nonnegative and sum to
one. Selection and mutation evolve it in dimensionless evolutionary time `τ`:

```text
∂p/∂τ = p(x) · (F(x, p) − mean(F)) + μ · ∂²p/∂x²

F(x, p) = −1/2 · ((x − e) / σK)²
          − ∫ C(x − y) p(y) dy

C(d) = exp(−d² / (2 σC²))
σC = σK / sqrt(B)
```

- `e` is the environmental optimum, moved by **ENV**.
- `σK` is the stabilizing niche width, set by **NICHE**.
- `μ` is mutation diffusion, set by **MUTATE**.
- `B` is the branching number, set by **COMPETE**.

The mean-fitness subtraction is the replicator term: successful traits take a
larger fraction of a fixed total population. **MASS L/R therefore report
fractions, not absolute animal counts.** Mutation diffuses mass only between
neighbouring traits. The ends use reflecting/no-flux boundaries; traits do not
wrap from +1 back to −1.

### What COMPETE means

COMPETE is not a generic population-level competition gain. It changes how
strongly an individual competes with *similar* phenotypes relative to the width
of the environmental niche. Turning it clockwise raises `B`, narrows `σC`, and
makes intermediate phenotypes more vulnerable to neighbours on either side.

For the idealized local analysis, the fitness curvature changes sign at `B = 1`:

- `B < 1`: competition is broad relative to the niche; one central cluster is
  favoured.
- `B > 1`: competition is local enough for disruptive selection and branching
  becomes possible.

The visible split is not a hard comparator at exactly `B = 1`. Finite mutation,
the 64-bin grid, environmental motion, and the time needed for two viable peaks
all shift the observed transition. High MUTATE can keep a broad population
unimodal even when COMPETE is high.

This is a **trait-structured mutation–selection model inspired by adaptive
dynamics**, not a literal simulation of genes, mating, or reproductive isolation.
Its two peaks are phenotypic clusters (or incipient species in the motivating
theory), not proof that complete biological speciation has occurred.

## Controls and inputs

| Control | Model range | Purpose |
| --- | --- | --- |
| **RATE** | `−8 … +4` oct | evolutionary speed; 0 = 8 `τ`/s, and each octave doubles or halves it |
| **MUTATE** | `μ = 0.00002 … 0.00012` | neighbouring-trait diffusion; low preserves distinct peaks, high smooths and can merge them |
| **COMPETE** | `B = 0.35 … 2.8` | relative competition locality; `B = 1` is near 0.505 on the knob |
| **NICHE** | `σK = 0.20 … 0.45` | width of the environmental opportunity around the optimum |
| **MUTANT** | `−1 … +1` | seed offset up to ±0.70 trait units from the current environment |
| **SEED** | button | mixes a narrow 2% mutant cohort into the population |

The default ecology is deliberately stabilizing: COMPETE 0.30 corresponds to
`B ≈ 0.65`, so it begins and remains one peak. Raise COMPETE past roughly noon to
enter the branching side. MUTATE defaults to 0.35 (`μ ≈ 0.000037`), NICHE to
about 0.50 (`σK ≈ 0.32`), and RATE to 0.

**MUTATE** and **COMPETE** have bipolar attenuverters and dedicated CV inputs;
10 V spans their full normalized knob range at full attenuation. The **RATE**
input is exponential too: +1 V doubles the evolutionary speed. **ENV** moves the
environmental optimum by 0.11 trait units per volt and clamps at `−0.55 … +0.55`;
it is not another NICHE-width control. Parameter and CV values are sanitized
before they reach the field solver. The combined RATE knob and CV is capped at
`+4.25` octaves, so positive CV stops tracking 1 V/oct near the knob maximum.

**SEED** accepts a rising trigger as well as the panel button. MUTANT is relative
to ENV, so moving the environment also moves the intended seed location; the
result is clamped to the interior trait range. A seed is only an invasion
attempt: it does not force a split and does not fire SPLIT immediately. Selection
may eliminate it. **RESET** returns to one narrow ancestor at the current
environment without generating a false merge event. Rack's context-menu
**Initialize** additionally returns all panel parameters to their defaults before
constructing that ancestor.

## Outputs

| Output | Meaning |
| --- | --- |
| **MASS L** | fraction of population in the lower-trait basin, scaled 0 … 10 V |
| **MASS R** | fraction in the higher-trait basin, scaled 0 … 10 V |
| **PITCH L** | lower-trait position at 2 V per trait unit (full domain −2 … +2 V) |
| **PITCH R** | higher-trait position at 2 V per trait unit (full domain −2 … +2 V) |
| **SPREAD** | standard deviation scaled as `min(spread / 0.5, 1) · 10 V` |
| **SPLIT** | 10 V / approximately 1 ms pulse when a viable two-cluster state is accepted |
| **MERGE** | 10 V / approximately 1 ms pulse when that state returns to one cluster |

Before a split, both PITCH outputs follow the same unresolved population mean,
MASS L carries the full 10 V population, and MASS R is 0 V. After a split, the
peak identities are always sorted by trait coordinate, so they cannot exchange
merely because one becomes taller. MASS is integrated on each side of the valley
rather than sampled only at the peak bin.

MASS L/R, PITCH L/R, and SPREAD use a 20 ms audio-rate one-pole presentation
smoother, removing the nominal 500 Hz field cadence without adding ecological
updates. SPLIT and MERGE follow the accepted detector state immediately and are
not smoothed.

SPLIT and MERGE use prominence, separation, basin-mass thresholds, persistence,
and different enter/leave thresholds. A broad or momentarily flat-topped hill
therefore does not chatter between states. These are ecological events, not a
periodic clock.

## Patching Finches

The most literal patch uses two oscillators and two VCAs:

```text
PITCH L/R → oscillator V/OCT
oscillators → VCA L/R audio inputs
MASS L/R  → VCA L/R CV inputs
```

It begins as one audible oscillator. When a viable second cluster establishes,
its mass opens the second VCA while the two pitch CVs diverge. Pan the voices
apart to hear branching as both interval and spatial separation.

Other useful patches:

- A slow LFO into ENV makes established clusters follow a moving opportunity.
- Gates into SEED repeatedly challenge the population with the phenotype chosen
  by MUTANT; only ecologically successful challenges become audible.
- Slow modulation of COMPETE crosses between stabilizing and disruptive
  regimes, producing SPLIT and MERGE events on accepted transitions.
- SPREAD can open a filter or increase modulation depth while one peak first
  broadens, before the discrete SPLIT event arrives.
- RATE can be slowed for a long-form transition or raised to audition the
  eventual outcome of a parameter change.

## Display and detection

The screen draws the current trait density, the environmental opportunity, and
the accepted cluster positions. The density is the current field, not an
oscilloscope history: its horizontal axis is phenotype and its vertical axis is
population mass.

Peak positions receive sub-bin refinement for smoother pitch CV. A candidate
split must have separated local maxima, a sufficiently deep intervening saddle,
meaningful population on both sides, and persistence in evolutionary time.
Merge detection is deliberately more permissive than split detection, providing
hysteresis around shallow valleys.

## State, performance, and limits

- **The complete musical state is saved with the patch.** Schema version 1 stores
  all 64 bin masses, the accepted split latch, and the pending split/merge
  persistence timers. A split held in the detector's weak hysteresis band and a
  transition partway toward its next event therefore resume as authored. Invalid
  versions, masses, latches, or timers are rejected transactionally. Older
  density-only patches still load; because they did not store detector history,
  their initial split state is inferred from the stronger entry threshold.
- SPLIT and MERGE pulses are transient signals, not saved state. Loading never
  manufactures an event or replays a pulse that occurred before the save. A
  pending detector timer continues after loading and fires only when its
  remaining persistence time completes. The lock-free save snapshot is refreshed
  on every 500 Hz ecological update; a save racing that publication can receive
  the previous complete field frame, at most roughly one 2 ms interval old.
- The solver is a fixed 64-bin field. Selection uses a bounded `O(64²)` Gaussian
  competition convolution; mutation uses a linear-time implicit tridiagonal
  solve. It allocates nothing in the audio callback.
- Reaction uses a positivity-preserving exponential update. Mutation uses
  backward-Euler diffusion, followed by normalization. The numerical workload is
  bounded even under hostile CV or corrupted patch input. Positive masses below
  `1e-30` are treated as numerically extinct before normalization, preventing
  float-subnormal slowdowns without affecting a biologically meaningful bin.
- Finches is intended as a slow CV instrument. The ecological field runs at a
  fixed 500 Hz control cadence and the display receives lock-free snapshots; it
  is not an audio-rate waveshaper. RATE is bounded so a control tick cannot ask
  the field core to advance more than 0.32 `τ` at once.
- The fixed trait interval is finite. An environment driven hard against its
  allowed edge has less room on one side and can produce visibly asymmetric
  clusters. The central ENV range is limited to reduce this boundary effect.
- The output contract intentionally recognizes at most one accepted pair. More
  complicated shoulders remain visible in the density, but outputs describe the
  strongest viable low/high partition rather than inventing extra voices.
- The discretized replicator–mutator field is a musically controllable rendering
  of evolutionary branching. It should not be used as a quantitative population
  genetics simulator.
- On the development i5-9600K, the standalone `-O2` field benchmark uses about
  0.3% of one core at the default RATE and 4.5% at the maximum bounded RATE
  request, assuming 500 updates/s. These are diagnostic hardware-dependent
  figures, not real-time deadlines; the stability suite uses a deliberately loose
  ceiling and separately verifies that reachable states contain no subnormals.
- Saved values are portable and reload exactly as finite IEEE-754 floats. Future
  evolution is deterministic for the same build and controls, but different
  compilers, `libm` implementations, or fast-math choices can eventually diverge
  by rounding around a split/merge threshold. This is not cross-platform
  bit-identical simulation.

`tools/stability/finches.cpp` exercises the defining one-to-two-to-one gesture,
the full exposed parameter box, no-flux boundaries, deterministic replay, state
restore including detector hysteresis and pending timers, hostile inputs, event
counts, the no-subnormal invariant, performance, positivity, normalization, and
timestep convergence. It runs in `make check`.

## Demo patches

`tools/make_patch_finches.py` writes four Core/Fundamental-only patches:

- **finches_1_branching** — cluster positions drive two VCO pitches while their
  masses open two VCAs: one audible phenotype becomes two.
- **finches_2_moving_niche** — a very slow, attenuated LFO moves ENV while two
  established clusters follow without losing their ecological spacing.
- **finches_3_mutant_invasion** — periodic SEED triggers repeatedly challenge a
  branching population with a nearby mutant. The cohort can bias the transition,
  but the aggregate density does not retain lineage labels or prove ancestry.
- **finches_4_split_merge** — slow COMPETE modulation crosses the transition;
  SPLIT and MERGE are auditioned on separate channels.

## References

- S. A. H. Geritz, É. Kisdi, G. Meszéna & J. A. J. Metz, [*Evolutionarily
  singular strategies and the adaptive growth and branching of the evolutionary
  tree*](https://doi.org/10.1023/A:1006554906681), Evolutionary Ecology 12
  (1998), 35–57.
- U. Dieckmann & M. Doebeli, [*On the origin of species by sympatric
  speciation*](https://doi.org/10.1038/22521), Nature 400 (1999), 354–357.
- M. Doebeli & U. Dieckmann, [*Evolutionary branching and sympatric speciation
  caused by different types of ecological interactions*](https://doi.org/10.1086/303417),
  American Naturalist 156 (2000), S77–S101.
