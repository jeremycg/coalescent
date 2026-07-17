# Archipelago (spatial local adaptation)

Archipelago is an eight-voice evolutionary CV instrument for VCV Rack 2. Each
voice is a habitat containing an abundance distribution over a continuous
trait. The environment favours a different trait along the archipelago;
selection pulls each population toward its local optimum, while mutation and
nearest-neighbour migration spread traits within and between habitats.

The result is a playable model of local adaptation. A gentle environmental
gradient forms a geographic pitch cline. Migration pulls the eight voices
together, strong local selection separates them, a central barrier can turn one
cline into two blocks, and a moving climate shifts the habitable range. TRAIT
and MASS leave as paired eight-channel polyphonic CVs, while global outputs
report the mean trait, left/right difference, and current migration activity.
COLONIZE and EXTINCT mark local occupancy changes.

Archipelago is part of the **Coalescent** plugin's *Fluctuations* series. It is
a deterministic spatial trait-density instrument, not an individual-based or
finite-population genetics simulator.

![Archipelago](img/archipelago.png)

## How it works

The state is a fixed `8 x 32` field. Habitat `i = 0 ... 7` contains non-negative
abundance `u[i][k]` in trait bin `k = 0 ... 31`. Trait bins are cell-centred on
`[-1, +1]` and habitat positions are evenly spaced on `[-1, +1]`:

```text
x[k] = -1 + (k + 1/2) * 2/32
z[i] = -1 + 2*i/7
n[i] = sum_k u[i][k]
theta[i] = climate + gradient * z[i]
```

`n[i]` is local population mass and `theta[i]` is the local environmental
optimum. The optimum is deliberately **not clamped** to the represented trait
interval. A large climate displacement can therefore leave edge habitats with
no well-adapted phenotype and cause a genuine range contraction rather than
pinning the opportunity to the display boundary.

In evolutionary time `tau`, the field follows

```text
du[i][k]/dtau =
    u[i][k] * (1 - n[i] - a * (x[k] - theta[i])^2)
  + D * L_trait(u[i])[k]
  + m * sum_(j neighbour i) w[i,j] * (u[j][k] - u[i][k])
```

where:

- `a` is local stabilizing-selection strength;
- `D` is mutation diffusion along trait space;
- `m` is nearest-neighbour migration rate;
- `L_trait` is the centred trait-space Laplacian with zero-flux boundaries;
- `w[i,j]` is an edge permeability: 1 for an open edge and
  `(1 - barrier)^2` for the central `3 <-> 4` edge.

The logistic term sets a unit local carrying mass in a perfectly adapted
habitat. The quadratic mismatch cost couples evolution to demography: a
population dragged far from its local optimum can decline, disappear, and later
be re-established by growth or migrants.

### Numerical update order

The module evaluates the field at the nearest integer sample-divider cadence to
500 Hz (exactly 500 Hz at 48 and 96 kHz). Each requested advance uses the exact
wall duration of that divider, divides model time into substeps `h <= 0.02 tau`,
and applies this symmetric positive split to every substep:

1. Apply an exact half-step of pairwise migration to active spatial edges in
   forward order.
2. Apply a backward-Euler half-step of no-flux mutation diffusion to each
   habitat.
3. Apply half of the local quadratic selection mortality.
4. Apply exact logistic growth for one full substep.
5. Apply the other half of selection mortality.
6. Apply the other backward-Euler mutation half-step.
7. Apply exact pairwise migration half-steps in reverse edge order.

Selection mortality multiplies each bin by
`exp(-h*a*(x - theta)^2/2)`. If mass before the logistic stage is `n`, exact
unit-capacity logistic growth gives

```text
n' = n * exp(h) / (1 + n * (exp(h) - 1))
```

and scales all bins in that habitat by `n'/n`. Each mutation half-step solves
the no-flux tridiagonal system

```text
(I - (h/2) * D * L_trait) u' = u.
```

For one spatial edge with rate `q = m*w`, exact pair mixing over duration `t`
maps bin values `A` and `B` to

```text
mean = (A + B) / 2
halfDifference = (A - B) * exp(-2*q*t) / 2
A' = mean + halfDifference
B' = mean - halfDifference
```

with `t = h/2` for each half-step. This edge map is positive and conserves
mass. Forward/reverse ordering makes the composed stepping-stone update
symmetric even though neighbouring edges share habitats. Backward-Euler
mutation is also positive and mass-conserving under the no-flux boundary.

After RESET, state restoration, and each solver substep, cell abundance below
`1e-280` is treated as numerical zero. This floor is roughly 80 orders of
magnitude below the smallest Gaussian tail seeded by RESET and 268 orders below
the `1e-12` mass used by public metrics. It prevents slow double-subnormal
arithmetic without deleting the standing variation used for climate adaptation.
An isolated habitat whose complete field reaches numerical zero cannot regrow
until migration supplies mass or RESET supplies a new initial condition.

The wrapper requests at most about `0.128 tau` on one nominal 500 Hz tick; the
field core also caps advances at its independent `0.25 tau` safety limit. No
model work is performed at audio rate.

## Controls

| Control | Range and default | Purpose |
| --- | --- | --- |
| **RATE** | `-8 ... +4` octaves, default `0` | evolutionary speed `4 * 2^rate tau/s`; RATE CV is added at 1 V/oct and the combined value is clamped to `-12 ... +4` |
| **SELECT** | `0 ... 1`, default `0.60` | maps to selection strength `a = 8 * u^2 / tau` |
| **MUTATE** | `0 ... 1`, default `0.45` | true zero at the bottom; otherwise `D = 1e-5 * 300^u / tau`, reaching `0.003 / tau` |
| **MIGRATE** | `0 ... 1`, default `0.40` | true zero at the bottom; otherwise `m = 0.002 * 1000^u / tau`, reaching `2 / tau` |
| **GRADIENT** | `-0.85 ... +0.85` trait, default `+0.55` | coefficient `g` in `theta = climate + g*z`; the habitat 1-to-8 optimum difference is `2*g` (1.10 by default), and positive values favour higher traits toward habitat 8 |
| **BARRIER** | `0 ... 1`, default `0` | closes the central `4 <-> 5` panel-labelled link with permeability `(1 - barrier)^2`; 0 is open and 1 is exactly closed |
| **CLIMATE** | `-0.85 ... +0.85` trait, default `0` | shared offset added to all eight environmental optima |
| **RESET** | button | restores the deterministic initial field and event state |
| **ROW / RING** | switch, default **ROW** | chooses the spatial neighbour graph |

The panel numbers habitats 1 to 8, while the equations use zero-based indices.
Thus the central barrier lies between panel habitats 4 and 5.

SELECT, MUTATE, MIGRATE, GRADIENT, and BARRIER have bipolar attenuverters. Their
CV contribution is

```text
0.1 * attenuverter * inputVoltage
```

in the corresponding raw knob coordinate, followed by clamping to the control's
range. At full positive attenuation, 10 V adds one raw unit. RATE uses its direct
1 V/oct input. CLIMATE has no attenuverter: its input adds a fixed
`0.1 trait/V` to the knob before the combined climate value is clamped to
`-0.85 ... +0.85`.

Control and CV changes alter future evolution. They do not redraw, transpose,
or renormalize the population already on screen.

### Row, ring, and barrier

**ROW** uses the seven edges `1-2, 2-3, ... 7-8`. End habitats each have one
neighbour. **RING** adds the `8-1` edge, giving every habitat two neighbours.
The environmental optimum remains a linear 1-to-8 gradient in both modes, so a
non-zero GRADIENT creates a deliberate environmental seam across the added ring
edge.

BARRIER always controls only the `4-5` edge. Closing it disconnects a row into
two four-habitat groups. It does **not** disconnect a ring, because migration can
still travel around the `8-1` side. In a ring it changes one route and can move
the strongest trait discontinuity, rather than creating two isolated halves.

## Inputs

| Input | Meaning |
| --- | --- |
| **RATE** | direct rate CV; 1 V doubles evolutionary speed |
| **SEL** | selection CV through the SELECT attenuverter |
| **MUT** | mutation CV through the MUTATE attenuverter |
| **MIG** | migration CV through the MIGRATE attenuverter |
| **GRAD** | environmental-gradient CV through the GRADIENT attenuverter |
| **BARR** | central-barrier CV through the BARRIER attenuverter |
| **CLIM** | direct shared-climate offset, `0.1 trait/V` |
| **RESET** | rising edge restores the deterministic initial field |

Archipelago has no STEP input. Its field is a continuous deterministic ecology,
not a sequence of discrete sampled generations.

## Outputs

| Output | Meaning |
| --- | --- |
| **TRAIT** | eight channels of local abundance-weighted mean trait at `2 V/trait`, full range `-2 ... +2 V` |
| **MASS** | eight channels of local abundance, with `0 ... 1` mapped to `0 ... 10 V` |
| **MEAN** | global abundance-weighted trait mean at `2 V/trait`, `-2 ... +2 V` |
| **DIFF** | signed right-half minus left-half trait mean at `2 V` per trait difference, full range `-4 ... +4 V` |
| **FLUX** | unsigned current migration activity, normalized to `0 ... 10 V` |
| **COLONIZE** | 10 V / approximately 1 ms pulse when any habitat enters occupied state |
| **EXTINCT** | 10 V / approximately 1 ms pulse when any habitat leaves occupied state |

TRAIT and MASS always carry eight polyphonic channels in habitat order. A raw
local mean is numerically valid above mass `1e-12`, but the reported TRAIT follows
the occupancy latch: it updates while the habitat is occupied, freezes when mass
crosses the EXTINCT threshold, and resumes from the newly established mean on
COLONIZE. This avoids a dying voice jumping to 0 V and producing a false pitch
gesture; MASS is the authoritative measure of whether the voice is present. On
RESET, each habitat starts with mass `0.75` in a narrow `sigma = 0.065` Gaussian
centred on its representable local optimum, so all eight voices begin occupied.

MEAN weights every trait bin by abundance across the entire field. If total mass
falls to `1e-12` or below it holds its last valid global mean. DIFF separately
computes abundance-weighted means for habitats 1-4 and 5-8, then subtracts left
from right. It is 0 V when either half has mass at or below `1e-12`. Its sign
therefore reports cline direction as well as magnitude.

FLUX is not a copy of the MIGRATE knob. For each active neighbour edge, the core
measures the total-variation difference between the two complete trait fields,
weights it by migration rate and edge permeability, averages over the active
graph, and normalizes the result to `0 ... 1`. For edge set `E`, it computes

```text
TV[i,j] = 0.5 * sum_k abs(u[i][k] - u[j][k])
flux = clamp(m * sum_((i,j) in E) w[i,j] * TV[i,j]
             / (2 * numberOfEdges), 0, 1)
FLUX voltage = 10 * flux
```

The `2` in the denominator is the maximum exposed migration rate; unit carrying
mass supplies the abundance reference. FLUX is 0 V when migration is off or all
neighbouring fields match, and rises only when migration is both enabled and has
different populations to mix.

TRAIT, MASS, MEAN, DIFF, and FLUX use a 20 ms one-pole presentation smoother at
audio rate. The smoother removes the nominal 500 Hz control cadence without
inventing extra ecological updates. COLONIZE and EXTINCT follow the occupancy
detector immediately and are not smoothed.

### Occupancy events

Occupancy uses hysteresis on local mass:

```text
unoccupied -> occupied at n >= 0.08  : COLONIZE
occupied   -> unoccupied at n <= 0.03: EXTINCT
```

The gap prevents a population near one threshold from firing alternating
events. The two outputs are monophonic OR events: simultaneous transitions in
several habitats produce one pulse on the corresponding jack. These are
operational local establishment and loss markers. They do not claim that a new
species evolved or that the last individual in a literal population died.

RESET and state restoration initialize the occupancy latches from the restored
field and deliberately emit no event pulse.

## Display

The display is a fixed space-by-trait map, not an oscilloscope. Habitats 1 to 8
run left to right and trait runs from `-1` to `+1` vertically. Each habitat
column shows its 32-bin abundance density. A strong line follows the local trait
means and a thinner line shows the eight environmental optima, making migration
load and local adaptation visible as separation between the two.

The central `4-5` seam shows barrier permeability. RING mode also marks the
`8-1` wrap link so the extra neighbour relation and its environmental seam are
not hidden. Low-mass habitats dim while retaining a faint marker for the held
TRAIT value. A brief highlight identifies the habitat responsible for a
COLONIZE or EXTINCT event.

The trait axis and abundance brightness are fixed rather than auto-scaled:
mutation broadening, climate-driven decline, and local extinction therefore
remain visually comparable over time. An arrow on the upper or lower trait rail
marks an environmental optimum outside the represented interval. The audio
thread publishes a bounded triple-buffered snapshot; the UI never reads mutable
solver state directly.

## Patching Archipelago

The literal eight-voice patch needs one polyphonic oscillator and VCA:

```text
TRAIT -> poly oscillator V/OCT
oscillator audio -> poly VCA audio
MASS -> poly VCA CV
```

At `2 V/trait`, the complete trait interval spans four octaves. Attenuate TRAIT
before V/OCT when a tighter chord is wanted. GRADIENT controls the chord's
direction and spread, SELECT tunes voices toward their local targets, and
MIGRATE pulls neighbouring voices into a smoother family.

Other useful patches:

- Raise MIGRATE and turn BARRIER clockwise. A smooth geographic chord separates
  into left and right groups as the central edge closes.
- Send a very slow triangle or random curve to CLIMATE. The optimum line moves
  first; populations lag, follow, contract, and recolonize. Use MASS for the
  eight voice levels and COLONIZE/EXTINCT for structural percussion.
- Patch DIFF to stereo width, delay time, or interval depth. It crosses through
  0 V when the cline reverses rather than reporting only unsigned spread.
- Use FLUX to open distortion, brightness, or diffusion while unlike neighbours
  exchange population, then close the effect after they converge.
- Modulate MIGRATE against SELECT. The same gradient moves between eight local
  identities and one more coherent travelling body.
- Switch to RING with non-zero GRADIENT. The `8-1` environmental seam becomes a
  persistent site of mixing pressure; FLUX exposes the resulting tension.
- Set MUTATE or MIGRATE exactly to zero to remove that process, rather than
  merely making it slow.

## How it differs from Islands and Finches

The three modules share evolutionary language but expose different models and
musical gestures.

| | **Archipelago** | **Islands** | **Finches** |
| --- | --- | --- | --- |
| State | continuous trait density and abundance in eight spatial habitats | frequency of one A/B allele in four finite demes | one continuous trait density in one habitat |
| Dynamics | deterministic selection, mutation diffusion, nearest-neighbour migration, and demography | stochastic Wright-Fisher binomial sampling with selection, mutation, and a common migrant pool | deterministic mutation-competition ecology with frequency-dependent branching |
| Spatial structure | row or ring; explicit gradient, climate offset, and central barrier | four demes coupled all-to-all through one migrant pool | none |
| Main musical output | poly-8 local trait pitches paired with poly-8 masses | four mono `0 ... 10 V` bounded allele-frequency lanes | two accepted peak pitches and their masses |
| Defining events | local colonization and extinction | global fixation, loss, sweep, and founder bottleneck | accepted split and merge |
| Primary gesture | a geographic cline follows, lags, contracts, or breaks across a barrier | finite populations wander and fix | one phenotype visibly branches into two |

In Islands, an "island" is a population-genetic deme and MIGRATE mixes expected
allele frequencies through a common pool before stochastic sampling. In
Archipelago, a habitat is a location on a graph and migration physically moves
the complete trait-density field only across graph edges. Archipelago has no
finite-population drift, allele fixation, founder sampling, or random-number
generator.

Archipelago and Finches both discretize a continuous trait, but Finches uses
frequency-dependent competition in one population to render evolutionary
branching. Archipelago has eight locally selected populations and renders
spatial local adaptation. A sharp cline or barrier-separated pair is not a
speciation claim, and Archipelago does not detect branches.

## State, performance, and limits

- **The complete evolved field is saved with the patch:** all 256 abundance
  cells, held local trait reports, the global fallback trait, and occupancy
  latches. Reloading continues the authored cline or range rather than silently
  rebuilding it from the knobs.
- State restoration validates dimensions and rejects non-finite or negative
  abundance. Derived metrics are recomputed from the current controls. One-shot
  pulse memory is cleared rather than serialized, so restoration emits neither
  COLONIZE nor EXTINCT.
- RESET is deterministic. It uses the current sanitized GRADIENT and CLIMATE,
  places mass `0.75` in each habitat as a narrow Gaussian around that local
  optimum, and clips only the reset centre to the representable trait domain.
  It is a new initial condition, not a random redraw.
- The model owns no RNG. Matching state, controls, CV, topology, and event timing
  produce the same deterministic evolution, subject to ordinary floating-point
  differences between builds and platforms.
- The solver runs at a nominal 500 Hz over a fixed `8 x 32` field. Reaction and migration
  are bounded linear passes; mutation uses eight fixed tridiagonal solves. It
  allocates nothing in the audio callback. Audio-rate work is limited to input
  handling, output smoothing, triggers, and scheduled snapshots.
- On the development i5-9600K, the standalone `-O2` field benchmark uses about
  `0.3%` of one core at default RATE and `1.5%` at maximum combined RATE under
  worst-case exposed controls. These are core-only lower bounds: Rack wrapper and
  UI work are excluded, and percentages vary with hardware and build flags.
- The finite `[-1, +1]` trait interval is intentional. Optima may move beyond it,
  but no unseen phenotype exists outside the field. Extreme climate and gradient
  settings can therefore cause edge loss.
- Mutation is deterministic diffusion, not a sampled mutation count. Migration
  is deterministic mass exchange, not individual dispersal. There is no genetic
  drift, recombination, dominance, linkage, sex, age structure, or evolving
  dispersal rate.
- All habitats share carrying capacity, selection, mutation, and migration
  controls. They differ through environmental optimum, graph position, current
  state, and the central edge permeability.
- MASS is continuous abundance. COLONIZE and EXTINCT are deliberately thresholded
  musical interpretations of that field, not extra discontinuities in the
  ecology.

The standalone Archipelago stability suite exercises positivity, bounded
hostile inputs, true-zero mutation and migration, row/ring conservation, exact
barrier closure, cline direction, climate-driven range changes, occupancy
hysteresis, numerical-extinction handling, standing-variation recovery,
deterministic replay, state restoration, output bounds, and timestep convergence.
It runs in `make check`.

## Demo patches

`tools/make_patch_archipelago.py` writes four Core/Fundamental-only patches:

- **archipelago_1_local_adaptation** - the paired TRAIT and MASS polyphonic
  outputs become an eight-voice locally adapted chord.
- **archipelago_2_migration** - slow migration modulation moves the voices
  between local identities and a coherent family.
- **archipelago_3_climate_shift** - a moving climate contracts and restores the
  audible range while COLONIZE and EXTINCT mark transitions.
- **archipelago_4_barrier** - high migration meets a nearly closed central edge,
  separating the left and right halves.

## References

- J. B. S. Haldane, [*The theory of a cline*](https://doi.org/10.1007/BF02986626),
  Journal of Genetics 48 (1948), 277-284.
- M. Slatkin, [*Gene flow and selection in a cline*](https://doi.org/10.1093/genetics/75.4.733),
  Genetics 75 (1973), 733-756.
- M. Kirkpatrick and N. H. Barton, [*Evolution of a species' range*](https://doi.org/10.1086/286054),
  American Naturalist 150 (1997), 1-23.
- J. Polechova and N. H. Barton,
  [*Limits to adaptation along environmental gradients*](https://doi.org/10.1073/pnas.1421515112),
  Proceedings of the National Academy of Sciences 112 (2015), 6401-6406.
