# Lineages (generative genealogy)

Lineages is a generative ancestral-tree sequencer for VCV Rack 2. It samples a
finite Kingman coalescent genealogy for 2 to 16 leaves, places neutral signed
mutations on its branches, and preserves that completed tree as a repeatable
musical composition.

In **ANCESTRY**, related polyphonic traits converge as playback moves from the
sampled present toward their most recent common ancestor (MRCA). In **DESCENT**,
one ancestral unison branches into a family of related voltages. The tree changes
only when **NEW** is requested; **RESET** replays the existing tree. A 0 to 10 V
**TIME** input can instead scrub its ancestral age directly.

Lineages is part of the **Coalescent** plugin. It is a finite neutral genealogy
and CV/event instrument, not an audio oscillator or a complete population-genetics
simulator.

![Lineages](img/lineages.png)

## How it works

### A tree viewed backward in time

The leaves along the bottom of the display are the sampled present. Trace any
two leaves upward and their branches eventually meet at a shared ancestor. Trace
all branches upward and the final meeting point is the **most recent common
ancestor**, or MRCA.

Lineages generates this tree with the standard Kingman coalescent. When `k`
ancestral lineages remain, the next waiting time is exponential with rate

```text
k * (k - 1) / 2
```

and one unordered pair is chosen uniformly to merge. With many lineages near
the present there are many possible pairs, so merge events tend to be crowded
together. With only two or three deep lineages, the next common ancestor tends
to take longer. This is why the display usually has dense branching near the
sampled present and longer isolated branches toward the root.

The sampled raw ages are normalized after generation:

```text
present = 0
MRCA    = 1
```

Normalization preserves the relative event spacing inside a tree while making
every generated tree take one complete traversal at the same RATE. The display
therefore shows relative ancestral age, not years or literal generations.

### Neutral mutations and musical traits

**MUTATE** sets the expected total number of mutations across the completed
tree:

```text
expected mutations = 64 * MUTATE^2
```

The actual count is a Poisson draw and can differ on every NEW tree. Each
mutation is placed along the total branch measure, so longer branches are more
likely to receive one, and it receives an independent `+1` or `-1` sign. The panel
control deliberately sets the expected *total* rather than a literal
per-generation mutation rate; this keeps the musical density comparable between
trees with different raw lengths.

The root trait is `0 V`. At any ancestral age, a leaf inherits every mutation
between that age and the sampled present on its path through the tree:

```text
TRAIT = signed inherited mutation sum * STEP / 12 V
```

STEP is measured in semitones, so the output can be used directly as a
V/OCT-style pitch offset. It changes the voltage spread immediately but does not
alter the genealogy. Final voltages are clamped to `-10 ... +10 V`; individual
trees are never normalized to an arbitrary common range.

The signs are an audible neutral trait, not nucleotide states. Opposite signs
can cancel to the same voltage even when two samples have different mutation
histories. DIVERSITY therefore counts inherited mutation differences directly
instead of measuring voltage variance.

### Ancestry and descent

- **ANCESTRY** starts at the sampled present (`0`) and moves toward the MRCA
  (`1`). Mutation differences disappear and related channels merge into exact
  equality. Every TRAITS channel is exactly `0 V` at the MRCA.
- **DESCENT** starts at the MRCA (`1`) and moves toward the sampled present
  (`0`). All channels begin in exact unison, then mutations split them into a
  related cluster.

TRAITS keeps one permanent channel per generated leaf throughout playback.
Channels do not disappear when ancestral lineages merge, and their left-to-right
order matches the leaves in the display. Equality between channels is therefore
part of the genealogy rather than a change in polyphonic cable shape.

## Controls

| Control | Range and default | Purpose |
| --- | --- | --- |
| **RATE** | `-6 ... +5` octaves, default `0` | internal traversal speed `0.1 * 2^RATE` complete trees/s; default is about 10 seconds per tree |
| **SAMPLES** | snapped `2 ... 16`, default `8` | leaf and TRAITS channel count for the next NEW tree; structural |
| **MUTATE** | `0 ... 1`, default `0.4` | expected `0 ... 64` total mutations using `64*u^2`; structural |
| **STEP** | `0 ... 12` semitones, default `2` | immediate voltage contribution of each signed mutation; non-structural |
| **DIRECTION** | **ANCESTRY / DESCENT**, default **DESCENT** | chooses internal movement direction without moving the current cursor |
| **LOOP** | off/on, default on | teleports to the direction's source after reaching its destination |
| **NEW** | button | generates a new tree from current SAMPLES/MUTATE and advances the local RNG |
| **RESET** | button | replays the current tree without changing it or consuming RNG |

RATE CV is added to the RATE knob at 1 V/oct. The combined exponent is clamped
to `-12 ... +5`, so the maximum internal speed is 3.2 complete trees/s and
extreme negative CV remains finite.

SAMPLES and MUTATE describe the *next* tree. Moving either away from the values
used to make the current tree lights the yellow indicator beside NEW and draws
an amber border around the display. The current topology, mutations, output
channel count, and playback continue unchanged. Returning both controls to the
generated values clears the indication; NEW installs a matching tree and also
clears it. RESET never clears or changes this pending state.

Mutation storage has a hard cap of 128. The exposed maximum has an expected
count of 64, so reaching the cap is an extreme Poisson tail rather than normal
operation.

### NEW, RESET, Initialize, and Randomize

These commands intentionally have different meanings:

- **NEW** creates a new topology, branch timing, mutation placement, and sign
  pattern. With TIME unpatched it moves to the selected direction's source.
- Panel/input **RESET** moves to that source on the existing tree. It does not
  consume RNG, change topology, or change mutations, so another replay is exact.
- Rack's context-menu **Initialize** first restores all parameter defaults, then
  generates a fresh default-parameter tree and resets its transport. It is not
  the same operation as the panel RESET button.
- Rack **Randomize** randomizes parameters, assigns a fresh local random stream,
  and generates a matching tree.

With TIME connected, NEW and RESET install or replay at the voltage-selected
age instead of forcing an endpoint, and do so without emitting historical event
pulses.

## Inputs and transport

| Input | Meaning |
| --- | --- |
| **RATE** | direct internal-rate CV at 1 V/oct |
| **TIME** | absolute ancestral age: `0 V` is the sampled present and `10 V` is the MRCA |
| **NEW** | rising edge performs the same operation as the NEW button |
| **RESET** | rising edge performs the same operation as the panel RESET button |

With TIME unpatched, RATE advances one scalar cursor. Changing DIRECTION midway
reverses movement from the current age; it does not reset the tree or fire an
event. With LOOP off, transport stops at its destination until NEW, RESET, or a
direction change makes movement possible again.

With TIME patched, its finite voltage is clamped to `0 ... 10 V` and controls the
cursor absolutely. RATE and LOOP no longer move it. DIRECTION still changes the
display wording but does not reinterpret the voltage. A non-finite TIME value
holds the last valid cursor.

Connecting TIME silently synchronizes to the current input voltage. Disconnecting
it silently hands the same cursor back to internal transport, which then continues
in the selected direction. Neither cable transition replays the event history.

### Event crossings and loops

NODE and MUTATION describe boundaries crossed by the cursor, not every rendered
frame. At an exact event age, equality belongs to the older/ancestral side. In
terms of an old cursor `a` and new cursor `b`, ancestry crosses events satisfying
`a < event <= b`, while descent crosses events satisfying
`b < event <= a`. Thus ancestry fires when it arrives exactly at a boundary;
descent assigns that exact point the ancestral state and fires as it leaves toward
the present. A complete crossing fires once, and crossing the same boundary in
reverse can fire it again.

A TIME jump may cross many nodes and mutations in one audio sample. Lineages
emits at most one pulse on each event-output type for that jump; pulse count does
not encode how many boundaries were skipped.

A loop wrap is a transport teleport, not a fast traversal through the whole
genealogy:

- An ancestry loop reaches the MRCA, emits its root events, then teleports to
  the present without firing every intermediate event.
- A descent loop reaches the present, teleports to the MRCA, and emits one MRCA
  start pulse without firing intermediate NODE or MUTATION events.

Overshoot is discarded at the boundary, so one update cannot race through
multiple complete trees.

## Outputs

| Output | Meaning |
| --- | --- |
| **TRAITS** | polyphonic signed mutation trait, one channel per leaf in the generated tree; `STEP / 12 V` per inherited mutation |
| **NODE** | 10 V / approximately 1 ms pulse when one or more coalescent nodes are crossed |
| **MUTATION** | 10 V / approximately 1 ms pulse when one or more mutation positions are crossed |
| **MRCA** | 10 V / approximately 1 ms pulse on arrival at the root or an explicit descent start there |
| **LINEAGES** | active ancestral-lineage count mapped from 0 to 10 V |
| **DIVERSITY** | mean pairwise inherited mutation distance at 2 V per mutation, clamped to 10 V |

For `n` generated samples and `k` active ancestral lineages:

```text
LINEAGES = 10 * (k - 1) / (n - 1) V
```

It is `10 V` at the sampled present and `0 V` at the MRCA. NODE also includes
the root node, so ancestry arrival at the MRCA can fire NODE and MRCA together.
Starting descent at the root fires MRCA but is not itself a NODE crossing.

For a currently inherited mutation carried by `a` of the `n` leaves, its
contribution to mean pairwise distance is

```text
a * (n - a) / choose(n, 2)
```

DIVERSITY sums those contributions, multiplies the result by the fixed
`2 V/mutation` scale, and clamps only the final voltage. It is not normalized per
tree. STEP has no effect on DIVERSITY, and signed TRAITS may coincide while
DIVERSITY remains nonzero. Both DIVERSITY and all TRAITS channels are exactly
zero at the mutation-free MRCA.

## Display

The display is a complete genealogy, not an oscilloscope history:

- sampled leaves are the small points along the bottom;
- ancestral age rises upward to the larger MRCA point;
- faint teal branches show the whole immutable tree;
- brighter mint branches show the portion traversed in the selected direction;
- the horizontal pale line is the current cursor;
- bright outlined dots at its branch intersections are the active ancestral
  lineages;
- internal points mark coalescent nodes;
- gold `/` ticks are positive mutations and red `\` ticks are negative
  mutations;
- bright mutation ticks are currently inherited at the cursor, while dim ticks
  are inactive at that ancestral age.

The complete tree appears immediately after NEW and its geometry remains stable
during playback. ANCESTRY/DESCENT appears at the top right. An amber frame and
dot, together with the panel's yellow NEW light, mean that structural controls
are waiting for NEW; they do not mean the current tree is invalid.

The audio thread publishes a self-contained display snapshot at about 45 Hz
through the plugin's shared triple buffer. The UI never follows mutable tree
arrays directly.

## Patching Lineages

### Polyphonic descent

The most literal patch lets one pitch become a related chord:

```text
TRAITS -> polyphonic oscillator V/OCT
oscillator -> polyphonic VCA or mixer
MRCA -> envelope or reset accent
```

Use DESCENT with LOOP off to hear exact unison spread into the sampled cluster.
STEP sets the interval size without redrawing the tree. Add a fixed base pitch
before TRAITS if the oscillator needs an absolute note rather than a bipolar
offset.

### Ancestral rhythm

Run slowly in ANCESTRY and route NODE and MUTATION to different envelopes,
drums, or sample changes. NODE follows mergers; MUTATION follows trait changes,
so their two irregular rhythms are related but not interchangeable. Use LINEAGES
to close a filter as branches merge and DIVERSITY to reduce delay spread or
stereo width on the way to the MRCA.

### Reversible CV scrub

Map a triangle, sequencer, or manually controlled CV into `0 ... 10 V` and patch
it to TIME. The tree becomes a deterministic reversible lookup: moving upward
removes derived mutations, moving downward restores them, and crossing the same
boundary in reverse fires its event again. A bipolar `-5 ... +5 V` LFO needs a
5 V offset before TIME to cover the full tree.

Other useful ideas:

- Trigger NEW only at larger compositional boundaries and use RESET for repeatable
  takes inside the current tree.
- Hold TIME near a dense group of young nodes, then sequence small movements
  across it for a related event pattern without generating another genealogy.
- Turn STEP to zero to mute pitch separation while retaining NODE, MUTATION,
  LINEAGES, and mutation-based DIVERSITY.
- Use fewer SAMPLES for exposed melodic counterpoint and 16 for denser clustered
  modulation. Remember to press NEW after changing SAMPLES.
- Duplicate an authored module for identical parallel processing. Randomize one
  copy before generating if the copies should have independent future trees.

### Demo patches

Three generated patches in `patches/` use only Lineages, Core, and Fundamental:

- `lineages_1_polyphonic_descent.vcv` turns TRAITS into a quiet eight-voice sine
  cluster that branches from unison;
- `lineages_2_ancestral_rhythm.vcv` sends NODE and MUTATION to separate envelope
  voices while LINEAGES and DIVERSITY shape their pitches;
- `lineages_3_cv_scrub.vcv` drives TIME with a slow unipolar triangle so one tree
  separates and reconverges reversibly.

Regenerate all three with `python3 tools/make_patch_lineages.py`.

## State, determinism, performance, and limits

- **The complete generated tree is saved with the patch:** topology, raw and
  normalized ages, planar channel order, mutations, generated SAMPLES/MUTATE,
  cursor, running state, direction, loop, pulse remainder, and the complete local
  RNG state. Reload restores the current outputs and composition and, on the
  same build and platform, also makes the next NEW produce the same replacement
  tree it would have produced before saving.
- Unsigned RNG fields are stored losslessly. Restored counts, indices, topology,
  times, masks, mutation positions/signs, cursor, pulses, and RNG state are all
  validated before installation. Invalid module state falls back atomically to a
  safe default tree rather than installing a partial genealogy.
- A duplicated module inherits the authored tree and next RNG state, so matching
  NEW requests can remain identical. Rack Randomize one copy to give it an
  independent random stream.
- The saved current tree is explicit data rather than a seed-only reconstruction.
  Exact next-NEW replay is guaranteed on the same build and platform. Future NEW
  generation uses ordinary floating-point `log`/`exp`, so bit-identical raw branch
  ages are not promised across different platform math libraries or compiler
  builds.
- Save snapshots are published at 500 Hz. A patch save can therefore capture
  transport and pulse presentation up to one control interval (about 2 ms)
  behind the wall-clock save instant, but the captured tree, playback state, and
  next RNG are internally consistent and restore exactly.
- The generator uses fixed arrays with hard limits of 16 leaves, 31 nodes, and
  128 mutations. NEW is bounded and infrequent. Normal audio-callback work is a
  scalar cursor advance, cached event-boundary updates, at most 16 output writes,
  three pulse generators, and scheduled snapshots; there is no allocation, file
  I/O, JSON work, or tree generation during ordinary playback.
- On the development i5-9600K, the standalone optimized playback benchmark takes
  about 21 ms for one million tiny reversible scrubs and about 441 ms per
  million-equivalent full-tree jumps with all 128 mutations. These are core-only,
  hardware-dependent figures; the latter is a deliberately hostile TIME pattern,
  not ordinary internal playback.
- TRAITS is stepped at mutation boundaries, LINEAGES is stepped at node
  boundaries, and the event outputs aggregate crossings. Lineages does not
  interpolate trait changes, turn event density into a clock grid, or report the
  number of events inside one large TIME jump.
- Tree age is normalized. RATE controls traversal of the finished musical object,
  not an inferred population size or physical number of generations.
- The Kingman model assumes a neutral, unstructured, random-mating ancestry in
  its large-population limit. Lineages has no selection, recombination, linkage,
  migration, population structure, changing population size, sex, or explicit
  individuals beyond the sampled leaves.
- MUTATE controls a musically calibrated expected total, and the signed trait is
  an audible abstraction. The module is suitable for demonstrating and hearing
  genealogical structure, not for quantitative inference from biological data.
- Signed steps can cancel, large STEP/dense mutation trees can reach the output
  clamp, and different NEW trees naturally have different spreads. This variation
  is preserved rather than hidden by per-tree normalization.
- There is no external clock mode in version 1. TIME is absolute ancestral age;
  use a sequencer or stepped voltage source when clocked access is needed.

The SDK-free Lineages stability suite exercises PCG and Poisson statistics,
Kingman waiting times and pair choice, tree topology, branch-proportional
mutation placement, fixed musical cases, reversible crossing semantics, loop
teleports, canonical output equivalence, hostile inputs, bounded worst-case
scrubbing, full persistence validation, and exact next-NEW restoration. The Rack
integration harness separately tests NEW versus RESET, Initialize, dirty-state
return, TIME connection changes and jumps, direction reversal, loop boundaries,
MRCA convergence, finite outputs, display snapshots, and JSON round trips. Both
are maintained with the module: the SDK-free suite runs in `make check`, while
the Rack integration harness is compiled and run separately against the Rack SDK.

## References

- J. F. C. Kingman, [*On the genealogy of large populations*](https://doi.org/10.2307/3213548),
  Journal of Applied Probability 19A (1982), 27-43.
- J. F. C. Kingman, [*The coalescent*](https://doi.org/10.1016/0304-4149(82)90011-4),
  Stochastic Processes and their Applications 13 (1982), 235-248.
