# Axon (FitzHugh–Nagumo)

Spiking-neuron oscillator for VCV Rack 2, and one half of Coalescent's **neuron
pair** — the twin of [Soma](soma.md), with which it shares an RK4 integration
strategy. Part of the **Coalescent** plugin's *Fluctuations* series — see the
[main README](../README.md).

Both neuron modules are **polyphonic** (up to 16 voices) — see [Polyphony](#polyphony).

![Coalescent neuron pair](img/neuron.png)

Axon is a voice built on the **FitzHugh–Nagumo** (FHN) model — a two-variable
reduction of the Hodgkin–Huxley neuron. The membrane voltage `v` is the audio
output. Above a current threshold the neuron free-runs as a relaxation
oscillator (a slow charge-up followed by a fast spike); below threshold it sits
at rest and fires exactly one spike each time you trigger it. One knob —
**CURRENT** — moves you across that boundary (a Hopf bifurcation), so the same
module is a drone oscillator at one setting and a percussion/transient voice at
another.

## How it works

The state is two coupled variables in dimensionless time:

```
dv/dt = v − v³/3 − w + I        (fast: membrane voltage)
dw/dt = ε·(v + a − b·w)          (slow: recovery)
```

`v` rushes along the cubic nullcline and snaps back (the spike); `w` is the slow
recovery that drags it down and sets up the next spike. **ε** (EPS) is the ratio
of the two timescales — small ε gives a sharp relaxation spike, large ε a
smoother, near-sinusoidal swing. **a** (SHAPE) shifts the asymmetry / threshold.
`b` is fixed at 0.8.

The system is **stiff** (a fast variable riding a slow one), so the integrator is
the real work: each sample takes a number of **RK4 substeps**, and that number
adapts to pitch to hold the substep size near `0.05` — up to a cap of 64 substeps,
above which (roughly the top octave) the step grows and pitch trades some
integration accuracy for bounded CPU. A finiteness reset + state clamp are the
backstop if forcing ever pushes it to run away. The `f()` derivative and the RK4
step are written so the **Hindmarsh–Rose** sibling (Soma) uses the same integration
strategy with one extra equation.

### Pitch is the simulation *speed*

The limit-cycle period is **emergent** — it depends on CURRENT, EPS and SHAPE —
so pitch is **open-loop calibrated, not phase-locked.** V/OCT maps to how fast
dimensionless time advances, calibrated (`RATE_CAL`) so the default voicing reads
C4 at 0 V. Tracking is within ~1 cent across the useful range at default params,
but **changing CURRENT / EPS / SHAPE detunes the pitch somewhat** — that coupling
is deliberate and part of the instrument's character, not a bug.

## Controls

| Control | Range | Purpose |
| --- | --- | --- |
| **PITCH** | ±4 oct | simulation speed (audio pitch); 0 = C4 |
| **CURRENT** | −0.2 … 1.6 | injected current `I`; the excitability / bifurcation control. Rests at both ends, oscillates in a middle band (~0.33–1.42 at default shape) |
| **EPS** | 0.01 … 0.30 | timescale ratio `ε`; small = sharp spike, large = smoother/near-sine |
| **SHAPE** | 0.4 … 1.0 | waveform asymmetry `a` (threshold position) |

**CURRENT** and **EPS** each have an attenuverter + CV input (±5 V). **V/OCT**
sums with PITCH (no attenuverter). **TRIG** injects a short decaying current pulse
on each rising edge — from rest that fires one spike (percussion); inside the
oscillating band it perturbs the phase. **SYNC** is a hard reset: a rising edge
re-seeds the orbit at the rest fixed point, so clocking it locks the cycle (and
sweeping a master against it gives the classic hard-sync timbre). There is no
mode switch: the regime is set purely by where CURRENT sits, and triggers are
honoured in both.

Outputs: **OUT** — the membrane voltage `v`, soft-clipped to ±5 V (`tanh`) and
internally DC-blocked at ~20 Hz (the limit cycle's mean is not zero). **SPIKE** —
a 10 V / ~1 ms pulse on each upward threshold crossing of `v` (one per spike, with
hysteresis so a noisy peak can't double-fire). **W** — the recovery variable as a
slow correlated ±5 V CV, intentionally *not* high-passed.

## Display

The screen traces the **phase portrait** in the `(v, w)` plane — the FHN limit
cycle is a glowing closed orbit, and a trigger from rest shows as a single loop
that jumps out and relaxes back. The two faint guides are the **nullclines**
(the cubic `v`-nullcline and the straight `w`-nullcline); their intersection is
the fixed point whose stability CURRENT controls, so you can watch the orbit grow
as CURRENT crosses into the oscillating band. The trail is read lock-free from a
~45 Hz snapshot — fine for a visualiser.

## Polyphony

Both neuron modules are polyphonic, up to **16 voices**, each an independent neuron
with its own integration state. The voice count is taken from the **V/OCT** cable's
channel count, falling back to **TRIG** — so a polyphonic gate/trigger (with no
pitch patched) gives you polyphonic percussion. Every CV input (CURRENT, EPS /
BURST, TRIG) is read per voice and normalled to channel 0 when it carries fewer
channels, and **OUT / SPIKE / W (Z)** are polyphonic. Knobs and attenuverters are
shared across all voices.

The display traces **every** active voice at once: each orbit is drawn on its own
hue, stepped across a narrow band around the module's accent colour (cyan for
Axon, amber for Soma), with a small `Nv` voice-count badge in the corner.

Feed them from a polyphonic source (a poly MIDI→CV, or `VCV Split`/`Merge`) and
sum the poly OUT with any mixer.

Polyphonic voices **keep their internal state** when the channel count drops and
rises again — a voice that goes silent isn't cleared, so it resumes its orbit
where it left off rather than restarting from rest. That's deliberate (continuity
across held chords), but if you want each new note to start clean, gate it with
the SYNC input.

## CPU

The integrator is the cost: per voice, Axon runs `oversample × K` RK4 substeps
per audio sample, where `K` (4…64) rises with pitch. So the worst case scales as
**voices × oversample × substeps** — at 16 voices, ×8 anti-aliasing, and the top
octave that's `16 × 8 × 64 ≈ 8000` steps/sample. That's a ceiling, not the normal
case (at moderate pitch `K` sits near its floor of 4), but it means a big
polyphonic patch at ×8 is genuinely heavy.

| Voices | Anti-aliasing | Cost |
| --- | --- | --- |
| 1 | Off / ×4 | light |
| 8–16 | ×4 (default) | moderate |
| 16 | ×8 | heavy, patch-dependent |

Rule of thumb: leave anti-aliasing at ×4 (or Off), and only reach for ×8 if you
hear aliasing on high notes. [Soma](soma.md) is a touch heavier still (a
three-variable system vs Axon's two).

## Patches

`tools/make_patches_neuron.py` writes seven Axon smoke-test patches into `patches/`
(and copies them into the Windows Rack patches folder if present):

- **axon_1_freerun** — default voicing → audio; play V/OCT
- **axon_2_blips** — sub-threshold CURRENT, an LFO square clocking TRIG (one spike
  per clock); take OUT for the spike voice
- **axon_3_selfevolving** — W self-patched into CURRENT CV: a slow wandering
  texture that rides its own recovery variable
- **axon_4_crossmod** — VCO SAW into CURRENT CV for FM-like sidebands
- **axon_5_sync** — a master VCO square clocking SYNC: hard sync. Sweep Axon's
  PITCH against the master for the classic sync-sweep timbre
- **axon_6_poly** — 4 voices spread across both pitch (a chord) and CURRENT, so
  each traces a differently-sized limit cycle: four coloured orbits on the scope
  (the portrait is pitch-invariant, so CURRENT is what separates them). OUT is
  summed back to mono through Sum
- **axon_7_midipoly** — a *playable* polyphonic voice: MIDI→CV drives V/OCT, its
  GATE opens an ADSR→VCA after Axon so notes start/stop (Axon free-runs, so the VCA
  is what you "play"). Pick your MIDI device and set the module's Polyphony channels
  in its right-click menu

## Notes / known limits

- **Pitch is emergent / approximate.** CURRENT, EPS and SHAPE pull the pitch a
  little (see above) — deliberate, not a bug. Calibration targets C4 at the
  default voicing.
- **Aliasing.** Spikes are sharp and the `tanh` soft-clip adds harmonics, so high
  notes can alias. The right-click **Anti-aliasing** option (Off / ×4 / ×8,
  default ×4) oversamples the whole output chain — DC-block + tanh — and decimates
  with a windowed-sinc FIR, which band-limits both the spike and the nonlinearity.
  Higher factors cost more CPU (scaling with voice count); turn it Off if you're
  CPU-bound and not playing high notes.
- **State is not saved.** `v`, `w` are transient and re-seed at rest on load;
  params persist. Deliberate.
- **DC.** OUT is DC-blocked (the limit-cycle mean ≠ 0). W is intentionally not
  blocked — it's a slow correlated CV.
- **SR changes** need no handler: the only SR-derived cached state (the DC-blocker
  cutoff) is refreshed when the sample rate changes; everything else recomputes
  per sample.

`tools/stability/axon.cpp` is a standalone replica of the kernel: it measures the
dimensionless period to set `RATE_CAL`, sweeps CURRENT × EPS × SHAPE × pitch and
asserts `v`,`w` stay finite/bounded, and checks V/OCT tracking.
`tools/render_wav_axon.cpp` auditions voicings offline (writes WAVs).

```bash
g++ -O2 -o /tmp/t tools/stability/axon.cpp && /tmp/t      # exit 0 = pass
python3 tools/panel_diagram_neuron.py                     # panel footprint check
```

## References

- R. FitzHugh, *Impulses and physiological states in theoretical models of nerve
  membrane*, Biophysical Journal 1 (6), 1961.
- J. Nagumo, S. Arimoto, S. Yoshizawa, *An active pulse transmission line
  simulating nerve axon*, Proc. IRE 50 (10), 1962.
- [FitzHugh–Nagumo model — Wikipedia](https://en.wikipedia.org/wiki/FitzHugh%E2%80%93Nagumo_model)
- [Scholarpedia: FitzHugh–Nagumo model](http://www.scholarpedia.org/article/FitzHugh-Nagumo_model)
