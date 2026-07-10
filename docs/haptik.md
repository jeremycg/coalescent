# Haptik

A scanned-synthesis oscillator for VCV Rack 2. Part of the **Coalescent**
plugin's *Fluctuations* series — see the [main README](../README.md).

![Haptik](img/haptik.png)

Haptik builds its waveform from a little physical system: a ring of masses joined
by springs. You set the ring moving (a pluck, a bump, noise, or an external
signal), and a phase pointer sweeps around the ring at audio rate, reading off
the displacement of each mass in turn. That sequence of displacements *is* the
output waveform. Because the sweep rate (pitch) is independent of the spring
physics, the same moving shape can be played at any pitch — or frozen and played
as a static wavetable.

## How it works

Scanned synthesis, invented by Max Mathews, Bill Verplank and Rob Shaw, separates
two things a normal oscillator fuses together:

- **A dynamical system** — here, a closed ring of `N` masses. Each mass is pulled
  toward its neighbours (the **COUPLE** stiffness) and toward its rest position
  (the **RATE**/centering force), and loses energy over time (**DAMP**). In Fast
  mode it is stepped once per audio sample with a symplectic (semi-implicit) Euler
  update (in Slow mode, every 256 samples and interpolated — see below); the step
  is cheap and, with the coupling capped at 0.9, the homogeneous lattice stays
  bounded (forced/lossless runaway is caught separately by a state clamp).
- **A scanning function** — a phase pointer that runs around the ring at the pitch
  frequency, linearly interpolating between adjacent masses to read out a sample.

The payoff is that **timbre and pitch are decoupled.** The shape evolves on its
own terms; the scan just decides how fast you read it. Freeze the shape and you
have a wavetable; let it ring and the wavetable morphs while it plays; drive it
from EXT IN and the ring resonates the input like a string body.

The **MODE** switch picks how fast the ring itself evolves:

- **Fast** (default) steps the lattice every sample, so its modes sit in the
  audible range and the voice behaves like a cross between a wavetable oscillator
  and a Karplus–Strong resonator. Here **COUPLE is the dominant timbral control**
  and **RATE** is subtle.
- **Slow** steps the lattice only every 256 samples, with the coefficients
  rescaled to that longer timestep. The modes drop into the 0–15 Hz "haptic"
  band — the original conception — so you scan a *slowly-morphing* wavetable,
  **RATE becomes a real evolution-rate control**, and FREEZE captures a clean
  table. Inter-frame interpolation keeps the morph smooth (no stepping), and the
  centering term is clamped so the divided step stays stable.

## Controls

| Control | Range | Purpose |
| --- | --- | --- |
| **N** | 8–128 | number of masses / table length |
| **PITCH** | ±4 oct | scan rate (audio pitch); 0 = C4 |
| **RATE** | 0.05–30 Hz | centering force — the ring's restoring/"breathing" rate (subtle vs COUPLE) |
| **COUPLE** | 0–0.9 | neighbour stiffness; the main timbral control (hard-capped at 0.9 for stability) |
| **DAMP** | 0–1 | energy loss; 0 = lossless drone |
| **INJECT** | 0–1 | excitation amount + EXT IN gain |
| **EXCITE** | impulse / bump / noise / drive | excitation shape applied on TRIG |
| **DRIVER** | 0–100% | where excitation / EXT IN enters the ring (default 25%) |
| **FREEZE** | Run / Freeze | hold the current shape as a static tone |
| **MODE** | Fast / Slow | audio-rate resonator vs haptic-rate morphing wavetable |

CV inputs (RATE, COUPLE, DAMP, INJECT) each have an attenuverter. **V/OCT** sums
with PITCH. **TRIG** re-excites with the current EXCITE shape on each rising edge.
**EXT IN** injects external audio as a force at the DRIVER mass (level scales with
INJECT). A polyphonic cable is summed to mono (Haptik has a single lattice).

Outputs: **OUT** — the scanned waveform, soft-clipped to ±5 V and internally
DC-blocked at ~20 Hz. **MOTION** — the displacement of **mass 0** (a fixed point
on the ring, *not* the driver mass) as a ±5 V CV; audio-rate in this build and
deliberately not high-passed.

## Display

The screen up top visualises the ring directly: the `N` masses are nodes on a
circle, each pushed in/out radially by its displacement and joined into the glowing
loop — that loop *is* the current waveform. The dim circle is zero displacement,
the orange node is the DRIVER mass (where excitation / EXT IN enters), and the
yellow dot is the scan read-head. RUN/FREEZE status and the mass count sit in the
corners. Freezing is the clearest demo: the shape locks while the scan dot keeps
sweeping — pitch reading a static table, made visible.

Two deliberate fudges, both for legibility:

- The scan dot sweeps at an illustrative ~0.5 Hz, **not** the true read rate. The
  real read-head runs at audio rate and can't be drawn at 60 fps (it would alias
  into several ghost dots), so this just conveys the concept.
- The ring is auto-scaled to a smoothed peak displacement so it stays a readable
  size instead of zooming with level. "Zero" is always the dim reference circle.

## Playing it

**Out of the box** Haptik defaults to **EXCITE = drive** (continuous excitation),
so patching `V/OCT` + `OUT → audio` gives an immediate, sustaining, evolving tone —
no trigger needed. Switch EXCITE to **impulse / bump / noise** (0–2) for the
struck/plucked behaviour instead.

In the pluck modes Haptik has **no built-in amplitude envelope** — a struck note
rings and decays according to DAMP, not according to how long you hold a key. A
typical playable patch:

```
MIDI→CV ─ V/OCT ─→ Haptik V/OCT
MIDI→CV ─ GATE  ─→ Haptik TRIG        (each note-on plucks)
Haptik OUT ─→ VCA ←─ ADSR (gated by GATE) ─→ Audio
```

The VCA+ADSR lets key-release cut the sound; without it, notes ring out per DAMP.
For reliable re-plucking on legato passages, drive TRIG from the MIDI→CV
**RETRIGGER** output. Haptik is **monophonic**, so set the MIDI module to one
voice; for chords, run one instance per voice. For drones, use DAMP=0 (one bump
rings forever) or EXCITE=3 (continuous drive).

## Patches

`tools/make_patches_haptik.py` writes four smoke-test patches into `patches/`
(and copies them into the Windows Rack patches folder if present):

- **haptik_1_plucked** — a single pluck (wire a clock/gate into TRIG to re-pluck)
- **haptik_2_drone** — lossless drone (DAMP = 0)
- **haptik_3_frozen** — frozen wavetable; play it with V/OCT
- **haptik_4_driven** — VCO → EXT IN, the ring as a resonator

## Notes

- Lattice state (displacement / velocity / scan phase) **is saved** with the patch,
  so an evolved or frozen sound reloads as itself. (Pre-2.2.1 patches, and any with
  malformed data, fall back to a fresh bump.)
- **FREEZE** holds the shape completely: while frozen, TRIG and excitation are
  ignored so plucking can't alter the held waveform (the scan keeps reading it).
- Excitations are zero-mean, so a pluck adds shape but no DC step — no click/thump
  through the output DC-blocker on each trigger.
- Stability has two parts. The **homogeneous** system is unconditionally stable
  with COUPLE ≤ 0.9 (ω_max = √(kCtr + 4·kSpr) < 2 for symplectic Euler) — **do not
  raise that cap.** But external forcing (TRIG, EXT IN) with DAMP=0 (lossless) can
  pump energy without bound, so `y[]`/`v[]` are additionally hard-clamped to ±16 —
  never reached in normal play, it just turns runaway into bounded saturation
  instead of NaN. `tools/stability/haptik.cpp` checks both (homogeneous sweep and a
  forced/lossless case) plus pitch:
  `g++ -O2 -o /tmp/t tools/stability/haptik.cpp && /tmp/t` (exit 0 = pass).

## References

- M. Mathews, W. Verplank, R. Shaw, *Scanned Synthesis*, Proc. International
  Computer Music Conference (ICMC), 2000 — the original technique, developed at
  Interval Research (1998–99).
- [Scanned synthesis — Wikipedia](https://en.wikipedia.org/wiki/Scanned_synthesis)
- [Csound: Scanned Synthesis (`scanu`/`scans`)](https://csound.com/docs/manual/SiggenScanTop.html)
  — the reference matrix-based implementation.
- [SuperCollider `Spring` UGen](https://doc.sccode.org/Classes/Spring.html) and a
  [scanned-synthesis example](http://sccode.org/1-4VP).
