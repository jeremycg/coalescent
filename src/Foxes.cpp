#include "plugin.hpp"
#include "dsp/rk4.hpp"
#include "dsp/display_snapshot.hpp"
#include "tanh_approx.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>

// FOXES — a three-trophic-level food chain: grass → bunnies → foxes. The
// nondimensional Hastings–Powell (1991) model, the three-species sibling of
// Bunnies. Unlike an autonomous two-species ecology (which can only rest or cycle),
// three levels give period doubling, multistability, and deterministic chaos.
//
//   f1(x) = a1·x/(1+b1·x)      f2(y) = a2·y/(1+b2·y)
//   dx/dt = x·(1 − x) − f1(x)·y                  (grass / resource)
//   dy/dt = f1(x)·y − f2(y)·z − d1·y             (bunnies / consumer)
//   dz/dt = f2(y)·z − d2·z                        (foxes / top predator)
//
// Canonical constants a1=5, a2=0.1, d1=0.4, d2=0.01 (fixed). WILD drives the
// bifurcation parameter b1 = 1 + 5.2·wild² (rest → chase → period doubling → chaos
// near b1=3); BALANCE sets the upper trophic saturation b2 = 1.75 + 0.5·balance.
// Neither is a monotonic "chaos amount" — both cross regular and chaotic windows.
//
// The autonomous system is dissipative on the nonnegative orthant (dS/dt =
// x(1−x) − d1·y − d2·z is bounded above), so it needs no amplitude servo: unforced,
// it never reaches the clamps. They are a backstop for corrupt input / numerics and
// for a sustained negative KICK (which legitimately drives grass to the floor). See
// tools/stability/foxes.cpp for the equilibrium/Hopf/Lyapunov/RATE_CAL proofs;
// its constants and maps mirror this file and the build fails if they drift.

static const int TRAIL_N = 2048;   // chronological phase-trail ring (long enough to reveal the chaotic attractor)

// Per-species output/display gains. Equalize the very different native ranges
// (measured centered RMS ≈ [0.071, 0.031, 0.248] at the periodic default) to
// ~matched output RMS. Fixed after calibration — no moving normalization. File
// scope (not a constexpr member) so a runtime-indexed read isn't an ODR-use that
// needs an out-of-line definition under -std=c++11.
static const float FOX_GAIN[3] = {8.f, 18.f, 2.5f};   // grass, bunny, fox

struct Foxes : Module {

    enum ParamId {
        RATE_PARAM,
        BALANCE_PARAM,
        WILD_PARAM,
        BALANCE_ATT_PARAM,
        WILD_ATT_PARAM,
        PARAMS_LEN
    };
    enum InputId {
        VOCT_INPUT,
        BALANCE_INPUT,
        WILD_INPUT,
        KICK_INPUT,
        INPUTS_LEN
    };
    enum OutputId {
        GRASS_OUTPUT,
        BUNNY_OUTPUT,
        FOX_OUTPUT,
        GRASS_PEAK_OUTPUT,
        BUNNY_PEAK_OUTPUT,
        FOX_PEAK_OUTPUT,
        OUTPUTS_LEN
    };
    enum LightId { LIGHTS_LEN };

    // ─── Canonical Hastings–Powell constants (published; not tunable) ─────────
    static constexpr float A1 = 5.f, A2 = 0.1f, D1 = 0.4f, D2 = 0.01f;

    // ─── State (transient; not saved) ─────────────────────────────────────────
    // A deterministic point on the periodic default attractor (b1=2.3, b2=2),
    // derived offline past the transient by tools/stability/foxes.cpp check [7].
    // Seeding here (not at the equilibrium) avoids a silent startup and is
    // deterministic across load/reset.
    static constexpr float SEED_X = 0.812781f, SEED_Y = 0.104883f, SEED_Z = 12.478951f;
    float x = SEED_X, y = SEED_Y, z = SEED_Z;
    float cx = SEED_X, cy = SEED_Y, cz = SEED_Z;    // active center (analytic equilibrium)
    float cachedB1 = -1.f, cachedB2 = -1.f;
    bool  centerValid = false;
    float prevPeak[3] = {}; bool rising[3] = {};    // per-species peak edge state
    dsp::PulseGenerator peakGen[3];

    // ─── Visualizer: live ring + lock-free snapshot ───────────────────────────
    // A chronological ring of fixed-scaled centered coordinates (NOT auto-normalized
    // — a strange attractor has a fixed shape and must not breathe for UI reasons).
    // The widget projects it isometrically. Captured at ~uniform dimensionless-time
    // spacing so the geometry is independent of audio pitch.
    //
    // These coordinate a clean UI hand-off (the display is a stable attractor
    // *portrait*, not a live scope — the audio-rate trajectory can't be followed at
    // frame rate). `resetGen` bumps on reset/recovery so the UI can identify and
    // eventually replace an old latch. `regimeGen` identifies the current control
    // regime, and `clean` says the ring is a FULL portrait of it (TRAIL_N points
    // captured since the last regime change). The UI only replaces an existing latch
    // with a clean snapshot, so a half-old/half-new mix is never promoted and the
    // previous portrait holds while a knob is mid-turn.
    struct Trail { float pt[TRAIL_N][3] = {}; int head = 0; int count = 0;
                   std::uint64_t resetGen = 0; std::uint64_t regimeGen = 0; bool clean = false; };
    Trail liveTrail;
    coalescent::DisplaySnapshot<Trail> displaySnapshot;
    dsp::ClockDivider pubDiv;
    float lastDisplayFs = 0.f;
    float capAccum = 0.f;   // dimensionless time accumulated since the last trail point
    std::uint64_t resetGen = 0, dispGen = 0;
    int   ptsSinceGen = 0;
    float regB1 = -1.f, regB2 = -1.f, regKick = 0.f;   // last display-regime signature

    // ─── Numerical policy (mirrors tools/stability/foxes.cpp; asserted there) ──
    static constexpr float RATE_CAL  = 61.387f;   // default period (tau) → C4; asserted within 5 cents
    static constexpr float HSUB_MAX  = 0.1f;      // validated 0-cent at the periodic default; bounded at max WILD
    static constexpr int   MIN_SUB   = 2;
    static constexpr int   MAX_SUB   = 64;
    static constexpr float PITCH_TOTAL_MIN = -8.f, PITCH_TOTAL_MAX = 8.f;
    static constexpr float STATE_MAX = 1e3f;
    static constexpr float POS_FLOOR = 1e-4f;
    // Full-scale Rack signals remain stable at the production h <= 0.1 policy.
    // The standalone sweep covers the complete sanitized +/-15 V input range.
    static constexpr float KICK_GAIN = 0.25f;
    static constexpr float BALANCE_CV_DEPTH = 0.1f, WILD_CV_DEPTH = 0.1f;
    static constexpr float POP_GATE = 0.15f;              // peak activity threshold on the gained centered value
    static constexpr float GATE_LEVEL = 10.f, GATE_TIME = 1e-3f;
    static constexpr float CAP_TAU = 0.225f;              // trail-capture spacing in dimensionless time

    Foxes() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);
        configParam(RATE_PARAM, -8.f, 4.f, 0.f, "Rate", " Hz", 2.f, dsp::FREQ_C4);
        configParam(BALANCE_PARAM, 0.f, 1.f, 0.5f, "Balance (upper trophic saturation)");
        configParam(WILD_PARAM, 0.f, 1.f, 0.5f, "Wild (rest → chase → period doubling → chaos)");
        configParam(BALANCE_ATT_PARAM, -1.f, 1.f, 0.f, "Balance CV");
        configParam(WILD_ATT_PARAM, -1.f, 1.f, 0.f, "Wild CV");

        configInput(VOCT_INPUT, "1V/oct rate");
        configInput(BALANCE_INPUT, "Balance CV");
        configInput(WILD_INPUT, "Wild CV");
        configInput(KICK_INPUT, "Kick / grass-force");

        configOutput(GRASS_OUTPUT, "Grass (resource)");
        configOutput(BUNNY_OUTPUT, "Bunnies (consumer)");
        configOutput(FOX_OUTPUT, "Foxes (top predator)");
        configOutput(GRASS_PEAK_OUTPUT, "Grass peak");
        configOutput(BUNNY_PEAK_OUTPUT, "Bunnies peak");
        configOutput(FOX_PEAK_OUTPUT, "Foxes peak");
    }

    static inline float f1(float X, float b1) { return A1 * X / (1.f + b1 * X); }
    static inline float f2(float Y, float b2) { return A2 * Y / (1.f + b2 * Y); }

    void reseed() { x = SEED_X; y = SEED_Y; z = SEED_Z; }
    void resetPeakMemory() { for (int i = 0; i < 3; ++i) { prevPeak[i] = 0.f; rising[i] = false; } }
    // Clears the ring, bumps resetGen (so the UI recognizes the replacement), and
    // resets the post-change point counter so a fresh portrait must be re-collected.
    void resetTrail() { liveTrail.head = liveTrail.count = 0; capAccum = 0.f; resetGen++; ptsSinceGen = 0; }

    void onReset(const ResetEvent& event) override {
        Module::onReset(event);
        reseed(); cx = SEED_X; cy = SEED_Y; cz = SEED_Z;
        cachedB1 = cachedB2 = -1.f; centerValid = false;
        resetPeakMemory();
        for (int i = 0; i < 3; ++i) peakGen[i].reset();
        resetTrail();
        regB1 = regB2 = -1.f; regKick = 0.f;
        liveTrail.resetGen = resetGen;
        liveTrail.regimeGen = dispGen;
        liveTrail.clean = false;
        displaySnapshot.writable() = liveTrail;
        displaySnapshot.publish();
        pubDiv.reset();
        lastDisplayFs = 0.f;
    }

    // Analytic interior coexistence point; caches on (b1,b2). Returns false and
    // leaves the previous valid center if any positivity condition fails (a
    // corrupt-CV backstop — never happens inside the exposed control box).
    bool updateCenter(float b1, float b2) {
        if (centerValid && b1 == cachedB1 && b2 == cachedB2) return true;
        float denom = A2 - b2 * D2;
        if (denom > 0.f) {
            float ys = D2 / denom;
            float disc = (b1 + 1.f) * (b1 + 1.f) - 4.f * A1 * b1 * ys;
            if (disc >= 0.f) {
                float xs = ((b1 - 1.f) + std::sqrt(disc)) / (2.f * b1);
                float fx = f1(xs, b1);
                float zs = ys * (fx - D1) / D2;
                if (xs > 0.f && ys > 0.f && fx > D1 && zs > 0.f
                    && std::isfinite(xs) && std::isfinite(zs)) {
                    cx = xs; cy = ys; cz = zs;
                    cachedB1 = b1; cachedB2 = b2; centerValid = true;
                    return true;
                }
            }
        }
        return false;
    }

    // rising→falling local max above POP_GATE on the gained centered value.
    void firePeak(int i, float gained, dsp::PulseGenerator& gen) {
        if (gained > prevPeak[i] + 1e-7f) rising[i] = true;
        else if (gained < prevPeak[i] - 1e-7f) {
            if (rising[i] && prevPeak[i] > POP_GATE) gen.trigger(GATE_TIME);
            rising[i] = false;
        }
        prevPeak[i] = gained;
    }

    void pushTrail(float sx, float sy, float sz) {
        liveTrail.pt[liveTrail.head][0] = sx;
        liveTrail.pt[liveTrail.head][1] = sy;
        liveTrail.pt[liveTrail.head][2] = sz;
        liveTrail.head = (liveTrail.head + 1) % TRAIL_N;
        if (liveTrail.count < TRAIL_N) liveTrail.count++;
        // Count to 2·TRAIL_N: the first ring after a regime change is a settling
        // allowance (overwritten), and the second ring becomes the portrait.
        if (ptsSinceGen < 2 * TRAIL_N) ptsSinceGen++;
    }

    void process(const ProcessArgs& args) override {
        const float fs = args.sampleRate;

        // ── read + map params (sanitize CV before use) ──
        float balCv  = inputs[BALANCE_INPUT].getVoltage();
        float wildCv = inputs[WILD_INPUT].getVoltage();
        if (!std::isfinite(balCv))  balCv = 0.f;
        if (!std::isfinite(wildCv)) wildCv = 0.f;
        float balance = clamp(params[BALANCE_PARAM].getValue()
            + balCv * params[BALANCE_ATT_PARAM].getValue() * BALANCE_CV_DEPTH, 0.f, 1.f);
        float wild = clamp(params[WILD_PARAM].getValue()
            + wildCv * params[WILD_ATT_PARAM].getValue() * WILD_CV_DEPTH, 0.f, 1.f);
        float kickV = inputs[KICK_INPUT].isConnected() ? inputs[KICK_INPUT].getVoltage() : 0.f;
        if (!std::isfinite(kickV)) kickV = 0.f;
        kickV = clamp(kickV, -15.f, 15.f);   // bound an extreme finite drive, not just non-finite
        float kick = kickV * KICK_GAIN;

        float b1 = 1.f + 5.2f * wild * wild;   // WILD → bifurcation parameter
        // BALANCE → upper trophic saturation, narrowed to [1.75, 2.25] (b2=2.0 at
        // centre). Above b2≈2.35 a fox-extinction basin captures the fixed seed and
        // rails FOX at −5 V, so the range is deliberately kept clear of it.
        float b2 = 1.75f + 0.5f * balance;

        // ── center (analytic); fall back to the last valid center on corrupt CV ──
        if (!updateCenter(b1, b2)) { reseed(); resetPeakMemory(); resetTrail(); }

        // ── display regime signature ──
        // Bump the regime generation when the portrait should refresh — a real change
        // in b1/b2 or a change in the sustained KICK level. ptsSinceGen (reset here)
        // then requires a FULL fresh trajectory before the ring counts as a complete
        // portrait, so mid-turn/half-old traces never reach the screen, and continuous
        // modulation simply holds the last settled portrait.
        if (std::fabs(b1 - regB1) > 0.01f || std::fabs(b2 - regB2) > 0.01f || std::fabs(kick - regKick) > 0.05f) {
            dispGen++; ptsSinceGen = 0; regB1 = b1; regB2 = b2; regKick = kick;
        }

        // ── pitch → dimensionless-time step, adaptive substepping ──
        // Sanitize V/OCT BEFORE the clamp: clamp() is fmin/fmax-based and maps a NaN
        // to a bound, so a corrupt V/OCT would otherwise pass the finite check and
        // drive full speed instead of resting.
        float voct = inputs[VOCT_INPUT].getVoltage();
        if (!std::isfinite(voct)) voct = 0.f;
        float pitchTotal = clamp(params[RATE_PARAM].getValue() + voct,
                                 PITCH_TOTAL_MIN, PITCH_TOTAL_MAX);
        float pitchHz = dsp::FREQ_C4 * dsp::approxExp2_taylor5(pitchTotal);
        float dtau = RATE_CAL * pitchHz / fs;
        dtau = std::min(dtau, HSUB_MAX * MAX_SUB);   // documented simulation-speed ceiling
        int   Ksub = clamp((int) std::ceil(dtau / HSUB_MAX), MIN_SUB, MAX_SUB);
        float h = dtau / Ksub;

        auto deriv = [&](const float* v, float* dv) {
            float X = std::max(v[0], POS_FLOOR), Y = std::max(v[1], POS_FLOOR), Z = std::max(v[2], POS_FLOOR);
            float g1 = f1(X, b1), g2 = f2(Y, b2);
            dv[0] = X * (1.f - X) - g1 * Y + kick;   // KICK: continuous grass-force, constant across RK stages
            dv[1] = g1 * Y - g2 * Z - D1 * Y;
            dv[2] = g2 * Z - D2 * Z;
        };

        float st[3] = {x, y, z};
        for (int s = 0; s < Ksub; ++s) {
            float pre[3] = {st[0], st[1], st[2]};   // for capture interpolation
            coalescent::rk4<3>(st, h, deriv);
            // isfinite BEFORE clamp: clamp() is fmin/fmax-based and maps NaN to a
            // bound, which would hide the failure from this backstop.
            if (!std::isfinite(st[0]) || !std::isfinite(st[1]) || !std::isfinite(st[2])) {
                reseed(); st[0] = x; st[1] = y; st[2] = z; resetPeakMemory(); resetTrail(); break;
            }
            st[0] = clamp(st[0], POS_FLOOR, STATE_MAX);
            st[1] = clamp(st[1], POS_FLOOR, STATE_MAX);
            st[2] = clamp(st[2], POS_FLOOR, STATE_MAX);

            // peak events from substep state (so fast grass/bunny peaks aren't missed
            // at the audio-sample boundary), on the gained centered value
            firePeak(0, (st[0] - cx) * FOX_GAIN[0], peakGen[0]);
            firePeak(1, (st[1] - cy) * FOX_GAIN[1], peakGen[1]);
            firePeak(2, (st[2] - cz) * FOX_GAIN[2], peakGen[2]);

            // Trail capture at uniform dimensionless time. Interpolate the point to
            // the exact CAP_TAU boundary inside this substep (h ≤ CAP_TAU, so at most
            // one boundary per substep — bounded) so the spacing is genuinely uniform
            // rather than quantized to the substep, independent of pitch.
            float before = capAccum;
            capAccum += h;
            if (capAccum >= CAP_TAU) {
                float f = clamp((CAP_TAU - before) / h, 0.f, 1.f);   // where in the substep the boundary fell
                capAccum -= CAP_TAU;
                float ix = pre[0] + (st[0] - pre[0]) * f;
                float iy = pre[1] + (st[1] - pre[1]) * f;
                float iz = pre[2] + (st[2] - pre[2]) * f;
                pushTrail((ix - cx) * FOX_GAIN[0], (iy - cy) * FOX_GAIN[1], (iz - cz) * FOX_GAIN[2]);
            }
        }
        x = st[0]; y = st[1]; z = st[2];

        // ── outputs: equilibrium-centered, fixed-gain, soft-clipped ──
        outputs[GRASS_OUTPUT].setVoltage(5.f * coalescent::fastTanh((x - cx) * FOX_GAIN[0]));
        outputs[BUNNY_OUTPUT].setVoltage(5.f * coalescent::fastTanh((y - cy) * FOX_GAIN[1]));
        outputs[FOX_OUTPUT].setVoltage(5.f * coalescent::fastTanh((z - cz) * FOX_GAIN[2]));
        outputs[GRASS_PEAK_OUTPUT].setVoltage(peakGen[0].process(args.sampleTime) ? GATE_LEVEL : 0.f);
        outputs[BUNNY_PEAK_OUTPUT].setVoltage(peakGen[1].process(args.sampleTime) ? GATE_LEVEL : 0.f);
        outputs[FOX_PEAK_OUTPUT].setVoltage(peakGen[2].process(args.sampleTime) ? GATE_LEVEL : 0.f);

        // ── publish the trail at ~45 Hz ──
        // `clean` = the ring is a complete post-settling portrait of regimeGen — one
        // ring discarded as a settling allowance, then one ring filled. Near a
        // bifurcation convergence can take longer, so this means ready/coherent, not a
        // mathematical proof that every transient has vanished. An existing portrait
        // is only replaced by a clean snapshot; initial startup may visibly fill before
        // its first clean portrait is ready.
        if (fs != lastDisplayFs) { lastDisplayFs = fs; pubDiv.setDivision(std::max(1, (int) std::round(fs / 45.f))); }
        if (pubDiv.process()) {
            liveTrail.resetGen = resetGen;
            liveTrail.regimeGen = dispGen;
            liveTrail.clean = (ptsSinceGen >= 2 * TRAIL_N);
            displaySnapshot.writable() = liveTrail;
            displaySnapshot.publish();
        }
    }
};

// ── Panel layout (mm) ──────────────────────────────────────────────────────────
namespace fpl {
    static const float KNOB_Y = 58.f, ATT_Y = 70.f, IN_Y = 84.f, OUT_Y = 100.f, PEAK_Y = 114.f;
    static const float KNOB_X[3] = {13.f, 30.48f, 47.96f};   // RATE, BALANCE, WILD
    static const float IN_X[4]   = {11.f, 25.5f, 40.f, 52.f}; // V/OCT, BAL, WILD, KICK (Bunnies-compatible)
    static const float OUT_X[3]  = {13.f, 30.48f, 47.96f};    // GRASS, BUNNY, FOX
}

// Isometric projection of the chronological three-species trail. Reads only the
// published snapshot; never touches the integrator.
struct FoxTrailView : widget::TransparentWidget {
    Foxes* module = nullptr;

    static constexpr float INV_SQRT2 = 0.70710678f, INV_SQRT6 = 0.40824829f;
    static constexpr double PLAYHEAD_SECS = 7.0;   // one-way sweep time of the slow playhead
    static constexpr float FOX_RADIUS = 4.2f;
    static constexpr float FOX_EXTENT = FOX_RADIUS * 1.6f;

    // Latched trail for the playhead. At audio rate the live ring churns almost
    // completely between frames, so a slow cursor over it would strobe; the dot
    // instead rides a stable copy, swapped only at a sweep endpoint.
    Foxes::Trail latched;
    bool haveLatch = false;
    bool latchedClean = false;   // was the latched snapshot a settled portrait?
    // mm per projected unit. Fixed (not moving) scale tuned so the canonical
    // chaotic teacup (b1=3) just fills the screen; the periodic default then reads
    // as a compact centered loop that grows into the attractor as WILD rises. Deep
    // WILD (b1>3) intentionally overshoots the bezel. Verified offline against the
    // exact isometric projection over both regimes.
    static constexpr float DISP_SCALE = 3.4f;

    // Dark "screen" + bezel (house style — the SVG is just panel + rails + screws).
    void draw(const DrawArgs& args) override {
        nvgBeginPath(args.vg);
        nvgRoundedRect(args.vg, 0, 0, box.size.x, box.size.y, mm2px(2.f));
        nvgFillColor(args.vg, nvgRGB(0x0a, 0x12, 0x0c));
        nvgFill(args.vg);
        nvgStrokeColor(args.vg, nvgRGB(0x2b, 0x40, 0x30));
        nvgStrokeWidth(args.vg, 1.2f);
        nvgStroke(args.vg);
        Widget::draw(args);
    }

    // Fixed soft-limit lens: linear until 0.75·lim (covers the default and chaos
    // range untouched), then a smooth tanh roll-off to lim. Deep-WILD overshoot
    // curves gracefully toward the bezel instead of hard-clamping into a flat line
    // along the edge. Fixed (frame-independent), so it does not make the attractor
    // breathe — it only bends the extreme excursions back on-screen.
    static float softLimit(float u, float lim) {
        float a = 0.75f * lim, au = std::fabs(u);
        if (au <= a) return u;
        float sign = u < 0.f ? -1.f : 1.f;
        return sign * (a + (lim - a) * std::tanh((au - a) / (lim - a)));
    }

    // (nx,ny,nz) fixed-scaled centered coords → screen point (isometric).
    Vec project(const Vec& ctr, const float* p) const {
        float px = (p[0] - p[1]) * INV_SQRT2;
        float py = (p[0] + p[1] - 2.f * p[2]) * INV_SQRT6;
        float s = mm2px(DISP_SCALE);
        // Reserve the fox's largest extent in the shared projection. The trace and
        // playhead therefore remain coincident even at deep-WILD screen edges.
        float ox = softLimit(px * s, std::max(1.f, box.size.x * 0.5f - FOX_EXTENT));
        float oy = softLimit(py * s, std::max(1.f, box.size.y * 0.5f - FOX_EXTENT));
        return Vec(ctr.x + ox, ctr.y - oy);
    }

    // Small fox glyph riding the playhead — matches the Bunnies bunny's minimal
    // single-colour construction: body ellipse, two pointed triangular ears, and a
    // swept bushy tail (the fox cues that a bunny's rounded ears don't).
    void drawFox(NVGcontext* vg, float px, float py, float r, NVGcolor col) {
        nvgFillColor(vg, col);
        // bushy tail, swept down-left
        nvgBeginPath(vg);
        nvgMoveTo(vg, px - r * 0.55f, py - r * 0.15f);
        nvgLineTo(vg, px - r * 1.60f, py + r * 0.35f);
        nvgLineTo(vg, px - r * 0.70f, py + r * 0.90f);
        nvgClosePath(vg); nvgFill(vg);
        // body
        nvgBeginPath(vg); nvgEllipse(vg, px, py, r * 0.85f, r * 0.70f); nvgFill(vg);
        // left pointed ear
        nvgBeginPath(vg);
        nvgMoveTo(vg, px - r * 0.70f, py - r * 0.40f);
        nvgLineTo(vg, px - r * 0.10f, py - r * 0.50f);
        nvgLineTo(vg, px - r * 0.55f, py - r * 1.40f);
        nvgClosePath(vg); nvgFill(vg);
        // right pointed ear
        nvgBeginPath(vg);
        nvgMoveTo(vg, px + r * 0.70f, py - r * 0.40f);
        nvgLineTo(vg, px + r * 0.10f, py - r * 0.50f);
        nvgLineTo(vg, px + r * 0.55f, py - r * 1.40f);
        nvgClosePath(vg); nvgFill(vg);
    }

    void drawLayer(const DrawArgs& args, int layer) override {
        if (layer != 1) return;
        Vec ctr = box.size.div(2);
        nvgScissor(args.vg, 0, 0, box.size.x, box.size.y);   // keep the trail/fox inside the bezel

        // origin marker (the coexistence equilibrium) — faint grass green
        nvgBeginPath(args.vg); nvgCircle(args.vg, ctr.x, ctr.y, 1.5f);
        nvgFillColor(args.vg, nvgRGBA(0x8a, 0xd0, 0x6a, 0x55)); nvgFill(args.vg);

        static const Foxes::Trail dummy;
        const Foxes::Trail& tr = module ? module->displaySnapshot.consume() : dummy;

        // The live ring churns almost completely between frames at audio rate, so a
        // slow playhead over it would strobe and would ride a different dataset than
        // the drawn trail. Instead LATCH one stable copy and draw BOTH the trail and
        // the fox from it (a stable attractor *portrait*, not a live scope). The latch
        // is swapped only for a good reason: initial fill, an explicit reset, or a
        // complete fresh portrait after recovery or a new regime (tr.clean). Initialize
        // clears the old portrait immediately; ordinary control changes keep it visible
        // while a replacement is collected. A completed swap is deferred toward
        // a sweep endpoint (fox near-stationary)
        // to soften it, but it is a deliberate portrait UPDATE, not seamless: the new
        // ring's endpoint is unrelated to the old, so the trace and fox do step to the
        // new shape. A half-old mid-turn ring is never clean, so it is never shown.
        double swept = system::getTime() / PLAYHEAD_SECS;
        double triNow = std::fmod(swept, 2.0); triNow = triNow < 1.0 ? triNow : 2.0 - triNow;
        bool atEndpoint = triNow < 0.03 || triNow > 0.97;
        // Swap to a settled portrait (tr.clean) once one is available for a NEW regime,
        // or once the current latch — which may have been an unclean fill/transient —
        // is superseded by a clean one of the same regime. Deferred to a sweep endpoint.
        const bool resetChanged = tr.resetGen != latched.resetGen;
        bool freshPortrait = tr.clean && (resetChanged
            || tr.regimeGen != latched.regimeGen || !latchedClean) && atEndpoint;
        if (!haveLatch || resetChanged || latched.count < TRAIL_N || freshPortrait) {
            latched = tr; latchedClean = tr.clean; haveLatch = true;
        }

        if (latched.count >= 2) {
            // Chronological trail, oldest → newest, fading in a bounded number of
            // alpha bands (no per-segment stroke, no closing chord across the ring seam).
            const int BANDS = 32;
            int start = (latched.head - latched.count + TRAIL_N) % TRAIL_N;
            nvgLineCap(args.vg, NVG_ROUND); nvgLineJoin(args.vg, NVG_ROUND);
            for (int band = 0; band < BANDS; ++band) {
                int lo = (int)((int64_t) band * (latched.count - 1) / BANDS);
                int hi = (int)((int64_t)(band + 1) * (latched.count - 1) / BANDS);
                if (hi <= lo) continue;
                float age = (float) band / BANDS;                 // 0 oldest … 1 newest
                int alpha = (int) (30 + 200 * age * age);
                nvgBeginPath(args.vg);
                for (int i = lo; i <= hi; ++i) {
                    Vec p = project(ctr, latched.pt[(start + i) % TRAIL_N]);
                    if (i == lo) nvgMoveTo(args.vg, p.x, p.y); else nvgLineTo(args.vg, p.x, p.y);
                }
                nvgStrokeColor(args.vg, nvgRGBA(0xf0, 0xa8, 0x3c, alpha));
                nvgStrokeWidth(args.vg, 1.3f); nvgStroke(args.vg);
            }
            // Slow playhead: sweep the fox along the SAME latched trail there-and-back
            // at a constant, watchable pace (reversing at the ends → no teleport seam).
            float fidx = (float) triNow * (latched.count - 1);
            int i0 = (int) fidx, i1 = std::min(i0 + 1, latched.count - 1);
            float fr = fidx - i0;
            Vec p0 = project(ctr, latched.pt[(start + i0) % TRAIL_N]);
            Vec p1 = project(ctr, latched.pt[(start + i1) % TRAIL_N]);
            Vec hp = p0.plus(p1.minus(p0).mult(fr));
            drawFox(args.vg, hp.x, hp.y, FOX_RADIUS, nvgRGB(0xff, 0xb4, 0x54));
        }

        // Title (top-left) + hint (bottom-right), DejaVuSans, amber accent.
        std::shared_ptr<Font> font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
        if (font) {
            nvgFontFaceId(args.vg, font->handle);
            nvgFontSize(args.vg, mm2px(3.2f));
            nvgTextLetterSpacing(args.vg, mm2px(0.5f));
            nvgFillColor(args.vg, nvgRGBA(0xf0, 0xc0, 0x78, 0xcc));
            nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
            nvgText(args.vg, mm2px(2.6f), mm2px(2.2f), "FOXES", NULL);
            nvgTextLetterSpacing(args.vg, 0.f);
            nvgFontSize(args.vg, mm2px(2.2f));
            nvgFillColor(args.vg, nvgRGBA(0xa8, 0xc0, 0x88, 0xaa));
            nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_BOTTOM);
            nvgText(args.vg, box.size.x - mm2px(2.4f), box.size.y - mm2px(1.8f), "food chain", NULL);
        }
    }
};

// Control labels (nanosvg ignores <text>, so labels live here on the panel).
struct FoxesLabels : widget::Widget {
    void draw(const DrawArgs& args) override {
        std::shared_ptr<Font> font = APP->window->loadFont(asset::system("res/fonts/Nunito-Bold.ttf"));
        if (!font) return;
        nvgFontFaceId(args.vg, font->handle);
        nvgFillColor(args.vg, nvgRGB(0xe6, 0xf2, 0xe6));   // near-white, house legibility spec
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
        auto label = [&](float xmm, float ymm, const char* s, float sz = 1.9f) {
            nvgFontSize(args.vg, mm2px(sz * 1.72f));
            nvgText(args.vg, mm2px(xmm), mm2px(ymm), s, nullptr);
        };
        using namespace fpl;
        const char* kl[3] = {"RATE", "BALANCE", "WILD"};
        for (int i = 0; i < 3; ++i) label(KNOB_X[i], KNOB_Y - 8.f, kl[i], 1.7f);
        const char* il[4] = {"V/OCT", "BAL", "WILD", "KICK"};
        for (int i = 0; i < 4; ++i) label(IN_X[i], IN_Y - 5.5f, il[i], 1.55f);
        const char* ol[3] = {"GRASS", "BUNNY", "FOX"};
        for (int i = 0; i < 3; ++i) label(OUT_X[i], OUT_Y - 5.5f, ol[i], 1.55f);
        const char* pl[3] = {"G↑", "B↑", "F↑"};   // peak events (↑ distinguishes them from the population outputs)
        for (int i = 0; i < 3; ++i) label(OUT_X[i], PEAK_Y - 5.5f, pl[i], 1.5f);
    }
};

struct FoxesWidget : ModuleWidget {

    FoxesWidget(Foxes* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Foxes.svg")));
        addPanelLabels<FoxesLabels>(this);
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f, 122.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 122.0f))));

        FoxTrailView* tv = new FoxTrailView();
        tv->module = module;
        tv->box.pos = mm2px(Vec(5.f, 8.f));
        tv->box.size = mm2px(Vec(50.96f, 38.f));
        addChild(tv);

        using namespace fpl;
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(KNOB_X[0], KNOB_Y)), module, Foxes::RATE_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(KNOB_X[1], KNOB_Y)), module, Foxes::BALANCE_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(KNOB_X[2], KNOB_Y)), module, Foxes::WILD_PARAM));
        addParam(createParamCentered<Trimpot>(mm2px(Vec(KNOB_X[1], ATT_Y)), module, Foxes::BALANCE_ATT_PARAM));
        addParam(createParamCentered<Trimpot>(mm2px(Vec(KNOB_X[2], ATT_Y)), module, Foxes::WILD_ATT_PARAM));

        const int inId[4] = {Foxes::VOCT_INPUT, Foxes::BALANCE_INPUT, Foxes::WILD_INPUT, Foxes::KICK_INPUT};
        for (int i = 0; i < 4; ++i)
            addInput(createInputCentered<PJ301MPort>(mm2px(Vec(IN_X[i], IN_Y)), module, inId[i]));
        const int outId[3] = {Foxes::GRASS_OUTPUT, Foxes::BUNNY_OUTPUT, Foxes::FOX_OUTPUT};
        const int popId[3] = {Foxes::GRASS_PEAK_OUTPUT, Foxes::BUNNY_PEAK_OUTPUT, Foxes::FOX_PEAK_OUTPUT};
        for (int i = 0; i < 3; ++i) {
            addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(OUT_X[i], OUT_Y)), module, outId[i]));
            addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(OUT_X[i], PEAK_Y)), module, popId[i]));
        }
    }
};

Model* modelFoxes = createModel<Foxes, FoxesWidget>("Foxes");
