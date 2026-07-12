#include "plugin.hpp"
#include "../dsp/rk4.hpp"
#include "../dsp/display_snapshot.hpp"
#include "../dsp/completed_path.hpp"
#include "tanh_approx.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>

// AXON — a spiking-neuron oscillator on the FitzHugh-Nagumo (FHN) system.
//
// The membrane voltage v is the audio output; pitch is the *speed* the
// dimensionless simulation is run at (open-loop calibrated, not phase-locked —
// the limit-cycle period is emergent, so I/eps/shape pull pitch a little; that
// is part of the instrument). Above a current threshold it free-runs as a
// relaxation oscillator; below threshold it rests and fires one spike per
// trigger.
//
// The integrator is the whole engineering job: RK4 with pitch-adaptive
// substepping so the stiff relaxation spike stays accurate/stable as pitch
// rises. f() and the RK4 step are factored so a Hindmarsh-Rose sibling (v2)
// can reuse them with one extra state variable.
//
// Polyphonic: up to 16 independent FHN voices. The voice count follows V/OCT
// (falling back to TRIG, so a poly trigger drives poly percussion even with no
// pitch cable). The display traces every active voice, each on its own hue
// stepped across a narrow band around the module's accent colour.

// Length of the (v,w) phase-portrait trail shown on the display.
static const int TRAIL = 512;
// Screen-space bins for the latest completed orbit the slow play head rides.
static const int   ORBIT_N    = 128;
static const float ORBIT_SECS = 7.f;   // wall-clock seconds per lap

struct Axon : Module {

    enum ParamId {
        PITCH_PARAM,
        CURRENT_PARAM,
        EPS_PARAM,
        SHAPE_PARAM,
        CURRENT_ATT_PARAM,
        EPS_ATT_PARAM,
        PARAMS_LEN
    };

    enum InputId {
        VOCT_INPUT,
        CURRENT_INPUT,
        EPS_INPUT,
        TRIG_INPUT,
        SYNC_INPUT,
        INPUTS_LEN
    };

    enum OutputId {
        OUT_OUTPUT,
        SPIKE_OUTPUT,
        W_OUTPUT,
        OUTPUTS_LEN
    };

    enum LightId {
        LIGHTS_LEN
    };

    // ─── Per-voice persistent state (polyphonic, up to 16 channels) ─────────
    // Each voice is an independent FHN cell, so the integration state and the
    // stateful DSP helpers (TRIG edge detect, spike detect/shaper, DC blocker)
    // are one-per-voice. The voice count is driven by V/OCT (or TRIG).
    static const int MAX_POLY = 16;
    static const int GROUPS   = MAX_POLY / 4;   // voices processed 4-wide (simd::float_4)
    // Poly voices run four at a time. Each group's integration state and stateful
    // DSP helpers hold one lane per voice; the whole oversampled output chain
    // (RK4 → clamp → spike detect → DC block → tanh → decimate) is vectorised.
    simd::float_4 vv4[GROUPS], ww4[GROUPS];      // FHN state per voice
    simd::float_4 trigPulse4[GROUPS] = {};       // decaying injected current from TRIG
    dsp::TSchmittTrigger<simd::float_4> trigIn4[GROUPS];   // TRIG edge detector (per lane)
    dsp::TSchmittTrigger<simd::float_4> syncIn4[GROUPS];   // SYNC edge detector (per lane)
    dsp::TSchmittTrigger<simd::float_4> spikeDet4[GROUPS]; // v crossing SPIKE_THRESH upward
    dsp::PulseGenerator spikeGen[MAX_POLY];      // SPIKE pulse shaper (scalar, one per voice)
    dsp::TRCFilter<simd::float_4> dcBlock4[GROUPS]; // DC blocker on OUT (limit-cycle mean ≠ 0)
    // Anti-aliasing decimators (one per group; only the active factor is used).
    dsp::Decimator<2, 8, simd::float_4> decim2_4[GROUPS];
    dsp::Decimator<4, 8, simd::float_4> decim4_4[GROUPS];
    dsp::Decimator<8, 8, simd::float_4> decim8_4[GROUPS];
    int   oversampleMode = 2;         // 0=off, 1=×2, 2=×4, 3=×8
    int   channels = 1;               // active voice count
    float lastFs = 0.f;               // detects SR change to refresh coefficients
    int   lastOs = 0;                 // detects oversample change to refresh coefficients
    float trigDecay = 0.f;            // cached per-sample TRIG pulse decay (fn of fs)
    float antiDenorm = 1e-18f;        // alternating sub-LSB dither, anti-denormal on dcBlock

    // ─── Display trail (phase portrait) ─────────────────────────────────────
    // The audio thread appends decimated (v,w) points to a per-voice ring, then
    // ~45 Hz publishes a coherent snapshot into a lock-free triple buffer; the UI
    // reads the latest complete frame. The head index, active-channel count and
    // effective CURRENT travel *with* the arrays, so the trails, head dots and
    // nullcline stay mutually consistent. See dsp/display_snapshot.hpp.
    float trailV[MAX_POLY][TRAIL] = {}, trailW[MAX_POLY][TRAIL] = {};
    int   trailIdx = 0, trailDecim = 0;
    // Audio-rate capture prevents the decimated trail from changing which side
    // of Axon's narrow spike a phase bin lands on from cycle to cycle. Completed
    // paths are promoted at most once per display frame; intermediate cycles are
    // cheap to capture and discard.
    using OrbitPath = coalescent::CompletedPath<ORBIT_N>;
    using OrbitPoint = OrbitPath::Point;
    OrbitPath orbitPath[MAX_POLY];
    int orbitCommitClock[MAX_POLY] = {};
    struct DisplayFrame {
        float v[MAX_POLY][TRAIL] = {}, w[MAX_POLY][TRAIL] = {};
        OrbitPoint orbit[MAX_POLY][ORBIT_N] = {}; // latest completed orbit the dot rides
        int   orbitCount[MAX_POLY] = {};
        bool  orbitValid[MAX_POLY] = {};
        bool  orbitClosed[MAX_POLY] = {};
        int   head = 0;
        int   channels = 1;
        float current = 0.6f;           // voice-0 effective CURRENT for the v-nullcline
    };
    coalescent::DisplaySnapshot<DisplayFrame> displaySnapshot;
    int   dispClock = 0;
    int   dispPeriod = 980;           // samples between snapshots (~45 Hz; refreshed on SR change)

    // ─── Tunable constants (set by ear/scope; RATE_CAL from tools/stability/) ─
    static constexpr float RATE_CAL    = 37.899004f; // dimensionless period at default → C4 at 0 V
    static constexpr float B_FIXED     = 0.8f;       // recovery linear coefficient b
    static constexpr float HSUB_MAX    = 0.05f;      // max dimensionless substep size
    static constexpr int   MIN_SUB     = 2;   // substep floor: 2 holds 0-cent accuracy at default (half the RK4 cost of a 4-floor)
    static constexpr int   MAX_SUB     = 64;
    static constexpr float TRIG_AMP    = 0.6f;       // injected current pulse height
    static constexpr float TRIG_TAU_MS = 3.f;        // pulse decay time constant (ms)
    static constexpr float SPIKE_THRESH = 1.0f;      // v level for a SPIKE pulse
    static constexpr float OUT_GAIN    = 1.0f;
    static constexpr float W_GAIN      = 2.5f;
    static constexpr float CV_DEPTH    = 0.1f;       // ±5V CV × att → ±0.5
    static constexpr float STATE_MAX   = 10.f;       // safety clamp on v,w

    Axon() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);

        configParam(PITCH_PARAM, -4.f, 4.f, 0.f, "Pitch", " Hz", 2.f, dsp::FREQ_C4);
        configParam(CURRENT_PARAM, -0.2f, 1.6f, 0.6f, "Current (excitability)");
        configParam(EPS_PARAM, 0.01f, 0.30f, 0.08f, "Timescale ε (spike sharpness)");
        configParam(SHAPE_PARAM, 0.4f, 1.0f, 0.7f, "Shape (asymmetry a)");
        configParam(CURRENT_ATT_PARAM, -1.f, 1.f, 0.f, "Current CV");
        configParam(EPS_ATT_PARAM, -1.f, 1.f, 0.f, "ε CV");

        configInput(VOCT_INPUT, "1V/oct pitch");
        configInput(CURRENT_INPUT, "Current CV");
        configInput(EPS_INPUT, "ε CV");
        configInput(TRIG_INPUT, "Trigger (fire a spike)");
        configInput(SYNC_INPUT, "Sync (hard reset to rest)");

        configOutput(OUT_OUTPUT, "Audio (membrane voltage v)");
        configOutput(SPIKE_OUTPUT, "Spike (trigger on each spike)");
        configOutput(W_OUTPUT, "Recovery w (slow CV)");

        for (int g = 0; g < GROUPS; g++) { vv4[g] = -1.2f; ww4[g] = -0.6f; }
        // A sub-pixel spacing stagger keeps identical poly voices from running
        // their bounded guide simplification on the same audio callback.
        for (int c = 0; c < MAX_POLY; ++c) {
            orbitPath[c].setMetric(50.f / 5.8f, 34.f / 3.8f);
            orbitPath[c].setMarkerAges(0.5f, 1024.f);
            orbitPath[c].setMinimumSpacing(0.15f + 0.0015f * c);
            orbitCommitClock[c] = 1 << 24;
        }
    }

    void onReset() override {
        oversampleMode = 2;   // restore default anti-aliasing (×4) on Initialize
        for (int g = 0; g < GROUPS; g++) {
            vv4[g] = -1.2f; ww4[g] = -0.6f; trigPulse4[g] = 0.f;
            trigIn4[g].reset(); syncIn4[g].reset(); spikeDet4[g].reset();
            dcBlock4[g].reset(); decim2_4[g].reset(); decim4_4[g].reset(); decim8_4[g].reset();
        }
        for (int c = 0; c < MAX_POLY; c++) {
            spikeGen[c].reset();
            orbitPath[c].reset();
            orbitCommitClock[c] = 1 << 24;
        }
    }

    json_t* dataToJson() override {
        json_t* root = json_object();
        // Save the actual factor (1/2/4/8) — immune to menu reordering. Keep the
        // legacy index key too so a 2.0.3 patch opened in 2.0.2 lands close.
        static const int factors[4] = {1, 2, 4, 8};
        json_object_set_new(root, "osFactor",
                            json_integer(factors[clamp(oversampleMode, 0, 3)]));
        json_object_set_new(root, "oversample", json_integer(oversampleMode));
        return root;
    }

    void dataFromJson(json_t* root) override {
        if (json_t* j = json_object_get(root, "osFactor")) {
            // 2.0.3+ patch: the factor is authoritative.
            switch ((int) json_integer_value(j)) {
                case 1:  oversampleMode = 0; break;
                case 2:  oversampleMode = 1; break;
                case 8:  oversampleMode = 3; break;
                default: oversampleMode = 2; break;   // ×4 (and anything unexpected)
            }
        }
        else if (json_t* j = json_object_get(root, "oversample")) {
            // Pre-2.0.3 patch: index meant 0=Off, 1=×4, 2=×8 — remap into the
            // new Off/×2/×4/×8 menu so old patches keep their anti-aliasing.
            int v = (int) json_integer_value(j);
            oversampleMode = v <= 0 ? 0 : (v == 1 ? 2 : 3);
        }
    }

    // FHN derivatives in dimensionless time, generic over the scalar type T so the
    // one expression serves scalar float and simd::float_4 (four poly voices at
    // once); the shared coalescent::rk4 step (integrator.hpp) is likewise templated.
    template <typename T>
    static inline void f(T v, T w, T Itot, T eps, T a, T& dv, T& dw) {
        dv = v - v * v * v / 3.f - w + Itot;
        dw = eps * (v + a - B_FIXED * w);
    }

    void process(const ProcessArgs& args) override {
        const float fs = args.sampleRate;
        const int os = oversampleMode == 0 ? 1 : (oversampleMode == 1 ? 2 : (oversampleMode == 2 ? 4 : 8));

        // SR/oversample-derived coefficients, recomputed only when fs or the
        // oversample factor changes. The DC-blocker runs in the oversampled
        // domain, so its cutoff is referenced to the oversampled rate (fs·os) to
        // keep the real-time corner at ~20 Hz. Caching trigDecay keeps a
        // transcendental off the per-sample path.
        if (fs != lastFs || os != lastOs) {
            for (int g = 0; g < GROUPS; g++) dcBlock4[g].setCutoffFreq(20.f / (fs * os));
            trigDecay = std::exp(-1.f / (TRIG_TAU_MS * 1e-3f * fs));
            dispPeriod = (int) (fs / 45.f);
            if (os != lastOs)
                for (int g = 0; g < GROUPS; g++) { decim2_4[g].reset(); decim4_4[g].reset(); decim8_4[g].reset(); }
            lastFs = fs;
            lastOs = os;
        }

        // Voice count follows V/OCT, but fall back to TRIG so a poly trigger
        // (percussion, no pitch cable) still spreads across voices.
        channels = std::max({inputs[VOCT_INPUT].getChannels(),
                             inputs[TRIG_INPUT].getChannels(), 1});
        outputs[OUT_OUTPUT].setChannels(channels);
        outputs[SPIKE_OUTPUT].setChannels(channels);
        outputs[W_OUTPUT].setChannels(channels);

        // Shared knob values; per-voice CV is added inside the loop.
        const float a         = params[SHAPE_PARAM].getValue();   // no CV by design
        const float epsBase   = params[EPS_PARAM].getValue();
        const float epsAtt    = params[EPS_ATT_PARAM].getValue();
        const float Ibase     = params[CURRENT_PARAM].getValue();
        const float Iatt      = params[CURRENT_ATT_PARAM].getValue();
        const float pitchKnob = params[PITCH_PARAM].getValue();

        // Alternating sub-LSB dither (anti-denormal), toggled once per sample;
        // applied to every voice's DC blocker so a parked rest can't denormalise.
        antiDenorm = -antiDenorm;

        float I0 = Ibase;               // voice-0 effective current, for the display nullcline
        const simd::float_4 aVec(a);    // SHAPE has no CV → same for every lane
        const int nGroups = (channels + 3) / 4;

        for (int g = 0; g < nGroups; g++) {
            const int base = g * 4;

            // Lanes at/above the channel count are inactive: SIMD lanes can't be
            // skipped, so they still compute, but their persistent state must not
            // evolve — preserving the documented "silent voices keep their state"
            // semantic from the scalar implementation. (Edge-detector state isn't
            // masked; it has no musical position to preserve.)
            const simd::float_4 activeMask =
                simd::float_4(0.f, 1.f, 2.f, 3.f) < (float) (channels - base);

            // ── params → physics (4 lanes) ──
            simd::float_4 eps = simd::clamp(epsBase
                + inputs[EPS_INPUT].getPolyVoltageSimd<simd::float_4>(base) * epsAtt * CV_DEPTH,
                0.01f, 0.30f);
            simd::float_4 I = Ibase
                + inputs[CURRENT_INPUT].getPolyVoltageSimd<simd::float_4>(base) * Iatt * CV_DEPTH;
            if (g == 0) I0 = I[0];

            simd::float_4 v = vv4[g], w = ww4[g];

            // ── hard sync: a rising edge re-seeds those lanes at the rest point ──
            simd::float_4 syncMask = syncIn4[g].process(
                inputs[SYNC_INPUT].getPolyVoltageSimd<simd::float_4>(base), 0.1f, 1.f);
            const int syncBits = simd::movemask(syncMask);
            v = simd::ifelse(syncMask, -1.2f, v);
            w = simd::ifelse(syncMask, -0.6f, w);

            // ── excitable trigger: rising edge injects a decaying current pulse ──
            simd::float_4 trigMask = trigIn4[g].process(
                inputs[TRIG_INPUT].getPolyVoltageSimd<simd::float_4>(base), 0.1f, 1.f);
            simd::float_4 tp = simd::ifelse(trigMask, TRIG_AMP, trigPulse4[g]);
            simd::float_4 Itot = I + tp;
            tp *= trigDecay;
            tp = simd::ifelse(simd::abs(tp) < 1e-30f, 0.f, tp);
            trigPulse4[g] = simd::ifelse(activeMask, tp, trigPulse4[g]);

            // ── pitch = simulation speed (open-loop). Run the group at its max
            // substep count K: a slower lane just gets a smaller-than-needed h,
            // so accuracy is preserved lane-by-lane (same reasoning as MIN_SUB). ──
            // Sanitize the pitch exponent BEFORE approxExp2: its internal float->int
            // shift is UB on a non-finite or absurd CV (NaN->0=C4; clamp bounds ±inf
            // and huge-finite). The subTau cap below then bounds the resulting step.
            simd::float_4 pexp = pitchKnob + inputs[VOCT_INPUT].getPolyVoltageSimd<simd::float_4>(base);
            pexp = simd::clamp(simd::ifelse(pexp == pexp, pexp, 0.f), -30.f, 30.f);
            simd::float_4 pitchHz = dsp::FREQ_C4 * dsp::approxExp2_taylor5(pexp);
            simd::float_4 subTau = RATE_CAL * pitchHz / fs / (float) os;
            // Hard guard (matches Operon/Bunnies): cap the per-sample sim-time so
            // h = subTau/K stays <= HSUB_MAX and Kf stays in range, even at extreme
            // V/OCT (an uncapped Kf can overflow the float->int below = UB). The
            // trade is that pitch goes flat above the cap with oversampling Off.
            subTau = simd::fmin(subTau, simd::float_4(HSUB_MAX * MAX_SUB));
            simd::float_4 Kf = simd::ceil(subTau / HSUB_MAX);
            int K = MIN_SUB;
            for (int l = 0; l < 4; l++) K = std::max(K, (int) Kf[l]);
            K = std::min(K, MAX_SUB);
            simd::float_4 h = subTau / (float) K;

            // ── render os oversampled samples, full chain each, then decimate ──
            simd::float_4 s[2] = {v, w};
            simd::float_4 osBuf[8];
            for (int l = 0; l < 4 && base + l < channels; ++l) {
                const int c = base + l;
                orbitCommitClock[c] = std::min(orbitCommitClock[c] + 1, dispPeriod);
            }
            for (int o = 0; o < os; o++) {
                for (int k = 0; k < K; k++)
                    coalescent::rk4<2>(s, h, [&](const simd::float_4* y, simd::float_4* d) {
                        f(y[0], y[1], Itot, eps, aVec, d[0], d[1]);
                    });

                // Backstop: reset any non-finite / runaway lane to rest, then clamp.
                simd::float_4 finite = (s[0] == s[0]) & (s[1] == s[1])
                                     & (simd::abs(s[0]) < 1e6f) & (simd::abs(s[1]) < 1e6f);
                const int orbitResetBits = (o == 0 ? syncBits : 0)
                    | ((~simd::movemask(finite)) & 0xf);
                s[0] = simd::ifelse(finite, s[0], -1.2f);
                s[1] = simd::ifelse(finite, s[1], -0.6f);
                s[0] = simd::clamp(s[0], -STATE_MAX, STATE_MAX);
                s[1] = simd::clamp(s[1], -STATE_MAX, STATE_MAX);

                // Spike detect per lane; fire the matching scalar pulse generators
                // (active lanes only — hidden lanes must not queue pulses).
                simd::float_4 spikeMask = spikeDet4[g].process(s[0], 0.f, SPIKE_THRESH);
                int sm = simd::movemask(spikeMask);
                for (int l = 0; l < 4 && base + l < channels; l++)
                    if (sm & (1 << l)) spikeGen[base + l].trigger(1e-3f);

                dcBlock4[g].process(s[0] + antiDenorm);
                osBuf[o] = 5.f * coalescent::fastTanh(dcBlock4[g].highpass() * OUT_GAIN);

                // Sample the guide at the actual oversampled state rate. This is
                // separate from the fading trail and resolves the narrow spike at
                // high pitch without publishing partial paths to the UI.
                for (int l = 0; l < 4 && base + l < channels; ++l) {
                    const int c = base + l;
                    const bool committed = orbitPath[c].push(
                        s[0][l], s[1][l], s[0][l], subTau[l],
                        orbitCommitClock[c] >= dispPeriod, orbitResetBits & (1 << l));
                    if (committed)
                        orbitCommitClock[c] = 0;
                }
            }
            // Write back active lanes only; inactive lanes keep their frozen state.
            vv4[g] = simd::ifelse(activeMask, s[0], vv4[g]);
            ww4[g] = simd::ifelse(activeMask, s[1], ww4[g]);

            simd::float_4 audio = os == 1 ? osBuf[0]
                : (os == 2 ? decim2_4[g].process(osBuf)
                :  os == 4 ? decim4_4[g].process(osBuf) : decim8_4[g].process(osBuf));

            // ── outputs (4 channels at once) ──
            outputs[OUT_OUTPUT].setVoltageSimd(audio, base);
            outputs[W_OUTPUT].setVoltageSimd(simd::clamp(s[1] * W_GAIN, -5.f, 5.f), base);
            simd::float_4 spikeOut;
            for (int l = 0; l < 4; l++)
                spikeOut[l] = spikeGen[base + l].process(args.sampleTime) ? 10.f : 0.f;
            outputs[SPIKE_OUTPUT].setVoltageSimd(spikeOut, base);
        }

        // ── feed the phase-portrait trails (all active voices, decimated) ──
        if (++trailDecim >= 4) {
            trailDecim = 0;
            for (int c = 0; c < channels; c++) {
                trailV[c][trailIdx] = vv4[c / 4][c % 4];
                trailW[c][trailIdx] = ww4[c / 4][c % 4];
            }
            trailIdx = (trailIdx + 1) % TRAIL;
        }
        if (++dispClock >= dispPeriod) {                 // publish display snapshot ~45 Hz
            dispClock = 0;
            DisplayFrame& fr = displaySnapshot.writable();
            for (int c = 0; c < channels; c++) {
                std::copy(trailV[c], trailV[c] + TRAIL, fr.v[c]);
                std::copy(trailW[c], trailW[c] + TRAIL, fr.w[c]);
                std::copy(orbitPath[c].path(), orbitPath[c].path() + ORBIT_N, fr.orbit[c]);
                fr.orbitCount[c] = orbitPath[c].pathSize();
                fr.orbitValid[c] = orbitPath[c].hasPath();
                fr.orbitClosed[c] = orbitPath[c].isClosed();
            }
            fr.head = trailIdx;
            fr.channels = channels;
            fr.current = I0;
            displaySnapshot.publish();
        }
    }
};


// ─── Phase-portrait display ────────────────────────────────────────────────
// The FHN limit cycle traced in the (v,w) plane: a glowing closed orbit when
// free-running, a point that jumps out and relaxes back on each trigger when
// excitable. The v-nullcline (the cubic) and w-nullcline (the line) are drawn
// as faint guides — their intersection is the fixed point whose stability the
// CURRENT knob controls (the Hopf bifurcation). State is read lock-free from
// the audio thread's snapshot. With several voices patched, each orbit is drawn
// on its own hue stepped across a narrow cyan band so the chord is legible.
struct PhaseDisplay : Widget {
    Axon* module = nullptr;
    // Adopt completed geometry only at a visual-lap boundary. Modulated cycles
    // can differ legitimately; replacing one under an unchanged phase would jump.
    Axon::OrbitPoint playOrbit[Axon::MAX_POLY][ORBIT_N] = {};
    int playCount[Axon::MAX_POLY] = {};
    long long playLap[Axon::MAX_POLY] = {};
    bool playValid[Axon::MAX_POLY] = {};
    bool playClosed[Axon::MAX_POLY] = {};

    // Plot ranges for v (x axis) and w (y axis). Sized with margin around the
    // orbit envelope (v∈[-2.0,2.0], w∈[-0.6,2.3] across the oscillating param
    // space) so the limit cycle sits comfortably inside the screen.
    static constexpr float VMIN = -2.9f, VMAX = 2.9f;
    static constexpr float WMIN = -1.1f, WMAX = 2.7f;

    // Voice → hue (0..1) within a band centred on Axon's cyan accent (~0.554).
    static float voiceHue(int v, int nv) {
        const float center = 0.554f, halfBand = 0.072f;
        if (nv <= 1) return center;
        return center - halfBand + 2.f * halfBand * v / (nv - 1);
    }

    void draw(const DrawArgs& args) override {
        // Screen background + bezel (base, non-illuminated layer).
        nvgBeginPath(args.vg);
        nvgRoundedRect(args.vg, 0, 0, box.size.x, box.size.y, mm2px(2.f));
        nvgFillColor(args.vg, nvgRGB(0x07, 0x07, 0x12));
        nvgFill(args.vg);
        nvgStrokeColor(args.vg, nvgRGB(0x2b, 0x2b, 0x4d));
        nvgStrokeWidth(args.vg, 1.2f);
        nvgStroke(args.vg);
        Widget::draw(args);
    }

    void drawLayer(const DrawArgs& args, int layer) override {
        if (layer == 1) {
            // Clip everything to the screen rectangle. The cubic v-nullcline (and
            // an orbit pushed past its envelope by CV / triggers) would otherwise
            // draw outside the box and spill onto the panel.
            nvgScissor(args.vg, 0.f, 0.f, box.size.x, box.size.y);

            auto X = [&](float v) { return (v - VMIN) / (VMAX - VMIN) * box.size.x; };
            auto Y = [&](float w) { return box.size.y - (w - WMIN) / (WMAX - WMIN) * box.size.y; };

            static const Axon::DisplayFrame dummy;
            const Axon::DisplayFrame& fr = module ? module->displaySnapshot.consume() : dummy;
            float a = module ? module->params[Axon::SHAPE_PARAM].getValue() : 0.7f;
            float I = module ? fr.current : 0.6f;   // CV-included, coherent with the trail
            int nv = module ? fr.channels : 1;
            const double playCycles = system::getTime() / (double) ORBIT_SECS;
            const double lapStart = std::floor(playCycles);
            const long long currentLap = (long long) lapStart;
            const float lapFraction = (float) (playCycles - lapStart);
            for (int v = 0; v < nv; ++v) {
                if (!fr.orbitValid[v]) {
                    playValid[v] = false;
                } else if (!playValid[v] || playLap[v] != currentLap) {
                    playCount[v] = clamp(fr.orbitCount[v], 2, ORBIT_N);
                    std::copy(fr.orbit[v], fr.orbit[v] + playCount[v], playOrbit[v]);
                    playClosed[v] = fr.orbitClosed[v];
                    playValid[v] = true;
                    playLap[v] = currentLap;
                }
            }

            // v-nullcline: w = v - v³/3 + I  (the cubic). Where dv/dt = 0.
            nvgBeginPath(args.vg);
            for (int i = 0; i <= 80; i++) {
                float v = VMIN + (VMAX - VMIN) * i / 80.f;
                float w = v - v * v * v / 3.f + I;
                float px = X(v), py = Y(w);
                if (i == 0) nvgMoveTo(args.vg, px, py); else nvgLineTo(args.vg, px, py);
            }
            nvgStrokeColor(args.vg, nvgRGBA(0x40, 0x70, 0x90, 0x55));
            nvgStrokeWidth(args.vg, 1.f);
            nvgStroke(args.vg);

            // w-nullcline: w = (v + a)/b  (the line). Where dw/dt = 0.
            nvgBeginPath(args.vg);
            {
                float w0 = (VMIN + a) / Axon::B_FIXED, w1 = (VMAX + a) / Axon::B_FIXED;
                nvgMoveTo(args.vg, X(VMIN), Y(w0));
                nvgLineTo(args.vg, X(VMAX), Y(w1));
            }
            nvgStrokeColor(args.vg, nvgRGBA(0x90, 0x50, 0x70, 0x55));
            nvgStrokeWidth(args.vg, 1.f);
            nvgStroke(args.vg);

            // The trajectory trails, one per active voice, oldest→newest brightening
            // along their length. More voices → dimmer trails so they don't overwhelm.
            // The fade is quantized into TRAIL_BANDS constant-alpha bands: each band is
            // one connected polyline stroked once, instead of a separate stroke per
            // TRAIL segment. At 16 voices that is ~32 strokes/voice rather than ~511 —
            // a ~16× cut in NanoVG calls for a visually indistinguishable (at panel scale) stepped gradient.
            if (module) {
                const int TRAIL_BANDS = 32;
                int idx = fr.head;   // newest just before idx (coherent with arrays)
                float trailA = clamp(204.f - (nv - 1) * 6.f, 112.f, 204.f);  // 0x70..0xcc
                nvgLineCap(args.vg, NVG_ROUND);
                nvgLineJoin(args.vg, NVG_ROUND);
                nvgStrokeWidth(args.vg, 1.6f);
                for (int v = 0; v < nv; v++) {
                    const float* dv = fr.v[v];
                    const float* dw = fr.w[v];
                    float hue = voiceHue(v, nv);
                    const Axon::OrbitPoint* orbit = playOrbit[v];

                    // The fading ring can cover less than one orbit at low pitch.
                    // Keep the completed guide faintly visible underneath it so the
                    // play head never appears to leave the trace it is following.
                    if (playValid[v]) {
                        const int orbitCount = playCount[v];
                        nvgBeginPath(args.vg);
                        nvgMoveTo(args.vg, X(orbit[0].x), Y(orbit[0].y));
                        for (int b = 1; b < orbitCount; ++b)
                            nvgLineTo(args.vg, X(orbit[b].x), Y(orbit[b].y));
                        if (playClosed[v])
                            nvgLineTo(args.vg, X(orbit[0].x), Y(orbit[0].y));
                        const int guideAlpha = clamp(64 - (nv - 1) * 2, 32, 64);
                        nvgStrokeColor(args.vg, nvgHSLA(hue, 0.85f, 0.62f, guideAlpha));
                        nvgStrokeWidth(args.vg, 1.f);
                        nvgStroke(args.vg);
                        nvgStrokeWidth(args.vg, 1.6f);
                    }

                    int curBand = -1;
                    auto strokeBand = [&](int band) {   // apply the band's alpha and flush
                        float alpha = (band + 0.5f) / TRAIL_BANDS;   // older = dimmer
                        nvgStrokeColor(args.vg, nvgHSLA(hue, 0.85f, 0.62f, (int)(alpha * trailA)));
                        nvgStroke(args.vg);
                    };
                    for (int k = 1; k < TRAIL; k++) {
                        int band = (k * TRAIL_BANDS) / TRAIL;
                        int i0 = (idx + k - 1) % TRAIL;
                        int i1 = (idx + k) % TRAIL;
                        if (band != curBand) {          // start a new band's polyline at the seam
                            if (curBand >= 0) strokeBand(curBand);
                            nvgBeginPath(args.vg);
                            nvgMoveTo(args.vg, X(dv[i0]), Y(dw[i0]));
                            curBand = band;
                        }
                        nvgLineTo(args.vg, X(dv[i1]), Y(dw[i1]));
                    }
                    if (curBand >= 0) strokeBand(curBand);
                    // Slow play head over the latest completed real orbit. Arc length is
                    // evaluated from this immutable display snapshot, so the dot moves
                    // steadily without feeding UI timing back into the audio thread.
                    float hx, hy;
                    if (playValid[v]) {
                        const int orbitCount = playCount[v];
                        const Axon::OrbitPoint head = coalescent::sampleCompletedPath(
                            orbit, orbitCount, playClosed[v], lapFraction, X, Y);
                        hx = X(head.x); hy = Y(head.y);
                    } else {
                        const int newest = (idx + TRAIL - 1) % TRAIL;
                        hx = X(dv[newest]); hy = Y(dw[newest]);
                    }
                    nvgBeginPath(args.vg); nvgCircle(args.vg, hx, hy, 4.f);
                    nvgFillColor(args.vg, nvgHSLA(hue, 0.85f, 0.72f, 0x55)); nvgFill(args.vg);
                    nvgBeginPath(args.vg); nvgCircle(args.vg, hx, hy, 2.f);
                    nvgFillColor(args.vg, nvgHSLA(hue, 0.55f, 0.88f, 0xff)); nvgFill(args.vg);
                }
            } else {
                // Browser preview: a static demo orbit.
                nvgBeginPath(args.vg);
                for (int i = 0; i <= 64; i++) {
                    float th = 2.f * (float)M_PI * i / 64.f;
                    float px = X(1.4f * std::cos(th)), py = Y(0.4f + 0.7f * std::sin(th));
                    if (i == 0) nvgMoveTo(args.vg, px, py); else nvgLineTo(args.vg, px, py);
                }
                nvgStrokeColor(args.vg, nvgRGBA(0x55, 0xc8, 0xff, 0x99));
                nvgStrokeWidth(args.vg, 1.6f);
                nvgStroke(args.vg);
            }

            // Screen title.
            std::shared_ptr<Font> font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
            if (font) {
                nvgFontFaceId(args.vg, font->handle);
                nvgFontSize(args.vg, mm2px(3.2f));
                nvgTextLetterSpacing(args.vg, mm2px(0.5f));
                nvgFillColor(args.vg, nvgRGBA(0x9a, 0xb0, 0xff, 0xcc));
                nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
                nvgText(args.vg, mm2px(2.6f), mm2px(2.2f), "AXON", NULL);
                nvgTextLetterSpacing(args.vg, 0.f);

                // Voice-count badge when polyphonic.
                if (nv > 1) {
                    char buf[16];
                    snprintf(buf, sizeof(buf), "%dv", nv);
                    nvgFontSize(args.vg, mm2px(2.2f));
                    nvgFillColor(args.vg, nvgRGBA(0x60, 0x80, 0xb0, 0xcc));
                    nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_TOP);
                    nvgText(args.vg, box.size.x - mm2px(2.4f), mm2px(2.2f), buf, NULL);
                }

                nvgFontSize(args.vg, mm2px(2.2f));
                nvgFillColor(args.vg, nvgRGBA(0x60, 0x80, 0xb0, 0xaa));
                nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_BOTTOM);
                nvgText(args.vg, box.size.x - mm2px(2.4f), box.size.y - mm2px(1.8f), "v–w", NULL);
            }

            nvgResetScissor(args.vg);
        }
        Widget::drawLayer(args, layer);
    }
};


// ─── Widget ─────────────────────────────────────────────────────────────────

struct AxonLabels : widget::Widget {

    void draw(const DrawArgs& args) override {
        std::shared_ptr<Font> font = APP->window->loadFont(asset::system("res/fonts/Nunito-Bold.ttf"));
        if (!font) return;

        nvgSave(args.vg);
        nvgFontFaceId(args.vg, font->handle);
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);

        const NVGcolor dim    = nvgRGB(0xe6, 0xe6, 0xf2);  // near-white for low-vision legibility (~13:1)
        const NVGcolor outclr = nvgRGB(0xf4, 0xf4, 0xfe);

        auto lbl = [&](float x, float y, float sz, NVGcolor col, const char* s) {
            nvgFontSize(args.vg, mm2px(sz * 1.72f));   // Nunito Bold, sized up for legibility
            nvgFillColor(args.vg, col);
            nvgText(args.vg, mm2px(x), mm2px(y), s, nullptr);
        };

        // Knob row labels (knobs at y=54)
        lbl( 7.5f,   46.f, 2.1f, dim, "PITCH");
        lbl(22.82f, 46.f, 2.1f, dim, "CURRENT");
        lbl(38.14f, 46.f, 2.1f, dim, "EPS");
        lbl(53.46f, 46.f, 2.1f, dim, "SHAPE");

        // CV jack labels (jacks at y=84)
        lbl(22.82f, 90.f, 2.1f, dim, "I CV");
        lbl(38.14f, 90.f, 2.1f, dim, "ε CV");

        // I/O row labels (row at y=112): three inputs | three outputs
        lbl( 6.5f, 118.5f, 2.1f, dim,    "V/OCT");
        lbl(15.5f, 118.5f, 2.1f, dim,    "TRIG");
        lbl(24.5f, 118.5f, 2.1f, dim,    "SYNC");
        lbl(36.5f, 118.5f, 2.1f, outclr, "OUT");
        lbl(45.5f, 118.5f, 2.1f, outclr, "SPIKE");
        lbl(54.5f, 118.5f, 2.1f, outclr, "W");

        nvgRestore(args.vg);
    }
};

struct AxonWidget : ModuleWidget {

    AxonWidget(Axon* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Axon.svg")));
        addPanelLabels<AxonLabels>(this);

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   122.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 122.0f))));

        // Phase-portrait scope across the top.
        PhaseDisplay* disp = new PhaseDisplay();
        disp->module = module;
        disp->box.pos  = mm2px(Vec(5.5f, 8.f));
        disp->box.size = mm2px(Vec(50.f, 34.f));
        addChild(disp);

        // Knob row (y=54): PITCH | CURRENT | EPS | SHAPE
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec( 7.5f,   54.f)), module, Axon::PITCH_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(22.82f, 54.f)), module, Axon::CURRENT_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(38.14f, 54.f)), module, Axon::EPS_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(53.46f, 54.f)), module, Axon::SHAPE_PARAM));

        // CV strips under CURRENT and EPS: attenuverter (y=72) + input jack (y=84)
        addParam(createParamCentered<Trimpot>(    mm2px(Vec(22.82f, 72.f)), module, Axon::CURRENT_ATT_PARAM));
        addParam(createParamCentered<Trimpot>(    mm2px(Vec(38.14f, 72.f)), module, Axon::EPS_ATT_PARAM));
        addInput(createInputCentered<PJ301MPort>( mm2px(Vec(22.82f, 84.f)), module, Axon::CURRENT_INPUT));
        addInput(createInputCentered<PJ301MPort>( mm2px(Vec(38.14f, 84.f)), module, Axon::EPS_INPUT));

        // I/O row (y=112): V/OCT, TRIG, SYNC (in) | OUT, SPIKE, W (out)
        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec( 6.5f, 112.f)), module, Axon::VOCT_INPUT));
        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec(15.5f, 112.f)), module, Axon::TRIG_INPUT));
        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec(24.5f, 112.f)), module, Axon::SYNC_INPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(36.5f, 112.f)), module, Axon::OUT_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(45.5f, 112.f)), module, Axon::SPIKE_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(54.5f, 112.f)), module, Axon::W_OUTPUT));
    }

    void appendContextMenu(Menu* menu) override {
        Axon* m = dynamic_cast<Axon*>(module);
        if (!m) return;
        menu->addChild(new MenuSeparator);
        menu->addChild(createIndexPtrSubmenuItem("Anti-aliasing",
            {"Off", "×2 oversampling", "×4 oversampling", "×8 oversampling"}, &m->oversampleMode));
    }
};

Model* modelAxon = createModel<Axon, AxonWidget>("Axon");
