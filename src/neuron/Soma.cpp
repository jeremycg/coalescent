#include "plugin.hpp"
#include "../dsp/rk4.hpp"
#include "../dsp/display_snapshot.hpp"
#include "../dsp/completed_path.hpp"
#include "tanh_approx.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>

// SOMA — a bursting/chaotic neuron oscillator on the Hindmarsh-Rose (HR) system.
//
// The sibling of Axon: HR adds a third, slow state variable `z` (adaptation) to
// the two-variable FitzHugh-Nagumo core, which is exactly what produces bursts
// (trains of spikes separated by quiescence) and, in a window of current, chaos.
// `x` (membrane potential) is the audio output; pitch is the *speed* the
// dimensionless simulation runs at, calibrated so the within-burst spike rate is
// C4 at 0 V (open-loop — the period is emergent, so CURRENT/BURST/ADAPT pull it).
//
//   dx/dt = y - a·x³ + b·x² - z + I      (fast: membrane potential, audio out)
//   dy/dt = c - d·x² - y                 (fast recovery / spiking)
//   dz/dt = r·( s·(x - x_R) - z )         (slow adaptation / bursting)
//
// a,b,c,d,x_R are fixed at the standard values; CURRENT=I, BURST=r, ADAPT=s are
// the exposed bifurcation controls. The integrator (RK4 + pitch-adaptive
// substepping) is the same engine as Axon, with one extra equation for the slow
// adaptation variable z.
//
// Polyphonic: up to 16 independent HR voices, voice count from V/OCT (falling
// back to TRIG). The display traces every active voice on its own hue, stepped
// across a narrow band around Soma's amber accent.

static const int TRAIL = 512;   // length of the (x,z) phase-portrait trail
// Screen-space bins for the latest completed burst/spike path the slow play head rides.
static const int   ORBIT_N    = 192;
static const float ORBIT_SECS = 7.f;   // wall-clock seconds per lap / round trip

struct Soma : Module {

    enum ParamId {
        PITCH_PARAM,
        CURRENT_PARAM,
        BURST_PARAM,
        ADAPT_PARAM,
        CURRENT_ATT_PARAM,
        BURST_ATT_PARAM,
        PARAMS_LEN
    };

    enum InputId {
        VOCT_INPUT,
        CURRENT_INPUT,
        BURST_INPUT,
        TRIG_INPUT,
        SYNC_INPUT,
        INPUTS_LEN
    };

    enum OutputId {
        OUT_OUTPUT,
        SPIKE_OUTPUT,
        Z_OUTPUT,
        OUTPUTS_LEN
    };

    enum LightId {
        LIGHTS_LEN
    };

    // ─── Per-voice persistent state (polyphonic, up to 16 channels) ─────────
    static const int MAX_POLY = 16;
    static const int GROUPS   = MAX_POLY / 4;   // voices processed 4-wide (simd::float_4)
    simd::float_4 xx4[GROUPS], yy4[GROUPS], zz4[GROUPS];   // HR state per voice
    simd::float_4 trigPulse4[GROUPS] = {};
    dsp::TSchmittTrigger<simd::float_4> trigIn4[GROUPS];
    dsp::TSchmittTrigger<simd::float_4> syncIn4[GROUPS];   // SYNC edge detector (hard reset to rest)
    dsp::TSchmittTrigger<simd::float_4> spikeDet4[GROUPS];
    dsp::PulseGenerator spikeGen[MAX_POLY];                // SPIKE pulse shaper (scalar, per voice)
    dsp::TRCFilter<simd::float_4> dcBlock4[GROUPS];
    // Anti-aliasing decimators (one per group; only the active factor is used).
    dsp::Decimator<2, 8, simd::float_4> decim2_4[GROUPS];
    dsp::Decimator<4, 8, simd::float_4> decim4_4[GROUPS];
    dsp::Decimator<8, 8, simd::float_4> decim8_4[GROUPS];
    int   oversampleMode = 2;         // 0=off, 1=×2, 2=×4, 3=×8
    int   channels = 1;               // active voice count
    float lastFs = 0.f;
    int   lastOs = 0;                 // detects oversample change to refresh coefficients
    float trigDecay = 0.f;            // cached per-sample TRIG pulse decay (fn of fs)
    float antiDenorm = 1e-18f;        // alternating sub-LSB dither, anti-denormal on dcBlock

    // ─── Display trail (x,z phase portrait) ─────────────────────────────────
    // Lock-free triple-buffer snapshot (see Axon / dsp/display_snapshot.hpp): the
    // audio thread appends to a per-voice ring and ~45 Hz publishes a complete
    // frame; the head index and active-voice count travel with the arrays so
    // trails and head dots stay coherent.
    float trailX[MAX_POLY][TRAIL] = {}, trailZ[MAX_POLY][TRAIL] = {};
    int   trailIdx = 0, trailDecim = 0;
    // The burst boundary is an upward crossing of a fast-minus-slow z envelope.
    // Both filters run in dimensionless model time, so the same detector follows
    // tonic spikes, bursts and low-pitch chaotic regimes. The completed path is
    // simplified incrementally, preserving short x tails without predicting its period.
    using OrbitPath = coalescent::CompletedPath<ORBIT_N>;
    using OrbitPoint = OrbitPath::Point;
    OrbitPath orbitPath[MAX_POLY];
    float zFast[MAX_POLY] = {}, zSlow[MAX_POLY] = {};
    int orbitCommitClock[MAX_POLY] = {};
    bool zFilterInit[MAX_POLY] = {};
    struct DisplayFrame {
        float x[MAX_POLY][TRAIL] = {}, z[MAX_POLY][TRAIL] = {};
        OrbitPoint orbit[MAX_POLY][ORBIT_N] = {}; // latest completed cycle the dot rides
        int   orbitCount[MAX_POLY] = {};
        bool  orbitValid[MAX_POLY] = {};
        bool  orbitClosed[MAX_POLY] = {};
        int   head = 0;
        int   channels = 1;
    };
    coalescent::DisplaySnapshot<DisplayFrame> displaySnapshot;
    int   dispClock = 0;
    int   dispPeriod = 980;           // samples between snapshots (~45 Hz; refreshed on SR change)

    // ─── Fixed HR parameters ────────────────────────────────────────────────
    static constexpr float A = 1.f, B = 3.f, C = 1.f, D = 5.f, XR = -1.6f;

    // ─── Tunable constants (RATE_CAL from tools/stability/soma.cpp) ──────────
    static constexpr float RATE_CAL    = 55.364003f; // tonic spike period (I=2.0, r=0.03) → C4 at 0 V
    static constexpr float HSUB_MAX    = 0.05f;
    static constexpr int   MIN_SUB     = 2;   // substep floor: 2 holds 0-cent accuracy at default (half the RK4 cost of a 4-floor)
    static constexpr int   MAX_SUB     = 64;
    static constexpr float TRIG_AMP    = 1.0f;        // injected current pulse (kicks a burst from rest)
    static constexpr float TRIG_TAU_MS = 5.f;
    static constexpr float SPIKE_THRESH = 1.0f;
    static constexpr float OUT_GAIN    = 1.5f;
    static constexpr float Z_CENTER    = 2.0f;        // z midpoint, so Z_OUT is bipolar around a burst
    static constexpr float Z_GAIN      = 2.5f;
    static constexpr float CV_DEPTH    = 0.2f;        // ±5V CV × att → ±1.0 (linear I; log₂ for r)
    static constexpr float STATE_MAX   = 25.f;        // backstop; normal |y| peaks ~12, max observed ~21
    static constexpr float R_MIN       = 0.001f, R_MAX = 0.05f;

    Soma() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);

        configParam(PITCH_PARAM, -4.f, 4.f, 0.f, "Pitch", " Hz", 2.f, dsp::FREQ_C4);
        configParam(CURRENT_PARAM, 0.4f, 4.0f, 2.0f, "Current (excitability)");
        // BURST stores log₂(r); display base 2 shows the adaptation rate r directly.
        configParam(BURST_PARAM, std::log2(R_MIN), std::log2(R_MAX), std::log2(0.03f),
                    "Burst / adaptation rate r", "", 2.f);
        configParam(ADAPT_PARAM, 1.0f, 5.0f, 4.0f, "Adaptation strength s");
        configParam(CURRENT_ATT_PARAM, -1.f, 1.f, 0.f, "Current CV");
        configParam(BURST_ATT_PARAM, -1.f, 1.f, 0.f, "Burst CV");

        configInput(VOCT_INPUT, "1V/oct pitch");
        configInput(CURRENT_INPUT, "Current CV");
        configInput(BURST_INPUT, "Burst rate CV");
        configInput(TRIG_INPUT, "Trigger (fire a burst)");
        configInput(SYNC_INPUT, "Sync (hard reset to rest)");

        configOutput(OUT_OUTPUT, "Audio (membrane potential x)");
        configOutput(SPIKE_OUTPUT, "Spike (trigger on each spike)");
        configOutput(Z_OUTPUT, "Adaptation z (burst-envelope CV)");

        for (int g = 0; g < GROUPS; g++) { xx4[g] = -1.6f; yy4[g] = -11.8f; zz4[g] = 2.f; }
        // A sub-pixel spacing stagger keeps identical poly voices from running
        // their bounded guide simplification on the same audio callback.
        for (int c = 0; c < MAX_POLY; ++c) {
            orbitPath[c].setMetric(50.f / 4.2f, 34.f / 3.7f);
            orbitPath[c].setMarkerAges(0.5f, 8192.f);
            orbitPath[c].setMinimumSpacing(0.15f + 0.0015f * c);
            orbitCommitClock[c] = 1 << 24;
        }
    }

    void onReset() override {
        oversampleMode = 2;   // restore default anti-aliasing (×4) on Initialize
        for (int g = 0; g < GROUPS; g++) {
            xx4[g] = -1.6f; yy4[g] = -11.8f; zz4[g] = 2.f; trigPulse4[g] = 0.f;
            trigIn4[g].reset(); syncIn4[g].reset(); spikeDet4[g].reset();
            dcBlock4[g].reset(); decim2_4[g].reset(); decim4_4[g].reset(); decim8_4[g].reset();
        }
        for (int c = 0; c < MAX_POLY; c++) {
            spikeGen[c].reset();
            orbitPath[c].reset();
            zFast[c] = zSlow[c] = 0.f;
            zFilterInit[c] = false;
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
            switch ((int) json_integer_value(j)) {
                case 1:  oversampleMode = 0; break;
                case 2:  oversampleMode = 1; break;
                case 8:  oversampleMode = 3; break;
                default: oversampleMode = 2; break;   // ×4 (and anything unexpected)
            }
        }
        else if (json_t* j = json_object_get(root, "oversample")) {
            // Pre-2.0.3 index meant 0=Off, 1=×4, 2=×8 — remap into Off/×2/×4/×8.
            int v = (int) json_integer_value(j);
            oversampleMode = v <= 0 ? 0 : (v == 1 ? 2 : 3);
        }
    }

    // HR derivatives — the FHN f() with one extra (slow) line. Generic over the
    // scalar type T so one expression serves scalar float and simd::float_4 (four
    // poly voices at once); the shared coalescent::rk4 step is likewise templated.
    template <typename T>
    static inline void f(T x, T y, T z, T Itot, T r, T s, T& dx, T& dy, T& dz) {
        dx = y - A*x*x*x + B*x*x - z + Itot;
        dy = C - D*x*x - y;
        dz = r * (s * (x - XR) - z);
    }

    // One-pole coefficient for a model-time step. This [2/2] Padé form closely
    // approximates 1-exp(-step/tau) over the bounded integrator step while
    // keeping transcendental functions out of the per-voice audio path.
    static inline float modelOnePole(float step, float tau) {
        const float x = step / tau;
        return clamp(x / (1.f + 0.5f * x + x * x / 12.f), 0.f, 1.f);
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

        // Voice count follows V/OCT, falling back to TRIG so a poly trigger
        // (percussive bursts, no pitch cable) still spreads across voices.
        channels = std::max({inputs[VOCT_INPUT].getChannels(),
                             inputs[TRIG_INPUT].getChannels(), 1});
        outputs[OUT_OUTPUT].setChannels(channels);
        outputs[SPIKE_OUTPUT].setChannels(channels);
        outputs[Z_OUTPUT].setChannels(channels);

        // Shared knob values; per-voice CV is added inside the loop.
        const float Ibase     = params[CURRENT_PARAM].getValue();
        const float Iatt      = params[CURRENT_ATT_PARAM].getValue();
        const float rLogBase  = params[BURST_PARAM].getValue();
        const float rAtt      = params[BURST_ATT_PARAM].getValue();
        const float s         = params[ADAPT_PARAM].getValue();
        const float pitchKnob = params[PITCH_PARAM].getValue();

        // Alternating sub-LSB dither (anti-denormal), toggled once per sample.
        antiDenorm = -antiDenorm;

        const simd::float_4 sVec(s);    // ADAPT has no CV → same for every lane
        const int nGroups = (channels + 3) / 4;

        for (int g = 0; g < nGroups; g++) {
            const int base = g * 4;

            // Lanes at/above the channel count are inactive: they still compute
            // (SIMD lanes can't be skipped) but must not evolve persistent state,
            // preserving the documented "silent voices keep their state" semantic.
            const simd::float_4 activeMask =
                simd::float_4(0.f, 1.f, 2.f, 3.f) < (float) (channels - base);

            // ── params → physics (4 lanes) ──
            simd::float_4 I = Ibase
                + inputs[CURRENT_INPUT].getPolyVoltageSimd<simd::float_4>(base) * Iatt * CV_DEPTH;
            // BURST is log₂(r); CV adds in the log domain (a multiplicative nudge on r).
            simd::float_4 r = simd::clamp(dsp::approxExp2_taylor5(rLogBase
                + inputs[BURST_INPUT].getPolyVoltageSimd<simd::float_4>(base) * rAtt * CV_DEPTH), R_MIN, R_MAX);

            simd::float_4 x = xx4[g], y = yy4[g], z = zz4[g];

            // ── hard sync: rising edge re-seeds those lanes to rest (all 3 vars) ──
            simd::float_4 syncMask = syncIn4[g].process(
                inputs[SYNC_INPUT].getPolyVoltageSimd<simd::float_4>(base), 0.1f, 1.f);
            const int syncBits = simd::movemask(syncMask);
            x = simd::ifelse(syncMask, -1.6f, x);
            y = simd::ifelse(syncMask, -11.8f, y);
            z = simd::ifelse(syncMask, 2.f, z);

            // ── excitable trigger: rising edge injects a decaying current pulse ──
            simd::float_4 trigMask = trigIn4[g].process(
                inputs[TRIG_INPUT].getPolyVoltageSimd<simd::float_4>(base), 0.1f, 1.f);
            simd::float_4 tp = simd::ifelse(trigMask, TRIG_AMP, trigPulse4[g]);
            simd::float_4 Itot = I + tp;
            tp *= trigDecay;
            tp = simd::ifelse(simd::abs(tp) < 1e-30f, 0.f, tp);
            trigPulse4[g] = simd::ifelse(activeMask, tp, trigPulse4[g]);

            // ── pitch = simulation speed (open-loop). Run the group at its max K
            // (a slower lane just gets a smaller-than-needed h; accuracy preserved). ──
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

            // ── render os oversampled samples (full output chain each), decimate ──
            simd::float_4 st[3] = {x, y, z};
            simd::float_4 osBuf[8];
            // Hoist the z-envelope one-pole coefficients: subTau[l] is fixed across the
            // oversample loop, so computing them once per lane here (instead of per lane
            // per oversample step) saves the divide/libm work at 16 voices ×8.
            float zFastCoef[4] = {}, zSlowCoef[4] = {};
            for (int l = 0; l < 4 && base + l < channels; ++l) {
                const int c = base + l;
                orbitCommitClock[c] = std::min(orbitCommitClock[c] + 1, dispPeriod);
                zFastCoef[l] = modelOnePole(subTau[l], 15.f);
                zSlowCoef[l] = modelOnePole(subTau[l], 300.f);
            }
            for (int o = 0; o < os; o++) {
                for (int k = 0; k < K; k++)
                    coalescent::rk4<3>(st, h, [&](const simd::float_4* Y, simd::float_4* d) {
                        f(Y[0], Y[1], Y[2], Itot, r, sVec, d[0], d[1], d[2]);
                    });

                // Backstop: reset any non-finite / runaway lane to rest, then clamp.
                simd::float_4 finite = (st[0] == st[0]) & (st[1] == st[1]) & (st[2] == st[2])
                    & (simd::abs(st[0]) < 1e6f) & (simd::abs(st[1]) < 1e6f) & (simd::abs(st[2]) < 1e6f);
                const int orbitResetBits = (o == 0 ? syncBits : 0)
                    | ((~simd::movemask(finite)) & 0xf);
                st[0] = simd::ifelse(finite, st[0], -1.6f);
                st[1] = simd::ifelse(finite, st[1], -11.8f);
                st[2] = simd::ifelse(finite, st[2], 2.f);
                st[0] = simd::clamp(st[0], -STATE_MAX, STATE_MAX);
                st[1] = simd::clamp(st[1], -STATE_MAX, STATE_MAX);
                st[2] = simd::clamp(st[2], -STATE_MAX, STATE_MAX);

                simd::float_4 spikeMask = spikeDet4[g].process(st[0], 0.f, SPIKE_THRESH);
                int sm = simd::movemask(spikeMask);
                for (int l = 0; l < 4 && base + l < channels; l++)
                    if (sm & (1 << l)) spikeGen[base + l].trigger(1e-3f);

                dcBlock4[g].process(st[0] + antiDenorm);
                osBuf[o] = 5.f * coalescent::fastTanh(dcBlock4[g].highpass() * OUT_GAIN);

                // The guide needs the oversampled states: Soma's positive-x peak
                // can fall entirely between host-rate samples at high pitch.
                for (int l = 0; l < 4 && base + l < channels; ++l) {
                    const int c = base + l;
                    const float xNow = st[0][l];
                    const float zNow = st[2][l];
                    const float modelStep = subTau[l];
                    const bool discontinuity = orbitResetBits & (1 << l);
                    if (!zFilterInit[c] || discontinuity) {
                        zFast[c] = zSlow[c] = zNow;
                        zFilterInit[c] = true;
                    } else {
                        zFast[c] += (zNow - zFast[c]) * zFastCoef[l];
                        zSlow[c] += (zNow - zSlow[c]) * zSlowCoef[l];
                    }

                    const bool committed = orbitPath[c].push(
                        xNow, zNow, zFast[c] - zSlow[c], modelStep,
                        orbitCommitClock[c] >= dispPeriod, discontinuity);
                    if (committed)
                        orbitCommitClock[c] = 0;
                }
            }
            // Write back active lanes only; inactive lanes keep their frozen state.
            xx4[g] = simd::ifelse(activeMask, st[0], xx4[g]);
            yy4[g] = simd::ifelse(activeMask, st[1], yy4[g]);
            zz4[g] = simd::ifelse(activeMask, st[2], zz4[g]);

            simd::float_4 audio = os == 1 ? osBuf[0]
                : (os == 2 ? decim2_4[g].process(osBuf)
                :  os == 4 ? decim4_4[g].process(osBuf) : decim8_4[g].process(osBuf));

            // ── outputs (4 channels at once) ──
            outputs[OUT_OUTPUT].setVoltageSimd(audio, base);
            outputs[Z_OUTPUT].setVoltageSimd(simd::clamp((st[2] - Z_CENTER) * Z_GAIN, -5.f, 5.f), base);
            simd::float_4 spikeOut;
            for (int l = 0; l < 4; l++)
                spikeOut[l] = spikeGen[base + l].process(args.sampleTime) ? 10.f : 0.f;
            outputs[SPIKE_OUTPUT].setVoltageSimd(spikeOut, base);
        }

        // ── feed the (x,z) phase-portrait trails (all voices; decimated, bursts are slow) ──
        if (++trailDecim >= 6) {
            trailDecim = 0;
            for (int c = 0; c < channels; c++) {
                trailX[c][trailIdx] = xx4[c / 4][c % 4];
                trailZ[c][trailIdx] = zz4[c / 4][c % 4];
            }
            trailIdx = (trailIdx + 1) % TRAIL;
        }
        if (++dispClock >= dispPeriod) {                 // publish display snapshot ~45 Hz
            dispClock = 0;
            DisplayFrame& fr = displaySnapshot.writable();
            for (int c = 0; c < channels; c++) {
                std::copy(trailX[c], trailX[c] + TRAIL, fr.x[c]);
                std::copy(trailZ[c], trailZ[c] + TRAIL, fr.z[c]);
                std::copy(orbitPath[c].path(), orbitPath[c].path() + ORBIT_N, fr.orbit[c]);
                fr.orbitCount[c] = orbitPath[c].pathSize();
                fr.orbitValid[c] = orbitPath[c].hasPath();
                fr.orbitClosed[c] = orbitPath[c].isClosed();
            }
            fr.head = trailIdx;
            fr.channels = channels;
            displaySnapshot.publish();
        }
    }
};


// ─── Phase-portrait display ────────────────────────────────────────────────
// The Hindmarsh-Rose attractor projected onto the (x, z) plane: fast spikes
// sweep horizontally while the slow adaptation z drifts up and down, so a burst
// shows as a cluster of spikes climbing along z and the quiescent gap as the
// slow return. In the chaotic window the trail never quite repeats. The faint
// diagonal is the z-nullcline z = s·(x − x_R) (where the slow drift reverses).
// With several voices patched, each attractor is drawn on its own hue stepped
// across a narrow amber band.
struct SomaDisplay : Widget {
    Soma* module = nullptr;
    // Hold one immutable guide for a complete visual lap. Chaotic cycles differ
    // legitimately; swapping guides mid-lap would teleport the slow play head.
    Soma::OrbitPoint playOrbit[Soma::MAX_POLY][ORBIT_N] = {};
    int playCount[Soma::MAX_POLY] = {};
    long long playLap[Soma::MAX_POLY] = {};
    bool playValid[Soma::MAX_POLY] = {};
    bool playClosed[Soma::MAX_POLY] = {};

    static constexpr float XMIN = -2.0f, XMAX = 2.2f;
    static constexpr float ZMIN = 0.3f, ZMAX = 4.0f;

    // Voice → hue (0..1) within a band centred on Soma's amber accent (~0.082).
    static float voiceHue(int v, int nv) {
        const float center = 0.082f, halfBand = 0.05f;
        if (nv <= 1) return center;
        return center - halfBand + 2.f * halfBand * v / (nv - 1);
    }

    void draw(const DrawArgs& args) override {
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
            nvgScissor(args.vg, 0.f, 0.f, box.size.x, box.size.y);

            auto X = [&](float x) { return (x - XMIN) / (XMAX - XMIN) * box.size.x; };
            auto Y = [&](float z) { return box.size.y - (z - ZMIN) / (ZMAX - ZMIN) * box.size.y; };

            float s = module ? module->params[Soma::ADAPT_PARAM].getValue() : 4.f;
            static const Soma::DisplayFrame dummy;
            const Soma::DisplayFrame& fr = module ? module->displaySnapshot.consume() : dummy;
            int nv = module ? fr.channels : 1;
            const double playCycles = system::getTime() / (double) ORBIT_SECS;
            const long long currentLap = (long long) std::floor(playCycles);
            const float lapFraction = (float) (playCycles - std::floor(playCycles));
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

            // z-nullcline: z = s·(x − x_R).
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, X(XMIN), Y(s * (XMIN - Soma::XR)));
            nvgLineTo(args.vg, X(XMAX), Y(s * (XMAX - Soma::XR)));
            nvgStrokeColor(args.vg, nvgRGBA(0x90, 0x50, 0x70, 0x55));
            nvgStrokeWidth(args.vg, 1.f);
            nvgStroke(args.vg);

            // Trails batched into TRAIL_BANDS constant-alpha bands (one polyline stroke
            // each) rather than a stroke per TRAIL segment — ~16× fewer NanoVG calls at
            // 16 voices for a stepped fade indistinguishable at panel scale. (See Axon for the rationale.)
            if (module) {
                const int TRAIL_BANDS = 32;
                int idx = fr.head;   // coherent with arrays
                float trailA = clamp(204.f - (nv - 1) * 6.f, 112.f, 204.f);  // 0x70..0xcc
                nvgLineCap(args.vg, NVG_ROUND);
                nvgLineJoin(args.vg, NVG_ROUND);
                nvgStrokeWidth(args.vg, 1.6f);
                for (int v = 0; v < nv; v++) {
                    const float* dx = fr.x[v];
                    const float* dz = fr.z[v];
                    float hue = voiceHue(v, nv);
                    const Soma::OrbitPoint* orbit = playOrbit[v];

                    // Show the completed guide beneath the fading ring. This keeps
                    // slow bursts and narrow spike tails visibly connected to the
                    // play head even when the decimated 512-point trail omits them.
                    if (playValid[v]) {
                        const int orbitCount = playCount[v];
                        nvgBeginPath(args.vg);
                        nvgMoveTo(args.vg, X(orbit[0].x), Y(orbit[0].y));
                        for (int b = 1; b < orbitCount; ++b)
                            nvgLineTo(args.vg, X(orbit[b].x), Y(orbit[b].y));
                        if (playClosed[v])
                            nvgLineTo(args.vg, X(orbit[0].x), Y(orbit[0].y));
                        const int guideAlpha = clamp(64 - (nv - 1) * 2, 32, 64);
                        nvgStrokeColor(args.vg, nvgHSLA(hue, 0.95f, 0.60f, guideAlpha));
                        nvgStrokeWidth(args.vg, 1.f);
                        nvgStroke(args.vg);
                        nvgStrokeWidth(args.vg, 1.6f);
                    }

                    int curBand = -1;
                    auto strokeBand = [&](int band) {
                        float alpha = (band + 0.5f) / TRAIL_BANDS;   // older = dimmer
                        nvgStrokeColor(args.vg, nvgHSLA(hue, 0.95f, 0.60f, (int)(alpha * trailA)));
                        nvgStroke(args.vg);
                    };
                    for (int k = 1; k < TRAIL; k++) {
                        int band = (k * TRAIL_BANDS) / TRAIL;
                        int i0 = (idx + k - 1) % TRAIL;
                        int i1 = (idx + k) % TRAIL;
                        if (band != curBand) {
                            if (curBand >= 0) strokeBand(curBand);
                            nvgBeginPath(args.vg);
                            nvgMoveTo(args.vg, X(dx[i0]), Y(dz[i0]));
                            curBand = band;
                        }
                        nvgLineTo(args.vg, X(dx[i1]), Y(dz[i1]));
                    }
                    if (curBand >= 0) strokeBand(curBand);
                    // Slow play head over the latest completed real burst/spike path.
                    // Uniform screen-space bins make the narrow horizontal tails visible.
                    // A chaotic path whose marker endpoints do not meet is traversed in
                    // both directions rather than closed with an invented off-trace chord.
                    float hx, hy;
                    if (playValid[v]) {
                        const float lap = lapFraction;
                        const int orbitCount = playCount[v];
                        const Soma::OrbitPoint head = coalescent::sampleCompletedPath(
                            orbit, orbitCount, playClosed[v], lap, X, Y);
                        hx = X(head.x); hy = Y(head.y);
                    } else {
                        const int newest = (idx + TRAIL - 1) % TRAIL;
                        hx = X(dx[newest]); hy = Y(dz[newest]);
                    }
                    nvgBeginPath(args.vg); nvgCircle(args.vg, hx, hy, 4.f);
                    nvgFillColor(args.vg, nvgHSLA(hue, 0.95f, 0.72f, 0x55)); nvgFill(args.vg);
                    nvgBeginPath(args.vg); nvgCircle(args.vg, hx, hy, 2.f);
                    nvgFillColor(args.vg, nvgHSLA(hue, 0.65f, 0.88f, 0xff)); nvgFill(args.vg);
                }
            } else {
                // Browser preview: an illustrative bursting loop.
                nvgBeginPath(args.vg);
                for (int i = 0; i <= 80; i++) {
                    float t = 2.f * (float)M_PI * i / 80.f;
                    float px = X(0.3f + 1.2f * std::cos(t) + 0.3f * std::sin(6*t));
                    float py = Y(2.0f + 1.1f * std::sin(t));
                    if (i == 0) nvgMoveTo(args.vg, px, py); else nvgLineTo(args.vg, px, py);
                }
                nvgStrokeColor(args.vg, nvgRGBA(0xff, 0x9b, 0x3a, 0x99));
                nvgStrokeWidth(args.vg, 1.6f);
                nvgStroke(args.vg);
            }

            std::shared_ptr<Font> font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
            if (font) {
                nvgFontFaceId(args.vg, font->handle);
                nvgFontSize(args.vg, mm2px(3.2f));
                nvgTextLetterSpacing(args.vg, mm2px(0.5f));
                nvgFillColor(args.vg, nvgRGBA(0xff, 0xc0, 0x9a, 0xcc));
                nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
                nvgText(args.vg, mm2px(2.6f), mm2px(2.2f), "SOMA", NULL);
                nvgTextLetterSpacing(args.vg, 0.f);

                // Voice-count badge when polyphonic.
                if (nv > 1) {
                    char buf[16];
                    snprintf(buf, sizeof(buf), "%dv", nv);
                    nvgFontSize(args.vg, mm2px(2.2f));
                    nvgFillColor(args.vg, nvgRGBA(0xb0, 0x80, 0x60, 0xcc));
                    nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_TOP);
                    nvgText(args.vg, box.size.x - mm2px(2.4f), mm2px(2.2f), buf, NULL);
                }

                nvgFontSize(args.vg, mm2px(2.2f));
                nvgFillColor(args.vg, nvgRGBA(0xb0, 0x80, 0x60, 0xaa));
                nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_BOTTOM);
                nvgText(args.vg, box.size.x - mm2px(2.4f), box.size.y - mm2px(1.8f), "x–z", NULL);
            }

            nvgResetScissor(args.vg);
        }
        Widget::drawLayer(args, layer);
    }
};


// ─── Widget ─────────────────────────────────────────────────────────────────

struct SomaLabels : widget::Widget {

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

        lbl( 7.5f,   46.f, 2.1f, dim, "PITCH");
        lbl(22.82f, 46.f, 2.1f, dim, "CURRENT");
        lbl(38.14f, 46.f, 2.1f, dim, "BURST");
        lbl(53.46f, 46.f, 2.1f, dim, "ADAPT");

        lbl(22.82f, 90.f, 2.1f, dim, "I CV");
        lbl(38.14f, 90.f, 2.1f, dim, "r CV");

        lbl( 6.5f, 118.5f, 2.1f, dim,    "V/OCT");
        lbl(15.5f, 118.5f, 2.1f, dim,    "TRIG");
        lbl(24.5f, 118.5f, 2.1f, dim,    "SYNC");
        lbl(36.5f, 118.5f, 2.1f, outclr, "OUT");
        lbl(45.5f, 118.5f, 2.1f, outclr, "SPIKE");
        lbl(54.5f, 118.5f, 2.1f, outclr, "Z");

        nvgRestore(args.vg);
    }
};

struct SomaWidget : ModuleWidget {

    SomaWidget(Soma* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Soma.svg")));
        addPanelLabels<SomaLabels>(this);

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   122.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 122.0f))));

        SomaDisplay* disp = new SomaDisplay();
        disp->module = module;
        disp->box.pos  = mm2px(Vec(5.5f, 8.f));
        disp->box.size = mm2px(Vec(50.f, 34.f));
        addChild(disp);

        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec( 7.5f,   54.f)), module, Soma::PITCH_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(22.82f, 54.f)), module, Soma::CURRENT_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(38.14f, 54.f)), module, Soma::BURST_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(53.46f, 54.f)), module, Soma::ADAPT_PARAM));

        addParam(createParamCentered<Trimpot>(    mm2px(Vec(22.82f, 72.f)), module, Soma::CURRENT_ATT_PARAM));
        addParam(createParamCentered<Trimpot>(    mm2px(Vec(38.14f, 72.f)), module, Soma::BURST_ATT_PARAM));
        addInput(createInputCentered<PJ301MPort>( mm2px(Vec(22.82f, 84.f)), module, Soma::CURRENT_INPUT));
        addInput(createInputCentered<PJ301MPort>( mm2px(Vec(38.14f, 84.f)), module, Soma::BURST_INPUT));

        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec( 6.5f, 112.f)), module, Soma::VOCT_INPUT));
        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec(15.5f, 112.f)), module, Soma::TRIG_INPUT));
        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec(24.5f, 112.f)), module, Soma::SYNC_INPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(36.5f, 112.f)), module, Soma::OUT_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(45.5f, 112.f)), module, Soma::SPIKE_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(54.5f, 112.f)), module, Soma::Z_OUTPUT));
    }

    void appendContextMenu(Menu* menu) override {
        Soma* m = dynamic_cast<Soma*>(module);
        if (!m) return;
        menu->addChild(new MenuSeparator);
        menu->addChild(createIndexPtrSubmenuItem("Anti-aliasing",
            {"Off", "×2 oversampling", "×4 oversampling", "×8 oversampling"}, &m->oversampleMode));
    }
};

Model* modelSoma = createModel<Soma, SomaWidget>("Soma");
