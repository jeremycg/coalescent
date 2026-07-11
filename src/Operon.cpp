#include "plugin.hpp"
#include "dsp/rk4.hpp"
#include "dsp/display_snapshot.hpp"
#include "tanh_approx.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>

// OPERON — a three-phase oscillator on the repressilator (Elowitz & Leibler,
// 2000): three genes in a ring, each repressing the next, whose delayed negative
// feedback sustains oscillation. Six ODEs (mRNA + protein per gene) integrated
// with the shared RK4 stepper; the three protein levels run ~120° apart, so the
// module is natively three-phase — stereo, chords, phase-locked modulation, and
// at sub-audio speed a three-phase LFO / clock.
//
//   dm_i/dt = -m_i + alpha / (1 + rep_i^n) + alpha0     (rep_i = protein of prev gene)
//   dp_i/dt = -beta · (p_i - m_i)
//
// ALPHA is the drive (crosses a Hopf threshold into oscillation), HILL the
// cooperativity (sine → pulse), BETA the decay ratio (character/phase skew),
// LEAK the basal floor. Pitch is the simulation *speed* (emergent period).
//
// Parametrisation follows the dimensionless Elowitz–Leibler form (BioModels
// BIOMD0000000012); defaults sit in the oscillating regime with the Hopf
// threshold reachable low on the ALPHA dial. The classic Elowitz set
// (alpha≈216, beta≈0.2, n=2, alpha0≈0.216) is a larger, slower voicing.

struct Operon : Module {

    enum ParamId {
        PITCH_PARAM,
        ALPHA_PARAM,
        HILL_PARAM,
        BETA_PARAM,
        LEAK_PARAM,
        ALPHA_ATT_PARAM,
        HILL_ATT_PARAM,
        BETA_ATT_PARAM,
        PARAMS_LEN
    };

    enum InputId {
        VOCT_INPUT,
        ALPHA_INPUT,
        HILL_INPUT,
        BETA_INPUT,
        PERTURB_INPUT,
        INPUTS_LEN
    };

    enum OutputId {
        OUT1_OUTPUT,
        OUT2_OUTPUT,
        OUT3_OUTPUT,
        GATE1_OUTPUT,
        GATE2_OUTPUT,
        GATE3_OUTPUT,
        OUTPUTS_LEN
    };

    enum LightId {
        LIGHTS_LEN
    };

    // ─── State (transient; not saved) ───────────────────────────────────────
    // Seeded ASYMMETRICALLY — with equal ICs the system stays on the symmetric
    // subspace and never oscillates.
    float m[3] = {0.2f, 0.1f, 0.3f};
    float p[3] = {0.1f, 0.4f, 0.15f};
    float lastCentered[3] = {};          // upward-crossing detection on p[k] - pStar
    float pStar = 0.f;                   // cached symmetric fixed point (= mStar)
    float pStarA = -1.f, pStarN = -1.f, pStarL = -1.f;  // params pStar was solved for
    bool  pStarValid = false;
    dsp::PulseGenerator gateGen[3];

    // Hill-response LUT: 1/(1+repⁿ) precomputed over rep, keyed on n — replaces the
    // per-sample pow with a lookup (~1.5× on the RK4 derivative), sub-cent even at
    // n=8 (8192 pts). Rebuilds are rate-limited (see process()) and fall back to
    // direct pow while n moves, so dialling HILL or CV'ing it can't thrash the
    // 8192-pow rebuild.
    static constexpr int   HILL_LUT_N     = 8192;
    static constexpr float HILL_LUT_XMAX  = 64.f;
    static constexpr int   HILL_LUT_SLICE = 256;   // entries filled per sample during a rebuild
    float hillLut[HILL_LUT_N + 1];
    float lutN = -1.f;      // n the active LUT was built for
    float nSettleRef = -1e9f; // anchor for the settle detector (see process())
    int   rebuildClock = 0; // consecutive samples n has stayed within N_MOVE_EPS of nSettleRef
    bool  lutValid = false;
    int   buildPos = -1;    // incremental rebuild cursor: -1 idle, else next entry to fill
    float buildN = -1.f;    // n the in-progress rebuild is filling for
    static constexpr float N_MOVE_EPS = 1e-4f;   // per-sample |Δn| above this = "actively modulated"

    static float hillEntry(int i, float n) {
        return 1.f / (1.f + std::pow(HILL_LUT_XMAX * i / HILL_LUT_N, n));
    }
    inline float hillLookup(float rep) const {   // rep already >= 0
        if (rep >= HILL_LUT_XMAX) return hillLut[HILL_LUT_N];
        float f = rep * (HILL_LUT_N / HILL_LUT_XMAX);
        int   i = (int) f;
        return hillLut[i] + (f - i) * (hillLut[i + 1] - hillLut[i]);
    }

    // ─── Display scope: the 3 protein levels over time (lock-free triple buffer) ──
    static const int HIST_N = 256;
    struct OperonScope {
        float p[3][HIST_N] = {};         // centered protein history (ring)
        int   head = 0;                  // next write index
        float peak = 1e-3f;              // smoothed y-scale
    };
    OperonScope liveScope;               // audio-thread only
    coalescent::DisplaySnapshot<OperonScope> displaySnapshot;
    dsp::ClockDivider dispDiv;
    float lastDisplayFs = 0.f;

    // ─── Tunable constants ──────────────────────────────────────────────────
    static constexpr float RATE_CAL     = 12.46f;  // measured default period (tools/stability/operon.cpp) → C4 at 0 V
    static constexpr float HSUB_MAX     = 0.05f;
    static constexpr int   MIN_SUB      = 2;   // substep floor: K=2 holds 0-cent and stays bounded at default
    static constexpr int   MAX_SUB      = 64;
    static constexpr float PITCH_TOTAL_MIN = -8.f, PITCH_TOTAL_MAX = 8.f;
    static constexpr float OUT_GAIN     = 0.3f;   // set offline: default RMS ~3.2V, peaks lightly saturate, high drive grits
    static constexpr float STATE_MAX    = 1e3f;
    static constexpr float PERTURB_GAIN = 0.5f;
    static constexpr float BIAS_EPS     = 1e-6f;   // zero-sum self-start symmetry breaker
    static constexpr float GATE_LEVEL   = 10.f, GATE_TIME = 1e-3f;
    static constexpr float GATE_LOW     = -0.05f, GATE_HIGH = 0.05f;  // centered units
    static constexpr int   CENTER_ITERS = 16;   // bisection fallback iterations
    static constexpr int   NEWTON_ITERS = 4;     // warm-started Newton (1 pow/iter)
    // Separate CV depths — the parameter scales differ a lot.
    static constexpr float ALPHA_CV_DEPTH = 6.f;
    static constexpr float HILL_CV_DEPTH  = 0.6f;
    static constexpr float BETA_CV_DEPTH  = 0.5f;

    Operon() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);

        configParam(PITCH_PARAM, -8.f, 4.f, 0.f, "Frequency", " Hz", 2.f, dsp::FREQ_C4);  // down to ~1 Hz for the LFO/clock use
        configParam(ALPHA_PARAM, 0.f, 60.f, 12.f, "Promoter strength (drive)");
        configParam(HILL_PARAM, 1.2f, 8.f, 2.5f, "Hill coefficient (cooperativity n)");
        configParam(BETA_PARAM, 0.2f, 5.f, 1.0f, "Protein/mRNA decay ratio");
        configParam(LEAK_PARAM, 0.f, 1.0f, 0.05f, "Basal leak");
        configParam(ALPHA_ATT_PARAM, -1.f, 1.f, 0.f, "Drive CV");
        configParam(HILL_ATT_PARAM, -1.f, 1.f, 0.f, "Hill CV");
        configParam(BETA_ATT_PARAM, -1.f, 1.f, 0.f, "Decay-ratio CV");

        configInput(VOCT_INPUT, "1V/oct pitch");
        configInput(ALPHA_INPUT, "Drive CV");
        configInput(HILL_INPUT, "Hill CV");
        configInput(BETA_INPUT, "Decay-ratio CV");
        configInput(PERTURB_INPUT, "Perturbation (into gene 1)");

        configOutput(OUT1_OUTPUT, "Protein 1");
        configOutput(OUT2_OUTPUT, "Protein 2");
        configOutput(OUT3_OUTPUT, "Protein 3");
        configOutput(GATE1_OUTPUT, "Gate 1");
        configOutput(GATE2_OUTPUT, "Gate 2");
        configOutput(GATE3_OUTPUT, "Gate 3");
    }

    void reseedAsymmetric() {
        m[0] = 0.2f; m[1] = 0.1f; m[2] = 0.3f;
        p[0] = 0.1f; p[1] = 0.4f; p[2] = 0.15f;
    }

    void resetGateMemory() {
        for (int k = 0; k < 3; ++k) lastCentered[k] = p[k] - pStar;
    }

    void onReset() override {
        reseedAsymmetric();
        pStarValid = false;
        pStar = 0.f; pStarA = pStarN = pStarL = -1.f;
        lutValid = false; lutN = -1.f; nSettleRef = -1e9f; rebuildClock = 0; buildPos = -1;
        for (int k = 0; k < 3; ++k) { lastCentered[k] = 0.f; gateGen[k].reset(); }
    }

    // Symmetric fixed point pStar (= mStar): the monotonic g(x) = x - alpha/(1+x^n)
    // - alpha0 has one root in (0, hi). g is strictly increasing (g' > 1), so Newton
    // warm-started from the previous root converges in 1–2 steps under per-sample
    // modulation — one pow() per step vs bisection's 16. Fall back to bisection on a
    // cold start (guess<0) or if a step leaves the bracket. max(x,eps) before pow —
    // n is non-integer, and x^(n-1) is reused from x^n as xn/x.
    static float solvePStar(float alpha, float n, float alpha0, float guess) {
        const float hi = std::max(1.f, alpha + alpha0 + 1.f);
        if (guess > 0.f && guess < hi) {
            float x = guess;
            for (int i = 0; i < NEWTON_ITERS; ++i) {
                float xn = std::pow(std::max(x, 1e-6f), n);
                float d  = 1.f + xn;
                float g  = x - alpha / d - alpha0;
                float gp = 1.f + alpha * n * (xn / std::max(x, 1e-6f)) / (d * d);   // g'(x) > 1
                float step = g / gp;
                x -= step;
                if (!(x > 0.f && x < hi)) break;                   // left the bracket → bisection
                if (std::fabs(step) < 1e-5f) return x;             // converged
            }
            // Did not converge within NEWTON_ITERS (a big jump) → fall through to bisection.
        }
        // Cold start / fallback: robust bisection.
        float lo = 0.f, h = hi;
        for (int i = 0; i < CENTER_ITERS; ++i) {
            float mid = 0.5f * (lo + h);
            float g = mid - alpha / (1.f + std::pow(std::max(mid, 0.f), n)) - alpha0;
            if (g > 0.f) h = mid; else lo = mid;
        }
        return 0.5f * (lo + h);
    }

    void process(const ProcessArgs& args) override {
        const float fs = args.sampleRate;

        // ── read + clamp params (per-sample; CV added with separate depths) ──
        float alpha  = clamp(params[ALPHA_PARAM].getValue()
                           + inputs[ALPHA_INPUT].getVoltage() * params[ALPHA_ATT_PARAM].getValue() * ALPHA_CV_DEPTH,
                             0.f, 80.f);
        float n      = clamp(params[HILL_PARAM].getValue()
                           + inputs[HILL_INPUT].getVoltage() * params[HILL_ATT_PARAM].getValue() * HILL_CV_DEPTH,
                             1.01f, 10.f);
        float beta   = clamp(params[BETA_PARAM].getValue()
                           + inputs[BETA_INPUT].getVoltage() * params[BETA_ATT_PARAM].getValue() * BETA_CV_DEPTH,
                             0.05f, 8.f);
        float alpha0 = clamp(params[LEAK_PARAM].getValue(), 0.f, 2.f);
        // Sanitize PERTURB before it enters the ODE: a non-finite CV would propagate
        // through rk4 into rep and hit hillLookup's array index (UB) *before* the
        // post-step finiteness guard below can catch it. Flush NaN/inf to 0, clamp range.
        float perturb = 0.f;
        if (inputs[PERTURB_INPUT].isConnected()) {
            float v = inputs[PERTURB_INPUT].getVoltage();
            if (std::isfinite(v)) perturb = clamp(v, -10.f, 10.f) * PERTURB_GAIN;
        }

        // ── symmetric fixed point for centering (cached on alpha/n/alpha0) ──
        if (!pStarValid || alpha != pStarA || n != pStarN || alpha0 != pStarL) {
            pStar = solvePStar(alpha, n, alpha0, pStarValid ? pStar : -1.f);   // warm-start Newton from the cached root
            pStarA = alpha; pStarN = n; pStarL = alpha0; pStarValid = true;
        }

        // Hill LUT vs direct pow: a STATIC n (knob at rest or a held HILL CV) uses the
        // 8192-point LUT; a moving n (audio-rate, or any drift) uses direct pow, since a
        // rebuilt-at-most-every-N-samples LUT would lag and detune.
        //
        // Settle detector: rebuildClock counts *consecutive* samples n has stayed within
        // N_MOVE_EPS of nSettleRef. ANY larger move — a fast sweep OR slow cumulative
        // drift past the window — resets the counter, the anchor, and any partial build.
        // A rebuild starts only after n has genuinely been still for ~2048 samples, so a
        // slow HILL LFO (whose per-sample delta is tiny but which never settles) stays on
        // direct pow instead of aborting-and-restarting a 256-entry fill every sample.
        if (std::fabs(n - nSettleRef) > N_MOVE_EPS) {
            nSettleRef = n; rebuildClock = 0; buildPos = -1;   // moved → reset window + drop partial build
        } else if (rebuildClock < 2048) {
            ++rebuildClock;                                    // saturates; only the >=2048 threshold matters
        }

        // Rebuild the LUT incrementally rather than in one sample: the full 8192-pow
        // table build measured ~60 µs (a real spike at high SR / many instances). Fill it
        // HILL_LUT_SLICE entries per sample and publish only when complete — the whole
        // rebuild spreads over ~33 samples (<1 ms). Safe to fill in place: the table is
        // only read when hillDirect is false, which requires a *matching* live LUT
        // (lutValid && n≈lutN) — impossible mid-build, since starting a build clears lutValid.
        //
        // Build for nSettleRef (the settle *anchor*), not the instantaneous n, and gate the
        // start on |nSettleRef - lutN|. Anchoring everything to the same band means an
        // in-band micro-oscillation (n wobbling within N_MOVE_EPS of the anchor) can't span
        // the 1e-4 rebuild threshold relative to a build-time n and thrash — it builds once
        // for the anchor and then the LUT covers the whole wobble. rebuildClock is also
        // cleared on completion so a fresh settle window is required before any next build.
        if (buildPos < 0 && rebuildClock >= 2048 && std::fabs(nSettleRef - lutN) > 1e-4f) {
            buildPos = 0; buildN = nSettleRef; lutValid = false;   // start; old LUT is now stale/unreadable
        }
        if (buildPos >= 0) {
            int end = std::min(buildPos + HILL_LUT_SLICE, HILL_LUT_N + 1);
            for (; buildPos < end; ++buildPos) hillLut[buildPos] = hillEntry(buildPos, buildN);
            if (buildPos > HILL_LUT_N) { lutN = buildN; lutValid = true; buildPos = -1; rebuildClock = 0; }  // go live
        }
        const bool hillDirect = !lutValid || std::fabs(n - lutN) > 1e-4f;

        // ── pitch = simulation speed; adaptive substepping ──
        float pitchTotal = clamp(params[PITCH_PARAM].getValue() + inputs[VOCT_INPUT].getVoltage(),
                                 PITCH_TOTAL_MIN, PITCH_TOTAL_MAX);
        float pitchHz = dsp::FREQ_C4 * dsp::approxExp2_taylor5(pitchTotal);
        float dtau = RATE_CAL * pitchHz / fs;
        dtau = std::min(dtau, HSUB_MAX * MAX_SUB);          // hard guard: keep h <= HSUB_MAX
        int   K = clamp((int) std::ceil(dtau / HSUB_MAX), MIN_SUB, MAX_SUB);
        float h = dtau / K;

        // ── derivative over y = [m0,m1,m2,p0,p1,p2]; rep_i = protein of prev gene ──
        auto deriv = [&](const float* Y, float* D) {
            for (int i = 0; i < 3; ++i) {
                float r0  = Y[3 + ((i + 2) % 3)];
                float rep = r0 > 0.f ? r0 : 0.f;   // required: n non-integer; also flushes NaN→0 (hillLookup index safety)
                float hr  = hillDirect ? 1.f / (1.f + std::pow(rep, n)) : hillLookup(rep);
                D[i]     = -Y[i] + alpha * hr + alpha0;
                D[3 + i] = -beta * (Y[3 + i] - Y[i]);
            }
            D[0] += perturb;                    // PERTURB into gene-0 transcription
            D[0] +=  BIAS_EPS;                  // zero-sum self-start symmetry breaker
            D[1] += -0.5f * BIAS_EPS;
            D[2] += -0.5f * BIAS_EPS;
        };

        float y[6] = {m[0], m[1], m[2], p[0], p[1], p[2]};
        for (int s = 0; s < K; ++s) {
            coalescent::rk4<6>(y, h, deriv);
            bool bad = false;
            for (int j = 0; j < 6; ++j) { if (!std::isfinite(y[j])) bad = true; y[j] = clamp(y[j], 0.f, STATE_MAX); }
            if (bad) {
                reseedAsymmetric();
                for (int j = 0; j < 3; ++j) { y[j] = m[j]; y[3 + j] = p[j]; }
                resetGateMemory();
                break;
            }
            for (int k = 0; k < 3; ++k) {                  // gates on CENTERED value, per substep
                float c = y[3 + k] - pStar;
                if (lastCentered[k] <= GATE_LOW && c >= GATE_HIGH) gateGen[k].trigger(GATE_TIME);
                lastCentered[k] = c;
            }
        }
        for (int i = 0; i < 3; ++i) { m[i] = y[i]; p[i] = y[3 + i]; }

        // ── outputs: centered, soft-clipped proteins + per-sample phase gates ──
        for (int k = 0; k < 3; ++k)
            outputs[OUT1_OUTPUT + k].setVoltage(5.f * coalescent::fastTanh((p[k] - pStar) * OUT_GAIN));
        for (int k = 0; k < 3; ++k)
            outputs[GATE1_OUTPUT + k].setVoltage(gateGen[k].process(args.sampleTime) ? GATE_LEVEL : 0.f);

        // ── display scope: record the 3 centered proteins into a history ring at
        // ~25 Hz and publish a full copy — a slowed time plot ~10 s wide
        // (clean at LFO/clock rates; an audio-rate tone naturally aliases here) ──
        if (fs != lastDisplayFs) {
            lastDisplayFs = fs;
            dispDiv.setDivision(std::max(1, (int) std::round(fs / 25.f)));  // ~25 Hz → ~10 s window, calmer scroll
        }
        if (dispDiv.process()) {
            float mx = 0.f;
            for (int k = 0; k < 3; ++k) {
                float c = p[k] - pStar;
                liveScope.p[k][liveScope.head] = c;
                mx = std::max(mx, std::fabs(c));
            }
            liveScope.head = (liveScope.head + 1) % HIST_N;
            liveScope.peak = std::max(mx, liveScope.peak * 0.999f + 1e-4f);
            displaySnapshot.writable() = liveScope;
            displaySnapshot.publish();
        }
    }
};

// ── Panel layout (mm), shared by the widget and its label overlay ────────────
namespace opl {
    static const float KNOB_Y = 58.f, ATT_Y = 70.f, IN_Y = 84.f, OUT_Y = 100.f, GATE_Y = 114.f;
    static const float COL5[5] = {10.f, 22.8f, 35.6f, 48.4f, 61.2f};  // 5-wide rows
    static const float COL3[3] = {18.f, 35.6f, 53.2f};                // 3-wide output rows
    static const float ATT3[3] = {22.8f, 35.6f, 48.4f};               // attenuverters
}

// Display: the three protein levels over time — three ~120°-staggered traces on
// one centreline, so the phase chase reads directly. Recorded at ~25 Hz (~10 s
// wide); reads only the published history, never the integrator.
struct OperonScope : widget::TransparentWidget {
    Operon* module = nullptr;

    // Dark "screen" + bezel (house style — the SVG is just panel + rails + screws).
    void draw(const DrawArgs& args) override {
        nvgBeginPath(args.vg);
        nvgRoundedRect(args.vg, 0, 0, box.size.x, box.size.y, mm2px(2.f));
        nvgFillColor(args.vg, nvgRGB(0x0a, 0x07, 0x12));
        nvgFill(args.vg);
        nvgStrokeColor(args.vg, nvgRGB(0x3a, 0x2b, 0x4d));
        nvgStrokeWidth(args.vg, 1.2f);
        nvgStroke(args.vg);
        Widget::draw(args);
    }

    void drawLayer(const DrawArgs& args, int layer) override {
        if (layer != 1) return;
        const int N = Operon::HIST_N;
        static const Operon::OperonScope dummy;
        const Operon::OperonScope& sc = module ? module->displaySnapshot.consume() : dummy;
        float peak = std::max(sc.peak, 1e-3f);
        const float W = box.size.x, H = box.size.y;
        // three stacked lanes (title occupies the top ~20%), one protein each
        const float laneMid[3] = { H * 0.36f, H * 0.58f, H * 0.80f };
        const float laneH = H * 0.10f;     // ±peak maps to ±laneH within a lane
        const NVGcolor col[3] = { nvgRGB(0xd0, 0xb8, 0xff), nvgRGB(0xa8, 0x82, 0xf0), nvgRGB(0x82, 0x60, 0xd8) };

        nvgScissor(args.vg, 0, 0, W, H);
        // faint vertical guides — eyeball phase alignment across the three lanes
        for (int g = 1; g < 8; ++g) {
            float x = W * g / 8.f;
            nvgBeginPath(args.vg); nvgMoveTo(args.vg, x, 0); nvgLineTo(args.vg, x, H);
            nvgStrokeColor(args.vg, nvgRGBA(0xc8, 0xb0, 0xff, 0x16)); nvgStrokeWidth(args.vg, 1.f); nvgStroke(args.vg);
        }
        for (int k = 0; k < 3; ++k) {
            // faint lane baseline
            nvgBeginPath(args.vg); nvgMoveTo(args.vg, 0, laneMid[k]); nvgLineTo(args.vg, W, laneMid[k]);
            nvgStrokeColor(args.vg, nvgRGBA(0x6a, 0x50, 0x9e, 0x35)); nvgStrokeWidth(args.vg, 1.f); nvgStroke(args.vg);
            // trace, oldest (left) → newest (right)
            nvgBeginPath(args.vg);
            for (int i = 0; i < N; ++i) {
                int idx = (sc.head + i) % N;    // sc.head = oldest sample
                float x = W * i / (N - 1);
                float y = laneMid[k] - clamp(sc.p[k][idx] / peak, -1.2f, 1.2f) * laneH;
                if (i == 0) nvgMoveTo(args.vg, x, y); else nvgLineTo(args.vg, x, y);
            }
            nvgStrokeColor(args.vg, col[k]);
            nvgStrokeWidth(args.vg, 1.3f);
            nvgLineJoin(args.vg, NVG_ROUND);
            nvgStroke(args.vg);
        }
        nvgResetScissor(args.vg);

        // Title (top-left, DejaVuSans, letter-spaced, violet accent).
        std::shared_ptr<Font> font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
        if (font) {
            nvgFontFaceId(args.vg, font->handle);
            nvgFontSize(args.vg, mm2px(3.2f));
            nvgTextLetterSpacing(args.vg, mm2px(0.5f));
            nvgFillColor(args.vg, nvgRGBA(0xc8, 0xb0, 0xff, 0xcc));
            nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
            nvgText(args.vg, mm2px(2.6f), mm2px(2.2f), "OPERON", NULL);
            nvgTextLetterSpacing(args.vg, 0.f);
        }
    }
};

// Control labels (nanosvg ignores <text>, so labels live here on the panel).
struct OperonLabels : widget::Widget {
    void draw(const DrawArgs& args) override {
        std::shared_ptr<Font> font = APP->window->loadFont(asset::system("res/fonts/Nunito-Bold.ttf"));
        if (!font) return;
        nvgFontFaceId(args.vg, font->handle);
        nvgFillColor(args.vg, nvgRGB(0xe6, 0xe6, 0xf2));   // near-white, ~13:1 (house legibility spec)
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
        auto label = [&](float xmm, float ymm, const char* s, float sz = 2.0f) {
            nvgFontSize(args.vg, mm2px(sz * 1.72f));       // Nunito Bold, sized up for legibility
            nvgText(args.vg, mm2px(xmm), mm2px(ymm), s, nullptr);
        };
        using namespace opl;
        const char* knobL[5] = {"PITCH", "DRIVE", "HILL", "DECAY", "LEAK"};
        for (int i = 0; i < 5; ++i) label(COL5[i], KNOB_Y - 8.f, knobL[i]);   // match the 8mm knob-label gap of the other modules
        const char* inL[5] = {"V/OCT", "DRIVE", "HILL", "DECAY", "PERT"};
        for (int i = 0; i < 5; ++i) label(COL5[i], IN_Y - 5.5f, inL[i], 1.7f);
        const char* outL[3] = {"OUT1", "OUT2", "OUT3"};
        const char* gateL[3] = {"GATE1", "GATE2", "GATE3"};
        for (int i = 0; i < 3; ++i) {
            label(COL3[i], OUT_Y - 5.5f, outL[i], 1.8f);
            label(COL3[i], GATE_Y - 5.5f, gateL[i], 1.8f);
        }
    }
};

struct OperonWidget : ModuleWidget {

    OperonWidget(Operon* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Operon.svg")));
        addPanelLabels<OperonLabels>(this);

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(65.12f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f, 122.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(65.12f, 122.0f))));

        OperonScope* scope = new OperonScope();
        scope->module = module;
        scope->box.pos = mm2px(Vec(6.f, 8.f));
        scope->box.size = mm2px(Vec(59.12f, 38.f));
        addChild(scope);

        using namespace opl;
        // Knob row: PITCH, ALPHA, HILL, BETA, LEAK
        const int knobId[5] = {Operon::PITCH_PARAM, Operon::ALPHA_PARAM, Operon::HILL_PARAM,
                               Operon::BETA_PARAM, Operon::LEAK_PARAM};
        for (int i = 0; i < 5; ++i)
            addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(COL5[i], KNOB_Y)), module, knobId[i]));
        // Attenuverters under ALPHA/HILL/BETA
        const int attId[3] = {Operon::ALPHA_ATT_PARAM, Operon::HILL_ATT_PARAM, Operon::BETA_ATT_PARAM};
        for (int i = 0; i < 3; ++i)
            addParam(createParamCentered<Trimpot>(mm2px(Vec(ATT3[i], ATT_Y)), module, attId[i]));
        // Input row: VOCT, ALPHA, HILL, BETA, PERTURB
        const int inId[5] = {Operon::VOCT_INPUT, Operon::ALPHA_INPUT, Operon::HILL_INPUT,
                             Operon::BETA_INPUT, Operon::PERTURB_INPUT};
        for (int i = 0; i < 5; ++i)
            addInput(createInputCentered<PJ301MPort>(mm2px(Vec(COL5[i], IN_Y)), module, inId[i]));
        // Output rows: OUT1/2/3 then GATE1/2/3
        const int outId[3] = {Operon::OUT1_OUTPUT, Operon::OUT2_OUTPUT, Operon::OUT3_OUTPUT};
        const int gateId[3] = {Operon::GATE1_OUTPUT, Operon::GATE2_OUTPUT, Operon::GATE3_OUTPUT};
        for (int i = 0; i < 3; ++i) {
            addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(COL3[i], OUT_Y)), module, outId[i]));
            addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(COL3[i], GATE_Y)), module, gateId[i]));
        }
    }
};

Model* modelOperon = createModel<Operon, OperonWidget>("Operon");
