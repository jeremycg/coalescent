#include "plugin.hpp"
#include "dsp/rk4.hpp"
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
    // subspace and never oscillates. See §1 of the plan.
    float m[3] = {0.2f, 0.1f, 0.3f};
    float p[3] = {0.1f, 0.4f, 0.15f};
    float lastCentered[3] = {};          // upward-crossing detection on p[k] - pStar
    float pStar = 0.f;                   // cached symmetric fixed point (= mStar)
    float pStarA = -1.f, pStarN = -1.f, pStarL = -1.f;  // params pStar was solved for
    bool  pStarValid = false;
    dsp::PulseGenerator gateGen[3];

    // Hill-response LUT: 1/(1+repⁿ) precomputed over rep, cached on n (rebuilt only
    // when HILL changes, like pStar). Replaces the per-sample pow with a lookup —
    // ~1.5× on the RK4 derivative. A live HILL_INPUT (audio-rate n would thrash the
    // rebuild) falls back to direct pow. Sub-cent accurate even at n=8 (8192 pts).
    static constexpr int   HILL_LUT_N    = 8192;
    static constexpr float HILL_LUT_XMAX = 64.f;
    float hillLut[HILL_LUT_N + 1];
    float lutN = -1.f;
    bool  lutValid = false;

    void buildHillLut(float n) {
        for (int i = 0; i <= HILL_LUT_N; ++i)
            hillLut[i] = 1.f / (1.f + std::pow(HILL_LUT_XMAX * i / HILL_LUT_N, n));
    }
    inline float hillLookup(float rep) const {   // rep already >= 0
        if (rep >= HILL_LUT_XMAX) return hillLut[HILL_LUT_N];
        float f = rep * (HILL_LUT_N / HILL_LUT_XMAX);
        int   i = (int) f;
        return hillLut[i] + (f - i) * (hillLut[i + 1] - hillLut[i]);
    }

    // ─── Display snapshot (lock-free double buffer; UI-only) ────────────────
    struct OperonSnapshot {
        float p[3] = {};                 // centered protein snapshot
        float peak = 1e-3f;              // smoothed normaliser
    };
    OperonSnapshot snapshots[2];
    std::atomic<int> snapshotIndex{0};
    dsp::ClockDivider dispDiv;
    float lastDisplayFs = 0.f;

    // ─── Tunable constants (final values set at M7) ─────────────────────────
    static constexpr float RATE_CAL     = 12.46f;  // measured default period (tools/stability/operon.cpp) → C4 at 0 V
    static constexpr float HSUB_MAX     = 0.05f;
    static constexpr int   MIN_SUB      = 2;   // profiled: K=2 holds 0-cent + bounded at default (adaptive wants 2; the 4-floor doubled the pow count)
    static constexpr int   MAX_SUB      = 64;
    static constexpr float PITCH_TOTAL_MIN = -8.f, PITCH_TOTAL_MAX = 8.f;
    static constexpr float OUT_GAIN     = 0.3f;   // set offline: default RMS ~3.2V, peaks lightly saturate, high drive grits
    static constexpr float STATE_MAX    = 1e3f;
    static constexpr float PERTURB_GAIN = 0.5f;
    static constexpr float BIAS_EPS     = 1e-6f;   // zero-sum self-start symmetry breaker
    static constexpr float GATE_LEVEL   = 10.f, GATE_TIME = 1e-3f;
    static constexpr float GATE_LOW     = -0.05f, GATE_HIGH = 0.05f;  // centered units
    static constexpr int   CENTER_ITERS = 16;
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
        lutValid = false; lutN = -1.f;
        for (int k = 0; k < 3; ++k) { lastCentered[k] = 0.f; gateGen[k].reset(); }
    }

    // Symmetric fixed point pStar (= mStar): the monotonic g(x) = x - alpha/(1+x^n)
    // - alpha0 has one root; bisection. max(mid,0) before pow — n is non-integer.
    static float solvePStar(float alpha, float n, float alpha0) {
        float lo = 0.f, hi = std::max(1.f, alpha + alpha0 + 1.f);
        for (int i = 0; i < CENTER_ITERS; ++i) {
            float mid = 0.5f * (lo + hi);
            float g = mid - alpha / (1.f + std::pow(std::max(mid, 0.f), n)) - alpha0;
            if (g > 0.f) hi = mid; else lo = mid;
        }
        return 0.5f * (lo + hi);
    }

    void process(const ProcessArgs& args) override {
        const float fs = args.sampleRate;

        // ── 3.1 read + clamp params (per-sample; CV added with separate depths) ──
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
        float perturb = inputs[PERTURB_INPUT].isConnected()
                      ? inputs[PERTURB_INPUT].getVoltage() * PERTURB_GAIN : 0.f;

        // ── 3.2 symmetric fixed point for centering (cached on alpha/n/alpha0) ──
        if (!pStarValid || alpha != pStarA || n != pStarN || alpha0 != pStarL) {
            pStar = solvePStar(alpha, n, alpha0);
            pStarA = alpha; pStarN = n; pStarL = alpha0; pStarValid = true;
        }

        // Hill LUT: rebuild only on n change (control-rate for a knob). With HILL
        // being CV'd, n can move audio-rate, so use direct pow instead of thrashing.
        const bool hillDirect = inputs[HILL_INPUT].isConnected();
        if (!hillDirect && (!lutValid || n != lutN)) {
            buildHillLut(n); lutN = n; lutValid = true;
        }

        // ── 3.4 pitch = simulation speed; adaptive substepping ──
        float pitchTotal = clamp(params[PITCH_PARAM].getValue() + inputs[VOCT_INPUT].getVoltage(),
                                 PITCH_TOTAL_MIN, PITCH_TOTAL_MAX);
        float pitchHz = dsp::FREQ_C4 * dsp::approxExp2_taylor5(pitchTotal);
        float dtau = RATE_CAL * pitchHz / fs;
        dtau = std::min(dtau, HSUB_MAX * MAX_SUB);          // hard guard: keep h <= HSUB_MAX
        int   K = clamp((int) std::ceil(dtau / HSUB_MAX), MIN_SUB, MAX_SUB);
        float h = dtau / K;

        // ── 3.3 derivative over y = [m0,m1,m2,p0,p1,p2]; rep_i = protein of prev gene ──
        auto deriv = [&](const float* Y, float* D) {
            for (int i = 0; i < 3; ++i) {
                float rep = std::max(Y[3 + ((i + 2) % 3)], 0.f);   // required: n non-integer
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

        // ── 3.5 outputs: centered, soft-clipped proteins + per-sample phase gates ──
        for (int k = 0; k < 3; ++k)
            outputs[OUT1_OUTPUT + k].setVoltage(5.f * coalescent::fastTanh((p[k] - pStar) * OUT_GAIN));
        for (int k = 0; k < 3; ++k)
            outputs[GATE1_OUTPUT + k].setVoltage(gateGen[k].process(args.sampleTime) ? GATE_LEVEL : 0.f);

        // ── display snapshot (~45 Hz, lock-free double buffer; UI-only) ──
        if (fs != lastDisplayFs) {
            lastDisplayFs = fs;
            dispDiv.setDivision(std::max(1, (int) std::round(fs / 45.f)));
        }
        if (dispDiv.process()) {
            int w = 1 - snapshotIndex.load(std::memory_order_relaxed);
            OperonSnapshot& s = snapshots[w];
            float mx = 0.f;
            for (int k = 0; k < 3; ++k) { s.p[k] = p[k] - pStar; mx = std::max(mx, std::fabs(s.p[k])); }
            const OperonSnapshot& prev = snapshots[1 - w];
            s.peak = std::max(mx, prev.peak * 0.95f + 1e-4f);   // slow decay, avoids /0
            snapshotIndex.store(w, std::memory_order_release);
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

// Gene-ring display. M1: three static dim nodes + repression arrows so the panel
// reads as a ring; the live, peak-normalised glow driven by the protein snapshot
// lands in M6. Reads only the published snapshot — never the integrator.
struct OperonRing : widget::TransparentWidget {
    Operon* module = nullptr;
    std::shared_ptr<Font> font;
    float glow[3] = {};   // per-node smoothed activity (breathe, not flicker at audio rate)

    // Dark "screen" + bezel, drawn by the widget (house style — the SVG is just
    // panel + rails + screws).
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
        Vec ctr = box.size.div(2);
        float R = std::min(box.size.x, box.size.y) * 0.32f;
        const float ang[3] = {-(float)M_PI / 2.f,
                              -(float)M_PI / 2.f + 2.f * (float)M_PI / 3.f,
                              -(float)M_PI / 2.f + 4.f * (float)M_PI / 3.f};
        Vec node[3];
        for (int k = 0; k < 3; ++k)
            node[k] = ctr.plus(Vec(std::cos(ang[k]), std::sin(ang[k])).mult(R));

        int read = module ? module->snapshotIndex.load(std::memory_order_acquire) : 0;
        const Operon::OperonSnapshot dummy;
        const Operon::OperonSnapshot& snap = module ? module->snapshots[read] : dummy;
        float peak = std::max(snap.peak, 1e-3f);

        // Repression arrows k -> (k+1)%3 (dim, static).
        for (int k = 0; k < 3; ++k) {
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, node[k].x, node[k].y);
            nvgLineTo(args.vg, node[(k + 1) % 3].x, node[(k + 1) % 3].y);
            nvgStrokeColor(args.vg, nvgRGBA(0x6a, 0x50, 0x9e, 0x60));
            nvgStrokeWidth(args.vg, 1.f);
            nvgStroke(args.vg);
        }
        // Nodes — activity (|centered| / peak) drives glow; rest reads dim (violet).
        for (int k = 0; k < 3; ++k) {
            float activity = clamp(std::fabs(snap.p[k]) / peak, 0.f, 1.f);
            glow[k] += (activity - glow[k]) * 0.06f;   // ~250 ms smoothing: breathe, not flicker
            float a = glow[k];
            float rad = mm2px(1.6f) + mm2px(2.2f) * a;
            nvgBeginPath(args.vg);
            nvgCircle(args.vg, node[k].x, node[k].y, rad);
            nvgFillColor(args.vg, nvgRGBA(0xa8, 0x78, 0xff, (int)(0x30 + 0xb0 * a)));
            nvgFill(args.vg);
        }

        // Title (top-left, DejaVuSans, letter-spaced, violet accent).
        if (!font) font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
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

struct OperonWidget : ModuleWidget {
    std::shared_ptr<Font> font;

    OperonWidget(Operon* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Operon.svg")));

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(65.12f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f, 122.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(65.12f, 122.0f))));

        OperonRing* ring = new OperonRing();
        ring->module = module;
        ring->box.pos = mm2px(Vec(6.f, 8.f));
        ring->box.size = mm2px(Vec(59.12f, 38.f));
        addChild(ring);

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

    // Control labels (nanosvg ignores <text>, so labels live here on the panel).
    void draw(const DrawArgs& args) override {
        ModuleWidget::draw(args);
        if (!font) font = APP->window->loadFont(asset::system("res/fonts/Nunito-Bold.ttf"));
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
        for (int i = 0; i < 5; ++i) label(COL5[i], KNOB_Y - 6.f, knobL[i]);
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

Model* modelOperon = createModel<Operon, OperonWidget>("Operon");
