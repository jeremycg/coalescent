#include "plugin.hpp"
#include "integrator.hpp"
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
// substepping) is the same engine as Axon, with one extra equation — the
// factoring the Axon plan called for.
//
// Polyphonic: up to 16 independent HR voices, voice count from V/OCT (falling
// back to TRIG). The display traces every active voice on its own hue, stepped
// across a narrow band around Soma's amber accent.

static const int TRAIL = 512;   // length of the (x,z) phase-portrait trail

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
    float xx[MAX_POLY], yy[MAX_POLY], zz[MAX_POLY];  // HR state per voice
    float trigPulse[MAX_POLY] = {};
    dsp::SchmittTrigger trigIn[MAX_POLY];
    dsp::SchmittTrigger syncIn[MAX_POLY];     // SYNC edge detector (hard reset to rest)
    dsp::SchmittTrigger spikeDet[MAX_POLY];
    dsp::PulseGenerator spikeGen[MAX_POLY];
    dsp::TRCFilter<float> dcBlock[MAX_POLY];
    // Anti-aliasing decimators (one per voice; only the active factor is used).
    // The whole output chain — DC-block then tanh soft-clip — runs oversampled
    // and is decimated here, so the sharp spikes and the tanh nonlinearity are
    // band-limited before reaching the output sample.
    dsp::Decimator<4, 8> decim4[MAX_POLY];
    dsp::Decimator<8, 8> decim8[MAX_POLY];
    int   oversampleMode = 1;         // 0 = off, 1 = ×4, 2 = ×8 (right-click menu)
    int   channels = 1;               // active voice count
    float lastFs = 0.f;
    int   lastOs = 0;                 // detects oversample change to refresh coefficients
    float trigDecay = 0.f;            // cached per-sample TRIG pulse decay (fn of fs)
    float antiDenorm = 1e-18f;        // alternating sub-LSB dither, anti-denormal on dcBlock

    // ─── Display trail (x,z phase portrait) ─────────────────────────────────
    // Double-buffered, lock-free snapshot (see Axon): the audio thread fills the
    // back buffer and flips dispBuf with a release store; the UI reads the front
    // buffer after an acquire load, with the head index and active-voice count
    // carried alongside so trails and head dots stay coherent.
    float trailX[MAX_POLY][TRAIL] = {}, trailZ[MAX_POLY][TRAIL] = {};
    int   trailIdx = 0, trailDecim = 0;
    float dispX[2][MAX_POLY][TRAIL] = {}, dispZ[2][MAX_POLY][TRAIL] = {};
    int   dispHead[2] = {};
    int   dispChannels[2] = {1, 1};
    std::atomic<int> dispBuf{0};
    int   dispClock = 0;

    // ─── Fixed HR parameters ────────────────────────────────────────────────
    static constexpr float A = 1.f, B = 3.f, C = 1.f, D = 5.f, XR = -1.6f;

    // ─── Tunable constants (RATE_CAL from tools/soma_stability_test.cpp) ─────
    static constexpr float RATE_CAL    = 14.925501f; // within-burst spike period at default → C4 at 0 V
    static constexpr float HSUB_MAX    = 0.05f;
    static constexpr int   MIN_SUB     = 4;
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
        configParam(BURST_PARAM, std::log2(R_MIN), std::log2(R_MAX), std::log2(0.006f),
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

        for (int c = 0; c < MAX_POLY; c++) { xx[c] = -1.6f; yy[c] = -11.8f; zz[c] = 2.f; }
    }

    void onReset() override {
        oversampleMode = 1;   // restore default anti-aliasing (×4) on Initialize
        for (int c = 0; c < MAX_POLY; c++) {
            xx[c] = -1.6f; yy[c] = -11.8f; zz[c] = 2.f; trigPulse[c] = 0.f;
            trigIn[c].reset(); syncIn[c].reset(); spikeDet[c].reset(); spikeGen[c].reset();
            dcBlock[c].reset(); decim4[c].reset(); decim8[c].reset();
        }
    }

    json_t* dataToJson() override {
        json_t* root = json_object();
        json_object_set_new(root, "oversample", json_integer(oversampleMode));
        return root;
    }

    void dataFromJson(json_t* root) override {
        if (json_t* j = json_object_get(root, "oversample"))
            oversampleMode = (int) json_integer_value(j);
    }

    // HR derivatives — the FHN f() with one extra (slow) line, as the Axon plan
    // anticipated. Same RK4 step below, now over three variables.
    static inline void f(float x, float y, float z, float Itot, float r, float s,
                         float& dx, float& dy, float& dz) {
        dx = y - A*x*x*x + B*x*x - z + Itot;
        dy = C - D*x*x - y;
        dz = r * (s * (x - XR) - z);
    }

    // Shared neuron::rk4 stepper over the three HR variables — arithmetic
    // identical to the original hand-written step (see integrator.hpp).
    static inline void rk4(float& x, float& y, float& z, float h,
                           float Itot, float r, float s) {
        float st[3] = {x, y, z};
        neuron::rk4<3>(st, h, [&](const float* Y, float* d) {
            f(Y[0], Y[1], Y[2], Itot, r, s, d[0], d[1], d[2]);
        });
        x = st[0];
        y = st[1];
        z = st[2];
    }

    void process(const ProcessArgs& args) override {
        const float fs = args.sampleRate;
        const int os = oversampleMode == 0 ? 1 : (oversampleMode == 1 ? 4 : 8);

        // SR/oversample-derived coefficients, recomputed only when fs or the
        // oversample factor changes. The DC-blocker runs in the oversampled
        // domain, so its cutoff is referenced to the oversampled rate (fs·os) to
        // keep the real-time corner at ~20 Hz. Caching trigDecay keeps a
        // transcendental off the per-sample path.
        if (fs != lastFs || os != lastOs) {
            for (int c = 0; c < MAX_POLY; c++) dcBlock[c].setCutoffFreq(20.f / (fs * os));
            trigDecay = std::exp(-1.f / (TRIG_TAU_MS * 1e-3f * fs));
            if (os != lastOs)
                for (int c = 0; c < MAX_POLY; c++) { decim4[c].reset(); decim8[c].reset(); }
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

        for (int c = 0; c < channels; c++) {
            // ── params → physics (per voice) ──
            float I = Ibase
                    + inputs[CURRENT_INPUT].getPolyVoltage(c) * Iatt * CV_DEPTH;
            // BURST is log₂(r); CV adds in the log domain (a multiplicative nudge on r).
            float rLog = rLogBase
                       + inputs[BURST_INPUT].getPolyVoltage(c) * rAtt * CV_DEPTH;
            float r = clamp(std::exp2(rLog), R_MIN, R_MAX);

            // ── hard sync: a rising edge re-seeds the voice at the rest fixed
            // point (all three state variables). The discontinuity is the sync. ──
            if (syncIn[c].process(inputs[SYNC_INPUT].getPolyVoltage(c), 0.1f, 1.f)) {
                xx[c] = -1.6f; yy[c] = -11.8f; zz[c] = 2.f;
            }

            // ── excitable trigger: rising edge injects a decaying current pulse ──
            if (trigIn[c].process(inputs[TRIG_INPUT].getPolyVoltage(c), 0.1f, 1.f))
                trigPulse[c] = TRIG_AMP;
            float Itot = I + trigPulse[c];
            trigPulse[c] *= trigDecay;
            if (std::fabs(trigPulse[c]) < 1e-30f) trigPulse[c] = 0.f;

            // ── pitch = simulation speed (open-loop; tracks within-burst spike rate) ──
            float pitchHz = dsp::FREQ_C4 * std::exp2(
                                pitchKnob + inputs[VOCT_INPUT].getPolyVoltage(c));
            float dtau = RATE_CAL * pitchHz / fs;
            float subTau = dtau / os;                             // advance per oversampled sample
            int   K = clamp((int) std::ceil(subTau / HSUB_MAX), MIN_SUB, MAX_SUB);
            float h = subTau / K;

            // ── render os oversampled samples (full output chain each), decimate ──
            float x = xx[c], y = yy[c], z = zz[c];
            float osBuf[8];
            for (int o = 0; o < os; o++) {
                for (int k = 0; k < K; k++)
                    rk4(x, y, z, h, Itot, r, s);

                // Backstop.
                if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
                    x = -1.6f; y = -11.8f; z = 2.f;
                }
                x = clamp(x, -STATE_MAX, STATE_MAX);
                y = clamp(y, -STATE_MAX, STATE_MAX);
                z = clamp(z, -STATE_MAX, STATE_MAX);

                if (spikeDet[c].process(x, 0.f, SPIKE_THRESH))
                    spikeGen[c].trigger(1e-3f);

                dcBlock[c].process(x + antiDenorm);
                osBuf[o] = 5.f * std::tanh(dcBlock[c].highpass() * OUT_GAIN);
            }
            xx[c] = x; yy[c] = y; zz[c] = z;

            float audio = os == 1 ? osBuf[0]
                        : (os == 4 ? decim4[c].process(osBuf) : decim8[c].process(osBuf));

            // ── outputs ──
            outputs[OUT_OUTPUT].setVoltage(audio, c);
            outputs[SPIKE_OUTPUT].setVoltage(spikeGen[c].process(args.sampleTime) ? 10.f : 0.f, c);
            // Z is the slow adaptation variable — a burst-envelope CV (not high-passed).
            outputs[Z_OUTPUT].setVoltage(clamp((z - Z_CENTER) * Z_GAIN, -5.f, 5.f), c);
        }

        // ── feed the (x,z) phase-portrait trails (all voices; decimated, bursts are slow) ──
        if (++trailDecim >= 6) {
            trailDecim = 0;
            for (int c = 0; c < channels; c++) {
                trailX[c][trailIdx] = xx[c];
                trailZ[c][trailIdx] = zz[c];
            }
            trailIdx = (trailIdx + 1) % TRAIL;
        }
        if (++dispClock >= (int)(fs / 45.f)) {       // publish display snapshot ~45 Hz
            dispClock = 0;
            int next = 1 - dispBuf.load(std::memory_order_relaxed);
            for (int c = 0; c < channels; c++) {
                std::copy(trailX[c], trailX[c] + TRAIL, dispX[next][c]);
                std::copy(trailZ[c], trailZ[c] + TRAIL, dispZ[next][c]);
            }
            dispHead[next] = trailIdx;
            dispChannels[next] = channels;
            dispBuf.store(next, std::memory_order_release);
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
    std::shared_ptr<Font> font;

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
            int b  = module ? module->dispBuf.load(std::memory_order_acquire) : 0;
            int nv = module ? module->dispChannels[b] : 1;

            // z-nullcline: z = s·(x − x_R).
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, X(XMIN), Y(s * (XMIN - Soma::XR)));
            nvgLineTo(args.vg, X(XMAX), Y(s * (XMAX - Soma::XR)));
            nvgStrokeColor(args.vg, nvgRGBA(0x90, 0x50, 0x70, 0x55));
            nvgStrokeWidth(args.vg, 1.f);
            nvgStroke(args.vg);

            if (module) {
                int idx = module->dispHead[b];   // coherent with arrays
                float trailA = clamp(204.f - (nv - 1) * 6.f, 112.f, 204.f);  // 0x70..0xcc
                nvgLineCap(args.vg, NVG_ROUND);
                for (int v = 0; v < nv; v++) {
                    const float* dx = module->dispX[b][v];
                    const float* dz = module->dispZ[b][v];
                    float hue = voiceHue(v, nv);
                    for (int k = 1; k < TRAIL; k++) {
                        int i0 = (idx + k - 1) % TRAIL;
                        int i1 = (idx + k) % TRAIL;
                        float alpha = (float) k / TRAIL;
                        nvgBeginPath(args.vg);
                        nvgMoveTo(args.vg, X(dx[i0]), Y(dz[i0]));
                        nvgLineTo(args.vg, X(dx[i1]), Y(dz[i1]));
                        nvgStrokeColor(args.vg, nvgHSLA(hue, 0.95f, 0.60f, (int)(alpha * trailA)));
                        nvgStrokeWidth(args.vg, 1.6f);
                        nvgStroke(args.vg);
                    }
                    int newest = (idx + TRAIL - 1) % TRAIL;
                    float hx = X(dx[newest]), hy = Y(dz[newest]);
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

            if (!font)
                font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
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
                    char buf[8];
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

struct SomaWidget : ModuleWidget {
    std::shared_ptr<Font> font;

    void draw(const DrawArgs& args) override {
        ModuleWidget::draw(args);
        if (!font)
            font = APP->window->loadFont(asset::system("res/fonts/Nunito-Bold.ttf"));
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

    SomaWidget(Soma* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Soma.svg")));

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
            {"Off", "×4 oversampling", "×8 oversampling"}, &m->oversampleMode));
    }
};

Model* modelSoma = createModel<Soma, SomaWidget>("Soma");
