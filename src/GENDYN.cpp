#include "plugin.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>

static const int MAX_N = 64;

struct GENDYN : Module {

    enum ParamId {
        N_PARAM,
        SCALE_AMP_PARAM,
        SCALE_DUR_PARAM,
        B_AMP_PARAM,
        B_DUR_CENTER_PARAM,
        B_DUR_WIDTH_PARAM,
        DISTRIBUTION_PARAM,
        SCALE_AMP_ATT_PARAM,
        SCALE_DUR_ATT_PARAM,
        B_AMP_ATT_PARAM,
        B_DUR_ATT_PARAM,
        PERSIST_PARAM,
        LOCK_PARAM,
        PARAMS_LEN
    };

    enum InputId {
        SCALE_AMP_INPUT,
        SCALE_DUR_INPUT,
        B_AMP_INPUT,
        B_DUR_INPUT,
        INPUTS_LEN
    };

    enum OutputId {
        AUDIO_OUTPUT,
        CYCLE_TRIG_OUTPUT,
        FREQ_CV_OUTPUT,
        OUTPUTS_LEN
    };

    enum LightId {
        LIGHTS_LEN
    };

    float amp[MAX_N];               // breakpoint amplitudes (volts, range ±bAmp)
    float dur[MAX_N];               // breakpoint durations  (samples)
    float step_amp[MAX_N];          // persistent primary walk: amplitude step
    float step_dur[MAX_N];          // persistent primary walk: duration step
    int   current_breakpoint = 0;
    int   current_sample     = 0;
    float current_amp        = 0.f; // amplitude at the start of the current segment
    float target_amp         = 0.f; // amplitude at the end   of the current segment
    int   current_dur        = 1;   // length of the current segment in samples
    float sum_dur            = 1.f; // cached sum of dur[0..N-1], updated per cycle
    float freq_cv            = 0.f; // cached FREQ output voltage, updated with sum_dur
    float norm_k             = 1.f; // playback duration scale; LOCK mode sets it so
                                    // the cycle length is exactly sampleRate/centerFreq
    float dur_err            = 0.f; // error-diffusion remainder for int segment lengths
    int   last_N             = -1;  // -1 forces reinit on the first process() call

    dsp::PulseGenerator cycleTrigger;

    // Display snapshot: ~45 Hz the audio thread publishes the breakpoints into a
    // lock-free double buffer (fill the back buffer, flip dispBuf with a release
    // store); the UI reads the front buffer after an acquire load. The walk drifts
    // slowly, so the polygon morphs smoothly, and the read no longer races the
    // per-cycle writer.
    float dispAmp[2][MAX_N] = {};
    float dispDur[2][MAX_N] = {};
    int   dispN[2]    = {13, 13};
    float dispBAmp[2] = {0.8f, 0.8f};
    std::atomic<int> dispBuf{0};
    int   dispClock = 0;

    GENDYN() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);

        configParam(N_PARAM,            2.f,   64.f,   13.f,  "Breakpoints (N)");
        getParamQuantity(N_PARAM)->snapEnabled = true;

        // Defaults retuned 2026-06 for the second-order walk: steps persist
        // across cycles, so net drift is ~8x faster at equal scale. 0.35x the
        // old first-order values preserves the same evolution timescale.
        configParam(SCALE_AMP_PARAM,    0.f,   1.f,    0.002f, "Amplitude step scale");
        configParam(SCALE_DUR_PARAM,    0.f,   1.f,    0.0035f,"Duration step scale");
        configParam(B_AMP_PARAM,        0.f,   1.f,    0.8f,  "Amplitude barrier");
        configParam(B_DUR_CENTER_PARAM, 20.f,  5000.f, 220.f, "Center frequency", " Hz");
        configParam(B_DUR_WIDTH_PARAM,  0.f,   1.f,    0.4f,  "Frequency barrier width");

        configParam(DISTRIBUTION_PARAM, 0.f,   3.f,    3.f,
            "Distribution (0=Cauchy  1=Gauss  2=Uniform  3=Logistic)");
        getParamQuantity(DISTRIBUTION_PARAM)->snapEnabled = true;

        configParam(PERSIST_PARAM,      0.f,   1.f,    0.3f,  "Glide persistence", "%", 0.f, 100.f);

        configSwitch(LOCK_PARAM,        0.f,   1.f,    0.f,   "Pitch lock (normalize durations)",
            {"Off", "On"});

        configParam(SCALE_AMP_ATT_PARAM,-1.f,  1.f,    0.f,   "Scale Amp attenuverter");
        configParam(SCALE_DUR_ATT_PARAM,-1.f,  1.f,    0.f,   "Scale Dur attenuverter");
        configParam(B_AMP_ATT_PARAM,    -1.f,  1.f,    0.f,   "B Amp attenuverter");
        configParam(B_DUR_ATT_PARAM,    -1.f,  1.f,    0.f,   "B Dur attenuverter");

        configInput(SCALE_AMP_INPUT, "Scale Amp CV");
        configInput(SCALE_DUR_INPUT, "Scale Dur CV");
        configInput(B_AMP_INPUT,     "Amplitude barrier CV");
        configInput(B_DUR_INPUT,     "Duration barrier CV");

        configOutput(AUDIO_OUTPUT,      "Audio");
        configOutput(CYCLE_TRIG_OUTPUT, "Cycle trigger");
        configOutput(FREQ_CV_OUTPUT,    "Frequency 1V/oct");
    }

    // ── DSP helpers ───────────────────────────────────────────────────────────

    // Reflect x into [lo, hi] by folding at each barrier (Serra eq. 3).
    // O(1) triangle fold, equivalent to mirroring repeatedly at the barriers;
    // constant-time even for far-out x (Cauchy draws have heavy tails).
    // Guards against a degenerate zero-width range.
    static float reflect(float x, float lo, float hi) {
        if (hi <= lo) return lo;
        const float range = hi - lo;
        float y = std::fmod(x - lo, 2.f * range);
        if (y < 0.f) y += 2.f * range;
        return lo + (y > range ? 2.f * range - y : y);
    }

    // Inverse-CDF samplers. `scale` is the distribution's spread parameter.
    static float cauchySample(float scale) {
        float u = clamp(rack::random::uniform(), 0.0001f, 0.9999f);
        return scale * std::tan(M_PI * (u - 0.5f));
    }

    static float gaussianSample(float scale) {
        // Box-Muller; clamp u1 away from 0 to avoid log(0).
        float u1 = clamp(rack::random::uniform(), 1e-6f, 1.f);
        float u2 = rack::random::uniform();
        return scale * std::sqrt(-2.f * std::log(u1)) * std::cos(2.f * M_PI * u2);
    }

    static float uniformSample(float scale) {
        return scale * (2.f * rack::random::uniform() - 1.f);
    }

    static float logisticSample(float scale) {
        float u = clamp(rack::random::uniform(), 0.0001f, 0.9999f);
        return scale * std::log(u / (1.f - u));
    }

    float drawSample(int dist, float scale) {
        switch (dist) {
            case 1:  return gaussianSample(scale);
            case 2:  return uniformSample(scale);
            case 3:  return logisticSample(scale);
            default: return cauchySample(scale);
        }
    }

    // ── Oscillator helpers ────────────────────────────────────────────────────

    // Randomise all breakpoints across the full barrier range and reset the
    // oscillator to breakpoint 0. Called on N change and at startup.
    void initBreakpoints(int N, float bAmp, float bDurMin, float bDurMax) {
        sum_dur = 0.f;
        for (int i = 0; i < N; i++) {
            amp[i] = (2.f * rack::random::uniform() - 1.f) * bAmp;
            dur[i] = bDurMin + rack::random::uniform() * (bDurMax - bDurMin);
            step_amp[i] = 0.f;
            step_dur[i] = 0.f;
            sum_dur += dur[i];
        }
        current_breakpoint = 0;
        current_sample     = 0;
        current_amp        = 0.f;
        target_amp         = amp[0];
        current_dur        = std::max(1, (int)dur[0]);
        dur_err            = 0.f;
    }

    // Apply one cycle of the second-order random walk to all breakpoints
    // (Serra 1993 eq. 2 + 4 with persistent steps; Hoffmann 2023; same
    // factoring as SC Gendy2). The persistent step is a normalized walk in
    // [-1, 1]: each cycle a draw of `spread` nudges it, and the step moves
    // the breakpoint scaled by scale*barrier. The spread sets the glide
    // persistence — how many cycles the step keeps its direction (~1 cycle
    // at spread 1, ~16 at 0.25, hundreds at 0.01). The step state is
    // scale-independent, making SCALE a pure gain: CV modulation rescales
    // motion instantly and reversibly instead of folding or starving the
    // accumulated steps. scale*barrier caps the per-cycle move: gainDur is
    // the maximum glissando rate, gainAmp the rate of timbral change.
    void runCycleUpdate(int N, float scaleAmp, float scaleDur,
                        float bAmp, float bDurMin, float bDurMax,
                        int dist, float spread) {
        const float pDurHW  = (bDurMax - bDurMin) * 0.5f;
        const float gainAmp = scaleAmp * bAmp;
        const float gainDur = scaleDur * pDurHW;
        sum_dur = 0.f;
        for (int i = 0; i < N; i++) {
            step_amp[i] = reflect(step_amp[i] + drawSample(dist, spread), -1.f, 1.f);
            amp[i]      = reflect(amp[i] + gainAmp * step_amp[i], -bAmp, bAmp);

            step_dur[i] = reflect(step_dur[i] + drawSample(dist, spread), -1.f, 1.f);
            dur[i]      = reflect(dur[i] + gainDur * step_dur[i], bDurMin, bDurMax);
            sum_dur += dur[i];
        }
    }

    // 1V/oct CV (C4 = 0V) for the frequency sampleRate / period.
    float computeFreqCV(float sampleRate, float period) const {
        if (period <= 0.f) return 0.f;
        return clamp(std::log2(sampleRate / period / dsp::FREQ_C4), -5.f, 5.f);
    }

    // LOCK mode (SC Gendy3-style): scale playback durations so each cycle
    // sums to exactly sampleRate/centerFreq. The duration *walk* stays in
    // its barriers untouched — relative durations keep evolving (timbre)
    // while pitch holds — so toggling LOCK is clean and stateless.
    void updateNormAndFreq(bool lockPitch, float sampleRate, float centerFreq) {
        norm_k = 1.f;
        if (lockPitch && sum_dur > 0.f)
            norm_k = (sampleRate / centerFreq) / sum_dur;
        freq_cv = computeFreqCV(sampleRate, sum_dur * norm_k);
    }

    // ── Main DSP loop ─────────────────────────────────────────────────────────

    void process(const ProcessArgs& args) override {

        // ── Read parameters ───────────────────────────────────────────────────
        int N = clamp((int)params[N_PARAM].getValue(), 2, MAX_N);

        float centerFreq = clamp(params[B_DUR_CENTER_PARAM].getValue(), 20.f, 5000.f);

        float scaleAmp = params[SCALE_AMP_PARAM].getValue();
        if (inputs[SCALE_AMP_INPUT].isConnected())
            scaleAmp += inputs[SCALE_AMP_INPUT].getVoltage()
                      * params[SCALE_AMP_ATT_PARAM].getValue() * 0.1f;
        scaleAmp = clamp(scaleAmp, 0.f, 1.f);

        float scaleDur = params[SCALE_DUR_PARAM].getValue();
        if (inputs[SCALE_DUR_INPUT].isConnected())
            scaleDur += inputs[SCALE_DUR_INPUT].getVoltage()
                      * params[SCALE_DUR_ATT_PARAM].getValue() * 0.1f;
        scaleDur = clamp(scaleDur, 0.f, 1.f);

        float bAmp = params[B_AMP_PARAM].getValue();
        if (inputs[B_AMP_INPUT].isConnected())
            bAmp += inputs[B_AMP_INPUT].getVoltage()
                  * params[B_AMP_ATT_PARAM].getValue() * 0.1f;
        bAmp = clamp(bAmp, 0.f, 1.f);

        float bDurWidth = params[B_DUR_WIDTH_PARAM].getValue();
        if (inputs[B_DUR_INPUT].isConnected())
            bDurWidth += inputs[B_DUR_INPUT].getVoltage()
                       * params[B_DUR_ATT_PARAM].getValue() * 0.1f;
        bDurWidth = clamp(bDurWidth, 0.f, 1.f);

        // Average samples-per-breakpoint at the target frequency, ±bDurWidth
        // fraction of that centre value. bDurWidth=0 collapses to a 1-sample
        // window, giving a fixed pitch.
        const float durCenter = args.sampleRate / (centerFreq * (float)N);
        const float halfWidth = bDurWidth * durCenter;
        float bDurMin = std::max(1.f, durCenter - halfWidth);
        float bDurMax = durCenter + halfWidth;
        if (bDurMax <= bDurMin) bDurMax = bDurMin + 1.f;

        int dist = clamp((int)params[DISTRIBUTION_PARAM].getValue(), 0, 3);

        const bool lockPitch = params[LOCK_PARAM].getValue() > 0.5f;

        // ── Reinit on N change ────────────────────────────────────────────────
        if (N != last_N) {
            initBreakpoints(N, bAmp, bDurMin, bDurMax);
            updateNormAndFreq(lockPitch, args.sampleRate, centerFreq);
            current_dur = std::max(1, (int)(dur[0] * norm_k));
            last_N = N;
        }

        // ── Generate output sample ────────────────────────────────────────────
        const float t      = (float)current_sample / (float)current_dur;
        const float output = current_amp + t * (target_amp - current_amp);
        outputs[AUDIO_OUTPUT].setVoltage(output * 5.f);

        // ── Auxiliary outputs ─────────────────────────────────────────────────
        outputs[CYCLE_TRIG_OUTPUT].setVoltage(
            cycleTrigger.process(args.sampleTime) ? 10.f : 0.f);

        // Instantaneous frequency (1V/oct, C4=0V). sum_dur — and with it
        // freq_cv — only changes at init and cycle updates, so the log2 is
        // computed there rather than per sample.
        outputs[FREQ_CV_OUTPUT].setVoltage(freq_cv);

        // ── Advance oscillator state ──────────────────────────────────────────
        if (++current_sample >= current_dur) {
            current_amp        = target_amp;
            current_sample     = 0;
            current_breakpoint = (current_breakpoint + 1) % N;

            if (current_breakpoint == 0) {
                // PERSIST knob -> step-walk draw spread, exponential: 0 -> 1.0
                // (near-white steps, first-order feel), 0.3 -> 0.25 (default),
                // 1 -> 0.01 (long steady glides). Needed once per cycle, so
                // the pow() lives here, off the per-sample path.
                const float persistSpread =
                    std::pow(10.f, -2.f * params[PERSIST_PARAM].getValue());
                runCycleUpdate(N, scaleAmp, scaleDur, bAmp, bDurMin, bDurMax,
                               dist, persistSpread);
                updateNormAndFreq(lockPitch, args.sampleRate, centerFreq);
                // Continuity at the cycle boundary is inherent: every segment
                // starts from current_amp (the value just reached), so the
                // wrap segment interpolates from the old amp[N-1] to the
                // freshly walked amp[0] with no jump (Serra eq. 1 reads the
                // two as the same point). All N segments carry the walk.
                cycleTrigger.trigger(1e-3f);
            }

            target_amp  = amp[current_breakpoint];
            // Error-diffused rounding of the (possibly normalized) duration:
            // plain truncation runs ~0.5 sample short per segment, which adds
            // up to a few percent of sharpness; carrying the remainder keeps
            // the long-run period exact (which LOCK mode depends on).
            const float fd = dur[current_breakpoint] * norm_k + dur_err;
            current_dur    = std::max(1, (int)(fd + 0.5f));
            dur_err        = clamp(fd - (float)current_dur, -4.f, 4.f);
        }

        // ── Refresh display snapshot (~45 Hz) ─────────────────────────────────
        if (++dispClock >= (int)(args.sampleRate / 45.f)) {       // publish snapshot ~45 Hz
            dispClock = 0;
            int next = 1 - dispBuf.load(std::memory_order_relaxed);
            std::copy(amp, amp + N, dispAmp[next]);
            std::copy(dur, dur + N, dispDur[next]);
            dispN[next]    = N;
            dispBAmp[next] = bAmp;
            dispBuf.store(next, std::memory_order_release);
        }
    }
};


// ─── Polygon display ───────────────────────────────────────────────────────
// The signature dynamic-stochastic-synthesis picture: the N breakpoints drawn
// as a piecewise-linear waveform (x = cumulative duration over one cycle,
// y = amplitude). As the two random walks evolve, the vertices drift — vertically
// (amplitude walk) and horizontally (duration walk) — so the whole polygon slowly
// morphs. The faint band marks the ±B AMP amplitude barriers the walk reflects in.
// Read lock-free from the audio thread's snapshot — fine for a display.
struct GENDYScope : Widget {
    GENDYN* module = nullptr;
    std::shared_ptr<Font> font;

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

            const float W = box.size.x, Hh = box.size.y;
            const float mid = Hh * 0.54f;        // waveform centreline (room for title)
            const float halfH = Hh * 0.40f;      // amp = ±1 maps to ±halfH
            auto Y = [&](float a) { return mid - clamp(a, -1.1f, 1.1f) * halfH; };

            int b = module ? module->dispBuf.load(std::memory_order_acquire) : 0;

            // Amplitude-barrier band (±B AMP): the reflecting walls of the amp walk.
            float bAmp = module ? module->dispBAmp[b] : 0.8f;
            for (float sgn : {-1.f, 1.f}) {
                nvgBeginPath(args.vg);
                nvgMoveTo(args.vg, 0, Y(sgn * bAmp));
                nvgLineTo(args.vg, W, Y(sgn * bAmp));
                nvgStrokeColor(args.vg, nvgRGBA(0x40, 0x90, 0x70, 0x40));
                nvgStrokeWidth(args.vg, 1.f);
                nvgStroke(args.vg);
            }

            // Build the polygon from the snapshot (or a demo shape in the browser).
            int N = (module && module->dispN[b] >= 2) ? module->dispN[b] : 13;
            float total = 0.f;
            for (int i = 0; i < N; i++)
                total += module ? std::max(1.f, module->dispDur[b][i]) : 1.f;
            if (total <= 0.f) total = 1.f;

            auto ampAt = [&](int i) {
                if (module) return module->dispAmp[b][i];
                return 0.7f * std::sin(2.f * (float)M_PI * 2.f * i / N);   // demo
            };
            auto durAt = [&](int i) {
                return module ? std::max(1.f, module->dispDur[b][i]) : 1.f;
            };

            // Vertices at cumulative segment starts; close the period at x = W.
            nvgBeginPath(args.vg);
            float cum = 0.f;
            for (int i = 0; i < N; i++) {
                float x = cum / total * W;
                if (i == 0) nvgMoveTo(args.vg, x, Y(ampAt(0)));
                else        nvgLineTo(args.vg, x, Y(ampAt(i)));
                cum += durAt(i);
            }
            nvgLineTo(args.vg, W, Y(ampAt(0)));     // wrap: period end == start value
            nvgLineCap(args.vg, NVG_ROUND);
            nvgLineJoin(args.vg, NVG_ROUND);
            nvgStrokeColor(args.vg, nvgRGBA(0x55, 0xe0, 0xa0, 0x55));   // green glow
            nvgStrokeWidth(args.vg, 3.2f);
            nvgStroke(args.vg);
            nvgStrokeColor(args.vg, nvgRGB(0xc0, 0xff, 0xd8));          // bright core
            nvgStrokeWidth(args.vg, 1.2f);
            nvgStroke(args.vg);

            // Breakpoint vertices (the walking points), if sparse enough to read.
            if (N <= 40) {
                cum = 0.f;
                for (int i = 0; i < N; i++) {
                    float x = cum / total * W;
                    nvgBeginPath(args.vg);
                    nvgCircle(args.vg, x, Y(ampAt(i)), 1.6f);
                    nvgFillColor(args.vg, nvgRGBA(0xc0, 0xff, 0xd8, 0xdd));
                    nvgFill(args.vg);
                    cum += durAt(i);
                }
            }

            if (!font)
                font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
            if (font) {
                nvgFontFaceId(args.vg, font->handle);
                nvgFontSize(args.vg, mm2px(3.2f));
                nvgTextLetterSpacing(args.vg, mm2px(0.4f));
                nvgFillColor(args.vg, nvgRGBA(0x9a, 0xff, 0xb8, 0xcc));
                nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
                nvgText(args.vg, mm2px(2.6f), mm2px(2.2f), "GENDYN", NULL);
                nvgTextLetterSpacing(args.vg, 0.f);

                char buf[16];
                std::snprintf(buf, sizeof(buf), "N=%d", N);
                nvgFontSize(args.vg, mm2px(2.3f));
                nvgFillColor(args.vg, nvgRGBA(0x60, 0xb0, 0x88, 0xaa));
                nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_BOTTOM);
                nvgText(args.vg, W - mm2px(2.4f), Hh - mm2px(1.8f), buf, NULL);
            }

            nvgResetScissor(args.vg);
        }
        Widget::drawLayer(args, layer);
    }
};


// ─── Widget ───────────────────────────────────────────────────────────────────

struct GENDYWidget : ModuleWidget {

    void draw(const DrawArgs& args) override {
        ModuleWidget::draw(args);

        std::shared_ptr<Font> font = APP->window->loadFont(
            asset::system("res/fonts/DejaVuSans.ttf"));
        if (!font) return;

        nvgSave(args.vg);
        nvgFontFaceId(args.vg, font->handle);
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);

        const NVGcolor dim    = nvgRGB(0x77, 0x77, 0x99);
        const NVGcolor outclr = nvgRGB(0xcc, 0xcc, 0xee);

        auto lbl = [&](float x, float y, float sz, NVGcolor col, const char* s) {
            nvgFontSize(args.vg, mm2px(sz));
            nvgFillColor(args.vg, col);
            nvgText(args.vg, mm2px(x), mm2px(y), s, nullptr);
        };

        // Top control row labels (knobs at y=54)
        lbl( 8.0f, 61.f, 1.8f, dim, "N");
        lbl(20.5f, 61.f, 1.8f, dim, "FREQ");
        lbl(30.5f, 61.f, 1.6f, dim, "LOCK");
        lbl(41.0f, 61.f, 1.8f, dim, "DIST");
        lbl(53.0f, 61.f, 1.7f, dim, "PERSIST");

        // CV channel-strip labels (above each knob at y=74)
        lbl( 9.00f, 68.f, 1.8f, dim, "S AMP");
        lbl(24.32f, 68.f, 1.8f, dim, "S DUR");
        lbl(39.64f, 68.f, 1.8f, dim, "B AMP");
        lbl(54.96f, 68.f, 1.8f, dim, "B WID");

        // Output labels (jacks at y=112)
        lbl(15.0f, 118.5f, 1.9f, outclr, "OUT");
        lbl(30.5f, 118.5f, 1.9f, outclr, "TRIG");
        lbl(46.0f, 118.5f, 1.9f, outclr, "FREQ");

        nvgRestore(args.vg);
    }

    GENDYWidget(GENDYN* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/GENDYN.svg")));

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   122.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 122.0f))));

        // Morphing-polygon scope across the top.
        GENDYScope* scope = new GENDYScope();
        scope->module = module;
        scope->box.pos  = mm2px(Vec(5.5f, 8.f));
        scope->box.size = mm2px(Vec(50.f, 36.f));
        addChild(scope);

        // Top control row (y=54): N | FREQ | LOCK | DIST | PERSIST
        addParam(createParamCentered<RoundBlackSnapKnob>(mm2px(Vec( 8.0f, 54.f)), module, GENDYN::N_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(    mm2px(Vec(20.5f, 54.f)), module, GENDYN::B_DUR_CENTER_PARAM));
        addParam(createParamCentered<CKSS>(              mm2px(Vec(30.5f, 54.f)), module, GENDYN::LOCK_PARAM));
        addParam(createParamCentered<RoundBlackSnapKnob>(mm2px(Vec(41.0f, 54.f)), module, GENDYN::DISTRIBUTION_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(    mm2px(Vec(53.0f, 54.f)), module, GENDYN::PERSIST_PARAM));

        // CV channel strips (knob y=74 / attenuverter y=86 / jack y=96)
        struct Strip { float x; int knob, att, in; };
        const Strip strips[] = {
            { 9.00f, GENDYN::SCALE_AMP_PARAM,  GENDYN::SCALE_AMP_ATT_PARAM, GENDYN::SCALE_AMP_INPUT},
            {24.32f, GENDYN::SCALE_DUR_PARAM,  GENDYN::SCALE_DUR_ATT_PARAM, GENDYN::SCALE_DUR_INPUT},
            {39.64f, GENDYN::B_AMP_PARAM,      GENDYN::B_AMP_ATT_PARAM,     GENDYN::B_AMP_INPUT},
            {54.96f, GENDYN::B_DUR_WIDTH_PARAM,GENDYN::B_DUR_ATT_PARAM,     GENDYN::B_DUR_INPUT},
        };
        for (const Strip& s : strips) {
            addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(s.x, 74.f)), module, s.knob));
            addParam(createParamCentered<Trimpot>(       mm2px(Vec(s.x, 86.f)), module, s.att));
            addInput(createInputCentered<PJ301MPort>(    mm2px(Vec(s.x, 96.f)), module, s.in));
        }

        // Outputs (y=112): OUT | TRIG | FREQ
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(15.0f, 112.f)), module, GENDYN::AUDIO_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(30.5f, 112.f)), module, GENDYN::CYCLE_TRIG_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(46.0f, 112.f)), module, GENDYN::FREQ_CV_OUTPUT));
    }
};

Model* modelGENDYN = createModel<GENDYN, GENDYWidget>("GENDYN");
