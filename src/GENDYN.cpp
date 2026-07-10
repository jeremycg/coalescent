#include "plugin.hpp"
#include "dsp/display_snapshot.hpp"
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

    // Seed shape (0=Sine 1=Triangle 2=Saw 3=Square 4=Random) + a menu re-seed
    // request. Sine default: the oscillator boots as a clean pitched tone that the
    // stochastic walk then evolves — far more discoverable than a random noise seed.
    int   initShape = 0;
    int   lastInitShape = -1;        // forces a reseed when the shape changes
    bool  reseedPending = false;
    bool  restoredPending = false;   // dataFromJson set the arrays; process() must recompute
                                     // the SR-dependent caches (norm_k / freq_cv / current_dur)
    dsp::TRCFilter<float> dcBlock;   // gentle ~5 Hz DC blocker on OUT (walk mean drifts)
    float lastFs = 0.f;

    // Display snapshot: ~45 Hz the audio thread publishes the breakpoints into a
    // lock-free triple buffer; the UI reads the latest complete frame, so it never
    // races the per-cycle writer. See dsp/display_snapshot.hpp.
    struct DisplayFrame {
        float amp[MAX_N] = {};
        float dur[MAX_N] = {};
        int   n = 13;
        float bAmp = 0.8f;
    };
    coalescent::DisplaySnapshot<DisplayFrame> displaySnapshot;
    int   dispClock = 0;
    int   dispPeriod = 980;         // samples between snapshots (~45 Hz; refreshed on SR change)

    GENDYN() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);

        configParam(N_PARAM,            2.f,   64.f,   13.f,  "Breakpoints (N)");
        getParamQuantity(N_PARAM)->snapEnabled = true;

        // Defaults for the second-order walk: steps persist across cycles, so net
        // drift is ~8x faster than an uncorrelated walk at equal scale; the 0.35
        // scaling holds a musical evolution timescale.
        configParam(SCALE_AMP_PARAM,    0.f,   1.f,    0.002f, "Amplitude step scale");
        configParam(SCALE_DUR_PARAM,    0.f,   1.f,    0.0035f,"Duration step scale");
        configParam(B_AMP_PARAM,        0.f,   1.f,    0.8f,  "Amplitude barrier");
        configParam(B_DUR_CENTER_PARAM, 20.f,  5000.f, 261.626f, "Center frequency", " Hz");
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

    json_t* dataToJson() override {
        json_t* root = json_object();
        json_object_set_new(root, "initShape", json_integer(initShape));
        // Persist the evolved waveform + walk velocities so a saved sound reloads as
        // itself, not a fresh reseed. Only once the arrays are real (last_N > 0).
        if (last_N > 0) {
            int n = last_N;
            json_object_set_new(root, "N", json_integer(n));
            json_t *ja = json_array(), *jd = json_array(), *jsa = json_array(), *jsd = json_array();
            for (int i = 0; i < n; i++) {
                json_array_append_new(ja,  json_real(amp[i]));
                json_array_append_new(jd,  json_real(dur[i]));
                json_array_append_new(jsa, json_real(step_amp[i]));
                json_array_append_new(jsd, json_real(step_dur[i]));
            }
            json_object_set_new(root, "amp", ja);
            json_object_set_new(root, "dur", jd);
            json_object_set_new(root, "stepAmp", jsa);
            json_object_set_new(root, "stepDur", jsd);
        }
        return root;
    }
    void dataFromJson(json_t* root) override {
        if (json_t* j = json_object_get(root, "initShape"))
            initShape = clamp((int) json_integer_value(j), 0, 4);
        // Restore the evolved waveform if present and well-formed (length matches the
        // saved N, durations >= 1, all finite); otherwise fall through to a reseed so
        // pre-serialization patches still load.
        json_t *jn = json_object_get(root, "N"),  *ja = json_object_get(root, "amp");
        json_t *jd = json_object_get(root, "dur"), *jsa = json_object_get(root, "stepAmp");
        json_t *jsd = json_object_get(root, "stepDur");
        if (!(jn && ja && jd && jsa && jsd)) return;
        int n = (int) json_integer_value(jn);
        if (n < 2 || n > MAX_N) return;
        if ((int) json_array_size(ja) != n || (int) json_array_size(jd) != n
            || (int) json_array_size(jsa) != n || (int) json_array_size(jsd) != n) return;
        // Durations are floored at 1 and capped well below INT_MAX so the float→int
        // conversions in playback (dur[0]*norm_k, error-diffused fd) can't overflow
        // on a corrupt/hand-edited patch — that would be undefined behaviour. Any
        // real duration is at most a few thousand samples.
        constexpr float DUR_MAX = 1e6f;
        float A[MAX_N], Dr[MAX_N], SA[MAX_N], SD[MAX_N];
        for (int i = 0; i < n; i++) {
            A[i]  = (float) json_number_value(json_array_get(ja,  i));
            Dr[i] = (float) json_number_value(json_array_get(jd,  i));
            SA[i] = (float) json_number_value(json_array_get(jsa, i));
            SD[i] = (float) json_number_value(json_array_get(jsd, i));
            if (!std::isfinite(A[i]) || !std::isfinite(Dr[i]) || Dr[i] < 1.f || Dr[i] > DUR_MAX
                || !std::isfinite(SA[i]) || !std::isfinite(SD[i])) return;   // malformed → reseed
        }
        float s = 0.f;
        for (int i = 0; i < n; i++) { amp[i] = A[i]; dur[i] = Dr[i]; step_amp[i] = SA[i]; step_dur[i] = SD[i]; s += Dr[i]; }
        sum_dur = s;
        current_breakpoint = 0; current_sample = 0;          // clean cycle start from the restored shape
        current_amp = amp[n - 1]; target_amp = amp[0];
        current_dur = std::max(1, (int) dur[0]); dur_err = 0.f;
        last_N = n; lastInitShape = initShape; reseedPending = false;   // suppress the reinit reseed
        restoredPending = true;   // first process() recomputes norm_k / freq_cv / current_dur (needs sampleRate)
    }

    void onReset() override {
        // Restore factory menu state and force a clean, unconditional re-seed on
        // the next process() (covers the case where N didn't change).
        initShape = 0;
        reseedPending = true;
        restoredPending = false;
        dcBlock.reset();
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

    // Seed all breakpoints and reset the oscillator to breakpoint 0. Deterministic
    // shapes (sine/triangle/saw/square) start from a clean waveform with equal
    // segment durations (durCenter) so the tone is pitched and recognisable; the
    // stochastic walk then evolves it. "Random" reproduces the classic noisy seed.
    void initBreakpoints(int N, float bAmp, float durCenter, float bDurMin, float bDurMax) {
        sum_dur = 0.f;
        for (int i = 0; i < N; i++) {
            const float ph = (float) i / (float) N;   // 0..1 over one cycle
            float a;
            switch (initShape) {
                case 1:  a = 1.f - 4.f * std::fabs(ph - 0.5f); break;    // triangle
                case 2:  a = 2.f * ph - 1.f; break;                      // saw
                case 3:  a = ph < 0.5f ? 1.f : -1.f; break;              // square
                case 4:  a = 2.f * rack::random::uniform() - 1.f; break; // random
                default: a = std::sin(2.f * (float) M_PI * ph); break;   // sine
            }
            amp[i] = a * bAmp;
            dur[i] = (initShape == 4)
                   ? bDurMin + rack::random::uniform() * (bDurMax - bDurMin)  // random pitch jitter
                   : std::max(1.f, durCenter);                                              // equal → stable pitch
            step_amp[i] = 0.f;
            step_dur[i] = 0.f;
            sum_dur += dur[i];
        }
        current_breakpoint = 0;
        current_sample     = 0;
        current_amp        = amp[N - 1];   // wrap value → the first cycle is exactly periodic
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
    void updateNormAndFreq(int N, bool lockPitch, float sampleRate, float centerFreq) {
        norm_k = 1.f;
        if (lockPitch && sum_dur > 0.f) {
            norm_k = (sampleRate / centerFreq) / sum_dur;
            // Playback floors each segment at 1 sample, so a cycle can't be
            // shorter than N samples. Clamp norm_k to keep the FREQ CV honest
            // when the LOCK target (fs/centerFreq) is below that reachable floor.
            norm_k = std::max(norm_k, (float) N / sum_dur);
        }
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
        // fraction of that centre value. bDurWidth=0 collapses the window to the
        // single value durCenter: gainDur→0 pins the duration walk and reflect()
        // returns durCenter, so the (error-diffused) playback holds the exact
        // centre pitch. Do NOT widen a zero-width window — that reintroduces a
        // duration walk and detunes (a +1-sample floor biased the pitch flat).
        const float durCenter = args.sampleRate / (centerFreq * (float)N);
        const float halfWidth = bDurWidth * durCenter;
        float bDurMin = std::max(1.f, durCenter - halfWidth);
        float bDurMax = std::max(bDurMin, durCenter + halfWidth);

        int dist = clamp((int)params[DISTRIBUTION_PARAM].getValue(), 0, 3);

        const bool lockPitch = params[LOCK_PARAM].getValue() > 0.5f;

        // ── Reinit on N change, seed-shape change, or a menu re-seed ───────────
        if (N != last_N || initShape != lastInitShape || reseedPending) {
            initBreakpoints(N, bAmp, durCenter, bDurMin, bDurMax);
            updateNormAndFreq(N, lockPitch, args.sampleRate, centerFreq);
            current_dur = std::max(1, (int)(dur[0] * norm_k));
            last_N = N; lastInitShape = initShape; reseedPending = false;
            restoredPending = false;                 // reinit already recomputed the caches
        } else if (restoredPending) {
            // Just restored from a patch: dataFromJson set the arrays but couldn't
            // compute the sample-rate-dependent caches. Do it now (without reseeding
            // the restored waveform) so the first cycle — pitch under LOCK and the
            // FREQ output — is correct from sample 0, not one cycle late.
            updateNormAndFreq(N, lockPitch, args.sampleRate, centerFreq);
            current_dur = std::max(1, (int)(dur[0] * norm_k));
            restoredPending = false;
        }

        // ── Generate output sample ────────────────────────────────────────────
        const float t      = (float)current_sample / (float)current_dur;
        const float output = current_amp + t * (target_amp - current_amp);
        if (args.sampleRate != lastFs) {          // refresh SR-derived caches on SR change
            dcBlock.setCutoffFreq(5.f / args.sampleRate);
            dispPeriod = (int) (args.sampleRate / 45.f);
            lastFs = args.sampleRate;
        }
        dcBlock.process(output);                  // strip the walk's slow DC drift
        outputs[AUDIO_OUTPUT].setVoltage(dcBlock.highpass() * 5.f);

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
                updateNormAndFreq(N, lockPitch, args.sampleRate, centerFreq);
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
        if (++dispClock >= dispPeriod) {                          // publish snapshot ~45 Hz
            dispClock = 0;
            DisplayFrame& fr = displaySnapshot.writable();
            std::copy(amp, amp + N, fr.amp);
            std::copy(dur, dur + N, fr.dur);
            fr.n    = N;
            fr.bAmp = bAmp;
            displaySnapshot.publish();
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

            static const GENDYN::DisplayFrame dummy;
            const GENDYN::DisplayFrame& fr = module ? module->displaySnapshot.consume() : dummy;

            // Amplitude-barrier band (±B AMP): the reflecting walls of the amp walk.
            float bAmp = module ? fr.bAmp : 0.8f;
            for (float sgn : {-1.f, 1.f}) {
                nvgBeginPath(args.vg);
                nvgMoveTo(args.vg, 0, Y(sgn * bAmp));
                nvgLineTo(args.vg, W, Y(sgn * bAmp));
                nvgStrokeColor(args.vg, nvgRGBA(0x40, 0x90, 0x70, 0x40));
                nvgStrokeWidth(args.vg, 1.f);
                nvgStroke(args.vg);
            }

            // Build the polygon from the snapshot (or a demo shape in the browser).
            int N = (module && fr.n >= 2) ? fr.n : 13;
            float total = 0.f;
            for (int i = 0; i < N; i++)
                total += module ? std::max(1.f, fr.dur[i]) : 1.f;
            if (total <= 0.f) total = 1.f;

            auto ampAt = [&](int i) {
                if (module) return fr.amp[i];
                return 0.7f * std::sin(2.f * (float)M_PI * 2.f * i / N);   // demo
            };
            auto durAt = [&](int i) {
                return module ? std::max(1.f, fr.dur[i]) : 1.f;
            };

            // Match the audio's segment/duration pairing: segment i interpolates
            // amp[i-1] -> amp[i] over dur[i], and a cycle begins from amp[N-1]. So
            // amp[i] sits at the cumulative END of dur[0..i], and the wrap value
            // amp[N-1] anchors both x=0 and x=W. (With equal durations this is the
            // old picture rotated one segment; with unequal durations it now matches
            // the waveform being heard.)
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, 0.f, Y(ampAt(N - 1)));
            float cum = 0.f;
            for (int i = 0; i < N; i++) {
                cum += durAt(i);
                nvgLineTo(args.vg, cum / total * W, Y(ampAt(i)));
            }
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
                    cum += durAt(i);                       // amp[i] at the end of dur[i] (see polygon above)
                    float x = cum / total * W;
                    nvgBeginPath(args.vg);
                    nvgCircle(args.vg, x, Y(ampAt(i)), 1.6f);
                    nvgFillColor(args.vg, nvgRGBA(0xc0, 0xff, 0xd8, 0xdd));
                    nvgFill(args.vg);
                }
            }

            std::shared_ptr<Font> font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
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

struct GENDYNLabels : widget::Widget {

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

        // Top control row labels (knobs at y=54)
        lbl( 8.0f, 46.f, 2.1f, dim, "N");
        lbl(20.5f, 46.f, 2.1f, dim, "FREQ");
        lbl(30.5f, 46.f, 2.1f, dim, "LOCK");
        lbl(41.0f, 46.f, 2.1f, dim, "DIST");
        lbl(53.0f, 46.f, 2.1f, dim, "PERSIST");

        // CV channel-strip labels (above each knob at y=74)
        lbl( 7.5f, 66.f, 2.1f, dim, "S AMP");
        lbl(22.82f, 66.f, 2.1f, dim, "S DUR");
        lbl(38.14f, 66.f, 2.1f, dim, "B AMP");
        lbl(53.46f, 66.f, 2.1f, dim, "B WID");

        // Output labels (jacks at y=112)
        lbl(15.0f, 118.5f, 2.1f, outclr, "OUT");
        lbl(30.5f, 118.5f, 2.1f, outclr, "TRIG");
        lbl(46.0f, 118.5f, 2.1f, outclr, "FREQ");

        nvgRestore(args.vg);
    }
};

struct GENDYWidget : ModuleWidget {

    GENDYWidget(GENDYN* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/GENDYN.svg")));
        addFramebufferedLabels<GENDYNLabels>(this);

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   122.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 122.0f))));

        // Morphing-polygon scope across the top.
        GENDYScope* scope = new GENDYScope();
        scope->module = module;
        scope->box.pos  = mm2px(Vec(5.5f, 8.f));
        scope->box.size = mm2px(Vec(50.f, 34.f));
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
            { 7.5f, GENDYN::SCALE_AMP_PARAM,  GENDYN::SCALE_AMP_ATT_PARAM, GENDYN::SCALE_AMP_INPUT},
            {22.82f, GENDYN::SCALE_DUR_PARAM,  GENDYN::SCALE_DUR_ATT_PARAM, GENDYN::SCALE_DUR_INPUT},
            {38.14f, GENDYN::B_AMP_PARAM,      GENDYN::B_AMP_ATT_PARAM,     GENDYN::B_AMP_INPUT},
            {53.46f, GENDYN::B_DUR_WIDTH_PARAM,GENDYN::B_DUR_ATT_PARAM,     GENDYN::B_DUR_INPUT},
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

    void appendContextMenu(Menu* menu) override {
        GENDYN* m = dynamic_cast<GENDYN*>(module);
        if (!m) return;
        menu->addChild(new MenuSeparator);
        menu->addChild(createIndexPtrSubmenuItem("Initial waveform",
            {"Sine", "Triangle", "Saw", "Square", "Random"}, &m->initShape));
        menu->addChild(createMenuItem("Re-seed waveform", "",
            [m]() { m->reseedPending = true; }));
    }
};

Model* modelGENDYN = createModel<GENDYN, GENDYWidget>("GENDYN");
