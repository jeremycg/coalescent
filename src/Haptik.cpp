#include "plugin.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>

// Array ceiling for the mass-spring ring. N (playable) tops out at 128; the
// arrays are sized larger so raising the N ceiling later needs no realloc.
static const int MAX_N = 256;

struct Haptik : Module {

    enum ParamId {
        N_PARAM,
        PITCH_PARAM,
        RATE_PARAM,
        COUPLE_PARAM,
        DAMP_PARAM,
        INJECT_PARAM,
        EXCITE_PARAM,
        FREEZE_PARAM,
        RATE_ATT_PARAM,
        COUPLE_ATT_PARAM,
        DAMP_ATT_PARAM,
        INJECT_ATT_PARAM,
        DRIVER_PARAM,        // appended last: keeps existing positional param IDs stable
        MODE_PARAM,          // Fast (audio-rate) / Slow (haptic divider)
        PARAMS_LEN
    };

    enum InputId {
        VOCT_INPUT,
        RATE_INPUT,
        COUPLE_INPUT,
        DAMP_INPUT,
        INJECT_INPUT,
        TRIG_INPUT,
        EXT_INPUT,
        INPUTS_LEN
    };

    enum OutputId {
        OUT_OUTPUT,
        MOTION_OUTPUT,
        OUTPUTS_LEN
    };

    enum LightId {
        LIGHTS_LEN
    };

    // ─── Persistent state ───────────────────────────────────────────────────
    float y[MAX_N] = {};            // displacement
    float v[MAX_N] = {};            // velocity
    float yPrev[MAX_N] = {};        // slow-mode: lattice state at the previous step (for interp)
    float scanPhase = 0.f;          // [0,1) scan pointer
    int   last_N    = -1;           // -1 forces reinit on first process()
    int   driverIdx = 0;            // mass that excitation / EXT IN drives
    int   divCounter = 0;           // slow-mode: samples since the last lattice step
    dsp::SchmittTrigger trig;
    dsp::TRCFilter<float> dcBlock;  // fixed internal DC blocker on OUT
    float lastFs = 0.f;             // detects SR change to refresh dcBlock cutoff

    // Display snapshot: ~45 Hz the audio thread publishes the lattice into a
    // lock-free double buffer (fill the back buffer, flip dispBuf with a release
    // store); the UI reads the front buffer after an acquire load instead of the
    // live y[] (which churns every sample). Lock-free and race-free.
    float dispY[2][MAX_N] = {};
    int   dispN[2] = {64, 64}, dispDriver[2] = {0, 0};
    std::atomic<int> dispBuf{0};
    int   dispClock = 0;

    // ─── Tunable constants ──────────────────────────────────────────────────
    static constexpr float DAMP_MAX_HZ = 800.f;
    static constexpr float OUT_GAIN    = 1.0f;
    static constexpr float MOTION_GAIN = 4.f;
    static constexpr float EXT_GAIN    = 0.005f;  // ±5V VCO at 0.2 saturated the ring; 0.005 lets it resonate
    static constexpr float BUMP_FRAC   = 0.125f;  // bump half-width = N * frac
    static constexpr float CV_DEPTH        = 0.1f;   // ±5V CV → ±0.5 at full attenuverter
    static constexpr float DRIVE_KEEPALIVE = 0.05f;  // continuous-drive noise amount (EXCITE=3)
    static constexpr float STATE_MAX       = 16.f;   // safety clamp on y[]/v[] (forced/lossless guard)
    static constexpr int   SLOW_DIV        = 256;    // slow-mode: step the lattice every N samples
    static constexpr float KCTR_MAX        = 0.35f;  // clamp centering so per-update ω_max < 2 (slow mode)

    Haptik() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);

        configParam(N_PARAM, 8.f, 128.f, 64.f, "Masses")->snapEnabled = true;

        configParam(PITCH_PARAM, -4.f, 4.f, 0.f, "Pitch", " Hz", 2.f, dsp::FREQ_C4);
        configParam(RATE_PARAM, std::log2(0.05f), std::log2(30.f), std::log2(3.f),
                    "Evolution rate", " Hz", 2.f);
        configParam(COUPLE_PARAM, 0.f, 0.9f, 0.3f, "Coupling");
        configParam(DAMP_PARAM, 0.f, 1.f, 0.35f, "Damping");
        configParam(INJECT_PARAM, 0.f, 1.f, 0.6f, "Inject");

        configParam(EXCITE_PARAM, 0.f, 3.f, 1.f, "Excitation")->snapEnabled = true;

        configSwitch(FREEZE_PARAM, 0.f, 1.f, 0.f, "Freeze", {"Run", "Freeze"});
        configSwitch(MODE_PARAM, 0.f, 1.f, 0.f, "Mode",
                     {"Fast (audio-rate resonator)", "Slow (haptic morph)"});

        configParam(RATE_ATT_PARAM, -1.f, 1.f, 0.f, "Rate CV");
        configParam(COUPLE_ATT_PARAM, -1.f, 1.f, 0.f, "Couple CV");
        configParam(DAMP_ATT_PARAM, -1.f, 1.f, 0.f, "Damp CV");
        configParam(INJECT_ATT_PARAM, -1.f, 1.f, 0.f, "Inject CV");

        // Driver position as a fraction of the ring; 0.25 == the old fixed N/4.
        configParam(DRIVER_PARAM, 0.f, 1.f, 0.25f, "Driver position", "%", 0.f, 100.f);

        configInput(VOCT_INPUT, "1V/oct pitch");
        configInput(RATE_INPUT, "Rate CV");
        configInput(COUPLE_INPUT, "Couple CV");
        configInput(DAMP_INPUT, "Damp CV");
        configInput(INJECT_INPUT, "Inject CV");
        configInput(TRIG_INPUT, "Trigger (re-excite)");
        configInput(EXT_INPUT, "External force");

        configOutput(OUT_OUTPUT, "Audio");
        configOutput(MOTION_OUTPUT, "Motion (mass 0 displacement, audio-rate CV)");
    }

    void onReset() override {
        last_N = -1;        // force a clean lattice reseed (y/v/scanPhase) on next process()
        trig.reset();       // clear trigger edge state
        dcBlock.reset();    // clear DC-blocker filter state
    }

    // Excitation writes displacement y[] (a "pluck"). amt is the inject amount.
    void applyExcite(int shape, float amt, int N) {
        float added = 0.f;   // total displacement injected, to remove its DC below
        switch (shape) {
            case 0:  // impulse: kick a single mass
                y[driverIdx] += amt;
                added = amt;
                break;
            case 1: {  // bump: Hann window, half-width BUMP_FRAC*N, wraps mod N
                int hw = std::max(1, (int)(BUMP_FRAC * N));
                for (int d = -hw; d <= hw; d++) {
                    float w = 0.5f * (1.f + std::cos((float)M_PI * (float)d / (float)hw));
                    int idx = ((driverIdx + d) % N + N) % N;
                    y[idx] += amt * w;
                    added += amt * w;
                }
                break;
            }
            case 2:  // noise: perturb every mass
                for (int i = 0; i < N; i++) {
                    float d = amt * (2.f * random::uniform() - 1.f);
                    y[i] += d;
                    added += d;
                }
                break;
            case 3:  // continuous drive: no impulse here (handled in kernel)
            default:
                break;
        }
        // Remove the DC the excitation injected, so a pluck adds shape but no net
        // offset — otherwise each pluck steps the scanned mean and thumps through the
        // output DC-blocker (audible as a click).
        if (added != 0.f) {
            float m = added / N;
            for (int i = 0; i < N; i++)
                y[i] -= m;
        }
    }

    void process(const ProcessArgs& args) override {
        const float fs = args.sampleRate;

        int N = (int) std::round(params[N_PARAM].getValue());
        N = clamp(N, 8, MAX_N);

        // Driver position is a live param (fraction of the ring), so recompute every
        // sample; set before reinit so the seed bump lands at the chosen mass.
        driverIdx = clamp((int)(params[DRIVER_PARAM].getValue() * N), 0, N - 1);

        // ── reinit on N change (also the first call) ──
        if (N != last_N) {
            std::fill(y, y + MAX_N, 0.f);
            std::fill(v, v + MAX_N, 0.f);
            std::fill(yPrev, yPrev + MAX_N, 0.f);
            scanPhase = 0.f;
            divCounter = 0;
            applyExcite(1, 1.f, N);            // always a full bump (ignores EXCITE/INJECT) so it sounds on load
            last_N = N;
        }

        // ── macros → physics ──
        // RATE_PARAM stores log2(Hz); CV adds in the log domain, att·CV_DEPTH scaled.
        float rateLog = params[RATE_PARAM].getValue()
                      + inputs[RATE_INPUT].getVoltage() * params[RATE_ATT_PARAM].getValue() * CV_DEPTH;
        bool  slow = params[MODE_PARAM].getValue() > 0.5f;
        int   D    = slow ? SLOW_DIV : 1;   // lattice steps every D samples

        float fEvo = clamp(std::exp2(rateLog), 0.05f, 30.f);
        // Effective timestep is D samples. In Fast mode (D=1) this is the original
        // behaviour: kCtr stays tiny so RATE is subtle. In Slow mode the D factor
        // makes RATE/centering meaningful; the clamp keeps per-update ω_max < 2.
        float wc   = 2.f * (float)M_PI * fEvo * D / fs;
        float kCtr = std::min(wc * wc, KCTR_MAX);

        // Coupling clamped hard at 0.9 keeps the *homogeneous* system stable
        // (ω_max = sqrt(kCtr + 4·kSpr) < 2 for symplectic Euler) — do not raise.
        // Forced/lossless boundedness is handled separately by the STATE_MAX clamp.
        float kSpr = clamp(params[COUPLE_PARAM].getValue()
                         + inputs[COUPLE_INPUT].getVoltage() * params[COUPLE_ATT_PARAM].getValue() * CV_DEPTH,
                         0.f, 0.9f);

        float damp = clamp(params[DAMP_PARAM].getValue()
                         + inputs[DAMP_INPUT].getVoltage() * params[DAMP_ATT_PARAM].getValue() * CV_DEPTH,
                         0.f, 1.f);
        float gamma = std::exp(-damp * DAMP_MAX_HZ * D / fs);   // velocity multiplier ≤ 1; D-scaled so
                                                                // wall-clock decay is divider-independent

        float amt = clamp(params[INJECT_PARAM].getValue()
                        + inputs[INJECT_INPUT].getVoltage() * params[INJECT_ATT_PARAM].getValue() * CV_DEPTH,
                        0.f, 1.f);

        bool freeze = params[FREEZE_PARAM].getValue() > 0.5f;
        int  shape  = (int) std::round(params[EXCITE_PARAM].getValue());

        // ── excitation ──
        // Hysteresis window (0.1..1V) so an offset/DC-coupled trigger source can't latch.
        // Gate on !freeze: a pluck writes y[] directly, so allowing it while frozen would
        // change the "held" waveform. trig.process() still runs to keep the edge in sync.
        if (trig.process(inputs[TRIG_INPUT].getVoltage(), 0.1f, 1.f) && !freeze)
            applyExcite(shape, amt, N);

        float drive = inputs[EXT_INPUT].isConnected()
                    ? inputs[EXT_INPUT].getVoltage() * EXT_GAIN * amt : 0.f;
        if (shape == 3 && !freeze)
            drive += amt * DRIVE_KEEPALIVE * (2.f * random::uniform() - 1.f);   // continuous keep-alive

        // ── dynamics: symplectic Euler, two passes (skip if frozen) ──
        // Pass 1 reads only y[] (one snapshot) to form all accelerations;
        // pass 2 reads only v[]. The ordering is required for correctness.
        // ── dynamics: symplectic Euler, two passes ──
        // Fast mode steps every sample; Slow mode steps every D samples (divCounter==0)
        // and keeps yPrev so the readout can interpolate between frames (no stepping).
        bool stepNow = !freeze && (!slow || divCounter == 0);
        if (stepNow) {
            if (slow) std::copy(y, y + N, yPrev);   // snapshot pre-step for inter-frame lerp

            // Velocity pass. Interior masses need no index wrap; the two ring-seam
            // masses are handled out of the hot loop to keep it modulo-free.
            // +1e-20f flushes denormals on long DAMP=0 tails (inaudible, DC-blocked).
            float a0 = kSpr * (y[N - 1] - 2.f * y[0] + y[1]) - kCtr * y[0];
            v[0] = (v[0] + a0) * gamma + 1e-20f;
            for (int i = 1; i < N - 1; i++) {
                float a = kSpr * (y[i - 1] - 2.f * y[i] + y[i + 1]) - kCtr * y[i];
                v[i] = (v[i] + a) * gamma + 1e-20f;
            }
            float aN = kSpr * (y[N - 2] - 2.f * y[N - 1] + y[0]) - kCtr * y[N - 1];
            v[N - 1] = (v[N - 1] + aN) * gamma + 1e-20f;

            // Driver force folded in here: (v+a)·gamma + drive·gamma == (v+a+drive)·gamma.
            v[driverIdx] += drive * gamma;

            // Bound the state itself. The ω_max<2 proof only covers the homogeneous
            // system; external forcing (TRIG, EXT IN) with DAMP=0 (lossless) can pump
            // energy without bound. This generous clamp is never reached in normal play
            // (|y|~1-2) and degrades runaway to bounded saturation instead of NaN.
            // Both must be clamped: clamping y alone lets v keep integrating and snap.
            for (int i = 0; i < N; i++) {
                y[i] = clamp(y[i] + v[i], -STATE_MAX, STATE_MAX);
                v[i] = clamp(v[i], -STATE_MAX, STATE_MAX);
            }
        }

        // ── scan readout (always) ──
        float pitchHz = dsp::FREQ_C4 * std::exp2(
                            params[PITCH_PARAM].getValue() + inputs[VOCT_INPUT].getVoltage());
        scanPhase += pitchHz / fs;
        scanPhase -= std::floor(scanPhase);
        float p   = scanPhase * N;
        int   i0  = std::min((int) p, N - 1);    // guard against p == N at phase rounding
        int   i1  = (i0 + 1) % N;
        float f   = p - i0;
        float s;
        if (slow && !freeze) {
            // Interpolate the two readout points between the previous and current
            // lattice frame (fr = samples since last step / D) so the held shape
            // ramps smoothly instead of stepping every D samples.
            float fr = (float) divCounter / (float) D;
            float a0 = yPrev[i0] + fr * (y[i0] - yPrev[i0]);
            float a1 = yPrev[i1] + fr * (y[i1] - yPrev[i1]);
            s = a0 + f * (a1 - a0);
        } else {
            s = y[i0] + f * (y[i1] - y[i0]);
        }
        if (slow && !freeze)                     // advance AFTER the readout (frame continuity)
            divCounter = (divCounter + 1) % D;

        // ── outputs ──
        if (fs != lastFs) {                      // recompute only on SR change (still
            dcBlock.setCutoffFreq(20.f / fs);    // graceful without onSampleRateChange)
            lastFs = fs;
        }
        dcBlock.process(s);                      // ~20 Hz high-pass; scanned mean wanders
        outputs[OUT_OUTPUT].setVoltage(5.f * std::tanh(dcBlock.highpass() * OUT_GAIN));
        // MOTION taps mass 0 directly (audio-rate in this design); not high-passed.
        outputs[MOTION_OUTPUT].setVoltage(clamp(y[0] * MOTION_GAIN, -5.f, 5.f));

        // ── refresh display snapshot (~45 Hz) ──
        if (++dispClock >= (int)(fs / 45.f)) {
            dispClock = 0;
            int next = 1 - dispBuf.load(std::memory_order_relaxed);
            std::copy(y, y + N, dispY[next]);
            dispN[next] = N;
            dispDriver[next] = driverIdx;
            dispBuf.store(next, std::memory_order_release);
        }
    }
};


// ─── Ring visualizer ─────────────────────────────────────────────────────────
// A scope/radar-style screen: the N masses are vertices on a circle, each pushed
// radially by its displacement y[i] and joined into a glowing closed loop (the
// live waveform). Concentric guides + radial spokes give a graticule; the DRIVER
// mass is an orange node; a single illustrative dot sweeps the ring as the scan
// read-head. State is read lock-free from the audio thread — fine for a display.
struct RingDisplay : Widget {
    Haptik* module = nullptr;
    float dispPeak = 0.6f;   // smoothed peak |y| for auto-scaling the radius
    double lastT  = 0.0;     // last frame time, for time-based smoothing
    std::shared_ptr<Font> font;   // cached once (loaded lazily in drawLayer)

    // Decorative shape shown in the module browser (no running module).
    static float demoShape(int i, int N) {
        return 0.6f * std::sin(2.f * (float)M_PI * 3.f * i / N);
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
            const float cx = box.size.x * 0.5f, cy = box.size.y * 0.5f;
            const float halfmin = std::min(cx, cy);
            const float R = halfmin * 0.66f, targetAmp = halfmin * 0.31f;
            const float rmin = halfmin * 0.10f, rmax = halfmin * 0.97f;

            int b = module ? module->dispBuf.load(std::memory_order_acquire) : 0;
            int N = (module && module->dispN[b] >= 2) ? module->dispN[b] : 64;
            bool freeze = module && module->params[Haptik::FREEZE_PARAM].getValue() > 0.5f;

            // Auto-scale: normalise the ring to a smoothed peak displacement so it stays
            // a readable size instead of zooming with absolute level. Fast attack lets a
            // new pluck fit quickly; slow (~1.2 s) release damps the shrink as it decays.
            float peak = 1e-4f;
            for (int i = 0; i < N; i++) {
                float yi = module ? module->dispY[b][i] : demoShape(i, N);
                peak = std::max(peak, std::fabs(yi));
            }
            double t = system::getTime();
            float dt = (lastT > 0.0) ? (float)(t - lastT) : 0.f;
            lastT = t;
            float tau = (peak > dispPeak) ? 0.10f : 1.2f;
            dispPeak += (peak - dispPeak) * (dt > 0.f ? 1.f - std::exp(-dt / tau) : 0.f);
            dispPeak = std::max(dispPeak, 0.05f);   // floor: silence stays small, not amplified

            auto radiusAt = [&](int i) {
                float yi = module ? module->dispY[b][i] : demoShape(i, N);
                float norm = clamp(yi / dispPeak, -1.4f, 1.4f);
                return clamp(R + norm * targetAmp, rmin, rmax);
            };
            auto ang = [&](float pos) {
                return 2.f * (float)M_PI * pos / N - (float)M_PI / 2.f;  // i=0 at top
            };

            // Graticule: concentric guide rings.
            for (float fr : {0.4f, 0.7f, 1.0f}) {
                nvgBeginPath(args.vg);
                nvgCircle(args.vg, cx, cy, R * fr);
                nvgStrokeColor(args.vg, nvgRGBA(0x40, 0x70, 0x90, fr == 1.0f ? 0x66 : 0x2c));
                nvgStrokeWidth(args.vg, 1.f);
                nvgStroke(args.vg);
            }
            // Graticule: radial spokes.
            for (int k = 0; k < 12; k++) {
                float th = 2.f * (float)M_PI * k / 12.f;
                nvgBeginPath(args.vg);
                nvgMoveTo(args.vg, cx, cy);
                nvgLineTo(args.vg, cx + R * std::cos(th), cy + R * std::sin(th));
                nvgStrokeColor(args.vg, nvgRGBA(0x30, 0x50, 0x70, 0x22));
                nvgStrokeWidth(args.vg, 0.75f);
                nvgStroke(args.vg);
            }

            // Live ring shape: translucent fill, soft glow stroke, bright core stroke.
            nvgBeginPath(args.vg);
            for (int i = 0; i < N; i++) {
                float th = ang((float)i), r = radiusAt(i);
                float px = cx + r * std::cos(th), py = cy + r * std::sin(th);
                if (i == 0) nvgMoveTo(args.vg, px, py); else nvgLineTo(args.vg, px, py);
            }
            nvgClosePath(args.vg);
            nvgFillColor(args.vg, nvgRGBA(0x40, 0xc0, 0xff, 0x22));
            nvgFill(args.vg);
            nvgStrokeColor(args.vg, nvgRGBA(0x55, 0xc8, 0xff, 0x40));
            nvgStrokeWidth(args.vg, 3.5f);
            nvgStroke(args.vg);
            nvgStrokeColor(args.vg, nvgRGB(0x9a, 0xe4, 0xff));
            nvgStrokeWidth(args.vg, 1.3f);
            nvgStroke(args.vg);

            // Mass nodes (only when sparse enough to read).
            if (N <= 48)
                for (int i = 0; i < N; i++) {
                    float th = ang((float)i), r = radiusAt(i);
                    nvgBeginPath(args.vg);
                    nvgCircle(args.vg, cx + r * std::cos(th), cy + r * std::sin(th), 1.3f);
                    nvgFillColor(args.vg, nvgRGBA(0xc0, 0xee, 0xff, 0xdd));
                    nvgFill(args.vg);
                }

            // Driver mass node (orange: glow + core).
            int drv = clamp(module ? module->dispDriver[b] : N / 4, 0, N - 1);
            {
                float th = ang((float)drv), r = radiusAt(drv);
                float px = cx + r * std::cos(th), py = cy + r * std::sin(th);
                nvgBeginPath(args.vg); nvgCircle(args.vg, px, py, 4.f);
                nvgFillColor(args.vg, nvgRGBA(0xff, 0x9b, 0x3a, 0x55)); nvgFill(args.vg);
                nvgBeginPath(args.vg); nvgCircle(args.vg, px, py, 2.f);
                nvgFillColor(args.vg, nvgRGB(0xff, 0xb0, 0x55)); nvgFill(args.vg);
            }

            // Scan read-head: a single illustrative dot sweeping ~0.5 Hz (NOT the true
            // read position — that runs at audio rate and can't be drawn at frame rate,
            // so no trail is shown). It keeps moving under FREEZE, demonstrating that
            // scanning is independent of the (now static) ring shape.
            float sp = (float) std::fmod(system::getTime() * 0.5, 1.0);
            float pos = sp * N;
            int i0 = std::min((int)pos, N - 1);
            float frac = pos - i0;
            float r = radiusAt(i0) + frac * (radiusAt((i0 + 1) % N) - radiusAt(i0));
            float th = ang(pos);
            float sx = cx + r * std::cos(th), sy = cy + r * std::sin(th);
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, cx, cy); nvgLineTo(args.vg, sx, sy);
            nvgStrokeColor(args.vg, nvgRGBA(0xff, 0xee, 0x88, 0x40));
            nvgStrokeWidth(args.vg, 1.f); nvgStroke(args.vg);
            nvgBeginPath(args.vg); nvgCircle(args.vg, sx, sy, 4.f);
            nvgFillColor(args.vg, nvgRGBA(0xff, 0xee, 0x88, 0x55)); nvgFill(args.vg);
            nvgBeginPath(args.vg); nvgCircle(args.vg, sx, sy, 2.2f);
            nvgFillColor(args.vg, nvgRGB(0xff, 0xf2, 0xa0)); nvgFill(args.vg);

            // Centre dot.
            nvgBeginPath(args.vg); nvgCircle(args.vg, cx, cy, 1.3f);
            nvgFillColor(args.vg, nvgRGBA(0x88, 0xcc, 0xff, 0x90)); nvgFill(args.vg);

            // Screen text: title (top-left), RUN/FREEZE status (top-right), N (bottom-left).
            if (!font)
                font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
            if (font) {
                nvgFontFaceId(args.vg, font->handle);
                nvgFontSize(args.vg, mm2px(3.4f));
                nvgTextLetterSpacing(args.vg, mm2px(0.5f));
                nvgFillColor(args.vg, nvgRGBA(0x9a, 0xb0, 0xff, 0xcc));
                nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
                nvgText(args.vg, mm2px(3.f), mm2px(2.6f), "HAPTIK", NULL);
                nvgTextLetterSpacing(args.vg, 0.f);

                nvgFontSize(args.vg, mm2px(2.4f));
                nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_TOP);
                if (freeze) { nvgFillColor(args.vg, nvgRGB(0x6a, 0xb0, 0xff)); }
                else        { nvgFillColor(args.vg, nvgRGBA(0x6a, 0xd0, 0x8a, 0xcc)); }
                nvgText(args.vg, box.size.x - mm2px(3.f), mm2px(2.9f), freeze ? "FREEZE" : "RUN", NULL);

                char buf[16];
                std::snprintf(buf, sizeof(buf), "N=%d", N);
                nvgFontSize(args.vg, mm2px(2.3f));
                nvgFillColor(args.vg, nvgRGBA(0x60, 0x80, 0xb0, 0xaa));
                nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_BOTTOM);
                nvgText(args.vg, mm2px(3.f), box.size.y - mm2px(2.6f), buf, NULL);
            }
        }
        Widget::drawLayer(args, layer);
    }
};


// ─── Widget ─────────────────────────────────────────────────────────────────

struct HaptikWidget : ModuleWidget {
    std::shared_ptr<Font> font;   // cached once (loaded lazily in draw)

    void draw(const DrawArgs& args) override {
        ModuleWidget::draw(args);

        if (!font)
            font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
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

        // Voice row labels (below each knob/switch at y=64)
        lbl(11.f, 70.f, 1.7f, dim, "N");
        lbl(26.f, 70.f, 1.7f, dim, "PITCH");
        lbl(41.f, 70.f, 1.7f, dim, "EXCITE");
        lbl(56.f, 70.f, 1.7f, dim, "DRIVER");
        lbl(70.f, 70.f, 1.7f, dim, "FREEZE");
        lbl(81.f, 70.f, 1.7f, dim, "MODE");

        // CV channel-strip labels (above each knob at y=82)
        lbl(13.f, 76.f, 1.9f, dim, "RATE");
        lbl(32.f, 76.f, 1.9f, dim, "COUPLE");
        lbl(51.f, 76.f, 1.9f, dim, "DAMP");
        lbl(70.f, 76.f, 1.9f, dim, "INJECT");

        // I/O labels (below each jack at y=113)
        lbl(13.f, 119.f, 1.8f, dim, "V/OCT");
        lbl(27.f, 119.f, 1.8f, dim, "TRIG");
        lbl(41.f, 119.f, 1.8f, dim, "EXT");
        lbl(64.f, 119.f, 1.9f, outclr, "OUT");
        lbl(78.f, 119.f, 1.9f, outclr, "MOTION");

        nvgRestore(args.vg);
    }

    HaptikWidget(Haptik* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Haptik.svg")));

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(85.42f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   122.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(85.42f, 122.0f))));

        // Ring visualizer — big scope screen across the top
        RingDisplay* ring = new RingDisplay();
        ring->module = module;
        ring->box.pos  = mm2px(Vec(6.5f, 8.f));
        ring->box.size = mm2px(Vec(78.f, 46.f));
        addChild(ring);

        // Voice row (y=64): N | PITCH | EXCITE | DRIVER | FREEZE | MODE
        addParam(createParamCentered<RoundBlackSnapKnob>(mm2px(Vec(11.f, 64.f)), module, Haptik::N_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(    mm2px(Vec(26.f, 64.f)), module, Haptik::PITCH_PARAM));
        addParam(createParamCentered<RoundBlackSnapKnob>(mm2px(Vec(41.f, 64.f)), module, Haptik::EXCITE_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(    mm2px(Vec(56.f, 64.f)), module, Haptik::DRIVER_PARAM));
        addParam(createParamCentered<CKSS>(              mm2px(Vec(70.f, 64.f)), module, Haptik::FREEZE_PARAM));
        addParam(createParamCentered<CKSS>(              mm2px(Vec(81.f, 64.f)), module, Haptik::MODE_PARAM));

        // CV channel strips (knob y=82 / attenuverter y=93 / jack y=102)
        struct Strip { float x; int knob, att, in; };
        const Strip strips[] = {
            {13.f, Haptik::RATE_PARAM,   Haptik::RATE_ATT_PARAM,   Haptik::RATE_INPUT},
            {32.f, Haptik::COUPLE_PARAM, Haptik::COUPLE_ATT_PARAM, Haptik::COUPLE_INPUT},
            {51.f, Haptik::DAMP_PARAM,   Haptik::DAMP_ATT_PARAM,   Haptik::DAMP_INPUT},
            {70.f, Haptik::INJECT_PARAM, Haptik::INJECT_ATT_PARAM, Haptik::INJECT_INPUT},
        };
        for (const Strip& s : strips) {
            addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(s.x, 82.f)), module, s.knob));
            addParam(createParamCentered<Trimpot>(       mm2px(Vec(s.x, 93.f)), module, s.att));
            addInput(createInputCentered<PJ301MPort>(    mm2px(Vec(s.x, 102.f)), module, s.in));
        }

        // I/O row (y=113): V/OCT | TRIG | EXT  ...  OUT | MOTION
        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec(13.f, 113.f)), module, Haptik::VOCT_INPUT));
        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec(27.f, 113.f)), module, Haptik::TRIG_INPUT));
        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec(41.f, 113.f)), module, Haptik::EXT_INPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(64.f, 113.f)), module, Haptik::OUT_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(78.f, 113.f)), module, Haptik::MOTION_OUTPUT));
    }
};

Model* modelHaptik = createModel<Haptik, HaptikWidget>("Haptik");
