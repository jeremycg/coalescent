#include "plugin.hpp"
#include "dsp/rk4.hpp"
#include "tanh_approx.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>

// BUNNIES — a predator–prey oscillator with two modes. PREY (bunnies) and PRED
// (foxes) chase each other by ~a quarter cycle: prey booms, predator booms a beat
// later, prey crashes, predator follows.
//
//   LV — Lotka–Volterra (conservative):
//     dx/dt = x·(1 − y)              dy/dt = gamma·y·(x − 1)
//   RM — Rosenzweig–MacArthur (self-correcting limit cycle):
//     g(x) = x/(1+b·x)
//     dx/dt = x·(1 − x/K) − y·g(x)   dy/dt = s·y·(g(x) − c)
//
// LV is *conservative* — neutral loops, no attractor — so a raw integration
// drifts. A closed-form conserved quantity V = gamma·(x−ln x) + (y−ln y) is used
// as both drift-killer and amplitude knob (WILD): each sample nudges the state
// down the gradient of (V−V0)², scaled by simulation time. RM self-corrects and
// needs no servo. Both integrate on the shared coalescent::rk4<2>.
//
// Servo constants (STAB_K, STAB_FLOOR, and the LV V0 range) were set by offline
// simulation: the plan's originals collapsed large orbits into the positivity
// floor; STAB_K=0.5 / STAB_FLOOR=0.2 / V0≤Vmin+3.5 is the stable, monotonic envelope.

static const int ORBIT_N = 512;   // phase-orbit trail length

struct Bunnies : Module {

    enum ParamId {
        RATE_PARAM,
        BALANCE_PARAM,
        WILD_PARAM,
        MODE_PARAM,
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
        PREY_OUTPUT,
        PRED_OUTPUT,
        PREY_POP_OUTPUT,
        PRED_POP_OUTPUT,
        OUTPUTS_LEN
    };
    enum LightId { LIGHTS_LEN };

    // ─── State (transient; not saved) ───────────────────────────────────────
    float x = 1.3f, y = 0.9f;        // seed OFF the fixed point (LV needs an orbit)
    float cx = 1.f, cy = 1.f;        // active center (fixed point)
    float rmCx = 1.f, rmCy = 1.f;    // cached RM fixed point (separate from active)
    float rmA = -1.f, rmB = -1.f, rmC = -1.f;
    bool  rmValid = false;
    int   lastMode = -1;
    float prevPrey = 0.f, prevPred = 0.f;
    bool  risingPrey = false, risingPred = false;
    dsp::PulseGenerator popGen[2];

    // ─── Visualizer: live ring + lock-free double-buffered snapshot ──────────
    struct Orbit { float pt[ORBIT_N][2] = {}; int head = 0; float peak = 1e-3f; };
    Orbit liveOrbit;                 // audio-thread only
    Orbit snap[2];
    std::atomic<int> snapIndex{0};
    dsp::ClockDivider trailDiv;
    float lastDisplayFs = 0.f;

    // ─── Tunable constants (final values at M8) ─────────────────────────────
    static constexpr float RATE_CAL  = 7.49f;  // measured LV default period × √gamma (tools/stability/bunnies.cpp) → C4
    static constexpr float HSUB_MAX  = 0.05f;
    static constexpr int   MIN_SUB   = 2;   // profiled: adaptive wants K=1 at default; 2 keeps 0-cent + floor margin, halves RK4
    static constexpr int   MAX_SUB   = 64;
    static constexpr float PITCH_TOTAL_MIN = -8.f, PITCH_TOTAL_MAX = 8.f;
    static constexpr float OUT_GAIN  = 0.9f;
    static constexpr float STATE_MAX = 1e3f;
    static constexpr float POS_FLOOR = 1e-4f;
    static constexpr float STAB_K       = 0.5f;   // servo gain (sim-time units) — retuned
    static constexpr float STAB_FLOOR   = 0.2f;   // reciprocal floor in the V-gradient — retuned
    static constexpr float MAX_STAB_STEP = 0.25f;
    static constexpr float LV_V0_RANGE  = 3.5f;   // WILD → V0 = Vmin + [0, 3.5]; keeps max WILD off the positivity floor (retuned from plan 8)
    static constexpr float KICK_GAIN = 0.5f;
    static constexpr float BALANCE_CV_DEPTH = 0.1f, WILD_CV_DEPTH = 0.1f;
    static constexpr float POP_MIN = 0.05f;
    static constexpr float GATE_LEVEL = 10.f, GATE_TIME = 1e-3f;
    static constexpr float RM_B = 0.5f, RM_S = 1.0f;
    static constexpr bool  LV_COMPENSATE_PITCH = true;

    enum { LV = 0, RM = 1 };

    Bunnies() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);
        configParam(RATE_PARAM, -8.f, 4.f, 0.f, "Rate", " Hz", 2.f, dsp::FREQ_C4);
        configParam(BALANCE_PARAM, 0.f, 1.f, 0.5f, "Balance (predator/prey asymmetry)");
        configParam(WILD_PARAM, 0.f, 1.f, 0.4f, "Wild (boom-bust size)");
        configSwitch(MODE_PARAM, 0.f, 1.f, 0.f, "Mode",
                     {"LV — Lotka–Volterra (conservative)", "RM — Rosenzweig–MacArthur (limit cycle)"});
        configParam(BALANCE_ATT_PARAM, -1.f, 1.f, 0.f, "Balance CV");
        configParam(WILD_ATT_PARAM, -1.f, 1.f, 0.f, "Wild CV");

        configInput(VOCT_INPUT, "1V/oct rate");
        configInput(BALANCE_INPUT, "Balance CV");
        configInput(WILD_INPUT, "Wild CV");
        configInput(KICK_INPUT, "Kick / prey-force");

        configOutput(PREY_OUTPUT, "Prey (bunnies)");
        configOutput(PRED_OUTPUT, "Predator (foxes)");
        configOutput(PREY_POP_OUTPUT, "Prey peak");
        configOutput(PRED_POP_OUTPUT, "Predator peak");
    }

    void reseed() { x = cx * 1.1f + 0.05f; y = cy * 0.9f + 0.02f; }
    // Clears peak-detector edge state. Deliberately does NOT reset popGen: an
    // in-flight ~1 ms pulse is allowed to finish across a mode change / reseed
    // (harmless, and avoids chopping a gate). onReset() clears popGen fully.
    void resetPeakMemory() { prevPrey = prevPred = 0.f; risingPrey = risingPred = false; }

    void onReset() override {
        x = 1.3f; y = 0.9f; cx = cy = 1.f; rmCx = rmCy = 1.f;
        rmValid = false; rmA = rmB = rmC = -1.f; lastMode = -1;
        resetPeakMemory();
        for (int i = 0; i < 2; ++i) popGen[i].reset();
        liveOrbit = Orbit{};
    }

    // rising→falling local max above POP_MIN on the centered signal.
    void firePeak(float c, float& prev, bool& rising, dsp::PulseGenerator& gen) {
        if (c > prev + 1e-7f) rising = true;
        else if (c < prev - 1e-7f) {
            if (rising && prev > POP_MIN) gen.trigger(GATE_TIME);
            rising = false;
        }
        prev = c;
    }

    void process(const ProcessArgs& args) override {
        const float fs = args.sampleRate;

        // ── 3.1 read + map params ──
        float balance = clamp(params[BALANCE_PARAM].getValue()
            + inputs[BALANCE_INPUT].getVoltage() * params[BALANCE_ATT_PARAM].getValue() * BALANCE_CV_DEPTH, 0.f, 1.f);
        float wild = clamp(params[WILD_PARAM].getValue()
            + inputs[WILD_INPUT].getVoltage() * params[WILD_ATT_PARAM].getValue() * WILD_CV_DEPTH, 0.f, 1.f);
        float kick = inputs[KICK_INPUT].isConnected() ? inputs[KICK_INPUT].getVoltage() * KICK_GAIN : 0.f;
        int   mode = params[MODE_PARAM].getValue() > 0.5f ? RM : LV;

        float gamma = 0.f, V0 = 0.f, K = 0.f, c = 0.f;
        if (mode == LV) {
            gamma = rack::math::rescale(balance, 0.f, 1.f, 0.2f, 5.f);
            V0 = (gamma + 1.f) + rack::math::rescale(wild, 0.f, 1.f, 0.f, LV_V0_RANGE);
        } else {
            c = clamp(rack::math::rescale(balance, 0.f, 1.f, 0.15f, 0.6f), 0.f, 0.95f / RM_B);  // keep b·c<1
            K = rack::math::rescale(wild, 0.f, 1.f, 1.2f, 12.f);
        }

        // ── 3.2 center (analytic) ──
        if (mode == LV) { cx = 1.f; cy = 1.f; }
        else {
            if (!rmValid || K != rmA || RM_B != rmB || c != rmC) {
                rmCx = c / (1.f - RM_B * c);
                rmCy = rmCx * (1.f - rmCx / K) / c;
                rmA = K; rmB = RM_B; rmC = c; rmValid = true;
            }
            cx = rmCx; cy = rmCy;
            if (!(cx > 0.f && cy > 0.f) || !std::isfinite(cx) || !std::isfinite(cy)) {
                cx = 1.f; cy = 1.f; reseed(); resetPeakMemory(); rmValid = false;
            }
        }

        // ── 3.3 MODE change → recenter + reseed ──
        if (mode != lastMode) { reseed(); resetPeakMemory(); lastMode = mode; }

        // ── 3.4 pitch (+ LV √gamma compensation), adaptive substepping ──
        float pitchTotal = clamp(params[RATE_PARAM].getValue() + inputs[VOCT_INPUT].getVoltage(),
                                 PITCH_TOTAL_MIN, PITCH_TOTAL_MAX);
        float pitchHz = dsp::FREQ_C4 * dsp::approxExp2_taylor5(pitchTotal);
        float dtau = RATE_CAL * pitchHz / fs;
        if (mode == LV && LV_COMPENSATE_PITCH) dtau /= std::sqrt(gamma);  // BALANCE = timbre, not pitch
        dtau = std::min(dtau, HSUB_MAX * MAX_SUB);
        int   Ksub = clamp((int) std::ceil(dtau / HSUB_MAX), MIN_SUB, MAX_SUB);
        float h = dtau / Ksub;

        auto deriv = [&](const float* v, float* dv) {
            float X = std::max(v[0], POS_FLOOR), Y = std::max(v[1], POS_FLOOR);
            if (mode == LV) { dv[0] = X * (1.f - Y);              dv[1] = gamma * Y * (X - 1.f); }
            else            { float g = X / (1.f + RM_B * X);
                              dv[0] = X * (1.f - X / K) - Y * g;  dv[1] = RM_S * Y * (g - c); }
            dv[0] += kick;                                        // continuous prey-force
        };

        float st[2] = {x, y};
        for (int s = 0; s < Ksub; ++s) {
            coalescent::rk4<2>(st, h, deriv);
            st[0] = clamp(st[0], POS_FLOOR, STATE_MAX);
            st[1] = clamp(st[1], POS_FLOOR, STATE_MAX);
            if (!std::isfinite(st[0]) || !std::isfinite(st[1])) { reseed(); st[0] = x; st[1] = y; resetPeakMemory(); break; }
        }

        // ── LV conserved-quantity servo — kills drift AND sets amplitude ──
        if (mode == LV) {
            float X = st[0], Y = st[1];
            float V   = gamma * (X - std::log(X)) + (Y - std::log(Y));
            float dVx = gamma * (1.f - 1.f / std::max(X, STAB_FLOOR));
            float dVy =         (1.f - 1.f / std::max(Y, STAB_FLOOR));
            float stab = STAB_K * dtau;
            float sx = clamp(-stab * (V - V0) * dVx, -MAX_STAB_STEP, MAX_STAB_STEP);
            float sy = clamp(-stab * (V - V0) * dVy, -MAX_STAB_STEP, MAX_STAB_STEP);
            st[0] = std::max(st[0] + sx, POS_FLOOR);
            st[1] = std::max(st[1] + sy, POS_FLOOR);
        }
        x = st[0]; y = st[1];

        // ── 3.5 outputs + peak gates ──
        float cX = x - cx, cY = y - cy;
        outputs[PREY_OUTPUT].setVoltage(5.f * coalescent::fastTanh(cX * OUT_GAIN));
        outputs[PRED_OUTPUT].setVoltage(5.f * coalescent::fastTanh(cY * OUT_GAIN));
        if (std::fabs(cX) < POP_MIN && std::fabs(cY) < POP_MIN) { risingPrey = risingPred = false; }
        firePeak(cX, prevPrey, risingPrey, popGen[0]);
        firePeak(cY, prevPred, risingPred, popGen[1]);
        outputs[PREY_POP_OUTPUT].setVoltage(popGen[0].process(args.sampleTime) ? GATE_LEVEL : 0.f);
        outputs[PRED_POP_OUTPUT].setVoltage(popGen[1].process(args.sampleTime) ? GATE_LEVEL : 0.f);

        // ── visualizer: append to live ring, then publish a full copy ──
        if (fs != lastDisplayFs) {
            lastDisplayFs = fs;
            // ~90 Hz: above the ~60 fps display (trail stays live) but far below the
            // old 300 Hz — less ring-copy on the audio thread and a proportionally
            // smaller double-buffer tear window.
            trailDiv.setDivision(std::max(1, (int) std::round(fs / 90.f)));
        }
        if (trailDiv.process()) {
            liveOrbit.pt[liveOrbit.head][0] = cX;
            liveOrbit.pt[liveOrbit.head][1] = cY;
            liveOrbit.head = (liveOrbit.head + 1) % ORBIT_N;
            float a = std::max(std::fabs(cX), std::fabs(cY));
            liveOrbit.peak = std::max(a, liveOrbit.peak * 0.999f + 1e-5f);
            int back = 1 - snapIndex.load(std::memory_order_relaxed);
            snap[back] = liveOrbit;                              // full ring copy
            snapIndex.store(back, std::memory_order_release);
        }
    }
};

// ── Panel layout (mm) ────────────────────────────────────────────────────────
namespace bpl {
    static const float KNOB_Y = 56.f, ATT_Y = 68.f, IN_Y = 84.f, OUT_Y = 100.f, POP_Y = 114.f;
    static const float KNOB_X[3] = {11.f, 25.5f, 40.f};   // RATE, BALANCE, WILD
    static const float MODE_X = 52.f;
    static const float IN_X[4] = {11.f, 25.5f, 40.f, 52.f};
    static const float OUT_X[2] = {18.f, 43.f};
}

// Phase-orbit display with a bunny riding the current point. Reads only the
// published snapshot; never touches the integrator.
struct OrbitView : widget::TransparentWidget {
    Bunnies* module = nullptr;
    std::shared_ptr<Font> font;

    // Dark "screen" + bezel (house style — the SVG is just panel + rails + screws).
    void draw(const DrawArgs& args) override {
        nvgBeginPath(args.vg);
        nvgRoundedRect(args.vg, 0, 0, box.size.x, box.size.y, mm2px(2.f));
        nvgFillColor(args.vg, nvgRGB(0x12, 0x0a, 0x0e));
        nvgFill(args.vg);
        nvgStrokeColor(args.vg, nvgRGB(0x4d, 0x2b, 0x38));
        nvgStrokeWidth(args.vg, 1.2f);
        nvgStroke(args.vg);
        Widget::draw(args);
    }

    void drawBunny(NVGcontext* vg, float px, float py, float r, NVGcolor col) {
        nvgFillColor(vg, col);
        nvgBeginPath(vg); nvgEllipse(vg, px, py, r * 0.9f, r * 0.7f); nvgFill(vg);          // body
        nvgBeginPath(vg); nvgEllipse(vg, px - r * 0.5f, py - r * 1.1f, r * 0.28f, r * 0.7f); nvgFill(vg);  // ear
        nvgBeginPath(vg); nvgEllipse(vg, px + r * 0.5f, py - r * 1.1f, r * 0.28f, r * 0.7f); nvgFill(vg);  // ear
    }

    void drawFox(NVGcontext* vg, float px, float py, float r, NVGcolor col) {
        nvgFillColor(vg, col);
        nvgBeginPath(vg);   // face: triangle, snout pointing down
        nvgMoveTo(vg, px - r * 0.85f, py - r * 0.5f); nvgLineTo(vg, px + r * 0.85f, py - r * 0.5f);
        nvgLineTo(vg, px, py + r * 0.95f); nvgClosePath(vg); nvgFill(vg);
        nvgBeginPath(vg);   // left ear
        nvgMoveTo(vg, px - r * 0.85f, py - r * 0.5f); nvgLineTo(vg, px - r * 0.45f, py - r * 1.55f);
        nvgLineTo(vg, px - r * 0.05f, py - r * 0.5f); nvgClosePath(vg); nvgFill(vg);
        nvgBeginPath(vg);   // right ear
        nvgMoveTo(vg, px + r * 0.85f, py - r * 0.5f); nvgLineTo(vg, px + r * 0.45f, py - r * 1.55f);
        nvgLineTo(vg, px + r * 0.05f, py - r * 0.5f); nvgClosePath(vg); nvgFill(vg);
    }

    void drawLayer(const DrawArgs& args, int layer) override {
        if (layer != 1) return;
        Vec ctr = box.size.div(2);
        float R = std::min(box.size.x, box.size.y) * 0.42f;

        int read = module ? module->snapIndex.load(std::memory_order_acquire) : 0;
        const Bunnies::Orbit dummy;
        const Bunnies::Orbit& o = module ? module->snap[read] : dummy;
        float peak = std::max(o.peak, 1e-3f);
        auto P = [&](float cx, float cy) {
            return ctr.plus(Vec(clamp(cx / peak, -1.1f, 1.1f) * R, -clamp(cy / peak, -1.1f, 1.1f) * R));
        };

        // origin marker (the fixed point)
        nvgBeginPath(args.vg); nvgCircle(args.vg, ctr.x, ctr.y, 1.5f);
        nvgFillColor(args.vg, nvgRGBA(0xff, 0x8a, 0x9a, 0x60)); nvgFill(args.vg);

        // fading trail oldest→newest
        int head = o.head;
        nvgLineCap(args.vg, NVG_ROUND); nvgLineJoin(args.vg, NVG_ROUND);
        for (int i = 1; i < ORBIT_N; ++i) {
            int a = (head + i - 1) % ORBIT_N, b = (head + i) % ORBIT_N;
            Vec pa = P(o.pt[a][0], o.pt[a][1]), pb = P(o.pt[b][0], o.pt[b][1]);
            float al = (float) i / ORBIT_N;
            nvgBeginPath(args.vg); nvgMoveTo(args.vg, pa.x, pa.y); nvgLineTo(args.vg, pb.x, pb.y);
            nvgStrokeColor(args.vg, nvgRGBA(0xff, 0x6b, 0x8a, (int)(0xd0 * al)));
            nvgStrokeWidth(args.vg, 1.4f); nvgStroke(args.vg);
        }
        // Slowed cursor: sweep the *recorded* orbit at a fixed ~0.2 Hz so the motion
        // is watchable at any pitch (an illustrative slow proxy — an audio-rate point
        // can't be shown at 60 fps; this rides the real loop, just slowly).
        float sweep = (float) std::fmod(system::getTime() * 0.2, 1.0);
        int cur = (head + (int) (sweep * ORBIT_N)) % ORBIT_N;
        Vec live = P(o.pt[cur][0], o.pt[cur][1]);
        float bunnyY = box.size.y - mm2px(3.2f);   // prey creature rides the bottom axis
        float foxX   = mm2px(3.2f);                 // predator creature rides the left axis

        // faint drop-lines: the live point is where prey (bunny, x) meets predator (fox, y)
        nvgStrokeColor(args.vg, nvgRGBA(0xff, 0x8a, 0x9a, 0x55));
        nvgStrokeWidth(args.vg, 1.f);
        nvgBeginPath(args.vg);
        nvgMoveTo(args.vg, live.x, bunnyY); nvgLineTo(args.vg, live.x, live.y);
        nvgMoveTo(args.vg, foxX, live.y);   nvgLineTo(args.vg, live.x, live.y);
        nvgStroke(args.vg);
        nvgBeginPath(args.vg); nvgCircle(args.vg, live.x, live.y, 2.f);
        nvgFillColor(args.vg, nvgRGBA(0xff, 0xe0, 0xe6, 0xdd)); nvgFill(args.vg);

        drawBunny(args.vg, live.x, bunnyY, 4.f, nvgRGB(0xff, 0xc0, 0xd0));   // prey (pink), slides ↔
        drawFox(args.vg, foxX, live.y, 4.f, nvgRGB(0xff, 0x9a, 0x5a));       // predator (orange), slides ↕

        // Title (top-left) + axis hint (bottom-right), DejaVuSans, coral accent.
        if (!font) font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
        if (font) {
            nvgFontFaceId(args.vg, font->handle);
            nvgFontSize(args.vg, mm2px(3.2f));
            nvgTextLetterSpacing(args.vg, mm2px(0.5f));
            nvgFillColor(args.vg, nvgRGBA(0xff, 0xb0, 0xc0, 0xcc));
            nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
            nvgText(args.vg, mm2px(2.6f), mm2px(2.2f), "BUNNIES", NULL);
            nvgTextLetterSpacing(args.vg, 0.f);
            nvgFontSize(args.vg, mm2px(2.2f));
            nvgFillColor(args.vg, nvgRGBA(0xd0, 0x80, 0x90, 0xaa));
            nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_BOTTOM);
            nvgText(args.vg, box.size.x - mm2px(2.4f), box.size.y - mm2px(1.8f), "prey–pred", NULL);
        }
    }
};

struct BunniesWidget : ModuleWidget {
    std::shared_ptr<Font> font;

    BunniesWidget(Bunnies* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Bunnies.svg")));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f, 122.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 122.0f))));

        OrbitView* ov = new OrbitView();
        ov->module = module;
        ov->box.pos = mm2px(Vec(5.f, 8.f));
        ov->box.size = mm2px(Vec(50.96f, 38.f));
        addChild(ov);

        using namespace bpl;
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(KNOB_X[0], KNOB_Y)), module, Bunnies::RATE_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(KNOB_X[1], KNOB_Y)), module, Bunnies::BALANCE_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(KNOB_X[2], KNOB_Y)), module, Bunnies::WILD_PARAM));
        addParam(createParamCentered<CKSS>(mm2px(Vec(MODE_X, KNOB_Y)), module, Bunnies::MODE_PARAM));
        addParam(createParamCentered<Trimpot>(mm2px(Vec(KNOB_X[1], ATT_Y)), module, Bunnies::BALANCE_ATT_PARAM));
        addParam(createParamCentered<Trimpot>(mm2px(Vec(KNOB_X[2], ATT_Y)), module, Bunnies::WILD_ATT_PARAM));

        const int inId[4] = {Bunnies::VOCT_INPUT, Bunnies::BALANCE_INPUT, Bunnies::WILD_INPUT, Bunnies::KICK_INPUT};
        for (int i = 0; i < 4; ++i)
            addInput(createInputCentered<PJ301MPort>(mm2px(Vec(IN_X[i], IN_Y)), module, inId[i]));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(OUT_X[0], OUT_Y)), module, Bunnies::PREY_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(OUT_X[1], OUT_Y)), module, Bunnies::PRED_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(OUT_X[0], POP_Y)), module, Bunnies::PREY_POP_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(OUT_X[1], POP_Y)), module, Bunnies::PRED_POP_OUTPUT));
    }

    void draw(const DrawArgs& args) override {
        ModuleWidget::draw(args);
        if (!font) font = APP->window->loadFont(asset::system("res/fonts/Nunito-Bold.ttf"));
        if (!font) return;
        nvgFontFaceId(args.vg, font->handle);
        nvgFillColor(args.vg, nvgRGB(0xe6, 0xe6, 0xf2));   // near-white, ~13:1 (house legibility spec)
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
        auto label = [&](float xmm, float ymm, const char* s, float sz = 1.9f) {
            nvgFontSize(args.vg, mm2px(sz * 1.72f));       // Nunito Bold, sized up for legibility
            nvgText(args.vg, mm2px(xmm), mm2px(ymm), s, nullptr);
        };
        using namespace bpl;
        const char* kl[3] = {"RATE", "BALANCE", "WILD"};
        for (int i = 0; i < 3; ++i) label(KNOB_X[i], KNOB_Y - 6.f, kl[i]);
        label(MODE_X, KNOB_Y - 6.f, "LV/RM", 1.5f);
        const char* il[4] = {"V/OCT", "BAL", "WILD", "KICK"};
        for (int i = 0; i < 4; ++i) label(IN_X[i], IN_Y - 5.5f, il[i], 1.55f);
        label(OUT_X[0], OUT_Y - 5.5f, "PREY", 1.6f);  label(OUT_X[1], OUT_Y - 5.5f, "PRED", 1.6f);
        label(OUT_X[0], POP_Y - 5.5f, "PREY↑", 1.5f); label(OUT_X[1], POP_Y - 5.5f, "PRED↑", 1.5f);
    }
};

Model* modelBunnies = createModel<Bunnies, BunniesWidget>("Bunnies");
