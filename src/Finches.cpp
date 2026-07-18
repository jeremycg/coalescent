#include "plugin.hpp"
#include "dsp/display_snapshot.hpp"
#include "dsp/finches_field.hpp"

#include <algorithm>
#include <cmath>

// FINCHES - evolutionary branching in a continuous trait distribution.
//
// The state is 64 nonnegative bin masses over trait x in [-1, 1]. A quadratic
// environmental niche pulls the population toward ENV, mutation diffuses mass
// between neighbouring traits, and a Gaussian competition kernel penalizes
// similar phenotypes. COMPETE controls the dimensionless branching number B:
// below B=1 the environmental optimum is locally stabilizing; above B=1 it is
// disruptive, although finite mutation shifts the visible two-peak transition.
//
// This is a slow, stateful CV source rather than an audio oscillator. The field
// runs at a fixed control rate; pitch, mass, and spread outputs are smoothed at
// audio rate. See dsp/finches_field.hpp for the positive, mass-conserving solver.

struct Finches : Module {
    enum ParamId {
        RATE_PARAM,
        MUTATE_PARAM,
        COMPETE_PARAM,
        NICHE_PARAM,
        MUTANT_PARAM,
        SEED_PARAM,
        MUTATE_ATT_PARAM,
        COMPETE_ATT_PARAM,
        PARAMS_LEN
    };
    enum InputId {
        RATE_INPUT,
        MUTATE_INPUT,
        COMPETE_INPUT,
        ENV_INPUT,
        SEED_INPUT,
        RESET_INPUT,
        INPUTS_LEN
    };
    enum OutputId {
        MASS_L_OUTPUT,
        MASS_R_OUTPUT,
        PITCH_L_OUTPUT,
        PITCH_R_OUTPUT,
        SPREAD_OUTPUT,
        SPLIT_OUTPUT,
        MERGE_OUTPUT,
        OUTPUTS_LEN
    };
    enum LightId { LIGHTS_LEN };

    static const int BINS = coalescent::FinchesField::kBins;

    coalescent::FinchesField field;
    dsp::ClockDivider fieldDiv;
    dsp::ClockDivider publishDiv;
    dsp::SchmittTrigger seedInputTrigger;
    dsp::SchmittTrigger seedButtonTrigger;
    dsp::SchmittTrigger resetTrigger;
    dsp::PulseGenerator splitPulse;
    dsp::PulseGenerator mergePulse;

    bool pendingSeed = false;
    bool pendingReset = false;
    float pendingSeedTrait = 0.f;
    float lastFs = 0.f;
    float fieldWallStep = 1.f / 500.f;
    float outputAlpha = 1.f;
    coalescent::FinchesField::Parameters currentParameters;

    enum SmoothId { MASS_L_SMOOTH, MASS_R_SMOOTH, PITCH_L_SMOOTH,
                    PITCH_R_SMOOTH, SPREAD_SMOOTH, SMOOTH_LEN };
    float target[SMOOTH_LEN] = {};
    float smooth[SMOOTH_LEN] = {};

    struct DisplayFrame {
        float mass[BINS] = {};
        float environment = 0.f;
        float niche = 0.32f;
        float lowTrait = 0.f;
        float highTrait = 0.f;
        bool split = false;

        // Give Rack's static module screenshot (where no engine is running) the
        // same deterministic ancestor users see at startup. Live snapshots
        // replace this frame as soon as process() begins.
        DisplayFrame() {
            float total = 0.f;
            for (int i = 0; i < BINS; ++i) {
                float z = coalescent::FinchesField::traitAt(i) / 0.055f;
                mass[i] = std::exp(-0.5f * z * z);
                total += mass[i];
            }
            for (int i = 0; i < BINS; ++i)
                mass[i] /= total;
        }
    };
    coalescent::DisplaySnapshot<DisplayFrame> displaySnapshot;

    struct SaveFrame {
        coalescent::FinchesField::State state;
        bool valid = false;
    };
    coalescent::DisplaySnapshot<SaveFrame> saveSnapshot;

    static constexpr float FIELD_HZ = 500.f;
    static constexpr float RATE_BASE = 8.f;       // evolutionary tau per wall second at RATE=0
    static constexpr float RATE_TOTAL_MIN = -12.f;
    static constexpr float RATE_TOTAL_MAX = 4.25f;
    static constexpr float PANEL_MUTATION_MAX = 0.00012f;
    static constexpr float CV_DEPTH = 0.1f;       // 10 V at full attenuverter spans a knob
    static constexpr float ENV_DEPTH = 0.11f;     // +/-5 V reaches almost the full ENV range
    static constexpr float MUTANT_RANGE = 0.70f;
    static constexpr float TRAIT_VOLTS = 2.f;     // trait position -> V/oct output
    static constexpr float SPREAD_FULL_SCALE = 0.50f;
    static constexpr float OUTPUT_TAU = 0.02f;
    static constexpr float OUTPUT_SNAP_EPSILON = 1e-9f;
    static constexpr float GATE_LEVEL = 10.f;
    static constexpr float GATE_TIME = 1e-3f;

    Finches() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);
        configParam(RATE_PARAM, -8.f, 4.f, 0.f, "Evolution rate", "x", 2.f, 1.f);
        configParam(MUTATE_PARAM, 0.f, 1.f, 0.35f, "Mutation diffusion", "",
                    PANEL_MUTATION_MAX / coalescent::FinchesField::mutationMin(),
                    coalescent::FinchesField::mutationMin());
        configParam(COMPETE_PARAM, 0.f, 1.f, 0.30f, "Competition locality (branching number)", "",
                    coalescent::FinchesField::branchingMax() / coalescent::FinchesField::branchingMin(),
                    coalescent::FinchesField::branchingMin());
        configParam(NICHE_PARAM, 0.f, 1.f, 0.48f, "Environmental niche width", " trait",
                    0.f, coalescent::FinchesField::nicheMax() - coalescent::FinchesField::nicheMin(),
                    coalescent::FinchesField::nicheMin());
        configParam(MUTANT_PARAM, -1.f, 1.f, 0.50f, "Mutant trait offset", " trait",
                    0.f, MUTANT_RANGE);
        configButton(SEED_PARAM, "Inject mutant");
        configParam(MUTATE_ATT_PARAM, -1.f, 1.f, 0.f, "Mutation CV");
        configParam(COMPETE_ATT_PARAM, -1.f, 1.f, 0.f, "Competition CV");

        configInput(RATE_INPUT, "Evolution rate (1 V doubles speed; combined limit +4.25 oct)");
        configInput(MUTATE_INPUT, "Mutation CV");
        configInput(COMPETE_INPUT, "Competition CV");
        configInput(ENV_INPUT, "Environmental optimum CV");
        configInput(SEED_INPUT, "Inject mutant trigger");
        configInput(RESET_INPUT, "Reset to one adapted cluster");

        configOutput(MASS_L_OUTPUT, "Lower-trait cluster abundance");
        configOutput(MASS_R_OUTPUT, "Higher-trait cluster abundance");
        configOutput(PITCH_L_OUTPUT, "Lower-trait position / unsplit mean (1 V/oct)");
        configOutput(PITCH_R_OUTPUT, "Higher-trait position / unsplit mean (1 V/oct)");
        configOutput(SPREAD_OUTPUT, "Trait spread");
        configOutput(SPLIT_OUTPUT, "Accepted split trigger");
        configOutput(MERGE_OUTPUT, "Accepted merge trigger");

        currentParameters = readParameters();
        setOutputTargets(true);
        publishSaveFrame();
    }

    static float safe(float x, float fallback = 0.f) {
        return std::isfinite(x) ? x : fallback;
    }

    static float clampf(float x, float lo, float hi) {
        return std::max(lo, std::min(x, hi));
    }

    static float exponentialMap(float lo, float hi, float unit) {
        unit = clampf(unit, 0.f, 1.f);
        return lo * std::pow(hi / lo, unit);
    }

    float inputVoltage(int id) {
        return inputs[id].isConnected() ? safe(inputs[id].getVoltage()) : 0.f;
    }

    float environmentFromInput() {
        return clampf(inputVoltage(ENV_INPUT) * ENV_DEPTH,
                      coalescent::FinchesField::environmentMin(),
                      coalescent::FinchesField::environmentMax());
    }

    coalescent::FinchesField::Parameters readParameters() {
        coalescent::FinchesField::Parameters p;
        float mutationUnit = clampf(safe(params[MUTATE_PARAM].getValue())
            + inputVoltage(MUTATE_INPUT) * safe(params[MUTATE_ATT_PARAM].getValue()) * CV_DEPTH,
            0.f, 1.f);
        float branchingUnit = clampf(safe(params[COMPETE_PARAM].getValue())
            + inputVoltage(COMPETE_INPUT) * safe(params[COMPETE_ATT_PARAM].getValue()) * CV_DEPTH,
            0.f, 1.f);
        float nicheUnit = clampf(safe(params[NICHE_PARAM].getValue(), 0.48f), 0.f, 1.f);

        p.mutation = exponentialMap(coalescent::FinchesField::mutationMin(),
                                    PANEL_MUTATION_MAX, mutationUnit);
        p.branching = exponentialMap(coalescent::FinchesField::branchingMin(),
                                     coalescent::FinchesField::branchingMax(), branchingUnit);
        p.niche = coalescent::FinchesField::nicheMin()
            + nicheUnit * (coalescent::FinchesField::nicheMax()
                          - coalescent::FinchesField::nicheMin());
        p.environment = environmentFromInput();
        return p;
    }

    float mutantTrait(float environment) {
        float offset = clampf(safe(params[MUTANT_PARAM].getValue()), -1.f, 1.f) * MUTANT_RANGE;
        return clampf(environment + offset,
                      coalescent::FinchesField::traitMin() + 0.05f,
                      coalescent::FinchesField::traitMax() - 0.05f);
    }

    void setOutputTargets(bool immediate = false) {
        const coalescent::FinchesField::Metrics& m = field.metrics();
        target[MASS_L_SMOOTH] = clampf(m.lowMass, 0.f, 1.f) * 10.f;
        target[MASS_R_SMOOTH] = clampf(m.highMass, 0.f, 1.f) * 10.f;
        target[PITCH_L_SMOOTH] = clampf(m.lowTrait, -1.f, 1.f) * TRAIT_VOLTS;
        target[PITCH_R_SMOOTH] = clampf(m.highTrait, -1.f, 1.f) * TRAIT_VOLTS;
        target[SPREAD_SMOOTH] = clampf(m.spread / SPREAD_FULL_SCALE, 0.f, 1.f) * 10.f;
        if (immediate)
            std::copy(target, target + SMOOTH_LEN, smooth);
    }

    void clearEventState() {
        splitPulse.reset();
        mergePulse.reset();
    }

    void onReset(const ResetEvent& event) override {
        Module::onReset(event);
        currentParameters = readParameters();
        field.reset(currentParameters.environment);
        pendingSeed = pendingReset = false;
        seedInputTrigger.reset();
        seedButtonTrigger.reset();
        resetTrigger.reset();
        clearEventState();
        setOutputTargets(true);
        publishSaveFrame();
        publishDisplayFrame(currentParameters);
    }

    void publishDisplayFrame(const coalescent::FinchesField::Parameters& p) {
        DisplayFrame& d = displaySnapshot.writable();
        field.copyMasses(d.mass, BINS);
        const coalescent::FinchesField::Metrics& m = field.metrics();
        d.environment = p.environment;
        d.niche = p.niche;
        d.lowTrait = m.lowTrait;
        d.highTrait = m.highTrait;
        d.split = m.split;
        displaySnapshot.publish();
    }

    void publishSaveFrame() {
        SaveFrame& s = saveSnapshot.writable();
        s.state = field.state();
        s.valid = true;
        saveSnapshot.publish();
    }

    json_t* dataToJson() override {
        json_t* root = json_object();
        const SaveFrame& s = saveSnapshot.consume();
        if (!s.valid)
            return root;

        json_object_set_new(root, "finchesVersion",
                            json_integer(coalescent::FinchesField::stateVersion()));
        json_t* density = json_array();
        for (int i = 0; i < BINS; ++i)
            json_array_append_new(density, json_real(s.state.mass[i]));
        json_object_set_new(root, "density", density);
        json_object_set_new(root, "split", json_boolean(s.state.split));
        json_object_set_new(root, "splitTimer", json_real(s.state.splitTimer));
        json_object_set_new(root, "mergeTimer", json_real(s.state.mergeTimer));
        return root;
    }

    void dataFromJson(json_t* root) override {
        json_t* density = json_object_get(root, "density");
        if (!density || !json_is_array(density) || json_array_size(density) != BINS)
            return;

        coalescent::FinchesField::State restored;
        for (int i = 0; i < BINS; ++i) {
            json_t* value = json_array_get(density, i);
            if (!value || !json_is_number(value))
                return;
            restored.mass[i] = static_cast<float>(json_number_value(value));
        }

        coalescent::FinchesField candidate;
        json_t* version = json_object_get(root, "finchesVersion");
        if (!version) {
            // Density-only patches predate persistence of the split latch and
            // detector timers. Preserve their historical strong-threshold load.
            if (!candidate.restoreMasses(restored.mass.data(), BINS))
                return;
        }
        else {
            json_t* split = json_object_get(root, "split");
            json_t* splitTimer = json_object_get(root, "splitTimer");
            json_t* mergeTimer = json_object_get(root, "mergeTimer");
            if (!json_is_integer(version)
                || json_integer_value(version) != coalescent::FinchesField::stateVersion()
                || !split || !json_is_boolean(split)
                || !splitTimer || !json_is_number(splitTimer)
                || !mergeTimer || !json_is_number(mergeTimer))
                return;
            restored.version = static_cast<int>(json_integer_value(version));
            restored.split = json_is_true(split);
            restored.splitTimer = static_cast<float>(json_number_value(splitTimer));
            restored.mergeTimer = static_cast<float>(json_number_value(mergeTimer));
            if (!candidate.restore(restored))
                return;
        }

        field = candidate;
        currentParameters = readParameters();
        pendingSeed = pendingReset = false;
        clearEventState();
        setOutputTargets(true);
        publishSaveFrame();
        publishDisplayFrame(currentParameters);
    }

    void process(const ProcessArgs& args) override {
        const float fs = args.sampleRate;
        if (fs != lastFs) {
            lastFs = fs;
            int fieldDivision = std::max(1, static_cast<int>(std::round(fs / FIELD_HZ)));
            fieldDiv.setDivision(fieldDivision);
            fieldDiv.reset();
            fieldWallStep = fieldDivision / fs;
            publishDiv.setDivision(std::max(1, static_cast<int>(std::round(fs / 45.f))));
            publishDiv.reset();
            outputAlpha = 1.f - std::exp(-args.sampleTime / OUTPUT_TAU);
        }

        bool seedEdge = seedInputTrigger.process(inputVoltage(SEED_INPUT), 0.1f, 1.f);
        seedEdge = seedButtonTrigger.process(safe(params[SEED_PARAM].getValue()), 0.1f, 0.9f) || seedEdge;
        if (seedEdge) {
            pendingSeed = true;
            pendingSeedTrait = mutantTrait(environmentFromInput());
        }
        if (resetTrigger.process(inputVoltage(RESET_INPUT), 0.1f, 1.f))
            pendingReset = true;

        if (fieldDiv.process()) {
            currentParameters = readParameters();
            bool resetThisTick = false;
            if (pendingReset) {
                field.reset(currentParameters.environment);
                pendingReset = false;
                pendingSeed = false;
                resetThisTick = true;
                clearEventState();
            }
            if (pendingSeed) {
                field.seed(pendingSeedTrait);
                pendingSeed = false;
            }

            if (!resetThisTick) {
                float rateOct = clampf(
                    safe(params[RATE_PARAM].getValue()) + inputVoltage(RATE_INPUT),
                    RATE_TOTAL_MIN, RATE_TOTAL_MAX);
                float deltaTau = RATE_BASE * std::exp2(rateOct) * fieldWallStep;
                field.advance(deltaTau, currentParameters);
                const coalescent::FinchesField::Metrics& metrics = field.metrics();
                if (metrics.splitEvent) splitPulse.trigger(GATE_TIME);
                if (metrics.mergeEvent) mergePulse.trigger(GATE_TIME);
            }
            setOutputTargets();
            publishSaveFrame();
        }

        for (int i = 0; i < SMOOTH_LEN; ++i) {
            const float delta = target[i] - smooth[i];
            if (std::fabs(delta) <= OUTPUT_SNAP_EPSILON)
                smooth[i] = target[i];
            else
                smooth[i] += delta * outputAlpha;
        }
        outputs[MASS_L_OUTPUT].setVoltage(smooth[MASS_L_SMOOTH]);
        outputs[MASS_R_OUTPUT].setVoltage(smooth[MASS_R_SMOOTH]);
        outputs[PITCH_L_OUTPUT].setVoltage(smooth[PITCH_L_SMOOTH]);
        outputs[PITCH_R_OUTPUT].setVoltage(smooth[PITCH_R_SMOOTH]);
        outputs[SPREAD_OUTPUT].setVoltage(smooth[SPREAD_SMOOTH]);
        outputs[SPLIT_OUTPUT].setVoltage(splitPulse.process(args.sampleTime) ? GATE_LEVEL : 0.f);
        outputs[MERGE_OUTPUT].setVoltage(mergePulse.process(args.sampleTime) ? GATE_LEVEL : 0.f);

        if (publishDiv.process())
            publishDisplayFrame(currentParameters);
    }
};

namespace finches_layout {
    static const float KNOB_Y = 58.f;
    static const float ATT_Y = 70.f;
    static const float IN_Y = 84.f;
    static const float OUT_Y = 100.f;
    static const float EVENT_Y = 114.f;
    static const float COL5[5] = {10.f, 22.8f, 35.6f, 48.4f, 61.2f};
    static const float COL3[3] = {18.f, 35.6f, 53.2f};
}

struct FinchesDensityView : widget::TransparentWidget {
    Finches* module = nullptr;
    static constexpr float VISIBLE_MASS = 0.24f;

    void draw(const DrawArgs& args) override {
        nvgBeginPath(args.vg);
        nvgRoundedRect(args.vg, 0.f, 0.f, box.size.x, box.size.y, mm2px(2.f));
        nvgFillColor(args.vg, nvgRGB(0x07, 0x13, 0x15));
        nvgFill(args.vg);
        nvgStrokeColor(args.vg, nvgRGB(0x2b, 0x4c, 0x48));
        nvgStrokeWidth(args.vg, 1.2f);
        nvgStroke(args.vg);
        Widget::draw(args);
    }

    static float traitX(float trait, float width) {
        return clamp((trait - coalescent::FinchesField::traitMin())
                         / (coalescent::FinchesField::traitMax()
                            - coalescent::FinchesField::traitMin()),
                     0.f, 1.f) * width;
    }

    static float massAt(const Finches::DisplayFrame& frame, float trait) {
        float bin = (trait - coalescent::FinchesField::traitMin())
                  / coalescent::FinchesField::binWidth() - 0.5f;
        int i0 = static_cast<int>(std::floor(bin));
        float frac = bin - i0;
        i0 = clamp(i0, 0, Finches::BINS - 1);
        int i1 = std::min(i0 + 1, Finches::BINS - 1);
        return frame.mass[i0] + (frame.mass[i1] - frame.mass[i0]) * clamp(frac, 0.f, 1.f);
    }

    void drawFinch(NVGcontext* vg, float x, float y, float r, int direction, NVGcolor colour) {
        float d = direction < 0 ? -1.f : 1.f;
        nvgFillColor(vg, colour);
        nvgBeginPath(vg);
        nvgEllipse(vg, x, y, r * 0.90f, r * 0.58f);
        nvgFill(vg);
        nvgBeginPath(vg);
        nvgCircle(vg, x + d * r * 0.62f, y - r * 0.40f, r * 0.43f);
        nvgFill(vg);
        nvgBeginPath(vg);
        nvgMoveTo(vg, x - d * r * 0.72f, y + r * 0.05f);
        nvgLineTo(vg, x - d * r * 1.35f, y - r * 0.18f);
        nvgLineTo(vg, x - d * r * 0.85f, y + r * 0.40f);
        nvgClosePath(vg);
        nvgFill(vg);
        nvgBeginPath(vg);
        nvgMoveTo(vg, x + d * r * 0.98f, y - r * 0.44f);
        nvgLineTo(vg, x + d * r * 1.42f, y - r * 0.28f);
        nvgLineTo(vg, x + d * r * 0.98f, y - r * 0.12f);
        nvgClosePath(vg);
        nvgFill(vg);
        nvgFillColor(vg, nvgRGBA(0x08, 0x14, 0x16, 0x90));
        nvgBeginPath(vg);
        nvgEllipse(vg, x - d * r * 0.12f, y + r * 0.02f, r * 0.48f, r * 0.30f);
        nvgFill(vg);
    }

    void drawLayer(const DrawArgs& args, int layer) override {
        if (layer != 1)
            return;
        static const Finches::DisplayFrame dummy;
        const Finches::DisplayFrame& frame = module ? module->displaySnapshot.consume() : dummy;
        const float W = box.size.x;
        const float H = box.size.y;
        const float top = H * 0.24f;
        const float base = H - mm2px(2.8f);
        const float plotHeight = base - top;
        auto densityY = [&](float mass) {
            return base - clamp(mass / VISIBLE_MASS, 0.f, 1.f) * plotHeight;
        };

        nvgScissor(args.vg, 0.f, 0.f, W, H);
        for (int guide = 1; guide < 4; ++guide) {
            float x = W * guide / 4.f;
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, x, top);
            nvgLineTo(args.vg, x, base);
            nvgStrokeColor(args.vg, nvgRGBA(0x7c, 0xb7, 0xaa, 0x18));
            nvgStrokeWidth(args.vg, 1.f);
            nvgStroke(args.vg);
        }

        nvgBeginPath(args.vg);
        for (int i = 0; i < Finches::BINS; ++i) {
            float trait = coalescent::FinchesField::traitAt(i);
            float z = (trait - frame.environment) / std::max(frame.niche, 1e-3f);
            float x = traitX(trait, W);
            float y = base - std::exp(-0.5f * z * z) * plotHeight * 0.56f;
            if (i == 0) nvgMoveTo(args.vg, x, y);
            else        nvgLineTo(args.vg, x, y);
        }
        nvgStrokeColor(args.vg, nvgRGBA(0xf0, 0xc8, 0x72, 0x45));
        nvgStrokeWidth(args.vg, 1.f);
        nvgStroke(args.vg);

        float envX = traitX(frame.environment, W);
        nvgBeginPath(args.vg);
        nvgMoveTo(args.vg, envX, top);
        nvgLineTo(args.vg, envX, base);
        nvgStrokeColor(args.vg, nvgRGBA(0xf0, 0xc8, 0x72, 0x70));
        nvgStrokeWidth(args.vg, 1.f);
        nvgStroke(args.vg);

        nvgBeginPath(args.vg);
        nvgMoveTo(args.vg, traitX(coalescent::FinchesField::traitAt(0), W), base);
        for (int i = 0; i < Finches::BINS; ++i)
            nvgLineTo(args.vg, traitX(coalescent::FinchesField::traitAt(i), W), densityY(frame.mass[i]));
        nvgLineTo(args.vg, traitX(coalescent::FinchesField::traitAt(Finches::BINS - 1), W), base);
        nvgClosePath(args.vg);
        nvgFillColor(args.vg, nvgRGBA(0x63, 0xc8, 0xb2, 0x35));
        nvgFill(args.vg);

        nvgBeginPath(args.vg);
        for (int i = 0; i < Finches::BINS; ++i) {
            float x = traitX(coalescent::FinchesField::traitAt(i), W);
            float y = densityY(frame.mass[i]);
            if (i == 0) nvgMoveTo(args.vg, x, y);
            else        nvgLineTo(args.vg, x, y);
        }
        nvgStrokeColor(args.vg, nvgRGB(0x83, 0xdd, 0xc8));
        nvgStrokeWidth(args.vg, 1.5f);
        nvgLineJoin(args.vg, NVG_ROUND);
        nvgStroke(args.vg);

        float birdR = mm2px(0.95f);
        auto birdY = [&](float trait) {
            return clamp(densityY(massAt(frame, trait)) - birdR * 0.75f,
                         top + birdR * 0.8f, base - birdR);
        };
        if (frame.split) {
            drawFinch(args.vg, traitX(frame.lowTrait, W), birdY(frame.lowTrait), birdR, -1,
                      nvgRGB(0x74, 0xc9, 0xf0));
            drawFinch(args.vg, traitX(frame.highTrait, W), birdY(frame.highTrait), birdR, 1,
                      nvgRGB(0xf0, 0xb8, 0x68));
        }
        else {
            drawFinch(args.vg, traitX(frame.lowTrait, W), birdY(frame.lowTrait), birdR, 1,
                      nvgRGB(0xe8, 0xd4, 0x82));
        }
        nvgResetScissor(args.vg);

        std::shared_ptr<Font> font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
        if (font) {
            nvgFontFaceId(args.vg, font->handle);
            nvgFontSize(args.vg, mm2px(3.2f));
            nvgTextLetterSpacing(args.vg, mm2px(0.5f));
            nvgFillColor(args.vg, nvgRGBA(0xe8, 0xd4, 0x82, 0xdd));
            nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
            nvgText(args.vg, mm2px(2.6f), mm2px(2.2f), "FINCHES", nullptr);
            nvgTextLetterSpacing(args.vg, 0.f);
            nvgFontSize(args.vg, mm2px(2.1f));
            nvgFillColor(args.vg, nvgRGBA(0x8a, 0xb7, 0xae, 0xa0));
            nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_BOTTOM);
            nvgText(args.vg, W - mm2px(2.4f), H - mm2px(1.8f), "trait ecology", nullptr);
        }
    }
};

struct FinchesLabels : widget::Widget {
    void draw(const DrawArgs& args) override {
        std::shared_ptr<Font> font = APP->window->loadFont(asset::system("res/fonts/Nunito-Bold.ttf"));
        if (!font)
            return;
        nvgFontFaceId(args.vg, font->handle);
        nvgTextLetterSpacing(args.vg, 0.f);
        nvgFillColor(args.vg, nvgRGB(0xe6, 0xf1, 0xef));
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
        auto label = [&](float xmm, float ymm, const char* text, float size = 1.65f) {
            nvgFontSize(args.vg, mm2px(size * 1.72f));
            nvgText(args.vg, mm2px(xmm), mm2px(ymm), text, nullptr);
        };
        using namespace finches_layout;
        const char* knobs[5] = {"RATE", "MUTATE", "COMPETE", "NICHE", "MUTANT"};
        for (int i = 0; i < 5; ++i)
            label(COL5[i], KNOB_Y - 8.f, knobs[i], i == 2 ? 1.50f : 1.62f);
        const char* inputs[5] = {"RATE", "MUTATE", "COMP", "ENV", "SEED"};
        for (int i = 0; i < 5; ++i)
            label(COL5[i], IN_Y - 5.5f, inputs[i], 1.45f);
        const char* outputs[5] = {"MASS L", "MASS R", "SPREAD", "PITCH L", "PITCH R"};
        for (int i = 0; i < 5; ++i)
            label(COL5[i], OUT_Y - 5.5f, outputs[i], 1.35f);
        const char* events[3] = {"RESET", "SPLIT", "MERGE"};
        for (int i = 0; i < 3; ++i)
            label(COL3[i], EVENT_Y - 5.5f, events[i], 1.50f);
    }
};

struct FinchesWidget : ModuleWidget {
    FinchesWidget(Finches* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Finches.svg")));
        addPanelLabels<FinchesLabels>(this);

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.f, 1.f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(65.12f, 1.f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.f, 122.f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(65.12f, 122.f))));

        FinchesDensityView* view = new FinchesDensityView();
        view->module = module;
        view->box.pos = mm2px(Vec(6.f, 8.f));
        view->box.size = mm2px(Vec(59.12f, 38.f));
        addChild(view);

        using namespace finches_layout;
        const int knobs[5] = {Finches::RATE_PARAM, Finches::MUTATE_PARAM, Finches::COMPETE_PARAM,
                              Finches::NICHE_PARAM, Finches::MUTANT_PARAM};
        for (int i = 0; i < 5; ++i)
            addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(COL5[i], KNOB_Y)), module, knobs[i]));
        addParam(createParamCentered<Trimpot>(mm2px(Vec(COL5[1], ATT_Y)), module, Finches::MUTATE_ATT_PARAM));
        addParam(createParamCentered<Trimpot>(mm2px(Vec(COL5[2], ATT_Y)), module, Finches::COMPETE_ATT_PARAM));
        addParam(createParamCentered<TL1105>(mm2px(Vec(COL5[4], ATT_Y)), module, Finches::SEED_PARAM));

        const int inputs[5] = {Finches::RATE_INPUT, Finches::MUTATE_INPUT, Finches::COMPETE_INPUT,
                               Finches::ENV_INPUT, Finches::SEED_INPUT};
        for (int i = 0; i < 5; ++i)
            addInput(createInputCentered<PJ301MPort>(mm2px(Vec(COL5[i], IN_Y)), module, inputs[i]));

        const int outputs[5] = {Finches::MASS_L_OUTPUT, Finches::MASS_R_OUTPUT, Finches::SPREAD_OUTPUT,
                                Finches::PITCH_L_OUTPUT, Finches::PITCH_R_OUTPUT};
        for (int i = 0; i < 5; ++i)
            addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(COL5[i], OUT_Y)), module, outputs[i]));
        addInput(createInputCentered<PJ301MPort>(mm2px(Vec(COL3[0], EVENT_Y)), module, Finches::RESET_INPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(COL3[1], EVENT_Y)), module, Finches::SPLIT_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(COL3[2], EVENT_Y)), module, Finches::MERGE_OUTPUT));
    }
};

Model* modelFinches = createModel<Finches, FinchesWidget>("Finches");
