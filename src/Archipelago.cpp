#include "plugin.hpp"
#include "dsp/archipelago_field.hpp"
#include "dsp/display_snapshot.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>

// ARCHIPELAGO - local adaptation across eight connected habitats.
//
// Each habitat carries a complete abundance density over a continuous trait.
// Selection follows a spatial environmental gradient, mutation diffuses in
// trait space, and migration couples neighbouring habitats. The deterministic
// field lives in dsp/archipelago_field.hpp; this wrapper supplies control-rate
// scheduling, CV presentation, event pulses, persistence, and visualization.

namespace {

float archipelagoClamp(float value, float low, float high) {
    return std::max(low, std::min(value, high));
}

float archipelagoSafe(float value, float fallback = 0.f) {
    return std::isfinite(value) ? value : fallback;
}

struct ArchipelagoSelectionQuantity : ParamQuantity {
    static float fromUnit(float unit) {
        unit = archipelagoClamp(unit, 0.f, 1.f);
        return 8.f * unit * unit;
    }

    float getDisplayValue() override {
        Param* parameter = getParam();
        return fromUnit(parameter ? parameter->getValue() : defaultValue);
    }

    void setDisplayValue(float displayValue) override {
        if (!std::isfinite(displayValue))
            return;
        displayValue = archipelagoClamp(displayValue, 0.f, 8.f);
        setImmediateValue(std::sqrt(displayValue / 8.f));
    }
};

struct ArchipelagoMutationQuantity : ParamQuantity {
    static float fromUnit(float unit) {
        unit = archipelagoClamp(unit, 0.f, 1.f);
        if (!(unit > 0.f))
            return 0.f;
        return 1e-5f * std::pow(300.f, unit);
    }

    float getDisplayValue() override {
        Param* parameter = getParam();
        return fromUnit(parameter ? parameter->getValue() : defaultValue);
    }

    void setDisplayValue(float displayValue) override {
        if (!std::isfinite(displayValue))
            return;
        if (!(displayValue > 0.f)) {
            setImmediateValue(0.f);
            return;
        }
        displayValue = archipelagoClamp(displayValue, 1e-5f, 0.003f);
        setImmediateValue(std::log(displayValue / 1e-5f) / std::log(300.f));
    }
};

struct ArchipelagoMigrationQuantity : ParamQuantity {
    static float fromUnit(float unit) {
        unit = archipelagoClamp(unit, 0.f, 1.f);
        if (!(unit > 0.f))
            return 0.f;
        return 0.002f * std::pow(1000.f, unit);
    }

    float getDisplayValue() override {
        Param* parameter = getParam();
        return fromUnit(parameter ? parameter->getValue() : defaultValue);
    }

    void setDisplayValue(float displayValue) override {
        if (!std::isfinite(displayValue))
            return;
        if (!(displayValue > 0.f)) {
            setImmediateValue(0.f);
            return;
        }
        displayValue = archipelagoClamp(displayValue, 0.002f, 2.f);
        setImmediateValue(std::log(displayValue / 0.002f) / std::log(1000.f));
    }
};

} // namespace

struct Archipelago : Module {
    enum ParamId {
        RATE_PARAM,
        SELECT_PARAM,
        MUTATE_PARAM,
        MIGRATE_PARAM,
        GRADIENT_PARAM,
        BARRIER_PARAM,
        CLIMATE_PARAM,
        RESET_PARAM,
        SELECT_ATT_PARAM,
        MUTATE_ATT_PARAM,
        MIGRATE_ATT_PARAM,
        GRADIENT_ATT_PARAM,
        BARRIER_ATT_PARAM,
        TOPOLOGY_PARAM,
        PARAMS_LEN
    };
    enum InputId {
        RATE_INPUT,
        SELECT_INPUT,
        MUTATE_INPUT,
        MIGRATE_INPUT,
        GRADIENT_INPUT,
        BARRIER_INPUT,
        CLIMATE_INPUT,
        RESET_INPUT,
        INPUTS_LEN
    };
    enum OutputId {
        TRAIT_OUTPUT,
        MASS_OUTPUT,
        MEAN_OUTPUT,
        DIFF_OUTPUT,
        FLUX_OUTPUT,
        COLONIZE_OUTPUT,
        EXTINCT_OUTPUT,
        OUTPUTS_LEN
    };
    enum LightId { LIGHTS_LEN };

    static const int HABITATS = coalescent::ArchipelagoField::kHabitats;
    static const int BINS = coalescent::ArchipelagoField::kBins;
    static const int VALUES = coalescent::ArchipelagoField::kValues;
    static const int TRAIT_SMOOTH = 0;
    static const int MASS_SMOOTH = HABITATS;
    static const int MEAN_SMOOTH = 2 * HABITATS;
    static const int DIFF_SMOOTH = MEAN_SMOOTH + 1;
    static const int FLUX_SMOOTH = DIFF_SMOOTH + 1;
    static const int SMOOTH_LEN = FLUX_SMOOTH + 1;

    static constexpr float FIELD_HZ = 500.f;
    static constexpr float RATE_BASE = 4.f;
    static constexpr float RATE_TOTAL_MIN = -12.f;
    static constexpr float RATE_TOTAL_MAX = 4.f;
    static constexpr float CV_DEPTH = 0.1f;
    static constexpr float CLIMATE_DEPTH = 0.1f;
    static constexpr float TRAIT_VOLTS = 2.f;
    static constexpr float OUTPUT_TAU = 0.02f;
    static constexpr float EVENT_TIME = 1e-3f;
    static constexpr float EVENT_LEVEL = 10.f;
    static constexpr float EVENT_GLOW_TIME = 1.1f;

    coalescent::ArchipelagoField field;
    dsp::ClockDivider fieldDivider;
    dsp::ClockDivider displayDivider;
    dsp::SchmittTrigger resetInputTrigger;
    dsp::SchmittTrigger resetButtonTrigger;
    dsp::PulseGenerator colonizePulse;
    dsp::PulseGenerator extinctPulse;

    bool pendingReset = false;
    float lastSampleRate = 0.f;
    float fieldWallStep = 1.f / FIELD_HZ;
    float outputAlpha = 1.f;
    float target[SMOOTH_LEN] = {};
    float smooth[SMOOTH_LEN] = {};
    float colonizeGlow[HABITATS] = {};
    float extinctGlow[HABITATS] = {};
    coalescent::ArchipelagoField::Parameters currentParameters;

    struct DisplayFrame {
        float density[VALUES] = {};
        float environment[HABITATS] = {};
        float trait[HABITATS] = {};
        float mass[HABITATS] = {};
        float colonizeGlow[HABITATS] = {};
        float extinctGlow[HABITATS] = {};
        float barrier = 0.f;
        std::uint32_t occupiedMask = 0u;
        bool ring = false;

        DisplayFrame() {
            coalescent::ArchipelagoField preview;
            coalescent::ArchipelagoField::Parameters p;
            preview.copyMasses(density, VALUES);
            const coalescent::ArchipelagoField::Metrics& m = preview.metrics();
            for (int i = 0; i < HABITATS; ++i) {
                environment[i] = static_cast<float>(
                    p.climate + p.gradient * coalescent::ArchipelagoField::habitatPosition(i));
                trait[i] = static_cast<float>(m.trait[i]);
                mass[i] = static_cast<float>(m.mass[i]);
            }
            barrier = static_cast<float>(p.barrier);
            occupiedMask = m.occupiedMask;
            ring = p.ring;
        }
    };
    coalescent::DisplaySnapshot<DisplayFrame> displaySnapshot;

    struct SaveFrame {
        coalescent::ArchipelagoField::State state;
        bool valid = false;
    };
    coalescent::DisplaySnapshot<SaveFrame> saveSnapshot;

    Archipelago() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);
        configParam(RATE_PARAM, -8.f, 4.f, 0.f,
                    "Evolution rate", "x", 2.f, 1.f);
        configParam<ArchipelagoSelectionQuantity>(
                    SELECT_PARAM, 0.f, 1.f, 0.60f,
                    "Local selection toward each habitat optimum", " / tau");
        configParam<ArchipelagoMutationQuantity>(
                    MUTATE_PARAM, 0.f, 1.f, 0.45f,
                    "Mutation diffusion in trait space", " / tau");
        configParam<ArchipelagoMigrationQuantity>(
                    MIGRATE_PARAM, 0.f, 1.f, 0.40f,
                    "Nearest-neighbor migration", " / tau");
        configParam(GRADIENT_PARAM, -0.85f, 0.85f, 0.55f,
                    "Environmental gradient coefficient (habitat 1 to 8 span = 2x)",
                    " trait");
        configParam(BARRIER_PARAM, 0.f, 1.f, 0.f,
                    "Central migration barrier (permeability = (1 - barrier)^2)");
        configParam(CLIMATE_PARAM, -0.85f, 0.85f, 0.f,
                    "Shared environmental optimum offset", " trait");
        configButton(RESET_PARAM,
                     "Reset population field at current gradient and climate");
        configParam(SELECT_ATT_PARAM, -1.f, 1.f, 0.f, "Selection CV");
        configParam(MUTATE_ATT_PARAM, -1.f, 1.f, 0.f, "Mutation CV");
        configParam(MIGRATE_ATT_PARAM, -1.f, 1.f, 0.f, "Migration CV");
        configParam(GRADIENT_ATT_PARAM, -1.f, 1.f, 0.f, "Gradient CV");
        configParam(BARRIER_ATT_PARAM, -1.f, 1.f, 0.f, "Barrier CV");
        configSwitch(TOPOLOGY_PARAM, 0.f, 1.f, 0.f, "Topology", {"Row", "Ring"});

        configInput(RATE_INPUT, "Evolution rate (1 V doubles speed)");
        configInput(SELECT_INPUT, "Selection CV");
        configInput(MUTATE_INPUT, "Mutation CV");
        configInput(MIGRATE_INPUT, "Migration CV");
        configInput(GRADIENT_INPUT, "Environmental gradient CV");
        configInput(BARRIER_INPUT, "Central barrier CV");
        configInput(CLIMATE_INPUT, "Climate offset CV (0.1 trait/V)");
        configInput(RESET_INPUT, "Reset trigger");

        configOutput(TRAIT_OUTPUT,
                     "Local mean trait (8 channels; 2 V/trait; holds while empty)");
        configOutput(MASS_OUTPUT, "Local population mass (8 channels; 0 to 10 V)");
        configOutput(MEAN_OUTPUT, "Global abundance-weighted trait mean (2 V/trait)");
        configOutput(DIFF_OUTPUT,
                     "Right-half minus left-half trait mean (2 V/trait difference)");
        configOutput(FLUX_OUTPUT,
                     "Current nearest-neighbor migration flux (0 to 10 V)");
        configOutput(COLONIZE_OUTPUT,
                     "Any population established trigger (10 V, approximately 1 ms)");
        configOutput(EXTINCT_OUTPUT,
                     "Any population lost trigger (10 V, approximately 1 ms)");

        currentParameters = readParameters();
        field.reset(currentParameters);
        setOutputTargets(true);
        publishDisplayFrame();
        publishSaveFrame();
    }

    float inputVoltage(int id) {
        return inputs[id].isConnected()
            ? archipelagoSafe(inputs[id].getVoltage()) : 0.f;
    }

    coalescent::ArchipelagoField::Parameters readParameters() {
        coalescent::ArchipelagoField::Parameters p;
        const float selectionUnit = archipelagoClamp(
            archipelagoSafe(params[SELECT_PARAM].getValue(), 0.60f)
                + inputVoltage(SELECT_INPUT)
                    * archipelagoSafe(params[SELECT_ATT_PARAM].getValue()) * CV_DEPTH,
            0.f, 1.f);
        const float mutationUnit = archipelagoClamp(
            archipelagoSafe(params[MUTATE_PARAM].getValue(), 0.45f)
                + inputVoltage(MUTATE_INPUT)
                    * archipelagoSafe(params[MUTATE_ATT_PARAM].getValue()) * CV_DEPTH,
            0.f, 1.f);
        const float migrationUnit = archipelagoClamp(
            archipelagoSafe(params[MIGRATE_PARAM].getValue(), 0.40f)
                + inputVoltage(MIGRATE_INPUT)
                    * archipelagoSafe(params[MIGRATE_ATT_PARAM].getValue()) * CV_DEPTH,
            0.f, 1.f);
        const float gradient = archipelagoClamp(
            archipelagoSafe(params[GRADIENT_PARAM].getValue(), 0.55f)
                + inputVoltage(GRADIENT_INPUT)
                    * archipelagoSafe(params[GRADIENT_ATT_PARAM].getValue()) * CV_DEPTH,
            -0.85f, 0.85f);
        const float barrier = archipelagoClamp(
            archipelagoSafe(params[BARRIER_PARAM].getValue())
                + inputVoltage(BARRIER_INPUT)
                    * archipelagoSafe(params[BARRIER_ATT_PARAM].getValue()) * CV_DEPTH,
            0.f, 1.f);
        const float climate = archipelagoClamp(
            archipelagoSafe(params[CLIMATE_PARAM].getValue())
                + inputVoltage(CLIMATE_INPUT) * CLIMATE_DEPTH,
            -0.85f, 0.85f);

        p.selection = ArchipelagoSelectionQuantity::fromUnit(selectionUnit);
        p.mutation = ArchipelagoMutationQuantity::fromUnit(mutationUnit);
        p.migration = ArchipelagoMigrationQuantity::fromUnit(migrationUnit);
        p.gradient = gradient;
        p.barrier = barrier;
        p.climate = climate;
        p.ring = archipelagoSafe(params[TOPOLOGY_PARAM].getValue()) >= 0.5f;
        return p;
    }

    void setOutputTargets(bool immediate = false) {
        const coalescent::ArchipelagoField::Metrics& m = field.metrics();
        for (int i = 0; i < HABITATS; ++i) {
            target[TRAIT_SMOOTH + i] = archipelagoClamp(
                static_cast<float>(m.trait[i]) * TRAIT_VOLTS, -2.f, 2.f);
            target[MASS_SMOOTH + i] = archipelagoClamp(
                static_cast<float>(m.mass[i]) * 10.f, 0.f, 10.f);
        }
        target[MEAN_SMOOTH] = archipelagoClamp(
            static_cast<float>(m.globalMean) * TRAIT_VOLTS, -2.f, 2.f);
        target[DIFF_SMOOTH] = archipelagoClamp(
            static_cast<float>(m.difference) * TRAIT_VOLTS, -4.f, 4.f);
        target[FLUX_SMOOTH] = archipelagoClamp(
            static_cast<float>(m.flux) * 10.f, 0.f, 10.f);
        if (immediate)
            std::copy(target, target + SMOOTH_LEN, smooth);
    }

    void clearEventState() {
        colonizePulse.reset();
        extinctPulse.reset();
        std::fill(colonizeGlow, colonizeGlow + HABITATS, 0.f);
        std::fill(extinctGlow, extinctGlow + HABITATS, 0.f);
    }

    void resetField(bool clearEdges) {
        currentParameters = readParameters();
        field.reset(currentParameters);
        pendingReset = false;
        clearEventState();
        if (clearEdges) {
            resetInputTrigger.reset();
            resetButtonTrigger.reset();
        }
        setOutputTargets(true);
        publishDisplayFrame();
        publishSaveFrame();
    }

    void onReset(const ResetEvent& event) override {
        Module::onReset(event);
        resetField(true);
    }

    void publishDisplayFrame() {
        DisplayFrame& d = displaySnapshot.writable();
        field.copyMasses(d.density, VALUES);
        const coalescent::ArchipelagoField::Metrics& m = field.metrics();
        for (int i = 0; i < HABITATS; ++i) {
            d.environment[i] = static_cast<float>(
                currentParameters.climate
                + currentParameters.gradient
                    * coalescent::ArchipelagoField::habitatPosition(i));
            d.trait[i] = static_cast<float>(m.trait[i]);
            d.mass[i] = static_cast<float>(m.mass[i]);
            d.colonizeGlow[i] = colonizeGlow[i];
            d.extinctGlow[i] = extinctGlow[i];
        }
        d.barrier = static_cast<float>(currentParameters.barrier);
        d.occupiedMask = m.occupiedMask;
        d.ring = currentParameters.ring;
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

        json_object_set_new(root, "archipelagoVersion", json_integer(1));
        json_t* density = json_array();
        for (int i = 0; i < VALUES; ++i)
            json_array_append_new(density, json_real(s.state.density[i]));
        json_object_set_new(root, "density", density);
        json_t* traits = json_array();
        for (int i = 0; i < HABITATS; ++i)
            json_array_append_new(traits, json_real(s.state.reportedTrait[i]));
        json_object_set_new(root, "reportedTrait", traits);
        json_object_set_new(root, "reportedGlobalMean",
                            json_real(s.state.reportedGlobalMean));
        json_object_set_new(root, "occupiedMask",
                            json_integer(s.state.occupiedMask));
        return root;
    }

    void dataFromJson(json_t* root) override {
        json_t* version = json_object_get(root, "archipelagoVersion");
        json_t* density = json_object_get(root, "density");
        json_t* traits = json_object_get(root, "reportedTrait");
        json_t* globalMean = json_object_get(root, "reportedGlobalMean");
        json_t* occupiedMask = json_object_get(root, "occupiedMask");
        if (!version || !json_is_integer(version)
            || json_integer_value(version) != 1
            || !density || !json_is_array(density)
            || json_array_size(density) != static_cast<std::size_t>(VALUES)
            || !traits || !json_is_array(traits)
            || json_array_size(traits) != static_cast<std::size_t>(HABITATS)
            || !globalMean || !json_is_number(globalMean)
            || !occupiedMask || !json_is_integer(occupiedMask))
            return;

        coalescent::ArchipelagoField::State restored;
        for (int i = 0; i < VALUES; ++i) {
            json_t* item = json_array_get(density, i);
            if (!item || !json_is_number(item))
                return;
            restored.density[i] = json_number_value(item);
        }
        for (int i = 0; i < HABITATS; ++i) {
            json_t* item = json_array_get(traits, i);
            if (!item || !json_is_number(item))
                return;
            restored.reportedTrait[i] = json_number_value(item);
        }
        restored.reportedGlobalMean = json_number_value(globalMean);
        const json_int_t mask = json_integer_value(occupiedMask);
        if (mask < 0 || mask > 0xff)
            return;
        restored.occupiedMask = static_cast<std::uint32_t>(mask);

        currentParameters = readParameters();
        coalescent::ArchipelagoField candidate;
        if (!candidate.restore(restored, currentParameters))
            return;

        field = candidate;
        pendingReset = false;
        resetInputTrigger.reset();
        resetButtonTrigger.reset();
        clearEventState();
        setOutputTargets(true);
        publishDisplayFrame();
        publishSaveFrame();
    }

    void collectEvents() {
        const coalescent::ArchipelagoField::Metrics& m = field.metrics();
        if (m.colonizeEvent)
            colonizePulse.trigger(EVENT_TIME);
        if (m.extinctEvent)
            extinctPulse.trigger(EVENT_TIME);
        for (int i = 0; i < HABITATS; ++i) {
            const std::uint32_t bit = UINT32_C(1) << i;
            if ((m.colonizeMask & bit) != 0u)
                colonizeGlow[i] = 1.f;
            if ((m.extinctMask & bit) != 0u)
                extinctGlow[i] = 1.f;
        }
    }

    void process(const ProcessArgs& args) override {
        if (args.sampleRate != lastSampleRate) {
            lastSampleRate = args.sampleRate;
            const int fieldDivision = std::max(
                1, static_cast<int>(std::round(args.sampleRate / FIELD_HZ)));
            fieldDivider.setDivision(fieldDivision);
            fieldDivider.reset();
            fieldWallStep = fieldDivision / args.sampleRate;
            displayDivider.setDivision(std::max(
                1, static_cast<int>(std::round(args.sampleRate / 45.f))));
            displayDivider.reset();
            outputAlpha = 1.f - std::exp(-args.sampleTime / OUTPUT_TAU);
        }

        const bool resetInputEdge = resetInputTrigger.process(
            inputVoltage(RESET_INPUT), 0.1f, 1.f);
        const bool resetButtonEdge = resetButtonTrigger.process(
            archipelagoSafe(params[RESET_PARAM].getValue()), 0.1f, 0.9f);
        if (resetInputEdge || resetButtonEdge)
            pendingReset = true;

        for (int i = 0; i < HABITATS; ++i) {
            colonizeGlow[i] = std::max(
                0.f, colonizeGlow[i] - args.sampleTime / EVENT_GLOW_TIME);
            extinctGlow[i] = std::max(
                0.f, extinctGlow[i] - args.sampleTime / EVENT_GLOW_TIME);
        }

        if (fieldDivider.process()) {
            currentParameters = readParameters();
            if (pendingReset) {
                field.reset(currentParameters);
                pendingReset = false;
                clearEventState();
            }
            else {
                const float rateOct = archipelagoClamp(
                    archipelagoSafe(params[RATE_PARAM].getValue())
                        + inputVoltage(RATE_INPUT),
                    RATE_TOTAL_MIN, RATE_TOTAL_MAX);
                const double deltaTau = static_cast<double>(RATE_BASE)
                    * std::exp2(static_cast<double>(rateOct))
                    * static_cast<double>(fieldWallStep);
                field.advance(deltaTau, currentParameters);
                collectEvents();
            }
            setOutputTargets();
            publishSaveFrame();
        }

        for (int i = 0; i < SMOOTH_LEN; ++i)
            smooth[i] += (target[i] - smooth[i]) * outputAlpha;

        outputs[TRAIT_OUTPUT].setChannels(HABITATS);
        outputs[MASS_OUTPUT].setChannels(HABITATS);
        for (int i = 0; i < HABITATS; ++i) {
            outputs[TRAIT_OUTPUT].setVoltage(smooth[TRAIT_SMOOTH + i], i);
            outputs[MASS_OUTPUT].setVoltage(smooth[MASS_SMOOTH + i], i);
        }
        outputs[MEAN_OUTPUT].setVoltage(smooth[MEAN_SMOOTH]);
        outputs[DIFF_OUTPUT].setVoltage(smooth[DIFF_SMOOTH]);
        outputs[FLUX_OUTPUT].setVoltage(smooth[FLUX_SMOOTH]);
        outputs[COLONIZE_OUTPUT].setVoltage(
            colonizePulse.process(args.sampleTime) ? EVENT_LEVEL : 0.f);
        outputs[EXTINCT_OUTPUT].setVoltage(
            extinctPulse.process(args.sampleTime) ? EVENT_LEVEL : 0.f);

        if (displayDivider.process())
            publishDisplayFrame();
    }
};

namespace archipelago_layout {
static const float KNOB_Y = 58.f;
static const float AUX_Y = 70.f;
static const float INPUT_Y = 84.5f;
static const float OUTPUT_Y = 101.f;
static const float EVENT_Y = 116.5f;
static const float KNOB_X[7] = {7.8f, 20.4f, 33.1f, 45.7f,
                                58.3f, 71.f, 83.6f};
static const float INPUT_X[8] = {6.4f, 17.65f, 28.9f, 40.15f,
                                 51.4f, 62.65f, 73.9f, 85.15f};
static const float OUTPUT_X[5] = {10.f, 27.85f, 45.72f, 63.58f, 81.44f};
static const float EVENT_X[2] = {36.f, 55.44f};
} // namespace archipelago_layout

struct ArchipelagoFieldView : widget::TransparentWidget {
    Archipelago* module = nullptr;

    static NVGcolor habitatColor(int habitat, unsigned alpha) {
        static const unsigned colors[Archipelago::HABITATS][3] = {
            {0x62, 0xd1, 0xbd}, {0x73, 0xd0, 0xb0},
            {0x8b, 0xcc, 0x9c}, {0xae, 0xc5, 0x86},
            {0xce, 0xb8, 0x78}, {0xdf, 0xa3, 0x78},
            {0xdc, 0x8e, 0x87}, {0xca, 0x82, 0xa3}
        };
        habitat = clamp(habitat, 0, Archipelago::HABITATS - 1);
        return nvgRGBA(colors[habitat][0], colors[habitat][1],
                       colors[habitat][2], alpha);
    }

    void draw(const DrawArgs& args) override {
        nvgBeginPath(args.vg);
        nvgRoundedRect(args.vg, 0.f, 0.f, box.size.x, box.size.y, mm2px(1.7f));
        nvgFillColor(args.vg, nvgRGB(0x06, 0x12, 0x14));
        nvgFill(args.vg);
        nvgStrokeColor(args.vg, nvgRGB(0x35, 0x54, 0x53));
        nvgStrokeWidth(args.vg, 1.2f);
        nvgStroke(args.vg);
        Widget::draw(args);
    }

    void drawLayer(const DrawArgs& args, int layer) override {
        if (layer != 1)
            return;
        static const Archipelago::DisplayFrame preview;
        const Archipelago::DisplayFrame& frame = module
            ? module->displaySnapshot.consume() : preview;
        const float width = box.size.x;
        const float height = box.size.y;
        const float left = mm2px(2.4f);
        const float right = width - mm2px(2.4f);
        const float top = mm2px(7.2f);
        const float bottom = height - mm2px(2.8f);
        const float plotWidth = right - left;
        const float plotHeight = bottom - top;
        const float habitatWidth = plotWidth / Archipelago::HABITATS;
        const float binHeight = plotHeight / Archipelago::BINS;
        const float visibleBinMass = 0.32f;
        auto yForTrait = [&](float trait) {
            float unit = (trait - static_cast<float>(
                coalescent::ArchipelagoField::traitMin()))
                / static_cast<float>(coalescent::ArchipelagoField::traitMax()
                                     - coalescent::ArchipelagoField::traitMin());
            return bottom - clamp(unit, 0.f, 1.f) * plotHeight;
        };

        nvgScissor(args.vg, 0.f, 0.f, width, height);
        for (int habitat = 0; habitat < Archipelago::HABITATS; ++habitat) {
            const float x = left + habitat * habitatWidth;
            for (int bin = 0; bin < Archipelago::BINS; ++bin) {
                const int index = habitat * Archipelago::BINS + bin;
                const float level = std::sqrt(clamp(
                    frame.density[index] / visibleBinMass, 0.f, 1.f));
                if (!(level > 0.008f))
                    continue;
                const unsigned alpha = static_cast<unsigned>(0x18 + level * 0xd0);
                nvgBeginPath(args.vg);
                nvgRect(args.vg, x + 0.5f,
                        bottom - (bin + 1) * binHeight + 0.25f,
                        habitatWidth - 1.f, std::max(0.5f, binHeight - 0.4f));
                nvgFillColor(args.vg, habitatColor(habitat, alpha));
                nvgFill(args.vg);
            }

            const float colonize = clamp(frame.colonizeGlow[habitat], 0.f, 1.f);
            const float extinct = clamp(frame.extinctGlow[habitat], 0.f, 1.f);
            if (colonize > 0.f || extinct > 0.f) {
                nvgBeginPath(args.vg);
                nvgRect(args.vg, x, top, habitatWidth, plotHeight);
                if (colonize >= extinct)
                    nvgFillColor(args.vg, nvgRGBA(
                        0x72, 0xe0, 0xa5, static_cast<unsigned>(colonize * 0x48)));
                else
                    nvgFillColor(args.vg, nvgRGBA(
                        0xf0, 0x78, 0x72, static_cast<unsigned>(extinct * 0x58)));
                nvgFill(args.vg);
            }

            const bool occupied = (frame.occupiedMask & (UINT32_C(1) << habitat)) != 0u;
            if (!occupied) {
                const float cx = x + habitatWidth * 0.5f;
                const float cy = top + plotHeight * 0.5f;
                nvgBeginPath(args.vg);
                nvgMoveTo(args.vg, cx - mm2px(0.8f), cy - mm2px(0.8f));
                nvgLineTo(args.vg, cx + mm2px(0.8f), cy + mm2px(0.8f));
                nvgMoveTo(args.vg, cx + mm2px(0.8f), cy - mm2px(0.8f));
                nvgLineTo(args.vg, cx - mm2px(0.8f), cy + mm2px(0.8f));
                nvgStrokeColor(args.vg, nvgRGBA(0xd9, 0xe3, 0xde, 0x38));
                nvgStrokeWidth(args.vg, 1.f);
                nvgStroke(args.vg);

                nvgBeginPath(args.vg);
                nvgCircle(args.vg, cx, yForTrait(frame.trait[habitat]),
                          mm2px(0.48f));
                nvgStrokeColor(args.vg, nvgRGBA(0xee, 0xf5, 0xec, 0x58));
                nvgStrokeWidth(args.vg, 1.f);
                nvgStroke(args.vg);
            }
        }

        nvgBeginPath(args.vg);
        for (int habitat = 0; habitat < Archipelago::HABITATS; ++habitat) {
            const float x = left + (habitat + 0.5f) * habitatWidth;
            const float y = yForTrait(frame.environment[habitat]);
            if (habitat == 0) nvgMoveTo(args.vg, x, y);
            else              nvgLineTo(args.vg, x, y);
        }
        nvgStrokeColor(args.vg, nvgRGBA(0xf0, 0xc5, 0x67, 0x82));
        nvgStrokeWidth(args.vg, 0.9f);
        nvgStroke(args.vg);
        for (int habitat = 0; habitat < Archipelago::HABITATS; ++habitat) {
            const float environment = frame.environment[habitat];
            if (environment >= coalescent::ArchipelagoField::traitMin()
                && environment <= coalescent::ArchipelagoField::traitMax())
                continue;
            const float x = left + (habitat + 0.5f) * habitatWidth;
            const bool above = environment > coalescent::ArchipelagoField::traitMax();
            const float tipY = above ? top : bottom;
            const float baseY = tipY + (above ? mm2px(1.1f) : -mm2px(1.1f));
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, x, tipY);
            nvgLineTo(args.vg, x - mm2px(0.7f), baseY);
            nvgLineTo(args.vg, x + mm2px(0.7f), baseY);
            nvgClosePath(args.vg);
            nvgFillColor(args.vg, nvgRGBA(0xf0, 0xc5, 0x67, 0x98));
            nvgFill(args.vg);
        }

        for (int habitat = 0; habitat < Archipelago::HABITATS; ++habitat) {
            const bool occupied = (frame.occupiedMask & (UINT32_C(1) << habitat)) != 0u;
            if (!occupied)
                continue;
            const float x = left + (habitat + 0.5f) * habitatWidth;
            const float y = yForTrait(frame.trait[habitat]);
            if (habitat + 1 < Archipelago::HABITATS
                && (frame.occupiedMask & (UINT32_C(1) << (habitat + 1))) != 0u) {
                const float nextX = left + (habitat + 1.5f) * habitatWidth;
                const float nextY = yForTrait(frame.trait[habitat + 1]);
                nvgBeginPath(args.vg);
                nvgMoveTo(args.vg, x, y);
                nvgLineTo(args.vg, nextX, nextY);
                nvgStrokeColor(args.vg, nvgRGBA(0xee, 0xf5, 0xec, 0xc8));
                nvgStrokeWidth(args.vg, 1.4f);
                nvgStroke(args.vg);
            }
            nvgBeginPath(args.vg);
            nvgCircle(args.vg, x, y, mm2px(0.55f));
            nvgFillColor(args.vg, nvgRGBA(0xf1, 0xf7, 0xf2, 0xe8));
            nvgFill(args.vg);
        }

        const float barrierX = left + 4.f * habitatWidth;
        nvgBeginPath(args.vg);
        nvgMoveTo(args.vg, barrierX, top);
        nvgLineTo(args.vg, barrierX, bottom);
        nvgStrokeColor(args.vg, nvgRGBA(
            0xf0, 0x8b, 0x73,
            static_cast<unsigned>(0x24 + clamp(frame.barrier, 0.f, 1.f) * 0xc0)));
        nvgStrokeWidth(args.vg, 1.f + 3.f * clamp(frame.barrier, 0.f, 1.f));
        nvgStroke(args.vg);

        if (frame.ring) {
            const float firstX = left + 0.5f * habitatWidth;
            const float lastX = right - 0.5f * habitatWidth;
            const float wrapY = bottom + mm2px(1.5f);
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, firstX, bottom - 0.5f);
            nvgBezierTo(args.vg, firstX, wrapY, lastX, wrapY,
                        lastX, bottom - 0.5f);
            nvgStrokeColor(args.vg, nvgRGBA(0x78, 0xca, 0xd5, 0x80));
            nvgStrokeWidth(args.vg, 1.1f);
            nvgStroke(args.vg);
        }
        nvgResetScissor(args.vg);

        std::shared_ptr<Font> font = APP->window->loadFont(
            asset::system("res/fonts/DejaVuSans.ttf"));
        if (font) {
            nvgFontFaceId(args.vg, font->handle);
            nvgFontSize(args.vg, mm2px(3.0f));
            nvgTextLetterSpacing(args.vg, mm2px(0.38f));
            nvgFillColor(args.vg, nvgRGBA(0xe9, 0xf1, 0xe9, 0xe8));
            nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
            nvgText(args.vg, mm2px(2.6f), mm2px(2.1f), "ARCHIPELAGO", nullptr);
            nvgTextLetterSpacing(args.vg, 0.f);
            nvgFontSize(args.vg, mm2px(1.9f));
            nvgFillColor(args.vg, nvgRGBA(0x98, 0xbe, 0xb5, 0xb0));
            nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_TOP);
            nvgText(args.vg, width - mm2px(2.5f), mm2px(2.6f),
                    frame.ring ? "RING" : "ROW", nullptr);
        }
    }
};

struct ArchipelagoLabels : widget::Widget {
    void draw(const DrawArgs& args) override {
        std::shared_ptr<Font> font = APP->window->loadFont(
            asset::system("res/fonts/Nunito-Bold.ttf"));
        if (!font)
            return;
        nvgFontFaceId(args.vg, font->handle);
        nvgTextLetterSpacing(args.vg, 0.f);
        nvgFillColor(args.vg, nvgRGB(0xe7, 0xf1, 0xed));
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
        auto label = [&](float x, float y, const char* text, float size = 1.3f) {
            nvgFontSize(args.vg, mm2px(size * 1.72f));
            nvgText(args.vg, mm2px(x), mm2px(y), text, nullptr);
        };
        using namespace archipelago_layout;
        const char* knobs[7] = {
            "RATE", "SELECT", "MUTATE", "MIGRATE", "GRADIENT", "BARRIER", "CLIMATE"
        };
        for (int i = 0; i < 7; ++i)
            label(KNOB_X[i], KNOB_Y - 8.f, knobs[i], i >= 4 ? 1.08f : 1.22f);
        label(KNOB_X[0], AUX_Y - 5.1f, "RESET", 1.15f);
        for (int i = 1; i <= 5; ++i)
            label(KNOB_X[i], AUX_Y - 5.1f, "CV", 1.08f);
        label(KNOB_X[6], AUX_Y - 5.1f, "ROW/RING", 0.92f);

        const char* inputs[8] = {
            "RATE", "SEL", "MUT", "MIG", "GRAD", "BARR", "CLIM", "RESET"
        };
        for (int i = 0; i < 8; ++i)
            label(INPUT_X[i], INPUT_Y - 5.5f, inputs[i], 1.08f);
        const char* outputs[5] = {"TRAIT", "MASS", "MEAN", "DIFF", "FLUX"};
        for (int i = 0; i < 5; ++i)
            label(OUTPUT_X[i], OUTPUT_Y - 5.5f, outputs[i], 1.25f);
        label(EVENT_X[0], EVENT_Y - 5.5f, "COLONIZE", 1.12f);
        label(EVENT_X[1], EVENT_Y - 5.5f, "EXTINCT", 1.15f);
    }
};

struct ArchipelagoWidget : ModuleWidget {
    ArchipelagoWidget(Archipelago* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Archipelago.svg")));
        addPanelLabels<ArchipelagoLabels>(this);

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.f, 1.f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(85.44f, 1.f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.f, 122.f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(85.44f, 122.f))));

        ArchipelagoFieldView* display = new ArchipelagoFieldView;
        display->module = module;
        display->box.pos = mm2px(Vec(6.f, 8.f));
        display->box.size = mm2px(Vec(79.44f, 38.f));
        addChild(display);

        using namespace archipelago_layout;
        const int knobs[7] = {
            Archipelago::RATE_PARAM, Archipelago::SELECT_PARAM,
            Archipelago::MUTATE_PARAM, Archipelago::MIGRATE_PARAM,
            Archipelago::GRADIENT_PARAM, Archipelago::BARRIER_PARAM,
            Archipelago::CLIMATE_PARAM
        };
        for (int i = 0; i < 7; ++i)
            addParam(createParamCentered<RoundBlackKnob>(
                mm2px(Vec(KNOB_X[i], KNOB_Y)), module, knobs[i]));
        addParam(createParamCentered<TL1105>(
            mm2px(Vec(KNOB_X[0], AUX_Y)), module, Archipelago::RESET_PARAM));
        const int attenuverters[5] = {
            Archipelago::SELECT_ATT_PARAM, Archipelago::MUTATE_ATT_PARAM,
            Archipelago::MIGRATE_ATT_PARAM, Archipelago::GRADIENT_ATT_PARAM,
            Archipelago::BARRIER_ATT_PARAM
        };
        for (int i = 0; i < 5; ++i)
            addParam(createParamCentered<Trimpot>(
                mm2px(Vec(KNOB_X[i + 1], AUX_Y)), module, attenuverters[i]));
        addParam(createParamCentered<CKSS>(
            mm2px(Vec(KNOB_X[6], AUX_Y)), module, Archipelago::TOPOLOGY_PARAM));

        const int inputs[8] = {
            Archipelago::RATE_INPUT, Archipelago::SELECT_INPUT,
            Archipelago::MUTATE_INPUT, Archipelago::MIGRATE_INPUT,
            Archipelago::GRADIENT_INPUT, Archipelago::BARRIER_INPUT,
            Archipelago::CLIMATE_INPUT, Archipelago::RESET_INPUT
        };
        for (int i = 0; i < 8; ++i)
            addInput(createInputCentered<PJ301MPort>(
                mm2px(Vec(INPUT_X[i], INPUT_Y)), module, inputs[i]));

        const int outputs[5] = {
            Archipelago::TRAIT_OUTPUT, Archipelago::MASS_OUTPUT,
            Archipelago::MEAN_OUTPUT, Archipelago::DIFF_OUTPUT,
            Archipelago::FLUX_OUTPUT
        };
        for (int i = 0; i < 5; ++i)
            addOutput(createOutputCentered<PJ301MPort>(
                mm2px(Vec(OUTPUT_X[i], OUTPUT_Y)), module, outputs[i]));
        addOutput(createOutputCentered<PJ301MPort>(
            mm2px(Vec(EVENT_X[0], EVENT_Y)), module, Archipelago::COLONIZE_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(
            mm2px(Vec(EVENT_X[1], EVENT_Y)), module, Archipelago::EXTINCT_OUTPUT));
    }
};

Model* modelArchipelago = createModel<Archipelago, ArchipelagoWidget>("Archipelago");
