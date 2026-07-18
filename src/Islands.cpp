#include "plugin.hpp"
#include "dsp/display_snapshot.hpp"
#include "dsp/hex64.hpp"
#include "dsp/islands_model.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>

// ISLANDS - four finite populations under Wright-Fisher drift.
//
// Each generation applies selection, symmetric mutation, and common-pool
// migration before independently sampling four new populations. The exact,
// deterministic population model lives in dsp/islands_model.hpp. This wrapper
// supplies clocking, audio-rate interpolation, trigger/gate semantics, Rack
// persistence, and the fixed-scale population-history display.

namespace {
constexpr std::uint64_t ISLANDS_FACTORY_SEED = UINT64_C(42);
constexpr std::uint64_t ISLANDS_FACTORY_STREAM = UINT64_C(54);
}

struct IslandsSizeQuantity : ParamQuantity {
    float getDisplayValue() override {
        Param* parameter = getParam();
        float value = parameter ? parameter->getValue() : defaultValue;
        return std::round(std::exp2(value));
    }

    void setDisplayValue(float displayValue) override {
        if (!std::isfinite(displayValue))
            return;
        displayValue = clamp(displayValue, 8.f, 4096.f);
        setImmediateValue(std::log2(displayValue));
    }
};

struct IslandsMutationQuantity : ParamQuantity {
    static float fromUnit(float unitValue) {
        unitValue = clamp(unitValue, 0.f, 1.f);
        if (!(unitValue > 0.f))
            return 0.f;
        static const float positiveMinimum = 1e-6f;
        static const float maximum = 0.05f;
        return positiveMinimum * std::pow(maximum / positiveMinimum, unitValue);
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
        static const float positiveMinimum = 1e-6f;
        static const float maximum = 0.05f;
        displayValue = clamp(displayValue, positiveMinimum, maximum);
        setImmediateValue(std::log(displayValue / positiveMinimum)
                          / std::log(maximum / positiveMinimum));
    }
};

struct IslandsGenerationQuantity : ParamQuantity {
    std::string getDisplayValueString() override {
        Param* parameter = getParam();
        if ((parameter ? parameter->getValue() : defaultValue) <= minValue)
            return "stopped";
        return ParamQuantity::getDisplayValueString();
    }

    std::string getUnit() override {
        Param* parameter = getParam();
        return (parameter ? parameter->getValue() : defaultValue) <= minValue
            ? "" : ParamQuantity::getUnit();
    }
};

struct Islands : Module {
    enum SeedAction {
        SEED_ACTION_NORMAL,
        SEED_ACTION_NEW_RANDOM,
        SEED_ACTION_FACTORY
    };

    enum ParamId {
        SIZE_PARAM,
        SELECT_PARAM,
        MUTATE_PARAM,
        MIGRATE_PARAM,
        GENERATIONS_PARAM,
        FOUNDER_PARAM,
        RESET_PARAM,
        SELECT_ATT_PARAM,
        MUTATE_ATT_PARAM,
        MIGRATE_ATT_PARAM,
        PARAMS_LEN
    };
    enum InputId {
        SIZE_INPUT,
        SELECT_INPUT,
        MUTATE_INPUT,
        MIGRATE_INPUT,
        GENERATIONS_INPUT,
        STEP_INPUT,
        FOUNDER_INPUT,
        RESET_INPUT,
        INPUTS_LEN
    };
    enum OutputId {
        I1_OUTPUT,
        I2_OUTPUT,
        I3_OUTPUT,
        I4_OUTPUT,
        MEAN_OUTPUT,
        HET_OUTPUT,
        FIX_A_OUTPUT,
        FIX_B_OUTPUT,
        LOSS_OUTPUT,
        SWEEP_OUTPUT,
        OUTPUTS_LEN
    };
    enum LightId { LIGHTS_LEN };

    enum ContinuousOutput {
        I1_CONTINUOUS,
        I2_CONTINUOUS,
        I3_CONTINUOUS,
        I4_CONTINUOUS,
        MEAN_CONTINUOUS,
        HET_CONTINUOUS,
        CONTINUOUS_LEN
    };

    static const int HISTORY_LEN = 128;

    coalescent::IslandsModel model;
    dsp::SchmittTrigger stepTrigger;
    dsp::SchmittTrigger founderInputTrigger;
    dsp::SchmittTrigger founderButtonTrigger;
    dsp::SchmittTrigger resetInputTrigger;
    dsp::SchmittTrigger resetButtonTrigger;
    dsp::PulseGenerator lossPulse;
    dsp::PulseGenerator sweepPulse;
    dsp::ClockDivider displayDivider;
    dsp::ClockDivider saveDivider;

    // Double precision matters here: at the slowest rate, a float phase can
    // stop changing once its per-sample increment falls below the local ULP.
    double generationPhase = 0.0;
    double rampPhase = 1.0;
    double rampDuration = 0.05;
    float rampFrom[CONTINUOUS_LEN] = {};
    float rampTo[CONTINUOUS_LEN] = {};
    float rampCurrent[CONTINUOUS_LEN] = {};
    float lastSampleRate = 0.f;
    int lastFounder = -1;
    float founderGlow = 0.f;
    std::uint64_t resetSeed = ISLANDS_FACTORY_SEED;
    std::uint64_t resetStream = ISLANDS_FACTORY_STREAM;
    SeedAction seedAction = SEED_ACTION_NORMAL;

    struct DisplayFrame {
        float history[coalescent::IslandsModel::kIslands][HISTORY_LEN] = {};
        float current[coalescent::IslandsModel::kIslands] = {};
        int head = 0;
        int count = 0;
        int nextFounder = 0;
        int lastFounder = -1;
        float founderGlow = 0.f;

        // Rack's module-browser render has no engine. Give it an informative,
        // deterministic population portrait rather than four empty traces.
        DisplayFrame() {
            static const float phase[4] = {0.2f, 1.4f, 2.7f, 4.3f};
            for (int j = 0; j < HISTORY_LEN; ++j) {
                float t = static_cast<float>(j) / static_cast<float>(HISTORY_LEN - 1);
                for (int i = 0; i < 4; ++i) {
                    float walk = 0.50f
                        + 0.19f * std::sin(5.4f * t + phase[i])
                        + 0.07f * std::sin(17.0f * t + phase[(i + 1) & 3]);
                    history[i][j] = std::max(0.04f, std::min(walk, 0.96f));
                }
            }
            head = 0;
            count = HISTORY_LEN;
            for (int i = 0; i < 4; ++i)
                current[i] = history[i][HISTORY_LEN - 1];
            nextFounder = 1;
            lastFounder = 0;
            founderGlow = 0.55f;
        }
    };

    DisplayFrame liveDisplay;
    coalescent::DisplaySnapshot<DisplayFrame> displaySnapshot;

    struct SaveFrame {
        coalescent::IslandsModel::State modelState;
        std::uint64_t resetSeed = ISLANDS_FACTORY_SEED;
        std::uint64_t resetStream = ISLANDS_FACTORY_STREAM;
        double generationPhase = 0.0;
        double rampPhase = 1.0;
        double rampDuration = 0.05;
        float rampFrom[CONTINUOUS_LEN] = {};
        float rampTo[CONTINUOUS_LEN] = {};
        float rampCurrent[CONTINUOUS_LEN] = {};
        int lastFounder = -1;
        float founderGlow = 0.f;
        bool valid = false;
    };
    coalescent::DisplaySnapshot<SaveFrame> saveSnapshot;

    static constexpr float SIZE_OCT_MIN = 3.f;
    static constexpr float SIZE_OCT_MAX = 12.f;
    static constexpr float GENERATION_OCT_MIN = -8.f;
    static constexpr float GENERATION_OCT_MAX = 7.64385619f; // log2(200)
    static constexpr float SELECT_CV_DEPTH = 0.05f;
    static constexpr float UNIT_CV_DEPTH = 0.1f;
    static constexpr float GATE_LEVEL = 10.f;
    static constexpr float PULSE_TIME = 1e-3f;
    static constexpr float MANUAL_RAMP_TIME = 0.05f;
    static constexpr float FOUNDER_GLOW_TIME = 1.25f;

    Islands() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);
        configParam<IslandsSizeQuantity>(SIZE_PARAM, SIZE_OCT_MIN, SIZE_OCT_MAX, 6.f,
                    "Population size", " copies", 2.f, 1.f);
        configParam(SELECT_PARAM,
                    static_cast<float>(coalescent::IslandsModel::selectionMin()),
                    static_cast<float>(coalescent::IslandsModel::selectionMax()),
                    0.f, "Selection (log-fitness advantage of A)");
        configParam<IslandsMutationQuantity>(MUTATE_PARAM, 0.f, 1.f, 0.35f,
                    "Symmetric mutation probability", " / generation");
        configParam(MIGRATE_PARAM, 0.f, 1.f, 0.15f,
                    "Common-pool migration fraction");
        configParam<IslandsGenerationQuantity>(
                    GENERATIONS_PARAM, GENERATION_OCT_MIN, GENERATION_OCT_MAX, 1.f,
                    "Generation rate", " /s", 2.f, 1.f);
        configButton(FOUNDER_PARAM, "Founder bottleneck");
        configButton(RESET_PARAM, "Reset populations");
        configParam(SELECT_ATT_PARAM, -1.f, 1.f, 0.f, "Selection CV");
        configParam(MUTATE_ATT_PARAM, -1.f, 1.f, 0.f, "Mutation CV");
        configParam(MIGRATE_ATT_PARAM, -1.f, 1.f, 0.f, "Migration CV");

        configInput(SIZE_INPUT, "Population size CV (1 V doubles copies)");
        configInput(SELECT_INPUT, "Selection CV");
        configInput(MUTATE_INPUT, "Mutation CV");
        configInput(MIGRATE_INPUT, "Migration CV");
        configInput(GENERATIONS_INPUT, "Generation rate CV (1 V doubles rate)");
        configInput(STEP_INPUT, "Advance one generation trigger");
        configInput(FOUNDER_INPUT, "Founder bottleneck trigger");
        configInput(RESET_INPUT, "Reset trigger");

        configOutput(I1_OUTPUT, "Island 1 allele-A frequency");
        configOutput(I2_OUTPUT, "Island 2 allele-A frequency");
        configOutput(I3_OUTPUT, "Island 3 allele-A frequency");
        configOutput(I4_OUTPUT, "Island 4 allele-A frequency");
        configOutput(MEAN_OUTPUT, "Mean allele-A frequency");
        configOutput(HET_OUTPUT, "Normalized within-island heterozygosity");
        configOutput(FIX_A_OUTPUT, "Global fixation at A gate");
        configOutput(FIX_B_OUTPUT, "Global fixation at B gate");
        configOutput(LOSS_OUTPUT, "Global allele-loss trigger");
        configOutput(SWEEP_OUTPUT, "Near-fixation sweep trigger");

        initializeState(true);
    }

    static float safe(float value, float fallback = 0.f) {
        return std::isfinite(value) ? value : fallback;
    }

    static float clampf(float value, float low, float high) {
        return std::max(low, std::min(value, high));
    }

    float inputVoltage(int id) {
        return inputs[id].isConnected() ? safe(inputs[id].getVoltage()) : 0.f;
    }

    float mutationFromUnit(float unit) const {
        unit = clampf(unit, 0.f, 1.f);
        if (!(unit > 0.f))
            return 0.f;
        return IslandsMutationQuantity::fromUnit(unit);
    }

    coalescent::IslandsModel::Parameters readGenerationParameters() {
        coalescent::IslandsModel::Parameters result;

        float sizeOct = safe(params[SIZE_PARAM].getValue(), 6.f) + inputVoltage(SIZE_INPUT);
        sizeOct = clampf(sizeOct, SIZE_OCT_MIN, SIZE_OCT_MAX);
        result.copies = static_cast<std::uint32_t>(std::lround(std::exp2(sizeOct)));

        float selection = safe(params[SELECT_PARAM].getValue())
            + inputVoltage(SELECT_INPUT) * safe(params[SELECT_ATT_PARAM].getValue())
                * SELECT_CV_DEPTH;
        result.selection = clampf(selection,
            static_cast<float>(coalescent::IslandsModel::selectionMin()),
            static_cast<float>(coalescent::IslandsModel::selectionMax()));

        float mutationUnit = safe(params[MUTATE_PARAM].getValue(), 0.35f)
            + inputVoltage(MUTATE_INPUT) * safe(params[MUTATE_ATT_PARAM].getValue())
                * UNIT_CV_DEPTH;
        result.mutation = mutationFromUnit(mutationUnit);

        float migration = safe(params[MIGRATE_PARAM].getValue(), 0.15f)
            + inputVoltage(MIGRATE_INPUT) * safe(params[MIGRATE_ATT_PARAM].getValue())
                * UNIT_CV_DEPTH;
        result.migration = clampf(migration,
            static_cast<float>(coalescent::IslandsModel::migrationMin()),
            static_cast<float>(coalescent::IslandsModel::migrationMax()));
        return result;
    }

    float generationRate() {
        float octave = safe(params[GENERATIONS_PARAM].getValue(), 1.f)
            + inputVoltage(GENERATIONS_INPUT);
        if (octave <= GENERATION_OCT_MIN)
            return 0.f;
        if (octave > GENERATION_OCT_MAX)
            octave = GENERATION_OCT_MAX;
        return std::exp2(octave);
    }

    void metricsToVoltages(float* destination) const {
        const coalescent::IslandsModel::Metrics& metrics = model.metrics();
        for (int i = 0; i < 4; ++i)
            destination[i] = clampf(static_cast<float>(metrics.frequency[i]), 0.f, 1.f) * 10.f;
        destination[MEAN_CONTINUOUS] = clampf(static_cast<float>(metrics.mean), 0.f, 1.f) * 10.f;
        destination[HET_CONTINUOUS] = clampf(static_cast<float>(metrics.diversity), 0.f, 1.f) * 10.f;
    }

    void resetHistory() {
        std::memset(liveDisplay.history, 0, sizeof(liveDisplay.history));
        liveDisplay.head = 0;
        liveDisplay.count = 0;
        liveDisplay.nextFounder = 0;
        liveDisplay.lastFounder = -1;
        liveDisplay.founderGlow = 0.f;
        pushHistory();
    }

    void pushHistory() {
        const coalescent::IslandsModel::Metrics& metrics = model.metrics();
        for (int i = 0; i < 4; ++i) {
            float value = clampf(static_cast<float>(metrics.frequency[i]), 0.f, 1.f);
            liveDisplay.history[i][liveDisplay.head] = value;
            liveDisplay.current[i] = value;
        }
        liveDisplay.head = (liveDisplay.head + 1) % HISTORY_LEN;
        if (liveDisplay.count < HISTORY_LEN)
            ++liveDisplay.count;
        liveDisplay.nextFounder = static_cast<int>(model.state().nextFounder);
        liveDisplay.lastFounder = lastFounder;
        liveDisplay.founderGlow = founderGlow;
    }

    void publishDisplayFrame() {
        liveDisplay.nextFounder = static_cast<int>(model.state().nextFounder);
        liveDisplay.lastFounder = lastFounder;
        liveDisplay.founderGlow = founderGlow;
        displaySnapshot.writable() = liveDisplay;
        displaySnapshot.publish();
    }

    void publishSaveFrame() {
        SaveFrame& saved = saveSnapshot.writable();
        saved.modelState = model.state();
        saved.resetSeed = resetSeed;
        saved.resetStream = resetStream;
        saved.generationPhase = generationPhase;
        saved.rampPhase = rampPhase;
        saved.rampDuration = rampDuration;
        std::copy(rampFrom, rampFrom + CONTINUOUS_LEN, saved.rampFrom);
        std::copy(rampTo, rampTo + CONTINUOUS_LEN, saved.rampTo);
        std::copy(rampCurrent, rampCurrent + CONTINUOUS_LEN, saved.rampCurrent);
        saved.lastFounder = lastFounder;
        saved.founderGlow = founderGlow;
        saved.valid = true;
        saveSnapshot.publish();
    }

    void initializeState(bool clearEdges) {
        model.reset(resetSeed, resetStream);
        generationPhase = 0.0;
        rampPhase = 1.0;
        rampDuration = MANUAL_RAMP_TIME;
        metricsToVoltages(rampTo);
        std::copy(rampTo, rampTo + CONTINUOUS_LEN, rampFrom);
        std::copy(rampTo, rampTo + CONTINUOUS_LEN, rampCurrent);
        lastFounder = -1;
        founderGlow = 0.f;
        lossPulse.reset();
        sweepPulse.reset();
        if (clearEdges) {
            stepTrigger.reset();
            founderInputTrigger.reset();
            founderButtonTrigger.reset();
            resetInputTrigger.reset();
            resetButtonTrigger.reset();
        }
        resetHistory();
        publishDisplayFrame();
        publishSaveFrame();
    }

    void onReset(const ResetEvent& e) override {
        const SeedAction action = seedAction;
        seedAction = SEED_ACTION_NORMAL;
        const bool customSeedAction = action != SEED_ACTION_NORMAL;
        if (action == SEED_ACTION_NEW_RANDOM) {
            resetSeed = random::u64();
            resetStream = random::u64() >> 1u;
        }
        else if (action == SEED_ACTION_FACTORY) {
            resetSeed = ISLANDS_FACTORY_SEED;
            resetStream = ISLANDS_FACTORY_STREAM;
        }
        else {
            Module::onReset(e);
        }
        // A context-menu reseed behaves like the panel RESET: held STEP and
        // FOUNDER inputs must not become fresh edges on the following sample.
        initializeState(!customSeedAction);
    }

    void onRandomize(const RandomizeEvent& e) override {
        Module::onRandomize(e);
        seedAction = SEED_ACTION_NORMAL;
        resetSeed = random::u64();
        resetStream = random::u64() >> 1u;
        initializeState(true);
    }

    void advanceRamp(float sampleTime) {
        if (rampPhase < 1.0) {
            const double duration = std::max(rampDuration, 1e-9);
            rampPhase = std::min(1.0,
                rampPhase + static_cast<double>(sampleTime) / duration);
            const float phase = static_cast<float>(rampPhase);
            for (int i = 0; i < CONTINUOUS_LEN; ++i)
                rampCurrent[i] = rampFrom[i] + (rampTo[i] - rampFrom[i]) * phase;
        }
        else {
            std::copy(rampTo, rampTo + CONTINUOUS_LEN, rampCurrent);
        }
    }

    void beginTransition(float rate) {
        std::copy(rampCurrent, rampCurrent + CONTINUOUS_LEN, rampFrom);
        metricsToVoltages(rampTo);
        rampPhase = 0.0;
        rampDuration = rate > 0.f ? 1.0 / static_cast<double>(rate)
                                  : static_cast<double>(MANUAL_RAMP_TIME);
        rampDuration = std::max(1.0 / 200.0, std::min(rampDuration, 512.0));
    }

    void collectEvents() {
        const coalescent::IslandsModel::Metrics& metrics = model.metrics();
        if (metrics.lossEvent)
            lossPulse.trigger(PULSE_TIME);
        if (metrics.sweepEvent)
            sweepPulse.trigger(PULSE_TIME);
    }

    void advanceGeneration(float rate) {
        model.advance(readGenerationParameters());
        collectEvents();
        beginTransition(rate);
        pushHistory();
        publishSaveFrame();
    }

    void applyFounder(float rate) {
        lastFounder = static_cast<int>(model.state().nextFounder);
        model.founder();
        founderGlow = 1.f;
        collectEvents();
        beginTransition(rate);
        pushHistory();
        publishSaveFrame();
    }

    static void appendFloatArray(json_t* root, const char* key,
                                 const float* values, int length) {
        json_t* array = json_array();
        for (int i = 0; i < length; ++i)
            json_array_append_new(array, json_real(values[i]));
        json_object_set_new(root, key, array);
    }

    static bool readFloatArray(json_t* root, const char* key,
                               float* values, int length, float low, float high) {
        json_t* array = json_object_get(root, key);
        if (!array || !json_is_array(array) || json_array_size(array) != static_cast<size_t>(length))
            return false;
        for (int i = 0; i < length; ++i) {
            json_t* item = json_array_get(array, i);
            if (!item || !json_is_number(item))
                return false;
            float value = static_cast<float>(json_number_value(item));
            if (!std::isfinite(value) || value < low || value > high)
                return false;
            values[i] = value;
        }
        return true;
    }

    static void appendUintArray(json_t* root, const char* key,
                                const std::array<std::uint32_t, 4>& values) {
        json_t* array = json_array();
        for (int i = 0; i < 4; ++i)
            json_array_append_new(array, json_integer(values[i]));
        json_object_set_new(root, key, array);
    }

    static bool readUintArray(json_t* root, const char* key,
                              std::array<std::uint32_t, 4>& values) {
        json_t* array = json_object_get(root, key);
        if (!array || !json_is_array(array) || json_array_size(array) != 4u)
            return false;
        for (int i = 0; i < 4; ++i) {
            json_t* item = json_array_get(array, i);
            if (!item || !json_is_integer(item))
                return false;
            json_int_t value = json_integer_value(item);
            if (value < 0 || value > 4096)
                return false;
            values[i] = static_cast<std::uint32_t>(value);
        }
        return true;
    }

    static bool readHex64(json_t* root, const char* key, std::uint64_t& value) {
        json_t* item = json_object_get(root, key);
        if (!item || !json_is_string(item))
            return false;
        const char* source = json_string_value(item);
        return source && coalescent::Hex64Codec::parse(
            source, json_string_length(item), value);
    }

    static bool readBoolean(json_t* root, const char* key, bool& value) {
        json_t* item = json_object_get(root, key);
        if (!item || !json_is_boolean(item))
            return false;
        value = json_is_true(item);
        return true;
    }

    static bool readFinite(json_t* root, const char* key, float& value,
                           float low, float high) {
        json_t* item = json_object_get(root, key);
        if (!item || !json_is_number(item))
            return false;
        float parsed = static_cast<float>(json_number_value(item));
        if (!std::isfinite(parsed) || parsed < low || parsed > high)
            return false;
        value = parsed;
        return true;
    }

    static bool readFiniteDouble(json_t* root, const char* key, double& value,
                                 double low, double high) {
        json_t* item = json_object_get(root, key);
        if (!item || !json_is_number(item))
            return false;
        double parsed = json_number_value(item);
        if (!std::isfinite(parsed) || parsed < low || parsed > high)
            return false;
        value = parsed;
        return true;
    }

    json_t* dataToJson() override {
        json_t* root = json_object();
        const SaveFrame& saved = saveSnapshot.consume();
        if (!saved.valid)
            return root;

        json_object_set_new(root, "islandsVersion", json_integer(2));
        appendUintArray(root, "counts", saved.modelState.counts);
        appendUintArray(root, "denominators", saved.modelState.denominators);
        char rngState[coalescent::Hex64Codec::TEXT_SIZE];
        char rngIncrement[coalescent::Hex64Codec::TEXT_SIZE];
        char resetSeedText[coalescent::Hex64Codec::TEXT_SIZE];
        char resetStreamText[coalescent::Hex64Codec::TEXT_SIZE];
        coalescent::Hex64Codec::format(saved.modelState.rngState, rngState);
        coalescent::Hex64Codec::format(saved.modelState.rngIncrement, rngIncrement);
        coalescent::Hex64Codec::format(saved.resetSeed, resetSeedText);
        coalescent::Hex64Codec::format(saved.resetStream, resetStreamText);
        json_object_set_new(root, "rngState", json_string(rngState));
        json_object_set_new(root, "rngIncrement", json_string(rngIncrement));
        json_object_set_new(root, "resetSeed", json_string(resetSeedText));
        json_object_set_new(root, "resetStream", json_string(resetStreamText));
        json_object_set_new(root, "nextFounder", json_integer(saved.modelState.nextFounder));
        json_object_set_new(root, "wasFixA", json_boolean(saved.modelState.wasFixA));
        json_object_set_new(root, "wasFixB", json_boolean(saved.modelState.wasFixB));
        json_object_set_new(root, "sweepAArmed", json_boolean(saved.modelState.sweepAArmed));
        json_object_set_new(root, "sweepBArmed", json_boolean(saved.modelState.sweepBArmed));
        json_object_set_new(root, "generationPhase", json_real(saved.generationPhase));
        json_object_set_new(root, "rampPhase", json_real(saved.rampPhase));
        json_object_set_new(root, "rampDuration", json_real(saved.rampDuration));
        appendFloatArray(root, "rampFrom", saved.rampFrom, CONTINUOUS_LEN);
        appendFloatArray(root, "rampTo", saved.rampTo, CONTINUOUS_LEN);
        appendFloatArray(root, "rampCurrent", saved.rampCurrent, CONTINUOUS_LEN);
        json_object_set_new(root, "lastFounder", json_integer(saved.lastFounder));
        json_object_set_new(root, "founderGlow", json_real(saved.founderGlow));
        // Keep the schema-2 shape readable by the immediately preceding build,
        // but never serialize an in-flight presentation pulse.
        json_object_set_new(root, "lossPulseRemaining", json_real(0.0));
        json_object_set_new(root, "sweepPulseRemaining", json_real(0.0));
        return root;
    }

    void dataFromJson(json_t* root) override {
        json_t* version = json_object_get(root, "islandsVersion");
        if (!version || !json_is_integer(version))
            return;
        const json_int_t versionValue = json_integer_value(version);
        if (versionValue != 1 && versionValue != 2)
            return;

        std::uint64_t restoredResetSeed = ISLANDS_FACTORY_SEED;
        std::uint64_t restoredResetStream = ISLANDS_FACTORY_STREAM;
        if (versionValue >= 2
            && (!readHex64(root, "resetSeed", restoredResetSeed)
                || !readHex64(root, "resetStream", restoredResetStream)
                || restoredResetStream > (UINT64_MAX >> 1u)))
            return;

        coalescent::IslandsModel::State restored;
        if (!readUintArray(root, "counts", restored.counts)
            || !readUintArray(root, "denominators", restored.denominators)
            || !readHex64(root, "rngState", restored.rngState)
            || !readHex64(root, "rngIncrement", restored.rngIncrement))
            return;
        if (versionValue >= 2
            && restored.rngIncrement != ((restoredResetStream << 1u) | 1u))
            return;

        json_t* nextFounderItem = json_object_get(root, "nextFounder");
        if (!nextFounderItem || !json_is_integer(nextFounderItem))
            return;
        json_int_t nextFounderValue = json_integer_value(nextFounderItem);
        if (nextFounderValue < 0 || nextFounderValue >= 4)
            return;
        restored.nextFounder = static_cast<std::uint32_t>(nextFounderValue);
        if (!readBoolean(root, "wasFixA", restored.wasFixA)
            || !readBoolean(root, "wasFixB", restored.wasFixB)
            || !readBoolean(root, "sweepAArmed", restored.sweepAArmed)
            || !readBoolean(root, "sweepBArmed", restored.sweepBArmed))
            return;

        double restoredGenerationPhase;
        double restoredRampPhase;
        double restoredRampDuration;
        float restoredFrom[CONTINUOUS_LEN];
        float restoredTo[CONTINUOUS_LEN];
        float restoredCurrent[CONTINUOUS_LEN];
        if (!readFiniteDouble(root, "generationPhase", restoredGenerationPhase, 0.0, 1.0)
            || !readFiniteDouble(root, "rampPhase", restoredRampPhase, 0.0, 1.0)
            || !readFiniteDouble(root, "rampDuration", restoredRampDuration, 1.0 / 200.0, 512.0)
            || !readFloatArray(root, "rampFrom", restoredFrom, CONTINUOUS_LEN, 0.f, 10.f)
            || !readFloatArray(root, "rampTo", restoredTo, CONTINUOUS_LEN, 0.f, 10.f)
            || !readFloatArray(root, "rampCurrent", restoredCurrent, CONTINUOUS_LEN, 0.f, 10.f))
            return;

        json_t* lastFounderItem = json_object_get(root, "lastFounder");
        if (!lastFounderItem || !json_is_integer(lastFounderItem))
            return;
        json_int_t restoredLastFounder = json_integer_value(lastFounderItem);
        if (restoredLastFounder < -1 || restoredLastFounder >= 4)
            return;
        float restoredFounderGlow;
        if (!readFinite(root, "founderGlow", restoredFounderGlow, 0.f, 1.f)
            || restoredGenerationPhase >= 1.0)
            return;

        // Validate both biological and presentation state before changing the
        // running module. A corrupt patch cannot manufacture an arbitrary CV
        // trajectory unrelated to the restored population.
        coalescent::IslandsModel candidate;
        if (!candidate.restore(restored))
            return;
        const coalescent::IslandsModel::Metrics& candidateMetrics = candidate.metrics();
        float expectedTarget[CONTINUOUS_LEN];
        for (int i = 0; i < 4; ++i)
            expectedTarget[i] = static_cast<float>(candidateMetrics.frequency[i]) * 10.f;
        expectedTarget[MEAN_CONTINUOUS] = static_cast<float>(candidateMetrics.mean) * 10.f;
        expectedTarget[HET_CONTINUOUS] = static_cast<float>(candidateMetrics.diversity) * 10.f;
        const float interpolationPhase = static_cast<float>(restoredRampPhase);
        for (int i = 0; i < CONTINUOUS_LEN; ++i) {
            float expectedCurrent = restoredFrom[i]
                + (restoredTo[i] - restoredFrom[i]) * interpolationPhase;
            if (std::fabs(restoredTo[i] - expectedTarget[i]) > 1e-4f
                || std::fabs(restoredCurrent[i] - expectedCurrent) > 1e-3f)
                return;
        }

        model = candidate;
        resetSeed = restoredResetSeed;
        resetStream = restoredResetStream;
        seedAction = SEED_ACTION_NORMAL;
        generationPhase = restoredGenerationPhase;
        rampPhase = restoredRampPhase;
        rampDuration = restoredRampDuration;
        std::copy(restoredFrom, restoredFrom + CONTINUOUS_LEN, rampFrom);
        std::copy(restoredTo, restoredTo + CONTINUOUS_LEN, rampTo);
        std::copy(restoredCurrent, restoredCurrent + CONTINUOUS_LEN, rampCurrent);
        lastFounder = static_cast<int>(restoredLastFounder);
        founderGlow = restoredFounderGlow;
        // LOSS and SWEEP are one-shot presentation signals, not authored state.
        // Schema-2 patches written by older versions may still contain their
        // former remainder fields; unknown JSON keys are deliberately ignored.
        lossPulse.reset();
        sweepPulse.reset();
        stepTrigger.reset();
        founderInputTrigger.reset();
        founderButtonTrigger.reset();
        resetInputTrigger.reset();
        resetButtonTrigger.reset();
        resetHistory();
        liveDisplay.lastFounder = lastFounder;
        liveDisplay.founderGlow = founderGlow;
        publishDisplayFrame();
        publishSaveFrame();
    }

    void process(const ProcessArgs& args) override {
        if (args.sampleRate != lastSampleRate) {
            lastSampleRate = args.sampleRate;
            displayDivider.setDivision(std::max(1,
                static_cast<int>(std::round(args.sampleRate / 45.f))));
            displayDivider.reset();
            saveDivider.setDivision(std::max(1,
                static_cast<int>(std::round(args.sampleRate / 500.f))));
            saveDivider.reset();
        }

        advanceRamp(args.sampleTime);
        founderGlow = std::max(0.f, founderGlow - args.sampleTime / FOUNDER_GLOW_TIME);

        const bool stepEdge = stepTrigger.process(inputVoltage(STEP_INPUT), 0.1f, 1.f);
        const bool founderInputEdge = founderInputTrigger.process(
            inputVoltage(FOUNDER_INPUT), 0.1f, 1.f);
        const bool founderButtonEdge = founderButtonTrigger.process(
            safe(params[FOUNDER_PARAM].getValue()), 0.1f, 0.9f);
        const bool resetInputEdge = resetInputTrigger.process(
            inputVoltage(RESET_INPUT), 0.1f, 1.f);
        const bool resetButtonEdge = resetButtonTrigger.process(
            safe(params[RESET_PARAM].getValue()), 0.1f, 0.9f);
        const float rate = generationRate();

        if (resetInputEdge || resetButtonEdge) {
            // Preserve Schmitt high states so a held reset produces one edge.
            initializeState(false);
        }
        else {
            if (founderInputEdge || founderButtonEdge)
                applyFounder(rate);
            if (stepEdge)
                advanceGeneration(rate);
        }

        if (rate > 0.f) {
            generationPhase += static_cast<double>(rate) * args.sampleTime;
            int guard = 0;
            while (generationPhase >= 1.0 && guard++ < 8) {
                generationPhase -= 1.0;
                advanceGeneration(rate);
            }
            if (generationPhase >= 1.0)
                generationPhase = std::fmod(generationPhase, 1.0);
        }

        for (int i = 0; i < CONTINUOUS_LEN; ++i)
            outputs[I1_OUTPUT + i].setVoltage(rampCurrent[i]);
        const coalescent::IslandsModel::Metrics& metrics = model.metrics();
        outputs[FIX_A_OUTPUT].setVoltage(metrics.fixA ? GATE_LEVEL : 0.f);
        outputs[FIX_B_OUTPUT].setVoltage(metrics.fixB ? GATE_LEVEL : 0.f);
        outputs[LOSS_OUTPUT].setVoltage(lossPulse.process(args.sampleTime) ? GATE_LEVEL : 0.f);
        outputs[SWEEP_OUTPUT].setVoltage(sweepPulse.process(args.sampleTime) ? GATE_LEVEL : 0.f);

        if (displayDivider.process())
            publishDisplayFrame();
        if (saveDivider.process())
            publishSaveFrame();
    }
};

namespace islands_layout {
    static const float KNOB_Y = 58.f;
    static const float ATT_Y = 70.f;
    static const float INPUT_Y = 84.5f;
    static const float OUTPUT_1_Y = 101.f;
    static const float OUTPUT_2_Y = 116.5f;
    static const float KNOB_X[5] = {8.5f, 24.55f, 40.64f, 56.73f, 72.78f};
    static const float INPUT_X[8] = {6.5f, 16.25f, 26.f, 35.75f,
                                     45.5f, 55.25f, 65.f, 74.75f};
}

struct IslandsHistoryView : widget::TransparentWidget {
    Islands* module = nullptr;

    static NVGcolor laneColor(int lane, unsigned alpha = 0xff) {
        static const unsigned color[4][3] = {
            {0x68, 0xd5, 0xc0}, {0xf0, 0xb8, 0x68},
            {0x78, 0xb8, 0xec}, {0xdf, 0x83, 0xa6}
        };
        return nvgRGBA(color[lane][0], color[lane][1], color[lane][2], alpha);
    }

    void draw(const DrawArgs& args) override {
        nvgBeginPath(args.vg);
        nvgRoundedRect(args.vg, 0.f, 0.f, box.size.x, box.size.y, mm2px(1.7f));
        nvgFillColor(args.vg, nvgRGB(0x06, 0x12, 0x13));
        nvgFill(args.vg);
        nvgStrokeColor(args.vg, nvgRGB(0x35, 0x56, 0x52));
        nvgStrokeWidth(args.vg, 1.2f);
        nvgStroke(args.vg);
        Widget::draw(args);
    }

    void drawLayer(const DrawArgs& args, int layer) override {
        if (layer != 1)
            return;
        static const Islands::DisplayFrame preview;
        const Islands::DisplayFrame& frame = module ? module->displaySnapshot.consume() : preview;
        const float width = box.size.x;
        const float height = box.size.y;
        const float left = mm2px(2.6f);
        const float right = width - mm2px(2.6f);
        const float top = mm2px(7.2f);
        const float bottom = height - mm2px(3.f);
        const float plotHeight = bottom - top;
        const float pointStep = (right - left) / static_cast<float>(Islands::HISTORY_LEN - 1);
        auto yFor = [&](float value) {
            return bottom - clamp(value, 0.f, 1.f) * plotHeight;
        };

        nvgScissor(args.vg, 0.f, 0.f, width, height);
        for (int guide = 0; guide <= 2; ++guide) {
            float y = top + plotHeight * guide / 2.f;
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, left, y);
            nvgLineTo(args.vg, right, y);
            nvgStrokeColor(args.vg, nvgRGBA(0xa0, 0xc4, 0xbc, guide == 1 ? 0x28 : 0x18));
            nvgStrokeWidth(args.vg, 1.f);
            nvgStroke(args.vg);
        }
        for (int guide = 1; guide < 4; ++guide) {
            float x = left + (right - left) * guide / 4.f;
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, x, top);
            nvgLineTo(args.vg, x, bottom);
            nvgStrokeColor(args.vg, nvgRGBA(0xa0, 0xc4, 0xbc, 0x12));
            nvgStrokeWidth(args.vg, 1.f);
            nvgStroke(args.vg);
        }

        int count = clamp(frame.count, 0, Islands::HISTORY_LEN);
        int oldest = (frame.head - count + Islands::HISTORY_LEN) % Islands::HISTORY_LEN;
        for (int lane = 0; lane < 4; ++lane) {
            if (count > 1) {
                nvgBeginPath(args.vg);
                for (int j = 0; j < count; ++j) {
                    int index = (oldest + j) % Islands::HISTORY_LEN;
                    float x = right - pointStep * static_cast<float>(count - 1 - j);
                    float y = yFor(frame.history[lane][index]);
                    if (j == 0) nvgMoveTo(args.vg, x, y);
                    else        nvgLineTo(args.vg, x, y);
                }
                nvgStrokeColor(args.vg, laneColor(lane, 0xc8));
                nvgStrokeWidth(args.vg, 1.35f);
                nvgLineJoin(args.vg, NVG_ROUND);
                nvgStroke(args.vg);
            }

            float markerY = yFor(frame.current[lane]);
            if (lane == frame.lastFounder && frame.founderGlow > 0.f) {
                nvgBeginPath(args.vg);
                nvgCircle(args.vg, right, markerY, mm2px(1.9f + frame.founderGlow));
                nvgFillColor(args.vg, laneColor(lane,
                    static_cast<unsigned>(0x18 + frame.founderGlow * 0x40)));
                nvgFill(args.vg);
            }
            if (lane == frame.nextFounder) {
                nvgBeginPath(args.vg);
                nvgCircle(args.vg, right, markerY, mm2px(1.45f));
                nvgStrokeColor(args.vg, nvgRGBA(0xee, 0xf5, 0xe8, 0xc0));
                nvgStrokeWidth(args.vg, 1.f);
                nvgStroke(args.vg);
            }
            nvgBeginPath(args.vg);
            nvgCircle(args.vg, right, markerY, mm2px(0.75f));
            nvgFillColor(args.vg, laneColor(lane));
            nvgFill(args.vg);
        }
        nvgResetScissor(args.vg);

        std::shared_ptr<Font> font = APP->window->loadFont(
            asset::system("res/fonts/DejaVuSans.ttf"));
        if (font) {
            nvgFontFaceId(args.vg, font->handle);
            nvgFontSize(args.vg, mm2px(3.2f));
            nvgTextLetterSpacing(args.vg, mm2px(0.45f));
            nvgFillColor(args.vg, nvgRGBA(0xe8, 0xef, 0xdd, 0xe8));
            nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
            nvgText(args.vg, mm2px(2.6f), mm2px(2.1f), "ISLANDS", nullptr);
            nvgTextLetterSpacing(args.vg, 0.f);
            nvgFontSize(args.vg, mm2px(2.f));
            nvgFillColor(args.vg, nvgRGBA(0x8d, 0xb7, 0xae, 0xa8));
            nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_TOP);
            nvgText(args.vg, width - mm2px(2.5f), mm2px(2.75f),
                    "Wright-Fisher drift", nullptr);
        }
    }
};

struct IslandsLabels : widget::Widget {
    void draw(const DrawArgs& args) override {
        std::shared_ptr<Font> font = APP->window->loadFont(
            asset::system("res/fonts/Nunito-Bold.ttf"));
        if (!font)
            return;
        nvgFontFaceId(args.vg, font->handle);
        nvgTextLetterSpacing(args.vg, 0.f);
        nvgFillColor(args.vg, nvgRGB(0xe6, 0xf1, 0xef));
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
        auto label = [&](float x, float y, const char* text, float size = 1.48f) {
            nvgFontSize(args.vg, mm2px(size * 1.72f));
            nvgText(args.vg, mm2px(x), mm2px(y), text, nullptr);
        };
        using namespace islands_layout;
        const char* knobs[5] = {"SIZE", "SELECT", "MUTATE", "MIGRATE", "GENERATIONS"};
        for (int i = 0; i < 5; ++i)
            label(KNOB_X[i], KNOB_Y - 8.f, knobs[i], i == 4 ? 1.25f : 1.48f);
        label(KNOB_X[0], ATT_Y - 5.2f, "FOUNDER", 1.25f);
        label(KNOB_X[4], ATT_Y - 5.2f, "RESET", 1.35f);

        const char* inputs[8] = {"SIZE", "SEL", "MUT", "MIG", "GEN", "STEP", "FOUND", "RESET"};
        for (int i = 0; i < 8; ++i)
            label(INPUT_X[i], INPUT_Y - 5.5f, inputs[i], 1.22f);
        const char* topOutputs[5] = {"I1", "I2", "I3", "I4", "MEAN"};
        const char* bottomOutputs[5] = {"HET", "FIX A", "FIX B", "LOSS", "SWEEP"};
        for (int i = 0; i < 5; ++i) {
            label(KNOB_X[i], OUTPUT_1_Y - 5.5f, topOutputs[i], 1.38f);
            label(KNOB_X[i], OUTPUT_2_Y - 5.5f, bottomOutputs[i], 1.30f);
        }
    }
};

struct IslandsWidget : ModuleWidget {
    IslandsWidget(Islands* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Islands.svg")));
        addPanelLabels<IslandsLabels>(this);

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.f, 1.f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(75.28f, 1.f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.f, 122.f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(75.28f, 122.f))));

        IslandsHistoryView* display = new IslandsHistoryView();
        display->module = module;
        display->box.pos = mm2px(Vec(6.f, 8.f));
        display->box.size = mm2px(Vec(69.28f, 38.f));
        addChild(display);

        using namespace islands_layout;
        const int knobs[5] = {Islands::SIZE_PARAM, Islands::SELECT_PARAM,
                              Islands::MUTATE_PARAM, Islands::MIGRATE_PARAM,
                              Islands::GENERATIONS_PARAM};
        for (int i = 0; i < 5; ++i)
            addParam(createParamCentered<RoundBlackKnob>(
                mm2px(Vec(KNOB_X[i], KNOB_Y)), module, knobs[i]));
        addParam(createParamCentered<TL1105>(mm2px(Vec(KNOB_X[0], ATT_Y)),
                                             module, Islands::FOUNDER_PARAM));
        addParam(createParamCentered<Trimpot>(mm2px(Vec(KNOB_X[1], ATT_Y)),
                                              module, Islands::SELECT_ATT_PARAM));
        addParam(createParamCentered<Trimpot>(mm2px(Vec(KNOB_X[2], ATT_Y)),
                                              module, Islands::MUTATE_ATT_PARAM));
        addParam(createParamCentered<Trimpot>(mm2px(Vec(KNOB_X[3], ATT_Y)),
                                              module, Islands::MIGRATE_ATT_PARAM));
        addParam(createParamCentered<TL1105>(mm2px(Vec(KNOB_X[4], ATT_Y)),
                                             module, Islands::RESET_PARAM));

        const int inputIds[8] = {Islands::SIZE_INPUT, Islands::SELECT_INPUT,
                                 Islands::MUTATE_INPUT, Islands::MIGRATE_INPUT,
                                 Islands::GENERATIONS_INPUT, Islands::STEP_INPUT,
                                 Islands::FOUNDER_INPUT, Islands::RESET_INPUT};
        for (int i = 0; i < 8; ++i)
            addInput(createInputCentered<PJ301MPort>(
                mm2px(Vec(INPUT_X[i], INPUT_Y)), module, inputIds[i]));

        const int topOutputs[5] = {Islands::I1_OUTPUT, Islands::I2_OUTPUT,
                                   Islands::I3_OUTPUT, Islands::I4_OUTPUT,
                                   Islands::MEAN_OUTPUT};
        const int bottomOutputs[5] = {Islands::HET_OUTPUT, Islands::FIX_A_OUTPUT,
                                      Islands::FIX_B_OUTPUT, Islands::LOSS_OUTPUT,
                                      Islands::SWEEP_OUTPUT};
        for (int i = 0; i < 5; ++i) {
            addOutput(createOutputCentered<PJ301MPort>(
                mm2px(Vec(KNOB_X[i], OUTPUT_1_Y)), module, topOutputs[i]));
            addOutput(createOutputCentered<PJ301MPort>(
                mm2px(Vec(KNOB_X[i], OUTPUT_2_Y)), module, bottomOutputs[i]));
        }
    }

    void applySeedAction(Islands::SeedAction seedAction, const std::string& historyName) {
        Islands* module = dynamic_cast<Islands*>(this->module);
        if (!module)
            return;

        history::ModuleChange* action = new history::ModuleChange;
        action->name = historyName;
        action->moduleId = module->id;
        action->oldModuleJ = APP->engine->moduleToJson(module);
        module->seedAction = seedAction;
        APP->engine->resetModule(module);
        action->newModuleJ = APP->engine->moduleToJson(module);
        APP->history->push(action);
    }

    void appendContextMenu(Menu* menu) override {
        Islands* module = dynamic_cast<Islands*>(this->module);
        menu->addChild(new MenuSeparator);
        menu->addChild(createMenuLabel("Random sequence"));
        menu->addChild(createMenuItem("New random seed", "", [this]() {
            applySeedAction(Islands::SEED_ACTION_NEW_RANDOM, "reseed Islands");
        }, !module));
        menu->addChild(createMenuItem("Restore factory seed", "", [this]() {
            applySeedAction(Islands::SEED_ACTION_FACTORY,
                            "restore Islands factory seed");
        }, !module));
    }
};

Model* modelIslands = createModel<Islands, IslandsWidget>("Islands");
