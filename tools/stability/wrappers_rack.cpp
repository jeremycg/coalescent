// Rack integration coverage for the production wrappers whose model mathematics
// are already exercised by the SDK-free shared-core contracts.
#include <engine/Engine.hpp>
#undef PRIVATE
#include <rack.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

using namespace rack;

Plugin* pluginInstance = nullptr;

// A few older wrappers use generic file-local constant names. Prefix them while
// compiling the production sources into this aggregate harness so their internal
// linkage names do not collide in one translation unit.
#define MAX_N GENDYN_MAX_N
#include "../../src/GENDYN.cpp"
#undef MAX_N

#define MAX_N HAPTIK_MAX_N
#include "../../src/Haptik.cpp"
#undef MAX_N

#define TRAIL AXON_TRAIL
#define ORBIT_N AXON_ORBIT_N
#define ORBIT_SECS AXON_ORBIT_SECS
#include "../../src/neuron/Axon.cpp"
#undef ORBIT_SECS
#undef ORBIT_N
#undef TRAIL

#define TRAIL SOMA_TRAIL
#define ORBIT_N SOMA_ORBIT_N
#define ORBIT_SECS SOMA_ORBIT_SECS
#include "../../src/neuron/Soma.cpp"
#undef ORBIT_SECS
#undef ORBIT_N
#undef TRAIL

#include "../../src/Operon.cpp"
#include "../../src/Bunnies.cpp"
#include "../../src/Foxes.cpp"
#include "../../src/Archipelago.cpp"

namespace {

static constexpr float SAMPLE_RATE = 48000.f;
static constexpr float SAMPLE_TIME = 1.f / SAMPLE_RATE;
static_assert(sizeof(float) == sizeof(std::uint32_t), "Rack requires 32-bit float");
static_assert(std::numeric_limits<float>::is_iec559, "finite check requires IEEE-754 float");

Module::ProcessArgs processArgs(std::int64_t frame = 0) {
    Module::ProcessArgs args;
    args.sampleRate = SAMPLE_RATE;
    args.sampleTime = SAMPLE_TIME;
    args.frame = frame;
    return args;
}

[[noreturn]] void fail(const char* module, const char* message) {
    std::fprintf(stderr, "%s Rack integration: FAIL: %s\n", module, message);
    std::exit(1);
}

void require(bool condition, const char* module, const char* message) {
    if (!condition)
        fail(module, message);
}

bool sameFloat(float a, float b) {
    std::uint32_t aBits = 0u;
    std::uint32_t bBits = 0u;
    std::memcpy(&aBits, &a, sizeof(aBits));
    std::memcpy(&bBits, &b, sizeof(bBits));
    return aBits == bBits;
}

// Keep the harness assertion valid under Rack's production
// -funsafe-math-optimizations flags. A library isfinite() call can otherwise be
// weakened by compiler assumptions about non-finite values.
bool finiteFloat(float value) {
    std::uint32_t bits = 0u;
    std::memcpy(&bits, &value, sizeof(bits));
    return (bits & UINT32_C(0x7f800000)) != UINT32_C(0x7f800000);
}

bool sameDouble(double a, double b) {
    std::uint64_t aBits = 0u;
    std::uint64_t bBits = 0u;
    std::memcpy(&aBits, &a, sizeof(aBits));
    std::memcpy(&bBits, &b, sizeof(bBits));
    return aBits == bBits;
}

template <typename ModuleType>
void requireFiniteOutputs(ModuleType& module, int outputCount, const char* name) {
    for (int outputId = 0; outputId < outputCount; ++outputId) {
        const int channels = std::max(1, module.outputs[outputId].getChannels());
        for (int channel = 0; channel < channels; ++channel) {
            if (!finiteFloat(module.outputs[outputId].getVoltage(channel)))
                fail(name, "non-finite output");
        }
    }
}

template <typename ModuleType>
void processFiniteSamples(ModuleType& module, int outputCount, const char* name,
                          int samples, std::int64_t firstFrame = 0) {
    for (int i = 0; i < samples; ++i) {
        module.process(processArgs(firstFrame + i));
        requireFiniteOutputs(module, outputCount, name);
    }
}

template <typename ModuleType>
void connectOutputs(ModuleType& module, int outputCount) {
    for (int outputId = 0; outputId < outputCount; ++outputId)
        module.outputs[outputId].channels = 1;
}

template <typename ModuleType>
void connectHostileInputs(ModuleType& module, int inputCount) {
    static const float values[4] = {
        std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::infinity(),
        -std::numeric_limits<float>::infinity(),
        1e30f
    };
    for (int inputId = 0; inputId < inputCount; ++inputId) {
        module.inputs[inputId].channels = 4;
        for (int channel = 0; channel < 4; ++channel)
            module.inputs[inputId].setVoltage(values[(inputId + channel) % 4], channel);
    }
}

template <typename ModuleType>
void initializeThroughRack(ModuleType& module) {
    Context context;
    contextSet(&context);
    engine::Engine engine;
    context.engine = &engine;
    engine.addModule(&module);
    engine.resetModule(&module);
    engine.removeModule(&module);
    context.engine = nullptr;
    contextSet(nullptr);
}

void testGendyn() {
    const char* name = "GENDYN";
    GENDYN module;
    connectOutputs(module, GENDYN::OUTPUTS_LEN);
    module.params[GENDYN::N_PARAM].setValue(8.f);
    module.params[GENDYN::LOCK_PARAM].setValue(1.f);
    module.initShape = 2;
    module.process(processArgs());
    require(module.last_N == 8 && module.lastInitShape == 2, name,
            "construction did not install the requested initial waveform");
    requireFiniteOutputs(module, GENDYN::OUTPUTS_LEN, name);

    std::array<float, 8> expectedAmp;
    std::array<float, 8> expectedDur;
    std::array<float, 8> expectedStepAmp;
    std::array<float, 8> expectedStepDur;
    for (int i = 0; i < 8; ++i) {
        expectedAmp[i] = module.amp[i];
        expectedDur[i] = module.dur[i];
        expectedStepAmp[i] = module.step_amp[i];
        expectedStepDur[i] = module.step_dur[i];
    }
    module.lockCorr = 1.125f;
    module.dur_err = 0.25f;
    module.publishSaveFrame(8, SAMPLE_RATE);
    json_t* saved = module.dataToJson();
    require(saved != nullptr, name, "dataToJson returned null");
    require(json_integer_value(json_object_get(saved, "N")) == 8
            && json_array_size(json_object_get(saved, "amp")) == 8u
            && json_array_size(json_object_get(saved, "dur")) == 8u
            && json_is_number(json_object_get(saved, "lockCorr"))
            && json_is_number(json_object_get(saved, "durErr"))
            && json_is_number(json_object_get(saved, "saveFs")),
            name, "persistence schema is incomplete");

    std::fill(module.amp, module.amp + 8, -0.75f);
    std::fill(module.dur, module.dur + 8, 1.f);
    std::fill(module.step_amp, module.step_amp + 8, 0.5f);
    std::fill(module.step_dur, module.step_dur + 8, -0.5f);
    module.dataFromJson(saved);
    require(module.restoredPending && module.last_N == 8 && module.initShape == 2,
            name, "valid persisted state was not accepted");
    for (int i = 0; i < 8; ++i) {
        require(sameFloat(module.amp[i], expectedAmp[i])
                && sameFloat(module.dur[i], expectedDur[i])
                && sameFloat(module.step_amp[i], expectedStepAmp[i])
                && sameFloat(module.step_dur[i], expectedStepDur[i]),
                name, "waveform or walk state did not restore exactly");
    }
    require(sameFloat(module.lockCorr, 1.125f)
            && sameFloat(module.restoredDurErr, 0.25f)
            && sameFloat(module.restoredFs, SAMPLE_RATE),
            name, "LOCK controller state did not restore exactly");
    module.process(processArgs(1));
    require(!module.restoredPending, name, "restored state was not installed on process");

    bool sawCycle = false;
    for (int frame = 2; frame < 4000; ++frame) {
        module.process(processArgs(frame));
        requireFiniteOutputs(module, GENDYN::OUTPUTS_LEN, name);
        sawCycle = sawCycle
            || module.outputs[GENDYN::CYCLE_TRIG_OUTPUT].getVoltage() > 1.f;
    }
    require(sawCycle, name, "cycle event output never fired");
    requireFiniteOutputs(module, GENDYN::OUTPUTS_LEN, name);
    json_decref(saved);

    module.params[GENDYN::N_PARAM].setValue(3.f);
    module.initShape = 4;
    initializeThroughRack(module);
    require(module.params[GENDYN::N_PARAM].getValue() == 13.f
            && module.initShape == 0 && module.reseedPending.load(),
            name, "Rack Initialize did not restore defaults and request a reseed");
    module.process(processArgs());
    require(module.last_N == 13, name, "Rack Initialize did not install the default waveform");

    module.params[GENDYN::SCALE_AMP_ATT_PARAM].setValue(1.f);
    module.params[GENDYN::SCALE_DUR_ATT_PARAM].setValue(1.f);
    module.params[GENDYN::B_AMP_ATT_PARAM].setValue(1.f);
    module.params[GENDYN::B_DUR_ATT_PARAM].setValue(1.f);
    connectHostileInputs(module, GENDYN::INPUTS_LEN);
    processFiniteSamples(module, GENDYN::OUTPUTS_LEN, name, 64);
    std::printf("GENDYN Rack persistence/lifecycle/events/hostile CV PASS\n");
}

void testHaptik() {
    const char* name = "Haptik";
    Haptik module;
    connectOutputs(module, Haptik::OUTPUTS_LEN);
    module.params[Haptik::N_PARAM].setValue(16.f);
    module.params[Haptik::EXCITE_PARAM].setValue(0.f);
    module.process(processArgs());
    require(module.last_N == 16, name, "construction did not seed the requested lattice");
    requireFiniteOutputs(module, Haptik::OUTPUTS_LEN, name);

    std::array<float, 16> expectedY;
    std::array<float, 16> expectedV;
    for (int i = 0; i < 16; ++i) {
        expectedY[i] = (i - 8) * 0.015625f;
        expectedV[i] = (7 - i) * 0.0078125f;
        module.y[i] = expectedY[i];
        module.v[i] = expectedV[i];
    }
    module.scanPhase = 0.375f;
    module.publishSaveFrame(16);
    json_t* saved = module.dataToJson();
    require(saved != nullptr, name, "dataToJson returned null");
    require(json_integer_value(json_object_get(saved, "N")) == 16
            && json_array_size(json_object_get(saved, "y")) == 16u
            && json_array_size(json_object_get(saved, "v")) == 16u
            && json_is_number(json_object_get(saved, "scanPhase")),
            name, "persistence schema is incomplete");

    std::fill(module.y, module.y + 16, 1.f);
    std::fill(module.v, module.v + 16, -1.f);
    module.scanPhase = 0.f;
    module.dataFromJson(saved);
    require(module.last_N == 16 && sameFloat(module.scanPhase, 0.375f)
            && module.divCounter == 0 && !module.wasFrozen,
            name, "valid persisted lattice was not accepted");
    for (int i = 0; i < 16; ++i) {
        require(sameFloat(module.y[i], expectedY[i])
                && sameFloat(module.v[i], expectedV[i])
                && sameFloat(module.yPrev[i], expectedY[i]),
                name, "lattice state did not restore exactly");
    }
    json_decref(saved);
    module.process(processArgs(1));
    requireFiniteOutputs(module, Haptik::OUTPUTS_LEN, name);

    module.params[Haptik::N_PARAM].setValue(128.f);
    module.params[Haptik::MODE_PARAM].setValue(1.f);
    initializeThroughRack(module);
    require(module.params[Haptik::N_PARAM].getValue() == 64.f
            && module.params[Haptik::MODE_PARAM].getValue() == 0.f
            && module.last_N == -1,
            name, "Rack Initialize did not restore defaults and clear the lattice");
    module.process(processArgs());
    require(module.last_N == 64, name, "Rack Initialize did not reseed the default lattice");

    // Prove the wrapper mapping itself, not merely the downstream state clamp:
    // a non-finite EXT sum is exactly neutral, while an ordinary finite drive
    // changes the same deterministic lattice.
    Haptik neutralDrive;
    Haptik hostileDrive;
    Haptik finiteDrive;
    Haptik* driveModules[3] = {&neutralDrive, &hostileDrive, &finiteDrive};
    for (int i = 0; i < 3; ++i) {
        driveModules[i]->params[Haptik::N_PARAM].setValue(16.f);
        driveModules[i]->params[Haptik::EXCITE_PARAM].setValue(0.f);
        driveModules[i]->process(processArgs());
        driveModules[i]->inputs[Haptik::EXT_INPUT].channels = 1;
    }
    neutralDrive.inputs[Haptik::EXT_INPUT].setVoltage(0.f);
    hostileDrive.inputs[Haptik::EXT_INPUT].setVoltage(
        std::numeric_limits<float>::quiet_NaN());
    finiteDrive.inputs[Haptik::EXT_INPUT].setVoltage(5.f);
    neutralDrive.process(processArgs(1));
    hostileDrive.process(processArgs(1));
    finiteDrive.process(processArgs(1));
    bool hostileMatchesNeutral = true;
    bool finiteDriveDiffers = false;
    for (int i = 0; i < 16; ++i) {
        hostileMatchesNeutral = hostileMatchesNeutral
            && sameFloat(hostileDrive.y[i], neutralDrive.y[i])
            && sameFloat(hostileDrive.v[i], neutralDrive.v[i]);
        finiteDriveDiffers = finiteDriveDiffers
            || !sameFloat(finiteDrive.y[i], neutralDrive.y[i])
            || !sameFloat(finiteDrive.v[i], neutralDrive.v[i]);
    }
    require(hostileMatchesNeutral && finiteDriveDiffers, name,
            "EXT non-finite fallback was not exactly neutral or finite drive had no effect");

    module.params[Haptik::RATE_ATT_PARAM].setValue(1.f);
    module.params[Haptik::COUPLE_ATT_PARAM].setValue(1.f);
    module.params[Haptik::DAMP_ATT_PARAM].setValue(1.f);
    module.params[Haptik::INJECT_ATT_PARAM].setValue(1.f);
    connectHostileInputs(module, Haptik::INPUTS_LEN);
    processFiniteSamples(module, Haptik::OUTPUTS_LEN, name, 64);
    std::printf("Haptik Rack persistence/lifecycle/hostile CV PASS\n");
}

template <typename Neuron>
void testNeuron(Neuron& module, const char* name,
                int pitchParam, int currentParam,
                int currentAttParam, int secondAttParam,
                int voctInput, int trigInput,
                int output, int spikeOutput, int outputCount, int inputCount,
                float currentDefault) {
    connectOutputs(module, outputCount);
    module.process(processArgs());
    require(module.lastOs == 4 && module.channels == 1, name,
            "factory-default x4 process path was not installed");
    requireFiniteOutputs(module, outputCount, name);

    module.oversampleMode = 0;
    module.params[pitchParam].setValue(1.f);
    module.inputs[voctInput].channels = 5;
    static const float offsets[5] = {-1.f, -0.5f, 0.f, 0.5f, 1.f};
    for (int channel = 0; channel < 5; ++channel)
        module.inputs[voctInput].setVoltage(offsets[channel], channel);

    bool sawSpike[5] = {};
    for (int frame = 0; frame < 12000; ++frame) {
        module.process(processArgs(frame));
        requireFiniteOutputs(module, outputCount, name);
        for (int channel = 0; channel < 5; ++channel)
            sawSpike[channel] = sawSpike[channel]
                || module.outputs[spikeOutput].getVoltage(channel) > 1.f;
    }
    require(module.outputs[output].getChannels() == 5
            && module.outputs[spikeOutput].getChannels() == 5,
            name, "V/OCT polyphony did not map to outputs");
    for (int channel = 0; channel < 5; ++channel)
        require(sawSpike[channel], name, "a polyphonic voice never emitted its event gate");
    requireFiniteOutputs(module, outputCount, name);

    module.inputs[voctInput].channels = 0;
    module.inputs[trigInput].channels = 3;
    for (int channel = 0; channel < 3; ++channel)
        module.inputs[trigInput].setVoltage(0.f, channel);
    module.process(processArgs(12000));
    requireFiniteOutputs(module, outputCount, name);
    require(module.outputs[output].getChannels() == 3
            && module.outputs[spikeOutput].getChannels() == 3,
            name, "TRIG channel fallback did not map to outputs");

    module.oversampleMode = 3;
    json_t* saved = module.dataToJson();
    require(saved != nullptr
            && json_integer_value(json_object_get(saved, "osFactor")) == 8,
            name, "anti-alias factor was not serialized losslessly");
    module.oversampleMode = 0;
    module.dataFromJson(saved);
    require(module.oversampleMode == 3, name, "anti-alias factor did not restore");
    json_decref(saved);

    module.params[currentParam].setValue(currentDefault + 0.1f);
    initializeThroughRack(module);
    require(module.params[pitchParam].getValue() == 0.f
            && sameFloat(module.params[currentParam].getValue(), currentDefault)
            && module.oversampleMode == 2 && module.channels == 1,
            name, "Rack Initialize did not restore defaults and voice state");
    module.process(processArgs());
    require(module.lastOs == 4, name,
            "Rack Initialize did not execute the default x4 process path");
    requireFiniteOutputs(module, outputCount, name);

    module.oversampleMode = 0;
    module.params[currentAttParam].setValue(1.f);
    module.params[secondAttParam].setValue(1.f);
    connectHostileInputs(module, inputCount);
    processFiniteSamples(module, outputCount, name, 64);
    std::printf("%s Rack polyphony/events/JSON/lifecycle/hostile CV PASS\n", name);
}

void testAxonCurrentFallback() {
    const char* name = "Axon";
    Axon neutral;
    Axon hostile;
    Axon driven;
    Axon* modules[3] = {&neutral, &hostile, &driven};
    for (int i = 0; i < 3; ++i) {
        connectOutputs(*modules[i], Axon::OUTPUTS_LEN);
        modules[i]->oversampleMode = 0;
        modules[i]->params[Axon::CURRENT_ATT_PARAM].setValue(1.f);
        modules[i]->inputs[Axon::VOCT_INPUT].channels = 4;
        modules[i]->inputs[Axon::CURRENT_INPUT].channels = 4;
        for (int lane = 0; lane < 4; ++lane) {
            modules[i]->inputs[Axon::VOCT_INPUT].setVoltage(0.f, lane);
            modules[i]->inputs[Axon::CURRENT_INPUT].setVoltage(0.f, lane);
        }
    }
    hostile.inputs[Axon::CURRENT_INPUT].setVoltage(
        std::numeric_limits<float>::quiet_NaN(), 1);
    hostile.inputs[Axon::CURRENT_INPUT].setVoltage(
        std::numeric_limits<float>::infinity(), 2);
    hostile.inputs[Axon::CURRENT_INPUT].setVoltage(
        -std::numeric_limits<float>::infinity(), 3);
    driven.inputs[Axon::CURRENT_INPUT].setVoltage(1.f, 1);
    neutral.process(processArgs());
    hostile.process(processArgs());
    driven.process(processArgs());

    bool hostileMatchesNeutral = true;
    for (int lane = 0; lane < 4; ++lane) {
        hostileMatchesNeutral = hostileMatchesNeutral
            && sameFloat(hostile.vv4[0][lane], neutral.vv4[0][lane])
            && sameFloat(hostile.ww4[0][lane], neutral.ww4[0][lane]);
    }
    const bool finiteDriveDiffers =
        !sameFloat(driven.vv4[0][1], neutral.vv4[0][1])
        || !sameFloat(driven.ww4[0][1], neutral.ww4[0][1]);
    require(hostileMatchesNeutral && finiteDriveDiffers, name,
            "CURRENT non-finite lanes were not neutral or finite CV had no effect");
}

void testSomaCurrentFallback() {
    const char* name = "Soma";
    Soma neutral;
    Soma hostile;
    Soma driven;
    Soma* modules[3] = {&neutral, &hostile, &driven};
    for (int i = 0; i < 3; ++i) {
        connectOutputs(*modules[i], Soma::OUTPUTS_LEN);
        modules[i]->oversampleMode = 0;
        modules[i]->params[Soma::CURRENT_ATT_PARAM].setValue(1.f);
        modules[i]->inputs[Soma::VOCT_INPUT].channels = 4;
        modules[i]->inputs[Soma::CURRENT_INPUT].channels = 4;
        for (int lane = 0; lane < 4; ++lane) {
            modules[i]->inputs[Soma::VOCT_INPUT].setVoltage(0.f, lane);
            modules[i]->inputs[Soma::CURRENT_INPUT].setVoltage(0.f, lane);
        }
    }
    hostile.inputs[Soma::CURRENT_INPUT].setVoltage(
        std::numeric_limits<float>::quiet_NaN(), 1);
    hostile.inputs[Soma::CURRENT_INPUT].setVoltage(
        std::numeric_limits<float>::infinity(), 2);
    hostile.inputs[Soma::CURRENT_INPUT].setVoltage(
        -std::numeric_limits<float>::infinity(), 3);
    driven.inputs[Soma::CURRENT_INPUT].setVoltage(1.f, 1);
    neutral.process(processArgs());
    hostile.process(processArgs());
    driven.process(processArgs());

    bool hostileMatchesNeutral = true;
    for (int lane = 0; lane < 4; ++lane) {
        hostileMatchesNeutral = hostileMatchesNeutral
            && sameFloat(hostile.xx4[0][lane], neutral.xx4[0][lane])
            && sameFloat(hostile.yy4[0][lane], neutral.yy4[0][lane])
            && sameFloat(hostile.zz4[0][lane], neutral.zz4[0][lane]);
    }
    const bool finiteDriveDiffers =
        !sameFloat(driven.xx4[0][1], neutral.xx4[0][1])
        || !sameFloat(driven.yy4[0][1], neutral.yy4[0][1])
        || !sameFloat(driven.zz4[0][1], neutral.zz4[0][1]);
    require(hostileMatchesNeutral && finiteDriveDiffers, name,
            "CURRENT non-finite lanes were not neutral or finite CV had no effect");
}

void testNeurons() {
    Axon axon;
    testNeuron(axon, "Axon",
               Axon::PITCH_PARAM, Axon::CURRENT_PARAM,
               Axon::CURRENT_ATT_PARAM, Axon::EPS_ATT_PARAM,
               Axon::VOCT_INPUT, Axon::TRIG_INPUT,
               Axon::OUT_OUTPUT, Axon::SPIKE_OUTPUT,
               Axon::OUTPUTS_LEN, Axon::INPUTS_LEN, 0.6f);

    Soma soma;
    testNeuron(soma, "Soma",
               Soma::PITCH_PARAM, Soma::CURRENT_PARAM,
               Soma::CURRENT_ATT_PARAM, Soma::BURST_ATT_PARAM,
               Soma::VOCT_INPUT, Soma::TRIG_INPUT,
               Soma::OUT_OUTPUT, Soma::SPIKE_OUTPUT,
               Soma::OUTPUTS_LEN, Soma::INPUTS_LEN, 2.f);

    testAxonCurrentFallback();
    testSomaCurrentFallback();
}

template <typename ModuleType>
void captureEvents(ModuleType& module, int firstOutput, int eventCount,
                   int outputCount, const char* name, int samples, bool* seen) {
    for (int frame = 0; frame < samples; ++frame) {
        module.process(processArgs(frame));
        requireFiniteOutputs(module, outputCount, name);
        for (int event = 0; event < eventCount; ++event)
            seen[event] = seen[event]
                || module.outputs[firstOutput + event].getVoltage() > 1.f;
    }
}

void testOperon() {
    const char* name = "Operon";
    Operon module;
    connectOutputs(module, Operon::OUTPUTS_LEN);
    bool gates[3] = {};
    captureEvents(module, Operon::GATE1_OUTPUT, 3, Operon::OUTPUTS_LEN,
                  name, 12000, gates);
    for (int i = 0; i < 3; ++i)
        require(gates[i], name, "a phase gate never fired");
    requireFiniteOutputs(module, Operon::OUTPUTS_LEN, name);

    module.params[Operon::ALPHA_PARAM].setValue(30.f);
    initializeThroughRack(module);
    require(module.params[Operon::ALPHA_PARAM].getValue() == 12.f
            && sameFloat(module.m[0], 0.2f) && sameFloat(module.p[1], 0.4f),
            name, "Rack Initialize did not restore defaults and seed");
    module.params[Operon::ALPHA_ATT_PARAM].setValue(1.f);
    module.params[Operon::HILL_ATT_PARAM].setValue(1.f);
    module.params[Operon::BETA_ATT_PARAM].setValue(1.f);
    connectHostileInputs(module, Operon::INPUTS_LEN);
    processFiniteSamples(module, Operon::OUTPUTS_LEN, name, 64);
    std::printf("Operon Rack events/lifecycle/hostile CV PASS\n");
}

void testBunnies() {
    const char* name = "Bunnies";
    Bunnies module;
    connectOutputs(module, Bunnies::OUTPUTS_LEN);
    bool events[2] = {};
    captureEvents(module, Bunnies::PREY_POP_OUTPUT, 2, Bunnies::OUTPUTS_LEN,
                  name, 12000, events);
    for (int i = 0; i < 2; ++i)
        require(events[i], name, "a population peak event never fired");
    requireFiniteOutputs(module, Bunnies::OUTPUTS_LEN, name);

    module.params[Bunnies::WILD_PARAM].setValue(1.f);
    initializeThroughRack(module);
    require(sameFloat(module.params[Bunnies::WILD_PARAM].getValue(), 0.4f)
            && sameFloat(module.x, 1.3f) && sameFloat(module.y, 0.9f),
            name, "Rack Initialize did not restore defaults and seed");
    module.params[Bunnies::BALANCE_ATT_PARAM].setValue(1.f);
    module.params[Bunnies::WILD_ATT_PARAM].setValue(1.f);
    connectHostileInputs(module, Bunnies::INPUTS_LEN);
    processFiniteSamples(module, Bunnies::OUTPUTS_LEN, name, 64);
    std::printf("Bunnies Rack events/lifecycle/hostile CV PASS\n");
}

void testFoxes() {
    const char* name = "Foxes";
    Foxes module;
    connectOutputs(module, Foxes::OUTPUTS_LEN);
    bool events[3] = {};
    captureEvents(module, Foxes::GRASS_PEAK_OUTPUT, 3, Foxes::OUTPUTS_LEN,
                  name, 16000, events);
    for (int i = 0; i < 3; ++i)
        require(events[i], name, "a trophic peak event never fired");
    requireFiniteOutputs(module, Foxes::OUTPUTS_LEN, name);

    module.params[Foxes::WILD_PARAM].setValue(1.f);
    initializeThroughRack(module);
    require(sameFloat(module.params[Foxes::WILD_PARAM].getValue(), 0.5f)
            && sameFloat(module.x, coalescent::foxes::SEED_X)
            && sameFloat(module.y, coalescent::foxes::SEED_Y)
            && sameFloat(module.z, coalescent::foxes::SEED_Z),
            name, "Rack Initialize did not restore defaults and seed");
    module.params[Foxes::BALANCE_ATT_PARAM].setValue(1.f);
    module.params[Foxes::WILD_ATT_PARAM].setValue(1.f);
    connectHostileInputs(module, Foxes::INPUTS_LEN);
    processFiniteSamples(module, Foxes::OUTPUTS_LEN, name, 64);
    std::printf("Foxes Rack events/lifecycle/hostile CV PASS\n");
}

bool sameArchipelagoState(const coalescent::ArchipelagoField::State& a,
                          const coalescent::ArchipelagoField::State& b) {
    if (a.occupiedMask != b.occupiedMask
        || !sameDouble(a.reportedGlobalMean, b.reportedGlobalMean))
        return false;
    for (int i = 0; i < coalescent::ArchipelagoField::kValues; ++i)
        if (!sameDouble(a.density[i], b.density[i]))
            return false;
    for (int i = 0; i < coalescent::ArchipelagoField::kHabitats; ++i)
        if (!sameDouble(a.reportedTrait[i], b.reportedTrait[i]))
            return false;
    return true;
}

void testArchipelago() {
    const char* name = "Archipelago";
    Archipelago module;
    connectOutputs(module, Archipelago::OUTPUTS_LEN);
    module.process(processArgs());
    require(module.outputs[Archipelago::TRAIT_OUTPUT].getChannels()
                == Archipelago::HABITATS
            && module.outputs[Archipelago::MASS_OUTPUT].getChannels()
                == Archipelago::HABITATS,
            name, "polyphonic outputs do not expose all habitats");
    requireFiniteOutputs(module, Archipelago::OUTPUTS_LEN, name);

    const coalescent::ArchipelagoField::Parameters parameters = module.readParameters();
    module.field.advance(0.4, parameters);
    module.publishSaveFrame();
    const coalescent::ArchipelagoField::State savedState = module.field.state();
    json_t* saved = module.dataToJson();
    require(saved != nullptr
            && json_integer_value(json_object_get(saved, "archipelagoVersion")) == 1
            && json_array_size(json_object_get(saved, "density"))
                == static_cast<std::size_t>(Archipelago::VALUES),
            name, "persistence schema is incomplete");
    module.field.reset(parameters);
    module.publishSaveFrame();
    module.dataFromJson(saved);
    require(sameArchipelagoState(module.field.state(), savedState), name,
            "field state did not round-trip exactly");
    json_decref(saved);

    // Align a panel RESET edge with the 500 Hz field tick. The field must equal a
    // fresh deterministic reset exactly; advancing on this same tick would make
    // even the first bin differ.
    module.field.advance(0.4, parameters);
    coalescent::ArchipelagoField expected;
    expected.reset(parameters);
    module.resetButtonTrigger.reset();
    module.params[Archipelago::RESET_PARAM].setValue(0.f);
    module.fieldDivider.reset();
    const int division = static_cast<int>(module.fieldDivider.getDivision());
    for (int i = 0; i < division - 1; ++i) {
        module.process(processArgs(i));
        requireFiniteOutputs(module, Archipelago::OUTPUTS_LEN, name);
    }
    module.params[Archipelago::RESET_PARAM].setValue(1.f);
    module.process(processArgs(division - 1));
    requireFiniteOutputs(module, Archipelago::OUTPUTS_LEN, name);
    require(sameArchipelagoState(module.field.state(), expected.state()), name,
            "RESET did not install the deterministic field or evolved on its reset tick");
    require(module.outputs[Archipelago::COLONIZE_OUTPUT].getVoltage() == 0.f
            && module.outputs[Archipelago::EXTINCT_OUTPUT].getVoltage() == 0.f,
            name, "RESET emitted an occupancy event");
    module.params[Archipelago::RESET_PARAM].setValue(0.f);

    module.params[Archipelago::GRADIENT_PARAM].setValue(-0.4f);
    initializeThroughRack(module);
    coalescent::ArchipelagoField initializedExpected;
    initializedExpected.reset(module.readParameters());
    require(sameFloat(module.params[Archipelago::GRADIENT_PARAM].getValue(), 0.55f)
            && sameArchipelagoState(module.field.state(), initializedExpected.state()),
            name, "Rack Initialize did not restore defaults and deterministic field");

    module.params[Archipelago::SELECT_ATT_PARAM].setValue(1.f);
    module.params[Archipelago::MUTATE_ATT_PARAM].setValue(1.f);
    module.params[Archipelago::MIGRATE_ATT_PARAM].setValue(1.f);
    module.params[Archipelago::GRADIENT_ATT_PARAM].setValue(1.f);
    module.params[Archipelago::BARRIER_ATT_PARAM].setValue(1.f);
    connectHostileInputs(module, Archipelago::INPUTS_LEN);
    processFiniteSamples(module, Archipelago::OUTPUTS_LEN, name,
        2 * static_cast<int>(module.fieldDivider.getDivision()));
    std::printf("Archipelago Rack persistence/reset-tick/lifecycle/hostile CV PASS\n");
}

} // namespace

int main() {
    rack::random::init();
    testGendyn();
    testHaptik();
    testNeurons();
    testOperon();
    testBunnies();
    testFoxes();
    testArchipelago();
    std::printf("Wrapper Rack integration: all 8 additional modules PASS\n");
    return 0;
}
