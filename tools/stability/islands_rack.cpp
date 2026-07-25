// Exercise the production Rack wrapper, including Rack's Initialize path.
#include <engine/Engine.hpp>
#undef PRIVATE
#include <rack.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

using namespace rack;

Plugin* pluginInstance = nullptr;

#include "../../src/Islands.cpp"

namespace {

Module::ProcessArgs zeroTimeArgs() {
    Module::ProcessArgs args;
    args.sampleRate = 48000.f;
    args.sampleTime = 0.f;
    args.frame = 0;
    return args;
}

void fail(const char* message) {
    std::fprintf(stderr, "Islands Rack integration: FAIL: %s\n", message);
    std::exit(1);
}

bool sameFloat(float a, float b) {
    std::uint32_t aBits = 0u;
    std::uint32_t bBits = 0u;
    std::memcpy(&aBits, &a, sizeof(aBits));
    std::memcpy(&bBits, &b, sizeof(bBits));
    return aBits == bBits;
}

bool sameState(const coalescent::IslandsModel::State& a,
               const coalescent::IslandsModel::State& b) {
    return a.counts == b.counts
        && a.denominators == b.denominators
        && a.rngState == b.rngState
        && a.rngIncrement == b.rngIncrement
        && a.nextFounder == b.nextFounder
        && a.wasFixA == b.wasFixA
        && a.wasFixB == b.wasFixB
        && a.sweepAArmed == b.sweepAArmed
        && a.sweepBArmed == b.sweepBArmed;
}

struct OutputFrame {
    float values[Islands::OUTPUTS_LEN] = {};
};

OutputFrame outputs(Islands& module) {
    OutputFrame result;
    for (int i = 0; i < Islands::OUTPUTS_LEN; ++i)
        result.values[i] = module.outputs[i].getVoltage();
    return result;
}

bool sameOutputs(const OutputFrame& a, const OutputFrame& b) {
    for (int i = 0; i < Islands::OUTPUTS_LEN; ++i)
        if (!sameFloat(a.values[i], b.values[i]))
            return false;
    return true;
}

void press(Islands& module, int paramId) {
    const Module::ProcessArgs args = zeroTimeArgs();
    module.params[paramId].setValue(1.f);
    module.process(args);
    module.params[paramId].setValue(0.f);
    module.process(args);
}

void stopClockAndArm(Islands& module) {
    module.params[Islands::GENERATIONS_PARAM].setValue(
        Islands::GENERATION_OCT_MIN);
    module.params[Islands::FOUNDER_PARAM].setValue(0.f);
    module.params[Islands::RESET_PARAM].setValue(0.f);
    module.inputs[Islands::STEP_INPUT].channels = 0;
    module.inputs[Islands::FOUNDER_INPUT].channels = 0;
    module.inputs[Islands::RESET_INPUT].channels = 0;
    module.process(zeroTimeArgs());
}

void testExactRoundTripAndNextGeneration() {
    Islands module;
    stopClockAndArm(module);
    module.params[Islands::SELECT_PARAM].setValue(0.13f);
    module.params[Islands::MUTATE_PARAM].setValue(0.62f);
    module.params[Islands::MIGRATE_PARAM].setValue(0.27f);
    for (int i = 0; i < 7; ++i)
        press(module, Islands::FOUNDER_PARAM);

    // Use STEP input for generations; panel FOUNDER above also exercises its
    // persisted round-robin index without relying on an automatic clock.
    for (int i = 0; i < 11; ++i) {
        module.inputs[Islands::STEP_INPUT].channels = 1;
        module.inputs[Islands::STEP_INPUT].setVoltage(10.f);
        module.process(zeroTimeArgs());
        module.inputs[Islands::STEP_INPUT].setVoltage(0.f);
        module.process(zeroTimeArgs());
    }
    module.publishSaveFrame();
    module.process(zeroTimeArgs());

    const coalescent::IslandsModel::State savedModel = module.model.state();
    const OutputFrame savedOutputs = outputs(module);
    json_t* saved = module.dataToJson();
    if (!saved)
        fail("dataToJson returned null");

    module.inputs[Islands::STEP_INPUT].setVoltage(10.f);
    module.process(zeroTimeArgs());
    module.inputs[Islands::STEP_INPUT].setVoltage(0.f);
    module.process(zeroTimeArgs());
    const coalescent::IslandsModel::State nextModel = module.model.state();
    if (sameState(nextModel, savedModel))
        fail("STEP did not advance the stochastic sequence");

    module.dataFromJson(saved);
    module.process(zeroTimeArgs());
    if (!sameState(module.model.state(), savedModel)
        || !sameOutputs(outputs(module), savedOutputs))
        fail("saved model or outputs did not restore exactly");

    module.inputs[Islands::STEP_INPUT].setVoltage(10.f);
    module.process(zeroTimeArgs());
    module.inputs[Islands::STEP_INPUT].setVoltage(0.f);
    module.process(zeroTimeArgs());
    if (!sameState(module.model.state(), nextModel))
        fail("post-restore generation did not reproduce exactly");

    json_t* rngState = json_object_get(saved, "rngState");
    json_t* rngIncrement = json_object_get(saved, "rngIncrement");
    if (!rngState || !rngIncrement || !json_is_string(rngState)
        || !json_is_string(rngIncrement)
        || json_string_length(rngState) != 16u
        || json_string_length(rngIncrement) != 16u)
        fail("RNG was not serialized as fixed-width hexadecimal");
    json_decref(saved);

    std::printf(
        "Islands persistence: saved RNG=%016llx next/reproduced=%016llx PASS\n",
        static_cast<unsigned long long>(savedModel.rngState),
        static_cast<unsigned long long>(nextModel.rngState));
}

void testTransientPulsesDoNotRoundTrip() {
    Islands module;
    stopClockAndArm(module);
    const coalescent::IslandsModel::State savedModel = module.model.state();

    module.lossPulse.trigger(Islands::PULSE_TIME);
    module.sweepPulse.trigger(Islands::PULSE_TIME);
    module.publishSaveFrame();
    json_t* saved = module.dataToJson();
    if (!saved)
        fail("pulse persistence fixture did not serialize");
    json_t* savedLossPulse = json_object_get(saved, "lossPulseRemaining");
    json_t* savedSweepPulse = json_object_get(saved, "sweepPulseRemaining");
    if (!json_is_number(savedLossPulse) || !json_is_number(savedSweepPulse)
        || json_number_value(savedLossPulse) != 0.0
        || json_number_value(savedSweepPulse) != 0.0)
        fail("schema-2 compatibility placeholders were absent or nonzero");

    // Loading the just-saved state while both source pulses are still live must
    // clear them rather than replaying either edge.
    module.dataFromJson(saved);
    if (!sameState(module.model.state(), savedModel)
        || module.lossPulse.remaining != 0.f
        || module.sweepPulse.remaining != 0.f)
        fail("live event pulse survived a save/load round trip");
    module.process(zeroTimeArgs());
    if (module.outputs[Islands::LOSS_OUTPUT].getVoltage() != 0.f
        || module.outputs[Islands::SWEEP_OUTPUT].getVoltage() != 0.f)
        fail("save/load replayed an event pulse at the outputs");

    // Older schema-2 JSON carried these fields. They remain accepted but have no
    // authority over the new module's transient pulse generators.
    json_object_set_new(saved, "lossPulseRemaining", json_real(Islands::PULSE_TIME));
    json_object_set_new(saved, "sweepPulseRemaining", json_real(Islands::PULSE_TIME));
    module.lossPulse.trigger(Islands::PULSE_TIME);
    module.sweepPulse.trigger(Islands::PULSE_TIME);
    module.dataFromJson(saved);
    if (!sameState(module.model.state(), savedModel)
        || module.lossPulse.remaining != 0.f
        || module.sweepPulse.remaining != 0.f)
        fail("legacy pulse fields were replayed or rejected the saved state");

    json_decref(saved);
    std::printf("Islands transient pulse persistence + legacy fields PASS\n");
}

void testSeedsResetInitializeAndFiniteOutputs() {
    Islands module;
    stopClockAndArm(module);
    const coalescent::IslandsModel::State factory = module.model.state();
    module.inputs[Islands::STEP_INPUT].channels = 1;
    module.inputs[Islands::STEP_INPUT].setVoltage(10.f);
    module.process(zeroTimeArgs());
    module.inputs[Islands::STEP_INPUT].setVoltage(0.f);
    module.process(zeroTimeArgs());
    press(module, Islands::RESET_PARAM);
    if (!sameState(module.model.state(), factory))
        fail("panel RESET did not replay the current seed");

    Context context;
    contextSet(&context);
    engine::Engine engine;
    context.engine = &engine;
    engine.addModule(&module);
    engine.randomizeModule(&module);
    const coalescent::IslandsModel::State randomized = module.model.state();
    if (randomized.rngState == factory.rngState
        && randomized.rngIncrement == factory.rngIncrement)
        fail("Rack Randomize did not install independent randomness");

    stopClockAndArm(module);
    module.inputs[Islands::STEP_INPUT].channels = 1;
    module.inputs[Islands::STEP_INPUT].setVoltage(10.f);
    module.process(zeroTimeArgs());
    module.inputs[Islands::STEP_INPUT].setVoltage(0.f);
    module.process(zeroTimeArgs());
    press(module, Islands::RESET_PARAM);
    if (!sameState(module.model.state(), randomized))
        fail("panel RESET did not replay the randomized seed");

    module.params[Islands::SIZE_PARAM].setValue(10.f);
    module.params[Islands::SELECT_PARAM].setValue(0.2f);
    engine.resetModule(&module);
    if (!sameFloat(module.params[Islands::SIZE_PARAM].getValue(), 6.f)
        || module.params[Islands::SELECT_PARAM].getValue() != 0.f
        || !sameState(module.model.state(), randomized))
        fail("Rack Initialize did not restore defaults and replay its seed");

    module.seedAction = Islands::SEED_ACTION_FACTORY;
    engine.resetModule(&module);
    if (!sameState(module.model.state(), factory))
        fail("factory seed action did not restore factory sequence");
    engine.removeModule(&module);
    context.engine = nullptr;
    contextSet(nullptr);

    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float inf = std::numeric_limits<float>::infinity();
    for (int param = 0; param < Islands::PARAMS_LEN; ++param)
        module.params[param].setValue(nan);
    for (int input = 0; input < Islands::INPUTS_LEN; ++input) {
        module.inputs[input].channels = 1;
        module.inputs[input].setVoltage((input & 1) ? inf : -inf);
    }
    module.process(zeroTimeArgs());
    for (int output = 0; output < Islands::OUTPUTS_LEN; ++output)
        if (!std::isfinite(module.outputs[output].getVoltage()))
            fail("non-finite output under hostile controls");

    std::printf("Islands seed/reset/Initialize + finite outputs PASS\n");
}

} // namespace

int main() {
    rack::random::init();
    testExactRoundTripAndNextGeneration();
    testTransientPulsesDoNotRoundTrip();
    testSeedsResetInitializeAndFiniteOutputs();
    std::printf("Islands Rack integration: PASS\n");
    return 0;
}
