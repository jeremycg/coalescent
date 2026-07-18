// Exercise the production Rack wrapper, including Rack's Initialize path.
#include <engine/Engine.hpp>
#undef PRIVATE
#include <rack.hpp>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

using namespace rack;

Plugin* pluginInstance = nullptr;

#include "../../src/Finches.cpp"

namespace {

Module::ProcessArgs zeroTimeArgs() {
    Module::ProcessArgs args;
    args.sampleRate = 48000.f;
    args.sampleTime = 0.f;
    args.frame = 0;
    return args;
}

void fail(const char* message) {
    std::fprintf(stderr, "Finches Rack integration: FAIL: %s\n", message);
    std::exit(1);
}

bool sameFloat(float a, float b) {
    std::uint32_t aBits = 0u;
    std::uint32_t bBits = 0u;
    std::memcpy(&aBits, &a, sizeof(aBits));
    std::memcpy(&bBits, &b, sizeof(bBits));
    return aBits == bBits;
}

bool sameState(const coalescent::FinchesField::State& a,
               const coalescent::FinchesField::State& b) {
    if (a.version != b.version || a.split != b.split
        || !sameFloat(a.splitTimer, b.splitTimer)
        || !sameFloat(a.mergeTimer, b.mergeTimer))
        return false;
    for (int i = 0; i < Finches::BINS; ++i)
        if (!sameFloat(a.mass[i], b.mass[i]))
            return false;
    return true;
}

void loadRoundTrip(Finches& module,
                   const coalescent::FinchesField::State& expected) {
    if (!module.field.restore(expected))
        fail("test state was rejected");
    module.publishSaveFrame();
    json_t* state = module.dataToJson();
    if (!state)
        fail("dataToJson returned null");

    module.field.reset(0.4f);
    module.publishSaveFrame();
    module.dataFromJson(state);
    if (!sameState(module.field.state(), expected))
        fail("complete detector state did not round-trip exactly");

    json_decref(state);
}

void testCompleteAndLegacyPersistence() {
    Finches module;

    coalescent::FinchesField::State pendingSplit = module.field.state();
    pendingSplit.split = false;
    pendingSplit.splitTimer = 0.125f;
    pendingSplit.mergeTimer = 0.f;
    loadRoundTrip(module, pendingSplit);

    coalescent::FinchesField::State pendingMerge = module.field.state();
    pendingMerge.split = true;
    pendingMerge.splitTimer = 0.f;
    pendingMerge.mergeTimer = 0.075f;
    loadRoundTrip(module, pendingMerge);

    json_t* legacy = json_object();
    json_t* density = json_array();
    for (int i = 0; i < Finches::BINS; ++i) {
        const float mass = (i == 15 || i == 48) ? 0.5f : 0.f;
        json_array_append_new(density, json_real(mass));
    }
    json_object_set_new(legacy, "density", density);
    module.dataFromJson(legacy);
    const coalescent::FinchesField::State restoredLegacy = module.field.state();
    if (!restoredLegacy.split || restoredLegacy.splitTimer != 0.f
        || restoredLegacy.mergeTimer != 0.f)
        fail("density-only state did not infer its historical split latch");
    json_decref(legacy);

    std::printf("Finches persistence: latch, timers, legacy density PASS\n");
}

void testContextInitializeAndFiniteOutputs() {
    Finches module;
    module.params[Finches::RATE_PARAM].setValue(3.f);
    module.params[Finches::MUTATE_PARAM].setValue(0.9f);
    module.params[Finches::COMPETE_PARAM].setValue(0.95f);
    module.params[Finches::NICHE_PARAM].setValue(0.1f);
    module.field.seed(0.7f, 0.2f);

    Context context;
    contextSet(&context);
    engine::Engine engine;
    context.engine = &engine;
    engine.addModule(&module);
    engine.resetModule(&module);
    engine.removeModule(&module);
    context.engine = nullptr;
    contextSet(nullptr);

    const coalescent::FinchesField::State initialized = module.field.state();
    if (module.params[Finches::RATE_PARAM].getValue() != 0.f
        || !sameFloat(module.params[Finches::MUTATE_PARAM].getValue(), 0.35f)
        || !sameFloat(module.params[Finches::COMPETE_PARAM].getValue(), 0.30f)
        || !sameFloat(module.params[Finches::NICHE_PARAM].getValue(), 0.48f)
        || initialized.split || initialized.splitTimer != 0.f
        || initialized.mergeTimer != 0.f)
        fail("Rack Initialize did not restore defaults and detector state");

    module.params[Finches::RATE_PARAM].setValue(
        std::numeric_limits<float>::quiet_NaN());
    module.inputs[Finches::RATE_INPUT].channels = 1;
    module.inputs[Finches::RATE_INPUT].setVoltage(
        std::numeric_limits<float>::infinity());
    module.inputs[Finches::ENV_INPUT].channels = 1;
    module.inputs[Finches::ENV_INPUT].setVoltage(
        -std::numeric_limits<float>::infinity());
    module.process(zeroTimeArgs());
    for (int output = 0; output < Finches::OUTPUTS_LEN; ++output)
        if (!std::isfinite(module.outputs[output].getVoltage()))
            fail("non-finite output under hostile controls");

    std::printf("Finches Rack Initialize + finite outputs PASS\n");
}

void testResetTickConvention() {
    Finches module;
    Module::ProcessArgs args = zeroTimeArgs();
    args.sampleRate = 500.f;
    args.sampleTime = 1.f / args.sampleRate;

    module.params[Finches::RATE_PARAM].setValue(4.f);
    module.inputs[Finches::ENV_INPUT].channels = 1;
    module.inputs[Finches::ENV_INPUT].setVoltage(2.f);
    module.inputs[Finches::RESET_INPUT].channels = 1;
    module.inputs[Finches::RESET_INPUT].setVoltage(0.f);
    module.inputs[Finches::SEED_INPUT].channels = 1;
    module.inputs[Finches::SEED_INPUT].setVoltage(0.f);
    module.process(args); // Prime both Schmitt triggers at the 500 Hz field cadence.

    module.field.seed(-0.7f, 0.2f);
    module.splitPulse.trigger(Finches::GATE_TIME);
    module.mergePulse.trigger(Finches::GATE_TIME);
    coalescent::FinchesField expected;
    expected.reset(module.environmentFromInput());

    module.inputs[Finches::RESET_INPUT].setVoltage(10.f);
    module.inputs[Finches::SEED_INPUT].setVoltage(10.f);
    module.process(args);
    if (!sameState(module.field.state(), expected.state()))
        fail("RESET tick advanced the field or admitted a simultaneous SEED");
    if (module.splitPulse.remaining != 0.f || module.mergePulse.remaining != 0.f)
        fail("RESET tick did not clear live event pulses");

    module.inputs[Finches::RESET_INPUT].setVoltage(0.f);
    module.inputs[Finches::SEED_INPUT].setVoltage(0.f);
    module.process(args);
    if (sameState(module.field.state(), expected.state()))
        fail("field did not resume evolution after the RESET tick");

    std::printf("Finches RESET tick: install now, evolve next tick PASS\n");
}

} // namespace

int main() {
    rack::random::init();
    testCompleteAndLegacyPersistence();
    testContextInitializeAndFiniteOutputs();
    testResetTickConvention();
    std::printf("Finches Rack integration: PASS\n");
    return 0;
}
