// The harness owns an Engine to exercise the same reset path as Rack's
// context-menu Initialize action; Engine construction is private to plugins.
#include <engine/Engine.hpp>
#undef PRIVATE
#include <rack.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

using namespace rack;

Plugin* pluginInstance = nullptr;

// Compile the production Rack wrapper into this focused integration harness.
#include "../../src/Lineages.cpp"

namespace {

Module::ProcessArgs zeroTimeArgs() {
    Module::ProcessArgs args;
    args.sampleRate = 48000.f;
    args.sampleTime = 0.f;
    args.frame = 0;
    return args;
}

void fail(const char* message) {
    std::fprintf(stderr, "Lineages Rack integration: FAIL: %s\n", message);
    std::exit(1);
}

bool sameFloat(float a, float b) {
    std::uint32_t aBits = 0u;
    std::uint32_t bBits = 0u;
    std::memcpy(&aBits, &a, sizeof(aBits));
    std::memcpy(&bBits, &b, sizeof(bBits));
    return aBits == bBits;
}

std::string treeText(const Tree& tree) {
    json_t* json = treeToJson(tree);
    char* text = json_dumps(json, JSON_COMPACT | JSON_SORT_KEYS);
    json_decref(json);
    if (!text)
        fail("tree JSON dump");
    const std::string result(text);
    std::free(text);
    return result;
}

std::uint64_t fnv1a(const std::string& text) {
    std::uint64_t hash = UINT64_C(1469598103934665603);
    for (std::size_t index = 0; index < text.size(); ++index) {
        hash ^= static_cast<unsigned char>(text[index]);
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

void armTriggers(Lineages& module) {
    const Module::ProcessArgs args = zeroTimeArgs();
    for (int output = 0; output < Lineages::OUTPUTS_LEN; ++output)
        if (!module.outputs[output].isConnected())
            module.outputs[output].channels = 1;
    module.params[Lineages::NEW_PARAM].setValue(0.f);
    module.params[Lineages::RESET_PARAM].setValue(0.f);
    module.inputs[Lineages::NEW_INPUT].channels = 0;
    module.inputs[Lineages::RESET_INPUT].channels = 0;
    module.process(args);
}

void press(Lineages& module, int paramId) {
    const Module::ProcessArgs args = zeroTimeArgs();
    module.params[paramId].setValue(1.f);
    module.process(args);
    module.params[paramId].setValue(0.f);
    module.process(args);
}

struct OutputFrame {
    int traitChannels = 0;
    float traits[MAX_LEAVES] = {};
    float node = 0.f;
    float mutation = 0.f;
    float mrca = 0.f;
    float lineages = 0.f;
    float diversity = 0.f;
};

OutputFrame captureOutputs(Lineages& module) {
    OutputFrame frame;
    frame.traitChannels = module.outputs[Lineages::TRAITS_OUTPUT].getChannels();
    for (int channel = 0; channel < frame.traitChannels; ++channel)
        frame.traits[channel] =
            module.outputs[Lineages::TRAITS_OUTPUT].getVoltage(channel);
    frame.node = module.outputs[Lineages::NODE_OUTPUT].getVoltage();
    frame.mutation = module.outputs[Lineages::MUTATION_OUTPUT].getVoltage();
    frame.mrca = module.outputs[Lineages::MRCA_OUTPUT].getVoltage();
    frame.lineages = module.outputs[Lineages::LINEAGES_OUTPUT].getVoltage();
    frame.diversity = module.outputs[Lineages::DIVERSITY_OUTPUT].getVoltage();
    return frame;
}

bool sameOutputs(const OutputFrame& a, const OutputFrame& b) {
    if (a.traitChannels != b.traitChannels ||
        !sameFloat(a.node, b.node) ||
        !sameFloat(a.mutation, b.mutation) ||
        !sameFloat(a.mrca, b.mrca) ||
        !sameFloat(a.lineages, b.lineages) ||
        !sameFloat(a.diversity, b.diversity))
        return false;
    for (int channel = 0; channel < a.traitChannels; ++channel)
        if (!sameFloat(a.traits[channel], b.traits[channel]))
            return false;
    return true;
}

bool sameDisplayGeometry(const Lineages::DisplayFrame& a,
                         const Lineages::DisplayFrame& b) {
    if (a.valid != b.valid || a.sampleCount != b.sampleCount ||
        a.nodeCount != b.nodeCount || a.mutationCount != b.mutationCount ||
        a.root != b.root)
        return false;
    for (int index = 0; index < a.nodeCount; ++index) {
        if (!sameFloat(a.nodes[index].x, b.nodes[index].x) ||
            !sameFloat(a.nodes[index].time, b.nodes[index].time) ||
            a.nodes[index].parent != b.nodes[index].parent ||
            a.nodes[index].descendantMask != b.nodes[index].descendantMask)
            return false;
    }
    for (int index = 0; index < a.mutationCount; ++index) {
        if (!sameFloat(a.mutations[index].time, b.mutations[index].time) ||
            a.mutations[index].branchChild != b.mutations[index].branchChild ||
            a.mutations[index].sign != b.mutations[index].sign)
            return false;
    }
    return true;
}

void requireDisplayMatchesPlayback(const Lineages& module,
                                   const Lineages::DisplayFrame& frame) {
    if (!frame.valid || frame.sampleCount != module.tree.leafCount ||
        frame.nodeCount != module.tree.nodeCount ||
        frame.mutationCount != module.tree.mutationCount ||
        frame.root != module.tree.root || frame.generation == 0u ||
        !std::isfinite(frame.cursor) || frame.cursor < 0.f || frame.cursor > 1.f)
        fail("display frame header");
    for (int index = 0; index < frame.nodeCount; ++index) {
        const Lineages::DisplayNode& node = frame.nodes[index];
        if (!std::isfinite(node.x) || !std::isfinite(node.time) ||
            node.x < 0.f || node.x > 1.f || node.time < 0.f || node.time > 1.f ||
            !sameFloat(node.x, static_cast<float>(
                module.tree.nodes[index].displayX)) ||
            !sameFloat(node.time, static_cast<float>(
                module.tree.nodes[index].normTime)) ||
            node.parent != module.tree.nodes[index].parent ||
            node.descendantMask != module.tree.nodes[index].descendantMask)
            fail("display node geometry");
    }
    for (int index = 0; index < frame.mutationCount; ++index) {
        const Lineages::DisplayMutation& mutation = frame.mutations[index];
        if (!std::isfinite(mutation.time) || mutation.time < 0.f ||
            mutation.time > 1.f || mutation.branchChild < 0 ||
            mutation.branchChild >= frame.nodeCount ||
            (mutation.sign != -1 && mutation.sign != 1) ||
            !sameFloat(mutation.time, static_cast<float>(
                module.tree.mutations[index].normTime)) ||
            mutation.branchChild != module.tree.mutations[index].branchChild ||
            mutation.sign != module.tree.mutations[index].sign)
            fail("display mutation geometry");
    }
    const coalescent::lineages::PlaybackValues& values = module.playback.values();
    if (frame.activeCount != values.activeLineageCount ||
        frame.activeCount != values.lineageCount ||
        frame.ancestry != (module.playback.direction() == Direction::Ancestry) ||
        !sameFloat(frame.cursor, static_cast<float>(module.playback.cursor())))
        fail("display transport state");
    for (int index = 0; index < frame.activeCount; ++index) {
        if (frame.activeMasks[index] != values.activeLineageMasks[index])
            fail("display active lineage mask");
    }
}

void requireFiniteOutputs(Lineages& module) {
    const OutputFrame frame = captureOutputs(module);
    for (int channel = 0; channel < frame.traitChannels; ++channel)
        if (!std::isfinite(frame.traits[channel]))
            fail("non-finite TRAITS output");
    if (!std::isfinite(frame.node) || !std::isfinite(frame.mutation) ||
        !std::isfinite(frame.mrca) || !std::isfinite(frame.lineages) ||
        !std::isfinite(frame.diversity))
        fail("non-finite monophonic output");
}

void testNewResetAndDirtyState() {
    Lineages module;
    armTriggers(module);
    if (!module.tree.validate() || module.tree.leafCount != 8 ||
        module.outputs[Lineages::TRAITS_OUTPUT].getChannels() != 8)
        fail("default tree/output contract");

    const std::string original = treeText(module.tree);
    const coalescent::lineages::Pcg32::State rngBeforeReset =
        module.generator.rngState();
    press(module, Lineages::RESET_PARAM);
    if (treeText(module.tree) != original ||
        module.generator.rngState().state != rngBeforeReset.state ||
        module.generator.rngState().increment != rngBeforeReset.increment)
        fail("panel RESET changed tree or RNG");

    module.params[Lineages::SAMPLES_PARAM].setValue(13.f);
    module.process(zeroTimeArgs());
    if (!module.isDirty() ||
        module.outputs[Lineages::TRAITS_OUTPUT].getChannels() != 8 ||
        module.tree.leafCount != 8)
        fail("dirty SAMPLES regenerated or changed polyphony");
    module.params[Lineages::SAMPLES_PARAM].setValue(8.f);
    module.process(zeroTimeArgs());
    if (module.isDirty())
        fail("returning SAMPLES did not clear dirty state");

    const float generatedMutate = static_cast<float>(
        mutateFromKey(module.generatedMutateKey));
    module.params[Lineages::MUTATE_PARAM].setValue(0.91f);
    module.process(zeroTimeArgs());
    if (!module.isDirty())
        fail("MUTATE did not mark tree dirty");
    module.params[Lineages::MUTATE_PARAM].setValue(generatedMutate);
    module.process(zeroTimeArgs());
    if (module.isDirty())
        fail("returning MUTATE did not clear dirty state");

    module.params[Lineages::SAMPLES_PARAM].setValue(11.f);
    module.params[Lineages::MUTATE_PARAM].setValue(0.7f);
    press(module, Lineages::NEW_PARAM);
    if (module.isDirty() || module.tree.leafCount != 11 ||
        module.outputs[Lineages::TRAITS_OUTPUT].getChannels() != 11 ||
        treeText(module.tree) == original)
        fail("NEW did not install matching replacement tree");

    std::printf("Rack lifecycle: NEW, panel RESET, dirty return PASS\n");
}

void testContextInitialize() {
    Lineages module;
    armTriggers(module);
    module.params[Lineages::RATE_PARAM].setValue(4.f);
    module.params[Lineages::SAMPLES_PARAM].setValue(14.f);
    module.params[Lineages::MUTATE_PARAM].setValue(0.9f);
    module.params[Lineages::STEP_PARAM].setValue(11.f);
    module.params[Lineages::DIRECTION_PARAM].setValue(0.f);
    module.params[Lineages::LOOP_PARAM].setValue(0.f);
    press(module, Lineages::NEW_PARAM);

    const std::string beforeTree = treeText(module.tree);
    const coalescent::lineages::Pcg32::State beforeRng =
        module.generator.rngState();
    const std::uint64_t beforeGeneration = module.displayGeneration;

    Context context;
    contextSet(&context);
    engine::Engine engine;
    context.engine = &engine;
    engine.addModule(&module);
    engine.resetModule(&module);
    engine.removeModule(&module);
    context.engine = nullptr;
    contextSet(nullptr);

    if (module.params[Lineages::RATE_PARAM].getValue() != 0.f ||
        module.params[Lineages::SAMPLES_PARAM].getValue() != 8.f ||
        !sameFloat(module.params[Lineages::MUTATE_PARAM].getValue(), 0.4f) ||
        module.params[Lineages::STEP_PARAM].getValue() != 2.f ||
        module.params[Lineages::DIRECTION_PARAM].getValue() != 1.f ||
        module.params[Lineages::LOOP_PARAM].getValue() != 1.f)
        fail("Rack Initialize did not restore default parameters");
    if (!module.tree.validate() || module.tree.leafCount != 8 ||
        module.generatedMutateKey != mutateKeyFromValue(0.4f) ||
        treeText(module.tree) == beforeTree)
        fail("Rack Initialize did not install a fresh default tree");
    const coalescent::lineages::Pcg32::State afterRng =
        module.generator.rngState();
    if (afterRng.state == beforeRng.state &&
        afterRng.increment == beforeRng.increment)
        fail("Rack Initialize did not advance the local RNG");
    if (module.playback.direction() != Direction::Descent ||
        !module.playback.loop() ||
        module.playback.cursor() != 1.0)
        fail("Rack Initialize did not reset default transport");
    if (module.displayGeneration != beforeGeneration + 1u)
        fail("Rack Initialize display generation");

    std::printf("Rack context Initialize: defaults, fresh tree, transport PASS\n");
}

void testDisplaySnapshot() {
    Lineages module;
    armTriggers(module);
    module.publishDisplayFrame();
    Lineages::DisplayFrame frame = module.displaySnapshot.consume();
    requireDisplayMatchesPlayback(module, frame);
    const std::uint64_t initialGeneration = frame.generation;
    const Lineages::DisplayFrame initialGeometry = frame;
    if (frame.activeCount != 1 ||
        frame.activeMasks[0] != module.tree.nodes[module.tree.root].descendantMask)
        fail("display MRCA active lineage");

    press(module, Lineages::RESET_PARAM);
    frame = module.displaySnapshot.consume();
    if (frame.generation != initialGeneration ||
        !sameDisplayGeometry(frame, initialGeometry))
        fail("panel RESET changed display generation or geometry");

    module.params[Lineages::SAMPLES_PARAM].setValue(12.f);
    module.process(zeroTimeArgs());
    module.publishDisplayFrame();
    frame = module.displaySnapshot.consume();
    if (!frame.dirty || frame.sampleCount != 8 ||
        frame.generation != initialGeneration)
        fail("display dirty state");
    module.params[Lineages::SAMPLES_PARAM].setValue(8.f);
    module.process(zeroTimeArgs());
    module.publishDisplayFrame();
    frame = module.displaySnapshot.consume();
    if (frame.dirty || frame.generation != initialGeneration)
        fail("display dirty state did not clear");

    module.params[Lineages::DIRECTION_PARAM].setValue(0.f);
    module.process(zeroTimeArgs());
    frame = module.displaySnapshot.consume();
    if (!frame.ancestry || frame.generation != initialGeneration ||
        !sameDisplayGeometry(frame, initialGeometry))
        fail("direction changed display geometry or generation");

    module.inputs[Lineages::TIME_INPUT].channels = 1;
    module.inputs[Lineages::TIME_INPUT].setVoltage(0.f);
    module.process(zeroTimeArgs());
    frame = module.displaySnapshot.consume();
    requireDisplayMatchesPlayback(module, frame);
    if (frame.activeCount != frame.sampleCount)
        fail("display present active lineage count");
    for (int index = 0; index < frame.activeCount; ++index) {
        const std::uint16_t mask = frame.activeMasks[index];
        if (mask == 0u || (mask & static_cast<std::uint16_t>(mask - 1u)) != 0u)
            fail("display present active mask was not a leaf");
    }

    module.inputs[Lineages::TIME_INPUT].setVoltage(4.25f);
    module.process(zeroTimeArgs());
    module.publishDisplayFrame();
    frame = module.displaySnapshot.consume();
    requireDisplayMatchesPlayback(module, frame);
    if (!sameFloat(frame.cursor, 0.425f))
        fail("display cursor did not follow TIME");

    press(module, Lineages::NEW_PARAM);
    frame = module.displaySnapshot.consume();
    requireDisplayMatchesPlayback(module, frame);
    if (frame.generation != initialGeneration + 1u || frame.dirty)
        fail("NEW display generation or dirty state");

    const Lineages::DisplayFrame preview = Lineages::previewDisplayFrame();
    if (!preview.valid || preview.sampleCount != 8 || preview.nodeCount != 15 ||
        preview.activeCount < 1 || preview.generation != 1u)
        fail("static display preview");

    std::printf("Rack display: geometry, cursor, active masks, dirty, generation PASS\n");
}

void testDirectionTimeLoopAndOutputs() {
    Lineages module;
    armTriggers(module);
    module.params[Lineages::MUTATE_PARAM].setValue(1.f);
    press(module, Lineages::NEW_PARAM);

    module.inputs[Lineages::TIME_INPUT].channels = 1;
    module.inputs[Lineages::TIME_INPUT].setVoltage(10.f);
    module.resetPulses();
    module.process(zeroTimeArgs());
    if (module.nodePulse.remaining != 0.f ||
        module.mutationPulse.remaining != 0.f ||
        module.mrcaPulse.remaining != 0.f)
        fail("TIME connection emitted historical pulses");
    const int samples = module.tree.leafCount;
    const float rootVoltage =
        module.outputs[Lineages::TRAITS_OUTPUT].getVoltage(0);
    for (int channel = 1; channel < samples; ++channel)
        if (!sameFloat(rootVoltage,
                       module.outputs[Lineages::TRAITS_OUTPUT].getVoltage(channel)))
            fail("TRAITS did not converge exactly at MRCA");

    module.inputs[Lineages::TIME_INPUT].setVoltage(0.f);
    module.process(zeroTimeArgs());
    if (module.nodePulse.remaining > PULSE_SECONDS ||
        module.mutationPulse.remaining > PULSE_SECONDS ||
        module.mrcaPulse.remaining > PULSE_SECONDS)
        fail("large TIME jump emitted an invalid pulse aggregate");
    if (module.outputs[Lineages::LINEAGES_OUTPUT].getVoltage() != 10.f)
        fail("present lineage voltage");

    module.resetPulses();
    module.inputs[Lineages::TIME_INPUT].channels = 0;
    const double beforeDisconnect = module.playback.cursor();
    module.process(zeroTimeArgs());
    if (module.playback.cursor() != beforeDisconnect ||
        module.nodePulse.remaining != 0.f ||
        module.mutationPulse.remaining != 0.f ||
        module.mrcaPulse.remaining != 0.f)
        fail("TIME disconnect moved transport or emitted events");

    module.params[Lineages::DIRECTION_PARAM].setValue(0.f);
    module.params[Lineages::LOOP_PARAM].setValue(0.f);
    press(module, Lineages::RESET_PARAM);
    if (module.playback.cursor() != 0.0)
        fail("ancestry RESET source");
    Module::ProcessArgs movement = zeroTimeArgs();
    movement.sampleTime = 0.1f;
    module.process(movement);
    const double afterAncestry = module.playback.cursor();
    if (!(afterAncestry > 0.0))
        fail("ancestry movement");
    module.params[Lineages::DIRECTION_PARAM].setValue(1.f);
    module.process(movement);
    if (module.playback.cursor() != afterAncestry)
        fail("direction reversal moved on switch sample");
    module.process(movement);
    if (!(module.playback.cursor() < afterAncestry))
        fail("descent did not reverse from current cursor");

    module.params[Lineages::LOOP_PARAM].setValue(1.f);
    module.params[Lineages::RATE_PARAM].setValue(5.f);
    module.process(zeroTimeArgs());
    press(module, Lineages::RESET_PARAM);
    movement.sampleTime = 1.f;
    module.process(movement);
    if (module.playback.cursor() != 1.0 || !module.playback.running())
        fail("descent loop did not teleport to source");

    module.params[Lineages::STEP_PARAM].setValue(
        std::numeric_limits<float>::quiet_NaN());
    module.params[Lineages::RATE_PARAM].setValue(
        std::numeric_limits<float>::infinity());
    module.inputs[Lineages::RATE_INPUT].channels = 1;
    module.inputs[Lineages::RATE_INPUT].setVoltage(
        -std::numeric_limits<float>::infinity());
    module.process(zeroTimeArgs());
    requireFiniteOutputs(module);

    std::printf("Rack transport: ancestry/descent, reversal, TIME, loop, finite PASS\n");
}

void testExactRoundTripAndNextNew() {
    Lineages module;
    armTriggers(module);
    module.params[Lineages::MUTATE_PARAM].setValue(1.f);
    module.params[Lineages::STEP_PARAM].setValue(5.f);
    press(module, Lineages::NEW_PARAM);

    module.inputs[Lineages::TIME_INPUT].channels = 1;
    module.inputs[Lineages::TIME_INPUT].setVoltage(10.f);
    module.process(zeroTimeArgs());
    module.inputs[Lineages::TIME_INPUT].setVoltage(0.f);
    module.process(zeroTimeArgs());
    module.publishSaveFrame();

    const std::string savedTree = treeText(module.tree);
    const OutputFrame savedOutputs = captureOutputs(module);
    json_t* patchState = module.dataToJson();
    if (!patchState)
        fail("dataToJson returned null");

    press(module, Lineages::NEW_PARAM);
    const std::string replacementTree = treeText(module.tree);
    if (replacementTree == savedTree)
        fail("replacement tree unexpectedly identical");

    const std::uint64_t generationBeforeRestore = module.displayGeneration;
    module.dataFromJson(patchState);
    module.process(zeroTimeArgs());
    if (treeText(module.tree) != savedTree ||
        !sameOutputs(captureOutputs(module), savedOutputs))
        fail("saved tree or outputs did not restore exactly");
    const Lineages::DisplayFrame restoredDisplay =
        module.displaySnapshot.consume();
    requireDisplayMatchesPlayback(module, restoredDisplay);
    if (module.displayGeneration != generationBeforeRestore + 1u ||
        restoredDisplay.generation != module.displayGeneration)
        fail("restored tree display generation");

    armTriggers(module);
    press(module, Lineages::NEW_PARAM);
    const std::string reproducedTree = treeText(module.tree);
    if (reproducedTree != replacementTree)
        fail("post-restore NEW did not reproduce exact next tree");

    json_t* nextRng = json_object_get(patchState, "nextRng");
    json_t* state = nextRng ? json_object_get(nextRng, "state") : nullptr;
    json_t* increment = nextRng ? json_object_get(nextRng, "increment") : nullptr;
    if (!state || !increment || !json_is_string(state) ||
        !json_is_string(increment) || json_string_length(state) != 16u ||
        json_string_length(increment) != 16u)
        fail("RNG was not serialized as fixed-width hexadecimal strings");

    std::printf(
        "Rack persistence: saved=%016llx replacement=%016llx reproduced=%016llx PASS\n",
        static_cast<unsigned long long>(fnv1a(savedTree)),
        static_cast<unsigned long long>(fnv1a(replacementTree)),
        static_cast<unsigned long long>(fnv1a(reproducedTree)));
    json_decref(patchState);
}

} // namespace

int main() {
    rack::random::init();
    testNewResetAndDirtyState();
    testContextInitialize();
    testDisplaySnapshot();
    testDirectionTimeLoopAndOutputs();
    testExactRoundTripAndNextNew();
    std::printf("Lineages Rack integration: PASS\n");
    return 0;
}
