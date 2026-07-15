#include "plugin.hpp"
#include "dsp/display_snapshot.hpp"
#include "lineages/state.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

// LINEAGES - a persistent Kingman genealogy played in either time direction.
//
// The SDK-free generator and playback engine own all biological and event-
// crossing rules. This file is intentionally limited to Rack transport,
// voltage mapping, persistence, and the genealogy display.

namespace {

using coalescent::lineages::Direction;
using coalescent::lineages::EventFacts;
using coalescent::lineages::Hex64Codec;
using coalescent::lineages::KingmanGenerator;
using coalescent::lineages::MAX_LEAVES;
using coalescent::lineages::MAX_MUTATIONS;
using coalescent::lineages::MAX_NODES;
using coalescent::lineages::Mutation;
using coalescent::lineages::Node;
using coalescent::lineages::Playback;
using coalescent::lineages::Snapshot;
using coalescent::lineages::Tree;

constexpr int MUTATE_QUANTIZATION = 65535;
constexpr float PULSE_SECONDS = 1e-3f;
constexpr float PULSE_VOLTS = 10.f;
constexpr float DIVERSITY_VOLTS_PER_MUTATION = 2.f;
constexpr float RATE_BASE_HZ = 0.1f;
constexpr float RATE_TOTAL_MIN = -12.f;
constexpr float RATE_TOTAL_MAX = 5.f;
constexpr float DISPLAY_HZ = 45.f;
constexpr float SAVE_HZ = 500.f;
constexpr float TRAIT_LIMIT = 10.f;

float finiteOr(float value, float fallback = 0.f) {
    return std::isfinite(value) ? value : fallback;
}

float clampFloat(float value, float low, float high) {
    if (value < low)
        return low;
    if (value > high)
        return high;
    return value;
}

double clampDouble(double value, double low, double high) {
    if (value < low)
        return low;
    if (value > high)
        return high;
    return value;
}

bool nearlyEqual(double a, double b) {
    if (!std::isfinite(a) || !std::isfinite(b))
        return false;
    const double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
    return std::fabs(a - b) <= 1e-12 * scale;
}

int mutateKeyFromValue(float value) {
    const float unit = clampFloat(finiteOr(value, 0.4f), 0.f, 1.f);
    return static_cast<int>(std::lround(unit * static_cast<float>(MUTATE_QUANTIZATION)));
}

double mutateFromKey(int key) {
    if (key < 0)
        key = 0;
    if (key > MUTATE_QUANTIZATION)
        key = MUTATE_QUANTIZATION;
    return static_cast<double>(key) / static_cast<double>(MUTATE_QUANTIZATION);
}

bool readInteger(json_t* object, const char* key, int low, int high,
                 int& destination) {
    json_t* item = json_object_get(object, key);
    if (!item || !json_is_integer(item))
        return false;
    const json_int_t value = json_integer_value(item);
    if (value < static_cast<json_int_t>(low) ||
        value > static_cast<json_int_t>(high))
        return false;
    destination = static_cast<int>(value);
    return true;
}

bool readBoolean(json_t* object, const char* key, bool& destination) {
    json_t* item = json_object_get(object, key);
    if (!item || !json_is_boolean(item))
        return false;
    destination = json_is_true(item);
    return true;
}

bool readFinite(json_t* object, const char* key, double& destination) {
    json_t* item = json_object_get(object, key);
    if (!item || !json_is_number(item))
        return false;
    const double value = json_number_value(item);
    if (!std::isfinite(value))
        return false;
    destination = value;
    return true;
}

bool readFiniteRange(json_t* object, const char* key, double low, double high,
                     double& destination) {
    double value = 0.0;
    if (!readFinite(object, key, value) || value < low || value > high)
        return false;
    destination = value;
    return true;
}

bool readHex(json_t* object, const char* key, std::uint64_t& destination) {
    json_t* item = json_object_get(object, key);
    if (!item || !json_is_string(item))
        return false;
    const char* text = json_string_value(item);
    return text && Hex64Codec::parse(text, json_string_length(item), destination);
}

void setHex(json_t* object, const char* key, std::uint64_t value) {
    char text[Hex64Codec::TEXT_SIZE];
    Hex64Codec::format(value, text);
    json_object_set_new(object, key, json_string(text));
}

json_t* nodeToJson(const Node& node) {
    json_t* object = json_object();
    json_object_set_new(object, "left", json_integer(node.left));
    json_object_set_new(object, "right", json_integer(node.right));
    json_object_set_new(object, "parent", json_integer(node.parent));
    json_object_set_new(object, "rawTime", json_real(node.rawTime));
    json_object_set_new(object, "normTime", json_real(node.normTime));
    json_object_set_new(object, "descendantMask", json_integer(node.descendantMask));
    json_object_set_new(object, "displayX", json_real(node.displayX));
    return object;
}

json_t* mutationToJson(const Mutation& mutation) {
    json_t* object = json_object();
    json_object_set_new(object, "branchChild", json_integer(mutation.branchChild));
    json_object_set_new(object, "rawTime", json_real(mutation.rawTime));
    json_object_set_new(object, "normTime", json_real(mutation.normTime));
    json_object_set_new(object, "descendantMask",
                        json_integer(mutation.descendantMask));
    json_object_set_new(object, "sign", json_integer(mutation.sign));
    return object;
}

json_t* treeToJson(const Tree& tree) {
    json_t* object = json_object();
    json_object_set_new(object, "leafCount", json_integer(tree.leafCount));
    json_object_set_new(object, "nodeCount", json_integer(tree.nodeCount));
    json_object_set_new(object, "root", json_integer(tree.root));
    json_object_set_new(object, "mutationCount", json_integer(tree.mutationCount));
    json_object_set_new(object, "generatedMutate", json_real(tree.generatedMutate));
    json_object_set_new(object, "mutationLambda", json_real(tree.mutationLambda));
    json_object_set_new(object, "totalRawBranchLength",
                        json_real(tree.totalRawBranchLength));
    json_object_set_new(object, "mutationsTruncated",
                        json_boolean(tree.mutationsTruncated));

    json_t* nodes = json_array();
    for (int index = 0; index < tree.nodeCount; ++index)
        json_array_append_new(nodes, nodeToJson(tree.nodes[index]));
    json_object_set_new(object, "nodes", nodes);

    json_t* leafOrder = json_array();
    for (int index = 0; index < tree.leafCount; ++index)
        json_array_append_new(leafOrder, json_integer(tree.leafOrder[index]));
    json_object_set_new(object, "leafOrder", leafOrder);

    json_t* mutations = json_array();
    for (int index = 0; index < tree.mutationCount; ++index)
        json_array_append_new(mutations, mutationToJson(tree.mutations[index]));
    json_object_set_new(object, "mutations", mutations);
    return object;
}

bool readNode(json_t* object, Node& destination) {
    if (!object || !json_is_object(object))
        return false;
    int left = -1;
    int right = -1;
    int parent = -1;
    int mask = 0;
    double rawTime = 0.0;
    double normTime = 0.0;
    double displayX = 0.0;
    if (!readInteger(object, "left", -1, MAX_NODES - 1, left) ||
        !readInteger(object, "right", -1, MAX_NODES - 1, right) ||
        !readInteger(object, "parent", -1, MAX_NODES - 1, parent) ||
        !readInteger(object, "descendantMask", 0, 0xffff, mask) ||
        !readFinite(object, "rawTime", rawTime) ||
        !readFiniteRange(object, "normTime", 0.0, 1.0, normTime) ||
        !readFiniteRange(object, "displayX", 0.0, 1.0, displayX))
        return false;
    destination.left = left;
    destination.right = right;
    destination.parent = parent;
    destination.rawTime = rawTime;
    destination.normTime = normTime;
    destination.descendantMask = static_cast<std::uint16_t>(mask);
    destination.displayX = displayX;
    return true;
}

bool readMutation(json_t* object, Mutation& destination) {
    if (!object || !json_is_object(object))
        return false;
    int branchChild = -1;
    int mask = 0;
    int sign = 0;
    double rawTime = 0.0;
    double normTime = 0.0;
    if (!readInteger(object, "branchChild", 0, MAX_NODES - 1, branchChild) ||
        !readInteger(object, "descendantMask", 0, 0xffff, mask) ||
        !readInteger(object, "sign", -1, 1, sign) ||
        !readFinite(object, "rawTime", rawTime) ||
        !readFiniteRange(object, "normTime", 0.0, 1.0, normTime) ||
        (sign != -1 && sign != 1))
        return false;
    destination.branchChild = branchChild;
    destination.rawTime = rawTime;
    destination.normTime = normTime;
    destination.descendantMask = static_cast<std::uint16_t>(mask);
    destination.sign = static_cast<std::int8_t>(sign);
    return true;
}

bool derivedTreeMatches(const Tree& serialized, const Tree& derived) {
    if (serialized.leafCount != derived.leafCount ||
        serialized.nodeCount != derived.nodeCount ||
        serialized.root != derived.root ||
        serialized.mutationCount != derived.mutationCount ||
        serialized.mutationsTruncated != derived.mutationsTruncated ||
        !nearlyEqual(serialized.generatedMutate, derived.generatedMutate) ||
        !nearlyEqual(serialized.mutationLambda, derived.mutationLambda) ||
        !nearlyEqual(serialized.totalRawBranchLength,
                     derived.totalRawBranchLength))
        return false;

    for (int index = 0; index < serialized.nodeCount; ++index) {
        const Node& expected = serialized.nodes[index];
        const Node& actual = derived.nodes[index];
        if (expected.left != actual.left || expected.right != actual.right ||
            expected.parent != actual.parent ||
            expected.descendantMask != actual.descendantMask ||
            !nearlyEqual(expected.rawTime, actual.rawTime) ||
            !nearlyEqual(expected.normTime, actual.normTime) ||
            !nearlyEqual(expected.displayX, actual.displayX))
            return false;
    }
    for (int index = 0; index < serialized.leafCount; ++index)
        if (serialized.leafOrder[index] != derived.leafOrder[index])
            return false;
    for (int index = 0; index < serialized.mutationCount; ++index) {
        const Mutation& expected = serialized.mutations[index];
        const Mutation& actual = derived.mutations[index];
        if (expected.branchChild != actual.branchChild ||
            expected.descendantMask != actual.descendantMask ||
            expected.sign != actual.sign ||
            !nearlyEqual(expected.rawTime, actual.rawTime) ||
            !nearlyEqual(expected.normTime, actual.normTime))
            return false;
    }
    return true;
}

bool readTree(json_t* object, Tree& destination) {
    if (!object || !json_is_object(object))
        return false;
    Tree serialized;
    if (!readInteger(object, "leafCount", 2, MAX_LEAVES,
                     serialized.leafCount) ||
        !readInteger(object, "nodeCount", 3, MAX_NODES,
                     serialized.nodeCount) ||
        !readInteger(object, "root", 0, MAX_NODES - 1, serialized.root) ||
        !readInteger(object, "mutationCount", 0, MAX_MUTATIONS,
                     serialized.mutationCount) ||
        !readFiniteRange(object, "generatedMutate", 0.0, 1.0,
                         serialized.generatedMutate) ||
        !readFinite(object, "mutationLambda", serialized.mutationLambda) ||
        !readFinite(object, "totalRawBranchLength",
                    serialized.totalRawBranchLength) ||
        !readBoolean(object, "mutationsTruncated",
                     serialized.mutationsTruncated))
        return false;

    if (serialized.nodeCount != 2 * serialized.leafCount - 1 ||
        serialized.root != serialized.nodeCount - 1 ||
        serialized.mutationLambda < 0.0 ||
        !(serialized.totalRawBranchLength > 0.0))
        return false;

    json_t* nodes = json_object_get(object, "nodes");
    json_t* leafOrder = json_object_get(object, "leafOrder");
    json_t* mutations = json_object_get(object, "mutations");
    if (!nodes || !json_is_array(nodes) ||
        json_array_size(nodes) != static_cast<std::size_t>(serialized.nodeCount) ||
        !leafOrder || !json_is_array(leafOrder) ||
        json_array_size(leafOrder) != static_cast<std::size_t>(serialized.leafCount) ||
        !mutations || !json_is_array(mutations) ||
        json_array_size(mutations) !=
            static_cast<std::size_t>(serialized.mutationCount))
        return false;

    for (int index = 0; index < serialized.nodeCount; ++index)
        if (!readNode(json_array_get(nodes, index), serialized.nodes[index]))
            return false;
    for (int index = 0; index < serialized.leafCount; ++index) {
        json_t* item = json_array_get(leafOrder, index);
        if (!item || !json_is_integer(item))
            return false;
        const json_int_t leaf = json_integer_value(item);
        if (leaf < 0 || leaf >= serialized.leafCount)
            return false;
        serialized.leafOrder[index] = static_cast<int>(leaf);
    }
    for (int index = 0; index < serialized.mutationCount; ++index)
        if (!readMutation(json_array_get(mutations, index),
                          serialized.mutations[index]))
            return false;

    Tree derived = serialized;
    if (!derived.validateAndDerive() || !derivedTreeMatches(serialized, derived))
        return false;
    destination = derived;
    return true;
}

json_t* snapshotToJson(const Snapshot& saved, int generatedMutateKey) {
    json_t* root = json_object();
    json_object_set_new(root, "lineagesStateVersion",
                        json_integer(coalescent::lineages::STATE_VERSION));
    json_object_set_new(root, "generatedSamples",
                        json_integer(saved.generated.tree.leafCount));
    json_object_set_new(root, "generatedMutateKey",
                        json_integer(generatedMutateKey));
    json_object_set_new(root, "tree", treeToJson(saved.generated.tree));

    json_t* rng = json_object();
    setHex(rng, "state", saved.generated.nextRng.state);
    setHex(rng, "increment", saved.generated.nextRng.increment);
    json_object_set_new(root, "nextRng", rng);

    json_t* transport = json_object();
    json_object_set_new(transport, "cursor", json_real(saved.cursor));
    json_object_set_new(transport, "direction",
                        json_integer(saved.direction == Direction::Descent ? 1 : 0));
    json_object_set_new(transport, "loop", json_boolean(saved.loop));
    json_object_set_new(transport, "running", json_boolean(saved.running));
    json_object_set_new(root, "transport", transport);

    json_t* pulses = json_object();
    json_object_set_new(pulses, "node", json_real(saved.nodePulseRemaining));
    json_object_set_new(pulses, "mutation",
                        json_real(saved.mutationPulseRemaining));
    json_object_set_new(pulses, "mrca", json_real(saved.mrcaPulseRemaining));
    json_object_set_new(root, "pulses", pulses);
    return root;
}

bool snapshotFromJson(json_t* root, Snapshot& destination,
                      int& generatedMutateKey) {
    if (!root || !json_is_object(root))
        return false;
    int version = 0;
    int generatedSamples = 0;
    int mutateKey = 0;
    if (!readInteger(root, "lineagesStateVersion",
                     coalescent::lineages::STATE_VERSION,
                     coalescent::lineages::STATE_VERSION, version) ||
        !readInteger(root, "generatedSamples", 2, MAX_LEAVES,
                     generatedSamples) ||
        !readInteger(root, "generatedMutateKey", 0, MUTATE_QUANTIZATION,
                     mutateKey))
        return false;

    Snapshot candidate;
    candidate.version = version;
    json_t* tree = json_object_get(root, "tree");
    if (!readTree(tree, candidate.generated.tree) ||
        candidate.generated.tree.leafCount != generatedSamples ||
        !nearlyEqual(candidate.generated.tree.generatedMutate,
                     mutateFromKey(mutateKey)))
        return false;

    json_t* rng = json_object_get(root, "nextRng");
    if (!rng || !json_is_object(rng) ||
        !readHex(rng, "state", candidate.generated.nextRng.state) ||
        !readHex(rng, "increment", candidate.generated.nextRng.increment))
        return false;

    json_t* transport = json_object_get(root, "transport");
    int directionValue = 0;
    if (!transport || !json_is_object(transport) ||
        !readFiniteRange(transport, "cursor", 0.0, 1.0,
                         candidate.cursor) ||
        !readInteger(transport, "direction", 0, 1, directionValue) ||
        !readBoolean(transport, "loop", candidate.loop) ||
        !readBoolean(transport, "running", candidate.running))
        return false;
    candidate.direction = directionValue == 0
        ? Direction::Ancestry : Direction::Descent;

    json_t* pulses = json_object_get(root, "pulses");
    double nodePulse = 0.0;
    double mutationPulse = 0.0;
    double mrcaPulse = 0.0;
    const double maximumPulse = Snapshot::maximumPulseRemaining();
    if (!pulses || !json_is_object(pulses) ||
        !readFiniteRange(pulses, "node", 0.0, maximumPulse, nodePulse) ||
        !readFiniteRange(pulses, "mutation", 0.0, maximumPulse,
                         mutationPulse) ||
        !readFiniteRange(pulses, "mrca", 0.0, maximumPulse, mrcaPulse))
        return false;
    candidate.nodePulseRemaining = static_cast<float>(nodePulse);
    candidate.mutationPulseRemaining = static_cast<float>(mutationPulse);
    candidate.mrcaPulseRemaining = static_cast<float>(mrcaPulse);

    if (!candidate.validateAndDerive())
        return false;
    destination = candidate;
    generatedMutateKey = mutateKey;
    return true;
}

} // namespace

struct LineagesMutationQuantity : ParamQuantity {
    float getDisplayValue() override {
        Param* parameter = getParam();
        const float value = parameter ? parameter->getValue() : defaultValue;
        const float unit = clampFloat(finiteOr(value, 0.4f), 0.f, 1.f);
        return 64.f * unit * unit;
    }

    void setDisplayValue(float displayValue) override {
        if (!std::isfinite(displayValue))
            return;
        const float expected = clampFloat(displayValue, 0.f, 64.f);
        setImmediateValue(std::sqrt(expected / 64.f));
    }
};

struct Lineages : Module {
    enum ParamIds {
        RATE_PARAM,
        SAMPLES_PARAM,
        MUTATE_PARAM,
        STEP_PARAM,
        DIRECTION_PARAM,
        LOOP_PARAM,
        NEW_PARAM,
        RESET_PARAM,
        PARAMS_LEN
    };
    enum InputIds {
        RATE_INPUT,
        TIME_INPUT,
        NEW_INPUT,
        RESET_INPUT,
        INPUTS_LEN
    };
    enum OutputIds {
        TRAITS_OUTPUT,
        NODE_OUTPUT,
        MUTATION_OUTPUT,
        MRCA_OUTPUT,
        LINEAGES_OUTPUT,
        DIVERSITY_OUTPUT,
        OUTPUTS_LEN
    };
    enum LightIds {
        DIRTY_LIGHT,
        LIGHTS_LEN
    };

    KingmanGenerator generator;
    Tree tree;
    Playback playback;
    int generatedMutateKey = mutateKeyFromValue(0.4f);

    dsp::SchmittTrigger newInputTrigger;
    dsp::SchmittTrigger newButtonTrigger;
    dsp::SchmittTrigger resetInputTrigger;
    dsp::SchmittTrigger resetButtonTrigger;
    dsp::PulseGenerator nodePulse;
    dsp::PulseGenerator mutationPulse;
    dsp::PulseGenerator mrcaPulse;

    bool timeConnectionKnown = false;
    bool timeWasConnected = false;
    dsp::ClockDivider displayDivider;
    dsp::ClockDivider saveDivider;

    struct DisplayNode {
        float x = 0.f;
        float time = 0.f;
        int parent = -1;
        std::uint16_t descendantMask = 0u;
    };

    struct DisplayMutation {
        float time = 0.f;
        int branchChild = -1;
        std::int8_t sign = 1;
    };

    struct DisplayFrame {
        DisplayNode nodes[MAX_NODES] = {};
        DisplayMutation mutations[MAX_MUTATIONS] = {};
        std::uint16_t activeMasks[MAX_LEAVES] = {};
        int sampleCount = 0;
        int nodeCount = 0;
        int mutationCount = 0;
        int root = -1;
        int activeCount = 0;
        float cursor = 1.f;
        bool ancestry = false;
        bool dirty = false;
        bool valid = false;
        std::uint64_t generation = 0u;
    };
    coalescent::DisplaySnapshot<DisplayFrame> displaySnapshot;
    std::uint64_t displayGeneration = 0u;

    struct SaveFrame {
        Snapshot state;
        int generatedMutateKey = mutateKeyFromValue(0.4f);
        bool valid = false;
    };
    coalescent::DisplaySnapshot<SaveFrame> saveSnapshot;

    Lineages() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);
        configParam(RATE_PARAM, -6.f, 5.f, 0.f, "Traversal rate", " trees/s",
                    2.f, RATE_BASE_HZ);
        configParam(SAMPLES_PARAM, 2.f, 16.f, 8.f, "Samples");
        getParamQuantity(SAMPLES_PARAM)->snapEnabled = true;
        configParam<LineagesMutationQuantity>(
            MUTATE_PARAM, 0.f, 1.f, 0.4f,
            "Mutation density", " expected mutations");
        configParam(STEP_PARAM, 0.f, 12.f, 2.f, "Mutation step", " semitones");
        configSwitch(DIRECTION_PARAM, 0.f, 1.f, 1.f, "Direction",
                     {"Ancestry", "Descent"});
        configSwitch(LOOP_PARAM, 0.f, 1.f, 1.f, "Loop", {"Off", "On"});
        configButton(NEW_PARAM, "Generate new tree");
        configButton(RESET_PARAM, "Replay current tree");

        configInput(RATE_INPUT, "Traversal rate (1 V/oct)");
        configInput(TIME_INPUT, "Absolute ancestral age (0 to 10 V)");
        configInput(NEW_INPUT, "Generate new tree trigger");
        configInput(RESET_INPUT, "Replay current tree trigger");

        configOutput(TRAITS_OUTPUT, "Sample traits (polyphonic)");
        configOutput(NODE_OUTPUT, "Node crossing trigger");
        configOutput(MUTATION_OUTPUT, "Mutation crossing trigger");
        configOutput(MRCA_OUTPUT, "Most recent common ancestor trigger");
        configOutput(LINEAGES_OUTPUT, "Active lineage count");
        configOutput(DIVERSITY_OUTPUT, "Mean pairwise mutation diversity");
        configLight(DIRTY_LIGHT, "New tree pending");

        displayDivider.setDivision(1024);
        saveDivider.setDivision(96);
        seedFromRack();
        if (!generateMatchingTree(false))
            installFactoryFallback();
        resetEdges();
        publishDisplayFrame();
        publishSaveFrame();
    }

    int currentSampleCount() {
        const float value = clampFloat(
            finiteOr(params[SAMPLES_PARAM].getValue(), 8.f), 2.f, 16.f);
        int count = static_cast<int>(std::lround(value));
        if (count < 2)
            count = 2;
        if (count > MAX_LEAVES)
            count = MAX_LEAVES;
        return count;
    }

    int currentMutateKey() {
        return mutateKeyFromValue(params[MUTATE_PARAM].getValue());
    }

    Direction currentDirection() {
        return finiteOr(params[DIRECTION_PARAM].getValue(), 1.f) >= 0.5f
            ? Direction::Descent : Direction::Ancestry;
    }

    bool currentLoop() {
        return finiteOr(params[LOOP_PARAM].getValue(), 1.f) >= 0.5f;
    }

    bool isDirty() {
        return currentSampleCount() != tree.leafCount ||
               currentMutateKey() != generatedMutateKey;
    }

    void markTreeInstalled() {
        ++displayGeneration;
        if (displayGeneration == 0u)
            ++displayGeneration;
    }

    static void fillDisplayFrame(DisplayFrame& frame, const Tree& sourceTree,
                                 const Playback& sourcePlayback,
                                 bool sourceDirty, std::uint64_t generation) {
        frame = DisplayFrame();
        frame.cursor = clampFloat(finiteOr(
            static_cast<float>(sourcePlayback.cursor()), 1.f), 0.f, 1.f);
        frame.ancestry = sourcePlayback.direction() == Direction::Ancestry;
        frame.dirty = sourceDirty;
        frame.generation = generation;

        if (sourceTree.leafCount < 2 || sourceTree.leafCount > MAX_LEAVES ||
            sourceTree.nodeCount != 2 * sourceTree.leafCount - 1 ||
            sourceTree.root < 0 || sourceTree.root >= sourceTree.nodeCount ||
            sourceTree.mutationCount < 0 ||
            sourceTree.mutationCount > MAX_MUTATIONS)
            return;

        frame.sampleCount = sourceTree.leafCount;
        frame.nodeCount = sourceTree.nodeCount;
        frame.mutationCount = sourceTree.mutationCount;
        frame.root = sourceTree.root;
        for (int index = 0; index < frame.nodeCount; ++index) {
            const Node& node = sourceTree.nodes[index];
            DisplayNode& displayNode = frame.nodes[index];
            displayNode.x = clampFloat(
                finiteOr(static_cast<float>(node.displayX), 0.5f), 0.f, 1.f);
            displayNode.time = clampFloat(
                finiteOr(static_cast<float>(node.normTime)), 0.f, 1.f);
            displayNode.parent = node.parent;
            displayNode.descendantMask = node.descendantMask;
        }
        for (int index = 0; index < frame.mutationCount; ++index) {
            const Mutation& mutation = sourceTree.mutations[index];
            DisplayMutation& displayMutation = frame.mutations[index];
            displayMutation.time = clampFloat(
                finiteOr(static_cast<float>(mutation.normTime)), 0.f, 1.f);
            displayMutation.branchChild = mutation.branchChild;
            displayMutation.sign = mutation.sign < 0 ? -1 : 1;
        }

        const coalescent::lineages::PlaybackValues& values =
            sourcePlayback.values();
        frame.activeCount = values.activeLineageCount;
        if (frame.activeCount < 0)
            frame.activeCount = 0;
        if (frame.activeCount > frame.sampleCount)
            frame.activeCount = frame.sampleCount;
        for (int index = 0; index < frame.activeCount; ++index)
            frame.activeMasks[index] = values.activeLineageMasks[index];
        frame.valid = true;
    }

    static DisplayFrame previewDisplayFrame() {
        KingmanGenerator previewGenerator;
        Tree previewTree;
        Playback previewPlayback;
        DisplayFrame frame;
        if (!previewGenerator.generate(8, 0.4, previewTree))
            return frame;
        previewPlayback.setDirection(Direction::Descent);
        previewPlayback.setLoop(true);
        previewPlayback.installTreeAtSource(previewTree);
        previewPlayback.scrub(previewTree, 0.56);
        fillDisplayFrame(frame, previewTree, previewPlayback, false, 1u);
        return frame;
    }

    void publishDisplayFrame() {
        DisplayFrame& frame = displaySnapshot.writable();
        fillDisplayFrame(frame, tree, playback, isDirty(), displayGeneration);
        displaySnapshot.publish();
    }

    void seedFromRack() {
        generator.seed(random::u64(), random::u64() >> 1u);
    }

    void resetEdges() {
        newInputTrigger.reset();
        newButtonTrigger.reset();
        resetInputTrigger.reset();
        resetButtonTrigger.reset();
    }

    void resetPulses() {
        nodePulse.reset();
        mutationPulse.reset();
        mrcaPulse.reset();
    }

    bool generateMatchingTree(bool emitMrca) {
        const int samples = currentSampleCount();
        const int mutateKey = currentMutateKey();
        Tree candidate;
        if (!generator.generate(samples, mutateFromKey(mutateKey), candidate))
            return false;

        tree = candidate;
        generatedMutateKey = mutateKey;
        playback.setDirection(currentDirection());
        playback.setLoop(currentLoop());
        const EventFacts facts = playback.installTreeAtSource(tree);
        markTreeInstalled();
        if (emitMrca)
            triggerFacts(facts);
        return true;
    }

    void installFactoryFallback() {
        Snapshot factory;
        if (!coalescent::lineages::makeFactorySnapshot(factory) ||
            !factory.install(generator, tree, playback)) {
            // The factory path is exhaustively tested and has no runtime input.
            // Keep a finite, silent module even if a future core change breaks it.
            tree = Tree();
            playback = Playback();
        }
        markTreeInstalled();
        generatedMutateKey = mutateKeyFromValue(
            static_cast<float>(coalescent::lineages::FACTORY_MUTATE));
        params[SAMPLES_PARAM].setValue(
            static_cast<float>(coalescent::lineages::FACTORY_LEAF_COUNT));
        params[MUTATE_PARAM].setValue(
            static_cast<float>(coalescent::lineages::FACTORY_MUTATE));
        params[DIRECTION_PARAM].setValue(1.f);
        params[LOOP_PARAM].setValue(1.f);
        resetPulses();
        resetEdges();
        timeConnectionKnown = false;
        timeWasConnected = false;
    }

    void triggerFacts(const EventFacts& facts) {
        if (facts.nodeCrossed)
            nodePulse.trigger(PULSE_SECONDS);
        if (facts.mutationCrossed)
            mutationPulse.trigger(PULSE_SECONDS);
        if (facts.requestsMrcaPulse())
            mrcaPulse.trigger(PULSE_SECONDS);
    }

    double timeTarget() {
        const float voltage = inputs[TIME_INPUT].getVoltage();
        if (!std::isfinite(voltage))
            return playback.cursor();
        return clampDouble(static_cast<double>(voltage) * 0.1, 0.0, 1.0);
    }

    float traversalRate() {
        float octave = finiteOr(params[RATE_PARAM].getValue(), 0.f);
        octave += finiteOr(inputs[RATE_INPUT].getVoltage(), 0.f);
        octave = clampFloat(octave, RATE_TOTAL_MIN, RATE_TOTAL_MAX);
        const float rate = RATE_BASE_HZ * std::exp2(octave);
        return std::isfinite(rate) && rate >= 0.f ? rate : 0.f;
    }

    void publishSaveFrame() {
        SaveFrame& frame = saveSnapshot.writable();
        frame.state.version = coalescent::lineages::STATE_VERSION;
        frame.state.generated = generator.capture(tree);
        frame.state.cursor = playback.cursor();
        frame.state.direction = playback.direction();
        frame.state.loop = playback.loop();
        frame.state.running = playback.running();
        frame.state.nodePulseRemaining = clampFloat(
            finiteOr(nodePulse.remaining), 0.f, Snapshot::maximumPulseRemaining());
        frame.state.mutationPulseRemaining = clampFloat(
            finiteOr(mutationPulse.remaining), 0.f,
            Snapshot::maximumPulseRemaining());
        frame.state.mrcaPulseRemaining = clampFloat(
            finiteOr(mrcaPulse.remaining), 0.f, Snapshot::maximumPulseRemaining());
        frame.generatedMutateKey = generatedMutateKey;
        frame.valid = tree.leafCount >= 2 && tree.leafCount <= MAX_LEAVES &&
                      tree.nodeCount == 2 * tree.leafCount - 1 &&
                      tree.root == tree.nodeCount - 1;
        saveSnapshot.publish();
    }

    void onSampleRateChange(const SampleRateChangeEvent& event) override {
        const float sampleRate = finiteOr(event.sampleRate, 48000.f);
        int displayDivision = static_cast<int>(
            std::lround(sampleRate / DISPLAY_HZ));
        if (displayDivision < 1)
            displayDivision = 1;
        displayDivider.setDivision(displayDivision);
        displayDivider.reset();
        int saveDivision = static_cast<int>(std::lround(sampleRate / SAVE_HZ));
        if (saveDivision < 1)
            saveDivision = 1;
        saveDivider.setDivision(saveDivision);
        saveDivider.reset();
    }

    void onReset(const ResetEvent& event) override {
        Module::onReset(event);
        resetPulses();
        if (!generateMatchingTree(false))
            installFactoryFallback();
        resetEdges();
        timeConnectionKnown = false;
        timeWasConnected = false;
        publishDisplayFrame();
        publishSaveFrame();
    }

    void onRandomize(const RandomizeEvent& event) override {
        Module::onRandomize(event);
        seedFromRack();
        resetPulses();
        if (!generateMatchingTree(false))
            installFactoryFallback();
        resetEdges();
        timeConnectionKnown = false;
        timeWasConnected = false;
        publishDisplayFrame();
        publishSaveFrame();
    }

    json_t* dataToJson() override {
        const SaveFrame& frame = saveSnapshot.consume();
        if (!frame.valid)
            return json_object();
        return snapshotToJson(frame.state, frame.generatedMutateKey);
    }

    void dataFromJson(json_t* root) override {
        Snapshot restored;
        int restoredMutateKey = 0;
        if (!snapshotFromJson(root, restored, restoredMutateKey) ||
            !restored.install(generator, tree, playback)) {
            installFactoryFallback();
            publishDisplayFrame();
            publishSaveFrame();
            return;
        }

        generatedMutateKey = restoredMutateKey;
        markTreeInstalled();
        params[DIRECTION_PARAM].setValue(
            restored.direction == Direction::Descent ? 1.f : 0.f);
        params[LOOP_PARAM].setValue(restored.loop ? 1.f : 0.f);
        nodePulse.remaining = restored.nodePulseRemaining;
        mutationPulse.remaining = restored.mutationPulseRemaining;
        mrcaPulse.remaining = restored.mrcaPulseRemaining;
        resetEdges();
        timeConnectionKnown = false;
        timeWasConnected = false;
        publishDisplayFrame();
        publishSaveFrame();
    }

    void writeContinuousOutputs() {
        const coalescent::lineages::PlaybackValues& values = playback.values();
        int samples = tree.leafCount;
        if (samples < 2 || samples > MAX_LEAVES)
            samples = 2;
        outputs[TRAITS_OUTPUT].setChannels(samples);

        const float step = clampFloat(
            finiteOr(params[STEP_PARAM].getValue(), 2.f), 0.f, 12.f);
        const float voltsPerMutation = step / 12.f;
        for (int channel = 0; channel < samples; ++channel) {
            float voltage = static_cast<float>(values.traitSteps[channel]) *
                            voltsPerMutation;
            voltage = clampFloat(finiteOr(voltage), -TRAIT_LIMIT, TRAIT_LIMIT);
            outputs[TRAITS_OUTPUT].setVoltage(voltage, channel);
        }

        int lineages = values.lineageCount;
        if (lineages < 1)
            lineages = 1;
        if (lineages > samples)
            lineages = samples;
        const float lineageVoltage = 10.f * static_cast<float>(lineages - 1) /
                                     static_cast<float>(samples - 1);
        outputs[LINEAGES_OUTPUT].setVoltage(
            clampFloat(finiteOr(lineageVoltage), 0.f, 10.f));

        float diversity = static_cast<float>(values.meanPairwiseMutationDistance) *
                          DIVERSITY_VOLTS_PER_MUTATION;
        diversity = clampFloat(finiteOr(diversity), 0.f, 10.f);
        outputs[DIVERSITY_OUTPUT].setVoltage(diversity);
    }

    void process(const ProcessArgs& args) override {
        const float sampleTime = std::max(0.f, finiteOr(args.sampleTime, 0.f));
        const Direction direction = currentDirection();
        const bool loop = currentLoop();
        const bool directionChanged = playback.direction() != direction;
        const bool loopChanged = playback.loop() != loop;
        playback.setDirection(direction);
        playback.setLoop(loop);

        const bool newFromInput = newInputTrigger.process(
            finiteOr(inputs[NEW_INPUT].getVoltage()), 0.1f, 1.f);
        const bool newFromButton = newButtonTrigger.process(
            finiteOr(params[NEW_PARAM].getValue()), 0.1f, 0.9f);
        const bool resetFromInput = resetInputTrigger.process(
            finiteOr(inputs[RESET_INPUT].getVoltage()), 0.1f, 1.f);
        const bool resetFromButton = resetButtonTrigger.process(
            finiteOr(params[RESET_PARAM].getValue()), 0.1f, 0.9f);
        const bool newRequested = newFromInput || newFromButton;
        const bool resetRequested = resetFromInput || resetFromButton;

        const bool timeConnected = inputs[TIME_INPUT].isConnected();
        const bool connectionChanged = timeConnectionKnown &&
                                       timeConnected != timeWasConnected;
        const bool initialTimeConnection = !timeConnectionKnown && timeConnected;
        bool immediateSave = directionChanged || loopChanged;

        if (newRequested) {
            Tree candidate;
            const int samples = currentSampleCount();
            const int mutateKey = currentMutateKey();
            if (generator.generate(samples, mutateFromKey(mutateKey), candidate)) {
                tree = candidate;
                generatedMutateKey = mutateKey;
                resetPulses();
                EventFacts facts;
                if (timeConnected)
                    playback.synchronize(tree, timeTarget());
                else
                    facts = playback.installTreeAtSource(tree);
                markTreeInstalled();
                triggerFacts(facts);
                immediateSave = true;
            }
        }
        else if (resetRequested) {
            resetPulses();
            EventFacts facts;
            if (timeConnected)
                playback.synchronize(tree, timeTarget());
            else
                facts = playback.resetToSource(tree);
            triggerFacts(facts);
            immediateSave = true;
        }
        else if (connectionChanged || initialTimeConnection) {
            if (timeConnected)
                playback.synchronize(tree, timeTarget());
            else
                playback.synchronize(tree, playback.cursor());
            immediateSave = true;
        }
        else if (timeConnected) {
            triggerFacts(playback.scrub(tree, timeTarget()));
        }
        else if (!directionChanged) {
            const double distance = static_cast<double>(traversalRate()) *
                                    static_cast<double>(sampleTime);
            triggerFacts(playback.advance(tree, distance));
        }

        timeConnectionKnown = true;
        timeWasConnected = timeConnected;

        writeContinuousOutputs();
        outputs[NODE_OUTPUT].setVoltage(
            nodePulse.process(sampleTime) ? PULSE_VOLTS : 0.f);
        outputs[MUTATION_OUTPUT].setVoltage(
            mutationPulse.process(sampleTime) ? PULSE_VOLTS : 0.f);
        outputs[MRCA_OUTPUT].setVoltage(
            mrcaPulse.process(sampleTime) ? PULSE_VOLTS : 0.f);
        lights[DIRTY_LIGHT].setBrightness(isDirty() ? 1.f : 0.f);

        const bool displayTick = displayDivider.process();
        if (immediateSave || displayTick)
            publishDisplayFrame();
        if (immediateSave || saveDivider.process())
            publishSaveFrame();
    }
};

namespace lineages_layout {
constexpr float COL4[4] = {10.f, 30.43f, 50.85f, 71.28f};
constexpr float COL3[3] = {14.f, 40.64f, 67.28f};
constexpr float KNOB_Y = 56.5f;
constexpr float ACTION_Y = 72.5f;
constexpr float INPUT_Y = 88.5f;
constexpr float OUTPUT_1_Y = 103.5f;
constexpr float OUTPUT_2_Y = 118.2f;
} // namespace lineages_layout

struct LineagesTreeView : widget::TransparentWidget {
    Lineages* module = nullptr;

    void draw(const DrawArgs& args) override {
        nvgBeginPath(args.vg);
        nvgRoundedRect(args.vg, 0.f, 0.f, box.size.x, box.size.y,
                       mm2px(1.7f));
        nvgFillColor(args.vg, nvgRGB(0x06, 0x12, 0x10));
        nvgFill(args.vg);
        nvgStrokeColor(args.vg, nvgRGB(0x3b, 0x57, 0x50));
        nvgStrokeWidth(args.vg, 1.2f);
        nvgStroke(args.vg);
        Widget::draw(args);
    }

    void drawLayer(const DrawArgs& args, int layer) override {
        if (layer != 1)
            return;

        static const Lineages::DisplayFrame preview =
            Lineages::previewDisplayFrame();
        const Lineages::DisplayFrame& frame = module
            ? module->displaySnapshot.consume() : preview;
        const float width = box.size.x;
        const float height = box.size.y;
        const float left = mm2px(2.8f);
        const float right = width - mm2px(2.8f);
        const float top = mm2px(7.2f);
        const float bottom = height - mm2px(4.5f);
        const float plotWidth = right - left;
        const float plotHeight = bottom - top;
        auto nodeX = [&](int index) {
            return left + frame.nodes[index].x * plotWidth;
        };
        auto timeY = [&](float time) {
            return bottom - clamp(time, 0.f, 1.f) * plotHeight;
        };

        nvgScissor(args.vg, 0.f, 0.f, width, height);
        if (frame.valid) {
            for (int guide = 1; guide < 4; ++guide) {
                const float y = timeY(static_cast<float>(guide) * 0.25f);
                nvgBeginPath(args.vg);
                nvgMoveTo(args.vg, left, y);
                nvgLineTo(args.vg, right, y);
                nvgStrokeColor(args.vg, nvgRGBA(0x79, 0xb2, 0xa4, 0x12));
                nvgStrokeWidth(args.vg, 1.f);
                nvgStroke(args.vg);
            }

            // The complete genealogy remains visible. Each child branch rises
            // vertically to its parent's age, then joins at the parent's X.
            nvgBeginPath(args.vg);
            for (int child = 0; child < frame.nodeCount; ++child) {
                const int parent = frame.nodes[child].parent;
                if (parent < 0 || parent >= frame.nodeCount)
                    continue;
                const float childX = nodeX(child);
                const float parentX = nodeX(parent);
                const float childY = timeY(frame.nodes[child].time);
                const float parentY = timeY(frame.nodes[parent].time);
                nvgMoveTo(args.vg, childX, childY);
                nvgLineTo(args.vg, childX, parentY);
                nvgLineTo(args.vg, parentX, parentY);
            }
            nvgStrokeColor(args.vg, nvgRGBA(0x75, 0xb4, 0xa5, 0x56));
            nvgStrokeWidth(args.vg, 1.05f);
            nvgLineCap(args.vg, NVG_ROUND);
            nvgLineJoin(args.vg, NVG_ROUND);
            nvgStroke(args.vg);

            // Highlight the portion already traversed in the selected direction.
            nvgBeginPath(args.vg);
            for (int child = 0; child < frame.nodeCount; ++child) {
                const int parent = frame.nodes[child].parent;
                if (parent < 0 || parent >= frame.nodeCount)
                    continue;
                const float childTime = frame.nodes[child].time;
                const float parentTime = frame.nodes[parent].time;
                const float childX = nodeX(child);
                const float parentX = nodeX(parent);
                const float parentY = timeY(parentTime);
                if (frame.ancestry) {
                    if (frame.cursor > childTime) {
                        const float endTime = std::min(parentTime, frame.cursor);
                        nvgMoveTo(args.vg, childX, timeY(childTime));
                        nvgLineTo(args.vg, childX, timeY(endTime));
                    }
                    if (frame.cursor >= parentTime) {
                        nvgMoveTo(args.vg, childX, parentY);
                        nvgLineTo(args.vg, parentX, parentY);
                    }
                }
                else {
                    if (frame.cursor < parentTime) {
                        const float startTime = std::max(childTime, frame.cursor);
                        nvgMoveTo(args.vg, childX, timeY(startTime));
                        nvgLineTo(args.vg, childX, parentY);
                    }
                    if (frame.cursor <= parentTime) {
                        nvgMoveTo(args.vg, childX, parentY);
                        nvgLineTo(args.vg, parentX, parentY);
                    }
                }
            }
            nvgStrokeColor(args.vg, nvgRGBA(0xa0, 0xe3, 0xcf, 0xd2));
            nvgStrokeWidth(args.vg, 1.35f);
            nvgLineCap(args.vg, NVG_ROUND);
            nvgLineJoin(args.vg, NVG_ROUND);
            nvgStroke(args.vg);

            // Signed mutations use both tint and slash orientation. Mutations
            // currently inherited by the leaves are brighter than future ones.
            for (int brightPass = 0; brightPass < 2; ++brightPass) {
                for (int signPass = 0; signPass < 2; ++signPass) {
                    nvgBeginPath(args.vg);
                    for (int index = 0; index < frame.mutationCount; ++index) {
                        const Lineages::DisplayMutation& mutation =
                            frame.mutations[index];
                        const int child = mutation.branchChild;
                        if (child < 0 || child >= frame.nodeCount)
                            continue;
                        const bool inherited = mutation.time > frame.cursor;
                        const bool positive = mutation.sign > 0;
                        if (inherited != (brightPass != 0) ||
                            positive != (signPass != 0))
                            continue;
                        const float x = nodeX(child);
                        const float y = timeY(mutation.time);
                        const float half = mm2px(0.72f);
                        if (positive) {
                            nvgMoveTo(args.vg, x - half, y + half);
                            nvgLineTo(args.vg, x + half, y - half);
                        }
                        else {
                            nvgMoveTo(args.vg, x - half, y - half);
                            nvgLineTo(args.vg, x + half, y + half);
                        }
                    }
                    const unsigned alpha = brightPass ? 0xe0 : 0x58;
                    nvgStrokeColor(args.vg, signPass
                        ? nvgRGBA(0xf0, 0xc5, 0x6d, alpha)
                        : nvgRGBA(0xe1, 0x86, 0x86, alpha));
                    nvgStrokeWidth(args.vg, brightPass ? 1.35f : 1.f);
                    nvgLineCap(args.vg, NVG_ROUND);
                    nvgStroke(args.vg);
                }
            }

            nvgBeginPath(args.vg);
            for (int index = 0; index < frame.sampleCount; ++index)
                nvgCircle(args.vg, nodeX(index), bottom, mm2px(0.35f));
            nvgFillColor(args.vg, nvgRGBA(0xc4, 0xe8, 0xde, 0xa8));
            nvgFill(args.vg);

            nvgBeginPath(args.vg);
            for (int index = frame.sampleCount; index < frame.nodeCount; ++index)
                nvgCircle(args.vg, nodeX(index), timeY(frame.nodes[index].time),
                          index == frame.root ? mm2px(0.55f) : mm2px(0.32f));
            nvgFillColor(args.vg, nvgRGBA(0xd7, 0xec, 0xe6, 0xa8));
            nvgFill(args.vg);

            const float cursorY = timeY(frame.cursor);
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, left, cursorY);
            nvgLineTo(args.vg, right, cursorY);
            nvgStrokeColor(args.vg, nvgRGBA(0xeb, 0xf4, 0xef, 0xb8));
            nvgStrokeWidth(args.vg, 1.15f);
            nvgStroke(args.vg);

            nvgBeginPath(args.vg);
            for (int active = 0; active < frame.activeCount; ++active) {
                const std::uint16_t mask = frame.activeMasks[active];
                int branch = -1;
                for (int node = 0; node < frame.nodeCount; ++node) {
                    if (frame.nodes[node].descendantMask == mask) {
                        branch = node;
                        break;
                    }
                }
                if (branch >= 0)
                    nvgCircle(args.vg, nodeX(branch), cursorY, mm2px(0.62f));
            }
            nvgFillColor(args.vg, nvgRGB(0xee, 0xf7, 0xf2));
            nvgFill(args.vg);
            nvgStrokeColor(args.vg, nvgRGBA(0x55, 0xae, 0x9c, 0xe8));
            nvgStrokeWidth(args.vg, 0.9f);
            nvgStroke(args.vg);
        }
        nvgResetScissor(args.vg);

        std::shared_ptr<Font> font = APP->window->loadFont(
            asset::system("res/fonts/DejaVuSans.ttf"));
        if (font) {
            nvgFontFaceId(args.vg, font->handle);
            nvgTextLetterSpacing(args.vg, 0.f);
            nvgFontSize(args.vg, mm2px(3.2f));
            nvgFillColor(args.vg, nvgRGBA(0xa8, 0xdb, 0xc9, 0xe8));
            nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
            nvgText(args.vg, mm2px(2.6f), mm2px(2.1f), "LINEAGES", nullptr);

            nvgFontSize(args.vg, mm2px(1.8f));
            nvgFillColor(args.vg, nvgRGBA(0x94, 0xbd, 0xb2, 0xb0));
            nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_TOP);
            nvgText(args.vg, width - mm2px(4.7f), mm2px(2.6f),
                    frame.ancestry ? "ANCESTRY" : "DESCENT", nullptr);

            nvgFontSize(args.vg, mm2px(1.9f));
            nvgFillColor(args.vg, nvgRGBA(0x83, 0xa9, 0x9c, 0xa0));
            nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_BOTTOM);
            nvgText(args.vg, width - mm2px(2.4f), height - mm2px(1.5f),
                    "ancestral tree", nullptr);
        }

        if (frame.dirty) {
            nvgBeginPath(args.vg);
            nvgRoundedRect(args.vg, 0.8f, 0.8f, width - 1.6f, height - 1.6f,
                           mm2px(1.55f));
            nvgStrokeColor(args.vg, nvgRGBA(0xef, 0xc6, 0x69, 0xb0));
            nvgStrokeWidth(args.vg, 1.2f);
            nvgStroke(args.vg);
            nvgBeginPath(args.vg);
            nvgCircle(args.vg, width - mm2px(2.7f), mm2px(3.2f),
                      mm2px(0.58f));
            nvgFillColor(args.vg, nvgRGB(0xef, 0xc6, 0x69));
            nvgFill(args.vg);
        }
    }
};

struct LineagesLabels : widget::Widget {
    void draw(const DrawArgs& args) override {
        std::shared_ptr<Font> font = APP->window->loadFont(
            asset::system("res/fonts/Nunito-Bold.ttf"));
        if (!font)
            return;
        nvgFontFaceId(args.vg, font->handle);
        nvgTextLetterSpacing(args.vg, 0.f);
        nvgFillColor(args.vg, nvgRGB(0xe6, 0xf1, 0xef));
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
        auto label = [&](float x, float y, const char* text, float size = 1.45f) {
            nvgFontSize(args.vg, mm2px(size * 1.72f));
            nvgText(args.vg, mm2px(x), mm2px(y), text, nullptr);
        };

        using namespace lineages_layout;
        const char* knobs[4] = {"RATE", "SAMPLES", "MUTATE", "STEP"};
        for (int index = 0; index < 4; ++index)
            label(COL4[index], KNOB_Y - 7.7f, knobs[index],
                  index == 1 ? 1.25f : 1.42f);
        const char* actions[4] = {"DIRECTION", "LOOP", "NEW", "RESET"};
        for (int index = 0; index < 4; ++index)
            label(COL4[index], ACTION_Y - 6.0f, actions[index],
                  index == 0 ? 1.08f : 1.32f);
        const char* inputs[4] = {"RATE", "TIME", "NEW", "RESET"};
        for (int index = 0; index < 4; ++index)
            label(COL4[index], INPUT_Y - 5.5f, inputs[index], 1.30f);
        const char* top[3] = {"TRAITS", "NODE", "MUTATION"};
        const char* bottom[3] = {"MRCA", "LINEAGES", "DIVERSITY"};
        for (int index = 0; index < 3; ++index) {
            label(COL3[index], OUTPUT_1_Y - 5.5f, top[index],
                  index == 2 ? 1.13f : 1.28f);
            label(COL3[index], OUTPUT_2_Y - 5.5f, bottom[index],
                  index == 0 ? 1.28f : 1.10f);
        }
    }
};

struct LineagesWidget : ModuleWidget {
    LineagesWidget(Lineages* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Lineages.svg")));
        addPanelLabels<LineagesLabels>(this);

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.f, 1.f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(75.28f, 1.f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.f, 122.f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(75.28f, 122.f))));

        LineagesTreeView* display = new LineagesTreeView();
        display->module = module;
        display->box.pos = mm2px(Vec(6.f, 8.f));
        display->box.size = mm2px(Vec(69.28f, 38.f));
        addChild(display);

        using namespace lineages_layout;
        const int knobs[4] = {Lineages::RATE_PARAM, Lineages::SAMPLES_PARAM,
                              Lineages::MUTATE_PARAM, Lineages::STEP_PARAM};
        for (int index = 0; index < 4; ++index)
            addParam(createParamCentered<RoundBlackKnob>(
                mm2px(Vec(COL4[index], KNOB_Y)), module, knobs[index]));

        addParam(createParamCentered<CKSS>(
            mm2px(Vec(COL4[0], ACTION_Y)), module, Lineages::DIRECTION_PARAM));
        addParam(createParamCentered<CKSS>(
            mm2px(Vec(COL4[1], ACTION_Y)), module, Lineages::LOOP_PARAM));
        addParam(createParamCentered<TL1105>(
            mm2px(Vec(COL4[2], ACTION_Y)), module, Lineages::NEW_PARAM));
        addParam(createParamCentered<TL1105>(
            mm2px(Vec(COL4[3], ACTION_Y)), module, Lineages::RESET_PARAM));
        addChild(createLightCentered<SmallLight<YellowLight>>(
            mm2px(Vec(COL4[2] + 6.1f, ACTION_Y)), module,
            Lineages::DIRTY_LIGHT));

        const int inputs[4] = {Lineages::RATE_INPUT, Lineages::TIME_INPUT,
                               Lineages::NEW_INPUT, Lineages::RESET_INPUT};
        for (int index = 0; index < 4; ++index)
            addInput(createInputCentered<PJ301MPort>(
                mm2px(Vec(COL4[index], INPUT_Y)), module, inputs[index]));

        const int topOutputs[3] = {Lineages::TRAITS_OUTPUT,
                                   Lineages::NODE_OUTPUT,
                                   Lineages::MUTATION_OUTPUT};
        const int bottomOutputs[3] = {Lineages::MRCA_OUTPUT,
                                      Lineages::LINEAGES_OUTPUT,
                                      Lineages::DIVERSITY_OUTPUT};
        for (int index = 0; index < 3; ++index) {
            addOutput(createOutputCentered<PJ301MPort>(
                mm2px(Vec(COL3[index], OUTPUT_1_Y)), module,
                topOutputs[index]));
            addOutput(createOutputCentered<PJ301MPort>(
                mm2px(Vec(COL3[index], OUTPUT_2_Y)), module,
                bottomOutputs[index]));
        }
    }
};

Model* modelLineages = createModel<Lineages, LineagesWidget>("Lineages");
