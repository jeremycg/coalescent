#pragma once

#include "kingman.hpp"

#include <array>
#include <cmath>
#include <cstdint>

namespace coalescent {
namespace lineages {

enum class Direction {
    Ancestry,
    Descent
};

// Canonical musical state at one normalized ancestral age. Trait values are
// integer mutation steps; the Rack wrapper applies STEP and output clamping.
struct PlaybackValues {
    int sampleCount = 0;
    std::array<int, MAX_LEAVES> traitSteps{};
    int lineageCount = 0;
    std::array<std::uint16_t, MAX_LEAVES> activeLineageMasks{};
    int activeLineageCount = 0;
    double meanPairwiseMutationDistance = 0.0;
};

// These are biological/transport facts, not pulse state. Keeping root arrival,
// departure, and explicit start separate lets the Rack wrapper implement the
// documented MRCA pulse policy without inferring it from cursor values.
struct EventFacts {
    bool nodeCrossed = false;
    bool mutationCrossed = false;
    bool arrivedAtMrca = false;
    bool departedMrca = false;
    bool startedAtMrca = false;
    bool wrapped = false;

    bool requestsMrcaPulse() const {
        return arrivedAtMrca || startedAtMrca;
    }

    void merge(const EventFacts& other) {
        nodeCrossed = nodeCrossed || other.nodeCrossed;
        mutationCrossed = mutationCrossed || other.mutationCrossed;
        arrivedAtMrca = arrivedAtMrca || other.arrivedAtMrca;
        departedMrca = departedMrca || other.departedMrca;
        startedAtMrca = startedAtMrca || other.startedAtMrca;
        wrapped = wrapped || other.wrapped;
    }
};

class Playback {
public:
    // The approved default opens at the MRCA and branches toward the samples.
    Playback() = default;

    double cursor() const {
        return cursor_;
    }

    Direction direction() const {
        return direction_;
    }

    bool loop() const {
        return loop_;
    }

    bool running() const {
        return running_;
    }

    const PlaybackValues& values() const {
        return values_;
    }

    void setDirection(Direction direction) {
        if (direction_ == direction)
            return;
        direction_ = direction;
        // Changing direction is one of the documented ways to resume stopped
        // transport. It never emits an event by itself.
        running_ = canAdvanceFromCursor();
    }

    void setLoop(bool loop) {
        loop_ = loop;
    }

    // NEW and RESET both use this endpoint behavior when TIME is unpatched.
    // Starting descent at the root requests MRCA, but is not a node crossing.
    EventFacts installTreeAtSource(const Tree& tree) {
        cursor_ = sourceCursor();
        running_ = true;
        rebuildMutationSuffixes(tree);
        canonicalize(tree);

        EventFacts facts;
        facts.startedAtMrca = direction_ == Direction::Descent;
        return facts;
    }

    EventFacts resetToSource(const Tree& tree) {
        cursor_ = sourceCursor();
        running_ = true;
        // RESET is infrequent and this keeps the public helper correct even if
        // a caller installs restored tree data immediately before resetting.
        rebuildMutationSuffixes(tree);
        canonicalize(tree);

        EventFacts facts;
        facts.startedAtMrca = direction_ == Direction::Descent;
        return facts;
    }

    // Use for NEW/RESET while TIME is patched, cable transitions, and load.
    // Non-finite TIME holds the previous valid cursor and emits no history.
    void synchronize(const Tree& tree, double targetCursor) {
        if (std::isfinite(targetCursor))
            cursor_ = clampCursor(targetCursor);
        running_ = canAdvanceFromCursor();
        rebuildMutationSuffixes(tree);
        canonicalize(tree);
    }

    // Restore is also pulse-free. The saved running bit is authoritative; the
    // wrapper validates persisted state before calling this method.
    void restore(const Tree& tree, double savedCursor, Direction savedDirection,
                 bool savedLoop, bool savedRunning) {
        direction_ = savedDirection;
        loop_ = savedLoop;
        cursor_ = std::isfinite(savedCursor)
            ? clampCursor(savedCursor)
            : sourceCursor();
        running_ = savedRunning;
        rebuildMutationSuffixes(tree);
        canonicalize(tree);
    }

    // Absolute TIME movement aggregates any number of crossed events into one
    // fact of each type. The sorted event caches produce the same destination
    // as canonicalValues() without rescanning every mutation for a slow LFO.
    EventFacts scrub(const Tree& tree, double targetCursor) {
        EventFacts facts;
        if (!std::isfinite(targetCursor))
            return facts;

        const double target = clampCursor(targetCursor);
        if (target == cursor_)
            return facts;
        facts = traverseTo(tree, target);
        running_ = canAdvanceFromCursor();
        return facts;
    }

    // Advance by a non-negative fraction of a complete tree. Loop overshoot is
    // deliberately discarded: one call can traverse at most one endpoint,
    // then teleports to the exact source with historical crossings suppressed.
    EventFacts advance(const Tree& tree, double normalizedDistance) {
        EventFacts facts;
        if (!running_ || !(normalizedDistance > 0.0) ||
            !std::isfinite(normalizedDistance))
            return facts;

        if (direction_ == Direction::Ancestry) {
            const double distanceToEnd = 1.0 - cursor_;
            if (normalizedDistance < distanceToEnd) {
                const double target = cursor_ + normalizedDistance;
                facts = traverseTo(tree, target);
            }
            else {
                facts = traverseTo(tree, 1.0);
                if (loop_) {
                    cursor_ = 0.0;
                    facts.wrapped = true;
                    canonicalize(tree);
                }
                else {
                    cursor_ = 1.0;
                    running_ = false;
                }
            }
        }
        else {
            const double distanceToEnd = cursor_;
            if (normalizedDistance < distanceToEnd) {
                const double target = cursor_ - normalizedDistance;
                facts = traverseTo(tree, target);
            }
            else {
                facts = traverseTo(tree, 0.0);
                if (loop_) {
                    cursor_ = 1.0;
                    facts.wrapped = true;
                    facts.startedAtMrca = true;
                    canonicalize(tree);
                }
                else {
                    cursor_ = 0.0;
                    running_ = false;
                }
            }
        }

        return facts;
    }

    // Equality is assigned to the ancestral side. Moving toward the root uses
    // old < event <= new; moving toward the samples uses new < event <= old.
    static EventFacts crossingFacts(const Tree& tree, double oldCursor,
                                    double newCursor) {
        EventFacts facts;
        if (!std::isfinite(oldCursor) || !std::isfinite(newCursor))
            return facts;

        const double oldValue = clampCursor(oldCursor);
        const double newValue = clampCursor(newCursor);
        if (oldValue == newValue)
            return facts;

        const bool upward = newValue > oldValue;
        const int leafCount = safeLeafCount(tree);
        const int nodeCount = safeNodeCount(tree);

        for (int i = leafCount; i < nodeCount; ++i) {
            const double eventTime = tree.nodes[i].normTime;
            if (!validEventTime(eventTime) ||
                !crossed(eventTime, oldValue, newValue, upward))
                continue;

            facts.nodeCrossed = true;
            if (i == tree.root) {
                if (upward)
                    facts.arrivedAtMrca = true;
                else
                    facts.departedMrca = true;
            }
        }

        const int mutationCount = safeMutationCount(tree);
        for (int i = 0; i < mutationCount; ++i) {
            const double eventTime = tree.mutations[i].normTime;
            if (validEventTime(eventTime) &&
                crossed(eventTime, oldValue, newValue, upward)) {
                facts.mutationCrossed = true;
                break;
            }
        }

        return facts;
    }

    // Reference implementation used after discontinuities and by tests. A
    // mutation is active only strictly below its age; at equality it has just
    // been removed and a node at equality has already merged its descendants.
    static PlaybackValues canonicalValues(const Tree& tree, double cursor) {
        PlaybackValues result;
        const double age = std::isfinite(cursor) ? clampCursor(cursor) : 0.0;
        const int n = safeLeafCount(tree);
        result.sampleCount = n;
        if (n == 0)
            return result;

        const int nodeCount = safeNodeCount(tree);
        int mergedNodes = 0;
        for (int i = n; i < nodeCount; ++i) {
            const double eventTime = tree.nodes[i].normTime;
            if (validEventTime(eventTime) && eventTime <= age)
                ++mergedNodes;
        }
        result.lineageCount = n - mergedNodes;
        if (result.lineageCount < 1)
            result.lineageCount = 1;
        if (result.lineageCount > n)
            result.lineageCount = n;

        rebuildActiveLineages(tree, age, result);

        const std::uint16_t allLeaves = n == MAX_LEAVES
            ? static_cast<std::uint16_t>(0xffffu)
            : static_cast<std::uint16_t>((1u << n) - 1u);
        int pairwiseSeparationSum = 0;
        const int mutationCount = safeMutationCount(tree);
        for (int i = 0; i < mutationCount; ++i) {
            const Mutation& mutation = tree.mutations[i];
            if (!validEventTime(mutation.normTime) ||
                !(mutation.normTime > age) ||
                (mutation.sign != -1 && mutation.sign != 1))
                continue;

            const std::uint16_t mask = static_cast<std::uint16_t>(
                mutation.descendantMask & allLeaves);
            for (int leaf = 0; leaf < n; ++leaf) {
                if ((mask & static_cast<std::uint16_t>(1u << leaf)) != 0u)
                    result.traitSteps[leaf] += mutation.sign;
            }

            const int descendants = popcount(mask);
            pairwiseSeparationSum += descendants * (n - descendants);
        }

        result.meanPairwiseMutationDistance = pairwiseDistance(
            pairwiseSeparationSum, n);

        return result;
    }

private:
    static double clampCursor(double cursor) {
        if (cursor <= 0.0)
            return 0.0;
        if (cursor >= 1.0)
            return 1.0;
        return cursor;
    }

    static bool validEventTime(double time) {
        return std::isfinite(time) && time >= 0.0 && time <= 1.0;
    }

    static bool crossed(double eventTime, double oldCursor,
                        double newCursor, bool upward) {
        return upward
            ? oldCursor < eventTime && eventTime <= newCursor
            : newCursor < eventTime && eventTime <= oldCursor;
    }

    static int safeLeafCount(const Tree& tree) {
        return tree.leafCount >= 2 && tree.leafCount <= MAX_LEAVES
            ? tree.leafCount
            : 0;
    }

    static int safeNodeCount(const Tree& tree) {
        if (tree.nodeCount < 0)
            return 0;
        return tree.nodeCount > MAX_NODES ? MAX_NODES : tree.nodeCount;
    }

    static int safeMutationCount(const Tree& tree) {
        if (tree.mutationCount < 0)
            return 0;
        return tree.mutationCount > MAX_MUTATIONS
            ? MAX_MUTATIONS
            : tree.mutationCount;
    }

    static int popcount(std::uint16_t value) {
        int count = 0;
        while (value != 0u) {
            value = static_cast<std::uint16_t>(value & (value - 1u));
            ++count;
        }
        return count;
    }

    static double pairwiseDistance(int separationSum, int leafCount) {
        const int pairCount = leafCount * (leafCount - 1) / 2;
        return pairCount > 0
            ? static_cast<double>(separationSum) /
                  static_cast<double>(pairCount)
            : 0.0;
    }

    // Each row is the complete mutation state for a boundary index: row i
    // contains mutations [i, mutationCount). The tree is immutable during
    // playback, so hostile full-range scrubs can copy one fixed row instead of
    // applying as many as 128 mutations to 16 leaves in each direction.
    void rebuildMutationSuffixes(const Tree& tree) {
        const int n = safeLeafCount(tree);
        const int mutationCount = safeMutationCount(tree);
        mutationSuffixTraits_[mutationCount].fill(0);
        mutationSuffixSeparation_[mutationCount] = 0;

        const std::uint16_t allLeaves = leafMask(n);
        for (int i = mutationCount - 1; i >= 0; --i) {
            mutationSuffixTraits_[i] = mutationSuffixTraits_[i + 1];
            mutationSuffixSeparation_[i] =
                mutationSuffixSeparation_[i + 1];

            const Mutation& mutation = tree.mutations[i];
            if (!validEventTime(mutation.normTime) ||
                (mutation.sign != -1 && mutation.sign != 1))
                continue;

            const std::uint16_t mask = static_cast<std::uint16_t>(
                mutation.descendantMask & allLeaves);
            for (int leaf = 0; leaf < n; ++leaf) {
                if ((mask & static_cast<std::uint16_t>(1u << leaf)) != 0u)
                    mutationSuffixTraits_[i][leaf] += mutation.sign;
            }
            const int descendants = popcount(mask);
            mutationSuffixSeparation_[i] += descendants * (n - descendants);
        }
    }

    void copyMutationSuffix(int leafCount, int mutationCount) {
        int boundary = nextMutation_;
        if (boundary < 0)
            boundary = 0;
        if (boundary > mutationCount)
            boundary = mutationCount;
        values_.traitSteps = mutationSuffixTraits_[boundary];
        pairwiseSeparationSum_ = mutationSuffixSeparation_[boundary];
        values_.meanPairwiseMutationDistance = pairwiseDistance(
            pairwiseSeparationSum_, leafCount);
    }

    static int nodeBoundary(const Tree& tree, double cursor) {
        int first = safeLeafCount(tree);
        int last = safeNodeCount(tree);
        while (first < last) {
            const int middle = first + (last - first) / 2;
            if (tree.nodes[middle].normTime <= cursor)
                first = middle + 1;
            else
                last = middle;
        }
        return first;
    }

    static int mutationBoundary(const Tree& tree, double cursor) {
        int first = 0;
        int last = safeMutationCount(tree);
        while (first < last) {
            const int middle = first + (last - first) / 2;
            if (tree.mutations[middle].normTime <= cursor)
                first = middle + 1;
            else
                last = middle;
        }
        return first;
    }

    // Valid installed trees guarantee chronological internal-node order and
    // sorted mutation order. Boundary lookup is fixed and logarithmic;
    // installation, cable transitions, load, and loop teleports still go
    // through canonical state.
    void canonicalize(const Tree& tree) {
        values_ = canonicalValues(tree, cursor_);

        const int n = safeLeafCount(tree);
        const int mutationCount = safeMutationCount(tree);
        nextNode_ = nodeBoundary(tree, cursor_);
        nextMutation_ = mutationBoundary(tree, cursor_);

        copyMutationSuffix(n, mutationCount);
    }

    EventFacts traverseTo(const Tree& tree, double targetCursor) {
        EventFacts facts;
        const double target = clampCursor(targetCursor);
        if (target == cursor_)
            return facts;

        const int n = safeLeafCount(tree);
        const int mutationCount = safeMutationCount(tree);
        const bool upward = target > cursor_;
        const int oldNodeBoundary = nextNode_;
        const int newNodeBoundary = nodeBoundary(tree, target);
        const bool activeLineagesChanged =
            newNodeBoundary != oldNodeBoundary;
        if (activeLineagesChanged) {
            facts.nodeCrossed = true;
            if (upward && oldNodeBoundary <= tree.root &&
                tree.root < newNodeBoundary)
                facts.arrivedAtMrca = true;
            if (!upward && newNodeBoundary <= tree.root &&
                tree.root < oldNodeBoundary)
                facts.departedMrca = true;
            nextNode_ = newNodeBoundary;
        }

        const int newMutationBoundary = mutationBoundary(tree, target);
        const bool mutationBoundaryChanged =
            newMutationBoundary != nextMutation_;
        if (mutationBoundaryChanged) {
            facts.mutationCrossed = true;
            nextMutation_ = newMutationBoundary;
        }

        cursor_ = target;
        values_.lineageCount = n - (nextNode_ - n);
        if (values_.lineageCount < 1)
            values_.lineageCount = 1;
        if (values_.lineageCount > n)
            values_.lineageCount = n;
        if (activeLineagesChanged)
            rebuildActiveLineages(tree, cursor_, values_);
        if (mutationBoundaryChanged)
            copyMutationSuffix(n, mutationCount);
        return facts;
    }

    static std::uint16_t leafMask(int leafCount) {
        if (leafCount <= 0)
            return 0u;
        return leafCount == MAX_LEAVES
            ? static_cast<std::uint16_t>(0xffffu)
            : static_cast<std::uint16_t>((1u << leafCount) - 1u);
    }

    // A lineage crosses the cursor when its node is no younger than the
    // cursor and its parent is strictly older. At equality the newly formed
    // internal node replaces its two child lineages.
    static void rebuildActiveLineages(const Tree& tree, double age,
                                      PlaybackValues& values) {
        values.activeLineageMasks.fill(0u);
        values.activeLineageCount = 0;
        const int nodeCount = safeNodeCount(tree);
        for (int i = 0; i < nodeCount &&
                        values.activeLineageCount < MAX_LEAVES; ++i) {
            const Node& node = tree.nodes[i];
            if (!validEventTime(node.normTime) || node.normTime > age)
                continue;

            bool parentIsOlder = node.parent < 0;
            if (node.parent >= 0 && node.parent < nodeCount) {
                const double parentTime = tree.nodes[node.parent].normTime;
                parentIsOlder = validEventTime(parentTime) && parentTime > age;
            }
            if (!parentIsOlder)
                continue;

            values.activeLineageMasks[values.activeLineageCount++] =
                node.descendantMask;
        }
        sortMasksByFirstLeaf(values.activeLineageMasks,
                             values.activeLineageCount);
    }

    static int firstLeaf(std::uint16_t mask) {
        for (int leaf = 0; leaf < MAX_LEAVES; ++leaf) {
            if ((mask & static_cast<std::uint16_t>(1u << leaf)) != 0u)
                return leaf;
        }
        return MAX_LEAVES;
    }

    static void sortMasksByFirstLeaf(
        std::array<std::uint16_t, MAX_LEAVES>& masks, int count) {
        for (int i = 1; i < count; ++i) {
            const std::uint16_t value = masks[i];
            const int key = firstLeaf(value);
            int j = i;
            while (j > 0 && firstLeaf(masks[j - 1]) > key) {
                masks[j] = masks[j - 1];
                --j;
            }
            masks[j] = value;
        }
    }

    double sourceCursor() const {
        return direction_ == Direction::Ancestry ? 0.0 : 1.0;
    }

    bool canAdvanceFromCursor() const {
        if (loop_)
            return true;
        return direction_ == Direction::Ancestry
            ? cursor_ < 1.0
            : cursor_ > 0.0;
    }

    double cursor_ = 1.0;
    Direction direction_ = Direction::Descent;
    bool loop_ = true;
    bool running_ = true;
    PlaybackValues values_{};
    int nextNode_ = 0;
    int nextMutation_ = 0;
    int pairwiseSeparationSum_ = 0;
    std::array<std::array<int, MAX_LEAVES>, MAX_MUTATIONS + 1>
        mutationSuffixTraits_{};
    std::array<int, MAX_MUTATIONS + 1> mutationSuffixSeparation_{};
};

} // namespace lineages
} // namespace coalescent
