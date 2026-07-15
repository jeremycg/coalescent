#pragma once

#include "random.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace coalescent {
namespace lineages {

static const int MAX_LEAVES = 16;
static const int MAX_NODES = 2 * MAX_LEAVES - 1;
static const int MAX_MUTATIONS = 128;

inline double mutationLambdaFor(double mutate) {
    if (!std::isfinite(mutate) || mutate < 0.0 || mutate > 1.0)
        return std::numeric_limits<double>::quiet_NaN();
    return 64.0 * mutate * mutate;
}

struct Node {
    int left;
    int right;
    int parent;
    double rawTime;
    double normTime;
    std::uint16_t descendantMask;
    double displayX;

    Node()
        : left(-1), right(-1), parent(-1), rawTime(0.0), normTime(0.0),
          descendantMask(0u), displayX(0.0) {}
};

struct Mutation {
    int branchChild;
    double rawTime;
    double normTime;
    std::uint16_t descendantMask;
    std::int8_t sign;

    Mutation()
        : branchChild(-1), rawTime(0.0), normTime(0.0),
          descendantMask(0u), sign(1) {}
};

struct Tree {
    int leafCount;
    int nodeCount;
    int root;
    std::array<Node, MAX_NODES> nodes;

    // Maps final planar output channel -> original leaf node index.
    std::array<int, MAX_LEAVES> leafOrder;

    int mutationCount;
    std::array<Mutation, MAX_MUTATIONS> mutations;
    double totalRawBranchLength;
    double generatedMutate;
    double mutationLambda;
    bool mutationsTruncated;

    Tree()
        : leafCount(0), nodeCount(0), root(-1), mutationCount(0),
          totalRawBranchLength(0.0), generatedMutate(0.0),
          mutationLambda(0.0), mutationsTruncated(false) {
        leafOrder.fill(-1);
    }

    // Validates canonical topology/raw event data and atomically reconstructs
    // every redundant field. Serialized parent links, masks, normalized times,
    // leaf order, display positions, and mutation masks are never trusted.
    bool validateAndDerive() {
        Tree candidate = *this;
        if (!candidate.deriveInPlace())
            return false;
        *this = candidate;
        return true;
    }

    // Tests whether the canonical data can form a valid tree. Derived fields
    // are deliberately ignored because restore rebuilds them.
    bool validate() const {
        Tree candidate = *this;
        return candidate.deriveInPlace();
    }

private:
    struct MutationLess {
        bool operator()(const Mutation& a, const Mutation& b) const {
            if (a.normTime != b.normTime)
                return a.normTime < b.normTime;
            if (a.rawTime != b.rawTime)
                return a.rawTime < b.rawTime;
            if (a.branchChild != b.branchChild)
                return a.branchChild < b.branchChild;
            return a.sign < b.sign;
        }
    };

    bool deriveNode(int index, int& nextChannel) {
        Node& node = nodes[index];
        if (index < leafCount) {
            if (nextChannel >= leafCount)
                return false;
            leafOrder[nextChannel] = index;
            node.descendantMask = static_cast<std::uint16_t>(1u << nextChannel);
            node.displayX = leafCount > 1
                ? static_cast<double>(nextChannel) / static_cast<double>(leafCount - 1)
                : 0.5;
            ++nextChannel;
            return true;
        }

        if (!deriveNode(node.left, nextChannel) || !deriveNode(node.right, nextChannel))
            return false;
        const std::uint16_t leftMask = nodes[node.left].descendantMask;
        const std::uint16_t rightMask = nodes[node.right].descendantMask;
        if (leftMask == 0u || rightMask == 0u || (leftMask & rightMask) != 0u)
            return false;
        node.descendantMask = static_cast<std::uint16_t>(leftMask | rightMask);
        node.displayX = 0.5 * (nodes[node.left].displayX + nodes[node.right].displayX);
        return true;
    }

    bool deriveInPlace() {
        if (leafCount < 2 || leafCount > MAX_LEAVES ||
            nodeCount != 2 * leafCount - 1 || root != nodeCount - 1 ||
            mutationCount < 0 || mutationCount > MAX_MUTATIONS ||
            !std::isfinite(generatedMutate) || generatedMutate < 0.0 || generatedMutate > 1.0)
            return false;

        const double expectedLambda = mutationLambdaFor(generatedMutate);
        mutationLambda = expectedLambda;
        if (mutationsTruncated && mutationCount != MAX_MUTATIONS)
            return false;

        leafOrder.fill(-1);
        for (int index = 0; index < nodeCount; ++index) {
            Node& node = nodes[index];
            node.parent = -1;
            node.normTime = 0.0;
            node.descendantMask = 0u;
            node.displayX = 0.0;

            if (!std::isfinite(node.rawTime))
                return false;
            if (index < leafCount) {
                if (node.left != -1 || node.right != -1 || node.rawTime != 0.0)
                    return false;
            }
            else {
                if (node.left < 0 || node.left >= index ||
                    node.right < 0 || node.right >= index || node.left == node.right ||
                    !(node.rawTime > nodes[node.left].rawTime) ||
                    !(node.rawTime > nodes[node.right].rawTime))
                    return false;
                if (index > leafCount && !(node.rawTime > nodes[index - 1].rawTime))
                    return false;
                if (nodes[node.left].parent != -1 || nodes[node.right].parent != -1)
                    return false;
                nodes[node.left].parent = index;
                nodes[node.right].parent = index;
            }
        }

        if (!(nodes[root].rawTime > 0.0) || !std::isfinite(nodes[root].rawTime) ||
            nodes[root].parent != -1)
            return false;
        for (int index = 0; index < nodeCount - 1; ++index)
            if (nodes[index].parent < 0)
                return false;

        const double rootRawTime = nodes[root].rawTime;
        for (int index = 0; index < nodeCount; ++index) {
            nodes[index].normTime = index < leafCount ? 0.0 : nodes[index].rawTime / rootRawTime;
            if (!std::isfinite(nodes[index].normTime) || nodes[index].normTime < 0.0 ||
                nodes[index].normTime > 1.0)
                return false;
            if (index >= leafCount &&
                !(nodes[index].normTime > nodes[index - 1].normTime))
                return false;
        }
        nodes[root].normTime = 1.0;

        int nextChannel = 0;
        if (!deriveNode(root, nextChannel) || nextChannel != leafCount)
            return false;
        const std::uint16_t expectedRootMask = static_cast<std::uint16_t>(
            (static_cast<std::uint32_t>(1u) << leafCount) - 1u);
        if (nodes[root].descendantMask != expectedRootMask)
            return false;

        totalRawBranchLength = 0.0;
        for (int index = 0; index < nodeCount; ++index) {
            if (index == root)
                continue;
            const double length = nodes[nodes[index].parent].rawTime - nodes[index].rawTime;
            if (!(length > 0.0) || !std::isfinite(length))
                return false;
            totalRawBranchLength += length;
        }
        if (!(totalRawBranchLength > 0.0) || !std::isfinite(totalRawBranchLength))
            return false;

        for (int index = 0; index < mutationCount; ++index) {
            Mutation& mutation = mutations[index];
            if (mutation.branchChild < 0 || mutation.branchChild >= nodeCount ||
                mutation.branchChild == root || !std::isfinite(mutation.rawTime) ||
                (mutation.sign != -1 && mutation.sign != 1))
                return false;
            const Node& child = nodes[mutation.branchChild];
            const Node& parent = nodes[child.parent];
            if (!(mutation.rawTime > child.rawTime) || !(mutation.rawTime < parent.rawTime))
                return false;
            mutation.normTime = mutation.rawTime / rootRawTime;
            if (!(mutation.normTime > child.normTime) ||
                !(mutation.normTime < parent.normTime))
                return false;
            mutation.descendantMask = child.descendantMask;
        }
        std::sort(mutations.begin(), mutations.begin() + mutationCount, MutationLess());
        return true;
    }
};

struct GeneratedState {
    Tree tree;
    Pcg32::State nextRng;

    bool validateAndDerive() {
        GeneratedState candidate = *this;
        if (!Pcg32::validState(candidate.nextRng) || !candidate.tree.validateAndDerive())
            return false;
        *this = candidate;
        return true;
    }
};

class KingmanGenerator {
public:
    KingmanGenerator(std::uint64_t seedValue = 42u, std::uint64_t stream = 54u)
        : rng_(seedValue, stream) {}

    void seed(std::uint64_t seedValue, std::uint64_t stream) {
        rng_.seed(seedValue, stream);
    }

    Pcg32::State rngState() const { return rng_.state(); }

    bool restoreRng(const Pcg32::State& saved) {
        Pcg32 candidate = rng_;
        if (!candidate.restore(saved))
            return false;
        rng_ = candidate;
        return true;
    }

    bool generateTopology(int leafCount, Tree& destination) {
        if (leafCount < 2 || leafCount > MAX_LEAVES)
            return false;

        Pcg32 trialRng = rng_;
        Tree candidate;
        candidate.leafCount = leafCount;
        candidate.nodeCount = 2 * leafCount - 1;
        candidate.root = candidate.nodeCount - 1;

        std::array<int, MAX_LEAVES> active;
        for (int index = 0; index < leafCount; ++index)
            active[index] = index;

        int activeCount = leafCount;
        double rawTime = 0.0;
        int nextNode = leafCount;
        while (activeCount > 1) {
            const double rate = 0.5 * static_cast<double>(activeCount) *
                                static_cast<double>(activeCount - 1);
            const double wait = trialRng.exponential(rate);
            const double nextTime = rawTime + wait;
            if (!(wait > 0.0) || !std::isfinite(wait) ||
                !(nextTime > rawTime) || !std::isfinite(nextTime))
                return false;

            const int firstPosition = static_cast<int>(
                trialRng.bounded(static_cast<std::uint32_t>(activeCount)));
            int secondPosition = static_cast<int>(
                trialRng.bounded(static_cast<std::uint32_t>(activeCount - 1)));
            if (secondPosition >= firstPosition)
                ++secondPosition;

            Node& parent = candidate.nodes[nextNode];
            parent.left = active[firstPosition];
            parent.right = active[secondPosition];
            parent.rawTime = nextTime;

            const int keep = std::min(firstPosition, secondPosition);
            const int remove = std::max(firstPosition, secondPosition);
            active[keep] = nextNode;
            active[remove] = active[activeCount - 1];
            --activeCount;
            ++nextNode;
            rawTime = nextTime;
        }

        if (active[0] != candidate.root || !candidate.validateAndDerive())
            return false;
        rng_ = trialRng;
        destination = candidate;
        return true;
    }

    bool generateMutations(double mutate, Tree& destination) {
        if (!std::isfinite(mutate) || mutate < 0.0 || mutate > 1.0)
            return false;

        Pcg32 trialRng = rng_;
        Tree candidate = destination;
        candidate.mutationCount = 0;
        candidate.generatedMutate = mutate;
        candidate.mutationLambda = mutationLambdaFor(mutate);
        candidate.mutationsTruncated = false;
        if (!candidate.validateAndDerive())
            return false;

        const Pcg32::PoissonResult count = trialRng.poissonCapped(
            candidate.mutationLambda, static_cast<std::uint32_t>(MAX_MUTATIONS));
        candidate.mutationCount = static_cast<int>(count.count);
        candidate.mutationsTruncated = count.truncated;

        std::array<int, MAX_NODES> branchChildren;
        int branchCount = 0;
        for (int index = 0; index < candidate.nodeCount; ++index)
            if (index != candidate.root)
                branchChildren[branchCount++] = index;

        for (int index = 0; index < candidate.mutationCount; ++index) {
            double distance = trialRng.uniformOpen() * candidate.totalRawBranchLength;
            int branchChild = branchChildren[branchCount - 1];
            for (int branch = 0; branch < branchCount; ++branch) {
                const int childIndex = branchChildren[branch];
                const Node& child = candidate.nodes[childIndex];
                const double length = candidate.nodes[child.parent].rawTime - child.rawTime;
                if (distance < length) {
                    branchChild = childIndex;
                    break;
                }
                distance -= length;
            }

            const Node& child = candidate.nodes[branchChild];
            const Node& parent = candidate.nodes[child.parent];
            const double fraction = trialRng.uniformOpen();
            double mutationTime = child.rawTime + fraction * (parent.rawTime - child.rawTime);
            if (!(mutationTime > child.rawTime))
                mutationTime = std::nextafter(child.rawTime, parent.rawTime);
            if (!(mutationTime < parent.rawTime))
                mutationTime = std::nextafter(parent.rawTime, child.rawTime);
            if (!(mutationTime > child.rawTime) || !(mutationTime < parent.rawTime))
                return false;

            Mutation& mutation = candidate.mutations[index];
            mutation.branchChild = branchChild;
            mutation.rawTime = mutationTime;
            mutation.sign = (trialRng.next() & 1u) != 0u ? 1 : -1;
        }

        if (!candidate.validateAndDerive())
            return false;
        rng_ = trialRng;
        destination = candidate;
        return true;
    }

    // Full NEW generation is atomic for both installed Tree and future RNG.
    bool generate(int leafCount, double mutate, Tree& destination) {
        KingmanGenerator trial = *this;
        Tree candidate;
        if (!trial.generateTopology(leafCount, candidate) ||
            !trial.generateMutations(mutate, candidate))
            return false;
        *this = trial;
        destination = candidate;
        return true;
    }

    GeneratedState capture(const Tree& tree) const {
        GeneratedState saved;
        saved.tree = tree;
        saved.nextRng = rng_.state();
        return saved;
    }

    // Validates the tree and post-generation RNG together. The generator is
    // installed only after callers accept the returned, derived tree copy.
    bool restore(const GeneratedState& saved, Tree& destination) {
        GeneratedState candidate = saved;
        if (!candidate.validateAndDerive())
            return false;
        Pcg32 candidateRng = rng_;
        if (!candidateRng.restore(candidate.nextRng))
            return false;
        rng_ = candidateRng;
        destination = candidate.tree;
        return true;
    }

private:
    Pcg32 rng_;
};

} // namespace lineages
} // namespace coalescent
