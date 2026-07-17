// SDK-free structural, statistical, and playback contract for Lineages.
//
//   g++ -std=c++11 -O2 tools/stability/lineages.cpp -o /tmp/lineages && /tmp/lineages
#include "../../src/lineages/random.hpp"
#include "../../src/lineages/kingman.hpp"
#include "../../src/lineages/playback.hpp"
#include "../../src/lineages/state.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

using coalescent::lineages::Direction;
using coalescent::lineages::EventFacts;
using coalescent::lineages::GeneratedState;
using coalescent::lineages::Hex64Codec;
using coalescent::lineages::KingmanGenerator;
using coalescent::lineages::MAX_LEAVES;
using coalescent::lineages::MAX_MUTATIONS;
using coalescent::lineages::MAX_NODES;
using coalescent::lineages::Mutation;
using coalescent::lineages::Node;
using coalescent::lineages::Pcg32;
using coalescent::lineages::Playback;
using coalescent::lineages::PlaybackValues;
using coalescent::lineages::Snapshot;
using coalescent::lineages::Tree;
using coalescent::lineages::makeFactorySnapshot;

static int failures = 0;

static void fail(const char* format, ...) {
    std::printf("FAIL: ");
    va_list args;
    va_start(args, format);
    std::vprintf(format, args);
    va_end(args);
    std::printf("\n");
    ++failures;
}

static int popcount(std::uint16_t value) {
    int result = 0;
    while (value != 0u) {
        value = static_cast<std::uint16_t>(value & (value - 1u));
        ++result;
    }
    return result;
}

static std::uint16_t allLeavesMask(int n) {
    return n == MAX_LEAVES
        ? static_cast<std::uint16_t>(0xffffu)
        : static_cast<std::uint16_t>((1u << n) - 1u);
}

static bool sameDoubleBits(double a, double b) {
    std::uint64_t aBits = 0u;
    std::uint64_t bBits = 0u;
    std::memcpy(&aBits, &a, sizeof(aBits));
    std::memcpy(&bBits, &b, sizeof(bBits));
    return aBits == bBits;
}

static bool sameFloatBits(float a, float b) {
    std::uint32_t aBits = 0u;
    std::uint32_t bBits = 0u;
    std::memcpy(&aBits, &a, sizeof(aBits));
    std::memcpy(&bBits, &b, sizeof(bBits));
    return aBits == bBits;
}

static bool sameNode(const Node& a, const Node& b) {
    return a.left == b.left && a.right == b.right && a.parent == b.parent
        && sameDoubleBits(a.rawTime, b.rawTime)
        && sameDoubleBits(a.normTime, b.normTime)
        && a.descendantMask == b.descendantMask
        && sameDoubleBits(a.displayX, b.displayX);
}

static bool sameMutation(const Mutation& a, const Mutation& b) {
    return a.branchChild == b.branchChild
        && sameDoubleBits(a.rawTime, b.rawTime)
        && sameDoubleBits(a.normTime, b.normTime)
        && a.descendantMask == b.descendantMask && a.sign == b.sign;
}

static bool sameTree(const Tree& a, const Tree& b) {
    if (a.leafCount != b.leafCount || a.nodeCount != b.nodeCount
        || a.root != b.root || a.mutationCount != b.mutationCount
        || !sameDoubleBits(a.totalRawBranchLength, b.totalRawBranchLength)
        || !sameDoubleBits(a.generatedMutate, b.generatedMutate)
        || !sameDoubleBits(a.mutationLambda, b.mutationLambda)
        || a.mutationsTruncated != b.mutationsTruncated)
        return false;
    for (int i = 0; i < MAX_NODES; ++i) {
        if (!sameNode(a.nodes[i], b.nodes[i]))
            return false;
    }
    for (int i = 0; i < MAX_LEAVES; ++i) {
        if (a.leafOrder[i] != b.leafOrder[i])
            return false;
    }
    for (int i = 0; i < MAX_MUTATIONS; ++i) {
        if (!sameMutation(a.mutations[i], b.mutations[i]))
            return false;
    }
    return true;
}

static bool sameRng(const Pcg32::State& a, const Pcg32::State& b) {
    return a.state == b.state && a.increment == b.increment;
}

static bool near(double observed, double expected, double tolerance = 1e-11) {
    return std::fabs(observed - expected) <= tolerance;
}

static bool contiguousMask(std::uint16_t mask) {
    if (mask == 0u)
        return false;
    int first = 0;
    while ((mask & static_cast<std::uint16_t>(1u << first)) == 0u)
        ++first;
    const int width = popcount(mask);
    const std::uint32_t expected = ((UINT32_C(1) << width) - 1u) << first;
    return mask == static_cast<std::uint16_t>(expected);
}

static bool finitePlaybackValues(const PlaybackValues& values) {
    if (values.sampleCount < 0 || values.sampleCount > MAX_LEAVES
        || values.lineageCount < 0 || values.lineageCount > MAX_LEAVES
        || values.activeLineageCount < 0
        || values.activeLineageCount > MAX_LEAVES
        || !std::isfinite(values.meanPairwiseMutationDistance)
        || values.meanPairwiseMutationDistance < 0.0)
        return false;
    for (int i = 0; i < MAX_LEAVES; ++i) {
        if (values.traitSteps[i] < -MAX_MUTATIONS
            || values.traitSteps[i] > MAX_MUTATIONS)
            return false;
    }
    return true;
}

static bool factsEmpty(const EventFacts& facts);

static void testRandom() {
    Pcg32 referenceRng(42u, 54u);
    const std::uint32_t reference[] = {
        UINT32_C(0xa15c02b7), UINT32_C(0x7b47f409), UINT32_C(0xba1d3330),
        UINT32_C(0x83d2f293), UINT32_C(0xbfa4784b), UINT32_C(0xcbed606e)
    };
    for (std::size_t i = 0; i < sizeof(reference) / sizeof(reference[0]); ++i) {
        const std::uint32_t observed = referenceRng.next();
        if (observed != reference[i])
            fail("PCG reference word %zu was %08x, expected %08x",
                 i, observed, reference[i]);
    }

    Pcg32 resumed(123u, 456u);
    for (int i = 0; i < 19; ++i)
        resumed.next();
    const Pcg32::State saved = resumed.state();
    std::uint32_t future[32];
    for (int i = 0; i < 32; ++i)
        future[i] = resumed.next();
    Pcg32 restored;
    if (!restored.restore(saved))
        fail("valid PCG state was rejected");
    for (int i = 0; i < 32; ++i) {
        if (restored.next() != future[i]) {
            fail("restored PCG future diverged at word %d", i);
            break;
        }
    }
    Pcg32::State invalid = saved;
    invalid.increment &= ~UINT64_C(1);
    const Pcg32::State beforeInvalid = restored.state();
    if (restored.restore(invalid)
        || restored.state().state != beforeInvalid.state
        || restored.state().increment != beforeInvalid.increment)
        fail("invalid even PCG increment was not rejected atomically");

    Pcg32 openRng(901u, 37u);
    for (int i = 0; i < 100000; ++i) {
        const double open = openRng.uniformOpen();
        const double closedOpen = openRng.uniformClosedOpen();
        if (!(open > 0.0 && open < 1.0)
            || !(closedOpen >= 0.0 && closedOpen < 1.0)) {
            fail("uniform endpoint contract failed");
            break;
        }
    }

    const double lambdas[] = {0.25, 4.0, 16.0, 64.0};
    const int draws = 60000;
    for (std::size_t li = 0; li < sizeof(lambdas) / sizeof(lambdas[0]); ++li) {
        const double lambda = lambdas[li];
        Pcg32 poissonRng(UINT64_C(0x12340000) + li, UINT64_C(0x56780000) + li);
        double sum = 0.0;
        double sumSquares = 0.0;
        int truncated = 0;
        for (int i = 0; i < draws; ++i) {
            const Pcg32::PoissonResult result = poissonRng.poissonCapped(lambda, 128u);
            sum += result.count;
            sumSquares += static_cast<double>(result.count) * result.count;
            truncated += result.truncated ? 1 : 0;
        }
        const double mean = sum / draws;
        const double variance = sumSquares / draws - mean * mean;
        const double meanTolerance = 7.0 * std::sqrt(lambda / draws) + 0.002;
        const double varianceTolerance = 0.035 * lambda + 0.025;
        std::printf("poisson lambda=%5.2f mean=%8.5f variance=%8.5f\n",
                    lambda, mean, variance);
        if (std::fabs(mean - lambda) > meanTolerance)
            fail("Poisson mean at lambda %.2f (error %.6f, tolerance %.6f)",
                 lambda, mean - lambda, meanTolerance);
        if (std::fabs(variance - lambda) > varianceTolerance)
            fail("Poisson variance at lambda %.2f (error %.6f, tolerance %.6f)",
                 lambda, variance - lambda, varianceTolerance);
        if (truncated != 0)
            fail("ordinary lambda %.2f unexpectedly reached mutation cap", lambda);
    }

    Pcg32 hostile(44u, 99u);
    if (hostile.poissonCapped(0.0, 128u).count != 0u
        || hostile.poissonCapped(-1.0, 128u).count != 0u
        || hostile.poissonCapped(NAN, 128u).count != 0u
        || hostile.poissonCapped(1.0, 0u).count != 0u)
        fail("hostile Poisson arguments did not return bounded zero");
    const Pcg32::PoissonResult capped = hostile.poissonCapped(1000000.0, 128u);
    if (capped.count != 128u || !capped.truncated)
        fail("large finite Poisson draw did not truncate at the hard cap");
}

static bool structuralInvariants(const Tree& tree, int expectedLeaves) {
    if (!tree.validate() || tree.leafCount != expectedLeaves
        || tree.nodeCount != 2 * expectedLeaves - 1
        || tree.root != tree.nodeCount - 1
        || tree.root < 0 || tree.root >= MAX_NODES
        || tree.nodes[tree.root].parent != -1
        || tree.nodes[tree.root].normTime != 1.0
        || tree.nodes[tree.root].descendantMask != allLeavesMask(expectedLeaves)
        || !(tree.totalRawBranchLength > 0.0)
        || !std::isfinite(tree.totalRawBranchLength))
        return false;

    bool originalSeen[MAX_LEAVES] = {};
    for (int channel = 0; channel < expectedLeaves; ++channel) {
        const int leaf = tree.leafOrder[channel];
        if (leaf < 0 || leaf >= expectedLeaves || originalSeen[leaf])
            return false;
        originalSeen[leaf] = true;
        const Node& node = tree.nodes[leaf];
        if (node.left != -1 || node.right != -1 || node.parent < expectedLeaves
            || node.parent >= tree.nodeCount || node.rawTime != 0.0
            || node.normTime != 0.0
            || node.descendantMask != static_cast<std::uint16_t>(1u << channel)
            || !std::isfinite(node.displayX) || node.displayX < 0.0
            || node.displayX > 1.0)
            return false;
        const double expectedX = static_cast<double>(channel)
                               / static_cast<double>(expectedLeaves - 1);
        if (!near(node.displayX, expectedX, 2e-15))
            return false;
    }

    double branchLength = 0.0;
    for (int i = 0; i < tree.nodeCount; ++i) {
        const Node& node = tree.nodes[i];
        if (!std::isfinite(node.rawTime) || !std::isfinite(node.normTime)
            || !std::isfinite(node.displayX) || node.rawTime < 0.0
            || node.normTime < 0.0 || node.normTime > 1.0
            || node.displayX < 0.0 || node.displayX > 1.0
            || !contiguousMask(node.descendantMask))
            return false;
        if (std::fabs(node.normTime
                      - node.rawTime / tree.nodes[tree.root].rawTime) > 2e-14)
            return false;
        if (i >= expectedLeaves) {
            if (i > expectedLeaves
                && !(node.normTime > tree.nodes[i - 1].normTime))
                return false;
            if (node.left < 0 || node.left >= i || node.right < 0
                || node.right >= i || node.left == node.right)
                return false;
            const Node& left = tree.nodes[node.left];
            const Node& right = tree.nodes[node.right];
            if (left.parent != i || right.parent != i
                || !(node.rawTime > left.rawTime)
                || !(node.rawTime > right.rawTime)
                || !(node.normTime > left.normTime)
                || !(node.normTime > right.normTime)
                || (left.descendantMask & right.descendantMask) != 0u
                || node.descendantMask
                    != static_cast<std::uint16_t>(left.descendantMask
                                                  | right.descendantMask))
                return false;
        }
        if (i != tree.root) {
            if (node.parent < expectedLeaves || node.parent >= tree.nodeCount)
                return false;
            branchLength += tree.nodes[node.parent].rawTime - node.rawTime;
        }
    }
    if (std::fabs(branchLength - tree.totalRawBranchLength)
        > 2e-13 * std::max(1.0, branchLength))
        return false;

    for (int leaf = 0; leaf < expectedLeaves; ++leaf) {
        int node = leaf;
        int hops = 0;
        while (node != tree.root && hops <= tree.nodeCount) {
            node = tree.nodes[node].parent;
            ++hops;
        }
        if (node != tree.root || hops > tree.nodeCount)
            return false;
    }

    if (tree.mutationCount < 0 || tree.mutationCount > MAX_MUTATIONS
        || !std::isfinite(tree.generatedMutate)
        || tree.generatedMutate < 0.0 || tree.generatedMutate > 1.0
        || !std::isfinite(tree.mutationLambda)
        || std::fabs(tree.mutationLambda
                     - 64.0 * tree.generatedMutate * tree.generatedMutate) > 1e-13)
        return false;
    for (int i = 0; i < tree.mutationCount; ++i) {
        const Mutation& mutation = tree.mutations[i];
        if (mutation.branchChild < 0 || mutation.branchChild >= tree.nodeCount
            || mutation.branchChild == tree.root
            || (mutation.sign != -1 && mutation.sign != 1)
            || !std::isfinite(mutation.rawTime)
            || !std::isfinite(mutation.normTime))
            return false;
        const Node& child = tree.nodes[mutation.branchChild];
        const Node& parent = tree.nodes[child.parent];
        if (!(mutation.rawTime > child.rawTime)
            || !(mutation.rawTime < parent.rawTime)
            || !(mutation.normTime > child.normTime)
            || !(mutation.normTime < parent.normTime)
            || mutation.descendantMask != child.descendantMask
            || std::fabs(mutation.normTime
                         - mutation.rawTime / tree.nodes[tree.root].rawTime) > 2e-14
            || (i > 0 && tree.mutations[i - 1].normTime > mutation.normTime))
            return false;
    }
    return true;
}

static void diagnoseStructural(const Tree& tree, int expectedLeaves) {
    const bool basicShape = tree.leafCount >= 2 && tree.leafCount <= MAX_LEAVES
        && tree.nodeCount == 2 * tree.leafCount - 1
        && tree.nodeCount <= MAX_NODES
        && tree.root >= 0 && tree.root < tree.nodeCount;
    if (!basicShape) {
        std::printf("  diagnostic: validate=%d leaves=%d/%d nodes=%d root=%d "
                    "(invalid basic shape)\n",
                    tree.validate() ? 1 : 0, tree.leafCount, expectedLeaves,
                    tree.nodeCount, tree.root);
        return;
    }

    double maxNormError = 0.0;
    double branchLength = 0.0;
    const double rootRawTime = tree.nodes[tree.root].rawTime;
    for (int i = 0; i < tree.nodeCount; ++i) {
        if (std::isfinite(rootRawTime) && rootRawTime != 0.0) {
            maxNormError = std::max(maxNormError, std::fabs(
                tree.nodes[i].normTime - tree.nodes[i].rawTime / rootRawTime));
        }
        if (i != tree.root) {
            const int parent = tree.nodes[i].parent;
            if (parent >= 0 && parent < tree.nodeCount) {
                branchLength += tree.nodes[parent].rawTime
                              - tree.nodes[i].rawTime;
            }
        }
    }
    std::printf("  diagnostic: validate=%d leaves=%d/%d nodes=%d root=%d "
                "rootNorm=%.17g branchErr=%.17g maxNormErr=%.17g\n",
                tree.validate() ? 1 : 0, tree.leafCount, expectedLeaves,
                tree.nodeCount, tree.root, tree.nodes[tree.root].normTime,
                branchLength - tree.totalRawBranchLength, maxNormError);
    for (int channel = 0; channel < tree.leafCount; ++channel) {
        const int leaf = tree.leafOrder[channel];
        if (leaf < 0 || leaf >= tree.leafCount) {
            std::printf("    ch%d invalid leaf=%d\n", channel, leaf);
            continue;
        }
        std::printf("    ch%d leaf%d mask=%04x x=%.17g expectedX=%.17g\n",
                    channel, leaf,
                    static_cast<unsigned>(tree.nodes[leaf].descendantMask),
                    tree.nodes[leaf].displayX,
                    static_cast<double>(channel) / (tree.leafCount - 1));
    }
}

static void testStructuralAndPersistence() {
    int generated = 0;
    for (int n = 2; n <= MAX_LEAVES; ++n) {
        for (int seed = 0; seed < 160; ++seed) {
            KingmanGenerator generator(UINT64_C(0x10000000) + seed,
                                        UINT64_C(0xabc000) + n * 257 + seed);
            Tree tree;
            const double mutate = static_cast<double>(seed % 9) / 8.0;
            if (!generator.generate(n, mutate, tree)
                || !structuralInvariants(tree, n)) {
                fail("tree invariant at n=%d seed=%d", n, seed);
                diagnoseStructural(tree, n);
                return;
            }
            ++generated;
        }
    }
    std::printf("structure: %d trees across n=2..16 valid and planar\n", generated);

    KingmanGenerator transactional(12345u, 67890u);
    Tree output;
    if (!transactional.generate(8, 0.4, output)) {
        fail("transactional test fixture generation failed");
        return;
    }
    const Tree baseline = output;
    const Pcg32::State rngBaseline = transactional.rngState();
    const double invalidMutate[] = {
        -0.001, 1.001, std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()
    };
    for (std::size_t i = 0; i < sizeof(invalidMutate) / sizeof(invalidMutate[0]); ++i) {
        if (transactional.generate(8, invalidMutate[i], output)
            || !sameTree(output, baseline)
            || !sameRng(transactional.rngState(), rngBaseline))
            fail("invalid mutate input %zu was not transactional", i);
    }
    if (transactional.generate(1, 0.4, output)
        || transactional.generate(MAX_LEAVES + 1, 0.4, output)
        || !sameTree(output, baseline)
        || !sameRng(transactional.rngState(), rngBaseline))
        fail("invalid sample count was not transactional");
    if (transactional.generateTopology(1, output)
        || transactional.generateTopology(MAX_LEAVES + 1, output)
        || !sameTree(output, baseline)
        || !sameRng(transactional.rngState(), rngBaseline))
        fail("invalid topology generation was not transactional");
    if (transactional.generateMutations(-1.0, output)
        || transactional.generateMutations(NAN, output)
        || !sameTree(output, baseline)
        || !sameRng(transactional.rngState(), rngBaseline))
        fail("invalid mutation generation was not transactional");

    Tree malformed = baseline;
    malformed.nodes[malformed.root].left = MAX_NODES;
    const Tree malformedBefore = malformed;
    if (malformed.validateAndDerive() || !sameTree(malformed, malformedBefore))
        fail("invalid derive was not rejected transactionally");

    // Patch restore trusts only canonical topology/raw event data. Every cached
    // parent, normalized time, mask, display coordinate, and aggregate must be
    // rebuilt rather than becoming a second competing source of truth.
    Tree poisoned = baseline;
    for (int i = 0; i < poisoned.nodeCount; ++i) {
        poisoned.nodes[i].parent = (i + 7) % poisoned.nodeCount;
        poisoned.nodes[i].normTime = NAN;
        poisoned.nodes[i].descendantMask = 0u;
        poisoned.nodes[i].displayX = -123.0;
    }
    poisoned.leafOrder.fill(-1);
    poisoned.totalRawBranchLength = NAN;
    poisoned.mutationLambda = NAN;
    for (int i = 0; i < poisoned.mutationCount; ++i) {
        poisoned.mutations[i].normTime = NAN;
        poisoned.mutations[i].descendantMask = 0u;
    }
    if (!poisoned.validateAndDerive() || !sameTree(poisoned, baseline))
        fail("restore did not canonically rebuild poisoned redundant fields");

    const GeneratedState saved = transactional.capture(baseline);
    if (!saved.tree.validate() || !Pcg32::validState(saved.nextRng))
        fail("captured generated state was invalid");
    Tree expectedNext;
    if (!transactional.generate(11, 0.61, expectedNext))
        fail("next-tree persistence fixture failed");
    KingmanGenerator restored(1u, 2u);
    Tree installed;
    if (!restored.restore(saved, installed) || !sameTree(installed, saved.tree))
        fail("valid generated state restore failed");
    Tree observedNext;
    if (!restored.generate(11, 0.61, observedNext)
        || !sameTree(expectedNext, observedNext)
        || !sameRng(transactional.rngState(), restored.rngState()))
        fail("restored RNG did not reproduce the next NEW tree");

    GeneratedState invalidSaved = saved;
    invalidSaved.nextRng.increment &= ~UINT64_C(1);
    const Pcg32::State beforeBadRestore = restored.rngState();
    const Tree beforeBadInstalled = installed;
    if (restored.restore(invalidSaved, installed)
        || !sameRng(beforeBadRestore, restored.rngState())
        || !sameTree(beforeBadInstalled, installed))
        fail("invalid generated-state restore changed RNG future");
}

static void testKingmanStatistics() {
    const int n = MAX_LEAVES;
    const int draws = 24000;
    double waitingSum[MAX_LEAVES + 1] = {};
    double mrcaSum = 0.0;
    double branchSum = 0.0;
    KingmanGenerator generator(UINT64_C(0xd00d1234), UINT64_C(0x9e3779b9));
    for (int draw = 0; draw < draws; ++draw) {
        Tree tree;
        if (!generator.generateTopology(n, tree)) {
            fail("Kingman statistics topology generation failed at draw %d", draw);
            return;
        }
        double previous = 0.0;
        for (int index = n; index < tree.nodeCount; ++index) {
            const int active = n - (index - n);
            waitingSum[active] += tree.nodes[index].rawTime - previous;
            previous = tree.nodes[index].rawTime;
        }
        mrcaSum += tree.nodes[tree.root].rawTime;
        branchSum += tree.totalRawBranchLength;
    }

    double worstWaitingRelativeError = 0.0;
    for (int active = 2; active <= n; ++active) {
        const double observed = waitingSum[active] / draws;
        const double expected = 2.0 / (static_cast<double>(active) * (active - 1));
        worstWaitingRelativeError = std::max(
            worstWaitingRelativeError, std::fabs(observed - expected) / expected);
    }
    const double observedMrca = mrcaSum / draws;
    const double expectedMrca = 2.0 * (1.0 - 1.0 / n);
    double harmonic = 0.0;
    for (int i = 1; i < n; ++i)
        harmonic += 1.0 / i;
    const double observedBranch = branchSum / draws;
    const double expectedBranch = 2.0 * harmonic;
    std::printf("Kingman n=%d: wait max-rel=%.4f MRCA=%.5f/%.5f branch=%.5f/%.5f\n",
                n, worstWaitingRelativeError, observedMrca, expectedMrca,
                observedBranch, expectedBranch);
    if (worstWaitingRelativeError > 0.035)
        fail("Kingman waiting-time mean error %.5f exceeded 0.035", worstWaitingRelativeError);
    if (std::fabs(observedMrca - expectedMrca) > 0.03)
        fail("Kingman MRCA mean error %.6f", observedMrca - expectedMrca);
    if (std::fabs(observedBranch - expectedBranch) > 0.07)
        fail("Kingman total branch-length mean error %.6f",
             observedBranch - expectedBranch);

    const int pairN = 8;
    const int pairDraws = 33600;
    int pairCounts[pairN][pairN] = {};
    int channelCounts[pairN][pairN] = {};
    KingmanGenerator pairGenerator(UINT64_C(0x77665544), UINT64_C(0x13579bdf));
    for (int draw = 0; draw < pairDraws; ++draw) {
        Tree tree;
        if (!pairGenerator.generateTopology(pairN, tree)) {
            fail("pair-uniformity topology generation failed");
            return;
        }
        int first = tree.nodes[pairN].left;
        int second = tree.nodes[pairN].right;
        if (first > second)
            std::swap(first, second);
        if (first < 0 || second >= pairN || first == second) {
            fail("first coalescing pair was not two original leaves");
            return;
        }
        ++pairCounts[first][second];
        for (int channel = 0; channel < pairN; ++channel)
            ++channelCounts[channel][tree.leafOrder[channel]];
    }
    const double pairExpected = static_cast<double>(pairDraws)
                              / (pairN * (pairN - 1) / 2);
    const double channelExpected = static_cast<double>(pairDraws) / pairN;
    double worstPairZ = 0.0;
    double worstChannelZ = 0.0;
    for (int first = 0; first < pairN; ++first) {
        for (int second = first + 1; second < pairN; ++second) {
            worstPairZ = std::max(worstPairZ,
                std::fabs(pairCounts[first][second] - pairExpected)
                    / std::sqrt(pairExpected));
        }
        for (int original = 0; original < pairN; ++original) {
            worstChannelZ = std::max(worstChannelZ,
                std::fabs(channelCounts[first][original] - channelExpected)
                    / std::sqrt(channelExpected));
        }
    }
    std::printf("uniformity n=8: first-pair max-z=%.3f channel max-z=%.3f\n",
                worstPairZ, worstChannelZ);
    if (worstPairZ > 6.0)
        fail("first coalescing pair selection is biased (max z %.3f)", worstPairZ);
    if (worstChannelZ > 6.0)
        fail("planar output order has original-leaf bias (max z %.3f)", worstChannelZ);
}

static void testMutations() {
    const int n = 8;
    KingmanGenerator topologyGenerator(UINT64_C(0x31415926), UINT64_C(0x27182818));
    Tree topology;
    if (!topologyGenerator.generateTopology(n, topology)) {
        fail("mutation test topology generation failed");
        return;
    }

    const double lambda = 12.0;
    const double mutate = std::sqrt(lambda / 64.0);
    const int draws = 20000;
    const int bins = 16;
    double countSum = 0.0;
    double countSquares = 0.0;
    long long totalMutations = 0;
    long long positiveSigns = 0;
    long long coordinateBins[bins] = {};
    long long observedBranch[MAX_NODES] = {};
    double expectedBranch[MAX_NODES] = {};
    int truncations = 0;

    KingmanGenerator mutationGenerator(UINT64_C(0xabcdef01), UINT64_C(0x10203040));
    for (int draw = 0; draw < draws; ++draw) {
        Tree tree = topology;
        if (!mutationGenerator.generateMutations(mutate, tree)
            || !structuralInvariants(tree, n)) {
            fail("mutation generation/invariant failure at draw %d", draw);
            return;
        }
        countSum += tree.mutationCount;
        countSquares += static_cast<double>(tree.mutationCount) * tree.mutationCount;
        totalMutations += tree.mutationCount;
        truncations += tree.mutationsTruncated ? 1 : 0;
        for (int child = 0; child < tree.nodeCount; ++child) {
            if (child == tree.root)
                continue;
            const double length = tree.nodes[tree.nodes[child].parent].rawTime
                                - tree.nodes[child].rawTime;
            expectedBranch[child] += tree.mutationCount * length
                                   / tree.totalRawBranchLength;
        }

        for (int m = 0; m < tree.mutationCount; ++m) {
            const Mutation& mutation = tree.mutations[m];
            ++observedBranch[mutation.branchChild];
            if (mutation.sign > 0)
                ++positiveSigns;

            double cumulative = 0.0;
            double coordinate = -1.0;
            for (int child = 0; child < tree.nodeCount; ++child) {
                if (child == tree.root)
                    continue;
                const double childTime = tree.nodes[child].rawTime;
                const double length = tree.nodes[tree.nodes[child].parent].rawTime
                                    - childTime;
                if (child == mutation.branchChild)
                    coordinate = (cumulative + mutation.rawTime - childTime)
                               / tree.totalRawBranchLength;
                cumulative += length;
            }
            if (!(coordinate > 0.0 && coordinate < 1.0)) {
                fail("mutation did not map strictly inside total branch measure");
                return;
            }
            int bin = static_cast<int>(coordinate * bins);
            if (bin < 0)
                bin = 0;
            if (bin >= bins)
                bin = bins - 1;
            ++coordinateBins[bin];
        }
    }

    const double countMean = countSum / draws;
    const double countVariance = countSquares / draws - countMean * countMean;
    const double signExpected = 0.5 * totalMutations;
    const double signZ = std::fabs(positiveSigns - signExpected)
                       / std::sqrt(0.25 * static_cast<double>(totalMutations));
    double worstBranchZ = 0.0;
    for (int child = 0; child < topology.nodeCount; ++child) {
        if (child == topology.root || expectedBranch[child] < 50.0)
            continue;
        const double probability = expectedBranch[child] / totalMutations;
        const double deviation = observedBranch[child] - expectedBranch[child];
        const double sigma = std::sqrt(expectedBranch[child]
                                      * std::max(0.0, 1.0 - probability));
        worstBranchZ = std::max(worstBranchZ, std::fabs(deviation) / sigma);
    }
    const double binExpected = static_cast<double>(totalMutations) / bins;
    double worstCoordinateZ = 0.0;
    for (int bin = 0; bin < bins; ++bin) {
        const double sigma = std::sqrt(binExpected * (1.0 - 1.0 / bins));
        worstCoordinateZ = std::max(worstCoordinateZ,
            std::fabs(coordinateBins[bin] - binExpected) / sigma);
    }
    std::printf("mutations: count %.4f/%.1f variance %.4f signs-z %.3f "
                "branch-z %.3f coordinate-z %.3f\n",
                countMean, lambda, countVariance, signZ,
                worstBranchZ, worstCoordinateZ);
    if (std::fabs(countMean - lambda) > 0.08)
        fail("generated mutation-count mean error %.6f", countMean - lambda);
    if (std::fabs(countVariance - lambda) > 0.25)
        fail("generated mutation-count variance error %.6f", countVariance - lambda);
    if (truncations != 0)
        fail("lambda 12 unexpectedly reached storage cap");
    if (signZ > 6.0)
        fail("mutation signs are imbalanced (z %.3f)", signZ);
    if (worstBranchZ > 6.0)
        fail("mutation branch selection is not length-weighted (z %.3f)", worstBranchZ);
    if (worstCoordinateZ > 6.0)
        fail("mutation positions are not uniform over branch measure (z %.3f)",
             worstCoordinateZ);

    Tree zero = topology;
    if (!mutationGenerator.generateMutations(0.0, zero)
        || zero.mutationCount != 0 || zero.mutationsTruncated
        || zero.mutationLambda != 0.0 || !zero.validate())
        fail("zero MUTATE did not create a valid mutation-free tree");
}

static void printTree(const char* label, const Tree& tree) {
    std::printf("%s: leaves=%d nodes=%d root=%d order=", label,
                tree.leafCount, tree.nodeCount, tree.root);
    for (int i = 0; i < tree.leafCount; ++i)
        std::printf("%s%d", i ? "," : "", tree.leafOrder[i]);
    std::printf("\n");
    for (int i = tree.leafCount; i < tree.nodeCount; ++i) {
        const Node& node = tree.nodes[i];
        std::printf("  N%d=(%d,%d) raw=%.12f norm=%.12f mask=%04x\n",
                    i, node.left, node.right, node.rawTime, node.normTime,
                    static_cast<unsigned>(node.descendantMask));
    }
    for (int i = 0; i < tree.mutationCount; ++i) {
        const Mutation& mutation = tree.mutations[i];
        std::printf("  M%d child=%d raw=%.12f norm=%.12f mask=%04x sign=%+d\n",
                    i, mutation.branchChild, mutation.rawTime,
                    mutation.normTime,
                    static_cast<unsigned>(mutation.descendantMask),
                    static_cast<int>(mutation.sign));
    }
}

static void testGoldenAndMusicalSnapshot() {
    KingmanGenerator goldenGenerator(42u, 54u);
    Tree golden;
    if (!goldenGenerator.generate(4, 0.2, golden)) {
        fail("golden n=4 tree generation failed");
        return;
    }
    printTree("golden n=4 seed=42/54 mutate=.2", golden);
    struct TraceEvent {
        double time;
        char kind;
        int index;
        TraceEvent(double eventTime, char eventKind, int eventIndex)
            : time(eventTime), kind(eventKind), index(eventIndex) {}
    };
    struct TraceLess {
        bool operator()(const TraceEvent& a, const TraceEvent& b) const {
            if (a.time != b.time)
                return a.time < b.time;
            if (a.kind != b.kind)
                return a.kind < b.kind;
            return a.index < b.index;
        }
    };
    std::vector<TraceEvent> trace;
    for (int i = golden.leafCount; i < golden.nodeCount; ++i)
        trace.push_back(TraceEvent(golden.nodes[i].normTime, 'N', i));
    for (int i = 0; i < golden.mutationCount; ++i)
        trace.push_back(TraceEvent(golden.mutations[i].normTime, 'M', i));
    std::sort(trace.begin(), trace.end(), TraceLess());
    std::printf("  ancestry trace: PRESENT");
    for (std::size_t i = 0; i < trace.size(); ++i)
        std::printf(" -> %c%d@%.6f", trace[i].kind, trace[i].index, trace[i].time);
    std::printf("\n  descent trace:  MRCA");
    for (std::size_t i = trace.size(); i-- > 0;)
        std::printf(" -> %c%d@%.6f", trace[i].kind, trace[i].index, trace[i].time);
    std::printf(" -> PRESENT\n");
    for (std::size_t i = 0; i < trace.size(); ++i) {
        const double lower = i == 0 ? 0.0
            : 0.5 * (trace[i - 1].time + trace[i].time);
        const double upper = i + 1 == trace.size() ? 1.0
            : 0.5 * (trace[i].time + trace[i + 1].time);
        const EventFacts upward = Playback::crossingFacts(
            golden, lower, trace[i].time);
        const EventFacts downward = Playback::crossingFacts(
            golden, trace[i].time, lower);
        if (upward.nodeCrossed != (trace[i].kind == 'N')
            || upward.mutationCrossed != (trace[i].kind == 'M')
            || downward.nodeCrossed != (trace[i].kind == 'N')
            || downward.mutationCrossed != (trace[i].kind == 'M'))
            fail("golden forward/reverse event trace mismatch at event %zu", i);
        if (upper > trace[i].time) {
            if (!factsEmpty(Playback::crossingFacts(
                    golden, trace[i].time, upper))
                || !factsEmpty(Playback::crossingFacts(
                    golden, upper, trace[i].time)))
                fail("golden exact event repeated before its next boundary at %zu", i);
        }
    }
    const int expectedOrder[] = {2, 1, 3, 0};
    const int expectedLeft[] = {1, 4, 2};
    const int expectedRight[] = {3, 0, 5};
    const double expectedNodeRaw[] = {
        0.076923861209, 0.298160407650, 0.587010166666
    };
    const double expectedNodeNorm[] = {
        0.131043490518, 0.507930568466, 1.0
    };
    const std::uint16_t expectedNodeMask[] = {0x0006u, 0x000eu, 0x000fu};
    const int expectedMutationChild[] = {0, 4, 5};
    const double expectedMutationRaw[] = {
        0.188914316254, 0.231403901114, 0.361497320077
    };
    const double expectedMutationNorm[] = {
        0.321824607106, 0.394207654747, 0.615828039453
    };
    const int expectedMutationSign[] = {1, -1, 1};
    const std::uint16_t expectedMutationMask[] = {0x0008u, 0x0006u, 0x000eu};
    bool goldenMatches = golden.mutationCount == 3;
    for (int i = 0; i < 4; ++i)
        goldenMatches = goldenMatches && golden.leafOrder[i] == expectedOrder[i];
    for (int i = 0; i < 3; ++i) {
        const Node& node = golden.nodes[4 + i];
        goldenMatches = goldenMatches && node.left == expectedLeft[i]
            && node.right == expectedRight[i]
            && near(node.rawTime, expectedNodeRaw[i])
            && near(node.normTime, expectedNodeNorm[i])
            && node.descendantMask == expectedNodeMask[i];
        const Mutation& mutation = golden.mutations[i];
        goldenMatches = goldenMatches
            && mutation.branchChild == expectedMutationChild[i]
            && near(mutation.rawTime, expectedMutationRaw[i])
            && near(mutation.normTime, expectedMutationNorm[i])
            && mutation.sign == expectedMutationSign[i]
            && mutation.descendantMask == expectedMutationMask[i];
    }
    const PlaybackValues goldenPresent = Playback::canonicalValues(golden, 0.0);
    const int expectedGoldenTraits[] = {0, 0, 0, 2};
    for (int i = 0; i < 4; ++i)
        goldenMatches = goldenMatches
            && goldenPresent.traitSteps[i] == expectedGoldenTraits[i];
    goldenMatches = goldenMatches
        && near(goldenPresent.meanPairwiseMutationDistance, 5.0 / 3.0, 1e-14);
    if (!goldenMatches)
        fail("fixed n=4 topology/mutation/playback golden changed");

    KingmanGenerator defaultGenerator(UINT64_C(0x4c494e4541474553),
                                      UINT64_C(0x434f414c45534345));
    Tree defaultTree;
    if (!defaultGenerator.generate(8, 0.4, defaultTree)) {
        fail("default musical snapshot generation failed");
        return;
    }
    const PlaybackValues present = Playback::canonicalValues(defaultTree, 0.0);
    const PlaybackValues root = Playback::canonicalValues(defaultTree, 1.0);
    std::printf("default n=8 mutate=.4 present steps:");
    for (int i = 0; i < present.sampleCount; ++i)
        std::printf(" %d", present.traitSteps[i]);
    std::printf("  diversity=%.5f (STEP=2 -> volts:",
                present.meanPairwiseMutationDistance);
    for (int i = 0; i < present.sampleCount; ++i)
        std::printf(" %+.3f", present.traitSteps[i] / 6.0);
    std::printf(")\n");
    if (present.sampleCount != 8 || present.lineageCount != 8
        || present.activeLineageCount != 8 || root.lineageCount != 1
        || root.activeLineageCount != 1 || root.sampleCount != 8
        || root.meanPairwiseMutationDistance != 0.0)
        fail("default snapshot endpoint state");
    for (int i = 0; i < 8; ++i) {
        if (root.traitSteps[i] != 0)
            fail("default snapshot did not converge to exact root unison");
    }

    const int ensemble = 512;
    double distinctSum = 0.0;
    double rangeSum = 0.0;
    double diversitySum = 0.0;
    int railRisk = 0;
    for (int draw = 0; draw < ensemble; ++draw) {
        Tree tree;
        if (!defaultGenerator.generate(8, 0.4, tree)) {
            fail("default musical ensemble generation failed");
            return;
        }
        const PlaybackValues values = Playback::canonicalValues(tree, 0.0);
        int minimum = values.traitSteps[0];
        int maximum = values.traitSteps[0];
        int unique[8];
        int uniqueCount = 0;
        for (int leaf = 0; leaf < 8; ++leaf) {
            minimum = std::min(minimum, values.traitSteps[leaf]);
            maximum = std::max(maximum, values.traitSteps[leaf]);
            bool seen = false;
            for (int i = 0; i < uniqueCount; ++i)
                seen = seen || unique[i] == values.traitSteps[leaf];
            if (!seen)
                unique[uniqueCount++] = values.traitSteps[leaf];
            if (std::abs(values.traitSteps[leaf]) >= 60)
                ++railRisk;
        }
        distinctSum += uniqueCount;
        rangeSum += maximum - minimum;
        diversitySum += values.meanPairwiseMutationDistance;
    }
    const double meanDistinct = distinctSum / ensemble;
    const double meanRange = rangeSum / ensemble;
    const double meanDiversity = diversitySum / ensemble;
    std::printf("default ensemble: distinct=%.3f range=%.3f steps (%.3f V at STEP=2) "
                "diversity=%.3f rail-risk=%d\n",
                meanDistinct, meanRange, meanRange / 6.0,
                meanDiversity, railRisk);
    if (meanDistinct < 2.0 || meanRange < 1.5 || meanDiversity < 0.5
        || railRisk != 0)
        fail("default mutation/STEP calibration did not form a useful bounded cluster");
}

static bool factsEmpty(const EventFacts& facts) {
    return !facts.nodeCrossed && !facts.mutationCrossed
        && !facts.arrivedAtMrca && !facts.departedMrca
        && !facts.startedAtMrca && !facts.wrapped;
}

static bool sameFacts(const EventFacts& a, const EventFacts& b) {
    return a.nodeCrossed == b.nodeCrossed
        && a.mutationCrossed == b.mutationCrossed
        && a.arrivedAtMrca == b.arrivedAtMrca
        && a.departedMrca == b.departedMrca
        && a.startedAtMrca == b.startedAtMrca
        && a.wrapped == b.wrapped;
}

static bool sameValues(const PlaybackValues& a, const PlaybackValues& b) {
    if (a.sampleCount != b.sampleCount || a.lineageCount != b.lineageCount
        || a.activeLineageCount != b.activeLineageCount
        || a.meanPairwiseMutationDistance != b.meanPairwiseMutationDistance)
        return false;
    for (int i = 0; i < MAX_LEAVES; ++i) {
        if (a.traitSteps[i] != b.traitSteps[i]
            || a.activeLineageMasks[i] != b.activeLineageMasks[i])
            return false;
    }
    return true;
}

static bool validActivePartition(const PlaybackValues& values) {
    if (!finitePlaybackValues(values)
        || values.lineageCount != values.activeLineageCount)
        return false;
    std::uint16_t combined = 0u;
    for (int i = 0; i < values.activeLineageCount; ++i) {
        const std::uint16_t mask = values.activeLineageMasks[i];
        if (mask == 0u || (combined & mask) != 0u)
            return false;
        combined = static_cast<std::uint16_t>(combined | mask);
        int firstTrait = 0;
        bool found = false;
        for (int leaf = 0; leaf < values.sampleCount; ++leaf) {
            if ((mask & static_cast<std::uint16_t>(1u << leaf)) == 0u)
                continue;
            if (!found) {
                firstTrait = values.traitSteps[leaf];
                found = true;
            }
            else if (values.traitSteps[leaf] != firstTrait) {
                return false;
            }
        }
        if (!found)
            return false;
    }
    return combined == allLeavesMask(values.sampleCount);
}

static bool traitsAndDiversityMatchBruteForce(const Tree& tree, double cursor,
                                               const PlaybackValues& values) {
    int traits[MAX_LEAVES] = {};
    double pairDistance = 0.0;
    for (int m = 0; m < tree.mutationCount; ++m) {
        const Mutation& mutation = tree.mutations[m];
        if (!(mutation.normTime > cursor))
            continue;
        for (int leaf = 0; leaf < tree.leafCount; ++leaf) {
            if ((mutation.descendantMask
                 & static_cast<std::uint16_t>(1u << leaf)) != 0u)
                traits[leaf] += mutation.sign;
        }
        for (int first = 0; first < tree.leafCount; ++first) {
            const bool firstHas = (mutation.descendantMask
                & static_cast<std::uint16_t>(1u << first)) != 0u;
            for (int second = first + 1; second < tree.leafCount; ++second) {
                const bool secondHas = (mutation.descendantMask
                    & static_cast<std::uint16_t>(1u << second)) != 0u;
                pairDistance += firstHas != secondHas ? 1.0 : 0.0;
            }
        }
    }
    for (int leaf = 0; leaf < tree.leafCount; ++leaf) {
        if (traits[leaf] != values.traitSteps[leaf])
            return false;
    }
    const double pairs = 0.5 * tree.leafCount * (tree.leafCount - 1);
    return near(values.meanPairwiseMutationDistance,
                pairs > 0.0 ? pairDistance / pairs : 0.0, 1e-13);
}

static Tree makeTwoLeafTree(bool cancellation) {
    Tree tree;
    tree.leafCount = 2;
    tree.nodeCount = 3;
    tree.root = 2;
    tree.nodes[2].left = 0;
    tree.nodes[2].right = 1;
    tree.nodes[2].rawTime = 1.0;
    tree.generatedMutate = 0.25;
    tree.mutationLambda = 4.0;
    tree.mutationCount = 2;
    tree.mutations[0].branchChild = 0;
    tree.mutations[0].rawTime = 0.25;
    tree.mutations[0].sign = 1;
    tree.mutations[1].branchChild = cancellation ? 0 : 1;
    tree.mutations[1].rawTime = 0.75;
    tree.mutations[1].sign = -1;
    if (!tree.validateAndDerive())
        fail("two-leaf playback fixture did not derive");
    return tree;
}

static void testExactPlaybackSemantics() {
    const Tree tree = makeTwoLeafTree(false);
    const PlaybackValues present = Playback::canonicalValues(tree, 0.0);
    const PlaybackValues firstBoundary = Playback::canonicalValues(tree, 0.25);
    const PlaybackValues secondBoundary = Playback::canonicalValues(tree, 0.75);
    const PlaybackValues root = Playback::canonicalValues(tree, 1.0);
    if (present.traitSteps[0] != 1 || present.traitSteps[1] != -1
        || present.meanPairwiseMutationDistance != 2.0
        || firstBoundary.traitSteps[0] != 0
        || firstBoundary.traitSteps[1] != -1
        || firstBoundary.meanPairwiseMutationDistance != 1.0
        || secondBoundary.traitSteps[0] != 0
        || secondBoundary.traitSteps[1] != 0
        || secondBoundary.meanPairwiseMutationDistance != 0.0
        || root.lineageCount != 1 || root.activeLineageCount != 1
        || root.traitSteps[0] != 0 || root.traitSteps[1] != 0)
        fail("canonical exact-event boundary ownership");

    EventFacts facts = Playback::crossingFacts(tree, 0.0, 0.25);
    if (!facts.mutationCrossed || facts.nodeCrossed)
        fail("upward arrival at mutation did not cross it");
    if (!factsEmpty(Playback::crossingFacts(tree, 0.25, 0.5)))
        fail("upward departure from exact mutation double-fired");
    if (!factsEmpty(Playback::crossingFacts(tree, 0.5, 0.25)))
        fail("downward arrival at exact mutation fired before state changed");
    facts = Playback::crossingFacts(tree, 0.25, 0.0);
    if (!facts.mutationCrossed || facts.nodeCrossed)
        fail("downward departure from exact mutation did not cross it");
    facts = Playback::crossingFacts(tree, 0.75, 1.0);
    if (!facts.nodeCrossed || !facts.arrivedAtMrca
        || facts.departedMrca || !facts.requestsMrcaPulse())
        fail("ancestry root-arrival facts");
    facts = Playback::crossingFacts(tree, 1.0, 0.75);
    if (!facts.nodeCrossed || !facts.departedMrca
        || facts.arrivedAtMrca || facts.requestsMrcaPulse())
        fail("descent root-departure facts");
    facts = Playback::crossingFacts(tree, 0.0, 1.0);
    if (!facts.nodeCrossed || !facts.mutationCrossed
        || !facts.arrivedAtMrca || !facts.requestsMrcaPulse())
        fail("large ancestry scrub did not aggregate each event type once");
    facts = Playback::crossingFacts(tree, 1.0, 0.0);
    if (!facts.nodeCrossed || !facts.mutationCrossed
        || !facts.departedMrca || facts.requestsMrcaPulse())
        fail("large descent scrub did not aggregate each event type once");
    if (!factsEmpty(Playback::crossingFacts(tree, 0.25, 0.25)))
        fail("stationary exact-event cursor fired");

    const Tree cancellation = makeTwoLeafTree(true);
    const PlaybackValues cancelled = Playback::canonicalValues(cancellation, 0.0);
    if (cancelled.traitSteps[0] != 0 || cancelled.traitSteps[1] != 0
        || cancelled.meanPairwiseMutationDistance != 2.0)
        fail("diversity incorrectly collapsed with signed trait cancellation");
}

static void testPlaybackTransportAndCache() {
    KingmanGenerator generator(1234u, 5678u);
    Tree tree;
    if (!generator.generate(16, 1.0, tree)) {
        fail("playback transport tree generation failed");
        return;
    }

    Playback playback;
    EventFacts facts = playback.installTreeAtSource(tree);
    if (!facts.startedAtMrca || facts.nodeCrossed
        || !facts.requestsMrcaPulse() || playback.cursor() != 1.0
        || playback.direction() != Direction::Descent || !playback.loop()
        || !playback.running()
        || !sameValues(playback.values(), Playback::canonicalValues(tree, 1.0)))
        fail("default descent install/root-start semantics");

    facts = playback.advance(tree, 0.1);
    if (!facts.departedMrca || !facts.nodeCrossed
        || facts.requestsMrcaPulse()
        || !sameValues(playback.values(),
                       Playback::canonicalValues(tree, playback.cursor())))
        fail("first descent movement/root departure semantics");

    for (int i = 0; i < 20000; ++i) {
        if (i % 13 == 0) {
            playback.setDirection(playback.direction() == Direction::Ancestry
                ? Direction::Descent : Direction::Ancestry);
        }
        playback.advance(tree, 0.000731);
        const PlaybackValues reference = Playback::canonicalValues(tree, playback.cursor());
        if (!sameValues(playback.values(), reference)
            || !validActivePartition(playback.values())
            || !traitsAndDiversityMatchBruteForce(
                tree, playback.cursor(), playback.values())) {
            fail("cached alternating-direction advance diverged at step %d cursor %.17g",
                 i, playback.cursor());
            return;
        }
    }

    for (int i = 0; i < 100000; ++i) {
        const int numerator = (i & 1) ? (i * 7919) % 10001
                                      : (i * 3571) % 10001;
        const double target = static_cast<double>(numerator) / 10000.0;
        playback.scrub(tree, target);
        const PlaybackValues reference = Playback::canonicalValues(tree, target);
        if (!sameValues(playback.values(), reference)
            || playback.values().lineageCount
                != playback.values().activeLineageCount) {
            fail("cached arbitrary scrub diverged at step %d target %.4f", i, target);
            return;
        }
        if (i % 997 == 0
            && (!validActivePartition(playback.values())
                || !traitsAndDiversityMatchBruteForce(
                    tree, target, playback.values()))) {
            fail("scrub canonical state failed brute-force check at step %d", i);
            return;
        }
    }
    std::printf("playback cache: 20k reversals + 100k arbitrary scrubs exact\n");

    Tree replacementTree;
    if (!generator.generate(3, 0.2, replacementTree)) {
        fail("playback replacement tree generation failed");
        return;
    }
    playback.synchronize(replacementTree, 0.37);
    const double replacementTargets[] = {0.37, 1.0, 0.0};
    for (std::size_t i = 0;
         i < sizeof(replacementTargets) / sizeof(replacementTargets[0]); ++i) {
        playback.scrub(replacementTree, replacementTargets[i]);
        if (!sameValues(playback.values(), Playback::canonicalValues(
                replacementTree, replacementTargets[i]))) {
            fail("high-to-low mutation cache replacement diverged at target %.2f",
                 replacementTargets[i]);
            return;
        }
    }

    Playback restoredPlayback;
    restoredPlayback.restore(
        tree, 0.4321, Direction::Ancestry, false, true);
    if (restoredPlayback.cursor() != 0.4321
        || restoredPlayback.direction() != Direction::Ancestry
        || restoredPlayback.loop() || !restoredPlayback.running()
        || !sameValues(restoredPlayback.values(),
                       Playback::canonicalValues(tree, 0.4321)))
        fail("valid playback restore did not reproduce saved transport/state");
    const EventFacts expectedRestoredAdvance = Playback::crossingFacts(
        tree, 0.4321, 0.5021);
    EventFacts restoredFacts = restoredPlayback.advance(tree, 0.07);
    if (!sameFacts(restoredFacts, expectedRestoredAdvance)
        || !near(restoredPlayback.cursor(), 0.5021, 1e-15)
        || !sameValues(restoredPlayback.values(),
                       Playback::canonicalValues(tree, restoredPlayback.cursor())))
        fail("restored playback diverged on internal advance");
    const EventFacts expectedRestoredScrub = Playback::crossingFacts(
        tree, restoredPlayback.cursor(), 0.81);
    restoredFacts = restoredPlayback.scrub(tree, 0.81);
    if (!sameFacts(restoredFacts, expectedRestoredScrub)
        || restoredPlayback.cursor() != 0.81
        || !sameValues(restoredPlayback.values(),
                       Playback::canonicalValues(tree, 0.81)))
        fail("restored playback diverged on absolute scrub");

    playback.setDirection(Direction::Ancestry);
    playback.setLoop(true);
    playback.synchronize(tree, 0.91);
    const EventFacts expectedAncestryTail = Playback::crossingFacts(tree, 0.91, 1.0);
    facts = playback.advance(tree, 99.0);
    if (!facts.wrapped || !facts.arrivedAtMrca || !facts.nodeCrossed
        || facts.mutationCrossed != expectedAncestryTail.mutationCrossed
        || !facts.requestsMrcaPulse() || playback.cursor() != 0.0
        || !sameValues(playback.values(), Playback::canonicalValues(tree, 0.0)))
        fail("ancestry wrap/overshoot discard/source state");

    playback.setDirection(Direction::Descent);
    playback.synchronize(tree, 0.09);
    const EventFacts expectedDescentTail = Playback::crossingFacts(tree, 0.09, 0.0);
    if (!expectedDescentTail.nodeCrossed || !expectedDescentTail.mutationCrossed)
        fail("descent wrap fixture lacks intermediate NODE/MUTATION events");
    facts = playback.advance(tree, 99.0);
    if (!facts.wrapped || !facts.startedAtMrca
        || facts.nodeCrossed != expectedDescentTail.nodeCrossed
        || facts.mutationCrossed != expectedDescentTail.mutationCrossed
        || !facts.requestsMrcaPulse() || playback.cursor() != 1.0
        || !sameValues(playback.values(), Playback::canonicalValues(tree, 1.0)))
        fail("descent wrap/overshoot discard/root state");

    playback.setDirection(Direction::Ancestry);
    playback.setLoop(false);
    playback.resetToSource(tree);
    playback.advance(tree, 99.0);
    if (playback.cursor() != 1.0 || playback.running())
        fail("non-loop ancestry did not stop at root");
    playback.setLoop(true);
    const double stoppedCursor = playback.cursor();
    facts = playback.advance(tree, 0.2);
    if (playback.running() || playback.cursor() != stoppedCursor || !factsEmpty(facts))
        fail("LOOP toggle incorrectly restarted stopped transport");
    playback.resetToSource(tree);
    if (!playback.running() || playback.cursor() != 0.0)
        fail("RESET did not restart newly enabled loop at source");

    playback.setLoop(false);
    playback.advance(tree, 99.0);
    playback.setDirection(Direction::Descent);
    if (!playback.running())
        fail("direction reversal did not resume stopped endpoint");
    facts = playback.advance(tree, 1e-6);
    if (!facts.nodeCrossed || !facts.departedMrca
        || facts.requestsMrcaPulse())
        fail("reversed descent did not cross root once on departure");

    playback.synchronize(tree, 0.4321);
    const double beforeNonFinite = playback.cursor();
    const PlaybackValues valuesBeforeNonFinite = playback.values();
    facts = playback.scrub(tree, std::numeric_limits<double>::quiet_NaN());
    if (!factsEmpty(facts) || playback.cursor() != beforeNonFinite
        || !sameValues(playback.values(), valuesBeforeNonFinite))
        fail("non-finite TIME did not hold previous cursor/state");
    const double invalidDistances[] = {
        0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()
    };
    for (std::size_t i = 0;
         i < sizeof(invalidDistances) / sizeof(invalidDistances[0]); ++i) {
        facts = playback.advance(tree, invalidDistances[i]);
        if (!factsEmpty(facts) || playback.cursor() != beforeNonFinite
            || !sameValues(playback.values(), valuesBeforeNonFinite))
            fail("hostile advance distance %zu changed transport", i);
    }
    playback.restore(tree, NAN, Direction::Descent, false, true);
    if (playback.cursor() != 1.0
        || !sameValues(playback.values(), Playback::canonicalValues(tree, 1.0)))
        fail("non-finite restored cursor did not select direction source");

    playback.setDirection(Direction::Ancestry);
    playback.setLoop(false);
    playback.synchronize(tree, 0.9);
    const double minimumDelta = 0.1 * std::pow(2.0, -12.0) / 48000.0;
    const double lowRateStart = playback.cursor();
    for (int i = 0; i < 1000; ++i)
        playback.advance(tree, minimumDelta);
    if (!(playback.cursor() > lowRateStart)
        || !near(playback.cursor(), lowRateStart + 1000.0 * minimumDelta, 1e-12))
        fail("minimum RATE cursor stalled or lost double precision "
             "(%.17g vs %.17g)", playback.cursor(),
             lowRateStart + 1000.0 * minimumDelta);

    playback.setLoop(true);
    playback.synchronize(tree, 0.5);
    const std::chrono::steady_clock::time_point begin
        = std::chrono::steady_clock::now();
    for (int i = 0; i < 1000000; ++i)
        playback.scrub(tree, (i & 1) ? 0.500001 : 0.5);
    const std::chrono::steady_clock::time_point end
        = std::chrono::steady_clock::now();
    const double micros = std::chrono::duration<double, std::micro>(end - begin).count();
    std::printf("playback benchmark: 1M tiny reversible scrubs %.1f us\n", micros);
    if (!sameValues(playback.values(),
                    Playback::canonicalValues(tree, playback.cursor())))
        fail("benchmark scrub left stale cached state");

    Tree cappedTree = tree;
    cappedTree.mutationCount = MAX_MUTATIONS;
    cappedTree.mutationsTruncated = true;
    cappedTree.generatedMutate = 1.0;
    for (int i = 0; i < MAX_MUTATIONS; ++i) {
        const int childIndex = i % cappedTree.root;
        const Node& child = cappedTree.nodes[childIndex];
        const Node& parent = cappedTree.nodes[child.parent];
        const double fraction = static_cast<double>(i + 1)
                              / static_cast<double>(MAX_MUTATIONS + 1);
        cappedTree.mutations[i].branchChild = childIndex;
        cappedTree.mutations[i].rawTime = child.rawTime
            + fraction * (parent.rawTime - child.rawTime);
        cappedTree.mutations[i].sign = (i & 1) ? -1 : 1;
    }
    if (!cappedTree.validateAndDerive()
        || !structuralInvariants(cappedTree, cappedTree.leafCount)) {
        fail("synthetic mutation-cap playback tree did not validate");
        return;
    }

    const int fullScrubCount = 200000;
    volatile unsigned fullScrubFacts = 0u;
    double fullMicros = std::numeric_limits<double>::infinity();
    for (int trial = 0; trial < 3; ++trial) {
        playback.synchronize(cappedTree, 0.0);
        const std::chrono::steady_clock::time_point fullBegin
            = std::chrono::steady_clock::now();
        for (int i = 0; i < fullScrubCount; ++i) {
            const EventFacts fullFacts = playback.scrub(
                cappedTree, (i & 1) ? 0.0 : 1.0);
            fullScrubFacts += static_cast<unsigned>(fullFacts.nodeCrossed)
                            + static_cast<unsigned>(fullFacts.mutationCrossed);
        }
        const std::chrono::steady_clock::time_point fullEnd
            = std::chrono::steady_clock::now();
        fullMicros = std::min(fullMicros,
            std::chrono::duration<double, std::micro>(
                fullEnd - fullBegin).count());
    }
    const double scaledFullMicros = fullMicros * (1000000.0 / fullScrubCount);
    const double fullToTinyRatio = scaledFullMicros / std::max(1.0, micros);
    std::printf("                    1M-equivalent capped full-tree scrubs %.1f us "
                "(%.2fx tiny)\n", scaledFullMicros, fullToTinyRatio);
    if (fullScrubFacts == 0u
        || !sameValues(playback.values(),
                       Playback::canonicalValues(cappedTree, 0.0)))
        fail("full-range benchmark missed events or left stale cached state");
    if (fullToTinyRatio > 64.0)
        fail("full-range TIME scrub regressed to per-mutation/leaf work (%.2fx)",
             fullToTinyRatio);

    Tree hostile;
    hostile.leafCount = std::numeric_limits<int>::max();
    hostile.nodeCount = std::numeric_limits<int>::max();
    hostile.mutationCount = std::numeric_limits<int>::max();
    const PlaybackValues hostileValues = Playback::canonicalValues(hostile, NAN);
    if (!finitePlaybackValues(hostileValues) || hostileValues.sampleCount != 0)
        fail("hostile malformed tree escaped bounded playback guards");
    if (!factsEmpty(Playback::crossingFacts(
            hostile, -std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity())))
        fail("non-finite crossing endpoints emitted facts");
}

static bool samePlaybackState(const Playback& a, const Playback& b) {
    return sameDoubleBits(a.cursor(), b.cursor())
        && a.direction() == b.direction()
        && a.loop() == b.loop()
        && a.running() == b.running()
        && sameValues(a.values(), b.values());
}

static bool sameSnapshot(const Snapshot& a, const Snapshot& b) {
    return a.version == b.version
        && sameTree(a.generated.tree, b.generated.tree)
        && sameRng(a.generated.nextRng, b.generated.nextRng)
        && sameDoubleBits(a.cursor, b.cursor)
        && a.direction == b.direction
        && a.loop == b.loop
        && a.running == b.running
        && sameFloatBits(a.nodePulseRemaining, b.nodePulseRemaining)
        && sameFloatBits(a.mutationPulseRemaining, b.mutationPulseRemaining)
        && sameFloatBits(a.mrcaPulseRemaining, b.mrcaPulseRemaining);
}

static void expectInvalidSnapshot(const Snapshot& input, const char* label) {
    Snapshot candidate = input;
    const Snapshot before = candidate;
    if (candidate.validateAndDerive())
        fail("invalid snapshot accepted: %s", label);
    else if (!sameSnapshot(candidate, before))
        fail("failed snapshot validation was not transactional: %s", label);
}

static bool buildSnapshotFixture(Snapshot& snapshot,
                                 KingmanGenerator& generator,
                                 Tree& tree, Playback& playback) {
    generator.seed(UINT64_C(0x5354415445464958),
                   UINT64_C(0x4c494e4541474553));
    if (!generator.generate(9, 0.53, tree))
        return false;
    playback = Playback();
    playback.setDirection(Direction::Ancestry);
    playback.setLoop(false);
    playback.installTreeAtSource(tree);
    playback.scrub(tree, 0.43125);

    snapshot = Snapshot();
    snapshot.generated = generator.capture(tree);
    snapshot.cursor = playback.cursor();
    snapshot.direction = playback.direction();
    snapshot.loop = playback.loop();
    snapshot.running = playback.running();
    snapshot.nodePulseRemaining = 0.000125f;
    snapshot.mutationPulseRemaining = 0.0005f;
    snapshot.mrcaPulseRemaining = Snapshot::maximumPulseRemaining();
    return snapshot.validateAndDerive();
}

static void testHex64Codec() {
    struct HexCase {
        std::uint64_t value;
        const char* canonical;
        const char* uppercase;
    };
    const HexCase cases[] = {
        {UINT64_C(0), "0000000000000000", "0000000000000000"},
        {UINT64_C(1), "0000000000000001", "0000000000000001"},
        {UINT64_C(0x0123456789abcdef), "0123456789abcdef", "0123456789ABCDEF"},
        {UINT64_C(0x8000000000000000), "8000000000000000", "8000000000000000"},
        {UINT64_MAX, "ffffffffffffffff", "FFFFFFFFFFFFFFFF"}
    };
    for (std::size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        char formatted[Hex64Codec::TEXT_SIZE];
        Hex64Codec::format(cases[i].value, formatted);
        if (std::strcmp(formatted, cases[i].canonical) != 0
            || formatted[Hex64Codec::DIGITS] != '\0')
            fail("Hex64 canonical formatting at case %zu", i);
        std::uint64_t parsed = UINT64_C(0xfeedfacecafebeef);
        if (!Hex64Codec::parse(formatted, Hex64Codec::DIGITS, parsed)
            || parsed != cases[i].value)
            fail("Hex64 canonical round trip at case %zu", i);
        parsed = 0u;
        if (!Hex64Codec::parseCString(cases[i].uppercase, parsed)
            || parsed != cases[i].value)
            fail("Hex64 uppercase parse at case %zu", i);
    }

    const std::uint64_t sentinel = UINT64_C(0x1122334455667788);
    const char embeddedNull[Hex64Codec::DIGITS] = {
        '0', '1', '2', '3', '\0', '5', '6', '7',
        '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'
    };
    const char tooLong[] = "0123456789abcdef0";
    const char invalidDigit[] = "0123456789abcdeg";
    std::uint64_t parsed = sentinel;
    if (Hex64Codec::parse(0, Hex64Codec::DIGITS, parsed)
        || parsed != sentinel)
        fail("Hex64 null parse was not transactional");
    if (Hex64Codec::parse("0123456789abcde", 15u, parsed)
        || parsed != sentinel)
        fail("Hex64 short parse was not transactional");
    if (Hex64Codec::parse(tooLong, 17u, parsed) || parsed != sentinel)
        fail("Hex64 long parse was not transactional");
    if (Hex64Codec::parse(invalidDigit, 16u, parsed) || parsed != sentinel)
        fail("Hex64 invalid-digit parse was not transactional");
    if (Hex64Codec::parse(embeddedNull, 16u, parsed) || parsed != sentinel)
        fail("Hex64 embedded-null parse was not transactional");
    if (Hex64Codec::parseCString(tooLong, parsed) || parsed != sentinel)
        fail("Hex64 overlong C string was not transactional");
    if (Hex64Codec::parseCString(0, parsed) || parsed != sentinel)
        fail("Hex64 null C string was not transactional");

    std::printf("state Hex64: canonical lowercase + exact case-insensitive parse\n");
}

static void testSnapshotState() {
    Snapshot snapshot;
    KingmanGenerator sourceGenerator;
    Tree sourceTree;
    Playback sourcePlayback;
    if (!buildSnapshotFixture(snapshot, sourceGenerator,
                              sourceTree, sourcePlayback)) {
        fail("snapshot fixture generation failed");
        return;
    }

    Snapshot poisoned = snapshot;
    Tree& poisonedTree = poisoned.generated.tree;
    for (int i = 0; i < poisonedTree.nodeCount; ++i) {
        poisonedTree.nodes[i].parent = (i + 3) % poisonedTree.nodeCount;
        poisonedTree.nodes[i].normTime = NAN;
        poisonedTree.nodes[i].descendantMask = 0u;
        poisonedTree.nodes[i].displayX = -42.0;
    }
    poisonedTree.leafOrder.fill(-1);
    poisonedTree.totalRawBranchLength = NAN;
    poisonedTree.mutationLambda = NAN;
    for (int i = 0; i < poisonedTree.mutationCount; ++i) {
        poisonedTree.mutations[i].normTime = NAN;
        poisonedTree.mutations[i].descendantMask = 0u;
    }
    if (!poisoned.validateAndDerive()
        || !sameTree(poisoned.generated.tree, snapshot.generated.tree))
        fail("Snapshot did not rebuild poisoned derived tree fields");

    Snapshot invalid = snapshot;
    invalid.version = 0;
    expectInvalidSnapshot(invalid, "version zero");
    invalid = snapshot;
    invalid.version = coalescent::lineages::STATE_VERSION + 1;
    expectInvalidSnapshot(invalid, "future version");
    invalid = snapshot;
    invalid.generated.nextRng.increment &= ~UINT64_C(1);
    expectInvalidSnapshot(invalid, "even PCG increment");
    invalid = snapshot;
    invalid.generated.tree.nodes[invalid.generated.tree.root].left = MAX_NODES;
    expectInvalidSnapshot(invalid, "invalid topology index");
    invalid = snapshot;
    invalid.generated.tree.nodes[invalid.generated.tree.root].rawTime = NAN;
    expectInvalidSnapshot(invalid, "non-finite root time");

    const double invalidCursors[] = {
        -std::numeric_limits<double>::epsilon(),
        1.0 + std::numeric_limits<double>::epsilon(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()
    };
    for (std::size_t i = 0;
         i < sizeof(invalidCursors) / sizeof(invalidCursors[0]); ++i) {
        invalid = snapshot;
        invalid.cursor = invalidCursors[i];
        expectInvalidSnapshot(invalid, "cursor");
    }
    invalid = snapshot;
    invalid.direction = static_cast<Direction>(-1);
    expectInvalidSnapshot(invalid, "negative direction");
    invalid = snapshot;
    invalid.direction = static_cast<Direction>(2);
    expectInvalidSnapshot(invalid, "unknown direction");

    // The running bit is authoritative on restore, but only states reachable
    // through Playback may be persisted. Stopped transport can only sit at the
    // selected direction's destination, and a non-looping transport cannot
    // remain running once it has reached that destination.
    invalid = snapshot;
    invalid.cursor = 0.5;
    invalid.direction = Direction::Ancestry;
    invalid.loop = false;
    invalid.running = false;
    expectInvalidSnapshot(invalid, "stopped ancestry away from destination");
    invalid = snapshot;
    invalid.cursor = 0.5;
    invalid.direction = Direction::Descent;
    invalid.loop = true;
    invalid.running = false;
    expectInvalidSnapshot(invalid, "stopped descent away from destination");
    invalid = snapshot;
    invalid.cursor = 1.0;
    invalid.direction = Direction::Ancestry;
    invalid.loop = false;
    invalid.running = true;
    expectInvalidSnapshot(invalid, "running non-loop ancestry destination");
    invalid = snapshot;
    invalid.cursor = 0.0;
    invalid.direction = Direction::Descent;
    invalid.loop = false;
    invalid.running = true;
    expectInvalidSnapshot(invalid, "running non-loop descent destination");

    // These endpoint states are reachable: transport can stop with LOOP off,
    // then LOOP can be enabled without implicitly restarting it.
    Snapshot stoppedEndpoint = snapshot;
    stoppedEndpoint.cursor = 1.0;
    stoppedEndpoint.direction = Direction::Ancestry;
    stoppedEndpoint.loop = true;
    stoppedEndpoint.running = false;
    if (!stoppedEndpoint.validateAndDerive())
        fail("reachable stopped ancestry endpoint was rejected");
    stoppedEndpoint = snapshot;
    stoppedEndpoint.cursor = 0.0;
    stoppedEndpoint.direction = Direction::Descent;
    stoppedEndpoint.loop = true;
    stoppedEndpoint.running = false;
    if (!stoppedEndpoint.validateAndDerive())
        fail("reachable stopped descent endpoint was rejected");

    const float invalidPulses[] = {
        -std::numeric_limits<float>::epsilon(),
        std::nextafter(Snapshot::maximumPulseRemaining(),
                       std::numeric_limits<float>::infinity()),
        std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::infinity()
    };
    for (int pulse = 0; pulse < 3; ++pulse) {
        for (std::size_t i = 0;
             i < sizeof(invalidPulses) / sizeof(invalidPulses[0]); ++i) {
            invalid = snapshot;
            if (pulse == 0)
                invalid.nodePulseRemaining = invalidPulses[i];
            else if (pulse == 1)
                invalid.mutationPulseRemaining = invalidPulses[i];
            else
                invalid.mrcaPulseRemaining = invalidPulses[i];
            expectInvalidSnapshot(invalid, "pulse remainder");
        }
    }

    KingmanGenerator installedGenerator(UINT64_C(0x99887766),
                                         UINT64_C(0x55443322));
    Tree installedTree;
    installedGenerator.generate(4, 0.1, installedTree);
    Playback installedPlayback;
    installedPlayback.setDirection(Direction::Descent);
    installedPlayback.setLoop(true);
    installedPlayback.installTreeAtSource(installedTree);
    if (!snapshot.install(installedGenerator, installedTree, installedPlayback)
        || !sameTree(installedTree, snapshot.generated.tree)
        || !sameRng(installedGenerator.rngState(), snapshot.generated.nextRng)
        || !samePlaybackState(installedPlayback, sourcePlayback))
        fail("valid Snapshot install did not reproduce composition/playback");

    Tree sourceNext;
    Tree installedNext;
    if (!sourceGenerator.generate(11, 0.73, sourceNext)
        || !installedGenerator.generate(11, 0.73, installedNext)
        || !sameTree(sourceNext, installedNext)
        || !sameRng(sourceGenerator.rngState(), installedGenerator.rngState()))
        fail("Snapshot install did not reproduce exact next NEW");

    Playback sourceContinuation = sourcePlayback;
    Playback installedContinuation = installedPlayback;
    for (int i = 0; i < 256; ++i) {
        if (i % 29 == 0) {
            const Direction direction = sourceContinuation.direction()
                == Direction::Ancestry ? Direction::Descent
                                       : Direction::Ancestry;
            sourceContinuation.setDirection(direction);
            installedContinuation.setDirection(direction);
        }
        const double target = static_cast<double>((i * 1879) % 1001) / 1000.0;
        const EventFacts sourceFacts = sourceContinuation.scrub(sourceTree, target);
        const EventFacts installedFacts = installedContinuation.scrub(installedTree, target);
        if (!sameFacts(sourceFacts, installedFacts)
            || !samePlaybackState(sourceContinuation, installedContinuation)) {
            fail("Snapshot playback continuation diverged at step %d", i);
            break;
        }
    }

    KingmanGenerator atomicGenerator(UINT64_C(0x123456), UINT64_C(0xabcdef));
    Tree atomicTree;
    atomicGenerator.generate(5, 0.2, atomicTree);
    Playback atomicPlayback;
    atomicPlayback.setDirection(Direction::Descent);
    atomicPlayback.setLoop(false);
    atomicPlayback.installTreeAtSource(atomicTree);
    atomicPlayback.scrub(atomicTree, 0.6);
    const Pcg32::State atomicRngBefore = atomicGenerator.rngState();
    const Tree atomicTreeBefore = atomicTree;
    const Playback atomicPlaybackBefore = atomicPlayback;
    invalid = snapshot;
    invalid.cursor = NAN;
    if (invalid.install(atomicGenerator, atomicTree, atomicPlayback)
        || !sameRng(atomicRngBefore, atomicGenerator.rngState())
        || !sameTree(atomicTreeBefore, atomicTree)
        || !samePlaybackState(atomicPlaybackBefore, atomicPlayback))
        fail("failed Snapshot install was not atomic");

    Snapshot factory;
    Snapshot factoryAgain;
    if (!makeFactorySnapshot(factory) || !makeFactorySnapshot(factoryAgain)
        || !factory.validateAndDerive() || !sameSnapshot(factory, factoryAgain)
        || factory.generated.tree.leafCount
            != coalescent::lineages::FACTORY_LEAF_COUNT
        || factory.generated.tree.generatedMutate
            != coalescent::lineages::FACTORY_MUTATE
        || factory.cursor != 1.0 || factory.direction != Direction::Descent
        || !factory.loop || !factory.running
        || factory.nodePulseRemaining != 0.f
        || factory.mutationPulseRemaining != 0.f
        || factory.mrcaPulseRemaining != 0.f)
        fail("factory fallback Snapshot invalid or nondeterministic");

    KingmanGenerator factoryGenerator(1u, 2u);
    Tree factoryTree;
    Playback factoryPlayback;
    if (!factory.install(factoryGenerator, factoryTree, factoryPlayback)
        || !structuralInvariants(factoryTree,
            coalescent::lineages::FACTORY_LEAF_COUNT)
        || factoryPlayback.cursor() != 1.0
        || factoryPlayback.values().lineageCount != 1
        || !sameValues(factoryPlayback.values(),
                       Playback::canonicalValues(factoryTree, 1.0)))
        fail("factory fallback Snapshot did not install a valid root state");

    std::printf("state Snapshot: validation, atomic install, factory, next NEW exact\n");
}

int main() {
    testRandom();
    testStructuralAndPersistence();
    testKingmanStatistics();
    testMutations();
    testGoldenAndMusicalSnapshot();
    testExactPlaybackSemantics();
    testPlaybackTransportAndCache();
    testHex64Codec();
    testSnapshotState();
    if (failures) {
        std::printf("Lineages stability: FAIL (%d findings)\n", failures);
        return 1;
    }
    std::printf("Lineages stability: PASS\n");
    return 0;
}
