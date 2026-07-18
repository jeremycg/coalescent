// SDK-free stability and behaviour contract for src/dsp/finches_field.hpp.
//
//   g++ -std=c++11 -O2 tools/stability/finches.cpp -o /tmp/finches && /tmp/finches
#include "../../src/dsp/finches_field.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>

using coalescent::FinchesField;

static bool valid(const FinchesField& field, float tolerance = 3e-6f) {
    double sum = 0.0;
    const std::array<float, FinchesField::kBins>& m = field.masses();
    for (int i = 0; i < FinchesField::kBins; ++i) {
        if (!std::isfinite(m[i]) || m[i] < 0.f || m[i] > 1.f
            || (m[i] > 0.f && m[i] < FinchesField::numericalExtinctionFloor()))
            return false;
        sum += m[i];
    }
    const FinchesField::Metrics& x = field.metrics();
    return std::fabs(sum - 1.0) <= tolerance && std::isfinite(x.mean) &&
           std::isfinite(x.spread) &&
           std::isfinite(x.lowTrait) && std::isfinite(x.highTrait) &&
           std::isfinite(x.divergence) && x.mean >= -1.f && x.mean <= 1.f &&
           x.spread >= 0.f && x.spread <= 1.f &&
           x.lowMass >= 0.f && x.lowMass <= 1.f &&
           x.highMass >= 0.f && x.highMass <= 1.f;
}

static bool sameState(const FinchesField::State& a,
                      const FinchesField::State& b) {
    return a.version == b.version
        && std::memcmp(a.mass.data(), b.mass.data(),
                       sizeof(float) * FinchesField::kBins) == 0
        && a.split == b.split
        && a.splitTimer == b.splitTimer
        && a.mergeTimer == b.mergeTimer;
}

static bool evolve(FinchesField& field, const FinchesField::Parameters& p,
                   float tau, float outerStep, int* splitEvents = 0, int* mergeEvents = 0) {
    int n = static_cast<int>(std::ceil(tau / outerStep));
    float h = tau / static_cast<float>(n);
    for (int i = 0; i < n; ++i) {
        field.advance(h, p);
        if (splitEvents && field.metrics().splitEvent) ++*splitEvents;
        if (mergeEvents && field.metrics().mergeEvent) ++*mergeEvents;
        if (!valid(field))
            return false;
    }
    return true;
}

static float maxMassDifference(const FinchesField& a, const FinchesField& b) {
    float error = 0.f;
    for (int i = 0; i < FinchesField::kBins; ++i)
        error = std::max(error, std::fabs(a.masses()[i] - b.masses()[i]));
    return error;
}

static float mirrorError(const FinchesField& field) {
    float error = 0.f;
    for (int i = 0; i < FinchesField::kBins / 2; ++i)
        error = std::max(error, std::fabs(field.masses()[i] -
                                         field.masses()[FinchesField::kBins - 1 - i]));
    return error;
}

static int significantPeakCount(const FinchesField& field) {
    const std::array<float, FinchesField::kBins>& m = field.masses();
    float maximum = *std::max_element(m.begin(), m.end());
    int count = 0;
    for (int i = 0; i < FinchesField::kBins; ++i) {
        float left = i ? m[i - 1] : -1.f;
        float right = i + 1 < FinchesField::kBins ? m[i + 1] : -1.f;
        if (m[i] >= 0.02f * maximum && m[i] >= left && m[i] >= right &&
            (m[i] > left || m[i] > right)) {
            // Treat the equal central pair produced by an even cell-centred grid
            // as one plateau rather than two biological peaks.
            if (i && m[i] == m[i - 1])
                continue;
            count++;
        }
    }
    return count;
}

int main() {
    int failures = 0;
    FinchesField::Parameters p;

    for (float mutation : {0.00002f, 0.00003f, 0.00004f, 0.00006f, 0.0001f}) {
        std::printf("calibration mutation=%g:", mutation);
        for (float branching : {2.1f, 2.8f}) {
            FinchesField f;
            p = FinchesField::Parameters(); p.mutation = mutation; p.branching = branching;
            evolve(f, p, 90.f, 0.16f);
            std::printf(" B%.1f=%c div%.3f spr%.3f", branching,
                        f.metrics().split ? 'S' : '.', f.metrics().divergence, f.metrics().spread);
        }
        std::printf("\n");
    }

    // A broad unresolved population remains one centered voice. In particular,
    // an equal-height shallow pair must not make both pitch identities choose
    // the lower candidate merely because the peak scan runs left to right.
    FinchesField unresolved;
    p = FinchesField::Parameters(); p.mutation = 0.00006f; p.branching = 2.1f;
    evolve(unresolved, p, 90.f, 0.16f);
    if (unresolved.metrics().split ||
        std::fabs(unresolved.metrics().lowTrait - unresolved.metrics().mean) > 1e-7f ||
        std::fabs(unresolved.metrics().highTrait - unresolved.metrics().mean) > 1e-7f) {
        std::printf("FAIL: unresolved pair did not retain one centered voice\n"); failures++;
    }

    // Calibrate the defining gesture: a stable ancestor, one accepted split,
    // then one accepted merge after crossing back below B=1.
    FinchesField gesture;
    int splitEvents = 0, mergeEvents = 0;
    p.mutation = 0.00002f; p.niche = 0.32f; p.environment = 0.f; p.branching = 0.65f;
    if (!evolve(gesture, p, 45.f, 0.08f, &splitEvents, &mergeEvents) || gesture.metrics().split) {
        std::printf("FAIL: stabilizing regime did not remain one peak\n"); failures++;
    }
    p.branching = 2.1f;
    if (!evolve(gesture, p, 90.f, 0.08f, &splitEvents, &mergeEvents) || !gesture.metrics().split) {
        std::printf("FAIL: disruptive regime did not branch (mean=%.4f spread=%.4f div=%.4f)\n",
                    gesture.metrics().mean, gesture.metrics().spread, gesture.metrics().divergence);
        for (int i = 0; i < FinchesField::kBins; ++i)
            std::printf("%s%.6f", i ? " " : "  density: ", gesture.masses()[i]);
        std::printf("\n");
        failures++;
    }
    float splitDivergence = gesture.metrics().divergence;
    float splitLowMass = gesture.metrics().lowMass;
    float splitHighMass = gesture.metrics().highMass;
    p.branching = 0.65f;
    if (!evolve(gesture, p, 90.f, 0.08f, &splitEvents, &mergeEvents) || gesture.metrics().split) {
        std::printf("FAIL: stabilizing return did not merge\n"); failures++;
    }
    std::printf("gesture: divergence=%.4f masses=%.3f/%.3f events split=%d merge=%d\n",
                splitDivergence, splitLowMass, splitHighMass, splitEvents, mergeEvents);
    if (splitEvents != 1 || mergeEvents != 1 || splitDivergence < 0.15f) {
        std::printf("FAIL: one-to-two-to-one event contract\n"); failures++;
    }

    // A centred ancestor under a centred environment must not pick a species
    // majority merely from convolution order or tridiagonal elimination order.
    FinchesField symmetric;
    p = FinchesField::Parameters(); p.mutation = 0.00002f; p.branching = 2.8f;
    evolve(symmetric, p, 100.f, 0.137f);
    float symmetryError = mirrorError(symmetric);
    std::printf("symmetry: max mirrored-bin error %.9f\n", symmetryError);
    if (!symmetric.metrics().split || symmetryError > 1e-7f ||
        std::fabs(symmetric.metrics().lowMass - symmetric.metrics().highMass) > 1e-5f) {
        std::printf("FAIL: centred branching acquired numerical directional bias\n"); failures++;
    }

    // Full exposed corners plus intermediate controls: finite, nonnegative,
    // normalized, bounded metrics. Seed each run to exercise asymmetric fields.
    int sweepRuns = 0;
    const float mutations[] = {FinchesField::mutationMin(), 0.00035f, FinchesField::mutationMax()};
    const float branchings[] = {FinchesField::branchingMin(), 0.9f, 1.1f, FinchesField::branchingMax()};
    const float niches[] = {FinchesField::nicheMin(), 0.32f, FinchesField::nicheMax()};
    const float environments[] = {FinchesField::environmentMin(), 0.f, FinchesField::environmentMax()};
    for (std::size_t mi = 0; mi < sizeof(mutations) / sizeof(*mutations); ++mi)
        for (std::size_t bi = 0; bi < sizeof(branchings) / sizeof(*branchings); ++bi)
            for (std::size_t ni = 0; ni < sizeof(niches) / sizeof(*niches); ++ni)
                for (std::size_t ei = 0; ei < sizeof(environments) / sizeof(*environments); ++ei) {
                    FinchesField f;
                    p.mutation = mutations[mi]; p.branching = branchings[bi];
                    p.niche = niches[ni]; p.environment = environments[ei];
                    f.reset(p.environment);
                    f.seed(-p.environment, 0.08f, 0.025f);
                    if (!evolve(f, p, 12.f, 0.16f)) {
                        std::printf("FAIL: sweep m=%g B=%g niche=%g env=%g\n",
                                    p.mutation, p.branching, p.niche, p.environment);
                        failures++;
                    }
                    sweepRuns++;
                }
    std::printf("stability: %d parameter/seed combinations finite + nonnegative + normalized\n", sweepRuns);

    // Long extreme runs remain within the intended one/two-morph vocabulary;
    // the exposed environment range leaves room for both peaks away from a wall.
    float closestEdge = 1.f;
    for (float niche : {FinchesField::nicheMin(), FinchesField::nicheMax()})
        for (float environment : {FinchesField::environmentMin(), 0.f, FinchesField::environmentMax()}) {
            FinchesField f;
            p = FinchesField::Parameters(); p.mutation = FinchesField::mutationMin();
            p.branching = FinchesField::branchingMax(); p.niche = niche; p.environment = environment;
            f.reset(environment);
            if (!evolve(f, p, 120.f, 0.16f) || significantPeakCount(f) > 2) {
                std::printf("FAIL: extreme regime produced invalid/three-plus morphology\n"); failures++;
            }
            if (f.metrics().split) {
                closestEdge = std::min(closestEdge, f.metrics().lowTrait - FinchesField::traitMin());
                closestEdge = std::min(closestEdge, FinchesField::traitMax() - f.metrics().highTrait);
            }
        }
    std::printf("extremes: closest accepted peak %.3f trait units from boundary\n", closestEdge);
    if (closestEdge < FinchesField::binWidth()) {
        std::printf("FAIL: exposed range pins an accepted peak to a boundary\n"); failures++;
    }

    // In a stabilizing regime the population follows the complete ENV range.
    float worstTracking = 0.f;
    for (float environment : {FinchesField::environmentMin(), -0.2f, 0.f, 0.2f,
                              FinchesField::environmentMax()}) {
        FinchesField f;
        p = FinchesField::Parameters(); p.branching = 0.65f; p.environment = environment;
        f.reset(-environment);
        evolve(f, p, 80.f, 0.13f);
        worstTracking = std::max(worstTracking, std::fabs(f.metrics().mean - environment));
    }
    std::printf("environment: worst settled mean error %.4f\n", worstTracking);
    if (worstTracking > 0.04f) {
        std::printf("FAIL: environmental optimum tracking\n"); failures++;
    }

    // No-flux boundaries retain injected mass and drive it inward rather than
    // wrapping it to the opposite edge.
    for (int side = 0; side < 2; ++side) {
        FinchesField edge;
        edge.reset(0.f);
        edge.seed(side ? 1.f : -1.f, 0.20f, 0.5f * FinchesField::binWidth());
        float oppositeBefore = edge.masses()[side ? 0 : FinchesField::kBins - 1];
        p = FinchesField::Parameters(); p.environment = side ? 0.5f : -0.5f;
        if (!evolve(edge, p, 4.f, 0.04f) ||
            edge.masses()[side ? 0 : FinchesField::kBins - 1] > oppositeBefore + 1e-6f) {
            std::printf("FAIL: boundary seed wrapped or destabilized\n"); failures++;
        }
    }

    // Outer call granularity must converge despite deliberately incommensurate
    // internal substeps.
    FinchesField coarse, fine;
    p = FinchesField::Parameters(); p.branching = 2.1f;
    evolve(coarse, p, 90.f, 0.137f);
    evolve(fine, p, 90.f, 0.0065f);
    float timestepError = maxMassDifference(coarse, fine);
    std::printf("timestep: max bin error %.7f (outer 0.137 vs 0.0065 tau)\n", timestepError);
    if (timestepError > 2e-4f || coarse.metrics().split != fine.metrics().split) {
        std::printf("FAIL: timestep convergence\n"); failures++;
    }

    // Event persistence is measured at the field's internal integration
    // cadence. One bounded outer call and the equivalent sequence of 0.02-tau
    // calls must therefore agree even when a qualifying morphology appears
    // partway through the interval. These fixtures begin eight internal steps
    // before the strong/weak detector boundary, leaving only 0.16 tau of the
    // 0.32-tau interval eligible for either event.
    FinchesField splitApproach;
    p = FinchesField::Parameters(); p.mutation = 0.00002f; p.branching = 2.8f;
    std::array<FinchesField::State, 9> splitHistory;
    int splitHistoryCount = 0;
    FinchesField::State beforeSplitCandidate;
    bool foundSplitCandidate = false;
    for (int step = 0; step < 20000 && !foundSplitCandidate; ++step) {
        const FinchesField::State before = splitApproach.state();
        if (splitHistoryCount < static_cast<int>(splitHistory.size()))
            splitHistory[splitHistoryCount++] = before;
        else {
            for (std::size_t i = 1; i < splitHistory.size(); ++i)
                splitHistory[i - 1] = splitHistory[i];
            splitHistory.back() = before;
        }
        splitApproach.advance(0.02f, p);
        if (splitApproach.state().splitTimer > 0.f
            && splitHistoryCount == static_cast<int>(splitHistory.size())) {
            beforeSplitCandidate = splitHistory.front();
            foundSplitCandidate = true;
        }
    }

    FinchesField coarseSplitCandidate, fineSplitCandidate;
    bool fineSplitEvent = false;
    bool splitGranularityOk = foundSplitCandidate
        && coarseSplitCandidate.restore(beforeSplitCandidate)
        && fineSplitCandidate.restore(beforeSplitCandidate);
    if (splitGranularityOk) {
        coarseSplitCandidate.advance(0.32f, p);
        for (int step = 0; step < 16; ++step) {
            fineSplitCandidate.advance(0.02f, p);
            fineSplitEvent = fineSplitEvent || fineSplitCandidate.metrics().splitEvent;
        }
        splitGranularityOk = !coarseSplitCandidate.metrics().splitEvent
            && !fineSplitEvent
            && sameState(coarseSplitCandidate.state(), fineSplitCandidate.state());
    }
    if (!splitGranularityOk) {
        std::printf("FAIL: coarse advance spuriously qualified SPLIT\n"); failures++;
    }
    else {
        // Continue from the identical pending state until the event occurs.
        // The coarse call must retain an event raised before its final substep.
        FinchesField coarseSplitEvent, fineSplitEventField;
        const FinchesField::State pending = coarseSplitCandidate.state();
        bool positiveFineSplitEvent = false;
        bool splitLatchOk = coarseSplitEvent.restore(pending)
            && fineSplitEventField.restore(pending);
        if (splitLatchOk) {
            coarseSplitEvent.advance(0.32f, p);
            for (int step = 0; step < 16; ++step) {
                fineSplitEventField.advance(0.02f, p);
                positiveFineSplitEvent = positiveFineSplitEvent
                    || fineSplitEventField.metrics().splitEvent;
            }
            splitLatchOk = coarseSplitEvent.metrics().splitEvent
                && positiveFineSplitEvent
                && sameState(coarseSplitEvent.state(), fineSplitEventField.state());
        }
        if (!splitLatchOk) {
            std::printf("FAIL: coarse advance did not latch SPLIT\n"); failures++;
        }
    }

    FinchesField mergeApproach;
    p = FinchesField::Parameters(); p.mutation = 0.00002f; p.branching = 2.8f;
    evolve(mergeApproach, p, 90.f, 0.13f);
    p.branching = 0.65f;
    std::array<FinchesField::State, 9> mergeHistory;
    int mergeHistoryCount = 0;
    FinchesField::State beforeMergeCandidate;
    bool foundMergeCandidate = false;
    for (int step = 0; step < 20000 && !foundMergeCandidate; ++step) {
        const FinchesField::State before = mergeApproach.state();
        if (mergeHistoryCount < static_cast<int>(mergeHistory.size()))
            mergeHistory[mergeHistoryCount++] = before;
        else {
            for (std::size_t i = 1; i < mergeHistory.size(); ++i)
                mergeHistory[i - 1] = mergeHistory[i];
            mergeHistory.back() = before;
        }
        mergeApproach.advance(0.02f, p);
        if (mergeApproach.state().mergeTimer > 0.f
            && mergeHistoryCount == static_cast<int>(mergeHistory.size())) {
            beforeMergeCandidate = mergeHistory.front();
            foundMergeCandidate = true;
        }
    }

    FinchesField coarseMergeCandidate, fineMergeCandidate;
    bool fineMergeEvent = false;
    bool mergeGranularityOk = foundMergeCandidate
        && coarseMergeCandidate.restore(beforeMergeCandidate)
        && fineMergeCandidate.restore(beforeMergeCandidate);
    if (mergeGranularityOk) {
        coarseMergeCandidate.advance(0.32f, p);
        for (int step = 0; step < 16; ++step) {
            fineMergeCandidate.advance(0.02f, p);
            fineMergeEvent = fineMergeEvent || fineMergeCandidate.metrics().mergeEvent;
        }
        mergeGranularityOk = !coarseMergeCandidate.metrics().mergeEvent
            && !fineMergeEvent
            && sameState(coarseMergeCandidate.state(), fineMergeCandidate.state());
    }
    if (!mergeGranularityOk) {
        std::printf("FAIL: coarse advance spuriously qualified MERGE\n"); failures++;
    }
    else {
        FinchesField coarseMergeEvent, fineMergeEventField;
        const FinchesField::State pending = coarseMergeCandidate.state();
        bool positiveFineMergeEvent = false;
        bool mergeLatchOk = coarseMergeEvent.restore(pending)
            && fineMergeEventField.restore(pending);
        if (mergeLatchOk) {
            coarseMergeEvent.advance(0.32f, p);
            for (int step = 0; step < 16; ++step) {
                fineMergeEventField.advance(0.02f, p);
                positiveFineMergeEvent = positiveFineMergeEvent
                    || fineMergeEventField.metrics().mergeEvent;
            }
            mergeLatchOk = coarseMergeEvent.metrics().mergeEvent
                && positiveFineMergeEvent
                && sameState(coarseMergeEvent.state(), fineMergeEventField.state());
        }
        if (!mergeLatchOk) {
            std::printf("FAIL: coarse advance did not latch MERGE\n"); failures++;
        }
    }

    // Reset and identical commands are deterministic bin-for-bin.
    FinchesField deterministicA, deterministicB;
    p = FinchesField::Parameters(); p.environment = 0.17f; p.branching = 2.3f;
    deterministicA.reset(p.environment); deterministicB.reset(p.environment);
    deterministicA.seed(-0.21f, 0.07f); deterministicB.seed(-0.21f, 0.07f);
    evolve(deterministicA, p, 18.f, 0.07f);
    evolve(deterministicB, p, 18.f, 0.07f);
    if (std::memcmp(deterministicA.masses().data(), deterministicB.masses().data(),
                    sizeof(float) * FinchesField::kBins) != 0) {
        std::printf("FAIL: deterministic replay differs\n"); failures++;
    }

    // The compatibility loader for pre-v1 density-only patches normalizes and
    // derives a strong split baseline without manufacturing an event.
    FinchesField savedSplit;
    p = FinchesField::Parameters(); p.mutation = 0.00002f; p.branching = 2.8f;
    evolve(savedSplit, p, 90.f, 0.13f);
    if (!savedSplit.metrics().split) {
        std::printf("FAIL: restore fixture did not split\n"); failures++;
    }
    FinchesField restored;
    float saved[FinchesField::kBins];
    savedSplit.copyMasses(saved, FinchesField::kBins);
    for (int i = 0; i < FinchesField::kBins; ++i) saved[i] *= 3.7f;
    if (!restored.restoreMasses(saved, FinchesField::kBins) || !valid(restored) ||
        !restored.metrics().split || restored.metrics().splitEvent || restored.metrics().mergeEvent) {
        std::printf("FAIL: valid restore contract\n"); failures++;
    }
    float beforeBad[FinchesField::kBins];
    restored.copyMasses(beforeBad, FinchesField::kBins);
    saved[7] = std::numeric_limits<float>::quiet_NaN();
    if (restored.restoreMasses(saved, FinchesField::kBins) ||
        std::memcmp(beforeBad, restored.masses().data(), sizeof(beforeBad)) != 0) {
        std::printf("FAIL: invalid restore was not rejected atomically\n"); failures++;
    }

    // A complete state round-trip preserves an accepted split while its peaks
    // are only in the weak hysteresis band. The old density-only inference is
    // deliberately used to discover that reachable band without duplicating
    // the detector thresholds in this test.
    FinchesField hysteresis;
    p = FinchesField::Parameters(); p.mutation = 0.00002f; p.branching = 2.8f;
    evolve(hysteresis, p, 90.f, 0.13f);
    p.branching = 0.65f;
    bool foundWeakBand = false;
    for (int step = 0; step < 2000 && hysteresis.metrics().split; ++step) {
        hysteresis.advance(0.01f, p);
        const FinchesField::State held = hysteresis.state();
        FinchesField densityOnly;
        if (densityOnly.restoreMasses(held.mass.data(), FinchesField::kBins)
            && !densityOnly.metrics().split && hysteresis.metrics().split) {
            FinchesField resumed;
            if (!resumed.restore(held) || !sameState(resumed.state(), held)
                || !resumed.metrics().split
                || resumed.metrics().splitEvent || resumed.metrics().mergeEvent
                || resumed.metrics().lowMass != hysteresis.metrics().lowMass
                || resumed.metrics().highMass != hysteresis.metrics().highMass
                || resumed.metrics().lowTrait != hysteresis.metrics().lowTrait
                || resumed.metrics().highTrait != hysteresis.metrics().highTrait) {
                std::printf("FAIL: weak-band split state round-trip\n"); failures++;
            }
            foundWeakBand = true;
            break;
        }
    }
    if (!foundWeakBand) {
        std::printf("FAIL: no reachable accepted weak-band split fixture\n"); failures++;
    }

    // Pending split and merge timers are musical state: a restored field must
    // emit the next event on the same future call as the uninterrupted field.
    FinchesField pendingSplit;
    p = FinchesField::Parameters(); p.mutation = 0.00002f; p.branching = 2.8f;
    FinchesField::State pendingSplitState;
    bool foundPendingSplit = false;
    for (int step = 0; step < 20000 && !pendingSplit.metrics().split; ++step) {
        pendingSplit.advance(0.01f, p);
        if (pendingSplit.state().splitTimer > 0.f) {
            pendingSplitState = pendingSplit.state();
            foundPendingSplit = true;
            break;
        }
    }
    FinchesField resumedSplit;
    if (!foundPendingSplit || !resumedSplit.restore(pendingSplitState)
        || !sameState(resumedSplit.state(), pendingSplitState)
        || resumedSplit.metrics().splitEvent || resumedSplit.metrics().mergeEvent) {
        std::printf("FAIL: pending split timer state round-trip\n"); failures++;
    }
    else {
        bool sawSplit = false;
        for (int step = 0; step < 100 && !sawSplit; ++step) {
            pendingSplit.advance(0.01f, p);
            resumedSplit.advance(0.01f, p);
            if (!sameState(pendingSplit.state(), resumedSplit.state())
                || pendingSplit.metrics().splitEvent != resumedSplit.metrics().splitEvent) {
                std::printf("FAIL: pending split future diverged\n"); failures++;
                break;
            }
            sawSplit = pendingSplit.metrics().splitEvent;
        }
        if (!sawSplit) {
            std::printf("FAIL: pending split fixture never emitted\n"); failures++;
        }
    }

    FinchesField pendingMerge;
    p = FinchesField::Parameters(); p.mutation = 0.00002f; p.branching = 2.8f;
    evolve(pendingMerge, p, 90.f, 0.13f);
    p.branching = 0.65f;
    FinchesField::State pendingMergeState;
    bool foundPendingMerge = false;
    for (int step = 0; step < 20000 && pendingMerge.metrics().split; ++step) {
        pendingMerge.advance(0.01f, p);
        if (pendingMerge.state().mergeTimer > 0.f) {
            pendingMergeState = pendingMerge.state();
            foundPendingMerge = true;
            break;
        }
    }
    FinchesField resumedMerge;
    if (!foundPendingMerge || !resumedMerge.restore(pendingMergeState)
        || !sameState(resumedMerge.state(), pendingMergeState)
        || resumedMerge.metrics().splitEvent || resumedMerge.metrics().mergeEvent) {
        std::printf("FAIL: pending merge timer state round-trip\n"); failures++;
    }
    else {
        bool sawMerge = false;
        for (int step = 0; step < 100 && !sawMerge; ++step) {
            pendingMerge.advance(0.01f, p);
            resumedMerge.advance(0.01f, p);
            if (!sameState(pendingMerge.state(), resumedMerge.state())
                || pendingMerge.metrics().mergeEvent != resumedMerge.metrics().mergeEvent) {
                std::printf("FAIL: pending merge future diverged\n"); failures++;
                break;
            }
            sawMerge = pendingMerge.metrics().mergeEvent;
        }
        if (!sawMerge) {
            std::printf("FAIL: pending merge fixture never emitted\n"); failures++;
        }
    }

    // Complete-state validation is transactional.
    const FinchesField::State beforeBadState = restored.state();
    FinchesField::State badState = beforeBadState;
    badState.version++;
    if (restored.restore(badState) || !sameState(restored.state(), beforeBadState)) {
        std::printf("FAIL: bad state version was not rejected atomically\n"); failures++;
    }
    badState = beforeBadState;
    badState.splitTimer = std::numeric_limits<float>::quiet_NaN();
    if (restored.restore(badState) || !sameState(restored.state(), beforeBadState)) {
        std::printf("FAIL: bad detector timer was not rejected atomically\n"); failures++;
    }
    badState = beforeBadState;
    const float displaced = badState.mass[0];
    badState.mass[0] = 1e-35f;
    badState.mass[FinchesField::kBins / 2] += displaced - 1e-35f;
    if (restored.restore(badState) || !sameState(restored.state(), beforeBadState)) {
        std::printf("FAIL: below-floor state was not rejected atomically\n"); failures++;
    }

    // An event is one-call state: seed and a zero-time advance cannot replay it.
    FinchesField eventLatch;
    p = FinchesField::Parameters(); p.mutation = 0.00002f; p.branching = 2.8f;
    int eventCount = 0;
    evolve(eventLatch, p, 90.f, 0.11f, &eventCount, 0);
    if (eventCount != 1 || !eventLatch.metrics().split) {
        std::printf("FAIL: event latch fixture\n"); failures++;
    }
    eventLatch.seed(0.8f);
    if (eventLatch.metrics().splitEvent || eventLatch.metrics().mergeEvent) {
        std::printf("FAIL: seed replayed stale event\n"); failures++;
    }
    eventLatch.advance(0.f, p);
    if (eventLatch.metrics().splitEvent || eventLatch.metrics().mergeEvent) {
        std::printf("FAIL: zero-time update replayed stale event\n"); failures++;
    }

    // NaN/inf/out-of-range inputs are sanitized before expensive operations.
    FinchesField hostile;
    p.mutation = std::numeric_limits<float>::quiet_NaN();
    p.branching = std::numeric_limits<float>::infinity();
    p.niche = -100.f; p.environment = 100.f;
    hostile.seed(std::numeric_limits<float>::infinity(), 100.f,
                 std::numeric_limits<float>::quiet_NaN());
    hostile.advance(std::numeric_limits<float>::infinity(), p);
    hostile.advance(-1.f, p);
    if (!valid(hostile)) {
        std::printf("FAIL: hostile input sanitization\n"); failures++;
    }

    // Edge-adapted populations used to retain float subnormals indefinitely.
    // This is reachable from reset alone and remains clean after a long run.
    FinchesField noSubnormal;
    p = FinchesField::Parameters();
    p.mutation = FinchesField::mutationMin();
    p.branching = FinchesField::branchingMax();
    p.niche = FinchesField::nicheMin();
    p.environment = FinchesField::environmentMax();
    noSubnormal.reset(p.environment);
    const bool resetHadNoTinyMass = valid(noSubnormal);
    evolve(noSubnormal, p, 300.f, 0.16f);
    int zeroBins = 0;
    float smallestPositive = 1.f;
    for (float value : noSubnormal.masses()) {
        if (value == 0.f)
            ++zeroBins;
        else
            smallestPositive = std::min(smallestPositive, value);
    }
    std::printf("numerics: %d extinct bins, smallest positive mass %.3g\n",
                zeroBins, smallestPositive);
    if (!resetHadNoTinyMass || !valid(noSubnormal) || zeroBins == 0
        || smallestPositive < std::numeric_limits<float>::min()) {
        std::printf("FAIL: numerical extinction floor / no-subnormal contract\n"); failures++;
    }

    // Diagnostic core-only cost at the default RATE and maximum bounded field
    // request. Hardware varies, so the threshold only catches gross regressions.
    const int defaultCalls = 2000;
    FinchesField defaultBenchmark;
    p = FinchesField::Parameters();
    const std::chrono::steady_clock::time_point defaultStart =
        std::chrono::steady_clock::now();
    for (int call = 0; call < defaultCalls; ++call)
        defaultBenchmark.advance(8.f / 500.f, p);
    const std::chrono::steady_clock::time_point defaultStop =
        std::chrono::steady_clock::now();
    const double defaultMicroseconds =
        std::chrono::duration<double, std::micro>(defaultStop - defaultStart).count()
        / defaultCalls;

    const int maximumCalls = 500;
    FinchesField maximumBenchmark;
    p.mutation = 0.00002f;
    p.branching = FinchesField::branchingMax();
    p.niche = FinchesField::nicheMin();
    p.environment = FinchesField::environmentMax();
    const float maximumAdvance = 8.f * std::exp2(4.25f) / 500.f;
    const std::chrono::steady_clock::time_point maximumStart =
        std::chrono::steady_clock::now();
    for (int call = 0; call < maximumCalls; ++call)
        maximumBenchmark.advance(maximumAdvance, p);
    const std::chrono::steady_clock::time_point maximumStop =
        std::chrono::steady_clock::now();
    const double maximumMicroseconds =
        std::chrono::duration<double, std::micro>(maximumStop - maximumStart).count()
        / maximumCalls;
    std::printf("benchmark: %.2f us default-rate / %.2f us maximum field advance "
                "(%.2f%% / %.2f%% of one core at 500 Hz)\n",
                defaultMicroseconds, maximumMicroseconds,
                defaultMicroseconds / 20.0, maximumMicroseconds / 20.0);
    if (!valid(defaultBenchmark) || !valid(maximumBenchmark)
        || defaultMicroseconds > 100.0 || maximumMicroseconds > 1000.0) {
        std::printf("FAIL: worst-case benchmark or terminal state\n"); failures++;
    }

    if (failures) {
        std::printf("FAIL: Finches stability contract (%d failures)\n", failures);
        return 1;
    }
    std::printf("PASS: Finches field stability and branching contract\n");
    return 0;
}
