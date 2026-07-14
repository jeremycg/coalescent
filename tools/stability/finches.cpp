// SDK-free stability and behaviour contract for src/dsp/finches_field.hpp.
//
//   g++ -std=c++11 -O2 tools/stability/finches.cpp -o /tmp/finches && /tmp/finches
#include "../../src/dsp/finches_field.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>

using coalescent::FinchesField;

static bool valid(const FinchesField& field, float tolerance = 3e-6f) {
    double sum = 0.0;
    const std::array<float, FinchesField::kBins>& m = field.masses();
    for (int i = 0; i < FinchesField::kBins; ++i) {
        if (!std::isfinite(m[i]) || m[i] < 0.f || m[i] > 1.f)
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

    // Valid restore is normalized, preserves a split baseline without an event,
    // and invalid payloads fail atomically.
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

    if (failures) {
        std::printf("FAIL: Finches stability contract (%d failures)\n", failures);
        return 1;
    }
    std::printf("PASS: Finches field stability and branching contract\n");
    return 0;
}
