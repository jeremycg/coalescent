// SDK-free stochastic and state contract for src/dsp/islands_model.hpp.
//
//   g++ -std=c++11 -O2 tools/stability/islands.cpp -o /tmp/islands && /tmp/islands
#include "../../src/dsp/hex64.hpp"
#include "../../src/dsp/islands_model.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <limits>

using coalescent::IslandsModel;
using coalescent::Hex64Codec;

static int failures = 0;

static void fail(const char* message) {
    std::printf("FAIL: %s\n", message);
    ++failures;
}

static bool sameState(const IslandsModel::State& a, const IslandsModel::State& b) {
    return a.counts == b.counts && a.denominators == b.denominators &&
           a.rngState == b.rngState && a.rngIncrement == b.rngIncrement &&
           a.nextFounder == b.nextFounder && a.wasFixA == b.wasFixA &&
           a.wasFixB == b.wasFixB && a.sweepAArmed == b.sweepAArmed &&
           a.sweepBArmed == b.sweepBArmed;
}

static bool valid(const IslandsModel& model) {
    const IslandsModel::State s = model.state();
    const IslandsModel::Metrics& m = model.metrics();
    double mean = 0.0;
    double diversity = 0.0;
    for (int i = 0; i < IslandsModel::kIslands; ++i) {
        if (s.denominators[i] < IslandsModel::copiesMin() ||
            s.denominators[i] > IslandsModel::copiesMax() ||
            s.counts[i] > s.denominators[i] ||
            !std::isfinite(m.frequency[i]) || m.frequency[i] < 0.0 ||
            m.frequency[i] > 1.0)
            return false;
        const double p = static_cast<double>(s.counts[i]) / s.denominators[i];
        if (std::fabs(m.frequency[i] - p) > 1e-15)
            return false;
        mean += p;
        diversity += p * (1.0 - p);
    }
    mean /= IslandsModel::kIslands;
    return std::isfinite(m.mean) && std::isfinite(m.diversity) &&
           std::fabs(m.mean - mean) < 1e-15 &&
           std::fabs(m.diversity - diversity) < 1e-15 &&
           m.diversity >= 0.0 && m.diversity <= 1.0;
}

static IslandsModel::State stateAt(const std::array<std::uint32_t, 4>& counts,
                                   const std::array<std::uint32_t, 4>& denominators,
                                   bool sweepAArmed, bool sweepBArmed) {
    IslandsModel source(987u, 123u);
    IslandsModel::State s = source.state();
    s.counts = counts;
    s.denominators = denominators;
    bool fixA = true;
    bool fixB = true;
    for (int i = 0; i < IslandsModel::kIslands; ++i) {
        fixA = fixA && counts[i] == denominators[i];
        fixB = fixB && counts[i] == 0u;
    }
    s.wasFixA = fixA;
    s.wasFixB = fixB;
    s.sweepAArmed = sweepAArmed;
    s.sweepBArmed = sweepBArmed;
    return s;
}

int main() {
    // Reference PCG-XSH-RR sequence from the canonical seed/stream example.
    IslandsModel::Pcg32 rng(42u, 54u);
    const std::uint32_t reference[] = {
        UINT32_C(0xa15c02b7), UINT32_C(0x7b47f409), UINT32_C(0xba1d3330),
        UINT32_C(0x83d2f293), UINT32_C(0xbfa4784b), UINT32_C(0xcbed606e)
    };
    for (std::size_t i = 0; i < sizeof(reference) / sizeof(reference[0]); ++i)
        if (rng.next() != reference[i])
            fail("PCG32 reference sequence changed");

    const std::uint64_t hexReference = UINT64_C(0xfedcba9876543210);
    char hexText[Hex64Codec::TEXT_SIZE];
    Hex64Codec::format(hexReference, hexText);
    std::uint64_t hexRoundTrip = 0u;
    if (!Hex64Codec::parseCString(hexText, hexRoundTrip) ||
        hexRoundTrip != hexReference)
        fail("shared Hex64 codec changed Islands persistence bits");

    IslandsModel model;
    if (!valid(model) || model.metrics().mean != 0.5 ||
        model.metrics().diversity != 1.0)
        fail("default 32/64 state or metrics");

    // Merely constructing and editing Parameters cannot redraw the model.
    const IslandsModel::State untouched = model.state();
    IslandsModel::Parameters p;
    p.copies = IslandsModel::copiesMax();
    p.selection = IslandsModel::selectionMax();
    p.mutation = IslandsModel::mutationMax();
    p.migration = IslandsModel::migrationMax();
    if (!sameState(untouched, model.state()))
        fail("parameter changes redrew population without an event");

    // Seed contract: equal seed/stream pairs reproduce every population draw.
    // Seed and stream are varied independently so neither input can silently
    // disappear from the implementation.
    const std::uint64_t replaySeed = UINT64_C(0x123456789abcdef0);
    const std::uint64_t replayStream = UINT64_C(0x13579bdf);
    const std::uint64_t alternateSeed = UINT64_C(0x0fedcba987654321);
    const std::uint64_t alternateStream = UINT64_C(0x2468ace0);
    IslandsModel seededA(replaySeed, replayStream);
    IslandsModel seededB(replaySeed, replayStream);
    IslandsModel differentSeed(alternateSeed, replayStream);
    IslandsModel differentStream(replaySeed, alternateStream);
    if (seededA.state().rngIncrement != ((replayStream << 1u) | 1u)
        || differentStream.state().rngIncrement != ((alternateStream << 1u) | 1u))
        fail("PCG stream selector did not set the expected increment");
    IslandsModel::Parameters seedParameters;
    seedParameters.copies = 257u;
    seedParameters.selection = 0.04;
    seedParameters.mutation = 0.012;
    seedParameters.migration = 0.31;
    bool seedDiverged = false;
    bool streamDiverged = false;
    for (int generation = 0; generation < 12; ++generation) {
        if (generation == 4) {
            seededA.founder();
            seededB.founder();
            differentSeed.founder();
            differentStream.founder();
        }
        seededA.advance(seedParameters);
        seededB.advance(seedParameters);
        differentSeed.advance(seedParameters);
        differentStream.advance(seedParameters);
        if (!sameState(seededA.state(), seededB.state()))
            fail("same seed/stream did not replay exactly");
        if (seededA.state().counts != differentSeed.state().counts)
            seedDiverged = true;
        if (seededA.state().counts != differentStream.state().counts)
            streamDiverged = true;
    }
    if (!seedDiverged)
        fail("distinct seed did not diverge in population state");
    if (!streamDiverged)
        fail("distinct stream did not diverge in population state");

    // Resetting to an explicitly chosen pair returns counts, founder cursor,
    // event latches, and RNG to the start of the same stochastic sequence.
    const std::uint64_t resetSeed = UINT64_C(0x1122334455667788);
    const std::uint64_t resetStream = UINT64_C(0x102030405060708);
    IslandsModel resetReplay(resetSeed, resetStream);
    std::array<IslandsModel::State, 14> expectedAfterReset;
    for (int pass = 0; pass < 2; ++pass) {
        if (pass == 1)
            resetReplay.reset(resetSeed, resetStream);
        for (int generation = 0; generation < 14; ++generation) {
            seedParameters.copies = 64u + static_cast<std::uint32_t>(generation * 37);
            seedParameters.selection = (generation % 3 - 1) * 0.08;
            seedParameters.mutation = 0.004 * (generation % 5);
            seedParameters.migration = 0.1 * (generation % 6);
            if (generation == 3 || generation == 10)
                resetReplay.founder();
            resetReplay.advance(seedParameters);
            if (pass == 0)
                expectedAfterReset[generation] = resetReplay.state();
            else if (!sameState(expectedAfterReset[generation], resetReplay.state()))
                fail("reset to chosen seed/stream did not replay sequence");
        }
    }

    // Founder targets rotate 0,1,2,3 and only the selected denominator changes.
    IslandsModel founder(77u, 91u);
    for (std::uint32_t target = 0u; target < 8u; ++target) {
        const IslandsModel::State before = founder.state();
        founder.founder();
        const IslandsModel::State after = founder.state();
        const std::uint32_t island = target % 4u;
        if (after.denominators[island] != IslandsModel::founderCopies() ||
            after.counts[island] > IslandsModel::founderCopies() ||
            after.nextFounder != (island + 1u) % 4u)
            fail("founder target, size, or cursor");
        for (std::uint32_t other = 0u; other < 4u; ++other)
            if (other != island &&
                (after.counts[other] != before.counts[other] ||
                 after.denominators[other] != before.denominators[other]))
                fail("founder altered a non-target island");
    }

    // Identical seeds and method order are bit-exact; restoring mid-sequence
    // resumes the same stochastic future and does not manufacture an event.
    IslandsModel replayA(123456u, 98765u), replayB(123456u, 98765u);
    for (int generation = 0; generation < 40; ++generation) {
        p.copies = 8u + static_cast<std::uint32_t>((generation * 97) % 4089);
        p.selection = -0.25 + 0.5 * (generation % 11) / 10.0;
        p.mutation = 0.05 * (generation % 7) / 6.0;
        p.migration = (generation % 9) / 8.0;
        if (generation == 9 || generation == 23) {
            replayA.founder();
            replayB.founder();
        }
        replayA.advance(p);
        replayB.advance(p);
        if (!sameState(replayA.state(), replayB.state()))
            fail("deterministic replay diverged");
    }
    IslandsModel resumed(1u, 1u);
    if (!resumed.restore(replayA.state()) || resumed.metrics().lossEvent ||
        resumed.metrics().sweepEvent)
        fail("valid restore failed or emitted an event");
    for (int generation = 0; generation < 20; ++generation) {
        p.copies = 257u + generation * 31u;
        p.selection = generation & 1 ? -0.11 : 0.07;
        p.mutation = 0.013;
        p.migration = 0.37;
        replayA.advance(p);
        resumed.advance(p);
        if (!sameState(replayA.state(), resumed.state()))
            fail("restored RNG/state did not resume exactly");
    }

    // Invalid state is rejected atomically.
    const IslandsModel::State beforeBadRestore = resumed.state();
    IslandsModel::State bad = beforeBadRestore;
    bad.counts[0] = bad.denominators[0] + 1u;
    if (resumed.restore(bad) || !sameState(beforeBadRestore, resumed.state()))
        fail("invalid count restore was not atomic");
    bad = beforeBadRestore;
    bad.rngIncrement &= ~UINT64_C(1);
    if (resumed.restore(bad) || !sameState(beforeBadRestore, resumed.state()))
        fail("invalid PCG stream restore was accepted");
    bad = beforeBadRestore;
    bad.nextFounder = 4u;
    if (resumed.restore(bad) || !sameState(beforeBadRestore, resumed.state()))
        fail("invalid founder cursor restore was accepted");

    // Genic log-fitness selection has the expected direction and magnitude.
    double selectedUp = 0.0;
    double selectedDown = 0.0;
    for (std::uint64_t seed = 0u; seed < 64u; ++seed) {
        IslandsModel up(1000u + seed, 2000u + seed);
        IslandsModel down(1000u + seed, 2000u + seed);
        p = IslandsModel::Parameters();
        p.copies = IslandsModel::copiesMax();
        p.selection = IslandsModel::selectionMax();
        up.advance(p);
        p.selection = IslandsModel::selectionMin();
        down.advance(p);
        selectedUp += up.metrics().mean;
        selectedDown += down.metrics().mean;
    }
    selectedUp /= 64.0;
    selectedDown /= 64.0;
    const double expectedUp = std::exp(0.25) / (1.0 + std::exp(0.25));
    std::printf("selection: means %.5f / %.5f (expected %.5f / %.5f)\n",
                selectedDown, selectedUp, 1.0 - expectedUp, expectedUp);
    if (std::fabs(selectedUp - expectedUp) > 0.003 ||
        std::fabs(selectedDown - (1.0 - expectedUp)) > 0.003)
        fail("selection bias does not match genic update");

    // Symmetric mutation recreates a lost allele; migration homogenizes the
    // parental probability but each island still receives an independent draw.
    IslandsModel mutation(22u, 33u);
    IslandsModel::State fixedB = stateAt({{0u, 0u, 0u, 0u}},
                                        {{4096u, 4096u, 4096u, 4096u}}, true, false);
    if (!mutation.restore(fixedB))
        fail("could not prepare fixed-B mutation state");
    p = IslandsModel::Parameters();
    p.copies = 4096u;
    p.mutation = 0.05;
    mutation.advance(p);
    if (mutation.metrics().fixB || mutation.metrics().mean < 0.035 ||
        mutation.metrics().mean > 0.065 || mutation.metrics().lossEvent)
        fail("mutation did not recover from exact fixation");

    const IslandsModel::State divided = stateAt({{0u, 4096u, 0u, 4096u}},
        {{4096u, 4096u, 4096u, 4096u}}, true, true);
    IslandsModel isolated(55u, 66u), mixed(55u, 66u);
    if (!isolated.restore(divided) || !mixed.restore(divided))
        fail("could not prepare divided-island state");
    p = IslandsModel::Parameters(); p.copies = 4096u; p.migration = 0.0;
    isolated.advance(p);
    if (isolated.state().counts != divided.counts)
        fail("zero migration failed to preserve fixed isolated islands");
    p.migration = 1.0;
    mixed.advance(p);
    double maxMixedError = 0.0;
    for (int i = 0; i < 4; ++i)
        maxMixedError = std::max(maxMixedError,
                                 std::fabs(mixed.metrics().frequency[i] - 0.5));
    if (maxMixedError > 0.03 ||
        (mixed.state().counts[0] == mixed.state().counts[1] &&
         mixed.state().counts[1] == mixed.state().counts[2] &&
         mixed.state().counts[2] == mixed.state().counts[3]))
        fail("common-pool migration or independent sampling");

    // One-generation neutral variance must be p(1-p)/N.
    const int replicates = 6000;
    double sum = 0.0;
    double sumSquares = 0.0;
    p = IslandsModel::Parameters(); p.copies = 64u;
    for (int replicate = 0; replicate < replicates; ++replicate) {
        IslandsModel neutral(100000u + replicate, 70000u + replicate);
        neutral.advance(p);
        const double x = neutral.metrics().frequency[0];
        sum += x;
        sumSquares += x * x;
    }
    const double neutralMean = sum / replicates;
    const double neutralVariance = (sumSquares - sum * sum / replicates) / (replicates - 1);
    const double expectedVariance = 0.25 / 64.0;
    std::printf("neutral: mean %.5f variance %.7f (expected %.7f)\n",
                neutralMean, neutralVariance, expectedVariance);
    if (std::fabs(neutralMean - 0.5) > 0.004 ||
        std::fabs(neutralVariance / expectedVariance - 1.0) > 0.08)
        fail("neutral mean/variance calibration");

    // Sweep entry is one-shot until the opposite side of the hysteresis band
    // rearms it. Test both A and B directions.
    IslandsModel sweep(91u, 29u);
    IslandsModel::State nearA = stateAt({{89u, 89u, 89u, 89u}},
                                       {{100u, 100u, 100u, 100u}}, true, true);
    if (!sweep.restore(nearA))
        fail("could not prepare A sweep state");
    p = IslandsModel::Parameters(); p.copies = 4096u; p.selection = 0.25;
    sweep.advance(p);
    if (!sweep.metrics().sweepEvent || sweep.metrics().mean < 0.9)
        fail("A sweep entry did not fire");
    p.selection = 0.0;
    sweep.advance(p);
    if (sweep.metrics().sweepEvent)
        fail("A sweep repeated while above threshold");
    p.selection = -0.25;
    for (int i = 0; i < 12 && sweep.metrics().mean > 0.8; ++i)
        sweep.advance(p);
    if (!sweep.state().sweepAArmed)
        fail("A sweep did not rearm below 0.8");
    p.selection = 0.25;
    bool secondASweep = false;
    for (int i = 0; i < 16 && !secondASweep; ++i) {
        sweep.advance(p);
        secondASweep = sweep.metrics().sweepEvent;
    }
    if (!secondASweep)
        fail("rearmed A sweep did not fire again");

    IslandsModel sweepB(92u, 30u);
    IslandsModel::State nearB = stateAt({{11u, 11u, 11u, 11u}},
                                       {{100u, 100u, 100u, 100u}}, true, true);
    if (!sweepB.restore(nearB))
        fail("could not prepare B sweep state");
    p = IslandsModel::Parameters(); p.copies = 4096u; p.selection = -0.25;
    sweepB.advance(p);
    if (!sweepB.metrics().sweepEvent || sweepB.metrics().mean > 0.1)
        fail("B sweep entry did not fire");

    // LOSS fires once on entry to exact global fixation, and FIX is an exact
    // count gate rather than a smoothed threshold.
    IslandsModel loss(314u, 159u);
    IslandsModel::State nearlyLost = stateAt({{0u, 0u, 0u, 1u}},
        {{8u, 8u, 8u, 8u}}, true, false);
    if (!loss.restore(nearlyLost))
        fail("could not prepare loss state");
    p = IslandsModel::Parameters(); p.copies = 8u;
    int lossEvents = 0;
    for (int i = 0; i < 100 && !loss.metrics().fixB; ++i) {
        loss.advance(p);
        if (loss.metrics().lossEvent) ++lossEvents;
    }
    if (!loss.metrics().fixB || lossEvents != 1)
        fail("global B loss/fixation entry contract");
    loss.advance(p);
    if (loss.metrics().lossEvent || !loss.metrics().fixB)
        fail("LOSS repeated while fixation persisted");

    // Hostile parameter values sanitize to a finite, bounded transition.
    IslandsModel hostile(1u, 2u);
    p = IslandsModel::Parameters();
    p.copies = 0u;
    p.selection = std::numeric_limits<double>::quiet_NaN();
    p.mutation = std::numeric_limits<double>::infinity();
    p.migration = -std::numeric_limits<double>::infinity();
    hostile.advance(p);
    if (!valid(hostile) || hostile.state().denominators[0] != 8u)
        fail("hostile low/non-finite parameters");
    p.copies = std::numeric_limits<std::uint32_t>::max();
    p.selection = 1e300;
    p.mutation = -1e300;
    p.migration = 1e300;
    hostile.advance(p);
    if (!valid(hostile) || hostile.state().denominators[0] != 4096u)
        fail("hostile high parameters");

    // Full exposed parameter corners remain bounded and restorable.
    int cornerRuns = 0;
    for (std::uint32_t copies : {8u, 4096u})
        for (double selection : {-0.25, 0.25})
            for (double mutationRate : {0.0, 0.05})
                for (double migrationRate : {0.0, 1.0}) {
                    IslandsModel corner(500u + cornerRuns, 900u + cornerRuns);
                    p = IslandsModel::Parameters();
                    p.copies = copies;
                    p.selection = selection;
                    p.mutation = mutationRate;
                    p.migration = migrationRate;
                    for (int generation = 0; generation < 24; ++generation) {
                        corner.advance(p);
                        if (!valid(corner))
                            fail("full-range corner left valid state space");
                    }
                    IslandsModel copy;
                    if (!copy.restore(corner.state()) || !sameState(copy.state(), corner.state()))
                        fail("full-range corner state was not restorable");
                    ++cornerRuns;
                }
    std::printf("stability: %d parameter corners bounded and restorable\n", cornerRuns);

    // Worst-case cost is deliberately transparent: exact sampling consumes
    // 4*4096 Bernoulli draws per generation, but remains event-rate work.
    IslandsModel benchmark(123u, 456u);
    p = IslandsModel::Parameters();
    p.copies = 4096u; p.selection = 0.03; p.mutation = 0.01; p.migration = 0.2;
    const int benchmarkGenerations = 500;
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    for (int i = 0; i < benchmarkGenerations; ++i)
        benchmark.advance(p);
    const std::chrono::steady_clock::time_point stop = std::chrono::steady_clock::now();
    const double microseconds =
        std::chrono::duration<double, std::micro>(stop - start).count() /
        benchmarkGenerations;
    std::printf("benchmark: %.2f us/generation at four x 4096 copies (%.2f%% at 200 gen/s)\n",
                microseconds, microseconds * 200.0 / 10000.0);
    // Leave ample room for slower CI hardware while still rejecting a change
    // that turns this event-rate operation into an audio-buffer-scale stall.
    if (!valid(benchmark) || microseconds > 250.0)
        fail("worst-case benchmark or terminal state");

    if (failures) {
        std::printf("islands stability: %d failure(s)\n", failures);
        return 1;
    }
    std::printf("islands stability: PASS\n");
    return 0;
}
