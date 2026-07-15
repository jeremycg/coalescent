// SDK-free spatial trait-field contract for src/dsp/archipelago_field.hpp.
//
//   g++ -std=c++11 -O2 tools/stability/archipelago.cpp -o /tmp/archipelago
//   /tmp/archipelago
#include "../../src/dsp/archipelago_field.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <limits>

using coalescent::ArchipelagoField;

static int failures = 0;

static void fail(const char* message) {
    std::printf("FAIL: %s\n", message);
    ++failures;
}

static bool sameState(const ArchipelagoField::State& a,
                      const ArchipelagoField::State& b) {
    return a.density == b.density &&
           a.reportedTrait == b.reportedTrait &&
           a.reportedGlobalMean == b.reportedGlobalMean &&
           a.occupiedMask == b.occupiedMask;
}

static double stateTotal(const ArchipelagoField::State& state) {
    double total = 0.0;
    for (int i = 0; i < ArchipelagoField::kValues; ++i)
        total += state.density[i];
    return total;
}

static double maxDensityDifference(const ArchipelagoField& a,
                                   const ArchipelagoField& b) {
    const ArchipelagoField::State sa = a.state();
    const ArchipelagoField::State sb = b.state();
    double worst = 0.0;
    for (int i = 0; i < ArchipelagoField::kValues; ++i)
        worst = std::max(worst, std::fabs(sa.density[i] - sb.density[i]));
    return worst;
}

static bool valid(const ArchipelagoField& field) {
    const ArchipelagoField::State state = field.state();
    const ArchipelagoField::Metrics& metrics = field.metrics();
    double total = 0.0;
    for (int habitat = 0; habitat < ArchipelagoField::kHabitats; ++habitat) {
        double mass = 0.0;
        for (int bin = 0; bin < ArchipelagoField::kBins; ++bin) {
            const double value = state.density[habitat * ArchipelagoField::kBins + bin];
            if (!std::isfinite(value) || value < 0.0 || value > 4.0 ||
                (value > 0.0 &&
                 value < ArchipelagoField::numericalExtinctionFloor()))
                return false;
            mass += value;
        }
        total += mass;
        if (!std::isfinite(metrics.mass[habitat]) ||
            std::fabs(metrics.mass[habitat] - mass) > 1e-10 ||
            !std::isfinite(metrics.trait[habitat]) ||
            metrics.trait[habitat] < -1.0 || metrics.trait[habitat] > 1.0)
            return false;
    }
    const std::uint32_t mask = (UINT32_C(1) << ArchipelagoField::kHabitats) - 1u;
    return std::isfinite(total) && std::isfinite(metrics.totalMass) &&
           std::fabs(metrics.totalMass - total) < 1e-10 &&
           std::isfinite(metrics.globalMean) && metrics.globalMean >= -1.0 &&
           metrics.globalMean <= 1.0 && std::isfinite(metrics.difference) &&
           metrics.difference >= -2.0 && metrics.difference <= 2.0 &&
           std::isfinite(metrics.flux) && metrics.flux >= 0.0 &&
           metrics.flux <= 1.0 && metrics.occupiedMask == state.occupiedMask &&
           (metrics.occupiedMask & ~mask) == 0u &&
           (metrics.colonizeMask & ~mask) == 0u &&
           (metrics.extinctMask & ~mask) == 0u &&
           metrics.colonizeEvent == (metrics.colonizeMask != 0u) &&
           metrics.extinctEvent == (metrics.extinctMask != 0u);
}

static ArchipelagoField::State fixtureState() {
    ArchipelagoField::State state;
    state.density.fill(0.0);
    state.reportedTrait.fill(0.0);
    state.reportedGlobalMean = 0.0;
    state.occupiedMask = 0u;
    return state;
}

static void putMass(ArchipelagoField::State& state, int habitat, int bin,
                    double mass) {
    state.density[habitat * ArchipelagoField::kBins + bin] += mass;
}

static void finishFixture(ArchipelagoField::State& state) {
    state.occupiedMask = 0u;
    double totalMass = 0.0;
    double totalMoment = 0.0;
    for (int habitat = 0; habitat < ArchipelagoField::kHabitats; ++habitat) {
        double mass = 0.0;
        double moment = 0.0;
        for (int bin = 0; bin < ArchipelagoField::kBins; ++bin) {
            const double value = state.density[habitat * ArchipelagoField::kBins + bin];
            mass += value;
            moment += value * ArchipelagoField::traitAt(bin);
        }
        if (mass > 1e-12)
            state.reportedTrait[habitat] = moment / mass;
        if (mass >= ArchipelagoField::colonizeThreshold())
            state.occupiedMask |= UINT32_C(1) << habitat;
        totalMass += mass;
        totalMoment += moment;
    }
    state.reportedGlobalMean = totalMass > 1e-12 ? totalMoment / totalMass : 0.0;
}

static void evolve(ArchipelagoField& field, const ArchipelagoField::Parameters& p,
                   double totalTau, double callTau,
                   std::uint32_t* colonizeMask = 0,
                   std::uint32_t* extinctMask = 0) {
    double elapsed = 0.0;
    while (elapsed + 1e-14 < totalTau) {
        const double delta = std::min(callTau, totalTau - elapsed);
        field.advance(delta, p);
        if (colonizeMask)
            *colonizeMask |= field.metrics().colonizeMask;
        if (extinctMask)
            *extinctMask |= field.metrics().extinctMask;
        elapsed += delta;
    }
}

int main() {
    ArchipelagoField::Parameters p;
    ArchipelagoField field;
    if (!valid(field) || field.metrics().occupiedMask != 0xffu)
        fail("default reset state or occupancy");
    for (int habitat = 0; habitat < ArchipelagoField::kHabitats; ++habitat) {
        if (std::fabs(field.metrics().mass[habitat] - 0.75) > 1e-12)
            fail("reset mass is not 0.75 per habitat");
        if (habitat > 0 &&
            !(field.metrics().trait[habitat] > field.metrics().trait[habitat - 1]))
            fail("default reset does not form an increasing cline");
    }
    if (field.metrics().colonizeEvent || field.metrics().extinctEvent)
        fail("reset emitted an occupancy event");

    // Public metrics are derived from abundance rather than panel parameters.
    // This two-block fixture also gives FLUX one known discontinuous edge.
    ArchipelagoField::State metricFixture = fixtureState();
    for (int habitat = 0; habitat < ArchipelagoField::kHabitats; ++habitat)
        putMass(metricFixture, habitat, habitat < 4 ? 8 : 24, 1.0);
    finishFixture(metricFixture);
    ArchipelagoField metricField;
    p = ArchipelagoField::Parameters();
    p.selection = 0.0; p.mutation = 0.0; p.migration = 1.0;
    if (!metricField.restore(metricFixture, p))
        fail("could not restore metric fixture");
    const double lowTrait = ArchipelagoField::traitAt(8);
    const double highTrait = ArchipelagoField::traitAt(24);
    const double expectedGlobal = 0.5 * (lowTrait + highTrait);
    const double expectedFlux = 1.0 / (2.0 * 7.0);
    if (std::fabs(metricField.metrics().globalMean - expectedGlobal) > 1e-15 ||
        std::fabs(metricField.metrics().difference - (highTrait - lowTrait)) > 1e-15 ||
        std::fabs(metricField.metrics().flux - expectedFlux) > 1e-15)
        fail("global mean, signed difference, or state-dependent flux calibration");
    for (int habitat = 0; habitat < ArchipelagoField::kHabitats; ++habitat) {
        const double expectedTrait = habitat < 4 ? lowTrait : highTrait;
        if (metricField.metrics().mass[habitat] != 1.0 ||
            metricField.metrics().trait[habitat] != expectedTrait)
            fail("local mass or trait metric calibration");
    }
    float copied[ArchipelagoField::kValues];
    metricField.copyMasses(copied, ArchipelagoField::kValues);
    for (int i = 0; i < ArchipelagoField::kValues; ++i)
        if (copied[i] != static_cast<float>(metricFixture.density[i]))
            fail("display density copy changed habitat-major state");

    // Restore canonicalizes numerical dust without erasing normal standing
    // variation. Exact zeros remain structural zeros through a positive update.
    ArchipelagoField::State numericalFixture = fixtureState();
    numericalFixture.density[0] = std::numeric_limits<double>::denorm_min();
    numericalFixture.density[1] = 1e-30;
    finishFixture(numericalFixture);
    ArchipelagoField numericalField;
    p = ArchipelagoField::Parameters();
    p.selection = 0.0; p.mutation = 0.0; p.migration = 0.0;
    if (!numericalField.restore(numericalFixture, p))
        fail("could not restore numerical-floor fixture");
    ArchipelagoField::State numericalState = numericalField.state();
    if (numericalState.density[0] != 0.0 ||
        numericalState.density[1] != 1e-30 ||
        numericalState.density[2] != 0.0)
        fail("restore did not canonicalize dust or preserve standing variation");
    numericalField.advance(0.02, p);
    numericalState = numericalField.state();
    if (numericalState.density[0] != 0.0 ||
        !(numericalState.density[1] > 1e-30) ||
        numericalState.density[2] != 0.0 ||
        numericalField.metrics().colonizeEvent ||
        numericalField.metrics().extinctEvent || !valid(numericalField))
        fail("numerical floor changed a structural zero or observable state");
    p.migration = 0.0;
    metricField.advance(0.0, p);
    if (metricField.metrics().flux != 0.0)
        fail("zero migration did not produce zero flux");
    p.migration = 1.0; p.barrier = 1.0;
    metricField.advance(0.0, p);
    if (metricField.metrics().flux != 0.0)
        fail("closed central edge still contributed flux");

    // True-zero process settings preserve an equilibrium state bit-for-bit.
    ArchipelagoField::State equilibrium = fixtureState();
    for (int habitat = 0; habitat < ArchipelagoField::kHabitats; ++habitat)
        putMass(equilibrium, habitat, 3 + habitat * 3, 1.0);
    finishFixture(equilibrium);
    p = ArchipelagoField::Parameters();
    p.selection = 0.0; p.mutation = 0.0; p.migration = 0.0;
    ArchipelagoField zeroProcesses;
    if (!zeroProcesses.restore(equilibrium, p))
        fail("could not restore equilibrium fixture");
    const ArchipelagoField::State equilibriumBefore = zeroProcesses.state();
    evolve(zeroProcesses, p, 2.0, 0.25);
    if (!sameState(equilibriumBefore, zeroProcesses.state()))
        fail("zero mutation/migration altered a unit-mass equilibrium");

    // Mutation and migration separately conserve mass. Unit mass in every
    // habitat also neutralizes logistic regulation, isolating those operators.
    ArchipelagoField mutationOnly;
    p.mutation = ArchipelagoField::mutationMax();
    if (!mutationOnly.restore(equilibrium, p))
        fail("could not restore mutation conservation fixture");
    evolve(mutationOnly, p, 8.0, 0.137);
    if (!valid(mutationOnly) ||
        std::fabs(stateTotal(mutationOnly.state()) - 8.0) > 2e-10)
        fail("no-flux mutation failed to conserve global mass");
    for (int habitat = 0; habitat < ArchipelagoField::kHabitats; ++habitat)
        if (std::fabs(mutationOnly.metrics().mass[habitat] - 1.0) > 3e-11)
            fail("no-flux mutation changed local abundance");

    ArchipelagoField migrationOnly;
    p.mutation = 0.0; p.migration = ArchipelagoField::migrationMax(); p.ring = true;
    if (!migrationOnly.restore(equilibrium, p))
        fail("could not restore migration conservation fixture");
    evolve(migrationOnly, p, 8.0, 0.137);
    if (!valid(migrationOnly) ||
        std::fabs(stateTotal(migrationOnly.state()) - 8.0) > 2e-10)
        fail("ring migration failed to conserve global mass");
    for (int habitat = 0; habitat < ArchipelagoField::kHabitats; ++habitat)
        if (std::fabs(migrationOnly.metrics().mass[habitat] - 1.0) > 3e-11)
            fail("migration changed equal local abundances");

    // The wrap edge must be absent in ROW and present in RING. Closing the
    // central edge must make the two halves exactly independent in ROW.
    ArchipelagoField::State atRightEnd = fixtureState();
    putMass(atRightEnd, 7, 25, 1.0);
    finishFixture(atRightEnd);
    p = ArchipelagoField::Parameters();
    p.selection = 0.0; p.mutation = 0.0; p.migration = 2.0; p.ring = false;
    ArchipelagoField row, ring;
    if (!row.restore(atRightEnd, p))
        fail("could not restore row fixture");
    ArchipelagoField::Parameters ringParameters = p;
    ringParameters.ring = true;
    if (!ring.restore(atRightEnd, ringParameters))
        fail("could not restore ring fixture");
    row.advance(0.02, p);
    ring.advance(0.02, ringParameters);
    const double rowWrap = row.metrics().mass[0];
    const double ringWrap = ring.metrics().mass[0];
    if (!(ringWrap > rowWrap + 0.01))
        fail("ring wrap edge is absent or row has an unintended wrap");

    ArchipelagoField::State centralSource = fixtureState();
    putMass(centralSource, 3, 20, 1.0);
    finishFixture(centralSource);
    ArchipelagoField closedBarrier, openBarrier;
    p.ring = false; p.barrier = 1.0;
    if (!closedBarrier.restore(centralSource, p))
        fail("could not restore closed-barrier fixture");
    ArchipelagoField::Parameters open = p;
    open.barrier = 0.0;
    if (!openBarrier.restore(centralSource, open))
        fail("could not restore open-barrier fixture");
    evolve(closedBarrier, p, 0.5, 0.02);
    evolve(openBarrier, open, 0.5, 0.02);
    double closedRight = 0.0;
    double openRight = 0.0;
    for (int habitat = 4; habitat < ArchipelagoField::kHabitats; ++habitat) {
        closedRight += closedBarrier.metrics().mass[habitat];
        openRight += openBarrier.metrics().mass[habitat];
    }
    if (closedRight != 0.0 || !(openRight > 0.05))
        fail("central barrier is not exactly closed/open as specified");

    // Local selection maintains a geographic cline; strong gene flow pulls
    // the same starting field toward a more homogeneous set of local traits.
    ArchipelagoField::Parameters adaptedParameters;
    adaptedParameters.selection = 6.0;
    adaptedParameters.mutation = 1.0e-4;
    adaptedParameters.migration = 0.004;
    adaptedParameters.gradient = 0.70;
    adaptedParameters.climate = 0.0;
    ArchipelagoField adapted, homogenized;
    adapted.reset(adaptedParameters);
    homogenized.reset(adaptedParameters);
    ArchipelagoField::Parameters homogenizedParameters = adaptedParameters;
    homogenizedParameters.migration = 2.0;
    evolve(adapted, adaptedParameters, 24.0, 0.113);
    evolve(homogenized, homogenizedParameters, 24.0, 0.113);
    const double adaptedSpan = adapted.metrics().trait[7] - adapted.metrics().trait[0];
    const double homogenizedSpan =
        homogenized.metrics().trait[7] - homogenized.metrics().trait[0];
    std::printf("adaptation: endpoint span %.4f isolated / %.4f high migration\n",
                adaptedSpan, homogenizedSpan);
    if (!(adaptedSpan > 1.30 && homogenizedSpan < 0.80 * adaptedSpan))
        fail("selection/migration balance did not distinguish adaptation and homogenization");

    // An extreme shared climate makes the high-gradient range edge hostile.
    // Returning the climate restores growth from residual density or migrants.
    ArchipelagoField::Parameters climateParameters;
    climateParameters.selection = 8.0;
    climateParameters.mutation = 1.0e-4;
    climateParameters.migration = 0.03;
    climateParameters.gradient = 0.85;
    climateParameters.climate = 0.0;
    ArchipelagoField climate;
    climate.reset(climateParameters);
    climateParameters.climate = 0.85;
    std::uint32_t climateExtinctions = 0u;
    evolve(climate, climateParameters, 28.0, 0.10, 0, &climateExtinctions);
    const std::uint32_t contractedMask = climate.metrics().occupiedMask;
    climateParameters.climate = 0.0;
    std::uint32_t climateColonizations = 0u;
    evolve(climate, climateParameters, 40.0, 0.10, &climateColonizations, 0);
    std::printf("climate: contracted mask 0x%02x, extinct 0x%02x, recolonized 0x%02x\n",
                static_cast<unsigned>(contractedMask),
                static_cast<unsigned>(climateExtinctions),
                static_cast<unsigned>(climateColonizations));
    if (climateExtinctions == 0u || contractedMask == 0xffu ||
        (climateColonizations & climateExtinctions) == 0u)
        fail("climate shift did not contract and recolonize the range");

    // A climate reversal with no mutation or migration must still adapt from
    // the tiny Gaussian tail seeded by RESET. This bounds the numerical floor
    // well below biologically meaningful standing variation.
    ArchipelagoField::Parameters standingParameters;
    standingParameters.selection = 8.0;
    standingParameters.mutation = 0.0;
    standingParameters.migration = 0.0;
    standingParameters.gradient = 0.0;
    standingParameters.climate = -0.85;
    standingParameters.barrier = 0.0;
    standingParameters.ring = false;
    ArchipelagoField standingVariation;
    standingVariation.reset(standingParameters);
    standingParameters.climate = 0.85;
    std::uint32_t standingExtinctions = 0u;
    std::uint32_t standingColonizations = 0u;
    evolve(standingVariation, standingParameters, 400.0, 0.128,
           &standingColonizations, &standingExtinctions);
    if (standingExtinctions != 0xffu || standingColonizations != 0xffu ||
        standingVariation.metrics().occupiedMask != 0xffu ||
        standingVariation.metrics().mass[0] < 0.9 ||
        standingVariation.metrics().trait[0] < 0.7)
        fail("numerical floor erased standing-variation climate recovery");

    // Hysteresis holds TRAIT while empty, enters at 0.08, exits at 0.03, and
    // reports each transition for one advance call only.
    ArchipelagoField::State nearColonization = fixtureState();
    putMass(nearColonization, 0, 24, 0.079);
    finishFixture(nearColonization);
    nearColonization.reportedTrait[0] = -0.6;
    ArchipelagoField occupancy;
    p = ArchipelagoField::Parameters();
    p.selection = 0.0; p.mutation = 0.0; p.migration = 0.0;
    if (!occupancy.restore(nearColonization, p) ||
        std::fabs(occupancy.metrics().trait[0] + 0.6) > 1e-15)
        fail("empty habitat did not restore its held trait");
    occupancy.advance(0.02, p);
    if ((occupancy.metrics().occupiedMask & 1u) == 0u ||
        (occupancy.metrics().colonizeMask & 1u) == 0u ||
        std::fabs(occupancy.metrics().trait[0] - ArchipelagoField::traitAt(24)) > 1e-12)
        fail("colonization threshold, event, or held-trait release");
    occupancy.advance(0.0, p);
    if (occupancy.metrics().colonizeEvent || occupancy.metrics().extinctEvent)
        fail("occupancy event replayed on a zero-time update");

    ArchipelagoField::State nearExtinction = fixtureState();
    putMass(nearExtinction, 7, 31, 0.031);
    finishFixture(nearExtinction);
    nearExtinction.occupiedMask |= 1u << 7;
    ArchipelagoField extinction;
    p.selection = 8.0; p.climate = -0.85; p.gradient = -0.85;
    if (!extinction.restore(nearExtinction, p))
        fail("could not restore extinction fixture");
    const double heldBeforeLoss = extinction.metrics().trait[7];
    extinction.advance(0.02, p);
    if ((extinction.metrics().occupiedMask & (1u << 7)) != 0u ||
        (extinction.metrics().extinctMask & (1u << 7)) == 0u ||
        extinction.metrics().trait[7] != heldBeforeLoss)
        fail("extinction threshold, event, or held trait");

    // Identical calls are bit-exact. A full state restore resumes the same
    // future and emits no occupancy transition merely because it was loaded.
    ArchipelagoField replayA, replayB;
    for (int step = 0; step < 80; ++step) {
        p = ArchipelagoField::Parameters();
        p.selection = 8.0 * (step % 11) / 10.0;
        p.mutation = 0.003 * (step % 7) / 6.0;
        p.migration = 2.0 * (step % 9) / 8.0;
        p.gradient = -0.85 + 1.7 * (step % 13) / 12.0;
        p.climate = -0.85 + 1.7 * (step % 5) / 4.0;
        p.barrier = (step % 6) / 5.0;
        p.ring = (step & 1) != 0;
        const double delta = 0.003 + 0.017 * (step % 8);
        replayA.advance(delta, p);
        replayB.advance(delta, p);
        if (!sameState(replayA.state(), replayB.state()))
            fail("deterministic replay diverged");
    }
    ArchipelagoField resumed;
    if (!resumed.restore(replayA.state(), p) ||
        resumed.metrics().colonizeEvent || resumed.metrics().extinctEvent ||
        !sameState(resumed.state(), replayA.state()))
        fail("full state restore changed state or emitted an event");
    for (int step = 0; step < 40; ++step) {
        p.climate = 0.75 * std::sin(0.31 * step);
        p.gradient = 0.75 * std::cos(0.23 * step);
        p.ring = step % 3 == 0;
        replayA.advance(0.071, p);
        resumed.advance(0.071, p);
        if (!sameState(replayA.state(), resumed.state()))
            fail("restored deterministic future diverged");
    }

    const ArchipelagoField::State beforeBadRestore = resumed.state();
    ArchipelagoField::State bad = beforeBadRestore;
    bad.density[17] = std::numeric_limits<double>::quiet_NaN();
    if (resumed.restore(bad, p) || !sameState(beforeBadRestore, resumed.state()))
        fail("non-finite restore was not rejected atomically");
    bad = beforeBadRestore;
    for (int bin = 0; bin < ArchipelagoField::kBins; ++bin)
        bad.density[bin] = 0.0;
    bad.occupiedMask |= 1u;
    if (resumed.restore(bad, p) || !sameState(beforeBadRestore, resumed.state()))
        fail("inconsistent occupancy restore was accepted");
    bad = beforeBadRestore;
    bad.reportedGlobalMean = 2.0;
    if (resumed.restore(bad, p) || !sameState(beforeBadRestore, resumed.state()))
        fail("out-of-range held metric restore was accepted");

    // Outer-call granularity converges around the same positive split flow.
    ArchipelagoField coarse, fine;
    p = ArchipelagoField::Parameters();
    p.selection = 5.3; p.mutation = 0.0011; p.migration = 0.73;
    p.gradient = -0.63; p.climate = 0.24; p.barrier = 0.37; p.ring = true;
    evolve(coarse, p, 12.0, 0.137);
    evolve(fine, p, 12.0, 0.0065);
    const double timestepError = maxDensityDifference(coarse, fine);
    std::printf("timestep: max cell error %.8f (0.137 vs 0.0065 tau calls)\n",
                timestepError);
    if (timestepError > 0.0001 ||
        coarse.metrics().occupiedMask != fine.metrics().occupiedMask)
        fail("outer-call timestep convergence");

    // Non-finite and extreme control inputs fail into the documented domain.
    ArchipelagoField hostile;
    const ArchipelagoField::State beforeInfiniteTime = hostile.state();
    p.selection = std::numeric_limits<double>::quiet_NaN();
    p.mutation = std::numeric_limits<double>::infinity();
    p.migration = -std::numeric_limits<double>::infinity();
    p.gradient = 1e300; p.climate = -1e300; p.barrier = 1e300; p.ring = true;
    hostile.advance(std::numeric_limits<double>::infinity(), p);
    if (!sameState(beforeInfiniteTime, hostile.state()))
        fail("non-finite time altered model state");
    hostile.advance(1e300, p);
    hostile.advance(-1.0, p);
    if (!valid(hostile))
        fail("hostile inputs left the valid state space");

    // Exposed parameter corners stay positive, finite, bounded, and restorable.
    int cornerRuns = 0;
    for (double selection : {0.0, 8.0})
        for (double mutation : {0.0, 0.003})
            for (double migration : {0.0, 2.0})
                for (double gradient : {-0.85, 0.85})
                    for (double climateValue : {-0.85, 0.85})
                        for (double barrier : {0.0, 1.0})
                            for (int topology = 0; topology < 2; ++topology) {
                                ArchipelagoField corner;
                                p = ArchipelagoField::Parameters();
                                p.selection = selection;
                                p.mutation = mutation;
                                p.migration = migration;
                                p.gradient = gradient;
                                p.climate = climateValue;
                                p.barrier = barrier;
                                p.ring = topology != 0;
                                corner.reset(p);
                                for (int call = 0; call < 16; ++call) {
                                    corner.advance(0.128, p);
                                    if (!valid(corner))
                                        fail("parameter corner left valid state space");
                                }
                                ArchipelagoField copy;
                                if (!copy.restore(corner.state(), p) ||
                                    !sameState(copy.state(), corner.state()))
                                    fail("parameter-corner state was not restorable");
                                ++cornerRuns;
                            }
    std::printf("stability: %d full-range parameter corners bounded\n", cornerRuns);

    // Report the core-only cost at default RATE and at the maximum one-tick
    // request. Timing is diagnostic and hardware-dependent; state invariants
    // are the portable regression for the former denormal cliff.
    const int benchmarkCalls = 2000;
    ArchipelagoField defaultBenchmark;
    p = ArchipelagoField::Parameters();
    const std::chrono::steady_clock::time_point defaultStart =
        std::chrono::steady_clock::now();
    for (int call = 0; call < benchmarkCalls; ++call)
        defaultBenchmark.advance(0.008, p);
    const std::chrono::steady_clock::time_point defaultStop =
        std::chrono::steady_clock::now();
    const double defaultMicroseconds =
        std::chrono::duration<double, std::micro>(defaultStop - defaultStart).count() /
        benchmarkCalls;

    ArchipelagoField benchmark;
    p.selection = 8.0; p.mutation = 0.003; p.migration = 2.0;
    p.gradient = 0.85; p.climate = -0.85; p.barrier = 0.31; p.ring = true;
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    for (int call = 0; call < benchmarkCalls; ++call) {
        p.climate = (call & 1) ? -0.85 : 0.85;
        benchmark.advance(0.128, p);
    }
    const std::chrono::steady_clock::time_point stop = std::chrono::steady_clock::now();
    const double microseconds =
        std::chrono::duration<double, std::micro>(stop - start).count() /
        benchmarkCalls;
    std::printf("benchmark: %.2f us default-rate / %.2f us maximum field advance "
                "(%.2f%% / %.2f%% of one core at 500 Hz)\n",
                defaultMicroseconds, microseconds,
                defaultMicroseconds / 20.0, microseconds / 20.0);
    if (!valid(defaultBenchmark) || !valid(benchmark) ||
        defaultMicroseconds > 1000.0 || microseconds > 1000.0)
        fail("worst-case benchmark or terminal state");

    if (failures) {
        std::printf("archipelago stability: %d failure(s)\n", failures);
        return 1;
    }
    std::printf("archipelago stability: PASS\n");
    return 0;
}
