#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace coalescent {

// Eight habitats carrying a non-normalized abundance density over a fixed
// cell-centred trait axis. The continuous model is
//
//   du_i(x)/dt = u_i(x) * (1 - n_i - s * (x - theta_i)^2)
//                 + D * d2u_i/dx2
//                 + m * sum_j w_ij * (u_j(x) - u_i(x)),
//   theta_i = climate + gradient * habitatPosition(i).
//
// Trait diffusion is no-flux. Migration uses either a row or ring graph; the
// edge between habitats 3 and 4 is attenuated by the barrier control. Positive
// split operators keep every intermediate density non-negative without an
// audio-thread allocation.
class ArchipelagoField {
public:
    static const int kHabitats = 8;
    static const int kBins = 32;
    static const int kValues = kHabitats * kBins;

    static constexpr double traitMin() { return -1.0; }
    static constexpr double traitMax() { return 1.0; }
    static constexpr double binWidth() {
        return (traitMax() - traitMin()) / static_cast<double>(kBins);
    }

    static constexpr double selectionMin() { return 0.0; }
    static constexpr double selectionMax() { return 8.0; }
    static constexpr double mutationMin() { return 0.0; }
    static constexpr double mutationMax() { return 0.003; }
    static constexpr double migrationMin() { return 0.0; }
    static constexpr double migrationMax() { return 2.0; }
    static constexpr double gradientMin() { return -0.85; }
    static constexpr double gradientMax() { return 0.85; }
    static constexpr double climateMin() { return -0.85; }
    static constexpr double climateMax() { return 0.85; }
    static constexpr double barrierMin() { return 0.0; }
    static constexpr double barrierMax() { return 1.0; }
    static constexpr double colonizeThreshold() { return 0.08; }
    static constexpr double extinctThreshold() { return 0.03; }
    static constexpr double maxAdvance() { return 0.25; }
    static constexpr double maxSubstep() { return 0.02; }
    static constexpr double numericalExtinctionFloor() { return 1e-280; }

    struct Parameters {
        double selection;
        double mutation;
        double migration;
        double gradient;
        double climate;
        double barrier;
        bool ring;

        Parameters()
            : selection(2.88), mutation(1.30e-4), migration(0.0317),
              gradient(0.55), climate(0.0), barrier(0.0), ring(false) {}
    };

    struct Metrics {
        std::array<double, kHabitats> mass;
        std::array<double, kHabitats> trait;
        double totalMass;
        double globalMean;
        double difference;
        double flux;
        std::uint32_t occupiedMask;
        std::uint32_t colonizeMask;
        std::uint32_t extinctMask;
        bool colonizeEvent;
        bool extinctEvent;

        Metrics()
            : mass(), trait(), totalMass(0.0), globalMean(0.0),
              difference(0.0), flux(0.0), occupiedMask(0u),
              colonizeMask(0u), extinctMask(0u),
              colonizeEvent(false), extinctEvent(false) {
            mass.fill(0.0);
            trait.fill(0.0);
        }
    };

    // Density is habitat-major: density[habitat * kBins + traitBin]. Held
    // traits and occupancy are part of state because empty habitats retain a
    // meaningful output until recolonization.
    struct State {
        std::array<double, kValues> density;
        std::array<double, kHabitats> reportedTrait;
        double reportedGlobalMean;
        std::uint32_t occupiedMask;

        State()
            : density(), reportedTrait(), reportedGlobalMean(0.0),
              occupiedMask(0u) {
            density.fill(0.0);
            reportedTrait.fill(0.0);
        }
    };

    ArchipelagoField() { reset(Parameters()); }

    static double traitAt(int bin) {
        bin = std::max(0, std::min(bin, kBins - 1));
        return traitMin() + (static_cast<double>(bin) + 0.5) * binWidth();
    }

    static double habitatPosition(int habitat) {
        habitat = std::max(0, std::min(habitat, kHabitats - 1));
        return -1.0 + 2.0 * static_cast<double>(habitat) /
                          static_cast<double>(kHabitats - 1);
    }

    const Metrics& metrics() const { return metrics_; }
    const std::array<double, kValues>& densities() const { return density_; }

    double densityAt(int habitat, int bin) const {
        if (habitat < 0 || habitat >= kHabitats || bin < 0 || bin >= kBins)
            return 0.0;
        return density_[index(habitat, bin)];
    }

    void copyMasses(float* destination, std::size_t count) const {
        if (!destination || count < static_cast<std::size_t>(kValues))
            return;
        for (int i = 0; i < kValues; ++i)
            destination[i] = static_cast<float>(density_[i]);
    }

    State state() const {
        State saved;
        saved.density = density_;
        saved.reportedTrait = reportedTrait_;
        saved.reportedGlobalMean = reportedGlobalMean_;
        saved.occupiedMask = occupiedMask_;
        return saved;
    }

    void reset() { reset(Parameters()); }

    // Starts each habitat at 75% of carrying capacity, locally adapted to the
    // nearest represented trait. The environmental optimum itself is not
    // clipped during evolution, so an extreme climate can remain hostile.
    void reset(const Parameters& requested) {
        const Parameters p = sanitized(requested);
        const double sigma = 0.065;
        for (int habitat = 0; habitat < kHabitats; ++habitat) {
            const double optimum = clamp(p.climate + p.gradient * habitatPosition(habitat),
                                         traitMin(), traitMax());
            double sum = 0.0;
            for (int bin = 0; bin < kBins; ++bin) {
                const double z = (traitAt(bin) - optimum) / sigma;
                const double value = std::exp(-0.5 * z * z);
                density_[index(habitat, bin)] = value;
                sum += value;
            }
            const double scale = sum > massEpsilon() ? 0.75 / sum : 0.0;
            for (int bin = 0; bin < kBins; ++bin)
                density_[index(habitat, bin)] *= scale;
        }
        repairDensity();

        occupiedMask_ = allHabitatsMask();
        reportedTrait_.fill(0.0);
        reportedGlobalMean_ = 0.0;
        metrics_ = Metrics();
        refreshMetrics(p, false);
        clearEvents();
    }

    // Invalid payloads leave the previous state untouched. Restore derives
    // the public metrics but deliberately does not emit occupancy events.
    bool restore(const State& saved, const Parameters& requested) {
        if (!validState(saved))
            return false;

        const Parameters p = sanitized(requested);
        density_ = saved.density;
        reportedTrait_ = saved.reportedTrait;
        reportedGlobalMean_ = saved.reportedGlobalMean;
        occupiedMask_ = saved.occupiedMask;
        repairDensity();
        clearEvents();
        refreshMetrics(p, false);
        clearEvents();
        return true;
    }

    // Advances in model time. Large or non-finite requests are bounded before
    // reaching transcendental functions; changing Parameters alone never
    // mutates the density.
    void advance(double deltaTau, const Parameters& requested) {
        const Parameters p = sanitized(requested);
        clearEvents();

        if (!std::isfinite(deltaTau) || !(deltaTau > 0.0)) {
            refreshMetrics(p, false);
            clearEvents();
            return;
        }

        deltaTau = std::min(deltaTau, maxAdvance());
        const int steps = std::max(1, static_cast<int>(
            std::ceil(deltaTau / maxSubstep())));
        const double h = deltaTau / static_cast<double>(steps);
        const double diffusionR = 0.5 * h * p.mutation /
                                  (binWidth() * binWidth());
        double diffusionC[kBins];
        double diffusionInverse[kBins];
        const bool hasDiffusion = prepareDiffusion(
            diffusionR, diffusionC, diffusionInverse);
        double selectionFactor[kValues];
        const bool hasSelection = prepareSelection(
            0.5 * h, p, selectionFactor);
        double migrationMix[kHabitats];
        const bool hasMigration = prepareMigration(
            0.5 * h, p, migrationMix);
        for (int step = 0; step < steps; ++step) {
            if (hasMigration)
                migrate(p, migrationMix, false);
            if (hasDiffusion)
                diffuse(diffusionR, diffusionC, diffusionInverse);
            if (hasSelection)
                select(selectionFactor);
            logistic(h);
            if (hasSelection)
                select(selectionFactor);
            if (hasDiffusion)
                diffuse(diffusionR, diffusionC, diffusionInverse);
            if (hasMigration)
                migrate(p, migrationMix, true);
            repairDensity();
            refreshMetrics(p, true);
        }
    }

private:
    std::array<double, kValues> density_;
    std::array<double, kHabitats> reportedTrait_;
    double reportedGlobalMean_ = 0.0;
    std::uint32_t occupiedMask_ = 0u;
    Metrics metrics_;

    static constexpr double massEpsilon() { return 1e-12; }
    static constexpr double stateBinMax() { return 4.0; }
    static constexpr double stateHabitatMax() { return 4.0; }
    static constexpr std::uint32_t allHabitatsMask() {
        return (UINT32_C(1) << kHabitats) - UINT32_C(1);
    }

    static int index(int habitat, int bin) { return habitat * kBins + bin; }

    static double clamp(double value, double lo, double hi) {
        return std::max(lo, std::min(value, hi));
    }

    static double finiteOr(double value, double fallback) {
        return std::isfinite(value) ? value : fallback;
    }

    static Parameters sanitized(const Parameters& requested) {
        const Parameters defaults;
        Parameters p;
        p.selection = clamp(finiteOr(requested.selection, defaults.selection),
                            selectionMin(), selectionMax());
        p.mutation = clamp(finiteOr(requested.mutation, defaults.mutation),
                           mutationMin(), mutationMax());
        p.migration = clamp(finiteOr(requested.migration, defaults.migration),
                            migrationMin(), migrationMax());
        p.gradient = clamp(finiteOr(requested.gradient, defaults.gradient),
                           gradientMin(), gradientMax());
        p.climate = clamp(finiteOr(requested.climate, defaults.climate),
                          climateMin(), climateMax());
        p.barrier = clamp(finiteOr(requested.barrier, defaults.barrier),
                          barrierMin(), barrierMax());
        p.ring = requested.ring;
        return p;
    }

    static double edgePermeability(int edge, const Parameters& p) {
        if (edge != 3)
            return 1.0;
        const double open = 1.0 - p.barrier;
        return open * open;
    }

    static void edgeHabitats(int edge, int& left, int& right) {
        if (edge == kHabitats - 1) {
            left = kHabitats - 1;
            right = 0;
        }
        else {
            left = edge;
            right = edge + 1;
        }
    }

    static bool prepareMigration(double duration, const Parameters& p,
                                 double* mix) {
        for (int edge = 0; edge < kHabitats; ++edge)
            mix[edge] = 0.0;
        if (!(duration > 0.0) || !(p.migration > 0.0))
            return false;
        const int edgeCount = p.ring ? kHabitats : kHabitats - 1;
        for (int edge = 0; edge < edgeCount; ++edge) {
            const double rate = p.migration * edgePermeability(edge, p);
            mix[edge] = 0.5 * (1.0 - std::exp(-2.0 * rate * duration));
        }
        return true;
    }

    void migrate(const Parameters& p, const double* mix, bool reverse) {
        const int edgeCount = p.ring ? kHabitats : kHabitats - 1;
        for (int order = 0; order < edgeCount; ++order) {
            const int edge = reverse ? edgeCount - 1 - order : order;
            int left = 0;
            int right = 0;
            edgeHabitats(edge, left, right);
            const double amount = mix[edge];
            if (!(amount > 0.0))
                continue;
            for (int bin = 0; bin < kBins; ++bin) {
                const int aIndex = index(left, bin);
                const int bIndex = index(right, bin);
                const double a = density_[aIndex];
                const double b = density_[bIndex];
                density_[aIndex] = a + amount * (b - a);
                density_[bIndex] = b + amount * (a - b);
            }
        }
    }

    // Backward Euler solve for a cell-centred finite-volume Laplacian with
    // Neumann boundary rows. The matrix is positive and has unit row sums.
    static bool prepareDiffusion(double r, double* cPrime,
                                 double* inverseDenominator) {
        if (!(r > 0.0) || !std::isfinite(r))
            return false;

        double denominator = 1.0 + r;
        inverseDenominator[0] = 1.0 / denominator;
        cPrime[0] = -r * inverseDenominator[0];
        for (int bin = 1; bin < kBins; ++bin) {
            const double upper = bin + 1 < kBins ? -r : 0.0;
            const double diagonal = bin + 1 < kBins ? 1.0 + 2.0 * r : 1.0 + r;
            denominator = diagonal + r * cPrime[bin - 1];
            inverseDenominator[bin] = 1.0 / denominator;
            cPrime[bin] = upper * inverseDenominator[bin];
        }
        return true;
    }

    void diffuse(double r, const double* cPrime,
                 const double* inverseDenominator) {

        double dPrime[kBins];
        double solution[kBins];
        for (int habitat = 0; habitat < kHabitats; ++habitat) {
            dPrime[0] = density_[index(habitat, 0)] * inverseDenominator[0];
            for (int bin = 1; bin < kBins; ++bin) {
                dPrime[bin] = (density_[index(habitat, bin)]
                               + r * dPrime[bin - 1]) * inverseDenominator[bin];
            }
            solution[kBins - 1] = dPrime[kBins - 1];
            for (int bin = kBins - 2; bin >= 0; --bin)
                solution[bin] = dPrime[bin] - cPrime[bin] * solution[bin + 1];
            for (int bin = 0; bin < kBins; ++bin)
                density_[index(habitat, bin)] = solution[bin];
        }
    }

    static bool prepareSelection(double duration, const Parameters& p,
                                 double* factor) {
        if (!(duration > 0.0) || !(p.selection > 0.0))
            return false;
        for (int habitat = 0; habitat < kHabitats; ++habitat) {
            const double optimum = p.climate + p.gradient * habitatPosition(habitat);
            for (int bin = 0; bin < kBins; ++bin) {
                const double distance = traitAt(bin) - optimum;
                factor[index(habitat, bin)] =
                    std::exp(-duration * p.selection * distance * distance);
            }
        }
        return true;
    }

    void select(const double* factor) {
        for (int i = 0; i < kValues; ++i)
            density_[i] *= factor[i];
    }

    // Exact scalar logistic flow, applied as a common scale to preserve each
    // habitat's trait composition during density regulation.
    void logistic(double duration) {
        if (!(duration > 0.0))
            return;
        const double growth = std::exp(duration);
        for (int habitat = 0; habitat < kHabitats; ++habitat) {
            double abundance = 0.0;
            for (int bin = 0; bin < kBins; ++bin)
                abundance += density_[index(habitat, bin)];
            if (!(abundance > 0.0))
                continue;
            const double denominator = 1.0 + abundance * (growth - 1.0);
            const double scale = growth / denominator;
            for (int bin = 0; bin < kBins; ++bin)
                density_[index(habitat, bin)] *= scale;
        }
    }

    void repairDensity() {
        for (int i = 0; i < kValues; ++i) {
            // This is far below seeded standing variation and observable mass,
            // but leaves ample headroom above slow double-subnormal arithmetic.
            if (!std::isfinite(density_[i]) ||
                density_[i] < numericalExtinctionFloor())
                density_[i] = 0.0;
            else if (density_[i] > stateBinMax())
                density_[i] = stateBinMax();
        }
    }

    void clearEvents() {
        metrics_.colonizeMask = 0u;
        metrics_.extinctMask = 0u;
        metrics_.colonizeEvent = false;
        metrics_.extinctEvent = false;
    }

    void refreshMetrics(const Parameters& p, bool emitEvents) {
        const std::uint32_t priorColonize = metrics_.colonizeMask;
        const std::uint32_t priorExtinct = metrics_.extinctMask;
        std::array<double, kHabitats> rawTrait;
        rawTrait.fill(0.0);
        metrics_.totalMass = 0.0;
        double totalMoment = 0.0;

        for (int habitat = 0; habitat < kHabitats; ++habitat) {
            double abundance = 0.0;
            double moment = 0.0;
            for (int bin = 0; bin < kBins; ++bin) {
                const double value = density_[index(habitat, bin)];
                abundance += value;
                moment += value * traitAt(bin);
            }
            metrics_.mass[habitat] = abundance;
            rawTrait[habitat] = abundance > massEpsilon()
                ? clamp(moment / abundance, traitMin(), traitMax())
                : reportedTrait_[habitat];
            metrics_.totalMass += abundance;
            totalMoment += moment;

            const std::uint32_t bit = UINT32_C(1) << habitat;
            const bool wasOccupied = (occupiedMask_ & bit) != 0u;
            bool occupied = wasOccupied;
            if (wasOccupied && abundance <= extinctThreshold())
                occupied = false;
            else if (!wasOccupied && abundance >= colonizeThreshold())
                occupied = true;

            if (occupied)
                occupiedMask_ |= bit;
            else
                occupiedMask_ &= ~bit;
            if (emitEvents && occupied != wasOccupied) {
                if (occupied)
                    metrics_.colonizeMask |= bit;
                else
                    metrics_.extinctMask |= bit;
            }
            if (occupied)
                reportedTrait_[habitat] = rawTrait[habitat];
            metrics_.trait[habitat] = reportedTrait_[habitat];
        }

        if (metrics_.totalMass > massEpsilon())
            reportedGlobalMean_ = clamp(totalMoment / metrics_.totalMass,
                                        traitMin(), traitMax());
        metrics_.globalMean = reportedGlobalMean_;

        double leftMass = 0.0;
        double rightMass = 0.0;
        double leftMoment = 0.0;
        double rightMoment = 0.0;
        for (int habitat = 0; habitat < kHabitats; ++habitat) {
            const double moment = metrics_.mass[habitat] * rawTrait[habitat];
            if (habitat < kHabitats / 2) {
                leftMass += metrics_.mass[habitat];
                leftMoment += moment;
            }
            else {
                rightMass += metrics_.mass[habitat];
                rightMoment += moment;
            }
        }
        metrics_.difference = leftMass > massEpsilon() && rightMass > massEpsilon()
            ? clamp(rightMoment / rightMass - leftMoment / leftMass,
                    traitMin() - traitMax(), traitMax() - traitMin())
            : 0.0;

        const int edgeCount = p.ring ? kHabitats : kHabitats - 1;
        double fluxNumerator = 0.0;
        for (int edge = 0; edge < edgeCount; ++edge) {
            int left = 0;
            int right = 0;
            edgeHabitats(edge, left, right);
            double difference = 0.0;
            for (int bin = 0; bin < kBins; ++bin)
                difference += std::fabs(density_[index(left, bin)]
                                        - density_[index(right, bin)]);
            fluxNumerator += edgePermeability(edge, p) * 0.5 * difference;
        }
        const double fluxScale = migrationMax() * static_cast<double>(edgeCount);
        metrics_.flux = fluxScale > 0.0
            ? clamp(p.migration * fluxNumerator / fluxScale, 0.0, 1.0)
            : 0.0;

        metrics_.occupiedMask = occupiedMask_;
        metrics_.colonizeMask |= priorColonize;
        metrics_.extinctMask |= priorExtinct;
        metrics_.colonizeEvent = metrics_.colonizeMask != 0u;
        metrics_.extinctEvent = metrics_.extinctMask != 0u;
    }

    static bool validState(const State& saved) {
        if (!std::isfinite(saved.reportedGlobalMean) ||
            saved.reportedGlobalMean < traitMin() ||
            saved.reportedGlobalMean > traitMax() ||
            (saved.occupiedMask & ~allHabitatsMask()) != 0u)
            return false;

        for (int habitat = 0; habitat < kHabitats; ++habitat) {
            if (!std::isfinite(saved.reportedTrait[habitat]) ||
                saved.reportedTrait[habitat] < traitMin() ||
                saved.reportedTrait[habitat] > traitMax())
                return false;
            double abundance = 0.0;
            for (int bin = 0; bin < kBins; ++bin) {
                const double value = saved.density[index(habitat, bin)];
                if (!std::isfinite(value) || value < 0.0 || value > stateBinMax())
                    return false;
                abundance += value;
            }
            if (!std::isfinite(abundance) || abundance > stateHabitatMax())
                return false;
            const bool occupied =
                (saved.occupiedMask & (UINT32_C(1) << habitat)) != 0u;
            if ((occupied && abundance <= extinctThreshold()) ||
                (!occupied && abundance >= colonizeThreshold()))
                return false;
        }
        return true;
    }
};

} // namespace coalescent
