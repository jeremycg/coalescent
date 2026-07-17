#pragma once

#include "pcg32.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace coalescent {

// Four-deme, two-allele Wright-Fisher model. Each generation applies genic
// selection, symmetric mutation, common-pool migration, then independent
// binomial sampling. Counts and denominators remain the authoritative state;
// changing Parameters alone therefore never redraws a population.
class IslandsModel {
public:
    static const int kIslands = 4;

    static constexpr std::uint32_t copiesMin() { return 8u; }
    static constexpr std::uint32_t copiesMax() { return 4096u; }
    static constexpr double selectionMin() { return -0.25; }
    static constexpr double selectionMax() { return 0.25; }
    static constexpr double mutationMin() { return 0.0; }
    static constexpr double mutationMax() { return 0.05; }
    static constexpr double migrationMin() { return 0.0; }
    static constexpr double migrationMax() { return 1.0; }
    static constexpr std::uint32_t founderCopies() { return 8u; }

    // Public alias retained for existing SDK-free callers and golden tests.
    using Pcg32 = coalescent::Pcg32;

    struct Parameters {
        std::uint32_t copies;
        double selection;
        double mutation;
        double migration;

        Parameters()
            : copies(64u), selection(0.0), mutation(0.0), migration(0.0) {}
    };

    struct Metrics {
        std::array<double, kIslands> frequency;
        double mean;
        double diversity;
        bool fixA;
        bool fixB;
        bool lossEvent;
        bool sweepEvent;

        Metrics()
            : frequency{{0.5, 0.5, 0.5, 0.5}}, mean(0.5), diversity(1.0),
              fixA(false), fixB(false), lossEvent(false), sweepEvent(false) {}
    };

    // All fields are fixed-width or booleans and can be serialized directly.
    // Encode the uint64 RNG fields losslessly (for example as hex strings), not
    // through a JSON floating-point number.
    struct State {
        std::array<std::uint32_t, kIslands> counts;
        std::array<std::uint32_t, kIslands> denominators;
        std::uint64_t rngState;
        std::uint64_t rngIncrement;
        std::uint32_t nextFounder;
        bool wasFixA;
        bool wasFixB;
        bool sweepAArmed;
        bool sweepBArmed;

        State()
            : counts{{32u, 32u, 32u, 32u}},
              denominators{{64u, 64u, 64u, 64u}},
              rngState(0u), rngIncrement(1u), nextFounder(0u),
              wasFixA(false), wasFixB(false),
              sweepAArmed(true), sweepBArmed(true) {}
    };

    IslandsModel(std::uint64_t seedValue = 42u, std::uint64_t stream = 54u) {
        reset(seedValue, stream);
    }

    void reset(std::uint64_t seedValue = 42u, std::uint64_t stream = 54u) {
        counts_.fill(32u);
        denominators_.fill(64u);
        rng_.seed(seedValue, stream);
        nextFounder_ = 0u;
        wasFixA_ = false;
        wasFixB_ = false;
        sweepAArmed_ = true;
        sweepBArmed_ = true;
        refreshMetrics(false);
    }

    const Metrics& metrics() const { return metrics_; }

    State state() const {
        State saved;
        saved.counts = counts_;
        saved.denominators = denominators_;
        const Pcg32::State rngState = rng_.state();
        saved.rngState = rngState.state;
        saved.rngIncrement = rngState.increment;
        saved.nextFounder = nextFounder_;
        saved.wasFixA = wasFixA_;
        saved.wasFixB = wasFixB_;
        saved.sweepAArmed = sweepAArmed_;
        saved.sweepBArmed = sweepBArmed_;
        return saved;
    }

    // Invalid payloads leave the previous state untouched. A valid restore
    // reconstructs metrics but deliberately emits no LOSS or SWEEP event.
    bool restore(const State& saved) {
        if (!validState(saved))
            return false;

        Pcg32 restoredRng;
        if (!restoredRng.restore(Pcg32::State(saved.rngState,
                                              saved.rngIncrement)))
            return false;

        counts_ = saved.counts;
        denominators_ = saved.denominators;
        rng_ = restoredRng;
        nextFounder_ = saved.nextFounder;
        wasFixA_ = saved.wasFixA;
        wasFixB_ = saved.wasFixB;
        sweepAArmed_ = saved.sweepAArmed;
        sweepBArmed_ = saved.sweepBArmed;
        refreshMetrics(false);
        return true;
    }

    void advance(const Parameters& requested) {
        metrics_.lossEvent = false;
        metrics_.sweepEvent = false;
        const Parameters p = sanitized(requested);

        std::array<double, kIslands> pool;
        double poolMean = 0.0;
        const double relativeFitness = std::exp(p.selection);
        for (int i = 0; i < kIslands; ++i) {
            const double frequency = static_cast<double>(counts_[i]) /
                                     static_cast<double>(denominators_[i]);
            const double denominator = 1.0 - frequency + frequency * relativeFitness;
            const double selected = denominator > 0.0
                ? frequency * relativeFitness / denominator : frequency;
            pool[i] = p.mutation + (1.0 - 2.0 * p.mutation) * selected;
            poolMean += pool[i];
        }
        poolMean /= static_cast<double>(kIslands);

        for (int i = 0; i < kIslands; ++i) {
            const double probability = clamp01((1.0 - p.migration) * pool[i]
                                               + p.migration * poolMean);
            counts_[i] = binomial(p.copies, probability);
            denominators_[i] = p.copies;
        }
        refreshMetrics(true);
    }

    // A founder event is a bottleneck, not a generation: the rotating target
    // island is resampled at eight copies from its own pre-trigger frequency.
    // Other islands and the generation parameters are untouched.
    void founder() {
        metrics_.lossEvent = false;
        metrics_.sweepEvent = false;
        const std::uint32_t island = nextFounder_;
        const double frequency = static_cast<double>(counts_[island]) /
                                 static_cast<double>(denominators_[island]);
        counts_[island] = binomial(founderCopies(), frequency);
        denominators_[island] = founderCopies();
        nextFounder_ = (nextFounder_ + 1u) % static_cast<std::uint32_t>(kIslands);
        refreshMetrics(true);
    }

private:
    std::array<std::uint32_t, kIslands> counts_;
    std::array<std::uint32_t, kIslands> denominators_;
    Pcg32 rng_;
    std::uint32_t nextFounder_ = 0u;
    bool wasFixA_ = false;
    bool wasFixB_ = false;
    bool sweepAArmed_ = true;
    bool sweepBArmed_ = true;
    Metrics metrics_;

    static double finiteOr(double value, double fallback) {
        return std::isfinite(value) ? value : fallback;
    }

    static double clamp(double value, double low, double high) {
        return std::max(low, std::min(value, high));
    }

    static double clamp01(double value) {
        return clamp(value, 0.0, 1.0);
    }

    static Parameters sanitized(const Parameters& requested) {
        const Parameters defaults;
        Parameters p;
        p.copies = std::max(copiesMin(), std::min(requested.copies, copiesMax()));
        p.selection = clamp(finiteOr(requested.selection, defaults.selection),
                            selectionMin(), selectionMax());
        p.mutation = clamp(finiteOr(requested.mutation, defaults.mutation),
                           mutationMin(), mutationMax());
        p.migration = clamp(finiteOr(requested.migration, defaults.migration),
                            migrationMin(), migrationMax());
        return p;
    }

    // Exact Bernoulli counting, rather than a Gaussian/Poisson approximation or
    // std::binomial_distribution. The 32-bit threshold makes the bit-for-bit
    // outcome part of this model's deterministic contract.
    std::uint32_t binomial(std::uint32_t trials, double probability) {
        if (!(probability > 0.0))
            return 0u;
        if (probability >= 1.0)
            return trials;

        const double scale = 4294967296.0; // 2^32
        const std::uint64_t threshold = static_cast<std::uint64_t>(probability * scale);
        std::uint32_t successes = 0u;
        for (std::uint32_t trial = 0u; trial < trials; ++trial)
            if (static_cast<std::uint64_t>(rng_.next()) < threshold)
                ++successes;
        return successes;
    }

    static bool validState(const State& saved) {
        if (saved.nextFounder >= static_cast<std::uint32_t>(kIslands) ||
            (saved.rngIncrement & 1u) == 0u)
            return false;

        bool fixA = true;
        bool fixB = true;
        double mean = 0.0;
        for (int i = 0; i < kIslands; ++i) {
            if (saved.denominators[i] < copiesMin() ||
                saved.denominators[i] > copiesMax() ||
                saved.counts[i] > saved.denominators[i])
                return false;
            fixA = fixA && saved.counts[i] == saved.denominators[i];
            fixB = fixB && saved.counts[i] == 0u;
            mean += static_cast<double>(saved.counts[i]) /
                    static_cast<double>(saved.denominators[i]);
        }
        mean /= static_cast<double>(kIslands);

        // These invariants are maintained by refreshMetrics(). Rejecting
        // inconsistent latch payloads prevents a corrupted patch manufacturing
        // an event on its first subsequent generation.
        if (saved.wasFixA != fixA || saved.wasFixB != fixB)
            return false;
        if ((mean >= 0.9 && saved.sweepAArmed) ||
            (mean <= 0.8 && !saved.sweepAArmed) ||
            (mean <= 0.1 && saved.sweepBArmed) ||
            (mean >= 0.2 && !saved.sweepBArmed))
            return false;
        return true;
    }

    void refreshMetrics(bool evaluateEvents) {
        metrics_.lossEvent = false;
        metrics_.sweepEvent = false;
        metrics_.mean = 0.0;
        double diversity = 0.0;
        bool fixA = true;
        bool fixB = true;
        for (int i = 0; i < kIslands; ++i) {
            const double frequency = static_cast<double>(counts_[i]) /
                                     static_cast<double>(denominators_[i]);
            metrics_.frequency[i] = frequency;
            metrics_.mean += frequency;
            diversity += frequency * (1.0 - frequency);
            fixA = fixA && counts_[i] == denominators_[i];
            fixB = fixB && counts_[i] == 0u;
        }
        metrics_.mean /= static_cast<double>(kIslands);
        metrics_.diversity = clamp01(diversity); // 4 * mean(p * (1-p))
        metrics_.fixA = fixA;
        metrics_.fixB = fixB;

        if (!evaluateEvents)
            return;

        metrics_.lossEvent = (fixA && !wasFixA_) || (fixB && !wasFixB_);
        wasFixA_ = fixA;
        wasFixB_ = fixB;

        if (metrics_.mean <= 0.8)
            sweepAArmed_ = true;
        if (metrics_.mean >= 0.2)
            sweepBArmed_ = true;
        if (sweepAArmed_ && metrics_.mean >= 0.9) {
            sweepAArmed_ = false;
            metrics_.sweepEvent = true;
        }
        if (sweepBArmed_ && metrics_.mean <= 0.1) {
            sweepBArmed_ = false;
            metrics_.sweepEvent = true;
        }
    }
};

} // namespace coalescent
