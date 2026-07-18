#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace coalescent {

// Fixed trait-density model for Finches. The state is a probability mass over
// cell-centred trait bins on [-1, 1]. Selection is a mass-conserving replicator
// step and mutation is no-flux diffusion:
//
//   dp/dt = p * (fitness - meanFitness) + mutation * d2p/dx2
//   fitness(x) = -0.5 * ((x - environment) / niche)^2
//                - integral competition(x-y) p(y) dy
//
// competition is Gaussian. Its width is niche / sqrt(branching), so the
// fitness curvature at a monomorphic environmental optimum changes sign at the
// dimensionless branching number B = 1. The split update uses a positive
// exponential replicator step followed by backward-Euler diffusion; both are
// allocation-free and the result is normalized after every substep.
class FinchesField {
public:
    static const int kBins = 64;
    // C++11 constexpr functions avoid the out-of-class definitions (and ODR
    // hazards) that header-only static float data members would require.
    static constexpr float traitMin() { return -1.f; }
    static constexpr float traitMax() { return 1.f; }
    static constexpr float binWidth() { return (traitMax() - traitMin()) / kBins; }

    // Exposed, tested numerical domain. Callers may pass values outside it;
    // advance() sanitizes and clamps before any indexing or transcendental use.
    static constexpr float mutationMin() { return 0.00002f; }
    static constexpr float mutationMax() { return 0.004f; }
    static constexpr float branchingMin() { return 0.35f; }
    static constexpr float branchingMax() { return 2.8f; }
    static constexpr float nicheMin() { return 0.20f; }
    static constexpr float nicheMax() { return 0.45f; }
    static constexpr float environmentMin() { return -0.55f; }
    static constexpr float environmentMax() { return 0.55f; }
    static constexpr float maxAdvance() { return 0.32f; }
    static constexpr float maxSubstep() { return 0.02f; }
    static constexpr float numericalExtinctionFloor() { return 1e-30f; }
    static constexpr float splitPersistence() { return 0.30f; }
    static constexpr float mergePersistence() { return 0.20f; }
    static constexpr int stateVersion() { return 1; }

    struct Parameters {
        float mutation;
        float branching;
        float niche;
        float environment;

        Parameters()
            : mutation(0.00004f), branching(2.1f), niche(0.32f), environment(0.f) {}
    };

    struct Metrics {
        float mean;
        float spread;
        float lowTrait;
        float highTrait;
        float lowMass;
        float highMass;
        float divergence;
        bool split;
        bool splitEvent;
        bool mergeEvent;

        Metrics()
            : mean(0.f), spread(0.f), lowTrait(0.f),
              highTrait(0.f), lowMass(1.f), highMass(0.f), divergence(0.f),
              split(false), splitEvent(false), mergeEvent(false) {}
    };

    // Complete authoritative state. Detector timers are evolutionary time, not
    // wall-clock time. Events themselves are intentionally one-call signals and
    // are never persisted.
    struct State {
        int version;
        std::array<float, kBins> mass;
        bool split;
        float splitTimer;
        float mergeTimer;

        State()
            : version(stateVersion()), split(false), splitTimer(0.f),
              mergeTimer(0.f) {
            mass.fill(0.f);
        }
    };

    FinchesField() { reset(); }

    static float traitAt(int index) {
        return traitMin() + (static_cast<float>(index) + 0.5f) * binWidth();
    }

    const std::array<float, kBins>& masses() const { return mass_; }
    const Metrics& metrics() const { return metrics_; }

    State state() const {
        State result;
        result.mass = mass_;
        result.split = split_;
        result.splitTimer = splitTimer_;
        result.mergeTimer = mergeTimer_;
        return result;
    }

    void copyMasses(float* destination, std::size_t count) const {
        if (!destination || count < static_cast<std::size_t>(kBins))
            return;
        std::copy(mass_.begin(), mass_.end(), destination);
    }

    // Deterministic narrow ancestor at the requested environmental optimum.
    void reset(float environment = 0.f) {
        environment = finiteOr(environment, 0.f);
        environment = clamp(environment, environmentMin(), environmentMax());
        const float sigma = 0.055f;
        for (int i = 0; i < kBins; ++i) {
            float z = (traitAt(i) - environment) / sigma;
            mass_[i] = std::exp(-0.5f * z * z);
        }
        normalize();
        split_ = false;
        splitTimer_ = 0.f;
        mergeTimer_ = 0.f;
        refreshMetrics(0.f, true);
    }

    // Mixes a narrow mutant cohort into the existing population. Amount is a
    // fraction of total post-injection mass, rather than an unbounded addition.
    void seed(float trait, float amount = 0.02f, float width = 0.035f) {
        metrics_.splitEvent = false;
        metrics_.mergeEvent = false;
        trait = clamp(finiteOr(trait, 0.f), traitMin(), traitMax());
        amount = clamp(finiteOr(amount, 0.02f), 0.f, 0.25f);
        width = clamp(finiteOr(width, 0.035f), 0.5f * binWidth(), 0.20f);
        if (!(amount > 0.f))
            return;

        float mutant[kBins];
        float mutantSum = 0.f;
        for (int i = 0; i < kBins; ++i) {
            float z = (traitAt(i) - trait) / width;
            mutant[i] = std::exp(-0.5f * z * z);
            mutantSum += mutant[i];
        }
        if (!(mutantSum > 0.f) || !std::isfinite(mutantSum))
            return;
        float inv = 1.f / mutantSum;
        for (int i = 0; i < kBins; ++i)
            mass_[i] = (1.f - amount) * mass_[i] + amount * mutant[i] * inv;
        normalize();
        refreshMetrics(0.f, false);
    }

    // Restores complete authored state without manufacturing an event. A bad
    // payload leaves the previous state untouched.
    bool restore(const State& source) {
        if (source.version != stateVersion()
            || !std::isfinite(source.splitTimer)
            || !std::isfinite(source.mergeTimer)
            || source.splitTimer < 0.f || source.splitTimer >= splitPersistence()
            || source.mergeTimer < 0.f || source.mergeTimer >= mergePersistence()
            || (source.split && source.splitTimer != 0.f)
            || (!source.split && source.mergeTimer != 0.f))
            return false;

        double sum = 0.0;
        for (int i = 0; i < kBins; ++i) {
            const float value = source.mass[i];
            if (!std::isfinite(value) || value < 0.f
                || (value > 0.f && value < numericalExtinctionFloor()))
                return false;
            sum += value;
        }
        if (!std::isfinite(sum) || std::fabs(sum - 1.0) > 1e-4)
            return false;

        mass_ = source.mass;
        split_ = source.split;
        splitTimer_ = source.splitTimer;
        mergeTimer_ = source.mergeTimer;
        metrics_.splitEvent = false;
        metrics_.mergeEvent = false;
        refreshMetrics(0.f, false, false);
        return true;
    }

    // Compatibility path for density-only state written before State v1. The
    // detector is initialized from the strong threshold because no historical
    // latch or persistence timers were stored by that format.
    bool restoreMasses(const float* source, std::size_t count) {
        if (!source || count != static_cast<std::size_t>(kBins))
            return false;
        std::array<float, kBins> restored;
        double sum = 0.0;
        for (int i = 0; i < kBins; ++i) {
            if (!std::isfinite(source[i]) || source[i] < 0.f)
                return false;
            restored[i] = source[i] < numericalExtinctionFloor() ? 0.f : source[i];
            sum += restored[i];
        }
        if (!(sum > 1e-12) || !std::isfinite(sum))
            return false;
        mass_ = restored;
        normalize();
        splitTimer_ = 0.f;
        mergeTimer_ = 0.f;
        refreshMetrics(0.f, true);
        return true;
    }

    // Advances by evolutionary time, not wall-clock seconds. Work is bounded:
    // excessively large deltas are clamped to maxAdvance(); non-finite deltas
    // fail closed and leave the state unchanged.
    void advance(float deltaTau, const Parameters& requested) {
        metrics_.splitEvent = false;
        metrics_.mergeEvent = false;

        deltaTau = finiteOr(deltaTau, 0.f);
        if (!(deltaTau > 0.f)) {
            refreshMetrics(0.f, false);
            return;
        }
        deltaTau = std::min(deltaTau, maxAdvance());
        Parameters p = sanitized(requested);
        buildCompetitionKernel(p.niche, p.branching);

        int steps = std::max(1, static_cast<int>(std::ceil(deltaTau / maxSubstep())));
        float h = deltaTau / static_cast<float>(steps);
        // Detector persistence follows the same bounded cadence as the field;
        // latch any crossing so later substeps cannot hide it from the caller.
        bool splitEvent = false;
        bool mergeEvent = false;
        for (int step = 0; step < steps; ++step) {
            reaction(h, p);
            diffuse(h * p.mutation / (binWidth() * binWidth()));
            normalize();
            refreshMetrics(h, false);
            splitEvent = splitEvent || metrics_.splitEvent;
            mergeEvent = mergeEvent || metrics_.mergeEvent;
        }
        metrics_.splitEvent = splitEvent;
        metrics_.mergeEvent = mergeEvent;
    }

private:
    struct Peaks {
        int lowIndex;
        int highIndex;
        int valleyIndex;
        float lowTrait;
        float highTrait;
        float lowHeight;
        float highHeight;
        float valleyHeight;
        float lowMass;
        float highMass;
        bool pair;

        Peaks()
            : lowIndex(0), highIndex(0), valleyIndex(0), lowTrait(0.f),
              highTrait(0.f), lowHeight(0.f), highHeight(0.f), valleyHeight(0.f),
              lowMass(1.f), highMass(0.f), pair(false) {}
    };

    std::array<float, kBins> mass_;
    std::array<float, kBins> kernel_;
    Metrics metrics_;
    bool split_ = false;
    float splitTimer_ = 0.f;
    float mergeTimer_ = 0.f;
    float kernelSigma_ = -1.f;

    static float clamp(float x, float lo, float hi) {
        return std::max(lo, std::min(x, hi));
    }

    static float finiteOr(float x, float fallback) {
        return std::isfinite(x) ? x : fallback;
    }

    static Parameters sanitized(const Parameters& requested) {
        Parameters defaults;
        Parameters p;
        p.mutation = clamp(finiteOr(requested.mutation, defaults.mutation),
                           mutationMin(), mutationMax());
        p.branching = clamp(finiteOr(requested.branching, defaults.branching),
                            branchingMin(), branchingMax());
        p.niche = clamp(finiteOr(requested.niche, defaults.niche),
                        nicheMin(), nicheMax());
        p.environment = clamp(finiteOr(requested.environment, defaults.environment),
                              environmentMin(), environmentMax());
        return p;
    }

    void buildCompetitionKernel(float niche, float branching) {
        float sigma = niche / std::sqrt(branching);
        if (sigma == kernelSigma_)
            return;
        kernelSigma_ = sigma;
        float invTwoSigma2 = 0.5f / (sigma * sigma);
        for (int d = 0; d < kBins; ++d) {
            float distance = static_cast<float>(d) * binWidth();
            kernel_[d] = std::exp(-distance * distance * invTwoSigma2);
        }
    }

    void reaction(float h, const Parameters& p) {
        float fitness[kBins];
        double meanFitness = 0.0;
        float invNiche = 1.f / p.niche;
        for (int i = 0; i < kBins; ++i) {
            // Accumulate equal-distance neighbours as a pair. Besides reducing
            // kernel lookups, this preserves exact mirror symmetry for a
            // symmetric population instead of selecting a side by summation
            // order at the unstable branching point.
            double competition = mass_[i];
            int maxDistance = std::max(i, kBins - 1 - i);
            for (int d = 1; d <= maxDistance; ++d) {
                float neighbours = 0.f;
                if (i - d >= 0) neighbours += mass_[i - d];
                if (i + d < kBins) neighbours += mass_[i + d];
                competition += static_cast<double>(neighbours) * kernel_[d];
            }
            float z = (traitAt(i) - p.environment) * invNiche;
            fitness[i] = -0.5f * z * z - static_cast<float>(competition);
            meanFitness += static_cast<double>(mass_[i]) * fitness[i];
        }

        for (int i = 0; i < kBins; ++i) {
            float exponent = h * (fitness[i] - static_cast<float>(meanFitness));
            exponent = clamp(exponent, -50.f, 50.f);
            mass_[i] *= std::exp(exponent);
        }
        normalize();
    }

    // Thomas solve for (I - r*L) pNew = p with finite-volume Neumann rows:
    // boundary diagonal 1+r, interior diagonal 1+2r. This symmetric matrix has
    // unit row sums, so its inverse conserves total mass up to roundoff.
    void diffuse(float r) {
        if (!(r > 0.f) || !std::isfinite(r))
            return;
        double fromLeft[kBins];
        double fromRight[kBins];
        solveDiffusion(static_cast<double>(r), false, fromLeft);
        solveDiffusion(static_cast<double>(r), true, fromRight);
        for (int i = 0; i < kBins; ++i)
            mass_[i] = static_cast<float>(0.5 * (fromLeft[i] + fromRight[i]));
    }

    // Solve once in either coordinate direction. Averaging the two equivalent
    // solves above removes the tiny directional bias of Thomas elimination;
    // that bias otherwise becomes audible after a perfectly symmetric ancestor
    // reaches a symmetry-breaking bifurcation.
    void solveDiffusion(double r, bool reverse, double* output) const {
        double cPrime[kBins];
        double dPrime[kBins];
        double solution[kBins];

        double diagonal = 1.0 + r;
        cPrime[0] = -r / diagonal;
        dPrime[0] = mass_[reverse ? kBins - 1 : 0] / diagonal;
        for (int i = 1; i < kBins; ++i) {
            double lower = -r;
            double upper = (i + 1 < kBins) ? -r : 0.0;
            diagonal = (i + 1 < kBins) ? (1.0 + 2.0 * r) : (1.0 + r);
            double denominator = diagonal - lower * cPrime[i - 1];
            cPrime[i] = upper / denominator;
            int source = reverse ? kBins - 1 - i : i;
            dPrime[i] = (mass_[source] - lower * dPrime[i - 1]) / denominator;
        }

        solution[kBins - 1] = dPrime[kBins - 1];
        for (int i = kBins - 2; i >= 0; --i)
            solution[i] = dPrime[i] - cPrime[i] * solution[i + 1];
        for (int i = 0; i < kBins; ++i)
            output[reverse ? kBins - 1 - i : i] = solution[i];
    }

    void normalize() {
        double sum = 0.0;
        for (int i = 0; i < kBins; ++i) {
            if (!std::isfinite(mass_[i])
                || mass_[i] < numericalExtinctionFloor())
                mass_[i] = 0.f;
            sum += mass_[i];
        }
        if (!(sum > 1e-20) || !std::isfinite(sum)) {
            mass_.fill(0.f);
            mass_[kBins / 2 - 1] = 0.5f;
            mass_[kBins / 2] = 0.5f;
            return;
        }
        // Scaling can move a value that was just above the floor below it. Each
        // repair pass removes at least one such bin, so kBins passes is a strict
        // bound; the next rescale only increases the remaining meaningful mass.
        for (int pass = 0; pass < kBins; ++pass) {
            const float inv = static_cast<float>(1.0 / sum);
            bool repaired = false;
            for (int i = 0; i < kBins; ++i) {
                mass_[i] *= inv;
                if (mass_[i] > 0.f && mass_[i] < numericalExtinctionFloor()) {
                    mass_[i] = 0.f;
                    repaired = true;
                }
            }
            if (!repaired)
                return;
            sum = 0.0;
            for (int i = 0; i < kBins; ++i)
                sum += mass_[i];
        }
    }

    float refinedTrait(int index) const {
        if (index <= 0 || index + 1 >= kBins)
            return traitAt(index);
        float left = mass_[index - 1];
        float center = mass_[index];
        float right = mass_[index + 1];
        float denominator = left - 2.f * center + right;
        float offset = 0.f;
        if (std::fabs(denominator) > 1e-12f)
            offset = 0.5f * (left - right) / denominator;
        return traitAt(index) + clamp(offset, -0.5f, 0.5f) * binWidth();
    }

    Peaks findPeaks() const {
        int maxima[kBins];
        int count = 0;
        for (int i = 0; i < kBins; ++i) {
            float left = i > 0 ? mass_[i - 1] : -1.f;
            float right = i + 1 < kBins ? mass_[i + 1] : -1.f;
            if (mass_[i] >= left && mass_[i] >= right &&
                (mass_[i] > left || mass_[i] > right))
                maxima[count++] = i;
        }

        Peaks result;
        int dominant = 0;
        for (int i = 1; i < kBins; ++i)
            if (mass_[i] > mass_[dominant]) dominant = i;
        result.lowIndex = result.highIndex = dominant;
        result.lowTrait = result.highTrait = refinedTrait(dominant);
        result.lowHeight = result.highHeight = mass_[dominant];

        float bestScore = -1.f;
        const int minBins = 4;
        for (int a = 0; a < count; ++a) {
            for (int b = a + 1; b < count; ++b) {
                int lo = maxima[a], hi = maxima[b];
                if (hi - lo < minBins)
                    continue;
                int valley = lo;
                for (int i = lo + 1; i < hi; ++i)
                    if (mass_[i] < mass_[valley]) valley = i;
                float smaller = std::min(mass_[lo], mass_[hi]);
                float depth = smaller > 0.f ? 1.f - mass_[valley] / smaller : 0.f;
                float score = smaller * std::max(depth, 0.f);
                if (score > bestScore) {
                    bestScore = score;
                    result.lowIndex = lo;
                    result.highIndex = hi;
                    result.valleyIndex = valley;
                    result.lowTrait = refinedTrait(lo);
                    result.highTrait = refinedTrait(hi);
                    result.lowHeight = mass_[lo];
                    result.highHeight = mass_[hi];
                    result.valleyHeight = mass_[valley];
                    result.pair = true;
                }
            }
        }

        if (result.pair) {
            result.lowMass = 0.f;
            result.highMass = 0.f;
            for (int i = 0; i < kBins; ++i) {
                if (i <= result.valleyIndex) result.lowMass += mass_[i];
                else                         result.highMass += mass_[i];
            }
        }
        return result;
    }

    void refreshMetrics(float elapsedTau, bool initializeDetector,
                        bool updateDetector = true) {
        bool oldSplitEvent = metrics_.splitEvent;
        bool oldMergeEvent = metrics_.mergeEvent;

        double mean = 0.0;
        for (int i = 0; i < kBins; ++i)
            mean += static_cast<double>(mass_[i]) * traitAt(i);
        double variance = 0.0;
        for (int i = 0; i < kBins; ++i) {
            double d = traitAt(i) - mean;
            variance += static_cast<double>(mass_[i]) * d * d;
        }

        Peaks peaks = findPeaks();
        float smaller = peaks.pair ? std::min(peaks.lowHeight, peaks.highHeight) : 0.f;
        float depth = smaller > 0.f ? 1.f - peaks.valleyHeight / smaller : 0.f;
        float separation = peaks.pair ? peaks.highTrait - peaks.lowTrait : 0.f;
        bool strongPair = peaks.pair && separation >= 4.f * binWidth() &&
                          depth >= 0.18f && peaks.lowMass >= 0.10f && peaks.highMass >= 0.10f;
        bool weakPair = peaks.pair && separation >= 3.5f * binWidth() &&
                        depth >= 0.08f && peaks.lowMass >= 0.07f && peaks.highMass >= 0.07f;

        if (!updateDetector) {
            oldSplitEvent = false;
            oldMergeEvent = false;
        }
        else if (initializeDetector) {
            split_ = strongPair;
            splitTimer_ = 0.f;
            mergeTimer_ = 0.f;
            oldSplitEvent = false;
            oldMergeEvent = false;
        }
        else if (!split_) {
            splitTimer_ = strongPair ? splitTimer_ + elapsedTau : 0.f;
            mergeTimer_ = 0.f;
            if (splitTimer_ >= splitPersistence()) {
                split_ = true;
                splitTimer_ = 0.f;
                oldSplitEvent = true;
            }
        }
        else {
            mergeTimer_ = weakPair ? 0.f : mergeTimer_ + elapsedTau;
            splitTimer_ = 0.f;
            if (mergeTimer_ >= mergePersistence()) {
                split_ = false;
                mergeTimer_ = 0.f;
                oldMergeEvent = true;
            }
        }

        metrics_.mean = static_cast<float>(mean);
        metrics_.spread = static_cast<float>(std::sqrt(std::max(variance, 0.0)));
        metrics_.split = split_;
        metrics_.splitEvent = oldSplitEvent;
        metrics_.mergeEvent = oldMergeEvent;

        if (split_ && peaks.pair) {
            metrics_.lowTrait = peaks.lowTrait;
            metrics_.highTrait = peaks.highTrait;
            metrics_.lowMass = peaks.lowMass;
            metrics_.highMass = peaks.highMass;
            metrics_.divergence = std::max(0.f, peaks.highTrait - peaks.lowTrait);
        }
        else {
            // Until a pair is accepted, both voices represent the unresolved
            // population as one. The mean stays continuous through a flat top
            // and does not choose the left member of an exactly symmetric,
            // still-too-shallow candidate pair.
            float trait = clamp(static_cast<float>(mean), traitMin(), traitMax());
            metrics_.lowTrait = trait;
            metrics_.highTrait = trait;
            metrics_.lowMass = 1.f;
            metrics_.highMass = 0.f;
            metrics_.divergence = 0.f;
        }
    }
};

} // namespace coalescent
