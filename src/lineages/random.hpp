#pragma once

#include <cmath>
#include <cstdint>

namespace coalescent {
namespace lineages {

// Module-local PCG-XSH-RR 64/32. This deliberately matches the engine used by
// Islands so generated state is independent of Rack and standard-library RNGs.
class Pcg32 {
public:
    struct State {
        std::uint64_t state;
        std::uint64_t increment;

        State(std::uint64_t stateValue = 0u, std::uint64_t incrementValue = 1u)
            : state(stateValue), increment(incrementValue) {}
    };

    struct PoissonResult {
        std::uint32_t count;
        bool truncated;

        PoissonResult(std::uint32_t countValue = 0u, bool wasTruncated = false)
            : count(countValue), truncated(wasTruncated) {}
    };

    Pcg32(std::uint64_t seedValue = 42u, std::uint64_t stream = 54u) {
        seed(seedValue, stream);
    }

    void seed(std::uint64_t seedValue, std::uint64_t stream) {
        state_ = 0u;
        increment_ = (stream << 1u) | 1u;
        next();
        state_ += seedValue;
        next();
    }

    std::uint32_t next() {
        const std::uint64_t oldState = state_;
        state_ = oldState * UINT64_C(6364136223846793005) + increment_;
        const std::uint32_t shifted = static_cast<std::uint32_t>(
            ((oldState >> 18u) ^ oldState) >> 27u);
        const std::uint32_t rotation = static_cast<std::uint32_t>(oldState >> 59u);
        return (shifted >> rotation) | (shifted << ((-rotation) & 31u));
    }

    State state() const {
        return State(state_, increment_);
    }

    static bool validState(const State& saved) {
        return (saved.increment & 1u) != 0u;
    }

    bool restore(const State& saved) {
        if (!validState(saved))
            return false;
        state_ = saved.state;
        increment_ = saved.increment;
        return true;
    }

    // Rejection avoids the modulo bias of next() % bound. Generation only uses
    // bounds <= 16, whose reject probability is at most 9 / 2^32. The fixed
    // retry budget keeps NEW bounded even under a corrupt/adversarial stream;
    // multiply-high is the deterministic last-resort mapping after 32 rejects.
    std::uint32_t bounded(std::uint32_t bound) {
        if (bound <= 1u)
            return 0u;
        const std::uint32_t threshold = static_cast<std::uint32_t>(-bound) % bound;
        for (int attempt = 0; attempt < 32; ++attempt) {
            const std::uint32_t value = next();
            if (value >= threshold)
                return value % bound;
        }
        return static_cast<std::uint32_t>(
            (static_cast<std::uint64_t>(next()) * bound) >> 32u);
    }

    // Exact 32-bit grid on [0, 1). This is suitable for categorical sampling.
    double uniformClosedOpen() {
        return std::ldexp(static_cast<double>(next()), -32);
    }

    // Midpoints of the same grid lie strictly inside (0, 1), making logarithms
    // and branch-interior sampling safe without retry loops.
    double uniformOpen() {
        return std::ldexp(static_cast<double>(next()) + 0.5, -32);
    }

    double exponential(double rate) {
        if (!(rate > 0.0) || !std::isfinite(rate))
            return 0.0;
        return -std::log(uniformOpen()) / rate;
    }

    // Knuth's exact product sampler is practical for Lineages (lambda <= 64).
    // Folding the tail into cap makes both work and storage strictly bounded.
    PoissonResult poissonCapped(double lambda, std::uint32_t cap) {
        if (!(lambda > 0.0) || !std::isfinite(lambda) || cap == 0u)
            return PoissonResult();

        const double limit = std::exp(-lambda);
        double product = 1.0;
        for (std::uint32_t count = 0u; count < cap; ++count) {
            product *= uniformOpen();
            if (product <= limit)
                return PoissonResult(count, false);
        }
        product *= uniformOpen();
        return PoissonResult(cap, product > limit);
    }

private:
    std::uint64_t state_ = 0u;
    std::uint64_t increment_ = 1u;
};

} // namespace lineages
} // namespace coalescent
