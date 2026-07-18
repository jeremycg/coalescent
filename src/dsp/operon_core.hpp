#pragma once

#include "fast_tanh.hpp"
#include "ode_policy.hpp"
#include "rk4.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace coalescent {
namespace operon {

static constexpr float RATE_CAL = 12.46f;
static constexpr float HSUB_MAX = 0.05f;
static constexpr int MIN_SUB = 2;
static constexpr int MAX_SUB = 64;
static constexpr float STATE_MAX = 1e3f;
static constexpr float BIAS_EPS = 1e-6f;
static constexpr float OUT_GAIN = 0.3f;
static constexpr float GATE_LOW = -0.05f;
static constexpr float GATE_HIGH = 0.05f;
static constexpr int CENTER_ITERS = 16;
static constexpr int NEWTON_ITERS = 4;

struct Parameters {
    float alpha;
    float n;
    float beta;
    float a0;
};

inline void seed(float state[6]) {
    state[0] = 0.2f;
    state[1] = 0.1f;
    state[2] = 0.3f;
    state[3] = 0.1f;
    state[4] = 0.4f;
    state[5] = 0.15f;
}

inline float solvePStar(float alpha, float n, float alpha0, float guess) {
    const float hi = std::max(1.f, alpha + alpha0 + 1.f);
    if (guess > 0.f && guess < hi) {
        float x = guess;
        for (int i = 0; i < NEWTON_ITERS; ++i) {
            const float safeX = std::max(x, 1e-6f);
            const float xn = std::pow(safeX, n);
            const float d = 1.f + xn;
            const float g = x - alpha / d - alpha0;
            const float gp = 1.f + alpha * n * (xn / safeX) / (d * d);
            const float step = g / gp;
            x -= step;
            if (!(x > 0.f && x < hi)) break;
            if (std::fabs(step) < 1e-5f) return x;
        }
    }

    float lo = 0.f;
    float high = hi;
    for (int i = 0; i < CENTER_ITERS; ++i) {
        const float mid = 0.5f * (lo + high);
        const float g = mid - alpha
            / (1.f + std::pow(std::max(mid, 0.f), n)) - alpha0;
        if (g > 0.f) high = mid;
        else lo = mid;
    }
    return 0.5f * (lo + high);
}

inline float directHillResponse(float repressor, float n) {
    return 1.f / (1.f + std::pow(repressor, n));
}

class HillLut {
public:
    static const int kSize = 8192;
    static const int kSlice = 256;
    static constexpr float xMax() { return 64.f; }
    static constexpr float moveEpsilon() { return 1e-4f; }
    static constexpr int settleSamples() { return 2048; }

    struct BuildWork {
        int entries;
        bool completed;
    };

    HillLut() { reset(); }

    void reset() {
        lutN_ = -1.f;
        settleRef_ = -1e9f;
        settleCount_ = 0;
        valid_ = false;
        buildPos_ = -1;
        buildN_ = -1.f;
    }

    static float entry(int index, float n) {
        return 1.f / (1.f + std::pow(xMax() * index / kSize, n));
    }

    BuildWork process(float n) {
        BuildWork work = {0, false};
        if (std::fabs(n - settleRef_) > moveEpsilon()) {
            settleRef_ = n;
            settleCount_ = 0;
            buildPos_ = -1;
        }
        else if (settleCount_ < settleSamples()) {
            ++settleCount_;
        }

        if (buildPos_ < 0 && settleCount_ >= settleSamples()
            && std::fabs(settleRef_ - lutN_) > moveEpsilon()) {
            buildPos_ = 0;
            buildN_ = settleRef_;
            valid_ = false;
        }

        if (buildPos_ >= 0) {
            const int end = std::min(buildPos_ + kSlice, kSize + 1);
            work.entries = end - buildPos_;
            for (; buildPos_ < end; ++buildPos_)
                values_[buildPos_] = entry(buildPos_, buildN_);
            if (buildPos_ > kSize) {
                lutN_ = buildN_;
                valid_ = true;
                buildPos_ = -1;
                settleCount_ = 0;
                work.completed = true;
            }
        }
        return work;
    }

    bool direct(float n) const {
        return !valid_ || std::fabs(n - lutN_) > moveEpsilon();
    }

    float lookup(float rep) const {
        if (rep >= xMax()) return values_[kSize];
        const float f = rep * (kSize / xMax());
        const int i = static_cast<int>(f);
        return values_[i] + (f - i) * (values_[i + 1] - values_[i]);
    }

    float response(float rep, float n, bool useDirect) const {
        return useDirect ? directHillResponse(rep, n) : lookup(rep);
    }

private:
    float values_[kSize + 1];
    float lutN_;
    float settleRef_;
    int settleCount_;
    bool valid_;
    int buildPos_;
    float buildN_;
};

template <typename HillResponse>
inline void derivative(const float state[6], float out[6],
                       const Parameters& p, float perturb,
                       HillResponse&& hillResponse) {
    for (int i = 0; i < 3; ++i) {
        const float rawRep = state[3 + ((i + 2) % 3)];
        const float rep = rawRep > 0.f ? rawRep : 0.f;
        const float hill = hillResponse(rep);
        out[i] = -state[i] + p.alpha * hill + p.a0;
        out[3 + i] = -p.beta * (state[3 + i] - state[i]);
    }
    out[0] += perturb;
    out[0] += BIAS_EPS;
    out[1] += -0.5f * BIAS_EPS;
    out[2] += -0.5f * BIAS_EPS;
}

template <typename HillResponse>
inline void step(float (&state)[6], float h, const Parameters& p,
                 float perturb, HillResponse&& hillResponse) {
    const auto deriv = [&](const float* value, float* out) {
        derivative(value, out, p, perturb, hillResponse);
    };
    coalescent::rk4<6>(state, h, deriv);
}

struct StateRepair {
    bool finite;
    bool clamped;
    bool exceededRange;

    explicit operator bool() const { return finite; }
};

// Validate and clamp one accepted RK4 substep. The diagnostic lets stability
// tests prove that normal trajectories do not lean on the safety ceiling.
inline StateRepair repairStateDetailed(float (&state)[6]) {
    StateRepair result = {true, false, false};
    for (int i = 0; i < 6; ++i) {
        if (!std::isfinite(state[i]))
            result.finite = false;
        else if (std::fabs(state[i]) > STATE_MAX)
            result.exceededRange = true;
        const float repaired = std::fmax(std::fmin(state[i], STATE_MAX), 0.f);
        result.clamped = result.clamped || repaired != state[i];
        state[i] = repaired;
    }
    return result;
}

inline bool repairState(float (&state)[6]) {
    return repairStateDetailed(state).finite;
}

template <typename HillResponse>
inline StateRepair advanceAcceptedSubstep(float (&state)[6], float h,
                                          const Parameters& p, float perturb,
                                          HillResponse&& hillResponse) {
    step(state, h, p, perturb, std::forward<HillResponse>(hillResponse));
    return repairStateDetailed(state);
}

inline OdeStepPlan stepPlan(float requestedDelta) {
    return makeOdeStepPlan(requestedDelta, HSUB_MAX, MIN_SUB, MAX_SUB);
}

inline bool gateStep(float centered, bool& armed) {
    if (centered <= GATE_LOW) {
        armed = true;
    }
    else if (armed && centered >= GATE_HIGH) {
        armed = false;
        return true;
    }
    return false;
}

inline float outputVoltage(float centered) {
    return 5.f * coalescent::fastTanh(centered * OUT_GAIN);
}

} // namespace operon
} // namespace coalescent
