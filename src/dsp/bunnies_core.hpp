#pragma once

#include "fast_tanh.hpp"
#include "ode_policy.hpp"
#include "rk4.hpp"

#include <algorithm>
#include <cmath>

namespace coalescent {
namespace bunnies {

enum Mode { LV = 0, RM = 1 };

static constexpr float RATE_CAL = 7.33f;
static constexpr float HSUB_MAX = 0.05f;
static constexpr int MIN_SUB = 2;
static constexpr int MAX_SUB = 64;
static constexpr float OUT_GAIN = 0.9f;
static constexpr float STATE_MAX = 1e3f;
static constexpr float POS_FLOOR = 1e-4f;
static constexpr float STAB_K = 0.5f;
static constexpr float STAB_FLOOR = 0.2f;
static constexpr float MAX_STAB_STEP = 0.25f;
static constexpr float LV_V0_RANGE = 3.5f;
static constexpr float KICK_GAIN = 0.5f;
static constexpr float POP_MIN = 0.05f;
static constexpr float RM_B = 0.5f;
static constexpr float RM_S = 1.f;

inline float gammaFromBalance(float balance) {
    return 0.2f + balance * 4.8f;
}

inline float targetInvariant(float gamma, float wild) {
    return gamma + 1.f + wild * LV_V0_RANGE;
}

inline float mortalityFromBalance(float balance) {
    return std::min(0.15f + balance * 0.45f, 0.95f / RM_B);
}

inline float capacityFromWild(float wild) {
    return 1.2f + wild * 10.8f;
}

inline bool rmCenter(float capacity, float mortality, float& x, float& y) {
    x = mortality / (1.f - RM_B * mortality);
    y = x * (1.f - x / capacity) / mortality;
    return x > 0.f && y > 0.f && std::isfinite(x) && std::isfinite(y);
}

inline void seed(float centerX, float centerY, float state[2]) {
    state[0] = centerX * 1.1f + 0.05f;
    state[1] = centerY * 0.9f + 0.02f;
}

struct Parameters {
    int mode;
    float gamma;
    float capacity;
    float mortality;
    float kick;
};

inline void derivative(const float state[2], float out[2],
                       const Parameters& p) {
    const float x = std::max(state[0], POS_FLOOR);
    const float y = std::max(state[1], POS_FLOOR);
    if (p.mode == LV) {
        out[0] = x * (1.f - y);
        out[1] = p.gamma * y * (x - 1.f);
    }
    else {
        const float response = x / (1.f + RM_B * x);
        out[0] = x * (1.f - x / p.capacity) - y * response;
        out[1] = RM_S * y * (response - p.mortality);
    }
    out[0] += p.kick;
}

inline void step(float (&state)[2], float h, const Parameters& p) {
    const auto deriv = [&](const float* value, float* out) {
        derivative(value, out, p);
    };
    coalescent::rk4<2>(state, h, deriv);
}

inline OdeStepPlan stepPlan(float requestedDelta, int mode, float gamma) {
    if (mode == LV) requestedDelta /= std::sqrt(gamma);
    return makeOdeStepPlan(requestedDelta, HSUB_MAX, MIN_SUB, MAX_SUB);
}

inline bool repairState(float (&state)[2]) {
    if (!std::isfinite(state[0]) || !std::isfinite(state[1])) return false;
    state[0] = std::fmax(std::fmin(state[0], STATE_MAX), POS_FLOOR);
    state[1] = std::fmax(std::fmin(state[1], STATE_MAX), POS_FLOOR);
    return true;
}

struct ServoResult {
    float invariant;
    float rawX;
    float rawY;
    float appliedX;
    float appliedY;
};

inline ServoResult applyLvServo(float (&state)[2], float gamma,
                                float target, float delta) {
    const float x = state[0];
    const float y = state[1];
    const float invariant = gamma * (x - std::log(x)) + (y - std::log(y));
    const float gradX = gamma * (1.f - 1.f / std::max(x, STAB_FLOOR));
    const float gradY = 1.f - 1.f / std::max(y, STAB_FLOOR);
    const float strength = STAB_K * delta;
    const float rawX = -strength * (invariant - target) * gradX;
    const float rawY = -strength * (invariant - target) * gradY;
    const float appliedX = std::fmax(std::fmin(rawX, MAX_STAB_STEP), -MAX_STAB_STEP);
    const float appliedY = std::fmax(std::fmin(rawY, MAX_STAB_STEP), -MAX_STAB_STEP);
    state[0] = std::max(state[0] + appliedX, POS_FLOOR);
    state[1] = std::max(state[1] + appliedY, POS_FLOOR);
    ServoResult result = {invariant, rawX, rawY, appliedX, appliedY};
    return result;
}

inline bool peakStep(float centered, float& previous, bool& rising) {
    bool fired = false;
    if (centered > previous + 1e-7f) {
        rising = true;
    }
    else if (centered < previous - 1e-7f) {
        fired = rising && previous > POP_MIN;
        rising = false;
    }
    previous = centered;
    return fired;
}

inline float outputVoltage(float centered) {
    return 5.f * coalescent::fastTanh(centered * OUT_GAIN);
}

} // namespace bunnies
} // namespace coalescent
