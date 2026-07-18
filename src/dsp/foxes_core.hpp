#pragma once

#include "fast_tanh.hpp"
#include "ode_policy.hpp"
#include "rk4.hpp"

#include <algorithm>
#include <cmath>

namespace coalescent {
namespace foxes {

static constexpr double A1 = 5.0;
static constexpr double A2 = 0.1;
static constexpr double D1 = 0.4;
static constexpr double D2 = 0.01;
static constexpr float RATE_CAL = 61.387f;
static constexpr float HSUB_MAX = 0.1f;
static constexpr int MIN_SUB = 2;
static constexpr int MAX_SUB = 64;
static constexpr float STATE_MAX = 1e3f;
static constexpr float POS_FLOOR = 1e-4f;
static constexpr float KICK_GAIN = 0.25f;
static constexpr float POP_GATE = 0.15f;
static constexpr int TRAIL_N = 2048;
static constexpr float CAP_TAU = 0.225f;
static constexpr float SEED_X = 0.812781f;
static constexpr float SEED_Y = 0.104883f;
static constexpr float SEED_Z = 12.478951f;

template <typename T>
inline T b1FromWild(T wild) {
    return T(1) + T(5.2) * wild * wild;
}

template <typename T>
inline T b2FromBalance(T balance) {
    return T(1.75) + T(0.5) * balance;
}

template <typename T>
inline T f1(T x, T b1) {
    return T(A1) * x / (T(1) + b1 * x);
}

template <typename T>
inline T f2(T y, T b2) {
    return T(A2) * y / (T(1) + b2 * y);
}

template <typename T>
inline void derivative(const T state[3], T out[3], T b1, T b2, T kick) {
    const T x = std::max(state[0], T(POS_FLOOR));
    const T y = std::max(state[1], T(POS_FLOOR));
    const T z = std::max(state[2], T(POS_FLOOR));
    const T response1 = f1(x, b1);
    const T response2 = f2(y, b2);
    out[0] = x * (T(1) - x) - response1 * y + kick;
    out[1] = response1 * y - response2 * z - T(D1) * y;
    out[2] = response2 * z - T(D2) * z;
}

template <typename T>
inline void step(T (&state)[3], T h, T b1, T b2, T kick) {
    const auto deriv = [&](const T* value, T* out) {
        derivative(value, out, b1, b2, kick);
    };
    coalescent::rk4<3>(state, h, deriv);
}

template <typename T>
inline bool equilibrium(T b1, T b2, T& x, T& y, T& z) {
    const T denominator = T(A2) - b2 * T(D2);
    if (!(denominator > T(0))) return false;
    y = T(D2) / denominator;
    const T discriminant = (b1 + T(1)) * (b1 + T(1))
        - T(4) * T(A1) * b1 * y;
    if (!(discriminant >= T(0))) return false;
    x = ((b1 - T(1)) + std::sqrt(discriminant)) / (T(2) * b1);
    const T response = f1(x, b1);
    z = y * (response - T(D1)) / T(D2);
    return x > T(0) && y > T(0) && response > T(D1) && z > T(0)
        && std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
}

inline void seed(float state[3]) {
    state[0] = SEED_X;
    state[1] = SEED_Y;
    state[2] = SEED_Z;
}

inline OdeStepPlan stepPlan(float requestedDelta) {
    return makeOdeStepPlan(requestedDelta, HSUB_MAX, MIN_SUB, MAX_SUB);
}

inline bool repairState(float (&state)[3]) {
    if (!std::isfinite(state[0]) || !std::isfinite(state[1])
        || !std::isfinite(state[2])) return false;
    for (int i = 0; i < 3; ++i)
        state[i] = std::fmax(std::fmin(state[i], STATE_MAX), POS_FLOOR);
    return true;
}

template <typename T>
inline bool peakStep(T gained, T& previous, bool& rising) {
    bool fired = false;
    if (gained > previous + T(1e-7f)) {
        rising = true;
    }
    else if (gained < previous - T(1e-7f)) {
        fired = rising && previous > T(POP_GATE);
        rising = false;
    }
    previous = gained;
    return fired;
}

inline float gain(int species) {
    return species == 0 ? 8.f : (species == 1 ? 18.f : 2.5f);
}

inline float outputVoltage(float centered, int species) {
    return 5.f * coalescent::fastTanh(centered * gain(species));
}

} // namespace foxes
} // namespace coalescent
