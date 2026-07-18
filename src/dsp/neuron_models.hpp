#pragma once

#include "ode_policy.hpp"
#include "rk4.hpp"

#include <cmath>

// SDK-free model core shared by Axon, Soma, and their standalone stability
// tests. Rack-specific trigger, filtering, output, and display plumbing stays in
// the module wrappers; this file owns the state transition and its scheduling
// and safety boundaries.

namespace coalescent {
namespace neuron {

// Rack's dsp::FREQ_C4 value, kept SDK-free so standalone tools use the exact
// reference frequency that the production wrappers feed into the model.
static constexpr float PITCH_REFERENCE_HZ = 261.6256f;

struct ScalarSchedule {
    float subTau;
    int substeps;
    float h;
};

struct ScalarMin {
    float operator()(float a, float b) const { return a < b ? a : b; }
};

struct ScalarCeil {
    float operator()(float x) const { return std::ceil(x); }
};

struct ScalarLane {
    float operator()(float value, int) const { return value; }
};

struct ScalarAbs {
    float operator()(float value) const { return std::fabs(value); }
};

struct ScalarSelect {
    float operator()(bool condition, float yes, float no) const {
        return condition ? yes : no;
    }
};

struct ScalarClamp {
    float operator()(float value, float low, float high) const {
        return value < low ? low : (value > high ? high : value);
    }
};

struct IgnoreObservation {
    template <typename T>
    inline void operator()(const T&) const {}
};

// Pitch is an exponent before Rack's approximate exp2. NaN is neutral (C4),
// while infinities and huge finite values saturate at the same safe bounds as
// the production wrappers.
template <typename T, typename Select, typename Clamp>
inline T sanitizePitchExponent(T exponent, Select select, Clamp clampValue) {
    exponent = select(exponent == exponent, exponent, T(0.f));
    return clampValue(exponent, -30.f, 30.f);
}

inline float sanitizePitchExponent(float exponent) {
    return sanitizePitchExponent(exponent, ScalarSelect(), ScalarClamp());
}

// Return the bounded amount of dimensionless model time represented by one
// oversampled output. MinValue is std::min-like for scalar tests and simd::fmin
// in Rack, preserving the wrapper's SIMD operation order.
template <typename Model, typename T, typename MinValue>
inline T boundedSubTau(T pitchHz, float sampleRate, int oversample, MinValue minValue) {
    const T unbounded = Model::RATE_CAL * pitchHz / sampleRate / (float) oversample;
    return minValue(unbounded, T(Model::HSUB_MAX * Model::MAX_SUB));
}

// A SIMD group runs at its largest required K. Slower lanes retain their own
// subTau and therefore use a smaller h; no lane loses scalar accuracy.
template <typename Model, typename T, typename CeilValue, typename LaneValue>
inline int groupSubstepsWithMinimum(T subTau, int lanes, int minimum,
                                    CeilValue ceilValue, LaneValue laneValue) {
    const T required = ceilValue(subTau / Model::HSUB_MAX);
    int substeps = minimum;
    for (int lane = 0; lane < lanes; ++lane) {
        const int laneSubsteps = (int) laneValue(required, lane);
        if (laneSubsteps > substeps)
            substeps = laneSubsteps;
    }
    return substeps < Model::MAX_SUB ? substeps : Model::MAX_SUB;
}

template <typename Model, typename T, typename CeilValue, typename LaneValue>
inline int groupSubsteps(T subTau, int lanes, CeilValue ceilValue, LaneValue laneValue) {
    return groupSubstepsWithMinimum<Model>(
        subTau, lanes, Model::MIN_SUB, ceilValue, laneValue);
}

template <typename Model>
inline ScalarSchedule scalarScheduleWithMinimum(float pitchHz, float sampleRate,
                                                int oversample, int minimum) {
    ScalarSchedule schedule;
    schedule.subTau = boundedSubTau<Model>(pitchHz, sampleRate, oversample, ScalarMin());
    schedule.substeps = groupSubstepsWithMinimum<Model>(
        schedule.subTau, 1, minimum, ScalarCeil(), ScalarLane());
    schedule.h = schedule.subTau / (float) schedule.substeps;
    return schedule;
}

template <typename Model>
inline ScalarSchedule scalarSchedule(float pitchHz, float sampleRate, int oversample = 1) {
    return scalarScheduleWithMinimum<Model>(pitchHz, sampleRate, oversample, Model::MIN_SUB);
}

template <typename Model>
inline ScalarSchedule scalarScheduleFromSubTau(float unboundedSubTau) {
    ScalarSchedule schedule;
    schedule.subTau = ScalarMin()(unboundedSubTau, Model::HSUB_MAX * Model::MAX_SUB);
    schedule.substeps = groupSubsteps<Model>(schedule.subTau, 1, ScalarCeil(), ScalarLane());
    schedule.h = schedule.subTau / (float) schedule.substeps;
    return schedule;
}

struct AxonCore {
    static const int STATE_COUNT = 2;
    static constexpr float RATE_CAL = 37.899004f;
    static constexpr float B_FIXED = 0.8f;
    static constexpr float HSUB_MAX = 0.05f;
    static const int MIN_SUB = 2;
    static const int MAX_SUB = 64;
    static constexpr float SPIKE_THRESH = 1.f;
    static constexpr float STATE_MAX = 10.f;
    static constexpr float RUNAWAY_MAX = 1e6f;
    static constexpr float REST_V = -1.2f;
    static constexpr float REST_W = -0.6f;

    template <typename T>
    static inline void derivative(const T* state, T* out, T current, T eps, T shape) {
        const T v = state[0];
        const T w = state[1];
        out[0] = v - v * v * v / 3.f - w + current;
        out[1] = eps * (v + shape - B_FIXED * w);
    }

    // One accepted observation interval. Axon currently observes events once
    // after this complete K-substep transition, not after each RK4 substep.
    template <typename T>
    static inline void advanceObservation(T (&state)[STATE_COUNT], T h, int substeps,
                                          T current, T eps, T shape) {
        for (int k = 0; k < substeps; ++k) {
            coalescent::rk4<STATE_COUNT>(state, h, [&](const T* y, T* d) {
                derivative(y, d, current, eps, shape);
            });
        }
    }

    template <typename T, typename Abs, typename Select, typename Clamp>
    static inline auto repair(T (&state)[STATE_COUNT], Abs absValue, Select select,
                              Clamp clampValue) -> decltype(state[0] == state[0]) {
        const auto finite = (state[0] == state[0]) & (state[1] == state[1])
            & (absValue(state[0]) < RUNAWAY_MAX) & (absValue(state[1]) < RUNAWAY_MAX);
        state[0] = select(finite, state[0], T(REST_V));
        state[1] = select(finite, state[1], T(REST_W));
        state[0] = clampValue(state[0], -STATE_MAX, STATE_MAX);
        state[1] = clampValue(state[1], -STATE_MAX, STATE_MAX);
        return finite;
    }

    static inline bool repair(float (&state)[STATE_COUNT]) {
        return repair(state, ScalarAbs(), ScalarSelect(), ScalarClamp());
    }
};

struct SomaCore {
    static const int STATE_COUNT = 3;
    static constexpr float A = 1.f;
    static constexpr float B = 3.f;
    static constexpr float C = 1.f;
    static constexpr float D = 5.f;
    static constexpr float XR = -1.6f;
    static constexpr float RATE_CAL = 55.364003f;
    static constexpr float HSUB_MAX = 0.05f;
    static const int MIN_SUB = 2;
    static const int MAX_SUB = 64;
    static constexpr float SPIKE_THRESH = 1.f;
    static constexpr float STATE_MAX = 25.f;
    static constexpr float RUNAWAY_MAX = 1e6f;
    static constexpr float R_MIN = 0.001f;
    static constexpr float R_MAX = 0.05f;
    static constexpr float R_DEFAULT = 0.03f;
    static constexpr float BURST_EXP_MIN = -30.f;
    static constexpr float BURST_EXP_MAX = 30.f;
    static constexpr float REST_X = -1.6f;
    static constexpr float REST_Y = -11.8f;
    static constexpr float REST_Z = 2.f;

    template <typename T>
    static inline void derivative(const T* state, T* out, T current, T r, T adapt) {
        const T x = state[0];
        const T y = state[1];
        const T z = state[2];
        out[0] = y - A * x * x * x + B * x * x - z + current;
        out[1] = C - D * x * x - y;
        out[2] = r * (adapt * (x - XR) - z);
    }

    // Advance one observation interval and expose every accepted RK4 substep to
    // event detectors. This keeps the transition single-sourced while allowing
    // narrow spikes to rise and fall between output-rate observations.
    template <typename T, typename Observer>
    static inline void advanceObservation(T (&state)[STATE_COUNT], T h, int substeps,
                                          T current, T r, T adapt, Observer observer) {
        for (int k = 0; k < substeps; ++k) {
            coalescent::rk4<STATE_COUNT>(state, h, [&](const T* y, T* d) {
                derivative(y, d, current, r, adapt);
            });
            observer(state[0]);
        }
    }

    template <typename T>
    static inline void advanceObservation(T (&state)[STATE_COUNT], T h, int substeps,
                                          T current, T r, T adapt) {
        advanceObservation(
            state, h, substeps, current, r, adapt, IgnoreObservation());
    }

    template <typename T, typename Abs, typename Select, typename Clamp>
    static inline auto repair(T (&state)[STATE_COUNT], Abs absValue, Select select,
                              Clamp clampValue) -> decltype(state[0] == state[0]) {
        const auto finite = (state[0] == state[0]) & (state[1] == state[1]) & (state[2] == state[2])
            & (absValue(state[0]) < RUNAWAY_MAX) & (absValue(state[1]) < RUNAWAY_MAX)
            & (absValue(state[2]) < RUNAWAY_MAX);
        state[0] = select(finite, state[0], T(REST_X));
        state[1] = select(finite, state[1], T(REST_Y));
        state[2] = select(finite, state[2], T(REST_Z));
        state[0] = clampValue(state[0], -STATE_MAX, STATE_MAX);
        state[1] = clampValue(state[1], -STATE_MAX, STATE_MAX);
        state[2] = clampValue(state[2], -STATE_MAX, STATE_MAX);
        return finite;
    }

    static inline bool repair(float (&state)[STATE_COUNT]) {
        return repair(state, ScalarAbs(), ScalarSelect(), ScalarClamp());
    }

    static inline float sanitizeBurstBase(float exponent) {
        const float fallback = std::log2(R_DEFAULT);
        if (!std::isfinite(exponent))
            return fallback;
        return ScalarClamp()(exponent, BURST_EXP_MIN, BURST_EXP_MAX);
    }

    template <typename T, typename Finite, typename Select, typename Clamp>
    static inline T sanitizeBurstExponent(T exponent, T fallback, Finite finite,
                                          Select select, Clamp clampValue) {
        return select(finite(exponent),
                      clampValue(exponent, BURST_EXP_MIN, BURST_EXP_MAX), fallback);
    }

    static inline float sanitizeBurstExponent(float exponent, float fallback) {
        struct ScalarFinite {
            bool operator()(float value) const { return std::isfinite(value); }
        };
        return sanitizeBurstExponent(exponent, fallback, ScalarFinite(), ScalarSelect(), ScalarClamp());
    }
};

} // namespace neuron
} // namespace coalescent
