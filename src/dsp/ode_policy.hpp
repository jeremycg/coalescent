#pragma once

#include "finite.hpp"

#include <algorithm>
#include <cmath>

namespace coalescent {

// Bounded fixed-step plan shared by the scalar ODE modules. Callers supply a
// non-negative requested simulation-time advance and their validated numerical
// limits; the returned step always satisfies h <= maxStep.
struct OdeStepPlan {
    float delta;
    int count;
    float h;
};

inline OdeStepPlan makeOdeStepPlan(float requested, float maxStep,
                                   int minCount, int maxCount) {
    const float delta = std::min(requested, maxStep * maxCount);
    const int count = std::min(maxCount, std::max(minCount,
        static_cast<int>(std::ceil(delta / maxStep))));
    OdeStepPlan plan = {delta, count, delta / count};
    return plan;
}

// Retained spelling for the float guard; coalescent::isFinite in finite.hpp is
// the shared predicate and documents why the exponent bits are inspected
// directly instead of calling the standard isfinite.
inline bool finiteFloat(float value) {
    return isFinite(value);
}

inline float finiteOr(float value, float fallback) {
    return finiteFloat(value) ? value : fallback;
}

// Generic form for Rack SIMD and other vector types. Callers provide the
// platform's finite predicate and lane-wise selector so this header stays
// SDK-free and non-finite lanes can fall back independently.
template <typename T, typename IsFinite, typename Select>
inline T finiteOr(T value, T fallback, IsFinite isFinite, Select select) {
    return select(isFinite(value), value, fallback);
}

inline float finiteClamp(float value, float fallback, float minimum, float maximum) {
    const float finite = finiteOr(value, fallback);
    return std::fmax(std::fmin(finite, maximum), minimum);
}

} // namespace coalescent
