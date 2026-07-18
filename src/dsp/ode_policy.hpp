#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

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

// Rack builds DSP with -funsafe-math-optimizations, under which compilers may
// fold std::isfinite() to true. Inspecting the IEEE-754 exponent bits keeps the
// hostile-input guard effective under the actual production flags.
inline bool finiteFloat(float value) {
    static_assert(sizeof(float) == sizeof(std::uint32_t), "32-bit float required");
    static_assert(std::numeric_limits<float>::is_iec559, "IEEE-754 float required");
    std::uint32_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return (bits & 0x7f800000u) != 0x7f800000u;
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
