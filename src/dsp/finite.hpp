#pragma once

#include <cstdint>
#include <cstring>
#include <limits>

namespace coalescent {

// Toolchain-independent replacement for the standard isfinite predicate.
//
// Rack builds plugin DSP with -O3 -funsafe-math-optimizations. That flag set
// alone does not enable -ffinite-math-only, and GCC 15 was measured to keep the
// standard predicate honest under exactly those flags. The guarantee is a
// property of the compiler rather than of the C++ standard, though: a toolchain
// that folds non-finite cases away (anything reaching -ffinite-math-only, which
// -ffast-math implies) is free to reduce that predicate to a constant true, and
// every hostile-input and state-repair path in this plugin would then stop
// rejecting NaN and infinity while still looking correct in review.
//
// Inspecting the IEEE-754 exponent field cannot be folded that way, costs a
// move and a mask, and returns exactly what the standard predicate returns for
// every input: true for zero, subnormals, and normal values; false for infinity
// and NaN. These overloads are therefore a drop-in replacement, and using them
// uniformly keeps the plugin's numerical trust boundaries independent of which
// compiler and flag set produced the binary.

namespace finite_detail {

inline std::uint32_t bitsOf(float value) {
    static_assert(sizeof(float) == sizeof(std::uint32_t), "32-bit float required");
    static_assert(std::numeric_limits<float>::is_iec559, "IEEE-754 float required");
    std::uint32_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

inline std::uint64_t bitsOf(double value) {
    static_assert(sizeof(double) == sizeof(std::uint64_t), "64-bit double required");
    static_assert(std::numeric_limits<double>::is_iec559, "IEEE-754 double required");
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

} // namespace finite_detail

inline bool isFinite(float value) {
    return (finite_detail::bitsOf(value) & UINT32_C(0x7f800000))
        != UINT32_C(0x7f800000);
}

inline bool isFinite(double value) {
    return (finite_detail::bitsOf(value) & UINT64_C(0x7ff0000000000000))
        != UINT64_C(0x7ff0000000000000);
}

// Pitch sanitizers distinguish NaN (neutral) from infinity (saturate by sign).
// Compare the sign-stripped representation with positive infinity so the test
// remains valid when the compiler assumes floating-point values are finite.
inline bool isNaN(float value) {
    return (finite_detail::bitsOf(value) & UINT32_C(0x7fffffff))
        > UINT32_C(0x7f800000);
}

inline bool isNaN(double value) {
    return (finite_detail::bitsOf(value) & UINT64_C(0x7fffffffffffffff))
        > UINT64_C(0x7ff0000000000000);
}

} // namespace coalescent
