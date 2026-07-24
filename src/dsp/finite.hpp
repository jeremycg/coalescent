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

inline bool isFinite(float value) {
    static_assert(sizeof(float) == sizeof(std::uint32_t), "32-bit float required");
    static_assert(std::numeric_limits<float>::is_iec559, "IEEE-754 float required");
    std::uint32_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return (bits & UINT32_C(0x7f800000)) != UINT32_C(0x7f800000);
}

inline bool isFinite(double value) {
    static_assert(sizeof(double) == sizeof(std::uint64_t), "64-bit double required");
    static_assert(std::numeric_limits<double>::is_iec559, "IEEE-754 double required");
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return (bits & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000);
}

} // namespace coalescent
