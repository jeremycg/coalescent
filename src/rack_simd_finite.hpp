#pragma once

#include <simd/Vector.hpp>

namespace coalescent {

// Rack's float_4 comparison operators currently map directly to SSE
// instructions, but classifying the IEEE-754 bits makes that contract explicit
// and independent of floating-point assumption flags.
inline rack::simd::float_4 simdFiniteMask(rack::simd::float_4 value) {
    const rack::simd::int32_4 magnitude =
        rack::simd::int32_4::cast(value) & rack::simd::int32_4(0x7fffffff);
    return rack::simd::float_4::cast(
        magnitude < rack::simd::int32_4(0x7f800000));
}

inline rack::simd::float_4 simdNaNMask(rack::simd::float_4 value) {
    const rack::simd::int32_4 magnitude =
        rack::simd::int32_4::cast(value) & rack::simd::int32_4(0x7fffffff);
    return rack::simd::float_4::cast(
        magnitude > rack::simd::int32_4(0x7f800000));
}

} // namespace coalescent
