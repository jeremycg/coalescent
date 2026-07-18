#pragma once

#include "dsp/fast_tanh.hpp"

namespace coalescent {

// Padé[7/6] rational approximation of tanh, saturating to ±1 beyond ±4.
//
// Used as the output soft-clip in Axon/Soma/Haptik/Bunnies/Foxes. Profiled ~14x faster than
// std::tanh with a max error of 6.7e-4 (-64 dB) over |x| <= 8 — inaudible — so
// it's a free substitute for the (transcendental) library tanh in the oversampled
// output chain. See tools/perf_minsub.cpp / the tanh profiling for the numbers.
// Four-lane version (poly path). Clamp the domain to ±4 as above (Padé at ±4 is
// 0.9993 vs the scalar's exact ±1 — a -72 dB difference, inaudible).
inline simd::float_4 fastTanh(simd::float_4 x) {
    x = simd::fmin(simd::fmax(x, -4.f), 4.f);
    const simd::float_4 x2 = x * x;
    const simd::float_4 n = x * (135135.f + x2 * (17325.f + x2 * (378.f + x2)));
    const simd::float_4 d = 135135.f + x2 * (62370.f + x2 * (3150.f + x2 * 28.f));
    return n / d;
}

} // namespace coalescent
