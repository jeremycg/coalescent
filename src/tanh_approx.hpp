#pragma once

namespace coalescent {

// Padé[7/6] rational approximation of tanh, saturating to ±1 beyond ±4.
//
// Used as the output soft-clip in Axon/Soma/Haptik. Profiled ~14x faster than
// std::tanh with a max error of 6.7e-4 (-64 dB) over |x| <= 8 — inaudible — so
// it's a free substitute for the (transcendental) library tanh in the oversampled
// output chain. See tools/perf_minsub.cpp / the tanh profiling for the numbers.
inline float fastTanh(float x) {
    if (x < -4.f) return -1.f;
    if (x >  4.f) return  1.f;
    const float x2 = x * x;
    const float n = x * (135135.f + x2 * (17325.f + x2 * (378.f + x2)));
    const float d = 135135.f + x2 * (62370.f + x2 * (3150.f + x2 * 28.f));
    return n / d;
}

} // namespace coalescent
