#pragma once

namespace coalescent {

// SDK-free scalar Pade[7/6] tanh approximation. The Rack-facing header adds a
// SIMD overload; standalone model tests include this scalar definition directly.
inline float fastTanh(float x) {
    if (x < -4.f) return -1.f;
    if (x >  4.f) return  1.f;
    const float x2 = x * x;
    const float n = x * (135135.f + x2 * (17325.f + x2 * (378.f + x2)));
    const float d = 135135.f + x2 * (62370.f + x2 * (3150.f + x2 * 28.f));
    return n / d;
}

} // namespace coalescent
