#pragma once

// Shared integrator for the Coalescent neuron pair (Axon / Soma).
//
// Both modules are stiff ODE oscillators integrated with the *same* scheme:
// classic explicit 4th-order Runge–Kutta, taken in pitch-adaptive substeps in
// dimensionless time. Only the dimensionality and the derivative differ —
// Axon is FitzHugh–Nagumo (N = 2: v, w), Soma is Hindmarsh–Rose (N = 3:
// x, y, z) — so the step itself lives here as a generic template and each
// module supplies its own f() via a small functor.
//
// The per-component arithmetic below is byte-for-byte the same expression that
// the two hand-written steppers used before extraction (k1..k4 with weights
// 1, 2, 2, 1 over 6), so audio/calibration are unchanged; see
// tools/integrator_equiv.cpp for the numeric proof.
//
// Templated on the scalar type T so the same code serves scalar float (the
// per-voice path) and rack::simd::float_4 (four poly voices at once). For float_4
// the step size h is a per-lane vector: a voice group runs at the group's max
// substep count K, so each lane's h = subTau_lane / K — never larger than its own
// scalar h, i.e. accuracy is preserved lane-by-lane.

namespace neuron {

// One RK4 step of size h over the N-dimensional state y, in place.
// deriv(const T* y, T* out) writes the time-derivative of state y.
template <int N, typename T, typename Deriv>
inline void rk4(T (&y)[N], T h, Deriv&& deriv) {
    T k1[N], k2[N], k3[N], k4[N], t[N];
    deriv(y, k1);
    for (int i = 0; i < N; ++i) t[i] = y[i] + 0.5f * h * k1[i];
    deriv(t, k2);
    for (int i = 0; i < N; ++i) t[i] = y[i] + 0.5f * h * k2[i];
    deriv(t, k3);
    for (int i = 0; i < N; ++i) t[i] = y[i] + h * k3[i];
    deriv(t, k4);
    for (int i = 0; i < N; ++i)
        y[i] += h / 6.f * (k1[i] + 2.f * k2[i] + 2.f * k3[i] + k4[i]);
}

} // namespace neuron
