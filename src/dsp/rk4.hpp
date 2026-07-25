#pragma once

// Generic explicit 4th-order Runge–Kutta step, shared by Coalescent's ODE
// modules — the neuron pair Axon/Soma (N = 2 / 3), Operon's repressilator (N = 6),
// and the ecological pair Bunnies/Foxes (N = 2 / 3). Each module supplies its own
// derivative via a small functor; only the dimensionality and f() differ, so the
// step itself lives here once.
//
// Templated on the scalar type T so one implementation serves scalar float (the
// per-voice path) and rack::simd::float_4 (four poly voices at once). For float_4
// the step size h is a per-lane vector: a voice group runs at the group's max
// substep count K, so each lane's h = subTau_lane / K — never larger than its own
// scalar h, i.e. accuracy is preserved lane-by-lane.
//
// `tools/integrator_equiv.cpp` checks this implementation against independent
// closed-form ODE solutions, so it does not preserve a second integrator copy.

namespace coalescent {

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

} // namespace coalescent
