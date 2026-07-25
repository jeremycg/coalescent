// SDK-free numerical contract for the production coalescent::rk4 template.
//
// This deliberately uses ODEs with closed-form solutions instead of retaining
// copies of any module's former handwritten stepper. It checks fourth-order
// convergence, multi-state operation, and long-run oscillator energy.
#include "../src/dsp/rk4.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>


template <typename T>
static T decayError(T step) {
    const T rate = T(0.7);
    const T endTime = T(8);
    const int count = static_cast<int>(std::round(endTime / step));
    step = endTime / T(count);
    T state[1] = {T(1.25)};
    for (int i = 0; i < count; ++i) {
        coalescent::rk4<1>(state, step, [&](const T* value, T* out) {
            out[0] = -rate * value[0];
        });
    }
    const T exact = T(1.25) * std::exp(-rate * endTime);
    return std::fabs(state[0] - exact);
}


static bool checkDecayConvergence() {
    const double coarse = decayError<double>(0.4);
    const double medium = decayError<double>(0.2);
    const double fine = decayError<double>(0.1);
    const double ratio1 = coarse / medium;
    const double ratio2 = medium / fine;
    const bool ok = fine < 1e-8 && ratio1 > 12.0 && ratio2 > 12.0;
    std::printf(
        "RK4 decay convergence: errors %.3e -> %.3e -> %.3e, ratios %.2f/%.2f: %s\n",
        coarse, medium, fine, ratio1, ratio2, ok ? "PASS" : "FAIL");
    return ok;
}


static bool checkIndependentStates() {
    const float rates[6] = {0.1f, 0.25f, 0.5f, 0.75f, 1.f, 1.5f};
    float state[6] = {1.f, 0.5f, 1.5f, 0.75f, 2.f, 1.25f};
    float initial[6];
    std::copy(state, state + 6, initial);
    const float step = 0.005f;
    const int count = 400;
    for (int i = 0; i < count; ++i) {
        coalescent::rk4<6>(state, step, [&](const float* value, float* out) {
            for (int lane = 0; lane < 6; ++lane)
                out[lane] = -rates[lane] * value[lane];
        });
    }

    float worst = 0.f;
    for (int lane = 0; lane < 6; ++lane) {
        const float exact = initial[lane] * std::exp(-rates[lane] * step * count);
        worst = std::max(worst, std::fabs(state[lane] - exact));
    }
    const bool ok = std::isfinite(worst) && worst < 2e-6f;
    std::printf("RK4 six-state analytic error %.3e: %s\n",
                worst, ok ? "PASS" : "FAIL");
    return ok;
}


static bool checkOscillatorEnergy() {
    const double omega = 1.7;
    const double step = 0.01;
    const int count = 200000;
    double state[2] = {1.2, -0.3};
    const double initialEnergy =
        0.5 * (state[1] * state[1] + omega * omega * state[0] * state[0]);
    double worstRelativeDrift = 0.0;

    for (int i = 0; i < count; ++i) {
        coalescent::rk4<2>(state, step, [&](const double* value, double* out) {
            out[0] = value[1];
            out[1] = -omega * omega * value[0];
        });
        const double energy =
            0.5 * (state[1] * state[1] + omega * omega * state[0] * state[0]);
        worstRelativeDrift = std::max(
            worstRelativeDrift, std::fabs(energy / initialEnergy - 1.0));
    }

    const bool ok = std::isfinite(state[0]) && std::isfinite(state[1])
                 && worstRelativeDrift < 2e-6;
    std::printf("RK4 oscillator max relative energy drift %.3e: %s\n",
                worstRelativeDrift, ok ? "PASS" : "FAIL");
    return ok;
}


int main() {
    const bool decay = checkDecayConvergence();
    const bool states = checkIndependentStates();
    const bool oscillator = checkOscillatorEnergy();
    if (!(decay && states && oscillator))
        return 1;
    std::printf("PASS: shared RK4 satisfies independent analytic contracts\n");
    return 0;
}
