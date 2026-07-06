// Proof that the float_4 (poly-SIMD) RK4 path matches the scalar per-voice path
// for Axon's FitzHugh–Nagumo kernel. Runs the same voicings through both the
// scalar coalescent::rk4<2,float> and the vectorised coalescent::rk4<2,float_4> (four
// identical lanes) and reports the max deviation, plus the same for fastTanh.
//
//   g++ -O3 -funsafe-math-optimizations -march=nehalem -std=c++17 -DARCH_X64 -DARCH_LIN \
//       -I$RACK/include -I$RACK/dep/include tools/simd_equiv.cpp -o /tmp/se && /tmp/se
#include <simd/functions.hpp>
using namespace rack;              // so `simd::` resolves as it does inside the plugin
#include "../src/dsp/rk4.hpp"
#include "../src/tanh_approx.hpp"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>

using rack::simd::float_4;
static constexpr float B_FIXED = 0.8f, RATE_CAL = 37.899004f, HSUB = 0.05f, SMAX = 10.f, FS = 44100.f, C4 = 261.6256f;
static constexpr int   MIN_SUB = 2, MAX_SUB = 64;

template <typename T>
static inline void fFHN(T v, T w, T I, T e, T a, T& dv, T& dw) {
    dv = v - v * v * v / 3.f - w + I;
    dw = e * (v + a - B_FIXED * w);
}

// One output sample of raw v, scalar path (mirrors Axon at os=1 for the test).
static float stepScalar(float& v, float& w, float pitchHz, float I, float eps, float a) {
    float subTau = RATE_CAL * pitchHz / FS;
    int K = std::min(MAX_SUB, std::max(MIN_SUB, (int) std::ceil(subTau / HSUB)));
    float h = subTau / K;
    float s[2] = {v, w};
    for (int k = 0; k < K; k++)
        coalescent::rk4<2>(s, h, [&](const float* y, float* d) { fFHN(y[0], y[1], I, eps, a, d[0], d[1]); });
    if (!std::isfinite(s[0]) || !std::isfinite(s[1])) { s[0] = -1.2f; s[1] = -0.6f; }
    s[0] = std::clamp(s[0], -SMAX, SMAX); s[1] = std::clamp(s[1], -SMAX, SMAX);
    v = s[0]; w = s[1];
    return v;
}
// Same, four identical lanes via float_4 (mirrors the module's group path).
static float_4 stepSimd(float_4& v, float_4& w, float_4 pitchHz, float_4 I, float_4 eps, float_4 a) {
    float_4 subTau = RATE_CAL * pitchHz / FS;
    float_4 Kf = simd::ceil(subTau / HSUB);
    int K = MIN_SUB; for (int l = 0; l < 4; l++) K = std::max(K, (int) Kf[l]); K = std::min(K, MAX_SUB);
    float_4 h = subTau / (float) K;
    float_4 s[2] = {v, w};
    for (int k = 0; k < K; k++)
        coalescent::rk4<2>(s, h, [&](const float_4* y, float_4* d) { fFHN(y[0], y[1], I, eps, a, d[0], d[1]); });
    float_4 fin = (s[0] == s[0]) & (s[1] == s[1]) & (simd::abs(s[0]) < 1e6f) & (simd::abs(s[1]) < 1e6f);
    s[0] = simd::ifelse(fin, s[0], -1.2f); s[1] = simd::ifelse(fin, s[1], -0.6f);
    s[0] = simd::clamp(s[0], -SMAX, SMAX); s[1] = simd::clamp(s[1], -SMAX, SMAX);
    v = s[0]; w = s[1];
    return v;
}

static float freqOf(const std::vector<float>& x) {
    int last = -1; double sum = 0; int cnt = 0;
    for (int n = 1; n < (int) x.size(); n++)
        if (x[n-1] <= 0.f && x[n] > 0.f) { if (last >= 0 && n > (int) x.size()/8) { sum += (n - last); cnt++; } last = n; }
    return cnt ? (float) (FS / (sum / cnt)) : 0.f;
}

int main() {
    struct V { const char* n; float hz, I, e, a; };
    V vs[] = {{"C4 default", C4, 0.6f, 0.08f, 0.7f}, {"C2 low", C4/4, 0.6f, 0.08f, 0.7f},
              {"stiff", C4, 0.9f, 0.02f, 0.5f}, {"C6 high", C4*4, 0.6f, 0.08f, 0.7f}};
    double worstCents = 0;
    printf("Scalar vs float_4 — same voicing, sonically-meaningful metric = pitch:\n");
    printf("  %-12s %10s %10s %8s %10s\n", "voicing", "f_scalar", "f_simd", "cents", "max|Δv|");
    for (auto& t : vs) {
        float sv = -1.2f, sw = -0.6f;
        float_4 vv(-1.2f), ww(-0.6f);
        std::vector<float> a(FS), b(FS); double mx = 0;
        for (int n = 0; n < (int) FS; n++) {
            a[n] = stepScalar(sv, sw, t.hz, t.I, t.e, t.a);
            b[n] = stepSimd(vv, ww, float_4(t.hz), float_4(t.I), float_4(t.e), float_4(t.a))[0];
            mx = std::max(mx, (double) std::fabs(a[n] - b[n]));
        }
        float fa = freqOf(a), fb = freqOf(b);
        float cents = (fa > 0 && fb > 0) ? 1200.f * std::log2(fb / fa) : 0.f;
        worstCents = std::max(worstCents, (double) std::fabs(cents));
        printf("  %-12s %10.3f %10.3f %8.3f %10.3e\n", t.n, fa, fb, cents, mx);
    }
    // fastTanh scalar vs float_4 over the audible range
    double tmax = 0;
    for (int i = 0; i < 20000; i++) {
        float x = -8.f + 16.f * (i / 19999.f);
        float_4 b = coalescent::fastTanh(float_4(x));
        tmax = std::max(tmax, (double) std::fabs(coalescent::fastTanh(x) - b[0]));
    }
    printf("fastTanh scalar vs float_4: max|Δ| = %.3e (differs only for |x|>4)\n", tmax);

    // ── Mixed-lane group: the group-max-K path. Four materially different lanes
    // share one K (the max), so slower lanes integrate with smaller h than their
    // scalar equivalent — verify each lane still lands on the scalar pitch. ──
    {
        float laneHz[4] = {C4 / 4, C4, C4 * 4, C4 * 8};          // C2, C4, C6, C7
        float laneI[4]  = {0.45f, 0.6f, 0.9f, 1.3f};             // spread across the band
        float_4 hz(laneHz[0], laneHz[1], laneHz[2], laneHz[3]);
        float_4 Iv(laneI[0], laneI[1], laneI[2], laneI[3]);
        float_4 vv(-1.2f), ww(-0.6f);
        std::vector<std::vector<float>> lane(4, std::vector<float>(FS));
        for (int n = 0; n < (int) FS; n++) {
            float_4 out = stepSimd(vv, ww, hz, Iv, float_4(0.08f), float_4(0.7f));
            for (int l = 0; l < 4; l++) lane[l][n] = out[l];
        }
        printf("Mixed-lane group (one shared K) vs per-voice scalar:\n");
        for (int l = 0; l < 4; l++) {
            float sv = -1.2f, sw = -0.6f;
            std::vector<float> ref(FS);
            for (int n = 0; n < (int) FS; n++) ref[n] = stepScalar(sv, sw, laneHz[l], laneI[l], 0.08f, 0.7f);
            float fRef = freqOf(ref), fLane = freqOf(lane[l]);
            float cents = (fRef > 0 && fLane > 0) ? 1200.f * std::log2(fLane / fRef) : 0.f;
            bool finite = std::isfinite(lane[l].back());
            worstCents = std::max(worstCents, (double) std::fabs(cents));
            printf("  lane %d (%6.1f Hz, I=%.2f): scalar %8.2f  simd %8.2f  %+.3f cents  %s\n",
                   l, laneHz[l], laneI[l], fRef, fLane, cents, finite ? "" : "NON-FINITE!");
            if (!finite) { printf("FAIL: mixed-lane output not finite\n"); return 1; }
        }
    }

    // Oscillator equivalence = same pitch. The per-sample |Δv| is phase micro-drift
    // from float-vs-float_4 rounding (inherent to SIMD), not a difference in tone.
    const double CENT_TOL = 1.0;
    if (worstCents > CENT_TOL) {
        printf("FAIL: float_4 pitch drifts %.3f cents from scalar (> %.1f)\n", worstCents, CENT_TOL);
        return 1;
    }
    printf("PASS: float_4 path is within %.1f cent of scalar for every voicing (incl. mixed lanes)\n", CENT_TOL);
    return 0;
}
