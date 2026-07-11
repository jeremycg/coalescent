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
#include <cstring>
#include <limits>
#include <utility>
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
    float subTau = std::min(RATE_CAL * pitchHz / FS, HSUB * MAX_SUB);   // mirror production cap
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
    float_4 subTau = simd::fmin(RATE_CAL * pitchHz / FS, float_4(HSUB * MAX_SUB));   // mirror production cap
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

// ── Soma / Hindmarsh–Rose: the 3-state path (grouping + shared-K + subTau cap) ──
static constexpr float HR_A = 1.f, HR_B = 3.f, HR_C = 1.f, HR_D = 5.f, HR_XR = -1.6f;
static constexpr float HR_RATE_CAL = 55.364003f, HR_SMAX = 25.f;
template <typename T>
static inline void fHR(T x, T y, T z, T I, T r, T s, T& dx, T& dy, T& dz) {
    dx = y - HR_A*x*x*x + HR_B*x*x - z + I;
    dy = HR_C - HR_D*x*x - y;
    dz = r * (s * (x - HR_XR) - z);
}
static float stepScalarHR(float& x, float& y, float& z, float pitchHz, float I, float r, float s) {
    float subTau = std::min(HR_RATE_CAL * pitchHz / FS, HSUB * MAX_SUB);
    int K = std::min(MAX_SUB, std::max(MIN_SUB, (int) std::ceil(subTau / HSUB)));
    float h = subTau / K;
    float st[3] = {x, y, z};
    for (int k = 0; k < K; k++)
        coalescent::rk4<3>(st, h, [&](const float* Y, float* d) { fHR(Y[0], Y[1], Y[2], I, r, s, d[0], d[1], d[2]); });
    for (int i = 0; i < 3; i++) { if (!std::isfinite(st[i])) { st[0]=-1.6f; st[1]=0.f; st[2]=0.f; break; } }
    for (int i = 0; i < 3; i++) st[i] = std::clamp(st[i], -HR_SMAX, HR_SMAX);
    x = st[0]; y = st[1]; z = st[2];
    return x;
}
static float_4 stepSimdHR(float_4& x, float_4& y, float_4& z, float_4 pitchHz, float_4 I, float_4 r, float_4 s) {
    float_4 subTau = simd::fmin(HR_RATE_CAL * pitchHz / FS, float_4(HSUB * MAX_SUB));
    float_4 Kf = simd::ceil(subTau / HSUB);
    int K = MIN_SUB; for (int l = 0; l < 4; l++) K = std::max(K, (int) Kf[l]); K = std::min(K, MAX_SUB);
    float_4 h = subTau / (float) K;
    float_4 st[3] = {x, y, z};
    for (int k = 0; k < K; k++)
        coalescent::rk4<3>(st, h, [&](const float_4* Y, float_4* d) { fHR(Y[0], Y[1], Y[2], I, r, s, d[0], d[1], d[2]); });
    float_4 fin = (st[0]==st[0]) & (st[1]==st[1]) & (st[2]==st[2])
        & (simd::abs(st[0])<1e6f) & (simd::abs(st[1])<1e6f) & (simd::abs(st[2])<1e6f);
    st[0] = simd::ifelse(fin, st[0], -1.6f); st[1] = simd::ifelse(fin, st[1], 0.f); st[2] = simd::ifelse(fin, st[2], 0.f);
    for (int i = 0; i < 3; i++) st[i] = simd::clamp(st[i], -HR_SMAX, HR_SMAX);
    x = st[0]; y = st[1]; z = st[2];
    return x;
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

    // ── Soma / HR: the 3-state path. Same idea — scalar vs float_4 (four identical
    // lanes) must land on the same pitch — plus a mixed-lane group for the shared-K
    // masking path, which the FHN test above cannot exercise for HR. ──
    {
        struct HRV { const char* n; float hz, I, r, s; };
        HRV vs[] = {{"tonic C4", C4, 2.0f, 0.03f, 4.f}, {"burst", C4, 2.6f, 0.008f, 4.f},
                    {"C2 low", C4/4, 2.0f, 0.03f, 4.f}, {"C6 high", C4*4, 2.0f, 0.03f, 4.f}};
        printf("Soma/HR (3-state) scalar vs float_4:\n");
        for (auto& t : vs) {
            float sx=-1.6f, sy=0.f, sz=0.f; float_4 vx(-1.6f), vy(0.f), vz(0.f);
            std::vector<float> a(FS), b(FS); double mx = 0;
            for (int n = 0; n < (int) FS; n++) {
                a[n] = stepScalarHR(sx, sy, sz, t.hz, t.I, t.r, t.s);
                b[n] = stepSimdHR(vx, vy, vz, float_4(t.hz), float_4(t.I), float_4(t.r), float_4(t.s))[0];
                mx = std::max(mx, (double) std::fabs(a[n] - b[n]));
            }
            float fa = freqOf(a), fb = freqOf(b);
            float cents = (fa > 0 && fb > 0) ? 1200.f * std::log2(fb / fa) : 0.f;
            worstCents = std::max(worstCents, (double) std::fabs(cents));
            printf("  %-10s scalar %8.2f  simd %8.2f  %+.3f cents  max|Δx|=%.3e\n", t.n, fa, fb, cents, mx);
        }
        // Mixed-lane HR group: four different voicings share one K (the max).
        float lHz[4] = {C4/4, C4, C4*2, C4*4}; float lI[4] = {2.0f, 2.4f, 2.0f, 2.2f};
        float_4 vx(-1.6f), vy(0.f), vz(0.f);
        std::vector<std::vector<float>> lane(4, std::vector<float>(FS));
        for (int n = 0; n < (int) FS; n++) {
            float_4 out = stepSimdHR(vx, vy, vz, float_4(lHz[0],lHz[1],lHz[2],lHz[3]),
                                     float_4(lI[0],lI[1],lI[2],lI[3]), float_4(0.03f), float_4(4.f));
            for (int l = 0; l < 4; l++) lane[l][n] = out[l];
        }
        printf("Mixed-lane HR group (one shared K) vs per-voice scalar:\n");
        for (int l = 0; l < 4; l++) {
            float sx=-1.6f, sy=0.f, sz=0.f; std::vector<float> ref(FS);
            for (int n = 0; n < (int) FS; n++) ref[n] = stepScalarHR(sx, sy, sz, lHz[l], lI[l], 0.03f, 4.f);
            float fRef = freqOf(ref), fLane = freqOf(lane[l]);
            float cents = (fRef > 0 && fLane > 0) ? 1200.f * std::log2(fLane / fRef) : 0.f;
            bool finite = std::isfinite(lane[l].back());
            worstCents = std::max(worstCents, (double) std::fabs(cents));
            printf("  lane %d (%6.1f Hz, I=%.2f): scalar %8.2f  simd %8.2f  %+.3f cents  %s\n",
                   l, lHz[l], lI[l], fRef, fLane, cents, finite ? "" : "NON-FINITE!");
            if (!finite) { printf("FAIL: mixed-lane HR output not finite\n"); return 1; }
        }
    }

    // ── Haptik lattice state/clamp pass: scalar vs float_4, must be BIT-identical.
    // Unlike the RK4 oscillators (where SIMD rounding causes benign phase drift), this
    // pass is a plain elementwise clamp with no accumulation, so the vectorized version
    // must match the scalar loop exactly. The module writes the clamp as fmax(fmin(x,hi),lo)
    // to mirror rack::clamp's op order (which also fixes _mm_min/max_ps's NaN and signed-
    // zero behavior to match std::fmin/fmax here). Swept over every production N=8..128
    // (all four tail residues) with a value set that includes the ±16 boundary, values far
    // past it, ±inf, NaN, ±0, and subnormals, compared with memcmp for true bit-equality. ──
    {
        const float SMAX_H = 16.f;          // Haptik STATE_MAX
        const float INF = std::numeric_limits<float>::infinity();
        const float NAN_ = std::numeric_limits<float>::quiet_NaN();
        const float SUB = std::numeric_limits<float>::denorm_min();
        const float extras[] = { 0.f, -0.f, SMAX_H, -SMAX_H, 100.f, -100.f, INF, -INF, NAN_, SUB, -SUB };
        const int   NE = (int) (sizeof(extras) / sizeof(extras[0]));

        long totalDiffs = 0;
        for (int N = 8; N <= 128; N++) {                     // every N, so every tail residue 0..3
            std::vector<float> y0(N), v0(N);
            for (int i = 0; i < N; i++) {
                // seed the first entries with the pathological values, the rest with a
                // spread that straddles the clamp boundary in both directions.
                if (i < NE) { y0[i] = extras[i]; v0[i] = extras[(i * 7 + 3) % NE]; }
                else {
                    float t = -1.f + 2.f * (i / (float)(N - 1));
                    y0[i] = t * 20.f * (1.f + 0.3f * std::sin(i * 1.7f));
                    v0[i] = -t * 25.f * (1.f + 0.2f * std::cos(i * 2.3f));
                }
            }
            std::vector<float> ys = y0, vs = v0;             // scalar reference (mirrors module)
            for (int i = 0; i < N; i++) {
                ys[i] = rack::math::clamp(ys[i] + vs[i], -SMAX_H, SMAX_H);
                vs[i] = rack::math::clamp(vs[i], -SMAX_H, SMAX_H);
            }
            std::vector<float> yv = y0, vv = v0;             // vectorized (mirrors module)
            const float_4 lo(-SMAX_H), hi(SMAX_H);
            int i = 0;
            for (; i + 4 <= N; i += 4) {
                float_4 Y = float_4::load(&yv[i]);
                float_4 V = float_4::load(&vv[i]);
                Y = simd::fmax(simd::fmin(Y + V, hi), lo);
                V = simd::fmax(simd::fmin(V, hi), lo);
                Y.store(&yv[i]); V.store(&vv[i]);
            }
            for (; i < N; i++) {
                yv[i] = rack::math::clamp(yv[i] + vv[i], -SMAX_H, SMAX_H);
                vv[i] = rack::math::clamp(vv[i], -SMAX_H, SMAX_H);
            }
            // memcmp: catches any bit difference (incl. signed zero), skipping NaN slots
            // where a raw bit-compare is not meaningful but both must still be NaN.
            for (int j = 0; j < N; j++) {
                for (auto pr : { std::make_pair(&ys[j], &yv[j]), std::make_pair(&vs[j], &vv[j]) }) {
                    if (std::isnan(*pr.first) || std::isnan(*pr.second)) {
                        if (std::isnan(*pr.first) != std::isnan(*pr.second)) totalDiffs++;
                    } else if (std::memcmp(pr.first, pr.second, sizeof(float)) != 0) {
                        totalDiffs++;
                    }
                }
            }
        }
        printf("Haptik clamp pass scalar vs float_4 (N=8..128, all tail residues, "
               "incl. inf/NaN/-0/subnormal): %ld bit diffs\n", totalDiffs);
        if (totalDiffs) { printf("FAIL: Haptik vectorized clamp differs from scalar\n"); return 1; }
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
