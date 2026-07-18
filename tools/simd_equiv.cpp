// Proof that the float_4 (poly-SIMD) RK4 path matches the scalar per-voice path
// for Axon's FitzHugh–Nagumo kernel. Runs the same voicings through both the
// scalar and vectorised forms of the exact shared production core (four
// identical lanes), and reports the max deviation, plus the same for fastTanh.
//
// Build and run through `make check-simd` with RACK_DIR set to the Rack SDK.
#include <simd/functions.hpp>
using namespace rack;              // so `simd::` resolves as it does inside the plugin
#include "../src/dsp/neuron_models.hpp"
#include "../src/dsp/haptik_core.hpp"
#include "../src/tanh_approx.hpp"
#include <cstdio>
#include <cmath>
#include <cstring>
#include <limits>
#include <utility>
#include <algorithm>
#include <vector>

using rack::simd::float_4;
using AxonCore = coalescent::neuron::AxonCore;
using SomaCore = coalescent::neuron::SomaCore;
static constexpr float FS = 44100.f;
static constexpr float C4 = coalescent::neuron::PITCH_REFERENCE_HZ;

// One accepted output observation through the exact scalar production core.
static float stepScalar(float& v, float& w, float pitchHz, float I, float eps, float a) {
    const coalescent::neuron::ScalarSchedule schedule =
        coalescent::neuron::scalarSchedule<AxonCore>(pitchHz, FS);
    float s[AxonCore::STATE_COUNT] = {v, w};
    AxonCore::advanceObservation(s, schedule.h, schedule.substeps, I, eps, a);
    AxonCore::repair(s);
    v = s[0]; w = s[1];
    return v;
}
// Same, four lanes via the exact group-max-K production core.
static float_4 stepSimd(float_4& v, float_4& w, float_4 pitchHz, float_4 I, float_4 eps, float_4 a) {
    const float_4 subTau = coalescent::neuron::boundedSubTau<AxonCore>(
        pitchHz, FS, 1,
        [](float_4 left, float_4 right) { return simd::fmin(left, right); });
    const int K = coalescent::neuron::groupSubsteps<AxonCore>(
        subTau, 4,
        [](float_4 value) { return simd::ceil(value); },
        [](float_4 value, int lane) { return value[lane]; });
    float_4 h = subTau / (float) K;
    float_4 s[AxonCore::STATE_COUNT] = {v, w};
    AxonCore::advanceObservation(s, h, K, I, eps, a);
    AxonCore::repair(
        s,
        [](float_4 value) { return simd::abs(value); },
        [](float_4 mask, float_4 yes, float_4 no) { return simd::ifelse(mask, yes, no); },
        [](float_4 value, float low, float high) { return simd::clamp(value, low, high); });
    v = s[0]; w = s[1];
    return v;
}

// ── Soma / Hindmarsh–Rose: the 3-state path (grouping + shared-K + subTau cap) ──
static float stepScalarHR(float& x, float& y, float& z, float pitchHz, float I, float r, float s) {
    const coalescent::neuron::ScalarSchedule schedule =
        coalescent::neuron::scalarSchedule<SomaCore>(pitchHz, FS);
    float st[SomaCore::STATE_COUNT] = {x, y, z};
    SomaCore::advanceObservation(st, schedule.h, schedule.substeps, I, r, s);
    SomaCore::repair(st);
    x = st[0]; y = st[1]; z = st[2];
    return x;
}
static float_4 stepSimdHR(float_4& x, float_4& y, float_4& z, float_4 pitchHz, float_4 I, float_4 r, float_4 s) {
    const float_4 subTau = coalescent::neuron::boundedSubTau<SomaCore>(
        pitchHz, FS, 1,
        [](float_4 left, float_4 right) { return simd::fmin(left, right); });
    const int K = coalescent::neuron::groupSubsteps<SomaCore>(
        subTau, 4,
        [](float_4 value) { return simd::ceil(value); },
        [](float_4 value, int lane) { return value[lane]; });
    float_4 h = subTau / (float) K;
    float_4 st[SomaCore::STATE_COUNT] = {x, y, z};
    SomaCore::advanceObservation(st, h, K, I, r, s);
    SomaCore::repair(
        st,
        [](float_4 value) { return simd::abs(value); },
        [](float_4 mask, float_4 yes, float_4 no) { return simd::ifelse(mask, yes, no); },
        [](float_4 value, float low, float high) { return simd::clamp(value, low, high); });
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
        float sv = AxonCore::REST_V, sw = AxonCore::REST_W;
        float_4 vv(AxonCore::REST_V), ww(AxonCore::REST_W);
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
        float_4 vv(AxonCore::REST_V), ww(AxonCore::REST_W);
        std::vector<std::vector<float>> lane(4, std::vector<float>(FS));
        for (int n = 0; n < (int) FS; n++) {
            float_4 out = stepSimd(vv, ww, hz, Iv, float_4(0.08f), float_4(0.7f));
            for (int l = 0; l < 4; l++) lane[l][n] = out[l];
        }
        printf("Mixed-lane group (one shared K) vs per-voice scalar:\n");
        for (int l = 0; l < 4; l++) {
            float sv = AxonCore::REST_V, sw = AxonCore::REST_W;
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
            float sx = SomaCore::REST_X, sy = SomaCore::REST_Y, sz = SomaCore::REST_Z;
            float_4 vx(SomaCore::REST_X), vy(SomaCore::REST_Y), vz(SomaCore::REST_Z);
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
        float_4 vx(SomaCore::REST_X), vy(SomaCore::REST_Y), vz(SomaCore::REST_Z);
        std::vector<std::vector<float>> lane(4, std::vector<float>(FS));
        for (int n = 0; n < (int) FS; n++) {
            float_4 out = stepSimdHR(vx, vy, vz, float_4(lHz[0],lHz[1],lHz[2],lHz[3]),
                                     float_4(lI[0],lI[1],lI[2],lI[3]), float_4(0.03f), float_4(4.f));
            for (int l = 0; l < 4; l++) lane[l][n] = out[l];
        }
        printf("Mixed-lane HR group (one shared K) vs per-voice scalar:\n");
        for (int l = 0; l < 4; l++) {
            float sx = SomaCore::REST_X, sy = SomaCore::REST_Y, sz = SomaCore::REST_Z;
            std::vector<float> ref(FS);
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

    // The production neuron wrappers use this exact SIMD policy for CURRENT
    // cable lanes. It must preserve finite lanes and neutralize hostile lanes
    // even under this tool's -funsafe-math-optimizations build.
    {
        const float floatMax = std::numeric_limits<float>::max();
        const float inf = std::numeric_limits<float>::infinity();
        const float nan = std::numeric_limits<float>::quiet_NaN();
        const float_4 raw(2.5f, nan, inf, -inf);
        const float_4 guarded = coalescent::finiteOr(
            raw, float_4(0.f),
            [floatMax](float_4 value) {
                return (value >= -floatMax) & (value <= floatMax);
            },
            [](float_4 mask, float_4 yes, float_4 no) {
                return simd::ifelse(mask, yes, no);
            });
        const bool finiteGuardOk = guarded[0] == 2.5f && guarded[1] == 0.f
            && guarded[2] == 0.f && guarded[3] == 0.f;
        printf("Neuron CURRENT SIMD finite-lane guard: %s\n",
               finiteGuardOk ? "PASS" : "FAIL");
        if (!finiteGuardOk)
            return 1;
    }

    {
        const float amount = 0.6f;
        const float finiteVoltage = -7.25f;
        const float expected = finiteVoltage * coalescent::haptik::kExternalGain * amount;
        const float finiteDrive =
            coalescent::haptik::externalDrive(finiteVoltage, amount);
        const float nanDrive = coalescent::haptik::externalDrive(
            std::numeric_limits<float>::quiet_NaN(), amount);
        const float positiveInfDrive = coalescent::haptik::externalDrive(
            std::numeric_limits<float>::infinity(), amount);
        const float negativeInfDrive = coalescent::haptik::externalDrive(
            -std::numeric_limits<float>::infinity(), amount);
        const bool externalGuardOk =
            std::fabs(finiteDrive - expected) < 1e-8f
            && nanDrive == 0.f && positiveInfDrive == 0.f && negativeInfDrive == 0.f;
        printf("Haptik EXT scalar finite guard under unsafe-math "
               "(finite=%.8f nan=%g +inf=%g -inf=%g): %s\n",
               finiteDrive, nanDrive, positiveInfDrive, negativeInfDrive,
               externalGuardOk ? "PASS" : "FAIL");
        if (!externalGuardOk)
            return 1;
    }

    // ── Haptik lattice state/clamp pass: scalar vs float_4, must be BIT-identical.
    // Unlike the RK4 oscillators (where SIMD rounding causes benign phase drift), this
    // pass is a plain elementwise clamp with no accumulation, so the vectorized version
    // must match the scalar loop exactly. Both paths call the shared production
    // advanceStateValue(), whose injected min/max operations preserve fmax(fmin(x,hi),lo)
    // ordering (which also fixes _mm_min/max_ps's NaN and signed-
    // zero behavior to match std::fmin/fmax here). Swept over every production N=8..128
    // (all four tail residues) with a value set that includes the ±16 boundary, values far
    // past it, ±inf, NaN, ±0, and subnormals, compared with memcmp for true bit-equality. ──
    {
        const float SMAX_H = coalescent::haptik::kStateMax;
        const float INF = std::numeric_limits<float>::infinity();
        const float NAN_ = std::numeric_limits<float>::quiet_NaN();
        const float SUB = std::numeric_limits<float>::denorm_min();
        const float extras[] = { 0.f, -0.f, SMAX_H, -SMAX_H, 100.f, -100.f, INF, -INF, NAN_, SUB, -SUB };
        const int   NE = (int) (sizeof(extras) / sizeof(extras[0]));

        long totalDiffs = 0;
        for (int N = coalescent::haptik::kMinPlayableMasses;
             N <= coalescent::haptik::kMaxPlayableMasses; N++) {
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
            std::vector<float> ys = y0, vs = v0;
            for (int i = 0; i < N; i++)
                coalescent::haptik::advanceStateValue(ys[i], vs[i]);
            std::vector<float> yv = y0, vv = v0;
            const float_4 lo(-SMAX_H), hi(SMAX_H);
            int i = 0;
            for (; i + 4 <= N; i += 4) {
                float_4 Y = float_4::load(&yv[i]);
                float_4 V = float_4::load(&vv[i]);
                coalescent::haptik::advanceStateValue(
                    Y, V, lo, hi,
                    [](const float_4& a, const float_4& b) {
                        return simd::fmin(a, b);
                    },
                    [](const float_4& a, const float_4& b) {
                        return simd::fmax(a, b);
                    });
                Y.store(&yv[i]); V.store(&vv[i]);
            }
            for (; i < N; i++)
                coalescent::haptik::advanceStateValue(yv[i], vv[i]);
            // memcmp: catches any bit difference (incl. signed zero), skipping NaN slots
            // where a raw bit-compare is not meaningful but both must still be NaN.
            for (int j = 0; j < N; j++) {
                const std::pair<const float*, const float*> comparisons[] = {
                    {&ys[j], &yv[j]}, {&vs[j], &vv[j]}
                };
                for (const auto& pr : comparisons) {
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
