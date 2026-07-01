// Regression test + calibration for Axon's FitzHugh-Nagumo kernel.
//
// Replicates the src/Axon.cpp integrator (factored f() + RK4 substepping in
// dimensionless time) in plain C++ so the safety invariant, pitch calibration,
// and stiffness behaviour can be checked without launching Rack. Keep this in
// sync with the kernel if the math changes.
//
//   Build & run:  g++ -O2 -o /tmp/axon_test stability_test.cpp && /tmp/axon_test
//   Exit 0 = all checks pass, 1 = a check failed.
//
// What it does (plan §5):
//   1. CALIBRATION: measure the dimensionless limit-cycle period T_dim at the
//      default params. RATE_CAL must equal T_dim so that, with
//      dtau = RATE_CAL * pitchHz / fs, the emergent audio fundamental comes out
//      at pitchHz (i.e. C4 at 0 V). Prints the value to bake into Axon.cpp.
//   2. STABILITY: sweep CURRENT × EPS × SHAPE, long run at several pitches,
//      assert v,w stay finite and bounded throughout.
//   3. PITCH: with the measured RATE_CAL, drive the full per-sample loop and
//      confirm the output fundamental tracks V/OCT within tolerance.
#include <cstdio>
#include <cmath>
#include <algorithm>

// ── kernel constants (mirror src/Axon.cpp) ──
static constexpr float B_FIXED   = 0.8f;
static constexpr float HSUB_MAX  = 0.05f;
static constexpr int   MIN_SUB   = 2;
static constexpr int   MAX_SUB   = 64;
static constexpr float STATE_MAX = 10.f;
static constexpr float FREQ_C4   = 261.6256f;   // rack dsp::FREQ_C4

static inline int clampi(int x, int lo, int hi) { return x < lo ? lo : (x > hi ? hi : x); }
static inline float clampf(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }

// FHN derivatives in dimensionless time (the factored f(); HR will add a 3rd line).
static inline void fhn(float v, float w, float Itot, float eps, float a,
                       float& dv, float& dw) {
    dv = v - v * v * v / 3.f - w + Itot;
    dw = eps * (v + a - B_FIXED * w);
}

// One RK4 substep of size h. Updates v,w in place.
static inline void rk4(float& v, float& w, float h, float Itot, float eps, float a) {
    float k1v, k1w, k2v, k2w, k3v, k3w, k4v, k4w;
    fhn(v,                 w,                 Itot, eps, a, k1v, k1w);
    fhn(v + 0.5f * h * k1v, w + 0.5f * h * k1w, Itot, eps, a, k2v, k2w);
    fhn(v + 0.5f * h * k2v, w + 0.5f * h * k2w, Itot, eps, a, k3v, k3w);
    fhn(v + h * k3v,        w + h * k3w,        Itot, eps, a, k4v, k4w);
    v += h / 6.f * (k1v + 2.f * k2v + 2.f * k3v + k4v);
    w += h / 6.f * (k1w + 2.f * k2w + 2.f * k3w + k4w);
}

// Advance dimensionless time by dtau using adaptive substepping; returns max|state|.
// Mirrors the per-sample integrate block in src/Axon.cpp.
static inline void advance(float& v, float& w, float dtau, float Itot, float eps, float a) {
    int K = clampi((int)std::ceil(dtau / HSUB_MAX), MIN_SUB, MAX_SUB);
    float h = dtau / K;
    for (int k = 0; k < K; k++) rk4(v, w, h, Itot, eps, a);
    if (!std::isfinite(v) || !std::isfinite(w)) { v = -1.2f; w = -0.6f; }
    v = clampf(v, -STATE_MAX, STATE_MAX);
    w = clampf(w, -STATE_MAX, STATE_MAX);
}

// Measure the dimensionless limit-cycle period by advancing tau in small fixed
// steps and timing upward zero(threshold)-crossings of v. Returns period in tau
// units, or -1 if it never oscillates (quiescent).
static double measurePeriodTau(float I, float eps, float a) {
    float v = -1.2f, w = -0.6f;
    const float thr = 1.0f;              // SPIKE_THRESH
    const float dstep = 1e-4f;           // fine tau step for accurate period
    // settle onto the limit cycle first
    double tau = 0.0;
    for (long i = 0; i < (long)(400.0 / dstep); i++) { advance(v, w, dstep, I, eps, a); tau += dstep; }
    // now time crossings
    float vprev = v;
    double firstCross = -1.0, lastCross = -1.0;
    int crossings = 0;
    long maxIter = (long)(2000.0 / dstep);
    for (long i = 0; i < maxIter && crossings < 21; i++) {
        advance(v, w, dstep, I, eps, a);
        tau += dstep;
        if (vprev < thr && v >= thr) {   // upward crossing
            if (firstCross < 0) firstCross = tau;
            lastCross = tau;
            crossings++;
        }
        vprev = v;
    }
    if (crossings < 2) return -1.0;
    return (lastCross - firstCross) / (crossings - 1);
}

int main() {
    int failures = 0;

    // ── 1. CALIBRATION at default params (I=0.6, eps=0.08, a=0.7) ──
    double Tdim = measurePeriodTau(0.6f, 0.08f, 0.7f);
    if (Tdim <= 0.0) {
        printf("FAIL: default params do not oscillate (T_dim=%.4f)\n", Tdim);
        failures++;
    } else {
        printf("CALIBRATION: dimensionless limit-cycle period T_dim = %.6f\n", Tdim);
        printf("             => set RATE_CAL = %.6ff   (gives C4 at 0 V)\n", Tdim);
    }
    float RATE_CAL = (Tdim > 0) ? (float)Tdim : 6.0f;

    // ── 2. STABILITY sweep ──
    const float fs = 48000.f;
    int sweepCount = 0;
    for (float I = -0.2f; I <= 1.6001f; I += 0.2f)
        for (float eps = 0.01f; eps <= 0.3001f; eps += 0.04f)
            for (float a = 0.4f; a <= 1.0001f; a += 0.15f)
                for (float oct = -4.f; oct <= 4.001f; oct += 4.f) {
                    float v = -1.2f, w = -0.6f;
                    float pitchHz = FREQ_C4 * std::exp2(oct);
                    float dtau = RATE_CAL * pitchHz / fs;
                    long n = (long)(fs * 0.5f);   // 0.5 s
                    for (long i = 0; i < n; i++) {
                        advance(v, w, dtau, I, eps, a);
                        if (!std::isfinite(v) || !std::isfinite(w) ||
                            std::fabs(v) > STATE_MAX + 1e-3f || std::fabs(w) > STATE_MAX + 1e-3f) {
                            printf("FAIL: state escaped at I=%.2f eps=%.3f a=%.2f oct=%.0f: v=%g w=%g\n",
                                   I, eps, a, oct, v, w);
                            failures++;
                            break;
                        }
                    }
                    sweepCount++;
                }
    printf("STABILITY: %d parameter/pitch combinations stayed finite & bounded\n", sweepCount);

    // ── 3. PITCH tracking: measure output period in samples at several octaves ──
    auto measureAudioHz = [&](float oct, float I, float eps, float a) -> double {
        float v = -1.2f, w = -0.6f;
        float pitchHz = FREQ_C4 * std::exp2(oct);
        float dtau = RATE_CAL * pitchHz / fs;
        const float thr = 1.0f;
        // settle
        for (long i = 0; i < (long)(fs * 0.3f); i++) advance(v, w, dtau, I, eps, a);
        float vprev = v;
        double firstCross = -1, lastCross = -1; int crossings = 0; long s = 0;
        long maxS = (long)(fs * 2.0f);
        for (; s < maxS && crossings < 21; s++) {
            advance(v, w, dtau, I, eps, a);
            if (vprev < thr && v >= thr) { if (firstCross < 0) firstCross = s; lastCross = s; crossings++; }
            vprev = v;
        }
        if (crossings < 2) return -1;
        double periodSamples = (lastCross - firstCross) / (crossings - 1);
        return fs / periodSamples;
    };
    for (float oct = -2.f; oct <= 2.001f; oct += 1.f) {
        double want = FREQ_C4 * std::exp2(oct);
        double got = measureAudioHz(oct, 0.6f, 0.08f, 0.7f);
        double cents = (got > 0) ? 1200.0 * std::log2(got / want) : 0;
        printf("PITCH: oct=%+.0f  want=%7.2f Hz  got=%7.2f Hz  (%+6.1f cents)\n", oct, want, got, cents);
        if (got < 0 || std::fabs(cents) > 50.0) {   // ½ semitone tolerance (open-loop)
            printf("  FAIL: pitch off by more than 50 cents\n");
            failures++;
        }
    }

    if (failures) { printf("\n%d CHECK(S) FAILED\n", failures); return 1; }
    printf("\nAll checks passed.\n");
    return 0;
}
