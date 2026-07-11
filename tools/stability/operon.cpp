// Standalone stability + calibration replica of Operon's repressilator kernel.
// Mirrors src/Operon.cpp process() math (dimensionless Elowitz–Leibler, RK4 with
// pitch-adaptive substepping). Sweeps ALPHA × HILL × BETA × LEAK at several
// pitches over long runs and asserts all six states stay finite and bounded;
// also reports the default-voicing dimensionless period (→ RATE_CAL) and where
// the ALPHA dial crosses the Hopf threshold. No Rack SDK needed.
//
//   g++ -O2 -o /tmp/t tools/stability/operon.cpp && /tmp/t     (exit 0 = pass)
#include <cstdio>
#include <cmath>
#include <algorithm>

static constexpr float HSUB_MAX = 0.05f, STATE_MAX = 1e3f, BIAS_EPS = 1e-6f;
static constexpr int   MIN_SUB = 2, MAX_SUB = 64;
static constexpr int   CENTER_ITERS = 16, NEWTON_ITERS = 4;   // mirror src/Operon.cpp

// pStar solvers. `bisect(iters)` is the reference; the module uses CENTER_ITERS as
// its fallback but a 40-iter bisection here gives a near-exact root to compare against.
static float pStarBisect(float alpha, float n, float a0, int iters) {
    float lo = 0.f, hi = std::max(1.f, alpha + a0 + 1.f);
    for (int i = 0; i < iters; ++i) { float mid = 0.5f * (lo + hi);
        float g = mid - alpha / (1.f + std::pow(std::max(mid, 0.f), n)) - a0; if (g > 0.f) hi = mid; else lo = mid; }
    return 0.5f * (lo + hi);
}
// Warm-started Newton — mirrors src/Operon.cpp exactly (falls back to bisection on a
// bracket exit or non-convergence within NEWTON_ITERS).
static float pStarNewton(float alpha, float n, float a0, float guess) {
    const float hi = std::max(1.f, alpha + a0 + 1.f);
    if (guess > 0.f && guess < hi) { float x = guess;
        for (int i = 0; i < NEWTON_ITERS; ++i) { float xn = std::pow(std::max(x, 1e-6f), n); float d = 1.f + xn;
            float g = x - alpha / d - a0; float gp = 1.f + alpha * n * (xn / std::max(x, 1e-6f)) / (d * d);
            float step = g / gp; x -= step; if (!(x > 0.f && x < hi)) break; if (std::fabs(step) < 1e-5f) return x; }
    }
    return pStarBisect(alpha, n, a0, CENTER_ITERS);
}

struct Par { float alpha, n, beta, a0; };

static inline void deriv(const float Y[6], float D[6], const Par& p, float perturb) {
    for (int i = 0; i < 3; ++i) {
        float rep = std::max(Y[3 + ((i + 2) % 3)], 0.f);
        D[i]     = -Y[i] + p.alpha / (1.f + std::pow(rep, p.n)) + p.a0;
        D[3 + i] = -p.beta * (Y[3 + i] - Y[i]);
    }
    D[0] += perturb;
    D[0] +=  BIAS_EPS; D[1] += -0.5f * BIAS_EPS; D[2] += -0.5f * BIAS_EPS;
}
static inline void rk4(float y[6], float h, const Par& p, float perturb) {
    float k1[6], k2[6], k3[6], k4[6], t[6];
    deriv(y, k1, p, perturb); for (int i = 0; i < 6; ++i) t[i] = y[i] + 0.5f * h * k1[i];
    deriv(t, k2, p, perturb); for (int i = 0; i < 6; ++i) t[i] = y[i] + 0.5f * h * k2[i];
    deriv(t, k3, p, perturb); for (int i = 0; i < 6; ++i) t[i] = y[i] + h * k3[i];
    deriv(t, k4, p, perturb);
    for (int i = 0; i < 6; ++i) y[i] += h / 6.f * (k1[i] + 2.f * k2[i] + 2.f * k3[i] + k4[i]);
}

static inline int subK(float dtau) {
    return std::min(MAX_SUB, std::max(MIN_SUB, (int) std::ceil(std::min(dtau, HSUB_MAX * MAX_SUB) / HSUB_MAX)));
}
// Fast finite/bounded check — single pass, no measurement.
static bool bounded(const Par& p, float dtau, int NS) {
    float y[6] = {0.2f, 0.1f, 0.3f, 0.1f, 0.4f, 0.15f};
    int K = subK(dtau); float h = std::min(dtau, HSUB_MAX * MAX_SUB) / K;
    for (int s = 0; s < NS; ++s) {
        for (int k = 0; k < K; ++k) rk4(y, h, p, 0.f);
        for (int j = 0; j < 6; ++j)
            if (!std::isfinite(y[j]) || std::fabs(y[j]) > STATE_MAX) return false;
    }
    return true;
}
// Centered-p0 amplitude + dimensionless period (for the two calibration voicings).
static void measure(const Par& p, float dtau, float& amp, float& periodTau) {
    float y[6] = {0.2f, 0.1f, 0.3f, 0.1f, 0.4f, 0.15f};
    int K = subK(dtau); float h = std::min(dtau, HSUB_MAX * MAX_SUB) / K;
    const int NS = 60000; double mean = 0; int cnt = 0;
    for (int s = 0; s < NS; ++s) { for (int k = 0; k < K; ++k) rk4(y, h, p, 0.f); if (s > NS / 2) { mean += y[3]; cnt++; } }
    mean /= std::max(cnt, 1);
    float mn = 1e9f, mx = -1e9f, prev = 0.f, lastCross = -1.f; double sumP = 0; int nP = 0;
    float y2[6] = {0.2f, 0.1f, 0.3f, 0.1f, 0.4f, 0.15f};
    for (int s = 0; s < NS; ++s) {
        for (int k = 0; k < K; ++k) rk4(y2, h, p, 0.f);
        if (s > NS / 2) {
            float v = y2[3]; mn = std::min(mn, v); mx = std::max(mx, v);
            float c = v - (float) mean; float tcur = s * K * h;
            if (prev <= 0.f && c > 0.f) { if (lastCross > 0) { sumP += tcur - lastCross; nP++; } lastCross = tcur; }
            prev = c;
        }
    }
    amp = mx - mn; periodTau = nP ? (float)(sumP / nP) : 0.f;
}

int main() {
    int fails = 0, runs = 0;
    // Bounded across the parameter box at three pitch regimes (incl. above the cap).
    for (float alpha : {0.f, 12.f, 40.f, 80.f})
    for (float n : {1.2f, 2.5f, 8.f})
    for (float beta : {0.2f, 1.f, 8.f})
    for (float a0 : {0.f, 0.5f})
    for (float dtau : {0.02f, 0.4f, 3.2f}) {
        Par p{alpha, n, beta, a0}; runs++;
        int NS = dtau > 1.f ? 8000 : 20000;   // high-K case runs fewer samples
        if (!bounded(p, dtau, NS)) { fails++; printf("  UNBOUNDED: a=%.0f n=%.1f b=%.1f a0=%.2f dtau=%.3f\n", alpha, n, beta, a0, dtau); }
    }
    printf("stability: %d/%d runs finite+bounded\n", runs - fails, runs);

    // Calibration: default voicing period (→ RATE_CAL) and Hopf crossing on ALPHA.
    float amp, per; Par def{12.f, 2.5f, 1.f, 0.05f};
    measure(def, 0.02f, amp, per);
    printf("default voicing (a=12 n=2.5 b=1 a0=0.05): amp=%.2f  period=%.2f tau  => RATE_CAL~%.2f\n", amp, per, per);
    printf("ALPHA Hopf crossing (n=2.5 b=1 a0=0.05): ");
    for (float a = 1.f; a <= 60.f; a += 1.f) { Par p{a, 2.5f, 1.f, 0.05f}; measure(p, 0.02f, amp, per); if (amp > 0.05f) { printf("oscillates from alpha~%.0f\n", a); break; } }

    // Hill LUT vs direct pow: max abs error of 1/(1+rep^n) over rep×n (guards the
    // module's 8192-point lookup, incl. the sharp n=8 knee).
    const int M = 8192; const float X = 64.f; static float lut[M + 1];
    double lutMax = 0;
    for (float nn : {1.2f, 2.5f, 5.f, 8.f}) {
        for (int i = 0; i <= M; ++i) lut[i] = 1.f / (1.f + std::pow(X * i / M, nn));
        for (float rep = 0.f; rep < X; rep += 0.013f) {
            float f = rep * (M / X); int i = (int) f;
            float lu = lut[i] + (f - i) * (lut[i + 1] - lut[i]);
            float ex = 1.f / (1.f + std::pow(rep, nn));
            lutMax = std::max(lutMax, (double) std::fabs(lu - ex));
        }
    }
    printf("Hill LUT max |error| vs pow (n up to 8): %.2e\n", lutMax);

    // Warm-started Newton must converge to the accurate root when seeded from a
    // nearby root (the per-sample-modulation case the module hits). Seed from the
    // root at a 1%-different ALPHA and require convergence to the 40-iter reference.
    double pStarMax = 0;
    for (float a0 : {0.f, 0.05f, 0.5f, 1.5f})
        for (float alpha = 0.5f; alpha <= 80.f; alpha *= 1.05f)
            for (float nn = 1.01f; nn <= 10.f; nn += 0.13f) {
                float ref    = pStarBisect(alpha, nn, a0, 40);
                float warm   = pStarBisect(alpha * 0.99f, nn, a0, 40);   // "previous sample" root
                float got    = pStarNewton(alpha, nn, a0, warm);
                pStarMax = std::max(pStarMax, (double) std::fabs(got - ref));
            }
    printf("pStar warm Newton vs accurate root: max |Δ| = %.2e\n", pStarMax);

    // Hill-LUT gate state machine (mirror of src/Operon.cpp): a build starts only
    // after n has been within N_MOVE_EPS of an anchor for 2048 consecutive samples;
    // any move — fast OR slow cumulative drift — resets the window and drops the
    // partial build. Assert (a) a slow sub-threshold HILL LFO does NOT thrash the
    // 256-entry fill (the regression the incremental rebuild could otherwise cause,
    // aborting+restarting every sample), and (b) a genuinely static n still completes
    // a build so the LUT goes live.
    {
        const int   LN = 8192, SLICE = 256;
        const float EPS = 1e-4f;   // N_MOVE_EPS
        auto gate = [&](float fs, float lfoHz, double secs, long& slices, long& completions) {
            float lutN = -1.f, nSettleRef = -1e9f, buildN = -1.f;
            int rebuildClock = 0, buildPos = -1; bool lutValid = false;
            slices = 0; completions = 0;
            long ns = (long) (secs * fs);
            for (long s = 0; s < ns; s++) {
                float n = 4.5f + 2.f * std::sin(2.0 * M_PI * lfoHz * s / fs);   // 2.5..6.5, in HILL range
                if (std::fabs(n - nSettleRef) > EPS) { nSettleRef = n; rebuildClock = 0; buildPos = -1; }
                else if (rebuildClock < 2048) ++rebuildClock;
                if (buildPos < 0 && rebuildClock >= 2048 && std::fabs(n - lutN) > 1e-4f) {
                    buildPos = 0; buildN = n; lutValid = false;
                }
                if (buildPos >= 0) {
                    int end = std::min(buildPos + SLICE, LN + 1);
                    if (end > buildPos) ++slices;
                    buildPos = end;
                    if (buildPos > LN) { lutN = buildN; lutValid = true; buildPos = -1; ++completions; }
                }
            }
        };
        long slices, comps;
        // (a) slow LFOs over 10 s: at most a few builds (turning-point settling), never
        //     the tens-of-thousands-of-slices-per-second thrash. Bound generously at the
        //     equivalent of ~20 full rebuilds (a full rebuild is ~33 slices).
        bool gateOk = true;
        for (float fs : {48000.f, 192000.f})
            for (float hz : {0.1f, 0.25f, 1.0f}) {
                gate(fs, hz, 10.0, slices, comps);
                if (slices > 20 * 33) { gateOk = false;
                    printf("  LUT gate THRASH: fs=%.0f lfo=%.2f -> %ld slices in 10 s\n", fs, hz, slices); }
            }
        // (b) a static n must still build the table.
        gate(48000.f, 0.0f, 1.0, slices, comps);
        if (comps < 1) { gateOk = false; printf("  LUT gate: static n never completed a build\n"); }
        printf("Hill LUT gate: slow-LFO no-thrash + static completes: %s\n", gateOk ? "PASS" : "FAIL");
        if (!gateOk) ++fails;
    }

    if (fails) { printf("FAIL: %d unbounded runs\n", fails); return 1; }
    if (lutMax > 1e-3) { printf("FAIL: Hill LUT error too large\n"); return 1; }
    if (pStarMax > 1e-4) { printf("FAIL: warm-started Newton pStar does not converge\n"); return 1; }
    printf("PASS: repressilator finite and bounded; Hill LUT within 1e-3 of pow; warm Newton pStar converges\n");
    return 0;
}
