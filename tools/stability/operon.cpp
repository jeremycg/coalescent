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

    if (fails) { printf("FAIL: %d unbounded runs\n", fails); return 1; }
    printf("PASS: repressilator finite and bounded across the parameter box\n");
    return 0;
}
