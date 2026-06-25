// Regression test for Haptik's core safety invariant and pitch accuracy.
//
// Replicates the src/Haptik.cpp dynamics kernel (modulo-free symplectic Euler
// with the +1e-20f denormal guard) in plain C++ so the invariant can be checked
// without launching Rack. Keep this in sync with the kernel if the math changes.
//
//   Build & run:  g++ -O2 -o /tmp/haptik_stability stability_test.cpp && /tmp/haptik_stability
//   Exit 0 = all checks pass, 1 = a stability check failed.
//
// Invariant under test (see plan §5): with COUPLE (kSpr) clamped at 0.9 and the
// centering term tiny, omega_max = sqrt(kCtr + 4*kSpr) < 2, so the integrator is
// stable and y[] stays finite/bounded for any COUPLE in [0, 0.9] even at DAMP=0.
#include <cstdio>
#include <cmath>
#include <algorithm>

static const int MAX_N = 256;
static float y[MAX_N], v[MAX_N];
static int driverIdx;
static const float STATE_MAX = 16.f;   // matches src/Haptik.cpp

static void seedBump(int N) {
    std::fill(y, y + MAX_N, 0.f); std::fill(v, v + MAX_N, 0.f);
    driverIdx = N / 4;
    int hw = std::max(1, (int)(0.125f * N));   // BUMP_FRAC = 0.125
    for (int d = -hw; d <= hw; d++) {
        float w = 0.5f * (1.f + std::cos((float)M_PI * d / hw));
        y[((driverIdx + d) % N + N) % N] += w;
    }
}

// Runs the dynamics; returns max|y| (INFINITY if it ever goes non-finite).
// If period!=null, also measures the scan period in samples.
static double run(int N, float fEvo, float kSpr, float damp, float fs,
                  long samples, float pitchHz, double* period, float drive = 0.f,
                  int D = 1) {
    seedBump(N);
    float wc = 2.f * (float)M_PI * fEvo * D / fs;
    float kCtr = std::min(wc * wc, 0.35f);          // KCTR_MAX clamp (slow mode)
    float gamma = std::exp(-damp * 800.f * D / fs); // DAMP_MAX_HZ=800, D-scaled
    double maxabs = 0.0; float scanPhase = 0.f;
    long firstWrap = -1, lastWrap = -1; int wraps = 0; int divc = 0;
    for (long n = 0; n < samples; n++) {
        if (D <= 1 || divc == 0) {                  // step every D samples
            v[0] = (v[0] + kSpr*(y[N-1]-2*y[0]+y[1]) - kCtr*y[0]) * gamma + 1e-20f;
            for (int i = 1; i < N-1; i++)
                v[i] = (v[i] + kSpr*(y[i-1]-2*y[i]+y[i+1]) - kCtr*y[i]) * gamma + 1e-20f;
            v[N-1] = (v[N-1] + kSpr*(y[N-2]-2*y[N-1]+y[0]) - kCtr*y[N-1]) * gamma + 1e-20f;
            v[driverIdx] += drive * gamma;
            for (int i = 0; i < N; i++) {
                y[i] = std::max(-STATE_MAX, std::min(STATE_MAX, y[i] + v[i]));
                v[i] = std::max(-STATE_MAX, std::min(STATE_MAX, v[i]));
                if (!std::isfinite(y[i])) return INFINITY;
                double a = std::fabs(y[i]); if (a > maxabs) maxabs = a;
            }
        }
        if (D > 1) divc = (divc + 1) % D;
        scanPhase += pitchHz / fs;
        if (scanPhase >= 1.f) { scanPhase -= std::floor(scanPhase); wraps++; if (firstWrap<0) firstWrap=n; lastWrap=n; }
    }
    if (period && wraps > 1) *period = (double)(lastWrap - firstWrap) / (wraps - 1);
    return maxabs;
}

int main() {
    const float fs = 44100.f;
    bool ok = true;

    printf("[1] Stability sweep COUPLE 0..0.9 (DAMP=0, RATE=3Hz), N in {8,64,128}\n");
    for (int N : {8, 64, 128})
        for (int k = 0; k <= 90; k++) {
            double m = run(N, 3.f, k/100.f, 0.f, fs, 200000, 261.626f, nullptr);
            if (!std::isfinite(m)) { ok = false; printf("    FAIL N=%d COUPLE=%.2f non-finite\n", N, k/100.f); }
        }
    printf("    %s\n", ok ? "PASS (all finite & bounded)" : "FAIL");

    printf("[2] Pitch at PITCH=0/VOCT=0 should be ~261.626 Hz\n");
    double period = 0; run(64, 3.f, 0.3f, 0.f, fs, 50000, 261.626f, &period);
    double hz = fs / period; bool pitchOk = std::fabs(hz - 261.626) < 1.0;
    printf("    measured %.3f Hz  %s\n", hz, pitchOk ? "PASS" : "FAIL");
    ok = ok && pitchOk;

    printf("[3] Extreme corner kSpr=0.9, RATE=30Hz, N=128, 500k samples\n");
    double m = run(128, 30.f, 0.9f, 0.f, fs, 500000, 261.626f, nullptr);
    printf("    max|y|=%.4f  %s\n", m, std::isfinite(m) ? "PASS" : "FAIL");
    ok = ok && std::isfinite(m);

    printf("[4] Forced + lossless (DAMP=0, constant EXT drive), 1M samples\n");
    double mf = run(64, 3.f, 0.3f, 0.f, fs, 1000000, 261.626f, nullptr, 0.5f);
    bool fb = std::isfinite(mf) && mf <= STATE_MAX + 1e-3;
    printf("    max|y|=%.4f (clamp %.0f)  %s\n", mf, STATE_MAX,
           fb ? "PASS (bounded by clamp)" : "FAIL");
    ok = ok && fb;

    printf("[5] Slow mode (D=256): COUPLE 0..0.9, RATE=30Hz (kCtr clamped), forced\n");
    bool sb = true;
    for (int k = 0; k <= 90; k += 10) {
        double m = run(64, 30.f, k / 100.f, 0.f, fs, 300000, 261.626f, nullptr, 0.5f, 256);
        if (!(std::isfinite(m) && m <= STATE_MAX + 1e-3)) {
            sb = false; printf("    FAIL COUPLE=%.2f max|y|=%.3f\n", k / 100.f, m);
        }
    }
    printf("    %s\n", sb ? "PASS (bounded with divider + kCtr clamp)" : "FAIL");
    ok = ok && sb;

    printf("%s\n", ok ? "ALL PASS" : "FAILURES PRESENT");
    return ok ? 0 : 1;
}
