// Regression test for Haptik's core safety invariant, pitch accuracy, and the
// Slow-mode FREEZE readout continuity (the captured interpolated frame).
//
// Replicates the src/Haptik.cpp dynamics kernel (modulo-free symplectic Euler
// with the +1e-20f denormal guard) in plain C++ so the invariant can be checked
// without launching Rack. Keep this in sync with the kernel if the math changes.
//
//   Build & run:  g++ -O2 -o /tmp/haptik_stability stability_test.cpp && /tmp/haptik_stability
//   Exit 0 = all checks pass, 1 = a stability check failed.
//
// Invariant under test: with COUPLE (kSpr) clamped at 0.9 and the
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
    float gamma = std::exp(-damp * damp * 250.f * D / fs); // DAMP_MAX_HZ=250, quadratic taper, D-scaled
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

    // [6] Slow-mode FREEZE captures the frame currently being *heard* (the
    //     interpolated shape between yPrev and y), so engaging FREEZE mid-frame
    //     doesn't jump the readout — the click the fix removed. Mirror the module:
    //     step every D samples keeping yPrev; the live readout lerps yPrev->y by
    //     fr = divCounter/D (src/Haptik.cpp:398). The FREEZE edge collapses that same
    //     interpolated frame into y (and yPrev) (src/Haptik.cpp:324), after which the
    //     frozen readout reads y directly — and must equal the live value. The old
    //     bug read the raw y endpoint instead, jumping by (1-fr)x the frame delta.
    {
        const int   N = 64, D = 256;                 // D = SLOW_DIV
        const float kSpr = 0.3f, fEvo = 3.f, drive = 0.5f, gamma = 1.f;  // DAMP=0
        const float wc = 2.f*(float)M_PI*fEvo*D/fs;
        const float kCtr = std::min(wc*wc, 0.35f);
        const int   targetDiv = D/2;                 // freeze halfway through a frame

        seedBump(N);
        static float yPrev[MAX_N]; std::fill(yPrev, yPrev + MAX_N, 0.f);
        int divc = 0;
        for (long n = 0; n < 5L*D + targetDiv; n++) { // evolve so yPrev and y differ
            if (divc == 0) {
                std::copy(y, y + N, yPrev);            // pre-step snapshot (module:352)
                v[0] = (v[0] + kSpr*(y[N-1]-2*y[0]+y[1]) - kCtr*y[0]) * gamma + 1e-20f;
                for (int i = 1; i < N-1; i++)
                    v[i] = (v[i] + kSpr*(y[i-1]-2*y[i]+y[i+1]) - kCtr*y[i]) * gamma + 1e-20f;
                v[N-1] = (v[N-1] + kSpr*(y[N-2]-2*y[N-1]+y[0]) - kCtr*y[N-1]) * gamma + 1e-20f;
                v[driverIdx] += drive * gamma;
                for (int i = 0; i < N; i++) {
                    y[i] = std::max(-STATE_MAX, std::min(STATE_MAX, y[i] + v[i]));
                    v[i] = std::max(-STATE_MAX, std::min(STATE_MAX, v[i]));
                }
            }
            divc = (divc + 1) % D;
        }
        // divc == targetDiv; yPrev = last pre-step frame, y = current frame.
        float fr = (float) divc / (float) D;
        float ycap[MAX_N];                            // FREEZE capture (module:324-329)
        for (int i = 0; i < N; i++) ycap[i] = yPrev[i] + fr*(y[i] - yPrev[i]);

        double maxJumpFixed = 0.0, maxJumpBug = 0.0;
        for (int i0 = 0; i0 < N; i0++) {
            int i1 = (i0 + 1) % N;
            for (float f : {0.0f, 0.37f, 0.5f, 0.83f}) {
                float a0 = yPrev[i0] + fr*(y[i0] - yPrev[i0]);   // live slow readout
                float a1 = yPrev[i1] + fr*(y[i1] - yPrev[i1]);
                float live        = a0 + f*(a1 - a0);
                float frozenFixed = ycap[i0] + f*(ycap[i1] - ycap[i0]);  // reads captured frame
                float frozenBug   = y[i0]    + f*(y[i1]    - y[i0]);     // old: reads raw endpoint
                maxJumpFixed = std::max(maxJumpFixed, (double) std::fabs(frozenFixed - live));
                maxJumpBug   = std::max(maxJumpBug,   (double) std::fabs(frozenBug   - live));
            }
        }
        bool contOk = maxJumpFixed < 1e-5;            // capture holds the heard shape
        bool discOk = maxJumpBug   > 1e-2;            // raw-endpoint read really would click
        printf("[6] Slow-mode FREEZE readout continuity (freeze at fr=%.2f)\n", fr);
        printf("    captured-frame jump=%.2e (raw-endpoint jump would be %.3f)\n",
               maxJumpFixed, maxJumpBug);
        if (!contOk) printf("    FAIL: FREEZE capture jumps the readout\n");
        if (!discOk) printf("    FAIL discrimination: raw-endpoint read wouldn't jump — test vacuous\n");
        printf("    %s\n", (contOk && discOk) ? "PASS" : "FAIL");
        ok = ok && contOk && discOk;
    }

    printf("%s\n", ok ? "ALL PASS" : "FAILURES PRESENT");
    return ok ? 0 : 1;
}
