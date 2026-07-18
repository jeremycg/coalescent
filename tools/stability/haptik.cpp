// SDK-free tests for Haptik's shared production lattice core: stability, pitch
// read-head accuracy, and Slow-mode FREEZE/readout continuity. The Rack wrapper
// and this binary call dsp/haptik_core.hpp; explicit raw-endpoint alternatives
// below are deliberate historical mutants used only for discrimination.
//
//   Build & run:  g++ -O2 -o /tmp/haptik_stability stability_test.cpp && /tmp/haptik_stability
//   Exit 0 = all checks pass, 1 = a stability check failed.
//
// Invariant under test: with COUPLE (kSpr) at its shared production maximum and the
// centering term tiny, omega_max = sqrt(kCtr + 4*kSpr) < 2, so the integrator is
// stable and y[] stays finite/bounded for any COUPLE in [0, 0.9] even at DAMP=0.
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <limits>
#include "../../src/dsp/haptik_core.hpp"

static const int MAX_N = coalescent::haptik::kMaxStorageMasses;
static float y[MAX_N], v[MAX_N];
static int driverIdx;
static const float STATE_MAX = coalescent::haptik::kStateMax;

static void seedBump(int N) {
    std::fill(y, y + MAX_N, 0.f); std::fill(v, v + MAX_N, 0.f);
    driverIdx = N / 4;
    const float added = coalescent::haptik::addHannBump(y, N, driverIdx, 1.f);
    coalescent::haptik::removeInjectedMean(y, N, added);
}

// Runs the dynamics; returns max|y| (INFINITY if it ever goes non-finite).
// If period!=null, also measures the scan period in samples.
static double run(int N, float fEvo, float kSpr, float damp, float fs,
                  long samples, float pitchHz, double* period, float drive = 0.f,
                  int D = 1) {
    seedBump(N);
    const float kCtr = coalescent::haptik::centeringCoefficient(fEvo, D, fs);
    // The production wrapper evaluates this shared exponent with Rack's exp2
    // approximation. The SDK-free stability oracle uses libm exp2 instead.
    const float gamma = std::exp2(
        coalescent::haptik::dampingExp2Argument(damp, D, fs));
    double maxabs = 0.0; float scanPhase = 0.f;
    long firstWrap = -1, lastWrap = -1; int wraps = 0; int divc = 0;
    for (long n = 0; n < samples; n++) {
        if (coalescent::haptik::shouldStep(false, D > 1, divc)) {
            coalescent::haptik::stepScalar(
                y, v, N, driverIdx, drive, kSpr, kCtr, gamma);
            for (int i = 0; i < N; i++) {
                if (!std::isfinite(y[i])) return INFINITY;
                double a = std::fabs(y[i]); if (a > maxabs) maxabs = a;
            }
        }
        if (D > 1)
            divc = coalescent::haptik::advanceDivider(divc, D);
        const float previousPhase = scanPhase;
        coalescent::haptik::advanceScanPhase(scanPhase, pitchHz, fs);
        if (scanPhase < previousPhase) {
            wraps++;
            if (firstWrap < 0) firstWrap = n;
            lastWrap = n;
        }
    }
    if (period && wraps > 1) *period = (double)(lastWrap - firstWrap) / (wraps - 1);
    return maxabs;
}

int main() {
    const float fs = 44100.f;
    bool ok = true;

    printf("[1] Stability sweep COUPLE 0..0.9 (DAMP=0, RATE=3Hz), N in {8,64,128}\n");
    for (int N : {coalescent::haptik::kMinPlayableMasses, 64,
                  coalescent::haptik::kMaxPlayableMasses})
        for (int k = 0; k <= static_cast<int>(
                 coalescent::haptik::kCouplingMax * 100.f); k++) {
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
    double m = run(coalescent::haptik::kMaxPlayableMasses,
                   coalescent::haptik::kEvolutionMaxHz,
                   coalescent::haptik::kCouplingMax,
                   0.f, fs, 500000, 261.626f, nullptr);
    printf("    max|y|=%.4f  %s\n", m, std::isfinite(m) ? "PASS" : "FAIL");
    ok = ok && std::isfinite(m);

    printf("[4] Forced + lossless (DAMP=0, constant EXT drive), 1M samples\n");
    double mf = run(64, 3.f, 0.3f, 0.f, fs, 1000000, 261.626f, nullptr, 0.5f);
    bool fb = std::isfinite(mf) && mf <= STATE_MAX + 1e-3;
    printf("    max|y|=%.4f (clamp %.0f)  %s\n", mf, STATE_MAX,
           fb ? "PASS (bounded by clamp)" : "FAIL");
    ok = ok && fb;

    printf("[5] Slow mode (D=%d): COUPLE 0..0.9, RATE=30Hz (kCtr clamped), forced\n",
           coalescent::haptik::kSlowDivider);
    bool sb = true;
    for (int k = 0; k <= static_cast<int>(
             coalescent::haptik::kCouplingMax * 100.f); k += 10) {
        double m = run(64, coalescent::haptik::kEvolutionMaxHz,
                       k / 100.f, 0.f, fs, 300000, 261.626f,
                       nullptr, 0.5f, coalescent::haptik::kSlowDivider);
        if (!(std::isfinite(m) && m <= STATE_MAX + 1e-3)) {
            sb = false; printf("    FAIL COUPLE=%.2f max|y|=%.3f\n", k / 100.f, m);
        }
    }
    printf("    %s\n", sb ? "PASS (bounded with divider + kCtr clamp)" : "FAIL");
    ok = ok && sb;

    // [6] Slow-mode FREEZE captures the frame currently being *heard* (the
    //     interpolated shape between yPrev and y), so engaging FREEZE mid-frame
    //     doesn't jump the readout — the click the fix removed. Step every D
    //     samples keeping yPrev; the live readout uses the shared
    //     interpolation primitive. The shared FREEZE capture collapses that same
    //     frame into y and yPrev, after which the frozen readout reads y directly.
    //     The old bug read the raw y endpoint instead, jumping by (1-fr)x the frame delta.
    {
        const int   N = 64, D = coalescent::haptik::kSlowDivider;
        const float kSpr = 0.3f, fEvo = 3.f, drive = 0.5f, gamma = 1.f;  // DAMP=0
        const float kCtr = coalescent::haptik::centeringCoefficient(fEvo, D, fs);
        const int   targetDiv = D/2;                 // freeze halfway through a frame

        seedBump(N);
        static float yPrev[MAX_N]; std::fill(yPrev, yPrev + MAX_N, 0.f);
        int divc = 0;
        for (long n = 0; n < 5L*D + targetDiv; n++) { // evolve so yPrev and y differ
            if (divc == 0) {
                std::copy(y, y + N, yPrev);
                coalescent::haptik::stepScalar(
                    y, v, N, driverIdx, drive, kSpr, kCtr, gamma);
            }
            divc = coalescent::haptik::advanceDivider(divc, D);
        }
        // divc == targetDiv; yPrev = last pre-step frame, y = current frame.
        float fr = (float) divc / (float) D;
        float ycap[MAX_N], capturedPrev[MAX_N];
        std::copy(y, y + N, ycap);
        std::copy(yPrev, yPrev + N, capturedPrev);
        coalescent::haptik::captureInterpolatedFrame(
            capturedPrev, ycap, N, fr);

        double maxJumpFixed = 0.0, maxJumpBug = 0.0, maxCaptureMismatch = 0.0;
        for (int i = 0; i < N; ++i)
            maxCaptureMismatch = std::max(
                maxCaptureMismatch,
                static_cast<double>(std::fabs(capturedPrev[i] - ycap[i])));
        for (int i0 = 0; i0 < N; i0++) {
            int i1 = (i0 + 1) % N;
            for (float f : {0.0f, 0.37f, 0.5f, 0.83f}) {
                float a0 = coalescent::haptik::interpolate(yPrev[i0], y[i0], fr);
                float a1 = coalescent::haptik::interpolate(yPrev[i1], y[i1], fr);
                float live = coalescent::haptik::interpolate(a0, a1, f);
                float frozenFixed = coalescent::haptik::interpolate(
                    ycap[i0], ycap[i1], f);
                float frozenBug   = y[i0]    + f*(y[i1]    - y[i0]);     // old: reads raw endpoint
                maxJumpFixed = std::max(maxJumpFixed, (double) std::fabs(frozenFixed - live));
                maxJumpBug   = std::max(maxJumpBug,   (double) std::fabs(frozenBug   - live));
            }
        }
        bool contOk = maxJumpFixed < 1e-5 && maxCaptureMismatch == 0.0;
        bool discOk = maxJumpBug   > 1e-2;            // raw-endpoint read really would click
        printf("[6] Slow-mode FREEZE readout continuity (freeze at fr=%.2f)\n", fr);
        printf("    captured-frame jump=%.2e (raw-endpoint jump would be %.3f)\n",
               maxJumpFixed, maxJumpBug);
        if (!contOk) printf("    FAIL: FREEZE capture jumps the readout\n");
        if (!discOk) printf("    FAIL discrimination: raw-endpoint read wouldn't jump — test vacuous\n");
        printf("    %s\n", (contOk && discOk) ? "PASS" : "FAIL");
        ok = ok && contOk && discOk;
    }

    // [7] MOTION follows the same yPrev->y interpolation as the scanned OUT in
    // Slow mode. Across a frame boundary the interpolated mass-0 trajectory is
    // continuous to one frame increment, whereas reading raw y[0] produces a
    // full-frame hold followed by the complete endpoint jump.
    {
        const int D = coalescent::haptik::kSlowDivider;
        const float frame[3] = {-0.8f, 0.8f, -0.2f};
        float prevLerp = frame[0], prevRaw = frame[1];
        double maxLerpStep = 0.0, maxRawStep = 0.0;
        for (int segment = 0; segment < 2; ++segment) {
            for (int divc = 0; divc < D; ++divc) {
                const float fr = (float) divc / (float) D;
                const float lerp = coalescent::haptik::interpolate(
                    frame[segment], frame[segment + 1], fr);
                const float raw = frame[segment + 1];
                maxLerpStep = std::max(maxLerpStep, (double) std::fabs(lerp - prevLerp));
                maxRawStep = std::max(maxRawStep, (double) std::fabs(raw - prevRaw));
                prevLerp = lerp;
                prevRaw = raw;
            }
        }
        const double largestFrameDelta = 1.6;
        bool motionOk = maxLerpStep <= largestFrameDelta / D + 1e-6;
        bool motionDiscriminates = maxRawStep > 0.5;
        printf("[7] Slow-mode MOTION interpolation continuity\n");
        printf("    interpolated max step=%.6f; raw endpoint max step=%.3f: %s\n",
               maxLerpStep, maxRawStep,
               (motionOk && motionDiscriminates) ? "PASS" : "FAIL");
        ok = ok && motionOk && motionDiscriminates;
    }

    {
        const float amount = 0.6f;
        const float finiteVoltage = -7.25f;
        const float expected = finiteVoltage * coalescent::haptik::kExternalGain * amount;
        const bool driveGuardOk =
            coalescent::haptik::externalDrive(finiteVoltage, amount) == expected
            && coalescent::haptik::externalDrive(
                std::numeric_limits<float>::quiet_NaN(), amount) == 0.f
            && coalescent::haptik::externalDrive(
                std::numeric_limits<float>::infinity(), amount) == 0.f
            && coalescent::haptik::externalDrive(
                -std::numeric_limits<float>::infinity(), amount) == 0.f;
        printf("[8] EXT IN non-finite poly-sum neutralization\n");
        printf("    finite mapping %.8f; hostile values -> 0: %s\n",
               expected, driveGuardOk ? "PASS" : "FAIL");
        ok = ok && driveGuardOk;
    }

    printf("%s\n", ok ? "ALL PASS" : "FAILURES PRESENT");
    return ok ? 0 : 1;
}
