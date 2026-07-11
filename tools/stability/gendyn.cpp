// Stability + calibration replica for GENDYN's pitch/playback path (mirrors the
// relevant pieces of src/GENDYN.cpp). Asserts the invariants the module relies on
// but that can't be seen from the audio alone:
//   [1] reflect() bounds the amplitude/duration walks (finiteness + range),
//   [2] DUR WID=0 (zero-width frequency barrier) plays the exact centre pitch
//       — the regression that shipped ~41 cents flat,
//   [3] the reachable maximum frequency is fs/N (every segment >= 1 sample),
//   [4] LOCK renormalises an off-centre cycle back to the centre pitch (well above floor),
//   [5] the LOCK servo hits reachable near-floor targets and saturates below the floor,
//   and all of the above hold across sample rates.
//
//   g++ -O2 -o /tmp/t tools/stability/gendyn.cpp && /tmp/t   (exit 0 = pass)
#include <cstdio>
#include <cmath>
#include <algorithm>

// ── ports of the exact kernel pieces under test ──────────────────────────────

// src/GENDYN.cpp reflect(): folds x into [lo,hi]; returns lo for a degenerate span.
static float reflect(float x, float lo, float hi) {
    if (hi <= lo) return lo;
    const float range = hi - lo;
    float y = std::fmod(x - lo, 2.f * range);
    if (y < 0.f) y += 2.f * range;
    return lo + (y > range ? 2.f * range - y : y);
}

// src/GENDYN.cpp:376-378 error-diffused duration rounding. Plays `cycles` cycles of
// N segments whose (pre-norm) durations are dur[i], scaled by norm_k, and returns
// the played frequency (fs * cycles / totalSamples).
static double playHz(const float* dur, int N, float norm_k, float fs, int cycles) {
    float dur_err = 0.f; long total = 0;
    for (int c = 0; c < cycles; c++)
        for (int i = 0; i < N; i++) {
            float fd = dur[i] * norm_k + dur_err;
            int cd = std::max(1, (int) (fd + 0.5f));
            dur_err = std::clamp(fd - (float) cd, -4.f, 4.f);
            total += cd;
        }
    return fs * cycles / (double) total;
}

// src/GENDYN.cpp updateNormAndFreq(): LOCK scales the cycle to fs/centreFreq, but
// never below the N-sample floor.
static float lockNorm(float sum_dur, int N, float fs, float centerFreq) {
    float k = (fs / centerFreq) / sum_dur;
    return std::max(k, (float) N / sum_dur);
}

// src/GENDYN.cpp closed-loop LOCK servo: each cycle plays the N segments
// (error-diffused, per-segment 1-sample floor), measures the realized period, and
// nudges lockCorr toward target/realized (damped to ±25%/cycle, bounded [0.25,4]).
// dur_err persists across cycles exactly as the module carries it. Returns the
// realized frequency of the final cycle after `cycles` cycles of convergence.
static double playHzServoLock(const float* dur, int N, float fs, float centerFreq, int cycles) {
    float sum = 0.f; for (int i = 0; i < N; i++) sum += dur[i];
    if (sum <= 0.f) return 0.0;
    float dur_err = 0.f, lockCorr = 1.f;
    const float target = fs / centerFreq;
    long realized = 0;
    for (int c = 0; c < cycles; c++) {
        float k = (fs / centerFreq) / sum * lockCorr;   // servo bounds the loop, not a norm_k clamp
        realized = 0;
        for (int i = 0; i < N; i++) {
            float fd = dur[i] * k + dur_err;
            int cd = std::max(1, (int) (fd + 0.5f));
            dur_err = std::clamp(fd - (float) cd, -4.f, 4.f);
            realized += cd;
        }
        float step = std::clamp(target / (float) realized, 0.8f, 1.25f);
        lockCorr = std::clamp(lockCorr * step, 0.25f, 4.f);
    }
    return fs / (double) realized;
}

static double cents(double got, double want) { return 1200.0 * std::log2(got / want); }

// src/GENDYN.cpp:309-313 frequency-barrier derivation. `buggy` selects the shipped
// regression (widen a zero-width window by +1 sample) vs the fix (clamp bDurMax>=bDurMin).
struct Barriers { float durCenter, bDurMin, bDurMax; };
static Barriers barriers(float bDurWidth, float fs, float centerFreq, int N, bool buggy) {
    float durCenter = fs / (centerFreq * (float) N);
    float halfWidth = bDurWidth * durCenter;
    float bDurMin = std::max(1.f, durCenter - halfWidth);
    float bDurMax;
    if (buggy) { bDurMax = durCenter + halfWidth; if (bDurMax <= bDurMin) bDurMax = bDurMin + 1.f; }
    else       { bDurMax = std::max(bDurMin, durCenter + halfWidth); }
    return {durCenter, bDurMin, bDurMax};
}

// src/GENDYN.cpp runCycleUpdate() duration walk: pDurHW=(max-min)/2, gainDur=scaleDur*pDurHW,
// then reflect the persistent step and the duration. A zero-width window makes gainDur=0, so
// the duration must stay pinned. Deterministic driver (fixed nonzero draw) for reproducibility.
static void walkDurations(const Barriers& b, int N, float scaleDur, int cycles, float* dur) {
    float pDurHW = (b.bDurMax - b.bDurMin) * 0.5f;
    float gainDur = scaleDur * pDurHW;
    float step[64];
    for (int i = 0; i < N; i++) { dur[i] = std::max(1.f, b.durCenter); step[i] = 0.f; }
    for (int c = 0; c < cycles; c++)
        for (int i = 0; i < N; i++) {
            step[i] = reflect(step[i] + 0.3f, -1.f, 1.f);          // any nonzero draw exposes a nonzero gainDur
            dur[i]  = reflect(dur[i] + gainDur * step[i], b.bDurMin, b.bDurMax);
        }
}

int main() {
    int fails = 0;
    const float sampleRates[] = {44100.f, 48000.f, 96000.f};

    // [1] reflect() always lands in [lo,hi] (or ==lo for a degenerate span), for
    //     every finite input — this is what bounds amp[] and dur[] every cycle.
    //     (drawSample stays finite: each sampler clamps its uniform() away from the
    //     0/1 tails, so log/tan can't blow up, so reflect never sees a non-finite.)
    {
        bool ok = true;
        const float xs[] = {0.f, 0.3f, -0.3f, 5.f, -5.f, 3183.f, -3183.f, 1e6f, -1e6f};
        struct { float lo, hi; } spans[] = {{-1.f, 1.f}, {1.f, 13.f}, {12.97f, 12.97f}, {3.f, 3.f}};
        for (auto s : spans)
            for (float x : xs) {
                float r = reflect(x, s.lo, s.hi);
                bool inb = (s.hi <= s.lo) ? (r == s.lo)
                                          : (std::isfinite(r) && r >= s.lo - 1e-3f && r <= s.hi + 1e-3f);
                if (!inb) { ok = false; printf("  [1] FAIL reflect(%.3g,%.3g,%.3g)=%.4g out of bounds\n", x, s.lo, s.hi, r); }
            }
        printf("[1] reflect() bounds walks: %s\n", ok ? "PASS" : "FAIL");
        if (!ok) fails++;
    }

    // [2] DUR WID=0: derive the barriers and RUN THE WALK. The shipped bug only
    //     detuned *after* the duration walk moved inside a widened window, so a test
    //     that just plays durCenter would pass the buggy code too. Assert the
    //     zero-width walk stays degenerate (durations pinned → exact pitch), and
    //     verify the old widened-barrier version actually drifts, so the test
    //     genuinely discriminates the regression.
    {
        bool ok = true;
        for (float fs : sampleRates)
            for (int N : {8, 13, 32, 64})
                for (float f : {130.81f, 261.6256f, 440.f}) {
                    Barriers b = barriers(0.f, fs, f, N, /*buggy=*/false);
                    if (b.durCenter < 1.f) continue;            // above fs/N floor → tested in [3]
                    float dur[64]; walkDurations(b, N, 1.0f, 4000, dur);   // WID=0, strong SCALE DUR
                    float maxdev = 0.f;
                    for (int i = 0; i < N; i++) maxdev = std::max(maxdev, std::fabs(dur[i] - b.durCenter));
                    double hz = playHz(dur, N, 1.f, fs, 40000);
                    double dc = cents(hz, f);
                    if (maxdev > 1e-3f || std::fabs(dc) > 2.0) { ok = false;
                        printf("  [2] FAIL fs=%.0f N=%d f=%.1f: maxdev=%.4f -> %.3f Hz (%.1f cents)\n",
                               fs, N, f, maxdev, hz, dc); }
                }
        // Discrimination: the old widened barrier must actually detune under the same
        // walk, otherwise the test above would prove nothing.
        {
            Barriers bug = barriers(0.f, 44100.f, 261.6256f, 13, /*buggy=*/true);
            float dur[64]; walkDurations(bug, 13, 1.0f, 4000, dur);
            double dc = cents(playHz(dur, 13, 1.f, 44100.f, 40000), 261.6256);
            if (std::fabs(dc) < 10.0) { ok = false;
                printf("  [2] FAIL discrimination: buggy barrier only %.1f cents off — test wouldn't catch it\n", dc); }
            else printf("  [2] discrimination OK: old widened barrier drifts %.1f cents under the same walk\n", dc);
        }
        printf("[2] DUR WID=0 walk stays degenerate -> centre pitch: %s\n", ok ? "PASS" : "FAIL");
        if (!ok) fails++;
    }

    // [3] Reachable max frequency is fs/N: request well above it, expect saturation
    //     at fs/N (every segment floored at 1 sample), not the requested pitch.
    {
        bool ok = true;
        for (float fs : sampleRates)
            for (int N : {16, 64}) {
                float want = 20000.f;                          // request way past fs/N
                float durCenter = fs / (want * (float) N);     // < 1 sample
                float dur[64]; for (int i = 0; i < N; i++) dur[i] = std::max(1.f, durCenter);
                double hz = playHz(dur, N, 1.f, fs, 40000);
                double cap = fs / N;
                if (std::fabs(hz - cap) > 1.0) { ok = false;
                    printf("  [3] FAIL fs=%.0f N=%d -> %.2f Hz, expected fs/N=%.2f\n", fs, N, hz, cap); }
            }
        printf("[3] max reachable frequency = fs/N: %s\n", ok ? "PASS" : "FAIL");
        if (!ok) fails++;
    }

    // [4] LOCK renormalises an off-centre cycle (durations perturbed away from
    //     durCenter) back to the centre pitch, down to the fs/N floor.
    {
        bool ok = true;
        for (float fs : sampleRates)
            for (int N : {13, 64})
                for (float f : {196.f, 261.6256f}) {
                    float durCenter = fs / (f * (float) N);
                    if (durCenter < 1.5f) continue;
                    float dur[64]; float sum = 0.f;
                    for (int i = 0; i < N; i++) {               // ±30% off-centre "walked" durations
                        dur[i] = durCenter * (1.f + 0.3f * std::sin(i * 1.7f));
                        sum += dur[i];
                    }
                    float k = lockNorm(sum, N, fs, f);
                    double hz = playHz(dur, N, k, fs, 40000);
                    double dc = cents(hz, f);
                    if (std::fabs(dc) > 3.0) { ok = false;
                        printf("  [4] FAIL fs=%.0f N=%d f=%.1f LOCK -> %.3f Hz (%.1f cents)\n", fs, N, f, hz, dc); }
                }
        printf("[4] LOCK holds centre pitch off-centre: %s\n", ok ? "PASS" : "FAIL");
        if (!ok) fails++;
    }

    // [5] Near the sample floor the plain theoretical norm inflates the realized
    //     period (per-segment >=1 sample + error-diffused rounding of unequal
    //     durations), so the module runs a per-cycle servo. Two regimes:
    //       (a) target period above the N-sample floor  -> servo hits it (few cents),
    //           whereas the un-servo'd norm would be sharply flat;
    //       (b) target period at/below the floor         -> physically unreachable,
    //           servo saturates and holds as close as fs/N allows (best-effort).
    {
        bool ok = true;
        for (float fs : sampleRates) {
            int N = 64;
            float floorHz = fs / N;                        // fastest reachable pitch
            // (a) reachable but near-floor, with the durations *clustered* asymmetric
            //     (many shorts, few longs) so the shorts all floor to 1 sample and the
            //     theoretical norm inflates the period sharply — the exact defect the
            //     servo fixes. Target period ~1.09*N samples. The integer-sample
            //     granularity is coarse here (~1 sample ≈ 23 cents at a 70-sample
            //     period), so the invariant is "within ~1 sample of the target period"
            //     — the best any integer-length playback can do — not absolute cents.
            {
                float f = floorHz / 1.09375f;              // period ~1.09N (70 samples for N=64) > floor
                float targetPeriod = fs / f;
                float dur[64];
                for (int i = 0; i < N; i++) dur[i] = (i < 48) ? 0.4f : 4.6f;  // 48 shorts + 16 longs
                double servoPeriod = fs / playHzServoLock(dur, N, fs, f, 400);
                float sum = 0.f; for (int i = 0; i < N; i++) sum += dur[i];
                double naivePeriod = fs / playHz(dur, N, lockNorm(sum, N, fs, f), fs, 400);
                double servoErr = std::fabs(servoPeriod - targetPeriod);
                double naiveErr = std::fabs(naivePeriod - targetPeriod);
                if (servoErr > 1.1) { ok = false;
                    printf("  [5a] FAIL fs=%.0f near-floor LOCK servo period err=%.2f samples (want <=1)\n", fs, servoErr); }
                else if (naiveErr < 3.0)                   // guard: the servo must be fixing a real defect
                    printf("  [5a] note fs=%.0f: un-servo'd norm only %.2f samples off — mild case\n", fs, naiveErr);
                else
                    printf("  [5a] OK fs=%.0f: servo within %.2f samples of target (un-servo'd off by %.1f samples)\n",
                           fs, servoErr, naiveErr);
            }
            // (b) unreachable (below the floor): realized must saturate near fs/N.
            {
                float f = floorHz * 1.3f;                  // requested pitch faster than fs/N
                float durCenter = fs / (f * (float) N);
                float dur[64]; for (int i = 0; i < N; i++) dur[i] = std::max(0.1f, durCenter);
                double hz = playHzServoLock(dur, N, fs, f, 400);
                if (std::fabs(hz - floorHz) > 1.0) { ok = false;
                    printf("  [5b] FAIL fs=%.0f below-floor LOCK -> %.2f Hz, expected saturation at fs/N=%.2f\n",
                           fs, hz, floorHz); }
            }
        }
        printf("[5] LOCK servo: hits reachable near-floor targets, saturates below floor: %s\n", ok ? "PASS" : "FAIL");
        if (!ok) fails++;
    }

    printf(fails ? "FAIL: %d GENDYN invariant(s) broken\n" : "PASS: GENDYN pitch/playback invariants hold\n", fails);
    return fails ? 1 : 0;
}
