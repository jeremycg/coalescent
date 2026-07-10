// Stability + calibration replica for GENDYN's pitch/playback path (mirrors the
// relevant pieces of src/GENDYN.cpp). Asserts the invariants the module relies on
// but that can't be seen from the audio alone:
//   [1] reflect() bounds the amplitude/duration walks (finiteness + range),
//   [2] DUR WID=0 (zero-width frequency barrier) plays the exact centre pitch
//       — the regression that shipped ~41 cents flat,
//   [3] the reachable maximum frequency is fs/N (every segment >= 1 sample),
//   [4] LOCK renormalises an off-centre cycle back to the centre pitch,
//   [5] all of the above hold across sample rates.
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

static double cents(double got, double want) { return 1200.0 * std::log2(got / want); }

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

    // [2] DUR WID=0 → durations pinned at durCenter → exact centre pitch, every SR.
    {
        bool ok = true;
        for (float fs : sampleRates)
            for (int N : {8, 13, 32, 64})
                for (float f : {130.81f, 261.6256f, 440.f}) {
                    float durCenter = fs / (f * (float) N);
                    if (durCenter < 1.f) continue;              // above fs/N floor → tested in [3]
                    float dur[64]; for (int i = 0; i < N; i++) dur[i] = durCenter;
                    double hz = playHz(dur, N, 1.f, fs, 40000);
                    double dc = cents(hz, f);
                    if (std::fabs(dc) > 2.0) { ok = false;
                        printf("  [2] FAIL fs=%.0f N=%d f=%.1f -> %.3f Hz (%.1f cents)\n", fs, N, f, hz, dc); }
                }
        printf("[2] DUR WID=0 plays centre pitch (was ~41 cents flat): %s\n", ok ? "PASS" : "FAIL");
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

    printf(fails ? "FAIL: %d GENDYN invariant(s) broken\n" : "PASS: GENDYN pitch/playback invariants hold\n", fails);
    return fails ? 1 : 0;
}
