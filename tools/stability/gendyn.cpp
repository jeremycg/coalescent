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
#include <vector>
#include <utility>

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

// src/GENDYN.cpp updateNormAndFreq() feed-forward scale: LOCK scales the cycle to
// fs/centreFreq (lockCorr = 1). Production no longer clamps norm_k to the floor —
// the per-cycle servo bounds the loop instead (exercised in test [5]); well above
// the floor (test [4]'s cases) the scale is simply fs/centreFreq/sum_dur.
static float lockNorm(float sum_dur, int N, float fs, float centerFreq) {
    (void) N;
    return (fs / centerFreq) / sum_dur;
}

// Mode selector for the servo mirror below: the shipped fix, and the two historical
// regressions the tests must still discriminate against.
enum ServoMode {
    SERVO_FIXED,        // cycleTarget feedback + anti-windup (>N) + reset-on-target-change
    BUG_FEEDBACK,       // pre-fix #1: feed the completed cycle back against the *next* target
    BUG_STALE,          // pre-fix #2: no anti-windup, no reset — carry lockCorr across targets
};

// Faithful mirror of the src/GENDYN.cpp LOCK servo (updateNormAndFreq + the cycle-
// boundary feedback). `sched` is the per-cycle (centerFreq, lockOn) the module sees;
// the shape `dur[]` is held fixed (the walk is exercised elsewhere). Returns the
// realized period of every cycle so transitions can be inspected; *finalCorr, if
// given, receives the converged lockCorr for restore tests.
//
// SERVO_FIXED mirrors the module exactly: feedback for a completed cycle uses the
// target THAT cycle was rendered against; the servo only integrates for reachable
// targets (period > the N-sample floor); and a material target change / unreachable
// target / LOCK toggle resets lockCorr to unity so a stale correction can't detune
// the next cycle. The two BUG_ modes drop one guard each so the tests can prove the
// fix discriminates.
static std::vector<long> servoRun(const float* dur, int N, float fs,
                                  const std::vector<std::pair<float,bool>>& sched,
                                  ServoMode mode = SERVO_FIXED, float lockCorr0 = 1.f,
                                  float* finalCorr = nullptr) {
    float sum = 0.f; for (int i = 0; i < N; i++) sum += dur[i];
    float lockCorr = lockCorr0, dur_err = 0.f, cycleTarget = 0.f, norm_k = 1.f;
    auto prime = [&](std::pair<float,bool> s) {             // reinit/restore branch primes cycle 0
        if (s.second) norm_k = (fs / s.first) / sum * lockCorr;
        else { norm_k = 1.f; lockCorr = 1.f; }
        cycleTarget = s.second ? fs / s.first : 0.f;
    };
    prime(sched[0]);
    std::vector<long> out; out.reserve(sched.size());
    for (size_t c = 0; c < sched.size(); c++) {
        long realized = 0;
        for (int i = 0; i < N; i++) {
            float fd = dur[i] * norm_k + dur_err;
            int   cd = std::max(1, (int) (fd + 0.5f));       // module's error-diffused round
            dur_err  = std::clamp(fd - (float) cd, -4.f, 4.f);
            realized += cd;
        }
        out.push_back(realized);

        // ── servo feedback for the just-finished cycle ──
        float fbT = cycleTarget;
        if (mode == BUG_FEEDBACK && c + 1 < sched.size())    // #1: feed back against the NEXT target
            fbT = sched[c+1].second ? fs / sched[c+1].first : 0.f;
        bool reachable = (mode == SERVO_FIXED) ? (fbT > (float) N) : (fbT > 0.f);  // anti-windup only when fixed
        if (fbT > 0.f && reachable && realized > 0)
            lockCorr = std::clamp(lockCorr * std::clamp(fbT / (float) realized, 0.8f, 1.25f), 0.25f, 4.f);

        // ── advance to the next cycle ──
        auto nx = (c + 1 < sched.size()) ? sched[c+1] : sched[c];
        float newTarget = nx.second ? fs / nx.first : 0.f;
        if (mode == SERVO_FIXED && (newTarget != cycleTarget || newTarget <= (float) N)) {
            lockCorr = 1.f; dur_err = 0.f;                   // reset-on-target-change / unreachable
        }
        if (nx.second) norm_k = (fs / nx.first) / sum * lockCorr;
        else { norm_k = 1.f; lockCorr = 1.f; }
        cycleTarget = newTarget;
    }
    if (finalCorr) *finalCorr = lockCorr;
    return out;
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

    // [5] LOCK servo, using only *production-reachable* state. Durations always come
    //     from reflect() into [bDurMin, bDurMax] with bDurMin >= 1, so every segment is
    //     >= 1 sample (unlike a synthetic sub-1 duration). Near the floor the error-
    //     diffused feed-forward is already close, so the servo's steady-state job is a
    //     ~1-sample trim; its important guarantees are that it (a) trims that residual,
    //     (b) saturates gracefully below the physical floor, and — the reason it exists
    //     as feedback rather than a one-shot — (c) does NOT double-correct when the
    //     target or LOCK state changes, and (d) can recall its converged value on load.
    {
        const int N = 64;

        // A production-valid clustered near-floor shape. At a 65-sample target the real
        // barrier range is [1, 2*65/64] = [1, 2.03125] (bDurMin=1, bDurMax=2*durCenter),
        // so use 57 minimum-length (1.0) + 7 near-max (2.03125) segments: the feed-forward
        // period realizes 66 samples and the servo must trim it to 65 (a real 1-sample
        // correction, not the vacuous 58/6 @2.1875 cluster whose feed-forward is already 65).
        float dNear[64];
        for (int i = 0; i < N; i++) dNear[i] = (i < 57) ? 1.0f : 2.03125f;

        bool ok = true;
        for (float fs : sampleRates) {
            float floorHz = fs / N;

            // (a) reachable near-floor: servo trims the residual to within a sample, and it
            //     is a real (not vacuous) correction — feed-forward alone is ~1 sample flat.
            {
                float f = floorHz / (65.f / N);            // target period 65 samples (> floor)
                std::vector<std::pair<float,bool>> hold(300, {f, true});
                long servo = servoRun(dNear, N, fs, hold).back();
                long ff    = servoRun(dNear, N, fs, {{f, true}}).front();   // feed-forward only (1 cycle)
                double servoErr = std::fabs(servo - fs / f), ffErr = std::fabs(ff - fs / f);
                if (servoErr > 1.1) { ok = false;
                    printf("  [5a] FAIL fs=%.0f near-floor servo period err=%.2f samples\n", fs, servoErr); }
                else if (ffErr < 0.9) { ok = false;        // fixture must actually exercise the servo
                    printf("  [5a] FAIL fs=%.0f feed-forward already exact (%.2f) — fixture vacuous\n", fs, ffErr); }
            }
            // (b) below the floor (target period < N): unreachable, must saturate at fs/N.
            {
                float f = floorHz * 1.3f;                  // asks for faster than fs/N
                std::vector<std::pair<float,bool>> hold(300, {f, true});
                long realized = servoRun(dNear, N, fs, hold).back();
                if (std::fabs(fs / realized - floorHz) > 1.0) { ok = false;
                    printf("  [5b] FAIL fs=%.0f below-floor -> %.2f Hz, expected fs/N=%.2f\n",
                           fs, fs / realized, floorHz); }
            }
        }

        // (c) target change must not double-correct. Well above the floor (durCenter ~11
        //     samples, period 704) so the regression is unmistakable: after converging at
        //     704, an octave up to 352 must land on 352, not 0.8*352 (the pre-fix servo
        //     fed the new target back against the old cycle → ~25% sharp, ~386 cents).
        {
            const float fs = 44100.f;
            float dHi[64]; for (int i = 0; i < N; i++) dHi[i] = 11.f;   // sum 704 = target
            float cf1 = fs / 704.f, cf2 = fs / 352.f;
            std::vector<std::pair<float,bool>> sched(80, {cf1, true});
            sched.push_back({cf2, true});                  // octave up on the 81st cycle
            long fixedPost = servoRun(dHi, N, fs, sched, SERVO_FIXED).at(80);
            long buggyPost = servoRun(dHi, N, fs, sched, BUG_FEEDBACK).at(80);
            double fixedC = cents(fs / fixedPost, fs / 352.f), buggyC = cents(fs / buggyPost, fs / 352.f);
            if (std::fabs(fixedC) > 8.0) { ok = false;
                printf("  [5c] FAIL octave change double-corrects: %.1f cents (want ~0)\n", fixedC); }
            if (std::fabs(buggyC) < 200.0) { ok = false;   // discrimination: the bug really was large
                printf("  [5c] FAIL discrimination: pre-fix servo only %.1f cents off — test vacuous\n", buggyC); }

            // (c') enabling LOCK on a cycle whose unlocked period differs from the target
            //      is the same hazard: unlocked natural period 704, LOCK on at target 352.
            std::vector<std::pair<float,bool>> eng(5, {cf1, false});
            for (int i = 0; i < 5; i++) eng.push_back({cf2, true});
            long fixedEng = servoRun(dHi, N, fs, eng, SERVO_FIXED).at(5);
            long buggyEng = servoRun(dHi, N, fs, eng, BUG_FEEDBACK).at(5);
            if (std::fabs(cents(fs / fixedEng, fs / 352.f)) > 8.0) { ok = false;
                printf("  [5c'] FAIL LOCK-engage double-corrects: %.1f cents\n", cents(fs / fixedEng, fs / 352.f)); }
            if (std::fabs(cents(fs / buggyEng, fs / 352.f)) < 200.0) { ok = false;
                printf("  [5c'] FAIL discrimination: pre-fix engage only %.1f cents off\n", cents(fs / buggyEng, fs / 352.f)); }
        }

        // (d) restore recalls the converged correction: a reloaded near-floor shape that
        //     saved its lockCorr lands on pitch immediately; a fresh reload (lockCorr=1)
        //     shows the one-sample transient the recall removes.
        {
            const float fs = 44100.f;
            float f = (fs / N) / (65.f / N);               // target 65 samples
            float converged = 1.f;
            std::vector<std::pair<float,bool>> hold(300, {f, true});
            servoRun(dNear, N, fs, hold, SERVO_FIXED, 1.f, &converged);
            long recalled = servoRun(dNear, N, fs, {{f, true}}, SERVO_FIXED, converged).front();
            long fresh    = servoRun(dNear, N, fs, {{f, true}}, SERVO_FIXED, 1.f).front();
            if (std::fabs(recalled - fs / f) > 1.1) { ok = false;
                printf("  [5d] FAIL restored lockCorr off by %.2f samples\n", std::fabs(recalled - fs / f)); }
            if (std::fabs(recalled - fs / f) > std::fabs(fresh - fs / f) + 1e-3) { ok = false;
                printf("  [5d] FAIL recall no better than a fresh reload\n"); }
        }

        // (e) a below-floor request winds no stale correction into a later reachable one.
        //     Hold LOCK at 5 kHz (target ~8.8 samples, far below the 64 floor) so the
        //     pre-fix servo would wind lockCorr to its 0.25 bound, then drop to C4: the
        //     first C4 cycle must land near C4, not stay floored at 64 samples (~+1677
        //     cents). Anti-windup + reset-on-target-change is what makes the fixed path
        //     land close; BUG_STALE (no guard) reproduces the blowup for discrimination.
        {
            const float fs = 44100.f;
            float dHi[64]; for (int i = 0; i < N; i++) dHi[i] = 11.f;
            float cf5k = 5000.f, cfC4 = fs / 168.56f;      // below floor, then a reachable C4
            std::vector<std::pair<float,bool>> sched(40, {cf5k, true});
            for (int i = 0; i < 8; i++) sched.push_back({cfC4, true});
            long fixedP = servoRun(dHi, N, fs, sched, SERVO_FIXED).at(40);
            double fixedC = cents(fs / fixedP, 261.626);
            double staleC = cents(fs / servoRun(dHi, N, fs, sched, BUG_STALE).at(40), 261.626);
            // reset-on-target-change clears both lockCorr AND dur_err, so the first C4
            // cycle lands within ~1 sample of 168.56 (≈±10 cents at this period).
            if (std::fabs(fixedP - 168.56) > 1.1) { ok = false;
                printf("  [5e] FAIL below-floor->C4 first cycle %ld samples (%.0f cents), want within 1 of 168.56\n", fixedP, fixedC); }
            if (std::fabs(staleC) < 1000.0) { ok = false;  // discrimination: the stale bug was huge
                printf("  [5e] FAIL discrimination: stale-correction only %.0f cents — test vacuous\n", staleC); }
        }

        printf("[5] LOCK servo: near-floor trim, below-floor saturation, no double-correct or stale-correction on target/LOCK change, restore recall: %s\n",
               ok ? "PASS" : "FAIL");
        if (!ok) fails++;
    }

    printf(fails ? "FAIL: %d GENDYN invariant(s) broken\n" : "PASS: GENDYN pitch/playback invariants hold\n", fails);
    return fails ? 1 : 0;
}
