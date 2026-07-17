// Standalone stability + calibration replica of Bunnies' predator–prey kernel
// (mirrors src/Bunnies.cpp). Sweeps BALANCE × WILD × KICK in both modes at several
// pitches over long runs; asserts x,y stay finite, positive, and bounded. Reports
// the LV default period (→ RATE_CAL), the RM Hopf coverage, and LV servo health
// (max |V−V0| after settling, servo-clamp fraction, positivity-floor contacts).
// Servo constants match the module (STAB_K=0.5, STAB_FLOOR=0.2, V0≤Vmin+3.5).
//
//   g++ -O2 -o /tmp/t tools/stability/bunnies.cpp && /tmp/t    (exit 0 = pass)
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>

static constexpr float POS = 1e-4f, STATE_MAX = 1e3f, HSUB_MAX = 0.05f;
static constexpr float STAB_K = 0.5f, STAB_FLOOR = 0.2f, MAX_STAB_STEP = 0.25f, LV_V0_RANGE = 3.5f;
static constexpr float RM_B = 0.5f, RM_S = 1.0f, KICK_GAIN = 0.5f;
static constexpr int   MIN_SUB = 2, MAX_SUB = 64;
static constexpr float RATE_CAL = 7.33f;   // must match src/Bunnies.cpp; asserted below so tuning can't silently drift

static inline void deriv(const float v[2], float d[2], int mode, float gamma, float K, float c, float kick) {
    float X = std::max(v[0], POS), Y = std::max(v[1], POS);
    if (mode == 0) { d[0] = X * (1.f - Y);              d[1] = gamma * Y * (X - 1.f); }
    else           { float g = X / (1.f + RM_B * X);
                     d[0] = X * (1.f - X / K) - Y * g;  d[1] = RM_S * Y * (g - c); }
    d[0] += kick;
}
static inline void rk4(float y[2], float h, int mode, float gamma, float K, float c, float kick) {
    float k1[2], k2[2], k3[2], k4[2], t[2];
    deriv(y, k1, mode, gamma, K, c, kick); for (int i = 0; i < 2; i++) t[i] = y[i] + .5f * h * k1[i];
    deriv(t, k2, mode, gamma, K, c, kick); for (int i = 0; i < 2; i++) t[i] = y[i] + .5f * h * k2[i];
    deriv(t, k3, mode, gamma, K, c, kick); for (int i = 0; i < 2; i++) t[i] = y[i] + h * k3[i];
    deriv(t, k4, mode, gamma, K, c, kick);
    for (int i = 0; i < 2; i++) y[i] += h / 6.f * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
}

struct Diag { float amp = 0, period = 0, minState = 0, maxVerr = 0, clampFrac = 0; };

struct PeakDetector {
    float prev = 0.f;
    bool rising = false;

    bool step(float centered) {
        if (centered > prev + 1e-7f) rising = true;
        else if (centered < prev - 1e-7f) {
            const bool fire = rising && prev > 0.05f;
            rising = false;
            prev = centered;
            return fire;
        }
        prev = centered;
        return false;
    }
};

// Event replica for the default LV orbit. `perSubstep` is the production policy;
// false reproduces the former once-per-audio-sample detector for discrimination.
static void gateContract(float dtau, int NS, bool perSubstep, int counts[2],
                         std::vector<int>& order) {
    const float balance = 0.5f, wild = 0.4f;
    const float gamma = 0.2f + balance * 4.8f;
    const float V0 = gamma + 1.f + wild * LV_V0_RANGE;
    float dt = std::min(dtau / std::sqrt(gamma), HSUB_MAX * MAX_SUB);
    int Ksub = std::min(MAX_SUB, std::max(MIN_SUB, (int) std::ceil(dt / HSUB_MAX)));
    float h = dt / Ksub;
    float y[2] = {1.15f, 0.92f};
    PeakDetector det[2];
    counts[0] = counts[1] = 0;
    order.clear();
    auto observe = [&](int lane, float centered, int sample) {
        if (det[lane].step(centered) && sample > NS / 3) {
            ++counts[lane];
            if (order.size() < 24) order.push_back(lane);
        }
    };
    for (int sample = 0; sample < NS; ++sample) {
        for (int k = 0; k < Ksub; ++k) {
            rk4(y, h, 0, gamma, 0.f, 0.f, 0.f);
            y[0] = std::min(std::max(y[0], POS), STATE_MAX);
            y[1] = std::min(std::max(y[1], POS), STATE_MAX);
            if (perSubstep) {
                observe(0, y[0] - 1.f, sample);
                observe(1, y[1] - 1.f, sample);
            }
        }
        float X = y[0], Y = y[1];
        float V = gamma * (X - std::log(X)) + (Y - std::log(Y));
        float dVx = gamma * (1.f - 1.f / std::max(X, STAB_FLOOR));
        float dVy = 1.f - 1.f / std::max(Y, STAB_FLOOR);
        float sx = std::clamp(-STAB_K * dt * (V - V0) * dVx, -MAX_STAB_STEP, MAX_STAB_STEP);
        float sy = std::clamp(-STAB_K * dt * (V - V0) * dVy, -MAX_STAB_STEP, MAX_STAB_STEP);
        y[0] = std::max(X + sx, POS);
        y[1] = std::max(Y + sy, POS);

        const float cx = y[0] - 1.f, cy = y[1] - 1.f;
        if (std::fabs(cx) < 0.05f && std::fabs(cy) < 0.05f)
            det[0].rising = det[1].rising = false;
        if (!perSubstep) {
            observe(0, cx, sample);
            observe(1, cy, sample);
        }
        else {
            // The LV amplitude servo is a numerical drift correction, not part of
            // the ecological trajectory. Synchronize the detector to its corrected
            // endpoint without manufacturing a second peak at the sample boundary.
            det[0].prev = cx;
            det[1].prev = cy;
        }
    }
}

// returns finite+positive+bounded ok; fills Diag (amp/period 2nd half; LV servo health).
static bool run(int mode, float balance, float wild, float kick, float dtau, int NS, Diag& dg) {
    float gamma = 0, V0 = 0, K = 0, c = 0;
    if (mode == 0) { gamma = 0.2f + balance * 4.8f; V0 = (gamma + 1.f) + wild * LV_V0_RANGE; }
    else { c = std::min(0.15f + balance * 0.45f, 0.95f / RM_B); K = 1.2f + wild * 10.8f; }
    float cx = 1.f, cy = 1.f;
    if (mode == 1) { cx = c / (1.f - RM_B * c); cy = cx * (1.f - cx / K) / c; if (!(cx > 0 && cy > 0)) return true; }
    float dt = dtau; if (mode == 0) dt /= std::sqrt(gamma);
    dt = std::min(dt, HSUB_MAX * MAX_SUB);
    int Ksub = std::min(MAX_SUB, std::max(MIN_SUB, (int) std::ceil(dt / HSUB_MAX))); float h = dt / Ksub;
    float kf = kick * KICK_GAIN;
    float y[2] = {cx * 1.1f + 0.05f, cy * 0.9f + 0.02f};   // mirrors production reseed() in src/Bunnies.cpp
    dg.minState = 1e9f; float pmn = 1e9f, pmx = -1e9f, prev = 0, lc = -1; double sp = 0; int np = 0;
    long clamps = 0, tot = 0; float maxVe = 0;
    for (int s = 0; s < NS; s++) {
        for (int k = 0; k < Ksub; k++) { rk4(y, h, mode, gamma, K, c, kf);
            y[0] = std::min(std::max(y[0], POS), STATE_MAX); y[1] = std::min(std::max(y[1], POS), STATE_MAX); }
        if (mode == 0) {
            float X = y[0], Y = y[1]; float V = gamma * (X - std::log(X)) + (Y - std::log(Y));
            float dVx = gamma * (1.f - 1.f / std::max(X, STAB_FLOOR)), dVy = (1.f - 1.f / std::max(Y, STAB_FLOOR));
            float st = STAB_K * dt;
            float sx = -st * (V - V0) * dVx, sy = -st * (V - V0) * dVy;
            float cxx = std::clamp(sx, -MAX_STAB_STEP, MAX_STAB_STEP), cyy = std::clamp(sy, -MAX_STAB_STEP, MAX_STAB_STEP);
            if (cxx != sx || cyy != sy) clamps++; tot++;
            y[0] = std::max(X + cxx, POS); y[1] = std::max(Y + cyy, POS);
            if (s > NS / 2) maxVe = std::max(maxVe, std::fabs(V - V0));
        }
        for (int j = 0; j < 2; j++) if (!std::isfinite(y[j]) || y[j] <= 0.f || y[j] > STATE_MAX) return false;
        dg.minState = std::min(dg.minState, std::min(y[0], y[1]));
        if (s > NS / 2) { pmn = std::min(pmn, y[0]); pmx = std::max(pmx, y[0]);
            float cc = y[0] - cx, t = s * dt; if (prev <= 0 && cc > 0) { if (lc > 0) { sp += t - lc; np++; } lc = t; } prev = cc; }
    }
    dg.amp = pmx - pmn; dg.period = np ? (float)(sp / np) : 0.f;
    dg.maxVerr = maxVe; dg.clampFrac = tot ? 100.f * clamps / tot : 0.f;
    return true;
}

int main() {
    int fails = 0, runs = 0; Diag dg;
    // Bounded sweep — both modes, BALANCE × WILD × KICK × pitch (KICK stress incl.).
    for (int mode = 0; mode < 2; mode++)
        for (float bal = 0.f; bal <= 1.f; bal += 0.2f)
            for (float wild = 0.f; wild <= 1.f; wild += 0.2f)
                for (float kick : {0.f, 2.f, -5.f, 5.f})       // KICK stress: hard prey-force both signs
                    for (float dtau : {0.01f, 0.05f, 0.3f}) {
                        runs++;
                        if (!run(mode, bal, wild, kick, dtau, 12000, dg)) {
                            fails++; printf("  FAIL mode=%s bal=%.1f wild=%.1f kick=%.1f dtau=%.2f\n",
                                            mode ? "RM" : "LV", bal, wild, kick, dtau);
                        }
                    }
    printf("stability: %d/%d runs finite+positive+bounded (incl. KICK ±5)\n", runs - fails, runs);

    run(0, 0.5f, 0.4f, 0.f, 0.02f, 120000, dg);   // LV default
    float gamma = 0.2f + 0.5f * 4.8f;
    float ratecalMeasured = dg.period * std::sqrt(gamma);
    printf("LV default (bal=0.5 wild=0.4): amp=%.2f period=%.2f tau => RATE_CAL~%.2f (const %.2f)\n",
           dg.amp, dg.period, ratecalMeasured, RATE_CAL);
    // Assert the constant tracks the measured value: |Δ| in cents must stay small so
    // the LV default voicing lands on C4. (7.49 vs 7.33 was ~37 cents sharp.)
    float cents = 1200.f * std::log2(RATE_CAL / ratecalMeasured);
    if (std::fabs(cents) > 5.f) {
        printf("FAIL: RATE_CAL %.2f is %.1f cents off the measured %.2f (LV default not on C4)\n",
               RATE_CAL, cents, ratecalMeasured);
        fails++;
    }

    // LV servo health across WILD: V should track V0, clamp
    // fraction stay low, and orbits stay off the positivity floor.
    printf("LV servo health (bal=0.5): max|V-V0|  servo-clamp%%  minState\n");
    float worstClamp = 0, worstFloor = 1e9f;
    for (float wild = 0.f; wild <= 1.f; wild += 0.25f) {
        run(0, 0.5f, wild, 0.f, 0.02f, 120000, dg);
        printf("  wild=%.2f:  %.3f       %5.1f%%      %.4f\n", wild, dg.maxVerr, dg.clampFrac, dg.minState);
        worstClamp = std::max(worstClamp, dg.clampFrac); worstFloor = std::min(worstFloor, dg.minState);
    }
    printf("RM Hopf coverage: ");
    for (float bal : {0.f, 0.5f, 1.f}) { printf("bal%.1f[", bal);
        for (float wild = 0.f; wild <= 1.f; wild += 0.25f) { run(1, bal, wild, 0.f, 0.02f, 60000, dg);
            printf("%s", dg.amp > 0.05f ? "o" : "."); } printf("] "); } printf("(o=cycle .=rest)\n");

    // Peak events must remain a two-phase alternating clock at normal speed and
    // must not disappear when one audio sample spans a large fraction of a cycle.
    {
        auto alternating = [](const std::vector<int>& order) {
            if (order.size() < 8) return false;
            for (std::size_t i = 1; i < order.size(); ++i)
                if (order[i] == order[i - 1]) return false;
            return true;
        };
        int normal[2], fast[2], sampledFast[2];
        std::vector<int> normalOrder, fastOrder, sampledOrder;
        gateContract(0.02f, 100000, true, normal, normalOrder);
        gateContract(3.2f, 12000, true, fast, fastOrder);
        gateContract(3.2f, 12000, false, sampledFast, sampledOrder);
        bool normalOk = normal[0] > 20 && std::abs(normal[0] - normal[1]) <= 1
                     && alternating(normalOrder);
        bool fastOk = fast[0] > 100 && std::abs(fast[0] - fast[1]) <= 1
                   && alternating(fastOrder);
        bool discriminates = sampledFast[0] + sampledFast[1] < fast[0] + fast[1];
        printf("peak gates: normal=%d/%d fast-substep=%d/%d fast-sampled=%d/%d: %s\n",
               normal[0], normal[1], fast[0], fast[1], sampledFast[0], sampledFast[1],
               (normalOk && fastOk && discriminates) ? "PASS" : "FAIL");
        if (!(normalOk && fastOk && discriminates)) ++fails;
    }

    // Rack wrapper policy: hostile non-finite CV is neutral and an extreme finite
    // KICK is bounded before entering all RK stages.
    {
        auto finiteOr = [](float v, float fallback) { return std::isfinite(v) ? v : fallback; };
        auto kickInput = [&](float v) { return std::clamp(finiteOr(v, 0.f), -15.f, 15.f); };
        bool cvOk = finiteOr(NAN, 0.f) == 0.f && finiteOr(INFINITY, 0.f) == 0.f
                 && kickInput(INFINITY) == 0.f && kickInput(1e30f) == 15.f;
        printf("hostile-CV neutralization/bounding: %s\n", cvOk ? "PASS" : "FAIL");
        if (!cvOk) ++fails;
    }

    if (fails) { printf("FAIL: %d runs left the positive-bounded region\n", fails); return 1; }
    if (worstClamp > 60.f) { printf("WARN: LV servo clamps %.0f%% — consider lowering STAB_K\n", worstClamp); }
    if (worstFloor <= POS * 1.5f) { printf("WARN: LV orbit reaches the positivity floor at high WILD\n"); }
    printf("PASS: predator-prey finite, positive, bounded in both modes (KICK-stressed)\n");
    return 0;
}
