// Standalone stability + calibration suite for Bunnies' predator–prey kernel.
// Calls the SDK-free production core while sweeping BALANCE × WILD × KICK in both modes at several
// pitches over long runs; asserts x,y stay finite, positive, and bounded. Reports
// the LV default period (→ RATE_CAL), the RM Hopf coverage, and LV servo health
// (max |V−V0| after settling, servo-clamp fraction, positivity-floor contacts).
// Servo and numerical policy come from the same production core used by Rack.
//
//   g++ -O2 -o /tmp/t tools/stability/bunnies.cpp && /tmp/t    (exit 0 = pass)
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include "../../src/dsp/bunnies_core.hpp"

namespace model = coalescent::bunnies;

struct Diag { float amp = 0, period = 0, minState = 0, maxVerr = 0, clampFrac = 0; };

// Shared-core event driver for the default LV orbit. `perSubstep` is the
// production policy; false is a deliberate sampled-boundary mutant that proves
// the regression distinguishes missed peaks.
static void gateContract(float dtau, int NS, bool perSubstep, int counts[2],
                         std::vector<int>& order) {
    const float balance = 0.5f, wild = 0.4f;
    const float gamma = model::gammaFromBalance(balance);
    const float V0 = model::targetInvariant(gamma, wild);
    const coalescent::OdeStepPlan plan = model::stepPlan(dtau, model::LV, gamma);
    const model::Parameters parameters = {model::LV, gamma, 0.f, 0.f, 0.f};
    float y[2];
    model::seed(1.f, 1.f, y);
    float previous[2] = {};
    bool rising[2] = {};
    counts[0] = counts[1] = 0;
    order.clear();
    auto observe = [&](int lane, float centered, int sample) {
        if (model::peakStep(centered, previous[lane], rising[lane]) && sample > NS / 3) {
            ++counts[lane];
            if (order.size() < 24) order.push_back(lane);
        }
    };
    for (int sample = 0; sample < NS; ++sample) {
        for (int k = 0; k < plan.count; ++k) {
            model::step(y, plan.h, parameters);
            if (!model::repairState(y)) return;
            if (perSubstep) {
                observe(0, y[0] - 1.f, sample);
                observe(1, y[1] - 1.f, sample);
            }
        }
        model::applyLvServo(y, gamma, V0, plan.delta);

        const float cx = y[0] - 1.f, cy = y[1] - 1.f;
        if (std::fabs(cx) < model::POP_MIN && std::fabs(cy) < model::POP_MIN)
            rising[0] = rising[1] = false;
        if (!perSubstep) {
            observe(0, cx, sample);
            observe(1, cy, sample);
        }
        else {
            // The LV amplitude servo is a numerical drift correction, not part of
            // the ecological trajectory. Synchronize the detector to its corrected
            // endpoint without manufacturing a second peak at the sample boundary.
            previous[0] = cx;
            previous[1] = cy;
        }
    }
}

// returns finite+positive+bounded ok; fills Diag (amp/period 2nd half; LV servo health).
static bool run(int mode, float balance, float wild, float kick, float dtau, int NS, Diag& dg) {
    float gamma = 0, V0 = 0, K = 0, c = 0;
    if (mode == model::LV) {
        gamma = model::gammaFromBalance(balance);
        V0 = model::targetInvariant(gamma, wild);
    }
    else {
        c = model::mortalityFromBalance(balance);
        K = model::capacityFromWild(wild);
    }
    float cx = 1.f, cy = 1.f;
    if (mode == model::RM && !model::rmCenter(K, c, cx, cy)) return true;
    const coalescent::OdeStepPlan plan = model::stepPlan(dtau, mode, gamma);
    const model::Parameters parameters = {mode, gamma, K, c, kick * model::KICK_GAIN};
    float y[2];
    model::seed(cx, cy, y);
    dg.minState = 1e9f; float pmn = 1e9f, pmx = -1e9f, prev = 0, lc = -1; double sp = 0; int np = 0;
    long clamps = 0, tot = 0; float maxVe = 0;
    for (int s = 0; s < NS; s++) {
        for (int k = 0; k < plan.count; k++) {
            model::step(y, plan.h, parameters);
            if (!model::repairState(y)) return false;
        }
        if (mode == model::LV) {
            const model::ServoResult servo = model::applyLvServo(y, gamma, V0, plan.delta);
            if (servo.appliedX != servo.rawX || servo.appliedY != servo.rawY) clamps++;
            tot++;
            if (s > NS / 2)
                maxVe = std::max(maxVe, std::fabs(servo.invariant - V0));
        }
        for (int j = 0; j < 2; j++)
            if (!std::isfinite(y[j]) || y[j] <= 0.f || y[j] > model::STATE_MAX) return false;
        dg.minState = std::min(dg.minState, std::min(y[0], y[1]));
        if (s > NS / 2) { pmn = std::min(pmn, y[0]); pmx = std::max(pmx, y[0]);
            float cc = y[0] - cx, t = s * plan.delta; if (prev <= 0 && cc > 0) { if (lc > 0) { sp += t - lc; np++; } lc = t; } prev = cc; }
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
    float gamma = model::gammaFromBalance(0.5f);
    float ratecalMeasured = dg.period * std::sqrt(gamma);
    printf("LV default (bal=0.5 wild=0.4): amp=%.2f period=%.2f tau => RATE_CAL~%.2f (const %.2f)\n",
           dg.amp, dg.period, ratecalMeasured, model::RATE_CAL);
    // Assert the constant tracks the measured value: |Δ| in cents must stay small so
    // the LV default voicing lands on C4. (7.49 vs 7.33 was ~37 cents sharp.)
    float cents = 1200.f * std::log2(model::RATE_CAL / ratecalMeasured);
    if (std::fabs(cents) > 5.f) {
        printf("FAIL: RATE_CAL %.2f is %.1f cents off the measured %.2f (LV default not on C4)\n",
               model::RATE_CAL, cents, ratecalMeasured);
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
        bool cvOk = coalescent::finiteOr(NAN, 0.f) == 0.f
                 && coalescent::finiteOr(INFINITY, 0.f) == 0.f
                 && coalescent::finiteClamp(INFINITY, 0.f, -15.f, 15.f) == 0.f
                 && coalescent::finiteClamp(1e30f, 0.f, -15.f, 15.f) == 15.f;
        printf("hostile-CV neutralization/bounding: %s\n", cvOk ? "PASS" : "FAIL");
        if (!cvOk) ++fails;
    }

    if (fails) { printf("FAIL: %d runs left the positive-bounded region\n", fails); return 1; }
    if (worstClamp > 60.f) { printf("WARN: LV servo clamps %.0f%% — consider lowering STAB_K\n", worstClamp); }
    if (worstFloor <= model::POS_FLOOR * 1.5f) { printf("WARN: LV orbit reaches the positivity floor at high WILD\n"); }
    printf("PASS: predator-prey finite, positive, bounded in both modes (KICK-stressed)\n");
    return 0;
}
