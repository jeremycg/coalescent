// Standalone stability + calibration suite for Operon's repressilator kernel.
// Calls the SDK-free production core (dimensionless Elowitz–Leibler, RK4 with
// pitch-adaptive substepping). Sweeps ALPHA × HILL × BETA × LEAK at several
// pitches over long runs and asserts all six states stay finite and bounded;
// also reports the default-voicing dimensionless period (→ RATE_CAL) and where
// the ALPHA dial crosses the Hopf threshold. No Rack SDK needed.
//
//   g++ -O2 -o /tmp/t tools/stability/operon.cpp && /tmp/t     (exit 0 = pass)
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include "../../src/dsp/operon_core.hpp"

namespace model = coalescent::operon;

// Independent pStar reference. A 40-iteration bisection gives a near-exact root
// to compare with the production core's bounded warm-Newton solver.
static float pStarBisect(float alpha, float n, float a0, int iters) {
    float lo = 0.f, hi = std::max(1.f, alpha + a0 + 1.f);
    for (int i = 0; i < iters; ++i) { float mid = 0.5f * (lo + hi);
        float g = mid - alpha / (1.f + std::pow(std::max(mid, 0.f), n)) - a0; if (g > 0.f) hi = mid; else lo = mid; }
    return 0.5f * (lo + hi);
}
using Par = model::Parameters;

static model::StateRepair advanceDirect(float (&state)[6], float h,
                                        const Par& parameters, float perturb = 0.f) {
    return model::advanceAcceptedSubstep(
        state, h, parameters, perturb, [&](float repressor) {
            return model::directHillResponse(repressor, parameters.n);
        });
}

// Gate contract from the production seed. Also counts the former one-step
// detector so this regression proves that it distinguishes the silent-gate bug.
static void gateContract(const Par& p, float dtau, int NS, int counts[3],
                         int legacyCounts[3], std::vector<int>& order) {
    float y[6];
    model::seed(y);
    const float center = pStarBisect(p.alpha, p.n, p.a0, 40);
    const coalescent::OdeStepPlan plan = model::stepPlan(dtau);
    bool armed[3] = {};
    float legacyPrev[3] = {y[3] - center, y[4] - center, y[5] - center};
    std::fill(counts, counts + 3, 0);
    std::fill(legacyCounts, legacyCounts + 3, 0);
    order.clear();
    for (int s = 0; s < NS; ++s) {
        for (int k = 0; k < plan.count; ++k) {
            if (!advanceDirect(y, plan.h, p))
                return;
            for (int lane = 0; lane < 3; ++lane) {
                const float c = y[3 + lane] - center;
                const bool fired = model::gateStep(c, armed[lane]);
                const bool legacyFired = legacyPrev[lane] <= model::GATE_LOW
                                      && c >= model::GATE_HIGH;
                legacyPrev[lane] = c;
                if (s > NS / 2 && fired) {
                    ++counts[lane];
                    if (order.size() < 24) order.push_back(lane);
                }
                if (s > NS / 2 && legacyFired) ++legacyCounts[lane];
            }
        }
    }
}
// Fast finite/bounded check — single pass, no measurement.
static bool bounded(const Par& p, float dtau, int NS) {
    float y[6];
    model::seed(y);
    const coalescent::OdeStepPlan plan = model::stepPlan(dtau);
    for (int s = 0; s < NS; ++s) {
        for (int k = 0; k < plan.count; ++k) {
            const model::StateRepair repair = advanceDirect(y, plan.h, p);
            if (!repair || repair.exceededRange) return false;
        }
        for (int j = 0; j < 6; ++j)
            if (!std::isfinite(y[j]) || std::fabs(y[j]) > model::STATE_MAX) return false;
    }
    return true;
}
// Centered-p0 amplitude + dimensionless period (for the two calibration voicings).
static void measure(const Par& p, float dtau, float& amp, float& periodTau) {
    float y[6];
    model::seed(y);
    const coalescent::OdeStepPlan plan = model::stepPlan(dtau);
    const int NS = 60000; double mean = 0; int cnt = 0;
    for (int s = 0; s < NS; ++s) {
        for (int k = 0; k < plan.count; ++k)
            if (!advanceDirect(y, plan.h, p)) return;
        if (s > NS / 2) { mean += y[3]; cnt++; }
    }
    mean /= std::max(cnt, 1);
    float mn = 1e9f, mx = -1e9f, prev = 0.f, lastCross = -1.f; double sumP = 0; int nP = 0;
    float y2[6];
    model::seed(y2);
    for (int s = 0; s < NS; ++s) {
        for (int k = 0; k < plan.count; ++k)
            if (!advanceDirect(y2, plan.h, p)) return;
        if (s > NS / 2) {
            float v = y2[3]; mn = std::min(mn, v); mx = std::max(mx, v);
            float c = v - (float) mean; float tcur = s * plan.count * plan.h;
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
    // production core's lookup, including the sharp n=8 knee).
    double lutMax = 0;
    for (float nn : {1.2f, 2.5f, 5.f, 8.f}) {
        model::HillLut lut;
        for (int sample = 0; sample < model::HillLut::settleSamples() + 64; ++sample)
            lut.process(nn);
        if (lut.direct(nn)) {
            printf("  Hill LUT failed to become active for n=%.1f\n", nn);
            ++fails;
            continue;
        }
        for (float rep = 0.f; rep < model::HillLut::xMax(); rep += 0.013f) {
            float lu = lut.lookup(rep);
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
                float got    = model::solvePStar(alpha, nn, a0, warm);
                pStarMax = std::max(pStarMax, (double) std::fabs(got - ref));
            }
    printf("pStar warm Newton vs accurate root: max |Δ| = %.2e\n", pStarMax);

    // Hill-LUT production state machine: a build starts only
    // after n has been within N_MOVE_EPS of an anchor for 2048 consecutive samples;
    // any move — fast OR slow cumulative drift — resets the window and drops the
    // partial build. Assert (a) a slow sub-threshold HILL LFO does NOT thrash the
    // 256-entry fill (the regression the incremental rebuild could otherwise cause,
    // aborting+restarting every sample), and (b) a genuinely static n still completes
    // a build so the LUT goes live.
    {
        auto gate = [&](float fs, float lfoHz, float amp, double secs, long& slices, long& completions) {
            model::HillLut lut;
            slices = 0; completions = 0;
            long ns = (long) (secs * fs);
            for (long s = 0; s < ns; s++) {
                float n = 4.5f + amp * std::sin(2.0 * M_PI * lfoHz * s / fs);   // centred in the HILL range
                const model::HillLut::BuildWork work = lut.process(n);
                if (work.entries > 0) ++slices;
                if (work.completed) ++completions;
            }
        };
        long slices, comps;
        bool gateOk = true;
        // (a) slow LFOs over 10 s: at most a few builds (turning-point settling), never
        //     the tens-of-thousands-of-slices-per-second thrash. Bound generously at the
        //     equivalent of ~20 full rebuilds (a full rebuild is ~33 slices).
        for (float fs : {48000.f, 192000.f})
            for (float hz : {0.1f, 0.25f, 1.0f}) {
                gate(fs, hz, 2.f, 10.0, slices, comps);
                if (slices > 20 * 33) { gateOk = false;
                    printf("  LUT gate THRASH: fs=%.0f lfo=%.2f -> %ld slices in 10 s\n", fs, hz, slices); }
            }
        // (a') an in-band micro-oscillation (±0.9·N_MOVE_EPS at 1 kHz — n never leaves the
        //      settle window) must build once and then let the LUT absorb the wobble, not
        //      rebuild every time n crosses the old build point. Bound at ~2 rebuilds.
        gate(48000.f, 1000.f, 0.9e-4f, 5.0, slices, comps);
        if (slices > 2 * 33) { gateOk = false;
            printf("  LUT gate micro-oscillation thrash: %ld slices in 5 s\n", slices); }
        // (b) a static n must still build the table.
        gate(48000.f, 0.0f, 0.f, 1.0, slices, comps);
        if (comps < 1) { gateOk = false; printf("  LUT gate: static n never completed a build\n"); }
        printf("Hill LUT gate: slow-LFO + micro-oscillation no-thrash, static completes: %s\n", gateOk ? "PASS" : "FAIL");
        if (!gateOk) ++fails;
    }

    // Exercise the exact LUT-backed callback and accepted-substep path used by
    // the Rack wrapper, not only the direct-pow reference path above.
    {
        model::HillLut lut;
        for (int sample = 0; sample < model::HillLut::settleSamples() + 64; ++sample)
            lut.process(def.n);
        float state[6];
        model::seed(state);
        const coalescent::OdeStepPlan plan = model::stepPlan(0.4f);
        bool lutRunOk = !lut.direct(def.n);
        for (int sample = 0; sample < 20000 && lutRunOk; ++sample) {
            lut.process(def.n);
            const bool direct = lut.direct(def.n);
            const auto response = [&](float repressor) {
                return lut.response(repressor, def.n, direct);
            };
            for (int substep = 0; substep < plan.count; ++substep) {
                const model::StateRepair repair = model::advanceAcceptedSubstep(
                    state, plan.h, def, 0.f, response);
                if (!repair || repair.exceededRange) {
                    lutRunOk = false;
                    break;
                }
            }
        }
        for (float value : state)
            lutRunOk = lutRunOk && std::isfinite(value)
                    && value >= 0.f && value <= model::STATE_MAX;
        printf("LUT-backed accepted-substep trajectory: %s\n",
               lutRunOk ? "PASS" : "FAIL");
        if (!lutRunOk) ++fails;
    }

    // Gate regression: default oscillation must produce three balanced lanes in
    // a stable cyclic order. The old detector required a single substep to jump
    // across the whole deadband and is expected to miss essentially everything.
    {
        int counts[3], legacy[3];
        std::vector<int> order;
        gateContract(def, 0.02f, 100000, counts, legacy, order);
        const int total = counts[0] + counts[1] + counts[2];
        const int legacyTotal = legacy[0] + legacy[1] + legacy[2];
        bool balanced = counts[0] > 20
                     && std::abs(counts[0] - counts[1]) <= 1
                     && std::abs(counts[1] - counts[2]) <= 1;
        bool cyclic = order.size() >= 12;
        int direction = 0;
        if (cyclic) {
            const int d = (order[1] - order[0] + 3) % 3;
            direction = d;
            cyclic = d == 1 || d == 2;
            for (std::size_t i = 1; cyclic && i < order.size(); ++i)
                cyclic = (order[i] - order[i - 1] + 3) % 3 == d;
        }
        const bool discriminates = legacyTotal * 4 < total;
        const char* orderName = direction == 1 ? "+1" : (direction == 2 ? "-1" : "invalid");
        printf("gate contract: armed=%d/%d/%d legacy=%d/%d/%d order=%s: %s\n",
               counts[0], counts[1], counts[2], legacy[0], legacy[1], legacy[2],
               orderName,
               (balanced && cyclic && discriminates) ? "PASS" : "FAIL");
        if (!(balanced && cyclic && discriminates)) ++fails;

        // Near-threshold chatter contract: reaching HIGH without first visiting
        // LOW is silent; one armed crossing fires once; hovering around HIGH and
        // returning only to the deadband cannot retrigger.
        bool armed = false;
        int nearFires = 0;
        for (float c : {0.051f, 0.049f, 0.052f, -0.049f, 0.051f,
                        -0.051f, -0.02f, 0.051f, 0.049f, 0.052f})
            if (model::gateStep(c, armed)) ++nearFires;
        printf("gate near-threshold hysteresis: fires=%d (want 1): %s\n",
               nearFires, nearFires == 1 ? "PASS" : "FAIL");
        if (nearFires != 1) ++fails;
    }

    // Hostile-CV policy used by the Rack wrapper: non-finite modulation is
    // neutral, while a finite value is bounded before entering the ODE/pitch map.
    {
        bool cvOk = coalescent::finiteOr(NAN, 0.f) == 0.f
                 && coalescent::finiteOr(INFINITY, 0.f) == 0.f
                 && coalescent::finiteOr(-INFINITY, 2.5f) == 2.5f
                 && coalescent::finiteOr(3.f, 0.f) == 3.f;
        printf("hostile-CV neutralization: %s\n", cvOk ? "PASS" : "FAIL");
        if (!cvOk) ++fails;

        float finite[6] = {-1.f, model::STATE_MAX + 1.f, 0.1f,
                           0.2f, 0.3f, 0.4f};
        const model::StateRepair finiteRepair = model::repairStateDetailed(finite);
        float corrupt[6];
        model::seed(corrupt);
        corrupt[2] = NAN;
        const model::StateRepair corruptRepair = model::repairStateDetailed(corrupt);
        if (!corruptRepair)
            model::seed(corrupt);  // Rack wrapper's atomic recovery action
        bool repairOk = finiteRepair && finiteRepair.clamped
                     && finiteRepair.exceededRange && finite[0] == 0.f
                     && finite[1] == model::STATE_MAX && !corruptRepair;
        for (float value : corrupt)
            repairOk = repairOk && std::isfinite(value)
                     && value >= 0.f && value <= model::STATE_MAX;
        printf("accepted-substep repair/reseed contract: %s\n",
               repairOk ? "PASS" : "FAIL");
        if (!repairOk) ++fails;
    }

    if (fails) { printf("FAIL: %d unbounded runs\n", fails); return 1; }
    if (lutMax > 1e-3) { printf("FAIL: Hill LUT error too large\n"); return 1; }
    if (pStarMax > 1e-4) { printf("FAIL: warm-started Newton pStar does not converge\n"); return 1; }
    printf("PASS: repressilator finite and bounded; Hill LUT within 1e-3 of pow; warm Newton pStar converges\n");
    return 0;
}
