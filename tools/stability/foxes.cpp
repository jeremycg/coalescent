// Standalone stability + calibration replica of Foxes' three-species food-chain
// kernel (mirrors src/Foxes.cpp). Implements the nondimensional Hastings–Powell
// (1991) tritrophic model and checks:
//   [1]  analytic-equilibrium residual over a dense BALANCE×WILD grid;
//   [2]  finite/positive/bounded over BALANCE×WILD×the full +/-15 V KICK range×pitch;
//        UNFORCED operation never touches the floor (negative KICK legitimately does);
//   [2b] a sweep from the FIXED production seed asserting no shaped output rails
//        (catches an alternate basin that the equilibrium-seeded [2] cannot see);
//   [3]  the canonical Hopf location (Routh–Hurwitz on the analytic Jacobian);
//   [4]  the default period → RATE_CAL within 5 cents;
//   [5]  the largest Lyapunov exponent by coupled state+tangent RK4 (≈0 at the
//        periodic default, robustly positive at canonical chaos) + fox-max spread;
//   [6]  a two-seed multistability scan near b1≈2.5 (reports coexisting attractors);
//   [2b] additionally reports FOX rail occupancy and fails only a truly dead output;
//   [7]  rest (WILD 0) settles near 0 V, and the default-attractor seed derivation;
//   [8]  input sanitization: NaN/inf V/OCT → rest and KICK bounded before the clamps;
//   [9]  timestep convergence: HSUB_MAX reproduces the fine-h period (≤5 cents) and
//        the chaotic fox-max spread (a statistical invariant);
//   [10] peak-gate contract: no chatter at rest, all fire when active, and the
//        default fires in grass→bunny→fox order;
//   [11] regime classification: the default is period-2 with a ~−20 dB subharmonic,
//        and the transition (WILD 0.57) is period-doubled beyond it but not chaos.
//   [12] the display's discarded first ring materially removes the WILD-to-rest
//        transient before the second-ring portrait is accepted.
//
// Constants below MUST match src/Foxes.cpp; RATE_CAL and the maps are asserted so
// calibration/mapping drift fails the build.
//
//   g++ -O2 -o /tmp/t tools/stability/foxes.cpp && /tmp/t    (exit 0 = pass)
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>

// ── canonical Hastings–Powell constants (published; do not tune) ──
static constexpr double A1 = 5.0, A2 = 0.1, D1 = 0.4, D2 = 0.01;

// ── musical maps (Coalescent product choices; asserted, not literature) ──
//   b1 = 1.0 + 5.2·wild²    (rest → cycle → period-doubling → chaos)
//   b2 = 1.75 + 0.5·balance (upper trophic saturation, narrowed to keep clear of
//                            the fox-extinction basin above b2≈2.35)
static inline double b1FromWild(double w) { return 1.0 + 5.2 * w * w; }
static inline double b2FromBal(double b)  { return 1.75 + 0.5 * b; }

// ── numerical policy (must match src/Foxes.cpp) ──
static constexpr float  POS_FLOOR = 1e-4f, STATE_MAX = 1e3f, HSUB_MAX = 0.1f;
static constexpr int    MIN_SUB = 2, MAX_SUB = 64;
static constexpr float  KICK_GAIN = 0.25f;
static constexpr int    TRAIL_N = 2048;
static constexpr double CAP_TAU = 0.225;
static constexpr float  RATE_CAL  = 61.387f;   // default period; asserted within 5 cents below

// Output shaping + the fixed production seed (must match src/Foxes.cpp) — used by
// the production-seed sweep and the output/rail contract checks.
static const double     FOX_GAIN[3] = {8.0, 18.0, 2.5};   // grass, bunny, fox
static constexpr double SEED_X = 0.812781, SEED_Y = 0.104883, SEED_Z = 12.478951;
static inline double    shaped(double centered, double gain) { return 5.0 * std::tanh(centered * gain); }

// Saturating functional responses. Templated so the exact float production kernel
// and the double-precision analysis share one definition.
template <typename T> static inline T f1(T x, T b1) { return T(A1) * x / (T(1) + b1 * x); }
template <typename T> static inline T f2(T y, T b2) { return T(A2) * y / (T(1) + b2 * y); }

// Derivative. KICK is a constant force into dx/dt across all RK stages of one
// sample; responses are evaluated on nonnegative local state so an intermediate
// negative RK stage cannot approach a rational-denominator singularity.
template <typename T>
static inline void deriv(const T v[3], T d[3], T b1, T b2, T kick) {
    T X = std::max(v[0], T(POS_FLOOR)), Y = std::max(v[1], T(POS_FLOOR)), Z = std::max(v[2], T(POS_FLOOR));
    T g1 = f1(X, b1), g2 = f2(Y, b2);
    d[0] = X * (T(1) - X) - g1 * Y + kick;
    d[1] = g1 * Y - g2 * Z - T(D1) * Y;
    d[2] = g2 * Z - T(D2) * Z;
}
template <typename T>
static inline void rk4(T y[3], T h, T b1, T b2, T kick) {
    T k1[3], k2[3], k3[3], k4[3], t[3];
    deriv(y, k1, b1, b2, kick); for (int i = 0; i < 3; i++) t[i] = y[i] + T(0.5) * h * k1[i];
    deriv(t, k2, b1, b2, kick); for (int i = 0; i < 3; i++) t[i] = y[i] + T(0.5) * h * k2[i];
    deriv(t, k3, b1, b2, kick); for (int i = 0; i < 3; i++) t[i] = y[i] + h * k3[i];
    deriv(t, k4, b1, b2, kick);
    for (int i = 0; i < 3; i++) y[i] += h / T(6) * (k1[i] + T(2) * k2[i] + T(2) * k3[i] + k4[i]);
}

// Analytic interior coexistence point. Returns false if any positivity condition
// fails (should never happen inside the exposed BALANCE×WILD box).
static bool equilibrium(double b1, double b2, double& xs, double& ys, double& zs) {
    double denom = A2 - b2 * D2;
    if (!(denom > 0)) return false;
    ys = D2 / denom;
    double disc = (b1 + 1) * (b1 + 1) - 4 * A1 * b1 * ys;
    if (!(disc >= 0)) return false;
    xs = ((b1 - 1) + std::sqrt(disc)) / (2 * b1);
    double fx = f1(xs, b1);
    zs = ys * (fx - D1) / D2;
    return xs > 0 && ys > 0 && fx > D1 && zs > 0;
}

// Jacobian at (x,y,z) for the b1,b2 system. Used by the Hopf and Lyapunov checks.
static void jacobian(double x, double y, double z, double b1, double b2, double J[3][3]) {
    double g1 = f1(x, b1), g2 = f2(y, b2);
    double g1p = A1 / ((1 + b1 * x) * (1 + b1 * x));
    double g2p = A2 / ((1 + b2 * y) * (1 + b2 * y));
    J[0][0] = 1 - 2 * x - g1p * y; J[0][1] = -g1;                 J[0][2] = 0;
    J[1][0] = g1p * y;             J[1][1] = g1 - g2p * z - D1;   J[1][2] = -g2;
    J[2][0] = 0;                   J[2][1] = g2p * z;             J[2][2] = g2 - D2;
}

struct Diag { float amp = 0, period = 0, minState = 0, maxState = 0; float floorFrac = 0, ceilFrac = 0; };

// Exact float production replica: adaptive substep, floor/ceiling backstop after
// each substep, KICK into dx/dt. Fills Diag over the 2nd half; returns
// finite+positive+bounded. Starts near the param's own equilibrium (perturbed).
static bool run(double bal, double wild, float kick, float dtau, int NS, Diag& dg) {
    float b1 = (float) b1FromWild(wild), b2 = (float) b2FromBal(bal);
    double xs, ys, zs;
    if (!equilibrium(b1FromWild(wild), b2FromBal(bal), xs, ys, zs)) return true; // guarded elsewhere
    float dt = std::min(dtau, HSUB_MAX * MAX_SUB);
    int Ksub = std::min(MAX_SUB, std::max(MIN_SUB, (int) std::ceil(dt / HSUB_MAX)));
    float h = dt / Ksub, kf = kick * KICK_GAIN;
    float y[3] = {(float)(xs * 1.2), (float) ys, (float) zs};
    dg.minState = 1e9f; dg.maxState = -1e9f;
    float pmn = 1e9f, pmx = -1e9f, prev = 0, lc = -1; double sp = 0; int np = 0;
    long floorHits = 0, ceilHits = 0, tot = 0;
    for (int s = 0; s < NS; s++) {
        for (int k = 0; k < Ksub; k++) {
            rk4<float>(y, h, b1, b2, kf);
            for (int j = 0; j < 3; j++) {
                if (!std::isfinite(y[j])) return false;           // check BEFORE clamp: don't let it mask NaN
                if (y[j] <= POS_FLOOR * 1.5f) floorHits++;
                if (y[j] >= STATE_MAX * 0.99f) ceilHits++;
                tot++;
                y[j] = std::min(std::max(y[j], POS_FLOOR), STATE_MAX);
            }
        }
        for (int j = 0; j < 3; j++) if (!std::isfinite(y[j]) || y[j] <= 0.f || y[j] > STATE_MAX) return false;
        dg.minState = std::min(dg.minState, std::min({y[0], y[1], y[2]}));
        dg.maxState = std::max(dg.maxState, std::max({y[0], y[1], y[2]}));
        // Period from upward x−xStar crossings. Track prev EVERY sample (so the
        // first 2nd-half crossing measures from a real phase, not a stale one) and
        // linearly interpolate the sub-sample crossing time (outer samples are a
        // coarse 0.02 tau; without interpolation a low-slope crossing quantizes to
        // ±1 sample ≈ ±1 tau and only averages out over even counts).
        float cc = y[0] - (float) xs;
        bool up = (prev <= 0.f && cc > 0.f);
        if (s > NS / 2) {
            pmn = std::min(pmn, y[0]); pmx = std::max(pmx, y[0]);
            if (up) { float frac = cc / (cc - prev), t = (s - frac) * dt;
                      if (lc > 0) { sp += t - lc; np++; } lc = t; }
        }
        prev = cc;
    }
    dg.amp = pmx - pmn; dg.period = np ? (float)(sp / np) : 0.f;
    dg.floorFrac = tot ? 100.f * floorHits / tot : 0.f;
    dg.ceilFrac  = tot ? 100.f * ceilHits  / tot : 0.f;
    return true;
}

// Largest Lyapunov exponent (double precision, offline). Co-evolves state and one
// tangent vector as a single coupled 6-D RK4 — the variational equation dδ/dt =
// J(state)·δ with J sampled at every RK stage's state, not frozen across the step —
// renormalizing every RENORM tau and accumulating log-growth after a burn-in.
static double lyapunov(double b1, double b2, double h, double burn, double span) {
    double xs, ys, zs;
    if (!equilibrium(b1, b2, xs, ys, zs)) return NAN;
    double y[3] = {xs * 1.2, ys, zs};
    for (double t = 0; t < burn; t += h) rk4<double>(y, h, b1, b2, 0.0);

    // coupled derivative of s = [state(0..2), tangent(3..5)]
    auto fc = [&](const double s[6], double o[6]) {
        deriv<double>(s, o, b1, b2, 0.0);                  // state part reuses the production derivative
        double J[3][3]; jacobian(s[0], s[1], s[2], b1, b2, J);
        for (int i = 0; i < 3; i++) o[3 + i] = J[i][0]*s[3] + J[i][1]*s[4] + J[i][2]*s[5];
    };
    double s[6] = {y[0], y[1], y[2], 1e-8, 0, 0}, sum = 0; int rn = 0;
    const double RENORM = 1.0; double acc = 0;
    for (double t = 0; t < span; t += h) {
        double k1[6], k2[6], k3[6], k4[6], tmp[6];
        fc(s, k1);   for (int i = 0; i < 6; i++) tmp[i] = s[i] + 0.5*h*k1[i];
        fc(tmp, k2); for (int i = 0; i < 6; i++) tmp[i] = s[i] + 0.5*h*k2[i];
        fc(tmp, k3); for (int i = 0; i < 6; i++) tmp[i] = s[i] + h*k3[i];
        fc(tmp, k4); for (int i = 0; i < 6; i++) s[i] += h/6.0*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
        acc += h;
        if (acc >= RENORM) {
            double n = std::sqrt(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
            sum += std::log(n / 1e-8); rn++;
            for (int i = 3; i < 6; i++) s[i] *= 1e-8 / n;   // renormalize the tangent to seed length
            acc = 0;
        }
    }
    return rn ? sum / (rn * RENORM) : NAN;
}

// Run the exact float kernel from the FIXED production seed (as the module does on
// load). Fills, per shaped output: mean, rail occupancy (fraction of time within
// 0.1 V of a ±5 V rail), and peak-to-peak. The main sweep seeds near each setting's
// own equilibrium, which HIDES an alternate basin; this exposes a BALANCE/WILD
// window that sticks the advertised FOX voice at a rail. A high mean-negative FOX is
// character (a burst-like signal with real positive excursions); a truly DEAD output
// is one that is railed almost always AND has negligible peak-to-peak.
// Optionally injects a constant KICK (already gain-scaled).
static void runFromSeed(double bal, double wild, float kick,
                        double meanV[3], double railFrac[3], double p2p[3]) {
    float b1 = (float) b1FromWild(wild), b2 = (float) b2FromBal(bal);
    double xs, ys, zs; equilibrium(b1FromWild(wild), b2FromBal(bal), xs, ys, zs);
    int Ksub = std::min(MAX_SUB, std::max(MIN_SUB, (int) std::ceil(0.025f / HSUB_MAX)));
    float h = 0.025f / Ksub;
    float y[3] = {(float) SEED_X, (float) SEED_Y, (float) SEED_Z};
    auto step = [&]() { for (int k = 0; k < Ksub; k++) { rk4<float>(y, h, b1, b2, kick);
        for (int j = 0; j < 3; j++) y[j] = std::min(std::max(y[j], POS_FLOOR), STATE_MAX); } };
    for (int s = 0; s < 160000; s++) step();                 // burn in past the transient (~4000 tau)
    double sum[3] = {0, 0, 0}, lo[3] = {9, 9, 9}, hi[3] = {-9, -9, -9}; long railed[3] = {0, 0, 0}, n = 0;
    const double cen[3] = {xs, ys, zs};
    for (int s = 0; s < 100000; s++) { step();
        for (int i = 0; i < 3; i++) {
            double v = shaped(y[i] - cen[i], FOX_GAIN[i]);
            sum[i] += v; lo[i] = std::min(lo[i], v); hi[i] = std::max(hi[i], v);
            if (std::fabs(v) > 4.9) railed[i]++;
        }
        n++;
    }
    for (int i = 0; i < 3; i++) { meanV[i] = sum[i] / n; railFrac[i] = (double) railed[i] / n; p2p[i] = hi[i] - lo[i]; }
}

// Peak-gate contract from the production seed. Mirrors src/Foxes.cpp firePeak
// exactly (rising→falling local max above POP_GATE on the gained centered value,
// evaluated per substep) and counts the events per species over a fixed window.
static const double POP_GATE = 0.15;   // mirrors src/Foxes.cpp
static void gateCounts(double bal, double wild, int counts[3]) {
    float b1 = (float) b1FromWild(wild), b2 = (float) b2FromBal(bal);
    double xs, ys, zs; equilibrium(b1FromWild(wild), b2FromBal(bal), xs, ys, zs);
    const double cen[3] = {xs, ys, zs};
    int Ksub = std::min(MAX_SUB, std::max(MIN_SUB, (int) std::ceil(0.02f / HSUB_MAX)));
    float h = 0.02f / Ksub;
    float y[3] = {(float) SEED_X, (float) SEED_Y, (float) SEED_Z};
    auto step = [&]() { for (int k = 0; k < Ksub; k++) { rk4<float>(y, h, b1, b2, 0.f);
        for (int j = 0; j < 3; j++) y[j] = std::min(std::max(y[j], POS_FLOOR), STATE_MAX); } };
    float prev[3] = {0, 0, 0}; bool rising[3] = {false, false, false};
    counts[0] = counts[1] = counts[2] = 0;
    for (int s = 0; s < 160000; s++) step();                 // burn in past the transient
    for (int s = 0; s < 40000; s++) {
        for (int k = 0; k < Ksub; k++) {
            rk4<float>(y, h, b1, b2, 0.f);
            for (int j = 0; j < 3; j++) y[j] = std::min(std::max(y[j], POS_FLOOR), STATE_MAX);
            for (int i = 0; i < 3; i++) {
                float g = (float) ((y[i] - cen[i]) * FOX_GAIN[i]);
                if (g > prev[i] + 1e-7f) rising[i] = true;
                else if (g < prev[i] - 1e-7f) { if (rising[i] && prev[i] > POP_GATE) counts[i]++; rising[i] = false; }
                prev[i] = g;
            }
        }
    }
}

// Fraction of default grass cycles whose peak order is grass→bunny→fox (the
// documented food-chain ordering: grass blooms, bunnies follow, foxes follow).
static double gateOrder() {
    double b1 = b1FromWild(0.5), b2 = 2.0, xs, ys, zs; equilibrium(b1, b2, xs, ys, zs);
    double cen[3] = {xs, ys, zs}, y[3] = {xs * 1.2, ys, zs}, h = 0.01;
    for (double t = 0; t < 5000; t += h) rk4<double>(y, h, b1, b2, 0.0);
    std::vector<double> pk[3]; double prev[3] = {0, 0, 0}; bool ris[3] = {false, false, false};
    double tau = 0;
    for (double t = 0; t < 3000; t += h) {
        rk4<double>(y, h, b1, b2, 0.0); tau += h;
        for (int i = 0; i < 3; i++) {
            double g = (y[i] - cen[i]) * FOX_GAIN[i];
            if (g > prev[i] + 1e-7) ris[i] = true;
            else if (g < prev[i] - 1e-7) { if (ris[i] && prev[i] > POP_GATE) pk[i].push_back(tau); ris[i] = false; }
            prev[i] = g;
        }
    }
    int ok = 0, tot = 0;
    for (size_t i = 0; i + 1 < pk[0].size(); i++) {
        double t0 = pk[0][i], t1 = pk[0][i + 1], tb = -1, tf = -1;
        for (double b : pk[1]) if (b > t0 && b < t1) { tb = b; break; }
        for (double f : pk[2]) if (f > t0 && f < t1) { tf = f; break; }
        if (tb > 0 && tf > 0) { tot++; if (t0 < tb && tb < tf) ok++; }
    }
    return tot ? (double) ok / tot : 0.0;
}

// C3/C4 magnitude ratio (dB) of the grass output at the default — a two-frequency
// DFT at the fundamental (1/P) and the period-2 subharmonic (1/2P). Guards the
// manual's "~−20 dB subharmonic" claim.
static double subharmonicDb() {
    const double PI = 3.14159265358979323846;
    double b1 = b1FromWild(0.5), b2 = 2.0, xs, ys, zs; equilibrium(b1, b2, xs, ys, zs);
    double y[3] = {xs * 1.2, ys, zs}, h = 0.01;
    for (double t = 0; t < 5000; t += h) rk4<double>(y, h, b1, b2, 0.0);
    double P = 61.387, f1 = 1.0 / P, f2 = 1.0 / (2 * P);
    double re1 = 0, im1 = 0, re2 = 0, im2 = 0, tau = 0; long N = 0;
    while (tau < 20 * P) {
        rk4<double>(y, h, b1, b2, 0.0); double g = y[0] - xs;
        re1 += g * std::cos(2 * PI * f1 * tau); im1 += g * std::sin(2 * PI * f1 * tau);
        re2 += g * std::cos(2 * PI * f2 * tau); im2 += g * std::sin(2 * PI * f2 * tau);
        tau += h; N++;
    }
    return 20.0 * std::log10(std::hypot(re2, im2) / std::hypot(re1, im1));
}

// Successive fox (z) local maxima after burn-in — tight cluster at the periodic
// default, broad return distribution in chaos.
static void foxMaxSpread(double b1, double b2, double h, double& lo, double& hi, int& count) {
    double xs, ys, zs; equilibrium(b1, b2, xs, ys, zs);
    double y[3] = {xs * 1.2, ys, zs};
    for (double t = 0; t < 4000; t += h) rk4<double>(y, h, b1, b2, 0.0);
    std::vector<double> mx; double zprev = y[2]; bool up = false;
    for (double t = 0; t < 12000; t += h) {
        rk4<double>(y, h, b1, b2, 0.0);
        if (y[2] < zprev && up) mx.push_back(zprev);
        up = y[2] > zprev; zprev = y[2];
    }
    lo = 1e9; hi = -1e9; int n = mx.size(), start = std::max(0, n - 80);
    for (int i = start; i < n; i++) { lo = std::min(lo, mx[i]); hi = std::max(hi, mx[i]); }
    count = n;
}

int main() {
    int fails = 0;

    // ── [1] analytic-equilibrium residual over a dense BALANCE×WILD grid ──
    double maxResid = 0;
    int invalid = 0;
    for (double bal = 0; bal <= 1.0001; bal += 0.05)
        for (double wild = 0; wild <= 1.0001; wild += 0.05) {
            double b1 = b1FromWild(wild), b2 = b2FromBal(bal), xs, ys, zs;
            if (!equilibrium(b1, b2, xs, ys, zs)) { invalid++; continue; }
            double s[3] = {xs, ys, zs}, d[3]; deriv<double>(s, d, b1, b2, 0.0);
            maxResid = std::max(maxResid, std::max({std::fabs(d[0]), std::fabs(d[1]), std::fabs(d[2])}));
        }
    printf("[1] equilibrium: grid residual max=%.2e  invalid centers=%d\n", maxResid, invalid);
    if (invalid) { printf("    FAIL: %d invalid equilibria inside the exposed box\n", invalid); fails++; }
    if (maxResid > 1e-9) { printf("    FAIL: equilibrium residual too large\n"); fails++; }

    // ── [2] full exposed-box stability (full KICK range, several pitches) ──
    int runs = 0, sf = 0; Diag dg; float worstFloorFrac = 0, worstCeilFrac = 0, worstUnforcedFloor = 0;
    for (double bal = 0; bal <= 1.0001; bal += 0.2)
        for (double wild = 0; wild <= 1.0001; wild += 0.1)
            for (float kick : {0.f, 2.f, 5.f, 10.f, 15.f, -5.f, -15.f})
                for (float dtau : {0.02f, 0.1f, 0.5f, 6.4f}) {     // incl. the integration-speed cap
                    runs++;
                    if (!run(bal, wild, kick, dtau, 8000, dg)) {
                        sf++; printf("    FAIL bal=%.1f wild=%.1f kick=%.1f dtau=%.2f\n", bal, wild, kick, dtau);
                    }
                    worstFloorFrac = std::max(worstFloorFrac, dg.floorFrac);
                    worstCeilFrac  = std::max(worstCeilFrac, dg.ceilFrac);
                    if (kick == 0.f) worstUnforcedFloor = std::max(worstUnforcedFloor, dg.floorFrac);
                }
    fails += sf;
    printf("[2] stability: %d/%d runs finite+positive+bounded  (worst floor%%=%.2f ceil%%=%.2f; unforced floor%%=%.2f)\n",
           runs - sf, runs, worstFloorFrac, worstCeilFrac, worstUnforcedFloor);
    if (worstCeilFrac > 0.f) { printf("    FAIL: state ceiling contacted — STATE_MAX is a corrupt-input backstop, not a limiter\n"); fails++; }
    if (worstUnforcedFloor > 0.5f) { printf("    FAIL: unforced operation leans on the positivity floor (%.2f%%)\n", worstUnforcedFloor); fails++; }

    // ── [2b] production-seed sweep: no advertised output may be DEAD ──
    // [2] seeds near each setting's own equilibrium, so it cannot see an alternate
    // basin. The module always starts from the fixed default-attractor seed, so sweep
    // from THAT across the whole box. A strongly negative mean is allowed (FOX is a
    // burst signal); the failure is a genuinely dead output — pinned to a rail (>98%)
    // with no usable swing (<0.5 V). This is the guard for the fox-extinction basin
    // that narrowing BALANCE to [1.75, 2.25] avoids. Rail occupancy is also reported.
    static const char* CH[3] = {"GRASS", "BUNNY", "FOX"};
    double worstMeanV = 0, worstRail = 0, worstP2pAtRail = 9; int deadHits = 0;
    double wmBal = 0, wmWild = 0; int wmCh = 0, wrCh = 0; double wrBal = 0, wrWild = 0;
    for (double bal = 0; bal <= 1.0001; bal += 0.1)
        for (double wild = 0; wild <= 1.0001; wild += 0.1) {
            double mv[3], rf[3], pp[3]; runFromSeed(bal, wild, 0.f, mv, rf, pp);
            for (int c = 0; c < 3; c++) {
                if (mv[c] < worstMeanV) { worstMeanV = mv[c]; wmBal = bal; wmWild = wild; wmCh = c; }
                if (rf[c] > worstRail)  { worstRail = rf[c]; wrBal = bal; wrWild = wild; wrCh = c; }
                // "dead" = essentially pinned to a rail with no usable swing
                if (rf[c] > 0.98 && pp[c] < 0.5) { deadHits++; worstP2pAtRail = std::min(worstP2pAtRail, pp[c]); }
            }
        }
    printf("[2b] production-seed sweep: most-negative mean = %+.2f V (%s @BAL %.1f WILD %.1f); "
           "peak rail occupancy = %.0f%% (%s @BAL %.1f WILD %.1f)\n",
           worstMeanV, CH[wmCh], wmBal, wmWild, 100.0 * worstRail, CH[wrCh], wrBal, wrWild);
    if (deadHits) { printf("    FAIL: %d setting(s) leave an output dead (railed >98%% with <0.5 V swing)\n", deadHits); fails++; }

    // ── [3] Hopf location: scan b1 for the Routh–Hurwitz crossing c2·c1−c0=0 ──
    double hopfB1 = NAN, prevRH = 0;
    for (double b1 = 1.0; b1 <= 3.0; b1 += 0.0005) {
        double b2 = 2.0, xs, ys, zs; if (!equilibrium(b1, b2, xs, ys, zs)) continue;
        double J[3][3]; jacobian(xs, ys, zs, b1, b2, J);
        double tr = J[0][0] + J[1][1] + J[2][2];
        double m01 = J[0][0]*J[1][1]-J[0][1]*J[1][0];
        double m02 = J[0][0]*J[2][2]-J[0][2]*J[2][0];
        double m12 = J[1][1]*J[2][2]-J[1][2]*J[2][1];
        double det = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
                   - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
                   + J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
        double c2 = -tr, c1 = m01 + m02 + m12, c0 = -det;
        double rh = c2 * c1 - c0;
        if (prevRH < 0 && rh >= 0) { /* leaving stability as b1 rises */ }
        if (prevRH > 0 && rh <= 0) hopfB1 = b1;
        prevRH = rh;
    }
    printf("[3] Hopf (b2=2): b1≈%.4f  (WILD≈%.4f)\n", hopfB1, std::sqrt((hopfB1 - 1.0) / 5.2));
    if (!(hopfB1 > 2.05 && hopfB1 < 2.20)) { printf("    FAIL: Hopf not near the expected b1≈2.1138\n"); fails++; }

    // ── [4] default calibration: measured period vs RATE_CAL (≤5 cents) ──
    run(0.5, 0.5, 0.f, 0.02f, 200000, dg);
    float cents = 1200.f * std::log2(RATE_CAL / dg.period);
    printf("[4] default cycle (bal=0.5 wild=0.5): amp=%.3f period=%.5f tau => RATE_CAL~%.5f (const %.3f, %+.2f cents)\n",
           dg.amp, dg.period, dg.period, RATE_CAL, cents);
    if (std::fabs(cents) > 5.f) { printf("    FAIL: RATE_CAL %.3f is %.1f cents off default period %.5f\n", RATE_CAL, cents, dg.period); fails++; }

    // ── [5] chaos verification: Lyapunov sign + fox-maxima return spread ──
    double lyDef = lyapunov(b1FromWild(0.5), 2.0, 0.01, 3000, 8000);
    double lyCha = lyapunov(3.0, 2.0, 0.01, 3000, 8000);
    double dLo, dHi, cLo, cHi; int dN, cN;
    foxMaxSpread(b1FromWild(0.5), 2.0, 0.005, dLo, dHi, dN);
    foxMaxSpread(3.0, 2.0, 0.005, cLo, cHi, cN);
    printf("[5] Lyapunov: default λ=%+.5f/tau (span %.3f)   chaos b1=3 λ=%+.5f/tau (span %.3f)\n",
           lyDef, dHi - dLo, lyCha, cHi - cLo);
    if (!(lyCha > 0.002)) { printf("    FAIL: canonical b1=3 not robustly chaotic (λ=%.5f)\n", lyCha); fails++; }
    if (!(lyDef < 0.001)) { printf("    FAIL: periodic default has positive λ (%.5f)\n", lyDef); fails++; }
    if (!((cHi - cLo) > 10 * (dHi - dLo))) { printf("    FAIL: fox-max return spread does not separate default from chaos\n"); fails++; }

    // ── [6] multistability scan (reported, not asserted) — two INDEPENDENT initial
    // conditions per b1 (not a continuation sweep). If a b1 settles to a different
    // grass amplitude from the two seeds, that column is multistable (coexisting
    // attractors); columns that agree are unique. This is the evidence that WILD is
    // not a monotonic chaos knob, not a CI assertion.
    { auto amp = [](double b1, double b2, double sx, double sy, double sz) {
          double xs, ys, zs; equilibrium(b1, b2, xs, ys, zs);
          double y[3] = {xs * sx, ys * sy, zs * sz};
          for (double t = 0; t < 6000; t += 0.01) rk4<double>(y, 0.01, b1, b2, 0.0);
          double mn = 1e9, mx = -1e9;
          for (double t = 0; t < 600; t += 0.01) { rk4<double>(y, 0.01, b1, b2, 0.0); mn = std::min(mn, y[0]); mx = std::max(mx, y[0]); }
          return mx - mn;
      };
      int multi = 0;
      printf("[6] multistability scan b1=2.40..2.60 (two seeds; *=seeds disagree → multistable): ");
      for (double b1 = 2.40; b1 <= 2.601; b1 += 0.02) {
          double aHi = amp(b1, 2.0, 1.2, 1.0, 1.0);   // seed above the equilibrium
          double aLo = amp(b1, 2.0, 0.6, 1.0, 1.6);   // a different basin probe
          bool disagree = std::fabs(aHi - aLo) > 0.05;
          if (disagree) multi++;
          printf("%s", disagree ? "*" : (aHi > 0.05 ? "o" : "."));
      }
      printf("  (%d multistable column%s)\n", multi, multi == 1 ? "" : "s"); }

    // ── [7] seed derivation + rest/output contract ──
    // A default-attractor seed for src/Foxes.cpp: integrate the default past its
    // transient and print a point on the cycle (deterministic, avoids silent startup).
    { double b1 = b1FromWild(0.5), b2 = 2.0, xs, ys, zs; equilibrium(b1, b2, xs, ys, zs);
      double y[3] = {xs * 1.2, ys, zs};
      for (double t = 0; t < 5000; t += 0.01) rk4<double>(y, 0.01, b1, b2, 0.0);
      printf("[7] default-attractor seed (paste into src/Foxes.cpp): {%.6ff, %.6ff, %.6ff}\n", y[0], y[1], y[2]); }
    { double mv[3], rf[3], pp[3]; runFromSeed(0.5, 0.0, 0.f, mv, rf, pp);   // WILD 0 → below Hopf → rest
      printf("    rest (WILD 0): mean outputs = (%+.2f, %+.2f, %+.2f) V\n", mv[0], mv[1], mv[2]);
      if (std::fabs(mv[0]) > 0.5 || std::fabs(mv[1]) > 0.5 || std::fabs(mv[2]) > 0.5) {
          printf("    FAIL: rest does not settle near 0 V\n"); fails++; } }

    // ── [8] input sanitization (mirrors src/Foxes.cpp guards) ──
    // Hostile V/OCT and KICK must be neutralized BEFORE the fmin/fmax clamps (which
    // map NaN to a bound → max speed). Confirm NaN/inf → rest/zero, huge finite →
    // bounded, and that the kernel stays finite under the sanitized hostile drive.
    { bool ok = true;
      auto sanVoct = [](float v){ if (!std::isfinite(v)) v = 0.f; return std::min(std::max(0.f + v, -8.f), 8.f); };
      auto sanKick = [](float v){ if (!std::isfinite(v)) v = 0.f; return std::min(std::max(v, -15.f), 15.f) * KICK_GAIN; };
      for (float bad : {(float) NAN, (float) INFINITY, -(float) INFINITY, 1e30f, -1e30f}) {
          float p = sanVoct(bad), k = sanKick(bad);
          if (!std::isfinite(p) || p < -8.f || p > 8.f) ok = false;
          if (!std::isfinite(k) || std::fabs(k) > 15.f * KICK_GAIN + 1e-3f) ok = false;
          double mv[3], rf[3], pp[3]; runFromSeed(0.5, 0.62, k, mv, rf, pp);
          for (int c = 0; c < 3; c++) if (!std::isfinite(mv[c])) ok = false;
      }
      // a NaN V/OCT must rest (0), NOT rail to max speed
      if (sanVoct((float) NAN) != 0.f) ok = false;
      printf("[8] sanitization: NaN/inf V·OCT→rest, KICK bounded, kernel finite under hostile drive: %s\n", ok ? "PASS" : "FAIL");
      if (!ok) fails++; }

    // ── [9] timestep convergence: production HSUB_MAX vs a fine reference ──
    // Periodic default → compare period (must land within 5 cents); chaos → compare
    // the fox-max spread, a statistical invariant (the trajectories must diverge, so
    // sample paths can't be compared). This is what freezes HSUB_MAX.
    { Diag dc; run(0.5, 0.5, 0.f, HSUB_MAX * 1.999f, 200000, dc);   // adaptive h → HSUB_MAX
      double pCoarse = dc.period;
      double b1 = b1FromWild(0.5), b2 = 2.0, xs, ys, zs; equilibrium(b1, b2, xs, ys, zs);
      double y[3] = {xs * 1.2, ys, zs};
      for (double t = 0; t < 4000; t += 0.005) rk4<double>(y, 0.005, b1, b2, 0.0);
      double prev = 0, lc = -1, sp = 0, T = 0; int np = 0;
      for (double t = 0; t < 8000; t += 0.005) {
          rk4<double>(y, 0.005, b1, b2, 0.0); T += 0.005; double cc = y[0] - xs;
          if (prev <= 0 && cc > 0) { double frac = cc / (cc - prev), tc = T - frac * 0.005;
              if (lc > 0) { sp += tc - lc; np++; } lc = tc; }
          prev = cc;
      }
      double pFine = np ? sp / np : 0, cents2 = 1200.0 * std::log2(pFine / pCoarse);
      double coLo, coHi, fiLo, fiHi; int ca, cb;
      foxMaxSpread(3.0, 2.0, 0.08, coLo, coHi, ca);    // coarse ~ HSUB_MAX
      foxMaxSpread(3.0, 2.0, 0.005, fiLo, fiHi, cb);   // fine reference
      double spCoarse = coHi - coLo, spFine = fiHi - fiLo;
      printf("[9] timestep convergence: default period coarse=%.4f fine=%.4f (%+.2f cents); "
             "chaos fox-max spread coarse=%.3f fine=%.3f\n", pCoarse, pFine, cents2, spCoarse, spFine);
      if (std::fabs(cents2) > 5.0) { printf("    FAIL: HSUB_MAX period differs from fine reference by >5 cents\n"); fails++; }
      if (spCoarse < 0.5 * spFine || spCoarse > 2.0 * spFine) { printf("    FAIL: chaos statistics not converged at HSUB_MAX\n"); fails++; } }

    // ── [10] peak-gate contract: no chatter at rest; all three fire when active ──
    { int gRest[3], gDef[3], gCha[3];
      gateCounts(0.5, 0.0, gRest);    // WILD 0 → rest
      gateCounts(0.5, 0.5, gDef);     // periodic default
      gateCounts(0.5, 0.62, gCha);    // canonical chaos
      printf("[10] peak gates: rest=(%d,%d,%d) default=(%d,%d,%d) chaos=(%d,%d,%d)  [grass,bunny,fox]\n",
             gRest[0], gRest[1], gRest[2], gDef[0], gDef[1], gDef[2], gCha[0], gCha[1], gCha[2]);
      if (gRest[0] + gRest[1] + gRest[2] > 0) { printf("    FAIL: gates chatter at rest\n"); fails++; }
      if (gDef[0] < 1 || gDef[1] < 1 || gDef[2] < 1) { printf("    FAIL: a species never fires at the default\n"); fails++; }
      if (gCha[0] < 1 || gCha[1] < 1 || gCha[2] < 1) { printf("    FAIL: a species never fires in chaos\n"); fails++; }
      double ord = gateOrder();
      printf("     default gate order grass→bunny→fox: %.0f%% of cycles\n", 100.0 * ord);
      if (ord < 0.8) { printf("    FAIL: default does not fire in grass→bunny→fox order\n"); fails++; } }

    // ── [11] regime classification guards the manual's claims ──
    // Default (WILD 0.5): period-2 (fox maxima alternate between two clusters) with a
    // small alternation, and a ~−20 dB grass subharmonic. Transition (WILD 0.57):
    // period-doubled well BEYOND the default (large alternation) but NOT chaos (λ ≤ 0).
    { auto maxSpread = [](double b1, double b2) {
          double xs, ys, zs; equilibrium(b1, b2, xs, ys, zs); double y[3] = {xs * 1.2, ys, zs};
          for (double t = 0; t < 5000; t += 0.005) rk4<double>(y, 0.005, b1, b2, 0.0);
          std::vector<double> m; double zp = y[2]; bool up = false;
          for (double t = 0; t < 3000; t += 0.005) {
              rk4<double>(y, 0.005, b1, b2, 0.0);
              if (y[2] < zp && up) m.push_back(zp);
              up = y[2] > zp; zp = y[2];
          }
          int n = (int) m.size(); double lo = 1e9, hi = -1e9;
          for (int i = std::max(0, n - 20); i < n; i++) { lo = std::min(lo, m[i]); hi = std::max(hi, m[i]); }
          // distinct clusters among the last maxima
          std::vector<double> c;
          for (int i = std::max(0, n - 20); i < n; i++) { bool nw = true;
              for (double u : c) if (std::fabs(m[i] - u) < 1e-3) nw = false;
              if (nw) c.push_back(m[i]); }
          return std::make_pair(hi - lo, (int) c.size());
      };
      auto dRes = maxSpread(b1FromWild(0.5), 2.0);
      auto tRes = maxSpread(b1FromWild(0.57), 2.0);
      double subDb = subharmonicDb();
      double lyTrans = lyapunov(b1FromWild(0.57), 2.0, 0.01, 3000, 6000);
      printf("[11] default: period-2=%s spread=%.4f subharmonic=%.1f dB;  transition: spread=%.3f λ=%+.5f/tau\n",
             dRes.second == 2 ? "yes" : "NO", dRes.first, subDb, tRes.first, lyTrans);
      if (dRes.second != 2) { printf("    FAIL: default is not the documented period-2 chase\n"); fails++; }
      if (subDb < -28.0 || subDb > -14.0) { printf("    FAIL: grass subharmonic %.1f dB is not the documented ~−20 dB\n", subDb); fails++; }
      if (tRes.first < 0.3) { printf("    FAIL: WILD 0.57 is not period-doubled beyond the default\n"); fails++; }
      if (lyTrans > 0.002) { printf("    FAIL: WILD 0.57 is chaotic, not the documented period-doubled transition\n"); fails++; } }

    // ── [12] display settling allowance: production seed -> WILD 0 rest ──
    // The UI discards one TRAIL_N-point window, then displays the next. Confirm that
    // this removes the large default-attractor decay before the resting portrait is
    // accepted. This guards the two-ring hand-off policy in src/Foxes.cpp.
    { double y[3] = {SEED_X, SEED_Y, SEED_Z};
      double ranges[2][3] = {};
      for (int window = 0; window < 2; window++) {
          double lo[3] = {1e9, 1e9, 1e9}, hi[3] = {-1e9, -1e9, -1e9};
          for (int p = 0; p < TRAIL_N; p++) {
              // CAP_TAU / 9 is exact enough here and mirrors the portrait's uniform
              // dimensionless-time spacing without involving audio sample rate.
              for (int k = 0; k < 9; k++) rk4<double>(y, CAP_TAU / 9.0, b1FromWild(0.0), 2.0, 0.0);
              for (int c = 0; c < 3; c++) { lo[c] = std::min(lo[c], y[c]); hi[c] = std::max(hi[c], y[c]); }
          }
          for (int c = 0; c < 3; c++) ranges[window][c] = hi[c] - lo[c];
      }
      printf("[12] portrait settling: first ranges=(%.3f,%.3f,%.3f) second=(%.4f,%.4f,%.4f)\n",
             ranges[0][0], ranges[0][1], ranges[0][2], ranges[1][0], ranges[1][1], ranges[1][2]);
      if (!(ranges[1][0] < 0.01 && ranges[1][1] < 0.01 && ranges[1][2] < 0.10
            && ranges[1][2] < 0.02 * ranges[0][2])) {
          printf("    FAIL: discarded first ring does not remove the WILD-to-rest display transient\n"); fails++; } }

    if (fails) { printf("FAIL: %d check(s) failed\n", fails); return 1; }
    printf("PASS: Hastings–Powell finite/positive/bounded; equilibrium exact; Hopf, RATE_CAL, chaos,\n");
    printf("      the rail guard, sanitization, gate contract, regimes, and portrait settling verified\n");
    return 0;
}
