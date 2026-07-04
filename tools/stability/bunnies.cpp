// Standalone stability + calibration replica of Bunnies' predator–prey kernel
// (mirrors src/Bunnies.cpp). Sweeps BALANCE × WILD in both modes at several
// pitches over long runs; asserts x,y stay finite, positive, and bounded. Also
// reports the LV default period (→ RATE_CAL) and the RM Hopf coverage. The LV
// servo constants here match the module (STAB_K=0.5, STAB_FLOOR=0.2, V0≤Vmin+4).
//
//   g++ -O2 -o /tmp/t tools/stability/bunnies.cpp && /tmp/t    (exit 0 = pass)
#include <cstdio>
#include <cmath>
#include <algorithm>

static constexpr float POS = 1e-4f, STATE_MAX = 1e3f, HSUB_MAX = 0.05f;
static constexpr float STAB_K = 0.5f, STAB_FLOOR = 0.2f, MAX_STAB_STEP = 0.25f, LV_V0_RANGE = 4.f;
static constexpr float RM_B = 0.5f, RM_S = 1.0f;
static constexpr int   MIN_SUB = 2, MAX_SUB = 64;

struct St { float x, y; };
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
// returns finite+positive+bounded ok; also prey amplitude & period (2nd half)
static bool run(int mode, float balance, float wild, float dtau, int NS, float& amp, float& period, float& minState) {
    float gamma = 0, V0 = 0, K = 0, c = 0;
    if (mode == 0) { gamma = 0.2f + balance * 4.8f; V0 = (gamma + 1.f) + wild * LV_V0_RANGE; }
    else { c = std::min(0.15f + balance * 0.45f, 0.95f / RM_B); K = 1.2f + wild * 10.8f; }
    float cx = 1.f, cy = 1.f;
    if (mode == 1) { cx = c / (1.f - RM_B * c); cy = cx * (1.f - cx / K) / c; if (!(cx > 0 && cy > 0)) return true; }
    float dt = dtau; if (mode == 0) dt /= std::sqrt(gamma);
    dt = std::min(dt, HSUB_MAX * MAX_SUB);
    int Ksub = std::min(MAX_SUB, std::max(MIN_SUB, (int) std::ceil(dt / HSUB_MAX))); float h = dt / Ksub;
    float y[2] = {cx * 1.3f + 0.05f, cy * 0.9f + 0.02f};
    minState = 1e9f; float pmn = 1e9f, pmx = -1e9f, prev = 0, lc = -1; double sp = 0; int np = 0;
    for (int s = 0; s < NS; s++) {
        for (int k = 0; k < Ksub; k++) { rk4(y, h, mode, gamma, K, c, 0.f);
            y[0] = std::min(std::max(y[0], POS), STATE_MAX); y[1] = std::min(std::max(y[1], POS), STATE_MAX); }
        if (mode == 0) {
            float X = y[0], Y = y[1]; float V = gamma * (X - std::log(X)) + (Y - std::log(Y));
            float dVx = gamma * (1.f - 1.f / std::max(X, STAB_FLOOR)), dVy = (1.f - 1.f / std::max(Y, STAB_FLOOR));
            float st = STAB_K * dt;
            y[0] = std::max(X + std::clamp(-st * (V - V0) * dVx, -MAX_STAB_STEP, MAX_STAB_STEP), POS);
            y[1] = std::max(Y + std::clamp(-st * (V - V0) * dVy, -MAX_STAB_STEP, MAX_STAB_STEP), POS);
        }
        for (int j = 0; j < 2; j++) if (!std::isfinite(y[j]) || y[j] <= 0.f || y[j] > STATE_MAX) return false;
        minState = std::min(minState, std::min(y[0], y[1]));
        if (s > NS / 2) { pmn = std::min(pmn, y[0]); pmx = std::max(pmx, y[0]);
            float cc = y[0] - cx, t = s * dt; if (prev <= 0 && cc > 0) { if (lc > 0) { sp += t - lc; np++; } lc = t; } prev = cc; }
    }
    amp = pmx - pmn; period = np ? (float)(sp / np) : 0.f;
    return true;
}

int main() {
    int fails = 0, runs = 0; float amp, per, mn;
    for (int mode = 0; mode < 2; mode++)
        for (float bal = 0.f; bal <= 1.f; bal += 0.2f)
            for (float wild = 0.f; wild <= 1.f; wild += 0.2f)
                for (float dtau : {0.01f, 0.05f, 0.3f}) {
                    runs++;
                    if (!run(mode, bal, wild, dtau, 12000, amp, per, mn)) {
                        fails++; printf("  FAIL mode=%s bal=%.1f wild=%.1f dtau=%.2f\n", mode ? "RM" : "LV", bal, wild, dtau);
                    }
                }
    printf("stability: %d/%d runs finite+positive+bounded\n", runs - fails, runs);

    run(0, 0.5f, 0.4f, 0.02f, 120000, amp, per, mn);   // LV default
    float gamma = 0.2f + 0.5f * 4.8f;
    printf("LV default (bal=0.5 wild=0.4): amp=%.2f period=%.2f tau => RATE_CAL~%.2f\n", amp, per, per * std::sqrt(gamma));
    printf("RM Hopf coverage (does WILD cross rest->cycle?):\n");
    for (float bal : {0.f, 0.5f, 1.f}) { printf("  bal=%.1f: ", bal);
        for (float wild = 0.f; wild <= 1.f; wild += 0.2f) { run(1, bal, wild, 0.02f, 60000, amp, per, mn);
            printf("w%.1f:%s ", wild, amp > 0.05f ? "cyc" : "rest"); } printf("\n"); }

    if (fails) { printf("FAIL: %d runs left the positive-bounded region\n", fails); return 1; }
    printf("PASS: predator-prey finite, positive, bounded in both modes\n");
    return 0;
}
