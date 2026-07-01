// Regression test + calibration for Soma's Hindmarsh-Rose kernel.
//
// Replicates the src/Soma.cpp integrator (factored f() + RK4 substepping in
// dimensionless time, the HR sibling of the FitzHugh-Nagumo kernel in Axon)
// so the safety invariant, pitch calibration, and high-pitch accuracy can be
// checked without launching Rack. Keep in sync with the kernel.
//
//   g++ -O2 -o /tmp/soma_test soma_stability_test.cpp && /tmp/soma_test
//   Exit 0 = all checks pass, 1 = a check failed.
//
// Hindmarsh-Rose:
//   dx/dt = y - a x³ + b x² - z + I      (fast: membrane potential, audio out)
//   dy/dt = c - d x² - y                 (fast recovery / spiking)
//   dz/dt = r ( s (x - x_R) - z )         (slow adaptation / bursting)
//
// Pitch tracks the *within-burst spike rate*: RATE_CAL = the within-burst spike
// period (in dimensionless τ) at the default voicing, so dtau = RATE_CAL·pitchHz/fs
// gives a spike buzz at pitchHz (C4 at 0 V). Regime changes (tonic↔bursting↔chaos)
// pull the period, so pitch is open-loop/approximate — deliberate, documented.
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>

static constexpr float A = 1.f, B = 3.f, C = 1.f, D = 5.f, XR = -1.6f;
static constexpr float HSUB_MAX = 0.05f;
static constexpr int   MIN_SUB = 2, MAX_SUB = 64;
static constexpr float STATE_MAX = 25.f;     // backstop; normal |y| peaks ~12 at rest
static constexpr float FREQ_C4 = 261.6256f;

// default voicing
static constexpr float DEF_I = 2.0f, DEF_R = 0.006f, DEF_S = 4.0f;

static inline int clampi(int x, int lo, int hi) { return x < lo ? lo : (x > hi ? hi : x); }
static inline float clampf(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }

static inline void f(float x, float y, float z, float I, float r, float s,
                     float& dx, float& dy, float& dz) {
    dx = y - A*x*x*x + B*x*x - z + I;
    dy = C - D*x*x - y;
    dz = r * (s * (x - XR) - z);
}
static inline void rk4(float& x, float& y, float& z, float h, float I, float r, float s) {
    float ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz;
    f(x,            y,            z,            I,r,s, ax,ay,az);
    f(x+.5f*h*ax,   y+.5f*h*ay,   z+.5f*h*az,   I,r,s, bx,by,bz);
    f(x+.5f*h*bx,   y+.5f*h*by,   z+.5f*h*bz,   I,r,s, cx,cy,cz);
    f(x+h*cx,       y+h*cy,       z+h*cz,       I,r,s, dx,dy,dz);
    x += h/6.f*(ax+2*bx+2*cx+dx);
    y += h/6.f*(ay+2*by+2*cy+dy);
    z += h/6.f*(az+2*bz+2*cz+dz);
}
// Advance dimensionless time by dtau with adaptive substepping (mirrors Soma.cpp).
static inline void advance(float& x, float& y, float& z, float dtau, float I, float r, float s) {
    int K = clampi((int)std::ceil(dtau / HSUB_MAX), MIN_SUB, MAX_SUB);
    float h = dtau / K;
    for (int k = 0; k < K; k++) rk4(x, y, z, h, I, r, s);
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) { x=-1.6f; y=-11.8f; z=2.f; }
    x = clampf(x, -STATE_MAX, STATE_MAX);
    y = clampf(y, -STATE_MAX, STATE_MAX);
    z = clampf(z, -STATE_MAX, STATE_MAX);
}

// Median within-burst spike period (τ): integrate at a fine fixed step, collect
// upward-crossing intervals, and take the median of those < 2.5× the min so the
// long inter-burst gaps don't pollute the spike-rate estimate.
static double measurePeriodTau(float I, float r, float s) {
    float x=-1.6f, y=-11.8f, z=2.f, step=5e-4f;
    for (long i = 0; i < (long)(800.0/step); i++) rk4(x,y,z,step,I,r,s);
    float xp = x; std::vector<double> isi; double last=-1, tau=0;
    for (long i = 0; i < (long)(3000.0/step); i++) {
        rk4(x,y,z,step,I,r,s); tau += step;
        if (xp < 1.0f && x >= 1.0f) { if (last >= 0) isi.push_back(tau-last); last = tau; }
        xp = x;
    }
    if (isi.size() < 3) return -1;
    std::sort(isi.begin(), isi.end());
    double mn = isi.front();
    std::vector<double> within;
    for (double v : isi) if (v < 2.5*mn) within.push_back(v);
    return within[within.size()/2];   // median within-burst ISI
}

int main() {
    int failures = 0;

    // ── 1. CALIBRATION: within-burst spike period at the default voicing ──
    double Tspk = measurePeriodTau(DEF_I, DEF_R, DEF_S);
    if (Tspk <= 0) { printf("FAIL: default voicing does not spike (T=%.3f)\n", Tspk); failures++; }
    else {
        printf("CALIBRATION: within-burst spike period at default = %.6f τ\n", Tspk);
        printf("             => set RATE_CAL = %.6ff   (C4 buzz at 0 V, default voicing)\n", Tspk);
    }
    float RATE_CAL = (Tspk > 0) ? (float)Tspk : 14.9f;

    // ── 2. STABILITY sweep over I × r × s × pitch ──
    const float fs = 48000.f;
    int n = 0; float maxAbs = 0.f;
    for (float I = 0.4f; I <= 4.001f; I += 0.4f)
        for (float r = 0.001f; r <= 0.05f; r *= 1.7f)
            for (float s = 1.f; s <= 5.001f; s += 1.f)
                for (float oct = -4.f; oct <= 4.001f; oct += 4.f) {
                    float x=-1.6f, y=-11.8f, z=2.f;
                    float pitchHz = FREQ_C4 * std::exp2(oct);
                    float dtau = RATE_CAL * pitchHz / fs;
                    long ns = (long)(fs * 0.4f);
                    bool bad = false;
                    for (long i = 0; i < ns; i++) {
                        advance(x,y,z,dtau,I,r,s);
                        float m = std::max(std::fabs(x), std::max(std::fabs(y), std::fabs(z)));
                        maxAbs = std::max(maxAbs, m);
                        if (!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(z) || m > STATE_MAX + 1e-3f) {
                            printf("FAIL: escaped I=%.2f r=%.4f s=%.1f oct=%.0f: x=%g y=%g z=%g\n",
                                   I,r,s,oct,x,y,z); failures++; bad = true; break;
                        }
                    }
                    if (bad) break;
                    n++;
                }
    printf("STABILITY: %d combos stayed finite & bounded; max|state| observed = %.2f (clamp %.0f)\n",
           n, maxAbs, STATE_MAX);

    // ── 3. PITCH/ACCURACY: within-burst spike rate vs V/OCT with adaptive substeps ──
    auto audioHz = [&](float oct, float I, float r, float s) -> double {
        float x=-1.6f, y=-11.8f, z=2.f;
        float pitchHz = FREQ_C4 * std::exp2(oct);
        float dtau = RATE_CAL * pitchHz / fs;
        for (long i = 0; i < (long)(fs*0.3f); i++) advance(x,y,z,dtau,I,r,s);
        float xp = x; std::vector<double> isi; double last=-1; long sIdx=0;
        long maxS = (long)(fs*2.0f);
        for (; sIdx < maxS; sIdx++) {
            advance(x,y,z,dtau,I,r,s);
            if (xp < 1.0f && x >= 1.0f) { if (last>=0) isi.push_back((sIdx-last)); last=sIdx; }
            xp = x;
        }
        if (isi.size() < 3) return -1;
        std::sort(isi.begin(), isi.end());
        double mn = isi.front(); std::vector<double> within;
        for (double v : isi) if (v < 2.5*mn) within.push_back(v);
        double per = within[within.size()/2];
        return fs / per;
    };
    for (float oct = -1.f; oct <= 2.001f; oct += 1.f) {
        double want = FREQ_C4 * std::exp2(oct);
        double got = audioHz(oct, DEF_I, DEF_R, DEF_S);
        double cents = (got > 0) ? 1200.0*std::log2(got/want) : 0;
        printf("PITCH: oct=%+.0f want=%7.2f Hz got=%7.2f Hz (%+6.1f cents)\n", oct, want, got, cents);
        if (got < 0 || std::fabs(cents) > 60.0) { printf("  FAIL: spike pitch off >60 cents\n"); failures++; }
    }

    if (failures) { printf("\n%d CHECK(S) FAILED\n", failures); return 1; }
    printf("\nAll checks passed.\n");
    return 0;
}
