// Numeric proof that the shared neuron::rk4 template (src/neuron/integrator.hpp)
// produces the same results as the original hand-written RK4 steppers it
// replaced, for both Axon (FitzHugh–Nagumo, N=2) and Soma (Hindmarsh–Rose,
// N=3).
//
// The per-element expressions are byte-for-byte the same, so under the plugin's
// own flags (-O3 -funsafe-math-optimizations -march=nehalem) and under strict
// IEEE math the two are *bit-identical*. Under some other fast-math flag combos
// the optimiser is free to reassociate the loop form differently from the scalar
// form (FMA/grouping), producing last-ULP differences (~1e-8 on state that
// swings ±2 — dynamically irrelevant and inaudible). So this checks a tight
// absolute tolerance rather than bit-exactness, which stays valid regardless of
// the flags CI happens to use while still catching any real algorithmic change
// (those differ by orders of magnitude more).
//
//   g++ -O2 -std=c++17 -o /tmp/eq tools/integrator_equiv.cpp && /tmp/eq
//   exit 0 = within tolerance across the sweep; 1 = a real divergence.
#include "../src/neuron/integrator.hpp"
#include <cstdio>
#include <cmath>
#include <cstdint>

// ── FHN derivative + ORIGINAL Axon rk4 (verbatim copy, pre-extraction) ──
static constexpr float B_FIXED = 0.8f;
static inline void fA(float v, float w, float Itot, float eps, float a, float& dv, float& dw) {
    dv = v - v * v * v / 3.f - w + Itot;
    dw = eps * (v + a - B_FIXED * w);
}
static inline void rk4A_old(float& v, float& w, float h, float Itot, float eps, float a) {
    float k1v, k1w, k2v, k2w, k3v, k3w, k4v, k4w;
    fA(v,                  w,                  Itot, eps, a, k1v, k1w);
    fA(v + 0.5f * h * k1v, w + 0.5f * h * k1w, Itot, eps, a, k2v, k2w);
    fA(v + 0.5f * h * k2v, w + 0.5f * h * k2w, Itot, eps, a, k3v, k3w);
    fA(v + h * k3v,        w + h * k3w,        Itot, eps, a, k4v, k4w);
    v += h / 6.f * (k1v + 2.f * k2v + 2.f * k3v + k4v);
    w += h / 6.f * (k1w + 2.f * k2w + 2.f * k3w + k4w);
}
static inline void rk4A_new(float& v, float& w, float h, float Itot, float eps, float a) {
    float s[2] = {v, w};
    neuron::rk4<2>(s, h, [&](const float* y, float* d) { fA(y[0], y[1], Itot, eps, a, d[0], d[1]); });
    v = s[0]; w = s[1];
}

// ── HR derivative + ORIGINAL Soma rk4 (verbatim copy, pre-extraction) ──
static constexpr float A = 1.f, B = 3.f, C = 1.f, D = 5.f, XR = -1.6f;
static inline void fS(float x, float y, float z, float Itot, float r, float s,
                      float& dx, float& dy, float& dz) {
    dx = y - A*x*x*x + B*x*x - z + Itot;
    dy = C - D*x*x - y;
    dz = r * (s * (x - XR) - z);
}
static inline void rk4S_old(float& x, float& y, float& z, float h, float Itot, float r, float s) {
    float ax,ay,az, bx,by,bz, cx,cy,cz, dx,dy,dz;
    fS(x,           y,           z,           Itot, r, s, ax,ay,az);
    fS(x+.5f*h*ax,  y+.5f*h*ay,  z+.5f*h*az,  Itot, r, s, bx,by,bz);
    fS(x+.5f*h*bx,  y+.5f*h*by,  z+.5f*h*bz,  Itot, r, s, cx,cy,cz);
    fS(x+h*cx,      y+h*cy,      z+h*cz,      Itot, r, s, dx,dy,dz);
    x += h/6.f*(ax + 2*bx + 2*cx + dx);
    y += h/6.f*(ay + 2*by + 2*cy + dy);
    z += h/6.f*(az + 2*bz + 2*cz + dz);
}
static inline void rk4S_new(float& x, float& y, float& z, float h, float Itot, float r, float s) {
    float st[3] = {x, y, z};
    neuron::rk4<3>(st, h, [&](const float* Y, float* d) { fS(Y[0], Y[1], Y[2], Itot, r, s, d[0], d[1], d[2]); });
    x = st[0]; y = st[1]; z = st[2];
}

int main() {
    // A real algorithmic change perturbs the RK4 result by O(state); pure
    // fast-math reassociation stays at the last ULP (~1e-8). 1e-5 sits two-plus
    // orders below a real change and three above the reassociation noise.
    const double TOL = 1e-5;

    // Deterministic LCG sweep across realistic state + parameter ranges.
    uint32_t st = 12345u;
    auto rnd = [&](float lo, float hi) {
        st = st * 1664525u + 1013904223u;
        return lo + (hi - lo) * ((st >> 8) / float(1 << 24));
    };
    double maxA = 0.0, maxS = 0.0;
    for (long i = 0; i < 2000000; ++i) {
        float h = rnd(1e-4f, 0.05f);
        // Axon
        {
            float v = rnd(-3.f, 3.f), w = rnd(-2.f, 2.f);
            float I = rnd(-0.2f, 1.6f), eps = rnd(0.01f, 0.30f), a = rnd(0.4f, 1.0f);
            float v1 = v, w1 = w, v2 = v, w2 = w;
            rk4A_old(v1, w1, h, I, eps, a);
            rk4A_new(v2, w2, h, I, eps, a);
            maxA = std::fmax(maxA, std::fmax(std::fabs(v1 - v2), std::fabs(w1 - w2)));
        }
        // Soma
        {
            float x = rnd(-2.f, 2.f), y = rnd(-12.f, 4.f), z = rnd(-2.f, 6.f);
            float I = rnd(0.4f, 4.0f), r = rnd(0.001f, 0.05f), s = rnd(1.0f, 5.0f);
            float x1=x,y1=y,z1=z, x2=x,y2=y,z2=z;
            rk4S_old(x1, y1, z1, h, I, r, s);
            rk4S_new(x2, y2, z2, h, I, r, s);
            maxS = std::fmax(maxS, std::fmax(std::fabs(x1 - x2),
                                   std::fmax(std::fabs(y1 - y2), std::fabs(z1 - z2))));
        }
    }
    printf("max abs diff over 2,000,000 samples — Axon=%.3e  Soma=%.3e  (tol %.0e)\n",
           maxA, maxS, TOL);
    if (maxA > TOL || maxS > TOL) {
        printf("FAIL: stepper extraction changed results beyond reassociation noise\n");
        return 1;
    }
    printf("PASS: shared neuron::rk4 matches the original steppers within tolerance\n");
    return 0;
}
