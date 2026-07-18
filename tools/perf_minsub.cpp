// Profiling: does lowering Axon's MIN_SUB substep floor from 4 to 2 (a) cut CPU
// and (b) change the sound / pitch? Uses Axon's exact shared production model,
// bounded schedule, accepted observation boundary, and safety repair at os=4;
// only the deliberately varied substep floor differs.
//
//   g++ -O3 -funsafe-math-optimizations -march=native -std=c++17 tools/perf_minsub.cpp -o /tmp/pm && /tmp/pm
//
// Reports, per voicing: emergent fundamental (Hz) at MIN_SUB=4 vs 2 (pitch
// drift in cents), max |Δ| between the two output streams, and wall-clock time.
#include <cstdio>
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>
#include "../src/dsp/neuron_models.hpp"

static constexpr float FS        = 44100.f;
using Core = coalescent::neuron::AxonCore;

// Render N samples of the raw membrane voltage v (pre-tanh) at pitchHz, os=4,
// given the substep floor. Returns the per-output-sample v[] trace.
static void render(std::vector<float>& out, int N, float pitchHz,
                   float I, float eps, float a, int os, int minsub) {
    float state[Core::STATE_COUNT] = {Core::REST_V, Core::REST_W};
    const coalescent::neuron::ScalarSchedule schedule =
        coalescent::neuron::scalarScheduleWithMinimum<Core>(pitchHz, FS, os, minsub);
    out.resize(N);
    for (int n = 0; n < N; n++) {
        for (int o = 0; o < os; o++) {
            Core::advanceObservation(
                state, schedule.h, schedule.substeps, I, eps, a);
            Core::repair(state);
        }
        out[n] = state[0];
    }
}

// Fundamental from mean rising-zero-crossing spacing (after a warmup).
static float freqHz(const std::vector<float>& x) {
    int last = -1; double sum = 0; int cnt = 0;
    for (int n = 1; n < (int)x.size(); n++)
        if (x[n-1] <= 0.f && x[n] > 0.f) {
            if (last >= 0 && n > (int)x.size()/8) { sum += (n - last); cnt++; }
            last = n;
        }
    return cnt ? (float)(FS / (sum / cnt)) : 0.f;
}

static double bench(float pitchHz, float I, float eps, float a, int os, int minsub, int N) {
    std::vector<float> o;
    auto t0 = std::chrono::high_resolution_clock::now();
    volatile float sink = 0;
    for (int r = 0; r < 20; r++) { render(o, N, pitchHz, I, eps, a, os, minsub); sink += o.back(); }
    auto t1 = std::chrono::high_resolution_clock::now();
    (void)sink;
    return std::chrono::duration<double, std::milli>(t1 - t0).count() / 20.0;
}

int main() {
    const int N = (int)FS; // 1 s
    struct V { const char* name; float pitchHz, I, eps, a; };
    V voicings[] = {
        {"default C4",   coalescent::neuron::PITCH_REFERENCE_HZ,     0.6f, 0.08f, 0.7f},
        {"low C2",       coalescent::neuron::PITCH_REFERENCE_HZ / 4, 0.6f, 0.08f, 0.7f},
        {"sharp spike",  coalescent::neuron::PITCH_REFERENCE_HZ,     0.9f, 0.02f, 0.5f},
        {"high C6",      coalescent::neuron::PITCH_REFERENCE_HZ * 4, 0.6f, 0.08f, 0.7f},
    };
    printf("os=4, comparing MIN_SUB 4 vs 2\n");
    printf("%-13s %10s %10s %8s %10s %9s %9s\n",
           "voicing","f@MIN4","f@MIN2","cents","max|dv|","ms@MIN4","ms@MIN2");
    for (auto& v : voicings) {
        std::vector<float> o4, o2;
        render(o4, N, v.pitchHz, v.I, v.eps, v.a, 4, 4);
        render(o2, N, v.pitchHz, v.I, v.eps, v.a, 4, 2);
        float f4 = freqHz(o4), f2 = freqHz(o2);
        float cents = (f4>0 && f2>0) ? 1200.f*std::log2(f2/f4) : 0.f;
        float mx = 0; for (int n=0;n<N;n++) mx = std::max(mx, std::fabs(o4[n]-o2[n]));
        double ms4 = bench(v.pitchHz, v.I, v.eps, v.a, 4, 4, N);
        double ms2 = bench(v.pitchHz, v.I, v.eps, v.a, 4, 2, N);
        printf("%-13s %10.2f %10.2f %8.2f %10.4f %9.2f %9.2f  (%.0f%% faster)\n",
               v.name, f4, f2, cents, mx, ms4, ms2, 100*(1-ms2/ms4));
    }
    return 0;
}
