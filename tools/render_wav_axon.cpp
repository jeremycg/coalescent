// Offline auditioning for Axon: uses the exact shared production FHN core and
// substep/safety policy with a SIMPLIFIED output stage — plain tanh
// soft-clip + 20 Hz DC blocker, no oversampling/FIR decimation — so voicings can
// be heard roughly without launching Rack. Not a bit-exact render of the plugin's
// output.
//
//   g++ -std=c++17 -O2 tools/render_wav_axon.cpp -o /tmp/axon_wav && /tmp/axon_wav
//
// Writes (in cwd):
//   axon_freerun.wav   free-run tone, CURRENT=0.6, plays a little V/OCT melody
//   axon_blips.wav     sub-threshold CURRENT=0.1, periodic triggers (one spike each)
//   axon_eps.wav       sweeps EPS to show spike sharpening
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include "../src/dsp/neuron_models.hpp"

using Core = coalescent::neuron::AxonCore;
static constexpr float TRIG_AMP = 0.6f, TRIG_TAU_MS = 3.f, OUT_GAIN = 1.0f;

static inline float clampf(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }

// Minimal one-pole high-pass (DC blocker), matching dsp::TRCFilter highpass at 20 Hz.
struct DCBlock {
    float xPrev = 0.f, yPrev = 0.f, a = 0.f;
    void setFreq(float fc, float fs) { float rc = 1.f/(2.f*(float)M_PI*fc); a = rc/(rc + 1.f/fs); }
    float process(float x) { float y = a*(yPrev + x - xPrev); xPrev = x; yPrev = y; return y; }
};

static void writeWav(const std::string& path, const std::vector<float>& s, int fs) {
    FILE* f = fopen(path.c_str(), "wb");
    int n = (int)s.size();
    auto u32 = [&](uint32_t v){ fwrite(&v,4,1,f); };
    auto u16 = [&](uint16_t v){ fwrite(&v,2,1,f); };
    fwrite("RIFF",1,4,f); u32(36 + n*2); fwrite("WAVE",1,4,f);
    fwrite("fmt ",1,4,f); u32(16); u16(1); u16(1); u32(fs); u32(fs*2); u16(2); u16(16);
    fwrite("data",1,4,f); u32(n*2);
    for (float x : s) { int16_t q = (int16_t)clampf(x*32767.f/5.f, -32768.f, 32767.f); fwrite(&q,2,1,f); }
    fclose(f);
    printf("wrote %s (%.2f s)\n", path.c_str(), (double)n/fs);
}

// Render one voice. trigPeriodS<=0 means no triggers (free-run).
static std::vector<float> render(float fs, double dur, float I, float eps, float a,
                                 float octStart, float octEnd, float trigPeriodS) {
    float state[Core::STATE_COUNT] = {Core::REST_V, Core::REST_W};
    float trigPulse = 0.f;
    DCBlock dc; dc.setFreq(20.f, fs);
    long N = (long)(dur*fs);
    std::vector<float> out; out.reserve(N);
    long trigSamp = trigPeriodS > 0 ? (long)(trigPeriodS*fs) : 0;
    for (long i = 0; i < N; i++) {
        if (trigSamp && (i % trigSamp == 0)) trigPulse = TRIG_AMP;
        float Itot = I + trigPulse;
        trigPulse *= std::exp(-1.f/(TRIG_TAU_MS*1e-3f*fs));
        float oct = octStart + (octEnd-octStart)*(float)i/N;
        float pitchHz = coalescent::neuron::PITCH_REFERENCE_HZ * std::exp2(oct);
        const coalescent::neuron::ScalarSchedule schedule =
            coalescent::neuron::scalarSchedule<Core>(pitchHz, fs);
        Core::advanceObservation(
            state, schedule.h, schedule.substeps, Itot, eps, a);
        Core::repair(state);
        out.push_back(5.f*std::tanh(dc.process(state[0])*OUT_GAIN));
    }
    return out;
}

int main() {
    const float fs = 48000.f;
    // Free-run melody: step through a few octaves on default voicing.
    std::vector<float> melody;
    float octs[] = {0.f, 0.25f, 0.5f, -0.25f, 0.f};
    for (float o : octs) { auto seg = render(fs, 0.5, 0.6f, 0.08f, 0.7f, o, o, 0); melody.insert(melody.end(), seg.begin(), seg.end()); }
    writeWav("axon_freerun.wav", melody, (int)fs);

    // Excitable blips: sub-threshold current, a trigger every 0.4 s.
    writeWav("axon_blips.wav", render(fs, 4.0, 0.1f, 0.08f, 0.7f, 0.f, 0.f, 0.4f), (int)fs);

    // EPS sweep: from sharp (0.02) to smooth (0.30) at fixed pitch, free-running.
    {
        std::vector<float> sweep; int steps = 6; float epsv[] = {0.02f,0.06f,0.10f,0.16f,0.24f,0.30f};
        for (int k = 0; k < steps; k++) { auto seg = render(fs, 0.7, 0.6f, epsv[k], 0.7f, 0.f, 0.f, 0); sweep.insert(sweep.end(), seg.begin(), seg.end()); }
        writeWav("axon_eps.wav", sweep, (int)fs);
    }
    return 0;
}
