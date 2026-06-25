// Offline auditioning for Axon: replicates the src/Axon.cpp kernel + output
// stage (tanh soft-clip + 20 Hz DC blocker) and writes WAV files so voicings
// can be heard without launching Rack. Keep in sync with the kernel.
//
//   g++ -O2 -o /tmp/axon_wav render_wav.cpp && /tmp/axon_wav
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

static constexpr float B_FIXED = 0.8f, HSUB_MAX = 0.05f, STATE_MAX = 10.f;
static constexpr int MIN_SUB = 4, MAX_SUB = 64;
static constexpr float RATE_CAL = 37.899004f, FREQ_C4 = 261.6256f;
static constexpr float TRIG_AMP = 0.6f, TRIG_TAU_MS = 3.f, OUT_GAIN = 1.0f;

static inline float clampf(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }
static inline int clampi(int x, int lo, int hi) { return x < lo ? lo : (x > hi ? hi : x); }

static inline void fhn(float v, float w, float I, float eps, float a, float& dv, float& dw) {
    dv = v - v * v * v / 3.f - w + I;
    dw = eps * (v + a - B_FIXED * w);
}
static inline void rk4(float& v, float& w, float h, float I, float eps, float a) {
    float k1v, k1w, k2v, k2w, k3v, k3w, k4v, k4w;
    fhn(v, w, I, eps, a, k1v, k1w);
    fhn(v + .5f*h*k1v, w + .5f*h*k1w, I, eps, a, k2v, k2w);
    fhn(v + .5f*h*k2v, w + .5f*h*k2w, I, eps, a, k3v, k3w);
    fhn(v + h*k3v, w + h*k3w, I, eps, a, k4v, k4w);
    v += h/6.f*(k1v + 2*k2v + 2*k3v + k4v);
    w += h/6.f*(k1w + 2*k2w + 2*k3w + k4w);
}

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
    float v = -1.2f, w = -0.6f, trigPulse = 0.f;
    DCBlock dc; dc.setFreq(20.f, fs);
    long N = (long)(dur*fs);
    std::vector<float> out; out.reserve(N);
    long trigSamp = trigPeriodS > 0 ? (long)(trigPeriodS*fs) : 0;
    for (long i = 0; i < N; i++) {
        if (trigSamp && (i % trigSamp == 0)) trigPulse = TRIG_AMP;
        float Itot = I + trigPulse;
        trigPulse *= std::exp(-1.f/(TRIG_TAU_MS*1e-3f*fs));
        float oct = octStart + (octEnd-octStart)*(float)i/N;
        float pitchHz = FREQ_C4*std::exp2(oct);
        float dtau = RATE_CAL*pitchHz/fs;
        int K = clampi((int)std::ceil(dtau/HSUB_MAX), MIN_SUB, MAX_SUB);
        float h = dtau/K;
        for (int k = 0; k < K; k++) rk4(v, w, h, Itot, eps, a);
        if (!std::isfinite(v)||!std::isfinite(w)) { v=-1.2f; w=-0.6f; }
        v = clampf(v, -STATE_MAX, STATE_MAX); w = clampf(w, -STATE_MAX, STATE_MAX);
        out.push_back(5.f*std::tanh(dc.process(v)*OUT_GAIN));
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
