// Offline auditioning for Soma: uses the exact shared production HR core and
// substep/safety policy with a SIMPLIFIED output stage — plain tanh
// soft-clip + DC block, no oversampling/FIR decimation — so HR voicings can be
// heard roughly without launching Rack. Not a bit-exact render of the plugin's
// output; for the exact chain, run it in Rack.
//
//   g++ -std=c++17 -O2 tools/render_wav_soma.cpp -o /tmp/soma_wav && /tmp/soma_wav
//
// Writes: soma_bursting.wav (default), soma_chaos.wav (I=3.25), soma_tonic.wav,
//         soma_blips.wav (sub-threshold, periodic triggers fire bursts).
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include "../src/dsp/neuron_models.hpp"

using Core = coalescent::neuron::SomaCore;
static constexpr float TRIG_AMP=1.0f, TRIG_TAU_MS=5.f, OUT_GAIN=1.5f;

static inline float clampf(float x,float lo,float hi){return x<lo?lo:(x>hi?hi:x);}
struct DCBlock { float xp=0,yp=0,a=0; void setFreq(float fc,float fs){float rc=1.f/(2.f*(float)M_PI*fc);a=rc/(rc+1.f/fs);} float p(float x){float y=a*(yp+x-xp);xp=x;yp=y;return y;} };

static void writeWav(const std::string&path,const std::vector<float>&s,int fs){
    FILE*f=fopen(path.c_str(),"wb"); int n=(int)s.size();
    auto u32=[&](uint32_t v){fwrite(&v,4,1,f);}; auto u16=[&](uint16_t v){fwrite(&v,2,1,f);};
    fwrite("RIFF",1,4,f);u32(36+n*2);fwrite("WAVE",1,4,f);fwrite("fmt ",1,4,f);u32(16);u16(1);u16(1);u32(fs);u32(fs*2);u16(2);u16(16);
    fwrite("data",1,4,f);u32(n*2);
    for(float x:s){int16_t q=(int16_t)clampf(x*32767.f/5.f,-32768.f,32767.f);fwrite(&q,2,1,f);}
    fclose(f); printf("wrote %s (%.2f s)\n",path.c_str(),(double)n/fs);
}
static std::vector<float> render(float fs,double dur,float I,float r,float s,float oct,float trigPeriodS){
    float state[Core::STATE_COUNT] = {Core::REST_X, Core::REST_Y, Core::REST_Z};
    float trig=0.f; DCBlock dc; dc.setFreq(20.f,fs);
    long N=(long)(dur*fs); std::vector<float> out; out.reserve(N);
    long ts=trigPeriodS>0?(long)(trigPeriodS*fs):0;
    float pitchHz=coalescent::neuron::PITCH_REFERENCE_HZ * std::exp2(oct);
    const coalescent::neuron::ScalarSchedule schedule =
        coalescent::neuron::scalarSchedule<Core>(pitchHz, fs);
    for(long i=0;i<N;i++){
        if(ts&&(i%ts==0)) trig=TRIG_AMP;
        float Itot=I+trig; trig*=std::exp(-1.f/(TRIG_TAU_MS*1e-3f*fs));
        Core::advanceObservation(
            state, schedule.h, schedule.substeps, Itot, r, s);
        Core::repair(state);
        out.push_back(5.f*std::tanh(dc.p(state[0])*OUT_GAIN));
    }
    return out;
}
int main(){
    const float fs=48000.f;
    writeWav("soma_bursting.wav", render(fs,4.0,2.0f,0.006f,4.0f,0.f,0), (int)fs);
    writeWav("soma_chaos.wav",    render(fs,4.0,3.25f,0.006f,4.0f,0.f,0), (int)fs);
    writeWav("soma_tonic.wav",    render(fs,3.0,3.0f,0.05f,1.0f,0.f,0), (int)fs);
    writeWav("soma_blips.wav",    render(fs,5.0,0.6f,0.006f,4.0f,0.f,0.6f), (int)fs);
    return 0;
}
