// Offline auditioning for Soma: replicates the src/Soma.cpp HR kernel (same
// constants, RK4, substep schedule) with a SIMPLIFIED output stage — plain tanh
// soft-clip + DC block, no oversampling/FIR decimation — so HR voicings can be
// heard roughly without launching Rack. Not a bit-exact render of the plugin's
// output; for the exact chain, run it in Rack.
//
//   g++ -O2 -o /tmp/soma_wav render_wav_soma.cpp && /tmp/soma_wav
//
// Writes: soma_bursting.wav (default), soma_chaos.wav (I=3.25), soma_tonic.wav,
//         soma_blips.wav (sub-threshold, periodic triggers fire bursts).
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

static constexpr float A=1,B=3,C=1,D=5,XR=-1.6f;
static constexpr float HSUB_MAX=0.05f, STATE_MAX=25.f;
static constexpr int MIN_SUB=2, MAX_SUB=64;                  // match src/neuron/Soma.cpp
static constexpr float RATE_CAL=55.364003f, FREQ_C4=261.6256f;  // match src/neuron/Soma.cpp (was a stale 14.925501)
static constexpr float TRIG_AMP=1.0f, TRIG_TAU_MS=5.f, OUT_GAIN=1.5f;

static inline float clampf(float x,float lo,float hi){return x<lo?lo:(x>hi?hi:x);}
static inline int clampi(int x,int lo,int hi){return x<lo?lo:(x>hi?hi:x);}

static inline void f(float x,float y,float z,float I,float r,float s,float&dx,float&dy,float&dz){
    dx=y-A*x*x*x+B*x*x-z+I; dy=C-D*x*x-y; dz=r*(s*(x-XR)-z);
}
static inline void rk(float&x,float&y,float&z,float h,float I,float r,float s){
    float ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz;
    f(x,y,z,I,r,s,ax,ay,az); f(x+.5f*h*ax,y+.5f*h*ay,z+.5f*h*az,I,r,s,bx,by,bz);
    f(x+.5f*h*bx,y+.5f*h*by,z+.5f*h*bz,I,r,s,cx,cy,cz); f(x+h*cx,y+h*cy,z+h*cz,I,r,s,dx,dy,dz);
    x+=h/6*(ax+2*bx+2*cx+dx); y+=h/6*(ay+2*by+2*cy+dy); z+=h/6*(az+2*bz+2*cz+dz);
}
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
    float x=-1.6f,y=-11.8f,z=2.f,trig=0.f; DCBlock dc; dc.setFreq(20.f,fs);
    long N=(long)(dur*fs); std::vector<float> out; out.reserve(N);
    long ts=trigPeriodS>0?(long)(trigPeriodS*fs):0;
    float pitchHz=FREQ_C4*std::exp2(oct); float dtau=RATE_CAL*pitchHz/fs;
    dtau=std::min(dtau,HSUB_MAX*MAX_SUB);   // mirror the production subTau cap: keeps h<=HSUB_MAX at extreme pitch
    int K=clampi((int)std::ceil(dtau/HSUB_MAX),MIN_SUB,MAX_SUB); float h=dtau/K;
    for(long i=0;i<N;i++){
        if(ts&&(i%ts==0)) trig=TRIG_AMP;
        float Itot=I+trig; trig*=std::exp(-1.f/(TRIG_TAU_MS*1e-3f*fs));
        for(int k=0;k<K;k++) rk(x,y,z,h,Itot,r,s);
        if(!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(z)){x=-1.6f;y=-11.8f;z=2.f;}
        x=clampf(x,-STATE_MAX,STATE_MAX);y=clampf(y,-STATE_MAX,STATE_MAX);z=clampf(z,-STATE_MAX,STATE_MAX);
        out.push_back(5.f*std::tanh(dc.p(x)*OUT_GAIN));
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
