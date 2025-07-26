#include "RtAudio.h"
#include <iostream>
#include <cstring>
#include <cmath>
#include "dsp.hpp"
#include "fifo.hpp"

using namespace std;

// type declarations
class scale_rotate {
private:
    sample_t factor; // for constellation scaling and rotating
public:
    void init(const sample_t&);
    void correct(sample_t&);
};

class local_oscillator {
private:
    double frq;
    double phi;
public:
    local_oscillator();
    sample_t mix(fifo<float>&);
    sample_t mix(float);
    void noise_measure();
    bool carrier_sync();
    double get_frq();
};

// global objects
#define RX_BUF_DEPTH    1024

#define PIx16 50.2654824574367

// fir lpf: kaiser win, fs=48k, fpass=750, fstop=1250, ripple=0.8db, attenuation=26.8db
static const float LPF[128] = {0.00248450157232583,0.00272126239724457,0.00292312563396990,0.00308355083689094,0.00319636939093471,0.00325592979788780,0.00325723853893578,0.00319609697908163,0.00306923151947558,0.00287441280670464,0.00261056004092097,0.00227783294394612,0.00187770312186331,0.00141300668474287,0.000887974747456610,0.000308240443700925,-0.000319177925121039,-0.000985919148661196,-0.00168234610464424,-0.00239763851277530,-0.00311990897171199,-0.00383634306490421,-0.00453336024656892,-0.00519679579883814,-0.00581209734082222,-0.00636453879997134,-0.00683944253250957,-0.00722241215407848,-0.00749956863000989,-0.00765778357163072,-0.00768491486087441,-0.00757003063336015,-0.00730362255126238,-0.00687780370935798,-0.00628649070858955,-0.00552555732429028,-0.00459296908229590,-0.00348888593725860,-0.00221573514863849,-0.000778252375312150,0.000816511455923319,0.00255921809002757,0.00443830434232950,0.00644007837399840,0.00854885391891003,0.0107471076771617,0.0130156790837646,0.0153339887037873,0.0176802836358547,0.0200319066643715,0.0223655831068754,0.0246577113866806,0.0268846843391657,0.0290231890976429,0.0310505200177431,0.0329448916018009,0.0346857346594334,0.0362539626657963,0.0376322567462921,0.0388052947819233,0.0397599637508392,0.0404855459928513,0.0409738682210445,0.0412194170057774,0.0412194170057774,0.0409738682210445,0.0404855459928513,0.0397599637508392,0.0388052947819233,0.0376322567462921,0.0362539626657963,0.0346857346594334,0.0329448916018009,0.0310505200177431,0.0290231890976429,0.0268846843391657,0.0246577113866806,0.0223655831068754,0.0200319066643715,0.0176802836358547,0.0153339887037873,0.0130156790837646,0.0107471076771617,0.00854885391891003,0.00644007837399840,0.00443830434232950,0.00255921809002757,0.000816511455923319,-0.000778252375312150,-0.00221573514863849,-0.00348888593725860,-0.00459296908229590,-0.00552555732429028,-0.00628649070858955,-0.00687780370935798,-0.00730362255126238,-0.00757003063336015,-0.00768491486087441,-0.00765778357163072,-0.00749956863000989,-0.00722241215407848,-0.00683944253250957,-0.00636453879997134,-0.00581209734082222,-0.00519679579883814,-0.00453336024656892,-0.00383634306490421,-0.00311990897171199,-0.00239763851277530,-0.00168234610464424,-0.000985919148661196,-0.000319177925121039,0.000308240443700925,0.000887974747456610,0.00141300668474287,0.00187770312186331,0.00227783294394612,0.00261056004092097,0.00287441280670464,0.00306923151947558,0.00319609697908163,0.00325723853893578,0.00325592979788780,0.00319636939093471,0.00308355083689094,0.00292312563396990,0.00272126239724457,0.00248450157232583};
#define LPF_LEN (sizeof(LPF)/sizeof(LPF[0]))

static sample_t tmp;

static fifo<float> q; // inter-thread data queue

static local_oscillator lo;
static fir_filter lpf;
static scale_rotate sr;

// functions
inline float amp(const sample_t& x) {
    return sqrt(x.I*x.I+x.Q*x.Q);
}

void scale_rotate::init(const sample_t& x) {
    float sra2 = x.I*x.I+x.Q*x.Q;
    factor.I = x.I/sra2;
    factor.Q = x.Q/sra2;
}

void scale_rotate::correct(sample_t& x) {
    sample_t xi = x;
    sample_t& xo = x;
    xo.I = xi.I*factor.I + xi.Q*factor.Q;
    xo.Q = xi.Q*factor.I - xi.I*factor.Q;
}

local_oscillator::local_oscillator() {
    frq = CARRIER_FRQ;
    phi = 0;
}

bool local_oscillator::carrier_sync() {
    return true;
}

sample_t local_oscillator::mix(fifo<float>& q) {
    float x = q.read();
    return mix(x);
}

sample_t local_oscillator::mix(float x) {
    sample_t mixed;
    mixed = (sample_t){(float)cos(phi)*x, (float)-sin(phi)*x};
    phi = phi + PIx16*frq/(SAMPLE_RATE*8);
    if (phi > PIx16) {
        phi -= PIx16;
    }
    return mixed;
}

double local_oscillator::get_frq() {
    return frq;
}

int rx_callback( void* /* out_buf */, void* in_buf, unsigned /* buf_samples */,  double /* timestamp */, RtAudioStreamStatus status, void* /* shared_data */) {
    if (status) cerr << "Overflow!" << endl;

    sample_t* buf = (sample_t*) in_buf;

    for (int i=0; i<RX_BUF_DEPTH; i++) {
        q.write(buf[i].L);
    }

    return 0;
}

float measure_noise(void) {
    // background noise measuring
    for (int j=0; j<(int)LPF_LEN; j++) {
        lpf.filter(lo.mix(q)); // make sure lpf is filled
    }

    float noise_level = 0;
    for (int i=0; i<NOISE_BODY; i++) {
        tmp = lpf.filter(lo.mix(q));
        noise_level += amp(tmp);
    }

    return noise_level;
}

void rx_demodulate(char* const data, unsigned& len_limit /* i/o */) {
    if (!data || !len_limit) {
        return;
    }

    // listen for preamble
    while (true) {
        // preamble
        for (unsigned short i=0; ; i++) {
            float x = q.read();
            cout << x << ' ';
            sample_t bb = lpf.filter(lo.mix(q));
            cout << bb.I << ' ' << bb.Q << endl;
        }
    }

    // protective margin
    for (int i=0; i<(int)LPF_LEN; i++) {
        lpf.filter(lo.mix(q));
    }
}

void ui(void) {
    lpf.init(LPF, LPF_LEN);

    cerr << "Silence! Measuring noise..." << endl;
    cout << "Noise level " << measure_noise() << endl;

    cerr << "Listening for data..." << endl;
    char s[256];
    unsigned len;

    while (true) {
        len = sizeof(s)-1;
        rx_demodulate(s, len);
        s[len] = '\0';
        if (len==0) {
            cerr << "Goodbye" << endl;
            break;
        } else {
            // cout << s << endl;
        }
    }
}

int main()
{
    RtAudio dev;
    if ( dev.getDeviceCount() < 1 ) {
        cerr << "No audio devices found!" << endl;
        return (-1);
    }
    RtAudio::StreamParameters parameters;
    parameters.deviceId = dev.getDefaultInputDevice();
    parameters.nChannels = 2;
    parameters.firstChannel = 0;
    unsigned int sampleRate = SAMPLE_RATE;
    unsigned int bufferFrames = RX_BUF_DEPTH;

    try {
        dev.openStream( NULL, &parameters, RTAUDIO_FLOAT32, sampleRate, &bufferFrames, &rx_callback, (void *)&q);
        dev.startStream();
    }
    catch ( RtAudioError& e ) {
        e.printMessage();
        return (-1);
    }

    ui();

    try {
        // Stop the stream
        dev.stopStream();
        dev.closeStream();
        return 0;
    }
    catch (RtAudioError& e) {
        e.printMessage();
        if (dev.isStreamOpen()) {
            dev.closeStream();
        }
        return (-1);
    }
}
