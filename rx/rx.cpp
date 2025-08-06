#include "RtAudio.h"
#include <iostream>
#include <cstring>
#include <cmath>
#include "dsp.hpp"
#include "fifo.hpp"
#include "sliwin.hpp"
#include <deque>

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
    bool carrier_sync();
    double get_frq();
};

class sliwin_peak_finder:protected sliwin {
private:
    bool _found;
    float _peak;
public:
    sliwin_peak_finder(size_t);
    bool slide(float);
    float get_peak(void);
};

// global objects
#define RX_BUF_DEPTH            1024
#define TH_COEF                 0.4f

#define PIx16 50.2654824574367

// fir lpf
static const float LPF[LPF_LEN] = LPF_COEF;

// chirp preamble match filter coef
static float* MF;

static sample_t tmp;

static fifo<float> q; // inter-thread data queue

static local_oscillator lo;
static fir_filter lpf;
static fir_filter mf;
static biquad_filter nbf;
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

sliwin_peak_finder::sliwin_peak_finder(size_t n):sliwin(n) {
    _found = false;
}

bool sliwin_peak_finder::slide(float v) {
    _found = false;
    slide(v);
    if (size()>=3) {
        v = operator[](1); // hope that [1] is the peak
        if ((v>operator[](0) && v>=operator[](2)) || (v>=operator[](0) && v>operator[](2))) { // neighborhood peak found
            for (size_t i=3; i<size(); i++) {
                if (operator[](i) > v) {
                    return false;
                }
            }
            _peak = v;
            _found = true; // regional peak found
        }
    }
    return _found;
}

float sliwin_peak_finder::get_peak(void) {
    if (_found)
        return _peak;
    else
        return 0.0f;
}

int rx_callback( void* /* out_buf */, void* in_buf, unsigned /* buf_samples */,  double /* timestamp */, RtAudioStreamStatus status, void* /* shared_data */) {
    if (status) cerr << "Overflow!" << endl;

    sample_t* buf = (sample_t*) in_buf;

    for (int i=0; i<RX_BUF_DEPTH; i++) {
        q.write(buf[i].L);
    }

    return 0;
}

void rx_demodulate(char* const data, unsigned& len_limit /* i/o */) {
    if (!data || !len_limit) {
        return;
    }

    // preamble
    {        
        sliwin_peak_finder win0(32); // regional peak finder
        sliwin_sum win1(32);

        for (unsigned short i=0; ; i++) {
            float sig = q.read();
            tmp = lpf.filter(lo.mix(sig));
            cout << sig << ' ' << tmp.I << ' ' << tmp.Q << ' ';
            if (i % 8 == 0) { // decimation
                tmp = mf.filter(tmp); // matching filtering
                cout << tmp.I << ' ' << tmp.Q << endl;

                if (win0.slide(amp(tmp))) { // regional peak found
                    float peak = win0.get_peak();
                    float th = win1.sum()*TH_COEF; // prior sum is used
                    if (win1.filled()) { // decisioning
                        if (peak>th) { // preamble got
                            break;
                        }
                    }
                    win1.slide(peak);
                }
            } else {
                cout << "0.0 0.0" << endl;
            }
        }
    }

    // carrier sync
    {
        for (unsigned short i=0; i<BUBBLE_BODY-8+1024; i++) { // skip the bubble and coarsely fill nbf (approx. 1024-point delay), match filter has 8-point delay
            float sig = q.read();
            tmp = lpf.filter(lo.mix(sig));
            cout << sig << ' ' << tmp.I << ' ' << tmp.Q << ' ';
            if (i % 4 == 0) { // decimation
                nbf.filter(tmp); // fill in the nbf biquad iir filter
            }
            cout << "0.0 0.0" << endl;
        }

        const size_t MAV0_SIZE = 32;
        const size_t MAV1_SIZE = 32;

        sample_t vco = {1.0f, 0.0f};
        deque<float> mav0; // phase lock checker
        deque<float> mav1; // freq stable checker

        for (size_t i=0; i<CARRIER_BODY-4096; i++) {
            float sig = q.read();
            tmp = lpf.filter(lo.mix(sig));
            cout << sig << ' ' << tmp.I << ' ' << tmp.Q << ' ';
            if (i % 4 == 0) { // decimation
                tmp = nbf.filter(tmp);
                cout << tmp.I << ' ' << tmp.Q << endl;
                double delta = atan2(tmp.Q*vco.I-tmp.I*vco.Q, tmp.I*vco.I+tmp.Q*vco.Q); // atan2(imag(tmp/vco), real(tmp/vco)), |vco|===1
                if (i % 256==0) {
                    mav0.push_front(delta);
                    if (mav0.size()>MAV0_SIZE) {
                        mav0.pop_back();
                        
                    }
                }
            } else {
                cout << "0.0 0.0" << endl;
            }
        }
    }

    len_limit = 0;
}

void init_filters() {
    // initialize filters
    lpf.init(LPF, LPF_LEN);

    const size_t MF_LEN = PREAM_BODY/8;
    MF = (float*)malloc(sizeof(float)*MF_LEN);
    for (size_t i=0; i<MF_LEN; i++) {
        MF[i] = chirp(PREAM_BODY-i*8);
    }
    mf.init(MF, MF_LEN);
    free(MF);

    const float IIR[6] = { // chebyshev type ii, 2nd order, fs 12k, fstop 40, astop 20db
        0.0993952155113220, -0.198703244328499, 0.0993952155113220, // B, numerator
        1, -1.98742473125458, 0.987511932849884, // A, denominator
    };
    nbf.init(IIR);
}

void ui(void) {
    init_filters();

    // protective margin, filling LPF
    for (int i=0; i<(int)LPF_LEN; i++) {
        lpf.filter(lo.mix(q));
    }

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
