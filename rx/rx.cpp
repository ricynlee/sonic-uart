#include "RtAudio.h"
#include <iostream>
#include <cmath>
#include "dsp.hpp"
#include "fifo.hpp"
#include "sliwin.hpp"

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
    double get_frq();
    void set_frq(double);
};

class sliwin_peak_finder:protected sliwin {
private:
    bool _found;
    float _peak;
public:
    sliwin_peak_finder(size_t);
    bool slide(float);
    float peak(void);
};

// global objects
#define RX_BUF_DEPTH            1024
#define TH_COEF                 0.45f
#define HARD_TH                 0.1f

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

void local_oscillator::set_frq(double offset) {
    frq += offset;
}

sliwin_peak_finder::sliwin_peak_finder(size_t n):sliwin(n) {
    _found = false;
}

bool sliwin_peak_finder::slide(float v) {
    _found = false;
    sliwin::slide(v);
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

float sliwin_peak_finder::peak(void) {
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
                    float peak = win0.peak();
                    float th = win1.sum()*TH_COEF; // prior sum is used
                    if (win1.filled() && peak>HARD_TH && peak>th) { // decisioning
                        break;
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


        const double beta = 0.001480778208934; // pll proportional coef
        const double alpha = 0.002193245422464; // pll integral coef
        double frq = 0.0;
        double phi = 0.0;
        sample_t vco = {1.0f, 0.0f};

        sliwin_stdd phi_win(32); // phase lock checker
        sliwin_stdd frq_win(32); // frq stable checker

        const int CC_UB = 96; // confidence coef, upper bound
        const int CC_LB = 32; // confidence coef, lower bound
        int cc = CC_UB; // confidence coef
        int sc = 0; // (frq) stable count
        float frq_offset;
        float confidence; // snr related

        for (size_t i=0; i<CARRIER_BODY-4096; i++) {
            float sig = q.read();
            tmp = lpf.filter(lo.mix(sig));
            cout << sig << ' ' << tmp.I << ' ' << tmp.Q << ' ';
            if (i % 4 == 0) { // 1/4 decimation
                tmp = nbf.filter(tmp);
                cout << tmp.I << ' ' << tmp.Q << endl;
                double delta = atan2((double)(tmp.Q*vco.I-tmp.I*vco.Q), (double)(tmp.I*vco.I+tmp.Q*vco.Q)); // atan2(imag(tmp/vco), real(tmp/vco)), |vco|===1
                if (i % 256==0) { // 1/256 decimation
                    phi_win.slide(delta);
                    if (phi_win.filled()) {
                        if (phi_win.stdd()<PI/9) {
                            cc = cc - 1;
                            if (cc<CC_LB) cc = CC_LB;
                        } else {
                            cc = cc + 1;
                            if (cc>CC_UB) cc = CC_UB;
                        }
                    }

                    frq_win.slide(delta);
                    if (frq_win.filled()) {
                        if (frq_win.stdd()<0.025) {
                            sc = sc + 1;
                            frq_offset = frq_win.mean();
                            confidence = 1-fast_exp((float)-sc/cc);
                        }
                    }
                }

                // pll
                frq += beta*delta;
                phi += alpha*delta + PIx16*frq/8.0/SAMPLE_RATE/4; // 4 stands for decimation rate
                if (phi>PIx16) phi-=PIx16;
                if (phi<-PIx16) phi+=PIx16;
                vco.I = cos(phi);
                vco.Q = sin(phi);
            } else {
                cout << "0.0 0.0" << endl;
            }
        }

        lo.set_frq(frq_offset);
        cerr << "Carrier @ " << lo.get_frq() << endl;
        cerr << "Confidence " << confidence << endl;
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
