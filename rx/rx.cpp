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
    bool carrier_sync();
    sample_t mix(const float&);
    double get_frq();
};

// global objects
#define RX_BUF_DEPTH    1024
#define CHIRP_DECIM     4

#define PIx16 50.2654824574367

static const float LPF[128] = LPF_COEF;
#define LPF_LEN (sizeof(LPF)/sizeof(LPF[0]))

#define MATCH_THRESH (1000./NOISE_BODY) // chirp capture decision threshold, multiple of noise_level

static float chirpseq[CHIRP_BODY/CHIRP_DECIM];
static float noise_level;
static sample_t tmp;

static fifo<float> q; // inter-thread data queue

static local_oscillator lo;
static fir_filter lpf;
static fir_filter mf; // match filter
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
    const double init_beta = 0.03;
    const double init_fflim = 10.0;
    const double adjcoef_beta = 0.90;
    const double adjcoef_fflim = 0.95;
    const double alpha = 0.0042;
    const int frqr_sample_period = 256;

    double beta = init_beta;
    double fflim = init_fflim;
    bool lock = false;
    double frqr_sample[8] = {0};

    frq = CARRIER_FRQ;
    phi = 0;

    for (int i=0; i<CARRIER_BODY; i++) {
        tmp = lpf.filter(lo.mix(q.read()));
        double phie = atan2(tmp.Q, tmp.I);

        // parameter adjusted for convergence
        if ((i % frqr_sample_period)==0) {
            double min = frq, max = frq;
            frqr_sample[(i / frqr_sample_period) & (sizeof(frqr_sample)/sizeof(double))] = frq;

            // lock check
            for (int j=0; j<(int)(sizeof(frqr_sample)/sizeof(double)); j++) {
                if (frqr_sample[j] < min) min = frqr_sample[j];
                if (frqr_sample[j] > max) max = frqr_sample[j];
            }
            lock = (max-min < fflim);

            // adjusting
            if (lock) {
                beta *= adjcoef_beta;
                fflim *= adjcoef_fflim;
            }
            else {
                beta /= adjcoef_beta;
                fflim /= adjcoef_fflim;
                if (beta > init_beta) beta = init_beta;
                if (fflim > init_fflim) fflim = init_fflim;
            }
        }

        // loop output
        phi = phi + alpha*phie;
        frq = frq + beta*phie;
    }

    return lock;
}

sample_t local_oscillator::mix(const float& x) {
    sample_t mixed = (sample_t){(float)cos(phi)*x, (float)-sin(phi)*x};
    phi = phi + 2*PI*frq/SAMPLE_RATE;
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

void measure_noise(void) {
    // background noise measuring
    for (int j=0; j<(int)LPF_LEN; j++) {
        lpf.filter(lo.mix(q.read())); // make sure lpf is filled
    }

    noise_level = 0;
    for (int i=0; i<NOISE_BODY; i++) {
        tmp = lpf.filter(lo.mix(q.read()));
        noise_level += amp(tmp);
    }
}

void rx_demodulate(char* const data, unsigned& len_limit /* i/o */) {
    if (!data || !len_limit) {
        return;
    }

    double doppler_coef;

    // listen for preamble
    while (true) {
        int latency;
        float mf_amp[3];
        
        // preamble0
        memset(mf_amp, 0, sizeof(mf_amp));
        for (unsigned short i=0; ; i++) {
            tmp = lpf.filter(lo.mix(q.read()));
            if (i % CHIRP_DECIM == 0) {
                mf_amp[0] = mf_amp[1];
                mf_amp[1] = mf_amp[2];
                mf_amp[2] = amp(mf.filter(tmp));

                if (mf_amp[1] > MATCH_THRESH*noise_level) {
                    if (mf_amp[1]>mf_amp[2] && mf_amp[1]>=mf_amp[0] && mf_amp[0]>=mf_amp[2]) {
                        float match_est = mf_amp[1] + (mf_amp[0]-mf_amp[2])/2;
                        latency = (int)(((match_est-mf_amp[0])/(mf_amp[1]-mf_amp[2])-1)*CHIRP_DECIM);
                    }
                    else if (mf_amp[1]>mf_amp[0] && mf_amp[1]>=mf_amp[2] && mf_amp[0]<=mf_amp[2]) {
                        float match_est = mf_amp[1] - (mf_amp[0]-mf_amp[2])/2;
                        latency = (int)((match_est-mf_amp[1])/(mf_amp[1]-mf_amp[0])*CHIRP_DECIM);
                    }
                    break;
                }
            }
        }
        mf.clear();

        cerr << "Preamble0  " << MATCH_THRESH*noise_level << endl;

        // skip bubble
        for (int i=0; i<BUBBLE_BODY+latency; i++) {
            lpf.filter(lo.mix(q.read()));
        }

        // carrier sync
        if (! lo.carrier_sync()) {
            cerr << "Carrier sync lost!" << endl;
            continue;
        }
        cerr << "Carrier @ " << lo.get_frq() << endl;
        doppler_coef = CARRIER_FRQ / lo.get_frq();

        // symbol sync
        bool sync = false;
        sample_t peak;
        float amp_peak = 0;
        memset(mf_amp, 0, sizeof(mf_amp));
        for (unsigned short i=0; i<CHIRP_BODY+BUBBLE_BODY*2; i++) {
            tmp = lpf.filter(lo.mix(q.read()));
            if (i % CHIRP_DECIM == 0) {
                tmp = mf.filter(tmp);
                mf_amp[0] = mf_amp[1];
                mf_amp[1] = mf_amp[2];
                mf_amp[2] = amp(tmp);

                if (mf_amp[2] > amp_peak) {
                    amp_peak = mf_amp[2];
                    peak = tmp;
                }

                if (mf_amp[1] > MATCH_THRESH*noise_level) {
                    if (mf_amp[1]>mf_amp[2] && mf_amp[1]>=mf_amp[0] && mf_amp[0]>=mf_amp[2]) {
                        float match_est = mf_amp[1] + (mf_amp[0]-mf_amp[2])/2;
                        latency = (int)(((match_est-mf_amp[0])/(mf_amp[1]-mf_amp[2])-1)*CHIRP_DECIM);
                    }
                    else if (mf_amp[1]>mf_amp[0] && mf_amp[1]>=mf_amp[2] && mf_amp[0]<=mf_amp[2]) {
                        float match_est = mf_amp[1] - (mf_amp[0]-mf_amp[2])/2;
                        latency = (int)((match_est-mf_amp[1])/(mf_amp[1]-mf_amp[0])*CHIRP_DECIM);
                    }
                    sync = true;
                    break;
                }
            }
        }
        mf.clear();

        if (! sync) {
            continue;
        }

        peak.I = peak.I * SYMBOL_BODY * CHIRP_DECIM / CHIRP_BODY;
        peak.Q = peak.Q * SYMBOL_BODY * CHIRP_DECIM / CHIRP_BODY;
        sr.init(peak);

        // skip bubble
        for (int i=0; i<(BUBBLE_BODY+latency)*doppler_coef; i++) {
            lpf.filter(lo.mix(q.read()));
        }
    }

    // packet length
    unsigned len = 0;
    for (int j=0; j<LENGTH_BITS; j++) {
        sample_t constel = (sample_t){0, 0};
        for (int i=0; i<SYMBOL_BODY; i++) {
            tmp = lpf.filter(lo.mix(q.read()));
            constel.I += tmp.I;
            constel.Q += tmp.Q;
        }
        constel.I /= SYMBOL_BODY;
        constel.Q /= SYMBOL_BODY;
        sr.correct(constel);

        cout << constel.I << ' ' << constel.Q << endl;

        len |= ((constel.I<=0)<<j);
    }

    if (len<len_limit) {
        len_limit=len;
    }

    // data processing
    for (unsigned k=0; k<len; k++) {
        char octet = 0;
        for (int j=0; j<8; j+=MODEM) {
            sample_t constel = (sample_t){0, 0};
            for (int i=0; i<SYMBOL_BODY; i++) {
                tmp = lpf.filter(lo.mix(q.read()));
                constel.I += tmp.I;
                constel.Q += tmp.Q;
            }
            constel.I /= SYMBOL_BODY;
            constel.Q /= SYMBOL_BODY;
            sr.correct(constel);

            cout << constel.I << ' ' << constel.Q << endl;

            unsigned char sym; // lsb first
            switch (MODEM) {
                default /* BPSK */:
                    sym = constel.I<0;
                    break;
                case QPSK:
                    sym = (constel.I<0) | ((constel.Q<=0)<<1);
                    break;
                case QAM16:
                    sym = (unsigned char)((0.75-constel.I)*2 + 0.5) | ((unsigned char)((0.75-constel.Q)*2 + 0.5)<<2);
            }
            octet |= (sym<<j);
        }

        if (k<len_limit) {
            data[k] = octet;
        }
    }

    // protective margin
    for (int i=0; i<(int)LPF_LEN; i++) {
        lpf.filter(lo.mix(q.read()));
    }
}

void ui(void) {
    lpf.init(LPF, LPF_LEN);

    for (int j=0; j<CHIRP_BODY; j+=CHIRP_DECIM) {
        chirpseq[j / CHIRP_DECIM] = chirp(j);
    }
    mf.init(chirpseq, CHIRP_BODY/CHIRP_DECIM);

    cerr << "Silence! Measuring noise..." << endl;
    measure_noise();

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
