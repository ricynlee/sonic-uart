#include "RtAudio.h"
#include <iostream>
#include <cmath>
#include "dsp.hpp"
#include "fifo.hpp"

using namespace std;

static const int RX_BUF_DEPTH = 1024;

inline float amp(sample_t x) {
    return sqrt(x.I*x.I+x.Q*x.Q);
}

sample_t sr; // for constellation scaling and rotating

void init_scale_rotate(sample_t x) {
    float sra2 = x.I*x.I+x.Q*x.Q;
    sr.I = x.I/sra2;
    sr.Q = x.Q/sra2;
}

inline void scale_rotate(sample_t& x) {
    sample_t xi = x;
    sample_t& xo = x;
    xo.I = xi.I*sr.I + xi.Q*sr.Q;
    xo.Q = xi.Q*sr.I - xi.I*sr.Q;
}

fifo<float> q; // inter-thread data queue

int rx_callback( void* /* out_buf */, void* in_buf, unsigned /* buf_samples */,  double /* timestamp */, RtAudioStreamStatus status, void* /* shared_data */) {
    if (status) cerr << "Overflow!" << endl;

    sample_t* buf = (sample_t*) in_buf;

    for (int i=0; i<RX_BUF_DEPTH; i++) {
        q.write(buf[i].L);
    }

    return 0;
}

float chirpseq[CHIRP_BODY/DECIMATION];

void rx_demodulate(char* const data, unsigned& len_limit /* i/o */) {
    if (!data || !len_limit) {
        return;
    }

    double lof = CARRIER_FRQ;
    double phi = 0;
    sample_t tmp;
    int latency;

    auto mix = [&phi, &lof](float x) {
        sample_t mixed = (sample_t){(float)cos(phi)*x, (float)-sin(phi)*x};
        phi = phi + 2*PI*lof/SAMPLE_RATE;
        return mixed;
    };

    // background noise measuring
    for (int j=0; j<512; j++) {
        q.read(); // make sure lpf is filled
    }

    float noise_level = 0;
    for (int i=0; i<NOISE_BODY; i++) {
        tmp = filter(mix(q.read()));
        noise_level += amp(tmp);
    }

    cerr << "noise level" << ' ' << noise_level << endl;

    // listening for preamble
    while (true) {
        float level = 0;
        for (int i=0; i<NOISE_BODY*4; i++) {
            tmp = filter(mix(q.read()));
            level += amp(tmp);
        }
        level /= 4;

        cerr << "level" << ' ' << level << endl;
    }

#if 0
    {
        sample_t win[ACCUMUL_TIMES*CHIPS] = {{0}};
        float peaka[3] = {0};

        while (true) {
            // shifting window (match filter)
            for (int i=0; i<ACCUMUL_TIMES*CHIPS-1; i++) {
                win[i] = win[i+1];
            }
            win[ACCUMUL_TIMES*CHIPS-1] = (sample_t){0, 0};

            for (int j=0; j<ACCUMUL_SAMPLES/16; j++) {
                sample_t accum = {0, 0};
                for (int i=0; i<16; i++) {
                    tmp = filter(mix(q.read()));
                    accum.I += tmp.I;
                    accum.Q += tmp.Q;
                }

                costas(accum);

                win[ACCUMUL_TIMES*CHIPS-1].I += accum.I;
                win[ACCUMUL_TIMES*CHIPS-1].Q += accum.Q;
            }

            extern float LO;
            cout << "LO=" << LO << endl;
            continue;

            win[ACCUMUL_TIMES*CHIPS-1].I /= ACCUMUL_SAMPLES;
            win[ACCUMUL_TIMES*CHIPS-1].Q /= ACCUMUL_SAMPLES;

            // capture threshold (by amplitude)
            float capture_thresh = 0;
            for (int i=0; i<ACCUMUL_TIMES*4; i++) {
                capture_thresh += amp(win[(CHIPS-4)*ACCUMUL_TIMES+i]);
            }
            capture_thresh = (capture_thresh-noise_level)*(CHIPS*0.5/4); // 4 past chips scanned

            // cout << capture_thresh << ' ';

            // capture
            tmp = (sample_t){0, 0};
            for (int i=0; i<CHIPS; i++) {
                for (int j=CHIP_PREFIX/ACCUMUL_SAMPLES; j<(CHIP_PREFIX+CHIP_BODY)/ACCUMUL_SAMPLES; j++) {
                    tmp.I += win[i*ACCUMUL_TIMES+j].I * (1-2*MSEQ[i]);
                    tmp.Q += win[i*ACCUMUL_TIMES+j].Q * (1-2*MSEQ[i]);
                }
            }
            peaka[0] = peaka[1];
            peaka[1] = peaka[2];
            peaka[2] = amp(tmp);

            // cout << peaka[2] << endl;
            // continue;

            if (capture_thresh < 16) {
                continue;
            }
            if (peaka[1]>peaka[2] && peaka[1]>=peaka[0] && peaka[0]>=peaka[2]) {
                float peaka_est = peaka[1] + (peaka[0]-peaka[2])/2;
                latency = (int)(((peaka_est-peaka[0])/(peaka[1]-peaka[2])-1)*ACCUMUL_SAMPLES);
                if (capture_thresh<peaka_est) {
                    break;
                }
            } else if (peaka[1]>peaka[0] && peaka[1]>=peaka[2] && peaka[0]<=peaka[2]) {
                float peaka_est = peaka[1] - (peaka[0]-peaka[2])/2;
                latency = (int)((peaka_est-peaka[1])/(peaka[1]-peaka[0])*ACCUMUL_SAMPLES);
                if (capture_thresh<peaka_est) {
                    break;
                }
            }
        }
    }

    // skip bubble
    for (int i=0; i<ACCUMUL_SAMPLES*(ACCUMUL_TIMES-1)+latency; i++) {
        filter(mix(q.read()));
    }

    // preamble1 & amplitude/phase aligning
    for (int i=0; i<SYMBOL_CYCLIC_PREFIX; i++) {
        filter(mix(q.read()));
    }
    sample_t sr = (sample_t){0, 0};
    for (int i=0; i<SYMBOL_BODY; i++) {
        tmp = filter(mix(q.read()));
        sr.I += tmp.I;
        sr.Q += tmp.Q;
    }
    sr.I /= (SYMBOL_BODY/128);
    sr.Q /= (SYMBOL_BODY/128);
    init_scale_rotate(sr);
    for (int i=0; i<SYMBOL_CYCLIC_SUFFIX; i++) {
        filter(mix(q.read()));
    }

    // packet length
    unsigned len = 0;
    for (int j=0; j<LENGTH_BITS; j++) {
        for (int i=0; i<SYMBOL_CYCLIC_PREFIX; i++) {
            filter(mix(q.read()));
        }
        sample_t constel = (sample_t){0, 0};
        for (int i=0; i<SYMBOL_BODY; i++) {
            tmp = filter(mix(q.read()));
            constel.I += tmp.I;
            constel.Q += tmp.Q;
        }
        constel.I /= (SYMBOL_BODY/128);
        constel.Q /= (SYMBOL_BODY/128);
        scale_rotate(constel);
        for (int i=0; i<SYMBOL_CYCLIC_SUFFIX; i++) {
            filter(mix(q.read()));
        }

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
            for (int i=0; i<SYMBOL_CYCLIC_PREFIX; i++) {
                filter(mix(q.read()));
            }
            sample_t constel = (sample_t){0, 0};
            for (int i=0; i<SYMBOL_BODY; i++) {
                tmp = filter(mix(q.read()));
                constel.I += tmp.I;
                constel.Q += tmp.Q;
            }
            constel.I /= (SYMBOL_BODY/128);
            constel.Q /= (SYMBOL_BODY/128);
            scale_rotate(constel);
            for (int i=0; i<SYMBOL_CYCLIC_SUFFIX; i++) {
                filter(mix(q.read()));
            }

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
    for (int i=0; i<TX_BUF_DEPTH; i++) {
        filter(mix(q.read()));
    }
#endif
}

void ui(void) {
    init_filter();

    for (int j=0; j<sizeof(chirpseq)/sizeof(chirpseq[0]); j++) {
        chirpseq[j] = chirp(j, true);
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
