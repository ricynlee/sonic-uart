#include "RtAudio.h"
#include <iostream>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include "dsp.hpp"
#include "fifo.hpp"

using namespace std;

typedef struct {
    int16_t left;
    int16_t right;
} sample_t;

inline int rx_filter(int x) { // iir lpf for demodulation
    int y;
    GEN_FILT(DEM_FILT_ORDER, y, x, DEM_B, DEM_A, DEM_B_FRAC_WID, DEM_A_FRAC_WID);
    return y;
}

inline void rx_filter_iq(int xi, int xq, int& yi, int& yq) { // iir lpf for demodulation
    yi = rx_filter(xi);
    GEN_FILT(DEM_FILT_ORDER, yq, xq, DEM_B, DEM_A, DEM_B_FRAC_WID, DEM_A_FRAC_WID);
}

fifo<int16_t> itdi; // inter-thread data interface

int rx_callback( void* /* out_buf */, void* in_buf, unsigned buf_samples,  double /* timestamp */, RtAudioStreamStatus status, void* /* shared_data */) {
    if (status) cerr << "Overflow!" << endl;

    sample_t* buf = (sample_t*) in_buf;

    int avg;
    for (int i=0; i<buf_samples; i++) {
        avg = buf[i].left + buf[i].right;
        avg = ((avg & 1) & ((avg >> 1) & 1)) + (avg >> 1); // improved round(mean(left, right))
        itdi.write(avg);
    }

    return 0;
}

inline int amp(int i, int q) {
#if 0
    // fast approx. sqrt(i*i+q*q)
    i = abs(i);
    q = abs(q);

    int s = i<q ? i : q;

    return (i+q-(s>>1)-(s>>2)+(s>>4));
#else
    return sqrt(i*i+q*q);
#endif
}

void rx_demodulate(char* const data, unsigned& len /* i/o */) {
    if (!data || !len) {
        return;
    }

    // listening for data & frame sync
    int A, I, Q;
    {
        int win[2][CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM] = {0}; // despreading: m-seq matching windows for I/Q
        int peak[3][3] = {0}; // peak finder, 0-i, 1-q, 2-amp
        while (true) {
            for (int i=0; i<CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1; i++) {
                win[0][i] = win[0][i+1];
                win[1][i] = win[1][i+1];
            }

            int& i0a = win[0][CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1];
            int& q0a = win[1][CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1];
            i0a = 0;
            q0a = 0;
            for (int i=0; i<ASYNC_ACCUM; i++) {
                int x = itdi.read();
                int i0, q0;
                rx_filter_iq((x*DEM_I_OSC[i & 7])>>OSC_FRAC_WID, (x*DEM_Q_OSC[i & 7])>>OSC_FRAC_WID, i0, q0);
                i0a += i0;
                q0a += q0;
            }
            i0a /= ASYNC_ACCUM;
            q0a /= ASYNC_ACCUM;

            cout << i0a << ' ' << q0a << ' ';

            peak[0][0] = peak[0][1];
            peak[0][1] = peak[0][2];
            peak[0][2] = 0;
            peak[1][0] = peak[1][1];
            peak[1][1] = peak[1][2];
            peak[1][2] = 0;
            peak[2][0] = peak[2][1];
            peak[2][1] = peak[2][2];
            for (int i=0; i<CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM; i++) {
                peak[0][2] += win[0][i] ^ (-MSEQ[i/(SAMPLES_PER_CHIP/ASYNC_ACCUM)]);
                peak[1][2] += win[1][i] ^ (-MSEQ[i/(SAMPLES_PER_CHIP/ASYNC_ACCUM)]);
            }
            peak[0][2] /= (CHIPS+1)*SAMPLES_PER_CHIP/ASYNC_ACCUM;
            peak[1][2] /= (CHIPS+1)*SAMPLES_PER_CHIP/ASYNC_ACCUM;
            peak[2][2] = amp(peak[0][2], peak[1][2]);

            cout << peak[2][2] << ' ' << win[0][CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1] << endl;
            continue;

            if (peak[2][1]>CAPTURE_THRESH && peak[2][1]>=peak[2][2] && peak[2][1]>=peak[2][0]) {
                // syncked/aligned
                A = peak[2][2];
                I = peak[0][2];
                Q = peak[1][2];
                break;
            }
        }
    }

    // carrier sync & skip bubble
    int16_t dem_i_osc[8];
    double phi = atan((double)Q/I);
    for (int i=0; i<8; i++) {
        dem_i_osc[i] = (int)(cos(6*i*PI/8 - phi)*16384*2);
        dem_i_osc[i] = (dem_i_osc[i]+1)/2;
    }
    for (int i=0; i<SAMPLES_PER_CHIP-ASYNC_ACCUM; i++) {
        rx_filter(itdi.read());
    }

    // data extraction & checking for ending
    {
        for (unsigned i=0; i<len; i++) {
            int i0a = 0;
            for (int j=0; j<SYNC_ACCUM; j++) {
                int x = itdi.read();
                int i0 = rx_filter(x*dem_i_osc[j & 7]);
                i0a += i0;
            }
            i0a /= SYNC_ACCUM;

            for (int j=0; j<8; j++) {
                if (i0a > DECISION_THRESH) {
                    data[i] &= ~(1<<j); // 1'b0
                } else if (i0a < -DECISION_THRESH) {
                    data[i] |= (1<<j); // 1'b1
                } else {
                    len = i;
                    break;
                }
            }
        }
    }
}

void ui(void) {
    cerr << "Listening for data..." << endl;
    char s[1024];
    unsigned len;

    while (true) {
        len = sizeof(s);
        rx_demodulate(s, len);
        cout << s << endl;
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
    unsigned int bufferFrames = 2048;

    try {
        dev.openStream( NULL, &parameters, RTAUDIO_SINT16, sampleRate, &bufferFrames, &rx_callback, (void *)&itdi);
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
