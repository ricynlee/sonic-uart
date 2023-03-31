#include "RtAudio.h"
#include <iostream>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>
#include "dsp.hpp"
#include "fifo.hpp"

using namespace std;

typedef struct {
    float left;
    float right;
} sample_t;

#if defined(__AVX__)
typedef union {
    __m256  d;
    __m256i di;
    float   f[8];
    int     i[8]; // only 32-bit int supported
} packed_t;

__attribute__((aligned(32))) float b[(ORDER+1)/8][8];
__attribute__((aligned(32))) float zi[(ORDER+1)/8][8];
__attribute__((aligned(32))) float zq[(ORDER+1)/8][8];
#else
float zi[ORDER+1];
float zq[ORDER+1];
#endif

void init_filter(void) {
#if defined(__AVX__)
    packed_t vindex = {
        .i = {0,16,32,48,64,80,96,112}
    };

    for (int i=0; i<(ORDER+1)/8; i++) {
        ((packed_t*)(zi+i))->d = _mm256_setzero_ps();
        ((packed_t*)(zq+i))->d = _mm256_setzero_ps();
        ((packed_t*)(b+i))->d = _mm256_i32gather_ps(B+i, vindex.d, 4); // reshape coef
    }
#else
    memset(zi, 0, sizeof(zi));
    memset(zq, 0, sizeof(zq));
#endif
}

inline float filteri(float x) {
    static int i = 0;
#if defined(__AVX__)
    packed_t s = {
        .f = {0}
    };

    ((packed_t*)(zi+i))->d = _mm256_permute_ps(((packed_t*)(zi+i))->d, 0x39);
    zi[i][3] = zi[i][7];
    zi[i][7] = x;

    i = (i+1)%(sizeof(zi)/sizeof(zi[0]));
    for (int j=0; j<(sizeof(zi)/sizeof(zi[0])); j++) {
        s.d = _mm256_fmadd_ps(((packed_t*)(zi+(j+i)%(sizeof(zi)/sizeof(zi[0]))))->d, ((packed_t*)(b+j))->d, s.d);
    }

    packed_t s_tmp;
    s_tmp.d = _mm256_hadd_ps(s.d, s.d);
    s.d = _mm256_hadd_ps(s_tmp.d, s_tmp.d); // or just return s_tmp.f[2] + [3] + [4] + [5]
    return (s.f[3] + s.f[4]);
#else // non-AVX

#endif
}

inline float filterq(float x) {
    static int i = 0;
#if defined(__AVX__)
    packed_t s = {
        .f = {0}
    };

    ((packed_t*)(zq+i))->d = _mm256_permute_ps(((packed_t*)(zq+i))->d, 0x39);
    zq[i][3] = zq[i][7];
    zq[i][7] = x;

    i = (i+1)%(sizeof(zq)/sizeof(zq[0]));
    for (int j=0; j<(sizeof(zq)/sizeof(zq[0])); j++) {
        s.d = _mm256_fmadd_ps(((packed_t*)(zq+(j+i)%(sizeof(zq)/sizeof(zq[0]))))->d, ((packed_t*)(b+j))->d, s.d);
    }

    packed_t s_tmp;
    s_tmp.d = _mm256_hadd_ps(s.d, s.d);
    s.d = _mm256_hadd_ps(s_tmp.d, s_tmp.d); // or just return s_tmp.f[2] + [3] + [4] + [5]
    return (s.f[3] + s.f[4]);
#else // non-AVX

#endif
}

fifo<float> q; // inter-thread data queue

int rx_callback( void* /* out_buf */, void* in_buf, unsigned /* buf_samples */,  double /* timestamp */, RtAudioStreamStatus status, void* /* shared_data */) {
    if (status) cerr << "Overflow!" << endl;

    sample_t* buf = (sample_t*) in_buf;

    for (int i=0; i<BUF_DEPTH; i++) {
        q.write((buf[i].left + buf[i].right)/2);
    }

    return 0;
}

inline float amp(float xi, float xq) {
    return sqrt(xi*xi+xq*xq);
}

void rx_demodulate(char* const data, unsigned& len /* i/o */) {
    if (!data || !len) {
        return;
    }

    // listening for preamble
    float wini[CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM];
    float winq[CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM];
    float peaki[3];
    float peakq[3];
    float peaka[3];
    float CARRIER[8];

    while (true) {
        for (int i=0; i<CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1; i++) {
            wini[i] = wini[i+1];
            winq[i] = winq[i+1];
        }

        wini[CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1] = 0;
        winq[CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1] = 0;
        for (int i=0; i<ASYNC_ACCUM; i++) {
            float x = q.read();
            wini[CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1] += filteri(x * COS[i & 7]);
            winq[CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1] += filterq(x * -SIN[i & 7]);
        }
        // wini[CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1] /= ASYNC_ACCUM;
        // winq[CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1] /= ASYNC_ACCUM;

        // cout << wini[CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1] << ' ' << winq[CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM-1] << endl;
        // continue;

        peaki[0] = peaki[1]; peaki[1] = peaki[2]; peaki[2] = 0;
        peakq[0] = peakq[1]; peakq[1] = peakq[2]; peakq[2] = 0;
        peaka[0] = peaka[1]; peaka[1] = peaka[2]; peaka[2] = 0;
        for (int i=0; i<CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM; i++) {
            peaki[2] += wini[i] * (1-2*MSEQ[i/(SAMPLES_PER_CHIP/ASYNC_ACCUM)]);
            peakq[2] += winq[i] * (1-2*MSEQ[i/(SAMPLES_PER_CHIP/ASYNC_ACCUM)]);
        }
        // peaki[2] /= CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM;
        // peakq[2] /= CHIPS*SAMPLES_PER_CHIP/ASYNC_ACCUM;
        peaka[2] = amp(peaki[2], peakq[2]);

        // cout << peaka[2] << endl;
        // continue;

        if (peaka[1]>CAPTURE_THRESH && peaka[1]>=peaka[2] && peaka[1]>=peaka[0]) {
            double phi = atan(peakq[1]/peaki[1]);
            for (int i=0; i<8; i++) {
                CARRIER[i] = cos(6*i*PI/8 - phi)/2;
            }
            break;
        }
    }

    // skip bubble
    for (int i=0; i<SAMPLES_PER_CHIP-ASYNC_ACCUM; i++) {
        q.read();
        filteri(0);
    }

    // data processing
    for (unsigned i=0; i<len; i++) {
        float y = 0;
        for (int j=0; j<SYNC_ACCUM; j++) {
            float x = q.read();
            y += filteri(x * CARRIER[j & 7]);
        }
        // y /= SYNC_ACCUM;

        // cout << y << endl;
        // continue;

        for (int j=0; j<8; j++) {
            if (y > DECISION_THRESH) {
                data[i] &= ~(1<<j); // 1'b0
            } else if (y < -DECISION_THRESH) {
                data[i] |= (1<<j); // 1'b1
            } else {
                len = i;
                break;
            }
        }
    }
}

void ui(void) {
    init_filter();

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
    unsigned int bufferFrames = BUF_DEPTH;

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
