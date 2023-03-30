#include "RtAudio.h"
#include <iostream>
#include <cstring>
#include <string>
#include <immintrin.h>
#include "fifo.hpp"
#include "dsp.hpp"

using namespace std;

#define SAMPLES_TX_ATOMIC_UNIT 1024 // common divisor of samples per DSSS code and samples per byte (2 symbols minimum can be a byte), >ORDER, and cannot be too small (e.g., 128) so as not to under run

typedef struct {
    float left;
    float right;
} sample_t;

#if defined(__AVX__)
typedef union {
    __m256 d;
    float  f[8];
    int    i[8]; // only 32-bit int supported
} packed_t;

__attribute__((aligned(4))) float b[(ORDER+1)/8][8];
__attribute__((aligned(4))) float z[(ORDER+1)/8][8];
#endif

void init_filter() {
    memset(z, 0, sizeof(z));
#if defined(__AVX__)

    packed_t vindex = {
        .i = {0,16,32,48,64,80,96,112}
    };
    for (int i=0; i<(ORDER+1)/8; i++) {
        ((packed_t*)(b+i*8))->d = _mm256_i32gather_ps(B+i, vindex.d, 4);
    }
#endif
}

inline float filter(float x) {
    static int i = 0;
#if defined(__AVX__)
    packed_t s = {
        .f = {0}
    };

    i = (i+sizeof(z)/sizeof(z[0])-1)%(sizeof(z)/sizeof(z[0]));
    for (int j=0; j<7; j++) {
        z[i][j] = z[i][j+1];
    }
    z[i][7] = x;

    i = (i+1)%(sizeof(z)/sizeof(z[0]));
    for (int j=0; j<(sizeof(z)/sizeof(z[0])); j++) {
        s.d = _mm256_fmadd_ps(((packed_t*)(z+(j+i)%(sizeof(z)/sizeof(z[0]))))->d, ((packed_t*)(b+j))->d, s.d);
    }
    float y = 0;
    for (int j=0; j<8; j++) {
        y += s.f[j];
    }

    i = (i+1)%(sizeof(z)/sizeof(z[0]));
    return y;
#else // non-AVX

#endif
}

fifo<float> q; // inter-thread data queue

int tx_callback( void* out_buf, void* /* in_buf */, unsigned /* buf_samples */, double /* timestamp */, RtAudioStreamStatus status, void* /* shared_data */) {
    if (status) cerr << "Underflow!" << endl;

    static int remaining_samples = SAMPLES_TX_ATOMIC_UNIT;

    sample_t* buffer = (sample_t*) out_buf;
    float x;
    for (int i=0; i<SAMPLES_TX_ATOMIC_UNIT; i++) {
        unsigned fifo_size = q.size();
        if ((fifo_size >= SAMPLES_TX_ATOMIC_UNIT-i) || (fifo_size >= remaining_samples)) {
            if (q.peek(x)) {
                if ((--remaining_samples)==0) {
                    remaining_samples = SAMPLES_TX_ATOMIC_UNIT;
                }
            } else {
                return 1; // should never reach this line
            }
        } else {
            x = 0;
        }

        buffer[i].left = x;
        buffer[i].right = x;
    }

    return 0;
}

void tx_modulate(const char* const data, unsigned len) {
    if (!data || !len) {
        return;
    }
    // frame header
    for (int i=0; i<CHIPS; i++) {
        for (int j=0; j<SAMPLES_PER_CHIP; j++) {
            q.write(SIN[j & 7] * filter(1-2*MSEQ[i]));
        }
    }
    // bubble for 1 chip
    for (int i=0; i<SAMPLES_PER_CHIP; i++) {
        q.write(SIN[i & 7] * filter(0));
    }
    // frame body
    int bit;
    for (const char* p=data; *p; p++) {
        for (int i=0; i<8; i++) {
            bit = ((*p)>>i) & 1; // lsb first
            for (int j=0; j<SAMPLES_PER_BIT; j++) {
                q.write(SIN[j & 7] * filter(1-2*bit));
            }
        }
    }
    // pick up remainders in the filter
    for (int i=0; i<SAMPLES_TX_ATOMIC_UNIT; i++) {
        q.write(SIN[i & 7] * filter(0));
    }
}

void ui(void) {
    string s;
    while (true) {
        cout << "> ";
        if (getline(cin, s).eof()) {
            break;
        }
        tx_modulate(s.c_str(), s.length());
    }
}

int main()
{
    RtAudio dev;
    if (dev.getDeviceCount() < 1) {
        cerr << "No audio devices found!" << endl;
        return false;
    }
    RtAudio::StreamParameters parameters;
    parameters.deviceId = dev.getDefaultOutputDevice();
    parameters.nChannels = 2;
    parameters.firstChannel = 0;
    unsigned int sampleRate = SAMPLE_RATE;
    unsigned int bufferFrames = SAMPLES_TX_ATOMIC_UNIT;

    try {
        dev.openStream(&parameters, NULL, RTAUDIO_SINT16, sampleRate, &bufferFrames, &tx_callback, (void *)&q);
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
