#include "RtAudio.h"
#include <iostream>
#include <cstring>
#include <immintrin.h>
#include "fifo.hpp"
#include "dsp.hpp"

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
__attribute__((aligned(32))) float z[(ORDER+1)/8][8];
#else
float z[ORDER+1];
#endif

void init_filter(void) {
#if defined(__AVX__)
    packed_t vindex = {
        .i = {0,16,32,48,64,80,96,112}
    };

    for (int i=0; i<(ORDER+1)/8; i++) {
        ((packed_t*)(z+i))->d = _mm256_setzero_ps();
        ((packed_t*)(b+i))->d = _mm256_i32gather_ps(B+i, vindex.d, 4); // reshape coef
    }
#else
    memset(z, 0, sizeof(z));
#endif
}

inline float filter(float x) {
    static int i = 0;
#if defined(__AVX__)
    packed_t s = {
        .f = {0}
    };

    ((packed_t*)(z+i))->d = _mm256_permute_ps(((packed_t*)(z+i))->d, 0x39);
    z[i][3] = z[i][7];
    z[i][7] = x;

    i = (i+1)%(sizeof(z)/sizeof(z[0]));
    for (int j=0; j<(sizeof(z)/sizeof(z[0])); j++) {
        s.d = _mm256_fmadd_ps(((packed_t*)(z+(j+i)%(sizeof(z)/sizeof(z[0]))))->d, ((packed_t*)(b+j))->d, s.d);
    }

    packed_t s_tmp;
    s_tmp.d = _mm256_hadd_ps(s.d, s.d);
    s.d = _mm256_hadd_ps(s_tmp.d, s_tmp.d); // or just return s_tmp.f[2] + [3] + [4] + [5]
    return (s.f[3] + s.f[4]);
#else // non-AVX

#endif
}

fifo<float> q; // inter-thread data queue

int tx_callback( void* out_buf, void* /* in_buf */, unsigned /* buf_samples */, double /* timestamp */, RtAudioStreamStatus status, void* /* shared_data */) {
    if (status) {
        cerr << "Underflow!" << endl;
    }

    static int remaining_samples = BUF_DEPTH;

    sample_t* buffer = (sample_t*) out_buf;
    float x;
    for (int i=0; i<BUF_DEPTH; i++) {
        unsigned fifo_size = q.size();
        if ((fifo_size >= BUF_DEPTH-i) || (fifo_size >= remaining_samples)) {
            if (q.peek(x)) {
                if ((--remaining_samples)==0) {
                    remaining_samples = BUF_DEPTH;
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
            q.write(COS[j & 7] * filter(1-2*MSEQ[i]));
        }
    }
    // bubble for 1 chip
    for (int i=0; i<SAMPLES_PER_CHIP; i++) {
        q.write(COS[i & 7] * filter(0));
    }
    // frame body
    int bit;
    for (const char* p=data; *p; p++) {
        for (int i=0; i<8; i++) {
            bit = ((*p)>>i) & 1; // lsb first
            for (int j=0; j<SAMPLES_PER_BIT; j++) {
                q.write(COS[j & 7] * filter(1-2*bit));
            }
        }
    }
    // pick up remainders in the filter
    for (int i=0; i<BUF_DEPTH; i++) {
        q.write(COS[i & 7] * filter(0));
    }
}

void ui(void) {
    init_filter();

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
        cerr << "No device!" << endl;
        return false;
    }
    RtAudio::StreamParameters parameters;
    parameters.deviceId = dev.getDefaultOutputDevice();
    parameters.nChannels = 2;
    parameters.firstChannel = 0;
    unsigned int sampleRate = SAMPLE_RATE;
    unsigned int bufferFrames = BUF_DEPTH;

    try {
        dev.openStream(&parameters, NULL, RTAUDIO_FLOAT32, sampleRate, &bufferFrames, &tx_callback, (void *)&q);
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
