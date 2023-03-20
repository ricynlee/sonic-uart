#include "RtAudio.h"
#include <cstdint>
#include <iostream>
#include <string>
#include "fifo.hpp"
#include "dsp.hpp"

using namespace std;

#define SAMPLES_TX_ATOMIC_UNIT 1024 // common divisor of samples per DSSS code and samples per byte (2 symbols minimum can be a byte), cannot be too small (e.g., 128) so as not to under run

typedef struct {
    int16_t left;
    int16_t right;
} sample_t;

inline int tx_filter(int x) { // iir bpf for modulation
    int y;
    GEN_FILT(MOD_FILT_ORDER, y, x, MOD_B, MOD_A, MOD_B_FRAC_WID, MOD_A_FRAC_WID);
    return y;
}

fifo<int16_t> itdi; // inter-thread data interface

int tx_callback( void* out_buf, void* /* in_buf */, unsigned /* buf_samples */, double /* timestamp */, RtAudioStreamStatus status, void* /* shared_data */) {
    if (status) cerr << "Underflow!" << endl;

    static int remaining_samples = SAMPLES_TX_ATOMIC_UNIT;

    sample_t* buffer = (sample_t*) out_buf;
    int16_t x;
    int y;
    for (int i=0; i<SAMPLES_TX_ATOMIC_UNIT; i++) {
        unsigned fifo_size = itdi.size();
        if ((fifo_size >= SAMPLES_TX_ATOMIC_UNIT-i) || (fifo_size >= remaining_samples)) {
            if (itdi.peek(x)) {
                if ((--remaining_samples)==0) {
                    remaining_samples = SAMPLES_TX_ATOMIC_UNIT;
                }
            } else {
                return 1; // should never reach this line
            }
        } else {
            x = 0;
        }

        y = tx_filter(x);

        buffer[i].left = y;
        buffer[i].right = y;
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
            itdi.write(MOD_OSC_2PSK[0][j & 7] ^ (-MSEQ[i]));
        }
    }
    // bubble for 1 chip
    for (int i=0; i<SAMPLES_PER_CHIP; i++) {
        itdi.write(0);
    }
    // frame body
    // int bit;
    // for (const char* p=data; *p; p++) {
    //     for (int i=0; i<8; i++) {
    //         bit = ((*p)>>i) & 1; // lsb first
    //         for (int j=0; j<SAMPLES_PER_BIT; j++) {
    //             itdi.write(MOD_OSC_2PSK[bit][j & 7]);
    //         }
    //     }
    // }
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
        dev.openStream(&parameters, NULL, RTAUDIO_SINT16, sampleRate, &bufferFrames, &tx_callback, (void *)&itdi);
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
