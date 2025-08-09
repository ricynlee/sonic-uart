#include "RtAudio.h"
#include <iostream>
#include <thread>
#include <cmath>
#include "fifo.hpp"
#include "dsp.hpp"

using namespace std;

#define TX_BUF_DEPTH 512    // common divisor of samples per chip and samples per symbol
                            // >ORDER
                            // cannot be too small (e.g., <256) in case of overflow/underflow

static float COS[8]; // {0.5, -0.353553390593274, 0, 0.353553390593274, -0.5, 0.353553390593275, 0, -0.353553390593274}; // 18kHz
static float SIN[8]; // {0, 0.353553390593274, -0.5, 0.353553390593274, 0, -0.353553390593273, 0.5, -0.353553390593274}; // 18kHz

// fir lpf
static const float LPF[LPF_LEN] = LPF_COEF;

static fifo<float> q; // inter-thread data queue

static fir_filter lpf;

int tx_callback( void* out_buf, void* /* in_buf */, unsigned /* buf_samples */, double /* timestamp */, RtAudioStreamStatus status, void* /* shared_data */) {
    if (status) {
        cerr << "Underflow!" << endl;
    }

    sample_t* buffer = (sample_t*) out_buf;

    static bool wearing = false;
    static bool RLn = false; // R/L channel wearing balancing
    unsigned fifo_size = q.size();

    if (fifo_size<TX_BUF_DEPTH) {
        if ( /*prior*/ wearing==true ) {
            RLn = !RLn;
        }
        wearing = false;
        for (int i=0; i<TX_BUF_DEPTH; i++) {
            buffer[i].R = 0;
            buffer[i].L = 0;
        }
    } else if (RLn) {
        wearing = true;
        for (int i=0; i<TX_BUF_DEPTH; i++) {
            buffer[i].R = q.read();
            buffer[i].L = 0;
        }
    } else /* !RLn */ {
        wearing = true;
        for (int i=0; i<TX_BUF_DEPTH; i++) {
            buffer[i].R = 0;
            buffer[i].L = q.read();
        }
    }

    return 0;
}

void tx_modulate(const char* const data, unsigned len) {
    if (len && (!data)) { // an exit sequence to rx side if len==0
        return;
    }

    sample_t constel, sample;

    // preamble: chirp for signal existence check, timing and gain control reference
    for (int i=0; i<PREAM_BODY; i++) {
        constel.I = chirp(i);
        constel.Q = 0;
        sample = lpf.filter(constel);
        q.write(COS[i&7]*sample.I - SIN[i&7]*sample.Q);
    }

    // bubble: avoid preamble-carrier interference
    constel.I = 0;
    constel.Q = 0;
    for (int i=0; i<BUBBLE_BODY; i++) {
        sample = lpf.filter(constel);
        q.write(COS[i&7]*sample.I - SIN[i&7]*sample.Q);
    }

    // carrier: for freq sync
    for (int i=0; i<CARRIER_BODY; i++) {
        constel.I = 1;
        constel.Q = 0;
        sample = lpf.filter(constel);
        q.write(COS[i&7]*sample.I - SIN[i&7]*sample.Q);
    }

    // bubble: avoid carrier-payload interference
    constel.I = 0;
    constel.Q = 0;
    for (int i=0; i<BUBBLE_BODY; i++) {
        sample = lpf.filter(constel);
        q.write(COS[i&7]*sample.I - SIN[i&7]*sample.Q);
    }

    // frame length in octets
    unsigned short header = 12345;
    for (int j=0; j<16; j++) {
        int bit = (header>>j) & 1; // lsb first
        constel.I = (1-2*bit);
        constel.Q = 0;
        for (int i=0; i<SYMBOL_BODY; i++) {
            sample = lpf.filter(constel);
            q.write(COS[i&7]*sample.I - SIN[i&7]*sample.Q);
        }
    }

    // // frame body
    // for (unsigned k=0; k<len; k++) {
    //     for (int j=0; j<8; j+=MODEM) {
    //         unsigned char sym = (data[k]>>j) & ((1<<MODEM)-1); // lsb first
    //         switch (MODEM) {
    //             default /* BPSK */:
    //                 constel.I = (1-2*sym);
    //                 constel.Q = 0;
    //                 break;
    //             case QPSK:
    //                 constel.I = 0.75 - 1.5*(sym & 1);
    //                 constel.Q = 0.75 - 1.5*(sym >> 1);
    //                 break;
    //             case QAM16:
    //                 constel.I = 0.75 - 0.5*(sym & 3);
    //                 constel.Q = 0.75 - 0.5*(sym >> 2);
    //         }
    //         for (int i=0; i<SYMBOL_BODY; i++) {
    //             sample = lpf.filter(constel);
    //             q.write(COS[i&7]*sample.I - SIN[i&7]*sample.Q);
    //         }
    //     }
    // }

    // pick up remainders in the filter & protective margin
    constel.I = 0;
    constel.Q = 0;
    for (int i=0; i<TX_BUF_DEPTH; i++) {
        sample = lpf.filter(constel);
        q.write(COS[i&7]*sample.I - SIN[i&7]*sample.Q);
    }
}

void ui(void) {
    // init lpf
    lpf.init(LPF, LPF_LEN);

    // init local oscillator lut
    for (int i=0; i<(int)(sizeof(COS)/sizeof(COS[0])); i++) {
        COS[i] = cos(2*PI*CARRIER_FRQ*i/SAMPLE_RATE)/2;
        SIN[i] = sin(2*PI*CARRIER_FRQ*i/SAMPLE_RATE)/2;
    }

    string s;
    while (true) {
        cout << "> ";
        if (getline(cin, s).eof()) {
            this_thread::sleep_for(chrono::milliseconds(200)); // avoid jamming of typing
            tx_modulate(NULL, 0);
            break;
        } else {
            this_thread::sleep_for(chrono::milliseconds(200)); // avoid jamming of typing
            tx_modulate(s.c_str(), s.length());
        }
    }

    while (q.size()) {
        this_thread::sleep_for(chrono::milliseconds(200)); // wait until all data is sent
    }
}

int main() {
    RtAudio dev;
    if (dev.getDeviceCount() < 1) {
        cerr << "No device!" << endl;
        return (-1);
    }
    RtAudio::StreamParameters parameters;
    parameters.deviceId = dev.getDefaultOutputDevice();
    parameters.nChannels = 2;
    parameters.firstChannel = 0;
    unsigned int sampleRate = SAMPLE_RATE;
    unsigned int bufferFrames = TX_BUF_DEPTH;

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
