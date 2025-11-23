#include "RtAudio.h"
#include <iostream>
#include <thread>
#include <cmath>
#include <cstdint>
#include <cstring>
#include "fifo.hpp"
#include "dsp.hpp"

using namespace std;

#define TX_BUF_DEPTH 512    // common divisor of samples per chip and samples per symbol
                            // >ORDER
                            // cannot be too small (e.g., <256) in case of overflow/underflow

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

    if (fifo_size<TX_BUF_DEPTH) { // idle, wait until buffer is filled (packet len must be n*TX_BUF_DEPTH)
        if ( /*prior*/ wearing==true ) {
            RLn = !RLn;
        }
        wearing = false;
        for (size_t i=0; i<TX_BUF_DEPTH; i++) {
            buffer[i].R = 0;
            buffer[i].L = 0;
        }
    } else if (RLn) {
        wearing = true;
        for (size_t i=0; i<TX_BUF_DEPTH; i++) {
            buffer[i].R = q.read();
            buffer[i].L = 0;
        }
    } else /* !RLn */ {
        wearing = true;
        for (size_t i=0; i<TX_BUF_DEPTH; i++) {
            buffer[i].R = 0;
            buffer[i].L = q.read();
        }
    }

    return 0;
}

void tx_octet(unsigned char c) {
    sample_t constel, sample;
    unsigned d = c;
    d = (d<<1) | 1U;
    for (int i=9; i>=0; i--) {
        constel.I = ((d & (1U<<i))>>i)*2;
        constel.I = 0.5*(constel.I - 1);
        for (int j=0; j<4; j++) {
            sample = lpf.filter(constel);
            q.write(sample.I);
            cout << sample.I << endl;
        }
    }
}

void tx_packet(unsigned len, unsigned char data[]) {
    size_t n = 0;
    sample_t constel, sample;

    n += 4;
    constel.I = 0.5;
    for (int j=0; j<4; j++) {
        sample = lpf.filter(constel);
        q.write(sample.I);
        cout << sample.I << endl;
    }

    for (unsigned i=0; i<len; i++) {
        n += 40;
        n %= TX_BUF_DEPTH;
        tx_octet(data[i]);
    }

    // pick up remainders in the filter & protective margin
    constel.I = 0;
    for (size_t i=0; i<TX_BUF_DEPTH *2-n; i++) {
        sample = lpf.filter(constel);
        q.write(sample.I);
    }
}

void ui(void) {
    // init lpf
    lpf.init(LPF, LPF_LEN);

    unsigned char txdata[4096];

    while (true) {
        cerr << "> ";
        cin.getline((char*)txdata, sizeof(txdata));
        if (cin.eof()) {
            break;
        }

        if (cin.fail()) { // too many chars in buffer
            cin.clear(); // leave it to next read
        }

        this_thread::sleep_for(chrono::milliseconds(200)); // avoid jamming of keyboard typing
        tx_packet(cin.gcount(), txdata);
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
