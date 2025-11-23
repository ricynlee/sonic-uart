#include "RtAudio.h"
#include <iostream>
#include <cmath>
#include "dsp.hpp"
#include "fifo.hpp"
#include "sliwin.hpp"

using namespace std;

// global objects
static fifo<float> q; // inter-thread data queue

#define RX_BUF_DEPTH            1024
#define TH                      0.5

int rx_callback( void* /* out_buf */, void* in_buf, unsigned /* buf_samples */,  double /* timestamp */, RtAudioStreamStatus status, void* /* shared_data */) {
    if (status) cerr << "Overflow!" << endl;

    sample_t* buf = (sample_t*) in_buf;

    for (int i=0; i<RX_BUF_DEPTH; i++) {
        q.write(buf[i].L);
    }

    return 0;
}

uint8_t rx_octet() {
    uint8_t c = 0;
    float sample;
    while (q.read()<=TH);
    while (q.read()>=-TH);
    for (int i=0; i<3; i++)
        q.read();
    for (int i=7; i>=0; i--) {
        q.read();
        sample = (q.read() + q.read()) / 2;
        // cout << sample << ' ' << 0 << endl;
        if (sample >= 0) {
            c |= (1U << i);
        }
        q.read();
    }
    q.read();
    sample = (q.read() + q.read()) / 2;
    if (sample >= 0) {
        // ERROR
    }
    // q.read();
    return c;
}

void ui(void) {
    while (true) {
        cerr << "Listening for data..." << endl;
        while (1) {
            cout << "=" << (int)rx_octet() << endl;
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
