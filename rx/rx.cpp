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
#define TH                      0.11f

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
    static float th = TH;

    // stop bit
    unsigned timeout = 0;
    do {
        sample = q.read();
        if (sample > 2*th) {
            th = sample/2;
            timeout = 0;
        } else {
            if (timeout==0x7fffu) {
                th = (th-TH)*0.9f + TH;
            }
            timeout = (timeout+1) & 0x7fffu;
        }
    } while(sample <= th);

    // start bit
    do {
        sample = q.read();
        if (sample > 2*th) {
            th = sample/2;
        }
    } while(sample > -th);

    for (int i=0; i<3; i++) {
        q.read();
    }

    // data bits
    for (int i=7; i>=0; i--) {
        q.read();
        sample = (q.read() + q.read())/2;
        if (sample >= 0) {
            c |= (1U << i);
        }
        q.read();
    }

    // skip possible ramping up of stop bits
    q.read();

    return c;
}

void ui(void) {
    cerr << "Listening for data..." << endl;    
    while (true) {
        char c = (char)rx_octet();
        if (c=='\0' || c=='\n')
            cout << endl;
        else
            cout << c;
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
