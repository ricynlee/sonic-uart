#include "RtAudio.h"
#include <iostream>
#include <thread>
#include "fifo.hpp"
#include "dsp.hpp"

using namespace std;

fifo<float> q; // inter-thread data queue

int tx_callback( void* out_buf, void* /* in_buf */, unsigned /* buf_samples */, double /* timestamp */, RtAudioStreamStatus status, void* /* shared_data */) {
    if (status) {
        cerr << "Underflow!" << endl;
    }

    sample_t* buffer = (sample_t*) out_buf;

    static bool wearing = false;
    static bool R = false;
    unsigned fifo_size = q.size();

    if (fifo_size<TX_BUF_DEPTH) {
        if ( /*prior*/ wearing==true ) {
            R = !R;
        }
        wearing = false;
        for (int i=0; i<TX_BUF_DEPTH; i++) {
            buffer[i].R = 0;
            buffer[i].L = 0;
        }
    } else if (R) {
        wearing = true;
        for (int i=0; i<TX_BUF_DEPTH; i++) {
            buffer[i].R = q.read();
            buffer[i].L = 0;
        }
    } else /* !R */ {
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

    // preamble0
    for (int j=0; j<CHIPS; j++) {
        constel.I = 0;
        constel.Q = 0;
        for (int i=0; i<CHIP_PREFIX; i++) {
            sample = filter(constel);
            q.write(COS[i&7] * sample.I - SIN[i&7] * sample.Q);
        }
        constel.I = 1-2*MSEQ[j];
        constel.Q = 0;
        for (int i=0; i<CHIP_BODY; i++) {
            sample= filter(constel);
            q.write(COS[i&7] * sample.I - SIN[i&7] * sample.Q);
        }
        constel.I = 0;
        constel.Q = 0;
        for (int i=0; i<CHIP_SUFFIX; i++) {
            sample= filter(constel);
            q.write(COS[i&7] * sample.I - SIN[i&7] * sample.Q);
        }
    }

    // bubble
    constel.I = 0;
    constel.Q = 0;
    for (int i=0; i<(CHIP_PREFIX+CHIP_BODY+CHIP_SUFFIX); i++) {
        sample= filter(constel);
        q.write(COS[i&7] * sample.I - SIN[i&7] * sample.Q);
    }

    // preamble1 (channel-measuring symbol)
    constel.I = 1;
    constel.Q = 0;
    for (int i=0; i<(SYMBOL_CYCLIC_PREFIX+SYMBOL_BODY+SYMBOL_CYCLIC_SUFFIX); i++) {
        sample= filter(constel);
        q.write(COS[i&7] * sample.I - SIN[i&7] * sample.Q);
    }

    // frame length in octets
    len &= ~((-1)<<LENGTH_BITS);
    for (int j=0; j<LENGTH_BITS; j++) {
        int bit = (len>>j) & 1; // lsb first
        constel.I = (1-2*bit);
        constel.Q = 0;
        for (int i=0; i<(SYMBOL_CYCLIC_PREFIX+SYMBOL_BODY+SYMBOL_CYCLIC_SUFFIX); i++) {
            sample= filter(constel);
            q.write(COS[i&7] * sample.I - SIN[i&7] * sample.Q);
        }
    }

    // frame body
    for (unsigned k=0; k<len; k++) {
        for (int j=0; j<8; j+=MODEM) {
            unsigned char sym = (data[k]>>j) & ((1<<MODEM)-1); // lsb first
            switch (MODEM) {
                default: // PSK2
                    constel.I = (1-2*sym);
                    constel.Q = 0;
                    break;
                case PSK4:
                    constel.I = 0.75 - 1.5*(sym & 1);
                    constel.Q = 0.75 - 1.5*(sym >> 1);
                    break;
                case QAM16:
                    constel.I = 0.75 - 0.5*(sym & 3);
                    constel.Q = 0.75 - 0.5*(sym >> 2);
            }
            for (int i=0; i<(SYMBOL_CYCLIC_PREFIX+SYMBOL_BODY+SYMBOL_CYCLIC_SUFFIX); i++) {
                sample= filter(constel);
                q.write(COS[i&7] * sample.I - SIN[i&7] * sample.Q);
            }
        }
    }

    // pick up remainders in the filter & protective margin
    constel.I = 0;
    constel.Q = 0;
    for (int i=0; i<TX_BUF_DEPTH; i++) {
        sample= filter(constel);
        q.write(COS[i&7] * sample.I - SIN[i&7] * sample.Q);
    }
}

void ui(void) {
    init_filter();

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

int main()
{
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
