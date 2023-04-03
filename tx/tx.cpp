#include "RtAudio.h"
#include <iostream>
#include <cstring>
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
    static bool right = false;
    unsigned fifo_size = q.size();

    if (fifo_size<TX_BUF_DEPTH) {
        if ( /*prior*/ wearing==true ) {
            right = !right;
        }
        wearing = false;
        for (int i=0; i<TX_BUF_DEPTH; i++) {
            buffer[i].right = 0;
            buffer[i].left = 0;
        }
    } else if (right) {
        wearing = true;
        for (int i=0; i<TX_BUF_DEPTH; i++) {
            buffer[i].right = q.read();
            buffer[i].left = 0;
        }
    } else /* !right */ {
        wearing = true;
        for (int i=0; i<TX_BUF_DEPTH; i++) {
            buffer[i].right = 0;
            buffer[i].left = q.read();
        }
    }

    return 0;
}

void tx_modulate(const char* const data, unsigned len) {
    if (len && (!data)) { // an exit sequence to rx side if len==0
        return;
    }

    // preamble
    for (int i=0; i<CHIPS; i++) {
        for (int j=0; j<MULTIPATH_MITIG; j++) {
            q.write(COS[j & 7] * filteri(0));
        }
        for (int j=0; j<SAMPLES_PER_CHIP-2*MULTIPATH_MITIG; j++) {
            q.write(COS[j & 7] * filteri(1-2*MSEQ[i]));
        }
        for (int j=0; j<MULTIPATH_MITIG; j++) {
            q.write(COS[j & 7] * filteri(0));
        }
    }

    // bubble
    for (int i=0; i<SAMPLES_PER_CHIP; i++) {
        q.write(COS[i & 7] * filteri(0));
    }

    // digital modem scheme - always in 2psk form
    char scheme = PSK_2;
    for (int i=0; i<SCHEME_BITS; i++) {
        int bit = (scheme>>i) & 1; // lsb first
        for (int j=0; j<MULTIPATH_MITIG; j++) {
            q.write(COS[j & 7] * filteri(0));
        }
        for (int j=0; j<SAMPLES_PER_SYM-2*MULTIPATH_MITIG; j++) {
            q.write(COS[j & 7] * filteri(1-2*bit));
        }
        for (int j=0; j<MULTIPATH_MITIG; j++) {
            q.write(COS[j & 7] * filteri(0));
        }
    }

    // frame length in octets - always in 2psk form
    len &= ~((-1)<<LENGTH_BITS);
    for (int i=0; i<LENGTH_BITS; i++) {
        int bit = (len>>i) & 1; // lsb first
        for (int j=0; j<MULTIPATH_MITIG; j++) {
            q.write(COS[j & 7] * filteri(0));
        }
        for (int j=0; j<SAMPLES_PER_SYM-2*MULTIPATH_MITIG; j++) {
            q.write(COS[j & 7] * filteri(1-2*bit));
        }
        for (int j=0; j<MULTIPATH_MITIG; j++) {
            q.write(COS[j & 7] * filteri(0));
        }
    }

    // frame body
    for (int i=0; i<len; i++) {
        for (int j=0; j<8; j++) {
            int bit = (data[i]>>j) & 1; // lsb first
            for (int k=0; k<MULTIPATH_MITIG; k++) {
                q.write(COS[k & 7] * filteri(0));
            }
            for (int k=0; k<SAMPLES_PER_SYM-2*MULTIPATH_MITIG; k++) {
                q.write(COS[k & 7] * filteri(1-2*bit));
            }
            for (int k=0; k<MULTIPATH_MITIG; k++) {
                q.write(COS[k & 7] * filteri(0));
            }
        }
    }

    // pick up remainders in the filter
    for (int i=0; i<TX_BUF_DEPTH; i++) {
        q.write(COS[i & 7] * filteri(0));
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
        return false;
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
