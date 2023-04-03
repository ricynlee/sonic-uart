#include "RtAudio.h"
#include <iostream>
#include <cmath>
#include "dsp.hpp"
#include "fifo.hpp"

using namespace std;

inline float amp(float xi, float xq) {
    return sqrt(xi*xi+xq*xq);
}

float srx, sry; // for constellation scaling and rotating

void init_scale_rotate(float sra, float sri, float srq) {
    sra *= sra;
    srx = sri/sra*(SAMPLES_PER_CHIP*CHIPS/DOWN_SAMPLE)/SAMPLES_PER_BIT;
    sry = srq/sra*(SAMPLES_PER_CHIP*CHIPS/DOWN_SAMPLE)/SAMPLES_PER_BIT;
}

inline void scale_rotate(float& xi, float& xq) {
    float ini = xi;
    float inq = xq;
    xi = ini*srx + inq*sry;
    xq = inq*srx - ini*sry;
}

fifo<float> q; // inter-thread data queue

int rx_callback( void* /* out_buf */, void* in_buf, unsigned /* buf_samples */,  double /* timestamp */, RtAudioStreamStatus status, void* /* shared_data */) {
    if (status) cerr << "Overflow!" << endl;

    sample_t* buf = (sample_t*) in_buf;

    for (int i=0; i<1024; i++) {
        q.write(buf[i].left);
    }

    return 0;
}

void rx_demodulate(char* const data, unsigned& len_limit /* i/o */) {
    if (!data || !len_limit) {
        return;
    }

    // listening for preamble
    {
        float wini[CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE] = {0};
        float winq[CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE] = {0};
        float peaki[3];
        float peakq[3];
        float peaka[3] = {0};

        while (true) {
            for (int i=0; i<CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE-1; i++) {
                wini[i] = wini[i+1];
                winq[i] = winq[i+1];
            }
            wini[CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE-1] = 0;
            winq[CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE-1] = 0;

            for (int i=0; i<DOWN_SAMPLE; i++) {
                float x = q.read();
                wini[CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE-1] += filteri(x * COS[i & 7] / 2);
                winq[CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE-1] += filterq(x * -SIN[i & 7] / 2);
            }
            wini[CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE-1] /= DOWN_SAMPLE;
            winq[CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE-1] /= DOWN_SAMPLE;

            // cout << wini[CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE-1] << ' ' << winq[CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE-1] << ' ';

            peaki[0] = peaki[1]; peaki[1] = peaki[2]; peaki[2] = 0;
            peakq[0] = peakq[1]; peakq[1] = peakq[2]; peakq[2] = 0;
            peaka[0] = peaka[1]; peaka[1] = peaka[2];
            for (int i=0; i<CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE; i++) {
                peaki[2] += wini[i] * (1-2*MSEQ[i/(SAMPLES_PER_CHIP/DOWN_SAMPLE)]);
                peakq[2] += winq[i] * (1-2*MSEQ[i/(SAMPLES_PER_CHIP/DOWN_SAMPLE)]);
            }
            peaka[2] = amp(peaki[2], peakq[2]);

            float capture_thresh = 0;
            for (int i=CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE/10*9; i<CHIPS*SAMPLES_PER_CHIP/DOWN_SAMPLE; i++) {
                capture_thresh += amp(wini[i], winq[i]);
            }

            // cout << peaka[2] << ' ' << capture_thresh << endl;
            // continue;

            if (capture_thresh>0.2) {
                capture_thresh *= 5; // fixed
                if (peaka[1]>capture_thresh && peaka[1]>=peaka[2] && peaka[1]>=peaka[0]) {
                    init_scale_rotate(peaka[1], peaki[1], peakq[1]);
                    break;
                }
            }
        }
    }

    // skip bubble
    for (int i=0; i<SAMPLES_PER_CHIP-DOWN_SAMPLE; i++) {
        float x = q.read();
        filteri(x * COS[i & 7] / 2);
        filterq(x * -SIN[i & 7] / 2);
    }

    // digital modem scheme
    char scheme = 0;
    for (int i=0; i<SCHEME_BITS; i++) {
        float consteli = 0;
        float constelq = 0;
        for (int j=0; j<SAMPLES_PER_BIT; j++) {
            float x = q.read();
            consteli += filteri(x * COS[j & 7] / 2);
            constelq += filterq(x * -SIN[j & 7] / 2);
        }
        scale_rotate(consteli, constelq);

        if (consteli > 0) {
            scheme &= ~(1<<i); // 1'b0
        } else {
            scheme |= (1<<i); // 1'b1
        }
    }

    // packet length
    unsigned len = 0;
    for (int i=0; i<LENGTH_BITS; i++) {
        float consteli = 0;
        float constelq = 0;
        for (int j=0; j<SAMPLES_PER_BIT; j++) {
            float x = q.read();
            consteli += filteri(x * COS[j & 7] / 2);
            constelq += filterq(x * -SIN[j & 7] / 2);
        }
        scale_rotate(consteli, constelq);

        if (consteli > 0) {
            len &= ~(1<<i); // 1'b0
        } else {
            len |= (1<<i); // 1'b1
        }
    }

    if (len<len_limit) {
        len_limit=len;
    }

    // data processing
    for (unsigned i=0; i<len; i++) {
        char octet = 0;

        // only 2psk is implemented
        for (int j=0; j<8; j++) {
            float consteli = 0;
            float constelq = 0;
            for (int k=0; k<SAMPLES_PER_BIT; k++) {
                float x = q.read();
                consteli += filteri(x * COS[k & 7] / 2);
                constelq += filterq(x * -SIN[k & 7] / 2);
            }
            scale_rotate(consteli, constelq);

            cout << consteli << ' ' << constelq << endl;

            if (consteli > 0) {
                octet &= ~(1<<j); // 1'b0
            } else {
                octet |= (1<<j); // 1'b1
            }
        }

        if (i<len_limit) {
            data[i] = octet;
        }
    }
}

void ui(void) {
    init_filter();

    cerr << "Listening for data..." << endl;
    char s[256];
    unsigned len;

    while (true) {
        len = sizeof(s)-1;
        rx_demodulate(s, len);
        s[len] = '\0';
        if (len==0) {
            cerr << "Goodbye" << endl;
            break;
        } else {
            cout << s << endl;
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
    unsigned int bufferFrames = 1024;

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
