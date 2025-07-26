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

// fir lpf: kaiser win, fs=48k, fpass=750, fstop=1250, ripple=0.8db, attenuation=26.8db
static const float LPF[128] = {0.00248450157232583,0.00272126239724457,0.00292312563396990,0.00308355083689094,0.00319636939093471,0.00325592979788780,0.00325723853893578,0.00319609697908163,0.00306923151947558,0.00287441280670464,0.00261056004092097,0.00227783294394612,0.00187770312186331,0.00141300668474287,0.000887974747456610,0.000308240443700925,-0.000319177925121039,-0.000985919148661196,-0.00168234610464424,-0.00239763851277530,-0.00311990897171199,-0.00383634306490421,-0.00453336024656892,-0.00519679579883814,-0.00581209734082222,-0.00636453879997134,-0.00683944253250957,-0.00722241215407848,-0.00749956863000989,-0.00765778357163072,-0.00768491486087441,-0.00757003063336015,-0.00730362255126238,-0.00687780370935798,-0.00628649070858955,-0.00552555732429028,-0.00459296908229590,-0.00348888593725860,-0.00221573514863849,-0.000778252375312150,0.000816511455923319,0.00255921809002757,0.00443830434232950,0.00644007837399840,0.00854885391891003,0.0107471076771617,0.0130156790837646,0.0153339887037873,0.0176802836358547,0.0200319066643715,0.0223655831068754,0.0246577113866806,0.0268846843391657,0.0290231890976429,0.0310505200177431,0.0329448916018009,0.0346857346594334,0.0362539626657963,0.0376322567462921,0.0388052947819233,0.0397599637508392,0.0404855459928513,0.0409738682210445,0.0412194170057774,0.0412194170057774,0.0409738682210445,0.0404855459928513,0.0397599637508392,0.0388052947819233,0.0376322567462921,0.0362539626657963,0.0346857346594334,0.0329448916018009,0.0310505200177431,0.0290231890976429,0.0268846843391657,0.0246577113866806,0.0223655831068754,0.0200319066643715,0.0176802836358547,0.0153339887037873,0.0130156790837646,0.0107471076771617,0.00854885391891003,0.00644007837399840,0.00443830434232950,0.00255921809002757,0.000816511455923319,-0.000778252375312150,-0.00221573514863849,-0.00348888593725860,-0.00459296908229590,-0.00552555732429028,-0.00628649070858955,-0.00687780370935798,-0.00730362255126238,-0.00757003063336015,-0.00768491486087441,-0.00765778357163072,-0.00749956863000989,-0.00722241215407848,-0.00683944253250957,-0.00636453879997134,-0.00581209734082222,-0.00519679579883814,-0.00453336024656892,-0.00383634306490421,-0.00311990897171199,-0.00239763851277530,-0.00168234610464424,-0.000985919148661196,-0.000319177925121039,0.000308240443700925,0.000887974747456610,0.00141300668474287,0.00187770312186331,0.00227783294394612,0.00261056004092097,0.00287441280670464,0.00306923151947558,0.00319609697908163,0.00325723853893578,0.00325592979788780,0.00319636939093471,0.00308355083689094,0.00292312563396990,0.00272126239724457,0.00248450157232583};
#define LPF_LEN (sizeof(LPF)/sizeof(LPF[0]))

static fifo<float> q; // inter-thread data queue

static fir_filter lpf;

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

    // preamble0: chirp for time & phase corrector
    for (int i=0; i<PREAM_BODY; i++) {
        constel.I = chirp(i);
        constel.Q = 0;
        sample = lpf.filter(constel);
        q.write(COS[i&7]*sample.I - SIN[i&7]*sample.Q);
    }

    // bubble: avoid payload-preamble interference
    constel.I = 0;
    constel.Q = 0;
    for (int i=0; i<BUBBLE_BODY; i++) {
        sample = lpf.filter(constel);
        q.write(COS[i&7]*sample.I - SIN[i&7]*sample.Q);
    }

    // // frame length in octets
    // len &= ~((-1)<<LENGTH_BITS);
    // for (int j=0; j<LENGTH_BITS; j++) {
    //     int bit = (len>>j) & 1; // lsb first
    //     constel.I = (1-2*bit);
    //     constel.Q = 0;
    //     for (int i=0; i<SYMBOL_BODY; i++) {
    //         sample = lpf.filter(constel);
    //         q.write(COS[i&7]*sample.I - SIN[i&7]*sample.Q);
    //     }
    // }

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
