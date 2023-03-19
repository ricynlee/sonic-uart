#include "RtAudio.h"
#include <cstdint>
#include <iostream>
#include "fifo.hpp"
#include "dsp.hpp"

using namespace std;

#define SAMPLES_TX_ATOMIC_UNIT 512 // common divisor of samples per DSSS code and samples per byte (2 symbols minimum can be a byte), cannot be too small (e.g., 128) so as not to under run

inline int tx_filter(int x) { // iir bpf for modulation
    int y;
    GEN_FILT(MOD_FILT_ORDER, y, x, MOD_B, MOD_A, MOD_B_FRAC_WID, MOD_A_FRAC_WID);
    return y;
}

int tx_callback( void* out_buf, void* /* in_buf */, unsigned samples, double /* timestamp */, RtAudioStreamStatus status, void* shared_data) {
  static int nTxFrameRemainder = SAMPLES_PER_CHIP*CHIPS_PER_BIT;

  if (status) cout << "Underflow!" << endl;

  fifo<int16_t>* pq = (fifo<int16_t>*)shared_data;
  signed short* buffer = (signed short*) outputBuffer;
  int X, Y;

  for (int i=0; i<nBufferFrames; i++) {
    if ((pq->size() >= nBufferFrames-i) || (pq->size() >= nTxFrameRemainder)) {
      X = pq->pop();
      if (nTxFrameRemainder==1) {
        nTxFrameRemainder = N_SAMPLES_PER_CHIP*N_CHIPS_PER_BIT;
      } else {
        nTxFrameRemainder--;
      }
    } else {
      X = 0;
    }

    y = filter();


    buffer[2*i+0] = y;
    buffer[2*i+1] = Y;
  }

  return 0;
}

int tx_process(fifo<signed short> q) {

}

int main()
{
    fifo<signed short> q;

    RtAudio dac;
    if (dev.getDeviceCount() < 1) {
        cerr << "No audio devices found!" << endl;
        return false;
    }
    RtAudio::StreamParameters parameters;
    parameters.deviceId = dac.getDefaultOutputDevice();
    parameters.nChannels = 2;
    parameters.firstChannel = 0;
    unsigned int sampleRate = 48000;
    unsigned int bufferFrames = SAMPLES_TX_ATOMIC_UNIT;

    try {
        dac.openStream(&parameters, NULL, RTAUDIO_SINT16, sampleRate, &bufferFrames, &tx_callback, (void *)&q);
        dac.startStream();
    }
    catch ( RtAudioError& e ) {
        e.printMessage();
        return (-1);
    }
  
  char c;
  cout << ">> ";
  while (true) {
    if (cin.get(c)) {
      for (int i=0; i<8; i++) {
        for (int j=0; j<N_CHIPS_PER_BIT; j++) {
          for (int k=0; k<N_SAMPLES_PER_CHIP; k++) {
            q.push((((c>>i) & 0x1) ^ MSEQ[j]) ? OSC[k & 0x7] : 0);
          }
        }
      }
      if (c=='\n') {
        cout << ">> ";
      }
    } else {
      break;
    }
  }

  try {
    // Stop the stream
    dac.stopStream();
  }
  catch (RtAudioError& e) {
    e.printMessage();
  }
  if ( dac.isStreamOpen() ) dac.closeStream();
  return 0;
}
