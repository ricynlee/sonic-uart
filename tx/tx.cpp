#include "RtAudio.h"
#include <iostream>
#include "ThreadSafeQueue.hpp"
#include "dsp.hpp"

using namespace std;

int callback( void* out_buf, void* /* in_buf */, unsigned samples, double /* timestamp */, RtAudioStreamStatus status, void* shared_data) {
  static unsigned int nTxFrameRemainder = SAMPLES_PER_CHIP*CHIPS_PER_BIT;
  static int z[ORDER] = {0};

  if (status) std::cout << "Underflow!" << std::endl;

  ThreadSafeQueue<signed short>* pq = (ThreadSafeQueue<signed short>*)userData;
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

    // iir filter
    Y = z[0] + ((X*B[0])>>KB);
    for (int n=1; n<ORDER; n++) {
      z[n-1] = z[n] + ((X*B[n])>>KB) - ((Y*A[n])>>KA);
    }
    z[ORDER-1] = ((X*B[ORDER])>>KB) - ((Y*A[ORDER])>>KA);

    buffer[2*i+0] = Y;
    buffer[2*i+1] = Y;
  }

  return 0;
}

int main()
{

  RtAudio dac;
  if ( dac.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 0 );
  }
  RtAudio::StreamParameters parameters;
  parameters.deviceId = dac.getDefaultOutputDevice();
  parameters.nChannels = 2;
  parameters.firstChannel = 0;
  unsigned int sampleRate = 48000;
  unsigned int bufferFrames = 2048;

  ThreadSafeQueue<signed short> q;

  try {
    dac.openStream( &parameters, NULL, RTAUDIO_SINT16,
                    sampleRate, &bufferFrames, &callback, (void *)&q);
    dac.startStream();
  }
  catch ( RtAudioError& e ) {
    e.printMessage();
    exit( 0 );
  }
  
  char c;
  std::cout << ">> ";
  while (true) {
    if (std::cin.get(c)) {
      for (int i=0; i<8; i++) {
        for (int j=0; j<N_CHIPS_PER_BIT; j++) {
          for (int k=0; k<N_SAMPLES_PER_CHIP; k++) {
            q.push((((c>>i) & 0x1) ^ MSEQ[j]) ? OSC[k & 0x7] : 0);
          }
        }
      }
      if (c=='\n') {
        std::cout << ">> ";
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
