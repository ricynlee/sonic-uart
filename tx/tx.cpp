#include "RtAudio.h"
#include <iostream>
#include "ThreadSafeQueue.hpp"

const signed short OSC[8] = {0, 11585, -16384, 11585, 0, -11585, 16384, -11585};
const signed char MSEQ[] = {1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0};

#define ORDER 6

const int B[ORDER+1] = {7408,0,-22223,0,22223,0,-7408};
const int KB = 22;

const int A[ORDER+1] ={2048,8000,15526,17884,13026,5629,1210};
const int KA = 11;

#define N_SAMPLES_PER_CHIP 128
#define N_CHIPS_PER_BIT 31

int callback( void* outputBuffer, void * /* inputBuffer */, unsigned int nBufferFrames,
         double /* streamTime */, RtAudioStreamStatus status, void *userData )
{
  static unsigned int nTxFrameRemainder = N_SAMPLES_PER_CHIP*N_CHIPS_PER_BIT;
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
