#include "RtAudio.h"
#include <iostream>
#include "ThreadSafeQueue.hpp"

const signed short OSC[8] = {0, 11585, -16384, 11585, 0, -11585, 16384, -11585};
const signed short MSEQ[] = {1, -1, -1, -1, -1, 1, -1, 1, -1, 1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1};

#define ORDER 6

const int B[ORDER+1] = {7408,0,-22223,0,22223,0,-7408};
const int KB = 22;

const int A[ORDER+1] ={2048,8000,15526,17884,13026,5629,1210};
const int KA = 11;

#define N_SAMPLES_PER_CHIP 128
#define N_CHIPS_PER_BIT 31

int callback( void* /* outputBuffer */, void* inputBuffer, unsigned int nBufferFrames,
         double /* streamTime */, RtAudioStreamStatus status, void* userData )
{
  static int z[ORDER] = {0};

  if (status) std::cout << "Overflow!" << std::endl;

  ThreadSafeQueue<signed short>* pq = (ThreadSafeQueue<signed short>*)userData;
  signed short* buffer = (signed short*) inputBuffer;
  int X, Y;

  for (int i=0; i<nBufferFrames; i++) {
    X = buffer[i];

    // iir filter
    Y = z[0] + ((X*B[0])>>KB);
    for (int n=1; n<ORDER; n++) {
      z[n-1] = z[n] + ((X*B[n])>>KB) - ((Y*A[n])>>KA);
    }
    z[ORDER-1] = ((X*B[ORDER])>>KB) - ((Y*A[ORDER])>>KA);

    pq->push(Y);
  }

  return 0;
}

int main()
{
  RtAudio adc;
  if ( adc.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 0 );
  }
  RtAudio::StreamParameters parameters;
  parameters.deviceId = adc.getDefaultInputDevice();
  parameters.nChannels = 1;
  parameters.firstChannel = 0;
  unsigned int sampleRate = 48000;
  unsigned int bufferFrames = 2048;

  ThreadSafeQueue<signed short> q;

  try {
    adc.openStream( NULL, &parameters, RTAUDIO_SINT16,
                    sampleRate, &bufferFrames, &callback, (void *)&q);
    adc.startStream();
  }
  catch ( RtAudioError& e ) {
    e.printMessage();
    exit( 0 );
  }
  
  int X;
  while (true) {
    while (q.empty());
    X = q.pop();
    // rectification
    if (X & 0x80000000) {
      X = -X;
    }
    // envelope detection


  }

  try {
    // Stop the stream
    adc.stopStream();
  }
  catch (RtAudioError& e) {
    e.printMessage();
  }
  if ( adc.isStreamOpen() ) adc.closeStream();
  return 0;
}
