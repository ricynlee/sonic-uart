#include "RtAudio.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "ThreadSafeQueue.hpp"

// oscillator
const signed short SIN[8] = {0, -11585, 16384, -11585, 0, 11585, -16384, 11585}; // -sin
const signed short COS[8] = {16384, -11585, 0, 11585, -16384, 11585, 0, -11585};
const int KOSC = 14;

// an extra 0 added in the front
const signed char MSEQ[] = {0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,1,1,1,0,1,1,1,1,0,0,0,1,0,1,1,0,0,1,1,0,1,1,0,0,0,0,1,1,1,1,0,0,1,1,1,0,0,0,0,1,0,1,0,1,1,1,1,1,1,1,1,0,0,1,0,1,1,1,1,0,1,0,0,1,0,1,0,0,0,0,1,1,0,1,1,1,0,1,1,0,1,1,1,1,1,0,1,0,1,1,1,0,1,0,0,0,0,0,1,1,0,0,1,0,1,0,1,0,1,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,1,1,0,1,1,0,1,0,1,0,0,1,1,0,1,0,0,1,1,1,1,1,1,0,1,1,1,0,0,1,1,0,0,1,1,1,1,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,1,0,0,1,0,0,1,1,0,0,0,1,0,0,1,1,1,0,1,0,1,0,1,1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,1,0,0,0,1,1,1,1,1};

// lpf
const int B[] = {5416, 21664, 32496, 21664, 5416};
const int KB = 26;
const int A[] ={4096, -14226, 18652, -10931, 2415};
const int KA = 12;

#define N_SAMPLES_PER_CHIP 64
#define N_CHIPS_PER_BIT (sizeof(MSEQ)/sizeof(MSEQ[0]))
#define ORDER (sizeof(B)/sizeof(B[0])-1)
#define N_ACC 16
#define DECISION_THRESH 160000

int callback( void* /* outputBuffer */, void* inputBuffer, unsigned int nBufferFrames,
         double /* streamTime */, RtAudioStreamStatus status, void* userData )
{
  static int z[ORDER] = {0};

  if (status) std::cerr << "Overflow!" << std::endl;

  ThreadSafeQueue<int>* pq = (ThreadSafeQueue<int>*)userData;
  signed short* buffer = (signed short*) inputBuffer;
  int x0, x1, x2;

  for (int i=0; i<nBufferFrames; i++) {
    x0 = buffer[i];

    // mix - I
    x1 = (x0*COS[i & 0x7])>>KOSC;

    // iir lpf - I
    x2 = z[0] + ((x1*B[0])>>KB);
    for (int n=1; n<ORDER; n++) {
      z[n-1] = z[n] + ((x1*B[n])>>KB) - ((x2*A[n])>>KA);
    }
    z[ORDER-1] = ((x1*B[ORDER])>>KB) - ((x2*A[ORDER])>>KA);

    // store - I
    pq->push(x2);

    // mix - Q
    x1 = (x0*SIN[i & 0x7])>>KOSC; // -sin

    // iir lpf - Q
    x2 = z[0] + ((x1*B[0])>>KB);
    for (int n=1; n<ORDER; n++) {
      z[n-1] = z[n] + ((x1*B[n])>>KB) - ((x2*A[n])>>KA);
    }
    z[ORDER-1] = ((x1*B[ORDER])>>KB) - ((x2*A[ORDER])>>KA);

    // store - Q
    pq->push(x2);
  }

  return 0;
}

inline int approx_dist(int x, int y) {
  x = abs(x);
  y = abs(y);
  int min = x<y ? x : y;
  return (x + y - (min >> 1) - (min >> 2) + (min >> 4));
}

int main()
{
  RtAudio adc;
  if ( adc.getDeviceCount() < 1 ) {
    std::cerr << "\nNo audio devices found!\n";
    exit( 0 );
  }
  RtAudio::StreamParameters parameters;
  parameters.deviceId = adc.getDefaultInputDevice();
  parameters.nChannels = 1;
  parameters.firstChannel = 0;
  unsigned int sampleRate = 48000;
  unsigned int bufferFrames = 2048;

  ThreadSafeQueue<int> q;

  try {
    adc.openStream( NULL, &parameters, RTAUDIO_SINT16 /*RTAUDIO_FLOAT32*/,
                    sampleRate, &bufferFrames, &callback, (void *)&q);
    adc.startStream();
  }
  catch ( RtAudioError& e ) {
    e.printMessage();
    exit( 0 );
  }

  std::cerr << "Listening..." << std::endl;

  int i0, q0, x0, y0;
  int y[N_CHIPS_PER_BIT*N_SAMPLES_PER_CHIP/N_ACC] = {0};
  int p[3];
  int d;
  char c;
  int bitIdx = 0;
  while (true) {
    y0 = 0;
    for (int i=0; i<N_ACC; i++) {
      while (q.empty());
      i0 = q.pop();
      while (q.empty());
      q0 = q.pop();

      x0 = sqrt(i0*i0+q0*q0);

      y0 += x0/N_ACC;
    }

    // std::cout << y0 << std::endl;
    // continue;

    for (int i=0; i<N_CHIPS_PER_BIT*N_SAMPLES_PER_CHIP/N_ACC-1; i++) {
      y[i] = y[i+1];
    }
    y[N_CHIPS_PER_BIT*N_SAMPLES_PER_CHIP/N_ACC-1] = y0;

    p[0] = p[1];
    p[1] = p[2];
    p[2] = 0;
    for (int i=0; i<N_CHIPS_PER_BIT*N_SAMPLES_PER_CHIP/N_ACC; i++) {
      p[2] += MSEQ[i/(N_SAMPLES_PER_CHIP/N_ACC)] ? y[i] : -y[i];
    }

    std::cout << p[2] << std::endl;
    continue;

    if (p[1]>DECISION_THRESH && p[1]>=p[2] && p[1]>=p[3]) { // bit 0
      c &= ~(1<<bitIdx);
      bitIdx = (bitIdx+1) & 0x7;
      d = 1;
    } else if (p[1]<-DECISION_THRESH && p[1]<=p[2] && p[1]<=p[3]) { // bit 1
      c |= (1<<bitIdx);
      bitIdx = (bitIdx+1) & 0x7;
      d = -1;
    } else {
      d = 0;
    }

    if (d) {
      std::cout << "Got bit " << ((d & 0x80000000)>>31) << std::endl;
      if (bitIdx==0) {
        std::cout << c;
      }
    }
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
