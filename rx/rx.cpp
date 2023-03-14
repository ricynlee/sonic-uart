#include "RtAudio.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "ThreadSafeQueue.hpp"

/*
const signed short SIN[8] = {0, 11585, -16384, 11585, 0, -11585, 16384, -11585};
const signed short COS[8] = {16384, -11585, 0, 11585, -16384, 11585, 0, -11585};
const int KOSC = 14;

const signed short MSEQ[] = {0, -1, -1, -1, -1, 0, -1, 0, -1, 0, 0, 0, -1, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, -1, -1, 0, 0, -1, 0, -1, -1};

#define ORDER 4

const int B[ORDER+1] = {5416, 21664, 32496, 21664, 5416};
const int KB = 26;

const int A[ORDER+1] ={4096, -14226, 18652, -10931, 2415};
const int KA = 12;

#define N_SAMPLES_PER_CHIP 128
#define N_CHIPS_PER_BIT 31

#define N_ACC 16

#define DECISION_THRESH 65536
*/

const float SIN[8] = {0,0.7071,-1.0000,0.7071,0,-0.7071,1.0000,-0.7071};
const float COS[8] = {1.0000,-0.7071,0,0.7071,-1.0000,0.7071,0,-0.7071};

const float MSEQ[] = {1, -1, -1, -1, -1, 1, -1, 1, -1, 1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1};

#define ORDER 4

const float B[ORDER+1] = {8.07046890258789e-05,0.000322818756103516,0.000484228134155273,0.000322818756103516,8.07046890258789e-05};
const float A[ORDER+1] ={1.0000,-3.4731,4.5537,-2.6687,0.5896};

#define N_SAMPLES_PER_CHIP 128
#define N_CHIPS_PER_BIT 31

#define N_ACC 16

#define DECISION_THRESH 65536


int callback( void* /* outputBuffer */, void* inputBuffer, unsigned int nBufferFrames,
         double /* streamTime */, RtAudioStreamStatus status, void* userData )
{
  static float z[ORDER] = {0};

  if (status) std::cout << "Overflow!" << std::endl;

  ThreadSafeQueue<float>* pq = (ThreadSafeQueue<float>*)userData;
  float* buffer = (float*) inputBuffer;
  float x0, x1, x2;

  for (int i=0; i<nBufferFrames; i++) {
    x0 = buffer[i];

    // mix - I
    x1 = x0*COS[i & 0x7];

    // iir lpf - I
    x2 = z[0] + (x1*B[0]);
    for (int n=1; n<ORDER; n++) {
      z[n-1] = z[n] + (x1*B[n]) - (x2*A[n]);
    }
    z[ORDER-1] = (x1*B[ORDER]) - (x2*A[ORDER]);

    // store - I
    pq->push(x2);

    // mix - Q
    x1 = x0*SIN[(i ^ 0x4) & 0x7]; // -sin

    // iir lpf - Q
    x2 = z[0] + (x1*B[0]);
    for (int n=1; n<ORDER; n++) {
      z[n-1] = z[n] + (x1*B[n]) - (x2*A[n]);
    }
    z[ORDER-1] = (x1*B[ORDER]) - (x2*A[ORDER]);

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
    std::cout << "\nNo audio devices found!\n";
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
    adc.openStream( NULL, &parameters, RTAUDIO_FLOAT32,
                    sampleRate, &bufferFrames, &callback, (void *)&q);
    adc.startStream();
  }
  catch ( RtAudioError& e ) {
    e.printMessage();
    exit( 0 );
  }
  
  float i0, q0, x0, y0;
  float y[N_CHIPS_PER_BIT*N_SAMPLES_PER_CHIP/N_ACC] = {0};
  float p[3];
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

      x0 = sqrt(i0*i0 + q0*q0);

      y0 += x0;
    }

    std::cout << y0 << std::endl;
    continue;

    for (int i=0; i<N_CHIPS_PER_BIT*N_SAMPLES_PER_CHIP/N_ACC-1; i++) {
      y[i] = y[i+1];
    }
    y[N_CHIPS_PER_BIT*N_SAMPLES_PER_CHIP/N_ACC-1] = y0;

    p[0] = p[1];
    p[1] = p[2];
    p[2] = 0;
    for (int i=0; i<N_CHIPS_PER_BIT*N_SAMPLES_PER_CHIP/N_ACC; i++) {
      p[2] += y[i] * MSEQ[i >> 3];
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
      std::cout << c;
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
