#include "RtAudio.h"
#include <iostream>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include "ThreadSafeQueue.hpp"

using namespace std;

inline int filter(int x) { // iir lpf for demodulation
    int y;
    GEN_FILT(DEM_FILT_ORDER, y, x, DEM_B, DEM_A, DEM_B_FRAC_WID, DEM_A_FRAC_WID);
    return y;
}


// an extra 0 added in the front
const int8_t MSEQ[] = {0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,1,1,1,0,1,1,1,1,0,0,0,1,0,1,1,0,0,1,1,0,1,1,0,0,0,0,1,1,1,1,0,0,1,1,1,0,0,0,0,1,0,1,0,1,1,1,1,1,1,1,1,0,0,1,0,1,1,1,1,0,1,0,0,1,0,1,0,0,0,0,1,1,0,1,1,1,0,1,1,0,1,1,1,1,1,0,1,0,1,1,1,0,1,0,0,0,0,0,1,1,0,0,1,0,1,0,1,0,1,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,1,1,0,1,1,0,1,0,1,0,0,1,1,0,1,0,0,1,1,1,1,1,1,0,1,1,1,0,0,1,1,0,0,1,1,1,1,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,1,0,0,1,0,0,1,1,0,0,0,1,0,0,1,1,1,0,1,0,1,0,1,1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,1,0,0,0,1,1,1,1,1};

// lpf
int32_t z[ORDER] = {0};


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









int callback( void* /* out_buf */, void* in_buf, unsigned samples,  double /* timestamp */, RtAudioStreamStatus status, void* shared_data) {
    if (status) std::cerr << "Overflow!" << std::endl;

    sample_t* buf = (sample_t*) in_buf;
    ThreadSafeQueue<int16_t>* pq = (ThreadSafeQueue<int16_t>*)shared_data;

    int avg;
    for (int i=0; i<samples; i++) {
        avg = buf[i].left + buf[i].right;
        avg = ((avg & 1) & ((avg >> 1) & 1)) + (avg >> 1); // improved round(mean(left, right))
        pq->push(avg);
    }

    return 0;
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
  parameters.nChannels = 2;
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
