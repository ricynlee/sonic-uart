#pragma once

#include <cstdint>

// parameter that can be overriden in makefile
#define SAMPLES_PER_CHIP 128 // for frequency spreading, 64~512
#define CAPTURE_THRESH (6400000/SAMPLES_PER_CHIP)
#define DECISION_THRESH (32768)

// sampling
static const int SAMPLE_RATE = 48000;
static const int CHANNELS = 2;

// oscillator for frequency mixing
static const int16_t MOD_OSC_2PSK[2][8] = {
    {0, 11585, -16384, 11585, 0, -11585, 16384, -11585}, // 1'b0
    {0, -11585, 16384, -11585, 0, 11585, -16384, 11585}, // 1'b1
};
static const int16_t DEM_I_OSC[8] = {16384, -11585, 0, 11585, -16384, 11585, 0, -11585}; // cos
static const int16_t DEM_Q_OSC[8] = {0, -11585, 16384, -11585, 0, 11585, -16384, 11585}; // -sin
static const int OSC_FRAC_WID = 14; // fraction width of fixed point

static const double PI = 3.1415926535897932384626433;

// band-pass filter for modulation
static const int MOD_B[] = {7408, 0, -22223, 0, 22223, 0, -7408};
static const int MOD_B_FRAC_WID = 22;
static const int MOD_A[] = {2048, 8000, 15526, 17884, 13026, 5629, 1210};
static const int MOD_A_FRAC_WID = 11;
static const int MOD_FILT_ORDER = sizeof(MOD_B)/sizeof(MOD_B[0])-1;

// low-pass filter for demodulation
// static const int16_t DEM_B[] = {5416, 21664, 32496, 21664, 5416};
// static const int DEM_B_FRAC_WID = 26;
// static const int16_t DEM_A[] ={4096, -14226, 18652, -10931, 2415};
// static const int DEM_A_FRAC_WID = 12;
static const int16_t DEM_B[] = {9699, 29098, 29098, 9699};
static const int DEM_B_FRAC_WID = 27;
static const int16_t DEM_A[] = {8192, -23173, 21887, -6902};
static const int DEM_A_FRAC_WID = 13;

static const int DEM_FILT_ORDER = sizeof(DEM_B)/sizeof(DEM_B[0])-1;

// maximum-length sequence for frequency spreading/despreading
static const int8_t MSEQ[] = {
    1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1
};
static const int CHIPS = sizeof(MSEQ)/sizeof(MSEQ[0]);
static const int SAMPLES_PER_BIT = SAMPLES_PER_CHIP*CHIPS/8;

// accumulation parameters
static const int ASYNC_ACCUM = SAMPLES_PER_CHIP/16;
static const int SYNC_ACCUM = SAMPLES_PER_BIT;

// filtering
#define GEN_FILT(ORDER, Y, X, B, A, KB, KA) do {              \
    static int Z[ORDER] = {0};                            \
    Y = Z[0] + ((X*B[0])>>KB);                            \
    for (int n=1; n<ORDER; n++) {                         \
        Z[n-1] = Z[n] + ((X*B[n])>>KB) - ((Y*A[n])>>KA);  \
    }                                                     \
    Z[ORDER-1] = ((X*B[ORDER])>>KB) - ((Y*A[ORDER])>>KA); \
} while (0)
