#pragma once

#define PI 3.1415926535897932384626433

// misc
#define SAMPLE_RATE     48000
#define CARRIER_FRQ     18000

// noise measuring
#define NOISE_BODY      16384

// chirp as preamble
#define PREAM_BODY      65536
#define BUBBLE_BODY     256
#define CARRIER_BODY    (65536-BUBBLE_BODY)

// symbol
#define SYMBOL_BODY     512

#define LENGTH_BITS     8

#define BPSK 1
#define QPSK 2
#define QAM16 4

#define MODEM BPSK

// typedefs
typedef struct {
    union {
        float L;
        float I;
    };
    union {
        float R;
        float Q;
    };
} sample_t;

// declarations
class fir_filter {
public:
    fir_filter();
    ~fir_filter();
    void init(const float* const, int); // initialize coefficients
    void clear(); // clear z
    sample_t filter(const sample_t&);
private:
    void* data;
};

class biquad_filter {
private:
    biquad_filter();
    ~biquad_filter();
    void init(const float*);
    void clear();
    sample_t filter(const sample_t&);
private:
    void* data;
};

float chirp(size_t);
