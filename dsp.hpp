#pragma once

#include <cstddef>
#include <cstdint>

#define PI 3.1415926535897932384626433

// misc
#define SAMPLE_RATE     48000
#define CARRIER_FRQ     18000

// chirp as preamble
#define PREAM_BODY      65536
#define BUBBLE_BODY     256
#define CARRIER_BODY    65536

// symbol
#define SYMBOL_BODY     512

// fir lpf: kaiser win, fs=48k, fpass=797, fstop=1385, ripple=0.5db, attenuation=30db
#define LPF_COEF {3.63542589950387e-05, -0.000730017551613079, 0.00133737595152190, -0.000245040263572821, -0.00297663896384212, 0.00540885213623617, -0.00232317554913813, -0.00699323216224549, 0.0147957976547783, -0.00930739549292851, -0.0122374567576044, 0.0338324208427928, -0.0283433012586116, -0.0172910505482885, 0.0784896146037067, -0.0933995502596032, -0.0204144320174365, 0.560360875376853, 0.560360875376853, -0.0204144320174365, -0.0933995502596032, 0.0784896146037067, -0.0172910505482885, -0.0283433012586116, 0.0338324208427928, -0.0122374567576044, -0.00930739549292851, 0.0147957976547783, -0.00699323216224549, -0.00232317554913813, 0.00540885213623617, -0.00297663896384212, -0.000245040263572821, 0.00133737595152190, -0.000730017551613079, 3.63542589950387e-05}
#define LPF_LEN  36

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

float fast_exp(float);
