#pragma once
#include <immintrin.h>

// oscillator for frequency mixing
static const float SIN[8] = {0, 0.353553390593274, -0.5, 0.353553390593274, 0, -0.353553390593273, 0.5, -0.353553390593274};
static const float COS[8] = {0.5, -0.353553390593274, 0, 0.353553390593274, -0.5, 0.353553390593275, 0, -0.353553390593274};
static const double PI = 3.1415926535897932384626433;

// fir lpf: eq ripple, fs=48k, fpass=650, fstop=1950, ripple=1db, attenuation=113db
static const float B[128] = {
    -6.37295807e-006, -9.64402807e-006, -1.66606351e-005, -2.67653431e-005, -4.07709485e-005, -5.96181817e-005, -8.43209491e-005, -1.15963798e-004, -1.55676418e-004, -2.04593045e-004, -2.63806258e-004, -3.34315671e-004, -4.16969153e-004, -5.12393366e-004, -6.20917999e-004, -7.42498960e-004, -8.76639446e-004, -1.02230883e-003, -1.17786461e-003, -1.34098053e-003, -1.50858436e-003, -1.67680834e-003, -1.84095348e-003, -1.99547131e-003, -2.13396898e-003, -2.24923855e-003, -2.33331113e-003, -2.37753754e-003, -2.37269676e-003, -2.30913283e-003, -2.17691949e-003, -1.96604966e-003, -1.66664540e-003, -1.26918685e-003, -7.64756172e-004, -1.45289450e-004, 5.96167054e-004, 1.46520650e-003, 2.46582087e-003, 3.60017805e-003, 4.86842263e-003, 6.26850734e-003, 7.79606588e-003, 9.44432337e-003, 1.12040592e-002, 1.30636226e-002, 1.50089944e-002, 1.70239136e-002, 1.90900564e-002, 2.11872607e-002, 2.32938118e-002, 2.53867656e-002, 2.74423212e-002, 2.94362064e-002, 3.13441083e-002, 3.31420936e-002, 3.48070748e-002, 3.63172106e-002, 3.76523510e-002, 3.87944020e-002, 3.97277139e-002, 4.04393598e-002, 4.09194119e-002, 4.11611348e-002, 4.11611348e-002, 4.09194119e-002, 4.04393598e-002, 3.97277139e-002, 3.87944020e-002, 3.76523510e-002, 3.63172106e-002, 3.48070748e-002, 3.31420936e-002, 3.13441083e-002, 2.94362064e-002, 2.74423212e-002, 2.53867656e-002, 2.32938118e-002, 2.11872607e-002, 1.90900564e-002, 1.70239136e-002, 1.50089944e-002, 1.30636226e-002, 1.12040592e-002, 9.44432337e-003, 7.79606588e-003, 6.26850734e-003, 4.86842263e-003, 3.60017805e-003, 2.46582087e-003, 1.46520650e-003, 5.96167054e-004, -1.45289450e-004, -7.64756172e-004, -1.26918685e-003, -1.66664540e-003, -1.96604966e-003, -2.17691949e-003, -2.30913283e-003, -2.37269676e-003, -2.37753754e-003, -2.33331113e-003, -2.24923855e-003, -2.13396898e-003, -1.99547131e-003, -1.84095348e-003, -1.67680834e-003, -1.50858436e-003, -1.34098053e-003, -1.17786461e-003, -1.02230883e-003, -8.76639446e-004, -7.42498960e-004, -6.20917999e-004, -5.12393366e-004, -4.16969153e-004, -3.34315671e-004, -2.63806258e-004, -2.04593045e-004, -1.55676418e-004, -1.15963798e-004, -8.43209491e-005, -5.96181817e-005, -4.07709485e-005, -2.67653431e-005, -1.66606351e-005, -9.64402807e-006, -6.37295807e-006, 
};
static const int ORDER = sizeof(B)/sizeof(B[0])-1;

// zc sequence

// maximum-length sequence for frequency spreading/despreading
static const char MSEQ[63] = {1,0,0,0,0,0,1,1,1,1,0,0,1,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,0,0,1,0,0,0,1,0,1,1,0,1,1,1,1,1,1,0,1,0,1,1,1,0,0,0,1,1,0,0,1,1,1,0,1};
static const int CHIPS = sizeof(MSEQ)/sizeof(MSEQ[0]);

// misc
static const int SAMPLE_RATE = 48000;
static const int CHANNELS = 2;
static const int SAMPLES_PER_CHIP = 256;
static const int DOWN_SAMPLE = SAMPLES_PER_CHIP/16;
static const int SAMPLES_PER_BIT = 256;
static const float TX_PA = 1.5;
static const int LENGTH_BITS = 7;
static const int SCHEME_BITS = 1;
static const int BUF_DEPTH = 512;   // common divisor of (SAMPLES_PER_CHIP*(CHIPS+1)) and samples per byte (2 symbols minimum can be a byte)
                                    // >ORDER
                                    // cannot be too small (e.g., <256) in case of overflow/underflow

// modem scheme
enum {
    PSK_2,
    QAM_16
};

// sample
typedef struct {
    float left;
    float right;
} sample_t;

// filter implementation
#if defined(__AVX__)
typedef union {
    __m256  d;
    __m256i di;
    float   f[8];
    int     i[8]; // only 32-bit int supported
} packed_t;

static __attribute__((aligned(32))) float b[(ORDER+1)/8][8];
static __attribute__((aligned(32))) float zi[(ORDER+1)/8][8];
static __attribute__((aligned(32))) float zq[(ORDER+1)/8][8];
#else
static float zi[ORDER+1];
static float zq[ORDER+1];
#endif

static inline void init_filter(void) {
#if defined(__AVX__)
    packed_t vindex = {
        .i = {0,16,32,48,64,80,96,112}
    };

    for (int i=0; i<(ORDER+1)/8; i++) {
        ((packed_t*)(zi+i))->d = _mm256_setzero_ps();
        ((packed_t*)(zq+i))->d = _mm256_setzero_ps();
        ((packed_t*)(b+i))->d = _mm256_i32gather_ps(B+i, vindex.di, 4); // reshape coef
    }
#else
    memset(zi, 0, sizeof(zi));
    memset(zq, 0, sizeof(zq));
#endif
}

static inline float filteri(float x) {
    static int i = 0;
#if defined(__AVX__)
    packed_t s = {
        .f = {0}
    };

    ((packed_t*)(zi+i))->d = _mm256_permute_ps(((packed_t*)(zi+i))->d, 0x39);
    zi[i][3] = zi[i][7];
    zi[i][7] = x;

    i = (i+1)%(sizeof(zi)/sizeof(zi[0]));
    for (int j=0; j<(sizeof(zi)/sizeof(zi[0])); j++) {
        s.d = _mm256_fmadd_ps(((packed_t*)(zi+(j+i)%(sizeof(zi)/sizeof(zi[0]))))->d, ((packed_t*)(b+j))->d, s.d);
    }

    packed_t s_tmp;
    s_tmp.d = _mm256_hadd_ps(s.d, s.d);
    s.d = _mm256_hadd_ps(s_tmp.d, s_tmp.d); // or just return s_tmp.f[2] + [3] + [4] + [5]
    return (s.f[3] + s.f[4]);
#else // non-AVX

#endif
}

static inline float filterq(float x) {
    static int i = 0;
#if defined(__AVX__)
    packed_t s = {
        .f = {0}
    };

    ((packed_t*)(zq+i))->d = _mm256_permute_ps(((packed_t*)(zq+i))->d, 0x39);
    zq[i][3] = zq[i][7];
    zq[i][7] = x;

    i = (i+1)%(sizeof(zq)/sizeof(zq[0]));
    for (int j=0; j<(sizeof(zq)/sizeof(zq[0])); j++) {
        s.d = _mm256_fmadd_ps(((packed_t*)(zq+(j+i)%(sizeof(zq)/sizeof(zq[0]))))->d, ((packed_t*)(b+j))->d, s.d);
    }

    packed_t s_tmp;
    s_tmp.d = _mm256_hadd_ps(s.d, s.d);
    s.d = _mm256_hadd_ps(s_tmp.d, s_tmp.d); // or just return s_tmp.f[2] + [3] + [4] + [5]
    return (s.f[3] + s.f[4]);
#else // non-AVX

#endif
}
