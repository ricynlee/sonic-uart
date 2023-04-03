#pragma once
#include <immintrin.h>

// oscillator for frequency mixing
static const float SIN[8] = {0, 0.353553390593274, -0.5, 0.353553390593274, 0, -0.353553390593273, 0.5, -0.353553390593274};
static const float COS[8] = {0.5, -0.353553390593274, 0, 0.353553390593274, -0.5, 0.353553390593275, 0, -0.353553390593274};
static const double PI = 3.1415926535897932384626433;

// fir lpf: eq ripple, fs=48k, fpass=750, fstop=2000, ripple=1db, attenuation=108db
static const float B[128] = {-7.33712295186706e-06, -1.19098694995046e-05, -2.12270824704319e-05, -3.49923466274049e-05, -5.45374241482932e-05, -8.14190061646514e-05, -0.000117392526590265, -0.000164412282174453, -0.000224599672947079, -0.000300198968034238, -0.000393526104744524, -0.000506901065818965, -0.000642562052235007, -0.000802567170467228, -0.000988686457276344, -0.00120228179730475, -0.00144417630508542, -0.00171451957430691, -0.00201265420764685, -0.00233698450028896, -0.00268485234118998, -0.00305242580361664, -0.00343461008742452, -0.00382498092949390, -0.00421574804931879, -0.00459775142371655, -0.00496049737557769, -0.00529223727062345, -0.00558009184896946, -0.00581021420657635, -0.00596800679340959, -0.00603837566450238, -0.00600602990016341, -0.00585581175982952, -0.00557305989786983, -0.00514399213716388, -0.00455610221251845, -0.00379856070503593, -0.00286260526627302, -0.00174191326368600, -0.000432942673796788, 0.00106476759538054, 0.00274835107848048, 0.00461143022403121, 0.00664400262758136, 0.00883240997791290, 0.0111593836918473, 0.0136041855439544, 0.0161428283900023, 0.0187483895570040, 0.0213913954794407, 0.0240402966737747, 0.0266619957983494, 0.0292224306613207, 0.0316872112452984, 0.0340222641825676, 0.0361945182085037, 0.0381725504994392, 0.0399272404611111, 0.0414323620498180, 0.0426651537418366, 0.0436067767441273, 0.0442427396774292, 0.0445632115006447, 0.0445632115006447, 0.0442427396774292, 0.0436067767441273, 0.0426651537418366, 0.0414323620498180, 0.0399272404611111, 0.0381725504994392, 0.0361945182085037, 0.0340222641825676, 0.0316872112452984, 0.0292224306613207, 0.0266619957983494, 0.0240402966737747, 0.0213913954794407, 0.0187483895570040, 0.0161428283900023, 0.0136041855439544, 0.0111593836918473, 0.00883240997791290, 0.00664400262758136, 0.00461143022403121, 0.00274835107848048, 0.00106476759538054, -0.000432942673796788, -0.00174191326368600, -0.00286260526627302, -0.00379856070503593, -0.00455610221251845, -0.00514399213716388, -0.00557305989786983, -0.00585581175982952, -0.00600602990016341, -0.00603837566450238, -0.00596800679340959, -0.00581021420657635, -0.00558009184896946, -0.00529223727062345, -0.00496049737557769, -0.00459775142371655, -0.00421574804931879, -0.00382498092949390, -0.00343461008742452, -0.00305242580361664, -0.00268485234118998, -0.00233698450028896, -0.00201265420764685, -0.00171451957430691, -0.00144417630508542, -0.00120228179730475, -0.000988686457276344, -0.000802567170467228, -0.000642562052235007, -0.000506901065818965, -0.000393526104744524, -0.000300198968034238, -0.000224599672947079, -0.000164412282174453, -0.000117392526590265, -8.14190061646514e-05, -5.45374241482932e-05, -3.49923466274049e-05, -2.12270824704319e-05, -1.19098694995046e-05, -7.33712295186706e-06};
static const int ORDER = sizeof(B)/sizeof(B[0])-1;

// maximum-length sequence for frequency spreading/despreading
static const char MSEQ[63] = {1,0,0,0,0,0,1,1,1,1,0,0,1,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,0,0,1,0,0,0,1,0,1,1,0,1,1,1,1,1,1,0,1,0,1,1,1,0,0,0,1,1,0,0,1,1,1,0,1};
static const int CHIPS = sizeof(MSEQ)/sizeof(MSEQ[0]);

// misc
static const int SAMPLE_RATE = 48000;
static const int MULTIPATH_MITIG = 128;  // 32 | 64 is allowed
static const int SAMPLES_PER_CHIP = 512;
static const int DOWN_SAMPLE = SAMPLES_PER_CHIP/16;
static const int SAMPLES_PER_SYM = 1024;
static const float TX_PA = 1.5;
static const int LENGTH_BITS = 7;
static const int SCHEME_BITS = 1;
static const int TX_BUF_DEPTH = 512;    // common divisor of (SAMPLES_PER_CHIP*(CHIPS+1)) and samples per byte (2 symbols minimum can be a byte)
                                        // >ORDER
                                        // cannot be too small (e.g., <256) in case of overflow/underflow

// modem scheme
enum {
    PSK_2,
    PSK_4
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
