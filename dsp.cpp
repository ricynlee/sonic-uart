#include <cstring>
#include <immintrin.h>
#include "dsp.hpp"

using namespace std;

// filter implementation
#if defined(__AVX__)
typedef union {
    __m256  d;
    __m256i di;
    float   f[8];
    int     i[8]; // only 32-bit int supported
} packed_t;

static __attribute__((aligned(32))) float b[(ORDER+1)/4][8];
static __attribute__((aligned(32))) float z[(ORDER+1)/4][8];
#endif

void init_filter(void) {
#if defined(__AVX__)
    packed_t vindex = {
        .i = {0,32,64,96,0,32,64,96}
    };

    for (int i=0; i<(ORDER+1)/4; i++) {
        ((packed_t*)(z+i))->d = _mm256_setzero_ps();
        ((packed_t*)(b+i))->d = _mm256_i32gather_ps(B+i, vindex.di, 4); // reshape coef
    }
#endif
}

sample_t filter(const sample_t& in) {
    sample_t out;

    static int i = 0;
#if defined(__AVX__)
    packed_t s = {
        .f = {0}
    };

    ((packed_t*)(z+i))->d = _mm256_permute_ps(((packed_t*)(z+i))->d, 0x39);
    z[i][3] = in.I;
    z[i][7] = in.Q;

    i = (i+1)%((ORDER+1)/4);
    for (int j=0; j<((ORDER+1)/4); j++) {
        s.d = _mm256_fmadd_ps(((packed_t*)(z+(j+i)%((ORDER+1)/4)))->d, ((packed_t*)(b+j))->d, s.d);
    }

    s.d = _mm256_hadd_ps(s.d, s.d);
    s.d = _mm256_hadd_ps(s.d, s.d);

    out.I = s.f[3];
    out.Q = s.f[4];

    return out;
#endif
}
