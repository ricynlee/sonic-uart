#include <cstring>
#include <immintrin.h>
#include <cmath>
#include "dsp.hpp"

using namespace std;

#if defined(__AVX__)
typedef union alignas(32) {
    __m256  d;
    __m256i di;
    float   f[8];
    int     i[8]; // only 32-bit int supported
} packed_t;

typedef struct {
    int       n;
    int       i;
    packed_t* b;
    packed_t* z;
} fir_filter_data_t;

#define Dn (((fir_filter_data_t*)data)->n)
#define Di (((fir_filter_data_t*)data)->i)
#define Db (((fir_filter_data_t*)data)->b)
#define Dz (((fir_filter_data_t*)data)->z)
#else
typedef struct {
    int   n;
    int   i;
    float (*b)[8];
    float (*z)[8];
} fir_filter_data_t;
#endif

fir_filter::fir_filter(const float* const coef, int order) {
    // coef length = 4n
    // order = 4n-1
#if defined(__AVX__)
    // assert((order+1) % 4 == 0);
    static const packed_t vindex = {
        .i = {0,32,64,96,0,32,64,96}
    };
    data = malloc(sizeof(fir_filter_data_t));
    Dn = order+1;
    Di = 0;
    Db = (packed_t*)_mm_malloc(Dn/4*8*2*sizeof(float), 32);
    Dz = Db + Dn/4;
    for (int i=0; i<Dn/4; i++) {
        Dz[i].d = _mm256_setzero_ps();
        Db[i].d = _mm256_i32gather_ps(coef+i, vindex.di, 4); // reshape coef
    }
#else
    // non-AVX routine not implemented
#endif
}

fir_filter::~fir_filter() {
    _mm_free(Db);
    free(data);
}

sample_t fir_filter::filter(const sample_t& in) {
#if defined(__AVX__)
    static packed_t summed;
    summed.d = _mm256_setzero_ps();

    Dz[Di].d = _mm256_permute_ps(Dz[Di].d, 0x39);
    Dz[Di].f[3] = in.I;
    Dz[Di].f[7] = in.Q;

    Di = (Di+1)%(Dn/4);
    for (int j=0; j<Dn/4; j++) {
        summed.d = _mm256_fmadd_ps(Dz[(Di+j)%(Dn/4)].d, Db[j].d, summed.d);
    }

    summed.d = _mm256_hadd_ps(summed.d, summed.d);
    summed.d = _mm256_hadd_ps(summed.d, summed.d);

    // out.I = summed.f[3];
    // out.Q = summed.f[4];

    return *(sample_t*)(summed.f+3);
#else
    // non-AVX routine not implemented
#endif
}

// chirp generator
float chirp(size_t index) {
    const double u = 1200. * SAMPLE_RATE / CHIRP_BODY / 2; 
    double t = (double)index/SAMPLE_RATE;
    return (float)cos(2*PI*(-600*t+u*t*t));
}
