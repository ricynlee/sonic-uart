#include <cstring>
#include <immintrin.h>
#include <cmath>
#include <cstdlib>
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
    packed_t* b;
    packed_t* z;
} fir_filter_data_t;
#else
typedef struct {
    int     n;
    int     i;
    float*  b;
    float*  z;
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
    data->n = order+1;
    data->i = 0;
    data->b = _mm_malloc(data->n * 2 * sizeof(float), sizeof(packed_t));
    data->z = data->b + data->n;
    for (int i=0; i<data->n/4; i++) {
        (data->z+i)->d = _mm256_setzero_ps();
        (data->b+i)->d = _mm256_i32gather_ps(coef+i, vindex.di, 4); // reshape coef
    }
#else
    // non-AVX routine not implemented
#endif
}

fir_filter::~fir_filter() {
    _mm_free(data->b);
    free(data);
}

sample_t fir_filter::filter(const sample_t& in) {
#if defined(__AVX__)
    static packed_t summed;
    summed.d = _mm256_setzero_ps();

    (data->z+data->i)->d = _mm256_permute_ps(data->z+data->i)->d, 0x39);
    z[data->i][3] = in.I;
    z[data->i][7] = in.Q;

    data->i = (data->i+1) % (data->n/4);
    for (int j=0; j<data->n/4; j++) {
        summed.d = _mm256_fmadd_ps((data->z + (data->i+j)%(data->n/4))->d, (data->b + j)->d, summed.d);
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
float chirp(size_t index, bool decimated) {
    const double u = 1200. * SAMPLE_RATE / CHIRP_BODY / 2; 
    double t = index*(double)(decimated ? DECIMATION : 1)/SAMPLE_RATE;
    return (float)cos(2*PI*(-600*t+u*t*t));
}
