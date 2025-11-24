#include <cstring>
#include <cstdint>
#include <immintrin.h>
#include <cmath>
#include "dsp.hpp"

// #include <iostream> // for dbg purposes

using namespace std;

#if defined(__AVX__)
typedef union alignas(32) {
    __m256  d;
    __m256i di;
    float   f[8];
    int     i[8]; // only 32-bit int supported
} packed_t;
#endif

///////////////////////////////////////////////////////////////////////////////////////////////
// FIR filter
///////////////////////////////////////////////////////////////////////////////////////////////
#if defined(__AVX__)
typedef struct {
    int       n;
    int       i;
    packed_t* b;
    packed_t* z;
} fir_filter_data_t;
#else
typedef struct {
    int      n;
    int      i;
    float    *b;
    sample_t *z;
} fir_filter_data_t;
#endif

#define FIR_Dn (((fir_filter_data_t*)data)->n)
#define FIR_Di (((fir_filter_data_t*)data)->i)
#define FIR_Db (((fir_filter_data_t*)data)->b)
#define FIR_Dz (((fir_filter_data_t*)data)->z)

fir_filter::fir_filter() {
    data = NULL;
}

void fir_filter::init(const float* const coef, int len /* order+1 */) {
    // coef length = 4n
    // order = 4n-1
    if (data)
        return;
    data = malloc(sizeof(fir_filter_data_t));
    FIR_Dn = len;
    FIR_Di = 0;
#if defined(__AVX__)
    packed_t vindex;
    for (int i=0; i<8; i++) {
        vindex.i[i] = (i%4)*len/4;
    }

    FIR_Db = (packed_t*)_mm_malloc(FIR_Dn*4*sizeof(float), 32);
    FIR_Dz = FIR_Db + FIR_Dn*2*sizeof(float)/sizeof(packed_t);

    for (int i=0; i<FIR_Dn/4; i++) {
        FIR_Db[i].d = _mm256_i32gather_ps(coef+i, vindex.di, 4); // reshape coef
        FIR_Dz[i].d = _mm256_setzero_ps();
    }

    // for (int i=0; i<FIR_Dn/4; i++) {
    //     for (int j=0; j<8; j++)
    //         cout << FIR_Db[i].f[j] << ' ';
    //     cout << endl;
    // }
#else
    FIR_Db = (float*)malloc(FIR_Dn*sizeof(float) + FIR_Dn*sizeof(sample_t));
    FIR_Dz = (sample_t*)(FIR_Db + FIR_Dn);
    memset(FIR_Dz, 0, len*sizeof(sample_t));
    memcpy(FIR_Db, coef, len*sizeof(float));
    // for (int i=0; i<FIR_Dn; i++)
    //     cout << FIR_Db[i] << ' ';
    // cout << endl;
#endif
}

void fir_filter::clear() {
#if defined(__AVX__)
    for (int i=0; i<FIR_Dn/4; i++) {
        FIR_Dz[i].d = _mm256_setzero_ps();
    }
#else
    memset(FIR_Dz, 0, FIR_Dn*sizeof(float));
#endif
}

fir_filter::~fir_filter() {
#if defined(__AVX__)
    if (data) {
        _mm_free(FIR_Db);
        free(data);
    }
#else
    if (data) {
        free(FIR_Db);
        free(data);
    }
#endif
}

sample_t fir_filter::filter(const sample_t& in) {
#if defined(__AVX__)
    static packed_t summed;
    summed.d = _mm256_setzero_ps();

    FIR_Dz[FIR_Di].d = _mm256_permute_ps(FIR_Dz[FIR_Di].d, 0x93);
    FIR_Dz[FIR_Di].f[0] = in.I;
    FIR_Dz[FIR_Di].f[4] = in.Q;

    // cout << "sum([" << endl;
    for (int j=0; j<FIR_Dn/4; j++) {
        summed.d = _mm256_fmadd_ps(FIR_Dz[(FIR_Di+j)%(FIR_Dn/4)].d, FIR_Db[j].d, summed.d);
        // cout << "    ";
        // for (int k=0; k<4; k++)
        //     cout << FIR_Db[j].f[k] << "*" << FIR_Dz[(FIR_Di+j)%(FIR_Dn/4)].f[k] << ' ';
        // cout << endl;
    }
    summed.d = _mm256_hadd_ps(summed.d, summed.d); // out.I = summed.f[3];
    summed.d = _mm256_hadd_ps(summed.d, summed.d); // out.Q = summed.f[4];
    // cout << "]) = " << summed.f[3] << endl << endl;

    FIR_Di = (FIR_Di+FIR_Dn/4-1) % (FIR_Dn/4);

    return *(sample_t*)(summed.f+3);
#else
    sample_t summed;
    summed.I = 0.0f;
    summed.Q = 0.0f;

    FIR_Dz[FIR_Di] = in;

    // cout << "sum([ ";
    for (int j=0; j<FIR_Dn; j++) {
        summed.I += FIR_Db[j]*FIR_Dz[(FIR_Di+j)%FIR_Dn].I;
        summed.Q += FIR_Db[j]*FIR_Dz[(FIR_Di+j)%FIR_Dn].Q;
        // cout << FIR_Db[j] << "*" << FIR_Dz[(FIR_Di+j)%FIR_Dn].Q << ' ';
    }
    // cout << "]) = " << summed.Q << endl;
    FIR_Di = (FIR_Di+FIR_Dn-1) % FIR_Dn;

    return summed;
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Fast approximate exponent
///////////////////////////////////////////////////////////////////////////////////////////////
float fast_exp(float v) {
    union {
        float f;
        int32_t i;
    } tmp;

    tmp.i = (int32_t)(v*12102203.0f) + (int32_t)1065353216;
    return tmp.f;
}
