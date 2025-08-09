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

#define Dn (((fir_filter_data_t*)data)->n)
#define Di (((fir_filter_data_t*)data)->i)
#define Db (((fir_filter_data_t*)data)->b)
#define Dz (((fir_filter_data_t*)data)->z)

fir_filter::fir_filter() {
    data = NULL;
}

void fir_filter::init(const float* const coef, int len /* order+1 */) {
    // coef length = 4n
    // order = 4n-1
    if (data)
        return;
    data = malloc(sizeof(fir_filter_data_t));
    Dn = len;
    Di = 0;
#if defined(__AVX__)
    packed_t vindex;
    for (int i=0; i<8; i++) {
        vindex.i[i] = (i%4)*len/4;
    }

    Db = (packed_t*)_mm_malloc(Dn*4*sizeof(float), 32);
    Dz = Db + Dn*2*sizeof(float)/sizeof(packed_t);

    for (int i=0; i<Dn/4; i++) {
        Db[i].d = _mm256_i32gather_ps(coef+i, vindex.di, 4); // reshape coef
        Dz[i].d = _mm256_setzero_ps();
    }

    // for (int i=0; i<Dn/4; i++) {
    //     for (int j=0; j<8; j++)
    //         cout << Db[i].f[j] << ' ';
    //     cout << endl;
    // }
#else
    Db = (float*)malloc(Dn*sizeof(float) + Dn*sizeof(sample_t));
    Dz = (sample_t*)(Db + Dn);
    memset(Dz, 0, len*sizeof(sample_t));
    memcpy(Db, coef, len*sizeof(float));
    // for (int i=0; i<Dn; i++)
    //     cout << Db[i] << ' ';
    // cout << endl;
#endif
}

void fir_filter::clear() {
#if defined(__AVX__)
    for (int i=0; i<Dn/4; i++) {
        Dz[i].d = _mm256_setzero_ps();
    }
#else
    memset(Dz, 0, Dn*sizeof(float));
#endif
}

fir_filter::~fir_filter() {
#if defined(__AVX__)
    if (data) {
        _mm_free(Db);
        free(data);
    }
#else
    if (data) {
        free(Db);
        free(data);
    }
#endif
}

sample_t fir_filter::filter(const sample_t& in) {
#if defined(__AVX__)
    static packed_t summed;
    summed.d = _mm256_setzero_ps();

    Dz[Di].d = _mm256_permute_ps(Dz[Di].d, 0x93);
    Dz[Di].f[0] = in.I;
    Dz[Di].f[4] = in.Q;

    // cout << "sum([" << endl;
    for (int j=0; j<Dn/4; j++) {
        summed.d = _mm256_fmadd_ps(Dz[(Di+j)%(Dn/4)].d, Db[j].d, summed.d);
        // cout << "    ";
        // for (int k=0; k<4; k++)
        //     cout << Db[j].f[k] << "*" << Dz[(Di+j)%(Dn/4)].f[k] << ' ';
        // cout << endl;
    }
    summed.d = _mm256_hadd_ps(summed.d, summed.d); // out.I = summed.f[3];
    summed.d = _mm256_hadd_ps(summed.d, summed.d); // out.Q = summed.f[4];
    // cout << "]) = " << summed.f[3] << endl << endl;

    Di = (Di+Dn/4-1) % (Dn/4);

    return *(sample_t*)(summed.f+3);
#else
    sample_t summed;
    summed.I = 0.0f;
    summed.Q = 0.0f;

    Dz[Di] = in;

    // cout << "sum([ ";
    for (int j=0; j<Dn; j++) {
        summed.I += Db[j]*Dz[(Di+j)%Dn].I;
        summed.Q += Db[j]*Dz[(Di+j)%Dn].Q;
        // cout << Db[j] << "*" << Dz[(Di+j)%Dn].Q << ' ';
    }
    // cout << "]) = " << summed.Q << endl;
    Di = (Di+Dn-1) % Dn;

    return summed;
#endif
}

// preamble generator
float chirp(size_t index) {
    const double u = 1280. * SAMPLE_RATE / PREAM_BODY / 2;
    double t = (double)index/SAMPLE_RATE;
    return (float)cos(2*PI*(-640*t+u*t*t));
}

// fast approximate exp
float fast_exp(float v) {
    union {
        float f;
        int32_t i;
    } tmp;

    tmp.i = (int32_t)(v*12102203.0f) + (int32_t)1065353216;
    return tmp.f;
}

typedef struct {
    float c[5]; // coef
    sample_t z[3];
} biquad_filter_data_t;

#undef Db
#undef Da
#undef Dz
#define Db (((biquad_filter_data_t*)data)->c)
#define Da (((biquad_filter_data_t*)data)->c + 2)
#define Dz (((biquad_filter_data_t*)data)->z)

biquad_filter::biquad_filter() {
    data = NULL;
}

biquad_filter::~biquad_filter() {
    if (data) {
        free(data);
    }
}

void biquad_filter::init(const float* const coef) {
    data = (biquad_filter_data_t*) malloc(sizeof(biquad_filter_data_t));
    memcpy(Db, coef, 3*sizeof(float));
    memcpy(Da+1, coef+4, 2*sizeof(float));
    clear();
}

void biquad_filter::clear() {
    memset(Dz, 0, sizeof Dz);
}

sample_t biquad_filter::filter(const sample_t& in) { // direct form ii
    sample_t out;
    Dz[0].I =   in.I -          Da[1]*Dz[1].I - Da[2]*Dz[2].I;
    out.I =     Db[0]*Dz[0].I + Db[1]*Dz[1].I + Db[2]*Dz[2].I;
    Dz[0].Q =   in.Q -          Da[1]*Dz[1].Q - Da[2]*Dz[2].Q;
    out.Q =     Db[0]*Dz[0].Q + Db[1]*Dz[1].Q + Db[2]*Dz[2].Q;
    Dz[2] = Dz[1];
    Dz[1] = Dz[0];
    return out;
}
