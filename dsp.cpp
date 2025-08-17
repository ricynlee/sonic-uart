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
// Chirp preamble
///////////////////////////////////////////////////////////////////////////////////////////////
float chirp(size_t index) {
    const double u = 1500. * SAMPLE_RATE / PREAM_BODY / 2;
    double t = (double)index/SAMPLE_RATE;
    return (float)cos(2*PI*(-750*t+u*t*t));
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

///////////////////////////////////////////////////////////////////////////////////////////////
// Biquad filter (2nd order IIR filter)
///////////////////////////////////////////////////////////////////////////////////////////////
typedef struct {
    float c[5]; // coef
    sample_t z[3];
} biquad_filter_data_t;

#define IIR_Db (((biquad_filter_data_t*)data)->c)
#define IIR_Da (((biquad_filter_data_t*)data)->c + 2)
#define IIR_Dz (((biquad_filter_data_t*)data)->z)

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
    memcpy(IIR_Db, coef, 3*sizeof(float));
    memcpy(IIR_Da+1, coef+4, 2*sizeof(float));
    clear();
}

void biquad_filter::clear() {
    memset(IIR_Dz, 0, sizeof IIR_Dz);
}

sample_t biquad_filter::filter(const sample_t& in) { // direct form ii
    sample_t out;
    IIR_Dz[0].I = in.I - IIR_Da[1]*IIR_Dz[1].I - IIR_Da[2]*IIR_Dz[2].I;
    out.I       = IIR_Db[0]*IIR_Dz[0].I + IIR_Db[1]*IIR_Dz[1].I + IIR_Db[2]*IIR_Dz[2].I;
    IIR_Dz[0].Q = in.Q - IIR_Da[1]*IIR_Dz[1].Q - IIR_Da[2]*IIR_Dz[2].Q;
    out.Q       = IIR_Db[0]*IIR_Dz[0].Q + IIR_Db[1]*IIR_Dz[1].Q + IIR_Db[2]*IIR_Dz[2].Q;
    IIR_Dz[2]   = IIR_Dz[1];
    IIR_Dz[1]   = IIR_Dz[0];
    return out;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// OFDM modem (FFT)
///////////////////////////////////////////////////////////////////////////////////////////////
// static sample_t weight_fft(int butterfly_index, int weight_index) {
//     sample_t result;
//     result.I = cos(2 * M_PI * butterfly_index * weight_index / 512);
//     result.Q = -sin(2 * M_PI * butterfly_index * weight_index / 512);
//     return result;
// }
// 
// static sample_t weight_x2_ifft(int butterfly_index, int weight_index) {
//     sample_t result;
//     // result.I = cos(2 * M_PI * butterfly_index * weight_index / 512) / 8; // unscaled
//     // result.Q = sin(2 * M_PI * butterfly_index * weight_index / 512) / 8; // unscaled
//     result.I = cos(2 * M_PI * butterfly_index * weight_index / 512) / 4; // scaled, (8/4)^3=8 times larger for 512-pt ifft
//     result.Q = sin(2 * M_PI * butterfly_index * weight_index / 512) / 4; // scaled, (8/4)^3=8 times larger for 512-pt ifft
//     return result;
// }
// 
// static sample_t cmplx_add(sample_t x, sample_t y) {
//     sample_t result;
//     result.I = x.I + y.I;
//     result.Q = x.Q + y.Q;
//     return result;
// }
// 
// static sample_t cmplx_sub(sample_t x, sample_t y) {
//     sample_t result;
//     result.I = x.I - y.I;
//     result.Q = x.Q - y.Q;
//     return result;
// }
// 
// static sample_t cmplx_mult(sample_t x, sample_t y) {
//     sample_t result;
//     result.I = x.I * y.I - x.Q * y.Q;
//     result.Q = x.Q * y.I + x.I * y.Q;
//     return result;
// }
// 
// void fft(sample_t* const data /* 512-pt inout */) { // natural order in, bit-reversed order out
//     sample_t butterfly[3][8]; // each stage of radix-8 calc requires 3 rounds of butterflies
//     int data_index[8];
//     int weight_index;
// 
//     for(int i = 0; i < 3; i++){
//         int group_num, group_size;
//         group_num = 1 << (3*i);
//         group_size = 1 << (6-3*i);
//         for(int j = 0; j < group_num; j++){
//             for(int k = 0; k < group_size; k++){
//                 for(int m = 0; m < 8; m++){
//                     data_index[m] = j * group_size * 8 + k + m * group_size;
//                 }
//                 weight_index = k * group_num;
// 
//                 butterfly[0][0] = cmplx_add(data[data_index[0]], data[data_index[4]]);
//                 butterfly[0][4] = cmplx_sub(data[data_index[0]], data[data_index[4]]);
//                 butterfly[0][2] = cmplx_add(data[data_index[2]], data[data_index[6]]);
//                 butterfly[0][6] = cmplx_sub(data[data_index[2]], data[data_index[6]]);
//                 butterfly[0][1] = cmplx_add(data[data_index[1]], data[data_index[5]]);
//                 butterfly[0][5] = cmplx_sub(data[data_index[1]], data[data_index[5]]);
//                 butterfly[0][3] = cmplx_add(data[data_index[3]], data[data_index[7]]);
//                 butterfly[0][7] = cmplx_sub(data[data_index[3]], data[data_index[7]]);
//                 butterfly[0][6] = cmplx_mult(butterfly[0][6], (sample_t){0,-1});
//                 butterfly[0][7] = cmplx_mult(butterfly[0][7], (sample_t){0,-1});
//                 butterfly[1][0] = cmplx_add(butterfly[0][0], butterfly[0][2]);
//                 butterfly[1][2] = cmplx_sub(butterfly[0][0], butterfly[0][2]);
//                 butterfly[1][4] = cmplx_add(butterfly[0][4], butterfly[0][6]);
//                 butterfly[1][6] = cmplx_sub(butterfly[0][4], butterfly[0][6]);
//                 butterfly[1][1] = cmplx_add(butterfly[0][1], butterfly[0][3]);
//                 butterfly[1][3] = cmplx_sub(butterfly[0][1], butterfly[0][3]);
//                 butterfly[1][5] = cmplx_add(butterfly[0][5], butterfly[0][7]);
//                 butterfly[1][7] = cmplx_sub(butterfly[0][5], butterfly[0][7]);
//                 butterfly[1][5] = cmplx_mult(butterfly[1][5], (sample_t){cos(M_PI/4),-sin(M_PI/4)});
//                 butterfly[1][3] = cmplx_mult(butterfly[1][3], (sample_t){0,-1});
//                 butterfly[1][7] = cmplx_mult(butterfly[1][7], (sample_t){cos(3*M_PI/4),-sin(3*M_PI/4)});
//                 butterfly[2][0] = cmplx_add(butterfly[1][0], butterfly[1][1]);
//                 butterfly[2][1] = cmplx_add(butterfly[1][4], butterfly[1][5]);
//                 butterfly[2][2] = cmplx_add(butterfly[1][2], butterfly[1][3]);
//                 butterfly[2][3] = cmplx_add(butterfly[1][6], butterfly[1][7]);
//                 butterfly[2][4] = cmplx_sub(butterfly[1][0], butterfly[1][1]);
//                 butterfly[2][5] = cmplx_sub(butterfly[1][4], butterfly[1][5]);
//                 butterfly[2][6] = cmplx_sub(butterfly[1][2], butterfly[1][3]);
//                 butterfly[2][7] = cmplx_sub(butterfly[1][6], butterfly[1][7]);
//                 data[data_index[0]] = cmplx_mult(butterfly[2][0], weight_fft(0,weight_index));
//                 data[data_index[1]] = cmplx_mult(butterfly[2][1], weight_fft(1,weight_index));
//                 data[data_index[2]] = cmplx_mult(butterfly[2][2], weight_fft(2,weight_index));
//                 data[data_index[3]] = cmplx_mult(butterfly[2][3], weight_fft(3,weight_index));
//                 data[data_index[4]] = cmplx_mult(butterfly[2][4], weight_fft(4,weight_index));
//                 data[data_index[5]] = cmplx_mult(butterfly[2][5], weight_fft(5,weight_index));
//                 data[data_index[6]] = cmplx_mult(butterfly[2][6], weight_fft(6,weight_index));
//                 data[data_index[7]] = cmplx_mult(butterfly[2][7], weight_fft(7,weight_index));
//             }
//         }
//     }
// }
// 
// void ifft_x8(sample_t* const data /* 512-pt inout */) { // bit-reversed order in, natural order out, output is x8 scaled
//     sample_t butterfly[3][8];
//     int data_index[8];
//     int weight_index;
// 
//     for(int i = 0; i < 3; i++){
//         int group_num, group_size;
//         group_num = 1 << (6-3*i);
//         group_size = 1 << (3*i);
//         for(int j = 0; j < group_num; j++){
//             for(int k = 0; k < group_size; k++){
//                 for(int m = 0; m < 8; m++){
//                     data_index[m] = j * group_size * 8 + k + m * group_size;
//                 }
//                 weight_index = k * group_num;
// 
//                 butterfly[0][0] = cmplx_mult(data[data_index[0]], weight_x2_ifft(0,weight_index));
//                 butterfly[0][1] = cmplx_mult(data[data_index[1]], weight_x2_ifft(1,weight_index));
//                 butterfly[0][2] = cmplx_mult(data[data_index[2]], weight_x2_ifft(2,weight_index));
//                 butterfly[0][3] = cmplx_mult(data[data_index[3]], weight_x2_ifft(3,weight_index));
//                 butterfly[0][4] = cmplx_mult(data[data_index[4]], weight_x2_ifft(4,weight_index));
//                 butterfly[0][5] = cmplx_mult(data[data_index[5]], weight_x2_ifft(5,weight_index));
//                 butterfly[0][6] = cmplx_mult(data[data_index[6]], weight_x2_ifft(6,weight_index));
//                 butterfly[0][7] = cmplx_mult(data[data_index[7]], weight_x2_ifft(7,weight_index));
//                 butterfly[1][0] = cmplx_add(butterfly[0][0], butterfly[0][4]);
//                 butterfly[1][4] = cmplx_sub(butterfly[0][0], butterfly[0][4]);
//                 butterfly[1][2] = cmplx_add(butterfly[0][2], butterfly[0][6]);
//                 butterfly[1][6] = cmplx_sub(butterfly[0][2], butterfly[0][6]);
//                 butterfly[1][1] = cmplx_add(butterfly[0][1], butterfly[0][5]);
//                 butterfly[1][5] = cmplx_sub(butterfly[0][1], butterfly[0][5]);
//                 butterfly[1][3] = cmplx_add(butterfly[0][3], butterfly[0][7]);
//                 butterfly[1][7] = cmplx_sub(butterfly[0][3], butterfly[0][7]);
//                 butterfly[1][6] = cmplx_mult(butterfly[1][6], (sample_t){0,1});
//                 butterfly[1][7] = cmplx_mult(butterfly[1][7], (sample_t){0,1});
//                 butterfly[2][0] = cmplx_add(butterfly[1][0], butterfly[1][2]);
//                 butterfly[2][2] = cmplx_sub(butterfly[1][0], butterfly[1][2]);
//                 butterfly[2][4] = cmplx_add(butterfly[1][4], butterfly[1][6]);
//                 butterfly[2][6] = cmplx_sub(butterfly[1][4], butterfly[1][6]);
//                 butterfly[2][1] = cmplx_add(butterfly[1][1], butterfly[1][3]);
//                 butterfly[2][3] = cmplx_sub(butterfly[1][1], butterfly[1][3]);
//                 butterfly[2][5] = cmplx_add(butterfly[1][5], butterfly[1][7]);
//                 butterfly[2][7] = cmplx_sub(butterfly[1][5], butterfly[1][7]);
//                 butterfly[2][5] = cmplx_mult(butterfly[2][5], (sample_t){cos(M_PI/4),sin(M_PI/4)});
//                 butterfly[2][3] = cmplx_mult(butterfly[2][3], (sample_t){0,1});
//                 butterfly[2][7] = cmplx_mult(butterfly[2][7], (sample_t){cos(3*M_PI/4),sin(3*M_PI/4)});
//                 data[data_index[0]] = cmplx_add(butterfly[2][0], butterfly[2][1]);
//                 data[data_index[1]] = cmplx_add(butterfly[2][4], butterfly[2][5]);
//                 data[data_index[2]] = cmplx_add(butterfly[2][2], butterfly[2][3]);
//                 data[data_index[3]] = cmplx_add(butterfly[2][6], butterfly[2][7]);
//                 data[data_index[4]] = cmplx_sub(butterfly[2][0], butterfly[2][1]);
//                 data[data_index[5]] = cmplx_sub(butterfly[2][4], butterfly[2][5]);
//                 data[data_index[6]] = cmplx_sub(butterfly[2][2], butterfly[2][3]);
//                 data[data_index[7]] = cmplx_sub(butterfly[2][6], butterfly[2][7]);
//             }
//         }
//     }
// }
// 
// static uint16_t reverse_bit(uint16_t n) { // radix-8, 512-pt only
//     n = ((n & 0b111000000) >> 6) | ((n & 0b000000111) << 6) | (n & 0b000111000);
//     return n;
// }
// 
// void reorder(sample_t* const data /* 512-pt inout */) { // in-place
//     sample_t tmp;
//     for (uint16_t i=0; i<512; i++) {
//         uint16_t j =  reverse_bit(i);
//         if (j>i) {
//             tmp = data[i];
//             data[i] = data[j];
//             data[j] = tmp;
//         }
//     }
// }
// 
// void ofdm_modulate(const sample_t* const constel /* 16-pt in */, sample_t* const data /* 512-pt out */) {
//     uint16_t j;
//     memset((void*)data, 0, sizeof(sample_t)*512);
// 
//     for(uint16_t i=0; i<8; i++) {
//         j = reverse_bit(i+1);
//         data[j].I = constel[i].I*4;
//         data[j].Q = constel[i].Q*4;
//     }
// 
//     for(uint16_t i=8; i<16; i++) {
//         j = reverse_bit(i+504-8);
//         data[j].I = constel[i].I*4;
//         data[j].Q = constel[i].Q*4;
//     }
// 
//     ifft_x8(data); // 16-pt data, constel x4 and ifft x8, so max abs val of time-domain signal is 1
// }
// 
// void ofdm_demodulate(sample_t* const data /* 512-pt in */, sample_t* const constel /* 16-pt out */) {
//     uint16_t j;
//     fft(data);
// 
//     for(uint16_t i=0; i<8; i++) {
//         j = reverse_bit(i+1);
//         constel[i] = data[j];
//     }
// 
//     for(uint16_t i=8; i<16; i++) {
//         j = reverse_bit(i+504-8);
//         constel[i] = data[j];
//     }
// }
// 