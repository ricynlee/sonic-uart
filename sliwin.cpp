#include "sliwin.hpp"
#include <immintrin.h>
#include <cmath>

///////////////////////////////////////////////////////////////////////
// sliwin
///////////////////////////////////////////////////////////////////////
sliwin::sliwin(size_t n) {
    if (n==0)
        throw "sliwin _size 0";
    _oldest = 0;
    _capacity = n;
    _size = 0;
#if defined(__AVX__)
    n = (n+7)/8*8; // ceiled to 8-float-aligned
    __m256* ptr = (__m256*)_mm_malloc(n*sizeof(float), 32);
    for (size_t i=0; i<n/8; i++) {
        ptr[i] = _mm256_setzero_ps();
    }
    _data = (float*)ptr;
#else
    _data = new float[n];
#endif
}

sliwin::~sliwin() {
#if defined(__AVX__)
    _mm_free(_data);
#else
    delete[] _data;
#endif
}

void sliwin::slide(const float e) {
    _data[_oldest] = e;
    _oldest = (_oldest+1) % _capacity;
    if (_size<_capacity)
        _size++;
}

float& sliwin::operator[](size_t i) { // [0] is latest
    if (i>=_size)
        throw "out of range";
    if (_size<_capacity)
        return _data[_size-1-i];
    else
        return _data[(_oldest+_capacity-1-i) % _capacity];
}

bool sliwin::filled(void) const {
    return (_size>=_capacity);
}

size_t sliwin::size(void) const {
    return _size;
}

///////////////////////////////////////////////////////////////////////
// sliwin_sum
///////////////////////////////////////////////////////////////////////
sliwin_sum::sliwin_sum(size_t n):sliwin(n) {
    _sum = 0.0f;
}

float sliwin_sum::slide(const float e) {
    _sum += e;
    if (_size>=_capacity)
        _sum -= _data[_oldest];
    sliwin::slide(e);
    return _sum;
}

float sliwin_sum::sum(void) const {
    return _sum;
}

///////////////////////////////////////////////////////////////////////
// sliwin_mean
///////////////////////////////////////////////////////////////////////
sliwin_mean::sliwin_mean(size_t n):sliwin_sum(n) {
    // nothing else
}

float sliwin_mean::slide(const float e) {
    sliwin_sum::slide(e);
    _mean = sum() / size();
    return _mean;
}

float sliwin_mean::mean(void) const {
    return _mean;
}

///////////////////////////////////////////////////////////////////////
// sliwin_stdd
///////////////////////////////////////////////////////////////////////
sliwin_stdd::sliwin_stdd(size_t n):sliwin_mean(n) {
    // nothing else
}

float sliwin_stdd::slide(const float e) {
    sliwin_mean::slide(e);
    _stdd = 0.0f;
#if defined(__AVX__)
    const __m256* ptr = (const __m256*)_data;
    __m256 tmp;
    for (size_t i=0; i<(size()+7)/8; i++) {
        tmp = _mm256_set1_ps(mean());
        tmp = _mm256_sub_ps(ptr[i], tmp);
        tmp = _mm256_mul_ps(tmp, tmp);
        for (size_t j=0; j<8 && i*8+j<size(); j++) {
            _stdd += ((const float*)(&tmp))[j];
        }
    }
#else
    for (size_t i=0; i<size(); i++) {
        _stdd += powf(_data[i]-mean(), 2);
    }
#endif
    if (size()>1) {
        _stdd = sqrtf(_stdd/(size()-1));
    } else {
        _stdd = 0;
    }
    return _stdd;
}

float sliwin_stdd::stdd(void) const {
    return _stdd;
}
