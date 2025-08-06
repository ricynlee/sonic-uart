#include "sliwin.hpp"

///////////////////////////////////////////////////////////////////////
// sliwin
///////////////////////////////////////////////////////////////////////
sliwin::sliwin(size_t n) {
    if (n==0)
        throw "sliwin _size 0";
    _data = new float[n];
    _oldest = 0;
    _capacity = n;
    _size = 0;
}

sliwin::~sliwin() {
    delete[] _data;
}

void sliwin::slide(const float e) {
    _data[_oldest] = e;
    _oldest = (_oldest+1) % _capacity;
    if (_size<_capacity)
        _size++;
}

float& sliwin::operator[](size_t i) {
    if (i>=_size)
        throw "out of range";
    if (_size<_capacity)
        return _data[i];
    else
        return _data[(_oldest+i) % _size];
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
