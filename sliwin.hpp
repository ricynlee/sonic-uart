#pragma once

#include <cstddef>

class sliwin {
protected:
    float* _data;
    size_t _oldest;
    size_t _capacity;
    size_t _size;
public:
    sliwin(size_t);
    ~sliwin();
public:
    void slide(const float);
    float& operator[](size_t);
    bool filled(void)const;
    size_t size(void)const;
};

class sliwin_sum: public sliwin {
protected:
    float _sum;
public:
    sliwin_sum(size_t);
    float slide(const float);
    float sum(void) const;
};

class sliwin_mean: public sliwin_sum {
protected:
    float _mean;
public:
    sliwin_mean(size_t);
    float slide(const float);
    float mean(void) const;
};

class sliwin_stdd: public sliwin_mean {
protected:
    float _stdd; // std deviation
public:
    sliwin_stdd(size_t);
    float slide(const float);
    float stdd(void) const;
};
