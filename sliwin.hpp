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
private:
    float _sum;
public:
    sliwin_sum(size_t);
    float slide(const float);
    float sum(void) const;
};
