#include "sliwin.hpp"
#include <iostream>
#include <random>

using namespace std;

class sliwin_peak_finder:public sliwin {
private:
    bool _found;
    float _peak;
public:
    sliwin_peak_finder(size_t);
    bool slide(float);
    float peak(void);
};

sliwin_peak_finder::sliwin_peak_finder(size_t n):sliwin(n) {
    _found = false;
}

bool sliwin_peak_finder::slide(float v) {
    _found = false;
    sliwin::slide(v);
    if (size()>=3) {
        v = operator[](1); // hope that [1] is the peak
        if ((v>operator[](0) && v>=operator[](2)) || (v>=operator[](0) && v>operator[](2))) { // neighborhood peak found
            for (size_t i=3; i<size(); i++) {
                if (operator[](i) > v) {
                    return false;
                }
            }
            _peak = v;
            _found = true; // regional peak found
        }
    }
    return _found;
}

float sliwin_peak_finder::peak(void) {
    if (_found)
        return _peak;
    else
        return 0.0f;
}

void test_sliwin_peak_finder() {
    uniform_real_distribution<float> dis(0.0, 10.0);
    minstd_rand gen;

    sliwin_peak_finder sw(32);
    for (int i=0; i<64; i++) {
        bool found = sw.slide(dis(gen));
        for (int j=0; j<sw.size(); j++) {
            cout << sw[j] << ' ';
        }
        if (found) {
            cout << "found " << sw.peak() << endl;
        } else {
            cout << endl;
        }
    }
}

void test_sliwin_stdd() {
    sliwin_stdd sw(32);
    for (int i=0; i<64; i++) {
        sw.slide((float)i);
        for (int j=0; j<sw.size(); j++) {
            cout << sw[j] << ' ';
        }
        cout << "stdd " << sw.stdd() << endl;
    }
}

int main() {
    test_sliwin_stdd();
    return 0;
}