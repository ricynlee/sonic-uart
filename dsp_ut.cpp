// g++ -std=c++11 dsp.cpp dsp_ut.cpp -o ut.exe

#include "dsp.hpp"
#include <cmath>
#include <iostream>

using namespace std;

#define TEST_CASE_SHARED_CODE \
    fir_filter fir; \
    fir.init(coefs, coefs_len); \
 \
    sample_t output; \
    for (unsigned i = 0; i<sizeof(input)/sizeof(input[0]); i++) { \
        output = fir.filter(input[i]); \
        if (expected_output[i].I==0 && fabs(output.I)>1e-4) { \
            return -1; \
        } else if (output.I/expected_output[i].I>(1+1e-4) || output.I/expected_output[i].I<(1-1e-4)) { \
            return -1; \
        } \
 \
        if (expected_output[i].Q==0 && fabs(output.Q)>1e-4) { \
            return -1; \
        } else if (output.Q/expected_output[i].Q>(1+1e-4) || output.Q/expected_output[i].Q<(1-1e-4)) { \
            return -1; \
        } \
    } \
 \
    return 0

int case1() {
    // 滤波器系数
    float coefs[] = {1.0f, 0.0f, 0.0f, 0.0f};
    int coefs_len = 4;

    // 测试输入
    sample_t input[] = {
        {1.0f, 0.0f},
        {0.0f, 1.0f},
        {-1.0f, 0.0f},
        {0.0f, -1.0f},
        {0.5f, 0.5f}
    };

    // 预期输出
    sample_t expected_output[] = {
        {1.0f, 0.0f},
        {0.0f, 1.0f},
        {-1.0f, 0.0f},
        {0.0f, -1.0f},
        {0.5f, 0.5f}
    };

    TEST_CASE_SHARED_CODE;
}

int case2() {
    float coefs[] = {0.5f, 0.5f, 0.0f, 0.0f};
    int coefs_len = 4;

    sample_t input[] = {
        {1.0f, 1.0f},
        {1.0f, 1.0f},
        {0.0f, 0.0f},
        {0.0f, 0.0f},
        {1.0f, -1.0f}
    };

    sample_t expected_output[] = {
        {0.5f, 0.5f},
        {1.0f, 1.0f},
        {0.5f, 0.5f},
        {0.0f, 0.0f},
        {0.5f, -0.5f}
    };

    TEST_CASE_SHARED_CODE;
}

int case3() {
    float coefs[] = {0.25f, 0.5f, 0.25f, 0.0f};
    int coefs_len = 4;

    sample_t input[] = {
        {1.0f, 1.0f},
        {1.0f, 1.0f},
        {1.0f, 1.0f},
        {-1.0f, -1.0f},
        {-1.0f, -1.0f},
        {-1.0f, -1.0f}
    };

    sample_t expected_output[] = {
        {0.25f, 0.25f},
        {0.75f, 0.75f},
        {1.0f, 1.0f},
        {0.5f, 0.5f},
        {-0.5f, -0.5f},
        {-1.0f, -1.0f}
    };

    TEST_CASE_SHARED_CODE;
}

int case4() {
    float coefs[] = {0.5f, -0.5f, 0.0f, 0.0f};
    int coefs_len = 4;

    sample_t input[] = {
        {1.0f, 1.0f},
        {1.0f, 1.0f},
        {0.0f, 0.0f},
        {0.0f, 0.0f},
        {-1.0f, -1.0f}
    };

    sample_t expected_output[] = {
        {0.5f, 0.5f},
        {0.0f, 0.0f},
        {-0.5f, -0.5f},
        {0.0f, 0.0f},
        {-0.5f, -0.5f}
    };

    TEST_CASE_SHARED_CODE;
}

int case5() {
    float coefs[] = {0.1f, 0.2f, 0.3f, 0.2f, 0.1f, 0.05f, 0.025f, 0.0125f};
    int coefs_len = 8;

    sample_t input[] = {
        {1.0f, 0.0f},
        {0.0f, 0.0f},
        {0.0f, 0.0f},
        {0.0f, 0.0f},
        {0.0f, 0.0f},
        {0.0f, 0.0f},
        {0.0f, 0.0f},
        {0.0f, 0.0f},
        {0.0f, 1.0f},
        {0.0f, 1.0f}
    };

    sample_t expected_output[] = {
        {0.1000f, 0.0000f},  // t=0
        {0.2000f, 0.0000f},  // t=1
        {0.3000f, 0.0000f},  // t=2
        {0.2000f, 0.0000f},  // t=3
        {0.1000f, 0.0000f},  // t=4
        {0.0500f, 0.0000f},  // t=5
        {0.0250f, 0.0000f},  // t=6
        {0.0125f, 0.0000f},  // t=7
        {0.0000f, 0.1000f},  // t=8
        {0.0000f, 0.3000f}   // t=9 (0.0125 + 0.025)
    };

    TEST_CASE_SHARED_CODE;
}

int case6() {
    float coefs[] = {0.25f, 0.25f, 0.25f, 0.25f};
    int coefs_len = 4;

    sample_t input[] = {
        {1.0f, 0.0f},
        {0.0f, 0.0f},
        {0.0f, 0.0f},
        {0.0f, 0.0f},
        {0.0f, 0.0f}
    };

    sample_t expected_output[] = {
        {0.25f, 0.0f},
        {0.25f, 0.0f},
        {0.25f, 0.0f},
        {0.25f, 0.0f},
        {0.0f, 0.0f}
    };

    TEST_CASE_SHARED_CODE;
}

int case7() {
    const float coefs[12] = {
        0.05f, 0.1f, 0.15f, 0.2f,      // 主瓣
        0.15f, 0.1f, 0.05f, 0.025f,     // 过渡带
        0.0125f, 0.00625f, 0.003125f, 0.0015625f  // 衰减尾
    };
    int coefs_len = 12;

    const sample_t input[26] = {
        {0.63f, 0.49f}, {0.81f, -0.22f}, {-0.75f, 0.31f}, {0.83f, -0.66f},
        {0.26f, 0.41f}, {-0.80f, -0.94f}, {-0.44f, -0.45f}, {0.09f, -0.91f},
        {0.92f, -0.81f}, {0.93f, 0.65f}, {-0.68f, 0.39f}, {0.94f, -0.37f},
        {0.91f, 0.90f}, {-0.03f, -0.93f}, {0.60f, -0.12f}, {-0.72f, -0.24f},
        {-0.16f, 0.53f}, {0.83f, 0.59f}, {0.58f, -0.63f}, {0.92f, -0.02f},
        {0.31f, -0.11f}, {-0.93f, 0.29f}, {0.70f, 0.42f}, {0.87f, 0.51f},
        {0.36f, -0.45f}, {0.52f, 0.36f}
    };

    const sample_t expected_output[26] = {
        {0.03150000f, 0.02450000f},   // t=0
        {0.1035000f, 0.03800000f},    // t=1
        {0.1380000f, 0.06700001f},    // t=2
        {0.2140000f, 0.06300001f},    // t=3
        {0.2400000f, 0.03050000f},   // t=4
        {0.1450000f, -0.02700000f},   // t=5
        {0.1030000f, -0.1380000f},    // t=6
        {-0.001750009f, -0.2162500f}, // t=7
        {-0.05837500f, -0.3753750f},  // t=8
        {0.006812509f, -0.3999375f},  // t=9
        {0.1004063f, -0.3724688f},   // t=10
        {0.2477031f, -0.2684844f},   // t=11
        {0.3358594f, -0.06037502f},  // t=12
        {0.3125468f, 0.01498437f},   // t=13
        {0.3823594f, 0.01100000f},   // t=14
        {0.3402812f, -0.004234400f}, // t=15
        {0.2331875f, -0.08643749f},  // t=16
        {0.1709687f, -0.03623439f},  // t=17
        {0.1003281f, -0.01126564f},  // t=18
        {0.1815937f, 0.05707812f},   // t=19
        {0.3103281f, 0.05342186f},   // t=20
        {0.3331875f, -0.005546877f},// t=21
        {0.3281250f, 0.009421861f}, // t=22
        {0.2655782f, 0.06049999f},  // t=23
        {0.2075781f, 0.1175469f},   // t=24
        {0.2925625f, 0.1590000f}    // t=25
    };

    TEST_CASE_SHARED_CODE;
}

#define TEST(n) \
    if (case##n()==0) { \
        cout << "CASE" #n " PASS" << endl; \
    } else { \
        cout << "CASE" #n " FAIL" << endl; \
    }

int main() {
#ifdef __AVX__
    cout << "AVX used" << endl;
#else
    cout << "AVX NOT used" << endl;
#endif
    TEST(1);
    TEST(2);
    TEST(3);
    TEST(4);
    TEST(5);
    TEST(6);
    TEST(7);

    return 0;
}