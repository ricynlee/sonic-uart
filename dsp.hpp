#pragma once
#include <immintrin.h>

// oscillator for frequency mixing
static const float SIN[8] = {0, 0.353553390593274, -0.5, 0.353553390593274, 0, -0.353553390593273, 0.5, -0.353553390593274}; // 18kHz
static const float COS[8] = {0.5, -0.353553390593274, 0, 0.353553390593274, -0.5, 0.353553390593275, 0, -0.353553390593274}; // 18kHz
// static const float SIN[8] = {0, 0.353553390593274, 0.5, 0.353553390593274, 0, -0.353553390593274, -0.5, -0.353553390593274}; // 6kHz
// static const float COS[8] = {0.5, 0.353553390593274, 0, -0.353553390593274, -0.5, -0.353553390593274, 0, 0.353553390593274}; // 6kHz
static const double PI = 3.1415926535897932384626433;

// fir lpf: kaiser win, fs=48k, fpass=750, fstop=1000, ripple=1db, attenuation=27db
static const float B[256] = {0.00120992911979556,0.00129564804956317,0.00136676570400596,0.00142174365464598,0.00145921355579048,0.00147800240665674,0.00147715606726706,0.00145595986396074,0.00141395791433752,0.00135096732992679,0.00126709148753434,0.00116272724699229,0.00103857077192515,0.000895617355126888,0.000735158042516559,0.000558772124350071,0.000368314736988395,0.000165901015861891,-4.61144918517675e-05,-0.000265162147115916,-0.000488485209643841,-0.000713171437382698,-0.000936187920160592,-0.00115441845264286,-0.00136470340657979,-0.00156388164032251,-0.00174883333966136,-0.00191652437206358,-0.00206404970958829,-0.00218867859803140,-0.00228789588436484,-0.00235944497399032,-0.00240136566571891,-0.00241203140467405,-0.00239018280990422,-0.00233495584689081,-0.00224590837024152,-0.00212303875014186,-0.00196680147200823,-0.00177811714820564,-0.00155837472993881,-0.00130942987743765,-0.00103359564673156,-0.000733627355657518,-0.000412701308960095,-7.43871496524662e-05,0.000277386046946049,0.000638369005173445,0.00100403709802777,0.00136964058037847,0.00173025834374130,0.00208085658960044,0.00241635111160576,0.00273166992701590,0.00302182114683092,0.00328195933252573,0.00350745348259807,0.00369395385496318,0.00383745739236474,0.00393437221646309,0.00398157583549619,0.00397647498175502,0.00391705287620425,0.00380191951990128,0.00363034987822175,0.00340231647714973,0.00311851664446294,0.00278038950636983,0.00239012623205781,0.00195067003369331,0.00146570859942585,0.000939657445997000,0.000377633259631693,-0.000214580635656603,-0.000830578152090311,-0.00146338273771107,-0.00210550706833601,-0.00274902139790356,-0.00338562787510455,-0.00400674063712359,-0.00460357265546918,-0.00516722677275538,-0.00568878883495927,-0.00615942478179932,-0.00657047890126705,-0.00691357301548123,-0.00718070333823562,-0.00736433872953057,-0.00745751056820154,-0.00745390541851521,-0.00734794698655605,-0.00713487528264523,-0.00681081600487232,-0.00637284573167563,-0.00581904314458370,-0.00514853606000543,-0.00436153262853622,-0.00345934880897403,-0.00244441628456116,-0.00132028630468994,-9.16171047720127e-05,0.00123584805987775,0.00265531521290541,0.00415898533537984,0.00573811167851090,0.00738306529819965,0.00908341445028782,0.0108280088752508,0.0126050785183907,0.0144023289903998,0.0162070598453283,0.0180062670260668,0.0197867769747973,0.0215353518724442,0.0232388228178024,0.0248842108994722,0.0264588501304388,0.0279505085200071,0.0293475054204464,0.0306388176977634,0.0318142026662827,0.0328642725944519,0.0337806157767773,0.0345558524131775,0.0351837202906609,0.0356591381132603,0.0359782464802265,0.0361384488642216,0.0361384488642216,0.0359782464802265,0.0356591381132603,0.0351837202906609,0.0345558524131775,0.0337806157767773,0.0328642725944519,0.0318142026662827,0.0306388176977634,0.0293475054204464,0.0279505085200071,0.0264588501304388,0.0248842108994722,0.0232388228178024,0.0215353518724442,0.0197867769747973,0.0180062670260668,0.0162070598453283,0.0144023289903998,0.0126050785183907,0.0108280088752508,0.00908341445028782,0.00738306529819965,0.00573811167851090,0.00415898533537984,0.00265531521290541,0.00123584805987775,-9.16171047720127e-05,-0.00132028630468994,-0.00244441628456116,-0.00345934880897403,-0.00436153262853622,-0.00514853606000543,-0.00581904314458370,-0.00637284573167563,-0.00681081600487232,-0.00713487528264523,-0.00734794698655605,-0.00745390541851521,-0.00745751056820154,-0.00736433872953057,-0.00718070333823562,-0.00691357301548123,-0.00657047890126705,-0.00615942478179932,-0.00568878883495927,-0.00516722677275538,-0.00460357265546918,-0.00400674063712359,-0.00338562787510455,-0.00274902139790356,-0.00210550706833601,-0.00146338273771107,-0.000830578152090311,-0.000214580635656603,0.000377633259631693,0.000939657445997000,0.00146570859942585,0.00195067003369331,0.00239012623205781,0.00278038950636983,0.00311851664446294,0.00340231647714973,0.00363034987822175,0.00380191951990128,0.00391705287620425,0.00397647498175502,0.00398157583549619,0.00393437221646309,0.00383745739236474,0.00369395385496318,0.00350745348259807,0.00328195933252573,0.00302182114683092,0.00273166992701590,0.00241635111160576,0.00208085658960044,0.00173025834374130,0.00136964058037847,0.00100403709802777,0.000638369005173445,0.000277386046946049,-7.43871496524662e-05,-0.000412701308960095,-0.000733627355657518,-0.00103359564673156,-0.00130942987743765,-0.00155837472993881,-0.00177811714820564,-0.00196680147200823,-0.00212303875014186,-0.00224590837024152,-0.00233495584689081,-0.00239018280990422,-0.00241203140467405,-0.00240136566571891,-0.00235944497399032,-0.00228789588436484,-0.00218867859803140,-0.00206404970958829,-0.00191652437206358,-0.00174883333966136,-0.00156388164032251,-0.00136470340657979,-0.00115441845264286,-0.000936187920160592,-0.000713171437382698,-0.000488485209643841,-0.000265162147115916,-4.61144918517675e-05,0.000165901015861891,0.000368314736988395,0.000558772124350071,0.000735158042516559,0.000895617355126888,0.00103857077192515,0.00116272724699229,0.00126709148753434,0.00135096732992679,0.00141395791433752,0.00145595986396074,0.00147715606726706,0.00147800240665674,0.00145921355579048,0.00142174365464598,0.00136676570400596,0.00129564804956317,0.00120992911979556};
static const int ORDER = sizeof(B)/sizeof(B[0])-1;

// maximum-length sequence for frequency spreading/despreading
static const char MSEQ[63] = {1,0,0,0,0,0,1,1,1,1,0,0,1,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,0,0,1,0,0,0,1,0,1,1,0,1,1,1,1,1,1,0,1,0,1,1,1,0,0,0,1,1,0,0,1,1,1,0,1};
static const int CHIPS = sizeof(MSEQ)/sizeof(MSEQ[0]);

// misc
static const int SAMPLE_RATE = 48000;

static const int CHIP_PREFIX = 1024;
static const int CHIP_BODY = 2048;
static const int CHIP_SUFFIX = CHIP_PREFIX;
static const int ACCUMUL_TIMES = 16;
static const int ACCUMUL_SAMPLES = (CHIP_PREFIX+CHIP_BODY+CHIP_SUFFIX)/ACCUMUL_TIMES;

static const int SYMBOL_CYCLIC_PREFIX = 1024;
static const int SYMBOL_BODY = 2048;
static const int SYMBOL_CYCLIC_SUFFIX = SYMBOL_CYCLIC_PREFIX;

static const int LENGTH_BITS = 7;

static const int TX_BUF_DEPTH = 512;    // common divisor of samples per chip and samples per symbol
                                        // >ORDER
                                        // cannot be too small (e.g., <256) in case of overflow/underflow

// typedefs
typedef struct {
    union {
        float L;
        float I;
    };
    union {
        float R;
        float Q;
    };
} sample_t;

// declarations
void init_filter(void);
sample_t filter(sample_t);

enum {
    PSK2 = 1,
    PSK4 = 2,
    QAM16 = 4
};

#define MODEM PSK2