#pragma once
#include "test_case_template.h"
#include <cmath>

template <>
struct TestCaseImpl<TestCase::BOOK_BRATU3> {                // discrete Bratu problem, x=(lambda,u_1..u_3)
    static constexpr size_t NDIM = 4;
    static constexpr size_t KDIM = 3;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-0.13100000000000001, -0.11700000000000001, -0.12000000000000001, -0.12300000000000001};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{3.613, 4.5999999999999996, 4.6069999999999993, 6.5140000000000002};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{24, 20, 20, 30};
    static constexpr std::array<double, NDIM> FIRST_POINT{1, 0.10544866558254885, 0.14144676492828928, 0.10544866558254885};
    static void func(const double *x, double *f) {
        const double h=1.0/4.0;
        double u[5]; u[0]=0; u[4]=0;
        for(int i=1;i<=3;++i) u[i]=x[i];
        for(int i=1;i<=3;++i) f[i-1]=(u[i-1]-2*u[i]+u[i+1])/(h*h)+x[0]*std::exp(u[i]);
    }
};

template <>
struct TestCaseImpl<TestCase::BOOK_BRATU5> {                // discrete Bratu problem, x=(lambda,u_1..u_5)
    static constexpr size_t NDIM = 6;
    static constexpr size_t KDIM = 5;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-0.13100000000000001, -0.11700000000000001, -0.12000000000000001, -0.12300000000000001, -0.126, -0.129};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{3.613, 3.3999999999999999, 5.0069999999999997, 6.5140000000000002, 5.0209999999999999, 3.4279999999999999};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{12, 10, 10, 12, 10, 10};
    static constexpr std::array<double, NDIM> FIRST_POINT{1, 0.077481307750760883, 0.12494700344221833, 0.14093802130893604, 0.12494700344221833, 0.077481307750760883};
    static void func(const double *x, double *f) {
        const double h=1.0/6.0;
        double u[7]; u[0]=0; u[6]=0;
        for(int i=1;i<=5;++i) u[i]=x[i];
        for(int i=1;i<=5;++i) f[i-1]=(u[i-1]-2*u[i]+u[i+1])/(h*h)+x[0]*std::exp(u[i]);
    }
};

template <>
struct TestCaseImpl<TestCase::BOOK_BLOBS3> {                 // 14 Gaussian blobs, level 1/2
    static constexpr size_t NDIM = 3;
    static constexpr size_t KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-2.0299999999999998, -2.0369999999999999, -2.0439999999999996};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{2.0099999999999998, 2.0149999999999997, 2.0199999999999996};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{40, 40, 40};
    static constexpr std::array<double, NDIM> FIRST_POINT{-0.55881730105692062, 1.4296094171501585, -0.99417476947276384};
    static void func(const double *x, double *f) {
        static constexpr double C[14][3] = {{-0.7973617340529725, 1.4296094171501585, -0.99417476947276384}, {-1.0262674866620558, -0.4803544300929361, -0.8622680109123011}, {0.54542637687291107, -1.2317459772049681, 1.2681899958549776}, {1.1460175650685085, -1.5909534970028159, 0.13269171750014164}, {-1.2580759231240322, -0.77454413196483118, -0.26593266997407139}, {-0.14842841006951524, -0.10193090897951795, 1.3680534427914441}, {-0.77193251384911865, -0.99875132548947865, 0.54563351697440066}, {1.4291798507263289, 1.3529948014200044, 1.2167999997551155}, {-1.3940617730444327, 1.3974275997732279, 0.47756918097432122}, {1.1889787056809484, -0.29408490222085026, -0.89795220957331778}, {0.93750424606301808, 0.51723030717370744, 0.89228829447063118}, {-0.95569696671433668, -1.170074440645144, 0.84360028684054855}, {-1.5352680841277231, 1.4259221801555393, -1.1677188399703504}, {0.32035291369472318, -0.25917744636330187, -0.56353583317095657}};
        double s=0;
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[0][j]; d+=t*t; } s+=std::exp(-d/(0.080683978875612808)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[1][j]; d+=t*t; } s+=std::exp(-d/(0.1648981727914734)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[2][j]; d+=t*t; } s+=std::exp(-d/(0.18740666755335969)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[3][j]; d+=t*t; } s+=std::exp(-d/(0.15663727783797129)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[4][j]; d+=t*t; } s+=std::exp(-d/(0.13694215180408226)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[5][j]; d+=t*t; } s+=std::exp(-d/(0.15389205141237255)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[6][j]; d+=t*t; } s+=std::exp(-d/(0.12760450041235866)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[7][j]; d+=t*t; } s+=std::exp(-d/(0.13079294185463161)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[8][j]; d+=t*t; } s+=std::exp(-d/(0.18594790620476695)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[9][j]; d+=t*t; } s+=std::exp(-d/(0.0939594651213395)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[10][j]; d+=t*t; } s+=std::exp(-d/(0.086696550492838115)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[11][j]; d+=t*t; } s+=std::exp(-d/(0.19306677170328798)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[12][j]; d+=t*t; } s+=std::exp(-d/(0.19239043027344713)); }
        { double d=0; for(int j=0;j<3;++j){ double t=x[j]-C[13][j]; d+=t*t; } s+=std::exp(-d/(0.11916264333964334)); }
        f[0]=0.5-s;
    }
};

template <>
struct TestCaseImpl<TestCase::BOOK_BLOBS4> {                 // 12 Gaussian blobs, level 1/2
    static constexpr size_t NDIM = 4;
    static constexpr size_t KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-2.0299999999999998, -2.0369999999999999, -2.0439999999999996, -2.0509999999999997};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{2.0099999999999998, 2.0149999999999997, 2.0199999999999996, 2.0249999999999999};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{24, 24, 24, 24};
    static constexpr std::array<double, NDIM> FIRST_POINT{1.1956410550289696, 0.77839992107648204, 1.4376515877023048, -1.3385494743626301};
    static void func(const double *x, double *f) {
        static constexpr double C[12][4] = {{0.96248960354593338, 0.77839992107648204, 1.4376515877023048, -1.3385494743626301}, {1.2741030320277731, 0.00039631253746241768, -0.16331717896571396, 0.59763236489395277}, {0.37272583020866845, -0.20308376974026099, -0.66845712466880591, 1.3397245639007718}, {1.02019378155188, -1.2721165392794274, -0.30568692094731742, 0.84288692409670585}, {1.2129899087306062, 1.4791341726500202, -0.82969315413048506, 1.2535579204320917}, {0.62358580380245687, -1.3817498177333314, -0.92598770349921189, -0.99618564394443476}, {-1.3523503942790935, 0.59971192137036144, -0.60754982927490786, 0.61367444579928776}, {-1.3385327601861114, 1.1878716108847533, 0.62583405984610652, 0.90292303730081613}, {0.34063334411044588, -0.13527803251877368, -0.98234499116830598, 1.4155610115651482}, {-1.1345809617675078, 0.068293141251335232, -1.211414017199864, -1.2531398531474527}, {0.62961662531331486, 1.2441259752916136, -0.12429028249006424, 0.9426924451598957}, {1.1506341376919917, 0.11337855340909142, 0.88158709054466566, -0.66257163972345556}};
        double s=0;
        { double d=0; for(int j=0;j<4;++j){ double t=x[j]-C[0][j]; d+=t*t; } s+=std::exp(-d/(0.078414520366283894)); }
        { double d=0; for(int j=0;j<4;++j){ double t=x[j]-C[1][j]; d+=t*t; } s+=std::exp(-d/(0.082840957704068782)); }
        { double d=0; for(int j=0;j<4;++j){ double t=x[j]-C[2][j]; d+=t*t; } s+=std::exp(-d/(0.06320999943337724)); }
        { double d=0; for(int j=0;j<4;++j){ double t=x[j]-C[3][j]; d+=t*t; } s+=std::exp(-d/(0.069823647464341332)); }
        { double d=0; for(int j=0;j<4;++j){ double t=x[j]-C[4][j]; d+=t*t; } s+=std::exp(-d/(0.16018225931235974)); }
        { double d=0; for(int j=0;j<4;++j){ double t=x[j]-C[5][j]; d+=t*t; } s+=std::exp(-d/(0.18547177359062864)); }
        { double d=0; for(int j=0;j<4;++j){ double t=x[j]-C[6][j]; d+=t*t; } s+=std::exp(-d/(0.068356579345542276)); }
        { double d=0; for(int j=0;j<4;++j){ double t=x[j]-C[7][j]; d+=t*t; } s+=std::exp(-d/(0.15538242294608051)); }
        { double d=0; for(int j=0;j<4;++j){ double t=x[j]-C[8][j]; d+=t*t; } s+=std::exp(-d/(0.16445622308743685)); }
        { double d=0; for(int j=0;j<4;++j){ double t=x[j]-C[9][j]; d+=t*t; } s+=std::exp(-d/(0.14918555751656359)); }
        { double d=0; for(int j=0;j<4;++j){ double t=x[j]-C[10][j]; d+=t*t; } s+=std::exp(-d/(0.12236053411517558)); }
        { double d=0; for(int j=0;j<4;++j){ double t=x[j]-C[11][j]; d+=t*t; } s+=std::exp(-d/(0.12048091894518118)); }
        f[0]=0.5-s;
    }
};
