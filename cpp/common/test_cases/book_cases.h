#pragma once
#include "test_case_template.h"
#include <cmath>

template <>
struct TestCaseImpl<TestCase::BOOK_SPHERE3> {
    static constexpr size_t NDIM = 3;
    static constexpr size_t KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.1, -1.1, -1.1};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.1, 1.1, 1.1};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{12, 12, 12};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.6, 0.48, 0.64};
    static void func(const double *x, double *f) {
        f[0] = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]-1.0;
    }
};

template <>
struct TestCaseImpl<TestCase::BOOK_SPHERE4> {
    static constexpr size_t NDIM = 4;
    static constexpr size_t KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.1, -1.1, -1.1, -1.1};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.1, 1.1, 1.1, 1.1};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{8, 8, 8, 8};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.5, 0.4, 0.3, 0.7071067811865476};
    static void func(const double *x, double *f) {
        f[0] = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]-1.0;
    }
};

template <>
struct TestCaseImpl<TestCase::BOOK_TORUS4> {
    static constexpr size_t NDIM = 4;
    static constexpr size_t KDIM = 2;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.1, -1.1, -0.51, -0.51};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.1, 1.1, 0.51, 0.51};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{10, 10, 10, 10};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.955336489125606, 0.29552020666133955, 0.38242109364224425, 0.3221088436188455};
    static void func(const double *x, double *f) {
        f[0] = x[0]*x[0]+x[1]*x[1]-1.0;
        f[1] = x[2]*x[2]+x[3]*x[3]-0.25;
    }
};

template <>
struct TestCaseImpl<TestCase::BOOK_ZCOSW> {
    static constexpr size_t NDIM = 4;
    static constexpr size_t KDIM = 2;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-3.17, -3.16, -3.18, -3.15};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{3.11, 3.12, 3.10, 3.13};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{10, 10, 10, 10};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.3, 0.2, 0.9745069929868755, -0.05949885707931209};
    static void func(const double *x, double *f) {
        f[0] = x[2]-std::cos(x[0])*std::cosh(x[1]);
        f[1] = x[3]+std::sin(x[0])*std::sinh(x[1]);
    }
};

template <>
struct TestCaseImpl<TestCase::BOOK_Z2COSW> {
    static constexpr size_t NDIM = 4;
    static constexpr size_t KDIM = 2;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-3.17, -3.16, -3.18, -3.15};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{3.11, 3.12, 3.10, 3.13};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{10, 10, 10, 10};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.3, 0.2, 0.9876306642408408, -0.03012201789271443};
    static void func(const double *x, double *f) {
        f[0] = x[2]*x[2]-x[3]*x[3]-std::cos(x[0])*std::cosh(x[1]);
        f[1] = 2.0*x[2]*x[3]+std::sin(x[0])*std::sinh(x[1]);
    }
};

template <>
struct TestCaseImpl<TestCase::BOOK_Z2W5> {
    static constexpr size_t NDIM = 4;
    static constexpr size_t KDIM = 2;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-4.03, -4.02, -4.04, -4.01};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{3.97, 3.98, 3.96, 3.99};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{12, 12, 12, 12};
    static constexpr std::array<double, NDIM> FIRST_POINT{1.096166836646079, 0.45557846060000695, 0.3, -0.2};
    static void func(const double *x, double *f) {
        const double u=x[2], v=x[3];
        const double re = u*u*u*u*u - 10.0*u*u*u*v*v + 5.0*u*v*v*v*v;
        const double im = 5.0*u*u*u*u*v - 10.0*u*u*v*v*v + v*v*v*v*v;
        f[0] = x[0]*x[0]-x[1]*x[1]-re-1.0;
        f[1] = 2.0*x[0]*x[1]-im-1.0;
    }
};

// S^2 x S^1 in R^5 (Introduction)
template <>
struct TestCaseImpl<TestCase::BOOK_S2XS1> {
    static constexpr size_t NDIM = 5;
    static constexpr size_t KDIM = 2;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.23, -1.21, -1.22, -0.61, -0.62};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.19, 1.17, 1.18, 0.59, 0.6};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{12, 12, 12, 8, 8};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.48, 0.6, 0.64, 0.3, 0.4};
    static void func(const double *x, double *f) {
        f[0] = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]-1.0;
        f[1] = x[3]*x[3]+x[4]*x[4]-0.25;
    }
};

// trilinear interpolant of corner values on [0,1]^3 (ambiguous cube of the Introduction)
inline double trilinear(const double *x, const double *c) {
    double s = 0;
    for (int b = 0; b < 8; ++b) {
        double w = 1;
        for (int i = 0; i < 3; ++i) w *= ((b >> i) & 1) ? x[i] : 1 - x[i];
        s += w * c[b];
    }
    return s;
}
// positive corners (1,0,0) and (0,1,1): not on the diagonal shared by the Kuhn simplices
template <>
struct TestCaseImpl<TestCase::BOOK_AMBIG_SPLIT> {
    static constexpr size_t NDIM = 3;
    static constexpr size_t KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{0, 0, 0};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1, 1, 1};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{1, 1, 1};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.8, 0.1, 0.1};
    static void func(const double *x, double *f) {
        static const double c[8] = {-1, 1, -1, -1, -1, -1, 1, -1};
        f[0] = trilinear(x, c);
    }
};
// positive corners (0,0,0) and (1,1,1): ends of the diagonal shared by all Kuhn simplices
template <>
struct TestCaseImpl<TestCase::BOOK_AMBIG_TUBE> {
    static constexpr size_t NDIM = 3;
    static constexpr size_t KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{0, 0, 0};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1, 1, 1};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{1, 1, 1};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.1, 0.1, 0.1};
    static void func(const double *x, double *f) {
        static const double c[8] = {1, -1, -1, -1, -1, -1, -1, 1};
        f[0] = trilinear(x, c);
    }
};

template <>
struct TestCaseImpl<TestCase::BOOK_S2XS1_COARSE> : TestCaseImpl<TestCase::BOOK_S2XS1> {
    static constexpr std::array<size_t, 5> DOMAIN_DIV{7, 7, 7, 5, 5};
};

// ---- J1a examples: M = F^{-1}(0) refined around N = M ∩ sub^{-1}(0)
template <>
struct TestCaseImpl<TestCase::BOOK_J1A_ELLIPSE> {          // ellipse and a line (Fig. 4.5)
    static constexpr size_t NDIM = 2, KDIM = 1, SUBDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{0, 0};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{8, 8};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{4, 4};
    static constexpr std::array<double, NDIM> FIRST_POINT{6.55, 3.95};
    static void func(const double *x, double *f) {
        f[0] = (x[0]-4.05)*(x[0]-4.05)/6.25 + (x[1]-3.95)*(x[1]-3.95)/4.0 - 1.0;
    }
    static void sub(const double *x, double *g) { g[0] = x[1] - 4.6 - 0.15*(x[0]-4.0); }
};
template <>
struct TestCaseImpl<TestCase::BOOK_J1A_TORUS> {            // torus refined around a wavy curve on it
    static constexpr size_t NDIM = 3, KDIM = 1, SUBDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.63, -1.61, -0.62};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.57, 1.59, 0.58};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{8, 8, 3};
    static constexpr std::array<double, NDIM> FIRST_POINT{1.4, 0.0, 0.0};
    static void func(const double *x, double *f) {
        const double r = std::sqrt(x[0]*x[0] + x[1]*x[1]) - 1.0;
        f[0] = r*r + x[2]*x[2] - 0.16;
    }
    static void sub(const double *x, double *g) { g[0] = x[2] - 0.25*std::sin(5.0*std::atan2(x[1], x[0])); }
};
template <>
struct TestCaseImpl<TestCase::BOOK_J1A_S3> {               // S^3 in R^4 refined around a 2-sphere
    static constexpr size_t NDIM = 4, KDIM = 1, SUBDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.2, -1.2, -1.2, -1.2};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.2, 1.2, 1.2, 1.2};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{4, 4, 4, 4};
    static constexpr std::array<double, NDIM> FIRST_POINT{1.0, 0.0, 0.0, 0.0};
    static void func(const double *x, double *f) { f[0] = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3] - 1.0; }
    static void sub(const double *x, double *g) { g[0] = x[3] - 0.3 - 0.1*x[0]; }
};

template <>
struct TestCaseImpl<TestCase::BOOK_J1A_ZCOSW> {            // z = cos w refined around the real curve Im w = 0.05
    static constexpr size_t NDIM = 4, KDIM = 2, SUBDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-3.17, -3.16, -3.18, -3.15};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{3.11, 3.12, 3.10, 3.13};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{5, 5, 5, 5};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.3, 0.2, 0.9745069929868755, -0.05949885707931209};
    static void func(const double *x, double *f) {
        f[0] = x[2]-std::cos(x[0])*std::cosh(x[1]);
        f[1] = x[3]+std::sin(x[0])*std::sinh(x[1]);
    }
    static void sub(const double *x, double *g) { g[0] = x[1] - 0.05; }
};

// ---------------- implicit modelling (Chapter 13) ----------------
namespace modeling {
inline double uminp(const double *f, int r, double p) { double s = 0; for (int i = 0; i < r; ++i) s += std::pow(f[i], -p); return std::pow(s, -1.0/p); }
inline double imaxp(const double *f, int r, double p) { double s = 0; for (int i = 0; i < r; ++i) s += std::pow(f[i], p); return std::pow(s, 1.0/p); }
inline double sq3(const double *x, double a, double b, double c) { return (x[0]-a)*(x[0]-a)+(x[1]-b)*(x[1]-b)+(x[2]-c)*(x[2]-c); }
}
template <>
struct TestCaseImpl<TestCase::MOD_TWO_SPHERES> {          // sharp union of two spheres
    static constexpr size_t NDIM = 3, KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.13, -1.12, -1.11};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{2.11, 1.14, 1.13};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{40, 28, 28};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.05, 0.594576344266849, 0.802483003457562};
    static void func(const double *x, double *f) { f[0] = std::min(modeling::sq3(x,0,0,0), modeling::sq3(x,1,0,0)) - 1.0; }
};
template <>
struct TestCaseImpl<TestCase::MOD_TWO_SPHERES_P8> {       // smooth union, p = 8
    static constexpr size_t NDIM = 3, KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.13, -1.12, -1.11};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{2.11, 1.14, 1.13};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{40, 28, 28};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.5, 0.88863072950733, 0.225484277147053};
    static void func(const double *x, double *f) { double g[2] = {modeling::sq3(x,0,0,0)+1e-12, modeling::sq3(x,1,0,0)+1e-12}; f[0] = modeling::uminp(g, 2, 8.0) - 1.0; }
};
template <>
struct TestCaseImpl<TestCase::MOD_EIGHT_TORI> {           // smooth union of eight tori, p = 10
    static constexpr size_t NDIM = 3, KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-9.13, -9.12, -2.11};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{9.11, 9.14, 2.13};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{80, 80, 20};
    static constexpr std::array<double, NDIM> FIRST_POINT{8.99964996624341, 0.03, 0.02};
    static void func(const double *x, double *f) {
        static const double c[8][2] = {{6,0},{-6,0},{0,6},{0,-6},{4,4},{4,-4},{-4,4},{-4,-4}};
        double g[8];
        for (int i = 0; i < 8; ++i) { double t = std::hypot(x[0]-c[i][0], x[1]-c[i][1]) - 2.0; g[i] = t*t + x[2]*x[2] + 1e-12; }
        f[0] = modeling::uminp(g, 8, 10.0) - 1.0;
    }
};
template <>
struct TestCaseImpl<TestCase::MOD_KNOT> {                 // smooth union of 151 spheres on a knot, p = 30
    static constexpr size_t NDIM = 3, KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.93, -1.92, -1.11};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.91, 1.94, 1.13};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{72, 72, 42};
    static constexpr std::array<double, NDIM> FIRST_POINT{1.41, 0.02, 0.356633311006689};
    static void func(const double *x, double *f) {
        static double C[151][3]; static bool init = false;
        if (!init) { for (int i = 0; i < 151; ++i) { double t = 4*M_PI*i/150; C[i][0] = (1+0.4*std::cos(1.5*t))*std::cos(t); C[i][1] = (1+0.4*std::cos(1.5*t))*std::sin(t); C[i][2] = 0.5*std::sin(1.5*t); } init = true; }
        double g[151];
        for (int i = 0; i < 151; ++i) g[i] = modeling::sq3(x, C[i][0], C[i][1], C[i][2]) - 0.24*0.24 + 1.0;
        f[0] = modeling::uminp(g, 151, 30.0) - 1.0;
    }
};
template <>
struct TestCaseImpl<TestCase::MOD_INTERSECT3> {           // smooth intersection of three spheres, p = 10
    static constexpr size_t NDIM = 3, KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.13, -1.12, -1.11};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{2.11, 2.14, 1.13};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{64, 64, 44};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.3, 0.4, 0.720568923735442};
    static void func(const double *x, double *f) { double g[3] = {modeling::sq3(x,0,0,0), modeling::sq3(x,0,1,0), modeling::sq3(x,std::sqrt(3.0)/2,0.5,0)}; f[0] = modeling::imaxp(g, 3, 10.0) - 1.0; }
};
template <>
struct TestCaseImpl<TestCase::MOD_R4_TWO_BALLS> {         // smooth union of two 3-spheres in R^4, p = 8
    static constexpr size_t NDIM = 4, KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.93, -1.12, -1.11, -1.13};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.91, 1.14, 1.13, 1.12};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{36, 22, 22, 22};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.8, 0.595280413843421, 0.803458819186495, 0.01};
    static void func(const double *x, double *f) {
        double a = (x[0]-0.8)*(x[0]-0.8), b = (x[0]+0.8)*(x[0]+0.8), r = x[1]*x[1]+x[2]*x[2]+x[3]*x[3];
        double g[2] = {a + r + 1e-12, b + r + 1e-12}; f[0] = modeling::uminp(g, 2, 8.0) - 1.0;
    }
};
template <>
struct TestCaseImpl<TestCase::MOD_R4_DIFF> {              // 4-ball minus the solid cylinder x1^2+x2^2 < 0.45^2, p-norm difference (p = 6)
    static constexpr size_t NDIM = 4, KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.13, -1.12, -1.11, -1.14};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.11, 1.14, 1.13, 1.12};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{32, 32, 32, 32};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.5, 0.6, 0.569938091651778, 0.254847774396861};
    static void func(const double *x, double *f) {
        double b = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3] + 1e-12, c = (x[0]*x[0]+x[1]*x[1])/0.2025 + 1e-12;
        double g[2] = {b, 1.0/c}; f[0] = modeling::imaxp(g, 2, 6.0) - 1.0;
    }
};
template <>
struct TestCaseImpl<TestCase::MOD_R4_TUBE> {              // smooth union of 200 3-spheres along a closed curve of R^4, p = 30
    static constexpr size_t NDIM = 4, KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.83, -1.82, -0.91, -0.96};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.81, 1.84, 0.93, 0.94};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{44, 44, 22, 22};
    static constexpr std::array<double, NDIM> FIRST_POINT{1.35, 0.01, 0.338163248211567, 0.513};
    static void func(const double *x, double *f) {
        static double C[200][4]; static bool init = false;
        if (!init) { for (int i = 0; i < 200; ++i) { double s = 2*M_PI*i/200; C[i][0] = (1+0.35*std::cos(3*s))*std::cos(2*s); C[i][1] = (1+0.35*std::cos(3*s))*std::sin(2*s); C[i][2] = 0.45*std::sin(3*s); C[i][3] = 0.5*std::cos(5*s); } init = true; }
        double g[200];
        for (int i = 0; i < 200; ++i) { double d = 0; for (int j = 0; j < 4; ++j) d += (x[j]-C[i][j])*(x[j]-C[i][j]); g[i] = d - 0.22*0.22 + 1.0; }
        f[0] = modeling::uminp(g, 200, 30.0) - 1.0;
    }
};

namespace modeling {
// pipe/box model after Moreira (2000): box and outer cylinder, minus three inner cylinders.
// primitives are written as positive functions g with g = 1 on the surface, g < 1 inside;
// union = min / p-norm min, intersection = max / p-norm max, complement = 1/g.
inline double pipebox(const double *x, double p) {
    auto U = [p](double a, double b) { return p <= 0 ? std::min(a, b) : std::pow(std::pow(a, -p) + std::pow(b, -p), -1.0/p); };
    auto I = [p](double a, double b) { return p <= 0 ? std::max(a, b) : std::pow(std::pow(a, p) + std::pow(b, p), 1.0/p); };
    auto sq = [](double t) { return t*t; };
    const double e = 1e-12;
    double box = I(I(sq(x[0]/1.25)+e, sq(x[1]/1.25)+e), sq(x[2]/1.25)+e);
    double cylo = I((sq(x[0])+sq(x[1]))/0.81+e, sq(x[2]/3.0)+e);
    double cz = I((sq(x[0])+sq(x[1]))/0.25+e, sq(x[2]/3.125)+e);
    double cx = I((sq(x[1])+sq(x[2]))/0.25+e, sq(x[0]/1.375)+e);
    double cy = I((sq(x[0])+sq(x[2]))/0.25+e, sq(x[1]/1.375)+e);
    double solid = U(box, cylo), holes = U(U(cz, cx), cy);
    return I(solid, 1.0/holes);
}
}
template <>
struct TestCaseImpl<TestCase::MOD_PIPEBOX> {
    static constexpr size_t NDIM = 3, KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.63, -1.62, -3.31};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.61, 1.64, 3.29};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{56, 56, 112};
    static constexpr std::array<double, NDIM> FIRST_POINT{1.25, 0.71, 0.73};
    static void func(const double *x, double *f) { f[0] = modeling::pipebox(x, 0.0) - 1.0; }
};
template <>
struct TestCaseImpl<TestCase::MOD_PIPEBOX_P8> {
    static constexpr size_t NDIM = 3, KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.63, -1.62, -3.31};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.61, 1.64, 3.29};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{56, 56, 112};
    static constexpr std::array<double, NDIM> FIRST_POINT{1.25001916, 0.71, 0.73};
    static void func(const double *x, double *f) { f[0] = modeling::pipebox(x, 8.0) - 1.0; }
};
template <>
struct TestCaseImpl<TestCase::BOOK_HYPERBOLA4> {          // generalized hyperbola x1 x2 x3 x4 = 1e-4 (GCMH article)
    static constexpr size_t NDIM = 4, KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1, -1, -1, -1};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1, 1, 1, 1};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{9, 9, 9, 9};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.5, 0.4, 0.3, 0.0016666666666666668};
    static void func(const double *x, double *f) { f[0] = x[0]*x[1]*x[2]*x[3] - 1e-4; }
};
// three-torus in R^6 on a coarser grid (benchmark of Chapter 12)
template <>
struct TestCaseImpl<TestCase::BOOK_T3_R6_B> {
    static constexpr size_t NDIM = 6, KDIM = 3;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-0.6, -0.6, -1.1, -1.1, -4.1, -4.1};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{0.6, 0.6, 1.1, 1.1, 4.1, 4.1};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{8, 8, 8, 8, 8, 8};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.4, 0.3, 0.6, 0.8, 1.2, 1.6};
    static void func(const double *x, double *f) {
        f[0] = x[0]*x[0] + x[1]*x[1] - 0.25;
        f[1] = x[2]*x[2] + x[3]*x[3] - 1.0;
        f[2] = x[4]*x[4] + x[5]*x[5] - 4.0;
    }
};
// convergence studies (Chapter 12): grid size CONV_N given at compile time
#ifndef CONV_N
#define CONV_N 16
#endif
template <>
struct TestCaseImpl<TestCase::BOOK_CONV_S2> {
    static constexpr size_t NDIM = 3, KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.137, -1.123, -1.111};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.113, 1.127, 1.139};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{CONV_N, CONV_N, CONV_N};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.48, 0.6, 0.64};
    static void func(const double *x, double *f) { f[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 1.0; }
};
template <>
struct TestCaseImpl<TestCase::BOOK_CONV_T2> {               // torus in R^4: x1^2+x2^2=1, x3^2+x4^2=1/4
    static constexpr size_t NDIM = 4, KDIM = 2;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.137, -1.123, -0.611, -0.597};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{1.113, 1.127, 0.589, 0.603};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{CONV_N, CONV_N, CONV_N/2, CONV_N/2};
    static constexpr std::array<double, NDIM> FIRST_POINT{0.6, 0.8, 0.3, 0.4};
    static void func(const double *x, double *f) { f[0] = x[0]*x[0] + x[1]*x[1] - 1.0; f[1] = x[2]*x[2] + x[3]*x[3] - 0.25; }
};
