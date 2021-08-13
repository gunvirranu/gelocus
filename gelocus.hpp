/*
 * gelocus -- Earth-centered reference frames in modern C++
 *
 * "Scientific code that shouldn't feel like it was written by scientists!"
 *
 * Licensed under the MIT License <http://opensource.org/licenses/MIT>.
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2021 Gunvir Ranu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef GELOCUS_HPP
#define GELOCUS_HPP

#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>

#ifndef GELOCUS_VEC3
#error Must define `GELOCUS_VEC3` type
#endif
#ifndef GELOCUS_MATRIX3
#error Must define `GELOCUS_MATRIX3` type
#endif
#ifndef GELOCUS_MATRIX3_TRANSPOSED
#error Must define `GELOCUS_MATRIX3_TRANSPOSED` method
#endif

namespace gelocus {

using std::size_t;
using Vec = GELOCUS_VEC3;
using Matrix = GELOCUS_MATRIX3;

constexpr double PI = 3.141592653589793238462643383279;
constexpr double RAD_PER_DEG = 0.017453292519943295769;
constexpr double DEG_PER_RAD = 57.29577951308232087680;

enum class Frame {
    J2000,
    MOD,
    TOD,
    TEME,
    PEF,
    ECEF
};

std::ostream& operator<<(std::ostream &os, Frame f);

struct EOPData {
    double xp, yp;  // [rad]
    double dPsi, dEps;
    // Not used for IAU-76/FK5 transformations
    double dX, dY;

    EOPData() : xp(0), yp(0), dPsi(0), dEps(0), dX(0), dY(0) {}
};

template <Frame F>
class Position {
public:
    Vec vec;

    explicit Position(Vec v) : vec(v) {}

    template <Frame To>
    Position<To> transform(double jd, EOPData eop = EOPData()) const;

    const auto& operator()(size_t i) const {
        return vec(i);
    }

    auto& operator()(size_t i) {
        return vec(i);
    }

    friend std::ostream& operator<<(std::ostream &os, const Position &p) {
        return os << "Position<" << F << "> { " << p(0) << ", " << p(1) << ", " << p(2) << " }";
    }
};

template <Frame From, Frame To>
class Transformation {
public:
    Matrix mat;

    explicit Transformation(Matrix mat) : mat(mat) {}

    explicit Transformation(double jd, EOPData eop = EOPData());

    Transformation<To, From> inverse() const {
        Matrix mat_T = this->mat.GELOCUS_MATRIX3_TRANSPOSED();
        return Transformation<To, From>(mat_T);
    }

    // TODO: Think about impling `operator*` as `Mat * Vec`
    Position<To> apply(const Position<From> &pos) const {
        return Position<To>(this->mat * pos.vec);
    }

    template <Frame FromR>
    Transformation<FromR, To> operator*(const Transformation<FromR, From> &rhs) const {
        return Transformation<FromR, To>(this->mat * rhs.mat);
    }

    friend std::ostream& operator<<(std::ostream &os, const Transformation &A) {
        return os << "Transformation<" << From << ", " << To << "> { " << A.mat << " }";
    }
};

std::ostream& operator<<(std::ostream &os, Frame f) {
    switch (f) {
        case Frame::J2000: return os << "J2000";
        case Frame::MOD: return os << "MOD";
        case Frame::TOD: return os << "TOD";
        case Frame::TEME: return os << "TEME";
        case Frame::PEF: return os << "PEF";
        case Frame::ECEF: return os << "ECEF";
    };
    throw std::invalid_argument("Invalid frame type");
}

template <Frame F>
template <Frame To>
Position<To> Position<F>::transform(const double jd, const EOPData eop) const {
    const auto A = Transformation<F, To>(jd, eop);
    return A.apply(*this);
}

constexpr double rad_to_deg(const double rad) {
    return rad * DEG_PER_RAD;
}

constexpr double deg_to_rad(const double deg) {
    return deg * RAD_PER_DEG;
}

// For pseudo-"private" implementation details
namespace detail {

// Common

constexpr double jd_to_jc(const double jd) {
    return (jd - 2451545) / 36525;
}

// FK5 / IAU-76 Theory

double greenwich_mean_sidereal_time(const double jc_ut1) {
    constexpr double cs[] = {
        67310.54841, 876600.0 * 3600 + 8640184.812866, 0.093104, -6.2e-6
    };
    // Horner's method for polynomial evaluation
    const double gmst_secs = cs[0] + jc_ut1 * (cs[1] + jc_ut1 * (cs[2] + jc_ut1 * cs[3]));
    double gmst = std::fmod(deg_to_rad(gmst_secs / 240), 2 * PI);
    if (gmst < 0) {
        gmst += 2 * PI;
    }
    return gmst;
}

// MOD to J2000
void iau76_precession(const double jc, Matrix &P) {
    constexpr double arcsec_to_rad = PI / (180.0 * 3600);
    // All in [arcsecond], uses Horner's method
    double zeta = jc * (2306.2181 + jc * (0.30188 + jc * 0.017998));
    double theta = jc * (2004.3109 + jc * (-0.42665 + jc * -0.041833));
    double z = jc * (2306.2181 + jc * (1.09468 + jc * 0.018203));
    zeta *= arcsec_to_rad;
    theta *= arcsec_to_rad;
    z *= arcsec_to_rad;

    const double cos_zeta = std::cos(zeta);
    const double sin_zeta = std::sin(zeta);
    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);
    const double cos_z = std::cos(z);
    const double sin_z = std::sin(z);
    P(0, 0) = cos_zeta * cos_theta * cos_z - sin_zeta * sin_z;
    P(0, 1) = cos_zeta * cos_theta * sin_z + sin_zeta * cos_z;
    P(0, 2) = cos_zeta * sin_theta;
    P(1, 0) = -sin_zeta * cos_theta * cos_z - cos_zeta * sin_z;
    P(1, 1) = -sin_zeta * cos_theta * sin_z + cos_zeta * cos_z;
    P(1, 2) = -sin_zeta * sin_theta;
    P(2, 0) = -sin_theta * cos_z;
    P(2, 1) = -sin_theta * sin_z;
    P(2, 2) = cos_theta;
}

struct Iau80NutCoeffSet {
    int l, l1, f, d, omg;   // Coefficients for l, l1, f, d, and omega
    double sp, spt;         // Longitude sine coefficients in 0.1 mas
    double ce, cet;         // Obliquity cosine coefficients in 0.1 mas
    int idx;                // Original index of coefficient set
};
constexpr size_t IAU80_TOTAL_NUT_TERMS = 106;
constexpr struct Iau80NutCoeffSet IAU80_NUT_COEFFS[IAU80_TOTAL_NUT_TERMS] = {
    // Ordered by largest coefficients first. Original index is `idx`.
    {  0,  0,  0,  0,  1, -171996.0, -174.2, 92025.0,  8.9,   1 },
    {  0,  0,  2, -2,  2,  -13187.0,   -1.6,  5736.0, -3.1,   9 },
    {  0,  0,  2,  0,  2,   -2274.0,   -0.2,   977.0, -0.5,  31 },
    {  0,  0,  0,  0,  2,    2062.0,    0.2,  -895.0,  0.5,   2 },
    {  0,  1,  0,  0,  0,    1426.0,   -3.4,    54.0, -0.1,  10 },
    {  1,  0,  0,  0,  0,     712.0,    0.1,    -7.0,  0.0,  32 },
    {  0,  1,  2, -2,  2,    -517.0,    1.2,   224.0, -0.6,  11 },
    {  0,  0,  2,  0,  1,    -386.0,   -0.4,   200.0,  0.0,  33 },
    {  1,  0,  2,  0,  2,    -301.0,    0.0,   129.0, -0.1,  34 },
    {  0, -1,  2, -2,  2,     217.0,   -0.5,   -95.0,  0.3,  12 },
    {  1,  0,  0, -2,  0,    -158.0,    0.0,    -1.0,  0.0,  35 },
    {  0,  0,  2, -2,  1,     129.0,    0.1,   -70.0,  0.0,  13 },
    { -1,  0,  2,  0,  2,     123.0,    0.0,   -53.0,  0.0,  36 },
    {  1,  0,  0,  0,  1,      63.0,    0.1,   -33.0,  0.0,  38 },
    {  0,  0,  0,  2,  0,      63.0,    0.0,    -2.0,  0.0,  37 },
    { -1,  0,  2,  2,  2,     -59.0,    0.0,    26.0,  0.0,  40 },
    { -1,  0,  0,  0,  1,     -58.0,   -0.1,    32.0,  0.0,  39 },
    {  1,  0,  2,  0,  1,     -51.0,    0.0,    27.0,  0.0,  41 },
    {  2,  0,  0, -2,  0,      48.0,    0.0,     1.0,  0.0,  14 },
    { -2,  0,  2,  0,  1,      46.0,    0.0,   -24.0,  0.0,   3 },
    {  0,  0,  2,  2,  2,     -38.0,    0.0,    16.0,  0.0,  42 },
    {  2,  0,  2,  0,  2,     -31.0,    0.0,    13.0,  0.0,  45 },
    {  2,  0,  0,  0,  0,      29.0,    0.0,    -1.0,  0.0,  43 },
    {  1,  0,  2, -2,  2,      29.0,    0.0,   -12.0,  0.0,  44 },
    {  0,  0,  2,  0,  0,      26.0,    0.0,    -1.0,  0.0,  46 },
    {  0,  0,  2, -2,  0,     -22.0,    0.0,     0.0,  0.0,  15 },
    { -1,  0,  2,  0,  1,      21.0,    0.0,   -10.0,  0.0,  47 },
    {  0,  2,  0,  0,  0,      17.0,   -0.1,     0.0,  0.0,  16 },
    {  0,  2,  2, -2,  2,     -16.0,    0.1,     7.0,  0.0,  18 },
    { -1,  0,  0,  2,  1,      16.0,    0.0,    -8.0,  0.0,  48 },
    {  0,  1,  0,  0,  1,     -15.0,    0.0,     9.0,  0.0,  17 },
    {  1,  0,  0, -2,  1,     -13.0,    0.0,     7.0,  0.0,  49 },
    {  0, -1,  0,  0,  1,     -12.0,    0.0,     6.0,  0.0,  19 },
    {  2,  0, -2,  0,  0,      11.0,    0.0,     0.0,  0.0,   4 },
    { -1,  0,  2,  2,  1,     -10.0,    0.0,     5.0,  0.0,  50 },
    {  1,  0,  2,  2,  2,      -8.0,    0.0,     3.0,  0.0,  54 },
    {  0, -1,  2,  0,  2,      -7.0,    0.0,     3.0,  0.0,  53 },
    {  0,  0,  2,  2,  1,      -7.0,    0.0,     3.0,  0.0,  58 },
    {  1,  1,  0, -2,  0,      -7.0,    0.0,     0.0,  0.0,  51 },
    {  0,  1,  2,  0,  2,       7.0,    0.0,    -3.0,  0.0,  52 },
    { -2,  0,  0,  2,  1,      -6.0,    0.0,     3.0,  0.0,  20 },
    {  0,  0,  0,  2,  1,      -6.0,    0.0,     3.0,  0.0,  57 },
    {  2,  0,  2, -2,  2,       6.0,    0.0,    -3.0,  0.0,  56 },
    {  1,  0,  0,  2,  0,       6.0,    0.0,     0.0,  0.0,  55 },
    {  1,  0,  2, -2,  1,       6.0,    0.0,    -3.0,  0.0,  58 },
    {  0,  0,  0, -2,  1,      -5.0,    0.0,     3.0,  0.0,  60 },
    {  0, -1,  2, -2,  1,      -5.0,    0.0,     3.0,  0.0,  21 },
    {  2,  0,  2,  0,  1,      -5.0,    0.0,     3.0,  0.0,  62 },
    {  1, -1,  0,  0,  0,       5.0,    0.0,     0.0,  0.0,  61 },
    {  1,  0,  0, -1,  0,      -4.0,    0.0,     0.0,  0.0,  24 },
    {  0,  0,  0,  1,  0,      -4.0,    0.0,     0.0,  0.0,  65 },
    {  0,  1,  0, -2,  0,      -4.0,    0.0,     0.0,  0.0,  63 },
    {  1,  0, -2,  0,  0,       4.0,    0.0,     0.0,  0.0,  64 },
    {  2,  0,  0, -2,  1,       4.0,    0.0,    -2.0,  0.0,  22 },
    {  0,  1,  2, -2,  1,       4.0,    0.0,    -2.0,  0.0,  23 },
    {  1,  1,  0,  0,  0,      -3.0,    0.0,     0.0,  0.0,  66 },
    {  1, -1,  0, -1,  0,      -3.0,    0.0,     0.0,  0.0,   6 },
    { -1, -1,  2,  2,  2,      -3.0,    0.0,     1.0,  0.0,  69 },
    {  0, -1,  2,  2,  2,      -3.0,    0.0,     1.0,  0.0,  72 },
    {  1, -1,  2,  0,  2,      -3.0,    0.0,     1.0,  0.0,  68 },
    {  3,  0,  2,  0,  2,      -3.0,    0.0,     1.0,  0.0,  71 },
    { -2,  0,  2,  0,  2,      -3.0,    0.0,     1.0,  0.0,   5 },
    {  1,  0,  2,  0,  0,       3.0,    0.0,     0.0,  0.0,  67 },
    { -1,  0,  2,  4,  2,      -2.0,    0.0,     1.0,  0.0,  82 },
    {  1,  0,  0,  0,  2,      -2.0,    0.0,     1.0,  0.0,  76 },
    { -1,  0,  2, -2,  1,      -2.0,    0.0,     1.0,  0.0,  74 },
    {  0, -2,  2, -2,  1,      -2.0,    0.0,     1.0,  0.0,   7 },
    { -2,  0,  0,  0,  1,      -2.0,    0.0,     1.0,  0.0,  70 },
    {  2,  0,  0,  0,  1,       2.0,    0.0,    -1.0,  0.0,  75 },
    {  3,  0,  0,  0,  0,       2.0,    0.0,     0.0,  0.0,  77 },
    {  1,  1,  2,  0,  2,       2.0,    0.0,    -1.0,  0.0,  73 },
    {  0,  0,  2,  1,  2,       2.0,    0.0,    -1.0,  0.0,  78 },
    {  1,  0,  0,  2,  1,      -1.0,    0.0,     0.0,  0.0,  91 },
    {  1,  0,  2,  2,  1,      -1.0,    0.0,     1.0,  0.0,  85 },
    {  1,  1,  0, -2,  1,      -1.0,    0.0,     0.0,  0.0, 102 },
    {  0,  1,  0,  2,  0,      -1.0,    0.0,     0.0,  0.0,  99 },
    {  0,  1,  2, -2,  0,      -1.0,    0.0,     0.0,  0.0,  30 },
    {  0,  1, -2,  2,  0,      -1.0,    0.0,     0.0,  0.0,  27 },
    {  1,  0, -2,  2,  0,      -1.0,    0.0,     0.0,  0.0, 103 },
    {  1,  0, -2, -2,  0,      -1.0,    0.0,     0.0,  0.0, 100 },
    {  1,  0,  2, -2,  0,      -1.0,    0.0,     0.0,  0.0,  94 },
    {  1,  0,  0, -4,  0,      -1.0,    0.0,     0.0,  0.0,  80 },
    {  2,  0,  0, -4,  0,      -1.0,    0.0,     0.0,  0.0,  83 },
    {  0,  0,  2,  4,  2,      -1.0,    0.0,     0.0,  0.0, 105 },
    {  0,  0,  2, -1,  2,      -1.0,    0.0,     0.0,  0.0,  98 },
    { -2,  0,  2,  4,  2,      -1.0,    0.0,     1.0,  0.0,  86 },
    {  2,  0,  2,  2,  2,      -1.0,    0.0,     0.0,  0.0,  90 },
    {  0, -1,  2,  0,  1,      -1.0,    0.0,     0.0,  0.0, 101 },
    {  0,  0, -2,  0,  1,      -1.0,    0.0,     0.0,  0.0,  97 },
    {  0,  0,  4, -2,  2,       1.0,    0.0,     0.0,  0.0,  92 },
    {  0,  1,  0,  0,  2,       1.0,    0.0,     0.0,  0.0,  28 },
    {  1,  1,  2, -2,  2,       1.0,    0.0,    -1.0,  0.0,  84 },
    {  3,  0,  2, -2,  2,       1.0,    0.0,     0.0,  0.0,  93 },
    { -2,  0,  2,  2,  2,       1.0,    0.0,    -1.0,  0.0,  81 },
    { -1,  0,  0,  0,  2,       1.0,    0.0,    -1.0,  0.0,  79 },
    {  0,  0, -2,  2,  1,       1.0,    0.0,     0.0,  0.0,  26 },
    {  0,  1,  2,  0,  1,       1.0,    0.0,     0.0,  0.0,  95 },
    { -1,  0,  4,  0,  2,       1.0,    0.0,     0.0,  0.0,  87 },
    {  2,  1,  0, -2,  0,       1.0,    0.0,     0.0,  0.0,  25 },
    {  2,  0,  0,  2,  0,       1.0,    0.0,     0.0,  0.0, 104 },
    {  2,  0,  2, -2,  1,       1.0,    0.0,    -1.0,  0.0,  89 },
    {  2,  0, -2,  0,  1,       1.0,    0.0,     0.0,  0.0,   8 },
    {  1, -1,  0, -2,  0,       1.0,    0.0,     0.0,  0.0,  88 },
    { -1,  0,  0,  1,  1,       1.0,    0.0,     0.0,  0.0,  29 },
    { -1, -1,  0,  2,  1,       1.0,    0.0,     0.0,  0.0,  96 },
    {  0,  1,  0,  1,  0,       1.0,    0.0,     0.0,  0.0, 106 }
};

// TOD to MOD
void iau80_nutation(
    const double jc, const EOPData eop,
    Matrix &N, double &mean_eps, double &omega, double &delta_psi
) {
    constexpr double P1MAS_TO_RADS = 1e-4 * PI / (180.0 * 3600);
    constexpr double inv3600 = 1.0 / 3600;
    mean_eps = 84381.448 + jc * (-46.8150 + jc * (-0.00059 + jc * 0.001813));
    mean_eps = deg_to_rad(std::fmod(mean_eps * inv3600, 360));
    // Delaunay fundamental arguments in degrees
    double l = 134.96298139 + jc * (1717915922.6330 + jc * (31.310 + jc * 0.064)) * inv3600;
    double l1 = 357.52772333 + jc * (129596581.2240 + jc * (-0.577 + jc * -0.012)) * inv3600;
    double f = 93.27191028 + jc * (1739527263.1370 + jc * (-13.257 + jc * 0.011)) * inv3600;
    double d = 297.85036306 + jc * (1602961601.3280 + jc * (-6.891 + jc * 0.019)) * inv3600;
    omega = 125.04452222 + jc * (-6962890.5390 + jc * (7.455 + jc * 0.008)) * inv3600;
    l = deg_to_rad(std::fmod(l, 360));
    l1 = deg_to_rad(std::fmod(l1, 360));
    f = deg_to_rad(std::fmod(f, 360));
    d = deg_to_rad(std::fmod(d, 360));
    omega = deg_to_rad(std::fmod(omega, 360));

    delta_psi = 0;
    double delta_eps = 0;
    // Sum in reverse order to preserve floating point accuracy
    for (int i = 106 - 1; i >= 0; --i) {
        constexpr auto& NUT_COEFFS = IAU80_NUT_COEFFS;
        const double arg = l * NUT_COEFFS[i].l
            + l1 * NUT_COEFFS[i].l1
            + f * NUT_COEFFS[i].f
            + d * NUT_COEFFS[i].d
            + omega * NUT_COEFFS[i].omg;
        delta_psi += (NUT_COEFFS[i].sp + jc * NUT_COEFFS[i].spt) * std::sin(arg);
        delta_eps += (NUT_COEFFS[i].ce + jc * NUT_COEFFS[i].cet) * std::cos(arg);
    }
    // Add in EOP corrections to GCRF
    delta_psi += eop.dPsi;
    delta_eps += eop.dEps;
    delta_psi = std::fmod(delta_psi * P1MAS_TO_RADS, 2 * PI);
    delta_eps = std::fmod(delta_eps * P1MAS_TO_RADS, 2 * PI);
    const double true_eps = mean_eps + delta_eps;

    const double sin_psi = std::sin(delta_psi);
    const double cos_psi = std::cos(delta_psi);
    const double sin_eps = std::sin(mean_eps);
    const double cos_eps = std::cos(mean_eps);
    const double sin_true_eps = std::sin(true_eps);
    const double cos_true_eps = std::cos(true_eps);
    N(0, 0) = cos_psi;
    N(0, 1) = cos_true_eps * sin_psi;
    N(0, 2) = sin_true_eps * sin_psi;
    N(1, 0) = -cos_eps * sin_psi;
    N(1, 1) = cos_true_eps * cos_eps * cos_psi + sin_true_eps * sin_eps;
    N(1, 2) = sin_true_eps * cos_eps * cos_psi - cos_true_eps * sin_eps;
    N(2, 0) = -sin_eps * sin_psi;
    N(2, 1) = cos_true_eps * sin_eps * cos_psi - sin_true_eps * cos_eps;
    N(2, 2) = sin_true_eps * sin_eps * cos_psi + cos_true_eps * cos_eps;
}

// PEF to TOD
void iau76_sidereal(
    const double jd, const double mean_eps, const double omega, const double delta_psi,
    Matrix &S
) {
    const double jc = jd_to_jc(jd);
    const double gmst = greenwich_mean_sidereal_time(jc);
    double ast = gmst + delta_psi * std::cos(mean_eps);
    if (jd > 2450449.5) {
        ast += 0.002640 * PI / (3600 * 180) * std::sin(omega);
        ast += 0.000063 * PI / (3600 * 180) * std::sin(2 * omega);
    }
    ast = std::fmod(ast, 2 * PI);
    const double sin_ast = std::sin(ast);
    const double cos_ast = std::cos(ast);
    S(0, 0) = cos_ast; S(0, 1) = -sin_ast; S(0, 2) = 0;
    S(1, 0) = sin_ast; S(1, 1) =  cos_ast; S(1, 2) = 0;
    S(2, 0) = 0      ; S(2, 1) =  0      ; S(2, 2) = 1;
}

// ECEF to PEF
void fk5_polar_motion(const EOPData eop, Matrix &PM) {
    const double sin_xp = std::sin(eop.xp);
    const double cos_xp = std::cos(eop.xp);
    const double sin_yp = std::sin(eop.yp);
    const double cos_yp = std::cos(eop.yp);
    PM(0, 0) = cos_xp         ; PM(0, 1) =  0     ; PM(0, 2) = -sin_xp         ;
    PM(1, 0) = sin_xp * sin_yp; PM(1, 1) =  cos_yp; PM(1, 2) =  cos_xp * sin_yp;
    PM(2, 0) = sin_xp * cos_yp; PM(2, 1) = -sin_yp; PM(2, 2) =  cos_xp * cos_yp;
}

void teme_to_pef(const double jd, Matrix &S) {
    const double jc = jd_to_jc(jd);
    const double gmst = greenwich_mean_sidereal_time(jc);
    const double omega_deg = 125.04452222 + jc * (
        -6962890.5390 + jc * (7.455 + jc * 0.008)
    ) / 3600;
    const double omega = deg_to_rad(std::fmod(omega_deg, 360));
    double gmstg = gmst;
    if (jd > 2450449.5) {
        gmstg += 0.002640 * PI / (3600 * 180) * sin(omega);
        gmstg += 0.000063 * PI / (3600 * 180) * sin(2 * omega);
    }
    gmstg = fmod(gmstg, 2 * PI);

    const double sg = std::sin(gmstg);
    const double cg = std::cos(gmstg);
    S(0, 0) =  cg; S(0, 1) = sg; S(0, 2) = 0;
    S(1, 0) = -sg; S(1, 1) = cg; S(1, 2) = 0;
    S(2, 0) =   0; S(2, 1) =  0; S(2, 2) = 1;
}

} // namespace detail

// Implemented "Basic" Transformations

template <Frame From, Frame To>
Transformation<From, To>::Transformation(const double jd, const EOPData eop) {
    static_assert(From == To, "This transformation is not implemented");
    (void) jd;  // To ignore unused value
    (void) eop;
    // Use  identity matrix since `From == To`
    mat(0, 0) = 1; mat(0, 1) = 0; mat(0, 2) = 0;
    mat(1, 0) = 0; mat(1, 1) = 1; mat(1, 2) = 0;
    mat(2, 0) = 0; mat(2, 1) = 0; mat(2, 2) = 1;
}

template <>
Transformation<Frame::MOD, Frame::J2000>::Transformation(const double jd, const EOPData eop) {
    const double jc = detail::jd_to_jc(jd);
    (void) eop;  // To ignore unused value
    detail::iau76_precession(jc, this->mat);
}

template <>
Transformation<Frame::TOD, Frame::MOD>::Transformation(const double jd, const EOPData eop) {
    const double jc = detail::jd_to_jc(jd);
    double mean_eps, omega, delta_psi;
    detail::iau80_nutation(jc, eop, this->mat, mean_eps, omega, delta_psi);
}

template <>
Transformation<Frame::PEF, Frame::TOD>::Transformation(const double jd, const EOPData eop) {
    const double jc = detail::jd_to_jc(jd);
    double mean_eps, omega, delta_psi;
    Matrix dummy;
    detail::iau80_nutation(jc, eop, dummy, mean_eps, omega, delta_psi);
    detail::iau76_sidereal(jd, mean_eps, omega, delta_psi, this->mat);
}

template <>
Transformation<Frame::ECEF, Frame::PEF>::Transformation(const double jd, const EOPData eop) {
    (void) jd;  // To ignore unused value
    detail::fk5_polar_motion(eop, this->mat);
}

template <>
Transformation<Frame::TEME, Frame::PEF>::Transformation(const double jd, const EOPData eop) {
    (void) eop;  // To ignore unused value
    detail::teme_to_pef(jd, this->mat);
}

// Compound transformations

template <>
Transformation<Frame::TOD, Frame::J2000>::Transformation(const double jd, const EOPData eop) {
    const auto tod_to_mod = Transformation<Frame::TOD, Frame::MOD>(jd, eop);
    const auto mod_to_j2000 = Transformation<Frame::MOD, Frame::J2000>(jd, eop);
    *this = mod_to_j2000 * tod_to_mod;
}

template <>
Transformation<Frame::PEF, Frame::J2000>::Transformation(const double jd, const EOPData eop) {
    const auto pef_to_tod = Transformation<Frame::PEF, Frame::TOD>(jd, eop);
    const auto tod_to_j2000 = Transformation<Frame::TOD, Frame::J2000>(jd, eop);
    *this = tod_to_j2000 * pef_to_tod;
}

template <>
Transformation<Frame::ECEF, Frame::J2000>::Transformation(const double jd, const EOPData eop) {
    const auto ecef_to_pef = Transformation<Frame::ECEF, Frame::PEF>(jd, eop);
    const auto pef_to_j2000 = Transformation<Frame::PEF, Frame::J2000>(jd, eop);
    *this = pef_to_j2000 * ecef_to_pef;
}

template <>
Transformation<Frame::TEME, Frame::J2000>::Transformation(const double jd, const EOPData eop) {
    const auto teme_to_pef = Transformation<Frame::TEME, Frame::PEF>(jd, eop);
    const auto pef_to_j2000 = Transformation<Frame::PEF, Frame::J2000>(jd, eop);
    *this = pef_to_j2000 * teme_to_pef;
}

template <>
Transformation<Frame::PEF, Frame::MOD>::Transformation(const double jd, const EOPData eop) {
    const auto pef_to_tod = Transformation<Frame::PEF, Frame::TOD>(jd, eop);
    const auto tod_to_mod = Transformation<Frame::TOD, Frame::MOD>(jd, eop);
    *this = tod_to_mod * pef_to_tod;
}

template <>
Transformation<Frame::ECEF, Frame::MOD>::Transformation(const double jd, const EOPData eop) {
    const auto ecef_to_pef = Transformation<Frame::ECEF, Frame::PEF>(jd, eop);
    const auto pef_to_mod = Transformation<Frame::PEF, Frame::MOD>(jd, eop);
    *this = pef_to_mod * ecef_to_pef;
}

template <>
Transformation<Frame::TEME, Frame::MOD>::Transformation(const double jd, const EOPData eop) {
    const auto teme_to_pef = Transformation<Frame::TEME, Frame::PEF>(jd, eop);
    const auto pef_to_mod = Transformation<Frame::PEF, Frame::MOD>(jd, eop);
    *this = pef_to_mod * teme_to_pef;
}

template <>
Transformation<Frame::ECEF, Frame::TOD>::Transformation(const double jd, const EOPData eop) {
    const auto ecef_to_pef = Transformation<Frame::ECEF, Frame::PEF>(jd, eop);
    const auto pef_to_tod = Transformation<Frame::PEF, Frame::TOD>(jd, eop);
    *this = pef_to_tod * ecef_to_pef;
}

template <>
Transformation<Frame::TEME, Frame::TOD>::Transformation(const double jd, const EOPData eop) {
    const auto teme_to_pef = Transformation<Frame::TEME, Frame::PEF>(jd, eop);
    const auto pef_to_tod = Transformation<Frame::PEF, Frame::TOD>(jd, eop);
    *this = pef_to_tod * teme_to_pef;
}

template <>
Transformation<Frame::TEME, Frame::ECEF>::Transformation(const double jd, const EOPData eop) {
    const auto teme_to_pef = Transformation<Frame::TEME, Frame::PEF>(jd, eop);
    const auto pef_to_ecef = Transformation<Frame::ECEF, Frame::PEF>(jd, eop).inverse();
    *this = pef_to_ecef * teme_to_pef;
}

// Inverse Transformations

template <>
Transformation<Frame::J2000, Frame::MOD>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::MOD, Frame::J2000>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::MOD, Frame::TOD>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::TOD, Frame::MOD>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::TOD, Frame::PEF>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::PEF, Frame::TOD>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::PEF, Frame::ECEF>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::ECEF, Frame::PEF>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::PEF, Frame::TEME>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::TEME, Frame::PEF>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::J2000, Frame::TOD>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::TOD, Frame::J2000>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::J2000, Frame::PEF>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::PEF, Frame::J2000>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::J2000, Frame::ECEF>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::ECEF, Frame::J2000>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::J2000, Frame::TEME>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::TEME, Frame::J2000>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::MOD, Frame::PEF>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::PEF, Frame::MOD>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::MOD, Frame::ECEF>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::ECEF, Frame::MOD>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::MOD, Frame::TEME>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::TEME, Frame::MOD>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::TOD, Frame::ECEF>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::ECEF, Frame::TOD>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::TOD, Frame::TEME>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::TEME, Frame::TOD>(jd, eop);
    *this = inverse.inverse();
}

template <>
Transformation<Frame::ECEF, Frame::TEME>::Transformation(const double jd, const EOPData eop) {
    auto const inverse = Transformation<Frame::TEME, Frame::ECEF>(jd, eop);
    *this = inverse.inverse();
}

} // namespace gelocus

#endif // GELOCUS_HPP
