/**
 * gelocus -- Earth-centered reference frame transformations library
 * "Scientific code that shouldn't feel like it was written by scientists!"
 * https://github.com/gunvirranu/gelocus
 *
 * SPDX-License-Identifier: MIT
 * Copyright (c) 2021 Gunvir Singh Ranu
 */

#ifndef GELOCUS_HPP
#define GELOCUS_HPP

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace gelocus {

using std::size_t;

enum class Frame {
    J2000,
    MOD,
    TOD,
    TEME,
    PEF,
    ECEF
};

template <Frame F>
class Position {
public:
    Vec vec;

    explicit Position(Vec v) : vec(v) {}

    template <Frame To>
    Position<To> transform(double jd) const;

    const auto& operator()(size_t i) const {
        return vec(i);
    }

    auto& operator()(size_t i) {
        return vec(i);
    }
};

template <Frame From, Frame To>
class Transformation {
public:
    Matrix mat;

    explicit Transformation(Matrix mat) : mat(mat) {}

    explicit Transformation(double jd);

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
};

template <Frame F>
template <Frame To>
Position<To> Position<F>::transform(const double jd) const {
    const auto A = Transformation<F, To>(jd, eop);
    return A.apply(*this);
}

// For pseudo-"private" implementation details
namespace detail {

// Common

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
