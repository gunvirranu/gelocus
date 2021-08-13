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
    PEF,
    ECEF
};

std::ostream& operator<<(std::ostream &os, Frame f);

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

    friend std::ostream& operator<<(std::ostream &os, const Position &p) {
        return os << "Position<" << F << "> { " << p(0) << ", " << p(1) << ", " << p(2) << " }";
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

    friend std::ostream& operator<<(std::ostream &os, const Transformation &A) {
        return os << "Transformation<" << From << ", " << To << "> { " << A.mat << " }";
    }
};

std::ostream& operator<<(std::ostream &os, Frame f) {
    switch (f) {
        case Frame::J2000: return os << "J2000";
        case Frame::MOD: return os << "MOD";
        case Frame::TOD: return os << "TOD";
        case Frame::PEF: return os << "PEF";
        case Frame::ECEF: return os << "ECEF";
    };
    throw std::invalid_argument("Invalid frame type");
}

template <Frame F>
template <Frame To>
Position<To> Position<F>::transform(const double jd) const {
    const auto A = Transformation<F, To>(jd);
    return A.apply(*this);
}

// Util Functions

constexpr double rad_to_deg(const double rad) {
    return rad * DEG_PER_RAD;
}

constexpr double deg_to_rad(const double deg) {
    return deg * RAD_PER_DEG;
}

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

// Implemented "Basic" Transformations

template <Frame From, Frame To>
Transformation<From, To>::Transformation(const double jd) {
    static_assert(From == To, "This transformation is not implemented");
    (void) jd;  // To ignore unused value
    // Use  identity matrix since `From == To`
    mat(0, 0) = 1; mat(0, 1) = 0; mat(0, 2) = 0;
    mat(1, 0) = 0; mat(1, 1) = 1; mat(1, 2) = 0;
    mat(2, 0) = 0; mat(2, 1) = 0; mat(2, 2) = 1;
}

} // namespace gelocus

#endif // GELOCUS_HPP
