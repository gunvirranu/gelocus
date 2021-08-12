#ifndef GELOCUS_HPP
#define GELOCUS_HPP

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

} // namespace gelocus

#endif // GELOCUS_HPP
