#ifndef GELOCUS_HPP
#define GELOCUS_HPP

#include <cstddef>
#include <iostream>
#include <stdexcept>

#ifndef GELOCUS_VEC3
#error Must define `GELOCUS_VEC3` type
#endif

namespace gelocus {

using std::size_t;
using Vec = GELOCUS_VEC3;

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

std::ostream& operator<<(std::ostream &os, Frame f) {
    switch (f) {
        case Frame::J2000: return os << "J2000";
        case Frame::MOD: return os << "MOD";
        case Frame::TOD: return os << "TOD";
        case Frame::PEF: return os << "PEF";
        case Frame::ECEF: return os << "ECEF";
    };
    throw std::runtime_error("Invalid frame type");
}

} // namespace gelocus

#endif // GELOCUS_HPP
