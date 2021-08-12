#ifndef GELOCUS_HPP
#define GELOCUS_HPP

namespace gelocus {

enum class Frame {
    J2000,
    MOD,
    TOD,
    PEF,
    ECEF
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
