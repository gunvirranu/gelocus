#include <iostream>
#include <array>
#include <stdexcept>

class Vec3 {
public:
    double x, y, z;
    Vec3() = default;
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
    const double& operator()(size_t i) const {
        switch (i) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: throw std::out_of_range("Index is out of range");
        }
    }
    double& operator()(size_t i) {
        switch (i) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: throw std::out_of_range("Index is out of range");
        }
    }
    friend std::ostream& operator<<(std::ostream &os, const Vec3 &v) {
        return os << "Vec3 { " << v.x << ", " << v.y << ", " << v.z << " }";
    }
};

class Mat3 {
public:
    std::array<std::array<double, 3>, 3> m;
    Mat3() = default;
    Mat3 transposed() const {
        Mat3 M;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                M(i, j) = (*this)(j, i);
            }
        }
        return M;
    }
    Mat3 operator*(const Mat3 &rhs) const {
        Mat3 out;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                out(i, j) = 0;
                for (size_t k = 0; k < 3; ++k) {
                    out(i, j) += (*this)(i, k) * rhs(k, j);
                }
            }
        }
        return out;
    }
    Vec3 operator*(const Vec3 &v) const {
        const auto &A = *this;
        Vec3 b;
        for (size_t i = 0; i < 3; ++i) {
            b(i) = A(i, 0) * v(0) + A(i, 1) * v(1) + A(i, 2) * v(2);
        }
        return b;
    }
    const double& operator()(size_t i, size_t j) const {
        return m[i][j];
    }
    double& operator()(size_t i, size_t j) {
        return m[i][j];
    }
    friend std::ostream& operator<<(std::ostream &os, const Mat3 &M) {
        return os << "Mat3 { " <<
            M(0, 0) << ", " << M(0, 1) << ", " << M(0, 2) << " ; " <<
            M(1, 0) << ", " << M(1, 1) << ", " << M(1, 2) << " ; " <<
            M(2, 0) << ", " << M(2, 1) << ", " << M(2, 2) << " }";
    }
};

#define GELOCUS_VEC3 Vec3
#define GELOCUS_MATRIX3 Mat3
#define GELOCUS_MATRIX3_TRANSPOSED transposed
#include "gelocus.hpp"

using gelocus::Frame;
template <Frame F> using Pos = gelocus::Position<F>;
template <Frame From, Frame To> using Trans = gelocus::Transformation<From, To>;

int main() {
    const auto v = Vec3(1, 2, 3);
    const auto p = Pos<Frame::TOD>(v);
    const auto jd = 2450000.0;
    const auto p2 = p.transform<Frame::MOD>(jd);

    std::cout << p << "\n";
    std::cout << p2 << "\n";
    std::cout << "Hello World!\n";
}
