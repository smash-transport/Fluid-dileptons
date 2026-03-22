#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include "setup.h"

namespace FluidDileptons {

// generic fixed-size vector class
template<typename T, size_t N>
class Vector {
  public:
    Vector() { data_.fill(T{}); }
    Vector(std::initializer_list<T> list) {
        std::copy(list.begin(), list.end(), data_.begin());
    }
    Vector(const Vector &other) = default;
    Vector &operator=(const Vector &other) = default;

    T operator[](size_t i) const { return data_[i]; }
    T &operator[](size_t i) { return data_[i]; }

    Vector operator+(const Vector &rhs) const {
        Vector res;
        for (size_t i = 0; i < N; ++i) res[i] = data_[i] + rhs[i];
        return res;
    }
    Vector operator-(const Vector &rhs) const {
        Vector res;
        for (size_t i = 0; i < N; ++i) res[i] = data_[i] - rhs[i];
        return res;
    }
    Vector operator-() const {
        Vector res;
        for (size_t i = 0; i < N; ++i) res[i] = -data_[i];
        return res;
    }
    Vector operator*(T scalar) const {
        Vector res;
        for (size_t i = 0; i < N; ++i) res[i] = data_[i] * scalar;
        return res;
    }
    Vector operator/(T scalar) const {
        Vector res;
        for (size_t i = 0; i < N; ++i) res[i] = data_[i] / scalar;
        return res;
    }
    Vector &operator+=(const Vector &rhs) {
        for (size_t i = 0; i < N; ++i) data_[i] += rhs[i];
        return *this;
    }
    Vector &operator-=(const Vector &rhs) {
        for (size_t i = 0; i < N; ++i) data_[i] -= rhs[i];
        return *this;
    }
    Vector &operator*=(T scalar) {
        for (size_t i = 0; i < N; ++i) data_[i] *= scalar;
        return *this;
    }
    Vector &operator/=(T scalar) {
        for (size_t i = 0; i < N; ++i) data_[i] /= scalar;
        return *this;
    }
    T operator*(const Vector &rhs) const {
        T sum{};
        for (size_t i = 0; i < N; ++i) sum += data_[i] * rhs[i];
        return sum;
    }
    T sqr() const {
        return (*this) * (*this);
    }
    T abs() const {
        return std::sqrt(this->sqr());
    }

  protected:
    std::array<T, N> data_;
};

using ThreeVector = Vector<double,3>;

class FourVector {
  private:
    Vector<double, 4> v_;
  public:
    FourVector() = default;
    FourVector(std::initializer_list<double> list) : v_{list} {}
    FourVector(double x0, ThreeVector xvec) : v_{ {x0, xvec[0], xvec[1], xvec[2]} } {}
    explicit FourVector(const Vector<double, 4>& vec) : v_{vec} {}
    FourVector &operator=(const FourVector &other) = default;

    FourVector operator+(const FourVector &rhs) const {
        return FourVector(v_ + rhs.v_);
    }
    FourVector operator-(const FourVector &rhs) const {
        return FourVector(v_ - rhs.v_);
    }
    FourVector operator-() const {
        return FourVector(-v_);
    }
    FourVector operator*(double s) const {
        return FourVector(v_ * s);
    }
    FourVector operator/(double s) const {
        return FourVector(v_ / s);
    }
    FourVector &operator+=(const FourVector &rhs) {
        v_ += rhs.v_;
        return *this;
    }
    FourVector &operator-=(const FourVector &rhs) {
        v_ -= rhs.v_;
        return *this;
    }
    FourVector &operator*=(double scalar) {
        v_ *= scalar;
        return *this;
    }
    FourVector &operator/=(double scalar) {
        v_ /= scalar;
        return *this;
    }

    double x0() const { return v_[0]; }
    double x1() const { return v_[1]; }
    double x2() const { return v_[2]; }
    double x3() const { return v_[3]; }
    ThreeVector threevec() const { return ThreeVector{v_[1],v_[2],v_[3]}; }

    double operator*(const FourVector &rhs) const {
        return this->x0() * rhs.x0() - this->threevec() * rhs.threevec();
    }
    double sqr() const { return (*this) * (*this); }
    double abs() const { return std::sqrt(sqr()); }
    double eta() const { return std::atanh(x3()/x0()); }
    double xT() const { return std::sqrt(x1()*x1() + x2()*x2()); }

    FourVector lorentz_boost(const ThreeVector& vel) const {
        const double v_sqr = vel.sqr();
        if (v_sqr >= 1) {
            throw std::runtime_error("Velocity for lorentz boost larger than one.");
        }
        const double gamma = 1. / std::sqrt(1. - v_sqr);
        const double xprime_0 = gamma * (x0() - threevec() * vel);
        const double constantpart = gamma / (gamma + 1) * (xprime_0 + x0());
        return FourVector(xprime_0, threevec() - vel * constantpart);
    }
};

inline std::ostream& operator<<(std::ostream& os, const ThreeVector& v) {
    os << v[0] << " " << v[1] << " " << v[2];
    return os;
}
inline std::ostream& operator<<(std::ostream& os, const FourVector& v) {
    os << v.x0() << " " << v.x1() << " " << v.x2() << " " << v.x3();
    return os;
}

// Marsaglia algorithm
inline ThreeVector sample_isotropic_3D_versor() {
    static thread_local std::ranlux24_base thread_random(std::random_device{}());
    static thread_local std::uniform_real_distribution<> thread_uniform(0.0, 1.0);
    double x1, x2, s;
    do {
        x1 = 2 * thread_uniform(thread_random) - 1;
        x2 = 2 * thread_uniform(thread_random) - 1;
        s = x1*x1 + x2*x2;
    } while (s >= 1);

    double factor = 2*sqrt(1 - s);

    double x = x1 * factor;
    double y = x2 * factor;
    double z = 1 - 2*s;
    return ThreeVector{x,y,z};
}

} // FluidDileptons
