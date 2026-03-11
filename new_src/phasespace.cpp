#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <cassert>
#include <algorithm>
#include <sstream>
#include <memory>

std::mt19937 rng;
std::uniform_real_distribution<double> uniform(0.0, 1.0);

constexpr double alpha_em = 1.0 / 137.0;
constexpr double alpha_em_sqr = alpha_em*alpha_em;
constexpr double CEM = 5.0 / 9.0; // u+d, no strange?
constexpr int Nc = 3;
constexpr double hbarc = 0.197327;
constexpr double m_electron = 0.000511;
constexpr double pi_cube = M_PI * M_PI * M_PI;

// generic fixed-size vector class
template<typename T, size_t N>
class Vector {
  public:
    Vector() { data_.fill(T{}); }
    Vector(std::initializer_list<T> list) {
        assert(list.size() == N && "initializer list must have N elements");
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

    // dot product
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
    FourVector &operator=(const FourVector &other) = default;

    FourVector operator+(const FourVector &rhs) const {
        return FourVector{v_ + rhs.v_};
    }
    FourVector operator-(const FourVector &rhs) const {
        return FourVector{v_ - rhs.v_};
    }
    FourVector operator*(double s) const {
        return FourVector{v_ * s};
    }
    FourVector operator/(double s) const {
        return FourVector{v_ / s};
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
    inline double sqr() const { return (*this) * (*this); }
    inline double abs() const { return std::sqrt(sqr()); }
    
/*
    ThreeVector velocity() const { return threevec() / x0(); }
    double tau() const { return std::sqrt(x0()*x0() - x3()*x3()); }
    double pT() const { return std::sqrt(x1()*x1() + x2()*x2()); }
    double phi() const { return std::atan2(x2(), x1()); }
    double eta() const { return std::atanh(x3() / x0()); }
*/

    FourVector lorentz_boost(const ThreeVector& vel) const {
        const double v_sqr = vel.sqr();
        const double gamma = v_sqr < 1. ? 1. / std::sqrt(1. - v_sqr) : 0;
        const double xprime_0 = gamma * (x0() - threevec() * vel);
        const double constantpart = gamma / (gamma + 1) * (xprime_0 + x0());
        return FourVector(xprime_0, threevec() - vel * constantpart);
    }
};
static const FourVector unit = FourVector{1,0,0,0};

struct HydroCell {
    FourVector position;
    FourVector landau_vel;
    double T, muB, nB, e_dens, nuc_dens;
};

enum class DileptonSource: int {
    rho = 9001,
    omega = 9002,
    phi = 9003,
    multipi = 9004,
    QGP = 9005
};

// Always evaluated in the rest frame of the cell. Position is not needed, 
// but could be useful for cuts (and SMASH integration)
class Dilepton {
  public:
    Dilepton() = default;
    // make constructor that receives cell
    Dilepton(double mass, double abs_mom, const DileptonSource& source) : mass_(mass) {

    }
    Dilepton(const FourVector& pos, const FourVector& mom,
            const DileptonSource& source) : position_(pos),
                                     momentum_(mom),
                                     emission_time_(pos.x0()),
                                     mass_(mom.abs()),
                                     source_(source) {}
    double m() { return mass_; }
    double q() const { return momentum_.threevec().abs(); }
    // return momentum of each lepton
    std::pair<FourVector, FourVector> single_lepton() const;
  private:
    double emission_time_ = 0, mass_ = 0;
    FourVector position_, momentum_;
    DileptonSource source_;
    std::shared_ptr<HydroCell> cell_ = nullptr;
};

std::pair<FourVector, FourVector> Dilepton::single_lepton() {
    FourVector p1{}, p2{};
    return std::make_pair<FourVector, FourVector>{p1,p2};
}

double dilepton_phasespace(double mee) {
    const double ratio_sqr = m_electron * m_electron / (mee * mee);
    if (1.0 - 4.0 * ratio_sqr < 0.0)
        return 0.0;
    return std::sqrt(1.0 - 4.0 * ratio_sqr) * (1.0 + 2.0 * ratio_sqr);
}

double fhat_2 (double m, double q, double T) {
    double fhat;
    const double q0 = std::sqrt(q * q + m * m);
    if (q==0) {
        fhat = std::tanh(q0/(4*T));
    } else {
        const double x_p = std::exp(-(q0 + q) / (2.0 * T));
        const double x_m = std::exp(-(q0 - q) / (2.0 * T));
        fhat = 1.0 + 2.0 * T / q * std::log((1.0 + x_p) / (1.0 + x_m));
    }
    return fhat;
}

double Qlat_total(double m, double q, double T, bool FF=true) {
    const double T_sqr = T*T;
    const double m_sqr = m*m;
    const double q0 = std::sqrt(q * q + m_sqr);
    double form_factor = 1;
    if (FF)
        form_factor = 4 * T_sqr / (4 * T_sqr + m_sqr);
    const double alpha_s = 6.0 * M_PI / (28.0 * std::log(T / 0.022));
    constexpr double K = 2.;
    const double Q_transv = 2 * M_PI * alpha_s * (T_sqr / m_sqr) * K * form_factor *
                                   std::log(1 + 2.912 * q0 / (4 * M_PI * alpha_s * T));
    const double Q_total = (2 + m_sqr / q0 / q0) * Q_transv / 3.;
    return Q_total;
}

// used for cross-check
double ImD_pQGP(double m, double q, double T) {
    return CEM * Nc * m * m * fhat_2(m, q, T) / (12 * M_PI);
}

double ImD_lat(double m, double q, double T) {
    const double fhat = fhat_2(m, q, T);
    const double Q_total = Qlat_total(m, q, T, false);
    // multiply by 4pi to match with fig 3 of arXiv:1304.2309
    const double ImD_total = (CEM * Nc / (12 * M_PI)) * m*m * (fhat + Q_total);
    return ImD_total;
}

// used for cross-check but is different from the paper
double ImD_rho_vacuum(double m, [[maybe_unused]] double q=0, [[maybe_unused]] double T=0) {
    constexpr double M0 = 0.778;
    constexpr double M0_sqr = M0*M0;
    constexpr double G0 = 0.149;
    constexpr double m_thr = 0.139;
    const double m_sqr = m*m;
    const double G = m >= m_thr ? G0 * std::sqrt((1-m_thr*m_thr/m_sqr) / (1-m_thr*m_thr/M0_sqr)) : 0;
    const double gamma = std::sqrt(M0_sqr*(M0_sqr+G*G));
    const double BW = m_sqr*M0*G*gamma / (std::pow(m_sqr-M0_sqr, 2) + m_sqr*G*G) / std::sqrt(M0_sqr+gamma);
    return BW/(2*M_PI);
}

double dR_dMd3q_without_ImD(double mee, double q, double T){
    if (T == 0.0)
        return 0.0;
    const double dil_ps = dilepton_phasespace(mee);
    if (dil_ps == 0.0)
        return 0.0;

    const double q0 = std::sqrt(q * q + mee*mee);
    const double f_BE = 1.0 / (std::exp(q0 / T) - 1.0);
    const double dR_d4q_coef = alpha_em_sqr * dil_ps  * f_BE / (pi_cube*mee*mee);
    // dR/dMd³q = M/q⁰ dR/d⁴q
    return (mee/q0) * dR_d4q_coef;
}

double dR_dMd3q_lat(double mee, double q, double T){
    const double coef = dR_dMd3q_without_ImD(mee, q, T);
    const double ImD = ImD_lat(mee, q, T);
    return coef * ImD;
}


double ImD_rho_medium(double mee, double q, double T) {
    std::cout << "Not implemented yet!"
    return 0;
}

void radiator(const HydroCell& cell) {

}

/*
double rho_lat(double q0, double T){
    constexpr double kappa = 0.0465;
    const double G = T*2.235;
    const double cte_BW = 1.098*G;

    const double BW = cte_BW*q0*(G/2)/(q0*q0+G*G/4.);
    const double cont = Nc*(1+kappa)*q0*q0*std::tanh(q0/4/T)/(2*M_PI);
    return (BW + cont)/4/M_PI;
}
*/

int main (int argc, char *argv[]) {
    const double T = 0.17;
    const double dq = 0.01;

    std::ofstream fout("dRdM_rhovac.txt");
    std::vector<std::string> buffer;
    buffer.reserve(300);
    for (double mee=0.01; mee<1.6; mee+=0.02) {
        double y = 0;
        for (double q=0;q<20;q+=dq) {
            y += q*q*dR_dMd3q_lat(mee, q, T);
        }
        std::ostringstream oss;
        // dR/dM² = (2M)⁻¹dR/dM and given in fm⁻⁴
        oss << mee << " " <<  4*M_PI*dq*y/std::pow(hbarc,4)/(2*mee)<< "\n";
        buffer.push_back(oss.str());
    }
    for (const auto& line : buffer) {
        fout << line;
    }
    fout.close();
    return 0;
}