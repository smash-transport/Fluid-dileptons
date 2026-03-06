#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include <sstream>

std::mt19937 rng;
std::uniform_real_distribution<double> uniform(0.0, 1.0);

constexpr double alpha_em = 1.0 / 137.0;
constexpr double alpha_em_sqr = alpha_em*alpha_em;
constexpr double CEM = 5.0 / 9.0; // u+d, no strange?
constexpr int Nc = 3;
constexpr double hbarc = 0.197327;
constexpr double m_electron = 0.000511;
constexpr double pi_cube = M_PI * M_PI * M_PI;

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
    const double Q_transv = 2 * M_PI * alpha_s * (T_sqr / m_sqr) * 2 * form_factor *
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

double ImD_rawa_rho()

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