#include <iostream>
#include <cmath>

#include "dilepton.h"
#include "setup.h"

namespace FluidDileptons {

std::pair<FourVector, FourVector> Dilepton::single_lepton() const {
    FourVector p1{}, p2{};
    notImplemented();
    return std::make_pair(p1, p2);
}

namespace Rates {

    static double dilepton_phasespace(double mee) {
        const double ratio_sqr = mass_electron * mass_electron / (mee * mee);
        if (1.0 - 4.0 * ratio_sqr < 0.0)
            return 0.0;
        return std::sqrt(1.0 - 4.0 * ratio_sqr) * (1.0 + 2.0 * ratio_sqr);
    }

    static double fhat_2 (const Parameters& par) {
        const double m = par.m, q = par.q, T = par.T;
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

    static double Qlat_total(const Parameters& par, bool FF=true) {
        const double m = par.m, q = par.q, T = par.T;
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
    double ImD_pQGP(const Parameters& par) {
        const double m = par.m;
        return CEM * Nc * m * m * fhat_2(par) / (12 * M_PI);
    }

    double ImD_latQGP(const Parameters& par) {
        const double m = par.m, q = par.q, T = par.T;
        // Calculation not valid below transition region (even there it's doubtful)
        if (T < 0.12)
            return 0;
        const double fhat = fhat_2(par);
        const double Q_total = Qlat_total(par, false);
        // multiply by 4pi to match with fig 3 of arXiv:1304.2309
        const double ImD_total = (CEM * Nc / (12 * M_PI)) * m*m * (fhat + Q_total);
        return ImD_total;
    }

    static double breit_wigner(double m, double m0, double G0, double m_thr) {
        const double m_sqr = m*m, m0_sqr = m0*m0;
        const double G = (m >= m_thr) ? G0 * std::sqrt((1 - m_thr*m_thr/m_sqr) / (1 - m_thr*m_thr/m0_sqr)) : 0;
        return m_sqr * m0 * G / (std::pow(m_sqr - m0_sqr, 2) + m_sqr * G * G) / std::sqrt(m0_sqr+G*G);
    }

    double ImD_rho_vacuum(const Parameters& par) {
        constexpr double M0 = 0.778;
        constexpr double G0 = 0.149;
        constexpr double m_thr = 2*mass_electron;
        return breit_wigner(par.m, M0, G0, m_thr);
    }
    double ImD_omega_vacuum(const Parameters& par) {
        constexpr double M0 = 0.783;
        constexpr double G0 = 0.00849;
        constexpr double m_thr = 2*mass_electron;
        return breit_wigner(par.m, M0, G0, m_thr);
    }
    double ImD_phi_vacuum(const Parameters& par) {
        constexpr double M0 = 1.019;
        constexpr double G0 = 0.00426;
        constexpr double m_thr = 2*mass_electron;
        return breit_wigner(par.m, M0, G0, m_thr);
    }

    double ImD_rho_medium(const Parameters& par) {
        std::cout << "Not implemented yet!";
        return 0;
    }

    double dR_dMd3q_without_ImD(const Parameters& par){
        const double mee = par.m, q = par.q, T = par.T;
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

    // In principle we could compute coef outside the function, but the performance is not
    // much worse, and this way we can parallelize the source loop together with m and q
    double choose_source(const Source s, const Parameters& par) {
        const double frac = par.QGP_fraction;
        const double coef = dR_dMd3q_without_ImD(par);
        double rate = 0;
/*
        static bool warned_vacuum = false;
        if (!warned_vacuum && (s == Source::rho || s == Source::omega || s == Source::phi)) {
            std::cerr << "Using vacuum spectral function for vector mesons, "
                      << "in-medium rates not yet implemented.\n";
            warned_vacuum = true;
        }
*/

        if (s == Source::QGP) {
            rate = frac * coef * ImD_latQGP(par);
        } else if (frac < 1) {
            if (s == Source::multipi) {
                rate =  (1-frac) * coef * ImD_multipi(par);
            } else if (s == Source::rho) {
                rate =  (1-frac) * coef * ImD_rho_vacuum(par);
            } else if (s == Source::omega) {
                rate =  (1-frac) * coef * ImD_omega_vacuum(par);
            } else if (s == Source::phi) {
                rate =  (1-frac) * coef * ImD_phi_vacuum(par);
            } else {
                notImplemented();
            }
        }
        return rate;
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
} // namespace Rates

}
