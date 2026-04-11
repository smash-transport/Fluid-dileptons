#include <iostream>
#include <cmath>
#include <filesystem>

#include "dilepton.h"
#include "setup.h"
#include "spectrum.h"
#include "VMtables.h"

namespace FluidDileptons {

std::pair<FourVector, FourVector> Dilepton::single_lepton() const {
    FourVector p1{}, p2{};
    notImplemented();
    return std::make_pair(p1, p2);
}

namespace Rates {

    namespace {

        std::size_t mq_index(std::size_t mass_index, std::size_t q_index) {
            return mass_index * MQGrid::mom_abs.size() + q_index;
        }

        const VectorMesonTable& rho_omega_table() {
            static const VectorMesonTable table = load_rho_omega_table(
                VectorMesonTables::rho_omega_path);
            return table;
        }

        const VectorMesonTable& phi_table() {
            static const VectorMesonTable table = load_phi_table(
                VectorMesonTables::phi_path);
            return table;
        }

        const VectorMesonGridCache& rho_omega_grid_cache() {
            static const VectorMesonGridCache cache =
                rho_omega_table().make_grid_cache(MQGrid::mom_abs, MQGrid::masses);
            return cache;
        }

        const VectorMesonGridCache& phi_grid_cache() {
            static const VectorMesonGridCache cache =
                phi_table().make_grid_cache(MQGrid::mom_abs, MQGrid::masses);
            return cache;
        }

        std::size_t rho_channel() {
            static const std::size_t channel = rho_omega_table().channel_index("rho");
            return channel;
        }

        std::size_t omega_channel() {
            static const std::size_t channel = rho_omega_table().channel_index("omega");
            return channel;
        }

        std::size_t phi_channel() {
            static const std::size_t channel = phi_table().channel_index("phi");
            return channel;
        }

        double interpolate_vector_meson(const VectorMesonTable& table,
                                        std::size_t channel,
                                        const Parameters& par,
                                        double pion_chemical_potential = 0.0,
                                        double kaon_chemical_potential = 0.0) {
            return table.interpolate(channel, par.T, par.nuc_dens, pion_chemical_potential,
                                     kaon_chemical_potential, par.q, par.m);
        }

    } // namespace

    double VectorMesonRateGrid::at(Source source,
                                   std::size_t mass_index,
                                   std::size_t q_index) const {
        const std::size_t index = mass_index * q_bins + q_index;
        if (source == Source::rho) {
            return rho[index];
        }
        if (source == Source::omega) {
            return omega[index];
        }
        if (source == Source::phi) {
            return phi[index];
        }
        return 0.0;
    }

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
        return interpolate_vector_meson(rho_omega_table(), rho_channel(), par);
    }

    double ImD_omega_medium(const Parameters& par) {
        return interpolate_vector_meson(rho_omega_table(), omega_channel(), par);
    }

    double ImD_phi_medium(const Parameters& par) {
        return interpolate_vector_meson(phi_table(), phi_channel(), par);
    }

    VectorMesonRateGrid precompute_vector_meson_rate_grid(double temperature,
                                                          double nucleon_density,
                                                          double qgp_fraction) {
        VectorMesonRateGrid grid;
        grid.q_bins = MQGrid::mom_abs.size();

        const std::size_t total_size = MQGrid::masses.size() * MQGrid::mom_abs.size();
        grid.rho.assign(total_size, 0.0);
        grid.omega.assign(total_size, 0.0);
        grid.phi.assign(total_size, 0.0);

        const double hadronic_fraction = 1.0 - qgp_fraction;
        if (hadronic_fraction <= 0.0) {
            return grid;
        }

        if (SpectralFunctionMode::vacuum_rho_omega && SpectralFunctionMode::vacuum_phi) {
            for (std::size_t mass_index = 0; mass_index < MQGrid::masses.size(); ++mass_index) {
                const double mass = MQGrid::masses[mass_index];
                if (mass <= 0.0001) {
                    continue;
                }
                for (std::size_t q_index = 0; q_index < MQGrid::mom_abs.size(); ++q_index) {
                    const double momentum = MQGrid::mom_abs[q_index];
                    const Parameters par{mass, momentum, temperature, nucleon_density, qgp_fraction};
                    const double coef = hadronic_fraction * dR_dMd3q_without_ImD(par);
                    if (coef == 0.0) {
                        continue;
                    }
                    const std::size_t index = mq_index(mass_index, q_index);
                    grid.rho[index] = coef * ImD_rho_vacuum(par);
                    grid.omega[index] = coef * ImD_omega_vacuum(par);
                    grid.phi[index] = coef * ImD_phi_vacuum(par);
                }
            }
            return grid;
        }

        VectorMesonCellContext rho_omega_context{};
        VectorMesonCellContext phi_context{};
        if (!SpectralFunctionMode::vacuum_rho_omega) {
            rho_omega_context = rho_omega_table().make_cell_context(temperature, nucleon_density, 0.0, 0.0);
        }
        if (!SpectralFunctionMode::vacuum_phi) {
            phi_context = phi_table().make_cell_context(temperature, nucleon_density, 0.0, 0.0);
        }

        for (std::size_t mass_index = 0; mass_index < MQGrid::masses.size(); ++mass_index) {
            const double mass = MQGrid::masses[mass_index];
            if (mass <= 0.0001) {
                continue;
            }
            const AxisBracket& rho_mass_bracket = rho_omega_grid_cache().masses[mass_index];
            const AxisBracket& phi_mass_bracket = phi_grid_cache().masses[mass_index];

            for (std::size_t q_index = 0; q_index < MQGrid::mom_abs.size(); ++q_index) {
                const double momentum = MQGrid::mom_abs[q_index];
                const Parameters par{mass, momentum, temperature, nucleon_density, qgp_fraction};
                const double coef = hadronic_fraction * dR_dMd3q_without_ImD(par);
                if (coef == 0.0) {
                    continue;
                }

                const std::size_t index = mq_index(mass_index, q_index);
                if (SpectralFunctionMode::vacuum_rho_omega) {
                    grid.rho[index] = coef * ImD_rho_vacuum(par);
                    grid.omega[index] = coef * ImD_omega_vacuum(par);
                } else {
                    grid.rho[index] = coef * rho_omega_table().interpolate_fast(
                        rho_channel(), rho_omega_context,
                        rho_omega_grid_cache().momenta[q_index], rho_mass_bracket);
                    grid.omega[index] = coef * rho_omega_table().interpolate_fast(
                        omega_channel(), rho_omega_context,
                        rho_omega_grid_cache().momenta[q_index], rho_mass_bracket);
                }

                if (SpectralFunctionMode::vacuum_phi) {
                    grid.phi[index] = coef * ImD_phi_vacuum(par);
                } else {
                    grid.phi[index] = coef * phi_table().interpolate_fast(
                        phi_channel(), phi_context,
                        phi_grid_cache().momenta[q_index], phi_mass_bracket);
                }
            }
        }

        return grid;
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
        if (s == Source::QGP) {
            rate = frac * coef * ImD_latQGP(par);
        } else if (frac < 1) {
            if (s == Source::multipi) {
                rate =  (1-frac) * coef * ImD_multipi(par);
            } else if (s == Source::rho) {
                rate =  (1-frac) * coef *
                        (SpectralFunctionMode::vacuum_rho_omega ?
                            ImD_rho_vacuum(par) : ImD_rho_medium(par));
            } else if (s == Source::omega) {
                rate =  (1-frac) * coef *
                        (SpectralFunctionMode::vacuum_rho_omega ?
                            ImD_omega_vacuum(par) : ImD_omega_medium(par));
            } else if (s == Source::phi) {
                rate =  (1-frac) * coef *
                        (SpectralFunctionMode::vacuum_phi ?
                            ImD_phi_vacuum(par) : ImD_phi_medium(par));
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
