#include <cmath>
#include <iterator>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "acceptance.h"
#include "hydrocell.h"
#include "spectrum.h"

namespace FluidDileptons {

double grand_canonical_density(double T, double mu, double mass, double spin, int degeneracy) {
    const int stat = static_cast<int>(round(2. * spin)) & 1 ? -1 : 1;
    double density = 0;
    if (mass / T > 10) {
        // Approximation for low T
        density = degeneracy * (2 * spin + 1) *
                   std::pow(mass * T / (2 * M_PI), 1.5) *
                    std::exp(-mass / T) * std::cosh(mu / T);
        return density;
    }
    for (int i = 1; i < 11; i++) {
        density += std::pow(stat, i + 1) *
                    std::cyl_bessel_k(2, i * mass / T) *
                    2 * std::cosh(mu / T) / i;
    }
    density *= degeneracy * (2 * spin + 1) * mass * mass * T /
                    (2 * M_PI * M_PI * std::pow(hbarc, -3));
    return density;
}

double grand_canonical_nucleon_density(double T, double mu) {
    constexpr double nuc_mass = 0.938;
    constexpr int g = 2;
    constexpr double nuc_J = 0.5;
    return grand_canonical_density(T, mu, nuc_mass, nuc_J, g);
}

const DileptonList& HydroCell::dileptons() const {
    return dileptons_;
}

void HydroCell::radiate() {
    if(radiated_) {
        throw std::runtime_error("Radiation already done for this cell, cannot radiate again.");
    }
    if (!AcceptanceCutter::in_spatial_range(position_)) {
        return;
    }

    if (OutputMode::dilepton) {
        dileptons_.reserve(N_oversample * all_sources.size() *
                             MQGrid::masses.size() * MQGrid::mom_abs.size());
    }

    const Rates::VectorMesonRateGrid vector_meson_rates =
        Rates::precompute_vector_meson_rate_grid(T_, nuc_dens_, QGP_fraction_);

    #pragma omp parallel for collapse(3)
    for (std::size_t source_index = 0; source_index < all_sources.size(); ++source_index) {
        for (std::size_t mass_index = 0; mass_index < MQGrid::masses.size(); ++mass_index) {
            for (std::size_t q_index = 0; q_index < MQGrid::mom_abs.size(); ++q_index) {
                const Source s = all_sources[source_index];
                const double m = MQGrid::masses[mass_index];
                const double q = MQGrid::mom_abs[q_index];
                if (m <= 0.0001) {
                    continue;
                }
                double weight = 0.0;
                if (s == Source::rho || s == Source::omega || s == Source::phi) {
                    weight = four_volume_ * vector_meson_rates.at(s, mass_index, q_index);
                } else {
                    const Rates::Parameters dilepton_par{m, q, T_, nuc_dens_, QGP_fraction_};
                    weight = four_volume_ * Rates::choose_source(s, dilepton_par);
                }
                if (weight == 0.0) {
                    continue;
                }
                for (int n_sample = 0; n_sample < N_oversample; ++n_sample) {
                    auto dil = Dilepton(m, q, s, this, weight);
                    dil.set_momentum(dil.momentum().lorentz_boost(-landau_vel_));
                    if (!AcceptanceCutter::in_momentum_range(dil.momentum()))
                        continue;
                    // TODO: acceptance cut on single lepton if needed
                    if (OutputMode::dilepton) {
                        #pragma omp critical
                        dileptons_.push_back(dil);
                    }
                    if (OutputMode::spectra)
                        Spectra::add(s, m, q, weight);
                }
            }
        }
    }
    radiated_ = true;
}

} // namespace FluidDileptons