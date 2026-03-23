#include <cmath>
#include <cassert>

#include "acceptance.h"
#include "hydrocell.h"
#include "spectrum.h"

namespace FluidDileptons {

double grand_canonical_density(double T, double mu, double mass, double spin, int degeneracy) {
    const int stat = static_cast<int>(round(2. * spin)) & 1 ? -1 : 1;
    double density = 0;
    for (int i = 1; i < 11; i++) {
        density += std::pow(stat, i + 1) *
                    std::cyl_bessel_k(2, i * mass / T) *
                    (std::exp(i * mu / T) + std::exp(-i * mu / T)) / i;
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
        dileptons_.reserve(N_oversample * MQGrid::masses.size() * MQGrid::qs.size() * all_sources.size());
    }

    for (double m: MQGrid::masses) {
        if (m <= 0.0001) {
            continue;
        }
        for (double q: MQGrid::qs) {
            const Rates::Parameters dilepton_par{m, q, T_, nuc_dens_};
            std::vector<double> source_weight;

            // Precompute weights before oversampling
            for (Source s: all_sources) {
                double weight = four_volume_ * Rates::choose_source(s, dilepton_par);
                weight *= (s == Source::QGP) ? QGP_fraction_ : (1 - QGP_fraction_);
                source_weight.push_back(weight / N_oversample); // can be 0;
            }

            for (int n_sample = 0; n_sample < N_oversample; ++n_sample) {
                for (size_t is = 0; is < all_sources.size(); ++is) {
                    Source s = all_sources[is];
                    const double weight = source_weight[is];
                    if (weight == 0)
                        continue;
                    auto dil = Dilepton(m, q, s, this, weight);
                    // boost to lab frame and do acceptance
                    dil.set_momentum(dil.momentum().lorentz_boost(-landau_vel_));
                    if (!AcceptanceCutter::in_momentum_range(dil.momentum()))
                        continue;
                    // TODO: acceptance cut on single lepton if needed
                    if (OutputMode::dilepton) {
                        dileptons_.push_back(dil);
                    }
                    if (OutputMode::spectra) {
                        Spectra::fill(dil);
                    }
                }
            }
        }
    }
    radiated_ = true;
}

} // FluidDileptons