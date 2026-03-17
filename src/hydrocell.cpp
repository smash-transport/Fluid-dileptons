#include <iostream>

#include "acceptance.h"
#include "hydrocell.h"
#include "dilepton.h"

namespace FluidDileptons {

int N_oversample = 1;

DileptonList HydroCell::dileptons() const {
    return dileptons_;
}

void HydroCell::radiate() {
    if (!AcceptanceCutter::in_spatial_range(position_))
        return;

    dileptons_.reserve(N_oversample * Grid::masses.size() * Grid::qs.size() * all_sources.size());
    for (double m: Grid::masses) {
        if (m <= 0.0001) {
            continue;
        }
        for (double q: Grid::qs) {
            const Rates::Parameters dilepton_par{m, q, T_, nuc_dens_};
            for (int n_sample = 0; n_sample < N_oversample; ++n_sample) {
                for (Source s: all_sources) {
                    double weight = four_volume_ * Rates::choose_source(s, dilepton_par);
                    weight *= (s == Source::QGP) ? QGP_fraction_ : (1 - QGP_fraction_);
                    if (weight == 0)
                        continue;
                    auto dil = Dilepton(m, q, s, this, weight / N_oversample);
                    // boost to lab frame and do acceptance
                    dil.set_momentum(dil.momentum().lorentz_boost(-landau_vel_));
                    if (!AcceptanceCutter::in_momentum_range(dil.momentum()))
                        continue;
                    // TODO: acceptance cut on single lepton if needed
                    dileptons_.push_back(dil);
                }
            }
        }
    }
}

}