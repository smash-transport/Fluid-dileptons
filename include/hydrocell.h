#pragma once

#include "vector.h"
#include "setup.h"
#include "dilepton.h"

namespace FluidDileptons {

// This assumes a non-interacting hadron gas
double grand_canonical_density(double T, double mu, double mass, double spin, int degeneracy = 1);

// includes n, p, nbar, pbar; but not baryonic resonances
double grand_canonical_nucleon_density(double T, double mu);

// Main interface with the hydro code
class HydroCell {
  public:
    HydroCell(double T, double muB, double QGP_frac,
            const FourVector& pos, const ThreeVector& vel,
            double four_volume = 1) :
              position_(pos),
              landau_vel_(vel),
              T_(T),
              muB_(muB),
              QGP_fraction_(QGP_frac),
              four_volume_(four_volume) {
        nuc_dens_ = grand_canonical_nucleon_density(T_, muB_);
    }
    FourVector position() const { return position_; }
    const DileptonList& dileptons() const;
    void radiate();

  private:
    bool radiated_ = false;
    FourVector position_;
    ThreeVector landau_vel_;
    double T_, muB_, nuc_dens_, QGP_fraction_;
    /*
     * Usually all cells have the same 4-volume, so it could be static,
     * but it's safer to take it from the hydro than from a config.
     * Alternatively, it can be unset and taken care of in the analysis.
     */
    double four_volume_;
    DileptonList dileptons_;
};

} // FluidDileptons
