#pragma once

#include "vector.h"
#include "setup.h"
#include "dilepton.h"

namespace FluidDileptons {

// This assumes a non-interacting hadron gas
static double grand_canonical_density(double T, double mu, double mass, double spin, int degeneracy);

// includes n, p, nbar, pbar; but not baryonic resonances
static double grand_canonical_nucleon_density(double T, double mu);

// Main interface with the hydro code
class HydroCell {
  public:
    HydroCell(double T, double muB, double QGP_frac,
            const FourVector& pos, const ThreeVector& vel,
            double four_volume = 1) :
                T_(T),
                muB_(muB),
                QGP_fraction_(QGP_frac),
                position_(pos),
                landau_vel_(vel),
                four_volume_(four_volume) {
        nuc_dens_ = grand_canonical_nucleon_density(T_, muB_);
    }
    FourVector position() const { return position_; }
    DileptonList dileptons() const;
    void radiate();

  private:
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


namespace Grid{
    inline std::vector<double> masses, qs;

    static void set_vector(std::vector<double>& vec,
                           double min, double max, double step) {
        vec.clear();
        for (size_t i = 0; min + i*step <= max; ++i)
            vec.push_back(min + i*step);
    }
    inline void set_masses(double min, double max, double step) {
        set_vector(masses, min, max, step);
    }
    inline void set_masses(std::vector<double> vec) {
        masses = vec;
    }
    inline void set_qs(double min, double max, double step) {
        set_vector(qs, min, max, step);
    }
    inline void set_qs(std::vector<double> vec) {
        qs = vec;
    }
}


static double grand_canonical_density(double T, double mu, double mass, double spin, int degeneracy=1) {
    const int stat = static_cast<int>(round(2. * spin)) & 1 ? -1 : 1;
    double density = 0;
    for (int i = 1; i < 11; i++) {
        density += std::pow(stat, i + 1) *
                    std::cyl_bessel_k(2, i * mass / T) *
                    (std::exp(i * mu / T) + std::exp(-i * mu / T)) / i;
                    // p,n                         pbar,nbar
    }
    density *= degeneracy * (2 * spin + 1) * mass * mass * T /
                    (2 * M_PI * M_PI * std::pow(hbarc, -3));
    return density; // in fm⁻³
}

static double grand_canonical_nucleon_density(double T, double mu) {
    constexpr double nuc_mass = 0.938;
    constexpr int g = 2; // isospin
    constexpr double nuc_J = 0.5;
    return grand_canonical_density(T, mu, nuc_mass, nuc_J, g);
}

} // FluidDileptons
