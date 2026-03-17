#pragma once

#include "vector.h"
#include "setup.h"

namespace FluidDileptons {

// The enum value is used in the output
enum Source {
    rho = 9001,
    omega = 9002,
    phi = 9003,
    multipi = 9004,
    QGP = 9005
};
static const std::vector<Source> all_sources {
    rho, omega, phi, multipi, QGP};

static int invalidArgument(std::string message = "") {
    throw std::invalid_argument(message);
    return -1;
}

// Always evaluated in the (Landau?) rest frame of the cell. Position is not needed
// explicitly, as it can be taken from the cell information
class Dilepton {
  public:
    Dilepton(double mass, double q, const Source source,
             const HydroCell* cell, double weight=0) :
        mass_(mass),
        q_(q),
        source_(source),
        cell_(cell),
        weight_((weight >= 0 && weight <= 1) ? weight : invalidArgument("weight = " + std::to_string(weight))) {
            ThreeVector vel = sample_isotropic_3D_versor()*q;
            momentum_ = FourVector{std::sqrt(mass*mass + q*q), vel};
        }

    double m() const { return mass_; }
    double q() const { return q_; }
    double weight() const { return weight_; }
    Source source() const { return source_; }
    FourVector momentum() const { return momentum_; }

    void set_weight(double w) { weight_ = w; }
    void set_momentum(const FourVector& mom) { momentum_ = mom; }

    // return 4-momentum of each lepton
    std::pair<FourVector, FourVector> single_lepton() const;
  private:
    double mass_, q_;
    double weight_;
    Source source_;
    FourVector momentum_;
    const HydroCell* cell_ = nullptr;
};

template <typename T>
using build_vector_ = std::vector<T, std::allocator<T>>;
using DileptonList = build_vector_<Dilepton>;

namespace Rates {
    struct Parameters {
        double m, q, T, nuc_dens;
    };

    double dilepton_phasespace(double mee);
    // used for cross-check
    double ImD_pQGP(const Parameters& par);
    double ImD_latQGP(const Parameters& par);
    // used for cross-check but leads to a different result from the paper
    double ImD_rho_vacuum(const Parameters& par);
    double ImD_rho_medium(const Parameters& par);
    double dR_dMd3q_without_ImD(const Parameters& par);
    double dR_dMd3q_QGP(const Parameters& par);
    double choose_source(const Source s,const Parameters& par);


    // Interpolation function for pion tadpole diagram
    double tadpole_interpolate(double temperature);

    // Imaginary part of vector correlator (4-pion contribution)
    double Im_PiV4pi(double s);

    // Imaginary part of 3-pion axial-vector correlator (a1 decay)
    double Im_PiA3pi(double s);

    // Imaginary part of 5-pion axial-vector correlator
    double Im_PiA5pi(double s);

    // Imaginary part of mixed vector correlator with pion effects
    double ImD_multipi(const Parameters& par);


} // namespace Rates

} // FluidDileptons
