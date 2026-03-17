#pragma once

#include <cmath>
#include <random>

namespace FluidDileptons {

constexpr double alpha_em = 1.0 / 137.0;
constexpr double alpha_em_sqr = alpha_em*alpha_em;
constexpr double CEM = 5.0 / 9.0; // u+d, no strange?
constexpr int Nc = 3;
constexpr double hbarc = 0.197327;
constexpr double m_electron = 0.000511;
constexpr double pi_cube = M_PI * M_PI * M_PI;
constexpr double pion_mass = 0.139;  // pion mass in GeV
constexpr double kaon_mass = 0.494;  // kaon mass in GeV

extern int N_oversample;

inline std::random_device rd;
inline std::ranlux24_base random(rd());
inline std::uniform_real_distribution<> uniform(0, 1);

// forward declarations
class Dilepton;
class HydroCell;
void output_cell_to_file(const std::string& filepath, const HydroCell& cell);

void notImplemented();

bool setup(const std::string& filepath);

}
