#pragma once

#include <cmath>
#include <random>

namespace FluidDileptons {

constexpr double alpha_em = 1.0 / 137.0;
constexpr double alpha_em_sqr = alpha_em*alpha_em;
constexpr double CEM = 5.0 / 9.0;
constexpr int Nc = 3;
constexpr double hbarc = 0.197327;
constexpr double mass_electron = 0.000511;
constexpr double pi_cube = M_PI * M_PI * M_PI;
constexpr double mass_pion = 0.139;

extern int N_oversample;
extern std::ranlux24_base random;
extern std::uniform_real_distribution<> uniform;

// forward declarations
class Dilepton;
class HydroCell;
void output_cell_to_file(const std::string& filepath, const HydroCell& cell);

void notImplemented();
bool setup(const std::string& filepath);

}
