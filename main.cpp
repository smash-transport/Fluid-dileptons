#include <iostream>
#include <fstream>
#include <sstream>

#include "fluiddileptons.h"

int main (int argc, char *argv[]) {
    const double T = 0.15, muB=0.1, lambdaQGP=0.4, fourVolume=0.5; // fake cell parameters
    FluidDileptons::FourVector pos = {10,0,0,3}; // fake cell 4-position
    FluidDileptons::ThreeVector vel = {0,0,0.7}; // fake cell Landau velocity

    /*
     * The radiation can be configured either inline or by config, or a combination of both,
     * which allows for dynamically adapting the parameters.
     */

    FluidDileptons::AcceptanceCutter::set_x_range(-0.1,0.1);
    FluidDileptons::setup("../dil_config.txt");

    FluidDileptons::HydroCell cell{T, muB, lambdaQGP, pos, vel, fourVolume};
    cell.radiate();
    FluidDileptons::output_cell_to_file("example.dat", cell);

    return 0;
}