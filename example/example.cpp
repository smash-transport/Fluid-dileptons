#include <iostream>

#include "fluiddileptons.h"

int main (int argc, char *argv[]) {
    const double T = 0.14, muB=0.2, lambdaQGP=0.5, fourVolume=1; // fake cell parameters
    const FluidDileptons::FourVector pos = {10,0,0,3}; // fake cell 4-position
    FluidDileptons::ThreeVector vel = {0,0.2,0}; // fake cell Landau velocity

    /*
     * The radiation can be configured either inline or by config, or a combination of both,
     * which allows for dynamically adapting the parameters.
     */
    bool setup = FluidDileptons::setup("../dilepton_config.txt");

    // For example:
    FluidDileptons::AcceptanceCutter::set_yrap_range(-1,1);

    /*
     * The main interface with the hydro code is the HydroCell class, which takes the cell
     * parameters and radiates dileptons, which can be accessed by the hydro or output to file.
     */
    FluidDileptons::HydroCell cell{T, muB, lambdaQGP, pos, vel, fourVolume};
    cell.radiate();
    FluidDileptons::output_cell_to_file("example_cell.dat", cell);
    FluidDileptons::output_spectra_to_file("example_spectra.dat");

    return 0;
}