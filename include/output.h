#pragma once

#include "setup.h"

namespace FluidDileptons {

void output_cell_to_file(const std::string& filepath, const HydroCell& cell);
void output_spectra_to_file(const std::string& filepath);
void output_source_spectrum_to_file(const std::string& filepath, Source source);

}
