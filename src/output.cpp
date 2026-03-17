#include <fstream>
#include <sstream>

#include "dilepton.h"
#include "hydrocell.h"
#include "output.h"
#include "vector.h"

namespace FluidDileptons {

void output_cell_to_file(const std::string& filepath, const HydroCell& cell) {
    std::ofstream fout(filepath);
    std::vector<std::string> buffer;
    buffer.reserve(300);
    for (const auto& dil: cell.dileptons() ) {
        std::ostringstream oss;
//        oss << cell.position().x0() << " " << dil.m() << " " <<  dil.momentum() << " " << dil.weight() << " " << dil.source() << "\n";
        oss << dil.m() << " " << dil.weight() << " " << dil.source() << "\n";
        buffer.push_back(oss.str());
    }
    for (const auto& line : buffer) {
        fout << line;
    }
}

}
