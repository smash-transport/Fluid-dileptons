#include <fstream>
#include <iomanip>
#include <sstream>

#include "dilepton.h"
#include "hydrocell.h"
#include "output.h"
#include "setup.h"
#include "spectrum.h"
#include "vector.h"

namespace FluidDileptons {

void output_cell_to_file(const std::string& filepath, const HydroCell& cell) {
    if (!OutputMode::dilepton) {
        throw std::invalid_argument("Dilepton output mode not enabled. Add it to the config to output dileptons.\n");
    }
    std::ofstream fout(filepath);
    std::vector<std::string> buffer;
    buffer.reserve(300);
    for (const auto& dil: cell.dileptons() ) {
        std::ostringstream oss;
//        oss << cell.position().x0() << " " << dil.m() << " " <<  dil.momentum() << " " << dil.weight() << " " << dil.source() << "\n";
        oss << dil.m() << " " << dil.q() << " " << dil.weight() << " " << dil.source() << "\n";
        buffer.push_back(oss.str());
    }
    for (const auto& line : buffer) {
        fout << line;
    }
}

void output_spectra_to_file(const std::string& filepath) {
    if (!OutputMode::spectra) {
        throw std::invalid_argument("Spectra output mode not enabled. Add it to the config to output spectra.\n");
    }
    std::ofstream fout(filepath);
    std::ostringstream buffer;

    buffer << "# Columns: q ";
    for (double q : MQGrid::qs) {
        buffer << q << " ";
    }
    buffer << "\n# Rows: mass ";
    for (double m : MQGrid::masses) {
        buffer << m << " ";
    }
    buffer << "\n" << std::scientific << std::setprecision(6);

    for (const Source source : all_sources) {
        const Spectrum& spectrum = Spectra::get(source);
        std::ostringstream block;
        block << "# source " << source << "\n";

        // One row per q-bin; each row has masses.size() entries.
        // This is faster (somehow) than using indices directly
        for (double q: MQGrid::qs) {
            for (double m: MQGrid::masses) {
                if (m > MQGrid::masses.front()) {
                    block << ' ';
                }
                block << spectrum.weight_at(m, q);
            }
            block << '\n';
        }
        buffer << block.str();
    }
    fout << buffer.str();
}

void output_source_spectrum_to_file(const std::string& filepath, Source source) {
    std::ofstream fout(filepath);
    std::ostringstream buffer;

    buffer << "# Columns: q_grid ";
    for (double q : MQGrid::qs) {
        buffer << q << " ";
    }
    buffer << "\n# Rows: mass_grid ";
    for (double m : MQGrid::masses) {
        buffer << m << " ";
    }
    buffer << "\n" << std::scientific << std::setprecision(6);

    const Spectrum& spectrum = Spectra::get(source);
    buffer << "# source " << source << "\n";

        // One row per q-bin; each row has masses.size() entries.
        // This is faster (somehow) than using indices directly
        for (double q: MQGrid::qs) {
            for (double m: MQGrid::masses) {
            if (m > MQGrid::masses.front()) {
                buffer << ' ';
            }
            buffer << spectrum.weight_at(m, q);
        }
        buffer << '\n';
    }

    fout << buffer.str();
}


} // namespace FluidDileptons
