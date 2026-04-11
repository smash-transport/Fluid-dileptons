#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>

#include "acceptance.h"
#include "hydrocell.h"
#include "setup.h"
#include "spectrum.h"

namespace FluidDileptons {

bool OutputMode::spectra = true;
bool OutputMode::dilepton = true;
bool SpectralFunctionMode::vacuum_rho_omega = false;
bool SpectralFunctionMode::vacuum_phi = false;
std::string VectorMesonTables::rho_omega_path;
std::string VectorMesonTables::phi_path;
int N_oversample = 1;

void notImplemented() {
    std::cerr << "Not implemented yet!" << std::endl;
}

// Trims whitespace from strings
static std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\r\n");
    size_t end = str.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

// Type alias for AcceptanceCutter and MQGrid function pointers
using RangeSetterFunc = void (*)(double, double);
using RangeStepSetterFunc = void (*)(double, double, double, bool);

static void parse_range(const std::string& key, const std::string& value,
                        RangeSetterFunc setter) {
    double min, max;
    std::stringstream ss(value);
    if (ss >> min >> max) {
        setter(min, max);
        std::cout << key << ": [" << min << ", " << max << "]\n";
    } else {
        std::cerr << "Could not parse " << key << ": " << value << "\n";
    }
}

static void parse_range_step(const std::string& key, const std::string& value,
                             RangeStepSetterFunc setter) {
    double min, max, step;
    std::stringstream ss(value);
    if (ss >> min >> max >> step) {
        setter(min, max, step, true);
        std::cout << key << ": [" << min << ", " << max << "] with step " << step << "\n";
    } else {
        std::cerr << "Could not parse " << key << ": " << value << "\n";
    }
}

static void suppress_output_mode(const std::string& mode) {
    if (mode == "spectra") {
        OutputMode::spectra = false;
        std::cout << "Disabled spectra output mode\n";
    } else if (mode == "dilepton") {
        OutputMode::dilepton = false;
        std::cout << "Disabled dilepton output mode\n";
    } else if (mode == "none") {
        std::cout << "All output modes enabled\n";
    } else {
        std::cerr << "Unknown output mode: " << mode << "\n";
    }
}

static void set_vector_meson_file_path(const std::string& key,
                                       const std::string& value,
                                       const std::filesystem::path& config_path,
                                       bool& vacuum_flag,
                                       std::string& target_path,
                                       std::ostringstream& oss_err) {
    if (value == "vacuum") {
        vacuum_flag = true;
        std::cout << key << ": vacuum\n";
        return;
    }

    std::filesystem::path table_path = value;
    if (table_path.empty()) {
        oss_err << key << " cannot be empty\n";
        return;
    }
    if (table_path.is_relative()) {
        table_path = config_path.parent_path() / table_path;
    }
    if (!std::filesystem::exists(table_path)) {
        oss_err << key << " does not exist: " << table_path.string() << "\n";
        return;
    }

    vacuum_flag = false;
    target_path = std::filesystem::weakly_canonical(table_path).string();
    std::cout << key << ": " << target_path << "\n";
}

static std::filesystem::path canonical_or_self(const std::string& filepath) {
    const std::filesystem::path path = filepath;
    std::error_code error;
    const auto canonical = std::filesystem::weakly_canonical(path, error);
    if (error) {
        return path;
    }
    return canonical;
}

bool setup(const std::string& filepath) {
    default_setup();
    std::ifstream config_file(filepath);
    std::ostringstream oss_err;
    const std::filesystem::path config_path = canonical_or_self(filepath);

    if (!config_file.is_open()) {
        throw std::invalid_argument("Could not open configuration file for fluid dileptons.\n");
    }

    std::string line;
    std::cout << "======== Reading fluid-dileptons configuration =========\n";
    while (std::getline(config_file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Parse key-value pair
        size_t delim_pos = line.find(':');
        if (delim_pos == std::string::npos) {
            oss_err << "'" << line << "' does not contain ':' separator\n";
        }

        std::string key = trim(line.substr(0, delim_pos));
        std::string value = trim(line.substr(delim_pos + 1));
        std::replace(value.begin(), value.end(), ',', ' '); // remove comma if needed (config looks nicer)

        if (key == "masses") {
            parse_range_step(key, value, MQGrid::set_masses);
        } else if (key == "abs_momenta") {
            parse_range_step(key, value, MQGrid::set_mom_abs);
        } else if (key == "x_range") {
            parse_range(key, value, AcceptanceCutter::set_x_range);
        } else if (key == "y_range") {
            parse_range(key, value, AcceptanceCutter::set_y_range);
        } else if (key == "eta_range") {
            parse_range(key, value, AcceptanceCutter::set_eta_range);
        } else if (key == "pt_range" || key == "pT_range") {
            parse_range("pt_range", value, AcceptanceCutter::set_pT_range);
        } else if (key == "yrap_range") {
            parse_range("yrap_range", value, AcceptanceCutter::set_yrap_range);
        } else if (key == "oversample") {
            // todo: make a generic parse_value?
            std::stringstream ss(value);
            int N;
            if (ss >> N) {
                N_oversample = N;
                std::cout << "N_oversample: " << N_oversample << "\n";
            } else {
                oss_err << "Could not parse N_oversample: " << value << "\n";
            }
        } else if (key == "suppress_output") {
            std::stringstream ss(value);
            std::string mode;
            while (ss >> mode) {
                suppress_output_mode(mode);
            }
        } else if (key == "rho_omega_table_path") {
            set_vector_meson_file_path(key, value, config_path,
                                       SpectralFunctionMode::vacuum_rho_omega,
                                       VectorMesonTables::rho_omega_path,
                                       oss_err);
        } else if (key == "phi_table_path") {
            set_vector_meson_file_path(key, value, config_path,
                                       SpectralFunctionMode::vacuum_phi,
                                       VectorMesonTables::phi_path,
                                       oss_err);
        } else {
            oss_err << key << " is not a configuration key for dileptons.\n";
        }
    }
    if (!oss_err.str().empty()) {
        std::cerr << "Errors while parsing fluid-dileptons configuration:\n" << oss_err.str();
        return false;
    }
    // Make this a config key?
    Spectra::initialize();

    config_file.close();
    // todo: validate values?

    std::cout << "======== Fluid-dileptons configuration read successfully ========\n";
    return true;
}

bool default_setup() {
    constexpr double large_cut = 100;

    MQGrid::set_masses(0.0, 2, 0.01);
    MQGrid::set_mom_abs(0.0, 3, 0.05);

    AcceptanceCutter::set_x_range(-large_cut, large_cut);
    AcceptanceCutter::set_y_range(-large_cut, large_cut);
    AcceptanceCutter::set_eta_range(-large_cut, large_cut);
    AcceptanceCutter::set_pT_range(0, large_cut);
    AcceptanceCutter::set_yrap_range(-large_cut, large_cut);

    OutputMode::spectra = true;
    OutputMode::dilepton = true;
    SpectralFunctionMode::vacuum_rho_omega = false;
    SpectralFunctionMode::vacuum_phi = false;
    VectorMesonTables::rho_omega_path = (std::filesystem::path(FLUID_DILEPTONS_SOURCE_DIR) /
                                         "old" / "diltables" / "rawa" /
                                         "ImDrho-or-f.dat").string();
    VectorMesonTables::phi_path = (std::filesystem::path(FLUID_DILEPTONS_SOURCE_DIR) /
                                   "old" / "diltables" / "rawa" /
                                   "ImDrho-phi.dat").string();

    N_oversample = 1;
    return true;
}

} // namespace FluidDileptons
