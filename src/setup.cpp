#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "acceptance.h"
#include "setup.h"
#include "hydrocell.h"

namespace FluidDileptons {

int N_oversample = 1;

std::random_device rd;
std::ranlux24_base random(rd());
std::uniform_real_distribution<> uniform(0, 1);

void notImplemented() {
    std::cerr << "Not implemented yet!" << std::endl;
}

// Trims whitespace from strings
static std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\r\n");
    size_t end = str.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

// Type alias for AcceptanceCutter and Grid function pointers
using RangeSetterFunc = void (*)(double, double);
using RangeStepSetterFunc = void (*)(double, double, double);

static void parse_range(const std::string& key, const std::string& value,
                        RangeSetterFunc setter) {
    double min, max;
    std::stringstream ss(value);
    if (ss >> min >> max) {
        setter(min, max);
        std::cout << "Set " << key << ": [" << min << ", " << max << "]\n";
    } else {
        std::cerr << "Could not parse " << key << ": " << value << "\n";
    }
}

static void parse_range_step(const std::string& key, const std::string& value,
                             RangeStepSetterFunc setter) {
    double min, max, step;
    std::stringstream ss(value);
    if (ss >> min >> max >> step) {
        setter(min, max, step);
        std::cout << "Set " << key << ": [" << min << ", " << max << "] with step " << step << "\n";
    } else {
        std::cerr << "Could not parse " << key << ": " << value << "\n";
    }
}

bool setup(const std::string& filepath) {
    default_setup();
    std::ifstream config_file(filepath);
    std::ostringstream oss_err;

    if (!config_file.is_open()) {
        throw std::invalid_argument("Could not open config file.\n");
    }

    std::string line;
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
            parse_range_step(key, value, Grid::set_masses);
        } else if (key == "abs_momenta") {
            parse_range_step(key, value, Grid::set_qs);
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
                std::cout << "Set N_oversample: " << N_oversample << "\n";
            } else {
                std::cerr << "Could not parse N_oversample: " << value << "\n";
            }
        } else {
            throw std::invalid_argument(key + " is not a configuration key for dileptons.\n");
        }
    }

    config_file.close();
    // todo: validate values?

    std::cout << "Configuration loaded successfully\n";
    return true;
}

void default_setup() {
    constexpr double large_cut = 100;

    Grid::set_masses(0.0, 2, 0.01);
    Grid::set_qs(0.0, 3, 0.1);

    AcceptanceCutter::set_x_range(-large_cut, large_cut);
    AcceptanceCutter::set_y_range(-large_cut, large_cut);
    AcceptanceCutter::set_eta_range(-large_cut, large_cut);
    AcceptanceCutter::set_pT_range(0, large_cut);
    AcceptanceCutter::set_yrap_range(-large_cut, large_cut);

    N_oversample = 1;
}

}