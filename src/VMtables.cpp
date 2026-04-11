#include "VMtables.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace FluidDileptons::Rates {

namespace {
constexpr double kTemperatureStart = 0.050;
constexpr double kTemperatureStep = 0.010;
constexpr double kDensityStep = 0.5;
constexpr double kMomentumStep = 0.25;
constexpr double kMassStep = 0.01;
constexpr double kMeVToGeV = 1e-3;
constexpr double kAxisTolerance = 1e-9;

std::size_t checked_flat_index(const std::array<std::size_t, 6>& shape,
                               std::size_t temperature_index,
                               std::size_t density_index,
                               std::size_t pion_index,
                               std::size_t kaon_index,
                               std::size_t momentum_index,
                               std::size_t mass_index) {
    if (temperature_index >= shape[0] || density_index >= shape[1] ||
        pion_index >= shape[2] || kaon_index >= shape[3] ||
        momentum_index >= shape[4] || mass_index >= shape[5]) {
        throw std::out_of_range("Rapp table index out of range");
    }

    return (((((temperature_index * shape[1] + density_index) * shape[2] + pion_index) *
              shape[3] + kaon_index) * shape[4] + momentum_index) *
            shape[5]) + mass_index;
}

void require_open(const std::ifstream& input, const std::string& filepath) {
    if (!input.is_open()) {
        throw std::runtime_error("Could not open Rapp table file: " + filepath);
    }
}

void read_block_header(std::ifstream& input, const std::string& filepath) {
    std::string header;
    input >> std::ws;
    std::getline(input, header);
    if (!input) {
        throw std::runtime_error("Unexpected end of file while reading block header from " + filepath);
    }
}

template <typename... Args>
void read_values(std::ifstream& input, const std::string& filepath, Args&... args) {
    if (!(input >> ... >> args)) {
        throw std::runtime_error("Unexpected end of file or malformed row while reading " + filepath);
    }
}

const std::vector<double>& get_axis(const VectorMesonTable& table, std::size_t axis_index) {
    switch (axis_index) {
        case 0: return table.axes.temperatures;
        case 1: return table.axes.baryon_densities;
        case 2: return table.axes.pion_chemical_potentials;
        case 3: return table.axes.kaon_chemical_potentials;
        case 4: return table.axes.momenta;
        case 5: return table.axes.masses;
        default: throw std::out_of_range("Vector meson axis index out of range");
    }
}

AxisInfo infer_axis_info(const std::vector<double>& axis, double expected_step = 0.0) {
    if (axis.empty()) {
        throw std::logic_error("Cannot infer metadata for an empty axis");
    }

    if (axis.size() == 1) {
        return {true, axis.front(), 0.0};
    }

    auto is_regular_with_step = [&](double step) {
        if (step <= 0.0) {
            return false;
        }
        for (std::size_t i = 1; i < axis.size(); ++i) {
            if (std::abs((axis[i] - axis[i - 1]) - step) >
                kAxisTolerance * std::max(1.0, std::abs(step))) {
                return false;
            }
        }
        return true;
    };

    if (is_regular_with_step(expected_step)) {
        return {true, axis.front(), expected_step};
    }

    const double inferred_step = axis[1] - axis[0];
    if (is_regular_with_step(inferred_step)) {
        return {true, axis.front(), inferred_step};
    }

    return {false, axis.front(), 0.0};
}

AxisBracket bracket_for_axis(const std::vector<double>& axis,
                             const AxisInfo& axis_info,
                             double value) {
    if (axis.empty()) {
        throw std::logic_error("Cannot interpolate on an empty axis");
    }

    if (axis.size() == 1 || value <= axis.front()) {
        return {0, 0, 0.0};
    }

    if (value >= axis.back()) {
        const std::size_t last = axis.size() - 1;
        return {last, last, 0.0};
    }

    if (axis_info.regular && axis_info.step > 0.0) {
        const double coordinate = (value - axis_info.min) / axis_info.step;
        std::size_t lower_index = static_cast<std::size_t>(std::floor(coordinate));
        if (lower_index >= axis.size() - 1) {
            lower_index = axis.size() - 2;
        }
        const std::size_t upper_index = lower_index + 1;
        const double lower_value = axis[lower_index];
        const double upper_value = axis[upper_index];
        const double delta = upper_value - lower_value;
        if (delta <= 0.0) {
            return {lower_index, lower_index, 0.0};
        }
        return {lower_index, upper_index, (value - lower_value) / delta};
    }

    const auto upper = std::upper_bound(axis.begin(), axis.end(), value);
    const std::size_t upper_index = static_cast<std::size_t>(upper - axis.begin());
    const std::size_t lower_index = upper_index - 1;
    const double lower_value = axis[lower_index];
    const double upper_value = axis[upper_index];
    const double delta = upper_value - lower_value;

    if (delta <= 0.0) {
        return {lower_index, lower_index, 0.0};
    }

    return {lower_index, upper_index, (value - lower_value) / delta};
}

double blend_weight(const AxisBracket& bracket, bool use_upper) {
    if (!use_upper || bracket.lower == bracket.upper) {
        return 1.0 - bracket.weight_upper;
    }
    return bracket.weight_upper;
}

std::size_t blend_index(const AxisBracket& bracket, bool use_upper) {
    return (use_upper && bracket.lower != bracket.upper) ? bracket.upper : bracket.lower;
}

} // namespace

void VectorMesonTable::finalize_metadata() {
    strides[5] = 1;
    for (int axis = 4; axis >= 0; --axis) {
        strides[axis] = strides[axis + 1] * layout.shape[axis + 1];
    }

    axis_info[0] = infer_axis_info(axes.temperatures, kTemperatureStep);
    axis_info[1] = infer_axis_info(axes.baryon_densities, kDensityStep);
    axis_info[2] = infer_axis_info(axes.pion_chemical_potentials);
    axis_info[3] = infer_axis_info(axes.kaon_chemical_potentials);
    axis_info[4] = infer_axis_info(axes.momenta, kMomentumStep);
    axis_info[5] = infer_axis_info(axes.masses, kMassStep);
}

std::size_t VectorMesonTable::channel_index(const std::string& channel_name) const {
    const auto it = std::find(channel_names.begin(), channel_names.end(), channel_name);
    if (it == channel_names.end()) {
        throw std::out_of_range("Unknown vector-meson channel: " + channel_name);
    }
    return static_cast<std::size_t>(it - channel_names.begin());
}

AxisBracket VectorMesonTable::bracket(std::size_t axis_index, double value) const {
    return bracket_for_axis(get_axis(*this, axis_index), axis_info[axis_index], value);
}

VectorMesonCellContext VectorMesonTable::make_cell_context(double temperature,
                                                           double density,
                                                           double pion_chemical_potential,
                                                           double kaon_chemical_potential) const {
    return {{
        bracket(0, temperature),
        bracket(1, density),
        bracket(2, pion_chemical_potential),
        bracket(3, kaon_chemical_potential)
    }};
}

VectorMesonGridCache VectorMesonTable::make_grid_cache(const std::vector<double>& momenta,
                                                       const std::vector<double>& masses) const {
    VectorMesonGridCache cache;
    cache.momenta.reserve(momenta.size());
    cache.masses.reserve(masses.size());
    for (double momentum : momenta) {
        cache.momenta.push_back(bracket(4, momentum));
    }
    for (double mass : masses) {
        cache.masses.push_back(bracket(5, mass));
    }
    return cache;
}

std::size_t VectorMesonTable::flat_index(std::size_t temperature_index,
                                         std::size_t density_index,
                                         std::size_t pion_index,
                                         std::size_t kaon_index,
                                         std::size_t momentum_index,
                                         std::size_t mass_index) const {
    return checked_flat_index(shape(), temperature_index, density_index, pion_index,
                              kaon_index, momentum_index, mass_index);
}

double VectorMesonTable::at(std::size_t channel,
                            std::size_t temperature_index,
                            std::size_t density_index,
                            std::size_t pion_index,
                            std::size_t kaon_index,
                            std::size_t momentum_index,
                            std::size_t mass_index) const {
    if (channel >= channel_values.size()) {
        throw std::out_of_range("Vector meson channel index out of range");
    }
    return channel_values[channel][flat_index(temperature_index, density_index, pion_index,
                                              kaon_index, momentum_index, mass_index)];
}

double VectorMesonTable::at(const std::string& channel_name,
                            std::size_t temperature_index,
                            std::size_t density_index,
                            std::size_t pion_index,
                            std::size_t kaon_index,
                            std::size_t momentum_index,
                            std::size_t mass_index) const {
    return at(channel_index(channel_name), temperature_index, density_index, pion_index,
              kaon_index, momentum_index, mass_index);
}

double VectorMesonTable::interpolate(std::size_t channel,
                                     double temperature,
                                     double density,
                                     double pion_chemical_potential,
                                     double kaon_chemical_potential,
                                     double momentum,
                                     double mass) const {
    const VectorMesonCellContext cell_context = make_cell_context(
        temperature, density, pion_chemical_potential, kaon_chemical_potential);
    return interpolate_fast(channel, cell_context, bracket(4, momentum), bracket(5, mass));
}

double VectorMesonTable::interpolate_fast(std::size_t channel,
                                          const VectorMesonCellContext& cell_context,
                                          const AxisBracket& momentum_bracket,
                                          const AxisBracket& mass_bracket) const {
    if (channel >= channel_values.size()) {
        throw std::out_of_range("Vector meson channel index out of range");
    }

    const std::array<AxisBracket, 6> brackets{
        cell_context.fixed_axes[0],
        cell_context.fixed_axes[1],
        cell_context.fixed_axes[2],
        cell_context.fixed_axes[3],
        momentum_bracket,
        mass_bracket
    };

    const double* values = channel_values[channel].data();
    std::array<std::size_t, 6> lower_indices{};
    std::array<std::size_t, 6> active_axes{};
    std::size_t active_count = 0;
    std::size_t base_index = 0;

    for (std::size_t axis = 0; axis < brackets.size(); ++axis) {
        lower_indices[axis] = brackets[axis].lower;
        base_index += lower_indices[axis] * strides[axis];
        if (!brackets[axis].fixed()) {
            active_axes[active_count++] = axis;
        }
    }

    if (active_count == 0) {
        return values[base_index];
    }

    double value = 0.0;
    for (unsigned int mask = 0; mask < (1u << active_count); ++mask) {
        double weight = 1.0;
        std::size_t index = base_index;

        for (std::size_t bit = 0; bit < active_count; ++bit) {
            const std::size_t axis = active_axes[bit];
            const bool use_upper = (mask & (1u << bit)) != 0;
            const AxisBracket& bracket_value = brackets[axis];
            weight *= blend_weight(bracket_value, use_upper);
            if (use_upper) {
                index += (bracket_value.upper - bracket_value.lower) * strides[axis];
            }
        }

        if (weight == 0.0) {
            continue;
        }

        value += weight * values[index];
    }

    return value;
}

double VectorMesonTable::interpolate(const std::string& channel_name,
                                     double temperature,
                                     double density,
                                     double pion_chemical_potential,
                                     double kaon_chemical_potential,
                                     double momentum,
                                     double mass) const {
    return interpolate(channel_index(channel_name), temperature, density,
                       pion_chemical_potential, kaon_chemical_potential,
                       momentum, mass);
}

VectorMesonTable load_rho_omega_table(const std::string& filepath) {
    std::ifstream input(filepath);
    require_open(input, filepath);

    VectorMesonTable table{rho_omega_layout};
    table.channel_names = {"rho", "omega"};
    table.axes.temperatures.resize(table.temperature_bins());
    table.axes.baryon_densities.resize(table.density_bins());
    table.axes.pion_chemical_potentials.resize(table.pion_chemical_potential_bins());
    table.axes.kaon_chemical_potentials.resize(table.kaon_chemical_potential_bins(), 0.0);
    table.axes.momenta.resize(table.momentum_bins());
    table.axes.masses.resize(table.mass_bins());
    table.channel_values.resize(table.channel_names.size(), std::vector<double>(table.total_size()));

    for (std::size_t i = 0; i < table.temperature_bins(); ++i) {
        read_block_header(input, filepath);
        table.axes.temperatures[i] = kTemperatureStart + static_cast<double>(i) * kTemperatureStep;

        for (std::size_t j = 0; j < table.density_bins(); ++j) {
            for (std::size_t k = 0; k < table.pion_chemical_potential_bins(); ++k) {
                for (std::size_t l = 0; l < table.kaon_chemical_potential_bins(); ++l) {
                    for (std::size_t n = 0; n < table.momentum_bins(); ++n) {
                        for (std::size_t m = 0; m < table.mass_bins(); ++m) {
                            double rho_density = 0.0;
                            double pion_mu = 0.0;
                            double mass = 0.0;
                            double momentum = 0.0;
                            double rho_value = 0.0;
                            double omega_value = 0.0;
                            read_values(input, filepath, rho_density, pion_mu, mass, momentum,
                                        rho_value, omega_value);

                            if (i == 0 && k == 0 && l == 0 && m == 0 && n == 0) {
                                table.axes.baryon_densities[j] = rho_density;
                            }
                            if (i == 0 && j == 0 && l == 0 && m == 0 && n == 0) {
                                table.axes.pion_chemical_potentials[k] = pion_mu * kMeVToGeV;
                            }
                            if (i == 0 && j == 0 && k == 0 && l == 0 && m == 0) {
                                table.axes.momenta[n] = momentum * kMeVToGeV;
                            }
                            if (i == 0 && j == 0 && k == 0 && l == 0 && n == 0) {
                                table.axes.masses[m] = mass;
                            }

                            const std::size_t index = table.flat_index(i, j, k, l, n, m);
                            table.channel_values[0][index] = rho_value;
                            table.channel_values[1][index] = omega_value;
                        }
                    }
                }
            }
        }
    }
    table.finalize_metadata();
    std::cout << "Loaded rho and omega spectral functions from " << filepath << "\n";
    return table;
}

VectorMesonTable load_phi_table(const std::string& filepath) {
    std::ifstream input(filepath);
    require_open(input, filepath);

    VectorMesonTable table{phi_layout};
    table.channel_names = {"phi"};
    table.axes.temperatures.resize(table.temperature_bins());
    table.axes.baryon_densities.resize(table.density_bins());
    table.axes.pion_chemical_potentials.resize(table.pion_chemical_potential_bins());
    table.axes.kaon_chemical_potentials.resize(table.kaon_chemical_potential_bins());
    table.axes.momenta.resize(table.momentum_bins());
    table.axes.masses.resize(table.mass_bins());
    table.channel_values.resize(table.channel_names.size(), std::vector<double>(table.total_size()));

    for (std::size_t i = 0; i < table.temperature_bins(); ++i) {
        read_block_header(input, filepath);
        table.axes.temperatures[i] = kTemperatureStart + static_cast<double>(i) * kTemperatureStep;

        for (std::size_t j = 0; j < table.density_bins(); ++j) {
            for (std::size_t k = 0; k < table.pion_chemical_potential_bins(); ++k) {
                for (std::size_t l = 0; l < table.kaon_chemical_potential_bins(); ++l) {
                    for (std::size_t n = 0; n < table.momentum_bins(); ++n) {
                        for (std::size_t m = 0; m < table.mass_bins(); ++m) {
                            double rho_density = 0.0;
                            double pion_mu = 0.0;
                            double kaon_mu = 0.0;
                            double mass = 0.0;
                            double momentum = 0.0;
                            double phi_value = 0.0;
                            read_values(input, filepath, rho_density, pion_mu, kaon_mu, mass,
                                        momentum, phi_value);

                            if (i == 0 && k == 0 && l == 0 && m == 0 && n == 0) {
                                table.axes.baryon_densities[j] = rho_density;
                            }
                            if (i == 0 && j == 0 && l == 0 && m == 0 && n == 0) {
                                table.axes.pion_chemical_potentials[k] = pion_mu * kMeVToGeV;
                            }
                            if (i == 0 && j == 0 && k == 0 && m == 0 && n == 0) {
                                table.axes.kaon_chemical_potentials[l] = kaon_mu * kMeVToGeV;
                            }
                            if (i == 0 && j == 0 && k == 0 && l == 0 && m == 0) {
                                table.axes.momenta[n] = momentum * kMeVToGeV;
                            }
                            if (i == 0 && j == 0 && k == 0 && l == 0 && n == 0) {
                                table.axes.masses[m] = mass;
                            }

                            const std::size_t index = table.flat_index(i, j, k, l, n, m);
                            table.channel_values[0][index] = phi_value;
                        }
                    }
                }
            }
        }
    }

    table.finalize_metadata();
    std::cout << "Loaded phi meson spectral function from " << filepath << "\n";
    return table;
}

} // namespace FluidDileptons::Rates