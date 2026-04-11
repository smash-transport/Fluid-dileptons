#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace FluidDileptons::Rates {

struct TableAxes {
    std::vector<double> temperatures;
    std::vector<double> baryon_densities;
    std::vector<double> pion_chemical_potentials;
    std::vector<double> kaon_chemical_potentials;
    std::vector<double> momenta;
    std::vector<double> masses;
};

struct TableLayout {
    std::array<std::size_t, 6> shape;
    bool has_secondary_channel = false;
};

inline constexpr TableLayout rho_omega_layout{{32, 13, 4, 1, 21, 150}, true};
inline constexpr TableLayout phi_layout{{32, 13, 4, 4, 21, 101}, false};

struct AxisBracket {
    std::size_t lower = 0;
    std::size_t upper = 0;
    double weight_upper = 0.0;

    bool fixed() const { return lower == upper; }
};

struct AxisInfo {
    bool regular = false;
    double min = 0.0;
    double step = 0.0;
};

struct VectorMesonCellContext {
    std::array<AxisBracket, 4> fixed_axes;
};

struct VectorMesonGridCache {
    std::vector<AxisBracket> momenta;
    std::vector<AxisBracket> masses;
};

struct VectorMesonTable {
    TableLayout layout;

    std::vector<std::string> channel_names;
    TableAxes axes;
    std::vector<std::vector<double>> channel_values;
    std::array<std::size_t, 6> strides{};
    std::array<AxisInfo, 6> axis_info{};

    std::size_t total_size() const {
        return layout.shape[0] * layout.shape[1] * layout.shape[2] *
               layout.shape[3] * layout.shape[4] * layout.shape[5];
    }

    std::size_t temperature_bins() const { return layout.shape[0]; }
    std::size_t density_bins() const { return layout.shape[1]; }
    std::size_t pion_chemical_potential_bins() const { return layout.shape[2]; }
    std::size_t kaon_chemical_potential_bins() const { return layout.shape[3]; }
    std::size_t momentum_bins() const { return layout.shape[4]; }
    std::size_t mass_bins() const { return layout.shape[5]; }

    bool has_secondary_channel() const { return layout.has_secondary_channel; }

    std::size_t channel_count() const { return channel_values.size(); }

    const std::array<std::size_t, 6>& shape() const {
        return layout.shape;
    }

    void finalize_metadata();
    std::size_t channel_index(const std::string& channel_name) const;
    AxisBracket bracket(std::size_t axis_index, double value) const;

    VectorMesonCellContext make_cell_context(double temperature,
                                             double density,
                                             double pion_chemical_potential = 0.0,
                                             double kaon_chemical_potential = 0.0) const;

    VectorMesonGridCache make_grid_cache(const std::vector<double>& momenta,
                                         const std::vector<double>& masses) const;

    std::size_t flat_index(std::size_t temperature_index,
                           std::size_t density_index,
                           std::size_t pion_index,
                           std::size_t kaon_index,
                           std::size_t momentum_index,
                           std::size_t mass_index) const;

    double at(std::size_t channel,
              std::size_t temperature_index,
              std::size_t density_index,
              std::size_t pion_index,
              std::size_t kaon_index,
              std::size_t momentum_index,
              std::size_t mass_index) const;

    double at(const std::string& channel_name,
              std::size_t temperature_index,
              std::size_t density_index,
              std::size_t pion_index,
              std::size_t kaon_index,
              std::size_t momentum_index,
              std::size_t mass_index) const;

    double interpolate(std::size_t channel,
                       double temperature,
                       double density,
                       double pion_chemical_potential,
                       double kaon_chemical_potential,
                       double momentum,
                       double mass) const;

    double interpolate_fast(std::size_t channel,
                            const VectorMesonCellContext& cell_context,
                            const AxisBracket& momentum_bracket,
                            const AxisBracket& mass_bracket) const;

    double interpolate(const std::string& channel_name,
                       double temperature,
                       double density,
                       double pion_chemical_potential,
                       double kaon_chemical_potential,
                       double momentum,
                       double mass) const;
};

VectorMesonTable load_rho_omega_table(const std::string& filepath);
VectorMesonTable load_phi_table(const std::string& filepath);

} // namespace FluidDileptons::Rates