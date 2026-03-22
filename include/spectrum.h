#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "dilepton.h"

namespace FluidDileptons {

namespace MQGrid {
    extern std::vector<double> masses, qs;
    void set_masses(double min, double max, double step, bool reset = false);
    void set_masses(std::vector<double> vec, bool reset = false);
    void set_qs(double min, double max, double step, bool reset = false);
    void set_qs(std::vector<double> vec, bool reset = false);
    std::pair<size_t, size_t> find_indices(double mass, double q);
}

// Histogram of dilepton weights on the fixed MQGrid (mass, q) for a single source.
class Spectrum {
  public:
    explicit Spectrum(Source source);

    Source source() const { return source_; }
    const std::vector<double>& masses() const { return MQGrid::masses; }
    const std::vector<double>& qs() const { return MQGrid::qs; }

    void reset() {
        std::fill(weights_.begin(), weights_.end(), 0.0);
    }

    void add(double mass, double q, double weight);

    void fill(const Dilepton& dilepton);
    void fill(const DileptonList& dileptons);
    void fill(const HydroCell& cell);

    double weight_at(double mass, double q) const;
    double weight_at(std::pair<size_t,size_t> indices) const;
    const std::vector<double>& data() const { return weights_; }

  private:

    size_t index(double mass, double q) const;
    size_t index(std::pair<size_t,size_t> indices) const;

    Source source_;
    std::vector<double> weights_;
};

// Global registry for all sources.
namespace Spectra {
    void initialize();
    void reset();

    const Spectrum& get(Source source);

    void fill(const Dilepton& dilepton);
    void fill(const DileptonList& dileptons);
    void fill(const HydroCell& cell);
}

} // namespace FluidDileptons
