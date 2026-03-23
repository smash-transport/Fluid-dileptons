#include <stdexcept>
#include <unordered_map>

#include "hydrocell.h"
#include "spectrum.h"

namespace FluidDileptons {

Spectrum::Spectrum(Source source) :
    source_(source) {
    weights_.assign(MQGrid::masses.size() * MQGrid::mom_abs.size(), 0.0);
}
void Spectrum::reset() {
    std::fill(weights_.begin(), weights_.end(), 0.0);
}
void Spectrum::fill(const Dilepton& dilepton) {
    if (dilepton.source() != source_)
        return;
    add(dilepton.m(), dilepton.q(), dilepton.weight());
}
void Spectrum::fill(const DileptonList& dileptons) {
    for (const Dilepton& d : dileptons)
        fill(d);
}
void Spectrum::fill(const HydroCell& cell) {
    fill(cell.dileptons());
}

void Spectrum::add(double mass, double q, double weight) {
    const size_t idx = index(mass, q);
    weights_[idx] += weight;
}
void Spectrum::add(std::pair<size_t,size_t> indices, double weight) {
    const size_t idx = index(indices);
    weights_[idx] += weight;
}

double Spectrum::weight_at(double mass, double q) const {
    return weights_[index(mass, q)];
}
size_t Spectrum::index(double mass, double q) const {
    return index(MQGrid::find_indices(mass, q));
}

// Using std::pair avoids having two overloads (of weight_at) with
// signatures that could be implicitly converted
double Spectrum::weight_at(std::pair<size_t,size_t> indices) const {
    const auto [mass_index, q_index] = indices;
    return weights_[mass_index * MQGrid::mom_abs.size() + q_index];
}
size_t Spectrum::index(std::pair<size_t,size_t> indices) const {
    const auto [mass_index, q_index] = indices;
    return mass_index * MQGrid::mom_abs.size() + q_index;
}

namespace Spectra {
    namespace {
        bool initialized = false;
        std::unordered_map<Source, Spectrum> all_spectra;

        // Hidden to forbid direct access, but it can be exposed if need be
        Spectrum& get_mutable(Source source) {
            auto it = all_spectra.find(source);
            if (it == all_spectra.end())
                throw std::invalid_argument("Unknown dilepton source.");
            return it->second;
        }
    }

    void initialize() {
        if (initialized)
            throw std::runtime_error("Spectra already initialized.");
        for (Source s : all_sources)
            all_spectra.emplace(s, Spectrum(s));
        initialized = true;
    }
    void reset() {
        for (auto& pair : all_spectra)
            pair.second.reset();
        initialized = false;
    }
    const Spectrum& get(Source source) {
        return get_mutable(source);
    }
    void add(Source source, double m, double q, double weight) {
        if (weight == 0.0)
            return;
        get_mutable(source).add(m, q, weight);
    }
    void add(Source source, std::pair<size_t,size_t> indices, double weight) {
        if (weight == 0.0)
            return;
        get_mutable(source).add(indices, weight);
    }
    void fill(const Dilepton& dilepton) {
        get_mutable(dilepton.source()).fill(dilepton);
    }
    void fill(const DileptonList& dileptons) {
        for (const Dilepton& d : dileptons)
            fill(d);
    }
    void fill(const HydroCell& cell) {
        fill(cell.dileptons());
    }
} // namespace Spectra

} // namespace FluidDileptons
