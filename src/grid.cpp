#include "spectrum.h"

namespace FluidDileptons {

/*
 The reset is a bit awkward but needed for the setup. Also allowing the user to have
 different grids is not bad but it should be done consciously, as it affects the output.
*/
namespace MQGrid {
    namespace {
        bool masses_set = false, qs_set = false;
    }
    std::vector<double> masses, qs;

    static void set_vector(std::vector<double>& vec, double min, double max, double step) {
        vec.clear();
        if (min + step >= max) {
            vec.push_back(min);
        } else {
            for (size_t i = 0; min + i*step <= max; ++i)
                vec.push_back(min + i*step);
        }
    }

    void set_masses(double min, double max, double step, bool reset) {
        if (masses_set && !reset) {
            std::cerr << "Warning: masses already set, ignoring.\n";
            return;
        }
        set_vector(masses, min, max, step);
        masses_set = true;
    }

    void set_masses(std::vector<double> vec, bool reset) {
        if (masses_set && !reset) {
            std::cerr << "Warning: masses already set, ignoring.\n";
            return;
        }
        masses = vec;
        masses_set = true;
    }

    void set_qs(double min, double max, double step, bool reset) {
        if (qs_set && !reset) {
            std::cerr << "Warning: qs already set, ignoring.\n";
            return;
        }
        set_vector(qs, min, max, step);
        qs_set = true;
    }

    void set_qs(std::vector<double> vec, bool reset) {
        if (qs_set && !reset) {
            std::cerr << "Warning: qs already set, ignoring.\n";
            return;
        }
        qs = vec;
        qs_set = true;
    }

    static size_t find_axis_index(const std::vector<double>& axis, double value) {
        const auto it = std::lower_bound(axis.begin(), axis.end(), value);

        size_t index = 0;
        if (it == axis.begin()) {
            index = 0;
        } else if (it == axis.end()) {
            index = axis.size() - 1;
        } else {
            const size_t right = static_cast<size_t>(it - axis.begin());
            const size_t left = right - 1;
            index = (std::abs(axis[left] - value) < std::abs(axis[right] - value)) ? left : right;
        }
        return index;
    }

    std::pair<size_t, size_t> find_indices(double mass, double q) {
        return {find_axis_index(masses, mass), find_axis_index(qs, q)};
    }
}

} // FluidDileptons
