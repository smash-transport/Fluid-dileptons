#pragma once

#include "vector.h"
#include "setup.h"

namespace FluidDileptons {

namespace AcceptanceCutter {
  static constexpr double large_cut = 100;

  void set_x_range(double min, double max);
  void set_y_range(double min, double max);
  void set_eta_range(double min, double max);
  void set_pT_range(double min, double max);
  void set_yrap_range(double min, double max);
  bool in_spatial_range(const FourVector& pos);
  bool in_momentum_range(const FourVector& mom);
}

}
