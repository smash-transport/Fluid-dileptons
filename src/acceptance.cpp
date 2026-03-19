#include "acceptance.h"

namespace FluidDileptons {

namespace AcceptanceCutter {

  namespace {
    double eta_s_min, eta_s_max, x_min, x_max, y_min, y_max;
    double pT_min, pT_max, yrap_min, yrap_max;
  }

  void set_x_range(double min, double max) {
    x_min = min;
    x_max = max;
  }
  void set_y_range(double min, double max) {
    y_min = min;
    y_max = max;
  }
  void set_eta_range(double min, double max) {
    eta_s_min = min;
    eta_s_max = max;
  }
  void set_pT_range(double min, double max) {
    pT_min = min;
    pT_max = max;
  }
  void set_yrap_range(double min, double max) {
    yrap_min = min;
    yrap_max = max;
  }

  bool in_spatial_range(const FourVector& pos) {
    const bool in_eta_s_range = (eta_s_min <= pos.eta()) && (pos.eta() < eta_s_max);
    const bool in_x_range = (x_min <= pos.x1()) && (pos.x1() < x_max);
    const bool in_y_range = (y_min <= pos.x2()) && (pos.x2() < y_max);
    return in_eta_s_range && in_x_range && in_y_range;
  }

  bool in_momentum_range(const FourVector& mom) {
    const bool in_yrap_range = (yrap_min <= mom.eta()) && (mom.eta() <= yrap_max);
    const bool in_pT_range = (pT_min <= mom.xT()) && (mom.xT() <= pT_max);
    return in_yrap_range && in_pT_range;
  }

}

}