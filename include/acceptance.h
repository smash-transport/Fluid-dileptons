#pragma once

#include "vector.h"
#include "setup.h"

namespace FluidDileptons {

namespace AcceptanceCutter {
  static constexpr double large_cut = 100;
  namespace {
    double eta_s_min_ = -large_cut, eta_s_max_ = large_cut;
    double x_min_ = -large_cut, x_max_ = large_cut;
    double y_min_ = -large_cut, y_max_ = large_cut;
    double pT_min_ = 0, pT_max_ = large_cut;
    double yrap_min_ = -large_cut, yrap_max_ = large_cut;
  }
  inline void set_x_range(double min, double max) {
    x_min_ = min;
    x_max_ = max;
  }
  inline void set_y_range(double min, double max) {
    y_min_ = min;
    y_max_ = max;
  }
  inline void set_eta_range(double min, double max) {
    eta_s_min_ = min;
    eta_s_max_ = max;
  }
  inline void set_pT_range(double min, double max) {
    pT_min_ = min;
    pT_max_ = max;
  }
  inline void set_yrap_range(double min, double max) {
    yrap_min_ = min;
    yrap_max_ = max;
  }
  inline bool in_spatial_range(const FourVector& pos) {
    const bool in_eta_s_range = (eta_s_min_<= pos.eta()) && (pos.eta() < eta_s_max_);
    const bool in_x_range = (x_min_<= pos.x1()) && (pos.x1() < x_max_);
    const bool in_y_range = (y_min_<= pos.x2()) && (pos.x2() < y_max_);
    return in_eta_s_range && in_x_range && in_y_range;
  }
  inline bool in_momentum_range(const FourVector& mom) {
    const bool in_yrap_range = (yrap_min_ <= mom.eta()) && (mom.eta() <= yrap_max_);
    const bool in_pT_range = (pT_min_ <= mom.xT()) && (mom.xT() <= pT_max_);
    return in_yrap_range && in_pT_range;
  }
}

}
