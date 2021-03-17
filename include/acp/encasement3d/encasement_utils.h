#pragma once

#include "acp/linmath/pv.h"
#include "acp/poly/poly3.h"

// Utilities for the encasement algorithm
namespace acp {
namespace encasement_utils {

// Convenince Double class
class Double {
 public:
  Double() : d(0.0) {}
  Double(const double d) : d(d) {}

  Double(const Double& e) : d(e.mid()) {}
  Double(const Parameter& e) : d(e.mid()) {}
  Double(const PParameter& e) : d(e.mid()) {}
  Double(const MParameter& e) : d(e.mid()) {}

  int size() const { return 1; }

  // Doubles can never be uninitialized
  bool uninitialized() const { return false; }

  Double& operator=(const Double& x) {
    d = x.d;
    return *this;
  }

  static Double constant(const double& d) { return Double(d); }

  bool operator==(const Double& x) const { return d == x.d; }

  Double operator-() const { return -d; }

  Double operator+(double x) const { return d + x; }
  Double operator-(double x) const { return d - x; }
  Double operator*(double x) const { return d * x; }
  Double operator/(double x) const { return d / x; }

  Double operator+(const Double& x) const { return d + x.d; }
  Double operator-(const Double& x) const { return d - x.d; }
  Double operator*(const Double& x) const { return d * x.d; }
  Double operator/(const Double& x) const { return d / x.d; }

  operator double() const { return d; }

  Double sqrt() const { return Double(std::sqrt(d)); }

  int sign(bool fail = true) const { return d == 0 ? 0 : (d < 0 ? -1 : 1); }

  double lb() const { return d; }
  double mid() const { return d; }
  double ub() const { return d; }

  Double lbP() const { return Double(d); }
  Double midP() const { return Double(d); }
  Double ubP() const { return Double(d); }

  Double abs() const { return std::abs(d); }
  bool operator<(const Double& b) const { return d < b.d; }
  bool operator>(const Double& b) const { return d > b.d; }

 private:
  double d;
};

std::map<std::pair<int, int>, std::vector<PV3<Parameter>>> FGCurves(
    const Poly3<Parameter>& f, const Poly3<Parameter>& g,
    const PV3<Parameter>& box,
    const std::vector<PV3<Parameter>>& boundary_verts, Double tolerance = 1e-6);

}  // namespace encasement_utils
}  // namespace acp
