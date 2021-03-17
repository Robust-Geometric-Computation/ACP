#ifndef PVD_H
#define PVD_H

#include "acp/core/acp.h"
#include "acp/linmath/pv.h"

namespace acp {

class PV2D {
 public:
  int size() { return 2; }
  double& operator[](int i) { return i == 0 ? x : y; }

  PV2D() {}
  PV2D(const double ix, const double iy) : x(ix), y(iy) {}
  PV2D(const PV2D& p) : x(p.x), y(p.y) {}
  template <class N>
  PV2D(const PV2<N>& p) : x(p.x.mid()), y(p.y.mid()) {}
  PV2D& operator=(const PV2D& p) {
    x = p.x;
    y = p.y;
    return *this;
  }
  double getX() const { return x; }
  double getY() const { return y; }
  double dot(const PV2D& b) const { return x * b.x + y * b.y; }
  double cross(const PV2D& b) const { return x * b.y - y * b.x; }
  PV2D operator+(const PV2D& b) const { return PV2D(x + b.x, y + b.y); }
  PV2D operator-(const PV2D& b) const { return PV2D(x - b.x, y - b.y); }
  PV2D operator-() const { return PV2D(-x, -y); }
  PV2D operator*(double k) const { return PV2D(x * k, y * k); }
  PV2D operator/(double k) const { return PV2D(x / k, y / k); }
  PV2D unit() const { return *this / sqrt((*this).dot(*this)); }
  double length() const { return sqrt((*this).dot(*this)); }

  // 2 times the signed area of the triangle this,b,c.
  double area(const PV2D& b, const PV2D& c) const {
    return (*this - c).cross(b - c);
  }

  double x, y;
};

inline PV2D operator*(double k, PV2D v) { return v * k; }

class PV3D {
 public:
  int size() { return 3; }

  PV3D() {}
  PV3D(const double ix, const double iy, const double iz)
      : x(ix), y(iy), z(iz) {}

  PV3D(const PV3D& p) : x(p.x), y(p.y), z(p.z) {}
  template <class N>
  PV3D(const PV3<N>& p) : x(p.x.mid()), y(p.y.mid()), z(p.z.mid()) {}
  PV3D& operator=(const PV3D& p) {
    x = p.x;
    y = p.y;
    z = p.z;
    return *this;
  }

  double getX() const { return x; }
  double getY() const { return y; }
  double getZ() const { return z; }

  double& operator[](int i) { return i == 0 ? x : i == 1 ? y : z; }

  double dot(const PV3D& b) const { return x * b.x + y * b.y + z * b.z; }

  PV3D operator+(const PV3D& b) const {
    return PV3D(x + b.x, y + b.y, z + b.z);
  }

  PV3D operator-() const { return PV3D(-x, -y, -z); }

  PV3D operator-(const PV3D& b) const {
    return PV3D(x - b.x, y - b.y, z - b.z);
  }

  PV3D operator*(double k) const { return PV3D(x * k, y * k, z * k); }

  PV3D operator/(double k) const { return PV3D(x / k, y / k, z / k); }

  PV3D cross(const PV3D& b) const {
    return PV3D(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
  }

  double tripleProduct(const PV3D& b, const PV3D& c) const {
    return dot(b.cross(c));
  }

  PV3D unit() const { return *this / sqrt((*this).dot(*this)); }

  double x, y, z;
};

inline PV3D operator*(double k, PV3D v) { return v * k; }

}  // namespace acp

#endif
