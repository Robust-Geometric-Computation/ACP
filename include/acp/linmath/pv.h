/*
  ACP (Adaptive Controlled Precision) Library
  for robust computational geometry

  Copyright (c) 2013-07-15
  Victor Joseph Milenkovic
  University of Miami
  vjm@miami.edu
  Elisha Sacks
  Purdue University
  eps@cs.purdue.edu

  This file contains code described in

Robust Complete Path Planning in the Plane
Victor Milenkovic, Elisha Sacks, and Steven Trac
Proceedings of the Workshop on the Algorithmic Foundations of Robotics (WAFR)
pages 37-52, 2012

   This code is under development.
   It is free for use for academic and research purposes.
   Permission from the authors is required for commercial use.
*/

#ifndef PV_H
#define PV_H

#include "acp/core/acp.h"

namespace acp {

template <class N>
class Scalar {
 public:
  int size() const { return 1; }
  const N& operator[](int i) const { return x; }
  N& operator[](int i) { return x; }
  operator const N&() const { return x; }
  operator N&() { return x; }

  Scalar() {}
  Scalar(const N& x) : x(x) {}
  template <class M>
  Scalar(const Scalar<M>& p) : x(p.x) {}
  Scalar& operator=(const Scalar& p) {
    x = p.x;
    return *this;
  }
  bool uninitialized() { return x.uninitialized(); }
  N getX() const { return x; }
  Scalar operator+(const Scalar& b) const { return Scalar(x + b.x); }
  Scalar operator-(const Scalar& b) const { return Scalar(x - b.x); }
  Scalar operator-() const { return Scalar(-x); }
  Scalar operator*(double k) const { return Scalar(x * k); }
  Scalar operator*(const N& k) const { return Scalar(x * k); }
  Scalar operator/(double k) const { return Scalar(x / k); }
  Scalar operator/(const N& k) const { return Scalar(x / k); }

  N x;
};

template <class N>
class PV2 {
 public:
  int size() const { return 2; }
  const N& operator[](int i) const { return i == 0 ? x : y; }
  N& operator[](int i) { return i == 0 ? x : y; }

  PV2() {}
  PV2(const N& x, const N& y) : x(x), y(y) {}
  static PV2<N> constant(double x, double y) {
    return PV2<N>(N::constant(x), N::constant(y));
  }
  static PV2<N> input(double x, double y) {
    return PV2<N>(N::input(x), N::input(y));
  }
  template <class M>
  PV2(const PV2<M>& p) : x(p.x), y(p.y) {}
  PV2& operator=(const PV2& p) {
    x = p.x;
    y = p.y;
    return *this;
  }
  bool uninitialized() { return x.uninitialized() && y.uninitialized(); }
  N getX() const { return x; }
  N getY() const { return y; }
  N dot(const PV2& b) const { return x * b.x + y * b.y; }
  N cross(const PV2& b) const { return x * b.y - y * b.x; }
  PV2 operator+(const PV2& b) const { return PV2(x + b.x, y + b.y); }
  PV2 operator-(const PV2& b) const { return PV2(x - b.x, y - b.y); }
  PV2 operator-() const { return PV2(-x, -y); }
  PV2 operator*(double k) const { return PV2(x * k, y * k); }
  PV2 operator*(const N& k) const { return PV2(x * k, y * k); }
  PV2 operator/(double k) const { return PV2(x / k, y / k); }
  PV2 operator/(const N& k) const { return PV2(x / k, y / k); }
  PV2 unit() const { return *this / (*this).dot(*this).sqrt(); }

  N length() const { return (*this).dot(*this).sqrt(); }

  // 2 times the signed area of the triangle this,b,c.
  N area(const PV2& b, const PV2& c) const { return (*this - c).cross(b - c); }

  // rotate by b
  PV2 rotate(const PV2& b) const {
    return PV2(b.x * x - b.y * y, b.y * x + b.x * y);
  }

  N x, y;
};

template <class N>
inline PV2<N> operator*(double k, PV2<N> v) {
  return v * k;
}

template <class N>
inline PV2<N> operator*(const N& k, PV2<N> v) {
  return v * k;
}

template <class N>
class PV3 {
 public:
  int size() const { return 3; }

  PV3() {}
  PV3(const N& ix, const N& iy, const N& iz) : x(ix), y(iy), z(iz) {}
  template <class M>
  PV3(const PV3<M>& p) : x(p.x), y(p.y), z(p.z) {}
  PV3& operator=(const PV3& p) {
    x = p.x;
    y = p.y;
    z = p.z;
    return *this;
  }

  static PV3<N> input(double x, double y, double z) {
    return PV3<N>(N::input(x), N::input(y), N::input(z));
  }

  static PV3<N> constant(double x, double y, double z) {
    return PV3<N>(N::constant(x), N::constant(y), N::constant(z));
  }

  bool uninitialized() {
    return x.uninitialized() && y.uninitialized() && z.uninitialized();
  }

  N getX() const { return x; }
  N getY() const { return y; }
  N getZ() const { return z; }

  const N& operator[](int i) const { return i == 0 ? x : i == 1 ? y : z; }
  N& operator[](int i) { return i == 0 ? x : i == 1 ? y : z; }

  N dot(const PV3& b) const { return x * b.x + y * b.y + z * b.z; }

  PV3 operator+(const PV3& b) const { return PV3(x + b.x, y + b.y, z + b.z); }

  PV3 operator-() const { return PV3(-x, -y, -z); }

  PV3 operator-(const PV3& b) const { return PV3(x - b.x, y - b.y, z - b.z); }

  PV3 operator*(double k) const { return PV3(x * k, y * k, z * k); }

  PV3 operator*(const N& k) const { return PV3(x * k, y * k, z * k); }

  PV3 operator/(double k) const { return PV3(x / k, y / k, z / k); }

  PV3 operator/(const N& k) const { return PV3(x / k, y / k, z / k); }

  PV3 cross(const PV3& b) const {
    return PV3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
  }

  N tripleProduct(const PV3& b, const PV3& c) const { return dot(b.cross(c)); }

  PV3 unit() const { return *this / (*this).dot(*this).sqrt(); }

  N length2() const { return (*this).dot(*this); }

  N length() const { return (*this).dot(*this).sqrt(); }

  // rotate by b around z axix
  PV3 rotateZ(const PV2<N>& b) const {
    return PV3(b.x * x - b.y * y, b.y * x + b.x * y, z);
  }

  N x, y, z;
};

template <class N>
inline PV3<N> operator*(double k, PV3<N> v) {
  return v * k;
}

template <class N>
inline PV3<N> operator*(const N& k, PV3<N> v) {
  return v * k;
}

}  // namespace acp

#endif
