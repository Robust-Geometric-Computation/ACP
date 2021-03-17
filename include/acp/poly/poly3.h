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

#ifndef POLY3_H
#define POLY3_H

#include <unordered_map>
#include <vector>
#include "acp/poly/poly2.h"
#include "acp/poly/poly3d.h"

using namespace std;
using namespace acp;

namespace acp {

namespace Poly3Inner {
class Index3 {
 public:
  int i[3];
  Index3(int i0, int i1, int i2) {
    i[0] = i0;
    i[1] = i1;
    i[2] = i2;
  }
};

class Hasher {
 public:
  size_t operator()(const Index3& a) const {
    return a.i[0] * 779230947 + a.i[1] * 247091631 + a.i[2] * 1194289623;
  }
};

class Equals {
 public:
  bool operator()(const Index3& a, const Index3& b) const {
    return a.i[0] == b.i[0] && a.i[1] == b.i[1] && a.i[2] == b.i[2];
  }
};
}  // namespace Poly3Inner

template <class N>
class Poly3 {
  typedef Poly3Inner::Index3 Index3;
  typedef Poly3Inner::Hasher Hasher;
  typedef Poly3Inner::Equals Equals;

 public:
  Poly3() : xd(-1), yd(-1), zd(-1), d(-1) {}

  void add(int i, int j, int k, const N& c) {
    auto it = ia.find(Index3(i, j, k));
    if (it == ia.end()) {
      ia[Index3(i, j, k)] = a.size();
      a.push_back(c);
    } else
      a[it->second] = a[it->second] + c;

    xd = i > xd ? i : xd;
    yd = j > yd ? j : yd;
    zd = k > zd ? k : zd;

    d = (i + j + k) > d ? (i + j + k) : d;
  }

  template <class M>
  Poly3(const Poly3<M>& p)
      : xd(p.xd), yd(p.yd), zd(p.zd), d(p.d), ia(p.ia), a(p.a.size()) {
    for (int i = 0; i < a.size(); i++) a[i] = p.a[i];
  }

  template <class M>
  Poly3(const Poly3D<M>& p)
      : xd(p.degX()),
        yd(p.degY()),
        zd(p.degZ()),
        d(p.degX() + p.degY() + p.degZ()) {
    for (int ix = 0; ix <= xd; ++ix)
      for (int iy = 0; iy <= yd; ++iy)
        for (int iz = 0; iz <= zd; ++iz) add(ix, iy, iz, p.get(ix, iy, iz));
  }

  N value(const PV3<N>& p) const { return value(p.x, p.y, p.z); }

  N value(const N& x, const N& y, const N& z) const {
    const N* xyz[3] = {&x, &y, &z};

    vector<N> powers[3];
    for (int i = 0; i < 3; i++) powers[i].push_back(N::constant(1));

    N sum = N::constant(0);

    for (auto pair = ia.begin(); pair != ia.end(); pair++) {
      Index3 ind = pair->first;
      N term = a[pair->second];

      for (int i = 0; i < 3; i++) {
        while (powers[i].size() <= ind.i[i])
          powers[i].push_back(*xyz[i] * powers[i].back());
        term = term * powers[i][ind.i[i]];
      }

      sum = sum + term;
    }

    return sum;
  }

  Poly3 der(int xyz) const {
    Poly3 that;

    Index3 ind(-1, -1, -1);
    int i = -1;

    for (auto [ind, i] : ia) {
      if (ind.i[xyz] > 0)
        that.add(ind.i[0] - !xyz, ind.i[1] - (xyz % 2), ind.i[2] - (xyz / 2),
                 ind.i[xyz] * a[i]);
    }
    return that;
  }

  PV3<N> gradient(const N& x, const N& y, const N& z) const {
    const N* xyz[3] = {&x, &y, &z};

    vector<N> powers[3];
    for (int i = 0; i < 3; i++) powers[i].push_back(N::constant(1));

    PV3<N> sum = PV3<N>(N(0), N(0), N(0));

    for (auto pair = ia.begin(); pair != ia.end(); pair++) {
      Index3 ind = pair->first;

      for (int i = 0; i < 3; i++)
        while (powers[i].size() <= ind.i[i])
          powers[i].push_back(*xyz[i] * powers[i].back());

      for (int i = 0; i < 3; i++) {
        if (ind.i[i] == 0) continue;
        N term = a[pair->second] * ind.i[i];
        for (int j = 0; j < 3; j++)
          term = term * powers[j /*i*/][ind.i[j] - (i == j)];
        sum[i] = sum[i] + term;
      }
    }

    return sum;
  }

  Poly3 operator+(const Poly3& that) const {
    Poly3 q;

    for (auto it = ia.begin(); it != ia.end(); it++) {
      Index3 i = it->first;
      int ind = it->second;
      q.add(i.i[0], i.i[1], i.i[2], a[ind]);
    }

    for (auto it = that.ia.begin(); it != that.ia.end(); it++) {
      Index3 i = it->first;
      int ind = it->second;
      q.add(i.i[0], i.i[1], i.i[2], that.a[ind]);
    }

    return q;
  }

  Poly3 operator-(const Poly3& that) const {
    Poly3 q;

    for (auto it = ia.begin(); it != ia.end(); it++) {
      Index3 i = it->first;
      int ind = it->second;
      q.add(i.i[0], i.i[1], i.i[2], a[ind]);
    }

    for (auto it = that.ia.begin(); it != that.ia.end(); it++) {
      Index3 i = it->first;
      int ind = it->second;
      q.add(i.i[0], i.i[1], i.i[2], -that.a[ind]);
    }

    return q;
  }

  Poly3 operator*(const Poly3& that) const {
    Poly3 q;

    for (auto it = ia.begin(); it != ia.end(); it++) {
      Index3 ai = it->first;
      int aai = it->second;

      for (auto it = that.ia.begin(); it != that.ia.end(); it++) {
        Index3 bi = it->first;
        int bai = it->second;

        int i = ai.i[0] + bi.i[0];
        int j = ai.i[1] + bi.i[1];
        int k = ai.i[2] + bi.i[2];

        q.add(i, j, k, a[aai] * that.a[bai]);
      }
    }

    return q;
  }

  // xyz = 0, 1, 2:  substituting for yz, xz, xy.
  Poly<N> substitute2(int ixyz, const N& x, const N& y, const N& z) const {
    const N* xyz[3] = {&x, &y, &z};

    vector<N> powers[3];
    for (int i = 0; i < 3; i++) powers[i].push_back(N::constant(1));

    Poly<N> ppoly;

    for (auto pair = ia.begin(); pair != ia.end(); pair++) {
      Index3 ind = pair->first;
      N term = a[pair->second];

      for (int i = 0; i < 3; i++)
        if (i != ixyz) {
          while (powers[i].size() <= ind.i[i])
            powers[i].push_back(*xyz[i] * powers[i].back());
          term = term * powers[i][ind.i[i]];
        }

      ppoly.add(term, ind.i[ixyz]);
    }

    return ppoly;
  }

  Poly2<N> substitute1(int ixyz, const N& xyz) const {
    int jxyz = (ixyz + 1) % 3;
    int kxyz = (jxyz + 1) % 3;

    vector<N> powers;
    powers.push_back(N(1));

    Poly2<N> poly2;

    for (auto pair = ia.begin(); pair != ia.end(); pair++) {
      Index3 ind = pair->first;
      N term = a[pair->second];

      while (powers.size() <= ind.i[ixyz])
        powers.push_back(xyz * powers.back());

      poly2.add(ind.i[jxyz], ind.i[/*jxyz*/ kxyz], term * powers[ind.i[ixyz]]);
    }

    return poly2;
  }

  // t is a flattened 4x4 matrix in column major form
  Poly3 transform(double t[16]) const {
    Poly3 x;
    Poly3 y;
    Poly3 z;
    Poly3 w;

    Poly3 unit;

    // First column
    x.add(1, 0, 0, N::constant(t[0]));
    y.add(1, 0, 0, N::constant(t[1]));
    z.add(1, 0, 0, N::constant(t[2]));
    w.add(1, 0, 0, N::constant(t[3]));

    // Second column
    x.add(0, 1, 0, N::constant(t[4]));
    y.add(0, 1, 0, N::constant(t[5]));
    z.add(0, 1, 0, N::constant(t[6]));
    w.add(0, 1, 0, N::constant(t[7]));

    // Third column
    x.add(0, 0, 1, N::constant(t[8]));
    y.add(0, 0, 1, N::constant(t[9]));
    z.add(0, 0, 1, N::constant(t[10]));
    w.add(0, 0, 1, N::constant(t[11]));

    // Fourth column
    x.add(0, 0, 0, N::constant(t[12]));
    y.add(0, 0, 0, N::constant(t[13]));
    z.add(0, 0, 0, N::constant(t[14]));
    w.add(0, 0, 0, N::constant(t[15]));

    unit.add(0, 0, 0, N::constant(1));

    std::vector<Poly3> xpow(xd + 1);
    std::vector<Poly3> ypow(yd + 1);
    std::vector<Poly3> zpow(zd + 1);
    std::vector<Poly3> wpow(d + 1);

    xpow[0] = ypow[0] = zpow[0] = wpow[0] = unit;

    for (int i = 1; i <= xd; i++) xpow[i] = x * xpow[i - 1];
    for (int i = 1; i <= yd; i++) ypow[i] = y * ypow[i - 1];
    for (int i = 1; i <= zd; i++) zpow[i] = z * zpow[i - 1];
    for (int i = 1; i <= d; i++) wpow[i] = w * wpow[i - 1];

    Poly3 sum;

    for (auto it = ia.begin(); it != ia.end(); it++) {
      Index3 ai = it->first;
      int aai = it->second;

      int i = ai.i[0];
      int j = ai.i[1];
      int k = ai.i[2];

      Poly3 coef;
      coef.add(0, 0, 0, a[aai]);

      sum = sum + (coef * xpow[i] * ypow[j] * zpow[k] * wpow[d - (i + j + k)]);
    }

    return sum;
  }

  unordered_map<Index3, int, Hasher, Equals> ia;
  vector<N> a;
  int size() const { return a.size(); }
  const N& operator[](int i) const { return a[i]; }
  N& operator[](int i) { return a[i]; }
  int xd, yd, zd, d;
};

}  // namespace acp
#endif
