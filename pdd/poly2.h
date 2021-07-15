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

#ifndef POLY2_H
#define POLY2_H

#include <unordered_map>
#include <vector>
//#include "acp/poly/poly.h"
#include "poly.h"
using namespace std;
using namespace acp;

namespace acp {

template <class N>
class Poly2 {
 public:
  Poly2() : d(-2), degx(-2), degy(-2) {}
  Poly2(int nt) : a(nt), m(2 * nt), d(-2), degx(-2), degy(-2) {}
  Poly2(int nt, N* a, int* m) : a(nt), m(2 * nt), d(-2), degx(-2), degy(-2) {
    int k = nt * 2;
    for (int i = 0; i < nt; ++i) this->a[i] = a[i];
    for (int i = 0; i < k; ++i) this->m[i] = m[i];
    setDegree();
    degx = degX();
    degy = degY();
  }
  template <class M>
  Poly2(const Poly2<M>& p)
      : a(p.a.size()), m(p.m), d(p.d), degx(p.degx), degy(p.degy) {
    for (int i = 0; i < a.size(); ++i) a[i] = N(p.a[i]);
  }
  void add(int i, int j, const N& c) {
    m.push_back(i);
    m.push_back(j);
    a.push_back(c);

    degx = i > degx ? i : degx;
    degy = j > degy ? j : degy;
    d = (i + j) > d ? (i + j) : d;
  }
  void setDegree() {
    d = 0;
    for (int i = 0; i < a.size(); ++i) {
      int di = 0;
      for (int j = 0; j < 2; ++j) di += m[2 * i + j];
      if (d < di) d = di;
    }
    degx = degX();
    degy = degY();
  }

  int degree() const { return d; }
  bool increased() const { return false; /* a[0].increased(); */ }
  bool decreased() const { return true; /* a[0].decreased(); */ }
  N value(N* x) const {
    N xp[degx + 1];
    N yp[degy + 1];

    N p = x[0];
    for (int i = 1; i <= degx; i++) {
      xp[i] = p;
      p = p * x[0];
    }

    p = x[1];
    for (int i = 1; i <= degy; i++) {
      yp[i] = p;
      p = p * x[1];
    }

    N y;
    auto mp = m.begin();
    for (int i = 0; i < a.size(); ++i) {
      N z = a[i];
      if (*mp) z = z * xp[*mp];
      mp++;
      if (*mp) z = z * yp[*mp];
      mp++;
      if (i == 0)
        y = z;
      else
        y = y + z;
    }

    return y;
  }

  N value(const N& x, const N& y) {
    N xy[] = {x, y};
    return value(xy);
  }

  N value(const std::initializer_list<N>& x) const {
    N xy[2] = {*(x.begin()), *(x.begin() + 1)};
    return value(xy);
  }

  N value(const PV2<N> &p) const {
    N xy[2] = {p.x, p.y};
    return value(xy);
  }

  static Poly2 one() {
    int nt = 1;
    N a[1];
    a[0] = N(1);
    int m[] = {0, 0};
    return Poly2(nt, a, m);
  }

  static Poly2 zero() {
    int nt = 1;
    N a[1];
    a[0] = N(0);
    int m[] = {0, 0};
    return Poly2(nt, a, m);
  }

  Poly<N> subX(N x) const {
    vector<N> coef(degy + 1);

    for (int i = 0; i < coef.size(); i++) {
      coef[i] = N(0);
    }

    vector<N> powx;
    powx.push_back(N(1));
    powx.push_back(x);
    for (int i = 2; i <= degx; i++) {
      powx.push_back(powx[powx.size() - 1] * x);
    }

    for (int i = 0; i < a.size(); i++) {
      int dx = m[2 * i];
      int dy = m[2 * i + 1];
      coef[dy] = coef[dy] + (a[i] * powx[dx]);
    }

    return Poly<N>(coef);
  }

  Poly<N> subY(N y) const {
    vector<N> coef(degx + 1);

    for (int i = 0; i < coef.size(); i++) {
      coef[i] = N(0);
    }

    vector<N> powy;
    powy.push_back(N(1));
    powy.push_back(y);
    for (int i = 2; i <= degy; i++) {
      powy.push_back(powy[powy.size() - 1] * y);
    }

    for (int i = 0; i < a.size(); i++) {
      int dx = m[2 * i];
      int dy = m[2 * i + 1];
      coef[dx] = coef[dx] + (a[i] * powy[dy]);
    }

    return Poly<N>(coef);
  }

  Poly2 derX() const {
    vector<N> coef;
    vector<int> pow;

    for (int i = 0; i < a.size(); i++) {
      if (m[2 * i] < 1) continue;
      coef.push_back(a[i] * m[2 * i]);
      pow.push_back(m[2 * i] - 1);
      pow.push_back(m[2 * i + 1]);
    }

    return Poly2(coef.size(), &coef[0], &pow[0]);
  }

  Poly2 derY() const {
    vector<N> coef;
    vector<int> pow;

    for (int i = 0; i < a.size(); i++) {
      if (m[2 * i + 1] < 1) continue;
      coef.push_back(a[i] * m[2 * i + 1]);
      pow.push_back(m[2 * i]);
      pow.push_back(m[2 * i + 1] - 1);
    }

    return Poly2(coef.size(), &coef[0], &pow[0]);
  }

  int degX() const {
    int d = 0;

    for (int i = 0; i < a.size(); i++) {
      int x = m[2 * i];
      d = std::max(d, x);
    }

    return d;
  }

  int degY() const {
    int d = 0;

    for (int i = 0; i < a.size(); i++) {
      int x = m[2 * i + 1];
      d = std::max(d, x);
    }

    return d;
  }

  void print_input() const {
    for (int i = 0; i < a.size(); i++) {
      printf("%.16f %d %d%s", a[i].mid(), m[2 * i], m[2 * i + 1],
             (i == a.size() - 1 ? "\n" : " "));
    }
  }

  void print() const {
    for (int i = 0; i < a.size(); i++) {
      printf("%.16f*x^(%d)*y^(%d) %s", a[i].mid(), m[2 * i], m[2 * i + 1],
             (i == a.size() - 1 ? "\n" : "+ "));
    }
  }

  PV2<N> polish(const PV2<N>& r) const {
    PV2<N> ret = r;
    return ret;
  }

  Poly2 operator*(const Poly2& b) const {
    int cnt = 0;
    // N ca[nt*b.nt];
    vector<N> ca(a.size() * b.a.size());
    int cm[2 * a.size() * b.a.size()];
    int exps[2];
    for (int i = 0; i < a.size(); i++) {
      auto mp = m.begin() + i * 2;
      for (int j = 0; j < b.a.size(); j++) {
        N coef = a[i] * b.a[j];
        auto bmp = b.m.begin() + j * 2;
        for (int k = 0; k < 2; k++) exps[k] = mp[k] + bmp[k];
        int l;
        for (l = 0; l < cnt; l++) {
          int* cmp = cm + l * 2;
          int h;
          for (h = 0; h < 2; h++)
            if (cmp[h] != exps[h]) break;
          if (h == 2) {
            ca[l] = ca[l] + coef;
            break;
          }
        }
        if (l == cnt) {
          int* cmp = cm + l * 2;
          for (int h = 0; h < 2; h++) cmp[h] = exps[h];
          ca[l] = coef;
          cnt++;
        }
      }
    }
    // return Poly2(cnt, ca, cm);
    return Poly2(cnt, &ca[0], cm);
  }

  Poly2 operator*(const N& p) const {
    vector<N> ca(a.size());
    int cm[2 * a.size()];

    for (int i = 0; i < a.size(); i++) {
      ca[i] = a[i] * p;
      cm[2 * i] = m[2 * i];
      cm[2 * i + 1] = m[2 * i + 1];
    }

    return Poly2(a.size(), &ca[0], cm);
  }

  Poly2 operator*(const double p) const {
    vector<N> ca(a.size());
    int cm[2 * a.size()];

    for (int i = 0; i < a.size(); i++) {
      ca[i] = a[i] * p;
      cm[2 * i] = m[2 * i];
      cm[2 * i + 1] = m[2 * i + 1];
    }

    return Poly2(a.size(), &ca[0], cm);
  }

  Poly2 plusCtimes(double c, const Poly2& b) const {
    int pnt = a.size();
    vector<N> pa(a.size() + b.a.size());
    int pm[2 * (a.size() + b.a.size())];

    for (int i = 0; i < a.size(); i++) pa[i] = a[i];
    int ntnv = a.size() * 2;
    for (int i = 0; i < ntnv; i++) pm[i] = m[i];

    auto bmp = b.m.begin();
    for (int i = 0; i < b.a.size(); i++) {
      int j;
      for (j = 0; j < a.size(); j++) {
        auto mp = m.begin() + j * 2;
        int k;
        for (k = 0; k < 2; k++)
          if (mp[k] != bmp[k]) break;
        if (k == 2) {
          pa[j] = pa[j] + c * b.a[i];
          break;
        }
      }
      if (j == a.size()) {
        pa[pnt] = c * b.a[i];
        int* pmp = pm + pnt * 2;
        for (int k = 0; k < 2; k++) pmp[k] = bmp[k];
        pnt++;
      }
      bmp += 2;
    }

    return Poly2(pnt, &pa[0], pm);
  }

  bool hasZeros () {
    for (N &ai : a)
      if (ai.sign() == 0)
	return true;
    return false;
  }

  Poly2 removeZeros () {
    int n = 0;
    for (int i = 0; i < a.size(); i++)
      if (a[i].sign() != 0) {
	if (n < i) {
	  a[n] = a[i];
	  m[2*n] = m[2*i];
	  m[2*n+1] = m[2*i+1];
	}
	++n;
      }
    a.resize(n);
    m.resize(2*n);
    setDegree();
    return *this;
  }

  Poly2 operator+(const Poly2& b) const { return plusCtimes(1.0, b); }
  Poly2 operator-(const Poly2& b) const { return plusCtimes(-1.0, b); }

  int size() const { return a.size(); }
  const N& operator[](int i) const { return a[i]; }
  N& operator[](int i) { return a[i]; }

  int d;
  vector<int> m;
  vector<N> a;
  int degx, degy;
};

}  // namespace acp
#endif
