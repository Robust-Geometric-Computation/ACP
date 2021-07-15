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

#ifndef POLY_H
#define POLY_H

#include <unordered_map>
#include "acp.h"
using namespace std;
using namespace acp;

namespace acp {

template<class T>
class Poly {
public:
  Poly () {}
  Poly (const T& c) { add(c, 0); }
  Poly (const vector<T> &v) {
    for (int i = 0; i < v.size(); i++)
      add(v[i], i);
  }

  template<class N>
  Poly (const Poly<N> &p) {
    for (int i = 0; i < p.a.size(); i++)
      a.push_back(p.a[i]);
  }
  
  unsigned int size () const { return a.size(); }
  const T &operator[] (unsigned int i) const { return a[i]; }
  T &operator[] (unsigned int i) { return a[i]; }
  T lc () const { return *a.rbegin(); }
  int deg () const { return a.size() - 1; }
  int degree () const { return deg(); }
  bool zero () const { return a.empty(); }

  void add (const T &c, int e) {
    int d = deg();
    if (e <= d)
      a[e] = a[e] + c;
    else {
      T zero = c;
      zero = 0;
      while (d + 1 < e) {
	a.push_back(zero);
	++d;
      }
      a.push_back(c);
    }
  }

  // make monic
  void monocize () {
    if (degree() < 0) {
      cerr << "Poly::monocize failed" << endl;
      exit(0);
    }
    T c = a.back().rcp();
    for (int i = 0; i < degree(); ++i)
      a[i] = a[i] * c;
    a[a.size() - 1] = 1;
  }

  Poly operator+ (const Poly &p) const {
    int d = deg(), e = p.deg();
    Poly q;
    if (d <= e) {
      for (int i = 0; i <= d; ++i)
	q.add(a[i] + p.a[i], i);
      for (int i = d + 1; i <= e; ++i)
	q.add(p.a[i], i);
    }
    else {
      for (int i = 0; i <= e; ++i)
	q.add(a[i] + p.a[i], i);
      for (int i = e + 1; i <= d; ++i)
	q.add(a[i], i);
    }
    return q;
  }

  Poly operator- (const Poly &p) const { return (*this) + (- p); }

  Poly operator- () const {
    Poly q;
    for (int i = 0; i <= deg(); ++i)
      q.add(- a[i], i);
    return q;
  }

  Poly operator* (const Poly &p) const {
    int d = deg(), e = p.deg();
    Poly q;
    for (int i = 0; i <= d; ++i)
      for (int j = 0; j <= e; ++j)
	q.add(a[i]*p.a[j], i + j);
    return q;
  }
  
  Poly plusCtimes(const T &c, const Poly &b) const {
  
    int d = a.size()-1;
    int bd = b.a.size()-1;
    int degree = max(d, bd);
    vector<T> coef(degree+1);
    for(int i = 0; i <= degree; i++) {
      coef[i]  = (d > bd ? a[i] : c * b.a[i]);
    }
    for(int i = 0; i <= min(d, bd); i++) {
      coef[i] = coef[i] + (d <= bd ? a[i] : c * b.a[i]);
    }
    return Poly(coef);
  }

  Poly rem (const Poly &p, Poly &q) const {
    T lp = p.lc();
    int dp = p.deg();
    Poly r(*this);
    while (dp <= r.deg()) {
      T k = r.lc()/lp;
      int dk = r.deg() - dp;
      q.add(k, dk);
      r.a.pop_back();
      for (int i = 0; i < dp; ++i)
	r.add(- k*p.a[i], i + dk);
      r.checkDegree();
    }

#ifdef RESIDUE_DEBUG
    Poly d = q * p + r - *this;
    d.checkDegree();
    if (d.deg() != -1) {
      cerr << "Poly::rem failed" << endl;
      exit(0);
    }
#endif

    return r;
  }

  Poly quotient (const Poly &p, Poly &r) const {
    Poly q;
    r = rem(p, q);
    return q;
  }

  void checkDegree () {
    while (!zero() && isZero(lc()))
      a.pop_back();
  }

  Poly gcd (const Poly &p, Poly &s, Poly &t) const
  {
    s = Poly(1);
    t = Poly(0);
    Poly r(*this), nr(p), ns, nt(T(1));
    r.checkDegree();
    nr.checkDegree();
    while (!nr.zero()) {
      Poly q, nnr = r.rem(nr, q), nns = s - q*ns, nnt = t - q*nt;
      r = nr;
      nr = nnr;
      s = ns;
      ns = nns;
      t = nt;
      nt = nnt;
    }
    return r;
  }

  Poly gcd (const Poly &p) const
  {
    Poly s, t;
    return gcd(p, s, t);
  }

  T value (const T &x) const
  {
    int d = deg();
    T y = a[d];
    for (int i = d - 1; i >= 0; --i)
      y = x*y + a[i];
    return y;
  }

  T der (const T &x) const
  {
    int d = deg();
    T y = a[d]*d;
    for (int i = d - 1; i > 0; --i)
      y = x*y + a[i]*i;
    return y;
  }

  Poly der () const
  {
    Poly g;
    for (int i = 1; i <= deg(); ++i)
      g.a.push_back(a[i]*i);
    return g;
  }
  
  Poly neg () const
  {
    Poly p;
    for (unsigned int i = 0u; i <= deg(); ++i)
      p.a.push_back(i%2 == 0 ? a[i] : - a[i]);
    return p;
  }

  Poly dual () const
  {
    Poly p;
    int d = deg();
    for (unsigned int i = 0u; i <= d; ++i)
      p.a.push_back(a[d-i]);
    return p;
  }

  Poly combine (const Poly &p, T s) const {
    Poly q;
    for (int i = 0; i <= deg(); ++i)
      q.add(p.a[i]*s - a[i]*(s - 1.0), i);
    return q;
  }

  T shrinkBracket (T x, int sub) const;

  T shrinkBracketB (const T &x, int sub) const;

  Poly moebius (const T &l, const T &u) const
  {
    return shift(l).dual().shift((u - l).rcp());
  }

  Poly shift (const T &s) const
  {
    Poly p = *this;
    if (s.sign() == 0)
      return p;
    int d = deg();
    for (int i = 0; i < d; ++i)
      for (int j = d - 1; j >= i; --j)
	p.a[j] = p.a[j] + s * p.a[j+1];
    return p;
  }

  vector<T> a;
};

vector<PTR<Object<Scalar>>> getRoots(PTR<Object<Poly>> poly);
vector<PTR<Object<Scalar>>> getRoots(PTR<Object<Poly>> poly,
                                     PTR<Object<Scalar>> lb,
                                     PTR<Object<Scalar>> ub);
}

#endif
