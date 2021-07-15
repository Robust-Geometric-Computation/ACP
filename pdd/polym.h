//  -*- Mode: c++ -*-
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

#ifndef POLYM_H
#define POLYM_H

#include "poly.h"

namespace acp {

// Multivariate polynomial
template<class N>
class PolyM {
public:
  int m;
  vector<N> a;

  // a[v2i[{4, 3, 7}]] is coefficient of x_0^4 x_1^3 x_2^2
  map<vector<int>, int> v2i;

  template<class T>
  PolyM (const PolyM<T> &p) : m(p.m), v2i(p.v2i) {
    for (int i = 0; i < p.a.size(); i++)
      a.push_back(p.a[i]);
  }
  unsigned int size () const { return a.size(); }
  const N &operator[] (unsigned int i) const { return a[i]; }
  N &operator[] (unsigned int i) { return a[i]; }

  // Create empty m-variate polynomial.
  PolyM (int m=0) : m(m) { if (m < 0) { cerr << "bad PolyM" << endl; exit(0); } }

  // Create 1-variate polynomial from univariate polynomial.
  template<class T>
  PolyM (const Poly<T> &p, int m=1, int index=0) : m(m) {
    if (index >= m) {
      cerr << "bad PolyM index" << endl;
      exit(0);
    }
    vector<int> v(m);
    for (int i = 0; i < p.a.size(); i++)
      if (!isZero(p.a[i])) {
        v[index] = i;
        add(v, p.a[i]);
      }
  }

  void add (const vector<int> &v, const N &c) {
    if (v.size() != m) {
      cerr << "bad PolyM::add" << endl;
      exit(0);
    }
    if (v2i.find(v) == v2i.end()) {
      v2i[v] = a.size();
      a.push_back(c);
    }
    else {
      int i = v2i[v];
      a[i] = a[i] + c;
    }
  }

  N value (vector<N> x) const {
    if (x.size() != m) {
      cerr << "bad PolyM::value" << endl;
      exit(0);
    }
    vector<vector<N>> pows;
    vector<N> one;
    one.push_back(N(1));
    for (int i = 0; i < m; ++i)
      pows.push_back(one);

    N v = 0;
    for (const pair<vector<int>, int> &vi : v2i) {
      const vector<int> &vec = vi.first;
      N c = a[vi.second];
      for (int i = 0; i < m; ++i) {
        while (pows[i].size() <= vec[i])
          pows[i].push_back(pows[i].back() * x[i]);
        c = c * pows[i][vec[i]];
      }
      v = v + c;
    }
    return v;
  }

  PolyM<N> operator+ (const PolyM<N> &b) const {
    if (m != b.m) {
      cerr << "bad PolyM::operator+" << endl;
      exit(0);
    }
    PolyM<N> p = *this;
    for (const pair<const vector<int>, int> vi : b.v2i)
      p.add(vi.first, b.a[vi.second]);
    return p;
  }

  PolyM<N> operator- (const PolyM<N> &b) const {
    if (m != b.m) {
      cerr << "bad PolyM::operator-" << endl;
      exit(0);
    }
    PolyM<N> p = *this;
    for (const pair<const vector<int>, int> vi : b.v2i)
      p.add(vi.first, -b.a[vi.second]);
    return p;
  }

  PolyM<N> operator- () const {
    PolyM<N> p = *this;
    for (int i = 0; i < p.a.size(); ++i)
      p.a[i] = -p.a[i];
    return p;
  }

  PolyM<N> operator* (const PolyM<N> &b) const {
    if (m != b.m) {
      cerr << "bad PolyM::operator*" << endl;
      exit(0);
    }
    PolyM<N> p(m);
    for (const pair<const vector<int>, int> v2iA : v2i)
      for (const pair<const vector<int>, int> v2iB : b.v2i) {
	vector<int> v = v2iA.first;
	for (int i = 0; i < v2iB.first.size(); ++i)
	  v[i] += v2iB.first[i];
	p.add(v, a[v2iA.second] * b.a[v2iB.second]);
      }
    return p;
  }

  PolyM<N> operator* (const N &k) const {
    PolyM<N> r(*this);
    for (int i = 0; i < a.size(); ++i)
      r.a[i] = r.a[i]*k;
    return r;
  }

  // Convert a 1-variate polynomial into a univariate polynomial
  Poly<N> toPoly () {
    if (m != 1) {
      cerr << "bad Pol<N> toPoly" << endl;
      exit(0);
    }
    Poly<N> p;
    for (const pair<vector<int>, int> &vi : v2i)
      p.add(a[vi.second], vi.first[0]);
    return p;
  }
};

class PolyM1Sign : public Primitive {
  Object<PolyM> *e1;
  Object<Scalar> *r0;
  Object<Poly> *m0;
  DeclareSign;
public:
  PolyM1Sign (Object<PolyM> *e1,
              Object<Scalar> *r0);
  PolyM1Sign (Object<PolyM> *e1,
              Object<Scalar> *r0,
              Object<Poly> *m0)
    : e1(e1), r0(r0), m0(m0) {}
};

class PolyM2Sign : public Primitive {
  Object<PolyM> *e2;
  Object<Scalar> *r0, *r1;
  Object<Poly> *m0, *m1;
  DeclareSign;
public:
  PolyM2Sign (Object<PolyM> *e2,
              Object<Scalar> *r0, Object<Scalar> *r1);
  PolyM2Sign (Object<PolyM> *e2,
              Object<Scalar> *r0, Object<Scalar> *r1,
              Object<Poly> *m0, Object<Poly> *m1)
    : e2(e2), r0(r0), r1(r1), m0(m0), m1(m1) {}
};

class PolyM3Sign : public Primitive {
  Object<PolyM> *e3;
  Object<Scalar> *r0, *r1, *r2;
  Object<Poly> *m0, *m1, *m2;
  DeclareSign;
public:
  PolyM3Sign (Object<PolyM> *e3,
              Object<Scalar> *r0, Object<Scalar> *r1, Object<Scalar> *r2);
};

class PolyM4Sign : public Primitive {
  Object<PolyM> *e4;
  Object<Scalar> *r0, *r1, *r2, *r3;
  Object<Poly> *m0, *m1, *m2, *m3;
  DeclareSign;
public:
  PolyM4Sign (Object<PolyM> *e4,
              Object<Scalar> *r0, Object<Scalar> *r1, Object<Scalar> *r2, Object<Scalar> *r3);
};

}
#endif

