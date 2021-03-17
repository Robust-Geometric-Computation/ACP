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

#ifndef ROOT_H
#define ROOT_H

#include "acp/core/object.h"
#include "acp/poly/poly.h"

namespace acp {

class SubKnot : public Object<Scalar> {
  PTR<Object<Scalar>> l, r;
  unsigned long int i, n;

  DeclareCalculate(Scalar) {
    return (l->get<N>().x * (n - i) + r->get<N>().x * i) / n;
  }

 public:
  SubKnot(Object<Scalar>* l, Object<Scalar>* r, unsigned long int i,
          unsigned long int n)
      : l(l), r(r), i(i), n(n) {}
};

class DInt {
 public:
  PTR<Object<Scalar>> l, u;
  unsigned long int n;

  DInt(Object<Scalar>* l, Object<Scalar>* u, unsigned long int n)
      : l(l), u(u), n(n) {}
};

class Root;

typedef vector<PTR<Object<Scalar>>> Roots;

class PolySolver {
  Roots linearRoots(Object<Scalar>* l, Object<Scalar>* u);
  Roots quadraticRoots();
  Roots quadraticRoots(Object<Scalar>* l, Object<Scalar>* u);
  Roots descartesRoots(Object<Scalar>* l, Object<Scalar>* u);
  bool descartes1(vector<DInt>& st, const DInt& i, int v, bool lflag);
  bool descartes2(vector<DInt>& st, const DInt& i, int v);
  unsigned long int descartes2int(const DInt& i, int v);
  PTR<Object<Poly>> poly;
  Poly<Parameter> get() { return poly->get<Parameter>(); }
  Poly<Parameter> getApprox(double e) { return poly->getApprox(e); }

  class CoeffSign : public Primitive {
    Object<Poly>* p;
    int i;

    DeclareSign { return p->get<N>().a[i]; }

   public:
    CoeffSign(Object<Poly>* p, int i) : p(p), i(i) {}
  };

 public:
  PolySolver(Object<Poly>* poly) : poly(poly) {}

  int deg() {
    int d = poly->getApprox(1.0).deg();
    while (d >= 0 && CoeffSign(poly, d) == 0) --d;
    return d;
  }

  Roots getRoots();
  Roots getRoots(PTR<Object<Scalar>> l, PTR<Object<Scalar>> u);
  Roots getRoots(double l, double u) {
    PTR<Object<Scalar>> pl = new Object<Scalar>(Parameter::constant(l), true);
    PTR<Object<Scalar>> pr = new Object<Scalar>(Parameter::constant(u), true);
    return getRoots(pl, pr);
  }
};

class Root : public Object<Scalar> {
 protected:
  PTR<Object<Poly>> p;

 public:
  Root(Object<Poly>* p) : p(p) {}
  PTR<Object<Poly>> getPoly() { return p; }
};

class LinearRoot : public Root {
  DeclareCalculate(Scalar) {
    Poly<N> q = p->get<N>();
    return -q[0] / q[1];
  }

 public:
  LinearRoot(Object<Poly>* p) : Root(p) {}
};

class QuadraticRoot : public Root {
  bool flag;

  DeclareCalculate(Scalar) {
    if (std::is_same<N, PParameter>::value) {
      if (!uninitialized()) {
        Parameter x = getCurrentP().x;
        return N(PParameter(x, true));
      }
      throw SignException(true);
    }
    Poly<N> q = p->get<N>();
    const N &a = q[2], &b = q[1], &c = q[0];
    N dd = (b * b - 4.0 * a * c);
    if (perI() != -1 && dd.sign() == -1) throw SignException(true);
    N d = dd.sqrt();
    if (b.sign() == -1) return flag ? 2.0 * c / (d - b) : 0.5 * (d - b) / a;
    return flag ? -0.5 * (d + b) / a : -2.0 * c / (d + b);
  }

 public:
  QuadraticRoot(Object<Poly>* p, bool flag) : Root(p), flag(flag) {}
};

class Descartes : public Object<Scalar> {
  Object<Poly>* p;
  Object<Scalar>*l, *u;

  DeclareCalculate(Scalar) {
    Poly<N> q = p->get<N>().moebius(l->get<N>().x, u->get<N>().x);
    int v = 0u, d = q.deg();
    int s = q.a[d].sign();
    for (int i = d - 1u; i >= 0; --i)
      if (q.a[i].sign() == -s) {
        ++v;
        s = -s;
      }
    return Scalar<N>(N(v));
  }

 public:
  Descartes(Object<Poly>* p, Object<Scalar>* l, Object<Scalar>* u)
      : p(p), l(l), u(u) {}
};

int descartes(Object<Poly>* p, Object<Scalar>* l, Object<Scalar>* u);

Scalar<Parameter> trace(const Poly<Parameter>& f, const Poly<Parameter>& g,
                        Parameter x);

class PolyRoot : public Root {
  PTR<Object<Scalar>> l, u;
  DeclareCalculate(Scalar) {
    if (std::is_same<N, PParameter>::value) {
      if (!uninitialized()) {
        Parameter x = getCurrentP().x;
        return N(PParameter(x, true));
      }
      throw SignException(true);
    }
    bool hom = perI() != -1;
    if (hom && Descartes(p, l, u).get<N>().x.lb() != 1)
      throw SignException(true);
    const Poly<N>& q = p->get<N>();
    int slb = q.value(l->get<N>().x).sign();
    int sub = q.value(u->get<N>().x).sign();
    N lu1(l->get<N>().x.ubP()), lu(lu1, lu1, true), ul1(u->get<N>().x.lbP()),
        ul(ul1, ul1, true);
    if (uninitialized() || hom) return q.shrinkBracket(N(lu, ul, true), sub);
    Parameter x = getCurrentP().x;
    N xl(x.lbP(), x.lbP(), true), xu(x.ubP(), x.ubP(), true);
    return q.shrinkBracket(N((xl - lu).sign(false) > 0 ? xl : lu,
                             (xu - ul).sign(false) < 0 ? xu : ul, true),
                           sub);
  }

 public:
  PolyRoot(Object<Poly>* p, Object<Scalar>* l, Object<Scalar>* u)
      : Root(p), l(l), u(u) {}
};

double polish(const Poly<Parameter>& f, double& x, int& k);

class CauchyBound : public Object<Scalar> {
  PTR<Object<Poly>> p;

  DeclareCalculate(Scalar) {
    Poly<N> q = p->get<N>();
    q.checkDegree();
    N b(1), qd = q.lc();
    for (int i = 0; i < q.deg(); ++i) {
      N bi = (q.a[i] / qd).abs();
      if (b < bi) b = bi;
    }
    return Scalar<N>(1.0 + b);
  }

 public:
  CauchyBound(Object<Poly>* p) : p(p) {}
};

class QuadraticRoots : public Primitive {
  Object<Poly>* p;

  DeclareSign {
    const Poly<N>& q = p->get<N>();
    return q[1] * q[1] - 4.0 * q[0] * q[2];
  }

 public:
  QuadraticRoots(Object<Poly>* p) : p(p) {}
};

class Order : public Primitive {
  Object<Scalar>*a, *b;

  DeclareSign { return b->get<N>().x - a->get<N>().x; }

 public:
  Order(Object<Scalar>* a, Object<Scalar>* b) : a(a), b(b) {}
};

template <class N>
bool zerop(const N& p) {
  return p.sign() == 0;
}

typedef PTR<Object<Scalar>> RootPTR;

}  // namespace acp
#endif
