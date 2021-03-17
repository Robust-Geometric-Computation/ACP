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

#include "acp/poly/root.h"
using namespace std;
using namespace acp;

namespace acp {

Roots PolySolver::getRoots() {
  int d = deg();
  if (d < 1) return Roots();
  if (d == 1) {
    Roots res;
    res.push_back(new LinearRoot(poly));
    return res;
  }
  if (d == 2) return quadraticRoots();
  double b = CauchyBound(poly).getApprox(1.0).x.ub();
  double lb = -randomNumber(1.0, 1.001) * b, ub = randomNumber(1.0, 1.001) * b;
  return getRoots(new Object<Scalar>(Parameter::constant(lb), true),
                  new Object<Scalar>(Parameter::constant(ub), true));
}

Roots PolySolver::getRoots(PTR<Object<Scalar>> l, PTR<Object<Scalar>> u) {
  int d = deg();
  if (d < 1) return Roots();
  if (d == 1) return linearRoots(l, u);
  if (d == 2) return quadraticRoots(l, u);
  return descartesRoots(l, u);
}

Roots PolySolver::linearRoots(Object<Scalar>* l, Object<Scalar>* u) {
  Roots res;
  RootPTR r = new LinearRoot(poly);
  if (Order(l, r) == 1 && Order(r, u) == 1) res.push_back(r);
  return res;
}

Roots PolySolver::quadraticRoots() {
  Roots res;
  if (QuadraticRoots(poly) == 1) {
    RootPTR r1 = new QuadraticRoot(poly, false),
            r2 = new QuadraticRoot(poly, true);
    if (Order(r1, r2) == 1) {
      res.push_back(r1);
      res.push_back(r2);
    } else {
      res.push_back(r2);
      res.push_back(r1);
    }
  }
  return res;
}

Roots PolySolver::quadraticRoots(Object<Scalar>* l, Object<Scalar>* u) {
  Roots res;
  if (QuadraticRoots(poly) == 1) {
    RootPTR r1 = new QuadraticRoot(poly, false),
            r2 = new QuadraticRoot(poly, true);
    if (Order(l, r1) == 1 && Order(r1, u) == 1) res.push_back(r1);
    if (Order(l, r2) == 1 && Order(r2, u) == 1) res.push_back(r2);
    if (res.size() == 2 && Order(res[0], res[1]) == -1)
      reverse(res.begin(), res.end());
  }
  return res;
}

Roots PolySolver::descartesRoots(Object<Scalar>* l, Object<Scalar>* u) {
  Roots res;
  vector<DInt> st;
  st.push_back(DInt(l, u, 4));
  int k = 0;
  while (!st.empty()) {
    ++k;
    DInt i = *st.rbegin();
    st.pop_back();
    int v = descartes(poly, i.l, i.u);
    if (v == 1)
      res.push_back(new PolyRoot(poly, i.l, i.u));
    else if (v > 1  // && !descartes1(st, i, v, true) && !descartes1(st, i, v,
                    // false)
             && !descartes2(st, i, v)) {
      PTR<Object<Scalar>> m = new SubKnot(i.l, i.u, 1, 2);
      unsigned long int n = max(4u, (unsigned int)sqrt(i.n));
      st.push_back(DInt(i.l, m, n));
      st.push_back(DInt(m, i.u, n));
    }
  }
  // cerr << "descartes iterations = " << k << endl;
  reverse(res.begin(), res.end());
  return res;
}

bool PolySolver::descartes1(vector<DInt>& st, const DInt& i, int v,
                            bool lflag) {
  PTR<Object<Scalar>> m = new SubKnot(i.l, i.u, lflag ? 1 : i.n - 1, i.n),
                      nl = lflag ? i.l : m, nu = lflag ? m : i.u;
  if (descartes(poly, nl, nu) != v) return false;
  unsigned long int ln = i.n, nn = min(ln * ln, 4294967295ul);
  st.push_back(DInt(nl, nu, nn));
  return true;
}

bool PolySolver::descartes2(vector<DInt>& st, const DInt& i, int v) {
  unsigned long int k = descartes2int(i, v);
  if (k == 0) return false;
  PTR<Object<Scalar>> nl = new SubKnot(i.l, i.u, k - 2, 4 * i.n),
                      nu = new SubKnot(i.l, i.u, k + 2, 4 * i.n);
  if (descartes(poly, nl, nu) != v) return false;
  unsigned long int ln = i.n, nn = min(ln * ln, 4294967295ul);
  st.push_back(DInt(nl, nu, nn));
  return true;
}

unsigned long int PolySolver::descartes2int(const DInt& i, int v) {
  Poly<Parameter> p = getApprox(1.0);
  Parameter il = i.l->getApprox(1.0).x, iu = i.u->getApprox(1.0).x,
            x = 0.5 * (il + iu), dp = p.der(x);
  if (dp.sign(false) == 0) return 0;
  Parameter nx = x - v * p.value(x) / dp, w = (nx - il) / (iu - il);
  return min(max(4 * i.n * w.mid(), 2.0), 4 * i.n - 2.0);
}

int descartes(Object<Poly>* p, Object<Scalar>* l, Object<Scalar>* u) {
  return Descartes(p, l, u).getApprox(1.0).x.lb();
}

Scalar<Parameter> trace(const Poly<Parameter>& f, const Poly<Parameter>& g,
                        Parameter x) {
  double s = 0.0, y = x.mid(), ds = 0.1, dsmin = 1e-10, emax = 1e-14;
  Poly<Parameter> r = f;
  while (s < 1.0) {
    if (ds < dsmin) throw SignException();
    double dy = r.der(Parameter::constant(y)).mid(),
           d = 1.0 / sqrt(1.0 + dy * dy), t = min(ds, (1.0 - s) * d),
           ns = ds == t ? s + ds : 1.0, ny = y + t * dy * d;
    r = ns == 1.0 ? g : f.combine(g, Parameter::constant(ns));
    double dr = fabs(r.lc().mid());
    int k;
    if (polish(r, ny, k) < dr * emax) {
      s = ns;
      y = ny;
      if (s == 1.0) break;
      if (k < 3) ds *= 2.0;
    } else if (fabs(dy) > dr * 1000)
      throw SignException();
    else
      ds *= 0.5;
  }
  return Scalar<Parameter>(Parameter::constant(y));
}

double polish(const Poly<Parameter>& f, double& x, int& k) {
  double v = f.value(Parameter::constant(x)).mid(), e = fabs(v);
  k = 0;
  while (e > 0.0) {
    ++k;
    double nx = x - v / f.der(Parameter::constant(x)).mid(),
           nv = f.value(Parameter::constant(nx)).mid();
    if (fabs(nv) > 0.5 * e) break;
    x = nx;
    v = nv;
    e = fabs(nv);
  }
  return e;
}

bool zerop(const Parameter& p) { return p.sign() == 0; }

}  // namespace acp

// debug

void pp(Object<Scalar>* p) { cerr << p->getApprox().x.mid() << endl; }

void pp(Object<Poly>* p) { p->getApprox(1.0).print(); }
