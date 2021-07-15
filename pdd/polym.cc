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

#include "polym.h"
#include "root.h"

#ifdef GCD_STATS 
extern int nGCD53, nGCD100, nGCD106, nGCD212, nGCD424;
#endif

namespace acp {

void throwSignException () {
  throw SignException();
}

template<class R>
class ResidueException {
};

template<class R>
void throwResidueException () {
  throw ResidueException<R>();
}

template<class N>
ostream& operator<< (ostream &out, const Poly<N> &p) {
  out << "Poly[";
  for (const N &c : p.a)
    out << " " << c;
  out << "]";
  return out;
}

// a mod m.  m must be monic.  Eliminates leading zero coefficients.
template<class N>
Poly<N> mod (const Poly<N>& a, const Poly<N> &m) {
  Poly<N> r = a;
  while(r.a.size() > 0 && isZero(r.a.back()))
    r.a.pop_back();
  if (r.degree() < m.degree())
    return r;

#ifdef RESIDUE_DEBUG
  Poly<N> q;
#endif

  while (r.degree() >= m.degree()) {
    int dd = r.degree() - m.degree();
    N c = -r.a.back();
#ifdef RESIDUE_DEBUG
    q.add(-c, dd);
#endif
    r.a.pop_back();
    for (int i = 0; i < m.degree(); ++i)
      r.add(c * m.a[i], i + dd);
    while(r.a.size() > 0 && isZero(r.a.back()))
      r.a.pop_back();
  }
  
#ifdef RESIDUE_DEBUG
  Poly<N> d = q * m + r - a;
  d.checkDegree();
  if (d.a.size() != 0)
    cout << "d " << d << endl;
  if (d.a.size() != 0) {
    cerr << "d.a not empty" << endl;
    exit(0);
  }
#endif
  
  return r;
}

// a must be monic
template<class N>
Poly<N> gcd (Poly<N> a, Poly<N> b) {
  while (a.a.size() > 1) {
    Poly<N> c = mod(b, a);
    if (c.a.size() == 0)
      break;
    // Multiply c by the reciprocal of its leading coefficent.
    c.monocize();
    b = a;
    a = c;
  }
  return a;
}

template<class N>
Poly<N> squareFree (Poly<N> f) {
  if (!isZero(f.a.back() - 1))
    f.monocize();
  Poly<N> fder = f.der();
  Poly<N> g = gcd(f, fder);
  if (g.deg() == 0)
    return f;
  Poly<N> rem;
  return f.quotient(g, rem);
}

// inverse of a modulo m (actually m/this->gcd(m).monocize())
// r is set to this->gcd(m).monocize(), so if it is a constant,
// then it is 1 and the output is the inverse of a modulo m
template<class N>
Poly<N> inverse (const Poly<N> &a, const Poly<N> &m, Poly<N> &r) {
  r = m;
  Poly<N> t, nt, nr = a;
#ifdef RESIDUE_DEBUG
  Poly<N> u, nu;
#endif
  nt.a.push_back(m.a[0].make(1)); // nt = 1
#ifdef RESIDUE_DEBUG
  u.a.push_back(m.a[0].make(1)); // u = 1
#endif
  while (nr.a.size() > 0) {
    Poly<N> rem;
    Poly<N> q = r.quotient(nr, rem), pt = t, pr = r;
#ifdef RESIDUE_DEBUG
    Poly<N> pu = u;
#endif
    t = nt;
    nt = pt - q*nt;
#ifdef RESIDUE_DEBUG
    u = nu;
    nu = pu - q*nu;
#endif
    r = nr;
    nr = pr - q*nr;
    while (nr.a.size() > 0 && isZero(nr.a.back()))
      nr.a.pop_back();
  }

  N rinv = r.a.back().rcp();
  for (int i = 0; i < t.a.size(); ++i)
    t.a[i] = t.a[i] * rinv;
#ifdef RESIDUE_DEBUG
  for (int i = 0; i < u.a.size(); ++i)
    u.a[i] = u.a[i] * rinv;
#endif

#ifdef RESIDUE_DEBUG
  //Poly<N> d = t * *this + u * m;
  Poly<N> d = t * a + u * m;
  cout << "d " << d << endl;
#endif

  r.monocize();

#ifdef RESIDUE_DEBUG
  if (r.a.size() == 1) {
    Poly<N> p = a * t;
    Poly<N> o = mod(p, m);
    cout << "inverse " << o << endl;
  }
#endif

  return t;
}

class CmpVecInt {
public:
  // Compare based on highest index element that is different
  // {4, 3, 7} < {2, 6, 7} because 3 < 6
  int operator() (const vector<int> &a, const vector<int> &b) const {
    if (a.size() != b.size()) {
      cerr << "error in CmpVecInt" << endl;
      exit(0);
    }
    for (int i = a.size(); --i >= 0;)
      if (a[i] != b[i])
        return a[i] < b[i];
    return false;
  }

  bool equals (const vector<int> &a, const vector<int> &b) const {
    return !(operator()(a, b) || operator()(b, a));
  }
};

PolyM<Parameter> Parameter::PolyMNtoPolyMRes
(const PolyM<Parameter> &p) const { return p; }
PolyM<PParameter> PParameter::PolyMNtoPolyMRes
(const PolyM<PParameter> &p) const { return p; }
PolyM<MParameter> MParameter::PolyMNtoPolyMRes
(const PolyM<MParameter> &p) const { return p; }

template<class N>
ostream& operator<< (ostream &out, const PolyM<N> &p) {
  out << "PolyM";
  for (const pair<const vector<int>, int> &vi : p.v2i) {
    for (int i : vi.first)
      out << " " << i;
    out << " " << p.a[vi.second];
  }
  return out;
}

// Evaluates a polynomial on r using the value method of the
// R (residue) coefficients.
template<class N, class R>
N value (const Poly<R> &h, const N &r) {
  N v = 0;
  for (int i = h.a.size(); --i >= 0;)
    v = v * r + h.a[i].value();
  return v;
}

template<class N, class R>
class Residue;

template<class N, class R>
inline bool isZero (const Residue<N, R> &a) { return a.isZeror(); }

template<class N, class R>
ostream& operator<< (ostream &out, const Residue<N, R> &r);

// N is base numerical class (Parameter or /*M*/Parameter)
// R is (Residue) class of coefficients
// Defines element of R[x]/m(x) with m(r)=0.
template<class N, class R>
class Residue {
public:
  Poly<R> *m; // modulus--must be monic
  N *r;       // root defined by m
  Poly<R> p;  // polynomial mod m, unique representation of residue element.
  N v;	      // Value p(r) with R evaluated recursively
  
  // Creates illegal residue with no modulus or root when a
  // zero-argument constructor is needed.
  Residue () : m(0), r(0) {}
  
  // Creates zero element.
  Residue (Poly<R> &m, N &r) : m(&m), r(&r), v(N(0)) {}
  
  // p mod m
  Residue (Poly<R> &m, N &r, const Poly<R> &p) 
    : m(&m), r(&r), p(mod(p, m)), v(::value(p, r)) {}

  Residue<N, R> operator+ (const Residue<N, R> &b) const { return Residue<N, R>(*m, *r, p + b.p); }
  Residue<N, R> operator- () const { return Residue<N, R>(*m, *r, -p); }
  Residue<N, R> operator- (const Residue<N, R> &b) const { return Residue<N, R>(*m, *r, p - b.p); }
  Residue<N, R> operator* (const Residue<N, R> &b) const { return Residue<N, R>(*m, *r, p * b.p); }

  Residue<N, R> rcp () {
    Poly<R> g;
    Poly<R> i = inverse(p, *m, g);
    if (g.a.size() == 1) {
      Residue<N, R> inv = Residue<N, R>(*m, *r, i);
#ifdef RESIDUE_DEBUG
      Residue<N, R> d = inv * *this - make(1);
      if (!isZero(d)) {
	cerr << "isZero(d) false" << endl;
	exit(0);
      }
#endif
      return inv;
    }
    // *m = *m / g;
    Poly<R> rem;
    *m = m->quotient(g, rem);
    if (rem.degree() != -1) {
      cerr << "bad rem degree in Residue" << endl;
      exit(0);
    }
    throwResidueException<R>(); // new modulus
    return Residue();
  }

  // Precondition:  !isZero(d)
  Residue<N, R> operator/ (const Residue<N, R> &b) {
    Poly<R> g;
    Poly<R> i = inverse(b.p, *m, g);
    if (g.a.size() == 1) {
      Residue<N, R> ret = Residue<N, R>(*m, *r, p * i);
#ifdef RESIDUE_DEBUG
      Residue<N, R> d = ret * b - *this;
      if (!isZero(d)) {
	cerr << "isZero(d) false" << endl;
	exit(0);
      }
#endif
      return ret;
    }
    // *m = *m / g;
    Poly<R> rem;
    *m = m->quotient(g, rem);
    if (rem.degree() != -1) {
      cerr << "bad rem degree in Residue" << endl;
      exit(0);
    }
    throwResidueException<R>(); // new modulus
    return Residue();
  }

  N value () const { return v; }

  // Convert a number into an element of this residue.
  Residue<N, R> make (const N &x) const {
    Residue<N, R> res(*m, *r);
    if (!isZero(x)) {
      res.p.a.push_back(m->a[0].make(x));
      res.v = x;
    }
    return res;
  }
  
  void operator= (const N &x) {
    p.a.clear();
    if (!isZero(x)) {
      p.a.push_back(m->a[0].make(x));
    }
    v = x;
  }

  // Convert the coefficients of p into elements of this residue
  Poly<Residue<N, R>> PolyNtoPolyRes (const Poly<N> &p) const {
    Poly<Residue<N, R>> pRes;
    for (int i = 0; i < p.a.size(); ++i)
      pRes.a.push_back(make(p.a[i]));
    return pRes;
  }
  
  // Convert PolyM<R>(x1,x2,x3,...) into PolyM<Res>(x2,x3,...) by
  // considering each polynomial in x1 as an element of Res.
  PolyM<Residue<N, R>> PolyMRtoPolyMRes (const PolyM<R> &p) const {
    PolyM<Residue<N, R>> polyM(p.m - 1);
    CmpVecInt cmp;
    vector<int> prev;
    Poly<R> poly;
    for (const pair<vector<int>, int> &vi : p.v2i) {
      vector<int> exp = vi.first;
      exp.erase(exp.begin(), exp.begin()+1);
      if (prev.size() > 0 && !cmp.equals(prev, exp)) {
        polyM.add(prev, Residue<N, R>(*m, *r, poly));
        poly.a.clear();
      }
      poly.add(p.a[vi.second], vi.first[0]);
      prev = exp;
    }
    if (prev.size() > 0)
      polyM.add(prev, Residue<N, R>(*m, *r, poly));
#ifdef RESIDUE_DEBUG
    cout << "PolyMRtoPolyMRes:" << endl;
    cout << p << endl;
    cout << polyM << endl;
#endif
    return polyM;
  }

  // If this Residue is ((Q[x1]/f(1))[x2]/f(x2))[x3]/f(x3), this
  // method converts an element of Q[x1,x2,x3,x4,x5,...] into an
  // element of Residue[x4,x5,...].  Recursive.
  PolyM<Residue<N, R>> PolyMNtoPolyMRes (const PolyM<N> &p) const {
    return PolyMRtoPolyMRes(m->a[0].PolyMNtoPolyMRes(p));
  }

  bool isZeror () const {
    if (p.a.size() == 0) // zero p
      return true;
    if (p.a.size() == 1) {// constant p
      if (isZero(p.a[0])) {
	cerr << "bad isZeror 1" << endl;
	exit(0);
      }
      return false;
    }
    if (v.sign(false) != 0) // p(r) != 0
      return false;
    Poly<R> g = gcd(*m, p);
    if (g.a.size() == 1) // trivial gcd
      return false;
    // Modulus m of residue must be composite!
    // Poly<R> h = *m / g;
    Poly<R> rem;
    Poly<R> h = m->quotient(g, rem);
    if (rem.degree() != -1) {
      cerr << "bad isZeror 2" << endl;
      exit(0);
    }
    // If g(r) != 0, then h(r) must be zero.
    N gr = ::value(g, *r);
    N hr = ::value(h, *r);
    if (gr.sign(false) != 0) {
      h.monocize();
      *m = h; // new modulus is h made monic
    }
    // If h(r) != 0, then g(r) must be zero.
    else if (hr.sign(false) != 0)
      *m = g; // new moduls is g
    // Cannot tell which (g(r) or h(r)) is zero.
    else
      throwSignException(); // 
    throwResidueException<R>(); // new modulus
    return false; // Keep the compiler happy.
  }
};

template<class N>
using Residue0 = Residue<N, N>;

template<class N>
using Residue1 = Residue<N, Residue0<N>>;

template<class N>
using Residue2 = Residue<N, Residue1<N>>;

template<class N, class R>
ostream& operator<< (ostream &out, const Residue<N, R> &r) {
  return out << "Res " << r.p;
}

PolyM1Sign::PolyM1Sign (Object<PolyM> *e1,
                        Object<Scalar> *r0)
  : e1(e1), r0(r0),
    m0(dynamic_cast<Root*>(r0)->getPoly()) {}

int countM1;

template<class N>
N PolyM1Sign::calculate () {
  if (++countM1 == 0)
    cout << "this is it" << endl;
  PolyM<N> e = e1->get<N>();
  if (e.m != 1) {
    cerr << "bad PolyM1Sign::calculate" << endl;
    exit(0);
  }
  vector<N> r = { r0->get<N>().x };
  N er = e.value(r);
#ifdef RESIDUE_DEBUG
  cout << "e " << e << endl;
  cout << "r[0] " << r[0] << endl;
  cout << "er " << er << endl;
#endif
  int s = er.sign(false);
  if (s != 0)
    return er;
  
  Poly<N> f0 = squareFree(m0->get<N>());
#ifdef RESIDUE_DEBUG
  cout << "f0 " << f0 << endl;
#endif

  Poly<N> ePoly = e.toPoly();
#ifdef RESIDUE_DEBUG
  cout << "e " << e << endl;
  cout << "ePoly " << ePoly << endl;
#endif

  Poly<N> g;

  g = gcd(f0, ePoly);
#ifdef RESIDUE_DEBUG
  N gr = value(g, r[0]);
  cout << "gr " << gr << endl;
  cout << "g " << g << endl;
#endif
  if (g.a.size() == 1) // trivial gcd
    throwSignException(); // not zero, unknown sign
  // Poly<Residue0<N>> h = f1 / g;
  Poly<N> rem;
  Poly<N> h = f0.quotient(g, rem);
  // g=gcd(e,f0) and f0=g*h
  // if h(r0)!=0, then g(r0)=0 so e(r0)=0
  N v = value(h, r[0]);
  if (v.sign(false) != 0) {
#ifdef GCD_STATS 
    if (curPrecision() == 53u)
      nGCD53++;
    else if (curPrecision() == 100u)
      nGCD100++;
    else if (curPrecision() == 106u)
      nGCD106++;
    else if (curPrecision() == 212u)
      nGCD212++;
    else
      cout << "curPrecision " << curPrecision() << endl;
#endif
    return N(0);
  }
  throwSignException();
  return N(0);
}

template Parameter PolyM1Sign::calculate<Parameter>();
template PParameter PolyM1Sign::calculate<PParameter>();
template MParameter PolyM1Sign::calculate<MParameter>();

PolyM2Sign::PolyM2Sign (Object<PolyM> *e2,
                        Object<Scalar> *r0, Object<Scalar> *r1)
  : e2(e2), r0(r0), r1(r1),
    m0(dynamic_cast<Root*>(r0)->getPoly()),
    m1(dynamic_cast<Root*>(r1)->getPoly()) {}

template<class N>
N PolyM2Sign::calculate () {
  PolyM<N> e = e2->get<N>();
  if (e.m != 2) {
    cerr << "bad PolyM2Sign::calculate" << endl;
    exit(0);
  }
  vector<N> r = { r0->get<N>().x, r1->get<N>().x };
  N er = e.value(r);
#ifdef RESIDUE_DEBUG
  cout << "e " << e << endl;
  cout << "r[0] " << r[0] << endl;
  cout << "r[1] " << r[1] << endl;
  cout << "er " << er << endl;
#endif
  int s = er.sign(false);
  if (s != 0)
    return er;
  
  Poly<N> f0 = squareFree(m0->get<N>());
#ifdef RESIDUE_DEBUG
  cout << "f0 " << f0 << endl;
#endif
  Residue0<N> res0 = Residue<N,N>(f0, r[0]);
  Poly<N> f1y = squareFree(m1->get<N>());
  Poly<N> f1yder = f1y.der();
  Poly<N> f1ygcd = gcd(f1y, f1yder);
  if (f1ygcd.deg() > 0) {
    // cout << "f1y.deg() " << f1y.deg() << " f1ygcd.deg() " << f1ygcd.deg() << endl;
    Poly<N> rem;
    f1y = f1y.quotient(f1ygcd, rem);
    if (rem.deg() != -1) {
      cerr << "bad rem degree in PolyM2Sign::calculate" << endl;
      exit(0);
    }
  }
  Poly<Residue0<N>> f1 = res0.PolyNtoPolyRes(f1y);
#ifdef RESIDUE_DEBUG
  cout << "f1y " << f1y << endl;
  cout << "f1 " << f1 << endl;
  N f1r = value(f1, r[1]);
  cout << "f1r " << f1r << endl;
#endif

  Poly<Residue0<N>> ePoly;
  Poly<Residue0<N>> g;

  while (true) {
    try {
      ePoly = res0.PolyMNtoPolyMRes(e).toPoly();
#ifdef RESIDUE_DEBUG
      cout << "e " << e << endl;
      cout << "ePoly " << ePoly << endl;
#endif
      g = gcd(f1, ePoly);
#ifdef RESIDUE_DEBUG
      N gr = value(g, r[1]);
      cout << "gr " << gr << endl;
#endif
      break;
    } catch (ResidueException<N> re) {
#ifdef RESIDUE_DEBUG
      cout << "f0 " << f0 << endl;
#endif
      res0 = Residue<N,N>(f0, r[0]);
      f1 = res0.PolyNtoPolyRes(f1y);
#ifdef RESIDUE_DEBUG
      cout << "f1y " << f1y << endl;
      cout << "f1 " << f1 << endl;
#endif
    }
  }
#ifdef RESIDUE_DEBUG
  cout << "g " << g << endl;
#endif
  if (g.a.size() == 1) // trivial gcd
    throwSignException(); // not zero, unknown sign
  // Poly<Residue0<N>> h = f1 / g;
  Poly<Residue0<N>> rem;
  Poly<Residue0<N>> h = f1.quotient(g, rem);
  // g=gcd(e,f1) and f1=g*h
  // if h(r1)!=0, then g(r1)=0 so e(r1)=0
  N v = value(h, r[1]);
  if (v.sign(false) != 0) {
#ifdef GCD_STATS 
    if (curPrecision() == 53u)
      nGCD53++;
    else if (curPrecision() == 100u)
      nGCD100++;
    else if (curPrecision() == 106u)
      nGCD106++;
    else if (curPrecision() == 212u)
      nGCD212++;
    else if (curPrecision() == 424u)
      nGCD424++;
    else
      cout << "curPrecision " << curPrecision() << endl;
#endif
    return N(0);
  }
  throwSignException();
  return N(0);
}

template Parameter PolyM2Sign::calculate<Parameter>();
template PParameter PolyM2Sign::calculate<PParameter>();
template MParameter PolyM2Sign::calculate<MParameter>();

PolyM3Sign::PolyM3Sign (Object<PolyM> *e3,
                        Object<Scalar> *r0, Object<Scalar> *r1, Object<Scalar> *r2)
  : e3(e3), r0(r0), r1(r1), r2(r2),
    m0(dynamic_cast<Root*>(r0)->getPoly()),
    m1(dynamic_cast<Root*>(r1)->getPoly()),
    m2(dynamic_cast<Root*>(r2)->getPoly()) {}

template<class N>
N PolyM3Sign::calculate () {
  PolyM<N> e = e3->get<N>();
  if (e.m != 3) {
    cerr << "bad PolyM3Sign::calculate" << endl;
    exit(0);
  }
  vector<N> r = { r0->get<N>().x, r1->get<N>().x, r2->get<N>().x };
  N er = e.value(r);
#ifdef RESIDUE_DEBUG
  cout << "e " << e << endl;
  cout << "r[0] " << r[0] << endl;
  cout << "r[1] " << r[1] << endl;
  cout << "r[2] " << r[2] << endl;
  cout << "er " << er << endl;
#endif
  int s = er.sign(false);
  if (s != 0)
    return er;
  Poly<N> f0 = squareFree(m0->get<N>());
#ifdef RESIDUE_DEBUG
  cout << "f0 " << f0 << endl;
#endif
  Residue0<N> res0 = Residue<N,N>(f0, r[0]);
  Poly<N> f1y = squareFree(m1->get<N>());
  Poly<Residue0<N>> f1 = res0.PolyNtoPolyRes(f1y);
#ifdef RESIDUE_DEBUG
  cout << "f1y " << f1y << endl;
  cout << "f1 " << f1 << endl;
  N f1r = value(f1, r[1]);
  cout << "f1r " << f1r << endl;
#endif
  
  Residue1<N> res1 = Residue1<N>(f1, r[1]);
  Poly<N> f2z = squareFree(m2->get<N>());
  Poly<Residue1<N>> f2 = res1.PolyNtoPolyRes(f2z);
#ifdef RESIDUE_DEBUG
  cout << "f2z " << f2z << endl;
  cout << "f2 " << f2 << endl;
  N f2r = value(f2, r[2]);
  cout << "f2r " << f2r << endl;
#endif
  Poly<Residue1<N>> ePoly;
  Poly<Residue1<N>> g;
  while (true) {
    try {
      ePoly = res1.PolyMNtoPolyMRes(e).toPoly();
#ifdef RESIDUE_DEBUG
      cout << "e " << e << endl;
      cout << "ePoly " << ePoly << endl;
#endif
      g = gcd(f2, ePoly);
#ifdef RESIDUE_DEBUG
      N gr = value(g, r[2]);
      cout << "gr " << gr << endl;
#endif
      break;
    } catch (ResidueException<N> re) {
#ifdef RESIDUE_DEBUG
      cout << "f0 " << f0 << endl;
#endif
      res0 = Residue0<N>(f0, r[0]);
      f1 = res0.PolyNtoPolyRes(f1y);
#ifdef RESIDUE_DEBUG
      cout << "f1y " << f1y << endl;
      cout << "f1 " << f1 << endl;
      N f1r = value(f1, r[1]);
      cout << "f1r " << f1r << endl;
#endif
      
      res1 = Residue1<N>(f1, r[1]);
      f2 = res1.PolyNtoPolyRes(f2z);
#ifdef RESIDUE_DEBUG
      cout << "f2z " << f2z << endl;
      cout << "f2 " << f2 << endl;
      N f2r = value(f2, r[2]);
      cout << "f2r " << f2r << endl;
#endif
    }
    catch (ResidueException<Residue0<N>> re) {
#ifdef RESIDUE_DEBUG
      cout << "f1y " << f1y << endl;
      cout << "f1 " << f1 << endl;
      N f1r = value(f1, r[1]);
      cout << "f1r " << f1r << endl;
#endif
      res1 = Residue1<N>(f1, r[1]);
      f2 = res1.PolyNtoPolyRes(f2z);
#ifdef RESIDUE_DEBUG
      cout << "f2z " << f2z << endl;
      cout << "f2 " << f2 << endl;
      N f2r = value(f2, r[2]);
      cout << "f2r " << f2r << endl;
#endif
    }
  }
#ifdef RESIDUE_DEBUG
  acp::disable();
  cout << "g " << g << endl;
  acp::enable();
#endif
  if (g.a.size() == 1) // trivial gcd
    throwSignException(); // not zero, unknown sign
  // Poly<Residue0<N>> h = f1 / g;
  Poly<Residue1<N>> rem;
  Poly<Residue1<N>> h = f2.quotient(g, rem);
  // g=gcd(e,f2) and f2=g*h
  // if h(r2)!=0, then g(r2)=0 so e(r2)=0
  N v = value(h, r[2]);
  if (v.sign(false) != 0)
    return N(0);
  throwSignException();
  return N(0);
}

template Parameter PolyM3Sign::calculate<Parameter>();
template PParameter PolyM3Sign::calculate<PParameter>();
template MParameter PolyM3Sign::calculate<MParameter>();


PolyM4Sign::PolyM4Sign (Object<PolyM> *e4,
                        Object<Scalar> *r0, Object<Scalar> *r1, Object<Scalar> *r2, Object<Scalar> *r3)
  : e4(e4), r0(r0), r1(r1), r2(r2), r3(r3),
    m0(dynamic_cast<Root*>(r0)->getPoly()),
    m1(dynamic_cast<Root*>(r1)->getPoly()),
    m2(dynamic_cast<Root*>(r2)->getPoly()),
    m3(dynamic_cast<Root*>(r3)->getPoly()) {}

template<class N>
N PolyM4Sign::calculate () {
  PolyM<N> e = e4->get<N>();
  if (e.m != 4) {
    cerr << "bad PolyM4Sign::calculate" << endl;
    exit(0);
  }
  vector<N> r = { r0->get<N>().x, r1->get<N>().x, r2->get<N>().x, r3->get<N>().x };
  N er = e.value(r);
#ifdef RESIDUE_DEBUG
  cout << "e " << e << endl;
  cout << "r[0] " << r[0] << endl;
  cout << "r[1] " << r[1] << endl;
  cout << "r[2] " << r[2] << endl;
  cout << "r[3] " << r[3] << endl;
  cout << "er " << er << endl;
#endif
  int s = er.sign(false);
  if (s != 0)
    return er;
  
  Poly<N> f0 = squareFree(m0->get<N>());
#ifdef RESIDUE_DEBUG
  cout << "f0 " << f0 << endl;
#endif
  Residue0<N> res0 = Residue<N,N>(f0, r[0]);
  Poly<N> f1y = squareFree(m1->get<N>());
  Poly<Residue0<N>> f1 = res0.PolyNtoPolyRes(f1y);
#ifdef RESIDUE_DEBUG
  cout << "f1y " << f1y << endl;
  cout << "f1 " << f1 << endl;
  N f1r = value(f1, r[1]);
  cout << "f1r " << f1r << endl;
#endif
  
  Residue1<N> res1 = Residue1<N>(f1, r[1]);
  Poly<N> f2z = squareFree(m2->get<N>());
  Poly<Residue1<N>> f2 = res1.PolyNtoPolyRes(f2z);
#ifdef RESIDUE_DEBUG
  cout << "f2z " << f2z << endl;
  cout << "f2 " << f2 << endl;
  N f2r = value(f2, r[2]);
  cout << "f2r " << f2r << endl;
#endif

  Residue2<N> res2 = Residue2<N>(f2, r[2]);
  Poly<N> f3w = squareFree(m3->get<N>());
  Poly<Residue2<N>> f3 = res2.PolyNtoPolyRes(f3w);
#ifdef RESIDUE_DEBUG
  cout << "f3w " << f3w << endl;
  cout << "f3 " << f3 << endl;
  N f3r = value(f3, r[3]);
  cout << "f3r " << f3r << endl;
#endif

  Poly<Residue2<N>> ePoly;
  Poly<Residue2<N>> g;

  while (true) {
    try {
      ePoly = res2.PolyMNtoPolyMRes(e).toPoly();
#ifdef RESIDUE_DEBUG
      cout << "e " << e << endl;
      cout << "ePoly " << ePoly << endl;
#endif
      g = gcd(f3, ePoly);
#ifdef RESIDUE_DEBUG
      N gr = value(g, r[3]);
      cout << "gr " << gr << endl;
#endif
      break;
    } catch (ResidueException<N> re) {
#ifdef RESIDUE_DEBUG
      cout << "f0 " << f0 << endl;
#endif
      res0 = Residue0<N>(f0, r[0]);
      f1 = res0.PolyNtoPolyRes(f1y);
#ifdef RESIDUE_DEBUG
      cout << "f1y " << f1y << endl;
      cout << "f1 " << f1 << endl;
      N f1r = value(f1, r[1]);
      cout << "f1r " << f1r << endl;
#endif
      
      res1 = Residue1<N>(f1, r[1]);
      f2 = res1.PolyNtoPolyRes(f2z);
#ifdef RESIDUE_DEBUG
      cout << "f2z " << f2z << endl;
      cout << "f2 " << f2 << endl;
      N f2r = value(f2, r[2]);
      cout << "f2r " << f2r << endl;
#endif
      
      res2 = Residue2<N>(f2, r[2]);
      f3 = res2.PolyNtoPolyRes(f3w);
#ifdef RESIDUE_DEBUG
      cout << "f3w " << f3w << endl;
      cout << "f3 " << f3 << endl;
      N f3r = value(f3, r[3]);
      cout << "f3r " << f3r << endl;
#endif
    }
    catch (ResidueException<Residue0<N>> re) {
#ifdef RESIDUE_DEBUG
      cout << "f1y " << f1y << endl;
      cout << "f1 " << f1 << endl;
      N f1r = value(f1, r[1]);
      cout << "f1r " << f1r << endl;
#endif

      res1 = Residue1<N>(f1, r[1]);
      f2 = res1.PolyNtoPolyRes(f2z);
#ifdef RESIDUE_DEBUG
      cout << "f2z " << f2z << endl;
      cout << "f2 " << f2 << endl;
      N f2r = value(f2, r[2]);
      cout << "f2r " << f2r << endl;
#endif
      
      res2 = Residue2<N>(f2, r[2]);
      f3 = res2.PolyNtoPolyRes(f3w);
#ifdef RESIDUE_DEBUG
      cout << "f3w " << f3w << endl;
      cout << "f3 " << f3 << endl;
      N f3r = value(f3, r[3]);
      cout << "f3r " << f3r << endl;
#endif
    }
    catch (ResidueException<Residue1<N>> re) {
#ifdef RESIDUE_DEBUG
      cout << "f2z " << f2z << endl;
      cout << "f2 " << f2 << endl;
      N f2r = value(f2, r[2]);
      cout << "f2r " << f2r << endl;
#endif
      
      res2 = Residue2<N>(f2, r[2]);
      f3 = res2.PolyNtoPolyRes(f3w);
#ifdef RESIDUE_DEBUG
      cout << "f3w " << f3w << endl;
      cout << "f3 " << f3 << endl;
      N f3r = value(f3, r[3]);
      cout << "f3r " << f3r << endl;
#endif
    }
  }
#ifdef RESIDUE_DEBUG
  acp::disable();
  cout << "g " << g << endl;
  acp::enable();
#endif
  if (g.a.size() == 1) // trivial gcd
    throwSignException(); // not zero, unknown sign
  Poly<Residue2<N>> rem;
  Poly<Residue2<N>> h = f3.quotient(g, rem);
  // g=gcd(e,f3) and f3=g*h
  // if h(r3)!=0, then g(r3)=0 so e(r3)=0
  N v = value(h, r[3]);
  if (v.sign(false) != 0) {
#ifdef RESIDUE_DEBUG
    cout << "sign 0" << endl;
#endif
    return N(0);
  }
  throwSignException();
  return N(0);
}

template Parameter PolyM4Sign::calculate<Parameter>();
template PParameter PolyM4Sign::calculate<PParameter>();
template MParameter PolyM4Sign::calculate<MParameter>();
}
