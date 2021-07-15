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

class PParameter;
class MParameter;

class Parameter {
  friend class PParameter;
  friend class MParameter;
  double l, u;
public:
  Parameter () : l(0.0), u(1.0) {}
  Parameter (double l, double u) : l(l), u(- u) {}
  Parameter (double x) : l(x), u(- x) {}
  Parameter (double l, double u, int i) : l(l), u(u) {}
  Parameter (const Parameter &il, const Parameter &iu, bool alg)
      : l(il.l), u(iu.u) {}
  Parameter (const PParameter &e);
  Parameter (const MParameter &e);
  bool uninitialized () const { return l == 0.0 && u == 1.0; }
  double lb () const { return l; }
  double ub () const { return - u; }

  // newacp
  double mid () const { return (lb() + ub()) / 2; }

  int sign (bool fail = true) const {
    if (l > 0.0)
      return 1;
    if (ub() < 0.0)
      return -1;
    if (!fail || (l == 0.0 && u == 0.0))
      return 0;
    throw SignException();
  }

  Parameter operator+ (const Parameter &b) const {
    return Parameter(l + b.l, u + b.u, 0);
  }

  Parameter operator+ (double b) const { return *this + Parameter(b); }

  Parameter operator- () const {
    return Parameter(u, l, 0);
  }

  Parameter operator- (const Parameter &b) const {
    return Parameter(l + b.u, u + b.l, 0);
  }

  Parameter operator- (double b) const { return *this - Parameter(b); }

  Parameter operator* (const Parameter &b) const {
    Parameter s = u >= 0.0 ? - *this : *this, t = u >= 0.0 ? - b : b;
    if (s.l > 0.0) {
      if (t.l >= 0.0)
        return Parameter(s.l*t.l, - s.u*t.u, 0);
      else if (t.u >= 0.0)
        return Parameter(- s.u*t.l, s.l*t.u, 0);
      else
        return Parameter(- s.u*t.l, - s.u*t.u, 0);
    }
    return Parameter(min(- s.l*t.u, - s.u*t.l), min(- s.l*t.l, - s.u*t.u), 0);
  }

  Parameter operator* (double b) const { return *this*Parameter(b); }

  Parameter rcp () const {
    if (sign() == 0) {
      cerr << "reciprocal of zero in Parameter" << endl;
      exit(0);
    }
    return Parameter(-1.0/u, -1.0/l, 0);
  }

  Parameter operator/ (const Parameter &b) const { return *this*b.rcp(); }
  Parameter operator/ (double b) const { return *this*Parameter(b).rcp(); }
  Parameter abs () const { return sign() == 1 ? *this : - *this; }

  Parameter sqrt () const {
    int s = -1;
    try {
      s = sign();
    }
    catch (SignException se) {}
    if (s == -1)
      throw SignException(true);
    return s == 0 ? Parameter(0.0)
      : Parameter(Parameter::sqrt(l).l, Parameter::sqrt(ub()).ub());
  }
      
  Parameter sqrt (double x) const {
    double s = ::sqrt(x);
    Parameter xx(x, x);
    Parameter lr = xx/s;
    if (s < lr.l)
      return Parameter(s, lr.ub());
    if (s > lr.ub())
      return Parameter(lr.l, s);
    return lr;
  }

  Parameter lbP () const { return Parameter(l, l); }
  Parameter ubP () const { return Parameter(ub(), ub()); }
  Parameter midP () const {
    double m = 0.5*(l + ub());
    return Parameter(m, m);
  }

  bool subset (const Parameter &b) const {
    return !(lb() < b.lb()) && !(b.ub() < ub()) &&
           (b.lb() < lb() || ub() < b.ub());
  }

  bool intersects (const Parameter &b) const {
    return !(ub() < b.lb() || b.ub() < lb());
  }

  Parameter interval (const Parameter &b) const { return Parameter(l, b.ub()); }

  Parameter intersect (const Parameter &b) const {
    if (ub() < b.lb() || b.ub() < lb()) {
      cerr << "error in Parameter::intersect" << endl;
      exit(0);
    }
    double il = lb() < b.lb() ? b.lb() : lb();
    double iu = ub() < b.ub() ? ub() : b.ub();
    return Parameter(il, iu);
  }

  double distance (const Parameter &b) const {
    if (ub() <= b.lb())
      return b.lb() - ub();
    if (b.ub() <= lb())
      return lb() - b.ub();
    return lb() == b.lb() && ub() == b.ub() ? 0.0 : -1.0;
  }

  double intervalWidth () const { return ub() - lb(); }
  // polym support
  Parameter value () const { return *this; }
  Parameter make (const Parameter &x) const { return x; }
  PolyM<Parameter> PolyMNtoPolyMRes (const PolyM<Parameter> &p) const;
};

