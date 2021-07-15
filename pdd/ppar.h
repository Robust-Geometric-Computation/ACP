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

class PParameter {
  friend class MParameter;
  Mods mods;
  bool alg;

public:
  static unsigned long int nModSign;
  static unsigned int nMixed, nRascal;
  Parameter p;

  PParameter () {}
  PParameter (double a) : p(a), mods(a), alg(false) {}
  PParameter (const Parameter &p, bool alg = false) : p(p), mods(p), alg(alg) {}

  PParameter (const Parameter &p, const Mods &mods, bool alg)
      : p(p), mods(mods), alg(alg) {
    if (sign(false)) ++nModSign;
  }

  PParameter (const MParameter &p);

  bool hasCurrentPrimes () { return mods.hasCurrentPrimes(); }

  PParameter operator+ (const PParameter &b) const {
    return PParameter(p + b.p, mods + b.mods, alg || b.alg);
  }

  PParameter operator+ (double b) const { return *this + PParameter(b); }

  PParameter operator- (const PParameter &b) const {
    return PParameter(p - b.p, mods - b.mods, alg || b.alg);
  }

  PParameter operator- (double b) const { return *this - PParameter(b); }

  PParameter operator- () const { return PParameter(-p, - mods, alg); }

  PParameter operator* (const PParameter &b) const {
    return PParameter(p*b.p, mods*b.mods, alg || b.alg);
  }

  PParameter operator* (double b) const { return *this*PParameter(b); }

  PParameter rcp () const { return PParameter(1) / *this; }

  PParameter operator/ (const PParameter &b) const {
    if (b.sign() == 0) {
      cerr << "reciprocal of zero in PParameter" << endl;
      exit(0);
    }
    try {
      return PParameter(p/b.p, mods/b.mods, alg || b.alg);
    }
    catch (SignException se) {
      throw SignException(alg || b.alg);
    }
  }

  PParameter operator/ (double b) const { return *this/PParameter(b); }
  PParameter abs() const { return sign() < 0 ? - *this : *this; }
  PParameter sqrt () const { return PParameter(p.sqrt(), true); }
  double lb () const { return mods.zero() ? 0.0 : p.lb(); }
  double ub () const { return mods.zero() ? 0.0 : p.ub(); }

  int sign (bool fail = true) const {
    if (p.lb() > 0.0 || p.ub() < 0.0) {
      if (mods.mixed()) {
        ++nMixed;
        for (int i = 0; i < NMods; ++i)
          if (mods.mod[i].a == 0)
	    throw mods.mod[i].p;
      }
      else if (mods.zero())
        ++nRascal;
      return p.lb() > 0.0 ? 1 : -1;
    }
    if (!fail || (p.lb() == 0 && p.ub() == 0) || mods.zero())
      return 0;
    throw SignException(alg);
  }

  PParameter (const Parameter &il, const Parameter &iu, bool alg) {
    cerr << "undefined in PParameter" << endl; exit(0);
  }
  PParameter lbP () const { cerr << "undefined in PParameter" << endl; exit(0); }
  PParameter ubP () const { cerr << "undefined in PParameter" << endl; exit(0); }
  PParameter midP () const { cerr << "undefined in PParameter" << endl; exit(0); }
  bool subset (const PParameter &b) const {
    cerr << "undefined in PParameter" << endl; exit(0);
  }
  bool intersects (const PParameter &b) const {
    cerr << "undefined in PParameter" << endl; exit(0);
  }
  PParameter interval (const PParameter &b) const {
    cerr << "undefined in PParameter" << endl; exit(0);
  }
  PParameter intersect (const PParameter &b) const {
    cerr << "undefined in PParameter" << endl; exit(0);
  }
  double intervalWidth () const { return p.intervalWidth(); }
  // polym support
  PParameter value () const { return *this; }
  PParameter make (const PParameter &x) const { return x; }
  PolyM<PParameter> PolyMNtoPolyMRes (const PolyM<PParameter> &p) const;
};

void pparReport ();

void resetPpar ();

