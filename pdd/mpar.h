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

class MValue {
  friend class Mods;
  mpfr_t m;
  unsigned int p;

  MValue (unsigned int p) : p(p) { mpfr_init2(m, p); }

public:
  MValue () : p(0) {}

  MValue (double x, unsigned int p) : p(p) {
    mpfr_init2(m, p);
    mpfr_set_d(m, x, GMP_RNDN);
  }

  void clear () { mpfr_clear(m); }

  MValue copy () const {
    MValue v(p);
    mpfr_set(v.m, m, GMP_RNDN);
    return v;
  }

  unsigned int getP () const { return p; }
  double value (mpfr_rnd_t round = GMP_RNDN) const {
    return mpfr_get_d(m, round);
  }

  MValue plus (const MValue &b, mpfr_rnd_t round) const {
    MValue res(max(p, b.p));
    mpfr_add(res.m, m, b.m, round);
    return res;
  }

  MValue minus () const {
    MValue res(p);
    mpfr_neg(res.m, m, GMP_RNDN);
    return res;
  }

  MValue minus (const MValue &b, mpfr_rnd_t round) const {
    MValue res(max(p, b.p));
    mpfr_sub(res.m, m, b.m, round);
    return res;
  }

  MValue times (const MValue &b, mpfr_rnd_t round) const {
    MValue res(max(p, b.p));
    mpfr_mul(res.m, m, b.m, round);
    return res;
  }

  MValue rcp (mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_d_div(res.m, 1.0, m, round);
    return res;
  }

  MValue sqrt (mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_sqrt(res.m, m, round);
    return res;
  }

  MValue mid (const MValue &b) const {
    MValue res = plus(b, GMP_RNDN);
    mpfr_mul_d(res.m, res.m, 0.5, GMP_RNDN);
    return res;
  }

  int sign () const { return mpfr_sgn(m); }
  bool operator== (const MValue &b) const { return mpfr_equal_p(m, b.m) != 0; }
  bool operator!= (const MValue &b) const { return mpfr_equal_p(m, b.m) == 0; }
  bool operator< (const MValue &b) const { return mpfr_less_p(m, b.m); }
};

class MParameter {
  friend class PParameter;
  template <template <typename> class P>
  friend class Object;
  MValue l, u;
  Mods mods;
  bool alg, clear;

  MParameter (const MValue &l, const MValue &u, const Mods &mods, bool alg)
      : l(l), u(u), mods(mods), alg(alg), clear(true) {}

public:
  static unsigned int nMixed, nRascal;

  MParameter () : clear(false) {}

  MParameter (double a)
    : l(a, curPrecision()), u(a, curPrecision()), mods(a), alg(false),
      clear(true) {}

  MParameter (const Parameter &a)
    : l(a.lb(), curPrecision()), u(a.ub(), curPrecision()), mods(a), alg(false),
      clear(true) {}
  
#ifdef EXACTRAT
  MParameter (const PParameter &a)
    : l(a.lb(), curPrecision()), u(a.ub(), curPrecision()), mods(a.lb()), alg(false),
      clear(true) {}
#else
  MParameter (const PParameter &a)
    : l(a.lb(), curPrecision()), u(a.ub(), curPrecision()), mods(a.mods), alg(a.alg),
      clear(true) {}
#endif
  MParameter (const MParameter &il, const MParameter &iu, bool ialg)
    : l(il.l.copy()), u(iu.u.copy()), mods(il.l, iu.u), alg(il.alg || iu.alg || ialg),
      clear(true) {
  }

  MParameter (const MParameter &a) : mods(a.mods), alg(a.alg) {
    if (!a.clear) {
      l = a.l;
      u = a.u;
      clear = false;
    }
    else {
      l = a.l.copy();
      u = a.u.copy();
      clear = true;
    }
  }

  ~MParameter () {
    if (clear) {
      l.clear();
      u.clear();
    }
  }

  bool hasCurrentPrimes () { return mods.hasCurrentPrimes(); }

  MParameter & operator= (const MParameter &a) {
    if (clear) {
      l.clear();
      u.clear();
    }
    if (!a.clear) {
      l = a.l;
      u = a.u;
      clear = false;
    }
    else {
      l = a.l.copy();
      u = a.u.copy();
      clear = true;
    }
    mods = a.mods;
    alg = a.alg;
    return *this;
  }

  MParameter operator+ (const MParameter &b) const {
    return MParameter(l.plus(b.l, GMP_RNDD), u.plus(b.u, GMP_RNDU),
                      mods + b.mods, alg || b.alg);
  }

  MParameter operator+ (double b) const { return *this + MParameter(b); }

  MParameter operator- (const MParameter &b) const {
    return MParameter(l.minus(b.u, GMP_RNDD), u.minus(b.l, GMP_RNDU),
                      mods - b.mods, alg || b.alg);
  }

  MParameter operator- (double b) const { return *this - MParameter(b); }

  MParameter operator- () const {
    return MParameter(u.minus(), l.minus(), - mods, alg);
  }

  MParameter operator* (const MParameter &b) const {
    return u.sign() == -1 ? (- *this).times(- b) : times(b);
  }

  MParameter times (const MParameter &b) const {
    bool nalg = alg || b.alg;
    Mods nmods = mods*b.mods;
    if (l.sign() == 1)
      return MParameter((b.l.sign() == 1 ? l : u).times(b.l, GMP_RNDD),
                        (b.u.sign() == 1 ? u : l).times(b.u, GMP_RNDU), nmods,
                        nalg);
    if (b.l.sign() == 1)
      return MParameter(l.times(b.u, GMP_RNDD), u.times(b.u, GMP_RNDU), nmods,
                        nalg);
    if (b.u.sign() == -1)
      return MParameter(u.times(b.l, GMP_RNDD), l.times(b.l, GMP_RNDU), nmods,
                        nalg);
    MValue cl1 = l.times(b.u, GMP_RNDD), cl2 = u.times(b.l, GMP_RNDD),
           cu1 = l.times(b.l, GMP_RNDU), cu2 = u.times(b.u, GMP_RNDU);
    bool lflag = cl1 < cl2, uflag = cu1 < cu2;
    MValue cl = lflag ? cl1 : cl2, cu = uflag ? cu2 : cu1;
    (lflag ? cl2 : cl1).clear();
    (uflag ? cu1 : cu2).clear();
    return MParameter(cl, cu, nmods, nalg);
  }

  MParameter operator* (double b) const { return *this*MParameter(b); }

  MParameter rcp () const {
    if (sign() == 0) {
      cerr << "reciprocal of zero in MParameter" << endl;
      exit(0);
    }
    return MParameter(u.rcp(GMP_RNDD), l.rcp(GMP_RNDU), Mods(1)/mods, alg);
  }

  MParameter operator/ (const MParameter &b) const { return *this*b.rcp(); }

  MParameter operator/ (double b) const { return *this*MParameter(b).rcp(); }

  MParameter abs () const { return sign() == 1 ? *this : - *this; }

  MParameter sqrt () const {
    if (sign() == -1) {
      cerr << "square root of negative number in MParameter" << endl;
      exit(0);
    }
    MValue nl = l.sqrt(MPFR_RNDD), nu = u.sqrt(MPFR_RNDU);
    return MParameter(nl, nu, Mods(nl, nu), true);
  }

  double lb () const { return mods.zero() ? 0.0 : l.value(GMP_RNDD); }
  double ub () const { return mods.zero() ? 0.0 : u.value(GMP_RNDU); }

  int sign (bool fail = true) const {
    int su = u.sign();
    if (su == -1) {
      if (mods.mixed())
        ++nMixed;
      else if (mods.zero())
        ++nRascal;
      return -1;
    }
    int sl = l.sign();
    if (sl == 1) {
      if (mods.mixed())
        ++nMixed;
      else if (mods.zero())
        ++nRascal;
      return 1;
    }
    if (!fail || (sl == 0 && su == 0) || mods.zero())
      return 0;
    throw SignException(alg);
  }

  MParameter lbP () const {
    return MParameter (l.copy(), l.copy(), Mods(l, l), alg);
  }
  MParameter ubP () const {
    return MParameter(u.copy(), u.copy(), Mods(u, u), alg);
  }

  MParameter midP () const {
    MValue x = l.mid(u);
    return MParameter(x, x.copy(), Mods(x, x), alg);
  }

  bool subset (const MParameter &b) const {
    return !(l < b.l) && !(b.u < u) && (b.l < l || u < b.u);
  }

  bool intersects (const MParameter &b) const { return !(u < b.l || b.u < l); }

  MParameter interval (const MParameter &b) const {
    return MParameter(l.copy(), b.u.copy(), Mods(l, b.u), alg || b.alg);
  }

  MParameter intersect (const MParameter &b) const {
    if (u < b.l || b.u < l) {
      cerr << "error in MParameter::intersect" << endl;
      exit(0);
    }
    MValue nl = l < b.l ? b.l : l, nu = u < b.u ? u : b.u;
    return MParameter(nl.copy(), nu.copy(), Mods(nl, nu), alg || b.alg);
  }

  double distance (const MParameter &b) const {
    if (u < b.l || u == b.l) {
      MValue x = b.l.minus(u, GMP_RNDN);
      double y = x.value();
      x.clear();
      return y;
    }
    if (b.u < l || b.u == l) {
      MValue x = l.minus(b.u, GMP_RNDN);
      double y = x.value();
      x.clear();
      return y;
    }
    return l == b.l && u == b.u ? 0.0 : -1.0;
  }

  double intervalWidth () const {
    MValue d = u.minus(l, GMP_RNDN);
    double w = d.value();
    d.clear();
    return w;
  }
  // polym support
  MParameter value () const { return *this; }
  MParameter make (const MParameter &x) const { return x; }
  PolyM<MParameter> PolyMNtoPolyMRes (const PolyM<MParameter> &p) const;
};

void mparReport ();

void resetMpar ();

