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

double nextDR (double x);

double prevDR (double x);

void rationalize (double x, double& u, double& v);

class PParameter {
  friend class MParameter;

  Parameter p;
  MP_RAT r;

  PParameter (double a, int dummy) {
    MP_INT i;
    mpz_init(&i);
    mpz_set_d(&i, a);
    mpq_init(&r);
    mpq_set_z(&r, &i);
    mpz_clear(&i);
  }

  void init (double a) {
    double u, v;
    rationalize(a, u, v);
    PParameter pu(u, 0), pv(v, 0);
    mpq_init(&r);
    mpq_div(&r, &pu.r, &pv.r);
  }

public:
  static unsigned long int nModSign;
  static unsigned int nMixed, nRascal;

  bool hasCurrentPrimes () { return true; }
  PParameter () { mpq_init(&r); }
  
  PParameter (const PParameter &a) {
    mpq_init(&r);
    if (a.alg())
      p = a.p;
    else
      mpq_set(&r, &a.r);
  }

  PParameter (double a) : PParameter(Parameter(a), false) {}

  PParameter (const Parameter &a, bool alg = false) {
    if (!alg) {
      if (a.lb() != a.ub()) {
	cerr << "internal error in PParameter" << endl;
	exit(0);
      }
      init(a.lb());
    }
    else {
      p = a;
      mpq_init(&r);
    }
  }
  
  ~PParameter () { mpq_clear(&r); }

  unsigned long int numBits () const {
    MP_INT n;
    mpz_init(&n);
    mpq_get_num(&n, &r);
    mpz_abs(&n, &n);
    unsigned long int l = 0;
    while (mpz_sgn(&n) == 1) {
      mpz_tdiv_q_2exp(&n, &n, 1);
      ++l;
    }
    mpz_clear(&n);
    return l;
  }

  void print () const {
    cout << &r << endl;
  }

  PParameter & operator= (const PParameter &x) {
    if (&x == this)
      return *this;
    mpq_clear(&r);
    mpq_init(&r);
    p = x.p;
    if (!x.alg())
      mpq_set(&r, &x.r);
    return *this;
  }

  bool alg () const { return !p.uninitialized(); }
  Parameter par () const { return alg() ? p : Parameter(lb(), ub()); }

  PParameter operator- () const {
    PParameter y;
    if (alg())
      y.p = -p;
    else
      mpq_neg(&y.r, &r);
    return y;
  }
  
  bool operator== (const PParameter &x) const {
    if (alg() || x.alg()) {
      cerr << "internal error in PParameter" << endl;
      exit(0);
    }
    return mpq_cmp(&r, &x.r) == 0;
  }
  
  PParameter operator+ (const PParameter &x) const {
    PParameter y;
    if (alg() || x.alg())
      y.p = par() + x.par();
    else
      mpq_add(&y.r, &r, &x.r);
    return y;
  }
  
  PParameter operator+ (double x) const { return *this + PParameter(x); }
  
  PParameter operator- (const PParameter &x) const {
    PParameter y;
    if (alg() || x.alg())
      y.p = par() - x.par();
    else
      mpq_sub(&y.r, &r, &x.r);
    return y;
  }
  
  PParameter operator- (double x) const { return *this - PParameter(x); }
  
  PParameter operator* (const PParameter &x) const {
    PParameter y;
    if (alg() || x.alg())
      y.p = par()*x.par();
    else
      mpq_mul(&y.r, &r, &x.r);
    return y;
  }
  
  PParameter operator* (double x) const { return *this*PParameter(x); }
  PParameter rcp () const { return PParameter(1) / *this; }

  PParameter operator/ (const PParameter &x) const {
    if (x.sign() == 0) {
      cerr << "reciprocal of zero in PParameter" << endl;
      exit(0);
    }
    PParameter y;
    if (alg() || x.alg())
      y.p = par()/x.par();
    else
      mpq_div(&y.r, &r, &x.r);
    return y;
  }
  
  PParameter operator/ (double x) const { return *this / PParameter(x); }
  PParameter abs () const { return sign() == 1 ? *this : - *this; }
  PParameter sqrt () const { throw SignException(true); }
  
  int sign (bool fail = true) const {
    if (alg())
      try {
        return p.sign(fail);
      }
      catch (SignException e) {
        throw SignException(true);
      }
    return mpq_sgn(&r);
  }
  
  double lb () const { return alg() ? p.lb() : prevDR(mpq_get_d(&r)); }
  double ub () const { return alg() ? p.ub() : nextDR(mpq_get_d(&r)); }
  PParameter (const Parameter &il, const Parameter &iu, bool alg) {
    cerr << "undefined in PParameter" << endl; exit(0);
  }
  PParameter (const MParameter &a);
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
  // polym support
  PParameter value () const { return *this; }
  PParameter make (const PParameter &x) const { return x; }
  PolyM<PParameter> PolyMNtoPolyMRes (const PolyM<PParameter> &p) const;
};

void pparReport ();

void resetPpar ();

