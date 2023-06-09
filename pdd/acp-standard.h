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

#ifndef ACP
#define ACP

#include "assert.h"
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <exception>
#include <iomanip>
#include <iostream>
#include <mpfr.h>
#include <limits.h>
#include <sys/time.h>

using namespace std;

namespace acp {

#define TMAX 128
extern pthread_key_t idkey;
extern unsigned int threadIds[128], curPrecisions[128];

unsigned int threadId();
unsigned int & curPrecision();

bool isEnabled ();
void enable ();
void disable ();
double randomNumber (double rmin, double rmax);
double getTime ();
double diffTime (const timeval &t2, const timeval &t1);

class SignException : public std::exception {
public:
  bool alg;

  SignException (bool alg = false) : alg(alg) {}
  virtual const char* what() const throw() { return "Not enough precision"; }
};

template<class N>
class PolyM;

#ifndef NOMODE
#include "par.h"
#else
#include "parnm.h"
#endif
#include "mod.h"
#ifndef EXACTRAT
#include "ppar.h"
#else
#include "ratpar.h"
#endif
#ifdef SEPBOUND
#include "mparsep.h"
#else
#include "mpar.h"
#endif

inline Parameter operator+ (double a, Parameter b) {
  return Parameter(a) + b;
}

inline Parameter operator- (double a, Parameter b) {
  return Parameter(a) - b;
}

inline Parameter operator* (double a, Parameter b) {
  return Parameter(a)*b;
}

inline Parameter operator/ (double a, Parameter b) {
  return Parameter(a)/b;
}

inline PParameter operator+ (double a, PParameter b) {
  return PParameter(a) + b;
}

inline PParameter operator- (double a, PParameter b) {
  return PParameter(a) - b;
}

inline PParameter operator* (double a, PParameter b) {
  return PParameter(a)*b;
}

inline PParameter operator/ (double a, PParameter b) {
  return PParameter(a)/b;
}

inline MParameter operator+ (double a, MParameter b) {
  return MParameter(a) + b;
}

inline MParameter operator- (double a, MParameter b) {
  return MParameter(a) - b;
}

inline MParameter operator* (double a, MParameter b) {
  return MParameter(a)*b;
}

inline MParameter operator/ (double a, MParameter b) {
  return MParameter(a)/b;
}

inline double psqrt (double a) { return sqrt(a); }
inline Parameter psqrt (Parameter a) { return a.sqrt(); }
inline PParameter psqrt (PParameter a) { return a.sqrt(); }
inline MParameter psqrt (MParameter a) { return a.sqrt(); }
inline bool isZero (double a) { return a == 0.0; }
inline bool isZero (const Parameter &a) { return a.sign() == 0; }
inline bool isZero (const PParameter &a) { return a.sign() == 0; }
inline bool isZero (const MParameter &a) { return a.sign() == 0; }

#include "pv.h"

class RefCnt {
  template<class T> friend class PTR;
  int refCnt;
  
  void incRef () {
    pthread_mutex_lock(&mutex);
    refCnt++;
    pthread_mutex_unlock(&mutex);
  }

  void decRef () { 
    pthread_mutex_lock(&mutex);
    bool flag = --refCnt == 0;
    pthread_mutex_unlock(&mutex);
    if (flag)
      delete this;
  }
protected:
  pthread_mutex_t mutex;
public:
  RefCnt () : refCnt(0) { mutex = PTHREAD_MUTEX_INITIALIZER; }
  virtual ~RefCnt () {
    if (refCnt > 0) {
      cerr << "invalid smart pointer" << endl;
      exit(0);
    }
  }
};

template<class T>
class PTR {
  T *t;
public:
  PTR () : t(0) {}
  PTR (T *t) : t(t) { incRef(); }
  PTR (const PTR &p) : t(p.t) { incRef(); }
  const PTR &operator= (const PTR &p) { 
    p.incRef(); decRef(); t = p.t; return *this; 
  }
  ~PTR () { decRef(); }
  void incRef () const { if (t != 0) t->incRef(); }
  void decRef () const { if (t != 0) t->decRef(); }
  operator T* () const { return t; }
  T *operator-> () const { return t; }
};

class BaseObject : public RefCnt {
protected:
  static unsigned int highPrecision, maxPrecision, homPrecision;
  static int perIs[TMAX];
  static double perTs[TMAX], maxP;
  static int & perI () { return perIs[threadId()]; }
  static double & perT () { return perTs[threadId()]; }
public:
  virtual ~BaseObject () {}
  static double delta;
  static void addThread (unsigned int i) {
    if (i >= TMAX) {
      cerr << "too many threads" << endl;
      exit(0);
    }
    pthread_setspecific(idkey, (void *) (threadIds + i));
    threadIds[i] = i;
    curPrecisions[i] = 53u;
    perIs[i] = -1;
  }
};

class PrecisionException : public std::exception {
public:
  virtual const char* what() const throw() {
    return "Maximum precision exceeded";
  }
};

class Interval {
public:
  double l, u;

  Interval () : l(0.0), u(0.0) {}
  Interval (double l, double u) : l(l), u(u) {}
  Interval (const Parameter &a) : l(a.lb()), u(a.ub()) {}
  double lb () const { return l; }
  double ub () const { return u; }
  double mid () const { return 0.5*(l + u); }
  bool intersects (const Interval &b) {
    return !(ub() < b.lb() || b.ub() < lb());
  }
  int sign () const {
    if (l > 0.0) return 1;
    if (u < 0.0) return -1;
    return 0;
  }
};

#define NHoms 1

template<template<typename> class P>
class Object : public BaseObject {
  friend class Parameter;
  class Elt {
  public:
    double t, *v;
    P<Parameter> p;
    P<MParameter> *mp;

    Elt (double t = 0.0) : t(t), v(0), mp(0) {}
    
    void setV (unsigned int n) {
      v = new double [n];
      for (int i = 0; i < n; ++i)
	v[i] = randomNumber(-1.0, 1.0);      
    }
  };
  
  typedef vector<Elt> Elts;
  
  P<Parameter> p;
  unsigned int pPrecision;
  P<PParameter> *pp;
  P<MParameter> *mp;

  bool checkAccuracy (double acc) {
    if (acc == 1.0) return true;
    for (int i = 0; i < p.size(); i++) {
      double l = p[i].lb(), u = p[i].ub();  
      if (!(acc > 1e-17 && u - l < acc*max(1.0,min(fabs(l), fabs(u)))) &&
	  u - l > acc*min(fabs(l), fabs(u))) {
	double ln = nextafter(l, 1.0 + u);
	if (ln < u) {
	  ln = nextafter(ln, 1.0 + u);
	  if (ln < u)
	    return false;
	}
      }
    }
    return true;
  }

 protected:
  Elts *e;
  bool input, clamp;

  virtual P<Parameter> calculateP () {
    cerr << "missing calculate method" << endl;
    exit(0);
    return P<Parameter>();
  }

  virtual P<PParameter> calculatePP () { 
    cerr << "missing calculate method" << endl; 
    exit(0); 
    return P<PParameter>();
  }

  virtual P<MParameter> calculateMP () { 
    cerr << "missing calculate method" << endl; 
    exit(0); 
    return P<MParameter>();
  }

  void noClear (P<MParameter> &x) const {
    for (int i = 0; i < x.size(); ++i) {
      if (!x[i].clear) {
	x[i].l = x[i].l.copy();
	x[i].u = x[i].u.copy();
      }
      x[i].clear = false;
    }
  }

  void clear (P<MParameter> *x) const {
    if (x) {
      for (int i = 0; i < x->size(); ++i)
	(*x)[i].clear = true;
      delete x;
    }
  }
  
 public:
  Object () : pPrecision(0), input(false), clamp(false), pp(0), mp(0), e(0) {}

  Object (const P<double> &d, bool perturb, bool clamp = false)
    : p(d), pPrecision(53u), input(true), clamp(clamp), pp(0), mp(0), e(0) {
    if (perturb)
      for (int i = 0; i < p.size(); ++i)
	p[i] = Parameter(d[i] + delta*(1.0 + fabs(d[i]))*randomNumber(-1.0, 1.0));
  }

  ~Object () {
    delete pp;
    clear(mp);
    if (e) {
      for (unsigned int i = 0u; i < NHoms; ++i)
	for (unsigned int j = 0u; j < e[i].size(); ++j) {
	  if (j == 0)
	    delete [] e[i][j].v;
	  clear(e[i][j].mp);
	}
      delete [] e; 
    }
  }

  bool uninitialized () const { return pPrecision == 0; }
  P<Parameter> getCurrentP () const { return p; }

  template<class N>
    const P<N> & get () {
    //if (input && std::is_same<N, Parameter>::value) // eps optimization
      //return get1(N());
    try {
      pthread_mutex_lock(&mutex);
      const P<N> &ne = get1(N());
      pthread_mutex_unlock(&mutex);
      return ne;
    }
    catch (SignException se) {
      pthread_mutex_unlock(&mutex);
      throw se;
    }
    catch (unsigned int p) {
      pthread_mutex_unlock(&mutex);
      throw p;
    }
  }

  const P<Parameter> & get1 (const Parameter &dummy) {
    if (uninitialized()) {
      p = calculateP();
      pPrecision = 53u;
    }
    return p;
  }

  const P<PParameter> & get1 (const PParameter &dummy) {
    if (pp) {
      if ((*pp)[0].hasCurrentPrimes()) {
	return *pp;
      }
      delete pp;
      pp = 0;
      if (mp == 0)
	pPrecision = 53u;
    }
#ifndef EXACTRAT
    if (mp) {
      if ((*mp)[0].hasCurrentPrimes()) {
	pp = new P<PParameter>(*mp);
	return *pp;
      }
      delete mp;
      mp = 0;
      pPrecision = 53u;
    }
#endif
    pp = new P<PParameter>(input ? P<PParameter>(p) : calculatePP());
    if (uninitialized() || !input)
      p = *pp;
    pPrecision = 100u;
    return *pp;
  }

  const P<MParameter> & get1 (const MParameter &dummy) {
    return perI() == -1 || clamp ? getM() : getMP();
  }
  
  const P<MParameter> & getM () {
    if (pPrecision < curPrecision() ||
        (mp && !(*mp)[0].hasCurrentPrimes())) { // added
      P<MParameter> *newMp = input ? new P<MParameter>(p) :
	new P<MParameter>(calculateMP());
      noClear(*newMp);
      clear(mp);
      mp = newMp;
      if (!input)
	p = *newMp;
      pPrecision = curPrecision();
    }
    return *mp;
  }

  const P<MParameter> & getMP () {
    if (!e)
      e = new Elts [NHoms];
    Elts &elts = e[perI()];
    double t = perT();
    for (unsigned int i = 0u; i < elts.size(); ++i)
      if (elts[i].t == t)
	return *elts[i].mp;
    Elt elt(t);
    if (input) {
      double *v;
      if (elts.empty()) {
	elt.setV(p.size());   
	v = elt.v;
      }
      else
	v = elts[0].v;
      elt.mp = new P<MParameter>(p);
      for (unsigned int i = 0u; i < p.size(); ++i)
	(*elt.mp)[i] = MParameter(p[i].lb()) + MParameter(t)*MParameter(v[i]);
    }
    else
      elt.mp = new P<MParameter>(calculateMP());
    noClear(*elt.mp);
    elts.push_back(elt);
    return *elt.mp;
  }

  P<Interval> getApprox (double acc = 1e-17) {
    int peri = perI(), cp = curPrecision();
    perI() = -1;
    curPrecision() = 53u;
    try {
      get<Parameter>();
      if (checkAccuracy(acc)) {
	perI() = peri;
	curPrecision() = cp;
	return p;
      }
    }
    catch (SignException se) {}
    catch (unsigned int p) { Mods::changePrime(p); }
    curPrecision() = highPrecision;
    while (true) {
      try {
	get<MParameter>();
        if (checkAccuracy(acc)) {
	  perI() = peri;
	  curPrecision() = cp;
	  return p;
        }
      }
      catch (SignException se) {}
      catch (unsigned int p) { Mods::changePrime(p); }
      curPrecision() *= 2u;
    }
  }

private:
  class IntervalOD : public Interval {
  public:
    IntervalOD () {}
    IntervalOD (const Interval &interval) : Interval(interval) {}
    operator double() const { return 0.5*(l + u); }
  };
  class IntervalOP : public Interval {
  public:
    IntervalOP () {}
    IntervalOP (const Interval &interval) : Interval(interval) {}
    operator Parameter() const { return Parameter(lb(), ub()); }
  };
public:
  P<double> getApproxMid (double acc = 1e-17) { return P<IntervalOD>(getApprox(acc)); }
  P<Parameter> getApproxParameter (double acc = 1e-17) { return P<IntervalOP>(getApprox(acc)); }
};

class Primitive : public BaseObject {
  virtual Parameter calculateP () = 0;
  virtual PParameter calculatePP () = 0;
  virtual MParameter calculateMP () = 0;
 public:
  static unsigned long int nPrim, nAmb, nDeg, nAlg, nHom, nId, maxBits;
  static double minPrim, minPR, tPar, tPPar, tMPar; 
  operator int () {
    ++nPrim;
    int cp = curPrecision();
    timeval t0, t1, t2, t3;
    gettimeofday(&t0, 0);
    bool alg = false;
    try {
      int s = calculateP().sign();
      gettimeofday(&t1, 0);
      tPar += diffTime(t1, t0);
      return s;
    }
    catch (SignException se) {}
    gettimeofday(&t1, 0);
    tPar += diffTime(t1, t0);
    ++nAmb;
    while (true) {
      try {
	PParameter x = calculatePP();
	int s = x.sign();
	if (s == 0) ++nDeg;
#ifdef EXACTRAT
	else
	  maxBits = max(maxBits, x.numBits());
#endif
	gettimeofday(&t2, 0);
	tPPar += diffTime(t2, t1);
	return s;
      }
      catch (unsigned int p) { Mods::changePrime(p); }
      catch (SignException se) {
	gettimeofday(&t2, 0);
	tPPar += diffTime(t2, t1);
	if (se.alg) {
	  ++nAlg;
	  alg = true;
	  int s = 2;
	  bool hom = homotopy(s);
	  if (s == 0) ++nDeg;
	  if (hom || s != 2) {
	    gettimeofday(&t3, 0);
	    tMPar += diffTime(t3, t2);
	    curPrecision() = cp;
	    return hom ? 0 : s;
	  }
	}
	break;
      }
    }
    curPrecision() = alg ? homPrecision : highPrecision;
    static unsigned int maxPrecisionUsed = highPrecision;
    while (!alg || curPrecision() <= maxPrecision) {
      try {
	MParameter x = calculateMP();
	int s = x.sign();
	if (s == 0)
	  ++nDeg;
	curPrecision() = cp;
	gettimeofday(&t3, 0);
	tMPar += diffTime(t3, t2);
	return s;
      }
      catch (SignException se) {
	curPrecision() *= 2u;
	if (curPrecision() > maxPrecisionUsed) {
	  maxPrecisionUsed = curPrecision();
	  cerr << "max precision used = " << maxPrecisionUsed << endl;
	}
      }
    }
#ifdef NO_HOMOTOPY
    cerr << "max precision exceeded!" << endl;
    curPrecision() = cp;
    return 0;
#endif
    throw PrecisionException();
  }

  bool homotopy (int &s) {
#ifdef NO_HOMOTOPY
    return false;
#endif
    s = 2;
    bool mixed = false;
    curPrecision() = homPrecision;
    MParameter x;
    while (true) {
      try {
	x = calculateMP();
	int sx = x.sign(false);
	if (sx || (x.lb() == 0.0 && x.ub() == 0.0))
	  s = sx;
	if (sx)
	  minPrim = min(minPrim,(s == 1 ? x.lb() : - x.ub()));
	break;
      }
      catch (SignException se) {
	curPrecision() *= 2u;
	if (curPrecision() > maxPrecision)
	  throw PrecisionException();
      }
    }
    if (s == 2)
      ++nHom;
    int n = 0;
    for (int i = 0; i < NHoms; ++i) {
      perI() = i;
      MParameter xp = homotopyAux();
      if (xp.sign(false) == 0)
	++n;
      if (s == 1 || s == -1) {
	double d = xp.distance(x);
	if (d > 0.0)
	  minPR = min(minPR, d/(xp - x).intervalWidth());
      }
    }
    perI() = -1;
    if ((n > 0 && n < NHoms) || ((s == 1 || s == -1) && n > 0))
      cerr << "mixed homotopy " << n << endl;
    if (n == 0)
      return false;
    if (s == 2)
      ++nId;
    return true;
  }

  MParameter homotopyAux () {
    double pert = perT();
    perT() = maxP;
    while (true) {
      try {
	MParameter x = calculateMP();
	perT() = pert;
	return x;
      }
      catch (SignException se) { perT() *= 0.1; }
    }
  }
};

#define DeclareCalculate(P)						\
  private:								\
    P<Parameter> calculateP () { return calculate<Parameter>(); }	\
    P<PParameter> calculatePP () { return calculate<PParameter>(); }	\
    P<MParameter> calculateMP () { return calculate<MParameter>(); }	\
    template<class N>							\
    P<N> calculate ()

#define DeclareSign							\
  private:								\
    Parameter calculateP () { return calculate<Parameter>(); }		\
    PParameter calculatePP () { return calculate<PParameter>(); }	\
    MParameter calculateMP () { return calculate<MParameter>(); }	\
    template<class N>							\
    N calculate ()

void primitiveReport ();

void resetReport ();
}
#endif
