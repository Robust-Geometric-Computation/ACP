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

#define ulong(a) ((unsigned long)(a))

unsigned int inverse (unsigned int a, unsigned int n);

class Mod {
public:
  Mod (unsigned int a = 0, unsigned int p = 0) : a(a), p(p) {}

  Mod operator- () const {
    long l = (-long(a))%p;
    return Mod(l < 0 ? l + p : l, p);
  }

  Mod operator+ (const Mod &x) const {
    return Mod((long(a) + long(x.a))%p, p);
  }

  Mod operator- (const Mod &x) const {
    long l = (long(a) - long(x.a))%p;
    return Mod(l < 0 ? l + p : l, p);
  }

  Mod operator* (const Mod &x) const {
    return Mod((ulong(a)*ulong(x.a))%p, p);
  }

  Mod operator/ (const Mod &x) const {
    if (x.a == 0)
      throw p;
    return Mod((ulong(a)*ulong(inverse(x.a, p)))%p, p);
  }

  unsigned int a, p;
};

class Modder {
  static const int eShift;
  static const int eMax;
  static const int eMin;

  vector<unsigned int> pow2v;
  unsigned int *pow2;
  unsigned int p;

public:
  Modder () : pow2v(0), pow2(0), p(0) {}
  Modder (unsigned int p) : p(p), pow2v(eMax - eMin + 1) {
    pow2 = &pow2v[0] - eMin;
    pow2[0] = 1;
    unsigned long x = 1;
    for (int e = 1; e <= eMax; e++) {
      x = (2 * x)%p;
      pow2[e] = x;
    }
    x = 1;
    unsigned int i2 = (p + 1)/2;
    for (int e = -1; e >= eMin; e--) {
      x = (i2*x)%p;
      pow2[e] = x;
    }
  }

  Mod mod (double x) {
    if (x == (long) x)
      return Mod((x >= 0 ? ((long) x)%p : p + ((long) x)%p), p);
    int e;
    double m = frexp(x, &e)*(1l << eShift);
    if (m != (long) m) {
      cerr << "internal error in Modder::mod" << endl;
      exit(0);
    }
    e -= eShift;
    long mp = ((long) m)%p;
    if (mp < 0) mp += p;
    if (!(eMin <= e && e <= eMax)) {
      cerr << "internal error in Modder::mod" << endl;
      exit(0);
    }
    return Mod(((unsigned long)mp*pow2[e])%p, p);
  }
};

#define NMods 2

unsigned int random32bitPrime ();

class MValue;

class Mods {
public:
  static unsigned int ps[NMods];
  static Modder modder[NMods];
  static pthread_mutex_t mutex;

  static void changePrime (unsigned int p);

  bool hasCurrentPrimes () {
    for (int i = 0; i < NMods; ++i)
      if (mod[i].p != ps[i]) {
        return false;
      }
    return true;
  }

  Mod mod[NMods];

  Mods () {}
  Mods (int start, int it, int up) {
    for (int i = 0; i < NMods; ++i) {
      if (ps[i] != 0) {
	cerr << "nonzero ps[i] in Mods::Mods" << endl;
	exit(0);
      }
      ps[i] = random32bitPrime();
      modder[i] = Modder(ps[i]);
    }
  }

  Mods (double r) {
    for (int i = 0; i < NMods; ++i) mod[i] = modder[i].mod(r);
  }

  Mods (const Parameter &p) {
    double l = p.lb(), u = p.ub();
    if (l == u)
      for (int i = 0; i < NMods; ++i)
	mod[i] = modder[i].mod(l);
    else
      for (int i = 0; i < NMods; ++i)
	mod[i] = Mod(random()%ps[i], ps[i]);
  }

  Mods (const MValue &l, const MValue &u);

  bool hasTheRightPrimes () const {
    for (int i = 0; i < NMods; ++i)
      if (mod[i].p != ps[i])
	return false;
    return true;
  }

  bool mixed () const {
    for (int i = 1; i < NMods; ++i)
      if ((mod[i].a == 0) != (mod[0].a == 0))
	return true;
    return false;
  }

  bool zero () const {
    for (int i = 0; i < NMods; ++i)
      if (mod[i].a != 0)
	return false;
    return true;
  }

  Mods operator- () const {
    Mods m;
    for (int i = 0; i < NMods; ++i)
      m.mod[i] = -mod[i];
    return m;
  }

  Mods operator+ (const Mods &x) const {
    Mods m;
    for (int i = 0; i < NMods; ++i)
      m.mod[i] = mod[i] + x.mod[i];
    return m;
  }

  Mods operator+ (double x) const { return *this + Mods(x); }

  Mods operator- (const Mods &x) const {
    Mods m;
    for (int i = 0; i < NMods; ++i)
      m.mod[i] = mod[i] - x.mod[i];
    return m;
  }

  Mods operator- (double x) const { return *this - Mods(x); }

  Mods operator* (const Mods &x) const {
    Mods m;
    for (int i = 0; i < NMods; ++i)
      m.mod[i] = mod[i]*x.mod[i];
    return m;
  }

  Mods operator* (double x) const { return *this*Mods(x); }

  Mods operator/ (const Mods &x) const {
    Mods m;
    for (int i = 0; i < NMods; ++i)
      m.mod[i] = mod[i]/x.mod[i];
    return m;
  }

  Mods operator/ (double x) const { return *this/Mods(x); }
};

