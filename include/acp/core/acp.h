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

//#define NO_MODE

#include <assert.h>
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <algorithm>
#include <exception>
#include <iomanip>
#include <iostream>
#include <vector>
#include <mpfr.h>

using namespace std;

namespace acp {

#define TMAX 128
extern pthread_key_t idkey;
extern unsigned int threadIds[128], curPrecisions[128];

unsigned int threadId();
unsigned int& curPrecision();

extern bool penabled, enabled;
void enable();
void disable();
bool isEnabled();

double randomNumber(double rmin, double rmax);

double nextD(double x);

double prevD(double x);

class SignException : public std::exception {
 public:
  bool alg;

  SignException(bool alg = false) : alg(alg) {}
  virtual const char* what() const throw() { return "Not enough precision"; }
};

class PParameter;
class MParameter;

class Parameter {
  friend class PParameter;
  friend class MParameter;
  double l, u;

#ifndef NO_MODE
  Parameter(double l, double u) : l(l), u(-u) {}
  Parameter(double l, double u, int i) : l(l), u(u) {}
#else
  Parameter(double l, double u) : l(l), u(u) {}
#endif

  Parameter(double x) : Parameter(x, x) {}
  Parameter sqrt(double x) const;

 public:
  static double delta;  // default is 2^{-26}
  void setSize(unsigned int s) {}
  int size() const { return 1; }
  Parameter& operator[](int i) { return *this; }
  const Parameter& elt(int i) const { return *this; }
  double lb() const { return l; }
#ifndef NO_MODE
  double ub() const { return -u; }
#else
  double ub() const { return u; }
#endif
  Parameter() : Parameter(0.0, -1.0) {}
  bool uninitialized() const { return l == 0.0 && ub() == -1.0; }
  Parameter(int x) : Parameter(x, x) {}

  static Parameter input(double x) {
    double y = x + delta * (1.0 + fabs(x)) * randomNumber(-1.0, 1.0);
    return Parameter(y, y);
  }

  static Parameter constant(double x) { return Parameter(x, x); }
  static Parameter interval(double x, double y) { return Parameter(x, y); }
  Parameter(const Parameter& il, const Parameter& iu, bool alg)
      : l(il.l), u(iu.u) {}
  void setAlg() {}
  Parameter(const PParameter& e);
  Parameter(const MParameter& e);

  int sign(bool fail = true) const {
    if (l > 0.0) return 1;
    if (ub() < 0.0) return -1;
    if (!fail || (l == 0.0 && u == 0.0)) return 0;
    throw SignException();
  }

  Parameter operator+(const Parameter& b) const {
#ifndef NO_MODE
    return Parameter(l + b.l, u + b.u, 0);
#else
    return Parameter(prevD(l + b.l), nextD(u + b.u));
#endif
  }

  Parameter operator+(double b) const { return *this + constant(b); }

  Parameter operator-() const {
#ifndef NO_MODE
    return Parameter(u, l, 0);
#else
    return Parameter(-u, -l);
#endif
  }

  Parameter operator-(const Parameter& b) const {
#ifndef NO_MODE
    return Parameter(l + b.u, u + b.l, 0);
#else
    return Parameter(prevD(l - b.u), nextD(u - b.l));
#endif
  }

  Parameter operator-(double b) const { return *this - constant(b); }

  Parameter operator*(const Parameter& b) const {
#ifndef NO_MODE
    Parameter s = u >= 0.0 ? -*this : *this, t = u >= 0.0 ? -b : b;
    if (s.l > 0.0) {
      if (t.l >= 0.0)
        return Parameter(s.l * t.l, -s.u * t.u, 0);
      else if (t.u >= 0.0)
        return Parameter(-s.u * t.l, s.l * t.u, 0);
      else
        return Parameter(-s.u * t.l, -s.u * t.u, 0);
    }
    return Parameter(min(-s.l * t.u, -s.u * t.l), min(-s.l * t.l, -s.u * t.u),
                     0);
#else
    Parameter s = u < 0.0 ? -*this : *this, t = u < 0.0 ? -b : b;
    if (s.l >= 0.0) {
      if (t.l >= 0.0)
        return Parameter(prevD(s.l * t.l), nextD(s.u * t.u));
      else if (t.u <= 0.0)
        return Parameter(prevD(s.u * t.l), nextD(s.l * t.u));
      else
        return Parameter(prevD(s.u * t.l), nextD(s.u * t.u));
    }
    if (t.l >= 0.0) return Parameter(prevD(s.l * t.u), nextD(s.u * t.u));
    if (t.u <= 0.0) return Parameter(prevD(s.u * t.l), nextD(s.l * t.l));
    double k1 = s.l * t.u, k2 = s.u * t.l, nl = k1 < k2 ? k1 : k2,
           k3 = s.l * t.l, k4 = s.u * t.u, nu = k3 < k4 ? k4 : k3;
    return Parameter(prevD(nl), nextD(nu));
#endif
  }

  Parameter operator*(double b) const { return *this * constant(b); }

  Parameter rcp() const {
    assert(sign());
#ifndef NO_MODE
    return Parameter(-1.0 / u, -1.0 / l, 0);
#else
    return Parameter(prevD(1.0 / u), nextD(1.0 / l));
#endif
  }

  Parameter operator/(const Parameter& b) const { return *this * b.rcp(); }
  Parameter operator/(double b) const { return *this * constant(b).rcp(); }
  Parameter abs() const { return sign() == 1 ? *this : -*this; }
  Parameter sqrt() const;
  double mid() const { return 0.5 * (l + ub()); }
  bool operator<(const Parameter& b) const { return (b - *this).sign() == 1; }
  bool operator>(const Parameter& b) const { return (b - *this).sign() == -1; }
  Parameter lbP() const { return Parameter(l, l); }
  Parameter ubP() const { return Parameter(ub(), ub()); }
  Parameter midP() const {
    double m = mid();
    return Parameter(m, m);
  }

  bool subset(const Parameter& b) const {
    return !(lb() < b.lb()) && !(b.ub() < ub()) &&
           (b.lb() < lb() || ub() < b.ub());
  }

  bool intersects(const Parameter& b) const {
    return !(ub() < b.lb() || b.ub() < lb());
  }

  Parameter interval(const Parameter& b) const { return Parameter(l, b.ub()); }

  Parameter intersect(const Parameter& b) const {
    assert(!(ub() < b.lb() || b.ub() < lb()));
    double il = lb() < b.lb() ? b.lb() : lb();
    double iu = ub() < b.ub() ? ub() : b.ub();
    return Parameter(il, iu);
  }

  double distance(const Parameter& b) const {
    if (ub() <= b.lb()) return b.lb() - ub();
    if (b.ub() <= lb()) return lb() - b.ub();
    return lb() == b.lb() && ub() == b.ub() ? 0.0 : -1.0;
  }

  double intervalWidth() const { return ub() - lb(); }
};

unsigned int inverse(unsigned int a, unsigned int n);

#define ulong(a) ((unsigned long)(a))

class Mod {
 public:
  Mod(unsigned int a = 0, unsigned int p = 0) : a(a), p(p) {}

  Mod operator-() const {
    long l = (-long(a)) % p;
    return Mod(l < 0 ? l + p : l, p);
  }

  Mod operator+(const Mod& x) const {
    return Mod((long(a) + long(x.a)) % p, p);
  }

  Mod operator-(const Mod& x) const {
    long l = (long(a) - long(x.a)) % p;
    return Mod(l < 0 ? l + p : l, p);
  }

  Mod operator*(const Mod& x) const {
    return Mod((ulong(a) * ulong(x.a)) % p, p);
  }

  Mod operator/(const Mod& x) const {
    if (x.a == 0) throw p;
    return Mod((ulong(a) * ulong(inverse(x.a, p))) % p, p);
  }

  unsigned int a, p;
};

class Modder {
  static const int eShift;
  static const int eMax;
  static const int eMin;

  vector<unsigned int> pow2v;
  unsigned int* pow2;
  unsigned int p;

 public:
  Modder() : pow2v(0), pow2(0), p(0) {}
  Modder(unsigned int p) : p(p), pow2v(eMax - eMin + 1) {
    pow2 = &pow2v[0] - eMin;
    pow2[0] = 1;
    unsigned long x = 1;
    for (int e = 1; e <= eMax; e++) {
      x = (2 * x) % p;
      pow2[e] = x;
    }
    x = 1;
    unsigned int i2 = (p + 1) / 2;
    for (int e = -1; e >= eMin; e--) {
      x = (i2 * x) % p;
      pow2[e] = x;
    }
  }

  Mod mod(double x) {
    if (x == (long)x)
      return Mod((x >= 0 ? ((long)x) % p : p + ((long)x) % p), p);
    int e;
    double m = frexp(x, &e) * (1l << eShift);
    assert(m == (long)m);
    e -= eShift;
    long mp = ((long)m) % p;
    if (mp < 0) mp += p;
    assert(eMin <= e && e <= eMax);
    return Mod(((unsigned long)mp * pow2[e]) % p, p);
  }
};

#define NMods 3
//#define DEGREE_BOUND
unsigned int random32bitPrime();

class MValue;

class Mods {
 public:
  static unsigned int ps[NMods];
  static Modder modder[NMods];
  static pthread_mutex_t mutex;

  static void changePrime(unsigned int p);

  bool hasCurrentPrimes() {
    for (int i = 0; i < NMods; i++)
      if (mod[i].p != ps[i]) {
        return false;
      }
    return true;
  }

  Mod mod[NMods];
#ifdef DEGREE_BOUND
  double degNu, degDe;
#endif

  void noBits() {
#ifdef DEGREE_BOUND
    degNu = 0;
    degDe = 0;
#endif
  }

  static double degNuMax;

  Mods()
#ifdef DEGREE_BOUND
      : degNu(0),
        degDe(0)
#endif
  {
  }

  Mods(double degNu, double degDe)
#ifdef DEGREE_BOUND
      : degNu(degNu),
        degDe(degDe)
#endif
  {
#ifdef DEGREE_BOUND
    if (degNu > degNuMax) {
      cerr << "degNuMax " << degNu << endl;
      degNuMax = degNu;
    }
#endif
  }

  Mods(int start, int it, int up)
#ifdef DEGREE_BOUND
      : degNu(1),
        degDe(0)
#endif
  {
    for (int i = 0; i < NMods; i++) {
      assert(ps[i] == 0);
      ps[i] = random32bitPrime();
      modder[i] = Modder(ps[i]);
    }
  }

  Mods(double r)
#ifdef DEGREE_BOUND
      : degNu(1),
        degDe(0)
#endif
  {
    for (int i = 0; i < NMods; i++) mod[i] = modder[i].mod(r);
  }

  Mods(const Parameter& p)
#ifdef DEGREE_BOUND
      : degNu(1),
        degDe(0)
#endif
  {
    double l = p.lb(), u = p.ub();
    if (l == u)
      for (int i = 0; i < NMods; i++) mod[i] = modder[i].mod(l);
    else
      for (int i = 0; i < NMods; i++) mod[i] = Mod(random() % ps[i], ps[i]);
  }

  Mods(const MValue& l, const MValue& u);

  bool hasTheRightPrimes() const {
    for (int i = 0; i < NMods; i++)
      if (mod[i].p != ps[i]) return false;
    return true;
  }

  bool mixed() const {
    for (int i = 1; i < NMods; i++)
      if ((mod[i].a == 0) != (mod[0].a == 0)) return true;
    return false;
  }

  bool zero() const {
    for (int i = 0; i < NMods; i++)
      if (mod[i].a != 0) return false;
    return true;
  }

  Mods operator-() const {
#ifdef DEGREE_BOUND
    Mods m(degNu, degDe);
#else
    Mods m;
#endif
    for (int i = 0; i < NMods; i++) m.mod[i] = -mod[i];
    return m;
  }

  Mods operator+(const Mods& x) const {
#ifdef DEGREE_BOUND
    Mods m(std::max(degNu + x.degDe, degDe + x.degNu), degDe + x.degDe);
#else
    Mods m;
#endif
    for (int i = 0; i < NMods; i++) m.mod[i] = mod[i] + x.mod[i];
    return m;
  }

  Mods operator+(double x) const { return *this + Mods(x); }

  Mods operator-(const Mods& x) const {
#ifdef DEGREE_BOUND
    Mods m(std::max(degNu + x.degDe, degDe + x.degNu), degDe + x.degDe);
#else
    Mods m;
#endif
    for (int i = 0; i < NMods; i++) m.mod[i] = mod[i] - x.mod[i];
    return m;
  }

  Mods operator-(double x) const { return *this - Mods(x); }

  Mods operator*(const Mods& x) const {
#ifdef DEGREE_BOUND
    Mods m(degNu + x.degNu, degDe + x.degDe);
#else
    Mods m;
#endif
    for (int i = 0; i < NMods; i++) m.mod[i] = mod[i] * x.mod[i];
    return m;
  }

  Mods operator*(double x) const { return *this * Mods(x); }

  Mods operator/(const Mods& x) const {
#ifdef DEGREE_BOUND
    Mods m(degNu + x.degDe, degDe + x.degNu);
#else
    Mods m;
#endif
    for (int i = 0; i < NMods; i++) m.mod[i] = mod[i] / x.mod[i];
    return m;
  }

  Mods operator/(double x) const { return *this / Mods(x); }
};

void mixedReport(int id, const Mods& mods);

void rationalize(double x, double& u, double& v);

//#define EXACTRAT

#ifdef EXACTRAT
class PParameter {
  friend class MParameter;

  Parameter p;
  MP_RAT r;

  PParameter(double a, int dummy) {
    MP_INT i;
    mpz_init(&i);
    mpz_set_d(&i, a);
    mpq_init(&r);
    mpq_set_z(&r, &i);
    mpz_clear(&i);
  }

  void init(double a) {
    double u, v;
    rationalize(a, u, v);
    PParameter pu(u, 0), pv(v, 0);
    mpq_init(&r);
    mpq_div(&r, &pu.r, &pv.r);
  }

 public:
  static unsigned long int nModSign;
  static unsigned int nMixed, nRascal;

  bool hasCurrentPrimes() { return true; }

  static PParameter constant(double a) { return PParameter(a); }
  PParameter() { mpq_init(&r); }
  PParameter(const PParameter& a) {
    mpq_init(&r);
    if (a.alg())
      p = a.p;
    else
      mpq_set(&r, &a.r);
  }
  PParameter(const Parameter& a, bool alg = false) {
    if (!alg) {
      assert(a.lb() == a.ub());
      init(a.lb());
    } else {
      p = a;
      mpq_init(&r);
    }
  }
  void setAlg() { alg = true; }

  ~PParameter() { mpq_clear(&r); }

  unsigned long int numBits() const {
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

  PParameter& operator=(const PParameter& x) {
    if (&x == this) return *this;
    mpq_clear(&r);
    mpq_init(&r);
    p = x.p;
    if (!x.alg()) mpq_set(&r, &x.r);
    return *this;
  }

  bool alg() const { return !p.uninitialized(); }
  Parameter par() const { return alg() ? p : Parameter(lb(), ub()); }

  PParameter operator-() const {
    PParameter y;
    if (alg())
      y.p = -p;
    else
      mpq_neg(&y.r, &r);
    return y;
  }
  bool operator==(const PParameter& x) const {
    assert(!alg() && !x.alg());
    return mpq_cmp(&r, &x.r) == 0;
  }
  PParameter operator+(const PParameter& x) const {
    PParameter y;
    if (alg() || x.alg())
      y.p = par() + x.par();
    else
      mpq_add(&y.r, &r, &x.r);
    return y;
  }
  PParameter operator+(double x) const { return *this + constant(x); }
  PParameter operator-(const PParameter& x) const {
    PParameter y;
    if (alg() || x.alg())
      y.p = par() - x.par();
    else
      mpq_sub(&y.r, &r, &x.r);
    return y;
  }
  PParameter operator-(double x) const { return *this - constant(x); }
  PParameter operator*(const PParameter& x) const {
    PParameter y;
    if (alg() || x.alg())
      y.p = par() * x.par();
    else
      mpq_mul(&y.r, &r, &x.r);
    return y;
  }
  PParameter operator*(double x) const { return *this * constant(x); }
  PParameter operator/(const PParameter& x) const {
    assert(x.sign());
    PParameter y;
    if (alg() || x.alg())
      y.p = par() / x.par();
    else
      mpq_div(&y.r, &r, &x.r);
    return y;
  }
  PParameter operator/(double x) const { return *this / constant(x); }
  PParameter sqrt() const { throw SignException(true); }
  int sign(bool fail = true) const {
    if (alg()) try {
        return p.sign();
      } catch (SignException e) {
        throw SignException(true);
      }
    return mpq_sgn(&r);
  }
  double lb() const { return alg() ? p.lb() : prevD(mpq_get_d(&r)); }
  double ub() const { return alg() ? p.ub() : nextD(mpq_get_d(&r)); }
  double mid() const { return alg() ? p.mid() : mpq_get_d(&r); }

  PParameter(const Parameter& il, const Parameter& iu, bool alg) { assert(0); }
  PParameter(const MParameter& a);
  PParameter abs() const { assert(0); }
  bool operator<(const PParameter& b) const { return (*this - b).sign() < 0; }
  bool operator>(const PParameter& b) const { return (*this - b).sign() > 0; }
  PParameter lbP() const { assert(0); }
  PParameter ubP() const { assert(0); }
  PParameter midP() const { assert(0); }
  bool subset(const PParameter& b) const { assert(0); }
  bool intersects(const PParameter& b) const { assert(0); }
  PParameter interval(const PParameter& b) const { assert(0); }
  PParameter intersect(const PParameter& b) const { assert(0); }
};

#else
class PParameter {
  friend class MParameter;
  Parameter p;
  Mods mods;
  bool alg;

 public:
  static unsigned long int nModSign;
  static unsigned int nMixed, nRascal;
  static PParameter constant(double a) { return PParameter(a); }
  PParameter() {}
  PParameter(const Parameter& p, bool alg = false) : p(p), mods(p), alg(alg) {}

  PParameter(const Parameter& p, const Mods& mods, bool alg)
      : p(p), mods(mods), alg(alg) {
    if (sign(false)) ++nModSign;
  }

  void setAlg() { alg = true; }

  PParameter(const MParameter& p);

  bool hasCurrentPrimes() { return mods.hasCurrentPrimes(); }

  PParameter operator+(const PParameter& b) const {
    return PParameter(p + b.p, mods + b.mods, alg || b.alg);
  }

  PParameter operator+(double b) const { return *this + constant(b); }

  PParameter operator-(const PParameter& b) const {
    return PParameter(p - b.p, mods - b.mods, alg || b.alg);
  }

  PParameter operator-(double b) const { return *this - constant(b); }

  PParameter operator-() const { return PParameter(-p, -mods, alg); }

  PParameter operator*(const PParameter& b) const {
    return PParameter(p * b.p, mods * b.mods, alg || b.alg);
  }

  PParameter operator*(double b) const { return *this * constant(b); }

  PParameter operator/(const PParameter& b) const {
    assert(b.sign());
    try {
      return PParameter(p / b.p, mods / b.mods, alg || b.alg);
    } catch (SignException se) {
      throw SignException(alg || b.alg);
    }
  }

  PParameter operator/(double b) const { return *this / constant(b); }

  PParameter sqrt() const { throw SignException(true); }
  double lb() const { return mods.zero() ? 0.0 : p.lb(); }
  double ub() const { return mods.zero() ? 0.0 : p.ub(); }
  double mid() const { return p.mid(); }

  int sign(bool fail = true) const {
    if (p.lb() > 0.0 || p.ub() < 0.0) {
      if (mods.mixed()) {
        ++nMixed;
        mixedReport(0, mods);
        for (int i = 0; i < NMods; i++)
          if (mods.mod[i].a == 0) throw mods.mod[i].p;
      } else if (mods.zero())
        ++nRascal;
      return p.lb() > 0.0 ? 1 : -1;
    }
    if (!fail || (p.lb() == 0 && p.ub() == 0) || mods.zero()) return 0;
    throw SignException(alg);
  }

  PParameter(const Parameter& il, const Parameter& iu, bool alg) { assert(0); }
  PParameter abs() const { return sign() < 0 ? -*this : *this; }
  bool operator<(const PParameter& b) const { return (*this - b).sign() < 0; }
  bool operator>(const PParameter& b) const { return (*this - b).sign() > 0; }
  PParameter lbP() const { assert(0); }
  PParameter ubP() const { assert(0); }
  PParameter midP() const { assert(0); }
  bool subset(const PParameter& b) const { assert(0); }
  bool intersects(const PParameter& b) const { assert(0); }
  PParameter interval(const PParameter& b) const { assert(0); }
  PParameter intersect(const PParameter& b) const { assert(0); }
  double intervalWidth() const { return p.intervalWidth(); }
};
#endif

class MValue {
  friend class Mods;
  mpfr_t m;
  unsigned int p;

  MValue(unsigned int p) : p(p) { mpfr_init2(m, p); }

 public:
  MValue() : p(0) {}

  MValue(double x, unsigned int p) : p(p) {
    mpfr_init2(m, p);
    mpfr_set_d(m, x, GMP_RNDN);
  }

  void clear() { mpfr_clear(m); }

  MValue copy() const {
    MValue v(p);
    mpfr_set(v.m, m, GMP_RNDN);
    return v;
  }

  unsigned int getP() const { return p; }
  double value(mpfr_rnd_t round = GMP_RNDN) const {
    return mpfr_get_d(m, round);
  }

  MValue plus(const MValue& b, mpfr_rnd_t round) const {
    MValue res(max(p, b.p));
    mpfr_add(res.m, m, b.m, round);
    return res;
  }

  MValue minus() const {
    MValue res(p);
    mpfr_neg(res.m, m, GMP_RNDN);
    return res;
  }

  MValue minus(const MValue& b, mpfr_rnd_t round) const {
    MValue res(max(p, b.p));
    mpfr_sub(res.m, m, b.m, round);
    return res;
  }

  MValue times(const MValue& b, mpfr_rnd_t round) const {
    MValue res(max(p, b.p));
    mpfr_mul(res.m, m, b.m, round);
    return res;
  }

  MValue rcp(mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_d_div(res.m, 1.0, m, round);
    return res;
  }

  MValue sqrt(mpfr_rnd_t round) const {
    MValue res(p);
    mpfr_sqrt(res.m, m, round);
    return res;
  }

  MValue mid(const MValue& b) const {
    MValue res = plus(b, GMP_RNDN);
    mpfr_mul_d(res.m, res.m, 0.5, GMP_RNDN);
    return res;
  }

  int sign() const { return mpfr_sgn(m); }
  bool operator==(const MValue& b) const { return mpfr_equal_p(m, b.m) != 0; }
  bool operator!=(const MValue& b) const { return mpfr_equal_p(m, b.m) == 0; }
  bool operator<(const MValue& b) const { return mpfr_less_p(m, b.m); }
};

class MParameter {
  friend class PParameter;
  template <template <typename> class P>
  friend class Object;
  MValue l, u;
  Mods mods;
  bool alg, clear, op;

  MParameter(const MValue& l, const MValue& u, const Mods& mods, bool alg)
      : l(l), u(u), mods(mods), alg(alg), clear(true), op(true) {}

 public:
  static unsigned int nMixed, nRascal;

  static MParameter constant(double a) {
    return MParameter(MValue(a, curPrecision()), MValue(a, curPrecision()),
                      Mods(a), false);
  }

  MParameter() : clear(false), op(false) {}

  MParameter(const Parameter& a)
      : l(a.lb(), curPrecision()),
        u(a.ub(), curPrecision()),
        mods(a),
        alg(false),
        clear(true),
        op(true) {}

#ifdef EXACTRAT
  MParameter(const PParameter& a)
      : l(a.lb(), curPrecision()),
        u(a.ub(), curPrecision()),
        mods(a.lb()),
        alg(false),
        clear(true),
        op(true) {}
#else
  MParameter(const PParameter& a)
      : l(a.lb(), curPrecision()),
        u(a.ub(), curPrecision()),
        mods(a.mods),
        alg(a.alg),
        clear(true),
        op(true) {}
#endif
  MParameter(const MParameter& il, const MParameter& iu, bool ialg)
      : l(il.l.copy()),
        u(iu.u.copy()),
        mods(il.l, iu.u),
        alg(il.alg || iu.alg || ialg),
        clear(true),
        op(true) {}

  void setAlg() { alg = true; }

  MParameter(const MParameter& a) : mods(a.mods), alg(a.alg), op(false) {
    if (a.op || !a.clear) {
      l = a.l;
      u = a.u;
      clear = a.clear;
      ((MParameter&)a).clear = ((MParameter&)a).op = false;
    } else {
      l = a.l.copy();
      u = a.u.copy();
      clear = true;
    }
  }

  ~MParameter() {
    if (clear) {
      l.clear();
      u.clear();
    }
  }

  bool hasCurrentPrimes() { return mods.hasCurrentPrimes(); }

  MParameter& operator=(const MParameter& a) {
    if (clear) {
      l.clear();
      u.clear();
    }
    if (a.op || !a.clear) {
      l = a.l;
      u = a.u;
      clear = a.clear;
      ((MParameter&)a).clear = ((MParameter&)a).op = false;
    } else {
      l = a.l.copy();
      u = a.u.copy();
      clear = true;
    }
    mods = a.mods;
    alg = a.alg;
    op = false;
    return *this;
  }

  MParameter operator+(const MParameter& b) const {
    return MParameter(l.plus(b.l, GMP_RNDD), u.plus(b.u, GMP_RNDU),
                      mods + b.mods, alg || b.alg);
  }

  MParameter operator+(double b) const { return *this + constant(b); }

  MParameter operator-(const MParameter& b) const {
    return MParameter(l.minus(b.u, GMP_RNDD), u.minus(b.l, GMP_RNDU),
                      mods - b.mods, alg || b.alg);
  }

  MParameter operator-(double b) const { return *this - constant(b); }

  MParameter operator-() const {
    return MParameter(u.minus(), l.minus(), -mods, alg);
  }

  MParameter operator*(const MParameter& b) const {
    return u.sign() == -1 ? (-*this).times(-b) : times(b);
  }

  MParameter times(const MParameter& b) const {
    bool nalg = alg || b.alg;
    Mods nmods = mods * b.mods;
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

  MParameter operator*(double b) const { return *this * constant(b); }

  MParameter rcp() const {
    assert(sign());
    return MParameter(u.rcp(GMP_RNDD), l.rcp(GMP_RNDU), Mods(1) / mods, alg);
  }

  MParameter operator/(const MParameter& b) const { return *this * b.rcp(); }

  MParameter operator/(double b) const { return *this * constant(b).rcp(); }

  MParameter sqrt() const {
    int s = sign();
    assert(s != -1);
    // assert(sign() != -1);
    MValue nl = l.sqrt(MPFR_RNDD), nu = u.sqrt(MPFR_RNDU);
    return MParameter(nl, nu, Mods(nl, nu), true);
  }

  MParameter abs() const { return sign() == 1 ? *this : -*this; }
  double lb() const { return mods.zero() ? 0.0 : l.value(GMP_RNDD); }
  double ub() const { return mods.zero() ? 0.0 : u.value(GMP_RNDU); }
  double mid() const { return 0.5 * (lb() + ub()); }

  int sign(bool fail = true) const {
    int su = u.sign();
    if (su == -1) {
      if (mods.mixed()) {
        ++nMixed;
        mixedReport(1, mods);
      } else if (mods.zero())
        ++nRascal;
      return -1;
    }
    int sl = l.sign();
    if (sl == 1) {
      if (mods.mixed()) {
        ++nMixed;
        mixedReport(2, mods);
      } else if (mods.zero())
        ++nRascal;
      return 1;
    }
    if (!fail || (sl == 0 && su == 0) || mods.zero()) return 0;
    throw SignException(alg);
  }

  bool operator<(const MParameter& b) const { return (b - *this).sign() == 1; }
  bool operator>(const MParameter& b) const { return (b - *this).sign() == -1; }
  MParameter lbP() const {
    return MParameter(l.copy(), l.copy(), Mods(l, l), alg);
  }
  MParameter ubP() const {
    return MParameter(u.copy(), u.copy(), Mods(u, u), alg);
  }

  MParameter midP() const {
    MValue x = l.mid(u);
    ;
    return MParameter(x, x.copy(), Mods(x, x), alg);
  }

  bool subset(const MParameter& b) const {
    return !(l < b.l) && !(b.u < u) && (b.l < l || u < b.u);
  }

  bool intersects(const MParameter& b) const { return !(u < b.l || b.u < l); }

  MParameter interval(const MParameter& b) const {
    return MParameter(l.copy(), b.u.copy(), Mods(l, b.u), alg || b.alg);
  }

  MParameter intersect(const MParameter& b) const {
    assert(!(u < b.l || b.u < l));
    MValue nl = l < b.l ? b.l : l, nu = u < b.u ? u : b.u;
    return MParameter(nl.copy(), nu.copy(), Mods(nl, nu), alg || b.alg);
  }

  double distance(const MParameter& b) const {
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

  double intervalWidth() const {
    MValue d = u.minus(l, GMP_RNDN);
    double w = d.value();
    d.clear();
    return w;
  }
};

inline Parameter operator+(double a, const Parameter& b) {
  return Parameter::constant(a) + b;
}

inline Parameter operator-(double a, const Parameter& b) {
  return Parameter::constant(a) - b;
}

inline Parameter operator*(double a, const Parameter& b) {
  return Parameter::constant(a) * b;
}

inline Parameter operator/(double a, const Parameter& b) {
  return Parameter::constant(a) / b;
}

inline bool operator<(double a, const Parameter& b) {
  return Parameter::constant(a) < b;
}

inline bool operator>(double a, const Parameter& b) {
  return Parameter::constant(a) > b;
}

inline PParameter operator+(double a, const PParameter& b) {
  return PParameter::constant(a) + b;
}

inline PParameter operator-(double a, const PParameter& b) {
  return PParameter::constant(a) - b;
}

inline PParameter operator*(double a, const PParameter& b) {
  return PParameter::constant(a) * b;
}

inline PParameter operator/(double a, const PParameter& b) {
  return PParameter::constant(a) / b;
}

inline bool operator<(double a, const PParameter& b) {
  return PParameter::constant(a) < b;
}

inline bool operator>(double a, const PParameter& b) {
  return PParameter::constant(a) > b;
}

inline MParameter operator+(double a, const MParameter& b) {
  return MParameter::constant(a) + b;
}

inline MParameter operator-(double a, const MParameter& b) {
  return MParameter::constant(a) - b;
}

inline MParameter operator*(double a, const MParameter& b) {
  return MParameter::constant(a) * b;
}

inline MParameter operator/(double a, const MParameter& b) {
  return MParameter::constant(a) / b;
}

inline bool operator<(double a, const MParameter& b) {
  return MParameter::constant(a) < b;
}

inline bool operator>(double a, const MParameter& b) {
  return MParameter::constant(a) > b;
}

}  // namespace acp
#endif
