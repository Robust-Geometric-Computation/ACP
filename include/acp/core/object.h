// -*- mode: c++ -*-

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

#ifndef OBJECT
#define OBJECT

#include <limits.h>
#include <sys/time.h>
#include <map>
#include <set>
#include "acp/core/acp.h"

double getTime();

#define ALG_PREC 265u
//#define ALG_PREC 212u

namespace acp {

class RefCnt {
  template <class T>
  friend class PTR;
  int refCnt;

  void incRef() {
    pthread_mutex_lock(&mutex);
    refCnt++;
    pthread_mutex_unlock(&mutex);
  }

  void decRef() {
    pthread_mutex_lock(&mutex);
    bool flag = --refCnt == 0;
    pthread_mutex_unlock(&mutex);
    if (flag) delete this;
  }

 protected:
  pthread_mutex_t mutex;

 public:
  RefCnt() : refCnt(0) { mutex = PTHREAD_MUTEX_INITIALIZER; }
  virtual ~RefCnt() { assert(refCnt == 0); }
};

template <class T>
class PTR {
 public:
  T* t;  // TODO MAKE PRIVATE AGAIN
  PTR() : t(0) {}
  PTR(T* t) : t(t) { incRef(); }
  PTR(const PTR& p) : t(p.t) { incRef(); }
  const PTR& operator=(const PTR& p) {
    p.incRef();
    decRef();
    t = p.t;
    return *this;
  }
  ~PTR() { decRef(); }
  void incRef() const {
    if (t != 0) t->incRef();
  }
  void decRef() const {
    if (t != 0) t->decRef();
  }
  operator T*() const { return t; }
  T* operator->() const { return t; }
};

class BaseObject : public RefCnt {
 protected:
  static unsigned int deltaPrecision, maxPrecision;
  static int perIs[TMAX];
  static double perTs[TMAX], maxP;
  static bool perInts[TMAX];
  static int& perI() { return perIs[threadId()]; }
  static double& perT() { return perTs[threadId()]; }
  static bool& perInt() { return perInts[threadId()]; }

 public:
  static double wMin;
  virtual ~BaseObject() {}
  static void addThread(unsigned int i) {
    assert(i < TMAX);
    pthread_setspecific(idkey, (void*)(threadIds + i));
    threadIds[i] = i;
    curPrecisions[i] = 53u;
    perIs[i] = -1;
    perInts[i] = false;
  }
};

class PrecisionException : public std::exception {
 public:
  virtual const char* what() const throw() {
    return "Maximum precision exceeded";
  }
};

#define NHoms 2

#ifdef INPUT_COMPLEXITY
extern vector<set<void*>> inputSets;
class InputSetPusher {
 public:
  int index;
  InputSetPusher() {
    index = inputSets.size();
    inputSets.push_back(set<void*>());
  }

  InputSetPusher(bool input, const set<void*>& inputs) {
    if (input || inputs.size() > 0)
      index = -1;
    else {
      index = inputSets.size();
      inputSets.push_back(inputs);
    }
  }
  ~InputSetPusher() {
    if (index != -1) inputSets.pop_back();
  }
};
#endif

template <template <typename> class P>
class Object : public BaseObject {
  class Elt {
   public:
    double t, ti, *v;
    P<Parameter> p;
    P<MParameter>*mp, *mpi;

    Elt() : t(0.0), ti(0.0), v(0), mp(0), mpi(0) {}
    ~Elt() {
      delete[] v;
      delete mp;
      delete mpi;
    }

    void setV(unsigned int n) {
      v = new double[n];
      for (int i = 0; i < n; ++i) v[i] = randomNumber(-1.0, 1.0);
    }
  };

  P<Parameter> p;
  unsigned int pPrecision;
  P<PParameter>* pp;
  P<MParameter>* mp;

  set<void*> inputs;

  bool checkAccuracy(double acc) {
    if (acc == 1.0) return true;
    for (int i = 0; i < p.size(); i++) {
      double l = p[i].lb(), u = p[i].ub();
      if (u - l > acc * min(fabs(l), fabs(u))) {
        double ln = nextafter(l, 1.0 + u);
        if (ln < u) {
          ln = nextafter(ln, 1.0 + u);
          if (ln < u) return false;
        }
      }
    }
    return true;
  }

 protected:
  Elt* e;
  bool input, clamp;

  virtual P<Parameter> calculateP() {
    cerr << "missing calculate Parameter" << endl;
    assert(0);
    return P<Parameter>();
  }

  virtual P<PParameter> calculatePP() {
    cerr << "missing calculate PParameter" << endl;
    assert(0);
    return P<PParameter>();
  }

  virtual P<MParameter> calculateMP() {
    cerr << "missing calculate MParameter" << endl;
    assert(0);
    return P<MParameter>();
  }

  void noClear(P<MParameter>& x) const {
    for (int i = 0; i < x.size(); ++i) {
      if (!x[i].clear) {
        x[i].l = x[i].l.copy();
        x[i].u = x[i].u.copy();
      }
      x[i].clear = false;
    }
  }

  void clear(P<MParameter>* x) const {
    if (x) {
      for (int i = 0; i < x->size(); ++i) (*x)[i].clear = true;
      delete x;
    }
  }

 public:
  Object() : pPrecision(0), input(false), clamp(false), pp(0), mp(0), e(0) {}

  Object(const P<Parameter>& p, bool clamp = false)
      : p(p), pPrecision(53u), input(true), clamp(clamp), pp(0), mp(0), e(0) {
    for (int i = 0; i < p.size(); i++)
      if (p[i].lb() != p[i].ub()) {
        std::cerr << "Input object has non-trivial interval." << std::endl;
        assert(0);
      }
  }

  ~Object() {
    delete pp;
    clear(mp);
    delete[] e;
  }
  bool uninitialized() const { return pPrecision == 0; }
  P<Parameter> getCurrentP() const { return p; }

  template <class N>
  P<N> get() {
#ifdef INPUT_COMPLEXITY
    InputSetPusher isp(input, inputs);
#endif

    try {
      pthread_mutex_lock(&mutex);
      P<N> ne = get1<N>();
      pthread_mutex_unlock(&mutex);
#ifdef INPUT_COMPLEXITY
      if (input) {
        for (int i = 0; i < inputSets.size(); i++) inputSets.at(i).insert(this);
      } else {
        if (isp.index != -1) inputs = inputSets.at(isp.index);
      }
#endif
      return ne;
    } catch (SignException se) {
      pthread_mutex_unlock(&mutex);
      throw se;
    } catch (unsigned int p) {
      pthread_mutex_unlock(&mutex);
      throw p;
    }
  }

  template <class N>
  P<N> get1() {
    int cp = curPrecision();
    if (cp == 53u) return get53();
    if (cp == 100u) {
      if (pp) {
        if ((*pp)[0].hasCurrentPrimes()) return *pp;
        delete pp;
        pp = 0;
        if (mp == 0) pPrecision = 53u;
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
      if (uninitialized() || !input) p = *pp;
      pPrecision = 100u;
      return *pp;
    }
    return perI() == -1 || clamp ? getM() : getMP();
  }

  P<Parameter> get53() {
    if (uninitialized()) {
      p = calculateP();
      pPrecision = 53u;
    }
    return p;
  }

  P<MParameter> getM() {
    if (pPrecision < curPrecision()) {
      P<MParameter>* newMp =
          input ? new P<MParameter>(p) : new P<MParameter>(calculateMP());
      noClear(*newMp);
      clear(mp);
      mp = newMp;
      if (!input) p = *newMp;
      pPrecision = curPrecision();
    }
    return *mp;
  }

  P<MParameter> getMP() {
    if (!e) e = new Elt[NHoms];
    Elt& ei = e[perI()];
    double t = perT();
    if (perInt()) {
      if (ei.mpi && ei.ti == t) return *ei.mpi;
      P<MParameter> mpi = input ? getMPI() : calculateMP();
      ei.ti = t;
      delete ei.mpi;
      ei.mpi = new P<MParameter>(mpi);
      return *ei.mpi;
    }
    if (ei.mp && ei.t == t) return *ei.mp;
    P<MParameter> mp = input ? getMPI() : calculateMP();
    if (!input) {
      int peri = perI();
      perI() = -1;
      P<MParameter> mp0 = get1<MParameter>();
      perI() = peri;
      for (int i = 0; i < mp.size(); ++i) {
        double d = mp[i].distance(mp0[i]);
        if (d > 0.0) {
          double iw = max(mp[i].intervalWidth(), mp0[i].intervalWidth());
          wMin = min(wMin, d / iw);
        }
      }
    }
    ei.t = t;
    delete ei.mp;
    ei.mp = new P<MParameter>(mp);
    return *ei.mp;
  }

  P<MParameter> getMPI() {
    Elt& ei = e[perI()];
    if (!ei.v) ei.setV(p.size());
    double t = perT();
    P<MParameter> q(p);  // q; // How could this have worked?
    if (perInt())
      for (int i = 0; i < p.size(); ++i) {
        MParameter a = MParameter::constant(p[i].lb()),
                   b = a + t * MParameter::constant(ei.v[i]);
        q[i] =
            ei.v[i] > 0.0 ? MParameter(a, b, false) : MParameter(b, a, false);
      }
    else
      for (int i = 0; i < p.size(); ++i)
        q[i] =
            MParameter::constant(p[i].lb()) + t * MParameter::constant(ei.v[i]);
    return q;
  }

  const P<Parameter>& getApprox(double acc = 1e-17) {
    unsigned int peri = perI();
    perI() = -1;
    try {
      get<Parameter>();
      if (checkAccuracy(acc)) {
        perI() = peri;
        return p;
      }
    } catch (SignException se) {
    }

#ifdef BLEEN
    curPrecision() = 100u;
    try {
      get<PParameter>();
      if (checkAccuracy(acc)) {
        perI() = peri;
        curPrecision() = 53u;
        return p;
      }
    } catch (SignException se) {
    }
#endif

    curPrecision() = 212u;
    while (curPrecision() <= maxPrecision) {
      try {
        get<MParameter>();
        if (checkAccuracy(acc)) {
          curPrecision() = 53u;
          perI() = peri;
          return p;
        }
      } catch (SignException se) {
      }
      curPrecision() += deltaPrecision;
    }
    cerr << "warning: checkAccuracy failure." << endl;
    curPrecision() = 53u;
    perI() = peri;
    return p;
  }
};

class Primitive : public BaseObject {
  virtual Parameter calculateP() = 0;
  virtual PParameter calculatePP() = 0;
  virtual MParameter calculateMP() = 0;

 public:
  static unsigned long int nPrim, nAmb, nDeg, nAlg, nHom, nId, maxBits;
  static double tPar, tPPar, tMPar, vMPar;
  operator int() {
#ifdef INPUT_COMPLEXITY
    InputSetPusher isp;
    static int maxInputs;
#endif
    ++nPrim;
    double t = getTime();
    try {
      int s = calculateP().sign();
      tPar += getTime() - t;
#ifdef INPUT_COMPLEXITY
      if (maxInputs < inputSets.at(isp.index).size()) {
        cerr << "maxInputs " << inputSets.at(isp.index).size() << endl;
        maxInputs = inputSets.at(isp.index).size();
      }
#endif
      return s;
    } catch (SignException se) {
    }
    double nt = getTime();
    tPar += nt - t;
    t = nt;
    curPrecision() = 100u;
    ++nAmb;
    bool alg = false;
    while (true) {
      try {
        PParameter x = calculatePP();
        int s = x.sign();
        if (s == 0) ++nDeg;
#ifdef EXACTRAT
        else {
          int xb = x.numBits();
          if (xb > maxBits) {
            maxBits = xb;
            cerr << "maxBits " << xb << endl;
          }
        }
        // maxBits = max(maxBits, x.numBits());
#endif
        curPrecision() = 53u;
        tPPar += getTime() - t;
#ifdef INPUT_COMPLEXITY
        if (maxInputs < inputSets.at(isp.index).size()) {
          cerr << "maxInputs " << inputSets.at(isp.index).size() << endl;
          maxInputs = inputSets.at(isp.index).size();
        }
#endif
        return s;
      } catch (unsigned int p) {
        Mods::changePrime(p);
      } catch (SignException se) {
#ifdef SET_ALG
        if (se.alg) {
#else
        if (false) {
#endif
          alg = true;
          ++nAlg;
          int s = 2;
          bool hom = homotopy(s);
          if (s == 0) ++nDeg;
          if (hom || s != 2) {
            curPrecision() = 53u;
            tMPar += getTime() - t;
#ifdef INPUT_COMPLEXITY
            if (maxInputs < inputSets.at(isp.index).size()) {
              cerr << "maxInputs " << inputSets.at(isp.index).size() << endl;
              maxInputs = inputSets.at(isp.index).size();
            }
#endif
            return hom ? 0 : s;
          }
        }
        break;
      }
    }
    curPrecision() = !alg ? 212u : ALG_PREC;
    static unsigned int maxPrecisionUsed = maxPrecision;
    while (!alg || curPrecision() <= maxPrecision) {
      try {
        MParameter x = calculateMP();
        int s = x.sign();
        if (s == 0)
          ++nDeg;
        else if (alg)
          vMPar = min(vMPar, (s == 1 ? x.lb() : -x.ub()) / x.intervalWidth());
        curPrecision() = 53u;
        tMPar += getTime() - t;
#ifdef INPUT_COMPLEXITY
        if (maxInputs < inputSets.at(isp.index).size()) {
          cerr << "maxInputs " << inputSets.at(isp.index).size() << endl;
          maxInputs = inputSets.at(isp.index).size();
        }
#endif
        return s;
      } catch (SignException se) {
#ifdef NO_HOMOTOPY
        if (curPrecision() > 424u) {
          curPrecision() = 53u;
          return 0;
        }
#endif
        curPrecision() += deltaPrecision;
        if (curPrecision() > maxPrecisionUsed) {
          maxPrecisionUsed = curPrecision();
          cerr << "max precision used = " << maxPrecisionUsed << endl;
        }
      }
    }
    throw PrecisionException();
    return 0;
  }

  bool homotopy(int& s) {
    s = 2;
    bool mixed = false;
    curPrecision() = ALG_PREC;
    while (true) {
      try {
        MParameter x = calculateMP();
        int sx = x.sign(false);
        if (sx || (x.lb() == 0.0 && x.ub() == 0.0)) {
          s = sx;
          if (s)
            vMPar = min(vMPar, (s == 1 ? x.lb() : -x.ub()) / x.intervalWidth());
          return false;
        }
        break;
      } catch (SignException se) {
        curPrecision() += deltaPrecision;
        if (curPrecision() > maxPrecision) throw PrecisionException();
      }
    }
    ++nHom;
    int n = 0;
    for (int i = 0; i < NHoms; ++i) {
      perI() = i;
      if (homotopyAux()) ++n;
    }
    perI() = -1;
    if (n > 0 && n < NHoms) cerr << "mixed homotopy" << endl;
    if (n == 0) return false;
    ++nId;
    return true;
  }

  bool homotopyAux() {
    perT() = maxP;
    while (true) {
      try {
        if (calculateMP().sign(false) == 0) return true;
        break;
      } catch (SignException se) {
        perT() *= 0.1;
      }
    }
    double ts = perT();
    int k = 0, kmax = round(deltaPrecision * 0.30102999566398117);
    while (true) {
      try {
        perInt() = true;
        calculateMP();
        perInt() = false;
        if (perT() < ts) {
          assert(perT() > 0.0);
          return calculateMP().sign(false) == 0;
        }
        return false;
      } catch (SignException se) {
        ++k;
        if (k < kmax)
          perT() *= 0.1;
        else {
          k = 0;
          curPrecision() += deltaPrecision;
        }
      }
    }
  }
};

#define DeclareCalculate(P)                                       \
 private:                                                         \
  P<Parameter> calculateP() { return calculate<Parameter>(); }    \
  P<PParameter> calculatePP() { return calculate<PParameter>(); } \
  P<MParameter> calculateMP() { return calculate<MParameter>(); } \
  template <class N>                                              \
  P<N> calculate()

#define DeclareSign                                            \
 private:                                                      \
  Parameter calculateP() { return calculate<Parameter>(); }    \
  PParameter calculatePP() { return calculate<PParameter>(); } \
  MParameter calculateMP() { return calculate<MParameter>(); } \
  template <class N>                                           \
  N calculate()

void primitiveReport();
void resetReport();
}  // namespace acp
#endif
