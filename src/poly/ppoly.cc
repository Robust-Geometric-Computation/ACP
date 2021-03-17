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

#include "acp/poly/ppoly.h"
#include <vector>
using namespace std;
using namespace acp;

namespace acp {

#ifdef PPOLY1
bool PPoly1::print_all = false;

Parameter PPoly1::polish(const Parameter& x) const {
  Parameter xubP = x.ubP();
  int sub = value(xubP).sign();
  return shrinkBracket(x, sub);
}

PPoly1::PPoly1(int d, Parameter* a) : d(d) {
  this->a = new Parameter[d + 1];
  if (a)
    for (int i = 0; i <= d; ++i) this->a[i] = a[i];
}

PPoly1::PPoly1(const PPoly1& f) : d(f.d) {
  a = new Parameter[d + 1];
  for (int i = 0; i <= d; ++i) a[i] = f.a[i];
}

PPoly1::PPoly1(double c, int d) : d(d) {
  this->a = new Parameter[d + 1];
  for (int i = 0; i < d; i++) a[i] = 0;
  a[d] = c;
}

PPoly1::PPoly1(const Parameter& c, int d) : d(d) {
  this->a = new Parameter[d + 1];
  for (int i = 0; i < d; i++) a[i] = 0;
  a[d] = c;
}

PPoly1& PPoly1::operator=(const PPoly1& f) {
  delete[] a;
  d = f.d;
  a = new Parameter[d + 1];
  for (int i = 0; i <= d; ++i) a[i] = f.a[i];
  return *this;
}

Parameter PPoly1::value(const Parameter& x) const {
  Parameter y = a[d];
  if (print_all) printf("a[%d]: [%.16e, %.16e]\n", d, a[d].lb(), a[d].ub());
  if (print_all)
    printf("yInterval / ylb : %.16e\n", y.intervalWidth() / y.lb());
  for (int i = d - 1; i >= 0; --i) {
    if (print_all) printf("a[%d]: [%.16e, %.16e]\n", i, a[i].lb(), a[i].ub());
    y = x * y + a[i];
    if (print_all) printf("y : [%.16e, %.16e]\n", y.lb(), y.ub());
    if (print_all)
      printf("yInterval / ylb : %.16e\n", y.intervalWidth() / y.lb());
  }
  return y;
}

Parameter PPoly1::der(const Parameter& x) const {
  Parameter y = d * a[d];
  for (int i = d - 1; i > 0; --i) y = x * y + i * a[i];
  return y;
}

PPoly1 PPoly1::der() const {
  PPoly1 g(d - 1);
  for (int i = 1; i <= d; ++i) g.a[i - 1] = i * a[i];
  return g;
}

// f * g = (F1 xd + F0) * (G1 xd + G0) = F1*G1 x2d + (F0*G1 + F1*G0)xd + F0*G0

PPoly1 PPoly1::karatsuba(const PPoly1& g) const {
  if (d != g.d) {
    printf("(karatsuba) defaulting to naive mult\n");
    return *this * g;
  }

  if (d <= 5) {
    return *this * g;
  }

  int fd = d / 2;

  int f1d = d - fd;

  int f0d = fd - 1;

  PPoly1 fp(1, fd);

  PPoly1 F1(f1d, &a[fd]);
  PPoly1 F0(f0d, &a[0]);

  PPoly1 G1(f1d, &g.a[fd]);
  PPoly1 G0(f0d, &g.a[0]);

  PPoly1 F1G1 = F1.karatsuba(G1);
  PPoly1 F0G0 = F0.karatsuba(G0);

  PPoly1 F1F0G1G0 = (F1 + F0).karatsuba(G1 + G0);

  return (F1.karatsuba(G1) * fp * fp) + ((F1F0G1G0 - F1G1 - F0G0) * fp) +
         (F0G0);
}

PPoly1 PPoly1::normalized() const {
  PPoly1 cop(*this);

  if (d >= 0) {
    Parameter max = cop.a[0];
    for (int i = 1; i <= d; i++) {
      max = (max < cop.a[i] ? cop.a[i] : max);
    }

    for (int i = 0; i <= d; i++) {
      cop.a[i] = cop.a[i] / max;
    }
  }

  return cop;
}

vector<Parameter> PPoly1::rootsR(const Parameter& lb,
                                 const Parameter& ub) const {
  if (d == 2) return roots2(lb, ub);
  if (d < 2) return roots1(lb, ub);
  PPoly1 df = der();
  vector<Parameter> dres = df.roots(lb, ub);
  int dn = dres.size();
  vector<Parameter> res;
  Parameter lbi = lb;
  int slbi = value(lbi).sign();
  for (int i = 0; i <= dn; ++i) {
    Parameter ubi = i < dn ? dres[i] : ub;
    int subi = value(ubi).sign();
    if (slbi == -subi) {
      Parameter x = lbi.innerInterval(ubi);
      res.push_back(shrinkBracket(x, subi));
    }
    lbi = ubi;
    slbi = subi;
  }
  return res;
}

vector<Parameter> PPoly1::roots1(const Parameter& lb,
                                 const Parameter& ub) const {
  vector<Parameter> res;
  if (d == 0) return res;
  if (d == 1) {
    Parameter r = -a[0] / a[1];
    if (lb < r && r < ub) res.push_back(r);
    return res;
  }
  assert(0);
  return res;
}

vector<Parameter> PPoly1::roots2(const Parameter& lb,
                                 const Parameter& ub) const {
  vector<Parameter> res;
  const Parameter &A = a[2], &B = a[1], &C = a[0];
  Parameter D2 = B * B - 4.0 * A * C;
  if (D2.sign() == -1) return res;
  Parameter D = D2.sqrt();
  Parameter r[2];

  // (-B - D)/2A < (-B + D)/2A if A>0
  bool Apos = (A.sign() == 1);
  if (B.sign() == -1) {
    Parameter DmB = D - B;
    r[!Apos] = 2.0 * C / DmB;
    r[Apos] = 0.5 * DmB / A;
  } else {
    Parameter DpB = D + B;
    r[!Apos] = -0.5 * DpB / A;
    r[Apos] = -2.0 * C / DpB;
  }

  if (lb < r[0] && r[0] < ub) res.push_back(r[0]);
  if (lb < r[1] && r[1] < ub) res.push_back(r[1]);

  return res;
}

Parameter PPoly1::shrinkBracket(Parameter x, int sub) const {
  static int count;
  count++;
  bool bflag = false;
  while (true) {
    bool flag = true;
    Parameter xm = x.midP(), fm = value(xm), df = der(x);
    if (df.sign(false) == 0)
      flag = false;
    else {
      Parameter nx = (xm - fm / df).intersect(x);
      if (nx.subset(x))
        x = nx;
      else
        flag = false;
    }
    if (!flag) {
      int sm = fm.sign(false);
      if (sm == sub) {
        Parameter nx = x.interval(xm);
        if (!nx.subset(x)) break;
        x = nx;
      } else if (sm == -sub) {
        Parameter nx = xm.interval(x);
        if (!nx.subset(x)) break;
        x = nx;
      } else if (!bflag) {
        bflag = true;
        x = shrinkBracketB(x, sub);
      } else
        break;
    }
  }
  return x;
}

Parameter PPoly1::shrinkBracketB(const Parameter& x, int sub) const {
  Parameter xlb = x.lbP(), xm = x.midP(), xub = x.ubP();
  while ((xlb - xm).sign(false) < 0) {
    Parameter nx = xlb.interval(xm).midP();
    if ((nx - xlb).sign(false) == 0 || (nx - xm).sign(false) == 0) break;
    int sm = value(nx).sign(false);
    if (sm == 0)
      xm = nx;
    else if (sm == sub)
      xub = nx;
    else
      xlb = nx;
  }
  xm = x.midP();
  while ((xm - xub).sign(false) < 0) {
    Parameter nx = xm.interval(xub).midP();
    if ((nx - xm).sign(false) == 0 || (nx - xub).sign(false) == 0) break;
    int sm = value(nx).sign(false);
    if (sm == 0)
      xm = nx;
    else if (sm == -sub)
      xlb = nx;
    else
      xub = nx;
  }
  Parameter nx = xlb.interval(xub);
  return nx.subset(x) ? nx : x;
}

PPoly1 PPoly1::operator*(const PPoly1& g) const {
  if (d == -1 || g.d == -1) {
    return PPoly1();
  }

  int degree = d + g.d;

  vector<Parameter> coef(degree + 1);

  for (int i = 0; i <= degree; i++) {
    coef[i] = 0;
  }

  for (int i = 0; i <= d; i++) {
    for (int j = 0; j <= g.d; j++) {
      coef[i + j] = coef[i + j] + a[i] * g.a[j];
    }
  }

  return PPoly1(degree, &coef[0]);
}

PPoly1 PPoly1::operator*(const double p) const {
  Parameter* ca = new Parameter[d + 1];

  for (int i = 0; i <= d; i++) {
    ca[i] = a[i] * p;
  }

  return PPoly1(d, ca);
}

PPoly1 PPoly1::operator*(const Parameter& p) const {
  Parameter* ca = new Parameter[d + 1];

  for (int i = 0; i <= d; i++) {
    ca[i] = a[i] * p;
  }

  return PPoly1(d, ca);
}

PPoly1 PPoly1::plusCtimes(double c, const PPoly1& b) const {
  int degree = max(d, b.d);

  // Parameter coef[degree+1];
  vector<Parameter> coef(degree + 1);

  for (int i = 0; i <= degree; i++) {
    coef[i] = (d > b.d ? a[i] : c * b.a[i]);
  }

  for (int i = 0; i <= min(d, b.d); i++) {
    coef[i] = coef[i] + (d <= b.d ? a[i] : c * b.a[i]);
  }

  return PPoly1(degree, &coef[0]);
}

PPoly1 PPoly1::plusCtimes(const Parameter& c, const PPoly1& b) const {
  int degree = max(d, b.d);

  vector<Parameter> coef(degree + 1);

  for (int i = 0; i <= degree; i++) {
    coef[i] = (d > b.d ? a[i] : c * b.a[i]);
  }
  /*
    printf("before: [");
    for(int i = 0; i <= degree; i++) {
      printf(" %.3f ", coef[i].mid());
    }
    printf("]\n");
  */

  for (int i = 0; i <= min(d, b.d); i++) {
    coef[i] = coef[i] + (d <= b.d ? a[i] : c * b.a[i]);
  }
  /*
    printf("after: [");
    for(int i = 0; i <= degree; i++) {
      printf(" %.3f ", coef[i].mid());
    }
    printf("]\n");
  */

  return PPoly1(degree, &coef[0]);
}

PPoly1 PPoly1::quotient(const PPoly1& divisor, PPoly1& remainder) {
  if (divisor.d > d) {
    remainder = *this;
    return PPoly1();
  }
  PPoly1 f = *this;
  PPoly1 g = divisor.monic();
  PPoly1 q(d - g.d);
  for (int i = q.d; i >= 0; i--) {
    q.a[i] = f.a[g.d + i];
    for (int j = 0; j < g.d; j++) f.a[i + j] = f.a[i + j] - q.a[i] * g.a[j];
  }
  f.d = g.d - 1;
  while (f.d >= 0 && f.a[f.d].sign() == 0) f.d--;
  remainder = f;
  return q;
}

PPoly1 PPoly1::gcd(const PPoly1& b) {
  if (d == -1) return b;
  if (b.d == -1) return *this;
  PPoly1 f = *this;
  PPoly1 g = b;
  while (g.d != -1) {
    PPoly1 r;
    f.quotient(g, r);
    f = g;
    g = r;
  }
  return f;
}

void PPoly1::print() {
  if (d < 0) {
    printf("0");
  } else {
    for (int i = d; i >= 0; i--) {
      printf("%s%.4f%s", i == d ? "" : " + ", a[i].mid(),
             (i == 0 ? ""
                     : (i == 1 ? " t"
                               : (" t^{" + std::to_string(i) + "}").c_str())));
    }
    // printf("\n");
  }
}

PPoly1 PPoly1::monic() const {
  PPoly1 p(d);
  for (unsigned int i = 0u; i < d; ++i) p.a[i] = a[i] / a[d];
  p.a[d] = 1;
  return p;
}

PPoly1 PPoly1::normal() const {
  PPoly1 p(d);

  int maxi = 0;
  for (unsigned i = 1u; i <= d; i++) {
    maxi = (a[i].abs() < a[maxi].abs() ? maxi : i);
  }

  double mid = a[maxi].mid();

  for (unsigned int i = 0u; i <= d; ++i) p.a[i] = a[i] / mid;

  return p;
}

PPoly1 PPoly1::shift(double s) const {
  PPoly1 p(*this);
  for (int i = 0; i < d; i++)
    for (int j = d - 1; j >= i; j--) p.a[j] = p.a[j] + p.a[j + 1] * s;
  return p;
}

PPoly1 PPoly1::shift(const Parameter& s) const {
  PPoly1 p(*this);
  for (int i = 0; i < d; i++)
    for (int j = d - 1; j >= i; j--) p.a[j] = p.a[j] + p.a[j + 1] * s;
  return p;
}

PPoly1 PPoly1::neg() const {
  PPoly1 p(d);
  for (unsigned int i = 0u; i <= d; ++i) p.a[i] = i % 2 == 0 ? a[i] : -a[i];
  return p;
}

PPoly1 PPoly1::dual() const {
  PPoly1 p(d);
  for (unsigned int i = 0u; i <= d; ++i) p.a[i] = a[d - i];
  return p;
}

template <class N>
bool zerop(const N& p) {
  return p.sign() == 0;
}

PPoly1 PPoly1::removeZeroRoots() const {
  unsigned int i = 0u;
  while (i < d && zerop<Parameter>(a[i])) ++i;
  PPoly1 p(d - i);
  for (unsigned int j = 0u; j <= d - i; ++j) p.a[j] = a[i + j];
  return p;
}

unsigned int PPoly1::descartes() const {
  bool pos = a[d] > 0.0;
  unsigned int s = 0u;
  for (int i = d - 1u; i >= 0; --i)
    if (!zerop(a[i]) && (a[i] > 0.0) != pos) {
      ++s;
      pos = !pos;
    }
  return s;
}

Parameter PPoly1::ubk() const {
  PPoly1 p = monic();
  Parameter u = 0;
  for (unsigned int k = 1u; k <= d; ++k)
    if (p.a[d - k] < 0.0) u = Parameter::max(u, (-p.a[d - k]).root(k));
  return 2.0 * u;
}

class Mobius {
  Mobius(double aa, double bb, double cc, double dd)
      : a(Parameter::constant(aa)),
        b(Parameter::constant(bb)),
        c(Parameter::constant(cc)),
        d(Parameter::constant(dd)) {}

  Mobius(const Parameter& a, const Parameter& b, const Parameter& c,
         const Parameter& d)
      : a(a), b(b), c(c), d(d) {}

 public:
  static Mobius identity() { return Mobius(1, 0, 0, 1); }

  static Mobius shift(double s) {
    return Mobius::shift(Parameter::constant(s));
  }

  static Mobius shift(const Parameter& s) { return Mobius(1, s, 0, 1); }

  static Mobius dual() { return Mobius(0, 1, 1, 0); }

  Mobius operator*(const Mobius& m) const;
  Parameter value(double x) const { return (a * x + b) / (c * x + d); }
  Parameter value(const Parameter& x) const {
    return (a * x + b) / (c * x + d);
  }
  Parameter inf() const { return a / c; }

  Parameter a, b, c, d;
};

Mobius Mobius::operator*(const Mobius& m) const {
  // (a (ma x + mb) / (mc x + md) + b) / (c (ma x + mb) / (mc x + md) + d)
  // (a (ma x + mb)  + b (mc x + md)) / (c (ma x + mb) + d (mc x + md))
  // ((a ma + b mc) x + (a mb + b md))) / ((c ma + d mc) x + (c mb + d md))
  Parameter na = a * m.a + b * m.c, nb = a * m.b + b * m.d,
            nc = c * m.a + d * m.c, nd = c * m.b + d * m.d;
  return Mobius(na, nb, nc, nd);
}

Ints PPoly1::rootInts(const Parameter& lb, const Parameter& ub) const {
  PPoly1 p = shift(lb);
  Mobius m = Mobius::shift(lb);
  Parameter s = 1 / (ub - lb);
  p = p.dual().shift(s);
  m = m * Mobius::dual() * Mobius::shift(s);
  Ints ints;
  vas(p, m, ints);
  return ints;
}

void PPoly1::vas(PPoly1 p, Mobius m, Ints& ints) const {
  if (zerop(p.a[0])) {
    Parameter r = m.value(0.0);
    ints.push_back(Int(r, r));
    p = p.removeZeroRoots();
  }
  unsigned int s = p.descartes();
  //// cout << "descartes " << s << endl;
  if (s <= 1u) {
    if (s == 1u) ints.push_back(vasInt(p, m));
    return;
  }
  // unsigned int lb = floor(1.0/p.dual().ubk());
  Parameter lb = (1.0 / p.dual().ubk()).lbP();
  //// cout << "shift " << lb.lb() << endl;
  if (lb > 0) {
    p = p.shift(lb);
    m = m * Mobius::shift(lb);
  }

  PPoly1 p1i = p.shift(1.0);
  Mobius m1i = m * Mobius::shift(1.0);
  //// cout << "try 1 inf (" << endl;
  vas(p1i, m1i, ints);
  //// cout << ") end 1 inf" << endl;

  if (p1i.descartes() < s) {
    PPoly1 p01 = p.dual().shift(1.0);
    Mobius m01 = m * Mobius::dual() * Mobius::shift(1.0);
    //// cout << "try 0 1 (" << endl;
    vas(p01, m01, ints);
    //// cout << ") end 1 inf" << endl;
  }
  if (zerop(p.value(1))) {
    Parameter r = m.value(1);
    ints.push_back(Int(r, r));
  }
}

Int PPoly1::vasInt(const PPoly1& p, const Mobius& m) const {
  if (zerop(m.d)) return Int(m.inf(), ubk());
  if (zerop(m.c)) return Int(m.value(0.0), ubk());
  Parameter m1 = m.value(0.0);
  Parameter m2 = m.inf();
  if (m1 < m2)
    return Int(m1, m2);
  else
    return Int(m2, m1);
}

vector<Parameter> PPoly1::rootsD(const Parameter& lb,
                                 const Parameter& ub) const {
  value(lb).sign();
  value(ub).sign();
  if (d == 2) return roots2(lb, ub);
  if (d < 2) return roots1(lb, ub);
  Ints ints = rootInts(lb, ub);
  vector<Parameter> res;
  for (Ints::iterator i = ints.begin(); i != ints.end(); ++i) {
    Parameter& lbi = i->first;
    Parameter& ubi = i->second;
    int slbi = value(lbi).sign();
    int subi = value(ubi).sign();
    Parameter x = lbi.innerInterval(ubi);
    res.push_back(shrinkBracket(x, subi));
  }
  // std::sort(res.begin(), res.end());
  for (int n = 1; n < res.size(); n++) {
    Parameter x = res[n];
    int i = n;
    while (i > 0 && res[i - 1] > x) {
      res[i] = res[i - 1];
      i--;
    }
    res[i] = x;
  }
  return res;
}
#endif
}  // namespace acp
