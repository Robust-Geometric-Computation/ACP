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

#ifndef PPOLY
#define PPOLY

#include <unordered_map>
#include <vector>
#include "acp/linmath/pv.h"
using namespace std;
using namespace acp;

namespace acp {

template <class T>
class PPoly {
 public:
  PPoly() {}
  PPoly(const T& c) { add(c, 0); }
  PPoly(const vector<T>& v) {
    for (int i = 0; i < v.size(); i++) add(v[i], i);
  }

  template <class N>
  PPoly(const PPoly<N>& p) {
    for (int i = 0; i < p.a.size(); i++) a.push_back(p.a[i]);
  }

  int size() const { return a.size(); }
  const T& operator[](int i) const { return a[i]; }
  T& operator[](int i) { return a[i]; }
  T lc() const { return *a.rbegin(); }
  int deg() const { return a.size() - 1; }
  int degree() const { return deg(); }
  bool zero() const { return a.empty(); }

  void add(const T& c, int e) {
    int d = deg();
    if (e <= d)
      a[e] = a[e] + c;
    else {
      while (d + 1 < e) {
        a.push_back(T(0));
        ++d;
      }
      a.push_back(c);
    }
  }

  PPoly operator+(const PPoly& p) const {
    int d = deg(), e = p.deg();
    PPoly q;
    if (d <= e) {
      for (int i = 0; i <= d; ++i) q.add(a[i] + p.a[i], i);
      for (int i = d + 1; i <= e; ++i) q.add(p.a[i], i);
    } else {
      for (int i = 0; i <= e; ++i) q.add(a[i] + p.a[i], i);
      for (int i = e + 1; i <= d; ++i) q.add(a[i], i);
    }
    return q;
  }

  PPoly operator-(const PPoly& p) const { return (*this) + (-p); }

  PPoly operator-() const {
    PPoly q;
    for (int i = 0; i <= deg(); ++i) q.add(-a[i], i);
    return q;
  }

  PPoly operator*(const PPoly& p) const {
    int d = deg(), e = p.deg();
    PPoly q;
    for (int i = 0; i <= d; ++i)
      for (int j = 0; j <= e; ++j) q.add(a[i] * p.a[j], i + j);
    return q;
  }

  PPoly rem(const PPoly& p, PPoly& q) const {
    T lp = p.lc();
    int dp = p.deg();
    PPoly r(*this);
    while (dp <= r.deg()) {
      T k = r.lc() / lp;
      int dk = r.deg() - dp;
      q.add(k, dk);
      r.a.pop_back();
      for (int i = 0; i < dp; ++i) r.add(-k * p.a[i], i + dk);
      r.checkDegree();
    }
    return r;
  }

  PPoly quotient(const PPoly& p, PPoly& r) const {
    PPoly q;
    r = rem(p, q);
    return q;
  }

  void checkDegree() {
    while (!zero() && zerop(lc())) a.pop_back();
  }

  PPoly gcd(const PPoly& p, PPoly& s, PPoly& t) const {
    s = PPoly(T(1));
    t = PPoly(T(0));
    PPoly r(*this), nr(p), ns, nt(T(1));
    r.checkDegree();
    nr.checkDegree();
    while (!nr.zero()) {
      PPoly q, nnr = r.rem(nr, q), nns = s - q * ns, nnt = t - q * nt;
      r = nr;
      nr = nnr;
      s = ns;
      ns = nns;
      t = nt;
      nt = nnt;
    }
    return r;
  }

  PPoly gcd(const PPoly& p) const {
    PPoly s, t;
    return gcd(p, s, t);
  }

  T value(const T& x) const {
    int d = deg();
    T y = a[d];
    for (int i = d - 1; i >= 0; --i) y = x * y + a[i];
    return y;
  }

  T der(const T& x) const {
    int d = deg();
    T y = d * a[d];
    for (int i = d - 1; i > 0; --i) y = x * y + i * a[i];
    return y;
  }

  PPoly der() const {
    PPoly g;
    for (int i = 1; i <= deg(); ++i) g.a.push_back(i * a[i]);
    return g;
  }

  PPoly neg() const {
    PPoly p;
    for (unsigned int i = 0u; i <= deg(); ++i)
      p.a.push_back(i % 2 == 0 ? a[i] : -a[i]);
    return p;
  }

  PPoly dual() const {
    PPoly p;
    int d = deg();
    for (unsigned int i = 0u; i <= d; ++i) p.a.push_back(a[d - i]);
    return p;
  }

  PPoly combine(const PPoly& p, T s) const {
    PPoly q;
    for (int i = 0; i <= deg(); ++i) q.add((1.0 - s) * a[i] + s * p.a[i], i);
    return q;
  }

  T shrinkBracket(T x, int sub) const {
    bool bflag = false;
    while (true) {
      bool flag = true;
      T xm = x.midP(), fm = value(xm), df = der(x);
      if (df.sign(false) == 0)
        flag = false;
      else {
        T nx = (xm - fm / df).intersect(x);
        if (nx.subset(x))
          x = nx;
        else
          flag = false;
      }
      if (!flag) {
        int sm = fm.sign(false);
        if (sm == sub) {
          T nx = x.interval(xm);
          if (!nx.subset(x)) break;
          x = nx;
        } else if (sm == -sub) {
          T nx = xm.interval(x);
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

  T shrinkBracketB(const T& x, int sub) const {
    T xlb = x.lbP(), xm = x.midP(), xub = x.ubP();
    while ((xlb - xm).sign(false) < 0) {
      T nx = xlb.interval(xm).midP();
      if ((nx - xlb).sign(false) == 0 || (nx - xm).sign(false) == 0) break;
      int sm = value(nx).sign(false);
      if (sm == 0)
        xm = nx;
      else if (sm == sub) {
        xub = nx;
        T nx = xlb.interval(xub);
        return nx.subset(x) ? nx : x;
      } else
        xlb = nx;
    }
    xm = x.midP();
    while ((xm - xub).sign(false) < 0) {
      T nx = xm.interval(xub).midP();
      if ((nx - xm).sign(false) == 0 || (nx - xub).sign(false) == 0) break;
      int sm = value(nx).sign(false);
      if (sm == 0)
        xm = nx;
      else if (sm == -sub) {
        xlb = nx;
        T nx = xlb.interval(xub);
        return nx.subset(x) ? nx : x;
      } else
        xub = nx;
    }
    T nx = xlb.interval(xub);
    return nx.subset(x) ? nx : x;
  }

  PPoly moebius(const T& l, const T& u) const {
    return shift(l).dual().shift(1.0 / (u - l));
  }

  PPoly shift(const T& s) const {
    PPoly p = *this;
    if (s.sign() == 0) return p;
    int d = deg();
    for (int i = 0; i < d; ++i)
      for (int j = d - 1; j >= i; --j) p.a[j] = p.a[j] + s * p.a[j + 1];
    return p;
  }

  void print() const {
    int d = deg();
    cerr << "( ";
    for (int i = 0; i <= d; ++i) cerr << a[i].mid() << " ";
    cerr << ")" << endl;
  }

  vector<T> a;
};

template <class T>
bool zerop(const PPoly<T>& p) {
  return p.zero() || (p.deg() == 0 && zerop(p.lc()));
}
}  // namespace acp

template <class N>
class PPoly2 {
 public:
  PPoly2() : d(-2), degx(-2), degy(-2) {}
  PPoly2(int nt) : a(nt), m(2 * nt), d(-2), degx(-2), degy(-2) {}
  PPoly2(int nt, N* a, int* m) : a(nt), m(2 * nt), d(-2), degx(-2), degy(-2) {
    int k = nt * 2;
    for (int i = 0; i < nt; ++i) this->a[i] = a[i];
    for (int i = 0; i < k; ++i) this->m[i] = m[i];
    setDegree();
    degx = degX();
    degy = degY();
  }
  template <class M>
  PPoly2(const PPoly2<M>& p)
      : a(p.a.size()), m(p.m), d(p.d), degx(p.degx), degy(p.degy) {
    for (int i = 0; i < a.size(); ++i) a[i] = p.a[i];
  }
  void setDegree() {
    d = 0;
    for (int i = 0; i < a.size(); ++i) {
      int di = 0;
      for (int j = 0; j < 2; ++j) di += m[2 * i + j];
      if (d < di) d = di;
    }
    degx = degX();
    degy = degY();
  }

  int degree() const { return d; }
  bool increased() const { return false; /* a[0].increased(); */ }
  bool decreased() const { return true; /* a[0].decreased(); */ }
  N value(N* x) const {
    N xp[degx + 1];
    N yp[degy + 1];

    N p = x[0];
    for (int i = 1; i <= degx; i++) {
      xp[i] = p;
      p = p * x[0];
    }

    p = x[1];
    for (int i = 1; i <= degy; i++) {
      yp[i] = p;
      p = p * x[1];
    }

    N y;
    auto mp = m.begin();
    for (int i = 0; i < a.size(); ++i) {
      N z = a[i];
      if (*mp) z = z * xp[*mp];
      mp++;
      if (*mp) z = z * yp[*mp];
      mp++;
      if (i == 0)
        y = z;
      else
        y = y + z;
    }

    return y;
  }

  N value(const std::initializer_list<N>& x) const {
    N xy[2] = {*(x.begin()), *(x.begin() + 1)};
    return value(xy);
  }

  N value(const PV2<N> p) {
    N xy[2] = {p.x, p.y};
    return value(xy);
  }

  static PPoly2 one() {
    int nt = 1;
    N a[1];
    a[0] = N(1);
    int m[] = {0, 0};
    return PPoly2(nt, a, m);
  }

  static PPoly2 zero() {
    int nt = 1;
    N a[1];
    a[0] = N(0);
    int m[] = {0, 0};
    return PPoly2(nt, a, m);
  }

#ifdef BLEEN
  PPoly1 subX(N x) const {
    vector<Parameter> coef(degy + 1);

    for (int i = 0; i < coef.size(); i++) {
      coef[i] = Parameter(0);
    }

    vector<N> powx;
    powx.push_back(Parameter(1));
    powx.push_back(x);
    for (int i = 2; i <= degx; i++) {
      powx.push_back(powx[powx.size() - 1] * x);
    }

    for (int i = 0; i < a.size(); i++) {
      int dx = m[2 * i];
      int dy = m[2 * i + 1];
      coef[dy] = coef[dy] + (a[i] * powx[dx]);
    }

    return PPoly1(degx, &coef[0]);
  }

  PPoly1 subY(N y) const {
    vector<Parameter> coef(degx + 1);

    for (int i = 0; i < coef.size(); i++) {
      coef[i] = Parameter(0);
    }

    vector<N> powy;
    powy.push_back(Parameter(1));
    powy.push_back(y);
    for (int i = 2; i <= degy; i++) {
      powy.push_back(powy[powy.size() - 1] * y);
    }

    for (int i = 0; i < a.size(); i++) {
      int dx = m[2 * i];
      int dy = m[2 * i + 1];
      coef[dx] = coef[dx] + (a[i] * powy[dy]);
    }

    return PPoly1(degy, &coef[0]);
  }
#endif

  PPoly2 derX() const {
    vector<N> coef;
    vector<int> pow;

    for (int i = 0; i < a.size(); i++) {
      if (m[2 * i] < 1) continue;
      coef.push_back(a[i] * m[2 * i]);
      pow.push_back(m[2 * i] - 1);
      pow.push_back(m[2 * i + 1]);
    }

    return PPoly2(coef.size(), &coef[0], &pow[0]);
  }

  PPoly2 derY() const {
    vector<N> coef;
    vector<int> pow;

    for (int i = 0; i < a.size(); i++) {
      if (m[2 * i + 1] < 1) continue;
      coef.push_back(a[i] * m[2 * i + 1]);
      pow.push_back(m[2 * i]);
      pow.push_back(m[2 * i + 1] - 1);
    }

    return PPoly2(coef.size(), &coef[0], &pow[0]);
  }

  int degX() const {
    int d = 0;

    for (int i = 0; i < a.size(); i++) {
      int x = m[2 * i];
      d = std::max(d, x);
    }

    return d;
  }

  int degY() const {
    int d = 0;

    for (int i = 0; i < a.size(); i++) {
      int x = m[2 * i + 1];
      d = std::max(d, x);
    }

    return d;
  }

  void print_input() const {
    for (int i = 0; i < a.size(); i++) {
      printf("%.16f %d %d%s", a[i].mid(), m[2 * i], m[2 * i + 1],
             (i == a.size() - 1 ? "\n" : " "));
    }
  }

  void print() const {
    for (int i = 0; i < a.size(); i++) {
      printf("%.16f*x^(%d)*y^(%d) %s", a[i].mid(), m[2 * i], m[2 * i + 1],
             (i == a.size() - 1 ? "\n" : "+ "));
    }
  }

  PV2<N> polish(const PV2<N>& r) const {
    PV2<N> ret = r;
    return ret;
  }

  PPoly2 operator*(const PPoly2& b) const {
    int cnt = 0;
    // N ca[nt*b.nt];
    vector<N> ca(a.size() * b.a.size());
    int cm[2 * a.size() * b.a.size()];
    int exps[2];
    for (int i = 0; i < a.size(); i++) {
      auto mp = m.begin() + i * 2;
      for (int j = 0; j < b.a.size(); j++) {
        N coef = a[i] * b.a[j];
        auto bmp = b.m.begin() + j * 2;
        for (int k = 0; k < 2; k++) exps[k] = mp[k] + bmp[k];
        int l;
        for (l = 0; l < cnt; l++) {
          int* cmp = cm + l * 2;
          int h;
          for (h = 0; h < 2; h++)
            if (cmp[h] != exps[h]) break;
          if (h == 2) {
            ca[l] = ca[l] + coef;
            break;
          }
        }
        if (l == cnt) {
          int* cmp = cm + l * 2;
          for (int h = 0; h < 2; h++) cmp[h] = exps[h];
          ca[l] = coef;
          cnt++;
        }
      }
    }
    // return PPoly2(cnt, ca, cm);
    return PPoly2(cnt, &ca[0], cm);
  }

  PPoly2 operator*(const N& p) const {
    vector<N> ca(a.size());
    int cm[2 * a.size()];

    for (int i = 0; i < a.size(); i++) {
      ca[i] = a[i] * p;
      cm[2 * i] = m[2 * i];
      cm[2 * i + 1] = m[2 * i + 1];
    }

    return PPoly2(a.size(), &ca[0], cm);
  }

  PPoly2 operator*(const double p) const {
    vector<N> ca(a.size());
    int cm[2 * a.size()];

    for (int i = 0; i < a.size(); i++) {
      ca[i] = a[i] * p;
      cm[2 * i] = m[2 * i];
      cm[2 * i + 1] = m[2 * i + 1];
    }

    return PPoly2(a.size(), &ca[0], cm);
  }

  PPoly2 plusCtimes(double c, const PPoly2& b) const {
    int pnt = a.size();
    vector<N> pa(a.size() + b.a.size());
    int pm[2 * (a.size() + b.a.size())];

    for (int i = 0; i < a.size(); i++) pa[i] = a[i];
    int ntnv = a.size() * 2;
    for (int i = 0; i < ntnv; i++) pm[i] = m[i];

    auto bmp = b.m.begin();
    for (int i = 0; i < b.a.size(); i++) {
      int j;
      for (j = 0; j < a.size(); j++) {
        auto mp = m.begin() + j * 2;
        int k;
        for (k = 0; k < 2; k++)
          if (mp[k] != bmp[k]) break;
        if (k == 2) {
          pa[j] = pa[j] + c * b.a[i];
          break;
        }
      }
      if (j == a.size()) {
        pa[pnt] = c * b.a[i];
        int* pmp = pm + pnt * 2;
        for (int k = 0; k < 2; k++) pmp[k] = bmp[k];
        pnt++;
      }
      bmp += 2;
    }

    return PPoly2(pnt, &pa[0], pm);
  }

  PPoly2 operator+(const PPoly2& b) const { return plusCtimes(1.0, b); }
  PPoly2 operator-(const PPoly2& b) const { return plusCtimes(-1.0, b); }

  int size() { return a.size(); }
  const N& operator[](int i) const { return a[i]; }
  N& operator[](int i) { return a[i]; }

  int d;
  vector<int> m;
  vector<N> a;
  int degx, degy;
};

namespace PPoly3Inner {
class Index3 {
 public:
  int i[3];
  Index3(int i0, int i1, int i2) {
    i[0] = i0;
    i[1] = i1;
    i[2] = i2;
  }
};

class Hasher {
 public:
  size_t operator()(const Index3& a) const {
    return a.i[0] * 779230947 + a.i[1] * 247091631 + a.i[2] * 1194289623;
  }
};

class Equals {
 public:
  bool operator()(const Index3& a, const Index3& b) const {
    return a.i[0] == b.i[0] && a.i[1] == b.i[1] && a.i[2] == b.i[2];
  }
};
}  // namespace PPoly3Inner

template <class N>
class PPoly3 {
  typedef PPoly3Inner::Index3 Index3;
  typedef PPoly3Inner::Hasher Hasher;
  typedef PPoly3Inner::Equals Equals;

 public:
  PPoly3() : xd(-1), yd(-1), zd(-1), d(-1) {}

  void add(int i, int j, int k, const N& c) {
    auto it = ia.find(Index3(i, j, k));
    if (it == ia.end()) {
      ia[Index3(i, j, k)] = a.size();
      a.push_back(c);
    } else
      a[it->second] = a[it->second] + c;

    xd = i > xd ? i : xd;
    yd = j > yd ? j : yd;
    zd = k > zd ? k : zd;

    d = (i + j + k) > d ? (i + j + k) : d;
  }

  template <class M>
  PPoly3(const PPoly3<M>& p)
      : xd(p.xd), yd(p.yd), zd(p.zd), d(p.d), ia(p.ia), a(p.a.size()) {
    for (int i = 0; i < a.size(); i++) a[i] = p.a[i];
  }

  N value(const N& x, const N& y, const N& z) const {
    const N* xyz[3] = {&x, &y, &z};

    vector<N> powers[3];
    for (int i = 0; i < 3; i++) powers[i].push_back(N::constant(1));

    N sum = N::constant(0);

    for (auto pair = ia.begin(); pair != ia.end(); pair++) {
      Index3 ind = pair->first;
      N term = a[pair->second];

      for (int i = 0; i < 3; i++) {
        while (powers[i].size() <= ind.i[i])
          powers[i].push_back(*xyz[i] * powers[i].back());
        term = term * powers[i][ind.i[i]];
      }

      sum = sum + term;
    }

    return sum;
  }

  PV3<N> gradient(const N& x, const N& y, const N& z) const {
    const N* xyz[3] = {&x, &y, &z};

    vector<N> powers[3];
    for (int i = 0; i < 3; i++) powers[i].push_back(N::constant(1));

    PV3<N> sum = PV3<N>(N(0), N(0), N(0));

    for (auto pair = ia.begin(); pair != ia.end(); pair++) {
      Index3 ind = pair->first;

      for (int i = 0; i < 3; i++)
        while (powers[i].size() <= ind.i[i])
          powers[i].push_back(*xyz[i] * powers[i].back());

      for (int i = 0; i < 3; i++) {
        if (ind.i[i] == 0) continue;
        N term = a[pair->second] * ind.i[i];
        for (int j = 0; j < 3; j++)
          term = term * powers[j /*i*/][ind.i[j] - (i == j)];
        sum[i] = sum[i] + term;
      }
    }

    return sum;
  }

  PPoly3 operator+(const PPoly3& that) const {
    PPoly3 q;

    for (auto it = ia.begin(); it != ia.end(); it++) {
      Index3 i = it->first;
      int ind = it->second;
      q.add(i.i[0], i.i[1], i.i[2], a[ind]);
    }

    for (auto it = that.ia.begin(); it != that.ia.end(); it++) {
      Index3 i = it->first;
      int ind = it->second;
      q.add(i.i[0], i.i[1], i.i[2], that.a[ind]);
    }

    return q;
  }

  PPoly3 operator-(const PPoly3& that) const {
    PPoly3 q;

    for (auto it = ia.begin(); it != ia.end(); it++) {
      Index3 i = it->first;
      int ind = it->second;
      q.add(i.i[0], i.i[1], i.i[2], a[ind]);
    }

    for (auto it = that.ia.begin(); it != that.ia.end(); it++) {
      Index3 i = it->first;
      int ind = it->second;
      q.add(i.i[0], i.i[1], i.i[2], -that.a[ind]);
    }

    return q;
  }

  PPoly3 operator*(const PPoly3& that) const {
    PPoly3 q;

    for (auto it = ia.begin(); it != ia.end(); it++) {
      Index3 ai = it->first;
      int aai = it->second;

      for (auto it = that.ia.begin(); it != that.ia.end(); it++) {
        Index3 bi = it->first;
        int bai = it->second;

        int i = ai.i[0] + bi.i[0];
        int j = ai.i[1] + bi.i[1];
        int k = ai.i[2] + bi.i[2];

        q.add(i, j, k, a[aai] * that.a[bai]);
      }
    }

    return q;
  }

  // xyz = 0, 1, 2:  substituting for yz, xz, xy.
  PPoly<N> substitute2(int ixyz, const N& x, const N& y, const N& z) const {
    const N* xyz[3] = {&x, &y, &z};

    vector<N> powers[3];
    for (int i = 0; i < 3; i++) powers[i].push_back(N::constant(1));

    PPoly<N> ppoly;

    for (auto pair = ia.begin(); pair != ia.end(); pair++) {
      Index3 ind = pair->first;
      N term = a[pair->second];

      for (int i = 0; i < 3; i++)
        if (i != ixyz) {
          while (powers[i].size() <= ind.i[i])
            powers[i].push_back(*xyz[i] * powers[i].back());
          term = term * powers[i][ind.i[i]];
        }

      ppoly.add(term, ind.i[ixyz]);
    }

    return ppoly;
  }

  // t is a flattened 4x4 matrix in column major form
  PPoly3 transform(double t[16]) const {
    PPoly3 x;
    PPoly3 y;
    PPoly3 z;
    PPoly3 w;

    PPoly3 unit;

    // First column
    x.add(1, 0, 0, N::constant(t[0]));
    y.add(1, 0, 0, N::constant(t[1]));
    z.add(1, 0, 0, N::constant(t[2]));
    w.add(1, 0, 0, N::constant(t[3]));

    // Second column
    x.add(0, 1, 0, N::constant(t[4]));
    y.add(0, 1, 0, N::constant(t[5]));
    z.add(0, 1, 0, N::constant(t[6]));
    w.add(0, 1, 0, N::constant(t[7]));

    // Third column
    x.add(0, 0, 1, N::constant(t[8]));
    y.add(0, 0, 1, N::constant(t[9]));
    z.add(0, 0, 1, N::constant(t[10]));
    w.add(0, 0, 1, N::constant(t[11]));

    // Fourth column
    x.add(0, 0, 0, N::constant(t[12]));
    y.add(0, 0, 0, N::constant(t[13]));
    z.add(0, 0, 0, N::constant(t[14]));
    w.add(0, 0, 0, N::constant(t[15]));

    unit.add(0, 0, 0, N::constant(1));

    std::vector<PPoly3> xpow(xd + 1);
    std::vector<PPoly3> ypow(yd + 1);
    std::vector<PPoly3> zpow(zd + 1);
    std::vector<PPoly3> wpow(d + 1);

    xpow[0] = ypow[0] = zpow[0] = wpow[0] = unit;

    for (int i = 1; i <= xd; i++) xpow[i] = x * xpow[i - 1];
    for (int i = 1; i <= yd; i++) ypow[i] = y * ypow[i - 1];
    for (int i = 1; i <= zd; i++) zpow[i] = z * zpow[i - 1];
    for (int i = 1; i <= d; i++) wpow[i] = w * wpow[i - 1];

    PPoly3 sum;

    for (auto it = ia.begin(); it != ia.end(); it++) {
      Index3 ai = it->first;
      int aai = it->second;

      int i = ai.i[0];
      int j = ai.i[1];
      int k = ai.i[2];

      PPoly3 coef;
      coef.add(0, 0, 0, a[aai]);

      sum = sum + (coef * xpow[i] * ypow[j] * zpow[k] * wpow[d - (i + j + k)]);
    }

    return sum;
  }

  unordered_map<Index3, int, Hasher, Equals> ia;
  vector<N> a;
  int size() const { return a.size(); }
  const N& operator[](int i) const { return a[i]; }
  N& operator[](int i) { return a[i]; }
  int xd, yd, zd, d;
};

#endif
