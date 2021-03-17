#ifndef PPOLYD
#define PPOLYD

#include <vector>
#include "acp/core/acp.h"
#include "acp/poly/poly2.h"

namespace acp {

class PPoly1D {
 public:
  PPoly1D(int d, double* a = 0);
  PPoly1D(const PPoly1D& f);
  template <class N>
  PPoly1D(const Poly<N>& f);
  ~PPoly1D() { delete[] a; }
  PPoly1D& operator=(const PPoly1D& f);

  double value(const double& x) const;
  double der(const double& x) const;
  PPoly1D der() const;

  PPoly1D operator*(const PPoly1D& g) const;
  PPoly1D plusCtimes(double c, const PPoly1D& b) const;
  PPoly1D operator+(const PPoly1D& b) const { return plusCtimes(1.0, b); }
  PPoly1D operator-(const PPoly1D& b) const { return plusCtimes(-1.0, b); }

  void print();

  int size() const { return d + 1; }
  double& operator[](int i) { return a[i]; }

  int d;
  double* a;
};

class PPoly2D {
 public:
  PPoly2D() : nt(0), a(0), m(0) {}
  PPoly2D(int nt) : nt(nt), a(new double[nt]), m(new int[2 * nt]) {}
  PPoly2D(int nt, double* a, int* m);
  PPoly2D(const PPoly2D& p);

  template <class N>
  PPoly2D(const Poly2<N>& p) : nt(p.size()) {
    int k = nt * 2;
    a = new double[nt];
    m = new int[k];
    for (int i = 0; i < nt; ++i) a[i] = p.a[i].mid();
    for (int i = 0; i < k; ++i) m[i] = p.m[i];

    d = p.d;
  }

  void setDegree();
  ~PPoly2D() {
    delete[] a;
    delete[] m;
  }
  PPoly2D& operator=(const PPoly2D& p);
  int degree() const { return d; }
  double value(double* x) const;

  PPoly2D derX() const;
  PPoly2D derY() const;

  PPoly2D operator*(const PPoly2D& g) const;
  PPoly2D plusCtimes(double c, const PPoly2D& b) const;
  PPoly2D operator+(const PPoly2D& b) const { return plusCtimes(1.0, b); }
  PPoly2D operator-(const PPoly2D& b) const { return plusCtimes(-1.0, b); }

  int size() const { return nt; }
  double& operator[](int i) { return a[i]; }

  int nt, d;
  double* a;
  int* m;
};

template <unsigned int N>
class PPolyD {
 public:
  PPolyD() : nt(0), a(0), m(0) {}
  PPolyD(int nt) : nt(nt), a(new double[nt]), m(new int[N * nt]) {}
  PPolyD(int nt, double* a, int* m);
  PPolyD(const PPolyD<N>& p);

  void setDegree();
  ~PPolyD() {
    delete[] a;
    delete[] m;
  }
  PPolyD<N>& operator=(const PPolyD<N>& p);
  int degree() const { return d; }
  double value(double* x) const;

  PPolyD<N> der(int i) const;

  PPolyD<N> operator*(const PPolyD<N>& g) const;
  PPolyD<N> plusCtimes(double c, const PPolyD<N>& b) const;
  PPolyD<N> operator+(const PPolyD<N>& b) const { return plusCtimes(1.0, b); }
  PPolyD<N> operator-(const PPolyD<N>& b) const { return plusCtimes(-1.0, b); }

  int size() { return nt; }
  double& operator[](int i) { return a[i]; }

  int nt, d;
  double* a;
  int* m;
};

}  // namespace acp
#endif
