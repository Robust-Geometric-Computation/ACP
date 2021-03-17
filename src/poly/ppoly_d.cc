#include "acp/poly/ppoly_d.h"
using namespace std;
using namespace acp;

namespace acp {

PPoly1D::PPoly1D(int d, double* a) : d(d) {
  this->a = new double[d + 1];
  if (a)
    for (int i = 0; i <= d; ++i) this->a[i] = a[i];
}

PPoly1D::PPoly1D(const PPoly1D& f) : d(f.d) {
  a = new double[d + 1];
  for (int i = 0; i <= d; ++i) a[i] = f.a[i];
}

#ifdef BLEEN
PPoly1D::PPoly1D(const PPoly& f) : d(f.size() - 1) {
  a = new double[d + 1];
  for (int i = 0; i <= d; ++i) a[i] = f.a[i].mid();
}
#endif

PPoly1D& PPoly1D::operator=(const PPoly1D& f) {
  delete[] a;
  d = f.d;
  a = new double[d + 1];
  for (int i = 0; i <= d; ++i) a[i] = f.a[i];
  return *this;
}

double PPoly1D::value(const double& x) const {
  double y = a[d];
  for (int i = d - 1; i >= 0; --i) y = x * y + a[i];
  return y;
}

double PPoly1D::der(const double& x) const {
  double y = d * a[d];
  for (int i = d - 1; i > 0; --i) y = x * y + i * a[i];
  return y;
}

PPoly1D PPoly1D::der() const {
  PPoly1D g(d - 1);
  for (int i = 1; i <= d; ++i) g.a[i - 1] = i * a[i];
  return g;
}

PPoly2D::PPoly2D(int nt, double* a, int* m) : nt(nt) {
  int k = nt * 2;
  this->a = new double[nt];
  this->m = new int[k];
  for (int i = 0; i < nt; ++i) this->a[i] = a[i];
  for (int i = 0; i < k; ++i) this->m[i] = m[i];
  setDegree();
}

PPoly2D::PPoly2D(const PPoly2D& p) : nt(p.size()) {
  if (p.a) {
    int k = nt * 2;
    a = new double[nt];
    m = new int[k];
    for (int i = 0; i < nt; ++i) a[i] = p.a[i];
    for (int i = 0; i < k; ++i) m[i] = p.m[i];
  } else {
    a = 0;
    m = 0;
  }
  d = p.d;
}

void PPoly2D::setDegree() {
  d = 0;
  for (int i = 0; i < nt; ++i) {
    int di = 0;
    for (int j = 0; j < 2; ++j) di += m[2 * i + j];
    if (d < di) d = di;
  }
}

PPoly2D& PPoly2D::operator=(const PPoly2D& p) {
  nt = p.nt;
  delete[] a;
  delete[] m;
  if (p.a) {
    int k = nt * 2;
    a = new double[nt];
    m = new int[k];
    for (int i = 0; i < nt; ++i) a[i] = p.a[i];
    for (int i = 0; i < k; ++i) m[i] = p.m[i];
  } else {
    a = 0;
    m = 0;
  }
  d = p.d;
  return *this;
}

double PPoly2D::value(double* x) const {
  // ??? vjm sentinal problem with y+z
  double y;
  int* mp = m;
  for (int i = 0; i < nt; ++i) {
    double z = a[i];
    for (int j = 0; j < 2; ++j, ++mp)
      for (int k = 0; k < *mp; ++k) z = z * x[j];
    if (i == 0)
      y = z;
    else
      y = y + z;
  }
  return y;
}

PPoly2D PPoly2D::operator*(const PPoly2D& b) const {
  int cnt = 0;
  // double ca[nt*b.nt];
  vector<double> ca(nt * b.nt);
  int cm[2 * nt * b.nt];
  int exps[2];
  for (int i = 0; i < nt; i++) {
    int* mp = m + i * 2;
    for (int j = 0; j < b.nt; j++) {
      double coef = a[i] * b.a[j];
      int* bmp = b.m + j * 2;
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
  // return PPoly2D(cnt, ca, cm);
  return PPoly2D(cnt, &ca[0], cm);
}

PPoly2D PPoly2D::plusCtimes(double c, const PPoly2D& b) const {
  int pnt = nt;
  vector<double> pa(nt + b.nt);
  int pm[2 * (nt + b.nt)];

  for (int i = 0; i < nt; i++) pa[i] = a[i];
  int ntnv = nt * 2;
  for (int i = 0; i < ntnv; i++) pm[i] = m[i];

  int* bmp = b.m;
  for (int i = 0; i < b.nt; i++) {
    int j;
    for (j = 0; j < nt; j++) {
      int* mp = m + j * 2;
      int k;
      for (k = 0; k < 2; k++)
        if (mp[k] != bmp[k]) break;
      if (k == 2) {
        pa[j] = pa[j] + c * b.a[i];
        break;
      }
    }
    if (j == nt) {
      pa[pnt] = c * b.a[i];
      int* pmp = pm + pnt * 2;
      for (int k = 0; k < 2; k++) pmp[k] = bmp[k];
      pnt++;
    }
    bmp += 2;
  }

  return PPoly2D(pnt, &pa[0], pm);
}

PPoly2D PPoly2D::derX() const {
  vector<double> coef;
  vector<int> pow;

  for (int i = 0; i < nt; i++) {
    if (m[2 * i] < 1) continue;
    coef.push_back(a[i] * m[2 * i]);
    pow.push_back(m[2 * i] - 1);
    pow.push_back(m[2 * i + 1]);
  }

  return PPoly2D(coef.size(), &coef[0], &pow[0]);
}

PPoly2D PPoly2D::derY() const {
  vector<double> coef;
  vector<int> pow;

  for (int i = 0; i < nt; i++) {
    if (m[2 * i + 1] < 1) continue;
    coef.push_back(a[i] * m[2 * i + 1]);
    pow.push_back(m[2 * i]);
    pow.push_back(m[2 * i + 1] - 1);
  }

  return PPoly2D(coef.size(), &coef[0], &pow[0]);
}

// PPoly2D PPoly2D::derX() const {
//
//  PPoly2D g(nt);
//
//  for(int i = 0; i < nt; i++) {
//    g.a[i] = a[i] * m[2*i];
//    g.m[2*i]   = max(m[2*i] - 1, 0);
//    g.m[2*i+1] = m[2*i + 1];
//  }
//
//  g.setDegree();
//
//  return g;
//
//}
//
// PPoly2D PPoly2D::derY() const {
//
//  PPoly2D g(nt);
//
//  for(int i = 0; i < nt; i++) {
//    g.a[i] = a[i] * m[2*i + 1];
//    g.m[2*i]   = m[2*i];
//    g.m[2*i+1] = max(m[2*i + 1] - 1, 0);
//  }
//
//  g.setDegree();
//
//  return g;
//
//}

//-----------------------------------------------------------

template <unsigned int N>
PPolyD<N>::PPolyD(int nt, double* a, int* m) : nt(nt) {
  int k = nt * N;
  this->a = new double[nt];
  this->m = new int[k];
  for (int i = 0; i < nt; ++i) this->a[i] = a[i];
  for (int i = 0; i < k; ++i) this->m[i] = m[i];
  setDegree();
}

template <unsigned int N>
PPolyD<N>::PPolyD(const PPolyD<N>& p) : nt(p.nt) {
  if (p.a) {
    int k = nt * N;
    a = new double[nt];
    m = new int[k];
    for (int i = 0; i < nt; ++i) a[i] = p.a[i];
    for (int i = 0; i < k; ++i) m[i] = p.m[i];
  } else {
    a = 0;
    m = 0;
  }
  d = p.d;
}

template <unsigned int N>
void PPolyD<N>::setDegree() {
  d = 0;
  for (int i = 0; i < nt; ++i) {
    int di = 0;
    for (int j = 0; j < N; ++j) di += m[N * i + j];
    if (d < di) d = di;
  }
}

template <unsigned int N>
PPolyD<N>& PPolyD<N>::operator=(const PPolyD<N>& p) {
  nt = p.nt;
  delete[] a;
  delete[] m;
  if (p.a) {
    int k = nt * N;
    a = new double[nt];
    m = new int[k];
    for (int i = 0; i < nt; ++i) a[i] = p.a[i];
    for (int i = 0; i < k; ++i) m[i] = p.m[i];
  } else {
    a = 0;
    m = 0;
  }
  d = p.d;
  return *this;
}

template <unsigned int N>
double PPolyD<N>::value(double* x) const {
  // ??? vjm sentinal problem with y+z
  double y;
  int* mp = m;
  for (int i = 0; i < nt; ++i) {
    double z = a[i];
    for (int j = 0; j < N; ++j, ++mp)
      for (int k = 0; k < *mp; ++k) z = z * x[j];
    if (i == 0)
      y = z;
    else
      y = y + z;
  }
  return y;
}

//???
template <unsigned int N>
PPolyD<N> PPolyD<N>::operator*(const PPolyD<N>& b) const {
  int cnt = 0;
  // double ca[nt*b.nt];
  vector<double> ca(nt * b.nt);
  int cm[N * nt * b.nt];
  int exps[N];
  for (int i = 0; i < nt; i++) {
    int* mp = m + i * N;
    for (int j = 0; j < b.nt; j++) {
      double coef = a[i] * b.a[j];
      int* bmp = b.m + j * N;
      for (int k = 0; k < N; k++) exps[k] = mp[k] + bmp[k];
      int l;
      for (l = 0; l < cnt; l++) {
        int* cmp = cm + l * N;
        int h;
        for (h = 0; h < N; h++)
          if (cmp[h] != exps[h]) break;
        if (h == N) {
          ca[l] = ca[l] + coef;
          break;
        }
      }
      if (l == cnt) {
        int* cmp = cm + l * N;
        for (int h = 0; h < N; h++) cmp[h] = exps[h];
        ca[l] = coef;
        cnt++;
      }
    }
  }
  // return PPoly2D(cnt, ca, cm);
  return PPoly2D(cnt, &ca[0], cm);
}

//???
template <unsigned int N>
PPolyD<N> PPolyD<N>::plusCtimes(double c, const PPolyD<N>& b) const {
  int pnt = nt;
  vector<double> pa(nt + b.nt);
  int pm[N * (nt + b.nt)];

  for (int i = 0; i < nt; i++) pa[i] = a[i];
  int ntnv = nt * N;
  for (int i = 0; i < ntnv; i++) pm[i] = m[i];

  int* bmp = b.m;
  for (int i = 0; i < b.nt; i++) {
    int j;
    for (j = 0; j < nt; j++) {
      int* mp = m + j * N;
      int k;
      for (k = 0; k < N; k++)
        if (mp[k] != bmp[k]) break;
      if (k == N) {
        pa[j] = pa[j] + c * b.a[i];
        break;
      }
    }
    if (j == nt) {
      pa[pnt] = c * b.a[i];
      int* pmp = pm + pnt * N;
      for (int k = 0; k < N; k++) pmp[k] = bmp[k];
      pnt++;
    }
    bmp += N;
  }

  return PPoly2D(pnt, &pa[0], pm);
}

template <unsigned int N>
PPolyD<N> PPolyD<N>::der(int j) const {
  assert(j < N && j >= 0);

  PPolyD<N> g(nt);

  for (int i = 0; i < nt; i++) {
    g.a[i] = a[i] * m[N * i + j];

    g.m[N * i + j] = max(m[N * i + j] - 1, 0);

    for (int k = 0; k < N; k++) {
      if (k == j) continue;
      g.m[N * i + k] = m[N * i + k];
    }
  }

  g.setDegree();

  return g;
}

}  // namespace acp
