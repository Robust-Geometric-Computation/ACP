#ifndef POLY2D_H
#define POLY2D_H

#include <vector>
#include "acp/poly/poly.h"

using namespace std;
using namespace acp;

template <class N>
class Poly2D {
  vector<N> a;
  int nx, ny;
  int ind(int ix, int iy) const { return ny * ix + iy; }
  void add(N& ai, const N& c) { ai = ai + c; }

 public:
  int size() const { return a.size(); }
  const N& operator[](int i) const { return a[i]; }
  N& operator[](int i) { return a[i]; }

  Poly2D() : nx(0), ny(0), a(0) {}
  Poly2D(int degx, int degy)
      : nx(degx + 1), ny(degy + 1), a((degx + 1) * (degy + 1)) {
    for (typename vector<N>::iterator it = a.begin(); it != a.end(); ++it)
      *it = N::constant(0);
  }
  template <class M>
  Poly2D(const Poly2D<M>& p) : nx(p.degX() + 1), ny(p.degY() + 1), a(p.size()) {
    for (int i = 0; i < a.size(); ++i) a[i] = p[i];
  }

  int degX() const { return nx - 1; }
  int degY() const { return ny - 1; }

  N get(int ix, int iy) const { return a[ind(ix, iy)]; }
  N set(int ix, int iy, const N& c) { return a[ind(ix, iy)] = c; }
  void add(int ix, int iy, const N& c) { return add(a[ind(ix, iy)], c); }

  N value(const N& x, const N& y) const {
    N sx = N::constant(0);
    for (int ix = nx; --ix >= 0;) {
      N sy = N::constant(0);
      for (int iy = ny; --iy >= 0;) sy = y * sy + get(ix, iy);
      sx = x * sx + sy;
    }
    return sx;
  }

  N value(const PV2<N>& xy) const { return value(xy.x, xy.y); }

  Poly<N> substitute1(int ixy, const N& v) const {
    int n[2] = {nx, ny};
    int i[2] = {0, 0};
    int iy = ixy, ix = (iy + 1) % 2;
    vector<N> a2;
    for (i[ix] = 0; i[ix] < n[ix]; ++i[ix]) {
      i[iy] = n[iy] - 1;
      N c = get(i[0], i[1]);
      for (i[iy] = n[iy] - 1; --i[iy] >= 0;) c = c * v + get(i[0], i[1]);
      a2.push_back(c);
    }

    Poly<N> p1 = Poly<N>(a2);

    // debug
    N xy[2];
    xy[ix] = N::constant(2.781828);
    xy[iy] = v;
    N p1v = p1.value(xy[ix]);
    N p2v = value(xy[0], xy[1]);
    assert((p1v - p2v).sign(false) == 0);

    return p1;
  }

  Poly2D<N> der(int xy) {
    int n[2] = {nx, ny};
    int iy = xy, ix = (iy + 1) % 2;
    --n[iy];
    Poly2D<N> poly(n[0] - 1, n[1] - 1);
    ++n[iy];
    int i[2] = {0, 0};
    for (i[ix] = 0; i[ix] < n[ix]; ++i[ix])
      for (i[iy] = 1; i[iy] < n[iy]; ++i[iy]) {
        N c = get(i[0], i[1]) * i[iy];
        --i[iy];
        poly.set(i[0], i[1], c);
        ++i[iy];
      }
    return poly;
  }

  N der(int xy, const PV2<N>& v) {
    int n[2] = {nx, ny};
    int iy = xy, ix = (iy + 1) % 2;
    int i[2] = {0, 0};
    N sx = N::constant(0);
    for (i[ix] = n[ix]; --i[ix] >= 0;) {
      N sy = N::constant(0);
      for (i[iy] = n[iy]; --i[iy] >= 1;)
        sy = v[iy] * sy + get(i[0], i[1]) * i[iy];
      sx = v[ix] * sx + sy;
    }

    Poly2D<N> df = der(xy);
    N d = df.value(v);
    assert((d - sx).sign(false) == 0);

    return sx;
  }

  PV2<N> gradient(const PV2<N>& v) { return PV2<N>(der(0, v), der(1, v)); }

#ifdef BLEEN
  Poly2D<N> operator+(const Poly2D<N>& b) const {
    Poly2D<N> p(std::max(degX(), b.degX()), std::max(degY(), b.degY()),
                std::max(degZ(), b.degZ()));

    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
        for (int iz = 0; iz < nz; ++iz) p.set(ix, iy, iz, get(ix, iy, iz));

    for (int ix = 0; ix < b.nx; ++ix)
      for (int iy = 0; iy < b.ny; ++iy)
        for (int iz = 0; iz < nz; ++iz)
          p.set(ix, iy, iz, p.get(ix, iy, iz) + b.get(ix, iy, iz));

    return p;
  }

  Poly2D<N> operator-(const Poly2D<N>& b) const {
    Poly2D<N> p(std::max(degX(), b.degX()), std::max(degY(), b.degY()),
                std::max(degZ(), b.degZ()));

    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
        for (int iz = 0; iz < nz; ++iz) p.set(ix, iy, iz, get(ix, iy, iz));

    for (int ix = 0; ix < b.nx; ++ix)
      for (int iy = 0; iy < b.ny; ++iy)
        for (int iz = 0; iz < nz; ++iz)
          p.set(ix, iy, iz, p.get(ix, iy, iz) - b.get(ix, iy, iz));

    return p;
  }

  Poly2D<N> operator*(const Poly2D<N>& b) const {
    Poly2D<N> p(degX() + b.degX(), degY() + b.degY(), degZ() + b.degZ());

    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
        for (int iz = 0; iz < nz; ++iz)
          for (int jx = 0; jx < b.nx; ++jx)
            for (int jy = 0; jy < b.ny; ++jy)
              for (int jz = 0; jz < b.nz; ++jz)
                p.add(ix + jx, iy + jy, iz + jz,
                      get(ix, iy, iz) * b.get(jx, jy, jz));

    return p;
  }
#endif

  void mobius(int xy, const N& lb, const N& ub) {
    int n[2] = {nx, ny}, i[2];
    int iy = xy, ix = (iy + 1) % 2;

    for (int m = 0; m < n[iy] - 1; m++)
      for (i[iy] = n[iy] - 1 - 1; i[iy] >= m; i[iy]--)
        for (i[ix] = 0; i[ix] < n[ix]; i[ix]++) {
          N c = get(i[0], i[1]);
          ++i[iy];
          c = c + lb * get(i[0], i[1]);
          --i[iy];
          set(i[0], i[1], c);
        }

    N s = 1.0 / (ub - lb);
    for (int m = 0; m < n[iy] - 1; m++)
      for (i[iy] = 1; i[iy] <= n[iy] - 1 - m; i[iy]++)
        for (i[ix] = 0; i[ix] < n[ix]; i[ix]++) {
          N c = get(i[0], i[1]);
          --i[iy];
          c = c + s * get(i[0], i[1]);
          ++i[iy];
          set(i[0], i[1], c);
        }
  }

  void mobius(const PV2<Parameter>& box) {
    // PV2<N> b(box);
    for (int xy = 0; xy < 2; ++xy)
      // mobius(xy, b[xy].lbP(), b[xy].ubP());
      mobius(xy, N::constant(box[xy].lb()), N::constant(box[xy].ub()));
  }
};

class Poly2Dval : public Primitive {
  PTR<Object<Poly2D>> f;
  PTR<Object<PV2>> xy;
  DeclareSign {
    Poly2D<N> f = this->f->get<N>();
    PV2<N> xy = this->xy->get<N>();
    return f.value(xy);
  }

 public:
  Poly2Dval(PTR<Object<Poly2D>> f, PTR<Object<PV2>> xy) : f(f), xy(xy) {}
};

vector<PTR<Object<PV2>>> getRoots(PTR<Object<Poly2D>> f, PTR<Object<Poly2D>> g,
                                  PV2<Parameter> box);

class Sub1Poly2D : public Object<Poly> {
  PTR<Object<Poly2D>> p2d;
  int ixy;
  PTR<Object<Scalar>> xy;
  DeclareCalculate(Poly) {
    return p2d->get<N>().substitute1(ixy, xy->get<N>());
  }

 public:
  // ixy = 0 or 1 means substituting for x or y
  Sub1Poly2D(PTR<Object<Poly2D>> p2d, int ixy, PTR<Object<Scalar>> xy)
      : p2d(p2d), ixy(ixy), xy(xy) {}
};

class P1to2 : public Object<PV2> {
  PTR<Object<Scalar>> p1;
  PTR<Object<Scalar>> s;
  int xy;
  DeclareCalculate(PV2) {
    PV2<N> p;
    p[xy] = s->get<N>();
    Scalar<N> p1 = this->p1->get<N>();
    p[(xy + 1) % 2] = p1.x;
    return p;
  }

 public:
  P1to2(PTR<Object<Scalar>> p1, PTR<Object<Scalar>> s, int xy)
      : p1(p1), s(s), xy(xy) {}
};

class PolyVal : public Primitive {
  PTR<Object<Poly>> f;
  PTR<Object<Scalar>> x;
  DeclareSign {
    Poly<N> f = this->f->get<N>();
    N x = this->x->get<N>().x;
    return f.value(x);
  }

 public:
  PolyVal(PTR<Object<Poly>> f, PTR<Object<Scalar>> x) : f(f), x(x) {}
};

static void swap(int& i, int& j) {
  int k = i;
  i = j;
  j = k;
}

template <class N>
bool lessMag(const N& a, const N& b) {
  if (a.sign(false) < 0)
    if (b.sign(false) < 0)
      return (-a - -b).sign(false) < 0;
    else
      return (-a - b).sign(false) < 0;
  else if (b.sign(false) < 0)
    return (a - -b).sign(false) < 0;
  else
    return (a - b).sign(false) < 0;
}

template <class N>
PV2<N> solve(N a[2][3]) {
  int r[2] = {0, 1};
  if (lessMag(a[r[0]][0], a[r[1]][0])) swap(r[0], r[1]);
  N c = a[r[1]][0] / a[r[0]][0];
  a[r[1]][0] = N::constant(0);
  for (int j = 1; j < 3; j++) a[r[1]][j] = a[r[1]][j] - c * a[r[0]][j];
  PV2<N> x;
  x[1] = -a[r[1]][2] / a[r[1]][1];
  x[0] = -(a[r[0]][2] + a[r[0]][1] * x[1]) / a[r[0]][0];
  return x;
}

template <class N>
PV2<N> newton(Poly2D<N> f, Poly2D<N> g, PV2<N> r) {
  PV2<N> m = PV2<N>(r.x.midP(), r.y.midP());
  N a[2][3] = {{f.der(0, r), f.der(1, r), f.value(m)},
               {g.der(0, r), g.der(1, r), g.value(m)}};
  PV2<N> s = m + solve(a);
  return PV2<N>(r.x.intersect(s.x), r.y.intersect(s.y));
}

PV2<Parameter> subdivide(PTR<Object<Poly2D>> f, PTR<Object<Poly2D>> g,
                         PV2<Parameter>& box);

class Root2D : public Object<PV2> {
  PTR<Object<Poly2D>> f, g;
  PV2<Parameter> box;
  DeclareCalculate(PV2) {
#ifdef SET_ALG
    if (std::is_same<N, PParameter>::value) {
      throw SignException(true);
    }
#else
    if (std::is_same<N, PParameter>::value) {
      throw SignException();
    }
#endif
    Poly2D<N> f = this->f->get<N>();
    Poly2D<N> g = this->g->get<N>();

    PV2<N> b = uninitialized() ? PV2<N>(box) : PV2<N>(getCurrentP());

    for (int i = 0; i < 100; i++) {
      PV2<N> bp;
      int count = 0;
      do {
        ++count;
        bp = b;
        try {
          b = newton(f, g, b);
        } catch (const SignException& e) {
          cout << e.what() << endl;
        }
      } while (b.x.subset(bp.x) || b.y.subset(bp.y));

      // instance was 0.0372 by 0.0041 with count = 55 at 212u
      if (curPrecision() == 53u || 
	  (curPrecision() > 53u && count > 1 &&
	   b.x.intervalWidth() < 1e-6 && b.y.intervalWidth() < 1e-6))
	break;

      if (curPrecision() > 53u) {
        unsigned int savePrecision = curPrecision();
        curPrecision() = 53u;
        box = subdivide(this->f, this->g, box);
        b = PV2<N>(box);
        curPrecision() = savePrecision;
      }
    }
    // assert(count > 1);
#ifdef SET_ALG
    b.x.setAlg();
    b.y.setAlg();
#endif
    return b;
  }

 public:
  Root2D(PTR<Object<Poly2D>> f, PTR<Object<Poly2D>> g, PV2<Parameter> box)
      : f(f), g(g), box(box) {}
};

class Descartes2D : public Primitive {
  PTR<Object<Poly2D>> f;
  PV2<Parameter> box;
  DeclareSign {
    Poly2D<N> f = this->f->get<N>();
    f.mobius(box);
    int s = f[0].sign();
    for (int i = 1; i < f.size(); ++i) {
      int si = f[i].sign();
      if (si != 0 && si != s) return N::constant(0);
    }
    return f[0];
  }

 public:
  Descartes2D(PTR<Object<Poly2D>> f, PV2<Parameter> box) : f(f), box(box) {}
};

// functions
vector<PTR<Object<PV2>>> boundaryRoots(PTR<Object<Poly2D>> f,
                                       const PV2<Parameter>& box, int xy,
                                       int lu);

vector<PTR<Object<PV2>>> boundaryRoots(PTR<Object<Poly2D>> f,
                                       const PV2<Parameter>& box);

vector<PTR<Object<PV2>>> getRoots(PTR<Object<Poly2D>> f, PTR<Object<Poly2D>> g,
                                  PV2<Parameter> box);

#endif
