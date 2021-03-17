#ifndef POLY3D_H
#define POLY3D_H

#include <vector>
#include "acp/encasement2d/encasement2d.h"
#include "acp/poly/poly2.h"
#include "acp/poly/poly2d.h"

using namespace std;
using namespace acp;

template <class N>
class Poly3D {
  int nx, ny, nz;
  int ind(int ix, int iy, int iz) const { return nz * (ny * ix + iy) + iz; }
  void add(N& ai, const N& c) { ai = ai + c; }

 public:
  vector<N> a;
  int size() const { return a.size(); }
  const N& operator[](int i) const { return a[i]; }
  N& operator[](int i) { return a[i]; }

  Poly3D() : nx(0), ny(0), nz(0), a(0) {}
  Poly3D(int degx, int degy, int degz)
      : nx(degx + 1),
        ny(degy + 1),
        nz(degz + 1),
        a((degx + 1) * (degy + 1) * (degz + 1)) {
    for (typename vector<N>::iterator it = a.begin(); it != a.end(); ++it)
      *it = N::constant(0);
  }
  template <class M>
  Poly3D(const Poly3D<M>& p)
      : nx(p.degX() + 1), ny(p.degY() + 1), nz(p.degZ() + 1), a(p.size()) {
    for (int i = 0; i < a.size(); ++i) a[i] = p[i];
  }

  int degX() const { return nx - 1; }
  int degY() const { return ny - 1; }
  int degZ() const { return nz - 1; }

  N get(int ix, int iy, int iz) const { return a[ind(ix, iy, iz)]; }
  N set(int ix, int iy, int iz, const N& c) { return a[ind(ix, iy, iz)] = c; }
  void add(int ix, int iy, int iz, const N& c) {
    return add(a[ind(ix, iy, iz)], c);
  }

  N value(const N& x, const N& y, const N& z) const {
    N sx = N::constant(0);
    for (int ix = nx; --ix >= 0;) {
      N sy = N::constant(0);
      for (int iy = ny; --iy >= 0;) {
        N sz = N::constant(0);
        for (int iz = nz; --iz >= 0;) sz = z * sz + get(ix, iy, iz);
        sy = y * sy + sz;
      }
      sx = x * sx + sy;
    }
    return sx;
  }

  N value(const PV3<N>& xyz) const { return value(xyz.x, xyz.y, xyz.z); }

  Poly<N> substitute2(int ix, const N& vx, int iy, const N& vy) {
    int iz = 3 - ix - iy;
    int n[3] = {nx, ny, nz};
    int i[3] = {0, 0, 0};

    vector<N> a1(n[iz]);

    for (i[iz] = 0; i[iz] < n[iz]; ++i[iz]) {
      i[ix] = n[ix];
      N cx(0);
      ;
      while (--i[ix] >= 0) {
        i[iy] = n[iy] - 1;
        N cy = get(i[0], i[1], i[2]);
        while (--i[iy] >= 0) cy = cy * vy + get(i[0], i[1], i[2]);

        cx = cx * vx + cy;
      }

      a1[i[iz]] = cx;
    }

    Poly<N> p1(a1);

    // debug
    N xyz[3];
    xyz[ix] = vx;
    xyz[iy] = vy;
    xyz[iz] = N::constant(1 /*-3.1415926*/);
    N p1v = p1.value(xyz[iz]);
    N p3v = value(xyz[0], xyz[1], xyz[2]);
    assert((p1v - p3v).sign(false) == 0);

    return p1;
  }

  Poly2<N> substitute1(int ixyz, const N& v) const {
    int n[3] = {nx, ny, nz};
    int i[3] = {0, 0, 0};

    int iz = ixyz, ix = (iz + 1) % 3, iy = (ix + 1) % 3;

    vector<N> a2;
    vector<int> m2;
    for (i[ix] = 0; i[ix] < n[ix]; ++i[ix])
      for (i[iy] = 0; i[iy] < n[iy]; ++i[iy]) {
        i[iz] = n[iz] - 1;
        N c = get(i[0], i[1], i[2]);
        while (--i[iz] >= 0) c = c * v + get(i[0], i[1], i[2]);

        if (c.sign() != 0) {
          m2.push_back(i[ix]);
          m2.push_back(i[iy]);
          a2.push_back(c);
        }
      }

    Poly2<N> p2 = Poly2<N>(a2.size(), &a2[0], &m2[0]);

#ifdef BLEEN
    // debug
    N xyz[3];
    xyz[ix] = N::constant(2.781828);
    xyz[iy] = N::constant(-3.1415926);
    xyz[iz] = v;
    N p2v = p2.value(xyz[ix], xyz[iy]);
    N p3v = value(xyz[0], xyz[1], xyz[2]);
    assert((p2v - p3v).sign(false) == 0);
#endif

    return p2;
  }

  Poly2D<N> substitute1d(int ixyz, const N& v) const {
    int n[3] = {nx, ny, nz};
    int i[3] = {0, 0, 0};

    int iz = ixyz, ix = (iz + 1) % 3, iy = (ix + 1) % 3;

    Poly2D<N> p2(n[ix]-1, n[iy]-1);

    for (i[ix] = 0; i[ix] < n[ix]; ++i[ix])
      for (i[iy] = 0; i[iy] < n[iy]; ++i[iy]) {
        i[iz] = n[iz] - 1;

        N c = get(i[0], i[1], i[2]);
        for (i[iz] = n[iz] - 1; --i[iz] >= 0;)
          c = c * v + get(i[0], i[1], i[2]);

        p2.set(i[ix], i[iy], c);
      }

    // debug
    N xyz[3];
    xyz[ix] = N::constant(2.781828);
    xyz[iy] = N::constant(-3.1415926);
    xyz[iz] = v;
    N p2v = p2.value(xyz[ix], xyz[iy]);
    N p3v = value(xyz[0], xyz[1], xyz[2]);
    assert((p2v - p3v).sign(false) == 0);

    return p2;
  }

  Poly3D<N> der(int xyz) {
    int n[3] = {nx, ny, nz};
    int iz = xyz, ix = (iz + 1) % 3, iy = (ix + 1) % 3;
    --n[iz];
    Poly3D<N> poly(n[0] - 1, n[1] - 1, n[2] - 1);
    ++n[iz];
    int i[3] = {0, 0, 0};
    for (i[ix] = 0; i[ix] < n[ix]; ++i[ix])
      for (i[iy] = 0; i[iy] < n[iy]; ++i[iy])
        for (i[iz] = 1; i[iz] < n[iz]; ++i[iz]) {
          N c = get(i[0], i[1], i[2]) * i[iz];
          --i[iz];
          poly.set(i[0], i[1], i[2], c);
          ++i[iz];
        }
    return poly;
  }

  N der(int xyz, const PV3<N>& v) {
    int n[3] = {nx, ny, nz};
    int iz = xyz, ix = (iz + 1) % 3, iy = (ix + 1) % 3;
    int i[3] = {0, 0, 0};
    N sx = N::constant(0);
    for (i[ix] = n[ix]; --i[ix] >= 0;) {
      N sy = N::constant(0);
      for (i[iy] = n[iy]; --i[iy] >= 0;) {
        N sz = N::constant(0);
        for (i[iz] = n[iz]; --i[iz] >= 1;)
          sz = v[iz] * sz + get(i[0], i[1], i[2]) * i[iz];
        sy = v[iy] * sy + sz;
      }
      sx = v[ix] * sx + sy;
    }

    // debug
    Poly3D<N> df = der(xyz);
    N d = df.value(v);
    assert((d - sx).sign(false) == 0);

    return sx;
  }

  PV3<N> gradient(const PV3<N>& v) {
    return PV3<N>(der(0, v), der(1, v), der(2, v));
  }

  Poly3D<N> operator+(const Poly3D<N>& b) const {
    Poly3D<N> p(std::max(degX(), b.degX()), std::max(degY(), b.degY()),
                std::max(degZ(), b.degZ()));

    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
        for (int iz = 0; iz < nz; ++iz) p.set(ix, iy, iz, get(ix, iy, iz));

    for (int ix = 0; ix < b.nx; ++ix)
      for (int iy = 0; iy < b.ny; ++iy)
        for (int iz = 0; iz < b.nz; ++iz)
          p.set(ix, iy, iz, p.get(ix, iy, iz) + b.get(ix, iy, iz));

    return p;
  }

  Poly3D<N> operator-(const Poly3D<N>& b) const {
    Poly3D<N> p(std::max(degX(), b.degX()), std::max(degY(), b.degY()),
                std::max(degZ(), b.degZ()));

    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
        for (int iz = 0; iz < nz; ++iz) p.set(ix, iy, iz, get(ix, iy, iz));

    for (int ix = 0; ix < b.nx; ++ix)
      for (int iy = 0; iy < b.ny; ++iy)
        for (int iz = 0; iz < b.nz; ++iz)
          p.set(ix, iy, iz, p.get(ix, iy, iz) - b.get(ix, iy, iz));

    return p;
  }

  Poly3D<N> operator*(const Poly3D<N>& b) const {
    Poly3D<N> p(degX() + b.degX(), degY() + b.degY(), degZ() + b.degZ());

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

  void mobius(int xyz, const N& lb, const N& ub) {
    int n[3] = {nx, ny, nz}, i[3];
    int iz = xyz, ix = (iz + 1) % 3, iy = (ix + 1) % 3;

    for (int m = 0; m < n[iz] - 1; m++)
      for (i[iz] = n[iz] - 1 - 1; i[iz] >= m; i[iz]--)
        for (i[ix] = 0; i[ix] < n[ix]; i[ix]++)
          for (i[iy] = 0; i[iy] < n[iy]; i[iy]++) {
            N c = get(i[0], i[1], i[2]);
            ++i[iz];
            c = c + lb * get(i[0], i[1], i[2]);
            --i[iz];
            set(i[0], i[1], i[2], c);
          }

    N s = 1.0 / (ub - lb);
    for (int m = 0; m < n[iz] - 1; m++)
      for (i[iz] = 1; i[iz] <= n[iz] - 1 - m; i[iz]++)
        for (i[ix] = 0; i[ix] < n[ix]; i[ix]++)
          for (i[iy] = 0; i[iy] < n[iy]; i[iy]++) {
            N c = get(i[0], i[1], i[2]);
            --i[iz];
            c = c + s * get(i[0], i[1], i[2]);
            ++i[iz];
            set(i[0], i[1], i[2], c);
          }
  }

  void mobius(const PV3<Parameter>& box) {
    PV3<N> b(box);
    // for (int xyz = 0; xyz < 3; ++xyz) mobius(xyz, b[xyz].lbP(),
    // b[xyz].ubP());
    for (int xyz = 0; xyz < 3; ++xyz)
      mobius(xyz, N::constant(b[xyz].lb()), N::constant(b[xyz].ub()));
  }

  // t is a flattened 3x3 matrix in column major form
  Poly3D<N> rotate(const double t[9]) const {
    Poly3D<N> x(1,1,1);
    Poly3D<N> y(1,1,1);
    Poly3D<N> z(1,1,1);

    Poly3D<N> unit(0,0,0);

    // First column
    x.add(1, 0, 0, N::constant(t[0]));
    y.add(1, 0, 0, N::constant(t[1]));
    z.add(1, 0, 0, N::constant(t[2]));

    // Second column
    x.add(0, 1, 0, N::constant(t[3]));
    y.add(0, 1, 0, N::constant(t[4]));
    z.add(0, 1, 0, N::constant(t[5]));

    // Third column
    x.add(0, 0, 1, N::constant(t[6]));
    y.add(0, 0, 1, N::constant(t[7]));
    z.add(0, 0, 1, N::constant(t[8]));

    unit.add(0, 0, 0, N::constant(1));

    std::vector<Poly3D<N>> xpow(nx);
    std::vector<Poly3D<N>> ypow(ny);
    std::vector<Poly3D<N>> zpow(nz);


    xpow[0] = ypow[0] = zpow[0] = unit;

    for (int i = 1; i < nx; i++) xpow[i] = x * xpow[i - 1];
    for (int i = 1; i < ny; i++) ypow[i] = y * ypow[i - 1];
    for (int i = 1; i < nz; i++) zpow[i] = z * zpow[i - 1];

    Poly3D<N> sum;

    for(int i = 0; i < nx; ++i) {
      const Poly3D<N>& xpoly = xpow[i];
      for(int j = 0; j < ny; ++j) {
        const Poly3D<N>& xypoly = xpoly * ypow[j];
        for(int k = 0; k < nz; ++k) {

          Poly3D<N> coef(0,0,0);
          coef.add(0, 0, 0, get(i, j, k));

          sum = sum + (coef * xypoly * zpow[k]);

        }
      }
    }

    return sum;
  }

};

vector<PTR<Object<PV3>>> getRoots(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g,
                                  PTR<Object<Poly3D>> h, PV3<Parameter> box);

class Poly3Dval : public Primitive {
  PTR<Object<Poly3D>> f;
  PTR<Object<PV3>> xyz;
  DeclareSign {
    Poly3D<N> f = this->f->get<N>();
    PV3<N> xyz = this->xyz->get<N>();
    return f.value(xyz);
  }

 public:
  Poly3Dval(PTR<Object<Poly3D>> f, PTR<Object<PV3>> xyz) : f(f), xyz(xyz) {}
};

class Sub1Poly3D : public Object<Poly2> {
  PTR<Object<Poly3D>> p3d;
  int ixyz;
  PTR<Object<Scalar>> xyz;
  DeclareCalculate(Poly2) {
    return p3d->get<N>().substitute1(ixyz, xyz->get<N>());
  }

 public:
  // ixyz = 0 or 1 or 2 means substituting for x or y or z
  Sub1Poly3D(PTR<Object<Poly3D>> p3d, int ixyz, PTR<Object<Scalar>> xyz)
      : p3d(p3d), ixyz(ixyz), xyz(xyz) {}
};

class Sub1dPoly3D : public Object<Poly2D> {
  PTR<Object<Poly3D>> p3d;
  int ixyz;
  PTR<Object<Scalar>> xyz;
  DeclareCalculate(Poly2D) {
    return p3d->get<N>().substitute1d(ixyz, xyz->get<N>());
  }

 public:
  // ixyz = 0 or 1 or 2 means substituting for x or y or z
  Sub1dPoly3D(PTR<Object<Poly3D>> p3d, int ixyz, PTR<Object<Scalar>> xyz)
      : p3d(p3d), ixyz(ixyz), xyz(xyz) {}
};

class P2to3 : public Object<PV3> {
  PTR<Object<PV2>> p2;
  PTR<Object<Scalar>> s;
  int xyz;
  DeclareCalculate(PV3) {
    PV3<N> p;
    p[xyz] = s->get<N>();
    PV2<N> p2 = this->p2->get<N>();
    p[(xyz + 1) % 3] = p2.x;
    p[(xyz + 2) % 3] = p2.y;
    return p;
  }

 public:
  P2to3(PTR<Object<PV2>> p2, PTR<Object<Scalar>> s, int xyz)
      : p2(p2), s(s), xyz(xyz) {}
};

template <class N>
PV3<N> solve(N a[3][4]) {
  int r[3] = {0, 1, 2};
  for (int i = 1; i < 3; i++)
    if (lessMag(a[r[0]][0], a[r[i]][0])) swap(r[0], r[i]);
  for (int i = 1; i < 3; i++) {
    N c = a[r[i]][0] / a[r[0]][0];
    a[r[i]][0] = N::constant(0);
    for (int j = 1; j < 4; j++) a[r[i]][j] = a[r[i]][j] - c * a[r[0]][j];
  }
  if (lessMag(a[r[1]][1], a[r[2]][1])) swap(r[1], r[2]);
  N c = a[r[2]][1] / a[r[1]][1];
  a[r[2]][1] = N::constant(0);
  for (int j = 2; j < 4; j++) a[r[2]][j] = a[r[2]][j] - c * a[r[1]][j];
  PV3<N> x;
  x[2] = -a[r[2]][3] / a[r[2]][2];
  x[1] = -(a[r[1]][3] + a[r[1]][2] * x[2]) / a[r[1]][1];
  x[0] = -(a[r[0]][3] + a[r[0]][2] * x[2] + a[r[0]][1] * x[1]) / a[r[0]][0];
  return x;
}

template <class N>
PV3<N> newton(Poly3D<N> f, Poly3D<N> g, Poly3D<N> h, PV3<N> r) {
  PV3<N> m = PV3<N>(r.x.midP(), r.y.midP(), r.z.midP());
  N a[3][4] = {{f.der(0, r), f.der(1, r), f.der(2, r), f.value(m)},
               {g.der(0, r), g.der(1, r), g.der(2, r), g.value(m)},
               {h.der(0, r), h.der(1, r), h.der(2, r), h.value(m)}};
  PV3<N> s = m + solve(a);
  return PV3<N>(r.x.intersect(s.x), r.y.intersect(s.y), r.z.intersect(s.z));
}

PV3<Parameter> subdivide(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g,
                         PTR<Object<Poly3D>> h, PV3<Parameter>& box);

class Root3D : public Object<PV3> {
  PTR<Object<Poly3D>> f, g, h;

 public:
  PV3<Parameter> box;
  DeclareCalculate(PV3) {
#ifdef SET_ALG
    if (std::is_same<N, PParameter>::value) {
      throw SignException(true);
    }
#else
    if (std::is_same<N, PParameter>::value) {
      throw SignException();
    }
#endif

    Poly3D<N> f = this->f->get<N>();
    Poly3D<N> g = this->g->get<N>();
    Poly3D<N> h = this->h->get<N>();
    PV3<N> b = uninitialized() ? PV3<N>(box) : PV3<N>(getCurrentP());

    for (int i = 0; i < 150; i++) {
      PV3<N> bp;
      int count = 0;
      do {
        ++count;
        bp = b;
        try {
          b = newton(f, g, h, b);
        } catch (const SignException& e) {
          cout << e.what() << endl;
        }
      } while (b.x.subset(bp.x) || b.y.subset(bp.y) || b.z.subset(bp.z));

      if (curPrecision() == 53u || (curPrecision() > 53u && count > 1)) break;

      if (curPrecision() > 53u) {
        unsigned int savePrecision = curPrecision();
        curPrecision() = 53u;
        box = subdivide(this->f, this->g, this->h, box);
        b = PV3<N>(box);
        curPrecision() = savePrecision;
      }
    }

#ifdef SET_ALG
    b.x.setAlg();
    b.y.setAlg();
    b.z.setAlg();
#endif
    return b;
  }

 public:
  Root3D(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g, PTR<Object<Poly3D>> h,
         PV3<Parameter> box)
      : f(f), g(g), h(h), box(box) {}
};

class Descartes3D : public Primitive {
  PTR<Object<Poly3D>> f;
  PV3<Parameter> box;
  DeclareSign {
    Poly3D<N> f = this->f->get<N>();
    f.mobius(box);
    int s = f[0].sign();
    for (int i = 1; i < f.size(); ++i) {
      int si = f[i].sign();
      if (si != 0 && si != s) return N::constant(0);
    }
    return f[0];
  }

 public:
  Descartes3D(PTR<Object<Poly3D>> f, PV3<Parameter> box) : f(f), box(box) {}
};

// functions

vector<PTR<Object<PV3>>> boundaryRoots(PTR<Object<Poly3D>> f,
                                       PTR<Object<Poly3D>> g,
                                       const PV3<Parameter>& box, int xyz,
                                       int lu);

vector<PTR<Object<PV3>>> boundaryRoots(PTR<Object<Poly3D>> f,
                                       PTR<Object<Poly3D>> g,
                                       const PV3<Parameter>& box);

vector<PTR<Object<PV3>>> getRoots(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g,
                                  PTR<Object<Poly3D>> h, PV3<Parameter> box);

#endif
