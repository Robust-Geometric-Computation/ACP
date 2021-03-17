#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

#include "acp/poly/poly3.h"
using namespace acp;

class Sub2Poly : public Object<Poly> {
  PTR<Object<Poly3>> ppoly3;
  int ixyz;
  PTR<Object<Scalar>> x, y, z;
  Poly<Parameter> calculate();

 public:
  // ixyz = 0 or 1 or 2 means substituting for y,z or x,z or x,y
  Sub2Poly(PTR<Object<Poly3>> ppoly3, int ixyz, PTR<Object<Scalar>> xyz[3])
      : ppoly3(ppoly3), ixyz(ixyz), x(xyz[0]), y(xyz[1]), z(xyz[2]) {}
  DeclareCalculate(Poly) {
    return ppoly3->get<N>().substitute2(ixyz, x->get<N>().x, y->get<N>().x,
                                        z->get<N>().x);
  }
};

class P3 {
 public:
  int xd, yd, zd;
  Parameter* a;  // size (xd+1)*(yd+1)*(zd+1)

  Parameter get(int i, int j, int k) {
    return a[i + (xd + 1) * (j + (yd + 1) * k)];
  }

  void set(int i, int j, int k, Parameter c) {
    a[i + (xd + 1) * (j + (yd + 1) * k)] = c;
  }

  P3(int xd, int yd, int zd) : xd(xd), yd(yd), zd(zd) {
    a = new Parameter[(xd + 1) * (yd + 1) * (zd + 1)];
    for (int i = 0; i <= xd; ++i) {
      for (int j = 0; j <= yd; ++j) {
        for (int k = 0; k <= zd; ++k) set(i, j, k, Parameter(0));
      }
    }
  }

  P3(Poly3<Parameter>& p) : P3(p.xd, p.yd, p.zd) {
    for (auto& term : p.ia)
      set(term.first.i[0], term.first.i[1], term.first.i[2], p.a[term.second]);
  }

  P3(P3& p3) : xd(p3.xd), yd(p3.yd), zd(p3.zd) {
    a = new Parameter[(xd + 1) * (yd + 1) * (zd + 1)];
    for (int i = 0; i <= xd; ++i) {
      for (int j = 0; j <= yd; ++j) {
        for (int k = 0; k <= zd; ++k) set(i, j, k, p3.get(i, j, k));
      }
    }
  }

  ~P3() { delete[] a; }

  void print() {
    for (int k = 0; k <= zd; ++k) {
      for (int i = 0; i <= xd; ++i)
        for (int j = 0; j <= yd; ++j) cout << get(i, j, k).mid() << " ";
      cout << endl;
    }
  }

  void mobiusX(double xld, double xud) {
    Parameter xl = Parameter::constant(xld);
    Parameter xu = Parameter::constant(xud);
    Parameter s = xl;
    for (int i = 0; i < xd; i++)
      for (int j = xd - 1; j >= i; j--)
        for (int k = 0; k <= yd; k++)
          for (int l = 0; l <= zd; l++)
            set(j, k, l, get(j, k, l) + s * get(j + 1, k, l));

    s = 1.0 / (xu - xl);
    for (int i = 0; i < xd; i++)
      for (int j = xd - 1; j >= i; j--)
        for (int k = 0; k <= yd; k++)
          for (int l = 0; l <= zd; l++)
            set(xd - j, k, l, get(xd - j, k, l) + s * get(xd - (j + 1), k, l));
  }

  void mobiusY(double yld, double yud) {
    Parameter yl = Parameter::constant(yld);
    Parameter yu = Parameter::constant(yud);
    Parameter s = yl;
    for (int i = 0; i < yd; i++)
      for (int j = yd - 1; j >= i; j--)
        for (int k = 0; k <= zd; k++)
          for (int l = 0; l <= xd; l++)
            set(l, j, k, get(l, j, k) + s * get(l, j + 1, k));

    s = 1.0 / (yu - yl);
    for (int i = 0; i < yd; i++)
      for (int j = yd - 1; j >= i; j--)
        for (int k = 0; k <= zd; k++)
          for (int l = 0; l <= xd; l++)
            set(l, yd - j, k, get(l, yd - j, k) + s * get(l, yd - (j + 1), k));
  }

  Parameter valueZ(double z) {
    Parameter v = get(0, 0, zd);
    for (int i = zd; --i >= 0;) v = v * z + get(0, 0, i);
    return v;
  }

  // random value in [a,b] with largest magnitude valueZ
  // n additional trials
  double bestZ(double a, double b, int n) {
    double m = randomNumber(a, b);
    Parameter pm = valueZ(m);
    for (int i = 0; i < n; i++) {
      double m2 = randomNumber(a, b);
      Parameter pm2 = valueZ(m2);
      if ((pm.lb() <= 0 && pm.ub() >= 0) ||
          (pm.lb() > 0 && pm2.lb() > pm.ub()) ||
          (pm.ub() < 0 && pm2.ub() < pm.lb())) {
        m = m2;
        pm = pm2;
      }
    }
    return m;
  }

  void mobiusZ(double zld, double zud) {
    Parameter zl = Parameter::constant(zld);
    Parameter zu = Parameter::constant(zud);
    Parameter s = zl;
    for (int i = 0; i < zd; i++)
      for (int j = zd - 1; j >= i; j--)
        for (int k = 0; k <= xd; k++)
          for (int l = 0; l <= yd; l++)
            set(k, l, j, get(k, l, j) + s * get(k, l, j + 1));

    s = 1.0 / (zu - zl);
    for (int i = 0; i < zd; i++)
      for (int j = zd - 1; j >= i; j--)
        for (int k = 0; k <= xd; k++)
          for (int l = 0; l <= yd; l++)
            set(k, l, zd - j, get(k, l, zd - j) + s * get(k, l, zd - (j + 1)));
  }

  int signXY(double z) {
    int s = 2;

    for (int i = 0; i <= xd; ++i) {
      for (int j = 0; j <= yd; ++j) {
        Parameter cij = get(i, j, zd);
        for (int k = zd; --k >= 0;) cij = cij * z + get(i, j, k);
        if (s == 2) {
          s = cij.sign(false);
          if (s == 0) return 0;
        } else if (s != cij.sign(false))
          return 0;
      }
    }

    return s;
  }

  int signZ(int iz) {
    int s = 2;

    for (int i = 0; i <= xd; ++i) {
      for (int j = 0; j <= yd; ++j) {
        if (s == 2) {
          if (get(i, j, iz).ub() < 0)
            s = -1;
          else if (get(i, j, iz).lb() > 0)
            s = 1;
          else if (get(i, j, iz).lb() < 0 && get(i, j, iz).ub() > 0)
            return 0;
        } else if (s == 1) {
          if (get(i, j, iz).lb() < 0) return 0;
        } else {
          if (get(i, j, iz).ub() > 0) return 0;
        }
      }
    }

    if (s == 2) return 0;

    return s;
  }

  int descartesZ() {
    int count = 0;
    int s = signZ(0);
    if (s == 0) return -1;
    for (int iz = 1; iz <= zd; iz++) {
      int si = signZ(iz);
      if (si == 0) return -1;
      if (si == -s) {
        s = -s;
        count++;
      }
    }
    return count;
  }
};

void probeZbug(const char* filename, int nsteps, Poly3<Parameter>& poly3) {
  ofstream out;
  out.open(filename);
  out << "P3" << endl;
  out << nsteps << " " << nsteps << endl;
  out << 255 << endl;

  cout << "probeZ " << filename << " " << nsteps << endl;
  PTR<Object<Poly3>> opoly3 = new Object<Poly3>(poly3);
  vector<double> xs(nsteps + 1);
  for (int ix = 0; ix <= nsteps; ++ix) xs[ix] = -1.0 + 2.0 * ix / nsteps;

  P3 p3(poly3);

  for (int ix = 0; ix < nsteps; ++ix) {
    cout << "ix " << ix << endl;
    P3 p3x(p3);
    p3x.mobiusX(xs[ix], xs[ix + 1]);
    for (int iy = 0; iy < nsteps; ++iy) {
      // cout << "iy " << iy << endl;
      P3 p3xy(p3x);
      p3xy.mobiusY(xs[iy], xs[iy + 1]);

      // if (ix == 256 && iy == 256)
      // p3xy.print();

      double zl = randomNumber(-1.1, -1);
      double zu = randomNumber(1, 1.1);
      vector<pair<double, double>> ints;
      ints.push_back(make_pair(zl, zu));

      Parameter root;

      while (ints.size() > 0) {
        pair<double, double> zz = ints.back();
        ints.pop_back();
        if (zz.first > 1.0 || zz.second < -1.0) continue;
        P3 p3xyz(p3xy);
        p3xyz.mobiusZ(zz.first, zz.second);
        int d = p3xyz.descartesZ();
        if (ix == 327 && iy == 499)
          cout << "zz " << zz.first << " " << zz.second << " d " << d << endl;
        if (d == 0) continue;
        if (d == 1) {
          PTR<Object<Scalar>> xyz[3];
          xyz[0] = new Object<Scalar>(
              Parameter::constant((xs[ix] + xs[ix + 1]) / 2));
          xyz[1] = new Object<Scalar>(
              Parameter::constant((xs[iy] + xs[iy + 1]) / 2));
          xyz[2] = xyz[0];
          PTR<Object<Poly>> poly = new Sub2Poly(opoly3, 2, xyz);
          PTR<Object<Scalar>> l =
              new Object<Scalar>(Parameter::constant(zz.first));
          PTR<Object<Scalar>> u =
              new Object<Scalar>(Parameter::constant(zz.second));
          vector<PTR<Object<Scalar>>> roots = getRoots(poly, l, u);
          assert(roots.size() == 1);
          Parameter r = roots[0]->getApprox();
          if (r.ub() < -1.0 || r.lb() > 1.0) continue;
          root = r;
          break;
        }
        if (zz.second - zz.first < 1.0 / 1024.0) {
          root = Parameter::constant(zz.first);
          // p3xyz.print();
          break;
        }
        double mid = (zz.first + zz.second) / 2;
        double dif = (zz.second - zz.first) / 8;
        double ma = mid - dif;
        double mb = mid + dif;
        double m = randomNumber(ma, mb);
        Parameter pm = p3xy.valueZ(m);
        for (int i = 0; i < 0; i++) {
          double m2 = randomNumber(ma, mb);
          Parameter pm2 = p3xy.valueZ(m2);
          if ((pm.lb() <= 0 && pm.ub() >= 0) ||
              (pm.lb() > 0 && pm2.lb() > pm.ub()) ||
              (pm.ub() < 0 && pm2.ub() < pm.lb())) {
            m = m2;
            pm = pm2;
          }
        }
        ints.push_back(make_pair(m, zz.second));
        ints.push_back(make_pair(zz.first, m));
      }

      int color[3] = {0, 0, 0};
      if (root.uninitialized()) {
        out << color[0] << " " << color[1] << " " << color[2] << endl;
        break;
      } else if (root.lb() < root.ub()) {
        Parameter x = Parameter::constant((xs[ix] + xs[ix + 1]) / 2);
        Parameter y = Parameter::constant((xs[iy] + xs[iy + 1]) / 2);
        Parameter z = root;
        PV3<Parameter> u = poly3.gradient(x, y, z).unit();
        double uz = u.z.mid();
        int c = (-uz + 1) * 255 / 2;
        assert(c <= 255 && c >= 0);
        color[0] = color[1] = color[2] = c;
        out << color[0] << " " << color[1] << " " << color[2] << endl;
        break;
      } else {
        color[0] = 255;
        int w = 20;
        // if (256 - w < ix && ix < 256 + w && 256 - w < iy && iy < 256 + w)
        // cout << "ambiguity " << ix << " " << iy << endl;

        Parameter x = Parameter::constant((xs[ix] + xs[ix + 1]) / 2);
        Parameter y = Parameter::constant((xs[iy] + xs[iy + 1]) / 2);
        Parameter z = root;
        PV3<Parameter> u = poly3.gradient(x, y, z).unit();
        double uz = u.z.mid();
        static double maxUZ;
        if (fabs(uz) > maxUZ) {
          cout << "maxUZ " << uz << " ix " << ix << " iy " << iy << endl;
          maxUZ = fabs(uz);
        }
        out << color[0] << " " << color[1] << " " << color[2] << endl;
        break;
      }
    }
  }

  out.close();
}

void probeZ(const char* filename, int nsteps, Poly3<Parameter>& poly3) {
  ofstream out;
  out.open(filename);
  out << "P3" << endl;
  out << nsteps << " " << nsteps << endl;
  out << 255 << endl;

  cout << "probeZ " << filename << " " << nsteps << endl;
  PTR<Object<Poly3>> opoly3 = new Object<Poly3>(poly3);
  vector<double> xs(nsteps + 1);
  for (int ix = 0; ix <= nsteps; ++ix) xs[ix] = -1.0 + 2.0 * ix / nsteps;

  P3 p3(poly3);

  int shortIntervals = 0;
  int numRoots = 0;
  double shortLength = 0;
  int badInterval = 0;

  for (int ix = 0; ix < nsteps; ++ix) {
    cout << "ix " << ix << endl;
    P3 p3x(p3);
    p3x.mobiusX(xs[ix], xs[ix + 1]);
    for (int iy = 0; iy < nsteps; ++iy) {
      // cout << "iy " << iy << endl;
      P3 p3xy(p3x);
      p3xy.mobiusY(xs[iy], xs[iy + 1]);

      // if (ix == 256 && iy == 256)
      // p3xy.print();

      for (int itrial = 0; itrial < 1; itrial++) {
        // double zl = randomNumber(-1.1, -1);
        // double zu = randomNumber(1, 1.1);
        double zl = p3xy.bestZ(-1.1, -1, 10);
        double zu = p3xy.bestZ(1, 1.1, 10);
        vector<pair<double, double>> ints;
        ints.push_back(make_pair(zl, zu));

        Parameter root;

        while (ints.size() > 0) {
          pair<double, double> zz = ints.back();
          ints.pop_back();
          if (zz.first > 1.0 || zz.second < -1.0) continue;
          P3 p3xyz(p3xy);
          p3xyz.mobiusZ(zz.first, zz.second);
          int d = p3xyz.descartesZ();
          if (ix == 327 && iy == 499)
            cout << "zz " << zz.first << " " << zz.second << " d " << d << endl;
          if (d == 0) continue;
          if (d == 1) {
            PTR<Object<Scalar>> xyz[3];
            xyz[0] = new Object<Scalar>(
                Parameter::constant((xs[ix] + xs[ix + 1]) / 2));
            xyz[1] = new Object<Scalar>(
                Parameter::constant((xs[iy] + xs[iy + 1]) / 2));
            xyz[2] = xyz[0];
            PTR<Object<Poly>> poly = new Sub2Poly(opoly3, 2, xyz);
            PTR<Object<Scalar>> l =
                new Object<Scalar>(Parameter::constant(zz.first));
            PTR<Object<Scalar>> u =
                new Object<Scalar>(Parameter::constant(zz.second));
            vector<PTR<Object<Scalar>>> roots = getRoots(poly, l, u);
            assert(roots.size() == 1);
            Parameter r = roots[0]->getApprox();
            if (r.ub() < -1.0 || r.lb() > 1.0) continue;
            if (root.uninitialized())  // multiroot
              root = r;
            numRoots++;

            PV3<Parameter> g = poly3.gradient(
                Parameter::constant((xs[ix] + xs[ix + 1]) / 2),
                Parameter::constant((xs[iy] + xs[iy + 1]) / 2), r);
            if (g.z.sign(false) == 0) {
              cerr << "Cannot calculate z's." << endl;
              continue;
            }

            Parameter dxdz = g.x / g.z;
            Parameter dydz = g.y / g.z;
            Parameter zzz[4] = {r + dxdz * 1 / nsteps + dydz * 1 / nsteps,
                                r - dxdz * 1 / nsteps + dydz * 1 / nsteps,
                                r + dxdz * 1 / nsteps - dydz * 1 / nsteps,
                                r - dxdz * 1 / nsteps - dydz * 1 / nsteps};
            double lo = r.lb();
            double up = r.ub();
            for (int i = 0; i < 4; i++) {
              if (lo > zzz[i].lb()) lo = zzz[i].lb();
              if (up < zzz[i].ub()) up = zzz[i].ub();
            }

            double width1 = up - lo;

            int slo;
            do {
              lo += lo - r.lb();
              if (lo < zz.first) lo = zz.first;
            } while ((slo = p3xy.signXY(lo)) == 0);

            int sup;
            do {
              up += up - r.ub();
              if (up > zz.second) up = zz.second;
            } while ((sup = p3xy.signXY(up)) == 0);

            assert(slo * sup == -1);

            int smid;
            double a = lo;
            double b = r.lb();
            while (b - a > 1 / 2.0 / nsteps)
              if ((smid = p3xy.signXY((a + b) / 2)) == 0)
                b = (a + b) / 2;
              else {
                assert(slo == smid);
                a = (a + b) / 2;
              }
            lo = a;

            a = r.ub();
            b = up;
            while (b - a > 1 / 2.0 / nsteps)
              if ((smid = p3xy.signXY((a + b) / 2)) == 0)
                a = (a + b) / 2;
              else {
                assert(sup == smid);
                b = (a + b) / 2;
              }
            up = b;

            double width2 = up - lo;
            static double widthRat;
            if (widthRat < width2 / width1) {
              // cerr << "widthRat " << width2 / width1 << endl;
              widthRat = width2 / width1;
            }
            cout << ((int)(0.5 + 10 * width2 / width1)) / 10.0 << " widthRat"
                 << endl;

            P3 p3xyz(p3xy);
            p3xyz.mobiusZ(lo, up);
            int d = p3xyz.descartesZ();
            if (d != 1) {
              // cerr << "Bad interval " << d << " " << up - lo << endl;
              badInterval++;
            }

            continue;
          }
          if (d < 0 && zz.second - zz.first < 1.0 / 8.0 /*512.0*/) {
            shortLength += zz.second - zz.first;
            root = Parameter::constant(zz.first);
            // p3xyz.print();
            shortIntervals++;
            continue;
          }

          double mid = (zz.first + zz.second) / 2;
          double dif = (zz.second - zz.first) / 8;
          double m = p3xy.bestZ(mid - dif, mid + dif, 10);

          ints.push_back(make_pair(m, zz.second));
          ints.push_back(make_pair(zz.first, m));
        }

        int color[3] = {0, 0, 0};
        if (root.uninitialized()) {
          out << color[0] << " " << color[1] << " " << color[2] << endl;
          break;
        } else if (root.lb() < root.ub()) {
          Parameter x = Parameter::constant((xs[ix] + xs[ix + 1]) / 2);
          Parameter y = Parameter::constant((xs[iy] + xs[iy + 1]) / 2);
          Parameter z = root;
          PV3<Parameter> u = poly3.gradient(x, y, z).unit();
          double uz = u.z.mid();
          int c = (-uz + 1) * 255 / 2;
          assert(c <= 255 && c >= 0);
          color[0] = color[1] = color[2] = c;
          out << color[0] << " " << color[1] << " " << color[2] << endl;
          break;
        } else if (itrial == 0) {
          color[0] = 255;
          int w = 20;
          // if (256 - w < ix && ix < 256 + w && 256 - w < iy && iy < 256 + w)
          // cout << "ambiguity " << ix << " " << iy << endl;

          Parameter x = Parameter::constant((xs[ix] + xs[ix + 1]) / 2);
          Parameter y = Parameter::constant((xs[iy] + xs[iy + 1]) / 2);
          Parameter z = root;
          PV3<Parameter> u = poly3.gradient(x, y, z).unit();
          double uz = u.z.mid();
          static double maxUZ;
          if (fabs(uz) > maxUZ) {
            cout << "maxUZ " << uz << " ix " << ix << " iy " << iy << endl;
            maxUZ = fabs(uz);
          }
          out << color[0] << " " << color[1] << " " << color[2] << endl;
          break;
        }
      }
    }
  }

  out.close();

  cout << "shortIntervals " << shortIntervals << endl;
  cout << "shortLength " << shortLength << endl;
  cout << "numRoots " << numRoots << endl;
  cout << "badInterval " << badInterval << endl;
}

void probeZbestz(const char* filename, int nsteps, Poly3<Parameter>& poly3) {
  ofstream out;
  out.open(filename);
  out << "P3" << endl;
  out << nsteps << " " << nsteps << endl;
  out << 255 << endl;

  cout << "probeZ " << filename << " " << nsteps << endl;
  PTR<Object<Poly3>> opoly3 = new Object<Poly3>(poly3);
  vector<double> xs(nsteps + 1);
  for (int ix = 0; ix <= nsteps; ++ix) xs[ix] = -1.0 + 2.0 * ix / nsteps;

  P3 p3(poly3);

  for (int ix = 0; ix < nsteps; ++ix) {
    cout << "ix " << ix << endl;
    P3 p3x(p3);
    p3x.mobiusX(xs[ix], xs[ix + 1]);
    for (int iy = 0; iy < nsteps; ++iy) {
      // cout << "iy " << iy << endl;
      P3 p3xy(p3x);
      p3xy.mobiusY(xs[iy], xs[iy + 1]);

      // if (ix == 256 && iy == 256)
      // p3xy.print();

      for (int itrial = 0; itrial < 1; itrial++) {
        // double zl = randomNumber(-1.1, -1);
        // double zu = randomNumber(1, 1.1);
        double zl = p3xy.bestZ(-1.1, -1, 10);
        double zu = p3xy.bestZ(1, 1.1, 10);
        vector<pair<double, double>> ints;
        ints.push_back(make_pair(zl, zu));

        Parameter root;

        while (ints.size() > 0) {
          pair<double, double> zz = ints.back();
          ints.pop_back();
          if (zz.first > 1.0 || zz.second < -1.0) continue;
          P3 p3xyz(p3xy);
          p3xyz.mobiusZ(zz.first, zz.second);
          int d = p3xyz.descartesZ();
          if (ix == 327 && iy == 499)
            cout << "zz " << zz.first << " " << zz.second << " d " << d << endl;
          if (d == 0) continue;
          if (d == 1) {
            PTR<Object<Scalar>> xyz[3];
            xyz[0] = new Object<Scalar>(
                Parameter::constant((xs[ix] + xs[ix + 1]) / 2));
            xyz[1] = new Object<Scalar>(
                Parameter::constant((xs[iy] + xs[iy + 1]) / 2));
            xyz[2] = xyz[0];
            PTR<Object<Poly>> poly = new Sub2Poly(opoly3, 2, xyz);
            PTR<Object<Scalar>> l =
                new Object<Scalar>(Parameter::constant(zz.first));
            PTR<Object<Scalar>> u =
                new Object<Scalar>(Parameter::constant(zz.second));
            vector<PTR<Object<Scalar>>> roots = getRoots(poly, l, u);
            assert(roots.size() == 1);
            Parameter r = roots[0]->getApprox();
            if (r.ub() < -1.0 || r.lb() > 1.0) continue;
            root = r;
            break;
          }
          if (zz.second - zz.first < 1.0 / 1024.0) {
            root = Parameter::constant(zz.first);
            // p3xyz.print();
            break;
          }

          double mid = (zz.first + zz.second) / 2;
          double dif = (zz.second - zz.first) / 8;
          double m = p3xy.bestZ(mid - dif, mid + dif, 10);

          ints.push_back(make_pair(m, zz.second));
          ints.push_back(make_pair(zz.first, m));
        }

        int color[3] = {0, 0, 0};
        if (root.uninitialized()) {
          out << color[0] << " " << color[1] << " " << color[2] << endl;
          break;
        } else if (root.lb() < root.ub()) {
          Parameter x = Parameter::constant((xs[ix] + xs[ix + 1]) / 2);
          Parameter y = Parameter::constant((xs[iy] + xs[iy + 1]) / 2);
          Parameter z = root;
          PV3<Parameter> u = poly3.gradient(x, y, z).unit();
          double uz = u.z.mid();
          int c = (-uz + 1) * 255 / 2;
          assert(c <= 255 && c >= 0);
          color[0] = color[1] = color[2] = c;
          out << color[0] << " " << color[1] << " " << color[2] << endl;
          break;
        } else if (itrial == 0) {
          color[0] = 255;
          int w = 20;
          // if (256 - w < ix && ix < 256 + w && 256 - w < iy && iy < 256 + w)
          // cout << "ambiguity " << ix << " " << iy << endl;

          Parameter x = Parameter::constant((xs[ix] + xs[ix + 1]) / 2);
          Parameter y = Parameter::constant((xs[iy] + xs[iy + 1]) / 2);
          Parameter z = root;
          PV3<Parameter> u = poly3.gradient(x, y, z).unit();
          double uz = u.z.mid();
          static double maxUZ;
          if (fabs(uz) > maxUZ) {
            cout << "maxUZ " << uz << " ix " << ix << " iy " << iy << endl;
            maxUZ = fabs(uz);
          }
          out << color[0] << " " << color[1] << " " << color[2] << endl;
          break;
        }
      }
    }
  }

  out.close();
}

void probeZrestart(const char* filename, int nsteps, Poly3<Parameter>& poly3) {
  ofstream out;
  out.open(filename);
  out << "P3" << endl;
  out << nsteps << " " << nsteps << endl;
  out << 255 << endl;

  cout << "probeZ " << filename << " " << nsteps << endl;
  PTR<Object<Poly3>> opoly3 = new Object<Poly3>(poly3);
  vector<double> xs(nsteps + 1);
  for (int ix = 0; ix <= nsteps; ++ix) xs[ix] = -1.0 + 2.0 * ix / nsteps;

  P3 p3(poly3);

  for (int ix = 0; ix < nsteps; ++ix) {
    cout << "ix " << ix << endl;
    P3 p3x(p3);
    p3x.mobiusX(xs[ix], xs[ix + 1]);
    for (int iy = 0; iy < nsteps; ++iy) {
      // cout << "iy " << iy << endl;
      P3 p3xy(p3x);
      p3xy.mobiusY(xs[iy], xs[iy + 1]);

      // if (ix == 256 && iy == 256)
      // p3xy.print();

      for (int itrial = 0; itrial < 1 /*3*/; itrial++) {
        double zl = randomNumber(-1.1, -1);
        double zu = randomNumber(1, 1.1);
        vector<pair<double, double>> ints;
        ints.push_back(make_pair(zl, zu));

        Parameter root;

        while (ints.size() > 0) {
          pair<double, double> zz = ints.back();
          ints.pop_back();
          if (zz.first > 1.0 || zz.second < -1.0) continue;
          P3 p3xyz(p3xy);
          p3xyz.mobiusZ(zz.first, zz.second);
          int d = p3xyz.descartesZ();
          if (ix == 327 && iy == 499)
            cout << "zz " << zz.first << " " << zz.second << " d " << d << endl;
          if (d == 0) continue;
          if (d == 1) {
            PTR<Object<Scalar>> xyz[3];
            xyz[0] = new Object<Scalar>(
                Parameter::constant((xs[ix] + xs[ix + 1]) / 2));
            xyz[1] = new Object<Scalar>(
                Parameter::constant((xs[iy] + xs[iy + 1]) / 2));
            xyz[2] = xyz[0];
            PTR<Object<Poly>> poly = new Sub2Poly(opoly3, 2, xyz);
            PTR<Object<Scalar>> l =
                new Object<Scalar>(Parameter::constant(zz.first));
            PTR<Object<Scalar>> u =
                new Object<Scalar>(Parameter::constant(zz.second));
            vector<PTR<Object<Scalar>>> roots = getRoots(poly, l, u);
            assert(roots.size() == 1);
            Parameter r = roots[0]->getApprox();
            if (r.ub() < -1.0 || r.lb() > 1.0) continue;
            root = r;
            break;
          }
          if (zz.second - zz.first < 1.0 / 1024.0) {
            root = Parameter::constant(zz.first);
            // p3xyz.print();
            break;
          }
          double m = (zz.first + zz.second) / 2;
          ints.push_back(make_pair(m, zz.second));
          ints.push_back(make_pair(zz.first, m));
        }

        int color[3] = {0, 0, 0};
        if (root.uninitialized()) {
          out << color[0] << " " << color[1] << " " << color[2] << endl;
          break;
        } else if (root.lb() < root.ub()) {
          Parameter x = Parameter::constant((xs[ix] + xs[ix + 1]) / 2);
          Parameter y = Parameter::constant((xs[iy] + xs[iy + 1]) / 2);
          Parameter z = root;
          PV3<Parameter> u = poly3.gradient(x, y, z).unit();
          double uz = u.z.mid();
          int c = (-uz + 1) * 255 / 2;
          assert(c <= 255 && c >= 0);
          color[0] = color[1] = color[2] = c;
          out << color[0] << " " << color[1] << " " << color[2] << endl;
          break;
        } else if (itrial == 0 /*2*/) {
          color[0] = 255;
          int w = 20;
          // if (256 - w < ix && ix < 256 + w && 256 - w < iy && iy < 256 + w)
          // cout << "ambiguity " << ix << " " << iy << endl;

          Parameter x = Parameter::constant((xs[ix] + xs[ix + 1]) / 2);
          Parameter y = Parameter::constant((xs[iy] + xs[iy + 1]) / 2);
          Parameter z = root;
          PV3<Parameter> u = poly3.gradient(x, y, z).unit();
          double uz = u.z.mid();
          static double maxUZ;
          if (fabs(uz) > maxUZ) {
            cout << "maxUZ " << uz << " ix " << ix << " iy " << iy << endl;
            maxUZ = fabs(uz);
          }
          out << color[0] << " " << color[1] << " " << color[2] << endl;
          break;
        }
      }
    }
  }

  out.close();
}

int descartes(Poly<Parameter>& p, double l, double u) {
  Poly<Parameter> q = p.moebius(Parameter::constant(l), Parameter::constant(u));
  int count = 0;
  int s = q.a[0].sign(false);
  if (s == 0) return -1;
  for (int i = 1; i < q.a.size(); i++) {
    int si = q.a[i].sign(false);
    if (si == 0) return -1;
    if (si == -s) {
      s = -s;
      count++;
    }
  }
  return count;
}

void probeZi(const char* filename, int nsteps, Poly3<Parameter>& poly3) {
  ofstream out;
  out.open(filename);
  out << "P3" << endl;
  out << nsteps << " " << nsteps << endl;
  out << 255 << endl;

  cout << "probeZ " << filename << " " << nsteps << endl;
  PTR<Object<Poly3>> opoly3 = new Object<Poly3>(poly3);
  vector<double> xs(nsteps + 1);
  for (int ix = 0; ix <= nsteps; ++ix) xs[ix] = -1.0 + 2.0 * ix / nsteps;

  for (int ix = 0; ix < nsteps; ++ix) {
    cout << "ix " << ix << endl;
    Poly2<Parameter> poly2 =
        poly3.substitute1(0, Parameter::interval(xs[ix], xs[ix + 1]));
    // Poly2<Parameter> poly2 = poly3.substitute1(0,
    // Parameter::constant(xs[ix]));
    for (int iy = 0; iy < nsteps; ++iy) {
      // cout << "iy " << iy << endl;
      Poly<Parameter> poly =
          poly2.subX(Parameter::interval(xs[iy], xs[iy + 1]));
      // Poly<Parameter> poly = poly2.subX(Parameter::constant(xs[iy]));

      if (ix == 256 && iy == 256) {
        cout << "poly2";
        for (int i = 0; i < poly2.a.size(); ++i)
          cout << " " << poly2.m[2 * i] << "," << poly2.m[2 * i + 1] << ","
               << poly2.a[i].lb() << "," << poly2.a[i].ub();
        cout << endl;
        cout << "poly";
        for (int i = 0; i < poly.a.size(); ++i)
          cout << " " << poly.a[i].lb() << "," << poly.a[i].ub();
        cout << endl;
      }

      for (int itrial = 0; itrial < 3; itrial++) {
        double zl = randomNumber(-1.1, -1);
        double zu = randomNumber(1, 1.1);
        vector<pair<double, double>> ints;
        ints.push_back(make_pair(zl, zu));

        Parameter root;

        while (ints.size() > 0) {
          pair<double, double> zz = ints.back();
          ints.pop_back();
          if (zz.first > 1.0 || zz.second < -1.0) continue;
          int d = descartes(poly, zz.first, zz.second);
          if (ix == 327 && iy == 499)
            cout << "zz " << zz.first << " " << zz.second << " d " << d << endl;
          if (d == 0) continue;
          if (d == 1) {
            PTR<Object<Scalar>> xyz[3];
            xyz[0] = new Object<Scalar>(
                Parameter::constant((xs[ix] + xs[ix + 1]) / 2));
            xyz[1] = new Object<Scalar>(
                Parameter::constant((xs[iy] + xs[iy + 1]) / 2));
            xyz[2] = xyz[0];
            PTR<Object<Poly>> poly = new Sub2Poly(opoly3, 2, xyz);
            PTR<Object<Scalar>> l =
                new Object<Scalar>(Parameter::constant(zz.first));
            PTR<Object<Scalar>> u =
                new Object<Scalar>(Parameter::constant(zz.second));
            vector<PTR<Object<Scalar>>> roots = getRoots(poly, l, u);
            if (roots.size() != 1) {
              // cout << "roots.size() " << roots.size() << endl;
              continue;
            }
            assert(roots.size() == 1);
            Parameter r = roots[0]->getApprox();
            if (r.ub() < -1.0 || r.lb() > 1.0) continue;
            root = r;
            break;
          }
          if (zz.second - zz.first < 1.0 / 1024.0) {
            root = Parameter::constant(zz.first);
            // p3xyz.print();
            break;
          }
          double m = (zz.first + zz.second) / 2;
          ints.push_back(make_pair(m, zz.second));
          ints.push_back(make_pair(zz.first, m));
        }

        int color[3] = {0, 0, 0};
        if (root.uninitialized()) {
          out << color[0] << " " << color[1] << " " << color[2] << endl;
          break;
        } else if (root.lb() < root.ub()) {
          Parameter x = Parameter::constant((xs[ix] + xs[ix + 1]) / 2);
          Parameter y = Parameter::constant((xs[iy] + xs[iy + 1]) / 2);
          Parameter z = root;
          PV3<Parameter> u = poly3.gradient(x, y, z).unit();
          double uz = u.z.mid();
          int c = (-uz + 1) * 255 / 2;
          assert(c <= 255 && c >= 0);
          color[0] = color[1] = color[2] = c;
          out << color[0] << " " << color[1] << " " << color[2] << endl;
          break;
        } else if (itrial == 2) {
          color[0] = 255;
          int w = 20;
          // if (256 - w < ix && ix < 256 + w && 256 - w < iy && iy < 256 + w)
          // cout << "ambiguity " << ix << " " << iy << endl;

          Parameter x = Parameter::constant((xs[ix] + xs[ix + 1]) / 2);
          Parameter y = Parameter::constant((xs[iy] + xs[iy + 1]) / 2);
          Parameter z = root;
          PV3<Parameter> u = poly3.gradient(x, y, z).unit();
          double uz = u.z.mid();
          static double maxUZ;
          if (fabs(uz) > maxUZ) {
            cout << "maxUZ " << uz << " ix " << ix << " iy " << iy << endl;
            maxUZ = fabs(uz);
          }
          out << color[0] << " " << color[1] << " " << color[2] << endl;
          break;
        }
      }
    }
  }

  out.close();
}

int main(int argc, char* argv[]) {
  Poly3<Parameter> poly3;

  if (argc > 1) {
    cout << "reading " << argv[1] << endl;
    ifstream in;
    in.open(argv[1]);
    int n;
    in >> n;
    cout << n << " terms" << endl;
    for (int i = 0; i < n; i++) {
      int ix, iy, iz;
      double c;
      in >> ix >> iy >> iz >> c;
      cout << ix << " " << iy << " " << iz << " " << c << endl;
      poly3.add(ix, iy, iz, Parameter::constant(c));
    }
  } else {
    // poly3.add(1, 0, 0, Parameter::constant(-2.0));
    poly3.add(0, 0, 0, Parameter::constant(-1.0));
    poly3.add(2, 0, 0, Parameter::constant(1.0));
    poly3.add(0, 2, 0, Parameter::constant(1.0));
    poly3.add(0, 0, 2, Parameter::constant(1.0));
  }

  enable();
  // probeZrestart("probe.ppm", 512, poly3);
  probeZ("probe.ppm", 512, poly3);
  // probeZ("probe.ppm", 256, poly3);
  disable();
}
