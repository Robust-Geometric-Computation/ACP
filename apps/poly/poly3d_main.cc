#include "acp/encasement2d/encasement2d.h"
#include "acp/poly/poly2d.h"
#include "acp/poly/poly3d.h"

using namespace std;
using namespace acp;

static PV3<Parameter> p(Parameter::input(1.0 / 3), Parameter::input(-1.0 / 3),
                        Parameter::input(2.0 / 3));

void solve_test() {
  double d[3][4] = {{0, 5, 6, 7}, {1, 2, 11, 13}, {2, 4, 6, 8}};
  Parameter a[3][4];
  Parameter b[3][4];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 4; ++j)
      b[i][j] = a[i][j] = Parameter::constant(randomNumber(-1, 1));
  // b[i][j] = a[i][j] = Parameter::constant(d[i][j]);
  PV3<Parameter> x = solve(a);
  for (int i = 0; i < 3; ++i) cout << x[i].mid() << " ";
  cout << endl;
  for (int i = 0; i < 3; ++i) {
    Parameter v = b[i][3];
    for (int j = 0; j < 3; j++) v = v + b[i][j] * x[j];
    cout << v.lb() << " " << v.ub() << " " << v.ub() - v.lb() << endl;
    assert(v.sign(false) == 0);
  }
}

Poly3D<Parameter> randomPoly(int degx, int degy, int degz, int deg) {
  Poly3D<Parameter> poly(degx, degy, degz);
  for (int ix = 0; ix <= degx; ix++)
    for (int iy = 0; iy <= degy; iy++)
      for (int iz = 0; iz <= degz; iz++)
        if (ix + iy + iz <= deg)
          poly.set(ix, iy, iz, Parameter::constant(randomNumber(-1, 1)));
  return poly;
}

void newton_test() {
  int deg = 5;
  Poly3D<Parameter> f = randomPoly(deg, deg, deg, deg);
  Poly3D<Parameter> g = randomPoly(deg, deg, deg, deg);
  Poly3D<Parameter> h = randomPoly(deg, deg, deg, deg);

  f.add(0, 0, 0, -f.value(p));
  g.add(0, 0, 0, -g.value(p));
  h.add(0, 0, 0, -h.value(p));

  assert(f.value(p).sign(false) == 0);
  assert(g.value(p).sign(false) == 0);
  assert(h.value(p).sign(false) == 0);

  double eps = 1e-3;
  Parameter eps_int =
      Parameter::constant(-eps).interval(Parameter::constant(2 * eps));
  PV3<Parameter> r(p.x + eps_int, p.y + eps_int, p.z + eps_int);
  for (int i = 0; i < 3; i++)
    cout << r[i].lb() << " " << r[i].ub() << " " << r[i].ub() - r[i].lb()
         << endl;

  PV3<Parameter> rPrev;
  do {
    rPrev = r;
    r = newton<Parameter>(f, g, h, r);
    for (int i = 0; i < 3; i++)
      cout << r[i].lb() << " " << r[i].ub() << " " << r[i].ub() - r[i].lb()
           << endl;
  } while (r.x.subset(rPrev.x) || r.y.subset(rPrev.y) || r.z.subset(rPrev.z));
}

void descartes_test() {
  int deg = 5;
  Poly3D<Parameter> fp = randomPoly(deg, deg, deg, deg);
  fp[0] = (fp[0] - fp.value(p)).midP();
  Parameter v = fp.value(p);
  cout << "f(p) " << v.lb() << " " << v.ub() << " " << v.ub() - v.lb() << endl;
  PTR<Object<Poly3D>> f = new Object<Poly3D>(fp);

  double eps = 1e-3;
  Parameter eps_int =
      Parameter::constant(-eps).interval(Parameter::constant(2 * eps));
  PV3<Parameter> b1(p.x + eps_int, p.y + eps_int, p.z + eps_int);
  int s1 = Descartes3D(f, b1);
  cout << "s1 " << s1 << endl;

  p = -p;
  PV3<Parameter> b2(p.x + eps_int, p.y + eps_int, p.z + eps_int);
  int s2 = Descartes3D(f, b2);
  cout << "s2 " << s2 << endl;
}

void boundaryRoots_test() {
  int deg = 2;
  Poly3D<Parameter> fp = randomPoly(deg, deg, deg, deg);
  fp[0] = (fp[0] - fp.value(p)).midP();
  PTR<Object<Poly3D>> f = new Object<Poly3D>(fp);
  Poly3D<Parameter> gp = randomPoly(deg, deg, deg, deg);
  gp[0] = (gp[0] - gp.value(p)).midP();
  PTR<Object<Poly3D>> g = new Object<Poly3D>(gp);

  double eps = 1e-3;
  Parameter eps_int =
      Parameter::constant(-eps).interval(Parameter::constant(2 * eps));
  PV3<Parameter> b1(p.x + eps_int, p.y + eps_int, p.z + eps_int);
  vector<PTR<Object<PV3>>> roots1 = boundaryRoots(f, g, b1);
  cout << "roots1.size() " << roots1.size() << endl;

  p = -p;
  PV3<Parameter> b2(p.x + eps_int, p.y + eps_int, p.z + eps_int);
  vector<PTR<Object<PV3>>> roots2 = boundaryRoots(f, g, b2);
  cout << "roots2.size() " << roots2.size() << endl;
}

void getRoots_test() {
  int deg = 10;
  Poly3D<Parameter> fp = randomPoly(deg, deg, deg, deg);
  fp[0] = (fp[0] - fp.value(p)).midP();
  PTR<Object<Poly3D>> f = new Object<Poly3D>(fp);
  Poly3D<Parameter> gp = randomPoly(deg, deg, deg, deg);
  gp[0] = (gp[0] - gp.value(p)).midP();
  PTR<Object<Poly3D>> g = new Object<Poly3D>(gp);
  Poly3D<Parameter> hp = randomPoly(deg, deg, deg, deg);
  hp[0] = (hp[0] - hp.value(p)).midP();
  PTR<Object<Poly3D>> h = new Object<Poly3D>(hp);

  Parameter m1p1 = Parameter::constant(-1).interval(Parameter::constant(1));
  PV3<Parameter> box(m1p1, m1p1, m1p1);

  getRoots(f, g, h, box);
}

int main(int argc, char* argv[]) {
  acp::enable();
  // solve_test();
  // newton_test();
  // descartes_test();
  // boundaryRoots_test();
  getRoots_test();
}
