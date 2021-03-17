#include "acp/poly/poly2d.h"

static PV2<Parameter> p(Parameter::input(1.0 / 3), Parameter::input(-2.0 / 3));

static void solve_test() {
  double d[2][3] = {{1, 2, 3}, {4, 5, 6}};
  Parameter a[2][3];
  Parameter b[2][3];
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      b[i][j] = a[i][j] = Parameter::constant(randomNumber(-1, 1));
  // b[i][j] = a[i][j] = Parameter::constant(d[i][j]);
  PV2<Parameter> x = solve(a);
  for (int i = 0; i < 2; ++i) cout << x[i].mid() << " ";
  cout << endl;
  for (int i = 0; i < 2; ++i) {
    Parameter v = b[i][2];
    for (int j = 0; j < 2; j++) v = v + b[i][j] * x[j];
    cout << v.lb() << " " << v.ub() << " " << v.ub() - v.lb() << endl;
    assert(v.sign(false) == 0);
  }
}

Poly2D<Parameter> randomPoly(int degx, int degy, int deg) {
  Poly2D<Parameter> poly(degx, degy);
  for (int ix = 0; ix <= degx; ix++)
    for (int iy = 0; iy <= degy; iy++)
      if (ix + iy <= deg)
        poly.set(ix, iy, Parameter::constant(randomNumber(-1, 1)));
  return poly;
}

void newton_test() {
  int deg = 5;
  Poly2D<Parameter> f = randomPoly(deg, deg, deg);
  Poly2D<Parameter> g = randomPoly(deg, deg, deg);

  f.add(0, 0, -f.value(p));
  g.add(0, 0, -g.value(p));

  assert(f.value(p).sign(false) == 0);
  assert(g.value(p).sign(false) == 0);

  double eps = 1e-3;
  Parameter eps_int =
      Parameter::constant(-eps).interval(Parameter::constant(2 * eps));
  PV2<Parameter> r(p.x + eps_int, p.y + eps_int);
  for (int i = 0; i < 2; i++)
    cout << r[i].lb() << " " << r[i].ub() << " " << r[i].ub() - r[i].lb()
         << endl;

  PV2<Parameter> rPrev;
  do {
    rPrev = r;
    r = newton<Parameter>(f, g, r);
    for (int i = 0; i < 2; i++)
      cout << r[i].lb() << " " << r[i].ub() << " " << r[i].ub() - r[i].lb()
           << endl;
  } while (r.x.subset(rPrev.x) || r.y.subset(rPrev.y));
}

void descartes_test() {
  int deg = 5;
  Poly2D<Parameter> fp = randomPoly(deg, deg, deg);
  fp[0] = (fp[0] - fp.value(p)).midP();
  Parameter v = fp.value(p);
  cout << "f(p) " << v.lb() << " " << v.ub() << " " << v.ub() - v.lb() << endl;
  PTR<Object<Poly2D>> f = new Object<Poly2D>(fp);

  double eps = 1e-3;
  Parameter eps_int =
      Parameter::constant(-eps).interval(Parameter::constant(2 * eps));
  PV2<Parameter> b1(p.x + eps_int, p.y + eps_int);
  int s1 = Descartes2D(f, b1);
  cout << "s1 " << s1 << endl;

  p = -p;
  PV2<Parameter> b2(p.x + eps_int, p.y + eps_int);
  int s2 = Descartes2D(f, b2);
  cout << "s2 " << s2 << endl;
}

void boundaryRoots_test() {
  int deg = 5;
  Poly2D<Parameter> fp = randomPoly(deg, deg, deg);
  fp[0] = (fp[0] - fp.value(p)).midP();
  PTR<Object<Poly2D>> f = new Object<Poly2D>(fp);

  double eps = 1e-3;
  Parameter eps_int =
      Parameter::constant(-eps).interval(Parameter::constant(2 * eps));
  PV2<Parameter> b1(p.x + eps_int, p.y + eps_int);
  vector<PTR<Object<PV2>>> roots1 = boundaryRoots(f, b1);
  cout << "roots1.size() " << roots1.size() << endl;

  p = -p;
  PV2<Parameter> b2(p.x + eps_int, p.y + eps_int);
  vector<PTR<Object<PV2>>> roots2 = boundaryRoots(f, b2);
  cout << "roots2.size() " << roots2.size() << endl;
}

void getRoots_test() {
  int deg = 5;
  Poly2D<Parameter> fp = randomPoly(deg, deg, deg);
  fp[0] = (fp[0] - fp.value(p)).midP();
  PTR<Object<Poly2D>> f = new Object<Poly2D>(fp);
  Poly2D<Parameter> gp = randomPoly(deg, deg, deg);
  gp[0] = (gp[0] - gp.value(p)).midP();
  PTR<Object<Poly2D>> g = new Object<Poly2D>(gp);

  Parameter m1p1 = Parameter::constant(-1).interval(Parameter::constant(1));
  PV2<Parameter> box(m1p1, m1p1);

  getRoots(f, g, box);
}

int main(int argc, char* argv[]) {
  enable();
  // solve_test();
  // newton_test();
  // descartes_test();
  // boundaryRoots_test();
  getRoots_test();
  return 0;
}
