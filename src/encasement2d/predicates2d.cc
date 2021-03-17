#include "acp/encasement2d/predicates2d.h"
#include <chrono>

template <class N>
N Line::Left::calculate() {
  return l->leftNeg<N>(q);
}

template Parameter Line::Left::calculate<Parameter>();
template PParameter Line::Left::calculate<PParameter>();
template MParameter Line::Left::calculate<MParameter>();

int Ellipse::count = 0;

vector<PTR<Object<Poly>>> NewtonPoly::bas;

bool LinePoly::print_all = false;

#ifdef BLEEN
PTR<Object<PV2>> XbetterApprox(PV2& p) {
  double lx = p.x.lb(), ly = p.y.lb();
  double ux = p.x.ub(), uy = p.y.ub();
  PV2 l = PV2::constant(lx, ly);
  PV2 r = PV2::constant((ux - lx) * randomNumber(0, 1),
                        (uy - ly) * randomNumber(0, 1));
  return new SumPoint(new InputPoint(l), new InputPoint(r));
}

PTR<Object<PV2>> betterApprox(PV2 p) {
  vector<PV2> lbs;
  while ((p.x.ub() - p.x.lb()) < 1e-8 * fabs(p.x.lb()) ||
         (p.y.ub() - p.y.lb()) < 1e-8 * fabs(p.y.lb())) {
    PV2 lb = PV2::constant(p.x.lb(), p.y.lb());
    lbs.push_back(lb);
    p = p - lb;
  }

  PV2 r = PV2::constant(randomNumber(p.x.lb(), p.x.ub()),
                        randomNumber(p.y.lb(), p.y.ub()));

  PTR<Object<PV2>> sum = new InputPoint(r);
  while (lbs.size() > 0) {
    sum = new SumPoint(new InputPoint(lbs.back()), sum);
    lbs.pop_back();
  }

  return sum;
}
#endif

class LineApproximator : public Primitive {
  PTR<Line> l;
  PTR<Object<PV2>> a, b;
  vector<double> ps;
  vector<double> vs;

  DeclareSign {
    PV2<N> p = l->getP<N>();
    PV2<N> v = l->getV<N>();
    bool between = ((a->get<N>() - p).cross(v).sign() !=
                    (b->get<N>() - p).cross(v).sign());

    while ((p.x.ub() - p.x.lb()) < 1e-8 * fabs(p.x.lb()) ||
           (p.y.ub() - p.y.lb()) < 1e-8 * fabs(p.y.lb())) {
      double x = p.x.lb(), y = p.y.lb();
      ps.push_back(y);
      ps.push_back(x);
      p = p - PV2<N>(N::constant(x), N::constant(y));
    }
    ps.push_back(randomNumber(p.y.lb(), p.y.ub()));
    ps.push_back(randomNumber(p.x.lb(), p.x.ub()));

    while ((v.x.ub() - v.x.lb()) < 1e-8 * fabs(v.x.lb()) ||
           (v.y.ub() - v.y.lb()) < 1e-8 * fabs(v.y.lb())) {
      double x = v.x.lb(), y = v.y.lb();
      vs.push_back(y);
      vs.push_back(x);
      v = v - PV2<N>(N::constant(x), N::constant(y));
    }
    vs.push_back(randomNumber(v.y.lb(), v.y.ub()));
    vs.push_back(randomNumber(v.x.lb(), v.x.ub()));

    return between ? N(-1) : N(1);
  }

 public:
  LineApproximator(PTR<Line> l, PTR<Object<PV2>> a, PTR<Object<PV2>> b)
      : l(l), a(a), b(b) {}

  PTR<Line> getApprox() {
    int sign = operator int();
    double x, y;

    x = ps.back();
    ps.pop_back();
    y = ps.back();
    ps.pop_back();
    PTR<Object<PV2>> p = new InputPoint(
        PV2<Parameter>(Parameter::constant(x), Parameter::constant(y)));
    while (ps.size() > 0) {
      x = ps.back();
      ps.pop_back();
      y = ps.back();
      ps.pop_back();

      p = new SumPoint(new InputPoint(PV2<Parameter>(Parameter::constant(x),
                                                     Parameter::constant(y))),
                       p);
    }

    x = vs.back();
    vs.pop_back();
    y = vs.back();
    vs.pop_back();
    PTR<Object<PV2>> v = new InputPoint(
        PV2<Parameter>(Parameter::constant(x), Parameter::constant(y)));
    while (vs.size() > 0) {
      x = vs.back();
      vs.pop_back();
      y = vs.back();
      vs.pop_back();

      v = new SumPoint(new InputPoint(PV2<Parameter>(Parameter::constant(x),
                                                     Parameter::constant(y))),
                       v);
    }

    return new Line(p, v);
  }
};

PTR<Line> Line::approximate(PTR<Object<PV2>> a, PTR<Object<PV2>> b) {
  assert(left(a) != left(b));

  PV2<Parameter> pp = p->getApprox(1.0);
  PV2<Parameter> vv = v->getApprox(1.0);

  PTR<Object<PV2>> pa = new InputPoint(pp.x.mid(), pp.y.mid());
  PTR<Object<PV2>> va = new InputPoint(vv.x.mid(), vv.y.mid());

  PTR<Line> l = new Line(pa, va);

  if (!(l->left(a) != l->left(b))) {
    l = LineApproximator(this, a, b).getApprox();
    assert(l->left(a) != l->left(b));
  }

  return l;
}

double dummy(double x) {
  x = x * x;
  return x;
}

template <class N>
Poly<N> substitute(PTR<Line> l, PTR<Object<Poly2>> f) {
  static int count;
  // substitute x = f(t) and y = g(t) into F(x, y) to get univariate F(t)
  //  if (++count == 2646739) {
  //    cout << "this is it" << endl;
  //    dummy(17);
  //  }

  Poly2<N> curve = f->get<N>();

  Poly<N> x_t;
  x_t.add(l->getP<N>().x, 0);
  // if vertical line, no linear term in x
  if (l->getAxis() != 1) x_t.add(l->getV<N>().x, 1);

  Poly<N> y_t;
  y_t.add(l->getP<N>().y, 0);
  // if horizontal line, no linear term in y
  if (l->getAxis() != 0) y_t.add(l->getV<N>().y, 1);

  vector<Poly<N>> xs;
  vector<Poly<N>> ys;

  xs.push_back(Poly<N>());
  xs.push_back(x_t);
  for (int i = 2; i <= curve.degx; i++) {
    xs.push_back(xs[xs.size() - 1] * x_t);
  }

  ys.push_back(Poly<N>());
  ys.push_back(y_t);
  for (int i = 2; i <= curve.degy; i++) {
    ys.push_back(ys[ys.size() - 1] * y_t);
  }

  // assert(coef_x[0].increased() == coef_y[0].increased());
  /*
    int degree = curve.d;

    vector<Parameter> coef(degree + 1);

    for(int i = 0; i <= degree; i++) {
      coef[i] = Parameter::constant(0);
    }

    for(int k = 0; k < curve.size(); k++) {

      if(curve.m[2*k] > 0 && curve.m[2*k+1] > 0) {

        int x = curve.m[2*k];
        int y = curve.m[2*k+1];

        if(xs[x].d == 0 && ys[y].d == 0) {

          coef[0] = coef[0] + (curve.a[k] * xs[x].a[0] * ys[y].a[0]);

        } else if(xs[x].d == 0) {

          for(int j = 0; j <= y; j++) {
            coef[j] = coef[j] + (curve.a[k] * xs[x].a[0] * ys[y].a[j]);
          }

        } else if(ys[y].d == 0) {

          for(int i = 0; i <= x; i++) {
            coef[i] = coef[i] + (curve.a[k] * xs[x].a[i] * ys[y].a[0]);
          }

        } else {

          for(int i = 0; i <= x; i++) {
            for(int j = 0; j <= y; j++) {
              coef[i+j] = coef[i+j] + (curve.a[k] * xs[x].a[i] * ys[y].a[j]);
            }
          }
        }

      } else if(curve.m[2*k] > 0) {

        int x = curve.m[2*k];

        if(xs[x].d == 0) {

          coef[0] = coef[0] + (curve.a[k] * xs[x].a[0]);

        } else {

          for(int i = 0; i <= x; i++) {
            coef[i] = coef[i] + (curve.a[k] * xs[x].a[i]);
          }

        }

      } else if(curve.m[2*k+1] > 0) {

        int y = curve.m[2*k+1];

        if(ys[y].d == 0) {

          coef[0] = coef[0] + (curve.a[k] * ys[y].a[0]);

        } else {

          for(int i = 0; i <= y; i++) {
            coef[i] = coef[i] + (curve.a[k] * ys[y].a[i]);
          }

        }

      } else {
        coef[0] = coef[0] + curve.a[k];
      }

    }

    return Poly(degree, &coef[0]);
  */

  // construct the base polynomial (starts as constant 0)
  N base[1];
  base[0] = N::constant(0);
  // TODO why did we do this? Bug in the Parameter::constant constructor
  // base[0].increasePrecision();
  // Poly<N> poly(base);
  Poly<N> poly(N(0));

  for (int i = 0; i < curve.size(); i++) {
    // temp polynomial to multiply against
    N temp[1];
    temp[0] = N::constant(1);
    // TODO why did we do this? Bug in the Parameter::constant constructor
    // temp[0].increasePrecision();
    // Poly<N> temp_poly(0, temp);
    Poly<N> temp_poly(N(1));

    // substitue x(t)
    // for(int j = 0; j < curve.m[2*i]; j++) {
    //  temp_poly = temp_poly * x_t;
    //}

    // substitute y(t)
    // for(int j = 0; j < curve.m[2*i + 1]; j++) {
    //  temp_poly = temp_poly * y_t;
    //}

    if (curve.m[2 * i] > 0) temp_poly = temp_poly * xs[curve.m[2 * i]];
    if (curve.m[2 * i + 1] > 0) temp_poly = temp_poly * ys[curve.m[2 * i + 1]];

    poly = poly.plusCtimes(curve.a[i], temp_poly);
  }

  return poly;
}

template <class N>
PV2<N> EllipseCriticalPoint::calculate() {
  Poly2<N> dx = e->get<N>().derX();
  Poly2<N> dy = e->get<N>().derY();

  N a, b, c, d, e, f;
  a = dx.a[0];
  b = dx.a[1];
  c = dx.a[3];

  d = dy.a[1];
  e = dy.a[2];
  f = dy.a[4];

  // Cramer's Rule for 2x2
  return PV2<N>((b * f - c * e) / (a * e - b * d),
                (c * d - a * f) / (a * e - b * d));
}

double max_error = -1;
double do_error = -1;

template <class N>
Poly<N> LinePoly::calculate() {
  /*
  if(print_all) {

    do_error = 1;
    max_error = 0;


    cout << endl;
    cout << "------------------------------------------------" << endl;
    cout << endl;

    auto t1 = std::chrono::high_resolution_clock::now();

    Poly sub = substitute(l, f);

    auto t2 = std::chrono::high_resolution_clock::now();

    sub.normalized().print();

    cout << "\nsubs: " <<
  std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count() * 0.001
  << " ms     \n";



    cout << endl;

    t1 = std::chrono::high_resolution_clock::now();

    Poly in = interpolate(l, f);

    t2 = std::chrono::high_resolution_clock::now();

    in.normalized().print();

   cout << "\ninterpolation: " <<
  std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count() * 0.001
  << " ms     \n";


    cout << endl;
    cout << "------------------------------------------------" << endl;
    cout << endl;

    do_error = -1;
    max_error = -1;

    if(Parameter::highPrecision == 53u) {
      cout << "double precision" << endl;
      Parameter::constant(0).sign();
    } else {
      cout << "higher precision" << endl;
    }

    return sub;

  } else {
  */

  return substitute<N>(l, f);

  if (l->getAxis() > -1) {
    return substitute<N>(l, f);
  } else {
    return interpolate<N>(l, f);
  }
  /*
  }

    //if y0 + t * y1 is just a constant, old substitution should be more
  computationally efficient
    //if(l->getV().x == 0 || l->getV().y == 0) {
    //  return substitute(l, f);
    //} else {
    //  return interpolate(l, f);
    //}

  */
}

template <class N>
Poly<N> interpolate(PTR<Line> l, PTR<Object<Poly2>> f) {
  // do interpolation on t = [0, 1, 2, 3, ... f.d] to obtain coefficients in
  // newton poly basis

  Poly2<N> ff = f->get<N>();

  int number_terms = ff.d + 1;

  vector<N> t;  // t values
  vector<N> b;  // values along line, t = [0, 1, 2, 3, ... ]

  vector<PV2<N>> p;  // points along line, t = [0, 1, 2, 3, ... ]

  for (int i = 0; i < number_terms; i++) {
    t.push_back(NewtonPoly::newtonRoot<N>(i));
    p.push_back(l->getP<N>(t[i]));
    b.push_back(ff.value(p[i]));
  }

  return interpolate<N>(b);
}

template <class N>
Poly<N> interpolate(vector<N> b) {
  vector<N> a;
  vector<N> t;  // t values

  vector<Poly<N>> newton;

  // set up points and values as well as precompute newton basis polys
  for (int i = 0; i < b.size(); i++) {
    t.push_back(NewtonPoly::newtonRoot<N>(i));
    newton.push_back(NewtonPoly::basis<N>(i));
  }

  // starting with the lowest terms, build up the ai's
  for (int i = 0; i < b.size(); i++) {
    N ai;

    N numer = b[i];
    N sum_term = N::constant(0);

    N denom = N::constant(1);

    // linear combination of all of the previous a[j] with previous basis polys
    // (newton[j]) evaluated at this t value (t[i])
    for (int j = 1; j <= i; j++) {
      denom = denom * (t[i] - t[i - j]);
      sum_term = (sum_term * (t[i] - t[i - j])) + a[i - j];
    }

    ai = (numer - sum_term) / denom;

    a.push_back(ai);
  }

  if (a.size() == 0) return Poly<N>();

  Poly<N> result = newton[0] * a[0];

  for (int i = 1; i < a.size(); i++) {
    result = result + (newton[i] * a[i]);
  }

  return result;
}
template <class N>
Poly<N> LinearCombinationPoly::calculate() {
  assert(a.size() == p.size());

  if (a.size() == 0) return Poly<N>();

  Poly<N> f = p[0]->get<N>() * a[0]->get<N>().x;

  for (int i = 1; i < a.size(); i++) {
    f = f + (p[i]->get<N>() * a[i]->get<N>().x);
  }

  return f;
}
template <class N>
Poly<N> ShiftOutPoly::calculate() {
  Poly<N> f = p->get<N>();

  vector<N> coef(f.size() - 1);
  for (int i = 0; i < coef.size(); i++) {
    coef[i] = f.a[i + 1];
  }

  return Poly<N>(coef);
}
template <class N>
Poly<N> MiddlePoly::calculate() {
  PTR<Line> l = new LineAB(p, q);

  // get each of their ppoly1s
  Poly<N> F = substitute<N>(l, f);
  Poly<N> G = substitute<N>(l, g);

  // Need to normalize the curves here by their gradient
  Poly2<N> fx = f->get<N>().derX();
  Poly2<N> fy = f->get<N>().derY();
  Poly2<N> gx = g->get<N>().derX();
  Poly2<N> gy = g->get<N>().derY();

  N xy[] = {p->get<N>().x, p->get<N>().y};

  PV2<N> grad_f(fx.value(xy), fy.value(xy));
  PV2<N> grad_g(gx.value(xy), gy.value(xy));

  N fc = (1.0 / grad_f.length());
  N gc = (1.0 / grad_g.length());

  Poly<N> fscale(fc);
  Poly<N> gscale(gc);

  F = F * fscale;
  G = G * gscale;

  // take the convention that G(p) > 0
  // and F(q) > 0
  if (PolySign(g, p) < 0) {
    // printf("need to flip G\n");
    G = -G;
  }

  if (PolySign(f, q) < 0) {
    // printf("need to flip F\n");
    F = -F;
  }

  /*delete l*/;

  return F - G;
}

class Approximator : public Primitive {
  PTR<Object<Scalar>> a, b, aver;
  vector<double> v;

  DeclareSign {
    N p = aver->get<N>();
    bool between = (a->get<N>().x < p && p < b->get<N>().x);

    while ((p.ub() - p.lb()) < 1e-8 * fabs(p.lb())) {
      double lb = p.lb();
      v.push_back(lb);
      p = p - lb;
    }
    v.push_back(randomNumber(p.lb(), p.ub()));

    return between ? N(-1) : N(1);
  }

 public:
  Approximator(PTR<Object<Scalar>> a, PTR<Object<Scalar>> b,
               PTR<Object<Scalar>> aver)
      : a(a), b(b), aver(aver) {}

  PTR<Object<Scalar>> getApprox() {
    int sign = operator int();
    // assert(sign == -1);
    PTR<Object<Scalar>> sum = new InputScalar(Parameter::constant(v.back()));
    v.pop_back();
    while (v.size() > 0) {
      sum = new SumAB(new InputScalar(Parameter::constant(v.back())), sum);
      v.pop_back();
    }
    return sum;
  }
};

PTR<Object<Scalar>> approxAverage(double t, PTR<Object<Scalar>> a,
                                  PTR<Object<Scalar>> b) {
  PTR<Object<Scalar>> aver = new AverageAB(t, a, b);
  assert((LessThan(a, aver) < 0 && LessThan(aver, b) < 0) ||
         (LessThan(b, aver) < 0 && LessThan(aver, a) < 0));
  double d = aver->getApprox(1.0).x.mid();
  PTR<Object<Scalar>> o = new InputScalar(Parameter::input(d));
  if (!((LessThan(a, o) < 0 && LessThan(o, b) < 0) ||
        (LessThan(b, o) < 0 && LessThan(o, a) < 0))) {
    /*
    double l = aver->get().lb();
    double u = aver->get().ub();
    assert(l < u);
    double r = (u - l) * randomNumber(0, 1);

    Parameter lp = Parameter::constant(l);
    Parameter rp = Parameter::constant(r);

    PTR<Object<Scalar>> lo = new InputScalar(lp);
    PTR<Object<Scalar>> ro = new InputScalar(rp);

    o = new SumAB(lo, ro);
    */

    o = Approximator(a, b, aver).getApprox();

    assert((LessThan(a, o) < 0 && LessThan(o, b) < 0) ||
           (LessThan(b, o) < 0 && LessThan(o, a) < 0));
  }
  return o;
}

PTR<Object<Scalar>> middleT(PTR<Object<PV2>> p, PTR<Object<Poly2>> f,
                            PTR<Object<PV2>> q, PTR<Object<Poly2>> g) {
  PTR<Object<Poly>> mp = new MiddlePoly(f, g, p, q);

  PTR<Object<Scalar>> zero = new InputScalar(Parameter::constant(0));
  PTR<Object<Scalar>> one = new InputScalar(Parameter::constant(1));

  std::vector<PTR<Object<Scalar>>> roots = getRoots(mp, zero, one);

  assert(roots.size() > 0);

  // Ascending asc;

  // return the middle root of the poly f(p + t*(q-p)) - g(p + t*(q-p)) = 0
  // std::sort(roots.begin(), roots.end(), asc);

  return roots[roots.size() / 2];
  /*
  PTR<Object<Scalar>> ret = new
  InputScalar(Parameter::input(root->get().mid())); if (!(LessThan(zero, ret) <
  0 && LessThan(ret, one) < 0)) { ret = Approximator(zero, one,
  root).getApprox(); assert(LessThan(zero, ret) < 0 && LessThan(ret, one) < 0);
  }
*/
  // assert(ret->get().sign() > 0 && (1 - ret->get()).sign() > 0);
  /*
  assert(Sign(ret) > 0 &&
         Sign(new DiffAB(
                  new InputScalar(Parameter::constant(1)), ret)) > 0);
  */
  /*
    return ret;
  */
}

PTR<Object<PV2>> middleV(PTR<Object<Poly2>> f, PTR<Object<Poly2>> g,
                         PTR<Object<PV2>> p, PTR<Object<PV2>> a,
                         PTR<Object<PV2>> b) {
  PTR<Object<PV2>> v = new NormalVector(new TangentVector(a, f));
  PTR<Object<PV2>> w = new NormalVector(new TangentVector(b, g));

  PTR<Object<PV2>> minusW = new NegativeVector(w);

  PTR<Object<PV2>> ab = new VectorAB(a, b);

  bool other = (Side(ab, v) == Side(ab, w));

  PTR<Object<PV2>> h = new HalfVector(v, other ? w : minusW);

  // PTR<Object<PV2>>  mv = new InputPoint(h->get().x.mid(), h->get().y.mid());

  assert(LeftOf(p, h, a) != LeftOf(p, h, b));

  return h;

  /*delete v*/;
  /*delete w*/;
  /*delete ab*/;
  /*delete minusW*/;
  /*delete h*/;

  // return mv;
}
template <class N>
Poly2<N> QuadApproxPoly2::calculate() {
  Poly2<N> ff = f->get<N>();

  if (i == 1) ff = ff * -1;

  Poly2<N> fx = ff.derX();
  Poly2<N> fy = ff.derY();
  Poly2<N> fxx = fx.derX();
  Poly2<N> fxy = fx.derY();
  Poly2<N> fyy = fy.derY();

  PV2<N> pp = p->get<N>();

  N a = pp.x.lbP();
  N b = pp.y.lbP();

  vector<N> coef;

  //

  // if(i == 0) {
  if (true) {
    coef.push_back(ff.value(a, b).lbP());
    coef.push_back(fx.value(a, b).lbP());
    coef.push_back(fy.value(a, b).lbP());
    coef.push_back(0.5 * fxx.value(pp.x, b).lbP());
    coef.push_back(fxy.value(pp.x, b).lbP());
    coef.push_back(0.5 * fyy.value(pp.x, pp.y).lbP());
  } else {
    coef.push_back(-(ff.value(a, b).ubP()));
    coef.push_back(-(fx.value(a, b).ubP()));
    coef.push_back(-(fy.value(a, b).ubP()));
    coef.push_back(-((0.5 * fxx.value(pp.x, b)).ubP()));
    coef.push_back(-(fxy.value(pp.x, b).ubP()));
    coef.push_back(-((0.5 * fyy.value(pp.x, pp.y)).ubP()));
  }

  //  if(i == 0) {
  //    coef.push_back(ff.value({a, b}).lbP());
  //    coef.push_back(fx.value({a, b}).lbP());
  //    coef.push_back(fy.value({a, b}).lbP());
  //    coef.push_back((0.5*fxx.value({pp.x, pp.y})).lbP());
  //    coef.push_back(fxy.value({pp.x, pp.y}).lbP());
  //    coef.push_back((0.5*fyy.value({pp.x, pp.y})).lbP());
  //  } else {
  //    coef.push_back(-ff.value({a, b}).ubP());
  //    coef.push_back(-fx.value({a, b}).ubP());
  //    coef.push_back(-fy.value({a, b}).ubP());
  //    coef.push_back((-0.5*fxx.value({pp.x, pp.y})).ubP());
  //    coef.push_back(-fxy.value({pp.x, pp.y}).ubP());
  //    coef.push_back((-0.5*fyy.value({pp.x, pp.y})).ubP());
  //  }

  int m[12] = {0, 0, 1, 0, 0, 1, 2, 0, 1, 1, 0, 2};

  return Poly2<N>(6, &coef[0], m);
}

template <class N>
PV2<N> QuadCriticalPoint::calculate() {
  Poly2<N> ff = f->get<N>();

  // f must be a quadratic with terms in this order
  assert(ff.size() == 6);
  assert(ff.m[0] == 0 && ff.m[1] == 0 && ff.m[2] == 1 && ff.m[3] == 0 &&
         ff.m[4] == 0 && ff.m[5] == 1 && ff.m[6] == 2 && ff.m[7] == 0 &&
         ff.m[8] == 1 && ff.m[9] == 1 && ff.m[10] == 0 && ff.m[11] == 2);

  N c1 = ff.a[1];
  N c2 = ff.a[2];

  N a1 = 2 * ff.a[3];
  N a2 = ff.a[4];
  N b1 = ff.a[4];
  N b2 = 2 * ff.a[5];

  N x = (b1 * c2 - c1 * b2) / (a1 * b2 - b1 * a2);
  N y = (c1 * a2 - a1 * c2) / (a1 * b2 - b1 * a2);

  return PV2<N>(x, y);
}
template <class N>
Poly2<N> CubicApproxPoly2::calculate() {
  Poly2<N> ff = f->get<N>();

  if (i == 1) ff = ff * -1;

  Poly2<N> fx = ff.derX();
  Poly2<N> fy = ff.derY();

  Poly2<N> fxx = fx.derX();
  Poly2<N> fxy = fx.derY();
  Poly2<N> fyy = fy.derY();

  Poly2<N> fxxx = fxx.derX();
  Poly2<N> fxxy = fxx.derY();
  Poly2<N> fxyy = fxy.derY();
  Poly2<N> fyyy = fyy.derY();

  PV2<N> pp = p->get<N>();

  N a = pp.x.lbP();
  N b = pp.y.lbP();

  vector<N> coef;

  // if(i == 0) {
  if (true) {
    coef.push_back(ff.value(a, b).lbP());

    coef.push_back(fx.value(a, b).lbP());
    coef.push_back(fy.value(a, b).lbP());

    coef.push_back((0.5 * fxx.value(a, b)).lbP());
    coef.push_back(fxy.value(a, b).lbP());
    coef.push_back((0.5 * fyy.value(a, b)).lbP());

    coef.push_back((fxxx.value(pp.x, b) / 6.0).lbP());
    coef.push_back(
        ((2 * fxxy.value(pp.x, b) + fxxy.value(pp.x, pp.y)) / 6.0).lbP());
    coef.push_back(
        ((2 * fxyy.value(pp.x, pp.y) + fxyy.value(pp.x, b)) / 6.0).lbP());
    coef.push_back((fyyy.value(pp.x, pp.y) / 6.0).lbP());

    // coef.push_back(((1.0/6.0)*fxxx.value(pp.x, pp.y)).lbP());
    // coef.push_back(((3.0/6.0)*(fxxy.value(pp.x, pp.y))).lbP());
    // coef.push_back(((3.0/6.0)*(fxyy.value(pp.x, pp.y))).lbP());
    // coef.push_back(((1.0/6.0)*fyyy.value(pp.x, pp.y)).lbP());

  } else {
    coef.push_back(-(ff.value(a, b).ubP()));

    coef.push_back(-(fx.value(a, b).ubP()));
    coef.push_back(-(fy.value(a, b).ubP()));

    coef.push_back(-((0.5 * fxx.value(a, b)).ubP()));
    coef.push_back(-((fxy.value(a, b)).ubP()));
    coef.push_back(-((0.5 * fyy.value(a, b)).ubP()));

    // coef.push_back(((-1.0/6.0)*fxxx.value(pp.x, b)).ubP());
    // coef.push_back(((-1.0/6.0)*(2*fxxy.value(pp.x, b) + fxxy.value(pp.x,
    // pp.y))).ubP()); coef.push_back(((-1.0/6.0)*(2*fxyy.value(pp.x, pp.y) +
    // fxyy.value(pp.x, b))).ubP()); coef.push_back(((-1.0/6.0)*fyyy.value(pp.x,
    // pp.y)).ubP());

    coef.push_back(-(((1.0 / 6.0) * fxxx.value(pp.x, pp.y)).ubP()));
    coef.push_back(-(((3.0 / 6.0) * (fxxy.value(pp.x, pp.y))).ubP()));
    coef.push_back(-(((3.0 / 6.0) * (fxyy.value(pp.x, pp.y))).ubP()));
    coef.push_back(-(((1.0 / 6.0) * fyyy.value(pp.x, pp.y)).ubP()));
  }

  int m[20] = {0, 0, 1, 0, 0, 1, 2, 0, 1, 1, 0, 2, 3, 0, 2, 1, 1, 2, 0, 3};

  Poly2<N> ppol = Poly2<N>(10, &coef[0], m);
  // ppol.print();
  return ppol;
}
template <class N>
Poly2<N> NApproxPoly2::calculate() {
  Poly2<N> ff = f->get<N>();

  if (i == 1) ff = ff * -1;

  PV2<N> pp = p->get<N>();
  N a = pp.x.lbP();
  N b = pp.y.lbP();

  vector<N> coef;
  vector<int> m;

  vector<Poly2<N>> functions;
  functions.push_back(ff);

  vector<N> constants;
  constants.push_back(Parameter::constant(1));

  N denom = Parameter::constant(1);

  coef.push_back(ff.value(a, b).lbP());
  m.push_back(0);
  m.push_back(0);

  for (int i = 1; i <= degree; i++) {
    denom = denom * Parameter::constant(i);

    vector<Poly2<N>> current_functions;
    vector<N> current_constants;

    current_functions.push_back(functions[0].derX());
    for (int j = 0; j < functions.size(); j++) {
      current_functions.push_back(functions[j].derY());
    }

    current_constants.push_back(Parameter::constant(1));

    for (int j = 0; j < constants.size() - 1; j++) {
      current_constants.push_back(constants[j] + constants[j + 1]);
    }

    current_constants.push_back(Parameter::constant(1));

    int x = i;
    int y = 0;

    for (int j = 0; j < current_functions.size(); j++) {
      N value;
      if (i == degree)
        value = (current_constants[j] / denom) *
                current_functions[j].value(pp.x, pp.y);
      else
        value =
            (current_constants[j] / denom) * current_functions[j].value(a, b);

      coef.push_back(value.lbP());
      m.push_back(x--);
      m.push_back(y++);
    }

    functions = current_functions;
    constants = current_constants;
  }

  return Poly2<N>(coef.size(), &coef[0], &m[0]);
}

template <class N>
N NewtonPoly::newtonRoot(int i) {
  // return Parameter::constant(i);
  int k = 0;
  while ((1 << (k + 1)) - 2 < i) {
    k++;
  }

  return N::constant((2.0 * i - 3.0 * ((1 << k) - 1)) / (1 << k));

  // return randomNumber(-1, 1);

  // return ((Parameter::constant(1).acos() * (2*Parameter::constant(i) - 1)) /
  // (2*Parameter::constant(n))).cos();
}

template Parameter NewtonPoly::newtonRoot<Parameter>(int i);
template PParameter NewtonPoly::newtonRoot<PParameter>(int i);
template MParameter NewtonPoly::newtonRoot<MParameter>(int i);

template <class N>
Poly<N> NewtonPoly::calculate() {
  // degree 0 basis poly, (1)
  if (!f) {
    return Poly<N>(N(1));
  }
  // degree > 0 basis poly
  else {
    vector<N> xy(2);
    xy[0] = -newtonRoot<N>(r);
    xy[1] = Parameter::constant(1);

    Poly<N> xr(xy);

    return f->get<N>() * xr;
  }
}

// i corresponds to the index of the basis poly, degree of that poly is ==
// itemplate<class N>
template <class N>
Poly<N> NewtonPoly::basis(int i) {
  // generate all of the basis polynomials up to the desired
  while (bas.size() <= i) {
    // Default Newton poly, 1
    if (bas.size() == 0) {
      bas.push_back(new NewtonPoly());
    }
    // Generic Newton poly, p18 = p17 * (x-18)
    else {
      bas.push_back(new NewtonPoly(bas.size() - 1, bas.back()));
    }
  }

  Poly<N> p = bas[i]->get<N>();

  return p;
}

template <class N>
Scalar<N> EquidistantParameter::calculate() {
  Poly2<N> ff = f->get<N>();
  Poly2<N> gg = g->get<N>();
  PV2<N> qq = q->get<N>();
  PV2<N> uu = u->get<N>();

  Poly2<N> ffx = ff.derX();
  Poly2<N> ffy = ff.derY();
  Poly2<N> ggx = gg.derX();
  Poly2<N> ggy = gg.derY();

  N fq = ff.value(&qq.x);
  N gq = gg.value(&qq.x);

  PV2<N> gradF = PV2<N>(ffx.value(&qq.x), ffy.value(&qq.x));
  PV2<N> gradG = PV2<N>(ggx.value(&qq.x), ggy.value(&qq.x));

  // use tangents instead of grad?
  // PV2 gradF = PV2(-ffy.value(&qq.x), ffx.value(&qq.x));
  // PV2 gradG = PV2(-ggy.value(&qq.x), ggx.value(&qq.x));

  N gradGLength = gradG.length();
  N gradFLength = gradF.length();

  N A = gradGLength * fq;
  N B = gradGLength * (gradF.dot(uu));

  N C = gradFLength * gq;
  N D = gradFLength * (gradG.dot(uu));

  //(f1(q) + grad f1(q) * u s) / |grad f1(q)| = (f2(q) + grad f2(q) * u s) /
  //|grad f2(q)|
  // A + B*S = C + D*S

  // s = (A - C) / (D - B)

  return (A - C) / (D - B);
}

template <class N>
PV2<N> LinearizationIntersection::calculate() {
  Poly2<N> ff = f->get<N>();
  Poly2<N> gg = g->get<N>();
  PV2<N> pp = p->get<N>();
  PV2<N> qq = (q ? q->get<N>() : pp);

  // To get v, we solve f1(p+v)=f2(p+v)=0 using linearization about p:
  // f1(p) + grad f1(p) * v = 0
  // f2(p) + grad f2(p) * v = 0

  PV2<N> fr(ff.derX().value(&pp.x), ff.derY().value(&pp.x));
  PV2<N> gr(gg.derX().value(&qq.x), gg.derY().value(&qq.x));
  PV2<N> b(-ff.value(&pp.x), -gg.value(&qq.x));

  N d = fr.x * gr.y - fr.y * gr.x;

  // pivot first row
  if ((fr.x.abs() - gr.x.abs()).sign(false) >= 0) {
    N c = gr.x / fr.x;

    gr.x = gr.x - c * fr.x;
    gr.y = gr.y - c * fr.y;
    b.y = b.y - c * b.x;

    N y = b.y / gr.y;

    N x = (b.x - fr.y * y) / fr.x;

    return PV2<N>(x, y);

  }
  // pivot on the second row
  else {
    N c = fr.x / gr.x;

    fr.x = fr.x - c * gr.x;
    fr.y = fr.y - c * gr.y;
    b.x = b.x - c * b.y;

    N y = b.x / fr.y;

    N x = (b.y - gr.y * y) / gr.x;

    return PV2<N>(x, y);
  }
}
template <class N>
Poly2<N> LinearizationPoly::calculate() {
  Poly2<N> ff = f->get<N>();
  PV2<N> pp = p->get<N>();

  Poly2<N> fx = ff.derX();
  Poly2<N> fy = ff.derY();

  N xy[3];
  xy[0] = ff.value(&pp.x);
  xy[1] = fx.value(&pp.x);
  xy[2] = fy.value(&pp.x);

  int m[6] = {0, 0, 1, 0, 0, 1};

  return Poly2<N>(3, xy, m);
}

PTR<Line> LinearizationPoly::getLine() {
  PTR<Object<Scalar>> fp = new PolyScalar(f, p);
  PTR<Object<PV2>> grad = new GradientVector(p, f);
  PTR<Object<Scalar>> gl = new DotScalar(grad, grad);

  PTR<Object<PV2>> q = new SumPoint(
      p, new NegativeVector(new ScaleVector(new DivAB(fp, gl), grad)));
  PTR<Object<PV2>> w = new Rot90(grad, nullptr);

  return new Line(q, w);
}

//--Predicates-------------------------------------------------------
template <class N>
N Intersects::calculate() {
  PV2<N> aa = a->get<N>();
  PV2<N> bb = b->get<N>();

  // REFACTOR
  if (aa.x.intersects(bb.x) && aa.y.intersects(bb.y)) {
    return N(1);
  } else {
    return N(-1);
  }
}

template <class N>
N LeftOf::calculate() {
  PV2<N> ap = p->get<N>();
  PV2<N> bp = ap + v->get<N>() * Parameter::constant(1);
  PV2<N> cp = o->get<N>();

  return (cp - ap).cross(bp - ap);
}

template Parameter LeftOf::calculate<Parameter>();
template PParameter LeftOf::calculate<PParameter>();
template MParameter LeftOf::calculate<MParameter>();

template <class N>
N LessThan::calculate() {
  return (a->get<N>() - b->get<N>());
}

template Parameter LessThan::calculate<Parameter>();
template PParameter LessThan::calculate<PParameter>();
template MParameter LessThan::calculate<MParameter>();

template <class N>
N LessThanAbs::calculate() {
  return (a->get<N>().x.abs() - b->get<N>().x.abs());
}

template <class N>
N Side::calculate() {
  return a->get<N>().cross(b->get<N>());
}

template Parameter Side::calculate<Parameter>();
template PParameter Side::calculate<PParameter>();
template MParameter Side::calculate<MParameter>();

template <class N>
N PolySign::calculate() {
  // static int poly_count = 0;
  // poly_count++;
  // printf("poly_count: %d\n", poly_count);
  N p_xy[2] = {p->get<N>().x, p->get<N>().y};
  return f->get<N>().value(p_xy);
}

template Parameter PolySign::calculate<Parameter>();
template PParameter PolySign::calculate<PParameter>();
template MParameter PolySign::calculate<MParameter>();

template <class N>
N Descartes2::calculate() {
  PV2<N> lp = l->get<N>();
  PV2<N> up = u->get<N>();
  Poly2<N> fp = f->get<N>();

  N xl = lp.x;
  N xu = up.x;
  N yl = lp.y;
  N yu = up.y;

  N a[fp.degx + 1][fp.degy + 1];
  for (int i = 0; i <= fp.degx; i++) {
    for (int j = 0; j <= fp.degy; j++) {
      a[i][j] = Parameter::constant(0);
    }
  }

  for (int i = 0; i < fp.size(); i++) {
    a[fp.m[2 * i]][fp.m[2 * i + 1]] = fp.a[i];
  }

  // To shift by s=a.  Let a[0..d] be the array of coefficients.  a[i] is the
  // coefficient of x^i.
  //
  // for i = 0 to d-1 (inclusive)
  //   for i' = d-1 to i
  //     a[i'] += s * a[i'+1]
  //
  // Rather than inverting, just remember the inversion.  To do the next shift
  // by s=1/(b-a)
  //
  // for i = 0 to d-1 (inclusive)
  //   for i' = d-1 to i
  //     a[d-i'] += s * a[d-(i'+1)]
  //
  // Now look for sign changes in a[0..d].
  //
  //
  // For a bivariate, let a[i][j] be the coefficient of x^i y^j.
  //
  // a[i'] += s * a[i'+1]
  //
  // means
  //
  // for j = 0 to d
  //   a[i'][j] += s * a[i'+1][j]
  //

  N s;

  // shift by xl
  s = xl;
  for (int i = 0; i < fp.degx; i++)
    for (int j = fp.degx - 1; j >= i; j--)
      for (int k = 0; k <= fp.degy; k++) a[j][k] = a[j][k] + s * a[j + 1][k];

  // shift by 1/(xu-xl)
  s = 1.0 / (xu - xl);
  for (int i = 0; i < fp.degx; i++)
    for (int j = fp.degx - 1; j >= i; j--)
      for (int k = 0; k <= fp.degy; k++)
        a[fp.degx - j][k] = a[fp.degx - j][k] + s * a[fp.degx - (j + 1)][k];

  // Suppose y goes from c to d.
  //
  // Check that each y polynomial g_i(y) = sum{j=0 to d) a[i][j] y^j has the
  // same sign at x=c using Horner.
  //
  // If yes, check that there are no zeros in [c,d] using a shift by c,
  // inverting, and by 1/(d-c).

  int sgn = 0;

  for (int i = 0; i <= fp.degx; i++) {
    N v = Parameter::constant(0);
    for (int j = fp.degy; j >= 0; j--) {
      v = v * yl + a[i][j];
    }

    int vsgn = v.sign();

    if (sgn == 0) {
      sgn = vsgn;
    } else if (sgn != vsgn) {
      return N(1);
    }
  }

  for (int i = 0; i <= fp.degx; i++) {
    // shift by yl
    s = yl;
    for (int j = 0; j < fp.degy; j++)
      for (int k = fp.degy - 1; k >= j; k--)
        a[i][k] = a[i][k] + s * a[i][k + 1];

    // shift by 1/(yu-yl)
    s = 1.0 / (yu - yl);
    for (int j = 0; j < fp.degy; j++)
      for (int k = fp.degy - 1; k >= j; k--)
        a[i][fp.degy - k] = a[i][fp.degy - k] + s * a[i][fp.degy - (k + 1)];

    // descartes
    sgn = a[i][0].sign();

    for (int j = 1; j <= fp.degy; j++) {
      if (sgn != a[i][j].sign()) return N(1);
    }
  }

  return N(0);
}

template Parameter Descartes2::calculate<acp::Parameter>();
template MParameter Descartes2::calculate<acp::MParameter>();
template PParameter Descartes2::calculate<acp::PParameter>();
template Parameter Intersects::calculate<acp::Parameter>();
template MParameter Intersects::calculate<acp::MParameter>();
template PParameter Intersects::calculate<acp::PParameter>();
template Parameter LessThanAbs::calculate<acp::Parameter>();
template MParameter LessThanAbs::calculate<acp::MParameter>();
template PParameter LessThanAbs::calculate<acp::PParameter>();
template Poly2<Parameter> NApproxPoly2::calculate<acp::Parameter>();
template Poly2<MParameter> NApproxPoly2::calculate<acp::MParameter>();
template Poly2<PParameter> NApproxPoly2::calculate<acp::PParameter>();
template Poly<Parameter> ShiftOutPoly::calculate<acp::Parameter>();
template Poly<MParameter> ShiftOutPoly::calculate<acp::MParameter>();
template Poly<PParameter> ShiftOutPoly::calculate<acp::PParameter>();
template Poly2<Parameter> QuadApproxPoly2::calculate<acp::Parameter>();
template Poly2<MParameter> QuadApproxPoly2::calculate<acp::MParameter>();
template Poly2<PParameter> QuadApproxPoly2::calculate<acp::PParameter>();
template Poly2<Parameter> CubicApproxPoly2::calculate<acp::Parameter>();
template Poly2<MParameter> CubicApproxPoly2::calculate<acp::MParameter>();
template Poly2<PParameter> CubicApproxPoly2::calculate<acp::PParameter>();
template Poly2<Parameter> LinearizationPoly::calculate<acp::Parameter>();
template Poly2<MParameter> LinearizationPoly::calculate<acp::MParameter>();
template Poly2<PParameter> LinearizationPoly::calculate<acp::PParameter>();
template PV2<Parameter> QuadCriticalPoint::calculate<acp::Parameter>();
template PV2<MParameter> QuadCriticalPoint::calculate<acp::MParameter>();
template PV2<PParameter> QuadCriticalPoint::calculate<acp::PParameter>();
template PV2<Parameter> EllipseCriticalPoint::calculate<acp::Parameter>();
template PV2<MParameter> EllipseCriticalPoint::calculate<acp::MParameter>();
template PV2<PParameter> EllipseCriticalPoint::calculate<acp::PParameter>();
template Scalar<Parameter> EquidistantParameter::calculate<acp::Parameter>();
template Scalar<MParameter> EquidistantParameter::calculate<acp::MParameter>();
template Scalar<PParameter> EquidistantParameter::calculate<acp::PParameter>();
template PV2<Parameter> LinearizationIntersection::calculate<acp::Parameter>();
template PV2<MParameter>
LinearizationIntersection::calculate<acp::MParameter>();
template PV2<PParameter>
LinearizationIntersection::calculate<acp::PParameter>();
template Poly<Parameter> LinePoly::calculate<acp::Parameter>();
template Poly<MParameter> LinePoly::calculate<acp::MParameter>();
template Poly<PParameter> LinePoly::calculate<acp::PParameter>();

template Poly<Parameter> interpolate(PTR<Line> l, PTR<Object<Poly2>> f);
template Poly<PParameter> interpolate(PTR<Line> l, PTR<Object<Poly2>> f);
template Poly<MParameter> interpolate(PTR<Line> l, PTR<Object<Poly2>> f);
