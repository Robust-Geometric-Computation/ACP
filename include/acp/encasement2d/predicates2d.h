#ifndef PREDICATES2D_H
#define PREDICATES2D_H

#include "acp/poly/poly2.h"

using namespace acp;
using namespace std;

// Macro to create one-argument primitive.
#define Primitive1(P, t1, v1)       \
  class P : public acp::Primitive { \
    t1 v1;                          \
    DeclareSign;                    \
                                    \
   public:                          \
    P(t1 v1) : v1(v1) {}            \
  };

// Macro to create two-argument primitive.
#define Primitive2(P, t1, v1, t2, v2)   \
  class P : public acp::Primitive {     \
    t1 v1;                              \
    t2 v2;                              \
    DeclareSign;                        \
                                        \
   public:                              \
    P(t1 v1, t2 v2) : v1(v1), v2(v2) {} \
  };

// Macro to create three-argument primitive.
#define Primitive3(P, t1, v1, t2, v2, t3, v3)          \
  class P : public acp::Primitive {                    \
    t1 v1;                                             \
    t2 v2;                                             \
    t3 v3;                                             \
    DeclareSign;                                       \
                                                       \
   public:                                             \
    P(t1 v1, t2 v2, t3 v3) : v1(v1), v2(v2), v3(v3) {} \
  };

// Macro to create four-argument primitive.
#define Primitive4(P, t1, v1, t2, v2, t3, v3, t4, v4)                 \
  class P : public acp::Primitive {                                   \
    t1 v1;                                                            \
    t2 v2;                                                            \
    t3 v3;                                                            \
    t4 v4;                                                            \
    DeclareSign;                                                      \
                                                                      \
   public:                                                            \
    P(t1 v1, t2 v2, t3 v3, t4 v4) : v1(v1), v2(v2), v3(v3), v4(v4) {} \
  };

// Macro to create five-argument primitive.
#define Primitive5(P, t1, v1, t2, v2, t3, v3, t4, v4, t5, v5) \
  class P : public acp::Primitive {                           \
    t1 v1;                                                    \
    t2 v2;                                                    \
    t3 v3;                                                    \
    t4 v4;                                                    \
    t5 v5;                                                    \
    DeclareSign;                                              \
                                                              \
   public:                                                    \
    P(t1 v1, t2 v2, t3 v3, t4 v4, t5 v5)                      \
        : v1(v1), v2(v2), v3(v3), v4(v4), v5(v5) {}           \
  };

class Line;

class InputPointConst : public Object<PV2> {
 public:
  InputPointConst(double x, double y)
      : Object<PV2>(PV2<Parameter>::constant(x, y)) {}
  InputPointConst(const PV2<Parameter>& x) : Object<PV2>(x) {}
};

class InputScalar : public Object<Scalar> {
 public:
  InputScalar(double r) : Object<Scalar>(Parameter::input(r)) {}
  InputScalar(const Parameter& x) : Object<Scalar>(x) {}
};

class PVScalar : public Object<Scalar> {
  PTR<Object<PV2>> point;
  int index;
  DeclareCalculate(Scalar) {
    return index < 1 ? point->get<N>().x : point->get<N>().y;
  }

 public:
  PVScalar(PTR<Object<PV2>> point, int index) : point(point), index(index) {}
};

class AverageAB : public Object<Scalar> {
  double t;
  PTR<Object<Scalar>> a;
  PTR<Object<Scalar>> b;
  DeclareCalculate(Scalar) {
    N aa = a->get<N>();
    return aa + t * (b->get<N>() - aa);
  }

 public:
  AverageAB(double t, PTR<Object<Scalar>> a, PTR<Object<Scalar>> b)
      : t(t), a(a), b(b) {}
};

class SumAB : public Object<Scalar> {
  PTR<Object<Scalar>> a;
  PTR<Object<Scalar>> b;
  DeclareCalculate(Scalar) { return a->get<N>() + b->get<N>(); }

 public:
  SumAB(PTR<Object<Scalar>> a, PTR<Object<Scalar>> b) : a(a), b(b) {}
};

class DiffAB : public Object<Scalar> {
  PTR<Object<Scalar>> a;
  PTR<Object<Scalar>> b;
  DeclareCalculate(Scalar) { return a->get<N>() - b->get<N>(); }

 public:
  DiffAB(PTR<Object<Scalar>> a, PTR<Object<Scalar>> b) : a(a), b(b) {}
};

class DivAB : public Object<Scalar> {
  PTR<Object<Scalar>> a;
  PTR<Object<Scalar>> b;
  DeclareCalculate(Scalar) { return a->get<N>() / b->get<N>(); }

 public:
  DivAB(PTR<Object<Scalar>> a, PTR<Object<Scalar>> b) : a(a), b(b) {}
};

class ProdAB : public Object<Scalar> {
  PTR<Object<Scalar>> a;
  PTR<Object<Scalar>> b;
  DeclareCalculate(Scalar) { return a->get<N>() * b->get<N>(); }

 public:
  ProdAB(PTR<Object<Scalar>> a, PTR<Object<Scalar>> b) : a(a), b(b) {}
};

class NegativeScalar : public Object<Scalar> {
  PTR<Object<Scalar>> a;
  DeclareCalculate(Scalar) { return -1 * (a->get<N>()); }

 public:
  NegativeScalar(PTR<Object<Scalar>> a) : a(a) {}
};

class Poly1Scalar : public Object<Scalar> {
  PTR<Object<Poly>> f;
  PTR<Object<Scalar>> p;
  DeclareCalculate(Scalar) { return f->get<N>().value(p->get<N>()); }

 public:
  Poly1Scalar(PTR<Object<Poly>> f, PTR<Object<Scalar>> p) : f(f), p(p) {}
};

class PolyScalar : public Object<Scalar> {
  PTR<Object<Poly2>> f;
  PTR<Object<PV2>> p;
  DeclareCalculate(Scalar) {
    PV2<N> pp = p->get<N>();
    N xy[2] = {pp.x, pp.y};
    return f->get<N>().value(xy);
  }

 public:
  PolyScalar(PTR<Object<Poly2>> f, PTR<Object<PV2>> p) : f(f), p(p) {}
};

class DotScalar : public Object<Scalar> {
  PTR<Object<PV2>> a;
  PTR<Object<PV2>> b;
  DeclareCalculate(Scalar) { return a->get<N>().dot(b->get<N>()); }

 public:
  DotScalar(PTR<Object<PV2>> a, PTR<Object<PV2>> b) : a(a), b(b) {}
};

class CrossScalar : public Object<Scalar> {
  PTR<Object<PV2>> a;
  PTR<Object<PV2>> b;
  DeclareCalculate(Scalar) { return a->get<N>().cross(b->get<N>()); }

 public:
  CrossScalar(PTR<Object<PV2>> a, PTR<Object<PV2>> b) : a(a), b(b) {}
};

class VectorLength : public Object<Scalar> {
  PTR<Object<PV2>> v;
  DeclareCalculate(Scalar) { return v->get<N>().length(); }

 public:
  VectorLength(PTR<Object<PV2>> v) : v(v) {}
};

class InputPoint : public Object<PV2> {
 public:
  InputPoint(const double x, const double y)
      : Object<PV2>(PV2<Parameter>(Parameter::input(x), Parameter::input(y))) {}
  InputPoint(const PV2<Parameter>& p) : Object<PV2>(p) {}
};

class DoubleInputInterval : public Object<Scalar> {
  double l, u;
  DeclareCalculate(Scalar) { return N(Parameter::interval(l, u)); }

 public:
  DoubleInputInterval(double l, double u) : l(l), u(u) {}
};

class Interval : public Object<Scalar> {
  PTR<Object<Scalar>> p;
  DeclareCalculate(Scalar) {
#ifdef BLEEN
    Parameter pp = p->get<N>();
    Parameter plb = pp.lbP();
    Parameter pub = pp.ubP();
    return Parameter(plb, pub);
#endif
    return p->get<N>();
  }

 public:
  Interval(const PTR<Object<Scalar>>& p) : p(p) {}
};

class ParameterInterval : public Object<Scalar> {
  PTR<Object<Scalar>> a;
  PTR<Object<Scalar>> b;
  DeclareCalculate(Scalar) {
    return N(Parameter(a->get<N>().x.interval(b->get<N>())));
  }

 public:
  ParameterInterval(PTR<Object<Scalar>> a, PTR<Object<Scalar>> b)
      : a(a), b(b) {}
};

class ParameterPoint : public Object<PV2> {
  PTR<Object<Scalar>> a;
  PTR<Object<Scalar>> b;
  DeclareCalculate(PV2) { return PV2<N>(a->get<N>(), b->get<N>()); }

 public:
  ParameterPoint(PTR<Object<Scalar>> a, PTR<Object<Scalar>> b) : a(a), b(b) {}
};

class SumPoint : public Object<PV2> {
  PTR<Object<PV2>> a;
  PTR<Object<PV2>> b;
  DeclareCalculate(PV2) { return a->get<N>() + b->get<N>(); }

 public:
  SumPoint(PTR<Object<PV2>> a, PTR<Object<PV2>> b) : a(a), b(b) {}
};

class DifferencePoint : public Object<PV2> {
  PTR<Object<PV2>> a;
  PTR<Object<PV2>> b;
  DeclareCalculate(PV2) { return a->get<N>() - b->get<N>(); }

 public:
  DifferencePoint(PTR<Object<PV2>> a, PTR<Object<PV2>> b) : a(a), b(b) {}
};

class CornerPoint : public Object<PV2> {
  PTR<Object<PV2>> p;
  int i, j;
  DeclareCalculate(PV2) {
    PV2<N> pp = p->get<N>();
    return PV2<N>((i == 0 ? pp.x.lbP() : pp.x.ubP()),
                  (j == 0 ? pp.y.lbP() : pp.y.ubP()));
  }

 public:
  CornerPoint(PTR<Object<PV2>> p, int i, int j) : p(p), i(i), j(j) {}
};

class VectorAB : public Object<PV2> {
  PTR<Object<PV2>> a;
  PTR<Object<PV2>> b;
  DeclareCalculate(PV2) { return b->get<N>() - a->get<N>(); }

 public:
  VectorAB(PTR<Object<PV2>> a, PTR<Object<PV2>> b) : a(a), b(b) {}
};

class CrossVector : public Object<Scalar> {
  PTR<Object<PV2>> a;
  PTR<Object<PV2>> b;
  DeclareCalculate(Scalar) { return b->get<N>().cross(a->get<N>()); }

 public:
  CrossVector(PTR<Object<PV2>> a, PTR<Object<PV2>> b) : a(a), b(b) {}
};

class MidVector : public Object<PV2> {
  PTR<Object<PV2>> a;
  PTR<Object<PV2>> b;
  DeclareCalculate(PV2) { return 0.5 * (a->get<N>() + b->get<N>()); }

 public:
  MidVector(PTR<Object<PV2>> a, PTR<Object<PV2>> b) : a(a), b(b) {}
};

class ScaleVector : public Object<PV2> {
  PTR<Object<PV2>> v;
  PTR<Object<Scalar>> s;
  double ds;
  bool param;
  DeclareCalculate(PV2) {
    if (param)
      return s->get<N>().x * v->get<N>();
    else
      return ds * v->get<N>();
  }

 public:
  ScaleVector(PTR<Object<Scalar>> s, PTR<Object<PV2>> v)
      : s(s), v(v), param(true) {}
  ScaleVector(double ds, PTR<Object<PV2>> v) : ds(ds), v(v), param(false) {}
};

class ScaleParameter : public Object<Scalar> {
  PTR<Object<Scalar>> v;
  PTR<Object<Scalar>> s;
  double ds;
  bool param;
  DeclareCalculate(Scalar) {
    if (param)
      return s->get<N>() * v->get<N>();
    else
      return ds * v->get<N>();
  }

 public:
  ScaleParameter(PTR<Object<Scalar>> s, PTR<Object<Scalar>> v)
      : s(s), v(v), param(true) {}
  ScaleParameter(double ds, PTR<Object<Scalar>> v)
      : ds(ds), v(v), param(false) {}
};

class TangentVector : public Object<PV2> {
  PTR<Object<PV2>> p;
  PTR<Object<Poly2>> f;
  DeclareCalculate(PV2) {
    Poly2<N> ff = f->get<N>();
    Poly2<N> dx = ff.derX();
    Poly2<N> dy = ff.derY();
    PV2<N> pp = p->get<N>();
    N xy[2] = {pp.x, pp.y};
    // gradient (dx, dy) rotated by -90: (dy, -dx)
    return PV2<N>(dy.value(xy), -dx.value(xy));
  }

 public:
  TangentVector(PTR<Object<PV2>> p, PTR<Object<Poly2>> f) : p(p), f(f) {}
};

class GradientVector : public Object<PV2> {
  PTR<Object<PV2>> p;
  PTR<Object<Poly2>> f;
  DeclareCalculate(PV2) {
    Poly2<N> ff = f->get<N>();
    Poly2<N> dx = ff.derX();
    Poly2<N> dy = ff.derY();
    PV2<N> pp = p->get<N>();
    N xy[2] = {pp.x, pp.y};
    return PV2<N>(dx.value(xy), dy.value(xy));
  }

 public:
  GradientVector(PTR<Object<PV2>> p, PTR<Object<Poly2>> f) : p(p), f(f) {}
};

class NormalVector : public Object<PV2> {
  PTR<Object<PV2>> v;
  DeclareCalculate(PV2) { return (v->get<N>()).unit(); }

 public:
  NormalVector(PTR<Object<PV2>> v) : v(v) {}
};

class AverageVector : public Object<PV2> {
  double t;
  PTR<Object<PV2>> v;
  PTR<Object<PV2>> w;
  DeclareCalculate(PV2) {
    PV2<N> vv = v->get<N>();
    return vv + (t * (w->get<N>() - vv));
  }

 public:
  AverageVector(double t, PTR<Object<PV2>> v, PTR<Object<PV2>> w)
      : t(t), v(v), w(w) {}
};

class AveragePoint : public Object<PV2> {
  PTR<Object<Scalar>> t;
  PTR<Object<PV2>> a;
  PTR<Object<PV2>> b;
  DeclareCalculate(PV2) {
    PV2<N> aa = a->get<N>();
    return aa + (t->get<N>().x * (b->get<N>() - aa));
  }

 public:
  AveragePoint(PTR<Object<Scalar>> t, PTR<Object<PV2>> a, PTR<Object<PV2>> b)
      : t(t), a(a), b(b) {}
};

// typedef PV2 (*transform_function)(PV2 p);

class TransformPoint : public Object<PV2> {
  PTR<Object<PV2>> p;
  double post_scale, draw_scale, x_off, y_off;
  DeclareCalculate(PV2) {
    PV2<N> pp = p->get<N>();
    return PV2<N>(post_scale * ((pp.x / draw_scale) - x_off),
                  post_scale * ((pp.y / draw_scale) - y_off));
  }

 public:
  TransformPoint(PTR<Object<PV2>> p, double post_scale, double draw_scale,
                 double x_off, double y_off)
      : p(p),
        post_scale(post_scale),
        draw_scale(draw_scale),
        x_off(x_off),
        y_off(y_off) {}
};

class HalfVector : public AverageVector {
 public:
  HalfVector(PTR<Object<PV2>> v, PTR<Object<PV2>> w)
      : AverageVector(0.5, v, w) {}
};

class Rot90 : public Object<PV2> {
  PTR<Object<PV2>> a;
  PTR<Object<PV2>> b;
  DeclareCalculate(PV2) {
    PV2<N> v = (b == nullptr ? a->get<N>() : b->get<N>() - a->get<N>());
    return PV2<N>(-v.y, v.x);
  }

 public:
  Rot90(PTR<Object<PV2>> a, PTR<Object<PV2>> b = nullptr) : a(a), b(b) {}
};

class Rot45 : public Object<PV2> {
  PTR<Object<PV2>> a;
  PTR<Object<PV2>> b;
  DeclareCalculate(PV2) {
    PV2<N> v = (b == nullptr ? a->get<N>() : b->get<N>() - a->get<N>());
    double c45 = cos(M_PI / 4.0);
    double s45 = sin(M_PI / 4.0);
    return PV2<N>(c45 * v.x - s45 * v.y, s45 * v.x + c45 * v.y);
  }

 public:
  Rot45(PTR<Object<PV2>> a, PTR<Object<PV2>> b) : a(a), b(b) {}
};

class NegativeVector : public Object<PV2> {
  PTR<Object<PV2>> v;
  DeclareCalculate(PV2) {
    PV2<N> w = v->get<N>();
    return PV2<N>(-w.x, -w.y);
  }

 public:
  NegativeVector(PTR<Object<PV2>> v) : v(v) {}
};

class InputPoly2 : public Object<Poly2> {
 public:
  InputPoly2(const Poly2<Parameter>& poly) : Object<Poly2>(poly) {}
};

class AddPoly2 : public Object<Poly2> {
  PTR<Object<Poly2>> f;
  PTR<Object<Poly2>> g;

  DeclareCalculate(Poly2) { return f->get<N>() + g->get<N>(); }

 public:
  AddPoly2(const PTR<Object<Poly2>> f, const PTR<Object<Poly2>> g)
      : f(f), g(g) {}
};

class AddNormalPoly2 : public Object<Poly2> {
  PTR<Object<Poly2>> f;
  PTR<Object<Poly2>> g;
  PTR<Object<PV2>> p;
  PTR<Object<PV2>> q;

  DeclareCalculate(Poly2) {
    PV2<N> pp = p->get<N>();
    PV2<N> qq = q->get<N>();
    Poly2<N> ff = f->get<N>();
    Poly2<N> gg = g->get<N>();
    Poly2<N> fx = ff.derX();
    Poly2<N> fy = ff.derY();
    Poly2<N> gx = gg.derX();
    Poly2<N> gy = gg.derY();
    PV2<N> grf(fx.value(&pp.x), fy.value(&pp.x));
    PV2<N> grg(gx.value(&qq.x), gy.value(&qq.x));

    return (f->get<N>() * (1.0 / grf.length())) +
           (g->get<N>() * (1.0 / grg.length()));
  }

 public:
  AddNormalPoly2(const PTR<Object<Poly2>> f, const PTR<Object<Poly2>> g,
                 PTR<Object<PV2>> p, PTR<Object<PV2>> q)
      : f(f), g(g), p(p), q(q) {}
};

class NegativePoly2 : public Object<Poly2> {
  PTR<Object<Poly2>> f;
  DeclareCalculate(Poly2) { return f->get<N>() * -1; }

 public:
  NegativePoly2(const PTR<Object<Poly2>> f) : f(f) {}
};

class OffsetPoly2 : public Object<Poly2> {
  PTR<Object<Poly2>> f;
  PTR<Object<Scalar>> a;
  DeclareCalculate(Poly2) {
    Poly2<N> pa;
    pa.add(0, 0, a->get<N>());
    // return f->get<N>() + (Poly2::one()*a->get<N>());
    return f->get<N>() + pa;
  }

 public:
  OffsetPoly2(const PTR<Object<Poly2>> f, const PTR<Object<Scalar>> a)
      : f(f), a(a) {}
  OffsetPoly2(const PTR<Object<Poly2>> f, const double a)
      : f(f), a(new InputScalar(Parameter::constant(a))) {}
};

class Ellipse : public Object<Poly2> {
  static int count;

 public:
  int id;
  Ellipse(const Poly2<Parameter>& poly) : id(count++), Object<Poly2>(poly) {}
};

class DirectionalDerivativePoly2 : public Object<Poly2> {
  PTR<Object<Poly2>> f;
  PTR<Object<PV2>> v;
  DeclareCalculate(Poly2) {
    Poly2<N> ff = f->get<N>();
    PV2<N> vv = v->get<N>();
    return (ff.derX() * vv.x) + (ff.derY() * vv.y);
  }

 public:
  DirectionalDerivativePoly2(const PTR<Object<Poly2>> f,
                             const PTR<Object<PV2>> v)
      : f(f) {
    this->v = new NormalVector(v);
  }
};

class DerXPoly2 : public Object<Poly2> {
  PTR<Object<Poly2>> f;
  DeclareCalculate(Poly2) { return f->get<N>().derX(); }

 public:
  DerXPoly2(const PTR<Object<Poly2>> f) : f(f) {}
};

class DerYPoly2 : public Object<Poly2> {
  PTR<Object<Poly2>> f;
  DeclareCalculate(Poly2) { return f->get<N>().derY(); }

 public:
  DerYPoly2(const PTR<Object<Poly2>> f) : f(f) {}
};

class HessianPoly2 : public Object<Poly2> {
  PTR<Object<Poly2>> f;
  DeclareCalculate(Poly2) {
    Poly2<N> ff = f->get<N>();
    Poly2<N> fx = ff.derX();
    Poly2<N> fy = ff.derY();
    Poly2<N> fxx = fx.derX();
    Poly2<N> fyy = fy.derY();
    Poly2<N> fxy = fx.derY();

    return (fxx * fyy) - (fxy * fxy);
  }

 public:
  HessianPoly2(const PTR<Object<Poly2>> f) : f(f) {}
};

class QuadApproxPoly2 : public Object<Poly2> {
  PTR<Object<Poly2>> f;
  PTR<Object<PV2>> p;
  int i;
  DeclareCalculate(Poly2);

 public:
  QuadApproxPoly2(PTR<Object<Poly2>> f, PTR<Object<PV2>> p, int i)
      : f(f), p(p), i(i) {}
};

class QuadCriticalPoint : public Object<PV2> {
  PTR<Object<Poly2>> f;
  DeclareCalculate(PV2);

 public:
  QuadCriticalPoint(PTR<Object<Poly2>> f) : f(f) {}
};

class CubicApproxPoly2 : public Object<Poly2> {
  PTR<Object<Poly2>> f;
  PTR<Object<PV2>> p;
  int i;
  DeclareCalculate(Poly2);

 public:
  CubicApproxPoly2(PTR<Object<Poly2>> f, PTR<Object<PV2>> p, int i)
      : f(f), p(p), i(i) {}
};

class NApproxPoly2 : public Object<Poly2> {
  PTR<Object<Poly2>> f;
  PTR<Object<PV2>> p;
  int degree;
  int i;
  DeclareCalculate(Poly2);

 public:
  NApproxPoly2(PTR<Object<Poly2>> f, PTR<Object<PV2>> p, int i, int degree)
      : f(f), p(p), i(i), degree(degree) {}
};

class NewtonPoly : public Object<Poly> {
  PTR<Object<Poly>> f;
  int r;

  static vector<PTR<Object<Poly>>> bas;

  DeclareCalculate(Poly);

 public:
  NewtonPoly(int r = -1, PTR<Object<Poly>> f = nullptr) : r(r), f(f) {}

  template <class N>
  static Poly<N> basis(int i);

  template <class N>
  static N newtonRoot(int i);
};

class LinearCombinationPoly : public PTR<Object<Poly>> {
  vector<PTR<Object<Scalar>>> a;
  vector<PTR<Object<Poly>>> p;

  DeclareCalculate(Poly);

 public:
  LinearCombinationPoly(vector<PTR<Object<Scalar>>>& a,
                        vector<PTR<Object<Poly>>>& p)
      : a(a), p(p) {}
};

class EquidistantParameter : public Object<Scalar> {
  PTR<Object<Poly2>> f;
  PTR<Object<Poly2>> g;
  PTR<Object<PV2>> q;
  PTR<Object<PV2>> u;
  DeclareCalculate(Scalar);

 public:
  EquidistantParameter(PTR<Object<Poly2>> f, PTR<Object<Poly2>> g,
                       PTR<Object<PV2>> q, PTR<Object<PV2>> u)
      : f(f), g(g), q(q), u(u) {}
};

class LinearizationIntersection : public Object<PV2> {
  PTR<Object<Poly2>> f;
  PTR<Object<Poly2>> g;
  PTR<Object<PV2>> p;
  PTR<Object<PV2>> q;
  DeclareCalculate(PV2);

 public:
  LinearizationIntersection(PTR<Object<Poly2>> f, PTR<Object<Poly2>> g,
                            PTR<Object<PV2>> p)
      : f(f), g(g), p(p), q(nullptr) {}
  LinearizationIntersection(PTR<Object<Poly2>> f, PTR<Object<Poly2>> g,
                            PTR<Object<PV2>> p, PTR<Object<PV2>> q)
      : f(f), g(g), p(p), q(q) {}
};

class LinearizationPoly : public Object<Poly2> {
  PTR<Object<Poly2>> f;
  PTR<Object<PV2>> p;

  DeclareCalculate(Poly2);

 public:
  LinearizationPoly(PTR<Object<Poly2>> f, PTR<Object<PV2>> p) : f(f), p(p) {}

  PTR<Line> getLine();
};

class LinearDistanceEstimate : public Object<Scalar> {
  PTR<Object<Poly2>> f;
  PTR<Object<PV2>> p;
  DeclareCalculate(Scalar) {
    PV2<N> pp = p->get<N>();
    Poly2<N> ff = f->get<N>();
    PV2<N> gradP(ff.derX().value(&(pp.x)), ff.derY().value(&(pp.x)));
    return ff.value(&(pp.x)).abs() / gradP.length();
  }

 public:
  LinearDistanceEstimate(PTR<Object<Poly2>> f, PTR<Object<PV2>> p)
      : f(f), p(p) {}
};

class EllipseCriticalPoint : public Object<PV2> {
  PTR<Object<Poly2>> e;
  DeclareCalculate(PV2);

 public:
  EllipseCriticalPoint(PTR<Object<Poly2>> e) : e(e) {}
};

class Line : public RefCnt {
  friend class Edge;

  template <class N>
  N leftNeg(PTR<Object<PV2>> q) {
    return (q->get<N>() - p->get<N>()).cross(v->get<N>());
  }

 public:
  PTR<Object<PV2>> p;
  PTR<Object<PV2>> v;
  int axis;  //-1 for not aligned, 0 for horiz, 1 for vert

  Line(PTR<Object<PV2>> p, PTR<Object<PV2>> v, int axis = -1)
      : p(p), v(v), axis(axis) {}
  template <class N>
  PV2<N> getP() {
    return p->get<N>();
  }
  template <class N>
  PV2<N> getV() {
    return v->get<N>();
  }
  template <class N>
  PV2<N> getP(PTR<Object<Scalar>> t) {
    return getP(t->get<N>().x);
  }
  template <class N>
  PV2<N> getP(const N& t) {
    return p->get<N>() + v->get<N>() * t;
  }
  int getAxis() const { return axis; }
  template <class N>
  N hit(PTR<Line> l) {
    // p + v t lies on l=(p',v')
    // (p + v t - p') x v' = 0
    // v x v' t = (p' - p) x v'
    // t = (p' - p) x v' / v x v';
    return (l->getP<N>() - getP<N>()).cross(l->getV<N>()) /
           getV<N>().cross(l->getV<N>());
  }

  class distance_object : public Object<Scalar> {
   public:
    distance_object(PTR<Object<PV2>> op, PTR<Object<PV2>> ov,
                    PTR<Object<PV2>> oq)
        : p(op), v(ov), q(oq) {}

    DeclareCalculate(Scalar) {
      PV2<N> pv = p->get<N>();
      PV2<N> vv = v->get<N>();
      PV2<N> qv = q->get<N>();
      PV2<N> pq = qv - pv;
      PV2<N> e = pq - pq.dot(vv) * vv;

      return e.length();
    }

   private:
    PTR<Object<PV2>> p, v, q;
  };

  PTR<Object<Scalar>> distance(PTR<Object<PV2>> q) {
    return new distance_object(p, v, q);
  }

  PTR<Line> approximate(PTR<Object<PV2>> a, PTR<Object<PV2>> b);

  Primitive2(Left, PTR<Line>, l, PTR<Object<PV2>>, q)

      int left(PTR<Object<PV2>> q) {
    return Left(this, q);
  }
};

class LineAB : public Line {
 public:
  LineAB(PTR<Object<PV2>> a, PTR<Object<PV2>> b, int axis = -1)
      : Line(a, new VectorAB(a, b), axis) {}
  // LineAB has to delete v because it creates it:
  //~LineAB () { delete dynamic_cast<VectorAB*>(v); }
};

class LLHit : public Object<Scalar> {
  PTR<Line> l1;
  PTR<Line> l2;
  DeclareCalculate(Scalar) { return l1->hit<N>(l2); }

 public:
  LLHit(PTR<Line> l1, PTR<Line> l2) : l1(l1), l2(l2) {}
};

class LinePoint : public Object<PV2> {
  PTR<Line> l;
  PTR<Object<Scalar>> t;
  DeclareCalculate(PV2) { return l->getP<N>(t); }

 public:
  LinePoint(PTR<Line> l, PTR<Object<Scalar>> t) : l(l), t(t) {}
};

// To intersect a line with an Ellipse, create a LinePoly and call
// getRoots().  For each root, create a LinePoint from the line and
// that root.  If you are only interested in roots in a cell, you
// should first calculate the the entrance and exit lines and compare
// the root to the two LLHit objects using a Different predicate.
class LinePoly : public Object<Poly> {
  PTR<Line> l;
  PTR<Object<Poly2>> f;
  DeclareCalculate(Poly);  // f(p + v t)

 public:
  LinePoly(PTR<Line> l, PTR<Object<Poly2>> f) : l(l), f(f) {}
  static bool print_all;
};

// Poly that shifts out the root at x=0 before returning the poly
class ShiftOutPoly : public Object<Poly> {
  PTR<Object<Poly>> p;
  DeclareCalculate(Poly);

 public:
  ShiftOutPoly(PTR<Object<Poly>> p) : p(p) {}
};

// Poly defined as f(p + t*(q-p)) - g(p + t*(q-p)) = 0
// To get the point for which to base a splitting line
// create a middle poly and call getRoots(). There should
// be at least one root. Choose the middle root if multiple.
class MiddlePoly : public Object<Poly> {
  PTR<Object<Poly2>> f;
  PTR<Object<Poly2>> g;
  PTR<Object<PV2>> p;
  PTR<Object<PV2>> q;
  DeclareCalculate(Poly);  // f(p + t*(q-p)) - g(p + t*(q-p)) = 0
 public:
  MiddlePoly(PTR<Object<Poly2>> f, PTR<Object<Poly2>> g, PTR<Object<PV2>> p,
             PTR<Object<PV2>> q)
      : f(f), g(g), p(p), q(q) {}
};

#ifdef BLEEN
class SubXPoly : public Object<Poly> {
  PTR<Object<Poly2>> f;
  PTR<Object<Scalar>> x;
  DeclareCalculate(Poly) { return f->get<N>().subX(x->get<N>()); }

 public:
  SubXPoly(PTR<Object<Poly2>> f, PTR<Object<Scalar>> x) : f(f), x(x) {}
};

class SubYPoly : public Object<Poly> {
  PTR<Object<Poly2>> f;
  PTR<Object<Scalar>> y;
  DeclareCalculate(Poly) { return f->get<N>().subY(y->get<N>()); }

 public:
  SubYPoly(PTR<Object<Poly2>> f, PTR<Object<Scalar>> y) : f(f), y(y) {}
};
#endif

class DerPoly : public Object<Poly> {
  PTR<Object<Poly>> f;
  DeclareCalculate(Poly) { return f->get<N>().der(); }

 public:
  DerPoly(PTR<Object<Poly>> f) : f(f) {}
};

class GradPoint : public Object<PV2> {
  PTR<Object<Poly2>> f;
  PTR<Object<PV2>> p;

  DeclareCalculate(PV2) {
    Poly2<N> fp = f->get<N>();
    PV2<N> pp = p->get<N>();
    return PV2<N>(fp.derX().value(pp.x, pp.y), fp.derY().value(pp.x, pp.y));
  }

 public:
  GradPoint(PTR<Object<Poly2>> f, PTR<Object<PV2>> p) : f(f), p(p) {}
  bool same_poly(PTR<GradPoint> that) { return f == that->f; }
};

template <class N>
Poly<N> substitute(PTR<Line> l, PTR<Object<Poly2>> f);

template <class N>
Poly<N> interpolate(vector<N> b);

template <class N>
Poly<N> interpolate(PTR<Line> l, PTR<Object<Poly2>> f);

PTR<Object<Scalar>> approxAverage(double t, PTR<Object<Scalar>> a,
                                  PTR<Object<Scalar>> b);

// f(p + (q - p) t) =approx= g(p + (q - p) t)
PTR<Object<Scalar>> middleT(PTR<Object<PV2>> p, PTR<Object<Poly2>> f,
                            PTR<Object<PV2>> q, PTR<Object<Poly2>> g);

// Approximate average of tangent vectors of f and g at p.
PTR<Object<PV2>> middleV(PTR<Object<Poly2>> f, PTR<Object<Poly2>> g,
                         PTR<Object<PV2>> p, PTR<Object<PV2>> a,
                         PTR<Object<PV2>> b);

// To create a splitting line
// t = new InputScalar(middleT(p, f, q, g))
// l = new LineAB(p, q);
// p = new LinePoint(l, t)
// v = new InputPoint(middleV(f, g, p))
// split = new Line(p, v)

//--Primitives-------------------------------------------------------

Primitive2(Intersects, PTR<Object<PV2>>, a, PTR<Object<PV2>>, b);

// sign() == -1 if counterclockwise
Primitive3(Orient2D, PTR<Object<PV2>>, a, PTR<Object<PV2>>, b, PTR<Object<PV2>>,
           c);
template <class N>
inline N Orient2D::calculate() {
  PV2<N> ap = a->get<N>();
  PV2<N> bp = b->get<N>();
  PV2<N> cp = c->get<N>();

  return (cp - ap).cross(bp - ap);
}

// sign() == -1 if o is left of p + t*v
Primitive3(LeftOf, PTR<Object<PV2>>, p, PTR<Object<PV2>>, v, PTR<Object<PV2>>,
           o);

// sign() == -1 if a.y < b.y
Primitive2(YOrder, PTR<Object<PV2>>, a, PTR<Object<PV2>>, b);
template <class N>
inline N YOrder::calculate() {
  return (a->get<N>().y - b->get<N>().y);
}

// sign() == -1 if a.x < b.x
Primitive2(XOrder, PTR<Object<PV2>>, a, PTR<Object<PV2>>, b);
template <class N>
inline N XOrder::calculate() {
  return (a->get<N>().x - b->get<N>().x);
}

// sign() == -1 if a < b
Primitive2(LessThan, PTR<Object<Scalar>>, a, PTR<Object<Scalar>>, b);

// sign() == -1 if |a| < |b|
Primitive2(LessThanAbs, PTR<Object<Scalar>>, a, PTR<Object<Scalar>>, b);

// sign() of a * b
Primitive2(Dot, PTR<Object<PV2>>, a, PTR<Object<PV2>>, b);
template <class N>
inline N Dot::calculate() {
  return (a->get<N>().dot(b->get<N>()));
}

Primitive2(Side, PTR<Object<PV2>>, a, PTR<Object<PV2>>, b);

Primitive2(Poly1Sign, PTR<Object<Poly>>, f, PTR<Object<Scalar>>, p);
template <class N>
inline N Poly1Sign::calculate() {
  Poly<N> ff = f->get<N>();
  N pp = p->get<N>();

  int i = -1;
  i++;

  return ff.value(pp);
}

Primitive2(PolySign, PTR<Object<Poly2>>, f, PTR<Object<PV2>>, p);

Primitive3(Descartes2, PTR<Object<Poly2>>, f, PTR<Object<PV2>>, l,
           PTR<Object<PV2>>, u);

Primitive1(Sign, PTR<Object<Scalar>>, a);
template <class N>
inline N Sign::calculate() {
  return a->get<N>();
}

Primitive4(Shorter, PTR<Object<PV2>>, a, PTR<Object<PV2>>, b, PTR<Object<PV2>>,
           c, PTR<Object<PV2>>, d);
template <class N>
inline N Shorter::calculate() {
  N l1 = (a->get<N>() - b->get<N>()).length();
  N l2 = (c->get<N>() - d->get<N>()).length();
  // return l1 < l2 ? -1 : 1;
  return l1 - l2;
}

Primitive2(CCWOrder, PTR<Object<PV2>>, a, PTR<Object<PV2>>, b);
template <class N>
inline N CCWOrder::calculate() {
  PV2<N> ap = a->get<N>();
  PV2<N> bp = b->get<N>();

  PV2<N> x = PV2<N>(N(1), N(0));

  int sa = x.cross(ap).sign();
  int sb = x.cross(bp).sign();

  if (sa * sb > 0) return ap.cross(bp);

  // REFACTOR
  // return sa;
  return x.cross(ap);
}

class Ascending {
 public:
  bool operator()(PTR<Object<Scalar>> const a, PTR<Object<Scalar>> const b) {
    return LessThan(a, b) == -1;
  }
};

class Descending {
 public:
  bool operator()(PTR<Object<Scalar>> const a, PTR<Object<Scalar>> const b) {
    return LessThan(a, b) == 1;
  }
};

class CCW {
 public:
  bool operator()(PTR<Object<PV2>> const a, PTR<Object<PV2>> const b) {
    return CCWOrder(a, b) == 1;
  }
};

#endif
