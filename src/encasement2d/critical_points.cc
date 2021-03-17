#include "acp/encasement2d/critical_points.h"
#include "acp/encasement2d/plot.h"
#include "acp/encasement2d/predicates2d.h"
#include "acp/poly/root.h"

vector<Rectangle*> Rectangle::rects;
int Rectangle::global_id = 0;

vector<Poly2<Parameter>> Rectangle::test_poly;

// With path compression
Rectangle* Rectangle::find() {
  if (this->parent == this)
    return this;
  else
    return (this->parent = this->parent->find());
}

// Union by rank
void Rectangle::unionize(Rectangle* that) {
  Rectangle* p1 = find();
  Rectangle* p2 = that->find();

  if (p1->rank < p2->rank) {
    p1->parent = p2;
    p2->rank = max(p1->rank + 1, p2->rank);
  } else {
    p2->parent = p1;
    p1->rank = max(p2->rank + 1, p1->rank);
  }
}

// Set the sign of the polynomials on this rect
//(f, fx, fy)
// 0   1  2
void Rectangle::set_sign(vector<PTR<Object<Poly2>>> v) {
  PTR<Object<PV2>> xy =
      new ParameterPoint(new DoubleInputInterval(box.x.lb(), box.x.ub()),
                         new DoubleInputInterval(box.y.lb(), box.y.ub()));

  PTR<Object<PV2>> l = new ParameterPoint(new InputScalar(box.x.lb()),
                                          new InputScalar(box.y.lb()));
  PTR<Object<PV2>> u = new ParameterPoint(new InputScalar(box.x.ub()),
                                          new InputScalar(box.y.ub()));

  this->sign_f = (new PolyScalar(v[0], xy))->getApprox(1.0).x.sign(false);
  this->sign_fx = (new PolyScalar(v[1], xy))->getApprox(1.0).x.sign(false);
  this->sign_fy = (new PolyScalar(v[2], xy))->getApprox(1.0).x.sign(false);
  // this->sign_hess = (new PolyScalar(v[3], xy))->get().sign(false);

  // if(this->sign_hess < 0) return;

  if (this->sign_f == 0 && this->sign_fx == 0 && this->sign_fy == 0) {
    // this->sign_f  = cubic_approx_sign(xy, v[0]);
    // this->sign_fx = cubic_approx_sign(xy, v[1]);
    // this->sign_fy = cubic_approx_sign(xy, v[2]);

    // this->sign_f   = (v[0]->get().d < 9 ? quad_approx_sign(xy, v[0]) :
    // cubic_approx_sign(xy, v[0])); this->sign_fx  = (v[1]->get().d < 9 ?
    // quad_approx_sign(xy, v[1]) : cubic_approx_sign(xy, v[1])); this->sign_fy
    // = (v[2]->get().d < 9 ? quad_approx_sign(xy, v[2]) : cubic_approx_sign(xy,
    // v[2]));

    // this->sign_f    = (v[0]->get().d < 9 ? quad_approx_sign(xy, v[0]) :
    // n_approx_sign(xy, v[0], v[0]->get().d / 2)); this->sign_fx   =
    // (v[1]->get().d < 9 ? quad_approx_sign(xy, v[1]) : n_approx_sign(xy, v[1],
    // v[1]->get().d / 2)); this->sign_fy   = (v[2]->get().d < 9 ?
    // quad_approx_sign(xy, v[2]) : n_approx_sign(xy, v[2], v[2]->get().d / 2));

    this->sign_f = (Descartes2(v[0], l, u) == 0 ? 1 : 0);
    this->sign_fx = (Descartes2(v[1], l, u) == 0 ? 1 : 0);
    this->sign_fy = (Descartes2(v[2], l, u) == 0 ? 1 : 0);

    // this->sign_hess = (v[3]->get().d < 9 ? quad_approx_sign(xy, v[3]) :
    // n_approx_sign(xy, v[3], v[3]->get().d / 2));

    // static bool printed = false;
    // if(!printed) {
    // printf("cubic_approx\n");
    // printed = true;
    //}

    // this->sign_f  = n_approx_sign(xy, v[0], 4);
    // this->sign_fx = n_approx_sign(xy, v[1], 4);
    // this->sign_fy = n_approx_sign(xy, v[2], 4);
    //	  static bool printed = false;
    //	  if(!printed) {
    //		printf("n_approx\n");
    //		printed = true;
    //	  }
  }

  // PV2 xy_p = xy->get();

  // this->sign_f =  v[0]->get().value(&xy_p.x).sign(false);
  // this->sign_fx = v[1]->get().value(&xy_p.x).sign(false);
  // this->sign_fy = v[2]->get().value(&xy_p.x).sign(false);
}

// Split a rectangle along major axis
Rectangle* Rectangle::split() {
  if (fabs(box.x.ub() - box.x.lb()) < fabs(box.y.ub() - box.y.lb())) {
    PV2<Parameter> box1(box.x, Parameter::interval(box.y.lb(), box.y.mid()));
    PV2<Parameter> box2(box.x, Parameter::interval(box.y.mid(), box.y.ub()));
    this->box = box1;
    this->num_splits++;
    return new Rectangle(box2, this->num_splits);
  } else {
    PV2<Parameter> box1(Parameter::interval(box.x.lb(), box.x.mid()), box.y);
    PV2<Parameter> box2(Parameter::interval(box.x.mid(), box.x.ub()), box.y);
    this->box = box1;
    this->num_splits++;
    return new Rectangle(box2, this->num_splits);
  }
}

int quad_approx_sign(PTR<Object<PV2>> xy, PTR<Object<Poly2>> f) {
  // get the four corners
  PTR<Object<PV2>> ll = new CornerPoint(xy, 0, 0);
  PTR<Object<PV2>> lr = new CornerPoint(xy, 1, 0);
  PTR<Object<PV2>> ul = new CornerPoint(xy, 0, 1);
  PTR<Object<PV2>> ur = new CornerPoint(xy, 1, 1);

  // evaluate the poly on the four corners
  PTR<Object<Scalar>> fll = new PolyScalar(f, ll);
  PTR<Object<Scalar>> flr = new PolyScalar(f, lr);
  PTR<Object<Scalar>> ful = new PolyScalar(f, ul);
  PTR<Object<Scalar>> fur = new PolyScalar(f, ur);

  // calculate the sign of f on the corners
  int sll = Sign(fll);
  int slr = Sign(flr);
  int sul = Sign(ful);
  int sur = Sign(fur);

  int sum = sll + slr + sul + sur;

  // if it's not consistent across all corners, fail
  if (abs(sum) != 4) return 0;

  if (f->getApprox(1.0).d < 2) return sll;

  PTR<Object<Poly2>> quad;

  // calculate the 2nd order Taylor approximation of f
  // function of (x - a) and (y - b)
  if (sll > 0) {
    quad = new QuadApproxPoly2(f, xy, 0);  // lower bound quad approx
  } else {
    quad = new QuadApproxPoly2(f, xy, 1);  // upper bound and negated approx
  }

  // int xysign = (new PolyScalar(quad, xy))->get().sign(false);
  // if(xysign != 0) return xysign;

  // Rectangle::test_poly.push_back(quad->get());

  // now the function is translated to (0, 0) to (w, h)
  PTR<Object<Scalar>> w = new VectorLength(new VectorAB(ll, lr));
  PTR<Object<Scalar>> h = new VectorLength(new VectorAB(ll, ul));

  PTR<Object<Scalar>> ZERO = new InputScalar(Parameter::constant(0));
  PTR<Object<Scalar>> ONE = new InputScalar(Parameter::constant(1));

  // translated box coordinates
  PTR<Object<PV2>> ll_z = new ParameterPoint(ZERO, ZERO);
  PTR<Object<PV2>> lr_z = new ParameterPoint(w, ZERO);
  PTR<Object<PV2>> ur_z = new ParameterPoint(w, h);
  PTR<Object<PV2>> ul_z = new ParameterPoint(ZERO, h);

  PTR<Object<Scalar>> qll = new PolyScalar(quad, ll_z);
  PTR<Object<Scalar>> qlr = new PolyScalar(quad, lr_z);
  PTR<Object<Scalar>> qur = new PolyScalar(quad, ur_z);
  PTR<Object<Scalar>> qul = new PolyScalar(quad, ul_z);

  int sll_z = Sign(qll);
  int slr_z = Sign(qlr);
  int sur_z = Sign(qur);
  int sul_z = Sign(qul);

  int qsum = sll_z + slr_z + sur_z + sul_z;

  if (qsum < 4) return 0;

  // get the polys whose roots might be intersections with the boundary
  // PTR<Object<Poly>> bot   = new InputPoly(new InputPoly1(substitute(new
  // LineAB(ll_z, lr_z), quad))); PTR<Object<Poly>> right = new InputPoly(new
  // InputPoly1(substitute(new LineAB(lr_z, ur_z), quad))); PTR<Object<Poly>>
  // top = new InputPoly(new InputPoly1(substitute(new LineAB(ul_z, ur_z),
  // quad))); PTR<Object<Poly>> left  = new InputPoly(new
  // InputPoly1(substitute(new LineAB(ll_z, ul_z), quad)));

  PTR<Object<Poly>> bot = new LinePoly(new LineAB(ll_z, lr_z), quad);
  PTR<Object<Poly>> right = new LinePoly(new LineAB(lr_z, ur_z), quad);
  PTR<Object<Poly>> top = new LinePoly(new LineAB(ul_z, ur_z), quad);
  PTR<Object<Poly>> left = new LinePoly(new LineAB(ll_z, ul_z), quad);

  // roots of the polys
  vector<PTR<Object<Scalar>>> bot_r = PolySolver(bot).getRoots(ZERO, ONE);
  vector<PTR<Object<Scalar>>> right_r = PolySolver(right).getRoots(ZERO, ONE);
  vector<PTR<Object<Scalar>>> top_r = PolySolver(top).getRoots(ZERO, ONE);
  vector<PTR<Object<Scalar>>> left_r = PolySolver(left).getRoots(ZERO, ONE);

  // if the quadratic intersects the boundary at all, fail
  if (bot_r.size() > 0) return 0;
  if (right_r.size() > 0) return 0;
  if (top_r.size() > 0) return 0;
  if (left_r.size() > 0) return 0;

  // Find the critical point
  PTR<Object<PV2>> crit = new QuadCriticalPoint(quad);

  PTR<Object<Scalar>> px = new PVScalar(crit, 0);
  PTR<Object<Scalar>> py = new PVScalar(crit, 1);

  PTR<Object<Scalar>> bx = new PVScalar(xy, 0);
  PTR<Object<Scalar>> by = new PVScalar(xy, 1);

  PTR<Object<Scalar>> box_x = new ParameterInterval(ZERO, w);
  PTR<Object<Scalar>> box_y = new ParameterInterval(ZERO, h);

  PTR<Object<PV2>> box_xy = new ParameterPoint(box_x, box_y);

  // if the critical point is in the box, evaluate f on the critical point
  if (Intersects(crit, box_xy) == 1) {
    PTR<Object<Scalar>> f_crit = new PolyScalar(quad, crit);

    // if the signs differ, there is a loop
    if (Sign(f_crit) > 0)
      return sll;
    else
      return 0;

  }
  // if the critical point is outside of the box, f is constant
  else {
    return sll;
  }
}

int cubic_approx_sign(PTR<Object<PV2>> xy, PTR<Object<Poly2>> f) {
  // get the four corners
  PTR<Object<PV2>> ll = new CornerPoint(xy, 0, 0);
  PTR<Object<PV2>> lr = new CornerPoint(xy, 1, 0);
  PTR<Object<PV2>> ul = new CornerPoint(xy, 0, 1);
  PTR<Object<PV2>> ur = new CornerPoint(xy, 1, 1);

  // evaluate the poly on the four corners
  PTR<Object<Scalar>> fll = new PolyScalar(f, ll);
  PTR<Object<Scalar>> flr = new PolyScalar(f, lr);
  PTR<Object<Scalar>> ful = new PolyScalar(f, ul);
  PTR<Object<Scalar>> fur = new PolyScalar(f, ur);

  // calculate the sign of f on the corners
  int sll = Sign(fll);
  int slr = Sign(flr);
  int sul = Sign(ful);
  int sur = Sign(fur);

  int sum = sll + slr + sul + sur;

  // if it's not consistent across all corners, fail
  if (abs(sum) != 4) return 0;

  if (f->getApprox(1.0).d < 2) return sll;

  PTR<Object<Poly2>> cubic;

  // calculate the 3rd order Taylor approximation of f
  // function of (x - a) and (y - b)
  if (sll > 0) {
    cubic = new CubicApproxPoly2(f, xy, 0);  // lower bound cubic approx
  } else {
    cubic = new CubicApproxPoly2(f, xy, 1);  // upper bound and negated approx
  }

  // Rectangle::test_poly.push_back(quad->get());

  // now the function is translated to (0, 0) to (w, h)
  PTR<Object<Scalar>> w = new VectorLength(new VectorAB(ll, lr));
  PTR<Object<Scalar>> h = new VectorLength(new VectorAB(ll, ul));

  PTR<Object<Scalar>> ZERO = new InputScalar(Parameter::constant(0));
  PTR<Object<Scalar>> ONE = new InputScalar(Parameter::constant(1));

  // translated box coordinates
  PTR<Object<PV2>> ll_z = new ParameterPoint(ZERO, ZERO);
  PTR<Object<PV2>> lr_z = new ParameterPoint(w, ZERO);
  PTR<Object<PV2>> ur_z = new ParameterPoint(w, h);
  PTR<Object<PV2>> ul_z = new ParameterPoint(ZERO, h);

  PTR<Object<Scalar>> qll = new PolyScalar(cubic, ll_z);
  PTR<Object<Scalar>> qlr = new PolyScalar(cubic, lr_z);
  PTR<Object<Scalar>> qur = new PolyScalar(cubic, ur_z);
  PTR<Object<Scalar>> qul = new PolyScalar(cubic, ul_z);

  int sll_z = Sign(qll);
  int slr_z = Sign(qlr);
  int sur_z = Sign(qur);
  int sul_z = Sign(qul);

  int qsum = sll_z + slr_z + sur_z + sul_z;

  if (qsum < 4) return 0;

  // get the polys whose roots might be intersections with the boundary
  // PTR<Object<Poly>> bot   = new InputPoly(new InputPoly1(substitute(new
  // LineAB(ll_z, lr_z), quad))); PTR<Object<Poly>> right = new InputPoly(new
  // InputPoly1(substitute(new LineAB(lr_z, ur_z), quad))); PTR<Object<Poly>>
  // top = new InputPoly(new InputPoly1(substitute(new LineAB(ul_z, ur_z),
  // quad))); PTR<Object<Poly>> left  = new InputPoly(new
  // InputPoly1(substitute(new LineAB(ll_z, ul_z), quad)));

  PTR<Object<Poly>> bot = new LinePoly(new LineAB(ll_z, lr_z), cubic);
  PTR<Object<Poly>> right = new LinePoly(new LineAB(lr_z, ur_z), cubic);
  PTR<Object<Poly>> top = new LinePoly(new LineAB(ul_z, ur_z), cubic);
  PTR<Object<Poly>> left = new LinePoly(new LineAB(ll_z, ul_z), cubic);

  // roots of the polys
  vector<PTR<Object<Scalar>>> bot_r = PolySolver(bot).getRoots(ZERO, ONE);
  vector<PTR<Object<Scalar>>> right_r = PolySolver(right).getRoots(ZERO, ONE);
  vector<PTR<Object<Scalar>>> top_r = PolySolver(top).getRoots(ZERO, ONE);
  vector<PTR<Object<Scalar>>> left_r = PolySolver(left).getRoots(ZERO, ONE);

  // if the quadratic intersects the boundary at all, fail
  if (bot_r.size() > 0) return 0;
  if (right_r.size() > 0) return 0;
  if (top_r.size() > 0) return 0;
  if (left_r.size() > 0) return 0;

  // Find the critical points

  PTR<Object<Scalar>> box_x = new ParameterInterval(ZERO, w);
  PTR<Object<Scalar>> box_y = new ParameterInterval(ZERO, h);

  PTR<Object<PV2>> box_xy = new ParameterPoint(box_x, box_y);

  PTR<Object<Poly2>> cubic_x = new DerXPoly2(cubic);
  PTR<Object<Poly2>> cubic_y = new DerYPoly2(cubic);

  vector<PTR<Object<PV2>>> roots;

  if (false) {
    PTR<Object<PV2>> bot_left = new ParameterPoint(ZERO, ZERO);
    PTR<Object<PV2>> top_right = new ParameterPoint(w, h);

    Encasement* enc = new Encasement(bot_left, top_right, 0, false);

    vector<PTR<Object<Poly2>>> curves;
    curves.push_back(cubic_x);
    curves.push_back(cubic_y);

    enc->init();

    enc->calculate(curves);

    roots = enc->getRoots();

  } else {
    ResultantArrangement result;

    roots = result.getIntersections(cubic_x, cubic_y);
  }

  for (int i = 0; i < roots.size(); i++) {
    if (Intersects(roots[i], box_xy) == 1) {
      PTR<Object<Scalar>> f_crit = new PolyScalar(cubic, roots[i]);
      if (Sign(f_crit) < 0) return 0;
    }
  }

  return sll;
}
int n_approx_sign(PTR<Object<PV2>> xy, PTR<Object<Poly2>> f, int n) {
  // get the four corners
  PTR<Object<PV2>> ll = new CornerPoint(xy, 0, 0);
  PTR<Object<PV2>> lr = new CornerPoint(xy, 1, 0);
  PTR<Object<PV2>> ul = new CornerPoint(xy, 0, 1);
  PTR<Object<PV2>> ur = new CornerPoint(xy, 1, 1);

  // evaluate the poly on the four corners
  PTR<Object<Scalar>> fll = new PolyScalar(f, ll);
  PTR<Object<Scalar>> flr = new PolyScalar(f, lr);
  PTR<Object<Scalar>> ful = new PolyScalar(f, ul);
  PTR<Object<Scalar>> fur = new PolyScalar(f, ur);

  // calculate the sign of f on the corners
  int sll = Sign(fll);
  int slr = Sign(flr);
  int sul = Sign(ful);
  int sur = Sign(fur);

  int sum = sll + slr + sul + sur;

  // if it's not consistent across all corners, fail
  if (abs(sum) != 4) return 0;

  if (f->getApprox(1.0).d < 2) return sll;

  PTR<Object<Poly2>> cubic;

  // calculate the 3rd order Taylor approximation of f
  // function of (x - a) and (y - b)
  if (sll > 0) {
    cubic = new NApproxPoly2(f, xy, 0, n);  // lower bound cubic approx
  } else {
    cubic = new NApproxPoly2(f, xy, 1, n);  // upper bound and negated approx
  }

  // Rectangle::test_poly.push_back(quad->get());

  // now the function is translated to (0, 0) to (w, h)
  PTR<Object<Scalar>> w = new VectorLength(new VectorAB(ll, lr));
  PTR<Object<Scalar>> h = new VectorLength(new VectorAB(ll, ul));

  PTR<Object<Scalar>> ZERO = new InputScalar(Parameter::constant(0));
  PTR<Object<Scalar>> ONE = new InputScalar(Parameter::constant(1));

  // translated box coordinates
  PTR<Object<PV2>> ll_z = new ParameterPoint(ZERO, ZERO);
  PTR<Object<PV2>> lr_z = new ParameterPoint(w, ZERO);
  PTR<Object<PV2>> ur_z = new ParameterPoint(w, h);
  PTR<Object<PV2>> ul_z = new ParameterPoint(ZERO, h);

  PTR<Object<Scalar>> qll = new PolyScalar(cubic, ll_z);
  PTR<Object<Scalar>> qlr = new PolyScalar(cubic, lr_z);
  PTR<Object<Scalar>> qur = new PolyScalar(cubic, ur_z);
  PTR<Object<Scalar>> qul = new PolyScalar(cubic, ul_z);

  int sll_z = Sign(qll);
  int slr_z = Sign(qlr);
  int sur_z = Sign(qur);
  int sul_z = Sign(qul);

  int qsum = sll_z + slr_z + sur_z + sul_z;

  if (qsum < 4) return 0;

  // get the polys whose roots might be intersections with the boundary
  // PTR<Object<Poly>> bot   = new InputPoly(new InputPoly1(substitute(new
  // LineAB(ll_z, lr_z), quad))); PTR<Object<Poly>> right = new InputPoly(new
  // InputPoly1(substitute(new LineAB(lr_z, ur_z), quad))); PTR<Object<Poly>>
  // top = new InputPoly(new InputPoly1(substitute(new LineAB(ul_z, ur_z),
  // quad))); PTR<Object<Poly>> left  = new InputPoly(new
  // InputPoly1(substitute(new LineAB(ll_z, ul_z), quad)));

  PTR<Object<Poly>> bot = new LinePoly(new LineAB(ll_z, lr_z), cubic);
  PTR<Object<Poly>> right = new LinePoly(new LineAB(lr_z, ur_z), cubic);
  PTR<Object<Poly>> top = new LinePoly(new LineAB(ul_z, ur_z), cubic);
  PTR<Object<Poly>> left = new LinePoly(new LineAB(ll_z, ul_z), cubic);

  // roots of the polys
  vector<PTR<Object<Scalar>>> bot_r = PolySolver(bot).getRoots(ZERO, ONE);
  vector<PTR<Object<Scalar>>> right_r = PolySolver(right).getRoots(ZERO, ONE);
  vector<PTR<Object<Scalar>>> top_r = PolySolver(top).getRoots(ZERO, ONE);
  vector<PTR<Object<Scalar>>> left_r = PolySolver(left).getRoots(ZERO, ONE);

  // if the quadratic intersects the boundary at all, fail
  if (bot_r.size() > 0) return 0;
  if (right_r.size() > 0) return 0;
  if (top_r.size() > 0) return 0;
  if (left_r.size() > 0) return 0;

  // Find the critical points

  PTR<Object<Poly2>> cubic_x = new DerXPoly2(cubic);
  PTR<Object<Poly2>> cubic_y = new DerYPoly2(cubic);

  PTR<Object<Scalar>> box_x = new ParameterInterval(ZERO, w);
  PTR<Object<Scalar>> box_y = new ParameterInterval(ZERO, h);

  PTR<Object<PV2>> box_xy = new ParameterPoint(box_x, box_y);

  vector<PTR<Object<PV2>>> roots;

  if (true) {
    PTR<Object<PV2>> bot_left = new ParameterPoint(ZERO, ZERO);
    PTR<Object<PV2>> top_right = new ParameterPoint(w, h);

    Encasement* enc = new Encasement(bot_left, top_right, 0, false);

    vector<PTR<Object<Poly2>>> curves;
    curves.push_back(cubic_x);
    curves.push_back(cubic_y);

    enc->init();

    enc->calculate(curves);

    roots = enc->getRoots();

  } else {
    ResultantArrangement result;

    roots = result.getIntersections(cubic_x, cubic_y);
  }

  for (int i = 0; i < roots.size(); i++) {
    if (Intersects(roots[i], box_xy) == 1) {
      PTR<Object<Scalar>> f_crit = new PolyScalar(cubic, roots[i]);
      if (Sign(f_crit) < 0) return 0;
    }
  }

  return sll;
}

void animate(vector<Rectangle*> good, queue<Rectangle*> q,
             Rectangle* bad = nullptr) {
  DrawItem item;
  vector<PTR<Object<PV2>>> draw_points;

  int s = q.size();

  vector<Rectangle*> qv;

  while (s-- > 0) {
    Rectangle* r = q.front();
    q.pop();
    qv.push_back(r);
    q.push(r);
  }

  if (bad) {
    PTR<Object<PV2>> rbox = new ParameterPoint(
        new DoubleInputInterval(bad->box.x.lb(), bad->box.x.ub()),
        new DoubleInputInterval(bad->box.y.lb(), bad->box.y.ub()));

    PTR<Object<PV2>> ll = new CornerPoint(rbox, 0, 0);
    PTR<Object<PV2>> lr = new CornerPoint(rbox, 1, 0);
    PTR<Object<PV2>> ul = new CornerPoint(rbox, 0, 1);
    PTR<Object<PV2>> ur = new CornerPoint(rbox, 1, 1);

    item.types.push_back(TRI);
    item.colors.push_back({1.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(ll);
    draw_points.push_back(lr);
    draw_points.push_back(ul);
    item.pvs.push_back(draw_points);

    item.types.push_back(TRI);
    item.colors.push_back({1.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(lr);
    draw_points.push_back(ur);
    draw_points.push_back(ul);
    item.pvs.push_back(draw_points);

    item.types.push_back(LINE);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(ll);
    draw_points.push_back(lr);
    item.pvs.push_back(draw_points);

    item.types.push_back(LINE);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(lr);
    draw_points.push_back(ur);
    item.pvs.push_back(draw_points);

    item.types.push_back(LINE);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(ur);
    draw_points.push_back(ul);
    item.pvs.push_back(draw_points);

    item.types.push_back(LINE);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(ul);
    draw_points.push_back(ll);
    item.pvs.push_back(draw_points);
  }

  for (vector<Rectangle*>::iterator it = good.begin(); it != good.end(); it++) {
    PTR<Object<PV2>> rbox = new ParameterPoint(
        new DoubleInputInterval((*it)->box.x.lb(), (*it)->box.x.ub()),
        new DoubleInputInterval((*it)->box.y.lb(), (*it)->box.y.ub()));

    PTR<Object<PV2>> ll = new CornerPoint(rbox, 0, 0);
    PTR<Object<PV2>> lr = new CornerPoint(rbox, 1, 0);
    PTR<Object<PV2>> ul = new CornerPoint(rbox, 0, 1);
    PTR<Object<PV2>> ur = new CornerPoint(rbox, 1, 1);

    item.types.push_back(TRI);
    item.colors.push_back({0.0f, 1.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(ll);
    draw_points.push_back(lr);
    draw_points.push_back(ul);
    item.pvs.push_back(draw_points);

    item.types.push_back(TRI);
    item.colors.push_back({0.0f, 1.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(lr);
    draw_points.push_back(ur);
    draw_points.push_back(ul);
    item.pvs.push_back(draw_points);

    item.types.push_back(LINE);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(ll);
    draw_points.push_back(lr);
    item.pvs.push_back(draw_points);

    item.types.push_back(LINE);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(lr);
    draw_points.push_back(ur);
    item.pvs.push_back(draw_points);

    item.types.push_back(LINE);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(ur);
    draw_points.push_back(ul);
    item.pvs.push_back(draw_points);

    item.types.push_back(LINE);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(ul);
    draw_points.push_back(ll);
    item.pvs.push_back(draw_points);
  }

  for (vector<Rectangle*>::iterator it = qv.begin(); it != qv.end(); it++) {
    PTR<Object<PV2>> rbox = new ParameterPoint(
        new DoubleInputInterval((*it)->box.x.lb(), (*it)->box.x.ub()),
        new DoubleInputInterval((*it)->box.y.lb(), (*it)->box.y.ub()));

    PTR<Object<PV2>> ll = new CornerPoint(rbox, 0, 0);
    PTR<Object<PV2>> lr = new CornerPoint(rbox, 1, 0);
    PTR<Object<PV2>> ul = new CornerPoint(rbox, 0, 1);
    PTR<Object<PV2>> ur = new CornerPoint(rbox, 1, 1);

    item.types.push_back(TRI);
    item.colors.push_back({0.3f, 0.3f, 0.3f});
    draw_points.clear();
    draw_points.push_back(ll);
    draw_points.push_back(lr);
    draw_points.push_back(ul);
    item.pvs.push_back(draw_points);

    item.types.push_back(TRI);
    item.colors.push_back({0.3f, 0.3f, 0.3f});
    draw_points.clear();
    draw_points.push_back(lr);
    draw_points.push_back(ur);
    draw_points.push_back(ul);
    item.pvs.push_back(draw_points);

    item.types.push_back(LINE);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(ll);
    draw_points.push_back(lr);
    item.pvs.push_back(draw_points);

    item.types.push_back(LINE);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(lr);
    draw_points.push_back(ur);
    item.pvs.push_back(draw_points);

    item.types.push_back(LINE);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(ur);
    draw_points.push_back(ul);
    item.pvs.push_back(draw_points);

    item.types.push_back(LINE);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(ul);
    draw_points.push_back(ll);
    item.pvs.push_back(draw_points);
  }

  // debugDrawCrit(item);
}

// Splits the starting rectangle until sign of f is constant and fx, fy
// are still unknown
vector<Rectangle*> critical_points(PTR<Object<Poly2>> f, Rectangle* start) {
  vector<Rectangle*> regions;
  queue<Rectangle*> q;

  q.push(start);

  Rectangle* s = new Rectangle(start->box);
  Rectangle* ss = new Rectangle(start->box);

  vector<PTR<Object<Poly2>>> v;

  v.push_back(f);                       // 0 f
  v.push_back(new DerXPoly2(v[0]));     // 1 fx
  v.push_back(new DerYPoly2(v[0]));     // 2 fy
  v.push_back(new HessianPoly2(v[0]));  // 2 fy

  vector<int> counts;

  while (!q.empty()) {
    // animate(regions, q);

    Rectangle* r = q.front();
    q.pop();

    // printf("[%f %f][%f %f] ", r->box.x.lb(), r->box.x.ub(), r->box.y.lb(),
    // r->box.y.ub()); fflush(stdin);

    r->set_sign(v);

    // printf("sign: %d %d %d\n", r->sign_f, r->sign_fx, r->sign_fy);

    // if(r->sign_hess >= 0 && r->sign_fx == 0 && r->sign_fy == 0) {
    if (r->sign_fx == 0 && r->sign_fy == 0) {
      if (r->sign_f != 0) {
        regions.push_back(r);
        int depth = r->num_splits;
        while (counts.size() < depth + 1) counts.push_back(0);
        counts[depth] = counts[depth] + 1;
      } else {
        // animate(regions, q, r);

        Rectangle* rr = new Rectangle(r->box);

        q.push(r->split());
        q.push(r);

        double amin = fmin(rr->box.x.ub() - rr->box.x.lb(),
                           rr->box.y.ub() - rr->box.y.lb());

        double smin =
            fmin(s->box.x.ub() - s->box.x.lb(), s->box.y.ub() - s->box.y.lb());

        if (amin < smin) {
          ss = s;
          s = rr;
          ss->set_sign(v);
          s->set_sign(v);
        }
      }
    } else {
      int depth = r->num_splits;
      while (counts.size() < depth + 1) counts.push_back(0);
      counts[depth] = counts[depth] + 1;
    }
  }

  //  printf("\nsmallest rectangle considered to split:\n");
  //  printf("[%.16f, %.16f] [%.16f, %.16f] w: %.16f h: %.16f\n\n",
  //  s->box.x.lb(), s->box.x.ub(), s->box.y.lb(), s->box.y.ub(), s->box.x.ub()
  //  - s->box.x.lb(), s->box.y.ub() - s->box.y.lb());
  //
  //  printf("\nsecond smallest rectangle considered to split:\n");
  //  printf("[%.16f, %.16f] [%.16f, %.16f] w: %.16f h: %.16f\n\n",
  //  ss->box.x.lb(), ss->box.x.ub(), ss->box.y.lb(), ss->box.y.ub(),
  //  ss->box.x.ub() - ss->box.x.lb(), ss->box.y.ub() - ss->box.y.lb());
  //
  //  printf("Number of terminating rects at depth:\n");
  //
  //  for(int i = 0; i < counts.size(); i++) {
  //    printf("  %d: %d\n", i, counts[i]);
  //  }

  // Parameter fs = f->get().value(&s->box.x);
  // Parameter fxs = f->get().derX().value(&s->box.x);
  // Parameter fys = f->get().derY().value(&s->box.x);

  // Parameter fss = f->get().value(&ss->box.x);
  // Parameter fxss = f->get().derX().value(&ss->box.x);
  // Parameter fyss = f->get().derY().value(&ss->box.x);

  // printf("f(s): [%g, %g]\n\n", fs.lb(), fs.ub());
  // printf("fx(s): [%g, %g]\n\n", fxs.lb(), fxs.ub());
  // printf("fy(s): [%g, %g]\n\n", fys.lb(), fys.ub());
  // printf("f(ss): [%g, %g]\n\n", fss.lb(), fss.ub());
  // printf("fx(ss): [%g, %g]\n\n", fxss.lb(), fxss.ub());
  // printf("fy(ss): [%g, %g]\n\n", fyss.lb(), fyss.ub());

  return regions;
}

// Returns set of rects representing all unique connected components
map<int, set<Rectangle*>> connected_components(vector<Rectangle*> v) {
  map<double, set<Rectangle*, CompareByLBY>> lbxMap;
  map<double, set<Rectangle*, CompareByLBY>> ubxMap;
  map<double, set<Rectangle*, CompareByLBX>> lbyMap;
  map<double, set<Rectangle*, CompareByLBX>> ubyMap;

  for (int i = 0; i < v.size(); i++) {
    lbxMap[v[i]->lbx()].insert(v[i]);
    ubxMap[v[i]->ubx()].insert(v[i]);
    lbyMap[v[i]->lby()].insert(v[i]);
    ubyMap[v[i]->uby()].insert(v[i]);
  }

  for (map<double, set<Rectangle*, CompareByLBY>>::iterator x = lbxMap.begin();
       x != lbxMap.end(); x++) {
    double xx = (*x).first;
    set<Rectangle*, CompareByLBY>::iterator ita = lbxMap[xx].begin();
    set<Rectangle*, CompareByLBY>::iterator itb = ubxMap[xx].begin();

    // printf("xx: %g\n", xx);

    while (ita != lbxMap[xx].end() && itb != ubxMap[xx].end()) {
      // printf("a lby: %g\n", (*ita)->lby());
      // printf("a uby: %g\n", (*ita)->uby());
      // printf("b lby: %g\n", (*itb)->lby());
      // printf("b uby: %g\n", (*itb)->uby());

      if ((*ita)->uby() < (*itb)->lby())
        ++ita;
      else if ((*ita)->lby() > (*itb)->uby())
        ++itb;
      else {
        Rectangle* r1 = *ita;
        Rectangle* r2 = *itb;

        if (r1->find() != r2->find()) {
          // printf("    join\n");
          r1->unionize(r2);
        } else {
          // printf("    already joined\n");
        }

        if (r1->uby() < r2->uby())
          ++ita;
        else
          ++itb;
      }
    }
  }

  for (map<double, set<Rectangle*, CompareByLBX>>::iterator y = lbyMap.begin();
       y != lbyMap.end(); y++) {
    double yy = (*y).first;
    set<Rectangle*, CompareByLBX>::iterator ita = lbyMap[yy].begin();
    set<Rectangle*, CompareByLBX>::iterator itb = ubyMap[yy].begin();

    // printf("yy: %g\n", yy);

    while (ita != lbyMap[yy].end() && itb != ubyMap[yy].end()) {
      // printf("a lbx: %g\n", (*ita)->lbx());
      // printf("a ubx: %g\n", (*ita)->ubx());
      // printf("b lbx: %g\n", (*itb)->lbx());
      // printf("b ubx: %g\n", (*itb)->ubx());

      if ((*ita)->ubx() < (*itb)->lbx())
        ++ita;
      else if ((*ita)->lbx() > (*itb)->ubx())
        ++itb;
      else {
        Rectangle* r1 = *ita;
        Rectangle* r2 = *itb;

        if (r1->find() != r2->find()) {
          // printf("    join\n");
          r1->unionize(r2);
        } else {
          // printf("    already joined\n");
        }

        if (r1->ubx() < r2->ubx())
          ++ita;
        else
          ++itb;
      }
    }
  }

  map<int, set<Rectangle*>> s;

  for (int i = 0; i < v.size(); i++) {
    s[v[i]->find()->id].insert(v[i]);
  }

  return s;
}
