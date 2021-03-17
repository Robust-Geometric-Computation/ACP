#include "acp/resultant/resultant_arrangement.h"
#include "acp/poly/root.h"

#include <chrono>

bool ResultantArrangement::alternates(Parameter& yp, Parameter& x_f,
                                      Parameter& x_g, PTR<Object<Poly2>> f,
                                      PTR<Object<Poly2>> g) {
  PTR<Object<Scalar>> yl = new InputParameter(yp.lbP());
  PTR<Object<Scalar>> yh = new InputParameter(yp.ubP());

  // get intersection with bottom of box
  PTR<Object<Poly>> f_bot = new SubstitutePoly1(yl, f, 1);
  PTR<Object<Poly>> g_bot = new SubstitutePoly1(yl, g, 1);

  PTR<Object<Scalar>> f_root = new PolyRoot(
      f_bot, new Object<Scalar>(x_f.lbP()), new Object<Scalar>(x_f.ubP()));
  PTR<Object<Scalar>> g_root = new PolyRoot(
      g_bot, new Object<Scalar>(x_g.lbP()), new Object<Scalar>(x_g.ubP()));

  // TEST TEST TEST

  PTR<Object<Scalar>> xl = new InputParameter(x_f.lbP());
  PTR<Object<Scalar>> xh = new InputParameter(x_f.ubP());

  vector<PTR<Object<Scalar>>> f_bot_r = PolySolver(f_bot).getRoots(xl, xh);

  xl = new InputParameter(x_g.lbP());
  xh = new InputParameter(x_g.ubP());

  vector<PTR<Object<Scalar>>> g_bot_r = PolySolver(g_bot).getRoots(xl, xh);

  assert(f_bot_r.size() == 1);
  assert(g_bot_r.size() == 1);

  // TEST TEST TEST

  bool f_to_g_bot = (LessThan(f_root, g_root) == -1);
  // bool f_to_g_bot = (LessThan(f_bot_r[0], g_bot_r[0]) == -1);

  // get intersection with top of box
  f_bot = new SubstitutePoly1(yh, f, 1);
  g_bot = new SubstitutePoly1(yh, g, 1);

  f_root = new PolyRoot(f_bot, new Object<Scalar>(x_f.lbP()),
                        new Object<Scalar>(x_f.ubP()));
  g_root = new PolyRoot(g_bot, new Object<Scalar>(x_g.lbP()),
                        new Object<Scalar>(x_g.ubP()));

  // TEST TEST TEST

  xl = new InputParameter(x_f.lbP());
  xh = new InputParameter(x_f.ubP());

  f_bot_r = PolySolver(f_bot).getRoots(xl, xh);

  xl = new InputParameter(x_g.lbP());
  xh = new InputParameter(x_g.ubP());

  g_bot_r = PolySolver(g_bot).getRoots(xl, xh);

  assert(f_bot_r.size() == 1);
  assert(g_bot_r.size() == 1);

  // TEST TEST TEST

  bool f_to_g_top = (LessThan(f_root, g_root) == -1);
  // bool f_to_g_top = (LessThan(f_bot_r[0], g_bot_r[0]) == -1);

  return f_to_g_bot != f_to_g_top;
}

void ResultantArrangement::alternates_double(vector<PV2<Parameter>>& boxes,
                                             PTR<Object<Poly2>> f,
                                             PTR<Object<Poly2>> g,
                                             vector<PV2<Parameter>>& verified) {
  // for each box, find the boxes with the shortest distance in y and in x take
  // half of this as the expansion size for the new larger box

  struct root_pair {
    PTR<Object<Scalar>> r;
    int i;
  };

  for (int i = 0; i < boxes.size(); i++) {
    int j = (i + 1) % boxes.size();

    // start with distance to the boundary as min for box expansion
    double xp =
        fmin((boxes[i].x - g_xl).abs().lb(), (boxes[i].x - g_xh).abs().lb());
    double yp =
        fmin((boxes[i].y - g_yl).abs().lb(), (boxes[i].y - g_yh).abs().lb());

    xp = fmin(xp, yp);
    yp = fmin(xp, yp);

    double d;

    while (j != i) {
      d = (boxes[i] - boxes[j]).length().lb();
      yp = fmin(d, yp);
      xp = fmin(d, xp);

      j = (j + 1) % boxes.size();
    }

    // go half of the shortest distance in either direction
    xp /= 2;
    yp /= 2;

    PTR<Object<Scalar>> xl =
        new InputParameter(Parameter::constant(boxes[i].x.lb() - xp));
    PTR<Object<Scalar>> xh =
        new InputParameter(Parameter::constant(boxes[i].x.ub() + xp));
    PTR<Object<Scalar>> yl =
        new InputParameter(Parameter::constant(boxes[i].y.lb() - yp));
    PTR<Object<Scalar>> yh =
        new InputParameter(Parameter::constant(boxes[i].y.ub() + yp));

    PTR<Object<Poly>> f_bot = new SubstitutePoly1(yl, f, 1);
    PTR<Object<Poly>> g_bot = new SubstitutePoly1(yl, g, 1);
    PTR<Object<Poly>> f_right = new SubstitutePoly1(xh, f, 0);
    PTR<Object<Poly>> g_right = new SubstitutePoly1(xh, g, 0);
    PTR<Object<Poly>> f_top = new SubstitutePoly1(yh, f, 1);
    PTR<Object<Poly>> g_top = new SubstitutePoly1(yh, g, 1);
    PTR<Object<Poly>> f_left = new SubstitutePoly1(xl, f, 0);
    PTR<Object<Poly>> g_left = new SubstitutePoly1(xl, g, 0);

    vector<PTR<Object<Scalar>>> f_bot_r = PolySolver(f_bot).getRoots(xl, xh);
    vector<PTR<Object<Scalar>>> g_bot_r = PolySolver(g_bot).getRoots(xl, xh);
    vector<PTR<Object<Scalar>>> f_right_r =
        PolySolver(f_right).getRoots(yl, yh);
    vector<PTR<Object<Scalar>>> g_right_r =
        PolySolver(g_right).getRoots(yl, yh);
    vector<PTR<Object<Scalar>>> f_top_r = PolySolver(f_top).getRoots(xl, xh);
    vector<PTR<Object<Scalar>>> g_top_r = PolySolver(g_top).getRoots(xl, xh);
    vector<PTR<Object<Scalar>>> f_left_r = PolySolver(f_left).getRoots(yl, yh);
    vector<PTR<Object<Scalar>>> g_left_r = PolySolver(g_left).getRoots(yl, yh);

    reverse(f_top_r.begin(), f_top_r.end());
    reverse(g_top_r.begin(), g_top_r.end());
    reverse(f_left_r.begin(), f_left_r.end());
    reverse(g_left_r.begin(), g_left_r.end());

    vector<struct root_pair> roots;

    // merge the list of roots on each side into a single list, ordered CCW

    int k = j = 0;

    while (j < f_bot_r.size() && k < g_bot_r.size()) {
      if (LessThan(f_bot_r[j], g_bot_r[k]) == -1) {
        roots.push_back({.r = f_bot_r[j++], .i = 0});
      } else {
        roots.push_back({.r = g_bot_r[k++], .i = 1});
      }
    }

    while (j < f_bot_r.size()) {
      roots.push_back({.r = f_bot_r[j++], .i = 0});
    }

    while (k < g_bot_r.size()) {
      roots.push_back({.r = g_bot_r[k++], .i = 1});
    }

    j = 0;
    k = 0;
    while (j < f_right_r.size() && k < g_right_r.size()) {
      if (LessThan(f_right_r[j], g_right_r[k]) == -1) {
        roots.push_back({.r = f_right_r[j++], .i = 0});
      } else {
        roots.push_back({.r = g_right_r[k++], .i = 1});
      }
    }

    while (j < f_right_r.size()) {
      roots.push_back({.r = f_right_r[j++], .i = 0});
    }

    while (k < g_right_r.size()) {
      roots.push_back({.r = g_right_r[k++], .i = 1});
    }

    // top and left are ordered descending because reversed
    // use LessThan == 1 to ensure larger first ordering on these sides

    j = 0;
    k = 0;
    while (j < f_top_r.size() && k < g_top_r.size()) {
      if (LessThan(f_top_r[j], g_top_r[k]) == 1) {
        roots.push_back({.r = f_top_r[j++], .i = 0});
      } else {
        roots.push_back({.r = g_top_r[k++], .i = 1});
      }
    }

    while (j < f_top_r.size()) {
      roots.push_back({.r = f_top_r[j++], .i = 0});
    }

    while (k < g_top_r.size()) {
      roots.push_back({.r = g_top_r[k++], .i = 1});
    }

    j = 0;
    k = 0;
    while (j < f_left_r.size() && k < g_left_r.size()) {
      if (LessThan(f_left_r[j], g_left_r[k]) == 1) {
        roots.push_back({.r = f_left_r[j++], .i = 0});
      } else {
        roots.push_back({.r = g_left_r[k++], .i = 1});
      }
    }

    while (j < f_left_r.size()) {
      roots.push_back({.r = f_left_r[j++], .i = 0});
    }

    while (k < g_left_r.size()) {
      roots.push_back({.r = g_left_r[k++], .i = 1});
    }

    // now calculate the parity of the # of intersections by using parity along
    // the border
    // https://books.google.com/books?id=KNYQmJf0V0cC&pg=PA50&lpg=PA50&dq=parity+of+curve+intersections&source=bl&ots=2HkYueE_QU&sig=kMQAJ-hTxlGLFZt0LLufDQ3A0Lc&hl=en&sa=X&ved=0ahUKEwitrsffu6jOAhUT6GMKHVdiCk8Q6AEIMzAE#v=onepage&q=parity%20of%20curve%20intersections&f=false
    // Weighted Ham-Sandwich Cuts (Prosenjit Bose and Stefan Langerman)

    // find the first f intersection
    k = j = 0;
    while (j < roots.size()) {
      if (roots[j].i == 0) break;
      j++;
    }

    // we need to find at least one f intersetion
    assert(j != roots.size());

    int count = 0;
    k = (j + 1) % roots.size();

    vector<int> par;

    // compute x1, x2, x3, x4, ...
    // count the number of g intersections inbetween pairs of f intersections
    while (k != j) {
      if (roots[k].i == 1) {
        count++;
      } else {
        par.push_back(count);
        count = 0;
      }

      k = (k + 1) % roots.size();
    }

    // don't forget the last pair
    par.push_back(count);

    // parity of intersections == parity of sum of odd indexed xj's
    int parity = 0;
    for (j = 1; j < par.size(); j += 2) {
      parity += par[j];
    }

    // if parity is 1 than there must be a single intersection in this larger
    // box containing a root
    if (parity % 2) {
      verified.push_back(boxes[i]);
    }
  }
}

vector<PTR<Object<PV2>>> ResultantArrangement::getIntersections(
    PTR<Object<Poly2>> f, PTR<Object<Poly2>> g) {
  PTR<Object<Scalar>> xl = new InputParameter(g_xl);
  PTR<Object<Scalar>> xh = new InputParameter(g_xh);
  PTR<Object<Scalar>> yl = new InputParameter(g_yl);
  PTR<Object<Scalar>> yh = new InputParameter(g_xh);

  // auto t1 = std::chrono::high_resolution_clock::now();

  // cout << "(resultant calcs): start" << endl;

  // solve for the resultant
  PTR<Object<Poly>> res_pp = new Resultant(f, g);
  Poly<Parameter> ress = res_pp->getApprox(1.0);
  ress.print();
  PTR<Object<Poly>> resultant = res_pp;
  // PTR<Object<Poly>> resultant_der = new DerPoly(res_pp);
  // PTR<Object<Poly>> resultant_qr = new Resultant(f, g);

  // int i = 0;

  // int phsign = Poly1Sign(resultant, yh);
  // int plsign = Poly1Sign(resultant, yl);

  // phsign = Poly1Sign(resultant, yh);
  // plsign = Poly1Sign(resultant, yl);

  // Poly r_e = resultant->getPoly();

  // Resultant::useEuclid = false;
  //
  // phsign = Poly1Sign(resultant_qr, yh);
  // plsign = Poly1Sign(resultant_qr, yl);

  // Poly r_q = resultant_qr->getPoly();

  // Resultant::useEuclid = true;
  ////Poly pp = resultant->getPoly();

  // Parameter inp = Parameter::constant(0.5);

  // auto t2 = std::chrono::high_resolution_clock::now();

  // cout << "(resultant calcs): " <<
  // std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count() * .001
  // << " ms\n";

  auto t1 = std::chrono::high_resolution_clock::now();

  // cout << "(resultant roots): start" << endl;

  // calculate the possible intersection intervals for y
  int k = 0;
  vector<PTR<Object<Scalar>>> y_int = PolySolver(resultant).getRoots(yl, yh);
  // vector<PTR<Object<Scalar>>> y_der_int =
  // PolySolver(resultant_der).getRoots(yl, yh);

  vector<PTR<Object<PV2>>> roots;

  // take the intersection of the resulting inteval roots
  vector<Parameter> intersection_f;
  vector<Parameter> intersection_g;

  vector<Parameter> yps;

  for (int i = 0; i < y_int.size(); i++) {
    PTR<Object<Scalar>> y = y_int[i];
    Parameter yp = y->getApprox(1.0);

    // substitute each candidate y interval into each poly
    PTR<Object<Poly>> f_sub = new SubstitutePoly1(y, f, 1);
    PTR<Object<Poly>> g_sub = new SubstitutePoly1(y, g, 1);

    // get the roots of each substitution and intersect those intervals
    vector<PTR<Object<Scalar>>> f_rootsR = PolySolver(f_sub).getRoots(xl, xh);
    vector<PTR<Object<Scalar>>> g_rootsR = PolySolver(g_sub).getRoots(xl, xh);

    vector<PTR<Object<Scalar>>> f_roots(f_rootsR.size());
    vector<PTR<Object<Scalar>>> g_roots(g_rootsR.size());

    for (int i = 0; i < f_roots.size(); i++) f_roots[i] = f_rootsR[i];
    for (int i = 0; i < g_roots.size(); i++) g_roots[i] = g_rootsR[i];

    sort(f_roots.begin(), f_roots.end(), Ascending());

    sort(g_roots.begin(), g_roots.end(), Ascending());

    int j = 0;
    int k = 0;

    while (j < f_roots.size() && k < g_roots.size()) {
      Parameter a = f_roots[j]->getApprox(1.0);
      Parameter b = g_roots[k]->getApprox(1.0);

      if (a.intersects(b)) {
        intersection_f.push_back(a);
        intersection_g.push_back(b);
        yps.push_back(yp);
      }

      // if(LessThan(f_roots[j], g_roots[k]) == -1) j++;
      if (a.ub() < b.ub())
        j++;
      else
        k++;
    }
  }

  bool high_precision_solver = false;

  if (high_precision_solver) {
    //(x_i, y) is a possible root, construct it and ensure f and g alternate
    for (int j = 0; j < intersection_f.size(); j++) {
      Parameter yp = yps[j];
      if (alternates(yp, intersection_f[j], intersection_g[j], f, g)) {
        PV2<Parameter> root(intersection_f[j].intersect(intersection_g[j]), yp);
        roots.push_back(new Root2(f, g, root));
      }
    }

  } else {
    // construct all of the intersection boxes
    vector<PV2<Parameter>> boxes;
    for (int i = 0; i < intersection_f.size(); i++) {
      boxes.push_back(PV2<Parameter>(
          intersection_f[i].intersect(intersection_g[i]), yps[i]));
    }

    vector<PV2<Parameter>> verified;

    alternates_double(boxes, f, g, verified);

    for (int i = 0; i < verified.size(); i++) {
      roots.push_back(new Root2(f, g, verified[i]));
    }
  }

  auto t2 = std::chrono::high_resolution_clock::now();

  // cout << "(resultant roots): " <<
  // std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count() * .001
  // << " ms\n";

  return roots;
}
