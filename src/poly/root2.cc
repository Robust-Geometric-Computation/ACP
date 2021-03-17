/*
  ACP (Adaptive Controlled Precision) Library
  for robust computational geometry

  Copyright (c) 2013-07-15
  Victor Joseph Milenkovic
  University of Miami
  vjm@miami.edu
  Elisha Sacks
  Purdue University
  eps@cs.purdue.edu

  This file contains code described in

Robust Complete Path Planning in the Plane
Victor Milenkovic, Elisha Sacks, and Steven Trac
Proceedings of the Workshop on the Algorithmic Foundations of Robotics (WAFR)
pages 37-52, 2012

   This code is under development.
   It is free for use for academic and research purposes.
   Permission from the authors is required for commercial use.
*/

#include "acp/poly/root2.h"
#include "acp/encasement2d/predicates2d.h"
#include "acp/poly/root.h"

using namespace std;
using namespace acp;

namespace acp {

bool Root2::polishing = false;

template <class N>
PV2<N> Root2::polish(PV2<Parameter> p) {
  /*
  vector<PV2> b;
  b.push_back(p);
  vector< PTR<Object<PPoly2>> > func;
  func.push_back(new Object<PPoly2>(ff));
  func.push_back(new Object<PPoly2>(gg));
  */
  // debugDraw(b, func);

  polishing = true;

  // printf("------------------------------------------------------------------------------------------\n\n");

  RootBoundary<N> rb = {.I = p, .v = vector<vector<PTR<Object<Scalar>>>>()};

  // run Newton's until it fails then try to subdivide
  // but if higher precision, just run Newton's
  if (curPrecision() == 53u) {
    do {
      newton<N>(rb);
    } while (subdivide<N>(rb));

  } else {
    newton<N>(rb);
  }

  polishing = false;

  // printf("I [%.16lf, %.16lf] [%.16lf, %.16lf]\n", rb.I.x.lb(), rb.I.x.ub(),
  // rb.I.y.lb(), rb.I.y.ub());

  // printf("------------------------------------------------------------------------------------------\n\n");

  return rb.I;
}

// sign() == -1 if a < b
// Primitive2(LessThan, Object<Parameter>* , a, Object<Parameter>* , b);
// int LessThan::sign() {
//  return (a->get() - b->get()).sign();
//}

template <class N>
vector<int> Root2::intersections(vector<vector<PTR<Object<Scalar>>>>& rbv) {
  vector<vector<PTR<Object<Scalar>>>> fv;
  vector<vector<PTR<Object<Scalar>>>> gv;

  for (int i = 0; i < 4; i++) fv.push_back(rbv[i]);
  for (int i = 4; i < 8; i++) gv.push_back(rbv[i]);

  vector<int> alt;

  for (int k = 0; k < 4; k++) {
    vector<PTR<Object<Scalar>>> frr = fv[k];
    vector<PTR<Object<Scalar>>> grr = gv[k];

    for (int i = 0, j = 0; i < frr.size() || j < grr.size();) {
      PTR<Object<Scalar>> fr = (i < frr.size() ? frr[i] : nullptr);
      PTR<Object<Scalar>> gr = (j < grr.size() ? grr[j] : nullptr);

      // k < 2 is bottom and right which are ascending roots, k >= 2 are top and
      // left which are descending
      if (fr != nullptr &&
          (gr == nullptr || (k < 2 ? fr->get<N>().x < gr->get<N>().x
                                   : gr->get<N>().x < fr->get<N>().x))) {
        alt.push_back(0);
        i++;
      } else {
        alt.push_back(1);
        j++;
      }
    }
  }

  return alt;
}

int Root2::parity(const vector<int>& alt) {
  // find the first f intersection
  int k = 0;
  int j = 0;
  while (j < alt.size()) {
    if (alt[j] == 0) break;
    j++;
  }

  // if f doesn't intersect this cell, keep going
  if (j == alt.size()) {
    return 0;
  }

  int count = 0;
  k = (j + 1) % alt.size();

  vector<int> par;

  // compute x1, x2, x3, x4, ...
  // count the number of g intersections inbetween pairs of f intersections
  while (k != j) {
    if (alt[k] == 1) {
      count++;
    } else {
      par.push_back(count);
      count = 0;
    }

    k = (k + 1) % alt.size();
  }

  // don't forget the last pair
  par.push_back(count);

  // parity of intersections == parity of sum of odd indexed xj's
  int parity = 0;
  for (j = 1; j < par.size(); j += 2) {
    parity += par[j];
  }

  return parity % 2;
}

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

template <class N>
void Root2::setRoots(RootBoundary<N>& rb) {
  PV2<N> b = rb.I;

  rb.v.clear();

  assert(curPrecision() == 53u);
  PTR<Object<Scalar>> xl = new Object<Scalar>(Parameter(b.x.lbP()));
  PTR<Object<Scalar>> xu = new Object<Scalar>(Parameter(b.x.ubP()));
  PTR<Object<Scalar>> yl = new Object<Scalar>(Parameter(b.y.lbP()));
  PTR<Object<Scalar>> yu = new Object<Scalar>(Parameter(b.y.ubP()));

  PTR<Object<Poly>> bot_f = new SubYPoly(f, yl);
  PTR<Object<Poly>> right_f = new SubXPoly(f, xu);
  PTR<Object<Poly>> top_f = new SubYPoly(f, yu);
  PTR<Object<Poly>> left_f = new SubXPoly(f, xl);

  rb.v.push_back(PolySolver(bot_f).getRoots(xl, xu));
  rb.v.push_back(PolySolver(right_f).getRoots(yl, yu));
  rb.v.push_back(PolySolver(top_f).getRoots(xl, xu));
  rb.v.push_back(PolySolver(left_f).getRoots(yl, yu));

  PTR<Object<Poly>> bot_g = new SubYPoly(g, yl);
  PTR<Object<Poly>> right_g = new SubXPoly(g, xu);
  PTR<Object<Poly>> top_g = new SubYPoly(g, yu);
  PTR<Object<Poly>> left_g = new SubXPoly(g, xl);

  rb.v.push_back(PolySolver(bot_g).getRoots(xl, xu));
  rb.v.push_back(PolySolver(right_g).getRoots(yl, yu));
  rb.v.push_back(PolySolver(top_g).getRoots(xl, xu));
  rb.v.push_back(PolySolver(left_g).getRoots(yl, yu));

  // top and left need to be reversed for parity check
  reverse(rb.v[2].begin(), rb.v[2].end());
  reverse(rb.v[3].begin(), rb.v[3].end());
  reverse(rb.v[6].begin(), rb.v[6].end());
  reverse(rb.v[7].begin(), rb.v[7].end());
}

template <class N>
Root2::RootBoundary<N> Root2::splitHoriz(RootBoundary<N>& rb, N s) {
  vector<vector<PTR<Object<Scalar>>>> bot(8);
  vector<vector<PTR<Object<Scalar>>>> top(8);

  // bot side for bot f and g roots are same
  bot[0] = rb.v[0];
  bot[4] = rb.v[4];

  // top side for top f and g roots are the same
  top[2] = rb.v[2];
  top[6] = rb.v[6];

  vector<vector<PTR<Object<Scalar>>>>& v = rb.v;

  // Separate roots around the split point, give left to left and right to right
  // for f [0, 4) and g [4, 8)
  for (int i = 0; i < 8; i += 4) {
    int li = 0;
    while (li < rb.v[i + 1].size() && v[i + 1][li]->get<N>().x < s) {
      li++;
    }
    //[0, li) are < s [li, size) are >= s

    int ui = 0;
    while (ui < rb.v[i + 3].size() && s < v[i + 3][ui]->get<N>().x) {
      ui++;
    }
    //[0, ui) are > s [ui, size) are <= s

    bot[i + 1] = vector<PTR<Object<Scalar>>>(rb.v[i + 1].begin(),
                                             rb.v[i + 1].begin() + li);
    top[i + 1] = vector<PTR<Object<Scalar>>>(rb.v[i + 1].begin() + li,
                                             rb.v[i + 1].end());

    top[i + 3] = vector<PTR<Object<Scalar>>>(rb.v[i + 3].begin(),
                                             rb.v[i + 3].begin() + ui);
    bot[i + 3] = vector<PTR<Object<Scalar>>>(rb.v[i + 3].begin() + ui,
                                             rb.v[i + 3].end());
  }

  PTR<Object<Scalar>> xl = new Object<Scalar>(Scalar<N>(rb.I.x.lbP()));
  PTR<Object<Scalar>> xu = new Object<Scalar>(Scalar<N>(rb.I.x.ubP()));

  PTR<Object<Poly>> mid_line_f =
      new SubYPoly(f, new Object<Scalar>(Scalar<N>(s)));
  PTR<Object<Poly>> mid_line_g =
      new SubYPoly(g, new Object<Scalar>(Scalar<N>(s)));

  vector<PTR<Object<Scalar>>> mid_roots_f_asc =
      PolySolver(mid_line_f).getRoots(xl, xu);
  vector<PTR<Object<Scalar>>> mid_roots_g_asc =
      PolySolver(mid_line_g).getRoots(xl, xu);

  vector<PTR<Object<Scalar>>> mid_roots_f_desc(mid_roots_f_asc.begin(),
                                               mid_roots_f_asc.end());
  vector<PTR<Object<Scalar>>> mid_roots_g_desc(mid_roots_g_asc.begin(),
                                               mid_roots_g_asc.end());

  reverse(mid_roots_f_desc.begin(), mid_roots_f_desc.end());
  reverse(mid_roots_g_desc.begin(), mid_roots_g_desc.end());

  top[0] = mid_roots_f_asc;
  top[4] = mid_roots_g_asc;

  bot[2] = mid_roots_f_desc;
  bot[6] = mid_roots_g_desc;

  PV2<N> I = rb.I;

  N yl = I.y.lbP();
  N yu = I.y.ubP();

  rb = {.I = PV2<N>(I.x, s.interval(yu)), .v = top};

  return {.I = PV2<N>(I.x, yl.interval(s)), .v = bot};
}

template <class N>
Root2::RootBoundary<N> Root2::splitVert(RootBoundary<N>& rb, N s) {
  vector<vector<PTR<Object<Scalar>>>> left(8);
  vector<vector<PTR<Object<Scalar>>>> right(8);

  // left side for left f and g roots are same
  left[3] = rb.v[3];
  left[7] = rb.v[7];

  // right side for right f and g roots are the same
  right[1] = rb.v[1];
  right[5] = rb.v[5];

  vector<vector<PTR<Object<Scalar>>>>& v = rb.v;

  // Separate roots around the split point, give left to left and right to right
  // for f [0, 4) and g [4, 8)
  for (int i = 0; i < 8; i += 4) {
    int li = 0;
    while (li < rb.v[i + 0].size() && v[i + 0][li]->get<N>().x < s) {
      li++;
    }
    //[0, li) are < s [li, size) are >= s

    int ui = 0;
    while (ui < rb.v[i + 2].size() && s < v[i + 2][ui]->get<N>().x) {
      ui++;
    }
    //[0, ui) are > s [ui, size) are <= s

    left[i + 0] = vector<PTR<Object<Scalar>>>(rb.v[i + 0].begin(),
                                              rb.v[i + 0].begin() + li);
    right[i + 0] = vector<PTR<Object<Scalar>>>(rb.v[i + 0].begin() + li,
                                               rb.v[i + 0].end());

    right[i + 2] = vector<PTR<Object<Scalar>>>(rb.v[i + 2].begin(),
                                               rb.v[i + 2].begin() + ui);
    left[i + 2] = vector<PTR<Object<Scalar>>>(rb.v[i + 2].begin() + ui,
                                              rb.v[i + 2].end());
  }

  assert(curPrecision() == 53u);
  PTR<Object<Scalar>> yl = new Object<Scalar>(Scalar<N>(rb.I.y.lbP()));
  PTR<Object<Scalar>> yu = new Object<Scalar>(Scalar<N>(rb.I.y.ubP()));

  PTR<Object<Poly>> mid_line_f =
      new SubXPoly(f, new Object<Scalar>(Scalar<N>(s)));
  PTR<Object<Poly>> mid_line_g =
      new SubXPoly(g, new Object<Scalar>(Scalar<N>(s)));

  vector<PTR<Object<Scalar>>> mid_roots_f_asc =
      PolySolver(mid_line_f).getRoots(yl, yu);
  vector<PTR<Object<Scalar>>> mid_roots_g_asc =
      PolySolver(mid_line_g).getRoots(yl, yu);

  vector<PTR<Object<Scalar>>> mid_roots_f_desc(mid_roots_f_asc.begin(),
                                               mid_roots_f_asc.end());
  vector<PTR<Object<Scalar>>> mid_roots_g_desc(mid_roots_g_asc.begin(),
                                               mid_roots_g_asc.end());

  reverse(mid_roots_f_desc.begin(), mid_roots_f_desc.end());
  reverse(mid_roots_g_desc.begin(), mid_roots_g_desc.end());

  left[1] = mid_roots_f_asc;
  left[5] = mid_roots_g_asc;

  right[3] = mid_roots_f_desc;
  right[7] = mid_roots_g_desc;

  PV2<N> I = rb.I;

  N xl = I.x.lbP();
  N xu = I.x.ubP();

  rb = {.I = PV2<N>(s.interval(xu), I.y), .v = right};

  return {.I = PV2<N>(xl.interval(s), I.y), .v = left};
  ;
  ;
  ;
}

// subdivide the cell once along major axis
template <class N>
bool Root2::subdivide(RootBoundary<N>& rb) {
  PV2<N> I = rb.I;

  // printf("SUBDIVIDE\n\n");

  int even_count = 0;
  int odd_count = 0;
  int odd_index = -1;

  vector<RootBoundary<N>> b;

  try {
    if (rb.v.size() == 0) {
      setRoots(rb);
    }

    if (I.x.intervalWidth() > I.y.intervalWidth()) {
      N xl = I.x.lbP();
      N xu = I.x.ubP();
      N xl75 = (0.75 * xl + 0.25 * xu).midP();
      N xl25 = (0.25 * xl + 0.75 * xu).midP();

      b.push_back(splitVert(rb, xl75));
      b.push_back(splitVert(rb, xl25));
      b.push_back(rb);

    } else {
      N yl = I.y.lbP();
      N yu = I.y.ubP();
      N yl75 = (0.75 * yl + 0.25 * yu).midP();
      N yl25 = (0.25 * yl + 0.75 * yu).midP();

      b.push_back(splitHoriz(rb, yl75));
      b.push_back(splitHoriz(rb, yl25));
      b.push_back(rb);
    }

    for (int it = 0; it < b.size(); it++) {
      // vector<PV2> list;
      // list.push_back(bb);

      // vector< PTR<Object<Poly2>> > fs;
      // fs.push_back(new Object<Poly2>(f));
      // fs.push_back(new Object<Poly2>(g));

      // debugDraw(list, fs);

      if (parity(intersections<N>(b[it].v))) {
        odd_count++;
        odd_index = it;
      } else {
        even_count++;
      }
    }

  } catch (SignException e) {
    rb.I = I;
    return false;
  }

  assert(odd_count >= 1);
  // if(odd_count < 1)
  //  return false;

  if (odd_count >= 1) {
    // I = b[odd_index];
    rb = b[odd_index];
    return true;
  } else {
    return false;
  }
}

template <class N>
void Root2::newton(RootBoundary<N>& rb) {
  PV2<N> I = rb.I;

  Poly2<N> ff = f->get<N>();
  Poly2<N> gg = g->get<N>();
  Poly2<N> fx = ff.derX();
  Poly2<N> fy = ff.derY();
  Poly2<N> gx = gg.derX();
  Poly2<N> gy = gg.derY();

  int count = 0;

  do {
    PV2<N> y(I.x.midP(), I.y.midP());

    PV2<N> foy(ff.value(&y.x), gg.value(&y.x));

    PV2<N> fr(fx.value(&I.x), fy.value(&I.x));
    PV2<N> gr(gx.value(&I.x), gy.value(&I.x));

    PV2<N> x;

    // let this be an exception in the end code
    if (!solve<N>(fr, gr, foy, x)) {
      // printf("NEWTON determinant 0, need to subdivide\n\n");
      break;
    }

    PV2<N> newI = y - x;

    // printf("Newton Result [%20.16g, %20.16g] [%20.16g, %20.16g]\n\n",
    // newI.x.lb(), newI.x.ub(), newI.y.lb(), newI.y.ub());

    // this should just be an assertion, not exception in the end code
    assert(newI.x.intersects(I.x) && newI.y.intersects(I.y));

    if (!newI.x.intersects(I.x) || !newI.y.intersects(I.y)) {
      // printf("NEWTON result disjoint\n\n");
      break;
    }

    /*
    vector<PV2> b;
    b.push_back(I);
    b.push_back(newI);
    vector< PTR<Object<Poly2>> > fs;
    fs.push_back(new Object<Poly2>(f));
    fs.push_back(new Object<Poly2n>(g));
    */
    // debugDraw(b, fs);

    PV2<N> inter(newI.x.intersect(I.x), newI.y.intersect(I.y));

    // this is a termination condition
    // if the new interval, (result of intersecting the result and old interval,
    // is not a proper subset of the old interval if(!inter.x.subset(I.x) &&
    // !inter.y.subset(I.y))

    // this is less efficient
    // if(inter.x.intervalWidth() >= I.x.intervalWidth() &&
    // inter.y.intervalWidth() >= I.y.intervalWidth()) {
    if (!inter.x.subset(I.x) && !inter.y.subset(I.y)) {
      // printf("NEWTON no progress\n\n");
      break;
    }
    /*
       int pI = parity(intersections(I));
       int pIn = parity(intersections(inter));

       if(pI == 0 || pIn == 0) {
         printf("here\n");
       }
   */
    I = inter;

    rb = {.I = I, .v = vector<vector<PTR<Object<Scalar>>>>()};

    // printf("      I       [%20.16g, %20.16g] [%20.16g, %20.16g]\n\n",
    // I.x.lb(), I.x.ub(), I.y.lb(), I.y.ub());

    // printf("NEWTON STEP\n\n");

  } while (true);
  //} while(++count < 100);

  // return I;
}

template <class N>
bool Root2::solve(PV2<N> fr, PV2<N> gr, PV2<N> b, PV2<N>& sol) {
  // pivot on the smallest element!
  // compare two things: sing(false) for comparison

  // if det contains zero, need to subdivide
  N d = fr.x * gr.y - fr.y * gr.x;
  if (!d.sign(false)) return false;

  if (!fr.x.sign(false) || !gr.x.sign(false)) return false;

  // pivot first row
  if ((fr.x.abs() - gr.x.abs()).sign(false) >= 0) {
    N c = gr.x / fr.x;

    gr.x = gr.x - c * fr.x;
    gr.y = gr.y - c * fr.y;
    b.y = b.y - c * b.x;

    if (!gr.y.sign(false)) return false;
    N y = b.y / gr.y;

    if (!fr.x.sign(false)) return false;
    N x = (b.x - fr.y * y) / fr.x;

    sol = PV2<N>(x, y);

  }
  // pivot on the second row
  else {
    N c = fr.x / gr.x;

    fr.x = fr.x - c * gr.x;
    fr.y = fr.y - c * gr.y;
    b.x = b.x - c * b.y;

    if (!fr.y.sign(false)) return false;
    N y = b.x / fr.y;

    if (!gr.x.sign(false)) return false;
    N x = (b.y - gr.y * y) / gr.x;

    sol = PV2<N>(x, y);
  }

  return true;
}

template PV2<Parameter> Root2::calculate<acp::Parameter>();
template PV2<MParameter> Root2::calculate<acp::MParameter>();
template PV2<PParameter> Root2::calculate<acp::PParameter>();
template PV2<Parameter> Root2::polish<acp::Parameter>(PV2<Parameter>);
template PV2<MParameter> Root2::polish<acp::MParameter>(PV2<Parameter>);
template PV2<PParameter> Root2::polish<acp::PParameter>(PV2<Parameter>);

}  // namespace acp
