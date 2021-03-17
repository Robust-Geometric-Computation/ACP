#ifndef CRITICAL_POINTS_H
#define CRITICAL_POINTS_H

#include "acp/core/object.h"
#include "acp/encasement2d/encasement2d.h"
#include "acp/encasement2d/predicates2d.h"
#include "acp/poly/root2.h"
#include "acp/resultant/resultant_arrangement.h"

#include <map>
#include <queue>
#include <set>
#include <vector>

using namespace acp;
using namespace std;

class Rectangle {
 public:
  PV2<Parameter> box;
  int sign_f;
  int sign_fx;
  int sign_fy;
  int sign_hess;
  int num_splits;

  static vector<Poly2<Parameter>> test_poly;

  Rectangle* parent;
  int rank;

  int id;

  Rectangle(double xl, double yl, double xu, double yu)
      : parent(this), rank(0), num_splits(0) {
    box = PV2<Parameter>(Parameter::interval(xl, xu),
                         Parameter::interval(yl, yu));
    rects.push_back(this);
    id = global_id++;
  }

  Rectangle(PV2<Parameter> box)
      : box(box), parent(this), rank(0), num_splits(0) {
    rects.push_back(this);
    id = global_id++;
  }

  Rectangle(const Rectangle& r)
      : box(r.box), parent(this), rank(0), num_splits(0) {
    rects.push_back(this);
    id = global_id++;
  }

  Rectangle(const Rectangle& r, int num_splits)
      : box(r.box), parent(this), rank(0), num_splits(num_splits) {
    rects.push_back(this);
    id = global_id++;
  }

  Rectangle& operator=(const Rectangle& r) {
    box = r.box;
    parent = this;
    rank = 0;
    id = r.id;
    return *this;
  }

  double lbx() { return box.x.lb(); }
  double ubx() { return box.x.ub(); }
  double lby() { return box.y.lb(); }
  double uby() { return box.y.ub(); }

  double midX() { return box.x.mid(); }
  double midY() { return box.y.mid(); }

  Rectangle* find();
  void unionize(Rectangle* that);

  void set_sign(vector<PTR<Object<Poly2>>> v);
  Rectangle* split();

  /*
    static void clean() {
      for(int i = 0; i < rects.size(); i++) {
        if(rects[i]) delete rects[i];
      }
      rects.clear();
      global_id = 0;
    }
  */

  static vector<Rectangle*> rects;
  static int global_id;
};

class CompareByLBY {
 public:
  bool operator()(Rectangle* const r, Rectangle* const s) const {
    return r->box.y.lb() < s->box.y.lb();
  }
};

class CompareByLBX {
 public:
  bool operator()(Rectangle* const r, Rectangle* const s) const {
    return r->box.x.lb() < s->box.x.lb();
  }
};

vector<Rectangle*> critical_points(PTR<Object<Poly2>> f, Rectangle* start);
map<int, set<Rectangle*>> connected_components(vector<Rectangle*> v);

int quad_approx_sign(PTR<Object<PV2>> xy, PTR<Object<Poly2>> f);
int cubic_approx_sign(PTR<Object<PV2>> xy, PTR<Object<Poly2>> f);
int n_approx_sign(PTR<Object<PV2>> xy, PTR<Object<Poly2>> f, int n);

#endif
