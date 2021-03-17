#include "acp/dcel/dcel.h"
#include "acp/encasement2d/plot.h"
#include "acp/poly/root.h"

bool test = false;

void Edge::split(PTR<Object<Scalar>> t_prime) {
  PTR<Edge> e_prev = prev;
  PTR<Edge> t_prev = twin->prev;
  PTR<Edge> e_next = next;
  PTR<Edge> t_next = twin->next;

  // disconnect the edge
  swap(t_prev);

  // create the new edge
  //(1) create with line instead
  // PTR<Edge>  new_edge = new Edge(line, t_prime, t2);
  // PTR<Edge>  new_edge = new Edge(new Vertex(new LinePoint(line, t_prime)),
  // head);

  PTR<Edge> new_edge = new Edge();
  new_edge->line = line;
  new_edge->t1 = t_prime;
  new_edge->t2 = t2;

  new_edge->tail = new Vertex(new LinePoint(line, t_prime));
  new_edge->head = head;

  new_edge->init();

  new_edge->type = EDGETYPE::ELINE;

  // shrink this edge
  //(2) Set the correct t1 and t2 respectively
  head = new_edge->tail;
  twin->tail = new_edge->tail;
  t2 = t_prime;
  twin->t1 = t_prime;

  // connect this to new edge
  swap(new_edge->twin);

  // connect new edge to prior connections
  new_edge->swap(t_prev);

  // set the faces
  new_edge->face = face;
  new_edge->twin->face = twin->face;

  // if there was no previous intersection, stop
  if (!intersection) return;

  PTR<Object<Scalar>> min = (LessThan(t1, t2) == -1 ? t1 : t2);
  PTR<Object<Scalar>> max = (LessThan(t1, t2) == -1 ? t2 : t1);

  /// printf/gc("\n  min: %g\n", min->get().mid());
  /// printf/gc("  min: %g\n", max->get().mid());
  /// printf/gc("  int: %g\n", intersection->get().mid());

  // The intersection is on the new edge
  if (!(LessThan(min, intersection) == -1 &&
        LessThan(intersection, max) == -1)) {
    /// printf/gc("  int goes to NEW edge\n");

    new_edge->intersection = intersection;
    new_edge->twin->intersection = intersection;
    new_edge->gradient_cross_point = new_edge->twin->gradient_cross_point =
        this->gradient_cross_point;

    this->gradient_cross_point = this->twin->gradient_cross_point = false;

    new_edge->poly = poly;
    new_edge->twin->poly = poly;

    this->poly = nullptr;
    this->twin->poly = nullptr;

    intersection = nullptr;
    twin->intersection = nullptr;
  } else {
    /// printf/gc("  int goes to edge\n\n");
  }
}

bool Edge::intersects(PTR<Line> l) {
  // return (LeftOf(l->p, l->v, tail->p)*LeftOf(l->p, l->v, head->p)) == -1;
  return l->left(tail->p) * l->left(head->p) == -1;
}

// intersect and split the edge between intersections
bool Edge::intersect(PTR<Object<Poly2>> f) {
  // if we've already intersected with this curve
  if (this->poly == f) return false;

  // Always ensure we're splitting with t1 < t2
  if (LessThan(this->t1, this->t2) == 1) {
    return twin->intersect(f);
  }

  PTR<Object<Poly>> lp = new LinePoly(this->line, f);

  vector<PTR<Object<Scalar>>> rootsR =
      PolySolver(lp).getRoots(this->t1, this->t2);

  if (rootsR.size() == 0) {
    // delete lp;
    return false;
  }

  vector<PTR<Object<Scalar>>> roots(rootsR.size());
  for (int i = 0; i < rootsR.size(); i++) roots[i] = rootsR[i];

  /// printf/gc("  roots size(): %ld\n", roots.size());

  Ascending asc;

  std::sort(roots.begin(), roots.end(), asc);

  // TEST TEST TEST
  //  if(test) {
  //
  //    PTR<Edge> joe = new Edge(line, t1, roots[0]);
  //    for(int i = 1; i < roots.size(); i++) {
  //      PTR<Edge>  n = new Edge(line, roots[i-1], roots[i]);
  //      joe->swap(n->twin);
  //      joe = n;
  //    }
  //    PTR<Edge> n = new Edge(line, roots[roots.size()-1], t2);
  //    joe->swap(n->twin);
  //    new Face(n->twin);
  //
  //    return false;
  //  }

  PTR<Edge> e = this;

  // remove the old intersection, we'll insert it later
  PTR<Object<Scalar>> old_intersection = this->intersection;
  PTR<Object<Poly2>> old_poly = this->poly;
  this->poly = nullptr;
  this->intersection = nullptr;
  twin->poly = nullptr;
  twin->intersection = nullptr;
  bool gcp = this->gradient_cross_point;
  this->gradient_cross_point = this->twin->gradient_cross_point = false;

  // split up this edge by the roots
  // If the root lies on this edge, then add the intersection to the edge
  // If the edge already has an intersection, then split at their midpoint
  // We've sorted the roots so we know that the intersection t is < roots[i]
  // So after we split, set roots[i] as the intersection point of e->next
  for (int i = 0; i < roots.size(); i++) {
    /////printf/gc("root[%d]: %g\n", i, roots[i]->get().mid());
    PTR<Object<Scalar>> root = roots[i];
    if (LessThan(e->t1, root) == -1 && LessThan(root, e->t2) == -1) {
      if (!e->intersection) {
        e->poly = f;
        e->twin->poly = f;
        e->intersection = root;
        e->twin->intersection = root;
        // This is where the face learns all of it's intersections from
        if (e->face) e->face->curves.insert(f);
        if (e->twin->face) e->twin->face->curves.insert(f);
      } else {
        // e->split(new AverageAB(0.5, e->intersection, roots[i]));
        assert(LessThan(e->intersection, root) == -1);
        assert(LessThan(root, e->t2) == -1);
        e->split(approxAverage(0.5, e->intersection, root));
        e = e->next;
        e->intersection = root;
        e->twin->intersection = root;
        e->poly = f;
        e->twin->poly = f;
      }

    } else {
      /*delete roots[i]*/;
    }
  }

  PTR<Edge> last = e->next;
  e = this;

  // insert the old intersection into the split up edges
  // We might have split up the edge into many edges when splitting
  // between roots on the polynomial. Find the edge that contains the t
  // value of the old intersection and either add it there, or split between
  // the existing intersection point (which will have to be a root above)
  // We don't know whether old_intersection < e->intersection though
  // so check before adding the intersection point
  if (old_intersection) {
    do {
      if (LessThan(e->t1, old_intersection) == -1 &&
          LessThan(old_intersection, e->t2) == -1) {
        if (!e->intersection) {
          e->poly = old_poly;
          e->twin->poly = old_poly;
          e->intersection = old_intersection;
          e->twin->intersection = old_intersection;
          e->gradient_cross_point = e->twin->gradient_cross_point = gcp;
        } else {
          // e->split(new AverageAB(0.5, e->intersection, old_intersection));
          e->split(approxAverage(0.5, e->intersection, old_intersection));
          if (e->intersection) {
            e = e->next;
          }
          e->intersection = old_intersection;
          e->twin->intersection = old_intersection;
          e->gradient_cross_point = e->twin->gradient_cross_point = gcp;
          e->poly = old_poly;
          e->twin->poly = old_poly;
          break;
        }
      }
      e = e->next;
    } while (e != last);
  }

  return true;
}

PTR<Edge> Edge::clone() {
  PTR<Edge> e = new Edge();
  e->tail = tail;
  e->head = head;

  e->line = line;
  e->t1 = t1;
  e->t2 = t2;

  PTR<Object<Poly2>> p = poly;
  PTR<Object<Scalar>> i = intersection;

  e->type = type;

  e->init();
  e->poly = e->twin->poly = p;
  e->intersection = e->twin->intersection = i;
  e->gradient_cross_point = e->twin->gradient_cross_point =
      this->gradient_cross_point;

  return e;
}

bool Face::intersect(PTR<Object<Poly2>> f) {
  vector<PTR<Edge>> edges;
  PTR<Edge> s = e;
  do {
    edges.push_back(s);
  } while ((s = s->next) != e);

  bool r = false;

  for (int i = 0; i < edges.size(); i++) {
    r = r || edges[i]->intersect(f);
  }

  return r;
}

//-------------------------------------------------------------------

PTR<Face> Face::make(PV2<Parameter>& p) {
  PV2<Parameter> ll(p.x.lbP(), p.y.lbP());
  PV2<Parameter> lr(p.x.ubP(), p.y.lbP());
  PV2<Parameter> ur(p.x.ubP(), p.y.ubP());
  PV2<Parameter> ul(p.x.lbP(), p.y.ubP());

  PTR<Vertex> ll_v = new Vertex(new InputPoint(ll));
  PTR<Vertex> lr_v = new Vertex(new InputPoint(lr));
  PTR<Vertex> ur_v = new Vertex(new InputPoint(ur));
  PTR<Vertex> ul_v = new Vertex(new InputPoint(ul));

  vector<PTR<Edge>> edges;

  edges.push_back(new Edge(ll_v, lr_v));
  edges.push_back(new Edge(lr_v, ur_v));
  edges.push_back(new Edge(ur_v, ul_v));
  edges.push_back(new Edge(ul_v, ll_v));

  for (int i = 0; i < 4; i++) {
    edges[i]->next = edges[(i + 1) % 4];
    edges[(i + 1) % 4]->prev = edges[i];
    edges[(i + 1) % 4]->twin->next = edges[i]->twin;
    edges[i]->twin->prev = edges[(i + 1) % 4]->twin;
  }

  return new Face(edges[0]);
}

//-------------------------------------------------------------------

vector<PTR<Edge>> Face::intersection_edges(PTR<Line> line) {
  vector<PTR<Edge>> i_edges;

  PTR<Edge> edge = this->e;

  do {
    // if you try to split the face withe a line that is already part of the
    // face, don't split
    if (edge->line == line) return vector<PTR<Edge>>();

    // don't check for splits with axis aligned lines that are the same
    // (identity)
    if ((!edge->axis_aligned || edge->line->axis != line->axis) &&
        edge->intersects(line))
      i_edges.push_back(edge);

    edge = edge->next;
  } while (edge != this->e);

  return i_edges;
}

PTR<Face> Face::split(PTR<Line> line) {
  vector<PTR<Edge>> i_edges = intersection_edges(line);

  if (i_edges.size() == 0) return nullptr;
  if (i_edges.size() != 2) {
    printf("FACE IS NOT CONVEX\n");
    return nullptr;
  }

  PTR<Edge> edge = this->e;

  PTR<Edge> new_edge = new Edge();

  i_edges[0]->split(new LLHit(i_edges[0]->line, line));
  i_edges[1]->split(new LLHit(i_edges[1]->line, line));

  new_edge->tail = i_edges[0]->head;
  new_edge->head = i_edges[1]->head;

  new_edge->line = line;
  new_edge->t1 = new LLHit(line, i_edges[0]->line);
  new_edge->t2 = new LLHit(line, i_edges[1]->line);
  new_edge->init();  // creates the twin etc...

  // set up connections
  new_edge->next = i_edges[1]->next;
  new_edge->prev = i_edges[0];
  new_edge->twin->next = i_edges[0]->next;
  new_edge->twin->prev = i_edges[1];

  i_edges[0]->next = new_edge;
  i_edges[1]->next = new_edge->twin;
  new_edge->next->prev = new_edge;
  new_edge->twin->next->prev = new_edge->twin;

  //(1) Create new Face. set this->e so that this face is left of line
  PTR<Face> f = new Face();

  if (LessThan(new_edge->t1, new_edge->t2) == -1) {
    this->e = new_edge;
    f->e = new_edge->twin;
    new_edge->face = this;
    new_edge->twin->face = f;
  } else {
    this->e = new_edge->twin;
    f->e = new_edge;
    new_edge->face = f;
    new_edge->twin->face = this;
  }
  /*
    ///printf/gc("tail: (%g, %g)\n", new_edge->tail->p->get().x.mid(),
                               new_edge->tail->p->get().y.mid());

    ///printf/gc("l(t1 = %g): (%g, %g)\n", new_edge->t1->get().mid(),
    new_edge->line->getP(new_edge->t1).x.mid(),
                                new_edge->line->getP(new_edge->t1).y.mid());

    ///printf/gc("head: (%g, %g)\n", new_edge->head->p->get().x.mid(),
                               new_edge->head->p->get().y.mid());

    ///printf/gc("l(t2 = %g): (%g, %g)\n", new_edge->t2->get().mid(),
    new_edge->line->getP(new_edge->t2).x.mid(),
                                new_edge->line->getP(new_edge->t2).y.mid());
  */

  // TEST TEST TEST
  test = false;

  //(2) Intersect new_edge with the polys in this face
  // Intersecting the new edge with a curve may split the
  // edge into N other edges. When checking the next curve
  // we must iterate through each of those individual edges
  // and call intersect on them, hence the start-end edges
  PTR<Edge> start =
      new_edge->prev;  // new_edge will always be start of the chain
                       // because splitting always keeps "this"'s tail
  PTR<Edge> end = new_edge->next;
  for (set<PTR<Object<Poly2>>>::iterator it = curves.begin();
       it != curves.end(); it++) {
    /// printf/gc("--------------------------\n");
    PTR<Edge> current = start->next;
    while (current != end) {
      PTR<Edge> c_prev = current->prev;  // To fix bug with missing intersection
      current->intersect(*it);
      current = c_prev->next->next;
    }
    /// printf/gc("--------------------------\n");
  }

  // TEST TEST TEST
  test = false;
  // return f;

  //(3) Establish all of the polys still in the face
  // Loop through each edge and collect the curves
  // that they intersect in the set
  this->curves.clear();
  PTR<Edge> current = this->e;
  do {
    if (current->poly) this->curves.insert(current->poly);
    current->face = this;
    current = current->next;
  } while (current != this->e);

  //(4) Establish all of the polys in the new face
  f->curves.clear();
  current = f->e;
  do {
    if (current->poly) f->curves.insert(current->poly);
    current->face = f;
    current = current->next;
  } while (current != f->e);

  // Always return the face to the right of the line
  return f;
}

bool Face::inside(PTR<Object<PV2>> p) {
  PTR<Edge> current = e;

  do {
    if (Orient2D(current->tail->p, current->head->p, p) >= 0) return false;

    current = current->next;

  } while (current != e);

  return true;
}

bool Face::contains(PTR<Object<Poly2>> f) {
  // return intersectionCount(f) != 0;
  // return curves.find(f) != curves.end();
  PTR<Edge> edge = e;
  do {
    if (edge->poly == f) return true;
  } while ((edge = edge->next) != e);

  return false;
}

// Rather expensive operation, tries to intersect
// each edge with the poly and checks for roots
// in the interval
bool Face::intersects(PTR<Object<Poly2>> f) {
  PTR<Edge> current = e;
  do {
    LinePoly* lp = new LinePoly(e->line, f);

    bool less = LessThan(e->t1, e->t2) == -1;

    vector<PTR<Object<Scalar>>> roots =
        PolySolver(lp).getRoots(less ? e->t1 : e->t2, less ? e->t2 : e->t1);
    bool inter = roots.size() > 0;
    for (int i = 0; i < roots.size(); i++) {
      /*delete roots[i]*/;
    }
    if (inter) return true;

  } while ((current = current->next) != e);

  return false;
}

int Face::intersectionCount(PTR<Object<Poly2>> f) {
  int count = 0;

  PTR<Edge> e = this->e;

  do {
    if (e->poly == f) count++;
    e = e->next;
  } while (e != this->e);

  return count;
}

// Test whether f and g intersect this face in an alternating manner
// If we find the first edge that f intersects, we don't know whether
// this is the entrance or exit of this branch of f, going counterclockwise

// Precondition: f and g both enter once and exit once
/*
bool Face::alternates(PTR<Object<Poly2>>  f, PTR<Object<Poly2>>  g) {
  assert(contains(f));
  assert(contains(g));

  assert(f != g);

  int fc = 0;
  int gc = 0;

  PTR<Edge> j = e;
  do {
    if(j->poly == f) fc++;
    if(j->poly == g) gc++;
    j = j->next;
  } while(j != e);

  //printf("\nfc: %d\n", fc);
  //printf("\ngc: %d\n", gc);

  PTR<Edge> start = e;

  while(start->poly != f) start = start->next;

  int count = 0;

  PTR<Edge> current = start;
  do {
    if(current->poly == g) count++;
  } while((current = current->next) != start && current->poly != f);

  return count % 2 == 1;
}
*/

bool Face::alternates(PTR<Object<Poly2>> f, PTR<Object<Poly2>> g) {
  assert(contains(f));
  assert(contains(g));

  assert(f != g);

  vector<int> alt;

  PTR<Edge> edge = e;
  do {
    if (edge->poly == f) alt.push_back(0);
    if (edge->poly == g) alt.push_back(1);
    edge = edge->next;
  } while (edge != e);

  int k = 0;
  int j = 0;
  while (j < alt.size()) {
    if (alt[j] == 0) break;
    j++;
  }

  // if f doesn't intersect this cell, keep going
  if (j == alt.size()) {
    return false;
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

  return parity % 2 == 1;
}

// Perform Newton's to try to closest approach
PTR<Object<PV2>> Face::closest_approach(PTR<Object<Poly2>> F,
                                        PTR<Object<Poly2>> G,
                                        PTR<Object<PV2>> start) {
  // turn on normal rounding mode
  assert(enabled);
  disable();

  PPoly2D f(F->getApprox(1.0));
  PPoly2D g(G->getApprox(1.0));

  // Our G function is grad(f) X grad(g) = 0

  PV2D delta_p(0, 0);

  //(1) Calculate the center of the cell
  PV2D p = centroid();
  // PV2D p(start->get().x.mid(), start->get().y.mid());
  PV2D cent = p;

  /////printf/gc("\nface: 0x%lx center (%f, %f)\n", (long)this, p.x, p.y);

  //(2) Solve for delta P and iterate

  double xy[2];

  PPoly2D fx = f.derX();
  PPoly2D fy = f.derY();
  PPoly2D gx = g.derX();
  PPoly2D gy = g.derY();

  PPoly2D fxx = fx.derX();
  PPoly2D fxy = fx.derY();
  PPoly2D fyy = fy.derY();

  PPoly2D gxx = gx.derX();
  PPoly2D gxy = gx.derY();
  PPoly2D gyy = gy.derY();

  xy[0] = p.x;
  xy[1] = p.y;

  double f_p = f.value(xy);
  double g_p = (fx.value(xy) * gy.value(xy) - fy.value(xy) * gx.value(xy));

  double maxf = fmax(fabs(f_p), fabs(g_p));

  int max_count = 0;
  int MAX_COUNT = 5;

  int iter = 0;
  int min_iter = 4;

  do {
    double fx_p = fx.value(xy);
    double fy_p = fy.value(xy);

    double fxv = fx_p;
    double fyv = fy_p;
    double gxv = gx.value(xy);
    double gyv = gy.value(xy);

    double fxxv = fxx.value(xy);
    double fxyv = fxy.value(xy);
    double fyyv = fyy.value(xy);

    double gxxv = gxx.value(xy);
    double gxyv = gxy.value(xy);
    double gyyv = gyy.value(xy);

    double gx_p = (fxxv * gyv + fxv * gxyv - fxyv * gxv - fyv * gxxv);
    double gy_p = (fxyv * gyv + fxv * gyyv - fyyv * gxv - fyv * gxyv);

    DrawItem item;
    vector<PTR<Object<PV2>>> draw_points;
    vector<PTR<Object<PV2>>> all_points;

    item.types.push_back(POINT);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.push_back(new InputPoint(PV2<Parameter>::constant(p.x, p.y)));
    all_points.push_back(draw_points[draw_points.size() - 1]);
    item.pvs.push_back(draw_points);

    PV2<Parameter> gradf = PV2<Parameter>::constant(fxv, fyv);
    PV2<Parameter> gradg = PV2<Parameter>::constant(gxv, gyv);
    gradf = gradf.unit();
    gradg = gradg.unit();
    PV2<Parameter> pp = PV2<Parameter>::constant(p.x, p.y);

    PV2<Parameter> pfg = pp + 0.3 * gradf;
    pfg = PV2<Parameter>::constant(pfg.x.mid(), pfg.y.mid());

    PV2<Parameter> pgg = pp + 0.3 * gradg;
    pgg = PV2<Parameter>::constant(pgg.x.mid(), pgg.y.mid());

    // item.types.push_back(VECTOR);
    // item.colors.push_back({0.7f, 0.0f, 0.0f});
    // draw_points.clear();
    // draw_points.push_back(prev);
    // draw_points.push_back(curr);
    // item.pvs.push_back(draw_points);

    item.types.push_back(POINT);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(new InputPoint(pp));
    item.pvs.push_back(draw_points);

    item.types.push_back(VECTOR);
    item.colors.push_back({0.7f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(new InputPoint(pp));
    draw_points.push_back(new InputPoint(pfg));
    item.pvs.push_back(draw_points);

    item.types.push_back(VECTOR);
    item.colors.push_back({0.0f, 0.7f, 0.0f});
    draw_points.clear();
    draw_points.push_back(new InputPoint(pp));
    draw_points.push_back(new InputPoint(pgg));
    item.pvs.push_back(draw_points);

    //    Parameter::enable();
    //
    //    double zoom = debugDrawFacePoint(this, F, G, 1, cent.x, cent.y, item,
    //    true); while(zoom < 0) {
    //      zoom = -zoom;
    //      zoom = debugDrawFacePoint(this, F, G, zoom, p.x, p.y, item, true);
    //    }
    //
    //    Parameter::disable();

    Mat2 M(fx_p, fy_p, gx_p, gy_p);
    PV2D b(-f_p, -g_p);

    delta_p = M.inverse() * b;

    p = p + delta_p;

    xy[0] = p.x;
    xy[1] = p.y;

    f_p = f.value(xy);

    fxv = fx.value(xy);
    fyv = fy.value(xy);
    gxv = gx.value(xy);
    gyv = gy.value(xy);
    g_p = (fxv * gyv - fyv * gxv);

    double lmaxf = fmax(fabs(f_p), fabs(g_p));
    if (lmaxf < maxf) {
      maxf = lmaxf;
      max_count = 0;
    } else if (iter >= min_iter) {
      max_count++;
    }

    if (!this->inside(new InputPoint(PV2<Parameter>::constant(p.x, p.y)))) {
      enable();
      return nullptr;
    }

    // PTR<Object<PV2>> curr = new InputPoint(PV2<Parameter>::constant(p.x,
    // p.y)); PTR<Object<PV2>> prev = all_points[all_points.size()-1];
    // all_points.push_back(curr);

    // PV2<Parameter> gradf = PV2<Parameter>::constant(fxv, fyv);
    // PV2<Parameter> gradg = PV2<Parameter>::constant(gxv, gyv);
    // gradf = gradf.unit();
    // gradg = gradg.unit();
    // PV2<Parameter> pp = PV2<Parameter>::constant(p.x, p.y);

    // PV2<Parameter> pfg = pp + 0.3*gradf;
    // pfg = PV2<Parameter>::constant(pfg.x.mid(), pfg.y.mid());

    // PV2<Parameter> pgg = pp + 0.3*gradg;
    // pgg = PV2<Parameter>::constant(pgg.x.mid(), pgg.y.mid());

    ////item.types.push_back(VECTOR);
    ////item.colors.push_back({0.7f, 0.0f, 0.0f});
    ////draw_points.clear();
    ////draw_points.push_back(prev);
    ////draw_points.push_back(curr);
    ////item.pvs.push_back(draw_points);
    //
    // item.types.push_back(POINT);
    // item.colors.push_back({0.0f, 0.0f, 0.0f});
    // draw_points.clear();
    // draw_points.push_back(new InputPoint(pp));
    // item.pvs.push_back(draw_points);

    // item.types.push_back(VECTOR);
    // item.colors.push_back({0.7f, 0.0f, 0.0f});
    // draw_points.clear();
    // draw_points.push_back(new InputPoint(pp));
    // draw_points.push_back(new InputPoint(pfg));
    // item.pvs.push_back(draw_points);

    // item.types.push_back(VECTOR);
    // item.colors.push_back({0.0f, 0.7f, 0.0f});
    // draw_points.clear();
    // draw_points.push_back(new InputPoint(pp));
    // draw_points.push_back(new InputPoint(pgg));
    // item.pvs.push_back(draw_points);

    // Parameter::enable();

    // zoom = debugDrawFacePoint(this, F, G, 1, p.x, p.y, item, true);
    // while(zoom < 0) {
    //  zoom = -zoom;
    //  zoom = debugDrawFacePoint(this, F, G, 1, p.x, p.y, item, true);
    //}
    //
    // Parameter::disable();

    iter++;

  } while (max_count < MAX_COUNT);

  enable();

  PTR<Object<PV2>> point = new InputPoint(PV2<Parameter>::constant(p.x, p.y));

  if (!this->inside(point)) {
    point = nullptr;
  }

  return point;
}
// Perform Newton's to try to find intersection of these curves inside
// this face
PTR<Object<PV2>> Face::intersection(PTR<Object<Poly2>> F,
                                    PTR<Object<Poly2>> G) {
  // turn on normal rounding mode
  assert(enabled);
  disable();

  PPoly2D f(F->getApprox(1.0));
  PPoly2D g(G->getApprox(1.0));

  double EPS = 1e-6;

  PV2D delta_p(0, 0);

  //(1) Calculate the center of the cell
  PV2D p = centroid();
  PV2D cent = p;

  DrawItem item;
  vector<PTR<Object<PV2>>> draw_points;
  vector<PTR<Object<PV2>>> all_points;

  item.types.push_back(POINT);
  item.colors.push_back({0.0f, 0.0f, 0.0f});
  draw_points.push_back(new InputPoint(PV2<Parameter>::constant(p.x, p.y)));
  all_points.push_back(draw_points[draw_points.size() - 1]);
  item.pvs.push_back(draw_points);

  enable();

  // double zoom = debugDrawFacePoint(this, F, G, 1, cent.x, cent.y, item,
  // true); while(zoom < 0) {
  //  zoom = -zoom;
  //  zoom = debugDrawFacePoint(this, F, G, zoom, p.x, p.y, item, true);
  //}

  disable();

  /////printf/gc("\nface: 0x%lx center (%f, %f)\n", (long)this, p.x, p.y);

  //(2) Solve for delta P and iterate

  double xy[2];

  PPoly2D fx = f.derX();
  PPoly2D fy = f.derY();
  PPoly2D gx = g.derX();
  PPoly2D gy = g.derY();

  xy[0] = p.x;
  xy[1] = p.y;

  double f_p = f.value(xy);
  double g_p = g.value(xy);

  double maxf = fmax(fabs(f_p), fabs(g_p));

  int max_count = 0;
  int MAX_COUNT = 5;

  int iter = 0;
  int min_iter = 4;

  do {
    double fx_p = fx.value(xy);
    double fy_p = fy.value(xy);
    double gx_p = gx.value(xy);
    double gy_p = gy.value(xy);

    // calculate the sign here for the new f/g minima finding
    Mat2 M(fx_p, fy_p, gx_p, gy_p);
    PV2D b(-f_p, -g_p);

    delta_p = M.inverse() * b;

    p = p + delta_p;

    xy[0] = p.x;
    xy[1] = p.y;

    f_p = f.value(xy);
    g_p = g.value(xy);

    double lmaxf = fmax(fabs(f_p), fabs(g_p));
    if (lmaxf < maxf) {
      maxf = lmaxf;
      max_count = 0;
    } else if (iter >= min_iter) {
      max_count++;
    }

    PTR<Object<PV2>> curr = new InputPoint(PV2<Parameter>::constant(p.x, p.y));
    PTR<Object<PV2>> prev = all_points[all_points.size() - 1];
    all_points.push_back(curr);

    item.types.push_back(VECTOR);
    item.colors.push_back({0.7f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(prev);
    draw_points.push_back(curr);
    item.pvs.push_back(draw_points);

    item.types.push_back(POINT);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(curr);
    item.pvs.push_back(draw_points);

    enable();

    // zoom = debugDrawFacePoint(this, F, G, fmax(2*lmaxf, 1e-15), p.x, p.y,
    // item, true); while(zoom < 0) {
    //  zoom = -zoom;
    //  zoom = debugDrawFacePoint(this, F, G, zoom, p.x, p.y, item, true);
    //}
    //
    // Parameter::disable();

    iter++;

  } while (max_count < MAX_COUNT);

  enable();

  PTR<Object<PV2>> point = new InputPoint(PV2<Parameter>::constant(p.x, p.y));

  if (!this->inside(point)) {
    point = nullptr;
  }

  return point;
}

bool Face::verify(PTR<Object<Poly2>> f, PTR<Object<Poly2>> g,
                  PTR<Object<PV2>> p) {
  if (!inside(p)) return false;

  if (!alternates(f, g)) return false;

  // printf("f: %d\n", static_cast<Ellipse*>((Object<Poly2>*)f)->id);
  // printf("g: %d\n", static_cast<Ellipse*>((Object<Poly2>*)g)->id);

  PTR<Edge> c = e;
  do {
    if (false && c->poly)
      printf("c->poly: %d\n",
             static_cast<Ellipse*>((Object<Poly2>*)c->poly)->id);
    if (c->poly != nullptr && c->poly != f && c->poly != g) return false;
  } while ((c = c->next) != e);

  intersection_point = p;
  this->f = f;
  this->g = g;

  return true;
}

PV2D Face::centroid() {
  PV2D p(0, 0);

  double total_weight = 0;

  PV2D o(e->tail->p->getApprox(1.0));

  PTR<Edge> current = e->next;

  do {
    PV2D a(current->tail->p->getApprox(1.0));
    PV2D b(current->head->p->getApprox(1.0));

    PV2D oa = a - o;
    PV2D ob = b - o;

    PV2D center = o + ((oa + ob) / 3.0);

    double weight = abs(oa.cross(ob));
    total_weight += weight;

    p = p + (weight * center);

  } while ((current = current->next) != e);

  return (p / total_weight);
}

PTR<Face> Face::clone() {
  PTR<Face> fc = new Face();

  fc->satisfied = satisfied;
  fc->narrowed = narrowed;
  fc->num_splits = num_splits;

  fc->curves = curves;
  fc->intersection_point = intersection_point;
  fc->f = f;
  fc->g = g;

  fc->e = e->clone();
  fc->e->face = fc;

  PTR<Edge> edge = e->next;
  ;
  PTR<Edge> n_edge = fc->e;

  do {
    PTR<Edge> nn_edge = edge->clone();
    n_edge->next = nn_edge;
    nn_edge->prev = n_edge;
    n_edge->twin->prev = nn_edge->twin;
    nn_edge->twin->next = n_edge->twin;

    nn_edge->face = fc;

    edge = edge->next;
    n_edge = nn_edge;
  } while (edge != e);

  n_edge->next = fc->e;
  fc->e->prev = n_edge;
  n_edge->twin->prev = fc->e->twin;
  fc->e->twin->next = n_edge->twin;

  return fc;
}

bool Face::rect() const {
  PTR<Edge> edge = e;

  int count = 0;
  do {
    if (!edge->axis_aligned) return false;
    edge = edge->next;
    // when we make a turn around the corner of the rect
    // if all previous edges are axis aligned, they should
    // share the same axis aligned line until the corner
    if (edge->line->axis != edge->prev->line->axis) count++;
  } while (edge != e);

  // printf("turn count: %d\n", count);

  return count == 4;
}

std::vector<PTR<Object<Scalar>>> Face::bounds() const {
  std::vector<PTR<Object<Scalar>>> list;

  PTR<Edge> edge = e;

  PTR<Object<Scalar>> xmin = new PVScalar(e->tail->p, 0);
  PTR<Object<Scalar>> xmax = new PVScalar(e->tail->p, 0);
  PTR<Object<Scalar>> ymin = new PVScalar(e->tail->p, 1);
  PTR<Object<Scalar>> ymax = new PVScalar(e->tail->p, 1);

  // TODO
  // maybe need to do this using the identity detection????

  while ((edge = edge->next) != e) {
    PTR<Object<Scalar>> x = new PVScalar(edge->tail->p, 0);
    PTR<Object<Scalar>> y = new PVScalar(edge->tail->p, 1);

    if ((x->getApprox(1.0).x - xmin->getApprox(1.0).x).sign(false) < 0)
      xmin = x;

    if ((y->getApprox(1.0).x - ymin->getApprox(1.0).x).sign(false) < 0)
      ymin = y;

    if ((xmax->getApprox(1.0).x - x->getApprox(1.0).x).sign(false) < 0)
      xmax = x;

    if ((ymax->getApprox(1.0).x - y->getApprox(1.0).x).sign(false) < 0)
      ymax = y;
  }

  list.push_back(xmin);
  list.push_back(ymin);
  list.push_back(xmax);
  list.push_back(ymax);

  return list;
}

template <class N>
PV2<N> BoundingBox::calculate() {
  N xl, xh, yl, yh;

  PTR<Edge> e = face->e;
  PV2<Parameter> p = e->tail->p->get<N>();
  xl = p.x.lbP();
  xh = p.x.ubP();
  yl = p.y.lbP();
  yh = p.y.ubP();

  while ((e = e->next) != face->e) {
    p = e->tail->p->get<N>();
    if (p.x.lbP() < xl) xl = p.x.lbP();
    if (p.x.ubP() > xh) xh = p.x.ubP();
    if (p.y.lbP() < yl) yl = p.y.lbP();
    if (p.y.ubP() > yh) yh = p.y.ubP();
  }

  return PV2<N>(xl.interval(xh), yl.interval(yh));
}
