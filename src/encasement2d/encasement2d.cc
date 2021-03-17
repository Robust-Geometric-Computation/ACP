#include "acp/encasement2d/encasement2d.h"
#include <unistd.h>
#include <chrono>
#include <stack>
#include <unordered_map>
#include "acp/encasement2d/plot.h"
#include "acp/poly/root.h"

void Encasement::calculate(vector<PTR<Object<Poly2>>> theCurves, bool verify) {
  sites.clear();

  for (int i = 0; i < theCurves.size(); i++) {
    sites.push_back(vector<PTR<Object<PV2>>>());
    // printf("%s (%s) poly: %p\n", (i ? "g" : "f"), (i ? "green" : "red"),
    // (Object<Poly2>*)theCurves[i]);
  }

  auto t1 = std::chrono::high_resolution_clock::now();

  if (verify)
    make_sites_verify(theCurves);
  else
    make_sites(theCurves);

  auto t2 = std::chrono::high_resolution_clock::now();

  // if(print)
  //   cout << "loops: " <<
  //   std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "
  //   ms\n";

  t1 = std::chrono::high_resolution_clock::now();

  bool aggregate = false;

  if (aggregate) {
    init();

    _calculate(theCurves, -1, -1);

    vector<PTR<Object<PV2>>> r = _getRoots();

    if (r.size() > 0) theRoots.insert(theRoots.end(), r.begin(), r.end());

  } else {
    for (int i = 0; i < theCurves.size(); i++) {
      for (int j = i + 1; j < theCurves.size(); j++) {
        vector<PTR<Object<Poly2>>> v;
        v.push_back(theCurves[i]);
        v.push_back(theCurves[j]);

        init();
        // for(int i = 0; i < v.size(); i++) {
        //  sites.push_back(vector< PTR<Object<PV2>> >());
        //}

        // make_sites(v);

        _calculate(v, i, j);

        vector<PTR<Object<PV2>>> r = _getRoots();

        // if(r.size() == 0)
        //  printf("curves %d and %d don't intersect\n", i, j);

        if (r.size() > 0) theRoots.insert(theRoots.end(), r.begin(), r.end());
      }
    }

    if (theCurves.size() == 1) {
      vector<PTR<Object<Poly2>>> v;
      v.push_back(theCurves[0]);
      init();
      // for(int i = 0; i < v.size(); i++) {
      //    sites.push_back(vector< PTR<Object<PV2>> >());
      //  }
      // make_sites(v);
      _calculate(v, 0, 0);
    }
  }

  t2 = std::chrono::high_resolution_clock::now();

  // if(print)
  //   cout << "calcs: " <<
  //   std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "
  //   ms\n";
}

void Encasement::_calculate(vector<PTR<Object<Poly2>>> theCurves, int i,
                            int j) {
  this->curves.insert(this->curves.end(), theCurves.begin(), theCurves.end());

  count = 1;

  // printf("STEPS: %d\n", STEPS);

  // srandom(this->seed);

  // printf("r: %g\n", randomNumber(0, 1));

  // first intersect with the bounding box, then eliminate
  // loops, intersecting the remaining unintersected curves
  // with the new boundaries
  intersect_edges(curves);  //(2)

  // printf("r: %g\n", randomNumber(0, 1));
  if (STEPS > 0 && count++ >= STEPS) {
    return;
  }

  if (i == -1 && j == -1) {
    for (int k = 0; k < curves.size(); k++) {
      eliminate_loops(k);
    }
  } else {
    eliminate_loops(i);  //(1)
    if (STEPS > 0 && count >= STEPS) return;
    if (i != j) eliminate_loops(j);
  }
  // printf("r: %g\n", randomNumber(0, 1));

  intersect_edges(curves);

  // printf("r: %g\n", randomNumber(0, 1));

  // TEST
  // refine_narrow(faces[0], *(faces[0]->curves.begin()));
  // refine_subdiv(faces[2], *(faces[2]->curves.begin()));
  // TEST

  // Two queues will be established, one for unsatisfied faces with
  // non-alternating curves, one for unsatisfied faces with alternating
  // curves. Solve each face appropriately until all faces are satisfied
  // and the queues are empty.

  // Set up queues

  stack<Item> non_inter_queue;
  stack<Item> inter_queue;
  stack<Item> self_queue;

  stack<Item> refine_narrow_queue;
  stack<Item> refine_subdiv_queue;

  for (int i = 0; i < faces.size(); i++) {
    add_to_queues(faces[i], non_inter_queue, inter_queue, self_queue);
  }

  //(4) Isolate non-intersecting segments
  //  If two curves intersect a face in a non-alternating manner
  //  creating a splitting line and split the face. Check the new
  //  face and edge for intersections. If this happens to a curve
  //  more than twice, refine the curve.

  //(5) Isolate intersections
  //  If two curves intersect a face in an alternating manner, that
  //  face contains an intersection. Use Newton's method to find the
  //  intersection point and then isolate it by bounding the gradient
  //  cones of the curves around the point. Check the resulting faces
  //  and edges for intersections. If Newton's fails, refine the face
  //  and deal with resulting faces.

  int split_count = 0;
  int isolation_count = 0;

  int parity = ((count + 1) % 2);
  // printf("parity: %d\n", parity);

  while (!non_inter_queue.empty() || !inter_queue.empty() ||
         !self_queue.empty() || !refine_narrow_queue.empty() ||
         !refine_subdiv_queue.empty()) {
    if (STEPS > 0) {
      count++;

      // if(count < STEPS && count%2 == parity) {

      //  if(count + 1 >= STEPS) {
      //    if(!self_queue.empty()) {
      //      curve_trace2(self_queue.front().face);
      //    } else if(!non_inter_queue.empty()) {
      //      curve_trace2(non_inter_queue.front().face);
      //    } else if(!inter_queue.empty()) {
      //      curve_trace2(inter_queue.front().face);
      //    }
      //    return;
      //  } else {
      //    continue;
      //  }
      //}

      if (count >= STEPS) return;

      // printf("count: %d STEPS: %d\t done\n\n", count, STEPS);
    }

    // printf("r: %g\n", randomNumber(0, 1));

    /*
    if(count % 10000 == 0) {
      sleep(1);
    }
    */

    if (!self_queue.empty()) {
      Item item = self_queue.top();
      self_queue.pop();

      PTR<Face> other = isolate_curves(item.face, item.f);

      add_to_queues(item.face, non_inter_queue, inter_queue, self_queue);
      add_to_queues(other, non_inter_queue, inter_queue, self_queue);

      continue;
    }

    else if (!refine_narrow_queue.empty()) {
      Item item = refine_narrow_queue.top();
      refine_narrow_queue.pop();

      refine_subdiv_queue.push(item);

      continue;

      // bool f_or_g = (rand() / (double)RAND_MAX) < 0.5;
      // bool f_or_g = (count % 2 == 0);
      bool f_or_g = randomNumber(0, 1) < 0.5;

      int face_count = refine_narrow(item.face, f_or_g ? item.f : item.g);

      /////printf/gc("face_count: %d\n", face_count);

      if (item.face->contains(item.f) && item.face->contains(item.g) &&
          item.face->intersectionCount(f_or_g ? item.f : item.g) == 2)
        refine_subdiv_queue.push({item.face, item.f, item.g, f_or_g});
      else
        add_to_queues(item.face, non_inter_queue, inter_queue, self_queue);

      if (face_count > 0)
        add_to_queues(faces[faces.size() - 1], non_inter_queue, inter_queue,
                      self_queue);

      if (face_count > 1)
        add_to_queues(faces[faces.size() - 2], non_inter_queue, inter_queue,
                      self_queue);

      continue;
    }

    else if (!refine_subdiv_queue.empty()) {
      Item item = refine_subdiv_queue.top();
      refine_subdiv_queue.pop();

      refine_subdiv(item.face, item.alternates ? item.f : item.g);

      add_to_queues(item.face, non_inter_queue, inter_queue, self_queue);
      add_to_queues(faces.back(), non_inter_queue, inter_queue, self_queue);

      continue;
    }

    // separate curves
    // if(inter_queue.empty() || (count % 2 == 0)) {
    else if (!non_inter_queue.empty()) {
      Item item = non_inter_queue.top();
      non_inter_queue.pop();

      // TODO testing separation
      // auto t1 = std::chrono::high_resolution_clock::now();
      // PTR<Face>  other = separate(item.face, item.f, item.g);
      PTR<Face> other = separate_new(item.face, item.f, item.g);
      // auto t2 = std::chrono::high_resolution_clock::now();
      // cout << "separate_new: " <<
      // 0.001*std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()
      // << " ms\n";

      split_count++;

      // assert(other != item.face);

      // if the face still contains the two curves, perform newton's iter
      // if it converges, isolate, if not after 2 splits, refine
      if (item.face->contains(item.f) && item.face->contains(item.g) &&
          item.face->alternates(item.f, item.g)) {
        int num_created = isolate_intersection(item.face, item.f, item.g);
        isolation_count++;
        bool to_queue = num_created < 0;

        if (to_queue) num_created = (-1 * num_created) - 1;

        // printf("num_created: %d\n", num_created);
        // printf("to_queue: %s\n", to_queue ? "yes" : "no");

        while (num_created > 0) {
          add_to_queues(faces[faces.size() - num_created], non_inter_queue,
                        inter_queue, self_queue);
          num_created--;
        }

        if (to_queue) {
          /////printf/gc("\nface split %d times\n", item.face->num_splits);
          if (item.face->num_splits < 1) {
            item.face->num_splits++;
            item.face->f = item.f;
            item.face->g = item.g;
            add_to_queues(item.face, non_inter_queue, inter_queue, self_queue);
          } else {
            item.face->num_splits = 0;
            item.face->f = nullptr;
            item.face->g = nullptr;
            if (item.face->intersectionCount(item.f) == 2 &&
                item.face->intersectionCount(item.g) == 2) {
              refine_narrow_queue.push({item.face, item.f, item.g, false});
            } else {
              add_to_queues(item.face, non_inter_queue, inter_queue,
                            self_queue);
            }
          }
        }

      } else {
        item.face->num_splits = 0;
        item.face->f = nullptr;
        item.face->g = nullptr;
        bool added =
            add_to_queues(item.face, non_inter_queue, inter_queue, self_queue);
        /////printf/gc("face %s added queues\n", added ? "" : "not");
      }

      if (!other) continue;

      if (other->contains(item.f) && other->contains(item.g) &&
          other->alternates(item.f, item.g)) {
        int num_created = isolate_intersection(other, item.f, item.g);
        isolation_count++;
        bool to_queue = num_created < 0;

        if (to_queue) num_created = (-1 * num_created) - 1;

        // printf("num_created: %d\n", num_created);
        // printf("to_queue: %s\n", to_queue ? "yes" : "no");

        while (num_created > 0) {
          add_to_queues(faces[faces.size() - num_created], non_inter_queue,
                        inter_queue, self_queue);
          num_created--;
        }

        if (to_queue) {
          if (other->num_splits < 1) {
            other->num_splits++;
            other->f = item.f;
            other->g = item.g;
            add_to_queues(other, non_inter_queue, inter_queue, self_queue);
          } else {
            other->num_splits = 0;
            other->f = nullptr;
            other->g = nullptr;
            if (other->intersectionCount(item.f) == 2 &&
                other->intersectionCount(item.g) == 2) {
              refine_narrow_queue.push({other, item.f, item.g, false});
            } else {
              add_to_queues(other, non_inter_queue, inter_queue, self_queue);
            }
          }
        }

      } else {
        other->num_splits = 0;
        other->f = nullptr;
        other->g = nullptr;
        bool added =
            add_to_queues(other, non_inter_queue, inter_queue, self_queue);
        /////printf/gc("other %s added queues\n", added ? "" : "not");
      }

      continue;

    }
    // isolate intersections
    else if (!inter_queue.empty()) {
      /////printf/gc("inter_queue\n");
      Item item = inter_queue.top();
      inter_queue.pop();

      int num_created = isolate_intersection(item.face, item.f, item.g);
      isolation_count++;
      bool to_queue = num_created < 0;

      if (to_queue) num_created = (-1 * num_created) - 1;

      // printf("num_created: %d\n", num_created);
      // printf("to_queue: %s\n", to_queue ? "yes" : "no");

      while (num_created > 0) {
        add_to_queues(faces[faces.size() - num_created], non_inter_queue,
                      inter_queue, self_queue);
        num_created--;
      }

      /////printf/gc("\n\nintersection isolation failed\n\n");
      if (to_queue) {
        refine_narrow_queue.push(item);
      }
    }
  }

  // printf("\nDONE\nencasement calculated\n\n");
  // printf("split count: %d\n", split_count);
  // printf("isolation count: %d\n", split_count);
}

#ifdef BLEEN
PTR<Encasement> Encasement::intersect(PTR<Encasement> that) {
  for (int i = 0; i < this->faces.size(); i++) {
    PTR<Face> f1 = this->faces[i];
    for (int j = 0; j < that->faces.size(); j++) {
      PTR<Face> f2 = that->faces[j];

      if (Intersects(f1->bb(), f2->bb()) == 1) {
      }
    }
  }

  // return nullptr;
}
#endif

// Encase a single curve and restrict it to cells
// where the curve is monotonic in all but one direction
//
// To begin with the directions are just fx and fy
void Encasement::calculate_single(PTR<Object<Poly2>> f) {
  init();

  PTR<Object<Poly2>> fx = new DerXPoly2(f);
  PTR<Object<Poly2>> fy = new DerYPoly2(f);

  PTR<Object<PV2>> dx =
      new NormalVector(new InputPoint(PV2<Parameter>::constant(1, 0)));
  PTR<Object<PV2>> dy =
      new NormalVector(new InputPoint(PV2<Parameter>::constant(0, 1)));
  PTR<Object<PV2>> ndy =
      new NormalVector(new InputPoint(PV2<Parameter>::constant(0, -1)));
  PTR<Object<PV2>> vxy =
      new NormalVector(new InputPoint(PV2<Parameter>::constant(1, 1)));
  PTR<Object<PV2>> vyx =
      new NormalVector(new InputPoint(PV2<Parameter>::constant(1, -1)));

  PTR<Object<Poly2>> fxfy = new DirectionalDerivativePoly2(f, vxy);
  PTR<Object<Poly2>> fyfx = new DirectionalDerivativePoly2(f, vyx);

  // PTR<Object<Poly2>> mid1 = new DirectionalDerivativePoly2(f, new
  // NormalVector(new SumPoint(ndy, vyx))); PTR<Object<Poly2>> mid2 = new
  // DirectionalDerivativePoly2(f, new NormalVector(new SumPoint(vyx, dx)));
  // PTR<Object<Poly2>> mid3 = new DirectionalDerivativePoly2(f, new
  // NormalVector(new SumPoint(dx, vxy))); PTR<Object<Poly2>> mid4 = new
  // DirectionalDerivativePoly2(f, new NormalVector(new SumPoint(vxy, dy)));

  vector<PTR<Object<Poly2>>> theCurves;
  theCurves.push_back(f);
  theCurves.push_back(fx);
  theCurves.push_back(fy);
  theCurves.push_back(fxfy);
  theCurves.push_back(fyfx);

  // theCurves.push_back(mid1);
  // theCurves.push_back(mid2);
  // theCurves.push_back(mid3);
  // theCurves.push_back(mid4);

  this->curves.insert(this->curves.end(), theCurves.begin(), theCurves.end());

  // (1) eliminate any possible loops
  sites.clear();

  for (int i = 0; i < curves.size(); i++) {
    sites.push_back(vector<PTR<Object<PV2>>>());
  }

  make_sites(theCurves);

  intersect_edges(curves);

  for (int i = 0; i < curves.size(); i++) eliminate_loops(i);

  intersect_edges(curves);

  printf("loops eliminated\n");

  // (2) Keep splitting with widest pair from set of all shortest pairs of
  //     boundary intersections, until there are <= 2 curves left in a cell
  //     e.g. {(f(p), fx(q)), (f(r), fy(s))}, if first pair of points has
  //     further distance apart (boundary distance?) use them to separate

  queue<PTR<Face>> q;

  for (int i = 0; i < faces.size(); i++) {
    if (faces[i]->contains(f) && faces[i]->curves.size() > 3) q.push(faces[i]);
  }

  while (q.size() > 0) {
    PTR<Face> face = q.front();
    q.pop();

    PTR<Face> other = separate(face);

    if (face->contains(f) && face->curves.size() > 3) q.push(face);
    if (other->contains(f) && other->curves.size() > 3) q.push(other);
  }
}

PTR<Face> Encasement::separate(PTR<Face> face) {
  // printf("Separate F and G\n");

  PTR<Object<Poly2>> f = curves[0];

  // look at all pairs of intersections along that boundary where the first
  // intersection is one with f. For each other curve, take the pair that has
  // the smallest distance. Out of all of those smallest pairs, choose the one
  // with largest distance around the boundary and split between those two
  vector<PTR<Edge>> f_edges;
  vector<PTR<Edge>> other_edges;
  vector<PTR<Object<Scalar>>> distances;

  for (int i = 0; i < curves.size(); i++) {
    f_edges.push_back(nullptr);
    other_edges.push_back(nullptr);
    distances.push_back(nullptr);
  }

  PTR<Edge> start = face->e;
  do {
    if (start->poly == f) break;
  } while ((start = start->next) != face->e);

  if (!start) return nullptr;

  for (int i = 1; i < curves.size(); i++) {
    PTR<Object<Poly2>> g = curves[i];

    PTR<Edge> e = start;

    vector<PTR<Edge>> pq_edges;

    do {
      if (e->poly == g || e->poly == f) pq_edges.push_back(e);
    } while ((e = e->next) != start);

    for (int j = 0; j < pq_edges.size(); j++) {
      int k = (j + 1) % pq_edges.size();

      PTR<Edge> p_e = pq_edges[j];
      PTR<Edge> q_e = pq_edges[k];

      if ((p_e->poly != q_e->poly)) {
        if (p_e->poly != f) {
          PTR<Edge> temp = p_e;
          p_e = q_e;
          q_e = temp;
        }

        PTR<Object<Scalar>> d = intersectionDistance(p_e, q_e);

        if (!distances[i] || LessThan(d, distances[i]) < 0) {
          distances[i] = d;
          f_edges[i] = p_e;
          other_edges[i] = q_e;
        }
      }
    }
  }

  int max_i = 1;

  while (max_i < distances.size() && !distances[max_i]) max_i++;

  if (max_i >= distances.size()) return nullptr;

  for (int i = max_i + 1; i < distances.size(); i++) {
    if (distances[i] && LessThan(distances[max_i], distances[i]) < 0) {
      max_i = i;
    }
  }

  PTR<Edge> e_p = f_edges[max_i];
  PTR<Edge> e_q = other_edges[max_i];

  PTR<Object<Poly2>> g = e_q->poly;

  PTR<Object<PV2>> p = e_p->intersection_point();
  PTR<Object<PV2>> q = e_q->intersection_point();

  // Flip signs so that f1(p2) > 0 and f2(p1) > 0.
  PTR<Object<Poly2>> FF = f;
  PTR<Object<Poly2>> GG = g;
  if (PolySign(f, q) < 0) FF = new NegativePoly2(f);
  if (PolySign(g, p) < 0) GG = new NegativePoly2(g);

  // u12 = (u2 - u1).unit(), u21 = -u12.
  PTR<Object<PV2>> u12 = new NormalVector(new DifferencePoint(q, p));
  PTR<Object<PV2>> u21 = new NormalVector(new DifferencePoint(p, q));

  // Unit normals n1 = (grad f1(p1)).unit() and n2 = (grad f2(p2)).unit().
  PTR<Object<PV2>> n1 = new NormalVector(new GradientVector(p, FF));
  PTR<Object<PV2>> n2 = new NormalVector(new GradientVector(q, GG));

  // Inward unit tangents v1 and v2.  (If p1 lies on edge a1b1 in CCW order,
  // (b1-a1) x v1 > 0....)
  PTR<Object<PV2>> v1 = new NormalVector(new TangentVector(p, FF));
  PTR<Object<PV2>> v2 = new NormalVector(new TangentVector(q, GG));

  PTR<Object<PV2>> a1b1 =
      new NormalVector(new DifferencePoint(e_p->head->p, e_p->tail->p));
  PTR<Object<PV2>> a2b2 =
      new NormalVector(new DifferencePoint(e_q->head->p, e_q->tail->p));

  if (Side(a1b1, v1) < 0) v1 = new NegativeVector(v1);
  if (Side(a2b2, v2) < 0) v2 = new NegativeVector(v2);

  PTR<Object<PV2>> P, Q, U, V;

  PTR<Object<Scalar>> u12DotV1 = new DotScalar(u12, v1);
  PTR<Object<Scalar>> u21DotV2 = new DotScalar(u21, v2);

  // Non-skewed
  if (Sign(u12DotV1) == Sign(u21DotV2)) {
    if (LessThan(u12DotV1, u21DotV2) < 0) {
      Q = p;
      U = u12;
    } else {
      Q = q;
      U = u21;
    }

    // printf("splitting: pq or qp\n");

  }
  // Skewed
  else {
    if (LessThan(u12DotV1, u21DotV2) < 0) {
      Q = p;
      U = n1;
      assert(Side(a1b1, n1) > 0);
    } else {
      Q = q;
      U = n2;
      assert(Side(a2b2, n2) > 0);
    }

    // printf("splitting: normal vector\n");
  }

  // Solve for for q+su equidistant from f1 and f2 using linearization about q:

  //(f1(q) + grad f1(q) * u s) / |grad f1(q)| = (f2(q) + grad f2(q) * u s) /
  //|grad f2(q)|

  // Solve for s and set p = q + s u.
  PTR<Object<Scalar>> S = new EquidistantParameter(FF, GG, Q, U);

  // I am assuming that u12 x v1 > 0 and u21 x v2 > 0, but if the curves are
  // really bent, this might not be true.  Put in an assert for now. Fallback
  // perpendicular bisector

  if (Side(u12, v1) != Side(u12, v2) ||
      LessThan(S, new InputScalar(Parameter::constant(0))) < 0) {
    // printf("fail perp\n");

    P = new MidVector(p, q);
    V = new NormalVector(new Rot90(new VectorAB(p, q), nullptr));

    PTR<Line> line_perf = new Line(P, V);

    PTR<Line> line = line_perf->approximate(p, q);

    PTR<Face> other = face->split(line);
    faces.push_back(other);
    other->num_splits = face->num_splits;

    return other;
  }

  Parameter s = S->getApprox(1.0);
  Parameter len = (new DifferencePoint(q, p))->getApprox(1.0).length();
  // printf("s: %g len: %g\n", s.mid(), len.mid());

  P = new SumPoint(Q, new ScaleVector(S, U));

  // Sanity check U * grad ( curve coming out of ) > 0
  if (U == u12) assert(Dot(U, n1) > 0);
  if (U == u21) assert(Dot(U, n2) > 0);

  //(1) Is P inside the cell (only if U isn't the edge of the cell
  //(2) FF(P) and GG(P) is positive
  //    if not div S by 2

  // If U is pointing along the edge between p and q (they're on the same edge)
  // make sure to avoid the inside the cell case

  PTR<Object<Scalar>> pqLength = new VectorLength(new DifferencePoint(q, p));

  bool inside, positive;
  int fail_count = 0;

  do {
    positive = true;
    inside = true;

    PTR<Edge> e = face->e;

    // insde the cell check
    if (e_p->line != e_q->line) {
      inside = face->inside(P);
    } else {
      inside = (LessThan(S, pqLength) < 0);
    }

    // positivity check
    if (PolySign(FF, P) < 0 || PolySign(GG, P) < 0) positive = false;

    Parameter s = S->getApprox(1.0);
    Parameter len = (new DifferencePoint(q, p))->getApprox(1.0).length();
    // printf("s: %g len: %g\n", s.mid(), len.mid());
    // printf("fail_count: %d positive: %s inside: %s\n", fail_count, (positive
    // ? "true" : "false"),
    //                                                                 (inside ?
    //                                                                 "true" :
    //                                                                 "false"));

    if (!positive || !inside) {
      S = new ScaleParameter(0.5, S);
      P = new SumPoint(Q, new ScaleVector(S, U));
      fail_count++;
    }

  } while (!positive || !inside);

  // printf("fail count: %d\n", fail_count);

  // To get v, we solve f1(p+v)=f2(p+v)=0 using linearization about p:
  // f1(p) + grad f1(p) * v = 0
  // f2(p) + grad f2(p) * v = 0
  V = new LinearizationIntersection(FF, GG, P);

  PTR<Line> line_perf = new Line(P, V);

  // If V goes wrong, make V point from P to the midpoint of p and q
  // If V goes wrong, make V the perpendicular bisector
  if (!(line_perf->left(p) != line_perf->left(q))) {
    V = new NormalVector(new Rot90(new VectorAB(p, q), nullptr));
    line_perf = new Line(P, V);
    // printf("fail perp no split\n");
  }

  // printf("succ split\n");

  PTR<Line> line = line_perf->approximate(p, q);

  PTR<Face> other = face->split(line);
  faces.push_back(other);
  other->num_splits = face->num_splits;

  return other;
}

PTR<Object<Scalar>> Encasement::intersectionDistance(PTR<Edge> e1,
                                                     PTR<Edge> e2) {
  PTR<Object<PV2>> p = e1->intersection_point();
  PTR<Object<PV2>> q = e2->intersection_point();

  // Flip signs so that f1(p2) > 0 and f2(p1) > 0.
  PTR<Object<Poly2>> FF = e1->poly;
  PTR<Object<Poly2>> GG = e2->poly;
  if (PolySign(e1->poly, q) < 0) FF = new NegativePoly2(e1->poly);
  if (PolySign(e2->poly, p) < 0) GG = new NegativePoly2(e2->poly);

  // u12 = (u2 - u1).unit(), u21 = -u12.
  PTR<Object<PV2>> u12 = new NormalVector(new DifferencePoint(q, p));
  PTR<Object<PV2>> u21 = new NormalVector(new DifferencePoint(p, q));

  // Unit normals n1 = (grad f1(p1)).unit() and n2 = (grad f2(p2)).unit().
  PTR<Object<PV2>> n1 = new NormalVector(new GradientVector(p, FF));
  PTR<Object<PV2>> n2 = new NormalVector(new GradientVector(q, GG));

  // Inward unit tangents v1 and v2.  (If p1 lies on edge a1b1 in CCW order,
  // (b1-a1) x v1 > 0....)
  PTR<Object<PV2>> v1 = new NormalVector(new TangentVector(p, FF));
  PTR<Object<PV2>> v2 = new NormalVector(new TangentVector(q, GG));

  PTR<Object<PV2>> a1b1 =
      new NormalVector(new DifferencePoint(e1->head->p, e1->tail->p));
  PTR<Object<PV2>> a2b2 =
      new NormalVector(new DifferencePoint(e2->head->p, e2->tail->p));

  if (Side(a1b1, v1) < 0) v1 = new NegativeVector(v1);
  if (Side(a2b2, v2) < 0) v2 = new NegativeVector(v2);

  PTR<Object<Scalar>> u12DotV1 = new DotScalar(u12, v1);
  PTR<Object<Scalar>> u21DotV2 = new DotScalar(u21, v2);

  // Non-skewed
  if (Sign(u12DotV1) == Sign(u21DotV2)) {
    return new VectorLength(new VectorAB(p, q));
  }
  // Skewed
  else {
    if (Sign(u12DotV1) < 0) {
      return new LinearDistanceEstimate(GG, p);
    } else {
      return new LinearDistanceEstimate(FF, q);
    }
  }
}

void Encasement::eliminate_loops(int i) {
  bool horiz = true;

  for (int x = 0; x < sites[i].size(); x++) {
    PTR<Object<PV2>> p = sites[i][x];

    for (int z = 0; z < faces.size(); z++) {
      PTR<Face> f = faces[z];

      // TODO I hope this is OK to remove
      // Parameter::handleIdentity = true;

      if (f->inside(p)) {
        if (STEPS > 0 && count++ >= STEPS) {
          /// printf/gc("count: %d STEPS: %d\t done\n\n", count, STEPS);
          // curve_trace2(f);
          // Parameter::handleIdentity = false;
          return;
        }

        // Parameter::handleIdentity = false;

        // PTR<Object<PV2>> w = new InputPoint(PV2::input(randomNumber(-1, 1),
        // randomNumber(-1, 1)));

        vector<PTR<Object<Scalar>>> bounds = f->bounds();

        double wid = min(
            abs(bounds[2]->getApprox(1.0).x.mid() - p->getApprox(1.0).x.mid()),
            abs(bounds[0]->getApprox(1.0).x.mid() - p->getApprox(1.0).x.mid()));
        double hei = min(
            abs(bounds[3]->getApprox(1.0).x.mid() - p->getApprox(1.0).y.mid()),
            abs(bounds[1]->getApprox(1.0).x.mid() - p->getApprox(1.0).y.mid()));

        PTR<Object<PV2>> w = new InputPoint(
            PV2<Parameter>::input(wid < hei ? 1 : 0, wid < hei ? 0 : 1));

        horiz = !horiz;

        PTR<Face> other = f->split(new Line(p, w));

        // assert(other != NULL);

        // TEST
        if (other) faces.push_back(other);

        break;
      }
    }
  }

  // Parameter::handleIdentity = false;
}

/*
void Encasement::eliminate_loops(int i, int j) {

  vector< vector< PTR<Object<PV2>> > > ssites;
  ssites.push_back(sites[i]);
  ssites.push_back(sites[j]);

  bool horiz = true;

  for(int x = 0; x < ssites.size(); x++) {
    vector< PTR<Object<PV2>> > sites_x = ssites[x];

    for(int y = 0; y < sites_x.size(); y++) {

      PTR<Object<PV2>> p = sites_x[y];

      for(int z = 0; z < faces.size(); z++) {

        PTR<Face> f = faces[z];

        if(f->inside(p)) {

          PTR<Object<PV2>> w = new InputPoint(PV2<Parameter>::input(horiz ? 1 :
0, horiz ? 0 : 1));

          horiz = !horiz;

          PTR<Face> other = f->split(new Line(p, w));

          //assert(other != NULL);

          //TEST
          if(other)
            faces.push_back(other);

          if(STEPS > 0 && count++ >= STEPS) {
            ///printf/gc("count: %d STEPS: %d\t done\n\n", count, STEPS);
            return;
          }

          break;

        }

      }
    }
  }

}
*/

void Encasement::make_sites_verify(vector<PTR<Object<Poly2>>> curves) {
  assert(faces.size() == 1);

  PTR<Object<PV2>> v = new InputPoint(PV2<Parameter>::constant(1, 0));

  PTR<Object<Scalar>> ZERO = new InputScalar(Parameter::constant(0));

  PTR<Face> face = faces[0];

  for (int i = 0; i < curves.size(); i++) {
    ResultantArrangement res;

    PTR<Object<Poly2>> fx = new DerXPoly2(curves[i]);
    PTR<Object<Poly2>> fy = new DerYPoly2(curves[i]);

    vector<PTR<Object<PV2>>> roots = res.getIntersections(fx, fy);

    for (int j = 0; j < roots.size(); j++) {
      PTR<Object<PV2>> p = roots[j];

      if (!face->inside(p)) continue;

      PTR<Line> line = new Line(p, v);

      vector<PTR<Edge>> i_edges = face->intersection_edges(line);

      PTR<Object<Scalar>> neg_r = new LLHit(line, i_edges[0]->line);
      PTR<Object<Scalar>> pos_r = new LLHit(line, i_edges[1]->line);

      if (LessThan(pos_r, neg_r) == -1) {
        PTR<Object<Scalar>> temp = neg_r;
        neg_r = pos_r;
        pos_r = temp;
      }

      assert(i_edges.size() == 2);
      assert(LessThan(neg_r, ZERO) == -1);
      assert(LessThan(ZERO, pos_r) == -1);

      PTR<Object<Poly>> line_poly = new LinePoly(line, curves[i]);

      vector<PTR<Object<Scalar>>> pos_roots =
          PolySolver(line_poly).getRoots(ZERO, pos_r);
      vector<PTR<Object<Scalar>>> neg_roots =
          PolySolver(line_poly).getRoots(neg_r, ZERO);

      if (pos_roots.size() <= 0 || neg_roots.size() <= 0) continue;

      PTR<Object<Scalar>> rr = pos_roots[0];
      PTR<Object<Scalar>> lr = neg_roots[neg_roots.size() - 1];

      PTR<Object<PV2>> l = new LinePoint(line, lr);
      PTR<Object<PV2>> r = new LinePoint(line, rr);

      PV2<Parameter> lp = l->getApprox(1.0);
      PV2<Parameter> rp = r->getApprox(1.0);

      PTR<Object<Scalar>> xc = new InputScalar((lp.x.ub() + rp.x.lb()) / 2);
      PTR<Object<Scalar>> yc =
          new InputScalar(Parameter::constant(p->getApprox(1.0).y.mid()));

      sites[i].push_back(new ParameterPoint(xc, yc));

      // verification

      PTR<Object<PV2>> mid = new ParameterPoint(xc, yc);
      PTR<Object<PV2>> vvv = new InputPoint(PV2<Parameter>::constant(0, 1));

      PTR<Object<PV2>> b = new SumPoint(mid, vvv);
      PTR<Object<PV2>> mid_line = new VectorAB(mid, b);

      PTR<Object<PV2>> l_vec = new VectorAB(mid, l);
      PTR<Object<PV2>> r_vec = new VectorAB(mid, r);

      assert(Side(mid_line, l_vec) != Side(mid_line, r_vec));
    }
  }
}

void Encasement::make_sites(vector<PTR<Object<Poly2>>> curves) {
  vector<map<int, set<Rectangle*>>> v;

  double xl, xu, yl, yu;

  xl = faces[0]->e->tail->p->getApprox(1.0).x.lb();
  xu = faces[0]->e->tail->p->getApprox(1.0).x.ub();
  yl = faces[0]->e->tail->p->getApprox(1.0).y.lb();
  yu = faces[0]->e->tail->p->getApprox(1.0).y.ub();

  double xl2, xu2, yl2, yu2;

  assert(faces.size() == 1);

  PTR<Edge> e = faces[0]->e->next;
  while (e != faces[0]->e) {
    xl2 = e->tail->p->getApprox(1.0).x.lb();
    xu2 = e->tail->p->getApprox(1.0).x.ub();
    yl2 = e->tail->p->getApprox(1.0).y.lb();
    yu2 = e->tail->p->getApprox(1.0).y.ub();

    xl = fmin(xl, xl2);
    xu = fmax(xu, xu2);
    yl = fmin(yl, yl2);
    yu = fmax(yu, yu2);

    e = e->next;
  }

  // Avoid identity for boxes
  xl -= Parameter::delta;
  xu += Parameter::delta;
  yl -= Parameter::delta;
  yu += Parameter::delta;

  // Rectangle::clean();

  for (int i = 0; i < curves.size(); i++) {
    Rectangle* start = new Rectangle(PV2<Parameter>(
        Parameter::interval(xl, xu), Parameter::interval(yl, yu)));

    v.push_back(connected_components(critical_points(curves[i], start)));
  }

  bool horiz = true;

  for (int i = 0; i < v.size(); i++) {
    map<int, set<Rectangle*>> s = v[i];

    for (map<int, set<Rectangle*>>::iterator it = s.begin(); it != s.end();
         it++) {
      set<Rectangle*> r = (*it).second;
      // Rectangle* rect = Rectangle::rects[(*it).first];
      Rectangle* rect = (*(*it).second.begin())->find();

      PTR<Object<PV2>> p =
          new InputPoint(PV2<Parameter>::constant(rect->midX(), rect->midY()));

      sites[i].push_back(p);
    }
  }

  // Rectangle::clean();
}

//(1) Eliminate loops
//  For each curve, verify that for each face does not contain an
//  entire loop. F([x0, x1], [y0, y1]) should not contain 0 where
//  x0, x1, y0, y1 describe the bounding box of the face
void Encasement::eliminate_loops_ellipse(vector<PTR<Object<Poly2>>> curves) {
  // At first there is only the bounding cell
  // Collect all of the ellipses that do not intersect the boundary
  vector<PTR<Object<Poly2>>> list1;
  vector<PTR<Object<Poly2>>> list2;

  list1.insert(list1.begin(), curves.begin(), curves.end());

  bool horiz = true;

  while (true) {
    for (int i = 0; i < list1.size(); i++) {
      bool found = false;
      for (int j = 0; j < faces.size(); j++) {
        if (faces[j]->contains(list1[i])) {
          found = true;
          break;
        }
      }

      if (!found) list2.push_back(list1[i]);
    }

    if (list2.size() == 0) return;

    // Get the first ellipses crit point
    for (int i = list2.size() - 1; i >= 0; i--) {
      PTR<Object<PV2>> cp = new EllipseCriticalPoint(list2[i]);

      int index = -1;
      for (int j = 0; j < faces.size(); j++) {
        if (faces[j]->inside(cp)) {
          index = j;
          break;
        }
      }

      // this ellipse is not within the bounding box at all
      // discard it and find another
      if (index == -1) {
        /*delete cp*/;
        list2.erase(list2.begin() + i);
        continue;
      }

      // PolySign on the critical point to see if sign is opposite of boundary
      // safe sign of PV2 of crit on Poly2 (call get) of each
      //  if not 0, safe to use mid
      //  if 0 do approx, passes between two verts of the critical point's PV2
      //  (horiz line)

      InputPoint* p = new InputPoint(cp->getApprox(1.0).x.mid(),
                                     cp->getApprox(1.0).y.mid());
      InputPoint* v =
          new InputPoint(PV2<Parameter>::input(horiz ? 1 : 0, horiz ? 0 : 1));

      PTR<Line> l = new Line(p, v);

      faces.push_back(faces[index]->split(l));

      break;
    }

    list1.clear();
    list1 = std::move(list2);

    assert(list2.size() == 0);

    // TODO
    // figure out why intersect_edges(list1) doesn't
    // catch every curve......
    intersect_edges(curves);

    horiz = !horiz;
  }
}

//(2) Intersect Edges
//  For each edge, check for intersections with each curve. If one
//  edge intersects multiple curves, split the edge between the
//  intersection point.
void Encasement::intersect_edges(vector<PTR<Object<Poly2>>> curves) {
  for (int i = 0; i < faces.size(); i++) {
    for (int j = 0; j < curves.size(); j++) {
      PTR<Object<Poly2>> f = curves[j];

      PTR<Edge> e = faces[i]->e;

      do {
        e->intersect(f);
      } while ((e = e->next) != faces[i]->e);
    }
  }
}

//(3) Isolate Segments
//  If a single curve intersects a face >= 4 times, split the face
//  between two of the curve segments, check the new splitting edge
//  for intersection with each curve in the cell and update the curves
//  contained in the two new cells. The only one guaranteed to still
//  intersect this face is the one we split between.
PTR<Face> Encasement::isolate_curves(PTR<Face> f, PTR<Object<Poly2>> curve) {
  // printf("Self Isolation\n");

  vector<PTR<Object<PV2>>> points;
  PTR<Edge> e = f->e;
  do {
    if (e->poly == curve) points.push_back(e->intersection_point());
  } while ((e = e->next) != f->e);

  assert(points.size() > 2);

  int index = 0;
  int jndex = 0;
  PTR<Object<Scalar>> smallest = nullptr;

  // for(int i = 1; i < points.size(); i++) {

  for (int i = 0; i < points.size(); i++) {
    for (int j = i + 2; j < points.size(); j += 2) {
      // if((i+2)%points.size() == index) continue;

      // PTR<Object<PV2>> l1_ab = new VectorAB(points[index],
      // points[(index+2)%points.size()]); PTR<Object<PV2>> l2_ab = new
      // VectorAB(points[i], points[(i+2)%points.size()]);

      // PTR<Object<Scalar>> l1 = new VectorLength(l1_ab);
      // PTR<Object<Scalar>> l2 = new VectorLength(l2_ab);

      PTR<Object<PV2>> a = points[i];
      PTR<Object<PV2>> b = points[j];

      PTR<Object<PV2>> ab = new VectorAB(a, b);

      PTR<Object<PV2>> v1 = new TangentVector(a, curve);
      PTR<Object<PV2>> v2 = new TangentVector(b, curve);

      // if(Side(ab, v1) != Side(ab, v2)) v2 = new NegativeVector(v2);

      PTR<Object<Scalar>> s = new DotScalar(v1, v2);

      if (Sign(s) < 0) s = new DotScalar(v1, new NegativeVector(v2));

      // if(i == 1) /////printf/gc("i: 0 length: %g\n", l1.mid());
      /////printf/gc("i: %d length: %g\n", i, l2.mid());
      if (!smallest || LessThan(s, smallest) < 0) {
        index = i;
        jndex = j;
        smallest = s;
      }
    }
  }

  /////printf/gc("index: %d\n", index);

  PTR<Object<PV2>> a = points[index];
  PTR<Object<PV2>> b = points[jndex];

  PTR<Object<PV2>> p = new AverageVector(0.5, a, b);
  // PTR<Object<PV2>>  v = new Rot90(a, b);
  PTR<Object<PV2>> ab = new VectorAB(a, b);
  PTR<Object<PV2>> v1 = new NormalVector(new TangentVector(a, curve));
  PTR<Object<PV2>> v2 = new NormalVector(new TangentVector(b, curve));
  PTR<Object<PV2>> v2_neg = new NegativeVector(v2);

  // TODO change ab

  static int bug_count = 0;
  bug_count++;

  // Don't  need to use ab as the vector to compare to here
  bool neg = (Side(ab, v1) == Side(ab, v2));
  // bool neg = (Dot(v1, v2) == 1);

  PTR<Object<PV2>> v = new AverageVector(0.5, v1, neg ? v2 : v2_neg);

  PTR<Line> l = new Line(p, v);

  PTR<Line> l_approx = l->approximate(a, b);

  PTR<Face> before = f;
  PTR<Face> other = f->split(l_approx);

  if (STEPS > 0 && count == STEPS) {
    // curve_trace2(other);
    // curve_trace2(f);
  }

  if (f == before)
    faces.push_back(other);
  else
    faces.push_back(f);

  return (f == before ? other : f);
}

/*
int Encasement::isolate_intersection(PTR<Face> face, PTR<Object<Poly2>> f,
PTR<Object<Poly2>> g) {

  assert(face->alternates(f, g));

  vector< PTR<Object<PV2>> > points;
  PTR<Object<Poly2>> poly;


  PTR<Edge> e = f->e;
  do {
    if(e->intersection) {
      points.push_back(e->intersection_point());
      poly = e->poly;
    }
    e = e->next;
  } while(e != f->e);

  assert(points.size() == 4);

  PTR<Object<PV2>> bb = new BoundingBox(face);

  //if the last poly seen was f then points[1], [3] are f points
  PTR<GradPoint> v0 = new GradPoint(f, (poly == f ? points[1] : points[0]));
  PTR<GradPoint> v1 = new GradPoint(f, (poly == f ? points[3] : points[2]));
  PTR<GradPoint> w0 = new GradPoint(g, (poly == f ? points[0] : points[1]));
  PTR<GradPoint> w1 = new GradPoint(g, (poly == f ? points[2] : points[3]));

  vector< PTR<GradPoint> > vecs;
  vecs.push_back(v0);
  vecs.push_back(v1);
  vecs.push_back(w0);
  vecs.push_back(w1);

  CCW ccw;

  //sort grad vectors around (1, 0)
  sort(vecs.begin(), vecs.end(), ccw);

  int match = -1;

  //if the cones already intersect trivially, then fail
  //some adjacent pair must have same poly
  for(int i = 0; i < vecs.size(); i++) {
    int j = (i+1)%vecs.size();
    if(vecs[i].same_poly(vecs[j])) {
      match = i;
    }
  }

  if(match < 0) {
    //TODO Fail!
  }

  v0 = vecs[match];
  v1 = vecs[(match+1)%vecs.size()];
  w0 = vecs[(match+2)%vecs.size()];
  w1 = vecs[(match+3)%vecs.size()];

  PTR<Object<PV2>> z0 = new MidVector(v1, w0);
  PTR<Object<PV2>> z0_perp = new Rot90(z0, nullptr);
  PTR<Object<PV2>> z1 = new MidVector(w1, v0);
  PTR<Object<PV2>> z1_perp = new Rot90(z1, nullptr);

  PTR<Object<Scalar>> p0 = new CrossVector(new GradPoint(f, bb), z0);
  PTR<Object<Scalar>> p1 = new CrossVector(new GradPoint(g, bb), z0);
  PTR<Object<Scalar>> q0 = new CrossVector(new GradPoint(f, bb), z0_perp);
  PTR<Object<Scalar>> q1 = new CrossVector(new GradPoint(g, bb), z0_perp);

  int sp0 = Sign(p0);
  int sp1 = Sign(p1);
  int sq0 = Sign(q0);
  int sq1 = Sign(q1);

  if(sp0 == 0 || sp1 == 0 || sq0 == 0 || sq1 == 0) {
    //TODO Fail!
  }

  if(sp0 != sp1 && sq0 == sq1) {
    //TODO success!
  }

  p0 = new CrossVector(new GradPoint(f, bb), z1);
  p1 = new CrossVector(new GradPoint(g, bb), z1);
  q0 = new CrossVector(new GradPoint(f, bb), z1_perp);
  q1 = new CrossVector(new GradPoint(g, bb), z1_perp);

  sp0 = Sign(p0);
  sp1 = Sign(p1);
  sq0 = Sign(q0);
  sq1 = Sign(q1);

  if(sp0 == 0 || sp1 == 0 || sq0 == 0 || sq1 == 0) {
    //TODO Fail!
  }

  if(sp0 != sp1 && sq0 == sq1) {
    //TODO success!
  }

}
*/

int Encasement::isolate_intersection_diamond(PTR<Face> face,
                                             PTR<Object<Poly2>> f,
                                             PTR<Object<Poly2>> g) {
  //(1) Get the point p from Newton's in double, if outside of face, return -1
  PTR<Object<PV2>> p = face->intersection(f, g);

  if (!p) return -1;

  //(2) Evaluate: f_eps = width(f(p)) g_eps = width(g(p))
  //    construct f+f_eps, f-f_eps, g+g_eps, g-g_eps
  PTR<Object<Scalar>> fval = new PolyScalar(f, p);
  PTR<Object<Scalar>> gval = new PolyScalar(g, p);
  double f_eps = 2 * fval->getApprox(1.0).x.intervalWidth();
  double g_eps = 2 * gval->getApprox(1.0).x.intervalWidth();

  int doubled = 0;

  DrawItem item;
  vector<PTR<Object<PV2>>> draw_points;

  PTR<Object<PV2>> p_p_p = nullptr;
  PTR<Object<PV2>> p_p_n = nullptr;
  PTR<Object<PV2>> p_n_p = nullptr;
  PTR<Object<PV2>> p_n_n = nullptr;

  // Run until found valid 4 points
  while (true) {
    if (doubled > 4) break;

    PTR<Object<Poly2>> f_p = new OffsetPoly2(f, f_eps);
    PTR<Object<Poly2>> f_n = new OffsetPoly2(f, -f_eps);
    PTR<Object<Poly2>> g_p = new OffsetPoly2(g, g_eps);
    PTR<Object<Poly2>> g_n = new OffsetPoly2(g, -g_eps);

    //(3) Run Newton's on all 4 combinations for p++, p--, p+-, p-+
    p_p_p = face->intersection(f_p, g_p);
    p_p_n = face->intersection(f_p, g_n);
    p_n_p = face->intersection(f_n, g_p);
    p_n_n = face->intersection(f_n, g_n);

    // failure, Newton's method left the face
    if (!p_p_p || !p_p_n || !p_n_p || !p_n_n) return -1;

    item.types.push_back(POINT);
    item.colors.push_back({0.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(p);
    item.pvs.push_back(draw_points);

    item.types.push_back(POINT);
    item.colors.push_back({1.0f, 0.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(p_p_p);
    item.pvs.push_back(draw_points);

    item.types.push_back(POINT);
    item.colors.push_back({0.0f, 1.0f, 0.0f});
    draw_points.clear();
    draw_points.push_back(p_p_n);
    item.pvs.push_back(draw_points);

    item.types.push_back(POINT);
    item.colors.push_back({0.0f, 0.0f, 1.0f});
    draw_points.clear();
    draw_points.push_back(p_n_p);
    item.pvs.push_back(draw_points);

    item.types.push_back(POINT);
    item.colors.push_back({1.0f, 0.0f, 1.0f});
    draw_points.clear();
    draw_points.push_back(p_n_n);
    item.pvs.push_back(draw_points);

    // double zoom = debugDrawFacePoint(face, f, g, fmax(f_eps, g_eps),
    // p->get().x.mid(), p->getApprox(1.0).y.mid(), item, true); while(zoom < 0)
    // {
    //  zoom = -zoom;
    //  zoom = debugDrawFacePoint(face, f, g, zoom, p->get().x.mid(),
    //  p->getApprox(1.0).y.mid(), item, true);
    //}

    //(4) Verify the sign on each is correct p++ > 0 > 0, p-- < 0 < 0, etc.
    //    If failure, double both epsilons and try again
    if (PolySign(f, p_p_p) != -1 || PolySign(g, p_p_p) != -1 ||
        PolySign(f, p_p_n) != -1 || PolySign(g, p_p_n) != 1 ||
        PolySign(f, p_n_p) != 1 || PolySign(g, p_n_p) != -1 ||
        PolySign(f, p_n_n) != 1 || PolySign(g, p_n_n) != 1) {
      f_eps *= 2;
      g_eps *= 2;
      doubled++;
      p_p_p = nullptr;
      p_p_n = nullptr;
      p_n_p = nullptr;
      p_n_n = nullptr;
      continue;
    }

    break;
  }

  // failure, unable to verify sign
  if (!p_p_p || !p_p_n || !p_n_p || !p_n_n) return -1;

  //(5) Construct a Root2 around the bounding box of the 4 points
  //    even if I passed an PTR<Object<PV2>> that calculates the bounding box of
  //    these 4 points they were all calculated in double and made objects by
  //    Parameter::constant anyway so this is all that can be done, unless we
  //    let the points go outside of the face or we do Newton's in higher
  //    precision

  double rxl = fmin(
      fmin(fmin(p_p_p->getApprox(1.0).x.lb(), p_p_n->getApprox(1.0).x.lb()),
           p_n_p->getApprox(1.0).x.lb()),
      p_n_n->getApprox(1.0).x.lb());
  double rxu = fmax(
      fmax(fmax(p_p_p->getApprox(1.0).x.ub(), p_p_n->getApprox(1.0).x.ub()),
           p_n_p->getApprox(1.0).x.ub()),
      p_n_n->getApprox(1.0).x.ub());
  double ryl = fmin(
      fmin(fmin(p_p_p->getApprox(1.0).y.lb(), p_p_n->getApprox(1.0).y.lb()),
           p_n_p->getApprox(1.0).y.lb()),
      p_n_n->getApprox(1.0).y.lb());
  double ryu = fmax(
      fmax(fmax(p_p_p->getApprox(1.0).y.ub(), p_p_n->getApprox(1.0).y.ub()),
           p_n_p->getApprox(1.0).y.ub()),
      p_n_n->getApprox(1.0).y.ub());

  // diamond check not completely bullet proof, change Root2
  // initialize Root2 off of PV2<Parameter>

  PV2<Parameter> r(Parameter::interval(rxl, rxu),
                   Parameter::interval(ryl, ryu));

  PTR<Object<PV2>> root = new Root2(f, g, r);

  return verify_intersection(face, f, g, root);
}

int Encasement::isolate_intersection_new(PTR<Face> face, PTR<Object<Poly2>> f,
                                         PTR<Object<Poly2>> g) {
  PTR<Object<PV2>> point = face->intersection(f, g);

  // TODO this is not optimal, just because Newton's doesn't
  //     converge to a point in the box doesn't mean we should
  //     stop, continue to try to isolate the gradient cones

  //     Need to implement the gradient cones!
  //     This only draws a small box around the intersection!

  if (!point) return -1;

  int face_count = 0;

  PTR<Object<PV2>> grad_f = new GradientVector(point, f);
  PTR<Object<PV2>> grad_g = new GradientVector(point, g);
  PTR<Object<PV2>> tan_f = new TangentVector(point, f);
  PTR<Object<PV2>> tan_g = new TangentVector(point, g);

  PTR<Object<PV2>> neg_grad_g = new NegativeVector(grad_g);
  PTR<Object<PV2>> neg_tan_g = new NegativeVector(tan_g);

  PTR<Object<PV2>> grad_f_n = new NormalVector(grad_f);
  PTR<Object<PV2>> grad_g_n = new NormalVector(grad_g);
  PTR<Object<PV2>> tan_f_n = new NormalVector(tan_f);
  PTR<Object<PV2>> tan_g_n = new NormalVector(tan_g);

  PTR<Object<PV2>> neg_grad_g_n = new NormalVector(neg_grad_g);
  PTR<Object<PV2>> neg_tan_g_n = new NormalVector(neg_tan_g);

  int sign_grad = Dot(grad_f, grad_g);
  int sign_tan = Dot(tan_f, tan_g);

  PTR<Object<PV2>> v =
      new HalfVector(grad_f_n, sign_grad == 1 ? grad_g_n : neg_grad_g_n);
  PTR<Object<PV2>> w =
      new HalfVector(tan_f_n, sign_tan == 1 ? tan_g_n : neg_tan_g_n);
  PTR<Object<Scalar>> f_p = new PolyScalar(f, point);
  PTR<Object<Scalar>> g_p = new PolyScalar(g, point);

  PTR<Object<Scalar>> grad_f_length = new VectorLength(grad_f);
  PTR<Object<Scalar>> grad_g_length = new VectorLength(grad_g);

  PTR<Object<Scalar>> vert_length = (LessThan(g_p, f_p) == -1 ? f_p : g_p);
  vert_length =
      new DivAB(vert_length,
                (LessThan(grad_f_length, grad_g_length) == -1 ? grad_f_length
                                                              : grad_g_length));

  PTR<Object<Scalar>> horiz_length = new DivAB(
      vert_length,
      new CrossVector(new NormalVector(grad_f), new NormalVector(grad_g)));

  /////printf/gc("\nhoriz_length: %g\n", horiz_length.mid());
  /////printf/gc("vert_length: %g\n\n", vert_length.mid());

  horiz_length = new ScaleParameter(1e15, horiz_length);
  vert_length = new ScaleParameter(1e15, vert_length);

  PTR<Object<PV2>> ap =
      new SumPoint(point, new ScaleVector(horiz_length, new NormalVector(w)));

  PTR<Object<PV2>> bp =
      new SumPoint(point, new ScaleVector(vert_length, new NormalVector(v)));

  PTR<Object<PV2>> cp = new SumPoint(
      point,
      new ScaleVector(new NegativeScalar(horiz_length), new NormalVector(w)));

  PTR<Object<PV2>> dp = new SumPoint(
      point,
      new ScaleVector(new NegativeScalar(vert_length), new NormalVector(v)));

  PTR<Object<PV2>> v_neg = new NegativeVector(v);
  PTR<Object<PV2>> w_neg = new NegativeVector(w);

  // assert(LeftOf(ap, v, point) != LeftOf(cp, v_neg, point));
  // assert(LeftOf(bp, w, point) != LeftOf(dp, w_neg, point));

  PTR<Object<PV2>> a =
      new InputPoint(ap->getApprox(1.0).x.mid(), ap->getApprox(1.0).y.mid());
  PTR<Object<PV2>> b =
      new InputPoint(bp->getApprox(1.0).x.mid(), bp->getApprox(1.0).y.mid());
  PTR<Object<PV2>> c =
      new InputPoint(cp->getApprox(1.0).x.mid(), cp->getApprox(1.0).y.mid());
  PTR<Object<PV2>> d =
      new InputPoint(dp->getApprox(1.0).x.mid(), dp->getApprox(1.0).y.mid());

  PTR<Object<PV2>> v_approx =
      new InputPoint(v->getApprox(1.0).x.mid(), v->getApprox(1.0).y.mid());
  PTR<Object<PV2>> w_approx =
      new InputPoint(w->getApprox(1.0).x.mid(), w->getApprox(1.0).y.mid());
  PTR<Object<PV2>> v_approx_neg = new NegativeVector(v_approx);
  PTR<Object<PV2>> w_approx_neg = new NegativeVector(w_approx);

  PTR<Line> l;

  if (LeftOf(a, v_approx, point) == -1)
    l = new Line(a, v_approx);
  else
    l = new Line(a, v_approx_neg);

  PTR<Face> other;

  other = face->split(l);
  if (other) {
    faces.push_back(other);
    face_count++;
  }

  assert(face->inside(point));

  if (LeftOf(c, v_approx, point) == -1)
    l = new Line(c, v_approx);
  else
    l = new Line(c, v_approx_neg);

  other = face->split(l);
  if (other) {
    faces.push_back(other);
    face_count++;
  }

  assert(face->inside(point));

  if (LeftOf(b, w_approx, point) == -1)
    l = new Line(b, w_approx);
  else
    l = new Line(b, w_approx_neg);

  other = face->split(l);
  if (other) {
    faces.push_back(other);
    face_count++;
  }

  assert(face->inside(point));

  if (LeftOf(d, w_approx, point) == -1)
    l = new Line(d, w_approx);
  else
    l = new Line(d, w_approx_neg);

  other = face->split(l);
  if (other) {
    faces.push_back(other);
    face_count++;
  }

  if (!face->verify(f, g, point)) face_count = (-1 * face_count) - 1;

  face->f = f;
  face->g = g;

  /////printf/gc("\nface_count: %d\n", face_count);
  // assert(face->inside(point));

  // May not alternate if it fails to isolate
  // assert(face->alternates(f, g));

  // Taken care of by the verify function
  // face->intersection_point = point;
  // face->f = f;
  // face->g = g;

  return face_count;
}

int Encasement::refine_narrow(PTR<Face> face, PTR<Object<Poly2>> f) {
  // printf("Refine Narrow\n");

  assert(face->contains(f));
  assert(face->intersectionCount(f) == 2);

  int ret = 0;

  // Top Line---------------------------------------------------------
  // Find a and b
  PTR<Object<PV2>> intersections[2];
  int index = 0;
  PTR<Edge> current = face->e;
  PTR<Edge> e1;
  PTR<Edge> e2;
  do {
    if (current->poly == f) {
      intersections[index++] = current->intersection_point();
      if (index == 1)
        e1 = current;
      else
        e2 = current;
    }
  } while ((current = current->next) != face->e);

  PTR<Object<PV2>> a = intersections[0];
  PTR<Object<PV2>> b = intersections[1];

  // get mid of a and b
  PTR<Object<PV2>> mid = new AverageVector(0.5, a, b);

  // get rot90 of b-a
  PTR<Object<PV2>> v = new VectorAB(a, b);
  PTR<Object<PV2>> perp_v = new Rot90(a, b);

  // Find ip1 and ip2 where rot90 b-a intersects the cell
  // t1 is the sign of ip1s t value and t2 for ip2
  // assert(t1.sign() != t2.sign()) <-- this we know for sure

  PTR<Line> perp_line_real = new Line(mid, perp_v);
  PTR<Line> perp_line = perp_line_real->approximate(a, b);

  PTR<Object<Scalar>> t_vals[2];
  current = face->e;
  index = 0;
  do {
    if (current->intersects(perp_line))
      t_vals[index++] = new LLHit(perp_line, current->line);
  } while ((current = current->next) != face->e);

  assert(index == 2);

  // assert t1 < t2
  bool swap = LessThan(t_vals[0], t_vals[1]) == -1;
  PTR<Object<Scalar>> t1 = swap ? t_vals[0] : t_vals[1];
  PTR<Object<Scalar>> t2 = swap ? t_vals[1] : t_vals[0];

  // assert(t1->get() < 0);
  // assert(t2->get() > 0);

  PTR<Object<Poly>> perp_poly = new LinePoly(perp_line, f);

  // get the roots between t1 and t2
  vector<PTR<Object<Scalar>>> roots = PolySolver(perp_poly).getRoots(t1, t2);

  PTR<Object<Scalar>> curve_hit = roots[roots.size() / 2];
  PTR<Object<PV2>> p = new LinePoint(perp_line, curve_hit);

  // ip1 always corresponds to p, swap if necessary
  // t1 can be zero if the two points lie on the same line
  swap = (t1->getApprox(1.0).x.sign(false) == 0 || Sign(t1) != Sign(curve_hit));
  PTR<Object<Scalar>> temp = t1;
  t1 = swap ? t2 : temp;
  t2 = swap ? temp : t2;

  // assert(t1->get().sign() == curve_hit->get().sign());
  // assert(t2->get().sign() != curve_hit->get().sign());

  PTR<Object<PV2>> ip1 = new LinePoint(perp_line, t1);
  PTR<Object<PV2>> ip2 = new LinePoint(perp_line, t2);

  // get tangent vector v at (p)
  PTR<Object<PV2>> tangent_v = new TangentVector(p, f);
  PTR<Object<PV2>> neg_tangent_v = new NegativeVector(tangent_v);

  // create line (10%(p, q), v) and approximate(p, q)
  PTR<Object<PV2>> p_ip1_10 = new AverageVector(0.5, p, ip1);

  bool neg = (LeftOf(p_ip1_10, tangent_v, mid) == -1);

  PTR<Line> tangent_line = new Line(p_ip1_10, neg ? tangent_v : neg_tangent_v);
  PTR<Line> tangent_line_approx = tangent_line->approximate(p, ip1);

  // split cell with line
  PTR<Face> other = face->split(tangent_line_approx);
  if (other) {
    faces.push_back(other);
    ret = 1;
  }

  if (e1->line != e2->line) {
    // Bottom Line---------------------------------------------------------
    // create line (10%(mid(a, b), ip2), b-a) and approximate(mid(a, b), ip2)
    PTR<Object<PV2>> mid_ip2_10 = new AverageVector(0.5, mid, ip2);

    PTR<Object<PV2>> neg_v = new NegativeVector(v);

    neg = (LeftOf(mid_ip2_10, v, mid) == -1);

    PTR<Line> ab_line = new Line(mid_ip2_10, neg ? v : neg_v);
    PTR<Line> ab_line_approx = ab_line->approximate(mid, ip2);

    // split the cell with line
    other = face->split(ab_line_approx);
    if (other) {
      faces.push_back(other);
      ret = 2;
    }

    /*delete mid_ip2_10*/;
    /*delete ab_line*/;
    /*delete neg_v*/;
  }

  // Clean up---------------------------------------------------------
  /*
  delete t1;
  delete t2;
  delete mid;
  delete v;
  delete perp_v;
  delete perp_line;
  delete a;
  delete b;
  delete perp_poly;
  delete p;
  delete ip1;
  delete ip2;
  delete tangent_v;
  delete p_ip1_10;
  delete tangent_line;
  delete neg_tangent_v;
  */

  return ret;
}

void Encasement::refine_subdiv(PTR<Face> face, PTR<Object<Poly2>> f) {
  // printf("Refine Subdivide\n");

  assert(face->contains(f));
  assert(face->intersectionCount(f) == 2);

  // Find a and b
  PTR<Object<PV2>> intersections[2];
  int index = 0;
  PTR<Edge> current = face->e;
  do {
    if (current->poly == f) {
      intersections[index++] = current->intersection_point();
    }
  } while ((current = current->next) != face->e);

  PTR<Object<PV2>> a = intersections[0];
  PTR<Object<PV2>> b = intersections[1];

  // get mid of a and b
  PTR<Object<PV2>> mid = new AverageVector(0.5, a, b);

  // get rot90 of b-a
  PTR<Object<PV2>> perp_v = new Rot90(a, b);

  // Find ip1 and ip2 where rot90 b-a intersects the cell
  // t1 is the sign of ip1s t value and t2 for ip2
  // assert(t1.sign() != t2.sign()) <-- this we know for sure

  PTR<Line> perp_line = new Line(mid, perp_v);
  PTR<Line> perp_line_approx = perp_line->approximate(a, b);

  PTR<Face> other = face->split(perp_line_approx);
  faces.push_back(other);
}

PTR<Face> Encasement::separate_new(PTR<Face> face, PTR<Object<Poly2>> f,
                                   PTR<Object<Poly2>> g) {
  // printf("Separate new\n");

  DrawItem item;

  if (!face->contains(f) || !face->contains(g)) return nullptr;

  vector<PTR<Object<PV2>>> points;
  vector<PTR<Edge>> intersection_edges;

  vector<int> f_or_g;

  PTR<Edge> e = face->e;
  do {
    if (e->poly == f) {
      points.push_back(e->intersection_point());
      intersection_edges.push_back(e);
      f_or_g.push_back(0);
    }
    if (e->poly == g) {
      points.push_back(e->intersection_point());
      intersection_edges.push_back(e);
      f_or_g.push_back(1);
    }
  } while ((e = e->next) != face->e);

  // printf("f_or_g [");
  // for(int i = 0; i < f_or_g.size(); i++) {
  // printf("%d ", f_or_g[i]);
  //}
  // printf("]\n");

  assert(f_or_g.size() == 4);

  int index = 0;

  while (f_or_g[index] == 1 || (f_or_g[index] == f_or_g[(index + 1) % 4])) {
    index++;
  }

  // printf("index: %d\n", index);

  PTR<Object<PV2>> a1 = points[index];
  PTR<Object<PV2>> b1 = points[(index + 1) % 4];
  PTR<Object<PV2>> b2 = points[(index + 2) % 4];
  PTR<Object<PV2>> a2 = points[(index + 3) % 4];

  PTR<Edge> e_a1 = intersection_edges[index];
  PTR<Edge> e_b1 = intersection_edges[(index + 1) % 4];
  PTR<Edge> e_b2 = intersection_edges[(index + 2) % 4];
  PTR<Edge> e_a2 = intersection_edges[(index + 3) % 4];

  PTR<Object<PV2>> p;
  PTR<Object<PV2>> q;
  PTR<Edge> e_p;
  PTR<Edge> e_q;

  // Flip signs so that f1(p2) > 0 and f2(p1) > 0.
  PTR<Object<Poly2>> FF = f;
  PTR<Object<Poly2>> GG = g;
  if (PolySign(f, b1) < 0) FF = new NegativePoly2(f);
  if (PolySign(g, a1) < 0) GG = new NegativePoly2(g);

  // do this test with f first, and if it fails check g?
  // if e_a1 OR e_a2 have sign 0 of df x dg (check the boolean on the edges)
  // skip

  // TODO
  // Instead of just giving up if the sign is 0 (boolean is true) use the second
  // derivative test in the paper

  if (!e_a1->gradient_cross_point && !e_a2->gradient_cross_point) {
    int sa1 = Side(new GradientVector(a1, FF), new GradientVector(a1, GG));
    int sa2 = Side(new GradientVector(a2, FF), new GradientVector(a2, GG));

    // If they have opposite signs then there is either a minimum or maximum on
    // f in between a and b find it and split at it in the direction of the
    // common gradient else just fall through to choosing the CLOSEST pair of
    // points to generate a splitting line

    // if it's a min
    if (sa1 == -1 && sa2 == 1) {
      PTR<Object<PV2>> gp = face->closest_approach(FF, GG, p);

      if (gp) {
        // cout << "Found min" << endl;

        PTR<Object<PV2>> grad_gp = new GradientVector(gp, FF);

        PTR<Line> line_perf = new Line(gp, grad_gp);

        PTR<Line> line;

        // If the line doesn't go between a1 and a2, use maximal vertices
        if (line_perf->left(a1) == line_perf->left(a2)) {
          vector<PTR<Vertex>> left_v;
          vector<PTR<Vertex>> right_v;

          vector<PTR<Vertex>> verts = face->get_verts();

          std::copy_if(
              verts.begin(), verts.end(), std::back_inserter(left_v),
              [&](PTR<Vertex> v) { return line_perf->left(v->p) < 0; });
          std::copy_if(
              verts.begin(), verts.end(), std::back_inserter(right_v),
              [&](PTR<Vertex> v) { return line_perf->left(v->p) > 0; });

          vector<PTR<Object<Scalar>>> left_d;
          vector<PTR<Object<Scalar>>> right_d;

          std::transform(
              left_v.begin(), left_v.end(), std::back_inserter(left_d),
              [&](PTR<Vertex> a) { return line_perf->distance(a->p); });
          std::transform(
              right_v.begin(), right_v.end(), std::back_inserter(right_d),
              [&](PTR<Vertex> a) { return line_perf->distance(a->p); });

          std::vector<PTR<Object<Scalar>>>::iterator lmax = std::max_element(
              left_d.begin(), left_d.end(),
              [&](PTR<Object<Scalar>> a, PTR<Object<Scalar>> b) {
                return LessThan(a, b) == 1;
              });
          std::vector<PTR<Object<Scalar>>>::iterator rmax = std::max_element(
              right_d.begin(), right_d.end(),
              [&](PTR<Object<Scalar>> a, PTR<Object<Scalar>> b) {
                return LessThan(a, b) == 1;
              });

          PTR<Vertex> lv = left_v[std::distance(left_d.begin(), lmax)];
          PTR<Vertex> rv = right_v[std::distance(right_d.begin(), rmax)];

          line = line_perf->approximate(lv->p, rv->p);

        } else {
          line = line_perf->approximate(a1, a2);
        }

        PTR<Face> other = face->split(line);
        faces.push_back(other);
        other->num_splits = face->num_splits;

        // TODO go through the edges generated by the split and find the one
        // with intersection with f with smallest t, mark this point on the edge
        // (and it's neighbor!) as having zero gradient cross product
        vector<PTR<Edge>> check_edges;
        PTR<Edge> edge = face->e;
        while (edge->line == line) {
          check_edges.push_back(edge);
          edge = edge->next;
        }

        int index = -1;
        for (int i = 0; i < check_edges.size(); i++) {
          if (check_edges[i]->intersection && check_edges[i]->poly == f) {
            if (index == -1 ||
                LessThanAbs(check_edges[i]->intersection,
                            check_edges[index]->intersection) < 0) {
              index = i;
            }
          }
        }

        assert(index != -1);

        // This is important, mark this point on the edge as having sign
        // gradient cross == 0
        check_edges[index]->gradient_cross_point =
            check_edges[index]->twin->gradient_cross_point = true;

        return other;

      }
      // FAILURE of Newton's (take note)
      else {
        // cout << "min finding failed" << endl;
      }
    }
  }
  // else if(!e_b1->gradient_cross_point && e_b2->gradient_cross_point) {
  //
  //}

  if (Shorter(a1, b1, a2, b2) == 1) {
    p = a2;
    q = b2;
    e_p = e_a2;
    e_q = e_b2;
  } else {
    p = a1;
    q = b1;
    e_p = e_a1;
    e_q = e_b1;
  }

  vector<PTR<Object<PV2>>> draw_points;

  //  item.types.push_back(LINE);
  //  item.colors.push_back({1.0f, 1.0f, 0.0f});
  //  draw_points.clear();
  //  draw_points.push_back(p);
  //  draw_points.push_back(q);
  //  item.pvs.push_back(draw_points);

  //  item.types.push_back(POINT);
  //  item.colors.push_back({1.0f, 0.0f, 0.0f});
  //  draw_points.clear();
  //  draw_points.push_back(p);
  //  item.pvs.push_back(draw_points);
  //
  //  item.types.push_back(POINT);
  //  item.colors.push_back({0.0f, 1.0f, 0.0f});
  //  draw_points.clear();
  //  draw_points.push_back(q);
  //  item.pvs.push_back(draw_points);

  // face->closest_approach(new AddPoly2(FF, new NegativePoly2(GG)), GG, p);
  // face->closest_approach(FF, GG, p);

  // u12 = (u2 - u1).unit(), u21 = -u12.
  PTR<Object<PV2>> u12 = new NormalVector(new DifferencePoint(q, p));
  PTR<Object<PV2>> u21 = new NormalVector(new DifferencePoint(p, q));

  // Unit normals n1 = (grad f1(p1)).unit() and n2 = (grad f2(p2)).unit().
  PTR<Object<PV2>> n1 = new NormalVector(new GradientVector(p, FF));
  PTR<Object<PV2>> n2 = new NormalVector(new GradientVector(q, GG));

  // Inward unit tangents v1 and v2.  (If p1 lies on edge a1b1 in CCW order,
  // (b1-a1) x v1 > 0....)
  PTR<Object<PV2>> v1 = new NormalVector(new TangentVector(p, FF));
  PTR<Object<PV2>> v2 = new NormalVector(new TangentVector(q, GG));

  PTR<Object<PV2>> a1b1 =
      new NormalVector(new DifferencePoint(e_p->head->p, e_p->tail->p));
  PTR<Object<PV2>> a2b2 =
      new NormalVector(new DifferencePoint(e_q->head->p, e_q->tail->p));

  if (Side(a1b1, v1) < 0) {
    v1 = new NegativeVector(v1);
  }

  if (Side(a2b2, v2) < 0) {
    v2 = new NegativeVector(v2);
  }

  PTR<Object<PV2>> P, Q, U, V;

  PTR<Object<Scalar>> u12DotV1 = new DotScalar(u12, v1);
  PTR<Object<Scalar>> u21DotV2 = new DotScalar(u21, v2);

  // Non-skewed
  if (Sign(u12DotV1) == Sign(u21DotV2)) {
    if (LessThan(u12DotV1, u21DotV2) < 0) {
      Q = p;
      U = u12;
    } else {
      Q = q;
      U = u21;
    }
  }
  // Skewed
  else {
    if (LessThan(u12DotV1, u21DotV2) < 0) {
      Q = p;
      U = n1;
    } else {
      Q = q;
      U = n2;
    }
  }

  //  item.types.push_back(VECTOR);
  //  item.colors.push_back({0.0f, 0.0f, 0.0f});
  //  draw_points.clear();
  //  draw_points.push_back(Q);
  //  draw_points.push_back(new SumPoint(Q, U));
  //  item.pvs.push_back(draw_points);
  //
  //  item.types.push_back(POINT);
  //  item.colors.push_back({1.0f, 0.0f, 0.0f});
  //  draw_points.clear();
  //  draw_points.push_back(Q);
  //  item.pvs.push_back(draw_points);
  //
  //  item.types.push_back(LINE);
  //  item.colors.push_back({1.0f, 1.0f, 0.0f});
  //  draw_points.clear();
  //  draw_points.push_back(p);
  //  draw_points.push_back(q);
  //  item.pvs.push_back(draw_points);
  //
  //  debugDrawFace(face, f, g, 1, 0, 0, item, true);

  // Solve for for q+su equidistant from f1 and f2 using linearization about q:

  //(f1(q) + grad f1(q) * u s) / |grad f1(q)| = (f2(q) + grad f2(q) * u s) /
  //|grad f2(q)|

  // Solve for s and set p = q + s u.
  PTR<Object<Scalar>> S = new EquidistantParameter(FF, GG, Q, U);

  // I am assuming that u12 x v1 > 0 and u21 x v2 > 0, but if the curves are
  // really bent, this might not be true.  Put in an assert for now. Fallback
  // perpendicular bisector

  if ((U == u12 && Dot(U, n1) < 0) || (U == u21 && Dot(U, n2) < 0) ||
      Side(u12, v1) != Side(u12, v2) ||
      LessThan(S, new InputScalar(Parameter::constant(0))) < 0) {
    // printf("fail perp\n");

    P = new MidVector(p, q);
    V = new NormalVector(new Rot90(new VectorAB(p, q), nullptr));

    PTR<Line> line_perf = new Line(P, V);

    PTR<Line> line = line_perf->approximate(p, q);

    PTR<Face> other = face->split(line);
    faces.push_back(other);
    other->num_splits = face->num_splits;

    if (STEPS > 0 && count == STEPS) {
      // curve_trace2(other);
      // curve_trace2(face);
    }

    return other;
  }

  Parameter s = S->getApprox(1.0);
  Parameter len = (new DifferencePoint(q, p))->getApprox(1.0).length();

  P = new SumPoint(Q, new ScaleVector(S, U));

  PTR<LinearizationPoly> flin = new LinearizationPoly(FF, Q);
  PTR<LinearizationPoly> glin = new LinearizationPoly(GG, Q);
  PTR<Line> fline = flin->getLine();
  PTR<Line> gline = glin->getLine();

  PTR<Object<Scalar>> NEG = new InputScalar(-10);
  PTR<Object<Scalar>> POS = new InputScalar(10);

  item.types.push_back(LINE);
  if (f == curves[0])
    item.colors.push_back({1.0f, 0.0f, 0.0f});
  else
    item.colors.push_back({0.0f, 1.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(new LinePoint(fline, NEG));
  draw_points.push_back(new LinePoint(fline, POS));
  item.pvs.push_back(draw_points);

  item.types.push_back(LINE);
  if (f == curves[0])
    item.colors.push_back({0.0f, 1.0f, 0.0f});
  else
    item.colors.push_back({1.0f, 0.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(new LinePoint(gline, NEG));
  draw_points.push_back(new LinePoint(gline, POS));
  item.pvs.push_back(draw_points);

  item.types.push_back(POINT);
  item.colors.push_back({0.0f, 0.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(p);
  item.pvs.push_back(draw_points);

  item.types.push_back(POINT);
  item.colors.push_back({0.0f, 0.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(q);
  item.pvs.push_back(draw_points);

  item.types.push_back(LINE);
  item.colors.push_back({1.0f, 1.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(p);
  draw_points.push_back(q);
  item.pvs.push_back(draw_points);

  item.types.push_back(VECTOR);
  item.colors.push_back({0.0f, 0.0f, 0.f});
  draw_points.clear();
  draw_points.push_back(Q);
  draw_points.push_back(P);
  item.pvs.push_back(draw_points);

  item.types.push_back(POINT);
  item.colors.push_back({1.0f, 0.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(P);
  item.pvs.push_back(draw_points);

  // printf("Equidistant Point\n");
  // debugDrawFace(face, FF, GG, 0.8, 0, 0, item, true);

  // Sanity check U * grad ( curve coming out of ) > 0
  if (U == u12) assert(Dot(U, n1) > 0);
  if (U == u21) assert(Dot(U, n2) > 0);

  //(1) Is P inside the cell (only if U isn't the edge of the cell
  //(2) FF(P) and GG(P) is positive
  //    if not div S by 2

  // If U is pointing along the edge between p and q (they're on the same edge)
  // make sure to avoid the inside the cell case

  PTR<Object<Scalar>> pqLength = new VectorLength(new DifferencePoint(q, p));

  bool inside, positive;
  int fail_count = 0;
  do {
    positive = true;
    inside = true;

    PTR<Edge> e = face->e;

    if (e_p->line != e_q->line) {
      inside = face->inside(P);
    } else {
      inside = (LessThan(S, pqLength) < 0);
    }

    // positivity check
    if (PolySign(FF, P) < 0 || PolySign(GG, P) < 0) positive = false;

    if (!positive || !inside) {
      S = new ScaleParameter(0.5, S);
      P = new SumPoint(Q, new ScaleVector(S, U));
      fail_count++;
    }

  } while (!positive || !inside);

  // printf("fail count: %d\n", fail_count);

  // To get v, we solve f1(p+v)=f2(p+v)=0 using linearization about p:
  // f1(p) + grad f1(p) * v = 0
  // f2(p) + grad f2(p) * v = 0
  V = new LinearizationIntersection(FF, GG, P);

  PTR<Line> line_perf = new Line(P, V);

  flin = new LinearizationPoly(FF, P);
  glin = new LinearizationPoly(GG, P);
  fline = flin->getLine();
  gline = glin->getLine();

  item.types.push_back(LINE);
  if (f == curves[0])
    item.colors.push_back({1.0f, 0.0f, 0.0f});
  else
    item.colors.push_back({0.0f, 1.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(new LinePoint(fline, NEG));
  draw_points.push_back(new LinePoint(fline, POS));
  item.pvs.push_back(draw_points);

  item.types.push_back(LINE);
  if (f == curves[0])
    item.colors.push_back({0.0f, 1.0f, 0.0f});
  else
    item.colors.push_back({1.0f, 0.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(new LinePoint(gline, NEG));
  draw_points.push_back(new LinePoint(gline, POS));
  item.pvs.push_back(draw_points);

  item.types.push_back(LINE);
  item.colors.push_back({1.0f, 1.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(p);
  draw_points.push_back(q);
  item.pvs.push_back(draw_points);

  item.types.push_back(VECTOR);
  item.colors.push_back({0.0f, 0.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(P);
  draw_points.push_back(new SumPoint(P, V));
  item.pvs.push_back(draw_points);

  item.types.push_back(POINT);
  item.colors.push_back({1.0f, 0.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(P);
  item.pvs.push_back(draw_points);

  // printf("Splitting Line\n");
  // debugDrawFace(face, f, g, 0.8, 0, 0, item, true);

  // If V goes wrong, make V point from P to the midpoint of p and q
  // If V goes wrong, make V the perpendicular bisector
  if (!(line_perf->left(p) != line_perf->left(q))) {
    V = new NormalVector(new Rot90(new VectorAB(p, q), nullptr));
    P = new MidVector(p, q);
    line_perf = new Line(P, V);

    // printf("fail perp no split\n");
  }

  item.types.push_back(VECTOR);
  item.colors.push_back({0.0f, 0.0f, 1.0f});
  draw_points.clear();
  draw_points.push_back(P);
  // draw_points.push_back(new SumPoint(P, new ScaleVector(0.5, new
  // NormalVector(V))));
  draw_points.push_back(new SumPoint(P, V));
  item.pvs.push_back(draw_points);

  item.types.push_back(VECTOR);
  item.colors.push_back({0.0f, 0.0f, 1.0f});
  draw_points.clear();
  draw_points.push_back(P);
  // draw_points.push_back(new SumPoint(P, new ScaleVector(0.5, new
  // NormalVector(new NegativeVector(V)))));
  draw_points.push_back(new SumPoint(P, new NegativeVector(V)));
  item.pvs.push_back(draw_points);

  item.types.push_back(POINT);
  item.colors.push_back({0.0f, 0.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(P);
  item.pvs.push_back(draw_points);
  // printf("succ split\n");
  // debugDrawFace(face, f, g, 1, 0, 0, item, true);

  PTR<Line> line = line_perf->approximate(p, q);

  PTR<Face> other = face->split(line);
  faces.push_back(other);
  other->num_splits = face->num_splits;

  return other;
}

PTR<Face> Encasement::separate(PTR<Face> face, PTR<Object<Poly2>> f,
                               PTR<Object<Poly2>> g) {
  if (!face->contains(f) || !face->contains(g)) return nullptr;

  vector<PTR<Object<PV2>>> points;

  vector<int> f_or_g;

  PTR<Edge> e = face->e;
  do {
    if (e->poly == f) {
      points.push_back(e->intersection_point());
      f_or_g.push_back(0);
    }
    if (e->poly == g) {
      points.push_back(e->intersection_point());
      f_or_g.push_back(1);
    }
  } while ((e = e->next) != face->e);

  // printf("f_or_g [");
  // for(int i = 0; i < f_or_g.size(); i++) {
  // printf("%d ", f_or_g[i]);
  //}
  // printf("]\n");

  assert(f_or_g.size() == 4);

  int index = 0;

  while (f_or_g[index] == 1 || (f_or_g[index] == f_or_g[(index + 1) % 4])) {
    index++;
  }

  // printf("index: %d\n", index);

  PTR<Object<PV2>> a1 = points[index];
  PTR<Object<PV2>> b1 = points[(index + 1) % 4];
  PTR<Object<PV2>> b2 = points[(index + 2) % 4];
  PTR<Object<PV2>> a2 = points[(index + 3) % 4];

  PTR<Object<PV2>> p;
  PTR<Object<PV2>> q;

  if (Shorter(a1, b1, a2, b2) == 1) {
    p = a2;
    q = b2;
  } else {
    p = a1;
    q = b1;
  }

  PTR<Object<Scalar>> t = middleT(p, f, q, g);
  PTR<Object<PV2>> r = new AveragePoint(t, p, q);
  PTR<Object<PV2>> v = middleV(f, g, r, p, q);

  PTR<Object<PV2>> rot_v = new Rot90(v, nullptr);

  PTR<Line> line_perf = new Line(r, v);

  PTR<Line> line = line_perf->approximate(p, q);

  /*
    PV2 pq = (q->get() - p->get()).unit();

    Parameter l1 = pq.dot(v->get().unit()).abs();
    Parameter l2 = pq.dot(rot_v->get().unit()).abs();

    ///printf/gc("v    : (%g, %g)\n", v->get().x.mid(),
    v->getApprox(1.0).y.mid());
    ///printf/gc("rot_v: (%g, %g)\n", rot_v->get().x.mid(),
    rot_v->getApprox(1.0).y.mid());

    ///printf/gc("q - p: (%g, %g)\n", pq.x.mid(), pq.y.mid());

    ///printf/gc("l1: %g\n", l1.mid());
    ///printf/gc("l2: %g\n", l2.mid());

    if(l1 < l2) {
      line = new Line(r, v);
      delete rot_v;
    }
    else {
      line = new Line(r, rot_v);
      delete v;
    }
  */
  PTR<Face> other = face->split(line);
  faces.push_back(other);
  other->num_splits = face->num_splits;
  return other;
}

bool Encasement::add_to_queues(PTR<Face> face, stack<Item>& nq, stack<Item>& q,
                               stack<Item>& sq) {
  // take same curve intersections first
  for (set<PTR<Object<Poly2>>>::iterator it = face->curves.begin();
       it != face->curves.end(); it++) {
    PTR<Object<Poly2>> f = *it;
    if (face->intersectionCount(f) > 2) {
      face->num_splits = 0;
      face->f = nullptr;
      face->g = nullptr;
      sq.push({face, f, nullptr, false});
      /////printf/gc("\n0x%lx added to self queue\n", (long)face);

      return true;
    }
  }

  if (face->num_splits > 0) {
    /// printf/gc("\nsplitting face one more time\n");

    if (!face->contains(face->f) || !face->contains(face->g)) {
      if (face->contains(face->f) && face->intersectionCount(face->f) > 2) {
        face->num_splits = 0;
        face->f = nullptr;
        face->g = nullptr;
        sq.push({face, face->f, nullptr, false});
        /////printf/gc("\n0x%lx added to self queue\n", (long)face);
        return true;
      } else if (face->contains(face->g) &&
                 face->intersectionCount(face->g) > 2) {
        face->num_splits = 0;
        face->f = nullptr;
        face->g = nullptr;
        sq.push({face, face->g, nullptr, false});
        /////printf/gc("\n0x%lx added to self queue\n", (long)face);
        return true;

      } else {
        return true;
      }
    }

    if (!face->alternates(face->f, face->g))
      nq.push({face, face->f, face->g, false});
    else
      q.push({face, face->f, face->g, true, "num_splits"});
    return true;
  }

  for (set<PTR<Object<Poly2>>>::iterator it = face->curves.begin();
       it != face->curves.end(); it++) {
    for (set<PTR<Object<Poly2>>>::iterator jt = it;
         ++jt != face->curves.end();) {
      PTR<Object<Poly2>> f = *it;
      PTR<Object<Poly2>> g = *jt;

      if (!face->contains(f) || !face->contains(g)) continue;

      if (!face->alternates(f, g)) {
        nq.push({face, f, g, false});
        /////printf/gc("\n0x%lx added to non inter queue\n", (long)face);
      } else {
        /////printf/gc("\n0x%lx added to inter queue\n", (long)face);
        q.push({face, f, g, face->alternates(f, g), "loop"});
      }

      // only add face to queue once, because if will be split
      return true;
    }
  }

  return false;
}

PTR<Face> Encasement::shrinkRoot(PTR<Face> f) {
  Parameter zero = Parameter::constant(0);
  Parameter one = Parameter::constant(1);

  PTR<Object<PV2>> vert = new InputPoint(PV2<Parameter>(zero, one));
  PTR<Object<PV2>> horiz = new InputPoint(PV2<Parameter>(one, zero));

  // find the x-y bounds of the box and split along longer axis
  std::vector<PTR<Object<Scalar>>> bounds = f->bounds();

  PTR<Object<Scalar>> xmag = new DiffAB(bounds[2], bounds[0]);
  PTR<Object<Scalar>> ymag = new DiffAB(bounds[3], bounds[1]);

  bool h;

  // split along the longer axis
  h = (LessThan(xmag, ymag) == -1);

  // printf("%s split\n", h ? "horiz" : "vert");

  // middle point of bounding rect, p for the splitting line regardless of horiz
  // or vert will be randomized

  // double r1 = (rand() * 0.6 / RAND_MAX) + 0.2;
  // double r2 = (rand() * 0.6 / RAND_MAX) + 0.2;

  double r1 = randomNumber(0.2, 0.8);
  double r2 = randomNumber(0.2, 0.8);

  PTR<Object<PV2>> p =
      new ParameterPoint(approxAverage(r1, bounds[0], bounds[2]),
                         approxAverage(r2, bounds[1], bounds[3]));

  PTR<Line> line = new Line(p, h ? horiz : vert);

  /*
  PTR<Line> lapprox = line->approximate(new ParameterPoint(bounds[0],
  bounds[1]), new ParameterPoint(bounds[2], bounds[3]));


  p = lapprox->p;
  */

  // p = new InputPoint(p->get().x.mid(), p->getApprox(1.0).y.mid());

  PTR<Line> split = new Line(p, h ? horiz : vert, h ? 0 : 1);

  PTR<Object<Poly2>> f_poly = f->f;
  PTR<Object<Poly2>> g_poly = f->g;

  PTR<Face> split_face = f->split(split);

  if (!split_face) {
    printf("error no split\n");
    printf("p (%g, %g)\n", p->getApprox(1.0).x.mid(),
           p->getApprox(1.0).y.mid());
    bad_faces.first = f;
    bad_faces.second = nullptr;
    return 0;
  }

  split_face->f = f->f = f_poly;
  split_face->g = f->g = g_poly;

  // box = box that still has the intersection
  bool falt = f->contains(f_poly) && f->contains(g_poly) &&
              f->alternates(f_poly, g_poly);
  bool salt = split_face->contains(f_poly) && split_face->contains(g_poly) &&
              split_face->alternates(f_poly, g_poly);

  if (falt == salt) {
    printf("both faces %s alternating\n", falt ? "" : "not");
    bad_faces.first = f;
    bad_faces.second = split_face;
    return 0;
  }

  assert(falt != salt);

  if (salt) f = split_face;

  return f;
}

// fix the faces vector thing where the wrong face is in the faces vector
// fix the no f and g in the face assertion error when getting roots?
// roots to g as well as f, whichever you hit first

int Encasement::isolate_intersection(PTR<Face> face, PTR<Object<Poly2>> f,
                                     PTR<Object<Poly2>> g) {
  return isolate_intersection_diamond(face, f, g);

  static int i_count = 0;

  i_count++;

  // printf("count: %d\n", i_count);

  face->f = f;
  face->g = g;

  if (!face->alternates(f, g)) {
    // printf("new isolation: face not alternating\n");
    return -1;
  }

  PTR<Face> original = face;

  // TODO if I take this out, then I have to check fxXfy on the BB of this face,
  //     and split the face until a Newton step can happen
  //  face = face->clone();
  //
  //  while(!face->rect()) {
  //    PTR<Face> other = shrinkRoot(face);
  //    if(other == nullptr) {
  //      //printf("other null here\n");
  //    }
  //    if(!other->alternates(f, g)) {
  //      //printf("other !alternates here\n");
  //    }
  //
  //    assert(other);
  //    face = other;
  //  }

  // create a Root2 out of the resulting PV2 and add it to roots
  std::vector<PTR<Object<Scalar>>> bounds = face->bounds();

  // TODO construct r with Parameter bounds

  // PV2 r(Parameter::interval(bounds[0]->get().lb(), bounds[2]->get().ub()),
  //      Parameter::interval(bounds[1]->get().lb(), bounds[3]->get().ub()));

  Parameter rxl = bounds[0]->getApprox(1.0).x.lbP();
  Parameter rxu = bounds[2]->getApprox(1.0).x.ubP();
  Parameter ryl = bounds[1]->getApprox(1.0).x.lbP();
  Parameter ryu = bounds[3]->getApprox(1.0).x.ubP();

  PV2<Parameter> r(rxl.interval(rxu), ryl.interval(ryu));

  PTR<Object<PV2>> root = new Root2(f, g, r);

  face = original;

  return verify_intersection(face, f, g, root);
}

// get average tangent vector at this root
// trace to edge of box in 4 directions
// substitute f and g into the line and make sure there are no roots
// other than the intersection point inside the box
// if there are take the two roots closest on each side of the root2
//   and cut with lines in between the root2 and the root of the line poly
int Encasement::verify_intersection(PTR<Face> face, PTR<Object<Poly2>> f,
                                    PTR<Object<Poly2>> g,
                                    PTR<Object<PV2>> root) {
  PTR<Object<Scalar>> ZERO = new InputScalar(Parameter::constant(0));

  PTR<Object<PV2>> f_t = new NormalVector(new TangentVector(root, f));
  PTR<Object<PV2>> g_t = new NormalVector(new TangentVector(root, g));

  if (Dot(f_t, g_t) == -1) {
    g_t = new NegativeVector(g_t);
    assert(Dot(f_t, g_t) == 1);
  }

  PTR<Object<PV2>> avg_t = new SumPoint(f_t, g_t);
  PTR<Object<PV2>> rot_t = new Rot90(avg_t, nullptr);

  PTR<Object<PV2>> avg_t_neg = new NegativeVector(avg_t);
  PTR<Object<PV2>> rot_t_neg = new NegativeVector(rot_t);

  static int assert_count = 0;

  // printf("assert: %d\n", ++assert_count);

  assert(CCWOrder(f_t, g_t) != 0);

  PTR<Line> avg_t_line = new Line(root, avg_t);
  PTR<Line> rot_t_line = new Line(root, rot_t);

  PTR<Object<Poly>> avg_t_poly_f =
      new ShiftOutPoly(new LinePoly(avg_t_line, f));
  PTR<Object<Poly>> rot_t_poly_f =
      new ShiftOutPoly(new LinePoly(rot_t_line, f));

  PTR<Object<Poly>> avg_t_poly_g =
      new ShiftOutPoly(new LinePoly(avg_t_line, g));
  PTR<Object<Poly>> rot_t_poly_g =
      new ShiftOutPoly(new LinePoly(rot_t_line, g));

  PTR<Object<Scalar>> avg_t_hit[2];  //[0] == negative hit, [1] == positive hit
  PTR<Object<Scalar>> rot_t_hit[2];

  PTR<Edge> e = face->e;
  int pos_count = 0;
  int neg_count = 0;
  do {
    if (e->intersects(avg_t_line)) {
      PTR<Object<Scalar>> hit = new LLHit(avg_t_line, e->line);
      int sign = Sign(hit);
      avg_t_hit[(sign == -1 ? 0 : 1)] = hit;
      (sign == -1 ? neg_count++ : pos_count++);
    }

    if (e->intersects(rot_t_line)) {
      PTR<Object<Scalar>> hit = new LLHit(rot_t_line, e->line);
      int sign = Sign(hit);
      rot_t_hit[(sign == -1 ? 0 : 1)] = hit;
      (sign == -1 ? neg_count++ : pos_count++);
    }

  } while ((e = e->next) != face->e);

  assert(pos_count == 2 && neg_count == 2);

  vector<PTR<Object<Scalar>>> avg_t_pos_roots_f =
      PolySolver(avg_t_poly_f).getRoots(ZERO, avg_t_hit[1]);
  vector<PTR<Object<Scalar>>> avg_t_neg_roots_f =
      PolySolver(avg_t_poly_f).getRoots(avg_t_hit[0], ZERO);
  vector<PTR<Object<Scalar>>> rot_t_pos_roots_f =
      PolySolver(rot_t_poly_f).getRoots(ZERO, rot_t_hit[1]);
  vector<PTR<Object<Scalar>>> rot_t_neg_roots_f =
      PolySolver(rot_t_poly_f).getRoots(rot_t_hit[0], ZERO);

  vector<PTR<Object<Scalar>>> avg_t_pos_roots_g =
      PolySolver(avg_t_poly_g).getRoots(ZERO, avg_t_hit[1]);
  vector<PTR<Object<Scalar>>> avg_t_neg_roots_g =
      PolySolver(avg_t_poly_g).getRoots(avg_t_hit[0], ZERO);
  vector<PTR<Object<Scalar>>> rot_t_pos_roots_g =
      PolySolver(rot_t_poly_g).getRoots(ZERO, rot_t_hit[1]);
  vector<PTR<Object<Scalar>>> rot_t_neg_roots_g =
      PolySolver(rot_t_poly_g).getRoots(rot_t_hit[0], ZERO);

  DrawItem item;
  vector<PTR<Object<PV2>>> draw_points;

  item.types.push_back(LINE);
  item.colors.push_back({1.0f, 1.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(root);
  draw_points.push_back(
      new SumPoint(root, new ScaleVector(avg_t_hit[1], avg_t)));
  item.pvs.push_back(draw_points);

  item.types.push_back(LINE);
  item.colors.push_back({1.0f, 1.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(root);
  draw_points.push_back(
      new SumPoint(root, new ScaleVector(avg_t_hit[0], avg_t)));
  item.pvs.push_back(draw_points);

  item.types.push_back(LINE);
  item.colors.push_back({0.0f, 0.0f, 1.0f});
  draw_points.clear();
  draw_points.push_back(root);
  draw_points.push_back(
      new SumPoint(root, new ScaleVector(rot_t_hit[1], rot_t)));
  item.pvs.push_back(draw_points);

  item.types.push_back(LINE);
  item.colors.push_back({0.0f, 0.0f, 1.0f});
  draw_points.clear();
  draw_points.push_back(root);
  draw_points.push_back(
      new SumPoint(root, new ScaleVector(rot_t_hit[0], rot_t)));
  item.pvs.push_back(draw_points);

  item.types.push_back(POINT);
  item.colors.push_back({0.0f, 0.0f, 0.0f});
  draw_points.clear();
  draw_points.push_back(root);
  item.pvs.push_back(draw_points);

  // debugDrawFace(face, f, g, 1, 0, 0, item, true);

  int face_count = 0;

  if (avg_t_pos_roots_f.size() > 0 || avg_t_pos_roots_g.size() > 0) {
    PTR<Object<Scalar>> t;

    if (avg_t_pos_roots_f.size() > 0 && avg_t_pos_roots_g.size() > 0) {
      PTR<Object<Scalar>> avrf = avg_t_pos_roots_f[0];
      PTR<Object<Scalar>> avrg = avg_t_pos_roots_g[0];

      t = (LessThan(avrf, avrg) == -1 ? avrf : avrg);

    } else if (avg_t_pos_roots_f.size() > 0) {
      t = avg_t_pos_roots_f[0];
    } else {
      t = avg_t_pos_roots_g[0];
    }

    int st = Sign(t);

    PTR<Object<PV2>> point = new LinePoint(avg_t_line, t);

    PTR<Object<PV2>> mid_p = new MidVector(root, point);
    PTR<Line> sep_line = new Line(mid_p, rot_t);
    sep_line = sep_line->approximate(root, point);

    PTR<Face> other = face->split(sep_line);

    if (other)
      face_count++;
    else
      assert(other != nullptr);

    faces.push_back(other);
  }

  face->f = f;
  face->g = g;
  assert(face->alternates(face->f, face->g));

  if (avg_t_neg_roots_f.size() > 0 || avg_t_neg_roots_g.size() > 0) {
    PTR<Object<Scalar>> t;

    int fi = avg_t_neg_roots_f.size() - 1;
    int gi = avg_t_neg_roots_g.size() - 1;

    if (avg_t_neg_roots_f.size() > 0 && avg_t_neg_roots_g.size() > 0) {
      PTR<Object<Scalar>> avrf = avg_t_neg_roots_f[fi];
      PTR<Object<Scalar>> avrg = avg_t_neg_roots_g[gi];

      t = (LessThan(avrf, avrg) == -1 ? avrf : avrg);

    } else if (avg_t_neg_roots_f.size() > 0) {
      t = avg_t_neg_roots_f[fi];
    } else {
      t = avg_t_neg_roots_g[gi];
    }

    int st = Sign(t);

    PTR<Object<PV2>> point = new LinePoint(avg_t_line, t);

    PTR<Object<PV2>> mid_p = new MidVector(root, point);
    PTR<Line> sep_line = new Line(mid_p, rot_t_neg);
    sep_line = sep_line->approximate(root, point);

    PTR<Face> other = face->split(sep_line);

    if (other)
      face_count++;
    else
      assert(other != nullptr);

    faces.push_back(other);
  }

  face->f = f;
  face->g = g;
  assert(face->alternates(face->f, face->g));

  if (rot_t_pos_roots_f.size() > 0 || rot_t_pos_roots_g.size() > 0) {
    PTR<Object<Scalar>> t;

    if (rot_t_pos_roots_f.size() > 0 && rot_t_pos_roots_g.size() > 0) {
      PTR<Object<Scalar>> rot_t_f = rot_t_pos_roots_f[0];
      PTR<Object<Scalar>> rot_t_g = rot_t_pos_roots_g[0];

      t = (LessThan(rot_t_f, rot_t_g) == -1 ? rot_t_f : rot_t_g);

    } else if (rot_t_pos_roots_f.size() > 0) {
      t = rot_t_pos_roots_f[0];
    } else {
      t = rot_t_pos_roots_g[0];
    }

    int st = Sign(t);

    PTR<Object<PV2>> point = new LinePoint(rot_t_line, t);

    PTR<Object<PV2>> mid_p = new MidVector(root, point);
    PTR<Line> sep_line = new Line(mid_p, avg_t_neg);
    sep_line = sep_line->approximate(root, point);

    PTR<Face> other = face->split(sep_line);

    if (other)
      face_count++;
    else
      assert(other != nullptr);

    faces.push_back(other);
  }

  face->f = f;
  face->g = g;
  assert(face->alternates(face->f, face->g));

  if (rot_t_neg_roots_f.size() > 0 || rot_t_neg_roots_g.size() > 0) {
    PTR<Object<Scalar>> t;

    int fi = rot_t_neg_roots_f.size() - 1;
    int gi = rot_t_neg_roots_g.size() - 1;

    if (rot_t_neg_roots_f.size() > 0 && rot_t_neg_roots_g.size() > 0) {
      PTR<Object<Scalar>> rot_t_f = rot_t_neg_roots_f[fi];
      PTR<Object<Scalar>> rot_t_g = rot_t_neg_roots_g[gi];

      t = (LessThan(rot_t_f, rot_t_g) == -1 ? rot_t_f : rot_t_g);

    } else if (rot_t_neg_roots_f.size() > 0) {
      t = rot_t_neg_roots_f[fi];
    } else {
      t = rot_t_neg_roots_g[gi];
    }

    int st = Sign(t);

    PTR<Object<PV2>> point = new LinePoint(rot_t_line, t);

    PTR<Object<PV2>> mid_p = new MidVector(root, point);
    PTR<Line> sep_line = new Line(mid_p, avg_t);
    sep_line = sep_line->approximate(root, point);

    PTR<Face> other = face->split(sep_line);

    if (other)
      face_count++;
    else
      assert(other != nullptr);

    faces.push_back(other);
  }

  assert(face->alternates(face->f, face->g));

  // printf("splits for isolation: %d\n", face_count);

  face->intersection_point = root;

  return face_count;
}

vector<PTR<Object<PV2>>> Encasement::getRoots() { return theRoots; }

vector<PTR<Object<PV2>>> Encasement::_getRoots() {
  vector<PTR<Object<PV2>>> roots;

  // for each face that has an intersection
  int fsize = faces.size();

  vector<PTR<Face>> nFaces;

  for (int i = 0; i < fsize; i++) {
    PTR<Face> f = faces[i];

    if (!f->intersection_point) {
      if (!(f->f != nullptr && f->g != nullptr && f->contains(f->f) &&
            f->contains(f->g) && f->alternates(f->f, f->g)))
        continue;
    }

    f = f->clone();

    // until the face has four axis aligned edges
    while (!f->rect()) {
      f = shrinkRoot(f);
      if (!f) return vector<PTR<Object<PV2>>>();
    }

    nFaces.push_back(f);
  }

  // TODO
  // here is where we will check thet the bounding circle of each box doesn't
  // contain any points of any other box so that the resulting triangulation
  // contains the roots

  fsize = nFaces.size();

  for (int i = 0; i < fsize; i++) {
    PTR<Face> f = nFaces[i];

    // create a Root2 out of the resulting PV2 and add it to roots
    std::vector<PTR<Object<Scalar>>> bounds = f->bounds();

    PV2<Parameter> r(Parameter::interval(bounds[0]->getApprox(1.0).x.lb(),
                                         bounds[2]->getApprox(1.0).x.ub()),
                     Parameter::interval(bounds[1]->getApprox(1.0).x.lb(),
                                         bounds[3]->getApprox(1.0).x.ub()));

    PTR<Object<PV2>> root = new Root2(f->f, f->g, r);

    roots.push_back(root);
  }

  return roots;
}

// Pre: getRoots was called, ensuring no vert of a root lies in the circumcircle
//     of any other root box
vector<PTR<Triangle>> Encasement::getTriangulation(
    vector<PTR<Object<PV2>>> roots) {
  vector<PTR<Object<PV2>>> bb;
  bb.push_back(bl->p);
  bb.push_back(br->p);
  bb.push_back(tr->p);
  bb.push_back(tl->p);

  Triangle::rootIndex.erase(Triangle::rootIndex.begin(),
                            Triangle::rootIndex.end());

  // set the rootIndexes for all of the root points
  for (int i = 0; i < bb.size(); i++) {
    Triangle::rootIndex[bb[i]] = -1;
  }

  vector<PTR<Object<PV2>>> points;

  vector<PTR<Object<PV2>>> original_points;

  // shrink by finding closest point
  /*
    for(int i = 0; i < roots.size(); i++) {
      PV2 r = roots[i]->get();
      Parameter dist;
      bool set = false;

      for(int j = 0; j < roots.size(); j++) {
        if(j == i) continue;
        PV2 s = roots[j]->get();
        Parameter l = (r - s).length();
        if(!set || l < dist) {
          dist = l;
          set = true;
        }
      }

      for(int j = 0; j < sites.size(); j++) {
        PV2 s = sites[j]->get();
        Parameter l = (r - s).length();
        if(!set || l < dist) {
          dist = l;
          set = true;
        }
      }

      Parameter w = dist * cos(M_PI * 0.25);
      Parameter l = dist * sin(M_PI * 0.25);

      PV2 p_bl = PV2(r.x.lbP() - w, r.y.lbP() - l);
      PV2 p_br = PV2(r.x.ubP() + w, r.y.lbP() - l);
      PV2 p_tr = PV2(r.x.ubP() + w, r.y.ubP() + l);
      PV2 p_tl = PV2(r.x.lbP() - w, r.y.ubP() + l);

      PTR<Object<PV2>> r_bl = new InputPoint(PV2(p_bl.x.lbP(), p_bl.y.lbP()));
      PTR<Object<PV2>> r_br = new InputPoint(PV2(p_br.x.lbP(), p_br.y.lbP()));
      PTR<Object<PV2>> r_tr = new InputPoint(PV2(p_tr.x.lbP(), p_tr.y.lbP()));
      PTR<Object<PV2>> r_tl = new InputPoint(PV2(p_tl.x.lbP(), p_tl.y.lbP()));

      points.push_back(r_bl);
      points.push_back(r_br);
      points.push_back(r_tr);
      points.push_back(r_tl);

    }

  */

  for (int i = 0; i < roots.size(); i++) {
    PV2<Parameter> r = roots[i]->getApprox(1.0);
    PTR<Object<PV2>> r_bl =
        new InputPoint(PV2<Parameter>(r.x.lbP(), r.y.lbP()));
    PTR<Object<PV2>> r_br =
        new InputPoint(PV2<Parameter>(r.x.ubP(), r.y.lbP()));
    PTR<Object<PV2>> r_tr =
        new InputPoint(PV2<Parameter>(r.x.ubP(), r.y.ubP()));
    PTR<Object<PV2>> r_tl =
        new InputPoint(PV2<Parameter>(r.x.lbP(), r.y.ubP()));

    // Triangle::rootIndex[r_bl] = i;
    // Triangle::rootIndex[r_br] = i;
    // Triangle::rootIndex[r_tr] = i;
    // Triangle::rootIndex[r_tl] = i;

    original_points.push_back(r_bl);
    original_points.push_back(r_br);
    original_points.push_back(r_tr);
    original_points.push_back(r_tl);

    PTR<Object<PV2>> p0 = new InputPoint(
        PV2<Parameter>(bl->p->getApprox(1.0).x.lbP() -
                           (1 + randomNumber(0, 1)) * Parameter::delta,
                       bl->p->getApprox(1.0).y.lbP() -
                           (1 + randomNumber(0, 1)) * Parameter::delta));
    PTR<Object<PV2>> p1 = new InputPoint(
        PV2<Parameter>(br->p->getApprox(1.0).x.ubP() +
                           (1 + randomNumber(0, 1)) * Parameter::delta,
                       br->p->getApprox(1.0).y.lbP() -
                           (1 + randomNumber(0, 1)) * Parameter::delta));
    PTR<Object<PV2>> p2 = new InputPoint(
        PV2<Parameter>(tr->p->getApprox(1.0).x.ubP() +
                           (1 + randomNumber(0, 1)) * Parameter::delta,
                       tr->p->getApprox(1.0).y.ubP() +
                           (1 + randomNumber(0, 1)) * Parameter::delta));
    PTR<Object<PV2>> p3 = new InputPoint(
        PV2<Parameter>(tl->p->getApprox(1.0).x.lbP() -
                           (1 + randomNumber(0, 1)) * Parameter::delta,
                       tl->p->getApprox(1.0).y.ubP() +
                           (1 + randomNumber(0, 1)) * Parameter::delta));

    // Triangle::rootIndex[p0] = i;
    // Triangle::rootIndex[p1] = i;
    // Triangle::rootIndex[p2] = i;
    // Triangle::rootIndex[p3] = i;

    points.push_back(p0);
    points.push_back(p1);
    points.push_back(p2);
    points.push_back(p3);
  }

  bool sat = false;

  while (!sat) {
    sat = true;

    for (int i = 0; i < points.size(); i += 4) {
      bool shrink = false;

      // two right triangle circumcircles
      // PTR<Object<Circle>> c1 = new Circle90(points[i], points[i+1],
      // points[i+2]);
      PTR<Object<Circle>> c1 =
          new Circle3(points[i], points[i + 1], points[i + 2]);

      for (int j = 0; j < points.size(); j += 4) {
        if (j == i) continue;
        for (int k = 0; k < 4; k++) {
          if (CircleDist(c1, points[j + k]) < 0) {
            sat = false;
            shrink = true;
            break;
          }
        }
      }

      bool done = false;

      for (int j = 0; j < sites.size(); j++) {
        for (int k = 0; k < sites[j].size(); k++) {
          if (CircleDist(c1, sites[j][k]) < 0) {
            sat = false;
            shrink = true;
            done = true;
            break;
          }
          if (done) break;
        }
      }

      if (shrink) {
        PV2<Parameter> p, q;
        p = points[i]->getApprox(1.0);
        q = original_points[i]->getApprox(1.0);
        points[i] = new InputPoint(PV2<Parameter>::constant(
            0.5 * (p.x.lb() + q.x.lb()), 0.5 * (p.y.lb() + q.y.lb())));
        p = points[i + 1]->getApprox(1.0);
        q = original_points[i + 1]->getApprox(1.0);
        points[i + 1] = new InputPoint(PV2<Parameter>::constant(
            0.5 * (p.x.ub() + q.x.ub()), 0.5 * (p.y.lb() + q.y.lb())));
        p = points[i + 2]->getApprox(1.0);
        q = original_points[i + 2]->getApprox(1.0);
        points[i + 2] = new InputPoint(PV2<Parameter>::constant(
            0.5 * (p.x.ub() + q.x.ub()), 0.5 * (p.y.ub() + q.y.ub())));
        p = points[i + 3]->getApprox(1.0);
        q = original_points[i + 3]->getApprox(1.0);
        points[i + 3] = new InputPoint(PV2<Parameter>::constant(
            0.5 * (p.x.lb() + q.x.lb()), 0.5 * (p.y.ub() + q.y.ub())));
      }
    }
  }

  for (int i = 0; i < points.size(); i++) {
    Triangle::rootIndex[points[i]] = (i / 4);
  }

  int count = 0;

  for (int i = 0; i < sites.size(); i++) {
    for (int k = 0; k < sites[i].size(); k++) {
      // don't add the sites for now
      points.push_back(sites[i][k]);
      count++;
      Triangle::rootIndex[sites[i][k]] = -count - 2;
    }
  }

  return delaunay(bb, points);
}

void Encasement::makeTriangulationEncasement(vector<PTR<Triangle>> root) {
  faces.clear();

  vector<PTR<Triangle>> tri;
  queue<PTR<Triangle>> q;

  for (int i = 0; i < root.size(); i++) {
    q.push(root[i]);
  }

  map<PTR<Triangle>, int> tri_map;

  while (q.size() > 0) {
    PTR<Triangle> t = q.front();
    q.pop();

    if (t->chld[0] == 0) {
      if (tri_map.find(t) == tri_map.end()) {
        tri.push_back(t);
        tri_map[t] = 1;
      }
    } else {
      for (int i = 0; i < 3; i++)
        if (t->chld[i]) q.push(t->chld[i]);
    }
  }

  for (int i = 0; i < tri.size(); i++) {
    PTR<Triangle> t = tri[i];

    // Go through each vertex and create and edge opposite of it
    // as well as the corresponding twin edge for the
    for (int j = 0; j < 3; j++) {
      int k = (j + 1) % 3;
      int l = (j + 2) % 3;

      if (!t->edges[j]) {
        PTR<Vertex> vk = new Vertex(t->vert[k]);
        PTR<Vertex> vl = new Vertex(t->vert[l]);
        t->edges[j] = new Edge(vk, vl);
        if (t->next[j]) {
          int j2 = t->next[j]->find(t);
          t->next[j]->edges[j2] = t->edges[j]->twin;
        }
      }
    }

    // Link the edges
    for (int j = 0; j < 3; j++) {
      int k = (j + 1) % 3;
      t->edges[j]->next = t->edges[k];
      t->edges[k]->prev = t->edges[j];
    }
  }

  map<int, int> root_map;

  for (int i = 0; i < tri.size(); i++) {
    PTR<Triangle> t1 = tri[i];

    if (t1->rootTriangle()) {
      int rootI = Triangle::rootIndex[t1->vert[0]];
      if (root_map.find(rootI) == root_map.end()) {
        bool done = false;

        for (int j = 0; j < 3; j++) {
          if (t1->next[j]->rootTriangle() &&
              t1->sameRoot(t1->next[j]->vert[0])) {
            assert(!done);

            PTR<Triangle> t2 = t1->next[j];

            int j2 = t2->find(t1);

            PTR<Edge> e1 = t1->edges[j];
            PTR<Edge> e2 = t2->edges[j2];

            assert(e1->twin == e2);

            e1->swap(e2->prev);
            e2->swap(e1->prev);

            faces.push_back(new Face(t1->edges[(j + 1) % 3]));

            done = true;
          }
        }

        root_map[rootI] = 1;
      }
    } else {
      faces.push_back(new Face(tri[i]->edges[0]));
    }
  }

  intersect_edges(curves);

  stack<Item> sq;

  for (int i = 0; i < faces.size(); i++) {
    PTR<Face> face = faces[i];

    add_to_queue(face, sq);
  }

  while (!sq.empty()) {
    Item item = sq.top();
    sq.pop();

    PTR<Face> other = isolate_curves(item.face, item.f);

    add_to_queue(item.face, sq);
    add_to_queue(other, sq);
  }
}

void Encasement::add_to_queue(PTR<Face> face, stack<Item>& sq) {
  // take same curve intersections first
  for (set<PTR<Object<Poly2>>>::iterator it = face->curves.begin();
       it != face->curves.end(); it++) {
    PTR<Object<Poly2>> f = *it;
    if (face->intersectionCount(f) > 2) {
      face->num_splits = 0;
      face->f = nullptr;
      face->g = nullptr;
      sq.push({face, f, nullptr, false});
      return;
    }
  }
}

// PRE: calculate has been called and all faces are valid encasement faces
void Encasement::refine(PTR<Object<Poly2>> f, double desired_length) {
  queue<PTR<Face>> refinement_faces;

  for (auto& face : faces) refinement_faces.push(face);

  PTR<Object<Scalar>> desired_length_p =
      new InputScalar(Parameter::constant(desired_length));

  // If a face is empty or satisfies ||b -a|| < desired_length for its two curve
  // boundary points then it can be removed. Otherwise, split the cell with the
  // perpendicular bisector of a and b.
  while (refinement_faces.size() > 0) {
    PTR<Face> face = refinement_faces.front();
    refinement_faces.pop();

    // grab all boundary intersections
    vector<PTR<Object<PV2>>> boundary_intersections;
    PTR<Edge> e = face->e;

    do {
      if (e->intersection && e->poly == f)
        boundary_intersections.push_back(e->intersection_point());
    } while ((e = e->next) != face->e);

    // if this is an emtpy face, skip it
    if (boundary_intersections.size() == 0) continue;

    // Exepect a valid encasement, so at most 2 boundary intersections with f
    // assert(boundary_intersections.size() == 2);

    PTR<Object<PV2>> p = boundary_intersections[0];
    PTR<Object<PV2>> q = boundary_intersections[1];

    PTR<Object<PV2>> pq = new VectorAB(p, q);

    // If pq is already shorter than the desired length, continue
    if (boundary_intersections.size() == 2 &&
        LessThan(new VectorLength(pq), desired_length_p) < 0)
      continue;

    // Otherwise split the cell with their perpendicular bisector and recurse
    PTR<Object<PV2>> P = new MidVector(p, q);
    PTR<Object<PV2>> V = new NormalVector(new Rot90(pq, nullptr));

    PTR<Line> line_perf = new Line(P, V);
    PTR<Line> line = line_perf->approximate(p, q);

    PTR<Face> other = face->split(line);

    faces.push_back(other);
    other->num_splits = face->num_splits;

    refinement_faces.push(face);
    refinement_faces.push(other);
  }
}
