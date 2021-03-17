#ifndef ENCASEMENT_H
#define ENCASEMENT_H

/**
 *
 *@author Joseph Masterjohn
 *
 */

#include "acp/core/acp.h"
#include "acp/core/object.h"
#include "acp/dcel/dcel.h"
#include "acp/delaunay/delaunay.h"
#include "acp/encasement2d/critical_points.h"
#include "acp/encasement2d/predicates2d.h"
#include "acp/poly/poly.h"

#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <iterator>
#include <map>
#include <queue>
#include <stack>
#include <vector>

using namespace std;

typedef struct {
  PTR<Face> face;
  PTR<Object<Poly2>> f;
  PTR<Object<Poly2>> g;
  bool alternates;
  string from;
} Item;

class Encasement {
 public:
  bool print;

  PTR<Object<PV2>> bottom_left;
  PTR<Object<PV2>> top_right;

  pair<PTR<Face>, PTR<Face>> bad_faces;

  vector<vector<PTR<Object<PV2>>>> sites;

  vector<PTR<Face>> faces;

  vector<PTR<Object<Poly2>>> curves;

  vector<PTR<Object<PV2>>> theRoots;

  PTR<Vertex> bl;
  PTR<Vertex> tl;
  PTR<Vertex> br;
  PTR<Vertex> tr;

  int seed;

  int STEPS;
  int count;

  // Initialize the encasement with a bounding box
  Encasement(PTR<Object<PV2>> bottom_left, PTR<Object<PV2>> top_right, int seed,
             bool print = true)
      : bottom_left(bottom_left), top_right(top_right), print(print) {
    this->seed = seed;
    init();
  }

  void clean() {
    // faces = vector< PTR<Face> >();
    faces.clear();
    // sites = vector< vector< PTR<Object<PV2>> > >();
    // curves = vector< PTR<Object<Poly2>> >();
    // theRoots = vector<PTR<Object<PV2>>>();
    // count = 1;
  }

  void init() {
    STEPS = -1;

    // printf("\n\n");

    clean();

    srand(this->seed);
    srandom(this->seed);
    // srand(1);
    // srandom(1);

    // TODO
    // put a print in random to make sure the same sequence is coming out

    PTR<Object<PV2>> X =
        new InputPoint(PV2<Parameter>(Parameter(1), Parameter(0)));
    PTR<Object<PV2>> Y =
        new InputPoint(PV2<Parameter>(Parameter(0), Parameter(1)));

    PTR<Object<PV2>> VEC = new VectorAB(bottom_left, top_right);

    PTR<Object<Scalar>> W = new DotScalar(VEC, X);
    PTR<Object<Scalar>> H = new DotScalar(VEC, Y);

    PTR<Object<PV2>> XV = new ScaleVector(W, X);
    PTR<Object<PV2>> YV = new ScaleVector(H, Y);

    PTR<Object<PV2>> top_left = new SumPoint(bottom_left, YV);
    PTR<Object<PV2>> bottom_right = new SumPoint(bottom_left, XV);

    PTR<Object<PV2>> b = new InputPoint(
        PV2<Parameter>(Parameter::constant(0), Parameter::constant(-2.3)));
    PTR<Object<PV2>> r = new InputPoint(
        PV2<Parameter>(Parameter::constant(2.4), Parameter::constant(0)));
    PTR<Object<PV2>> t = new InputPoint(
        PV2<Parameter>(Parameter::constant(-0.3), Parameter::constant(2.4)));
    PTR<Object<PV2>> l = new InputPoint(
        PV2<Parameter>(Parameter::constant(-2.4), Parameter::constant(-0.2)));

    // PV2 top_left = PV2::input(bottom_left.getX().mid(),
    //				   bottom_left.getY().mid() + length.mid());
    // PV2 bottom_right = PV2::input(bottom_left.getX().mid() + width.mid(),
    //	       		       bottom_left.getY().mid());
    // PV2 top_right = PV2::input(bottom_left.getX().mid() + width.mid(),
    //		    		    bottom_left.getY().mid() + length.mid());

    bl = new Vertex(bottom_left);
    tl = new Vertex(top_left);
    br = new Vertex(bottom_right);
    tr = new Vertex(top_right);

    // PTR<Vertex> bv = new Vertex(b);
    // PTR<Vertex> rv = new Vertex(r);
    // PTR<Vertex> tv = new Vertex(t);
    // PTR<Vertex> lv = new Vertex(l);

    vector<PTR<Edge>> edges;

    edges.push_back(new Edge(bl, br));
    edges.push_back(new Edge(br, tr));
    edges.push_back(new Edge(tr, tl));
    edges.push_back(new Edge(tl, bl));

    // edges.push_back(new Edge(bl, bv));
    // edges.push_back(new Edge(bv, br));
    // edges.push_back(new Edge(br, rv));
    // edges.push_back(new Edge(rv, tr));
    // edges.push_back(new Edge(tr, tv));
    // edges.push_back(new Edge(tv, tl));
    // edges.push_back(new Edge(tl, lv));
    // edges.push_back(new Edge(lv, bl));

    // connect every other edge because of how they're stored in edges[]
    for (int i = 0; i < 4; i++) {
      edges[i]->next = edges[(i + 1) % 4];
      edges[(i + 1) % 4]->prev = edges[i];
      edges[(i + 1) % 4]->twin->next = edges[i]->twin;
      edges[i]->twin->prev = edges[(i + 1) % 4]->twin;
    }

    faces.push_back(new Face(edges[0]));
  }

  void calculate(vector<PTR<Object<Poly2>>> curves, bool verify = false);
  void _calculate(vector<PTR<Object<Poly2>>> curves, int i, int j);

  void calculate_single(PTR<Object<Poly2>> f);

  void refine(PTR<Object<Poly2>> f, double desired_length);

  void make_sites_verify(vector<PTR<Object<Poly2>>> curves);
  void make_sites(vector<PTR<Object<Poly2>>> curves);

  PTR<Encasement> intersect(PTR<Encasement> that);

  // Step 1: No cell contains a loop of a curve
  void eliminate_loops_ellipse(vector<PTR<Object<Poly2>>> curves);
  // void eliminate_loops(int i, int j);
  void eliminate_loops(int i);
  // Step 2: Intersect all cells with all curves
  void intersect_edges(vector<PTR<Object<Poly2>>> curves);
  // Step 3: All curves intersect a cell once
  PTR<Face> isolate_curves(PTR<Face> f, PTR<Object<Poly2>> curve);

  // Find an intersection using Newtons and then box it in creating
  // possibly 4 new faces. Return the number of faces created.
  int isolate_intersection(PTR<Face> face, PTR<Object<Poly2>> f,
                           PTR<Object<Poly2>> g);
  int isolate_intersection_new(PTR<Face> face, PTR<Object<Poly2>> f,
                               PTR<Object<Poly2>> g);
  int isolate_intersection_diamond(PTR<Face> face, PTR<Object<Poly2>> f,
                                   PTR<Object<Poly2>> g);

  int verify_intersection(PTR<Face> face, PTR<Object<Poly2>> f,
                          PTR<Object<Poly2>> g, PTR<Object<PV2>> r);

  int refine_narrow(PTR<Face> face, PTR<Object<Poly2>> f);
  void refine_subdiv(PTR<Face> face, PTR<Object<Poly2>> f);

  PTR<Face> separate(PTR<Face> face, PTR<Object<Poly2>> f,
                     PTR<Object<Poly2>> g);
  PTR<Face> separate_new(PTR<Face> face, PTR<Object<Poly2>> f,
                         PTR<Object<Poly2>> g);
  PTR<Face> separate(PTR<Face> face);

  PTR<Object<Scalar>> intersectionDistance(PTR<Edge> e1, PTR<Edge> e2);

  bool add_to_queues(PTR<Face> face, stack<Item>& nq, stack<Item>& q,
                     stack<Item>& sq);
  void add_to_queue(PTR<Face> face, stack<Item>& sq);

  PTR<Face> shrinkRoot(PTR<Face> f);

  vector<PTR<Object<PV2>>> getRoots();
  vector<PTR<Object<PV2>>> _getRoots();

  vector<PTR<Triangle>> getTriangulation(vector<PTR<Object<PV2>>>);

  void makeTriangulationEncasement(vector<PTR<Triangle>>);
};

#endif
