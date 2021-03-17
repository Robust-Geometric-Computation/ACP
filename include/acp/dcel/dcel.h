#ifndef DCEL_H
#define DCEL_H

#include "acp/core/acp.h"
#include "acp/core/object.h"
#include "acp/encasement2d/predicates2d.h"
#include "acp/linmath/mat.h"
#include "acp/linmath/pv_d.h"
#include "acp/poly/ppoly_d.h"

#include <set>
#include <unordered_set>
#include <vector>

using namespace acp;
using namespace std;

class Vertex;
class Edge;
class Face;
class BoundingBox;

class BoundingBox : public Object<PV2> {
  PTR<Face> face;
  DeclareCalculate(PV2);

 public:
  BoundingBox(PTR<Face> face) : face(face) {}
};

class Vertex : public RefCnt {
 public:
  PTR<Object<PV2>> p;  // coordiantes
  PTR<Edge> e;         // incident edge

  Vertex(PTR<Object<PV2>> p) : p(p) {}
};

//--Enum for the type of edge for destruction purposes
typedef enum { V2V, ELINE, TWIN } EDGETYPE;

class Edge : public RefCnt {
 public:
  //--Standard DCEL fields---------------------------------------------
  PTR<Vertex> tail;  // tail
  PTR<Vertex> head;  // head

  PTR<Edge> twin;  // twin

  PTR<Face> face;  // incident face

  PTR<Edge> next;  // next
  PTR<Edge> prev;  // prev

  //--Line parameterization--------------------------------------------
  PTR<Line> line;
  PTR<Object<Scalar>> t1;
  PTR<Object<Scalar>> t2;

  //--Reference to intersecting poly-----------------------------------
  PTR<Object<Poly2>> poly;
  PTR<Object<Scalar>> intersection;

  bool gradient_cross_point;

  EDGETYPE type;

  //--Axis aligned indicator for constructing root2s
  bool axis_aligned;

  //--Constructors-----------------------------------------------------
  // Connecting two vertices with no existing line between them
  Edge(PTR<Vertex> tail, PTR<Vertex> head) : tail(tail), head(head) {
    line = new LineAB(tail->p, head->p);
    t1 = new InputScalar(Parameter::constant(0));
    t2 = new InputScalar(Parameter::constant(1));

    init();

    type = EDGETYPE::V2V;
  }

  // Construct an edge from a line and two points on the line
  Edge(PTR<Line> line, PTR<Object<Scalar>> t1, PTR<Object<Scalar>> t2)
      : line(line), t1(t1), t2(t2) {
    tail = new Vertex(new LinePoint(line, t1));
    head = new Vertex(new LinePoint(line, t2));

    init();

    type = EDGETYPE::ELINE;
  }

  // Empty constructor for face splitting, must call init() after setting up
  // vertices and line
  Edge() {}

  //--Utility methods--------------------------------------------------
  PTR<Edge> vprev() { return twin->next; }

  void swap(PTR<Edge> that) {
    PTR<Edge> temp = this->next;
    this->next = that->next;
    that->next = temp;
    this->next->prev = this;
    that->next->prev = that;
  }

  PTR<Object<PV2>> intersection_point() {
    assert(poly != nullptr && intersection != nullptr);
    return new LinePoint(line, intersection);
  }

  bool operator==(const Edge& other) {
    return this->tail == other.tail && this->head == other.head;
  }

  bool operator!=(const Edge& other) { return !operator==(other); }

  PTR<Edge> clone();

  //--Methods implemented in .cc file----------------------------------
  void split(PTR<Object<Scalar>> t_prime);  // splits an edge at l(t_prime)

  // states whether this Edge intersects line l
  bool intersects(PTR<Line> l);

  // intersect this edge with curve f, split edge between intersections
  bool intersect(PTR<Object<Poly2>> f);

  // Pre:
  //  Edge is partially initialized with 2 verts, a line, and endpoints
  //
  void init() {
    twin = new Edge();

    twin->tail = head;
    twin->head = tail;

    twin->line = line;
    twin->t1 = t2;
    twin->t2 = t1;

    twin->twin = this;
    twin->type = EDGETYPE::TWIN;

    next = twin;
    prev = twin;
    twin->next = this;
    twin->prev = this;

    twin->face = face = nullptr;
    twin->poly = poly = nullptr;
    twin->intersection = intersection = nullptr;
    twin->gradient_cross_point = gradient_cross_point = false;

    axis_aligned = line->getAxis() != -1;
    twin->axis_aligned = axis_aligned;
  }
};

class Face : public RefCnt {
 public:
  PTR<Edge> e;  // incident edge
  bool satisfied;

  bool narrowed;

  int num_splits;

  std::set<PTR<Object<Poly2>>> curves;  // when satisfied, size <= 2
  PTR<Object<PV2>> intersection_point;  // when satisfied and size == 2
  PTR<Object<Poly2>> f;
  PTR<Object<Poly2>> g;

  PTR<Face> parent;
  int rank;

  Face()
      : e(nullptr),
        satisfied(false),
        narrowed(false),
        intersection_point(nullptr),
        f(nullptr),
        g(nullptr),
        num_splits(0),
        _bb(nullptr),
        parent(this),
        rank(0) {}

  Face(PTR<Edge> e)
      : e(e),
        satisfied(false),
        narrowed(false),
        intersection_point(nullptr),
        f(nullptr),
        g(nullptr),
        num_splits(0),
        _bb(nullptr),
        parent(this),
        rank(0) {
    PTR<Edge> current = e;
    do {
      current->face = this;
    } while ((current = current->next) != e);
  }

  void clear_parent() { parent = this; }

  PTR<Face> find() {
    if (this->parent == this)
      return this;
    else
      return (this->parent = this->parent->find());
  }

  void unionize(PTR<Face> that) {
    PTR<Face> p1 = find();
    PTR<Face> p2 = that->find();

    if (p1->rank < p2->rank) {
      p1->parent = p2;
      p2->rank = max(p1->rank + 1, p2->rank);
    } else {
      p2->parent = p1;
      p1->rank = max(p2->rank + 1, p1->rank);
    }
  }

  std::vector<PTR<Vertex>> get_verts() const {
    std::vector<PTR<Vertex>> v;
    PTR<Edge> curr = e;
    do {
      v.push_back(curr->head);
    } while ((curr = curr->next) != e);

    return v;
  }

  //--Methods to be implemented in .cc file----------------------------
  vector<PTR<Edge>> intersection_edges(PTR<Line> line);

  PTR<Face> split(PTR<Line> line);  // split a face by a line

  PTR<Face> clone();  // true to remove from list of faces

  bool intersect(PTR<Object<Poly2>> f);

  static PTR<Face> make(PV2<Parameter>& p);

  bool inside(PTR<Object<PV2>> p);

  // Does this face already know it intersects with f
  bool contains(PTR<Object<Poly2>> f);

  // Check if any edge of this face intersects f
  bool intersects(PTR<Object<Poly2>> f);

  // count the # of times f intersects this face
  int intersectionCount(PTR<Object<Poly2>> f);

  // true if f and g intersect this Face alternating
  bool alternates(PTR<Object<Poly2>> f, PTR<Object<Poly2>> g);

  bool verify(PTR<Object<Poly2>> f, PTR<Object<Poly2>> g, PTR<Object<PV2>> p);

  PTR<Object<PV2>> intersection(PTR<Object<Poly2>> F, PTR<Object<Poly2>> G);

  PTR<Object<PV2>> closest_approach(PTR<Object<Poly2>> F, PTR<Object<Poly2>> G,
                                    PTR<Object<PV2>> p);

  PV2D centroid();

  bool rect() const;

  std::vector<PTR<Object<Scalar>>> bounds() const;

  PTR<Object<PV2>> _bb;

  PTR<Object<PV2>> bb() {
    if (!_bb) _bb = new BoundingBox(this);
    return _bb;
  }
};

#endif
