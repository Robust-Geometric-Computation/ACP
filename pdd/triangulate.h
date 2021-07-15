#ifndef	TRIANGULATE
#define TRIANGULATE

#include "polyhedron.h"

class Edge2;

class Vertex2 {
 public:
  Vertex2 () : p(0) {}
  Vertex2 (Point *p, int coord) : p(p), coord(coord), edge(0), left(0) {}
  void addEdge (Edge2 *e);
  Edge2 * findEdge (Vertex2 *h) const;
  int yOrder (Vertex2 *b) const;
  int leftTurn (Vertex2 *b, Vertex2 *c) const;

  Point *p;
  int coord;
  Edge2 *edge, *left;
  bool flag;
};

typedef vector<Vertex2 *> Vertices2;

class Vertex2yOrder {
 public:
  bool operator() (Vertex2 *v, Vertex2 *w) const {
    return v != w && v->yOrder(w) == -1;
  }
};

class Edge2 {
 public:
  Edge2 (Vertex2 *tail, Edge2 *twin, bool flag)
    : tail(tail), helper(0), twin(twin), next(0), flag(flag) {}
  Vertex2 * head () const { return twin->tail; }
  bool increasingY ();
  bool clockwise (Edge2 *e);

  Vertex2 *tail, *helper;
  Edge2 *twin, *next;
  bool flag;
};

class Edge2xOrder {
 public:
  bool operator() (Edge2 *e, Edge2 *f) const {
    if (e == f)
      return false;
    if (e->head() == f->head())
      return f->tail->leftTurn(f->head(), e->tail) == 1;
    return e->head()->yOrder(f->head()) == 1 ?
      f->tail->leftTurn(f->head(), e->head()) == 1 :
      e->tail->leftTurn(e->head(), f->head()) == -1;
  }
};

typedef set<Edge2 *, Edge2xOrder> Edge2Set;

typedef map<Point *, Vertex2 *> PV2Map;

Edge2 * getEdge (Vertex2 *tail, Vertex2 *head, bool flag, bool tflag);

void triangulate (const vector<Points *> &reg, int coord, Triangles &tr);

void triangulateInit (const vector<Points *> &reg, int coord, Vertices2 &vertices);

Vertex2 * getVertex2 (Point *p, int coord, PV2Map &vmap);

void triangulate (Vertices2 &vertices, Triangles &tr);

void makeMonotone (Vertices2 &vertices);

void makeMonotone (Vertex2 *v, Edge2Set &sweep);

Edge2 * findLeft (Vertex2 *v, Edge2Set &sweep);

void triangulateMonotone (const Vertices2 &vertices, Triangles &tr);

void monotoneLoop (Edge2 *e, Vertices2 &u);

void triangulateLoop (const Vertices2 &u, Triangles &tr);

void triangulateSame (Vertex2 *v, Vertices2 &stack, Triangles &tr);

void triangulateOther (Vertex2 *v, Vertices2 &stack, Triangles &tr);

Triangle triangle (Vertex2 *a, Vertex2 *b, Vertex2 *c, bool flag);

void deleteVertices (Vertices2 &vertices);

#endif
