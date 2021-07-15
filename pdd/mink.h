#ifndef MINK
#define MINK

#include "polyhedron.h"

Polyhedron * minkowskiSum (Polyhedron *a, Polyhedron *b, bool check = false);

Shells minkowskiShells (Polyhedron *a, Polyhedron *b, const Shells &sh);

Shells minkowskiShellsAux (Polyhedron *a, Polyhedron *b, const Shells &sh);

void minkowskiShellsT (void *ptr);

class MSData {
 public:
  unsigned int i, is, ie;
  Polyhedron *a;
  Octree<Face *> *octree;
  const Shells *sh;
  Shells res;
};

bool minkowskiShell (Polyhedron *a, Octree<Face *> *octree, Shell *s);

class UnitVector : public Point {
  PTR<Point> a;

  DeclareCalculate (PV3) { return a->get<N>().unit(); }
 public:
  UnitVector (Point *a) : a(a) {}
};

class GreatCircleN : public Point {
  Point *t, *h;

  DeclareCalculate (PV3) { return t->get<N>().cross(h->get<N>()); }
 public:
  GreatCircleN (Point *t, Point *h) : t(t), h(h) {}
};

class GreatCircleMinX : public Point {
  Point *n;
  
  DeclareCalculate (PV3) {
    PV3<N> p = n->get<N>();
    return PV3<N>(- p.y*p.y - p.z*p.z, p.y*p.x, p.z*p.x);
  }
 public:
  GreatCircleMinX (Point *n) : n(n) {}
};

class GreatCircleMinY : public Point {
  Point *n;
  
  DeclareCalculate (PV3) {
    PV3<N> p = n->get<N>();
    return PV3<N>(p.x*p.y, - p.x*p.x - p.z*p.z, p.z*p.y);
  }
 public:
  GreatCircleMinY (Point *n) : n(n) {}
};

class GreatCircleMinZ : public Point {
  Point *n;
  
  DeclareCalculate (PV3) {
    PV3<N> p = n->get<N>();
    return PV3<N>(p.x*p.z, p.y*p.z, - p.x*p.x - p.y*p.y);
  }
 public:
  GreatCircleMinZ (Point *n) : n(n) {}
};

class EENormal : public Point {
  HEdge *e, *f;

  DeclareCalculate (PV3) { return e->getU<N>().cross(f->getU<N>()); }
 public:
  EENormal (HEdge *e, HEdge *f) : e(e), f(f) {}
};

class FNormal : public Point {
  Face *f;
  
  DeclareCalculate (PV3) { return f->getP()->get<N>().n; }
 public:
  FNormal (Face *f) : f(f) {}
};

class TripleProduct :public Primitive {
  Point *a, *b, *c;

  DeclareSign {
    return a->get<N>().tripleProduct(b->get<N>(), c->get<N>());
  }
 public:
  TripleProduct (Point *a, Point *b, Point *c) : a(a), b(b), c(c) {}
};

void sphereBBox (Edge *e, double *bbox);

void sphereBBox (Point *t, Point *h, double *bbox);

void sphereBBox (const HEdges &ed, double *bbox);

int coordinate (Point *a, int c);

class BSPElt {
 public:
  BSPElt (Vertex *v, const HEdges &ed, bool convex = true);
  BSPElt (Edge *e, bool convex = true);
  BSPElt (Face *f);
  BSPElt (const BSPElt &x);
  void setPbb ();
  bool compatible (const BSPElt &e) const { return (l & e.l) == 0u; }
  int side (Point *r) const;

  union {
    Vertex *v;
    Edge *e;
    Face *f;
  } d;
  HEdges ed;
  unsigned int l;
  double bbox[6];
  PV3<Parameter> pbb;
  bool convex;
};

typedef vector<BSPElt> BSPElts;

void BSPTree (BSPElts &aelts, BSPElts &belts, BSPElts &ea, BSPElts &eb,
	      int nmax = 40, int dmax = 20, unsigned int c = 0u);

void BSPPartition (BSPElts &elts, Point *r, unsigned int c, BSPElts &elts1,
		   BSPElts &elts2);

void BSPLeaf (const BSPElts &aelts, const BSPElts &belts, BSPElts &ea, 
	      BSPElts &eb);

class MinkHullFace {
 public:
  MinkHullFace (HEdge *e, MinkHullFace *prev, MinkHullFace *next)
    : e(e), prev(prev), next(next) {}
  bool inCset (HEdge *f) const
  { return find(cset.begin(), cset.end(), f) != cset.end(); }
  void updateCset (MinkHullFace *h, HEdge *f);
  bool conflict (HEdge *f) const;
  void cone (HEdges &hedges) const;

  HEdges cset;
  HEdge *e;
  MinkHullFace *prev, *next;
};

class DegenerateConflict1 : public Primitive {
  Point *a, *b, *c, *d;

  DeclareSign {
    PV3<N> u = b->get<N>() - a->get<N>(), v = c->get<N>() - a->get<N>(),
      w = d->get<N>() - a->get<N>(), x = u.cross(v);
    return x.tripleProduct(v, w);
  }
 public:
  DegenerateConflict1 (Point *a, Point *b, Point *c, Point *d)
    : a(a), b(b), c(c), d(d) {}
};

class DegenerateConflict2 : public Primitive {
  Point *a, *b, *c, *d;

  DeclareSign {
    PV3<N> u = b->get<N>() - a->get<N>(), v = c->get<N>() - a->get<N>(),
      w = d->get<N>() - a->get<N>(), x = u.cross(v);
    return x.tripleProduct(w, u);
  }
 public:
  DegenerateConflict2 (Point *a, Point *b, Point *c, Point *d)
    : a(a), b(b), c(c), d(d) {}
};

bool convexCone (Vertex *v, HEdges &hedges);

MinkHullFace * initHull (HEdges &hedges);

int circulationEEE (HEdge *e, HEdge *f, HEdge *g);

MinkHullFace * updateHull (MinkHullFace *hull, HEdge *e, bool &flag);

MinkHullFace * updateHullAux (MinkHullFace *fs, MinkHullFace *fe, HEdge *e);

void deleteHull (MinkHullFace *hull);

bool convexOrder (const HEdges &hedges);

Polyhedron * convolution (Polyhedron *a, Polyhedron *b, bool check);

typedef map<pair<Point *, Point *>, Point *> PPPMap;

class SF {
 public:
  PTR<Point> a, b, c, d;

  SF (Point *a, Point *b, Point *c, Point *d = 0) : a(a), b(b), c(c), d(d) {}
};

typedef vector<SF> SFs;

void sumVF (Polyhedron *a, Polyhedron *b, bool avflag, PPPMap &pmap, SFs &sfs);

void convexVertices (Polyhedron *a, BSPElts &elts);

bool compatibleVF (HEdges &ed, Face *f);

void sumVF (Vertex *v, Face *f, bool avflag, PPPMap &pmap, SFs &sfs);

class InnerProductEF : public Primitive {
  HEdge *e;
  Face *f;

  DeclareSign {
    return e->getU<N>().dot(f->getP()->get<N>().n);
  }
 public:
  InnerProductEF (HEdge *e, Face *f) : e(e), f(f) {}
};

Point * sumPP (Point *a, Point *b, bool aflag, PPPMap &pmap);

void sumEE (Polyhedron *a, Polyhedron *b, PPPMap &pmap, SFs &sfs);

void convexEdges (Polyhedron *a, BSPElts &elts);

int convexEdge (Edge *e);

class ConvexEdge : public Primitive {
  HEdge *e1, *e2;

  DeclareSign {
    return e1->getU<N>().tripleProduct(e1->getF()->getP()->get<N>().n,
				       e2->getF()->getP()->get<N>().n);
  }
 public:
  ConvexEdge (HEdge *e1, HEdge *e2) : e1(e1), e2(e2) {}
};

bool compatibleEE (Edge *e, Edge *f, bool &aflag);

void sumEE (Edge *e, Edge *f, bool aflag, PPPMap &pmap, SFs &sfs);

Polyhedron * convolution (const SFs &sfs, bool check);

PVMap pvMap (const SFs &sfs, Polyhedron *c);

Vertices vertices (const SF &s, Polyhedron *a, PVMap &pvmap);

#endif
