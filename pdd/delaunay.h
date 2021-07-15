//  -*- Mode: c++ -*-
#ifndef DELAUNAY
#define DELAUNAY
#include "io.h"
#include <unordered_map>
#include <unordered_set>

class Facet {
public:
  Point *a, *b, *c;
  unsigned int id;
  
  Facet () : a(0), b(0), c(0), id(0u) {}

  Facet (Point *ai, Point *bi, Point *ci, unsigned int id = 0u) : id(id) {
    Point *pts[] = {ai, bi, ci};
    int i = pts[0] < pts[1] ? 0 : 1;
    if (pts[2] < pts[i])
      i = 2;
    a = pts[i];
    b = pts[(i+1)%3];
    c = pts[(i+2)%3];
  }

  bool ghost () const { return !a || !b || !c; }
  
  Point * vertex (unsigned int i) const {
    switch (i) {
    case 0: return a;
    case 1: return b;
    case 2: return c;
    }
    return 0;
  }

  Point * oppositePoint (Point *u, Point *v) const {
    for (unsigned int i = 0u; i < 3u; ++i)
      if (vertex(i) == u && vertex((i+1)%3) == v)
	return vertex((i+2)%3);
    return 0;
  }

  bool operator== (const Facet &x) const { return a == x.a && b == x.b && c == x.c; }
};

typedef vector<Facet> Facets;

class FacetHash {
public:
  size_t operator() (const Facet &f) const {
    return size_t(f.a) + size_t(f.b) + size_t(f.c);
  }
};

typedef unordered_set<Facet, FacetHash> FacetSet;

template<class N>
N inSphere (const PV3<N> &a, const PV3<N> &b, const PV3<N> &c, const PV3<N> &d,
	    const PV3<N> &e)
{
  PV3<N> u1 = a - e, u2 = b - e, u3 = c - e, u4 = d - e;
  return - u1.dot(u1)*u2.tripleProduct(u3, u4) + u2.dot(u2)*u1.tripleProduct(u3, u4)
    - u3.dot(u3)*u1.tripleProduct(u2, u4) + u4.dot(u4)*u1.tripleProduct(u2, u3);
}

class InSphere : public Primitive {
  Point *a, *b, *c, *d, *e;

  DeclareSign {
    PV3<N> aa = a->get<N>(), bb = b->get<N>(), cc = c->get<N>(), dd = d->get<N>(),
		 ee = e->get<N>();
    return inSphere(aa, bb, cc, dd, ee);
  }
public:
  InSphere (Point *a, Point *b, Point *c, Point *d, Point *e)
    : a(a), b(b), c(c), d(d), e(e) {}
};

class TriangleLinePoint : public Point {
  TrianglePlane *p;
  Point *t, *h;

  DeclareCalculate (PV3) {
    PlaneData<N> pd = p->get<N>();
    PV3<N> a = p->getA()->get<N>(), b = p->getB()->get<N>(), c = p->getC()->get<N>(),
      tt = t->get<N>(), u = h->get<N>() - tt;
    return tt - ((pd.k + pd.n.dot(tt))/pd.n.dot(u))*u;
  }
public:
  TriangleLinePoint (TrianglePlane *p, Point *t, Point *h) : p(p), t(t), h(h) {}
};

int intersectTriangleLine (TrianglePlane *p, Point *t, Point *h);

template<class N>
PV3<N> circleCenter (const PV3<N> a, const PV3<N> &b, const PV3<N> &c)
{
  PV3<N> ab = b - a, bc = c - b, n = ab.cross(bc), p = (a + b)/2,
    q = (b + c)/2, u = n.cross(ab), v = n.cross(bc);
  N k = (q - p).tripleProduct(v, n)/u.tripleProduct(v, n);
  return p + k*u;
}

class CircleCenter : public Point {
  PTR<Point> a, b, c;
  unsigned int id;
private:
  DeclareCalculate (PV3) {
    return circleCenter(a->get<N>(), b->get<N>(), c->get<N>());
  }
public:
  CircleCenter (Point *a, Point *b, Point *c, unsigned int id = 0u)
    : a(a), b(b), c(c), id(id) {}
  Point * getA () const { return a; }
  Point * getB () const { return b; }
  Point * getC () const { return c; }
  unsigned int getId () const { return id; }
};

bool inDiametricBall (Point *p, Point *a, Point *b, Point *c);

class Tetrahedron {
  Point * pts[4];
  bool flag;
  
public:
  Tetrahedron (Point *a, Point *b, Point *c, Point *d) : flag(true) {
    pts[0] = a; pts[1] = b; pts[2] = c; pts[3] = d;
  }
  
  bool valid () const { return flag; }
  void setValid (bool f) { flag = f; }
  Point * vertex (unsigned int i) const { return pts[i]; }
  bool isVertex (Point *p) const {
    return p == pts[0] || p == pts[1] || p == pts[2] || p == pts[3];
  }
  bool ghost () const { return !(pts[0] && pts[1] && pts[2] && pts[3]); }

  Facet facet (unsigned int i) const { 
    switch (i) {
    case 0: return Facet(pts[0], pts[2], pts[3]);
    case 1: return Facet(pts[2], pts[1], pts[3]);
    case 2: return Facet(pts[1], pts[0], pts[3]);
    case 3: return Facet(pts[0], pts[1], pts[2]);
    }
    cerr << "illegal facet index in Delaunay" << endl;
    exit(0);
  }

  Facet exitFacet (Point *t, Point *h) const {
    for (unsigned int i = 0u; i < 4u; ++i) {
      Facet f = facet(i);
      PTR<TrianglePlane> pl = new TrianglePlane(f.a, f.b, f.c);
      if (Side(pl, t) == -1 && Side(pl, h) == 1 &&
	  intersectTriangleLine(pl, t, h) > -1)
	return f;
    }
    cerr << "no exit facet in Delaunay" << endl;
    exit(0);
  }

  bool inCircumsphere (Point *p) const {
    if (ghost())
      for (unsigned int i = 0u; i < 4u; ++i) {
	Facet f = facet(i);
	if (!f.ghost()) {
	  PTR<TrianglePlane> pl = new TrianglePlane(f.a, f.b, f.c);
	  int s = Side(pl, p);
	  return s == -1 || s == 0 && inDiametricBall(p, f.a, f.b, f.c);
	}
      }
    return InSphere(pts[0], pts[1], pts[2], pts[3], p) == 1;
  }

  Point * oppositeVertex (const Facet &f) {
    for (unsigned int i = 0u; i < 4u; ++i)
      if (pts[i] != f.a && pts[i] != f.b && pts[i] != f.c)
	return pts[i];
    return 0;
  }    
};

typedef vector<Tetrahedron *> Tetrahedrons;
typedef pair<Facet, Tetrahedron *> FTPair;
typedef unordered_map<Facet, Tetrahedron *, FacetHash> FTMap;
typedef vector<Point*> PPoints;

Tetrahedrons delaunay (const PPoints &pts, double *bbox, FTMap &ftmap);

Tetrahedrons delaunay (PPoints &pts, FTMap &ftmap);

void delaunayInit (PPoints &pts, Tetrahedrons &ts, FTMap &ftmap);

bool firstTetrahedron (PPoints &pts);

void swap (PPoints &pts, unsigned int i, unsigned int j);

void addTetrahedron (Tetrahedrons &ts, FTMap &ftmap, Point *a, Point *b,
		     Point* c, Point *d);

void addTetrahedron (Tetrahedrons &ts, FTMap &ftmap, Tetrahedron *t);

void insertPoint (Point *p, Tetrahedrons &ts, FTMap &ftmap);

Tetrahedron * findTetrahedron (Point *p, Tetrahedron *t, const FTMap &ftmap);

Tetrahedron * getOppositeTetrahedron (const Facet &f, const FTMap &ftmap);

void insertPoint (Point *p, Tetrahedron *t, Tetrahedrons &ts, FTMap &ftmap);

Tetrahedrons cavity (Point *p, Tetrahedron *t, FTMap &ftmap);

void fillCavity (Tetrahedrons &ts, FTMap &ftmap, Point *p, const Tetrahedrons &cav);

Facets cavityBoundary (const Tetrahedrons &cav, const FTMap &ftmap);

void removeTetrahedrons (const Tetrahedrons &ts, FTMap &ftmap);

void removeTetrahedron (Tetrahedron *t, FTMap &ftmap);

void deleteTetrahedrons (const Tetrahedrons &ts);

Polyhedron * polyhedron (const Tetrahedrons &ts);

bool tetrahedralMesh (Polyhedron *a);

#endif
