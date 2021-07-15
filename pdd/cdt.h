//  -*- Mode: c++ -*-
#ifndef CDT
#define CDT
#include "delaunay.h"
#include "heap.h"

Polyhedron * cdt (Polyhedron *a);

typedef pair<Point *, Tetrahedron *> PTPair;
typedef map<Point *, Tetrahedron *> PTMap;

void update (const Tetrahedrons &ts, unsigned int s, PTMap &ptmap);

class Segment {
 public:
  Point *a, *b;
  bool input, aflag, bflag;

  Segment () : a(0), b(0) {}
  Segment (Point *a, Point *b, bool input = false, bool aflag = false, bool bflag = false) 
    : a(a), b(b), input(input), aflag(aflag), bflag(bflag) {}
  bool operator== (const Segment &e) const { return a == e.a && b == e.b; }
  bool operator< (const Segment &e) const { return a < e.a || a == e.a && b < e.b; }
};

typedef vector<Segment> Segments;
typedef set<Segment> SegmentSet;
typedef pair<Segment, Facet> SFPair;
typedef map<Segment, Facet> SFMap;

SFMap facets (Polyhedron *a);

bool acute (Vertex *a);

class SmallAngle : public Primitive {
  Point *a, *b, *c;
  double k;
  
  DeclareSign {
    PV3<N> p = b->get<N>(), u = c->get<N>() - p, v = a->get<N>() - p;
    N uu = u.dot(u), vv = v.dot(v), uv = u.dot(v);
    return uv*uv - uu*vv*k*k;
  }
public:
  SmallAngle (Point *a, Point *b, Point *c, double k) : a(a), b(b), c(c), k(k) {}
};

void recoverSegments (Tetrahedrons &ts, FTMap &ftmap, PTMap &ptmap, SFMap &sfmap,
		      Points &steiner, Octree<Point *> *octree);

Segments missingSegments (const Tetrahedrons &ts, const SFMap &sfmap);

Segments segments (Tetrahedron *t);

void splitSegment (Tetrahedrons &ts, FTMap &ftmap, PTMap &ptmap, SFMap &sfmap,
		   Points &steiner, Octree<Point *> *octree, Segments &sm);

Point * splitSegment (Points &steiner, Octree<Point *> *octree, const Segment &s,
		      Segment &s1, Segment &s2);

Point * splitPoint (Octree<Point *> *octree, const Segment &s);

Point * referencePoint (Octree<Point *> *octree, const Segment &s);

void encroachedBBox (const Segment &s, double *bb);

class EncroachedSegment : public Primitive {
  Point *a, *b, *c;

  DeclareSign {
    PV3<N> p = a->get<N>(), q = b->get<N>(), r = c->get<N>();
    return (p - r).dot(r - q);
  }
public:
  EncroachedSegment (Point *a, Point *b, Point *c) : a(a), b(b), c(c) {}
};

double pointDistance (Point *a, Point *b);

class Distance : public Object<Scalar> {
  Point *a, *b;

  DeclareCalculate (Scalar) {
    PV3<N> u = b->get<N>() - a->get<N>();
    return u.dot(u).sqrt();
  }
public:
  Distance (Point *a, Point *b) : a(a), b(b) {}
};

class AffinePoint : public Point {
  PTR<Point> a, b;
  double s;

  DeclareCalculate (PV3) {
    PV3<N> p = a->get<N>(), u = b->get<N>() - p;
    return p + u*s;
  }
 public:
  AffinePoint (Point *a, Point *b, double s) : a(a), b(b), s(s) {}
};

void splitFacet (SFMap &sfmap, const Segment &s, Point *p);

Segments cavityBoundary (SFMap &sfmap, Point *p, const Facet &f);

void fillCavity (SFMap &sfmap, const Segment &s, Point *p, const Segments &cb,
		 unsigned int id);

void missingSegments (const Tetrahedrons &ts, unsigned int s, const SFMap &sfmap,
		      SegmentSet &sc, const Tetrahedrons &cav, Segments &sm);

void recoverFacets (Tetrahedrons &ts, FTMap &ftmap, PTMap &ptmap,
		    const SFMap &sfmap);

Facets missingFacets (const FTMap &ftmap, const SFMap &sfmap);

void recoverFacets (Tetrahedrons &ts, FTMap &ftmap, PTMap &ptmap,
		    const SFMap &sfmap, Facets &fm);

Facets facetPolygon (const FTMap &ftmap, const SFMap &sfmap, Facets &fm);

bool inputSegment (const SFMap &sfmap, Point *a, Point *b, Segment &s);

bool inputFacet (const SFMap &sfmap, const Facet &f, Facet &g);

bool cavityFacets (const FTMap &ftmap, const PTMap &ptmap, const Facets &fp,
		   Facets &fs);

Tetrahedron * containingTetrahedron (const FTMap &ftmap, const PTMap &ptmap,
				     const Facet &f);

bool contains (Tetrahedron *t, Point *p);

int intersects (const Facet &f, const Facet &g);

int intersects (TrianglePlane *p, TrianglePlane *q);

int intersectsAux (TrianglePlane *p, TrianglePlane *q);

bool overlaps (TrianglePlane *p, TrianglePlane *q);

int contains (TrianglePlane *p, Point *a);

bool intersectsEdges (TrianglePlane *p, Point *t, Point *h);

bool intersectsEE (Point *a, Point *b, Point *c, Point *d, int pc);

template<class N>
N orient3d (const PV3<N> &a, const PV3<N> &b, const PV3<N> &c, const PV3<N> &d)
{
  PV3<N> u = a - d, v = b - d, w = c - d;
  return u.tripleProduct(v, w);
}

template<class N>
class PV4 {
public:
  N x, y, z, w;
  
  PV4 () {}
  PV4 (const N &ix, const N &iy, const N &iz, const N &iw)
    : x(ix), y(iy), z(iz), w(iw) {}
  template<class M>
  PV4 (const PV4<M> &p) : x(p.x), y(p.y), z(p.z), w(p.w) {}
  PV4 & operator= (const PV4 &p) {
     x = p.x; y = p.y; z = p.z; w = p.w; return *this;
  }
  unsigned int size () const { return 4; }
  const N &operator[] (unsigned int i) const { return i == 0 ? x : i == 1 ? y : i == 2 ? z : w; }
  N &operator[] (unsigned int i) { return i == 0 ? x : i == 1 ? y : i == 2 ? z : w; }
  PV4 operator- (const PV4 &b) const { return PV4(x - b.x, y - b.y, z - b.z, w - b.w); }
};

template<class N>
N orient4d (const PV4<N> &a, const PV4<N> &b, const PV4<N> &c, const PV4<N> &d,
	    const PV4<N> &e)
{
  PV3<N> u(a.x - e.x, a.y - e.y, a.z - e.z), v(b.x - e.x, b.y - e.y, b.z - e.z),
    w(c.x - e.x, c.y - e.y, c.z - e.z), x(d.x - e.x, d.y - e.y, d.z - e.z);
  N uw = a.w - e.w, vw = b.w - e.w, ww = c.w - e.w, xw = d.w - e.w;
  return uw*v.tripleProduct(w, x) - vw*u.tripleProduct(w, x) +
    ww*u.tripleProduct(v, x) - xw*u.tripleProduct(v, w);
}

class FlipLifter : public Object<PV4> {
  Point *v, *a, *b, *c;

  DeclareCalculate (PV4) {
    PV3<N> vv = v->get<N>();
    if (!a)
      return PV4<N>(vv.x, vv.y, vv.z, N(0));
    N k = orient3d(a->get<N>(), b->get<N>(), c->get<N>(), vv);
    return PV4<N>(vv.x, vv.y, vv.z, k);
  }
public:
  FlipLifter (Point *v) : v(v), a(0), b(0), c(0) {} 
  FlipLifter (Point *v, Point *a, Point *b, Point *c) : v(v), a(a), b(b), c(c) {}
};

typedef vector<PTR<FlipLifter>> FlipLifters;

class FlipTime : public Object<Scalar> {
  PTR<FlipLifter> x, y, z, v, w;

  DeclareCalculate (Scalar) {
    PV4<N> xx = x->get<N>(), yy = y->get<N>(), zz = z->get<N>(), vv = v->get<N>(),
      ww = w->get<N>();
    PV3<N> x3(xx.x, xx.y, xx.z), y3(yy.x, yy.y, yy.z), z3(zz.x, zz.y, zz.z),
      v3(vv.x, vv.y, vv.z), w3(ww.x, ww.y, ww.z);
    N a = inSphere(x3, y3, z3, v3, w3), b = orient4d(xx, yy, zz, vv, ww);
    return b.sign() == 0 ? N(-1) : - a/b;
  }
public:
  FlipTime (FlipLifter *x, FlipLifter *y, FlipLifter *z, FlipLifter *v, FlipLifter *w)
    : x(x), y(y), z(z), v(v), w(w) {}
};

class ScalarOrder : public Primitive {
  Object<Scalar> *a, *b;

  DeclareSign { return b->get<N>().x - a->get<N>().x; }
public:
  ScalarOrder (Object<Scalar> *a, Object<Scalar> *b) : a(a), b(b) {}
};    

class FlipElt {
public:
  Facet f;
  Tetrahedron *t, *u;
  PTR<FlipTime> s;

  bool valid () const { return t->valid() && u->valid(); }

  bool operator== (const FlipElt &x) const {
    return t == x.t && u == x.u || t == x.u && u == x.t;
  }
  
  bool operator< (const FlipElt &x) const {
    return !(*this == x) && ScalarOrder(s, x.s) == 1;
  }

  FlipElt () : t(0), u(0), s(0) {}
  FlipElt (const Facet &f, Tetrahedron *t, Tetrahedron *u, FlipTime *s)
    : f(f), t(t), u(u), s(s) {}
};

typedef vector<FlipElt> FlipElts;

typedef Heap<FlipElt> FlipQ;

void flipInsertPolygon (Tetrahedrons &ts, FTMap &ftmap, PTMap &ptmap,
			const Facets &fp, const Facets &fs);

FlipQ flipInit (const FTMap &ftmap, const Facets &fs, TrianglePlane *p);

void certify (const FTMap &ftmap, const Facet &f, TrianglePlane *p,
	      Object<Scalar> *ts, FlipQ &q);

FlipElts nextFlips (FlipQ &q);

bool flip (Tetrahedrons &ts, FTMap &ftmap, PTMap &ptmap, const FlipElt &fe,
	   Tetrahedrons &cav);

void certify (const FTMap &ftmap, Tetrahedrons &cav, TrianglePlane *p,
	      Object<Scalar> *ts, FlipQ &q);

Tetrahedrons innerTetrahedrons (const Tetrahedrons &ts, const FTMap &ftmap,
				const FacetSet &fs);

Tetrahedron * getTetrahedron (const Facet &f, const FTMap &ftmap);

void setValid (const Tetrahedrons &ts, bool flag);

#endif
