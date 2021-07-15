//  -*- Mode: c++ -*-
#ifndef CDT
#define CDT
#include "delaunay.h"
#include "heap.h"

Polyhedron * mesh (Polyhedron *a, double r, double cs);

typedef set<Point *> PPointSet;
typedef pair<Point *, Tetrahedron *> PTPair;
typedef unordered_map<Point *, Tetrahedron *> PTMap;

class AffinePoint : public Point {
public:
  Point *a, *b;
  double s;
  unsigned int id;
private:
  DeclareCalculate (PV3) {
    PV3<N> p = a->get<N>(), u = b->get<N>() - p;
    return p + u*s;
  }
 public:
  AffinePoint (Point *a, Point *b, double s, unsigned int id = 0u)
    : a(a), b(b), s(s), id(id) {}
  unsigned int getId () const { return id; }
};

double distance (Point *a, Point *b);

class Distance : public Object<Scalar> {
  Point *a, *b;

  DeclareCalculate (Scalar) {
    PV3<N> u = b->get<N>() - a->get<N>();
    return u.dot(u).sqrt();
  }
public:
  Distance (Point *a, Point *b) : a(a), b(b) {}
};

class Segment {
 public:
  Point *a, *b;
  unsigned int id;
  bool aflag, bflag;

  Segment () : a(0), b(0) {}
  Segment (Point *a, Point *b, unsigned int id = 0u, bool aflag = false, bool bflag = false) 
    : a(a), b(b), id(id), aflag(aflag), bflag(bflag) {}
  bool operator== (const Segment &e) const { return a == e.a && b == e.b; }
  bool operator< (const Segment &e) const { return a < e.a || a == e.a && b < e.b; }
};

class SegmentCC : public Segment, public RefCnt {
public:
  double bbox[6];
  
  SegmentCC (const Segment &s) : Segment(s) {
    PTR<Point> p = new AffinePoint(a, b, 0.5);
    double r = distance(p, a);
    PV3<double> q = p->getApproxMid();
    for (unsigned int i = 0u; i < 3u; ++i) {
      bbox[2*i] = q[i] - r;
      bbox[2*i+1] = q[i] + r;
    }
  }

  void getBBox (double *bb) const { copyBBox(bbox, bb); }
};

class SegmentHash {
public:
  size_t operator() (const Segment &s) const {
    return size_t(s.a) + size_t(s.b);
  }
};

typedef vector<Segment> Segments;
typedef set<Segment> SegmentSet;
typedef vector<PTR<SegmentCC>> SegmentCCs;
typedef pair<Segment, Facet> SFPair;
typedef vector<SFPair> SFPairs;
typedef unordered_map<Segment, Facet, SegmentHash> SFMap;
typedef pair<Point *, double> PRPair;
typedef unordered_map<Point *, double> PRMap;
typedef pair<Facet, Tetrahedron *> FTPair;
typedef vector<FTPair> FTPairs;

class Skinny {
 public:
  Tetrahedron *t;
  PTR<Point> c, p;
  double r;

  Skinny () : t(0), c(0), p(0), r(0) {}
  Skinny (Tetrahedron *t, Point *c, Point *p, double r)
    : t(t), c(c), p(p), r(r) {}
  bool operator< (const Skinny &x) const { return x.r < r; }
};

typedef Heap<Skinny> SkinnyQ;

class IIPair {
public:
  unsigned int m, n;

  IIPair (unsigned int i, unsigned int j) : m(min(i, j)), n(max(i, j)) {}
  bool operator< (const IIPair &x) const { return m < x.m || m == x.m && n < x.n; }
};

class FacetPlane : public TrianglePlane {
public:
  unsigned int id;
  double bbox[6];

  FacetPlane (Point *a, Point *b, Point *c, unsigned int id)
    : TrianglePlane(a, b, c), id(id) {
    a->getBBox(bbox);
    double bb[6];
    b->getBBox(bb);
    mergeBBox(bb, bbox);
    c->getBBox(bb);
    mergeBBox(bb, bbox);
  }

  void getBBox (double *b) const { copyBBox(bbox, b); }
};

typedef vector<PTR<FacetPlane>> FacetPlanes;

class Split {
public:
  Segment s;
  Point *par;
  bool ent;

  Split (const Segment &s, Point *par, bool ent) : s(s), par(par), ent(ent) {}
};

typedef vector<Split> Splits;

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

class Mesh {
  Polyhedron *a;
  double r, cs, rm;
  SFMap sfmap;
  Octree<PTR<FacetPlane>> *foctree;
  PRMap prmap;
  set<IIPair> ent;
  Octree<PTR<Point>> *poctree;
  Points ppars;
  Tetrahedrons ts;
  FTMap ftmap;
  PTMap ptmap;
  Octree<PTR<SegmentCC>> *soctree;
  set<Tetrahedron *> fenced;

  void init ();
  void initSFMap ();
  unsigned int segmentId (Point *a, Point *b) const;
  bool inputSegment (Point *a, Point *b, Segment &s) const;
  double getRR (Point *p) const { return prmap.find(p)->second; }
  void initFoctree ();
  bool entwined (unsigned int i, unsigned int j) const {
    return ent.find(IIPair(i, j)) != ent.end();
  }
  void initEntwined ();
  unsigned int facetId (Point *a, Point *b, Point *c) const;
  bool inputFacet (const Facet &f, Facet &g) const;
  bool entwinedFE (Face *f, HEdge *e) const;
  void initPRMap ();
  PPoints visibleVertices (Point *p, double dmax) const;
  bool visible (Point *p, Point *q) const;
  void updatePRMap (AffinePoint *p, bool cd = false, Point *par = 0,
		    bool entw = false);
  PPoints neighbors (AffinePoint *p) const;
  PPoints neighborsCD (Point *p) const;
  Tetrahedrons vertexTetrahedrons (Point *p) const;
  double relaxedDistance (AffinePoint *p, Point *q) const;
  void updatePRMap (CircleCenter *p);
  double calculateRR (CircleCenter *p) const;
  bool eligible (Point *p, const Facet &f) const;
  void updatePTMap (unsigned int s);
  SFPairs encroachedSegments () const;
  void encroachedSegments (unsigned int st, SegmentSet &cand, SFPairs &enc) const;
  void splitSegment (SFPairs &enc);
  AffinePoint * splitSegment (const Segment &s, Segment &s1, Segment &s2);
  AffinePoint * splitPoint (const Segment &s) const;
  void splitSegmentFacet (const Segment &s, Point *p, Facets *fm = 0);
  Segments cavityBoundary (Point *p, const Facet &f, Facets &fs);
  void eraseFacets (const Facets &fs);
  void fillCavity (const Segment &s, Point *p, const Segments &cb,
		   unsigned int id);
  void boundaryFacets (const Segments &cb, Facets &fm) const;
  SegmentSet inputSegments (const Tetrahedrons &cav) const;
  void initSoctree ();
  void refine ();
  Facets missingFacets () const;
  FTPairs encroachedFacets () const;
  void encroachedFacets (unsigned int s, FTPairs &enc) const;
  void encroachedFacet (Point *p, const Facet &f, FTPairs &enc) const;
  Facet orthogonalProjection (Point *p, const Facet &f) const;  
  SkinnyQ skinnyTetrahedrons () const;
  Tetrahedrons innerTetrahedrons () const;
  void skinnyTetrahedrons (const Tetrahedrons &tc, unsigned int s,
			   SkinnyQ &sk) const;
  bool checkSkinny (Tetrahedron *t, SkinnyQ &sk) const;
  bool sliver (Tetrahedron *t) const;
  void skinnyTetrahedrons (unsigned int st, SkinnyQ &sk) const;
  void splitSegment (Splits &sp, Facets &fm, FTPairs &fenc);
  Tetrahedrons segmentTetrahedrons (const Segment &s) const;
  Tetrahedron * intersectingTetrahedron (const Segment &s) const;
  Tetrahedron * findValid () const;
  Tetrahedron * containingTetrahedron (Point *p, Tetrahedron *t) const;
  void encroachedSegments (Splits &sp, const Segment &s1, const Segment &s2,
			   Point *p, Point *par, bool ent) const;
  Segments encroachedSegments (Point *p) const;
  void encroachedSegment (Splits &sp, const Segment &s) const;
  void recoverFacet (Facets &fm, FTPairs &fenc);
  Facets facetPolygon (Facets &fm) const;
  bool cavityFacets (const Facets &fp, Facets &fs) const;
  Tetrahedron * containingTetrahedron (const Facet &f) const;
  void flipInsertPolygon (const Facets &fp, const Facets &fs);
  FlipQ flipInit (const Facets &fs, TrianglePlane *p) const;
  void certify (const Facet &f, TrianglePlane *p, Object<Scalar> *ts, FlipQ &q) const;
  FlipElts nextFlips (FlipQ &q) const;
  bool flip (const FlipElt &fe, Tetrahedrons &cav);
  void certify (Tetrahedrons &cav, TrianglePlane *p, Object<Scalar> *ts, FlipQ &q) const;
  void splitFacet (Splits &sp, Facets &fm, FTPairs &fenc);
  void splitFacet (Splits &sp, Facets &fm, FTPairs &fenc, const Facet &f,
		   CircleCenter *p, Tetrahedron *t);
  Tetrahedrons facetTetrahedrons (const Facet &f, Point *p, Tetrahedron *t) const;
  void splitFacet (const Facet &f, Point *p, Facets &fm);
  void removeType3Neighbors (Point *p, Facets &fm);
  void removeVertex (Point *p, Facets &fm);
  void fillCavity (Facets &fm, PPoints &pts, const Facets &fs);
  bool checkCavity (Facets &fm, PPoints &pts, FacetSet &fs,
		    const FTMap &ftmapc);
  void splitTetrahedron (Splits &sp, Facets &fm, FTPairs &fenc, SkinnyQ &sk);
  Point * occluder (Tetrahedron *t, Point *c, Facet &f) const;
  void splitTetrahedron (Splits &sp, FTPairs &fenc, SkinnyQ &sk);
  bool visible (Point *p, const Segment &s) const;
  Facets encroachedFacets (Tetrahedron *t, Point *p);
  bool partiallyVisible (Point *p, const Facet &f) const;
  bool entwinedOccluders (Point *c, const Segment &s) const;
  void splitTetrahedron (const Facet &f, Point *cp, FTPairs &fenc, SkinnyQ &sk);
  void insertVertex (Point *p, Tetrahedrons &tp, Splits &sp, Facets &fm);
  Tetrahedrons cavityCD (Point *p, const Tetrahedrons &pts) const;
  bool starShaped (Point *p, const FacetSet &fs, Facet &fb) const;
  void missingSF (const Tetrahedrons &cav, const FacetSet &fs, Point *par,
		  Splits &sp, Facets &fm) const;
  void insertVertexT (SkinnyQ &sk);
  void updatePRMap (Tetrahedron *t, Point *c, Point *p);
  bool checkCD (unsigned int s) const;
public:
  Mesh (Polyhedron *a, double r, double cs, double rm)
    : a(a), r(r), cs(cs), rm(rm), foctree(0), poctree(0), soctree(0) {}
  ~Mesh () {
    deleteTetrahedrons(ts); delete foctree; delete poctree; delete soctree;
  }
  Tetrahedrons mesh ();
};

class Acute : public Primitive {
  Point *a, *b, *c;
  
  DeclareSign {
    PV3<N> p = b->get<N>(), u = c->get<N>() - p, v = a->get<N>() - p;
    return u.dot(v);
  }
public:
  Acute (Point *a, Point *b, Point *c) : a(a), b(b), c(c) {}
};

class AcuteEdge : public Primitive {
  Edge *e;

  DeclareSign {
    PV3<N> m = e->getHEdge(0)->getF()->getP()->get<N>().n,
      n = e->getHEdge(1)->getF()->getP()->get<N>().n;
    return - m.dot(n);
  }
public:
  AcuteEdge (Edge *e) : e(e) {}
};

class EncroachedSegment : public Primitive {
  Point *a, *b, *c;

  DeclareSign {
    PV3<N> p = a->get<N>(), q = b->get<N>(), r = c->get<N>();
    return (p - r).dot(r - q);
  }
public:
  EncroachedSegment (Point *a, Point *b, Point *c) : a(a), b(b), c(c) {}
};

class InDiametricEllipse : public Primitive {
  Point *p, *a, *b, *c;
  
  DeclareSign {
    PV3<N> ap = a->get<N>(), bp = b->get<N>(), cp = c->get<N>(),
      o = circleCenter(ap, bp, cp), n = (cp - bp).cross(ap - bp), ao = ap - o,
      po = p->get<N>() - o, z = po.dot(n)/n.dot(n)*n, xy = po - z;
    return ao.dot(ao) - xy.dot(xy) - 3*z.dot(z);
  }
public:
  InDiametricEllipse (Point *p, Point *a, Point *b, Point *c) : p(p), a(a), b(b), c(c) {}
};

class PointNearPlane : public Primitive {
  Point *a;
  Plane *p;
  double d;

  DeclareSign {
    PlaneData<N> pd = p->get<N>();
    PV3<N> n = pd.n;
    N k = n.dot(a->get<N>()) + pd.k;
    return d*d*n.dot(n) - k*k;
  }
public:
  PointNearPlane (Point *a, Plane *p, double d) : a(a), p(p), d(d) {}
};

void setValid (const Tetrahedrons &ts, bool flag);

Segments segments (Tetrahedron *t);

bool contains (Tetrahedron *t, Point *p);

bool intersects (Tetrahedron *t, const Segment &s);

int intersects (TrianglePlane *p, Point *t, Point *h, bool strict = true);

bool overlaps (TrianglePlane *p, Point *t, Point *h);

int intersects (const Facet &f, const Facet &g);

int intersects (TrianglePlane *p, TrianglePlane *q);

int intersectsAux (TrianglePlane *p, TrianglePlane *q);

bool overlaps (TrianglePlane *p, TrianglePlane *q);

int contains (TrianglePlane *p, Point *a);

bool intersectsEdges (TrianglePlane *p, Point *t, Point *h);

bool intersectsEE (Point *a, Point *b, Point *c, Point *d, int pc);

PPoints points (const Facets &fs);

Tetrahedron * getTetrahedron (const Facet &f, const FTMap &ftmap);

Tetrahedrons innerTetrahedrons (const FTMap &ftmap, const FacetSet &fs);

class OrthogonalSide : public Primitive {
  Point *p, *a, *b;
  TrianglePlane *pl;

  DeclareSign {
    PV3<N>  pp = p->get<N>(), aa = a->get<N>(), bb = b->get<N>(),
      n = pl->get<N>().n.cross(bb - aa);
    return (pp - aa).dot(n);
  }
public:
  OrthogonalSide (Point *p, Point *a, Point *b, TrianglePlane *pl)
    : p(p), a(a), b(b), pl(pl) {}
};

class SphereCenter : public Point {
public:
  Point *a, *b, *c, *d;
private:
  DeclareCalculate (PV3) {
    PV3<N> aa = a->get<N>(), bb = b->get<N>(), cc = c->get<N>(), dd = d->get<N>(),
      n = (cc - bb).cross(aa - bb), o = circleCenter(aa, bb, cc), oa = o - aa,
      od = o - dd, ad = aa - dd;
    N k = (oa.dot(oa) - od.dot(od))/(2*n.dot(ad));
    return o + k*n;
  }
public:
  SphereCenter (Point *a, Point *b, Point *c, Point *d) : a(a), b(b), c(c), d(d) {}
  SphereCenter (double x, double y, double z) : Point(x, y, z, false) {}
};

class RayFacetPoint : public Point {
  PTR<Point> t, h;
  Facet f;

  DeclareCalculate (PV3) {
    PV3<N> a = t->get<N>(), u = h->get<N>() - a, fa = f.a->get<N>(),
      fb = f.b->get<N>(), fc = f.c->get<N>(), n = (fc - fb).cross(fa - fb);
    N k = n.dot(fa - a)/n.dot(u);
    return a + k*u;
  }
public:
  RayFacetPoint (Point *t, Point *h, const Facet &f) : t(t), h(h), f(f) {}
};

class ProjectionPoint : public Point {
  TrianglePlane *pl;
  Point *p;

  DeclareCalculate (PV3) {
    PV3<N> pp = p->get<N>(), pa = pp - pl->getA()->get<N>(), n = pl->get<N>().n;
    return pp - (pa.dot(n)/n.dot(n))*n;
  }
public:
  ProjectionPoint (TrianglePlane *pl, Point *p) : pl(pl), p(p) {}
};

class SmallAngle : public Primitive {
  Point *a, *b, *c, *d;
  double cs;

  DeclareSign {
    PV3<N> p = a->get<N>(), q = b->get<N>(),
      m = (c->get<N>() - p).cross(q - p), n = (d->get<N>() - q).cross(p - q);
    N mn = m.dot(n), mm = m.dot(m), nn = n.dot(n);
    return mn*mn - (cs*cs)*mm*nn;
  }
public:
  SmallAngle (Point *a, Point *b, Point *c, Point *d, double cs)
    : a(a), b(b), c(c), d(d), cs(cs) {}
};

bool containsProjection (Point *a, Point *b, Point *c, Point *q);

Point * randomPointInSphere (Point *o, double r);

#endif
