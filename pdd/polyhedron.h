#ifndef POLY
#define POLY

#include "acp.h"
#include "octree.h"

using namespace std;
using namespace acp;

//#define USETHREADS

unsigned nthreads (unsigned int n);

template<class N>
N cross (const PV3<N> &a, const PV3<N> &b, const int coord)
{
  int s = coord > 0 ? 1 : -1;
  switch (s*coord) {
  case 1:
    return s*(a.getY()*b.getZ() - a.getZ()*b.getY());
  case 2:
    return s*(a.getZ()*b.getX() - a.getX()*b.getZ());
  case 3:
    return s*(a.getX()*b.getY() - a.getY()*b.getX());
  default:
    cerr << "bad coord in cross" << endl;
    exit(0);
  }
}

class Plane;

class Point : public Object<PV3> {
  friend class Plane;
  friend class TrianglePlane;
  friend class EPPoint;
  friend class MidPoint;
  friend class CentroidPoint;
  friend class Vertex;
  friend class Edge;
  friend class Face;
  friend class Polyhedron;
  friend class Convolution;
 public:
  Point () {}
  Point (double x, double y, double z, bool perturb = true)
    : Object(PV3<double>(x, y, z), perturb) {}
  void getBBox (double *bbox);
  bool identical (Point *p);
  int order (Point *p);
  bool identicalI (Point *p);
  bool onLine (Point *a, Point *b);
  int side (Plane *a);
};

typedef vector<PTR<Point>> Points;

extern PTR<Point> Rdir;

class OnLine : public Primitive {
  Point *a, *t, *h;
  
  DeclareSign {
    PV3<N> u = a->get<N>() - t->get<N>(), v = h->get<N>() - t->get<N>();
    return Rdir->get<N>().tripleProduct(u, v);
  }
 public:
  OnLine (Point *a, Point *t, Point *h) : a(a), t(t), h(h) {}
};

class PointOrderE : public Primitive {
  Point *a, *t, *h;

  DeclareSign {
    return (a->get<N>() - t->get<N>()).dot(h->get<N>() - t->get<N>());
  }
 public:
  PointOrderE (Point *a, Point *t, Point *h) : a(a), t(t), h(h) {}
};

bool onEdge (Point *a, Point *t, Point *h, bool strict);

class PointOrder : public Primitive {
  Point *a, *b, *r;

  DeclareSign {
    return r->get<N>().dot(b->get<N>() - a->get<N>());
  }
 public:
  PointOrder (Point *a, Point *b, Point *r) : a(a), b(b), r(r) {}
};

class PointOrderR : public Primitive {
  Point *a, *b;

  DeclareSign {
    return Rdir->get<N>().dot(b->get<N>() - a->get<N>());
  }
 public:
  PointOrderR (Point *a, Point *b) : a(a), b(b) {}
};

class PointOrderPP : public Primitive {
  Point *a, *b, *t, *h;

  DeclareSign {
    return (h->get<N>() - t->get<N>()).dot(b->get<N>() - a->get<N>());
  }
 public:
  PointOrderPP (Point *a, Point *b, Point *t, Point *h)
    : a(a), b(b), t(t), h(h) {}
};

class LeftTurn : public Primitive {
  Point *a, *b, *c;
  int pc;

  DeclareSign {
    return cross(c->get<N>() - b->get<N>(), a->get<N>() - b->get<N>(), pc);
  }
 public:
  LeftTurn (Point *a, Point *b, Point *c, int pc) 
    : a(a), b(b), c(c), pc(pc) {}
};

class Orientation :public Primitive {
  Point *a, *b, *c, *d;

  DeclareSign {
    PV3<N> u = d->get<N>() - a->get<N>(), v = b->get<N>() - a->get<N>(),
      w = c->get<N>() - a->get<N>();
    return u.tripleProduct(w, v);
  }
 public:
  Orientation (Point *a, Point *b, Point *c, Point *d)
    : a(a), b(b), c(c), d(d) {}
};

namespace acp {
class CloserPair : public Primitive {
  Point *a, *b, *c, *d;

  DeclareSign {
    PV3<N> ab = a->get<N>() - b->get<N>(), cd = c->get<N>() - d->get<N>();
    return cd.dot(cd) - ab.dot(ab);
  }
 public:
  CloserPair (Point *a, Point *b, Point *c, Point *d)
    : a(a), b(b), c(c), d(d) {}
};
}
bool closerPair (Point *a, Point *b, Point *c, Point *d);

template<class N>
class PlaneData {
 public:
  PV3<N> n;
  N k;
  unsigned int size () const { return 4; }
  const N & operator[] (unsigned int i) const { return i < 3 ? n[i] : k; }
  N & operator[] (unsigned int i) { return i < 3 ? n[i] : k; }
  PlaneData () {}
  PlaneData (const PV3<N> &n, const N &k) : n(n), k(k) {}
  template<class M>
    PlaneData (const PlaneData<M> &p) : n(p.n), k(p.k) {}
};

class Plane : public Object<PlaneData> {
  friend class Point;
  friend class Face;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  int pc;
 public:
  int getPC () const { return pc; }
  bool coplanar (Plane *p);
};

class PlaneSideO : public Object<Scalar> {
  Plane *p;
  Point *a;

  DeclareCalculate (Scalar) {
    PlaneData<N> pp = p->get<N>();
    PV3<N> aa = a->get<N>();
    return pp.n.dot(aa) + pp.k;
  }
 public:
  PlaneSideO (Plane *p, Point *a) : p(p), a(a) {}
};

int projectionCoordinate (Plane *p);

class PlaneNormalCoord : public Primitive {
  Plane *p;
  unsigned int i;

  DeclareSign { return p->get<N>().n[i]; }
 public:
  PlaneNormalCoord (Plane *p, unsigned int i) : p(p), i(i) {}
};

class Coplanar : public Primitive {
  Plane *p, *q;

  DeclareSign {
    return Rdir->get<N>().tripleProduct(p->get<N>().n, q->get<N>().n);
  }
 public:
  Coplanar (Plane *p, Plane *q) : p(p), q(q) {}
};

class Side : public Primitive {
  Plane *p;
  Point *a;

  DeclareSign {
    return a->get<N>().dot(p->get<N>().n) + p->get<N>().k;
  }
 public:
  Side (Plane *p, Point *a) : p(p), a(a) {}
};

class PlaneRayAlignment : public Primitive {
  Plane *p;
  Point *u;

  DeclareSign {
    return p->get<N>().n.dot(u->get<N>());
  }
 public:
  PlaneRayAlignment (Plane *p, Point *u) : p(p), u(u) {}
};

double bboxSize (double *bb);

void copyBBox (const double *bbf, double *bbt);

void mergeBBox (const double *bbf, double *bbt);

bool bboxOverlap (const double *a, const double *b, double s = 0.0);

bool bboxOverlap (Point *a, const double *bbox);

class TrianglePlane : public Plane {
 protected:
  PTR<Point> a, b, c;

  DeclareCalculate(PlaneData) {
    PV3<N> ap = a->get<N>(), bp = b->get<N>(), cp = c->get<N>(),
      n = (cp - bp).cross(ap - bp);
    N k = - n.dot(ap);
    return PlaneData<N>(n, k);
  }
 public:
  TrianglePlane () : a(0), b(0), c(0) {}
  TrianglePlane (Point *a, Point *b, Point *c) : a(a), b(b), c(c) {
    pc = projectionCoordinate(this);
  }
  
  Point * getA () const { return a; }
  Point * getB () const { return b; }
  Point * getC () const { return c; }
};

class NegPoint : public Point {
  friend class Point;
  PTR<Point> a;
  
  DeclareCalculate(PV3) { return - a->get<N>(); }
 public:
  NegPoint (Point *a) : a(a) {}
};

class SumPoint : public Point {
  friend class Point;
  PTR<Point> a, b;
  
  DeclareCalculate (PV3) {
    return a->get<N>() + b->get<N>(); }
 public:
  SumPoint (Point *a, Point *b) : a(a), b(b) {}
};

class DiffPoint : public Point {
  friend class Point;
  PTR<Point> a, b;
  
  DeclareCalculate (PV3) { return a->get<N>() - b->get<N>(); }
 public:
  DiffPoint (Point *a, Point *b) : a(a), b(b) {}
};

class MidPoint : public Point {
  friend class Point;
  PTR<Point> a, b;
  double t;
  DeclareCalculate (PV3) { return a->get<N>()*(1-t) + b->get<N>()*t; }
 public:
  MidPoint (Point *p, Point *q) {
    MidPoint *pm = dynamic_cast<MidPoint*>(p);
    MidPoint *qm = dynamic_cast<MidPoint*>(q);
    if (pm != 0 && qm == 0) {
      if (q == pm->a) {
	a = pm->a;
	b = pm->b;
	t = pm->t / 2;
	return;
      }
      if (q == pm->b) {
	a = pm->a;
	b = pm->b;
	t = (1+pm->t) / 2;
	return;
      }
    } else if (pm == 0 && qm != 0) {
      if (p == qm->a) {
	a = qm->a;
	b = qm->b;
	t = qm->t / 2;
	return;
      }
      if (p == qm->b) {
	a = qm->a;
	b = qm->b;
	t = (1+qm->t) / 2;
	return;
      }
    }
    else if (pm != 0 && qm != 0 && pm->a == qm->a && pm->b == qm->b) {
      a = pm->a;
      b = pm->b;
      t = (pm->t + qm->t) / 2;
    }

    a = p;
    b = q;
    t = 0.5;
  }
};

class CentroidPoint : public Point {
  friend class Point;
  Points pts;

  DeclareCalculate (PV3) {
    PV3<N> a = pts[0]->get<N>();
    for (int i = 1; i < pts.size(); ++i)
      a = a + pts[i]->get<N>();
    return a/pts.size();
  }
 public:
  CentroidPoint (const Points &ipts) {
    for (int i = 0; i < ipts.size(); ++i)
      pts.push_back(ipts[i]);
  }
  
  CentroidPoint (Point *a, Point *b) {
    pts.push_back(a);
    pts.push_back(b);
  }

  CentroidPoint (Point *a, Point *b, Point *c) {
    pts.push_back(a);
    pts.push_back(b);
    pts.push_back(c);
  }

  CentroidPoint (Point *a, Point *b, Point *c, Point *d) {
    pts.push_back(a);
    pts.push_back(b);
    pts.push_back(c);
    pts.push_back(d);
  }  
};

class ScalePoint : public Point {
  PTR<Point> p;
  double unit;

  DeclareCalculate (PV3) { return p->get<N>()*unit; }
 public:
  ScalePoint(Point *p, double unit) : p(p), unit(unit) {}
};

class EPPoint : public Point {
  PTR<Point> t, h;
  PTR<TrianglePlane> p;
  
  DeclareCalculate (PV3) {
    PV3<N> a = p->getA()->get<N>(), n = p->get<N>().n, tp = t->get<N>(),
      u = h->get<N>() - tp;
    N k = n.dot(a - tp)/n.dot(u);
    return tp + k*u;
  }
 public:
  EPPoint (Point *t, Point *h, TrianglePlane *p) : t(t), h(h), p(p) {}
  TrianglePlane * getP () { return p; }
};

class RayPlanePoint : public Point {
  PTR<Point> t, r;
  TrianglePlane *p;

  DeclareCalculate (PV3) {
    PV3<N> a = t->get<N>(), u = r->get<N>(), n = p->get<N>().n;
    N k = - (n.dot(a) + p->get<N>().k)/n.dot(u);
    return a + k*u;
  }
 public:
  RayPlanePoint (Point *t, Point *r, TrianglePlane *p) : t(t), r(r), p(p) {}
};

class RayZPlanePoint : public Point {
  PTR<Point> t, r;
  double z;

  DeclareCalculate (PV3) {
    PV3<N> a = t->get<N>(), u = r->get<N>();
    N k = (z - a.z)/u.z;
    return a + k*u;
  }
 public:
  RayZPlanePoint (Point *t, Point *r, double z)
    : t(t), r(r), z(z) {}
};

class Edge;
class HEdge;
class Face;

class Vertex {
  friend class Edge;
  friend class HEdge;
  friend class Face;
  friend class Shell;
  friend class Cell;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  PTR<Point> p;
  vector<Edge *> edges;
  double bbox[6];
 public:
  Vertex (Point *p) : p(p) { p->getBBox(bbox); }
  Point * getP () { return p; }
  int EdgesN () const { return edges.size(); }
  Edge * getEdge (int i) const { return edges[i]; }
  double * getBBox () { return bbox; }
  void getBBox (double *b) const { copyBBox(bbox, b); }
  vector<HEdge *> outgoingHEdges () const;
  vector<Face *> incidentFaces () const;
  HEdge * connected (Vertex *a) const;
};

typedef vector<Vertex *> Vertices;

class Triangle {
 public:
  PTR<Point> a, b, c;
  Triangle (Point *a, Point *b, Point *c) : a(a), b(b), c(c) {}
};

typedef vector<Triangle> Triangles;

class Vertex;

class HEdge {
  friend class Edge;
  friend class Face;
  friend class HFace;
  friend class Shell;
  friend class Cell;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  Edge *e;
  bool forward, flag;
  Face *f;
  HEdge *next;
 public:
  HEdge () : e(0), forward(false), next(0), f(0) {}
  HEdge (Edge *e, bool forward) 
    : e(e), forward(forward), next(0), f(0) {}
  Vertex * tail () const;
  Vertex * head () const;
  Edge * getE () const { return e; }
  bool getForward () const { return forward; }
  Face * getF () const { return f; }
  void setF (Face *g) { f = g; }
  HEdge * getNext () const { return next; }
  void setNext (HEdge *h) { next = h; }
  template<class N> PV3<N> getU ();
  template<class N> PV3<N> getN ();
  HEdge * cw () const;
  HEdge * ccw () const;
  Vertices loop () const;
  Points pointLoop () const;
  vector<HEdge *> edgeLoop ();
};

typedef vector<HEdge *> HEdges;

class Edge {
  friend class EEPoint;
  friend class Vertex;
  friend class HEdge;
  friend class Face;
  friend class Shell;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  Vertex *t, *h;
  double bbox[6];
  vector<HEdge *> hedges;
  Vertices vertices;
  vector <Edge *> dedges;
 public:
  Edge (Vertex *t, Vertex *h);
  void setBBox ();
  ~Edge ();
  HEdge * addHEdge (bool forward);
  void removeHEdge (HEdge *e);
  Vertex * getT () const { return t; }
  Vertex * getH () const { return h; }
  template<class N> PV3<N> getU () { return h->p->get<N>() - t->p->get<N>(); }
  int HEdgesN () const { return hedges.size(); }
  HEdge * getHEdge (int i) const { return i < hedges.size() ? hedges[i] : 0; }
  double * getBBox () { return bbox; }
  void sortHEdges ();
};

typedef vector<Edge *> Edges;

template<class N>
PV3<N> HEdge::getU ()
{
  return forward ? e->getU<N>() : - e->getU<N>();
}

class EEPoint : public Point {
  PTR<Point> et, eh, ft, fh;
  int pc;
  vector<Face *> fa;
  
  DeclareCalculate (PV3) {
    PV3<N> a = et->get<N>(), u = eh->get<N>() - a, b = ft->get<N>(),
      v = fh->get<N>() - b;
    N k = cross(b - a, v, pc)/cross(u, v, pc);
    return a + k*u;
  }
 public:
  EEPoint (Edge *e, Edge *f, int pc) : et(e->t->getP()), eh(e->h->getP()), 
    ft(f->t->getP()), fh(f->h->getP()), pc(pc) {}

  EEPoint (Point *et, Point *eh, Point *ft, Point *fh, int pc)
    : et(et), eh(eh), ft(ft), fh(fh), pc(pc) {}

  void addFace (Face *f) { fa.push_back(f); }
  vector<Face *> getFaces () { return fa; }
};

class Shell;

class HFace {
  friend class Face;
  friend class Shell;
  friend class Cell;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  Face *f;
  Shell *s;
 public:
  HFace () : f(0), s(0) {}
  Face * getF () const { return f; }
  Shell * getS () const { return s; }
  bool pos () const;
  HFace * twin () const;
  template<class N> PV3<N> getN () const;
  vector<HFace *> neighbors () const;
  HFace * neighbor (HEdge *e) const;
};

typedef vector<HFace *> HFaces;

class Face {
  friend class HEdge;
  friend class HFace;
  friend class Shell;
  friend class Cell;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  HEdge *h;
  PTR<TrianglePlane> p;
  HFace hfaces[2];
  double bbox[6];
 public:
  Face (HEdge *h, TrianglePlane *p = 0);
  Face (Point *a, Point *b, Point *c);
  void addLoop (HEdge *h, bool flag);
  void update ();
  TrianglePlane * getP () { return p; }
  HEdge * getBoundary () const { return h; }
  double * getBBox () { return bbox; }
  void getBBox (double *b) const { copyBBox(bbox, b); }
  HFace * getHFace (int i) { return hfaces + i; }
  bool boundaryVertex (Point *a) const;
  bool boundaryVertex (Vertex *v) const;
  bool boundaryEdge (Edge *e) const;
  bool sharedEdge (Face *f) const;
  void sharedVertices (Face *f, Points &pts) const;
  Points boundaryPoints () const;
  virtual HEdges boundaryHEdges () const { return h->edgeLoop(); }
  bool intersects (Face *g, bool strict);
  bool intersectsFP (Face *g, int *sg);
  bool checkFP (int *s, int n) const;
  bool verifyFP (Face *g, int *sg);
  bool intersectsFE (Face *g, int *sg, bool strict);
  bool intersectsFEP (Point *et, Point *eh, bool strict);
  bool intersectsEE (Point *et, Point *eh, Point *ft, Point *fh, bool strict);
  bool contains (Point *a, bool strict, int *ie = 0);
  Point * centroid () const;
  PTR<Point> rayIntersection (Point *a, Point *r);
  virtual void triangulate (Triangles &tr);
};

typedef vector<Face *> Faces;

template<class N>
PV3<N> HEdge::getN ()
{
  return forward ? f->p->get<N>().n : - f->p->get<N>().n; 
}

template<class N>
PV3<N> HFace::getN () const
{
  return pos() ? f->p->get<N>().n : - f->p->get<N>().n;
}

class EdgeOrder1 : public Primitive {
  HEdge *e;

  DeclareSign { return Rdir->get<N>().dot(e->getN<N>()); }
 public:
  EdgeOrder1 (HEdge *e) : e(e) {}
};

class EdgeOrder2 : public Primitive {
  Edge *e;
  HEdge *f, *g;

  DeclareSign { return e->getU<N>().tripleProduct(g->getN<N>(), f->getN<N>()); }
 public:
  EdgeOrder2 (Edge *e, HEdge *f, HEdge *g) : e(e), f(f), g(g) {}
};

class EdgeOrder {
  Edge *e;
 public:
  EdgeOrder (Edge *e) : e(e) {}
  bool operator() (HEdge *f, HEdge *g) const {
    if (f == g) return false;
    int sf = EdgeOrder1(f), sg = EdgeOrder1(g);
    return sf != sg ? sf == 1 : EdgeOrder2(e, f, g) == 1;
  }
};

class Cross2 : public Primitive {
  Point *a, *b;
  int c;

  DeclareSign { return cross(a->get<N>(), b->get<N>(), c); }
 public:
  Cross2 (Point *a, Point *b, int c) : a(a), b(b), c(c) {}
};

bool rayEdgeIntersection (Point *a, Point *r, Point *t, Point *h, int c);

class Cell;

class Shell {
  friend class Cell;
  friend class Polyhedron;
  friend class Convolution;
 protected:
  HFaces hfaces;
  double bbox[6];
  Octree<Face *> *octreef;
  Cell *c;
  void setBBox ();
  void setOctree ();
 public:
  Shell () : octreef(0), c(0) {};
  ~Shell ();
  Shell (const HFaces &hf);
  HFaces & getHFaces () { return hfaces; }
  double * getBBox () { return bbox; }
  Cell * getC() const { return c; }
  bool outer () const;
  Vertex * vmax (Point *r) const;
  Vertex * vmax (Face *f, Point *r) const;
  bool contains (Shell *s) const;
  bool subset (const Shell *s) const;
  int contains (Point *a) const;
  void rayBBox (Point *a, Point *r, double *rb) const;
  int euler () const;
  bool pos () const;
};

typedef vector<Shell *> Shells;

void deleteShells (const Shells &sh);

class SlopeOrder : public Primitive {
  Edge *e, *f;
  Point *x;

  DeclareSign {
    PV3<N> u = e->getU<N>(), v = f->getU<N>(), r = x->get<N>();
    N ur = u.dot(r), vr = v.dot(r);
    return vr*vr*u.dot(u) - ur*ur*v.dot(v);
  }
 public:
  SlopeOrder (Edge *e, Edge *f, Point *x) : e(e), f(f), x(x) {}
};

class Convex : public Primitive {
  HEdge *e;
  HFace *f;

  DeclareSign {
    HFace *g = f->neighbor(e);
    Face *ff = f->getF(), *gf = g->getF();
    PV3<N> nf = ff->getP()->get<N>().n, ng = g->getN<N>(), u = e->getU<N>();
    return u.tripleProduct(ng, nf);
  }
 public:
  Convex (HEdge *e, HFace *f) : e(e), f(f) {}
};

class HFaceNormal : public Point {
  PTR<TrianglePlane> p;
  bool pos;

  DeclareCalculate (PV3) { return pos ? p->get<N>().n : - p->get<N>().n; }
 public:
  HFaceNormal (HFace *f) : p(f->getF()->getP()), pos(f->pos()) {}
};

class Cell {
  friend class Polyhedron;
  friend class Convolution;
 protected:
  Shell *outer;
  Shells inner;
  int wn;
 public:
  Cell (Shell *outer) : outer(outer), wn(0) { if (outer) outer->c = this; }
  ~Cell () { delete outer; deleteShells(inner); }
  int nShells () const { return inner.size() + (outer ? 1 : 0); }
  
  Shell * getShell (int i) const {
    if (outer)
      return i == 0 ? outer : inner[i-1];
    return inner[i];
  }
    
  double * getBBox () { return outer ? outer->getBBox() : 0; }
  void getBBox (double *b) const { copyBBox(outer->bbox, b); }
  void addInner (Shell *s) { inner.push_back(s); s->c = this; }
  bool contains (Point *p) const;
  Point * interiorPoint () const;
  int getWN () const { return wn; }
  PTR<Point> interiorPoint (Vertex *a, double d) const;
};

typedef vector<Cell *> Cells;

HFace * largestHFace (const HFaces &hf);

double areaSquared (Face *f);

typedef pair<Point *, Vertex*> PVPair;

typedef map<Point *, Vertex *> PVMap;

enum SetOp { Union, Intersection, Complement };

class Polyhedron : public RefCnt {
 public:
  Vertices vertices;
  Edges edges;
  Faces faces;
  Cells cells;
  double bbox[6];

  ~Polyhedron ();
  Vertex * getVertex (Point *p);
  Vertex * getVertex (double x, double y, double z, bool perturbed = true) {
    return getVertex(new Point(x, y, z, perturbed));
  }
  Vertex * getVertex (Point *p, PVMap &pvmap);
  Edge * getEdge (Vertex *a, Vertex *b);
  HEdge * addHEdge (Vertex *a, Vertex *b);
  HEdge * getHEdge (Vertex *a, Vertex *b);
  Face * addTriangle (Vertex *a, Vertex *b, Vertex *c, TrianglePlane *p = 0);
  Face * addRectangle (Vertex *a, Vertex *b, Vertex *c, Vertex *d);
  Face * addFace (HEdge *h);
  void formCells ();
  void formCellsAux (const Shells &shells);
  void removeShell (Shell *s);
  void formShells (Shells &shells);
  Shell * formShell (HFace *f) const;
  Cell * enclosingCell (Shell *s, Octree<Cell *> *octreec) const;
  void clearCells ();
  Face * addTriangle (Point *a, Point *b, Point *c, PVMap &pvmap, TrianglePlane *p = 0);
  Face * addTriangle (Face *f, PVMap &pvmap, TrianglePlane *p = 0);
  Polyhedron * copy () const;
  Polyhedron * scale (double unit) const;
  Polyhedron * negative () const;
  Polyhedron * translate (Point *t) const;
  Polyhedron * negativeTranslate (Point *t) const;
  bool intersects (Polyhedron *a, bool strict, Octree<Face *> *aoctree = 0) const;
  bool contains (Point *p) const;
  int containingCell (Point *p) const;
  bool intersectsEdges (Octree<Face *> *octree, bool strict) const;
  Polyhedron * boolean (Polyhedron *a, SetOp op);
  Polyhedron * mergeCells (const set<Cell *> &cin) const;
  Polyhedron * selfUnion ();
  Polyhedron * cellPolyhedron (int i) const;
  void addHFaces (const HFaces &hf, PVMap &pvmap);
  void replaceVertex (Face *f, Vertex *v, Vertex *w);
  void removeLoop (HEdge *e);
  void removeHEdge (HEdge *e);
  void moveVertex (Vertex *v, Point *p);
  void removeNullFaces ();
  Octree<Face *> * faceOctree (double s = 0.0) const;
  Octree<Cell *> * cellOctree () const;
  void computeWindingNumbers ();
  void updateWN (Cell *c, HFaces &st) const;
  void describe ();
};

typedef vector<Polyhedron *> Polyhedrons;

Face * faceVertices (Vertex *a, Vertex *b, Vertex *c);

Face * faceVertices (Vertex *a, Vertex *b, Vertex *c, Vertex *d);

Polyhedron * subdivide (Polyhedron *a, bool oneway);

class FFE {
 public:
  Face *f, *g;
  Points *pts;
  double bbox[6];

  FFE () : f(0), g(0), pts(0) {}

  FFE (Face *f, Face *g, Point *p, Point *q)
    : f(f), g(g), pts(new Points) {
    pts->push_back(p); pts->push_back(q);
    p->getBBox(bbox);
    double bbq[6];
    q->getBBox(bbq);
    mergeBBox(bbq, bbox);
  }

  FFE (Face *f, Face *g, Points *pts) : f(f), g(g), pts(pts) {
    pts->at(0)->getBBox(bbox);
    double bbq[6];
    pts->back()->getBBox(bbq);
    mergeBBox(bbq, bbox);
  }
};

void intersectFF (Polyhedron *a,
		  map<Edge *, Points *> &epsmap,
		  map<pair<Face *, Face *>, FFE> &ffemap,
		  vector<pair<Face *, Edge *>> &fe,
		  vector<pair<Edge *, Edge *>> &ee);

void intersectFE (const vector<pair<Face *, Face *>> &ff1,
		  map<pair<Face *, Edge *>, PTR<Point>> &fepmap,
		  map<Edge *, Points *> &epsmap,
		  map<pair<Face *, Face *>, Points> &ffpsmap,
		  vector<pair<Face *, Face *>> &ff,
		  vector<pair<Face *, Edge *>> &fe,
		  vector<pair<Edge *, Edge *>> &ee);

void intersectFET (void *ptr);

class FEData {
 public:
  unsigned int i, is, ie;
  const vector<pair<Face *, Face *>> *ff1;
  vector<pair<Edge *, PTR<Point>>> ep;
  vector<pair<pair<Face *, Edge *>, PTR<Point>>> fep;
  vector<pair<pair<Face *, Face *>, PTR<Point>>> ffp;
  vector<pair<Face *, Face *>> ff;
  vector<pair<Face *, Edge *>> fe;
  vector<pair<Edge *, Edge *>> ee;
};

void intersectFE (Face *f, Face *g,
		  vector<pair<Edge *, PTR<Point>>> &ep,
		  vector<pair<pair<Face *, Edge *>, PTR<Point>>> &fep,
		  vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp,
		  vector<pair<Face *, Face *>> &ff,
		  vector<pair<Face *, Edge *>> &fe,
		  vector<pair<Edge *, Edge *>> &ee);		  

void intersectFE (Face *f, Face *g, int *sg,
		  vector<pair<Edge *, PTR<Point>>> &ep,
		  vector<pair<pair<Face *, Edge *>, PTR<Point>>> &fep,
		  vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp,
		  vector<pair<Face *, Edge *>> &fe,
		  vector<pair<Edge *, Edge *>> &ee);

void intersectFEP (Face *f, HEdge *h,
		   vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp,
		   vector<pair<Edge *, PTR<Point>>> &ep,
		   vector<pair<Face *, Edge *>> &fe,
		   vector<pair<Edge *, Edge *>> &ee);

bool intersectEE (Edge *e, Edge *f, int c,
		   vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp,
		  vector<pair<Edge *, PTR<Point>>> &ep);

bool intersectEE (Point *et, Point *eh, Point *ft, Point *fh,
		  int c, Points &pe, Points &pf);

void intersectFV (Face *f, HEdge *h,
		  vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp,
		  vector<pair<Edge *, PTR<Point>>> &ep);

void updateFFP (Face *f, Vertex *v,
		vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp);

void updateFFP (Edge *e1, Edge *e2, Point *p,
		vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp);

void intersectFEG (Face *f, HEdge *h,
		   vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp,
		   vector<pair<Edge *, PTR<Point>>> &ep,
		   vector<pair<pair<Face *, Edge *>, PTR<Point>>> &fep);

bool first (HEdge *h);

bool signChange (int *s);

void update (Edge *e, Point *p, map<Edge *, Points *> &epsmap);

void update (const pair<Face *, Face *> &ff, Point *p,
	     map<pair<Face *, Face *>, Points> &fppmap);

void intersectFF (const vector<pair<Face *, Face *>> &ff,
		  const map<pair<Face *, Edge *>, PTR<Point>> &fepmap,
		  const map<Edge *, Points *> &epsmap,
		  map<pair<Face *, Face *>, Points> &ffpsmap,
		  map<pair<Face *, Face *>, FFE> &ffemap);

class FFData {
 public:
  unsigned int i, is, ie;
  const vector<pair<Face *, Face *>> *ff;
  const map<pair<Face *, Edge *>, PTR<Point>> *fepmap;
  const map<Edge *, Points *> *epsmap;
  const map<pair<Face *, Face *>, Points> *ffpsmap;
  vector<pair< pair<Face *, Face *>, pair<PTR<Point>, PTR<Point>>>> ffpp;
};

void intersectFFT (void *ptr);

void intersectFF (Face *f, Face *g,
		  const map<pair<Face *, Edge *>, PTR<Point>> &fepmap,
		  const map<Edge *, Points *> &epsmap,
		  const map<pair<Face *, Face *>, Points> &ffpsmap,
		  vector<pair< pair<Face *, Face *>, pair<PTR<Point>, PTR<Point>>>> &ffpp);

void findFE (Face *f, Face *g,
	     const map<pair<Face *, Edge *>, PTR<Point>> &fepmap,
	     Points &pts);

void sharedVertices (Face *f, Face *g,
		     const map<pair<Face *, Face *>, Points> &ffpsmap, Points &pts);

void removeDuplicates (Points &pts);

class FFOrder : public Primitive {
  Face *f, *g;
  Point *p, *q;

  DeclareSign {
    PV3<N> u = f->getP()->get<N>().n.cross(g->getP()->get<N>().n),
      pq = p->get<N>() - q->get<N>();
    return u.dot(pq);
  }
 public:
  FFOrder (Face *f, Face *g, Point *p, Point *q) : f(f), g(g), p(p), q(q) {}
};

void update (Face *f, Face *g, Point *a, Point *b,
	     map<pair<Face *, Face *>, FFE> &ffemap);

map<Face *, vector<FFE>> feMap (map<Edge *, Points *> &epsmap,
				map<pair<Face *, Face *>, FFE> &ffemap,
				const vector<pair<Face *, Edge *>> &fe);

void feMapFE (map<Edge *, Points *> &epsmap,
	      map<pair<Face *, Face *>, FFE> &ffemap,
	      Face *f, Edge *e, map<Face *, vector<FFE>> &femap);

class PointOrderRP {
 public:
  bool operator() (Point *a, Point *b) const {
    return a != b && PointOrderR(a, b) == 1;
  }
};

void update (Face *f, const FFE &ffe, map<Face *, vector<FFE>> &femap);

void intersectFFF (const map<pair<Face *, Face *>, FFE> &ffemap,
		   const map<Face *, vector<FFE>> &femap,
		   vector<pair<Points *, Points *>> &col);

void intersectFFFAux (const map<pair<Face *, Face *>, FFE> &ffemap,
		      const map<Face *, vector<FFE>> &femap,
		      vector<pair<Points *, Points *>> &col,
		      vector<pair<Points *, PTR<Point>>> &psps);

class FFFData {
 public:
  unsigned int i, is, ie;
  const map<pair<Face *, Face *>, FFE> *ffemap;
  const map<Face *, vector<FFE>> *femap;
  vector<pair<Points *, Points *>> col;
  vector<pair<Points *, PTR<Point>>> psps;
};

void intersectFFFT (void *ptr);

void intersectFFF (Face *f, const vector<FFE> &ed,
		   const map<pair<Face *, Face *>, FFE> &ffemap,
		   vector<pair<Points *, Points *>> &col,
		   vector<pair<Points *, PTR<Point>>> &psps);

class PointOrderPPP {
 public:
  Point *t, *h;
    
  PointOrderPPP (Point *t, Point *h) : t(t), h(h) {}
    
  bool operator() (Point *a, Point *b) const {
    return a != b && PointOrderPP(a, b, t, h) == 1;
  }
};

void subedgesE (map<Edge *, Points *> &epsmap, bool oneway, set<Point *> &bpts,
		vector<pair<PTR<Point>, PTR<Point>>> &equiv);

class SEEData {
 public:
  unsigned int i, is, ie;
  const map<Edge *, Points *> *epsmap;
  bool oneway;
  set<Point *> bpts;
  vector<pair<PTR<Point>, PTR<Point>>> equiv;
};

void subedgesET (void *ptr);

void subedgesE (Edge *e, Points &pts, bool oneway, set<Point *> &bpts,
		vector<pair<PTR<Point>, PTR<Point>>> &equiv);

void subedgesE (Points &pts, bool oneway, set<Point *> &bpts,
		vector<pair<PTR<Point>, PTR<Point>>> &equiv);

void subedgesE1 (Points &pts, set<Point *> &bpts,
		 vector<pair<PTR<Point>, PTR<Point>>> &equiv);

void subedgesP (Points &pts, vector<pair<PTR<Point>, PTR<Point>>> &equiv);

void subedgesFF (map<pair<Face *, Face *>, FFE> &ffemap, bool oneway,
		 const set<Point *> &bpts,
		 vector<pair<PTR<Point>, PTR<Point>>> &equiv);

class SEFData {
 public:
  unsigned int i, is, ie;
  map<pair<Face *, Face *>, FFE> *ffemap;
  bool oneway;
  const set<Point *> *bpts;
  vector<pair<PTR<Point>, PTR<Point>>> equiv;
};

void subedgesFFT (void *ptr);

void subedgesFF (FFE &ffe, bool oneway, const set<Point *> &bpts,
		 vector<pair<PTR<Point>, PTR<Point>>> &equiv);

void subedgesFF (FFE &ffe, const set<Point *> &bpts,
		 vector<pair<PTR<Point>, PTR<Point>>> &equiv);

Face * otherFace (Face *f, Face *g, const Faces &fa);

map<PTR<Point>, PTR<Point>> pMap (const map<Edge *, Points *> &epsmap,
				  const vector<pair<Edge *, Edge *>> &ee,
				  const vector<pair<Points *, Points *>> &col);

void pMapEE (Edge *e, Edge *f, const map<Edge *, Points *> &epsmap,
	     map<PTR<Point>, PTR<Point>> &pmap);

void pMapCol (const Points &pts1, const Points &pts2, 
	      map<PTR<Point>, PTR<Point>> &pmap);

void pMapEquiv (Point *a, Point *b, map<PTR<Point>, PTR<Point>> &pmap);

Point * getPoint (Point *p, const map<PTR<Point>, PTR<Point>> &pmap);

void update (const vector<pair<PTR<Point>, PTR<Point>>> &equiv,
	      map<PTR<Point>, PTR<Point>> &pmap);

class SEdge {
 public:
  PTR<Point> tail;
  bool forward, flag;
  SEdge *twin, *cw;

  SEdge (Point *tail, bool forward, SEdge *twin = 0)
    : tail(tail), forward(forward), flag(false), twin(twin), cw(0) {}
  Point * head () const { return twin->tail; }
  SEdge * next () const { return twin->cw; }
  Points loop () const;
};  

typedef vector<SEdge *> SEdges;

class SFace {
 public:
  Face *f;
  vector<Points> b;

  SFace (Face *f, const Points &pts) : f(f) { b.push_back(pts); }
};

typedef vector<SFace> SFaces;

Polyhedron * subfaces (Polyhedron *a, bool oneway,
		       const map<Edge *, Points *> &epsmap,
		       const map<Face *, vector<FFE>> &femap,
		       const map<PTR<Point>, PTR<Point>> &pmap);

void subfacesT (void *ptr);

class SUData {
 public:
  unsigned int i, is, ie;
  Polyhedron *a;
  bool oneway;
  const map<Edge *, Points *> *epsmap;
  const map<Face *, vector<FFE>> *femap;
  const map<PTR<Point>, PTR<Point>> *pmap;
  SFaces sf;
};

void subfaces (Face *f, bool oneway,
	       const map<Edge *, Points *> &epsmap,
	       const map<Face *, vector<FFE>> &femap,
	       const map<PTR<Point>, PTR<Point>> &pmap, SFaces &sf);

SEdges subedges (Face *f, bool oneway,
		 const map<Edge *, Points *> &epsmap,
		 const map<Face *, vector<FFE>> &femap,
		 const map<PTR<Point>, PTR<Point>> &pmap);

void subedgesH (HEdge *h, const map<Edge *, Points *> &epsmap,
		vector<pair<PTR<Point>, PTR<Point>>> &pp);

void subedges (const Points &pts, bool fflag, bool bflag,  
	       vector<pair<PTR<Point>, PTR<Point>>> &pp);

void subedgesFFE (Face *f, const FFE &ffe, bool oneway,
		  vector<pair<PTR<Point>, PTR<Point>>> &pp);

class SEdgeCWOrder {
  int c;
 public:
  SEdgeCWOrder (int c) : c(c) {}
  bool operator() (SEdge *e, SEdge *f) const {
    Point *et = e->tail, *ft = f->tail;
    if (et != ft)
      return et < ft;
    Point *eh = e->head(), *fh = f->head();
    if (eh == fh)
      return f->forward < e->forward;
    int s1 = PointOrderR(et, eh), s2 = PointOrderR(et, fh);
    return s1 != s2 ? s1 == 1 : LeftTurn(eh, et, fh, c) == 1;
  }
};

class SEdgeHeadOrder {
 public:
  bool operator() (SEdge *e, SEdge *f) const {
    Point *p = e->head(), *q = f->head();
    return p != q && PointOrderR(p, q) == 1;
  }
};

SEdges subedgesPP (const vector<pair<PTR<Point>, PTR<Point>>> &pp, int c,
		   const map<PTR<Point>, PTR<Point>> &pmap);

void sedges (Point *t, Point *h, SEdges &ed);

void subfaces (Face *f, bool oneway, SEdges &ed, SFaces &sf);

void linkSEdges (SEdges &ed, int c);

SEdge * findLoop (SEdge *h);

void addInner (SEdge *e, bool oneway, int c, SFaces &sf);

bool contains (SEdge *e, int c, const Points &pts);

class MFace : public Face {
  HEdges inner;
  
 public:
  MFace (HEdge *h, TrianglePlane *p) : Face(h, p) {}
  const HEdges & getInner () const { return inner; }
  void addInner (HEdge *e) { inner.push_back(e); addLoop(e, false); }
  void update (HEdge *e);
  virtual HEdges boundaryHEdges () const;
  virtual void triangulate (Triangles &tr);
};

void triangulate (const vector<Points *> &reg, int coord, Triangles &tr);

void addFace (const SFace &f, Polyhedron *a, PVMap &pvmap);

Vertices loop (const Points &pts, Polyhedron *a, PVMap &pvmap);

bool currentFace (const Vertices &ve, Polyhedron *a);

HEdge * addLoop (const Vertices &ve, Polyhedron *a);

bool outerLoop (HEdge *h);

void deleteMaps (map<Edge *, Points *> &epsmap,
		 map<pair<Face *, Face *>, FFE> &ffemap);

Polyhedron * triangulate (Polyhedron *a);

Triangles triangulate (const Faces &fa);

void triangulateT (void *ptr);

class TRData {
 public:
  unsigned int i, is, ie;
  const Faces *fa;
  Triangles tr;
};

class PointProjection : public Object<Scalar> {
  Point *a, *r;

  DeclareCalculate (Scalar) { return a->get<N>().dot(r->get<N>()); }
 public:
  PointProjection (Point *a, Point *r) : a(a), r(r) {}
};

class PointR {
 public:
  Point *a;
  double l, u;

  PointR (Point *a) : a(a) {
    PointProjection r(a, Rdir);
    Scalar<Interval> k = r.getApprox(1.0);
    l = k.x.l;
    u = k.x.u;
  }

  bool operator< (const PointR &x) const {
    if (a == x.a || x.u < l)
      return false;
    return u < x.l || a->order(x.a) == 1;
  }

  bool operator== (const PointR &x) const {
    if (a == x.a)
      return true;
    if (x.u < l || u < x.l)
      return false;
    return a->order(x.a) == 0;
  }
};

PVMap pvMap (const Points &pts, Polyhedron *a);

bool inSet (bool ina, bool inb, SetOp op);

Polyhedron * overlay (Polyhedron **poly, int n);

Polyhedron * multiUnion (Polyhedron **poly, int n);

Polyhedron * coalesce (Polyhedron *a);

Polyhedron * coalesceFaces (Polyhedron *a);

vector<set<Face *>> groupFaces (Polyhedron *a);

void coalesceFace (const set<Face *> &fs, SFaces &sf);

void coalesceEdges (Polyhedron *a); 

void coalesceLoop (HEdge *h, vector<pair<HEdge *, HEdge *>> &hh);

bool coalesceEdge (HEdge *e, HEdge *f);

void coalesceEdges (Polyhedron *a, HEdge *s, HEdge *e);

Polyhedron * box (double *b, bool perturb = true);

Polyhedron * sphere (double ox, double oy, double oz, double r, double err);

Polyhedron * sphere (double err);

Polyhedron * octohedron ();

double sphereError (Polyhedron *a);

Polyhedron * sphereRefine (Polyhedron *a);

Polyhedron * randomTets (int n, double u, double v);

void randomTet (Polyhedron *a, double u, double v);

Point * randomPoint (double d);

// chloe


class RayPoint : public Point {
  PTR<Point> a, b;
  double s;

  DeclareCalculate (PV3) {
    PV3<N> p = a->get<N>(), u = b->get<N>() - p;
    return p + u*s;
  }
 public:
  RayPoint (Point *a, Point *b, double d) : a(a), b(b) {
    PV3<double> u = b->getApproxMid() - a->getApproxMid();
    double uu = sqrt(u.dot(u)), k = min(0.5*uu, d);
    s = k*randomNumber(31.0/32.0, 33.0/32.0);
  }
};

static PV3<Parameter> PV3ItoPV3P (const PV3<Interval> &pv3i) {
  return PV3<Parameter>(Parameter(pv3i.x.l, pv3i.x.u),
                        Parameter(pv3i.y.l, pv3i.y.u),
                        Parameter(pv3i.z.l, pv3i.z.u));
}

class BisectorPoint : public Point {
  HFace *f1, *f2;
  HEdge *e;
  double k1, k2;

  DeclareCalculate (PV3) {
    PV3<N> n1 = f1->getN<N>(), n2 = f2->getN<N>(), u = e->getU<N>();
    if (u.tripleProduct(n1, n2).sign() == 0)
      return n1;
    return n1*k1 + n2*k2;
  }
 public:
  BisectorPoint (HFace *f1, HFace *f2, HEdge *e) : f1(f1), f2(f2), e(e) {
    PV3<Parameter> n1 = PV3ItoPV3P(f1->getF()->getP()->getApprox().n);
    PV3<Parameter> n2 = PV3ItoPV3P(f2->getF()->getP()->getApprox().n);
    double s1 = f1->pos() ? 1.0 : -1.0, s2 = f2->pos() ? 1.0 : -1.0,
      d1 = ::sqrt(n1.dot(n1).ub()), d2 = ::sqrt(n2.dot(n2).ub());
    k1 = 0.5*s1/max(1e-16, d1);
    k2 = 0.5*s2/max(1e-16, d2);
  }
};

#endif
