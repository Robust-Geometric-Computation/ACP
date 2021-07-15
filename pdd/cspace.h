//  -*- Mode: c++ -*-
#ifndef CSPACE
#define CSPACE

#include "mink.h"
#include "root.h"

#define PI 3.141592653589793

typedef set<Vertex *> VertexSet;

VertexSet intersection (const VertexSet &a, const VertexSet &b);

class Angle : public Object<PV2> {
  VertexSet vs;
protected:
  Angle () {}
  Angle (const VertexSet &vs) : vs(vs) {}
  Angle (double c, double s) : Object(PV2<double>(c, s), false, true) {}
public:
  static PTR<Angle> mpi, ppi;
  VertexSet getVS () const { return vs; }  
  bool order (Angle *b);

  template<class N> N dt () {
    if (this == Angle::mpi)
      return N(1);
    PV2<N> p = get<N>();
    N s = p.y.sign() == -1 ? p.x/(1.0 - p.y) : - p.x/(1.0 + p.y);
    return 2.0/(1.0 + s*s);
  }

  double theta () {
    PV2<double> p = getApproxMid();
    return atan2(p.y, p.x);
  }
};

typedef vector<PTR<Angle>> Angles;

class AngleUpper : public Primitive {
  Angle *a;

  DeclareSign {
    if (a == Angle::mpi) return N(-1);
    if (a == Angle::ppi) return N(1);
    return a->get<N>().y;
  }
public:
  AngleUpper (Angle *a) : a(a) {}
};

class AngleOrder : public Primitive {
  Angle *a, *b;

  DeclareSign {
    PV2<N> p = a->get<N>(), q = b->get<N>();
    return p.y.sign() != q.y.sign() ? q.y : p.cross(q);
  }
public:
  AngleOrder (Angle *a, Angle *b) : a(a), b(b) {}
};

bool inInterval (Angle *a, Angle *s, Angle *e);

bool intervalOverlap (Angle *s1, Angle *e1, Angle *s2, Angle *e2);

typedef pair<Angle *, Angle *> AngleInterval;

typedef vector<AngleInterval> AngleIntervals;

AngleIntervals intersect (const AngleIntervals &a, const AngleIntervals &b);

class Rot90Angle : public Angle {
  PTR<Angle> a;

  DeclareCalculate (PV2) {
    PV2<N> p = a->get<N>();
    return PV2<N>(- p.y, p.x);
  }
public:
  Rot90Angle (Angle *a) : a(a) {}
};

class BisectorAngle : public Angle {
  PTR<Angle> a, b;

  DeclareCalculate (PV2) {
    PV2<N> pa = a->get<N>(), pb = b->get<N>();
    N dc = pb.x - pa.x, ds = pb.y - pa.y, d = (dc*dc + ds*ds).sqrt();
    return PV2<N>(ds/d, - dc/d);
  }  
  public:
  BisectorAngle (Angle *ai, Angle *bi) : a(ai), b(bi) {
    if (b->order(a)) {
      a = bi; b = ai;
    }
  }
};

class ParametricAngle : public Angle {
  bool upper;
  PTR<Object<Scalar>> t;

  DeclareCalculate (PV2) {
    N s = t->get<N>(), d = 1.0 + s*s;
    return upper ? PV2<N>(-2.0*s/d, (1.0 - s*s)/d)
      : PV2<N>(2.0*s/d, (s*s - 1.0)/d);
  }
public:
  ParametricAngle (bool upper, double ti)
    : upper(upper), t(new Object<Scalar>(Scalar<double>(ti), false, true)) {}
  ParametricAngle (bool upper, Object<Scalar> *ti, VertexSet vs)
    : upper(upper), t(ti), Angle(vs) {}
  bool getUpper () const { return upper; }
  PTR<Object<Scalar>> getT () const { return t; }
};

class LinearPoly : public Object<Poly> {
  double k;

  DeclareCalculate (Poly) {
    Poly<N> p;
    p.a.push_back(- k);
    p.a.push_back(1);
    return p;
  }
public:
  LinearPoly (double k) : k(k) {}
};

class AngleParameter : public Object<Scalar> {
  Angle *a;

  DeclareCalculate (Scalar) {
    if (a == Angle::mpi)
      return N(-1);
    PV2<N> p = a->get<N>();
    return p.y.sign() == -1 ? p.x/(1.0 - p.y) : - p.x/(1.0 + p.y);
  }
public:
  AngleParameter (Angle *a) : a(a) {}
};

// k1*sin^2(x) + k2*cos^2(x) + k3*sin(x)*cos(x) + k4*sin(x) + k5*cos(x) + k6
template<class N>
class SinCosData {
 public:
  N k[6];
  bool linear;
  
  unsigned int size () const { return 6; }
  const N & operator[] (unsigned int i) const { return k[i]; }
  N & operator[] (unsigned int i) { return k[i]; }
  SinCosData () {}
  
  SinCosData (const N &k1, const N &k2, const N &k3) : linear(true) {
    k[0] = k1; k[1] = k2; k[2] = k3; k[3] = k[4] = k[5] = N(0);
  }
  
  SinCosData (const N &k1, const N &k2, const N &k3, const N &k4, const N &k5,
	      const N &k6) : linear(false) {
    k[0] = k1; k[1] = k2; k[2] = k3; k[3] = k4; k[4] = k5; k[5] = k6;
  }

  template<class M>
  SinCosData (const SinCosData<M> &x) : linear(x.linear) {
    for (int i = 0; i < 6; ++i) k[i] = x[i];
  }
};

class SinCosPoly : public Object<Poly> {
  PTR<Object<SinCosData>> x;
  bool flag;

  DeclareCalculate (Poly) {
    SinCosData<N> k = x->get<N>();
    Poly<N> p;
    double u = flag ? 1.0 : -1.0;
    if (k.linear) {
      p.a.push_back(k[0]*u + k[2]);
      p.a.push_back(-2.0*k[1]*u);
      p.a.push_back(k[2] - k[0]*u);
    }
    else {
      p.a.push_back(k[0] + k[3]*u + k[5]);
      p.a.push_back(-2.0*(k[2] + k[4]*u));
      p.a.push_back(2.0*(- k[0] + 2.0*k[1] + k[5]));
      p.a.push_back(2.0*(k[2] - k[4]*u));
      p.a.push_back(k[0] - k[3]*u + k[5]);
      p.checkDegree();
    }
    return p;
  }
public:
  SinCosPoly (Object<SinCosData> *x, bool flag) : x(x), flag(flag) {}
};

Angles sinCosAngles (Object<SinCosData> *x, Angle *s, Angle *e, const VertexSet &vs);

void sinCosAngles (Object<SinCosData> *x, bool upper, double l, double u,
		   const VertexSet &vs, Angles &as);

class PolyEF : public Object<SinCosData> {
  Edge *e;
  Face *f;
  bool aflag;

  DeclareCalculate (SinCosData) {
    PV3<N> u =  e->getU<N>(), n = f->getP()->get<N>().n;
    N k1 = aflag ? u.x*n.y - u.y*n.x : n.x*u.y - n.y*u.x,
      k2 = u.x*n.x + u.y*n.y, k3 = u.z*n.z;
    return SinCosData<N>(k1, k2, k3);
  }
public:
  PolyEF (Edge *e, Face *f, bool aflag)
     : e(e), f(f), aflag(aflag) {}
};

class PolyEFMin : public Primitive {
  Edge *e;
  Face *f;

  DeclareSign {
    PV3<N> u =  e->getU<N>(), n = f->getP()->get<N>().n;
    N k2 = u.x*n.x + u.y*n.y, k3 = u.z*n.z;
    return k3 - k2;
  }
public:
  PolyEFMin (Edge *e, Face *f) : e(e), f(f) {}
};

class EFKey {
public:
  Vertex *v[5];
  EFKey (Edge *e, Face *f) {
    v[0] = e->getT();
    v[1] = e->getH();
    sort(v, v + 2);
    HEdge *b = f->getBoundary();
    v[2] = b->tail();
    v[3] = b->head();
    v[4] = b->getNext()->head();
    sort(v + 2, v + 5);
  }

  EFKey (Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4, Vertex *v5) {
    v[0] = v1; v[1] = v2; v[2] = v3; v[3] = v4; v[4] = v5;
    sort(v, v + 2);
    sort(v + 2, v + 5);
  }

  bool operator< (const EFKey &x) const {
    for (int i = 0; i < 5; ++i)
      if (v[i] < x.v[i])
	return true;
      else if (x.v[i] < v[i])
	return false;
    return false;
  }
};  

class PolySFL : public Object<SinCosData> {
  Vertex *v[5];
  bool aflag;

  DeclareCalculate (SinCosData) {
    PV3<N> a[5];
    for (int i = 0; i < 5; ++i)
      a[i] = v[i]->getP()->get<N>();
    PV3<N> u = a[1] - a[0], n = (a[4] - a[3]).cross(a[2] - a[3]);
    N k1 = aflag ? u.x*n.y - u.y*n.x : n.x*u.y - n.y*u.x,
      k2 = u.x*n.x + u.y*n.y, k3 = u.z*n.z;
    return SinCosData<N>(k1, k2, k3);
  }
public:
  PolySFL (const EFKey &k, bool aflag) : aflag(aflag) {
    for (int i = 0; i < 5; ++i)
      v[i] = k.v[i];
  }
};

Angles anglesEF (Edge *e, Face *f, bool aflag);

Angles anglesSFL (const EFKey &k, bool aflag);

class CPoint : public Object<PV2> {
  friend class CFace;
  friend class CEdge;
  friend class Cspace;
protected:
  PTR<Angle> a;
public:
  CPoint (PTR<Angle> a) : a(a) {}
  PTR<Angle> getA () const { return a; }
  
  void getBBox (double *bbox) {
    PV2<Interval> xy = getApprox(1e-6);
    int s = AngleUpper(a) == 1 ? 1 : -1;
    Interval z = AngleParameter(a).getApprox(1e-6).x;
    bbox[0] = xy.x.l;
    bbox[1] = xy.x.u;
    bbox[2] = xy.y.l;
    bbox[3] = xy.y.u;
    bbox[4] = z.l + s;
    bbox[5] = z.u + s;
  }
};

typedef vector<PTR<CPoint>> CPoints;

class COrder :public Primitive {
  CPoint *a, *t, *h;

  DeclareSign {
    return (a->get<N>() - t->get<N>()).dot(h->get<N>() - t->get<N>());
  }
public:
  COrder (CPoint *a, CPoint *t, CPoint *h) : a(a), t(t), h(h) {}
};

bool onEdge (CPoint *a, CPoint *t, CPoint *h);

class CPointOrderX :public Primitive {
  CPoint *a, *b;

  DeclareSign { return b->get<N>().x - a->get<N>().x; }
public:
  CPointOrderX (CPoint *a, CPoint *b) : a(a), b(b) {}
};

class CPointOrderXO {
public:
  bool operator() (CPoint *a, CPoint *b) const {
    return a != b && CPointOrderX(a, b) == 1;
  }
};

class CEdge;

template<class N>
class SpiralData {
public:
  PV2<N> a, b;
  unsigned int size () const { return 4; }
  const N & operator[] (unsigned int i) const { return i < 2 ? a[i] : b[i-2]; }
  N & operator[] (unsigned int i) { return i < 2 ? a[i] : b[i-2]; }
  SpiralData () {}
  SpiralData (const PV2<N> &a, const PV2<N> &b) : a(a), b(b) {}
  template<class M>
    SpiralData (const SpiralData<M> &p) : a(p.a), b(p.b) {}
};

class Spiral : public Object<SpiralData> {
public:
  Vertex *v;
  Edge *e;
  bool aflag;
  VertexSet vs;
  vector<CEdge *> edges;

  Spiral (Vertex *v, Edge *e, bool aflag) : v(v), e(e), aflag(aflag) {
    vs.insert(v);
    vs.insert(e->getT());
    vs.insert(e->getH());
  }

  void vertices (VertexSet &vsa, VertexSet &vsb) const {
    if (aflag) {
      vsa.insert(v);
      vsb.insert(e->getT());
      vsb.insert(e->getH());
    }
    else {
      vsa.insert(e->getT());
      vsa.insert(e->getH());
      vsb.insert(v);
    }
  }

  int sharedVertices (Angle *a) const { return intersection(vs, a->getVS()).size(); }

  void bboxTP (Angle *as, Angle *ae, double *bbox);

  template<class N>
  PV2<N> xy (Angle *a) {
    SpiralData<N> sd = get<N>();
    return sd.a.rotate(a->get<N>()) + sd.b;
  }

  PV2<double> xyApprox (Angle *a) {
    SpiralData<double> sd = getApproxMid();
    PV2<double> aa = a->getApproxMid();
    return sd.a.rotate(aa) + sd.b;
  }

  DeclareCalculate (SpiralData) {
    PV3<N> p = v->getP()->get<N>(), t = e->getT()->getP()->get<N>(),
    h = e->getH()->getP()->get<N>(), u = h - t, w = t.cross(u);
    if (aflag)
      return SpiralData<N>(PV2<N>(- p.x, - p.y),
			   PV2<N>((p.z*u.x - w.y)/u.z, (p.z*u.y + w.x)/u.z));
    return SpiralData<N>(PV2<N>((w.y - p.z*u.x)/u.z, - (w.x + p.z*u.y)/u.z),
			 PV2<N>(p.x, p.y));
  }
};

class SpiralXMaxAngle : public Angle {
  Spiral *s;

  DeclareCalculate (PV2) {
    PV2<N> a = s->get<N>().a;
    return PV2<N>(a.x, - a.y);
  }
public:
  SpiralXMaxAngle (Spiral *s) : s(s) {}
};

class SpiralOrder {
public:
  bool operator() (Spiral *a, Spiral *b) const {
    if (a == b) return false;
    return a->v < b->v || a->v == b->v && a->e < b->e;
  }
};

typedef set<Spiral *, SpiralOrder> SpiralSet;

class SpiralPoint : public CPoint {
  friend class CEdge;
  Spiral *s;
  DeclareCalculate(PV2) { return s->xy<N>(a); }
public:
  SpiralPoint (PTR<Angle> a, Spiral *s) : CPoint(a), s(s) {}
  Spiral * getS () const { return s; }
};

class CEdge;

class HCEdge;

class CVertex {
  friend class CEdge;
  friend class CFace;
  friend class Cspace;
protected:
  PTR<CPoint> p;
  vector<CEdge *> edges;
  double bbox[6];
public:
  CVertex (CPoint *p) : p(p) { p->getBBox(bbox); }
  CPoint * getP () { return p; }
  int EdgesN () const { return edges.size(); }
  CEdge * getEdge (int i) const { return edges[i]; }
  double * getBBox () { return bbox; }
  void outgoingHEdges (vector<HCEdge *> &ed) const;
};

typedef vector<CVertex *> CVertices;

typedef set<CVertex *> CVertexSet;

CVertexSet intersection (const CVertexSet &a, const CVertexSet &b);

class VertexAngleOrder {
public:
  bool operator() (CVertex *a, CVertex *b) const {
    return a != b && a->getP()->getA()->order(b->getP()->getA());
  }
};

class CFace;

class CEdge {
  friend class CPoint;
  friend class CVertex;
  friend class HCEdge;
  friend class CFace;
  friend class Cspace;
protected:
  CVertex *t, *h;
  bool flag, inc, dec;
  Spiral *s;
  double bbox[6];
  vector<HCEdge *> hedges;
  CVertices vertices;
  map<Angle *, PTR<CPoint>> apmap;
public:
  CEdge (CVertex *t, CVertex *h, bool flag);
  ~CEdge ();
  CVertex * getT () const { return t; }
  CVertex * getH () const { return h; }
  Spiral * getSpiral () const { return s; }
  int HEdgesN () const { return hedges.size(); }
  HCEdge * getHEdge (int i) const { return hedges[i]; }
  double * getBBox () { return bbox; }
  void addVertex (CVertex *v);
  HCEdge * addHEdge (bool forward);
  void removeHEdge (HCEdge *e);
  void edgeVertices (CVertices &ve) const;
  CFace * otherFace (CFace *f) const;
  bool piEdge () const;
  bool horizontal () const { return t->p->a == h->p->a; }
  bool increasing () const { return inc; }
  bool decreasing () const { return dec; }
  AngleInterval angleInterval () const;
  bool contains (Angle *a) const { return inInterval(a, t->p->a, h->p->a); }
  bool contains (CPoint *p) const { return onEdge(p, t->p, h->p); }
  CPoint * point (Angle *a);
  CPoint * pointAux (Angle *a) const;
  template<class N> PV3<N> getU (CPoint *p) const;
  void sortHEdges ();
};

typedef vector<CEdge *> CEdges;

typedef set<CEdge *> CEdgeSet;

class HCEdge {
  friend class CEdge;
  friend class CFace;
  friend class Cspace;
protected:
  CEdge *e;
  bool forward, flag;
  CFace *f;
  HCEdge *next, *hp;
public:
  HCEdge (CEdge *e, bool forward) 
    : e(e), forward(forward), flag(false), next(0), hp(0), f(0) {}
  CVertex * tail () const { return forward ? e->t : e->h; }
  CVertex * head () const { return forward ? e->h : e->t; }
  CEdge * getE () const { return e; }
  bool getForward () const { return forward; }
  CFace * getF () const { return f; }
  HCEdge * getNext () const { return next; }
  void setNext (HCEdge *h) { next = h; }
  HCEdge * getHP () const { return hp; }
  HCEdge * cw () const;
  HCEdge * ccw () const;
  void loop (vector<HCEdge *> &ed);
  bool onLoop ();
  bool increasing () const { return forward ? e->increasing() : e->decreasing(); }
  bool decreasing () const { return forward ? e->decreasing() : e->increasing(); }
  template<class N> PV3<N> getU (CPoint *p) const;
  template<class N> PV3<N> getN (CPoint *p) const;
};

typedef vector<HCEdge *> HCEdges;

typedef set<HCEdge *> HCEdgeSet;

class CEdgeOrder1 : public Primitive {
  HCEdge *e;
  CPoint *p;

  DeclareSign { return Rdir->get<N>().dot(e->getN<N>(p)); }
public:
  CEdgeOrder1 (HCEdge *e, CPoint *p) : e(e), p(p) {}
};

class CEdgeOrder2 : public Primitive {
  CEdge *e;
  HCEdge *f, *g;
  CPoint *p;

  DeclareSign { return e->getU<N>(p).tripleProduct(g->getN<N>(p), f->getN<N>(p)); }
public:
  CEdgeOrder2 (CEdge *e, HCEdge *f, HCEdge *g, CPoint *p) : e(e), f(f), g(g), p(p) {}
};

class CEdgeOrder {
  CEdge *e;
public:
  CEdgeOrder (CEdge *e) : e(e) {}
  bool operator() (HCEdge *f, HCEdge *g) const {
    if (f == g) return false;
    CPoint *p = e->getT()->getP();
    int sf = CEdgeOrder1(f, p), sg = CEdgeOrder1(g, p);
    return sf != sg ? sf == 1 : CEdgeOrder2(e, f, g, p) == 1;
  }
};

enum HHCEdgeType { RightEnd, LeftEnd, LeftStart, RightStart};

class HHCEdge {
public:
  HCEdge *e;
  bool f;
  HHCEdge (HCEdge *e, bool f) : e(e), f(f) {}
  CVertex * tail () const { return f ? e->tail() : e->head(); }
  CVertex * head () const { return f ? e->head() : e->tail(); }
  
  HHCEdgeType type () const {
    if (f)
      return e->increasing() ? RightStart : LeftEnd;
    return e->increasing() ? RightEnd : LeftStart;
  }
  
  bool operator< (const HHCEdge &x) const {
    Angle *a = tail()->getP()->getA(), *xa = x.tail()->getP()->getA();
    return a == xa ? type() < x.type() : a->order(xa);
  }
};

typedef vector<HHCEdge> HHCEdges;

class Slab {
  double bbox[6];
public:
  CEdge *l, *r;
  PTR<Angle> s, e;
  CPoints pb, pt;
  CFace *f;

  Slab (CEdge *l, CEdge *r, Angle *s, Angle *e, CFace *f)
    : l(l), r(r), s(s), e(e), f(f) { bbox[0] = 1.0; bbox[1] = 0.0; }
  double * getBBox ();
  void addPoint (CPoint *p);
};

typedef vector<Slab *> Slabs;

typedef set<Slab *> SlabSet;

bool inSlab (CPoint *p, Slab *s);

class InSlab : public Primitive {
  CPoint *p;
  Slab *s;

  DeclareSign {
    Angle *a = p->getA();
    PTR<CPoint> pl = s->l->point(a), pr = s->r->point(a);
    PV2<N> ppl = pl->get<N>(), ppr = pr->get<N>(), u = ppr - ppl;
    N k = u.dot(p->get<N>() - ppl);
    return k*(u.dot(u) - k);
  }
public:
  InSlab (CPoint *p, Slab *s) : p(p), s(s) {}
};

class CShell;

enum CFaceType { CFaceVF, CFaceFV, CFaceEEP, CFaceEEM };

class CFace {
  friend class HCEdge;
  friend class CShell;
  friend class Cspace;
protected:
  CFaceType type;
  int *feature1, *feature2;
  HCEdges boundary;
  double bbox[6];
  CEdges edges;
  VertexSet vs;
  CShell *s;
  Slabs slabs;
public:
  CFace (Vertex *v, Face *f, bool aflag);
  CFace (Edge *e, Edge *f, bool pflag);
  CFace (CFaceType type, int *feature1, int *feature2)
    : type(type), feature1(feature1), feature2(feature2), s(0) {}
  ~CFace () {
    for (Slabs::iterator s = slabs.begin(); s != slabs.end(); ++s)
      delete *s;
  }
  CFaceType getType () const { return type; }
  int * getFeature1 () const { return feature1; }
  int * getFeature2 () const { return feature2; }
  HCEdge * getBoundary (int i) const { return boundary[i]; }
  const HCEdges &getBoundary () { return boundary; } 
  double * getBBox () { return bbox; }
  void getBBox (double *b) const { copyBBox(bbox, b); }
  const CEdges & getEdges () const { return edges; }
  const Slabs & getSlabs () const { return slabs; }
  Angle * startAngle () const { return boundary[0]->tail()->p->a; }
  Angle * endAngle () const { return boundary[0]->next->head()->p->a; }
  CEdge * leftEdge () const { return boundary[0]->next->next->next->e; }
  CEdge * rightEdge () const { return boundary[0]->next->e; }
  void addBoundary (HCEdge *e);
  void vertices (VertexSet &vsa, VertexSet &vsb) const;
  Spiral * sharedSpiral (CFace *f) const;
  void boundaryVertices (CVertices &ve) const;
  bool boundaryEdge (CEdge *e) const;
  void boundaryHEdges (HCEdges &ed) const;
  void sharedBoundaryVertices (CFace *f, CVertices &vfg) const;
  void edgeVertices (CVertexSet &vs) const;
  int sharedVertices (Angle *a) const { return intersection(vs, a->getVS()).size(); }
  bool sameBoundary (CVertex *a, CVertex *b) const;
  bool sharedBoundaryVertex (CFace *f, CFace *g) const;
  bool bboxOverlap (CFace *f, CFace *g) const;
  bool angleOverlap (CFace *f, CFace *g) const;
  bool inInterval (Angle *a) const;
  bool contains (CPoint *p) const;
  bool contains2(CPoint *p) const;
  template<class N> PV3<N> getN (CPoint *p) const;
  template<class N> void coeffs (PV2<N> &a1, PV2<N> &a2, N &k1, N &k2, N &k3) const;
  template<class N> void coeffsVF (PV2<N> &a1, PV2<N> &a2, N &k1, N &k2, N &k3) const;
  template<class N> void coeffsFV (PV2<N> &a1, PV2<N> &a2, N &k1, N &k2, N &k3) const;
  template<class N> void coeffsEE (PV2<N> &a1, PV2<N> &a2, N &k1, N &k2, N &k3) const;
  template<class N> void line (Angle *a, N &u, N &v, N &w) const;
  CFace * neighbor (HCEdge *e) const;
  void formSlabs();
  void initSlabs (HHCEdges &hhe, CPoints &hpts) const;
  void updateRE (Slabs &sl, const HHCEdge &h);
  void updateLE (Slabs &sl, const HHCEdge &h);
  void updateLS (Slabs &sl, const HHCEdge &h);
  void updateRS (Slabs &sl, const HHCEdge &h);
  CEdge * slabNextIncreasing (Slab *s) const;
  void splitSlab (Slab *s, const Angles &an, Slabs &sl);
};

typedef vector<CFace *> CFaces;

class FFPoint : public CPoint {
  CFace *f1, *f2;
  DeclareCalculate(PV2) {
    N a1, b1, c1, a2, b2, c2;
    f1->line(a, a1, b1, c1);
    f2->line(a, a2, b2, c2);
    N den = a1*b2 - a2*b1, x = (b1*c2 - b2*c1)/den, y = (a2*c1 - a1*c2)/den;
    return PV2<N>(x, y);
  }
public:
  FFPoint (Angle *a, CFace *f1, CFace *f2)
    : CPoint(a), f1(f1), f2(f2) {}
  CFace * getF1 () const { return f1; }
  CFace * getF2 () const { return f2; }
};    

class IndependentFaces :public Primitive {
  CFace *f, *g;
  Angle *a;

  DeclareSign {
    N a1, b1, c1, a2, b2, c2;
    f->line<N>(a, a1, b1, c1);
    g->line<N>(a, a2, b2, c2);
    return a1*b2 - a2*b1;
  }
public:
  IndependentFaces (CFace *f, CFace *g, Angle *a) : f(f), g(g), a(a) {}
};

class EdgeOrderL : public Primitive {
  CPoint *t, *h, *v, *w;

  DeclareSign {
    PV2<N> u = h->get<N>() - t->get<N>(), vw = w->get<N>() - v->get<N>();
    return u.dot(vw);
  }
public:
  EdgeOrderL (CPoint *t, CPoint *h, CPoint *v, CPoint *w)
    : t(t), h(h), v(v), w(w) {}
};

class EdgePointOrderL {
public:
  EdgePointOrderL (CPoint *t, CPoint *h) : t(t), h(h) {}
  bool operator() (CPoint *v, CPoint *w) const {
    return v != w && EdgeOrderL(t, h, v, w) == 1;
  }

  CPoint *t, *h;
};

class EdgeVertexOrderL {
public:
  EdgeVertexOrderL (CEdge *e) : t(e->getT()), h(e->getH()) {}
  bool operator() (CVertex *v, CVertex *w) const {
    return v != w && EdgeOrderL(t->getP(), h->getP(), v->getP(), w->getP()) == 1;
  }

  CVertex *t, *h;
};

typedef pair<CEdge *, CEdge *> CEEPair;

typedef set<CEEPair> CEEPairSet;

typedef pair<EFKey, Angles> EFAPair;

typedef map<EFKey, Angles> EFAMap;

typedef pair<Angle *, Spiral *> ASPair;

typedef pair<ASPair, CVertex *> ASVPair;

typedef map<ASPair, CVertex *> ASVMap;

typedef pair<CEdge *, CFace *> CEFPair;

typedef pair<CEFPair, CVertices> CEFVPair;

typedef map<CEFPair, CVertices> CEFVMap;

typedef pair<Spiral *, CFace *> SFPair;

typedef pair<Angle *, CFace *> AFPair;

typedef pair<AFPair, CVertex *> AFVPair;

typedef map<AFPair, CVertex *> AFVMap;

typedef pair<SFPair, Angles> SFAPair;

typedef map<SFPair, Angles> SFAMap;

class LinePatchPoint : public CPoint {
protected:
  CEdge *e;
  CFace *f;
  DeclareCalculate(PV2) {
    CPoint *et = e->getT()->getP(), *eh = e->getH()->getP();
    Spiral *l = f->leftEdge()->getSpiral(), *r = f->rightEdge()->getSpiral();
    PV2<N> pl = l->xy<N>(a), v = r->xy<N>(a) - pl, u = eh->get<N>() - et->get<N>(),
      w = et->get<N>() - pl;
    N k = u.cross(w)/u.cross(v);
    return pl + k*v;    
  }
public:
  LinePatchPoint (CEdge *e, CFace *f)
    : CPoint(e->getT()->getP()->getA()), e(e), f(f) {}
};

typedef pair<CFace *, CFace *> CFFPair;

typedef pair<CFFPair, CEdges> CFFEPair;

typedef map<CFFPair, CEdges> CFFEMap;

class PolySF : public Object<SinCosData> {
protected:
  Spiral *s;
  CFace *f;

  DeclareCalculate (SinCosData) {
    SpiralData<N> sd1 = s->get<N>(), sd2 = f->leftEdge()->getSpiral()->get<N>(),
      sd3 = f->rightEdge()->getSpiral()->get<N>();
    PV2<N> a12 = sd1.a - sd2.a, b12 = sd1.b - sd2.b, a23 = sd2.a - sd3.a,
      b23 = sd2.b - sd3.b;
    N k1 = a23.dot(b12) - a12.dot(b23),  k2 = a12.cross(b23) - a23.cross(b12),
      k3 = a12.cross(a23) + b12.cross(b23);
    return SinCosData<N>(k1, k2, k3);
  }
public:
  PolySF (Spiral *s, CFace *f) : s(s), f(f) {}
};

Angles anglesSF (Spiral *s, CFace *f);

class PolyFFF : public Object<SinCosData> {
  CFace *f1, *f2, *f3;
    
  DeclareCalculate (SinCosData) {
    N k1, k2, k3, k4, k5, k6;
    coeffs(k1, k2, k3, k4, k5, k6);
    return SinCosData<N>(k1, k2, k3, k4, k5, k6);
  };

  template<class N>
  void coeffs (N &k1, N &k2, N &k3, N &k4, N &k5, N &k6) {
    PV2<N> a11, a21, a12, a22, a13, a23;
    N k11, k21, k31, k12, k22, k32, k13, k23, k33;
    f1->coeffs<N>(a11, a21, k11, k21, k31);
    f2->coeffs<N>(a12, a22, k12, k22, k32);
    f3->coeffs<N>(a13, a23, k13, k23, k33);
    k1 = k2 = k3 = k4 = k5 = k6 = N(0);
    coeffs<N>(k11, k21, k31, a12, a22, a13, a23, k1, k2, k3, k4, k5, k6);
    coeffs<N>(k12, k22, k32, a13, a23, a11, a21, k1, k2, k3, k4, k5, k6);
    coeffs<N>(k13, k23, k33, a11, a21, a12, a22, k1, k2, k3, k4, k5, k6);
  }

  template<class N>
  void coeffs (N &k1, N &k2, N &k3, PV2<N> &a1, PV2<N> &b1, PV2<N> &a2,
	       PV2<N> &b2, N &x1, N &x2, N &x3, N &x4, N &x5, N &x6) {
    N k4 = a2.dot(b1) - a1.dot(b2), k5 = a1.cross(b2) - a2.cross(b1),
      k6 = a1.cross(a2) + b1.cross(b2);
    x1 = x1 + k1*k4;
    x2 = x2 + k2*k5;
    x3 = x3 + k1*k5 + k2*k4;
    x4 = x4 + k1*k6 + k3*k4;
    x5 = x5 + k2*k6 + k3*k5;
    x6 = x6 + k3*k6;
  }

public:
  PolyFFF (CFace *f1, CFace *f2, CFace *f3) : f1(f1), f2(f2), f3(f3) {}
};

Angles anglesFFF (CFace *f1, CFace *f2, CFace *f3);

typedef pair<CVertex *, CVertex *> CVVPair;

typedef map<CVertex *, CVertex *> CVVMap;

class CCell;

class CShell {
public:
  CFaces faces;
  CVertex *vm;
  bool outer;
  CCell *c;

  CShell () : vm(0), outer(false), c(0) {}
  void init ();
};

typedef vector<CShell *> CShells;
  
class CShellOrder {
public:
  bool operator() (CShell *a, CShell *b) const {
    return a != b && CPointOrderX(a->vm->getP(), b->vm->getP()) == 1;
  }
};

class CCell {
public:
  CShells shells;
  ~CCell ();
};

typedef vector<CCell *> CCells;

class Cspace : public RefCnt {
public:
  CVertices vertices;
  CEdges edges;
  CFaces faces;
  CCells cells;
  double bbox[6];
  EFAMap efamap;
  SpiralSet ss;
  ASVMap asvmap;
  CEFVMap efvmap;
  AFVMap afvmap;
  SFAMap sfamap;
  Cspace *par;

  Cspace (Cspace *par = 0) : par(par) {}
  ~Cspace ();
  Angles angles (Edge *e, Face *f, bool aflag);
  Angles angles (const EFKey &k, bool aflag);
  void angleIntervals (Edge *e, Face *f, bool aflag, AngleIntervals &pos,
		       AngleIntervals &neg);
  void angleIntervals (HEdge *e, Face *f, bool aflag, AngleIntervals &pos,
		       AngleIntervals &neg);
  Spiral * getSpiral (Vertex *v, Edge *e, bool aflag);
  CVertex * getVertex (CPoint *p);
  CVertex * getVertex (Angle *a, Spiral *s);
  HCEdge * addHEdge (CVertex *a, CVertex *b, bool flag, HCEdge *hp = 0);
  CEdge * getEdge (CVertex *a, CVertex *b, bool flag, HCEdge *hp = 0);
  HCEdge * addLoop (const CVertices &ve);
  void patches (Polyhedron *a, Polyhedron *b);
  void patchesVF (Polyhedron *a, Polyhedron *b, bool aflag);
  void patchVF (Vertex *v, Face *f, bool aflag);
  void patchVF (Vertex *v, Face *f, bool aflag, Angle *s, Angle *e);
  void patch (Angle *s, Angle *e, Spiral *l, Spiral *r, CFace *f);
  void patchesEE (Polyhedron *a, Polyhedron *b);
  void patchEE (Edge *e, Edge *f);
  void patchEE (Edge *e1, Edge *e2, Angle *s, Angle *e, bool flag);
  void intersectFF ();
  void intersectEE (CFace *f, CFace *g, CEEPairSet &eps);
  void intersectEE (CEdge *e, CEdge *f);
  void intersectSS (CEdge *e, CEdge *f) const;
  void intersectSL (CEdge *e, CEdge *f);
  void intersectLL (CEdge *e, CEdge *f);
  void intersectFF (CFace *f, CFace *g);
  void intersectFF (CFace *f, CFace *g, CVertices &vfg);
  void intersectEF (HCEdge *e, CFace *f, CVertices &vfg);
  void intersectEF (CEdge *e, CFace *f, CVertices &ve);
  CVertex * intersectEFLine (CEdge *e, CFace *f);
  Angles intersectSF (Spiral *s, CFace *f);
  Angles intersectSFLine (Spiral *s, CFace *f);
  void formFF (CFace *f, CFace *g, CVertices &vfg);
  bool intersectsFFLine (CFace *f, CFace *g) const;
  void formFFLine (CFace *f, CFace *g, CVertices &vfg);
  CVertex * spiralVertex (Angle *a, Spiral *s);
  void formFFCurve (CFace *f, CFace *g, CVertices &vfg);
  void formFF (CFace *f, CFace *g, CVertex *v, CVertex *w);
  void intersectFFF ();
  void formFFEMap (CFFEMap &ffemap) const;
  void intersectFFF (CFace *f, const CFFEMap &ffemap);
  void intersectFFF (CFace *f1, CFace *f2, CFace *f3, CEdge *e12, CEdge *e13,
		     const CFFEMap &ffemap);
  void intersectFFFH (CFace *f1, CFace *f2, CFace *f3, CEdge *eh,
		      CEdge *e1, CEdge *e2);
  void intersectFFFG (CFace *f1, CFace *f2, CFace *f3, CEdge *e12, CEdge *e13,
		      const CEdges &ed);
  void sortVertices ();
  Cspace * subfaces ();
  void subfaces (CFace *f, Cspace *a, CVVMap &vvmap) const;
  void subedges (CFace *f, Cspace *a, CVVMap &vvmap, HCEdges &he) const;
  void subedges (HCEdge *e, Cspace *a, CVVMap &vvmap, HCEdges &he) const;
  CVertex * getVertex (CVertex *v, CVVMap &vmap);
  void setNext (CFace *f, const HCEdges &he) const;
  CFace * addFace (HCEdge *e, CFace *f);
  void removeBad (const HCEdges &ed);
  void formCells (Polyhedron *a, Polyhedron *b);
  void formShells (Polyhedron *a, Polyhedron *b, CShells &sh);
  void formShells (CShells &sh) const;
  CShell * enclosingShell (CShell *s, Octree<CFace *> *octree);
  void removeFace (CFace *f) const;
  void describe () const;
};

class PointOrderZ :public Primitive {
  Point *a, *b;

  DeclareSign { return b->get<N>().z - a->get<N>().z; }
public:
  PointOrderZ (Point *a, Point *b) : a(a), b(b) {}
};

bool inIntervalZ (Vertex *v, Edge *e);

class OrientationVF : public Primitive {
  bool aflag;
  Face *f;
  Angle *a;
  Spiral *l, *r;

  DeclareSign {
    PV3<N> n3 = aflag ? f->getP()->get<N>().n :
      - f->getP()->get<N>().n.rotateZ(a->get<N>());
    PV2<N> n(n3.x, n3.y), t = r->xy<N>(a) - l->xy<N>(a);
    return n.cross(t);
  }
public:
  OrientationVF (bool aflag, Face *f, Angle *a, Spiral *l, Spiral *r)
    : aflag(aflag), f(f), a(a), l(l), r(r) {}
};

class OrientationEE : public Primitive {
  bool flag;
  Edge *e, *f;
  Angle *a;
  Spiral *l, *r;

  DeclareSign {
    PV3<N> n3 = e->getU<N>().rotateZ(a->get<N>()).cross(f->getU<N>());
    PV2<N> n(n3.x, n3.y), t = r->xy<N>(a) - l->xy<N>(a);
    return flag ? t.cross(n) : n.cross(t);
  }
public:
  OrientationEE (bool flag, Edge *e, Edge *f, Angle *a, Spiral *l, Spiral *r)
    : flag(flag), e(e), f(f), a(a), l(l), r(r) {}
};

class CFFOrder :public Primitive {
  CFace *f, *g;
  CVertex *v, *w;

  DeclareSign {
    CPoint *vp = v->getP(), *wp = w->getP();
    Angle *va = vp->getA(), *wa = wp->getA();
    PV3<N> u = f->getN<N>(vp).cross(g->getN<N>(vp));
    if (va == wa) {
      PV2<N> vw = vp->get<N>() - wp->get<N>(), u2(u.x, u.y);
      return u2.dot(vw);
    }
    return - u.z;
  }
public:
  CFFOrder (CFace *f, CFace *g, CVertex *v, CVertex *w)
    : f(f), g(g), v(v), w(w) {}
};

class EdgeVertexOrderS {
public:
  EdgeVertexOrderS (CEdge *e) : inc(e->increasing()) {}
  bool operator() (CVertex *v, CVertex *w) const {
    return v != w && inc == v->getP()->getA()->order(w->getP()->getA());
  }
  bool inc;
};

class CVertexHHEdgeOrder1 :public Primitive {
  CVertex *v;
  HCEdge *e;

  DeclareSign {
    N k = Rdir->get<N>().dot(e->getU<N>(v->getP()));
    return v == e->tail() ? k : - k;
  }
public:
  CVertexHHEdgeOrder1 (CVertex *v, HCEdge *e) : v(v), e(e) {}
};

class CVertexHHEdgeOrder2 :public Primitive {
  CVertex *v;
  HCEdge *e, *f;
  CFace *g;

  DeclareSign {
    CPoint *p = v->getP();
    N k = g->getN<N>(p).tripleProduct(f->getU<N>(p), e->getU<N>(p));
    return (v == e->tail()) == (v == f->tail()) ? k : - k;
  }
public:
  CVertexHHEdgeOrder2 (CVertex *v, HCEdge *e, HCEdge *f, CFace *g)
    : v(v), e(e), f(f), g(g) {}
};  

class HHCEdgeOrder {
  CFace *g;
public:
  HHCEdgeOrder (CFace *g) : g(g) {}
  bool operator() (const HHCEdge &e, const HHCEdge &f) const {
    CVertex *u = e.tail(), *v = f.tail();
    if (u != v)
      return u < v;
    if (e.e == f.e)
      return false;
    int se = CVertexHHEdgeOrder1(u, e.e), sf = CVertexHHEdgeOrder1(u, f.e);
    return se != sf ? se == 1 : CVertexHHEdgeOrder2(u, e.e, f.e, g) == 1;
  }
};

class RayPatchPointX : public CPoint {
protected:
  CPoint *p;
  CFace *f;
  DeclareCalculate(PV2) {
    HCEdge *e = f->getBoundary(0);
    CFace *g = e->getHP()->getF();

    Spiral *l = g->leftEdge()->getSpiral(), *r = g->rightEdge()->getSpiral();
    PV2<N> pl = l->xy<N>(a), v = r->xy<N>(a) - pl, w = p->get<N>() - pl;
    N k = w.y/v.y;
    return pl + k*v;
  }
public:
  RayPatchPointX (CPoint *p, CFace *f) : CPoint(p->getA()), p(p), f(f) {}
};

Cspace * cspace (Polyhedron *a, Polyhedron *b);

class TransformedPoint : public Point {
  CPoint *c;
  Point *a;
  
  DeclareCalculate(PV3) {
    PV3<N> p = a->get<N>();
    PV2<N> q = PV2<N>(p.x, p.y).rotate(c->getA()->get<N>()) + c->get<N>();
    return PV3<N>(q.x, q.y, p.z);
  }

public:
  TransformedPoint (CPoint *c, Point *a) : c(c), a(a) {}
};

Polyhedron * transform (Polyhedron *p, CPoint *c);

typedef pair<CPoint *, Vertex *> CPVPair;

typedef map<CPoint *, Vertex *> CPVMap;

typedef pair<CEdge *, CPoints> EPPair;

typedef map<CEdge *, CPoints> EPMap;

class CTriangle {
public:
  PTR<CPoint> a, b, c;
  double bbox[6];
  
  CTriangle (CPoint *a, CPoint *b, CPoint *c) : a(a), b(b), c(c) {
    bbox[0] = 1.0; bbox[1] = 0.0;
  }

  CPoint * point (int i) const {
    if (i == 0) return a;
    return i == 1 ? b : c;
  }
};

typedef vector<CTriangle> CTriangles;

typedef pair<Slab *, CTriangles> STPair;

typedef map<Slab *, CTriangles> STMap;

class DistancePL : public Primitive {
  Point *a, *t, *h;
  double d;

  DeclareSign {
    PV3<N> tp = t->get<N>(), u = h->get<N>() - tp,
      p = tp + ((a->get<N>() - tp).dot(u)/u.dot(u))*u, w = a->get<N>() - p;
    return d*d - w.dot(w);
  }
 public:
  DistancePL (Point *a, Point *t, Point *h, double d) : a(a), t(t), h(h), d(d) {}
};

Polyhedron * discretize (Cspace *c, double d);

Polyhedron * discretize (const CFaces &fa, double d, CPVMap &pvmap);

EPMap discretizeEdges (const CFaces &fa, double d);

CPoints discretize (CEdge *e, double d);

void discretize (CEdge *e, Angle *as, Angle *ae, double d, CPoints &pts);

bool close (CPoint *a, CPoint *t, CPoint *h, double d);

void delentil (const CFaces &fa, EPMap &epmap);

void delentil (Slab *s, EPMap &epmap);

void addPoint (CEdge *e, Angle *a, EPMap &epmap);

CTriangles discretize (Slab *s, const EPMap &epmap);

CPoints getPoints (CEdge *f, Angle *s, Angle *e, const EPMap &epmap);

CTriangles discretize (const CPoints &l, const CPoints &r, Slab *s);

void ctriangle (CPoint *a, CPoint *b, CPoint *c, CTriangles &tr);

void splitB (const CTriangle &t, const CPoints &pb, bool lflag, CTriangles &tr);

void splitT (const CTriangle &t, const CPoints &pt, bool lflag, CTriangles &tr);

void addTriangles (const STMap &stmap, CPVMap &pvmap, Polyhedron *a);

Vertex * getVertex (CPoint *p, CPVMap &cpvmap, Polyhedron *a);

class PointCPoint : public Point {
  PTR<CPoint> p;
  DeclareCalculate(PV3) {
    PV2<N> xy = p->get<N>();
    Angle *a = p->getA();
    if (a == Angle::mpi)
      return PV3<N>(xy.x, xy.y, N(-2));
    PV2<N> q = a->get<N>();
    N z = q.y.sign() == -1 ? q.x/(1.0 - q.y) - 1 : - q.x/(1.0 + q.y) + 1;
    return PV3<N>(xy.x, xy.y, z);
  }

public:
  PointCPoint (CPoint *p) : p(p) {}
  CPoint * getCP () const { return p; }
};

void facesPI (Cspace *c, Polyhedron *a, CPVMap &pvmap, bool flag);

Points * loopPI (Polyhedron *a, CPVMap &cpvmap, CEdgeSet &es, CEdge *e,
		 PVMap &pvmap);

#endif
