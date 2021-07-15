//  -*- Mode: c++ -*-

#include "cspace.h"
#include"io.h"
#include "polym.h"

template<class N>
void rotateParametric (const PV2<N> &a, const PV2<N> &b, int u, Poly<N> &x,
		       Poly<N> &y)
{
  x.add(u*b.x + a.y, 0);
  x.add(2*a.x, 1);
  x.add(u*b.x - a.y, 2);
  y.add(u*b.y - a.x, 0);
  y.add(2*a.y, 1);
  y.add(u*b.y + a.x, 2);
}

template<class N>
PV3<Poly<N>> pointParametric (Spiral *s, int u, Poly<N> &w)
{
  SpiralData<N> sd = s->get<N>();
  Poly<N> z(- u);
  z.add(1, 1);
  PV3<Poly<N>> p;
  rotateParametric(sd.a, sd.b, u, p.x, p.y);
  w.add(u, 0);
  w.add(u, 2);
  p.z = z*w;
  return p;
}

template<class N>
PV3<Poly<N>> lineParametric (CFace *f, int u)
{
  PV2<N> a, b;
  N k1, k2, k3;
  f->coeffs<N>(a, b, k1, k2, k3);
  PV3<Poly<N>> p;
  rotateParametric(a, b, u, p.x, p.y);
  p.z.add(u*k3 - k1, 0);
  p.z.add(2*k2, 1);
  p.z.add(u*k3 + k1, 2);
  return p;
}

template<class N>
PV3<Poly<N>> pointParametric (CFace *f1, CFace *f2, int u, Poly<N> &w)
{
  PV3<Poly<N>> l1 = lineParametric<N>(f1, u), l2 = lineParametric<N>(f2, u);
  w = l1.x*l2.y - l1.y*l2.x;
  Poly<N> z(- u);
  z.add(1, 1);
  return PV3<Poly<N>>(l1.y*l2.z - l1.z*l2.y, l1.z*l2.x - l1.x*l2.z, z*w);
}

class CPointRecord {
public:
  Spiral *s;
  CFace *f1, *f2;
  ParametricAngle *a;

  CPointRecord () : s(0), f1(0), f2(0), a(0) {}
  
  template<class N>
  PV3<Poly<N>> point (Poly<N> &w) {
    int u = a->getUpper() ? -1 : 1;
    return s ? pointParametric(s, u, w) : pointParametric(f1, f2, u, w);
  }

  template<class N>
  PV3<PolyM<N>> point (PolyM<N> &w, int m, int i) {
    Poly<N> p2;
    PV3<Poly<N>> p1 = point(p2);
    w = PolyM<N>(p2, m, i);
    return PV3<PolyM<N>>(PolyM<N>(p1.x, m, i), PolyM<N>(p1.y, m, i),
			 PolyM<N>(p1.z, m, i));
  }
};

class PointOrderPolyM : public Object<PolyM> {
  CPointRecord a, b;

  DeclareCalculate (PolyM) {
    PolyM<N> a2, b2;
    PV3<PolyM<N>> a1 = a.point(a2, 2, 0), b1 = b.point(b2, 2, 1),
      w = b1*a2 - a1*b2;
    PV3<N> r = Rdir->get<N>();
    return (w.x*r.x + w.y*r.y + w.z*r.z)*a2*b2;
  }
public:
  PointOrderPolyM (const CPointRecord &a, const CPointRecord &b)
    : a(a), b(b) {}
};

class LinePolyM : public Object<PolyM> {
  CPointRecord a, t, h;

  DeclareCalculate (PolyM) {
    PolyM<N> a2, t2, h2;
    PV3<PolyM<N>> a1 = a.point(a2, 3, 0), t1 = t.point(t2, 3, 1),
      h1 = h.point(h2, 3, 2),
      w = a1.cross(h1)*t2 + t1.cross(a1)*h2 + h1.cross(t1)*a2;
    PV3<N> r = Rdir->get<N>();
    return (w.x*r.x + w.y*r.y + w.z*r.z)*a2*t2*h2;
  }
public:
  LinePolyM (const CPointRecord &a, const CPointRecord &t, const CPointRecord &h)
    : a(a), t(t), h(h) {}
};

class SidePolyM : public Object<PolyM> {
  CPointRecord p, a, b, c;

  DeclareCalculate (PolyM) {
    PolyM<N> p2, a2, b2, c2;
    PV3<PolyM<N>> p1 = p.point(p2, 4, 0), a1 = a.point(a2, 4, 1),
      b1 = b.point(b2, 4, 2), c1 = c.point(c2, 4, 3);
    return (p1.tripleProduct(b1, c1)*a2 - p1.tripleProduct(b1, a1)*c2 -
	    p1.tripleProduct(a1, c1)*b2 - a1.tripleProduct(b1, c1)*p2)*p2*a2*b2*c2;

  }
public:
  SidePolyM (const CPointRecord &p, const CPointRecord &a, const CPointRecord &b, 
	     const CPointRecord &c) : p(p), a(a), b(b), c(c) {}
};

bool getCPoint (Point *a, CPointRecord &r)
{
  PointCPoint *b = dynamic_cast<PointCPoint *>(a);
  if (!b)
    return false;
  CPoint *c = b->getCP();
  r.a = dynamic_cast<ParametricAngle *>((Angle *) c->getA());
  if (!r.a)
    return false;
  SpiralPoint *s = dynamic_cast<SpiralPoint *>(c);
  if (s) {
    r.s = s->getS();
    return true;
  }
  FFPoint *ff = dynamic_cast<FFPoint *>(c);
  if (ff) {
    r.f1 = ff->getF1();
    r.f2 = ff->getF2();
    return true;
  }
  return false;
}

class AngleOrderO {
public:
  bool operator() (Angle *a, Angle *b) const {
    return a != b && a->order(b);
  }
};

bool rational (ParametricAngle *a)
{
  Root *r = (Root *) (Object<Scalar> *) a->getT();
  return degree(r->getPoly()) == 1;
}

int pointOrderParametric (Point *a, Point *b)
{
 CPointRecord ra, rb;
 if (!(getCPoint(a, ra) && getCPoint(b, rb)))
   return -2;
 PTR<Object<PolyM>> f = new PointOrderPolyM(ra, rb);
 return PolyM2Sign(f, ra.a->getT(), rb.a->getT());
}

int onLineParametric (Point *a, Point *t, Point *h)
{
  CPointRecord ra, rt, rh;
  if (!(getCPoint(a, ra) && getCPoint(t, rt) && getCPoint(h, rh)))
    return -2;
  /* optimization removed for residue testing
  set<ParametricAngle *, AngleOrderO> as = {ra.a, rt.a, rh.a};
  if (as.size() < 3)
    return -2;
  */
  PTR<Object<PolyM>> f = new LinePolyM(ra, rt, rh);
  return PolyM3Sign(f, ra.a->getT(), rt.a->getT(), rh.a->getT());
}

int sideParametric (Point *p, Point *a, Point *b, Point *c)
{
  CPointRecord rp, ra, rb, rc;
  if (!(getCPoint(p, rp) && getCPoint(a, ra) && getCPoint(b, rb) &&
	getCPoint(c, rc)))
    return -2;
  /* optimization removed for residue testing
  set<ParametricAngle *, AngleOrderO> as = {rp.a, ra.a, rb.a, rc.a};
  if (as.size() < 4)
    return -2;
  */
  PTR<Object<PolyM>> f = new SidePolyM(rp, ra, rb, rc);
  int s = PolyM4Sign(f, rp.a->getT(), ra.a->getT(), rb.a->getT(), rc.a->getT());
  return s;
}

unsigned long int nnn = 0;
double tAcp = 0.0, tRes = 0.0;

int main (int argc, char *argv[])
{
  if (argc < 3)
    return 0;
  ifstream astr(argv[1]), bstr(argv[2]);
  if (!(astr.good() && bstr.good()))
    return 0;
  acp::enable();
  double t0 = getTime(), dd = 0.01;
  PTR<Polyhedron> a0 = readPolyhedronVTK(astr), b = readPolyhedronVTK(bstr);
  PTR<Point> t = new Point(-0.5*(a0->bbox[0] + a0->bbox[1]),
			   -0.5*(a0->bbox[2] + a0->bbox[3]),
			   b->bbox[5] - a0->bbox[5] - 0.1);
  PTR<Polyhedron> a = a0->translate(t);
  PTR<Cspace> c = cspace(a, b);
  PTR<Polyhedron> p0 = discretize(c, dd), p = p0->selfUnion();
  t0 = getTime() - t0;
  double ta = tAcp/double(nnn), tr = tRes/double(nnn),
    ra = ta == 0.0 ? 0.0 : tr/ta;
  cerr << "faces: " << p->faces.size() << "; cpu time: " << t0 
       << " predicates: " << nnn << endl << "time per predicate: acp " << ta
       << " residue " << tr << " ratio " << ra << endl;
  acp::disable();
}
