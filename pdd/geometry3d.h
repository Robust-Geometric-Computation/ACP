#ifndef GEOMETRY3D
#define GEOMETRY3D

#include "polyhedron.h"

class PlaneSide : public Primitive {
  Plane *plane;
  Point *point;

  DeclareSign {
    PV3<N> v = plane->get<N>().n, p = point->get<N>();
    N d = plane->get<N>().k;
    return v.dot(p) - d;
  }
 public:
  PlaneSide (Plane *plane, Point *point) : plane(plane), point(point) {}
};

class DiffLength : public Primitive {
  Point *i, *j;

  DeclareSign {
    return i->get<N>().dot(i->get<N>()) - j->get<N>().dot(j->get<N>());
  }
 public:
  DiffLength (Point *i, Point *j) : i(i), j(j) {}
};

class IsTangent : public Primitive {
  Point *v0, *v1, *v2;

  DeclareSign {
    PV3<N> a = v0->get<N>(), b = v1->get<N>(), c = v2->get<N>(),
      n = (b - a).cross(c - a), nsplit(n.y, - n.x, N(0));
    N k = n.dot(a);
    int psa = nsplit.dot(a).sign(), psb = nsplit.dot(b).sign(),
      psc = nsplit.dot(c).sign();
    return psa == psb && psb == psc ? N(0) : N(1);
  }
 public:
  IsTangent (Point *v0, Point *v1, Point *v2) : v0(v0), v1(v1), v2(v2) {}
};

class Orient3D : public Primitive {
  Point *a, *b, *c, *d;

  DeclareSign {
    PV3<N> dd = d->get<N>(), da = a->get<N>() - dd, db = b->get<N>() - dd,
      dc = c->get<N>() - dd;
    return da.cross(db).dot(dc);
  }
 public:
  Orient3D (Point *a, Point *b, Point *c, Point *d) : a(a), b(b), c(c), d(d) {}
};

class InSphere : public Primitive {
  Point *a, *b, *c, *d, *e;

  DeclareSign {
    PV3<N> ee = e->get<N>(), ea = a->get<N>() - ee, eb = b->get<N>() - ee,
      ec = c->get<N>() - ee, ed = d->get<N>() - ee;
    return ea.cross(eb).dot(ec)*ed.dot(ed) -
      eb.cross(ec).dot(ed)*ea.dot(ea) +
      ec.cross(ed).dot(ea)*eb.dot(eb) -
      ed.cross(ea).dot(eb)*ec.dot(ec);
  }
 public:
  InSphere (Point *a, Point *b, Point *c, Point *d, Point *e)
    : a(a), b(b), c(c), d(d), e(e) {}
};

class RotationPoint : public Point {
  PTR<Point> point, sca;

  DeclareCalculate (PV3) {
    N sint = sca->get<N>().x, cost = sca->get<N>().y;
    PV3<N> p = point->get<N>();
    return PV3<N>(cost*p.x - sint*p.y, sint*p.x + cost*p.y, p.z);
  }
  
 public:
  RotationPoint (Point *point, Point *sca) : point(point), sca(sca) {} 
};

class TangentIntersectionPoint : public Point {
  PTR<Point> point, sca;

  DeclareCalculate (PV3) {
    N alpha = sca->get<N>().z;
    PV3<N> p = point->get<N>();
    return PV3<N>(p.x - alpha*p.y, p.y + alpha*p.x, p.z);
  }

 public:
  TangentIntersectionPoint (Point *point, Point *sca) : point(point), sca(sca) {}
};

class ZIntercectPoint : public Point {
  PTR<Point> a, b, c;

  DeclareCalculate (PV3) {
    N t = (c->get<N>().z - a->get<N>().z)/(b->get<N>().z - a->get<N>().z);
    return a->get<N>() + t*(b->get<N>() - a->get<N>());
  }

 public:
 ZIntercectPoint (Point *a, Point *b, Point *c) : a(a), b(b), c(c) {}
};

class SplitPlane : public Plane {
  PTR<Point> v0, v1, v2;
  
  DeclareCalculate (PlaneData) {
    PV3<N> a = v0->get<N>(), b = v1->get<N>(), c = v2->get<N>(),
      n_tri = (b - a).cross(c - a); 
    N k_tri = n_tri.dot(a); 
    PV3<N> n(n_tri.y, - n_tri.x, N(0));
    return PlaneData<N>(n, N(0));
  }

 public:
  SplitPlane (Point *v0, Point *v1, Point *v2) : v0(v0), v1(v1), v2(v2) {}
};

class PointNormalPlane : public Plane {
  PTR<Object<PV3>> point, normal;
  
  DeclareCalculate (PlaneData) {
    PV3<N> p = point->get<N>(), n = normal->get<N>();
    N d = n.dot(p);
  return PlaneData<N>(n, d);
  }

 public:
  PointNormalPlane (Point *point, Point *normal) : point(point), normal(normal) {}
};

class IntersectionPoint : public Point {
  PTR<Point> tail, head;
  PTR<Plane> plane;

  DeclareCalculate (PV3) {
    PV3<N> t = tail->get<N>(), h = head->get<N>(), v = plane->get<N>().n;
    N d = plane->get<N>().k, s = (d - t.dot(v))/(h - t).dot(v);
    return t + (h - t)*s;
  }
  
 public:
  IntersectionPoint (Point *tail, Point *head, Plane *plane)
    : tail(tail), head(head), plane(plane) {}
};

class NdotV : public Primitive {
  PTR<Point> head, tail, a, b, c;
  DeclareSign {
    PV3<N> t = tail->get<N>(), v = head->get<N>() - t,
      n = (c->get<N>() - b->get<N>()).cross(a->get<N>() - b->get<N>());
    return n.dot(v);
  }
 public:
  NdotV (Point *tail, Point *head, Point *a, Point *b, Point *c)
    : tail(tail), head(head), a(a), b(b), c(c) {}
};

class FaceIntersectionPoint : public Point {
  PTR<Point> head, tail, a, b, c;
  HFace * hface;
  
  DeclareCalculate (PV3) {
    PV3<N> t = tail->get<N>(), v = head->get<N>() - t,
      n = (c->get<N>() - b->get<N>()).cross(a->get<N>() - b->get<N>());
    N kt = - n.dot(a->get<N>()), k = - (n.dot(t) + kt)/n.dot(v);
    return t + k*v;
  }

 public:
  FaceIntersectionPoint (Point *tail, Point *head, HFace *hface)
    : tail(tail), head(head), hface(hface)
  { a = hface->getF()->getP()->getA();
    b = hface->getF()->getP()->getB();
    c = hface->getF()->getP()->getC();
    if (NdotV(tail, head, a, b, c) == 0) {
      cerr << "FaceIntersectionPoint failed" << endl;
      exit(0);
    }
  }
  
  HFace * getHFace() { return hface; }
};

#define vol(a,b,c,d) (a-d).cross(b-d).dot(c-d)

inline bool operator< (const Parameter &a, const Parameter &b) {
  return (a - b).sign() < 0;
}

inline bool operator> (const Parameter &a, const Parameter &b) {
  return (a - b).sign() > 0;
}

inline bool operator< (const PParameter &a, const PParameter &b) {
  return (a - b).sign() < 0;
}

inline bool operator> (const PParameter &a, const PParameter &b) {
  return (a - b).sign() > 0;
}

inline bool operator< (const MParameter &a, const MParameter &b) {
  return (a - b).sign() < 0;
}

inline bool operator> (const MParameter &a, const MParameter &b) {
  return (a - b).sign() > 0;
}

class FaceNearestPoint : public Point {
  PTR<Point> point, pa, pb, pc;

  DeclareCalculate(PV3) {
    PV3<N> p = point->get<N>(), a = pa->get<N>(), b = pb->get<N>(),
      c = pc->get<N>(), n = (a - b).cross(c - b);
    N d = n.dot(b), s = (d - p.dot(n))/n.dot(n);
    PV3<N> q = p + n*s;
    bool sideab = vol(q, a, b, b + n).sign() == vol(c, a, b, b + n).sign();
    bool sideac = vol(q, a, c, c + n).sign() == vol(b, a, c, c + n).sign();
    bool sidebc = vol(q, b, c, c + n).sign() == vol(a, b, c, c + n).sign();
    if (sideac && sideab && sidebc) 
      return q;
    if (sideab && sideac) {
      N t = (q-c).dot(b-c);
      if (t.sign() < 0)
	return c;
      if ((t - N(1)).sign() > 0)
	return b;
      return c + t*(b-c);
    }
    if (sideab && sidebc) {
      N t = (q-c).dot(a-c);
      if (t.sign() < 0)
	return c;
      if ((t - N(1)).sign() > 0)
	return a;
      return c + t*(a-c);
    }
    if (sideac && sidebc) {
      N t = (q-a).dot(b-a);
      if (t.sign() < 0)
	return a;
      if ((t - N(1)).sign() > 0)
	return b;
      return a + t*(b-a);    
    }
    if (!sideab && !sideac)
      return a;
    if (!sideab && !sidebc)
      return b;
    if (!sideac && !sidebc)
      return c;
    cerr << "FaceNearestPoint failed" << endl;
    exit(0);
    return p;
  }
 
public: 
  FaceNearestPoint (Point *point, HFace *hf) : point(point) {
    pa = hf->getF()->getP()->getA();
    pb = hf->getF()->getP()->getB();
    pc = hf->getF()->getP()->getC();
  }
};



#endif
