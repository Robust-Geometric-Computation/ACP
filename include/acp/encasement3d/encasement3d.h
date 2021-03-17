#include "acp/encasement3d/disjointset.h"
#include "acp/poly/poly.h"
#include "acp/poly/poly3d.h"
#include "acp/core/timer.h"

// mesh
template <class N>
Poly<N> remainder(const Poly<N>& a, const Poly<N>& b, Poly<N>& q) {
  N lb = b.lc();
  int db = b.deg();
  Poly<N> r = a;
  while (db <= r.deg()) {
    N k = r.lc() / lb;
    int dk = r.deg() - db;
    q.add(k, dk);
    r.a.pop_back();
    for (int i = 0; i < db; ++i) r.add(-k * b.a[i], i + dk);
  }
  return r;
}

// mesh
// Sign is zero if poly has a double root.
class NoDoubleZero : public Primitive {
  PTR<Object<Poly>> poly;
  DeclareSign {
    Poly<N> a = poly->get<N>();
    Poly<N> b = a.der();
    while (b.deg() > 0) {
      Poly<N> q;
      Poly<N> r = remainder(a, b, q);
      a = b;
      b = r;
    }
    return b.a[0];
  }

 public:
  NoDoubleZero(PTR<Object<Poly>> poly) : poly(poly) {}
};

// mesh
// poly has a double root.  Removes that root (both of them).
class SimplePoly : public Object<Poly> {
  PTR<Object<Poly>> poly;
  DeclareCalculate(Poly) {
    Poly<N> f = poly->get<N>();
    Poly<N> a = f;
    Poly<N> b = a.der();
    while (b.deg() > 1) {
      Poly<N> q;
      Poly<N> r = remainder(a, b, q);
      a = b;
      b = r;
    }
    Poly<N> q;
    Poly<N> r = remainder(f, b * b, q);
    return q;
  }

 public:
  SimplePoly(PTR<Object<Poly>> poly) : poly(poly) {}
};

class V3I : public Object<Scalar> {
  PTR<Object<PV3>> p;
  int i;
  DeclareCalculate(Scalar) { return p->get<N>()[i]; }

 public:
  V3I(PTR<Object<PV3>> p, int i) : p(p), i(i) {}
};

class P3ISIS : public Object<Poly> {
  PTR<Object<Poly3D>> p3d;
  int ix, iy;
  PTR<Object<Scalar>> sx, sy;
  DeclareCalculate(Poly) {
    return p3d->get<N>().substitute2(ix, sx->get<N>().x, iy, sy->get<N>().x);
  }

 public:
  P3ISIS(PTR<Object<Poly3D>> p3d, int ix, PTR<Object<Scalar>> sx, int iy,
         PTR<Object<Scalar>> sy)
      : p3d(p3d), ix(ix), sx(sx), iy(iy), sy(sy) {}
};

class P3V3xy : public Object<Poly> {
  PTR<Object<Poly3D>> p3d;
  PTR<Object<PV3>> p;
  DeclareCalculate(Poly) {
    PV3<N> p = this->p->get<N>();
    return p3d->get<N>().substitute2(0, p.x, 1, p.y);
  }

 public:
  P3V3xy(PTR<Object<Poly3D>> p3d, PTR<Object<PV3>> p) : p3d(p3d), p(p) {}
};

class SS : public Object<PV2> {
  PTR<Object<Scalar>> sx, sy;
  DeclareCalculate(PV2) {
    PV2<N> p;
    p.x = sx->get<N>().x;
    p.y = sy->get<N>().x;
    return p;
  }

 public:
  SS(PTR<Object<Scalar>> sx, PTR<Object<Scalar>> sy) : sx(sx), sy(sy) {}
};

class SSS : public Object<PV3> {
  PTR<Object<Scalar>> sx, sy, sz;
  DeclareCalculate(PV3) {
    PV3<N> p;
    p.x = sx->get<N>().x;
    p.y = sy->get<N>().x;
    p.z = sz->get<N>().x;
    return p;
  }

 public:
  SSS(PTR<Object<Scalar>> sx, PTR<Object<Scalar>> sy, PTR<Object<Scalar>> sz)
      : sx(sx), sy(sy), sz(sz) {}
};

class V3xyS : public Object<PV3> {
  PTR<Object<PV3>> p3;
  PTR<Object<Scalar>> z;
  DeclareCalculate(PV3) {
    PV3<N> p3 = this->p3->get<N>();
    PV3<N> p;
    p.x = p3.x;
    p.y = p3.y;
    p.z = z->get<N>().x;
    return p;
  }

 public:
  V3xyS(PTR<Object<PV3>> p3, PTR<Object<Scalar>> z) : p3(p3), z(z) {}
};

class V3ItoV2 : public Object<PV2> {
  PTR<Object<PV3>> p3;
  int xyz;
  DeclareCalculate(PV2) {
    PV3<N> p3 = this->p3->get<N>();
    PV2<N> p;
    p.x = p3[(xyz + 1) % 3];
    p.y = p3[(xyz + 2) % 3];
    return p;
  }

 public:
  V3ItoV2(PTR<Object<PV3>> p3, int xyz) : p3(p3), xyz(xyz) {}
};

class P3derI : public Object<Poly3D> {
  PTR<Object<Poly3D>> p;
  int xyz;
  DeclareCalculate(Poly3D) { return p->get<N>().der(xyz); }

 public:
  P3derI(PTR<Object<Poly3D>> p, int xyz) : p(p), xyz(xyz) {}
};

class P3derIV3 : public Primitive {
  PTR<Object<Poly3D>> p;
  int xyz;
  PTR<Object<PV3>> v;
  DeclareSign { return p->get<N>().der(xyz, v->get<N>()); }

 public:
  P3derIV3(PTR<Object<Poly3D>> p, int xyz, PTR<Object<PV3>> v)
      : p(p), xyz(xyz), v(v) {}
};

class P2derI : public Object<Poly2D> {
  PTR<Object<Poly2D>> p;
  int xy;
  DeclareCalculate(Poly2D) { return p->get<N>().der(xy); }

 public:
  P2derI(PTR<Object<Poly2D>> p, int xy) : p(p), xy(xy) {}
};

class P2derIV2 : public Primitive {
  PTR<Object<Poly2D>> p;
  int xy;
  PTR<Object<PV2>> v;
  DeclareSign { return p->get<N>().der(xy, v->get<N>()); }

 public:
  P2derIV2(PTR<Object<Poly2D>> p, int xy, PTR<Object<PV2>> v)
      : p(p), xy(xy), v(v) {}
};

class DFxDGI : public Object<Poly3D> {
  PTR<Object<Poly3D>> f, g;
  int xyz;
  DeclareCalculate(Poly3D) {
    Poly3D<N> f = this->f->get<N>();
    Poly3D<N> g = this->g->get<N>();
    int ix = (xyz + 1) % 3, iy = (xyz + 2) % 3;
    return f.der(ix) * g.der(iy) - f.der(iy) * g.der(ix);
  }

 public:
  DFxDGI(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g, int xyz)
      : f(f), g(g), xyz(xyz) {}
};

class DFxDGIV3 : public Primitive {
  PTR<Object<Poly3D>> f, g;
  int xyz;
  PTR<Object<PV3>> v;
  DeclareSign {
    Poly3D<N> f = this->f->get<N>();
    Poly3D<N> g = this->g->get<N>();
    int ix = (xyz + 1) % 3, iy = (xyz + 2) % 3;
    PV3<N> v = this->v->get<N>();
    return f.der(ix, v) * g.der(iy, v) - f.der(iy, v) * g.der(ix, v);
  }

 public:
  DFxDGIV3(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g, int xyz,
           PTR<Object<PV3>> v)
      : f(f), g(g), xyz(xyz), v(v) {}
};

class SminusS : public Primitive {
  PTR<Object<Scalar>> s, t;
  DeclareSign { return s->get<N>().x - t->get<N>().x; }

 public:
  SminusS(PTR<Object<Scalar>> s, PTR<Object<Scalar>> t) : s(s), t(t) {}
};

class V3IminusD : public Primitive {
  PTR<Object<PV3>> p;
  int xyz;
  double d;
  DeclareSign { return p->get<N>()[xyz] - d; }

 public:
  V3IminusD(PTR<Object<PV3>> p, int xyz, double d) : p(p), xyz(xyz), d(d) {}
};

class V3IminusS : public Primitive {
  PTR<Object<PV3>> p;
  int xyz;
  PTR<Object<Scalar>> s;
  DeclareSign { return p->get<N>()[xyz] - s->get<N>().x; }

 public:
  V3IminusS(PTR<Object<PV3>> p, int xyz, PTR<Object<Scalar>> s)
      : p(p), xyz(xyz), s(s) {}
};

class V3minusV3I : public Primitive {
  PTR<Object<PV3>> p, q;
  int xyz;
  DeclareSign { return p->get<N>()[xyz] - q->get<N>()[xyz]; }

 public:
  V3minusV3I(PTR<Object<PV3>> p, PTR<Object<PV3>> q, int xyz)
      : p(p), q(q), xyz(xyz) {}
};

class V2IminusD : public Primitive {
  PTR<Object<PV2>> p;
  int xy;
  double d;
  DeclareSign { return p->get<N>()[xy] - d; }

 public:
  V2IminusD(PTR<Object<PV2>> p, int xy, double d) : p(p), xy(xy), d(d) {}
};

// mesh
static bool contains(const PV2<Parameter>& box, PTR<Object<PV2>> p) {
  return (
      V2IminusD(p, 0, box[0].lb()) >= 0 && V2IminusD(p, 0, box[0].ub()) <= 0 &&
      V2IminusD(p, 1, box[1].lb()) >= 0 && V2IminusD(p, 1, box[1].ub()) <= 0);
}

class V2minusV2I : public Primitive {
  PTR<Object<PV2>> p, q;
  int xy;
  DeclareSign { return p->get<N>()[xy] - q->get<N>()[xy]; }

 public:
  V2minusV2I(PTR<Object<PV2>> p, PTR<Object<PV2>> q, int xy)
      : p(p), q(q), xy(xy) {}
};

class V3mV3ImV3mV3I : public Primitive {
  PTR<Object<PV3>> pi, qi, pj, qj;
  int i, j;
  DeclareSign {
    return (pi->get<N>()[i] - qi->get<N>()[i]) -
           (pj->get<N>()[j] - qj->get<N>()[j]);
  }

 public:
  V3mV3ImV3mV3I(PTR<Object<PV3>> pi, PTR<Object<PV3>> qi, int i,
                PTR<Object<PV3>> pj, PTR<Object<PV3>> qj, int j)
      : pi(pi), qi(qi), i(i), pj(pj), qj(qj), j(j) {}
};

class XmXmYmY : public Primitive {
  PTR<Object<PV2>> px, qx, py, qy;
  DeclareSign {
    return (px->get<N>().x - qx->get<N>().x) -
           (py->get<N>().y - qy->get<N>().y);
  }

 public:
  XmXmYmY(PTR<Object<PV2>> px, PTR<Object<PV2>> qx, PTR<Object<PV2>> py,
          PTR<Object<PV2>> qy)
      : px(px), qx(qx), py(py), qy(qy) {}
};

// Solve for root of f=0 for f in fs and x_{xyz}=d for xyz,d in xyzs,ds.
vector<PTR<Object<PV3>>> getRoots(vector<PTR<Object<Poly3D>>> fs,
                                  vector<int> xyzs, vector<double> ds,
                                  const PV3<Parameter>& box) {
  assert(fs.size() + xyzs.size() == 3);
  assert(xyzs.size() == ds.size());

  if (xyzs.size() == 0) {
    return getRoots(fs[0], fs[1], fs[2], box);
  }

  vector<PTR<Object<PV3>>> roots;

  if (xyzs.size() == 1) {
    int iz = xyzs[0];
    int ix = (iz + 1) % 3, iy = (iz + 2) % 3;
    PTR<Object<Scalar>> s = new Object<Scalar>(Parameter::constant(ds[0]));
    PTR<Object<Poly2D>> f2 = new Sub1dPoly3D(fs[0], iz, s);
    PTR<Object<Poly2D>> g2 = new Sub1dPoly3D(fs[1], iz, s);
    PV2<Parameter> box2(box[ix], box[iy]);
    vector<PTR<Object<PV2>>> roots2d = getRoots(f2, g2, box2);
    for (auto it = roots2d.begin(); it != roots2d.end(); ++it)
      roots.push_back(new P2to3(*it, s, iz));
  } else if (xyzs.size() == 2) {
    int ix = xyzs[0];
    int iy = xyzs[1];
    assert(ix != iy);
    int iz = 3 - ix - iy;
    PTR<Object<Scalar>> sx = new Object<Scalar>(Parameter::constant(ds[0]));
    PTR<Object<Scalar>> sy = new Object<Scalar>(Parameter::constant(ds[1]));
    assert(0 <= ix && ix < 3);
    assert(0 <= iy && iy < 3);
    PTR<Object<Poly>> f1 = new P3ISIS(fs[0], ix, sx, iy, sy);
    PTR<Object<Scalar>> lbz = new Object<Scalar>(box[iz].lbP());
    PTR<Object<Scalar>> ubz = new Object<Scalar>(box[iz].ubP());
    vector<PTR<Object<Scalar>>> roots1 = acp::getRoots(f1, lbz, ubz);
    PTR<Object<Scalar>> s[3];
    s[ix] = sx;
    s[iy] = sy;
    for (int h = 0; h < roots1.size(); ++h) {
      s[iz] = roots1[h];
      roots.push_back(new SSS(s[0], s[1], s[2]));
    }
  } else {
    assert(xyzs[0] != xyzs[1] && xyzs[1] != xyzs[2] && xyzs[2] != xyzs[0]);
    PV3<Parameter> r;
    for (int i = 0; i < 3; ++i) r[xyzs[i]] = Parameter::constant(ds[i]);
    roots.push_back(new Object<PV3>(r));
  }

  return roots;
}

// polynomial x - c
Poly3D<Parameter> slicePoly(double c) {
  Poly3D<Parameter> p(1, 0, 0);
  p[0] = Parameter::constant(-c);
  p[1] = 1;
  return p;
}

// polynomial -x_{xyz} + box[xyz].lb() or x - box[xyz].ub()
// side of box with outward normal
Poly3D<Parameter> boxPoly(const PV3<Parameter>& box, int xyz, int lu) {
  int deg[3] = {0, 0, 0};
  deg[xyz] = 1;
  Poly3D<Parameter> p(deg[0], deg[1], deg[2]);
  int s = -1 + 2 * lu;
  p[0] = Parameter::constant((lu == 0 ? box[xyz].lb() : box[xyz].ub()) * -s);
  p[1] = s;
  return p;
}

class Encasement3D;

// mesh
class Vert3DInfo {
 public:
  int type, jf, jg, jh;
  Vert3DInfo(int type, int jf, int jg, int jh)
      : type(type), jf(jf), jg(jg), jh(jh) {}
  Vert3DInfo() { type = jf = jg = jh = -1; }
};

// Encasement of polynomial f in an x-slices
class FXEncasement {
 public:
  // fIndex:  index of f in e3
  // xValue:  x value of slice
  // box:  yz box of 2D encasement
  FXEncasement(Encasement3D* e3, int fIndex, double xValue, PV2<Parameter> box);

  // mesh
  FXEncasement(Encasement3D* e3, int fIndex, PTR<Object<PV3>> v3,
               Vert3DInfo v3info);
  PTR<Object<Poly3D>> f3;
  PTR<Object<PV3>> v3;
  PV3<Parameter> v3p;
  Vert3DInfo v3info;
  PTR<Object<PV2>> identicalRoot;
  class Vert;
  PTR<Vert> identicalVert;

  Encasement3D* e3;
  int fIndex;
  double xValue;
  PTR<Object<Scalar>> xObject;

  // boundary and input polynomials parallel to e3->bf
  // bf[0] and bf[1] (x-boundary polys in e3) are null
  // for all others, substitute x=xValue
  vector<PTR<Object<Poly2D>>> bf;
  PV2<Parameter> box;

  // fyz[0] and fyz[1] are really fx and fy because we are in yz plane.
  PTR<Object<Poly2D>> f, fyz[2];

  class Edge;

  enum VertType {
    FG = 0,  // intersection of f with an input or boundary surface (index jg)
    FD = 1,  // intersection of f with fyz[0] or fyz[1] (jg = 0 or 1)
    FB = 2   // intersection of f with cel boundary y=c or z=c (jg = 0 or 1)
  };

  class Cel {
   public:
    Cel() : yz(-1) {
      spv.s = 0;
      lb.c = rb.c = 0;
    }
    ~Cel() {
      if (0 <= yz && yz <= 2) {
        delete lb.c;
        delete rb.c;
      }
    }

    // debug
    PV2<Parameter> cbox;

    // 0,1 if internal and split on y,z.
    // 2 if split on plane
    // -1 if leaf
    int yz;

    // split value or pointer to splitting line or internal vertex
    union {
      double s;
      Object<Poly2D>* p;
      Vert* v;
    } spv;

    // left cel or first boundary vertex
    union {
      Cel* c;
      Vert* b;
    } lb;

    // right cel or second boundary vertex
    union {
      Cel* c;
      Vert* b;
    } rb;

    void setbverts(vector<Vert*>& bverts) {
      assert(bverts.size() <= 2);
      assert(lb.b == 0);
      if (bverts.size() > 0) lb.b = bverts[0];
      assert(rb.b == 0);
      if (bverts.size() > 1) rb.b = bverts[1];
    }

    // add turning vertex f=fyz[jg]=0, recursive
    void addT(Vert* v);

    // add all boundary vertices, recursive
    void addB(FXEncasement* e, vector<Vert*>& bverts,
              const PV2<Parameter>& box);

    // add f=g=0 vertices, recursive
    // p and side are set if ancestor split on line
    void addI(FXEncasement* e, Vert* v, const PV2<Parameter>& box,
              Object<Poly2D>* p = 0, int side = 0);

    // Set the Cel fields of all vertices inside or bounding leaf Cels.
    void setVertCels(FXEncasement* e);

    // f(p)=0.  find the Edge that contains p.
    Edge* getEdge(PTR<Object<PV2>> p);
  } top;

  class Vert : public RefCnt {
   public:
    // see VertType
    Vert(VertType type, PTR<Object<PV2>> p, int jg)
        : type(type), p(p), jg(jg), boxvert(false), twin(0) {
      c[0] = c[1] = 0;
      e[0] = e[1] = 0;
    }
    VertType type;
    PTR<Object<PV2>> p;
    int jg;
    Cel* c[2];
    bool boxvert;
    // for boundary Vert, e[0] is arrangement Edge
    // else e[i]->v[i] == v[i], i=0,1
    Edge* e[2];  // arrangement Edge of boundary Vert
    // mesh
    Vert* twin;  // identical vertex

    // inform this boundary Vert of one of its two leaf cells
    void setBCel(Cel* bc) {
      if (c[0] == 0)
        c[0] = bc;
      else if (c[1] == 0)
        c[1] = bc;
      else
        assert(0);
    }

    // for disjoint-set tree
    int& size() { return jg; }
    Vert*& parent() { return *(Vert**)&e[0]; }
    Edge*& edge() { return e[0]; }
  };

  // fg[ig] is set of solutions to f(p)=g(p)=0 for g=bf[ig]
  vector<vector<PTR<Vert>>> fg;

  // solutions to f(p)=fyz[i](p)=0, for i=0,1
  vector<PTR<Vert>> fd[2];

  // solve f(y,z)=0 and y=d or z=d (yz=0 or 1) inside box
  vector<PTR<Object<PV2>>> getRoots(int yz, double d,
                                    const PV2<Parameter>& box) {
    PTR<Object<Scalar>> s = new Object<Scalar>(Parameter::constant(d));
    PTR<Object<Poly>> f1 = new Sub1Poly2D(f, yz, s);

    // mesh
    if (NoDoubleZero(f1) == 0) f1 = new SimplePoly(f1);

    PTR<Object<Scalar>> lb =
        new Object<Scalar>(Parameter::constant(box[!yz].lb()));
    PTR<Object<Scalar>> ub =
        new Object<Scalar>(Parameter::constant(box[!yz].ub()));
    vector<PTR<Object<Scalar>>> roots1 = acp::getRoots(f1, lb, ub);
    vector<PTR<Object<PV2>>> roots(roots1.size());
    for (int i = 0; i < roots.size(); i++)
      roots[i] = new P1to2(roots1[i], s, yz);
    return roots;
  }

  // solve f=g=0 for g=bf[ig]
  vector<PTR<Object<PV2>>> getRoots(int ig) {
    if (ig < 2 || ig == fIndex) {
      vector<PTR<Object<PV2>>> roots;
      return roots;
    }

    PV2<Parameter> bbox = box;

    // mesh
    bool hasIdentity =
        (v3info.type == 0 /*FGH*/ &&
         (ig == v3info.jf || ig == v3info.jg || ig == v3info.jh)) ||
        (v3info.type == 1 /*FGD*/ && (v3info.jf == ig || v3info.jg == ig) &&
         v3info.jh != 0) ||
        (v3info.type == 2 /*FGB*/ && (ig == v3info.jf || ig == v3info.jg)) ||
        (v3info.type == 3 /*FFzH*/ &&
         ((v3info.jf == fIndex && v3info.jh == ig) ||
          (v3info.jh == fIndex && v3info.jf == ig))) ||
        (v3info.type == 6 /*FSH*/ && v3info.jf == ig);
    bool expand = hasIdentity && v3info.type == 0 &&
                  ((ig >= 6 && v3info.jf < 6) || v3info.jg < 6);
    int je = ig == v3info.jf ? v3info.jg : v3info.jf;
    if (expand) {
      int yz = (je - 2) / 2;
      if (je % 2 == 0)
        bbox[yz] = Parameter::constant(box[yz].lb() - 1e-6)
                       .interval(Parameter::constant(box[yz].ub()));
      else
        bbox[yz] = Parameter::constant(box[yz].lb())
                       .interval(Parameter::constant(box[yz].ub() + 1e-6));
    }
    bool hasTangency =
        (v3info.type == 1 /*FGD*/ && (v3info.jf == ig || v3info.jg == ig) &&
         v3info.jh == 0 && ig >= 6);
    if (hasTangency) {
      extern PTR<Object<PV2>> degenerateIntersection;
      degenerateIntersection = new V3ItoV2(v3, 0);
    }

    vector<PTR<Object<PV2>>> roots;
    if (ig < 6) {
      int yz = ig / 2 - 1;
      int lu = ig % 2;
      double d = lu == 0 ? box[yz].lb() : box[yz].ub();
      roots = getRoots(yz, d, bbox);
    } else {
      if (v3info.type == 2 /*FGB*/ && v3p.x.lb() < v3p.x.ub()) {
	extern int yzAvoidSplit;
	extern double sAvoidSplit;
	if (v3p.y.lb() == v3p.y.ub()) {
	  yzAvoidSplit = 1;
	  sAvoidSplit = v3p.y.lb();
	}
	else if (v3p.z.lb() == v3p.z.ub()) {
	  yzAvoidSplit = 2;
	  sAvoidSplit = v3p.z.lb();
	}
	else
	  assert(0);
      }

      roots = ::getRoots(bf[ig], f, bbox);

      {
	extern int yzAvoidSplit;
	yzAvoidSplit = 0;
      }
    }

    // mesh
    if (hasTangency) {
      extern PTR<Object<PV2>> degenerateIntersection;
      roots.push_back(degenerateIntersection);
      identicalRoot = degenerateIntersection;
      degenerateIntersection = 0;
      return roots;
    }
    if (hasIdentity) {
      PV2<Parameter> p55(Parameter::constant(-0.5), Parameter::constant(0.5));
      Poly2D<Parameter> f2p = f->getApprox(1);
      Parameter f2pofp55 = f2p.value(p55);

      bool failed = false;
      for (PTR<Object<PV2>> r : roots) {
        PV2<Parameter> rp = r->getApprox(1);
        if (rp.x.intersects(v3p.y) && rp.y.intersects(v3p.z)) {
          if (identicalRoot != 0) {
            failed = true;
            break;
          }
          identicalRoot = r;
        }
      }
      if (failed) {
        identicalRoot = 0;
        for (PTR<Object<PV2>> r : roots) {
          PV2<Parameter> rp = r->getApprox();
          if (rp.x.intersects(v3p.y) && rp.y.intersects(v3p.z)) {
            assert(identicalRoot == 0);
            identicalRoot = r;
          }
        }
      }
      assert(identicalRoot != 0);
    }
    if (expand)
      for (int i = roots.size(); --i > 0;)
        if (roots[i] != identicalRoot && Poly2Dval(bf[je], roots[i]) > 0)
          roots.erase(roots.begin() + i);

    return roots;
  }

  vector<PTR<Object<PV2>>> getRootsYZ(int yz) {
    PV2<Parameter> bbox = box;

    // mesh
    bool hasIdentity =
        (v3info.type == 1 /*FGD*/ && v3info.jf < 6 && v3info.jh == 0 &&
         (v3info.jf - 2) / 2 == !yz) ||
        (v3info.type == 3 /*FFzH*/ && v3info.jf == fIndex && yz == 1) ||
        (v3info.type == 4 /*FFzD*/ &&
         (yz == 1 ||
          (yz == 0 && v3info.jh == 0 && P3derIV3(f3, 1, v3) == 0))) ||
        (v3info.type == 5 /*FFzB*/ && yz == 1);
    bool expand = (hasIdentity && (v3info.type == 1 ||
                                   (v3info.type == 3 && v3info.jh < 6)));  // ||
    // (v3info.type == 0 && v3info.jf < 6 && v3info.jg < 6
    int je = v3info.type == 1 ? v3info.jf : v3info.jh;
    if (expand) {
      int eyz = (je - 2) / 2;
      if (je % 2 == 0)
        bbox[eyz] = Parameter::constant(box[eyz].lb() - 1e-6)
                        .interval(Parameter::constant(box[eyz].ub()));
      else
        bbox[eyz] = Parameter::constant(box[eyz].lb())
                        .interval(Parameter::constant(box[eyz].ub() + 1e-6));
    }
    bool hasSingularity = (v3info.type == 4 /*FFzD*/ && v3info.jh == 0 &&
                           (P3derIV3(f3, 1, v3) == 0 || yz == 1));
    if (hasSingularity) {
      extern PTR<Object<PV2>> degenerateIntersection;
      degenerateIntersection = new V3ItoV2(v3, 0);
    }

    if (v3info.type == 2 /*FGB*/ && v3p.x.lb() < v3p.x.ub()) {
      extern int yzAvoidSplit;
      extern double sAvoidSplit;
      if (v3p.y.lb() == v3p.y.ub()) {
	yzAvoidSplit = 1;
	sAvoidSplit = v3p.y.lb();
      }
      else if (v3p.z.lb() == v3p.z.ub()) {
	yzAvoidSplit = 2;
	sAvoidSplit = v3p.z.lb();
      }
      else
	assert(0);
    }

    vector<PTR<Object<PV2>>> roots = ::getRoots(fyz[yz], f, bbox);

    {
      extern int yzAvoidSplit;
      yzAvoidSplit = 0;
    }

    // mesh
    if (hasSingularity) {
      extern PTR<Object<PV2>> degenerateIntersection;
      roots.push_back(degenerateIntersection);
      identicalRoot = degenerateIntersection;
      degenerateIntersection = 0;
      return roots;
    }
    if (hasIdentity) {
      bool failed = false;
      for (PTR<Object<PV2>> r : roots) {
        PV2<Parameter> rp = r->getApprox(1);
        if (rp.x.intersects(v3p.y) && rp.y.intersects(v3p.z)) {
          if (identicalRoot != 0) {
            failed = true;
            break;
          }
          identicalRoot = r;
        }
      }
      if (failed) {
        identicalRoot = 0;
        for (PTR<Object<PV2>> r : roots) {
          PV2<Parameter> rp = r->getApprox();
          if (rp.x.intersects(v3p.y) && rp.y.intersects(v3p.z)) {
            assert(identicalRoot == 0);
            identicalRoot = r;
          }
        }
      }
      assert(identicalRoot != 0);
    }
    if (expand)
      for (int i = roots.size(); --i > 0;)
        if (roots[i] != identicalRoot && Poly2Dval(bf[je], roots[i]) > 0)
          roots.erase(roots.begin() + i);

    return roots;
  }

  void initVerts() {
    // mesh
    // Swapped the two for-loops.  I hope this doesn't break anything!

    for (int yz = 0; yz < 2; ++yz) {
      vector<PTR<Object<PV2>>> roots = getRootsYZ(yz);
      fd[yz].resize(roots.size());
      for (int h = 0; h < roots.size(); ++h) {
        fd[yz][h] = new Vert(FD, roots[h], yz);

        // mesh
        if (roots[h] == identicalRoot) {
          if (identicalVert == 0)
            identicalVert = fd[yz][h];
          else {
            assert(identicalVert->twin == 0);
            identicalVert->twin = fd[yz][h];
          }
          identicalRoot = 0;
        }
      }
    }

    fg.resize(bf.size());
    for (int ig = 2; ig < bf.size(); ++ig) {
      vector<PTR<Object<PV2>>> roots = getRoots(ig);
      fg[ig].resize(roots.size());
      for (int h = 0; h < roots.size(); ++h) {
        fg[ig][h] = new Vert(FG, roots[h], ig);

        // mesh
        if (roots[h] == identicalRoot) {
          if (identicalVert == 0)
            identicalVert = fg[ig][h];
          else {
            assert(identicalVert->twin == 0);
            identicalVert->twin = fg[ig][h];
          }
          identicalRoot = 0;
        }
      }
    }
  }

  // Cel boundary vertices: f=0 and y=c or z=c
  vector<PTR<Vert>> fb;

  // splitting lines used in tree, saved in one place for memory management
  // purposes
  vector<PTR<Object<Poly2D>>> lines;

  void initCels() {
    for (int yz = 0; yz < 2; ++yz)
      for (Vert* v : fd[yz]) top.addT(v);

    // boxvert.txt
    assert(fIndex >= 6);  // no "flat" cases to worry about?
    vector<Vert*> bverts;
#ifdef BLEEN
    for (int yz = 0; yz < 2; yz++) {
      for (int lu = 0; lu < 2; lu++) {
        double d = lu == 0 ? box[yz].lb() : box[yz].ub();
        vector<PTR<Object<PV2>>> broots = getRoots(yz, d, box);
        for (int h = 0; h < broots.size(); ++h) {
          Vert* bvert = new Vert(FB, broots[h], yz);
          fb.push_back(bvert);
          bverts.push_back(bvert);
          bvert->boxvert = true;
        }
      }
    }
#endif
    for (int ig = 2; ig < 6; ++ig) {
      int yz = (ig - 2) / 2;
      for (Vert* vert : fg[ig]) {
        if (identicalVert != 0 && vert == identicalVert->twin) continue;
        Vert* bvert = new Vert(FB, vert->p, yz);
        fb.push_back(bvert);
        bverts.push_back(bvert);
        bvert->boxvert = true;
      }
    }
    int nboxbverts = bverts.size();

    top.addB(this, bverts, box);

    for (int ig = 0; ig < bf.size(); ++ig)
      for (Vert* v : fg[ig]) top.addI(this, v, box);

    top.setVertCels(this);

    fb.erase(fb.begin(), fb.begin() + nboxbverts);
  }

  class Edge : public RefCnt {
   public:
    Edge(Vert* b) : b(b) { v[0] = v[1] = 0; }
    Vert *b, *v[2];
  };

  vector<PTR<Edge>> edges;

  // find all arrangement edges
  void findEdges();

  static void expand(Edge* edge, Vert* v[2]);
};

class Encasement3D {
 public:
  Encasement3D(const PV3<Parameter>& box, vector<PTR<Object<Poly3D>>> polys);

  PV3<Parameter> box;

  // box polys and f polys.  bf.size() == 6 + polys.size()
  vector<PTR<Object<Poly3D>>> bf;

  int ind(int i, int j, int k) {
    if (k != bf.size())
      assert(!((j < 6 && i / 2 == j / 2) || (k < 6 && j / 2 == k / 2)));
    return i + j * (j - 1) / 2 + k * (k - 1) * (k - 2) / 6;
  }

  class Vert;
  class Edge;
  class Loop;
  class Face;
  class Shell;
  class Cell;

  enum VertType {
    FGH = 0,   // three input or boundary surfaces (jf < jg < jh)
    FGD = 1,   // two input or boundary (jf < jg) and
               // (dF x dG) x, y or z (jh = 0, 1, or 2)
    FGB = 2,   // two input or boundary (jf < jg) and
               // a Cel boundary (jh == 0,1,2 for x,y,z=c)
    FFzH = 3,  // input (jf >= 6), fz (jg == -1), and
               // input or boundary (jh != jf)
    FFzD = 4,  // input (jf >= 6), fz (jg == -1), and
               // (dF x dFz) x, y or z (jh = 0, 1, or 2)
    FFzB = 5,  // input (jf >= 6), fz (jg == -1), and
               // a Cel boundary (jh == 0,1,2 for x,y,z=c)
    FSH = 6,   // input (jf >= 6), slice (jg), and
               // input or boundary (jh != jf, jh >= 0) or Fz (jh == -1)
    FSD = 7,   // input (jf >= 6), slice (jg), and
               // (dF x dS) x, y or z (jh = 0, 1, or 2) (never 0)
    FSB = 8,   // input (jf >= 6), slice (jg), and
               // a Cel boundary (jh == 0,1,2 for x,y,z=c)
    FSFz = 9,  // input (jf >= 6), slice (jg), and
               // Fz (jh == -1)  DEPRECATED
    FXY = 10
  };  // input or boundary (jf) intersecting vertical line (jg=jh=-1)

  // encasement subdivision tree Cel
  class Cel {
   public:
    Cel() : xyz(-1) {
      spv.s = 0;
      lb.c = rb.c = 0;
    }
    ~Cel() {
      if (0 <= xyz && xyz <= 3) {
        delete lb.c;
        delete rb.c;
      }
    }

    // 0,1,2 if internal and split on x,y,z.
    // 3 if split on plane
    // -1 if leaf
    int xyz;

    // split value or pointer to splitting plane or internal vertex
    union {
      double s;
      Object<Poly3D>* p;
      Vert* v;
    } spv;

    // left Cel or first boundary vertex
    union {
      Cel* c;
      Vert* b;
    } lb;

    // right Cel or second boundary vertex
    union {
      Cel* c;
      Vert* b;
    } rb;

    // store contents of bverts in lb.b and rb.b
    void setbverts(vector<Vert*>& bverts) {
      assert(bverts.size() <= 2);
      assert(lb.b == 0);
      if (bverts.size() > 0) lb.b = bverts[0];
      assert(rb.b == 0);
      if (bverts.size() > 1) rb.b = bverts[1];
    }

    // add turning points
    void addT(Vert* v);

    // add boundary vertices, intersections of f and g with boundaries of Cels
    void addB(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g, VertType bType,
              int jf, int jg, vector<Vert*>& bverts,
              vector<PTR<Vert>>& bvertsAll, vector<PTR<Object<Poly3D>>>& planes,
              const PV3<Parameter>& box);

    // add intersection vertices, f=g=h=0
    void addI(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g, VertType bType,
              int jf, int jg, Vert* v, vector<PTR<Vert>>& bvertsAll,
              int xyzAvoid, const PV3<Parameter>& box, Object<Poly3D>* p = 0,
              int side = 0);

    // inform Verts of incidence with leaf Cels
    void setVertCels(Encasement3D* e, Object<Poly3D>* f, Object<Poly3D>* g);

    // f(p)=g(p)=0.  get arrangement Edge that contains p
    Edge* getEdge(PTR<Object<PV3>> p, int xyzAvoid);
  };

  // arrangement vertex
  class Vert : public RefCnt {
   public:
    // see VertType
    Vert(VertType type, PTR<Object<PV3>> p, int jf, int jg, int jh)
        : type(type), p(p), jf(jf), jg(jg), jh(jh) {
      c[0] = c[1] = c[2] = 0;
    }
    VertType type;
    PTR<Object<PV3>> p;
    int jf, jg, jh;

    // if intersection Vert,
    // this Vert is inside c[0] of g,h encasement
    // this Vert is inside c[1] of f,h encasement
    // this Vert is inside c[2] of f,g encasement
    // if boundary Vert
    // this Vert touches leaf Cels c[0] and c[1] in f,g encasement
    Cel* c[3];

    // put surface polynomials f,g,h into fgh
    void getFGH(Encasement3D* e, Object<Poly3D>* fgh[3]);

    // store incident Cel into bounding vertex
    // use first available of c[0] and c[1]
    void setBCel(Cel* bc) {
      if (c[0] == 0)
        c[0] = bc;
      else if (c[1] == 0)
        c[1] = bc;
      else
        assert(0);
    }

    // disjoint-set tree
    // c[2] is unused by boundary Vert
    // but we will have to save jh because I am bloody-minded
    // ultimately, c[2] will point to the arrangement Edge containing
    // this (boundary) Vert
    int& size() { return jh; }
    Vert*& parent() { return *(Vert**)&c[2]; }
    Edge*& edge() { return *(Edge**)&c[2]; }
  };

  // coordinate with maximum difference
  // if xyzAvoid is 0,1 or 2, avoid it and also extreme coordinate of turning
  // points
  static int maxDiffCoord(Vert* u, Vert* v, int xyzAvoid = -1,
                          bool ffz_only = false);

  // bf[i], bf[j], bf[k]
  vector<PTR<Object<PV3>>> getRoots(int i, int j, int k);

  // bf[i] and two other polys
  vector<PTR<Object<PV3>>> getRoots(int i, PTR<Object<Poly3D>> f,
                                    PTR<Object<Poly3D>> g);

  // x=cx and bf[i] and bf[j]
  vector<PTR<Object<PV3>>> getRoots(double cx, int i, int j);

  // x=cx and bf[i] and g (actually fz)
  vector<PTR<Object<PV3>>> getRoots(double cx, int jf, PTR<Object<Poly3D>> g);

  // create bType intersections (see VertType)
  // x_{xyz}=m, f, and g
  static vector<PTR<Object<PV3>>> getRoots(PTR<Object<Poly3D>> f,
                                           PTR<Object<Poly3D>> g,
                                           VertType bType, int jf, int jg,
                                           int xyz, double m,
                                           const PV3<Parameter>& box);

  // fgh[ind(i, j, k)] (i<j<k) is set of intersections of bf[i], bf[j], bf[k]
  // fgh.size() ==  ind(0, 0, bf.size())
  vector<vector<PTR<Vert>>> fgh;

  void initFGH() {
    fgh.resize(ind(0, 0, bf.size()));
    for (int i = 0; i < bf.size(); ++i)
      for (int j = i + 1; j < bf.size(); ++j)
        for (int k = j + 1; k < bf.size(); ++k) {
          if ((i == 0 && j == 1) || (i == 2 && j == 3) || (i == 4 && j == 5))
            continue;
          if ((j == 2 && k == 3) || (j == 4 && k == 5)) continue;

          vector<PTR<Object<PV3>>> roots = getRoots(i, j, k);
          fgh[ind(i, j, k)].resize(roots.size());
          for (int h = 0; h < roots.size(); ++h)
            fgh[ind(i, j, k)][h] = new Vert(FGH, roots[h], i, j, k);
        }
  }

  // f,g encasement
  class FG : public RefCnt {
   public:
    FG(Encasement3D* e, int jf, int jg);
    FG(Encasement3D* e, int jf, int jg /* -1 for fz, js for slice */,
       PTR<Object<Poly3D>> g /* fz or slice */);
    Encasement3D* e;
    int jf, jg;  // if jg == -1, use g
    PTR<Object<Poly3D>> g;
    PTR<Object<Poly3D>> dFxdG[3];  // at most two, default 0 and 1 (x and y)
    vector<PTR<Vert>> fgd[3];      // f,g,dFxdG[i] intersections
    vector<PTR<Vert>> fgb;         // Cel boundary vertices
    vector<PTR<Object<Poly3D>>> planes;  // Cel splitting planes
    Cel top;                             // top of subdivision tree
    int nboxbverts;

    // add turning points and boundary Verts
    void addTB();

    int getAvoid() {
      int nz = (dFxdG[0] == 0) + (dFxdG[1] == 0) + (dFxdG[2] == 0);
      assert(nz == 3 || nz == 1);
      int xyzAvoid = -1;
      if (nz == 1) {
        xyzAvoid = dFxdG[0] == 0 ? 0 : dFxdG[1] == 0 ? 1 : 2;
        assert(dFxdG[xyzAvoid] == 0);
      }
      return xyzAvoid;
    }

    PTR<Object<Poly3D>> getF() { return e->bf[jf]; }
    PTR<Object<Poly3D>> getG() { return g == 0 ? e->bf[jg] : g; }

    // add set of intersection vertices
    void addIs(vector<PTR<Vert>>& verts) {
      PTR<Object<Poly3D>> f = getF();
      PTR<Object<Poly3D>> gg = getG();
      VertType bType = g == 0 ? FGB : jg == -1 ? FFzB : FSB;
      for (auto iv = verts.begin(); iv != verts.end(); ++iv)
        top.addI(f, gg, bType, jf, jg, *iv, fgb, getAvoid(), e->box);
    }

    // inform Verts of the Cels they touch
    void setVertCels() {
      top.setVertCels(e, getF(), getG());
      // get rid of box boundary boundary points
      fgb.erase(fgb.begin(), fgb.begin() + nboxbverts);
    }

    vector<PTR<Edge>> edges;

    // find f,g arrangement Edges
    void findEdges();

    // f(p)=g(p)=0.  get arrangement Edge that contains p
    Edge* getEdge(PTR<Object<PV3>> p) { return top.getEdge(p, getAvoid()); }
  };

  vector<PTR<FG>> fg;  // all f,g encasements
  int ind(int i, int j) {
    assert(!((j < 6 && i / 2 == j / 2)));
    return i + j * (j - 1) / 2;
  }

  void initFG() {
    fg.resize(ind(0, bf.size()));
    for (int i = 0; i < bf.size(); ++i)
      for (int j = i + 1; j < bf.size(); ++j) {
        if (j < 6 && i / 2 == j / 2) continue;
        fg[ind(i, j)] = new FG(this, i, j);
      }
  }

  // directed f,g edge
  class DEdge {
    int fgfbnp;

   public:
    // fg: edge as element of f or g
    // fb: forward (v[0] to v[1]) or backward
    // np: the cell with f<0 or f>0 (g<0 or g>0) for the face to the left
    // --left as viewed from the positive side of f (or g)
    DEdge() : e(0), fgfbnp(0) {}
    DEdge(Edge* e, int fg, int fb, int np = 0)
        : e(e), fgfbnp(fg + 2 * fb + 4 * np) {}
    Edge* e;

    bool operator!=(const DEdge& d) const {
      return e != d.e || fgfbnp != d.fgfbnp;
    }
    bool isNull() { return e == 0; }
    int fg() { return fgfbnp % 2; }
    int fb() { return fgfbnp / 2 % 2; }
    int np() { return fgfbnp / 4; }
    Vert* tail() { return e->v[fb()]; }
    Vert* head() { return e->v[!fb()]; }
    DEdge next();
    Loop*& loop();
    DEdge twin() { return DEdge(e, fg(), !fb(), np()); }
    DEdge twin(int fb) { return DEdge(e, fg(), fb, np()); }
    /*  fg  fb  np      fg  fb  np
         0   0   0       1   1   0
         0   0   1       1   0   0
         0   1   0       1   1   1
         0   1   1       1   0   1  */
    DEdge sameShell(int np) {
      // Have to check both ways because it could be a boundary edge
      // with no shell in one direction.
      if (e->next[!fg()][0].isNull() && e->next[!fg()][1].isNull())
        // FS or FFz edge
        return twin();
      else
        return DEdge(e, !fg(), !fg() ^ np, fg() ^ fb());
    }
  };

  class Loop : public RefCnt {
   public:
    Loop(DEdge d) : dedge(d), next(this), comp(-1), face(0) {}
    DEdge dedge;
    Loop* next;  // circular linked list of loops of a face
    Face* face;
    bool sameFace(Loop* that) {
      Loop* loop = this;
      do {
        if (loop == that) return true;
        loop = loop->next;
      } while (loop != this);
      return false;
    }
    void joinFace(Loop* that) {
      Loop* thisNext = this->next;
      Loop* thatNext = that->next;
      this->next = thatNext;
      that->next = thisNext;
    }
    void setFace(Face* face) {
      Loop* loop = this;
      do {
        loop->face = face;
        loop = loop->next;
      } while (loop != this);
    }

    int comp;

    // set the connected component number of this Loop
    // get the farthest Verts in the xyz direction
    void setComp(int comp, Vert*& vxmin, Vert*& vxmax, int xyz) {
      this->comp = comp;
      DEdge d = dedge;
      do {
        if (vxmin == 0)
          vxmin = vxmax = d.tail();
        else if (V3minusV3I(d.tail()->p, vxmin->p, xyz) < 0)
          vxmin = d.tail();
        else if (V3minusV3I(d.tail()->p, vxmax->p, xyz) > 0)
          vxmax = d.tail();
        if (d.twin().loop() == 0)
          ;  // boundary edge, no loop on other side
        else if (d.twin().loop()->comp == -1)
          d.twin().loop()->setComp(comp, vxmin, vxmax, xyz);
        else
          assert(d.twin().loop()->comp == comp);
        d = d.next();
      } while (d != dedge);
    }
  };

  // information about a boundary or input surface
  class F : public RefCnt {
   public:
    F(Encasement3D* e, int jf);
    Encasement3D* e;
    int jf;
    PTR<Object<Poly3D>> fz;

    // ffzh[ih] are the intersections of f,fz, and bf[h]
    vector<vector<PTR<Vert>>> ffzh;

    PTR<FG> ffz;  // f,fz encasement

    // slices to cut holes in f=0
    vector<PTR<Object<Poly3D>>> slices;

    // fsh[js][jh] are intersections of f, slices[js], and bf[jh]
    vector<vector<vector<PTR<Vert>>>> fsh;

    // fsfz[js] are the intersections of f, slices[js], and fz
    // vector<vector<PTR<Vert>>> fsfz;  // DEPRECATED

    // fs[js] is the f,slices[js] encasement
    vector<PTR<FG>> fs;

    vector<DEdge> dedges;  // one DEdge per edge with fg set

    // set the next fields of each Edge
    void setEdgeNexts();

    vector<PTR<Loop>> loops;
    vector<PTR<Face>> faces;
    void findLoops();
    void findFaces();

    // each Shell involves multiple Fs, but it is a convenient method
    // to find the Shells associated with the Faces of this F
    void findShells();

    // get a point that lies on Face face of this F
    PTR<Object<PV3>> getPointOnFace(Face* face);

    // f(p)=0.  Find the Face that contains p.
    Face* getFace(PTR<Object<PV3>> p);
  };

  vector<PTR<F>> f;  // all F information

  void initF() {
    f.resize(bf.size());
    for (int i = 0; i < bf.size(); ++i) /*f[i] = */
      new F(this, i);
  }

  // add intersections to all top level FGs
  // FGH (to fg, fh, and gh), FFzH (to fh), FSH (to fh)
  void addIs();

  // get x values hitting all the holes in f
  vector<double> getSlices(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> fz);

  class Edge : public RefCnt {
    // get the next DEdge around the loop to the left
    // viewing from above
    // fg: f or g
    // fb: forward or backward
    DEdge getNext(int fg, int fb) {
      Object<Poly3D>* fgb[3];
      b->getFGH(e, fgb);  // get f and g for the edge
      Object<Poly3D>* fgh[3];
      v[!fb]->getFGH(e, fgh);  // get f, g, and h for the head

      // This is the case that we're on a boundary/g edge tracing out
      // a face of g that winds out of the box. Don't have any next
      // so that no loop can be formed.
      if (b->jf < 6 && b->jg >= 6 && fg == 1 && fb == 0) {
        return DEdge();
      }

      Cel* c = nullptr;

      // See https://github.com/joemasterjohn/acp/issues/2
      // IF we're at a turning point
      //   IF there's a non-null c[0] (a dual vertex)
      //       AND [(we're tracing FG (Boundary / Surface) and fg == 1 (Surface)
      //            OR we're tracing FS (Surface / Slice) and fg == 0 (Surface))
      //       OR we're tracing FFz]
      //
      //     set c properly (c[0] or c[2])
      //
      //   ELSE we're either a normal turning point OR we're tracing a loop on a
      //        boundary/slice face and hit turning point
      //
      //     keep going forward

      // If we're at a FGB intersection, F is the boundary and G is the non
      // boundary surface because boundaries have lower indices than
      // non-boundaries. If fg == 0 in that case, we should cruise through the
      // turning point, because we're tracing a face on the boundary surface.

      // If we're at an FSB intersection, F is the non slice and G is the slice.
      // If fg == 1 (can this ever happen?) we're tracing on the slice (??)
      // We should never encounter an edge with FSB endpoint and fg == 1, right?

      if (v[!fb]->type % 3 == 1) {  // turning point

        // b is an internal boundary point on the edge, and an XXB vertex

        // Edge traces a face of non-boundary or
        // Edge traces a face of non-slice
        bool nonFlatSurface =
            (b->type == FGB && fg == 1) || (b->type == FSB & fg == 0);

        // Edge traces a face of F
        bool ffzSurface = b->type == FFzB;

        // Dual vertex, AND (we're tracing FFz around F OR
        //   we're tracing FB or FS around the face of a non-boundary/non-slice)
        if (v[!fb]->c[0] != nullptr && (nonFlatSurface || ffzSurface)) {
          if (nonFlatSurface) {
            c = v[!fb]->c[0];
          } else if (ffzSurface) {
            c = v[!fb]->c[2];
          }
        } else {
          // Just a normal old turning point, continue to the next edge
          Edge* edgeL = v[!fb]->c[2]->lb.b->edge();
          Edge* edgeR = v[!fb]->c[2]->rb.b->edge();
          assert(this == edgeL || this == edgeR);
          DEdge n = this == edgeR ? DEdge(edgeL, fg, fb) : DEdge(edgeR, fg, fb);
          assert(n.tail() == v[!fb]);
          return n;
        }

      } else {
        int i = 0;
        for (; i < 3; i++) {
          if (fgh[i] == fgb[!fg]) {
            break;  // i is index of surface "g" (fg) in fgh
          }
        }

        // alternatively
        // int i = static_cast<int>(std::find(fgh, fgh+3, fgb[!fg]) - fgh);

        c = v[!fb]->c[i];  // c[i] is in "f"h arrangement
      }

      if (c == 0) {
        // No encasement cell, so cannot turn.
        // Go straight.
        // get the index in fgh that is not f or g
        int notForG = 0;
        for (; notForG < 3; ++notForG)
          if (fgh[notForG] != fgb[0] && fgh[notForG] != fgb[1]) break;

        // node is c[notForG]
        Edge* edgeL = v[!fb]->c[notForG]->lb.b->edge();
        Edge* edgeR = v[!fb]->c[notForG]->rb.b->edge();
        assert(this == edgeL || this == edgeR);
        DEdge n = this == edgeR ? DEdge(edgeL, fg, fb) : DEdge(edgeR, fg, fb);
        assert(n.tail() == v[!fb]);
        return n;
      }

      int s = 2 * (fg ^ fb) - 1;  // sign of "g" for left turn

      assert(c->lb.b != 0);
      Object<Poly3D>* fgb2[3];
      c->lb.b->getFGH(e, fgb2);  // ???

      // debug
      // confirm that c is in the "f"h encasement
      assert(fgb2[0] != fgb[!fg] && fgb2[1] != fgb[!fg]);
      assert(fgb[fg] == fgb2[0] || fgb[fg] == fgb2[1]);

      Edge* edge = 0;
      // look for sign of "g" equal s
      if (Poly3Dval(fgb[!fg], c->lb.b->p) == s)
        edge = c->lb.b->edge();
      else if (c->rb.b != 0)
        edge = c->rb.b->edge();
      else
        return DEdge();

      return DEdge(edge, fgb[fg] == fgb2[1], v[!fb] == edge->v[fb] ? fb : !fb);
    }

   public:
    Edge(Encasement3D* e, Vert* b) : e(e), b(b) {
      v[0] = v[1] = 0;
      loop[0][0] = loop[0][1] = loop[1][0] = loop[1][1] = 0;
    }
    Encasement3D* e;
    Vert *b, *v[2];
    DEdge next[2 /*f g*/][2 /*forward backward*/];
    Loop* loop[2 /*f g*/][2 /*forward backward*/];
    void setNexts(int fg) {
      for (int fb = 0; fb < 2; ++fb) {
        next[fg][fb] = getNext(fg, fb);
      }
    }
  };

  void setVertCels();
  void initEdges();
  void findLoops();

  // mesh
  // v is edge->v[iv] for an Edge of the FXEncasement
  // v will lie on an edge of the Encasement3D
  DEdge getDEdge(FXEncasement::Vert* v, PTR<Object<Scalar>> x, int jf, int iv);

  // v is edge->v[iv] for an Edge of the FXEncasement
  // v will lie on an edge of the Encasement3D
  Loop* getLoop(FXEncasement::Vert* v, PTR<Object<Scalar>> x, int jf, int iv) {
    return getDEdge(v, x, jf, iv).loop();
  }

  // v is the intersection of x_{ind}=d with a curve on a face of the box
  // We consider pairs of consecutive vertices.
  // iv=0 means lower, iv=1 means upper
  // lu = 0, 1 means lb(), ub() face.
  // return the Loop that is to that side of the curve (above for iv=0, below
  // for iv=1)
  Loop* getLoop(Vert* v, int iv, int xind, int lu);

  class Face : public RefCnt {
   public:
    Face(Encasement3D* enc, int jf, Loop* loop) : enc(enc), jf(jf), loop(loop) {
      shell[0] = shell[1] = 0;
    }
    Encasement3D* enc;
    int jf;
    Loop* loop;
    Shell* shell[2];
    void setShell(int np, Shell* shell);

    // false if shell would be outside of box
    bool hasShell(int np) { return !loop->dedge.sameShell(np).next().isNull(); }
  };

  PTR<Object<PV3>> getPointOnFace(Face* face) {
    return f[face->jf]->getPointOnFace(face);
  }

  class Compare {
    int xyz;

   public:
    Compare(int xyz) : xyz(xyz) {}
    bool operator()(Vert* u, Vert* v) const {
      return V3minusV3I(u->p, v->p, xyz) < 0;
    }
  };

  class XCompare {
   public:
    bool operator()(pair<Object<PV3>*, bool> u,
                    pair<Object<PV3>*, bool> v) const {
      return V3minusV3I(u.first, v.first, 0) < 0;
    }
  };

  void findFaces() {
    for (int jf = 0; jf < bf.size(); jf++) f[jf]->findFaces();
  }

  class Shell : public RefCnt {
   public:
    Shell(Encasement3D* enc) : enc(enc), next(this), cell(0), comp(-1) {}
    Encasement3D* enc;
    vector<Face*> faces;
    Shell* next;
    Cell* cell;
    bool sameCell(Shell* that) {
      Shell* shell = this;
      do {
        if (shell == that) return true;
        shell = shell->next;
      } while (shell != this);
      return false;
    }
    void joinCell(Shell* that) {
      Shell* thisNext = this->next;
      Shell* thatNext = that->next;
      this->next = thatNext;
      that->next = thisNext;
    }
    void setCell(Cell* cell) {
      Shell* shell = this;
      do {
        shell->cell = cell;
        cell->shells.push_back(shell);
        shell = shell->next;
      } while (shell != this);
    }

    int comp;
    void setComp(int comp) {
      this->comp = comp;
      for (Face* face : faces) {
        int ishell = this == face->shell[1];
        assert(face->shell[ishell] == this);
        if (face->shell[!ishell] == 0)
          ;  // box boundary face
        else if (face->shell[!ishell]->comp == -1)
          face->shell[!ishell]->setComp(comp);
        else
          assert(face->shell[!ishell]->comp == comp);
      }
    }
  };

  vector<PTR<Shell>> shells;

  void findShells() {
    for (int jf = 0; jf < bf.size(); ++jf) f[jf]->findShells();
  }

  class Cell : public RefCnt {
   public:
    Cell(Encasement3D* enc) : enc(enc) {}
    Encasement3D* enc;
    vector<Shell*> shells;
  };

  vector<PTR<Cell>> cells;

  void findCells();

  Cell* getCell(PTR<Object<PV3>> p);
};

FXEncasement::FXEncasement(Encasement3D* e3, int fIndex, double xValue,
                           PV2<Parameter> box)
    : f3(e3->bf[fIndex]),
      e3(e3),
      fIndex(fIndex),
      xValue(xValue),
      box(box),
      bf(e3->bf.size()) {
  xObject = new Object<Scalar>(Parameter::constant(xValue));
  // e3->bf[0,1] are the x-bounding surfaces and are parallel to x=xValue
  for (int i = 2; i < bf.size(); ++i)
    bf[i] = new Sub1dPoly3D(e3->bf[i], 0, xObject);
  f = bf[fIndex];
  for (int yz = 0; yz < 2; ++yz) fyz[yz] = new P2derI(f, yz);
  initVerts();
  initCels();
  findEdges();
}

// mesh
FXEncasement::FXEncasement(Encasement3D* e3, int fIndex, PTR<Object<PV3>> v3,
                           Vert3DInfo v3info)
    : f3(e3->bf[fIndex]),
      e3(e3),
      fIndex(fIndex),
      xValue(0),
      box(e3->box[1], e3->box[2]),
      bf(e3->bf.size()),
      v3(v3),
      v3info(v3info) {
  v3p = v3->getApprox();
  xObject = new V3I(v3, 0);
  for (int i = 2; i < bf.size(); ++i)
    bf[i] = new Sub1dPoly3D(e3->bf[i], 0, xObject);
  f = bf[fIndex];
  for (int yz = 0; yz < 2; ++yz) fyz[yz] = new P2derI(f, yz);
  initVerts();
  initCels();
  findEdges();
}

int maxDiffCoord(FXEncasement::Vert* u, FXEncasement::Vert* v,
                 bool checkForTurning = true) {
  if (checkForTurning) {
    // Remember, we are in yz plane so yz=1 is really the z direction.
    // Not monotonic in y or z if it contains f(p)=fyz[1 or 0](p) = 0
    if (u->type == FXEncasement::FD) return u->jg;
    if (v->type == FXEncasement::FD) return v->jg;
  }

  FXEncasement::Vert* vymax = V2minusV2I(u->p, v->p, 0) > 0 ? u : v;
  FXEncasement::Vert* vymin = v == vymax ? u : v;
  FXEncasement::Vert* vzmax = V2minusV2I(u->p, v->p, 1) > 0 ? u : v;
  FXEncasement::Vert* vzmin = v == vzmax ? u : v;

  return XmXmYmY(vymax->p, vymin->p, vzmax->p, vzmin->p) > 0 ? 0 : 1;
}

// return a double value between the yz coordinates of u and v
// try a cheap approximation first
double getMiddle(FXEncasement::Vert* u, FXEncasement::Vert* v, int yz) {
  assert(V2minusV2I(u->p, v->p, yz) < 0);
  double m =
      (u->p->getApprox(1.0)[yz].ub() + v->p->getApprox(1.0)[yz].lb()) / 2;
  if (V2IminusD(u->p, yz, m) < 0 && V2IminusD(v->p, yz, m) > 0) return m;
  m = (u->p->getApprox()[yz].ub() + v->p->getApprox()[yz].lb()) / 2;
  assert(V2IminusD(u->p, yz, m) < 0 && V2IminusD(v->p, yz, m) > 0);
  return m;
}

// add a f=fyz[i]=0 vertex to the tree of Cels
void FXEncasement::Cel::addT(Vert* v) {
  // internal split on yz
  if (0 <= yz && yz < 2) {
    if (V2IminusD(v->p, yz, spv.s) < 0)
      lb.c->addT(v);
    else
      rb.c->addT(v);
    return;
  }

  // no splitting by lines yet
  assert(yz != 2);

  // empty leaf Cel
  if (spv.v == 0) {
    spv.v = v;
    return;
  }

  // mesh
  if (v == spv.v->twin) return;

  // two vertices in leaf
  // split current Cel between verts and insert them into two new leaf Cels
  yz = maxDiffCoord(v, spv.v, false);
  Vert* vL = V2minusV2I(v->p, spv.v->p, yz) < 0 ? v : spv.v;
  Vert* vR = vL == spv.v ? v : spv.v;
  spv.s = getMiddle(vL, vR, yz);
  lb.c = new Cel();
  rb.c = new Cel();
  lb.c->spv.v = vL;
  rb.c->spv.v = vR;
}

// x_{xy} < m half of box
PV2<Parameter> leftBox(PV2<Parameter> box, int xy, double m) {
  box[xy] = Parameter::constant(box[xy].lb()).interval(Parameter::constant(m));
  return box;
}

// x_{xy} > m half of box
PV2<Parameter> rightBox(PV2<Parameter> box, int xy, double m) {
  box[xy] = Parameter::constant(m).interval(Parameter::constant(box[xy].ub()));
  return box;
}

PTR<Object<Poly2D>> separatingLine(PTR<Object<Poly2D>> f,
                                   vector<FXEncasement::Vert*>& bverts,
                                   const PV2<Parameter>& box) {
  return 0;
}

// largest coordinate of box
int maxCoord(const PV2<Parameter>& box) {
  return box[0].ub() - box[0].lb() < box[1].ub() - box[1].lb();
}

// add boundary vertices
void FXEncasement::Cel::addB(FXEncasement* e, vector<Vert*>& bverts,
                             const PV2<Parameter>& box) {
  // mesh
  cbox = box;

  if (0 <= yz && yz < 2) {
    vector<Vert*> bvertsL, bvertsR;
    for (int h = 0; h < bverts.size(); ++h)
      if (V2IminusD(bverts[h]->p, yz, spv.s) < 0)
        bvertsL.push_back(bverts[h]);
      else
        bvertsR.push_back(bverts[h]);

    vector<PTR<Object<PV2>>> broots = e->getRoots(yz, spv.s, box);
    for (int h = 0; h < broots.size(); ++h) {
      Vert* bvert = new Vert(FB, broots[h], yz);
      e->fb.push_back(bvert);  // global set of boundary vertices
      bvertsL.push_back(bvert);
      bvertsR.push_back(bvert);
    }

    lb.c->addB(e, bvertsL, leftBox(box, yz, spv.s));
    rb.c->addB(e, bvertsR, rightBox(box, yz, spv.s));
    return;
  }

  assert(yz != 2);

  if (bverts.size() <= 2) {
    setbverts(bverts);
    return;
  }

  // mesh
  bool containsFG = false;
  if (spv.v != 0 && spv.v->twin != 0 && bverts.size() == 4) {
    for (int ig = 0; !containsFG && ig < e->bf.size(); ++ig)
      for (Vert* v : e->fg[ig])
        if (contains(box, v->p)) {
          containsFG = true;
          break;
        }
    if (!containsFG) {
      lb.b = 0;
      rb.b = 0;
      for (Vert* b : bverts)
        if (P2derIV2(e->f, 1, b->p) > 0) {
          if (lb.b == 0)
            lb.b = b;
          else {
            assert(rb.b == 0);
            rb.b = b;
          }
        }
      for (Vert* b : bverts)
        if (b != lb.b && b != rb.b) {
          if (lb.b->twin == 0)
            lb.b->twin = b;
          else
            rb.b->twin = b;
        }
      return;
    }
  }

  // take the vertex out of the Cel because we will using the spv
  // union to store the splitting line
  Vert* v = spv.v;
  spv.v = 0;

  if (bverts.size() <= 4 && !containsFG) {
    // for splitting line purposes, an odd number of boundary points
    // means an internal vertex that is really a boundary point
    vector<Vert*> bvertsSL = bverts;
    if (bvertsSL.size() == 3) {
      assert(v != 0);
      bvertsSL.push_back(v);
    }

    PTR<Object<Poly2D>> p = separatingLine(e->f, bvertsSL, box);
    if (p != 0) {             // found a splitting line!
      e->lines.push_back(p);  // save it in global list
      yz = 2;
      spv.p = p;  // this Cel now split by line

      vector<Vert*> bvertsL, bvertsR;
      for (int h = 0; h < bverts.size(); ++h)
        if (Poly2Dval(spv.p, bverts[h]->p) < 0)
          bvertsL.push_back(bverts[h]);
        else
          bvertsR.push_back(bverts[h]);

      assert(bvertsL.size() <= 2);
      assert(bvertsR.size() <= 2);

      lb.c = new Cel();
      rb.c = new Cel();

      // move interval vertex down to correct new leaf
      if (v != 0) {
        if (Poly2Dval(spv.p, v->p) < 0)
          lb.c->spv.v = v;
        else
          rb.c->spv.v = v;
      }

      // no need to recurse
      // just put boundary vertices into leaf Cels
      lb.c->setbverts(bvertsL);
      rb.c->setbverts(bvertsR);
      return;
    }
  }

  // more than 4 boundary vertices or line could not be found
  // split the box down the middle and recurse
  yz = maxCoord(box);
  spv.s = (box[yz].lb() + box[yz].ub()) / 2;
  vector<Vert*> bvertsL, bvertsR;
  for (int h = 0; h < bverts.size(); ++h)
    if (V2IminusD(bverts[h]->p, yz, spv.s) < 0)
      bvertsL.push_back(bverts[h]);
    else
      bvertsR.push_back(bverts[h]);

  vector<PTR<Object<PV2>>> broots = e->getRoots(yz, spv.s, box);
  for (int h = 0; h < broots.size(); ++h) {
    Vert* bvert = new Vert(FB, broots[h], yz);
    e->fb.push_back(bvert);
    bvertsL.push_back(bvert);
    bvertsR.push_back(bvert);
  }

  lb.c = new Cel();
  rb.c = new Cel();

  if (v != 0) {
    if (V2IminusD(v->p, yz, spv.s) < 0)
      lb.c->spv.v = v;
    else
      rb.c->spv.v = v;
  }

  lb.c->addB(e, bvertsL, leftBox(box, yz, spv.s));
  rb.c->addB(e, bvertsR, rightBox(box, yz, spv.s));
}

// add intersection f=g=0 to subtree
// p and side are set if ancestor was split by line
void FXEncasement::Cel::addI(FXEncasement* e, Vert* v,
                             const PV2<Parameter>& box, Object<Poly2D>* p,
                             int side) {
  // mesh
  cbox = box;

  if (0 <= yz && yz < 2) {
    if (V2IminusD(v->p, yz, spv.s) < 0)
      lb.c->addI(e, v, leftBox(box, yz, spv.s), p, side);
    else
      rb.c->addI(e, v, rightBox(box, yz, spv.s), p, side);
    return;
  }

  if (yz == 2) {
    assert(p == 0);  // only one line split along each path to leaf
    if (Poly2Dval(spv.p, v->p) < 0)
      lb.c->addI(e, v, box, spv.p, -1);  // everything on -1 side
    else
      rb.c->addI(e, v, box, spv.p, 1);  // everything to +1 side
    return;
  }

  if (spv.v == 0) {
    spv.v = v;
    return;
  }

  // mesh
  if (v == spv.v->twin) return;

  // two vertices in cell, split between them
  yz = maxDiffCoord(v, spv.v);
  Vert* vL = V2minusV2I(v->p, spv.v->p, yz) < 0 ? v : spv.v;
  Vert* vR = vL == spv.v ? v : spv.v;
  spv.s = getMiddle(vL, vR, yz);

  Vert* bverts[2] = {lb.b, rb.b};
  lb.c = new Cel();
  rb.c = new Cel();

  // debug
  lb.c->cbox = leftBox(box, yz, spv.s);
  rb.c->cbox = rightBox(box, yz, spv.s);

  lb.c->spv.v = vL;  // put interval verts
  rb.c->spv.v = vR;  // into new leaf Cels
  vector<Vert*> bvertsL, bvertsR;
  for (int h = 0; h < 2; ++h)
    if (bverts[h] != 0) {
      if (V2IminusD(bverts[h]->p, yz, spv.s) < 0)
        bvertsL.push_back(bverts[h]);
      else
        bvertsR.push_back(bverts[h]);
    }

  vector<PTR<Object<PV2>>> broots = e->getRoots(yz, spv.s, box);
  for (int h = 0; h < broots.size(); ++h)
    // only consider split vertices on correct side of line
    if (p == 0 || Poly2Dval(p, broots[h]) == side) {
      Vert* bvert = new Vert(FB, broots[h], yz);
      e->fb.push_back(bvert);
      bvertsL.push_back(bvert);
      bvertsR.push_back(bvert);
    }

  // this will set off an alarm if there were more than two vertices
  // in the middle hence more than two into a leaf
  lb.c->setbverts(bvertsL);
  rb.c->setbverts(bvertsR);
}

// inform all vertices of the leaf cells they touch
void FXEncasement::Cel::setVertCels(FXEncasement* e) {
  if (0 <= yz && yz <= 2) {
    lb.c->setVertCels(e);
    rb.c->setVertCels(e);
    return;
  }

  if (spv.v != 0) {
    spv.v->c[0] = this;

    // mesh
    if (spv.v->twin != 0) spv.v->twin->c[0] = this;
  }

  if (lb.b != 0) {
    if (lb.b->boxvert)  // box boundary point
      lb.b = 0;         // remove it
    else {
      lb.b->setBCel(this);

      // mesh
      if (spv.v != 0 && spv.v->twin != 0 && spv.v->type == FD &&
          spv.v->twin->type == FD && lb.b->twin != 0)
        lb.b->twin->setBCel(this);
    }
  }
  if (rb.b != 0) {
    if (rb.b->boxvert)  // box boundary point
      rb.b = 0;         // remove it
    else {
      rb.b->setBCel(this);

      // mesh
      if (spv.v != 0 && spv.v->twin != 0 && spv.v->type == FD &&
          spv.v->twin->type == FD && rb.b->twin != 0)
        rb.b->twin->setBCel(this);
    }
  }
  if (lb.b == 0 && rb.b != 0) {
    lb.b = rb.b;
    rb.b = 0;
  }
}

// find all the arrangement edges
void FXEncasement::findEdges() {
  // save the jg field of each boundary vertex
  vector<int> yzs(fb.size());
  for (int ib = 0; ib < fb.size(); ++ib) yzs[ib] = fb[ib]->jg;

  for (Vert* b : fb) MakeSet<Vert>(b);
  // for each boundary vertex
  for (Vert* b : fb)
    // for each cell it touches
    for (int ic = 0; ic < 2; ++ic)
      // if the cell is empty, meaning it is spanned by a single edge
      if (b->c[ic]->spv.v == 0) {
        Cel* c = b->c[ic];
        if (c->lb.b != 0 && c->rb.b != 0)
          // its two boundary vertices belong to the same edge
          Union<Vert>(c->lb.b, c->rb.b);
      }
  // make sure each Vert is a root of a disjoint-set tree or one below
  // the root
  for (Vert* b : fb) Find<Vert>(b);
  for (Vert* b : fb)
    // each root is a new Edge
    if (b->parent() == b) {
      b->edge() = new Edge(b);
      edges.push_back(b->edge());  // save in global list
      b->size() = 0;               // mark them as visited
    }

  for (Vert* b : fb)
    if (b->size() != 0) {  // not visited
      b->edge() = b->parent()->edge();
      b->size() = 0;
    }

  // put back jg (actually yz) information
  for (int ib = 0; ib < fb.size(); ++ib) fb[ib]->jg = yzs[ib];

  // set v[0,1] tails and head of edges
  for (Vert* b : fb)
    for (int ic = 0; ic < 2; ++ic)
      // b incident on non-empty cell
      // its internal vertex is an Edge endpoint
      if (b->c[ic]->spv.v != 0) {
        Vert* v = b->c[ic]->spv.v;
        // set edge->v[0,1] from left to right (increasing y actually)
        // if fz is positive
        int svb = V2minusV2I(b->p, v->p, 0);
        int sfz = P2derIV2(f, 1, b->p);
#ifdef BLEEN
        if (svb * sfz > 0) {
          assert(b->edge()->v[0] == 0);
          b->edge()->v[0] = v;
          v->e[0] = b->edge();
        } else {
          assert(b->edge()->v[1] == 0);
          b->edge()->v[1] = v;
          v->e[1] = b->edge();
        }
#endif
        int vInd = !(svb * sfz > 0);

        assert(b->edge()->v[vInd] == 0);

        // mesh
        if (v->twin != 0 && v->type == FD && v->twin->type == FD) {
          if (b->twin != 0) {
            b->edge()->v[vInd] = v;
            v->e[vInd] = b->edge();
          } else {
            b->edge()->v[vInd] = v->twin;
            v->twin->e[vInd] = b->edge();
          }
        } else

        {
          b->edge()->v[vInd] = v;
          v->e[vInd] = b->edge();
        }
      }
}

// f(p)=0.  find Edge that contains p.
FXEncasement::Edge* FXEncasement::Cel::getEdge(PTR<Object<PV2>> p) {
  if (0 <= yz && yz < 2) {
    if (V2IminusD(p, yz, spv.s) < 0)
      return lb.c->getEdge(p);
    else
      return rb.c->getEdge(p);
  }

  if (yz == 2) {
    if (Poly2Dval(spv.p, p) < 0)
      return lb.c->getEdge(p);
    else
      return rb.c->getEdge(p);
  }

  // only one Edge in cell
  if (spv.v == 0) return lb.b->edge();

  // each half-edge is monotonic in yzmon.
  int yzmon = spv.v->type == FD ? spv.v->jg : 0;
  // pick the Edge of the boundary vertex to the same side of the
  // internal vertex
  int s1 = V2minusV2I(p, spv.v->p, yzmon);        // side of p
  int s2 = V2minusV2I(lb.b->p, spv.v->p, yzmon);  // side of lb.b

  // debug
  if (rb.b != 0)
    assert(V2minusV2I(rb.b->p, spv.v->p, yzmon) == -s2);  // side of lb.b

  if (s1 == s2)
    return lb.b->edge();
  else
    return rb.b->edge();
}

// Puts together edges with with a common x turning point (fx = 0)
void FXEncasement::expand(Edge* e, Vert* v[2]) {
  Edge *first = e, *last = e;
  while (first->v[0]->type == FD && first->v[0]->jg == 0)
    first = first->v[0]->e[1];
  while (last->v[1]->type == FD && last->v[1]->jg == 0) last = last->v[1]->e[0];
  v[0] = first->v[0];
  v[1] = last->v[1];
}

Encasement3D::Encasement3D(const PV3<Parameter>& box,
                           vector<PTR<Object<Poly3D>>> polys) {
  for (int xyz = 0; xyz < 3; ++xyz)
    for (int lu = 0; lu < 2; ++lu)
      bf.push_back(new Object<Poly3D>(boxPoly(box, xyz, lu)));
  for (auto it = polys.begin(); it != polys.end(); ++it) bf.push_back(*it);

  this->box = box;

  Timer timer;
  timer.Mark("start");

  initFGH();
  timer.Mark("initFGH()");

  initFG();
  timer.Mark("initFG()");

  initF();
  timer.Mark("initF()");

  addIs();
  timer.Mark("addIs()");

  setVertCels();
  timer.Mark("setVertCels()");

  initEdges();
  timer.Mark("initEdges()");

  findLoops();
  timer.Mark("findLoops()");

  findFaces();
  timer.Mark("findFaces()");

  findShells();
  timer.Mark("findShells()");

  findCells();
  timer.Mark("findCells()");

  timer.Report();
}

vector<PTR<Object<PV3>>> Encasement3D::getRoots(int i, int j, int k) {
  vector<PTR<Object<Poly3D>>> fs;
  vector<int> xyzs;
  vector<double> ds;

  int is[3] = {i, j, k};
  for (int h = 0; h < 3; ++h)
    if (is[h] >= 6)
      fs.push_back(bf[is[h]]);
    else {
      int xyz = is[h] / 2;
      int lu = is[h] % 2;
      double d = lu == 0 ? box[xyz].lb() : box[xyz].ub();
      xyzs.push_back(xyz);
      ds.push_back(d);
    }

  vector<PTR<Object<PV3>>> roots = ::getRoots(fs, xyzs, ds, box);

  return roots;
}

vector<PTR<Object<PV3>>> Encasement3D::getRoots(int i, PTR<Object<Poly3D>> f,
                                                PTR<Object<Poly3D>> g) {
  if (i >= 6) return ::getRoots(bf[i], f, g, box);

  int xyz = i / 2;
  int lu = i % 2;
  double d = lu == 0 ? box[xyz].lb() : box[xyz].ub();

  vector<PTR<Object<Poly3D>>> fs = {f, g};
  vector<int> xyzs = {xyz};
  vector<double> ds = {d};

  vector<PTR<Object<PV3>>> roots = ::getRoots(fs, xyzs, ds, box);

  return roots;
}

vector<PTR<Object<PV3>>> Encasement3D::getRoots(double x, int i, int j) {
  assert(j >= 6);
  vector<PTR<Object<Poly3D>>> fs = {bf[j]};
  vector<int> xyzs = {0};
  vector<double> ds = {x};

  if (i >= 6)
    fs.push_back(bf[i]);
  else {
    int xyz = i / 2;
    int lu = i % 2;
    double d = lu == 0 ? box[xyz].lb() : box[xyz].ub();
    xyzs.push_back(xyz);
    ds.push_back(d);
  }

  vector<PTR<Object<PV3>>> roots = ::getRoots(fs, xyzs, ds, box);

  return roots;
}

vector<PTR<Object<PV3>>> Encasement3D::getRoots(double cx, int jf,
                                                PTR<Object<Poly3D>> g) {
  assert(jf >= 6);
  vector<PTR<Object<Poly3D>>> fs = {bf[jf], g};
  vector<int> xyzs = {0};
  vector<double> ds = {cx};

  vector<PTR<Object<PV3>>> roots = ::getRoots(fs, xyzs, ds, box);

  return roots;
}

Encasement3D::FG::FG(Encasement3D* e, int jf, int jg) : e(e), jf(jf), jg(jg) {
  int nt = 0;  // only add turning points for two directions
  if (jg >= 6) {
    for (int xyz = 0; xyz < 3 && nt < 2; ++xyz)
      if (xyz != jf / 2) {
        dFxdG[xyz] = new DFxDGI(e->bf[jf], e->bf[jg], xyz);
        ++nt;
      }
    assert(nt == 2);
  }

  for (int xyz = 0; xyz < 3; ++xyz)
    if (dFxdG[xyz] != 0) {
      vector<PTR<Object<PV3>>> roots = e->getRoots(jf, e->bf[jg], dFxdG[xyz]);
      for (auto r : roots) {
        Root3D* r3 = dynamic_cast<Root3D*>(r.t);
        if (r3 != nullptr) {
          printf("box(%8lf %8lf %8lf)\n", r3->box.x.mid(), r3->box.y.mid(),
                 r3->box.z.mid());
        }
      }
      fgd[xyz].resize(roots.size());
      for (int h = 0; h < roots.size(); ++h)
        fgd[xyz][h] = new Vert(FGD, roots[h], jf, jg, xyz);
    }
  addTB();
}

Encasement3D::FG::FG(Encasement3D* e, int jf,
                     int jg /* -1 for fz, js for slice */,
                     PTR<Object<Poly3D>> g /* fz or slice */)
    : e(e), jf(jf), jg(jg), g(g) {
  assert(jf >= 6);
  int nt = 0;  // only add turning points for two directions
  // Don't compute dFXdFs_x for F and x slices nor for F and Fz
  // dFxdFz_x == FyFzz and dFxdFz_y == -FxFzz (Identity)
  for (int xyz = 0; xyz < 3 && nt < 2; ++xyz)
    if (!(jg >= 0 && xyz == 0) && !(jg == -1 && xyz == 1)) {
      dFxdG[xyz] = new DFxDGI(e->bf[jf], g, xyz);
      ++nt;
    }
  assert(nt == 2);

  for (int xyz = 0; xyz < 3; ++xyz)
    if (dFxdG[xyz] != 0) {
      vector<PTR<Object<PV3>>> roots = e->getRoots(jf, g, dFxdG[xyz]);
      fgd[xyz].resize(roots.size());
      for (int h = 0; h < roots.size(); ++h)
        fgd[xyz][h] = new Vert(jg < 0 ? FFzD : FSD, roots[h], jf, jg, xyz);
    }
  addTB();
}

void Encasement3D::FG::addTB() {
  for (int xyz = 0; xyz < 3; ++xyz)
    for (int i = 0; i < fgd[xyz].size(); ++i) top.addT(fgd[xyz][i]);

  PTR<Object<Poly3D>> f = e->bf[jf];
  PTR<Object<Poly3D>> g = this->g == 0 ? e->bf[jg] : this->g;
  VertType bType = this->g == 0 ? FGB : jg == -1 ? FFzB : FSB;

  vector<Vert*> bverts;

  for (int xyz = 0; xyz < 3; xyz++) {
    if (jf / 2 == xyz || (this->g == 0 && jg / 2 == xyz) ||
        (this->g != 0 && jg >= 0 && xyz == 0))
      continue;
    for (int lu = 0; lu < 2; lu++) {
      double d = lu == 0 ? e->box[xyz].lb() : e->box[xyz].ub();
      vector<PTR<Object<PV3>>> broots =
          getRoots(f, g, bType, jf, jg, xyz, d, e->box);
      for (int h = 0; h < broots.size(); ++h) {
        Vert* bvert = new Vert(bType, broots[h], jf, jg, xyz);
        fgb.push_back(bvert);
        bverts.push_back(bvert);
        bvert->c[2] = &top;
      }
    }
  }
  nboxbverts = bverts.size();

  top.addB(f, g, bType, jf, jg, bverts, fgb, planes, e->box);
}

vector<double> Encasement3D::getSlices(PTR<Object<Poly3D>> f,
                                       PTR<Object<Poly3D>> fz) {
  vector<PTR<Object<PV3>>> splits, joins;

  PTR<Object<Poly3D>> fy = new P3derI(f, 1);
  vector<PTR<Object<PV3>>> ffyfz = ::getRoots(f, fy, fz, box);
  for (auto it = ffyfz.begin(); it != ffyfz.end(); ++it) {
    int sfyyfzz_fyzfyz = DFxDGIV3(fy, fz, 0, *it);
    if (sfyyfzz_fyzfyz > 0) continue;
    int sfx = P3derIV3(f, 0, *it);
    int sfzz = P3derIV3(fz, 2, *it);
    if (sfx == sfzz)
      splits.push_back(*it);
    else
      joins.push_back(*it);
  }

  PTR<Object<Poly3D>> fzz = new P3derI(fz, 2);
  vector<PTR<Object<PV3>>> ffzfzz = ::getRoots(f, fz, fzz, box);
  for (auto it = ffzfzz.begin(); it != ffzfzz.end(); ++it) {
    int sfy = Poly3Dval(fy, *it);
    int sfzzz = P3derIV3(fzz, 2, *it);
    int sfxfyz_fyfxz = DFxDGIV3(f, fz, 2, *it);
    if (sfy * sfzzz * sfxfyz_fyfxz > 0)
      splits.push_back(*it);
    else
      joins.push_back(*it);
  }

  vector<pair<Object<PV3>*, bool>> all;

  for (Object<PV3>* split : splits)
    all.push_back(pair<Object<PV3>*, bool>(split, 0));

  for (Object<PV3>* join : joins)
    all.push_back(pair<Object<PV3>*, bool>(join, 1));

  sort(all.begin(), all.end(), XCompare());

  vector<double> slices;

  for (int i = 1; i < all.size(); ++i)
    if (!all[i - 1].second && all[i].second) {
      double d = (all[i - 1].first->getApprox().x.ub() +
                  all[i].first->getApprox().x.lb()) /
                 2;
      assert(V3IminusD(all[i - 1].first, 0, d) < 0);
      assert(V3IminusD(all[i].first, 0, d) > 0);
      slices.push_back(d);
    }

  // DEBUG
  if (false && slices.size() == 0 && ffyfz.size() == 2) {
    cout << "DEBUG in getSlices!!" << endl;
    double d =
        (ffyfz[0]->getApprox().x.ub() + ffyfz[1]->getApprox().x.lb()) / 2;
    slices.push_back(d);
  }

  return slices;

#ifdef BLEEN
  double d = box[0].lb();
  while (true) {
    PTR<Object<PV3>> sMin;
    for (auto it = splits.begin(); it != splits.end(); ++it)
      if (V3IminusD(*it, 0, d) > 0 &&
          (sMin == 0 || V3minusV3I(*it, sMin, 0) < 0))
        sMin = *it;

    if (sMin == 0) return slices;

    PTR<Object<PV3>> jMin;
    for (auto it = joins.begin(); it != joins.end(); ++it)
      if (V3minusV3I(*it, sMin, 0) > 0 &&
          (jMin == 0 || V3minusV3I(*it, jMin, 0) < 0))
        jMin = *it;

    if (jMin == 0) return slices;

    for (auto it = splits.begin(); it != splits.end(); ++it)
      if (*it != sMin && V3minusV3I(*it, sMin, 0) > 0 &&
          V3minusV3I(*it, jMin, 0) < 0)
        sMin = *it;

    d = (sMin->getApprox().x.ub() + jMin->getApprox().x.lb()) / 2;
    assert(V3IminusD(sMin, 0, d) < 0);
    assert(V3IminusD(jMin, 0, d) > 0);
    slices.push_back(d);
  }
#endif
}

Encasement3D::F::F(Encasement3D* e, int jf) : e(e), jf(jf) {
  e->f[jf] = this;

  if (jf < 6) return;

  Timer timer("    ");
  timer.Mark("start");

  fz = new P3derI(e->bf[jf], 2);

  ffzh.resize(e->bf.size());
  // intersection of f, fz, and x=c is a y turning point of f,x=c
  // intersection of f, fz, and y=c is a x turning point of f,y=c
  // so skip them
  for (int jh = 4; jh < e->bf.size(); jh++)
    if (jh != jf) {
      vector<PTR<Object<PV3>>> roots = e->getRoots(jh, e->bf[jf], fz);
      ffzh[jh].resize(roots.size());
      for (int h = 0; h < roots.size(); ++h)
        ffzh[jh][h] = new Vert(FFzH, roots[h], jf, -1, jh);
    }

  timer.Mark("ffzh getRoots()");

  ffz = new FG(e, jf, -1, fz);
  for (int jh = 0; jh < e->bf.size(); jh++) ffz->addIs(ffzh[jh]);
  for (int jh = 0; jh < 4; ++jh) {
    FG* hf = e->fg[e->ind(jh, jf)];
    // add y turning point for x=c and vice versa
    ffz->addIs(hf->fgd[1 - jh / 2]);
  }

  timer.Mark("ffz addIs()");

  vector<double> slicesD = e->getSlices(e->bf[jf], fz);
  slices.resize(slicesD.size());
  for (int i = 0; i < slicesD.size(); ++i)
    slices[i] = new Object<Poly3D>(slicePoly(slicesD[i]));

  timer.Mark("ffz getSlices()");

  fsh.resize(slices.size());
  // fsfz.resize(slices.size()); NO These are y turning points in fs encasement.
  for (int is = 0; is < slices.size(); ++is) {
    fsh[is].resize(e->bf.size());
    for (int jh = 2; jh < e->bf.size(); ++jh) {
      vector<PTR<Object<PV3>>> roots;
      if (jh < jf)
        roots = e->getRoots(slicesD[is], jh, jf);
      else if (jh > jf)
        roots = e->getRoots(slicesD[is], jf, jh);
      else
        continue;
      fsh[is][jh].resize(roots.size());
      for (int h = 0; h < roots.size(); ++h)
        fsh[is][jh][h] = new Vert(FSH, roots[h], jf, is, jh);
    }

    /* DEPRECATED
    vector<PTR<Object<PV3>>> roots = e->getRoots(slicesD[is], jf, fz);
    fsfz[is].resize(roots.size());
    for (int h = 0; h < roots.size(); ++h)
      fsfz[is][h] = new Vert(FSFz, roots[h], jf, is, -1);
    ffz->addIs(fsfz[is]);
    */
  }

  timer.Mark("fsh getRoots()");

  fs.resize(slices.size());
  for (int is = 0; is < slices.size(); ++is) {
    fs[is] = new FG(e, jf, is, slices[is]);
    for (int jh = 2; jh < e->bf.size(); ++jh) fs[is]->addIs(fsh[is][jh]);
    // fs[is]->addIs(fsfz[is]); NO These are y turning points.
    fs[is]->setVertCels();

    ffz->addIs(fs[is]->fgd[1]);  // adding y turning points to ffz
  }

  timer.Mark("fsh addIs()");

  ffz->setVertCels();
  timer.Mark("ffz setVertCels()");

  timer.Report();
}

int Encasement3D::maxDiffCoord(Vert* u, Vert* v, int xyzAvoid, bool ffz_only) {
  int uxyz, vxyz;

  if (ffz_only) {
    uxyz = u->type == FFzD ? u->jh : -1;
    vxyz = v->type == FFzD ? v->jh : -1;
  } else {
    uxyz = u->type % 3 == 1 ? u->jh : -1;
    vxyz = v->type % 3 == 1 ? v->jh : -1;
  }

  int xyzMax = -1;
  Vert *uMax = 0, *vMax = 0;
  for (int xyz = 0; xyz < 3; ++xyz) {
    if (xyzAvoid != -1 && (xyz == xyzAvoid || xyz == uxyz || xyz == vxyz))
      continue;
    if (V3minusV3I(u->p, v->p, xyz) > 0) {
      if (xyzMax == -1 ||
          V3mV3ImV3mV3I(uMax->p, vMax->p, xyzMax, u->p, v->p, xyz) < 0) {
        xyzMax = xyz;
        uMax = u;
        vMax = v;
      }
    } else {
      if (xyzMax == -1 ||
          V3mV3ImV3mV3I(uMax->p, vMax->p, xyzMax, v->p, u->p, xyz) < 0) {
        xyzMax = xyz;
        uMax = v;
        vMax = u;
      }
    }
  }
  assert(xyzMax != -1);
  return xyzMax;
}

double getMiddle(Encasement3D::Vert* u, Encasement3D::Vert* v, int xyz) {
  assert(V3minusV3I(u->p, v->p, xyz) < 0);
  double m =
      (u->p->getApprox(1.0)[xyz].ub() + v->p->getApprox(1.0)[xyz].lb()) / 2;
  if (V3IminusD(u->p, xyz, m) < 0 && V3IminusD(v->p, xyz, m) > 0) return m;
  m = (u->p->getApprox()[xyz].ub() + v->p->getApprox()[xyz].lb()) / 2;
  double uu = u->p->getApprox()[xyz].ub();
  double ul = u->p->getApprox()[xyz].lb();
  double vu = v->p->getApprox()[xyz].ub();
  double vl = v->p->getApprox()[xyz].lb();
  double uum = uu - m;
  double ulm = ul - m;
  double vum = vu - m;
  double vlm = uu - m;
  assert(V3IminusD(u->p, xyz, m) < 0 && V3IminusD(v->p, xyz, m) > 0);
  return m;
}

void Encasement3D::Cel::addT(Vert* v) {
  if (0 <= xyz && xyz < 3) {
    if (V3IminusD(v->p, xyz, spv.s) < 0)
      lb.c->addT(v);
    else
      rb.c->addT(v);
    return;
  }

  assert(xyz != 3);

  if (spv.v == 0) {
    spv.v = v;
    return;
  }

  xyz = maxDiffCoord(v, spv.v);
  Vert* vL = V3minusV3I(v->p, spv.v->p, xyz) < 0 ? v : spv.v;
  Vert* vR = vL == spv.v ? v : spv.v;
  spv.s = getMiddle(vL, vR, xyz);
  lb.c = new Cel();
  rb.c = new Cel();
  lb.c->spv.v = vL;
  rb.c->spv.v = vR;
}

vector<PTR<Object<PV3>>> Encasement3D::getRoots(PTR<Object<Poly3D>> f,
                                                PTR<Object<Poly3D>> g,
                                                VertType bType, int jf, int jg,
                                                int xyz, double m,
                                                const PV3<Parameter>& box) {
  vector<PTR<Object<Poly3D>>> fs;
  vector<int> xyzs = {xyz};
  vector<double> ds = {m};

  if (bType == FGB) {
    int jfg[2] = {jf, jg};
    PTR<Object<Poly3D>> fg[2] = {f, g};
    for (int i = 0; i < 2; ++i)
      if (jfg[i] < 6) {
        int xyz = jfg[i] / 2;
        int lu = jfg[i] % 2;
        double d = lu == 0 ? box[xyz].lb() : box[xyz].ub();
        xyzs.push_back(xyz);
        ds.push_back(d);
      } else
        fs.push_back(fg[i]);
  } else if (bType == FFzB) {
    fs.push_back(f);
    fs.push_back(g);
  } else if (bType == FSB) {
    fs.push_back(f);
    Poly3D<Parameter> s = g->getApprox(1.0);
    int xyz = 0 * s.degX() + 1 * s.degY() + 2 * s.degZ();
    assert(s[0].lb() == s[0].ub());
    double d = -s[0].lb();
    xyzs.push_back(xyz);
    ds.push_back(d);
  }

  vector<PTR<Object<PV3>>> roots = ::getRoots(fs, xyzs, ds, box);

  return roots;
}

PV3<Parameter> leftBox(PV3<Parameter> box, int xyz, double m) {
  box[xyz] =
      Parameter::constant(box[xyz].lb()).interval(Parameter::constant(m));
  return box;
}

PV3<Parameter> rightBox(PV3<Parameter> box, int xyz, double m) {
  box[xyz] =
      Parameter::constant(m).interval(Parameter::constant(box[xyz].ub()));
  return box;
}

PTR<Object<Poly3D>> separatingPlane(PTR<Object<Poly3D>> f,
                                    PTR<Object<Poly3D>> g,
                                    vector<Encasement3D::Vert*>& bverts,
                                    const PV3<Parameter>& box) {
  return 0;
}

int maxCoord(const PV3<Parameter>& box, int xyzAvoid = -1) {
  std::vector<int> ind{0, 1, 2};

  int ret =
      static_cast<int>(std::max_element(ind.begin(), ind.end(),
                                        [&](const int& i, const int& j) {
                                          if (i == xyzAvoid) {
                                            return true;
                                          }
                                          if (j == xyzAvoid) {
                                            return false;
                                          }
                                          return (box[i].ub() - box[i].lb()) <
                                                 (box[j].ub() - box[j].lb());
                                        }) -
                       ind.begin());

  return ret;
  //  int xyzMax = 0;
  //  double dMax = box[0].ub() - box[0].lb();
  //  for (int xyz = 1; xyz < 3; xyz++) {
  //    double d = box[xyz].ub() - box[xyz].lb();
  //    if (dMax < d) {
  //      dMax = d;
  //      xyzMax = xyz;
  //    }
  //  }
  //  return xyzMax;
}

void Encasement3D::Cel::addB(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g,
                             VertType bType, int jf, int jg,
                             vector<Vert*>& bverts,
                             vector<PTR<Vert>>& bvertsAll,
                             vector<PTR<Object<Poly3D>>>& planes,
                             const PV3<Parameter>& box) {
  if (0 <= xyz && xyz < 3) {
    vector<Vert*> bvertsL, bvertsR;
    for (int h = 0; h < bverts.size(); ++h)
      if (V3IminusD(bverts[h]->p, xyz, spv.s) < 0)
        bvertsL.push_back(bverts[h]);
      else
        bvertsR.push_back(bverts[h]);

    vector<PTR<Object<PV3>>> broots =
        getRoots(f, g, bType, jf, jg, xyz, spv.s, box);
    for (int h = 0; h < broots.size(); ++h) {
      Vert* bvert = new Vert(bType, broots[h], jf, jg, xyz);
      bvertsAll.push_back(bvert);
      bvertsL.push_back(bvert);
      bvertsR.push_back(bvert);
    }

    lb.c->addB(f, g, bType, jf, jg, bvertsL, bvertsAll, planes,
               leftBox(box, xyz, spv.s));
    rb.c->addB(f, g, bType, jf, jg, bvertsR, bvertsAll, planes,
               rightBox(box, xyz, spv.s));
    return;
  }

  assert(xyz != 3);

  if (bverts.size() <= 2) {
    setbverts(bverts);
    return;
  }

  Vert* v = spv.v;
  spv.v = 0;

  if (bverts.size() <= 4) {
    vector<Vert*> bverts4 = bverts;
    if (bverts4.size() == 3) {  // if an odd number of boundary vertices
      assert(v != 0);           // "internal" vert actually lies on boundary
      bverts4.push_back(v);
    }
    PTR<Object<Poly3D>> p = separatingPlane(f, g, bverts4, box);
    if (p != 0) {
      planes.push_back(p);
      xyz = 3;
      spv.p = p;

      vector<Vert*> bvertsL, bvertsR;
      for (int h = 0; h < bverts.size(); ++h)
        if (Poly3Dval(spv.p, bverts[h]->p) < 0)
          bvertsL.push_back(bverts[h]);
        else
          bvertsR.push_back(bverts[h]);

      assert(bvertsL.size() <= 2);
      assert(bvertsR.size() <= 2);

      lb.c = new Cel();
      rb.c = new Cel();

      if (v != 0) {
        if (Poly3Dval(spv.p, v->p) < 0)
          lb.c->spv.v = v;
        else
          rb.c->spv.v = v;
      }

      lb.c->setbverts(bvertsL);
      rb.c->setbverts(bvertsR);
      return;
    }
  }

  // If f is a boundary, avoid splitting on that boundaries dimension
  // Should never be the case that jg < 6
  if (jf < 6) {
    assert(jg >= 6);
  }

  int avoid = -1;
  if (jf < 6) avoid = jf / 2;
  if (bType == FSB) avoid = 0;
  xyz = maxCoord(box, avoid);

  spv.s = (box[xyz].lb() + box[xyz].ub()) / 2;
  assert(box[xyz].lb() < spv.s && spv.s < box[xyz].ub());
  vector<Vert*> bvertsL, bvertsR;
  for (int h = 0; h < bverts.size(); ++h)
    if (V3IminusD(bverts[h]->p, xyz, spv.s) < 0)
      bvertsL.push_back(bverts[h]);
    else
      bvertsR.push_back(bverts[h]);

  vector<PTR<Object<PV3>>> broots =
      getRoots(f, g, bType, jf, jg, xyz, spv.s, box);
  for (int h = 0; h < broots.size(); ++h) {
    Vert* bvert = new Vert(bType, broots[h], jf, jg, xyz);
    bvertsAll.push_back(bvert);
    bvertsL.push_back(bvert);
    bvertsR.push_back(bvert);
  }

  lb.c = new Cel();
  rb.c = new Cel();

  if (v != 0) {
    if (V3IminusD(v->p, xyz, spv.s) < 0)
      lb.c->spv.v = v;
    else
      rb.c->spv.v = v;
  }

  lb.c->addB(f, g, bType, jf, jg, bvertsL, bvertsAll, planes,
             leftBox(box, xyz, spv.s));
  rb.c->addB(f, g, bType, jf, jg, bvertsR, bvertsAll, planes,
             rightBox(box, xyz, spv.s));
}

void Encasement3D::Cel::addI(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g,
                             VertType bType, int jf, int jg, Vert* v,
                             vector<PTR<Vert>>& bvertsAll, int xyzAvoid,
                             const PV3<Parameter>& box, Object<Poly3D>* p,
                             int side) {
  if (0 <= xyz && xyz < 3) {
    if (V3IminusD(v->p, xyz, spv.s) < 0)
      lb.c->addI(f, g, bType, jf, jg, v, bvertsAll, xyzAvoid,
                 leftBox(box, xyz, spv.s), p, side);
    else
      rb.c->addI(f, g, bType, jf, jg, v, bvertsAll, xyzAvoid,
                 rightBox(box, xyz, spv.s), p, side);
    return;
  }

  if (xyz == 3) {
    if (Poly3Dval(spv.p, v->p) < 0)
      lb.c->addI(f, g, bType, jf, jg, v, bvertsAll, xyzAvoid, box, spv.p, -1);
    else
      rb.c->addI(f, g, bType, jf, jg, v, bvertsAll, xyzAvoid, box, spv.p, 1);
    return;
  }

  if (spv.v == 0) {
    spv.v = v;
    return;
  }

  xyz = maxDiffCoord(v, spv.v, xyzAvoid, bType == FFzB);
  Vert* vL = V3minusV3I(v->p, spv.v->p, xyz) < 0 ? v : spv.v;
  Vert* vR = vL == spv.v ? v : spv.v;
  spv.s = getMiddle(vL, vR, xyz);

  Vert* bverts[2] = {lb.b, rb.b};
  lb.c = new Cel();
  rb.c = new Cel();
  lb.c->spv.v = vL;
  rb.c->spv.v = vR;
  vector<Vert*> bvertsL, bvertsR;
  for (int h = 0; h < 2; ++h)
    if (bverts[h] != 0) {
      if (V3IminusD(bverts[h]->p, xyz, spv.s) < 0)
        bvertsL.push_back(bverts[h]);
      else
        bvertsR.push_back(bverts[h]);
    }

  vector<PTR<Object<PV3>>> broots =
      getRoots(f, g, bType, jf, jg, xyz, spv.s, box);
  for (int h = 0; h < broots.size(); ++h)
    if (p == 0 || Poly3Dval(p, broots[h]) == side) {
      Vert* bvert = new Vert(bType, broots[h], jf, jg, xyz);
      bvertsAll.push_back(bvert);
      bvertsL.push_back(bvert);
      bvertsR.push_back(bvert);
    }

  lb.c->setbverts(bvertsL);
  rb.c->setbverts(bvertsR);
}

Encasement3D::Edge* Encasement3D::Cel::getEdge(PTR<Object<PV3>> p,
                                               int xyzAvoid) {
  if (0 <= xyz && xyz < 3) {
    if (V3IminusD(p, xyz, spv.s) < 0)
      return lb.c->getEdge(p, xyzAvoid);
    else
      return rb.c->getEdge(p, xyzAvoid);
  }

  if (xyz == 3) {
    if (Poly3Dval(spv.p, p) < 0)
      return lb.c->getEdge(p, xyzAvoid);
    else
      return rb.c->getEdge(p, xyzAvoid);
  }

  if (spv.v == 0) return lb.b->edge();

  int xyzvv = 0;
  for (int i = 0; i < 2; ++i)
    if (xyzAvoid == xyzvv || (spv.v->type % 3 == 1 && spv.v->jh == xyzvv))
      ++xyzvv;
  assert(!(xyzAvoid == xyzvv || (spv.v->type % 3 == 1 && spv.v->jh == xyzvv)));
  int s1 = V3minusV3I(p, spv.v->p, xyzvv);
  int s2 = V3minusV3I(lb.b->p, spv.v->p, xyzvv);

  // debug
  if (rb.b != 0) assert(V3minusV3I(rb.b->p, spv.v->p, xyzvv) == -s2);

  if (s1 == 0) {
    // cout << "getEdge identity" << endl;
    return lb.b->edge();
  }
  if (s1 == s2)
    return lb.b->edge();
  else
    return rb.b->edge();
}

void Encasement3D::Cel::setVertCels(Encasement3D* e, Object<Poly3D>* f,
                                    Object<Poly3D>* g) {
  if (0 <= xyz && xyz <= 3) {
    lb.c->setVertCels(e, f, g);
    rb.c->setVertCels(e, f, g);
    return;
  }

  if (spv.v != 0) {
    // See https://github.com/joemasterjohn/acp/issues/2
    // The below logic doesn't work for FGD or FSD verts that are also
    // added to FFz encasements. First check that G (which is Fz) is
    // not in the list of returned fgh. If not put it in c[0].
    Object<Poly3D>* fgh[3];
    spv.v->getFGH(e, fgh);

    // This should always be true. Sanity check.
    assert(std::find(fgh, fgh + 3, f) != fgh + 3);

    // Look for g, if not there put this Cel at c[0]
    if (std::find(fgh, fgh + 3, g) == fgh + 3) {
      assert(spv.v->c[0] == nullptr);  // No one else should have set c[0]
      assert(spv.v->type % 3 ==
             1);  // must be a turning point inserted into FFz
      spv.v->c[0] = this;

    } else {
      // Find  the index that ISN'T f or g

      // alternatively:
      // int i = static_cast<int>(std::find(fgh, fgh+3, f) - fgh)
      // int j = static_cast<int>(std::find(fgh, fgh+3, g) - fgh)
      // spv.v->c[3 - i - j] = this;

      for (int i = 0; i < 3; i++) {
        if (fgh[i] != f && fgh[i] != g) {
          assert(spv.v->c[i] == 0);
          spv.v->c[i] = this;
        }
      }
    }
  }

  if (lb.b != 0) {
    if (lb.b->c[2] != 0)  // box boundary point
      lb.b = 0;           // remove it
    else
      lb.b->setBCel(this);
  }
  if (rb.b != 0) {
    if (rb.b->c[2] != 0)  // box boundary point
      rb.b = 0;           // remove it
    else
      rb.b->setBCel(this);
  }
  if (lb.b == 0 && rb.b != 0) {
    lb.b = rb.b;
    rb.b = 0;
  }
}

void Encasement3D::Vert::getFGH(Encasement3D* e, Object<Poly3D>* fgh[3]) {
  switch (type) {
    case FGH:
      fgh[0] = e->bf[jf];
      fgh[1] = e->bf[jg];
      fgh[2] = e->bf[jh];
      break;
    case FGD:
      fgh[0] = e->bf[jf];
      fgh[1] = e->bf[jg];
      fgh[2] = e->fg[e->ind(jf, jg)]->dFxdG[jh];
      break;
    case FGB:
      fgh[0] = e->bf[jf];
      fgh[1] = e->bf[jg];
      fgh[2] = 0;
      break;
    case FFzH:
      fgh[0] = e->bf[jf];
      fgh[1] = e->f[jf]->fz;
      fgh[2] = e->bf[jh];
      break;
    case FFzD:
      fgh[0] = e->bf[jf];
      fgh[1] = e->f[jf]->fz;
      fgh[2] = e->f[jf]->ffz->dFxdG[jh];
      break;
    case FFzB:
      fgh[0] = e->bf[jf];
      fgh[1] = e->f[jf]->fz;
      fgh[2] = 0;
      break;
    case FSH:
      fgh[0] = e->bf[jf];
      fgh[1] = e->f[jf]->slices[jg];
      fgh[2] = e->bf[jh];
      break;
    case FSD:
      fgh[0] = e->bf[jf];
      fgh[1] = e->f[jf]->slices[jg];
      fgh[2] = e->f[jf]->fs[jg]->dFxdG[jh];
      break;
    case FSB:
      fgh[0] = e->bf[jf];
      fgh[1] = e->f[jf]->slices[jg];
      fgh[2] = 0;
      break;
    case FSFz:
      fgh[0] = e->bf[jf];
      fgh[1] = e->f[jf]->slices[jg];
      fgh[2] = e->f[jf]->fz;
      break;
    case FXY:
      fgh[0] = e->bf[jf];
      fgh[1] = 0;
      fgh[2] = 0;
      break;
    default:
      assert(0);
  }
}

void Encasement3D::addIs() {
  for (int jf = 0; jf < bf.size(); ++jf) {
    for (int jg = jf + 1; jg < bf.size(); ++jg) {
      if (jg < 6 && jf / 2 == jg / 2) continue;
      for (int jh = jg + 1; jh < bf.size(); ++jh) {
        if (jh < 6 && jg / 2 == jh / 2) continue;
        vector<PTR<Vert>> verts = fgh[ind(jf, jg, jh)];
        fg[ind(jf, jg)]->addIs(verts);
        fg[ind(jf, jh)]->addIs(verts);
        fg[ind(jg, jh)]->addIs(verts);
      }
    }

    if (jf < 6) continue;
    F* ff = f[jf];
    for (int jh = 0; jh < bf.size(); jh++)
      if (jf < jh) {
        if (jh < 6 && jf / 2 == jh / 2) continue;
        fg[ind(jf, jh)]->addIs(ff->ffzh[jh]);
        for (int is = 0; is < ff->slices.size(); ++is)
          fg[ind(jf, jh)]->addIs(ff->fsh[is][jh]);
      } else if (jh < jf) {
        if (jf < 6 && jf / 2 == jh / 2) continue;
        fg[ind(jh, jf)]->addIs(ff->ffzh[jh]);
        for (int is = 0; is < ff->slices.size(); ++is)
          fg[ind(jh, jf)]->addIs(ff->fsh[is][jh]);
      }
  }
}

void Encasement3D::setVertCels() {
  for (int jf = 0; jf < bf.size(); ++jf)
    for (int jg = jf + 1; jg < bf.size(); ++jg) {
      if (jg < 6 && jf / 2 == jg / 2) continue;
      fg[ind(jf, jg)]->setVertCels();
    }
}

void Encasement3D::FG::findEdges() {
  vector<int> xyzs(fgb.size());
  for (int ib = 0; ib < fgb.size(); ++ib) xyzs[ib] = fgb[ib]->jh;

  for (Vert* b : fgb) MakeSet<Vert>(b);
  for (Vert* b : fgb)
    for (int ic = 0; ic < 2; ++ic)
      if (b->c[ic]->spv.v == 0) {
        Cel* c = b->c[ic];
        if (c->lb.b != 0 && c->rb.b != 0) Union<Vert>(c->lb.b, c->rb.b);
      }
  for (Vert* b : fgb) Find<Vert>(b);
  for (Vert* b : fgb)
    if (b->parent() == b) {
      b->edge() = new Edge(e, b);
      edges.push_back(b->edge());
      e->f[jf]->dedges.push_back(DEdge(b->edge(), 0, 0));
      if (g == 0) e->f[jg]->dedges.push_back(DEdge(b->edge(), 1, 0));
      b->size() = 0;
    }
  for (Vert* b : fgb)
    if (b->size() != 0) {
      b->edge() = b->parent()->edge();
      b->size() = 0;
    }

  for (int ib = 0; ib < fgb.size(); ++ib) fgb[ib]->jh = xyzs[ib];

  for (Vert* b : fgb)
    for (int ic = 0; ic < 2; ++ic)
      if (b->c[ic]->spv.v != 0) {
        int xyz = b->jh;
        Vert* v = b->c[ic]->spv.v;
        int svb = V3minusV3I(b->p, v->p, xyz);
        int sdfxdg = 0;
        sdfxdg = DFxDGIV3(getF(), getG(), xyz, b->p);
        /*
        if (this->g != 0 || jg >= 6)
          //sdfxdg = Poly3Dval(dFxdG[xyz], b->p);
          sdfxdg = DFxDGIV3(getF(), getG(), xyz, b->p);
        else {
          int df[] = { 0, 0, 0 };
          df[jf/2] = jf % 2 == 0 ? -1 : 1;
          int dg[] = { 0, 0, 0 };
          dg[jg/2] = jg % 2 == 0 ? -1 : 1;
          sdfxdg = (df[(xyz+1)%3] * dg[(xyz+2)%3] -
                    df[(xyz+2)%3] * dg[(xyz+1)%3]);
        }
        */
        assert(sdfxdg != 0);
        if (svb * sdfxdg > 0) {
          assert(b->edge()->v[0] == 0);
          b->edge()->v[0] = v;
        } else {
          assert(b->edge()->v[1] == 0);
          b->edge()->v[1] = v;
        }
      }
}

void Encasement3D::F::setEdgeNexts() {
  for (DEdge dedge : dedges) {
    dedge.e->setNexts(dedge.fg());
  }
}

void Encasement3D::initEdges() {
  for (int jf = 0; jf < bf.size(); ++jf) {
    for (int jg = jf + 1; jg < bf.size(); ++jg) {
      if (jg < 6 && jf / 2 == jg / 2) continue;
      fg[ind(jf, jg)]->findEdges();
    }

    if (jf < 6) continue;
    F* ff = f[jf];
    ff->ffz->findEdges();

    for (int is = 0; is < ff->fs.size(); ++is) ff->fs[is]->findEdges();
  }

  for (int jf = 0; jf < bf.size(); ++jf) f[jf]->setEdgeNexts();
}

Encasement3D::DEdge Encasement3D::DEdge::next() {
  DEdge n = e->next[fg()][fb()];
  return DEdge(n.e, n.fg(), n.fb(), np());
}

Encasement3D::Loop*& Encasement3D::DEdge::loop() { return e->loop[fg()][fb()]; }

void Encasement3D::findLoops() {
  for (int jf = 0; jf < bf.size(); ++jf) {
    f[jf]->findLoops();
  }
}

void Encasement3D::F::findLoops() {
  for (DEdge dedge : dedges) {
    for (int fb = 0; fb < 2; ++fb) {
      DEdge d = dedge.twin(fb);
    }
  }

  auto canLoop = [](const DEdge& dedge) {
    DEdge d = dedge;
    do {
      if (d.isNull()) return false;
    } while ((d = d.next()) != dedge);
    return true;
  };

  for (DEdge dedge : dedges)
    for (int fb = 0; fb < 2; ++fb) {
      DEdge d = dedge.twin(fb);

      // if (d.loop() != 0 || d.next().isNull()) continue;
      if (d.loop() != 0 || !canLoop(d)) continue;

      Loop* loop = new Loop(d);
      loops.push_back(loop);
      DEdge d0 = d;
      do {
        d.loop() = loop;
        d = d.next();
      } while (d != d0);
    }
}

void Encasement3D::F::findFaces() {
  int xind = jf >= 6 ? 0 : ((jf / 2) + 1) % 3;

  int nComps = 0;
  vector<Vert*> vxmins, vxmaxs;
  for (Loop* loop : loops)
    if (loop->comp == -1) {
      Vert *vxmin = 0, *vxmax = 0;
      loop->setComp(nComps++, vxmin, vxmax, xind);
      vxmins.push_back(vxmin);
      vxmaxs.push_back(vxmax);
    }

  vector<double> xs;
  double x = e->box[xind].lb() - 1;
  while (true) {
    Vert* vr = 0;
    for (int i = 0; i < vxmins.size(); ++i)
      if (V3IminusD(vxmins[i]->p, xind, x) > 0 &&
          (vr == 0 || V3minusV3I(vxmaxs[i]->p, vr->p, xind) < 0))
        vr = vxmaxs[i];
    if (vr == 0) break;
    Vert* vl = 0;
    for (int i = 0; i < vxmins.size(); ++i)
      if (V3IminusD(vxmins[i]->p, xind, x) > 0 &&
          V3minusV3I(vxmins[i]->p, vr->p, xind) < 0 &&
          (vl == 0 || V3minusV3I(vxmins[i]->p, vl->p, xind) > 0))
        vl = vxmins[i];
    x = getMiddle(vl, vr, xind);
    xs.push_back(x);
  }

  if (nComps > 1) {
    for (double d : xs) {
      if (jf >= 6) {
        PV2<Parameter> box2(e->box[1], e->box[2]);
        FXEncasement e2(e, jf, d, box2);
        for (FXEncasement::Edge* edgeFX : e2.edges) {
          if (edgeFX->v[0]->type == FXEncasement::FD && edgeFX->v[0]->jg == 0)
            continue;
          FXEncasement::Vert* v[2];
          FXEncasement::expand(edgeFX, v);
          Loop* loop[2];
          for (int iv = 0; iv < 2; ++iv)
            loop[iv] = e->getLoop(v[iv], e2.xObject, jf, iv);
          if (!(loop[0] == loop[1] || loop[0]->sameFace(loop[1])))
            loop[0]->joinFace(loop[1]);
        }
      } else {
        vector<PTR<Vert>> verts;
        for (int jg = 0; jg < e->bf.size(); ++jg) {
          if (jg / 2 == jf / 2) continue;
          if (jg / 2 == xind) continue;
          vector<PTR<Object<PV3>>> roots =
              e->getRoots(e->bf[jf], e->bf[jg], FGB, jf, jg, xind, d, e->box);
          for (PTR<Object<PV3>> root : roots)
            verts.push_back(new Vert(FGB, root, jf, jg, xind));
        }
        sort(verts.begin(), verts.end(), Compare((xind + 1) % 3));
        for (int h = 0; h < verts.size() - 1; ++h) {
          Loop* loop[2];
          for (int iv = 0; iv < 2; ++iv)
            loop[iv] = e->getLoop(verts[h + iv], iv, xind, jf % 2);
          if (!(loop[0] == loop[1] || loop[0]->sameFace(loop[1])))
            loop[0]->joinFace(loop[1]);
        }
      }
    }
  }

  for (PTR<Loop> loop : loops)
    if (loop->face == 0) {
      Face* face = new Face(e, jf, loop);
      faces.push_back(face);
      loop->setFace(face);
    }
}

// v is edge->v[iv] for an Edge of the FXEncasement
Encasement3D::DEdge Encasement3D::getDEdge(FXEncasement::Vert* v,
                                           PTR<Object<Scalar>> x, int jf,
                                           int iv) {
  PTR<Object<PV3>> p3 = new P2to3(v->p, x, 0);
  int fg = 0;
  Edge* edge = 0;

  if (v->type == FXEncasement::FD) {
    assert(v->jg == 1);  // should be z turning point
    edge = f[jf]->ffz->getEdge(p3);
  } else if (jf < v->jg)
    edge = this->fg[ind(jf, v->jg)]->getEdge(p3);
  else {
    edge = this->fg[ind(v->jg, jf)]->getEdge(p3);
    fg = 1;
  }

  int sv01x = V3minusV3I(edge->v[0]->p, edge->v[1]->p, 0);
  int fb = (sv01x > 0) ^ iv;

  return DEdge(edge, fg, fb);
}

// v is the intersection of x_{xyz}=d with a curve on a face of the box
// We consider pairs of consecutive vertices.
// iv=0 means lower, iv=1 means upper
// xyz is the prevailing left right direction
// lu = 0 or 1 is the lb() or ub() face
// return the Loop that is to that side of the curve (above for iv=0, below for
// iv=1)
Encasement3D::Loop* Encasement3D::getLoop(Vert* v, int iv, int xyz, int lu) {
  Edge* edge = 0;

  int fg = 0;

  if (v->jf < v->jg)
    edge = this->fg[ind(v->jf, v->jg)]->getEdge(v->p);
  else {
    edge = this->fg[ind(v->jg, v->jf)]->getEdge(v->p);
    fg = 1;
  }

  int sv01 = V3minusV3I(edge->v[0]->p, edge->v[1]->p, xyz);

  int fb = iv ^ !lu ^ (sv01 > 0);

  return DEdge(edge, fg, fb).loop();
}

void Encasement3D::F::findShells() {
  for (Face* face : faces) {
    for (int np = 0; np < 2; ++np)
      if (face->shell[np] == 0 && face->hasShell(np)) {
        Shell* shell = new Shell(e);
        e->shells.push_back(shell);
        face->setShell(np, shell);
      }
  }
}

void Encasement3D::Face::setShell(int np, Shell* shell) {
  this->shell[np] = shell;
  shell->faces.push_back(this);
  Loop* l = loop;
  do {
    DEdge d = l->dedge;
    do {
      DEdge d2 = d.sameShell(np);
      int np2 = d2.np();
      if (d2.loop()->face->shell[np2] == 0)
        d2.loop()->face->setShell(np2, shell);
      d = d.next();
    } while (d != l->dedge);
    l = l->next;
  } while (l != loop);
}

void Encasement3D::findCells() {
  int ncomps = 0;
  for (Shell* shell : shells)
    if (shell->comp == -1) shell->setComp(ncomps++);

  // Because of the order shells are created, the shell with comp=0
  // must touch the left boundary.

  PTR<Object<Scalar>> lbz =
      new Object<Scalar>(Parameter::constant(box[2].lb()));
  PTR<Object<Scalar>> ubz =
      new Object<Scalar>(Parameter::constant(box[2].ub()));

  for (int icomp = 1; icomp < ncomps; ++icomp) {
    Shell* shell = 0;
    for (Shell* s : shells)
      if (s->comp == icomp) {
        shell = s;
        break;
      }
    PTR<Object<PV3>> pof = getPointOnFace(shell->faces[0]);
    vector<PTR<Vert>> verts;
    for (int jf = 4; jf < bf.size(); ++jf) {
      PTR<Object<Poly>> fxy = new P3V3xy(bf[jf], pof);

      // DEBUG
      Poly3D<Parameter> fp = bf[jf]->getApprox();
      PV3<Parameter> pofp = pof->getApprox();
      Parameter val = fp.value(pofp);

      vector<PTR<Object<Scalar>>> roots;
      if (jf == 4)
        roots.push_back(lbz);
      else if (jf == 5)
        roots.push_back(ubz);
      else
        roots = acp::getRoots(fxy, lbz, ubz);
      for (PTR<Object<Scalar>> z : roots) {
        Vert* vxy = new Vert(FXY, new V3xyS(pof, z), jf, -1, -1);
        verts.push_back(vxy);
      }
    }
    sort(verts.begin(), verts.end(), Compare(2));
    vector<int> sfz(verts.size());
    vector<Face*> faces(verts.size());
    for (int h = 0; h < verts.size(); ++h) {
      sfz[h] = P3derIV3(bf[verts[h]->jf], 2, verts[h]->p);
      if (h == 0) {
        assert(verts[h]->jf == 4);
        faces[h] = f[4]->faces[0];
      } else if (h == verts.size() - 1) {
        assert(verts[h]->jf == 5);
        faces[h] = f[5]->faces[0];
      } else
        faces[h] = f[verts[h]->jf]->getFace(verts[h]->p);
    }
    for (int h = 1; h < verts.size(); ++h) {
      Shell* sl = faces[h - 1]->shell[sfz[h - 1] > 0];
      Shell* su = faces[h]->shell[sfz[h] < 0];
      if (!sl->sameCell(su)) sl->joinCell(su);
    }
  }

  for (Shell* shell : shells)
    if (shell->cell == 0) {
      Cell* cell = new Cell(this);
      cells.push_back(cell);
      shell->setCell(cell);
    }
}

PTR<Object<PV3>> Encasement3D::F::getPointOnFace(Face* face) {
  Vert *vxmin = 0, *vxmax = 0;
  Loop* loop = face->loop;
  do {
    DEdge dedge = loop->dedge;
    do {
      Vert* tail = dedge.tail();
      if (vxmin == 0)
        vxmin = vxmax = tail;
      else if (V3minusV3I(tail->p, vxmin->p, 0) < 0)
        vxmin = tail;
      else if (V3minusV3I(tail->p, vxmax->p, 0) > 0)
        vxmax = tail;
      dedge = dedge.next();
    } while (dedge != loop->dedge);
    loop = loop->next;
  } while (loop != face->loop);

  double x = getMiddle(vxmin, vxmax, 0);
  PV2<Parameter> box2(e->box[1], e->box[2]);
  FXEncasement e2(e, jf, x, box2);

  for (FXEncasement::Edge* edgeFX : e2.edges) {
    if (edgeFX->v[0]->type == FXEncasement::FD && edgeFX->v[0]->jg == 0)
      continue;
    FXEncasement::Vert* v[2];
    FXEncasement::expand(edgeFX, v);
    Loop* loop[2];
    for (int iv = 0; iv < 2; ++iv)
      loop[iv] = e->getLoop(v[iv], e2.xObject, jf, iv);
    assert(loop[0]->face == loop[1]->face);
    if (loop[0]->face != face) continue;
    return new P2to3(edgeFX->b->p, e2.xObject, 0);
  }

  assert(0);
  return 0;
}

Encasement3D::Face* Encasement3D::F::getFace(PTR<Object<PV3>> p) {
  PV3<Parameter> pp = p->getApprox(1.0);
  assert(pp.x.lb() == pp.x.ub());
  double x = pp.x.lb();
  PTR<Object<PV2>> p2 = new V3ItoV2(p, 0);
  PV2<Parameter> box2(e->box[1], e->box[2]);
  FXEncasement e2(e, jf, x, box2);
  FXEncasement::Edge* edgeFX = e2.top.getEdge(p2);
  FXEncasement::Vert* v[2];
  FXEncasement::expand(edgeFX, v);
  Loop* loop[2];
  for (int iv = 0; iv < 2; ++iv)
    loop[iv] = e->getLoop(v[iv], e2.xObject, jf, iv);
  assert(loop[0]->face == loop[1]->face);
  return loop[0]->face;
}

Encasement3D::Cell* Encasement3D::getCell(PTR<Object<PV3>> p) {
  for (int xyz = 0; xyz < 3; ++xyz)
    for (int lu = 0; lu < 2; ++lu) {
      double d = lu == 0 ? box[xyz].lb() : box[xyz].ub();
      int s = V3IminusD(p, xyz, d);
      if (s != 1 - 2 * lu) return 0;
    }

  PTR<Object<Scalar>> lbz = new Object<Scalar>(box[2].lbP());
  PTR<Object<Scalar>> ubz = new Object<Scalar>(box[2].ubP());
  int jfMin = -1;
  PTR<Object<Scalar>> zMin;
  for (int jf = 4; jf < bf.size(); ++jf) {
    PTR<Object<Poly>> fxy = new P3V3xy(bf[jf], p);
    vector<PTR<Object<Scalar>>> roots = acp::getRoots(fxy, lbz, ubz);
    for (PTR<Object<Scalar>> z : roots) {
      if (V3IminusS(p, 2, z) > 0) continue;
      if (jfMin == -1 || SminusS(z, zMin) < 0) {
        jfMin = jf;
        zMin = z;
      }
    }
  }

  PTR<Object<PV3>> q = new V3xyS(p, zMin);
  int sfz = P3derIV3(bf[jfMin], 2, q);
  Face* face = f[jfMin]->getFace(q);
  return face->shell[sfz < 0]->cell;
}
