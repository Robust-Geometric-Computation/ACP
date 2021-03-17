#include <vector>
#include "acp/core/object.h"
#include "acp/linmath/pv.h"
using namespace acp;
using namespace std;

typedef PTR<Object<PV3>> Vert;

map<Object<PV3>*, int> vindex;
vector<Vert> verts;
Vert indexVert(Vert v) {
  vindex[v] = verts.size();
  verts.push_back(v);
  return v;
}

class ZxN : public Object<PV2> {
  vector<Vert> abc;
  DeclareCalculate(PV2) {
    PV3<N> a = abc[0]->get<N>();
    PV3<N> b = abc[1]->get<N>();
    PV3<N> c = abc[2]->get<N>();
    PV3<N> n = (b - a).cross(c - a);
    return PV2<N>(-n.y, n.x);
  }

 public:
  ZxN(vector<Vert> abc) : abc(abc) { assert(abc.size() == 3); }
};

class ZxNVert : public Object<PV3> {
  vector<Vert> abc;
  int i;
  DeclareCalculate(PV3) {
    PV3<N> a = abc[0]->get<N>();
    PV3<N> b = abc[1]->get<N>();
    PV3<N> c = abc[2]->get<N>();
    PV3<N> n = (b - a).cross(c - a);
    PV3<N> zxn = PV3<N>(-n.y, n.x, N(0));
    PV3<N> p = abc[i]->get<N>();
    PV3<N> q = abc[(i + 1) % 3]->get<N>();
    N pd = p.dot(zxn);
    N qd = q.dot(zxn);
    return p * (qd / (qd - pd)) - q * (pd / (qd - pd));
  }

 public:
  ZxNVert(vector<Vert> abc, int i) : abc(abc), i(i) {
    assert(abc.size() == 3 && 0 <= i && i < 3);
  }
};

class Dot23 : public Primitive {
  Object<PV2>* p2;
  Object<PV3>* p3;
  DeclareSign {
    PV2<N> p2p = p2->get<N>();
    PV3<N> p3p = p3->get<N>();
    return p2p.x * p3p.x + p2p.y * p3p.y;
  }

 public:
  Dot23(Object<PV2>* p2, Object<PV3>* p3) : p2(p2), p3(p3) {}
};

class DiffZ : public Primitive {
  Object<PV3>*a, *b;
  DeclareSign {
    if (a == b) return N(0);
    return a->get<N>().z - b->get<N>().z;
  }

 public:
  DiffZ(Object<PV3>* a, Object<PV3>* b) : a(a), b(b) {}
};

class SplitABatCz : public Object<PV3> {
  Object<PV3>*a, *b, *c;
  DeclareCalculate(PV3) {
    PV3<N> ap = a->get<N>();
    PV3<N> bp = b->get<N>();
    N cz = c->get<N>().z;
    N s = (cz - ap.z) / (bp.z - ap.z);
    return PV3<N>(ap.x + (bp.x - ap.x) * s, ap.y + (bp.y - ap.y) * s, cz);
  }

 public:
  SplitABatCz(Object<PV3>* a, Object<PV3>* b, Object<PV3>* c)
      : a(a), b(b), c(c) {}
};

vector<vector<Vert>> splitZxN(vector<Vert>& tri) {
  PTR<Object<PV2>> zxn = new ZxN(tri);
  vector<vector<Vert>> polys;
  int firstPlus = -1;
  int firstMinus = -1;
  int signs[6];
  bool noNeg = true, noPos = true;
  for (int i = 0; i < 3; ++i) {
    signs[i] = Dot23(zxn, tri[i]);
    signs[i + 3] = signs[i];
    if (signs[i] > 0) noNeg = false;
    if (signs[i] < 0) noPos = false;
  }
  if (noNeg || noPos) {
    polys.push_back(tri);
    return polys;
  }
  Vert p2m, m2p;
  for (int s = -1; s <= 1; s += 2) {
    vector<Vert> poly;
    for (int i = 1; i < 6; ++i)
      if (poly.size() == 0 && signs[i - 1] * s <= 0 && signs[i] * s > 0) {
        if (signs[i - 1] == 0)
          poly.push_back(tri[(i - 1) % 3]);
        else if (s == -1) {
          poly.push_back(indexVert(new ZxNVert(tri, (i - 1) % 3)));
          m2p = poly.back();
        } else
          poly.push_back(p2m);
        poly.push_back(tri[i % 3]);
      } else if (poly.size() > 0 && signs[i] * s > 0) {
        assert(signs[i - 1] * s > 0);
        poly.push_back(tri[i % 3]);
      } else if (poly.size() > 0 && signs[i] * s <= 0) {
        assert(signs[i - 1] * s > 0);
        if (signs[i] == 0)
          poly.push_back(tri[i % 3]);
        else if (s == -1) {
          poly.push_back(indexVert(new ZxNVert(tri, (i - 1) % 3)));
          p2m = poly.back();
        } else
          poly.push_back(m2p);
        break;
      }
    polys.push_back(poly);
  }

  cerr << "splitZxN";
  for (int i = 0; i < 3; i++) cerr << " " << vindex[tri[i]];
  cerr << " into " << polys.size() << endl;
  for (vector<Vert>& poly : polys) {
    for (Vert& v : poly) cerr << vindex[v] << " ";
    cerr << endl;
  }
  return polys;
}

void shiftTo0(vector<Vert>& poly, int i) {
  assert(0 <= i && i < poly.size());
  vector<Vert> copy = poly;
  for (int j = i; j < poly.size(); ++j) poly[j - i] = copy[j];
  for (int j = 0; j < i; ++j) poly[j + poly.size() - i] = copy[j];
}

void duplicate(vector<Vert>& poly, int i) {
  assert(0 <= i && i < poly.size());
  poly.insert(poly.begin() + i, poly[i]);
}

vector<Vert> splitLeft(vector<Vert>& poly) {
  Vert v = indexVert(new SplitABatCz(poly[1], poly[2], poly.back()));
  vector<Vert> trap(4);
  trap[0] = poly[0];
  trap[1] = poly[1];
  trap[2] = v;
  trap[3] = poly.back();
  poly.pop_back();
  poly[0] = trap[3];
  poly[1] = trap[2];
  return trap;
}

vector<Vert> splitRight(vector<Vert>& poly) {
  Vert v = indexVert(new SplitABatCz(poly.back(), poly[0], poly[2]));
  vector<Vert> trap(4);
  trap[0] = poly[0];
  trap[1] = poly[1];
  trap[2] = poly[2];
  trap[3] = v;
  poly.erase(poly.begin() + 2);
  poly[0] = trap[3];
  poly[1] = trap[2];
  return trap;
}

vector<Vert> splitBoth(vector<Vert>& poly) {
  vector<Vert> trap(4);
  trap[0] = poly[0];
  trap[1] = poly[1];
  trap[2] = poly[2];
  trap[3] = poly.back();
  poly.pop_back();
  poly.erase(poly.begin() + 2);
  poly[0] = trap[3];
  poly[1] = trap[2];
  return trap;
}

void print(const char* s, vector<Vert>& poly) {
  cerr << s;
  for (int j = 0; j < poly.size(); j++) cerr << " " << vindex[poly[j]];
  cerr << endl;
}

void trapezoidalize(vector<Vert>& poly, vector<vector<Vert>>& traps) {
  print("trapezoidalize", poly);
  int n = poly.size();
  int i = 0;
  while (DiffZ(poly[i], poly[(i + 1) % n]) <= 0) i = (i + 1) % n;
  while (DiffZ(poly[i], poly[(i + 1) % n]) > 0) i = (i + 1) % n;
  shiftTo0(poly, i);
  print("shift", poly);
  if (DiffZ(poly[0], poly[1]) < 0) duplicate(poly, 0);
  print("dup0", poly);
  n = poly.size();
  i = 1;
  while (DiffZ(poly[i], poly[(i + 1) % n]) < 0) i = (i + 1) % n;
  if (DiffZ(poly[i], poly[(i + 1) % n]) > 0) duplicate(poly, i);
  print("dupi", poly);
  while (poly.size() > 4) {
    assert(DiffZ(poly.back(), poly[i]) < 0 || DiffZ(poly[2], poly[i]) < 0);
    if (DiffZ(poly.back(), poly[2]) < 0)
      traps.push_back(splitLeft(poly));
    else if (DiffZ(poly.back(), poly[2]) > 0)
      traps.push_back(splitRight(poly));
    else
      traps.push_back(splitBoth(poly));
    print("trap", traps.back());
    print("poly", poly);
  }
  traps.push_back(poly);
}

void triangulate(int i0, vector<Vert>& poly, vector<vector<Vert>>& tris) {
  for (int i = 2; i < poly.size(); ++i) {
    vector<Vert> tri(3);
    tri[0] = poly[i0];
    tri[1] = poly[(i0 + i - 1) % poly.size()];
    tri[2] = poly[(i0 + i) % poly.size()];
    if (tri[0] == tri[1] || tri[0] == tri[2] || tri[1] == tri[2]) continue;
    tris.push_back(tri);
  }
}

class A2mB2 : public Primitive {
  Vert a, b;
  DeclareSign {
    return a->get<N>().dot(a->get<N>()) - b->get<N>().dot(b->get<N>());
  }

 public:
  A2mB2(Vert a, Vert b) : a(a), b(b) {}
};

class Rot : public Object<PV3> {
  PTR<Object<PV2>> u;
  Vert p;
  DeclareCalculate(PV3) {
    PV2<N> up = u->get<N>();
    PV3<N> pp = p->get<N>();
    return PV3<N>(up.x * pp.x - up.y * pp.y, up.x * pp.y + up.y * pp.x, pp.z);
  }

 public:
  Rot(PTR<Object<PV2>> u, Vert p) : u(u), p(p) {}
};

class MidTan : public Object<PV3> {
  PTR<Object<PV2>> u;
  Vert p;
  DeclareCalculate(PV3) {
    PV2<N> up = u->get<N>();
    N t = (1 - up.x) / up.y;
    PV3<N> pp = p->get<N>();
    return PV3<N>(pp.x - t * pp.y, pp.y + t * pp.x, pp.z);
  }

 public:
  MidTan(PTR<Object<PV2>> u, Vert p) : u(u), p(p) {}
};

class SignXYZ : public Primitive {
  Vert v;
  int i;
  DeclareSign { return v->get<N>()[i]; }

 public:
  SignXYZ(Vert v, int i) : v(v), i(i) {}
};

map<Object<PV3>*, Object<PV3>*> saveMidTan;
Vert getMidTan(PTR<Object<PV2>> u, Vert v) {
  if (SignXYZ(v, 0) == 0 && SignXYZ(v, 1) == 0) return v;
  if (saveMidTan.find(v) == saveMidTan.end())
    return indexVert(saveMidTan[v] = new MidTan(u, v));
  return saveMidTan[v];
}

map<Object<PV3>*, Object<PV3>*> saveRot;
Vert getRot(PTR<Object<PV2>> u, Vert v) {
  if (SignXYZ(v, 0) == 0 && SignXYZ(v, 1) == 0) return v;
  if (saveRot.find(v) == saveRot.end())
    return indexVert(saveRot[v] = new Rot(u, v));
  return saveRot[v];
}

// Eventually need to avoid duplicate Rot and MidTan!
vector<Vert> pentagon(PTR<Object<PV2>> u, Vert a, Vert b) {
  vector<Vert> pent;
  pent.push_back(a);
  pent.push_back(b);
  pent.push_back(getMidTan(u, b));
  pent.push_back(getRot(u, b));
  pent.push_back(getRot(u, a));
  return pent;
}

class ABxCD : public Primitive {
  Vert a, b, c, d;
  DeclareSign {
    return (b->get<N>() - a->get<N>()).cross(d->get<N>() - c->get<N>()).z;
  }

 public:
  ABxCD(Vert a, Vert b, Vert c, Vert d) : a(a), b(b), c(c), d(d) {}
};

void outerTriangles(int i, vector<Vert>& bot, vector<Vert>& top,
                    vector<vector<Vert>>& triangles) {
  int j = (i + 1) % 5;

  if (bot[i] == bot[j]) {
    vector<Vert> tri;
    tri.push_back(bot[i]);
    tri.push_back(top[j]);
    tri.push_back(top[i]);
    triangles.push_back(tri);
    return;
  }

  if (top[i] == top[j]) {
    vector<Vert> tri;
    tri.push_back(bot[i]);
    tri.push_back(bot[j]);
    tri.push_back(top[i]);
    triangles.push_back(tri);
    return;
  }

  if (i == 0 || i == 3 || ABxCD(bot[i], bot[j], top[i], top[j]) > 0) {
    vector<Vert> tri;
    tri.push_back(bot[i]);
    tri.push_back(bot[j]);
    tri.push_back(top[i]);
    triangles.push_back(tri);
    tri.clear();
    tri.push_back(bot[j]);
    tri.push_back(top[j]);
    tri.push_back(top[i]);
    triangles.push_back(tri);
    return;
  }

  {
    vector<Vert> tri;
    tri.push_back(bot[i]);
    tri.push_back(bot[j]);
    tri.push_back(top[j]);
    triangles.push_back(tri);
    tri.clear();
    tri.push_back(bot[i]);
    tri.push_back(top[j]);
    tri.push_back(top[i]);
    triangles.push_back(tri);
    return;
  }
}

void sweep(PTR<Object<PV2>> u, vector<Vert>& trap, vector<vector<Vert>>& tris) {
  cerr << "sweep trap";
  for (int i = 0; i < 4; i++) cerr << " " << vindex[trap[i]];
  cerr << endl;
  int sgnBot = A2mB2(trap[0], trap[1]);
  int sgnTop = A2mB2(trap[3], trap[2]);
  assert(sgnBot * sgnTop != -1);
  bool pos = sgnBot < 0 || sgnTop < 0;
  vector<Vert> bot =
      pos ? pentagon(u, trap[0], trap[1]) : pentagon(u, trap[1], trap[0]);
  vector<Vert> top =
      pos ? pentagon(u, trap[3], trap[2]) : pentagon(u, trap[2], trap[3]);
  for (int i = 0; i < 5; i++) outerTriangles(i, bot, top, tris);
  reverse(bot.begin(), bot.end());
  triangulate(pos ? 0 : 4, bot, tris);
  triangulate(pos ? 4 : 0, top, tris);
}

vector<vector<Vert>> sweep(PTR<Object<PV2>> u, vector<Vert>& triangle) {
  vector<vector<Vert>> triangles;
  if (DiffZ(triangle[0], triangle[1]) == 0 &&
      DiffZ(triangle[0], triangle[2]) == 0)
    return triangles;
  // split by (z x n) * p = 0
  vector<vector<Vert>> split = splitZxN(triangle);
  vector<vector<Vert>> traps;
  for (vector<Vert>& poly : split)
    // split by z = v.z for each Vert
    trapezoidalize(poly, traps);
  for (vector<Vert>& trap : traps)
    // sweep each trapezoid
    sweep(u, trap, triangles);
  return triangles;
}

vector<Vert> readTriangle(istream& in) {
  vector<Vert> tri;
  for (int i = 0; i < 3; i++) {
    PV3<Parameter> p;
    for (int j = 0; j < 3; j++) {
      double x;
      in >> x;
      p[j] = Parameter::constant(x);
    }
    tri.push_back(indexVert(new Object<PV3>(p)));
  }
  return tri;
}

void writeVTK(ostream& out, vector<vector<Vert>>& triangles) {
  out << "# vtk DataFile Version 3.0" << endl
      << "vtk output" << endl
      << "ASCII" << endl
      << "DATASET POLYDATA" << endl
      << "POINTS " << verts.size() << " double" << endl;

  for (Vert& p : verts) {
    PV3<Parameter> pp = p->getApprox();
    out << pp.x.mid() << " " << pp.y.mid() << " " << pp.z.mid() << endl;
  }

  out << endl
      << "POLYGONS " << triangles.size() << " " << triangles.size() * 4 << endl;
  for (vector<Vert>& tri : triangles) {
    out << 3;
    for (Vert& v : tri) {
      assert(vindex.find(v) != vindex.end());
      out << " " << vindex[v];
    }
    out << endl;
  }
}

class U2D : public Object<PV2> {
  PTR<Object<Scalar>> s;
  DeclareCalculate(PV2) {
    N t = s->get<N>().x;
    return PV2<N>((1 - t * t) / (1 + t * t), 2 * t / (1 + t * t));
  }

 public:
  U2D(PTR<Object<Scalar>> s) : s(s) {}
};

int main(int argc, char* argv[]) {
  acp::enable();
  PTR<Object<PV2>> u = new U2D(new Object<Scalar>(Parameter::constant(0.25)));
  vector<Vert> triangle = readTriangle(cin);
  vector<vector<Vert>> triangles = sweep(u, triangle);
  writeVTK(cout, triangles);
  acp::disable();
}
