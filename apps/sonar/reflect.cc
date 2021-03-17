#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

#include "acp/core/object.h"
#include "acp/linmath/pv.h"
using namespace acp;

class TriangleNormal : public Object<PV3> {
  DeclareCalculate(PV3) {
    PV3<N> p0 = p[0]->get<N>();
    return (p[1]->get<N>() - p0).cross(p[2]->get<N>() - p0).unit();
  }

 public:
  PTR<Object<PV3>> p[3];
  TriangleNormal(PTR<Object<PV3>> a, PTR<Object<PV3>> b, PTR<Object<PV3>> c) {
    p[0] = a;
    p[1] = b;
    p[2] = c;
  }
};

class DoubleReflect : public Object<PV3> {
 public:
  PTR<TriangleNormal> t[2];
  DoubleReflect(PTR<TriangleNormal> t0, PTR<TriangleNormal> t1) {
    t[0] = t0;
    t[1] = t1;
  }

  DeclareCalculate(PV3) {
    // reflection of p is p - 2 (p-a)*n n
    // reflection of o is o - 2 (o-a)*n n = 2 a*n n
    PV3<N> n0 = t[0]->get<N>();
    PV3<N> a0 = t[0]->p[0]->get<N>();
    PV3<N> n1 = t[1]->get<N>();
    PV3<N> a1 = t[1]->p[0]->get<N>();

    PV3<N> p = n0 * (a0.dot(n0) * 2);
    return p - n1 * ((p - a1).dot(n1) * 2);
  }
};

class Dot : public Primitive {
  Object<PV3>*u, *v;
  DeclareSign { return u->get<N>().dot(v->get<N>()); }

 public:
  Dot(Object<PV3>* u, Object<PV3>* v) : u(u), v(v) {}
};

class TripleProduct : public Primitive {
  Object<PV3>*u, *v, *w;
  DeclareSign { return u->get<N>().cross(v->get<N>()).dot(w->get<N>()); }

 public:
  TripleProduct(Object<PV3>* u, Object<PV3>* v, Object<PV3>* w)
      : u(u), v(v), w(w) {}
};

class Orient3D : public Primitive {
  Object<PV3>*a, *b, *c, *d;
  DeclareSign {
    PV3<N> aN = a->get<N>();
    return (b->get<N>() - aN).cross(c->get<N>() - aN).dot(d->get<N>() - aN);
  }

 public:
  Orient3D(Object<PV3>* a, Object<PV3>* b, Object<PV3>* c, Object<PV3>* d)
      : a(a), b(b), c(c), d(d) {}
};

class Triangle : public RefCnt {
 public:
  PTR<TriangleNormal> n;

  Triangle(PTR<Object<PV3>> a, PTR<Object<PV3>> b, PTR<Object<PV3>> c)
      : n(new TriangleNormal(a, b, c)) {
    assert((Object<PV3>*)b != 0);
  }

  bool pierces(Object<PV3>* v) {
    return (Dot(v, n) < 0 && TripleProduct(n->p[0], n->p[1], v) < 0 &&
            TripleProduct(n->p[1], n->p[2], v) < 0 &&
            TripleProduct(n->p[2], n->p[0], v) < 0);
  }

  class PiercePoint : public Object<PV3> {
    PTR<TriangleNormal> n;
    PTR<Object<PV3>> v;
    DeclareCalculate(PV3) {
      // (t v - a) * n = 0
      // t v*n - a*n = 0
      // t = a*n / v*n
      PV3<N> nN = n->get<N>();
      PV3<N> vN = v->get<N>();
      return vN * (n->p[0]->get<N>().dot(nN) / vN.dot(nN));
    }

   public:
    PiercePoint(PTR<TriangleNormal> n, PTR<Object<PV3>> v) : n(n), v(v) {}
  };

  PTR<Object<PV3>> piercePoint(PTR<Object<PV3>> v) {
    return new PiercePoint(n, v);
  }

  class Side : public Primitive {
    Object<PV3>* p;
    TriangleNormal* n;
    DeclareSign { return (p->get<N>() - n->p[0]->get<N>()).dot(n->get<N>()); }

   public:
    Side(Object<PV3>* p, TriangleNormal* n) : n(n), p(p) {}
  };

  bool intersectsSegment(Object<PV3>* a, Object<PV3>* b) {
    int sa = Side(a, n);
    int sb = Side(b, n);
    if ((sa > 0) == (sb > 0)) return false;

    int s01 = Orient3D(a, b, n->p[0], n->p[1]);
    int s12 = Orient3D(a, b, n->p[1], n->p[2]);

    if (s01 != s12) return false;

    int s20 = Orient3D(a, b, n->p[2], n->p[0]);
    return s12 == s20;
  }
};

class TrianglePair {
 public:
  PTR<Triangle> t[2];
  PTR<Object<PV3>> p[2];  // pierce points
  TrianglePair(PTR<Triangle> t0, PTR<Triangle> t1, PTR<Object<PV3>> p0,
               PTR<Object<PV3>> p1) {
    t[0] = t0;
    t[1] = t1;
    p[0] = p0;
    p[1] = p1;
  }
};

#define CHECK_PIERCES
#define CHECK_BLOCKED

vector<TrianglePair> reflects(vector<PTR<Triangle>> ts) {
  vector<TrianglePair> pairs;

  PTR<Object<PV3>> o =
      new Object<PV3>(PV3<Parameter>(Parameter(0), Parameter(0), Parameter(0)));

  for (int i = 0; i < ts.size(); ++i) {
    PTR<Triangle> ti = ts[i];
    for (int j = i + 1; j < ts.size(); ++j) {
      PTR<Triangle> tj = ts[j];

      PTR<Object<PV3>> vij = new DoubleReflect(ti->n, tj->n);
#ifdef CHECK_PIERCES
      if (!tj->pierces(vij)) continue;
#endif

      PTR<Object<PV3>> vji = new DoubleReflect(tj->n, ti->n);
#ifdef CHECK_PIERCES
      if (!ti->pierces(vji)) continue;
#endif

      PTR<Object<PV3>> pi = ti->piercePoint(vji);
      PTR<Object<PV3>> pj = tj->piercePoint(vij);

#ifdef CHECK_BLOCKED
      bool blocked = false;
      for (int k = 0; k < ts.size(); ++k)
        if (k != i && k != j && ts[k]->intersectsSegment(pi, pj)) {
          blocked = true;
          break;
        }
      if (blocked) continue;

      for (int k = 0; k < ts.size(); ++k)
        if (k != i && k != j && ts[k]->intersectsSegment(o, pi)) {
          blocked = true;
          break;
        }
      if (blocked) continue;

      for (int k = 0; k < ts.size(); ++k)
        if (k != i && k != j && ts[k]->intersectsSegment(o, pj)) {
          blocked = true;
          break;
        }
      if (blocked) continue;
#endif

      pairs.push_back(TrianglePair(ti, tj, pi, pj));
    }
  }

  return pairs;
}

void rotate(int ixyz, double a, double xyz[3]) {
  disable();
  double c = cos(a);
  double s = sin(a);
  int j = (ixyz + 1) % 3;
  int k = (ixyz + 2) % 3;
  double x = c * xyz[j] - s * xyz[k];
  double y = s * xyz[j] + c * xyz[k];
  xyz[j] = x;
  xyz[k] = y;
}

vector<PTR<Triangle>> readTriangles(char* vertexFile, char* triangleFile) {
  vector<PTR<Triangle>> ts;

  ifstream in;
  in.open(vertexFile);

  vector<PTR<Object<PV3>>> vs;
  double w, xyz[3];
  while (in >> w >> xyz[0] >> xyz[1] >> xyz[2]) {
    rotate(2, 0.01, xyz);
    rotate(0, 0.1, xyz);
    rotate(1, 1, xyz);
    vs.push_back(new Object<PV3>(PV3<Parameter>(Parameter::input(xyz[0]),
                                                Parameter::input(xyz[1]),
                                                Parameter::input(xyz[2]))));
    assert(w == vs.size());
  }

  cout << vs.size() << " vertices" << endl;

  in.close();

  in.open(triangleFile);

  int i, j, k;
  while (in >> i >> j >> k)
    ts.push_back(new Triangle(vs[i - 1], vs[j - 1], vs[k - 1]));

  cout << ts.size() << " triangles" << endl;

  return ts;
}

class AngleLessThan : public Primitive {
  double width;
  PTR<Object<PV3>> u, v;
  DeclareSign {
    PV3<N> uN = u->get<N>();
    PV3<N> vN = v->get<N>();
    PV3<N> uxv = uN.cross(vN);
    N uDu = uN.dot(uN);
    N vDv = vN.dot(vN);
    return uxv.dot(uxv) / uDu / vDv / width / width - 1;
  }

 public:
  AngleLessThan(PTR<Object<PV3>> u, PTR<Object<PV3>> v, double width)
      : u(u), v(v), width(width) {}
};

void display(const char* file, vector<PTR<Triangle>>& ts,
             vector<TrianglePair>& pairs) {
  disable();

  double maxAngle = 14.0 * M_PI / 180;
  int nx = 1024;
  double reflectWidth = maxAngle / 128;

  vector<double> xs(nx);
  for (int i = 0; i < nx; ++i) xs[i] = tan(-maxAngle + i * 2 * maxAngle / nx);

  enable();

  cout << pairs.size() << " pairs" << endl;

  ofstream out;
  out.open(file);
  out << "P3" << endl;
  out << nx << " " << nx << endl;
  out << 255 << endl;

  PTR<Object<PV3>> o =
      new Object<PV3>(PV3<Parameter>(Parameter(0), Parameter(0), Parameter(0)));

  vector<double> minX(ts.size()), maxX(ts.size()), minZ(ts.size()),
      maxZ(ts.size());
  for (int it = 0; it < ts.size(); ++it)
    for (int ip = 0; ip < 3; ++ip) {
      if (ip == 0 || minX[it] > (ts[it]->n->p[ip]->getApprox(1.0).x /
                                 ts[it]->n->p[ip]->getApprox(1.0).y)
                                    .lb())
        minX[it] = (ts[it]->n->p[ip]->getApprox(1.0).x /
                    ts[it]->n->p[ip]->getApprox(1.0).y)
                       .lb();
      if (ip == 0 || maxX[it] < (ts[it]->n->p[ip]->getApprox(1.0).x /
                                 ts[it]->n->p[ip]->getApprox(1.0).y)
                                    .ub())
        maxX[it] = (ts[it]->n->p[ip]->getApprox(1.0).x /
                    ts[it]->n->p[ip]->getApprox(1.0).y)
                       .ub();
      if (ip == 0 || minZ[it] > (ts[it]->n->p[ip]->getApprox(1.0).z /
                                 ts[it]->n->p[ip]->getApprox(1.0).y)
                                    .lb())
        minZ[it] = (ts[it]->n->p[ip]->getApprox(1.0).z /
                    ts[it]->n->p[ip]->getApprox(1.0).y)
                       .lb();
      if (ip == 0 || maxZ[it] < (ts[it]->n->p[ip]->getApprox(1.0).z /
                                 ts[it]->n->p[ip]->getApprox(1.0).y)
                                    .ub())
        maxZ[it] = (ts[it]->n->p[ip]->getApprox(1.0).z /
                    ts[it]->n->p[ip]->getApprox(1.0).y)
                       .ub();
    }

  for (int j = nx - 1; j >= 0; --j) {
    for (int i = 0; i < nx; ++i) {
      PTR<Object<PV3>> v = new Object<PV3>(
          PV3<Parameter>(Parameter::constant(xs[i]), Parameter(1),
                         Parameter::constant(xs[j])));

      bool reflectSpot = false;
      for (TrianglePair pair : pairs)
        if (AngleLessThan(v, pair.p[0], reflectWidth) < 0 ||
            AngleLessThan(v, pair.p[1], reflectWidth) < 0) {
          reflectSpot = true;
          break;
        }
      if (reflectSpot) {
        out << "255 0 0" << endl;
        continue;
      }

      PTR<Object<PV3>> p;
      PTR<Triangle> nearT;
      for (int it = 0; it < ts.size(); ++it) {
        if (xs[i] < minX[it] || xs[i] > maxX[it] || xs[j] < minZ[it] ||
            xs[j] > maxZ[it])
          continue;
        PTR<Triangle> t = ts[it];
        if (t->pierces(v) && (p == 0 || t->intersectsSegment(o, p))) {
          p = t->piercePoint(v);
          nearT = t;
        }
      }

      if (nearT == 0)
        out << "0 0 0" << endl;
      else {
        PV3<Parameter> n = nearT->n->getApprox();
        if (false)
          out << (int)(255 * ((n.x.mid() + 1) / 2)) << " "
              << (int)(255 * ((n.y.mid() + 1) / 2)) << " "
              << (int)(255 * ((n.z.mid() + 1) / 2)) << endl;
        else
          out << (int)(255 * (-n.y.mid() + 0.25) / 1.25) << " "
              << (int)(255 * (-n.y.mid() + 0.25) / 1.25) << " "
              << (int)(255 * (-n.y.mid() + 0.25) / 1.25) << endl;
      }
    }
  }

  out.close();
}

int main(int argc, char* argv[]) {
  acp::enable();
  vector<PTR<Triangle>> ts = readTriangles(argv[1], argv[2]);
  vector<TrianglePair> pairs = reflects(ts);
  display("ncolor.ppm", ts, pairs);
  acp::disable();
  return 0;
}
