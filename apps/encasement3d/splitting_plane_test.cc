#ifdef OSX
#include <Accelerate/Accelerate.h>
#else
extern "C" {
#include <clapack.h>
}
#endif

#include <algorithm>
#include <vector>

#include "acp/encasement2d/encasement2d.h"
#include "acp/encasement3d/encasement3d.h"
#include "acp/encasement3d/encasement_utils.h"
#include "acp/linmath/pv.h"
#include "acp/poly/poly3d.h"

namespace acp {
namespace splitting_plane_test {

// random point in a box
PV3<double> RandomPoint(const PV3<double>& box_l, const PV3<double>& box_u) {
  return PV3<double>(randomNumber(box_l.x, box_u.x),
                     randomNumber(box_l.y, box_u.y),
                     randomNumber(box_l.z, box_u.z));
}

// random points in given box
std::vector<PV3<double>> RandomPoints(const PV3<double>& box_l,
                                      const PV3<double>& box_u, const int n) {
  std::vector<PV3<double>> v(n);
  std::generate(v.begin(), v.end(),
                [&]() { return RandomPoint(box_l, box_u); });
  return v;
}

// random surface of given degree
Poly3D<Parameter> RandomSurface(int degree, const PV3<double>& box_l,
                                const PV3<double>& box_u) {
  int num_points = ((degree + 1) * (degree + 2) * (degree + 3) / 6) - 1;

  std::vector<PV3<double>> points = RandomPoints(box_l, box_u, num_points);

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  std::transform(points.begin(), points.end(), std::back_inserter(x),
                 [&x](const PV3<double>& p) { return p.x; });
  std::transform(points.begin(), points.end(), std::back_inserter(y),
                 [&y](const PV3<double>& p) { return p.y; });
  std::transform(points.begin(), points.end(), std::back_inserter(z),
                 [&z](const PV3<double>& p) { return p.z; });

  double m[num_points * num_points];
  double b[num_points];
  int ipiv[num_points];

  // setup the matrix and rhs for interpolation
  for (int col = 0; col < num_points; col++) {
    int c = 0;

    std::vector<double> x_pow(degree + 1);
    std::vector<double> y_pow(degree + 1);
    std::vector<double> z_pow(degree + 1);
    x_pow[0] = y_pow[0] = z_pow[0] = 1;

    for (int i = 1; i <= degree; i++) {
      x_pow[i] = x_pow[i - 1] * x[col];
      y_pow[i] = y_pow[i - 1] * y[col];
      z_pow[i] = z_pow[i - 1] * z[col];
    }

    for (int i = 0; i <= degree; i++) {
      for (int j = 0; j <= degree - i; j++) {
        for (int k = 0; k <= degree - i - j; k++) {
          // std::cout << i << " " << j << " " << k;
          if (i == 0 && j == 0 && k == 0) {
            b[col] = -1;
            // std::cout << " -1" << std::endl;
          } else {
            m[c++ * num_points + col] = x_pow[i] * y_pow[j] * z_pow[k];
            // std::cout << " m[" << (c - 1) << "*" << num_points << " + " << k
            //          << "] ";
            // std::cout << x_pow[i] * y_pow[j] * z_pow[k] << std::endl;
          }
        }
      }
    }
  }

#ifdef OSX
  __CLPK_integer info;
#else
  int info;
#endif

#ifdef OSX
  __CLPK_integer n = num_points;
  __CLPK_integer nrhs = 1;
  __CLPK_integer lda = num_points;
  __CLPK_integer ipivm[num_points];
  __CLPK_integer ldb = num_points;
  dgesv_(&n, &nrhs, m, &lda, ipivm, b, &ldb, &info);
#else
  info = clapack_dgesv(CblasColMajor, num_points, 1, m, num_points, ipiv, b,
                       num_points);
#endif

  for (int i = 0; i < num_points; i++) {
    // std::cout << "b[" << i << "]: " << b[i] << std::endl;
  }

  Poly3D<Parameter> poly(degree, degree, degree);

  int c = 0;

  poly.set(0, 0, 0, Parameter::constant(1));

  for (int i = 0; i <= degree; i++) {
    for (int j = 0; j <= degree - i; j++) {
      for (int k = 0; k <= degree - i - j; k++) {
        if (i != 0 || j != 0 || k != 0) {
          poly.set(i, j, k, Parameter::constant(b[c++]));
        }
      }
    }
  }

  return poly;
}

vector<PTR<Object<PV2>>> boundary_intersections(PTR<Object<Poly2>> f,
                                                PTR<Object<Poly2>> g,
                                                PTR<Object<PV2>> bottom_left,
                                                PTR<Object<PV2>> top_right) {
  Encasement* enc = new Encasement(bottom_left, top_right, rand());
  enc->init();

  std::vector<PTR<Object<Poly2>>> curves = {f, g};
  enc->calculate(curves, false);

  return enc->getRoots();
}

// f3 and g3 are the surfaces
// 0 <= xyz <= 2 is the selection of which axis is constant
// 0 <= low_or_high <= 1 selects the lower or upper box boundary on the xyz axis
// values contains [x_low, x_high, y_low, y_high, z_low, z_high]
vector<PV3<Parameter>> boundary_intersections(
    PTR<Object<Poly3D>> f3, PTR<Object<Poly3D>> g3, int xyz, int low_or_high,
    std::vector<PTR<Object<Scalar>>> values) {
  PTR<Object<Scalar>> value = values[2 * xyz + low_or_high];

  PTR<Object<Poly2>> f2 = new Sub1Poly3D(f3, xyz, value);
  PTR<Object<Poly2>> g2 = new Sub1Poly3D(g3, xyz, value);

  int j = (xyz + 1) % 3;
  int k = (xyz + 2) % 3;

  PTR<Object<PV2>> minXY = new SS(values[2 * j], values[2 * k]);
  PTR<Object<PV2>> maxXY = new SS(values[2 * j + 1], values[2 * k + 1]);

  vector<PTR<Object<PV2>>> b = boundary_intersections(f2, g2, minXY, maxXY);

  vector<PV3<Parameter>> intersections;

  std::transform(b.begin(), b.end(), std::back_inserter(intersections),
                 [&](const PTR<Object<PV2>>& p) {
                   PV2<Parameter> pp = p->getApprox();
                   PV3<Parameter> v;
                   v[xyz] = value->getApprox().x;
                   v[j] = pp.x;
                   v[k] = pp.y;
                   return v;
                 });

  return intersections;
}

class ABdotN : public Primitive {
  Object<PV3> *a, *b, *n;
  DeclareSign {
    return (b->get<N>() - a->get<N>()).dot(n->get<N>());
  }
public:
  ABdotN (Object<PV3> *a, Object<PV3> *b, Object<PV3> *n)
    : a(a), b(b), n(n) {}
};

class ABdotM_CDdotN : public Primitive {
  Object<PV3> *a, *b, *m, *c, *d, *n;
  DeclareSign {
    PV3<N> mv = m->get<N>();
    PV3<N> nv = n->get<N>();
    N abDm = (b->get<N>() - a->get<N>()).dot(mv);
    N cdDn = (c->get<N>() - d->get<N>()).dot(nv);
    // abDm^2/m^2 - cdDn^2/n^2
    // abDm^2 * n^2 - cdDn^2 * m^2
    return (abDm * abDm) * nv.dot(nv) - (cdDn * cdDn) * mv.dot(mv);
  }
public:
  ABdotM_CDdotN (Object<PV3> *a, Object<PV3> *b, Object<PV3> *m,
		 Object<PV3> *c, Object<PV3> *d, Object<PV3> *n)
    : a(a), b(b), m(m), c(c), d(d), n(n) {}
};

class ABx : public Primitive {
  Object<PV3> *a, *b;
  DeclareSign { return b->get<N>().x - a->get<N>().x; }
public:
  ABx (Object<PV3> *a, Object<PV3> *b) 
    : a(a), b(b) {}
};

class ABCz : public Primitive {
  Object<PV3> *a, *b, *c;
  DeclareSign { 
    PV3<N> ap = a->get<N>();
    PV3<N> bp = b->get<N>();
    PV3<N> cp = c->get<N>();
    return (bp.x - ap.x) * (cp.y - ap.y) - (bp.y - ap.y) * (cp.x - ap.x);
  }
public:
  ABCz (Object<PV3> *a, Object<PV3> *b, Object<PV3> *c) 
    : a(a), b(b), c(c) {}
};

class ABCD : public Primitive {
  Object<PV3> *a, *b, *c, *d;
  DeclareSign { 
    PV3<N> ap = a->get<N>();
    PV3<N> bp = b->get<N>();
    PV3<N> cp = c->get<N>();
    PV3<N> dp = d->get<N>();
    return (bp - ap).cross(cp - ap).dot(dp - ap);
  }
public:
  ABCD (Object<PV3> *a, Object<PV3> *b, Object<PV3> *c, Object<PV3> *d) 
    : a(a), b(b), c(c), d(d) {}
};

class ABCdotN : public Primitive {
  Object<PV3> *a, *b, *c, *n;
  DeclareSign { 
    PV3<N> ap = a->get<N>();
    PV3<N> bp = b->get<N>();
    PV3<N> cp = c->get<N>();
    PV3<N> np = c->get<N>();
    return (bp - ap).cross(cp - ap).dot(np);
  }
public:
  ABCdotN (Object<PV3> *a, Object<PV3> *b, Object<PV3> *c, Object<PV3> *n) 
    : a(a), b(b), c(c), n(n) {}
};

class AB2_CD2 : public Primitive {
  Object<PV3> *a, *b, *c, *d;
  DeclareSign { 
    PV3<N> ab = b->get<N>() - a->get<N>();;
    PV3<N> cd = d->get<N>() - c->get<N>();;
    return ab.dot(ab) - cd.dot(cd);
  }
public:
  AB2_CD2 (Object<PV3> *a, Object<PV3> *b, Object<PV3> *c, Object<PV3> *d) 
    : a(a), b(b), c(c), d(d) {}
};

class ABdotAC : public Primitive {
  Object<PV3> *a, *b, *c;
  DeclareSign { 
    PV3<N> ap = a->get<N>();
    PV3<N> bp = b->get<N>();
    PV3<N> cp = c->get<N>();
    return (bp - ap).dot(cp - ap);
  }
public:
  ABdotAC (Object<PV3> *a, Object<PV3> *b, Object<PV3> *c) 
    : a(a), b(b), c(c) {}
};

class AB : public Object<PV3> {
  PTR<Object<PV3>> a, b;
  DeclareCalculate(PV3) {
    return b->get<N>() - a->get<N>();
  }
public:
  AB (PTR<Object<PV3>> a, PTR<Object<PV3>> b)
    : a(a), b(b) {}
};

class ABC : public Object<PV3> {
  PTR<Object<PV3>> a, b, c;
  DeclareCalculate(PV3) {
    PV3<N> ap = a->get<N>();
    PV3<N> bp = b->get<N>();
    PV3<N> cp = c->get<N>();
    return (bp - ap).cross(cp - ap);
  }
public:
  ABC (PTR<Object<PV3>> a, PTR<Object<PV3>> b, PTR<Object<PV3>> c)
    : a(a), b(b), c(c) {}
};

class ABxCD : public Object<PV3> {
  PTR<Object<PV3>> a, b, c, d;
  DeclareCalculate(PV3) {
    PV3<N> ap = a->get<N>();
    PV3<N> bp = b->get<N>();
    PV3<N> cp = c->get<N>();
    PV3<N> dp = d->get<N>();
    return (bp - ap).cross(dp - cp);
  }
public:
  ABxCD (PTR<Object<PV3>> a, PTR<Object<PV3>> b, PTR<Object<PV3>> c, PTR<Object<PV3>> d)
    : a(a), b(b), c(c), d(d) {}
};

// Shortest vector connecting line AB to point C
class ABtoC : public Object<PV3> {
  PTR<Object<PV3>> a, b, c;
  DeclareCalculate(PV3) {
    PV3<N> ap = a->get<N>();
    PV3<N> ab = b->get<N>() - ap;
    PV3<N> ac = c->get<N>() - ap;
    // (c - (a + t ab)) * ab = 0
    // (ac - t ab) * ab = 0
    // ac * ab = t ab * ab
    // t = (ac * ab) / (ab * ab)
    // n = ac - ab t
    // n = ac - ab ((ac * ab) / (ab * ab))
    return  ac - ab * (ac.dot(ab)) / (ab.dot(ab));
  }
public:
  ABtoC (PTR<Object<PV3>> a, PTR<Object<PV3>> b, PTR<Object<PV3>> c)
    : a(a), b(b), c(c) {}
};

class MidPoint : public Object<PV3> {
  PTR<Object<PV3>> a, b;
  DeclareCalculate(PV3) {
    return (a->get<N>() + b->get<N>()) / 2;
  }
public:
  MidPoint (PTR<Object<PV3>> a, PTR<Object<PV3>> b)
    : a(a), b(b) {}
};
  
class Triangle {
public:
  int p[3], t[3];
  PTR<Object<PV3>> n;

  Triangle () { p[0] = p[1] = p[2] = -1; t[0] = t[1] = t[2] = -1; }
  Triangle (const vector<PTR<Object<PV3>>> &ps, int a, int b, int c) { 
    p[0] = a; p[1] = b; p[2] = c;
    t[0] = t[1] = t[2] = -1;
    n = new ABC(ps[a], ps[b], ps[c]);
  }
};
  
void convexHull (const std::vector<PTR<Object<PV3>>> &ps, vector<Triangle> &ts, 
		 map<pair<int, int>, int> &e2t, int a, int b) {
  if (e2t.find(pair<int,int>(a,b)) != e2t.end())
    return;

  int c = 0;
  while (c == a || c == b)
    c++;

  for (int i = 0; i < ps.size(); ++i)
    if (i != a && i != b && i != c && ABCD(ps[a], ps[b], ps[c], ps[i]) > 0)
      c = i;

  for (int i = 0; i < ps.size(); ++i)
    if (i != a && i != b && i != c)
      assert(ABCD(ps[a], ps[b], ps[c], ps[i]) < 0);

  int t = ts.size();
  ts.push_back(Triangle(ps, a, b, c));
  e2t[pair<int,int>(a,b)] = t;
  e2t[pair<int,int>(b,c)] = t;
  e2t[pair<int,int>(c,a)] = t;

  convexHull(ps, ts, e2t, c, b);
  convexHull(ps, ts, e2t, a, c);
}

vector<Triangle> convexHull (const std::vector<PTR<Object<PV3>>> &ps) {
  vector<Triangle> ts;
  map<pair<int, int>, int> e2t;

  int a = 0;
  for (int i = 1; i < ps.size(); ++i)
    if (ABx(ps[a], ps[i]) > 0)
      a = i;

  for (int i = 0; i < ps.size(); ++i)
    if (i != a)
      assert(ABx(ps[a], ps[i]) < 0);

  int b = a == 0 ? 1 : 0;
  for (int i = 0; i < ps.size(); ++i)
    if (i != a && i != b && ABCz(ps[a], ps[b], ps[i]) > 0)
      b = i;

  for (int i = 0; i < ps.size(); ++i)
    if (i != a && i != b)
      assert(ABCz(ps[a], ps[b], ps[i]) < 0);

  convexHull(ps, ts, e2t, a, b);

  for (Triangle &t : ts) {
    for (int i = 0; i < 3; ++i)
      t.t[i] = e2t[pair<int,int>(t.p[(i+2)%3],t.p[(i+1)%3])];
  }

  return ts;
}

// Is the nearest point to v on the hull on edge ab?
// nLeft and nRight are outward normals left and right of ab.
// s = 1 means a,b from fs and v from gs.
// If so, check if distance is larger than from fMax to gMax along nMax,
// and update fMax, etc.
void checkEV (PTR<Object<PV3>> a, PTR<Object<PV3>> b, PTR<Object<PV3>> v, 
	      PTR<Object<PV3>> nLeft, PTR<Object<PV3>> nRight, int s,
	      PTR<Object<PV3>> &nMax, PTR<Object<PV3>> &fMax, PTR<Object<PV3>> &gMax, int &sMax) {
  if (ABdotAC(a, b, v) > 0 && ABdotAC(b, a, v) > 0 &&
      ABCdotN(a, v, b, nLeft) > 0 && ABCdotN(a, v, b, nRight) < 0) {
    PTR<Object<PV3>> n = new ABtoC(a, b, v);
    if (fMax == 0 || ABdotM_CDdotN(a, v, n, fMax, gMax, nMax) > 0) {
      nMax = n;
      sMax = s;
      if (s == 1) {
	fMax = a;
	gMax = v;
	cout << "EV" << endl;
      }
      else {
	fMax = v;
	gMax = a;
	cout << "VE" << endl;
      }
    }
  }
}
  
class PolyPlane : public Object<Poly3D> {
  PTR<Object<PV3>> p, n;
  int s;
  DeclareCalculate(Poly3D) {
    Poly3D<N> poly(1, 1, 1);
    PV3<N> nv = n->get<N>() * s;
    PV3<N> pp = p->get<N>();
    poly.set(0, 0, 0, -nv.dot(pp));
    poly.set(1, 0, 0, nv.x);
    poly.set(0, 1, 0, nv.y);
    poly.set(0, 0, 1, nv.z);
    return poly;
   }
public:
  PolyPlane (PTR<Object<PV3>> p, PTR<Object<PV3>> n, int s)
    : p(p), n(n), s(s) {}
};

// TODO
Poly3D<Parameter> splitting_plane(const std::vector<PV3<Parameter>>& f_points,
                                  const std::vector<PV3<Parameter>>& g_points) {
  // Compute a splitting plane that separates all points in f_points from
  // all points in g_points

  vector<PTR<Object<PV3>>> fs;
  for (int i = 0; i < 160 && i < f_points.size(); i += 32)
    fs.push_back(new Object<PV3>(f_points[i]));

  cout << "fs.size() " << fs.size() << endl;
  // fs.erase(fs.begin() + 80, fs.end());

  vector<PTR<Object<PV3>>> gs;
  for (int i = 0; i < 160 && i < g_points.size(); i += 32)
    gs.push_back(new Object<PV3>(g_points[i]));

  cout << "gs.size() " << gs.size() << endl;
  // gs.erase(gs.begin() + 80, gs.end());

  PTR<Object<PV3>> nMax, fMax, gMax;
  int sMax = 1;

  PTR<Object<PV3>> fMin, gMin;
  for (PTR<Object<PV3>> &f : fs)
    for (PTR<Object<PV3>> &g : gs)
      if (fMin == 0 || AB2_CD2(f, g, fMin, gMin) < 0) {
	fMin = f;
	gMin = g;
      }
	
#ifndef BLEEN	
  PV3<Parameter> fp = fMin->getApprox();
  PV3<Parameter> gp = gMin->getApprox();
  cout << "fg " << (gp - fp).length().mid() << endl;
#endif
  nMax = new AB(fMin, gMin);
  bool separated = true;
  for (PTR<Object<PV3>> &f : fs)
    if (f != fMin && ABdotN(fMin, f, nMax) > 0) {
      cout << "not separated f" << endl;
      separated = false;
      break;
    }
  if (separated)
    for (PTR<Object<PV3>> &g : gs)
      if (g != gMin && ABdotN(gMin, g, nMax) < 0) {
	cout << "not separated g" << endl;
	separated = false;
	break;
      }
  if (false && separated) {
    fMax = fMin;
    gMax = gMin;
    sMax = 1;
    cout << "VV" << endl;
  }
  
  vector<Triangle> fhull = convexHull(fs);
  vector<Triangle> ghull = convexHull(gs);

  for (Triangle &t : fhull) {
    PTR<Object<PV3>> n = t.n;
    PTR<Object<PV3>> f = fs[t.p[0]];

    PTR<Object<PV3>> g;
    for (PTR<Object<PV3>> p : gs)
      if (g == 0 || ABdotN(g, p, n) < 0)
	g = p;

    if (ABdotN(f, g, n) > 0 &&
	(fMax == 0 ||
	 ABdotM_CDdotN(f, g, n, fMax, gMax, nMax) > 0)) {
      fMax = f;
      gMax = g;
      nMax = n;
      sMax = 1;
      cout << "FV" << endl;

      PV3<Parameter> fMaxp = fMax->getApprox();
      PV3<Parameter> gMaxp = gMax->getApprox();
      PV3<Parameter> nMaxp = nMax->getApprox() * sMax;
      Parameter fgdotn = (gMaxp - fMaxp).dot(nMaxp);
      Parameter l = ((fgdotn * fgdotn) / nMaxp.dot(nMaxp)).sqrt();
      cout << "l " << l.mid() << endl;
    }

    for (int i = 0; i < 3; i++)
      checkEV(fs[t.p[(i+1)%3]], fs[t.p[(i+2)%3]], g,
	      n, fhull[t.t[i]].n, 1,
	      nMax, fMax, gMax, sMax);
  }
	
  for (Triangle &t : ghull) {
    PTR<Object<PV3>> n = t.n;
    PTR<Object<PV3>> g = gs[t.p[0]];

    for (PTR<Object<PV3>> g2 : gs)
      assert(ABdotN(g, g2, n) <= 0);

    PTR<Object<PV3>> f;
    for (PTR<Object<PV3>> p : fs)
      if (f == 0 || ABdotN(f, p, n) < 0)
	f = p;

    for (PTR<Object<PV3>> f2 : fs)
      assert(ABdotN(f, f2, n) >= 0);

    if (ABdotN(f, g, n) < 0 &&
	(fMax == 0 ||
	 ABdotM_CDdotN(f, g, n, fMax, gMax, nMax) > 0)) {
      fMax = f;
      gMax = g;
      nMax = n;
      sMax = -1;
      cout << "VF" << endl;

      PV3<Parameter> fMaxp = fMax->getApprox();
      PV3<Parameter> gMaxp = gMax->getApprox();
      PV3<Parameter> nMaxp = nMax->getApprox() * sMax;
      Parameter fgdotn = (gMaxp - fMaxp).dot(nMaxp);
      Parameter l = ((fgdotn * fgdotn) / nMaxp.dot(nMaxp)).sqrt();
      cout << "l " << l.mid() << endl;
    }

    for (int i = 0; i < 3; i++)
      checkEV(gs[t.p[(i+1)%3]], gs[t.p[(i+2)%3]], f,
	      n, ghull[t.t[i]].n, -1,
	      nMax, fMax, gMax, sMax);
  }
	
  for (Triangle &tf : fhull) {
    for (int ef = 0; ef < 3; ++ef) {
      if (tf.p[(ef+1)%3] > tf.p[(ef+2)%3])
	continue;
      PTR<Object<PV3>> af = fs[tf.p[(ef+1)%3]];
      PTR<Object<PV3>> bf = fs[tf.p[(ef+2)%3]];
      PTR<Object<PV3>> nfLeft = tf.n;
      PTR<Object<PV3>> nfRight = fhull[tf.t[ef]].n;
      
      for (Triangle &tg : ghull) {
	for (int eg = 0; eg < 3; ++eg) {
	  PTR<Object<PV3>> ag = gs[tg.p[(eg+1)%3]];
	  PTR<Object<PV3>> bg = gs[tg.p[(eg+2)%3]];
	  PTR<Object<PV3>> ngLeft = tg.n;
	  PTR<Object<PV3>> ngRight = ghull[tg.t[eg]].n;
      
	  if (ABdotN(af, bf, ngRight) < 0 &&
	      ABdotN(af, bf, ngLeft) > 0 &&
	      ABdotN(ag, bg, nfRight) < 0 &&
	      ABdotN(ag, bg, nfLeft) > 0) {
	    if (ABCD(af, bf, ag, bg) < 0) {
	      PTR<Object<PV3>> n = new ABxCD(af, bf, ag, bg);
	      assert(ABdotN(af, ag, n) > 0);
	    
	      if (fMax == 0 ||
		  ABdotM_CDdotN(af, ag, n, fMax, gMax, nMax) > 0) {
		fMax = af;
		gMax = ag;
		nMax = n;
		sMax = 1;

		cout << "EE" << endl;

		PV3<Parameter> fMaxp = fMax->getApprox();
		PV3<Parameter> gMaxp = gMax->getApprox();
		PV3<Parameter> nMaxp = nMax->getApprox() * sMax;
		Parameter fgdotn = (gMaxp - fMaxp).dot(nMaxp);
		Parameter l = ((fgdotn * fgdotn) / nMaxp.dot(nMaxp)).sqrt();
		cout << "l " << l.mid() << endl;
	      }
	    }

	    checkEV(af, bf, ag,
		    nfLeft, nfRight, 1,
		    nMax, fMax, gMax, sMax);
	    checkEV(af, bf, bg,
		    nfLeft, nfRight, 1,
		    nMax, fMax, gMax, sMax);
	    checkEV(ag, bg, af,
		    ngLeft, ngRight, -1,
		    nMax, fMax, gMax, sMax);
	    checkEV(ag, bg, bf,
		    ngLeft, ngRight, -1,
		    nMax, fMax, gMax, sMax);
	  }
	}
      }
    }
  }

  if (fMax != 0) {
    cout << "separated" << endl;
    PTR<Object<PV3>> pMax = new MidPoint(fMax, gMax);
    
    PV3<Parameter> fMaxp = fMax->getApprox();
    PV3<Parameter> gMaxp = gMax->getApprox();

    cout << "fg " << (fMaxp - gMaxp).length().mid() << endl;

    PV3<Parameter> nMaxp = nMax->getApprox() * sMax;
    Parameter fgdotn = (gMaxp - fMaxp).dot(nMaxp);
    cout << "fgdotn " << fgdotn.mid() << endl;
    Parameter l = ((fgdotn * fgdotn) / nMaxp.dot(nMaxp)).sqrt();
    cout << "l " << l.mid() << endl;

    for (PTR<Object<PV3>> &f : fs)
      assert(sMax * ABdotN(pMax, f, nMax) < 0);
    for (PTR<Object<PV3>> &g : gs)
      assert(sMax * ABdotN(pMax, g, nMax) > 0);

    PTR<Object<Poly3D>> poly = new PolyPlane(pMax, nMax, sMax);
    for (PTR<Object<PV3>> &f : fs)
      assert(Poly3Dval(poly, f) < 0);
    for (PTR<Object<PV3>> &g : gs)
      assert(Poly3Dval(poly, g) > 0);
  }

  return Poly3D<Parameter>();
}

}  // namespace splitting_plane_test
}  // namespace acp

int main(int argc, char** argv) {
  enable();

  if (argc > 2) {
    srand(atoi(argv[2]));
  }

  PV3<double> box_l(-0.5, -0.5, -0.5);
  PV3<double> box_u(0.5, 0.5, 0.5);

  std::vector<PTR<Object<Scalar>>> values = {
      new Object<Scalar>(Parameter::constant(box_l.x)),
      new Object<Scalar>(Parameter::constant(box_u.x)),
      new Object<Scalar>(Parameter::constant(box_l.y)),
      new Object<Scalar>(Parameter::constant(box_u.y)),
      new Object<Scalar>(Parameter::constant(box_l.z)),
      new Object<Scalar>(Parameter::constant(box_u.z))};

  int degree = (argc > 1 ? atoi(argv[1]) : 2);  // default degree 2

  int num_terms = (degree + 1) * (degree + 2) * (degree + 3) / 6;

  int npolys = (argc > 3 ? atoi(argv[3]) : 2);  // default 2 surfaces

  std::vector<PTR<Object<Poly3D>>> polys;
  for (int i = 0; i < npolys; i++) {
    Poly3D<Parameter> poly =
        acp::splitting_plane_test::RandomSurface(degree, box_l, box_u);
    polys.push_back(new Object<Poly3D>(poly));
  }

  // You can use the ^^above^^ to generate random Poly3s, or hard code your own
  vector<PV3<Parameter>> all_intersections;

  PTR<Object<Poly3D>> f = polys[0];
  PTR<Object<Poly3D>> g = polys[1];

  for (int xyz = 0; xyz < 3; xyz++) {
    for (int low_or_high = 0; low_or_high < 2; low_or_high++) {
      vector<PV3<Parameter>> intersections =
          acp::splitting_plane_test::boundary_intersections(
              f, g, xyz, low_or_high, values);

      all_intersections.insert(all_intersections.begin(), intersections.begin(),
                               intersections.end());
    }
  }

  // Only compute splitting plane when there are two FG curve segments
  assert(all_intersections.size() == 4);

  Poly3<Parameter> fp(f->getApprox());
  Poly3<Parameter> gp(g->getApprox());

  PV3<Parameter> box(Parameter::interval(box_l.x, box_u.x),
                     Parameter::interval(box_l.y, box_u.y),
                     Parameter::interval(box_l.z, box_u.z));

  // If all goes well, there should be two things in the map, whose
  // .first in the <key, value> pair has (i, j) where j is non-negative
  // If j is -1, that means FGCurves didn't find a matching boundary
  // endpoint when tracing. This will result in more than 2 curves
  // being traced out.
  std::map<std::pair<int, int>, vector<PV3<Parameter>>> curves =
      acp::encasement_utils::FGCurves(fp, gp, box, all_intersections);

  // An assert here about how many elements in curves? In case something goes
  // wrong
  assert(curves.size() == 2);

  // Ok if you got here, then there are two curves in curves that connect up
  // unique pairs of boundary points. The ideal case for creating a splitting
  // plane.

  // Compute the splitting plane.
  // TODO

  vector<PV3<Parameter>> points[2];
  int ipoints = 0;
  for (std::pair<std::pair<int, int>, vector<PV3<Parameter>>> curve : curves)
    points[ipoints++] = curve.second;

  acp::splitting_plane_test::splitting_plane(points[0], points[1]);

  std::cout << "fg curves size: " << curves.size() << std::endl;

  disable();
}

int mainDR (int degree, int seed) {
  enable();

  srand(seed);

  PV3<double> box_l(-0.5, -0.5, -0.5);
  PV3<double> box_u(0.5, 0.5, 0.5);

  std::vector<PTR<Object<Scalar>>> values = {
      new Object<Scalar>(Parameter::constant(box_l.x)),
      new Object<Scalar>(Parameter::constant(box_u.x)),
      new Object<Scalar>(Parameter::constant(box_l.y)),
      new Object<Scalar>(Parameter::constant(box_u.y)),
      new Object<Scalar>(Parameter::constant(box_l.z)),
      new Object<Scalar>(Parameter::constant(box_u.z))};

  int num_terms = (degree + 1) * (degree + 2) * (degree + 3) / 6;

  int npolys = 2;

  std::vector<PTR<Object<Poly3D>>> polys;
  for (int i = 0; i < npolys; i++) {
    Poly3D<Parameter> poly =
        acp::splitting_plane_test::RandomSurface(degree, box_l, box_u);
    polys.push_back(new Object<Poly3D>(poly));
  }

  // You can use the ^^above^^ to generate random Poly3s, or hard code your own
  vector<PV3<Parameter>> all_intersections;

  PTR<Object<Poly3D>> f = polys[0];
  PTR<Object<Poly3D>> g = polys[1];

  for (int xyz = 0; xyz < 3; xyz++) {
    for (int low_or_high = 0; low_or_high < 2; low_or_high++) {
      vector<PV3<Parameter>> intersections =
          acp::splitting_plane_test::boundary_intersections(
              f, g, xyz, low_or_high, values);

      all_intersections.insert(all_intersections.begin(), intersections.begin(),
                               intersections.end());
    }
  }

  disable();
  return all_intersections.size();
}

int mainY(int argc, char** argv) {
  int degree = (argc > 1 ? atoi(argv[1]) : 2);  // default degree 2

  for (int seed = 1; seed <= 10; seed++) {
    int nints = mainDR(degree, seed);
    cout << seed << " " << nints << endl;
  }

  return 0;
}
