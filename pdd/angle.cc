#include "angle.h"
#include "polym.h"
#include "root.h"
#include <iomanip>
using namespace std;
using namespace acp;

extern Contact current_contact;
extern Angle *current_angle;
extern int nIdTests, nIdentities;
extern bool commonRoot (Contact contact, RootAngle *root);

bool inGCD;

#ifndef NOTIME
#ifdef BLEEN
double getTimeX ()
{
  timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1e-6*tv.tv_usec;
}
#endif
double rootTime;
double rootHighTime;
double predTime;
double highTime;
double idenTime;
#endif

namespace acp {
  ZeroAngle *angle0 = new ZeroAngle();

template<class N>
N Angle::CompareGCD::calculate () {
  int cp = curPrecision();
  inGCD = false;
  PV2<N> aU = a->get<N>();
  PV2<N> bU = b->get<N>();
  if (aU.y.sign() > 0 && bU.y.sign() < 0)
    return N(-1);
  if (aU.y.sign() < 0 && bU.y.sign() > 0)
    return N(1);
  N bxa = bU.cross(aU);
  int bxas = bxa.sign(false);
  if (bxas != 0)
    return bxa;

  //if (curPrecision() == 53u)
  // nIdTests++;

  bool flag = false;
  N s = b->getS<N>(flag);
  Poly<N> f = a->getPoly()->univariate<N>();
  Poly<N> g = b->getPoly()->univariate<N>();
  
  inGCD = true;
  Poly<N> fder = f.der();
  Poly<N> f2 = f.gcd(fder);
  if (f2.deg() > 0) {
    Poly<N> r;
    f = f.quotient(f2, r);
    assert(r.deg() == -1);
  }

  Poly<N> fg = f.gcd(g);
  if (fg.deg() > 0) {
    Poly<N> r;
    f = f.quotient(fg, r);
    assert(r.deg() == -1);
  }

  inGCD = false;
  N fs = f.value(s);

  nIdTests++;

  if (curPrecision() == 53u)
    nGCD53++;
  else if (curPrecision() == 100u)
    nGCD100++;
  else if (curPrecision() == 106u)
    nGCD106++;
  else if (curPrecision() == 212u)
    nGCD212++;
  else
    cout << "curPrecision " << curPrecision() << endl;

  if (fs.sign(false) != 0) {
    nIdentities++;
    return N(0);
  }

  throw SignException();
  assert(0);
  return N(0);
}

class RootAngleScalar : public Object<Scalar> {
  RootAngle *a;
  DeclareCalculate(Scalar) {
    if (std::is_same<N, PParameter>::value) {
      if (!uninitialized()) {
        Parameter x = getCurrentP().x;
        return N(PParameter(x, true));
      }
      throw SignException(true);
    }
    bool flag = false;
    return a->getS<N>(flag);
  }
public:
  RootAngleScalar (RootAngle *a) : a(a) {}
};

class RootAnglePoly : public Object<Poly> {
  RootAngle *a;
  DeclareCalculate(Poly) {
    return a->getPoly()->univariate<N>();
  }
public:
  RootAnglePoly (RootAngle *a) : a(a) {}
};

class RootAnglePolyM : public Object<PolyM> {
  RootAngle *a;
  DeclareCalculate(PolyM) {
    return PolyM<N>(a->getPoly()->univariate<N>());
  }
public:
  RootAnglePolyM (RootAngle *a) : a(a) {}
};

int Angle::compare (Angle *a, Angle *b) {
  if (a == b)
    return 0;
  if (a == angle0)
    return -1;
  if (b == angle0)
    return 1;

#ifdef GCD_IDENTITY
  RootAngle *ra = dynamic_cast<RootAngle*>(a);
  RootAngle *rb = dynamic_cast<RootAngle*>(b);
  if (ra != 0 && rb != 0)
    return CompareGCD(ra, rb);
#endif

  if (false) {
    RootAngle *ra = dynamic_cast<RootAngle*>(a);
    RootAngle *rb = dynamic_cast<RootAngle*>(b);
    RootAngle *other = ra == current_angle ? rb : ra;
    if (curPrecision() == 53u)
      nIdTests++;
    if (commonRoot(current_contact, other))
      nIdentities++;
  }

  int ya = Y(a);
  int yb = Y(b);
  assert(ya != 0 && yb != 0);
  if (ya > 0 && yb < 0)
    return -1;
  if (ya < 0 && yb > 0)
    return 1;
  int s = a->get<Parameter>().cross(b->get<Parameter>()).sign(false);
  if (s != 0)
    return -s;
  nIdTests++;
#ifdef RESIDUE_IDENTITY1
  RootAngle *ra = dynamic_cast<RootAngle*>(a);
  RootAngle *rb = dynamic_cast<RootAngle*>(b);
  if (ra != 0 && rb != 0) {
    RootAngleScalar sa(ra);
    RootAnglePoly pa(ra);
    RootAnglePolyM eb(rb);
    int sRES = PolyM1Sign(&eb, &sa, &pa);
    // assert(a->identical(b) == (sRES == 0));
    if (sRES == 0) {
      nIdentities++;
      return 0;
    }
  }
#endif
#ifdef RESIDUE_IDENTITY2
  RootAngle *ra = dynamic_cast<RootAngle*>(a);
  RootAngle *rb = dynamic_cast<RootAngle*>(b);
  if (ra != 0 && rb != 0) {
    RootAngleScalar sa(ra), sb(rb);
    RootAnglePoly pa(ra), pb(rb);
    PolyM<double> ep(2);
    vector<int> v(2);
#ifdef BLEEN    
    // x0 - x1
    v[0] = 1;
    v[1] = 0;
    ep.add(v, 1);
    v[0] = 0;
    v[1] = 1;
    ep.add(v, -1);
#endif
    // (1 - x0^2, 2 x0) x (1 - x1^2, 2 x1)
    // x1 - x0^2 x1 - x0 + x0x1^2
    v[0] = 0;
    v[1] = 1;
    ep.add(v, 1);
    v[0] = 2;
    v[1] = 1;
    ep.add(v, -1);
    v[0] = 1;
    v[1] = 0;
    ep.add(v, -1);
    v[0] = 1;
    v[1] = 2;
    ep.add(v, 1);
    Object<PolyM> e(ep, false);
    RootAnglePolyM eb(rb);
    int sRES2 = PolyM2Sign(&e, &sa, &sb, &pa, &pb);
    // int sHOM2 = Cross(a, b);
    // assert(sRES2 == sHOM2);
    if (sRES2 == 0)
      nIdentities++;
    return -sRES2;
  }
#endif
#ifdef TABLE_IDENTITY
  if (a->identical(b)) {
    nIdentities++;
    return 0;
  }
#endif
#ifdef HOMOTOPY_IDENTITY
  RootAngle *ra = dynamic_cast<RootAngle*>(a);
  RootAngle *rb = dynamic_cast<RootAngle*>(b);
  if (ra != 0 && rb != 0) {
    RootAngleScalar sa(ra);
    RootAnglePoly pb(rb);
    int sHOM1 = PolySign(&pb, &sa);
    if (sHOM1 == 0) {
      nIdentities++;
      return 0;
    }
  }
#endif
  s = Cross(a, b);

  if (s == 0)
    nIdentities++;
  return -s;
}

int Angle::circulation (Angle *a, Angle *b, Angle *c) {
  if (a->identical(b) || b->identical(c))
    return 0;
  if (a->identical(c))
    assert(0);
  if (Cross(a, c) > 0) {
    if (Cross(a, b) > 0 && Cross(c, b) < 0)
      return 1;
  }
  else {
    if (Cross(a, b) > 0 || Cross(c, b) < 0)
      return 1;
  }
  return -1;
}

void Angle::foo3 () {
  Angle *a, *b, *c;
  Angle::circulation(a, b, c);
}

template<class N>
PV2<N> ExplicitAngle::calculate () {
  Poly2<N> pp2 = getPoly()->get<N>();
  N &a = pp2.a[1];
  N &b = pp2.a[2];
  N &c = pp2.a[0];
  N r = -c / (a*a + b*b);
  PV2<N> p(a * r, b * r);
  PV2<N> q = PV2<N>(-p.y, p.x) * (1 / p.dot(p) - 1).sqrt();
  if (index == 0)
    return p + q;
  else
    return p - q;
#ifdef BLEEN
  /*
  if (index != (p.x.sign() < 0 && 
                (p.y.sign() > 0 ? 
                 ((p.y + q.y).sign() < 0) : 
                 ((p.y - q.y).sign() > 0))))
    return p + q;
  else
    return p - q;
  */
  PV2<N> ret;
  if (index != (p.x.sign() < 0 && 
                (p.y.sign() > 0 ? 
                 ((p.y + q.y).sign() < 0) : 
                 ((p.y - q.y).sign() > 0))))
    ret = p + q;
  else
    ret = p - q;
  PV2<Parameter> cp = getCurrentP();
  if (!cp.x.uninitialized()) {
    assert(!(cp.x.ub() < ret.x.lb() || ret.x.ub() < cp.x.lb()));
    assert(!(cp.y.ub() < ret.y.lb() || ret.y.ub() < cp.y.lb()));
  }
  return ret;
#endif
}

template<class N>
N AnglePoly::Discriminant::calculate () {
  Poly2<N> pp2 = poly->get<N>();
  N &a = pp2.a[1];
  N &b = pp2.a[2];
  N &c = pp2.a[0];
  return (a*a + b*b - c*c);
}

vector<PTR<Angle>> AnglePoly::getExplicitAngles () {
  vector<PTR<Angle>> angles;
  if (Discriminant(this) > 0) {
    angles.push_back(new ExplicitAngle(this, 0));
    angles.push_back(new ExplicitAngle(this, 1));
  }
  return angles;
}

vector<PTR<Angle>> AnglePoly::getImplicitAngles () {
  vector<PTR<Angle>> angles;
  PTR<AnglePoly1> poly = new AnglePoly1(this);
  Roots roots = PolySolver(poly).getRoots();
  for (int i = 0; i < roots.size(); i++) {
#ifdef BLEEN
    Parameter s = roots[i]->get<Parameter>().x;
    PV2<Parameter> u((1-s*s)/(1+s*s), 2*s/(1+s*s));
    cout << "poly(root) "
	 << s.lb() << " "
	 << u.x.lb() << " " << u.y.lb() << " "
	 << ppoly.value(s).lb() << endl;
#endif
    ImplicitAngle *angle = new ImplicitAngle(this, roots[i], i);
    angles.push_back(angle);
  }
  return angles;
}

int choose (int k, int n)
{
  int res = 1;
  for (int i = 0; i < k; ++i)
    res *= n - i;
  for (int i = 2; i <= k; ++i)
    res /= i;
  return res;
}

int power (int b, int d)
{
  int p = 1;
  for (int i = 0; i < d; ++i)
    p *= b;
  return p;
}

template<class N>
static void term (const N &u, int c, int s, int d, vector<N> &a) {
  int e = d - c - s;
  int power2s = power(2, s);
  for (int k = 0; k <= c + e; ++k) {
    int l = 0;
    int minus1toi = 1;
    for (int i = 0; i <= k; ++i) {
      l += minus1toi * choose(i, c) * choose(k - i, e);
      minus1toi = -minus1toi;
    }
    l *= power2s;
    int p = 2*k + s;
    a[p] = a[p] + u*l;
  }
}

template<class N>
Poly<N> AnglePoly::univariate () {
  Poly2<N> poly = get<N>();
  int d = poly.degree(), dg = 0;
  vector<N> a(2*d+1);
  for (int i = 0; i < 2*d + 1; ++i)
    a[i] = N(0);
  for (int i = 0; i < poly.a.size(); ++i) {
    int c = poly.m[2*i], s = poly.m[2*i+1];
    term(poly.a[i], c, s, d, a);
    dg = max(dg, 2*d - s);
  }
  while (a.size() > dg + 1)
    a.pop_back();
#ifdef BLEEN
  Poly<N> p(a);
  N s = N(-7.47455);
  PV2<N> v((1-s*s)/(1+s*s), 2*s/(1+s*s));
  N polyval = poly.value(v);
  for (int i = 0; i < d; i++)
    polyval = polyval * (1+s*s);
  cout << "univariate "
       << polyval.lb() << " "
       << p.value(s).lb() << endl;
  return p;
#endif
  return Poly<N>(a);
}

template<class N>
Poly<N> AnglePoly1::calculate () {
  Poly<N> g= poly2->univariate<N>();
#ifndef TABLE_IDENTITY
  inGCD = true;
  Poly<N> gder = g.der();
  Poly<N> rem;
  Poly<N> f = g.gcd(gder);
  if (f.degree() > 0) {
    Poly<N> r;
    g = g.quotient(f, r);
  }
  inGCD = false;
#endif
  return g;
}

template<class N>
N RootAngle::Sign::calculate () { 
  return poly->get<N>().value(root->get<N>()); 
}
}

#ifdef BLEEN
Poly2 Unitizer::getPoly () {
  Poly2 poly(2);
  poly.m[0] = 1; poly.m[1] = 0;
  poly.m[2] = 0; poly.m[3] = 1;
  poly.setDegree();
  PV2 xy = v->getV();
  poly.a[0] = -xy.y;
  poly.a[1] = xy.x;
  return poly;
}

int Unitizer::Dot::sign () {
  return v->getV().dot(r->getU()).sign();
}

Angle *Unitizer::unitize (VectorObject *v) {
  Unitizer *u = new Unitizer(v);
  vector<Angle*> roots = u->getRoots();
  if (Unitizer::Dot(v, roots[0]) > 0)
    return roots[0];
  else
    return roots[1];
}

Angle *MeanAngle::mean (Angle *u, Angle *v, double t) {
  MeanAngle *m = new MeanAngle(u, v, t);
  return Unitizer::unitize(m);
}
#endif // BLEEN

