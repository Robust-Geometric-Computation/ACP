#ifndef ANGLE_H
#define ANGLE_H
//#include "object.h"
//#include "pv.h"
#include "poly2.h"
#include "contact.h"
#include "assert.h"

#ifndef NOTIME
#include <sys/time.h>
double getTimeX ();
extern double rootTime;
extern double rootHighTime;
extern double predTime;
extern double highTime;
extern double idenTime;
#endif

extern int nGCD53, nGCD106, nGCD100, nGCD212, nRootCmps;

//#define TABLE_IDENTITY
//#define GCD_IDENTITY
//#define RESIDUE_IDENTITY1
//#define RESIDUE_IDENTITY2
//#define HOMOTOPY_IDENTITY

#ifdef GCD_IDENTITY
    extern vector<Contact> gcdPairs;
#endif

namespace acp {

class RootAngle;

// An angle.  Equivalent to point on unit circle.
class Angle : public Object<PV2> {
public:
  enum Type { ZERO, INPUT, EXPLICIT, IMPLICIT, MEAN };
  virtual Type getType () = 0;
  virtual ~Angle () {}
  virtual bool identical (Angle *that) { return this == that; }

  template<class N>
  N getS (bool &flag) { PV2<N> u = get<N>(); return u.y / (1 + u.x); }

  template<class N>
  static PV2<N> s2u (const Scalar<N> &s) {
    N d = 1 + s.x * s.x;
    return PV2<N>((1 - s.x * s.x) / d, s.x * 2 / d);
  }

  class Y : public Primitive {
    Angle *a;
    DeclareSign { return a->get<N>().y; }
  public:
    Y (Angle *a) : a(a) {}
  };

  class Cross : public Primitive {
    Angle *a, *b;
    DeclareSign { return a->get<N>().cross(b->get<N>()); }
  public:
    Cross (Angle *a, Angle *b) : a(a), b(b) {}
  };

  // -1 if 0 <= a < b < 2 PI
  // 0 if equal or identical
  static int compare (Angle *a, Angle *b);

  // 1 if going from a to c through b is counterclockwise
  // 0 if a same as b or b same as c
  static int circulation (Angle *a, Angle *b, Angle *c);

  class CompareGCD : public Primitive {
    RootAngle* a;
    RootAngle* b;
    DeclareSign;
  public:
    CompareGCD (RootAngle* a, RootAngle* b) : a(a), b(b) {}
  };

private:
  void foo2 ();
  void foo3 ();
};

class ZeroAngle : public Angle {
  DeclareCalculate(PV2) {
    return PV2<N>(N(1), N(0));
  }
public:
  Type getType () { return ZERO; }
  bool identical (Angle *that) { return that->getType() == ZERO; }
};

extern ZeroAngle *angle0;

class InputAngle : public Angle {
  static double tanhalf (double x) {
    if (!isEnabled())
      return tan(x/2);
    disable();
    double th = tan(x/2);
    enable();
    return th;
  }
  PTR<Object<Scalar>> s;
  DeclareCalculate(PV2) { return s2u(s->get<N>()); }
public:
  virtual Type getType () { return INPUT; }
  // Create from angle in radians (approximately).
  InputAngle (double angle) : s(new Object<Scalar>(tanhalf(angle), false)) {}
};

class AnglePoly : public Object<Poly2> {
  std::vector<PTR<Angle>> getAngles () {
    return getApprox(1).d == 1 ? getExplicitAngles() : getImplicitAngles();
  }
public:
  template<class N>
  N value (Angle *angle) { return get<N>().value(angle->get<N>()); }
  virtual bool identical (AnglePoly *that) { return this == that; }
  std::vector<PTR<Angle>> getAnglesSorted () {
    std::vector<PTR<Angle>> angles = getAngles();
    for (int i = 0; i < angles.size(); i++) {
      PTR<Angle> angle = angles[i];
      int j = i;
      while (j > 0 && Angle::compare(angle, angles[j-1]) < 0) {
	angles[j] = angles[j-1];
	j--;
      }
      angles[j] = angle;
    }
    for (int i = 1; i < angles.size(); i++)
      assert(Angle::compare(angles[i-1], angles[i]) < 0);
    return angles;
  }
  enum Type { EF, EF2, VF, EE, EEE, EFF, FFF, NNN, FFFF, U };
  virtual Type getType () = 0;
  virtual void getContact (Contact &contact) = 0;
  Contact getContact () {
    Contact contact;
    getContact(contact);
    contact.normalize();
    return contact;
  }
  template<class N>
  Poly<N> univariate ();
private:
  class Discriminant : public Primitive {
    AnglePoly* poly;
    DeclareSign;
  public:
    Discriminant (AnglePoly* poly) : poly(poly) {}
  };
  std::vector<PTR<Angle>> getExplicitAngles ();
  std::vector<PTR<Angle>> getImplicitAngles ();
};

class AnglePoly1 : public Object<Poly> {
  PTR<AnglePoly> poly2;
public:
  AnglePoly1 (AnglePoly *poly2) : poly2(poly2) {}

  DeclareCalculate(Poly);
};

class RootAngle : public Angle {
  friend class AnglePoly;
  PTR<AnglePoly> poly;
protected:
  const int index;
public:
  RootAngle (AnglePoly *poly, int index) : poly(poly), index(index) {}

  AnglePoly *getPoly () { return poly; }

  bool identical (Angle *thatAngle) {
    if (Angle::identical(thatAngle)) return true;
    RootAngle *that = dynamic_cast<RootAngle*>(thatAngle);
    if (that == 0) return false;
#ifdef BLEEN
    gcdPairs.push_back(this->poly->getContact());
    gcdPairs.push_back(that->poly->getContact());
#endif
    if (this->index != that->index) return false;
    if (this->poly->identical(that->poly)) return true;
    return this->poly->getContact() == that->poly->getContact();
  }

  class Sign : public Primitive {
    AnglePoly* poly;
    Angle* root;
    DeclareSign;
  public:
    Sign (AnglePoly* poly, Angle* root) : poly(poly), root(root) {}
  };
};

class ExplicitAngle : public RootAngle {
  friend class AnglePoly;
  // For a x + b y + c = 0, flag false/true means smaller/larger s.
  ExplicitAngle (AnglePoly *poly, bool flag) : RootAngle(poly, flag) {}
  ~ExplicitAngle () {}
  DeclareCalculate(PV2);
public:
  virtual Type getType () { return EXPLICIT; }
};

class ImplicitAngle : public RootAngle {
  friend class AnglePoly;
  PTR<Object<Scalar>> root;
  ImplicitAngle (AnglePoly *poly, Object<Scalar> *root, int index) : RootAngle(poly, index), root(root) {}
  DeclareCalculate(PV2) { return s2u(root->get<N>()); }
public:
  virtual Type getType () { return IMPLICIT; }
};

class MeanAngle : public Angle {
  PTR<Angle> u, v;
  double t;
  DeclareCalculate(PV2) {
    PV2<N> m = u->get<N>() * t + v->get<N>() * (1-t);
    return(m / m.dot(m).sqrt());
  }
public:
  Type getType () { return MEAN; }
  MeanAngle (Angle *u, Angle *v, double t) : u(u), v(v), t(t) {}
};

}

#endif
