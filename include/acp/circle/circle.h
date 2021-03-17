#ifndef CIRCLE_H
#define CIRCLE_H

#include "acp/encasement2d/predicates2d.h"

using namespace std;
using namespace acp;

// Pure parameter circle class
template <class N>
class Circle {
 public:
  PV2<N> center;
  N radius2;

  PV2<N> getCenter() { return center; }
  N getRadius2() { return radius2; }

  Circle() {}
  Circle(const PV2<N>& center, const N& radius2)
      : center(center), radius2(radius2) {}

  template <class T>
  Circle(const Circle<T>& c) : center(c.center), radius2(c.radius2) {}

  int size() const { return 3; }
  N& operator[](int i) { return i < 2 ? center[i] : radius2; }
  const N& operator[](int i) const { return i < 2 ? center[i] : radius2; }
};

Primitive2(CircleDist, PTR<Object<Circle>>, circle, PTR<Object<PV2>>, point)

    // Circle whose diameter is ab.
    class Circle2 : public Object<Circle> {
  PTR<Object<PV2>> a, b;
  DeclareCalculate(Circle);

 public:
  Circle2(PTR<Object<PV2>> a, PTR<Object<PV2>> b) : a(a), b(b) {}
};

// Circumcircle of triangle abc.
class Circle3 : public Object<Circle> {
  PTR<Object<PV2>> a, b, c;
  DeclareCalculate(Circle);

 public:
  Circle3(PTR<Object<PV2>> a, PTR<Object<PV2>> b, PTR<Object<PV2>> c)
      : a(a), b(b), c(c) {}
};

// Circumcircle of right triangle abc. Angle abc needs to be 90 degrees.
class Circle90 : public Object<Circle> {
  PTR<Object<PV2>> a, b, c;
  DeclareCalculate(Circle);

 public:
  Circle90(PTR<Object<PV2>> a, PTR<Object<PV2>> b, PTR<Object<PV2>> c)
      : a(a), b(b), c(c) {}
};

extern vector<PTR<Object<Circle>>> circles;

PTR<Object<Circle>> smallestCircle(vector<PTR<Object<PV2>>>& points, int n);

PTR<Object<Circle>> smallestCircle(PTR<Object<PV2>> p,
                                   vector<PTR<Object<PV2>>>& points, int n);

PTR<Object<Circle>> smallestCircle(PTR<Object<PV2>> p, PTR<Object<PV2>> q,
                                   vector<PTR<Object<PV2>>>& points, int n);

void randomize(vector<PTR<Object<PV2>>>& points);

#endif
