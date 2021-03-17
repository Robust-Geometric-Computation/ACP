#include "acp/circle/circle.h"

template <class N>
N CircleDist::calculate() {
  Circle<N> cir = circle->get<N>();
  PV2<N> c = cir.getCenter();
  PV2<N> p = point->get<N>();
  PV2<N> cp = p - c;
  N r2 = cir.getRadius2();
  return (cp.dot(cp) - r2);
}

template Parameter CircleDist::calculate<Parameter>();
template PParameter CircleDist::calculate<PParameter>();
template MParameter CircleDist::calculate<MParameter>();

template <class N>
Circle<N> Circle2::calculate() {
  PV2<N> ap = a->get<N>();
  PV2<N> bp = b->get<N>();
  PV2<N> c = (ap + bp) * 0.5;
  PV2<N> ab = bp - ap;
  N r2 = ab.dot(ab) * 0.25;
  return Circle<N>(c, r2);
}

template <class N>
Circle<N> Circle3::calculate() {
  PV2<N> ap = a->get<N>();
  PV2<N> bp = b->get<N>();
  PV2<N> cp = c->get<N>();
  PV2<N> abm = (ap + bp) * 0.5;
  PV2<N> ab = bp - ap;
  PV2<N> ab90(-ab.y, ab.x);
  PV2<N> bcm = (bp + cp) * 0.5;
  PV2<N> bc = cp - bp;
  // (abm + ab90 t - bcm) * bc = 0
  // abm - bcm = (a + b)/2 - (b + c)/2 = (a - c)/2
  PV2<N> ca2 = (ap - cp) * 0.5;
  // (ca2 + ab90 t) * bc = 0
  // ab90 * bc t = -ca2 * bc
  // t = -(ca2 * bc) / (ab90 * bc)
  N t = -ca2.dot(bc) / ab90.dot(bc);
  PV2<N> center = abm + ab90 * t;
  PV2<N> center2a = ap - center;
  N radius2 = center2a.dot(center2a);

  return Circle<N>(center, radius2);
}

template <class N>
Circle<N> Circle90::calculate() {
  PV2<N> ap = a->get<N>();
  PV2<N> cp = c->get<N>();

  PV2<N> center = (ap + cp) * 0.5;
  PV2<N> ac2 = (cp - ap) * 0.5;

  N radius2 = ac2.dot(ac2);

  return Circle<N>(center, radius2);
}

// List of all circles created by the smallestCircle functions.
vector<PTR<Object<Circle>>> circles;

// Smallest circle containing the first n points in points.
PTR<Object<Circle>> smallestCircle(vector<PTR<Object<PV2>>>& points, int n) {
  PTR<Object<Circle>> circle = new Circle2(points[0], points[1]);
  circles.push_back(circle);

  for (int i = 2; i < n; i++) {
    if (CircleDist(circle, points[i]) > 0)
      circle = smallestCircle(points[i], points, i);
  }

  return circle;
}

// Smallest circle through p containing the first n points in points.
PTR<Object<Circle>> smallestCircle(PTR<Object<PV2>> p,
                                   vector<PTR<Object<PV2>>>& points, int n) {
  PTR<Object<Circle>> circle = new Circle2(p, points[0]);
  circles.push_back(circle);

  for (int i = 1; i < n; i++) {
    if (CircleDist(circle, points[i]) > 0)
      circle = smallestCircle(p, points[i], points, i);
  }

  return circle;
}

// Smallest circle through p and q containing the first n points in points.
PTR<Object<Circle>> smallestCircle(PTR<Object<PV2>> p, PTR<Object<PV2>> q,
                                   vector<PTR<Object<PV2>>>& points, int n) {
  PTR<Object<Circle>> circle = new Circle2(p, q);
  circles.push_back(circle);

  for (int i = 0; i < n; i++) {
    if (CircleDist(circle, points[i]) > 0) {
      circle = new Circle3(p, q, points[i]);
      circles.push_back(circle);
    }
  }

  return circle;
}

// random integer >= 0 and < n.
int randomn(int n) {
  int r = random() % n;
  if (r < 0) r += n;
  return r;
}

void randomize(vector<PTR<Object<PV2>>>& points) {
  for (int i = points.size() - 1; i > 0; i--) {
    int j = randomn(i + 1);
    PTR<Object<PV2>> temp = points[i];
    points[i] = points[j];
    points[j] = temp;
  }
}
