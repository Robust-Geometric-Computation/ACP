#include "hull.h"
#include "geometry3d.h"
#include "mink.h"
#include "io.h"
#include <cstring>

class InputParameter : public Object<Scalar> {
public:
  InputParameter (double x) : Object(Scalar<double>(x), true) {}
};

class OuterApproxFace;

class SinCosAlpha : public Point {
  PTR<Object<Scalar>> tan_theta;
  DeclareCalculate (PV3) {
    N t = tan_theta->get<N>().x, sint = 2*t/(1+t*t), cost = (1-t*t)/(1+t*t),
      alpha = (1-cost)/sint;
    return PV3<N>(sint, cost, alpha);
  }
public:
  SinCosAlpha (Object<Scalar> *t) : tan_theta(t) {}
};

PTR<Polyhedron> loadPoly_VTK(const char * filename, bool perturbed=false);
void savePoly_VTK(PTR<Polyhedron> p, const char * filename);

class FreeSpace {
public:
  PTR<Polyhedron> robot, obstacle;
  std::vector<PTR<Polyhedron>> blockspaces_close, blockspaces_rough;
  int numRotations;
  bool inner_approximation;
  FreeSpace(PTR<Polyhedron> robot, PTR<Polyhedron> obstacle, PTR<Object<Scalar> > tan_half_angle, int numRotations, bool inner_approximation = true);
  PTR<Polyhedron> generateSweep(std::vector<OuterApproxFace> & trapezoids);
  void generateFreeSpaces(std::vector<PTR<Polyhedron>> & spaces, PTR<Polyhedron> sweep, const char * basename);
};

double tanHalfAngle (int n);
