#ifndef RESULTANT_H
#define RESULTANT_H

#include "acp/encasement2d/predicates2d.h"
#include "acp/poly/poly.h"

using namespace acp;
using namespace std;

/*
 * Exception for when someone creates a Root2 that is not an actual intersection
 * of f and g
 */

class NoIntersectionException : public std::exception {
 public:
  virtual const char* what() const throw() {
    return "Not an actual intersection point";
  }
};

extern NoIntersectionException noIntersectionException;

/*
 * Basic Input parameter
 */

class InputParameter : public Object<Scalar> {
 public:
  InputParameter(const Parameter& p) : Object<Scalar>(p) {}
};

template <class N>
class SylvesterMatrix {
 public:
  // give the determinant of this matrix
  Poly<N> det();

  Poly<N> det4x4();

  Poly<N> get(int i, int j);

  SylvesterMatrix(Poly2<N> f, Poly2<N> g) : f(f), g(g) {
    construct();
    // print();
  }

  void print();

 private:
  int degX(Poly2<N> f);
  vector<Poly<N>> fillVector(Poly2<N> f);

  void construct();
  void createMatrix();

  // give row index of poly of min degree in column j (restricted to rows <= j)
  int findMinDeg(int j);

  // swap row i and row j
  void swapRow(int i, int j);

  // reduce row i2 with i1, (i2, j)'s degree will be reduced by 1 to
  // symbolically reduce the degree
  void reduce(int j, int i1, int i2);

  // count the number of deg -1 polys in this column (restricted to rows <= j)
  int countNegDeg(int j);

  Poly2<N> f;
  Poly2<N> g;
  vector<Poly<N>> fv;
  vector<Poly<N>> gv;
  int fd;                     // degX(f)
  int gd;                     // degX(g)
  int n;                      // deg(f) + deg(g)
  vector<vector<Poly<N>>> m;  // the matrix
};

template <class N>
class ResultantEuclid {
  Poly2<N> f, g;

 public:
  ResultantEuclid(Poly2<N> f, Poly2<N> g) : f(f), g(g) {}

  N getRi(int i);
  int degree(vector<Poly<N>>& f);
  N resultant(vector<N>& f, vector<N>& g);
  Poly<N> resultant(vector<Poly<N>>& f, vector<Poly<N>>& g);

  vector<Poly<N>> fillVector(Poly2<N> f);

  Poly<N> det();
};

/*
 * Returns determinant of the sylvester matrix where coeficients are polys in y
 * (second component)
 */

class Resultant : public Object<Poly> {
 public:
  Resultant(PTR<Object<Poly2>> f, PTR<Object<Poly2>> g) : f(f), g(g) {}
  DeclareCalculate(Poly);

  static bool useEuclid;

 private:
  PTR<Object<Poly2>> f;
  PTR<Object<Poly2>> g;
};

/*
 * Substitute a value of x into a Poly2<N> to induce a Poly<N> in y
 */

class SubstitutePoly1 : public Object<Poly> {
 public:
  SubstitutePoly1(PTR<Object<Scalar>> x, PTR<Object<Poly2>> f, int i)
      : x(x), f(f) {
    this->i = (i == 0 ? 0 : 1);
  }
  DeclareCalculate(Poly);

 private:
  PTR<Object<Scalar>> x;
  PTR<Object<Poly2>> f;
  int i;
};

/*
 * A result of the resultant root finding procedure
 */

class ResultantRoot : public Object<PV2> {
 public:
  ResultantRoot(PTR<Object<Scalar>> x, PTR<Object<Scalar>> y,
                PTR<Object<Poly2>> f, PTR<Object<Poly2>> g)
      : x(x), y(y), f(f), g(g), verified(-1) {}

  DeclareCalculate(PV2);

 private:
  PTR<Object<Scalar>> x;
  PTR<Object<Scalar>> y;
  PTR<Object<Poly2>> f;
  PTR<Object<Poly2>> g;
  int verified;
};
/*
//sign() == -1 if a < b
Primitive2(LessThan, PTR<Object<Scalar>> , a, PTR<Object<Scalar>> , b);

class Ascending {
public:
  bool operator()(PTR<Object<Scalar>>  const a, PTR<Object<Scalar>>  const b) {
    return LessThan(a, b) == -1;
  }
};

class Descending {
public:
  bool operator()(PTR<Object<Scalar>>  const a, PTR<Object<Scalar>>  const b) {
    return LessThan(a, b) == 1;
  }
};
*/

#endif
