#include <cstdio>
#include <ctime>
#include <iomanip>
#include <map>
#include <set>
#include <unordered_set>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <utility>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <unistd.h>
#include <glm/gtc/quaternion.hpp> 
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/matrix_decompose.hpp>

/**
#include <ft2build.h>
#include FT_FREETYPE_H
**/
#include "acp/encasement2d/encasement2d.h"
#include "acp/encasement3d/encasement_utils.h"
#include "acp/poly/poly3.h"
#include "acp/poly/root.h"
#include "acp/visualization/stb_image.h"
#include "acp/visualization/Playback.h"

/**
#include "acp/ftype/Shader.h"
#include "acp/ftype/TextRenderer.h"
#include "acp/ftype/ResourceManager.h"
**/


typedef PTR<Object<Scalar>> Real;
typedef PTR<Object<Poly3>> POLY3;

bool verbose = false;

bool drawSurfaces = true;
bool drawLines = true;
bool drawBox = true;
bool drawFG = true;

void displayChains(PTR<Object<Poly3>> pol3);
void displayFGCurves(PTR<Object<Poly3>> f, PTR<Object<Poly3>> g);
void setCubeFaces();

int current_curve;
vector<double> curve_min;
vector<double> curve_max;

const char* interpolate_file = "poses.txt";
Playback *playback_object = new Playback(interpolate_file);
void main_interpolate();
void recordCurrentPose(const char* filename);
void vregrid();
void restorePoly();
void updateModel();
void assignInverses();

    std::clock_t bench;
    double currentbench;

int previous_width;
int previous_height;


glm::mat4 trs;
glm::mat4 trsinv;
glm::mat4 pose_model;
glm::mat4 pose_model_inverse;
glm::mat4 P_Inverse;
glm::mat4 V_Inverse;
//void displayRegridding();
//bool dregrid = false;
//bool regridding_on_screen = false;
//TextRenderer *RText;


class vertex {
 public:
  float x, y, z;
  vertex() {}
  vertex(float x0, float y0, float z0) {
    x = x0;
    y = y0;
    z = z0;
  }
  bool operator==(const vertex& v) const {
    if (fabs(x - v.x) < 1e-5 && fabs(y - v.y) < 1e-5 && fabs(z - v.z) < 1e-5)
      return true;

    return false;
  }
  float& operator[](int i) { return i == 0 ? x : i == 1 ? y : z; }
};

vertex vmin;
vertex vmax;

class FXYZ : public Primitive {
  Object<Poly3>* poly3;
  Object<Scalar>* x;
  Object<Scalar>* y;
  Object<Scalar>* z;

 public:
  FXYZ(Object<Poly3>* poly3, Object<Scalar>* x, Object<Scalar>* y,
       Object<Scalar>* z)
      : poly3(poly3), x(x), y(y), z(z) {}
  DeclareSign {
    return poly3->get<N>().value(x->get<N>().x, y->get<N>().x, z->get<N>().x);
  }
};

class TransformedPoly3 : public Object<Poly3> {
  PTR<Object<Poly3>> p;
  double t[16];
  DeclareCalculate(Poly3) { return p->get<N>().transform(t); }

 public:
  TransformedPoly3(PTR<Object<Poly3>> p, glm::mat4& m) : p(p) {
    float* f = (float*)&m;
    for (int i = 0; i < 16; i++) t[i] = f[i];
  }
};

class UnitGradient : public Object<PV3> {
  POLY3 p;
  Real x, y, z;

 public:
  UnitGradient(POLY3 p, Real x, Real y, Real z) : p(p), x(x), y(y), z(z) {}
  DeclareCalculate(PV3) {
    return p->get<N>()
        .gradient(x->get<N>().x, y->get<N>().x, z->get<N>().x)
        .unit();
  }
};

// give this function a matrix and cast it to a float star
// t[16] is a matric in column major form
// making a new ppolu 3 using the current one... make the transformation matrix
// and cast it in as float star

// The univariate polynomial resulting from substituting for two
// variables in a trivariate polynomial
class Sub2Poly : public Object<Poly> {
  POLY3 ppoly3;
  int ixyz;
  Real x, y, z;
  Poly<Parameter> calculate();

 public:
  // ixyz = 0 or 1 or 2 means substituting for y,z or x,z or x,y
  Sub2Poly(POLY3 ppoly3, int ixyz, Real xyz[3])
      : ppoly3(ppoly3), ixyz(ixyz), x(xyz[0]), y(xyz[1]), z(xyz[2]) {}
  DeclareCalculate(Poly) {
    return ppoly3->get<N>().substitute2(ixyz, x->get<N>().x, y->get<N>().x,
                                        z->get<N>().x);
  }
};

// a + (i/n) (b - a)
class IofN : public Object<Scalar> {
  int i, n;
  Real a, b;

 public:
  IofN(int i, int n, Real a, Real b) : i(i), n(n), a(a), b(b) {}
  DeclareCalculate(Scalar) {
    return a->get<N>().x + i * ((b->get<N>().x - a->get<N>().x) / n);
  }
};

Roots edgeRoots(POLY3 poly3, Real xyz[3], int istep, Real step) {
  //   cout << "edgeRoots...";
  PTR<Object<Poly>> poly = new Sub2Poly(poly3, istep, xyz);
  Roots roots = PolySolver(poly).getRoots(xyz[istep], step);
  //   cout << "done" << endl;
  return roots;
}

class Square {
 public:
  // this is the minimum coord of all xyz
  int xyz[3];
  // this is the direction that you do not step in the sqaure (because a square
  // goes in two)
  int istep;
  Square() {}
  Square(int x, int y, int z, int istep) {
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
    this->istep = istep;
  }
};

bool operator==(const Square& lhs, const Square& rhs) {
  return lhs.xyz[0] == rhs.xyz[0] && lhs.xyz[1] == rhs.xyz[1] &&
         lhs.xyz[2] == rhs.xyz[2] && lhs.istep == rhs.istep;
}

bool operator<(const Square& lhs, const Square& rhs) {
  for (int i = 0; i < 3; i++)
    if (lhs.xyz[i] != rhs.xyz[i]) return lhs.xyz[i] < rhs.xyz[i];
  return lhs.istep < rhs.istep;
}

class GridRoots {
 public:
  int nSteps[3];
  Real minXYZ[3];

  void setMinXYZ(Real xyz[3]) {
    for (int i = 0; i < 3; i++) minXYZ[i] = xyz[i];
  }

  Real maxXYZ[3];

  void setMaxXYZ(Real xyz[3]) {
    for (int i = 0; i < 3; i++) maxXYZ[i] = xyz[i];
  }

  vector<Real> allXYZ[3];

  void setSteps(int nSteps[3]) {
    for (int i = 0; i < 3; i++) {
      this->nSteps[i] = nSteps[i];
      allXYZ[i].push_back(minXYZ[i]);
      for (int j = 1; j < nSteps[i]; j++)
        allXYZ[i].push_back(new IofN(j, nSteps[i], minXYZ[i], maxXYZ[i]));
      allXYZ[i].push_back(maxXYZ[i]);
    }
  }

  class EdgeRoots {
   public:
    int xyz[3];  // grid "point" such as 7, 3, 4
    // actual grid point (allXYZ[0][xyz[0]], allXYZ[1][xyz[1]],
    // allXYZ[2][xyz[2]])
    int istep;    // istep = 0, 1, 2 means edge in x, y, z direction
    Roots roots;  // intersections on that edge

    EdgeRoots(int xyz[3], int istep, Roots roots) : istep(istep), roots(roots) {
      for (int i = 0; i < 3; i++) this->xyz[i] = xyz[i];
    }
  };

  vector<EdgeRoots> allRoots;

  void calcEdgeRoots(POLY3 poly3) {
    allRoots.clear();
    int xyz[3], istep;
    Real rxyz[3];
    for (xyz[0] = 0; xyz[0] < allXYZ[0].size(); xyz[0]++) {
      rxyz[0] = allXYZ[0][xyz[0]];
      for (xyz[1] = 0; xyz[1] < allXYZ[1].size(); xyz[1]++) {
        rxyz[1] = allXYZ[1][xyz[1]];
        for (xyz[2] = 0; xyz[2] < allXYZ[2].size(); xyz[2]++) {
          rxyz[2] = allXYZ[2][xyz[2]];
          for (int istep = 0; istep < 3; istep++)
            if (xyz[istep] < allXYZ[istep].size() - 1) {
              Roots roots =
                  edgeRoots(poly3, rxyz, istep, allXYZ[istep][xyz[istep] + 1]);
              if (roots.size() > 0)
                allRoots.push_back(EdgeRoots(xyz, istep, roots));
            }
        }
      }
    }
  }
};

#define BOUNDARY_PAIRS
#ifdef BOUNDARY_PAIRS
std::vector<int> boundary_pairs(PTR<Object<Poly2>> f,
                                PTR<Object<PV2>> bottom_left,
                                PTR<Object<PV2>> top_right) {
  Encasement* enc = new Encasement(bottom_left, top_right, rand());
  enc->init();

  std::vector<PTR<Object<Poly2>>> curves = {f};
  enc->calculate(curves, false);

  int count = 0;
  vector<PTR<Face>>& faces = enc->faces;

#ifdef BLEEN
  for (int i = 0; i < faces.size(); ++i) {
    cout << "face " << i;
    Edge* start = faces[i]->e;
    Edge* e = start;
    do {
      cout << " " << i;
      cout << (e->intersection ? "i" : "e");
      if (e->twin->face == 0)
        cout << -1;
      else
        for (int j = 0; j < faces.size(); ++j)
          if (e->twin->face == faces[j]) cout << j;
    } while ((e = e->next) != start);
    cout << endl;
  }
#endif

  // get the boundary edge whose tail is bl
  Edge* edge0 = 0;
  for (int i = 0; i < faces.size(); ++i) {
    Edge* start = faces[i]->e;
    Edge* e = start;
    do {
      assert(++count < 1000000);
      if (!e->twin->face) break;
    } while ((e = e->next) != start);
    if (!e->twin->face) {
      edge0 = e;
      break;
    }
  }
  assert(edge0 != 0);
  while (edge0->tail != enc->bl) {
    assert(++count < 1000000);
    edge0 = edge0->next->twin->face ? edge0->next->twin->next : edge0->next;
  }

  // find all boundary edges that intersect f
  vector<PTR<Edge>> intersection_edges;
  Edge* e = edge0;
  do {
    assert(++count < 1000000);
    assert(!e->twin->face);
    if (e->intersection) intersection_edges.push_back(e);
    e = e->next->twin->face ? e->next->twin->next : e->next;
  } while (e != edge0);

  // boundary edge i connects to edge pairs[i]
  vector<int> pairs(intersection_edges.size());
  for (int i = 0; i < intersection_edges.size(); ++i) {
    Edge* start = intersection_edges[i];
    do {
      Edge* e = start;
      do {
        assert(++count < 1000000);
        e = e->next;
      } while (!e->intersection);
      assert(e != start);
      start = e->twin;
    } while (start->face);
    start = start->twin;
    pairs[i] = -1;
    for (int j = 0; j < intersection_edges.size(); ++j)
      if (intersection_edges[j] == start) pairs[i] = j;
    assert(pairs[i] != -1);
  }
  return pairs;
}

std::vector<std::vector<PTR<Object<PV2>>>> boundary_curves(PTR<Object<Poly2>> f,
                                                           Encasement* enc) {
  std::vector<std::vector<PTR<Object<PV2>>>> curves;

  for (const auto& face : enc->faces) {
    PTR<Edge> e = face->e;
    std::vector<PTR<Object<PV2>>> edge_intersections;

    do {
      if (e->intersection && e->poly == f)
        edge_intersections.push_back(e->intersection_point());
    } while ((e = e->next) != face->e);

    if (edge_intersections.size() > 0) {
      assert(edge_intersections.size() == 2);
      curves.push_back(edge_intersections);
    }
  }
  return curves;
}

// The bivariate polynomial resulting from substituting for one
// variable in a trivariate polynomial
class Sub1Poly3 : public Object<Poly2> {
  POLY3 poly3;
  int ixyz;
  Real xyz;
  DeclareCalculate(Poly2) {
    return poly3->get<N>().substitute1(ixyz, xyz->get<N>());
  }

 public:
  // ixyz = 0 or 1 or 2 means substituting for x or y or z
  Sub1Poly3(POLY3 poly3, int ixyz, Real xyz)
      : poly3(poly3), ixyz(ixyz), xyz(xyz) {}
};

class XYPoint : public Object<PV2> {
  Real x, y;
  DeclareCalculate(PV2) { return PV2<N>(x->get<N>(), y->get<N>()); }

 public:
  XYPoint(Real x, Real y) : x(x), y(y) {}
};

typedef GridRoots::EdgeRoots EdgeRoots;

int sideIndex(Square& sq, EdgeRoots& r) {
  int j = (sq.istep + 1) % 3;
  int k = (sq.istep + 2) % 3;
  assert(r.istep == j || r.istep == k);
  assert(r.xyz[sq.istep] == sq.xyz[sq.istep]);
  assert(r.xyz[j] == sq.xyz[j] || r.xyz[j] == sq.xyz[j] + 1);
  assert(r.xyz[k] == sq.xyz[k] || r.xyz[k] == sq.xyz[k] + 1);
  bool x = r.xyz[j] != sq.xyz[j];
  bool y = r.xyz[k] != sq.xyz[k];
  bool z = r.istep == k && !x;
  assert(x + y + z <= 1);
  return x + 2 * y + 3 * z;
}

bool lessThan(Square& sq, pair<EdgeRoots*, int>& a, pair<EdgeRoots*, int>& b) {
  int aSide = sideIndex(sq, *a.first);
  int bSide = sideIndex(sq, *b.first);
  if (aSide != bSide) return aSide < bSide;
  if (aSide < 2) return a.second < b.second;
  return b.second < a.second;
}

// sorted is a list of boundary intersections for square sq in
// counter-clockwise order starting at the bottom-left corner.
// inserts r into the sorted list.
void insertSorted(Square& sq, pair<EdgeRoots*, int>& r,
                  vector<pair<EdgeRoots*, int>>& sorted) {
  int i = sorted.size();
  sorted.push_back(r);
  while (i > 0 && lessThan(sq, r, sorted[i - 1])) {
    sorted[i] = sorted[i - 1];
    --i;
  }
  sorted[i] = r;
}

// returns pairs[] such that boundary intersection sorted[i] is
// connected to sorted[pairs[i]]
vector<int> boundary_pairs(PTR<Object<Poly3>> f3, GridRoots& gr,
                           const Square& sq) {
  PTR<Object<Poly2>> f2 =
      new Sub1Poly3(f3, sq.istep, gr.allXYZ[sq.istep][sq.xyz[sq.istep]]);
  int j = (sq.istep + 1) % 3;
  int k = (sq.istep + 2) % 3;
  PTR<Object<PV2>> minXY =
      new XYPoint(gr.allXYZ[j][sq.xyz[j]], gr.allXYZ[k][sq.xyz[k]]);
  PTR<Object<PV2>> maxXY =
      new XYPoint(gr.allXYZ[j][sq.xyz[j] + 1], gr.allXYZ[k][sq.xyz[k] + 1]);


  //vector<int> b = boundary_pairs(f2, minXY, maxXY);

  return boundary_pairs(f2, minXY, maxXY);
}

vector<vector<PTR<Object<PV2>>>> boundary_chains(PTR<Object<Poly2>> f,
                                                 PTR<Object<PV2>> bottom_left,
                                                 PTR<Object<PV2>> top_right) {
  Encasement* enc = new Encasement(bottom_left, top_right, rand());
  enc->init();

  std::vector<PTR<Object<Poly2>>> curves = {f};
  enc->calculate(curves, false);
  enc->refine(f, 0.1);

  vector<vector<PTR<Object<PV2>>>> b_curve = boundary_curves(f, enc);

  return b_curve;
}

// returns pairs[] such that boundary intersection sorted[i] is
// connected to sorted[pairs[i]]
vector<vector<vertex>> boundary_chains(PTR<Object<Poly3>> f3, GridRoots& gr,
                                       const Square& sq) {
  PTR<Object<Poly2>> f2 =
      new Sub1Poly3(f3, sq.istep, gr.allXYZ[sq.istep][sq.xyz[sq.istep]]);
  int j = (sq.istep + 1) % 3;
  int k = (sq.istep + 2) % 3;

  PTR<Object<PV2>> minXY =
      new XYPoint(gr.allXYZ[j][sq.xyz[j]], gr.allXYZ[k][sq.xyz[k]]);
  PTR<Object<PV2>> maxXY =
      new XYPoint(gr.allXYZ[j][sq.xyz[j] + 1], gr.allXYZ[k][sq.xyz[k] + 1]);
  vector<vector<PTR<Object<PV2>>>> b = boundary_chains(f2, minXY, maxXY);
  vector<vector<vertex>> chains;

  for (vector<PTR<Object<PV2>>> i : b) {
    vector<vertex> inner_chain;
    for (PTR<Object<PV2>> p : i) {
      PV2<Parameter> pp = p->getApprox();
      vertex v;
      v[sq.istep] = gr.allXYZ[sq.istep][sq.xyz[sq.istep]]->getApprox().x.mid();
      v[j] = pp.x.mid();
      v[k] = pp.y.mid();
      inner_chain.push_back(v);
    }
    chains.push_back(inner_chain);
  }
  return chains;
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

vector<PV3<Parameter>> boundary_intersections(PTR<Object<Poly3>> f3,
                                              PTR<Object<Poly3>> g3,
                                              GridRoots& gr, const Square& sq) {
  PTR<Object<Poly2>> f2 =
      new Sub1Poly3(f3, sq.istep, gr.allXYZ[sq.istep][sq.xyz[sq.istep]]);
  PTR<Object<Poly2>> g2 =
      new Sub1Poly3(g3, sq.istep, gr.allXYZ[sq.istep][sq.xyz[sq.istep]]);

  int j = (sq.istep + 1) % 3;
  int k = (sq.istep + 2) % 3;

  PTR<Object<PV2>> minXY =
      new XYPoint(gr.allXYZ[j][sq.xyz[j]], gr.allXYZ[k][sq.xyz[k]]);
  PTR<Object<PV2>> maxXY =
      new XYPoint(gr.allXYZ[j][sq.xyz[j] + 1], gr.allXYZ[k][sq.xyz[k] + 1]);

  vector<PTR<Object<PV2>>> b = boundary_intersections(f2, g2, minXY, maxXY);

  vector<PV3<Parameter>> intersections;

  std::transform(b.begin(), b.end(), std::back_inserter(intersections),
                 [&](const PTR<Object<PV2>>& p) {
                   PV2<Parameter> pp = p->getApprox();
                   PV3<Parameter> v;
                   v[sq.istep] =
                       gr.allXYZ[sq.istep][sq.xyz[sq.istep]]->getApprox().x;
                   v[j] = pp.x;
                   v[k] = pp.y;
                   return v;
                 });

  return intersections;
}

#endif

void makeSquare(const GridRoots::EdgeRoots& r, Square s[4]) {
  if (r.istep == 0) {
    s[0] = Square(r.xyz[0], r.xyz[1], r.xyz[2], 1);
    s[1] = Square(r.xyz[0], r.xyz[1], r.xyz[2], 2);
    s[2] = Square(r.xyz[0], r.xyz[1] - 1, r.xyz[2], 2);  // was a 1, now a 2
    s[3] = Square(r.xyz[0], r.xyz[1], r.xyz[2] - 1, 1);  // was a 2 now a 1
  } else if (r.istep == 1) {
    s[0] = Square(r.xyz[0], r.xyz[1], r.xyz[2], 0);
    s[1] = Square(r.xyz[0], r.xyz[1], r.xyz[2], 2);
    s[2] = Square(r.xyz[0] - 1, r.xyz[1], r.xyz[2], 2);  // was a 0 now a 2
    s[3] = Square(r.xyz[0], r.xyz[1], r.xyz[2] - 1, 0);  // was a 2 now a 0
  } else {
    s[0] = Square(r.xyz[0], r.xyz[1], r.xyz[2], 0);
    s[1] = Square(r.xyz[0], r.xyz[1], r.xyz[2], 1);
    s[2] = Square(r.xyz[0] - 1, r.xyz[1], r.xyz[2], 1);  // was a 0 now a 1
    s[3] = Square(r.xyz[0], r.xyz[1] - 1, r.xyz[2], 0);  // was a 1 now a 0
  }
}

bool neighbors(const GridRoots::EdgeRoots& a, const GridRoots::EdgeRoots& b) {
  Square sa[4];
  Square sb[4];
  makeSquare(a, sa);
  makeSquare(b, sb);

  return sa[0] == sb[0] || sa[0] == sb[1] || sa[0] == sb[2] || sa[0] == sb[3] ||
         sa[1] == sb[0] || sa[1] == sb[1] || sa[1] == sb[2] || sa[1] == sb[3] ||
         sa[2] == sb[0] || sa[2] == sb[1] || sa[2] == sb[2] || sa[2] == sb[3] ||
         sa[3] == sb[0] || sa[3] == sb[1] || sa[3] == sb[2] || sa[3] == sb[3];
}
/*class vertex {
public:
    float x, y, z;
};*/
class Cube {
 public:
  bool operator<(const Cube& rhs) const {
    if (this->x < rhs.x) {
      return true;
    } else if (this->x > rhs.x) {
      return false;
    } else {
      if (this->y < rhs.y) {
        return true;
      } else if (this->y > rhs.y) {
        return false;
      } else {
        return this->z < rhs.z;
      }
    }
  }

  int x, y, z;
  Cube(int x1, int y1, int z1) : x{x1}, y{y1}, z{z1} {}
};

// makes 6 squares of the Cube
void makeSquaresFromCube(Cube r, Square s[6]) {
  s[0] = Square(r.x, r.y, r.z, 0);
  s[1] = Square(r.x, r.y, r.z, 1);
  s[2] = Square(r.x, r.y, r.z, 2);
  s[3] = Square(r.x + 1, r.y, r.z, 0);
  s[4] = Square(r.x, r.y + 1, r.z, 1);
  s[5] = Square(r.x, r.y, r.z + 1, 2);
}

static void error_callback(int error, const char* description) {
  fprintf(stderr, "Error: %s\n", description);
}

float deg_x = 0.01;
float deg_y = 0.01;
float deg_z = 0.01;

float scale_factor = 1.0;

// float inc_deg = 0.000001;
static void key_callback(GLFWwindow* window, int key, int scancode, int action,
                         int mods);
static void rotate(float dx, float dy, float dz);
static void scale(float fScale);
static void move(float dx, float dy, float dz);
static void mouse(GLFWwindow* window, int button, int action, int mods);
static void motion(GLFWwindow* window, double x, double y);
static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

class DrawVerticesInCube {
 public:
  GridRoots grid;
  set<Cube> verticesInCube;
  // map<Cube, vector<GridRoots::EdgeRoots>> verticesInCube;
  map<Cube, GLuint> vboForCube;
  GLuint vertex_shader, fragment_shader, program, vertex_shader_sphere,
      fragment_shader_sphere, program_sphere, vertex_shader_line,
      fragment_shader_line, program_line, texture_vertex_shader,
      texture_fragment_shader, texture_program;
  GLuint texture_VBO, texture_EBO;
  GLuint texture;
  GLint mvp_location;
  GLint mit_location;
  GLint m_location;
  GLint v_location;
  GLint p_location;
  GLint mvp_location_lines;

  GLint color_location;
  GLint vpos_location;
  GLint vcol_location;
  GLuint texture_vpos_location;
  GLuint texture_texCoord_location;
  GLint vpos_location_sphere;
  GLint vcol_location_sphere;
  GLint vpos_location_line;
  GLint vcol_location_line;
  GLint light_location;
  GLFWwindow* window;

  glm::mat4 uMIT;
  glm::mat4 transform;
  glm::mat4 M;
  glm::mat4 V;
  glm::mat4 P;
  glm::mat4 uMITLines;
  glm::mat4 transformLines;
  glm::mat4 MLines;
  // glm::mat4 V;
  // glm::mat4 P;
  glm::vec3 light_world;
  glm::vec3 light;
  POLY3 poly3;

  std::vector<vector<vertex>> vertex_all_vect;
  pair<vector<vertex>, vector<vector<vertex>>> lineVertices;
  vector<vector<vertex>> fgVertices;

  std::vector<std::vector<double>> cubePlanes;
  std::vector<vertex> vertex_all;
  std::vector<vertex> vertex_all_intersect_points;
  map<POLY3, std::vector<vertex>> surfaceVerteces;
  // map<Poly1*,std::vector<vertex>> surfaceVerteces_poly1;

  // even is position
  // odd is color
  const char* vertex_shader_text =
      "uniform mat4 MVP;\n"
      "uniform mat4 V;\n"
      "uniform mat4 M;\n"
      "uniform mat4 MIT;\n"
      "attribute vec3 vCol;\n"
      "attribute vec3 vPos;\n"
      "varying float intensity;\n"
      "varying vec3 color;\n"
      "varying vec3 grad;\n"
      "varying vec3 position;\n"
      "varying vec3 camera_position;\n"

      "void main()\n"
      "{\n"
      "    gl_Position = MVP * vec4(vPos, 1.0);\n"
      "    color = (MVP*vec4(vCol*0.5+0.5,0.0)).xyz;//vec3(1.f,0.f,0.f);\n"
      "    grad = normalize((V*MIT*vec4(vCol,0.0)).xyz);//vec3(1.f,0.f,0.f);\n"
      "    if(grad.z < 0.0) { grad = -grad; }\n"
      "    position = gl_Position.xyz;\n"
      "    camera_position = (V * M * vec4(vPos, 1.0)).xyz;\n"
      "}\n";

  // this needs to have another variable that can store position
  // ligth can be passed in as a uniform variable
  // fragment == pixel
  // unit normal can be the unitized gradient
  // it needs to know where the light is and where the fragment/pixel
  const char* fragment_shader_text =
      "varying float intensity;\n"
      "uniform vec3 uRGB;\n"
      "uniform vec3 light;\n"  // light in view coordinates
      "vec3 myRGB = uRGB;\n"
      "varying vec3 color;\n"
      "varying vec3 grad;\n"
      "varying vec3 position;\n"
      "varying vec3 camera_position;\n"

      "float ambientCoefficient;\n"
      "vec3 color_of_light;\n"
      "vec3 diffuse;\n"

      "void main()\n"
      "{\n"

      "ambientCoefficient = 0.1;\n"
      "color_of_light = vec3(1, 1, 1);\n"
      "float materialShininess = 16.0;\n"
      "float specularStrength  = 0.5;\n"

      "vec3 surfaceToLight = -normalize(light);\n"
      "vec3 surfaceToCamera = -normalize(camera_position);\n"

      // ambient
      "vec3 ambient = ambientCoefficient * myRGB * color_of_light;\n"

      // diffuse
      "float intensity = dot(normalize(grad), surfaceToLight);\n"
      "if (intensity < 0.0) { intensity = -intensity; myRGB = 0.75 * "
      "myRGB; }\n"  // making back of surface a darker color
      // "if (intensity < 0.0) { intensity = -intensity; myRGB = 1.0 -
      // myRGB; }\n" //making back of surface a complementary color
      "diffuse = intensity*myRGB;\n"

      // specular
      "float specularCoefficient = 0.0;\n"
      "if (intensity > 0.0) { specularCoefficient = pow(max(0.0, "
      "dot(surfaceToCamera, reflect(-surfaceToLight, grad))), "
      "materialShininess); }\n"
      "vec3 specular = specularStrength * specularCoefficient * "
      "color_of_light;\n"

      // no light attenuation since directional lighting
      // attenuation
      //"float light_attenuation = 0.2;\n"
      //"float distanceToLight = length(light - position);\n"
      //"float attenuation = 1.0 / (1.0 + light_attenuation *
      // pow(distanceToLight, 2.0));\n"

      // linear color (color before gamma correction)
      "vec3 linearColor = ambient + diffuse + specular;\n"

      // final color (after gamma correction)
      //"vec3 gamma = vec3(1.0/2.2);\n"
      //"vec3 finalColor = pow(linearColor, gamma);\n"
      "vec3 finalColor = linearColor;\n"

      "gl_FragColor = vec4(finalColor,1.0);\n"
      // "gl_FragColor = vec4(diffuse,1.0);\n" //this is how we had it
      // before (but with the darker color instead of complementary)

      "}\n";

  const char* texture_vertex_shader_text =
      "attribute vec3 aPos;\n"
      "attribute vec3 aColor;\n"
      "attribute vec2 aTexCoord;\n"
      "varying vec3 ourColor;\n"
      "varying vec2 TexCoord;\n"
      "void main()\n"
      "{\n"
      "    gl_Position = vec4(aPos, 1.0);\n"
      "    ourColor = aColor;\n"
      "    TexCoord = aTexCoord;\n"
      "}\n";

  const char* texture_fragment_shader_text =
      "varying vec3 ourColor;\n"
      "varying vec2 TexCoord;\n"
      "uniform sampler2D ourTexture;\n"
      "void main()\n"
      "{\n"
      "    gl_FragColor = texture2D(ourTexture, TexCoord);\n"
      "}\n";

  const char* vertex_shader_text_sphere =
      "attribute vec3 vCol;\n"
      "attribute vec3 vPos;\n"
      "varying float intensity;\n"
      "varying vec3 color;\n"
      "varying vec3 grad;\n"
      "void main()\n"
      "{\n"
      "    gl_Position = vec4(vPos, 1.0);\n"
      "    color = (vec4(vCol*0.5+0.5,0.0)).xyz;//vec3(1.f,0.f,0.f);\n"
      "    grad = normalize((vec4(vCol,0.0)).xyz);//vec3(1.f,0.f,0.f);\n"
      "}\n";

  const char* fragment_shader_text_sphere =
      "void main()\n"
      "{\n"
      "vec3 RGB = vec3(0.8, 0, 0);\n"
      "gl_FragColor = vec4(RGB,1.0);\n"
      "}\n";

  const char* vertex_shader_text_line =
      "uniform mat4 MVP;\n"
      "attribute vec3 vPos;\n"
      "void main()\n"
      "{\n"
      "    gl_Position = MVP*vec4(vPos, 1.0);\n"
      "    gl_Position.z -= 0.05;\n"
      "}\n";

  const char* fragment_shader_text_line =
      "uniform vec3 colorLine;\n"
      "vec3 RGB = colorLine;\n"
      "void main()\n"
      "{\n"
      //"vec3 RGB = vec3(0.8, 0.5, 0.2);\n"
      "gl_FragColor = vec4(RGB,1.0);\n"
      "}\n";

  DrawVerticesInCube() = default;

  void setGrid(const GridRoots& grid) {
    this->grid = grid;
    //        grid.setMaxXYZ()
  }

  void setMaxXYZ(Real xyz[3]) { grid.setMaxXYZ(xyz); }
  void setMinXYZ(Real xyz[3]) { grid.setMinXYZ(xyz); }

  void setSteps(int nSteps[3]) { this->grid.setSteps(nSteps); }

  void calcEdgeRoots() { this->grid.calcEdgeRoots(poly3); }

  void updateVerticesInCube() {
    assert(0);
    verticesInCube = findCubes();
  }

  void setPoly3(POLY3 poly3) { this->poly3 = poly3; }

  vector<POLY3> poly3s;
  vector<POLY3> poly3s_original;
  void addPoly3(POLY3 poly3) { poly3s.push_back(poly3); }
  void addPoly3Original(POLY3 poly3) { poly3s_original.push_back(poly3); }

  pair<int, int> loopsNotClosed(
      pair<vector<pair<EdgeRoots*, int>>, vector<int>> squareInfoForCube[6]) {
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < squareInfoForCube[i].first.size(); j++) {
        if (squareInfoForCube[i].first[j].second != -1) {
          pair<int, int> p;
          p.first = i;
          p.second = j;

          return p;
        }
      }
    }
    pair<int, int> p;
    p.first = -1;
    p.second = -1;

    return p;
  }

  pair<int, int> findDuplicateEdgeRoot(
      pair<EdgeRoots*, int> EdgeIndexPair,
      pair<vector<pair<EdgeRoots*, int>>, vector<int>> squareInfoForCube[6]) {
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < squareInfoForCube[i].first.size(); j++) {
        if (EdgeIndexPair.first == squareInfoForCube[i].first[j].first &&
            EdgeIndexPair.second == squareInfoForCube[i].first[j].second) {
          pair<int, int> p;
          p.first = i;
          p.second = j;

          return p;
        }
      }
    }
    pair<int, int> p;
    p.first = -1;
    p.second = -1;

    return p;
  }
  void AddSurface(const GridRoots& grid, const POLY3& poly3) {
    

    vertex_all.clear();


    this->grid = grid;
    this->poly3 = poly3;
    verticesInCube.clear();

    this->grid.calcEdgeRoots(poly3);

    map<Square, pair<vector<pair<EdgeRoots*, int>>, vector<int>>> squareInfo =
        makeSquareInfoMap();

    set<Cube> cubeSet = findCubes();
    for (auto& cube : cubeSet) {
      vector<vector<pair<EdgeRoots*, int>>> verticesInCube;
      Cube c = cube;
      Square sq[6];
      makeSquaresFromCube(c, sq);
      pair<vector<pair<EdgeRoots*, int>>, vector<int>> squareInfoForCube[6];
      vector<pair<vector<pair<EdgeRoots*, int>>, vector<int>>> infoDebug;
      for (int i = 0; i < 6; ++i) {
        squareInfoForCube[i] = squareInfo[sq[i]];
        infoDebug.push_back(squareInfoForCube[i]);
      }
      pair<int, int> p;
      p.first = -1;
      p.second = -1;

      static int count;
      static int count2;
      count2++;

      while (loopsNotClosed(squareInfoForCube) != p) {
        vector<pair<EdgeRoots*, int>> verticesInLoop;

        // first vertice in loop
        pair<int, int> squareIndexVertexIndex =
            loopsNotClosed(squareInfoForCube);
        do {
          assert(count++ < 1000000);
          int squareIndex = squareIndexVertexIndex.first;
          int vertexIndex = squareIndexVertexIndex.second;
          verticesInLoop.push_back(
              squareInfoForCube[squareIndex].first[vertexIndex]);
          squareInfoForCube[squareIndex].first[vertexIndex].second = -1;

          vertexIndex = squareInfoForCube[squareIndex].second[vertexIndex];
          pair<EdgeRoots*, int> vertex =
              squareInfoForCube[squareIndex].first[vertexIndex];
          squareInfoForCube[squareIndex].first[vertexIndex].second = -1;
          squareIndexVertexIndex =
              findDuplicateEdgeRoot(vertex, squareInfoForCube);
          if (squareIndexVertexIndex == p) assert(vertex == verticesInLoop[0]);

        } while (squareIndexVertexIndex != p);

        verticesInCube.push_back(verticesInLoop);
      }

      // triangulation
      for (int a = 0; a < verticesInCube.size(); a++) {
        vector<pair<EdgeRoots*, int>> eroots = verticesInCube[a];

        int num = eroots.size();
        int nTri = num - 3 + 1;

        for (int i = 1; i < num - 1; i++) {
          int ms[] = {0, i, i + 1};
          for (int j = 0; j < 3; j++) {
            int m = ms[j];

            vertex color;
            vertex v;
            auto r = eroots[m % num].first;  // gives you EdgeRoots object

            Real vx = grid.allXYZ[0][r->xyz[0]];  // arrow instead of dot  //
                                                  // the edge is here
            Real vy = grid.allXYZ[1][r->xyz[1]];
            Real vz = grid.allXYZ[2][r->xyz[2]];
            v.x = grid.allXYZ[0][r->xyz[0]]->getApprox().x.mid();
            v.y = grid.allXYZ[1][r->xyz[1]]->getApprox().x.mid();
            v.z = grid.allXYZ[2][r->xyz[2]]->getApprox().x.mid();

            int k = eroots[m % num].second;  // int in the pair is k
            if (r->istep == 0) {
              vx = r->roots[k];
              v.x = r->roots[k]->getApprox().x.mid();
            } else if (r->istep == 1) {
              vy = r->roots[k];
              v.y = r->roots[k]->getApprox().x.mid();
            } else {
              vz = r->roots[k];
              v.z = r->roots[k]->getApprox().x.mid();
            }

            PTR<Object<PV3>> unitGrad = new UnitGradient(poly3, vx, vy, vz);
            PV3<Parameter> grad_unit = unitGrad->getApprox();

            color.x = grad_unit.x.mid();
            color.y = grad_unit.y.mid();
            color.z = grad_unit.z.mid();
            vertex_all.push_back(color);
            vertex_all.push_back(v);

            surfaceVerteces[poly3].push_back(v);
          }
        }
      }
    }
    vertex_all_vect.push_back(vertex_all);
  }

  explicit DrawVerticesInCube(const GridRoots& grid, const POLY3& poly3) {
    this->grid = grid;
    // iterating over all the roots
    // edge is shared by four cubes, unless its on a boundry
    verticesInCube = findCubes();

    initGLFW();

    bufferData();
    initShader();
    initShaderSphere();

    runGLFW();
  }

  void runGLFW() {
    double previousTime = glfwGetTime();
    int frameCount = 0;

    double fpsTime = glfwGetTime();

    double TARGET_FPS = 60.0;

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    glfwWindowHint(GLFW_REFRESH_RATE, 60);

    while (!glfwWindowShouldClose(window)) {
      float ratio;
      int width, height;

      int change;

      glfwGetFramebufferSize(window, &width, &height);
      ratio = width / (float)height;

      P[0][0] = 1 / ratio;



      glViewport(0, 0, width, height);

      if(width != previous_width || height != previous_height)
      {
        assignInverses();
        previous_width = width;
        previous_height = height;
      }

      glClearColor(0.0, 0.0, 0.0, 1.0);

      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      drawTexture();

      /**
     if(dregrid == true && regridding_on_screen == false)
      {   

       

        
          RText = new TextRenderer(width, height);
          RText->Load("nakula.ttf", 100);
          RText->RenderText("Regridding...", 5.0f, 5.0f, 1.0f, glm::vec3(0.0f, 1.0f, 0.0f));

          regridding_on_screen = true;
      }
      else
      {
        if(regridding_on_screen == true)
        { 
          key_callback(window, GLFW_KEY_V, glfwGetKeyScancode(GLFW_KEY_V), GLFW_PRESS, 0);
        }
      }
     **/
      glClear(GL_DEPTH_BUFFER_BIT);

      if (drawSurfaces) draw2();

      drawLine(lineVertices, fgVertices);

      glfwSwapBuffers(window);


      if(playback_object->keep_interpolating == true)
      {
        if(playback_object->isEmpty == true)
        {
          playback_object->keep_interpolating = false;
          main_interpolate();
        }
        else
        {
          main_interpolate();
        }
        
      }
      else
      {

        if(playback_object->interpolating == true)
        {
          main_interpolate();
        }
      }

      glfwPollEvents();

    

      // Measure speed
      double currentTime = glfwGetTime();
      frameCount++;
      // If a second has passed.
      if (currentTime - previousTime >= 1.0) {
        // Display the frame count here any way you want.
        double fps = (frameCount / (currentTime - previousTime));

        glfwSetWindowTitle(window, (std::string("3D Encasement ") +
                                    std::to_string(fps) + std::string(" fps"))
                                       .c_str());

        frameCount = 0;
        previousTime = currentTime;
      }

      while (glfwGetTime() < fpsTime + 1.0 / TARGET_FPS) {
        std::chrono::duration<double> sleepTime((fpsTime + 1.0 / TARGET_FPS) -
                                                glfwGetTime());
        std::this_thread::sleep_for(sleepTime);
      }
      fpsTime += 1.0 / TARGET_FPS;

      

    }

    glfwDestroyWindow(window);

    glfwTerminate();
  }

  void initGLFW() {
    if (!glfwInit()) exit(EXIT_FAILURE);
    glfwSetErrorCallback(error_callback);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

    window = glfwCreateWindow(640, 480, "3D Encasement", NULL, NULL);
    glfwSetMouseButtonCallback(window,
                               mouse);  // callback for mouse click inputs
    glfwSetCursorPosCallback(window, motion);  // callback for mouse movements
    glfwSetScrollCallback(window, scroll_callback);  // scole call back

    if (!window) {
      glfwTerminate();
      exit(EXIT_FAILURE);
    }

    glfwSetKeyCallback(window, key_callback);

    glfwMakeContextCurrent(window);
    if (glewInit() != GLEW_OK) {
      fprintf(stderr, "Failed to initialize GLEW\n");
      exit(EXIT_FAILURE);
    }

    glfwSwapInterval(1);

    float ratio;
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    previous_width = width;
    previous_height = height;

    ratio = width / (float)height;

    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDepthFunc(GL_LEQUAL);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }

  void initShader() {
    vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, &vertex_shader_text, NULL);
    glCompileShader(vertex_shader);

    GLint isCompiled = 0;
    glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE) {
      GLint maxLength = 0;
      glGetShaderiv(vertex_shader, GL_INFO_LOG_LENGTH, &maxLength);

      // The maxLength includes the NULL character
      std::vector<GLchar> errorLog(maxLength);
      glGetShaderInfoLog(vertex_shader, maxLength, &maxLength, &errorLog[0]);

      printf("<vertex shader compilation error>: %s\n", &errorLog[0]);

      // Provide the infolog in whatever manor you deem best.
      // Exit with failure.
      glDeleteShader(vertex_shader);  // Don't leak the shader.
      exit(1);
    }

    fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, &fragment_shader_text, NULL);
    glCompileShader(fragment_shader);

    isCompiled = 0;
    glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE) {
      GLint maxLength = 0;
      glGetShaderiv(fragment_shader, GL_INFO_LOG_LENGTH, &maxLength);

      // The maxLength includes the NULL character
      std::vector<GLchar> errorLog(maxLength);
      glGetShaderInfoLog(fragment_shader, maxLength, &maxLength, &errorLog[0]);

      printf("<fragment shader compilation error>: %s\n", &errorLog[0]);

      // Provide the infolog in whatever manor you deem best.
      // Exit with failure.
      glDeleteShader(fragment_shader);  // Don't leak the shader.
      exit(1);
    }

    program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);

    mvp_location = glGetUniformLocation(program, "MVP");
    mit_location = glGetUniformLocation(program, "MIT");
    color_location = glGetUniformLocation(program, "uRGB");
    m_location = glGetUniformLocation(program, "M");
    v_location = glGetUniformLocation(program, "V");
    p_location = glGetUniformLocation(program, "P");

    vpos_location = glGetAttribLocation(program, "vPos");
    vcol_location = glGetAttribLocation(program, "vCol");

    light_location = glGetUniformLocation(program, "light");

    transform = glm::mat4(1.0f);
    // cout<<glm::mat4(1.0f)<<endl;
    M = glm::mat4(1.0f);
    V = glm::mat4(1.0f);
    P = glm::ortho(-1.0f, 1.0f, -1.0f, 1.0f);

    light_world = glm::vec3(200, 300, -400);
    light = glm::vec3(
        V * glm::vec4(light_world, 0));  // view is not changing for now

    transform = P * V * M;
    glUseProgram(program);
    glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(transform));
    uMIT = glm::inverse(glm::transpose(M));
    glUniformMatrix4fv(mit_location, 1, GL_FALSE, glm::value_ptr(uMIT));
    glUniformMatrix4fv(v_location, 1, GL_FALSE, glm::value_ptr(V));
    glUniformMatrix4fv(m_location, 1, GL_FALSE, glm::value_ptr(M));

    glUniform3fv(light_location, 1, glm::value_ptr(light));

    P_Inverse = glm::inverse(P);
    V_Inverse = glm::inverse(V);
  }

  void initTexture() {
    // build and compile our shader zprogram
    // ------------------------------------
    texture_vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(texture_vertex_shader, 1, &texture_vertex_shader_text, NULL);
    glCompileShader(texture_vertex_shader);

    GLint isCompiled = 0;
    glGetShaderiv(texture_vertex_shader, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE) {
      GLint maxLength = 0;
      glGetShaderiv(texture_vertex_shader, GL_INFO_LOG_LENGTH, &maxLength);

      // The maxLength includes the NULL character
      std::vector<GLchar> errorLog(maxLength);
      glGetShaderInfoLog(texture_vertex_shader, maxLength, &maxLength,
                         &errorLog[0]);

      printf("<vertex shader compilation error>: %s\n", &errorLog[0]);

      // Provide the infolog in whatever manor you deem best.
      // Exit with failure.
      glDeleteShader(texture_vertex_shader);  // Don't leak the shader.
      exit(1);
    }

    texture_fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(texture_fragment_shader, 1, &texture_fragment_shader_text,
                   NULL);
    glCompileShader(texture_fragment_shader);

    isCompiled = 0;
    glGetShaderiv(texture_fragment_shader, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE) {
      GLint maxLength = 0;
      glGetShaderiv(texture_fragment_shader, GL_INFO_LOG_LENGTH, &maxLength);

      // The maxLength includes the NULL character
      std::vector<GLchar> errorLog(maxLength);
      glGetShaderInfoLog(texture_fragment_shader, maxLength, &maxLength,
                         &errorLog[0]);

      printf("<fragment shader compilation error>: %s\n", &errorLog[0]);

      // Provide the infolog in whatever manor you deem best.
      // Exit with failure.
      glDeleteShader(texture_fragment_shader);  // Don't leak the shader.
      exit(1);
    }

    texture_program = glCreateProgram();
    glAttachShader(texture_program, texture_vertex_shader);
    glAttachShader(texture_program, texture_fragment_shader);
    glLinkProgram(texture_program);

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------
    float vertices[] = {
        // positions          // colors           // texture coords
        1.0f,  1.0f,  0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f,  // top right
        1.0f,  -1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f,  // bottom right
        -1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,  // bottom left
        -1.0f, 1.0f,  0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f   // top left
    };
    unsigned int indices[] = {
        0, 1, 3,  // first triangle
        1, 2, 3   // second triangle
    };

    glGenBuffers(1, &texture_VBO);
    glGenBuffers(1, &texture_EBO);

    glBindBuffer(GL_ARRAY_BUFFER, texture_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, texture_EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices,
                 GL_STATIC_DRAW);

    texture_vpos_location = glGetAttribLocation(texture_program, "aPos");
    texture_texCoord_location =
        glGetAttribLocation(texture_program, "aTexCoord");

    // position attribute
    glVertexAttribPointer(texture_vpos_location, 3, GL_FLOAT, GL_FALSE,
                          8 * sizeof(float), nullptr);
    glEnableVertexAttribArray(texture_vpos_location);

    // texture coord attribute
    glVertexAttribPointer(texture_texCoord_location, 2, GL_FLOAT, GL_FALSE,
                          8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(texture_texCoord_location);

    // load and create a texture
    // -------------------------
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D,
                  texture);  // all upcoming GL_TEXTURE_2D operations now have
                             // effect on this texture object
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_REPEAT);  // set texture wrapping to GL_REPEAT (default
                                 // wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    // load image, create texture and generate mipmaps
    int width, height, nrChannels;
    // The FileSystem::getPath(...) is part of the GitHub repository so we can
    // find files on any IDE/platform; replace it with your own image path.
    stbi_set_flip_vertically_on_load(true);
    unsigned char* data =
        stbi_load("gradient.jpg", &width, &height, &nrChannels, 0);
    if (data) {
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB,
                   GL_UNSIGNED_BYTE, data);
      glGenerateMipmap(GL_TEXTURE_2D);
    } else {
      std::cout << "Failed to load texture" << std::endl;
    }
    stbi_image_free(data);
  }

  void initShaderSphere() {
    vertex_shader_sphere = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader_sphere, 1, &vertex_shader_text_sphere, NULL);
    glCompileShader(vertex_shader_sphere);

    fragment_shader_sphere = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader_sphere, 1, &fragment_shader_text_sphere,
                   NULL);
    glCompileShader(fragment_shader_sphere);

    program_sphere = glCreateProgram();
    glAttachShader(program_sphere, vertex_shader_sphere);
    glAttachShader(program_sphere, fragment_shader_sphere);
    glLinkProgram(program_sphere);

    vpos_location_sphere = glGetAttribLocation(program_sphere, "vPos");
    vcol_location_sphere = glGetAttribLocation(program_sphere, "vCol");
  }

  void initShaderLine() {
    vertex_shader_line = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader_line, 1, &vertex_shader_text_line, NULL);
    glCompileShader(vertex_shader_line);

    fragment_shader_line = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader_line, 1, &fragment_shader_text_line, NULL);
    glCompileShader(fragment_shader_line);

    program_line = glCreateProgram();
    glLineWidth(2.5);
    glAttachShader(program_line, vertex_shader_line);
    glAttachShader(program_line, fragment_shader_line);
    glLinkProgram(program_line);

    glUseProgram(program_line);
    mvp_location_lines = glGetUniformLocation(program_line, "MVP");

    vpos_location_line = glGetAttribLocation(program_line, "vPos");
    vcol_location_line = glGetUniformLocation(program_line, "colorLine");

    transformLines = glm::mat4(1.0f);
    MLines = glm::mat4(1.0f);
    transformLines = MLines * V * P;

    glUniformMatrix4fv(mvp_location_lines, 1, GL_FALSE,
                       glm::value_ptr(transformLines));
  }

  GLuint all_buffer;

  vector<vertex> all_vertices;

  GLuint abuffer;
  void bufferData() {
    assert(0);
#ifdef BLEEN
    // glGenVertexArrays(1, &abuffer);
    // glBindVertexArray(abuffer);
    Poly3 ppoly3 = poly3->getCurrentValue();
    for (auto& cube_vertices : verticesInCube) {
      auto& cube = cube_vertices.first;
      vector<GridRoots::EdgeRoots> eroots = cube_vertices.second;
      vertex color;
      vertex v;
      vector<vertex> vertices;
      // SPLITTING INTO TRIANGLES???
      int num = eroots.size();
      int nTri = num - 3 + 1;
      for (int i = 0; i < nTri; ++i)
        for (int m = i; m < i + 3; ++m)
        // for (auto r: eroots)
        {
          auto r = eroots[m];
          v.x = grid.allXYZ[0][r.xyz[0]]->getApprox().x.mid();
          v.y = grid.allXYZ[1][r.xyz[1]]->getApprox().x.mid();
          v.z = grid.allXYZ[2][r.xyz[2]]->getApprox().x.mid();

          for (int k = 0; k < r.roots.size(); ++k) {
            if (r.istep == 0) {
              v.x = r.roots[k]->getApprox().x.mid();
            } else if (r.istep == 1) {
              v.y = r.roots[k]->getApprox().x.mid();
            } else {
              v.z = r.roots[k]->getApprox().x.mid();
            }
          }
          // printf("%d,%d,%d: ",r.xyz[0],r.xyz[1],r.xyz[2]);
          // printf("%f,%f,%f\n",v.x,v.y,v.z);

          // smooth shader
          // calc gradient
          PV3 grad = ppoly3.gradient(v.x, v.y, v.z);
          PV3 grad_unit = grad.unit();
          color.x = grad_unit.x.mid();
          color.y = grad_unit.y.mid();
          color.z = grad_unit.z.mid();

          vertices.push_back(color);
          vertices.push_back(v);
        }

      GLuint buffer;
      glGenBuffers(1, &buffer);
      glBindBuffer(GL_ARRAY_BUFFER, buffer);
      glBufferData(GL_ARRAY_BUFFER, sizeof(vertex) * vertices.size(),
                   vertices.data(), GL_STATIC_DRAW);

      //            glEnableVertexAttribArray(vpos_location);
      //            glVertexAttribPointer(vpos_location, 3, GL_FLOAT, GL_FALSE,
      //            sizeof(vertex), nullptr);
      //
      vboForCube[cube] = buffer;
    }

    //        cout << sizeof(vertex) << endl;
    //        cout << all_vertices.size() << endl;
    //        glGenBuffers(1, &all_buffer);
    //        glBindBuffer(GL_ARRAY_BUFFER, all_buffer);
    //        glBufferData(GL_ARRAY_BUFFER, sizeof(vertex)*all_vertices.size(),
    //        all_vertices.data(), GL_STATIC_DRAW);

#endif
  }
  vector<GLuint> vbo_vect;
  // int iterationForDraw2;
  void draw2() {
    for (int i = 0; i < vertex_all_vect.size(); i++) {
      glUseProgram(program);

      glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(transform));
      glUniformMatrix4fv(v_location, 1, GL_FALSE, glm::value_ptr(V));
      glUniformMatrix4fv(m_location, 1, GL_FALSE, glm::value_ptr(M));
      uMIT = glm::inverse(glm::transpose(M));
      glUniformMatrix4fv(mit_location, 1, GL_FALSE, glm::value_ptr(uMIT));
      glUniform3fv(light_location, 1, glm::value_ptr(light));

      if (i == 0) {
        glUniform3f(color_location, 0.3f, 0.7f, 0.7f);
      }

      if (i == 2) {
        glUniform3f(color_location, 0.7f, 0.2f, 0.4f);
      }

      if (i == 1) {
        glUniform3f(color_location, 0.7f, 0.6f, 0.1f);
      }
      if (i== 3) {
        glUniform3f(color_location, 0.1f, 0.4f, 0.8f);
      }

      size_t offset = sizeof(vertex);
      glBindBuffer(GL_ARRAY_BUFFER, vbo_vect[i]);

      glEnableVertexAttribArray(vpos_location);
      glEnableVertexAttribArray(vpos_location_sphere);
      glVertexAttribPointer(vpos_location, 3, GL_FLOAT, GL_FALSE,
                            2 * sizeof(vertex), (GLvoid*)offset);
      glVertexAttribPointer(vpos_location_sphere, 3, GL_FLOAT, GL_FALSE,
                            2 * sizeof(vertex), (GLvoid*)offset);

      glBindBuffer(GL_ARRAY_BUFFER, vbo_vect[i]);
      glEnableVertexAttribArray(vcol_location);
      glEnableVertexAttribArray(vcol_location_sphere);
      glVertexAttribPointer(vcol_location, 3, GL_FLOAT, GL_FALSE,
                            2 * sizeof(vertex), nullptr);
      glVertexAttribPointer(vcol_location_sphere, 3, GL_FLOAT, GL_FALSE,
                            2 * sizeof(vertex), nullptr);
      glDrawArrays(GL_TRIANGLES, 0, vertex_all_vect[i].size() / 2);
      glUseProgram(0);
    }
  }

  void drawTexture() {
    // bind textures on corresponding texture units
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texture);

    // render container
    glUseProgram(texture_program);

    glBindBuffer(GL_ARRAY_BUFFER, texture_VBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, texture_EBO);

    texture_vpos_location = glGetAttribLocation(texture_program, "aPos");
    texture_texCoord_location =
        glGetAttribLocation(texture_program, "aTexCoord");

    // position attribute
    glVertexAttribPointer(texture_vpos_location, 3, GL_FLOAT, GL_FALSE,
                          8 * sizeof(float), nullptr);
    glEnableVertexAttribArray(texture_vpos_location);

    // texture coord attribute
    glVertexAttribPointer(texture_texCoord_location, 2, GL_FLOAT, GL_FALSE,
                          8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(texture_texCoord_location);

    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
  }

  void drawSphere() {
    glUseProgram(program_sphere);

    size_t offset = sizeof(vertex);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_vect[vertex_all_vect.size() - 1]);

    glEnableVertexAttribArray(vpos_location_sphere);
    glVertexAttribPointer(vpos_location_sphere, 3, GL_FLOAT, GL_FALSE,
                          2 * sizeof(vertex), (GLvoid*)offset);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_vect[vertex_all_vect.size() - 1]);
    glEnableVertexAttribArray(vcol_location_sphere);
    glVertexAttribPointer(vcol_location_sphere, 3, GL_FLOAT, GL_FALSE,
                          2 * sizeof(vertex), nullptr);
    glDrawArrays(GL_TRIANGLES, 0,
                 vertex_all_vect[vertex_all_vect.size() - 1].size() / 2);
    glUseProgram(0);
  }

  void drawLine(pair<vector<vertex>, vector<vector<vertex>>> lineVertices,
                vector<vector<vertex>> fgVertices) {
    if (drawBox) {
      GLuint vbo_cube;
      glGenBuffers(1, &vbo_cube);
      glBindBuffer(GL_ARRAY_BUFFER, vbo_cube);
      glBufferData(GL_ARRAY_BUFFER, sizeof(vertex) * lineVertices.first.size(),
                   &lineVertices.first[0], GL_STATIC_DRAW);

      glUseProgram(program_line);

      // this is for the box
      glUniform3f(vcol_location_line, 1.0f, 1.0f, 1.0f);
      glLineWidth(2.0);

      glBindBuffer(GL_ARRAY_BUFFER, vbo_cube);

      glUniformMatrix4fv(mvp_location_lines, 1, GL_FALSE,
                         glm::value_ptr(transformLines));
      glEnableVertexAttribArray(vpos_location_line);
      glVertexAttribPointer(vpos_location_line, 3, GL_FLOAT, GL_FALSE,
                            sizeof(vertex), nullptr);

      glDrawArrays(GL_LINES, 0, lineVertices.first.size());
      glUseProgram(0);
    }

    glLineWidth(4.0);

    if (drawLines) {
      for (int i = 0; i < lineVertices.second.size(); i++) {
        GLuint vbo_line;
        glGenBuffers(1, &vbo_line);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_line);
        glBufferData(GL_ARRAY_BUFFER,
                     sizeof(vertex) * lineVertices.second[i].size(),
                     &lineVertices.second[i][0], GL_STATIC_DRAW);

        glUseProgram(program_line);

        if (i % 4 == 0) {
          glUniform3f(vcol_location_line, 0.5f, 0.9f, 0.9f);
        }

        if (i % 4 == 2) {
          glUniform3f(vcol_location_line, 0.9f, 0.4f, 0.6f);
        }

        if (i % 4 == 1) {
          glUniform3f(vcol_location_line, 0.9f, 0.9f, 0.1f);
        }
        if (i % 4 == 3) {
          glUniform3f(vcol_location_line, 0.1f, 0.5f, 0.9f);
        }

        glBindBuffer(GL_ARRAY_BUFFER, vbo_line);

        glUniformMatrix4fv(mvp_location_lines, 1, GL_FALSE,
                           glm::value_ptr(transformLines));
        glEnableVertexAttribArray(vpos_location_line);
        glVertexAttribPointer(vpos_location_line, 3, GL_FLOAT, GL_FALSE,
                              sizeof(vertex), nullptr);

        glDrawArrays(GL_LINES, 0, lineVertices.second[i].size());
        glUseProgram(0);
      }
    }

    if (drawFG) {
      for (int i = 0; i < fgVertices.size(); i++) {
        double mn = curve_min[i];
        double mx = curve_max[i];

        // std::cout << "[" << i << "] min: " << mn << "  max: " << mx <<
        // std::endl;

        if (mx < mn) continue;

        int num_verts = static_cast<int>((mx - mn) * fgVertices[i].size());
        int vert_start = static_cast<int>(mn * (fgVertices[i].size() - 1));
        if (vert_start % 2 == 1) {
          if (vert_start == fgVertices.size() - 1)
            vert_start--;
          else
            vert_start++;
        }

        if (vert_start + num_verts > fgVertices[i].size()) {
          num_verts = (fgVertices[i].size() - vert_start);
        }

        // std::cout << num_verts << " out of " << fgVertices[i].size()
        //          << " vertices" << std::endl;
        // std::cout << "starting at index " << vert_start << std::endl;

        GLuint vbo_line;
        glGenBuffers(1, &vbo_line);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_line);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertex) * num_verts,
                     &fgVertices[i][vert_start], GL_STATIC_DRAW);

        glUseProgram(program_line);

        glUniform3f(vcol_location_line, 0.9f, 0.0f, 0.0f);

        glBindBuffer(GL_ARRAY_BUFFER, vbo_line);

        glUniformMatrix4fv(mvp_location_lines, 1, GL_FALSE,
                           glm::value_ptr(transformLines));
        glEnableVertexAttribArray(vpos_location_line);
        glVertexAttribPointer(vpos_location_line, 3, GL_FLOAT, GL_FALSE,
                              sizeof(vertex), nullptr);

        glDrawArrays(GL_LINES, 0, num_verts);
        glUseProgram(0);
      }
    }
  }

  void bufferData2() {
    for (int i = 0; i < vertex_all_vect.size(); i++) {
      GLuint buffer;
      glGenBuffers(1, &buffer);
      glBindBuffer(GL_ARRAY_BUFFER, buffer);
      int si = vertex_all_vect[i].size();
      glBufferData(GL_ARRAY_BUFFER, sizeof(vertex) * vertex_all_vect[i].size(),
                   &vertex_all_vect[i][0], GL_STATIC_DRAW);


      vbo_vect.push_back(buffer);  


      //vbo_vect[i] = buffer;
    }
  }

  // filling in squareInfo map
  map<Square, pair<vector<pair<EdgeRoots*, int>>, vector<int>>>
  makeSquareInfoMap() {

    map<Square, pair<vector<pair<EdgeRoots*, int>>, vector<int>>> squareInfo;
    for (int e = 0; e < grid.allRoots.size(); e++) {
      for (int i = 0; i < grid.allRoots[e].roots.size(); i++) {
        pair<EdgeRoots*, int> r;
        r.first = &grid.allRoots[e];
        r.second = i;
        Square se[4];
        makeSquare(grid.allRoots[e], se);

        
        for (int j = 0; j <= 3; j++) {
          int istep = se[j].istep;
          int jstep = (istep + 1) % 3;
          int kstep = (jstep + 1) % 3;
          if (se[j].xyz[0] == -1 || se[j].xyz[1] == -1 || se[j].xyz[2] == -1 ||
              se[j].xyz[istep] > grid.nSteps[istep] ||
              se[j].xyz[jstep] >= grid.nSteps[jstep] ||
              se[j].xyz[kstep] >= grid.nSteps[kstep])
            continue;
          insertSorted(se[j], r, squareInfo[se[j]].first);
        }
      }

    }

    int count = 0;
    for (auto& p : squareInfo) {
      if (count == 4) {
        vector<int> c = boundary_pairs(poly3, grid, p.first);
      }
      p.second.second = boundary_pairs(poly3, grid, p.first);
      count++;
    }
    return squareInfo;
  }

  set<Cube> findCubes() {
    set<Cube> cubeSet;

    // map from cube to edge root (edge with all the intersections on it)
    for (int i = 0; i < grid.allRoots.size(); i++) {
      // r is of type edge roots
      auto& r = grid.allRoots[i];

      int s1 = FXYZ(poly3, grid.allXYZ[0][r.xyz[0]], grid.allXYZ[1][r.xyz[1]],
                    grid.allXYZ[2][r.xyz[2]]);
      int step[3] = {0, 0, 0};
      step[r.istep] = 1;
      int s2 = FXYZ(poly3, grid.allXYZ[0][r.xyz[0] + step[0]],
                    grid.allXYZ[1][r.xyz[1] + step[1]],
                    grid.allXYZ[2][r.xyz[2] + step[2]]);

      if (s1 == s2) {
        if (verbose) cout << "same sign!" << endl;
      }

      // the outer three ifs are x,y,z
      if (r.istep == 0) {
        for (int k = 0; k < r.roots.size(); ++k) {
          // four cases, somewhere in the middle, along an edge, on a corner
          if (r.xyz[0] < grid.nSteps[0] && r.xyz[1] < grid.nSteps[1] &&
              r.xyz[2] < grid.nSteps[2]) {
            auto c1 = Cube(r.xyz[0], r.xyz[1], r.xyz[2]);
            cubeSet.insert(c1);
          }

          if (r.xyz[0] < grid.nSteps[0] && r.xyz[1] > 0 &&
              r.xyz[2] < grid.nSteps[2]) {
            auto c2 = Cube(r.xyz[0], r.xyz[1] - 1, r.xyz[2]);
            cubeSet.insert(c2);
          }

          if (r.xyz[0] < grid.nSteps[0] && r.xyz[1] < grid.nSteps[1] &&
              r.xyz[2] > 0) {
            auto c3 = Cube(r.xyz[0], r.xyz[1], r.xyz[2] - 1);
            cubeSet.insert(c3);
          }

          if (r.xyz[0] < grid.nSteps[0] && r.xyz[1] > 0 && r.xyz[2] > 0) {
            auto c4 = Cube(r.xyz[0], r.xyz[1] - 1, r.xyz[2] - 1);
            cubeSet.insert(c4);
          }
        }
      } else if (r.istep == 1) {
        for (int k = 0; k < r.roots.size(); ++k) {
          if (r.xyz[0] < grid.nSteps[0] && r.xyz[1] < grid.nSteps[1] &&
              r.xyz[2] < grid.nSteps[2]) {
            auto c1 = Cube(r.xyz[0], r.xyz[1], r.xyz[2]);
            cubeSet.insert(c1);
          }

          if (r.xyz[0] > 0 && r.xyz[1] < grid.nSteps[1] &&
              r.xyz[2] < grid.nSteps[2]) {
            auto c2 = Cube(r.xyz[0] - 1, r.xyz[1], r.xyz[2]);
            cubeSet.insert(c2);
          }

          if (r.xyz[0] < grid.nSteps[0] && r.xyz[1] < grid.nSteps[1] &&
              r.xyz[2] > 0) {
            auto c3 = Cube(r.xyz[0], r.xyz[1], r.xyz[2] - 1);
            cubeSet.insert(c3);
          }
          if (r.xyz[0] > 0 && r.xyz[1] < grid.nSteps[1] && r.xyz[2] > 0) {
            auto c4 = Cube(r.xyz[0] - 1, r.xyz[1], r.xyz[2] - 1);
            cubeSet.insert(c4);
          }
        }
      } else {
        for (int k = 0; k < r.roots.size(); ++k) {
          if (r.xyz[0] < grid.nSteps[0] && r.xyz[1] < grid.nSteps[1] &&
              r.xyz[2] < grid.nSteps[2]) {
            auto c1 = Cube(r.xyz[0], r.xyz[1], r.xyz[2]);
            cubeSet.insert(c1);
          }

          if (r.xyz[0] > 0 && r.xyz[1] < grid.nSteps[1] &&
              r.xyz[2] < grid.nSteps[2]) {
            auto c2 = Cube(r.xyz[0] - 1, r.xyz[1], r.xyz[2]);
            cubeSet.insert(c2);
          }

          if (r.xyz[0] < grid.nSteps[0] && r.xyz[1] > 0 &&
              r.xyz[2] < grid.nSteps[2]) {
            auto c3 = Cube(r.xyz[0], r.xyz[1] - 1, r.xyz[2]);
            cubeSet.insert(c3);
          }
          if (r.xyz[0] > 0 && r.xyz[1] > 0 && r.xyz[2] < grid.nSteps[2]) {
            auto c4 = Cube(r.xyz[0] - 1, r.xyz[1] - 1, r.xyz[2]);
            cubeSet.insert(c4);
          }
        }
      }
    }
    return cubeSet;
  };
};

auto drawer = DrawVerticesInCube();
auto draw_sphere = DrawVerticesInCube();
// Set number of steps to 10 for each dimension.
int nSteps[3] = {15, 15, 15};
// Set lower bounds to (-0.5, -0.5, -0.5)
Real lo = new Object<Scalar>(Parameter::constant(-1.0));
Real los[3] = {lo, lo, lo};

// Set upper bounds to (0.5, 0.5, 0.5)
Real hi = new Object<Scalar>(Parameter::constant(1.0));
Real his[3] = {hi, hi, hi};

double posx;
double posy;
double posz;
int button_mode = -1;



static void key_callback(GLFWwindow* window, int key, int scancode, int action,
                         int mods) {
  float ratio;
  int width, height;
  glfwGetFramebufferSize(window, &width, &height);
  ratio = width / (float)height;

  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
    glfwSetWindowShouldClose(window, GL_TRUE);
  } else if (key == GLFW_KEY_UP &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    if (mods & GLFW_MOD_SHIFT) {
      curve_max[current_curve] = std::min(1.0, curve_max[current_curve] + 0.01);
      std::cout << "max[" << current_curve << "]: " << curve_max[current_curve]
                << std::endl;
    } else {
      curve_min[current_curve] = std::min(1.0, curve_min[current_curve] + 0.01);
      std::cout << "min[" << current_curve << "]: " << curve_min[current_curve]
                << std::endl;
    }
  } else if (key == GLFW_KEY_DOWN &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    if (mods & GLFW_MOD_SHIFT) {
      curve_max[current_curve] = std::max(0.0, curve_max[current_curve] - 0.01);
      std::cout << "max[" << current_curve << "]: " << curve_max[current_curve]
                << std::endl;
    } else {
      curve_min[current_curve] = std::max(0.0, curve_min[current_curve] - 0.01);
      std::cout << "min[" << current_curve << "]: " << curve_min[current_curve]
                << std::endl;
    }
  } else if (key == GLFW_KEY_LEFT &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    current_curve = (current_curve - 1 + curve_min.size()) % curve_min.size();

    std::cout << "current curve: " << current_curve << std::endl;
  } else if (key == GLFW_KEY_RIGHT &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    current_curve = (current_curve + 1) % curve_min.size();
    std::cout << "current curve: " << current_curve << std::endl;
  } else if (key == GLFW_KEY_O &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    rotate(0, 0, deg_z);
  } else if (key == GLFW_KEY_P &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    rotate(0, 0, -deg_z);
  } else if (key == GLFW_KEY_A &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    move(-deg_x, 0, 0);
  } else if (key == GLFW_KEY_D &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    move(deg_x, 0, 0);
  } else if (key == GLFW_KEY_S &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    move(0, -deg_y, 0);
  } else if (key == GLFW_KEY_W &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    move(0, deg_y, 0);
  }

  else if (key == GLFW_KEY_Q &&
           (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    move((float)0, 0, -deg_z);
  } else if (key == GLFW_KEY_E &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    move((float)0, 0, deg_z);
  }

  else if (key == GLFW_KEY_J &&
           (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    drawSurfaces = !drawSurfaces;
  } else if (key == GLFW_KEY_K &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    drawLines = !drawLines;
  } else if (key == GLFW_KEY_L &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    drawBox = !drawBox;
  } else if (key == GLFW_KEY_H &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {
    drawFG = !drawFG;
  } else if (key == GLFW_KEY_M &&
             (action == GLFW_REPEAT || action == GLFW_PRESS)) {

      if(playback_object->isEmpty == true && playback_object->interpolating == false)
      {
        std::cout << "Creating new playback object. \n\n";
        playback_object = new Playback(interpolate_file);
      }
      main_interpolate();
        
  } else if(key == GLFW_KEY_N && 
            (action == GLFW_REPEAT || action == GLFW_PRESS)){
    

    playback_object->recordCurrentPose("poses2.txt", drawer.M);


  } else if(key == GLFW_KEY_V && 
            (action == GLFW_REPEAT || action == GLFW_PRESS)){

    vregrid();

  } else if(key == GLFW_KEY_G && (action == GLFW_REPEAT || action == GLFW_PRESS)){


    restorePoly();


  } else if(key == GLFW_KEY_B && (action == GLFW_REPEAT || action == GLFW_PRESS)){
    
    
      if(playback_object->isEmpty == true && playback_object->interpolating == false)
      {
        std::cout << "Creating new playback object. \n\n";
        playback_object = new Playback(interpolate_file);
      }
    
    playback_object->keep_interpolating = true;
    //main_interpolate();
    
  } else if (key == GLFW_KEY_U && action == GLFW_PRESS) {
    cout << "You pressed u!!" << endl;

    std::clock_t start;
    double duration;
    start = std::clock();

    drawer.vbo_vect.clear();
   // drawer.vertex_all.clear();
    drawer.vertex_all_vect.clear();

    glm::mat4 PV_Inv = glm::inverse(drawer.P * drawer.V);
    glm::mat4 PVM_Inv = glm::inverse(drawer.transform);
    drawer.M = PV_Inv;
    drawer.transform = drawer.P * drawer.V * drawer.M;

    drawer.MLines = PV_Inv;
    drawer.transformLines = drawer.P * drawer.V * drawer.MLines;

    for (int i = 0; i < drawer.poly3s.size(); i++) {
      drawer.poly3s[i] = new TransformedPoly3(drawer.poly3s[i], PVM_Inv);
      if (verbose) cout << "transformed poly " << i << endl;
    }

    for (int i = 0; i < drawer.poly3s.size(); i++) {
      drawer.AddSurface(drawer.grid, drawer.poly3s[i]);
    }
    drawer.bufferData2();

    // Get new cell vertices, put them into lineVertices
    // set vmin and vmax
    // glm::vec4 vmin_gl = PVM_Inv * glm::vec4(-1, -1, -1, 1);
    // glm::vec4 vmax_gl = PVM_Inv * glm::vec4(1, 1, 1, 1);

    // cout << "TEST w: " << vmin_gl.w << " " << vmax_gl.w << endl;

    // vmin.x = vmin_gl.x;
    // vmin.y = vmin_gl.y;
    // vmin.z = vmin_gl.z;
    // vmax.x = vmax_gl.x;
    // vmax.y = vmax_gl.y;
    // vmax.z = vmax_gl.z;

    // cout << "vmin: " << vmin.x << " " << vmin.y << " " << vmin.z << endl;
    // cout << "vmax: " << vmax.x << " " << vmax.y << " " << vmax.z << endl;

    // setCubeFaces();

    drawer.lineVertices.second.clear();

    for (int i = 0; i < drawer.poly3s.size(); i++) {
      displayChains(drawer.poly3s[i]);
    }

    current_curve = 0;
    curve_min.clear();
    curve_max.clear();
    drawer.fgVertices.clear();

    for (int i = 0; i < drawer.poly3s.size(); i++) {
      for (int j = i + 1; j < drawer.poly3s.size(); j++) {
        displayFGCurves(drawer.poly3s[i], drawer.poly3s[j]);
      }
    }

    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    cout << "regrid time " << duration << endl;
  }
}

static void rotate(float dx, float dy, float dz) {
  
  glm::mat4 R(1.f);
  R = glm::rotate(R, -dx, glm::vec3(0.0f, 1.0f, 0.0f));
  R = glm::rotate(R, -dy, glm::vec3(1.0f, 0.0f, 0.0f));
  R = glm::rotate(R, -dz, glm::vec3(0.0f, 0.0f, 1.0f));
  drawer.M = R * drawer.M;
  drawer.MLines = R * drawer.MLines;

  updateModel();
}

static void scale(float fScale) {

  glm::mat4 S(1.f);
  S = glm::scale(S, glm::vec3(fScale, fScale, fScale));
  drawer.M = S * drawer.M;
  drawer.MLines = S * drawer.MLines;
 
  updateModel();
                     
}

static void move(float dx, float dy, float dz) {
  glm::mat4 T(1.f);
  T = glm::translate(T, glm::vec3(dx, dy, dz));
  drawer.M = T * drawer.M;
  drawer.MLines = T * drawer.MLines;
  
  updateModel();
}

static void scroll_callback(GLFWwindow* window, double xoffset,
                            double yoffset) {
  if (yoffset < 0 && (button_mode == -1)) {
    scale(1.1);
  }

  else if (button_mode == -1) {
    scale(0.9);
  }

  else if (yoffset < 0 && (button_mode == 2)) {
    // double dz=-2*(yoffset - posz)/height;
    rotate(0, 0, -deg_z);
    // posz=yoffset;
  } else if (yoffset > 0 && (button_mode == 2)) {
    // double dz=-2*(yoffset - posz)/height;
    rotate(0, 0, deg_z);
    // posz=yoffset;
  }

  else if (yoffset < 0 && (button_mode == 1)) {
    // double dz=-2*(yoffset - posz)/height;
    move((float)0, 0, -deg_z);
    // posz=yoffset;
  } else if (yoffset > 0 && (button_mode == 1)) {
    // double dz=-2*(yoffset - posz)/height;
    move((float)0, 0, deg_z);
    // posz=yoffset;
  }
}
static void mouse(GLFWwindow* window, int button, int action, int mods) {
  if (action == GLFW_PRESS) {
    glfwGetCursorPos(window, &posx, &posy);
    if (button == GLFW_MOUSE_BUTTON_LEFT)
      button_mode = 1;
    else if (button == GLFW_MOUSE_BUTTON_RIGHT)
      button_mode = 2;
    // else if(button == GLFW_MOUSE_BUTTON_LEFT && )

  } else if (action == GLFW_RELEASE) {
    button_mode = -1;
  }
}

static void motion(GLFWwindow* window, double x, double y) {
  int width, height;
  glfwGetFramebufferSize(drawer.window, &width, &height);
  if (button_mode == 1) {
    float dx = 2 * (x - posx) / width;
    float dy = -2 * (y - posy) / height;
    move(dx, dy, (float)0);
    posx = x;
    posy = y;
  } else if (button_mode == 2) {
    double dx = -2 * (x - posx) / width;
    double dy = -2 * (y - posy) / height;
    // double dz=-2*(y-posy)/height;
    rotate(dx, dy, 0);
    posx = x;
    posy = y;
  }
}

void ReadPolyText(const char* filename, Poly3<Parameter>* ppoly3) {
  std::ifstream fin;
  fin.open(filename);
  string str;

  getline(fin, str);
  getline(fin, str);
  getline(fin, str);
  getline(fin, str);
  getline(fin, str);
  getline(fin, str);
  int num = atoi(str.c_str());
  getline(fin, str);

  const char* d = " ";
  for (int i = 0; i < num; ++i) {
    getline(fin, str);
    char* p;
    p = strtok((char*)str.c_str(), d);
    double tmp[4];
    int m = 0;
    while (p) {
      tmp[m++] = atof(p);
      p = strtok(NULL, d);
    }
    ppoly3->add(tmp[0], tmp[1], tmp[2], Parameter::input(tmp[3]));
  }
  return;
}

void ReadBondingBox(const char* filename) {
  std::ifstream fin;
  fin.open(filename);
  string str;

  getline(fin, str);
  getline(fin, str);

  int num = atoi(str.c_str());
  const char* d = " ";

  getline(fin, str);

  for (int i = 0; i < num; ++i) {
    getline(fin, str);
    char* p;
    p = strtok((char*)str.c_str(), d);
    double tmp[3];
    int m = 0;
    while (p) {
      tmp[m++] = atof(p);
      p = strtok(NULL, d);
    }
    // vertex v;
    if (i == 0) {
      vmin.x = tmp[0];
      vmin.y = tmp[1];
      vmin.z = tmp[2];
    }
    if (i == 1) {
      vmax.x = tmp[0];
      vmax.y = tmp[1];
      vmax.z = tmp[2];
    }
  }

  cout << "vmin: " << vmin.x << " " << vmin.y << " " << vmin.z << endl;
  cout << "vmax: " << vmax.x << " " << vmax.y << " " << vmax.z << endl;

  los[0] = new Object<Scalar>(Parameter::input(vmin.x));
  los[1] = new Object<Scalar>(Parameter::input(vmin.y));
  los[2] = new Object<Scalar>(Parameter::input(vmin.z));
  his[0] = new Object<Scalar>(Parameter::input(vmax.x));
  his[1] = new Object<Scalar>(Parameter::input(vmax.y));
  his[2] = new Object<Scalar>(Parameter::input(vmax.z));

  setCubeFaces();
}

void setCubeFaces() {
  //     v3 + + + + + + + + + + + v2
  //      +\                      +\
  //      + \                     + \
  //      +  \                    +  \
  //      +   \                   +   \
  //      +   v7 + + + + + + + + +++ + v6
  //      +    +                  +    +
  //      +    +                  +    +
  //      +    +                  +    +   ^
  //      +    +                  +    +   |
  //      +    +                  +    +   |
  //      +    +                  +    +   Z
  //     v0 + +++ + + + + + + + + v1   +   |
  // ^    \   +                    \   +   |
  //  \    \  +                     \  +   v
  //   Y    \ +                      \ +
  //    \    \+                       \+
  //     v    v4 + + + + + + + + + + + v5
  //
  //                 <---X--->
  // v0 = vmin
  // v6 = vmax

  // 6 Outward Facing Faces:
  //
  //  {v0, v4, v7, v3} Left
  //  {v4, v5, v6, v7} Front
  //  {v5, v1, v2, v6} Right
  //  {v1, v0, v3, v2} Back
  //  {v3, v7, v6, v2} Top
  //  {v0, v1, v5, v4} Bottom

  vertex v0(vmin.x, vmin.y, vmin.z);
  vertex v1(vmax.x, vmin.y, vmin.z);
  vertex v2(vmax.x, vmax.y, vmin.z);
  vertex v3(vmin.x, vmax.y, vmin.z);

  vertex v4(vmin.x, vmin.y, vmax.z);
  vertex v5(vmax.x, vmin.y, vmax.z);
  vertex v6(vmax.x, vmax.y, vmax.z);
  vertex v7(vmin.x, vmax.y, vmax.z);

  vector<vertex> cubeLines;

  cubeLines.insert(cubeLines.end(),
                   {v0, v1, v1, v2, v2, v3, v3, v0});  // Back Face
  cubeLines.insert(cubeLines.end(),
                   {v4, v5, v5, v6, v6, v7, v7, v4});  // Front Face
  cubeLines.insert(cubeLines.end(),
                   {v0, v4, v1, v5, v2, v6, v3, v7});  // Connecting Edges

  drawer.lineVertices.first = cubeLines;

  std::vector<std::vector<double>> cubePlanes;

  cubePlanes.push_back({1, 0, 0, vmin.x});
  cubePlanes.push_back({-1, 0, 0, vmax.x});
  cubePlanes.push_back({0, 1, 0, vmin.y});
  cubePlanes.push_back({0, -1, 0, vmax.y});
  cubePlanes.push_back({0, 0, 1, vmin.z});
  cubePlanes.push_back({0, 0, -1, vmax.z});

  drawer.cubePlanes = cubePlanes;

  return;
}

void displayChains(PTR<Object<Poly3>> pol3) {
  auto drawerforCube = DrawVerticesInCube();
  // GridRoots grid;

  // TODO read from file???
  int nSteps[3] = {1, 1, 1};

  drawerforCube.grid.setMinXYZ(los);
  drawerforCube.grid.setMaxXYZ(his);
  drawerforCube.grid.setSteps(nSteps);

  Square sq[6];

  Cube cube = Cube(0, 0, 0);
  makeSquaresFromCube(cube, sq);
  vector<vertex> surfaceLines;
  vector<vector<vertex>> surfaceLinesVect;

  for (int i = 0; i < 6; i++) {
    vector<vector<vertex>> chainsOnSquare =
        boundary_chains(pol3, drawerforCube.grid, sq[i]);

    for (vector<vertex> chain : chainsOnSquare) {
      for (vertex v : chain) {
        surfaceLines.push_back(v);
      }
    }
  }
  surfaceLinesVect.push_back(surfaceLines);
  drawer.lineVertices.second.insert(drawer.lineVertices.second.end(),
                                    surfaceLinesVect.begin(),
                                    surfaceLinesVect.end());
}

void displayFGCurves(PTR<Object<Poly3>> f, PTR<Object<Poly3>> g) {
  auto drawerforCube = DrawVerticesInCube();

  int nSteps[3] = {1, 1, 1};

  drawerforCube.grid.setMinXYZ(los);
  drawerforCube.grid.setMaxXYZ(his);
  drawerforCube.grid.setSteps(nSteps);

  Square sq[6];
  Cube cube = Cube(0, 0, 0);
  makeSquaresFromCube(cube, sq);

  vector<PV3<Parameter>> all_intersections;

  for (int i = 0; i < 6; i++) {
    vector<PV3<Parameter>> intersections =
        boundary_intersections(f, g, drawerforCube.grid, sq[i]);

    all_intersections.insert(all_intersections.begin(), intersections.begin(),
                             intersections.end());
  }

  std::cout << "Total boundary intersections of f and g: "
            << all_intersections.size() << std::endl;
  /*
    if (all_intersections.size() == 2) {
      surfaceLines.push_back(vertex(all_intersections[0].x.mid(),
                                    all_intersections[0].y.mid(),
                                    all_intersections[0].z.mid()));
      surfaceLines.push_back(vertex(all_intersections[1].x.mid(),
                                    all_intersections[1].y.mid(),
                                    all_intersections[1].z.mid()));
      drawer.lineVertices.second.insert(drawer.lineVertices.second.end(),
                                        surfaceLines);
      std::cout << "second lines size: " << drawer.lineVertices.second.size()
                << std::endl;
      return;
    }
  */
  Poly3<Parameter> fp = f->getApprox(1.0);
  Poly3<Parameter> gp = g->getApprox(1.0);

  PV3<Parameter> box(Parameter::interval(los[0]->getApprox(1.0).x.lb(),
                                         his[0]->getApprox(1.0).x.ub()),
                     Parameter::interval(los[1]->getApprox(1.0).x.lb(),
                                         his[1]->getApprox(1.0).x.ub()),
                     Parameter::interval(los[2]->getApprox(1.0).x.lb(),
                                         his[2]->getApprox(1.0).x.ub()));

  std::map<std::pair<int, int>, vector<PV3<Parameter>>> curves =
      acp::encasement_utils::FGCurves(fp, gp, box, all_intersections);

  std::cout << "fg curves size: " << curves.size() << std::endl;

  for (auto& [ind, segment] : curves) {
    vector<vertex> surfaceLines;
    std::transform(segment.begin(), segment.end(),
                   std::back_inserter(surfaceLines),
                   [](const PV3<Parameter>& p) {
                     return vertex(p.x.mid(), p.y.mid(), p.z.mid());
                   });
    drawer.fgVertices.insert(drawer.fgVertices.end(), surfaceLines);
    curve_min.push_back(0);
    curve_max.push_back(1);
  }
}


void vregrid()
{
    cout << "You pressed v!!" << endl;
    std::clock_t start;

    double duration;
    start = std::clock();

    drawer.vbo_vect.clear();

    //drawer.vertex_all.clear();
    drawer.vertex_all_vect.clear();

    //T-INV
    //glm::mat4 PVM_INV = glm::inverse(drawer.transform);
    glm::mat4 PVM_INV =  glm::transpose(drawer.uMIT)* V_Inverse *P_Inverse;

    glm::vec4 triangle_vertex_data = glm::vec4(0, 0, 0, 1.0f);
    glm::vec4 line_vertex_data = glm::vec4(0, 0, 0, 1.0f);
    glm::vec4 FG_vertex_data = glm::vec4(0, 0, 0, 1.0f);
    //Change f to the original f AND apply T-INV
    for (int i = 0; i < drawer.poly3s.size(); i++)
    {
      drawer.poly3s[i] = new TransformedPoly3(drawer.poly3s_original[i], PVM_INV);
      if (verbose) cout << "transformed poly " << i << endl;
    }


    //Set triangles vertices

    
    for (int i = 0; i < drawer.poly3s.size(); i++) 
    {
      drawer.AddSurface(drawer.grid, drawer.poly3s[i]);
    }


    //Apply T-INV to all triangle vertices

    for(int i = 0; i < drawer.vertex_all_vect.size(); i++)
    {
      for(int j = 0; j < drawer.vertex_all_vect[i].size(); j++)
      {
        vertex &triangle_vertex = drawer.vertex_all_vect[i][j];

        triangle_vertex_data[0] = triangle_vertex.x;
        triangle_vertex_data[1] = triangle_vertex.y;
        triangle_vertex_data[2] = triangle_vertex.z;
        triangle_vertex_data[3] = 1.0f;

        triangle_vertex_data = PVM_INV * triangle_vertex_data;

        drawer.vertex_all_vect[i][j].x = triangle_vertex_data[0];
        drawer.vertex_all_vect[i][j].y = triangle_vertex_data[1];
        drawer.vertex_all_vect[i][j].z = triangle_vertex_data[2];
      }
    }

    //Update buffers 
    drawer.bufferData2();

    
    drawer.lineVertices.second.clear();

    for (int i = 0; i < drawer.poly3s.size(); i++) 
    {
      displayChains(drawer.poly3s[i]);
    }

    for(int i = 0; i < drawer.lineVertices.second.size(); i++)
    {
      for(int j = 0; j < drawer.lineVertices.second[i].size(); j++)
      {
        vertex line_vertex = drawer.lineVertices.second[i][j];

        line_vertex_data[0] = line_vertex.x;
        line_vertex_data[1]= line_vertex.y;
        line_vertex_data[2] = line_vertex.z;
        line_vertex_data[3] = 1.0f;


        line_vertex_data = PVM_INV * line_vertex_data;

        drawer.lineVertices.second[i][j].x = line_vertex_data[0];
        drawer.lineVertices.second[i][j].y = line_vertex_data[1];
        drawer.lineVertices.second[i][j].z = line_vertex_data[2];
      }
    }

    //clear gridbox coordinates and set them back to {+-1, +-1, +-1}
    drawer.lineVertices.first.clear();
    setCubeFaces();

    //apply T-INV to gridbox
    for(int i = 0; i < drawer.lineVertices.first.size(); i++)
    {
      vertex line_vertex = drawer.lineVertices.first[i];
      
      line_vertex_data[0] = line_vertex.x;
      line_vertex_data[1]= line_vertex.y;
      line_vertex_data[2] = line_vertex.z;
      line_vertex_data[3] = 1.0f;

      line_vertex_data = PVM_INV * line_vertex_data;

      drawer.lineVertices.first[i].x = line_vertex_data[0];
      drawer.lineVertices.first[i].y = line_vertex_data[1];
      drawer.lineVertices.first[i].z = line_vertex_data[2];
    }

    current_curve = 0;
    curve_min.clear();
    curve_max.clear();
    drawer.fgVertices.clear();

    for (int i = 0; i < drawer.poly3s.size(); i++) {
      for (int j = i + 1; j < drawer.poly3s.size(); j++) {
        displayFGCurves(drawer.poly3s[i], drawer.poly3s[j]);
      }
    }

    for(int i = 0; i < drawer.fgVertices.size(); i++)
    {
      for(int j = 0; j < drawer.fgVertices[i].size(); j++)
      {
        vertex FG_vertex = drawer.fgVertices[i][j];

        FG_vertex_data[0] = FG_vertex.x;
        FG_vertex_data[1] = FG_vertex.y;
        FG_vertex_data[2] =  FG_vertex.z;
        FG_vertex_data[3] = 1.0f;

        FG_vertex_data = PVM_INV * FG_vertex_data;

        drawer.fgVertices[i][j].x = FG_vertex_data[0];
        drawer.fgVertices[i][j].y = FG_vertex_data[1];
        drawer.fgVertices[i][j].z = FG_vertex_data[2];
      }
    }

  duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
  cout << "regrid time " << duration << endl;
      
  return;

}


void assignInverses()
{
  P_Inverse = glm::inverse(drawer.P);
  V_Inverse = glm::inverse(drawer.V);
  return;
}
void restorePoly()
{
    drawer.vbo_vect.clear();
    //drawer.vertex_all.clear();
    drawer.vertex_all_vect.clear();

    glm::mat4 identity_mat = glm::mat4(1.f);

    drawer.M = identity_mat;
    drawer.MLines = identity_mat;
    drawer.transform = drawer.P * drawer.V * drawer.M;


    drawer.transformLines = drawer.P * drawer.V * drawer.MLines;

    for (int i = 0; i < drawer.poly3s.size(); i++) {
      drawer.poly3s[i] = new TransformedPoly3(drawer.poly3s_original[i], identity_mat);
      if (verbose) cout << "transformed poly " << i << endl;
    }
    
      for (int i = 0; i < drawer.poly3s.size(); i++) {
      drawer.AddSurface(drawer.grid, drawer.poly3s[i]);
    }
    
        drawer.bufferData2();
         drawer.lineVertices.second.clear();
         
         for (int i = 0; i < drawer.poly3s.size(); i++) {
      displayChains(drawer.poly3s[i]);
    }


    drawer.lineVertices.first.clear();
    setCubeFaces();

    current_curve = 0;
    curve_min.clear();
    curve_max.clear();
    drawer.fgVertices.clear();

    for (int i = 0; i < drawer.poly3s.size(); i++) {
      for (int j = i + 1; j < drawer.poly3s.size(); j++) {
        displayFGCurves(drawer.poly3s[i], drawer.poly3s[j]);
      }
    }

  updateModel();
  return;
}


/*
void interpolate(){

  if(playback_object->started == false)
  {
    std::cout << "Model has been set to the first pose.\n\n";
    drawer.M = playback_object->getFirstPose();
    drawer.MLines = playback_object->getFirstPose();

    updateModel();
    return;

  }
  else
  {
    if(playback_object->interpolating== false)
    {
      std::cout << "Interpolating...\n\n";
      playback_object->updateM1(drawer.M);
      playback_object->setTransformation();

      playback_object->interpolating = true;

      trs = playback_object->getTRS();
      trsinv = playback_object->getTRSinv();
      pose_model = playback_object->M1;
      pose_model_inverse = drawer.uMIT;
      
      return;
    }
    else
    {
      if(playback_object->counter > playback_object->steps - 1)
      {
        std::cout << "Finished.\n\n";
        std::cout << "M translation before update coords : " << drawer.M[3][0] << " " << drawer.M[3][1] << " " << drawer.M[3][2] << "\n\n";
        drawer.M = playback_object->M2;
        drawer.MLines = playback_object->M2;
        updateModel();

        playback_object->interpolating = false;

        std::cout << "M translation  after update coords : " << drawer.M[3][0] << " " << drawer.M[3][1] << " " << drawer.M[3][2] << "\n\n";

        if(playback_object->isRegrid == true)
        {
          vregrid();
        }
        return;
      }
      else
      {
       
        usleep(30000);

        pose_model = trs * pose_model;

        drawer.M = pose_model;

        drawer.MLines = pose_model;

        drawer.transform = drawer.P * drawer.V * drawer.M;
        drawer.transformLines = drawer.P * drawer.V * drawer.MLines;

        glUniformMatrix4fv(drawer.mvp_location, 1, GL_FALSE, glm::value_ptr(drawer.transform));
        glUniformMatrix4fv(drawer.mvp_location_lines, 1, GL_FALSE, glm::value_ptr(drawer.transformLines));
        glUniformMatrix4fv(drawer.v_location, 1, GL_FALSE, glm::value_ptr(drawer.V));
        glUniformMatrix4fv(drawer.m_location, 1, GL_FALSE, glm::value_ptr(drawer.M));

        //drawer.uMIT = glm::inverse(glm::transpose(drawer.M));


        /*
         * 
         * trsinv = Transpose(trs_inverse)
         *
         * drawer.uMIT = Transpose(M_inverse)
         * 
         * trs inv * drawer.uMIT = Transpose(inverse(trs * M)) 
         *
           */
        /*
        pose_model_inverse = trsinv * pose_model_inverse;
        
        drawer.uMIT = pose_model_inverse;

        glUniformMatrix4fv(drawer.mit_location, 1, GL_FALSE, glm::value_ptr(drawer.uMIT));
        playback_object->counter += 1.0;
        return;
      }
      return;
    }
    
    return;
  }
    return;
}
*/

void main_interpolate(){



  playback_object->interpolate(drawer.M, drawer.MLines, drawer.uMIT);

  if(playback_object->update_model == true){

    updateModel();

    if(playback_object->interpolating == false){

      if(playback_object->isRegrid == true){
        vregrid();
      }

/*
       if(playback_object->keep_interpolating == true){
      //std::cout << "set interp to keep interp\n\n";
      std::cout << "value of keep_interpolating = " << playback_object->keep_interpolating << " \n\n";
      playback_object->interpolate(drawer.M, drawer.MLines, drawer.uMIT);
      return;
    }*/

   }

    playback_object->update_model = false;
    return;
  }
  else{
      
    drawer.transform = drawer.P * drawer.V * drawer.M;
    drawer.transformLines = drawer.P * drawer.V * drawer.MLines;

    glUniformMatrix4fv(drawer.mvp_location, 1, GL_FALSE, glm::value_ptr(drawer.transform));
    glUniformMatrix4fv(drawer.mvp_location_lines, 1, GL_FALSE, glm::value_ptr(drawer.transformLines));
    glUniformMatrix4fv(drawer.v_location, 1, GL_FALSE, glm::value_ptr(drawer.V));
    glUniformMatrix4fv(drawer.m_location, 1, GL_FALSE, glm::value_ptr(drawer.M));
    glUniformMatrix4fv(drawer.mit_location, 1, GL_FALSE, glm::value_ptr(drawer.uMIT));
    return;
  }

     
  
}

void updateModel(){
  
  drawer.transform = drawer.P * drawer.V * drawer.M;
  drawer.transformLines = drawer.P * drawer.V * drawer.MLines;

  glUniformMatrix4fv(drawer.mvp_location, 1, GL_FALSE, glm::value_ptr(drawer.transform));
  glUniformMatrix4fv(drawer.mvp_location_lines, 1, GL_FALSE, glm::value_ptr(drawer.transformLines));
  glUniformMatrix4fv(drawer.v_location, 1, GL_FALSE, glm::value_ptr(drawer.V));
  glUniformMatrix4fv(drawer.m_location, 1, GL_FALSE, glm::value_ptr(drawer.M));

  drawer.uMIT = glm::inverse(glm::transpose(drawer.M));

  glUniformMatrix4fv(drawer.mit_location, 1, GL_FALSE, glm::value_ptr(drawer.uMIT));
}


int main(int argc, char* argv[]) {
  enable();

  // Create a trivariate f(x,y,z)
  // x^4 + y^3 - z - 0.25
  // test: x2+y2+z2=0.5

  // char delim[] = "=";
  // char* word = argv[1];
  // char* ptr = strtok(word, delim);
  // ptr = strtok(NULL, delim);

  //   string a = "../";
  //   a += ptr;

  // const char *cc = a.c_str();
  // const char* cc = ptr;

  ReadBondingBox("boundingBox.txt");

  // The grid and (eventually) its intersections.
  GridRoots grid;

  const char* cmd_input;

  grid.setMinXYZ(los);
  grid.setMaxXYZ(his);
  grid.setSteps(nSteps);

  std::vector<Poly3<Parameter>> polys;
  std::vector<POLY3> poly_objects;

  for (int i = 1; i < argc; i++) {

    cmd_input = argv[i];

    if(cmd_input[0] == '-')
    {
      const char* flag_word = "pose";

      for(int i = 0; i < strlen(flag_word); i++)
      {
        if(cmd_input[i+1] != flag_word[i])
        {
          std::cout << "Unrecognized flag. Usage: -<flag>=<filename>\n\n";
          exit(0);
        }
      }

      interpolate_file = &cmd_input[strlen("-pose=")];
      playback_object = new Playback(interpolate_file);
      continue;
    }

    Poly3<Parameter> p;
    ReadPolyText(argv[i], &p);
    POLY3 p_object = new Object<Poly3>(p);

    polys.push_back(p);
    poly_objects.push_back(p_object);

    displayChains(p_object);
    drawer.addPoly3Original(p_object);
    drawer.addPoly3(p_object);
    drawer.AddSurface(grid, p_object);
  }

  current_curve = 0;
  curve_min.clear();
  curve_min.clear();

  for (int i = 0; i < poly_objects.size(); i++) {
    for (int j = i + 1; j < poly_objects.size(); j++) {
      displayFGCurves(poly_objects[i], poly_objects[j]);
    }
  }

  drawer.initGLFW();
  drawer.bufferData2();
  drawer.initShader();
  drawer.initTexture();
  drawer.initShaderSphere();
  drawer.initShaderLine();
  drawer.runGLFW();
}

