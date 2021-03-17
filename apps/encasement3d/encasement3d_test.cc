#ifdef OSX
#include <Accelerate/Accelerate.h>
#else
extern "C" {
#include <clapack.h>
}
#endif

#include <algorithm>
#include <vector>

#include "acp/encasement3d/encasement3d.h"
#include "acp/encasement3d/encasement_utils.h"
#include "acp/linmath/pv.h"
#include "acp/poly/poly3d.h"

namespace acp {
namespace encasement3d_test {

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

}  // namespace encasement3d_test
}  // namespace acp

#include "acp/encasement3d/mesh3d.h"

class UnitNormal : public Object<PV3> {
  PTR<Object<Poly3D>> p;
  PTR<Object<PV3>> v;
  DeclareCalculate(PV3) {
    return p->get<N>().gradient(v->get<N>()).unit();
  }
public:
  UnitNormal (PTR<Object<Poly3D>> p, PTR<Object<PV3>> v)
    : p(p), v(v) {}
};

class PV3withZeros : public Object<PV3> {
  class V3I : public Primitive {
    PTR<Object<PV3>> p;
    int i;
    DeclareSign {
      return p->get<N>()[i];
    }
  public:
    V3I (PTR<Object<PV3>> p, int i) : p(p), i(i) {}
  };

  PTR<Object<PV3>> pv3;
  bool zero[3];

  DeclareCalculate(PV3) {
    PV3<N> p = pv3->get<N>();
    for (int i = 0; i < 3; ++i)
      if (zero[i])
        p[i] = N(0);
    return p;
  }
    
public:
  PV3withZeros (PTR<Object<PV3>> p) : pv3(p) {
    for (int i = 0; i < 3; ++i)
      zero[i] = (V3I(p, i) == 0);
  }
};

PTR<Object<Poly3D>> makeTorus () {
#ifdef BLEEN
  Poly3D<Parameter> ppoly(2, 2, 2);
  ppoly.set(0, 0, 0, Parameter::input(-0.25));
  ppoly.set(1, 0, 0, Parameter::input(0.0001));
  ppoly.set(0, 1, 0, Parameter::input(0.0002));
  ppoly.set(0, 0, 1, Parameter::input(0.0003));
  ppoly.set(0, 1, 1, Parameter::input(0.0004));
  ppoly.set(1, 0, 1, Parameter::input(0.0005));
  ppoly.set(1, 1, 0, Parameter::input(0.0006));
  ppoly.set(2, 0, 0, Parameter::input(1));
  ppoly.set(0, 2, 0, Parameter::input(1));
  ppoly.set(0, 0, 2, Parameter::input(1));
#endif

  acp::disable();
  double rot = 0.1;
  double len = sqrt(1 + rot * rot);
  double c = 1 / len;
  double s = rot / len;
  acp::enable();

  // (𝑥2+𝑦2+𝑧2+𝑅2−𝑟2)2−4𝑅2(𝑥2+𝑦2)
  Parameter r = Parameter::input(0.1);
  Parameter R = Parameter::input(0.3);

  Poly3D<Parameter> pxyzRr(2, 2, 2);
  pxyzRr.set(0, 0, 0, R * R - r * r);
  pxyzRr.set(2, 0, 0, Parameter::input(1));
  pxyzRr.set(0, 2, 0, Parameter::input(1));
  // pxyzRr.set(0, 0, 2, Parameter::input(1));
  // r -> r (1 - z)
  // -r^2 = -r^2 + 2 r^2 z - r^2 z^2
  pxyzRr.set(0, 0, 1, r * r * 2);
  pxyzRr.set(0, 0, 2, Parameter::input(1) - r * r);
  
  Poly3D<Parameter> torus = pxyzRr * pxyzRr;
  Parameter m4R2 = R * R * -4;
  // torus.add(2, 0, 0, m4R2);
  // x^2 -> c^2 x^2 + 2 c s x y + s^2 y^2
  torus.add(2, 0, 0, m4R2 * c * c);
  torus.add(1, 1, 0, m4R2 * 2 * c * s);
  torus.add(0, 2, 0, m4R2 * s * s);
  torus.add(0, 0, 2, m4R2);

  int n = 20;
  for (int i = 0; i < n; ++i) {
    double z = -0.5 + i * 1.0 / n;
    PV3<Parameter> p(0, 0, Parameter::constant(z));
    Parameter v = torus.value(p);
    cout << z << " " << v.mid() << endl;
  }

  for (int ix = 0; ix <= 4; ++ix)
    for (int iy = 0; iy <= 4; ++iy)
      for (int iz = 0; iz <= 4; ++iz)
        if (ix + iy + iz <= 4)
          torus.add(ix, iy, iz, Parameter::input(0.0001 * ix + 0.0002 * iy + 0.0004 * iz));

  for (Parameter &p : torus.a)
    p = Parameter::constant(p.mid());

  return new Object<Poly3D>(torus);
}

int main(int argc, char** argv) {
  enable();

  if (argc > 2) {
    srand(atoi(argv[2]));
  }

  PV3<double> box_l(-0.5, -0.5, -0.5);
  PV3<double> box_u(0.5, 0.5, 0.5);

  int degree = (argc > 1 ? atoi(argv[1]) : 3);

  //  Poly3D<Parameter> poly =
  //      acp::encasement3d_test::RandomSurface(degree, box_l, box_u);

  int num_terms = (degree + 1) * (degree + 2) * (degree + 3) / 6;

  //  std::cout << "# num of polys" << std::endl;
  //  std::cout << "1" << std::endl;
  //  std::cout << std::endl;
  //  std::cout << "Poly 2" << std::endl;
  //  std::cout << "# num of terms" << std::endl;
  //  std::cout << num_terms << std::endl;
  //  std::cout << "# terms" << std::endl;
  //
  //  for (int i = 0; i <= degree; i++) {
  //    for (int j = 0; j <= degree - i; j++) {
  //      for (int k = 0; k <= degree - i - j; k++) {
  //        std::cout << i << " " << j << " " << k << " " << poly.get(i, j,
  //        k).mid()
  //                  << std::endl;
  //      }
  //    }
  //  }

  PV3<Parameter> box(
      Parameter::constant(box_l.x).interval(Parameter::constant(box_u.x)),
      Parameter::constant(box_l.y).interval(Parameter::constant(box_u.y)),
      Parameter::constant(box_l.z).interval(Parameter::constant(box_u.z)));

  int npolys = (argc > 3 ? atoi(argv[3]) : 3);

  std::vector<PTR<Object<Poly3D>>> polys;
  for (int i = 0; i < npolys; i++) {
    Poly3D<Parameter> poly =
        acp::encasement3d_test::RandomSurface(degree, box_l, box_u);
    polys.push_back(new Object<Poly3D>(poly));
    //    std::cout << "# num of polys" << std::endl;
    //    std::cout << "1" << std::endl;
    //    std::cout << std::endl;
    //    std::cout << "Poly 2" << std::endl;
    //    std::cout << "# num of terms" << std::endl;
    //    std::cout << num_terms << std::endl;
    //    std::cout << "# terms" << std::endl;
    //
    //    for (int i = 0; i <= degree; i++) {
    //      for (int j = 0; j <= degree - i; j++) {
    //        for (int k = 0; k <= degree - i - j; k++) {
    //          std::cout << i << " " << j << " " << k << " "
    //                    << poly.get(i, j, k).mid() << std::endl;
    //        }
    //      }
    //    }
  }

  if (npolys == 0)
    polys.push_back(makeTorus());

  // polys.clear();

  // test max coord
  //  PV3<Parameter> bbox(5, 4, 3);
  //  assert(maxCoord(bbox, -1) == 0);
  //  assert(maxCoord(bbox, 0) == 1);
  //  assert(maxCoord(bbox, 1) == 0);
  //  bbox = PV3<Parameter>(5, 5, 5);
  //  assert(maxCoord(bbox, -1) == 0);
  //  assert(maxCoord(bbox, 0) == 1);
  //  assert(maxCoord(bbox, 2) == 0);

  //  Poly3D<Parameter> plane1(1, 1, 1);
  //  //  assert(plane1.a.size() == 4);
  //  plane1.a[0] = Parameter::input(0.1);
  //  plane1.a[1] = Parameter::input(1.1);
  //  plane1.a[2] = Parameter::input(0.1);
  //  plane1.a[4] = Parameter::input(0.1);
  //  plane1.a[3] = Parameter::input(0.01);
  //  plane1.a[5] = Parameter::input(0.01);
  //  plane1.a[6] = Parameter::input(0.01);
  //  plane1.a[7] = Parameter::input(0.01);
  //  polys.push_back(new Object<Poly3D>(plane1));

  //  Poly3D<Parameter> plane1b(1, 1, 1);
  //  //  assert(plane1b.a.size() == 4);
  //  plane1b.a[0] = Parameter::input(-0.1);
  //  plane1b.a[1] = Parameter::input(1.1);
  //  plane1b.a[2] = Parameter::input(0.1);
  //  plane1b.a[4] = Parameter::input(0.1);
  //  plane1b.a[3] = Parameter::input(0.01);
  //  plane1b.a[5] = Parameter::input(0.01);
  //  plane1b.a[6] = Parameter::input(0.01);
  //  plane1b.a[7] = Parameter::input(0.01);
  //  polys.push_back(new Object<Poly3D>(plane1b));
  //
  //  Poly3D<Parameter> plane2(1, 1, 1);
  //  //  assert(plane2.a.size() == 4);
  //  plane2.a[0] = Parameter::input(0.1);
  //  plane2.a[1] = Parameter::input(0.1);
  //  plane2.a[2] = Parameter::input(1.1);
  //  plane2.a[4] = Parameter::input(0.1);
  //  plane2.a[3] = Parameter::input(0.01);
  //  plane2.a[5] = Parameter::input(0.01);
  //  plane2.a[6] = Parameter::input(0.01);
  //  plane2.a[7] = Parameter::input(0.01);
  //  polys.push_back(new Object<Poly3D>(plane2));
  //
  //  Poly3D<Parameter> plane2b(1, 1, 1);
  //  //  assert(plane2b.a.size() == 4);
  //  plane2b.a[0] = Parameter::input(-0.1);
  //  plane2b.a[1] = Parameter::input(0.1);
  //  plane2b.a[2] = Parameter::input(1.1);
  //  plane2b.a[4] = Parameter::input(0.1);
  //  plane2b.a[3] = Parameter::input(0.01);
  //  plane2b.a[5] = Parameter::input(0.01);
  //  plane2b.a[6] = Parameter::input(0.01);
  //  plane2b.a[7] = Parameter::input(0.01);
  //  polys.push_back(new Object<Poly3D>(plane2b));
  //
  //  Poly3D<Parameter> plane3(1, 1, 1);
  //  //  assert(plane3.a.size() == 4);
  //  plane3.a[0] = Parameter::input(0.1);
  //  plane3.a[1] = Parameter::input(0.1);
  //  plane3.a[2] = Parameter::input(0.1);
  //  plane3.a[4] = Parameter::input(1.1);
  //  plane3.a[3] = Parameter::input(0.01);
  //  plane3.a[5] = Parameter::input(0.01);
  //  plane3.a[6] = Parameter::input(0.01);
  //  plane3.a[7] = Parameter::input(0.01);
  //  polys.push_back(new Object<Poly3D>(plane3));
  //
  //  Poly3D<Parameter> plane3b(1, 1, 1);
  //  //  assert(plane3b.a.size() == 4);
  //  plane3b.a[0] = Parameter::input(-0.1);
  //  plane3b.a[1] = Parameter::input(0.1);
  //  plane3b.a[2] = Parameter::input(0.1);
  //  plane3b.a[4] = Parameter::input(1.1);
  //  plane3b.a[3] = Parameter::input(0.01);
  //  plane3b.a[5] = Parameter::input(0.01);
  //  plane3b.a[6] = Parameter::input(0.01);
  //  plane3b.a[7] = Parameter::input(0.01);
  //  polys.push_back(new Object<Poly3D>(plane3b));
  //
  //  Poly3D<Parameter> sphere1(2, 2, 2);
  //  for (int ix = 0; ix < 2; ix++)
  //    for (int iy = 0; iy < 2; iy++)
  //      for (int iz = 0; iz < 2; iz++)
  //	sphere1.set(ix, iy, iz, Parameter::input(0.01));
  //  sphere1.set(0, 0, 0, -Parameter::input(0.25));
  //  sphere1.set(2, 0, 0, Parameter::input(1.0));
  //  sphere1.set(0, 2, 0, Parameter::input(1.0));
  //  sphere1.set(0, 0, 2, Parameter::input(1.0));
  //  polys.push_back(new Object<Poly3D>(sphere1));

  //  Poly3D<Parameter> sphere2(2, 2, 2);
  //  for (int ix = 0; ix < 2; ix++)
  //    for (int iy = 0; iy < 2; iy++)
  //      for (int iz = 0; iz < 2; iz++)
  //        sphere2.set(ix, iy, iz, Parameter::input(0.01));
  //
  //  sphere2.set(0, 0, 0, -Parameter::input(0.1));
  //  sphere2.set(2, 0, 0, Parameter::input(1.0));
  //  sphere2.set(0, 2, 0, Parameter::input(1.0));
  //  sphere2.set(0, 0, 2, Parameter::input(1.0));
  //  polys.push_back(new Object<Poly3D>(sphere2));

  Encasement3D* encasement = new Encasement3D(box, polys);
  cout << encasement->cells.size() << endl;
  cout << encasement->f[6]->faces.size() << endl;;
  class Mesh3D mesh3d(encasement, 1e-3);
  
  cout << "meshing" << endl;
  for (int h = 6; h < encasement->f.size(); ++h) {
    for (int i = 0; i < encasement->f[h]->faces.size(); ++i) {
      if (i == 8 || i == 9 || i == 10 || i == 11)
	/* continue */;
      cout << "i = " << i << endl;
      const vector<Mesh3D::Triangle> &mesh = 
        mesh3d.getMesh(encasement->f[h]->faces[i]);
      cout << mesh.size() << endl;
    }
  }

  cout << "premeshing" << endl;
  for (int h = 6; h < encasement->f.size(); ++h) {
    for (int i = 0; i < encasement->f[h]->faces.size(); ++i) {
      if (i == 8 || i == 9 || i == 10 || i == 11)
	/* continue */;
      const vector<Mesh3D::Triangle> &mesh = 
        mesh3d.getMesh(encasement->f[h]->faces[i]);
      cout << mesh.size() << endl;
    }
  }

  int nzeros = 0;

  cout << "remeshing" << endl;
  for (int h = 6; h < encasement->f.size(); ++h) {
    for (int i = 0; i < encasement->f[h]->faces.size(); ++i) {
      if (i == 8 || i == 9 || i == 10 || i == 11)
	/* continue */;
      const vector<Mesh3D::Triangle> &mesh = 
        mesh3d.getMesh(encasement->f[h]->faces[i]);
      cout << mesh.size() << endl;
      std::stringstream ss;
      ss << "face" 
         << std::setfill('0') << std::setw(3) << h-6
         << "-"
         << std::setfill('0') << std::setw(3) << i
         << ".vtk";
      mesh3d.writeMesh(ss.str().c_str(), mesh);
      
      for (const Mesh3D::Triangle &t : mesh)
        for (int j = 0; j < 3; ++j) {
          PTR<Object<PV3>> u = new UnitNormal(encasement->bf[h], t.p[j]);
          PTR<Object<PV3>> u0 = new PV3withZeros(u);
          PV3<Parameter> up = u0->getApprox();
          for (int k = 0; k < 3; ++k)
            if (up[k].lb() == 0 && up[k].ub() == 0)
              ++nzeros;
        }
    }
  }
  
  cout << "nzeros " << nzeros << endl;

#ifdef BLEEN
  PV3<Parameter> pp = PV3<Parameter>::input(0, 0, 0);
  cout << pp.x.lb() << " " << pp.y.lb() << " " << pp.z.lb() << endl;
  cout << pp.x.ub() << " " << pp.y.ub() << " " << pp.z.ub() << endl;
  PTR<Object<PV3>> point = new Object<PV3>(pp);

  Encasement3D::Cell *cell = encasement->getCell(point);

  for (Encasement3D::Shell *shell : cell->shells) {
    for (Encasement3D::Face *face : shell->faces) {
      if (face->jf < 6)
        continue;
      const vector<Mesh3D::Triangle> &mesh = 
        mesh3d.getMesh(face);
      std::stringstream ss;
      ss << "cell" 
         << std::setfill('0') << std::setw(3) << face->jf - 6
         << ".vtk";
      mesh3d.writeMesh(ss.str().c_str(), mesh);
    }
  }
#endif

  disable();
}
