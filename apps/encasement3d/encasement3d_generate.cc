#ifdef OSX
#include <Accelerate/Accelerate.h>
#else
extern "C" {
#include <clapack.h>
}
#endif

#include <algorithm>
#include <chrono>
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
                                const PV3<double>& box_u,
                                std::vector<PV3<double>> points) {
  int num_points = ((degree + 1) * (degree + 2) * (degree + 3) / 6) - 1;

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
            // std::cout << " m[" << (c - 1) << "*" << num_points << " + " <<
            // k
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
          // poly.set(i, j, k, Parameter::constant(b[c++]));
          poly.set(i, j, k, Parameter::input(b[c++]));
        }
      }
    }
  }

  return poly;
}

Poly3D<Parameter> RandomSurface(int degree, const PV3<double>& box_l,
                                const PV3<double>& box_u) {
  int num_points = ((degree + 1) * (degree + 2) * (degree + 3) / 6) - 1;

  std::vector<PV3<double>> points = RandomPoints(box_l, box_u, num_points);

  return RandomSurface(degree, box_l, box_u, points);
}

}  // namespace encasement3d_test
}  // namespace acp

int main(int argc, char** argv) {
  enable();

  if (argc > 3) {
    srand(atoi(argv[3]));
  }

  PV3<double> box_l(-1, -1, -1);
  PV3<double> box_u(1, 1, 1);

  int degree = (argc > 1 ? atoi(argv[1]) : 4);

  int number_of_surfaces = (argc > 2 ? atoi(argv[2]) : 3);

  Poly3D<Parameter> poly =
      acp::encasement3d_test::RandomSurface(degree, box_l, box_u);

  int num_terms = (degree + 1) * (degree + 2) * (degree + 3) / 6;

  std::cout << "# num of polys" << std::endl;
  std::cout << "1" << std::endl;
  std::cout << std::endl;
  std::cout << "Poly 2" << std::endl;
  std::cout << "# num of terms" << std::endl;
  std::cout << num_terms << std::endl;
  std::cout << "# terms" << std::endl;

  for (int i = 0; i <= degree; i++) {
    for (int j = 0; j <= degree - i; j++) {
      for (int k = 0; k <= degree - i - j; k++) {
        std::cout << i << " " << j << " " << k << " " << poly.get(i, j, k).mid()
                  << std::endl;
      }
    }
  }

  // Random surfaces in a box

  std::cout << number_of_surfaces << " surfaces of degree " << degree
            << std::endl;

  PV3<Parameter> box(
      Parameter::constant(box_l.x).interval(Parameter::constant(box_u.x)),
      Parameter::constant(box_l.y).interval(Parameter::constant(box_u.y)),
      Parameter::constant(box_l.z).interval(Parameter::constant(box_u.z)));

  std::vector<PTR<Object<Poly3D>>> polys;
  for (int i = 0; i < number_of_surfaces; i++) {
    polys.push_back(new Object<Poly3D>(
        acp::encasement3d_test::RandomSurface(degree, box_l, box_u)));
    std::cout << "# num of polys" << std::endl;
    std::cout << "1" << std::endl;
    std::cout << std::endl;
    std::cout << "Poly 2" << std::endl;
    std::cout << "# num of terms" << std::endl;
    std::cout << num_terms << std::endl;
    std::cout << "# terms" << std::endl;

    for (int i = 0; i <= degree; i++) {
      for (int j = 0; j <= degree - i; j++) {
        for (int k = 0; k <= degree - i - j; k++) {
          std::cout << i << " " << j << " " << k << " "
                    << poly.get(i, j, k).mid() << std::endl;
        }
      }
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
  }

  if (argc > 5) {
    double l = atof(argv[4]);
    double u = atof(argv[5]);
    box_l = PV3<double>(l, l, l);
    box_u = PV3<double>(u, u, u);

    box = PV3<Parameter>(
        Parameter::constant(box_l.x).interval(Parameter::constant(box_u.x)),
        Parameter::constant(box_l.y).interval(Parameter::constant(box_u.y)),
        Parameter::constant(box_l.z).interval(Parameter::constant(box_u.z)));

    std::cout << "box(" << l << ", " << u << ")" << std::endl;
  }

  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < number_of_surfaces; i++) {
    for (int j = i + 1; j < number_of_surfaces; j++) {
      for (int k = j + 1; k < number_of_surfaces; k++) {
        std::vector<PTR<Object<PV3>>> roots =
            getRoots(polys[i], polys[j], polys[k], box);

        std::cout << "Intersections(" << i << ", " << j << ", " << k
                  << "): " << roots.size() << std::endl;
      }
    }
  }

  auto stop = std::chrono::high_resolution_clock::now();

  auto duration =
      std::chrono::duration_cast<std::chrono::seconds>(stop - start);

  std::cout << duration.count() << std::endl;

  disable();
}
