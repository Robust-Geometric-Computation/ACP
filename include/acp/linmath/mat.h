#ifndef MAT_H
#define MAT_H

#include <iostream>
#include "acp/linmath/pv_d.h"

using namespace std;

namespace acp {

class Mat2 {
 public:
  Mat2() {}

  Mat2(const double& ia, const double& ib, const double& ic, const double& id)
      : a(ia), b(ib), c(ic), d(id) {}

  Mat2(const Mat2& m) : a(m.a), b(m.b), c(m.c), d(m.d) {}

  Mat2& operator=(const Mat2& m) {
    a = m.a;
    b = m.b, c = m.c, d = m.d;
    return *this;
  }

  Mat2 operator+(const double k) { return Mat2(a + k, b + k, c + k, d + k); }
  Mat2 operator-(const double k) { return Mat2(a - k, b - k, c - k, d - k); }
  Mat2 operator*(const double k) { return Mat2(a * k, b * k, c * k, d * k); }
  Mat2 operator/(const double k) { return Mat2(a / k, b / k, c / k, d / k); }

  PV2D operator*(const PV2D& v) {
    return PV2D(a * v.x + b * v.y, c * v.x + d * v.y);
  }

  Mat2 operator+(const Mat2& B) {
    return Mat2(a + B.a, b + B.b, c + B.c, d + B.d);
  }

  Mat2 operator-(const Mat2& B) {
    return Mat2(a - B.a, b - B.b, c - B.c, d - B.d);
  }

  Mat2 operator/(const Mat2& B) {
    return Mat2(a / B.a, b / B.b, c / B.c, d / B.d);
  }

  Mat2 operator*(const Mat2& B) {
    double aa, bb, cc, dd;

    aa = a * B.a + b * B.c;
    bb = a * B.b + b * B.d;
    cc = c * B.a + d * B.c;
    dd = c * B.b + d * B.d;

    return Mat2(aa, bb, cc, dd);
  }

  Mat2 inverse() {
    double det = (a * d) - (b * c);

    return Mat2(d, -b, -c, a) * (1 / det);
  }

  void print() {
    cout << "| " << a << " " << b << " |\n"
         << "| " << c << " " << d << " |\n";
  }

 protected:
  double a, b, c, d;
};

inline Mat2 operator*(double k, Mat2 m) { return m * k; }

}  // namespace acp

#endif
