#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/remove.h>
#include <thrust/sequence.h>

#include <stdio.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "cuda_interval.h"
#include "cuda_poly.h"

#include "../critical_points.h"
#include "../dcel.h"
#include "../encasement.h"
#include "../ppoly.h"

using namespace std;

typedef double Number;

template <typename T>
struct i_pair {
  interval_gpu<T> x;
  interval_gpu<T> y;
};

// this functor returns true if the argument is constant sign on f and g, and
// false otherwise
template <typename T>
struct is_ambiguous : public thrust::unary_function<T, bool> {
  const cuda_poly2<T>* f;
  const cuda_poly2<T>* g;

  is_ambiguous(const cuda_poly2<T>* _f, const cuda_poly2<T>* _g)
      : f(_f), g(_g) {}

  __host__ __device__ bool operator()(i_pair<T> I) {
    // evaluate f and g on I and return true if both are non constant
    interval_gpu<T> fv = f->value(I.x, I.y);
    interval_gpu<T> gv = g->value(I.x, I.y);

    return fv.lower() <= 0 && fv.upper() >= 0 && gv.lower() <= 0 &&
           gv.upper() >= 0;
  }
};

// this functor returns true if the argument is constant sign on f and g, and
// false otherwise
template <typename T>
struct is_ambiguous_quad : public thrust::unary_function<T, bool> {
  const cuda_poly2<T>* f;
  const cuda_poly2<T>* fx;
  const cuda_poly2<T>* fy;
  const cuda_poly2<T>* fxx;
  const cuda_poly2<T>* fxy;
  const cuda_poly2<T>* fyy;

  const cuda_poly2<T>* g;
  const cuda_poly2<T>* gx;
  const cuda_poly2<T>* gy;
  const cuda_poly2<T>* gxx;
  const cuda_poly2<T>* gxy;
  const cuda_poly2<T>* gyy;

  is_ambiguous_quad(const cuda_poly2<T>* _f, const cuda_poly2<T>* _fx,
                    const cuda_poly2<T>* _fy, const cuda_poly2<T>* _fxx,
                    const cuda_poly2<T>* _fxy, const cuda_poly2<T>* _fyy,
                    const cuda_poly2<T>* _g, const cuda_poly2<T>* _gx,
                    const cuda_poly2<T>* _gy, const cuda_poly2<T>* _gxx,
                    const cuda_poly2<T>* _gxy, const cuda_poly2<T>* _gyy)
      : f(_f),
        fx(_fx),
        fy(_fy),
        fxx(_fxx),
        fxy(_fxy),
        fyy(_fyy),
        g(_g),
        gx(_gx),
        gy(_gy),
        gxx(_gxx),
        gxy(_gxy),
        gyy(_gyy) {}

  __host__ __device__ int sign(interval_gpu<T> p) {
    return p.lower() < 0 ? (p.upper() > 0 ? 0 : -1) : 1;
  }

  __host__ __device__ bool sign(int index, i_pair<T> I) {
    const cuda_poly2<T>* f;
    const cuda_poly2<T>* fx;
    const cuda_poly2<T>* fy;
    const cuda_poly2<T>* fxx;
    const cuda_poly2<T>* fxy;
    const cuda_poly2<T>* fyy;

    interval_gpu<T> a[6];
    int m[12] = {0, 0, 1, 0, 0, 1, 2, 0, 1, 1, 0, 2};
    int nt = 6;
    int degx = 2;
    int degy = 2;

    // storage for substituted univariates and derivatives
    interval_gpu<T> sub[3];

    interval_gpu<T> der_sub[2];

    if (index == 0) {
      f = this->f;
      fx = this->fx;
      fy = this->fy;
      fxx = this->fxx;
      fxy = this->fxy;
      fyy = this->fyy;
    } else {
      f = this->g;
      fx = this->gx;
      fy = this->gy;
      fxx = this->gxx;
      fxy = this->gxy;
      fyy = this->gyy;
    }

    // evaluate the poly on the four corners
    interval_gpu<T> fll = f->value(I.x.lower(), I.y.lower());
    interval_gpu<T> flr = f->value(I.x.upper(), I.y.lower());
    interval_gpu<T> ful = f->value(I.x.lower(), I.y.upper());
    interval_gpu<T> fur = f->value(I.x.upper(), I.y.upper());

    // calculate the sign of f on the corners
    int sll = sign(fll);
    int slr = sign(flr);
    int sul = sign(ful);
    int sur = sign(fur);

    int sum = sll + slr + sul + sur;

    // if it's not consistent across all corners, fail
    if (abs(sum) != 4) return true;

    // now the function is translated to (0, 0) to (w, h)
    double w = I.x.upper() - I.x.lower();
    double h = I.y.upper() - I.y.lower();

    double lx = I.x.lower();
    double ly = I.y.lower();

    // sll is either 1 or -1 and tells us whether to use f or -f for minimum of
    // quad

    a[0] = (T(sll) * f->value(lx, ly)).lower();
    a[1] = (T(sll) * fx->value(lx, ly)).lower();
    a[2] = (T(sll) * fy->value(lx, ly)).lower();
    a[3] = (T(sll) * 0.5 * fxx->value(I.x, ly)).lower();
    a[4] = (T(sll) * fxy->value(I.x, ly)).lower();
    a[5] = (T(sll) * 0.5 * fyy->value(I.x, I.y)).lower();

    // evaluate sign of quadratic on corners of translated box

    interval_gpu<T> qll = cuda_poly2<T>::value(a, m, nt, degx, degy, 0, 0);
    interval_gpu<T> qlr = cuda_poly2<T>::value(a, m, nt, degx, degy, w, 0);
    interval_gpu<T> qur = cuda_poly2<T>::value(a, m, nt, degx, degy, w, h);
    interval_gpu<T> qul = cuda_poly2<T>::value(a, m, nt, degx, degy, 0, h);

    int sll_z = sign(qll);
    int slr_z = sign(qlr);
    int sur_z = sign(qur);
    int sul_z = sign(qul);

    int qsum = sll_z + slr_z + sur_z + sul_z;

    if (qsum < 4) return true;

    // check that substituting each side of the box [0, 0]x[w, h] doesn't have
    // any roots (check sign of derivative on endpoints)

    // quad   = a0 + a1x + a2y + a3 x^2 + a4 xy + a5 y^2

    //------------------------------------------------------------------------------------------------------------------------

    // sub_y = (a0 + a2(y) + a5 (y)^2) + (a1 + a4(y)) x + a3 x^2

    // sub lower y = 0
    sub[0] = a[0];
    sub[1] = a[1];
    sub[2] = a[3];

    // compute derivative from substituted
    der_sub[0] = sub[1];
    der_sub[1] = 2.0 * sub[2];

    // evaluate derivative on endpoints to see if min in the interval
    interval_gpu<T> left_d_val = der_sub[0];                      // x = 0
    interval_gpu<T> right_d_val = der_sub[0] + (w * der_sub[1]);  // x = w

    // min exists within interval, evaluate the sub quadratic at the min
    if (sign(left_d_val) != sign(right_d_val)) {
      interval_gpu<T> min_x = (-1.0 * der_sub[0]) / der_sub[1];
      interval_gpu<T> min_val =
          sub[0] + (min_x * sub[1]) + (min_x * min_x * sub[2]);
      if (sign(min_val) <= 0) return true;
    }

    //------------------------------------------------------------------------------------------------------------------------

    // sub_y = (a0 + a2(y) + a5 (y)^2) + (a1 + a4(y)) x + a3 x^2

    // sub upper y = h
    sub[0] = a[0] + (a[2] * h) + (a[5] * h * h);
    sub[1] = a[1] + (a[4] * h);
    sub[2] = a[3];

    // compute derivative from substituted
    der_sub[0] = sub[1];
    der_sub[1] = 2.0 * sub[2];

    // evaluate derivative on endpoints to see if min in the interval
    left_d_val = der_sub[0];                      // x = 0
    right_d_val = der_sub[0] + (w * der_sub[1]);  // x = w

    // min exists within interval, evaluate the sub quadratic at the min
    if (sign(left_d_val) != sign(right_d_val)) {
      interval_gpu<T> min_x = (-1.0 * der_sub[0]) / der_sub[1];
      interval_gpu<T> min_val =
          sub[0] + (min_x * sub[1]) + (min_x * min_x * sub[2]);
      if (sign(min_val) <= 0) return true;
    }

    //------------------------------------------------------------------------------------------------------------------------

    // sub_x = (a0 + a1(x) + a3 (x)^2) + (a2 + a4(x)) y + a5 y^2

    // sub left  x = 0
    sub[0] = a[0];
    sub[1] = a[2];
    sub[2] = a[5];

    // compute derivative from substituted
    der_sub[0] = sub[1];
    der_sub[1] = 2.0 * sub[2];

    // evaluate derivative on endpoints to see if min in the interval
    left_d_val = der_sub[0];                      // y = 0
    right_d_val = der_sub[0] + (h * der_sub[1]);  // y = h

    // min exists within interval, evaluate the sub quadratic at the min
    if (sign(left_d_val) != sign(right_d_val)) {
      interval_gpu<T> min_y = (-1.0 * der_sub[0]) / der_sub[1];
      interval_gpu<T> min_val =
          sub[0] + (min_y * sub[1]) + (min_y * min_y * sub[2]);
      if (sign(min_val) <= 0) return true;
    }

    //------------------------------------------------------------------------------------------------------------------------

    // sub_x = (a0 + a1(x) + a3 (x)^2) + (a2 + a4(x)) y + a5 y^2

    // sub right x = w
    sub[0] = a[0] + (a[1] * w) + (a[3] * w * w);
    sub[1] = a[2] + (a[4] * w);
    sub[2] = a[5];

    // compute derivative from substituted
    der_sub[0] = sub[1];
    der_sub[1] = 2.0 * sub[2];

    // evaluate derivative on endpoints to see if min in the interval
    left_d_val = der_sub[0];                      // y = 0
    right_d_val = der_sub[0] + (h * der_sub[1]);  // y = h

    // min exists within interval, evaluate the sub quadratic at the min
    if (sign(left_d_val) != sign(right_d_val)) {
      interval_gpu<T> min_y = (-1.0 * der_sub[0]) / der_sub[1];
      interval_gpu<T> min_val =
          sub[0] + (min_y * sub[1]) + (min_y * min_y * sub[2]);
      if (sign(min_val) <= 0) return true;
    }

    //-------------------------------------------------------------------------------------------------------------------------

    // quad   = a0 + a1x + a2y + a3 x^2 + a4 xy + a5 y^2
    // quad_x = a1 + 2 a3 x +   a4 y
    // quad_y = a2 +   a4 x + 2 a5 y

    // constant coefficients
    interval_gpu<T> c1 = a[1];
    interval_gpu<T> c2 = a[2];

    // coefficients of x
    interval_gpu<T> a1 = 2.0 * a[3];
    interval_gpu<T> a2 = a[4];

    // coefficients of y
    interval_gpu<T> b1 = a[4];
    interval_gpu<T> b2 = 2.0 * a[5];

    interval_gpu<T> critx = (b1 * c2 - c1 * b2) / (a1 * b2 - b1 * a2);
    interval_gpu<T> crity = (c1 * a2 - a1 * c2) / (a1 * b2 - b1 * a2);

    if (0.0 < critx && critx < w && 0.0 < crity && crity < h) {
      interval_gpu<T> crit_val = f->value(critx, crity);
      if (sign(crit_val) < 0) return true;
    }

    return false;
  }

  __host__ __device__ bool operator()(i_pair<T> I) {
    // evaluate f and g on I and return true if both are non constant
    interval_gpu<T> fv = f->value(I.x, I.y);
    interval_gpu<T> gv = g->value(I.x, I.y);

    // if (fv.lower() > 0 || fv.upper() < 0 || gv.lower() > 0
    //		|| gv.upper() < 0) {
    //	return false;
    //} else {
    return sign(0, I) && sign(1, I);
    //}
  }
};

template <typename T>
struct function_value : public thrust::unary_function<T, T> {
  const cuda_poly2<T>* f;
  const cuda_poly2<T>* g;

  function_value(const cuda_poly2<T>* _f, const cuda_poly2<T>* _g)
      : f(_f), g(_g) {}

  __host__ __device__ i_pair<T> operator()(i_pair<T> I) {
    // evaluate f and g on I and return true if both are non constant
    interval_gpu<T> fv = f->value(I.x, I.y);
    interval_gpu<T> gv = g->value(I.x, I.y);

    i_pair<T> p;
    p.x = fv;
    p.y = gv;

    return p;
  }
};

template <typename T>
struct subdivide : public thrust::unary_function<int, i_pair<T>> {
  int total;

  i_pair<T>* intervals;

  int nintervals;

  subdivide(int _total, i_pair<T>* _intervals, int _nintervals)
      : total(_total), intervals(_intervals), nintervals(_nintervals) {}

  __device__ i_pair<T> operator()(const int& i) const {
    int sub_per_int = total / nintervals;
    int rows = int(sqrtf(sub_per_int));
    int cols = sub_per_int / rows;

    sub_per_int = rows * cols;

    int interval_index = i / sub_per_int;

    if (interval_index >= nintervals)
      return {interval_gpu<T>(T(0), T(0)), interval_gpu<T>(T(0), T(0))};

    i_pair<T> ip = intervals[interval_index];

    // split more on the larger axis
    if ((ip.x.upper() - ip.x.lower()) < (ip.y.upper() - ip.y.lower())) {
      int temp = rows;
      rows = cols;
      cols = temp;
    }

    interval_gpu<T> xi = ip.x;
    interval_gpu<T> yi = ip.y;

    int j = i % sub_per_int;

    int x = j / cols;
    int y = j % cols;

    T xil = xi.lower();
    T xiu = xi.upper();
    T xiw = xiu - xil;

    T yil = yi.lower();
    T yiu = yi.upper();
    T yiw = yiu - yil;

    T xl = xil + (xiw * (T(x) / rows));
    T yl = yil + (yiw * (T(y) / cols));

    T xu = xil + (xiw * (T(x + 1) / rows));
    T yu = yil + (yiw * (T(y + 1) / cols));

    if (x == rows - 1) xu = xiu;
    if (y == cols + 1) yu = yiu;

    i_pair<T> ret;
    ret.x = interval_gpu<T>(xl, xu);
    ret.y = interval_gpu<T>(yl, yu);

    return ret;
  }
};

vector<PTR<Object<PPoly2>>> loadCurves(const char* fname);

thrust::host_vector<i_pair<Number>> get_regions_gpu(PTR<Object<PPoly2>> f_in,
                                                    PTR<Object<PPoly2>> g_in,
                                                    double xl, double yl,
                                                    double xu, double yu,
                                                    int sub, int iterations);

vector<Rectangle*> get_all_rects(PTR<Object<PPoly2>> f, PTR<Object<PPoly2>> g,
                                 double xl, double yl, double xu, double yu,
                                 int sub, int iterations);

vector<Rectangle*> get_regions(PTR<Object<PPoly2>> f, PTR<Object<PPoly2>> g,
                               double xl, double yl, double xu, double yu,
                               int sub, int iterations);
