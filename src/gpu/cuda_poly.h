#ifndef CUDA_POLY_H
#define CUDA_POLY_H

#include <vector>
#include "cuda_interval_lib.h"

template <class T>
class cuda_poly2 {
 public:
  __host__ cuda_poly2(interval_gpu<T>* coef, int* deg, int nt);

  __device__ cuda_poly2();

  __device__ interval_gpu<T> value(interval_gpu<T> const& x,
                                   interval_gpu<T> const& y) const;

  __host__ cuda_poly2<T>* copy_device();

  __host__ __device__ static interval_gpu<T> value(interval_gpu<T>* a, int* m,
                                                   int nt, int degx, int degy,
                                                   interval_gpu<T> const& x,
                                                   interval_gpu<T> const& y) {
    interval_gpu<T> xp[400];
    interval_gpu<T> yp[400];

    interval_gpu<T> p = x;
    for (int i = 1; i <= degx; i++) {
      xp[i] = p;
      p = p * x;
    }

    p = y;
    for (int i = 1; i <= degy; i++) {
      yp[i] = p;
      p = p * y;
    }

    interval_gpu<T> v;

    for (int i = 0; i < nt; i++) {
      interval_gpu<T> z = a[i];
      int xd = m[2 * i];
      int yd = m[2 * i + 1];
      if (xd > 0) z = z * xp[xd];
      if (yd > 0) z = z * yp[yd];

      if (i == 0)
        v = z;
      else
        v = v + z;
    }

    return v;
  }

  interval_gpu<T>* a;
  int* m;
  int degx;
  int degy;
  int nt;
  int d;
};

template <class T>
inline __host__ cuda_poly2<T>* cuda_poly2<T>::copy_device() {
  cuda_poly2<T>* dev_f;

  interval_gpu<T>* aptr;
  int* mptr;

  cudaMalloc((void**)&dev_f, sizeof(cuda_poly2<T>));

  cudaMalloc((void**)&aptr, sizeof(interval_gpu<T>) * nt);
  cudaMalloc((void**)&mptr, sizeof(int) * 2 * nt);

  cudaMemcpy(aptr, a, sizeof(interval_gpu<T>) * nt, cudaMemcpyHostToDevice);
  cudaMemcpy(mptr, m, sizeof(int) * 2 * nt, cudaMemcpyHostToDevice);

  cudaMemcpy(&(dev_f->a), &aptr, sizeof(interval_gpu<T>*),
             cudaMemcpyHostToDevice);
  cudaMemcpy(&(dev_f->m), &mptr, sizeof(int*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(dev_f->nt), &nt, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(&(dev_f->degx), &degx, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(&(dev_f->degy), &degy, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(&(dev_f->d), &d, sizeof(int), cudaMemcpyHostToDevice);

  return dev_f;
}

template <class T>
inline __host__ cuda_poly2<T>::cuda_poly2(interval_gpu<T>* coef, int* deg,
                                          int num) {
  nt = num;

  a = new interval_gpu<T>[nt];
  m = new int[2 * nt];

  for (int i = 0; i < nt; i++) {
    a[i] = coef[i];
    m[2 * i] = deg[2 * i];
    m[2 * i + 1] = deg[2 * i + 1];
  }

  d = degx = degy = -1;

  for (int i = 0; i < nt; i++) {
    degx = max(m[2 * i], degx);
    degy = max(m[2 * i + 1], degy);
    d = max(d, degx + degy);
  }
}

template <class T>
inline __device__ interval_gpu<T> cuda_poly2<T>::value(
    interval_gpu<T> const& x, interval_gpu<T> const& y) const {
  interval_gpu<T> xp[400];
  interval_gpu<T> yp[400];

  interval_gpu<T> p = x;
  for (int i = 1; i <= degx; i++) {
    xp[i] = p;
    p = p * x;
  }

  p = y;
  for (int i = 1; i <= degy; i++) {
    yp[i] = p;
    p = p * y;
  }

  interval_gpu<T> v;

  for (int i = 0; i < nt; i++) {
    interval_gpu<T> z = a[i];
    int xd = m[2 * i];
    int yd = m[2 * i + 1];
    if (xd > 0) z = z * xp[xd];
    if (yd > 0) z = z * yp[yd];

    if (i == 0)
      v = z;
    else
      v = v + z;
  }

  return v;
}

#endif
