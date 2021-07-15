/*
  ACP (Adaptive Controlled Precision) Library
  for robust computational geometry

  Copyright (c) 2013-07-15
  Victor Joseph Milenkovic
  University of Miami
  vjm@miami.edu
  Elisha Sacks
  Purdue University
  eps@cs.purdue.edu

  This file contains code described in

  Robust Complete Path Planning in the Plane
  Victor Milenkovic, Elisha Sacks, and Steven Trac
  Proceedings of the Workshop on the Algorithmic Foundations of Robotics (WAFR)
  pages 37-52, 2012

  This code is under development.
  It is free for use for academic and research purposes.
  Permission from the authors is required for commercial use.
*/

#include "root.h"

namespace acp {

vector<PTR<Object<Scalar>>> getRoots(PTR<Object<Poly>> poly) {
  return PolySolver(poly).getRoots();
}

vector<PTR<Object<Scalar>>> getRoots(PTR<Object<Poly>> poly,
                                     PTR<Object<Scalar>> lb,
                                     PTR<Object<Scalar>> ub) {
  return PolySolver(poly).getRoots(lb, ub);
}

template<class T>
T Poly<T>::shrinkBracket (T x, int sub) const {
  bool bflag = false;
  while (true) {
    bool flag = true;
    T xm = x.midP(), fm = value(xm), df = der(x);
    if (df.sign(false) == 0)
      flag = false;
    else {
      T nx = (xm - fm/df).intersect(x);
      if (nx.subset(x))
        x = nx;
      else
        flag = false;
    }
    if (!flag) {
      int sm = fm.sign(false);
      if (sm == sub) {
        T nx = x.interval(xm);
        if (!nx.subset(x))
          break;
        x = nx;
      }
      else if (sm == - sub) {
        T nx = xm.interval(x);
        if (!nx.subset(x))
          break;
        x = nx;
      }
      else if (!bflag) {
        bflag = true;
        x = shrinkBracketB(x, sub);
      }
      else
        break;
    }
  }
#ifdef SEPBOUND
  if (std::is_same<T, MParameter>::value)
    MParameter::polySep(a, x);
#endif
  return x;
}

template<class T>
T Poly<T>::shrinkBracketB (const T &x, int sub) const {
  T xlb = x.lbP(), xm = x.midP(), xub = x.ubP();
  while ((xlb - xm).sign(false) < 0) {
    T nx = xlb.interval(xm).midP();
    if ((nx - xlb).sign(false) == 0 || (nx - xm).sign(false) == 0)
      break;
    int sm = value(nx).sign(false);
    if (sm == 0)
      xm = nx;
    else if (sm == sub) {
      xub = nx;
      T nx = xlb.interval(xub);
      return nx.subset(x) ? nx : x;
    }
    else
      xlb = nx;
  }
  xm = x.midP();
  while ((xm - xub).sign(false) < 0) {
    T nx = xm.interval(xub).midP();
    if ((nx - xm).sign(false) == 0 || (nx - xub).sign(false) == 0)
      break;
    int sm = value(nx).sign(false);
    if (sm == 0)
      xm = nx;
    else if (sm == -sub) {
      xlb = nx;
      T nx = xlb.interval(xub);
      return nx.subset(x) ? nx : x;
    }
    else
      xub = nx;
  }
  T nx = xlb.interval(xub);
  return nx.subset(x) ? nx : x;
}

template class Poly<Parameter>;
template class Poly<PParameter>;
template class Poly<MParameter>;
}
