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

#ifndef ROOT2_H
#define ROOT2_H

#include "acp/core/object.h"
#include "acp/poly/poly2.h"

namespace acp {

class Root2;

class Root2PTR : public PTR<Root2> {
 public:
  Root2PTR() {}
  Root2PTR(Root2* r) : PTR<Root2>(r) {}
  operator PTR<Object<PV2>>() const;
};

class Root2 : public Object<PV2> {
  PV2<Parameter> r;
  PTR<Object<Poly2>> f, g;

 public:
  Root2(PTR<Object<Poly2>> f, PTR<Object<Poly2>> g, const PV2<Parameter>& r)
      : r(r), f(f), g(g) {}
  PTR<Object<Poly2>> getF() { return f; }
  PTR<Object<Poly2>> getG() { return g; }

 private:
  DeclareCalculate(PV2) {
    if (uninitialized())
      return polish<N>(r);
    else
      return polish<N>(getCurrentP());
  }

  // order of roots in v: bot, right, top, left
  // 0-3 for f, 4-7 for g
  template <class N>
  class RootBoundary {
   public:
    PV2<N> I;
    vector<vector<PTR<Object<Scalar>>>> v;
  };

  /*
  //order of roots in v: bot, right, top, left
  typedef struct {
    PV2 I;
    vector< vector<RootPTR> > v;
  } RootBoundary;
  */

  static bool polishing;

  template <class N>
  PV2<N> polish(PV2<Parameter> p);

  template <class N>
  void newton(RootBoundary<N>& I);

  template <class N>
  bool subdivide(RootBoundary<N>& I);

  template <class N>
  std::vector<int> intersections(vector<vector<PTR<Object<Scalar>>>>& rbv);

  int parity(const std::vector<int>& alt);

  template <class N>
  void setRoots(RootBoundary<N>& rb);

  template <class N>
  RootBoundary<N> splitHoriz(RootBoundary<N>& rb, N c);

  template <class N>
  RootBoundary<N> splitVert(RootBoundary<N>& rb, N c);

  template <class N>
  bool solve(PV2<N>, PV2<N>, PV2<N>, PV2<N>&);
};

inline Root2PTR::operator PTR<Object<PV2>>() const {
  return PTR<Object<PV2>>(operator Root2*());
}

}  // namespace acp
#endif
