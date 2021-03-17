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

#include "acp/core/object.h"

double getTime() {
  timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1e-6 * tv.tv_usec;
}

namespace acp {
bool inInt = false;
unsigned int BaseObject::deltaPrecision = 53u, BaseObject::maxPrecision = 848u;
int BaseObject::perIs[TMAX] = {-1};
double BaseObject::perTs[TMAX], BaseObject::maxP = 1e-8,
                                BaseObject::wMin = 1e100;
bool BaseObject::perInts[TMAX] = {false};
unsigned long int Primitive::nPrim = 0u, Primitive::nAmb = 0u,
                  Primitive::nDeg = 0u, Primitive::nAlg = 0u,
                  Primitive::nHom = 0u, Primitive::nId = 0u,
                  Primitive::maxBits = 0u;
double Primitive::tPar = 0.0, Primitive::tPPar = 0.0, Primitive::tMPar = 0.0,
       Primitive::vMPar = 1e100;

vector<set<void*>> inputSets;

void primitiveReport() {
  unsigned int nm = PParameter::nMixed + MParameter::nMixed,
               nr = PParameter::nRascal + MParameter::nRascal;
  cerr << setprecision(6);
  cerr << "prim " << Primitive::nPrim << " amb " << Primitive::nAmb << " deg "
       << Primitive::nDeg << " unamb mod " << PParameter::nModSign;
  if (nm || nr) cerr << " mixed mod " << nm << " rascal " << nr;
  if (Primitive::maxBits) cerr << " max bits " << Primitive::maxBits;
  cerr << endl;
  if (Primitive::nAlg)
    cerr << "alg " << Primitive::nAlg << " hom " << Primitive::nHom
         << " deg hom " << Primitive::nId << " min hom ratio "
         << BaseObject::wMin << " min nonzero prim " << Primitive::vMPar
         << endl;
  double tPrim = Primitive::tPar + Primitive::tPPar + Primitive::tMPar,
         pp = tPrim == 0.0 ? 0.0 : 100.0 * Primitive::tPar / tPrim;
  cerr << "prim time: par " << Primitive::tPar << " ppar " << Primitive::tPPar
       << " mpar " << Primitive::tMPar << " total " << tPrim << " par % " << pp
       << endl;
}

void resetReport() {
  Primitive::nPrim = 0;
  Primitive::nAmb = 0;
  Primitive::nDeg = 0;
  PParameter::nModSign = 0;
  Primitive::nAlg = 0;
  Primitive::nHom = 0;
  Primitive::tPar = 0;
  Primitive::tPPar = 0;
  Primitive::tMPar = 0;
}
}  // namespace acp
