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

#include "acp.h"
using namespace acp;

namespace acp {

pthread_key_t idkey;
unsigned int threadIds[TMAX], curPrecisions[TMAX] = {53u};

unsigned int threadId ()
{
  void* i = pthread_getspecific(idkey);
  return i ? *(unsigned int*)i : 0;
}

unsigned int & curPrecision ()
{
  return curPrecisions[threadId()];
}

bool enabled = false;

bool isEnabled () { return enabled; }

void enable ()
{
  if (!enabled) {
    enabled = true;
    pthread_key_create(&idkey, NULL);
  }
#ifndef NOMODE
  fesetround(FE_DOWNWARD);
#endif
}

void disable ()
{
#ifndef NOMODE
  fesetround(FE_TONEAREST);
#endif
}

double randomNumber (double rmin, double rmax)
{
  return rmin + (rmax - rmin)*random()/double(RAND_MAX);
}

double getTime ()
{
  timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1e-6*tv.tv_usec;
}

double diffTime (const timeval &t2, const timeval &t1)
{
  return t2.tv_sec - t1.tv_sec + 1e-6*(t2.tv_usec - t1.tv_usec);
}

#ifndef NOMODE
#include "par.cc"
#else
#include "parnm.cc"
#endif
#include "mod.cc"
#ifndef EXACTRAT
#include "ppar.cc"
#else
#include "ratpar.cc"
#endif
#include "mpar.cc"

unsigned int BaseObject::highPrecision = 106u,
  BaseObject::maxPrecision = 424u, BaseObject::homPrecision = 265u;
  int BaseObject::perIs[TMAX] = {-1};
  double BaseObject::delta = pow(2.0, -27), BaseObject::perTs[TMAX],
    BaseObject::maxP = 1e-8;
  unsigned long int Primitive::nPrim = 0u, Primitive::nAmb = 0u,
    Primitive::nDeg = 0u, Primitive::nAlg = 0u, Primitive::nHom = 0u,
    Primitive::nId = 0u, Primitive::maxBits = 0u;
  double Primitive::minPrim = 1e100, Primitive::minPR = 1e100,
    Primitive::tPar = 0.0, Primitive::tPPar = 0.0, Primitive::tMPar = 0.0;
  
void primitiveReport ()
{
  double a = 100.0*double(Primitive::nAmb)/double(Primitive::nPrim),
    d = 100.0*double(Primitive::nDeg)/double(Primitive::nAmb),
    t = Primitive::tPar + Primitive::tPPar + Primitive::tMPar,
    f = t == 0.0 ? 0.0 :100.0*Primitive::tPar/t,
    m = t == 0.0 ? 0.0 : 100.0*Primitive::tPPar/t,
    e = t == 0.0 ? 0.0 : 100.0*Primitive::tMPar/t;
  cout << setprecision(7) << endl << "ACP statistics" << endl
       << "p = " << Primitive::nPrim << " a = " << a << " d = " << d
       << " t = " << t << " f = " << f << " m = " << m << endl
       << "e = " << e;
  pparReport();
  if (Primitive::maxBits)
    cout << " b = " << Primitive::maxBits;
  mparReport();
  if (Primitive::nAlg) {
    double l = 100.0*double(Primitive::nAlg)/double(Primitive::nAmb),
      i = 100.0*double(Primitive::nId)/double(Primitive::nAmb),
      r = log10(Primitive::minPR);
    cout << " l = " << l << " i = " << i << " r = " << r;
  }
  cout << endl;
}

void resetReport ()
{
  resetPpar();
  resetMpar();
  Primitive::nPrim = 0u;
  Primitive::nAmb = 0u;
  Primitive::nDeg = 0u;
  Primitive::maxBits = 0u;
  Primitive::nAlg = 0u;
  Primitive::nHom = 0u;
  Primitive::nId = 0u;
  Primitive::minPrim = 1e100;
  Primitive::minPR = 1e100;
  Primitive::tPar = 0;
  Primitive::tPPar = 0;
  Primitive::tMPar = 0;
}

}  // namespace acp
