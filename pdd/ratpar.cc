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

double nextDR (double x)
{
  static double p = 1.0 / (1 << 26) / (1 << 26);
  double x2 = x + p * fabs(x);
  if (x2 == x)
    return nextafter(x, 1.0);
  double x3 = 0.5 * (x + x2);
  return x3 > x ? x3 : x2;
}

double prevDR (double x)
{
  static double p = 1.0 / (1 << 26) / (1 << 26);
  double x2 = x - p * fabs(x);
  if (x2 == x)
    return nextafter(x, -1.0);
  double x3 = 0.5 * (x + x2);
  return x3 < x ? x3 : x2;
}

void rationalize (double x, double& u, double& v)
{
  if (x == long(x)) {
    u = x;
    v = 1.0;
  }
  else {
    int exp;
    u = frexp(x, &exp) * (1l << 53);
    exp -= 53;
    while (long(u) % 2lu == 0lu) {
      u *= 0.5;
      ++exp;
    }
    v = exp2(-exp);
  }
}

void pparReport ()
{
}

void resetPpar ()
{
}

unsigned long int PParameter::nModSign = 0u;
unsigned int PParameter::nMixed = 0u;
unsigned int PParameter::nRascal = 0u;
