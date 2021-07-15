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

unsigned long int PParameter::nModSign = 0u;
unsigned int PParameter::nMixed = 0u;
unsigned int PParameter::nRascal = 0u;

PParameter::PParameter (const MParameter& a)
  : p(a.lb(), a.ub()), mods(a.mods), alg(a.alg) {}

void pparReport ()
{
  if (PParameter::nMixed > 0)
    cerr << " mixed mod " << PParameter::nMixed;
  if (PParameter::nRascal > 0)
    cerr << " rascal mod " << PParameter::nRascal;
}

void resetPpar ()
{
  PParameter::nMixed = 0u;
  PParameter::nRascal = 0u;
  PParameter::nModSign = 0u;
}
