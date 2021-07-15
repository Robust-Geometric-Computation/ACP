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

unsigned int MParameter::nMixed = 0u;
unsigned int MParameter::nRascal = 0u;

void mparReport ()
{
  if (MParameter::nMixed > 0)
    cerr << " mixed hom " << MParameter::nMixed;
  if (MParameter::nRascal > 0)
    cerr << " rascal hom " << MParameter::nRascal;
}

void resetMpar ()
{
  MParameter::nMixed = 0u;
  MParameter::nRascal = 0u;
}
