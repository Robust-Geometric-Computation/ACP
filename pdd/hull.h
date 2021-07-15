#ifndef HULL
#define HULL

#include "polyhedron.h"

Polyhedron * convexHull (Polyhedron *a);

Polyhedron * convexHull (Points &pts);

int * randomPermutation (int n);

int randomInteger (int lb, int ub);

Polyhedron * hullInit (Points &pts);

class ConflictGraph {
  map<Vertex *, set<Face *> *> vcon;
  map<Face *, set<Vertex *> *> fcon;
 public:
  ~ConflictGraph ();
  void insert (Vertex *v, Face *f);
  void erase (Face *f);
  set<Face *> * conflicts (Vertex *v);
  set<Vertex *> * conflicts (Face *f);
  void update (Vertex *v, Face *f, set<Vertex *> *vs);
};

ConflictGraph conflictGraphInit (Polyhedron *a);

bool visible (Vertex *v, Face *f);

void expandHull (Polyhedron *a, Vertex *v, ConflictGraph &cg);

HEdges horizon (const set<Face *> &fs);

void expandHull (Polyhedron *a, Vertex *v, ConflictGraph &cg, HEdge *h);

#endif
