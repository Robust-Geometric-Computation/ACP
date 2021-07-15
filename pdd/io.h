#ifndef IO
#define IO

#include "polyhedron.h"
#include <fstream>
#include <sstream>

Polyhedron * readPolyhedron (istream &istr, bool perturb = true);

void readCells (istream &istr, Polyhedron *a);

Cell * readCell (istream &istr, Polyhedron *a);

Shell * readShell (istream &istr, Polyhedron *a);

typedef pair<Vertex *, unsigned int> VIPair;

typedef map<Vertex *, unsigned int> VIMap;

typedef pair<Face *, unsigned int> FIPair;

typedef map<Face *, unsigned int> FIMap;

void writePolyhedron (Polyhedron *a, ostream &ostr);

void writeCells (Polyhedron *a, ostream &ostr);

void writeShell (Shell *s, FIMap &fimap, ostream &ostr);

Polyhedron * readPolyhedronVTK (istream &istr, bool perturb = true);

void skipComments (istream &istr);

void writePolyhedronVTK (Polyhedron *a, ostream &ostr);

void writePolyhedronVTK (const Faces &fa, ostream &ostr);

int getPoint (VIMap &vimap, vector<PV3<double>> &pts, Vertex *v);

void outputVTK (const vector<PV3<double>> &pts, const vector<unsigned int> &data,
		int ptype, ostream &ostr);

Polyhedron * readPolyhedronOBJ (istream &istr, bool perturb = true);

bool readPointOBJ (istream &istr, double &x, double &y, double &z);

bool readTriangleOBJ (istream &istr, int n, int &i, int &j, int &k);

void readIndexOBJ (istream &istr, int n, int &i);

void addTriangleOBJ (Polyhedron *a, Vertex *u, Vertex *v, Vertex *w);

#endif
