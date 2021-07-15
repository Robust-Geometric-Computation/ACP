#ifndef PACK
#define PACK

#include "io.h"
#include "mink.h"

double pack3 (Polyhedron *p1, Polyhedron *p2, Polyhedron *p3, bool check,
	      double minsep, PTR<Point> t[3]);

bool pack3 (Polyhedron *p1, Polyhedron *p2, Polyhedron *p3, PTR<Point> t[3],
	    double s, Polyhedron *u12, Polyhedron *u13, Polyhedron *u23);

bool pack3 (Polyhedron *p1, Polyhedron *p2, Polyhedron *p3, PTR<Point> t[3],
	    double s, bool flag, bool &fail,
	    Polyhedron *u12, Polyhedron *u13, Polyhedron *u23);

void freeSpace (Polyhedron *a, double *box, double *res);

void freeSpace (double *a, double *b, double *res);

Polyhedron * expandedBox (double *b);

Polyhedron * coalesceD (Polyhedron *a);

PTR<Point> interiorPoint (Polyhedron *a);

bool badv123 (double *v012, Polyhedron *u12, Polyhedron *&v12,
	      double *v023, Polyhedron *u23, Polyhedron *&v23, Point *t13);

void pack3output (Polyhedron *a, Polyhedron *b, Polyhedron *c,
		  PTR<Point> t[3], ostream &ostr);

#endif
