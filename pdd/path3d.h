#include "geometry3d.h"
#include "io.h"
#include <cstring>
#include <queue>
#include <map>

void findPath(Polyhedron * blockspace, int cell_index, PTR<Point> start, PTR<Point> end, Points &path);
void findPath(Polyhedron * blockspace, PTR<Point> start, PTR<Point> end, Points &path);

void testBoundaryCondition(Polyhedron * poly, int mode);