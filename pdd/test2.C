#include "cdt.h"

int main (int argc, char *argv[])
{
  if (argc < 2)
    return 0;
  ifstream str(argv[1]);
  if (!str.good())
    return 0;
  acp::enable();
  double t = getTime();
  PTR<Polyhedron> a = readPolyhedronVTK(str, false), b = cdt(a);
  tetrahedralMesh(b);
  t = getTime() - t;
  cerr << "faces " << b->faces.size() << " cpu time " << t << endl;
  acp::primitiveReport();
  acp::disable();
}

