#include "cdt.h"
#include "mink.h"

int main (int argc, char *argv[])
{
  if (argc < 3)
    return 0;
  ifstream astr(argv[1]), bstr(argv[2]);
  if (!(astr.good() && bstr.good()))
    return 0;
  acp::enable();
  double t = getTime();
  PTR<Polyhedron> a = readPolyhedronVTK(astr, false),
    b = readPolyhedronVTK(bstr, false), c = minkowskiSum(a, b, true),
    d = cdt(c);
  t = getTime() - t;
  cerr << "faces " << d->faces.size() << " cpu time " << t << endl;
  acp::primitiveReport();
  acp::disable();
}

