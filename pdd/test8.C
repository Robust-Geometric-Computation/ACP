#include "cspace.h"
#include "io.h"

int main (int argc, char *argv[])
{
  if (argc < 3)
    return 0;
  ifstream astr(argv[1]), bstr(argv[2]);
  if (!(astr.good() && bstr.good()))
    return 0;
  acp::enable();
  double t0 = getTime();
  PTR<Polyhedron> a0 = readPolyhedronVTK(astr), b = readPolyhedronVTK(bstr);
  PTR<Point> t = new Point(-0.5*(a0->bbox[0] + a0->bbox[1]),
			   -0.5*(a0->bbox[2] + a0->bbox[3]),
			   b->bbox[5] - a0->bbox[5] - 0.1);
  PTR<Polyhedron> a = a0->translate(t);
  PTR<Cspace> c = cspace(a, b);
  double dd = 0.01, ds = 1e-15*bboxSize(c->bbox);
  PTR<Polyhedron> p0 = discretize(c, dd), p = p0->selfUnion();
  t0 = getTime() - t0;
  cerr << "faces: " << p->faces.size() << "; cpu time: " << t0 << endl;
  acp::primitiveReport();
  acp::disable();
}

