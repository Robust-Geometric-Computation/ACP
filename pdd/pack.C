#include "pack.h"

double pack3 (Polyhedron *p1, Polyhedron *p2, Polyhedron *p3, bool check,
	      double minsep, PTR<Point> t[3])
{
  Polyhedron *p2m = p2->negative(), *p3m = p3->negative(),
    *u12 = minkowskiSum(p1, p2m, check), *u13 = minkowskiSum(p1, p3m, check),
    *u23 = minkowskiSum(p2, p3m, check);
  delete p2m;
  delete p3m;
  double d[] = {bboxSize(p1->bbox), bboxSize(p2->bbox), bboxSize(p3->bbox)};
  sort(d, d + 3);
  double lb = d[2], ub = d[1] + d[2] + 0.1, tt = 0.0;
  int k = 0;
  //cerr << "initial interval: [" << lb << ", " << ub << "]" << endl;
  t[0] = t[1] = t[2] = 0;
  while (ub - lb > minsep) {
    ++k;
    double s = 0.5*(lb + ub), t0 = getTime();
    bool flag = pack3(p1, p2, p3, t, s, u12, u13, u23);
    double t1 = getTime() - t0;
    //cerr << setprecision(17);
    //cerr << k << (flag ? ". yes " : ". no ") << s << "; time: " << t1 << endl;
    tt += t1;
    if (flag)
      ub = s;
    else
      lb = s;
  }
  cerr << "minimum size: " << ub << "; cpu time: " << tt << endl;
  delete u12;
  delete u13;
  delete u23;
  return ub;
}

bool pack3 (Polyhedron *p1, Polyhedron *p2, Polyhedron *p3, PTR<Point> t[3],
	    double s, Polyhedron *u12, Polyhedron *u13, Polyhedron *u23)
{
  bool fail;
  if (pack3(p1, p2, p3, t, s, true, fail, u12, u13, u23))
    return true;
  if (!fail)
    return false;
  cerr << "heuristic failed" << endl;
  return pack3(p1, p2, p3, t, s, false, fail, u12, u13, u23);
}

bool pack3 (Polyhedron *p1, Polyhedron *p2, Polyhedron *p3, PTR<Point> t[3],
	    double s, bool flag, bool &fail,
	    Polyhedron *u12, Polyhedron *u13, Polyhedron *u23)
{
  fail = false;
  double u01[6], u02[6], u03[6], v012[6], v013[6], v023[6],
    bb[6] = {0.0, s, 0.0, s, 0.0, s};
  freeSpace(p1, bb, u01);
  freeSpace(p2, bb, u02);
  freeSpace(p3, bb, u03);
  freeSpace(u01, u02, v012);
  freeSpace(u01, u03, v013);
  freeSpace(u02, u03, v023);
  Polyhedron *b012 = flag ? expandedBox(v012) : box(v012, false),
    *b013 = box(v013, false),
    *b023 = flag ? expandedBox(v023) : box(v023, false),
    *v12 = coalesceD(b012->boolean(u12, Complement)),
    *v13 = coalesceD(b013->boolean(u13, Complement)),
    *v23 = coalesceD(b023->boolean(u23, Complement)),
    *v123 = coalesceD(minkowskiSum(v12, v23, true)),
    *w13 = v13->boolean(v123, Intersection);
  bool res = w13->faces.size() > 0;
  if (res) {
    PTR<Point> t13 = interiorPoint(w13);
    if (flag && badv123(v012, u12, v12, v023, u23, v23, t13))
      fail = true;
    else {
      PTR<Point> t13m = new NegPoint(t13);
      Polyhedron *u01b = box(u01, false), *u02b = box(u02, false),
	*u03b = box(u03, false), *u03t = u03b->translate(t13m),
	*v01 = coalesceD(u01b->boolean(u03t, Intersection)),
	*v01m = v01->negative(), 
	*v01u02 = coalesceD(minkowskiSum(u02b, v01m, true)),
	*v23t = v23->negativeTranslate(t13),
	*temp = coalesceD(v23t->boolean(v01u02, Intersection)),
	*w12 = coalesceD(v12->boolean(temp, Intersection));
      if (w12->faces.empty())
	res = false;
      else {
	PTR<Point> t12 = interiorPoint(w12), t12m = new NegPoint(t12);
	Polyhedron *u02t = u02b->translate(t12m),
	  *w01 = v01->boolean(u02t, Intersection);
	t[0] = interiorPoint(w01);
	t[1] = new SumPoint(t[0], t12);
	t[2] = new SumPoint(t[0], t13);
	delete u02t; delete w01;
      }
      delete u01b; delete u02b; delete u03b; delete u03t; delete v01; delete v01m;
      delete v01u02; delete v23t; delete temp; delete w12;
    }
  }
  delete b012; delete b013; delete b023; delete v12; delete v13; delete v23;
  delete v123; delete w13;
  return res && !fail;
}

void freeSpace (Polyhedron *a, double *box, double *res)
{
  double *bboxa = a->bbox;
  for (int i = 0; i < 6; ++i)
    res[i] = box[i] - bboxa[i];
}

// a moves with respect to b
void freeSpace (double *a, double *b, double *res)
{
  for (int i = 0; i < 3; ++i) {
    res[2*i] = b[2*i] - a[2*i+1];
    res[2*i+1] = b[2*i+1] - a[2*i];
  }
}

Polyhedron * expandedBox (double *b)
{
  double e[6], s = BaseObject::delta*bboxSize(b);
  for (int i = 0; i < 3; ++i) {
    e[2*i] = b[2*i] - s;
    e[2*i+1] = b[2*i+1] + s;
  }
  return box(e, true);
}

Polyhedron * coalesceD (Polyhedron *a)
{
  Polyhedron *b = coalesce(a);
  delete a;
  return b;
}

PTR<Point> interiorPoint (Polyhedron *a)
{
  a->computeWindingNumbers();
  for (int i = 1; i < a->cells.size(); ++i)
    if (a->cells[i]->getWN() == 1)
      return a->cells[i]->interiorPoint();
  return 0;
}

bool badv123 (double *v012, Polyhedron *u12, Polyhedron *&v12,
	      double *v023, Polyhedron *u23, Polyhedron *&v23, Point *t13)
{
  Polyhedron *b012 = box(v012, false), *b023 = box(v023, false);
  delete v12;
  v12 = b012->boolean(u12, Complement);
  delete v23;
  v23 = b023->boolean(u23, Complement);
  Polyhedron *v12m = v12->negativeTranslate(t13);
  bool res = !v12m->intersects(v23, true, 0);
  delete b012;
  delete b023;
  delete v12m;
  return res;
}

void pack3output (Polyhedron *a, Polyhedron *b, Polyhedron *c,
		  PTR<Point> t[3], ostream &ostr)
{
  Polyhedron *p[] = {a, b, c}, *q[3];
  Faces fa;
  for (int i = 0; i < 3; ++i) {
    q[i] = p[i]->translate(t[i]);
    fa.insert(fa.end(), q[i]->faces.begin(), q[i]->faces.end());
  }
  writePolyhedronVTK(fa, ostr);
  for (int i = 0; i < 3; ++i)
    delete q[i];
}

int main (int argc, char *argv[])
{
  if (argc < 4)
    return 0;
  ifstream astr(argv[1]), bstr(argv[2]), cstr(argv[3]);
  if (!(astr.good() && bstr.good() && cstr.good()))
    return 0;
  acp::enable();
  PTR<Polyhedron> a = readPolyhedronVTK(astr, false),
    b = readPolyhedronVTK(bstr, false), c = readPolyhedronVTK(cstr, false);
  double minsep = 1e-6;
  PTR<Point> t[3];
  pack3(a, b, c, true, minsep, t);
  //if (t[0])
  //pack3output(a, b, c, t, cout);
  acp::primitiveReport();
  acp::disable();
}
