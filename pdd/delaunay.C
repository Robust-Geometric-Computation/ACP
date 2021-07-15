#include "delaunay.h"

Tetrahedrons delaunay (const PPoints &pts, double *bbox, FTMap &ftmap)
{
  Octree<Point *> *octree = Octree<Point *>::octree(pts, bbox, 0.0, 1, 100);
  PPoints spts;
  octree->list(spts);
  delete octree;
  return delaunay(spts, ftmap);
}

Tetrahedrons delaunay (PPoints &pts, FTMap &ftmap)
{
  Tetrahedrons ts;
  delaunayInit(pts, ts, ftmap);
  for (unsigned int i = 4u; i < pts.size(); ++i)
    insertPoint(pts[i], ts, ftmap);
  return ts;
}

void delaunayInit (PPoints &pts, Tetrahedrons &ts, FTMap &ftmap)
{
  if (!firstTetrahedron(pts)) {
    cerr << "coplanar points in Delaunay" << endl;
    exit(0);
  }
  addTetrahedron(ts, ftmap, pts[0], pts[1], pts[2], pts[3]);
  Tetrahedron *t = ts[0];
  for (unsigned int i = 0u; i < 4u; ++i) {
    Facet f = t->facet(i);
    addTetrahedron(ts, ftmap, f.c, f.b, f.a, 0);
  }
}

bool firstTetrahedron (PPoints &pts)
{
  unsigned int n = pts.size();
  for (unsigned int i = 0u; i + 3u < n; ++i)
    for (unsigned int j = i + 1u; j + 2u < n; ++j)
      for (unsigned int k = j + 1u; k + 1u < n; ++k)
	for (unsigned int l = k + 1u; l < n; ++l) {
	  int s = Orientation(pts[i], pts[j], pts[k], pts[l]);
	  if (s) {
	    swap(pts, 0, i);
	    swap(pts, 1, j);
	    swap(pts, 2, k);
	    swap(pts, 3, l);
	    if (s == -1) 
	      swap(pts, 0, 1);
	    return true;
	  }
	}
  return false;
}

void swap (PPoints &pts, unsigned int i, unsigned int j)
{
  if (i != j) {
    Point *p = pts[i];
    pts[i] = pts[j];
    pts[j] = p;
  }
}

void addTetrahedron (Tetrahedrons &ts, FTMap &ftmap, Point *a, Point *b,
		     Point *c, Point *d)
{
  Tetrahedron *t = new Tetrahedron(a, b, c, d);
  addTetrahedron(ts, ftmap, t);
}

void addTetrahedron (Tetrahedrons &ts, FTMap &ftmap, Tetrahedron *t)
{
  ts.push_back(t);
  for (unsigned int i = 0u; i < 4u; ++i) {
    Facet f = t->facet(i);
    FTMap::iterator iter = ftmap.find(f);
    if (iter == ftmap.end())
      ftmap.insert(FTPair(f, t));
    else
      iter->second = t;
  }
}

void insertPoint (Point *p, Tetrahedrons &ts, FTMap &ftmap)
{
  Tetrahedrons::reverse_iterator r = ts.rbegin();
  while ((*r)->ghost())
    ++r;
  Tetrahedron *t = findTetrahedron(p, *r, ftmap);
  insertPoint(p, t, ts, ftmap);
}

Tetrahedron * findTetrahedron (Point *p, Tetrahedron *t, const FTMap &ftmap)
{
  PTR<Point> o = new CentroidPoint(t->vertex(0), t->vertex(1), t->vertex(2),
				   t->vertex(3));
  while (!t->inCircumsphere(p)) {
    Facet f = t->exitFacet(o, p);
    t = getOppositeTetrahedron(f, ftmap);
  }
  return t;
}

Tetrahedron * getOppositeTetrahedron (const Facet &f, const FTMap &ftmap)
{
  Facet g(f.c, f.b, f.a);
  FTMap::const_iterator i = ftmap.find(g);
  return i == ftmap.end() ? 0 : i->second;
}

void insertPoint (Point *p, Tetrahedron *t, Tetrahedrons &ts, FTMap &ftmap)
{
  Tetrahedrons cav = cavity(p, t, ftmap);
  fillCavity(ts, ftmap, p, cav);
}

Tetrahedrons cavity (Point *p, Tetrahedron *t, FTMap &ftmap)
{
  Tetrahedrons st, cav;
  st.push_back(t);
  while (!st.empty()) {
    Tetrahedron *s = st.back();
    st.pop_back();
    if (s->valid()) {
      s->setValid(false);
      cav.push_back(s);
      for (unsigned int i = 0u; i < 4u; ++i) {
	Tetrahedron *u = getOppositeTetrahedron(s->facet(i), ftmap);
	if (u && u->valid() && u->inCircumsphere(p))
	  st.push_back(u);
      }
    }
  }
  return cav;
}

bool inDiametricBall (Point *p, Point *a, Point *b, Point *c)
{
  PTR<Point> o = new CircleCenter(a, b, c);
  return closerPair(p, o, a, o);
}

void fillCavity (Tetrahedrons &ts, FTMap &ftmap, Point *p, const Tetrahedrons &cav)
{
  Facets fs = cavityBoundary(cav, ftmap);
  removeTetrahedrons(cav, ftmap);
  for (Facets::iterator f = fs.begin(); f != fs.end(); ++f)
    addTetrahedron(ts, ftmap, f->a, f->b, f->c, p);
}

Facets cavityBoundary (const Tetrahedrons &cav, const FTMap &ftmap)
{
  Facets fs;
  for (Tetrahedrons::const_iterator t = cav.begin(); t != cav.end(); ++t)
    for (unsigned int i = 0u; i < 4u; ++i) {
      Facet f = (*t)->facet(i);
      Tetrahedron *u = getOppositeTetrahedron(f, ftmap);
      if (!u || u->valid())
	fs.push_back(f);
    }
  return fs;
}

void removeTetrahedrons (const Tetrahedrons &ts, FTMap &ftmap)
{
  for (Tetrahedrons::const_iterator t = ts.begin(); t != ts.end(); ++t)
    removeTetrahedron(*t, ftmap);
}

void removeTetrahedron (Tetrahedron *t, FTMap &ftmap)
{
  for (unsigned int i = 0u; i < 4u; ++i)
    ftmap.erase(t->facet(i));
}

void deleteTetrahedrons (const Tetrahedrons &ts)
{
  for (Tetrahedrons::const_iterator t = ts.begin(); t != ts.end(); ++t)
    delete *t;
}

int intersectTriangleLine (TrianglePlane *p, Point *t, Point *h)
{
  PTR<Point> a = new TriangleLinePoint(p, t, h);
  int pc = p->getPC(), s1 = LeftTurn(p->getA(), p->getB(), a, pc);
  if (s1 == -1)
    return -1;
  int s2 = LeftTurn(p->getB(), p->getC(), a, pc);
  if (s2 == -1)
    return s2;
  int s3 = LeftTurn(p->getC(), p->getA(), a, pc);
  return min(s1, min(s2, s3));
}

Polyhedron * polyhedron (const Tetrahedrons &ts)
{
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (Tetrahedrons::const_iterator t = ts.begin(); t != ts.end(); ++t)
    if ((*t)->valid() && !(*t)->ghost())
      for (unsigned int i = 0u; i < 4u; ++i) {
	Facet f = (*t)->facet(i);
	Vertex *u = a->getVertex(f.a, pvmap), *v = a->getVertex(f.b, pvmap),
	  *w = a->getVertex(f.c, pvmap);
	HEdges ue = u->outgoingHEdges();
	bool flag = false;
	for (HEdges::iterator h = ue.begin(); !flag && h != ue.end(); ++h)
	  flag = (*h)->head() == w && (*h)->getNext()->head() == v;
	if (!flag)
	  a->addTriangle(u, v, w);
      }
  return a;
}

bool tetrahedralMesh (Polyhedron *a)
{
  if (a->cells.empty())
    a->formCells();
  vector<int> bad;
  for (unsigned int i = 1u; i < a->cells.size(); ++i) {
    Cell *c = a->cells[i];
    if (!(c->nShells() == 1 && c->getShell(0)->euler() == 2 &&
	  c->getShell(0)->getHFaces().size() == 4))
      bad.push_back(i);
  }
  if (bad.empty())
    return true;
  cerr << "bad cells: ";
  for (vector<int>::iterator i = bad.begin(); i != bad.end(); ++i)
    cerr << *i << " ";
  cerr << endl;
  return false;
}

// debug

int find (Tetrahedron *t, const Tetrahedrons &ts)
{
  for (unsigned int i = 0u; i < ts.size(); ++i)
    if (ts[i] == t)
      return i;
  return -1;
}

bool check (const Tetrahedrons &ts, const FTMap &ftmap)
{
  bool res = true;
  for (unsigned int i = 0u; i < ts.size(); ++i) {
    Tetrahedron *t = ts[i];
    if (t->valid() && !t->ghost())
      for (unsigned int j = 0u; j < 4u; ++j)
	if (!getOppositeTetrahedron(t->facet(j), ftmap)) {
	  res = false;
	  cerr << "missing tetrahedron " << i << " opposite facet " << j << endl;
	}
  }
  return res;
}

string outi (int i);

void pfacets (const Facets &fs, int i)
{
  vector<PV3<double>> pts;
  vector<unsigned int> idata;
  int j = 0;
  for (Facets::const_iterator f = fs.begin(); f != fs.end(); ++f)
    if (!f->ghost()) {
      pts.push_back(f->a->getApproxMid());
      pts.push_back(f->b->getApproxMid());
      pts.push_back(f->c->getApproxMid());
      idata.push_back(3);
      idata.push_back(3*j);
      idata.push_back(3*j+1);
      idata.push_back(3*j+2);
      ++j;
    }
  ofstream ostr(outi(i).c_str());
  outputVTK(pts, idata, true, ostr);  
}

void pfacet (const Facet &f, int i)
{
  Facets fs;
  fs.push_back(f);
  pfacets(fs, i);
}

void pfacets (const FacetSet &fs, int i)
{
  Facets fa(fs.begin(), fs.end());
  pfacets(fa, i);
}

void ptetra (const Tetrahedrons &ts, int i, bool all, int s, int e)
{
  Facets fs;
  for (unsigned int i = s; i < e; ++i)
    if (ts[i] && (all || ts[i]->valid()))
      for (unsigned int j = 0u; j < 4u; ++j) {
	Facet f = ts[i]->facet(j);
	if (!f.ghost())
	  fs.push_back(f);
      }
  pfacets(fs, i);
}

void ptetra (Tetrahedrons &ts, int i)
{
  ptetra(ts, i, false, 0, ts.size());
}

void ptetra1 (Tetrahedrons &ts, int i)
{
  ptetra(ts, i, true, 0, ts.size());
}

void ptetra (Tetrahedron *te, int i)
{
  Tetrahedrons ts;
  ts.push_back(te);
  ptetra(ts, i);
}

void pps (const Points &pts);

void pf (const Facet &f)
{
  Points pts;
  pts.push_back(f.a);
  pts.push_back(f.b);
  pts.push_back(f.c);
  pps(pts);
}

void pfs (const Facets &fs)
{
  cerr << "(";
  for (Facets::const_iterator f = fs.begin(); f != fs.end(); ++f)
    pf(*f);
  cerr << ")" << endl;
}

