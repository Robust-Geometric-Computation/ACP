#include "hull.h"

Polyhedron * convexHull (Polyhedron *a)
{
  Points pts;
  for (Vertices::iterator v = a->vertices.begin(); v != a->vertices.end(); ++v)
    pts.push_back((*v)->getP());
  return convexHull(pts);
}

Polyhedron * convexHull (Points &pts)
{
  int n = pts.size(), *p = randomPermutation(n);
  Points ppts;
  for (int i = 0; i < n; ++i)
    ppts.push_back(pts[p[i]]);
  delete [] p;
  Polyhedron *a = hullInit(ppts);
  if (!a)
    return 0;
  ConflictGraph cg = conflictGraphInit(a);
  for (Vertices::iterator v = a->vertices.begin() + 4; v != a->vertices.end(); ++v)
    expandHull(a, *v, cg);
  a->removeNullFaces();
  return a;
}

int * randomPermutation (int n)
{
  int *p = new int [n];
  for (int i = 0; i < n; ++i)
    p[i] = i;
  for (int i = 0; i < n - 1; ++i) {
    int j = randomInteger(i, n - 1), ti = p[i];
    p[i] = p[j];
    p[j] = ti;
  }
  return p;
}

int randomInteger (int lb, int ub)
{
  double r = random()/(double) RAND_MAX;
  int d = (int) floor(r*(ub - lb + 1));
  return lb + d;
}

Polyhedron * hullInit (Points &pts)
{
  int n = pts.size();
  if (n < 4)
    return 0;
  int i = 2;
  while (i < n && pts[i]->onLine(pts[0], pts[1]))
    ++i;
  if (i == n)
    return 0;
  if (i > 2) {
    PTR<Point> p = pts[2];
    pts[2] = pts[i];
    pts[i] = p;
  }
  i = 3;
  int s;
  while (i < n) {
    s = Orientation(pts[0], pts[1], pts[2], pts[i]);
    if (s)
      break;
    ++i;
  }
  if (i == n)
    return 0;
  if (i > 3) {
    PTR<Point> p = pts[3];
    pts[3] = pts[i];
    pts[i] = p;
  }
  if (s == -1)
    reverse(pts.begin(), pts.begin() + 2);
  Polyhedron *a = new Polyhedron;
  for (Points::iterator p = pts.begin(); p != pts.end(); ++p)
    a->getVertex(*p);
  a->addTriangle(a->vertices[0], a->vertices[2], a->vertices[3]);
  a->addTriangle(a->vertices[2], a->vertices[1], a->vertices[3]);
  a->addTriangle(a->vertices[1], a->vertices[0], a->vertices[3]);
  a->addTriangle(a->vertices[0], a->vertices[1], a->vertices[2]);
  return a;
}

ConflictGraph::~ConflictGraph ()
{
  for (map<Vertex *, set<Face *> *>::iterator i = vcon.begin(); i != vcon.end(); ++i)
    delete i->second;
  for (map<Face *, set<Vertex *> *>::iterator i = fcon.begin(); i != fcon.end(); ++i)
    delete i->second;
}

void ConflictGraph::insert (Vertex *v, Face *f)
{
  map<Vertex *, set<Face *> *>::iterator vi = vcon.find(v);
  if (vi == vcon.end()) {
    set<Face *> *fs = new set<Face *>;
    fs->insert(f);
    vcon.insert(pair<Vertex *, set<Face *> *>(v, fs));
  }
  else
    vi->second->insert(f);
  map<Face *, set<Vertex *> *>::iterator fi = fcon.find(f);
  if (fi == fcon.end()) {
    set<Vertex *> *vs = new set<Vertex *>;
    vs->insert(v);
    fcon.insert(pair<Face *, set<Vertex *> *>(f, vs));
  }
  else
    fi->second->insert(v);
}

void ConflictGraph::erase (Face *f)
{
  set<Vertex *> *vs = fcon.find(f)->second;
  for (set<Vertex *>::iterator v = vs->begin(); v != vs->end(); ++v) 
    vcon.find(*v)->second->erase(f);
}

set<Face *> * ConflictGraph::conflicts (Vertex *v)
{
  map<Vertex *, set<Face *> *>::iterator i = vcon.find(v);
  return i == vcon.end() || i->second->empty() ? 0 : i->second;
}

set<Vertex *> * ConflictGraph::conflicts (Face *f)
{
  map<Face *, set<Vertex *> *>::iterator i = fcon.find(f);
  return i == fcon.end() || i->second->empty() ? 0 : i->second;
}

void ConflictGraph::update (Vertex *v, Face *f, set<Vertex *> *vs)
{
  if (vs)
    for (set<Vertex *>::iterator w = vs->begin(); w != vs->end(); ++w)
      if (*w != v && visible(*w, f))
	insert(*w, f);
}

ConflictGraph conflictGraphInit (Polyhedron *a)
{
  ConflictGraph cg;
  for (Vertices::iterator v = a->vertices.begin() + 4; v != a->vertices.end(); ++v)
    for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f)
      if (visible(*v, *f))
	cg.insert(*v, *f);
  return cg;
}

bool visible (Vertex *v, Face *f)
{
  int s = v->getP()->side(f->getP());
  return s == 1 || s == 0 && !f->contains(v->getP(), false);
}

void expandHull (Polyhedron *a, Vertex *v, ConflictGraph &cg)
{
  set<Face *> *fs = cg.conflicts(v);
  if (!fs)
    return;
  HEdges he = horizon(*fs);
  for (HEdges::iterator h = he.begin(); h != he.end(); ++h)
    expandHull(a, v, cg, *h);
  Faces fa(fs->begin(), fs->end());
  for (Faces::iterator f = fa.begin(); f != fa.end(); ++f) {
    cg.erase(*f);
    a->removeLoop((*f)->getBoundary());
  }
}

HEdges horizon (const set<Face *> &fs)
{
  HEdges he;
  for (set<Face *>::const_iterator f = fs.begin(); f != fs.end(); ++f) {
    HEdges fe = (*f)->getBoundary()->edgeLoop();
    for (HEdges::iterator g = fe.begin(); g != fe.end(); ++g)
      if (fs.find((*g)->ccw()->getF()) == fs.end())
	he.push_back(*g);
  }
  return he;
}

void expandHull (Polyhedron *a, Vertex *v, ConflictGraph &cg, HEdge *h)
{
  set<Vertex *> *vs1 = cg.conflicts(h->getF()), *vs2 = cg.conflicts(h->ccw()->getF());
  Face *f = a->addTriangle(v, h->tail(), h->head());
  cg.update(v, f, vs1);
  cg.update(v, f, vs2);  
}

/*
#include "io.h"
int main (int argc, char *argv[])
{
  if (argc < 2)
    return 0;
  ifstream astr(argv[1]);
  if (!astr.good())
    return 0;
  bool perturbed = argc == 2 ? true : *argv[2] != 't';
  acp::enable();
  PTR<Polyhedron> a = readPolyhedronVTK(astr, perturbed),
    b = convexHull(a);
  b->formCells();
  b->describe();
  writePolyhedronVTK(b->faces, cout);
  acp::disable();
}
*/
