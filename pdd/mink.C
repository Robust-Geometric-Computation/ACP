#include "mink.h"

Polyhedron * minkowskiSum (Polyhedron *a, Polyhedron *b, bool check)
{
  Polyhedron *con = convolution(a, b, check),
    *c = subdivide(con, true), *d = new Polyhedron;
  Shells sh;
  c->formShells(sh);
  Shells msh = minkowskiShells(a, b, sh);
  Faces fa;
  for (Shells::iterator s = msh.begin(); s != msh.end(); ++s) {
    const HFaces &hf = (*s)->getHFaces();
    for (HFaces::const_iterator f = hf.begin(); f != hf.end(); ++f)
      fa.push_back((*f)->getF());
  }
  Triangles tr = triangulate(fa);
  PVMap pvmap;
  for (Triangles::iterator t = tr.begin(); t != tr.end(); ++t)
    d->addTriangle(t->a, t->b, t->c, pvmap);
  delete c;
  delete con;
  deleteShells(sh);
  return d;
}

Shells minkowskiShells (Polyhedron *a, Polyhedron *b, const Shells &sh)
{
  Shells psh;
  for (Shells::const_iterator s = sh.begin(); s != sh.end(); ++s)
    if ((*s)->pos())
      psh.push_back(*s);
  if (psh.size() == 1)
    return psh;
  return a->faces.size() < b->faces.size() ? minkowskiShellsAux(a, b, psh)
    : minkowskiShellsAux(b, a, psh);
}

Shells minkowskiShellsAux (Polyhedron *a, Polyhedron *b, const Shells &sh)
{
  Octree<Face *> *octree = b->faceOctree();
  const unsigned int n = nthreads(8);
  unsigned int k = sh.size(), m = k/n, is = 0;
  MSData msd[n];
  for (int i = 0; i < n; ++i) {
    msd[i].i = i;
    msd[i].is = is;
    is = i + 1 == n ? k : is + m;
    msd[i].ie = is;
    msd[i].a = a;
    msd[i].octree = octree;
    msd[i].sh = &sh;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) minkowskiShellsT, (void *) (msd + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  Shells res;
  for (int i = 0; i < n; ++i)
    res.insert(res.end(), msd[i].res.begin(), msd[i].res.end());
  delete octree;
  return res;
}

void minkowskiShellsT (void *ptr)
{
  MSData *msd = (MSData *) ptr;
  BaseObject::addThread(msd->i);
  for (int i = msd->is; i < msd->ie; ++i) {
    Shell *s = msd->sh->at(i);
    if (minkowskiShell(msd->a, msd->octree, s))
      msd->res.push_back(s);
  }
}

bool minkowskiShell (Polyhedron *a, Octree<Face *> *octree, Shell *s)
{
  const HFaces &hf = s->getHFaces();
  Point *p = s->getHFaces()[0]->getF()->getBoundary()->tail()->getP();
  Polyhedron *c = a->negativeTranslate(p);
  bool res = !c->intersectsEdges(octree, true);
  delete c;
  return res;
}

void sphereBBox (Edge *e, double *bbox)
{
  PTR<Point> en1 = new FNormal(e->getHEdge(0)->getF()),
    en2 = new FNormal(e->getHEdge(1)->getF());
  sphereBBox(en1, en2, bbox);
}

void sphereBBox (Point *t, Point *h, double *bbox)
{
  PTR<Point> n = new GreatCircleN(t, h), px = new GreatCircleMinX(n),
    py = new GreatCircleMinY(n), pz = new GreatCircleMinZ(n);
  int tx = TripleProduct(n, t, px), hx = TripleProduct(n, px, h),
    ty = TripleProduct(n, t, py), hy = TripleProduct(n, py, h),
    tz = TripleProduct(n, t, pz), hz = TripleProduct(n, pz, h);
  PTR<Point> tu = new UnitVector(t), hu = new UnitVector(h);
  PV3<Interval> tp = tu->getApprox(1.0), hp = hu->getApprox(1.0);
  for (unsigned int i = 0u; i < 3u; ++i) {
    bbox[2*i] = min(tp[i].l, hp[i].l);
    bbox[2*i+1] = max(tp[i].u, hp[i].u);
  }
  if (tx == 1 && hx == 1) {
    UnitVector u(px);
    bbox[0] = min(bbox[0], u.getApprox(1.0).x.l);
  }
  else if (tx == -1 && hx == -1) {
    UnitVector u(px);
    bbox[1] = max(bbox[1], - u.getApprox(1.0).x.l);
  }
  if (ty == 1 && hy == 1) {
    UnitVector u(py);
    bbox[2] = min(bbox[2], u.getApprox(1.0).y.l);
  }
  else if (ty == -1 && hy == -1) {
    UnitVector u(py);
    bbox[3] = max(bbox[3], - u.getApprox(1.0).y.l);
  }
  if (tz == 1 && hz == 1) {
    UnitVector u(pz);
    bbox[4] = min(bbox[4], u.getApprox(1.0).z.l);
  }
  else if (tz == -1 && hz == -1) {
    UnitVector u(pz);
    bbox[5] = max(bbox[5], - u.getApprox(1.0).z.l);
  }
}

void sphereBBox (const HEdges &ed, double *bbox)
{
  bbox[0] = bbox[2] = bbox[4] = 1.0;
  bbox[1] = bbox[3] = bbox[5] = -1.0;
  bool flags[6] = {false, false, false, false, false, false};
  int n = ed.size();
  Points pts;
  for (int i = 0; i < n; ++i)
    pts.push_back(new EENormal(ed[i], ed[(i+1)%n]));
  for (int i = 0; i < n; ++i) {
    Point *t = pts[i], *h = pts[(i+1)%n];
    double eb[6];
    sphereBBox(t, h, eb);
    mergeBBox(eb, bbox); 
    PTR<Point> n = new GreatCircleN(t, h);
    for (int j = 0; j < 3; ++j) {
      int s = coordinate(n, j);
      if (s == 0) {
	bbox[2*j] = -1.0;
	bbox[2*j+1] = 1.0;
      }
      else if (s == -1)
	flags[2*j] = true;
      else
	flags[2*j+1] = true;
    }
  }
  for (int i = 0; i < 3; ++i) {
    if (flags[2*i] && !flags[2*i+1])
      bbox[2*i] = -1.0;
    if (!flags[2*i] && flags[2*i+1])
      bbox[2*i+1] = 1.0;
  }
}

int coordinate (Point *a, int c)
{
  PV3<Interval> ap = a->getApprox(1.0);
  if (ap[c].l > 0.0)
    return 1;
  if (ap[c].u < 0.0)
    return -1;
  return 0;
}

BSPElt::BSPElt (Vertex *v, const HEdges &ed, bool convex) : l(0u), ed(ed), convex(convex)
{
  d.v = v;
  if (ed.empty()) {
    bbox[0] = bbox[2] = bbox[4] = -1.0;
    bbox[1] = bbox[3] = bbox[5] = 1.0;
  }
  else
    sphereBBox(ed, bbox);
  setPbb();
}

BSPElt::BSPElt (Edge *e, bool convex) : l(0), convex(convex)
{
  d.e = e;
  sphereBBox(e, bbox);
  setPbb();
}

BSPElt::BSPElt (Face *f) : l(0), convex(true)
{
  d.f = f;
  PTR<Point> n = new FNormal(f);
  UnitVector u(n);
  u.getBBox(bbox);
  setPbb();
}

BSPElt::BSPElt (const BSPElt &x)
{
  d.v = x.d.v;
  ed = x.ed;
  l = x.l;
  copyBBox(x.bbox, bbox);
  pbb = x.pbb;
  convex = x.convex;
}

void BSPElt::setPbb ()
{
  Parameter xl(bbox[0]), xu(bbox[1]), yl(bbox[2]), yu(bbox[3]),
    zl(bbox[4]), zu(bbox[5]);
  pbb.x = Parameter(xl, xu, false);
  pbb.y = Parameter(yl, yu, false);
  pbb.z = Parameter(zl, zu, false);
}

int BSPElt::side (Point *r) const
{
  PV3<Interval> rr = r->getApprox(1.0);
  PV3<Parameter> rp(Parameter(Parameter(rr.x.l), Parameter(rr.x.u), false),
		    Parameter(Parameter(rr.y.l), Parameter(rr.y.u), false),
		    Parameter(Parameter(rr.z.l), Parameter(rr.z.u), false));
  Parameter k = rp.dot(pbb);
  return k.sign(false);
}

void BSPTree (BSPElts &aelts, BSPElts &belts, BSPElts &ea, BSPElts &eb, 
	      int nmax, int dmax, unsigned int c)
{
  if (dmax == 0 || aelts.size() < nmax && belts.size() < nmax)
    BSPLeaf(aelts, belts, ea, eb);
  else {
    Point r(randomNumber(-1.0, 1.0), randomNumber(-1.0, 1.0), randomNumber(-1.0, 1.0));
    BSPElts aelts1, aelts2, belts1, belts2;
    BSPPartition(aelts, &r, c, aelts1, aelts2);
    BSPPartition(belts, &r, c, belts1, belts2);
    if (aelts1.size() == aelts.size() || aelts2.size() == aelts.size() ||
	belts1.size() == belts.size() || belts2.size() == belts.size())
      BSPLeaf(aelts, belts, ea, eb);
    else {
      BSPTree(aelts1, belts1, ea, eb, nmax, dmax - 1, c + 1u);
      BSPTree(aelts2, belts2, ea, eb, nmax, dmax - 1, c + 1u);
    }
  }
}

void BSPPartition (BSPElts &elts, Point *r, unsigned int c, BSPElts &elts1, BSPElts &elts2)
{
  for (BSPElts::iterator e = elts.begin(); e != elts.end(); ++e) {
    int s = e->side(r);
    if (s == -1)
      elts1.push_back(*e);
    else if (s == 1)
      elts2.push_back(*e);
    else {
      elts1.push_back(*e);
      BSPElt f(*e);
      f.l += 1u << c;
      elts2.push_back(f);
    }
  }
}

void BSPLeaf (const BSPElts &aelts, const BSPElts &belts, BSPElts &ea, 
	      BSPElts &eb)
{
  for (BSPElts::const_iterator e = aelts.begin(); e != aelts.end(); ++e)    
    for (BSPElts::const_iterator f = belts.begin(); f != belts.end(); ++f)
      if (e->compatible(*f) && bboxOverlap(e->bbox, f->bbox)) {
	ea.push_back(*e);
	eb.push_back(*f);
      }
}

bool MinkHullFace::conflict (HEdge *f) const
{
  int s = circulationEEE(f, e, next->e);
  if (s != 0) return s == 1;
  Point *a = e->tail()->getP(), *b = e->head()->getP(),
    *c = next->e->head()->getP(), *d = f->head()->getP();
  return d->onLine(a, b) || d->onLine(a, c) ||
    DegenerateConflict1(a, b, c, d) == 1 || DegenerateConflict2(a, b, c, d) == 1;
}

void MinkHullFace::updateCset (MinkHullFace *h, HEdge *f)
{
  for (HEdges::iterator g = h->cset.begin(); g != h->cset.end(); ++g)
    if (*g != f && conflict(*g))
      cset.push_back(*g);
}

void MinkHullFace::cone (HEdges &hedges) const
{
  const MinkHullFace *h = this;
  do {
    hedges.push_back(h->e);
    h = h->next;
  }
  while (h != this);
}

bool convexCone (Vertex *v, HEdges &hedges)
{
  HEdges he;
  for (int i = 0; i < v->EdgesN(); ++i) {
    HEdge *e = v->getEdge(i)->getHEdge(0);
    he.push_back(e->tail() == v ? e : e->ccw());
  }
  MinkHullFace *hull = initHull(he);
  if (!hull)
    return false;
  for (int i = 3; i < he.size(); ++i) {
    bool flag;
    hull = updateHull(hull, he[i], flag);
    if (!flag) {
      deleteHull(hull);
      return false;
    }
  }
  hull->cone(hedges);
  deleteHull(hull);
  return convexOrder(hedges);
}

MinkHullFace * initHull (HEdges &hedges)
{
  HEdge *e = hedges[0], *f = 0, *g = 0;
  for (int i = 2; i < hedges.size(); ++i) {
    int s = circulationEEE(e, hedges[1], hedges[i]);
    if (s) {
      if (s == -1) {
	f = hedges[1];
	g = hedges[i];
	hedges[i] = hedges[2];
	hedges[2] = g;
      }
      else {
	f = hedges[i];
	g = hedges[1];
	hedges[i] = hedges[2];
	hedges[2] = f;
      }
      break;
    }
  }
  if (!f)
    return 0;
  MinkHullFace *he = new MinkHullFace(e, 0, 0), *hf = new MinkHullFace(f, he, 0), 
    *hg = new MinkHullFace(g, hf, he);
  he->next = hf;
  he->prev = hg;
  hf->next = hg;
  MinkHullFace *faces[3] = {he, hf, hg};
  for (int i = 0; i < 3; ++i)
    for (int j = 3; j < hedges.size(); ++j)
      if (faces[i]->conflict(hedges[j]))
	faces[i]->cset.push_back(hedges[j]);
  return he;
}

int circulationEEE (HEdge *e, HEdge *f, HEdge *g)
{
  return Orientation(e->tail()->getP(), e->head()->getP(), 
		     g->head()->getP(), f->head()->getP());
}

MinkHullFace * updateHull (MinkHullFace *hull, HEdge *e, bool &flag)
{
  flag = true;
  MinkHullFace *fs = hull;
  if (fs->inCset(e)) {
    while (fs->prev != hull && fs->prev->inCset(e))
      fs = fs->prev;
    if (fs->prev == hull) {
      flag = false;
      return hull;
    }
  }
  else {
    do
      fs = fs->next;
    while (fs != hull && !fs->inCset(e));
    if (fs == hull)
      return hull;
  }
  MinkHullFace *fe = fs;
  while (fe->next->inCset(e))
    fe = fe->next;
  return updateHullAux(fs, fe, e);
}

MinkHullFace * updateHullAux (MinkHullFace *fs, MinkHullFace *fe, HEdge *e)
{
  MinkHullFace *f1 = new MinkHullFace(fs->e, fs->prev, 0),
    *f2 = new MinkHullFace(e, f1, fe->next);
  f1->next = f2;
  f1->updateCset(fs->prev, e);
  f1->updateCset(fs, e);
  f2->updateCset(fe, e);
  f2->updateCset(fe->next, e);
  f1->prev->next = f1;
  f2->next->prev = f2;
  fe->next = 0;
  while (fs) {
    MinkHullFace *ptr = fs->next;
    delete fs;
    fs = ptr;
  }
  return f1;
}

void deleteHull (MinkHullFace *hull)
{
  hull->prev->next = 0;
  while (hull) {
    MinkHullFace *nhull = hull->next;
    delete hull;
    hull = nhull;
  }
}

bool convexOrder (const HEdges &hedges)
{
  HEdge *e1 = hedges[0], *e2 = hedges[1], *e3 = hedges[2];
  do {
    e1 = e1->getNext()->getNext()->ccw();
    if (e1 == e2) 
      return true;
    if (e1 == e3)
      return false;
  }
  while (e1 != hedges[0]);
  return false;
}

Polyhedron * convolution (Polyhedron *a, Polyhedron *b, bool check)
{
  SFs sfs;
  PPPMap pmap;
  sumVF(a, b, true, pmap, sfs);
  sumVF(b, a, false, pmap, sfs);
  sumEE(a, b, pmap, sfs);
  return convolution(sfs, check);
}

void sumVF (Polyhedron *a, Polyhedron *b, bool avflag, PPPMap &pmap, SFs &sfs)
{
  BSPElts aelts, belts, ea, eb;
  convexVertices(a, aelts);
  for (Faces::iterator f = b->faces.begin(); f != b->faces.end(); ++f)
    belts.push_back(BSPElt(*f));
  BSPTree(aelts, belts, ea, eb);
  for (int i = 0; i < ea.size(); ++i)
    if (compatibleVF(ea[i].ed, eb[i].d.f))
      sumVF(ea[i].d.v, eb[i].d.f, avflag, pmap, sfs);
}

void convexVertices (Polyhedron *a, BSPElts &elts)
{
  for (Vertices::iterator v = a->vertices.begin(); v != a->vertices.end(); ++v) {
    HEdges ed;
    if (convexCone(*v, ed))
      elts.push_back(BSPElt(*v, ed));
  }
}

bool compatibleVF (HEdges &ed, Face *f)
{
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    if (InnerProductEF(*e, f) == 1)
      return false;
  return true;
}

void sumVF (Vertex *v, Face *f, bool avflag, PPPMap &pmap, SFs &sfs)
{
  Points pf = f->getBoundary()->pointLoop();
  sfs.push_back(SF(sumPP(v->getP(), pf[0], avflag, pmap),
		   sumPP(v->getP(), pf[1], avflag, pmap),
		   sumPP(v->getP(), pf[2], avflag, pmap)));
}

Point * sumPP (Point *a, Point *b, bool aflag, PPPMap &pmap)
{
  Point *aa = aflag ? a : b, *bb = aflag ? b : a;
  pair<Point *, Point *> pp(aa, bb);
  PPPMap::iterator iter = pmap.find(pp);
  if (iter != pmap.end())
    return iter->second;
  Point *p = new SumPoint(aa, bb);
  pmap.insert(pair<pair<Point *, Point *>, Point *>(pp, p));
  return p;
}

void sumEE (Polyhedron *a, Polyhedron *b, PPPMap &pmap, SFs &sfs)
{
  BSPElts aelts, belts, ea, eb;
  convexEdges(a, aelts);
  convexEdges(b, belts);
  BSPTree(aelts, belts, ea, eb);
  for (int i = 0; i < ea.size(); ++i) {
    Edge *eai = ea[i].d.e, *ebi = eb[i].d.e;
    bool aflag;
    if (compatibleEE(eai, ebi, aflag))
      sumEE(eai, ebi, aflag, pmap, sfs);
  }
}

void convexEdges (Polyhedron *a, BSPElts &elts)
{
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if (convexEdge(*e) == 1)
      elts.push_back(BSPElt(*e));
}

int convexEdge (Edge *e)
{
  HEdge *e1 = e->getHEdge(0), *e2 = e->getHEdge(1);
  if (e1->getF()->getP() == e2->getF()->getP())
    return 0;
  return ConvexEdge(e1, e2);
}

bool compatibleEE (Edge *e, Edge *f, bool &aflag)
{
  HEdge *e1 = e->getHEdge(0)->getForward() ? e->getHEdge(0) : e->getHEdge(1),
    *f1 = f->getHEdge(0)->getForward() ? f->getHEdge(0) : f->getHEdge(1),
    *e2 = e1->ccw(), *f2 = f1->ccw();
  int s1 = InnerProductEF(f1, e1->getF()), s2 = InnerProductEF(f1, e2->getF());
  aflag = s1 == 1 || s1 == 0 && s2 == -1;
  if (s1 == s2)
    return false;
  int s3 = InnerProductEF(e1, f1->getF());
  if (s1 != 0 && s1 == s3 || s2 != 0 && s2 == - s3)
    return false;
  int s4 = InnerProductEF(e1, f2->getF());
  if (s3 == s4)
    return false;
  return s4 == 0 || (s1 == 0 ? s2 == - s4 : s1 == s4);
}

void sumEE (Edge *e, Edge *f, bool aflag, PPPMap &pmap, SFs &sfs)
{
  Edge *ee = aflag ? e : f, *ff = aflag ? f : e;
  sfs.push_back(SF(sumPP(ee->getT()->getP(), ff->getT()->getP(), aflag, pmap),
		   sumPP(ee->getH()->getP(), ff->getT()->getP(), aflag, pmap),
		   sumPP(ee->getH()->getP(), ff->getH()->getP(), aflag, pmap),
		   sumPP(ee->getT()->getP(), ff->getH()->getP(), aflag, pmap)));
}

Polyhedron * convolution (const SFs &sfs, bool check)
{
  Polyhedron *c = new Polyhedron;
  PVMap pvmap = check ? pvMap(sfs, c) : PVMap();
  for (SFs::const_iterator s = sfs.begin(); s != sfs.end(); ++s) {
    Vertices ve = vertices(*s, c, pvmap);
    if (ve.size() == 3 && (!check || !faceVertices(ve[0], ve[1], ve[2])))
      c->addTriangle(ve[0], ve[1], ve[2]);
    else if (ve.size() == 4 && (!check || !faceVertices(ve[0], ve[1], ve[2], ve[3]) &&
				!ve[0]->getP()->onLine(ve[1]->getP(), ve[2]->getP())))
      c->addRectangle(ve[0], ve[1], ve[2], ve[3]);
  }
  return c;
}

PVMap pvMap (const SFs &sfs, Polyhedron *c)
{
  Points pts;
  for (SFs::const_iterator s = sfs.begin(); s != sfs.end(); ++s) {
    pts.push_back(s->a);
    pts.push_back(s->b);
    pts.push_back(s->c);
    if (s->d)
      pts.push_back(s->d);
  }
  return pvMap(pts, c);
}

Vertices vertices (const SF &s, Polyhedron *a, PVMap &pvmap)
{
  Vertices ve;
  ve.push_back(a->getVertex(s.a, pvmap));
  ve.push_back(a->getVertex(s.b, pvmap));
  ve.push_back(a->getVertex(s.c, pvmap));
  if (s.d)
    ve.push_back(a->getVertex(s.d, pvmap));
  return ve;
}

