#include "cdt.h"

Polyhedron * cdt (Polyhedron *a)
{
  PPoints pts;
  for (Vertices::iterator v = a->vertices.begin(); v != a->vertices.end(); ++v)
    pts.push_back((*v)->getP());
  FTMap ftmap;
  Tetrahedrons ts = delaunay(pts, a->bbox, ftmap);
  PTMap ptmap;
  update(ts, 0u, ptmap);
  SFMap sfmap = facets(a);
  Points steiner;
  Octree<Point *> *octree = Octree<Point *>::octree(pts, a->bbox, 0.0, 1, 100);
  recoverSegments(ts, ftmap, ptmap, sfmap, steiner, octree);
  delete octree;
  recoverFacets(ts, ftmap, ptmap, sfmap);
  FacetSet fs;
  for (SFMap::iterator i = sfmap.begin(); i != sfmap.end(); ++i)
    fs.insert(i->second);
  Tetrahedrons ti = innerTetrahedrons(ts, ftmap, fs);
  Polyhedron *b = polyhedron(ti);
  deleteTetrahedrons(ts);
  return b;
}

void update (const Tetrahedrons &ts, unsigned int s, PTMap &ptmap)
{
  for (unsigned int i = s; i < ts.size(); ++i) {
    Tetrahedron *t = ts[i];
    if (t->valid() && !t->ghost())
      for (unsigned int j = 0u; j < 4u; ++j) {
	Point *p = t->vertex(j);
	PTMap::iterator iter = ptmap.find(p);
	if (iter == ptmap.end())
	  ptmap.insert(PTPair(p, t));
	else
	  iter->second = t;
      }
  }
}

SFMap facets (Polyhedron *a)
{
  set<Point*> ac;
  for (Vertices::iterator v = a->vertices.begin(); v != a->vertices.end(); ++v)
    if (acute(*v))
      ac.insert((*v)->getP());
  SFMap sfmap;
  for (unsigned int i = 0u; i < a->faces.size(); ++i) {
    Points pts = a->faces[i]->boundaryPoints();
    bool f1 = ac.find(pts[0]) != ac.end(), f2 = ac.find(pts[1]) != ac.end(),
      f3 = ac.find(pts[2]) != ac.end();
    Segment s1(pts[0], pts[1], true, f1, f2), s2(pts[1], pts[2], true, f2, f3),
      s3(pts[2], pts[0], true, f3, f1);
    Facet fa(pts[0], pts[1], pts[2], i + 1u);
    sfmap.insert(SFPair(s1, fa));
    sfmap.insert(SFPair(s2, fa));
    sfmap.insert(SFPair(s3, fa));
  }
  return sfmap;
}

bool acute (Vertex *a)
{
  HEdges ed = a->outgoingHEdges();
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    if (SmallAngle((*e)->getNext()->head()->getP(), (*e)->tail()->getP(),
		   (*e)->head()->getP(), 0.5) == 1)
      return true;
  return false;
}

void recoverSegments (Tetrahedrons &ts, FTMap &ftmap, PTMap &ptmap,
		      SFMap &sfmap, Points &steiner, Octree<Point *> *octree)
{
  Segments sm = missingSegments(ts, sfmap);
  while (!sm.empty())
    splitSegment(ts, ftmap, ptmap, sfmap, steiner, octree, sm);
}

Segments missingSegments (const Tetrahedrons &ts, const SFMap &sfmap)
{
  SegmentSet ss;
  for (Tetrahedrons::const_iterator t = ts.begin(); t != ts.end(); ++t)
    if ((*t)->valid() && !(*t)->ghost()) {
      Segments st = segments(*t);
      ss.insert(st.begin(), st.end());
    }
  Segments sm;
  for (SFMap::const_iterator s = sfmap.begin(); s != sfmap.end(); ++s)
    if (s->first.input && ss.find(s->first) == ss.end())
      sm.push_back(s->first);
  return sm;
}

Segments segments (Tetrahedron *t)
{
  Segments segs;
  segs.push_back(Segment(t->vertex(0), t->vertex(2)));
  segs.push_back(Segment(t->vertex(2), t->vertex(3)));
  segs.push_back(Segment(t->vertex(3), t->vertex(0)));
  segs.push_back(Segment(t->vertex(2), t->vertex(1)));
  segs.push_back(Segment(t->vertex(1), t->vertex(3)));
  segs.push_back(Segment(t->vertex(3), t->vertex(2)));
  segs.push_back(Segment(t->vertex(1), t->vertex(0)));
  segs.push_back(Segment(t->vertex(0), t->vertex(3)));
  segs.push_back(Segment(t->vertex(3), t->vertex(1)));
  segs.push_back(Segment(t->vertex(0), t->vertex(1)));
  segs.push_back(Segment(t->vertex(1), t->vertex(2)));
  segs.push_back(Segment(t->vertex(2), t->vertex(0)));
  return segs;
}

void splitSegment (Tetrahedrons &ts, FTMap &ftmap, PTMap &ptmap, SFMap &sfmap,
		   Points &steiner, Octree<Point *> *octree, Segments &sm)
{
  Segment s = sm.back(), s1, s2;
  sm.pop_back();
  if (sfmap.find(s) == sfmap.end())
    return;
  Point *p = splitSegment(steiner, octree, s, s1, s2);
  octree->insert(p);
  Tetrahedron *tp = findTetrahedron(p, ptmap.find(s.a)->second, ftmap);
  unsigned int st = ts.size();
  Tetrahedrons cav = cavity(p, tp, ftmap);
  fillCavity(ts, ftmap, p, cav);
  update(ts, st, ptmap);
  splitFacet(sfmap, s, p);
  splitFacet(sfmap, Segment(s.b, s.a), p);
  SegmentSet sc;
  sc.insert(s1);
  sc.insert(s2);
  missingSegments(ts, st, sfmap, sc, cav, sm);
}

Point * splitSegment (Points &steiner, Octree<Point *> *octree, const Segment &s,
		      Segment &s1, Segment &s2)
{
  Point *p = splitPoint(octree, s);
  steiner.push_back(p);
  s1 = Segment(s.a, p, s.input, s.aflag, false);
  s2 = Segment(p, s.b, s.input, false, s.bflag);
  return p;
}

Point * splitPoint (Octree<Point *> *octree, const Segment &s)
{
  if (!s.aflag && !s.bflag) // eps: added this case
    return new CentroidPoint(s.a, s.b);
  Point *p = referencePoint(octree, s);
  double ab = pointDistance(s.a, s.b);
  if (s.aflag == s.bflag) {
    double ap = pointDistance(s.a, p);
    if (ap < 0.5*ab)
      return new AffinePoint(s.a, s.b, ap/ab);
    double pb = pointDistance(p, s.b);
    if (pb < 0.5*ab)
      return new AffinePoint(s.b, s.a, pb/ab);
    return new CentroidPoint(s.a, s.b);
  }
  Point *a = s.aflag ? s.a : s.b, *b = s.aflag ? s.b : s.a;
  double ap = pointDistance(a, p);
  Point *v = new AffinePoint(a, b, ap/ab);
  if (closerPair(v, p, v, b))
    return v;
  double pv = pointDistance(p, v), av = pointDistance(a, v),
    r = pv < 0.5*av ? av - pv : 0.5*av;
  delete v;
  return new AffinePoint(a, b, r/ab);
}

Point * referencePoint (Octree<Point *> *octree, const Segment &s)
{
  double bb[6];
  encroachedBBox(s, bb);
  PPoints pts;
  octree->find(bb, pts);
  Point *p = 0;
  PTR<Point> o = new CentroidPoint(s.a, s.b);
  for (unsigned int i = 0u; i < pts.size(); ++i) {
    Point *pi = pts[i];
    if (pi != s.a && pi != s.b && EncroachedSegment(s.a, s.b, pi) > -1 &&
	(!p || closerPair(pi, o, p, o)))
	p = pi;
  }
  return p;
}

void encroachedBBox (const Segment &s, double *bb)
{
  PV3<double> a = s.a->getApproxMid(), b = s.b->getApproxMid(),
    m = 0.5*(a + b);
  double r = 0.5*pointDistance(s.a, s.b);
  bb[0] = m.x - r;
  bb[1] = m.x + r;
  bb[2] = m.y - r;
  bb[3] = m.y + r;
  bb[4] = m.z - r;
  bb[5] = m.z + r;
}

double pointDistance (Point *a, Point *b)
{
  return Distance(a, b).getApproxMid();
}

void splitFacet (SFMap &sfmap, const Segment &s, Point *p)
{
  Facet f = sfmap.find(s)->second;
  Segments cb = cavityBoundary(sfmap, p, f);
  fillCavity(sfmap, s, p, cb, f.id);
}

Segments cavityBoundary (SFMap &sfmap, Point *p, const Facet &f)
{
  Segments cb;
  Facets st;
  st.push_back(f);
  while (!st.empty()) {
    Facet f = st.back();
    st.pop_back();
    for (unsigned int i = 0u; i < 3u; ++i) {
      SFMap::iterator si = sfmap.find(Segment(f.vertex(i), f.vertex((i+1)%3)));
      Segment s = si->first;
      sfmap.erase(si);
      if (s.input)
	cb.push_back(s);
      else {
	SFMap::iterator gi = sfmap.find(Segment(s.b, s.a));
	if (gi != sfmap.end()) {
	  Facet g = gi->second;
	  if (inDiametricBall(p, g.vertex(0), g.vertex(1), g.vertex(2)))
	    st.push_back(g);
	  else
	    cb.push_back(s);
	}
      }
    }
  }
  return cb;
}

void fillCavity (SFMap &sfmap, const Segment &s, Point *p, const Segments &cb,
		 unsigned int id)
{
  for (Segments::const_iterator c = cb.begin(); c != cb.end(); ++c)
    if (!(*c == s)) {
      Facet f(p, c->a, c->b, id);
      Segment s1(p, c->a, s.b == c->a, false, c->aflag),
	s3(c->b, p, s.a == c->b, c->bflag, false);
      sfmap.insert(SFPair(s1, f));
      sfmap.insert(SFPair(*c, f));
      sfmap.insert(SFPair(s3, f));
    }
}

void missingSegments (const Tetrahedrons &ts, unsigned int s, const SFMap &sfmap,
		      SegmentSet &sc, const Tetrahedrons &cav, Segments &sm)
{
  for (Tetrahedrons::const_iterator t = cav.begin(); t != cav.end(); ++t) {
    Segments se = segments(*t);
    for (Segments::iterator s = se.begin(); s != se.end(); ++s) {
      SFMap::const_iterator si = sfmap.find(*s);
      if (si != sfmap.end() && si->first.input)
	sc.insert(si->first);
    }
  }
  for (unsigned int i = s; i < ts.size(); ++i) {
    Segments se = segments(ts[i]);
    for (Segments::iterator s = se.begin(); s != se.end(); ++s)
      sc.erase(*s);
  }
  sm.insert(sm.end(), sc.begin(), sc.end());
}

void recoverFacets (Tetrahedrons &ts, FTMap &ftmap, PTMap &ptmap,
		    const SFMap &sfmap)
{
  Facets fm = missingFacets(ftmap, sfmap);
  while (!fm.empty())
    recoverFacets(ts, ftmap, ptmap, sfmap, fm);
}

Facets missingFacets (const FTMap &ftmap, const SFMap &sfmap)
{
  Facets fm;
  for (SFMap::const_iterator i = sfmap.begin(); i != sfmap.end(); ++i)
    if (i->first.a == i->second.a & i->first.b == i->second.b &&
	!getTetrahedron(i->second, ftmap))
      fm.push_back(i->second);
  return fm;
}

void recoverFacets (Tetrahedrons &ts, FTMap &ftmap, PTMap &ptmap,
		    const SFMap &sfmap, Facets &fm)
{
  Facets fp = facetPolygon(ftmap, sfmap, fm), fs;
  if (fp.empty())
    return;
  if (!cavityFacets(ftmap, ptmap, fp, fs)) {
    cerr << "cavityFacets false in recoverFacets" << endl;
    exit(0);
  }
  flipInsertPolygon(ts, ftmap, ptmap, fp, fs);
}

Facets facetPolygon (const FTMap &ftmap, const SFMap &sfmap, Facets &fm)
{
  Facets st, fs;
  FacetSet done;
  st.push_back(fm.back());
  fm.pop_back();
  while (!st.empty()) {
    Facet f = st.back();
    st.pop_back();
    if (inputFacet(sfmap, f, f) && done.insert(f).second &&
	!getTetrahedron(f, ftmap)) {
      fs.push_back(f);
      for (unsigned int i = 0u; i < 3u; ++i) {
	Point *a = f.vertex(i), *b = f.vertex((i+1)%3);
	Segment s;
	if (!inputSegment(sfmap, a, b, s)) {
	  Facet g = sfmap.find(Segment(b, a))->second;
	  if (done.find(g) == done.end())
	    st.push_back(g);
	}
      }
    }
  }
  return fs;
}

bool inputSegment (const SFMap &sfmap, Point *a, Point *b, Segment &s)
{
  SFMap::const_iterator si = sfmap.find(Segment(a, b));
  if (si == sfmap.end())
    return false;
  s = si->first;
  return s.input;
}

bool inputFacet (const SFMap &sfmap, const Facet &f, Facet &g)
{
  SFMap::const_iterator i = sfmap.find(Segment(f.a, f.b));
  if (i == sfmap.end() || i->second.oppositePoint(f.a, f.b) != f.c)
    return false;
  g = i->second;
  return true;
}

bool cavityFacets (const FTMap &ftmap, const PTMap &ptmap, const Facets &fp,
		   Facets &fs)
{
  Tetrahedrons st, cav;
  st.push_back(containingTetrahedron(ftmap, ptmap, fp[0]));
  while (!st.empty()) {
    Tetrahedron *t = st.back();
    st.pop_back();
    if (t->valid()) {
      t->setValid(false);
      cav.push_back(t);
      for (unsigned int i = 0u; i < 4u; ++i) {
	Facet f = t->facet(i);
	Tetrahedron *u = getOppositeTetrahedron(f, ftmap);
	if (u->valid() && !u->ghost()) {
	  bool flag = false;
	  for (Facets::const_iterator g = fp.begin(); !flag && g != fp.end(); ++g) {
	    int s = intersects(f, *g);
	    if (s == 0) {
	      setValid(cav, true);
	      return false;
	    }
	    flag = s == 1;
	  }
	  if (flag) {
	    st.push_back(u);
	    fs.push_back(f);
	  }
	}
      }
    }
  }
  setValid(cav, true);
  return true;
}

Tetrahedron * containingTetrahedron (const FTMap &ftmap, const PTMap &ptmap,
				     const Facet &f)
{
  Tetrahedron *t = ptmap.find(f.a)->second;
  PTR<Point> p = new CentroidPoint(f.a, f.b, f.c),
    o = new CentroidPoint(t->vertex(0), t->vertex(1), t->vertex(2), t->vertex(3));
  while (!contains(t, p)) {
    Facet f = t->exitFacet(o, p);
    t = getOppositeTetrahedron(f, ftmap);
  }
  return t;
}

bool contains (Tetrahedron *t, Point *p)
{
  for (unsigned int i = 0u; i < 4u; ++i) {
    Facet f = t->facet(i);
    PTR<TrianglePlane> pl = new TrianglePlane(f.a, f.b, f.c);
    if (Side(pl, p) == 1)
      return false;
  }
  return true;
}

int intersects (const Facet &f, const Facet &g)
{
  PTR<TrianglePlane> p = new TrianglePlane(f.a, f.b, f.c),
    q = new TrianglePlane(g.a, g.b, g.c);
  return intersects(p, q);
}

int intersects (TrianglePlane *p, TrianglePlane *q)
{
  int s = intersectsAux(p, q);
  return s == -1 ? intersectsAux(q, p) : s;
}

int intersectsAux (TrianglePlane *p, TrianglePlane *q)
{
  Point *a = q->getA(), *b = q->getB(), *c = q->getC();
  int sa = Side(p, a), sb = Side(p, b);
  if (sa*sb == -1 && intersectTriangleLine(p, a, b) == 1)
    return 1;
  int sc = Side(p, c);
  if (sb*sc == -1 && intersectTriangleLine(p, b, c) == 1 ||
      sc*sa == -1 && intersectTriangleLine(p, c, a) == 1)
    return 1;
 return sa == 0 && sb == 0 && sc == 0 && overlaps(p, q) ? 0 : -1;
}

bool overlaps (TrianglePlane *p, TrianglePlane *q)
{
  Point *v[] = {q->getA(), q->getB(), q->getC()};
  for (unsigned int i = 0u; i < 3u; ++i)
    if (contains(p, v[i]) == 1 ||
	intersectsEdges(p, v[i], v[(i+1u)%3u]))
      return true;
  return false;
}

int contains (TrianglePlane *p, Point *a)
{
  if (a == p->getA() || a == p->getB() || a == p->getC())
    return -2;
  int pc = p->getPC(),
    lt1 = LeftTurn(p->getA(), p->getB(), a, pc),
    lt2 = LeftTurn(p->getB(), p->getC(), a, pc),
    lt3 = LeftTurn(p->getC(), p->getA(), a, pc);
  return min(lt1, min(lt2, lt3));
}

bool intersectsEdges (TrianglePlane *p, Point *t, Point *h)
{
  int pc = p->getPC();
  return intersectsEE(p->getA(), p->getB(), t, h, pc) ||
    intersectsEE(p->getB(), p->getC(), t, h, pc) ||
    intersectsEE(p->getC(), p->getA(), t, h, pc);
}

bool intersectsEE (Point *a, Point *b, Point *c, Point *d, int pc)
{
  return LeftTurn(a, b, c, pc)*LeftTurn(a, b, d, pc) == -1 &&
    LeftTurn(c, d, a, pc)*LeftTurn(c, d, b, pc) == -1;
}

void flipInsertPolygon (Tetrahedrons &ts, FTMap &ftmap, PTMap &ptmap,
			const Facets &fp, const Facets &fs)
{
  PTR<TrianglePlane> p = new TrianglePlane(fp[0].a, fp[0].b, fp[0].c);
  FlipQ q = flipInit(ftmap, fs, p);
  while (true) {
    FlipElts fes = nextFlips(q);
    if (fes.empty())
      break;
    bool flag = false;
    int st = ts.size();
    Tetrahedrons cav;
    for (FlipElts::iterator f = fes.begin(); !flag && f != fes.end(); ++f)
      flag = flip(ts, ftmap, ptmap, *f, cav);
    if (!(flag && fes.size() <= 3)) {
      cerr << "bad flip in flipInsertPolyhedron" << endl;
      exit(0);
    }
    certify(ftmap, cav, p, fes[0].s, q);
  }
}

FlipQ flipInit (const FTMap &ftmap, const Facets &fs, TrianglePlane *p)
{
  FlipQ q;
  PTR<Object<Scalar>> ts = new Object<Scalar>(Scalar<double>(0.0), false);
  for (Facets::const_iterator f = fs.begin(); f != fs.end(); ++f)
    certify(ftmap, *f, p, ts, q);
  return q;
}

void certify (const FTMap &ftmap, const Facet &f, TrianglePlane *p,
	      Object<Scalar> *ts, FlipQ &q)
{
  Tetrahedron *t = getTetrahedron(f, ftmap), *u = getOppositeTetrahedron(f, ftmap);
  PPoints pts;
  for (unsigned int i = 0u; i < 3u; ++i)
    pts.push_back(f.vertex(i));
  pts.push_back(t->oppositeVertex(f));
  pts.push_back(u->oppositeVertex(Facet(f.c, f.b, f.a)));
  FlipLifters fl;
  for (PPoints::iterator v = pts.begin(); v != pts.end(); ++v)
    if (Side(p, *v) == 1)
      fl.push_back(new FlipLifter(*v, p->getA(), p->getB(), p->getC()));
    else
      fl.push_back(new FlipLifter(*v));
  PTR<FlipTime> s = new FlipTime(fl[0], fl[1], fl[2], fl[3], fl[4]);
  if (ScalarOrder(ts, s) > -1)
    q.insert(FlipElt(f, t, u, s));
}

FlipElts nextFlips (FlipQ &q)
{
  FlipElts fes;
  while (!q.empty()) {
    FlipElt x = q.min();
    q.removeMin();
    if (x.valid()) {
      fes.push_back(x);
      while (!q.empty()) {
	FlipElt y = q.min();
	if (x == y || ScalarOrder(x.s, q.min().s) == 0) {
	  q.removeMin();
	  if (y.valid() && find(fes.begin(), fes.end(), y) == fes.end())
	    fes.push_back(y);
	}
	else
	  break;
      }
      break;
    }
  }
  return fes;
}

bool flip (Tetrahedrons &ts, FTMap &ftmap, PTMap &ptmap, const FlipElt &fe,
	   Tetrahedrons &cav)
{
  Point *x = fe.t->oppositeVertex(fe.f), *y = fe.u->oppositeVertex(fe.f);
  Tetrahedrons td, tc;
  td.push_back(fe.t);
  td.push_back(fe.u);
  for (unsigned int i = 0u; i < 3u; ++i) {
    Point *a = fe.f.vertex(i), *b = fe.f.vertex((i+1)%3);
    if (Orientation(a, b, x, y) == 1) {
      td.push_back(getTetrahedron(Facet(a, b, x), ftmap));
      if (td.back()->oppositeVertex(Facet(a, b, x)) != y) {
	for (Tetrahedrons::iterator i = tc.begin(); i != tc.end(); ++i)
	  delete *i;
      	return false;
      }
    }
    else
      tc.push_back(new Tetrahedron(a, b, y, x));
  }
  setValid(td, false);
  removeTetrahedrons(td, ftmap);
  for (Tetrahedrons::iterator i = tc.begin(); i != tc.end(); ++i)
    addTetrahedron(ts, ftmap, *i);
  update(ts, ts.size() - tc.size(), ptmap);
  cav = tc;
  return true;
}

void certify (const FTMap &ftmap, Tetrahedrons &cav, TrianglePlane *p,
	      Object<Scalar> *ts, FlipQ &q)
{
  setValid(cav, false);
  Facets fs = ::cavityBoundary(cav, ftmap);
  setValid(cav, true);
  for (Facets::iterator f = fs.begin(); f != fs.end(); ++f) {
    int sa = Side(p, f->a), sb = Side(p, f->b), sc = Side(p, f->c);
    if (sa*sb == -1 || sa*sc == -1 || sb*sc == -1)
      certify(ftmap, *f, p, ts, q);
  }
}

Tetrahedrons innerTetrahedrons (const Tetrahedrons &ts, const FTMap &ftmap,
				const FacetSet &fs)
{
  Tetrahedrons st, res;
  for (FacetSet::const_iterator f = fs.begin(); f != fs.end(); ++f)
    st.push_back(getTetrahedron(*f, ftmap));
  while (!st.empty()) {
    Tetrahedron *t = st.back();
    st.pop_back();
    if (t->valid()) {
      t->setValid(false);
      res.push_back(t);
      for (unsigned int i = 0u; i < 4u; ++i) {
	Facet f = t->facet(i);
	if (fs.find(f) == fs.end()) {
	  Tetrahedron *u = getOppositeTetrahedron(f, ftmap);
	  if (!u->ghost())
	    st.push_back(u);
	}
      }
    }
  }
  for (Tetrahedrons::iterator t = res.begin(); t != res.end(); ++t)
    (*t)->setValid(true);
  return res;
}

Tetrahedron * getTetrahedron (const Facet &f, const FTMap &ftmap)
{
  FTMap::const_iterator i = ftmap.find(f);
  return i == ftmap.end() ? 0 : i->second;
}

void setValid (const Tetrahedrons &ts, bool flag)
{
  for (Tetrahedrons::const_iterator t = ts.begin(); t != ts.end(); ++t)
    (*t)->setValid(flag);
}

// debug

Polyhedron * polyhedron (const SFMap &sfmap)
{
  FacetSet fs;
  for (SFMap::const_iterator iter = sfmap.begin(); iter != sfmap.end(); ++iter)
    fs.insert(iter->second);
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (FacetSet::iterator f = fs.begin(); f != fs.end(); ++f)
    a->addTriangle(f->vertex(0), f->vertex(1), f->vertex(2), pvmap);
  return a;
}

Polyhedron * polyhedron (const SFMap &sfmap, unsigned int id)
{
  FacetSet fs;
  for (SFMap::const_iterator iter = sfmap.begin(); iter != sfmap.end(); ++iter)
    if (iter->second.id == id)
      fs.insert(iter->second);
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (FacetSet::iterator f = fs.begin(); f != fs.end(); ++f)
    a->addTriangle(f->vertex(0), f->vertex(1), f->vertex(2), pvmap);
  return a;
}
