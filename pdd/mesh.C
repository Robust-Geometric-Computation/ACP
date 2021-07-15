#include "mesh.h"

Polyhedron * mesh (Polyhedron *a, double r, double cs, double rm)
{
  Mesh m(a, r, cs, rm);
  Tetrahedrons ts = m.mesh();
  return polyhedron(ts);
}

Tetrahedrons Mesh::mesh ()
{
  double t = getTime();
  init();
  a->computeWindingNumbers();
  refine();
  Tetrahedrons res = innerTetrahedrons();
  t = getTime() - t;
  unsigned int nf = 0, ns = 0u;
  for (Tetrahedrons::iterator t = res.begin(); t != res.end(); ++t)
    if (fenced.find(*t) != fenced.end())
      ++nf;
    else if (sliver(*t))
      ++ns;
  cerr << "tetrahedrons " << res.size() << " fenced " << nf << " sliver " << ns
       << " cpu time " << t << endl;
  acp::primitiveReport();
  return res;
}

void Mesh::init ()
{
  initSFMap();
  Points pts;
  for (Vertices::iterator v = a->vertices.begin(); v != a->vertices.end(); ++v)
    pts.push_back((*v)->getP());
  poctree = Octree<PTR<Point>>::octree(pts, a->bbox, 0.0, 1, 100);
  initFoctree();
  initEntwined();
  initPRMap();
  PPoints ppts;
  for (Points::iterator p = pts.begin(); p != pts.end(); ++p)
    ppts.push_back(*p);
  ts = delaunay(ppts, a->bbox, ftmap);
  updatePTMap(0u);
  SFPairs enc = encroachedSegments();
  while (!enc.empty())
    splitSegment(enc);
  initSoctree();
}

void Mesh::initSFMap ()
{
  unsigned int id = 1u;
  set<Segment> ss;
  for (Edges::const_iterator e = a->edges.begin(); e != a->edges.end(); ++e) {
    Point *t = (*e)->getT()->getP(), *h = (*e)->getH()->getP();
    ss.insert(Segment(t, h, id));
    ss.insert(Segment(h, t, id));
    ++id;
  }
  Facets fs;
  PPointSet acute;
  for (Faces::const_iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    Points pts = (*f)->boundaryPoints();
    fs.push_back(Facet(pts[0], pts[1], pts[2], id));
    ++id;
    for (unsigned int i = 0u; i < 3u; ++i) {
      Point *a = pts[i], *b = pts[(i+1u)%3u], *c = pts[(i+2u)%3u];
      if (Acute(a, b, c) == 1) {
	acute.insert(b);
	ent.insert(IIPair(ss.find(Segment(a, b))->id,
			  ss.find(Segment(b, c))->id));
      }
    }
  }
  for (Facets::iterator f = fs.begin(); f != fs.end(); ++f)
    for (unsigned int i = 0u; i < 3u; ++i) {
      Segment s = *ss.find(Segment(f->vertex(i), f->vertex((i+1u)%3u)));
      s.aflag = acute.find(s.a) != acute.end();
      s.bflag = acute.find(s.b) != acute.end();
      sfmap.insert(SFPair(s, *f));
    }
}

unsigned int Mesh::segmentId (Point *a, Point *b) const
{
  Segment s;
  return inputSegment(a, b, s) ? s.id : 0u;
}

bool Mesh::inputSegment (Point *a, Point *b, Segment &s) const
{
  SFMap::const_iterator si = sfmap.find(Segment(a, b));
  if (si == sfmap.end())
    return false;
  s = si->first;
  return s.id != 0u;
}

void Mesh::initFoctree ()
{
  FacetPlanes trs;
  for (SFMap::iterator si = sfmap.begin(); si != sfmap.end(); ++si) {
    Facet &f = si->second;
    if (si->first.a == f.a && si->first.b == f.b)
      trs.push_back(new FacetPlane(f.a, f.b, f.c, f.id));
  }
  foctree = Octree<PTR<FacetPlane>>::octree(trs, a->bbox, 0.0, 1, 100);
}

void Mesh::initEntwined ()
{
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if (AcuteEdge(*e) == 1) {
      Points b1 = (*e)->getHEdge(0)->pointLoop(),
	b2 = (*e)->getHEdge(1)->pointLoop();
      ent.insert(IIPair(facetId(b1[0], b1[1], b1[2]),
			facetId(b2[0], b2[1], b2[2])));
    }
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    HEdges ed = (*f)->boundaryHEdges();
    for (unsigned int i = 0u; i < 3u; ++i) {
      HEdge *e = ed[i]->ccw()->getNext()->getNext()->ccw();
      while (e->getF() != *f) {
	if (entwinedFE(*f, e))
	  ent.insert(IIPair(segmentId(e->tail()->getP(), e->head()->getP()),
			    facetId(ed[0]->tail()->getP(), ed[1]->tail()->getP(),
				    ed[2]->tail()->getP())));
	e = e->getNext()->getNext()->ccw();
      }
    }
  }
}

unsigned int Mesh::facetId (Point *a, Point *b, Point *c) const
{
  Facet f;
  return inputFacet(Facet(a, b, c), f) ? f.id : 0u;
}

bool Mesh::inputFacet (const Facet &f, Facet &g) const
{
  SFMap::const_iterator i = sfmap.find(Segment(f.a, f.b));
  if (i == sfmap.end() || i->second.oppositePoint(f.a, f.b) != f.c)
    return false;
  g = i->second;
  return true;
}

bool Mesh::entwinedFE (Face *f, HEdge *e) const
{
  Point *t = e->tail()->getP();
  TrianglePlane *pl = f->getP();
  Points pts = f->boundaryPoints();
  int pc = pl->getPC();
  PTR<Point> h = new ProjectionPoint(pl, e->head()->getP());
  for (unsigned int i = 0u; i < 3u; ++i)
    if (pts[i] == t) {
      Point *p = pts[(i+1u)%3u], *q = pts[(i+2u)%3u];
      return LeftTurn(h, t, p, pc) == 1 && LeftTurn(q, t, h, pc) == 1;
    }
  return false;
}

void Mesh::initPRMap ()
{
  for (Vertices::iterator v = a->vertices.begin(); v != a->vertices.end(); ++v) {
    Point *p = (*v)->getP();
    HEdges ed = (*v)->outgoingHEdges();
    double d = 1e100;
    for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
      d = min(d, distance(p, (*e)->head()->getP()));
    PPoints pv = visibleVertices(p, d);
    for (PPoints::iterator q = pv.begin(); q != pv.end(); ++q)
      d = min(d, distance(p, *q));
    prmap.insert(PRPair(p, d));
  }
}

PPoints Mesh::visibleVertices (Point *p, double dmax) const
{
  double bb[6];
  p->getBBox(bb);
  for (unsigned int i = 0u; i < 3u; ++i) {
    bb[2*i] -= dmax;
    bb[2*i+1] += dmax;
  }
  Points cand;
  poctree->find(bb, cand);
  PPoints pv;
  for (Points::iterator q = cand.begin(); q != cand.end(); ++q)
    if (*q != p && distance(*q, p) < dmax && visible(p, *q))
      pv.push_back(*q);
  return pv;
}

bool Mesh::visible (Point *p, Point *q) const
{
  double bb[6], bbq[6];
  p->getBBox(bb);
  q->getBBox(bbq);
  mergeBBox(bbq, bb);
  FacetPlanes cand;
  foctree->find(bb, cand);
  for (FacetPlanes::iterator f = cand.begin(); f != cand.end(); ++f)
    if (intersects(*f, p, q) == 1)
      return false;
  return true;
}

void Mesh::updatePRMap (AffinePoint *p, bool cd, Point *par, bool entw)
{
  bool midpoint = p->s == 0.5;
  CircleCenter *parc = dynamic_cast<CircleCenter *>(par);
  if (parc)
    entw = entwined(p->getId(), parc->getId());
  double d = 1e100;
  if (par && midpoint && !entw)
    d = distance(p, p->a);
  else if (par && !midpoint)
    d = min(distance(p, p->a), distance(p, p->b));
  else {
    PPoints pts = cd ? neighborsCD(p) : neighbors(p);
    for (PPoints::iterator q = pts.begin(); q != pts.end(); ++q)
      d = min(d, relaxedDistance(p, *q));
    if (entw) {
      double dpar;
      if (parc)
	dpar = calculateRR(parc);
      else {
	SphereCenter *pars = dynamic_cast<SphereCenter *>(par);
	dpar = distance(pars, pars->a);
      }
      d = max(min(distance(p, p->a), distance(p, p->b)),
	      min(dpar*sqrt(0.5), d));
    }
  }
  prmap.insert(PRPair(p, d));
}

PPoints Mesh::neighbors (AffinePoint *p) const
{
  PPoints pts = visibleVertices(p, distance(p, p->a));
  if (find(pts.begin(), pts.end(), p->a) == pts.end())
    pts.push_back(p->a);
  if (find(pts.begin(), pts.end(), p->b) == pts.end())
    pts.push_back(p->b);
  return pts;
}

PPoints Mesh::neighborsCD (Point *p) const
{
  Tetrahedrons tp = vertexTetrahedrons(p);
  setValid(tp, true);
  PPoints pts;
  PPointSet done;
  for (Tetrahedrons::iterator t = tp.begin(); t != tp.end(); ++t)
    for (unsigned int i = 0u; i < 4u; ++i) {
      Point *q = (*t)->vertex(i);
      if (q != p && done.insert(q).second)
	pts.push_back(q);
    }
  return pts;
}

Tetrahedrons Mesh::vertexTetrahedrons (Point *p) const
{
  Tetrahedrons st, tn;
  st.push_back(ptmap.find(p)->second);
  while (!st.empty()) {
    Tetrahedron *t = st.back();
    st.pop_back();
    if (t->valid()) {
      t->setValid(false);
      tn.push_back(t);
      for (unsigned int i = 0u; i < 4u; ++i) {
	Tetrahedron *u = getOppositeTetrahedron(t->facet(i), ftmap);
	if (u->valid() && !u->ghost() && u->isVertex(p))
	  st.push_back(u);
      }
    }
  }
  return tn;
}

double Mesh::relaxedDistance (AffinePoint *p, Point *q) const
{
  double d = distance(p, q);
  bool flag = p->s == 0.5 && (q == p->a || q == p->b);
  if (!flag) {
    AffinePoint *qa = dynamic_cast<AffinePoint *>(q);
    if (qa)
      flag = entwined(p->getId(), qa->getId());
    else {
      CircleCenter *qc = dynamic_cast<CircleCenter *>(q);
      flag = qc && entwined(p->getId(), qc->getId());
    }
  }
  return flag ? max(d, getRR(q)) : d;
}

void Mesh::updatePRMap (CircleCenter *p)
{
  double d = calculateRR(p);
  prmap.insert(PRPair(p, d));
}

double Mesh::calculateRR (CircleCenter *p) const
{
  Facet f(p->getA(), p->getB(), p->getC());
  double d = distance(p, f.a);
  PPoints pn = visibleVertices(p, d);
  for (PPoints::iterator q = pn.begin(); q != pn.end(); ++q)
    if (eligible(*q, f) && inDiametricBall(*q, f.a, f.b, f.c))
      d = min(d, distance(p, *q));
  return d;
}

bool Mesh::eligible (Point *p, const Facet &f) const
{
  PTR<Plane> pl = new TrianglePlane(f.a, f.b, f.c);
  if (Side(pl, p) == 0)
    return false;
  AffinePoint *pa = dynamic_cast<AffinePoint *>(p);
  if (pa)
    return !(entwined(pa->getId(), f.id) &&
	     PointNearPlane(p, pl, sqrt(2.0)*getRR(p)) == 1);
  CircleCenter *pc = dynamic_cast<CircleCenter *>(p);
  return !(pc && entwined(pc->getId(), f.id)
	   && PointNearPlane(p, pl, getRR(p)) == 1);
}

void Mesh::updatePTMap (unsigned int s)
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

SFPairs Mesh::encroachedSegments () const
{
  SegmentSet cand;
  for (SFMap::const_iterator s = sfmap.begin(); s != sfmap.end(); ++s)
    if (s->first.id)
      cand.insert(s->first);
  SFPairs enc;
  encroachedSegments(0u, cand, enc);
  return enc;
}

void Mesh::encroachedSegments (unsigned int st, SegmentSet &cand, SFPairs &enc) const
{
  for (Tetrahedrons::const_iterator t = ts.begin() + st; t != ts.end(); ++t)
    if ((*t)->valid() && !(*t)->ghost())
      for (unsigned int i = 0u; i < 4u; ++i) {
	Facet f = (*t)->facet(i);
	for (unsigned int j = 0u; j < 3u; ++j) {
	  Segment s;
	  if (inputSegment(f.vertex(j), f.vertex((j+1)%3), s)) {
	    cand.erase(s);
	    if (EncroachedSegment(s.a, s.b, f.vertex((j+2)%3)) == 1)
	      enc.push_back(SFPair(s, f));
	  }
	}
      }
  for (SegmentSet::iterator s = cand.begin(); s != cand.end(); ++s)
    enc.push_back(SFPair(*s, Facet()));
}

void Mesh::splitSegment (SFPairs &enc)
{
  Segment s = enc.back().first;
  Facet f = enc.back().second;
  enc.pop_back();
  SFMap::const_iterator si = sfmap.find(s);
  if (si == sfmap.end() ||
      f.a && (f.id ? !(f == si->second) : !getTetrahedron(f, ftmap)))
    return;
  Segment sr, s1, s2;
  inputSegment(s.b, s.a, sr);
  AffinePoint *p = splitSegment(s, s1, s2);
  Tetrahedron *tp = findTetrahedron(p, ptmap.find(s.a)->second, ftmap);
  unsigned int st = ts.size();
  Tetrahedrons cav = cavity(p, tp, ftmap);
  ::fillCavity(ts, ftmap, p, cav);
  updatePTMap(st);
  updatePRMap(p);
  splitSegmentFacet(s, p);
  splitSegmentFacet(sr, p);
  SegmentSet cand = inputSegments(cav);
  cand.insert(s1);
  cand.insert(s2);
  encroachedSegments(st, cand, enc);
}

AffinePoint * Mesh::splitSegment (const Segment &s, Segment &s1, Segment &s2)
{
  AffinePoint *p = splitPoint(s);
  poctree->insert(PTR<Point>(p));
  s1 = Segment(s.a, p, s.id, false, false);
  s2 = Segment(p, s.b, s.id, false, s.aflag && s.bflag);
  return p;
}

AffinePoint * Mesh::splitPoint (const Segment &s) const
{
  double ab = distance(s.a, s.b);
  if (s.aflag && s.bflag) {
    double x = pow(2.0, floor(log2(ab)) - 1)/ab;
    return new AffinePoint(s.a, s.b, x, s.id);
  }
  if (s.aflag) {
    double x = pow(2.0, ceil(log2(ab/3.0)))/ab;
    return new AffinePoint(s.a, s.b, x, s.id);
  }
  if (s.bflag) {
    double x = pow(2.0, ceil(log2(ab/3.0)))/ab;
    return new AffinePoint(s.b, s.a, x, s.id);
  }
  return new AffinePoint(s.a, s.b, 0.5, s.id);
}

void Mesh::splitSegmentFacet (const Segment &s, Point *p, Facets *fm)
{
  Facet f = sfmap.find(s)->second;
  Facets fs;
  Segments cb = cavityBoundary(p, f, fs);
  eraseFacets(fs);
  fillCavity(s, p, cb, f.id);
  if (fm)
    boundaryFacets(cb, *fm);
}

Segments Mesh::cavityBoundary (Point *p, const Facet &f, Facets &fs)
{
  Segments cb;
  Facets st;
  FacetSet done;
  st.push_back(f);
  while (!st.empty()) {
    Facet f = st.back();
    st.pop_back();
    if (done.insert(f).second) {
      fs.push_back(f);
      for (unsigned int i = 0u; i < 3u; ++i) {
	Segment s = sfmap.find(Segment(f.vertex(i), f.vertex((i+1)%3)))->first;
	if (s.id)
	  cb.push_back(s);
	else {
	  Facet g = sfmap.find(Segment(s.b, s.a))->second;
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

void Mesh::eraseFacets (const Facets &fs)
{
  for (Facets::const_iterator f = fs.begin(); f != fs.end(); ++f)
    for (unsigned int i = 0u; i < 3u; ++i)
      sfmap.erase(Segment(f->vertex(i), f->vertex((i+1)%3)));
}

void Mesh::fillCavity (const Segment &s, Point *p, const Segments &cb,
		       unsigned int id)
{
  for (Segments::const_iterator c = cb.begin(); c != cb.end(); ++c)
    if (!(*c == s)) {
      Facet f(p, c->a, c->b, id);
      Segment s1(p, c->a, s.b == c->a ? s.id : 0u, false, c->aflag),
	s3(c->b, p, s.a == c->b ? s.id : 0u, c->bflag, false);
      sfmap.insert(SFPair(s1, f));
      sfmap.insert(SFPair(*c, f));
      sfmap.insert(SFPair(s3, f));
    }
}

void Mesh::boundaryFacets (const Segments &cb, Facets &fm) const
{
  for (Segments::const_iterator s = cb.begin(); s != cb.end(); ++s) {
    SFMap::const_iterator si = sfmap.find(*s);
    if (si != sfmap.end() && !getTetrahedron(si->second, ftmap))
      fm.push_back(si->second);
  }
}

SegmentSet Mesh::inputSegments (const Tetrahedrons &cav) const
{
  SegmentSet ss;
  for (Tetrahedrons::const_iterator t = cav.begin(); t != cav.end(); ++t) {
    Segments se = segments(*t);
    for (Segments::iterator s = se.begin(); s != se.end(); ++s) {
      Segment si;
      if (inputSegment(s->a, s->b, si))
	ss.insert(si);
    }
  }
  return ss;
}

void Mesh::initSoctree ()
{
  set<Segment> ss;
  for (SFMap::iterator i = sfmap.begin(); i != sfmap.end(); ++i) {
    const Segment &s = i->first;
    if (s.id != 0u && s.a < s.b)
      ss.insert(s);
  }
  SegmentCCs sv;
  for (set<Segment>::iterator s = ss.begin(); s != ss.end(); ++s)
    sv.push_back(new SegmentCC(*s));
  soctree = Octree<PTR<SegmentCC>>::octree(sv, a->bbox, 0.0, 1, 100);
}

void Mesh::refine ()
{
  Splits sp;
  Facets fm = missingFacets();
  FTPairs fenc = encroachedFacets();
  SkinnyQ sk;
  bool skflag = false;
  while (true) {
    bool check = true;
    int nt = ts.size();
    if (!sp.empty())
      splitSegment(sp, fm, fenc);
    else if (!fm.empty())
      recoverFacet(fm, fenc);
    else if (!fenc.empty())
      splitFacet(sp, fm, fenc);
    else if (!skflag) {
      skflag = true;
      sk = skinnyTetrahedrons();
    }
    else if (!sk.empty()) {
      check = false;
      splitTetrahedron(sp, fm, fenc, sk);
    }
    else
      break;
    if (skflag && check)
      skinnyTetrahedrons(nt, sk);
  }
}

Facets Mesh::missingFacets () const
{
  Facets fm;
  for (SFMap::const_iterator i = sfmap.begin(); i != sfmap.end(); ++i)
    if (i->first.a == i->second.a & i->first.b == i->second.b &&
	!getTetrahedron(i->second, ftmap))
      fm.push_back(i->second);
  return fm;
}

FTPairs Mesh::encroachedFacets () const
{
  FTPairs fenc;
  encroachedFacets(0u, fenc);
  return fenc;
}

void Mesh::encroachedFacets (unsigned int s, FTPairs &enc) const
{
  for (Tetrahedrons::const_iterator t = ts.begin() + s; t != ts.end(); ++t)
    if ((*t)->valid() && !(*t)->ghost())
      for (unsigned int i = 0u; i < 4u; ++i) {
	Facet f = (*t)->facet(i);
	encroachedFacet((*t)->oppositeVertex(f), f, enc);
      }
}

void Mesh::encroachedFacet (Point *p, const Facet &f, FTPairs &enc) const
{
  Facet g;
  if (inputFacet(f, g) && eligible(p, g) &&
      InDiametricEllipse(p, g.a, g.b, g.c) == 1) {
    Facet h = orthogonalProjection(p, g);
    Tetrahedron *u = getTetrahedron(h, ftmap);
    enc.push_back(FTPair(h, u));
  }
}

Facet Mesh::orthogonalProjection (Point *p, const Facet &f) const
{
  Facets st;
  FacetSet done;
  st.push_back(f);
  while (!st.empty()) {
    Facet g = st.back();
    st.pop_back();
    if (done.insert(g).second) {
      if (containsProjection(g.a, g.b, g.c, p))
	return g;
      for (unsigned int i = 0u; i < 3u; ++i) {
	Segment s = sfmap.find(Segment(g.vertex(i), g.vertex((i+1)%3)))->first;
	if (s.id == 0u)
	  st.push_back(sfmap.find(Segment(s.b, s.a))->second);
      }
    }
  }
  return Facet();
}

SkinnyQ Mesh::skinnyTetrahedrons () const
{
  Tetrahedrons ti = innerTetrahedrons();
  SkinnyQ sk;
  for (Tetrahedrons::const_iterator t = ti.begin(); t != ti.end(); ++t)
    if ((*t)->valid() && !(*t)->ghost() && fenced.find(*t) == fenced.end())
      checkSkinny(*t, sk);
  return sk;
}

Tetrahedrons Mesh::innerTetrahedrons () const
{
  FacetSet fs;
  for (SFMap::const_iterator i = sfmap.begin(); i != sfmap.end(); ++i)
    fs.insert(i->second);
  return ::innerTetrahedrons(ftmap, fs);
}

bool Mesh::checkSkinny (Tetrahedron *t, SkinnyQ &sk) const
{
  PTR<Point> c = new SphereCenter(t->vertex(0), t->vertex(1),
				  t->vertex(2), t->vertex(3));
  PV3<double> q = c->getApproxMid(1e-8);
  c = new SphereCenter(q.x, q.y, q.z);
  double lt = getRR(t->vertex(0));
  for (unsigned int i = 1u; i < 4u; ++i)
    lt = min(lt, getRR(t->vertex(i)));
  Segments ss = segments(t);
  Segment sm = ss[0];
  for (unsigned int i = 1u; i < 6u; ++i)
    if (closerPair(ss[i].a, ss[i].b, sm.a, sm.b))
      sm = ss[i];
  double dc = distance(c, t->vertex(0)), dsm = distance(sm.a, sm.b),
    tr = dc/max(lt, dsm);
  if (tr > r) {
    sk.insert(Skinny(t, c, c, tr));
    return true;
  }
  if (sliver(t) && dc > rm) {
    sk.insert(Skinny(t, c, c, 0.0));
    return true;
  }
  return false;
}

bool Mesh::sliver (Tetrahedron *t) const
{
  const unsigned int k[] = {0, 2, 3, 1, 2, 1, 3, 0, 1, 0, 3, 2, 0, 1, 2, 3};
  if (cs == 1.0)
    return false;
  for (unsigned int i = 0u; i < 4u; ++i)
    if (SmallAngle(t->vertex(k[4u*i]), t->vertex(k[4u*i+1]),
		   t->vertex(k[4u*i+2]), t->vertex(k[4u*i+3]), cs) == 1)
      return true;
  return false;
}

void Mesh::skinnyTetrahedrons (unsigned int st, SkinnyQ &sk) const
{
  for (Tetrahedrons::const_iterator t = ts.begin() + st; t != ts.end(); ++t)
    if ((*t)->valid() && !(*t)->ghost() && fenced.find(*t) == fenced.end() &&
	a->contains(PTR<Point>(new CentroidPoint((*t)->vertex(0), (*t)->vertex(1),
						 (*t)->vertex(2), (*t)->vertex(3)))))
      checkSkinny(*t, sk);
}

void Mesh::splitSegment (Splits &sp, Facets &fm, FTPairs &fenc)
{
  Segment s = sp.back().s;
  Point *par = sp.back().par;
  bool ent = sp.back().ent;
  sp.pop_back();
  if (sfmap.find(s) == sfmap.end())
    return;
  Tetrahedrons tp = segmentTetrahedrons(s);
  Segment sr, s1, s2;
  inputSegment(s.b, s.a, sr);
  AffinePoint *p = splitSegment(s, s1, s2);
  soctree->insert(PTR<SegmentCC>(new SegmentCC(s1)));
  soctree->insert(PTR<SegmentCC>(new SegmentCC(s2)));
  insertVertex(p, tp, sp, fm);
  splitSegmentFacet(s, p, &fm);
  splitSegmentFacet(sr, p, &fm);
  encroachedSegments(sp, s1, s2, p, par, ent);
  bool cd = sp.empty() && fm.empty();
  updatePRMap(p, cd, par, ent);
}

Tetrahedrons Mesh::segmentTetrahedrons (const Segment &s) const
{
  Tetrahedrons st, cav;
  st.push_back(intersectingTetrahedron(s));
  while (!st.empty()) {
    Tetrahedron *t = st.back();
    st.pop_back();
    if (t->valid()) {
      t->setValid(false);
      cav.push_back(t);
      for (unsigned int i = 0u; i < 4u; ++i) {
	Tetrahedron *u = getOppositeTetrahedron(t->facet(i), ftmap);
	if (u->valid() && !u->ghost() && intersects(u, s))
	  st.push_back(u);
      }
    }
  }
  setValid(cav, true);
  return cav;
}

Tetrahedron * Mesh::intersectingTetrahedron (const Segment &s) const
{
  Tetrahedron *t = getTetrahedron(sfmap.find(s)->second, ftmap);
  if (t)
    return t;
  PTR<Point> p = new CentroidPoint(s.a, s.b);
  t = ptmap.find(s.a)->second;
  if (!t->valid())
    t = ptmap.find(s.b)->second;
  if (!t->valid())
    t = findValid();
  return containingTetrahedron(p, t);
}

Tetrahedron * Mesh::findValid () const
{
  for (Tetrahedrons::const_iterator t = ts.begin(); t != ts.end(); ++t)
    if ((*t)->valid() && !(*t)->ghost())
      return *t;
  return 0;
}

Tetrahedron * Mesh::containingTetrahedron (Point *p, Tetrahedron *t) const
{
  PTR<Point> o = new CentroidPoint(t->vertex(0), t->vertex(1), t->vertex(2),
				   t->vertex(3));
  while (!contains(t, p)) {
    Facet f = t->exitFacet(o, p);
    t = getOppositeTetrahedron(f, ftmap);
  }
  return t;
}

void Mesh::encroachedSegments (Splits &sp, const Segment &s1, const Segment &s2,
			       Point *p, Point *par, bool ent) const
{
  Segments se = encroachedSegments(p);
  for (Segments::iterator s = se.begin(); s != se.end(); ++s)
    sp.push_back(Split(*s, par, ent));
  encroachedSegment(sp, s1);
  encroachedSegment(sp, s2);
}

Segments Mesh::encroachedSegments (Point *p) const
{
  double bb[6];
  p->getBBox(bb);
  SegmentCCs cand;
  soctree->find(bb, cand);
  Segments enc;
  for (SegmentCCs::iterator s = cand.begin(); s != cand.end(); ++s)
    if (sfmap.find(**s) != sfmap.end() &&
	EncroachedSegment((*s)->a, (*s)->b, p) == 1)
      enc.push_back(**s);
  return enc;
}

void Mesh::encroachedSegment (Splits &sp, const Segment &s) const
{
  SegmentCC scc(s);
  Points pts;
  poctree->find(scc.bbox, pts);
  for (Points::const_iterator p = pts.begin(); p != pts.end(); ++p)
    if (EncroachedSegment(s.a, s.b, *p) == 1)
      sp.push_back(Split(s, 0, false));
}

void Mesh::recoverFacet (Facets &fm, FTPairs &fenc)
{
  Facets fp = facetPolygon(fm);
  if (fp.empty())
    return;
  Facets fs;
  if (cavityFacets(fp, fs))
    flipInsertPolygon(fp, fs);
  else {
    Facets nfm;
    for (Facets::iterator f = fm.begin(); f != fm.end(); ++f)
      if (find(fp.begin(), fp.end(), *f) == fp.end())
	nfm.push_back(*f);
    fm = nfm;
    for (Facets::iterator f = fp.begin(); f != fp.end(); ++f)
      fenc.push_back(FTPair(*f, 0));
  }
}

Facets Mesh::facetPolygon (Facets &fm) const
{
  Facets st, fs;
  FacetSet done;
  st.push_back(fm.back());
  fm.pop_back();
  while (!st.empty()) {
    Facet f = st.back();
    st.pop_back();
    if (inputFacet(f, f) && done.insert(f).second && !getTetrahedron(f, ftmap)) {
      fs.push_back(f);
      for (unsigned int i = 0u; i < 3u; ++i) {
	Point *a = f.vertex(i), *b = f.vertex((i+1)%3);
	Segment s;
	if (!inputSegment(a, b, s)) {
	  Facet g = sfmap.find(Segment(b, a))->second;
	  if (done.find(g) == done.end())
	    st.push_back(g);
	}
      }
    }
  }
  return fs;
}

bool Mesh::cavityFacets (const Facets &fp, Facets &fs) const
{
  Tetrahedrons st, cav;
  st.push_back(containingTetrahedron(fp[0]));
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

Tetrahedron * Mesh::containingTetrahedron (const Facet &f) const
{
  Tetrahedron *t = ptmap.find(f.a)->second;
  if (!t->valid())
    t = ptmap.find(f.b)->second;
  if (!t->valid())
    t = ptmap.find(f.c)->second;
  if (!t->valid())
    t = findValid();
  PTR<Point> p = new CentroidPoint(f.a, f.b, f.c);
  return containingTetrahedron(p, t);
}

void Mesh::flipInsertPolygon (const Facets &fp, const Facets &fs)
{
  PTR<TrianglePlane> p = new TrianglePlane(fp[0].a, fp[0].b, fp[0].c);
  FlipQ q = flipInit(fs, p);
  while (true) {
    FlipElts fes = nextFlips(q);
    if (fes.empty())
      break;
    bool flag = false;
    Tetrahedrons cav;
    for (FlipElts::iterator f = fes.begin(); !flag && f != fes.end(); ++f)
      flag = flip(*f, cav);
    if (!(flag && fes.size() <= 3)) {
      cerr << "bad flip in Mesh::flipInsertPolygon" << endl;
      exit(0);
    }
    certify(cav, p, fes[0].s, q);
  }
}

FlipQ Mesh::flipInit (const Facets &fs, TrianglePlane *p) const
{
  FlipQ q;
  PTR<Object<Scalar>> ts = new Object<Scalar>(Scalar<double>(0), false);
  for (Facets::const_iterator f = fs.begin(); f != fs.end(); ++f)
    certify(*f, p, ts, q);
  return q;
}

void Mesh::certify (const Facet &f, TrianglePlane *p, Object<Scalar> *ts, FlipQ &q) const
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

FlipElts Mesh::nextFlips (FlipQ &q) const
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

bool Mesh::flip (const FlipElt &fe, Tetrahedrons &cav)
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
  updatePTMap(ts.size() - tc.size());
  cav = tc;
  return true;
}

void Mesh::certify (Tetrahedrons &cav, TrianglePlane *p, Object<Scalar> *ts, FlipQ &q) const
{
  setValid(cav, false);
  Facets fs = ::cavityBoundary(cav, ftmap);
  setValid(cav, true);
  for (Facets::iterator f = fs.begin(); f != fs.end(); ++f) {
    int sa = Side(p, f->a), sb = Side(p, f->b), sc = Side(p, f->c);
    if (sa*sb == -1 || sa*sc == -1 || sb*sc == -1)
      certify(*f, p, ts, q);
  }
}

void Mesh::splitFacet (Splits &sp, Facets &fm, FTPairs &fenc)
{
  Facet f = fenc.back().first;
  Tetrahedron *t = fenc.back().second;
  fenc.pop_back();
  if (!inputFacet(f, f) || t && !t->valid())
    return;
  CircleCenter *p = new CircleCenter(f.a, f.b, f.c, f.id);
  Segments senc = encroachedSegments(p);
  if (senc.empty()) {
    poctree->insert(PTR<Point>(p));
    splitFacet(sp, fm, fenc, f, p, t);
  }
  else {
    ppars.push_back(p);
    sp.push_back(Split(senc[0], p, false));
    fenc.push_back(FTPair(f, t));
  }
}

void Mesh::splitFacet (Splits &sp, Facets &fm, FTPairs &fenc, const Facet &f,
		       CircleCenter *p, Tetrahedron *t)
{
  Tetrahedron *u = t == 0 ? containingTetrahedron(f) : t;
  Tetrahedrons tp = facetTetrahedrons(f, p, u);
  int nt = ts.size();
  insertVertex(p, tp, sp, fm);
  splitFacet(f, p, fm);
  if (t)
    removeType3Neighbors(p, fm);
  updatePRMap(p);
  encroachedFacets(nt, fenc);
}

Tetrahedrons Mesh::facetTetrahedrons (const Facet &f, Point *p, Tetrahedron *t) const
{
  Tetrahedrons tp;
  tp.push_back(t);
  if (t == getTetrahedron(f, ftmap))
    tp.push_back(getOppositeTetrahedron(f, ftmap));
  else
    for (unsigned int i = 0u; i < 4u; ++i) {
      Facet g = t->facet(i);
      if (Orientation(g.a, g.b, g.c, p) == 0) {
	tp.push_back(getOppositeTetrahedron(g, ftmap));
	break;
      }
    }
  return tp;
}

void Mesh::splitFacet (const Facet &f, Point *p, Facets &fm)
{
  Facets fs;
  Segments cb = cavityBoundary(p, f, fs);
  eraseFacets(fs);
  fillCavity(Segment(), p, cb, f.id);
  boundaryFacets(cb, fm);
}

void Mesh::removeType3Neighbors (Point *p, Facets &fm)
{
  return;
  PPoints pn = neighborsCD(p);
  for (PPoints::iterator q = pn.begin(); q != pn.end(); ++q)
    if (dynamic_cast<SphereCenter *>(*q))
      removeVertex(*q, fm);
}

void Mesh::removeVertex (Point *p, Facets &fm)
{
  Tetrahedrons cav = vertexTetrahedrons(p);
  Facets fs = ::cavityBoundary(cav, ftmap);
  removeTetrahedrons(cav, ftmap);
  PPoints ps = points(fs);
  int st = ts.size();
  fillCavity(fm, ps, fs);
  ptmap.erase(p);
  prmap.erase(p);
  poctree->remove(PTR<Point>(p));
  updatePTMap(st);
}

void Mesh::fillCavity (Facets &fm, PPoints &pts, const Facets &fs)
{
  FacetSet fss(fs.begin(), fs.end());
  while (true) {
    FTMap ftmapc;
    Tetrahedrons cav = delaunay(pts, ftmapc);
    if (checkCavity(fm, pts, fss, ftmapc)) {
      Tetrahedrons ti = ::innerTetrahedrons(ftmapc, fss);
      for (Tetrahedrons::iterator t = ti.begin(); t != ti.end(); ++t)
	addTetrahedron(ts, ftmap, (*t)->vertex(0), (*t)->vertex(1),
		       (*t)->vertex(2), (*t)->vertex(3));
      deleteTetrahedrons(cav);
      return;
    }
   deleteTetrahedrons(cav);
  }
}

bool Mesh::checkCavity (Facets &fm, PPoints &pts, FacetSet &fs,
			const FTMap &ftmapc)
{
  for (FacetSet::iterator f = fs.begin(); f != fs.end(); ++f)
    if (ftmapc.find(*f) == ftmapc.end()) {
      Facet fi;
      if (inputFacet(*f, fi))
	fm.push_back(fi);
      Tetrahedron *t = getOppositeTetrahedron(*f, ftmap);
      Point *p = t->oppositeVertex(*f);
      for (unsigned int i = 0u; i < 4u; ++i) {
	Facet g = t->facet(i);
	if (!fs.erase(Facet(g.c, g.b, g.a)))
	  fs.insert(g);
      }
      t->setValid(false);
      removeTetrahedron(t, ftmap);
      if (find(pts.begin(), pts.end(), p) == pts.end())
	pts.push_back(p);
    return false;
    }
  return true;
}

void Mesh::splitTetrahedron (Splits &sp, Facets &fm, FTPairs &fenc, SkinnyQ &sk)
{
  Tetrahedron *t = sk.min().t;
  if (!(t->valid() && fenced.find(t) == fenced.end())) {
    sk.removeMin();
    return;
  }
  Facet f;
  PTR<Point> cp = occluder(t, sk.min().p, f);
  if (cp)
    splitTetrahedron(f, cp, fenc, sk);
  else
    splitTetrahedron(sp, fenc, sk);
}

Point * Mesh::occluder (Tetrahedron *t, Point *c, Facet &f) const
{
  PTR<Point> q = new CentroidPoint(t->vertex(0), t->vertex(1), t->vertex(2),
				   t->vertex(3));
  while (!contains(t, c)) {
    f = t->exitFacet(q, c);
    if (inputFacet(f, f))
      return new RayFacetPoint(c, q, f);
    t = getOppositeTetrahedron(f, ftmap);
  }
  return 0;
}

void Mesh::splitTetrahedron (Splits &sp, FTPairs &fenc, SkinnyQ &sk)
{
  Tetrahedron *t = sk.min().t;
  PTR<Point> p = sk.min().p;
  Segments senc = encroachedSegments(p);
  for (Segments::iterator s = senc.begin(); s != senc.end(); ++s)
    if (visible(p, *s)) {
      sp.push_back(Split(*s, p, false));
      return;
    }
  Facets fec = encroachedFacets(t, p);
  for (Facets::iterator f = fec.begin(); f != fec.end(); ++f) {
    Face ff(f->a, f->b, f->c);
    PTR<Point> q = new ProjectionPoint(ff.getP(), p);
    if (ff.contains(q, true) && visible(q, p)) {
      fenc.push_back(FTPair(*f, getTetrahedron(*f, ftmap)));
      return;
    }
  }
  for (Facets::iterator f = fec.begin(); f != fec.end(); ++f)
    if (partiallyVisible(p, *f)) {
      Facet g = orthogonalProjection(p, *f);
      fenc.push_back(FTPair(g, getTetrahedron(g, ftmap)));
      return;
    }
  if (!senc.empty()) {
    bool flag = entwinedOccluders(p, senc[0]);
    sp.push_back(Split(senc[0], p, flag));
    return;
  }
  insertVertexT(sk);
}

bool Mesh::visible (Point *p, const Segment &s) const
{
  PTR<FacetPlane> f = new FacetPlane(p, s.a, s.b, 0u);
  FacetPlanes cand;
  foctree->find(f->bbox, cand);
  for (FacetPlanes::iterator g = cand.begin(); g != cand.end(); ++g)
    if (intersects(f, *g) == 1)
      return false;
  return true;
}

Facets Mesh::encroachedFacets (Tetrahedron *t, Point *p)
{
  Facets res;
  Tetrahedrons cav = cavity(p, t, ftmap);
  for (Tetrahedrons::iterator c = cav.begin(); c != cav.end(); ++c)
    if (!(*c)->ghost())
      for (unsigned int i = 0u; i < 4u; ++i) {
	Facet f = (*c)->facet(i), g;
	if (inputFacet(f, g) && InDiametricEllipse(p, f.a, f.b, f.c) == 1)
	  res.push_back(g);
      }
  setValid(cav, true);
  return res;
}

// to do: exact algorithm
bool Mesh::partiallyVisible (Point *p, const Facet &f) const
{
  if (visible(p, f.a) || visible(p, f.b) || visible(p, f.c))
    return true;
  PTR<Point> m = new CentroidPoint(f.a, f.b, f.c);
  return visible(p, m);
}

bool Mesh::entwinedOccluders (Point *c, const Segment &s) const
{
  PTR<FacetPlane> f = new FacetPlane(c, s.a, s.b, 0u);
  FacetPlanes cand;
  foctree->find(f->bbox, cand);
  for (FacetPlanes::iterator g = cand.begin(); g != cand.end(); ++g)
    if (intersects(f, *g) == 1 && !entwined((*g)->id, s.id))
      return false;
  return true;
}

void Mesh::splitTetrahedron (const Facet &f, Point *cp, FTPairs &fenc, SkinnyQ &sk)
{
  Tetrahedron *t = sk.min().t;
  PTR<Point> p = sk.min().p;
  for (unsigned int i = 0u; i < 4u; ++i) {
    Point *w = t->vertex(i);
    if (eligible(w, f) && InDiametricEllipse(w, f.a, f.b, f.c) == 1) {
      Facet g = orthogonalProjection(w, f);
      fenc.push_back(FTPair(g, getTetrahedron(g, ftmap)));
      return;
    }
  }
  Facet g = orthogonalProjection(cp, f);
  if (InDiametricEllipse(p, g.a, g.b, g.c) == 1) {
    fenc.push_back(FTPair(g, getTetrahedron(g, ftmap)));
    return;
  }
  fenced.insert(t);
  sk.removeMin();
}

void Mesh::insertVertex (Point *p, Tetrahedrons &tp, Splits &sp, Facets &fm)
{
  Tetrahedrons cav = cavityCD(p, tp);
  Facets fs = ::cavityBoundary(cav, ftmap);
  FacetSet fss(fs.begin(), fs.end());
  Facet fb;
  while (!starShaped(p, fss, fb)) {
    Tetrahedron *t = getOppositeTetrahedron(fb, ftmap);
    cav.push_back(t);
    t->setValid(false);
    for (unsigned int i = 0u; i < 4u; ++i) {
      Facet g = t->facet(i), h(g.c, g.b, g.a);
      if (!fss.erase(h))
	fss.insert(g);
    }
  }
  removeTetrahedrons(cav, ftmap);
  int st = ts.size();
  for (FacetSet::iterator f = fss.begin(); f != fss.end(); ++f)
    addTetrahedron(ts, ftmap, f->a, f->b, f->c, p);
  updatePTMap(st);
  missingSF(cav, fss, p, sp, fm);
}

Tetrahedrons Mesh::cavityCD (Point *p, const Tetrahedrons &pts) const
{
  Tetrahedrons st = pts, cav;
  while (!st.empty()) {
    Tetrahedron *s = st.back();
    st.pop_back();
    if (s->valid()) {
      s->setValid(false);
      cav.push_back(s);
      for (unsigned int i = 0u; i < 4u; ++i) {
	Facet f = s->facet(i), g(f.c, f.b, f.a), h;
	if (!inputFacet(f, g) && !inputFacet(g, h)) {
	  Tetrahedron *t = getOppositeTetrahedron(f, ftmap);
	  if (t && t->valid() && t->inCircumsphere(p))
	    st.push_back(t);
	}
      }
    }
  }
  return cav;
}

bool Mesh::starShaped (Point *p, const FacetSet &fs, Facet &fb) const
{
  for (FacetSet::const_iterator f = fs.begin(); f != fs.end(); ++f)
    if (!f->ghost() && Orientation(f->a, f->b, f->c, p) != 1) {
      fb = *f;
      return false;
    }
  return true;
}

void Mesh::missingSF (const Tetrahedrons &cav, const FacetSet &fs, Point *par,
		      Splits &sp, Facets &fm) const
{
  SegmentSet ss;
  for (FacetSet::const_iterator f = fs.begin(); f != fs.end(); ++f)
    for (unsigned int i = 0u; i < 3u; ++i)
      ss.insert(Segment(f->vertex(i), f->vertex((i+1u)%3u)));
  for (Tetrahedrons::const_iterator t = cav.begin(); t != cav.end(); ++t)
    for (unsigned int i = 0u; i < 4u; ++i) {
      Facet f = (*t)->facet(i);
      if (fs.find(f) == fs.end() && inputFacet(f, f) &&
	  ftmap.find(f) == ftmap.end()) {
	fm.push_back(f);
	Segment s;
	for (unsigned int j = 0u; j < 3u; ++j) {
	  Point *a = f.vertex(j), *b = f.vertex((j+1u)%3u);
	  if (a < b && inputSegment(a, b, s) && ss.find(s) == ss.end())
	    sp.push_back(Split(s, par, false));
	}
      }
    }
}

void Mesh::insertVertexT (SkinnyQ &sk)
{
  Tetrahedron *t = sk.min().t;
  PTR<Point> c = sk.min().c, p = sk.min().p;
  double r = sk.min().r;
  sk.removeMin();
  Tetrahedrons tp;
  tp.push_back(containingTetrahedron(p, t));
  Tetrahedrons cav = cavityCD(p, tp);
  Facets fs = ::cavityBoundary(cav, ftmap);
  updatePRMap(t, c, p);
  bool bad = false;
  Tetrahedrons nts;
  SkinnyQ nsk;
  for (Facets::iterator f = fs.begin(); !bad && f != fs.end(); ++f) {
    Tetrahedron *u = new Tetrahedron(f->a, f->b, f->c, p);
    nts.push_back(u);
    bad = !checkSkinny(u, nsk) && sliver(u);
  }
  if (bad) {
    setValid(cav, true);
    prmap.erase(p);
    deleteTetrahedrons(nts);
    double d = 0.5*distance(c, t->vertex(0));
    p = randomPointInSphere(c, d);
    sk.insert(Skinny(t, c, p, r));
  }
  else {
    poctree->insert(p);
    removeTetrahedrons(cav, ftmap);
    int st = ts.size();
    for (Tetrahedrons::iterator t = nts.begin(); t != nts.end(); ++t)
      addTetrahedron(ts, ftmap, *t);
    updatePTMap(st);
    while (!nsk.empty()) {
      sk.insert(nsk.min());
      nsk.removeMin();
    }
  }
}

void Mesh::updatePRMap (Tetrahedron *t, Point *c, Point *p)
{
  double d = distance(p, t->vertex(0));
  if (p != c)
    for (unsigned int i = 1u; i < 4u; ++i)
      d = min(d, distance(p, t->vertex(i)));
  prmap.insert(PRPair(p, d));
}

bool Mesh::checkCD (unsigned int s) const // debug
{
  bool res = true;
  for (unsigned int i = s; i < ts.size(); ++i) {
    Tetrahedron *t = ts[i];
    if (t->valid() && !t->ghost())
      for (unsigned int j = 0u; j < 4u; ++j) {
	Facet f = t->facet(j), g(f.c, f.b, f.a);
	if (!inputFacet(f, f) && !inputFacet(g, g)) {
	  Tetrahedron *u = getTetrahedron(g, ftmap);
	  if (t < u && !u->ghost()) {
	    Point *y = u->oppositeVertex(g);
	    int s = InSphere(t->vertex(0), t->vertex(1), t->vertex(2),
			     t->vertex(3), y);
	    if (s > -1) {
	      res = false;
	      int nu = find(ts.begin(), ts.end(), u) - ts.begin();
	      cerr << "bad tets " << i << " and " << nu << endl;
	    }
	  }
	}
      }
  }
  return res;
}

double distance (Point *a, Point *b)
{
  return Distance(a, b).getApproxMid(1e-8).x;
}

void setValid (const Tetrahedrons &ts, bool flag)
{
  for (Tetrahedrons::const_iterator t = ts.begin(); t != ts.end(); ++t)
    (*t)->setValid(flag);
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

bool intersects (Tetrahedron *t, const Segment &s)
{
  bool va = t->isVertex(s.a), vb = t->isVertex(s.b);
  if (va && vb || !va && contains(t, s.a) || !vb && contains(t, s.b))
    return true;
  for (unsigned int i = 0u; i < 4u; ++i) {
    Facet f = t->facet(i);
    PTR<TrianglePlane> p = new TrianglePlane(f.a, f.b, f.c);
    if (intersects(p, s.a, s.b, true) == 1)
      return true;
  }
  return false;
}

int intersects (TrianglePlane *p, Point *t, Point *h, bool strict)
{
  int st = Side(p, t), sh = Side(p, h);
  if (!strict && st == 0 && sh == 0 && overlaps(p, t, h))
    return 0;
  return st*sh == -1 && intersectTriangleLine(p, t, h) == 1 ? 1 : -1;
}

bool overlaps (TrianglePlane *p, Point *t, Point *h)
{
  int ct = contains(p, t), ch = contains(p, h);
  return ct >= 0 || ch >= 0 || ct == -2 && ch == -2 ||
    intersectsEdges(p, t, h);
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

PPoints points (const Facets &fs)
{
  PPoints pts;
  PPointSet done;
  for (Facets::const_iterator f = fs.begin(); f != fs.end(); ++f)
    for (unsigned int i = 0u; i < 3u; ++i) {
      Point *p = f->vertex(i);
      if (done.insert(p).second)
	pts.push_back(p);
    }
  return pts;
}

Tetrahedron * getTetrahedron (const Facet &f, const FTMap &ftmap)
{
  FTMap::const_iterator i = ftmap.find(f);
  return i == ftmap.end() ? 0 : i->second;
}

Tetrahedrons innerTetrahedrons (const FTMap &ftmap, const FacetSet &fs)
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
  setValid(res, true);
  return res;
}

bool containsProjection (Point *a, Point *b, Point *c, Point *q)
{
  PTR<TrianglePlane> p = new TrianglePlane(a, b, c);
  return OrthogonalSide(q, a, b, p) > -1 && OrthogonalSide(q, b, c, p) > -1 &&
    OrthogonalSide(q, c, a, p) > -1;
}

Point * randomPointInSphere (Point *o, double r)
{
  double dx = randomNumber(-1.0, 1.0), dy = randomNumber(-1.0, 1.0),
    dz = randomNumber(-1.0, 1.0),
    d = randomNumber(0.0, r)/sqrt(dx*dx + dy*dy + dz*dz);
  PV3<double> p = o->getApproxMid(1e-8);
  return new SphereCenter(p.x + dx*d, p.y + dy*d, p.z + dz*d);
}

// debug

AffinePoint * afp (Point *p)
{
  return dynamic_cast<AffinePoint *>(p);
}

CircleCenter *ccp (Point *p)
{
  return dynamic_cast<CircleCenter *>(p);
}

SphereCenter *scp (Point *p)
{
  return dynamic_cast<SphereCenter *>(p);
}

void pps (const Points &pts);

void ps (const Segment &s)
{
  Points pts;
  pts.push_back(s.a);
  pts.push_back(s.b);
  pps(pts);
}

void pss (const Segments &ss)
{
  cerr << "(";
  for (Segments::const_iterator s = ss.begin(); s != ss.end(); ++s)
    ps(*s);
  cerr << ")" << endl;
}

void pedge (Edge *, int);

void psegment (const Segment &s, int i)
{
  Vertex t(s.a), h(s.b);
  Edge e(&t, &h);
  pedge(&e, i);
}

Polyhedron * polyhedron (const FacetSet &fs)
{
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (FacetSet::const_iterator f = fs.begin(); f != fs.end(); ++f)
    a->addTriangle(f->vertex(0), f->vertex(1), f->vertex(2), pvmap);
  return a;
}

Polyhedron * polyhedron (const SFMap &sfmap)
{
  FacetSet fs;
  for (SFMap::const_iterator iter = sfmap.begin(); iter != sfmap.end(); ++iter)
    fs.insert(iter->second);
  return polyhedron(fs);
}

Polyhedron * polyhedron (const SFMap &sfmap, unsigned int id)
{
  FacetSet fs;
  for (SFMap::const_iterator iter = sfmap.begin(); iter != sfmap.end(); ++iter)
    if (iter->second.id == id)
      fs.insert(iter->second);
  return polyhedron(fs);
}

int ori (Point *a, Point *b, Point *c, Point *d)
{
  return Orientation(a, b, c, d);
}

SphereCenter * sc (Point *a, Point *b, Point *c, Point *d)
{
  return new SphereCenter(a, b, c, d);
}

int main (int argc, char *argv[])
{
  if (argc < 2)
    return 0;
  ifstream str(argv[1]);
  if (!str.good())
    return 0;
  acp::enable();
  PTR<Polyhedron> a = readPolyhedronVTK(str, false);
  double r = 2.0, cs = cos(10.0*3.14159/180.0), rm = 0.01;
  PTR<Polyhedron> b = mesh(a, r, cs, rm);
  acp::disable();
}

