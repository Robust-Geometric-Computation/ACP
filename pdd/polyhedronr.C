#include "polyhedron.h"

unsigned int nthreads (unsigned int n)
{
#ifdef USETHREADS
  return n;
#else
  return 1;
#endif
}

void Point::getBBox (double *bbox)
{
  PV3<Interval> p = getApprox(1.0);
  for (unsigned int i = 0u; i < 3u; ++i) {
    bbox[2*i] = p[i].l;
    bbox[2*i+1] = p[i].u;
  }
}

bool Point::identical (Point *a)
{
  return order(a) == 0;
}

int Point::order (Point *a)
{
  if (this == a)
    return 0;
  SumPoint *s1 = dynamic_cast<SumPoint *>(this);
  if (s1) {
    SumPoint *s2 = dynamic_cast<SumPoint *>(a);
    if (s2 &&
	(s1->a->identicalI(s2->a) && s1->b->identicalI(s2->b) ||
	 s1->a->identicalI(s2->b) && s1->b->identicalI(s2->a)))
      return 0;
  }
#ifdef RESTEST2
  extern double tAcp, tRes;
  extern unsigned long int nnn;
  timeval t0, t1, t2;
  gettimeofday(&t0, 0);
  int s = PointOrderR(this, a);
  gettimeofday(&t1, 0);
  double ta = diffTime(t1, t0);
  int pointOrderParametric (Point *a, Point *b);
  int sr = pointOrderParametric(this, a);
  gettimeofday(&t2, 0);
  double tr = diffTime(t2, t1);
  if (sr == -2)
    return s;
  if (sr != s) {
    cerr << "bad residue" << endl;
    exit(0);
  }
  tAcp += ta;
  tRes += tr;
  ++nnn;
  return s;
#endif
  return PointOrderR(this, a);
}

bool Point::identicalI (Point *a)
{
  if (!(input && a->input))
    return false;
  PV3<Interval> p = getApprox(1.0), q = a->getApprox(1.0);
  return p.x.l == q.x.l && p.y.l == q.y.l && p.z.l == q.z.l;
}

bool Point::onLine (Point *a, Point *b)
{
  if (this == a || this == b)
    return true;
#ifdef RESTEST3
  extern double tAcp, tRes;
  extern unsigned long int nnn;
  timeval t0, t1, t2;
  gettimeofday(&t0, 0);
  bool res = OnLine(this, a, b) == 0;
  if (!res) return res;
  gettimeofday(&t1, 0);
  double ta = diffTime(t1, t0);
  int onLineParametric (Point *a, Point *b, Point *c);
  int sr = onLineParametric(this, a, b);
  gettimeofday(&t2, 0);
  double tr = diffTime(t2, t1);
  if (sr == -2)
    return res;
  if ((sr == 0) != res) {
    cerr << "bad residue" << endl;
    exit(0);
  }
  tAcp += ta;
  tRes += tr;
  ++nnn;
  return res;
#endif
  return OnLine(this, a, b) == 0;
}

int Point::side (Plane *a)
{
  TrianglePlane *t = dynamic_cast<TrianglePlane *>(a);
  if (t->getA() == this || t->getB() == this || t->getC() == this)
    return 0;
#ifdef RESTEST4
  extern double tAcp, tRes;
  extern unsigned long int nnn;
  timeval t0, t1, t2;
  gettimeofday(&t0, 0);
  int s = Side(a, this);
  if (s) return s;
  gettimeofday(&t1, 0);
  double ta = diffTime(t1, t0);
  int sideParametric (Point *a, Point *b, Point *c, Point *d);
  int sr = sideParametric(this, t->getA(), t->getB(), t->getC());
  gettimeofday(&t2, 0);
  double tr = diffTime(t2, t1);
  if (sr == -2)
    return s;
  if (sr != s) {
    cerr << "bad residue" << endl;
    exit(0);
  }
  tAcp += ta;
  tRes += tr;
  ++nnn;

  if (nnn == 5000) {
    double ta = tAcp/double(nnn), tr = tRes/double(nnn),
      ra = ta == 0.0 ? 0.0 : tr/ta;
  cerr << " predicates: " << nnn << endl << "time per predicate: acp " << ta
       << " residue " << tr << " ratio " << ra << endl;
  exit(0);
  }
  return s;
#endif
  return Side(a, this);
}

PTR<Point> Rdir = new Point(0.8401877171547095, 0.394382926819093, 0.7830992237586059);

bool onEdge (Point *a, Point *t, Point *h, bool strict)
{
  if (a->identical(t) || a->identical(h))
    return !strict;
  return PointOrderE(a, t, h) == 1 && PointOrderE(a, h, t) == 1;
}

bool closerPair (Point *a, Point *b, Point *c, Point *d)
{
  if (a == c && b == d || a == d && b == c)
    return false;
  return CloserPair(a, b, c, d) == 1;
}

bool Plane::coplanar (Plane *p)
{
  if (this == p)
    return true;
  return Coplanar(this, p) == 0;
}

int projectionCoordinate (Plane *p)
{
  PV3<Interval> n = p->getApprox(1.0).n;
  int s[] = {n.x.sign(), n.y.sign(), n.z.sign()};
  if (s[0] && s[1] && s[2]) {
    for (unsigned int i = 0u; i < 3u; ++i)
      if (s[i] == -1) {
	double t = - n[i].l;
	n[i].l = - n[i].u;
	n[i].u = t;
      }
    for (int i = 0; i < 3; ++i)
      if (n[i].l > n[(i+1)%3].u && n[i].l > n[(i+2)%3].u)
	return (i+1)*s[i];
  }
  for (int i = 0; i < 3; ++i)
    if (s[i])
      return (i+1)*s[i];
  for (int i = 0; i < 3; ++i) {
    int si = PlaneNormalCoord(p, i);
    if (si)
      return (i+1)*si;
  }
  exit(0);
}

double bboxSize (double *bb)
{
  return max(bb[1] - bb[0], max(bb[3] - bb[2], bb[5] - bb[4]));
}

void copyBBox (const double *bbf, double *bbt)
{
  for (int i = 0; i < 6; ++i)
    bbt[i] = bbf[i];
}

void mergeBBox (const double *bbf, double *bbt)
{
  for (int i = 0; i < 3; ++i) {
    bbt[2*i] = min(bbt[2*i], bbf[2*i]);
    bbt[2*i+1] = max(bbt[2*i+1], bbf[2*i+1]);
  }
}

bool bboxOverlap (const double *a, const double *b, double s)
{
  for (int i = 0; i < 3; ++i)
    if (a[2*i+1] + s < b[2*i] || b[2*i+1] + s < a[2*i])
      return false;
  return true;
}

bool bboxOverlap (Point *a, const double *bbox)
{
  double abox[6];
  a->getBBox(abox);
  return bboxOverlap(abox, bbox);
}

HEdges Vertex::outgoingHEdges () const
{
  HEdges ed;
  for (Edges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    for (HEdges::iterator f = (*e)->hedges.begin(); f != (*e)->hedges.end(); ++f)
      if ((*f)->tail() == this)
	ed.push_back(*f);
  return ed;
}

Faces Vertex::incidentFaces () const
{
  Faces fa;
  for (Edges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    for (HEdges::iterator f = (*e)->hedges.begin(); f != (*e)->hedges.end(); ++f)
      if ((*f)->tail() == this)
	fa.push_back((*f)->getF());
  return fa;
}

HEdge * Vertex::connected (Vertex *a) const
{
  HEdges ed = outgoingHEdges();
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    if ((*e)->head() == a)
      return *e;
  return 0;
}

Vertex * HEdge::tail () const
{
  return forward ? e->t : e->h;
}

Vertex * HEdge::head () const
{
  return forward ? e->h: e->t;
}

HEdge * HEdge::cw () const
{
  int n = e->hedges.size();
  for (int i = 0; i < n; ++i)
    if (e->hedges[i] == this)
      return e->hedges[(i+1)%n];
  return 0;
}

HEdge * HEdge::ccw () const
{
  int n = e->hedges.size();
  for (int i = 0; i < n; ++i)
    if (e->hedges[i] == this)
      return e->hedges[i == 0 ? n - 1 : i - 1];
  return 0;
}

Vertices HEdge::loop () const
{
  Vertices ve;
  const HEdge *e = this;
  do {
    ve.push_back(e->tail());
    e = e->next;
  }
  while (e != this);
  return ve;
}

Points HEdge::pointLoop () const
{
  Vertices ve = loop();
  Points pts;
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
    pts.push_back((*v)->p);
  return pts;
}

HEdges HEdge::edgeLoop ()
{
  HEdges ed;
  HEdge *e = this;
  do {
    ed.push_back(e);
    e = e->next;
  }
  while (e != this);
  return ed;
}

Edge::Edge (Vertex *t, Vertex *h) : t(t), h(h)
{ 
  hedges.reserve(2);
  setBBox();
}

void Edge::setBBox ()
{
  copyBBox(t->bbox, bbox);
  mergeBBox(h->bbox, bbox);
}

Edge::~Edge ()
{
  for (HEdges::iterator h = hedges.begin(); h != hedges.end(); ++h)
    delete *h;
}

HEdge * Edge::addHEdge (bool forward)
{
  HEdge *e = new HEdge(this, forward);
  hedges.push_back(e);
  return e;
}

void Edge::removeHEdge (HEdge *e)
{
  HEdges::iterator j = remove(hedges.begin(), hedges.end(), e);
  hedges.erase(j, hedges.end());
  delete e;
}

void Edge::sortHEdges ()
{
  if (hedges.size() > 2)
    sort(hedges.begin(), hedges.end(), EdgeOrder(this));
}

bool HFace::pos () const 
{ 
  return this == f->hfaces; 
}
  
HFace * HFace::twin () const 
{
  return this == f->hfaces ? f->hfaces + 1 : f->hfaces;
}

HFaces HFace::neighbors () const
{
  HEdges ed = f->boundaryHEdges();
  HFaces hf;
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    hf.push_back(neighbor(*e));
  return hf;
}

HFace * HFace::neighbor (HEdge *e) const
{
  bool ccw = pos() == e->forward;
  HEdge *f = ccw ? e->ccw() : e->cw();
  bool flag = e->forward == f->forward ? pos() : !pos();
  Face *g = f->f;
  return flag ? g->hfaces + 1 : g->hfaces;
}

Face::Face (HEdge *h, TrianglePlane *pi) : h(h)
{
  p = pi ? pi : new TrianglePlane(h->tail()->p, h->head()->p, h->next->head()->p);
  hfaces[0].f = hfaces[1].f = this;
  addLoop(h, true);
}

Face::Face (Point *a, Point *b, Point *c) : h(0)
{
  p = new TrianglePlane(a, b, c);
  hfaces[0].f = hfaces[1].f = this;
  a->getBBox(bbox);
  double bb[6];
  b->getBBox(bb);
  mergeBBox(bb, bbox);
  c->getBBox(bb);
  mergeBBox(bb, bbox);
}

void Face::addLoop (HEdge *h, bool flag)
{
  HEdges ed = h->edgeLoop();
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    (*e)->f = this;
  if (flag) {
    copyBBox(ed[0]->tail()->bbox, bbox);
    for (int i = 1; i < ed.size(); ++i)
      mergeBBox(ed[i]->tail()->bbox, bbox);
  }
}

void Face::update ()
{
  Vertices ve = h->loop();
  p = new TrianglePlane(ve[0]->p, ve[1]->p, ve[2]->p);
  copyBBox(ve[0]->bbox, bbox);
  for (int i = 1; i < ve.size(); ++i)
    mergeBBox(ve[i]->bbox, bbox);
}

bool Face::boundaryVertex (Point *a) const
{
  Vertices ve = h->loop();
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
    if ((*v)->p->identical(a))
      return true;
  return false;
}

bool Face::boundaryVertex (Vertex *v) const
{
  Vertices ve = h->loop();
  return find(ve.begin(), ve.end(), v) != ve.end();
}

bool Face::boundaryEdge (Edge *e) const
{
  for (HEdges::iterator f = e->hedges.begin(); f != e->hedges.end(); ++f)
    if ((*f)->f == this)
      return true;
  return false;
}

bool Face::sharedEdge (Face *f) const
{
  if (h && f->h) {
    HEdges ed = h->edgeLoop();
    for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
      for (HEdges::iterator i = (*e)->e->hedges.begin(); 
	   i != (*e)->e->hedges.end(); ++i)
	if ((*i)->f == f)
	  return true;
    return false;
  }
  Points p1 = boundaryPoints(), p2 = f->boundaryPoints();
  int n = 0;
  for (Points::iterator p = p1.begin(); p != p1.end(); ++p)
    if (find(p2.begin(), p2.end(), *p) != p2.end())
      ++n;
  return n > 1;
}

void Face::sharedVertices (Face *f, Points &pts) const
{
  Vertices ve = h->loop(), vef = f->h->loop();
  for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
    if (find(vef.begin(), vef.end(), *v) != vef.end())
      pts.push_back((*v)->p);
}

Points Face::boundaryPoints () const
{
  if (h)
    return h->pointLoop();
  Points pts;
  pts.push_back(p->getA());
  pts.push_back(p->getB());
  pts.push_back(p->getC());
  return pts;
}

bool Face::intersects (Face *g, bool strict)
{
  int sf[4], sg[4];
  if (intersectsFP(g, sg) || g->intersectsFP(this, sf) ||
      verifyFP(g, sg) || g->verifyFP(this, sf))
    return false;
  return intersectsFE(g, sg, strict) || g->intersectsFE(this, sf, strict);  
}

bool Face::intersectsFP (Face *g, int *sg)
{
  Points pf = boundaryPoints(), pg = g->boundaryPoints();
  int n = pg.size();
  for (int i = 0; i < n; ++i)
    if (find(pf.begin(), pf.end(), pg[i]) != pf.end())
      sg[i] = 2;
    else {
      PlaneSideO pgi(getP(), pg[i]);
      int s =pgi.getApprox(1.0).x.sign();
      sg[i] = s ? s : 3;
    }
  return checkFP(sg, n);
}

bool Face::checkFP (int *s, int n) const
{
  bool pflag = false, nflag = false;
  int nv = 0;
  for (int i = 0; i < n; ++i)
    if (s[i] == 0 || s[i] == 3)
      return false;
    else if (s[i] == 1)
      pflag = true;
    else if (s[i] == -1)
      nflag = true;
    else
      ++nv;
  return nv < 2 && !(pflag && nflag);
}

bool Face::verifyFP (Face *g, int *sg)
{
  Points pg = g->boundaryPoints();
  int n = pg.size();
  for (int i = 0; i < n; ++i)
    if (sg[i] == 3)
      sg[i] = pg[i]->side(getP());
  return checkFP(sg, n);
}

bool Face::intersectsFE (Face *g, int *sg, bool strict)
{
  Points pg = g->boundaryPoints();
  int n = pg.size();
  if (!strict) {
    bool flag = false;
    for (int i = 0; i < n; ++i)
      if ((sg[i] == 0 || sg[i] == 2) && (sg[(i+1)%n] == 0 || sg[(i+1)%n] == 2)) {
	flag = true;
	if (intersectsFEP(pg[i], pg[(i+1)%n], strict))
	  return true;
      }
    if (flag || sharedEdge(g))
      return false;
    for (int i = 0; i < n; ++i)
      if (sg[i] == 0 && contains(pg[i], strict))
	return true;
  }
  for (int i = 0; i < n; ++i)
    if (sg[i] && sg[i] == - sg[(i+1)%n]) {
      PTR<Point> p = new EPPoint(pg[i], pg[(i+1)%n], getP());
      if (bboxOverlap(p, bbox) && contains(p, strict))
	return true;
    }
  return false;
}

bool Face::intersectsFEP (Point *et, Point *eh, bool strict)
{
  if (contains(et, strict) || contains(eh, strict))
    return true;
  Points pf = boundaryPoints();
  int n = pf.size();
  for (int i = 0; i < n; ++i)
    if (intersectsEE(pf[i], pf[(i+1)%n], et, eh, strict))
      return true;
  return false;
}

bool Face::intersectsEE (Point *et, Point *eh, Point *ft, Point *fh, bool strict)
{
  if (et == ft || et == fh || eh == ft || eh == fh)
    return false;
  int c = p->getPC(), tp1 = LeftTurn(et, ft, fh, c);
  if (!strict && tp1 == 0 && onEdge(et, ft, fh, true))
    return true;
  int tp2 = LeftTurn(eh, ft, fh, c);
  if (!strict && tp2 == 0 && onEdge(eh, ft, fh, true))
    return true;
  if (tp1*tp2 > -1)
    return false;
  int tp3 = LeftTurn(ft, et, eh, c);
  if (!strict && tp3 == 0 && onEdge(ft, et, eh, true))
    return true;
  int tp4 = LeftTurn(fh, et, eh, c);
  return !strict && tp4 == 0 && onEdge(fh, et, eh, true) ||
    tp3*tp4 == -1;
}

bool Face::contains (Point *a, bool strict, int *ie)
{
  if (!bboxOverlap(a, bbox))
    return false;
  Points pts = boundaryPoints();
  int n = pts.size(), c = p->getPC();
  Point *et = pts.back();
  for (int i = 0; i < n; ++i) {
    Point *eh = pts[i];
    if (a == et || a == eh)
      return false;
    int s = LeftTurn(a, et, eh, c);
    if (s == -1)
      return false;
    if (s == 0)
      if (strict)
	return false;
      else if (onEdge(a, et, eh, true)) {
	if (ie)
	  *ie = i == 0 ? n - 1 : i - 1;
      	return true;
      }
      else
	return false;
    et = eh;
  }
  return true;
}

Point * Face::centroid () const
{
  Point *a = h->tail()->p, *b = h->head()->p, *c = h->next->head()->p;
  return new CentroidPoint(a, b, c);
}

PTR<Point> Face::rayIntersection (Point *a, Point *r)
{
  if (a->side(p) == 0)
    return contains(a, false) ? a : 0;
  if (PlaneRayAlignment(p, r) == 0)
    return 0;
  PTR<Point> q = new RayPlanePoint(a, r, p);
  if (PointOrder(a, q, r) == 1 && (contains(q, false) || boundaryVertex(q)))
    return q;
  return 0;
}

void Face::triangulate (Triangles &tr)
{
  Vertices ve = h->loop();
  int n = ve.size();
  for (int i = 1; i + 1 < n; ++i)
    tr.push_back(Triangle(ve[0]->getP(), ve[i]->getP(), ve[i+1]->getP()));
}

Shell::~Shell ()
{
  if (octreef)
    delete octreef;
}

Shell::Shell (const HFaces &hf) : hfaces(hf), c(0)
{
  for (HFaces::iterator h = hfaces.begin(); h != hfaces.end(); ++h)
    (*h)->s = this;
  setBBox();
  setOctree();
}

void Shell::setBBox ()
{
  bbox[0] = bbox[2] = bbox[4] = 1e20;
  bbox[1] = bbox[3] = bbox[5] = -1e20;
  for (HFaces::iterator f = hfaces.begin(); f != hfaces.end(); ++f)
    mergeBBox((*f)->f->bbox, bbox);
}

void Shell::setOctree ()
{
  Faces fa;
  for (HFaces::iterator f = hfaces.begin(); f != hfaces.end(); ++f)
    fa.push_back((*f)->f);
  octreef = Octree<Face *>::octree(fa, bbox);
}

bool Shell::outer () const
{
  PTR<Point> r = new Point(0.0, 0.0, 1.0);
  Vertex *vm = vmax(r);
  HEdge *em = 0;
  HFace *fm = 0;
  for (Edges::iterator e = vm->edges.begin(); e != vm->edges.end(); ++e)
    for (HEdges::iterator h = (*e)->hedges.begin(); h != (*e)->hedges.end(); ++h)
      for (int i = 0; i < 2; ++i)
	if ((*h)->f->hfaces[i].s == this) {
	  if (!em || em->e != *e && SlopeOrder(*e, em->e, r) == 1) {
	    HEdge *nem = *h;
	    HFace *nfm = nem->f->hfaces + i;
	    if (!nfm->f->p->coplanar(nfm->neighbor(nem)->f->p)) {
	      em = nem;
	      fm = nfm;
	      }
	  }
	  break;
	}
  return Convex(em, fm) == 1;
}

Vertex * Shell::vmax (Point *r) const
{
  Vertex *v = hfaces[0]->f->h->tail();
  while (true) {
    Vertex *w = v;
    for (Edges::iterator e = v->edges.begin(); e != v->edges.end(); ++e)
      for (HEdges::iterator h = (*e)->hedges.begin(); h != (*e)->hedges.end(); ++h)
	if ((*h)->f->hfaces[0].s == this || (*h)->f->hfaces[1].s == this) {
	  Vertex *nw = vmax((*h)->f, r);
	  if (w != nw && PointOrder(w->p, nw->p, r) == 1)
	    w = nw;
	}
    if (v == w) {
      double rb[6];
      rayBBox(w->p, r, rb);
      Faces fa;
      octreef->find(rb, fa);
      for (Faces::iterator f = fa.begin(); v == w && f != fa.end(); ++f)
	if ((*f)->rayIntersection(w->p, r) != 0)
	  w = vmax(*f, r);
    }
    if (v == w)
      break;
    v = w;
  }
  return v;
}

Vertex * Shell::vmax (Face *f, Point *r) const
{
  Vertices ve = f->h->loop();
  Vertex *v = ve[0];
  for (int i = 1; i < ve.size(); ++i)
    if (PointOrder(v->p, ve[i]->p, r) == 1)
      v = ve[i];
  return v;
}

bool Shell::contains (Shell *s) const
{
  if (!bboxOverlap(bbox, s->bbox) || s->subset(this))
    return false;
  for (HFaces::const_iterator h = s->hfaces.begin(); h != s->hfaces.end(); ++h) {
    Vertices ve = (*h)->f->h->loop();
    for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v) {
      int i = contains((*v)->p);
      if (i == 1)
	return true;
      if (i == -1)
	return false;
    }
  }
  return false;
}

bool Shell::subset (const Shell *s) const
{
  if (s->hfaces.size() < hfaces.size())
    return false;
  for (HFaces::const_iterator h = hfaces.begin(); h != hfaces.end(); ++h)
    if ((*h)->twin()->s != s)
      return false;
  return true;
}

int Shell::contains (Point *a) const
{
  if (!bboxOverlap(a, bbox))
    return -1;
  PTR<Point> r = new Point(0.0, 0.0, 1.0);
  double rb[6];
  rayBBox(a, r, rb);
  Faces fa;
  octreef->find(rb, fa);
  bool res = false;
  for (Faces::iterator f = fa.begin(); f != fa.end(); ++f)
    if ((*f)->boundaryVertex(a) ||
	a->side((*f)->getP()) == 0 && (*f)->contains(a, false))
      return 0;
    else if ((*f)->rayIntersection(a, r) != 0)
      res = !res;
  return res ? 1 : -1;
}

void Shell::rayBBox (Point *a, Point *r, double *rb) const
{
  a->getBBox(rb);
  RayZPlanePoint q(a, r, bbox[5]);
  double qb[6];
  q.getBBox(qb);
  mergeBBox(qb, rb);
}

int Shell::euler () const
{
  set<Vertex *> vs;
  set<Edge *> es;
  set<Face *> fs;
  for (HFaces::const_iterator f = hfaces.begin(); f != hfaces.end(); ++f)
    fs.insert((*f)->f);
  for (set<Face *>::iterator f = fs.begin(); f != fs.end(); ++f) {
    HEdge *e = (*f)->h;
    do {
      vs.insert(e->tail());
      es.insert(e->e);
      e = e->next;
    }
    while (e != (*f)->h);
  }
  int nv = vs.size(), ne = es.size(), nf = fs.size();
  return nv - ne + nf;
}

bool Shell::pos () const
{
  for (HFaces::const_iterator h = hfaces.begin(); h != hfaces.end(); ++h)
    if (!(*h)->pos())
      return false;
  return true;
}

void deleteShells (const Shells &sh)
{
  for (Shells::const_iterator s = sh.begin(); s != sh.end(); ++s)
    delete *s;
}

bool Cell::contains (Point *p) const
{
  if (outer && outer->contains(p) < 1)
    return false;
  for (Shells::const_iterator s = inner.begin(); s != inner.end(); ++s)
    if ((*s)->contains(p) > -1)
      return false;
  return true;
}

Point * Cell::interiorPoint () const
{
  Shell *s = getShell(0);
  HFace *f = largestHFace(outer->hfaces);
  PTR<Point> p = f->f->centroid(), n = new HFaceNormal(f), qmin = 0;
  for (int i = 0; i < nShells(); ++i) {
    Shell *s = getShell(i);
    for (HFaces::iterator h = s->hfaces.begin(); h != s->hfaces.end(); ++h) {
      Face *g = (*h)->f;
      if (g != f->f) {
	PTR<Point> q = g->rayIntersection(p, n);
	if (q && p != q && (!qmin || CloserPair(p, q, p, qmin) == 1))
	  qmin = q;
      }
    }
  }
  if (qmin)
    return new CentroidPoint(p, qmin);
  return 0;
}

// for chloe
PTR<Point> Cell::interiorPoint (Vertex *a, double d) const
{
  HEdges ed = a->outgoingHEdges();
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    for (int i = 0; i < 2; ++i) {
      HFace *f1 = (*e)->f->hfaces + i;
      if (f1->s->c == this) {
	HFace *f2 = f1->neighbor(*e);
	PTR<Point> p = new RayPoint(a->p, (*e)->head()->p, d),
	  n = new BisectorPoint(f1, f2, *e), qmin = 0;
	for (int j = 0; j < nShells(); ++j) {
	  Shell *s = getShell(j);
	  for (HFaces::iterator h = s->hfaces.begin(); h != s->hfaces.end(); ++h) {
	    Face *g = (*h)->f;
	    if (g != f1->f && g != f2->f) {
	      PTR<Point> q = g->rayIntersection(p, n);
	      if (q && (!qmin || CloserPair(p, q, p, qmin) == 1))
		qmin = q;
	    }
	  }
	}
	if (qmin)
	  return new RayPoint(p, qmin, d);
	else
	  return new RayPoint(p, new SumPoint(p, n), d);
      }
    }
  return 0;
}

HFace * largestHFace (const HFaces &hf)
{
  HFace *f = hf[0];
  double a = areaSquared(f->getF());
  for (int i = 1; i < hf.size(); ++i) {
    double ai = areaSquared(hf[i]->getF());
    if (ai > a) {
      f = hf[i];
      a = ai;
    }
  }
  return f;
}

double areaSquared (Face *f)
{
  Points pts = f->getBoundary()->pointLoop();
  PV3<double> a = pts[0]->getApproxMid(1.0), b = pts[1]->getApproxMid(1.0),
    c = pts[2]->getApproxMid(1.0), d = (c - b).cross(a - b);
  return d.dot(d);
}

Polyhedron::~Polyhedron ()
{
  for (Vertices::iterator v = vertices.begin(); v != vertices.end(); ++v)
    delete *v;
  for (Edges::iterator e = edges.begin(); e != edges.end(); ++e)
    delete *e;
  for (Faces::iterator f = faces.begin(); f != faces.end(); ++f)
    delete *f;
  for (Cells::iterator c = cells.begin(); c != cells.end(); ++c)
    delete *c;
}

Vertex * Polyhedron::getVertex (Point *p)
{
  Vertex *v = new Vertex(p);
  if (vertices.empty())
    copyBBox(v->bbox, bbox);
  else
    mergeBBox(v->bbox, bbox);
  vertices.push_back(v);
  return v;
}

Vertex * Polyhedron::getVertex (Point *p, PVMap &pvmap)
{
  PVMap::iterator iter = pvmap.find(p);
  if (iter != pvmap.end())
    return iter->second;
  Vertex *w = getVertex(p);
  pvmap.insert(PVPair(p, w));
  return w;
}

Edge * Polyhedron::getEdge (Vertex *a, Vertex *b)
{
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if ((*e)->t == a && (*e)->h == b || (*e)->t == b && (*e)->h == a)
      return *e;
  Edge *e = new Edge(a, b);
  edges.push_back(e);
  a->edges.push_back(e);
  b->edges.push_back(e);
  return e;
}

HEdge * Polyhedron::addHEdge (Vertex *a, Vertex *b)
{
  Edge *e = getEdge(a, b);
  bool forward = e->t == a;
  return e->addHEdge(forward);
}

HEdge * Polyhedron::getHEdge (Vertex *a, Vertex *b)
{
  Edge *e = getEdge(a, b);
  bool forward = e->t == a;
  for (HEdges::iterator f = e->hedges.begin(); f != e->hedges.end(); ++f)
    if ((*f)->forward == forward && !(*f)->f)
      return *f;
  return e->addHEdge(forward);
}

Face * Polyhedron::addTriangle (Vertex *a, Vertex *b, Vertex *c, TrianglePlane *p)
{
  HEdge *u = addHEdge(a, b), *v = addHEdge(b, c), *w = addHEdge(c, a);
  u->setNext(v);
  v->setNext(w);
  w->setNext(u);
  Face *f = new Face(u, p);
  faces.push_back(f);
  return f;
}

Face * Polyhedron::addRectangle (Vertex *a, Vertex *b, Vertex *c, Vertex *d)
{
  HEdge *u = addHEdge(a, b), *v = addHEdge(b, c), *w = addHEdge(c, d), 
    *x = addHEdge(d, a);
  u->setNext(v);
  v->setNext(w);
  w->setNext(x);
  x->setNext(u);
  Face *f = new Face(u);
  faces.push_back(f);
  return f;
}

Face * Polyhedron::addFace (HEdge *h)
{
  Face *f = new Face(h);
  faces.push_back(f);
  return f;
}

void Polyhedron::formCells () {
  Shells shells;
  formShells(shells);
  formCellsAux(shells);
}

void Polyhedron::formCellsAux (const Shells &shells)
{
  Shells inner;
  cells.push_back(new Cell(0));
  bool flag = false;
  for (Shells::const_iterator s = shells.begin(); s != shells.end(); ++s)
    if ((*s)->hfaces.size() < 3) {
      flag = true;
      removeShell(*s);
    }
    else {
      (*s)->setBBox();
      (*s)->setOctree();
      if ((*s)->outer())
	cells.push_back(new Cell(*s));
      else
	inner.push_back(*s);
    }
  if (flag)
    removeNullFaces();
  Octree<Cell *> *octreec = cellOctree();
  for (Shells::iterator s = inner.begin(); s != inner.end(); ++s)
    enclosingCell(*s, octreec)->addInner(*s);
  delete octreec;
}

void Polyhedron::removeShell (Shell *s)
{
  for (HFaces::iterator h = s->hfaces.begin(); h != s->hfaces.end(); ++h) {
    HEdge *e = (*h)->f->h;
    if (e)
      removeLoop(e);
  }
  delete s;
}

void Polyhedron::formShells (Shells &shells)
{
  for (Edges::iterator e = edges.begin(); e != edges.end(); ++e)
    (*e)->sortHEdges();
  for (Faces::iterator f = faces.begin(); f != faces.end(); ++f)
    for (int i = 0; i < 2; ++i)
      if (!(*f)->hfaces[i].s)
	shells.push_back(formShell((*f)->hfaces + i));
}

Shell * Polyhedron::formShell (HFace *f) const
{
  Shell *s = new Shell;
  HFaces st;
  st.push_back(f);
  while (!st.empty()) {
    HFace *g = *(st.end()-1);
    st.pop_back();
    if (!g->s) {
      g->s = s;
      s->hfaces.push_back(g);
      HFaces gn = g->neighbors();
      st.insert(st.end(), gn.begin(), gn.end());
    }
  }
  return s;
}

Cell * Polyhedron::enclosingCell (Shell *s, Octree<Cell *> *octreec) const
{
  Cells ce;
  octreec->find(s->hfaces[0]->f->h->tail()->bbox, ce);
  Cell *c = 0;
  for (Cells::iterator d = ce.begin(); d != ce.end(); ++d)
    if ((*d)->outer->contains(s) &&
      (!c || c->outer->contains((*d)->outer)))
      c = *d;
  return c ? c : cells[0];
}

void Polyhedron::clearCells ()
{
  for (Faces::iterator f = faces.begin(); f != faces.end(); ++f)
    for (int i = 0; i < 2; ++i)
      (*f)->hfaces[i].s = 0;
  for (Cells::iterator c = cells.begin(); c != cells.end(); ++c)
    delete *c;
  cells.clear();
}

Face * Polyhedron::addTriangle (Point *a, Point *b, Point *c, PVMap &pvmap,
				TrianglePlane *p)
{
  Vertex *u = getVertex(a, pvmap), *v = getVertex(b, pvmap),
    *w = getVertex(c, pvmap);
  return addTriangle(u, v, w, p);
}

Face * Polyhedron::addTriangle (Face *f, PVMap &pvmap, TrianglePlane *p)
{
  Points pf = f->h->pointLoop();
  return addTriangle(pf[0], pf[1], pf[2], pvmap, p);
}

Polyhedron * Polyhedron::copy () const
{
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f)
    a->addTriangle(*f, pvmap, (*f)->p);
  return a;
}

Polyhedron * Polyhedron::scale (double unit) const
{
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (Vertices::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
    pvmap.insert(PVPair((*v)->p, a->getVertex(new ScalePoint((*v)->p, unit))));
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f)
    a->addTriangle(*f, pvmap);
  return a;
}

Polyhedron * Polyhedron::negative () const
{
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (Vertices::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
    pvmap.insert(PVPair((*v)->p, a->getVertex(new NegPoint((*v)->p))));
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f) {
    Points pf = (*f)->h->pointLoop();
    a->addTriangle(pf[2], pf[1], pf[0], pvmap);
  }
  return a;
}

Polyhedron * Polyhedron::translate (Point *t) const
{
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (Vertices::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
    pvmap.insert(PVPair((*v)->p, a->getVertex(new SumPoint(t, (*v)->p))));
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f)
    a->addTriangle(*f, pvmap);
  return a;
}

Polyhedron * Polyhedron::negativeTranslate (Point *t) const
{
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (Vertices::const_iterator v = vertices.begin(); v != vertices.end(); ++v)
    pvmap.insert(PVPair((*v)->p, a->getVertex(new DiffPoint(t, (*v)->p))));
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f) {
    Points pf = (*f)->h->pointLoop();
    a->addTriangle(pf[2], pf[1], pf[0], pvmap);
  }
  return a;
}

bool Polyhedron::intersects (Polyhedron *a, bool strict, Octree<Face *> *aoctree) const
{
  Octree<Face *> *octree = aoctree ? aoctree : a->faceOctree();
  bool res = contains(a->vertices[0]->p) || a->contains(vertices[0]->p) ||
    intersectsEdges(octree, strict);
  if (!aoctree)
    delete octree;
  return res;
}

bool Polyhedron::contains (Point *p) const
{
  if (!bboxOverlap(p, bbox))
    return false;
  for (int i = 1; i < cells.size(); ++i)
    if (cells[i]->wn == 1 && cells[i]->contains(p))
      return true;
  return false;
}

int Polyhedron::containingCell (Point *p) const
{
  for (int i = 0; i< cells.size(); i++)
    if (cells[i]->contains(p))
      return i;
  return -1;
}

bool Polyhedron::intersectsEdges (Octree<Face *> *octree, bool strict) const
{
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f) {
    Faces fa;
    octree->find((*f)->bbox, fa);
    for (Faces::iterator g = fa.begin(); g != fa.end(); ++g)
      if ((*f)->intersects(*g, strict))
	return true;
  }
  return false;
}

Polyhedron * Polyhedron::boolean (Polyhedron *a, SetOp op)
{
  Polyhedron *poly[] = {this, a}, *c = overlay(poly, 2);
  computeWindingNumbers();
  a->computeWindingNumbers();
  c->formCells();
  set<Cell *> cin;
  for (int i = 1; i < c->cells.size(); ++i) {
    Cell *ci = c->cells[i];
    PTR<Point> p = ci->interiorPoint();
    bool ina = contains(p), inb = a->contains(p);
    if (inSet(ina, inb, op))
      cin.insert(ci);
  }
  Polyhedron *d = c->mergeCells(cin);
  delete c;
  return d;
}

Polyhedron * Polyhedron::mergeCells (const set<Cell *> &cin) const
{
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  for (Faces::const_iterator f = faces.begin(); f != faces.end(); ++f) {
    bool in1 = cin.find((*f)->hfaces[0].s->c) != cin.end(),
      in2 = cin.find((*f)->hfaces[1].s->c) != cin.end();
    if (in1 != in2) {
      Points pf = (*f)->h->pointLoop();
      if (in2)
	a->addTriangle(pf[0], pf[1], pf[2], pvmap, (*f)->p);
      else
	a->addTriangle(pf[2], pf[1], pf[0], pvmap);
    }
  }
  return a;
}

Polyhedron * Polyhedron::selfUnion ()
{
  Polyhedron *a = subdivide(this, false), *b = triangulate(a);
  b->computeWindingNumbers();
  set<Cell *> bin;
  for (int i = 1; i < b->cells.size(); ++i)
    if (b->cells[i]->wn > 0)
      bin.insert(b->cells[i]);
  Polyhedron *c = b->mergeCells(bin);
  delete a;
  delete b;
  return c;
}

// used in simplify from here on

Polyhedron * Polyhedron::cellPolyhedron (int i) const
{
  Cell *c = cells[i];
  Polyhedron *a = new Polyhedron;
  PVMap pvmap;
  if (c->outer)
    a->addHFaces(c->outer->hfaces, pvmap);
  for (Shells::const_iterator s = c->inner.begin(); s != c->inner.end(); ++s)
    a->addHFaces((*s)->hfaces, pvmap);
  a->formCells();
  return a;  
}

void Polyhedron::addHFaces (const HFaces &hf, PVMap &pvmap)
{
  for (HFaces::const_iterator h = hf.begin(); h != hf.end(); ++h)
    addTriangle((*h)->f, pvmap, (*h)->f->p);
}

void Polyhedron::replaceVertex (Face *f, Vertex *v, Vertex *w)
{
  HEdge *e = f->h;
  while (e->next->head() != v)
    e = e->next;
  removeHEdge(e->next->next);
  removeHEdge(e->next);
  Vertex *t = e->tail(), *h = e->head();
  HEdge *en = getHEdge(h, w), *enn = getHEdge(w, t);
  en->f = enn->f = f;
  e->next = en;
  en->next = enn;
  enn->next = e;
  f->h = e;
  f->p = new TrianglePlane(t->p, h->p, w->p);
  copyBBox(e->e->bbox, f->bbox);
  mergeBBox(en->e->bbox, f->bbox);
  mergeBBox(enn->e->bbox, f->bbox);
}

void Polyhedron::removeLoop (HEdge *e)
{
  if (e->f)
    e->f->h = 0;
  HEdges ed = e->edgeLoop();
  for (HEdges::iterator h = ed.begin(); h != ed.end(); ++h)
    removeHEdge(*h);
}

void Polyhedron::removeHEdge (HEdge *h)
{
  Edge *e = h->e;
  e->removeHEdge(h);
  if (e->hedges.empty()) {
    Edges::iterator k = remove(e->t->edges.begin(), e->t->edges.end(), e);
    e->t->edges.erase(k, e->t->edges.end());
    k = remove(e->h->edges.begin(), e->h->edges.end(), e);
    e->h->edges.erase(k, e->h->edges.end());
  }
}

void Polyhedron::moveVertex (Vertex *v, Point *p)
{
  v->p = p;
  p->getBBox(v->bbox);
  mergeBBox(v->bbox, bbox);
  HEdges ed = v->outgoingHEdges();
  for (HEdges::iterator h = ed.begin(); h != ed.end(); ++h) {
    (*h)->e->setBBox();
    if ((*h)->f)
      (*h)->f->update();
    else
      addFace(*h);
  }
}

void Polyhedron::removeNullFaces ()
{
  int i = 0;
  while (i < vertices.size())
    if (vertices[i]->edges.empty()) {
      delete vertices[i];
      vertices[i] = vertices.back();
      vertices.pop_back();
    }
    else
      ++i;
  i = 0;
  while (i < edges.size())
    if (edges[i]->hedges.empty()) {
      delete edges[i];
      edges[i] = edges.back();
      edges.pop_back();
    }
    else
      ++i;
  i = 0;
  while (i < faces.size())
    if (!faces[i]->h) {
      delete faces[i];
      faces[i] = faces.back();
      faces.pop_back();
    }
    else
      ++i;
}

Octree<Face *> * Polyhedron::faceOctree (double s) const
{
  return Octree<Face *>::octree(faces, bbox, s);
}

Octree<Cell *> * Polyhedron::cellOctree () const
{
  Cells ce;
  for (int i = 1; i < cells.size(); ++i)
    ce.push_back(cells[i]);
  return Octree<Cell *>::octree(ce, bbox);
}

void Polyhedron::computeWindingNumbers ()
{
  if (cells.empty())
    formCells();
  cells[0]->wn = 0;
  HFaces st;
  updateWN(cells[0], st);
  set<Cell *> done;
  while (!st.empty()) {
    HFace *f = st.back();
    st.pop_back();
    if (done.insert(f->s->c).second) {
      f->s->c->wn = f->twin()->s->c->wn + (f->pos() ? -1 : 1);
      updateWN(f->s->c, st);
    }
  }
}

void Polyhedron::updateWN (Cell *c, HFaces &st) const
{
  for (int i = 0; i < c->nShells(); ++i) {
    Shell *s = c->getShell(i);
    for (HFaces::iterator h = s->hfaces.begin(); h != s->hfaces.end(); ++h)
      st.push_back((*h)->twin());
  }
}

void Polyhedron::describe ()
{
  if (cells.empty())
    formCells();
  bool wflag = false;
  for (int i = 1; !wflag && i < cells.size(); ++i)
    wflag = cells[i]->wn;
  Cell *c = cells[0];
  cerr << "unbounded cell: " << c->inner.size() << " shells: ";
  for (Shells::iterator s = c->inner.begin(); s != c->inner.end(); ++s)
    cerr << (*s)->hfaces.size() << " hfaces, euler = " << (*s)->euler() << "; ";
  cerr << endl;
  for (int i = 1; i < cells.size(); ++i) {
    Cell *c = cells[i];
    cerr << "bounded cell " << i << ": ";
    if (wflag)
      cerr << "wn = " << c->wn << "; ";
    cerr << c->nShells() << " shells; ";
    for (int  j = 0; j < c->nShells(); ++j) {
      Shell *s = c->getShell(j);
      cerr << s->hfaces.size() << " hfaces, euler = " << s->euler() << "; ";
    }
    cerr << endl;
  }
}

Face * faceVertices (Vertex *a, Vertex *b, Vertex *c)
{
  HEdges ea = a->outgoingHEdges();
  for (HEdges::iterator u = ea.begin(); u != ea.end(); ++u)
    if ((*u)->head() == b) {
      Face *f = (*u)->getF();
      if (f) {
	HEdge *v = (*u)->getNext();
	if (v->head() == c && v->getF() == f) {
	  HEdge *w = v->getNext();
	  if (w->head() == a && w->getF() == f)
	    return f;
	}
      }
    }
  return 0;
}

Face * faceVertices (Vertex *a, Vertex *b, Vertex *c, Vertex *d)
{
  HEdges ea = a->outgoingHEdges();
  for (HEdges::iterator u = ea.begin(); u != ea.end(); ++u)
    if ((*u)->head() == b) {
      Face *f = (*u)->getF();
      if (f) {
	HEdge *v = (*u)->getNext();
	if (v->head() == c && v->getF() == f) {
	  HEdge *w = v->getNext();
	  if (w->head() == d && w->getF() == f) {
	    HEdge *x = w->getNext();
	    if (x->head() == a && x->getF() == f)
	      return f;
	  }
	}
      }
    }
  return 0;
}

Polyhedron * subdivide (Polyhedron *a, bool oneway)
{
  map<Edge *, Points * > epsmap;
  map<pair<Face *, Face *>, FFE> ffemap;
  vector<pair<Face *, Edge *>> fe;
  vector<pair<Edge *, Edge *>> ee;
  intersectFF(a, epsmap, ffemap, fe, ee);
  map<Face *, vector<FFE>> femap = feMap(epsmap, ffemap, fe);
  vector<pair<Points *, Points *>> col;
  intersectFFF(ffemap, femap, col);
  map<PTR<Point>, PTR<Point>> pmap = pMap(epsmap, ee, col);
  set<Point *> bpts;
  vector<pair<PTR<Point>, PTR<Point>>> equiv;
  subedgesE(epsmap, oneway, bpts, equiv);
  subedgesFF(ffemap, oneway, bpts, equiv);
  update(equiv, pmap);
  Polyhedron *b = subfaces(a, oneway, epsmap, femap, pmap);
  deleteMaps(epsmap, ffemap);
  return b;
}

void intersectFF (Polyhedron *a, map<Edge *, Points *> &epsmap,
		  map<pair<Face *, Face *>, FFE> &ffemap,
		  vector<pair<Face *, Edge *>> &fe,
		  vector<pair<Edge *, Edge *>> &ee)
{
  Octree<Face *> *octree = a->faceOctree();
  vector<pair<Face *, Face *>> ff1, ff;
  octree->pairs(ff1);
  delete octree;
  map<pair<Face *, Edge *>, PTR<Point>> fepmap;
  map<pair<Face *, Face *>, Points> ffpsmap;
  intersectFE(ff1, fepmap, epsmap, ffpsmap, ff, fe, ee);
  intersectFF(ff, fepmap, epsmap, ffpsmap, ffemap);
}

void intersectFE (const vector<pair<Face *, Face *>> &ff1,
		  map<pair<Face *, Edge *>, PTR<Point>> &fepmap,
		  map<Edge *, Points *> &epsmap,
		  map<pair<Face *, Face *>, Points> &ffpsmap,
		  vector<pair<Face *, Face *>> &ff,
		  vector<pair<Face *, Edge *>> &fe,
		  vector<pair<Edge *, Edge *>> &ee)
{
  const unsigned int n = nthreads(8);
  unsigned int k = ff1.size(), m = k/n, is = 0;
  FEData fed[n];
  for (int i = 0; i < n; ++i) {
    fed[i].i = i;
    fed[i].is = is;
    is = i + 1 == n ? k : is + m;
    fed[i].ie = is;
    fed[i].ff1 = &ff1;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) intersectFET, (void *) (fed + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  for (int i = 0; i < n; ++i) {
    for (vector<pair<Edge *, PTR<Point>>>::iterator x = fed[i].ep.begin();
	 x != fed[i].ep.end(); ++x)
      update(x->first, x->second, epsmap);
    for (vector<pair<pair<Face *, Edge *>, PTR<Point>>>::iterator
	   x = fed[i].fep.begin(); x != fed[i].fep.end(); ++x)
      fepmap.insert(*x);
    for (vector<pair<pair<Face *, Face *>, PTR<Point>>>::iterator x = fed[i].ffp.begin();
	 x != fed[i].ffp.end(); ++x)
      update(x->first, x->second, ffpsmap);
    ff.insert(ff.end(), fed[i].ff.begin(), fed[i].ff.end());
    fe.insert(fe.end(), fed[i].fe.begin(), fed[i].fe.end());
    ee.insert(ee.end(), fed[i].ee.begin(), fed[i].ee.end());
  }
}

void intersectFET (void *ptr)
{
  FEData *fed = (FEData *) ptr;
  BaseObject::addThread(fed->i);
  for (int i = fed->is; i < fed->ie; ++i) {
    const pair<Face *, Face *> &ff = fed->ff1->at(i);
    intersectFE(ff.first, ff.second, fed->ep, fed->fep, fed->ffp, fed->ff,
		fed->fe, fed->ee);
  }
}

void intersectFE (Face *f, Face *g,
		  vector<pair<Edge *, PTR<Point>>> &ep,
		  vector<pair<pair<Face *, Edge *>, PTR<Point>>> &fep,
		  vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp,
		  vector<pair<Face *, Face *>> &ff,
		  vector<pair<Face *, Edge *>> &fe,
		  vector<pair<Edge *, Edge *>> &ee)
{
  int sf[] = {0, 0, 0, 0}, sg[] = {0, 0, 0, 0};
  if (f->intersectsFP(g, sg) || g->intersectsFP(f, sf) ||
      f->verifyFP(g, sg) || g->verifyFP(f, sf))
    return;
  intersectFE(f, g, sg, ep, fep, ffp, fe, ee);
  intersectFE(g, f, sf, ep, fep, ffp, fe, ee);
  if (signChange(sf) && signChange(sg))
    ff.push_back(pair<Face *, Face *>(f, g));
}

void intersectFE (Face *f, Face *g, int *sg,
		  vector<pair<Edge *, PTR<Point>>> &ep,
		  vector<pair<pair<Face *, Edge *>, PTR<Point>>> &fep,
		  vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp,
		  vector<pair<Face *, Edge *>> &fe,
		  vector<pair<Edge *, Edge *>> &ee)
{
  HEdges eg = g->getBoundary()->edgeLoop();
  int n = eg.size();
  bool flag = false;
  for (int i = 0; i < n; ++i)
    if ((sg[i] == 0 || sg[i] == 2) && (sg[(i+1)%n] == 0 || sg[(i+1)%n] == 2)) {
      intersectFEP(f, eg[i], ffp, ep, fe, ee);
      flag = true;
    }
  if (flag || f->sharedEdge(g))
    return;
  for (int i = 0; i < n; ++i)
    if (sg[i] == 0)
      intersectFV(f, eg[i], ffp, ep);
    else if (sg[i] == - sg[(i+1)%n])
      intersectFEG(f, eg[i], ffp, ep, fep);
}

void intersectFEP (Face *f, HEdge *h,
		   vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp,
		   vector<pair<Edge *, PTR<Point>>> &ep,
		   vector<pair<Face *, Edge *>> &fe,
		   vector<pair<Edge *, Edge *>> &ee)
{
  if (!first(h))
    return;
  Edge *e = h->getE();
  if (f->boundaryEdge(e) || !bboxOverlap(f->getBBox(), e->getBBox()))
    return;
  HEdges ef = f->getBoundary()->edgeLoop();
  for (HEdges::iterator i = ef.begin(); i != ef.end(); ++i)
    if (first(*i) && bboxOverlap(e->getBBox(), (*i)->getE()->getBBox()) &&
	intersectEE(e, (*i)->getE(), f->getP()->getPC(), ffp, ep))
      ee.push_back(pair<Edge *, Edge *>(e, (*i)->getE()));
  fe.push_back(pair<Face *, Edge *>(f, e));
}

bool intersectEE (Edge *e, Edge *f, int c,
		  vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp,
		  vector<pair<Edge *, PTR<Point>>> &ep)
{
  Points pe, pf;
  bool res = intersectEE(e->getT()->getP(), e->getH()->getP(),
			 f->getT()->getP(), f->getH()->getP(), c, pe, pf);
  for (Points::iterator p = pe.begin(); p != pe.end(); ++p) {
    ep.push_back(pair<Edge *, PTR<Point>>(e, *p));
    updateFFP(e, f, *p, ffp);
  }
  for (Points::iterator p = pf.begin(); p != pf.end(); ++p) {
    ep.push_back(pair<Edge *, PTR<Point>>(f, *p));
    updateFFP(e, f, *p, ffp);
  }
  return res;
}

bool intersectEE (Point *et, Point *eh, Point *ft, Point *fh,
		  int c, Points &pe, Points &pf)
{
  int tp1 = et == ft || et == fh ? 0 : LeftTurn(et, ft, fh, c),
     tp2 = eh == ft || eh == fh ? 0 : LeftTurn(eh, ft, fh, c);
  if (tp1*tp2 == 1)
    return false;
  if (tp1 == 0 && onEdge(et, ft, fh, true))
    pf.push_back(et);
  if (tp2 == 0 && onEdge(eh, ft, fh, true))
    pf.push_back(eh);
  int tp3 = ft == et || ft == eh ? 0 : LeftTurn(ft, et, eh, c),
    tp4 = fh == et || fh == eh ? 0 : LeftTurn(fh, et, eh, c);
  if (tp3 == 0 && onEdge(ft, et, eh, true))
    pe.push_back(ft);
  if (tp4 == 0 && onEdge(fh, et, eh, true))
    pe.push_back(fh);
  if (tp1 == 0 && tp2 == 0)
    return !(pe.empty() && pf.empty());
  if (tp1*tp2 == -1 && tp3*tp4 == -1) {
    PTR<Point> p = new EEPoint(et, eh, ft, fh, c);
    pe.push_back(p);
    pf.push_back(p);
    }
  return false;
}

void intersectFV (Face *f, HEdge *h,
		  vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp,
		  vector<pair<Edge *, PTR<Point>>> &ep)
{
  int ie = -1;
  Point *p = h->tail()->getP();
  if (f->contains(p, false, &ie)) {
    if (ie == -1)
      updateFFP(f, h->tail(), ffp);
    else {
      HEdges hef = f->getBoundary()->edgeLoop();
      Edge *ef = hef[ie]->getE();
      ep.push_back(pair<Edge *, PTR<Point>>(ef, p));
      updateFFP(ef, h->getE(), p, ffp);
    }
  }
}

void updateFFP (Face *f, Vertex *v,
		vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp)
{
  HEdges ed = v->outgoingHEdges();
  for (HEdges::iterator e = ed.begin(); e != ed.end(); ++e) {
    Face *g = (*e)->getF();
    pair<Face *, Face *> ff(f < g ? f : g, f < g ? g : f);
    ffp.push_back(pair<pair<Face *, Face *>, PTR<Point>>(ff, v->getP()));
  }
}

void updateFFP (Edge *e1, Edge *e2, Point *p,
		vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp)
{
  for (int i = 0; i < e1->HEdgesN(); ++i) {
    Face *f = e1->getHEdge(i)->getF();
    for (int j = 0; j < e2->HEdgesN(); ++j) {
      Face *g = e2->getHEdge(j)->getF();
      pair<Face *, Face *> ff(f < g ? f : g, f < g ? g : f);
      ffp.push_back(pair<pair<Face *, Face *>, PTR<Point>>(ff, p));
    }
  }
}

void intersectFEG (Face *f, HEdge *h,
		   vector<pair<pair<Face *, Face *>, PTR<Point>>> &ffp,
		   vector<pair<Edge *, PTR<Point>>> &ep,
		   vector<pair<pair<Face *, Edge *>, PTR<Point>>> &fep)
{
  if (!first(h))
    return;
  Edge *e = h->getE();
  if (!bboxOverlap(f->getBBox(), e->getBBox()))
    return;
  PTR<Point> p = new EPPoint(e->getT()->getP(), e->getH()->getP(), f->getP());
  int ie = -1;
  if (f->contains(p, false, &ie)) {
    ep.push_back(pair<Edge *, PTR<Point>>(e, p));
    if (ie == -1) {
      pair<Face *, Edge *> fe(f, e);
      fep.push_back(pair<pair<Face *, Edge *>, PTR<Point>>(fe, p));
    }
    else {
      HEdges hef = f->getBoundary()->edgeLoop();
      Edge *ef = hef[ie]->getE();
      ep.push_back(pair<Edge *, PTR<Point>>(ef, p));
      updateFFP(ef, h->getE(), p, ffp);
    }
  }
}

bool first (HEdge *h)
{
  return h == h->getE()->getHEdge(0);
}

bool signChange (int *s)
{
  bool pflag = false, nflag = false;
  for (int i = 0; i < 4; ++i)
    if (s[i] == 1)
      pflag = true;
    else if (s[i] == -1)
      nflag = true;
  return pflag && nflag;
}

void update (Edge *e, Point *p, map<Edge *, Points *> &epsmap)
{
  map<Edge *, Points *>::iterator i = epsmap.find(e);
  if (i == epsmap.end()) {
    Points *pts = new Points;
    pts->push_back(p);
    epsmap.insert(pair<Edge *, Points *>(e, pts));
  }
  else
    i->second->push_back(p);
}

void update (const pair<Face *, Face *> &ff, Point *p,
	     map<pair<Face *, Face *>, Points> &fppmap)
{
  map<pair<Face *, Face *>, Points>::iterator i = fppmap.find(ff);
  if (i == fppmap.end()) {
    Points pts;
    pts.push_back(p);
    fppmap.insert(pair<pair<Face *, Face *>, Points>(ff, pts));
  }
  else
    i->second.push_back(p);
}

void intersectFF (const vector<pair<Face *, Face *>> &ff,
		  const map<pair<Face *, Edge *>, PTR<Point>> &fepmap,
		  const map<Edge *, Points *> &epsmap,
		  map<pair<Face *, Face *>, Points> &ffpsmap,
		  map<pair<Face *, Face *>, FFE> &ffemap)
{
  const unsigned int n = nthreads(8);
  unsigned int k = ff.size(), m = k/n, is = 0;
  FFData ffd[n];
  for (int i = 0; i < n; ++i) {
    ffd[i].i = i;
    ffd[i].is = is;
    is = i + 1 == n ? k : is + m;
    ffd[i].ie = is;
    ffd[i].ff = &ff;
    ffd[i].fepmap = &fepmap;
    ffd[i].epsmap = &epsmap;
    ffd[i].ffpsmap = &ffpsmap;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) intersectFFT, (void *) (ffd + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  vector<pair< pair<Face *, Face *>, pair<PTR<Point>, PTR<Point>>>> ffpp;
  for (int i = 0; i < n; ++i)
    for (vector<pair< pair<Face *, Face *>, pair<PTR<Point>, PTR<Point>>>>::iterator
	 j = ffd[i].ffpp.begin(); j != ffd[i].ffpp.end(); ++j)
    update(j->first.first, j->first.second, j->second.first, j->second.second,
	   ffemap);
}

void intersectFFT (void *ptr)
{
  FFData *ffd = (FFData *) ptr;
  BaseObject::addThread(ffd->i);
  for (int i = ffd->is; i < ffd->ie; ++i) {
    const pair<Face *, Face *> &ff = ffd->ff->at(i);
    intersectFF(ff.first, ff.second, *ffd->fepmap, *ffd->epsmap,
		*ffd->ffpsmap, ffd->ffpp);
  }
}

void intersectFF (Face *f, Face *g,
		  const map<pair<Face *, Edge *>, PTR<Point>> &fepmap,
		  const map<Edge *, Points *> &epsmap,
		  const map<pair<Face *, Face *>, Points> &ffpsmap,
		  vector<pair< pair<Face *, Face *>, pair<PTR<Point>, PTR<Point>>>> &ffpp)
{
  Points pts;
  findFE(f, g, fepmap, pts);
  findFE(g, f, fepmap, pts);
  if (pts.size() < 2)
    f->sharedVertices(g, pts);
  if (pts.size() < 2)
    sharedVertices(f, g, ffpsmap, pts);
  if (pts.size() != 2)
    return;
  Face *ff = f < g ? f : g, *gg = f < g ? g : f;
  bool flag = FFOrder(ff, gg, pts[0], pts[1]) == 1;
  pair<Face *, Face *> fg(ff, gg);
  pair<PTR<Point>, PTR<Point>> pq(pts[flag ? 0 : 1], pts[flag ? 1 : 0]);
  ffpp.push_back(pair<pair<Face *, Face *>, pair<PTR<Point>, PTR<Point>>>(fg, pq));
}

void findFE (Face *f, Face *g,
	     const map<pair<Face *, Edge *>, PTR<Point>> &fepmap,
	     Points &pts)
{
  HEdges eg = g->getBoundary()->edgeLoop();
  for (HEdges::iterator h = eg.begin(); h != eg.end(); ++h) {
    pair<Face *, Edge *> fe(f, (*h)->getE());
    map<pair<Face *, Edge *>, PTR<Point>>::const_iterator i = fepmap.find(fe);
    if (i != fepmap.end())
      pts.push_back(i->second);
  }
}

void sharedVertices (Face *f, Face *g,
		     const map<pair<Face *, Face *>, Points> &ffpsmap, Points &pts)
{
  pair<Face *, Face *> ff(f < g ? f : g, f < g ? g : f);
  map<pair<Face *, Face *>, Points>::const_iterator i = ffpsmap.find(ff);
  if (i != ffpsmap.end()) {
    pts.insert(pts.end(), i->second.begin(), i->second.end());
    removeDuplicates(pts);
  }
}

void removeDuplicates (Points &pts)
{
  for (int i = 0; i + 1 < pts.size(); ++i) {
    int j = i + 1;
    while (j < pts.size()) {
      if (pts[i]->identical(pts[j])) {
	pts[j] = pts.back();
	pts.pop_back();
      }
      else
	++j;
    }
  }
}

void update (Face *f, Face *g, Point *a, Point *b,
	     map<pair<Face *, Face *>, FFE> &ffemap)
{
  pair<Face *, Face *> fg(f, g);
  FFE ffe(f, g, a, b); 
  ffemap.insert(pair<pair<Face *, Face *>, FFE>(fg, ffe));
}

map<Face *, vector<FFE>> feMap (map<Edge *, Points *> &epsmap,
				map<pair<Face *, Face *>, FFE> &ffemap,
				const vector<pair<Face *, Edge *>> &fe)
{
  map<Face *, vector<FFE>> femap;
  for (map<pair<Face *, Face *>, FFE>::iterator i = ffemap.begin();
       i != ffemap.end(); ++i) {
    update(i->first.first, i->second, femap);
    update(i->first.second, i->second, femap);
  }
  for (vector<pair<Face *, Edge *>>::const_iterator i = fe.begin(); i != fe.end(); ++i)
    feMapFE(epsmap, ffemap, i->first, i->second, femap);
  return femap;
}

void feMapFE (map<Edge *, Points *> &epsmap,
	      map<pair<Face *, Face *>, FFE> &ffemap,
	      Face *f, Edge *e, map<Face *, vector<FFE>> &femap)
{
  HEdges ef = f->getBoundary()->edgeLoop();
  for (HEdges::iterator h = ef.begin(); h != ef.end(); ++h)
    if ((*h)->tail()->getP()->onLine(e->getT()->getP(), e->getH()->getP()) &&
	(*h)->head()->getP()->onLine(e->getT()->getP(), e->getH()->getP()))
      return;
  Points *pts;
  map<Edge *, Points *>::const_iterator i = epsmap.find(e);
  if (i == epsmap.end()) {
    pts = new Points;
    epsmap.insert(pair<Edge *, Points *>(e, pts));
  }
  else
    pts = i->second;
  pts->insert(pts->begin(), (e->getT()->getP()));
  pts->push_back(e->getH()->getP());
  FFE ffe(f, 0, pts);
  update(f, ffe, femap);
  for (int i = 0; i < e->HEdgesN(); ++i) {
    Face *g = e->getHEdge(i)->getF();
    pair<Face *, Face *> fg(f < g ? f : g, f < g ? g : f);
    ffemap.insert(pair<pair<Face *, Face *>, FFE>(fg, FFE()));
  }
}

void update (Face *f, const FFE &ffe, map<Face *, vector<FFE>> &femap)
{
  map<Face *, vector<FFE>>::iterator i = femap.find(f);
  if (i == femap.end()) {
    vector<FFE> ffes;
    ffes.push_back(ffe);
    femap.insert(pair<Face *, vector<FFE>>(f, ffes));
  }
  else
    i->second.push_back(ffe);
}

void intersectFFF (const map<pair<Face *, Face *>, FFE> &ffemap,
		   const map<Face *, vector<FFE>> &femap,
		   vector<pair<Points *, Points *>> &col)
{
  vector<pair<Points *, PTR<Point>>> psps;
  intersectFFFAux(ffemap, femap, col, psps);
  for (vector<pair<Points *, PTR<Point>>>::iterator i = psps.begin();
       i != psps.end(); ++i)
    i->first->push_back(i->second);
}

void intersectFFFAux (const map<pair<Face *, Face *>, FFE> &ffemap,
		      const map<Face *, vector<FFE>> &femap,
		      vector<pair<Points *, Points *>> &col,
		      vector<pair<Points *, PTR<Point>>> &psps)
{
  const unsigned int n = nthreads(8);
  unsigned int k = femap.size(), m = k/n, is = 0;
  FFFData fd[n];
  for (int i = 0; i < n; ++i) {
    fd[i].i = i;
    fd[i].is = is;
    is = i + 1 == n ? k : is + m;
    fd[i].ie = is;
    fd[i].ffemap = &ffemap;
    fd[i].femap = &femap;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) intersectFFFT, (void *) (fd + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  for (int i = 0; i < n; ++i) {
    psps.insert(psps.end(), fd[i].psps.begin(), fd[i].psps.end());
    col.insert(col.end(), fd[i].col.begin(), fd[i].col.end());
  }
}

void intersectFFFT (void *ptr)
{
  FFFData *fd = (FFFData *) ptr;
  BaseObject::addThread(fd->i);
  map<Face *, vector<FFE>>::const_iterator x = fd->femap->begin();
  for (int i = 0; i < fd->is; ++i, ++x)
    ;
  for (int i = fd->is; i < fd->ie; ++i, ++x)
    intersectFFF(x->first, x->second, *fd->ffemap, fd->col, fd->psps);
}

void intersectFFF (Face *f, const vector<FFE> &ed,
		   const map<pair<Face *, Face *>, FFE> &ffemap,
		   vector<pair<Points *, Points *>> &col,
		   vector<pair<Points *, PTR<Point>>> &psps)
{
  int n = ed.size(), pc = f->getP()->getPC();
  for (int i = 0; i + 1 < n; ++i) {
    const FFE &ei = ed[i];
    Face *fi = ei.f == f ? ei.g : ei.f;
    for (int j = i + 1; j < n; ++j) {
      const FFE &ej = ed[j];
      Face *fj = ej.f == f ? ej.g : ej.f;
      if ((fi || fj) && bboxOverlap(ei.bbox, ej.bbox)) {
	bool flag = !fi || !fj;
	map<pair<Face *, Face *>, FFE>::const_iterator y = ffemap.end();
	if (!flag) {
	  pair<Face *, Face *> fifj(fi < fj ? fi : fj, fi < fj ? fj : fi);
	  y = ffemap.find(fifj);
	  if (y != ffemap.end())
	    flag = !y->second.pts || f < fi && f < fj;
	}
	if (flag) {
	  Points pi, pj;
	  if (intersectEE(ei.pts->at(0), ei.pts->at(1), ej.pts->at(0), ej.pts->at(1),
			  pc, pi, pj))
	    col.push_back(pair<Points *, Points *>(ei.pts, ej.pts));
	  if (fi)
	    for (Points::iterator p = pi.begin(); p != pi.end(); ++p)
	      psps.push_back(pair<Points *, PTR<Point>>(ei.pts, *p));
	  if (fj)
	    for (Points::iterator p = pj.begin(); p != pj.end(); ++p)
	      psps.push_back(pair<Points *, PTR<Point>>(ej.pts, *p));
	  if (fi && fj && pi.size() == 1 && pj.size() == 1 &&
	      y != ffemap.end() && y->second.pts) {
	    psps.push_back(pair<Points *, PTR<Point>>(y->second.pts, pi[0]));
	    EEPoint *pee = dynamic_cast<EEPoint *>((Point *) pi[0]);
	    if (pee) {
	      pee->addFace(f);
	      pee->addFace(fi);
	      pee->addFace(fj);
	    }
	  }
	}
      }
    }
  }
}

void subedgesE (map<Edge *, Points *> &epsmap, bool oneway, set<Point *> &bpts,
		vector<pair<PTR<Point>, PTR<Point>>> &equiv)
{
  const unsigned int n = nthreads(8);
  unsigned int k = epsmap.size(), m = k/n, is = 0;
  SEEData see[n];
  for (int i = 0; i < n; ++i) {
    see[i].i = i;
    see[i].is = is;
    is = i + 1 == n ? k : is + m;
    see[i].ie = is;
    see[i].epsmap = &epsmap;
    see[i].oneway = oneway;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) subedgesET, (void *) (see + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  for (int i = 0; i < n; ++i) {
    bpts.insert(see[i].bpts.begin(), see[i].bpts.end());
    equiv.insert(equiv.end(), see[i].equiv.begin(), see[i].equiv.end());
  }
}

void subedgesET (void *ptr)
{
  SEEData *see = (SEEData *) ptr;
  BaseObject::addThread(see->i);
  map<Edge *, Points *>::const_iterator x = see->epsmap->begin();
  for (int i = 0; i < see->is; ++i, ++x)
    ;
  for (int i = see->is; i < see->ie; ++i, ++x)
    subedgesE(x->first, *x->second, see->oneway, see->bpts, see->equiv);
}

void subedgesE (Edge *e, Points &pts, bool oneway, set<Point *> &bpts,
		vector<pair<PTR<Point>, PTR<Point>>> &equiv)
{
  Point *t = e->getT()->getP(), *h = e->getH()->getP();
  if (pts.size() > 1)
    sort(pts.begin(), pts.end(), PointOrderPPP(t, h));
  pts.insert(pts.begin(), t);
  pts.push_back(h);
  subedgesE(pts, oneway, bpts, equiv);
}

void subedgesE (Points &pts, bool oneway, set<Point *> &bpts,
		vector<pair<PTR<Point>, PTR<Point>>> &equiv)
{
  if (oneway)
    subedgesE1(pts, bpts, equiv);
  else
    subedgesP(pts, equiv);
}

void subedgesE1 (Points &pts, set<Point *> &bpts,
		 vector<pair<PTR<Point>, PTR<Point>>> &equiv)
{
  Points pp;
  PTR<Point> u = new DiffPoint(pts.back(), pts[0]);
  int i = 0, n = pts.size();
  bool lfree = true, pblocked = false;
  while (i + 1 < n) {
    int j = i + 1;
    while (j + 1 < n && pts[j]->identical(pts[j+1]))
      ++j;
    for (int k = i + 1; k < j; ++k)
      equiv.push_back(pair<PTR<Point>, PTR<Point>>(pts[j], pts[k]));
    int s = 0;
    if (i + 1 == j) {
      EPPoint *epp = dynamic_cast<EPPoint *>((Point *) pts[j]);
      if (epp)
	s = PlaneRayAlignment(epp->getP(), u);
    }
    if (lfree && s <= 0) {
      pp.push_back(pts[i]);
      pp.push_back(pts[j]);
      pblocked = false;
    }
    else {
      if (pblocked)
	bpts.insert(pts[i]);
      pblocked = true;
    }
    i = j;
    lfree = s >= 0;
  }
  pts = pp;
}

void subedgesP (Points &pts, vector<pair<PTR<Point>, PTR<Point>>> &equiv)
{
  Points pp;
  int i = 0, n = pts.size();
  while (i + 1 < n) {
    int j = i + 1;
    while (j + 1 < n && pts[j]->identical(pts[j+1]))
      ++j;
    for (int k = i + 1; k < j; ++k)
      equiv.push_back(pair<PTR<Point>, PTR<Point>>(pts[j], pts[k]));
    pp.push_back(pts[i]);
    pp.push_back(pts[j]);
    i = j;
  }
  pts = pp;
}

void subedgesFF (map<pair<Face *, Face *>, FFE> &ffemap, bool oneway,
		 const set<Point *> &bpts,
		 vector<pair<PTR<Point>, PTR<Point>>> &equiv)
{
  const unsigned int n = nthreads(8);
  unsigned int k = ffemap.size(), m = k/n, is = 0;
  SEFData sef[n];
  for (int i = 0; i < n; ++i) {
    sef[i].i = i;
    sef[i].is = is;
    is = i + 1 == n ? k : is + m;
    sef[i].ie = is;
    sef[i].ffemap = &ffemap;
    sef[i].oneway = oneway;
    sef[i].bpts = &bpts;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) subedgesFFT, (void *) (sef + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  for (int i = 0; i < n; ++i)
    equiv.insert(equiv.end(), sef[i].equiv.begin(), sef[i].equiv.end());  
}

void subedgesFFT (void *ptr)
{
  SEFData *sef = (SEFData *) ptr;
  BaseObject::addThread(sef->i);
  map<pair<Face *, Face *>, FFE>::iterator x = sef->ffemap->begin();
  for (int i = 0; i < sef->is; ++i, ++x)
    ;
  for (int i = sef->is; i < sef->ie; ++i, ++x)
    subedgesFF(x->second, sef->oneway, *sef->bpts, sef->equiv);
}

void subedgesFF (FFE &ffe, bool oneway, const set<Point *> &bpts,
		 vector<pair<PTR<Point>, PTR<Point>>> &equiv)
{
  if (ffe.pts) {
    if (ffe.pts->size() > 1)
      sort(ffe.pts->begin(), ffe.pts->end(), 
	   PointOrderPPP(ffe.pts->at(0), ffe.pts->at(1)));
    if (oneway)
      subedgesFF(ffe, bpts, equiv);
    else
      subedgesP(*ffe.pts, equiv);
  }
}

void subedgesFF (FFE &ffe, const set<Point *> &bpts,
		 vector<pair<PTR<Point>, PTR<Point>>> &equiv)
{
  Points &pts = *ffe.pts, pp;
  PTR<Point> u = new DiffPoint(pts[1], pts[0]);
  int i = 0, n = pts.size();
  bool lfree = true;
  while (i + 1 < n) {
    int j = i + 1;
    while (j + 1 < n && pts[j]->identical(pts[j+1]))
      ++j;
    for (int k = i + 1; k < j; ++k)
      equiv.push_back(pair<PTR<Point>, PTR<Point>>(pts[j], pts[k]));
    int s = 0;
    if (i + 1 == j && i + 2 < n) {
      EEPoint *eh = dynamic_cast<EEPoint *>((Point *) pts[i+1]);
      if (eh) {
	Face *g = otherFace(ffe.f, ffe.g, eh->getFaces());
	if (g)
	  s = PlaneRayAlignment(g->getP(), u);
      }
    }
    if ((i > 0 || bpts.find(pts[i]) == bpts.end()) &&
	(i + 2 < n || bpts.find(pts[i+1]) == bpts.end()) &&
	lfree && s <= 0) {
      pp.push_back(pts[i]);
      pp.push_back(pts[j]);
    }
    lfree = s >= 0;
    i = j;
  }
  pts = pp;
}

Face * otherFace (Face *f, Face *g, const Faces &fa)
{
  for (Faces::const_iterator h = fa.begin(); h != fa.end(); ++h)
    if (*h != f && *h != g)
      return *h;
  return 0;
}

map<PTR<Point>, PTR<Point>> pMap (const map<Edge *, Points *> &epsmap,
				  const vector<pair<Edge *, Edge *>> &ee,
				  const vector<pair<Points *, Points *>> &col)
{
  map<PTR<Point>, PTR<Point>> pmap;
  for (vector<pair<Edge *, Edge *>>::const_iterator i = ee.begin();
       i != ee.end(); ++i)
    pMapEE(i->first, i->second, epsmap, pmap);
  for (vector<pair<Points *, Points *>>::const_iterator i = col.begin();
       i != col.end(); ++i)
    pMapCol(*i->first, *i->second, pmap);
  return pmap;
}

void pMapEE (Edge *e, Edge *f, const map<Edge *, Points *> &epsmap,
	     map<PTR<Point>, PTR<Point>> &pmap)
{
  map<Edge *, Points *>::const_iterator ei = epsmap.find(e), fi = epsmap.find(f);
  if (ei != epsmap.end() && fi != epsmap.end())
    pMapCol(*ei->second, *fi->second, pmap);
}

void pMapCol (const Points &pts1, const Points &pts2, 
	      map<PTR<Point>, PTR<Point>> &pmap)
{
  Points pts;
  for (Points::const_iterator p = pts1.begin(); p != pts1.end(); ++p)
    pts.push_back(getPoint(*p, pmap));
  for (Points::const_iterator p = pts2.begin(); p != pts2.end(); ++p)
    pts.push_back(getPoint(*p, pmap));
  sort(pts.begin(), pts.end(), PointOrderRP());
  int i = 0, n = pts.size();
  while (i + 1 < n) {
    int j = i + 1;
    while (j < n && pts[i]->identical(pts[j])) {
      pMapEquiv(pts[i], pts[j], pmap);
      ++j;
    }
    i = j;
  }		  
}

void pMapEquiv (Point *a, Point *b, map<PTR<Point>, PTR<Point>> &pmap)
{
  if (a != b) {
    Point *c = getPoint(a, pmap), *d = getPoint(b, pmap);
    if (c != d)
      pmap.insert(pair<PTR<Point>, PTR<Point>>(c, d));
  }
}

Point * getPoint (Point *p, const map<PTR<Point>, PTR<Point>> &pmap)
{
  while (true) {
    map<PTR<Point>, PTR<Point>>::const_iterator i = pmap.find(p);
    if (i == pmap.end())
      break;
    p = i->second;
  }
  return p;
}

void update (const vector<pair<PTR<Point>, PTR<Point>>> &equiv,
	      map<PTR<Point>, PTR<Point>> &pmap)
{
  for (vector<pair<PTR<Point>, PTR<Point>>>::const_iterator i = equiv.begin();
       i != equiv.end(); ++i)
    pMapEquiv(i->first, i->second, pmap);
}

Points SEdge::loop () const
{
  Points pts;
  const SEdge *e = this;
  do {
    pts.push_back(e->tail);
    e = e->next();
  }
  while (e != this);
  return pts;
}

Polyhedron * subfaces (Polyhedron *a, bool oneway,
		       const map<Edge *, Points *> &epsmap,
		       const map<Face *, vector<FFE>> &femap,
		       const map<PTR<Point>, PTR<Point>> &pmap)
{
  const unsigned int n = nthreads(16);
  unsigned int k = a->faces.size(), m = k/n, is = 0;
  SUData sud[n];
  for (int i = 0; i < n; ++i) {
    sud[i].i = i;
    sud[i].is = is;
    is = i + 1 == n ? k : is + m;
    sud[i].ie = is;
    sud[i].a = a;
    sud[i].oneway = oneway;
    sud[i].epsmap = &epsmap;
    sud[i].femap = &femap;
    sud[i].pmap = &pmap;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) subfacesT, (void *) (sud + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  Polyhedron *b = new Polyhedron;
  PVMap pvmap;
  for (int i = 0; i < n; ++i)
    for (SFaces::iterator s = sud[i].sf.begin(); s != sud[i].sf.end(); ++s)
      addFace(*s, b, pvmap);
  return b;
}

void subfacesT (void *ptr)
{
  SUData *sud = (SUData *) ptr;
  BaseObject::addThread(sud->i);
  for (int i = sud->is; i < sud->ie; ++i)
    subfaces(sud->a->faces[i], sud->oneway, *sud->epsmap,
	     *sud->femap, *sud->pmap, sud->sf);
}

void subfaces (Face *f, bool oneway, const map<Edge *, Points *> &epsmap,
	       const map<Face *, vector<FFE>> &femap,
	       const map<PTR<Point>, PTR<Point>> &pmap, SFaces &sf)
{
  SEdges ed = subedges(f, oneway, epsmap, femap, pmap);
  subfaces(f, oneway, ed, sf);
}

SEdges subedges (Face *f, bool oneway,
		 const map<Edge *, Points *> &epsmap,
		 const map<Face *, vector<FFE>> &femap,
		 const map<PTR<Point>, PTR<Point>> &pmap)
{
  vector<pair<PTR<Point>, PTR<Point>>> pp;
  HEdges ef = f->getBoundary()->edgeLoop();
  for (HEdges::iterator h = ef.begin(); h != ef.end(); ++h)
    subedgesH(*h, epsmap, pp);
  map<Face *, vector<FFE>>::const_iterator i = femap.find(f);
  if (i != femap.end())
    for (vector<FFE>::const_iterator j = i->second.begin();
	 j != i->second.end(); ++j)
      subedgesFFE(f, *j, oneway, pp);
  return subedgesPP(pp, f->getP()->getPC(), pmap);
}

void subedgesH (HEdge *h, const map<Edge *, Points *> &epsmap,
		vector<pair<PTR<Point>, PTR<Point>>> &pp)
{
  map<Edge *, Points *>::const_iterator i = epsmap.find(h->getE());
  if (i == epsmap.end()) {
    Points pts;
    pts.push_back(h->tail()->getP());
    pts.push_back(h->head()->getP());
    subedges(pts, true, false, pp);
  }
  else
    subedges(*i->second, h->getForward(), !h->getForward(), pp);
}

void subedges (const Points &pts, bool fflag, bool bflag,  
	       vector<pair<PTR<Point>, PTR<Point>>> &pp)
{
  for (int i = 0; i + 1 < pts.size(); i += 2) {
    if (fflag)
      pp.push_back(pair<PTR<Point>, PTR<Point>>(pts[i], pts[i+1]));
    if (bflag)
      pp.push_back(pair<PTR<Point>, PTR<Point>>(pts[i+1], pts[i]));
  }
}

void subedgesFFE (Face *f, const FFE &ffe, bool oneway,
		  vector<pair<PTR<Point>, PTR<Point>>> &pp)
{
  Points &pts = *ffe.pts;
  if (ffe.g) {
    bool fflag = !oneway || ffe.f == f, gflag = !oneway || ffe.f != f;
    subedges(pts, fflag, gflag, pp);
  }
  else
    for (int i = 0; i + 1 < pts.size(); i += 2) {
      CentroidPoint c(pts[i], pts[i+1]);
      if (f->contains(&c, true)) {
	pp.push_back(pair<PTR<Point>, PTR<Point>>(pts[i], pts[i+1]));
	pp.push_back(pair<PTR<Point>, PTR<Point>>(pts[i+1], pts[i]));
      }
    }
}

SEdges subedgesPP (const vector<pair<PTR<Point>, PTR<Point>>> &pp, int c,
		   const map<PTR<Point>, PTR<Point>> &pmap)
{
  bool flag = pmap.empty();
  set<pair<PTR<Point>, PTR<Point>>> pps;
  for (vector<pair<PTR<Point>, PTR<Point>>>::const_iterator i = pp.begin();
       i != pp.end(); ++i) {
    Point *t = flag ? (Point *) i->first : getPoint(i->first, pmap),
      *h = flag ? (Point *) i->second : getPoint(i->second, pmap);
    if (t != h)
      pps.insert(pair<PTR<Point>, PTR<Point>>(t, h));
  }
  SEdges ed;
  for (set<pair<PTR<Point>, PTR<Point>>>::const_iterator i = pps.begin();
       i != pps.end(); ++i)
    sedges(i->first, i->second, ed);
  return ed;
}

void sedges (Point *t, Point *h, SEdges &ed)
{
  SEdge *f = new SEdge(t, true), *b = new SEdge(h, false, f);
  f->twin = b;
  ed.push_back(f);
  ed.push_back(b);
}

void subfaces (Face *f, bool oneway, SEdges &ed, SFaces &sf)
{
  int c = f->getP()->getPC();
  linkSEdges(ed, c);
  SEdges outer, inner;
  for (SEdges::const_iterator e = ed.begin(); e != ed.end(); ++e) {
    SEdge *l = findLoop(*e);
    if (l)
      if (LeftTurn(l->tail, l->head(), l->next()->head(), c) == 1)
	outer.push_back(l);
      else
	inner.push_back(l);
  }
  sort(outer.begin(), outer.end(), SEdgeHeadOrder());
  sort(inner.begin(), inner.end(), SEdgeHeadOrder());
  SFaces nf;
  for (SEdges::iterator e = outer.begin(); e != outer.end(); ++e)
    nf.push_back(SFace(f, (*e)->loop()));
  for (SEdges::iterator e = inner.begin(); e != inner.end(); ++e)
    addInner(*e, oneway, c, nf);
  sf.insert(sf.end(), nf.begin(), nf.end());
  for (SEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    delete *e;
}

void linkSEdges (SEdges &ed, int c)
{
  sort(ed.begin(), ed.end(), SEdgeCWOrder(c));
  int i = 0, n = ed.size();
  while (i < n) {
    int j = i + 1;
    while (j < n && ed[i]->tail == ed[j]->tail) {
      ed[j-1]->cw = ed[j];
      ++j;
    }
    ed[j-1]->cw = ed[i];
    i = j;
  }
}

SEdge * findLoop (SEdge *e)
{
  if (e->flag || !e->forward)
    return 0;
  SEdge *l = e, *f = e;
  do {
    f->flag = true;
    if (SEdgeHeadOrder()(f, l))
      l = f;
    f = f->next();
    if (!(f && f->forward))
      return 0;
  }
  while (f != e);
  return l;
}

void addInner (SEdge *e, bool oneway, int c, SFaces &sf)
{
  for (int i = sf.size() - 1; i >= 0; --i)
    if (contains(e, c, sf[i].b[0])) {
      if (oneway)
	for (int j = 1; j < sf[i].b.size(); ++j)
	  if (contains(e, c, sf[i].b[j]))
	    return;
      sf[i].b.push_back(e->loop());
      return;
    }
}

bool contains (SEdge *e, int c, const Points &pts)
{
  int n = pts.size();
  Point *a = 0;
  SEdge *e0 = e;
  while (true) {
    a = e->head();
    bool flag = true;
    for (int i = 0; flag && i < n; ++i)
      flag = a != pts[i];
    if (flag)
      break;
    e = e->next();
    if (e == e0)
      return false;
  }
  Point *et = pts.back();
  bool res = false;
  for (int i = 0; i < n; ++i) {
    Point *eh = pts[i];
    if (rayEdgeIntersection(a, Rdir, et, eh, c))
      res = !res;
    et = eh;
  }
  return res;
}

bool rayEdgeIntersection (Point *a, Point *r, Point *t, Point *h, int c)
{
  PTR<Point> at = new DiffPoint(t, a), ah = new DiffPoint(h, a);
  if (Cross2(at, r, c) == Cross2(ah, r, c))
    return false;
  PTR<Point> th = new DiffPoint(h, t);
  return Cross2(at, th, c) == Cross2(r, th, c);
}

void MFace::update (HEdge *e)
{
  if (!h->getF())
    h = e;
  else
    for (int i = 0; i < inner.size(); ++i)
      if (!inner[i]->getF()) {
	inner[i] = e;
	break;
      }
}

HEdges MFace::boundaryHEdges () const
{
  HEdges ed = h->edgeLoop();
  for (HEdges::const_iterator e = inner.begin(); e != inner.end(); ++e) {
    HEdges eed = (*e)->edgeLoop();
    ed.insert(ed.end(), eed.begin(), eed.end());
  }
  return ed;
}

void MFace::triangulate (Triangles &tr)
{
  vector<Points> pp;
  pp.push_back(h->pointLoop());
  for (HEdges::iterator e = inner.begin(); e != inner.end(); ++e)
    pp.push_back((*e)->pointLoop());
  if (pp.size() == 1 && pp[0].size() == 3)
    tr.push_back(Triangle(pp[0][0], pp[0][1], pp[0][2]));
  else {
    vector<Points *> reg;
    for (vector<Points>::iterator i = pp.begin(); i != pp.end(); ++i)
      reg.push_back(&*i);
    ::triangulate(reg, p->getPC(), tr);
  }
}

void addFace (const SFace &sf, Polyhedron *a, PVMap &pvmap)
{
  Vertices ve = loop(sf.b[0], a, pvmap);
  if (currentFace(ve, a))
    return;
  MFace *f = new MFace(addLoop(ve, a), sf.f->getP());
  for (int i = 1; i < sf.b.size(); ++i)
    f->addInner(addLoop(loop(sf.b[i], a, pvmap), a));
  a->faces.push_back(f);
}

Vertices loop (const Points &pts, Polyhedron *a, PVMap &pvmap)
{
  Vertices ve;
  for (Points::const_iterator p = pts.begin(); p != pts.end(); ++p)
    ve.push_back(a->getVertex(*p, pvmap));
  return ve;
}

bool currentFace (const Vertices &ve, Polyhedron *a)
{
  int n = ve.size();
  Edge *e = a->getEdge(ve[0], ve[1]);
  for (int i = 0; i < e->HEdgesN(); ++i) {
    HEdge *h = e->getHEdge(i);
    if (outerLoop(h)) {
      Vertices vh = h->loop();
      if (vh.size() == n) {
	bool flag = true;
	if (h->tail() == ve[0])
	  for (int j = 0; flag && j < n; ++j)
	    flag = ve[j] == vh[j];
	else
	  for (int j = 0; flag && j < n; ++j)
	    flag = ve[j] == vh[(n+1-j)%n];
	if (flag)
	  return true;
      }
    }
  }
  return false;
}

bool outerLoop (HEdge *h)
{
  HEdges ed = h->getF()->getBoundary()->edgeLoop();
  return find(ed.begin(), ed.end(), h) != ed.end();
}

HEdge * addLoop (const Vertices &ve, Polyhedron *a)
{
  Vertex *t = ve.back();
  HEdges ed;
  for (Vertices::const_iterator h = ve.begin(); h != ve.end(); ++h) {
    HEdge *f = a->getHEdge(t, *h);
    ed.push_back(f);
    t = *h;
  }
  for (int i = 0; i + 1 < ed.size(); ++i)
    ed[i]->setNext(ed[i+1]);
  ed.back()->setNext(ed[0]);
  return ed[0];
}

void deleteMaps (map<Edge *, Points *> &epsmap,
		 map<pair<Face *, Face *>, FFE> &ffemap)
{
  for (map<Edge *, Points *>::iterator i = epsmap.begin(); i != epsmap.end(); ++i)
    delete i->second;
  for (map<pair<Face *, Face *>, FFE>::iterator i = ffemap.begin();
       i != ffemap.end(); ++i)
    delete i->second.pts;
}

Polyhedron * triangulate (Polyhedron *a)
{
  Triangles tr = triangulate(a->faces);
  Polyhedron *b = new Polyhedron;
  PVMap pvmap;
  int k = 0;
  for (Triangles::iterator t = tr.begin(); t != tr.end(); ++t, ++k)
    b->addTriangle(t->a, t->b, t->c, pvmap);
  return b;
}

Triangles triangulate (const Faces &fa)
{
  const unsigned int n = nthreads(8);
  unsigned int k = fa.size(), m = k/n, is = 0;
  TRData trd[n];
  for (int i = 0; i < n; ++i) {
    trd[i].i = i;
    trd[i].is = is;
    is = i + 1 == n ? k : is + m;
    trd[i].ie = is;
    trd[i].fa = &fa;
  }
  pthread_t threads[n];
  for (int i = 0; i < n; ++i)
    pthread_create(threads + i, 0,
		   (void* (*)(void*)) triangulateT, (void *) (trd + i));
  for (int i = 0; i < n; ++i)
    pthread_join(threads[i], 0);
  Triangles tr;
  for (int i = 0; i < n; ++i)
    tr.insert(tr.end(), trd[i].tr.begin(), trd[i].tr.end());
  return tr;
}

void triangulateT (void *ptr)
{
  TRData *trd = (TRData *) ptr;
  BaseObject::addThread(trd->i);
  for (int i = trd->is; i < trd->ie; ++i)
    trd->fa->at(i)->triangulate(trd->tr);
}

PVMap pvMap (const Points &pts, Polyhedron *a)
{
  vector<PointR> pr;
  for (Points::const_iterator p = pts.begin(); p != pts.end(); ++p)
    pr.push_back(PointR(*p));
  sort(pr.begin(), pr.end());
  int i = 0, n = pr.size();
  PVMap pvmap;
  while (i < n) {
    Vertex *v = a->getVertex(pr[i].a);
    pvmap.insert(pair<Point *, Vertex *>(pr[i].a, v));
    int j = i + 1;
    while (j < n && pr[i] == pr[j]) {
      pvmap.insert(pair<Point *, Vertex *>(pr[j].a, v));
      ++j;
    }
    i = j;
  }
  return pvmap;
}

bool inSet (bool ina, bool inb, SetOp op)
{
  switch (op) {
  case Union:
    return ina || inb;
  case Intersection:
    return ina && inb;
  case Complement:
    return ina && !inb;
  }
}

Polyhedron * overlay (Polyhedron **poly, int n)
{
  Points pts;
  for (int i = 0; i < n; ++i)
    for (Vertices::iterator v = poly[i]->vertices.begin();
	   v != poly[i]->vertices.end(); ++v)
      pts.push_back((*v)->getP());
  Polyhedron *c = new Polyhedron;
  PVMap pvmap = pvMap(pts, c);
  for (int i = 0; i < n; ++i) {
    const Faces &fa = poly[i]->faces;
    for (Faces::const_iterator f = fa.begin(); f != fa.end(); ++f)
      c->addTriangle(*f, pvmap, (*f)->getP());
  }
  Polyhedron *d = subdivide(c, false), *e = triangulate(d);
  delete c;
  delete d;
  return e;
}

Polyhedron * multiUnion (Polyhedron **poly, int n)
{
  Polyhedron *c = overlay(poly, n), *d = new Polyhedron;
  for (int i = 0; i < n; ++i)
    poly[i]->computeWindingNumbers();
  c->formCells();
  set<Cell *> cin;
  for (int i = 1; i < c->cells.size(); ++i) {
    Cell *ci = c->cells[i];
    PTR<Point> p = ci->interiorPoint();
    bool flag = false;
    for (int j = 0; !flag && j < n; ++j)
      flag = poly[j]->contains(p);
    if (flag)
      cin.insert(ci);
  }
  PVMap pvmap;
  for (Faces::iterator f = c->faces.begin(); f != c->faces.end(); ++f) {
    bool in1 = cin.find((*f)->getHFace(0)->getS()->getC()) != cin.end(),
      in2 = cin.find((*f)->getHFace(1)->getS()->getC()) != cin.end();
    if (in1 != in2) {
      Points pf = (*f)->getBoundary()->pointLoop();
      if (in2)
	d->addTriangle(pf[0], pf[1], pf[2], pvmap, (*f)->getP());
      else
	d->addTriangle(pf[2], pf[1], pf[0], pvmap);
    }
  }
  delete c;
  return d;
}

Polyhedron * coalesce (Polyhedron *a)
{
  Polyhedron *b = coalesceFaces(a);
  coalesceEdges(b);
  Polyhedron *c = triangulate(b);
  delete b;
  return c;
}

Polyhedron * coalesceFaces (Polyhedron *a)
{
  vector<set<Face *>> fg = groupFaces(a);
  SFaces sf;
  for (vector<set<Face *>>::iterator f = fg.begin(); f != fg.end(); ++f)
    coalesceFace(*f, sf);
  Polyhedron *b = new Polyhedron;
  PVMap pvmap;
  for (SFaces::iterator f = sf.begin(); f != sf.end(); ++f)
    addFace(*f, b, pvmap);
  return b;
}

vector<set<Face *>> groupFaces (Polyhedron *a)
{
  vector<set<Face *>> res;
  set<Face *> done;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f)
    if (done.insert(*f).second) {
      Faces st;
      set<Face *> fg;
      st.push_back(*f);
      while (!st.empty()) {
	Face *g = st.back();
	st.pop_back();
	if (g->getP()->coplanar((*f)->getP()) && fg.insert(g).second) {
	  done.insert(g);
	  HEdges ed = g->getBoundary()->edgeLoop();
	  for (HEdges::iterator h = ed.begin(); h != ed.end(); ++h)
	    if ((*h)->getE()->HEdgesN() == 2 &&
		(*h)->getForward() != (*h)->ccw()->getForward())
	      st.push_back((*h)->ccw()->getF());
	}
      }
      res.push_back(fg);
    }
  return res;
}

void coalesceFace (const set<Face *> &fs, SFaces &sf)
{
  if (fs.size() == 1) {
    Face *f = *fs.begin();
    sf.push_back(SFace(f, f->getBoundary()->pointLoop()));
    return;
  }
  set<pair<PTR<Point>, PTR<Point>>> pp;
  for (set<Face *>::const_iterator f = fs.begin(); f != fs.end(); ++f) {
    HEdges ef = (*f)->getBoundary()->edgeLoop();
    for (HEdges::iterator h = ef.begin(); h != ef.end(); ++h) {
      pair<PTR<Point>, PTR<Point>> pf((*h)->tail()->getP(), (*h)->head()->getP()),
	pb(pf.second, pf.first);
      set<pair<PTR<Point>, PTR<Point>>>::iterator ib = pp.find(pb);
      if (ib == pp.end())
	pp.insert(pf);
      else
	pp.erase(ib);
    }
  }
  SEdges ed;
  for (set<pair<PTR<Point>, PTR<Point>>>::iterator i = pp.begin(); i != pp.end(); ++i)
    sedges(i->first, i->second, ed);
  subfaces(*fs.begin(), true, ed, sf);
}

void coalesceEdges (Polyhedron *a)
{
  vector<pair<HEdge *, HEdge *>> hh;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    MFace *m = dynamic_cast<MFace *>(*f);
    coalesceLoop(m->getBoundary(), hh);
    const HEdges &he = m->getInner();
    for (HEdges::const_iterator h = he.begin(); h != he.end(); ++h)
      coalesceLoop(*h, hh);
  }
  for (vector<pair<HEdge *, HEdge *>>::iterator i = hh.begin(); i != hh.end(); ++i)
    coalesceEdges(a, i->first, i->second);
}

void coalesceLoop (HEdge *h, vector<pair<HEdge *, HEdge *>> &hh)
{
  while (h->tail()->getP()->onLine(h->head()->getP(),
				   h->getNext()->head()->getP()))
    h = h->getNext();
  h = h->getNext();
  HEdge *e = h;
  do {
    HEdge *f = e;
    while (coalesceEdge(e, f->getNext()))
      f = f->getNext();
    if (e != f)
      hh.push_back(pair<HEdge *, HEdge *>(e, f));
    e = f->getNext();
  }
  while (e != h);
}

bool coalesceEdge (HEdge *e, HEdge *f)
{
  int n = e->getE()->HEdgesN();
  if (f->getE()->HEdgesN() != n)
    return false;
  HEdge *ee = e->ccw(), *ff = f->ccw();
  for (int i = 0; i + 1 < n; ++i) {
    if (ee->getF() != ff->getF())
      return false;
    ee = ee->ccw();
    ff = ff->ccw();
  }
  return e->tail()->getP()->onLine(e->head()->getP(), f->head()->getP());
}

void coalesceEdges (Polyhedron *a, HEdge *s, HEdge *e)
{
  MFace *f = dynamic_cast<MFace *>(s->getF());
  HEdge *en = e->getNext(), *sp = en;
  while (sp->getNext() != s)
    sp = sp->getNext();
  HEdge *se = a->getHEdge(s->tail(), e->head());
  se->setF(f);
  sp->setNext(se);
  se->setNext(en);
  while (s != en) {
    s->setF(0);
    s = s->getNext();
  }
  f->update(se);
}

// shapes

Polyhedron * box (double *b, bool perturb)
{
  Polyhedron *p = new Polyhedron;
  p->getVertex(b[0], b[2], b[4], perturb);
  p->getVertex(b[1], b[2], b[4], perturb);
  p->getVertex(b[1], b[3], b[4], perturb);
  p->getVertex(b[0], b[3], b[4], perturb);
  p->getVertex(b[0], b[2], b[5], perturb);
  p->getVertex(b[1], b[2], b[5], perturb);
  p->getVertex(b[1], b[3], b[5], perturb);
  p->getVertex(b[0], b[3], b[5], perturb);
  p->addTriangle(p->vertices[0], p->vertices[2], p->vertices[1]);
  p->addTriangle(p->vertices[0], p->vertices[3], p->vertices[2]);
  p->addTriangle(p->vertices[0], p->vertices[1], p->vertices[5]);
  p->addTriangle(p->vertices[0], p->vertices[5], p->vertices[4]);
  p->addTriangle(p->vertices[1], p->vertices[2], p->vertices[6]);
  p->addTriangle(p->vertices[1], p->vertices[6], p->vertices[5]);
  p->addTriangle(p->vertices[2], p->vertices[3], p->vertices[7]);
  p->addTriangle(p->vertices[2], p->vertices[7], p->vertices[6]);
  p->addTriangle(p->vertices[3], p->vertices[0], p->vertices[4]);
  p->addTriangle(p->vertices[3], p->vertices[4], p->vertices[7]);
  p->addTriangle(p->vertices[4], p->vertices[5], p->vertices[6]);
  p->addTriangle(p->vertices[4], p->vertices[6], p->vertices[7]);
  return p;
}

Polyhedron * sphere (double ox, double oy, double oz, double r, double err)
{
  Polyhedron *a = sphere(err/r);
  Polyhedron *b = new Polyhedron;
  PVMap pvmap;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    Vertices va = (*f)->getBoundary()->loop(), vb;
    for (Vertices::iterator v = va.begin(); v != va.end(); ++v) {
      PVMap::iterator i = pvmap.find((*v)->getP());
      if (i == pvmap.end()) {
	PV3<double> vi = (*v)->getP()->getApproxMid(1.0);
	Vertex *w = b->getVertex(ox + r*vi.x, oy + r*vi.y, oz + r*vi.z);
	pvmap.insert(PVPair((*v)->getP(), w));
	vb.push_back(w);
      }
      else
	vb.push_back(i->second);
    }
    b->addTriangle(vb[0], vb[1], vb[2]);
  }
  delete a;
  return b;
}

Polyhedron * sphere (double err)
{
  Polyhedron *a = octohedron();
  while (sphereError(a) > err) {
    Polyhedron *b = sphereRefine(a);
    delete a;
    a = b;
  }
  return a;
}

Polyhedron * octohedron ()
{
  Polyhedron *a = new Polyhedron;
  a->getVertex(1.0, 0.0, 0.0);
  a->getVertex(-1.0, 0.0, 0.0);
  a->getVertex(0.0, 1.0, 0.0);
  a->getVertex(0.0, -1.0, 0.0);
  a->getVertex(0.0, 0.0, 1.0);
  a->getVertex(0.0, 0.0, -1.0);
  a->addTriangle(a->vertices[2], a->vertices[4], a->vertices[0]);
  a->addTriangle(a->vertices[1], a->vertices[4], a->vertices[2]);
  a->addTriangle(a->vertices[3], a->vertices[4], a->vertices[1]);
  a->addTriangle(a->vertices[0], a->vertices[4], a->vertices[3]);
  a->addTriangle(a->vertices[5], a->vertices[2], a->vertices[0]);
  a->addTriangle(a->vertices[5], a->vertices[1], a->vertices[2]);
  a->addTriangle(a->vertices[5], a->vertices[3], a->vertices[1]);
  a->addTriangle(a->vertices[5], a->vertices[0], a->vertices[3]);
  return a;
}

double sphereError (Polyhedron *a)
{
  double d = 1.0;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    Vertices ve = (*f)->getBoundary()->loop();
    PV3<double> u = ve[0]->getP()->getApproxMid(1.0),
      v = ve[1]->getP()->getApproxMid(1.0), w = ve[2]->getP()->getApproxMid(1.0),
      p = (u + v + w)/3.0;
    d = min(d, p.dot(p));
  }
  return 1.0 - sqrt(d);
}

Polyhedron * sphereRefine (Polyhedron *a)
{
  Polyhedron *b = new Polyhedron;
  map<Edge *, Vertex *> evmap;
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e) {
    PV3<double> pm = (*e)->getT()->getP()->getApproxMid(1.0) +
      (*e)->getH()->getP()->getApproxMid(1.0);
    double k = 1.0/pm.length();
    Vertex *vm = b->getVertex(k*pm.x, k*pm.y, k*pm.z);
    evmap.insert(pair<Edge *, Vertex *>(*e, vm));
  }
  PVMap pvmap;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    HEdge *e = (*f)->getBoundary();
    Vertices ve = e->loop();
    Vertex *u = b->getVertex(ve[0]->getP(), pvmap), *v = b->getVertex(ve[1]->getP(), pvmap),
      *w = b->getVertex(ve[2]->getP(), pvmap), *uv = evmap.find(e->getE())->second,
      *vw = evmap.find(e->getNext()->getE())->second,
      *wu = evmap.find(e->getNext()->getNext()->getE())->second;
    b->addTriangle(u, uv, wu);
    b->addTriangle(uv, v, vw);
    b->addTriangle(vw, w, wu);
    b->addTriangle(uv, vw, wu);
  }
  return b;
}

// n tetrahedra (t, t + a, t + b, t + c) with the coordinates of t in [0, u]
// and the coordinates of a, b, and c in [0, v]
Polyhedron * randomTets (int n, double u, double v)
{
  Polyhedron *a = new Polyhedron;
  for (int i = 0; i < n; ++i)
    randomTet(a, u, v);
  Polyhedron *b = subdivide(a, false), *c = triangulate(b);
  delete b;
  delete a;
  return c;
}

void randomTet (Polyhedron *a, double u, double v)
{
  PTR<Point> p0 = randomPoint(u),
    p[] = {p0, new SumPoint(p0, randomPoint(v)),
	   new SumPoint(p0, randomPoint(v)), new SumPoint(p0, randomPoint(v))};
  if (Orientation(p[0], p[1], p[2], p[3]) == -1)
    reverse(p, p + 4);
  Vertex *v1 = a->getVertex(p[0]), *v2 = a->getVertex(p[1]),
    *v3 = a->getVertex(p[2]), *v4 = a->getVertex(p[3]);
  a->addTriangle(v1, v2, v4);
  a->addTriangle(v2, v3, v4);
  a->addTriangle(v3, v1, v4);
  a->addTriangle(v3, v2, v1);
}

Point * randomPoint (double d)
{
  return new Point(randomNumber(0.0, d), randomNumber(0.0, d),
		   randomNumber(0.0, d), false);
}

// debug

void find (Point *p, const Vertices &vertices)
{
  for (int i = 0; i < vertices.size(); ++i)
    if (p == vertices[i]->getP())
      cerr << i << " ";
  cerr << endl;
}

int find (Vertex *v, const Vertices &vertices)
{
  for (int i = 0; i < vertices.size(); ++i)
    if (v == vertices[i])
      return i;
  return -1;
}

int find (Edge *e, const Edges &edges)
{
  for (int i = 0; i < edges.size(); ++i)
    if (e == edges[i])
      return i;
  return -1;
}

int find (Face *f, const Faces &faces)
{
  for (int i = 0; i < faces.size(); ++i)
    if (f == faces[i])
      return i;
  return -1;
}

int find (HFace *f, const HFaces &faces)
{
  for (int i = 0; i < faces.size(); ++i)
    if (f == faces[i])
      return i;
  return -1;
}

void find (Plane *p, const Faces &faces)
{
  for (int i = 0; i < faces.size(); ++i)
    if (faces[i]->getP()->coplanar(p))
      cerr << i << " ";
  cerr << endl;
}

int find (Shell *s, const Shells &shells)
{
  for (int i = 0; i < shells.size(); ++i)
    if (s == shells[i])
      return i;
  return -1;
}

int find (Cell *c, const Cells &cells)
{
  for (int i = 0; i < cells.size(); ++i)
    if (c == cells[i])
      return i;
  return -1;
}

void pp1 (const PV3<double> &p)
{
  cerr << setprecision(16);
  cerr << "(" << p.x << " " << p.y << " " << p.z << ")";
}

void pp (const PV3<double> &p)
{
  pp1(p);
  cerr << endl;
}

void pp (Point *p)
{
  pp(p->getApproxMid());
}

void pps (const Points &pts)
{
  cerr << "(" << endl;
  for (Points::const_iterator p = pts.begin(); p != pts.end(); ++p)
    pp(*p);
  cerr << ")" << endl;
}

void pv (Vertex *v)
{
  pp(v->getP());
}

void pvs (const Vertices &ve)
{
  cerr << "(" << endl;
  for (Vertices::const_iterator v = ve.begin(); v != ve.end(); ++v)
    pv(*v);
  cerr << ")" << endl;
}

void ptr (const Triangle &t)
{
  cerr << "(";
  pp(t.a);
  pp(t.b);
  pp(t.c);
  cerr << ")";
}

void ptrs (const Triangles &tr)
{
  cerr << "(";
  for (Triangles::const_iterator t = tr.begin(); t != tr.end(); ++t)
    ptr(*t);
  cerr << ")" << endl;
}

void pe (Edge *e)
{
  cerr << "(";
  pp1(e->getT()->getP()->getApproxMid());
  pp1(e->getH()->getP()->getApproxMid());
  cerr << ")" << endl;
}

void pes (const Edges &ed)
{
  cerr << "(";
  for (Edges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e);
  cerr << ")" << endl;
}

void pe (HEdge *e)
{
  cerr << "(";
  pp1(e->tail()->getP()->getApproxMid());
  pp1(e->head()->getP()->getApproxMid());
  cerr << ")" << endl;
}

void pes (const HEdges &ed)
{
  cerr << "(";
  for (HEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e);
  cerr << ")" << endl;
}

void pl (HEdge *e)
{
  Vertices ve = e->loop();
  pvs(ve);
}

void pe (SEdge *e)
{
  cerr << "(";
  pp1(e->tail->getApproxMid());
  pp1(e->head()->getApproxMid());
  cerr << ")" << endl;
}

void pes (const SEdges &ed)
{
  cerr << "(";
  for (SEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e);
  cerr << ")" << endl;
}

void pl (SEdge *e)
{
  Points pts = e->loop();
  pps(pts);
}

void pf (Face *f)
{
  pvs(f->getBoundary()->loop());
}

void pfs (const Faces &fa)
{
  cerr << "(";
  for (Faces::const_iterator f = fa.begin(); f != fa.end(); ++f)
    pf(*f);
  cerr << ")" << endl;
}

void pfs (const HFaces &fa)
{
  cerr << "(";
  for (HFaces::const_iterator f = fa.begin(); f != fa.end(); ++f)
    pf((*f)->getF());
  cerr << ")" << endl;
}

TrianglePlane * tpl (Plane *p)
{
  return dynamic_cast<TrianglePlane *>(p);
}

SumPoint * spt (Point *v)
{
  return dynamic_cast<SumPoint *>(v);
}

EPPoint * eppt (Point *v)
{
  return dynamic_cast<EPPoint *>(v);
}		      

EEPoint * eept (Point *v)
{
  return dynamic_cast<EEPoint *>(v);
}

double distancePP (Point *v, Point *w)
{
  PV3<double> u = v->getApproxMid() - w->getApproxMid();
  return u.length();
}

double distance (Point *a, Point *t, Point *h)
{
  PV3<double> ap = a->getApproxMid(), tp = t->getApproxMid(),
    u = h->getApproxMid() - tp;
  double k = (ap - tp).dot(u)/u.dot(u);
  PV3<double> p = tp + k*u, w = ap - p;
  return w.length();
}

double distance (Point *v, Plane *p)
{
  PlaneData<double> pd = p->getApproxMid();
  double k = pd.n.length(), d = pd.n.dot(v->getApproxMid()) + pd.k;
  return d/k;
}

double distance (Point *v, Point *a, Point *b, Point *c)
{
  PV3<double> aa = a->getApproxMid(), bb = b->getApproxMid(),
    cc = c->getApproxMid(), n = (cc - bb).cross(aa - bb);
  double k = n.length(), d = n.dot(v->getApproxMid() - bb);
  return d/k;
}

double distanceEE (Point *a, Point *b, Point *c, Point *d)
{
  PV3<double> aa = a->getApproxMid(), bb = b->getApproxMid(),
    cc = c->getApproxMid(), dd = d->getApproxMid(), u = bb - aa, v = dd - cc,
    w = u.cross(v);
  double k = w.length();
  return w.dot(cc - aa)/k;
}

double distance (Edge *e, Edge *f)
{
  return distanceEE(e->getT()->getP(), e->getH()->getP(),
		    f->getT()->getP(), f->getH()->getP());
}

void edgesNM (Shell *s)
{
  set<Edge *> es;
  const HFaces &hf = s->getHFaces();
  for (HFaces::const_iterator f = hf.begin(); f != hf.end(); ++f) {
    HEdges he = (*f)->getF()->getBoundary()->edgeLoop();
    for (HEdges::iterator e = he.begin(); e != he.end(); ++e)
      es.insert((*e)->getE());
  }
  for (set<Edge *>::iterator e = es.begin(); e != es.end(); ++e) {
    int n = 0;
    for (int i = 0; i < (*e)->HEdgesN(); ++i) {
      Face *f = (*e)->getHEdge(i)->getF();
      for (int j = 0; j < 2; ++j)
	if (f->getHFace(j)->getS() == s)
	  ++n;
    }
    if (n != 2)
      cerr << "non-manifold edge " << n << endl;
  }
}

void edgesNM (Polyhedron *a, Edges &ed)
{
  for (int i = 0; i < a->edges.size(); ++i) {
    Edge *e = a->edges[i];
    int n = e->HEdgesN();
    if (n == 0) {
      cerr << "no hedges " << i << endl;
      ed.push_back(e);
    }
    else if (n == 1) {
      cerr << "dangling edge " << i << endl;
      ed.push_back(e);
    }
    else if (n > 2) {
      cerr << "non-manifold edge " << i << " hedges " << n << endl;
      ed.push_back(e);
    }
    else if (e->getHEdge(0)->getForward() == e->getHEdge(1)->getForward()) {
      cerr << "same sign hedges" << endl;
      ed.push_back(e);
    }
  }
}

void edgesNM (Polyhedron *a)
{
  Edges ed;
  edgesNM(a, ed);
}

bool vertexNM (Vertex *v)
{
  HEdges ed = v->outgoingHEdges();
  int k = 0;
  HEdge *e = ed[0];
  do {
    e = e->getNext()->getNext()->ccw();
    ++k;
  }
  while (e != ed[0]);
  return k < ed.size();
}

void verticesNM (Polyhedron *a)
{
  for (int i = 0; i < a->vertices.size(); ++i)
    if (vertexNM(a->vertices[i]))
      cerr << "non manifold vertex " << i << endl;
}

double * rotationMatrix (double t, double *u)
{
  double c = cos(t), d = 1.0 - c, s = sin(t), ux = u[0], uy = u[1], uz = u[2],
    *m = new double[9];
  m[0] = ux*ux*d + c;
  m[1] = ux*uy*d + uz*s;
  m[2] = ux*uz*d - uy*s;
  m[3] = ux*uy*d - uz*s;
  m[4] = uy*uy*d + c;
  m[5] = uy*uz*d + ux*s;
  m[6] = ux*uz*d + uy*s;
  m[7] = uy*uz*d - ux*s;
  m[8] = uz*uz*d + c;
  return m;
}

Point * rotate (Point *p, double *m)
{
  PV3<double> q = p->getApproxMid();
  double x = m[0]*q.x + m[3]*q.y + m[6]*q.z,
    y = m[1]*q.x + m[4]*q.y + m[7]*q.z,
    z = m[2]*q.z + m[5]*q.y + m[8]*q.z;
  return new Point(x, y, z, false);
}

Polyhedron * rotate (Polyhedron *a, double t,  double *u)
{
  Polyhedron *b = new Polyhedron;
  PVMap pvmap;
  double *m = rotationMatrix(t, u);
  for (Vertices::const_iterator v = a->vertices.begin(); v != a->vertices.end(); ++v)
    pvmap.insert(PVPair((*v)->getP(), b->getVertex(rotate((*v)->getP(), m))));
  for (Faces::const_iterator f = a->faces.begin(); f != a->faces.end(); ++f)
    b->addTriangle(*f, pvmap);
  delete [] m;
  return b;
}
