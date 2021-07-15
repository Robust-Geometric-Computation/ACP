#include "cspace.h"

VertexSet intersection (const VertexSet &a, const VertexSet &b)
{
  VertexSet c;
  set_intersection(a.begin(), a.end(), b.begin(), b.end(),
		   inserter(c, c.begin()));
  return c;
}

PTR<Angle> Angle::mpi = new Angle(-1.0, 0.0);
PTR<Angle> Angle::ppi = new Angle(-1.0, 0.0);

bool Angle::order (Angle *b)
{
  if (this == b || this == ppi || b == mpi)
    return false;
  if (this == mpi || b == ppi)
    return true;
  return AngleOrder(this, b) == 1;
}

bool inInterval (Angle *a, Angle *s, Angle *e)
{
  return s->order(e) ? s->order(a) && a->order(e)
    : e->order(a) && a->order(s);
}

bool intervalOverlap (Angle *s1, Angle *e1, Angle *s2, Angle *e2)
{
  return s2->order(e1) && s1->order(e2);
}

void intersect (const AngleInterval &x, const AngleInterval &y,
		AngleIntervals &res)
{
  Angle *ns = x.first->order(y.first) ? y.first : x.first,
    *ne = x.second->order(y.second) ? x.second : y.second;
  if (ns->order(ne))
    res.push_back(AngleInterval(ns, ne));
}

AngleIntervals intersect (const AngleIntervals &a, const AngleIntervals &b)
{
  AngleIntervals res;
  for (AngleIntervals::const_iterator ai = a.begin(); ai != a.end(); ++ai)
    for (AngleIntervals::const_iterator bi = b.begin(); bi != b.end(); ++bi)
      intersect(*ai, *bi, res);
  return res;
}

Angles sinCosAngles (Object<SinCosData> *x, Angle *s, Angle *e, const VertexSet &vs)
{
  Angles as;
  double l = AngleParameter(s).getApprox(1.0).x.l,
    u = AngleParameter(e).getApprox(1.0).x.u;
  if (AngleUpper(e) < 1)
    sinCosAngles(x, false, l, u, vs, as);
  else if (AngleUpper(s) > -1)
    sinCosAngles(x, true, l, u, vs, as);
  else {
    sinCosAngles(x, false, l, 1.0, vs, as);
    sinCosAngles(x, true, -1.0, u, vs, as);
  }
  return as;
}

void sinCosAngles (Object<SinCosData> *x, bool upper, double l, double u,
		   const VertexSet &vs, Angles &as)
{
  PTR<Object<Poly>> p = new SinCosPoly(x, upper);
  vector<PTR<Object<Scalar>>> rts = getRoots(p, l, u);
  for (int i = 0; i < rts.size(); ++i)
    as.push_back(new ParametricAngle(upper, rts[i], vs));
}

Angles anglesEF (Edge *e, Face *f, bool aflag)
{
  PTR<Object<SinCosData>> x = new PolyEF(e, f, aflag);
  VertexSet vs;
  vs.insert(e->getT());
  vs.insert(e->getH());
  Vertices ve = f->getBoundary()->loop();
  vs.insert(ve.begin(), ve.end());
  return sinCosAngles(x, Angle::mpi, Angle::ppi, vs);
}

Angles anglesSFL (const EFKey &k, bool aflag)
{
  PTR<Object<SinCosData>> x = new PolySFL(k, aflag);
  VertexSet vs(k.v, k.v + 5);
  return sinCosAngles(x, Angle::mpi, Angle::ppi, vs);
}

bool onEdge (CPoint *a, CPoint *t, CPoint *h)
{
  if ( a == t || a == h)
    return false;
  Angle *ta = t->getA(), *ha = h->getA();
  return ta == ha ? COrder(a, t, h) == 1 && COrder(a, h, t) == 1
    : inInterval(a->getA(), ta, ha);
}

void Spiral::bboxTP (Angle *as, Angle *ae, double *bbox)
{
  PTR<Angle> a1 = new SpiralXMaxAngle(this), a2 = new Rot90Angle(a1),
    a3 = new Rot90Angle(a2), a4 = new Rot90Angle(a3);
  SpiralData<Interval> sd = getApprox();
  double rax = max(fabs(sd.a.x.l), fabs(sd.a.x.u)),
    ray = max(fabs(sd.a.y.l), fabs(sd.a.y.u)),
    r = sqrt(rax*rax + ray*ray) + 1e-10;
  if (inInterval(a1, as, ae))
    bbox[1] = sd.b.x.u + r;
  if (inInterval(a2, as, ae))
    bbox[3] = sd.b.y.u + r;
  if (inInterval(a3, as, ae))
    bbox[0] = sd.b.x.l - r;
  if (inInterval(a4, as, ae))
    bbox[2] = sd.b.y.l - r;
}

void CVertex::outgoingHEdges (HCEdges &ed) const
{
  for (CEdges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    for (HCEdges::iterator f = (*e)->hedges.begin(); f != (*e)->hedges.end(); ++f)
      if ((*f)->tail() == this)
	ed.push_back(*f);
}

CVertexSet intersection (const CVertexSet &a, const CVertexSet &b)
{
  CVertexSet c;
  set_intersection(a.begin(), a.end(), b.begin(), b.end(),
		   inserter(c, c.begin()));
  return c;
}

CEdge::CEdge (CVertex *t, CVertex *h, bool flag) : t(t), h(h), flag(flag), s(0)
{
  inc = t->p->a->order(h->p->a);
  dec = !inc && h->p->a->order(t->p->a);
  copyBBox(t->bbox, bbox);
  mergeBBox(h->bbox, bbox);
  if (!flag)
    return;
  SpiralPoint *tp = dynamic_cast<SpiralPoint *>(t->getP()),
    *hp = dynamic_cast<SpiralPoint *>(h->getP());
  s = (tp && hp && tp->getS() == hp->getS()) ? tp->getS() : 0;
  if (s)
    if (tp->a->order(hp->a))
      s->bboxTP(tp->a, hp->a, bbox);
    else
      s->bboxTP(hp->a, tp->a, bbox);
}

CEdge::~CEdge ()
{
  for (HCEdges::iterator e = hedges.begin(); e != hedges.end(); ++e)
    delete *e;
}

void CEdge::addVertex (CVertex *v)
{
  if (find(vertices.begin(), vertices.end(), v) == vertices.end())
    vertices.push_back(v);
}

HCEdge * CEdge::addHEdge (bool forward)
{
  HCEdge *e = new HCEdge(this, forward);
  hedges.push_back(e);
  return e;
}

void CEdge::removeHEdge (HCEdge *e)
{
  HCEdges::iterator j = remove(hedges.begin(), hedges.end(), e);
  hedges.erase(j, hedges.end());
  delete e;
}

void CEdge::edgeVertices (CVertices &ve) const
{
  ve.push_back(t);
  ve.insert(ve.end(), vertices.begin(), vertices.end());
  ve.push_back(h);
}

CFace * CEdge::otherFace (CFace *f) const
{
  return hedges[0]->f == f ? hedges[1]->f : hedges[0]->f;
}

bool CEdge::piEdge () const
{
  return t->p->a == h->p->a && (t->p->a == Angle::mpi || t->p->a == Angle::ppi);
}

AngleInterval CEdge::angleInterval () const
{
  Angle *at = t->p->getA(), *ah = h->p->getA();
  if (at->order(ah))
    return AngleInterval(at, ah);
  return AngleInterval(ah, at);
}

CPoint * CEdge::point (Angle *a)
{
  map<Angle *, PTR<CPoint>>::iterator i = apmap.find(a);
  if (i != apmap.end())
    return i->second;
  CPoint *p = pointAux(a);
  apmap.insert(pair<Angle *, PTR<CPoint>>(a, p));
  return p;
}

CPoint * CEdge::pointAux (Angle *a) const
{
  if (s)
    return new SpiralPoint(a, s);
  HCEdge *hp = hedges[0]->hp;
  if (!hp)
    return new FFPoint(a, hedges[0]->f, hedges[1]->f);
  CEdge *ep = hp->e;
  if (ep->s)
    return new SpiralPoint(a, ep->s);
  return new FFPoint(a, ep->hedges[0]->f, ep->hedges[1]->f);
}

template<class N>
PV3<N> CEdge::getU (CPoint *p) const
{
  if (horizontal()) {
    PV2<N> uxy = h->getP()->get<N>() - t->getP()->get<N>();
    return PV3<N>(uxy.x, uxy.y, N(0));
  }
  if (s) {
    PV2<N> u = s->get<N>().a.rotate(p->a->get<N>());
    N dt = p->a->dt<N>(), vz = increasing() ? N(1) : N(-1);
    return PV3<N>(- vz*dt*u.y, vz*dt*u.x, vz);
  }
  return hedges[1]->f->getN<N>(p).cross(hedges[0]->f->getN<N>(p));
}

void CEdge::sortHEdges ()
{
  if (hedges.size() > 2)
    sort(hedges.begin(), hedges.end(), CEdgeOrder(this));
}

HCEdge * HCEdge::cw () const
{
  int n = e->hedges.size();
  for (int i = 0; i < n; ++i)
    if (e->hedges[i] == this)
      return e->hedges[(i+1)%n];
}

HCEdge * HCEdge::ccw () const
{
  int n = e->hedges.size();
  for (int i = 0; i < n; ++i)
    if (e->hedges[i] == this)
      return e->hedges[i == 0 ? n - 1 : i - 1];
}

void HCEdge::loop (HCEdges &ed)
{
  HCEdge *h = this;
  do {
    ed.push_back(h);
    h = h->next;
  }
  while (h != this);
}

bool HCEdge::onLoop ()
{
  HCEdge *h = this;
  do {
    h->flag = true;
    h = h->next;
    if (!h)
      return false;
  }
  while (h != this);
  return true;
}

template<class N>
PV3<N> HCEdge::getU (CPoint *p) const
{
  if (hp)
    return hp->getU<N>(p);
  return forward ? e->getU<N>(p) : - e->getU<N>(p);
}

template<class N>
PV3<N> HCEdge::getN (CPoint *p) const
{
  return forward ? f->getN<N>(p) : - f->getN<N>(p);
}

double * Slab::getBBox ()
{
  if (bbox[1] < bbox[0]) {
    PTR<CPoint> ls = l->point(s), rs = r->point(s), le = l->point(e), re = r->point(e);
    ls->getBBox(bbox);
    double bb[6];
    rs->getBBox(bb);
    mergeBBox(bb, bbox);
    le->getBBox(bb);
    mergeBBox(bb, bbox);
    re->getBBox(bb);
    mergeBBox(bb, bbox);
  }
  return bbox;
}

void Slab::addPoint (CPoint *p)
{
  if ((p->getA() == s || p->getA() == e) && inSlab(p, this))
    if (p->getA() == s)
      pb.push_back(p);
    else
      pt.push_back(p);
}

bool inSlab (CPoint *p, Slab *s)
{
  return p != s->l->getT()->getP() && p != s->l->getH()->getP() &&
    p != s->r->getT()->getP() && p != s->r->getH()->getP() &&
    InSlab(p, s) == 1;
}

CFace::CFace (Vertex *v, Face *f, bool aflag) : s(0)
{
  type = aflag ? CFaceVF : CFaceFV;
  feature1 = (int *) v;
  feature2 = (int *) f;
  HEdge *e = f->getBoundary();
  vs.insert(v);
  vs.insert(e->tail());
  vs.insert(e->head());
  vs.insert(e->getNext()->head());
}

CFace::CFace (Edge *e, Edge *f, bool pflag) : s(0)
{
  type = pflag ? CFaceEEP : CFaceEEM;
  feature1 = (int *) e;
  feature2 = (int *) f;
  vs.insert(e->getT());
  vs.insert(e->getH());
  vs.insert(f->getT());
  vs.insert(f->getH());
}

void CFace::addBoundary (HCEdge *e)
{
  HCEdges ed;
  e->loop(ed);
  if (boundary.empty()) {
    copyBBox(ed[0]->e->bbox, bbox);
    for (int i = 1; i < ed.size(); ++i)
      mergeBBox(ed[i]->e->bbox, bbox);
  }
  boundary.push_back(e);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    (*e)->f = this;
}

void CFace::vertices (VertexSet &vsa, VertexSet &vsb) const
{
  if (type == CFaceVF) {
    vsa.insert((Vertex *) feature1);
    HEdge *b = ((Face *) feature2)->getBoundary();
    vsb.insert(b->tail());
    vsb.insert(b->head());
    vsb.insert(b->getNext()->head());
  }
  else if (type == CFaceFV) {
    vsb.insert((Vertex *) feature1);
    HEdge *b = ((Face *) feature2)->getBoundary();
    vsa.insert(b->tail());
    vsa.insert(b->head());
    vsa.insert(b->getNext()->head());
  }
  else {
    Edge *ea = (Edge *) feature1, *eb = (Edge *) feature2;
    vsa.insert(ea->getT());
    vsa.insert(ea->getH());
    vsb.insert(eb->getT());
    vsb.insert(eb->getH());
  }
}

Spiral * CFace::sharedSpiral (CFace *f) const
{
  Spiral *sl = leftEdge()->s, *sr = rightEdge()->s,
    *fsl = f->leftEdge()->s, *fsr = f->rightEdge()->s;
  if (sl == fsl || sl == fsr)
    return sl;
  if (sr == fsl || sr == fsr)
    return sr;
  return 0;
}

void CFace::boundaryVertices (CVertices &ve) const
{
  HCEdges ed;
  boundaryHEdges(ed);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    ve.push_back((*e)->tail());
}

bool CFace::boundaryEdge (CEdge *e) const
{
  for (HCEdges::iterator f = e->hedges.begin(); f != e->hedges.end(); ++f)
    if ((*f)->f == this)
      return true;
  return false;
}

void CFace::boundaryHEdges (HCEdges &ed) const
{
  for (HCEdges::const_iterator b = boundary.begin(); b != boundary.end(); ++b)
    (*b)->loop(ed);
}

void CFace::sharedBoundaryVertices (CFace *f, CVertices &vfg) const
{
  CVertexSet vs1, vs2;
  edgeVertices(vs1);
  f->edgeVertices(vs2);
  CVertexSet vs = intersection(vs1, vs2);
  vfg.insert(vfg.end(), vs.begin(), vs.end());
}

void CFace::edgeVertices (CVertexSet &vs) const
{
  HCEdges ed;
  boundaryHEdges(ed);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e) {
    CVertices ve;
    (*e)->e->edgeVertices(ve);
    vs.insert(ve.begin(), ve.end());
  }
}

bool CFace::sameBoundary (CVertex *a, CVertex *b) const
{
  HCEdges ed;
  boundaryHEdges(ed);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e) {
    CVertices ve;
    (*e)->e->edgeVertices(ve);
    if (find(ve.begin(), ve.end(), a) != ve.end() &&
	find(ve.begin(), ve.end(), b) != ve.end())
      return true;
  }
  return false;
}

bool CFace::sharedBoundaryVertex (CFace *f, CFace *g) const
{
  CVertexSet vs1, vs2,
  edgeVertices(vs1);
  f->edgeVertices(vs2);
  CVertexSet vs12 = intersection(vs1, vs2);
  if (vs12.empty())
    return false;
  CVertexSet vs3;
  g->edgeVertices(vs3);
  CVertexSet vs123 = intersection(vs3, vs12);
  return !vs123.empty();
}

bool CFace::bboxOverlap (CFace *f, CFace *g) const
{
  double bb[6];
  for (int i = 0; i < 3; ++i) {
    bb[2*i] = max(bbox[2*i], f->bbox[2*i]);
    bb[2*i+1] = min(bbox[2*i+1], f->bbox[2*i+1]);
  }
  return ::bboxOverlap(bb, g->bbox);
}

bool CFace::angleOverlap (CFace *f, CFace *g) const
{
  Angle *s = startAngle(), *e = endAngle(), *fs = f->startAngle(),
    *fe = f->endAngle(), *gs = g->startAngle(), *ge = g->endAngle();
  if (s->order(fs))
    s = fs;
  if (s->order(gs))
    s = gs;
  if (fe->order(e))
    e = fe;
  if (ge->order(e))
    e = ge;
  return s->order(e);
}

bool CFace::inInterval (Angle *a) const
{
  return startAngle()->order(a) && a->order(endAngle());
}

// patch
bool CFace::contains (CPoint *p) const
{
  double bbp[6];
  p->getBBox(bbp);
  if (!::bboxOverlap(bbox, bbp))
    return false;
  CEdge *l = leftEdge(), *r = rightEdge();
  if (!inInterval(p->a))
    return false;
  SpiralPoint pl(p->a, l->getSpiral()), pr(p->a, r->getSpiral());
  return onEdge(p, &pl, &pr);
}

// general face
bool CFace::contains2 (CPoint *p) const
{
  HCEdges ed;
  boundaryHEdges(ed);
  CPoints pts;
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    if ((*e)->e->contains(p->a))
      pts.push_back((*e)->e->point(p->a));
  sort(pts.begin(), pts.end(), CPointOrderXO());
  bool res = false;
  for (int i = 0; !res && i + 1 < pts.size(); i += 2)
    res = onEdge(p, pts[i], pts[i+1]);
  return res;
}

template<class N>
PV3<N> CFace::getN (CPoint *p) const
{
  Angle *th = p->getA();
  N dth = th->dt<N>();
  PV2<N> pn = p->get<N>(), thn = th->get<N>();
  if (type == CFaceVF || type == CFaceFV) {
    Vertex *v = (Vertex *) feature1;
    Face *f = (Face *) feature2;
    PV3<N> vp = v->getP()->get<N>(), n = f->getP()->get<N>().n;
    PV2<N> v2(vp.x, vp.y), n2(n.x, n.y);
    if (type == CFaceVF)
      return PV3<N>(n.x, n.y, dth*v2.rotate(thn).cross(n2));
    PV2<N> dd = - n2.rotate(th->get<N>());
    N dt = dth*n2.rotate(thn).cross(v2 - pn);
    return PV3<N>(dd.x, dd.y, dt);
  }
  Edge *e = (Edge *) feature1, *f = (Edge *) feature2;
  PV3<N> a = e->getT()->getP()->get<N>(), u = e->getU<N>(),
    b = f->getT()->getP()->get<N>(), v = f->getU<N>(),
    dd = u.rotateZ(th->get<N>()).cross(v), bd(b.x - pn.x, b.y - pn.y, b.z),
    bdv = bd.cross(v), ua = u.cross(a);
  PV2<N> dtu = dth*PV2<N>(- u.y, u.x).rotate(thn),
    dtua = dth*PV2<N>(- ua.y, ua.x).rotate(thn);
  N dt = dtu.dot(PV2<N>(bdv.x, bdv.y)) - dtua.dot(PV2<N>(v.x, v.y));
  PV3<N> df(dd.x, dd.y, dt);
  return type == CFaceEEP ? - df : df;
}

// ((th(a1) + a2).(x, y) + k1*s + k2*c + k3 = 0
template<class N>
void CFace::coeffs (PV2<N> &a1, PV2<N> &a2, N &k1, N &k2, N &k3) const
{
  if (type == CFaceVF)
    coeffsVF(a1, a2, k1, k2, k3);
  else if (type == CFaceFV)
    coeffsFV(a1, a2, k1, k2, k3);
  else
    coeffsEE(a1, a2, k1, k2, k3);
}

template<class N>
void CFace::coeffsVF (PV2<N> &a1, PV2<N> &a2, N &k1, N &k2, N &k3) const
{
  Vertex *v = (Vertex *) feature1;
  Face *f = (Face *) feature2;
  PV3<N> p = v->getP()->get<N>(), a = f->getBoundary()->tail()->getP()->get<N>(),
    n = f->getP()->get<N>().n;
  a1.x = a1.y = N(0);
  a2.x = n.x;
  a2.y = n.y;
  k1 = p.x*n.y - p.y*n.x;
  k2 = p.x*n.x + p.y*n.y;
  k3 = p.z*n.z - a.dot(n);
}

template<class N>
void CFace::coeffsFV (PV2<N> &a1, PV2<N> &a2, N &k1, N &k2, N &k3) const
{
  Vertex *v = (Vertex *) feature1;
  Face *f = (Face *) feature2;
  PV3<N> p = v->getP()->get<N>(), a = f->getBoundary()->tail()->getP()->get<N>(),
    n = f->getP()->get<N>().n;
  a1.x = n.x;
  a1.y = n.y;
  a2.x = a2.y = N(0);
  k1 = p.x*n.y - p.y*n.x;
  k2 = - p.x*n.x - p.y*n.y;
  k3 = a.dot(n) - p.z*n.z;
}

template<class N>
void CFace::coeffsEE (PV2<N> &a1, PV2<N> &a2, N &k1, N &k2, N &k3) const
{
  Edge *e = (Edge *) feature1, *f = (Edge *) feature2;
  PV3<N> a = e->getT()->getP()->get<N>(), u = e->getH()->getP()->get<N>() - a,
    b = f->getT()->getP()->get<N>(), v = f->getH()->getP()->get<N>() - b,
    w = u.cross(a);
  a1.x = u.y*v.z;
  a1.y = - u.x*v.z;
  a2.x = - u.z*v.y;
  a2.y = u.z*v.x;
  k1 = v.x*w.y - v.y*w.x - v.z*(b.x*u.x + b.y*u.y) + b.z*(u.x*v.x + u.y*v.y);
  k2 = - v.x*w.x - v.y*w.y + v.z*(b.y*u.x - b.x*u.y) - b.z*(u.x*v.y - u.y*v.x);
  k3 = - v.z*w.z + u.z*(b.x*v.y - b.y*v.x);
}

template<class N>
void CFace::line (Angle *a, N &u, N &v, N &w) const
{
  PV2<N> a1, a2, q = a->get<N>();
  N k1, k2, k3;
  coeffs<N>(a1, a2, k1, k2, k3);
  PV2<N> n = a1.rotate(q) + a2;
  u = n.x;
  v = n.y;
  w = k1*q.y + k2*q.x + k3;
}

CFace * CFace::neighbor (HCEdge *e) const
{
  CEdge *ee = e->e;
  HCEdge *f = e->forward ? e->ccw() : e->cw();
  return f->forward == e->forward ? 0 : f->f;
}

void CFace::formSlabs ()
{
  HHCEdges hhe;
  CPoints hpts;
  initSlabs(hhe, hpts);
  Slabs sl;
  for (HHCEdges::iterator e = hhe.begin(); e != hhe.end(); ++e)
    switch (e->type()) {
    case RightEnd:
      updateRE(sl, *e);
      break;
    case LeftEnd:
      updateLE(sl, *e);
      break;
    case LeftStart:
      updateLS(sl, *e);
      break;
    case RightStart:
      updateRS(sl, *e);
    }
  for (CPoints::iterator p = hpts.begin(); p != hpts.end(); ++p)
    for (Slabs::iterator s = slabs.begin(); s != slabs.end(); ++s)
      (*s)->addPoint(*p);
  for (Slabs::iterator s = slabs.begin(); s != slabs.end(); ++s) {
    PTR<CPoint> bl = (*s)->l->point((*s)->s), br = (*s)->r->point((*s)->s),
      tl = (*s)->l->point((*s)->e), tr = (*s)->r->point((*s)->e);
    sort((*s)->pb.begin(), (*s)->pb.end(), EdgePointOrderL(bl, br));
    sort((*s)->pt.begin(), (*s)->pt.end(), EdgePointOrderL(tl, tr));
  } 
  for (Slabs::iterator s = sl.begin(); s != sl.end(); ++s)
    if (find(slabs.begin(), slabs.end(), *s) == slabs.end())
      delete *s;
}

void CFace::initSlabs (HHCEdges &hhe, CPoints &hpts) const
{
  HCEdges ed;
  boundary[0]->loop(ed);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e) {
    hpts.push_back((*e)->tail()->getP());
    if (!(*e)->e->horizontal()) {
      hhe.push_back(HHCEdge(*e, true));
      hhe.push_back(HHCEdge(*e, false));
    }
  }
  sort(hhe.begin(), hhe.end());
}

void CFace::updateRE (Slabs &sl, const HHCEdge &h)
{
  for (Slabs::iterator s = sl.begin(); s != sl.end(); ++s)
    if ((*s)->r == h.e->e) {
      if ((*s)->l) {
	CPoint *p = h.tail()->p;
	slabs.push_back(new Slab((*s)->l, (*s)->r, (*s)->s, p->a, this));
      }
      (*s)->r = 0;
      return;
    } 
}

void CFace::updateLE (Slabs &sl, const HHCEdge &h)
{
  for (Slabs::iterator s = sl.begin(); s != sl.end(); ++s)
    if ((*s)->l == h.e->e) {
      if ((*s)->r) {
	CPoint *p = h.tail()->p;
	slabs.push_back(new Slab((*s)->l, (*s)->r, (*s)->s, p->a, this));
	for (Slabs::iterator t = sl.begin(); t != sl.end(); ++t)
	  if ((*t)->l && !(*t)->r) {
	    (*t)->r = (*s)->r;
	    (*t)->s = p->a;
	    (*t)->pb = (*s)->pb;
	    (*t)->pt = (*s)->pt;
	    (*s)->r = 0;
	    break;
	  }
      }
      (*s)->l = 0;
      return;
    }
}

void CFace::updateLS (Slabs &sl, const HHCEdge &h)
{
  CPoint *p = h.tail()->p;
  Angle *a = p->a;
  for (Slabs::iterator s = sl.begin(); s != sl.end(); ++s)
    if (!(*s)->l && (*s)->r) {
      (*s)->l = h.e->e;
      (*s)->s = a;
      return;
    }
  for (Slabs::iterator s = sl.begin(); s != sl.end(); ++s)
    if ((*s)->l && (*s)->r && inSlab(p, *s)) {
      slabs.push_back(new Slab((*s)->l, (*s)->r, (*s)->s, a, this));
      CEdge *osl = (*s)->l;
      (*s)->l = h.e->e;
      (*s)->s = a;
      sl.push_back(new Slab(osl, 0, a, 0, this));
      return;
    }
  sl.push_back(new Slab(h.e->e, 0, a, 0, this));
}

void CFace::updateRS (Slabs &sl, const HHCEdge &h)
{
  CPoint *p = h.tail()->p;
  Angle *a = p->a;
  Slabs cand;
  for (Slabs::iterator s = sl.begin(); s != sl.end(); ++s)
    if ((*s)->l && (!(*s)->r || inSlab(p, *s)))
      cand.push_back(*s);
  Slab *sn = 0;
  if (cand.size() == 1)
    sn = cand[0];
  else
    for (Slabs::iterator s = cand.begin(); s != cand.end(); ++s)
      if (slabNextIncreasing(*s) == h.e->e) {
	sn = *s;
	break;
      }
  if (sn->r)
    slabs.push_back(new Slab(sn->l, sn->r, sn->s, a, this));
  sn->r = h.e->e;
  sn->s = a;
}

CEdge * CFace::slabNextIncreasing (Slab *s) const
{
  HCEdge *e = s->l->hedges[0]->f == this ? s->l->hedges[0] : s->l->hedges[1];
  e = e->next;
  while (e->e->horizontal())
    e = e->next;
  return e->e;
}

void CFace::splitSlab (Slab *s, const Angles &an, Slabs &sl)
{
  Angle *as = s->s;
  for (int i = 0; i <= an.size(); ++i)
    if (i < an.size()) {
      sl.push_back(new Slab(s->l, s->r, as, an[i], this));
      as = an[i];
    }
    else
      sl.push_back(new Slab(s->l, s->r, as, s->e, this));
  sl[0]->pb = s->pb;
  (*sl.rbegin())->pt = s->pt;
  Slabs::iterator iter = remove(slabs.begin(), slabs.end(), s);
  slabs.erase(iter, slabs.end());
  delete s;
  slabs.insert(slabs.end(), sl.begin(), sl.end());
}

Angles anglesSF (Spiral *s, CFace *f)
{
  PTR<Object<SinCosData>> x = new PolySF(s, f);
  return sinCosAngles(x, f->startAngle(), f->endAngle(), VertexSet());
}

Angles anglesFFF (CFace *f1, CFace *f2, CFace *f3)
{
  Angle *s1 = f1->startAngle(), *s2 = f2->startAngle(), *s3 = f3->startAngle(),
    *e1 = f1->endAngle(), *e2 = f2->endAngle(), *e3 = f3->endAngle(),
    *s = s1->order(s2) ? s2 : s1, *e = e1->order(e2) ? e1 : e2;
  if (s->order(s3))
    s = s3;
  if (e3->order(e))
    e = e3;
  if (s->order(e)) {
    PTR<Object<SinCosData>> x = new PolyFFF(f1, f2, f3);
    return sinCosAngles(x, s, e, VertexSet());
  }
}

void CShell::init ()
{
  CVertexSet vs;
  for (CFaces::iterator f = faces.begin(); f != faces.end(); ++f) {
    (*f)->s = this;
    CVertices vf;
    (*f)->boundaryVertices(vf);
    vs.insert(vf.begin(), vf.end());
  }
  CVertexSet::iterator v = vs.begin();
  vm = *v;
  ++v;
  while (v != vs.end()) {
    if (CPointOrderX((*v)->getP(), vm->getP()) == 1)
      vm = *v;
    ++v;
  }
}

CCell::~CCell ()
{
  for (CShells::iterator s = shells.begin(); s != shells.end(); ++s)
    delete *s;
}

Cspace::~Cspace ()
{
  if (par)
    delete par;
  for (CVertices::iterator v = vertices.begin(); v != vertices.end(); ++v)
    delete *v;
  for (CEdges::iterator e = edges.begin(); e != edges.end(); ++e)
    delete *e;
  for (CFaces::iterator f = faces.begin(); f != faces.end(); ++f)
    delete *f;
  for (CCells::iterator c = cells.begin(); c != cells.end(); ++c)
    delete *c;
  for (SpiralSet::iterator s = ss.begin(); s != ss.end(); ++s)
    delete *s;
}

Angles Cspace::angles (Edge *e, Face *f, bool aflag)
{
  EFKey ef(e, f);
  EFAMap::iterator i = efamap.find(ef);
  if (i != efamap.end())
    return i->second;
  Angles a = anglesEF(e, f, aflag);
  efamap.insert(EFAPair(ef, a));
  return a;
}

Angles Cspace::angles (const EFKey &k, bool aflag)
{
  EFAMap::iterator i = efamap.find(k);
  if (i != efamap.end())
    return i->second;
  Angles a = anglesSFL(k, aflag);
  efamap.insert(EFAPair(k, a));
  return a;
}

void Cspace::angleIntervals (Edge *e, Face *f, bool aflag, AngleIntervals &pos,
			     AngleIntervals &neg)
{
  Angles a = angles(e, f, aflag);
  bool amin = PolyEFMin(e, f) == 1;
  if (!a.empty()) {
    if (amin) {
      pos.push_back(AngleInterval(Angle::mpi, a[0]));
      pos.push_back(AngleInterval(a[1], Angle::ppi));
      neg.push_back(AngleInterval(a[0], a[1]));
    }
    else {
      pos.push_back(AngleInterval(a[0], a[1]));
      neg.push_back(AngleInterval(Angle::mpi, a[0]));
      neg.push_back(AngleInterval(a[1], Angle::ppi));
    }
  }
  else if (amin)
    pos.push_back(AngleInterval(Angle::mpi, Angle::ppi));
  else
    neg.push_back(AngleInterval(Angle::mpi, Angle::ppi));
}

void Cspace::angleIntervals (HEdge *e, Face *f, bool aflag, AngleIntervals &pos,
			     AngleIntervals &neg)
{
  if (e->getForward())
    angleIntervals(e->getE(), f, aflag, pos, neg);
  else
    angleIntervals(e->getE(), f, aflag, neg, pos);
}

Spiral * Cspace::getSpiral (Vertex *v, Edge *e, bool aflag)
{
  Spiral *s = new Spiral(v, e, aflag);
  pair<SpiralSet::iterator, bool> x = ss.insert(s);
  if (x.second)
    return s;
  delete s;
  return *x.first;
}

CVertex * Cspace::getVertex (CPoint *p)
{
 CVertex *v = new CVertex(p);
  if (vertices.empty())
    copyBBox(v->bbox, bbox);
  else
    mergeBBox(v->bbox, bbox);
  vertices.push_back(v);
  return v;
}

CVertex * Cspace::getVertex (Angle *a, Spiral *s)
{
  ASPair as(a, s);
  ASVMap::iterator i = asvmap.find(as);
  if (i != asvmap.end())
    return i->second;
  SpiralPoint *p = new SpiralPoint(a, s);
  CVertex *v = getVertex(p);
  asvmap.insert(ASVPair(as, v));
  return v;
}

HCEdge * Cspace::addHEdge (CVertex *a, CVertex *b, bool flag, HCEdge *hp)
{
  CEdge *e = getEdge(a, b, flag, hp);
  bool forward = e->t == a;
  HCEdge *he = e->addHEdge(forward);
  he->hp = hp;
  return he;
}

CEdge * Cspace::getEdge (CVertex *a, CVertex *b, bool flag, HCEdge *hp)
{
  for (CEdges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if (((*e)->t == a && (*e)->h == b || (*e)->t == b && (*e)->h == a) &&
	(*e)->flag == flag && !(*e)->hedges.empty() &&
	(flag || !hp ||hp->e == (*e)->hedges[0]->hp->e))
      return *e;
  for (CEdges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if (((*e)->t == a && (*e)->h == b || (*e)->t == b && (*e)->h == a) &&
	(*e)->hedges.empty()) {
      (*e)->flag = flag;
      return *e;
    }
  CEdge *e = new CEdge(a, b, flag);
  edges.push_back(e);
  a->edges.push_back(e);
  b->edges.push_back(e);
  return e;
}

HCEdge * Cspace::addLoop (const CVertices &ve)
{
 int n = ve.size();
 HCEdge *ep = addHEdge(ve[n-1], ve[0], true), *e0 = ep;
 for (int i = 0; i + 1 < n; ++i) {
   HCEdge *e = addHEdge(ve[i], ve[i+1], true);
   ep->next = e;
   ep = e;
  }
 ep->next = e0;
 return e0->next;
}

void Cspace::patches (Polyhedron *a, Polyhedron *b)
{
  patchesVF(a, b, true);
  patchesVF(b, a, false);
  patchesEE(a, b);
}

void Cspace::patchesVF (Polyhedron *a, Polyhedron *b, bool aflag)
{
  for (Vertices::iterator v = a->vertices.begin(); v != a->vertices.end(); ++v) {
    HEdges ed;
    if (convexCone(*v, ed))
      for (Faces::iterator f = b->faces.begin(); f != b->faces.end(); ++f)
	patchVF(*v, *f, aflag);
  }
}

void Cspace::patchVF (Vertex *v, Face *f, bool aflag)
{
  double *vb = v->getBBox(), *fb = f->getBBox();
  if (vb[5] <= fb[4] || fb[5] <= vb[4])
    return;
  HEdges he = v->outgoingHEdges();
  AngleIntervals ais, aisd;
  angleIntervals(he[0], f, aflag, ais, aisd);
  for (int i = 1; !ais.empty() && i < he.size(); ++i) {
    AngleIntervals epos, eneg;
    angleIntervals(he[i], f, aflag, epos, eneg);
    ais = intersect(ais, epos);
  }
  for (AngleIntervals::iterator i = ais.begin(); i != ais.end(); ++i)
    patchVF(v, f, aflag, i->first, i->second);
}

void Cspace::patchVF (Vertex *v, Face *f, bool aflag, Angle *s, Angle *e)
{
  Spiral *l = 0, *r = 0;
  HEdges he = f->getBoundary()->edgeLoop();
  for (HEdges::iterator h = he.begin(); !r && h != he.end(); ++h) {
    Edge *hh = (*h)->getE();
    if (inIntervalZ(v, hh)) {
      Spiral *sp = getSpiral(v, hh, aflag);
      if (!l)
	l = sp;
      else
	r = sp;
    }
  }
  if (OrientationVF(aflag, f, s, l, r) == -1) {
    Spiral *temp = l;
    l = r;
    r = temp;
  }
  CFace *p = new CFace(v, f, aflag);
  faces.push_back(p);
  patch(s, e, l, r, p);
}

void Cspace::patch (Angle *s, Angle *e, Spiral *l, Spiral *r, CFace *f)
{
  CVertices ve;
  ve.push_back(getVertex(s, l));
  ve.push_back(getVertex(s, r));
  ve.push_back(getVertex(e, r));
  ve.push_back(getVertex(e, l));
  HCEdge *h = addLoop(ve);
  f->addBoundary(h);
  l->edges.push_back(f->leftEdge());
  r->edges.push_back(f->rightEdge());
}

void Cspace::patchesEE (Polyhedron *a, Polyhedron *b)  
{
  Edges aedges, bedges;
  for (Edges::iterator e = a->edges.begin(); e != a->edges.end(); ++e)
    if (convexEdge(*e) == 1)
      aedges.push_back(*e);
  for (Edges::iterator e = b->edges.begin(); e != b->edges.end(); ++e)
    if (convexEdge(*e) == 1)
      bedges.push_back(*e);
  for (Edges::iterator e = aedges.begin(); e != aedges.end(); ++e)
    for (Edges::iterator f = bedges.begin(); f != bedges.end(); ++f)
      patchEE(*e, *f);
}

void Cspace::patchEE (Edge *e, Edge *f)
{
  double *eb = e->getBBox(), *fb = f->getBBox();
  if (eb[5] <= fb[4] || fb[5] <= eb[4])
    return;
  AngleIntervals pos1, neg1, pos2, neg2, pos3, neg3, pos4, neg4;
  angleIntervals(f, e->getHEdge(0)->getF(), false, pos1, neg1);
  angleIntervals(f, e->getHEdge(1)->getF(), false, pos2, neg2);
  angleIntervals(e, f->getHEdge(0)->getF(), true, pos3, neg3);
  angleIntervals(e, f->getHEdge(1)->getF(), true, pos4, neg4);
  AngleIntervals pos = intersect(intersect(intersect(pos1, neg2), pos3), neg4),
    neg = intersect(intersect(intersect(neg1, pos2), neg3), pos4);
  for (AngleIntervals::iterator i = pos.begin(); i != pos.end(); ++i)
    patchEE(e, f, i->first, i->second, true);
  for (AngleIntervals::iterator i = neg.begin(); i != neg.end(); ++i)
    patchEE(e, f, i->first, i->second, false);
}

void Cspace::patchEE (Edge *e1, Edge *e2, Angle *s, Angle *e, bool pflag)
{
  Spiral *l = 0, *r = 0;
  Vertex *vs[] = {e1->getT(), e1->getH(), e2->getT(), e2->getH()};
  Edge *es[] = {e2, e2, e1, e1};
  bool fs[] = {true, true, false, false};
  for (int i = 0; i < 4; ++i)
    if (inIntervalZ(vs[i], es[i])) {
      Spiral *sp = getSpiral(vs[i], es[i], fs[i]);
      if (!l)
	l = sp;
      else
	r = sp;
    }
  if (OrientationEE(pflag, e1, e2, s, l, r) == -1) {
    Spiral *temp = l;
    l = r;
    r = temp;
  }
  CFace *p = new CFace(e1, e2, pflag);
  faces.push_back(p);
  patch(s, e, l, r, p);
}

void Cspace::intersectFF ()
{
  Octree<CFace *> *octree = Octree<CFace *>::octree(faces, bbox);
  vector<pair<CFace *, CFace *>> ff;
  octree->pairs(ff);
  CEEPairSet eps;
  for (vector<pair<CFace *, CFace *>>::iterator i = ff.begin(); i != ff.end(); ++i) {
    intersectEE(i->first, i->second, eps);
    intersectFF(i->first, i->second);
  }
  delete octree;
}

void Cspace::intersectEE (CFace *f, CFace *g, CEEPairSet &eps)
{
  HCEdges ef, eg;
  f->boundaryHEdges(ef);
  g->boundaryHEdges(eg);
  for (HCEdges::iterator i = ef.begin(); i != ef.end(); ++i) {
    CEdge *fi = (*i)->e;
    for (HCEdges::iterator j = eg.begin(); j != eg.end(); ++j) {
      CEdge *gj = (*j)->e;
      if (bboxOverlap(fi->bbox, gj->bbox) &&
	  eps.insert(CEEPair(fi < gj ? fi : gj, fi < gj ? gj : fi)).second)
	intersectEE(fi, gj);
    }
  }
}

void Cspace::intersectEE (CEdge *e, CEdge *f)
{
  Spiral *se = e->getSpiral(), *sf = f->getSpiral();
  if (se) {
    if (se == sf)
      intersectSS(e, f);
    else if (!sf)
      intersectSL(e, f);
  }
  else if (sf)
    intersectSL(f, e);
  else
    intersectLL(e, f);
}

void Cspace::intersectSS (CEdge *e, CEdge *f) const
{
  if (e->contains(f->t->p->a))
    e->addVertex(f->t);
  if (e->contains(f->h->p->a))
    e->addVertex(f->h);
  if (f->contains(e->t->p->a))
    f->addVertex(e->t);
  if (f->contains(e->h->p->a))
    f->addVertex(e->h);
}

void Cspace::intersectSL (CEdge *e, CEdge *f)
{
  if (f->piEdge())
    return;
  Angle *eta = e->t->p->a, *eha = e->h->p->a, *fa = f->t->p->a;
  if (fa == eta && f->contains(e->t->p))
    f->addVertex(e->t);
  else if (fa == eha && f->contains(e->h->p))
    f->addVertex(e->h);
  else if (inInterval(fa, eta, eha) &&
	   e->getSpiral()->sharedVertices(fa) == 3) {
    CVertex *v = getVertex(fa, e->getSpiral());
    if (f->contains(v->p)) {
      e->addVertex(v);
      f->addVertex(v);
    }
  }
}

void Cspace::intersectLL (CEdge *e, CEdge *f)
{
  Angle *a = e->t->p->a;
  if (a != f->t->p->a)
    return;
  if (a == Angle::mpi || a == Angle::ppi) {
    CFace *ef = e->hedges[0]->f, *ff = f->hedges[0]->f;
    if (IndependentFaces(ef, ff, a) == 0)
      return;
    PTR<CPoint> p = new FFPoint(a, ef, ff);
    if (e->contains(p) && f->contains(p)) {
      CVertex *v = getVertex(p);
      e->addVertex(v);
      f->addVertex(v);
    }
    return;
  }
  if (e->contains(f->t->p))
    e->addVertex(f->t);
  if (e->contains(f->h->p))
    e->addVertex(f->h);
  if (f->contains(e->t->p))
    f->addVertex(e->t);
  if (f->contains(e->h->p))
    f->addVertex(e->h);
}

void Cspace::intersectFF (CFace *f, CFace *g)
{
  CVertices vfg;
  intersectFF(f, g, vfg);
  intersectFF(g, f, vfg);
  formFF(f, g, vfg);
}

void Cspace::intersectFF (CFace *f, CFace *g, CVertices &vfg)
{
  HCEdges eg;
  g->boundaryHEdges(eg);
  for (HCEdges::iterator e = eg.begin(); e != eg.end(); ++e)
    intersectEF(*e, f, vfg);
}

void Cspace::intersectEF (HCEdge *e, CFace *f, CVertices &vfg)
{
  CEdge *ee = e->e;
  if (!bboxOverlap(ee->bbox, f->bbox))
    return;
  AngleInterval eai = ee->angleInterval();
  if (!intervalOverlap(eai.first, eai.second, f->startAngle(), f->endAngle()))
    return;
  CEFPair ef(ee, f);
  CEFVMap::iterator i = efvmap.find(ef);
  if (i != efvmap.end()) {
    vfg.insert(vfg.end(), i->second.begin(), i->second.end());
    return;
  }
  CVertices ve;
  intersectEF(ee, f, ve);
  efvmap.insert(CEFVPair(ef, ve));
  vfg.insert(vfg.end(), ve.begin(), ve.end());
}

void Cspace::intersectEF (CEdge *e, CFace *f, CVertices &ve)
{
  Spiral *s = e->getSpiral();
  if (!s) {
    CVertex *v = intersectEFLine(e, f);
    if (v)
      ve.push_back(v);
    return;
  }
  if (s == f->leftEdge()->s || s == f->rightEdge()->s)
    return;
  Angles asf = intersectSF(s, f);
  if (!asf.empty()) {
    AngleInterval eai = e->angleInterval();
    Angle *fs = f->startAngle(), *fe = f->endAngle(),
      *as = eai.first->order(fs) ? fs : eai.first,
      *ae = eai.second->order(fe) ? eai.second : fe;
    for (int i = 0; i < asf.size(); ++i) {
      Angle *a = asf[i];
      if (inInterval(a, as, ae)) {
	SpiralPoint p(a, s);
	if (f->contains(&p)) {
	  CVertex *v = getVertex(a, s);
	  e->vertices.push_back(v);
	  ve.push_back(v);
	}
      }
    }
  }
}

CVertex * Cspace::intersectEFLine (CEdge *e, CFace *f)
{
  Angle *a = e->t->p->a;
  AFPair af(a, f);
  CVertex *v = 0;
  AFVMap::iterator i = afvmap.find(af);
  if (i != afvmap.end())
    v = i->second;
  else if (f->sharedVertices(a) < 3) {
    PTR<CPoint> p = new LinePatchPoint(e, f);
    if (f->contains(p))
      v = getVertex(p);
    afvmap.insert(AFVPair(af, v));
  }
  if (v && e->contains(v->p)) {
    e->addVertex(v);
    return v;
  }
  return 0;
}

Angles Cspace::intersectSF (Spiral *s, CFace *f)
{
  Angles al = intersectSFLine(s, f);
  if (!al.empty())
    return al;
  SFPair sf(s, f);
  SFAMap::iterator i = sfamap.find(sf);
  if (i != sfamap.end())
    return i->second;
  Angles a = anglesSF(s, f);
  sfamap.insert(SFAPair(sf, a));
  return a;
}

Angles Cspace::intersectSFLine (Spiral *s, CFace *f)
{
  VertexSet vsa, vsb;;
  s->vertices(vsa, vsb);
  f->vertices(vsa, vsb);
  if (vsa.size() == 2 && vsb.size() == 3) {
    VertexSet::iterator i = vsb.begin();
    Vertex *v1 = *vsa.begin(), *v2 = *vsa.rbegin(), *v3 = *i, *v5 = *vsb.rbegin();
    ++i;
    EFKey k(v1, v2, v3, *i, v5);
    return angles(k, true);
  }
  if (vsa.size() == 3 && vsb.size() == 2) {
    VertexSet::iterator i = vsa.begin();
    Vertex *v1 = *vsb.begin(), *v2 = *vsb.rbegin(), *v3 = *i, *v5 = *vsa.rbegin();
    ++i;
    EFKey k(v1, v2, v3, *i, v5);
    return angles(k, false);
  }
  return Angles();
}

void Cspace::formFF (CFace *f, CFace *g, CVertices &vfg)
{
  if (intersectsFFLine(f, g))
    formFFLine(f, g, vfg);
  else
    formFFCurve(f, g, vfg);
}

bool Cspace::intersectsFFLine (CFace *f, CFace *g) const
{
  VertexSet vsa, vsb;
  f->vertices(vsa, vsb);
  g->vertices(vsa, vsb);
  return vsa.size() == 2 && vsb.size() == 3 || vsa.size() == 3 && vsb.size() == 2;
}

void Cspace::formFFLine (CFace *f, CFace *g, CVertices &vfg)
{
  sort(vfg.begin(), vfg.end(), VertexAngleOrder());
  int i = 0, n = vfg.size();
  while (i < n)
    if (i + 1 < n && vfg[i]->p->a == vfg[i+1]->p->a) {
      formFF(f, g, vfg[i], vfg[i+1]);
      i += 2;
    }
    else {
      Spiral *s = f->sharedSpiral(g);
      formFF(f, g, vfg[i], spiralVertex(vfg[i]->p->a, s));
      ++i;
    }
}

CVertex * Cspace::spiralVertex (Angle *a, Spiral *s)
{
  int nv = vertices.size();
  CVertex *v = getVertex(a, s);
  if (nv < vertices.size())
    for (CEdges::iterator e = s->edges.begin(); e != s->edges.end(); ++e)
      if ((*e)->contains(v->p))
	(*e)->addVertex(v);
  return v;
}

void Cspace::formFFCurve (CFace *f, CFace *g, CVertices &vfg)
{
  int n = vfg.size();
  if (n != 2 && n != 4)
    f->sharedBoundaryVertices(g, vfg);
  sort(vfg.begin(), vfg.end(), VertexAngleOrder());
  for (int i = 0; i + 1 < vfg.size(); i += 2)
    if (!(n == 0 && f->sameBoundary(vfg[i], vfg[i+1]) &&
	  g->sameBoundary(vfg[i], vfg[i+1])))
      formFF(f, g, vfg[i], vfg[i+1]);
}

void Cspace::formFF (CFace *f, CFace *g, CVertex *v, CVertex *w)
{
  CEdge *e = CFFOrder(f, g, v, w) == 1 ? getEdge(v, w, false)
    : getEdge(w, v, false);
  f->edges.push_back(e);
  e->addHEdge(true)->f = f;
  g->edges.push_back(e);
  e->addHEdge(false)->f = g;
} 

void Cspace::intersectFFF ()
{
  CFFEMap ffemap;
  formFFEMap(ffemap);
  for (CFaces::iterator f = faces.begin(); f != faces.end(); ++f)
    intersectFFF(*f, ffemap);
}

void Cspace::formFFEMap (CFFEMap &ffemap) const
{

  CEdges::const_iterator e = edges.begin();
  while (e != edges.end() && (*e)->flag)
    ++e;
  while (e != edges.end()) {
    CFace *f = (*e)->hedges[0]->f, *g = (*e)->hedges[1]->f;
    CFFPair fg(f < g ? f : g, f < g ? g : f);
    CFFEMap::iterator i = ffemap.find(fg);
    if (i == ffemap.end()) {
      CEdges ed;
      ed.push_back(*e);
      ffemap.insert(CFFEPair(fg, ed));
    }
    else
      i->second.push_back(*e);
    ++e;
  }
}

void Cspace::intersectFFF (CFace *f, const CFFEMap &ffemap)
{
  int n = f->edges.size();
  for (int i = 0; i + 1 < n; ++i) {
    CEdge *ei = f->edges[i];
    CFace *fi = ei->otherFace(f);
    if (f < fi)
      for (int j = i + 1; j < n; ++j) {
	CEdge *ej = f->edges[j];
	CFace *fj = ej->otherFace(f);
	if (f < fj && fi != fj)
	  intersectFFF(f, fi, fj, ei, ej, ffemap);
      }
  }
}

void Cspace::intersectFFF (CFace *f1, CFace *f2, CFace *f3, CEdge *e12, CEdge *e13,
			   const CFFEMap &ffemap)
{
  if (!(f1->bboxOverlap(f2, f3) && f1->angleOverlap(f2, f3)))
    return;
  CFFPair f23(f2 < f3 ? f2 : f3, f2 < f3 ? f3 : f2);
  CFFEMap::const_iterator iter = ffemap.find(f23);
  if (iter == ffemap.end())
    return;
  const CEdges ed = iter->second;
  if (e12->horizontal())
      for (CEdges::const_iterator e23 = ed.begin(); e23 != ed.end(); ++e23)
	intersectFFFH(f1, f2, f3, e12, e13, *e23);
  else if (e13->horizontal())
    for (CEdges::const_iterator e23 = ed.begin(); e23 != ed.end(); ++e23)
      intersectFFFH(f1, f3, f2, e13, e12, *e23);
  else if (ed[0]->horizontal())
    for (CEdges::const_iterator e23 = ed.begin(); e23 != ed.end(); ++e23)
      intersectFFFH(f2, f3, f1, *e23, e12, e13);
  else
    intersectFFFG(f1, f2, f3, e12, e13, ed);
}

void Cspace::intersectFFFH (CFace *f1, CFace *f2, CFace *f3, CEdge *eh,
			    CEdge *e1, CEdge *e2)
{
  Angle *a = eh->t->p->a;
  if (e1->contains(a) && e2->contains(a)) {
    CVertex *v = intersectEFLine(eh, f3);
    if (v && onEdge(v->p, eh->t->p, eh->h->p)) {
      eh->addVertex(v);
      e1->addVertex(v);
      e2->addVertex(v);
    }
  }
}

void Cspace::intersectFFFG (CFace *f1, CFace *f2, CFace *f3, CEdge *e12,
			    CEdge *e13, const CEdges &ed)
{
  Angles a = anglesFFF(f1, f2, f3);
  if (a.size() == 0)
    return;
  for (int i = 0; i < a.size(); ++i) {
    Angle *ai = a[i];
    if (e12->contains(ai) && e13->contains(ai))
      for (CEdges::const_iterator e23 = ed.begin(); e23 != ed.end(); ++e23)
	if ((*e23)->contains(ai)) {
	  PTR<CPoint> p = new FFPoint(ai, f1, f2);
	  CVertex *v = getVertex(p);
	  e12->addVertex(v);
	  e13->addVertex(v);
	  (*e23)->addVertex(v);
	  break;
	}
  }
}

void Cspace::sortVertices ()
{
  for (CEdges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    if (!(*e)->vertices.empty())
      if ((*e)->horizontal())
	sort((*e)->vertices.begin(), (*e)->vertices.end(), EdgeVertexOrderL(*e));
      else
	sort((*e)->vertices.begin(), (*e)->vertices.end(), EdgeVertexOrderS(*e));
}

Cspace * Cspace::subfaces ()
{
  Cspace *a = new Cspace(this);
  CVVMap vvmap;
  for (CFaces::const_iterator f = faces.begin(); f != faces.end(); ++f)
    subfaces(*f, a, vvmap);
  return a;
}

void Cspace::subfaces (CFace *f, Cspace *a, CVVMap &vvmap) const
{
  HCEdges he, outer, bad;
  subedges(f, a, vvmap, he);
  for (HCEdges::iterator e = he.begin(); e != he.end(); ++e)
    if (!(*e)->flag)
      if ((*e)->onLoop())
	outer.push_back(*e);
      else
	bad.push_back(*e);
  for (HCEdges::iterator e = outer.begin(); e != outer.end(); ++e)
    a->addFace(*e, f);
  a->removeBad(bad);
}

void Cspace::subedges (CFace *f, Cspace *a, CVVMap &vvmap, HCEdges &he) const
{
  HCEdges ed;
  f->boundaryHEdges(ed);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    subedges(*e, a, vvmap, he);
  for (CEdges::iterator e = f->edges.begin(); e != f->edges.end(); ++e) {
    HCEdge *h = (*e)->hedges[0]->f == f ? (*e)->hedges[0] : (*e)->hedges[1];
    subedges(h, a, vvmap, he);
  }
  setNext(f, he);
}

void Cspace::subedges (HCEdge *e, Cspace *a, CVVMap &vvmap, HCEdges &he) const
{
  CVertices ve;
  e->e->edgeVertices(ve);
  CVertex *t = a->getVertex(ve[0], vvmap);
  for (int i = 1; i < ve.size(); ++i) {
    CVertex *h = a->getVertex(ve[i], vvmap);
    if (e->forward)
      he.push_back(a->addHEdge(t, h, e->e->flag, e));
    else
      he.push_back(a->addHEdge(h, t, e->e->flag, e));
    t = h;
  }
}

CVertex * Cspace::getVertex (CVertex *v, CVVMap &vmap)
{
  CVVMap::iterator iter = vmap.find(v);
  if (iter != vmap.end())
    return iter->second;
  CVertex *w = getVertex(v->p);
  vmap.insert(CVVPair(v, w));
  return w;
}

void Cspace::setNext (CFace *f, const HCEdges &he) const
{
  HHCEdges hhe;
  for (HCEdges::const_iterator e = he.begin(); e != he.end(); ++e) {
    hhe.push_back(HHCEdge(*e, true));
    hhe.push_back(HHCEdge(*e, false));
  }
  sort(hhe.begin(), hhe.end(), HHCEdgeOrder(f));
  int i = 0, n = hhe.size();
  while (i < n) {
    int m = 1;
    while (i + m < n && hhe[i].tail() == hhe[i+m].tail())
      ++m;
    int j = 0;
    while (j < m) {
      int k = i + (j + 1)%m;
      if (!hhe[i+j].f && hhe[k].f)
	hhe[i+j].e->next = hhe[k].e;
      ++j;
    }
    i += m;
  }
}

CFace * Cspace::addFace (HCEdge *e, CFace *f)
{
  CFace *g = new CFace(f->type, f->feature1, f->feature2);
  faces.push_back(g);
  g->addBoundary(e);
  return g;
}

void Cspace::removeBad (const HCEdges &ed)
{
  HCEdgeSet es;
  for (HCEdges::const_iterator e = ed.begin(); e != ed.end(); ++e) {
    HCEdge *f = *e;
    while (f) {
      es.insert(f);
      f = f->next;
      if (f == *e)
	break;
    }
  }
  for (HCEdgeSet::iterator e = es.begin(); e != es.end(); ++e)
    (*e)->e->removeHEdge(*e);
}

void Cspace::formCells (Polyhedron *a, Polyhedron *b)
{
  CShells sh;
  formShells(a, b, sh);
  sort(sh.begin(), sh.end(), CShellOrder());
  Octree<CFace *> *octree = Octree<CFace *>::octree(faces, bbox);
  cells.push_back(new CCell);
  for (int i = 0; i < sh.size(); ++i) {
    CShell *s = sh[i];
    CShell *t = i == 0 ? 0 : enclosingShell(s, octree);
    if (!t) {
      s->c = cells[0];
      cells[0]->shells.push_back(s);
    }
    else if (t->outer) {
      s->c = t->c;
      t->c->shells.push_back(s);
    }
    else {
      s->outer = true;
      s->c = new CCell;
      s->c->shells.push_back(s);
      cells.push_back(s->c);
    }
  }
  delete octree;
}

void Cspace::formShells (Polyhedron *a, Polyhedron *b, CShells &sh)
{
  CShells sh0;
  formShells(sh0);
  for (CShells::iterator s = sh0.begin(); s != sh0.end(); ++s) {
    CPoint *c = (*s)->faces[0]->boundary[0]->tail()->p;
    Polyhedron *ac = transform(a, c);
    bool bad = ac->intersects(b, true);
    delete ac;
    if (!bad) {
      (*s)->init();
      sh.push_back(*s);
    }
    else
      delete *s;
  }
  int i = 0;
  while (i < faces.size())
    if (faces[i]->s)
      ++i;
    else {
      removeFace(faces[i]);
      faces[i] = *faces.rbegin();
      faces.pop_back();
    }
}

void Cspace::formShells (CShells &sh) const
{
  for (CEdges::const_iterator e = edges.begin(); e != edges.end(); ++e)
    (*e)->sortHEdges();
  set<CFace *> done;
  for (CFaces::const_iterator f = faces.begin(); f != faces.end(); ++f)
    if (done.insert(*f).second) {
      CShell *s = new CShell;
      CFaces st;
      st.push_back(*f);
      bool flag = true;
      while (!st.empty()) {
	CFace *g = *st.rbegin();
	st.pop_back();
	s->faces.push_back(g);
	HCEdges ed;
	g->boundaryHEdges(ed);
	for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
	  if (!(*e)->e->piEdge()) {
	    CFace *h = g->neighbor(*e);
	    if (!h)
	      flag = false;
	    else if (done.insert(h).second)
	      st.push_back(h);
	  }
      }
      if (flag)
	sh.push_back(s);
      else
	delete s;
    }
}

CShell * Cspace::enclosingShell (CShell *s, Octree<CFace *> *octree)
{
  CPoint *a = s->vm->p;
  double bb[6];
  a->getBBox(bb);
  bb[0] = bbox[0];
  CFaces fa;
  octree->find(bb, fa);
  CFace *f = 0;
  PTR<CPoint> p = 0;
  for (CFaces::iterator g = fa.begin(); g != fa.end(); ++g)
    if ((*g)->s != s) {
      PTR<CPoint> q = new RayPatchPointX(a, *g);
      double qb[6];
      q->getBBox(qb);
      if (bboxOverlap(qb, (*g)->bbox) &&
	  CPointOrderX(q, a) == 1 && (*g)->contains2(q) &&
	  (!p || CPointOrderX(q, p) == 1)) {
	p = q;
	f = *g;
      }
    }
  return f ? f->s : 0;
}

void Cspace::removeFace (CFace *f) const
{
  HCEdges ed;
  f->boundaryHEdges(ed);
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    (*e)->e->removeHEdge(*e);
  delete f;
}

void Cspace::describe () const
{
  for (int i = 0; i < cells.size(); ++i) {
    CCell *c = cells[i];
    cerr << "cell " << i << ": " << c->shells.size() << " shells; ";
    for (CShells::iterator s = c->shells.begin(); s != c->shells.end(); ++s)
      cerr << (*s)->faces.size() << " hfaces ";
    cerr << endl;
  }
}

bool inIntervalZ (Vertex *v, Edge *e)
{
  Point *p = v->getP(), *t = e->getT()->getP(), *h = e->getH()->getP();
  if (p == t || p == h)
    return false;
  int s = PointOrderZ(t, h);
  return PointOrderZ(t, p) == s && PointOrderZ(p, h) == s;
}

Cspace * cspace (Polyhedron *a, Polyhedron *b)
{
  Cspace *c = new Cspace;
  c->patches(a, b);
  c->intersectFF();
  c->intersectFFF();
  c->sortVertices();
  Cspace *d = c->subfaces();
  d->formCells(a, b);
  return d;
}

Polyhedron * transform (Polyhedron *p, CPoint *c)
{
  Polyhedron *q = new Polyhedron;
  PVMap pvmap;
  for (Vertices::const_iterator v = p->vertices.begin(); v != p->vertices.end(); ++v)
    pvmap.insert(PVPair((*v)->getP(), q->getVertex(new TransformedPoint(c, (*v)->getP()))));
  for (Faces::const_iterator f = p->faces.begin(); f != p->faces.end(); ++f)
    q->addTriangle(*f, pvmap);
  return q;
}

Polyhedron * discretize (Cspace *c, double d)
{
  for (CFaces::iterator f = c->faces.begin(); f != c->faces.end(); ++f)
    (*f)->formSlabs();
  CPVMap pvmap;
  Polyhedron *a = discretize(c->faces, d, pvmap);
  facesPI(c, a, pvmap, true);
  facesPI(c, a, pvmap, false);
  return a;
}

Polyhedron * discretize (const CFaces &fa, double d, CPVMap &pvmap)
{
  EPMap epmap = discretizeEdges(fa, d);
  delentil(fa, epmap);
  STMap stm;
  for (CFaces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
    const Slabs &sl = (*f)->getSlabs();
    for (Slabs::const_iterator s = sl.begin(); s != sl.end(); ++s)
      stm.insert(STPair(*s, discretize(*s, epmap)));
  }
  Polyhedron *a = new Polyhedron;
  addTriangles(stm, pvmap, a);
  return a;
}

EPMap discretizeEdges (const CFaces &fa, double d)
{
  EPMap epmap;
  set<CEdge *> es;
  for (CFaces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
    HCEdges ed;
    (*f)->boundaryHEdges(ed);
    for (HCEdges::iterator h = ed.begin(); h != ed.end(); ++h) {
      CEdge *e = (*h)->getE();
      if (!e->piEdge())
	es.insert(e);
    }
  }
  for (set<CEdge *>::iterator e = es.begin(); e != es.end(); ++e)
    epmap.insert(EPPair(*e, discretize(*e, d)));
  return epmap;
}

CPoints discretize (CEdge *e, double d)
{
  CVertices ve;
  for (int i = 0; i < e->HEdgesN(); ++i) {
    CVertices vei;
    e->getHEdge(i)->getF()->boundaryVertices(vei);
    for (CVertices::iterator v = vei.begin(); v != vei.end(); ++v)
      if (e->contains((*v)->getP()->getA()))
	ve.push_back(*v);
  }
  ve.push_back(e->getT());
  ve.push_back(e->getH());
  sort(ve.begin(), ve.end(), VertexAngleOrder());
  int n = ve.size();
  CPoints pts;
  pts.push_back(ve[0]->getP());
  for (int i = 0; i + 1 < n; ++i) {
    Angle *a1 = ve[i]->getP()->getA(), *a2 = ve[i+1]->getP()->getA();
    if (a1 != a2) {
      discretize(e, a1, a2, d, pts);
      if (i + 2 < n)
	pts.push_back(e->point(a2));
    }
  }
  pts.push_back(ve[n-1]->getP());
  return pts;
}

void discretize (CEdge *e, Angle *as, Angle *ae, double d, CPoints &pts)
{
  PTR<CPoint> ps = e->point(as), pe = e->point(ae);
  CPoints st, res;
  res.push_back(ps);
  st.push_back(pe);
  while (!st.empty()) {
    CPoint *p = *res.rbegin(), *q = *st.rbegin();
    PTR<Angle> ba = new BisectorAngle(p->getA(), q->getA());
    double da = AngleParameter(ba).getApprox().x.l;
    PTR<LinearPoly> pa = new LinearPoly(da);
    PTR<Angle> a = new ParametricAngle(AngleUpper(ba) == 1, new LinearRoot(pa),
				       VertexSet());
    PTR<CPoint> m = e->point(a);
    if (!inInterval(a, p->getA(), q->getA()) || close(m, p, q, d)) {
      res.push_back(q);
      st.pop_back();
    }
    else
      st.push_back(m);
  }
  for (int i = 1; i + 1 < res.size(); ++i)
    pts.push_back(res[i]);
}

bool close (CPoint *a, CPoint *t, CPoint *h, double d)
{
  PointCPoint pa(a), pt(t), ph(h);
  return DistancePL(&pa, &pt, &ph, d) == 1;
}

void delentil (const CFaces &fa, EPMap &epmap)
{
  for (CFaces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
    const Slabs &sl = (*f)->getSlabs();
    for (Slabs::const_iterator s = sl.begin(); s != sl.end(); ++s)
      delentil(*s, epmap);
  }
}

void delentil (Slab *s, EPMap &epmap)
{
  CPoints l = getPoints(s->l, s->s, s->e, epmap),
    r = getPoints(s->r, s->s, s->e, epmap);
  if (l.size() == 2 && r.size() == 2 && l[0] == r[0] && l[1] == r[1]) {
    PTR<Angle> a = new BisectorAngle(s->s, s->e);
    addPoint(s->l, a, epmap);
  }
}

void addPoint (CEdge *e, Angle *a, EPMap &epmap)
{
  CPoints &pts = epmap.find(e)->second;
  for (CPoints::iterator i = pts.begin(); i + 1 != pts.end(); ++i)
    if (inInterval(a, (*i)->getA(), (*(i+1))->getA())) {
      pts.insert(i+1, e->point(a));
      break;
    }
}

CTriangles discretize (Slab *s, const EPMap &epmap)
{
  CPoints l = getPoints(s->l, s->s, s->e, epmap),
    r = getPoints(s->r, s->s, s->e, epmap);
  CTriangles trs = discretize(l, r, s), tr;
  int is = s->pb.empty() ? 0 : 1, n = trs.size(),
    ie = s->pt.empty() ? n : n - 1;
  if (is == 1)
    splitB(trs[0], s->pb, trs[0].a == l[1], tr);
  for (int i = is; i < ie; ++i)
    tr.push_back(trs[i]);
  if (ie < n)
    splitT(trs[n-1], s->pt, trs[n-1].a == *l.rbegin(), tr);
  return tr;
}

CPoints getPoints (CEdge *f, Angle *s, Angle *e, const EPMap &epmap)
{
  const CPoints &fpts = epmap.find(f)->second;
  CPoints pts;
  CPoints::const_iterator p = fpts.begin();
  while ((*p)->getA()->order(s))
    ++p;
  while (true) {
    pts.push_back(*p);
    if (!(*p)->getA()->order(e))
      break;
    ++p;
  }
  return pts;
}

CTriangles discretize (const CPoints &l, const CPoints &r, Slab *s)
{
  CTriangles tr;
  int il = 0, ir = 0, nl = l.size(), nr = r.size();
  while (ir + 1 < nr)
    if (r[ir]->getA()->order(l[il]->getA())) {
      int i2 = il + 1;
      while (i2 < nl && !r[ir+1]->getA()->order(l[i2]->getA()))
	++i2;
      double mid = 0.5*(r[ir]->getA()->theta() + r[ir+1]->getA()->theta());
      int m = il;
      for (int j = il + 1; j < i2; ++j)
        if (fabs(l[j]->getA()->theta() - mid) < fabs(l[m]->getA()->theta() - mid))
          m = j;
      for (int j = il + 1; j <= m; ++j)
	ctriangle(l[j], l[j-1], r[ir], tr);
      ctriangle(r[ir], r[ir+1], l[m], tr);
      for (int j = m + 1; j < i2; ++j)
	ctriangle(l[j], l[j-1], r[ir+1], tr);
      il = i2 - 1;
      ++ir;
    }
    else {
      int i2 = ir + 1;
      while (r[i2]->getA()->order(l[il+1]->getA()))
	++i2;
      double mid = 0.5*(l[il]->getA()->theta() + l[il+1]->getA()->theta());
      int m = ir;
      for (int j = ir + 1; j < i2; ++j)
        if (fabs(r[j]->getA()->theta() - mid) < fabs(r[m]->getA()->theta() - mid))
          m = j;
      for (int j = ir + 1; j <= m; ++j)
	ctriangle(r[j-1], r[j], l[il], tr);
      ctriangle(l[il+1], l[il], r[m], tr);
      for (int j = m + 1; j < i2; ++j)
	ctriangle(r[j-1], r[j], l[il+1], tr);
      ir = i2 - 1;
      ++il;
    }
  return tr;
}

void ctriangle (CPoint *a, CPoint *b, CPoint *c, CTriangles &tr)
{
  if (a == b || a == c || b == c)
    return;
  tr.push_back(CTriangle(a, b, c));
}

void splitB (const CTriangle &t, const CPoints &pb, bool lflag, CTriangles &tr)
{
  int n = pb.size();
  if (lflag) {
    ctriangle(t.a, t.b, pb[0], tr);
    for (int i = 0; i < n; ++i)
      ctriangle(t.a, pb[i], i + 1 < n ? pb[i+1] : t.c, tr);
  }
  else {
    ctriangle(t.a, t.b, pb[n-1], tr);
    for (int i = 0; i < n; ++i)
      ctriangle(t.b, i == 0 ? t.c : pb[i-1], pb[i], tr);
  }
}

void splitT (const CTriangle &t, const CPoints &pt, bool lflag, CTriangles &tr)
{
  int n = pt.size();
  if (lflag) {
    ctriangle(t.a, t.b, pt[0], tr);
    for (int i = 0; i < n; ++i)
      ctriangle(t.b, i + 1 < n ? pt[i+1] : t.c, pt[i], tr);
  }
  else {
    ctriangle(t.a, t.b, pt[n-1], tr);
    for (int i = 0; i < n; ++i)
      ctriangle(t.a, pt[i], i == 0 ? t.c : pt[i-1], tr);
  }
}

void addTriangles (const STMap &stm, CPVMap &pvmap, Polyhedron *a)
{
  for (STMap::const_iterator i = stm.begin(); i != stm.end(); ++i)
    for (CTriangles::const_iterator t = i->second.begin(); t != i->second.end(); ++t)
      a->addTriangle(getVertex(t->a, pvmap, a),
		     getVertex(t->b, pvmap, a),
		     getVertex(t->c, pvmap, a));
}

Vertex * getVertex (CPoint *p, CPVMap &pvmap, Polyhedron *a)
{
  CPVMap::iterator i = pvmap.find(p);
  if (i != pvmap.end())
    return i->second;
  Vertex *v = a->getVertex(new PointCPoint(p));
  pvmap.insert(CPVPair(p, v));
  return v;
}

void facesPI (Cspace *c, Polyhedron *a, CPVMap &cpvmap, bool flag)
{
  CEdgeSet es;
  PVMap pvmap;
  vector<Points *> reg;
  for (CEdges::const_iterator e = c->edges.begin(); e != c->edges.end(); ++e)
    if ((*e)->HEdgesN() == 1 && (*e)->piEdge() &&
	(*e)->getT()->getP()->getA() == (flag ? Angle::ppi : Angle::mpi) &&
	es.find(*e) == es.end())
      reg.push_back(loopPI(a, cpvmap, es, *e, pvmap));
  Triangles tr;
  triangulate(reg, flag ? 3 : -3, tr);
  for (Triangles::iterator t = tr.begin(); t != tr.end(); ++t)
    a->addTriangle(pvmap.find(t->a)->second, pvmap.find(t->b)->second,
		   pvmap.find(t->c)->second);
  for (vector<Points *>::const_iterator r = reg.begin(); r != reg.end(); ++r)
    delete(*r);
}

Points * loopPI (Polyhedron *a, CPVMap &cpvmap, CEdgeSet &es, CEdge *e,
		 PVMap &pvmap)
{
  Points *pts = new Points;
  CEdge *e0 = e;
  do {
    es.insert(e);
    Vertex *v = getVertex(e->getH()->getP(), cpvmap, a);
    pvmap.insert(PVPair(v->getP(), v));
    pts->push_back(v->getP());
    CVertex *t = e->getT();
    for (int i = 0; i < t->EdgesN(); ++i) {
      CEdge *f = t->getEdge(i);
      if (e != f && f->HEdgesN() == 1) {
	e = f;
	break;
      }
    }
  }
  while (e != e0);
  return pts;
}

// debug

PointCPoint * pcpt (Point *p)
{
  return dynamic_cast<PointCPoint *>(p);
}

SpiralPoint * spt (CPoint *p)
{
  return dynamic_cast<SpiralPoint *>(p);
}

LinePatchPoint * lppt (CPoint *p)
{
  return dynamic_cast<LinePatchPoint *>(p);
}

ParametricAngle * pangle (Angle *a)
{
  return dynamic_cast<ParametricAngle *>(a);
}

FFPoint * ffp (CPoint *p)
{
  return dynamic_cast<FFPoint *>(p);
}

void pa1 (Angle *a)
{
  if (a == Angle::mpi)
    cerr << "-pi";
  else if (a == Angle::ppi)
    cerr << "pi";
  else
    cerr << a->theta();
}

void pa (Angle *a)
{
  pa1(a);
  cerr << endl;
}

void pai1 (AngleInterval &ai)
{
  cerr << "(";
  pa1(ai.first);
  cerr << " ";
  pa1(ai.second);
  cerr << ")";
}

void pai (AngleInterval &ai)
{
  pai1(ai);
  cerr << endl;
}

void pais (AngleIntervals &ais)
{
  cerr << "(";
  for (AngleIntervals::iterator i = ais.begin(); i != ais.end(); ++i)
    pai1(*i);
  cerr << ")" << endl;
}

void pp1 (CPoint *p)
{
  PV2<double> xy = p->getApproxMid();
  cerr << "(" << xy.x << " " << xy.y << " ";
  pa1(p->getA());
  cerr << ")";
}

void pp (CPoint *p)
{
  cerr << setprecision(16);
  pp1(p);
  cerr << endl;
}

void pv (CVertex *v)
{
  pp(v->getP());
}

void pvs (const CVertices &ve)
{
  cerr << "(" << endl;
  for (CVertices::const_iterator v = ve.begin(); v != ve.end(); ++v)
    pv(*v);
  cerr << ")" << endl;
}

void pe (CEdge *e)
{
  cerr << "(";
  pp1(e->getT()->getP());
  cerr << " ";
  pp1(e->getH()->getP());
  cerr << ")" << endl;
}

void pes (const CEdges &ed)
{
  cerr << "(";
  for (CEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e);
  cerr << ")" << endl;
}

void pe (HCEdge *e)
{
  cerr << "(";
  pp1(e->tail()->getP());
  cerr << " ";
  pp1(e->head()->getP());
  cerr << ")" << endl;
}

void pes (const HCEdges &ed)
{
  cerr << "(";
  for (HCEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e);
  cerr << ")" << endl;
}

typedef vector<PV3<double>> PV3s;

void pedge (CEdge *e, double d, PV3s &pts)
{
  CPoints cpts;
  cpts.push_back(e->getT()->getP());
  if (!e->horizontal())
    discretize(e, e->getT()->getP()->getA(), e->getH()->getP()->getA(), d, cpts);
  cpts.push_back(e->getH()->getP());
  for (CPoints::iterator i = cpts.begin(); i != cpts.end(); ++i) {
    PointCPoint p(*i);
    pts.push_back(p.getApproxMid());
  }
}

void pp1 (const PV3<double> &p);

void ppts (const PV3s &pts)
{
  cerr << "(";
  for (int i = 0; i < pts.size(); ++i)
    pp1(pts[i]);
  cerr << ")" << endl;
}

void pe (CEdge *e, double d)
{
  PV3s pts;
  pedge(e, d, pts);
  ppts(pts);
}

void pes (const CEdges &ed, double d)
{
  cerr << "(";
  for (CEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e, d);
  cerr << ")" << endl;
}

void pe (CEdge *e, Angle *as, Angle *ae, double d)
{
  CPoints cpts;
  discretize(e, as, ae, d, cpts);
  PV3s pts;
  for (CPoints::iterator i = cpts.begin(); i != cpts.end(); ++i) {
    PointCPoint p(*i);
    pts.push_back(p.getApproxMid());
  }
  ppts(pts);
}

void pe (HCEdge *e, double d)
{
  PV3s pts;
  pedge(e->getE(), d, pts);
  if (!e->getForward())
    reverse(pts.begin(), pts.end());
  ppts(pts);
}

void pes (const HCEdges &ed, double d)
{
  cerr << "(";
  for (HCEdges::const_iterator e = ed.begin(); e != ed.end(); ++e)
    pe(*e, d);
  cerr << ")" << endl;
}

void pl (HCEdge *e)
{
  HCEdges ed;
  e->loop(ed);
  CVertices ve;
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    ve.push_back((*e)->head());
  pvs(ve);
}

void pl (HCEdge *e, double d)
{
  HCEdges ed;
  e->loop(ed);
  cerr << "(";
  for (HCEdges::iterator f = ed.begin(); f != ed.end(); ++f)
    pe(*f, d);
  cerr << ")" << endl;
}

void pf (CFace *f, Angle *a)
{
  PTR<CPoint> p = f->leftEdge()->point(a), q = f->rightEdge()->point(a);
  cerr << "(";
  pp1(p);
  cerr << " ";
  pp1(q);
  cerr << ")" << endl;
}

bool find (Spiral *s, const SpiralSet &ss)
{
  return ss.find(s) != ss.end();
}

int find (CPoint *p, const CPoints &pts)
{
  for (int i = 0; i < pts.size(); ++i)
    if (p == pts[i])
      return i;
  return -1;
}

int find (CVertex *v, const CVertices &vertices)
{
  for (int i = 0; i < vertices.size(); ++i)
    if (v == vertices[i])
      return i;
  return -1;
}

int find (CPoint *p, const CVertices &vertices)
{
  for (int i = 0; i < vertices.size(); ++i)
    if (p == vertices[i]->getP())
      return i;
  return -1;
}

int find (CEdge *e, const CEdges &edges)
{
  for (int i = 0; i < edges.size(); ++i)
    if (e == edges[i])
      return i;
  return -1;
}

int find (CFace *f, const CFaces &faces)
{
  for (int i = 0; i < faces.size(); ++i)
    if (f == faces[i])
      return i;
  return -1;
}

int find (Slab *s, const Slabs &sl)
{
  for (int i = 0; i < sl.size(); ++i)
    if (s == sl[i])
      return i;
  return -1;
}

void plines (const vector<PV3s> &lines, int i);;

void pedge (CEdge *e, double d, int i)
{
  PV3s pts;
  pedge(e, d, pts);
  vector<PV3s> lines;
  lines.push_back(pts);
  plines(lines, i);
}

void pedges (const CEdges &edges, double d, int i)
{
  vector<PV3s> lines;
  for (CEdges::const_iterator e = edges.begin(); e != edges.end(); ++e) {
    PV3s pts;
    pedge(*e, d, pts);
    lines.push_back(pts);
  }
  plines(lines, i);
}

void pedges (const HCEdges &edges, double d, int i)
{
  vector<PV3s> lines;
  for (HCEdges::const_iterator e = edges.begin(); e != edges.end(); ++e) {
    PV3s pts;
    pedge((*e)->getE(), d, pts);
    lines.push_back(pts);
  }
  plines(lines, i);
}

void pboundary (CFace *f, double d, int i)
{
  HCEdges ed;
  f->boundaryHEdges(ed);
  pedges(ed, d, i);
}

void pfaces (const Faces &fa, int i);

void ptriangles (const CTriangles &tr, int i)
{
  CPVMap pvmap;
  Polyhedron *a = new Polyhedron;
  for (CTriangles::const_iterator t = tr.begin(); t != tr.end(); ++t)
    a->addTriangle(getVertex(t->a, pvmap, a),
		   getVertex(t->b, pvmap, a),
		   getVertex(t->c, pvmap, a));
  pfaces(a->faces, i);
  delete a;
}

void ptriangle (const CTriangle &t, int i)
{
  CTriangles tr;
  tr.push_back(t);
  ptriangles(tr, i);
}

void pfaces (const CFaces &fa, double d, int i)
{
  for (CFaces::const_iterator f = fa.begin(); f != fa.end(); ++f)
    (*f)->formSlabs();
  EPMap epmap = discretizeEdges(fa, d);
  delentil(fa, epmap);
  STMap stm;
  for (CFaces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
    const Slabs &sl = (*f)->getSlabs();
    for (Slabs::const_iterator s = sl.begin(); s != sl.end(); ++s)
      stm.insert(STPair(*s, discretize(*s, epmap)));
  }
  CTriangles tr;
  for (STMap::iterator j = stm.begin(); j != stm.end(); ++j)
    tr.insert(tr.end(), j->second.begin(), j->second.end());
  ptriangles(tr, i);
}

void pface (CFace *f, double d, int i)
{
  CFaces fa;
  fa.push_back(f);
  pfaces(fa, d, i);
}

void pfaces (const CFaces &fa, double d, int i, int is, int ie)
{
  CFaces fad;
  for (int j = is; j < ie; ++j)
    fad.push_back(fa[j]);
  pfaces(fad, d, i);
}

void pslab (Slab *s, double d, int i)
{
  CFaces fa;
  fa.push_back(s->f);
  EPMap epmap = discretizeEdges(fa, d);
  CTriangles tr = discretize(s, epmap);
  ptriangles(tr, i);
}

double distance (CPoint *a, CPoint *b)
{
  PV2<double> u = a->getApproxMid() - b->getApproxMid();
  double uu = u.dot(u),
    da = fabs(a->getA()->theta() - b->getA()->theta());
  return uu + da;
}

double distance (CPoint *p, CFace *f)
{
  Angle *pa = p->getA();
  Spiral *l = f->leftEdge()->getSpiral(), *r = f->rightEdge()->getSpiral();
  PV2<double> a = p->getApproxMid(),
    t = l->xyApprox(pa), h = r->xyApprox(pa),
    u = (h - t).unit(), v(a.x - t.x, a.y - t.y);
  return u.cross(v);
}

void findSpiral (Cspace *a, Spiral *s)
{
  for (int i = 0; i < a->faces.size(); ++i)
    if (a->faces[i]->leftEdge()->getSpiral() == s ||
	a->faces[i]->rightEdge()->getSpiral() == s)
      cerr << i << " ";
  cerr << endl;
}

void findAngle (Cspace *a, Angle *x)
{
  for (int i = 0; i < a->faces.size(); ++i)
    if (a->faces[i]->startAngle() == x || a->faces[i]->endAngle() == x)
      cerr << i << " ";
  cerr << endl;
}

CPoints crossSection (CFace *f, Angle *a)
{
  HCEdges ed;
  f->boundaryHEdges(ed);
  CPoints pts;
  for (HCEdges::iterator e = ed.begin(); e != ed.end(); ++e)
    if ((*e)->getE()->contains(a))
      pts.push_back((*e)->getE()->point(a));
  sort(pts.begin(), pts.end(), CPointOrderXO());
  return pts;
}

class CIntersectsP : public Primitive {
  CPoint *a, *b, *c, *d;

  DeclareSign {
    PV2<N> pa = a->get<N>(), pb = b->get<N>(), pc = c->get<N>(),
      pd = d->get<N>(), u = pb - pa;
    N k1 = u.cross(pc - pa), k2 = u.cross(pd - pa);
    if (k1.sign()*k2.sign() != -1)
      return N(-1);
    PV2<N> v = pd - pc;
    N k3 = v.cross(pa - pc), k4 = v.cross(pb - pc);
    return - k3*k4;
  }
 public:
  CIntersectsP (CPoint *a, CPoint *b, CPoint *c, CPoint *d)
    : a(a), b(b), c(c), d(d) {}
};

bool intersects (CFace *f, CFace *g, Angle *a)
{
  CPoints fpts = crossSection(f, a), gpts = crossSection(g, a);
  bool res = false;
  for (int i = 0; !res && i < fpts.size(); i += 2)
    for (int j = 0; !res && j < gpts.size(); j += 2)
      res = CIntersectsP(fpts[i], fpts[i+1], gpts[i], gpts[i+1]) == 1;
  return res;
}

