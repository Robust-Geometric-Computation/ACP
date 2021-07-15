#include "triangulate.h"

void Vertex2::addEdge (Edge2 *e)
{
  if (!edge) {
    edge = e;
    e->next = e;
  }
  else if (e->clockwise(edge)) {
    e->next = edge;
    Edge2 *f = edge;
    while (f->next != edge)
      f = f->next;
    f->next = e;
    edge = e;
  }
  else {
    Edge2 *f = edge;
    while (f->next != edge && f->next->clockwise(e))
      f = f->next;
    e->next = f->next;
    f->next = e;
  }
}

Edge2 * Vertex2::findEdge (Vertex2 *h) const
{
  if (!edge)
    return 0;
  Edge2 *e = edge;
  do {
    if (e->head() == h)
      return e;
    e = e->next;
  }
  while (e != edge);
  return 0;
}

int Vertex2::yOrder (Vertex2 *b) const
{
  return PointOrderR(p, b->p);
}

int Vertex2::leftTurn (Vertex2 *b, Vertex2 *c) const
{
  return LeftTurn(p, b->p, c->p, coord);
}

bool Edge2::increasingY ()
{
  return tail->yOrder(head()) == 1;
}

bool Edge2::clockwise (Edge2 *e)
{
  bool inc = increasingY(), einc = e->increasingY();
  return inc != einc ? inc : head()->leftTurn(tail, e->head()) == 1;
}

Edge2 * getEdge (Vertex2 *tail, Vertex2 *head, bool flag, bool tflag)
{
  Edge2 *oe = tail->findEdge(head);
  if (oe)
    return oe;
  Edge2 *e = new Edge2(tail, 0, flag), *et = new Edge2(head, e, tflag);
  e->twin = et;
  tail->addEdge(e);
  head->addEdge(et);
  return e;
}

void triangulate (const vector<Points *> &reg, int coord, Triangles &tr)
{
  Vertices2 vertices;
  triangulateInit(reg, coord, vertices);
  triangulate(vertices, tr);
  deleteVertices(vertices);
}

void triangulateInit (const vector<Points *> &reg, int coord, Vertices2 &vertices)
{
  PV2Map vmap;
  for (int i = 0; i < reg.size(); ++i) {
    const Points &b = *reg[i];
    int n = b.size();
    for (int ie = 0; ie < n; ++ie) {
      int ien = ie + 1 == n ? 0 : ie + 1;
      Vertex2 *t = getVertex2(b[ie], coord, vmap),
	*h = getVertex2(b[ien], coord, vmap);
      Edge2 *e = getEdge(t, h, false, true);
      e->flag = false;
    }
  }
  for (PV2Map::iterator x = vmap.begin(); x != vmap.end(); ++x)
    vertices.push_back(x->second);
}

Vertex2 * getVertex2 (Point *p, int coord, PV2Map &vmap)
{
  PV2Map::iterator x = vmap.find(p);
  if (x != vmap.end())
    return x->second;
  Vertex2 *w = new Vertex2(p, coord);
  vmap.insert(pair<Point *, Vertex2 *>(p, w));
  return w;
}

void triangulate (Vertices2 &vertices, Triangles &tr)
{
  makeMonotone(vertices);
  triangulateMonotone(vertices, tr);
}

void makeMonotone (Vertices2 &vertices)
{
  sort(vertices.begin(), vertices.end(), Vertex2yOrder());
  Edge2Set sweep;
  for (Vertices2::iterator v = vertices.begin(); v != vertices.end(); ++v)
    makeMonotone(*v, sweep);
}

void makeMonotone (Vertex2 *v, Edge2Set &sweep)
{
  Vertices2 d;
  bool dflag = true;
  Edge2 *e = v->edge;
  while (e->increasingY()) {
    dflag = false;
    if (e->helper->flag)
      d.push_back(e->helper);
    sweep.erase(e);
    e = e->next;
    if (e == v->edge)
      break;
  }
  Edge2 *f = findLeft(v, sweep);
  if (f) {
    if (f->helper->flag || dflag)
      d.push_back(f->helper);
    f->helper = v;
  }
  v->flag = true;
  while (!e->increasingY()) {
    v->flag = false;
    e->twin->helper = v;
    sweep.insert(e->twin);
    e = e->next;
    if (e == v->edge)
      break;
  }
  for (Vertices2::iterator w = d.begin(); w != d.end(); ++w)
    getEdge(v, *w, false, false);
}

Edge2 * findLeft (Vertex2 *v, Edge2Set &sweep)
{
  if (sweep.empty())
    return 0;
  Edge2 *e = v->edge->increasingY() ? v->edge : v->edge->twin;
  Edge2Set::iterator f = sweep.lower_bound(e);
  if (f == sweep.begin())
    return 0;
  --f;
  return (*f)->twin->flag ? 0 : *f;
}

void triangulateMonotone (const Vertices2 &vertices, Triangles &tr)
{
  for (Vertices2::const_iterator v = vertices.begin(); v != vertices.end(); ++v) {
    Edge2 *e0 = (*v)->edge, *e = e0;
    do {
      if (!e->flag && !e->increasingY()) {
	Vertices2 u;
	monotoneLoop(e, u);
	triangulateLoop(u, tr);
      }
      e = e->next;
    }
    while (e != e0);
  }
}

void monotoneLoop (Edge2 *e, Vertices2 &u)
{
  bool rflag = false;
  while (!e->flag) {
    rflag = rflag || e->increasingY();
    e->tail->flag = rflag;
    u.push_back(e->tail);
    e->flag = true;
    e = e->twin->next;
  }
  sort(u.begin(), u.end(), Vertex2yOrder());
}

void triangulateLoop (const Vertices2 &u, Triangles &tr)
{
  Vertices2 stack;
  stack.push_back(u[0]);
  stack.push_back(u[1]);
  int n = u.size();
  for (int i = 2; i + 1 < n; ++i)
    if (u[i]->flag == stack.back()->flag)
      triangulateSame(u[i], stack, tr);
    else
      triangulateOther(u[i], stack, tr);
  int m = stack.size() - 1;
  if (u[n-1]->p->onLine(stack[m]->p, stack[m-1]->p)) {
    stack.push_back(u[n-1]);
    for (int i = 1; i <= m; ++i)
      tr.push_back(triangle(stack[0], stack[i], stack[i+1], stack[0]->flag));
  }
  else
    triangulateOther(u[n-1], stack, tr);
}

void triangulateSame (Vertex2 *v, Vertices2 &stack, Triangles &tr)
{
  int n = stack.size() - 1;
  while (n > 0) {
    if (v->p->onLine(stack[n-1]->p, stack[n]->p))
      break;
    bool aflag = stack[n-1]->leftTurn(stack[n], v) == 1;
    if (v->flag ? aflag : !aflag)
      break;
    tr.push_back(triangle(v, stack[n], stack[n-1], v->flag));
    stack.pop_back();
    --n;
  }
  stack.push_back(v);
}

void triangulateOther (Vertex2 *v, Vertices2 &stack, Triangles &tr)
{
  int n = stack.size() - 1;
  Vertex2 *q = stack[n];
  bool flag = q->flag;
  while (n > 0) {
    tr.push_back(triangle(v, stack[n], stack[n-1], flag));
    --n;
  }
  stack.clear();
  stack.push_back(q);
  stack.push_back(v);
}

Triangle triangle (Vertex2 *a, Vertex2 *b, Vertex2 *c, bool flag)
{
  return flag ? Triangle(a->p, b->p, c->p) : Triangle(a->p, c->p, b->p);
}

void deleteVertices (Vertices2 &vertices)
{
  for (unsigned int i = 0; i < vertices.size(); ++i) {
    vector<Edge2 *> edges;
    Edge2 *e = vertices[i]->edge, *e0 = e;
    do {
      edges.push_back(e);
      e = e->next;
    }
    while (e != e0);
    for (unsigned int j = 0; j < edges.size(); ++j)
      delete edges[j];
    delete vertices[i];
  }
}
