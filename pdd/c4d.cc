#include "c4d.h"
#include <fstream>
using namespace std;

#ifndef NOTTIME
int objCount;
#endif

bool noSweep = false;

C4D::Event *currentEvent;

#ifdef BLEEN
void C4D::EventSet::approximateUpdate (Event *event, double crossStep) {
  if (isAll()) return;

  for (int i = 0; i < events.size(); i += 2) {
    PTR<Event> a = events[i];
    PTR<Event> b = events[i+1];
    if (a == b) {
      double cross = event->getU().cross(a->getU()).mid();
      if (cross > 0 && cross < crossStep) {
        events[i] = event;
        return;
      }
      if (cross < 0 && cross > -crossStep) {
        events[i+1] = event;
        return;
      }
      continue;
    }

    if (event == a || event == b || event->between(a, b) > 0)
      return;
    if (fabs(event->getU().cross(a->getU()).mid()) < crossStep) {
      events[i] = event;
      return;
    }
    if (fabs(event->getU().cross(b->getU()).mid()) < crossStep) {
      events[i+1] = event;
      return;
    }
  }    

  events.push_back(event);
  events.push_back(event);
}
#endif

#ifdef BLEEN
int C4D::Circulation::sign () {
  if (a->getU().cross(c->getU()).sign() > 0) {
    if (a->getU().cross(b->getU()).sign() > 0 &&
        c->getU().cross(b->getU()).sign() < 0)
      return 1;
  }
  else {
    if (a->getU().cross(b->getU()).sign() > 0 ||
        c->getU().cross(b->getU()).sign() < 0)
      return 1;
  }
  return -1;
}
#endif

// End angles are sorted on [0, 2 PI).
void C4D::EventSet::sort (vector<PTR<Event> > &events) {
  for (int n = 2; n < events.size(); n += 2) {
    PTR<Event> a = events[n];
    PTR<Event> b = events[n+1];
    int i = n;
    while (i > 0 && b->lessThan(events[i-1])) {
      events[i] = events[i-2];
      events[i+1] = events[i-1];
      i -= 2;
    }
    events[i] = a;
    events[i+1] = b;
  }
}

void check (const vector<PTR<C4D::Event>> &a) {
  assert(a.size() % 2 == 0);
  if (a.size() <= 2) return;

  for (int i = 2; i < a.size(); i++)
    assert(Angle::compare(a[i-1]->getAngle(), a[i]->getAngle()) < 0);

  assert(Angle::compare(a[0]->getAngle(), a[1]->getAngle()) < 0 ||
	 Angle::compare(a[0]->getAngle(), a.back()->getAngle()) > 0);
}

C4D::EventSet C4D::EventSet::intersect (const EventSet &that) const {
  if (isNull()) return *this;
  if (that.isNull()) return that;
  if (isAll()) return that;
  if (that.isAll()) return *this;

  EventSet out;
  const vector<PTR<Event> > &a = events;
  const vector<PTR<Event> > &b = that.events;
  vector<PTR<Event> > &c = out.events;
  
  static int count;
  if (++count == 0)
    cout << "this is it" << endl;

  check(a);
  check(b);
    
  for (int i = 0; i < a.size(); i++)
    assert(!a[i]->getAngle()->identical(a[(i+1)%a.size()]->getAngle())); 
  for (int i = 0; i < b.size(); i++)
    assert(!b[i]->getAngle()->identical(b[(i+1)%b.size()]->getAngle())); 

  for (int i = 0; i < a.size(); i += 2) {
    for (int j = 0; j < b.size(); j += 2) {
      bool x2 = a[i]->identical(b[j]);

      if (a[i]->between(b[j], b[j+1]) || a[i]->identical(b[j])) {
        c.push_back(a[i]);
        if (b[j+1]->between(a[i], a[i+1]))
          c.push_back(b[j+1]);
        else {
          // assert(a[i+1]->between(b[j], b[j+1]));
          c.push_back(a[i+1]);
        }
	assert(!c[c.size()-2]->getAngle()->identical(c[c.size()-1]->getAngle()));
      }
      if (b[j]->between(a[i], a[i+1])) {
        c.push_back(b[j]);
        if (a[i+1]->between(b[j], b[j+1]))
          c.push_back(a[i+1]);
        else {
          // assert(b[j+1]->between(a[i], a[i+1]));
          c.push_back(b[j+1]);
        }
	assert(!c[c.size()-2]->getAngle()->identical(c[c.size()-1]->getAngle()));
      }

      int cmpaa = Angle::compare(a[i+0]->getAngle(), a[i+1]->getAngle());
      int cmpbb = Angle::compare(b[j+0]->getAngle(), b[j+1]->getAngle());
      int cmp00 = Angle::compare(a[i+0]->getAngle(), b[j+0]->getAngle());
      int cmp01 = Angle::compare(a[i+0]->getAngle(), b[j+1]->getAngle());
      int cmp10 = Angle::compare(a[i+1]->getAngle(), b[j+0]->getAngle());
      int cmp11 = Angle::compare(a[i+1]->getAngle(), b[j+1]->getAngle());
      if (count == -2)
	cout << count << endl;
    }
  }
  for (int i = 0; i < c.size(); i++)
    assert(!c[i]->getAngle()->identical(c[(i+1)%c.size()]->getAngle())); 


  sort(c);
  check(c);
  for (int i = 0; i < c.size(); i++)
    assert(!c[i]->getAngle()->identical(c[(i+1)%c.size()]->getAngle())); 
  return out;
}

C4D::EventSet C4D::EventSet::complement () const 
{
  if (isNull()) return EventSet(true);
  if (isAll()) return EventSet();
  EventSet out;
  if (Angle::compare(events[0]->getAngle(), events[1]->getAngle()) < 0) {
    out.events.push_back(events.back());
    for (int i = 0; i < events.size()-1; i++)
      out.events.push_back(events[i]);
  }
  else {
    for (int i = 1; i < events.size(); i++)
      out.events.push_back(events[i]);
    out.events.push_back(events[0]);
  }
  return out;
}

bool C4D::EventSet::contains (Angle *angle) const {
  if (isNull()) return false;
  if (isAll()) return true;
  for (int i = 0; i < events.size(); i += 2)
    if (Angle::circulation(events[i]->getAngle(), angle, events[i+1]->getAngle()) > 0)
      return true;
  return false;
}

bool C4D::EventSet::contains (Event *event) const {
  if (isNull()) return false;
  if (isAll()) return true;
  for (int i = 0; i < events.size(); i += 2)
    if (event->between(events[i], events[i+1]))
      return true;
  return false;
}

void C4D::Obj::informEvents() {
  vector<PTR<Event> > &events = life.events;
  set<Event*, Compare> &eventList = getC4D()->eventList;
  for (int i = 0; i < events.size(); i++) {
    PTR<Event> event = events[i];
    set<Event*,Compare>::iterator it = eventList.find(event);
    if (it == eventList.end())
      eventList.insert(event);
    else if (event != *it)
      // eliminate duplicate event
      // space leak, use smart pointer for event
      events[i] = event = *it;
    if (i%2 == 0)
      event->appear.insert(this);
    else
      event->vanish.insert(this);
  }
}

int C4D::Edge::insertIndex (Event *event, Vert *v) {
  if (v == (Vert *) 0xa079b0 || v == (Vert *) 0xa06da0)
    cout << "insertIndex this is it" << endl;
  int a = 0, b = verts.size()-1;
  while (a <= b) {
    int c = (a+b)/2;
    if (OrderUVE(event, v, dynamic_cast<Vert*>(verts[c]), this) < 0)
      b = c - 1;
    else
      a = c + 1;
  }
  return a;
}

// Add Vert to verts.
void C4D::Edge::addVert (Event *event, Vert *v) {
  if (!smallestTwin ()) {
    twin->addVert(event, v);
    return;
  }        
  for (int i = 0; i < verts.size(); i++)
    if (verts[i] == v)
      return;
  int a = insertIndex(event, v);
  verts.insert(verts.begin() + a, v);
}

int C4D::Edge::findVert (Vert *v) {
#ifdef BLEEN1
  if (!smallestTwin())
    return twin->findVert(v);
#endif
  for (int i = 0; i < verts.size(); i++)
    if (verts[i] == v)
      return i;
  return -1;
}

vector<C4D::SubEdge*> C4D::Edge::makeSubEdges (set<Obj*, Compare> &objs) {
  vector<SubEdge*> subs;
  if (!smallestTwin()) {
    assert(0);
    return subs;
  }

  Vert *t = 0, *h = getTail();
  for (int j = 0; j <= verts.size(); j++) {
    t = h;
    h = j < verts.size() ? dynamic_cast<Vert*>(verts[j]) : getHead();
    SubEdge *sub = SubEdge::make(t, h, this);
    objs.insert(sub);
    objs.insert(sub->twin);
    subs.push_back(sub);
  }
  return subs;
}

vector<C4D::SubEdge*> C4D::Edge::makeSubEdges (unordered_set<Obj*, Equal, Equal> &objs) {
  vector<SubEdge*> subs;
  if (!smallestTwin()) {
    assert(0);
    return subs;
  }

  Vert *t = 0, *h = getTail();
  for (int j = 0; j <= verts.size(); j++) {
    t = h;
    h = j < verts.size() ? dynamic_cast<Vert*>(verts[j]) : getHead();
    SubEdge *sub = SubEdge::make(t, h, this);
    objs.insert(sub);
    objs.insert(sub->twin);
    subs.push_back(sub);
  }
  return subs;
}

vector<C4D::SubEdge*> C4D::Edge::makeSubEdges (unordered_set<Obj*, Equal, Equal> *objs) {
  vector<SubEdge*> subs;
  if (!smallestTwin()) {
    assert(0);
    return subs;
  }

  Vert *t = 0, *h = getTail();
  for (int j = 0; j <= verts.size(); j++) {
    t = h;
    h = j < verts.size() ? dynamic_cast<Vert*>(verts[j]) : getHead();
    SubEdge *sub = SubEdge::make(t, h, this);
    objs[sub->getType()].insert(sub);
    objs[sub->twin->getType()].insert(sub->twin);
    subs.push_back(sub);
  }
  return subs;
}

template<class N>
N C4D::CircFFF::calculate () {
  return a->get<N>(event).cross(b->get<N>(event)).dot(c->get<N>(event));
}

template<class N>
N C4D::LeftVEVE::calculate () {
  return ((v1->get<N>(event) - v2->get<N>(event)).dot
          (e1->get<N>(event).cross(e2->get<N>(event))));    
}

C4D::EdgeAB *C4D::VertXYZ::getEdge (Vert *head) {
  for (int i = 0; i < edges.size(); i++) {
    if (edges[i]->getHead() == head)
      return dynamic_cast<EdgeAB*>(edges[i]);
  }
  return 0;
}

C4D::FaceABC *C4D::VertXYZ::getFace (Vert *b, Vert *c) {
  EdgeAB *ab = getEdge(b);
  if (ab == 0)
    return 0;
  if (ab->face->contains(c))
    return ab->getFace();
  if (ab->twin->face->contains(c))
    return ab->getTwin()->getFace();
  return 0;
}

void C4D::VertXYZ::addEdge (Edge *edge) {
  assert(edge->getTail() == this);
  if (getEdge(edge->getHead()) == 0)
    edges.push_back(edge);
}

void C4D::VertXYZ::makeCone () {
  assert(edges.size() >= 3);
  vector<Edge *> e, extra;
  e.push_back(edges[0]);
  if (CircEEE0(edges[0], edges[1], edges[2]) < 0) {
    e.push_back(edges[1]);
    e.push_back(edges[2]);
  } else {
    e.push_back(edges[2]);
    e.push_back(edges[1]);
  }

  for (int j = 3; j < edges.size(); j++) {
    Edge *edge = edges[j];
    int lastCirc = CircEEE0(e.back(), e[0], edge);
    int prevCirc = lastCirc;
    int iStart = -1;
    int iEnd = -1;
    for (int i = 0; i < e.size(); i++) {
      int currCirc = i == e.size()-1 ? lastCirc : CircEEE0(e[i], e[i+1], edge);
      if (prevCirc > 0 && currCirc < 0)
        iStart = i;
      else if (prevCirc < 0 && currCirc > 0)
        iEnd = i;
      prevCirc = currCirc;
    }
    if (iStart == -1) {
      if (lastCirc < 0)
        continue;
      else
        return;
    }

    extra.clear();
    if (iEnd < iStart) {
      for (int i = iStart; i < e.size(); i++)
        extra.push_back(e[i]);
      for (int i = 0; i <= iEnd; i++)
        extra.push_back(e[i]);
    }
    else
      for (int i = iStart; i <= iEnd; i++)
        extra.push_back(e[i]);
    e.clear();
    e = extra;
    e.push_back(edge);
  }

  Edge *first = e[0];
  edges.clear();
  Edge *curr = first;
  do {
    edges.push_back(curr);
    const vector<Edge *> &fedges = curr->face->getOuterLoop();
    for (int i = 0; i < fedges.size(); i++)
      if (fedges[i]->getHead() == this)
        curr = fedges[i]->getTwin();
  } while (curr != first);

  int n = 0;
  int inds[3];
  for (int i = 0; n < 3 && i < edges.size(); i++) {
    for (int j = 0; j < e.size(); j++)
      if (edges[i] == e[j]) {
        inds[n++] = j;
        break;
      }
  }
  // This may fail.  What you need to do is construct a cone for a
  // single face-connected component and check it has the right
  // orientation.
  assert(n == 3);
  if ((inds[1] - inds[0]) * (inds[2] - inds[1]) * (inds[0] - inds[2]) > 0)
    return;

  cone = e;
  assert(cone.size() != 2);
}

bool C4D::VertXYZ::isCompatible (Event *event, FaceABC *f) {
  if (cone.size() == 0)
    return false;
  for (int i = 0; i < cone.size(); i++)
    if (DotEF(event, cone[i], f) > 0)
      return false;
  return true;
}

bool C4D::readObjFile (const char *file, bool moving) {
  assert(SUBF == 14);
  ifstream in;
  in.open(file);
  if (!in)
    return false;
  return readObjFile(in, moving);
}

bool C4D::readObjFile (istream &in, bool moving) {
  int msign = moving ? -1 : 1;
  enable();
  string s;
  while (in >> s) {
    if (s == "v") {
      double x, y, z;
      if (!(in >> x >> y >> z)) {
        disable();
        return false;
      }
#ifdef BLEEN
      extern void sf_rotate (double &x, double &y, double &z);
      sf_rotate(x, y, z);
      cout << "v " << x << " " << y << " " << z << endl;
#endif
      new VertXYZ(this, msign * x, msign * y, msign * z, moving);
    }
    else if (s == "f") {
      int a, b, c;
      if (!(in >> a >> b >> c)) {
        disable();
        return false;
      }
      // new FaceABC(verts[moving][a-1], verts[moving][b-1], verts[moving][c-1]);
      if (!moving)
        new FaceABC(verts[moving][a], verts[moving][b], verts[moving][c]);
      else
        new FaceABC(verts[moving][a], verts[moving][c], verts[moving][b]);
    }
    else {
      disable();
      return false;
    }
  }
  disable();
  return true;
}

void C4D::initSums (bool noCone) {
  enable();
  sumVVs.resize(verts[0].size() * verts[1].size());
  if (!noCone)
    for (int moving = 0; moving < 2; moving++) {
      for (vector<VertXYZ*>::iterator it = verts[moving].begin();
           it != verts[moving].end(); it++)
        (*it)->makeCone();
    }
#ifdef BLEEN
  sumVEs.resize(verts[0].size() * edges[1].size());
  sumEVs.resize(verts[1].size() * edges[0].size());
  sumVFs.resize(verts[0].size() * faces[1].size());
  sumFVs.resize(verts[1].size() * faces[0].size());
  sumEEs.resize(edges[0].size() * edges[1].size());
#endif
  disable();
}

C4D::C4D (C4D *that)  : event0(this, angle0), EVENT0(&event0) {
  for (int m = 0; m < 2; m++) {
    for (int i = 0; i < that->verts[m].size(); i++)
      new VertXYZ(this, that->verts[m][i]);
    for (int i = 0; i < that->faces[m].size(); i++) {
      Face *face = that->faces[m][i];
      const vector<Edge*> &edges = face->getOuterLoop();
      new FaceABC(verts[m][dynamic_cast<VertXYZ*>(edges[0]->tail)->id],
                  verts[m][dynamic_cast<VertXYZ*>(edges[1]->tail)->id],
                  verts[m][dynamic_cast<VertXYZ*>(edges[2]->tail)->id]);
    }
  }
  initSums();
  enable();
}

C4D::Obj *C4D::findObj (Obj *obj, set<C4D::Obj *, C4D::Compare> objs) {
  set<C4D::Obj *, C4D::Compare>::iterator it = objs.find(obj);
  if (it != objs.end())
    return *it;
  return 0;
}

bool C4D::compare (Event *event, set<C4D::Obj *, C4D::Compare> &objs0,
                   set<C4D::Obj *, C4D::Compare> &objs1) {
  static int count = 0;
  if (++count == 11)
    cout << "this is it" << endl;

  if (count == 0) {
    SumFace *faces[] = { (SumFace*) 0x6ab060, (SumFace*) 0x6d1d60, (SumFace*) 0x6b0960};
    for (int i = 0; i < 3; i++) {
      SumFace *face = faces[i];
      const vector<Edge*> &edges = face->getOuterLoop();
      for (int j = 0; j < edges.size(); j++) {
	Edge *edge = edges[j];
	cout << edge->tail;
	if (edge->smallestTwin())
	  for (int k = 0; k < edge->verts.size(); k++) {
	    IntEF *vert = dynamic_cast<IntEF*>(edge->verts[k]);
	    cout << " " << vert << " " << vert->f;
	    if (vert->f == faces[0] || vert->f == faces[1] || vert->f == faces[2])
	      cout << "*";
	  }
	else
	  for (int k = edge->twin->verts.size()-1; k >= 0; k--) {
	    IntEF *vert = dynamic_cast<IntEF*>(edge->twin->verts[k]);
	    cout << " " << vert << " " << vert->f;
	    if (vert->f == faces[0] || vert->f == faces[1] || vert->f == faces[2])
	      cout << "*";
	  }
	cout << endl;
	if (!edge->smallestTwin())
	  edge = edge->twin;
	for (int k = 0; k <= edge->verts.size(); k++)
	  assert(OrderUVE(event, edge->getVert(k-1), edge->getVert(k), edge) == -1);
      }
      cout << endl;
    }
    IntEdge *ff[] = {faces[0]->findLiveEdge(faces[1]),
		     faces[1]->findLiveEdge(faces[2]),
		     faces[2]->findLiveEdge(faces[0]) };
    IntEdge *ff01 = IntEdge::get(faces[0], faces[1]);
    int sTail = faces[2]->side(event, ff[0]->tail);
    int sHead = faces[2]->side(event, ff[0]->twin->tail);
    assert(sTail != sHead);
    if (sTail > 0)
      ff[0] = ff[0]->getTwin();
    bool inside = faces[2]->intersectsInside(event, ff[0]);
    assert(inside);
  }

  set<C4D::Obj *, C4D::Compare>::iterator it0 = objs0.begin();
  set<C4D::Obj *, C4D::Compare>::iterator it1 = objs1.begin();
  Compare compare;
  int n0 = objs0.size();
  int n1 = objs1.size();
  while (it0 != objs0.end() && it1 != objs1.end()) {
    Obj *obj0 = *it0;
    Obj *obj1 = *it1;
    if (obj1 == (Obj*) 0x2c3ff10)
      cout << "this is it" << endl;
    if (obj1->getType() == SUBE) {
      SubEdge *sub1 = dynamic_cast<SubEdge*>(obj1);
      if (sub1->e == (C4D::SumEdge *) 0x2a090a0)
	cout << "this is it" << endl;
    }
    if (compare(obj0, obj1)) {
      cout << "obj0 < obj1" << endl;

      if (obj0->getType() == IFFF) {
	IntFFF *v = dynamic_cast<IntFFF*>(obj0);
	IntFFF *v2 = IntFFF::make(v->faces[0], v->faces[2], v->faces[1]);
	IntEdge *es[3];
	int tailSide[3], headSide[3], dir[3], dirtv[3], dirvh[3];
	for (int i = 0; i < 3; i++) {
	  es[i] = v->faces[(i+1)%3]->findLiveEdge(v->faces[(i+2)%3]);
	  dir[i] = OrderUVE(currentEvent,
			    es[i]->tail, es[i]->twin->tail, es[i]);
	  dirtv[i] = OrderUVE((i == 0 || i == 2) ? event : currentEvent,
			      es[i]->tail, v, es[i]);
	  dirvh[i] = OrderUVE(i == 1 ? event : currentEvent,
			      v, es[i]->twin->tail, es[i]);
	  tailSide[i] = v->faces[i]->side(event, es[i]->tail);
	  headSide[i] = v->faces[i]->side(event, es[i]->twin->tail);
	}
	cout << "IFFF this is it" << endl;
      }

      if (obj0->getType() == SUBE) {
	{
          SubEdge *sub0 = (SubEdge*)obj0;
          Edge *e = sub0->e;
          assert(e->verts.size() == 2);
          Vert *v0 = e->getVert(0);
          Vert *v1 = e->getVert(1);
          int order = OrderUVE(event, v0, v1, e);
          cout << "order " << order << endl;
	}
	{
          SubEdge *sub1 = (SubEdge*)obj1;
          Edge *e = sub1->e;
          assert(e->verts.size() == 2);
          Vert *v0 = e->getVert(0);
          Vert *v1 = e->getVert(1);
          int order = OrderUVE(event, v0, v1, e);
          cout << "order " << order << endl;
	}
	cout << "order this is it" << endl;
      }
      if (obj1->getType() == SUBE) {
	SubEdge *sub1 = dynamic_cast<SubEdge*>(obj1);
	if (sub1->e->verts.size() > 0) {
	  Obj *vert1 = sub1->e->verts[0];
	  Obj *vert0 = findObj(vert1, objs0);
	  cout << "this is it" << endl;
	}
      }
      compare(obj0, obj1);
      if (obj0->getType() == IFFF) {
	IntFFF *fff = dynamic_cast<IntFFF*>(obj0);
	for (int i = 0; i < 3; i++) {
	  SumFace *face = fff->faces[i];
	  const vector<Edge*> &edges = face->getOuterLoop();
	  for (int j = 0; j < edges.size(); j++) {
	    Edge *edge = edges[j];
	    cout << edge->tail;
	    for (int k = 0; k < edge->verts.size(); k++) {
	      IntEF *vert = dynamic_cast<IntEF*>(edge->verts[k]);
	      cout << " " << vert << " " << vert->f;
	    }
	    cout << endl;
	    if (!edge->smallestTwin())
	      edge = edge->twin;
	    for (int k = 0; k <= edge->verts.size(); k++)
	      assert(OrderUVE(event, edge->getVert(k-1), edge->getVert(k), edge) == -1);
	  }
	}
	IntEdge *ff[] = {fff->faces[0]->findLiveEdge(fff->faces[1]),
			 fff->faces[1]->findLiveEdge(fff->faces[2]),
			 fff->faces[2]->findLiveEdge(fff->faces[0]) };
	IntEdge *ff01 = IntEdge::get(fff->faces[0], fff->faces[1]);
	int sTail = fff->faces[2]->side(event, ff[0]->tail);
	int sHead = fff->faces[2]->side(event, ff[0]->twin->tail);
	assert(sTail != sHead);
	if (sTail > 0)
	  ff[0] = ff[0]->getTwin();
	bool inside = fff->faces[2]->intersectsInside(event, ff[0]);
	assert(inside);
      }
    }
    else if (compare(obj1, obj0)) {
      cout << "obj1 < obj0" << endl;
      if (obj0->getType() == SUBE) {
	SubEdge *sub0 = dynamic_cast<SubEdge*>(obj0);

	int sign = OrderUVE(event, sub0->tail, sub0->twin->tail, sub0->e);
	assert(sub0->e->verts.size() >= 2);
	SubEdge *edge = SubEdge::make(sub0->e->getVert(0),
				      sub0->e->getVert(1),
				      sub0->e);
	edge->addEvents();
      }
      if (obj1->getType() == IFFF) {
        IntFFF *fff = dynamic_cast<IntFFF*>(obj1);
	if (false) {
	  IntEdge *ff = fff->faces[0]->findLiveEdge(fff->faces[1]);
	  int sTail = fff->faces[2]->side(event, ff->tail);
	  int sHead = fff->faces[2]->side(event, ff->twin->tail);
	  assert(sTail != sHead);
	  if (sTail > 0)
	    ff = ff->getTwin();
	  bool inside = fff->faces[2]->intersectsInside(event, ff);
	  assert(inside);
	}
	else {
          // (SumFace*) 0x6d0d60, (SumFace*) 0x6af960, (SumFace*) 0x6aa060};
	  SumFace *faces[] =
            { dynamic_cast<SumFace*>(findObj(fff->faces[0], objs0)),
              dynamic_cast<SumFace*>(findObj(fff->faces[1], objs0)),
              dynamic_cast<SumFace*>(findObj(fff->faces[2], objs0)) };
	  for (int i = 0; i < 3; i++) {
	    SumFace *face = faces[i];
	    const vector<Edge*> &edges = face->getOuterLoop();
	    for (int j = 0; j < edges.size(); j++) {
	      SumEdge *edge = dynamic_cast<SumEdge*>(edges[j]);
	      if (!edge->smallestTwin())
		edge = edge->getTwin();
	      for (int k = 1; k < edge->verts.size(); k++) {
		IntEF *v0 = dynamic_cast<IntEF*>(edge->verts[k-1]);
		IntEF *v1 = dynamic_cast<IntEF*>(edge->verts[k]);
		if ((v0->f == faces[(i+1)%3] && v1->f == faces[(i+2)%3]) ||
		    (v1->f == faces[(i+1)%3] && v0->f == faces[(i+2)%3]))
		  cout << "IntFFF this is it" << endl;
	      }
	    }
	  }
	  IntEdge *ff[] = {faces[0]->findLiveEdge(faces[1]),
			   faces[1]->findLiveEdge(faces[2]),
			   faces[2]->findLiveEdge(faces[0]) };
	  IntEdge *ff01 = IntEdge::get(faces[0], faces[1]);
	  int sTail = faces[2]->side(event, ff[0]->tail);
	  int sHead = faces[2]->side(event, ff[0]->twin->tail);
	  assert(sTail != sHead);
	  if (sTail > 0)
	    ff[0] = ff[0]->getTwin();
	  bool inside = faces[2]->intersectsInside(event, ff[0]);
	  assert(inside);
	}
      }
    }

    it0++;
    it1++;
  }
  if (it0 != objs0.end())
    cout << "it0 not at end" << endl;
  if (it1 != objs1.end())
    cout << "it1 not at end" << endl;
  return true;
}

string C4D::Face::toString () {
  string s;
  vector<Edge *> edges = getOuterLoop();
  for (int i = 0; i < edges.size(); i++) {
    s = s + edges[i]->getTail()->toString() + " ";
  }
  return s;
}

bool C4D::Face::contains (Vert *v) {
  const vector<Edge *> &outer = getOuterLoop();
  for (int i = 0; i < outer.size(); i++)
    if (v == outer[i]->getTail())
      return true;

  for (int h = 0; h < numInnerLoops(); h++) {
    const vector<Edge *> &inner = getInnerLoop(h);
    for (int i = 0; i < inner.size(); i++)
      if (v == inner[i]->getTail())
        return true;
  }

  return false;
}

bool C4D::Face::contains (Edge *e) {
  const vector<Edge *> &outer = getOuterLoop();
  for (int i = 0; i < outer.size(); i++)
    if (e == outer[i])
      return true;

  for (int h = 0; h < numInnerLoops(); h++) {
    const vector<Edge *> &inner = getInnerLoop(h);
    for (int i = 0; i < inner.size(); i++)
      if (e == inner[i])
        return true;
  }

  return false;
}

template<class N>
Poly2<N> C4D::PolyEF::calculate () {
  Poly2<N> poly(3);
  poly.m[0] = 0; poly.m[1] = 0;
  poly.m[2] = 1; poly.m[3] = 0;
  poly.m[4] = 0; poly.m[5] = 1;
  poly.setDegree();
  PV3<N> v = e[1]->get<N>() - e[0]->get<N>();
  PV3<N> n = (f[1]->get<N>() - f[0]->get<N>()).cross(f[2]->get<N>() - f[0]->get<N>());
  poly.a[0] = v.z * n.z;
  poly.a[1] = v.x * n.x + v.y * n.y;
  poly.a[2] = e[0]->isMoving() ?
    v.x * n.y - v.y * n.x :
    n.x * v.y - n.y * v.x;
  return poly;
}

C4D::EventSet C4D::PolyEF::getNegative (EdgeAB* e, FaceABC* f) {
  PTR<PolyEF> polyEF = new PolyEF(e, f);
  return EventSet(e->getC4D(), polyEF->getAnglesSorted(), Sign(e, f) < 0);
}

bool C4D::PolyEF::identical (AnglePoly *thatAnglePoly) {
  if (this == thatAnglePoly) return true;
  PolyEF *that = dynamic_cast<PolyEF*>(thatAnglePoly);
  if (that == 0) return false;
  if (this->e[0] != that->e[0]) return false;
  if (this->e[1] != that->e[1]) return false;
  if (this->f[0] != that->f[0]) return false;
  if (this->f[1] != that->f[1]) return false;
  if (this->f[2] != that->f[2]) return false;
  return true;
}

template<class N>
N C4D::PolyEF::Sign::calculate () {
  PV3<N> v = e->Edge::get<N>();
  PV3<N> n = f->Face::get<N>();
  return v.dot(n);
}

template<class N>
N C4D::PolyEF::HasRoots::calculate () {
  PV3<N> v = e->Edge::get<N>();
  PV3<N> n = f->Face::get<N>();
  return (n.z*n.z)*(v.z*v.z) - (n.x*n.x + n.y*n.y)*(v.x*v.x + v.y*v.y);
}

bool C4D::PolyEF::contains (Vert *obj) {
  if (obj->getType() != SVV)
    return false;
  SumVV *vert = (SumVV*) obj;
  int m = e[0]->isMoving();
  if (vert->vv[m] != e[0] && vert->vv[m] != e[1])
    return false;
  if (vert->vv[!m] != f[0] && vert->vv[!m] != f[1] && vert->vv[!m] != f[2])
    return false;
  return true;
}

bool C4D::PolyEF::contains (SumFace *face) {
  const vector<Edge*> &edges = face->getOuterLoop();
  for (int i = 0; i < edges.size(); i++)
    if (!contains(edges[i]->tail))
      return false;
  return true;
}

template<class N>
N C4D::PolyVF::Sign::calculate () {
  return (v->get<N>() - f->getFP<N>()).dot(f->get<N>());
}

template<class N>
Poly2<N> C4D::PolyVF::calculate () {
  SumVV *vFace = dynamic_cast<SumVV*>(f->getOuterLoop()[0]->getTail());
  PV3<N> a = v->vv[0]->get<N>() - vFace->vv[0]->get<N>();
  PV3<N> b = v->vv[1]->get<N>() - vFace->vv[1]->get<N>();
  return C4D::getPoly(a, b, f);
}

C4D::EventSet C4D::PolyVF::getNegative (SumVV *v, SumFace *f) {
  if (f->getType() == SVF) {
    SumVF* svf = dynamic_cast<SumVF*>(f);
    if (svf->f->contains(v->vv[1])) {
      PTR<PolyEF> polyEF = new PolyEF(v->vv[0], svf->v, svf->f);
      return EventSet(v->getC4D(), polyEF->getAnglesSorted(), Sign(v, f) < 0);
    }
  }
  else if (f->getType() == SFV) {
    SumFV* sfv = dynamic_cast<SumFV*>(f);
    if (sfv->f->contains(v->vv[0])) {
      PTR<PolyEF> polyEF = new PolyEF(v->vv[1], sfv->v, sfv->f);
      return EventSet(v->getC4D(), polyEF->getAnglesSorted(), Sign(v, f) < 0);
    }
  }
  else {
    SumEE* see = dynamic_cast<SumEE*>(f);
    for (int i = 0; i <= 1; i++)
      if (see->ee[i]->contains(v->vv[i])) {
        PTR<PolyEF> polyEF =
          new PolyEF(see->ee[i],
                     see->ee[!i]->getTail(), see->ee[!i]->getHead(), v->vv[!i]);
        return EventSet(v->getC4D(), polyEF->getAnglesSorted(), Sign(v, f) < 0);
      }
  }
  
  PTR<PolyVF> polyVF = new PolyVF(v, f);
  return EventSet(v->getC4D(), polyVF->getAnglesSorted(), Sign(v, f) < 0);
}

bool C4D::PolyVF::sameEvent (Event *event, SumVV *v, SumFace *f) {
  AnglePoly *poly = dynamic_cast<RootAngle*>(event->getAngle())->getPoly();

  if (f->getType() == SVF) {
    SumVF* svf = dynamic_cast<SumVF*>(f);
    if (svf->f->contains(v->vv[1])) {
      PolyEF polyEF(v->vv[0], svf->v, svf->f);
      return polyEF.identical(poly);
    }
  }
  else if (f->getType() == SFV) {
    SumFV* sfv = dynamic_cast<SumFV*>(f);
    if (sfv->f->contains(v->vv[0])) {
      PolyEF polyEF(v->vv[1], sfv->v, sfv->f);
      return polyEF.identical(poly);
    }
  }
  else {
    SumEE* see = dynamic_cast<SumEE*>(f);
    for (int i = 0; i <= 1; i++)
      if (see->ee[i]->contains(v->vv[i])) {
        PolyEF polyEF(see->ee[i], see->ee[!i]->getTail(),
                      see->ee[!i]->getHead(), v->vv[!i]);
        return polyEF.identical(poly);
      }
  }
  
  PolyVF polyVF(v, f);
  return polyVF.identical(poly);
}

bool C4D::PolyVF::identical (AnglePoly *thatAnglePoly) {
  if (this == thatAnglePoly) return true;
  PolyVF *that = dynamic_cast<PolyVF*>(thatAnglePoly);
  if (that == 0) return false;
  if (this->v == that->v && this->f == that->f)
    return true;

  return false;
}

template<class N>
N C4D::PolyEE::Sign::calculate () {
  return (e1->getTail()->get<N>() - e0->getTail()->get<N>()).dot(e0->get<N>().cross(e1->get<N>()));
}

template<class N>
Poly2<N> C4D::PolyEE::calculate () {
  Poly2<N> poly(3);
  poly.m[0] = 0; poly.m[1] = 0;
  poly.m[2] = 1; poly.m[3] = 0;
  poly.m[4] = 0; poly.m[5] = 1;
  poly.setDegree();

  SumVV *p0 = e0->getTail();
  SumVV *p1 = e1->getTail();
  PV3<N> a = p1->vv[0]->get<N>() - p0->vv[0]->get<N>();
  PV3<N> b = p1->vv[1]->get<N>() - p0->vv[1]->get<N>();
  PV3<N> v0 = e0->get<N>();
  PV3<N> v1 = e1->get<N>();
  
  if (e0->getType() == SEV)
    if (e1->getType() == SEV) {
      // (a + q b) * (v0 x v1)
      // a * (v0 x v1) + q b * (v0 x v1)
      // (v0 x v1) * q b + a * (v0 x v1)
      PV3<N> v0xv1 = v0.cross(v1);
      setABC(poly, v0xv1, b, a.dot(v0xv1));
    }
    else {
      // (a + q b) * (v0 x q v1)
      // a * (v0 x q v1) + q b * (v0 x q v1)
      // (a x v0) * q v1 + v0 * q (v1 x b)
      setAB(poly, a.cross(v0), v1);
      addAB(poly, v0, v1.cross(b));
    }
  else
    if (e1->getType() == SEV) {
      // (a + q b) * (q v0 x v1)
      // a * (q v0 x v1) + q b * (q v0 x v1)
      // (v1 x a) * q v0 + v1 * q (b x v0)
      setAB(poly, v1.cross(a), v0);
      addAB(poly, v1, b.cross(v0));
    }
    else {
      // (a + q b) * (q v0 x q v1)
      // a * q (v0 x v1) + q b * q (v0 x v1)
      // a * q (v0 x v1) + b * (v0 x v1)
      PV3<N> v0xv1 = v0.cross(v1);
      setABC(poly, a, v0xv1, b.dot(v0xv1));
    }
  return poly;
}

C4D::EventSet C4D::PolyEE::getNegative (SumEdge *e0, SumEdge *e1) {
  if (e0->v->isMoving() == e1->v->isMoving()) {
    VertXYZ* v = 0;
    if (e0->e->contains(e1->e->getTail()))
      v = e1->e->getHead();
    else if (e0->e->contains(e1->e->getHead()))
      v = e1->e->getTail();
    if (v != 0) {
      PTR<PolyEF> polyEF = new PolyEF(e0->v, e1->v,
                                  e0->e->getTail(), e0->e->getHead(), v);
      return EventSet(v->getC4D(), polyEF->getAnglesSorted(), Sign(e0, e1) < 0);
    }
  }
  else {
    SumEdge *ee[] = { e0, e1 };
    for (int i = 0; i <= 1; i++)
      if (ee[i]->e->contains(ee[!i]->v)) {
        PTR<PolyEF> polyEF =
          new PolyEF(ee[i]->e,
                     ee[!i]->e->getTail(), ee[!i]->e->getHead(), ee[i]->v);
        return EventSet(e0->getC4D(), polyEF->getAnglesSorted(), Sign(e0, e1) < 0);
      }
  }
  
  PTR<PolyEE> polyEE = new PolyEE(e0, e1);
  return EventSet(e0->getC4D(), polyEE->getAnglesSorted(), Sign(e0, e1) < 0);
};

bool C4D::PolyEE::sameEvent (Event *event, SumEdge *e0, SumEdge *e1) {
  AnglePoly *poly = dynamic_cast<RootAngle*>(event->getAngle())->getPoly();

  if (e0->v->isMoving() == e1->v->isMoving()) {
    VertXYZ* v = 0;
    if (e0->e->contains(e1->e->getTail()))
      v = e1->e->getHead();
    else if (e0->e->contains(e1->e->getHead()))
      v = e1->e->getTail();
    if (v != 0) {
      PolyEF polyEF(e0->v, e1->v, e0->e->getTail(), e0->e->getHead(), v);
      return polyEF.identical(poly);
    }
  }
  else {
    SumEdge *ee[] = { e0, e1 };
    for (int i = 0; i <= 1; i++)
      if (ee[i]->e->contains(ee[!i]->v)) {
        PolyEF polyEF(ee[i]->e, ee[!i]->e->getTail(),
                      ee[!i]->e->getHead(), ee[i]->v);
        return polyEF.identical(poly);
      }
  }
  
  PolyEE polyEE(e0, e1);
  return polyEE.identical(poly);
};

bool C4D::PolyEE::identical (AnglePoly *thatAnglePoly) {
  if (this == thatAnglePoly) return true;
  PolyEE *that = dynamic_cast<PolyEE*>(thatAnglePoly);
  if (that == 0) return false;
  if ((this->e0 == that->e0 || this->e0 == that->e0->twin)
      &&
      (this->e1 == that->e1 || this->e1 == that->e1->twin))
    return true;

  return false;
}

void C4D::SumEdge::handleEE (Event *event, SumEdge *that) {
  if (!smallestTwin ()) {
    getTwin()->handleEE(event, that);
    return;
  }

  static int count;
  if (++count == 0)
    cout << "handleEE this is it" << endl;

  set<Obj*> *objs[] = { &event->vanish, &event->appear };
  vector<IntEF*> vVanish, vAppear;
  vector<IntEF*> *vs[] = { &vVanish, &vAppear };
  vector<int> sign;
  // i = 0: vanish.  i = 1: appear.
  for (int i = 0; i < 2; i++)
    // for (int j = 0; j < (*objs[i]).size(); j++)
    for (Obj* obj : *objs[i])
      if (obj->getType() == IEF) {
        IntEF *v = dynamic_cast<IntEF*>(obj);
        if ((v->e == this || v->e == getTwin())) {
          if (that->faceIndex(v->f) != -1) {
            vs[i]->push_back(v);
	    if (i == 1)
	      sign.push_back(1);
          }
          else if (that->getTwin()->faceIndex(v->f) != -1) {
            vs[i]->push_back(v);
	    if (i == 1)
	      sign.push_back(-1);
          }
        }
      }

  int minVanish = verts.size(), maxVanish = -1;
  for (int i = 0; i < vVanish.size(); i++) {
    int vIndex = findVert(vVanish[i]);
    assert(vIndex != -1);
    if (minVanish > vIndex)
      minVanish = vIndex;
    if (maxVanish < vIndex)
      maxVanish = vIndex;
  }
  if (minVanish <= maxVanish)
    assert(maxVanish - minVanish + 1 == vVanish.size());
  else
    assert(vVanish.size() == 0);

  for (int j = 1; j < vAppear.size(); j++) {
    IntEF *v = vAppear[j];
    int i = j;
    while (i > 0 &&
           CircEFF(event, that, v->f, vAppear[i-1]->f) ==
           sign[j] * DotEF(event, this, vAppear[i-1]->f)) {
      // int circEFF = CircEFF(event, that, v->f, vAppear[i-1]->f);
      // int dotEF = DotEF(event, this, vAppear[i-1]->f);
      vAppear[i] = vAppear[i-1];
      i--;
    }
    vAppear[i] = v;
  }           

  vector<Vert*> appear;
  for (int i = 0; i < vAppear.size(); i++)
    appear.push_back(vAppear[i]);

  handleReplace(event, minVanish, maxVanish, appear);
}

bool C4D::SumEdge::intersects (Event *event, SumFace *face) {
  IntEF *ef = IntEF::get(this, face);
  if (ef != 0 && ef->life.contains(event))
    return true;
  ef = IntEF::get(getTwin(), face);
  return ef != 0 && ef->life.contains(event);
}  

#ifndef HOMOTOPY_IDENTITY
// positive if v is to the outer side of this face (same as normal)
// zero of v is a vertex of this face
// negative otherwise.
template<class N>
N C4D::SumFace::Side::calculate () {
  Type vType = v->getType();
  if (vType == SVV) {
    for (int i = 0; i < face->edges.size(); i++)
      if (v == face->edges[i]->getTail())
        return N(0);
  }
  else if (vType == IEF) {
    IntEF *ef = dynamic_cast<IntEF*>(v);
    if (ef->f == face)
      return N(0);
    for (int i = 0; i < face->edges.size(); i++)
      if (ef->e == face->edges[i] || ef->e == face->edges[i]->twin)
        return N(0);
  }
  else
    assert(0);

  bool isIdentity = false;
  Angle *angle = event->getAngle();
  if (angle->getType() == Angle::IMPLICIT || 
      angle->getType() == Angle::EXPLICIT) {
    Contact contact = face->getC4D()->getContact(v, face);
    RootAngle *root = dynamic_cast<RootAngle*>(angle);

    if (curPrecision() == 100u)
      nGCD100++;
    else
      nGCD212++;

    isIdentity = face->getC4D()->commonRoot(contact, root);
  }

  if (isIdentity) {
    PV3<N> p = face->getFP<N>(event);
    PV3<N> n = face->get<N>(event);
    PV3<N> q = v->get<N>(event);
    PV3<N> pd = face->getFPd<N>(event);
    PV3<N> nd = face->getd<N>(event);
    PV3<N> qd = v->getd<N>(event);
    return (nd.dot(q - p) + n.dot(qd - pd));
  }

  PV3<N> p = face->getFP<N>(event);
  PV3<N> n = face->get<N>(event);
  PV3<N> q = v->get<N>(event);
  N nDOTpq = n.dot(q - p);
  /*
  int nDOTpqs = nDOTpq.sign(false);
  if (nDOTpqs != 0)
    return nDOTpqs;
  */
  return nDOTpq;
}
#else
template<class N>
N C4D::SumFace::SideH::calculate () {
  PV3<N> p = face->getFP<N>(event);
  PV3<N> n = face->get<N>(event);
  PV3<N> q = v->get<N>(event);
  N nDOTpq = n.dot(q - p);
  /*
  int nDOTpqs = nDOTpq.sign(false);
  if (nDOTpqs != 0)
    return nDOTpqs;
  */
  return nDOTpq;
}

template<class N>
N C4D::SumFace::SideD::calculate () {
  PV3<N> p = face->getFP<N>(event);
  PV3<N> n = face->get<N>(event);
  PV3<N> q = v->get<N>(event);
  PV3<N> pd = face->getFPd<N>(event);
  PV3<N> nd = face->getd<N>(event);
  PV3<N> qd = v->getd<N>(event);
  return (nd.dot(q - p) + n.dot(qd - pd));
}
#endif

// Given e->intersects(this), does it intersect the interior of this?
bool C4D::SumFace::intersectsInside (Event *event, Edge *e) {
  assert(side(event, e->tail) < 0 && side(event, e->twin->tail) > 0);
  for (int i = 0; i < edges.size(); i++)
    if (LeftOf(event, e, edges[i]) < 0)
      return false;
  return true;
}

// True if line through v parallel to e intersects inside face.
// No condition on direction of e.
bool C4D::SumFace::intersectsInside (Event *event, Vert *v, Edge *e) {
  int signEF = DotEF(event, e, this);
  for (int i = 0; i < edges.size(); i++)
    if (LeftVEVE(event, v, e, edges[i]->getTail(), edges[i]) != signEF)
      return false;
  return true;
}

C4D::IntEdge *C4D::SumFace::findLiveEdge (SumFace *that) {
  IntEdge dummy(that);
  set<IntEdge *, Compare>::iterator it = liveEdges.find(&dummy);
  if (it == liveEdges.end())
    return 0;
  return *it;
}

C4D::EventSet C4D::SumFace::getIntLife (SumFace *that) {
  if ((this == (C4D::SumFace *) 0x6e02d0 && that == (C4D::SumFace *) 0x730d80)
      ||
      (that == (C4D::SumFace *) 0x6e02d0 && this == (C4D::SumFace *) 0x730d80))
    cout << "getIntLife this is it" << endl;

  EventSet intLife;
  // cout << "size " << intLife.events.size() << endl;
  IntEF *intEF = 0;
  for (int i = 0; i < this->edges.size(); i++) {
    SumEdge *e = dynamic_cast<SumEdge*>(this->edges[i]);
    if ((intEF = IntEF::get(e, that)) != 0)
      intLife = intLife.combine(intEF->life);
    // cout << "size " << intLife.events.size() << endl;
    e = e->getTwin();
    if ((intEF = IntEF::get(e, that)) != 0)
      intLife = intLife.combine(intEF->life);
    // cout << "size " << intLife.events.size() << endl;
  }
  for (int i = 0; i < that->edges.size(); i++) {
    SumEdge *e = dynamic_cast<SumEdge*>(that->edges[i]);
    if ((intEF = IntEF::get(e, this)) != 0)
      intLife = intLife.combine(intEF->life);
    e = e->getTwin();
    if ((intEF = IntEF::get(e, this)) != 0)
      intLife = intLife.combine(intEF->life);
  }
    
  intLife = intLife.intersect(this->life);
  if (intLife.isNull()) return intLife;
  intLife = intLife.intersect(that->life);
  if (intLife.isNull()) return intLife;

  return intLife;
}

C4D::Vert *C4D::Event::resolveVert (Vert *v) {
  if (v->getType() != IEF)
    return v;
  IntEF *ef = (IntEF *) v;
  RootAngle *root = dynamic_cast<RootAngle*>((Angle *) a);
  if (root == 0)
    return v;
  if (root->getPoly()->getType() != AnglePoly::EF2)
    return v;
  PolyEF *poly = dynamic_cast<PolyEF*>(root->getPoly());
  if (!poly->contains(ef->f))
    return v;
  Vert *vr = v;
  if (poly->contains(ef->e->tail))
    vr = ef->e->tail;
  if (poly->contains(ef->e->twin->tail))
    vr = ef->e->twin->tail;
  static int count;
  if (v != vr && ++count == 0) {
    int sign = OrderUVE(this, v, vr, ef->e);
    cout << sign << endl;
  }
  return vr;
}  

void C4D::SumFace::intersectEdges (Event *event, bool doHandle) {
  RootAngle *root = dynamic_cast<RootAngle*>(event->getAngle());
  PolyEF *poly = (root != 0 && root->getPoly()->getType() == AnglePoly::EF2) ?
    dynamic_cast<PolyEF*>(root->getPoly()) : 0;
  bool polyContainsThis = (poly != 0 && poly->contains(this));

  for (set<IntEdge *, Compare>::iterator ie = liveEdges.begin();
       ie != liveEdges.end(); ++ie) {
    IntEdge *ab = *ie;
    Vert *ar = event->resolveVert(ab->getTail());
    Vert *br = event->resolveVert(ab->getHead());
    SumFace *f = ab->getFace();

    if (polyContainsThis && f->hasCoplanarEdge(poly))
      continue;

    for (set<IntEdge *, Compare>::iterator je = ie;
         ++je != liveEdges.end();) {
      IntEdge *cd = *je;
      Vert *cr = event->resolveVert(cd->getTail());
      Vert *dr = event->resolveVert(cd->getHead());
      SumFace *g = cd->getFace();

      if (polyContainsThis && g->hasCoplanarEdge(poly))
	continue;

      if (f->neighborOf(g))
	continue;

      int gSa = g->side(event, ar);
      int gSb = g->side(event, br);
      if (gSa * gSb != -1)
        continue;

      int fSc = f->side(event, cr);
      int fSd = f->side(event, dr);
      if (fSc * fSd != -1)
        continue;

      assert(gSa * fSc < 0);
      IntFFF *v =
        gSa < 0 ? IntFFF::make(this, f, g) : IntFFF::make(this, g, f);
      ab->addVert(event, v);
      cd->addVert(event, v);

      IntEdge *fg = f->findLiveEdge(g);

      assert(findLiveEdge(f));
      assert(findLiveEdge(g));
      assert(f->findLiveEdge(this));
      assert(g->findLiveEdge(this));
      assert(f->findLiveEdge(g));
      assert(g->findLiveEdge(f));

      if (!doHandle)
	fg->addVert(event, v);
      else
	fg->handleInsert(event, v);
    }
  }
}

C4D::SumFace *C4D::VertXYZ::getSum (FaceABC *f) {
  static int count;
  count++;
  if (cone.size() == 0)
    return 0;

  EventSet events(true);

  for (int i = 0; i < cone.size(); i++) {
    EdgeAB *e = dynamic_cast<EdgeAB*>(cone[i]);
    if (PolyEF::hasRoots(e, f)) {
      EventSet eventsI = PolyEF::getNegative(e, f);
      events = events.intersect(eventsI);
      if (events.isNull())
        return 0;
    }
    else if (PolyEF::Sign(e, f) > 0)
      return 0;
  }
  
  SumFace *sum = 0;
  if (isMoving())
    sum = SumFV::make(f, this);
  else
    sum = SumVF::make(this, f);
  sum->life = events;
  sum->updateEdgeLives();
  return sum;
}  

C4D::SumFace *C4D::EdgeAB::getSum (EdgeAB *that) {
  EdgeAB *edges[] = { that, that, getTwin(), getTwin() };
  FaceABC *faces[] = { getFace(), getTwin()->getFace(), that->getFace(), that->getTwin()->getFace() };
  int signs[] = { 1, -1, 1, -1 };

  EventSet events(true);
  for (int i = 0; i < 4; i++) {
    EdgeAB *e = edges[i];
    FaceABC *f = faces[i];
    int sign = signs[i];

    if (PolyEF::hasRoots(e, f)) {
      // PolyEF *polyEF = new PolyEF(e, f);
      // EventSet eventsI = polyEF->getNegative();
      EventSet eventsI = PolyEF::getNegative(e, f);
      if (sign > 0)
        eventsI = eventsI.complement();
      events = events.intersect(eventsI);
      if (events.isNull())
        return 0;
    }
    else if ((PolyEF::Sign(e, f) > 0) != (sign > 0))
      return 0;
  }

  SumFace *sum = SumEE::make(this, that);
  sum->life = events;
  sum->updateEdgeLives();
  return sum;
}

template<class N>
Poly2<N> C4D::getPoly (PV3<N> a, PV3<N> b, SumFace *f) {
  Poly2<N> poly(3);
  poly.m[0] = 0; poly.m[1] = 0;
  poly.m[2] = 1; poly.m[3] = 0;
  poly.m[4] = 0; poly.m[5] = 1;
  poly.setDegree();
  
  if (f->getType() == SFV) {
    PV3<N> n = f->get<N>();
    // (a + q b) * n = n * q b + a * n
    setABC(poly, n, b, a.dot(n));
  }
  else if (f->getType() == SVF) {
    PV3<N> n = f->get<N>();
    // (a + q b) * q n = a * (q n) + b * n
    setABC(poly, a, n, b.dot(n));
  }
  else if (f->getType() == SEE) {
    SumEE *see = dynamic_cast<SumEE*>(f);
    PV3<N> v0 = see->ee[0]->get<N>();
    PV3<N> v1 = see->ee[1]->get<N>();
    // (a + q b) * (v0 x q v1)
    // a * (v0 x q v1) + (q b) * (v0 x q v1)
    // (a x v0) * q v1 + v0 * q (v1 x b)
    setAB(poly, a.cross(v0), v1);
    addAB(poly, v0, v1.cross(b));
  }
  
  return poly;
}

template<class N>
Poly2<N> C4D::getPoly (PV3<N> a, PV3<N> b) {
  Poly2<N> poly(3);
  poly.m[0] = 0; poly.m[1] = 0;
  poly.m[2] = 1; poly.m[3] = 0;
  poly.m[4] = 0; poly.m[5] = 1;
  poly.setDegree();
  setAB(poly, a, b);
  return poly;
}

template<class N>
Poly2<N> C4D::getPolyA (PV3<N> a, SumFace *f) {
  if (f->getType() == SFV) {
    Poly2<N> poly(1);
    poly.m[0] = 0; poly.m[1] = 0;
    poly.setDegree();
    // a * n
    poly.a[0] = a.dot(f->get<N>());
    return poly;
  }

  Poly2<N> poly(3);
  poly.m[0] = 0; poly.m[1] = 0;
  poly.m[2] = 1; poly.m[3] = 0;
  poly.m[4] = 0; poly.m[5] = 1;
  poly.setDegree();
  
  if (f->getType() == SVF) {
    // a * q n
    setAB(poly, a, f->get<N>());
  }
  else if (f->getType() == SEE) {
    SumEE *see = dynamic_cast<SumEE*>(f);
    PV3<N> v0 = see->ee[0]->get<N>();
    PV3<N> v1 = see->ee[1]->get<N>();
    // a * (v0 x q v1)
    // (a x v0) * q v1
    setAB(poly, a.cross(v0), v1);
  }
  
  return poly;
}

template<class N>
Poly2<N> C4D::getPolyB (PV3<N> b, SumFace *f) {
  if (f->getType() == SVF) {
    Poly2<N> poly(1);
    poly.m[0] = 0; poly.m[1] = 0;
    poly.setDegree();
    // q b * q n = b * n
    poly.a[0] = b.dot(f->get<N>());
    return poly;
  }

  Poly2<N> poly(3);
  poly.m[0] = 0; poly.m[1] = 0;
  poly.m[2] = 1; poly.m[3] = 0;
  poly.m[4] = 0; poly.m[5] = 1;
  poly.setDegree();
  
  if (f->getType() == SFV) {
    // q b * n = n * q b
    setAB(poly, f->get<N>(), b);
  }
  else if (f->getType() == SEE) {
    SumEE *see = dynamic_cast<SumEE*>(f);
    PV3<N> v0 = see->ee[0]->get<N>();
    PV3<N> v1 = see->ee[1]->get<N>();
    // q b * (v0 x q v1)
    // v0 * q (v1 x b)
    setAB(poly, v0, v1.cross(b));
  }
  
  return poly;
}

template<class N>
Poly2<N> C4D::getPoly (SumEdge *e, SumFace *f) {
  PV3<N> v = e->get<N>();
  if (e->getType() == SEV)
    return getPolyA(v, f);
  else
    return getPolyB(v, f);
}

template<class N>
Poly2<N> C4D::PolyEEE::calculate () {
  for (int i = 0; i < 3; i++) {
    int j = (i+1)%3;
    int k = (j+1)%3;
    if (!e[i][0]->isMoving() && !e[j][0]->isMoving() && e[k][0]->isMoving())
      return C4D::getPoly(getV<N>(i).cross(getV<N>(j)), getV<N>(k));
    if (e[i][0]->isMoving() && e[j][0]->isMoving() && !e[k][0]->isMoving())
      return C4D::getPoly(getV<N>(k), getV<N>(i).cross(getV<N>(j)));
  }

  assert(0);
  return Poly2<N>();
}
  
template<class N>
Poly2<N> C4D::PolyFFF::calculate () {
  for (int i = 0; i < 3; i++) {
    int j = (i+1)%3;
    int k = (j+1)%3;
    if (nMoving(i) == 0 && nMoving(j) == 0)
      return getPolyA(getN<N>(i).cross(getN<N>(j)), k);
  }

  for (int i = 0; i < 3; i++) {
    int j = (i+1)%3;
    int k = (j+1)%3;
    if (nMoving(i) == 3 && nMoving(j) == 3)
      return getPolyB(getN<N>(i).cross(getN<N>(j)), k);
  }

  for (int i = 0; i < 3; i++) {
    int j = (i+1)%3;
    int k = (j+1)%3;
    if (nMoving(i) == 2) {
      PV3<N> v = getV<N>(i, 0);
      PV3<N> w = getV<N>(i, 1);
      return (getPolyA(v, j) * getPolyB(w, k) -
              getPolyA(v, k) * getPolyB(w, j));
    }
  }
  
  assert(0);
  return Poly2<N>();
}

template<class N>
Poly2<N> C4D::PolyFFF::getPolyA (PV3<N> a, int i) {
  if (nMoving(i) == 0) {
    Poly2<N> poly(1);
    poly.m[0] = 0; poly.m[1] = 0;
    poly.setDegree();
    // a * n
    poly.a[0] = a.dot(getN<N>(i));
    return poly;
  }

  Poly2<N> poly(3);
  poly.m[0] = 0; poly.m[1] = 0;
  poly.m[2] = 1; poly.m[3] = 0;
  poly.m[4] = 0; poly.m[5] = 1;
  poly.setDegree();
  
  if (nMoving(i) == 3) {
    // a * q n
    setAB(poly, a, getN<N>(i));
  }
  else if (nMoving(i) == 2) {
    PV3<N> v0 = getV<N>(i, 0);
    PV3<N> v1 = getV<N>(i, 1);
    // a * (v0 x q v1)
    // (a x v0) * q v1
    setAB(poly, a.cross(v0), v1);
  }
  
  return poly;
}

template<class N>
Poly2<N> C4D::PolyFFF::getPolyB (PV3<N> b, int i) {
  if (nMoving(i) == 3) {
    Poly2<N> poly(1);
    poly.m[0] = 0; poly.m[1] = 0;
    poly.setDegree();
    // q b * q n = b * n
    poly.a[0] = b.dot(getN<N>(i));
    return poly;
  }

  Poly2<N> poly(3);
  poly.m[0] = 0; poly.m[1] = 0;
  poly.m[2] = 1; poly.m[3] = 0;
  poly.m[4] = 0; poly.m[5] = 1;
  poly.setDegree();
  
  if (nMoving(i) == 0) {
    // q b * n = n * q b
    setAB(poly, getN<N>(i), b);
  }
  else if (nMoving(i) == 2) {
    PV3<N> v0 = getV<N>(i, 0);
    PV3<N> v1 = getV<N>(i, 1);
    // q b * (v0 x q v1)
    // v0 * q (v1 x b)
    setAB(poly, v0, v1.cross(b));
  }
  
  return poly;
}

template<class N>
Poly2<N> C4D::getPoly (SumEdge *e0, SumEdge *e1, SumEdge *e2) {
  SumEdge *es[] = { e0, e1, e2 };

  for (int i = 0; i < 3; i++) {
    int j = (i+1)%3;
    int k = (j+1)%3;
    if (es[i]->getType() == SEV && 
        es[j]->getType() == SEV && 
        es[k]->getType() == SVE)
      return getPoly(es[i]->get<N>().cross(es[j]->get<N>()), es[k]->get<N>());
    if (es[i]->getType() == SVE && 
        es[j]->getType() == SVE && 
        es[k]->getType() == SEV)
      return getPoly(es[k]->get<N>(), es[i]->get<N>().cross(es[j]->get<N>()));
  }

  assert(0);
  return Poly2<N>();
}

template<class N>
Poly2<N> C4D::getPoly (SumFace *f0, SumFace *f1, SumFace *f2) {
  SumFace *fs[] = { f0, f1, f2 };

  for (int i = 0; i < 3; i++) {
    int j = (i+1)%3;
    int k = (j+1)%3;
    if (fs[i]->getType() == SFV && fs[j]->getType() == SFV)
      return getPolyA(fs[i]->get<N>().cross(fs[j]->get<N>()), fs[k]);
  }

  for (int i = 0; i < 3; i++) {
    int j = (i+1)%3;
    int k = (j+1)%3;
    if (fs[i]->getType() == SVF && fs[j]->getType() == SVF)
      return getPolyB(fs[i]->get<N>().cross(fs[j]->get<N>()), fs[k]);
  }

  for (int i = 0; i < 3; i++) {
    int j = (i+1)%3;
    int k = (j+1)%3;
    if (fs[i]->getType() == SEE) {
      SumEE *see = dynamic_cast<SumEE*>(fs[i]);
      PV3<N> v = see->ee[0]->get<N>();
      PV3<N> w = see->ee[1]->get<N>();
      return (getPolyA(v, fs[j]) * getPolyB(w, fs[k]) -
              getPolyA(v, fs[k]) * getPolyB(w, fs[j]));
    }
  }

  assert(0);
  return Poly2<N>();
}

#ifdef BLEEN
Poly2 C4D::PolyEFF::getPoly () {
  /*
    SumEdge p0, V0
    SumFace p1, N1
    SumFace p2, N2
    p = p0 - p1
    q = p - p2
    
    (p + t V0) * N1 = 0
    p * N1 = -V0 * N1 t
    t = -(p * N1) / (V0 * N1)
    p - (p * N1) / (V0 * N1) V0
    (q - (p * N1) / (V0 * N1) V0) * N2
    
    (q * N2) (V0 * N1) - (p * N1) (V0 * N2)
    over
    V0 * N1
  */

  PV3<N> p0 = ef->e->getTail()->getP();
  PV3<N> p1 = ef->f->getP();
  PV3<N> p2 = f->getP();
  PV3<N> p = p0 - p1;
  PV3<N> q = p - p2;

  return (C4D::getPolyA(q, f) * C4D::getPoly(ef->e, ef->f) -
          C4D::getPolyA(p, ef->f) * C4D::getPoly(ef->e, f));
}
#endif

template<class N>
Poly2<N> C4D::PolyEFF::calculate () {
  /*
    SumEdge p0, V0
    SumFace p1, N1
    SumFace p2, N2
    p10 = p0 - p1
    p20 = p0 - p2
    
    p0 + t V0
    (p0 + t V0 - p1) * N1 = 0
    (p10 + t V0) * N1 = 0
    p10 * N1 = -V0 * N1 t
    t = -(p10 * N1) / (V0 * N1)
    p0 - (p10 * N1) / (V0 * N1) V0
    (p0 - (p10 * N1) / (V0 * N1) V0 - p2) * N2
    (p20 - (p10 * N1) / (V0 * N1) V0) * N2
    (p20 * N2) (V0 * N1) - (p10 * N1) (V0 * N2)
    over
    V0 * N1
  */

  SumVV *ev = dynamic_cast<SumVV*>(e->getTail());
  SumVV *f1v = dynamic_cast<SumVV*>(f1->getOuterLoop()[0]->getTail());
  SumVV *f2v = dynamic_cast<SumVV*>(f2->getOuterLoop()[0]->getTail());
  PV3<N> p0a = ev->vv[0]->get<N>();
  PV3<N> p0b = ev->vv[1]->get<N>();
  PV3<N> p1a = f1v->vv[0]->get<N>();
  PV3<N> p1b = f1v->vv[1]->get<N>();
  PV3<N> p2a = f2v->vv[0]->get<N>();
  PV3<N> p2b = f2v->vv[1]->get<N>();
  PV3<N> p10a = p0a - p1a;
  PV3<N> p10b = p0b - p1b;
  PV3<N> p20a = p0a - p2a;
  PV3<N> p20b = p0b - p2b;

  Poly2<N> polyA = C4D::getPoly<N>(p20a, p20b, f2);
  Poly2<N> polyB = C4D::getPoly<N>(e, f1);
  Poly2<N> polyC = C4D::getPoly<N>(p10a, p10b, f1);
  Poly2<N> polyD = C4D::getPoly<N>(e, f2);
  
  return polyA * polyB - polyC * polyD;

  /*
  return (C4D::getPoly(p20a, p20b, f2) * C4D::getPoly(e, f1) -
          C4D::getPoly(p10a, p10b, f1) * C4D::getPoly(e, f2));
  */
}

bool C4D::PolyEFF::isEpEEE (EdgeAB *&ab, EdgeAB *cd[3]) {
#ifndef NOTIME
  idenTime -= getTime();
#endif
  bool ret = isEpEEE_(ab, cd);
#ifndef NOTIME
  idenTime += getTime();
#endif
  return ret;
}

bool C4D::PolyEFF::isEpEEE_ (EdgeAB *&ab, EdgeAB *cd[3]) {
  if (!(e->getType() == SVE &&
        f1->getType() == SEE &&
        f2->getType() == SEE))
    return false;

  EdgeAB *ab1 = dynamic_cast<SumEE*>(f1)->ee[0];
  if (!ab1->contains(dynamic_cast<SumVE*>(e)->v))
    return false;

  EdgeAB *ab2 = dynamic_cast<SumEE*>(f2)->ee[0];
  if (ab1 != ab2 && ab1 != ab2->twin)
    return false;

  ab = ab1->smallestTwin() ? ab1 : ab1->getTwin();

  cd[0] = dynamic_cast<SumVE*>(e)->e;
  cd[1] = dynamic_cast<SumEE*>(f1)->ee[1];
  cd[2] = dynamic_cast<SumEE*>(f2)->ee[1];
  for (int i = 0; i < 3; i++)
    if (!cd[i]->smallestTwin())
      cd[i] = cd[i]->getTwin();

  for (int i = 0; i < 2; i++)
    for (int j = i+1; j < 3; j++)
      if (cd[j]->lessThan(cd[i])) {
        EdgeAB *temp = cd[i];
        cd[i] = cd[j];
        cd[j] = temp;
      }
  
  return true;
}

bool C4D::PolyEFF::isEEEpE (EdgeAB *cd[3], EdgeAB *&ab) {
#ifndef NOTIME
  idenTime -= getTime();
#endif
  bool ret = isEEEpE_(cd, ab);
#ifndef NOTIME
  idenTime += getTime();
#endif
  return ret;
}

bool C4D::PolyEFF::isEEEpE_ (EdgeAB *cd[3], EdgeAB *&ab) {
  if (!(e->getType() == SEV &&
        f1->getType() == SEE &&
        f2->getType() == SEE))
    return false;

  EdgeAB *ab1 = dynamic_cast<SumEE*>(f1)->ee[1];
  if (!ab1->contains(dynamic_cast<SumEV*>(e)->v))
    return false;

  EdgeAB *ab2 = dynamic_cast<SumEE*>(f2)->ee[1];
  if (ab1 != ab2 && ab1 != ab2->twin)
    return false;

  ab = ab1->smallestTwin() ? ab1 : ab1->getTwin();

  cd[0] = dynamic_cast<SumEV*>(e)->e;
  cd[1] = dynamic_cast<SumEE*>(f1)->ee[0];
  cd[2] = dynamic_cast<SumEE*>(f2)->ee[0];
  for (int i = 0; i < 3; i++)
    if (!cd[i]->smallestTwin())
      cd[i] = cd[i]->getTwin();

  for (int i = 0; i < 2; i++)
    for (int j = i+1; j < 3; j++)
      if (cd[j]->lessThan(cd[i])) {
        EdgeAB *temp = cd[i];
        cd[i] = cd[j];
        cd[j] = temp;
      }
  
  return true;
}

bool C4D::PolyFFFF::isEpEEE (EdgeAB *&ab, EdgeAB *cd[3]) {
#ifndef NOTIME
  idenTime -= getTime();
#endif
  bool ret = isEpEEE_(ab, cd);
#ifndef NOTIME
  idenTime += getTime();
#endif
  return ret;
}

bool C4D::PolyFFFF::isEpEEE_ (EdgeAB *&ab, EdgeAB *cd[3]) {
  if (!(f[1]->getType() == SEE &&
        f[2]->getType() == SEE &&
        f[3]->getType() == SEE))
    return false;

  EdgeAB *ab0123[4];
  ab0123[0] = f[0]->getType() == SEE ? dynamic_cast<SumEE*>(f[0])->ee[0] : 0;
  for (int i = 1; i < 4; i++)
    ab0123[i] = dynamic_cast<SumEE*>(f[i])->ee[0];
  for (int i = 0; i < 4; i++)
    if (ab0123[i] != 0 && !ab0123[i]->smallestTwin())
      ab0123[i] = ab0123[i]->getTwin();

  int i0;
  for (i0 = 0; i0 < 4; i0++)
    if (ab0123[i0] == ab0123[(i0+1)%4] &&
        ab0123[i0] == ab0123[(i0+2)%4])
      break;
  if (i0 == 4)
    return false;

  ab = ab0123[i0];

  for (int i = 0; i < 3; i++) {
    cd[i] = dynamic_cast<SumEE*>(f[(i0+i)%4])->ee[1];
    if (!cd[i]->smallestTwin())
      cd[i] = cd[i]->getTwin();
  }

  for (int i = 0; i < 2; i++)
    for (int j = i+1; j < 3; j++)
      if (cd[j]->lessThan(cd[i])) {
        EdgeAB *temp = cd[i];
        cd[i] = cd[j];
        cd[j] = temp;
      }
  
  return true;
}

bool C4D::PolyFFFF::isEEEpE (EdgeAB *cd[3], EdgeAB *&ab) {
#ifndef NOTIME
  idenTime -= getTime();
#endif
  bool ret = isEEEpE_(cd, ab);
#ifndef NOTIME
  idenTime += getTime();
#endif
  return ret;
}

bool C4D::PolyFFFF::isEEEpE_ (EdgeAB *cd[3], EdgeAB *&ab) {
  if (!(f[1]->getType() == SEE &&
        f[2]->getType() == SEE &&
        f[3]->getType() == SEE))
    return false;

  EdgeAB *ab0123[4];
  ab0123[0] = f[0]->getType() == SEE ? dynamic_cast<SumEE*>(f[0])->ee[1] : 0;
  for (int i = 1; i < 4; i++)
    ab0123[i] = dynamic_cast<SumEE*>(f[i])->ee[1];
  for (int i = 0; i < 4; i++)
    if (ab0123[i] != 0 && !ab0123[i]->smallestTwin())
      ab0123[i] = ab0123[i]->getTwin();

  int i0;
  for (i0 = 0; i0 < 4; i0++)
    if (ab0123[i0] == ab0123[(i0+1)%4] &&
        ab0123[i0] == ab0123[(i0+2)%4])
      break;
  if (i0 == 4)
    return false;

  ab = ab0123[i0];

  for (int i = 0; i < 3; i++) {
    cd[i] = dynamic_cast<SumEE*>(f[(i0+i)%4])->ee[0];
    if (!cd[i]->smallestTwin())
      cd[i] = cd[i]->getTwin();
  }

  for (int i = 0; i < 2; i++)
    for (int j = i+1; j < 3; j++)
      if (cd[j]->lessThan(cd[i])) {
        EdgeAB *temp = cd[i];
        cd[i] = cd[j];
        cd[j] = temp;
      }
  
  return true;
}

void C4D::addEvents (SumEdge *e, SumFace *f1, SumFace *f2) {
  static int count;
  count++;
  PolyEFF temp(e, f1, f2);
#ifndef NOTIME
  idenTime -= getTime();
#endif
  bool found = (polyEFFs.find(&temp) != polyEFFs.end());
#ifndef NOTIME
  idenTime += getTime();
#endif
  if (found)
    return;
  PolyEFF *polyEFF = new PolyEFF(e, f1, f2);
  polyEFFs.insert(polyEFF);
  // polyEFF->incRef();

  // get 2 x 2 IntEFs...
  IntEF *ef1s[] = { IntEF::get(e, f1), IntEF::get(e->getTwin(), f1) };
  IntEF *ef2s[] = { IntEF::get(e, f2), IntEF::get(e->getTwin(), f2) };

#ifndef NOTIME
  idenTime -= getTime();
#endif

  PolyEF *polyEF = f1->related(f2);

#ifndef NOTIME
  idenTime += getTime();
#endif

  if (polyEF != 0) {
    vector<PTR<Angle>> roots = polyEF->getAnglesSorted();
    for (int i = 0; i < roots.size(); i++) {
      Angle *root = roots[i];
      IntEF *ef1 = 0;
      if (ef1s[0] && ef1s[0]->life.contains(root))
        ef1 = ef1s[0];
      else if (ef1s[1] && ef1s[1]->life.contains(root))
        ef1 = ef1s[1];
      else {
        delete root;
        continue;
      }
      IntEF *ef2 = 0;
      if (ef2s[0] && ef2s[0]->life.contains(root))
        ef2 = ef2s[0];
      else if (ef2s[1] && ef2s[1]->life.contains(root))
        ef2 = ef2s[1];
      else {
        delete root;
        continue;
      }
      Event temp(this, root);
      Event *event = 0;
      set<Event*, Compare>::iterator it = eventList.find(&temp);
      if (it != eventList.end()) {
        delete root;
        event = *it;
      }
      else {
        event = new Event(this, root);
        eventList.insert(event);
      }
      if (event == (C4D::Event *) 0x1ca1b40)
	cout << "swap this is it" << endl;
      event->swap.insert(pair<Obj*,Obj*>(ef1,ef2));
    }
    return;
  }
 
#ifndef NOTIME
  idenTime -= getTime();
#endif
  bool f1nf2 = f1->neighborOf(f2);
#ifndef NOTIME
  idenTime += getTime();
#endif
  if (f1nf2)
    return;

#ifndef NOTIME
  idenTime -= getTime();
#endif
  SumVV *svv = f1->commonVert(f2);
  if (svv != 0) {
    SumVV *tail = e->getTail();
    SumVV *head = e->getTwin()->getTail();
    if (tail->vv[0] == svv->vv[0] || tail->vv[1] == svv->vv[1] ||
        head->vv[0] == svv->vv[0] || head->vv[1] == svv->vv[1]) {
#ifndef NOTIME
      idenTime += getTime();
#endif
      PolyNNN(f1, f2, e).addEvents();
      return;
    }
  }

#ifndef NOTIME
  idenTime += getTime();
#endif

  BigPoly *poly = polyEFF;

  {
    EdgeAB *ab, *cd[3];
    if (polyEFF->isEpEEE(ab, cd)) {
      cout << "isEpEEE " << count << endl;
      set<BigPoly*, BigPoly::CompareEpEEE>::iterator
        it = polyEpEEEs.find(polyEFF);
      if (it == polyEpEEEs.end())
        polyEpEEEs.insert(polyEFF);
      else
        poly = *it;
    }
  }

  {
    EdgeAB *ab[3], *cd;
    if (polyEFF->isEEEpE(ab, cd)) {
      cout << "isEEEpE " << count << endl;
      set<BigPoly*, BigPoly::CompareEEEpE>::iterator
        it = polyEEEpEs.find(polyEFF);
      if (it == polyEEEpEs.end())
        polyEEEpEs.insert(polyEFF);
      else
        poly = *it;
    }
  }

  vector<PTR<Angle>> roots = poly->getAnglesSorted();

  for (int i = 0; i < roots.size(); i++) {
    Angle *root = roots[i];
    /*
      Event temp(this, root, true);
      int s = OrderUVE(&temp,
      ef1s[0] ? ef1s[0] : ef1s[1],
      ef2s[0] ? ef2s[0] : ef2s[1],
      e);
    */
    // int s2 = RootAngle::Sign(polyEFF, root);
    if (((ef1s[0] && ef1s[0]->life.contains(root)) ||
         (ef1s[1] && ef1s[1]->life.contains(root))) &&
        ((ef2s[0] && ef2s[0]->life.contains(root)) ||
         (ef2s[1] && ef2s[1]->life.contains(root)))) {
      Event temp(this, root);
      Event *event = 0;
      if (eventList.find(&temp) != eventList.end()) {
	event = *eventList.find(&temp);
        delete root;
      }
      else {
	event = new Event(this, root);
	eventList.insert(event);
      }
      event->appear.insert(polyEFF->e);
      event->appear.insert(polyEFF->f1);
      event->appear.insert(polyEFF->f2);
      event->vanish.insert(polyEFF->e);
    }
  }
}

C4D::IntEF *C4D::findLiveEF (SumEdge *e, SumFace *f) {
  for (int i = 0; i < 2; i++) {
    IntEF *v = IntEF::get(e, f, sweepSets[IEF]);
    if (v != 0)
      return v;
    e = e->getTwin();
  }
  assert(0);
  return 0;
}

#ifdef BLEEN
int C4D::Edge::findVert (Obj *vert) {
  for (int i = 0; i < verts.size(); i++)
    if (verts[i] == vert)
      return i;
  assert(0);
  return -1;
}
#endif

// Sub from verts[i] to verts[i+1] vanishes. i can be -1 or size-1.
void C4D::Edge::vanishSub (Event *event, int i) {
  if (!smallestTwin()) {
    assert(i == -1);
    twin->vanishSub(event, twin->verts.size()-1);
    return;
  }
  SubEdge *sub = SubEdge::make(getVert(i), getVert(i+1), this);
  assert(sub->life.events.size() % 2 == 1);
  sub->life.events.push_back(event);
  event->vanish.insert(sub);
  assert(sub->getTwin()->life.events.size() % 2 == 1);
  sub->getTwin()->life.events.push_back(event);
  event->vanish.insert(sub->getTwin());
}

// Sub from verts[i] to verts[i+1] appears. i can be -1 or size-1.
void C4D::Edge::appearSub (Event *event, int i) {
  if (!smallestTwin()) {
    assert(i == -1);
    twin->appearSub(event, twin->verts.size()-1);
    return;
  }
  SubEdge *sub = SubEdge::make(getVert(i), getVert(i+1), this);
  assert(sub->life.events.size() % 2 == 0);
  sub->life.events.push_back(event);
  event->appear.insert(sub);
  assert(sub->getTwin()->life.events.size() % 2 == 0);
  sub->getTwin()->life.events.push_back(event);
  event->appear.insert(sub->getTwin());
  if (0 <= i && i+1 < verts.size())
    sub->addEvents();
}

void C4D::Edge::vanishSubs (Event *event) {
  for (int i = -1; i < 0 || i < verts.size(); i++)
    vanishSub(event, i);
}

void C4D::Edge::appearSubs (Event *event) {
  for (int i = -1; i < 0 || i < verts.size(); i++)
    appearSub(event, i);
}

// Swap verts[i1] and verts[i1+1]
void C4D::Edge::handleSwap (Event *event, int i1) {
  int i2 = i1+1;
  for (int i = i1-1; i < i2+1; i++)
    vanishSub(event, i);
    
  Obj *temp = verts[i1];
  verts[i1] = verts[i2];
  verts[i2] = temp;
  for (int i = i1-1; i < i2+1; i++)
    appearSub(event, i);
}

void C4D::Edge::handleSwap (Event *event, Vert *v1, Vert *v2) {
  if (!smallestTwin()) {
    getTwin()->handleSwap(event, v1, v2);
    return;
  }
  int i1 = findVert(v1);
  assert(i1 != -1);
  int i2 = findVert(v2);
  assert(i2 != -1);
  if (i1+1 == i2)
    handleSwap(event, i1);
  else if (i2+1 == i1)
    handleSwap(event, i2);
  else
    assert(0);
}  

void C4D::Event::makeVanish (Vert *v) {
  if (v->getType() == IFFF && v->life.events.back() != this) {
    if (v == (IntFFF*) 0x2aa4980)
      cout << "this is it" << endl;
    v->life.events.push_back(this);
    if (v == (C4D::Obj *) 0xdcf7080)
      cout << "vVanish this is it" << endl;
    vanish.insert(v);
  }
}

void C4D::Event::makeAppear (Vert *v) {
  if (v == (Vert*) 0x5e414fa0)
    cout << "makeAppear this is it" << endl;
  if (v->getType() == IFFF &&
      (v->life.events.size() == 0 || v->life.events.back() != this)) {
    v->life.events.push_back(this);
    appear.insert(v);
  }
}

// Add or remove v at tail end.
void C4D::Edge::handleTail (Event *event, Vert *v) {
  if (!smallestTwin())
    twin->handleHead(event, v);
  else if (verts.size() > 0 && verts[0] == v) {
    event->makeVanish(v);
    vanishSub(event, -1);
    vanishSub(event, 0);
    verts.erase(verts.begin());
    appearSub(event, -1);
  }
  else {
    event->makeAppear(v);
    vanishSub(event, -1);
    verts.insert(verts.begin(), v);
    appearSub(event, -1);
    appearSub(event, 0);
  }
}

// Add or remove v at head end.
void C4D::Edge::handleHead (Event *event, Vert *v) {
  if (!smallestTwin())
    twin->handleTail(event, v);
  else if (verts.size() > 0 && verts.back() == v) {
    event->makeVanish(v);
    vanishSub(event, verts.size()-2);
    vanishSub(event, verts.size()-1);
    verts.pop_back();
    appearSub(event, verts.size()-1);
  }
  else {
    event->makeAppear(v);
    vanishSub(event, verts.size()-1);
    verts.push_back(v);
    appearSub(event, verts.size()-2);
    appearSub(event, verts.size()-1);
  }
}

void C4D::Edge::handleRemove (Event *event, Vert *v) {
  if (!smallestTwin()) {
    twin->handleRemove(event, v);
    return;
  }
  event->makeVanish(v);
  int index = findVert(v);
  assert(index != -1);
  vanishSub(event, index-1);
  vanishSub(event, index);
  verts.erase(verts.begin()+index);
  appearSub(event, index-1);
}

void C4D::Edge::handleInsert (Event *event, Vert *v) {
  if (!smallestTwin()) {
    twin->handleInsert(event, v);
    return;
  }
  event->makeAppear(v);
  int index = findVert(v);
  assert(index == -1);
  int a = 0, b = verts.size()-1;
  while (a <= b) {
    int c = (a+b)/2;
    if (OrderUVE(event, v, dynamic_cast<Vert*>(verts[c]), this) < 0)
      b = c - 1;
    else
      a = c + 1;
  }
  vanishSub(event, a-1);
  verts.insert(verts.begin() + a, v);
  appearSub(event, a-1);
  appearSub(event, a);
}

void C4D::Edge::handleReplace (Event *event, Vert *v0, Vert *v1) {
  if (!smallestTwin()) {
    twin->handleReplace(event, v0, v1);
    return;
  }
  event->makeVanish(v0);
  event->makeAppear(v1);
  int index = findVert(v0);
  assert(index != -1);
  assert(findVert(v1) == -1);
    
  vanishSub(event, index-1);
  vanishSub(event, index);
  verts[index] = v1;
  appearSub(event, index-1);
  appearSub(event, index);
}

void C4D::Edge::handleReplace (Event *event, int minVanish, int maxVanish,
                               vector<Vert*> &appear) {
  if (minVanish > maxVanish && appear.size() == 0)
    return;

  if (!smallestTwin()) {
    vector<Vert*> reverse;
    for (int i = appear.size(); --i >= 0;)
      reverse.push_back(appear[i]);
    twin->handleReplace(event, minVanish, maxVanish, reverse);
    return;
  }

  if (minVanish <= maxVanish) {
    for (int i = minVanish; i <= maxVanish; i++)
      event->makeVanish(getVert(i));
    for (int i = minVanish-1; i <= maxVanish; i++)
      vanishSub(event, i);
    verts.erase(verts.begin()+minVanish, verts.begin()+maxVanish+1);
  }

  int index = minVanish;
  if (minVanish > maxVanish) {
    assert(appear.size() > 0);
    // set index
    Vert *v = appear[0];
    int a = 0, b = verts.size()-1;
    while (a <= b) {
      int c = (a+b)/2;
      if (OrderUVE(event, v, dynamic_cast<Vert*>(verts[c]), this) < 0)
	b = c - 1;
      else
	a = c + 1;
    }
    index = a;
    vanishSub(event, index-1);
  }
  for (int i = 0; i < appear.size(); i++)
    event->makeAppear(appear[i]);
  verts.insert(verts.begin()+index, appear.begin(), appear.end());
  int last = index + appear.size();
  for (int i = index-1; i < last; i++)
    appearSub(event, i);
}

void C4D::Edge::handleRemove (Event *event, Vert *v0, Vert *v1) {
  if (!smallestTwin()) {
    twin->handleRemove(event, v0, v1);
    return;
  }
  event->makeVanish(v0);
  event->makeVanish(v1);
  int index = findVert(v0);
  assert(index != -1);
  if (index+1 == verts.size() || verts[index+1] != v1) {
    assert(index > 0 && verts[index-1] == v1);
    index--;
  }

  vanishSub(event, index-1);
  vanishSub(event, index);
  vanishSub(event, index+1);
  verts.erase(verts.begin()+index);
  verts.erase(verts.begin()+index);
  appearSub(event, index-1);
}

void C4D::Edge::handleInsert (Event *event, Vert *v0, Vert *v1) {
  if (!smallestTwin()) {
    twin->handleInsert(event, v1, v0);
    return;
  }
  assert(findVert(v0) == -1);
  assert(findVert(v1) == -1);
  event->makeAppear(v0);
  event->makeAppear(v1);
  int index = insertIndex(event, v0);
  assert(index != -1);
    
  vanishSub(event, index-1);
  verts.insert(verts.begin()+index, v0);
  verts.insert(verts.begin()+index+1, v1);
  appearSub(event, index-1);
  appearSub(event, index);
  appearSub(event, index+1);
}

#ifdef BLEEN
void C4D::SumEdge::handleEFF (Event *event) {
  assert(smallestTwin());
  assert(event->vanish[0] == this);
  assert(event->appear[0] == this);
  SumFace *f[] = { dynamic_cast<SumFace*>(event->appear[1]),
                   dynamic_cast<SumFace*>(event->appear[2]) };
  event->vanish.erase(event->vanish.begin());
  event->appear.erase(event->appear.begin());
  event->appear.erase(event->appear.begin());
  event->appear.erase(event->appear.begin());

  C4D *c4d = getC4D();
  IntEF *ef[] = { c4d->findLiveEF(this, f[0]), c4d->findLiveEF(this, f[1]) };
  // int order = C4D::OrderUVE(event, ef[0], ef[1], this);
  int i[] = { findVert(ef[0]), findVert(ef[1]) };
  assert(i[0] != -1 && i[1] != -1);
  if (i[0] < i[1]) {
    assert(i[0]+1 == i[1]);
    handleSwap(event, i[0]);
  }
  else {
    assert(i[1]+1 == i[0]);
    handleSwap(event, i[1]);
  }

  // Two end events for each face incident on this.
  IntFFF *v[2] = { 0, 0 };
  int s[2] = { 0, 0 };
  int n = 0;
  SumEdge *e = this;
  for (int i = 0; i < 2; i++) {
    if (e->face != 0) {
      s[n] = (CircFFF(event, e->face, f[0], f[1]) == 1);
      v[n] = IntFFF::make(e->getFace(), f[!s[n]], f[s[n]]);

      for (int j = 0; j < 2; j++) {
        IntEdge *ff = e->getFace()->findLiveEdge(f[j]);
        if (ff->getTail() == ef[j])
          ff->handleTail(event, v[n]);
        else {
          assert(ff->getTwin()->getTail() == ef[j]);
          ff->handleHead(event, v[n]);
        }
      }

      n++;
    }
    e = e->getTwin();
  }

  IntEdge *ff = f[0]->findLiveEdge(f[1]);
  if (ff->containsVert(v[0])) {
    if (v[1] == 0)
      ff->handleRemove(event, v[0]);
    else if (ff->containsVert(v[1]))
      ff->handleRemove(event, v[0], v[1]);
    else
      ff->handleReplace(event, v[0], v[1]);
  }
  else {
    if (v[1] == 0)
      ff->handleInsert(event, v[0]);
    else if (ff->containsVert(v[1]))
      ff->handleReplace(event, v[1], v[0]);
    else if (DotEF(event, ff, face) == -CircEFF(event, this, face, twin->face))
      ff->handleInsert(event, v[0], v[1]);
    else
      ff->handleInsert(event, v[1], v[0]);
  }
}
#endif

void C4D::SumEdge::handleEFF (Event *event) {
  assert(0);
#ifdef SWEEP_ALL_FACES
  assert(smallestTwin());
  assert(event->vanish[0] == this);
  assert(event->appear[0] == this);
  SumFace *f[] = { dynamic_cast<SumFace*>(event->appear[1]),
                   dynamic_cast<SumFace*>(event->appear[2]) };
  event->vanish.erase(event->vanish.begin());
  event->appear.erase(event->appear.begin());
  event->appear.erase(event->appear.begin());
  event->appear.erase(event->appear.begin());

  handleEFF(event, f, false);
#endif
}

void C4D::SumEdge::handleEFF (Event *event, SumFace *f[2], bool f0parallel) {
  C4D *c4d = getC4D();
  IntEF *ef[] = { 0, c4d->findLiveEF(this, f[1]) };
  if (!f0parallel)
    ef[0] = c4d->findLiveEF(this, f[0]);

  if (!f0parallel) {
    int i[] = { findVert(ef[0]), findVert(ef[1]) };
    assert(i[0] != -1 && i[1] != -1);
    if (i[0] < i[1]) {
      assert(i[0]+1 == i[1]);
      handleSwap(event, i[0]);
    }
    else {
      assert(i[1]+1 == i[0]);
      handleSwap(event, i[1]);
    }
  }

  PolyEF *poly = f0parallel ?
    dynamic_cast<PolyEF*>(dynamic_cast<RootAngle*>(event->getAngle())->getPoly()) : 0;

  IntEdge *ff = f[0]->findLiveEdge(f[1]);
  // Two end events for each face incident on this.

  SumVV *commonV = 0;
  if (ff->tail->getType() == SVV)
    commonV = dynamic_cast<SumVV*>(ff->tail);
  else if (ff->twin->tail->getType() == SVV)
    commonV = dynamic_cast<SumVV*>(ff->twin->tail);
  

  vector<IntFFF*> vVanish, vAppear;
  vector<int> sign;
  SumEdge *e = this;
  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < e->faces.size(); k++) {
      SumFace *face = e->faces[k];
      
      if (poly != 0 &&
	  (// face is parallel to f[0].
	   poly->contains(face) ||
	   // Intersection with f[0] vanished not appeared.
	   face->findLiveEdge(f[0]) == 0))
	continue;
      
      if (commonV != 0 && face->contains(commonV))
	continue;
      
      // identity
      // Three SEEs with the same fixed (or moving) edge.
      if (face->getType() == SEE &&
	  f[0]->getType() == SEE &&
	  f[1]->getType() == SEE) {
	SumEE *see = dynamic_cast<SumEE*>(face);
	SumEE *see0 = dynamic_cast<SumEE*>(f[0]);
	SumEE *see1 = dynamic_cast<SumEE*>(f[1]);
	if ((see->ee[0] == see0->ee[0] || see->ee[0] == see0->ee[0]->twin) &&
	    (see->ee[0] == see1->ee[0] || see->ee[0] == see1->ee[0]->twin))
	  continue;
	if ((see->ee[1] == see0->ee[1] || see->ee[1] == see0->ee[1]->twin) &&
	    (see->ee[1] == see1->ee[1] || see->ee[1] == see1->ee[1]->twin))
	  continue;
      }

      int s = (CircFFF(event, face, f[0], f[1]) == 1);
      IntFFF *v = IntFFF::make(face, f[!s], f[s]);
      if (ff->containsVert(v)) {
	assert(!f0parallel);
        vVanish.push_back(v);
      }
      else if (!(f0parallel && f[0]->findLiveEdge(face) == 0)) {
        vAppear.push_back(v);
        sign.push_back(1 - 2 * i);
      }

      for (int j = 0; j < 2; j++) {
        IntEdge *facef = face->findLiveEdge(f[j]);
        if (j == 0 && f0parallel)
	  facef->handleInsert(event, v);
        else if (facef->getTail() == ef[j])
          facef->handleTail(event, v);
        else {
          assert(facef->getTwin()->getTail() == ef[j]);
          facef->handleHead(event, v);
        }
      }
    }
    e = e->getTwin();
  }

  int minVanish = ff->verts.size(), maxVanish = -1;
  if (minVanish < ff->twin->verts.size())
    minVanish = ff->twin->verts.size();
  for (int i = 0; i < vVanish.size(); i++) {
    int vIndex = ff->findVert(vVanish[i]);
    assert(vIndex != -1);
    if (minVanish > vIndex)
      minVanish = vIndex;
    if (maxVanish < vIndex)
      maxVanish = vIndex;
  }
  if (minVanish <= maxVanish)
    /* assert(maxVanish - minVanish + 1 == vVanish.size()) */;
  else
    assert(vVanish.size() == 0);

  for (int j = 1; j < vAppear.size(); j++) {
    IntFFF *v = vAppear[j];
    int i = j;
    while (i > 0 &&
           CircEFF(event, this, v->otherFace(f[0], f[1]), vAppear[i-1]->otherFace(f[0], f[1])) ==
           sign[j] * DotEF(event, ff, vAppear[i-1]->otherFace(f[0], f[1]))) {
      vAppear[i] = vAppear[i-1];
      i--;
    }
    vAppear[i] = v;
  }           

  vector<Vert*> vAppear2;
  for (int i = 0; i < vAppear.size(); i++)
    vAppear2.push_back(vAppear[i]);

  ff->handleReplace(event, minVanish, maxVanish, vAppear2);
}

#ifdef BLEEN
C4D::PolyEEE::PolyEEE (SumEdge *e0, SumEdge *e1, SumEdge *e2) {
  e[0] = e0;
  int i;
  for (i = 1; i > 0 && C4D::Compare()(e1, e[i-1]); i--)
    e[i] = e[i-1];
  e[i] = e1;
  for (i = 2; i > 0 && C4D::Compare()(e2, e[i-1]); i--)
    e[i] = e[i-1];
  e[i] = e2;
}

C4D::PolyFFF::PolyFFF (SumFace *f0, SumFace *f1, SumFace *f2) {
  f[0] = f0;
  int i;
  for (i = 1; i > 0 && C4D::Compare()(f1, f[i-1]); i--)
    f[i] = f[i-1];
  f[i] = f1;
  for (i = 2; i > 0 && C4D::Compare()(f2, f[i-1]); i--)
    f[i] = f[i-1];
  f[i] = f2;
}
#endif

C4D::PolyNNN::Pair::Pair (SumVV *v, SumFace *f) {
  const vector<Edge *> &edges = f->getOuterLoop();
  int i = 0;
  while (edges[i]->tail != v) i++;
  SumVV *a = dynamic_cast<SumVV*>(edges[i]->twin->tail);
  SumVV *b = dynamic_cast<SumVV*>(edges[(i + edges.size() - 1) % edges.size()]->tail);
  setVV(a->vv[a->vv[0] == v->vv[0]], b->vv[b->vv[0] == v->vv[0]]);
}

C4D::PolyNNN::Pair::Pair (SumVV *v, SumEdge *e) {
  SumVV *t = e->getTail();
  SumVV *h = e->getTwin()->getTail();

  if (t->vv[0] == v->vv[0] && h->vv[0] == v->vv[0])
    setVV(t->vv[1], h->vv[1]);
  else if (t->vv[1] == v->vv[1] && h->vv[1] == v->vv[1])
    setVV(t->vv[0], h->vv[0]);
  else if (t->vv[0] == v->vv[0] || t->vv[1] == v->vv[1])
    setVV(h->vv[0], h->vv[1]);
  else if (h->vv[0] == v->vv[0] || h->vv[1] == v->vv[1])
    setVV(t->vv[0], t->vv[1]);
  else
    assert(0);
}      

C4D::PolyNNN::PolyNNN (SumFace *f0, SumFace *f1, SumFace *f2) {
  vert = f0->commonVert(f1);
  assert(vert != 0);
  assert(f2->contains(vert));

  pairs[0] = Pair(vert, f0);
  Pair pair1(vert, f1);
  if (pairs[0] < pair1)
    pairs[1] = pair1;
  else {
    pairs[1] = pairs[0];
    pairs[0] = pair1;
  }
  Pair pair2(vert, f2);
  if (pairs[1] < pair2)
    pairs[2] = pair2;
  else {
    pairs[2] = pairs[1];
    if (pairs[0] < pair2)
      pairs[1] = pair2;
    else {
      pairs[1] = pairs[0];
      pairs[0] = pair2;
    }
  }
}

C4D::PolyNNN::PolyNNN (SumFace *f0, SumFace *f1, SumEdge *e) {
  vert = f0->commonVert(f1);
  assert(vert != 0);

  pairs[0] = Pair(vert, f0);
  Pair pair1(vert, f1);
  if (pairs[0] < pair1)
    pairs[1] = pair1;
  else {
    pairs[1] = pairs[0];
    pairs[0] = pair1;
  }
  Pair pair2(vert, e);
  if (pairs[1] < pair2)
    pairs[2] = pair2;
  else {
    pairs[2] = pairs[1];
    if (pairs[0] < pair2)
      pairs[1] = pair2;
    else {
      pairs[1] = pairs[0];
      pairs[0] = pair2;
    }
  }
}

bool C4D::PolyNNN::identical (AnglePoly *thatAP) {
  if (thatAP->getType() != NNN)
    return false;
  PolyNNN *that = dynamic_cast<PolyNNN*>(thatAP);
  return (this->vert->vv[0] == that->vert->vv[0] &&
          this->vert->vv[1] == that->vert->vv[1] &&
          this->pairs[0] == that->pairs[0] &&
          this->pairs[1] == that->pairs[1] &&
          this->pairs[2] == that->pairs[2]);
}

template<class N>
Poly2<N> C4D::PolyNNN::getPolyA (PV3<N> a, Pair &pair) {
  PV3<N> p0 = vert->vv[0]->get<N>();
  PV3<N> p1 = vert->vv[1]->get<N>();

  if (!pair.vv[1]->isMoving()) {
    PV3<N> n = (pair.vv[0]->get<N>() - p0).cross(pair.vv[1]->get<N>() - p0);
    Poly2<N> poly(1);
    poly.m[0] = 0; poly.m[1] = 0;
    poly.setDegree();
    // a * n
    poly.a[0] = a.dot(n);
    return poly;
  }

  Poly2<N> poly(3);
  poly.m[0] = 0; poly.m[1] = 0;
  poly.m[2] = 1; poly.m[3] = 0;
  poly.m[4] = 0; poly.m[5] = 1;
  poly.setDegree();
  
  if (pair.vv[0]->isMoving()) {
    PV3<N> n = (pair.vv[0]->get<N>() - p1).cross(pair.vv[1]->get<N>() - p1);
    // a * q n
    setAB(poly, a, n);
  }
  else {
    PV3<N> v0 = pair.vv[0]->get<N>() - p0;
    PV3<N> v1 = pair.vv[1]->get<N>() - p1;
    // a * (v0 x q v1)
    // (a x v0) * q v1
    setAB(poly, a.cross(v0), v1);
  }
  
  return poly;
}

template<class N>
Poly2<N> C4D::PolyNNN::getPolyB (PV3<N> b, Pair &pair) {
  PV3<N> p0 = vert->vv[0]->get<N>();
  PV3<N> p1 = vert->vv[1]->get<N>();

  if (pair.vv[0]->isMoving()) {
    PV3<N> n = (pair.vv[0]->get<N>() - p1).cross(pair.vv[1]->get<N>() - p1);
    Poly2<N> poly(1);
    poly.m[0] = 0; poly.m[1] = 0;
    poly.setDegree();
    // q b * q n = b * n
    poly.a[0] = b.dot(n);
    return poly;
  }

  Poly2<N> poly(3);
  poly.m[0] = 0; poly.m[1] = 0;
  poly.m[2] = 1; poly.m[3] = 0;
  poly.m[4] = 0; poly.m[5] = 1;
  poly.setDegree();
  
  if (!pair.vv[1]->isMoving()) {
    PV3<N> n = (pair.vv[0]->get<N>() - p0).cross(pair.vv[1]->get<N>() - p0);
    // q b * n = n * q b
    setAB(poly, n, b);
  }
  else {
    PV3<N> v0 = pair.vv[0]->get<N>() - p0;
    PV3<N> v1 = pair.vv[1]->get<N>() - p1;
    // q b * (v0 x q v1)
    // v0 * q (v1 x b)
    setAB(poly, v0, v1.cross(b));
  }
  
  return poly;
}

template<class N>
Poly2<N> C4D::PolyNNN::calculate () {
  C4D *c4d = vert->getC4D();
  PV3<N> p0 = vert->vv[0]->get<N>();
  PV3<N> p1 = vert->vv[1]->get<N>();
  if (!pairs[1].vv[1]->isMoving()) {
    assert(pairs[2].vv[1]->isMoving());
    PV3<N> n0 = (pairs[0].vv[0]->get<N>() - p0).cross(pairs[0].vv[1]->get<N>() - p0);
    PV3<N> n1 = (pairs[1].vv[0]->get<N>() - p0).cross(pairs[1].vv[1]->get<N>() - p0);
    PV3<N> n = n0.cross(n1);
    // n * pairs[2]
    return getPolyA(n, pairs[2]);
  }
  else if (pairs[1].vv[0]->isMoving()) {
    assert(!pairs[0].vv[0]->isMoving());
    PV3<N> n1 = (pairs[1].vv[0]->get<N>() - p1).cross(pairs[1].vv[1]->get<N>() - p1);
    PV3<N> n2 = (pairs[2].vv[0]->get<N>() - p1).cross(pairs[2].vv[1]->get<N>() - p1);
    PV3<N> n = n1.cross(n2);
    // q n * pairs[0]
    return getPolyB(n, pairs[0]);
  }
  else {
    PV3<N> v0 = pairs[1].vv[0]->get<N>() - p0;
    PV3<N> v1 = pairs[1].vv[1]->get<N>() - p1;
    // (n0 x (v0 x q v1)) * n2
    // ((n0 * q v1) v0 - (n0 * v0) q v1) * n2
    // (q v1 * n0) (v0 * n2) - (v0 * n0) (q v1 * n2)
    return (getPolyB(v1, pairs[0]) * getPolyA(v0, pairs[2]) -
            getPolyA(v0, pairs[0]) * getPolyB(v1, pairs[2]));
  }
}

void C4D::PolyNNN::getFaces (Pair &pair, SumFace *faces[2]) {
  faces[0] = faces[1] = 0;
  if (!pair.vv[0]->isMoving() && !pair.vv[1]->isMoving()) {
    VertXYZ *a = vert->vv[0];
    VertXYZ *b = pair.vv[0];
    VertXYZ *c = pair.vv[1];
    VertXYZ *d = vert->vv[1];
    // abc + d
    FaceABC *abc = a->getFace(b, c);
    if (abc != 0)
      faces[0] = SumFV::get(abc, d);
  }
  else if (pair.vv[0]->isMoving() && pair.vv[1]->isMoving()) {
    VertXYZ *a = vert->vv[0];
    VertXYZ *b = vert->vv[1];
    VertXYZ *c = pair.vv[0];
    VertXYZ *d = pair.vv[1];
    // a + bcd
    FaceABC *bcd = b->getFace(c, d);
    if (bcd != 0)
      faces[0] = SumVF::get(a, bcd);
  }
  else {
    VertXYZ *a = vert->vv[0];
    VertXYZ *b = pair.vv[0];
    VertXYZ *c = vert->vv[1];
    VertXYZ *d = pair.vv[1];
    // ab + cd
    EdgeAB *ab = a->id < b->id ? a->getEdge(b) : b->getEdge(a);
    if (ab == 0)
      return;
    EdgeAB *cd = c->getEdge(d);
    if (cd == 0)
      return;
    int n = 0;
    for (int i = 0; i < 2; i++) {
      if ((faces[n] = SumEE::get(ab, cd)) != 0)
        n++;
      cd = cd->getTwin();
    }
  }
}

void C4D::PolyNNN::getEdges (Pair &pair, SumEdge *edges[2]) {
  edges[0] = edges[1] = 0;
  if (!pair.vv[0]->isMoving() && !pair.vv[1]->isMoving()) {
    VertXYZ *a = pair.vv[0];
    VertXYZ *b = pair.vv[1];
    VertXYZ *c = vert->vv[1];
    // ab + c
    EdgeAB *ab = a->getEdge(b);
    if (ab != 0)
      edges[0] = SumEV::get(ab, c);
  }
  else if (pair.vv[0]->isMoving() && pair.vv[1]->isMoving()) {
    VertXYZ *a = vert->vv[0];
    VertXYZ *b = pair.vv[0];
    VertXYZ *c = pair.vv[1];
    // a + bc
    EdgeAB *bc = b->getEdge(c);
    if (bc != 0)
      edges[0] = SumVE::get(a, bc);
  }
  else {
    VertXYZ *a = vert->vv[0];
    VertXYZ *b = pair.vv[0];
    VertXYZ *c = vert->vv[1];
    VertXYZ *d = pair.vv[1];
    EdgeAB *ab = a->getEdge(b);
    // ab + d
    if (ab != 0)
      edges[0] = SumEV::get(ab, d);
    EdgeAB *cd = c->getEdge(d);
    // b + cd
    if (cd != 0)
      edges[edges[0] != 0] = SumVE::get(b, cd);
  }
}

void C4D::PolyNNN::addEvents () {
  static int count;
  if (++count == 0)
    cout << "PolyNNN::addEvents this is it" << endl;

  C4D *c4d = vert->getC4D();
  set<PTR<PolyNNN>, PolyNNN::Compare> &polyNNNs = c4d->polyNNNs;
  if (polyNNNs.find(this) != polyNNNs.end())
    return;
  PolyNNN *poly = new PolyNNN(*this);
  polyNNNs.insert(poly);
  // poly->incRef();

  SumFace *faces[3][2];
  SumEdge *edges[3][2];
  for (int i = 0; i < 3; i++) {
    getFaces(pairs[i], faces[i]);
    getEdges(pairs[i], edges[i]);
  }

  vector<PTR<Angle>> roots = poly->getAnglesSorted();
  for (int h = 0; h < roots.size(); h++) {
    Angle *root = roots[h];
    Event temp(c4d, root);
    bool valid = true;

    SumFace *face[3] = { 0, 0, 0 };
    int nfaces = 0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 2; j++)
        if (faces[i][j] != 0 && faces[i][j]->life.contains(root))
          face[i] = faces[i][j];
      if (face[i] != 0)
        nfaces++;
    }

    if (nfaces < 2) {
      delete root;
      continue;
    }

    if (nfaces == 2)
      for (int i = 0; i < 3; i++) {
        if (face[i] != 0)
          continue;
        SumFace *face1 = face[(i+1)%3];
        SumFace *face2 = face[(i+2)%3];
        IntEdge *ff = IntEdge::get(face1, face2);
        if (ff == 0 || !ff->life.contains(root))
          break;
        for (int j = 0; j < 2; j++) {
          SumEdge *edge = edges[i][j];
          if (edge != 0 && edge->life.contains(root) &&
              edge->intersects(&temp, face1) &&
              edge->intersects(&temp, face2))
            nfaces++;
        }
      }


    if (nfaces == 2) {
      delete root;
      continue;
    }
      
    Event *event = new Event(c4d, root);
    c4d->eventList.insert(event);
  }    
}

void C4D::PolyNNN::handleEvent (Event *event) {
  int nfaces = 0;
  SumFace *face[3] = { 0, 0, 0 };
  for (int i = 0; i < 3; i++) {
    SumFace *faces[2];
    getFaces(pairs[i], faces);
    for (int j = 0; j < 2; j++)
      if (faces[j] != 0 && faces[j]->life.contains(event))
        face[i] = faces[j];
    if (face[i] != 0)
      nfaces++;
  }
  assert(nfaces >= 2);

  if (nfaces == 3) {
    // for each pair of FF ff1 and ff2
    // for each pair of FFF fff1 on ff1 and fff2 on ff2
    // if they have the same other plane
    // get their common FF ff
    // swap fff1 and fff2 on ff
    for (int i = 0; i < 3; i++) {
      int i1 = (i+1)%3;
      int i2 = (i+2)%3;
      IntEdge *ff1 = face[i]->findLiveEdge(face[i1]);
      IntEdge *ff2 = face[i]->findLiveEdge(face[i2]);
      if (ff1 == 0 || ff2 == 0)
        continue;
      if (!ff1->smallestTwin())
        ff1 = ff1->getTwin();
      if (!ff2->smallestTwin())
        ff2 = ff2->getTwin();
      for (int j1 = 0; j1 < ff1->verts.size(); j1++) {
        IntFFF *fff1 = dynamic_cast<IntFFF*>(ff1->verts[j1]);
        SumFace *other = fff1->otherFace(face[i], face[i1]);
        for (int j2 = 0; j2 < ff2->verts.size(); j2++) {
          IntFFF *fff2 = dynamic_cast<IntFFF*>(ff2->verts[j2]);
          if (fff2->otherFace(face[i], face[i2]) != other)
            continue;
          IntEdge *ff = face[i]->findLiveEdge(other);
          ff->handleSwap(event, fff1, fff2);
        }
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    if (!(nfaces == 3 || face[i] == 0))
      continue;
    int i1 = (i+1)%3;
    int i2 = (i+2)%3;
    IntEdge *ff = face[i1]->findLiveEdge(face[i2]);
    if (ff == 0)
      continue;
    SumFace *f12[2] = { face[i1], face[i2] };
    SumEdge *edges[2];
    getEdges(pairs[i], edges);
    for (int j = 0; j < 2; j++) {
      SumEdge *edge = edges[j];
      if (edge != 0 && !edge->smallestTwin())
	edge = edge->getTwin();
      if (edge != 0 && edge->life.contains(event) &&
          edge->intersects(event, face[i1]) &&
          edge->intersects(event, face[i2]))
        edge->handleEFF(event, f12, false);
    }
  }
}

C4D::PolyFFFF::PolyFFFF (SumFace *f0, SumFace *f1, SumFace *f2, SumFace *f3) {
  f[0] = f0;
  int i;
  for (i = 1; i > 0 && C4D::Compare()(f1, f[i-1]); i--)
    f[i] = f[i-1];
  f[i] = f1;
  for (i = 2; i > 0 && C4D::Compare()(f2, f[i-1]); i--)
    f[i] = f[i-1];
  f[i] = f2;
  for (i = 3; i > 0 && C4D::Compare()(f3, f[i-1]); i--)
    f[i] = f[i-1];
  f[i] = f3;
}

template<class N>
Poly2<N> C4D::PolyFFFF::calculate () {
  SumVV *f0v = dynamic_cast<SumVV*>(f[0]->getOuterLoop()[0]->getTail());
  SumVV *f1v = dynamic_cast<SumVV*>(f[1]->getOuterLoop()[0]->getTail());
  SumVV *f2v = dynamic_cast<SumVV*>(f[2]->getOuterLoop()[0]->getTail());
  SumVV *f3v = dynamic_cast<SumVV*>(f[3]->getOuterLoop()[0]->getTail());

  PV3<N> p0a = f0v->vv[0]->get<N>();
  PV3<N> p0b = f0v->vv[1]->get<N>();
  PV3<N> p1a = f1v->vv[0]->get<N>();
  PV3<N> p1b = f1v->vv[1]->get<N>();
  PV3<N> p2a = f2v->vv[0]->get<N>();
  PV3<N> p2b = f2v->vv[1]->get<N>();
  PV3<N> p3a = f3v->vv[0]->get<N>();
  PV3<N> p3b = f3v->vv[1]->get<N>();

#ifdef BLEEN
  return ((C4D::getPoly(p0a, p0b, f[0]) * C4D::getPoly(f[1], f[2], f[3]) -
           C4D::getPoly(p1a, p1b, f[1]) * C4D::getPoly(f[0], f[2], f[3]) +
           C4D::getPoly(p2a, p2b, f[2]) * C4D::getPoly(f[0], f[1], f[3]) -
           C4D::getPoly(p3a, p3b, f[3]) * C4D::getPoly(f[0], f[1], f[2])));
#endif
  Poly2<N> poly0 = C4D::getPoly<N>(p0a, p0b, f[0]), poly123 = C4D::getPoly<N>(f[1], f[2], f[3]);
  Poly2<N> poly1 = C4D::getPoly<N>(p1a, p1b, f[1]), poly023 = C4D::getPoly<N>(f[0], f[2], f[3]);
  Poly2<N> poly2 = C4D::getPoly<N>(p2a, p2b, f[2]), poly013 = C4D::getPoly<N>(f[0], f[1], f[3]);
  Poly2<N> poly3 = C4D::getPoly<N>(p3a, p3b, f[3]), poly012 = C4D::getPoly<N>(f[0], f[1], f[2]);

  Poly2<N> ret = poly0 * poly123 - poly1 * poly023 + poly2 * poly013 - poly3 * poly012;

  static int count;
  ++count;
  // cout << "PolyFFFF " << count << endl;
  bool hasZeros = ret.hasZeros();

  return ret.removeZeros();

    // return poly0 * poly123 - poly1 * poly023 + poly2 * poly013 - poly3 * poly012;

#ifdef BLEEN
  SumEdge *e0 = dynamic_cast<SumEdge*>(f0->getOuterLoop()[0]);
  PV3<N> p0 = e0->getTail()->getP();
  PV3<N> p1 = f1->getP();
  PV3<N> p2 = f2->getP();
  PV3<N> p3 = f2->getP();
  PV3<N> p = p0 - p1;
  PV3<N> q = p - p2;
  PV3<N> r = q - p3;

  // ((r * N3) (V0 * N1) - (p * N1) (V0 * N3)) [N0 N1 N2]
  // -
  // ((q * N2) (V0 * N1) - (p * N1) (V0 * N2)) [N0 N1 N3]

  return ((C4D::getPolyA<N>(r, f3) * C4D::getPoly<N>(e0, f1) -
           C4D::getPolyA<N>(p, f1) * C4D::getPoly<N>(e0, f3)) *
          C4D::getPoly<N>(f0, f1, f2) -
          (C4D::getPolyA<N>(q, f2) * C4D::getPoly<N>(e0, f1) -
           C4D::getPolyA<N>(p, f1) * C4D::getPoly<N>(e0, f2)) *
          C4D::getPoly<N>(f0, f1, f3));
#endif
}

void C4D::addEvents (SumFace *fe1, SumFace *fe2, SumFace *fv1, SumFace *fv2) {
  static int count;
  count++;
  PolyFFFF temp(fe1, fe2, fv1, fv2);
  if (temp.f[0] == (SumFace*) 0x786a40 &&
      temp.f[1] == (SumFace*) 0xa6ade0 &&
      temp.f[2] == (SumFace*) 0xac8200 &&
      temp.f[3] == (SumFace*) 0xaf4b50)
    cout << "addEvents this is it" << endl;
#ifndef NOTIME
  idenTime -= getTime();
#endif
  bool found = (polyFFFFs.find(&temp) != polyFFFFs.end());
#ifndef NOTIME
  idenTime += getTime();
#endif
  if (found)
    return;
  PolyFFFF *polyFFFF = new PolyFFFF(fe1, fe2, fv1, fv2);
  polyFFFFs.insert(polyFFFF);
  //polyFFFF->incRef();

  SumFace **f = temp.f;

#ifndef NOTIME
  idenTime -= getTime();
#endif

  PolyEF *polyEF = fv1->related(fv2);
  SumFace *fe[] = { fe1, fe2 };
  vector<IntEF*> v;
  for (int i = 0; i < 2; i++) {
    const vector<Edge*> &edges = fe[i]->getOuterLoop();
    for (int j = 0; j < edges.size(); j++) {
      SumEdge *edge = dynamic_cast<SumEdge*>(edges[j]);
      // identity
      if (polyEF != 0 && polyEF->contains(edge)) {
#ifndef NOTIME
	idenTime += getTime();
#endif
	return;
      }
#ifdef NOTIME
      idenTime -= getTime();
#endif
      IntEF *vert0 = IntEF::get(edge, fe[1-i]);
      if (vert0 != 0)
        v.push_back(vert0);
      IntEF *vert1 = IntEF::get(edge->getTwin(), fe[1-i]);
      if (vert1 != 0)
        v.push_back(vert1);
#ifdef NOTIME
      idenTime += getTime();
#endif
    }
  }

#ifdef NOTIME
  idenTime += getTime();
#endif

#ifdef NOTIME
  idenTime -= getTime();
#endif

  Vert *commonV = 0;
  const vector<Edge*> &edges = fe[0]->getOuterLoop();
  for (int i = 0; i < edges.size(); i++)
    if (fe[1]->contains(edges[i]->tail)) {
      assert(commonV == 0);
      commonV = edges[i]->tail;
    }

#ifndef NOTIME
  idenTime += getTime();
#endif

  if (polyEF != 0) {
    IntEdge *edge = IntEdge::get(fe1, fe2);
    if (!edge->smallestTwin())
      edge = edge->getTwin();

    vector<PTR<Angle>> roots = polyEF->getAnglesSorted();
    for (int i = 0; i < roots.size(); i++) {
      Angle *root = roots[i];
      Event temp(this, root);

      if (!edge->life.contains(&temp)) {
        delete root;
        continue;
      }

      if (!fv1->life.contains(&temp)) {
        delete root;
        continue;
      }

      if (!fv2->life.contains(&temp)) {
        delete root;
        continue;
      }

      Vert *tail = 0, *head = 0;
      int count = 0;
      for (int k = 0; k < v.size(); k++)
        if (v[k]->life.contains(root)) {
          count++;
          if (edge->face == v[k]->f) {
            if (edge->twin->face->contains(v[k]->e))
              tail = v[k];
            else
              head = v[k];
          }
          else if (edge->twin->face == v[k]->f) {
            if (edge->face->contains(v[k]->e))
              head = v[k];
            else
              tail = v[k];
          }
        }
      assert(tail != 0 || head != 0);
      if (tail == 0) {
	tail = commonV;
	count++;
      }
      else if (head == 0) {
	head = commonV;
	count++;
      }
      assert(count == 2 && tail != 0 && head != 0);

      IntFFF *vv[2];
      bool valid = true;
      SumFace *fv[] = { fv1, fv2 };
      for (int j = 0; j < 2; j++) {
        SumFace *face = fv[j];
        int sTail = face->side(&temp, tail);
        int sHead = face->side(&temp, head);
        if (sTail * sHead > 0) {
          valid = false;
          break;
        }
        if (sTail == -1 &&
            face->intersectsInside(&temp, tail, edge))
          vv[j] = IntFFF::make(edge->getTwin()->getFace(), edge->getFace(),
                               face);
        else if (sHead == -1 &&
                 face->intersectsInside(&temp, head, edge->getTwin()))
          vv[j] = IntFFF::make(edge->getFace(), edge->getTwin()->getFace(),
                               face);
        else {
          valid = false;
          break;
        }
      }
      if (!valid) {
        delete root;
        continue;
      }

      Event *event = 0;
      set<Event*, Compare>::iterator it = eventList.find(&temp);
      if (it != eventList.end()) {
        delete root;
        event = *it;
      }
      else {
        event = new Event(this, root);
        eventList.insert(event);
      }
      if (event == (C4D::Event *) 0x1ca1b40)
	cout << "swap this is it" << endl;
      event->swap.insert(pair<Obj*,Obj*>(vv[0], vv[1]));
    }
    return;
  }

#ifndef NOTIME
  idenTime -= getTime();
#endif

  for (int i = 0; i < 4; i++) {
    int j = (i+1)%4;
    int k = (j+1)%4;
    SumVV *svv = f[i]->commonVert(f[j]);
    if (svv != 0 && f[k]->contains(svv)) {
#ifndef NOTIME
      idenTime += getTime();
#endif
      PolyNNN(f[i], f[j], f[k]).addEvents();
      return;
    }
  }
      
#ifndef NOTIME
  idenTime += getTime();
#endif

#ifdef BLEEN
  {
    const vector<Edge*> &edges = fv1->getOuterLoop();
    for (int i = 0; i < edges.size(); i++)
      if (fv2->contains(edges[i]->tail) &&
	  (fe1->contains(edges[i]->tail) || fe2->contains(edges[i]->tail)))
	return;
  }
#endif

#ifndef NOTIME
  idenTime -= getTime();
#endif
  bool neighbors = fv1->neighborOf(fv2);
#ifndef NOTIME
  idenTime += getTime();
#endif

  if (neighbors)
    return;

  IntEdge *e[3][4];
  for (int i = 0; i < 3; i++)
    for (int j = i+1; j < 4; j++)
      if ((e[i][j] = IntEdge::get(f[i], f[j])) == 0)
        return;

  BigPoly *poly = polyFFFF;

  {
    EdgeAB *ab, *cd[3];
    if (poly->isEpEEE(ab, cd)) {
      cout << "FFFF isEpEEE " << count << endl;
      set<BigPoly*, BigPoly::CompareEpEEE>::iterator
        it = polyEpEEEs.find(poly);
      if (it == polyEpEEEs.end())
        polyEpEEEs.insert(poly);
      else
        poly = *it;
    }
  }

  {
    EdgeAB *ab[3], *cd;
    if (poly->isEEEpE(ab, cd)) {
      cout << "FFFF isEEEpE " << count << endl;
      set<BigPoly*, BigPoly::CompareEEEpE>::iterator
        it = polyEEEpEs.find(poly);
      if (it == polyEEEpEs.end())
        polyEEEpEs.insert(poly);
      else
        poly = *it;
    }
  }

  vector<PTR<Angle>> roots = poly->getAnglesSorted();

  for (int h = 0; h < roots.size(); h++) {
    Angle *root = roots[h];
    Event temp(this, root);
    bool valid = true;

    for (int i = 0; i < 4; i++)
      if (!f[i]->life.contains(root))
        valid = false;
    if (!valid) {
      delete root;
      continue;
    }

    for (int i = 0; i < 3; i++)
      for (int j = i+1; j < 4; j++)
        if (!e[i][j]->life.contains(root))
          valid = false;
    if (!valid) {
      delete root;
      continue;
    }

    IntEdge *e12 = IntEdge::get(fe1, fe2);
    SumFace *fv[] = { fv1, fv2 };
    for (int j = 0; j < 2; j++) {
      Vert *tail = 0, *head = 0;
      int count = 0;
      for (int k = 0; k < v.size(); k++)
        if (v[k]->life.contains(root)) {
          count++;
          if (fv[j]->side(&temp, v[k]) < 0)
            tail = v[k];
          else
            head = v[k];
        }
      if (commonV != 0) {
	count++;
	if (fv[j]->side(&temp, commonV) < 0)
	  tail = commonV;
	else
	  head = commonV;
      }
      assert(count == 2);
      if (tail == 0 || head == 0 ||
          !fv[j]->intersectsInside(&temp, tail, e12)) {
        valid = false;
        break;
      }
    }
    if (!valid) {
      delete root;
      continue;
    }

    Event *event = 0;
    if (eventList.find(&temp) != eventList.end()) {
      event = *eventList.find(&temp);
      delete root;
    }
    else {
      event = new Event(this, root);
      eventList.insert(event);
    }
    for (int i = 0; i < 4; i++)
      event->appear.insert(f[i]);
    event->vanish.insert(f[0]);
  }
}

// SumEdge sumE collides with this face at EF event, forming IntEdge intE.
void C4D::SumFace::handleEF (Event *event, IntEdge *intE, SumEdge *sumE) {
  RootAngle *root = dynamic_cast<RootAngle*>(event->getAngle());
  PolyEF *poly = dynamic_cast<PolyEF*>(root->getPoly());

  vector<IntEdge*> newEdges;

  for (int i = 0; i < 2; i++) {
    const vector<SumFace*> &faces = sumE->getFaces();
    for (int k = 0; k < faces.size(); k++) {
      SumFace *face = faces[k];
      IntEdge *edge = findLiveEdge(face);
      if (edge == 0)
        continue;
      newEdges.push_back(edge);
      /*
        if (edge->smallestTwin())
        edge->appearSubs(event);
        else
        edge->twin->appearSubs(event);
      */
    }
    sumE = sumE->getTwin();
  }

  Vert *th[] = { intE->tail, intE->twin->tail };
  for (int i = 0; i < 2; i++) {
    if (th[i]->getType() == IEF) {
      IntEF *v = dynamic_cast<IntEF*>(th[i]);
      if (v->f != intE->twin->face)
	continue;
      if (poly->contains(v->e->tail)) {
        th[i] = v->e->tail;
	assert(!poly->contains(v->e->twin->tail));
      }
      else if (poly->contains(v->e->twin->tail))
        th[i] = v->e->twin->tail;
      else
	assert(0);
    }
  }

  for (set<IntEdge*, Compare>::iterator it = liveEdges.begin();
       it != liveEdges.end(); ++it) {
    IntEdge *e = *it;
    if (find(newEdges.begin(), newEdges.end(), e) != newEdges.end())
      continue;
    int tailSign = e->getFace()->side(event, th[0]);
    int headSign = e->getFace()->side(event, th[1]);
    if (tailSign * headSign != -1)
      continue;
    tailSign = intE->getFace()->side(event, e->tail);
    headSign = intE->getFace()->side(event, e->twin->tail);
    if (tailSign * headSign != -1)
      continue;

    tailSign = e->getFace()->side(event, sumE->tail);
    headSign = e->getFace()->side(event, sumE->twin->tail);
    SumEdge *se = tailSign == -1 ? sumE : sumE->getTwin();
    assert(tailSign * headSign == -1);
    assert(e->getFace()->intersectsInside(event, se));

    SumFace *f[] = { this, e->getFace() };
    sumE->handleEFF(event, f, true);
  }
}

void C4D::SumVV::handleVF (Event *event, SumFace *fSweep, SumEdge *eHit) {
  static int count;
  ++count;
  RootAngle *root = dynamic_cast<RootAngle*>(event->getAngle());
  PolyEF *poly = root->getPoly()->getType() == AnglePoly::EF2 ?
    dynamic_cast<PolyEF*>(root->getPoly()) : 0;
  if (poly != 0)
    assert(poly->contains(fSweep));
           
  // Edges swept into by sweep plane.
  vector<IntEdge*> edges;
  int sweepSign = 0;

  for (int i = 0; i < intEdges.size(); i++) {
    IntEdge *ie = intEdges[i];
    
    // IntEdge parallel to sweep face.
    if (poly != 0 &&
        (poly->contains(ie->getFace()) ||
         poly->contains(ie->getTwin()->getFace())))
      continue;

    // Face appears but is sweeping OFF ie instead of ON it.
    if (poly != 0 &&
        event->findAppear(fSweep) != -1 &&
        DotEF(event, eHit, fSweep) != DotEF(event, ie, fSweep))
      continue;
    
    IntEdge *ie1 = fSweep->findLiveEdge(ie->getFace());
    IntEdge *ie2 = fSweep->findLiveEdge(ie->getTwin()->getFace());
    
    // One of the other IntEdges no longer exists.
    if (ie1 == 0 || ie2 == 0)
      continue;

    // If FFF exists, remove it.  Done.
    if (ie->numVerts() != 0) {
      IntFFF *v0 = dynamic_cast<IntFFF*>(ie->getVert(0));
      if (v0->hasFace(fSweep)) {
        ie->handleRemove(event, v0);
        ie1->handleRemove(event, v0);
        ie2->handleRemove(event, v0);
        continue;
      }
    }

    // Determine if sweep direction matches normal to face.
    // Should be the same for all intEdges.
    if (sweepSign == 0)
      sweepSign = DotEF(event, ie, fSweep);
    else
      assert(DotEF(event, ie, fSweep) == sweepSign);

    IntFFF *v = (sweepSign < 0) ?
      IntFFF::make(fSweep, ie->getFace(), ie->getTwin()->getFace()) :
      IntFFF::make(fSweep, ie->getTwin()->getFace(), ie->getFace());

    // Add new vertex to IntEdge.
    ie->handleTail(event, v);

    // Save IntEdge to handle other two edges meeting at vertex.
    edges.push_back(ie);
  }
  
  // Collect faces other than the sweep face meeting at a vertex.
  vector<SumFace*> faces;
  for (int i = 0; i < edges.size(); i++) {
    IntFFF *v = dynamic_cast<IntFFF*>(edges[i]->getVert(0));
    assert(v->hasFace(fSweep));
    for (int j = 0; j < 3; j++) {
      SumFace *face = v->faces[j];
      if (face == fSweep)
        continue;
      int k = 0;
      while (k < faces.size() && faces[k] != face) k++;
      if (k == faces.size())
        faces.push_back(face);
    }
  }
  
  // Sort vertices along edge defined by each face.
  for (int i = 0; i < faces.size(); i++) {
    SumFace *face = faces[i];
    Edge *edge = (sweepSign < 0) ?
      face->findLiveEdge(fSweep) : fSweep->findLiveEdge(face);
    vector<IntEdge*> ordered;
    for (int j = 0; j < edges.size(); j++) {
      IntEdge *ie = edges[j];
      IntFFF *v = dynamic_cast<IntFFF*>(ie->getVert(0));
      if (!v->hasFace(face))
        continue;
      ordered.push_back(ie);
      int k = ordered.size()-1;
      while (k > 0 && CircEEF(event, ie, ordered[k-1], face) < 0) {
        ordered[k] = ordered[k-1];
        k--;
      }
      ordered[k] = ie;
    }

    // Add vertices to edge.
    vector<Vert*> verts;
    for (int k = 0; k < ordered.size(); k++)
      verts.push_back(ordered[k]->getVert(0));
    edge->handleReplace(event, -2, -3, verts);
  }
}

void C4D::handleFFFF (Event *event) {
  assert(0);
#ifdef SWEEP_ALL_FACES
  SumFace *f[4];
  for (int i = 0; i < 4; i++)
    f[i] = dynamic_cast<SumFace*>(event->appear[i]);
  assert(event->vanish[0] == f[0]);
  event->vanish.erase(event->vanish.begin());
  event->appear.erase(event->appear.begin());
  event->appear.erase(event->appear.begin());
  event->appear.erase(event->appear.begin());
  event->appear.erase(event->appear.begin());

  int k, l;
  for (int i = 0; i < 3; i++) {
    for (int j = i+1; j < 4; j++) {
      for (k = 0; k == i || k == j; k++);
      for (l = 0; l == i || l == j || l == k; l++);
      IntEdge *edge = f[i]->findLiveEdge(f[j]);
      if (!edge->smallestTwin())
	edge = edge->getTwin();
      IntFFF *vk = IntFFF::get(f[i], f[j], f[k]);
      if (vk == 0 || edge->findVert(vk) == -1)
        vk = IntFFF::get(f[j], f[i], f[k]);
      IntFFF *vl = IntFFF::get(f[i], f[j], f[l]);
      if (vl == 0 || edge->findVert(vl) == -1)
        vl = IntFFF::get(f[j], f[i], f[l]);
      int index[] = { edge->findVert(vk), edge->findVert(vl) };
      if (index[0] == -1 || index[1] == -1) {
	RootAngle *root = dynamic_cast<RootAngle*>(event->getAngle());
	BigPoly *poly = dynamic_cast<BigPoly*>(root->getPoly());
	EdgeAB *ab, *cd[3];
	assert(poly->isEpEEE(ab, cd) || poly->isEEEpE(cd, ab));
	continue;
      }
      if (index[0] < index[1]) {
        assert(index[0]+1 == index[1]);
        edge->handleSwap(event, index[0]);
      }
      else {
        assert(index[1]+1 == index[0]);
        edge->handleSwap(event, index[1]);
      }
    }
  }
#endif
}

void C4D::SubEdge::addEvents () {
  if (e->getType() != IFF) {
    assert(e->getType() == SVE || e->getType() == SEV);
    SumEdge *edge = dynamic_cast<SumEdge*>(e);
    IntEF *tail = dynamic_cast<IntEF*>(getTail());
    IntEF *head = dynamic_cast<IntEF*>(getHead());
    getC4D()->addEvents(edge, tail->f, head->f);
  }
  else {
    IntEdge *edge = dynamic_cast<IntEdge*>(e);
    SumFace *f0 = edge->getFace();
    SumFace *f1 = edge->getTwin()->getFace();
    IntFFF *tail = dynamic_cast<IntFFF*>(getTail());
    IntFFF *head = dynamic_cast<IntFFF*>(getHead());
    SumFace *tailf = tail->otherFace(f0, f1);
    SumFace *headf = head->otherFace(f0, f1);
    getC4D()->addEvents(f0, f1, tailf, headf);
  }
}

void printLife (C4D::EventSet &life) {
  cout << "life " << life.events.size()/2;
  for (int i = 0; i < life.events.size(); i++) {
    PV2<Parameter> u = life.events[i]->getAngle()->get<Parameter>();
    double x = u.x.mid();
    double y = u.y.mid();
    // cout << " " << u.x.mid() << " " << u.y.mid();
    disable();
    cout << " " << x << " " << y;
    enable();
  }
  cout << endl;
}

bool printFaceLife = false;

vector<C4D::Face *> C4D::getSumFaces () {
  static int count;
  vector<Face *> sumFaces;
  for (int i = 0; i < verts[0].size(); i++) {
    VertXYZ *v = verts[0][i];
    if (v->cone.size() == 0)
      continue;
    for (int j = 0; j < faces[!0].size(); j++) {
      FaceABC *f = faces[!0][j];
      count++;
      SumFace *sum = v->getSum(f);
      if (sum != 0) {
        sumFaces.push_back(sum);
	if (printFaceLife) {
	  cout << "face "
	       << "1 " << i << " "
	       << "3 " 
	       << f->getVert(0)->id << " "
	       << f->getVert(1)->id << " "
	       << f->getVert(2)->id << endl;
	  printLife(sum->life);
	}
     }
    }
  }
  for (int i = 0; i < verts[1].size(); i++) {
    VertXYZ *v = verts[1][i];
    if (v->cone.size() == 0)
      continue;
    //cout << "testing robot vertex " << i << endl;
    for (int j = 0; j < faces[!1].size(); j++) {
      FaceABC *f = faces[!1][j];
      SumFace *sum = v->getSum(f);
      if (sum != 0) {
        sumFaces.push_back(sum);
	if (printFaceLife) {
	  cout << "face "
	       << "3 " 
	       << f->getVert(0)->id << " "
	       << f->getVert(1)->id << " "
	       << f->getVert(2)->id << " "
	       << "1 " << i << endl;
	  printLife(sum->life);
	}
      }
    }
  }
  for (int i = 0; i < edges[0].size(); i++) {
    EdgeAB *a = edges[0][i];
    if (!a->isConvex())
      continue;
    if (a->getTail()->id > a->getHead()->id)
      continue;
    for (int j = 0; j < edges[!0].size(); j++) {
      EdgeAB *b = edges[!0][j];
      if (!b->isConvex())
        continue;
      SumFace *sum = a->getSum(b);
      if (sum != 0) {
        sumFaces.push_back(sum);
	if (printFaceLife) {
	  cout << "face "
	       << "2 " << a->getTail()->id << " " << a->getHead()->id << " "
	       << "2 " << b->getTail()->id << " " << b->getHead()->id << endl;
	  printLife(sum->life);
	}
      }
    }
  }
  return sumFaces;
}

vector<C4D::Face *> C4D::getSumFaces (Event *event) {
  vector<Face *> sumFaces;
  for (int i = 0; i < verts[0].size(); i++) {
    VertXYZ *v = verts[0][i];
    if (v->cone.size() == 0)
      continue;
    for (int j = 0; j < faces[!0].size(); j++) {
      FaceABC *f = faces[!0][j];
      if (v->isCompatible(event, f)) {
        SumVF *sumVF = SumVF::make(v, f);
        sumFaces.push_back(sumVF);
	if (!sumVF->life.isNull())
	  assert(sumVF->life.contains(event));
        // string s = sumFaces.back()->toString();
        // cout << s << endl;
      }
      else {
        SumFace *sum = SumVF::get(v, f);
        if (sum != 0)
          assert(!sum->life.contains(event));
      }
    }
  }
  for (int i = 0; i < verts[1].size(); i++) {
    VertXYZ *v = verts[1][i];
    if (v->cone.size() == 0)
      continue;
    for (int j = 0; j < faces[!1].size(); j++) {
      FaceABC *f = faces[!1][j];
      if (v->isCompatible(event, f)) {
        sumFaces.push_back(SumFV::make(f, v));
	if (!sumFaces.back()->life.isNull())
	  assert(sumFaces.back()->life.contains(event));
        // string s = sumFaces.back()->toString();
        // cout << s << endl;
      }
      else {
        SumFace *sum = SumFV::get(f, v);
        if (sum != 0)
          assert(!sum->life.contains(event));
      }
    }
  }
  for (int i = 0; i < edges[0].size(); i++) {
    EdgeAB *a = edges[0][i];
    if (!a->isConvex())
      continue;
    if (a->getTail()->id > a->getHead()->id)
      continue;
    for (int j = 0; j < edges[!0].size(); j++) {
      EdgeAB *b = edges[!0][j];
      if (!b->isConvex())
        continue;
      if (a->isCompatible(event, b) &&
          b->isCompatible(event, dynamic_cast<EdgeAB*>(a->getTwin()))) {
        sumFaces.push_back(SumEE::make(a, b));
        
        Face *f = sumFaces.back();
        // assert(a->isOutside(event, f));
        // assert(dynamic_cast<EdgeAB*>(a->getTwin())->isOutside(event, f));
        // assert(b->isOutside(event, f));
        // assert(dynamic_cast<EdgeAB*>(b->getTwin())->isOutside(event, f));
	if (!f->life.isNull())
	  assert(f->life.contains(event));        
        // cout << sumFaces.back()->toString() << endl;
      }
      else {
        SumFace *sum = SumEE::get(a, b);
        if (sum != 0)
          assert(!sum->life.contains(event));
      }
    }
  }
  
  return sumFaces;
}

C4D::EventSet C4D::IntEF::getLife (SumEdge *e, SumFace *f) {
  if (e == (C4D::SumEdge *) 0x68a1d0 &&
      f == (C4D::SumFace *) 0x695b70)
    cout << "this is it IntEF::getLife" << endl;

  EventSet life = e->life.intersect(f->life);
  if (life.isNull())
    return life;

  SumVV *t = e->getTail();
  // PolyVF *tPoly = new PolyVF(t, f);
  // EventSet tNeg = tPoly->getNegative();
  EventSet tNeg = PolyVF::getNegative(t, f);
  life = life.intersect(tNeg);
  if (life.isNull())
    return life;

  SumVV *h = e->getHead();
  // PolyVF *hPoly = new PolyVF(h, f);
  // EventSet hPos = hPoly->getNegative().complement();
  EventSet hPos = PolyVF::getNegative(h, f).complement();
  life = life.intersect(hPos);
  if (life.isNull())
    return life;

  const vector<Edge *> &edges = f->getOuterLoop();
  for (int i = 0; i < edges.size(); i++) {
    SumEdge *side = dynamic_cast<SumEdge*>(edges[i]);
    // PolyEE *poly = new PolyEE(side, e);
    // EventSet eNeg = poly->getNegative();
    EventSet eNeg = PolyEE::getNegative(side, e);
    life = life.intersect(eNeg);
    if (life.isNull())
      return life;
  }
  
  return life;
}

void C4D::IntEdge::checkForIntersections (Event *event) {
  static int count = 0;
  assert(smallestTwin());
  SumFace *sface = getFace();
  SumFace *tface = getTwin()->getFace();

  assert(event->getAngle()->getType() == Angle::IMPLICIT ||
         event->getAngle()->getType() == Angle::EXPLICIT);
  RootAngle *root = dynamic_cast<RootAngle*>(event->getAngle());
  assert(root->getPoly()->getType() == AnglePoly::EF2);
  PolyEF *poly = dynamic_cast<PolyEF*>(root->getPoly());
  Vert *th[] = { tail, twin->tail };
  for (int i = 0; i < 2; i++)
    if (th[i]->getType() == IEF) {
      IntEF *ef = dynamic_cast<IntEF*>(th[i]);
      if (poly->contains(ef->f)) {
        if (poly->contains(ef->e->tail)) {
          assert(!poly->contains(ef->e->twin->tail));
          th[i] = ef->e->tail;
        }
        else if (poly->contains(ef->e->twin->tail)) {
          th[i] = ef->e->twin->tail;
        }
      }
    }

  SumFace *f_ef = 0;
  if (poly->contains(sface))
    f_ef = sface;
  else if (poly->contains(tface))
    f_ef = tface;
  else
    return;

  if (f_ef != sface)
    assert(poly->contains(tface));
  SumFace *f_nef = f_ef != sface ? sface : tface;
  SumEdge *e_ef = 0;
  const vector<Edge*> &edges = f_nef->getOuterLoop();
  for (int i = 0; i < edges.size(); i++)
    if (poly->contains(edges[i]))
      e_ef = dynamic_cast<SumEdge*>(edges[i]);
  if (e_ef == 0)
    return;

  for (set<IntEdge *, Compare>::iterator ie = sface->liveEdges.begin();
       ie != sface->liveEdges.end(); ++ie) {
    if (event->findAppear(*ie) != -1)
      continue;
    SumFace *f = (*ie)->getFace();
    if (f == tface)
      continue;
    if (f->neighborOf(tface))
      continue;

    int sTail = f->side(event, th[0]);
    int sHead = f->side(event, th[1]);
    if (sTail * sHead != -1)
      continue;
    IntEdge *e = sTail == -1 ? this : getTwin();
    bool inside = f->intersectsInside(event, e);
    if (!inside)
      continue;
    IntFFF *v = IntFFF::make(e->getTwin()->getFace(), e->getFace(), f);
    addVert(event, v);

    IntEdge *ff_nef = f->findLiveEdge(f_nef);
    // ??? really should be handleTail or handleHead
    ff_nef->handleInsert(event, v);

    // SumFace *f_nef2 = e_ef->getTwin()->getFace();

    IntEdge *ff_ef = f->findLiveEdge(f_ef);

    vector<SumFace*> faces;
    int n0;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < e_ef->getFaces().size(); j++) {
        SumFace *face = e_ef->getFaces()[j];
        if (f_ef->findLiveEdge(face) != 0)
          faces.push_back(face);
      }
      e_ef = e_ef->getTwin();
      if (i == 0)
        n0 = faces.size();
    }
    assert(faces.size() >= 1);
    
    if (faces.size() == 1) {
      ff_ef->handleInsert(event, v);
      continue;
    }

    int iSmaller = -1;
    for (int i = 1; i < faces.size(); i++)
      if (Compare()(faces[i], f_nef))
        iSmaller = i;
    if (iSmaller != -1)
      continue;

    vector<IntFFF*> vs;
    for (int i = 0; i < faces.size(); i++) {
      IntFFF *v = (e->getFace() == f_ef) == (i >= n0) ?
        IntFFF::make(f_ef, faces[i], f) : IntFFF::make(faces[i], f_ef, f);
      vs.push_back(v);
      int j = i;
      int sign = i > 0 ? -1 : 1;
      while (j > 0 &&
             CircEFF(event, e_ef, faces[i], vs[j-1]->otherFace(f_ef, f)) ==
             sign * DotEF(event, ff_ef, vs[j-1]->otherFace(f_ef, f))) {
        vs[j] = vs[j-1];
	j--;
      }
      vs[j] = v;
    }

    vector<Vert*> vs2;
    for (int i = 0; i < vs.size(); i++)
      vs2.push_back(vs[i]);

    ff_ef->handleReplace(event, 1, 0, vs2);
  }
}

set<C4D::Obj *, C4D::Compare> C4D::getObjects$ () {
  enable();
  set<C4D::Obj *, C4D::Compare> objects = getObjects();
  disable();
  return objects;
}

// Create SumFace and SumEdge objects (and their lifetimes);
set<C4D::Obj *, C4D::Compare> C4D::getSums (vector<Face *> &sumFaces) {
  set<Obj *, Compare> objs;

  for (int i = 0; i < sumFaces.size(); i++) {
    SumFace *f = dynamic_cast<SumFace*>(sumFaces[i]);
    objs.insert(f);
    const vector<Edge *> &fedges = f->getOuterLoop();
    int fn = fedges.size();
    Vert *fverts[fn];
    for (int j = 0; j < fn; j++) {
      fverts[j] = fedges[j]->getTail();

      SumEdge *edge = dynamic_cast<SumEdge*>(fedges[j]);
      edge->life = edge->life.combine(f->life);
      edge->getTwin()->life = edge->life;
      objs.insert(edge);
      objs.insert(edge->getTwin());
    }
  }

  return objs;
}

// Create SumFace, SumEdge, IntEF, IntEdge and their lifetimes.
set<C4D::Obj *, C4D::Compare> C4D::getObjects () {
  static int count;
  vector<Face *> sumFaces = getSumFaces();

  set<Obj *, Compare> objs;
  for (int i = 0; i < sumFaces.size(); i++) {
    SumFace *f = dynamic_cast<SumFace*>(sumFaces[i]);
    objs.insert(f);
    const vector<Edge *> &fedges = f->getOuterLoop();
    int fn = fedges.size();
    Vert *fverts[fn];
    for (int j = 0; j < fn; j++) {
      fverts[j] = fedges[j]->getTail();

      SumEdge *edge = dynamic_cast<SumEdge*>(fedges[j]);
      edge->life = edge->life.combine(f->life);
      edge->getTwin()->life = edge->life;
      objs.insert(edge);
      objs.insert(edge->getTwin());
    }
  }

  for (int i = 0; i < sumFaces.size(); i++) {
    SumFace *f = dynamic_cast<SumFace*>(sumFaces[i]);
    const vector<Edge *> &fedges = f->getOuterLoop();
    int fn = fedges.size();
    Vert *fverts[fn];
    for (int j = 0; j < fn; j++)
      fverts[j] = fedges[j]->getTail();

    for (int j = 0; j < sumVEs.size(); j++) {
      SumEdge *e = sumVEs[j];
      if (f == (C4D::SumFace *) 0x695b70 &&
	  e == (C4D::Edge *) 0x68a1d0)
	cout << "this is it sumVEs" << endl;

      if (e == 0)
        continue;

      bool neighbor = false;
      for (int k = 0; k < fn; k++)
        if (fverts[k] == e->getTail() || fverts[k] == e->getHead())
          neighbor = true;
      if (neighbor)
        continue;

      EventSet efLife = IntEF::getLife(e, f);
      if (efLife.isNull())
        continue;

      IntEF *intEF = IntEF::make(e, f);
      intEF->life = efLife;
      objs.insert(intEF);
    }

    for (int j = 0; j < sumEVs.size(); j++) {
      SumEdge *e = sumEVs[j];
      if (e == 0)
        continue;

      bool neighbor = false;
      for (int k = 0; k < fn; k++)
        if (fverts[k] == e->getTail() || fverts[k] == e->getHead())
          neighbor = true;
      if (neighbor)
        continue;

      EventSet efLife = IntEF::getLife(e, f);
      if (efLife.isNull())
        continue;
      
      IntEF *intEF = IntEF::make(e, f);
      intEF->life = efLife;
      objs.insert(intEF);
    }
  }

  for (int i = 0; i < sumFaces.size(); i++) {
    SumFace *f = dynamic_cast<SumFace*>(sumFaces[i]);
    for (int j = i+1; j < sumFaces.size(); j++) {
      if (j == i) continue;
      SumFace *g = dynamic_cast<SumFace*>(sumFaces[j]);

      if ((f == (C4D::SumFace *) 0x695b70 ||
	   f == (C4D::SumFace *) 0x6e3cc0) &&
	  (g == (C4D::SumFace *) 0x695b70 ||
	   g == (C4D::SumFace *) 0x6e3cc0)) {
	cout << "this is it fg" << endl;
	const vector<Edge*> &fedges = f->getOuterLoop();
	for (int k = 0; k < fedges.size(); k++) {
	  Edge *e = fedges[k];
	  int sTail = g->side(EVENT0, e->tail);
	  int sHead = g->side(EVENT0, e->twin->tail);
	  if (sTail == -1 && sHead == 1 && g->intersectsInside(EVENT0, e))
	    cout << "this is it fg intersection" << endl;
	  e = e->twin;
	  sTail = g->side(EVENT0, e->tail);
	  sHead = g->side(EVENT0, e->twin->tail);
	  if (sTail == -1 && sHead == 1 && g->intersectsInside(EVENT0, e))
	    cout << "this is it fg intersection" << endl;
	}
	const vector<Edge*> &gedges = g->getOuterLoop();
	for (int k = 0; k < gedges.size(); k++) {
	  Edge *e = gedges[k];
	  int sTail = f->side(EVENT0, e->tail);
	  int sHead = f->side(EVENT0, e->twin->tail);
	  if (sTail == -1 && sHead == 1 && f->intersectsInside(EVENT0, e))
	    cout << "this is it fg intersection" << endl;
	  e = e->twin;
	  sTail = f->side(EVENT0, e->tail);
	  sHead = f->side(EVENT0, e->twin->tail);
	  if (sTail == -1 && sHead == 1 && f->intersectsInside(EVENT0, e))
	    cout << "this is it fg intersection" << endl;
	}
      }

      EventSet intLife = f->getIntLife(g);
      if (intLife.isNull()) continue;

      IntEdge *e = IntEdge::make(f, g);
      if (e == (IntEdge *) 0x298fce0)
	cout << "this is it IntEdge" << endl;
      e->life = intLife;
      e->getTwin()->life = intLife;
      objs.insert(e);
      objs.insert(e->getTwin());
    }
  }

  return objs;
}

// Create IntEF and IntEdge objects (and their lifetimes) that
// intersect SumFace f and add to objs.
void  C4D::getInts (vector<Face *> &sumFaces, SumFace *f, 
                    set<Obj *, Compare> &objs) {
  {
    const vector<Edge *> &fedges = f->getOuterLoop();
    int fn = fedges.size();
    Vert *fverts[fn];
    for (int j = 0; j < fn; j++)
      fverts[j] = fedges[j]->getTail();

    for (int j = 0; j < sumVEs.size(); j++) {
      SumEdge *e = sumVEs[j];
      if (e == 0)
        continue;

      bool neighbor = false;
      for (int k = 0; k < fn; k++)
        if (fverts[k] == e->getTail() || fverts[k] == e->getHead())
          neighbor = true;
      if (neighbor)
        continue;

      EventSet efLife = IntEF::getLife(e, f);
      if (efLife.isNull())
        continue;

      IntEF *intEF = IntEF::make(e, f);
      intEF->life = efLife;
      objs.insert(intEF);
    }

    for (int j = 0; j < sumEVs.size(); j++) {
      SumEdge *e = sumEVs[j];
      if (e == 0)
        continue;

      bool neighbor = false;
      for (int k = 0; k < fn; k++)
        if (fverts[k] == e->getTail() || fverts[k] == e->getHead())
          neighbor = true;
      if (neighbor)
        continue;

      EventSet efLife = IntEF::getLife(e, f);
      if (efLife.isNull())
        continue;
      
      IntEF *intEF = IntEF::make(e, f);
      intEF->life = efLife;
      objs.insert(intEF);
    }
  }

  for (int i = 0; i < sumFaces.size(); i++) {
    SumFace *g = dynamic_cast<SumFace*>(sumFaces[i]);
    if (g == f)
      continue;
    const vector<Edge *> &gedges = g->getOuterLoop();
    int gn = gedges.size();
    Vert *gverts[gn];
    for (int j = 0; j < gn; j++)
      gverts[j] = gedges[j]->getTail();

    for (int j = 0; j < f->numEdges(); j++) {
      SumEdge *e = f->getEdge(j);

      bool neighbor = false;
      for (int k = 0; k < gn; k++)
        if (gverts[k] == e->getTail() || gverts[k] == e->getHead())
          neighbor = true;
      if (neighbor)
        continue;
      
      {
        EventSet egLife = IntEF::getLife(e, g);
        if (!egLife.isNull()) {
          IntEF *intEF = IntEF::make(e, g);
          intEF->life = egLife;
          objs.insert(intEF);
        }
      }
      {
        EventSet egLife = IntEF::getLife(e->getTwin(), g);
        if (!egLife.isNull()) {
          IntEF *intEF = IntEF::make(e->getTwin(), g);
          intEF->life = egLife;
          objs.insert(intEF);
        }
      }
    }
  }

  {
    for (int j = 0; j < sumFaces.size(); j++) {
      SumFace *g = dynamic_cast<SumFace*>(sumFaces[j]);
      if (g == f)
        continue;

      EventSet intLife = f->getIntLife(g);
      if (intLife.isNull()) continue;

      IntEdge *e = IntEdge::make(f, g);
      e->life = intLife;
      e->getTwin()->life = intLife;
      objs.insert(e);
      objs.insert(e->getTwin());
    }
  }
}

set<C4D::Obj *, C4D::Compare> C4D::getObjects (Event *event) {
  SumFace *bad0 = (SumFace*)0x2a06090;
  SumFace *bad1 = (SumFace*)0x2a08f20;
  SumFace *bad2 = (SumFace*)0x2a090c0;

  enable();

  for (int i = 0; i < sumVEs.size(); i++)
    if (sumVEs[i] != 0)
      sumVEs[i]->clearFaces();

  for (int i = 0; i < sumEVs.size(); i++)
    if (sumEVs[i] != 0)
      sumEVs[i]->clearFaces();

  vector<Face *> sumFaces = getSumFaces(event);
  set<Obj *, Compare> objs;
  for (vector<Face *>::iterator it = sumFaces.begin();
       it != sumFaces.end(); ++it) {
    SumFace *face = dynamic_cast<SumFace*>(*it);
    objs.insert(face);
    for (set<IntEdge *, Compare>::iterator ie = face->liveEdges.begin();
         ie != face->liveEdges.end(); ++ie)
      (*ie)->verts.clear();
    face->liveEdges.clear();
    const vector<Edge *> &edges = face->getOuterLoop();
    for (int i = 0; i < edges.size(); i++) {
      edges[i]->verts.clear();
      edges[i]->twin->verts.clear();
      dynamic_cast<SumEdge*>(edges[i])->addFace(dynamic_cast<SumFace*>(*it));
    }
  }

  for (int i = 0; i < sumFaces.size(); i++) {
    SumFace *f = dynamic_cast<SumFace*>(sumFaces[i]);
    const vector<Edge *> &fedges = f->getOuterLoop();
    int fn = fedges.size();

    for (int j = 0; j < fn; j++) {
      objs.insert(fedges[j]);
      objs.insert(fedges[j]->twin);
    }

    for (int j = i+1; j < sumFaces.size(); j++) {
      SumFace *g = dynamic_cast<SumFace*>(sumFaces[j]);
      const vector<Edge *> &gedges = g->getOuterLoop();
      int gn = gedges.size();

      //if ((f == (SumFace*)0x100308f10 || f == (SumFace*)0x10033e220) &&
      //          (g == (SumFace*)0x100308f10 || g == (SumFace*)0x10033e220))
      // cout << "this is it" << endl;

#ifdef BLEEN
      for (int k = 0; k < fn; k++)
        cout << g->intersects(event, fedges[k]) << " "
             << g->intersects(event, fedges[k]->twin) << " "
             << (fedges[k]->tail->getP(event) - f->getP(event)).dot(f->getN(event)).mid() << endl;

      for (int k = 0; k < gn; k++)
        cout << f->intersects(event, gedges[k]) << " "
             << f->intersects(event, gedges[k]->twin) << " "
             << (gedges[k]->tail->getP(event) - g->getP(event)).dot(g->getN(event)).mid() << endl;
#endif

      bool neighbors = false;
      for (int k = 0; k < fn; k++)
        if (dynamic_cast<SumEdge*>(fedges[k])->hasFace(g))
          neighbors = true;
      if (neighbors)
        continue;

      int fnz = 0, fnn = 0, fnp = 0;
      int fside[4];
      for (int k = 0; k < fn; k++) {
        fside[k] = g->side(event, fedges[k]->tail);
        if (fside[k] < 0)
          fnn++;
        else if (fside[k] > 0)
          fnp++;
        else
          fnz++;
      }

#ifdef BLEEN
      if (fnz == 2) {
        PV3<N> vz = gedges[1]->get<N>(event);
        PV3<N> vx = g->getN(event);
        PV3<N> vy = vz.cross(vx);

        PV3<N> p2[2] = { verts[0][2]->getP(event), verts[1][2]->getP(event) };

        for (int fm = 0; fm < 2; fm++) {
          for (int i = 0; i < 4; i++) {
            PV3<N> p = verts[fm][i]->getP(event) - p2[fm];
            cout << i << " "
                 << vx.dot(p).mid() << " "
                 << vy.dot(p).mid() << " "
                 << vz.dot(p).mid() << endl;
          }
        }
      }
#endif

      assert(fnz < 2);
      if (fnn == 0 || fnp == 0)
        continue;

      int gnz = 0, gnn = 0, gnp = 0;
      int gside[4];
      for (int k = 0; k < gn; k++) {
        gside[k] = f->side(event, gedges[k]->tail);
        if (gside[k] < 0)
          gnn++;
        else if (gside[k] > 0)
          gnp++;
        else
          gnz++;
      }

      assert(gnz < 2);
      if (gnn == 0 || gnp == 0)
        continue;

      IntEdge *edge = 0;

      for (int k = 0; k < fn; k++)
        if (fside[k] < 0 && fside[(k+1)%fn] > 0 &&
            g->intersectsInside(event, fedges[k])) {
          if (edge == 0)
            edge = IntEdge::make(f, g);
          assert(edge->twin->tail == 0);
          edge->twin->tail = IntEF::make(dynamic_cast<SumEdge*>(fedges[k]), g);
          objs.insert(edge->twin->tail);
          fedges[k]->addVert(event, edge->twin->tail);
          if (edge != 0 && edge->tail != 0 && edge->twin->tail != 0)
            break;
        }
        else if (fside[k] > 0 && fside[(k+1)%fn] < 0 &&
                 g->intersectsInside(event, fedges[k]->twin)) {
          if (edge == 0)
            edge = IntEdge::make(f, g);
          assert(edge->tail == 0);
          edge->tail = IntEF::make(dynamic_cast<SumEdge*>(fedges[k]->twin), g);
          objs.insert(edge->tail);
          fedges[k]->twin->addVert(event, edge->tail);
          if (edge != 0 && edge->tail != 0 && edge->twin->tail != 0)
            break;
        }

      if (edge != 0 && fnz > 0)
        for (int k = 0; k < fn; k++)
          if (fside[k] == 0) {
            if (fside[(k+1)%fn] > 0) {
              assert(edge->twin->tail == 0);
              edge->twin->tail = fedges[k]->tail;
            }
            else {
              assert(edge->tail == 0);
              edge->tail = fedges[k]->tail;
            }
          }
      
      if (edge != 0 && edge->tail != 0 && edge->twin->tail != 0) {
        f->liveEdges.insert(edge->getTwin());
        g->liveEdges.insert(edge);
	if (edge == (IntEdge *) 0x2a0cfd0 ||
	    edge->getTwin() == (IntEdge *) 0x2a0cfd0)
	  cout << "this is it IntEdge3" << endl;
        objs.insert(edge);
        objs.insert(edge->getTwin());
        continue;
      }
        
      for (int k = 0; k < gn; k++)
        if (gside[k] < 0 && gside[(k+1)%gn] > 0 &&
            f->intersectsInside(event, gedges[k])) {
          if (edge == 0)
            edge = IntEdge::make(f, g);
          assert(edge->tail == 0);
          edge->tail = IntEF::make(dynamic_cast<SumEdge*>(gedges[k]), f);
          objs.insert(edge->tail);
          gedges[k]->addVert(event, edge->tail);
          if (edge != 0 && edge->tail != 0 && edge->twin->tail != 0)
            break;
        }
        else if (gside[k] > 0 && gside[(k+1)%gn] < 0 &&
                 f->intersectsInside(event, gedges[k]->twin)) {
          if (edge == 0)
            edge = IntEdge::make(f, g);
          assert(edge->twin->tail == 0);
          edge->twin->tail = IntEF::make(dynamic_cast<SumEdge*>(gedges[k]->twin), f);
          objs.insert(edge->twin->tail);
          gedges[k]->twin->addVert(event, edge->twin->tail);
          if (edge != 0 && edge->tail != 0 && edge->twin->tail != 0)
            break;
        }

      if (edge == 0)
        continue;

      if (edge != 0 && edge->tail != 0 && edge->twin->tail != 0) {
        f->liveEdges.insert(edge->getTwin());
        g->liveEdges.insert(edge);
	if (edge == (IntEdge *) 0x2a0cfd0 ||
	    edge->getTwin() == (IntEdge *) 0x2a0cfd0)
	  cout << "this is it IntEdge2" << endl;
        objs.insert(edge);
        objs.insert(edge->getTwin());
        continue;
      }
        
      assert(gnz > 0);
      for (int k = 0; k < gn; k++)
        if (gside[k] == 0) {
          if (gside[(k+1)%gn] > 0) {
            assert(edge->tail == 0);
            edge->tail = gedges[k]->tail;
          }
          else {
            assert(edge->twin->tail == 0);
            edge->twin->tail = gedges[k]->tail;
          }
        }

      assert(edge->tail != 0 && edge->twin->tail != 0);
      f->liveEdges.insert(edge->getTwin());
      g->liveEdges.insert(edge);
      if (edge == (IntEdge *) 0x2a0cfd0 ||
	  edge->getTwin() == (IntEdge *) 0x2a0cfd0)
	cout << "this is it IntEdge1" << endl;
      objs.insert(edge);
      objs.insert(edge->getTwin());
    }
  } 

  for (vector<Face *>::iterator jt = sumFaces.begin();
       jt != sumFaces.end(); ++jt) {
    const vector<Edge *> &edges = (*jt)->getOuterLoop();
    for (int i = 0; i < edges.size(); i++) {
      // I think I forgot to put IntEF's into verts....
      // assert(edges[i]->verts.size() == 0);
      if (edges[i]->smallestTwin())
	edges[i]->makeSubEdges(objs);
      else if (dynamic_cast<SumEdge*>(edges[i]->twin)->getFaces().size() == 0)
	edges[i]->twin->makeSubEdges(objs);
    }
  }

  for (vector<Face *>::iterator jt = sumFaces.begin();
       jt != sumFaces.end(); ++jt) {
    SumFace *face = dynamic_cast<SumFace*>(*jt);
    set<IntEdge *, Compare> &edges = face->liveEdges;
    for (set<IntEdge *, Compare>::iterator ie = edges.begin();
         ie != edges.end(); ++ie) {
      Edge *ab = *ie;
      Vert *a = ab->getTail();
      Vert *b = ab->getHead();
      SumFace *f = dynamic_cast<SumFace*>(ab->face);
      for (set<IntEdge *, Compare>::iterator je = ie;
           ++je != edges.end();) {
        Edge *cd = *je;
        Vert *c = cd->getTail();
        Vert *d = cd->getHead();
        SumFace *g = dynamic_cast<SumFace*>(cd->face);

	if ((face == bad0 || face == bad1 || face == bad2) &&
	    (f == bad0 || f == bad1 || f == bad2) &&
	    (g == bad0 || g == bad1 || g == bad2))
	  cout << "this is it bad" << endl;

        int gSa = g->side(event, a);
        int gSb = g->side(event, b);
        if (gSa * gSb != -1)
          continue;

        int fSc = f->side(event, c);
        int fSd = f->side(event, d);
        if (fSc * fSd != -1)
          continue;

        assert(gSa * fSc < 0);
        IntFFF *v =
          gSa < 0 ? IntFFF::make(face, f, g) : IntFFF::make(face, g, f);
        objs.insert(v);
        ab->addVert(event, v);
        cd->addVert(event, v);

        IntEdge *fg = f->findLiveEdge(g);
        fg->addVert(event, v);
      }
    }
  }

  for (vector<Face *>::iterator jt = sumFaces.begin();
       jt != sumFaces.end(); ++jt) {
    SumFace *face = dynamic_cast<SumFace*>(*jt);
    set<IntEdge *, Compare> &edges = face->liveEdges;
    for (set<IntEdge *, Compare>::iterator ie = edges.begin();
         ie != edges.end(); ++ie) {
      if (*ie == (IntEdge *) 0x2a0cfd0)
	cout << "this is it IntEdge" << endl;
      if ((*ie)->smallestTwin())
	(*ie)->makeSubEdges(objs);
    }
  }

  disable();
  return objs;
}

C4D::Event *global_event;

void C4D::sweep () {
  C4D c4d(this);
  PTR<Event> EVENT0PTR = EVENT0;

  // Create SumFace, SumEdge, IntEF, IntEdge and their lifetimes.
  set<Obj*, Compare> objs = getObjects();

  // Put Events of each Obj into eventList.
  // Inform events of who appears and vanishes.
  // Eliminate duplicate Events.
  // Put initial objects into sweepSet.
  for (set<Obj*, Compare>::iterator it = objs.begin();
       it != objs.end(); ++it) {
    (*it)->informEvents();
    if (*it == (Obj *) 0x69f0f0)
      cout << "this is it Obj" << endl;
    if ((*it)->contains(EVENT0)) {
      sweepSets[(*it)->getType()].insert(*it);
      if (*it == (Obj *) 0x2990ce0)
	cout << "this is it" << endl;
    }
  }

  // At angle 0.  Set faces of SumEdges. Add IntEFs to SumEdges. Add
  // IntEdges to SumFaces.
  for (int ei = ERROR; ei <= SUBF; ei++)
    for (unordered_set<Obj*, Equal, Equal>::iterator it = sweepSets[ei].begin();
         it != sweepSets[ei].end(); ++it) {
      switch ((*it)->getType()) {
      case SVE:
      case SEV:
        // Prepare edges to receive vertices.
        assert(dynamic_cast<Edge*>(*it)->verts.size() == 0);
        dynamic_cast<Edge*>(*it)->verts.clear();
        break;
      case SVF:
      case SFV:
      case SEE:
        {
          SumFace *f = dynamic_cast<SumFace*>(*it);
          // Inform SumEdges whose face they border.
          f->initEdges();
          // Prepare faces to receive IntEdges
          assert(f->liveEdges.size() == 0);
          f->liveEdges.clear();
        }
        break;
      case IEF:
        {
          // Add IntEF to SumEdge.
          IntEF *p = dynamic_cast<IntEF*>(*it);
          p->e->addVert(EVENT0, p);
          break;
        }
      case IFF:
        {
          // Add IntEdge to SumFace.
          IntEdge *e = dynamic_cast<IntEdge*>(*it);
          assert(e->tail == 0);
          e->getTwin()->getFace()->liveEdges.insert(e);
          break;
        }
      default:
        assert(0);
        break;
      }
    }
  
  // Create SubEdges of SumEdges and enqueue their events.  Set tails
  // of IntEdges.
  for (int ei = ERROR; ei <= SUBF; ei++)
    for (unordered_set<Obj*, Equal, Equal>::iterator it = sweepSets[ei].begin();
         it != sweepSets[ei].end(); ++it) {
      switch ((*it)->getType()) {
      case SVE:
      case SEV:
        {
          SumEdge *e = dynamic_cast<SumEdge*>(*it);
          if (e->smallestTwin()) {
            if (e == (C4D::SumEV *) 0x69c860)
              cout << "SEV this is it";
            vector<SubEdge*> subs = e->makeSubEdges(sweepSets);
            for (int i = 0; i < subs.size(); i++)
              sweepSets[subs[i]->getType()].insert(subs[i]);
            for (int i = 1; i < subs.size()-1; i++)
              subs[i]->addEvents();
          }
        }
        break;
      case SVF:
      case SFV:
      case SEE:
        break;
      case IEF:
        {
          IntEF *p = dynamic_cast<IntEF*>(*it);
          for (int i = 0; i < p->e->getFaces().size(); i++) {
            IntEdge *e = p->e->getFaces()[i]->findLiveEdge(p->f);
            assert(e->tail == 0);
            e->setTail(p);
          }
          for (int i = 0; i < p->e->getTwin()->getFaces().size(); i++) {
            IntEdge *e = p->f->findLiveEdge(p->e->getTwin()->getFaces()[i]);
            assert(e->tail == 0);
            e->setTail(p);
          }
          break;
        }
      case IFF:
        {
          IntEdge *e = dynamic_cast<IntEdge*>(*it);
          if (e == (Obj *) 0x2990ce0)
            cout << "this is it IntEdge" << endl;
          assert(e->tail != 0 || e->twin->tail != 0);
          if (e->tail == 0 || e->twin->tail == 0) {
            const vector<Edge*> &edges = e->face->getOuterLoop();
            for (int i = 0; i < edges.size(); i++)
              if (e->twin->face->contains(edges[i]->tail)) {
                if (e->tail == 0)
                  e->setTail(edges[i]->tail);
                else
                  e->getTwin()->setTail(edges[i]->tail);
                break;
              }
            assert(e->tail != 0 && e->twin->tail != 0);
          }
        }
        break;
      case IFFF:
        break;
      case SUBE:
        break;
      default:
        assert(0);
        break;
      }
    }

  double sweepStart = getTime();
  rootTime = 0;
  rootHighTime = 0;
  predTime = 0;
  highTime = 0;

  // Intersect IntEdges and create IntFFFs.  Create SubEdges of
  // IntEdges.  Start IntFFFs at EVENT0.  Start SubEdges at EVENT0.
  for (int ei = ERROR; ei <= SUBF; ei++)
    for (unordered_set<Obj*, Equal, Equal>::iterator it = sweepSets[ei].begin();
         it != sweepSets[ei].end(); ++it) {
      switch ((*it)->getType()) {
      case SVE:
      case SEV:
        break;
      case SVF:
      case SFV:
      case SEE:
        {
          SumFace *f = dynamic_cast<SumFace*>(*it);
          f->intersectEdges(EVENT0);
        }
        break;
      case IEF:
        break;
      case IFF:
        {
          IntEdge *e = dynamic_cast<IntEdge*>(*it);
          if (e->smallestTwin()) {
            vector<SubEdge*> subs = e->makeSubEdges(sweepSets);
            for (int i = 0; i < subs.size(); i++)
              sweepSets[subs[i]->getType()].insert(subs[i]);
            for (int i = 1; i < subs.size(); i++)
              sweepSets[subs[i]->getTail()->getType()].insert(subs[i]->getTail());
            for (int i = 1; i < subs.size()-1; i++)
              subs[i]->addEvents();
          }
        }
        break;
      case IFFF:
        // Live IntFFFs are given EVENT0 as temporary starting event.
        (*it)->life.events.push_back(EVENT0);
        break;
      case SUBE:
        // Live SubEdges are given EVENT0 as temporary starting event.
        (*it)->life.events.push_back(EVENT0);
        break;
      default:
        assert(0);
        break;
      }
    }

  set<C4D::Obj *, C4D::Compare> objs1 = c4d.getObjects(c4d.EVENT0);
  enable();
  // compare(c4d.EVENT0, sweepSet, objs1);

  if (true) {
    set<Event*, Compare>::iterator it = eventList.begin();
    MeanAngle meanAngle(c4d.EVENT0->getAngle(), (*it)->getAngle(), 0.5);
    Event meanEvent(this, &meanAngle);
    // compare(&meanEvent, sweepSet, objs1);
  }

  // Sweep.
  int nEvents = 0;
  for (set<Event*, Compare>::iterator it = eventList.begin();
       it != eventList.end(); ++it) {
    Event *event = *it;
    // vjm debug
    global_event = event;
    nEvents++;
    if (nEvents % 10 == 0)
#ifdef NOTIME
      cout << nEvents << endl;
#else
    cout << nEvents << " " << objCount << " "
         << rootTime << " " << rootHighTime << " "
         << predTime << " " << highTime << " "
         << (getTime() - sweepStart) << " " << idenTime << endl;
#endif

    if (nEvents == 1634)
      cout << "this is it Event" << endl;

    assert(0);
#ifdef SWEEP_ALL_FACES
    // Raw SubEdge Event:
    if (event->polyType() == AnglePoly::NNN)
      dynamic_cast<PolyNNN*>(event->getPoly())->handleEvent(event);
    else while (event->vanish.size() > 0 && event->appear.size() > 0 &&
                event->vanish[0] == event->appear[0]) {
        switch (event->vanish[0]->getType()) {
        case SVE:
        case SEV:
          dynamic_cast<SumEdge*>(event->vanish[0])->handleEFF(event);
          break;
        case SVF:
        case SFV:
        case SEE:
          handleFFFF(event);
          break;
        default:
          assert(0);
          break;
        }
      }
#endif

    // Do swaps.
    for (pair<Obj*,Obj*> oo : event->swap) {
      if (oo.first->getType() == IEF) {
        IntEF *v1 = dynamic_cast<IntEF*>(oo.first);
        IntEF *v2 = dynamic_cast<IntEF*>(oo.second);
        assert(v1->e == v2->e);
        SumEdge *e = v1->e;
        if (!e->smallestTwin())
          e = e->getTwin();
        int i1 = e->findVert(v1);
        int i2 = e->findVert(v2);
        assert(i1 != -1 && i2 != -1);
        if (i1+1 == i2)
          e->handleSwap(event, i1);
        else if (i2+1 == i1)
          e->handleSwap(event, i2);
        else
          assert(0);
      }
      else {
        assert(oo.first->getType() == IFFF);
        IntFFF *v1 = dynamic_cast<IntFFF*>(oo.first);
        IntFFF *v2 = dynamic_cast<IntFFF*>(oo.second);
        IntEdge *e = v1->commonEdge(v2);
        assert(e != 0);
        if (!e->smallestTwin())
          e = e->getTwin();
        int i1 = e->findVert(v1);
        int i2 = e->findVert(v2);
        assert(i1 != -1 && i2 != -1);
        if (i1+1 == i2)
          e->handleSwap(event, i1);
        else if (i2+1 == i1)
          e->handleSwap(event, i2);
        else
          assert(0);
      }
    }    

    // Update their lives below?

    // In case the event is the collision of two SumEdges.
    // SumEdge *ee[] = { 0, 0 };
    vector<SumEdge*> ee[2];

    for (Obj *obj : event->vanish) {
      switch (obj->getType()) {
      case SVE:
      case SEV:
        {
          // If a SumEdge vanishes, so do all its SubEdges.
          SumEdge *edge = dynamic_cast<SumEdge*>(obj);
          if (edge->smallestTwin()) {
            edge->vanishSubs(event);
            edge->verts.clear();
          }
        }
        break;
      case SVF:
      case SFV:
      case SEE:
        {
	  SumFace *face = dynamic_cast<SumFace*>(obj);
          const vector<Edge*> &edges = face->getOuterLoop();
	  for (int j = 0; j < edges.size(); j++) {
	    SumEdge *edge = dynamic_cast<SumEdge*>(edges[j]);
	    edge->removeFace(face);
          }
        }
        break;
      case IEF:
        {
          bool do_break = false;
          IntEF *vert = dynamic_cast<IntEF*>(obj);
          // Check if the IntEF is vanishing because its SumEdge is.
	  if (event->vanish.find(vert->e) != event->vanish.end())
	    break;
          // Vanishing because its SumFace is.
	  if (event->vanish.find(vert->f) != event->vanish.end()) {
	    vert->e->handleRemove(event, vert);
	    break;
	  }
#ifdef BLEEN
          for (int i = 0; i < event->vanish.size(); i++)
            if (event->vanish[i] == vert->e)
              do_break = true;
          if (do_break)
            break;
          // Vanishing because its SumFace is.
          for (int i = 0; i < event->vanish.size(); i++)
            if (event->vanish[i] == vert->f) {
              vert->e->handleRemove(event, vert);
              do_break = true;
            }
          if (do_break)
            break;
#endif
          // Result of VF event.
          RootAngle *root = dynamic_cast<RootAngle*>(event->getAngle());
          SumEdge *edge = vert->e;
          for (int i = 0; i < 2; i++) {
            if (PolyVF::sameEvent(event, edge->getTail(), vert->f)) {
              edge->handleTail(event, vert);
              do_break = true;
              break;
            }
            edge = edge->getTwin();
          }
          if (do_break)
            break;

          // Result of EE event.

          // Already seen:
          for (int i = 0; !do_break && i < ee[0].size(); i++)
            if (((vert->e == ee[0][i] ||
                  vert->e == ee[0][i]->getTwin()) &&
                 ee[1][i]->hasFace(vert->f)) ||
                (((vert->e == ee[1][i] ||
                   vert->e == ee[1][i]->getTwin()) &&
                  ee[0][i]->hasFace(vert->f))))
              do_break = true;
          if (do_break)
            break;

          const vector<Edge*> &edges = vert->f->getOuterLoop();
          for (int i = 0; i < edges.size(); i++) {
            SumEdge *edge2 = dynamic_cast<SumEdge*>(edges[i]);
            if (PolyEE::sameEvent(event, edge, edge2)) {
              // Should not have seen an ee collision already.
              ee[0].push_back(edge);
              ee[1].push_back(edge2);

              // ee[0].back()->handleEE(event, ee[1].back());
              // ee[1].back()->handleEE(event, ee[0].back());

              do_break = true;
              break;
            }
          }
          if (do_break)
            break;

          assert(0);
        }
        break;
      case IFF:
        {
          // If a IntEdge vanishes, so do all its SubEdges and IntFFFs.
          // vjm If its endpoint is a sum vertex, remove the edge from its list.
          IntEdge *edge = dynamic_cast<IntEdge*>(obj);
	  if (edge == (Obj *) 0x2991070)
	    cout << "this is it IntEdge" << endl;
          SumFace *face = edge->getTwin()->getFace();
          set<IntEdge*, Compare>::iterator it = face->liveEdges.find(edge);
          assert(it != face->liveEdges.end());
          face->liveEdges.erase(it);
          if (!edge->smallestTwin())
            break;
          edge->vanishSubs(event);
          edge->setTail(0);
          edge->getTwin()->setTail(0);
          for (int i = 0; i < edge->verts.size(); i++)
            event->makeVanish(edge->getVert(i));
	  edge->verts.clear();
        }
        break;
      case IFFF:
        // event has already been added to life.
	{
	  IntFFF *v = (IntFFF*) obj;
	  for (int i = 0; i < 3; i++) {
	    int j = (i+1)%3;
	    IntEdge *e = v->faces[i]->findLiveEdge(v->faces[j]);
	    if (e == 0)
	      continue;
	    if (!e->smallestTwin())
	      e = e->getTwin();
	    if (e->containsVert(v))
	      e->handleRemove(event, v);
	  }
	}
        break;
      case SUBE:
        // event has already been added to life.
        break;
      default:
        assert(0);
        break;
      }
    }
    
    // Add new intersection edges to faces so they are available to look up.
    for (Obj *obj : event->appear)
      if (obj->getType() == IFF) {
	IntEdge *e = dynamic_cast<IntEdge*>(obj);
	assert(e->tail == 0);
        e->getTwin()->getFace()->liveEdges.insert(e);
      }
    
    // List of vert/face incidences.
    vector<SumVV*> vfVerts;
    vector<SumFace*> vfFaces;
    vector<SumEdge*> vfEdges;

    // List of edge/face incidences.
    vector<IntEdge*> efIntEdges;
    vector<SumEdge*> efSumEdges;
    vector<SumFace*> efFaces;

    for (Obj *obj : event->appear) {
      switch (obj->getType()) {
      case SVE:
      case SEV:
	{
          SumEdge *edge = dynamic_cast<SumEdge*>(obj);
          assert(edge->verts.size() == 0);
	}
        break;
      case SVF:
      case SFV:
      case SEE:
	{
	  SumFace *face = dynamic_cast<SumFace*>(obj);
	  const vector<Edge*> &edges = face->getOuterLoop();
	  for (int j = 0; j < edges.size(); j++) {
	    SumEdge *edge = dynamic_cast<SumEdge*>(edges[j]);
	    edge->addFace(face);
	    if (edge->verts.size() != 0 || edge->twin->verts.size() != 0) {
	      vector<Obj*> &verts = edge->verts.size() != 0 ?
		edge->verts : edge->twin->verts;
	      for (int k = 0; k < verts.size(); k++) {
		IntEF *ef = dynamic_cast<IntEF*>(verts[k]);
		if (event->findVanish(ef) != -1 || event->findAppear(ef) != -1)
		  continue;
		IntEdge *e = face->findLiveEdge(ef->f);
		if (ef->e == edge)
		  e->setTail(ef);
		else if (ef->e == edge->twin)
		  e->getTwin()->setTail(ef);
		else
		  assert(0);
	      }
	    }
	  }
	}
	break;
      case IEF:
        {
          IntEF *vert = dynamic_cast<IntEF*>(obj);
          SumEdge *edge = vert->e;
          SumFace *face = vert->f;

          // Update tail of IFF
          for (int i = 0; i < vert->e->getFaces().size(); i++) {
            IntEdge *e = vert->e->getFaces()[i]->findLiveEdge(face);
            if (event->findAppear(e) == -1) {
              e->vanishSub(event, -1);
              e->setTail(vert);
              e->appearSub(event, -1);
            }
            else
              e->setTail(vert);
          }
          for (int i = 0; i < vert->e->getTwin()->getFaces().size(); i++) {
            IntEdge *e = face->findLiveEdge(vert->e->getTwin()->getFaces()[i]);
            if (event->findAppear(e) == -1) {
              e->vanishSub(event, -1);
              e->setTail(vert);
              e->appearSub(event, -1);
            }
            else
              e->setTail(vert);
          }

          bool do_break = false;

          // Result of EE event?
          // Already seen:
          for (int i = 0; !do_break && i < ee[0].size(); i++)
            if (((vert->e == ee[0][i] ||
                  vert->e == ee[0][i]->getTwin()) &&
                 ee[1][i]->hasFace(face)) ||
                (((vert->e == ee[1][i] ||
                   vert->e == ee[1][i]->getTwin()) &&
                  ee[0][i]->hasFace(face))))
              do_break = true;
          if (do_break)
            break;

          // Not seen yet:
          const vector<Edge*> &edges = face->getOuterLoop();
          for (int i = 0; i < edges.size(); i++) {
            SumEdge *edge2 = dynamic_cast<SumEdge*>(edges[i]);
            if (PolyEE::sameEvent(event, edge, edge2)) {
              ee[0].push_back(edge);
              ee[1].push_back(edge2);

              // ee[0].back()->handleEE(event, ee[1].back());
              // ee[1].back()->handleEE(event, ee[0].back());

              do_break = true;
              break;
            }
          }
          if (do_break)
            break;

          // IEF appears because edge appears.
          // Edge is new, so add all the vertices before constructing subedges.
	  if (event->appear.find(vert->e) != event->appear.end()) {
	    vert->e->addVert(event, vert);
            break;
	  }

	  bool face_appear = false;
          // IEF appears because face appears.
          // Insert vertex into existing edge.
	  if (event->appear.find(face) != event->appear.end()) {
	    vert->e->handleInsert(event, vert);
	    face_appear = true;
	  }

          // Result of VF event.
          // Could be PolyVF or PolyEF.
	  RootAngle *root = dynamic_cast<RootAngle*>(event->getAngle());
	  PolyEF *poly = root->getPoly()->getType() == AnglePoly::EF2 ?
	    dynamic_cast<PolyEF*>(root->getPoly()) : 0;
	  bool polyContainsF = (poly != 0 && poly->contains(face));

	  // Face can appear and sweep through vertex.
	  if (face_appear && !polyContainsF)
	    break;

          for (int i = 0; i < 2; i++) {
            if (PolyVF::sameEvent(event, edge->getTail(), face)) {
	      if (!face_appear)
		edge->handleTail(event, vert);

              SumVV *tail = edge->getTail();
              if (tail->intEdges.size() > 0) {
                int j;
                for (j = 0; j < vfVerts.size(); j++)
                  if (tail == vfVerts[j] && face == vfFaces[j])
                    break;
                if (j == vfVerts.size()) {
                  vfVerts.push_back(tail);
                  vfFaces.push_back(face);
                  vfEdges.push_back(edge);
                }
              }
              do_break = true;
              break;
            }
            edge = edge->getTwin();
          }
          if (do_break || face_appear)
            break;

          assert(0);
        }
        break;
      case IFF:
        {
          IntEdge *e = dynamic_cast<IntEdge*>(obj);
          assert(e->tail != 0 || e->twin->tail != 0);
          if (e->tail == 0 || e->twin->tail == 0) {
            const vector<Edge*> &edges = e->face->getOuterLoop();
            for (int i = 0; i < edges.size(); i++)
              if (e->twin->face->contains(edges[i]->tail)) {
                if (e->tail == 0)
                  e->setTail(edges[i]->tail);
                else
                  e->getTwin()->setTail(edges[i]->tail);
                break;
              }

            // DEBUGGING CODE:
	    if (e->tail == 0 || e->twin->tail == 0) {
	      int n = 0;
	      const vector<Edge*> &edges = e->face->getOuterLoop();
	      for (int i = 0; i < edges.size(); i++) {
		cout << edges[i] << " " << edges[i]->tail << endl;
		if (e->twin->face->contains(edges[i]->tail))
		  n++;
	      }
	      cout << endl;
	      const vector<Edge*> &edges2 = e->twin->face->getOuterLoop();
	      for (int i = 0; i < edges2.size(); i++)
		cout << edges2[i] << " " << edges2[i]->tail << endl;
	      cout << "n " << n << endl;

	      // for (unordered_set<IntEF *, Equal, Equal>::iterator it = intEFs.begin();
              for (auto it = intEFs.begin();
		   it != intEFs.end(); ++it) {
		IntEF *p = *it;
		if ((p->e->hasFace(e->getFace()) ||
                     p->e->hasFace(e->getTwin()->getFace())) &&
		    (p->f == e->face ||
		     p->f == e->twin->face))
		  cout << "this is it" << endl;
	      }

	      IntEdge *edge = e; 
	      SumFace *face01[] = { e->getFace(), e->getTwin()->getFace() };
	      const vector<Edge*> *edges01[] =
		{ &face01[0]->getOuterLoop(), &face01[1]->getOuterLoop() };
	      for (int h = 0; h < 2; h++) {
		for (int j = 0; j < edges01[h]->size(); j++) {
		  SumEdge *e = dynamic_cast<SumEdge*>((*edges01[h])[j]);
		  for (int k = 0; k < 2; k++) {
		    int tside = face01[!h]->side(event, e->tail);
		    int hside = face01[!h]->side(event, e->twin->tail);
		    assert(tside != 0 && hside != 0);
		    if (tside < 0 && hside > 0 &&
			face01[!h]->intersectsInside(event, e)) {
		      IntEF *ef = IntEF::get(e, face01[!h]);
		      cout << "this is it" << endl;
		    }
		    e = e->getTwin();
		  }
		}
	      }
	    }
            assert(e->tail != 0 && e->twin->tail != 0);
          }
        }
        break;
      case IFFF:
        break;
      case SUBE:
        break;
      default:
        assert(0);
        break;
      }
    }

    for (Obj *obj : event->appear) {
      switch (obj->getType()) {
      case SVE:
      case SEV:
        {
          SumEdge *edge = dynamic_cast<SumEdge*>(obj);
          if (!edge->smallestTwin())
            break;
          edge->appearSubs(event);
        }
        break;
      case SVF:
      case SFV:
      case SEE:
        dynamic_cast<SumFace*>(obj)->intersectEdges(event, true);
        break;
      case IEF:
        break;
      case IFF:
        {
          // If a IntEdge appears, so do all its SubEdges and IntFFFs.
          IntEdge *edge = dynamic_cast<IntEdge*>(obj);
	  if (edge == (IntEdge *) 0x4c5eb80)
	    cout << "IFF this is it" << endl;

	  if (edge->smallestTwin())
	    edge->appearSubs(event);

#ifdef BLEEN
          // If it appears because one of its sum faces just appeared,
          // then all the intersections on that sum face were
          // calculated already.  Create its subs.
          if (event->findAppear(edge->face) != -1 ||
              event->findAppear(edge->twin->face) != -1)
            break;
#endif

          RootAngle *root = dynamic_cast<RootAngle*>(event->getAngle());
          // VF and EE events also handled already.
          if (root->getPoly()->getType() == AnglePoly::VF ||
              root->getPoly()->getType() == AnglePoly::EE)
            break;

          PolyEF *poly = dynamic_cast<PolyEF*>(root->getPoly());
          SumFace *sface = edge->getFace();
          SumFace *tface = edge->getTwin()->getFace();
          bool polyContainsFace = poly->contains(sface);
          bool polyContainsTwin = poly->contains(tface);

          // Prefer to deal with edge which has non-special face to its left.
          if (polyContainsFace)
            break;

	  // Not the result of an edge/face incidence.
          if (!polyContainsTwin)
	    break;

          const vector<Edge*> &edges = sface->getOuterLoop();
          SumEdge *sumE = 0;
          for (int i = 0; i < edges.size(); i++)
            if (poly->contains(edges[i]))
              sumE = dynamic_cast<SumEdge*>(edges[i]);

	  // Only a vertex/face coincidence.
          if (sumE == 0)
	    break;
          
          SumEdge *smallE = sumE->smallestTwin() ? sumE : sumE->getTwin();
          bool seenEF = false;
          for (int i = 0; i < efSumEdges.size(); i++)
            if (efSumEdges[i] == smallE && efFaces[i] == tface)
              seenEF = true;
          if (seenEF)
            break;

          efSumEdges.push_back(smallE);
	  efIntEdges.push_back(edge);
          efFaces.push_back(tface);

          // tface->handleEF(event, edge, smallE);
        }
        break;
      case IFFF:
        break;
      case SUBE:
        break;
      default:
        assert(0);
        break;
      }
    }

    for (int i = 0; i < ee[0].size(); i++) {
      ee[0][i]->handleEE(event, ee[1][i]);
      ee[1][i]->handleEE(event, ee[0][i]);
    }

    for (int i = 0; i < vfVerts.size(); i++)
      vfVerts[i]->handleVF(event, vfFaces[i], vfEdges[i]);

    for (int i = 0; i < efFaces.size(); i++)
      efFaces[i]->handleEF(event, efIntEdges[i], efSumEdges[i]);

#ifdef BLEEN    
    for (int iv = 0; iv < newFFFs.size(); iv++) {
      IntFFF *v = newFFFs[iv];
      for (int i = 0; i < 3; i++) {
        SumFace *fi = v->faces[i];
        SumFace *fj = v->faces[(i+1)%3];
        IntEdge *eij = fi->findLiveEdge(fj);
        eij->handleInsert(event, v);
      }
    }
#endif

    // Eliminate interim subs.
    for (auto it = event->vanish.begin(); it != event->vanish.end();) {
      auto jt = it++;
      vector<PTR<Event> > &objEvents = (*jt)->life.events;
      int s = objEvents.size();
      assert(s >= 2);
      if (objEvents[s-2] == objEvents[s-1])
	event->vanish.erase(jt);
    }
#ifdef BLEEN      
    int n = 0;
    for (int i = 0; i < event->vanish.size(); i++) {
      Obj *obj = event->vanish[i];
      vector<PTR<Event> > &objEvents = obj->life.events;
      int s = objEvents.size();
      assert(s >= 2);
      if (objEvents[s-2] != objEvents[s-1])
	event->vanish[n++] = obj;
    }
    event->vanish.resize(n);	
#endif

    for (auto it = event->appear.begin(); it != event->appear.end();) {
      auto jt = it++;
      vector<PTR<Event> > &objEvents = (*jt)->life.events;
      int s = objEvents.size();
      if (s >= 2 && objEvents[s-2] == objEvents[s-1])
	event->appear.erase(jt);
    }
#ifdef BLEEN
    n = 0;
    for (int i = 0; i < event->appear.size(); i++) {
      Obj *obj = event->appear[i];
      vector<PTR<Event> > &objEvents = obj->life.events;
      int s = objEvents.size();
      if (s < 2 || objEvents[s-2] != objEvents[s-1])
	event->appear[n++] = obj;
    }
    event->appear.resize(n);	
#endif

    for (Obj *obj : event->vanish) {
      // remove from sweepSet
      // DEBUG!!!
      // assert(sweepSets[obj->getType()].find(obj) != sweepSets[obj->getType()].end());
      sweepSets[obj->getType()].erase(sweepSets[obj->getType()].find(obj));
    }

    for (Obj *obj : event->appear) {
      if (obj == (C4D::SubEdge *) 0x29c81a0)
	cout << "this is it SubEdge" << endl;

      // add to sweepSet
      // DEBUG!!!
      // assert(sweepSets[obj->getType()].find(obj) == sweepSets[obj->getType()].end());
      sweepSets[obj->getType()].insert(obj);
    }

    set<Event*, Compare>::iterator jt = it;
    jt++;
    if (jt != eventList.end()) {
      // assert((*jt)->id == event->id + 1);
      MeanAngle meanAngle(event->getAngle(), (*jt)->getAngle(), 0.5);
      Event meanEvent(this, &meanAngle);
      // if (nEvents >= 1000000000)
      // objs1 = c4d.getObjects(&meanEvent);
      enable();
      currentEvent = event;
      // if (nEvents >= 1000000000)
      // compare(&meanEvent, sweepSet, objs1);
    }
    else {
      objs1 = c4d.getObjects(EVENT0);
      enable();
      // compare(EVENT0, sweepSet, objs1);
    }
  }

  // Replace EVENT0 for live IFFFs and SubEdges with actual starting event.
  for (int ei = ERROR; ei <= SUBF; ei++)
    for (unordered_set<Obj*, Equal, Equal>::iterator it = sweepSets[ei].begin();
         it != sweepSets[ei].end(); ++it)
      if ((*it)->getType() == IFFF || (*it)->getType() == SUBE) {
        assert((*it)->life.events[0] == EVENT0);
        (*it)->life.events[0] = (*it)->life.events.back();
        (*it)->life.events.pop_back();
      }
}

int sweep_count = 0;
void C4D::checkEdgeList (Event *event, SumFace *base, IntEdge *e,
                         unordered_set<IntEdge*, Equal, Equal> &eSet) {
  vector<IntFFF*> verts;
  for (auto it = eSet.begin(); it != eSet.end(); it++) {
    IntEdge *d = *it;
    int det = d->getTwin()->getFace()->side(event, e->tail);
    int deh = d->getTwin()->getFace()->side(event, e->twin->tail);
    int edt = e->getTwin()->getFace()->side(event, d->tail);
    int edh = e->getTwin()->getFace()->side(event, d->twin->tail);

    if (det * deh != -1)
      continue;

    if (edt * edh != -1)
      continue;

    assert(edt * det == -1);
    IntFFF *v = edt > 0 ?
      IntFFF::make(base, d->getTwin()->getFace(), e->getTwin()->getFace()) :
      IntFFF::make(base, e->getTwin()->getFace(), d->getTwin()->getFace());
      
    int i = verts.size();
    verts.push_back(0);
    while (i > 0 && OrderUVE(event, verts[i-1], v, e) > 0) {
      verts[i] = verts[i-1];
      i--;
    }
    verts[i] = v;
  }

  if (!(verts.size() == e->verts.size()))
    cout << "sweep_count " << sweep_count << endl;
  assert(verts.size() == e->verts.size());
  for (int i = 0; i < verts.size(); i++) {
    if (!(e->verts[i] == verts[i]))
      cout << "sweep_count " << sweep_count << endl;
    assert(e->verts[i] == verts[i]);
  }
}

void C4D::intersectEdges (SumFace *base, Event *event, IntEdge *d, IntEdge *e, 
                          int &nFFFs, bool insert) {
  int di = -1, ei;
  for (ei = 0; ei < e->verts.size(); ei++)
    if ((di = d->findVert(dynamic_cast<Vert*>(e->verts[ei]))) != -1)
      break;
          
  if (!insert && di == -1)
    return;

  bool intersecting = false;
  int det = d->getTwin()->getFace()->side(event, e->tail);
  int deh = d->getTwin()->getFace()->side(event, e->twin->tail);
  int edt = 0, edh = 0;
  if (det * deh == -1) {
    edt = e->getTwin()->getFace()->side(event, d->tail);
    edh = e->getTwin()->getFace()->side(event, d->twin->tail);
    if (edt * edh == -1)
      intersecting = true;
  }

  if ((di != -1) == intersecting)
    return;

  if (!insert && di != -1) {
    e->verts[ei]->life.events.push_back(event);
    e->verts.erase(e->verts.begin() + ei);
    generateSwaps(e, ei-1);
    d->verts.erase(d->verts.begin() + di);
    generateSwaps(d, di-1);
    return;
  }

  assert(edt * det == -1);
  IntFFF *fff = edt > 0 ?
    IntFFF::make(base, d->getTwin()->getFace(), e->getTwin()->getFace()) :
    IntFFF::make(base, e->getTwin()->getFace(), d->getTwin()->getFace());
  fff->life.events.push_back(event);
  nFFFs++;
          
  di = d->insertIndex(event, fff);
  d->verts.insert(d->verts.begin() + di, fff);
  generateSwaps(d, di-1);
  generateSwaps(d, di);

  ei = e->insertIndex(event, fff);
  e->verts.insert(e->verts.begin() + ei, fff);
  generateSwaps(e, ei-1);
  generateSwaps(e, ei);
}

void C4D::sweep2 (int iFace) {
  // C4D c4d(this);
  vector<Face *> sumFaces = getSumFaces();
  cout << "sumFaces " << sumFaces.size() << endl;
  set<Obj*, Compare> objs = getSums(sumFaces);
  sweep2(iFace, sumFaces, objs);
}

void C4D::sweep2 (vector<int> &iFaces, int from, int to) {
  // C4D c4d(this);
  vector<Face *> sumFaces = getSumFaces();
  cout << "sumFaces " << sumFaces.size() << endl;
  set<Obj*, Compare> objs = getSums(sumFaces);
  for (int i = from; i < to; i++)
    sweep2(iFaces[i], sumFaces, objs);
}

void C4D::sweep2 (int iFace, vector<Face *> & sumFaces, set<Obj*, Compare> &objs) {
  SumFace *base = dynamic_cast<SumFace*>(sumFaces[iFace]);
  getInts(sumFaces, base, objs);

  // Put Events of each IntEF and IntEdge into eventList.
  // Inform events of who appears and vanishes.
  // Eliminate duplicate Events.
  // Objects live at angle 0 "appear" at EVENT0.
  int nBaseV = 0, nEdgeV = 0, nEdges = 0;

  PTR<Event> eVENT0 = new Event(this, angle0);

  for (set<Obj*, Compare>::iterator it = objs.begin();
       it != objs.end(); ++it)
    if ((*it)->getType() == IEF || (*it)->getType() == IFF) {
      Obj* obj = *it;
      if ((*it)->getType() == IEF) {
	IntEF *e = dynamic_cast<IntEF*>(*it);
	if (e->f == base)
	  nBaseV++;
	else {
	  nEdgeV++;
          // WARNING WILL ROBINSON!!!!
          (*it)->life = base->life.intersect((*it)->life);
        }
      }
      else
	nEdges++;
      (*it)->informEvents();
      if ((*it)->contains(EVENT0))
        eVENT0->appear.insert(*it);
    }

  cout << "nBaseV " << nBaseV << " nEdgeV " << nEdgeV << " nEdges " << nEdges << endl;

  // Add eVENT0 as an imaginary event.
  eventList.insert(eVENT0);

  if (noSweep)
    return;

  unordered_set<IntEF*, Equal, Equal> vSet;
  unordered_set<IntEdge*, Equal, Equal> eSet;
  map<Vert*, set<IntEdge*, Compare>, Compare> vv2ff;

  const vector<Edge*> &edges = base->getOuterLoop();
  for (int i = 0; i < edges.size(); i++)
    edges[i]->face = base;

  int nSEswaps = 0, nFFswaps = 0, nFFFs = 0;

  Event *prevent = 0;
  for (auto it = eventList.begin(); it != eventList.end(); ++it) {
    sweep_count++;
    Event *event = *it;
    currentEvent = event;
    if (false && sweep_count > 0 && prevent != 0) {
      if (sweep_count % 100 == 0)
	cout << "sweep_count " << sweep_count << endl;
      MeanAngle meanAngle(prevent->getAngle(), event->getAngle(), /*0.99999*/ /*0.00001*/ 0.5);
      Event meanEvent(this, &meanAngle);
      for (auto eit = eSet.begin(); eit != eSet.end(); eit++) {
	IntEdge *e = *eit;
        checkEdgeList(&meanEvent, base, e, eSet);
      }
    }
    prevent = event;

    for (Obj *obj : event->vanish) {
      Type type = obj->getType();
      if (type == IEF) {
        IntEF *v = dynamic_cast<IntEF*>(obj);
        if (v->f == base) {
          auto it = vSet.find(v);
          assert(it != vSet.end());
          vSet.erase(it);
          for (int j = 0; j < v->e->faces.size(); j++) {
            IntEdge *d = v->e->faces[j]->intEdge;
            assert(d->tail == v);
            d->tail = 0;
          }
          for (int j = 0; j < v->e->getTwin()->faces.size(); j++) {
            IntEdge *d = v->e->getTwin()->faces[j]->intEdge;
            assert(d->twin->tail == v);
            d->twin->tail = 0;
          }
          assert(v->e->verts.size() == 1);
          v->e->verts.clear();
        }
        else {
          SumEdge *e = v->e;
          if (e->face != base) {
            e = e->getTwin();
            assert(e->face == base);
          }

          if (e == v->e) {
            assert(v->f->intEdge->twin->tail == v);
            v->f->intEdge->twin->tail = 0;
          }
          else {
           assert(v->f->intEdge->tail == v);
            v->f->intEdge->tail = 0;
          }

          int i = e->findVert(v);
          assert(i != -1);
          e->verts.erase(e->verts.begin() + i);
          generateSwaps(e, i-1);
        }
      }
      else if (type == IFF) {
        IntEdge *e = dynamic_cast<IntEdge*>(obj);
	if (e->face != base)
	  continue;

        auto eit = eSet.find(e);
        assert(eit != eSet.end());
        eSet.erase(eit);

        SumFace *f = e->getTwin()->getFace();
        for (int j = 0; j < f->numEdges(); j++) {
	  SumEdge *d = f->getEdge(j);
	  if (d->verts.size() == 1) {
	    IntEF *v = dynamic_cast<IntEF*>(d->verts[0]);
	    assert(e->tail == v);
	    e->tail = 0;
	  }
	  else if (d->twin->verts.size() == 1) {
	    IntEF *v = dynamic_cast<IntEF*>(d->twin->verts[0]);
	    assert(e->twin->tail == v);
	    e->twin->tail = 0;
	  }
          d->removeFace(f);
	}
        assert(f->intEdge == e);
        f->intEdge = 0;

        for (int h = 0; h < e->verts.size(); h++) {
          IntFFF *v = dynamic_cast<IntFFF*>(e->verts[h]);
          v->life.events.push_back(event);
          SumFace *o = v->otherFace(base, e->getTwin()->getFace());
          IntEdge *d = o->intEdge;
          int j = d->findVert(v);
          assert(j != -1);
          d->verts.erase(d->verts.begin() + j);
          generateSwaps(d, j-1);
        }

        e->verts.clear();

        assert(e->tail == 0 || e->twin->tail == 0);
        if (e->tail != 0 || e->twin->tail != 0) {
          Vert *&v = e->tail != 0 ? e->tail : e->twin->tail;
          assert(v->getType() == SVV);
          auto it = vv2ff[v].find(e);
          assert(it != vv2ff[v].end());
          vv2ff[v].erase(it);
	  v = 0;
        }
      }
    }

    for (Obj *obj : event->appear) {
      Type type = obj->getType();
      if (type == IFF) {
        IntEdge *e = dynamic_cast<IntEdge*>(obj);
	if (e->face != base)
	  continue;

        SumFace *f = e->getTwin()->getFace();
        f->intEdge = e;
      }
    }

    set<IntEF*, Compare> bSet;

    for (Obj *obj : event->appear) {
      Type type = obj->getType();
      if (type != IEF)
        continue;
      IntEF *v = dynamic_cast<IntEF*>(obj);

      if (v->f == base) {
        vSet.insert(v);
        assert(v->e->verts.size() == 0);
        v->e->verts.push_back(v);

        for (int fi = 0; fi < v->e->faces.size(); fi++) {
          IntEdge *d = v->e->faces[fi]->intEdge;
          assert(d->tail == 0);
          d->tail = v;
        }
        for (int fi = 0; fi < v->e->getTwin()->faces.size(); fi++) {
          IntEdge *d = v->e->getTwin()->faces[fi]->intEdge;
          assert(d->twin->tail == 0);
          d->twin->tail = v;
        }

        for (auto eit = eSet.begin(); eit != eSet.end(); eit++) {
          IntEdge *e = *eit;
          if (e->tail != v && e->twin->tail != v)
            generateSwaps(v, e->getTwin()->getFace());
        }
      }
      else {
        bSet.insert(v);
        SumEdge *e = v->e;
        if (e->face != base)
          e = e->getTwin();
        assert(e->face == base);
        
        int vi = e->insertIndex(event, v);
        e->verts.insert(e->verts.begin() + vi, v);
        generateSwaps(e, vi-1);
        generateSwaps(e, vi);
        
        IntEdge *d = v->f->intEdge;
        if (v->e == e) {
          assert(d->twin->tail == 0);
          d->twin->tail = v;
        }
        else {
          assert(d->tail == 0);
          d->tail = v;
        }
      }
    }

    for (Obj *obj : event->appear) {
      Type type = obj->getType();
      if (type == IFF) {
        IntEdge *e = dynamic_cast<IntEdge*>(obj);
	if (e->face != base)
	  continue;

        SumFace *f = e->getTwin()->getFace();
        for (int j = 0; j < f->numEdges(); j++) {
          SumEdge *d = f->getEdge(j);
          d->addFace(f);
          if (d->verts.size() == 1) {
            IntEF *v = dynamic_cast<IntEF*>(d->verts[0]);
	    if (e->tail != 0) {
	      IntEF *t = dynamic_cast<IntEF*>(e->tail);
	      assert(vSet.find(t) == vSet.end());
	    }
            assert(e->tail == 0);
            e->tail = v;
          }
          else if (d->twin->verts.size() == 1) {
            IntEF *v = dynamic_cast<IntEF*>(d->twin->verts[0]);
            assert(e->twin->tail == 0);
            e->twin->tail = v;
          }
        }

        if (e->tail == 0 || e->twin->tail == 0) {
          assert(e->tail != 0 || e->twin->tail != 0);
          SumVV *v = e->getFace()->commonVert(e->getTwin()->getFace());
	  if (v == 0) {
	    for (int k = 0; k < base->numEdges(); k++) {
	      SumEdge *ek = base->getEdge(k);
	      for (int l = 0; l < ek->verts.size(); l++) {
		IntEF *vl = dynamic_cast<IntEF*>(ek->verts[l]);
		if (vl->f == e->twin->face)
		  cout << "this is it" << endl;
	      }
	    }

	    SumFace *face = e->getTwin()->getFace();
	    for (int k = 0; k < face->numEdges(); k++) {
	      SumEdge *ek = face->getEdge(k);
	      for (auto vit = vSet.begin(); vit != vSet.end(); vit++) {
		IntEF *vl = *vit;
		if (vl->e == ek || vl->e == ek->getTwin())
		  cout << "this is it" << endl;
	      }
	    }
	  }
          assert(v != 0);
          assert(e->tail != v && e->twin->tail != v);
          if (e->tail == 0)
            e->tail = v;
          else
            e->twin->tail = v;
          assert(vv2ff[v].find(e) == vv2ff[v].end());
          vv2ff[v].insert(e);
        }
      }
    }

    for (Obj *obj : event->appear) {
      Type type = obj->getType();
      if (type == IFF) {
        IntEdge *e = dynamic_cast<IntEdge*>(obj);
	if (e->face != base)
	  continue;

        for (auto dit = eSet.begin(); dit != eSet.end(); dit++) {
          IntEdge *d = *dit;
          intersectEdges(base, event, d, e, nFFFs, true);
        }

        eSet.insert(e);

        for (auto vit = vSet.begin(); vit != vSet.end(); vit++)
          generateSwaps(*vit, e->getTwin()->getFace());
      }
    }

    for (Obj *obj : event->vanish) {
      if (obj->getType() != IEF)
        continue;
      IntEF *u = dynamic_cast<IntEF*>(obj);
      if (u->f == base)
        continue;
      if (u->f->intEdge == 0)
        continue;
      SumEdge *ue = u->e;
      if (ue->face != base)
        ue = ue->getTwin();
      assert(ue->face == base);

      IntEdge *e = u->f->intEdge;
      Vert *vs[2] = { e->tail, e->twin->tail };
      for (int th = 0; th <= 1; th++) {
        if (vs[th]->getType() != IEF)
          continue;
        IntEF *v = dynamic_cast<IntEF*>(vs[th]);
        if (v->f == base)
          continue;
        SumEdge *ve = v->e;
        if (ve->face != base)
          ve = ve->getTwin();
        assert(ve->face == base);
        if (ue == ve)
          continue;
        SumVV *vv = 0;
        if (ue->twin->tail == ve->tail)
          vv = dynamic_cast<SumVV*>(ve->tail);
        else if (ve->twin->tail == ue->tail)
          vv = dynamic_cast<SumVV*>(ue->tail);
        else
          continue;
        assert(vv != 0);

        set<IntEdge*, Compare> &vvSet = vv2ff[vv];
        for (auto dit = vvSet.begin(); dit != vvSet.end(); dit++) {
          IntEdge *d = *dit;
          intersectEdges(base, event, d, e, nFFFs, false);
	}
        for (auto dit = vvSet.begin(); dit != vvSet.end(); dit++) {
          IntEdge *d = *dit;
          intersectEdges(base, event, d, e, nFFFs, true);
        }
      }
    }

    for (pair<Obj*,Obj*> oo : event->swap) {
      Obj *obj0 = oo.first;
      Obj *obj1 = oo.second;
      Type type0 = obj0->getType();
      Type type1 = obj1->getType();
      if (type1 == IEF) {
        assert(type0 == IEF);
        IntEF *v0 = dynamic_cast<IntEF*>(obj0);
        IntEF *v1 = dynamic_cast<IntEF*>(obj1);

        SumEdge *e = v0->e;
        if (e->face != base)
          e = e->getTwin();
        assert(e->face == base);

        int i0 = e->findVert(v0);
        if (i0 == -1)
          continue;
        int i1 = e->findVert(v1);
        if (i1 == -1)
          continue;

        if (bSet.find(v0) != bSet.end() || bSet.find(v1) != bSet.end()) {
          assert((OrderUVE(event, v0, v1, e) < 0) == (i0 < i1));
          continue;
        }

	nSEswaps++;
        e->verts[i0] = v1;
        e->verts[i1] = v0;

        if (i0 < i1) {
          generateSwaps(e, i0-1);
          generateSwaps(e, i1);
          assert(OrderUVE(event, v1, v0, e) < 0);
        }
        else if (i1 < i0) {
          generateSwaps(e, i1-1);
          generateSwaps(e, i0);
          assert(OrderUVE(event, v0, v1, e) < 0);
        }
        else
          assert(0);

        IntEdge *e0 = v0->f->intEdge;
        IntEdge *e1 = v1->f->intEdge;
        
        int h0 = e0->tail == v0 ? 0 : e0->verts.size()-1;
        int h1 = e1->tail == v1 ? 0 : e1->verts.size()-1;
        if (0 <= h0 && h0 < e0->verts.size() &&
            0 <= h1 && h1 < e1->verts.size() &&
            e0->verts[h0] == e1->verts[h1]) {
          IntFFF *vFFF = dynamic_cast<IntFFF*>(e0->verts[h0]);
          vFFF->life.events.push_back(event);
          e0->verts.erase(e0->verts.begin() + h0);
          e1->verts.erase(e1->verts.begin() + h1);
        }
        else if (true) {
          intersectEdges(base, event, e0, e1, nFFFs, false);
          intersectEdges(base, event, e0, e1, nFFFs, true);
          continue;
        }
        else {
          int e0e1t = e0->getTwin()->getFace()->side(event, e1->tail);
          int e0e1h = e0->getTwin()->getFace()->side(event, e1->twin->tail);
          if (e0e1t * e0e1h != -1)
            continue;
          int e1e0t = e1->getTwin()->getFace()->side(event, e0->tail);
          int e1e0h = e1->getTwin()->getFace()->side(event, e0->twin->tail);
          if (e1e0t * e1e0h != -1)
            continue;
          
          assert(e0e1t * e1e0t == -1);
          IntFFF *fff = e1e0t > 0 ?
            IntFFF::make(base, e0->getTwin()->getFace(), e1->getTwin()->getFace()) :
            IntFFF::make(base, e1->getTwin()->getFace(), e0->getTwin()->getFace());
          fff->life.events.push_back(event);
	  nFFFs++;

          if (e0->tail == v0) {
            e0->verts.insert(e0->verts.begin(), fff);
            generateSwaps(e0, 0);
          }
          else {
            e0->verts.insert(e0->verts.end(), fff);
            generateSwaps(e0, e0->verts.size()-2);
          }

          if (e1->tail == v1) {
            e1->verts.insert(e1->verts.begin(), fff);
            generateSwaps(e1, 0);
          }
          else {
            e1->verts.insert(e1->verts.end(), fff);
            generateSwaps(e1, e1->verts.size()-2);
          }
        }
      }
      else if (SVF <= type1 && type1 <= SEE) {
        assert(type0 == IEF);
        IntEF *v = dynamic_cast<IntEF*>(obj0);
        SumFace *f = dynamic_cast<SumFace*>(obj1);
        if (vSet.find(v) == vSet.end())
          continue;
        if (f->intEdge == 0)
          continue;
        IntEdge *e = f->intEdge;
        assert(eSet.find(e) != eSet.end());

        for (int h0 = 0, h1 = 0; 
             h0 < v->e->faces.size() || h1 < v->e->getTwin()->faces.size();
             h0 < v->e->faces.size() ? h0++ : h1++) {
          SumFace *g = h0 < v->e->faces.size() ? 
            v->e->faces[h0] : v->e->getTwin()->faces[h1];
          IntEdge *d = g->intEdge;
          intersectEdges(base, event, d, e, nFFFs, false);
        }
        for (int h0 = 0, h1 = 0; 
             h0 < v->e->faces.size() || h1 < v->e->getTwin()->faces.size();
             h0 < v->e->faces.size() ? h0++ : h1++) {
          SumFace *g = h0 < v->e->faces.size() ? 
            v->e->faces[h0] : v->e->getTwin()->faces[h1];
          IntEdge *d = g->intEdge;
          intersectEdges(base, event, d, e, nFFFs, true);
        }
      }
      else if (type1 == IFFF) {
        assert(type0 == IFFF);
        IntFFF *v0 = dynamic_cast<IntFFF*>(obj0);
        IntFFF *v1 = dynamic_cast<IntFFF*>(obj1);

        if (v0->life.events.back() == event ||
            v1->life.events.back() == event)
          continue;

        IntEdge *e = v0->commonEdge(v1);
        if (e->face != base)
          e = e->getTwin();
        assert(e->face == base);

        int i0 = e->findVert(v0);
        if (i0 == -1)
          continue;
        int i1 = e->findVert(v1);
        if (i1 == -1)
          continue;
        // assert(i0 == i1+1 || i0 == i1-1);
        
	nFFswaps++;
        e->verts[i0] = v1;
        e->verts[i1] = v0;
        
        if (i0 < i1) {
          generateSwaps(e, i0-1);
          generateSwaps(e, i1);
        }
        else if (i1 < i0) {
          generateSwaps(e, i1-1);
          generateSwaps(e, i0);
        }
        else
          assert(0);
      }
      else
        assert(0);
    }
  }

  cout << "nSEswaps " << nSEswaps << " nFFswaps " << nFFswaps << endl;
  cout << "nFFFs " << nFFFs << endl;
  cout << "nGCD100 " << nGCD100 << endl;
  cout << "nGCD212 " << nGCD212 << endl;
  cout << "sweep2 completed." << endl;
}

// Select m random integers in the range 0..n-1 and put them into
// inds[0..m-1] in increasing order.
void randomInds (int m, int n, int *inds) {
  int i = 0;
  for (int j = 0; j < n; j++)
    if (random() % (n - j) < (m - i))
      inds[i++] = j;
#ifdef BLEEN
  cout << "randomInds";
  for (int i = 0; i < m; i++)
    cout << " " << inds[i];
  cout << endl;
#endif
}
#ifdef BLEEN
void randomInds (int m, int n, int *inds) {
  int i;
  int k = 0;
  while (k < m) {
    int ind = random() % n;
    for (i = 0; i < k && inds[i] < ind; i++);
    if (i < k && inds[i] == ind)
      continue;
    for (int j = k; j > i; j--)
      inds[j] = inds[j-1];
    inds[i] = ind;
    k++;
  }
}
#endif

void permute (int n, int *is) {
  while (n > 0) {
    int i = random() % n;
    int temp = is[i];
    is[i] = is[--n];
    is[n] = temp;
  }
}

void permutation (int n, int *is) {
  for (int i = 0; i < n; i++)
    is[i] = i;
  permute(n, is);
}

C4D::EdgeAB *C4D::makeEdge (VertXYZ *a, VertXYZ *b) {
  EdgeAB *ab = a->getEdge(b);
  if (ab != 0)
    return ab;
  return new EdgeAB(a, b);
}

C4D::FaceABC *C4D::makeFace (VertXYZ *a, VertXYZ *b, VertXYZ *c) {
  FaceABC *abc = a->getFace(b, c);
  if (abc != 0)
    return abc;
  return new FaceABC(a, b, c);
}

// intersect im[0<m] with in[0<n], output io[0<o] with o returned.
// !!! not efficient for large m and n
int intersect (int *im, int m, int *in, int n, int *io) {
  int o = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (im[i] == in[j])
        io[o++] = im[i];
  return o;
}

bool contains (int ind, int *inds, int n) {
  for (int i = 0; i < n; i++)
    if (ind == inds[i])
      return true;
  return false;
}

bool Face4::matches (int nF[4], int indsF[4][3]) {
  for (int i = 0; i < 4; i++) {
    if (ns[0][i] != nF[i])
      return false;
    for (int j = 0; j < nF[i]; j++)
      if (is[0][i][j] != indsF[i][j])
	return false;
  }
  return true;
}

Face4::Face4 (int nF[4], int indsF[4][3],
              int nM[4], int indsM[4][3]) {
  for (int i = 0; i < 4; i++) {
    ns[0][i] = nF[i];
    ns[1][i] = nM[i];
    for (int j = 0; j < 3; j++) {
      is[0][i][j] = indsF[i][j];
      is[1][i][j] = indsM[i][j];
    }
  }
}

Face4::Face4 (const Face4 &f, bool mf) {
  if (!mf) {
    *this = f;
    return;
  }
  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 4; i++) {
      ns[fm][i] = f.ns[!fm][i];
      for (int j = 0; j < 3; j++)
        is[fm][i][j] = f.is[!fm][i][j];
    }
}

void Face4::reindex () {
  int iNew[2][100];
  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < ns[fm][i]; j++) {
        assert(0 <= is[fm][i][j] && is[fm][i][j] < 100);
        iNew[fm][is[fm][i][j]] = -1;
      }
  for (int fm = 0; fm <= 1; fm++) {
    int n = 0;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < ns[fm][i]; j++)
        if (iNew[fm][is[fm][i][j]] == -1)
          iNew[fm][is[fm][i][j]] = n++;
  }
  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < ns[fm][i]; j++)
        is[fm][i][j] = iNew[fm][is[fm][i][j]];

  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 4; i++)
      for (int j0 = 0; j0 < ns[fm][i]; j0++)
	for (int j1 = j0+1; j1 < ns[fm][i]; j1++)
	  if (is[fm][i][j0] > is[fm][i][j1]) {
	    int temp = is[fm][i][j0];
	    is[fm][i][j0] = is[fm][i][j1];
	    is[fm][i][j1] = temp;
	  }  
}

void Face4::swap (int &a, int &b) {
  int temp = a;
  a = b;
  b = temp;
}

void Face4::swap (int a[3], int b[3]) {
  for (int i = 0; i < 3; i++) {
    int temp = a[i];
    a[i] = b[i];
    b[i] = temp;
  }
}

void Face4::swapFaces (int i, int j) {
  swap(ns[0][i], ns[0][j]);
  swap(ns[1][i], ns[1][j]);
  swap(is[0][i], is[0][j]);
  swap(is[1][i], is[1][j]);
}

void sortArray (int n, int *a) {
  for (int i = 0; i < n; i++)
    for (int j = i+1; j < n; j++)
      if (a[i] > a[j]) {
	int temp = a[i];
	a[i] = a[j];
	a[j] = temp;
      }
}

void Face4::replace (int r[2][12]) {
  for (int fm = 0; fm <= 1; fm++) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < ns[fm][i]; j++)
	is[fm][i][j] = r[fm][is[fm][i][j]];
      sortArray(ns[fm][i], is[fm][i]);
    }
  }
  for (int i = 0; i < 4 && ns[0][i] > 0; i++)
    for (int j = i+1; j < 4 && ns[0][j] > 0; j++) {
      if (ns[0][i] != ns[0][j] || ns[1][i] != ns[1][j])
	continue;
      int cmp = 0;
      for (int fm = 0; cmp == 0 && fm <= 1; fm++)
	for (int k = 0; cmp == 0 && k < ns[fm][i]; k++)
	  cmp = is[fm][i][k] - is[fm][j][k];
      if (cmp > 0)
	swapFaces(i, j);
    }
}

void Face4::replace2 (int r[2][12]) {
  for (int fm = 0; fm <= 1; fm++) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < ns[fm][i]; j++)
	is[fm][i][j] = r[fm][is[fm][i][j]];
      sortArray(ns[fm][i], is[fm][i]);
    }
  }
}

typedef int SWAPS[2];
int n0012 = 1, sw0012[][2] = {{0, 1}};
int n0112 = 1, sw0112[][2] = {{1, 2}};
int n0122 = 1, sw0122[][2] = {{2, 3}};
int n0011 = 3, sw0011[][2] = {{0, 1}, {2, 3}, {0, 1}};
int n0001 = 5, sw0001[][2] = {{0, 1}, {1, 2}, {0, 1}, {1, 2}, {0, 1}};
int n0111 = 5, sw0111[][2] = {{1, 2}, {2, 3}, {1, 2}, {2, 3}, {1, 2}};
int n0000 = 23, sw0000[][2] = {{0, 1}, {1, 2}, {0, 1}, {1, 2}, {0, 1}, {1, 3},
                               {0, 1}, {1, 2}, {0, 1}, {1, 2}, {0, 1}, {1, 3},
                               {0, 1}, {1, 2}, {0, 1}, {1, 2}, {0, 1}, {0, 3},
                               {0, 1}, {1, 2}, {0, 1}, {1, 2}, {0, 1}};
void Face4::minimize1 () {
  for (int i = 0; i < 4; i++)
    for (int j = i+1; j < 4; j++)
      if (ns[0][i] > ns[0][j])
        swapFaces(i, j);

  reindex();

  int values[4];
  values[0] = 0;
  for (int i = 1; i < 4; i++)
    values[i] = values[i-1] + (ns[0][i] > ns[0][i-1]);

  int pow10[] = { 1000, 100, 10, 1 };
  int rank[2][12];

  int values3;
  Face4 old = *this;
  do {
    values3 = values[3];
    old = *this;

    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < 12; i++)
        rank[fm][i] = 0;
    
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < ns[fm][i]; j++)
          rank[fm][is[fm][i][j]] += pow10[values[i]];
    
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < 4; i++)
        for (int j0 = 0; j0 < ns[fm][i]; j0++)
          for (int j1 = j0+1; j1 < ns[fm][i]; j1++)
            if (rank[fm][is[fm][i][j0]] < rank[fm][is[fm][i][j1]])
              swap(is[fm][i][j0], is[fm][i][j1]);
    
    for (int i0 = 0; i0 < 4; i0++)
      for (int i1 = i0+1; i1 < 4; i1++) {
        if (ns[0][i0] != ns[0][i1])
          continue;
        bool doSwap = false;
        for (int fm = 0; fm <= 1; fm++)
          for (int j = 0; j < ns[fm][i0]; j++)
            if (rank[fm][is[fm][i0][j]] > rank[fm][is[fm][i1][j]]) {
              fm = 2;
              break;
            }
            else if (rank[fm][is[fm][i0][j]] < rank[fm][is[fm][i1][j]]) {
              doSwap = true;
              fm = 2;
              break;
            }
        if (doSwap)
          swapFaces(i0, i1);
      }
    
    for (int i = 1; i < 4; i++) {
      if (ns[0][i] != ns[0][i-1]) {
        values[i] = values[i-1] + 1;
        continue;
      }
      bool greater = false;
      for (int fm = 0; fm <= 1; fm++)
        for (int j = 0; j < ns[fm][i]; j++)
          if (rank[fm][is[fm][i-1][j]] > rank[fm][is[fm][i][j]]) {
            greater = true;
            fm = 2;
            break;
          }
      values[i] = values[i-1] + greater;
    }

    reindex();
  } while (values[3] > values3 || !(*this == old));

  if (values[3] == 3)
    return;

  int nswaps;
  SWAPS *swaps;
  if (values[3] == 2) {
    if (values[0] == values[1]) {
      nswaps = n0012;
      swaps = sw0012;
    }
    else if (values[1] == values[2]) {
      nswaps = n0112;
      swaps = sw0112;
    }
    else if (values[2] == values[3]) {
      nswaps = n0122;
      swaps = sw0122;
    }
    else
      assert(0);
  }
  else if (values[3] == 1) {
    if (values[0] == values[1] && values[2] == values[3]) {
      nswaps = n0011;
      swaps = sw0011;
    }
    else if (values[0] == values[1]) {
      nswaps = n0001;
      swaps = sw0001;
    }
    else {
      nswaps = n0111;
      swaps = sw0111;
    }
  }
  else {
    nswaps = n0000;
    swaps = sw0000;
  }

  Face4 minsofar = *this;
  int values2[] = { 0, 1, 2, 3 };
  Compare compare;

  for (int i = 0; i < nswaps; i++) {
    swapFaces(swaps[i][0], swaps[i][1]);

    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < 12; i++)
        rank[fm][i] = 0;
    
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < ns[fm][i]; j++)
          rank[fm][is[fm][i][j]] += pow10[values2[i]];
    
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < 4; i++)
        for (int j0 = 0; j0 < ns[fm][i]; j0++)
          for (int j1 = j0+1; j1 < ns[fm][i]; j1++)
            if (rank[fm][is[fm][i][j0]] < rank[fm][is[fm][i][j1]])
              swap(is[fm][i][j0], is[fm][i][j1]);

    reindex();

    if (compare(*this, minsofar))
      minsofar = *this;
  }

  *this = minsofar;
  // if (values[3] < 3)
  // cout << "this is it" << endl;

#ifdef BLEEN
  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 12; i++)
      rank[fm][i] = 0;

  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < ns[fm][i]; j++)
        rank[fm][is[fm][i][j]] += pow10[ns[0][i]];

  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 4; i++)
      for (int j0 = 0; j0 < ns[fm][i]; j0++)
        for (int j1 = j0+1; j1 < ns[fm][i]; j1++)
          if (rank[fm][is[fm][i][j0]] < rank[fm][is[fm][i][j1]])
            swap(is[fm][i][j0], is[fm][i][j1]);

  reindex();

  for (int i0 = 0; i0 < 4; i0++)
    for (int i1 = i0+1; i1 < 4; i1++) {
      if (ns[0][i0] != ns[0][i1])
        continue;
      bool doSwap = false;
      for (int fm = 0; fm <= 1; fm++)
        for (int j = 0; j < ns[fm][i0]; j++)
          if (rank[fm][is[fm][i0][j]] > rank[fm][is[fm][i1][j]]) {
            fm = 2;
            break;
          }
          else if (rank[fm][is[fm][i0][j]] < rank[fm][is[fm][i1][j]]) {
            doSwap = true;
            fm = 2;
            break;
          }
      if (doSwap)
        swapFaces(i0, i1);
    }

  reindex();

  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 12; i++)
      rank[fm][i] = 0;

  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < ns[fm][i]; j++)
        rank[fm][is[fm][i][j]] += pow10[i];

  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 4; i++)
      for (int j0 = 0; j0 < ns[fm][i]; j0++)
        for (int j1 = j0+1; j1 < ns[fm][i]; j1++)
          if (rank[fm][is[fm][i][j0]] < rank[fm][is[fm][i][j1]])
            swap(is[fm][i][j0], is[fm][i][j1]);

  reindex();

  for (int i0 = 0; i0 < 4; i0++)
    for (int i1 = i0+1; i1 < 4; i1++) {
      if (ns[0][i0] != ns[0][i1])
        continue;
      bool doSwap = false;
      for (int fm = 0; fm <= 1; fm++)
        for (int j = 0; j < ns[fm][i0]; j++)
          if (rank[fm][is[fm][i0][j]] > rank[fm][is[fm][i1][j]]) {
            fm = 2;
            break;
          }
          else if (rank[fm][is[fm][i0][j]] < rank[fm][is[fm][i1][j]]) {
            doSwap = true;
            fm = 2;
            break;
          }
      if (doSwap)
        swapFaces(i0, i1);
    }

  reindex();

  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 12; i++)
      rank[fm][i] = 0;

  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < ns[fm][i]; j++)
        rank[fm][is[fm][i][j]] += pow10[i];

  for (int fm = 0; fm <= 1; fm++)
    for (int i = 0; i < 4; i++)
      for (int j0 = 0; j0 < ns[fm][i]; j0++)
        for (int j1 = j0+1; j1 < ns[fm][i]; j1++)
          if (rank[fm][is[fm][i][j0]] < rank[fm][is[fm][i][j1]])
            swap(is[fm][i][j0], is[fm][i][j1]);

  reindex();
#endif
}

void Face4::minimize () {
  minimize1();
#ifdef BLEEN
  Face4 mf(*this, true);
  mf.minimize1();
  if (Compare()(mf, *this))
    *this = mf;
#endif
}

int numUnique (int iMax, int ns[4], int is[4][3]) {
  bool appears[iMax];
  for (int i = 0; i < iMax; i++)
    appears[i] = false;
  int n = 0;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < ns[i]; j++)
      if (!appears[is[i][j]]) {
        n++;
        appears[is[i][j]] = true;
      }
  return n;
}

int randomInt (int lo, int hi) {
  return lo + (random() % (hi - lo));
}

AnglePoly *C4D::getPoly (Face4 face4r) {
    int *nF = face4r.ns[0];
    int *nM = face4r.ns[1];
    typedef int int43[4][3];
    int43 &indsF = face4r.is[0];
    int43 &indsM = face4r.is[1];

    if (nF[0] == 0 && nM[0] == 3 &&
	indsM[0][0] == 3 &&
	indsM[0][1] == 8 &&
	indsM[0][2] == 9 &&
	nF[1] == 2 && nM[1] == 2 &&
	indsF[1][0] == 3 &&
	indsF[1][1] == 27 &&
	indsM[1][0] == 9 &&
	indsM[1][1] == 24 &&
	nF[2] == 2 && nM[2] == 2 &&
	indsF[2][0] == 27 &&
	indsF[2][1] == 34 &&
	indsM[2][0] == 9 &&
	indsM[2][1] == 10)
      cout << "this is it getPoly" << endl;

    // Extraneous Polynomials
#ifdef BLEEN
    for (int i = 0; i < 4; i++)
      if (nF[i] == 0 && nM[i] != 0)
	indsF[i][nF[i]++] = 0;
      else if (nF[i] != 0 && nM[i] == 0)
	indsM[i][nM[i]++] = 0;
#endif

    // EdgeAB *ab = 0;
    VertXYZ *ab[2] = { 0, 0 };
    // FaceABC *cde = 0;
    VertXYZ *cde[3] = { 0, 0, 0 };
    VertXYZ *edgesV[3][2] { { 0, 0 }, { 0, 0 }, { 0, 0 } };
    SumVV *vert = 0;
    SumEdge *edges[4] = { 0, 0, 0, 0 };
    SumFace *faces[4] = { 0, 0, 0, 0 };
    AnglePoly *poly = 0;

    if (nF[0] == 2 && nM[0] == 3) {
      // ab = i2edge[0][indsF[0]];
      ab[0] = verts[0][indsF[0][0]];
      ab[1] = verts[0][indsF[0][1]];
      // cde = i3face[1][indsM[0]];
      cde[0] = verts[1][indsM[0][0]];
      cde[1] = verts[1][indsM[0][1]];
      cde[2] = verts[1][indsM[0][2]];
      return new PolyEF(ab[0], ab[1], cde[0], cde[1], cde[2]);
    }

    if (nF[0] == 3 && nM[0] == 2) {
      // ab = i2edge[1][indsM[0]];
      ab[0] = verts[1][indsM[0][0]];
      ab[1] = verts[1][indsM[0][1]];
      // cde = i3face[0][indsF[0]];
      cde[0] = verts[0][indsF[0][0]];
      cde[1] = verts[0][indsF[0][1]];
      cde[2] = verts[0][indsF[0][2]];
      return new PolyEF(ab[0], ab[1], cde[0], cde[1], cde[2]);
    }

    if (nF[2] + nM[2] == 2 && (nF[2] == 0 || nM[2] == 0)) {
      VertXYZ *e[3][2];
      for (int i = 0; i < 3; i++)
	if (nM[i] == 0) {
	  e[i][0] = verts[0][indsF[i][0]];
	  e[i][1] = verts[0][indsF[i][1]];
	}
	else {
	  e[i][0] = verts[1][indsM[i][0]];
	  e[i][1] = verts[1][indsM[i][1]];
	}
      return new PolyEEE(e);
    }
	  
    if (nF[0] == 0 || nM[0] == 0 ||
	(nF[0] + nM[0] == 4 && nF[3] + nM[3] == 0)) {
      VertXYZ *v[3][4];
      for (int i = 0; i < 3; i++)
	if (nM[i] == 0) {
	  for (int j = 0; j < 3; j++)
	    v[i][j] = verts[0][indsF[i][j]];
	  v[i][3] = 0;
	}
	else if (nF[i] == 0) {
	  for (int j = 0; j < 3; j++)
	    v[i][j] = verts[1][indsM[i][j]];
	  v[i][3] = 0;
	}
	else {
	  v[i][0] = verts[0][indsF[i][0]];
	  v[i][1] = verts[0][indsF[i][1]];
	  v[i][2] = verts[1][indsM[i][0]];
	  v[i][3] = verts[1][indsM[i][1]];
	}
      return new PolyFFF(v);
    }    

    if (nF[0] == 1 && nM[0] == 1) {
      VertXYZ *vf = verts[0][indsF[0][0]];
      VertXYZ *vm = verts[1][indsM[0][0]];
      vert = SumVV::make(vf, vm);
    }

    for (int i = 0; i < 4 && nF[i] != 0; i++)
      if (nF[i] + nM[i] == 3) {
        if (nF[i] == 1) {
          VertXYZ *a = verts[0][indsF[i][0]];
          EdgeAB *bc = EdgeAB::make(this, 1, indsM[i]);
          edges[i] = SumVE::make(a, bc);
        }
        else {
          EdgeAB *ab = EdgeAB::make(this, 0, indsF[i]);
          VertXYZ *c = verts[1][indsM[i][0]];
          edges[i] = SumEV::make(ab, c);
        }
      }
      else if (nF[i] + nM[i] == 4) {
        if (nF[i] == 1) {
          VertXYZ *a = verts[0][indsF[i][0]];
          FaceABC *bcd = FaceABC::make(this, 1, indsM[i]);
          faces[i] = SumVF::make(a, bcd);
        }
        else if (nF[i] == 2) {
          EdgeAB *ab = EdgeAB::make(this, 0, indsF[i]);
          EdgeAB *cd = EdgeAB::make(this, 1, indsM[i]);
          faces[i] = SumEE::make(ab, cd);
        }
        else {
          FaceABC *abc = FaceABC::make(this, 0, indsF[i]);
          VertXYZ *d = verts[1][indsM[i][0]];
          faces[i] = SumFV::make(abc, d);
        }
      }
    
    if (nF[2] == 0) {
      if (nF[0] + nM[0] == 2)
        poly = new PolyVF(vert, faces[1]);
      else
        poly = new PolyEE(edges[0], edges[1]);
    }
    else if (nF[3] == 0)
      poly = new PolyEFF(edges[0], faces[1], faces[2]);
    else
      poly = new PolyFFFF(faces[0], faces[1], faces[2], faces[3]);

    return poly;
}

#ifdef BLEEN
void C4D::getContact (VertXYZ *v, int i, Contact &contact) {
  int fm = v->isMoving();
  contact.is[fm][i][contact.ns[fm][i]++] = v->id;
}

void C4D::getContact (EdgeAB *e, int i, Contact &contact) {
  getContact(e->getTail(), i, contact);
  getContact(e->getHead(), i, contact);
}

void C4D::getContact (FaceABC *f, int i, Contact &contact) {
  for (int j = 0; j < 3; j++)
    getContact(f->getVert(j), i, contact);
}

void C4D::getContact (SumVV *v, int i, Contact &contact) {
  getContact(v->vv[0], i, contact);
  getContact(v->vv[1], i, contact);
}

void C4D::getContact (SumEdge *e, int i, Contact &contact) {
  if (e->getType() == SVE) {
    SumVE *sve = dynamic_cast<SumVE*>(e);
    getContact(sve->v, i, contact);
    getContact(sve->e, i, contact);
  }
  else {
    assert(e->getType() == SEV);
    SumEV *sev = dynamic_cast<SumEV*>(e);
    getContact(sev->v, i, contact);
    getContact(sev->e, i, contact);
  }
}

void C4D::getContact (SumFace *f, int i, Contact &contact) {
  if (f->getType() == SVF) {
    SumVF *svf = dynamic_cast<SumVF*>(f);
    getContact(svf->v, i, contact);
    getContact(svf->f, i, contact);
  }
  else if (f->getType() == SFV) {
    SumFV *sfv = dynamic_cast<SumFV*>(f);
    getContact(sfv->v, i, contact);
    getContact(sfv->f, i, contact);
  }
  else {
    assert(f->getType() == SEE);
    SumEE *see = dynamic_cast<SumEE*>(f);
    getContact(see->ee[0], i, contact);
    getContact(see->ee[1], i, contact);
  }
}

Contact C4D::getContact (AnglePoly *poly) {
  Contact contact;
  switch (poly->getType()) {
  case AnglePoly::EF2:
    {
      PolyEF *polyEF = dynamic_cast<PolyEF*>(poly);
      for (int i = 0; i < 2; i++)
        getContact(polyEF->e[i], 0, contact);
      for (int i = 0; i < 3; i++)
        getContact(polyEF->f[i], 0, contact);
    }
    break;
  case AnglePoly::VF:
    {
      PolyVF *polyVF = dynamic_cast<PolyVF*>(poly);
      getContact(polyVF->v, 0, contact);
      getContact(polyVF->f, 1, contact);
    }
    break;
  case AnglePoly::EE:
    {
      PolyEE *polyEE = dynamic_cast<PolyEE*>(poly);
      getContact(polyEE->ee[0], 0, contact);
      getContact(polyEE->ee[1], 1, contact);
    }
    break;
  case AnglePoly::EFF:
    {
      PolyEFF *polyEFF = dynamic_cast<PolyEFF*>(poly);
      getContact(polyEFF->e, 0, contact);
      getContact(polyEFF->f1, 1, contact);
      getContact(polyEFF->f2, 2, contact);
    }
    break;
  case AnglePoly::FFFF:
    {
      PolyFFFF *polyFFFF = dynamic_cast<PolyFFFF*>(poly);
      for (int i = 0; i < 4; i++)
        getContact(polyFFFF->f[i], i, contact);
    }
  default:
    assert(0);
  }
  contact.normalize();
  return contact;
}
#endif

template<class N>
N C4D::ContactSign::calculate () { 
  N value = poly->value<N>(root);

  int svalue = value.sign(false);
  if (svalue != 0)
    return value;

  bool flag = true;
  N s = root->getS<N>(flag);
  Poly<N> f = root->getPoly()->univariate<N>();
  Poly<N> g = poly->univariate<N>();

  Poly<N> fder = f.der();
  Poly<N> f2 = f.gcd(fder);
  if (f2.degree() > 0) {
    Poly<N> r;
    f = f.quotient(f2, r);
    assert(r.degree() == -1);
  }
  Poly<N> fg = f.gcd(g);
  if (fg.degree() > 0) {
    Poly<N> r;
    f = f.quotient(fg, r);
    assert(r.degree() == -1);
  }
  N fs = f.value(s);

  if (fs.sign(false) != 0)
    return N(0);

  value.sign();
  assert(0);
  return N(0);
}

bool C4D::commonRoot (Contact &contact, RootAngle *root) {
#ifdef TABLE_IDENTITY
  AnglePoly *poly = root->getPoly();
  Contact factor = getContact(poly);
  vector<Contact> factors = contact.getFactors();
  for (int i = 0; i < factors.size(); i++)
    if (factors[i] == factor)
      return true;
  return false;
#else
  PTR<AnglePoly> poly = getPoly(contact);
  return ContactSign(poly, root) == 0;
#ifdef BLEEN
  AnglePoly *poly = getPoly(contact);
  Parameter value = poly->value(root);

  if (value.sign(false) != 0)
    return false;

  bool flag = false;
  Parameter s = root->getS(flag);
  Poly1 f = root->getPoly()->univariate(flag);
  Poly1 g = poly->univariate(flag);
  delete poly;

  Poly1 fder = f.der();
  Poly1 f2 = f.gcd(fder);
  if (f2.degree() > 0) {
    Poly1 r;
    f = f.quotient(f2, r);
    assert(r.degree() == -1);
  }
  Poly1 fg = f.gcd(g);
  if (fg.degree() > 0) {
    Poly1 r;
    f = f.quotient(fg, r);
    assert(r.degree() == -1);
  }
  Parameter fs = f.value(s);

  if (fs.sign(false) != 0)
    return true;

  value.sign();
  assert(0);
  return true;
#endif
#endif
}

vector<C4D::Event*> C4D::solve (Contact &contact) {
  static int count = 0;
  ++count;
  vector<Contact> factors = contact.getFactors();
#ifndef TABLE_IDENTITY
    if (factors.size() > 1) {
      factors.clear();
      factors.push_back(contact);
    }
#endif    
  vector<Event*> events;
  for (int i = 0; i < factors.size(); i++) {
    Contact &factor = factors[i];
    if (contact2events.find(factor) == contact2events.end()) {
      AnglePoly *poly = getPoly(factor);
      vector<PTR<Angle>> angles = poly->getAnglesSorted();
      vector<Event*> fevents;
      for (int j = 0; j < angles.size(); j++) {
        Event *event = new Event(this, angles[j]);
	if (event->lessThan(currentEvent)) {
	  delete event;
	  continue;
	}
        auto it = eventList.find(event);
        if (it == eventList.end())
          eventList.insert(event);
        else {
          delete event;
          event = *it;
        }
        fevents.push_back(event);
      }
      contact2events[factor] = fevents;
    }
    vector<Event*> &fevents = contact2events[factor];
    for (int j = 0; j < fevents.size(); j++)
      events.push_back(fevents[j]);
  }
  return events;
}

class Compare {
public:
  bool operator() (const vector<int> &a, const vector<int> &b) const {
    if (a.size() != b.size())
      return a.size() < b.size();
    for (int i = 0; i < a.size(); i++)
      if (a[i] != b[i])
	return a[i] < b[i];
    return false;
  }
};

int alts[3334];

class Alists {
public:
  int n, nF, nM;
  int a[16];

  Alists (int n, int nF, int nM) : n(n), nF(nF), nM(nM) {
    for (int i = 0; i < 16; i++)
      a[i] = 0;
  }

  void setF (int i, int b) { a[i] |= (1 << (4 - 1 - b)); }
  void setM (int i, int b) { a[nF + i] |= (1 << (4 - 1 - b)); }

  bool getF (int i, int b) { return (a[i] & (1 << (4 - 1 - b))) != 0; }
  bool getM (int i, int b) { return (a[nF + i] & (1 << (4 - 1 - b))) != 0; }

  static void sort (int *a, int l) {
    for (int i = 1; i < l; i++) {
      int c = a[i];
      int j = i;
      while (--j >= 0 && a[j] < c)
        a[j+1] = a[j];
      a[j+1] = c;
    }
  }      

  static void sort (int *a, int l, int perm[12]) {
    for (int i = 1; i < l; i++) {
      int c = a[i];
      int c2 = perm[i];
      int j = i;
      while (--j >= 0 && a[j] < c) {
        a[j+1] = a[j];
	perm[j+1] = perm[j];
      }
      a[j+1] = c;
      perm[j+1] = c2;
    }
  }      

  void sort () {
    sort(a, nF);
    sort(a + nF, nM);
  }

  void sort (int perm[2][12]) {
    sort(a, nF, perm[0]);
    sort(a + nF, nM, perm[1]);
  }

#ifdef BLEEN
  void swap (int i, int j) {
    int *swap = swapBits[i][j];
    for (int k = 0; k < nF + nM; k++)
      a[k] = swap[a[k]];
    sort();
  }
#endif

  bool operator< (const Alists &that) const {
    for (int i = 0; i < 16; i++)
      if (a[i] != that.a[i])
        return a[i] > that.a[i];
    return false;
  }

#ifdef BLEEN  
  Alists minimum (int nSwaps, SWAPS *swaps, Perm &pMin) {
    Alists alMin = *this;
    Perm perm;
    for (int i = 0; i < nSwaps; i++) {
      swap(swaps[i][0], swaps[i][1]);
      perm.swap(swaps[i][0], swaps[i][1]);
      if (*this < alMin) {
        alMin = *this;
        pMin = perm;
        // cout << " swapped ";
      }
    }
    return alMin;
  }
#endif

  void print (ostream &out) {
    for (int i = 0; i < nF + nM; i++) {
      if (i != 0)
        out << " ";
      for (int j = 0; j < n; j++)
        if ((a[i] & (1 << (4 - j - 1))) != 0)
          out << j;
    }
  }

  class Compare {
  public:
    bool operator() (const Alists &a, const Alists &b) const {
      if (a.n != b.n)
        return a.n < b.n;
      if (a.nF != b.nF)
        return a.nF < b.nF;
      if (a.nM != b.nM)
        return a.nM < b.nM;
      return a < b;
    }
  };
};

class Graph;

map<Alists, Graph*, Alists::Compare> alist2graph;

#ifdef BLEEN
class Graph {
public:
  int n;
  int nF[4];
  int nM[4];
  int f[4][4];
  int m[4][4];

  Graph () {}

  Graph (const Face4 &face4) {
    n = 0;
    for (int i = 0; i < 4; i++) {
      if (face4.ns[0][i] > 0) {
	n++;
	nF[i] = face4.ns[0][i];
	nM[i] = face4.ns[1][i];
      }
      else {
      	nF[i] = 0;
	nM[i] = 0;
      }
      for (int j = 0; j < 4; j++) {
	f[i][j] = -1;
	m[i][j] = -1;
      }
    }

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nF[i]; j++)
	f[i][j] = face4.is[0][i][j];
      for (int j = 0; j < nM[i]; j++)
	m[i][j] = face4.is[1][i][j];
    }
  }

  Face4 getFace4 () {
    for (int i = n; i < 4; i++) {
      nF[i] = nM[i] = 0;
    }

    int indsF[4][3];
    int indsM[4][3];
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
	indsF[i][j] = indsM[i][j] = -1;
    
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nF[i]; j++)
	indsF[i][j] = f[i][j];
      for (int j = 0; j < nM[i]; j++)
	indsM[i][j] = m[i][j];
    }

    Face4 face4(nF, indsF, nM, indsM);
    return face4;
  }

  Graph (int n, int nF[4], int nM[4]) : n(n) {
    for (int i = 0; i < 4; i++) {
      this->nF[i] = 0;
      this->nM[i] = 0;
      for (int j = 0; j < 4; j++) {
	this->f[i][j] = -1;
	this->m[i][j] = -1;
      }
    }

    for (int i = 0; i < n; i++) {
      this->nF[i] = nF[i];
      this->nM[i] = nM[i];
    }
  }

  bool equal (int a, int b) {
    if (nF[a] != nF[b])
      return false;
    if (nM[a] != nM[b])
      return false;
    for (int i = 0; i < nF[a]; i++)
      if (f[a][i] != f[b][i])
        return false;
    for (int i = 0; i < nM[a]; i++)
      if (m[a][i] != m[b][i])
        return false;
    return true;
  }

  int getNFV () {
    int nFV = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < nF[i]; j++)
        if (nFV <= f[i][j])
          nFV = f[i][j] + 1;
    return nFV;
  }
    
  int getNMV () {
    int nMV = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < nM[i]; j++)
        if (nMV <= m[i][j])
          nMV = m[i][j] + 1;
    return nMV;
  }
    
  void sort(int *list, int l) {
    for (int i = 1; i < l; i++) {
      int c = list[i];
      int j = i;
      while (--j >= 0 && list[j] < c)
        list[j+1] = list[j];
      list[j+1] = c;
    }
  }      

  Alists getLists () {
    Alists lists(n, getNFV(), getNMV());
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nF[i]; j++)
        lists.setF(f[i][j], i);
      for (int j = 0; j < nM[i]; j++)
        lists.setM(m[i][j], i);
    }
    // lists.sort();
    return lists;
  }

  Alists getLists (int perm[2][12]) {
    Alists lists(n, getNFV(), getNMV());
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nF[i]; j++)
        lists.setF(f[i][j], i);
      for (int j = 0; j < nM[i]; j++)
        lists.setM(m[i][j], i);
    }
    lists.sort(perm);
    return lists;
  }

  int aLists (int *lists) {
    int *fLists = lists;
    int nFV = getNFV();
    for (int v = 0; v < nFV; v++)
      fLists[v] = 0;
    int bit = 1 << n;
    for (int i = 0; i < n; i++) {
      bit /= 2;
      for (int j = 0; j < nF[i]; j++)
        fLists[f[i][j]] |= bit;
    }
    sort(fLists, nFV);

    int *mLists = lists + nFV;
    int nMV = getNMV();
    for (int v = 0; v < nMV; v++)
      mLists[v] = 0;
    bit = 1 << n;
    for (int i = 0; i < n; i++) {
      bit /= 2;
      for (int j = 0; j < nM[i]; j++)
        mLists[m[i][j]] |= bit;
    }
    sort(mLists, nMV);

    return nFV + nMV;
  }

  void apply (Alists lists) {
    for (int i = 0; i < n; i++)
      nF[i] = nM[i] = 0;
    for (int v = 0; v < lists.nF; v++)
      for (int i = 0; i < n; i++)
        if (lists.getF(v, i))
          f[i][nF[i]++] = v;
    for (int v = 0; v < lists.nM; v++)
      for (int i = 0; i < n; i++)
        if (lists.getM(v, i))
          m[i][nM[i]++] = v;
    Alists lists2 = getLists();
    assert(!(lists < lists2) && !(lists2 < lists));
  }

#ifdef BLEEN
  Graph getMin (Perm &pMin) {
    int eqNext = 0;
    int bit = 1;
    for (int i = 1; i < n; i++) {
      if (nF[i-1] == nF[i] && nM[i-1] == nM[i])
        eqNext |= bit;
      bit *= 2;
    }
    Alists lists = getLists();
    Alists aMin = lists.minimum(nswapss[eqNext], swapss[eqNext], pMin);

    Graph gMin = *this;
    gMin.apply(aMin);
    return gMin;
  }
#endif

  void print (ostream &out) {
    if (true) {
      out << n << "   ";
      for (int i = 0; i < n; i++) {
        out << " " << nF[i];
        for (int j = 0; j < nF[i]; j++)
          out << " " << f[i][j];
        out << "  " << nM[i];
        for (int j = 0; j < nM[i]; j++)
          out << " " << m[i][j];
	out << "  ";
      }
    }
    else {
      for (int i = 0; i < n; i++) {
        if (i != 0)
          out << " ";
        for (int j = 0; j < nF[i]; j++)
          out << f[i][j];
        out << "+";
        for (int j = 0; j < nM[i]; j++)
          out << m[i][j];
      }
    }
#ifdef BLEEN
    int lists[16];
    int l = aLists(lists);
    for (int i = 0; i < l; i++) {
      out << " ";
      for (int j = 0; j < n; j++)
        if ((lists[i] & (1 << (4 - j - 1))) != 0)
          out << j;
    }
    Alists lists = getLists();

    if (false) {
      int l = getNFV() + getNMV();
      for (int i = 0; i < l; i++) {
        out << " ";
        for (int j = 0; j < n; j++)
          if ((lists.a[i] & (1 << (4 - j - 1))) != 0)
            out << j;
      }
    }

    Graph gmin = *this;
    gmin.apply(lists);
    gmin.printFM(out);

    out << endl;
#endif
  }

  bool parallelFaces () {
    for (int i = 0; i < n; i++)
      for (int j = i+1; j < n; j++)
	if ((nF[i] == 3 && nF[j] == 3 &&
	     f[i][0] == f[j][0] &&
	     f[i][1] == f[j][1] &&
	     f[i][2] == f[j][2]) ||
	    (nM[i] == 3 && nM[j] == 3 &&
	     m[i][0] == m[j][0] &&
	     m[i][1] == m[j][1] &&
	     m[i][2] == m[j][2]))
	  return true;
    return false;
  }
  
#ifdef BLEEN
  void gen (bool fm, int nV, int i, int j) {
    if (!fm) {
      if (i < n) {
        if (j < nF[i]) {
          int s = j == 0 ? 0 : f[i][j-1]+1;
          for (int v = s; v < nV; v++) {
            f[i][j] = v;
            gen(fm, nV, i, j+1);
          }
          f[i][j] = nV;
          gen(fm, nV+1, i, j+1);
        }
        else
          gen(fm, nV, i+1, 0);
      }
      else
        gen(true, 0, 0, 0);
    }
    else {
      if (i < n) {
        if (j < nM[i]) {
          int s = j == 0 ? 0 : m[i][j-1]+1;
          for (int v = s; v < nV; v++) {
            m[i][j] = v;
            gen(fm, nV, i, j+1);
          }
          m[i][j] = nV;
          gen(fm, nV+1, i, j+1);
        }
        else
          gen(fm, nV, i+1, 0);
      }
      else {
        for (int i = 0; i < n; i++)
          for (int j = i+1; j < n; j++)
            if (equal(i, j))
              return;

        static int count;
        // if (++count < 1006559)
        // return;
        // Alists lists = getLists();
        // apply(lists);
        // print(cout);
        // cout << " ";
        // lists.print(cout);
        Perm pMin;
        Graph gMin = getMin(pMin);
        gMin.print(cout);
        cout << endl;
        return;
      }
    }
  }

  void gen () {
    gen(false, 0, 0, 0);
  }
#endif
};
#endif

class Graph {
public:
  int n;
  int nF[4];
  int nM[4];
  int f[4][4];
  int m[4][4];

  bool isContact () {
    if (n == 1)
      return false;
    for (int i = 0; i < n; i++)
      if (nF[i] == 0 || nM[i] == 0)
	return false;
    if (n == 3 && nF[0] + nM[0] == nF[1] + nM[1] && nF[0] + nM[0] == nF[2] + nM[2])
      return false;
    return true;
  }

  long key () {
    long s = 1;
    long k = 0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nF[i]; j++) {
	k += s * f[i][j];
	s *= 16;
      }
      for (int j = 0; j < nM[i]; j++) {
	k += s * m[i][j];
	s *= 16;
      }
    }
    return k;
  }

  int compareC (int i, int j) {
    if (nF[i] + nM[i] != nF[j] + nM[j])
      return (nF[i] + nM[i]) - (nF[j] + nM[j]);
    if (nF[i] != nF[j])
      return nF[i] - nF[j];
    for (int k = 0; k < nF[i]; k++)
      if (f[i][k] != f[j][k])
        return f[i][k] - f[j][k];
    for (int k = 0; k < nM[i]; k++)
      if (m[i][k] != m[j][k])
        return m[i][k] - m[j][k];
    return 0;
  }

  void swap (int &a, int &b) {
    int tmp = a;
    a = b;
    b = tmp;
  }

  void swapC (int i, int j) {
    for (int k = 0; k < nF[i]; k++)
      swap(f[i][k], f[j][k]);
    for (int k = 0; k < nM[i]; k++)
      swap(m[i][k], m[j][k]);
  }

  void sort2 (int a[4], int n) {
    for (int i = 0; i < n; i++)
      for (int j = i+1; j < n; j++)
        if (a[i] > a[j])
          swap(a[i], a[j]);
  }

  void normalize () {
    for (int i = 0; i < 4; i++) {
      sort2(f[i], nF[i]);
      sort2(m[i], nM[i]);
    }
    for (int i = 0; i < 4; i++)
      for (int j = i+1; j < 4 && (nF[j] > 0 || nM[j] > 0); j++)
        if (compareC(i, j) > 0)
          swapC(i, j);
  }

  bool operator== (Graph &g) {
    for (int i = 0; i < 4; i++)
      if (nF[i] != g.nF[i] || nM[i] != g.nM[i])
        return false;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < nF[i]; j++)
        if (f[i][j] != g.f[i][j])
          return false;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < nM[i]; j++)
        if (m[i][j] != g.m[i][j])
          return false;
    return true;
  }

  Graph () {}

  Graph (const Face4 &face4) {
    n = 0;
    for (int i = 0; i < 4; i++) {
      if (face4.ns[0][i] > 0 || face4.ns[1][i] > 0) {
	n++;
	nF[i] = face4.ns[0][i];
	nM[i] = face4.ns[1][i];
      }
      else {
      	nF[i] = 0;
	nM[i] = 0;
      }
      for (int j = 0; j < 4; j++) {
	f[i][j] = -1;
	m[i][j] = -1;
      }
    }

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nF[i]; j++)
	f[i][j] = face4.is[0][i][j];
      for (int j = 0; j < nM[i]; j++)
	m[i][j] = face4.is[1][i][j];
    }
  }

  Face4 getFace4 () {
    for (int i = n; i < 4; i++) {
      nF[i] = nM[i] = 0;
    }

    int indsF[4][3];
    int indsM[4][3];
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++)
	indsF[i][j] = indsM[i][j] = -1;
    
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nF[i]; j++)
	indsF[i][j] = f[i][j];
      for (int j = 0; j < nM[i]; j++)
	indsM[i][j] = m[i][j];
    }

    Face4 face4(nF, indsF, nM, indsM);
    return face4;
  }

  Graph (int n, int nF[4], int nM[4]) : n(n) {
    for (int i = 0; i < 4; i++) {
      this->nF[i] = 0;
      this->nM[i] = 0;
      for (int j = 0; j < 4; j++) {
	this->f[i][j] = -1;
	this->m[i][j] = -1;
      }
    }

    for (int i = 0; i < n; i++) {
      this->nF[i] = nF[i];
      this->nM[i] = nM[i];
    }
  }

  bool equal (int a, int b) {
    if (nF[a] != nF[b])
      return false;
    if (nM[a] != nM[b])
      return false;
    for (int i = 0; i < nF[a]; i++)
      if (f[a][i] != f[b][i])
        return false;
    for (int i = 0; i < nM[a]; i++)
      if (m[a][i] != m[b][i])
        return false;
    return true;
  }

  int getNFV () {
    int nFV = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < nF[i]; j++)
        if (nFV <= f[i][j])
          nFV = f[i][j] + 1;
    return nFV;
  }
    
  int getNMV () {
    int nMV = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < nM[i]; j++)
        if (nMV <= m[i][j])
          nMV = m[i][j] + 1;
    return nMV;
  }
    
  void sort(int *list, int l) {
    for (int i = 1; i < l; i++) {
      int c = list[i];
      int j = i;
      while (--j >= 0 && list[j] < c)
        list[j+1] = list[j];
      list[j+1] = c;
    }
  }      

  Alists getLists () {
    Alists lists(n, getNFV(), getNMV());
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < nF[i]; j++)
        lists.setF(f[i][j], i);
      for (int j = 0; j < nM[i]; j++)
        lists.setM(m[i][j], i);
    }
    // lists.sort();
    return lists;
  }

  int aLists (int *lists) {
    int *fLists = lists;
    int nFV = getNFV();
    for (int v = 0; v < nFV; v++)
      fLists[v] = 0;
    int bit = 1 << n;
    for (int i = 0; i < n; i++) {
      bit /= 2;
      for (int j = 0; j < nF[i]; j++)
        fLists[f[i][j]] |= bit;
    }
    sort(fLists, nFV);

    int *mLists = lists + nFV;
    int nMV = getNMV();
    for (int v = 0; v < nMV; v++)
      mLists[v] = 0;
    bit = 1 << n;
    for (int i = 0; i < n; i++) {
      bit /= 2;
      for (int j = 0; j < nM[i]; j++)
        mLists[m[i][j]] |= bit;
    }
    sort(mLists, nMV);

    return nFV + nMV;
  }

  void apply (Alists lists) {
    for (int i = 0; i < n; i++)
      nF[i] = nM[i] = 0;
    for (int v = 0; v < lists.nF; v++)
      for (int i = 0; i < n; i++)
        if (lists.getF(v, i))
          f[i][nF[i]++] = v;
    for (int v = 0; v < lists.nM; v++)
      for (int i = 0; i < n; i++)
        if (lists.getM(v, i))
          m[i][nM[i]++] = v;
    Alists lists2 = getLists();
    assert(!(lists < lists2) && !(lists2 < lists));
  }

#ifdef BLEEN
  Graph getMin (Perm &pMin) {
    int eqNext = 0;
    int bit = 1;
    for (int i = 1; i < n; i++) {
      if (nF[i-1] == nF[i] && nM[i-1] == nM[i])
        eqNext |= bit;
      bit *= 2;
    }
    Alists lists = getLists();
    Alists aMin = lists.minimum(nswapss[eqNext], swapss[eqNext], pMin);

    Graph gMin = *this;
    gMin.apply(aMin);
    return gMin;
  }
#endif

  void print (ostream &out) {
    if (true) {
      out << n << "   ";
      for (int i = 0; i < n; i++) {
        out << " " << nF[i];
        for (int j = 0; j < nF[i]; j++)
          out << " " << f[i][j];
        out << "  " << nM[i];
        for (int j = 0; j < nM[i]; j++)
          out << " " << m[i][j];
	out << "  ";
      }
    }
    else {
      for (int i = 0; i < n; i++) {
        if (i != 0)
          out << " ";
        for (int j = 0; j < nF[i]; j++)
          out << f[i][j];
        out << "+";
        for (int j = 0; j < nM[i]; j++)
          out << m[i][j];
      }
    }
#ifdef BLEEN
    int lists[16];
    int l = aLists(lists);
    for (int i = 0; i < l; i++) {
      out << " ";
      for (int j = 0; j < n; j++)
        if ((lists[i] & (1 << (4 - j - 1))) != 0)
          out << j;
    }
    Alists lists = getLists();

    if (false) {
      int l = getNFV() + getNMV();
      for (int i = 0; i < l; i++) {
        out << " ";
        for (int j = 0; j < n; j++)
          if ((lists.a[i] & (1 << (4 - j - 1))) != 0)
            out << j;
      }
    }

    Graph gmin = *this;
    gmin.apply(lists);
    gmin.printFM(out);

    out << endl;
#endif
  }

  int depends (int i) {
    int d = 0;
    for (int j = 0; j < nF[i]; j++)
      d |= (1 << f[i][j]);
    for (int j = 0; j < nM[i]; j++)
      d |= (1 << (m[i][j] + 16));
    return d;
  }

  static int subset (int *a, int m, int *b, int n) {
    int count = 0;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
	if (a[i] == b[j])
	  count++;
    return count == m;
  }

  int parallelEdgeFace () {
    for (int i = 0; i < n; i++)
      if (nF[i] + nM[i] == 3)
	for (int j = i+1; j < n; j++)
	  if (nF[j] + nM[j] == 4 &&
	      ((nF[i] == 2 && subset(f[i], nF[i], f[j], nF[j])) ||
               (nM[i] == 2 && subset(m[i], nM[i], m[j], nM[j]))))
	    return depends(i) | depends(j);
    return 0;
  }

  int parallelFaces () {
    for (int i = 0; i < n; i++)
      for (int j = i+1; j < n; j++)
	if ((nF[i] == 3 && nF[j] == 3 &&
	     f[i][0] == f[j][0] &&
	     f[i][1] == f[j][1] &&
	     f[i][2] == f[j][2]) ||
	    (nM[i] == 3 && nM[j] == 3 &&
	     m[i][0] == m[j][0] &&
	     m[i][1] == m[j][1] &&
	     m[i][2] == m[j][2]))
	  return depends(i) | depends(j);
    return 0;
  }

  static int intersect (int *a, int m, int *b, int n, int *c) {
    int count = 0;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
	if (a[i] == b[j])
	  c[count++] = a[i];
    return count;
  }

  int parallelFaces3 () {
    int a[3];
    for (int i = 0; i < n; i++)
      if (nF[i] + nM[i] == 4) {
	for (int j = i+1; j < n; j++)
	  if (intersect(f[i], nF[i], f[j], nF[j], a) == 2) {
	    for (int k = j+1; k < n; k++)
	      if (subset(a, 2, f[k], nF[k]))
		return depends(i) | depends(j) | depends(k);
	  }
	  else
            if (intersect(m[i], nM[i], m[j], nM[j], a) == 2)
              for (int k = j+1; k < n; k++)
                if (subset(a, 2, m[k], nM[k]))
                  return depends(i) | depends(j) | depends(k);
      }
    return 0;
  }

  int commonEdge () {
    int a[3], b[3];
    for (int i = 0; i < n; i++)
      if (nF[i] + nM[i] == 4) {
	for (int j = i+1; j < n; j++)
	  if (intersect(f[i], nF[i], f[j], nF[j], a) == 2 &&
	      intersect(m[i], nM[i], m[j], nM[j], b) == 1) {
	    for (int k = 0; k < n; k++)
	      if (k != i && k != j && subset(a, 2, f[k], nF[k]))
		return ((1 << a[0]) | (1 << a[1]) | (1 << (b[0] + 16)) |
			depends(k));
	  }
	  else if (intersect(m[i], nM[i], m[j], nM[j], a) == 2 &&
		   intersect(f[i], nF[i], f[j], nF[j], b) == 1) {
	    for (int k = 0; k < n; k++)
	      if (k != i && k != j && subset(a, 2, m[k], nM[k]))
		return ((1 << (a[0] + 16)) | (1 << (a[1] + 16)) | (1 << b[0]) |
			depends(k));
	  }
      }
    return 0;
  }

  int parallelDepends () {
    int par;
    par = parallelEdgeFace();
    if (par != 0)
      return par;
    par = parallelFaces();
    if (par != 0)
      return par;
    par = parallelFaces3();
    if (par != 0) {
      int com = commonEdge();
      if (com != 0)
	return com;
      return par;
    }
    return 0;
    // return parallelEdgeFace() | parallelFaces() | parallelFaces3();
  }

#ifdef BLEEN
  void gen (bool fm, int nV, int i, int j) {
    if (!fm) {
      if (i < n) {
        if (j < nF[i]) {
          int s = j == 0 ? 0 : f[i][j-1]+1;
          for (int v = s; v < nV; v++) {
            f[i][j] = v;
            gen(fm, nV, i, j+1);
          }
          f[i][j] = nV;
          gen(fm, nV+1, i, j+1);
        }
        else
          gen(fm, nV, i+1, 0);
      }
      else
        gen(true, 0, 0, 0);
    }
    else {
      if (i < n) {
        if (j < nM[i]) {
          int s = j == 0 ? 0 : m[i][j-1]+1;
          for (int v = s; v < nV; v++) {
            m[i][j] = v;
            gen(fm, nV, i, j+1);
          }
          m[i][j] = nV;
          gen(fm, nV+1, i, j+1);
        }
        else
          gen(fm, nV, i+1, 0);
      }
      else {
        for (int i = 0; i < n; i++)
          for (int j = i+1; j < n; j++)
            if (equal(i, j))
              return;

        static int count;
        // if (++count < 1006559)
        // return;
        // Alists lists = getLists();
        // apply(lists);
        // print(cout);
        // cout << " ";
        // lists.print(cout);
        Perm pMin;
        Graph gMin = getMin(pMin);
        gMin.print(cout);
        cout << endl;
        return;
      }
    }
  }

  void gen () {
    gen(false, 0, 0, 0);
  }
#endif
};

#ifdef BLEEN
map<vector<int>, Graph, Compare> vec2graph;
C4D *c4d;

void createPerms (Alists &lists, Alists &perm, int fm, int v) {
  int nF = lists.nF;
  int nM = lists.nM;
  if (fm == 0) {
    if (v < nF) {
      int start = nF;
      while (start > 0 && perm.a[start-1] != lists.a[v])
	start--;
      for (int v2 = start; v2 < nF; v2++)
	if (perm.a[v2] == 0) {
	  perm.a[v2] = lists.a[v];
	  createPerms(lists, perm, fm, v+1);
	  perm.a[v2] = 0;
	}
    }
    else
      createPerms(lists, perm, 1, 0);
  }
  else {
    if (v < nM) {
      int start = nM;
      while (nF + start > 0 && perm.a[nF + start-1] != lists.a[nF + v])
	start--;
      for (int v2 = start; v2 < nM; v2++)
	if (perm.a[nF + v2] == 0) {
	  perm.a[nF + v2] = lists.a[nF + v];
	  createPerms(lists, perm, fm, v+1);
	  perm.a[nF + v2] = 0;
	}
    }
    else {
      int nF[4], nM[4];
      Graph graph(lists.n, nF, nM);
      graph.apply(perm);

      Face4 face4 = graph.getFace4();
      AnglePoly *poly = c4d->getPoly(face4);
      extern vector<int> coeffs26;
      vector<PTR<Angle>> roots = poly->getAnglesSorted();
      int coeffs[7] = { 0, 0, 0, 0, 0, 0, 0 };
      for (int i = 0; i < coeffs26.size(); i++)
	coeffs[i] = coeffs26[i];
      vec2graph[coeffs26] = graph;
    }
  }
}
#endif

#ifdef BLEEN
void read (istream &in, vector< vector<int> > &factors) {
  factors.clear();

  int nf;
  in >> nf;
  for (int i = 0; i < nf; i++) {
    int nc;
    in >> nc;
    vector<int> factor;
    for (int j = 0; j < nc; j++) {
      int c;
      in >> c;
      factor.push_back(c);
    }
    factors.push_back(factor);
  }
}
#endif

class Perm4 {
public:
  int perm[4];
};

int nGCD53, nGCD106, nGCD100, nGCD212, nGCD424, nRootCmps;

void C4D::contacts2polys () {
  vector<Face4> face4s;
  istream &in = std::cin;
  // in.open("../identity/min-class.txt");
  int n;
  while (in >> n) {
    int nF[4] = { 0, 0, 0, 0 };
    int nM[4] = { 0, 0, 0, 0 };
    int indsF[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
    int indsM[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
    
    for (int i = 0; i < n; i++) {
      in >> nF[i];
      for (int j = 0; j < nF[i]; j++)
	in >> indsF[i][j];

      in >> nM[i];
      for (int j = 0; j < nM[i]; j++)
	in >> indsM[i][j];
    }
    
    Face4 four(nF, indsF, nM, indsM);
    face4s.push_back(four);
  }

  int nFM[] = { 13, 13 };
  for (int fm = 0; fm <= 1; fm++) {
    for (int i = 0; i < nFM[fm]; i++)
      if (false)
	new VertXYZ(this,
		    randomNumber(-1, 1),
		    randomNumber(-1, 1),
		    randomNumber(-1, 1),
		    fm);
      else
	new VertXYZ(this,
#ifdef RANDOM10BITS
		    randomInt(0, 1024),
		    randomInt(0, 1024),
		    randomInt(0, 1024),
#else
		    randomInt(0, 1024*1024),
		    randomInt(0, 1024*1024),
		    randomInt(0, 1024*1024),
#endif
		    fm, true);

    Int2 i2;
    for (i2.i[0] = 0; i2.i[0] < nFM[fm]; i2.i[0]++)
      for (i2.i[1] = i2.i[0]+1; i2.i[1] < nFM[fm]; i2.i[1]++)
	i2edge[fm][i2] = new EdgeAB(verts[fm][i2.i[0]], verts[fm][i2.i[1]]);

    Int3 i3;
    for (i3.i[0] = 0; i3.i[0] < nFM[fm]; i3.i[0]++)
      for (i3.i[1] = i3.i[0]+1; i3.i[1] < nFM[fm]; i3.i[1]++)
	for (i3.i[2] = i3.i[1]+1; i3.i[2] < nFM[fm]; i3.i[2]++)
	  i3face[fm][i3] = new FaceABC(verts[fm][i3.i[0]], verts[fm][i3.i[1]], verts[fm][i3.i[2]]);
  }

  initSums(true);
  enable();
  assert(0); // need to make this function a friend
  // EParameter::curPrecision() = 848u;

#ifdef BLEEN
  // DEBUG
  int iPoly = 0;
  for (int iRep = 0; iRep < face4s.size(); iRep++) {
    Face4 face4r = face4s[iRep];

    AnglePoly *poly = getPoly(face4r);

    iPoly++;

    Poly2<N> f = poly->getPoly();
    Poly1 g = AnglePoly::univariate(f, true);
    for (int i = 0; i <= g.d; i++) {
      if (i > 0)
	cout << " ";
      g.a[i].print();
    }
    for (int i = g.d+1; i <= 6; i++)
      cout << " 0";

    if (iRep % 2 == 0)
      cout << "    ";
    else
      cout << endl;

    delete poly;
  }
#endif
}  

const int prime26a = 66986537;
const int prime26b = 66986519;
int prime26 = prime26a;

bool divisibleBy1plusS2 (vector<int> c) {
  vector<int> q(c.size()-2);
  for (int i = q.size()-1; i >= 0; i--) {
    q[i] = c[i+2];
    c[i] = (c[i] - c[i+2] + prime26) % prime26;
  }
  return(c[0] == 0 && c[1] == 0);
}

vector<int> divideBy1plusS2 (vector<int> c) {
  vector<int> q(c.size()-2);
  for (int i = q.size()-1; i >= 0; i--) {
    q[i] = c[i+2];
    c[i] = (c[i] - c[i+2] + prime26) % prime26;
  }
  assert(c[0] == 0 && c[1] == 0);
  return q;
}

map<vector<int>, Graph, Compare> vec2graph;
map<vector<int>, vector<Graph>, Compare> vec2graphs;
C4D *c4d;
Graph basic;
vector<int> basic_coeffs;
unordered_set<long> graph_keys;
bool really_basic;

long powerX (long a, long e, long p) {
  if (e == 0)
    return 1;
  if (e == 1)
    return a;
  long ae2 = powerX(a, e/2, p);
  if (e % 2 == 0)
    return (ae2 * ae2) % p;
  else
    return (((ae2 * ae2) % p) * a) % p;
}

long inverseX (long a, long p) {
  return powerX(a, p-2, p);
}

bool operator== (const vector<int> &c1, const vector<int> &c2) {
  if (c1.size() != c2.size())
    return false;
  for (int i = 0; i < c1.size(); i++)
    if (c1[i] != c2[i])
      return false;
  return true;
}

bool eqVector (const vector<int> &c1, const vector<int> &c2) {
  return c1 == c2;
}

vector<int> getCoeffs (Graph &graph, double prime=prime26) {
  double save = prime26;
  prime26 = prime;
  Face4 face4 = graph.getFace4();
  AnglePoly *poly = c4d->getPoly(face4);

  vector<int> coeffs26X;
  Poly<MParameter> g = poly->univariate<MParameter>();
  int l = g.deg();
  while (l >= 0 && g.a[l].lb() == 0)
    l--;
  long ainv = l != -1 ? inverseX(g.a[l].lb(), prime26) : 1;
  for (int i = 0; i <= g.deg(); i++) {
    long a = g.a[i].lb();
    a *= ainv;
    a %= prime26;
    if (a < 0)
      a += prime26;
    coeffs26X.push_back(a);
    // cout << " " << a;
  }

  delete poly;
  prime26 = save;

  vector<int> coeffs = coeffs26X;
  while (divisibleBy1plusS2(coeffs))
    coeffs = divideBy1plusS2(coeffs);
  return coeffs;
}

void createPerms (Alists &lists, Alists &perm, int fm, int v) {
  int nF = lists.nF;
  int nM = lists.nM;
  if (fm == 0) {
    if (v < nF) {
      int start = nF;
      while (start > 0 && perm.a[start-1] != lists.a[v])
	start--;
      for (int v2 = start; v2 < nF; v2++)
	if (perm.a[v2] == 0) {
	  perm.a[v2] = lists.a[v];
	  createPerms(lists, perm, fm, v+1);
	  perm.a[v2] = 0;
	}
    }
    else
      createPerms(lists, perm, 1, 0);
  }
  else {
    if (v < nM) {
      int start = nM;
      while (nF + start > 0 && perm.a[nF + start-1] != lists.a[nF + v])
	start--;
      for (int v2 = start; v2 < nM; v2++)
	if (perm.a[nF + v2] == 0) {
	  perm.a[nF + v2] = lists.a[nF + v];
	  createPerms(lists, perm, fm, v+1);
	  perm.a[nF + v2] = 0;
	}
    }
    else {
      int nF[4], nM[4];
      Graph graph(lists.n, nF, nM);
      graph.apply(perm);
      graph.normalize();

      long key = graph.key();
      if (graph_keys.find(key) != graph_keys.end())
	return;
      graph_keys.insert(key);
      assert(graph_keys.find(key) != graph_keys.end());

      vector<int> coeffs = getCoeffs(graph);
#ifdef BLEEN
      Face4 face4 = graph.getFace4();
      AnglePoly *poly = c4d->getPoly(face4);
      extern vector<int> coeffs26;
      vector<PTR<Angle>> roots = poly->getAnglesSorted();
      assert(roots.size() == 0);
      delete poly;
      extern vector<int> coeffs26;
      vector<int> coeffs;
      if (coeffs26.size() == basic_coeffs.size())
	coeffs = coeffs26;
      else if (coeffs26.size() == basic_coeffs.size()+2)
	coeffs = divideBy1plusS2(coeffs26);
      else
	assert(0);
#endif
      assert(coeffs.size() == basic_coeffs.size());

      bool eq_coeffs = (coeffs == basic_coeffs);
#ifdef BLEEN
      bool eq_coeffs = true;
      for (int i = 0; i < coeffs.size(); i++)
	if (coeffs[i] != basic_coeffs[i])
	  eq_coeffs = false;
#endif

      bool is_auto = false;
      if (eq_coeffs && !(graph == basic)) {
	cout << "basic ";
	graph.normalize();
	basic.print(cout);
	cout << endl;
	cout << "auto ";
	graph.print(cout);
	cout << endl;
        is_auto = true;

	Graph basic2 = basic;
	for (int i = 0; i < basic2.n; i++) {
	  for (int j = 0; j < basic2.nF[i]; j++)
	    basic2.f[i][j]++;
	  for (int j = 0; j < basic2.nM[i]; j++)
	    basic2.m[i][j]++;
	}
	vector<int> basic_coeffs2 = getCoeffs(basic2);

	Graph graph2 = graph;
	for (int i = 0; i < graph2.n; i++) {
	  for (int j = 0; j < graph2.nF[i]; j++)
	    graph2.f[i][j]++;
	  for (int j = 0; j < graph2.nM[i]; j++)
	    graph2.m[i][j]++;
	}
	vector<int> graph_coeffs2 = getCoeffs(graph2);

	assert(basic_coeffs2 == graph_coeffs2);
      }

      auto it = vec2graph.find(coeffs);
      // if (coeffs.size() < 7) {
      if (graph.n < 4) {
	if (it == vec2graph.end())
	  vec2graph[coeffs] = graph;
        else {
          auto its = vec2graphs.find(coeffs);
          if (its == vec2graphs.end()) {
            vector<Graph> graphs;
            graphs.push_back(it->second);
            graphs.push_back(graph);
            vec2graphs[coeffs] = graphs;
          }
          else
            its->second.push_back(graph);
        }
      }
      else if (is_auto) {
        auto its = vec2graphs.find(coeffs);
        if (its == vec2graphs.end()) {
          vector<Graph> graphs;
          graphs.push_back(basic);
          graphs.push_back(graph);
          vec2graphs[coeffs] = graphs;
        }
        else
          its->second.push_back(graph);
      }
#ifdef CHECK_MATCHES
      else {
	if (it != vec2graph.end()) {
	  vector<int> c1b = getCoeffs(graph, prime26b);
	  vector<int> c2b = getCoeffs(it->second, prime26b);
	  if (c1b == c2b) {
	    cout << "perm ";
	    it->second.print(cout);
	    cout << endl;
	    cout << "matches ";
	    graph.print(cout);
	    cout << endl;

	    really_basic = false;
	  }
	}
      }
#endif
    }
  }
}

void read (istream &in, vector< vector<int> > &factors) {
  factors.clear();

  int nf;
  in >> nf;
  for (int i = 0; i < nf; i++) {
    int nc;
    in >> nc;
    vector<int> factor;
    for (int j = 0; j < nc; j++) {
      int c;
      in >> c;
      factor.push_back(c);
    }
    factors.push_back(factor);
  }
}

int numOnes (int x) {
  int n = 0;
  while (x > 0) {
    if (x % 2 == 1)
      n++;
    x /= 2;
  }
  return n;
}

#ifdef GCD_IDENTITY
vector<Contact> gcdPairs;
#endif

C4D *current_c4d;
Contact current_contact;
Angle *current_angle;
int nIdTests, nIdentities;
bool commonRoot (Contact contact, RootAngle *root) {
  return current_c4d->commonRoot(contact, root);
}

void C4D::findIdentitiesT () {
  c4d = this;
  alts[1111] = 24;
  alts[1112] = 6;
  alts[1113] = 6;
  alts[1122] = 4;
  alts[1123] = 2;
  alts[1133] = 4;
  alts[1222] = 6;
  alts[1223] = 2;
  alts[1233] = 2;
  alts[1333] = 6;
  alts[2222] = 24;
  alts[2223] = 6;
  alts[2233] = 4;
  alts[2333] = 6;
  alts[3333] = 24;
  int nalts = 0;

  vector<Face4> face4s;
  ifstream in;
  in.open("min-class3.txt");

  int n;
  while (in >> n) {
    int nF[4] = { 0, 0, 0, 0 };
    int nM[4] = { 0, 0, 0, 0 };
    int indsF[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
    int indsM[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
    
    for (int i = 0; i < n; i++) {
      in >> nF[i];
      for (int j = 0; j < nF[i]; j++)
	in >> indsF[i][j];
      in >> nM[i];
      for (int j = 0; j < nM[i]; j++)
	in >> indsM[i][j];
    }
    
    Face4 four(nF, indsF, nM, indsM);
    face4s.push_back(four);
  }

#define NPOOL 12
  // int nFM[] = { 32, 32 };
  int nFM[] = { NPOOL, NPOOL };
  double x, y, z;
  bool duplicates = false;//true;
  for (int fm = 0; fm <= 1; fm++) {
    for (int i = 0; i < nFM[fm]; i++)
      if (true) {
	if (!duplicates || i % 2 == 0) {
	  x = randomNumber(-1, 1);
	  y = randomNumber(-1, 1);
	  z = randomNumber(-1, 1);
	}
	new VertXYZ(this, x, y, z, fm);
      }
      else {
	new VertXYZ(this,
		    randomInt(0, 1024),
		    randomInt(0, 1024),
		    randomInt(0, 1024),
		    fm, true);
      }
    Int2 i2;
    for (i2.i[0] = 0; i2.i[0] < nFM[fm]; i2.i[0]++)
      for (i2.i[1] = i2.i[0]+1; i2.i[1] < nFM[fm]; i2.i[1]++)
	i2edge[fm][i2] = new EdgeAB(verts[fm][i2.i[0]], verts[fm][i2.i[1]]);

    Int3 i3;
    for (i3.i[0] = 0; i3.i[0] < nFM[fm]; i3.i[0]++)
      for (i3.i[1] = i3.i[0]+1; i3.i[1] < nFM[fm]; i3.i[1]++)
	for (i3.i[2] = i3.i[1]+1; i3.i[2] < nFM[fm]; i3.i[2]++)
	  i3face[fm][i3] = new FaceABC(verts[fm][i3.i[0]], verts[fm][i3.i[1]], verts[fm][i3.i[2]]);
  }

  initSums(true);
  enable();
  Contact::readFactors();

  int nResolved = 0;

  current_c4d = this;

  for (int iTest = 0; iTest < 100000; iTest++) {
    //if (iTest % 1000 == 0)
    //cout << "iTest " << iTest << endl;
    // int iClass = 4564;
    int iClass = randomInt(0, face4s.size());
    int r[2][NPOOL];
    permutation(NPOOL, r[0]);
    permutation(NPOOL, r[1]);

    int r12[2][12];
    for (int i = 0; i < 12; i++) {
      r12[0][i] = r[0][i];
      r12[1][i] = r[1][i];
    }

    Face4 face4rep = face4s[iClass];
    Face4 face4 = face4rep;

    face4.replace(r12);
    Contact contact(face4.ns, face4.is);
    contact.normalize();
    current_contact = contact;
    vector<Contact> factors = contact.getFactors();

    if (factors.size() < 0) {
      cout << "contact " << iClass << " angles:" << endl;
      AnglePoly *poly = getPoly(contact);
      Poly2<Parameter> ppoly2 = poly->getApproxMid(1.0);
      vector<PTR<Angle>> angles = poly->getAnglesSorted();
      for (int j = 0; j < angles.size(); j++) {
	PV2<Parameter> xy = angles[j]->getApproxMid(1.0);
	cout << xy.x.lb() << " " << xy.x.ub() << " "
	     << xy.y.lb() << " " << xy.y.ub() << " "
	     << poly->value<Parameter>(angles[j]).lb() << " "
	     << poly->value<Parameter>(angles[j]).ub() << endl;
      }
      for (int i = 0; i < factors.size(); i++) {
	cout << "factor " << i << " angles:" << endl;
	AnglePoly *poly = getPoly(factors[i]);
	Poly2<Parameter> ppoly2 = poly->getApproxMid(1.0);
	vector<PTR<Angle>> angles = poly->getAnglesSorted();
	for (int j = 0; j < angles.size(); j++) {
	  PV2<Parameter> xy = angles[j]->getApproxMid(1.0);
	  cout << xy.x.lb() << " " << xy.x.ub() << " "
	       << xy.y.lb() << " " << xy.y.ub() << " "
	       << poly->value<Parameter>(angles[j]).lb() << " "
	       << poly->value<Parameter>(angles[j]).ub() << endl;
	}
      }
      cout << endl;
    }

#ifndef TABLE_IDENTITY
    if (factors.size() > 1) {
      factors.clear();
      factors.push_back(contact);
    }
#endif    

    for (int i = 0; i < factors.size(); i++) {
      Contact &factor = factors[i];
      if (contact2events.find(factor) == contact2events.end()) {
        AnglePoly *poly = getPoly(factor);
        vector<PTR<Angle>> angles = poly->getAnglesSorted();
        vector<Event*> fevents;
        for (int j = 0; j < angles.size(); j++) {
          Event *event = new Event(this, angles[j]);
          current_angle = angles[j];
          auto it = eventList.find(event);
          if (it == eventList.end())
            eventList.insert(event);
          else {
            delete event;
            event = *it;
          }
          fevents.push_back(event);
        }
        contact2events[factor] = fevents;
      }
      else
        nResolved++;
    }
  }

  cout << "nResolved " << nResolved << endl;
  // cout << "nPermsFound " << nPermsFound << endl;
  cout << "nGCD53 " << nGCD53 << endl;
  cout << "nGCD100 " << nGCD100 << endl;
  cout << "nGCD106 " << nGCD106 << endl;
  cout << "nGCD212 " << nGCD212 << endl;
  cout << "nGCD424 " << nGCD424 << endl;
  cout << "nIdTests " << nIdTests << endl;
  cout << "nIdentities " << nIdentities << endl;
  cout << "nRootCmps " << nRootCmps << endl;

#ifdef ROOT_SEPARATION2
  ofstream out;
  out.open("gcdpairs.txt");
  for (int i = 0; i < gcdPairs.size(); i++) {
    gcdPairs[i].write(out);
    out << endl;
  }
  out.close();
#endif
}

istream &operator>> (istream &in, Graph &graph) {
  int n;
  if (in >> n) {
    int nF[4] = { 0, 0, 0, 0 };
    int nM[4] = { 0, 0, 0, 0 };
    int indsF[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
    int indsM[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
    for (int i = 0; i < n; i++) {
      in >> nF[i];
      for (int j = 0; j < nF[i]; j++)
	in >> indsF[i][j];
      in >> nM[i];
      for (int j = 0; j < nM[i]; j++)
	in >> indsM[i][j];
    }
    
    Face4 four(nF, indsF, nM, indsM);
    graph = Graph(four);
  }
  return in;
}

vector<int> multiply (vector<int> &a, vector<int> &b) {
  vector<int> c(a.size() + b.size() - 1);
  for (int d = 0; d < c.size(); d++)
    c[d] = 0;
  for (int i = 0; i < a.size(); i++) {
    long ai = a[i];
    for (int j = 0; j < b.size(); j++) {
      long bj = b[j];
      c[i + j] = (c[i + j] + (ai * bj) % prime26) % prime26;
    }
  }
  return c;
}

void C4D::checkFactors () {
  c4d = this;

  int nFM[] = { 13, 13 };
  for (int fm = 0; fm <= 1; fm++) {
    for (int i = 0; i < nFM[fm]; i++)
      if (false)
	new VertXYZ(this,
		    randomNumber(-1, 1),
		    randomNumber(-1, 1),
		    randomNumber(-1, 1),
		    fm);
      else
	new VertXYZ(this,
		    randomInt(0, 1024*1024),
		    randomInt(0, 1024*1024),
		    randomInt(0, 1024*1024),
		    fm, true);

    Int2 i2;
    for (i2.i[0] = 0; i2.i[0] < nFM[fm]; i2.i[0]++)
      for (i2.i[1] = i2.i[0]+1; i2.i[1] < nFM[fm]; i2.i[1]++)
	i2edge[fm][i2] = new EdgeAB(verts[fm][i2.i[0]], verts[fm][i2.i[1]]);

    Int3 i3;
    for (i3.i[0] = 0; i3.i[0] < nFM[fm]; i3.i[0]++)
      for (i3.i[1] = i3.i[0]+1; i3.i[1] < nFM[fm]; i3.i[1]++)
	for (i3.i[2] = i3.i[1]+1; i3.i[2] < nFM[fm]; i3.i[2]++)
	  i3face[fm][i3] = new FaceABC(verts[fm][i3.i[0]], verts[fm][i3.i[1]], verts[fm][i3.i[2]]);
  }

  initSums(true);
  enable();

  ifstream in;
  in.open("min-class3.txt");
  ifstream in2;
  in2.open("eqsub-classes3.txt");

  Graph graph;
  while (in >> graph) {
    vector<int> coeffs = getCoeffs(graph);
    
    int m;
    in2 >> m;
    
    vector<int> product;
    for (int i = 0; i < m; i++) {
      Graph graph2;
      in2 >> graph2;
      if (i == 0)
        product = getCoeffs(graph2);
      else {
        vector<int> factor = getCoeffs(graph2);
        product = multiply(product, factor);
      }
    }

    for (int i = 0; i < coeffs.size() || i < product.size(); i++)
      if (i >= coeffs.size())
        assert(product[i] == 0);
      else if (i >= product.size())
        assert(coeffs[i] == 0);
      else
        assert(coeffs[i] == product[i]);
  }
}

void C4D::findIdentitiesModP () {
  c4d = this;
  alts[1111] = 24;
  alts[1112] = 6;
  alts[1113] = 6;
  alts[1122] = 4;
  alts[1123] = 2;
  alts[1133] = 4;
  alts[1222] = 6;
  alts[1223] = 2;
  alts[1233] = 2;
  alts[1333] = 6;
  alts[2222] = 24;
  alts[2223] = 6;
  alts[2233] = 4;
  alts[2333] = 6;
  alts[3333] = 24;
  int nalts = 0;

  vector<Face4> face4s;
  ifstream in;
  in.open("min-class3.txt");
  ifstream in2;
  in2.open("eq-factor3.txt");
  ifstream in3;
  in3.open("depends3.txt");
  ofstream out;
  out.open("eq-classes3.txt");
  int n;
  while (in >> n) {
    int nF[4] = { 0, 0, 0, 0 };
    int nM[4] = { 0, 0, 0, 0 };
    int indsF[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
    int indsM[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
    
    for (int i = 0; i < n; i++) {
      in >> nF[i];
      for (int j = 0; j < nF[i]; j++)
	in >> indsF[i][j];
      in >> nM[i];
      for (int j = 0; j < nM[i]; j++)
	in >> indsM[i][j];
    }
    
    Face4 four(nF, indsF, nM, indsM);
    face4s.push_back(four);
  }

  int nFM[] = { 13, 13 };
  for (int fm = 0; fm <= 1; fm++) {
    for (int i = 0; i < nFM[fm]; i++)
      if (false)
	new VertXYZ(this,
		    randomNumber(-1, 1),
		    randomNumber(-1, 1),
		    randomNumber(-1, 1),
		    fm);
      else
	new VertXYZ(this,
		    randomInt(0, 1024*1024),
		    randomInt(0, 1024*1024),
		    randomInt(0, 1024*1024),
		    fm, true);

    Int2 i2;
    for (i2.i[0] = 0; i2.i[0] < nFM[fm]; i2.i[0]++)
      for (i2.i[1] = i2.i[0]+1; i2.i[1] < nFM[fm]; i2.i[1]++)
	i2edge[fm][i2] = new EdgeAB(verts[fm][i2.i[0]], verts[fm][i2.i[1]]);

    Int3 i3;
    for (i3.i[0] = 0; i3.i[0] < nFM[fm]; i3.i[0]++)
      for (i3.i[1] = i3.i[0]+1; i3.i[1] < nFM[fm]; i3.i[1]++)
	for (i3.i[2] = i3.i[1]+1; i3.i[2] < nFM[fm]; i3.i[2]++)
	  i3face[fm][i3] = new FaceABC(verts[fm][i3.i[0]], verts[fm][i3.i[1]], verts[fm][i3.i[2]]);
  }

  initSums(true);
  enable();

  int iA = 0, iB = 0, iC = 0;

  int nF[] = { 0, 0, 0, 0 };
  int nM[] = { 0, 0, 0, 0 };
  Graph g(0, nF, nM);

  vector<int> bad;

#ifdef EDGE_FACE_PARALLEL
  bad.clear();
  bad.push_back(6480144);
  bad.push_back(56007641);
  bad.push_back(1);
  g.n = -1;
  vec2graph[bad] = g;

  bad.clear();
  bad.push_back(54772632);
  bad.push_back(8028415);
  bad.push_back(1);
  g.n = -2;
  vec2graph[bad] = g;

  bad.clear(); // 01 r12 r03 coplanar
  bad.push_back(60055967);
  bad.push_back(1319169);
  bad.push_back(1);
  g.n = -3;
  vec2graph[bad] = g;
#endif

#ifdef BLEEN
  bad.clear();
  bad.push_back(54803762);
  bad.push_back(61517002);
  bad.push_back(1);
  g.n = -4;
  vec2graph[bad] = g;

  bad.clear();
  bad.push_back(18493610);
  bad.push_back(56602397);
  bad.push_back(1);
  g.n = -4;
  vec2graph[bad] = g;

  bad.clear();
  bad.push_back(47403059);
  bad.push_back(65172823);
  bad.push_back(1);
  g.n = -5;
  vec2graph[bad] = g;

  bad.clear();
  bad.push_back(4273458);
  bad.push_back(13588952);
  bad.push_back(1);
  g.n = -6;
  vec2graph[bad] = g;

  bad.clear();
  bad.push_back(49122595);
  bad.push_back(10608799);
  bad.push_back(1);
  g.n = -7;
  vec2graph[bad] = g;

  bad.clear();
  bad.push_back(43497994);
  bad.push_back(27821812);
  bad.push_back(1);
  g.n = -8;
  vec2graph[bad] = g;

  bad.clear();
  bad.push_back(58819928);
  bad.push_back(55273705);
  bad.push_back(1);
  g.n = -9;
  vec2graph[bad] = g;

  bad.clear();
  bad.push_back(61242463);
  bad.push_back(272065);
  bad.push_back(1);
  g.n = -10;
  vec2graph[bad] = g;

  bad.clear();
  bad.push_back(54339061);
  bad.push_back(25420561);
  bad.push_back(1);
  g.n = -9;
  vec2graph[bad] = g;
#endif

  // DEBUG
  int iPoly = 0;
  for (int iRep = 0 /* 10319 */; iRep < face4s.size(); iRep++) {
    if (iRep % 100 == 0)
      cout << "iRep " << iRep << endl;
    Face4 face4r = face4s[iRep];

    Graph graph(face4r);
      
    vector< vector<int> > factors;
    read(in2, factors);

    int nfactors = factors.size();
    out << factors.size() << "  ";

    int par = 0;
    /*
    if (nfactors > 1)
      par = graph.parallelDepends();
    */

    int ndepends;
    in3 >> ndepends;
    assert(ndepends == nfactors);

    int nfound = 0, nContactFactors = 0;
    for (int i = 0; i < factors.size(); i++) {
      int fdep, mdep;
      in3 >> fdep >> mdep;
      int dep = fdep + (mdep << 16);
      if (par != 0 && numOnes(dep) > 5 && (par != dep)) {
	out << " " << -4 << " ";
	continue;
      }

      map<vector<int>, Graph, Compare>::iterator it =
	vec2graph.find(factors[i]);
      if (it != vec2graph.end()) {
	nfound++;

        auto its = vec2graphs.find(factors[i]);
        if (its == vec2graphs.end()) {
          out << " ";
          it->second.print(out);
	  if (it->second.isContact())
	    nContactFactors++;
        }
        else {
          out << " -999 " << its->second.size();
          for (int i = 0; i < its->second.size(); i++) {
            out << " ";
            its->second[i].print(out);
          }
        }
	// if (factors.size() == 1)
	// assert(getCoeffs(graph, prime26b) == getCoeffs(it->second, prime26b));
      }
      else
	assert(factors.size() == 1);
    }
      
    int dummy;
    if (ndepends == 2)
      for (int i = 0; i < 2; i++)
	in3 >> dummy;
    else if (ndepends == 3)
      for (int i = 0; i < 6; i++)
	in3 >> dummy;

    if (nfound != 0) {
      assert(nfound == factors.size());
#ifdef BLEEN
      if (nfound != factors.size() && !graph.parallelFaces()) {
	// cout << "missing factor " << (iRep+1) << endl;
	cout << "missing factor ";
	graph.print(cout);
	cout << endl;
      }
#endif
      if (nfound > 1 && nfound == nContactFactors) {
	cout << "interesting ";
	graph.print(cout);
	cout << endl;
      }
    }
    else if (factors.size() == 1) {
      // cout << endl;
      if (graph.n < 5) {
	Alists lists = graph.getLists();
	Alists perm = lists;
	for (int i = 0; i < 16; i++)
	  perm.a[i] = 0;
        basic = graph;
	basic.normalize();
	assert(basic == graph);
	basic_coeffs = factors[0];

#ifdef BLEEN
	Face4 face4 = graph.getFace4();
	AnglePoly *poly = c4d->getPoly(face4);
	extern vector<int> coeffs26;
	vector<PTR<Angle>> roots = poly->getAnglesSorted();
	assert(roots.size() == 0);
	delete poly;
	for (int i = 0; i < coeffs26.size(); i++)
	  basic_coeffs[i] = coeffs26[i];

	if (coeffs26.size() == factors[0].size()) {
	  for (int i = 0; i < coeffs26.size(); i++)
	    assert(coeffs26[i] != 0 && factors[0][i] == coeffs26[i]);
	  static int biggestD = 0;
	  if (coeffs26.size() > biggestD) {
	    biggestD = coeffs26.size();
	    cout << "biggestD " << biggestD << endl;
	  }
	}
	else {
	  static int biggestN = 0;
	  static int biggestD = 0;
	  if (graph.n > biggestN || coeffs26.size() > biggestD) {
	    cout << "graph.n " << graph.n << " "
		 << "coeffs26.size() " << coeffs26.size() << " "
		 << "factors[0].size() " << factors[0].size() << endl;
	    biggestN = graph.n;
	    biggestD = coeffs26.size();
	  }

	  assert(coeffs26.size() == factors[0].size()+2);
	  vector<int> coeffs26d = divideBy1plusS2(coeffs26);
	  for (int i = 0; i < coeffs26d.size(); i++)
	    assert(coeffs26d[i] != 0 && factors[0][i] == coeffs26d[i]);
	}
#endif

	graph_keys.clear();
	really_basic = true;
	createPerms(lists, perm, 0, 0);

        auto its = vec2graphs.find(basic_coeffs);
        if (its == vec2graphs.end()) {
          out << " ";
          graph.print(out);
        }
        else {
          out << " -99 " << its->second.size();
          for (int i = 0; i < its->second.size(); i++) {
            out << " ";
            its->second[i].print(out);
          }
        }

#ifdef CHECK_FOR_MATCHES
	if (really_basic && graph.n == 4) {
	  auto it = vec2graph.find(basic_coeffs);
	  assert(it == vec2graph.end());
	  vec2graph[basic_coeffs] = basic;
	}
#endif
      }
    }

    out << endl;
  }

  out.close();
}  

#ifdef BLEEN
void C4D::findIdentities () {
  c4d = this;
  alts[1111] = 24;
  alts[1112] = 6;
  alts[1113] = 6;
  alts[1122] = 4;
  alts[1123] = 2;
  alts[1133] = 4;
  alts[1222] = 6;
  alts[1223] = 2;
  alts[1233] = 2;
  alts[1333] = 6;
  alts[2222] = 24;
  alts[2223] = 6;
  alts[2233] = 4;
  alts[2333] = 6;
  alts[3333] = 24;
  int nalts = 0;

  vector<Face4> face4s;
  map<Face4, int, Face4::Compare> face4index;
  vector< vector<Face4> > factorss;
  ifstream in;
  in.open("../identity/min-class.txt");
  ifstream in2;
  in2.open("../identity/eqsub-classes.txt");
  int n;
  while (in >> n) {
    int nF[4] = { 0, 0, 0, 0 };
    int nM[4] = { 0, 0, 0, 0 };
    int indsF[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
    int indsM[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
    
    for (int i = 0; i < n; i++) {
      in >> nF[i];
      for (int j = 0; j < nF[i]; j++)
	in >> indsF[i][j];
      in >> nM[i];
      for (int j = 0; j < nM[i]; j++)
	in >> indsM[i][j];
    }
    
    Face4 four(nF, indsF, nM, indsM);
    face4index[four] = face4s.size();
    // int index = face4index[four];
    // assert(index == face4s.size());
    face4s.push_back(four);

    int m;
    in2 >> m;
    vector<Face4> fs;
    
    for (int i = 0; i < m; i++) {
      int n2;
      in2 >> n2;

      if (n2 <= 0)
	continue;

      int nF[4] = { 0, 0, 0, 0 };
      int nM[4] = { 0, 0, 0, 0 };
      int indsF[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
      int indsM[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
      
      for (int i = 0; i < n2; i++) {
	in2 >> nF[i];
	for (int j = 0; j < nF[i]; j++)
	  in2 >> indsF[i][j];
	in2 >> nM[i];
	for (int j = 0; j < nM[i]; j++)
	  in2 >> indsM[i][j];
      }
      
      Face4 four(nF, indsF, nM, indsM);
      fs.push_back(four);
    }

    factorss.push_back(fs);
  }

  map<Face4, Perm4, Face4::Compare> face4perm;
  ifstream inPerms;
  inPerms.open("../identity/perms.txt");
  while (inPerms >> n) {
        int nF[4] = { 0, 0, 0, 0 };
    int nM[4] = { 0, 0, 0, 0 };
    int indsF[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
    int indsM[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, };
    
    for (int i = 0; i < n; i++) {
      inPerms >> nF[i];
      for (int j = 0; j < nF[i]; j++)
	inPerms >> indsF[i][j];
      inPerms >> nM[i];
      for (int j = 0; j < nM[i]; j++)
	inPerms >> indsM[i][j];
    }
    
    Perm4 perm;
    for (int i = 0; i < 4; i++)
      inPerms >> perm.perm[i];

    Face4 four(nF, indsF, nM, indsM);
    face4perm[four] = perm;
  }

#define NPOOL 32
  // int nFM[] = { 32, 32 };
  int nFM[] = { NPOOL, NPOOL };
  for (int fm = 0; fm <= 1; fm++) {
    for (int i = 0; i < nFM[fm]; i++)
      if (true)
	new VertXYZ(this,
		    randomNumber(-1, 1),
		    randomNumber(-1, 1),
		    randomNumber(-1, 1),
		    fm);
      else
	new VertXYZ(this,
		    randomInt(0, 1024),
		    randomInt(0, 1024),
		    randomInt(0, 1024),
		    fm, true);

    Int2 i2;
    for (i2.i[0] = 0; i2.i[0] < nFM[fm]; i2.i[0]++)
      for (i2.i[1] = i2.i[0]+1; i2.i[1] < nFM[fm]; i2.i[1]++)
	i2edge[fm][i2] = new EdgeAB(verts[fm][i2.i[0]], verts[fm][i2.i[1]]);

    Int3 i3;
    for (i3.i[0] = 0; i3.i[0] < nFM[fm]; i3.i[0]++)
      for (i3.i[1] = i3.i[0]+1; i3.i[1] < nFM[fm]; i3.i[1]++)
	for (i3.i[2] = i3.i[1]+1; i3.i[2] < nFM[fm]; i3.i[2]++)
	  i3face[fm][i3] = new FaceABC(verts[fm][i3.i[0]], verts[fm][i3.i[1]], verts[fm][i3.i[2]]);
  }

  initSums(true);
  enable();

  set<Face4, Face4::Compare> solved;

  int nResolved = 0;
  int nPermsFound = 0;

  for (int iTest = 0; iTest < 1000000; iTest++) {
    // int iClass = 4564;
    int iClass = randomInt(0, face4s.size());
    int r[2][12];
    permutation(12, r[0]);
    permutation(12, r[1]);

    Face4 face4rep = face4s[iClass];
#ifdef ONLINE_LOOKUP
    int face4repI = face4index[face4rep];
    assert(face4repI == iClass);

    if (false) {
    Face4 face4ins = face4rep;
    face4ins.replace(r);
    Graph graph(face4ins);

    int iTo012[2][NPOOL];
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < NPOOL; i++)
	iTo012[fm][i] = -1;
    int iFrom012[2][12];
    int n012[2] = { 0, 0 };
    
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < 4; i++)
	for (int j = 0; j < face4ins.ns[fm][i]; j++) {
	  if (iTo012[fm][face4ins.is[fm][i][j]] == -1) {
	    iTo012[fm][face4ins.is[fm][i][j]] = n012[fm];
	    iFrom012[fm][n012[fm]] = face4ins.is[fm][i][j];
	    n012[fm]++;
	  }
	  if (fm == 0)
	    graph.f[i][j] = iTo012[fm][face4ins.is[fm][i][j]];
	  else
	    graph.m[i][j] = iTo012[fm][face4ins.is[fm][i][j]];
	}

    int perm[2][12];
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < n012[fm]; i++)
	perm[fm][i] = i;

    Alists lists = graph.getLists(perm);
    graph.apply(lists);

    Face4 face4lup = graph.getFace4();
    map<Face4, Perm4, Face4::Compare>::iterator it = face4perm.find(face4lup);
    if (it != face4perm.end())
      nPermsFound++;
    }

    Face4 face4ins = face4rep;
    face4ins.replace2(r);
    Graph graph(face4ins);

    int iTo012[2][NPOOL];
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < NPOOL; i++)
	iTo012[fm][i] = -1;
    int iFrom012[2][12];
    int n012[2] = { 0, 0 };
    
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < 4; i++)
	for (int j = 0; j < face4ins.ns[fm][i]; j++) {
	  if (iTo012[fm][face4ins.is[fm][i][j]] == -1) {
	    iTo012[fm][face4ins.is[fm][i][j]] = n012[fm];
	    iFrom012[fm][n012[fm]] = face4ins.is[fm][i][j];
	    n012[fm]++;
	  }
	  if (fm == 0)
	    graph.f[i][j] = iTo012[fm][face4ins.is[fm][i][j]];
	  else
	    graph.m[i][j] = iTo012[fm][face4ins.is[fm][i][j]];
	}

    int perm[2][12];
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < n012[fm]; i++)
	perm[fm][i] = i;

    Alists lists = graph.getLists(perm);
    graph.apply(lists);

    Face4 face4lup = graph.getFace4();

    assert(!(face4lup < face4rep) && !(face4rep < face4lup));
#endif

#ifdef BLEEN
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < n012[fm]; i++)
	assert(r[fm][i] == iFrom012[fm][perm[fm][i]]);
#endif

    vector<Face4> factors = factorss[iClass];

#ifndef TABLE_IDENTITY
    if (factors.size() > 1) {
      factors.clear();
      factors.push_back(face4rep);
    }
#endif    

    for (int j = 0; j < factors.size(); j++) {
      Face4 factor = factors[j];
      factor.replace(r);
      if (solved.find(factor) == solved.end()) {
	if (solved.find(factor) != solved.end() && factor.ns[0][2] != 0)
	  nResolved++;
	solved.insert(factor);
	AnglePoly *poly = getPoly(factor);
	vector<PTR<Angle>> roots = poly->getAnglesSorted();
	// delete them?

        Contact contact = getContact(poly);
        AnglePoly *poly2 = getPoly(contact);
        Contact contact2 = getContact(poly2);

	assert(contact == contact2);
	if (poly->getType() == AnglePoly::EFF) {
	  assert(PolyEFF::Equal()(dynamic_cast<PolyEFF*>(poly),
				  dynamic_cast<PolyEFF*>(poly2)));
	}
	else if (poly->getType() == AnglePoly::FFFF) {
	  assert(PolyFFFF::Equal()(dynamic_cast<PolyFFFF*>(poly),
				   dynamic_cast<PolyFFFF*>(poly2)));
	}
	else
	  assert(poly->identical(poly2));

	for (int h = 0; h < roots.size(); h++) {
	  Angle *root = roots[h];
	  Event *event = new Event(this, root);
	  eventList.insert(event);
	}
      }
      else
	nResolved++;
    }
  }

  cout << "nResolved " << nResolved << endl;
  cout << "nPermsFound " << nPermsFound << endl;
  cout << "nRootSeps " << nRootSeps << endl;
}  
#endif

void C4D::findIdentities (int nFixed, int nMoving) {
  int nF[4] = { 2, 2, 2, 2 };
  int nM[4] = { 2, 2, 2, 2 };
  int indsF[4][3] = {{0, 1, -1}, {0, 2, -1}, {1, 2, -1}, {3, 4, -1}};
  int indsM[4][3] = {{0, 1, -1}, {0, 2, -1}, {1, 2, -1}, {3, 4, -1}};
#ifdef BLEEN
  int nF[4] = {1, 2, 2, 2};
  int nM[4] = {3, 2, 2, 2};
  int indsF[4][3] = {{0, -1, -1}, {0, 1, -1}, {2, 3, -1}, {2, 4, -1}};
  int indsM[4][3] = {{0, 1, 2}, {3, 4, -1}, {0, 1, -1}, {0, 3, -1}};
#endif

  Face4 four(nF, indsF, nM, indsM);
  Face4 fmin(four);
  for (int i = 0; i < 10; i++)
    fmin.minimize();

  int perm5[] = { 0, 1, 2, 3, 4 };

  set<Face4, Face4::Compare> fours, foursu;

  for (int p = 0; p < 000000; p++) {
    if (p % 10000 == 0)
      cout << "p " << p << " fours " << fours.size() << endl;
    for (int fm = 0; fm <= 1; fm++) {
      permute(5, perm5);
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < four.ns[fm][i]; j++)
          four.is[fm][i][j] = perm5[four.is[fm][i][j]];
	permute(four.ns[fm][i], four.is[fm][i]);
      }
    }
    int ij[2];
    randomInds(2, 4, ij);
    four.swapFaces(ij[0], ij[1]);
    // if (random() % 2 == 0)
    // four = Face4(four, true);

    foursu.insert(four);

    for (int fm = 10; fm <= 1; fm++)
      for (int i = 0; i < 4; i++)
	for (int j0 = 0; j0 < four.ns[fm][i]; j0++)
	  for (int j1 = j0+1; j1 < four.ns[fm][i]; j1++)
	    if (four.is[fm][i][j0] > four.is[fm][i][j1]) {
	      int temp = four.is[fm][i][j0];
	      four.is[fm][i][j0] = four.is[fm][i][j1];
	      four.is[fm][i][j1] = temp;
	    }  

    fours.insert(four);

    Face4 face4(four);
    face4.minimize();
    assert(face4 == fmin);
  }

  int itn = 0;
  set<Face4, Face4::Compare>::iterator itu = foursu.begin();
  for (set<Face4, Face4::Compare>::iterator it = fours.begin();
       it != fours.end(); ++it, ++itu) {
    bool sinu = (foursu.find(*it) != foursu.end());
    bool uins = (fours.find(*itu) != fours.end());
    cout << sinu << " " << uins << endl;
    if (itn++ % 1971 == 0) {
      for (int j = 0; j < (*it).ns[0][0]; j++)
	cout << (*it).is[0][0][j] << " ";
      cout << endl;
    }
  }

  for (long iTrial = 0; iTrial < 0000000; iTrial++) {
    if (iTrial % 100000 == 0)
      cout << "iTrial " << iTrial << " fours " << fours.size() << endl;
    SumFace *faces[4] = { 0, 0, 0, 0 };
    int iFace;
    int nF[4] = { 0, 0, 0, 0};
    int indsF[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}};
    int nM[4] = { 0, 0, 0, 0};
    int indsM[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}};
    for (iFace = 0; iFace < 4; iFace++) {
      int inds[4];
      do {
        randomInds(4, nFixed + nMoving, inds);
      } while (inds[0] >= nFixed || inds[3] < nFixed);
      // !(inds[1] < nFixed && inds[2] >= nFixed));
      while (inds[nF[iFace]] < nFixed) {
        indsF[iFace][nF[iFace]] = inds[nF[iFace]];
        nF[iFace]++;
      }
      while (nF[iFace] + nM[iFace] < 4) {
        indsM[iFace][nM[iFace]] = inds[nF[iFace] + nM[iFace]] - nFixed;
        nM[iFace]++;
      }
      // assert(nF[iFace] == 2 && nM[iFace] == 2);
    }

    Face4 face4(nF, indsF, nM, indsM);
    Face4 fmin2(face4);
    fmin2.minimize();
    if (fmin2 == fmin)
      fours.insert(face4);
  }

  cout << "num fours " << fours.size() << endl;

  if (false)
    return;

  int nFM[] = { nFixed, nMoving };
  for (int fm = 0; fm <= 1; fm++) {
    for (int i = 0; i < nFM[fm]; i++)
      new VertXYZ(this,
		  randomNumber(-1, 1), randomNumber(-1, 1), randomNumber(-1, 1),
		  fm);

    Int2 i2;
    for (i2.i[0] = 0; i2.i[0] < nFM[fm]; i2.i[0]++)
      for (i2.i[1] = i2.i[0]+1; i2.i[1] < nFM[fm]; i2.i[1]++)
	i2edge[fm][i2] = new EdgeAB(verts[fm][i2.i[0]], verts[fm][i2.i[1]]);

    Int3 i3;
    for (i3.i[0] = 0; i3.i[0] < nFM[fm]; i3.i[0]++)
      for (i3.i[1] = i3.i[0]+1; i3.i[1] < nFM[fm]; i3.i[1]++)
	for (i3.i[2] = i3.i[1]+1; i3.i[2] < nFM[fm]; i3.i[2]++)
	  i3face[fm][i3] = new FaceABC(verts[fm][i3.i[0]], verts[fm][i3.i[1]], verts[fm][i3.i[2]]);
  }

  enable();
  initSums(true);

  map<Face4, int, Face4::Compare> face4num;
  int nface4num = 0;
  int nfour = 0;

  for (int iTrial = 0; face4num.size() < 674; iTrial++) {
    if (iTrial % 100000 == 0)
      // cout << "iTrial " << iTrial << " nfour " << nfour << endl;
      cout << "iTrial " << iTrial << " face4num " << face4num.size() << endl;
    SumFace *faces[4] = { 0, 0, 0, 0 };
    int iFace;
    int nF[4] = { 0, 0, 0, 0};
    int indsF[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}};
    int nM[4] = { 0, 0, 0, 0};
    int indsM[4][3] = {{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}};
    for (iFace = 0; iFace < 4; iFace++) {
      int inds[4];
      do {
        randomInds(4, nFixed + nMoving, inds);
      } while (inds[0] >= nFixed || inds[3] < nFixed);
      while (inds[nF[iFace]] < nFixed) {
        indsF[iFace][nF[iFace]] = inds[nF[iFace]];
        nF[iFace]++;
      }
      while (nF[iFace] + nM[iFace] < 4) {
        indsM[iFace][nM[iFace]] = inds[nF[iFace] + nM[iFace]] - nFixed;
        nM[iFace]++;
      }
#ifdef BLEEN
      nF[iFace] = 1 + random() % 3;
      randomInds(nF[iFace], nFixed, indsF[iFace]);
      nM[iFace] = 4 - nF[iFace];
      randomInds(nM[iFace], nMoving, indsM[iFace]);
#endif
      switch (nF[iFace]) {
      case 1: {
        VertXYZ *a = verts[0][indsF[iFace][0]];
        FaceABC *bcd = i3face[1][indsM[iFace]];
        faces[iFace] = SumVF::make(a, bcd);
        break;
      }
      case 2: {
        EdgeAB *ab = i2edge[0][indsF[iFace]];
        EdgeAB *cd = i2edge[1][indsM[iFace]];
        faces[iFace] = SumEE::make(ab, cd);
        break;
      }
      case 3: {
        FaceABC *abc = i3face[0][indsF[iFace]];
        VertXYZ *d = verts[1][indsM[iFace][0]];
        faces[iFace] = SumFV::make(abc, d);
        break;
      }
      default:
        assert(0);
        break;
      }

      int jFace;
      for (jFace = 0; jFace < iFace; jFace++) {
        if (faces[iFace]->neighborOf2(faces[jFace]))
          break;
        PolyEF *polyEF = faces[iFace]->related(faces[jFace]);
        if (polyEF != 0) {
          delete polyEF;
          break;
        }
      }
      if (jFace < iFace)
        break;
    }
    if (iFace < 4)
      continue;

    if (numUnique(nFixed, nF, indsF) != nFixed ||
        numUnique(nMoving, nM, indsM) != nMoving)
      continue;

    // Check if all four faces share a (sum) vertex.
    const vector<Edge*> &edges0 = faces[0]->getOuterLoop();
    int i;
    for (i = 0; i < edges0.size(); i++) {
      Vert *vert = edges0[i]->tail;
      int j;
      for (j = 1; j < 4; j++)
        if (!faces[j]->contains(vert))
          break;
      if (j == 4)
        break;
    }
    if (i != edges0.size())
      continue;

    // Check if three faces share a (sum) vertex.
    for (i = 0; i < edges0.size(); i++) {
      Vert *vert = edges0[i]->tail;
      int n = 1;
      int j;
      for (j = 1; j < 4; j++)
        if (faces[j]->contains(vert))
          n++;
      if (n >= 3)
        break;
    }
    if (i != edges0.size())
      continue;

    const vector<Edge*> &edges1 = faces[1]->getOuterLoop();
    for (i = 0; i < edges1.size(); i++) {
      Vert *vert = edges1[i]->tail;
      int n = 1;
      int j;
      for (j = 2; j < 4; j++)
        if (faces[j]->contains(vert))
          n++;
      if (n >= 3)
        break;
    }
    if (i != edges1.size())
      continue;

    // Check if three fixed edges or triangles share an edge.
    for (i = 0; i < 4; i++) {
      int i0 = i;
      int i1 = (i+1)%4;
      int i2 = (i+2)%4;
      int io[3], o = 0;
      o = intersect(indsF[i0], nF[i0], indsF[i1], nF[i1], io);
      o = intersect(io, o, indsF[i2], nF[i2], io);
      if (o >= 2)
        break;
    }
    if (i != 4)
      continue;

    // Check if three moving edges or triangles share an edge.
    for (i = 0; i < 4; i++) {
      int i0 = i;
      int i1 = (i+1)%4;
      int i2 = (i+2)%4;
      int io[3], o = 0;
      o = intersect(indsM[i0], nM[i0], indsM[i1], nM[i1], io);
      o = intersect(io, o, indsM[i2], nM[i2], io);
      if (o >= 2)
        break;
    }
    if (i != 4)
      continue;

    Face4 face4(nF, indsF, nM, indsM);
    face4.minimize();

    if (false) {
      Face4 face4(nF, indsF, nM, indsM);
      Face4 fmin2(face4);
      fmin2.minimize();
      if (fmin2 == fmin)
	nfour++;
    }

    // fixed are 0, 0, 0 or contains 01, contains 01
    if (face4.ns[0][0] == 1 && face4.is[0][0][0] == 0 &&
	face4.ns[0][1] == 1 && face4.is[0][1][0] == 0 &&
	((face4.ns[0][2] == 1 && face4.is[0][2][0] == 0) ||
	 (face4.is[0][2][0] == 0 && face4.is[0][2][1] == 1)) &&
	face4.is[0][3][0] == 0 && face4.is[0][3][1] == 1)
      continue;

    // 0+0xx, 0+0xx, 1+0xx, 1+0xx
    if (face4.ns[0][0] == 1 && face4.is[0][0][0] == 0 &&
	face4.ns[0][1] == 1 && face4.is[0][1][0] == 0 &&
	face4.ns[0][2] == 1 && face4.is[0][2][0] == 1 &&
	face4.ns[0][3] == 1 && face4.is[0][3][0] == 1 &&
	face4.is[1][0][0] == 0 &&
	face4.is[1][1][0] == 0 &&
	face4.is[1][2][0] == 0 &&
	face4.is[1][3][0] == 0)
      continue;

    int nsI[2][4] = {{1, 1, 1, 2}, {3, 3, 3, 2}};
    int isI[2][4][3] = {{{0, -1, -1}, {1, -1, -1}, {2, -1, -1}, {0, 1, -1}}, {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {0, 4, -1}}};
    Face4 face4I(nsI[0], isI[0], nsI[1], isI[1]);
    if (face4 == face4I)
      continue;

    int nsA[2][4] = {{1, 1, 2, 2}, {3, 3, 2, 2}};
    int isA[2][4][3] = {{{0, -1, -1}, {1, -1, -1}, {0, 2, -1}, {1, 2, -1}}, {{0, 1, 2}, {0, 1, 3}, {0, 3, -1}, {0, 2, -1}}};
    Face4 face4A(nsA[0], isA[0], nsA[1], isA[1]);
    if (face4 == face4A)
      continue;

    int nsL[2][4] = {{1, 1, 2, 2}, {3, 3, 2, 2}};
    int isL[2][4][3] = {{{0, -1, -1}, {1, -1, -1}, {0, 2, -1}, {1, 2, -1}}, {{0, 1, 2}, {0, 3, 4}, {0, 3, -1}, {0, 1, -1}}};
    Face4 face4L(nsL[0], isL[0], nsL[1], isL[1]);
    if (face4 == face4L)
      continue;

    int nsD[2][4] = {{1, 1, 3, 3}, {3, 3, 1, 1}};
    int isD[2][4][3] = {{{0, -1, -1}, {1, -1, -1}, {0, 1, 2}, {0, 1, 3}}, {{0, 1, 2}, {0, 1, 3}, {0, -1, -1}, {1, -1, -1}}};
    Face4 face4D(nsD[0], isD[0], nsD[1], isD[1]);
    if (face4 == face4D)
      continue;

    int nsC[2][4] = {{1, 2, 2, 3}, {3, 2, 2, 1}};
    int isC[2][4][3] = {{{0, -1, -1}, {0, 1, -1}, {0, 2, -1}, {0, 1, 2}}, {{0, 1, 2}, {0, 3, -1}, {1, 3, -1}, {2, -1, -1}}};
    Face4 face4C(nsC[0], isC[0], nsC[1], isC[1]);
    if (face4 == face4C)
      continue;

    int nsB[2][4] = {{2, 2, 2, 2}, {2, 2, 2, 2}};
    int isB[2][4][3] = {{{0, 1, -1}, {0, 1, -1}, {0, 2, -1}, {0, 2, -1}}, {{0, 1, -1}, {2, 3, -1}, {0, 2, -1}, {1, 3, -1}}};
    Face4 face4B(nsB[0], isB[0], nsB[1], isB[1]);
    if (face4 == face4B)
      continue;
    
    face4num[face4]++;
    nface4num++;

    /*
      PolyFFFF temp(faces[0], faces[1], faces[2], faces[3]);
      if (polyFFFFs.find(&temp) != polyFFFFs.end())
      continue;
      PolyFFFF *poly = new PolyFFFF(faces[0], faces[1], faces[2], faces[3]);
      polyFFFFs.insert(poly);
      poly->incRef();

      vector<PTR<Angle>> roots = poly->getAnglesSorted();
      for (int h = 0; h < roots.size(); h++) {
      Angle *root = roots[h];
      Event *event = new Event(this, root);
      eventList.insert(event);
      }
    */
  }

  Face4 minface4 = (*face4num.begin()).first;
  int minnum = (*face4num.begin()).second;
  Face4 maxface4 = (*face4num.begin()).first;
  int maxnum = (*face4num.begin()).second;
  for (map<Face4, int, Face4::Compare>::iterator it = face4num.begin();
       it != face4num.end(); ++it) {
    cout << (*it).second << endl;
    if ((*it).second < minnum) {
      minface4 = (*it).first;
      minnum = (*it).second;
    }
    if ((*it).second > maxnum) {
      maxface4 = (*it).first;
      maxnum = (*it).second;
    }
  }
  cout << face4num.size() << " " << minnum << endl;
  cout << face4num.size() << " " << maxnum << endl;

  vector<Face4> face4s(face4num.size());
  int nface4 = 0;
  for (map<Face4, int, Face4::Compare>::iterator it = face4num.begin();
       it != face4num.end(); ++it)
    face4s[nface4++] = (*it).first;
    
  nface4num = 0;
  face4num.clear();

  int permF[nFixed];
  for (int i = 0; i < nFixed; i++)
    permF[i] = i;
  int permM[nMoving];
  for (int i = 0; i < nMoving; i++)
    permM[i] = i;

  for (int iTrial = 0; nface4num < 100000000; iTrial++) {
    if (iTrial % 100000 == 0)
      cout << "iTrial " << iTrial << " nface4num " << nface4num << endl;

    Face4 face4r = face4s[random() % face4s.size()];
    permute(nFixed, permF);
    permute(nMoving, permM);

    permute(nFixed, permF);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < four.ns[0][i]; j++)
        four.is[0][i][j] = permF[four.is[0][i][j]];
      
    permute(nMoving, permM);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < four.ns[1][i]; j++)
        four.is[1][i][j] = permM[four.is[1][i][j]];
      
    for (int fm = 0; fm <= 1; fm++)
      for (int i = 0; i < 4; i++)
        for (int j0 = 0; j0 < four.ns[fm][i]; j0++)
          for (int j1 = j0+1; j1 < four.ns[fm][i]; j1++)
            if (four.is[fm][i][j0] > four.is[fm][i][j1]) {
              int temp = four.is[fm][i][j0];
              four.is[fm][i][j0] = four.is[fm][i][j1];
              four.is[fm][i][j1] = temp;
            }

    int *nF = face4r.ns[0];
    int *nM = face4r.ns[1];
    
    typedef int int43[4][3];
    int43 &indsF = face4r.is[0];
    int43 &indsM = face4r.is[1];
    SumFace *faces[4] = { 0, 0, 0, 0 };

    int iFace;
    for (iFace = 0; iFace < 4; iFace++) {
      switch (nF[iFace]) {
      case 1: {
        VertXYZ *a = verts[0][indsF[iFace][0]];
        FaceABC *bcd = i3face[1][indsM[iFace]];
        faces[iFace] = SumVF::make(a, bcd);
        break;
      }
      case 2: {
        EdgeAB *ab = i2edge[0][indsF[iFace]];
        EdgeAB *cd = i2edge[1][indsM[iFace]];
        faces[iFace] = SumEE::make(ab, cd);
        break;
      }
      case 3: {
        FaceABC *abc = i3face[0][indsF[iFace]];
        VertXYZ *d = verts[1][indsM[iFace][0]];
        faces[iFace] = SumFV::make(abc, d);
        break;
      }
      default:
        assert(0);
        break;
      }

      int jFace;
      for (jFace = 0; jFace < iFace; jFace++) {
        if (faces[iFace]->neighborOf2(faces[jFace]))
          break;
        PolyEF *polyEF = faces[iFace]->related(faces[jFace]);
        if (polyEF != 0) {
          delete polyEF;
          break;
        }
      }
      if (jFace < iFace)
        break;
    }
    if (iFace < 4)
      continue;

    if (numUnique(nFixed, nF, indsF) != nFixed ||
        numUnique(nMoving, nM, indsM) != nMoving)
      continue;

    // Check if all four faces share a (sum) vertex.
    const vector<Edge*> &edges0 = faces[0]->getOuterLoop();
    int i;
    for (i = 0; i < edges0.size(); i++) {
      Vert *vert = edges0[i]->tail;
      int j;
      for (j = 1; j < 4; j++)
        if (!faces[j]->contains(vert))
          break;
      if (j == 4)
        break;
    }
    if (i != edges0.size())
      continue;

    // Check if three faces share a (sum) vertex.
    for (i = 0; i < edges0.size(); i++) {
      Vert *vert = edges0[i]->tail;
      int n = 1;
      int j;
      for (j = 1; j < 4; j++)
        if (faces[j]->contains(vert))
          n++;
      if (n >= 3)
        break;
    }
    if (i != edges0.size())
      continue;

    const vector<Edge*> &edges1 = faces[1]->getOuterLoop();
    for (i = 0; i < edges1.size(); i++) {
      Vert *vert = edges1[i]->tail;
      int n = 1;
      int j;
      for (j = 2; j < 4; j++)
        if (faces[j]->contains(vert))
          n++;
      if (n >= 3)
        break;
    }
    if (i != edges1.size())
      continue;

    // Check if three fixed edges or triangles share an edge.
    for (i = 0; i < 4; i++) {
      int i0 = i;
      int i1 = (i+1)%4;
      int i2 = (i+2)%4;
      int io[3], o = 0;
      o = intersect(indsF[i0], nF[i0], indsF[i1], nF[i1], io);
      o = intersect(io, o, indsF[i2], nF[i2], io);
      if (o >= 2)
        break;
    }
    if (i != 4)
      continue;

    // Check if three moving edges or triangles share an edge.
    for (i = 0; i < 4; i++) {
      int i0 = i;
      int i1 = (i+1)%4;
      int i2 = (i+2)%4;
      int io[3], o = 0;
      o = intersect(indsM[i0], nM[i0], indsM[i1], nM[i1], io);
      o = intersect(io, o, indsM[i2], nM[i2], io);
      if (o >= 2)
        break;
    }
    if (i != 4)
      continue;

    Face4 face4(nF, indsF, nM, indsM);
    face4.minimize();

    if (false) {
      Face4 face4(nF, indsF, nM, indsM);
      Face4 fmin2(face4);
      fmin2.minimize();
      if (fmin2 == fmin)
	nfour++;
    }

    // fixed are 0, 0, 0 or contains 01, contains 01
    if (face4.ns[0][0] == 1 && face4.is[0][0][0] == 0 &&
	face4.ns[0][1] == 1 && face4.is[0][1][0] == 0 &&
	((face4.ns[0][2] == 1 && face4.is[0][2][0] == 0) ||
	 (face4.is[0][2][0] == 0 && face4.is[0][2][1] == 1)) &&
	face4.is[0][3][0] == 0 && face4.is[0][3][1] == 1)
      continue;

    // 0+0xx, 0+0xx, 1+0xx, 1+0xx
    if (face4.ns[0][0] == 1 && face4.is[0][0][0] == 0 &&
	face4.ns[0][1] == 1 && face4.is[0][1][0] == 0 &&
	face4.ns[0][2] == 1 && face4.is[0][2][0] == 1 &&
	face4.ns[0][3] == 1 && face4.is[0][3][0] == 1 &&
	face4.is[1][0][0] == 0 &&
	face4.is[1][1][0] == 0 &&
	face4.is[1][2][0] == 0 &&
	face4.is[1][3][0] == 0)
      continue;

    int nsI[2][4] = {{1, 1, 1, 2}, {3, 3, 3, 2}};
    int isI[2][4][3] = {{{0, -1, -1}, {1, -1, -1}, {2, -1, -1}, {0, 1, -1}}, {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {0, 4, -1}}};
    Face4 face4I(nsI[0], isI[0], nsI[1], isI[1]);
    if (face4 == face4I)
      continue;

    int nsA[2][4] = {{1, 1, 2, 2}, {3, 3, 2, 2}};
    int isA[2][4][3] = {{{0, -1, -1}, {1, -1, -1}, {0, 2, -1}, {1, 2, -1}}, {{0, 1, 2}, {0, 1, 3}, {0, 3, -1}, {0, 2, -1}}};
    Face4 face4A(nsA[0], isA[0], nsA[1], isA[1]);
    if (face4 == face4A)
      continue;

    int nsL[2][4] = {{1, 1, 2, 2}, {3, 3, 2, 2}};
    int isL[2][4][3] = {{{0, -1, -1}, {1, -1, -1}, {0, 2, -1}, {1, 2, -1}}, {{0, 1, 2}, {0, 3, 4}, {0, 3, -1}, {0, 1, -1}}};
    Face4 face4L(nsL[0], isL[0], nsL[1], isL[1]);
    if (face4 == face4L)
      continue;

    int nsD[2][4] = {{1, 1, 3, 3}, {3, 3, 1, 1}};
    int isD[2][4][3] = {{{0, -1, -1}, {1, -1, -1}, {0, 1, 2}, {0, 1, 3}}, {{0, 1, 2}, {0, 1, 3}, {0, -1, -1}, {1, -1, -1}}};
    Face4 face4D(nsD[0], isD[0], nsD[1], isD[1]);
    if (face4 == face4D)
      continue;

    int nsC[2][4] = {{1, 2, 2, 3}, {3, 2, 2, 1}};
    int isC[2][4][3] = {{{0, -1, -1}, {0, 1, -1}, {0, 2, -1}, {0, 1, 2}}, {{0, 1, 2}, {0, 3, -1}, {1, 3, -1}, {2, -1, -1}}};
    Face4 face4C(nsC[0], isC[0], nsC[1], isC[1]);
    if (face4 == face4C)
      continue;

    int nsB[2][4] = {{2, 2, 2, 2}, {2, 2, 2, 2}};
    int isB[2][4][3] = {{{0, 1, -1}, {0, 1, -1}, {0, 2, -1}, {0, 2, -1}}, {{0, 1, -1}, {2, 3, -1}, {0, 2, -1}, {1, 3, -1}}};
    Face4 face4B(nsB[0], isB[0], nsB[1], isB[1]);
    if (face4 == face4B)
      continue;
    
    face4num[face4]++;
    nface4num++;

    PolyFFFF temp(faces[0], faces[1], faces[2], faces[3]);
    if (polyFFFFs.find(&temp) != polyFFFFs.end())
      continue;
    PolyFFFF *poly = new PolyFFFF(faces[0], faces[1], faces[2], faces[3]);
    polyFFFFs.insert(poly);
    //poly->incRef();

    vector<PTR<Angle>> roots = poly->getAnglesSorted();
    for (int h = 0; h < roots.size(); h++) {
      Angle *root = roots[h];
      Event *event = new Event(this, root);
      eventList.insert(event);
    }
  }

  minface4 = (*face4num.begin()).first;
  minnum = (*face4num.begin()).second;
  maxface4 = (*face4num.begin()).first;
  maxnum = (*face4num.begin()).second;
  for (map<Face4, int, Face4::Compare>::iterator it = face4num.begin();
       it != face4num.end(); ++it) {
    cout << (*it).second << endl;
    if ((*it).second < minnum) {
      minface4 = (*it).first;
      minnum = (*it).second;
    }
    if ((*it).second > maxnum) {
      maxface4 = (*it).first;
      maxnum = (*it).second;
    }
  }
  cout << face4num.size() << " " << minnum << endl;
  cout << face4num.size() << " " << maxnum << endl;
}

#ifdef C4D_MAIN
int main (int argc, char *argv[]) {
  if (argc > 1) {
    int seed = atoi(argv[1]);
    cout << "seed " << seed << endl;
    srandom(seed);
  }

  C4D c4d;

#ifdef BLEEN
  // c4d.findIdentities(6, 6);
  c4d.findIdentities();
  if (true)
    return 0;
#endif
  
  Parameter::delta = 1e-4;

  const char *input = "frame.obj";

  cout << "reading fixed" << endl;
  c4d.readObjFile(input, false);
  
  if (false) {
    C4D::VertXYZ *a = c4d.verts[0][4];
    C4D::VertXYZ *b = c4d.verts[1][7];
    C4D::EdgeAB *ab = c4d.edges[1][6];
    C4D::FaceABC *abc = c4d.faces[0][9];
    C4D::SumVV *svv = C4D::SumVV::get(a, b);
    C4D::SumVE *sve = C4D::SumVE::get(a, ab);
    C4D::SumFV *sfv = C4D::SumFV::get(abc, b);
  }

  cout << "reading moving" << endl;
  c4d.readObjFile(input, true);
  cout << "init sums" << endl;
  c4d.initSums();

  if (true) {
    enable();
    Contact::readFactors();

    int iFace = atoi(argv[2]);
    cout << "iFace " << iFace << endl;
    c4d.sweep2(iFace);
    return 0;
  }

  set<C4D::Obj*, C4D::Compare> allobjs = c4d.getObjects$();

  cout << "creating event ranges" << endl;

  //C4D::EventQ event(&c4d, 0.1);
  double crossStep = 1.00001 * sin(2*M_PI/100);
  int n = 0;
  for (double angle = 0; angle < 2*M_PI; angle += 2*M_PI/100) {
    n++;
    // C4D::EventQ *event = new C4D::EventQ(&c4d, angle);
    C4D::Event *event = new C4D::Event(&c4d, new InputAngle(angle));
    c4d.current = event;
    set<C4D::Obj*, C4D::Compare> objects = c4d.getObjects$(event);
    for (set<C4D::Obj*, C4D::Compare>::iterator it = objects.begin();
         it != objects.end(); ++it)
      if ((*it)->getType() != C4D::SVF &&
          (*it)->getType() != C4D::SFV &&
          (*it)->getType() != C4D::SEE &&
          (*it)->getType() != C4D::IEF)
        (*it)->updateEvents(event, crossStep);

    set<C4D::Obj*, C4D::Compare>::iterator it = allobjs.begin();
    C4D::Obj *prev = 0;
    while (it != allobjs.end()) {
      assert ((*it)->getType() == C4D::IEF ||
              (*it)->getType() == C4D::IFF);
      if (prev != 0)
        assert(C4D::Compare()(prev, *it));
      prev = *it;
      it++;
    }
    it = allobjs.begin();
    set<C4D::Obj*, C4D::Compare>::iterator jt = objects.begin();
    prev = 0;
    while (jt != objects.end()) {
      if ((*jt)->getType() == C4D::IEF ||
          (*jt)->getType() == C4D::IFF)
        assert(allobjs.find(*jt) != allobjs.end());
      if (prev != 0)
        assert(C4D::Compare()(prev, *jt));
      prev = *jt;
      if (*jt == (C4D::Obj*)0x102968510)
        cout << "this is it" << endl;
      jt++;
    }
    jt = objects.begin();
    static int count;
    while (it != allobjs.end() || jt != objects.end())
      if (it == allobjs.end()) {
        assert((*jt)->getType() != C4D::IEF &&
               (*jt)->getType() != C4D::IFF);
        jt++;
      }
      else if (jt == objects.end()) {
        assert(!(*it)->life.contains(event));
        it++;
      }
      else if ((*jt)->getType() != C4D::IEF &&
               (*jt)->getType() != C4D::IFF) {
        if (*jt == (C4D::Obj*)0x102968510)
          cout << "this is it" << endl;
        jt++;
      }
      else if (*it != *jt) {
        set<C4D::Obj*, C4D::Compare>::iterator kt = objects.find(*it);
        assert(objects.find(*it) == objects.end());
        assert(allobjs.find(*jt) != allobjs.end());
        assert(!(*it)->life.contains(event));
        it++;
      }
      else {
        assert((*it)->life.contains(event));
        if (*jt == (C4D::Obj*)0x102968510)
          cout << "this is it" << endl;
        it++;
        jt++;
      }
  }  

  cout << "done" << endl;

  C4D::Event event(&c4d, new InputAngle(3.2672563598));

#ifdef JUSTFACES
  cout << "getting sum faces" << endl;
  vector<C4D::Face*> sumFaces = c4d.getSumFaces(&event);
  cout << "size " << sumFaces.size() << endl;
#endif
  cout << "getting objects" << endl;
  set<C4D::Obj*, C4D::Compare> objects = c4d.getObjects(&event);
  objects = c4d.getObjects(&event);
  cout << "size " << objects.size() << endl;
  for (set<C4D::Obj*, C4D::Compare>::iterator it = objects.begin();
       it != objects.end(); ++it) {
    C4D::Type type = (*it)->getType();
    if (type == C4D::SUBE) {
      C4D::SubEdge *edge = dynamic_cast<C4D::SubEdge*>(*it);
      cout << "edge " << edge << " ";
    }
    else {
      C4D::SumFace *face = dynamic_cast<C4D::SumFace*>(*it);
      cout << "face " << face << " ";
    }
    cout << (*it)->life.size() << endl;
  }
}
#endif
