#ifndef C4D_H
#define C4D_H

#include "angle.h"
#include <fstream>
#include <sstream>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <math.h>
using std::string;
using std::stringstream;
using std::vector;
using std::set;
using std::unordered_set;
using std::unordered_map;
using std::map;
using std::ifstream;
using std::istream;
using std::cout;
using std::endl;
using namespace acp;

// Macro to create one-argument primitive.
#define Primitive1(P, t1, v1)                                    \
  class P : public acp::Primitive {                              \
    t1 v1;                                                       \
  DeclareSign;                                                   \
  public:                                                        \
  P (t1 v1) : v1(v1) {}                                          \
  };

// Macro to create two-argument primitive.
#define Primitive2(P, t1, v1, t2, v2)                            \
  class P : public acp::Primitive {                              \
    t1 v1; t2 v2;                                                \
  DeclareSign;                                                   \
  public:                                                        \
  P (t1 v1, t2 v2) : v1(v1), v2(v2) {}                           \
  };

// Macro to create three-argument primitive.
#define Primitive3(P, t1, v1, t2, v2, t3, v3)                       \
  class P : public acp::Primitive {                                 \
    t1 v1; t2 v2; t3 v3;                                            \
  DeclareSign;                                                      \
  public:                                                           \
  P (t1 v1, t2 v2, t3 v3) : v1(v1), v2(v2), v3(v3) {}               \
  };

// Macro to create four-argument primitive.
#define Primitive4(P, t1, v1, t2, v2, t3, v3, t4, v4)                 \
  class P : public acp::Primitive {                                   \
    t1 v1; t2 v2; t3 v3; t4 v4;                                       \
    DeclareSign;                                                      \
  public:                                                             \
  P (t1 v1, t2 v2, t3 v3, t4 v4) : v1(v1), v2(v2), v3(v3), v4(v4) {}  \
  };

// Macro to create five-argument primitive.
#define Primitive5(P, t1, v1, t2, v2, t3, v3, t4, v4, t5, v5)           \
  class P : public acp::Primitive {                                     \
    t1 v1; t2 v2; t3 v3; t4 v4; t5 v5;                                  \
    DeclareSign;                                                        \
  public:                                                               \
  P (t1 v1, t2 v2, t3 v3, t4 v4, t5 v5)                                 \
    : v1(v1), v2(v2), v3(v3), v4(v4), v5(v5) {}                         \
  };

inline string int2string(int number)
{
  stringstream ss;//create a stringstream
  ss << number;//add number to the stream
  return ss.str();//return a string with the contents of the stream
}

static PV3<double> PV3PtoPV3d (const PV3<Parameter> &PV3P) {
  assert(PV3P.x.lb() == PV3P.x.ub() &&
         PV3P.y.lb() == PV3P.y.ub() &&
         PV3P.z.lb() == PV3P.z.ub());
  return PV3<double>(PV3P.x.lb(), PV3P.y.lb(), PV3P.z.lb());
}

#ifdef BLEEN
template<class T>
class PTR {
  T *t;
public:
  PTR () : t(0) {}
  PTR (T *t) : t(t) { incRef(); }
  PTR (const PTR &p) : t(p.t) { incRef(); }
  const PTR &operator= (const PTR &p) { p.incRef(); decRef(); t = p.t; return *this; }
  ~PTR () { decRef(); }
  void incRef () const { if (t != 0) t->incRef(); }
  void decRef () const { if (t != 0) t->decRef(); }
  operator T* () const { return t; }
  T *operator-> () const { return t; }
};
#endif

#ifndef NOTTIME
extern int objCount;
#endif

class Face4 {
  void minimize1 ();
public:
  int ns[2][4], is[2][4][3];

  Face4 () {}

  Face4 (int nF[4], int indsF[4][3],
         int nM[4], int indsM[4][3]);

  Face4 (const Face4 &f, bool mf);

  class Compare {
  public:
    bool operator() (const Face4 &a, const Face4 &b) const {
      for (int fm = 0; fm <= 1; fm++)
	for (int i = 0; i < 4; i++)
	  if (a.ns[fm][i] != b.ns[fm][i])
	    return a.ns[fm][i] < b.ns[fm][i];
      for (int fm = 0; fm <= 1; fm++)
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < a.ns[fm][i]; j++)
            if (a.is[fm][i][j] != b.is[fm][i][j])
              return a.is[fm][i][j] < b.is[fm][i][j];
      return false;
    }
  };

  bool operator== (const Face4 &b) const {
    Compare compare;
    return !compare(*this, b) && !compare(b, *this);
  }

  bool operator< (const Face4 &b) const {
    Compare compare;
    return compare(*this, b);
  }

  void reindex ();

  void minimize ();

  void swap(int &a, int &b);
  void swap(int a[3], int b[3]);
  void swapFaces(int i, int j);

  void replace (int r[2][12]);
  void replace2 (int r[2][12]);

  bool matches (int nF[4], int indsF[4][3]);
};

class C4D {
public:
  // Rotate p about the z-axis according q = (cos theta, sin theta).
  template<class N>
  static PV3<N> rotate (PV2<N> q, PV3<N> p) {
    return PV3<N>(q.x * p.x - q.y * p.y, q.x * p.y + q.y * p.x, p.z);
  }
  
  // Derivative with respect to theta.
  template<class N>
  static PV3<N> rotated (PV2<N> q, PV3<N> p) {
    return PV3<N>(-q.y * p.x - q.x * p.y, -q.y * p.y + q.x * p.x, N(0));
  }
  
  enum Type { ERROR,
              V, E, F,
              SVV, SVE, SEV, SVF, SFV, SEE,
              IEF, IFF, IFFF,
              SUBE, SUBF };

  class Obj;
  class Vert;
  class Event : public RefCnt {
    C4D *c4d;
    PTR<Angle> a;
    int refCount;
  public:
    Angle *getAngle () { return a; }
    RootAngle *getRootAngle () { return dynamic_cast<RootAngle*>((Angle *) a); }
    AnglePoly *getPoly () { return getRootAngle()->getPoly(); }
    AnglePoly::Type polyType () { return getPoly()->getType(); }

    Event (C4D *c4d, Angle *a) : c4d(c4d), a(a), refCount(0) {}
    ~Event () {}

    template<class N>
    PV2<N> getU () { return a->get<N>(); }

    bool lessThan (Event *that) {
      //#ifdef ROOT_SEPARATION
      nRootCmps++;
      //#endif
      return Angle::compare(this->a, that->a) < 0;
    }
    bool between (Event *a, Event *b) {
      return Angle::circulation(a->a, this->a, b->a) > 0;
    }
    bool identical (Event *that) { return this->a->identical(that->a); }

    // List of objects that vanish and appear at this event.
    set<Obj*> vanish, appear;

    // List of objects that need to be swapped, by twos.
    set<pair<Obj*,Obj*>> swap;
    
    int findVanish (Obj *obj) {
      if (vanish.find(obj) == vanish.end())
	return -1;
      return 0;
    }

    int findAppear (Obj *obj) {
      if (appear.find(obj) == appear.end())
	return -1;
      return 0;
    }


    // For IFFF vertices because their events are determined during the sweep.
    // Add this vanish event to v's life and add v to vanish list.
    void makeVanish (Vert *v);
    // Add this appear event to v's life and add v to appear list.
    void makeAppear (Vert *v);

    // If an IntEF is identical to a SumVV at a PolyEF root, return
    // the SumVV.  Identity!
    Vert *resolveVert (Vert *v);
  } event0, *EVENT0;
  
  class EventSet {
    bool all;
    static void sort (vector<PTR<Event> > &events);

  public:
    vector<PTR<Event> > events;
    EventSet () : events(0), all(false) {}
    EventSet (bool all) : events(0), all(all) {}
    EventSet (vector<PTR<Event> > events) : events(events), all(false) {}
    EventSet (C4D *c4d, vector<PTR<Angle>> angles)
      : events(angles.size()), all(false) {
      if (angles.size() == 2)
	Angle::compare(angles[0], angles[1]);
      for (int i = 0; i < angles.size(); i++)
        events[i] = new Event(c4d, angles[i]);
    }

    EventSet (C4D *c4d, vector<PTR<Angle>> angles, bool containsZero)
      : events(angles.size()), all(false) {
      if (angles.size() == 2)
	Angle::compare(angles[0], angles[1]);
#ifdef BLEEN
      for (int j = 0; j < angles.size(); j++) {
	bool flag = true;
	PV2<Parameter> u = Angle::s2u<Parameter>(angles[j]->getS<Parameter>(flag));
	if (j == 7)
	  cout << j << endl;
      }
#endif
      if (containsZero)
        if (angles.size() == 0)
          all = true;
        else {
          for (int i = 1; i < angles.size(); i++)
            events[i-1] = new Event(c4d, angles[i]);
          events.back() = new Event(c4d, angles.front());
        }
      else
        for (int i = 0; i < angles.size(); i++)
          events[i] = new Event(c4d, angles[i]);
    }

    bool isNull () const { return events.size() == 0 && !all; }
    bool isAll () const { return all; }
    int size () const { return events.size(); }
    bool contains (Angle *angle) const;
    bool contains (Event *event) const;
    EventSet complement () const;
    EventSet intersect (const EventSet &that) const;
    EventSet combine (const EventSet &that) const {
      return complement().intersect(that.complement()).complement();
    }

    // Expand live interval to include event if event lies within crossStep.
    // Otherwise, make a zero-length live interval at event.
    void approximateUpdate (Event *event, double crossStep);
  };    

  class Obj {
  public:
    virtual C4D *getC4D () { assert(0); return 0; }
    virtual string toString () { assert(0); return ""; }
    virtual Type getType () { assert(0); return ERROR; }
    virtual bool lessThan (Obj *o) { assert(0); return false; }
    virtual size_t hashCode () const { assert(0); return 0; }
    EventSet life;
    void updateEvents (Event *event, double crossStep) {
      life.approximateUpdate(event, crossStep);
    }
    bool contains (Event *event) {
      return life.contains(event);
    }

#ifndef NOTIME
    Obj () { objCount++; }
#endif

    // Add this object's events to the C4D eventList.
    // Add this object to vanish or appear for the event.
    // Consolidate duplicate events.
    void informEvents();
  
    // position, tangent, or normal at event angle
    template<class N>
    PV3<N> get (Event *event) {
      return std::is_same<N, Parameter>::value ? 
        (PV3<N>)calculateP(event) : std::is_same<N, PParameter>::value ? 
        (PV3<N>)calculatePP(event) : (PV3<N>)calculateMP(event);
    }
    // derivative with respect to event angle
    template<class N>
    PV3<N> getd (Event *event) {
      return std::is_same<N, Parameter>::value ? 
        (PV3<N>)calculatedP(event) : std::is_same<N, PParameter>::value ? 
        (PV3<N>)calculatedPP(event) : (PV3<N>)calculatedMP(event);
    }

    virtual PV3<Parameter> calculateP (Event *event) = 0;
    virtual PV3<PParameter> calculatePP (Event *event) = 0;
    virtual PV3<MParameter> calculateMP (Event *event) = 0;

    virtual PV3<Parameter> calculatedP (Event *event) = 0;
    virtual PV3<PParameter> calculatedPP (Event *event) = 0;
    virtual PV3<MParameter> calculatedMP (Event *event) = 0;

    // get at angle zero
    template<class N>
    PV3<N> get () { return get<N>(getC4D()->EVENT0); }

    // Get from ``outside'' ACP.
    PV3<Parameter> get$ (Event *event) {
      enable();
      PV3<Parameter> p = get<Parameter>(event);
      disable();
      return p;
    }
  };

#define Declare(C)								\
    PV3<Parameter> C ## P (Event *event) { return C<Parameter>(event); }	\
    PV3<PParameter> C ## PP (Event *event) { return C<PParameter>(event); }	\
    PV3<MParameter> C ## MP (Event *event) { return C<MParameter>(event); }	\
    template<class N>								\
    PV3<N> C (Event *event)

  class IntEdge;
  class IntEF;
  class Compare {
  public:
    bool operator() (Obj *a, Obj *b) const {
      if (a->getType() != b->getType())
        return a->getType() < b->getType();
      return a->lessThan(b);
    }

    bool operator() (IntEdge *a, IntEdge *b) const {
      return operator()((Obj *) a, (Obj *)b);
    }

    bool operator() (IntEF *a, IntEF *b) const {
      return operator()((Obj *) a, (Obj *)b);
    }

    bool operator() (Event *a, Event *b) const {
      static int count;
      if (++count == 0) {
	Contact ca = dynamic_cast<RootAngle*>(a->getAngle())->getPoly()->getContact();
	Contact::VMap vma;
	Contact ka = ca.getKey(vma);
	Contact cb = dynamic_cast<RootAngle*>(b->getAngle())->getPoly()->getContact();
	Contact::VMap vmb;
	Contact kb = cb.getKey(vmb);
	cout << "this is it" << endl;
      }
      return a->lessThan(b);
    }
  };

  class Equal {
  public:
    // bool operator() (IntEdge *const a, IntEdge *const b) const;

    bool operator() (Obj *const a, Obj *const b) const {
      Compare compare;
      return !(compare(a, b) || compare(b, a));
    }

    size_t operator() (Obj *const obj) const {
      return obj->hashCode();
    }
  };

  class Edge;

  class Vert : public Obj {
  public:
  };

  class Face;
  class SubEdge;

  class Edge : public Obj {
  protected:
  public:
    Vert *tail;
    Edge *twin;

    /*
      vector<Face *> faces;
      Edge (Vert *tail, Edge *twin) : tail(tail), twin(twin), faces(1) {
      faces.clear();
      }

      Face *&face () { assert(faces.size() == 1); return faces[0]; }

      int faceIndex (Face *face) {
      for (int i = 0; i < faces.size(); i++)
      if (faces[i] == face)
      return i;
      return -1;
      }
    */
    // Face to the left of this edge.
    Face *face;
    Edge (Vert *tail, Edge *twin) : tail(tail), twin(twin), face(0) {}

    Vert *getTail () const { return tail; }
    Edge *getTwin () { return twin; }
    Vert *getHead () { return twin->tail; }
    C4D *getC4D () { return tail ? tail->getC4D() : face->getC4D(); }
    bool contains (Vert *v) { return getTail() == v || getHead() == v; }

    // Uniquely selects either this or twin.
    bool smallestTwin () {
      // return face != 0 && (twin->face == 0 || lessThan(twin));
      return lessThan(twin);
    }

    // Verts that lie on this Edge ordered by direction of getV(event).
    // Only in smallest twin.
    vector<Obj *> verts;

    int numVerts () {
      return smallestTwin() ? verts.size() : twin->verts.size();
    }

    // Tail is -1 and head is verts.size().
    Vert *getVert (int i) {
      if (!smallestTwin())
        return twin->getVert(twin->verts.size() - 1 - i);
      if (i == -1)
        return tail;
      if (i == verts.size())
        return twin->tail;
      return dynamic_cast<Vert*>(verts[i]);
    }

    // Index at which new vertex should be inserted.
    int insertIndex (Event *event, Vert *v);

    // Add Vert to verts 
    void addVert (Event *event, Vert *v);

    // Index of v in verts or -1.
    int findVert (Vert *v);

    bool containsVert (Vert *v) { return findVert(v) != -1; }

    // Create a pair of subedges for each interval between
    // tail, elements of verts, and head.
    // Insert them into objs.
    vector<SubEdge*> makeSubEdges (set<Obj*, Compare> &objs);
    vector<SubEdge*> makeSubEdges (unordered_set<Obj*, Equal, Equal> &objs);
    vector<SubEdge*> makeSubEdges (unordered_set<Obj*, Equal, Equal> *objs);

    // Swap vertices at i1 and i2=i1+1.
    // Update vanish and appear of event to reflect subedges.
    void handleSwap (Event *event, int i1);

    // Swap v1 and v2, who must be neighbors.
    void handleSwap (Event *event, Vert *v1, Vert *v2);

    // Add or remove v from tail end.
    void handleTail(Event *event, Vert *v);

    void handleHead (Event *event, Vert *v);

    void handleRemove (Event *event, Vert *v);
    void handleInsert (Event *event, Vert *v);
    void handleReplace (Event *event, Vert *v0, Vert *v1);
    void handleReplace (Event *event, int minVanish, int MaxVanish,
                        vector<Vert*> &appear);

    // v0 and v1 must be neighbors, but not necessarily in that order.
    void handleRemove (Event *event, Vert *v0, Vert *v1);

    // Insert v0 and v1 as neighbors in that order.
    void handleInsert (Event *event, Vert *v0, Vert *v1);

    // Subedge from i to i+1 vanishes.  Update its life and event->vanish.
    void vanishSub (Event *event, int i);

    // Call vanishSubs for all subedges.
    void vanishSubs (Event *event);

    // Subedge from i to i+1 appears.  Update its life and event->vanish.
    // Also calculate all events associated with this subedge.
    void appearSub (Event *event, int i);

    // Call appearSubs for all subedges.
    void appearSubs (Event *event);
  };
  
  class Face : public Obj {
  public:
    Vert *getVert () { return getOuterLoop()[0]->getTail(); }

    bool contains (Vert *v);

    bool contains (Edge *e);

    // Point in face.
    template<class N>
    PV3<N> getFP (Event *event) { return getVert()->get<N>(event); }
    template<class N>
    PV3<N> getFPd (Event *event) { return getVert()->getd<N>(event); }
    // Get point from ``outside'' ACP.
    PV3<Parameter> getFP$ (Event *event) {
      enable();
      PV3<Parameter> p = getFP<Parameter>(event);
      disable();
      return p;
    }
    template<class N>
    PV3<N> getFP () { return getFP<N>(getC4D()->EVENT0); }

    // loop[i]->getHead() == loop[i+1]->getTail()
    // Counterclockwise when viewed from outside according to N.
    virtual const vector<Edge *> &getOuterLoop () = 0;

    // number of holes
    virtual int numInnerLoops () = 0;

    // i'th hole.  head to tail, but clockwise.
    virtual const vector<Edge *> &getInnerLoop (int i) = 0;

    C4D *getC4D () { return getVert()->getC4D(); }

    string toString ();

    // Intersection of (plane of) face with (line of) edge.
    template<class N>
    PV3<N> intersection (Event *event, Edge *e) {
      PV3<N> t = e->getTail()->get<N>(event);
      PV3<N> v = e->get<N>(event);
      PV3<N> p = getFP<N>(event);
      PV3<N> n = get<N>(event);
      
      // n * (t + v s - p) = 0
      // (n * v) s = n * (p - t)
      // s = (n * (p - t)) / (n * v)
      N s = n.dot(p - t) / n.dot(v);
      PV3<N> q = t + v * s;
      return q;
    }

      // Intersection of (plane of) face with (line of) edge.
    template<class N>
    PV3<N> intersectiond (Event *event, Edge *e) {
      PV3<N> t = e->getTail()->get<N>(event);
      PV3<N> v = e->get<N>(event);
      PV3<N> p = getFP<N>(event);
      PV3<N> n = get<N>(event);

      PV3<N> td = e->getTail()->getd<N>(event);
      PV3<N> vd = e->getd<N>(event);
      PV3<N> pd = getFPd<N>(event);
      PV3<N> nd = getd<N>(event);
      
      // n * (t + v s - p) = 0
      // (n * v) s = n * (p - t)
      // s = (n * (p - t)) / (n * v)
      // d (a/b) = d a / b - a d b / b^2
      N a = n.dot(p - t);
      N ad = nd.dot(p - t) + n.dot(pd - td);
      N b = n.dot(v);
      N bd = nd.dot(v) + n.dot(vd);
      N s = a / b;
      N sd = ad / b - a * bd / (b * b);
      // PV3 q = t + v * s;
      PV3<N> qd = td + vd * s + v * sd;
      return qd;
    }
};

  // v * n at angle 0
  class DotEF0 : public Primitive {
    Edge *e; Face *f;
  public:
    DotEF0 (Edge *e, Face *f) : e(e), f(f) {}
    DeclareSign { return e->get<N>().dot(f->get<N>()); }
  };

  // v * n
  class DotEF : public Primitive {
    Event *event; Edge *e; Face *f;
  public:
    DotEF (Event *event, Edge *e, Face *f) : event(event), e(e), f(f) {}
    DeclareSign { return e->get<N>(event).dot(f->get<N>(event)); }
  };

  // Positive if a,b,c are right-handed and angle=0.
  class CircEEE0 : public Primitive {
    Edge *a, *b, *c;
  public:
    CircEEE0 (Edge *a, Edge *b, Edge *c) : a(a), b(b), c(c) {}
    DeclareSign { return a->get<N>().cross(b->get<N>()).dot(c->get<N>()); }
  };

  // Positive if a and normals to b and c are right-handed.
  class CircEFF0 : public Primitive {
    Edge *a; Face *b; Face *c;
  public:
    CircEFF0 (Edge *a, Face *b, Face *c) : a(a), b(b), c(c) {}
    DeclareSign { return a->get<N>().cross(b->get<N>()).dot(c->get<N>()); }
  };

  class CircEFF : public Primitive {
    Event *e; Edge *a; Face *b; Face *c;
  public:
    CircEFF (Event *e, Edge *a, Face *b, Face *c) : e(e), a(a), b(b), c(c) {}
    DeclareSign { return a->get<N>(e).cross(b->get<N>(e)).dot(c->get<N>(e)); }
  };

  class CircEEF : public Primitive {
    Event *e; Edge *a; Edge *b; Face *c;
  public:
    CircEEF (Event *e, Edge *a, Edge *b, Face *c) : e(e), a(a), b(b), c(c) {}
    DeclareSign { return a->get<N>(e).cross(b->get<N>(e)).dot(c->get<N>(e)); }
  };

  // Given ab intersects a plane through e from below to above,
  // does it intersect it to the left of e's line?
  class LeftOf : public Primitive {
    Event *event;
    Edge *ab, *e;
  public:
    LeftOf (Event *event, Edge *ab, Edge *e) : event(event), ab(ab), e(e) {}
    DeclareSign {
      PV3<N> a = ab->getTail()->get<N>(event);
      PV3<N> abv = ab->get<N>(event);
      PV3<N> c = e->getTail()->get<N>(event);
      PV3<N> v = e->get<N>(event);
      return v.cross(a - c).dot(abv);
    }
  };

  // Sign of [na, nb, nc].
  Primitive4(CircFFF, Event*, event, Face*, a, Face*, b, Face*, c);

  // Sign of (v1 - v2) * (e1 x e2)
  Primitive5(LeftVEVE, Event*, event, Vert*, v1, Edge*, e1, Vert*, v2, Edge*, e2);

#ifndef HOMOTOPY_IDENTITY
  // u and v lie on e.
  // Negative if u precedes v in the direction of e->getV().
  class OrderUVE : public Primitive {
    Event *event; Vert *u; Vert *v; Edge *e;
  public:
    OrderUVE (Event *event, Vert *u, Vert *v, Edge *e)
      : event(event), u(u), v(v), e(e) {}
    DeclareSign {
      static int count = 0;
      ++count;
      C4D *c4d = 0;
      bool isIdentity = false;
      Angle *angle = event->getAngle();
      if (angle->getType() == Angle::IMPLICIT || 
          angle->getType() == Angle::EXPLICIT) {
        Contact contact;
        if (u->getType() == IEF && v->getType() == IEF && 
            SVE <= e->getType() && e->getType() <= SEV) {
          IntEF *u2 = dynamic_cast<IntEF*>(u);
          IntEF *v2 = dynamic_cast<IntEF*>(v);
          SumEdge *e2 = dynamic_cast<SumEdge*>(e);
          assert(u2->e == e2 || u2->e == e2->getTwin());
          assert(v2->e == e2 || v2->e == e2->getTwin());
          c4d = u2->getC4D();
          contact = c4d->getContact(u2, v2);
        }
        else if (u->getType() == IFFF && v->getType() == IFFF && 
                 e->getType() == IFF) {
          IntFFF *u2 = dynamic_cast<IntFFF*>(u);
          IntFFF *v2 = dynamic_cast<IntFFF*>(v);
          IntEdge *e2 = dynamic_cast<IntEdge*>(e);
          assert(u2->hasFace(e2->getFace()) && 
                 u2->hasFace(e2->getTwin()->getFace()));
          assert(v2->hasFace(e2->getFace()) && 
                 v2->hasFace(e2->getTwin()->getFace()));
          c4d = u2->getC4D();
          contact = c4d->getContact(u2, v2);
        }
        else
          assert(0);

        RootAngle *root = dynamic_cast<RootAngle*>(angle);
        isIdentity = c4d->commonRoot(contact, root);
      }

      if (isIdentity)
        return ((u->getd<N>(event) - v->getd<N>(event)).dot(e->get<N>(event)) +
                (u->get<N>(event) - v->get<N>(event)).dot(e->getd<N>(event)));
                                                    
      return (u->get<N>(event) - v->get<N>(event)).dot(e->get<N>(event));
    }
  };
#else
  class OrderUVE {
    Event *event; Vert *u; Vert *v; Edge *e;
  public:
    OrderUVE (Event *event, Vert *u, Vert *v, Edge *e)
      : event(event), u(u), v(v), e(e) {}
    operator int () {
      int s = OrderUVEH(event, u, v, e);
      if (s != 0)
	return s;
      return OrderUVED(event, u, v, e);
    }
  };
  class OrderUVEH : public Primitive {
    Event *event; Vert *u; Vert *v; Edge *e;
  public:
    OrderUVEH (Event *event, Vert *u, Vert *v, Edge *e)
      : event(event), u(u), v(v), e(e) {}
    DeclareSign {
      return (u->get<N>(event) - v->get<N>(event)).dot(e->get<N>(event));
    }
  };
  class OrderUVED : public Primitive {
    Event *event; Vert *u; Vert *v; Edge *e;
  public:
    OrderUVED (Event *event, Vert *u, Vert *v, Edge *e)
      : event(event), u(u), v(v), e(e) {}
    DeclareSign {
      return ((u->getd<N>(event) - v->getd<N>(event)).dot(e->get<N>(event)) +
	      (u->get<N>(event) - v->get<N>(event)).dot(e->getd<N>(event)));
    }
  };
#endif

  class EdgeAB;
  class FaceABC;
  class SumFace;

  // Input vertex belonging to obstacle or robot.
  class VertXYZ;
  // verts[0] is obstacle vertices, verts[1] is robot vertices.
  vector<VertXYZ *> verts[2];
  class VertXYZ : public Vert, public Object<PV3> {
    C4D *c4d;
    bool moving;
  public:
    VertXYZ (C4D *c4d, double x, double y, double z, bool moving)
      : c4d(c4d), id(c4d->verts[moving].size()), moving(moving),
	Object<PV3>(PV3<double>(x, y, z), true) {
      c4d->verts[moving].push_back(this);
    }

    VertXYZ (C4D *c4d, double x, double y, double z, bool moving, bool constant)
      : c4d(c4d), id(c4d->verts[moving].size()), moving(moving), 
        Object<PV3>(PV3<double>(x, y, z), true) {
      c4d->verts[moving].push_back(this);
    }

    VertXYZ (C4D *c4d, VertXYZ *v)
      : c4d(c4d), moving(v->moving), id(c4d->verts[moving].size()),
        // Object<PV3>(PV3PtoPV3d(v->get<Parameter>()), false) {
        Object<PV3>(v->getApproxMid(), false) {
      c4d->verts[moving].push_back(this);
    }
      
    C4D *getC4D () { return c4d; }
    Type getType () { return V; }
    const int id;
    int getID () const { return id; }
    int getIDM () const { return 2 * id + isMoving(); }

    void getContact (int i, Contact &contact) {
      int fm = isMoving();
      contact.is[fm][i][contact.ns[fm][i]++] = id;
    }

    bool lessThan (Obj *o) {
      return getIDM() < dynamic_cast<VertXYZ*>(o)->getIDM();
    }
    size_t hashCode () const { return getIDM(); }
    bool isMoving () const { return moving; }

    template<class N>
    PV3<N> get () { return Object<PV3>::get<N>(); }

    template<class N>
    PV3<N> get (Event *event) { return Vert::get<N>(event); }

    Declare(calculate) {
      return event != getC4D()->EVENT0 && moving ?
        rotate(event->getU<N>(), get<N>()) : get<N>();
    }

    Declare(calculated) {
      return event != getC4D()->EVENT0 && moving ?
        rotated(event->getU<N>(), get<N>()) : PV3<N>(N(0), N(0), N(0));
    }

    // Edges with e->tail == this
    // Sorted counterclockwise by call to makeCone.
    vector<Edge *> edges;
    EdgeAB *getEdge (Vert *head);
    void addEdge (Edge *edge);

    FaceABC *getFace (Vert *b, Vert *c);

    // Edges of convex hull in a neighborhood of this vertex.
    // In counterclockwise order as viewed from outside.
    // Empty if vertex is concave.
    vector<Edge *> cone;
    // Initialize cone.
    void makeCone ();
    // Use cone to tell if the vertex and face f are compatible.
    bool isCompatible (Event *event, FaceABC *f);

    string toString () { return int2string(id); }

    SumFace *getSum (FaceABC *f);
  };

  class Int2 {
  public:
    int i[2];
    Int2 () { i[0] = i[1] = 0; }
    Int2 (const int *ii) { i[0] = ii[0]; i[1] = ii[1]; }
    Int2 (int t, int h) {
      i[0] = t;
      i[1] = h;
    }
    Int2 (int t, int h, bool sort) {
      if (t < h) {
        i[0] = t;
        i[1] = h;
      }
      else {
        i[0] = h;
        i[1] = t;
      }
    }
    bool operator() (const Int2 &a, const Int2 &b) const {
      if (a.i[0] != b.i[0])
	return a.i[0] < b.i[0];
      else
	return a.i[1] < b.i[1];
    }
  };

  // Input edge.
  vector<EdgeAB *> edges[2];
  map<Int2, EdgeAB*, Int2> i2edge[2];

  class EdgeAB : public Edge, public Object<PV3> {
    DeclareCalculate(PV3) {
      return(getHead()->Object<PV3>::get<N>() - getTail()->Object<PV3>::get<N>());
    }
    EdgeAB (VertXYZ *tail, EdgeAB *twin)
      : Edge(tail, twin), id(getC4D()->edges[isMoving()].size()) {
      getC4D()->edges[isMoving()].push_back(this);
    }
  public:
    Type getType () { return E; }
    void getContact (int i, Contact &contact) {
      getTail()->getContact(i, contact);
      getHead()->getContact(i, contact);
    }

    const int id;
    int getID () const { return id; }
    int getIDM () const { return 2 * id + isMoving(); }
    bool lessThan (Obj *o) {
      return getIDM() < dynamic_cast<EdgeAB*>(o)->getIDM();
    }
    size_t hashCode () const { return getIDM(); }

    EdgeAB (VertXYZ *tail, VertXYZ *head)
      : Edge(tail, new EdgeAB(head, this)),
        id(getC4D()->edges[isMoving()].size()) {
      assert(this->tail != 0);
      getC4D()->edges[isMoving()].push_back(this);
      getC4D()->i2edge[isMoving()][Int2(tail->id, head->id, true)] = this;
      tail->addEdge(this);
      getHead()->addEdge(getTwin());
    }

    static EdgeAB *make (C4D *c4d, int fm, Int2 i2) {
      EdgeAB *&e = c4d->i2edge[fm][i2];
      if (e != 0)
        return e;
      return new EdgeAB(c4d->verts[fm][i2.i[0]], 
                        c4d->verts[fm][i2.i[1]]);
    }

    template<class N>
    PV3<N> get () { return Object<PV3>::get<N>(); }

    template<class N>
    PV3<N> get (Event *event) { return Edge::get<N>(event); }

    Declare(calculate) {
      return event != getC4D()->EVENT0 && isMoving() ?
        rotate(event->getU<N>(), Object<PV3>::get<N>()) : Object<PV3>::get<N>();
    }

    Declare(calculated) {
      return event != getC4D()->EVENT0 && isMoving() ?
        rotated(event->getU<N>(), Object<PV3>::get<N>()) : PV3<N>(N(0), N(0), N(0));
    }

    bool isMoving () const { return getTail()->isMoving(); }

    bool isConvex () { return CircEFF0(this, face, twin->face) > 0; }

    bool isCompatible (Event *event, EdgeAB *that) {
      // return (that->getV(event).dot(face->getN(event)) > 0 &&
      // that->getV(event).dot(twin->face->getN(event)) < 0);
      return (DotEF(event, that, face) > 0 &&
              DotEF(event, that, twin->face) < 0);
    }

    SumFace *getSum (EdgeAB *that);

    Primitive3(IsOutside, Event*, event, Face*, f, EdgeAB*, e);

    // Face f is externally tangent to this edge with respect to its
    // left face.  Just for checking.
    bool isOutside (Event *event, Face *f) {
      return IsOutside(event, f, this) > 0;
      // return f->getN(event).cross(getV(event)).dot(face->getN(event)) > 0;
    }

    VertXYZ *getTail () const {
      return dynamic_cast<VertXYZ*>(Edge::getTail());
    }
    VertXYZ *getHead () {
      return dynamic_cast<VertXYZ*>(Edge::getHead());
    }
    EdgeAB *getTwin () {
      return dynamic_cast<EdgeAB*>(Edge::getTwin());
    }
    FaceABC *getFace () {
      return dynamic_cast<FaceABC*>(face);
    }

    FaceABC *getFace (VertXYZ *v) {
      if (getFace()->contains(v))
        return getFace();
      if (getTwin()->getFace()->contains(v))
        return getTwin()->getFace();
      return 0;
    }
  };

  class PolyEF;

  class Int3 {
  public:
    int i[3];
    Int3 () { i[0] = i[1] = i[2] = 0; }
    Int3 (const int *iii) { i[0] = iii[0]; i[1] = iii[1]; i[2] = iii[2]; }
    Int3 (int a, int b, int c) {
      i[0] = a;
      i[1] = b;
      i[2] = c;
    }
    Int3 (int a, int b, int c, bool sort) {
      if (a < b) {
        i[0] = a;
        i[1] = b;
      }
      else {
        i[0] = b;
        i[1] = a;
      }
      int j = 2;
      while (j > 0 && i[j-1] > c) {
        i[j] = i[j-1];
        j--;
      }
      i[j] = c;
    }

    bool operator() (const Int3 &a, const Int3 &b) const {
      if (a.i[0] != b.i[0])
	return a.i[0] < b.i[0];
      else if (a.i[1] != b.i[1])
	return a.i[1] < b.i[1];
      else
	return a.i[2] < b.i[2];
    }
  };

  // Input face.
  vector<FaceABC *> faces[2];
  map<Int3, FaceABC*, Int3> i3face[2];

  class FaceABC : public Face, public Object<PV3> {
    vector<Edge *> edges;
    friend class PolyEF;
    DeclareCalculate(PV3) {
      return(edges[0]->get<N>().cross(edges[1]->get<N>()));
    }
  public:
    Type getType () { return F; }
    void getContact (int i, Contact &contact) {
      for (int j = 0; j < 3; j++)
        getVert(j)->getContact(i, contact);
    }

    const int id;
    int getID () const { return id; }
    int getIDM () const { return 2 * id + isMoving(); }
    bool lessThan (Obj *o) {
      return getIDM() < dynamic_cast<FaceABC*>(o)->getIDM();
    }
    size_t hashCode () const { return getIDM(); }

    FaceABC (VertXYZ *a, VertXYZ *b, VertXYZ *c)
      : id(a->getC4D()->faces[a->isMoving()].size()) {
      VertXYZ *abc[] = { a, b, c };
      for (int i = 0; i < 3; i++) {
        int j = (i+1)%3;
        Edge *eij = abc[i]->getEdge(abc[j]);
        if (eij == 0)
          eij = new EdgeAB(abc[i], abc[j]);
        eij->face = this;
        edges.push_back(eij);
      }
      a->getC4D()->faces[isMoving()].push_back(this);
      a->getC4D()->i3face[isMoving()][Int3(a->id, b->id, c->id, true)] = this;
    }

    static FaceABC *make (C4D *c4d, int fm, Int3 i3) {
      FaceABC *&f = c4d->i3face[fm][i3];
      if (f != 0)
        return f;
      return new FaceABC(c4d->verts[fm][i3.i[0]], 
                         c4d->verts[fm][i3.i[1]],
                         c4d->verts[fm][i3.i[2]]);
    }

    bool isMoving () const {
      return dynamic_cast<VertXYZ*>(edges[0]->getTail())->isMoving();
    }

    template<class N>
    PV3<N> get () { return Object<PV3>::get<N>(); }

    template<class N>
    PV3<N> get (Event *event) { return Face::get<N>(event); }

    Declare(calculate) {
      return event != getC4D()->EVENT0 && isMoving() ?
        rotate(event->getU<N>(), Object<PV3>::get<N>()) : Object<PV3>::get<N>();
    }

    Declare(calculated) {
      return event != getC4D()->EVENT0 && isMoving() ?
        rotated(event->getU<N>(), Object<PV3>::get<N>()) : PV3<N>(N(0), N(0), N(0));
    }

    const vector<Edge *> &getOuterLoop () { return edges; }
    int numInnerLoops () { return 0; }
    const vector<Edge *> &getInnerLoop (int i) { assert(0); return edges; }

    VertXYZ *getVert (int i) { return dynamic_cast<VertXYZ*>(edges[i]->tail); }
  };

  class PolyVF;

  static void sort (VertXYZ *vs[], int n) {
    for (int i = 1; i < n; i++) {
      VertXYZ* v = vs[i];
      int j = i;
      while (j > 0 && v->lessThan(vs[j-1])) {
        vs[j] = vs[j-1];
        j--;
      }
      vs[j] = v;
    }
  }

  class PolyEF : public AnglePoly {
    VertXYZ *e[2], *f[3];
    DeclareCalculate(Poly2);
  public:
    Type getType () { return EF2; }

    void getContact (Contact &contact) {
      for (int i = 0; i < 2; i++)
        e[i]->getContact(0, contact);
      for (int i = 0; i < 3; i++)
        f[i]->getContact(0, contact);
    }

    PolyEF (EdgeAB *edge, FaceABC *face) {
      e[0] = edge->getTail();
      e[1] = edge->getHead();
      sort(e, 2);
      f[0] = dynamic_cast<VertXYZ*>(face->edges[0]->getTail());
      f[1] = dynamic_cast<VertXYZ*>(face->edges[1]->getTail());
      f[2] = dynamic_cast<VertXYZ*>(face->edges[2]->getTail());
      sort(f, 3);
    }

    PolyEF (VertXYZ* e0, VertXYZ* e1, FaceABC *face) {
      e[0] = e0;
      e[1] = e1;
      sort(e, 2);
      f[0] = dynamic_cast<VertXYZ*>(face->edges[0]->getTail());
      f[1] = dynamic_cast<VertXYZ*>(face->edges[1]->getTail());
      f[2] = dynamic_cast<VertXYZ*>(face->edges[2]->getTail());
      sort(f, 3);
    }

    PolyEF (EdgeAB *edge, VertXYZ* f0, VertXYZ* f1, VertXYZ* f2) {
      e[0] = edge->getTail();
      e[1] = edge->getHead();
      sort(e, 2);
      f[0] = f0;
      f[1] = f1;
      f[2] = f2;
      sort(f, 3);
    }

    PolyEF (VertXYZ* e0, VertXYZ* e1, VertXYZ* f0, VertXYZ* f1, VertXYZ* f2) {
      e[0] = e0;
      e[1] = e1;
      sort(e, 2);
      f[0] = f0;
      f[1] = f1;
      f[2] = f2;
      sort(f, 3);
    }

    bool identical (AnglePoly *thatAnglePoly);

    friend class PolyVF;

    static EventSet getNegative (EdgeAB* e, FaceABC* f);

    Primitive2(Sign, EdgeAB*, e, FaceABC*, f);
    //int sign () { return Sign(e, f); }

    Primitive2(HasRoots, EdgeAB*, e, FaceABC*, f);
    static bool hasRoots (EdgeAB *e, FaceABC *f) {
      return HasRoots(e, f) < 0;
    }

    bool contains (Vert *vert);
    bool contains (Edge *edge) {
      return contains(edge->tail) && contains(edge->twin->tail);
    }
    bool contains (SumFace *face);
  };

  class SumEdge;
  //class IntEdge;

  class SumVV;
  vector<SumVV *> sumVVs;
  map<Int2, SumVV*, Int2> i2sumVV;
  class SumVV : public Vert {
    SumVV (VertXYZ *a, VertXYZ *b) {
      vv[0] = a;
      vv[1] = b;
      a->getC4D()->sumVVs.push_back(this);
      get(a, b) = this;
    }

  public:
    Type getType () { return SVV; }
    void getContact (int i, Contact &contact) {
      vv[0]->getContact(i, contact);
      vv[1]->getContact(i, contact);
    }

    bool lessThan (Obj *o) {
      SumVV *that = dynamic_cast<SumVV*>(o);
      // if (vv[0] != that->vv[0])
      if (vv[0]->lessThan(that->vv[0]) || that->vv[0]->lessThan(vv[0]))
        return vv[0]->lessThan(that->vv[0]);
      return vv[1]->lessThan(that->vv[1]);
    }
    size_t hashCode () const {
      return vv[0]->hashCode() + 1804289383 * vv[1]->hashCode();
    }
      
    /*
    static int getID (VertXYZ *a, VertXYZ *b) {
      return a->id + b->id * a->getC4D()->verts[0].size();
    }
    int getID () { return getID(vv[0], vv[1]); }
    */

    static SumVV *&get (VertXYZ *a, VertXYZ *b) {
      assert(!a->isMoving());
      assert(b->isMoving());
      return a->getC4D()->i2sumVV[Int2(a->id, b->id)];
    }

    static SumVV *make (VertXYZ *a, VertXYZ *b) {
      SumVV *&vv = get(a, b);
      if (vv != 0)
        return vv;
      return new SumVV(a, b);
    }

    VertXYZ *vv[2];

    Declare(calculate) {
      return vv[0]->get<N>(event) + vv[1]->get<N>(event);
    }

    Declare(calculated) {
      return vv[0]->getd<N>(event) + vv[1]->getd<N>(event);
    }

    template<class N>
    PV3<N> get () {
      return vv[0]->get<N>() + vv[1]->get<N>();
    }

    string toString () { return vv[0]->toString() + "+" + vv[1]->toString(); }

    C4D *getC4D () { return vv[0]->getC4D(); }

    vector<IntEdge*> intEdges;

    void addIntEdge (IntEdge *e) {
      assert(this == e->tail);
      for (int i = 0; i < intEdges.size(); i++)
        if (intEdges[i] == e)
          assert(0);
      intEdges.push_back(e);
    }

    void removeIntEdge (IntEdge *e) {
      assert(this == e->tail);
      for (int i = 0; i < intEdges.size(); i++)
        if (intEdges[i] == e) {
          intEdges.erase(intEdges.begin() + i);
          return;
        }
      assert(0);
    }

    void handleVF (Event *event, SumFace *fSweep, SumEdge *eHit);
  };

  class PolyEE;

  class SumEdge : public Edge {
  protected:
  public:
    VertXYZ *v;
    EdgeAB *e;
    friend class PolyEE;
    vector<SumFace*> faces;

    void getContact (int i, Contact &contact) {
      v->getContact(i, contact);
      e->getContact(i, contact);
    }

    bool lessThan (Obj *o) {
      SumEdge *that = dynamic_cast<SumEdge*>(o);
      // if (v != that->v)
      if (v->lessThan(that->v) || that->v->lessThan(v))
        return v->lessThan(that->v);
      return e->lessThan(that->e);
    }
    size_t hashCode () const {
      return v->hashCode() + 846930886 * e->hashCode();
    }

    SumEdge (Vert *tail, Edge *twin, VertXYZ *v, EdgeAB *e)
      : Edge(tail, twin), v(v), e(e) {}

    SumVV *getTail () { return dynamic_cast<SumVV*>(Edge::getTail()); }
    SumVV *getHead () { return dynamic_cast<SumVV*>(Edge::getHead()); }
    SumEdge *getTwin () { return dynamic_cast<SumEdge*>(Edge::getTwin()); }
    SumFace *getFace () { return dynamic_cast<SumFace*>(face); }

    int faceIndex (SumFace *face) {
      for (int i = 0; i < faces.size(); i++)
        if (faces[i] == face)
          return i;
      return -1;
    }

    void addFace (SumFace *face) {
      assert(faceIndex(face) == -1);
      faces.push_back(face);
    }

    void removeFace (SumFace *face) {
      int i = faceIndex(face);
      assert(i != -1);
      faces.erase(faces.begin() + i);
    }

    void clearFaces () {
      faces.clear();
    }

    const vector<SumFace*> &getFaces () { return faces; }

    bool hasFace (SumFace *face) {
      return faceIndex(face) != -1 || getTwin()->faceIndex(face) != -1;
    }

    Declare(calculate) { return e->get<N>(event); }
    Declare(calculated) { return e->getd<N>(event); }
    template<class N>
    PV3<N> get () { return e->get<N>(); }

    void handleEFF (Event *event);
    void handleEFF (Event *event, SumFace *f[2], bool f0parallel);
    void handleEE (Event *event, SumEdge *that);

    bool intersects (Event *event, SumFace *face);
  };
    
  class SumVE;
  vector<SumVE *> sumVEs;
  map<Int2, SumVE*, Int2> i2sumVE;
  class SumVE : public SumEdge {
    SumVE (VertXYZ *v, EdgeAB *e)
      : SumEdge(SumVV::make(v, e->getTail()), new SumVE(v, e->getTwin(), this),
                v, e) {
      v->getC4D()->sumVEs.push_back(this);
      get(v, e) = this;
    }

    SumVE (VertXYZ *v, EdgeAB *e, Edge *twin)
      : SumEdge(SumVV::make(v, e->getTail()), twin, v, e) {
      v->getC4D()->sumVEs.push_back(this);
      get(v, e) = this;
    }

  public:
    Type getType () { return SVE; }
    /*
    static int getID (VertXYZ *v, EdgeAB *e) {
      return v->id + e->id * v->getC4D()->verts[0].size();
    }
    int getID () { return getID(v, e); }
    */
    static SumVE *&get (VertXYZ *v, EdgeAB *e) {
      // return v->getC4D()->sumVEs[getID(v, e)];
      return v->getC4D()->i2sumVE[Int2(v->id, e->id)];
    }

    static SumVE *make (VertXYZ *v, EdgeAB *e) {
      SumVE *&ve = get(v, e);
      if (ve != 0)
        return ve;
      return new SumVE(v, e);
    }
  };
  
  class SumEV;
  vector<SumEV *> sumEVs;
  map<Int2, SumEV*, Int2> i2sumEV;
  class SumEV : public SumEdge {
    SumEV (EdgeAB *e, VertXYZ *v)
      : SumEdge(SumVV::make(e->getTail(), v), new SumEV(e->getTwin(), v, this),
                v, e) {
      if (this == (C4D::SumEV *) 0x69c860)
	cout << "SumEV this is it";
      v->getC4D()->sumEVs.push_back(this);
      get(e, v) = this;
    }

    SumEV (EdgeAB *e, VertXYZ *v, Edge *twin)
      : SumEdge(SumVV::make(e->getTail(), v), twin, v, e) {
      if (this == (C4D::SumEV *) 0x69c860)
	cout << "SumEV this is it";
      v->getC4D()->sumEVs.push_back(this);
      get(e, v) = this;
    }

  public:
    Type getType () { return SEV; }
    /*
    static int getID (EdgeAB *e, VertXYZ *v) {
      return v->id + e->id * v->getC4D()->verts[1].size();
    }
    int getID () { return getID(e, v); }
    */
    static SumEV *&get (EdgeAB *e, VertXYZ *v) {
      // return v->getC4D()->sumEVs[getID(e, v)];
      return v->getC4D()->i2sumEV[Int2(e->id, v->id)];
    }

    static SumEV *make (EdgeAB *e, VertXYZ *v) {
      SumEV *&ev = get(e, v);
      if (ev != 0)
        return ev;
      if (ev == (C4D::Obj *) 0x69c860)
	cout << "SumEV this is it" << endl;
      return new SumEV(e, v);
    }
  };
  
  class SumFace : public Face {
  protected:
    vector<Edge *> edges;
  public:
    virtual void getContact (int i, Contact &contact) = 0;

    const vector<Edge *> &getOuterLoop () { return edges; }
    int numInnerLoops () { return 0; }
    const vector<Edge *> &getInnerLoop (int i) { assert(0); return edges; }

    SumVV *getVert () { return dynamic_cast<SumVV*>(Face::getVert()); }

    void initEdges () {
      for (int i = 0; i < edges.size(); i++)
        getEdge(i)->addFace(this);
    }

    void updateEdgeLives () {
      for (int i = 0; i < edges.size(); i++)
        edges[i]->life = life.combine(edges[i]->life);
    }          

#ifndef HOMOTOPY_IDENTITY
    Primitive3(Side, SumFace*, face, Event*, event, Vert*, v);
#else
    Primitive3(SideH, SumFace*, face, Event*, event, Vert*, v);
    Primitive3(SideD, SumFace*, face, Event*, event, Vert*, v);

    class Side {
      SumFace *face; Event *event; Vert *v;
    public:
      Side (SumFace *face, Event *event, Vert *v) :
	face(face), event(event), v(v) {}
      operator int () {
	Type vType = v->getType();
	if (vType == SVV) {
	  for (int i = 0; i < face->edges.size(); i++)
	    if (v == face->edges[i]->getTail())
	      return 0;
	}
	else if (vType == IEF) {
	  IntEF *ef = dynamic_cast<IntEF*>(v);
	  if (ef->f == face)
	    return 0;
	  for (int i = 0; i < face->edges.size(); i++)
	    if (ef->e == face->edges[i] || ef->e == face->edges[i]->twin)
	      return 0;
	}
	else
	  assert(0);

	int s = SideH(face, event, v);
	if (s != 0)
	  return s;
	return SideD(face, event, v);
      }
    };
#endif

    // positive if v is to the outer side of this face (same as normal)
    // zero of v is a vertex of this face
    // negative otherwise.
    int side (Event *event, Vert *v) { return Side(this, event, v); }

#ifdef BLEEN
    // True if interior e transversely intersects plane of this face
    // from inside to outside.
    bool intersectsPlane (Event *event, Edge *e) {
      return side(event, e->getTail()) < 0 && side(event, e->getHead()) > 0;
    }
#endif

    // True if that intersection is in interior of this face.
    bool intersectsInside (Event *event, Edge *e);

    // True if line through v parallel to e intersects inside face.
    // No condition on direction of e.
    bool intersectsInside (Event *event, Vert *v, Edge *e);

#ifdef BLEEN
    // True if interior e transversely intersects interior of this
    // face from inside to outside.
    bool intersectsX (Event *event, Edge *e) {
      return intersectsPlaneX(event, e) && intersectsInside(event, e);
    }
#endif

    // These are the twins of the intersection edges that belong to
    // this face so they can be looked up by the other face.
    // intEdges:  all that are known.
    // liveEdges:  ones that are live for the current event.
    set<IntEdge *, Compare> liveEdges;

    // Find live intersection edge with that face.
    IntEdge *findLiveEdge (SumFace *that);

    // Get the lifetime of this's intersection with that.
    EventSet getIntLife (SumFace *that);

    // Find all intersections of IntEdges
    void intersectEdges (Event *event, bool doHandle=false);

    int numEdges () { return edges.size(); }

    SumEdge *getEdge (int i) { return dynamic_cast<SumEdge*>(edges[i]); }

    SumVV *commonVert (SumFace *that) {
      for (int i = 0; i < edges.size(); i++)
        if (that->contains(edges[i]->tail))
          return dynamic_cast<SumVV*>(edges[i]->tail);
      return 0;
    }

    bool neighborOf (SumFace *that) {
      for (int i = 0; i < edges.size(); i++) {
        if (getEdge(i)->faceIndex(that) != -1)
          return true;
        if (getEdge(i)->getTwin()->faceIndex(that) != -1)
          return true;
      }
      return false;
    }

    bool neighborOf2 (SumFace *that) {
      for (int i = 0; i < this->edges.size(); i++)
	for (int j = 0; j < that->edges.size(); j++)
	  if (this->edges[i] == that->edges[j] ||
	      this->edges[i] == that->edges[j]->twin)
	    return true;
      return false;
    }

    // Add this SumFace's generators to fixed and moving.  Check for
    // duplicates.  Stop at 3.  Return false if about to add 4th.
    virtual bool generators (VertXYZ *fixed[3], VertXYZ *moving[3]) = 0;

    bool addGenerator (VertXYZ *gen, VertXYZ *gens[3]) {
      return ((gens[0] == 0 && (gens[0] = gen) != 0) ||
              (gens[1] == 0 &&
               (gens[0] == gen || (gens[1] = gen) != 0)) ||
              (gens[2] == 0 &&
               (gens[0] == gen || gens[1] == gen || (gens[2] = gen) != 0)) ||
              (gens[0] == gen || gens[1] == gen || gens[2] == gen));
    }

    PolyEF *related (SumFace *that) {
      VertXYZ *fixed[] = { 0, 0, 0 };
      VertXYZ *moving[] = { 0, 0, 0 };
      if (!this->generators(fixed, moving) ||
          !that->generators(fixed, moving) ||
          (fixed[2] != 0 && moving[2] != 0))
        return 0;
      if (fixed[2] == 0)
        return new PolyEF(fixed[0], fixed[1], moving[0], moving[1], moving[2]);
      else
        return new PolyEF(moving[0], moving[1], fixed[0], fixed[1], fixed[2]);
    }

    void handleEF (Event *event, IntEdge *intE, SumEdge *sumE);

    bool hasCoplanarEdge (PolyEF *poly) {
      for (int i = 0; i < edges.size(); i++)
	if (poly->contains(edges[i]))
	  return true;
      return false;
    }

    // IntEFs that have been "next" to it on the same SumFace.
    set<IntEF *, Compare> nexts;
    bool solved (IntEF *v) {
      if (nexts.find(v) != nexts.end())
        return true;
      nexts.insert(v);
      return false;
    }

    SumFace () { intEdge = 0; }

    IntEdge *intEdge;
  };
    
  class SumVF;
  vector<SumVF *> sumVFs;
  map<Int2, SumVF*, Int2> i2sumVF;
  class SumVF : public SumFace {
    SumVF (VertXYZ *v, FaceABC *f) : v(v), f(f) {
      const vector<Edge *> &fedges = f->getOuterLoop();
      for (int i = 0; i < fedges.size(); i++) {
        edges.push_back(SumVE::make(v, dynamic_cast<EdgeAB*>(fedges[i])));
        // edges.back()->face = this;
      }
      v->getC4D()->sumVFs.push_back(this);
      get(v, f) = this;
    }
  public:
    VertXYZ *v;
    FaceABC *f;

    void getContact (int i, Contact &contact) {
      v->getContact(i, contact);
      f->getContact(i, contact);
    }

    Type getType () { return SVF; }
    bool lessThan (Obj *o) {
      SumVF *that = dynamic_cast<SumVF*>(o);
      // if (v != that->v)
      if (v->lessThan(that->v) || that->v->lessThan(v))
        return v->lessThan(that->v);
      return f->lessThan(that->f);
    }
    size_t hashCode () const {
      return v->hashCode() + 1681692777 * f->hashCode();
    }

    /*
    static int getID (VertXYZ *v, FaceABC *f) {
      return v->id + f->id * v->getC4D()->verts[0].size();
    }
    int getID () { return getID(v, f); }
    */

    static SumVF *&get (VertXYZ *v, FaceABC *f) {
      assert(!v->isMoving());
      // return v->getC4D()->sumVFs[getID(v, f)];
      return v->getC4D()->i2sumVF[Int2(v->id, f->id)];
    }

    static SumVF *make (VertXYZ *v, FaceABC *f) {
      SumVF *vf = get(v, f);
      if (vf != 0)
        return vf;
      return new SumVF(v, f);
    }

    Declare(calculate) { return f->get<N>(event); }
    Declare(calculated) { return f->getd<N>(event); }

    friend class PolyVF;

    VertXYZ *getV () { return v; }
    FaceABC *getF () { return f; }

    bool generators (VertXYZ *fixed[3], VertXYZ *moving[3]) {
      return (addGenerator(v, fixed) &&
              addGenerator(f->getVert(0), moving) &&
              addGenerator(f->getVert(1), moving) &&
              addGenerator(f->getVert(2), moving));
    }
  };
  
  class SumFV;
  vector<SumFV *> sumFVs;
  map<Int2, SumFV*, Int2> i2sumFV;
  class SumFV : public SumFace {
    SumFV (FaceABC *f, VertXYZ *v) : v(v), f(f) {
      const vector<Edge *> &fedges = f->getOuterLoop();
      for (int i = 0; i < fedges.size(); i++) {
        edges.push_back(SumEV::make(dynamic_cast<EdgeAB*>(fedges[i]), v));
        // edges.back()->face = this;
      }
      v->getC4D()->sumFVs.push_back(this);
      get(f, v) = this;
    }
  public:
    VertXYZ *v;
    FaceABC *f;

    void getContact (int i, Contact &contact) {
      v->getContact(i, contact);
      f->getContact(i, contact);
    }

    Type getType () { return SFV; }
    bool lessThan (Obj *o) {
      SumFV *that = dynamic_cast<SumFV*>(o);
      // if (v != that->v)
      if (v->lessThan(that->v) || that->v->lessThan(v))
        return v->lessThan(that->v);
      return f->lessThan(that->f);
    }
    size_t hashCode () const {
      return v->hashCode() + 1714636915 * f->hashCode();
    }

    /*
    static int getID (FaceABC *f, VertXYZ *v) {
      return v->id + f->id * v->getC4D()->verts[1].size();
    }
    int getID () { return getID(f, v); }
    */

    static SumFV *&get (FaceABC *f, VertXYZ *v) {
      return v->getC4D()->i2sumFV[Int2(f->id, v->id)];
    }
    
    static SumFV *make (FaceABC *f, VertXYZ *v) {
      SumFV *fv = get(f, v);
      if (fv != 0)
        return fv;
      return new SumFV(f, v);
    }

    Declare(calculate) { return f->get<N>(event); }
    Declare(calculated) { return f->getd<N>(event); }

    friend class PolyVF;

    VertXYZ *getV () { return v; }
    FaceABC *getF () { return f; }

    bool generators (VertXYZ *fixed[3], VertXYZ *moving[3]) {
      return (addGenerator(v, moving) &&
              addGenerator(f->getVert(0), fixed) &&
              addGenerator(f->getVert(1), fixed) &&
              addGenerator(f->getVert(2), fixed));
    }
  };
  
  class SumEE;
  vector<SumEE *> sumEEs;
  map<Int2, SumEE*, Int2> i2sumEE;
  class SumEE : public SumFace {
    SumEE (EdgeAB *a, EdgeAB *b) {
      ee[0] = a;
      ee[1] = b;
      edges.push_back(SumEV::make(a, b->getTail()));
      edges.push_back(SumVE::make(a->getHead(), b));
      edges.push_back(SumEV::make(a->getTwin(), b->getHead()));
      edges.push_back(SumVE::make(a->getTail(), b->getTwin()));
      // for (int i = 0; i < 4; i++)
      // edges[i]->face = this;
      a->getC4D()->sumEEs.push_back(this);
      get(a, b) = this;
    }
  public:
    EdgeAB *ee[2];

    void getContact (int i, Contact &contact) {
      ee[0]->getContact(i, contact);
      ee[1]->getContact(i, contact);
    }

    Type getType () { return SEE; }
    bool lessThan (Obj *o) {
      SumEE *that = dynamic_cast<SumEE*>(o);
      // if (ee[0] != that->ee[0])
      if (ee[0]->lessThan(that->ee[0]) || that->ee[0]->lessThan(ee[0]))
        return ee[0]->lessThan(that->ee[0]);
      return ee[1]->lessThan(that->ee[1]);
    }
    size_t hashCode () const {
      return ee[0]->hashCode() +  1957747793 * ee[1]->hashCode();
    }
      
    /*
    static int getID (EdgeAB *a, EdgeAB *b) {
      return a->id + b->id * a->getC4D()->edges[0].size();
    }
    int getID () { return getID(ee[0], ee[1]); }
    */

    static SumEE *&get (EdgeAB *a, EdgeAB *b) {
      // return a->getC4D()->sumEEs[getID(a, b)];
      return a->getC4D()->i2sumEE[Int2(a->id, b->id)];
    }

    static SumEE *make (EdgeAB *a, EdgeAB *b) {
      SumEE *ab = get(a, b);
      if (ab != 0)
        return ab;
      return new SumEE(a, b);
    }
    
    Declare(calculate) {
      return ee[0]->get<N>(event).cross(ee[1]->get<N>(event));
    }

    Declare(calculated) {
      return (ee[0]->getd<N>(event).cross(ee[1]->get<N>(event)) +
              ee[0]->get<N>(event).cross(ee[1]->getd<N>(event)));
    }

    bool generators (VertXYZ *fixed[3], VertXYZ *moving[3]) {
      return (addGenerator(ee[0]->getTail(), fixed) &&
              addGenerator(ee[0]->getTwin()->getTail(), fixed) &&
              addGenerator(ee[1]->getTail(), moving) &&
              addGenerator(ee[1]->getTwin()->getTail(), moving));
    }
  };

  class PolyVVVV : public AnglePoly {
  public:
    virtual C4D *getC4D () = 0;

    virtual int sign () = 0;

    EventSet getNegative () {
      vector<PTR<Angle>> angles = getAnglesSorted();
      // assert(angles.size() == 2);
      EventSet events(getC4D(), angles);
      if (sign() > 0)
        return events;
      else
        return events.complement();
    }
  };

  template<class N>
  static void setABC (Poly2<N> &poly, const PV3<N> &a, const PV3<N> &b, const N &c) {
    poly.a[0] = a.z * b.z + c;
    poly.a[1] = a.x * b.x + a.y * b.y;
    poly.a[2] = a.y * b.x - a.x * b.y;
  }
  
  template<class N>
  static void setAB (Poly2<N> &poly, const PV3<N> &a, const PV3<N> &b) {
    poly.a[0] = a.z * b.z;
    poly.a[1] = a.x * b.x + a.y * b.y;
    poly.a[2] = a.y * b.x - a.x * b.y;
  }
  
  template<class N>
  static void addAB (Poly2<N> &poly, const PV3<N> &a, const PV3<N> &b) {
    poly.a[0] = poly.a[0] + a.z * b.z;
    poly.a[1] = poly.a[1] + a.x * b.x + a.y * b.y;
    poly.a[2] = poly.a[2] + a.y * b.x - a.x * b.y;
  }

  template<class N>
  static Poly2<N> getPoly (PV3<N> a, PV3<N> b);
  template<class N>
  static Poly2<N> getPoly (PV3<N> a, PV3<N> b, SumFace *f);
  template<class N>
  static Poly2<N> getPolyA (PV3<N> a, SumFace *f);
  template<class N>
  static Poly2<N> getPolyB (PV3<N> b, SumFace *f);
  template<class N>
  static Poly2<N> getPoly (SumEdge *e, SumFace *f);
  template<class N>
  static Poly2<N> getPoly (SumEdge *e0, SumEdge *e1, SumEdge *e2);
  template<class N>
  static Poly2<N> getPoly (SumFace *f0, SumFace *f1, SumFace *f2);

  class PolyVF : public PolyVVVV {
    SumVV *v;
    SumFace *f;
    DeclareCalculate(Poly2);
  public:
    Type getType () { return VF; }
    void getContact (Contact &contact) {
      v->getContact(0, contact);
      f->getContact(1, contact);
    }

    PolyVF (SumVV *v, SumFace *f) : v(v), f(f) {}
    virtual C4D *getC4D () { return v->getC4D(); }

    Primitive2(Sign, SumVV*, v, SumFace*, f);
    int sign () { return Sign(v, f); }
    bool identical (AnglePoly *that);

    static EventSet getNegative (SumVV *v, SumFace *f);
    static bool sameEvent (Event *event, SumVV *v, SumFace *f);
  };

  class PolyEE : public PolyVVVV {
    SumEdge *e0, *e1;
    DeclareCalculate(Poly2);
  public:
    Type getType () { return EE; }
    void getContact (Contact &contact) {
      e0->getContact(0, contact);
      e1->getContact(1, contact);
    }

    virtual C4D *getC4D () { return e0->getC4D(); }
    PolyEE (SumEdge *e0, SumEdge *e1) {
      if (e0->getTwin()->lessThan(e0))
        e0 = e0->getTwin();
      if (e1->getTwin()->lessThan(e1))
        e1 = e1->getTwin();
      if (Compare()(e0, e1)) {
        this->e0 = e0;
        this->e1 = e1;
      }
      else {
        this->e0 = e1;
        this->e1 = e0;
      }
    }

    Primitive2(Sign, SumEdge*, e0, SumEdge*, e1);
    int sign () { return Sign(e0, e1); }
    bool identical (AnglePoly *thatAnglePoly);

    static EventSet getNegative (SumEdge *e0, SumEdge *e1);
    static bool sameEvent (Event *event, SumEdge *e0, SumEdge *e1);
  };

  class IntEF : public Vert {
    IntEF (SumEdge *e, SumFace *f) : e(e), f(f) {}
  public:
    SumEdge *e;
    SumFace *f;
    Type getType () { return IEF; }
    bool lessThan (Obj *o) {
      IntEF *that = dynamic_cast<IntEF*>(o);
      // if (f != that->f)
      if (Compare()(f, that->f) || Compare()(that->f, f))
        return Compare()(f, that->f);
      return Compare()(e, that->e);
    }
    size_t hashCode () const {
      return e->hashCode() + 424238335 * f->hashCode();
    }

    static IntEF *get (SumEdge *e, SumFace *f) {
      IntEF dummy(e, f);
      unordered_set<IntEF *, Equal, Equal> &intEFs = e->getC4D()->intEFs;
      unordered_set<IntEF *, Equal, Equal>::iterator it = intEFs.find(&dummy);
      if (it != intEFs.end())
        return *it;
      return 0;
    }      

    static IntEF *get (SumEdge *e, SumFace *f,
                       unordered_set<Obj*, Equal, Equal> &objSet) {
      IntEF dummy(e, f);
      unordered_set<Obj*, Equal, Equal>::iterator it = objSet.find(&dummy);
      if (it != objSet.end())
        return dynamic_cast<IntEF*>(*it);
      return 0;
    }      

    static IntEF *make (SumEdge *e, SumFace *f) {
      IntEF dummy(e, f);
      unordered_set<IntEF *, Equal, Equal> &intEFs = e->getC4D()->intEFs;
      unordered_set<IntEF *, Equal, Equal>::iterator it = intEFs.find(&dummy);
      IntEF *v = 0;
      if (it != intEFs.end())
        v = *it;
      else {
        v = new IntEF(e, f);
        intEFs.insert(v);
      }
      return v;
    }      

    Declare(calculate) { return f->intersection<N>(event, e); }
    Declare(calculated) { return f->intersectiond<N>(event, e); }
    C4D *getC4D () { return e->getC4D(); }

    static EventSet getLife (SumEdge *e, SumFace *f);

    // IntEFs that have been next to it on the same SumEdge.
    // Stored in lesser IntEF.
    unordered_set<IntEF *, Equal, Equal> nexts;
    bool solved (IntEF *that) {
      if (that->lessThan(this))
        return that->solved(this);
      if (nexts.find(that) != nexts.end())
        return true;
      nexts.insert(that);
      return false;
    }

    set<IntEdge*, Compare> eSet;
  };
  unordered_set<IntEF*, Equal, Equal> intEFs;

  IntEF *findLiveEF (SumEdge *e, SumFace *f);

  class IntEdge : public Edge {
    IntEdge (IntEdge *twin) : Edge(0, twin) {}

    IntEdge (SumFace *f, SumFace *g) : Edge(0, new IntEdge(this)) {
      face = f;
      twin->face = g;
      f->getC4D()->intEdges.insert(this);
      g->getC4D()->intEdges.insert(getTwin());
    }

  public:
    IntEdge (SumFace *f) : Edge(0, 0) { face = f; }

    static IntEdge *get (SumFace *f, SumFace *g) {
      IntEdge dummyF(f);
      IntEdge dummyG(g);
      dummyF.twin = &dummyG;
      dummyG.twin = &dummyF;
      unordered_set<IntEdge*, Equal, Equal> &intEdges = f->getC4D()->intEdges;
      // unordered_set<IntEdge*, Equal, Equal> &intEdgesG = g->getC4D()->intEdges;
      unordered_set<IntEdge*, Equal, Equal>::iterator it = intEdges.find(&dummyF);
      if (it != intEdges.end())
        return *it;
      return 0;
    }      

    static IntEdge *make (SumFace *f, SumFace *g) {
      IntEdge dummyF(f);
      IntEdge dummyG(g);
      dummyF.twin = &dummyG;
      dummyG.twin = &dummyF;
      unordered_set<IntEdge*, Equal, Equal> &intEdges = f->getC4D()->intEdges;
      unordered_set<IntEdge*, Equal, Equal>::iterator it = intEdges.find(&dummyF);
      IntEdge *e = 0;
      if (it != intEdges.end()) {
        e = *it;
        e->tail = 0;
        e->twin->tail = 0;
      }
      else {
        e = new IntEdge(f, g);
	if (e == (IntEdge *) 0x2990ce0)
	  cout << "this is it IntEdge" << endl;
      }
      return e;
    }      

    Type getType () { return IFF; }
    bool lessThan (Obj *o) {
      IntEdge *that = dynamic_cast<IntEdge*>(o);
      // if (face != that->face || twin == 0 || that->twin == 0)
      if (Compare()(face, that->face) || Compare()(that->face, face) ||
          twin == 0 || that->twin == 0)
        return Compare()(face, that->face);
      else
        return Compare()(twin->face, that->twin->face);
    }
    size_t hashCode () const {
      assert(twin != 0);
      // return face->hashCode() + twin == 0 ? 0 : 719885386 * twin->face->hashCode();
      return face->hashCode() + 719885386 * twin->face->hashCode();
    }

    IntEdge *getTwin () { return dynamic_cast<IntEdge*>(twin); }
    SumFace *getFace () { return dynamic_cast<SumFace*>(face); }

    Declare(calculate) {
      return twin->face->get<N>(event).cross(face->get<N>(event));
    }

    Declare(calculated) {
      return (twin->face->getd<N>(event).cross(face->get<N>(event)) +
              twin->face->get<N>(event).cross(face->getd<N>(event)));
    }

    void checkForIntersections (Event *event);

    void setTail (Vert *v) {
      assert(tail != v);
      if (this == (IntEdge *) 0x9e9460 && v == (Vert *) 0x6c5cd0)
	cout << "setTail this is it" << endl;
      if (tail != 0 && tail->getType() == SVV)
        dynamic_cast<SumVV*>(tail)->removeIntEdge(this);
      tail = v;
      if (tail != 0 && tail->getType() == SVV)
        dynamic_cast<SumVV*>(tail)->addIntEdge(this);
    }
  };    
  unordered_set<IntEdge*, Equal, Equal> intEdges;

  class SubEdge : public Edge {
    SubEdge (Vert *v, Edge *e) : Edge(v, 0), e(e) {}

    SubEdge (Vert *t, Vert *h, Edge *e)
      : Edge(t, new SubEdge(h, e->twin, this)), e(e) {}

    SubEdge (Vert *v, Edge *e, SubEdge *twin) : Edge(v, twin), e(e) {}

  public:
    Edge *e;

    static SubEdge *make (Vert *t, Vert *h, Edge *e) {
      SubEdge dummy0(t, e);
      SubEdge dummy1(h, e->twin);
      dummy0.twin = &dummy1;
      dummy1.twin = &dummy0;
      unordered_set<SubEdge *, Equal, Equal> &subEdges = t->getC4D()->subEdges;
      unordered_set<SubEdge *, Equal, Equal>::iterator it = subEdges.find(&dummy0);
      if (it != subEdges.end()) {
	if (*it == (C4D::SubEdge *) 0x29c81a0)
	  cout << "this is it SubEdge1" << endl;
        return *it;
      }
      SubEdge *sub = new SubEdge(t, h, e);
      subEdges.insert(sub);
      subEdges.insert(sub->getTwin());
      if (sub == (C4D::Obj *) 0x2b265e50 || sub == (C4D::Edge *) 0x2b265ec0)
	cout << "this is it SubEdge2" << endl;
      if (e == (Edge *) 0x28ccf280 &&
	  t == (Vert *) 0x29256530 &&
	  h == (Vert *) 0x2927cb80) {
	extern  Event *global_event;
	int uve = OrderUVE(global_event, t, h, e);
	cout << "this is it SubEdge2" << endl;
      }
      return sub;
    }

    Type getType () { return SUBE; }
    bool lessThan (Obj *o) {
      SubEdge *that = dynamic_cast<SubEdge*>(o);
      // if (e != that->e)
      if (Compare()(e, that->e) || Compare()(that->e, e))
        return Compare()(e, that->e);
      // if (tail != that->tail)
      if (Compare()(tail, that->tail) || Compare()(that->tail, tail))
        return Compare()(tail, that->tail);
      return Compare()(twin->tail, that->twin->tail);
    }
    size_t hashCode () const {
      return e->hashCode() + 1649760492 * tail->hashCode() + 596516649 * twin->tail->hashCode();
    }

    Declare(calculate) { return e->get<N>(event); }
    Declare(calculated) { return e->getd<N>(event); }

    SubEdge *getTwin () { return dynamic_cast<SubEdge*>(twin); }

    void addEvents();
  };      
  unordered_set<SubEdge *, Equal, Equal> subEdges;

  class BigPoly : public AnglePoly {
  public:
    BigPoly () {}

    virtual bool isEpEEE (EdgeAB *&ab, EdgeAB *cd[3]) = 0;

    class CompareEpEEE {
    public:
      bool operator() (BigPoly *poly1, BigPoly *poly2) const {
        EdgeAB *ab1, *cd1[3];
        bool poly1test = poly1->isEpEEE(ab1, cd1);
        assert(poly1test);
        EdgeAB *ab2, *cd2[3];
        bool poly2test = poly2->isEpEEE(ab2, cd2);
        assert(poly2test);
        if (ab1 != ab2)
          return ab1->lessThan(ab2);
        for (int i = 0; i < 3; i++)
          if (cd1[i] != cd2[i])
            return cd1[i]->lessThan(cd2[i]);
        return false;
      }
    };

    virtual bool isEEEpE (EdgeAB *cd[3], EdgeAB *&ab) = 0;

    class CompareEEEpE {
    public:
      bool operator() (BigPoly *poly1, BigPoly *poly2) const {
        EdgeAB *ab1, *cd1[3];
        bool poly1test = poly1->isEEEpE(cd1, ab1);
        assert(poly1test);
        EdgeAB *ab2, *cd2[3];
        bool poly2test = poly2->isEEEpE(cd2, ab2);
        assert(poly2test);
        if (ab1 != ab2)
          return ab1->lessThan(ab2);
        for (int i = 0; i < 3; i++)
          if (cd1[i] != cd2[i])
            return cd1[i]->lessThan(cd2[i]);
        return false;
      }
    };
  };

  set<BigPoly*, BigPoly::CompareEpEEE> polyEpEEEs;
  set<BigPoly*, BigPoly::CompareEEEpE> polyEEEpEs;

  class PolyEFF : public BigPoly {
    DeclareCalculate(Poly2);
  public:
    SumEdge *e;
    SumFace *f1, *f2;

    Type getType () { return EFF; }
    void getContact (Contact &contact) {
      e->getContact(0, contact);
      f1->getContact(1, contact);
      f2->getContact(2, contact);
    }

    PolyEFF (SumEdge *e, SumFace *f1, SumFace *f2) {
      if (e->lessThan(e->getTwin()))
        this->e = e;
      else
        this->e = e->getTwin();
      if (C4D::Compare()(f1, f2)) {
        this->f1 = f1;
        this->f2 = f2;
      }
      else {
        this->f1 = f2;
        this->f2 = f1;
      }
    }

    class Compare {
    public:
      bool operator() (PolyEFF *a, PolyEFF *b) const {
        if (a->e != b->e)
          return C4D::Compare()(a->e, b->e);
        if (a->f1 != b->f1)
          return C4D::Compare()(a->f1, b->f1);
        if (a->f2 != b->f2)
          return C4D::Compare()(a->f2, b->f2);
        return false;
      }
    };

    size_t hashCode () const {
      return e->hashCode() + 783368690 * f1->hashCode() + 1102520059 * f2->hashCode();
    }

    class Equal {
    public:
      bool operator() (PolyEFF *a, PolyEFF *b) const {
        return a->e == b->e && a->f1 == b->f1 && a->f2 == b->f2;
      }

      size_t operator() (PolyEFF *a) const {
        return a->hashCode();
      }
    };

    bool isEpEEE (EdgeAB *&ab, EdgeAB *cd[3]);
    bool isEEEpE (EdgeAB *cd[3], EdgeAB *&ab);
    bool isEpEEE_ (EdgeAB *&ab, EdgeAB *cd[3]);
    bool isEEEpE_ (EdgeAB *cd[3], EdgeAB *&ab);
  };

  // set<PolyEFF*, PolyEFF::Compare> polyEFFs;
  unordered_set<PTR<PolyEFF>, PolyEFF::Equal, PolyEFF::Equal> polyEFFs;
    
  void addEvents (SumEdge *e, SumFace *f1, SumFace *f2);

#ifdef BLEEN
  class PolyEEE : public AnglePoly {
    SumEdge *e[3];
  public:
    Type getType () { return EEE; }
    PolyEEE (SumEdge *e0, SumEdge *e1, SumEdge *e2);
    Poly2 getPoly () { return C4D::getPoly(e[0], e[1], e[2]); }

    class Compare {
    public:
      bool operator() (PolyEEE *a, PolyEEE *b) const {
        for (int i = 0; i < 3; i++)
          if (a->e[i] != b->e[i])
            return C4D::Compare()(a->e[i], b->e[i]);
        return false;
      }
    };
  };

  class PolyFFF : public AnglePoly {
    SumFace *f[3];
  public:
    Type getType () { return FFF; }
    PolyFFF (SumFace *f0, SumFace *f1, SumFace *f2);
    Poly2 getPoly () { return C4D::getPoly(f[0], f[1], f[2]); }

    class Compare {
    public:
      bool operator() (PolyFFF *a, PolyFFF *b) const {
        for (int i = 0; i < 3; i++)
          if (a->f[i] != b->f[i])
            return C4D::Compare()(a->f[i], b->f[i]);
        return false;
      }
    };
  };
#endif

  class PolyEEE : public AnglePoly {
    //SumEdge *e[3];
    VertXYZ *e[3][2];
    DeclareCalculate(Poly2);
  public:
    Type getType () { return EEE; }

    void getContact (Contact &contact) {
      for(int i = 0; i < 3; i++)
	for (int j = 0; j < 2; j++)
	  e[i][j]->getContact(i, contact);
    }

    //PolyEEE (SumEdge *e0, SumEdge *e1, SumEdge *e2)
    PolyEEE (VertXYZ *e[3][2]) {
      for (int i = 0; i < 3; i++)
	for (int j = 0; j < 2; j++)
	  this->e[i][j] = e[i][j];
    }
      
    template<class N>
    PV3<N> getV (int i) {
      return e[i][1]->get<N>() - e[i][0]->get<N>();
    }

    // Poly2 getPoly () { return C4D::getPoly(e[0], e[1], e[2]); }

    class Compare {
    public:
      bool operator() (PolyEEE *a, PolyEEE *b) const {
	assert(0);
	/*
        for (int i = 0; i < 3; i++)
          if (a->e[i] != b->e[i])
            return C4D::Compare()(a->e[i], b->e[i]);
	*/
        return false;
      }
    };
  };

  class PolyFFF : public AnglePoly {
    VertXYZ *v[3][4];
    // SumFace *f[3];
    DeclareCalculate(Poly2);
  public:
    Type getType () { return FFF; }

    void getContact (Contact &contact) {
      for(int i = 0; i < 3; i++)
	for (int j = 0; j < 4; j++)
	  if (v[i][j] != 0)
	    v[i][j]->getContact(i, contact);
    }

    // PolyFFF (SumFace *f0, SumFace *f1, SumFace *f2);
    PolyFFF (VertXYZ *v[3][4]) {
      for (int i = 0; i < 3; i++)
	for (int j = 0; j < 4; j++)
	  this->v[i][j] = v[i][j];
    }

    int nMoving (int i) {
      if (v[i][3] != 0)
	return 2;
      if (v[i][0]->isMoving())
	return 3;
      return 0;
    }	

    template<class N>
    PV3<N> getN (int i) {
      return (v[i][1]->get<N>() - v[i][0]->get<N>()).cross(v[i][2]->get<N>() - v[i][0]->get<N>());
    }

    template<class N>
    PV3<N> getV (int i, int fm) {
      return v[i][2*fm + 1]->get<N>() - v[i][2*fm]->get<N>();
    }

    // Poly2 getPoly () { return C4D::getPoly(f[0], f[1], f[2]); }

    template<class N>
    Poly2<N> getPolyA (PV3<N> a, int i);
    template<class N>
    Poly2<N> getPolyB (PV3<N> b, int i);

    class Compare {
    public:
      bool operator() (PolyFFF *a, PolyFFF *b) const {
	assert(0);
#ifdef BLEEN
        for (int i = 0; i < 3; i++)
          if (a->f[i] != b->f[i])
            return C4D::Compare()(a->f[i], b->f[i]);
#endif
        return false;
      }
    };
  };

  class PolyNNN : public AnglePoly {
    void getContact (Contact &contact) { assert(0); }

    SumVV *vert;

    class Pair {
    public:
      VertXYZ *vv[2];

      void setVV (VertXYZ *a, VertXYZ *b) {
        if (a->isMoving() < b->isMoving() ||
            (a->isMoving() ==  b->isMoving() && a->id < b->id)) {
          vv[0] = a;
          vv[1] = b;
        }
        else {
          vv[0] = b;
          vv[1] = a;
        }
      }

      Pair () { vv[0] = vv[1] = 0; }
      Pair (VertXYZ *a, VertXYZ *b) { setVV(a, b); }
      Pair (SumVV *v, SumFace *f);
      Pair (SumVV *v, SumEdge *e);

      bool operator== (const Pair &that) const {
        return vv[0] == that.vv[0] && vv[1] == that.vv[1];
      }

      bool operator!= (const Pair &that) const {
        return vv[0] != that.vv[0] || vv[1] != that.vv[1];
      }

      bool operator< (const Pair &that) const {
        if (vv[0]->isMoving() < that.vv[0]->isMoving())
          return true;
        if (vv[0]->isMoving() > that.vv[0]->isMoving())
          return false;
        if (vv[1]->isMoving() < that.vv[1]->isMoving())
          return true;
        if (vv[1]->isMoving() > that.vv[1]->isMoving())
          return false;
        if (vv[0]->id < that.vv[0]->id)
          return true;
        if (vv[0]->id > that.vv[0]->id)
          return false;
        if (vv[1]->id < that.vv[1]->id)
          return true;
        if (vv[1]->id > that.vv[1]->id)
          return false;
        return false;
      }
    } pairs[3];

    template<class N>
    Poly2<N> getPolyA (PV3<N> a, Pair &pair);
    template<class N>
    Poly2<N> getPolyB (PV3<N> b, Pair &pair);

    DeclareCalculate(Poly2);

  public:
    Type getType () { return NNN; }
    PolyNNN (SumFace *f0, SumFace *f1, SumFace *f2);
    PolyNNN (SumFace *f0, SumFace *f1, SumEdge *e);
    virtual bool identical (AnglePoly *that);
    class Compare {
    public:
      bool operator() (PolyNNN *a, PolyNNN *b) const {
        if (a->vert != b->vert)
          return a->vert->lessThan(b->vert);
        for (int i = 0; i < 3; i++)
          if (a->pairs[i] != b->pairs[i])
            return a->pairs[i]< b->pairs[i];
        return false;
      }
    };

    void getFaces (Pair &pair, SumFace *faces[2]);
    void getEdges (Pair &pair, SumEdge *edges[2]);

    void addEvents ();
    void handleEvent (Event *event);
  };

  set<PTR<PolyNNN>, PolyNNN::Compare> polyNNNs;

  class PolyFFFF : public BigPoly {
    DeclareCalculate(Poly2);
  public:
    SumFace *f[4];
    Type getType () { return FFFF; }
    void getContact (Contact &contact) {
      for (int i = 0; i < 4; i++)
        f[i]->getContact(i, contact);
    }

    PolyFFFF (SumFace *f0, SumFace *f1, SumFace *f2, SumFace *f3);

    class Compare {
    public:
      bool operator() (PolyFFFF *a, PolyFFFF *b) const {
        for (int i = 0; i < 4; i++)
          if (a->f[i] != b->f[i])
            return C4D::Compare()(a->f[i], b->f[i]);
        return false;
      }
    };

    size_t hashCode () const {
      return f[0]->hashCode() +
        2044897763 * f[1]->hashCode() +
        1967513926 * f[2]->hashCode() +
        1365180540 * f[3]->hashCode();
    }
    class Equal {
    public:
      bool operator() (PolyFFFF *a, PolyFFFF *b) const {
        for (int i = 0; i < 4; i++)
          if (a->f[i] != b->f[i])
            return false;
        return true;
      }

      size_t operator() (PolyFFFF *a) const { return a->hashCode(); }
    };

    bool isEpEEE (EdgeAB *&ab, EdgeAB *cd[3]);
    bool isEEEpE (EdgeAB *cd[3], EdgeAB *&ab);
    bool isEpEEE_ (EdgeAB *&ab, EdgeAB *cd[3]);
    bool isEEEpE_ (EdgeAB *cd[3], EdgeAB *&ab);
  };

  // set<PolyFFFF*, PolyFFFF::Compare> polyFFFFs;
  unordered_set<PTR<PolyFFFF>, PolyFFFF::Equal, PolyFFFF::Equal> polyFFFFs;
    
  void addEvents (SumFace *f0, SumFace *f1, SumFace *f2, SumFace *f3);

  void handleFFFF (Event *event);

  class IntFFF : public Vert {
    // Faces are counterclockwise as viewed from octant that is
    // outside all faces.

    IntFFF (SumFace *f, SumFace *g, SumFace *h) {
      SumFace *fgh[] = { f, g, h };
      int i = 0;
      for (int j = 1; j < 3; j++)
        if (Compare()(fgh[j], fgh[i]))
          i = j;
      for (int j = 0; j < 3; j++)
        faces[j] = fgh[(i+j)%3];
    }

  public:
    SumFace *faces[3];
    Type getType () { return IFFF; }
    static IntFFF *make (SumFace *f, SumFace *g, SumFace *h) {
      IntFFF dummy(f, g, h);
      unordered_set<IntFFF *, Equal, Equal> &intFFFs = f->getC4D()->intFFFs;
      unordered_set<IntFFF *, Equal, Equal>::iterator it = intFFFs.find(&dummy);
      if (it != intFFFs.end()) {
        return *it;
      }
      IntFFF *v = new IntFFF(f, g, h);
      intFFFs.insert(v);
      return v;
    }
    static IntFFF *get (SumFace *f, SumFace *g, SumFace *h) {
      IntFFF dummy(f, g, h);
      unordered_set<IntFFF *, Equal, Equal> &intFFFs = f->getC4D()->intFFFs;
      unordered_set<IntFFF *, Equal, Equal>::iterator it = intFFFs.find(&dummy);
      if (it != intFFFs.end()) {
        return *it;
      }
      return 0;
    }
    bool lessThan (Obj *o) {
      IntFFF *that = dynamic_cast<IntFFF*>(o);
      for (int i = 0; i < 3; i++)
        // if (faces[i] != that->faces[i])
        if (Compare()(faces[i], that->faces[i]) ||
            Compare()(that->faces[i], faces[i]))
          return Compare()(faces[i], that->faces[i]);
      return false;
    }
    size_t hashCode () const {
      return faces[0]->hashCode() +
        1189641421 * faces[1]->hashCode() +
        1025202362 * faces[2]->hashCode();
    }

    C4D *getC4D () { return faces[0]->getC4D(); }

    Declare(calculate) {
      // return faces[2]->intersection(event, faces[0]->findLiveEdge(faces[1]));
      for (int i = 0; i < 3; i++) {
	IntEdge *e = IntEdge::get(faces[(i+1)%3], faces[(i+2)%3]);
	if (e != 0)
	  return faces[i]->intersection<N>(event, e);
      }
      assert(0);
      return PV3<N>(N(0), N(0), N(0));
    }

    Declare(calculated) {
      // return faces[2]->intersectiond(event, faces[0]->findLiveEdge(faces[1]));
      for (int i = 0; i < 3; i++) {
	IntEdge *e = IntEdge::get(faces[(i+1)%3], faces[(i+2)%3]);
	if (e != 0)
	  return faces[i]->intersectiond<N>(event, e);
      }
      assert(0);
      return PV3<N>(N(0), N(0), N(0));
    }

    SumFace *otherFace (SumFace *f0, SumFace *f1) {
      for (int i = 0; i < 3; i++)
        if (faces[i] != f0 && faces[i] != f1)
          return faces[i];
      assert(0);
      return 0;
    }

    bool hasFace (SumFace *f) {
      return faces[0] == f || faces[1] == f || faces[2] == f;
    }

    IntEdge *commonEdge (IntFFF *that) {
      for (int i = 0; i < 3; i++)
        if (!that->hasFace(faces[i])) {
          assert(that->hasFace(faces[(i+1)%3]));
          assert(that->hasFace(faces[(i+2)%3]));
          return IntEdge::get(faces[(i+1)%3], faces[(i+2)%3]);
        }
      return 0;
    }

    // IntFFFs that have been next to it on the same IntEdge.
    // Stored in lesser IntFFF.
    unordered_set<IntFFF *, Equal, Equal> nexts;
    bool solved (IntFFF *that) {
      if (that->lessThan(this))
        return that->solved(this);
      if (nexts.find(that) != nexts.end())
        return true;
      nexts.insert(that);
      return false;
    }
  };
  unordered_set<IntFFF *, Equal, Equal> intFFFs;

  bool readObjFile (const char *file, bool moving);
  bool readObjFile (istream &in, bool moving);

  void initSums (bool noCone=false);

  vector<Face *> getSumFaces ();
  vector<Face *> getSumFaces (Event *event);
  set<Obj *, Compare> getObjects ();
  set<Obj *, Compare> getObjects$ ();
  set<Obj *, Compare> getObjects (Event *event);
  set<Obj *, Compare> getObjects$ (Event *event) {
    assert(0);
    return getObjects(event);
  }

  set<Obj *, Compare> getSums (vector<Face *> &sumFaces);
  void  getInts (vector<Face *> &sumFaces, SumFace *f, 
                 set<Obj *, Compare> &objs);

  set<Event*, Compare> eventList;
  // set<Obj*, Compare> sweepSet;
  unordered_set<Obj*, Equal, Equal> sweepSets[16];

  void sweep ();

  Event *currentEvent;
  void sweep2 (int iFace);
  void sweep2 (vector<int> &iFaces, int from, int to);
  void sweep2 (int iFace, vector<Face *> & sumFaces, set<Obj*, Compare> &objs);

  void intersectEdges (SumFace *base, Event *event, IntEdge *d, IntEdge *e, 
                       int &nFFFs, bool insert);

  C4D () : event0(this, angle0), EVENT0(&event0) {}

  C4D (C4D *that);

  Obj *findObj (Obj *obj, set<C4D::Obj *, C4D::Compare> objs);
  bool compare (Event *event, set<C4D::Obj *, C4D::Compare> &objs0,
                set<C4D::Obj *, C4D::Compare> &objs1);

  Event *current;

  EdgeAB *makeEdge (VertXYZ *a, VertXYZ *b);
  FaceABC *makeFace (VertXYZ *a, VertXYZ *b, VertXYZ *c);
  void findIdentities (int nFixed, int nMoving);
  void contacts2polys ();
  void findIdentitiesModP ();
  void findIdentities ();
  void findIdentitiesT ();

  AnglePoly *getPoly (Face4 face4r);
  AnglePoly *getPoly (Contact &contact) {
    Face4 face4(contact.ns[0], contact.is[0],
		contact.ns[1], contact.is[1]);
    return getPoly(face4);
  }

  Contact getContact (AnglePoly *poly) {
    Contact contact;
    poly->getContact(contact);
    contact.normalize();
    return contact;
  }

  Primitive2(ContactSign, AnglePoly*, poly, RootAngle*, root);
  bool commonRoot (Contact &contact, RootAngle *root);

  unordered_map<Contact, vector<Event*>, Contact::Equal, Contact::Equal> contact2events;

  vector<Event*> solve (Contact &contact);

  Contact getContact (IntEF *a, IntEF *b) {
    assert(a->e == b->e || a->e == b->e->getTwin());
    Contact contact;
    a->e->getContact(0, contact);
    a->f->getContact(1, contact);
    b->f->getContact(2, contact);
    contact.normalize();
    return contact;
  }

  Contact getContact (SumVV *v, SumFace *f) {
    Contact contact;
    v->getContact(0, contact);
    f->getContact(1, contact);
    contact.normalize();
    return contact;
  }

  Contact getContact (IntEF *v, SumFace *f) {
    Contact contact;
    v->e->getContact(0, contact);
    v->f->getContact(1, contact);
    f->getContact(2, contact);
    contact.normalize();
    return contact;
  }

  Contact getContact (Vert *v, SumFace *f) {
    Type vType = v->getType();
    if (vType == SVV)
      return getContact(dynamic_cast<SumVV*>(v), f);
    else if (vType == IEF)
      return getContact(dynamic_cast<IntEF*>(v), f);
    else
      assert(0);
    return Contact();
  }

  Contact getContact (IntFFF *a, IntFFF *b) {
    Contact contact;
    SumFace *faces[6] = { a->faces[0], a->faces[1], a->faces[2] };
    int n = 3;
    for (int i = 0; i < 3; i++) {
      bool seen = false;
      for (int j = 0; j < n; j++)
        if (b->faces[i] == faces[j])
          seen = true;
      if (!seen)
        faces[n++] = b->faces[i];
    }
    assert(n == 4);
    for (int i = 0; i < 4; i++)
      faces[i]->getContact(i, contact);
    contact.normalize();
    return contact;
  }

  void generateSwaps (IntEF *a, IntEF *b) {
    if (a->solved(b))
      return;
    Contact contact = getContact(a, b);
    vector<Event*> events = solve(contact);
    for (int i = 0; i < events.size(); i++) {
      events[i]->swap.insert(pair<Obj*,Obj*>(a,b));
    }
  }

  void generateSwaps (SumEdge *e, int i) {
    if (i < 0 || i+1 >= e->verts.size())
      return;
    generateSwaps(dynamic_cast<IntEF*>(e->verts[i]),
                  dynamic_cast<IntEF*>(e->verts[i+1]));
  }

  void generateSwaps (IntEF *v, SumFace *f) {
    if (f->solved(v))
      return;
    Contact contact = getContact(v, f);
    vector<Event*> events = solve(contact);
    for (int i = 0; i < events.size(); i++) {
      events[i]->swap.insert(pair<Obj*,Obj*>(v,f));
    }
  }

  void generateSwaps (IntFFF *a, IntFFF *b) {
    if (a->solved(b))
      return;
    Contact contact = getContact(a, b);
    vector<Event*> events = solve(contact);
    for (int i = 0; i < events.size(); i++) {
      if ((a->life.events.size() > 0 && a->life.events.back() == events[i]) ||
	  (b->life.events.size() > 0 && b->life.events.back() == events[i]))
	continue;
      events[i]->swap.insert(pair<Obj*,Obj*>(a,b));
    }
  }

  void generateSwaps (IntEdge *e, int i) {
    if (i < 0 || i+1 >= e->verts.size())
      return;
    generateSwaps(dynamic_cast<IntFFF*>(e->verts[i]),
                  dynamic_cast<IntFFF*>(e->verts[i+1]));
  }

  void checkEdgeList (Event *event, SumFace *base, IntEdge *e,
                      unordered_set<IntEdge*, Equal, Equal> &eSet);

  void checkFactors ();
};

template<class N>
inline N C4D::EdgeAB::IsOutside::calculate () {
  return f->get<N>(event).cross(e->get<N>(event)).dot(e->getFace()->get<N>(event));
}

/*
inline bool C4D::Equal::operator() (IntEdge *const a, IntEdge *const b) const {
  return a->face != b->face ? a->face < b->face : a->twin->face < b->twin->face;
}
*/

extern bool noSweep;
#endif
