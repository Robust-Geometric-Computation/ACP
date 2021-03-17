//#include "acp/encasement3d/encasement3d.h"
#include <fstream>
#include <iostream>
#include <map>
using namespace std;

class V2distV2V2mS : public Primitive {
  PTR<Object<PV2>> p, a, b;
  double epsilon;
  DeclareSign {
    PV2<N> ab = b->get<N>() - a->get<N>();
    PV2<N> ap = p->get<N>() - a->get<N>();

    // (ap - ab t) * ab = 0
    // ap * ab = ab * ab t
    // t = (ap * ab) / (ab * ab)
    // u = ap - ab t
    // u * u - eps * eps

    N t = ap.dot(ab) / ab.dot(ab);
    PV2<N> u = ap - ab * t;
    return u.dot(u) - epsilon * epsilon;
  }

 public:
  V2distV2V2mS(PTR<Object<PV2>> p, PTR<Object<PV2>> a, PTR<Object<PV2>> b,
               double epsilon)
      : p(p), a(a), b(b), epsilon(epsilon) {}
};

class V2V2ImJ : public Primitive {
  PTR<Object<PV2>> a, b;
  int i, j;
  DeclareSign {
    PV2<N> ab = b->get<N>() - a->get<N>();
    return ab[i] * ab[i] - ab[j] * ab[j];
  }

 public:
  V2V2ImJ(PTR<Object<PV2>> a, PTR<Object<PV2>> b, int i, int j)
      : a(a), b(b), i(i), j(j) {}
};

class MidPoint2 : public Object<PV2> {
  PTR<Object<PV2>> a, b;
  DeclareCalculate(PV2) { return (a->get<N>() + b->get<N>()) / 2; }

 public:
  MidPoint2(PTR<Object<PV2>> a, PTR<Object<PV2>> b) : a(a), b(b) {}
};

class V3V3mV3V3 : public Primitive {
  PTR<Object<PV3>> a, b, c, d;
  DeclareSign {
    PV3<N> ab = b->get<N>() - a->get<N>();
    PV3<N> cd = d->get<N>() - c->get<N>();
    return ab.dot(ab) - cd.dot(cd);
  }

 public:
  V3V3mV3V3(PTR<Object<PV3>> a, PTR<Object<PV3>> b, PTR<Object<PV3>> c,
            PTR<Object<PV3>> d)
      : a(a), b(b), c(c), d(d) {}
};

class V3distV3V3mS : public Primitive {
  PTR<Object<PV3>> p, a, b;
  double epsilon;
  DeclareSign {
    PV3<N> ab = b->get<N>() - a->get<N>();
    PV3<N> ap = p->get<N>() - a->get<N>();

    // (ap - ab t) * ab = 0
    // ap * ab = ab * ab t
    // t = (ap * ab) / (ab * ab)
    // u = ap - ab t
    // u * u - eps * eps

    N t = ap.dot(ab) / ab.dot(ab);
    PV3<N> u = ap - ab * t;
    return u.dot(u) - epsilon * epsilon;
  }

 public:
  V3distV3V3mS(PTR<Object<PV3>> p, PTR<Object<PV3>> a, PTR<Object<PV3>> b,
               double epsilon)
      : p(p), a(a), b(b), epsilon(epsilon) {}
};

class V3V3ImJ : public Primitive {
  PTR<Object<PV3>> a, b;
  int i, j;
  DeclareSign {
    PV3<N> ab = b->get<N>() - a->get<N>();
    return ab[i] * ab[i] - ab[j] * ab[j];
  }

 public:
  V3V3ImJ(PTR<Object<PV3>> a, PTR<Object<PV3>> b, int i, int j)
      : a(a), b(b), i(i), j(j) {}
};

class MidPoint3 : public Object<PV3> {
  PTR<Object<PV3>> a, b;
  DeclareCalculate(PV3) { return (a->get<N>() + b->get<N>()) / 2; }

 public:
  MidPoint3(PTR<Object<PV3>> a, PTR<Object<PV3>> b) : a(a), b(b) {}
};

class V3V3atV3I : public Object<PV3> {
  PTR<Object<PV3>> a, b, p;
  int i;
  DeclareCalculate(PV3) {
    PV3<N> ap = a->get<N>();
    PV3<N> ab = b->get<N>() - ap;
    return ap + ab * ((p->get<N>()[i] - ap[i]) / ab[i]);
  }

 public:
  V3V3atV3I(PTR<Object<PV3>> a, PTR<Object<PV3>> b, PTR<Object<PV3>> p, int i)
      : a(a), b(b), p(p), i(i) {}
};

class V3V3ofT : public Object<PV3> {
  PTR<Object<PV3>> a, b;
  double t;
  DeclareCalculate(PV3) {
    PV3<N> ap = a->get<N>();
    PV3<N> ab = b->get<N>() - ap;
    return ap + ab * t;
  }

 public:
  V3V3ofT(PTR<Object<PV3>> a, PTR<Object<PV3>> b, double t)
      : a(a), b(b), t(t) {}
};

class XV2 : public Object<PV3> {
  PTR<Object<Scalar>> x;
  PTR<Object<PV2>> v2;
  DeclareCalculate(PV3) {
    PV2<N> p2 = v2->get<N>();
    PV3<N> p;
    p.x = x->get<N>().x;
    p.y = p2.x;
    p.z = p2.y;
    return p;
  }

 public:
  XV2(PTR<Object<Scalar>> x, PTR<Object<PV2>> v2) : x(x), v2(v2) {}
};

class V3V3atSI : public Object<PV3> {
  PTR<Object<PV3>> a, b;
  PTR<Object<Scalar>> s;
  int i;
  DeclareCalculate(PV3) {
    PV3<N> a_ = a->get<N>();
    PV3<N> ab = b->get<N>() - a_;
    N x = s->get<N>().x;
    return a_ + ab * ((x - a_[i]) / ab[i]);
  }

 public:
  V3V3atSI(PTR<Object<PV3>> a, PTR<Object<PV3>> b, PTR<Object<Scalar>> s, int i)
      : a(a), b(b), s(s), i(i) {}
};

class P3yyzz_yzyz : public Primitive {
  PTR<Object<Poly3D>> poly3d;
  PTR<Object<PV3>> pv3;
  DeclareSign {
    Poly3D<N> f = poly3d->get<N>();
    PV3<N> p = pv3->get<N>();
    Poly3D<N> fy = f.der(1);
    Poly3D<N> fz = f.der(2);
    N fyy = fy.der(1, p);
    N fyz = fy.der(2, p);
    N fzz = fz.der(2, p);
    return fyy * fzz - fyz * fyz;
  }

 public:
  P3yyzz_yzyz(PTR<Object<Poly3D>> poly3d, PTR<Object<PV3>> pv3)
      : poly3d(poly3d), pv3(pv3) {}
};

class P3zzz : public Primitive {
  PTR<Object<Poly3D>> poly3d;
  PTR<Object<PV3>> pv3;
  DeclareSign {
    Poly3D<N> f = poly3d->get<N>();
    PV3<N> p = pv3->get<N>();
    Poly3D<N> fz = f.der(2);
    Poly3D<N> fzz = fz.der(2);
    N fzzz = fzz.der(2, p);
    return fzzz;
  }

 public:
  P3zzz(PTR<Object<Poly3D>> poly3d, PTR<Object<PV3>> pv3)
      : poly3d(poly3d), pv3(pv3) {}
};

class Mesh3D {
 public:
  Mesh3D(Encasement3D* enc3d, double epsilon)
      : enc3d(enc3d), epsilon(epsilon) {}

  Encasement3D* enc3d;
  double epsilon;

  // edges[e] is polygonal chain approximating edge e
  map<Encasement3D::Edge*, vector<PTR<Object<PV3>>>> edge2chain;

  class Triangle {
   public:
    PTR<Object<PV3>> p[3];

    Triangle(PTR<Object<PV3>> a, PTR<Object<PV3>> b, PTR<Object<PV3>> c) {
      p[0] = a;
      p[1] = b;
      p[2] = c;
    }
  };

  static bool isYZ(Encasement3D::Edge* edge) {
    return edge->b->jf == 0 || edge->b->jf == 1 ||
           edge->b->type == Encasement3D::FSB;
  }

  static bool isTail(Encasement3D::DEdge dedge) {
    assert(!isYZ(dedge.e));
    return V3minusV3I(dedge.tail()->p, dedge.head()->p, 0) < 0;
  }

  class Rib {
   public:
    Encasement3D::Face* face;
    Encasement3D::Edge *tail[2], *head[2];
    Rib *prev, *next;
    vector<PTR<Object<PV3>>> verts;
    // debug
    int id;

    Rib(Encasement3D::Face* face, Encasement3D::Edge* tail,
        Encasement3D::Edge* head)
        : face(face) {
      this->tail[0] = tail;
      this->tail[1] = 0;
      this->head[0] = head;
      this->head[1] = 0;
      prev = next = 0;
    }

    Rib(Encasement3D::Face* face, Encasement3D::DEdge dtail,
        Encasement3D::DEdge dhead)
        : face(face) {
      assert(isTail(dtail));
      assert(!isTail(dhead));
      assert(!isYZ(dtail.e));
      assert(!isYZ(dhead.e));
      this->tail[0] = dtail.e;
      this->tail[1] = 0;
      this->head[0] = dhead.e;
      this->head[1] = 0;
      prev = next = 0;
    }
  };

  class FaceMesh {
   public:
    FaceMesh(Mesh3D* mesh3d, Encasement3D::Face* face)
        : mesh3d(mesh3d), face(face), signFz(0) {
      poly = face->enc->bf[face->jf];
      approximateEdges();
      createRibs();
      findNeighbors();
      generateTriangles();
    }

    Mesh3D* mesh3d;
    Encasement3D::Face* face;
    PTR<Object<Poly3D>> poly;
    int signFz;
    vector<Rib*> ribs;
    vector<Triangle> triangles;

    // Save rib if other endpoint is chain vertex instead of extra point.
    map<Object<PV3>*, Rib*> duplicates;

    void approximateEdges() {
      Encasement3D::Loop* loop = face->loop;
      do {
        Encasement3D::DEdge dedge = loop->dedge;
        do {
          Encasement3D::Edge* edge = dedge.e;
          if (signFz == 0 && edge->b->type != Encasement3D::FFzB)
            signFz = P3derIV3(poly, 2, edge->b->p);
          if (mesh3d->edge2chain.find(edge) == mesh3d->edge2chain.end())
            mesh3d->edge2chain[edge] = mesh3d->approximate(edge);
          dedge = dedge.next();
        } while (dedge != loop->dedge);
        loop = loop->next;
      } while (loop != face->loop);
      if (signFz == 0) {
        PTR<Object<PV3>> pof = face->enc->getPointOnFace(face);
        signFz = P3derIV3(poly, 2, pof);
      }
      assert(signFz != 0);
    }

    void createRibs() {
      Encasement3D::Loop* loop = face->loop;
      do {
        Encasement3D::DEdge start = loop->dedge;
        while (isYZ(start.e)) {
          start = start.next();
          assert(start != loop->dedge);
        }
        assert(!isYZ(start.e));
        Encasement3D::DEdge dedge = start;
        do {
          if (!isYZ(dedge.e)) {
            vector<PTR<Object<PV3>>>& chain = mesh3d->edge2chain[dedge.e];
            if (chain[0] == dedge.tail()->p)
              for (int i = 1; i < chain.size() - 1; ++i) {
                PTR<Object<PV3>> v = chain[i];
                if (isExtra(v)) continue;
                if (duplicates.find(v) != duplicates.end()) continue;
                makeRib(dedge, v, vert3DInfo(dedge.e->b));
              }
            else
              for (int i = chain.size() - 1; --i >= 1;) {
                PTR<Object<PV3>> v = chain[i];
                if (isExtra(v)) continue;
                if (duplicates.find(v) != duplicates.end()) continue;
                makeRib(dedge, v, vert3DInfo(dedge.e->b));
              }
            makeHeadRibs(dedge);
          }
          dedge = dedge.next();
        } while (dedge != start);
        loop = loop->next;
      } while (loop != face->loop);
    }

    bool isExtra(Object<PV3>* v) { return dynamic_cast<V3V3atSI*>(v) != 0; }

    Vert3DInfo vert3DInfo(Encasement3D::Vert* v) {
      return Vert3DInfo(v->type, v->jf, v->jg, v->jh);
    }

    Rib* pushRib(Rib* rib) {
      int nribs = ribs.size();
      // debug
      rib->id = nribs;
      ribs.push_back(rib);
      return rib;
    }

    void makeRib(Encasement3D::DEdge dedge, PTR<Object<PV3>> v,
                 Vert3DInfo vert3DInfo) {
      // Make an FXEncasement with x=v.x.

      // debug
      static int count;
      ++count;

      FXEncasement fxe(face->enc, face->jf, v, vert3DInfo);

      PTR<Object<Scalar>> x = new V3I(v, 0);

      // Find the cel c that contains (v.y, v.z).
      FXEncasement::Cel* c = fxe.identicalVert->c[0];

      Encasement3D::DEdge tail, head;
      bool vIsTail = isTail(dedge);
      FXEncasement::Edge *first = 0, *last = 0;
      if (vIsTail) {
        tail = dedge;
        first = last = c->spv.v->e[0];
        if (P2derIV2(fxe.f, 1, first->b->p) != signFz) {
          assert(c->spv.v->twin != 0);
          first = last = c->spv.v->twin->e[0];
          assert(P2derIV2(fxe.f, 1, first->b->p) == signFz);
        }
        while (last->v[1]->type == FXEncasement::FD && last->v[1]->jg == 0)
          last = last->v[1]->e[0];
        head = face->enc->getDEdge(last->v[1], x, face->jf, 1);
      } else {
        head = dedge;
        first = last = c->spv.v->e[1];
        if (P2derIV2(fxe.f, 1, first->b->p) != signFz) {
          if (c->spv.v->twin == 0) {
            PV3<Parameter> vp = v->getApprox();
            cout << "PANIC!"
                 << " " << vp.x.mid() << " " << vp.y.mid() << " " << vp.z.mid()
                 << endl;
          } else {
            assert(c->spv.v->twin != 0);
            first = last = c->spv.v->twin->e[1];
            assert(P2derIV2(fxe.f, 1, first->b->p) == signFz);
          }
        }
        while (first->v[0]->type == FXEncasement::FD && first->v[0]->jg == 0)
          first = first->v[0]->e[1];
        tail = face->enc->getDEdge(first->v[0], x, face->jf, 0);
      }

      Rib* rib = pushRib(new Rib(face, tail, head));

      PTR<Object<PV3>> extra = getExtra(vIsTail ? head.e : tail.e, x);
      if (!isExtra(extra)) duplicates[extra] = rib;

      rib->verts.push_back(vIsTail ? v : extra);
      FXEncasement::Edge* edge = first;
      for (FXEncasement::Edge* edge = first; edge->v[0] != last->v[1];
           edge = edge->v[1]->e[0]) {
        vector<PTR<Object<PV2>>> chain = mesh3d->approximate(&fxe, edge);
        for (int i = 1; i < chain.size() - (edge == last); ++i)
          rib->verts.push_back(new XV2(x, chain[i]));
        if (edge == last) break;
      }
      rib->verts.push_back(vIsTail ? extra : v);
    }

    PTR<Object<PV3>> getExtra(Encasement3D::Edge* edge, PTR<Object<Scalar>> x) {
      vector<PTR<Object<PV3>>>& chain = mesh3d->edge2chain[edge];
      int a = -1, b = chain.size();
      int dir = V3minusV3I(chain[b - 1], chain[a + 1], 0);
      while (a + 1 < b) {
        int m = (a + b) / 2;
        int s = V3IminusS(chain[m], 0, x) * dir;
        if (s < 0)
          a = m;
        else if (s > 0)
          b = m;
        else
          return chain[m];
      }
      assert(0 <= a && b < chain.size());
      PTR<Object<PV3>> p = new V3V3atSI(chain[a], chain[b], x, 0);
      chain.insert(chain.begin() + b, p);
      return p;
    }

    // FGD!!!!
    void makeHeadRibs(Encasement3D::DEdge dedge) {
      Encasement3D::Vert* vert = dedge.head();
      PTR<Object<PV3>> v = vert->p;
      bool vIsTail = isTail(dedge);
      Encasement3D::DEdge dnext = dedge.next();

      if (!isXturning(dedge)) {
        if (!isYZ(dedge.next().e)) {
          Rib* rib = 0;
          if (duplicates.find(dedge.head()->p) != duplicates.end())
            rib = duplicates[dedge.head()->p];
          else {
            makeRib(dedge, dedge.head()->p, vert3DInfo(vert));
            rib = ribs.back();
          }
          assert(isTail(dnext) == vIsTail);
          if (vIsTail) {
            rib->tail[0] = dedge.e;
            rib->tail[1] = dnext.e;
          } else {
            rib->head[0] = dedge.e;
            rib->head[1] = dnext.e;
          }
          return;
        }

        // dnext is parallel to yz plane
        vector<PTR<Object<PV3>>> yzchain;
        yzchain.push_back(v);
        do {
          vector<PTR<Object<PV3>>>& chain = mesh3d->edge2chain[dnext.e];
          if (chain[0] == yzchain.back())
            yzchain.insert(yzchain.end(), chain.begin() + 1, chain.end());
          else {
            assert(chain.back() == yzchain.back());
            for (int i = chain.size() - 1; --i >= 0;)
              yzchain.push_back(chain[i]);
          }
          dnext = dnext.next();
        } while (isYZ(dnext.e));

        Encasement3D::DEdge tail = vIsTail ? dedge : dnext;
        Encasement3D::DEdge head = vIsTail ? dnext : dedge;
        Rib* rib = pushRib(new Rib(face, tail, head));
        if (vIsTail)
          ribs.back()->tail[1] = tail.e;
        else
          ribs.back()->head[1] = head.e;
        if (vIsTail)
          rib->verts = yzchain;
        else
          for (int i = yzchain.size(); --i >= 0;)
            rib->verts.push_back(yzchain[i]);
        return;
      }

      Encasement3D::DEdge tail = vIsTail ? dedge : dnext;
      Encasement3D::DEdge head = vIsTail ? dnext : dedge;
      assert(!isYZ(tail.e) && isTail(tail));
      assert(!isYZ(head.e) && !isTail(head));

      if (!isSplitOrJoin(dedge)) {
        // local minimum or maximum in x
        // zero-length rib
        Rib* rib = pushRib(new Rib(face, tail, head));
        rib->verts.push_back(v);
        if (vIsTail)
          rib->tail[1] = tail.e;
        else
          rib->head[1] = head.e;
        return;
      }

      // split or join
      makeRib(head, v, vert3DInfo(vert));
      Rib* ribt = ribs.back();
      makeRib(tail, v, vert3DInfo(vert));
      Rib* ribh = ribs.back();

      Rib* rib = pushRib(new Rib(face, ribt->tail[0], ribh->head[0]));
      rib->verts = ribt->verts;
      for (int i = 1; i < ribh->verts.size(); ++i)
        rib->verts.push_back(ribh->verts[i]);
    }

    bool isXturning(Encasement3D::DEdge dedge) {
      return !isYZ(dedge.e) && !isYZ(dedge.next().e) &&
             isTail(dedge) != isTail(dedge.next());
#ifdef BLEEN
      return ((v->type == Encasement3D::FGD && v->jf >= 4) ||
              v->type == Encasement3D::FFzD) &&
             v->jh == 0;
#endif
    }

    bool isSplitOrJoin(Encasement3D::DEdge dedge) {
      assert(isXturning(dedge));
      Encasement3D::Vert* vert = dedge.head();
      if (vert->type == Encasement3D::FFzD) {
        if (P3derIV3(poly, 1, vert->p) == 0)
          return P3yyzz_yzyz(poly, vert->p) < 0;
        else
          return P3zzz(poly, vert->p) == signFz;
      }

      if (vert->type == Encasement3D::FGD) {
        if (vert->jf > 3) {
          Object<Poly3D>* fgh[3];
          dedge.e->b->getFGH(mesh3d->enc3d, fgh);
          return isTail(dedge) ^ (signFz > 0) ^ (dedge.fb() == 0) ^
                 (DFxDGIV3(fgh[0], fgh[1], 1, vert->p) > 0);
        } else {
          assert(vert->jf > 1);
          int lu = (vert->jf - 2) % 2;
          assert(vert->jg >= 6);
          PTR<Object<Poly3D>> f = mesh3d->enc3d->bf[vert->jg];
          int sfy = P3derIV3(f, 1, vert->p);
          PTR<Object<Poly3D>> fz = new P3derI(f, 2);
          int sfzz = P3derIV3(fz, 2, vert->p);
          return lu ^ (sfy > 0) ^ (sfzz > 0);
        }
      }
      return false;
    }

    bool isXextreme(Encasement3D::DEdge dedge) {
      assert(0);
      Encasement3D::Vert* vert = dedge.head();
      assert(isXturning(dedge));
      if (vert->type == Encasement3D::FFzD)
        return P3derIV3(poly, 1, vert->p) == 0 &&
               P3yyzz_yzyz(poly, vert->p) > 0;
      if (vert->jf < 4) return false;
      Object<Poly3D>* fgh[3];
      dedge.e->b->getFGH(mesh3d->enc3d, fgh);
      return !(isTail(dedge) ^ (signFz > 0) ^ (dedge.fb() == 0) ^
               (DFxDGIV3(fgh[0], fgh[1], 1, vert->p) > 0));
    }

    void findNeighbors() {
      std::sort(ribs.begin(), ribs.end(), Compare(signFz));
      for (int i = 0; i < ribs.size(); ++i) {
        Rib* rib = ribs[i];
        if (rib->tail[0] == rib->tail[1]) {
          assert(rib->prev != 0);
          continue;
        }
        if (rib->prev != 0 && i >= 2 &&
            isSplitJoin(ribs[i - 2], ribs[i - 1], ribs[i]))
          continue;
        if (rib->prev != 0 && i < ribs.size() - 2 &&
            isSplitJoin(ribs[i], ribs[i + 1], ribs[i + 2])) {
          assert(ribs[i + 1]->prev != 0);
          ++i;
          continue;
        }
        for (int j = i + 1; j < ribs.size(); ++j) {
          Rib* rib2 = ribs[j];
          if (commonTail(rib, rib2) && commonHead(rib, rib2)) {
            rib->next = rib2;
            rib2->prev = rib;
            break;
          }
        }
        assert(rib->next);
      }
    }

    bool isSplitJoin(Rib* t, Rib* h, Rib* th) {
      return (t->tail[0] == th->tail[0] && t->head[0] != th->head[0] &&
              h->tail[0] != th->tail[0] && h->head[0] == th->head[0]);
    }

    class Compare {
      int signFz;

     public:
      Compare(int signFz) : signFz(signFz) {}
      bool operator()(Rib* a, Rib* b) {
        if (a == b) return false;
        int s = V3minusV3I(a->verts[0], b->verts[0], 0);
        if (s < 0) return true;
        if (s > 0) return false;
        s = signFz * V3minusV3I(a->verts.back(), b->verts.back(), 1);
        if (s < 0) return true;
        if (s > 0) return false;
        return signFz * V3minusV3I(a->verts[0], b->verts[0], 1) > 0;
      }
    };

    Encasement3D::Edge* commonTail(Rib* a, Rib* b) {
      if (a->tail[0] == b->tail[0] || a->tail[0] == b->tail[1])
        return a->tail[0];
      if (a->tail[1] == b->tail[0] || a->tail[1] == b->tail[1])
        return a->tail[1];
      return 0;
    }

    Encasement3D::Edge* commonHead(Rib* a, Rib* b) {
      if (a->head[0] == b->head[0] || a->head[0] == b->head[1])
        return a->head[0];
      if (a->head[1] == b->head[0] || a->head[1] == b->head[1])
        return a->head[1];
      return 0;
    }

    void generateTriangles() {
      // generateTriangles(ribs[167], ribs[169]);
      // if (true) return;
      for (Rib* rib : ribs) {
        if (rib->next != 0) generateTriangles(rib, rib->next);
      }
    }

    class EdgeTriangle {
     public:
      int triangleIndex;         // index of triangle in triangles
      Encasement3D::Edge* edge;  // edge containing side of triangle
      int vertexIndex;  // (lowest) index of a vertex of triangle in chain[edge]

      EdgeTriangle(int triangleIndex, Encasement3D::Edge* edge, int vertexIndex)
          : triangleIndex(triangleIndex),
            edge(edge),
            vertexIndex(vertexIndex) {}
    };

    vector<EdgeTriangle> edgeTriangles;

    // Find lowest index of a vertex of t in chain[e] starting search at
    // startingIndex.
    int findIndex(Triangle& t, Encasement3D::Edge* e, int startingIndex = 0) {
      vector<PTR<Object<PV3>>>& chain = mesh3d->edge2chain[e];
      for (int i = startingIndex; i < chain.size(); ++i)
        for (int j = 0; j < 3; ++j)
          if (chain[i] == t.p[j]) return i;
      assert(0);
      return -1;
    }

    void generateTriangles(Rib* l, Rib* r) {
      vector<PTR<Object<PV3>>>& u = l->verts;
      vector<PTR<Object<PV3>>>& v = r->verts;

      int dir = u.size() > 1 ? V3minusV3I(u.back(), u[0], 1)
                             : V3minusV3I(v.back(), v[0], 1);
      assert(dir != 0);

      int firstTriangle = triangles.size();
      int i = 0, j = 0;

      while (i < u.size() - 1 || j < v.size() - 1)
        if (i < u.size() - 1 && (j == v.size() - 1 ||
#ifdef IGNORE_Z
                                 dir * V3minusV3I(u[i], v[j], 1) < 0
#else
                                 V3V3mV3V3(u[i + 1], v[j], u[i], v[j + 1]) < 0
#endif
                                 )) {
          triangles.push_back(Triangle(u[i], v[j], u[i + 1]));
          ++i;
        } else {
          triangles.push_back(Triangle(u[i], v[j], v[j + 1]));
          ++j;
        }

      Encasement3D::Edge* tail = commonTail(l, r);
      edgeTriangles.push_back(EdgeTriangle(
          firstTriangle, tail, findIndex(triangles[firstTriangle], tail)));
      Encasement3D::Edge* head = commonHead(l, r);
      edgeTriangles.push_back(EdgeTriangle(triangles.size() - 1, head,
                                           findIndex(triangles.back(), head)));
    }

    int findIndex(Triangle& t, PTR<Object<PV3>> p) {
      for (int i = 0; i < 3; ++i)
        if (t.p[i] == p) return i;
      assert(0);
      return -1;
    }

    void remesh() {
      for (int i = 0; i < edgeTriangles.size(); ++i) {
        EdgeTriangle& et = edgeTriangles[i];
        // refresh value of et.vertexIndex, which may have changed
        et.vertexIndex =
            findIndex(triangles[et.triangleIndex], et.edge, et.vertexIndex);
        int nextIndex =
            findIndex(triangles[et.triangleIndex], et.edge, et.vertexIndex + 1);
        if (nextIndex == et.vertexIndex + 1) continue;
        triangles.push_back(triangles[et.triangleIndex]);
        Triangle& t1 = triangles[et.triangleIndex];
        Triangle& t2 = triangles[triangles.size() - 1];
        vector<PTR<Object<PV3>>>& chain = mesh3d->edge2chain[et.edge];
        t1.p[findIndex(t1, chain[nextIndex])] = chain[et.vertexIndex + 1];
        t2.p[findIndex(t2, chain[et.vertexIndex])] = chain[et.vertexIndex + 1];
      }
    }
  };

  map<Encasement3D::Face*, FaceMesh*> face2mesh;

  const vector<Triangle>& getMesh(Encasement3D::Face* face) {
    if (face2mesh.find(face) == face2mesh.end())
      face2mesh[face] = new FaceMesh(this, face);
    face2mesh[face]->remesh();
    return face2mesh[face]->triangles;
  }

  void writeMesh(const char* filename, const vector<Triangle>& triangles) {
    vector<Object<PV3>*> points;
    map<Object<PV3>*, int> index;

    for (Triangle t : triangles)
      for (int i = 0; i < 3; ++i) {
        Object<PV3>* p = t.p[i];
        if (index.find(p) == index.end()) {
          index[p] = points.size();
          points.push_back(p);
        }
      }

    ofstream out;
    out.open(filename);

    out << "# vtk DataFile Version 3.0" << endl;
    out << "vtk output" << endl;
    out << "ASCII" << endl;
    out << "DATASET POLYDATA" << endl;
    out << "POINTS " << points.size() << " double" << endl;

    for (Object<PV3>* p : points) {
      PV3<Parameter> pp = p->getApprox(1);
      out << pp.x.mid() << " " << pp.y.mid() << " " << pp.z.mid() << endl;
    }

    out << endl;
    out << "POLYGONS " << triangles.size() << " " << 4 * triangles.size()
        << endl;

    for (Triangle t : triangles) {
      out << 3;
      for (int i = 0; i < 3; ++i) out << " " << index[t.p[i]];
      out << endl;
    }

    out.close();
  }

  // FXEncasement Chains

  // Create a polygonal chain that approximates FXEncasement edge.
  vector<PTR<Object<PV2>>> approximate(FXEncasement* encfx,
                                       FXEncasement::Edge* edge) {
    vector<FXEncasement::Vert*> verts = initialChain(edge);
    assert(verts[0] == edge->v[0]);
    assert(verts.back() == edge->v[1]);

    vector<PTR<Object<PV2>>> chain;
    chain.push_back(verts[0]->p);

    for (int i = 1; i < verts.size(); ++i)
      addPoint(encfx, edge, chain, verts[i]->p);

    return chain;
  }

  // edge->v[0], boundary verts, edge->v[1]
  vector<FXEncasement::Vert*> initialChain(FXEncasement::Edge* edge) {
    vector<FXEncasement::Vert*> v;
    vector<FXEncasement::Vert*> v0 = halfChain(edge->b, edge->b->c[0]);
    vector<FXEncasement::Vert*> v1 = halfChain(edge->b, edge->b->c[1]);

    vector<FXEncasement::Vert*>* pv0 =
        v0.back() == edge->v[0] || v0.back()->twin == edge->v[0] ? &v0 : &v1;
    vector<FXEncasement::Vert*>* pv1 =
        v0.back() == edge->v[0] || v0.back()->twin == edge->v[0] ? &v1 : &v0;
    assert(pv0->back() == edge->v[0] || pv0->back()->twin == edge->v[0]);
    (*pv0)[pv0->size() - 1] = edge->v[0];
    assert(pv1->back() == edge->v[1] || pv1->back()->twin == edge->v[1]);
    (*pv1)[pv1->size() - 1] = edge->v[1];

    for (int i = pv0->size(); --i >= 0;) v.push_back((*pv0)[i]);
    v.push_back(edge->b);
    for (int i = 0; i < pv1->size(); ++i) v.push_back((*pv1)[i]);

    return v;
  }

  // Returns chain of boundary Vert plus one arrangement Vert on Edge
  // starting with boundary vert b and going in direction of Cel c.
  // b is not in the output.
  vector<FXEncasement::Vert*> halfChain(FXEncasement::Vert* b,
                                        FXEncasement::Cel* c) {
    vector<FXEncasement::Vert*> verts;
    while (c->spv.v == 0) {
      b = c->lb.b == b ? c->rb.b : c->lb.b;
      c = b->c[0] == c ? b->c[1] : b->c[0];
      verts.push_back(b);
    }
    verts.push_back(c->spv.v);
    return verts;
  }

  // Add q to a chain that lies within epsilon of the curve.
  void addPoint(FXEncasement* encfx, FXEncasement::Edge* edge,
                vector<PTR<Object<PV2>>>& chain, PTR<Object<PV2>> q) {
    PTR<Object<PV2>> p = chain.back();
    PTR<Object<PV2>> r = sample(encfx, edge, p, q);
    if (V2distV2V2mS(r, p, q, epsilon) < 0) {
      chain.push_back(r);
      chain.push_back(q);
    } else {
      addPoint(encfx, edge, chain, r);
      addPoint(encfx, edge, chain, q);
    }
  }

  // Sample curve "halfway" between p and q.
  PTR<Object<PV2>> sample(FXEncasement* encfx, FXEncasement::Edge* edge,
                          PTR<Object<PV2>> p, PTR<Object<PV2>> q) {
    Object<Poly2D>* line = 0;
    int side = 0;
    FXEncasement::Cel* cel = 0;
    PV2<Parameter> bbox = getBBox(encfx, edge, p, q, line, side, cel);
    int yzMax = maxDiffCoord(p, q);
    if (cel->spv.v != 0 && cel->spv.v->type == FXEncasement::FD) {
      // box contains a y or z turning point
      if (cel->spv.v->jg == !yzMax) yzMax = !yzMax;
      // box contain a singular point
      if (cel->spv.v->twin != 0 && cel->spv.v->twin->type == FXEncasement::FD)
        yzMax = 0;
    }
    double mid = getMiddle(p, q, yzMax);
    vector<PTR<Object<PV2>>> roots = encfx->getRoots(yzMax, mid, bbox);
    assert(roots.size() > 0);
    int edgeSignFz = P2derIV2(encfx->f, 1, edge->b->p);  // calculate just once?
    PTR<Object<PV2>> r;
    int n = 0;
    for (PTR<Object<PV2>> root : roots)
      if (P2derIV2(encfx->f, 1, root) == edgeSignFz &&
          (line == 0 || Poly2Dval(line, root) == side)) {
        r = root;
        ++n;
      }
    assert(n == 1);
    return r;
  }

  PV2<Parameter> getBBox(FXEncasement* encfx, FXEncasement::Edge* edge,
                         PTR<Object<PV2>> p, PTR<Object<PV2>> q,
                         Object<Poly2D>*& line, int& side,
                         FXEncasement::Cel*& cel) {
    PTR<Object<PV2>> m = new MidPoint2(p, q);
    PV2<Parameter> box = encfx->box;
    return getBBox(&encfx->top, m, box, line, side, cel);
  }

  PV2<Parameter> getBBox(FXEncasement::Cel* c, PTR<Object<PV2>> m,
                         PV2<Parameter> box, Object<Poly2D>*& line, int& side,
                         FXEncasement::Cel*& cel) {
    if (0 <= c->yz && c->yz < 2) {
      if (V2IminusD(m, c->yz, c->spv.s) < 0)
        return getBBox(c->lb.c, m, leftBox(box, c->yz, c->spv.s), line, side,
                       cel);
      else
        return getBBox(c->rb.c, m, rightBox(box, c->yz, c->spv.s), line, side,
                       cel);
    }

    if (c->yz == 2) {
      assert(line == 0);
      line = c->spv.p;
      side = Poly2Dval(line, m);
      if (side < 0)
        return getBBox(c->lb.c, m, box, line, side, cel);
      else
        return getBBox(c->rb.c, m, box, line, side, cel);
    }

    assert(c->yz == -1);
    cel = c;
    return box;
  }

  int maxDiffCoord(PTR<Object<PV2>> a, PTR<Object<PV2>> b) {
    return V2V2ImJ(a, b, 0, 1) > 0 ? 0 : 1;
  }

  double getMiddle(double a, double b) { return a + (b - a) * 0.500123456789; }

  double getMiddle(PTR<Object<PV2>> u, PTR<Object<PV2>> v, int yz) {
    double m =
        getMiddle(u->getApprox(1.0)[yz].ub(), v->getApprox(1.0)[yz].lb());
    if (V2IminusD(u, yz, m) * V2IminusD(v, yz, m) == -1) return m;
    m = getMiddle(u->getApprox()[yz].ub(), v->getApprox()[yz].lb());
    assert(V2IminusD(u, yz, m) * V2IminusD(v, yz, m) == -1);
    return m;
  }

  // Encasement3D Chains

  // Create a polygonal chain that approximates Encasement3D edge.
  vector<PTR<Object<PV3>>> approximate(Encasement3D::Edge* edge) {
    vector<Encasement3D::Vert*> verts = initialChain(edge);
    assert(verts[0] == edge->v[0]);
    assert(verts.back() == edge->v[1]);

    vector<PTR<Object<PV3>>> chain;
    chain.push_back(verts[0]->p);

    // Use verts if it already is an epsilon-chain.
    if (verts.size() == 3 &&
        V3distV3V3mS(verts[1]->p, verts[0]->p, verts[2]->p, epsilon) < 0) {
      chain.push_back(verts[1]->p);
      chain.push_back(verts[2]->p);
      return chain;
    }

#ifdef START_WITH_BOUNDARY_POINTS
    for (int i = 1; i < verts.size(); ++i)
      addPoint(edge, chain, verts[i]->p, verts);
#else
    addPoint(edge, chain, verts.back()->p, verts);
#endif

#ifndef START_WITH_BOUNDARY_POINTS
    // removing boundary points to avoid identities (sigh)
    int j = 1;
    int n = 0;
    for (int i = 0; i < chain.size(); ++i)
      if (chain[i] != verts[j]->p || j == verts.size() - 1)
        chain[n++] = chain[i];
      else
        ++j;
    chain.resize(n);
#endif
    return chain;
  }

  // edge->v[0], boundary verts, edge->v[1]
  vector<Encasement3D::Vert*> initialChain(Encasement3D::Edge* edge) {
    vector<Encasement3D::Vert*> v0 = halfChain(edge->b, edge->b->c[0]);
    vector<Encasement3D::Vert*> v1 = halfChain(edge->b, edge->b->c[1]);

    vector<Encasement3D::Vert*>* pv0 = edge->v[0] == v0.back() ? &v0 : &v1;
    vector<Encasement3D::Vert*>* pv1 = edge->v[0] == v0.back() ? &v1 : &v0;
    assert(pv0->back() == edge->v[0]);
    assert(pv1->back() == edge->v[1]);
    vector<Encasement3D::Vert*> v;

    for (int i = pv0->size(); --i >= 0;) v.push_back((*pv0)[i]);
    v.push_back(edge->b);
    for (int i = 0; i < pv1->size(); ++i) v.push_back((*pv1)[i]);

    return v;
  }

  // Returns chain of boundary Vert plus one arrangement Vert on Edge
  // starting with boundary vert b and going in direction of b->c[i].
  // b is not in the output.
  vector<Encasement3D::Vert*> halfChain(Encasement3D::Vert* b,
                                        Encasement3D::Cel* c) {
    vector<Encasement3D::Vert*> verts;
    while (c->spv.v == 0) {
      b = c->lb.b == b ? c->rb.b : c->lb.b;
      c = b->c[0] == c ? b->c[1] : b->c[0];
      verts.push_back(b);
    }
    verts.push_back(c->spv.v);
    return verts;
  }

  // Add q to a chain that lies within epsilon of the curve.
  void addPoint(Encasement3D::Edge* edge, vector<PTR<Object<PV3>>>& chain,
                PTR<Object<PV3>> q, vector<Encasement3D::Vert*>& verts) {
    PTR<Object<PV3>> p = chain.back();
    PTR<Object<PV3>> r = sample(edge, p, q, verts);
    if (V3distV3V3mS(r, p, q, epsilon) < 0) {
      chain.push_back(r);
      chain.push_back(q);
    } else {
      addPoint(edge, chain, r, verts);
      addPoint(edge, chain, q, verts);
    }
  }

  // Sample curve "halfway" between p and q.
  PTR<Object<PV3>> sample(Encasement3D::Edge* edge, PTR<Object<PV3>> p,
                          PTR<Object<PV3>> q) {
    Object<Poly3D>* plane = 0;
    int side = 0;
    Encasement3D::Cel* cel = 0;
    PV3<Parameter> bbox = getBBox(edge, p, q, plane, side, cel);

    int xyzAvoid = -1;
    if (cel->spv.v != 0 && (cel->spv.v->type == Encasement3D::FGD ||
                            cel->spv.v->type == Encasement3D::FFzD ||
                            cel->spv.v->type == Encasement3D::FSD))
      xyzAvoid = cel->spv.v->jh;
    Encasement3D::FG* fg = getFG(edge);
    if (cel->spv.v != 0 && cel->spv.v->type == Encasement3D::FFzH &&
        fg->jg == -1)
      xyzAvoid = 1;
    int xyzMax = maxDiffCoord(p, q, xyzAvoid);
    if (cel->spv.v != 0 && cel->spv.v->type == Encasement3D::FFzD &&
        fg->jg == -1) {
      if (xyzAvoid == 0)
        xyzMax = 2;
      else if (xyzAvoid == 2)
        xyzMax = 0;
    }
    if (cel->spv.v != 0 && cel->spv.v->type == Encasement3D::FGD &&
        fg->dFxdG[xyzMax] == 0) {
      assert(xyzMax != xyzAvoid);
      xyzMax = 3 - xyzMax - xyzAvoid;
    }
    double mid = getMiddle(p, q, xyzMax);
    Encasement3D::VertType bType = edge->b->type;
    vector<PTR<Object<PV3>>> roots = Encasement3D::getRoots(
        fg->getF(), fg->getG(), bType, fg->jf, fg->jg, xyzMax, mid, bbox);
    if(roots.size() == 0) {
      cout << "roots 0";
    }
    assert(roots.size() > 0);
    PTR<Object<PV3>> r;
    if (plane == 0) {
      assert(roots.size() == 1);
      r = roots[0];
    } else {
      int n = 0;
      for (PTR<Object<PV3>> root : roots)
        if (Poly3Dval(plane, root) == side) {
          r = root;
          ++n;
        }
      assert(n == 1);
    }
    return r;
  }

  // verts is a polygonal chain, and p and q are points on this chain.
  // Return the point on the chain halfway between p and q in the longest
  // dimension. Set pp and qq to the nearest points in p and q and chain
  // vertices.
  PTR<Object<PV3>> getMidPoint(PTR<Object<PV3>> p, PTR<Object<PV3>> q,
                               vector<Encasement3D::Vert*>& verts,
                               PTR<Object<PV3>>& pp, PTR<Object<PV3>>& qq) {
    // x-monotonic (xy=0) unless it lies in an x plane, then use y (xy=1)
    int xy = (V3minusV3I(p, q, 0) == 0);

    int dir = V3minusV3I(verts.back()->p, verts[0]->p, xy);

    int ip = 1;
    while (ip < verts.size() - 1 && dir * V3minusV3I(p, verts[ip]->p, xy) > 0)
      ++ip;
    if (!(dir * V3minusV3I(p, verts[ip]->p, xy) < 0)) {
      cout << "p[xy] " << p->getApprox(1.0)[xy].lb() << " "
           << p->getApprox(1.0)[xy].ub() << endl;
      cout << "verts[ip].p[xy] " << verts[ip]->p->getApprox(1.0)[xy].lb() << " "
           << verts[ip]->p->getApprox(1.0)[xy].ub() << endl;
      cout << "verts[0].p[xy] " << verts[0]->p->getApprox(1.0)[xy].lb() << " "
           << verts[0]->p->getApprox(1.0)[xy].ub() << endl;
      cout << "verts[2].p[xy] " << verts[2]->p->getApprox(1.0)[xy].lb() << " "
           << verts[2]->p->getApprox(1.0)[xy].ub() << endl;
      int s = dir * V3minusV3I(p, verts[ip]->p, 0);
      cout << "s " << s << endl;
      cout << "verts.size() " << verts.size() << endl;
      cout << "ip " << ip << endl;
    }
    assert(dir * V3minusV3I(p, verts[ip]->p, xy) < 0);
    int iq = verts.size() - 2;
    while (iq > 0 && dir * V3minusV3I(verts[iq]->p, q, xy) > 0) --iq;
    assert(dir * V3minusV3I(verts[iq]->p, q, xy) < 0);
    PTR<Object<PV3>> m = new V3V3ofT(p, q, 0.500102030405);
    if (ip > iq) {
      pp = p;
      qq = q;
      return m;
    }
    int xyz = maxDiffCoord(p, q, -1);
    // cout << "m[xyz] " << m->getApprox(1.0)[xyz].mid() << endl;
    for (int i = ip - 1; i <= iq; ++i) {
      pp = i < ip ? p : verts[i]->p;
      qq = i < iq ? verts[i + 1]->p : q;
      // cout << "pp[xyz] " << pp->getApprox(1.0)[xyz].mid() << endl;
      // cout << "qq[xyz] " << qq->getApprox(1.0)[xyz].mid() << endl;
      if (V3minusV3I(pp, m, xyz) * V3minusV3I(qq, m, xyz) < 0)
        return new V3V3atV3I(pp, qq, m, xyz);
    }
    assert(0);
    return 0;
  }

  // Sample curve "halfway" between p and q.
  PTR<Object<PV3>> sample(Encasement3D::Edge* edge, PTR<Object<PV3>> p,
                          PTR<Object<PV3>> q,
                          vector<Encasement3D::Vert*>& verts) {
    // if (V3minusV3I(p, q, 0) == 0) return sample(edge, p, q);

    PTR<Object<PV3>> pp, qq;
    PTR<Object<PV3>> m = getMidPoint(p, q, verts, pp, qq);

    Object<Poly3D>* plane = 0;
    int side = 0;
    Encasement3D::Cel* cel = 0;
    PV3<Parameter> bbox = getBBox(edge, m, plane, side, cel);

    int xyzAvoid = -1;
    if (cel->spv.v != 0 && (cel->spv.v->type == Encasement3D::FGD ||
                            cel->spv.v->type == Encasement3D::FFzD ||
                            cel->spv.v->type == Encasement3D::FSD))
      xyzAvoid = cel->spv.v->jh;
    Encasement3D::FG* fg = getFG(edge);
    if (cel->spv.v != 0 && cel->spv.v->type == Encasement3D::FFzH &&
        fg->jg == -1)
      xyzAvoid = 1;
    int xyzMax = maxDiffCoord(p, q, xyzAvoid);
    if (cel->spv.v != 0 && cel->spv.v->type == Encasement3D::FFzD &&
        fg->jg == -1) {
      if (xyzAvoid == 0)
        xyzMax = 2;
      else if (xyzAvoid == 2)
        xyzMax = 0;
    }
    if (cel->spv.v != 0 && cel->spv.v->type == Encasement3D::FGD &&
        fg->dFxdG[xyzMax] == 0) {
      assert(xyzMax != xyzAvoid);
      xyzMax = 3 - xyzMax - xyzAvoid;
    }
    double mid = getMiddle(m, pp, qq, xyzMax);
    Encasement3D::VertType bType = edge->b->type;
    vector<PTR<Object<PV3>>> roots = Encasement3D::getRoots(
        fg->getF(), fg->getG(), bType, fg->jf, fg->jg, xyzMax, mid, bbox);
    assert(roots.size() > 0);
    PTR<Object<PV3>> r;
    if (plane == 0) {
      assert(roots.size() == 1);
      r = roots[0];
    } else {
      int n = 0;
      for (PTR<Object<PV3>> root : roots)
        if (Poly3Dval(plane, root) == side) {
          r = root;
          ++n;
        }
      assert(n == 1);
    }
    return r;
  }

  PV3<Parameter> getBBox(Encasement3D::Edge* edge, PTR<Object<PV3>> p,
                         PTR<Object<PV3>> q, Object<Poly3D>*& plane, int& side,
                         Encasement3D::Cel*& cel) {
    PTR<Object<PV3>> m = new MidPoint3(p, q);
    PV3<Parameter> box = enc3d->box;
    Encasement3D::FG* fg = getFG(edge);
    return getBBox(&fg->top, m, box, plane, side, cel);
  }

  PV3<Parameter> getBBox(Encasement3D::Edge* edge, PTR<Object<PV3>> m,
                         Object<Poly3D>*& plane, int& side,
                         Encasement3D::Cel*& cel) {
    PV3<Parameter> box = enc3d->box;
    Encasement3D::FG* fg = getFG(edge);
    return getBBox(&fg->top, m, box, plane, side, cel);
  }

  Encasement3D::FG* getFG(Encasement3D::Edge* edge) {
    Encasement3D::Vert* b = edge->b;
    switch (b->type) {
      case Encasement3D::FGB:
        return enc3d->fg[enc3d->ind(b->jf, b->jg)];
      case Encasement3D::FFzB:
        return enc3d->f[b->jf]->ffz;
      case Encasement3D::FSB:
        return enc3d->f[b->jf]->fs[b->jg];
      default:
        assert(0);
    };
    return 0;
  }

  PV3<Parameter> getBBox(Encasement3D::Cel* c, PTR<Object<PV3>> m,
                         PV3<Parameter> box, Object<Poly3D>*& plane, int& side,
                         Encasement3D::Cel*& cel) {
    if (0 <= c->xyz && c->xyz < 3) {
      if (V3IminusD(m, c->xyz, c->spv.s) < 0)
        return getBBox(c->lb.c, m, leftBox(box, c->xyz, c->spv.s), plane, side,
                       cel);
      else
        return getBBox(c->rb.c, m, rightBox(box, c->xyz, c->spv.s), plane, side,
                       cel);
    }

    if (c->xyz == 3) {
      assert(plane == 0);
      plane = c->spv.p;
      side = Poly3Dval(plane, m);
      if (side < 0)
        return getBBox(c->lb.c, m, box, plane, side, cel);
      else
        return getBBox(c->rb.c, m, box, plane, side, cel);
    }

    assert(c->xyz == -1);
    cel = c;
    return box;
  }

  int maxDiffCoord(PTR<Object<PV3>> a, PTR<Object<PV3>> b, int xyzAvoid) {
    if (xyzAvoid != -1) {
      int i = (xyzAvoid + 1) % 3;
      int j = (xyzAvoid + 2) % 3;
      return V3V3ImJ(a, b, i, j) > 0 ? i : j;
    }
    int i = V3V3ImJ(a, b, 0, 1) > 0 ? 0 : 1;
    return V3V3ImJ(a, b, i, 2) > 0 ? i : 2;
  }

#ifdef BLEEN
  double getMiddle(PTR<Object<PV3>> u, PTR<Object<PV3>> v, int xyz) {
    double m = (u->getApprox(1.0)[xyz].ub() + v->getApprox(1.0)[xyz].lb()) / 2;
    if (V3IminusD(u, xyz, m) * V3IminusD(v, xyz, m) == -1) return m;
    m = (u->getApprox()[xyz].ub() + v->getApprox()[xyz].lb()) / 2;
    assert(V3IminusD(u, xyz, m) * V3IminusD(v, xyz, m) == -1);
    return m;
  }
#endif

  double getMiddle(PTR<Object<PV3>> u, PTR<Object<PV3>> v, int xyz) {
    double m =
        getMiddle(u->getApprox(1.0)[xyz].ub(), v->getApprox(1.0)[xyz].lb());
    if (V3IminusD(u, xyz, m) * V3IminusD(v, xyz, m) == -1) return m;
    m = getMiddle(u->getApprox()[xyz].ub(), v->getApprox()[xyz].lb());
    assert(V3IminusD(u, xyz, m) * V3IminusD(v, xyz, m) == -1);
    return m;
  }

  double getMiddle(PTR<Object<PV3>> mid, PTR<Object<PV3>> u, PTR<Object<PV3>> v,
                   int xyz) {
    assert(V3minusV3I(u, mid, xyz) * V3minusV3I(v, mid, xyz) < 0);
    Parameter midA = mid->getApprox(1.0)[xyz];
    double m = getMiddle(midA.lb(), midA.ub());
    assert(V3IminusD(u, xyz, m) * V3IminusD(v, xyz, m) == -1);
    return m;
  }
};
