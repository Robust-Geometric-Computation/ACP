#include "path3d.h"

const bool DEBUG = false;
const bool VERBOSE = false;

void save(Points path, const char * filename) {
  ofstream ostr;
  ostr.open(filename);
  if (ostr.is_open()) { 
    ostr << setprecision(20) << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl << "ASCII" << endl
         << "DATASET POLYDATA" << endl 
         << "POINTS " << path.size() << " double" << endl;
    for (int i=0; i<path.size(); i++)
      ostr << path[i]->getApprox(1e-16).getX().mid() << " " << path[i]->getApprox(1e-16).getY().mid() << " " << path[i]->getApprox(1e-16).getZ().mid() << endl;
    ostr<<endl<<"LINES 1 "<<path.size()+1<<endl<<path.size()<<" ";
    for (int i=0; i<path.size(); i++)
      ostr<<i<<" ";
    ostr.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
}

template<class N>
class TransformationData {
 public:
  PV3<N> u, v, w, t;
  TransformationData () {}
  TransformationData (PV3<N> u, PV3<N> v, PV3<N> w, PV3<N> t) : u(u), v(v), w(w), t(t) {}

  template<class M>
  TransformationData (const TransformationData<M> &p) : u(p.u), v(p.v), w(p.w), t(p.t) {}

  N & operator[](int i) { 
    if (i < 3) return u[i];
    if (i < 6) return v[i%3];
    if (i < 9) return w[i%3];
    return t[i%3];
  }
  const N & operator[](int i) const { 
    if (i < 3) return u[i];
    if (i < 6) return v[i%3];
    if (i < 9) return w[i%3];
    return t[i%3];
  }
  int size () const { return 12; }
};

class Transformation : public Object<TransformationData> {
 public:
  template<class N>
  PV3<N> getU () { return get<N>().u; }
  template<class N>
  PV3<N> getV () { return get<N>().v; }
  template<class N>
  PV3<N> getW () { return get<N>().w; }
  template<class N>
  PV3<N> getT () { return get<N>().t; }
};

//transformation which rotates d about bc to lie in the same plane as abc
class UnfoldTriangleTransformation : public Transformation {
  PTR<Point> pa, pb, pc, pd;
  DeclareCalculate(TransformationData) {
    PV3<N> a = pa->get<N>();
    PV3<N> b = pb->get<N>();
    PV3<N> c = pc->get<N>();
    PV3<N> d = pd->get<N>();

    PV3<N> u = (c-b).unit();
    PV3<N> uprime = u;
    PV3<N> w = (d-b).cross(u).unit();
    PV3<N> wprime = u.cross(a-b).unit();
    PV3<N> v = w.cross(u);
    PV3<N> vprime = wprime.cross(u);

    PV3<N> t = b - b.dot(u)*(uprime) - b.dot(v)*(vprime) - b.dot(w)*(wprime);
    PV3<N> x = uprime.getX()*u + vprime.getX()*v + wprime.getX()*w;
    PV3<N> y = uprime.getY()*u + vprime.getY()*v + wprime.getY()*w;
    PV3<N> z = uprime.getZ()*u + vprime.getZ()*v + wprime.getZ()*w;
    return TransformationData<N>(x,y,z,t);
  }
 public:
  UnfoldTriangleTransformation(PTR<Point> pa, PTR<Point> pb, PTR<Point> pc, PTR<Point> pd) : pa(pa), pb(pb), pc(pc), pd(pd) {}
  UnfoldTriangleTransformation(HFace * hf1, HFace * hf2) {
    PTR<Point> p1[3];
    p1[0] = hf1->getF()->getBoundary()->tail()->getP();
    p1[1] = hf1->getF()->getBoundary()->getNext()->tail()->getP();
    p1[2] = hf1->getF()->getBoundary()->getNext()->getNext()->tail()->getP();
  
    PTR<Point> p2[3];
    p2[0] = hf2->getF()->getBoundary()->tail()->getP();
    p2[1] = hf2->getF()->getBoundary()->getNext()->tail()->getP();
    p2[2] = hf2->getF()->getBoundary()->getNext()->getNext()->tail()->getP();

    int j1, j2;
    for (int i=0; i<3; i++)
      if (p1[i] != p2[0] && p1[i] != p2[1] && p1[i] != p2[2])
        j1 = i;
    
    for (int i=0; i<3; i++)
      if (p1[0] != p2[i] && p1[1] != p2[i] && p1[2] != p2[i])
        j2 = i;
   
    this->pa = p1[j1];
    this->pb = p1[(j1+1)%3];      
    this->pc = p1[(j1+2)%3];
    this->pd = p2[j2];      
  }
};

class NxNdotR : public Primitive {
  PTR<Point> pa, pb, pc, pd;
  static PTR<Object<PV3>> pr;

  DeclareSign {
    if (pr == 0)
      pr = new Object<PV3>(PV3<double>(randomNumber(-1, 1), randomNumber(-1, 1),
				       randomNumber(-1, 1)), false);
    PV3<N> a = pa->get<N>();
    PV3<N> b = pb->get<N>();
    PV3<N> c = pc->get<N>();
    PV3<N> d = pd->get<N>();

    PV3<N> u = (c-b).unit();
    PV3<N> uprime = u;
    PV3<N> w = (d-b).cross(u).unit();
    PV3<N> wprime = u.cross(a-b).unit();

    return w.cross(wprime).dot(pr->get<N>());
  }
 public:
  NxNdotR (HFace * hf1, HFace * hf2) {
    PTR<Point> p1[3];
    p1[0] = hf1->getF()->getBoundary()->tail()->getP();
    p1[1] = hf1->getF()->getBoundary()->getNext()->tail()->getP();
    p1[2] = hf1->getF()->getBoundary()->getNext()->getNext()->tail()->getP();
  
    PTR<Point> p2[3];
    p2[0] = hf2->getF()->getBoundary()->tail()->getP();
    p2[1] = hf2->getF()->getBoundary()->getNext()->tail()->getP();
    p2[2] = hf2->getF()->getBoundary()->getNext()->getNext()->tail()->getP();

    int j1, j2;
    for (int i=0; i<3; i++)
      if (p1[i] != p2[0] && p1[i] != p2[1] && p1[i] != p2[2])
        j1 = i;
    
    for (int i=0; i<3; i++)
      if (p1[0] != p2[i] && p1[1] != p2[i] && p1[2] != p2[i])
        j2 = i;
   
    this->pa = p1[j1];
    this->pb = p1[(j1+1)%3];      
    this->pc = p1[(j1+2)%3];
    this->pd = p2[j2];      
  }
};

PTR<Object<PV3>> NxNdotR::pr;

//transformation which rotates triangle abc into the xy-plane centered around b
class XYPlaneTriangleTransfromation : public Transformation {
  PTR<Point> pa, pb, pc, pd;
  DeclareCalculate(TransformationData) {
    PV3<N> a = pa->get<N>();
    PV3<N> b = pb->get<N>();
    PV3<N> c = pc->get<N>();

    PV3<N> u = (c-b).unit();
    PV3<N> w = (a-b).cross(u).unit();
    PV3<N> v = u.cross(w);
    PV3<N> t = b - PV3<N>(b.dot(u), b.dot(v), b.dot(w));

    return TransformationData<N>(u,v,w,t);
  }
 public:
  XYPlaneTriangleTransfromation(PTR<Point> pa, PTR<Point> pb, PTR<Point> pc) : pa(pa), pb(pb), pc(pc) {}
};

class CompositeTransformation : public Transformation {
  PTR<Transformation> f, g;
  DeclareCalculate(TransformationData) {
    PV3<N> u1 = g->getU<N>();
    PV3<N> v1 = g->getV<N>();
    PV3<N> w1 = g->getW<N>();
    PV3<N> t1 = g->getT<N>();
    PV3<N> u2 = f->getU<N>();
    PV3<N> v2 = f->getV<N>();
    PV3<N> w2 = f->getW<N>();
    PV3<N> t2 = f->getT<N>();

    PV3<N> x1 = PV3<N>(u1.getX(), v1.getX(), w1.getX());
    PV3<N> y1 = PV3<N>(u1.getY(), v1.getY(), w1.getY());
    PV3<N> z1 = PV3<N>(u1.getZ(), v1.getZ(), w1.getZ());

    PV3<N> u = PV3<N>(x1.dot(u2), y1.dot(u2), z1.dot(u2)); 
    PV3<N> v = PV3<N>(x1.dot(v2), y1.dot(v2), z1.dot(v2)); 
    PV3<N> w = PV3<N>(x1.dot(w2), y1.dot(w2), z1.dot(w2)); 
    PV3<N> t = PV3<N>(u2.dot(t1), v2.dot(t1), w2.dot(t1)) + t2;

    return TransformationData<N>(u,v,w,t);
  }
 public:
  CompositeTransformation(PTR<Transformation> f, PTR<Transformation> g) : f(f), g(g) {}
};

//apply transformation t to point p
class TransformedPoint : public Point {
  PTR<Point> point;
  PTR<Transformation> t;
  DeclareCalculate(PV3) {
    PV3<N> p = point->get<N>();
    PV3<N> u = t->getU<N>();
    PV3<N> v = t->getV<N>();
    PV3<N> w = t->getW<N>();
    PV3<N> tt = t->getT<N>();
    return PV3<N>(u.dot(p), v.dot(p), w.dot(p)) + tt;
  }
 public:
  TransformedPoint(PTR<Point> point, PTR<Transformation> t) : point(point), t(t) {}
};

struct ComparePointOrder {
  PTR<Point> r;
  ComparePointOrder(PTR<Point> r) : r(r) {}
  bool operator()(PTR<FaceIntersectionPoint> i, PTR<FaceIntersectionPoint> j) {
    return (PointOrder(i, j, r) > 0); 
  }
};

class XYComponents : public Object<PV2> {
 protected:
  PTR<Point> p;
  DeclareCalculate(PV2) { return PV2<N>(p->get<N>().getX(), p->get<N>().getY()); }
 public:
  XYComponents (PTR<Point> p) : p(p) {} 
};

class PathVertex {
 public:
  PTR<Point> original;
  PTR<Object<PV2> > transformed2d;
  PathVertex(PTR<Point> p) {
    original = p; transformed2d = 0;
  }
  PathVertex(PTR<Point> original, PTR<Object<PV2> > transformed2d) : original(original), transformed2d(transformed2d) {}
};

void save(std::vector<PathVertex*> & vertPath, const char * filename) {
  ofstream ostr;
  ostr.open(filename);
  if (ostr.is_open()) { 
    ostr << setprecision(20) << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl << "ASCII" << endl
         << "DATASET POLYDATA" << endl 
         << "POINTS " << vertPath.size() << " double" << endl;
    for (int i=0; i<vertPath.size(); i++) {
      ostr << vertPath[i]->transformed2d->getApprox(1e-6).getX().mid()<<" "<<vertPath[i]->transformed2d->getApprox(1e-6).getY().mid()<<" 0"<<endl;
    }
    ostr<<endl<<"LINES 1 "<<vertPath.size()+1<<endl<<vertPath.size()<<" ";
    for (int i=0; i<vertPath.size(); i++)
      ostr<<i<<" ";
    ostr.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
}

class AreaABC : public Primitive {
  PathVertex * pva, * pvb, *pvc;

  DeclareSign {
    PV2<N> a = pva->transformed2d->get<N>();
    PV2<N> b = pvb->transformed2d->get<N>();
    PV2<N> c = pvc->transformed2d->get<N>();
    return (b-a).cross(c-a);
  }

 public:
  AreaABC (PathVertex* pva, PathVertex* pvb, PathVertex* pvc) : pva(pva), pvb(pvb), pvc(pvc) {}
};

//find the 3d equivalent of the 2d intersection between triangle edge p->q and line segment c->d
class ABintersectCDto3D : public Point {
 protected:
  PathVertex *vp, *vq, *vc, *vd;
 public:
  ABintersectCDto3D(PathVertex * vp, PathVertex * vq, PathVertex * vc, PathVertex * vd)
   : vp(vp), vq(vq), vc(vc), vd(vd) {}
  DeclareCalculate(PV3) {
    PV2<N> p = vp->transformed2d->get<N>();
    PV2<N> q = vq->transformed2d->get<N>();
    PV2<N> c = vc->transformed2d->get<N>();
    PV2<N> d = vd->transformed2d->get<N>();
    PV3<N> p3d = vp->original->get<N>();
    PV3<N> q3d = vq->original->get<N>();

    N t = -((p-c).cross(d-c)) / ((q-p).cross(d-c));
    assert(t.sign() > 0);
    assert((t - N(1)).sign() < 0);
    return p3d + t *(q3d-p3d);
  }
};

void savePivots(Points path1, const char * filename) {
  Points path;
  for (int i=0; i<path1.size(); i++) {
    Point * p = path1[i];
    if ((dynamic_cast<ABintersectCDto3D *> (p)) == NULL)
      path.push_back(path1[i]);
  }

  ofstream ostr;
  ostr.open(filename);
  if (ostr.is_open()) { 
    ostr << setprecision(20) << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl << "ASCII" << endl
         << "DATASET POLYDATA" << endl 
         << "POINTS " << path.size() << " double" << endl;
    for (int i=0; i<path.size(); i++)
      ostr << path[i]->getApprox(1e-16).getX().mid() << " " << path[i]->getApprox(1e-16).getY().mid() << " " << path[i]->getApprox(1e-16).getZ().mid() << endl;
    ostr<<endl<<"LINES 1 "<<path.size()+1<<endl<<path.size()<<" ";
    for (int i=0; i<path.size(); i++)
      ostr<<i<<" ";
    ostr.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
}

class PathTriangle {
 public:
  PathVertex * p[3];
  HFace * hface;
  PTR<Transformation> cumulative;
  PathVertex * getVertex(PTR<Point> q) {
    for (int i=0; i<3; i++)
      if (p[i]->original == q)
        return p[i];
    return NULL;
  }
  PathTriangle(PathVertex * p0, PathVertex * p1, PathVertex * p2, HFace * hface) {
    p[0] = p0; p[1] = p1; p[2] = p2; this->hface = hface;
  }
  PathTriangle(HFace * hf,  PathTriangle prev) {
    PTR<Point> ps[3];
    ps[0] = hf->getF()->getBoundary()->tail()->getP();
    ps[1] = hf->getF()->getBoundary()->getNext()->tail()->getP();
    ps[2] = hf->getF()->getBoundary()->getNext()->getNext()->tail()->getP();
    
    int uncommon = -1;
    for (int i=0; i<3; i++)
      if (ps[i] != prev.p[0]->original  &&  ps[i] != prev.p[1]->original  &&  ps[i] != prev.p[2]->original)
        uncommon = i;
    assert(uncommon > -1);

    PTR<Point> p = ps[uncommon];
    ps[uncommon] = ps[2];
    ps[2] = p;

    this->p[0] = prev.getVertex(ps[0]); 
    this->p[1] = prev.getVertex(ps[1]); 
    this->p[2] = new PathVertex(ps[2]);
    this->hface = hf;
  }
};

void savePathTriangles(std::vector<PathTriangle> & ts, const char * filename, bool saveTranformed) {
  ofstream ostr;
  ostr.open(filename);
  if (ostr.is_open()) {
    ostr << setprecision(20) << "# vtk DataFile Version 3.0" << endl
         << "vtk output" << endl << "ASCII" << endl
         << "DATASET POLYDATA" << endl 
         << "POINTS " << ts.size()*3 << " double" << endl;
    for (int i=0; i<ts.size(); i++)
      for (int j=0; j<3; j++)
        if (saveTranformed)
          ostr << ts[i].p[j]->transformed2d->getApprox().getX().mid() << " " << ts[i].p[j]->transformed2d->getApprox().getY().mid() << " 0.0" << endl;
        else
          ostr << ts[i].p[j]->original->getApprox().getX().mid() << " " << ts[i].p[j]->original->getApprox().getY().mid() << " " << ts[i].p[j]->original->getApprox().getZ().mid() << endl;

    ostr << endl << "POLYGONS " << ts.size() << " " << 4*ts.size() << endl;
    for (int i=0; i<ts.size(); i++)
      ostr << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << endl;
    ostr.close();
  } else { cout<<"could not open file"<<endl; return; }
}

void flattenTriangles(std::vector<PathTriangle> & triangles, int n) {
  PTR<Transformation> cumulative = triangles[n-1].cumulative;
  for (int i=n; i<triangles.size(); i++) {
    if (triangles[i].p[0]->transformed2d == 0)
      triangles[i].p[0]->transformed2d = new XYComponents(new TransformedPoint(triangles[i].p[0]->original, cumulative));
    if (triangles[i].p[1]->transformed2d == 0)
      triangles[i].p[1]->transformed2d = new XYComponents(new TransformedPoint(triangles[i].p[1]->original, cumulative));

    if (NxNdotR(triangles[i-1].hface, triangles[i].hface) != 0) {
      PTR<Transformation> t = new UnfoldTriangleTransformation(triangles[i-1].hface, triangles[i].hface);
      cumulative = new CompositeTransformation(cumulative, t);
    }
    triangles[i].p[2]->transformed2d = new XYComponents(new TransformedPoint(triangles[i].p[2]->original, cumulative));
    triangles[i].cumulative = cumulative;
  }
}

void commonEdge(PathTriangle t1, PathTriangle t2, PathVertex *& p, PathVertex *& q) {
  int uncommon = -1;
  int commonCount = 0;
  for (int i=0; i<3; i++)
    if (t1.p[i] != t2.p[0] && t1.p[i] != t2.p[1] && t1.p[i] != t2.p[2])
      uncommon = i;
    else
      commonCount++;
  assert(uncommon > -1); assert(commonCount == 2);
  p = t1.p[(uncommon+1)%3];
  q = t1.p[(uncommon+2)%3];
}

class PathEdge {
 public:
  PathVertex * left;
  PathVertex * right;
  PathEdge(PathVertex * left, PathVertex * right) : left(left), right(right) {}
};

bool hasVertexPoint(HFace * hf, PTR<Point> p) {
  PTR<Point> p0 = hf->getF()->getBoundary()->tail()->getP();
  PTR<Point> p1 = hf->getF()->getBoundary()->getNext()->tail()->getP();
  PTR<Point> p2 = hf->getF()->getBoundary()->getNext()->getNext()->tail()->getP();
  return ((p0 == p) || (p1 == p) || (p2 == p));
}

void path2Dto3D(std::vector<PathVertex*> &vertPath, std::vector<int> &vertPathIndices, std::vector<PathEdge> &edges,  Points &path) {
  path.push_back(vertPath[0]->original);
  int j = 1;
  for (int i=0; i<edges.size(); i++) {
    if (j == vertPath.size() && (i == edges.size()-1) && edges[i].right == vertPath[vertPath.size()-1]) 
      break;
    assert(j < vertPath.size());
    if (vertPathIndices[j] == i) {
      path.push_back(vertPath[j]->original);
      if (vertPath[j] == edges[edges.size()-1].right)
        break;
      j++;
    } else if (AreaABC(edges[i].left, edges[i].right, vertPath[j])   != 0
            && AreaABC(edges[i].left, edges[i].right, vertPath[j-1]) != 0 ) {
      path.push_back(new ABintersectCDto3D(edges[i].left, edges[i].right, vertPath[j-1], vertPath[j]));
      path[path.size()-1]->getApprox();
    }
  }
}

class State {
 public:
  std::vector<PathEdge> edges;
  std::vector<PathVertex*> left;
  std::vector<PathVertex*> right;
  int nPath;
  State(std::vector<PathEdge> & e, int nedges, std::vector<PathVertex*> & l, std::vector<PathVertex*> r, int n) {
    if (nedges > e.size())
      nedges = e.size();
    edges.insert(edges.begin(), e.begin(), e.begin()+nedges);
    right.insert(right.begin(), r.begin(), r.end());
    left.insert(left.begin(), l.begin(), l.end());
    nPath = n;
  }
  void save(const char * filename) {
    ofstream ostr;
    ostr.open(filename);
    if (ostr.is_open()) {
      ostr << setprecision(20)<< edges.size()<<" "<<nPath<<endl;
      ostr << left.size()<<" "<<right.size()<<endl;
      for (int i=0; i<edges.size(); i++) {
        ostr<<edges[i].left->transformed2d->getApprox(1e-16).getX().mid()<<" "<<edges[i].left->transformed2d->getApprox(1e-16).getY().mid()<<endl;
        ostr<<edges[i].right->transformed2d->getApprox(1e-16).getX().mid()<<" "<<edges[i].right->transformed2d->getApprox(1e-16).getY().mid()<<endl;
      }
      for (int i=0; i<left.size(); i++)
        ostr<<left[i]->transformed2d->getApprox(1e-16).getX().mid()<<" "<<left[i]->transformed2d->getApprox(1e-16).getY().mid()<<endl;
      for (int i=0; i<right.size(); i++)
        ostr<<right[i]->transformed2d->getApprox(1e-16).getX().mid()<<" "<<right[i]->transformed2d->getApprox(1e-16).getY().mid()<<endl;
      ostr.close();
    } else { cout<<"could not open file"<<endl; }
  }
};

void shortestPath(std::vector<PathTriangle> & triangles, PathVertex * start, PathVertex * end, std::vector<PathVertex*> & left, std::vector<int> & leftIndices, std::vector<PathEdge> & edges) {
  for (int i=0; i<triangles.size()-1; i++) {
    PathVertex * p, * q;
    commonEdge(triangles[i], triangles[i+1], p, q);
    PathVertex * from;
    for (int j=0; j<3; j++)
      if (triangles[i].p[j] != p && triangles[i].p[j] != q)
        from = triangles[i].p[j];
    PathVertex * l = ((AreaABC(from, p, q) > 0)? q : p);
    PathVertex * r = ((l != p)? p : q);
    edges.push_back(PathEdge(l, r)); 
  }
  if (edges[0].left->original == start->original || edges[0].left->original->identical(start->original))
    start = edges[0].left;
  if (edges[0].right->original == start->original || edges[0].right->original->identical(start->original))
    start = edges[0].right;
  if (edges[edges.size()-1].left->original == end->original || edges[edges.size()-1].left->original->identical(end->original))
    end = edges[edges.size()-1].left;
  if (edges[edges.size()-1].right->original == end->original || edges[edges.size()-1].right->original->identical(end->original))
    end = edges[edges.size()-1].right;
  edges.push_back(PathEdge(end, end));

  std::vector<PathVertex*> right;
  std::vector<int> rightIndices;
  left.push_back(start); leftIndices.push_back(-1);
  right.push_back(start); rightIndices.push_back(-1);

  if (edges[0].left != start) { left.push_back(edges[0].left); leftIndices.push_back(0); }
  if (edges[0].right != start) { right.push_back(edges[0].right); rightIndices.push_back(0); }
  char s[50]; int state_num = 0;

  if (DEBUG) {
    sprintf(s, "state/%d.txt", state_num++);
    State(edges, 2, left, right, 1).save(s);
  }

  int nPath = 1; //left[0] and right[0] are the same, the rest of the path isn't
  for (int i=1; i<edges.size(); i++) {
    if (right[right.size()-1] != edges[i].right) {
      while(right.size()>nPath && AreaABC(right[right.size()-2], right[right.size()-1], edges[i].right) > 0)
      {  right.pop_back(); rightIndices.pop_back();  }

      if (right.size() == nPath) {
        while (left.size()>nPath && AreaABC(right[right.size()-1], edges[i].right, left[nPath]) < 0) {
          right.push_back(left[nPath]);
          rightIndices.push_back(leftIndices[nPath]);
          nPath++;
        }
      }
      right.push_back(edges[i].right);
      rightIndices.push_back(i);

      if (DEBUG) {
        sprintf(s, "state/%d.txt", state_num++);
        State(edges, i+2, left, right, nPath).save(s);
      }
    }

    if (left[left.size()-1] != edges[i].left) {
      while(left.size()>nPath && AreaABC(left[left.size()-2], left[left.size()-1], edges[i].left) < 0)
      {  left.pop_back(); leftIndices.pop_back();  }

      if (left.size() == nPath) {
        while (right.size() > nPath && AreaABC(left[left.size()-1], edges[i].left, right[nPath]) > 0) {
          left.push_back(right[nPath]);
          leftIndices.push_back(rightIndices[nPath]);
          nPath++;
        }
      }
      left.push_back(edges[i].left);
      leftIndices.push_back(i);

      if (DEBUG) {
        sprintf(s, "state/%d.txt", state_num++);
        State(edges, i+2, left, right, nPath).save(s);
      }
    }
  }
}

//populate newPath with a sequence of PathTriangles from oldPath[startIndex] to oldPathpendIndex] with goes
//the other way around pv
void otherWay(vector<PathTriangle> & oldPath, vector<PathTriangle> & newPath, int startIndex, int endIndex, PathVertex * pv) {
  HFace * nextHF = NULL;
  HFace * hfaceStart = oldPath[startIndex].hface;
  std::vector<HFace*> neighbors = hfaceStart->neighbors();
  assert(neighbors.size() == 3);
  
  if (endIndex != startIndex+1 && (neighbors[0] == oldPath[endIndex].hface || neighbors[1] == oldPath[endIndex].hface || neighbors[2] == oldPath[endIndex].hface))
    return;

  for (int i=0; i< neighbors.size(); i++)
    if (neighbors[i] != oldPath[startIndex-1].hface && neighbors[i] != oldPath[startIndex+1].hface && hasVertexPoint(neighbors[i],pv->original))
      nextHF = neighbors[i];
  assert(nextHF != NULL);
  PathTriangle next(nextHF, oldPath[startIndex]);
  newPath.push_back(next);
  while (true) {
    neighbors.clear();
    HFace * prev = ((newPath.size()==1)? oldPath[startIndex].hface : newPath[newPath.size()-2].hface);
    neighbors = newPath[newPath.size()-1].hface->neighbors();
    assert(neighbors.size() == 3);
    for (int i=0; i<neighbors.size(); i++)
      if (neighbors[i] != prev && hasVertexPoint(neighbors[i],pv->original))
        nextHF = neighbors[i];
    if (nextHF == oldPath[endIndex].hface)
      return;
    PathTriangle next(nextHF, newPath[newPath.size()-1]);
    newPath.push_back(next); 
  }
}



void localPath(PTR<Point> a, PTR<Point> b, HFaces & pathfaces, Points & path) {
  vector<long> randomLongs;
  set<long> pathHashes;
  std::vector<PathTriangle> triangles;
  std::vector<PTR<Transformation> > transformations;
  HFace * hf = pathfaces[0];
  PathVertex * v0 = new PathVertex(hf->getF()->getBoundary()->tail()->getP());
  PathVertex * v1 = new PathVertex(hf->getF()->getBoundary()->getNext()->tail()->getP());
  PathVertex * v2 = new PathVertex(hf->getF()->getBoundary()->getNext()->getNext()->tail()->getP());
  triangles.push_back(PathTriangle(v0,v1,v2,hf));
  for (int i=1; i<pathfaces.size(); i++) {
    PathTriangle prev= triangles[i-1];
    hf = pathfaces[i];

    PTR<Point> p[3];
    p[0] = hf->getF()->getBoundary()->tail()->getP();
    p[1] = hf->getF()->getBoundary()->getNext()->tail()->getP();
    p[2] = hf->getF()->getBoundary()->getNext()->getNext()->tail()->getP();

    int uncommon1 = -1;
    for (int j=0; j<3; j++)
      if (prev.p[j]->original != p[0] && prev.p[j]->original != p[1] && prev.p[j]->original != p[2])
        uncommon1 = j;
    assert(uncommon1 > -1);

    int uncommon2 = -1;
    for (int j=0; j<3; j++)
      if (p[j] != prev.p[0]->original && p[j] != prev.p[1]->original && p[j] != prev.p[2]->original)
        uncommon2 = j;
    assert(uncommon2 > -1);

    PathTriangle t(prev.p[(uncommon1+1)%3], prev.p[(uncommon1+2)%3], new PathVertex(p[uncommon2]), hf);
    triangles.push_back(t);
    // transformations.push_back(new UnfoldTriangleTransformation(prev.p[uncommon1]->original, t.p[0]->original, t.p[1]->original, t.p[2]->original));  
  } 
  PTR<Transformation> xyplane = new XYPlaneTriangleTransfromation(triangles[0].p[0]->original, triangles[0].p[1]->original, triangles[0].p[2]->original);
  for (int i=0; i<3; i++)
    triangles[0].p[i]->transformed2d = new XYComponents(new TransformedPoint(triangles[0].p[i]->original, xyplane));
  triangles[0].cumulative = xyplane;
  PathVertex * start = new PathVertex((PTR<Point>) a);
  start->transformed2d = new XYComponents(new TransformedPoint(start->original, xyplane));
  flattenTriangles(triangles, 1);
  PathVertex * end = new PathVertex((PTR<Point>) b);
  end->transformed2d = new XYComponents(new TransformedPoint(end->original, triangles[triangles.size()-1].cumulative));

  std::vector<PathVertex*> vertPath; 
  std::vector<int> vertPathIndices; 
  std::vector<PathEdge> edges;

  shortestPath(triangles, start, end, vertPath, vertPathIndices, edges);

  if (DEBUG) {
    path2Dto3D(vertPath, vertPathIndices, edges, path);
    save(path, "path0.vtk");
    save(vertPath, "path0-2d.vtk");
    savePathTriangles(triangles, "triangles.vtk", false);
    savePathTriangles(triangles, "triangles-2d.vtk", true);
  }  

  bool changedThisIteration;
  int iter_num  = 0;
  char s[50];
  do {
    iter_num++;
    changedThisIteration = false;
    if (VERBOSE) cout<<"new iteration------------------ "<<iter_num<<endl;
    for (int i=0; i<triangles.size(); i++) {
      
      if (vertPathIndices[vertPathIndices.size()-2] < i)
        continue;
      int vertPathPos;
      for (int j=0; j<vertPath.size(); j++) {
        if (vertPathIndices[j] >= i) {
          i = vertPathIndices[j];
          vertPathPos = j;
          break;
        }
      }

      int startIndex = i;
      int endIndex = startIndex+1;
      while(edges[endIndex].left == vertPath[vertPathPos] || edges[endIndex].right == vertPath[vertPathPos])
        endIndex++;

      std:vector<PathTriangle> newPath;
      if (VERBOSE) cout<<"startIndex: "<<startIndex<<" endIndex: "<<endIndex<<endl;

      PathVertex * replacingVertex = vertPath[vertPathPos];
      otherWay(triangles, newPath, startIndex, endIndex, replacingVertex);

      if (DEBUG) {
        std::vector<PathTriangle> tmp;
        tmp.insert(tmp.begin(), triangles.begin()+startIndex, triangles.begin()+endIndex+1);
        sprintf(s, "old%d-%02d-2d.vtk", iter_num, i);
        savePathTriangles(tmp, s, true);
        sprintf(s, "old%d-%02d.vtk", iter_num, i);
        savePathTriangles(tmp, s, false);
      }

      triangles.erase(triangles.begin()+startIndex+1, triangles.begin()+endIndex);
      triangles.insert(triangles.begin()+startIndex+1, newPath.begin(), newPath.end());
      endIndex = startIndex + newPath.size() + 1;

      //give triangles[endIndex] new PathVertices to match the new intermediate path
      //also replace the one that doesn't match with a new vertex, it may match something earlier in the path
      int uncommon = -1;
      int commonCount = 0;
      for (int j=0; j<3; j++) {
        int matchingIndex = -1;
        for (int k=0; k<3; k++)
          if (triangles[endIndex].p[j]->original == triangles[endIndex-1].p[k]->original)
            matchingIndex = k;
        if (matchingIndex == -1) {
          uncommon = j;
          triangles[endIndex].p[j] = new PathVertex(triangles[endIndex].p[j]->original);
        } else {
          triangles[endIndex].p[j] = triangles[endIndex-1].p[matchingIndex];
          commonCount++;
        }
      }
      assert(uncommon > -1); assert(commonCount == 2);
      //swap the the triangle[endIndex]'s vertices so that the uncommon 1 is at position 2
      PathVertex * v_tmp = triangles[endIndex].p[uncommon];
      triangles[endIndex].p[uncommon] = triangles[endIndex].p[2];
      triangles[endIndex].p[2] = v_tmp;
      //replace these vertices in an subsequent triangles 
      for (int j=0; j<3; j++) {
        PathVertex * v = triangles[endIndex].p[j];
        for (int k=endIndex+1; k<triangles.size(); k++) {
          if (triangles[k].p[0]->original == v->original)
            triangles[k].p[0] = v;
          else if (triangles[k].p[1]->original == v->original)
            triangles[k].p[1] = v;
          else if (triangles[k].p[2]->original == v->original)
            triangles[k].p[2] = v;
          else
            break;
        }
      }

      flattenTriangles(triangles, startIndex+1);
      end->transformed2d = new XYComponents(new TransformedPoint(end->original, triangles[triangles.size()-1].cumulative));

      if (DEBUG) {
        std::vector<PathTriangle> tmp;
        tmp.insert(tmp.begin(), triangles.begin()+startIndex, triangles.begin()+endIndex+1);
        sprintf(s, "new%d-%02d-2d.vtk", iter_num, i);
        savePathTriangles(tmp, s, true);
        sprintf(s, "new%d-%02d.vtk", iter_num, i);
        savePathTriangles(tmp, s, false);
      }

      std::vector<PathVertex*> oldPath;
      oldPath.insert(oldPath.begin(), vertPath.begin(), vertPath.end());

      vertPath.clear();
      vertPathIndices.clear();
      edges.clear();

      shortestPath(triangles, start, end, vertPath, vertPathIndices, edges);

      if (DEBUG) {
        sprintf(s, "triangles%d-%02d-2d.vtk", iter_num, i);
        savePathTriangles(triangles, s, true);
        sprintf(s, "triangles%d-%02d.vtk", iter_num, i);
        savePathTriangles(triangles, s, false);
        sprintf(s, "path%d-%02d-2d.vtk", iter_num, i);
        save(vertPath, s);
      }

      //to do: simplify this
      bool pathChanged = false;
      long hash = 0;
      for (int j = 0; j < vertPath.size(); j++) {
	while (j >= randomLongs.size())
	  randomLongs.push_back(random());
	hash += randomLongs[j] * (vertPath[j]->original - (Point*) 0);
      }

      if (pathHashes.find(hash) == pathHashes.end()) {
	pathChanged = true;
	pathHashes.insert(hash);
      }

      /*
      if (oldPath.size() != vertPath.size())
        pathChanged = true;
      else 
        for (int j=0; j<vertPath.size(); j++)
          if (vertPath[j]->original != oldPath[j]->original && !vertPath[j]->original->identical(oldPath[j]->original))
            pathChanged = true;
      */
      changedThisIteration = changedThisIteration || pathChanged;

      if (!pathChanged && VERBOSE)
        cout<<"no change"<<endl;

      if (DEBUG && pathChanged) {
        path.clear();
        path2Dto3D(vertPath, vertPathIndices, edges, path);
        sprintf(s, "path%d-%02d.vtk", iter_num, i);
        save(path, s);
        sprintf(s, "pivot%d-%02d.vtk", iter_num, i);
        savePivots(path, s);
      }

    }

    if (iter_num > 90) {
      cout << iter_num << " triangles";
      for (int i = 0; i < triangles.size(); i++)
	cout << " " << triangles[i].hface;
      cout << endl;
    }
    if (iter_num > 100) {
      cout<<"WARN: path3d not terminating"<<endl;
      break;
    }
  } while(changedThisIteration);

  path.clear();
  path2Dto3D(vertPath, vertPathIndices, edges, path);
}

void bfs(HFace * fa, HFace * fb, HFaces & pathfaces) {
  assert(fa->getS() == fb->getS());
  HFaces pathfaces_rev;

  std::map<HFace *, HFace *> parents;
  std::queue<HFace *> q;
  q.push(fa);
  parents[fa] = NULL;

  while (!q.empty()) {
    HFace * current = q.front(); q.pop();
    if (current == fb) {
      while(current != NULL) {
        pathfaces_rev.push_back(current);
        current = parents[current];
      }
      break;
    }

    HFaces hfaces = current->neighbors();
    for (int i=(hfaces.size()-1); i>=0; i--) {
      HFace * hf = hfaces[i];
      if (parents.find(hf) == parents.end()) {
        parents[hf] = current;
        q.push(hf);
      }
    }
  }
  if (pathfaces_rev.size() == 0) {
    cout<<"no path found by bfs"<<endl;
    return;
  }

  for (int i=pathfaces_rev.size()-1; i>=0; i--)
    pathfaces.push_back(pathfaces_rev[i]);
}

void testBoundaryCondition(Polyhedron * poly, int mode) {
  int START_ON_VERTEX = 0;
  int START_ON_EDGE = 1;
  int END_ON_VERTEX = 2;
  int END_ON_EDGE = 3;
  poly->computeWindingNumbers();

  PTR<FaceIntersectionPoint> start  = 0;
  PTR<FaceIntersectionPoint> end  = 0;
  while(start == 0) {
    std::vector<PTR<FaceIntersectionPoint> > points;
    //choose a random 2 points, they define a ray
    PTR<Point> p = new Point(random()/double(RAND_MAX), random()/double(RAND_MAX),random()/double(RAND_MAX));
    PTR<Point> q = new Point(random()/double(RAND_MAX), random()/double(RAND_MAX),random()/double(RAND_MAX));
    PTR<Point> r =  new DiffPoint(p, q);
    //choose a random and shell
    int cell_index = floor(random()/double(RAND_MAX) * poly->cells.size());
    Cell * cell = poly->cells[cell_index];
    int shell_index = floor(random()/double(RAND_MAX) * cell->nShells());


    //do ray casting
    Shell * shell = cell->getShell(shell_index);
    for (int j=0; j<shell->getHFaces().size(); j++) {
      HFace * hface = shell->getHFaces()[j];
      PTR<FaceIntersectionPoint> s = new FaceIntersectionPoint(p, q, hface);
      bool onFace = hface->getF()->contains(s, false);
      if (onFace)
        points.push_back(s);
    }

    //sort all the intersections, set start and end to the first 2 if they exist
    if (points.size() > 1) {
      std::sort(points.begin(), points.end(), ComparePointOrder(r));
      start = points[0];
      end = points[1];
    }
  }

  //get bfs path
  HFaces bfsPath;
  bfs(start->getHFace(), end->getHFace(), bfsPath);
  assert(bfsPath.size() !=0);
  Points path;

  if (mode == START_ON_VERTEX) {
    //replace start with one of its face's vertices
    HFace * hface = start->getHFace();
    PTR<Point> newStart = hface->getF()->getBoundary()->tail()->getP();
    localPath(newStart, (PTR<Point>) end, bfsPath, path);
  } else if (mode  == START_ON_EDGE) {
    //replace start with a point on one of its face's edges
    HFace * hface = start->getHFace();
    PTR<Point> a = hface->getF()->getBoundary()->tail()->getP();
    PTR<Point> b = hface->getF()->getBoundary()->head()->getP();
    PTR<Point> newStart = new MidPoint(a,b);
    localPath(newStart, (PTR<Point>) end, bfsPath, path);
  } else if (mode == END_ON_VERTEX) {
    //replace end with one of its face's vertices
    HFace * hface = end->getHFace();
    PTR<Point> newEnd = hface->getF()->getBoundary()->tail()->getP();
    localPath((PTR<Point>) start, newEnd, bfsPath, path);
  } else if (mode  == END_ON_EDGE) {
    //replace end with a point on one of its face's edges
    HFace * hface = end->getHFace();
    PTR<Point> a = hface->getF()->getBoundary()->tail()->getP();
    PTR<Point> b = hface->getF()->getBoundary()->head()->getP();
    PTR<Point> newEnd = new MidPoint(a,b);
    localPath((PTR<Point>) start, newEnd, bfsPath, path);
  }
}

void findPath(Polyhedron * blockspace, int cell_index, PTR<Point> start, PTR<Point> end, Points &path) {
  //TO DO: make shortestPath algorithm can handle when start or end is a vertex

  assert(blockspace->cells.size() != 0);
  Cell * cell = blockspace->cells[cell_index];
  bool startOutside = blockspace->containingCell(start) != cell_index;
  bool endOutside = blockspace->containingCell(end) != cell_index;

  //find nearest vertex to any surface points
  PTR<Point> p_start = NULL;
  PTR<Point> p_end = NULL;
  HFace * hf_start = NULL;
  HFace * hf_end = NULL;
  if (startOutside || endOutside) {
    for (int i=0; i<cell->nShells(); i++) {
      Shell * shell = cell->getShell(i);
      for (int j=0; j<shell->getHFaces().size(); j++) {
        HFace * hface = shell->getHFaces()[j];
        if (startOutside) {
          PTR<Point> q = new FaceNearestPoint(start, hface);
          if (p_start == NULL || CloserPair(start, q, start, p_start) > 0 ) {
            p_start = q;
            hf_start = hface;
          }
        }
        if (endOutside) {
          PTR<Point> q = new FaceNearestPoint(end, hface);
          if (p_end == NULL || CloserPair(end, q, end, p_end) > 0 ) {
            p_end = q;
            hf_end = hface;
          }
        }
      }
    }
  }

  if (startOutside)
    start = p_start;
  if (endOutside)
    end = p_end;

  std::vector<PTR<FaceIntersectionPoint> > points;
  PTR<Point> r =  new DiffPoint(end, start);

  for (int i=0; i<cell->nShells(); i++) {
    Shell * shell = cell->getShell(i);
    for (int j=0; j<shell->getHFaces().size(); j++) {
      HFace * hface = shell->getHFaces()[j];

      if (hface == hf_start || hface == hf_end)
        continue;
      PTR<FaceIntersectionPoint> q = new FaceIntersectionPoint(start, end, hface);
      bool onFace = hface->getF()->contains(q, false);
      bool afterStart = (PointOrder(start, q, r) == 1);
      bool beforeEnd = (PointOrder(q, end, r) == 1);
      if (onFace && afterStart && beforeEnd)
        points.push_back(q);
    }
  }

  std::sort(points.begin(), points.end(), ComparePointOrder(r));
  if (VERBOSE) cout<<"intersections: "<<points.size()<<endl;
  
  path.push_back(start);
  if (points.size() > 0) {
    PTR<Point> p = new MidPoint(start, points[0]);
    if (startOutside && blockspace->containingCell(p) != cell_index) {
      HFaces subHfaces;
      Points subPath;
      bfs(hf_start, points[0]->getHFace(), subHfaces);
      localPath(start, (PTR<Point>)points[0], subHfaces, subPath);
      path.insert(path.end(), subPath.begin()+1, subPath.end());
    } 
    for (int i=1; i<points.size(); i++) {
      p = new MidPoint(points[i-1], points[i]);
      if (blockspace->containingCell(p) != cell_index) {
        //findpath points[i-1] to points[i]
        HFaces subHfaces;
        Points subPath;
        bfs(points[i-1]->getHFace(), points[i]->getHFace(), subHfaces);
        localPath((PTR<Point>)points[i-1], (PTR<Point>)points[i], subHfaces, subPath);
        path.insert(path.end(), subPath.begin()+1, subPath.end());
      } else
        path.push_back((PTR<Point>)points[i-1]);
    }
    p = new MidPoint(points[points.size()-1], end);
    if (endOutside && blockspace->containingCell(p) != cell_index) {
      //findpath start to points[0]
      HFaces subHfaces;
      Points subPath;
      bfs(points[points.size()-1]->getHFace(), hf_end, subHfaces);
      localPath((PTR<Point>)points[points.size()-1], end, subHfaces, subPath);
      path.insert(path.end(), subPath.begin()+1, subPath.end());
    } else 
      path.push_back((PTR<Point>)points[points.size()-1]);
  } else if (startOutside && endOutside) {
    PTR<Point> p = new MidPoint(start, end);
    if (blockspace->containingCell(p) != cell_index) {
      HFaces subHfaces; Points subPath;
      bfs(hf_start, hf_end, subHfaces);
      localPath(start, end, subHfaces, subPath);
      path.insert(path.end(), subPath.begin()+1, subPath.end());
    }
  }

  path.push_back(end);  
}

void findPath(Polyhedron * blockspace, PTR<Point> start, PTR<Point> end, Points &path) {
  //make sure start and end are not in blockspace
  if (blockspace->contains(start)) {
    cout<<"invalid start"<<endl;
    return;
  }
  if (blockspace->contains(end)) {
    cout<<"invalid end"<<endl;
    return;
  }
  
  //make sure start and end are in the same component of blockspace
  blockspace->computeWindingNumbers();
  int cell_index = blockspace->containingCell(start);
  if (blockspace->containingCell(end) != cell_index) {
    cout<<"start and end are not in the same component"<<endl;
    return;
  }

  //to do: triangulate? can we count on blockspaces already being triangulated?

  //do raycasting
  Cell * cell = blockspace->cells[cell_index];
  std::vector<PTR<FaceIntersectionPoint> > points;
  PTR<Point> r =  new DiffPoint(end, start);
  for (int i=0; i<cell->nShells(); i++) {
    Shell * shell = cell->getShell(i);
    for (int j=0; j<shell->getHFaces().size(); j++) {
      HFace * hface = shell->getHFaces()[j];
      PTR<FaceIntersectionPoint> q = new FaceIntersectionPoint(start, end, hface);
      bool onFace = hface->getF()->contains(q, false);
      bool afterStart = (PointOrder(start, q, r) == 1);
      bool beforeEnd = (PointOrder(q, end, r) == 1);
      if (onFace && afterStart && beforeEnd)
        points.push_back(q);
    }
  }
  
  if (VERBOSE) cout<<"intersections: "<<points.size()<<endl;
  assert(points.size() % 2 == 0);

  std::sort(points.begin(), points.end(), ComparePointOrder(r));

  path.push_back(start);
  for (int i=0; i<points.size(); i+=2) {
    HFaces subHfaces;
    Points subPath;
    bfs(points[i]->getHFace(), points[i+1]->getHFace(), subHfaces);
    localPath((PTR<Point>)points[i], (PTR<Point>)points[i+1], subHfaces, subPath);
    path.insert(path.end(), subPath.begin(), subPath.end());
  }
  if (points.size()>0)
    path.push_back((PTR<Point>) points[points.size()-1]);
  path.push_back(end);
}
