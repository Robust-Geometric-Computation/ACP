#include "io.h"

Polyhedron * readPolyhedron (istream &istr, bool perturb)
{
  Polyhedron *a = new Polyhedron;
  int nv, nf;
  istr >> nv >> nf;
  acp::disable();
  for (int i = 0; i < nv; ++i) {
    double x, y, z;
    istr >> x >> y >> z;
    a->getVertex(new Point(x, y, z, perturb));
  }
  acp::enable();
  for (int i = 0; i < nf; ++i) {
    int u, v, w;
    istr >> u >> v >> w;
    Vertex *uv = a->vertices[u], *vv = a->vertices[v], *wv = a->vertices[w];
    a->addTriangle(uv, vv, wv);
  }
  istr >> ws;
  if (istr.peek() != EOF)
    readCells(istr, a);
  return a;
}

void readCells (istream &istr, Polyhedron *a)
{
  int nc;
  istr >> nc;
  for (int i = 0; i < nc; ++i)
    a->cells.push_back(readCell(istr, a));
}

Cell * readCell (istream &istr, Polyhedron *a)
{
  int ns;
  istr >> ns;
  Shells sh;
  for (int i = 0; i < ns; ++i)
    sh.push_back(readShell(istr, a));
  int i0 = a->cells.empty() ? 0 : 1;
  Cell *c = new Cell(i0 == 0 ? 0 : sh[0]);
  for (int i = i0; i < ns; ++i)
    c->addInner(sh[i]);
  return c;
}

Shell * readShell (istream &istr, Polyhedron *a)
{
  int nf;
  istr >> nf;
  HFaces hf;
  for (int i = 0; i < nf; ++i) {
    int f, h = 0;
    istr >> f;
    if (f < 0) {
      f = - f;
      h = 1;
    }    
    hf.push_back(a->faces[f-1]->getHFace(h));
  }
  return new Shell(hf);
}

void writePolyhedron (Polyhedron *a, ostream &ostr)
{
  VIMap vimap;
  ostr << setprecision(17);
  ostr << a->vertices.size() << " " << a->faces.size() << endl << endl;
  acp::disable();
  for (int i = 0; i < a->vertices.size(); ++i) {
    Vertex *v = a->vertices[i];
    vimap.insert(VIPair(v, i));
    PV3<double> p = v->getP()->getApproxMid();
    ostr << p.x << " " << p.y << " " << p.z << endl;
  }
  acp::enable();
  ostr << endl;
  for (Faces::iterator f = a->faces.begin(); f != a->faces.end(); ++f) {
    Vertices ve = (*f)->getBoundary()->loop();
    ostr << vimap.find(ve[0])->second << " " << vimap.find(ve[1])->second
	 << " " << vimap.find(ve[2])->second << endl;
  }
  ostr << endl;
  if (!a->cells.empty())
    writeCells(a, ostr);  
}

void writeCells (Polyhedron *a, ostream &ostr)
{
  FIMap fimap;
  for (int i = 0; i < a->faces.size(); ++i)
    fimap.insert(FIPair(a->faces[i], i + 1));
  int nc = a->cells.size();
  ostr << nc << endl << endl;
  for (Cells::iterator c = a->cells.begin(); c != a->cells.end(); ++c) {
    int ns = (*c)->nShells();
    ostr << ns << endl;
    for (int j = 0; j < ns; ++j)
      writeShell((*c)->getShell(j), fimap, ostr);
    ostr << endl;
  }
}

void writeShell (Shell *s, FIMap &fimap, ostream &ostr)
{
  const HFaces &hf = s->getHFaces();
  ostr << hf.size() << endl;
  for (HFaces::const_iterator h = hf.begin(); h != hf.end(); ++h) {
    int s = (*h)->pos() ? 1 : -1, f = fimap.find((*h)->getF())->second;
    ostr << s*f << " ";
  }
  ostr << endl;
}

Polyhedron * readPolyhedronVTK (istream &istr, bool perturb)
{
  Polyhedron *a = new Polyhedron;
  skipComments(istr);
  string dummy;
  do {
    istr >> dummy;
  }
  while (!(dummy == "POINTS" || dummy == "points"));
  int nv;
  istr >> nv >> dummy;
  acp::disable();
  for (int i = 0; i < nv; ++i) {
    double x, y, z;
    istr >> x >> y >> z;
    a->getVertex(new Point(x, y, z, perturb));
  }
  acp::enable();
  int nf;
  istr >> dummy >> nf >> dummy;
  for (int i = 0; i < nf; ++i) {
    int k, u, v, w;
    istr >> k >> u >> v >> w;
    if (k != 3) {
      cerr << "non triangle in readPolyhedronVTK" << endl;
      exit(0);
    }
    Vertex *uv = a->vertices[u], *vv = a->vertices[v], *wv = a->vertices[w];
    a->addTriangle(uv, vv, wv);
  }
  return a;
}

void skipComments (istream &istr)
{
  char s[10000];
  ws(istr);
  while (!istr.eof() && (istr.peek() == '#' || istr.peek() == '\n'))
    istr.getline(s, 1000);
}

void writePolyhedronVTK (Polyhedron *a, ostream &ostr)
{
  writePolyhedronVTK(a->faces, ostr);
}
  
void writePolyhedronVTK (const Faces &fa, ostream &ostr)
{
  vector<PV3<double>> pts;
  vector<unsigned int> data;
  VIMap vimap;
  for (Faces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
    Vertices ve = (*f)->getBoundary()->loop();
    data.push_back(ve.size());
    for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
      data.push_back(getPoint(vimap, pts, *v));
  }
  outputVTK(pts, data, 1, ostr);
}

int getPoint (VIMap &vimap, vector<PV3<double>> &pts, Vertex *v)
{
  VIMap::iterator iter = vimap.find(v);
  if (iter != vimap.end())
    return iter->second;
  int k = pts.size();
  pts.push_back(v->getP()->getApproxMid());
  vimap.insert(VIPair(v, k));
  return k;
}

void outputVTK (const vector<PV3<double>> &pts, const vector<unsigned int> &data,
		int ptype, ostream &ostr)
{
  int nv = pts.size();
  ostr << setprecision(20) << "# vtk DataFile Version 3.0" << endl
       << "vtk output" << endl << "ASCII" << endl
       << "DATASET POLYDATA" << endl 
       << "POINTS " << nv << " double" << endl;
  acp::disable();
  for (vector<PV3<double>>::const_iterator p = pts.begin(); p != pts.end(); ++p) 
    ostr << p->x << " " << p->y << " " << p->z << endl;
  acp::enable();
  int np = 0, i = 0;
  while (i < data.size()) {
    ++np;
    i += data[i] + 1;
  }
  ostr << endl;
  if (ptype == 1) ostr << "POLYGONS";
  else if (ptype == 2) ostr << "LINES";
  else ostr << "VERTICES";
  ostr << " " << np << " " << data.size() << endl;
  i = 0;
  while (i < data.size()) {
    ostr << data[i] << " ";
    for (int j = 0; j < data[i]; ++j)
      ostr << data[i+j+1] << " ";
    ostr << endl;
    i += data[i] + 1;
  }
}

Polyhedron * readPolyhedronOBJ (istream &istr, bool perturb)
{
  Polyhedron *a = new Polyhedron;
  Vertices ve;
  acp::disable();
  while (istr.peek() != EOF) {
    double x, y, z;
    if (readPointOBJ(istr, x, y, z))
      ve.push_back(a->getVertex(new Point(x, y, z, perturb)));
    else
      break;
  }
  acp::enable();
  int n = ve.size();
  while (istr.peek() != EOF) {
    int i, j, k;
    if (readTriangleOBJ(istr, n, i,  j, k))
      addTriangleOBJ(a, ve[i], ve[j], ve[k]);
  }
  return a;
}

bool readPointOBJ (istream &istr, double &x, double &y, double &z)
{
  while (!istr.eof()) {
    if (istr.peek() == 'f')
      return false;
    if (istr.peek() == 'v') {
      istr.get();
      if (istr.get() == ' ')
	break;
    }
    istr.ignore(1000, '\n');
  }
  if (istr.eof())
    return false;
  istr >> x >> y >> z;
  istr.ignore(1000, '\n');
  return true;
}

bool readTriangleOBJ (istream &istr, int n, int &i, int &j, int &k)
{
  while (!istr.eof() && istr.peek() != 'f')
    istr.ignore(1000, '\n');
  if (istr.eof())
    return false;
  istr.get();
  readIndexOBJ(istr, n, i);
  readIndexOBJ(istr, n, j);
  readIndexOBJ(istr, n, k);
  istr.ignore(1000, '\n');
  return true;
}

void readIndexOBJ (istream &istr, int n, int &i)
{
  istr >> i;
  if (i > 0)
    --i;
  else
    i = n - i;
  while (istr.peek() != ' ' && istr.peek() != '\n')
    istr.get();
}

void addTriangleOBJ (Polyhedron *a, Vertex *u, Vertex *v, Vertex *w)
{
  if (u == v || u == w || v == w || u->getP()->onLine(v->getP(), w->getP()))
    return;
  HEdges ue = u->outgoingHEdges();
  for (HEdges::iterator h = ue.begin(); h != ue.end(); ++h)
    if ((*h)->head() == w && (*h)->getNext()->head() == v)
      return;
  a->addTriangle(u, v, w);
}

// debug

void writePolyhedronVTK (const HFaces &fa, ostream &ostr)
{
  vector<PV3<double>> pts;
  vector<unsigned int> data;
  VIMap vimap;
  for (HFaces::const_iterator f = fa.begin(); f != fa.end(); ++f) {
    Face *ff = (*f)->getF();
    Vertices ve = ff->getBoundary()->loop();
    if (!(*f)->pos())
      reverse(ve.begin(), ve.end());
    data.push_back(ve.size());
    for (Vertices::iterator v = ve.begin(); v != ve.end(); ++v)
      data.push_back(getPoint(vimap, pts, *v));
  }
  outputVTK(pts, data, 1, ostr);
}

string outi (int i)
{
  return string("out") + std::to_string(i) + ".vtk";
}

void plines (const vector<vector<PV3<double>>> &lines, int i)
{
  vector<PV3<double>> pts;
  vector<unsigned int> data;
  int k = 0;
  for (vector<vector<PV3<double>>>::const_iterator l = lines.begin();
       l != lines.end(); ++l) {
    data.push_back(l->size());
    for (vector<PV3<double>>::const_iterator p = l->begin(); p != l->end();
	 ++p, ++k) {
      pts.push_back(*p);
      data.push_back(k);
    }
  }
  ofstream ostr(outi(i).c_str());
  outputVTK(pts, data, 2, ostr);
}

void pfaces (const Faces &faces, int i)
{
  ofstream ostr(outi(i).c_str());
  writePolyhedronVTK(faces, ostr);
}

void pfaces (const Faces &faces, int i, int is, int ie)
{
  Faces fa;
  for (int j = is; j < ie; ++j)
    fa.push_back(faces[j]);
  pfaces(fa, i);
}

void pfaces (Polyhedron *a, int i)
{
  pfaces(a->faces, i);
}

void pface (Face *f, int i)
{
  Faces fa;
  fa.push_back(f);
  pfaces(fa, i);
}

void pfaces (const set<Face *> &fs, int i)
{
  Faces fa(fs.begin(), fs.end());
  pfaces(fa, i);
}

void pfaces (const HFaces &hfaces, int i)
{
  ofstream ostr(outi(i).c_str());
  writePolyhedronVTK(hfaces, ostr);
}

void pcells (Polyhedron *a)
{
  if (a->cells.empty())
    a->formCells();
  for (int i = 1; i < a->cells.size(); ++i)
    pfaces(a->cells[i]->getShell(0)->getHFaces(), i);
}

void pfaces (const Faces &faces)
{
  writePolyhedronVTK(faces, cout);
}

void pfaces (const HFaces &hfaces)
{
  Faces faces;
  for (HFaces::const_iterator f = hfaces.begin(); f != hfaces.end(); ++f)
    faces.push_back((*f)->getF());
  writePolyhedronVTK(faces, cout);
}

void pface (Face *f)
{
  Faces faces;
  faces.push_back(f);
  writePolyhedronVTK(faces, cout);
}

void pedges (const Edges &edges, ostream &ostr)
{
  VIMap vidmap;
  vector<PV3<double>> pts;
  vector<unsigned int> data;
  for (Edges::const_iterator e = edges.begin(); e != edges.end(); ++e) {
    data.push_back(2);
    data.push_back(getPoint(vidmap, pts, (*e)->getT()));
    data.push_back(getPoint(vidmap, pts, (*e)->getH()));
  }
  outputVTK(pts, data, 2, ostr);
}

void pedges (const Edges &edges)
{
  pedges(edges, cout);
}

void pedges (const Edges &edges, int i)
{
  ofstream ostr(outi(i).c_str());
  pedges(edges, ostr);
}

void pedges (const HEdges &hedges)
{
  Edges edges;
  for (HEdges::const_iterator e = hedges.begin(); e != hedges.end(); ++e)
    edges.push_back((*e)->getE());
  pedges(edges);
}

void pedge (Edge *e)
{
  Edges ed;
  ed.push_back(e);
  pedges(ed);
}

void pedge (HEdge *e)
{
  pedge(e->getE());
}

void pedge (Edge *e, int i)
{
  Edges ed;
  ed.push_back(e);
  pedges(ed, i);
}

int getPoint (map<Point *, unsigned int> &pimap, vector<PV3<double>> &pts, Point *p)
{
  map<Point *, unsigned int>::iterator iter = pimap.find(p);
  if (iter != pimap.end())
    return iter->second;
  int k = pts.size();
  pts.push_back(p->getApproxMid());
  pimap.insert(pair<Point *, unsigned int>(p, k));
  return k;
}

void pedges (const SEdges &edges, ostream &ostr)
{
  map<Point *, unsigned int> pidmap;
  vector<PV3<double>> pts;
  vector<unsigned int> data;
  for (SEdges::const_iterator e = edges.begin(); e != edges.end(); ++e) {
    data.push_back(2);
    data.push_back(getPoint(pidmap, pts, (*e)->tail));
    data.push_back(getPoint(pidmap, pts, (*e)->head()));
  }
  outputVTK(pts, data, 2, ostr);
}

void pedges (const SEdges &edges)
{
  pedges(edges, cout);
}

void pvertices (const Vertices &ve, int i)
{
  VIMap vidmap;
  vector<PV3<double>> pts;
  vector<unsigned int> data;
  data.push_back(ve.size());
  for (Vertices::const_iterator v = ve.begin(); v != ve.end(); ++v)
    data.push_back(getPoint(vidmap, pts, *v));
  ofstream ostr(outi(i).c_str());
  outputVTK(pts, data, 3, ostr);
}
