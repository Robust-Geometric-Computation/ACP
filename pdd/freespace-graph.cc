#include "freespace-graph.h"
bool DEBUG = true;
const double SHRINKING_RADIUS = 1e-5;

map<string, PTRPolyhedron> polyhedra;

PTRPolyhedron loadPoly(const char * filename) {
  string name(filename);
  cout << "loadPoly " << name << endl;
  if (polyhedra.find(name) != polyhedra.end()) {
    PTRPolyhedron p = polyhedra[name];
#ifdef BLEEN
    cout << "cached" << endl;
    cout << "p->cells.size() " << p->cells.size() << endl;
    if (p->cells.size() > 0)
      cout << "p->cells[0]->inner.size() " << p->cells[0]->inner.size() << endl;
#endif
    return polyhedra[name];
  }

  PTRPolyhedron poly;
  ifstream infile (filename);
  if (infile.is_open()) {
    poly = readPolyhedron (infile, false);
    infile.close();
  } else {
    cout<<"could not read from file"<<endl;
    return NULL;
  }

  return poly;
}

PTRPolyhedron loadPolyVTK(const char * filename) {
  cout << "loadPolyVTK " << filename << endl;
  string name(filename);
  if (polyhedra.find(name) != polyhedra.end()) {
    PTRPolyhedron p = polyhedra[name];
#ifdef BLEEN
    cout << "cached" << endl;
    cout << "p->cells.size() " << p->cells.size() << endl;
    if (p->cells.size() > 0)
      cout << "p->cells[0]->inner.size() " << p->cells[0]->inner.size() << endl;
#endif
    return polyhedra[name];
  }
  PTRPolyhedron poly;
  ifstream infile (filename);
  if (infile.is_open()) {
    poly = readPolyhedronVTK (infile);
    infile.close();
  } else {
    cout<<"could not read from file"<<endl;
    return NULL;
  }

  return poly;
}


void savePoly(PTRPolyhedron p, const char * filename) {
  string name(filename);
  cout << "savePoly " << name << endl;
  polyhedra[name] = p;
#ifdef BLEEN
  string name2(filename);
  assert(polyhedra.find(name2) != polyhedra.end());
  assert(polyhedra[name2] == p);
  cout << "p->cells.size() " << p->cells.size() << endl;
  if (p->cells.size() > 0)
    cout << "p->cells[0]->inner.size() " << p->cells[0]->inner.size() << endl;
  string name3("frustum-tworooms-graph/0-0.tri");
  if (polyhedra.find(name3) != polyhedra.end()) {
    p = loadPoly("frustum-tworooms-graph/0-0.tri");
  }
#endif

#ifdef BLEEN
  ofstream out;
  out.open(filename);
  if (out.is_open()) {
    writePolyhedronVTK (p->faces, out);
    out.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
#endif

#ifdef BLEEN
  if (true)
    return;

  ofstream out;
  out.open(filename);
  if (out.is_open()) {
    writePolyhedron (p, out);
    out.close();
  } else {
    cout<<"could not write to file"<<endl;
  }
#endif
}


pair <PTR<Point>, PTR<Point> > nearestPointPair(PTRPolyhedron poly, int i, int j) {
  PTR<Point> a = NULL;
  PTR<Point> b = NULL;

  //populate HFace lists
  HFaces hfs1, hfs2;
  assert(poly->cells.at(i)->nShells() > 0);
  for (int k=0; k<poly->cells[i]->nShells(); k++) {
    Shell * s = poly->cells[i]->getShell(k);
    HFaces hfs = s->getHFaces();
    hfs1.insert(hfs1.end(), hfs.begin(), hfs.end());
  }
  for (int k=0; k<poly->cells[j]->nShells(); k++) {
    Shell * s = poly->cells[j]->getShell(k);
    HFaces hfs = s->getHFaces();
    hfs2.insert(hfs2.end(), hfs.begin(), hfs.end());
  }

  //populate vertex lists
  Vertices v1, v2;
  for (int k=0; k<hfs1.size(); k++) {
    HEdge * heFirst = hfs1[k]->getF()->getBoundary();
    HEdge * he = heFirst;
    v1.push_back(he->tail());
    do {
      if (std::find(v1.begin(), v1.end(), he->tail()) == v1.end())
        v1.push_back(he->tail());
      he = he->getNext();
    } while (he != heFirst);
  }
  for (int k=0; k<hfs2.size(); k++) {
    HEdge * heFirst = hfs2[k]->getF()->getBoundary();
    HEdge * he = heFirst;
    v2.push_back(he->tail());
    do {
      if (std::find(v2.begin(), v2.end(), he->tail()) == v2.end())
        v2.push_back(he->tail());
      he = he->getNext();
    } while (he != heFirst);
  }

  Vertex *av = 0;
  Vertex *bv = 0;

  //find closest pair
  for (int k=0; k<v1.size(); k++)
    for (int l=0; l<v2.size(); l++)
      if (av == 0 || CloserPair(v1[k]->getP(), v2[l]->getP(), av->getP(), bv->getP()) > 0) {
        av = v1[k];
        bv = v2[l];
      }

  a = poly->cells[i]->interiorPoint(av, 1e-6);
  b = poly->cells[j]->interiorPoint(bv, 1e-6);
  assert(a != 0);
  assert(b != 0);
  // a = av->getP();
  // b = bv->getP();

  return std::make_pair(a,b);
}

PTR<Point> nearestPoint(PTRPolyhedron poly, int cell_index, PTR<Point> p) {
  assert(!poly->cells.empty());
  
  int p_index = poly->containingCell(p);
  if (p_index == cell_index)
    return p;

  HFaces hfs1;
  for (int k=0; k<poly->cells[cell_index]->nShells(); k++) {
    Shell * s = poly->cells[cell_index]->getShell(k);
    HFaces hfs = s->getHFaces();
    hfs1.insert(hfs1.end(), hfs.begin(), hfs.end());
  }
  Vertices v;
  for (int k=0; k<hfs1.size(); k++) {
    HEdge * heFirst = hfs1[k]->getF()->getBoundary();
    HEdge * he = heFirst;
    do {
      if (std::find(v.begin(), v.end(), he->tail()) == v.end())
        v.push_back(he->tail());
      he = he->getNext();
    } while (he != heFirst);
  }

  Vertex *q = v[0];
  for (int i=1; i<v.size(); i++) {
    if (CloserPair(v[i]->getP(), p, q->getP(), p) > 0)
      q = v[i];
  }
  return poly->cells[cell_index]->interiorPoint(q, 1e-6);
  // return q->getP();
}

PTR<Point> nearestPoint(PTRPolyhedron poly, std::vector<int> & cells, PTR<Point> p) {
  PTR<Point> closest = 0;
  for (int i=0; i<cells.size(); i++) {
    PTR<Point> q = nearestPoint(poly, cells[i], p);
    if (closest == NULL || CloserPair(q, p, closest, p) > 0)
      closest = q;
  }
  assert(closest != NULL);
  return closest;
}

void FreeSpaceGraph::deepestPath(PTRFreeSpaceGraphNode n, PTR<Point> p, PTR<Point> q, std::vector<std::pair<PTRFreeSpaceGraphNode , PTR<Point> > > & path) {
  // invariant: (n,p) is already on the path, (n,q) will be added after return
  if (n->children.size() == 0) {
    return;
  }

  PTRFreeSpaceGraphNode P;
  PTRFreeSpaceGraphNode Q;
  PTR<Point> p1=0;
  PTR<Point> q1=0;
  std::string s = string(dir) + "/" + std::to_string(n->level+1) + "-" + std::to_string(n->blockspace_index) + ".tri";
  PTRPolyhedron childspace = loadPoly(s.c_str());
  assert(childspace != NULL);
  assert(childspace->cells.size() != 0);
  for (int i=0; i<childspace->cells.size(); i++)
    if (childspace->cells[i]->getWN() == 0) {
      PTR<Point> r = nearestPoint(childspace, i, p);
      if (p1 == 0 || r == p || CloserPair(r, p, p1, p) > 0) {
        p1 = r;
        P = graph[n->level+1][n->blockspace_index]->get(i);
      }
      r = nearestPoint(childspace, i, q);
      if (q1 == 0 || r == q || CloserPair(r, q, q1, q) > 0) {
        q1 = r;
        Q = graph[n->level+1][n->blockspace_index]->get(i);
      }
    }
  /*delete*/ childspace;
  assert(p != NULL);
  assert(q != NULL);
  assert(p1 != NULL);
  assert(q1 != NULL);
  
  if (P == Q) {
    if (p != p1) {
      path.push_back(make_pair(n, p1));
      path.push_back(make_pair(P, p1));
    } else {
      path.push_back(make_pair(P, p));
    }
    deepestPath(P, p1, q1, path);
    if (q != q1) {
      path.push_back(make_pair(Q, q1));
      path.push_back(make_pair(n, q1));
    } else {
      path.push_back(make_pair(Q, q));
    }
  } else {
    int i = std::find(P->siblings.begin(), P->siblings.end(), Q) - P->siblings.begin();
    assert(i < P->siblings.size());
    PTR<Point> a = P->siblingPoints[i].first;
    PTR<Point> b = P->siblingPoints[i].second;
    assert(a != NULL);
    assert(b != NULL);

    if (p != p1) {
      path.push_back(make_pair(n, p1));
      path.push_back(make_pair(P, p1));
    } else {
      path.push_back(make_pair(P, p));
    }
    deepestPath(P, p1, a, path);
    path.push_back(make_pair(P, a));
    path.push_back(make_pair(n, a));
    path.push_back(make_pair(n, b));
    path.push_back(make_pair(Q, b));
    deepestPath(Q, b, q1, path);
    if (q != q1) {
      path.push_back(make_pair(Q, q1));
      path.push_back(make_pair(n, q1));
    } else {
      path.push_back(make_pair(Q, q));
    }
  }
}

void FreeSpaceGraph::nodePointPath(std::vector<PTRFreeSpaceGraphNode> & nodes, PTR<Point> a, PTR<Point> b, std::vector<std::pair<PTRFreeSpaceGraphNode , PTR<Point> > > & path) {
  //working backwards, find the points we should go to in each node in the path
  std::vector<std::pair<PTRFreeSpaceGraphNode , PTR<Point> > > rev;
  rev.push_back(make_pair(nodes[nodes.size()-1], b));
  for (int i=nodes.size()-2; i>=0; i--) {
    PTRFreeSpaceGraphNode from = rev[rev.size()-1].first;
    PTRFreeSpaceGraphNode to = nodes[i];
    PTR<Point> p = rev[rev.size()-1].second;

    if (from->parent == to) {
      assert(p != NULL);
      rev.push_back(make_pair(to, p));
      continue;
    }

    if (to->parent == from) {
      PTRPolyhedron to_blockspace = loadPoly((dir + "/" + std::to_string(to->level) + "-" + std::to_string(to->blockspace_index) + ".tri").c_str());
      assert(to_blockspace != NULL);
      if (to_blockspace->containingCell(p) == to->cell_index) {
        rev.push_back(make_pair(to, p));
        assert(p != NULL);
      } else {
        PTR<Point> q = nearestPoint(to_blockspace, to->cell_index, p);
        // deepestPath(from, p, q, rev); //I don't think it makes sense to use deepest path when going from parent to child
        rev.push_back(make_pair(from, q));
        rev.push_back(make_pair(to, q));
        assert(q != NULL);
      }
      /*delete*/ to_blockspace;
      continue;
    } 

    //neighbor edge
    int k = std::find(to->neighbors.begin(), to->neighbors.end(), from) - to->neighbors.begin();
    assert(k < to->neighbors.size());
    assert(from->level == to->level);
    int bi = (to->blockspace_index == (from->blockspace_index+1)%blockspaces_per_level)? from->blockspace_index : to->blockspace_index;
    int bj = (bi == to->blockspace_index)? from->blockspace_index : to->blockspace_index;
    std::string s = std::string(dir) + "/" + std::to_string(to->level) + "-" + std::to_string(bi) + "-" + std::to_string(bj) + ".tri";
    PTRPolyhedron intersection = loadPoly(s.c_str() );
    int cell_index = intersection->containingCell(p);
    assert(intersection != NULL);
    if ( std::find(to->neighborIntersectionIndices[k].begin(), to->neighborIntersectionIndices[k].end(), cell_index) != to->neighborIntersectionIndices[k].end()) {
      rev.push_back(make_pair(to, p));
      assert(p != NULL);
    } else {
      PTR<Point> q = nearestPoint(intersection, to->neighborIntersectionIndices[k], p);
      assert(p != NULL);
      assert(q != NULL);
      deepestPath(from, p, q, rev);
      rev.push_back(make_pair(from, q));
      rev.push_back(make_pair(to, q));
      assert(q != NULL);
    }
    /*delete*/ intersection; 
  }
  rev.push_back(make_pair(nodes[0], a));
  assert(a != NULL);

  for (int i=rev.size()-1; i>=0; i--)
    path.push_back(rev[i]);
}

void FreeSpaceGraph::bfsPath(PTRFreeSpaceGraphNode start, PTRFreeSpaceGraphNode end, std::vector<PTRFreeSpaceGraphNode> & nodes) {
  //first use bfs to find the initial node path
  std::map<PTRFreeSpaceGraphNode, PTRFreeSpaceGraphNode> prev;
  std::queue<PTRFreeSpaceGraphNode> q;
  q.push(start);
  prev[start] = NULL;
  bool done = false;
  while (!(q.empty() || done)) { 
    PTRFreeSpaceGraphNode n = q.front(); q.pop();
    //enqueue neighbors
    for (int i=0; i<n->neighbors.size(); i++)
      if (n->neighbors[i]->enabled && prev.find(n->neighbors[i]) == prev.end()) {
        prev[n->neighbors[i]] = n;
        q.push(n->neighbors[i]);
        done = done || (n->neighbors[i] == end);
      }
    //enqueue children
    for (int i=0; i<n->children.size(); i++)
      if (n->children[i]->enabled && prev.find(n->children[i]) == prev.end()) {
        prev[n->children[i]] = n;
        q.push(n->children[i]);
        done = done || (n->children[i] == end);
      }
    //enqueue parent
    if (n->parent != NULL && n->parent->enabled && prev.find(n->parent) == prev.end()) {
      prev[n->parent] = n;
      q.push(n->parent);
      done = done || (n->parent == end);
    }
  }

  assert (prev.find(end) != prev.end());

  std::vector<PTRFreeSpaceGraphNode> rev;
  PTRFreeSpaceGraphNode n = end;
  while (n != NULL) {
    rev.push_back(n);
    n = prev[n];
  }
  for (int i=rev.size()-1; i>=0; i--)
    nodes.push_back(rev[i]);
}

bool FreeSpaceGraph::isConnected(PTRFreeSpaceGraphNode start, PTRFreeSpaceGraphNode end) {
  assert(start->enabled); assert(end->enabled);
  if (start == end)
    return true;

  std::queue<PTRFreeSpaceGraphNode> q;
  std::set<PTRFreeSpaceGraphNode> visited;
  q.push(start);
  while (!q.empty()) {
    PTRFreeSpaceGraphNode n = q.front(); q.pop();
    if (n == end)
      assert(false);
    
    //enqueue neighbors
    for (int i=0; i<n->neighbors.size(); i++) 
      if (n->neighbors[i]->enabled) {
        if (n->neighbors[i] == end)
          return true;
        if (visited.find(n->neighbors[i]) == visited.end()) {
          q.push(n->neighbors[i]); 
          visited.insert(n->neighbors[i]);
        }
      }
    
    //enqueue children
    for (int i=0; i<n->children.size(); i++) 
      if (n->children[i]->enabled) {
        if (n->children[i] == end)
          return true;
        if (visited.find(n->children[i]) == visited.end()) {
          q.push(n->children[i]); 
          visited.insert(n->children[i]);
        }
      }

    //enqueue parent
    if (n->parent == end)
      return true;
    if (n->parent != NULL && n->parent->enabled && visited.find(n->parent) == visited.end())
      q.push(n->parent);
  }
  return false;
}

double FreeSpaceGraph::angle(int blockspace_index) {
  return blockspace_index*theta + theta/2;
}

void FreeSpaceGraph::getPath(PTRNode start, PTRNode end, PTR<Point> a, PTR<Point> b, std::vector<std::pair<PTR<Point>, double > > & path) {
  if (!isConnected(start, end)) {
    cout<<"no path exists"<<endl;
    return;
  }

  //deepen start and end as much as possible
  for (int level= start->level+1; level <num_levels; level++) {
    PTRPolyhedron space = loadPoly((dir + "/" + std::to_string(level) + "-" + std::to_string(start->blockspace_index) + ".tri").c_str());
    int cell_index = space->containingCell(a);
    /*delete*/ space;
    if (graph[level][start->blockspace_index]->get(cell_index) != NULL)
      start = graph[level][start->blockspace_index]->get(cell_index);
    else
      break;
  }
  for (int level= end->level+1; level <num_levels; level++) {
    PTRPolyhedron space = loadPoly((dir + "/" + std::to_string(level) + "-" + std::to_string(end->blockspace_index) + ".tri").c_str());
    int cell_index = space->containingCell(b);
    /*delete*/ space;
    if (graph[level][end->blockspace_index]->get(cell_index) != NULL)
      end = graph[level][end->blockspace_index]->get(cell_index);
    else
      break;
  }

  for (int i=0; i<num_levels; i++)
    for (int j=0; j<blockspaces_per_level; j++)
      for (int k=0; k<graph[i][j]->nodes.size(); k++) 
        if (graph[i][j]->nodes[k] != start && graph[i][j]->nodes[k] != end) {
          graph[i][j]->nodes[k]->enabled = false;
          if (!isConnected(start, end))
            graph[i][j]->nodes[k]->enabled = true;
        }

  std::vector<PTRFreeSpaceGraphNode> nodes;
  std::vector<std::pair<PTRFreeSpaceGraphNode , PTR<Point> > > node_point_path;
  bfsPath(start, end, nodes);
  
  if (DEBUG) {
    for (int i=0; i<nodes.size(); i++)
      cout<<nodes[i]->level<<" "<<nodes[i]->blockspace_index<<" "<<nodes[i]->cell_index<<endl;
    cout<<endl;
  }

  nodePointPath(nodes, a, b, node_point_path);

  if (DEBUG) {
    for (int i=0; i< node_point_path.size(); i++) {
      PTRNode n = node_point_path[i].first;
      PTR<Point> p = node_point_path[i].second;
      cout<<n->level<<" "<<n->blockspace_index<<" "<<n->cell_index<<" "<<p<<endl;
    }
    cout<<endl;
  }

  path.push_back(make_pair(node_point_path[0].second, angle(node_point_path[0].first->blockspace_index)));
  if (DEBUG) cout<<node_point_path[0].second->getApprox().getX().mid()<<" "<<node_point_path[0].second->getApprox().getY().mid()<<" "<<node_point_path[0].second->getApprox().getZ().mid()<<" "<<angle(node_point_path[0].first->blockspace_index)<<" "<<node_point_path[0].first->level<<" "<<node_point_path[0].first->blockspace_index<<" "<<node_point_path[0].first->cell_index<<endl;
  for (int i=1; i<node_point_path.size(); i++) {
    PTRNode n1 = node_point_path[i-1].first;
    PTRNode n2 = node_point_path[i].first;
    PTR<Point> p1 = node_point_path[i-1].second;
    PTR<Point> p2 = node_point_path[i].second;
    if (n1 == n2 && p1 == p2)
      continue;
    if (n1 == n2) {
      //load blockspace associated with n1
      std::string s = std::string(dir) + "/" + std::to_string(n1->level) + "-" + std::to_string(n1->blockspace_index) + ".tri";
      PTRPolyhedron poly = loadPoly(s.c_str());
      //use path3d to find a path from p1 to p2
      Points subPath;
      findPath(poly, n1->cell_index, p1, p2, subPath);
      /*delete*/ poly;
      for (int i=1; i< subPath.size(); i++) {
        path.push_back(make_pair(subPath[i],angle(n1->blockspace_index)));
        if (DEBUG) cout<<subPath[i]->getApprox().getX().mid()<<" "<<subPath[i]->getApprox().getY().mid()<<" "<<subPath[i]->getApprox().getZ().mid()<<" "<<angle(n1->blockspace_index)<<" "<<n1->level<<" "<<n1->blockspace_index<<" "<<n1->cell_index<<((subPath.size() > 2)?" path3d" : "")<<endl;
      }
    } else {
      assert(p1 == p2);
      //angle change if needed
      if (n1->blockspace_index != n2->blockspace_index) {
        path.push_back(make_pair(p2,angle(n2->blockspace_index)));
      }
      if (DEBUG) cout<<p1->getApprox().getX().mid()<<" "<<p1->getApprox().getY().mid()<<" "<<p1->getApprox().getZ().mid()<<" "<<angle(n2->blockspace_index)<<" "<<n2->level<<" "<<n2->blockspace_index<<" "<<n2->cell_index<<endl;
    }
  }
}

FreeSpaceGraph::FreeSpaceGraph(std::vector<PTRPolyhedron> & close_blockspaces, std::vector<PTRPolyhedron> & rough_blockspaces, double theta, double clearance_unit, int num_levels, const char * dir) {
  double unit = clearance_unit;
  this->theta = theta;
  this->num_levels = num_levels;
  this->blockspaces_per_level = close_blockspaces.size();
  this->dir = dir;
  std::vector<PTRPolyhedron> blockspaces;
  blockspaces.insert(blockspaces.begin(), close_blockspaces.begin(), close_blockspaces.end());

  for (int i = 0; i<blockspaces.size(); i++)
    //eps removed: simplify(blockspaces[i], 1e-6);

  //create directory
  if (mkdir(dir, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) == -1 && errno != EEXIST)
    cout<<"could not create directry"<<endl;

  //initialize the graph
  for(int i=0; i<num_levels; i++) 
    graph.push_back(std::vector<PTRBlockSpaceNode>());
  for (int i=0; i<num_levels; i++) 
    for (int j=0; j< blockspaces.size(); j++)
      graph[i].push_back(new FreeSpaceGraph::BlockSpaceNode(i, j));

  std::vector<PTRPolyhedron> prev_blockspaces;
  time_t start,end, start2, end2;
  for (int level = 0; level<num_levels; level++) {
    cout<<"level "<<level<<endl;
    time(&start);
    double simplify_time = 0;

    if (level > 0) {
      if (level > 2)
        for (int i=0; i<prev_blockspaces.size(); i++)
          /*delete*/ prev_blockspaces[i];
      prev_blockspaces.clear();

      prev_blockspaces.insert(prev_blockspaces.begin(), blockspaces.begin(), blockspaces.end());
      blockspaces.clear();

      PTRPolyhedron ball = loadPolyVTK("ex/sphere.vtk");
      PTRPolyhedron unit_ball = ball->scale(unit);
      unit = unit * 2;
      for (int i=0; i<prev_blockspaces.size(); i++) {
        blockspaces.push_back(minkowskiSum(rough_blockspaces[i], unit_ball));
        time(&start2);
        //eps removed: simplify(blockspaces[i], 1e-6);
        time(&end2);
        simplify_time += difftime(end2, start2);

      }
      
      /*delete*/ unit_ball;
    }

    cout<<"creating nodes"<<endl;
    //create nodes
    for (int i=0; i<blockspaces.size(); i++) {
      blockspaces[i]->computeWindingNumbers();
      for (int j=0; j<blockspaces[i]->cells.size(); j++) {
	assert(blockspaces[i]->cells.at(j)->nShells() > 0);
        if (blockspaces[i]->cells[j]->getWN() == 0) {
          if (level > 0) {
            //find the parent of this cell
            // PTR<Point> p = pointInCell(blockspaces[i], j);
            cout<<i<<" "<<j<<endl;
            int comp;
            if (j == 0)
              comp = 0;
            else {
              PTR<Point> p = blockspaces[i]->cells[j]->interiorPoint();
              comp = prev_blockspaces[i]->containingCell(p);
            }
            PTRFreeSpaceGraphNode parent = graph[level-1][i]->get(comp);
            assert(parent != NULL);
            PTRFreeSpaceGraphNode newNode = new FreeSpaceGraph::Node(level, i, j, parent);
            graph[level][i]->nodes.push_back(newNode);
            parent->children.push_back(newNode);
          } else 
            graph[level][i]->nodes.push_back(new FreeSpaceGraph::Node(level, i, j, NULL));
        }
      }
    }

    PTRPolyhedron tmp = loadPolyVTK("ex/sphere.vtk");
    PTRPolyhedron shrinking_sphere = tmp->scale(SHRINKING_RADIUS);

    cout<<"creating edges"<<endl;
    //create neighbor edges
    for (int i=0; i<blockspaces.size(); i++) {
      int j = (i+1)%blockspaces.size();
      cout << "computing block_union" << endl;
      PTRPolyhedron block_union = blockspaces[i]->boolean(blockspaces[j], Union);

#ifdef SHRINK_BLOCK_UNION
      cout << "shrinking block_union" << endl;
      tmp = minkowskiSum(block_union, shrinking_sphere);
      /*delete*/ block_union; block_union = tmp;

      cout <<"start simplify block_union"<<endl;
      time(&start2);
      //eps removed: simplify(block_union, 1e-6);
      time(&end2);
      simplify_time += difftime(end2, start2);
      cout <<"end simplify block_union"<<endl;
#endif
      cout << "computing winding numbers" << endl;
      block_union->computeWindingNumbers();

      cout << "saving block_union" <<endl;
      std::string s = std::string(dir) + "/" + std::to_string(level) + "-" + std::to_string(i) + "-" + std::to_string(j) + ".tri";
      savePoly(block_union, s.c_str());

      for (int k=0; k< block_union->cells.size(); k++)
        if (block_union->cells[k]->getWN() == 0) {
          cout<<i<<" "<<j<<" "<<k<<endl;
          int ci, cj;
          if (k == 0) {
            ci = 0; cj = 0;
          } else {
            PTR<Point> p = block_union->cells[k]->interiorPoint();
            ci = blockspaces[i]->containingCell(p);
            cj = blockspaces[j]->containingCell(p); 
          }
          PTRNode ni = graph[level][i]->get(ci);
          PTRNode nj = graph[level][j]->get(cj);
	  if (nj == 0) // vjm
	    cout << "cj " << cj << endl;
          assert(nj != NULL); assert(ni != NULL);
          int pos_i = std::find(ni->neighbors.begin(), ni->neighbors.end(), nj) - ni->neighbors.begin();
          if (pos_i == ni->neighbors.size()) {
            ni->neighbors.push_back(nj);
            ni->neighborIntersectionIndices.push_back(std::vector<int>());
          }
          int pos_j = std::find(nj->neighbors.begin(), nj->neighbors.end(), ni) - nj->neighbors.begin();
          if (pos_j == nj->neighbors.size()) {
            nj->neighbors.push_back(ni);
            nj->neighborIntersectionIndices.push_back(std::vector<int>());
          }
          ni->neighborIntersectionIndices[pos_i].push_back(k);
          nj->neighborIntersectionIndices[pos_j].push_back(k);
        }
      /*delete*/ block_union;  
    }

    cout<<"finding siblings"<<endl;
    for (int i=0; i<blockspaces.size(); i++) {
      int nCells = graph[level][i]->nodes.size();
      if (nCells > 1)
        for (int s1 = 0; s1 <nCells-1; s1++)
          for (int s2=s1+1; s2<nCells; s2++) {
            PTRFreeSpaceGraphNode sibling1 = graph[level][i]->nodes[s1];
            PTRFreeSpaceGraphNode sibling2 = graph[level][i]->nodes[s2];
            sibling1->siblings.push_back(sibling2);
            sibling2->siblings.push_back(sibling1);
            std::pair< PTR<Point>, PTR<Point> > ab = nearestPointPair(blockspaces[i], sibling1->cell_index, sibling2->cell_index);
            sibling1->siblingPoints.push_back(ab);
            sibling2->siblingPoints.push_back(std::make_pair(ab.second,ab.first));
          }
    }

    for (int i=0; i<blockspaces.size(); i++) {
      std::string s = std::string(dir) + "/" + std::to_string(level) + "-" + std::to_string(i) + ".tri";
      cout<<"saving blockspace "<<i<<endl;
      savePoly(blockspaces[i], s.c_str());
    }
    time(&end);
    cout << difftime(end, start)/60.0 << " minutes elapsed on level " << level << endl;
    cout << simplify_time/60.0 << " minutes spent simplifying on level " << level << endl;
  }

  for(int i=0; i<prev_blockspaces.size(); i++)
    /*delete*/ prev_blockspaces[i];
  prev_blockspaces.clear();

  if (num_levels > 1)
    for (int i=0; i<blockspaces.size(); i++)
      /*delete*/ blockspaces[i];

  blockspaces.clear();

  cout<<"saving"<<endl;

  std::string s = string(dir) + "/graph.txt";
  ofstream out;
  out.open(s.c_str());
  out << setprecision(20);
  out << num_levels << " " << blockspaces_per_level << " " << theta << " " << clearance_unit <<endl;
  if (out.is_open()) {
    for (int i=0; i<num_levels; i++)
      for (int j=0; j<blockspaces_per_level; j++)
        for (int k=0; k<graph[i][j]->nodes.size(); k++) {
          PTRFreeSpaceGraphNode node = graph[i][j]->nodes[k];
          out << node->level << " " << node->blockspace_index << " " << node->cell_index << ",";
          if (node->parent != NULL) out << " " << node->parent->level << " " << node->parent->blockspace_index << " " << node->parent->cell_index << ",";
          else out << " NULL,";
          for (int l=0; l<node->children.size(); l++)
            out<< " " << node->children[l]->level << " " << node->children[l]->blockspace_index << " " << node->children[l]->cell_index;
          out<<",";
          for (int l=0; l<node->neighbors.size(); l++) {
            out<< " " << node->neighbors[l]->level << " " << node->neighbors[l]->blockspace_index << " " << node->neighbors[l]->cell_index << " [";
            for (int m=0; m<node->neighborIntersectionIndices[l].size(); m++)
              out<<(m>0?" ":"")<<node->neighborIntersectionIndices[l][m];
            out<<"]";
          }
          if (node->siblings.size()>0) {
          out<<",";
            for (int l=0; l<node->siblings.size(); l++)
              out<< " " << node->siblings[l]->level << " " << node->siblings[l]->blockspace_index << " " << node->siblings[l]->cell_index;
            out<<",";
            for (int l=0; l<node->siblingPoints.size(); l++)
              out << " " << node->siblingPoints[l].first->getApprox().getX().mid()
                  << " " << node->siblingPoints[l].first->getApprox().getY().mid()
                  << " " << node->siblingPoints[l].first->getApprox().getZ().mid()
                  << " " << node->siblingPoints[l].second->getApprox().getX().mid()
                  << " " << node->siblingPoints[l].second->getApprox().getY().mid()
                  << " " << node->siblingPoints[l].second->getApprox().getZ().mid();
          }
          out<<endl;
        }
    out.close();
  }
}

FreeSpaceGraph::FreeSpaceGraph(const char * dir) {
  this->dir = dir;
  std::string s = string(dir) + "/graph.txt";
  ifstream infile (s.c_str());
  double clearance_unit;
  if (infile.is_open()) {
    infile >> num_levels >> blockspaces_per_level >> theta >> clearance_unit;

    //initialize graph
    for(int i=0; i<num_levels; i++) 
      graph.push_back(std::vector<PTRBlockSpaceNode>());
    for (int i=0; i<num_levels; i++)
      for (int j=0; j< blockspaces_per_level; j++)
        graph[i].push_back(new FreeSpaceGraph::BlockSpaceNode(i, j));

    std::string line;
    std::getline(infile, line); //eat new line
    while (std::getline(infile, line)) {
      istringstream ss(line);
      string id, parent, children, neighbors, siblings, siblingPoints;
      PTRFreeSpaceGraphNode n;
      if (getline(ss, id, ',')) {
        istringstream ss2(id);
        int i,j,k;
        if (ss2 >> i >> j >> k)
          n = graph[i][j]->getOrCreate(k);
      } else { cout << "cannot read line ("<<line<<")"<<endl; return; }
      if (getline(ss, parent, ',') && parent.find("NULL") == string::npos) {
        istringstream ss2(parent);
        int i,j,k;
        if (ss2 >> i >> j >> k)
          n->parent = graph[i][j]->getOrCreate(k);
      }
      if (getline(ss, children, ',')) {
        istringstream ss2(children);
        int i,j,k;
        while (ss2 >> i >> j >> k)
          n->children.push_back(graph[i][j]->getOrCreate(k));
      }
      if (getline(ss, neighbors, ',')) {
        istringstream ss2(neighbors);
        std::string s1, s2;
        while(getline(ss2, s1, '[')) {
          getline(ss2, s2, ']');
          istringstream ss3(s1);
          int i, j, k;
          ss3 >> i >> j >> k;
          PTRNode neighbor = graph[i][j]->getOrCreate(k);
          n->neighbors.push_back(neighbor);
          n->neighborIntersectionIndices.push_back(std::vector<int>());
          istringstream ss4(s2);
          while (ss4 >> i)
            n->neighborIntersectionIndices[n->neighbors.size()-1].push_back(i);
        }
      }
      if (getline(ss, siblings, ',')) {
        istringstream ss2(siblings);
        int i,j,k;
        while (ss2 >> i >> j >> k)
          n->siblings.push_back(graph[i][j]->getOrCreate(k));
      }
      if (getline(ss, siblingPoints, ',')) {
        istringstream ss2(siblingPoints);
        double x1, y1, z1, x2, y2, z2;
        while (ss2 >> x1 >> y1 >> z1 >> x2 >> y2 >> z2) {
          PTR<Point> a = new Point(x1, y1, z1);
          PTR<Point> b = new Point(x2, y2, z2);
          n->siblingPoints.push_back(std::make_pair(a,b));
        }
      }
    }
    infile.close();
  } else {
    cout<<"could not read from file"<<endl;
  }
}
