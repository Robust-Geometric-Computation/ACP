#include "path3d.h"
#include "mink.h"
#include "io.h"
#include <sys/stat.h>
#include <errno.h>
#include <string>
#include <iostream>

typedef PTR<Polyhedron> PTRPolyhedron;

PTRPolyhedron loadPoly(const char * filename);
PTRPolyhedron loadPolyVTK(const char * filename);
void savePoly(PTRPolyhedron p, char * filename);

class FreeSpaceGraph {
 public: 
  double theta;

  class Node;
  typedef Node* PTRNode;

  class Node : public RefCnt {
   public:
    int level;
    int blockspace_index;
    int cell_index;
    bool enabled;
    PTRNode parent;  
    std::vector<PTRNode> children; 
    std::vector<PTRNode> neighbors;
    std::vector<std::vector<int> > neighborIntersectionIndices; 
    std::vector<PTRNode> siblings;
    std::vector< std::pair< PTR<Point>, PTR<Point> > > siblingPoints;
    Node(int level, int blockspace_index, int cell_index, PTRNode parent)
     : level(level), blockspace_index(blockspace_index), cell_index(cell_index), parent(parent), enabled(true) {}
  };

  class BlockSpaceNode;
  typedef BlockSpaceNode* PTRBlockSpaceNode;

  class BlockSpaceNode {
   public:
    int level;
    int blockspace_index;
    std::vector<PTRNode> nodes;
    BlockSpaceNode(int level, int blockspace_index) : level(level), blockspace_index(blockspace_index) {}
    PTRNode get(int index) {
      for (int i=0; i<nodes.size(); i++)
        if (nodes[i]->cell_index == index)
          return nodes[i];
      return NULL;
    }
    PTRNode getOrCreate(int index) {
      PTRNode n = get(index);
      if (n == NULL) {
        n = new Node(level, blockspace_index, index, NULL);
        nodes.push_back(n);
      }
      return n;
    }
  };

  int num_levels;
  int blockspaces_per_level;
  std::string dir;
  std::vector<std::vector<PTRBlockSpaceNode> > graph;

  typedef FreeSpaceGraph::Node* PTRFreeSpaceGraphNode;

  void deepestPath(PTRFreeSpaceGraphNode n, PTR<Point> p, PTR<Point> q, std::vector<std::pair<PTRFreeSpaceGraphNode , PTR<Point> > > & path);
  void nodePointPath(std::vector<PTRFreeSpaceGraphNode> & nodes, PTR<Point> a, PTR<Point> b, std::vector<std::pair<PTRFreeSpaceGraphNode , PTR<Point> > > & path);
  void bfsPath(PTRFreeSpaceGraphNode start, PTRFreeSpaceGraphNode end, std::vector<PTRFreeSpaceGraphNode> & nodes);
  bool isConnected(PTRFreeSpaceGraphNode start, PTRFreeSpaceGraphNode end);
  void getPath(PTRNode start, PTRNode end, PTR<Point> a, PTR<Point> b, std::vector<std::pair<PTR<Point>, double > > & path);
  double angle(int blockspace_index);

  FreeSpaceGraph(std::vector<PTRPolyhedron> & close_blockspaces, std::vector<PTRPolyhedron> & rough_blockspaces, double theta, double clearance_unit, int num_levels, const char * dir);
  FreeSpaceGraph(const char * dir);
};

typedef FreeSpaceGraph::Node* PTRFreeSpaceGraphNode;
