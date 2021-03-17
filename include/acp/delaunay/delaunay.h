#include <map>
#include "acp/circle/circle.h"
#include "acp/dcel/dcel.h"

using namespace std;

class Triangle : public RefCnt {
 public:
  static int n;
  int index;
  PTR<Object<PV2>> vert[3];
  PTR<Triangle> next[3];
  PTR<Triangle> chld[3];

  PTR<Edge> edges[3];

  static map<PTR<Object<PV2>>, int> rootIndex;

  // Create a triangle from three points.
  Triangle(PTR<Object<PV2>> a, PTR<Object<PV2>> b, PTR<Object<PV2>> c) {
    vert[0] = a;
    vert[1] = b;
    vert[2] = c;
    next[0] = next[1] = next[2] = 0;
    chld[0] = chld[1] = chld[2] = 0;
    index = n++;
    edges[0] = NULL;
    edges[1] = NULL;
    edges[2] = NULL;
  }

  // Create a large triangle that contains a set of points in its interior.
  Triangle(vector<PTR<Object<PV2>>> points);

  // Copy vertices and next triangles of that.  Children are set to null.
  Triangle(PTR<Triangle> that) {
    *this = *that;
    chld[0] = chld[1] = chld[2] = 0;
    index = n++;
  }

  // Is point one of the vertices of this triangle?
  bool hasVertex(PTR<Object<PV2>> point) {
    return point == vert[0] || point == vert[1] || point == vert[2];
  }

  // Is point contained inside this triangle?
  bool contains(PTR<Object<PV2>> point) {
    return (Orient2D(point, vert[0], vert[1]) == -1 &&
            Orient2D(point, vert[1], vert[2]) == -1 &&
            Orient2D(point, vert[2], vert[0]) == -1);
  }

  // Locate the leaf triangle (no children) under this one which contains p.
  // Pre: this triangle contains p.
  PTR<Triangle> locate(PTR<Object<PV2>> p);

  // Split this triangle into three child triangles at point p.
  // Set child triangles.
  void split(PTR<Object<PV2>> p);

  // Find the index of that among the next (neighboring) triangles of this.
  int find(PTR<Triangle> that);

  // Debugging. Check to make sure that all the neighbors of this
  // triangle think that this triangle is a neighbor.
  void checkNeighbors();

  // This triangle replaces oldT along the edge shared with neighbor next[i].
  // Update the neighbor to point to this triangle instead.
  void updateNeighbor(int i, PTR<Triangle> oldT);

  // Flip this triangle with next[i].  They share the common edge
  // opposite vert[i].
  void flip(int i);

  // Check if this triangle needs to be flipped with next[i].  They
  // share the common edge opposite vert[i].
  // If so, do the flip and then recursively check the two new
  // triangles if they need to be flipped with the neighbor opposite
  // this same vertex.
  void checkForFlip(int i);

  bool sameRoot(PTR<Object<PV2>>);

  bool rootTriangle();
};

// Compute the Delaunay triangulation as a tree rooted at the bounding
// triangle of the points.  The actual triangulation is the leaves of
// the tree.
vector<PTR<Triangle>> delaunay(vector<PTR<Object<PV2>>>& bb,
                               vector<PTR<Object<PV2>>>& points);
