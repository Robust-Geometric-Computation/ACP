#include "acp/delaunay/delaunay.h"

int Triangle::n;

map<PTR<Object<PV2>>, int> Triangle::rootIndex;

// bounding triangle of points.......
// change to two triangle bounding box
Triangle::Triangle(vector<PTR<Object<PV2>>> points) {
  PTR<Object<PV2>> pMinX = points[0];
  for (int i = 1; i < points.size(); i++)
    if (XOrder(pMinX, points[i]) == -1) pMinX = points[i];
  double minX = pMinX->getApprox(1.0).x.mid();
  PTR<Object<PV2>> pMinY = points[0];
  for (int i = 1; i < points.size(); i++)
    if (YOrder(pMinY, points[i]) == -1) pMinY = points[i];
  double minY = pMinY->getApprox(1.0).y.mid();
  PTR<Object<PV2>> pMaxX = points[0];
  for (int i = 1; i < points.size(); i++)
    if (XOrder(pMaxX, points[i]) == 1) pMaxX = points[i];
  double maxX = pMaxX->getApprox(1.0).x.mid();
  PTR<Object<PV2>> pMaxY = points[0];
  for (int i = 1; i < points.size(); i++)
    if (YOrder(pMaxY, points[i]) == 1) pMaxY = points[i];
  double maxY = pMaxY->getApprox(1.0).y.mid();

  vert[0] = new InputPoint(minX - 0.1, minY - 0.1);
  vert[1] = new InputPoint(maxX * 2 - minX + 0.2, minY - 0.1);
  vert[2] = new InputPoint(minX - 0.1, 2 * maxY - minY + 0.2);

  next[0] = next[1] = next[2] = 0;
  chld[0] = chld[1] = chld[2] = 0;
  index = n++;
}

// Locate the leaf triangle (no children) under this one which contains p.
// Pre: this triangle contains p.
PTR<Triangle> Triangle::locate(PTR<Object<PV2>> p) {
  if (chld[0] == 0) return this;

  for (int i = 0; i < 3; i++)
    if (chld[i] != 0 && chld[i]->contains(p)) return chld[i]->locate(p);

  return 0;
}

// Split this triangle into three child triangles at point p.
// Set child triangles.
void Triangle::split(PTR<Object<PV2>> p) {
  // Each child starts as a copy of this
  // Set vert[i] of child i to p
  // Update the neighbor across the opposite edge to point to the
  // child instead of this.
  for (int i = 0; i < 3; i++) {
    chld[i] = new Triangle(this);
    chld[i]->vert[i] = p;
    chld[i]->updateNeighbor(i, this);
  }

  // Now next[i] for child i is set correctly, but the other two next
  // fields should be set to the other two children instead of the
  // neighbors of this.
  for (int i = 0; i < 3; i++) {
    int j = (i + 1) % 3;
    int k = (i + 2) % 3;
    chld[i]->next[j] = chld[j];
    chld[i]->next[k] = chld[k];
  }

  // Debugging code.  Only a partial check of correctness, though.
  for (int i = 0; i < 3; i++) chld[i]->checkNeighbors();
}

// Find the index of that among the next (neighboring) triangles of this.
int Triangle::find(PTR<Triangle> that) {
  for (int i = 0; i < 3; i++)
    if (next[i] == that) return i;
  assert(false);
  return -1;
}

// Debugging. Check to make sure that all the neighbors of this
// triangle think that this triangle is a neighbor.
void Triangle::checkNeighbors() {
  // cout << "checking neighbors" << endl;
  for (int i = 0; i < 3; i++)
    if (next[i] != 0) next[i]->find(this);
}

// This triangle replaces oldT along the edge shared with neighbor next[i].
// Update the neighbor to point to this triangle instead.
void Triangle::updateNeighbor(int i, PTR<Triangle> oldT) {
  if (next[i] != 0) next[i]->next[next[i]->find(oldT)] = this;
}

// Flip this triangle with next[i].  They share the common edge
// opposite vert[i].
void Triangle::flip(int i) {
  // EXERCISE

  // Determine i, j, k in order for this triangle and i2, j2, k2 for
  // the neighboring triangle where its i2 is across from i.
  int j = (i + 1) % 3;
  int k = (i + 2) % 3;

  int i2 = next[i]->find(this);
  int j2 = (i2 + 1) % 3;
  int k2 = (i2 + 2) % 3;

  if (vert[k] != next[i]->vert[j2]) {
    int temp = j2;
    j2 = k2;
    k2 = temp;
  }

  // Child 0 is a copy of this triangle.
  // Change its vert[k].
  // Set its next[i].
  // Update that neighbor.
  // It has the correct next[i], but that neighbor doesn't know it.
  chld[0] = new Triangle(this);
  chld[0]->vert[j] = next[i]->vert[i2];
  chld[0]->next[i] = next[i]->next[k2];
  chld[0]->updateNeighbor(i, next[i]);
  chld[0]->updateNeighbor(j, this);

  // Similarly child 1.
  chld[1] = new Triangle(this);
  chld[1]->vert[k] = next[i]->vert[i2];
  chld[1]->next[i] = next[i]->next[j2];
  chld[1]->updateNeighbor(i, next[i]);
  chld[1]->updateNeighbor(k, this);

  // Child 0 and child 1 are next to each other.
  chld[0]->next[k] = chld[1];
  chld[1]->next[j] = chld[0];

  // next[i] has the same two children.
  next[i]->chld[0] = chld[0];
  next[i]->chld[1] = chld[1];

  // debugging
  chld[0]->checkNeighbors();
  chld[1]->checkNeighbors();
}

// check if all
bool Triangle::rootTriangle() {
  return rootIndex[vert[0]] != -1 && rootIndex[vert[0]] == rootIndex[vert[1]] &&
         rootIndex[vert[1]] == rootIndex[vert[2]];
}

bool Triangle::sameRoot(PTR<Object<PV2>> d) {
  if (!rootTriangle()) return false;

  return rootIndex[vert[0]] == rootIndex[d];
}

// Check if this triangle needs to be flipped with next[i].  They
// share the common edge opposite vert[i].
// If so, do the flip and then recursively check the two new
// triangles if they need to be flipped with the neighbor opposite
// this same vertex.
void Triangle::checkForFlip(int i) {
  if (next[i] == 0) return;

  // EXERCISE 3

  // Set j,k and i2 as in flip.
  // Set a, b, c as verts i, j, k and d as i2 in next[i].
  int j = (i + 1) % 3;
  int k = (i + 2) % 3;

  int i2 = next[i]->find(this);
  int j2 = (i2 + 1) % 3;
  int k2 = (i2 + 2) % 3;

  if (vert[k] != next[i]->vert[j2]) {
    int temp = j2;
    j2 = k2;
    k2 = temp;
  }

  PTR<Object<PV2>> a = vert[i];
  PTR<Object<PV2>> b = vert[j];
  PTR<Object<PV2>> c = vert[k];
  PTR<Object<PV2>> d = next[i]->vert[i2];

  // do not flip with triangles on the same root as you
  // avoids identity with Circ::Dist
  if (sameRoot(d)) {
    printf("(%g, %g) [%d]\n", a->getApprox(1.0).x.mid(),
           a->getApprox(1.0).y.mid(), rootIndex[a]);
    printf("(%g, %g) [%d]\n", b->getApprox(1.0).x.mid(),
           b->getApprox(1.0).y.mid(), rootIndex[b]);
    printf("(%g, %g) [%d]\n", c->getApprox(1.0).x.mid(),
           c->getApprox(1.0).y.mid(), rootIndex[c]);
    printf("(%g, %g) [%d]\n", d->getApprox(1.0).x.mid(),
           d->getApprox(1.0).y.mid(), rootIndex[d]);
    printf("\n");
    return;
  }

  // If ad is not inside abc and dcb, forget it.
  // If d is not inside circle abc, forget it.
  PTR<Object<Circle>> circle = new Circle3(a, b, c);
  if (Orient2D(a, b, d) == -1 && Orient2D(a, d, c) == -1 &&
      CircleDist(circle, d) < 0) {
    // Do the flip
    this->flip(i);

    // Check the two new children for flip.
    chld[0]->checkForFlip(i);
    chld[1]->checkForFlip(i);
  }
}

vector<PTR<Triangle>> delaunay(vector<PTR<Object<PV2>>>& bb,
                               vector<PTR<Object<PV2>>>& points) {
  if (points.size() == 0) return vector<PTR<Triangle>>();

  vector<PTR<Triangle>> roots;

  roots.push_back(new Triangle(bb[0], bb[1], bb[2]));
  roots.push_back(new Triangle(bb[0], bb[2], bb[3]));

  roots[0]->next[1] = roots[1];
  roots[1]->next[2] = roots[0];

  for (int i = 0; i < points.size(); i++) {
    PTR<Object<PV2>> p = points[i];

    if (roots[0]->contains(p)) {
      PTR<Triangle> t1 = roots[0]->locate(p);
      assert(t1);

      t1->split(p);
      for (int j = 0; j < 3; j++) t1->chld[j]->checkForFlip(j);

    } else if (roots[1]->contains(p)) {
      PTR<Triangle> t2 = roots[1]->locate(p);
      assert(t2);

      t2->split(p);
      for (int j = 0; j < 3; j++) t2->chld[j]->checkForFlip(j);

    } else {
      fprintf(stderr, "Point not in triangulation\n");
      exit(1);
    }
  }

  return roots;
}
