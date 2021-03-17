#ifndef INTERSECT_H
#define INTERSECT_H
#include <stdio.h>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class vertex {
 public:
  float x, y, z;
  vertex() {}
  vertex(float x0, float y0, float z0) {
    x = x0;
    y = y0;
    z = z0;
  }
  bool operator==(const vertex& v) const {
    if (fabs(x - v.x) < 1e-5 && fabs(y - v.y) < 1e-5 && fabs(z - v.z) < 1e-5)
      return true;

    return false;
  }
};

class Triangle {
 public:
  vertex v1;
  vertex v2;
  vertex v3;

  bool operator==(const Triangle& tri) const {
    if (v1 == tri.v1 && v2 == tri.v2 && v3 == tri.v3) return true;

    return false;
  }
};

static void VSub(float res[3], float a[3], float b[3]) {
  res[0] = a[0] - b[0];
  res[1] = a[1] - b[1];
  res[2] = a[2] - b[2];
}

static void Vcross(float res[3], float a[3], float b[3]) {
  res[0] = a[1] * b[2] - a[2] * b[1];
  res[1] = a[0] * b[2] - a[2] * b[0];
  res[2] = a[1] * b[0] - a[0] * b[1];
}

static float Vdot(float a[3], float b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static bool EdgeIntersection(vertex& v1, vertex& v2, Triangle& tri,
                             vertex& intersect_v) {
  float sp[3] = {v1.x, v1.y, v1.z};
  float sq[3] = {v2.x, v2.y, v2.z};

  float a[3] = {tri.v1.x, tri.v1.y, tri.v1.z};
  float b[3] = {tri.v2.x, tri.v2.y, tri.v2.z};
  float c[3] = {tri.v3.x, tri.v3.y, tri.v3.z};

  float v, w;
  float ab[3], ac[3], qp[3], ap[3], norm[3], e[3];
  VSub(ab, b, a);
  VSub(ac, c, a);
  VSub(qp, sp, sq);
  // Compute triangle normal. Can be precalculated or cached if
  // intersecting multiple segments against the same triangle
  Vcross(norm, ab, ac);
  // Compute denominator d. If d <= 0, segment is parallel to or points
  // away from triangle, so exit early

  float d = Vdot(qp, norm);
  if (d <= 0.0f) return false;
  // Compute intersection t value of pq with plane of triangle. A ray
  // intersects if 0 <= t. Segment intersects if 0 <= t <= 1. Delay
  // dividing by d until intersection has been found to pierce triangle
  VSub(ap, sp, a);
  float t = Vdot(ap, norm);
  if (t < 1e-9) return false;
  if (t > d) return false;
  // For segment; exclude this code line for a ray test
  // Compute barycentric coordinate components and test if within bounds
  Vcross(e, qp, ap);
  v = Vdot(ac, e);
  if (v < 1e-9 || v > d) return false;
  w = -Vdot(ab, e);
  if (w < 1e-10 || v + w > d) return false;
  // Segment/ray intersects triangle. Perform delayed division
  t /= d;

  float fIntersect[3];
  fIntersect[0] = sp[0] - qp[0] * t;
  fIntersect[1] = sp[1] - qp[1] * t;
  fIntersect[2] = sp[2] - qp[2] * t;
  intersect_v.x = fIntersect[0];
  intersect_v.y = fIntersect[1];
  intersect_v.z = fIntersect[2];
  return true;
}

static int TriangleIntersection(Triangle& tri1, Triangle& tri2,
                                vertex& intersect_v) {
  int res = 0;
  res = EdgeIntersection(tri1.v1, tri1.v2, tri2, intersect_v);
  if (res > 0) return res;
  res = EdgeIntersection(tri1.v1, tri1.v3, tri2, intersect_v);
  if (res > 0) return res;
  res = EdgeIntersection(tri1.v2, tri1.v3, tri2, intersect_v);
  if (res > 0) return res;

  res = EdgeIntersection(tri2.v1, tri2.v2, tri1, intersect_v);
  if (res > 0) return res;
  res = EdgeIntersection(tri2.v1, tri2.v3, tri1, intersect_v);
  if (res > 0) return res;
  res = EdgeIntersection(tri2.v2, tri2.v3, tri1, intersect_v);
  if (res > 0) return res;

  return 0;
}

static void SurfaceIntersection(std::vector<Triangle>& surface1,
                                std::vector<Triangle>& surface2,
                                std::vector<Triangle>& vecIntersectTri,
                                std::vector<vertex>& vecIntersPoints) {
  int num1 = surface1.size();
  int num2 = surface2.size();

  for (int i = 0; i < num1; ++i) {
    for (int j = 0; j < num2; ++j) {
      vertex v;
      int tmp = TriangleIntersection(surface1[i], surface2[j], v);
      if (tmp > 0) {
        vecIntersectTri.push_back(surface1[i]);
        vecIntersectTri.push_back(surface2[j]);
        vecIntersPoints.push_back(v);
      }
    }
  }

  // unique
  vecIntersectTri.erase(
      std::unique(vecIntersectTri.begin(), vecIntersectTri.end()),
      vecIntersectTri.end());
  vecIntersPoints.erase(
      std::unique(vecIntersPoints.begin(), vecIntersPoints.end()),
      vecIntersPoints.end());
  return;
}

static void SurfaceIntersection(std::vector<vertex>& surface1,
                                std::vector<vertex>& surface2,
                                std::vector<vertex>& vecIntersectTri,
                                std::vector<vertex>& vecIntersectPoints) {
  std::vector<Triangle> surface_tri1;
  std::vector<Triangle> surface_tri2;
  int num1 = surface1.size();

  num1 = num1 / 3;
  for (int i = 0; i < num1; ++i) {
    Triangle tri;
    tri.v1 = surface1[3 * i];
    tri.v2 = surface1[3 * i + 1];
    tri.v3 = surface1[3 * i + 2];
    surface_tri1.push_back(tri);
  }

  int num2 = surface2.size();

  num2 = num2 / 3;
  for (int i = 0; i < num2; ++i) {
    Triangle tri;
    tri.v1 = surface2[3 * i];
    tri.v2 = surface2[3 * i + 1];
    tri.v3 = surface2[3 * i + 2];
    surface_tri2.push_back(tri);
  }

  std::vector<Triangle> intersect_tri;
  SurfaceIntersection(surface_tri1, surface_tri2, intersect_tri,
                      vecIntersectPoints);

  vecIntersectTri.clear();
  int num = intersect_tri.size();
  for (int i = 0; i < num; ++i) {
    vecIntersectTri.push_back(intersect_tri[i].v1);
    vecIntersectTri.push_back(intersect_tri[i].v2);
    vecIntersectTri.push_back(intersect_tri[i].v3);
  }
  return;
}

// zbuffer
vertex project(glm::mat4 transform, vertex& v) {
  vertex res;
  res.x = transform[0][0] * v.x + transform[0][1] * v.y + transform[0][2] * v.z;
  res.y = transform[1][0] * v.x + transform[1][1] * v.y + transform[1][2] * v.z;
  res.z = v.z;
  return res;
}

float Area(vertex& v, vertex& v1, vertex& v2) {
  float farea = v1.x * v2.y + v2.x * v.y + v.x * v1.y - v1.x * v.y -
                v2.x * v1.y - v.x * v2.y;
  return farea * 0.5;
}

float Angle(vertex& v, vertex& v1, vertex& v2) {
  float v_v1[3] = {v.x - v1.x, v.y - v1.y};
  float v_v2[3] = {v.x - v2.x, v.y - v2.y};
  float len1 = sqrt(v_v1[0] * v_v1[0] + v_v1[1] * v_v1[1]);
  float len2 = sqrt(v_v2[0] * v_v2[0] + v_v2[1] * v_v2[1]);
  float fdot = v_v1[0] * v_v2[0] + v_v1[1] * v_v2[1];
  if (len1 * len2 == 0) return 0;

  return acos(fdot / (len1 * len2));
}

// angle sum = 360
bool IsVertexInTriangle(vertex& v, vertex& v1, vertex& v2, vertex& v3) {
  /*float fangle1=Angle(v,v1,v2);
  float fangle2=Angle(v,v3,v2);
  float fangle3=Angle(v,v1,v3);
  float pi2=2*3.1415927;
  if(fabs(fangle1+fangle2+fangle3-pi2)<0.1)
      return 1;
  return 0;*/
  float farea1 = Area(v1, v2, v);
  float farea2 = Area(v1, v3, v);
  float farea3 = Area(v2, v3, v);

  float farea = Area(v1, v2, v3);
  if (fabs(farea1 + farea2 + farea3 - farea) < 0.001) return true;

  return false;
}

void barycentric(vertex& v, vertex& v1, vertex& v2, vertex& v3, float uvz[3]) {
  float s1[3];
  float s2[3];
  s2[0] = v3.y - v1.y;
  s2[1] = v2.y - v1.y;
  s2[2] = v1.y - v.y;
  s1[0] = v3.x - v1.x;
  s1[1] = v2.x - v1.x;
  s1[2] = v1.x - v.x;
  float u[3];
  Vcross(u, s1, s2);
  uvz[0] = -1;
  uvz[1] = 1;
  uvz[2] = 1;
  if (std::abs(u[2]) > 1e-2)  // dont forget that u[2] is integer. If it is zero
                              // then triangle ABC is degenerate
  {
    uvz[0] = 1.0 - (u[0] + u[1]) / u[2];
    uvz[1] = u[1] / u[2];
    uvz[2] = u[0] / u[2];
  }
  return;
}

void zbuffer(glm::mat4 transform, std::vector<vertex>& allvertex,
             std::vector<vertex>& intersecVertex,
             std::vector<vertex>& hideVertex) {
  // glm::mat transform;
  int width = 512;
  int height = 512;
  double* pZbuffer = (double*)malloc(sizeof(double) * width * height);
  for (int i = 0; i < width; ++i)
    for (int j = 0; j < height; ++j) pZbuffer[j * width + i] = 1e10;

  float spcx = 2.0 / width;  //[-1,1]
  float spcy = 2.0 / height;

  int num1 = allvertex.size() / 3;
  for (int i = 0; i < num1; ++i) {
    // triangle projection
    vertex p_v1 = project(transform, allvertex[3 * i]);
    vertex p_v2 = project(transform, allvertex[3 * i + 1]);
    vertex p_v3 = project(transform, allvertex[3 * i + 2]);

    // x/y range in screen,left corner (-1,-1)
    float x0 = p_v1.x < p_v2.x ? p_v1.x : p_v2.x;
    x0 = x0 < p_v2.x ? x0 : p_v2.x;
    float x1 = p_v1.x > p_v2.x ? p_v1.x : p_v2.x;
    x1 = x1 > p_v2.x ? x1 : p_v2.x;
    float y0 = p_v1.y < p_v2.y ? p_v1.y : p_v2.y;
    y0 = y0 < p_v2.y ? y0 : p_v2.y;
    float y1 = p_v1.y > p_v2.y ? p_v1.y : p_v2.y;
    y1 = y1 > p_v2.y ? y1 : p_v2.y;

    int c0 = (x0 + 1) / spcx, c1 = (x1 + 1) / spcx;
    int r0 = (y0 + 1) / spcy, r1 = (y1 + 1) / spcy;

    c0 = c0 < 0 ? 0 : c0;
    c1 = c1 > width - 1 ? width - 1 : c1;
    r0 = r0 < 0 ? 0 : r0;
    r1 = r1 > height - 1 ? height - 1 : r1;
    for (int r = r0; r <= r1; ++r) {
      for (int c = c0; c <= c1; ++c) {
        float xx = -1 + c0 * spcx;
        float yy = -1 + r0 * spcy;
        vertex v;
        v.x = xx;
        v.y = yy;
        bool bInTri = IsVertexInTriangle(v, p_v1, p_v2, p_v3);
        if (bInTri) {
          float uvz[3];
          barycentric(v, p_v1, p_v2, p_v3, uvz);
          float zz = p_v1.z * uvz[0] + p_v2.z * uvz[1] + p_v3.z * uvz[2];

          if (pZbuffer[r * width + c] > zz) pZbuffer[r * width + c] = zz;
        }
      }
    }
  }

  //
  int numIntersect = intersecVertex.size();
  for (int i = 0; i < numIntersect; ++i) {
    vertex v = project(transform, intersecVertex[i]);
    int r = (v.y + 1) / spcy;
    int c = (v.x + 1) / spcx;
    if (r < 0 || r > height - 1 || c < 0 || c > width - 1) continue;

    double zz = pZbuffer[r * width + c];
    if (pZbuffer[r * width + c] < v.z) hideVertex.push_back(v);
  }
  free(pZbuffer);
  pZbuffer = NULL;
  return;
}

#endif  // INTERSECT_H
