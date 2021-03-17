#include <algorithm>
#include <set>
#include <vector>

#include "acp/encasement3d/encasement_utils.h"

// Utilities for the encasement algorithm
namespace acp {
namespace encasement_utils {

// Determinant of 2x2 matrix of Doubles
//
// @param[in] m
//     2x2 matrix, row major
//
// @returns The determinant of m.
inline Double Determinant(const Double m[2][2]) {
  return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

// Determinant of 3x3 matrix of Doubles
//
// @param[in] m
//     3x3 matrix, row major
//
// @returns The determinant of m.
inline Double Determinant(const Double m[3][3]) {
  const Double c1[2][2] = {{m[1][1], m[1][2]}, {m[2][1], m[2][2]}};
  const Double c2[2][2] = {{m[1][2], m[1][0]}, {m[2][2], m[2][0]}};
  const Double c3[2][2] = {{m[1][0], m[1][1]}, {m[2][0], m[2][1]}};
  return m[0][0] * Determinant(c1) + m[0][1] * Determinant(c2) +
         m[0][2] * Determinant(c3);
}

// Solve 3x3 system of equations in Double using Cramer's rule
//
// @pre m != 0
//
// @param[in] m
//     The 3x3 matrix in the system.
// @param[in] b
//     The RHS of the system.
// @param[out] x
//     The solution to the system Ax = b. In MATLAB speak: x = A \ b.
void Cramer(const Double m[3][3], const Double b[3], Double* x) {
  Double s = Determinant(m);

  const Double c1[3][3] = {{b[0], m[0][1], m[0][2]},
                           {b[1], m[1][1], m[1][2]},
                           {b[2], m[2][1], m[2][2]}};

  const Double c2[3][3] = {{m[0][0], b[0], m[0][2]},
                           {m[1][0], b[1], m[1][2]},
                           {m[2][0], b[2], m[2][2]}};

  const Double c3[3][3] = {{m[0][0], m[0][1], b[0]},
                           {m[1][0], m[1][1], b[1]},
                           {m[2][0], m[2][1], b[2]}};

  Double s1 = Determinant(c1);

  Double s2 = Determinant(c2);

  Double s3 = Determinant(c3);
  x[0] = s1 / s;
  x[1] = s2 / s;
  x[2] = s3 / s;
}

// Finds the first instace of a PV3<Double> in points whose
// distance from p is less than tolerance
//
// @param[in] p
//     The point we're searching for.
// @param[in] points
//     The set of points in which we're searching.
// @param[in] tolerance (optional)
//     The distance threshold for considering a point close enough to p.
//
// @returns The index of the first PV3<Double> in points that matches the
// criteria
//         OR -1 if no such point exists.
int FindClosest(PV3<Double> p, std::vector<PV3<Double>> points,
                Double tolerance = 1e-6) {
  auto q = std::find_if(
      points.begin(), points.end(),
      [&](const PV3<Double>& a) { return (a - p).length2() <= tolerance; });

  if (q == points.end()) return -1;

  return static_cast<int>(q - points.begin());
}

// Trace the FG curve starting at the given starting vertex until it leave the
// box
//
// @param[in] f
//     The poly3 f(x,y,z).
// @param[in] g
//     The poly3 g(x,y,z).
// @param[in] box_l
//     The lower boundary corner of the bounding box.
// @param[in] box_u
//     The upper boundary corner of the bounding box.
// @param[in] starting_vert
//     The point at which to start the tracing.
//
// @returns A vector of PV3<Double> representing the polyline approximation of
// the fg
//          curve in the box starting at starting_vert.
std::vector<PV3<Double>> CurveTrace(const Poly3<Double>& f,
                                    const Poly3<Double>& g,
                                    const PV3<Double>& box_l,
                                    const PV3<Double>& box_u,
                                    const PV3<Double>& starting_vert) {
  // The polyline curve segment.
  std::vector<PV3<Double>> segment;

  // Establish a step size (start with just 0.01 of the smallest dimension)
  const Double min_dimension = std::min(
      box_u.x - box_l.x, std::min(box_u.y - box_l.y, box_u.z - box_l.z));

  const Double step_size = min_dimension * 1e-2;

  auto Inside = [](const PV3<Double>& p, const PV3<Double>& box_l,
                   const PV3<Double>& box_u) {
    return box_l.x <= p.x && p.x <= box_u.x && box_l.y <= p.y &&
           p.y <= box_u.y && box_l.z <= p.z && p.z <= box_u.z;
  };

  Poly3<Double> fx = f.der(0);
  Poly3<Double> fy = f.der(1);
  Poly3<Double> fz = f.der(2);
  Poly3<Double> gx = g.der(0);
  Poly3<Double> gy = g.der(1);
  Poly3<Double> gz = g.der(2);

  auto GradF = [&](const PV3<Double>& p) {
    return PV3<Double>(fx.value(p), fy.value(p), fz.value(p));
  };

  auto GradG = [&](const PV3<Double>& p) {
    return PV3<Double>(gx.value(p), gy.value(p), gz.value(p));
  };

  auto GradFCrossGradG = [&](const PV3<Double>& p) {
    return PV3<Double>(fy.value(p) * gz.value(p) - fz.value(p) * gy.value(p),
                       fz.value(p) * gx.value(p) - fx.value(p) * gz.value(p),
                       fx.value(p) * gy.value(p) - fy.value(p) * gx.value(p));
  };

  // Start from starting_vert and trace along Tan(f, g) == ∇F × ∇G
  PV3<Double> current_point = starting_vert;

  do {
    // Take a step
    PV3<Double> v = GradFCrossGradG(current_point).unit();
    PV3<Double> next_point = current_point + v * step_size;

    // After taking a step, take a Newton step back towards the FG curve
    // Step in the direction of Δp where:
    //   Δp * (∇F × ∇G)(p) = 0
    //   Δp * ∇F(p) = -F(p)
    //   Δp * ∇G(p) = -G(p)

    int newton_steps = 2;

    for (int i = 0; i < newton_steps; i++) {
      PV3<Double> c1 = GradFCrossGradG(next_point);
      PV3<Double> c2 = GradF(next_point);
      PV3<Double> c3 = GradG(next_point);
      Double fp = f.value(next_point);
      Double gp = g.value(next_point);

      Double m[3][3] = {
          {c1.x, c1.y, c1.z}, {c2.x, c2.y, c2.z}, {c3.x, c3.y, c3.z}};

      Double b[3] = {0.0, -fp, -gp};

      Double x[3];

      Cramer(m, b, &x[0]);

      next_point = next_point + PV3<Double>(x[0], x[1], x[2]);
    }

    // Detect if we've left the box OR if we've made no progress (should never
    // happen, so we'll ignore for now)
    if (Inside(next_point, box_l, box_u)) {
      segment.push_back(current_point);
      segment.push_back(next_point);
      current_point = next_point;
    } else {
      break;
    }

  } while (true);

  return segment;
}

// Computes poly chains of FG curves inside a box
//
// @param[in] f
//     The poly3 f(x,y,z).
// @param[in] g
//     The poly3 g(x,y,z).
// @param[in] box
//     The bounding region of the cell we're tracing within.
// @param[in] boundary_verts
//     List of boundary vertices of the cell we're tracing.
//     These serve as starting and ending points for the curve segments.
// @param[in] tolerance (optional)
//     Distance to consider two points on a curve segment to be equivalent.
//     Used to determine segment endpoints from boundary intersections
//
// @returns A map from index pair (i,j) to an open polyline.
//          i and j are the indices of the two input boundary vertices who
//          bound a single curve polyline.
//          Each polyline is represented by a vector of PV3.
//
//          If we trace a curve from vertex i that fails to reach another
//          boundary vertex, the pair (i, -1) is the key for that polyline.
std::map<std::pair<int, int>, std::vector<PV3<Parameter>>> FGCurves(
    const Poly3<Parameter>& f, const Poly3<Parameter>& g,
    const PV3<Parameter>& box,
    const std::vector<PV3<Parameter>>& boundary_verts, Double tolerance) {
  // Shadow the inputs with Double number type versions
  Poly3<Double> f_d(f);
  Poly3<Double> g_d(g);

  PV3<Double> box_l(box.x.lb(), box.y.lb(), box.z.lb());
  PV3<Double> box_u(box.x.ub(), box.y.ub(), box.z.ub());

  std::vector<PV3<Double>> boundary_verts_d;
  std::copy(boundary_verts.begin(), boundary_verts.end(),
            std::back_inserter(boundary_verts_d));

  // Create an empty map from boundary pairs to curve segment polylines
  std::map<std::pair<int, int>, std::vector<PV3<Parameter>>> curve_segments;

  // Create an empy set of vertices to mark processed verts
  std::set<int> processed;

  // Loop through each unprocessed boundary vert and curve trace.
  for (int i = 0; i < boundary_verts_d.size(); i++) {
    // Skip already processed vertices
    if (processed.find(i) != processed.end()) continue;

    // Trace the FG curve segment starting from boundary_verts[i].
    std::vector<PV3<Double>> segment =
        CurveTrace(f_d, g_d, box_l, box_u, boundary_verts_d[i]);

    if (segment.size() == 0) continue;

    // Find a boundary vertex close enough to the end of the chain.
    // If we can't find one we'll return a chain with unknown ending
    // vertex (vi, -1)
    int end_index = FindClosest(segment.back(), boundary_verts_d, tolerance);

    std::vector<PV3<Parameter>> segment_N;
    std::transform(segment.begin(), segment.end(),
                   std::back_inserter(segment_N),
                   [](const PV3<Double>& p) -> PV3<Parameter> {
                     return PV3<Parameter>::constant(p.x, p.y, p.z);
                   });

    curve_segments.insert(
        std::make_pair(std::make_pair(i, end_index), segment_N));

    // Keep track of processed boundary verts
    processed.insert(i);
    if (end_index != -1) processed.insert(end_index);
  }

  return curve_segments;
}

}  // namespace encasement_utils
}  // namespace acp
