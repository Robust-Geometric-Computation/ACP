#include "acp/resultant/resultant.h"
#include "acp/poly/root2.h"

NoIntersectionException noIntersectionException;

bool Resultant::useEuclid = true;

template <class N>
Poly<N> SylvesterMatrix<N>::get(int i, int j) {
  return m[i][j];
}

template <class N>
void SylvesterMatrix<N>::createMatrix() {
  m = vector<vector<Poly<N>>>(n);

  for (int i = 0; i < n; i++) {
    vector<Poly<N>> v(n);

    for (int j = 0; j < n; j++) {
      int ii = i;
      int jj = j;

      // coef comes from f
      if (i < gd) {
        jj = (n + j - i) % n;
        v[j] = fv[jj];
      }
      // coef comes from g
      else {
        ii -= gd;
        jj = (n + j - ii) % n;
        v[j] = gv[jj];
      }
    }

    m[i] = v;
  }
}

template <class N>
int SylvesterMatrix<N>::findMinDeg(int j) {
  int d;
  int x = -1;
  for (int i = 0; i <= j; i++) {
    if (m[i][j].deg() != -1) {
      if (x == -1 || m[i][j].deg() < d) {
        d = m[i][j].deg();
        x = i;
      }
    }
  }
  return x;
}

template <class N>
void SylvesterMatrix<N>::swapRow(int i, int j) {
  vector<Poly<N>> temp = m[i];
  m[i] = m[j];
  m[j] = temp;
}

template <class N>
void SylvesterMatrix<N>::reduce(int j, int i1, int i2) {
  Poly<N> f = m[i1][j];
  Poly<N> g = m[i2][j];

  int d = g.deg() - f.deg();
  N a = g.a[g.deg()] / f.a[f.deg()];

  // Poly<N> h(a, d); //a * (x) ^ d
  Poly<N> h;  // a * (x) ^ d
  h.add(a, d);
  for (int c = j; c >= 0; c--) {
    m[i2][c] = m[i2][c] - (h * m[i1][c]);
  }

  // We know we've eliminated the leading term of this poly
  // m[i2][j].d--;
  m[i2][j].a.pop_back();

  // While leading coefficient is exactly 0 (N::constant(0)) reduce the degree
  // TODO
  // IS this right???
  while (m[i2][j].deg() >= 0 && m[i2][j].a[m[i2][j].deg()].sign() == 0) {
    // m[i2][j].deg()--;
    m[i2][j].a.pop_back();
  }
}

/*
 * Fill the vector with n - deg(f) trailing 0 polynomials
 * Make sure n is set before calling
 */

template <class N>
vector<Poly<N>> SylvesterMatrix<N>::fillVector(Poly2<N> f) {
  vector<Poly<N>> v(n);
#ifdef BLEEN
  for (int i = 0; i < n; i++) {
    v[i] = Poly<N>();
  }
#endif

  // TODO optimize?

  // Fill the coef in decreasing degree
  for (int i = 0; i < f.size(); i++) {
    int j = f.m[2 * i];
    // v[f.d - j] = v[f.d - j] + Poly<N>(f.a[i], f.m[2*i+1]);
    v[f.d - j].add(f.a[i], f.m[2 * i + 1]);
  }

  return v;
}

template <class N>
void SylvesterMatrix<N>::construct() {
  fd = f.degX();
  gd = g.degX();
  n = fd + gd;

  fv = fillVector(f);
  gv = fillVector(g);
  /*
    for(int i = 0; i < fd; i++) {

      for(int j = i+1; j < fd; j++) {
        printf("Poly<N>nomialGCD[");
        fv[i].print();
        printf(", ");
        fv[j].print();
        printf("];\n\n");
      }

      for(int j = 0; j < gd; j++) {
        printf("Poly<N>nomialGCD[");
        fv[i].print();
        printf(", ");
        gv[j].print();
        printf("];\n\n");
      }
    }

    for(int i = 0; i < gd; i++) {
      for(int j = i+1; j < gd; j++) {
        printf("Poly<N>nomialGCD[");
        gv[i].print();
        printf(", ");
        gv[j].print();
        printf("];\n\n");
      }
    }
  */

  createMatrix();
}

template <class N>
int SylvesterMatrix<N>::countNegDeg(int j) {
  int c = 0;
  for (int i = 0; i <= j; i++) {
    if (m[i][j].deg() == -1) c++;
  }
  return c;
}

// M =
//
// [ a, b, c, 0]
// [ 0, a, b, c]
// [ d, e, f, 0]
// [ 0, d, e, f]
//
// >> det(M)
//
// ans =
//
// a^2*f^2 - a*b*e*f - 2*a*c*d*f + a*c*e^2 + b^2*d*f - b*c*d*e + c^2*d^2

// A^2*F^2 - A*B*E*F - 2*A*C*D*F + A*C*E^2 + B^2*D*F - B*C*D*E + C^2*D^2

template <class N>
Poly<N> SylvesterMatrix<N>::det4x4() {
  assert(get(0, 3).deg() < 0 & get(1, 0).deg() < 0 && get(2, 3).deg() < 0 &&
         get(3, 0).deg() < 0);
  Poly<N> a = get(0, 0);
  Poly<N> b = get(0, 1);
  Poly<N> c = get(0, 2);
  Poly<N> d = get(2, 0);
  Poly<N> e = get(2, 1);
  Poly<N> f = get(2, 2);

  // printf("\na\n");
  // a.print();
  // printf("\nb\n");
  // b.print();
  // printf("\nc\n");
  // c.print();
  // printf("\nd\n");
  // d.print();
  // printf("\ne\n");
  // e.print();
  // printf("\nf\n");
  // f.print();
  // printf("\n\n");

  // a^2*f^2 - a*b*e*f - 2*a*c*d*f + a*c*e^2 + b^2*d*f - b*c*d*e + c^2*d^2
  return (((a * a * f * f) - (a * b * e * f))
              .plusCtimes(N(-2), (a * c * d * f))) +
         (a * c * e * e) + (b * b * d * f) - (b * c * d * e) + (c * c * d * d);
}

template <class N>
Poly<N> SylvesterMatrix<N>::det() {
  //  do_error = 1;
  //  max_error = 0;

  if (f.degree() == 2 && g.degree() == 2) {
    return det4x4();
  }

  // QR decomposition algorithm

  // do elimination on each column until a single constant poly is the only term
  // left
  for (int j = n - 1; j > 0; j--) {
    // is there only non-zero term remaining
    int c = countNegDeg(j);

    while (c < j) {
      // find the row with the poly of min degree
      int i = findMinDeg(j);

      // go through all rows that are not the reducing row and do not have
      // degree -1 and reduce
      for (int x = 0; x <= j; x++) {
        if (x == i || m[x][j].deg() == -1) continue;
        reduce(j, i, x);
      }

      c = countNegDeg(j);
    }

    int r = 0;

    // find the remaining row
    for (int i = 0; i <= j; i++) {
      if (m[i][j].deg() != -1) r = i;
    }

    // the one left over needs to be degree 0
    assert(m[r][j].deg() == 0);

    // swap them so the term is on the jth diagonal
    if (r != j) swapRow(r, j);
  }

  // take the determinant of the lower triangular matrix
  Poly<N> res(N(1));
  for (int i = 0; i < n; i++) {
    res = res * m[i][i];
  }

  //  do_error = -1;
  //  max_error = -1;

  // return res;
  return m[0][0];
}

template <class N>
void SylvesterMatrix<N>::print() {
  printf("\\begin{matrix}\n");

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("(");
      get(i, j).print();
      printf(") %s ", j == n - 1 ? "\\\\" : " & ");
    }
    printf("\n");
  }
  printf("\\end{matrix}\n");
  printf("\n\n");
  fflush(stdout);

  // printf("fv: ");
  // for(int i = 0; i < fv.size(); i++) {
  //  fv[i].print();
  //  printf("\t");
  //}
  // printf("\n");

  // printf("gv: ");
  // for(int i = 0; i < gv.size(); i++) {
  //  gv[i].print();
  //  printf("\t");
  //}
  // printf("\n\n");
}

template <class N>
N ResultantEuclid<N>::getRi(int i) {
  return N::constant(i);
}

// Total degree of f[0] + f[1] y + f[2] y^2 + ....
template <class N>
int ResultantEuclid<N>::degree(vector<Poly<N>>& f) {
  int max = -1;
  for (int i = 0; i < f.size(); i++)
    if (max < f[i].deg() + i) max = f[i].deg() + i;
  return max;
}

template <class N>
vector<Poly<N>> ResultantEuclid<N>::fillVector(Poly2<N> f) {
  vector<Poly<N>> v(f.degX() + 1);

#ifdef BLEEN
  for (int i = 0; i < f.degX() + 1; i++) {
    v[i] = Poly<N>(0);
  }

  // Fill the coef in decreasing degree
  for (int i = 0; i < f.size(); i++) {
    int j = f.m[2 * i];
    v[j] = v[j] + Poly<N>(f.a[i], f.m[2 * i + 1]);
  }
#endif

  for (int i = 0; i < f.size(); i++) {
    int j = f.m[2 * i];
    v[j].add(f.a[i], f.m[2 * i + 1]);
  }

  return v;
}

// Resultant value after plugging in for x.
// Essentially gcd of sum f[i] y^i and sum g[i] y^i.
template <class N>
N ResultantEuclid<N>::resultant(vector<N>& f, vector<N>& g) {
  if (f.size() == 1) return f[0];

  N prod = N::constant(1);

  // Reduce g to g mod f as polynomials
  for (int i = g.size() - f.size(); i >= 0; i--) {
    prod = prod * f.back();
    if (f.back().sign() == 0) cout << "Oops: special case 2." << endl;

    N s = g.back() / f.back();
    g.pop_back();  // Leading term of g becomes zero.
    for (int j = 0; j < f.size() - 1; j++) g[i + j] = g[i + j] - s * f[j];
  }

  return prod * resultant(g, f);
}

// f[i] is coefficient of y^i, etc.
template <class N>
Poly<N> ResultantEuclid<N>::resultant(vector<Poly<N>>& f, vector<Poly<N>>& g) {
  if (f.size() > g.size()) return resultant(g, f);

  vector<N> fv(f.size());
  vector<N> gv(g.size());
  int bezout = degree(f) * degree(g);
  vector<N> values(bezout + 1);

  for (int i = 0; i <= bezout; i++) {
    N ri = NewtonPoly::newtonRoot<N>(i);

    fv.clear();
    for (int j = 0; j < f.size(); j++) fv.push_back(f[j].value(ri));
    if (fv.back().sign() == 0) cout << "Oops: special case." << endl;

    gv.clear();
    for (int j = 0; j < g.size(); j++) gv.push_back(g[j].value(ri));
    if (gv.back().sign() == 0) cout << "Oops: special case" << endl;

    values[i] = resultant(fv, gv);
  }

  return interpolate(values);
}

template <class N>
Poly<N> ResultantEuclid<N>::det() {
  vector<Poly<N>> fcoef = fillVector(f);
  vector<Poly<N>> gcoef = fillVector(g);

  return resultant(fcoef, gcoef);
}

template <class N>
Poly<N> normal(Poly<N> poly) {
  Poly<N> p(poly);

  int maxi = 0;
  for (unsigned i = 1u; i < p.a.size(); i++) {
    maxi = (p.a[i].abs() < p.a[maxi].abs() ? maxi : i);
  }

  double mid = p.a[maxi].mid();

  for (unsigned int i = 0u; i < p.a.size(); ++i) p.a[i] = p.a[i] / mid;

  return p;
}

template <class N>
Poly<N> Resultant::calculate() {
  // printf("BEGIN resultant calc precision [%d]\n", N::highPrecision);

  if (useEuclid) {
    //    do_error = 1;
    //    max_error = 0;

    ResultantEuclid<N> s(f->get<N>(), g->get<N>());

    Poly<N> pr = normal(s.det());

    //    do_error = -1;
    //    max_error = -1;

    // printf("END resultant calc precision [%d]\n", N::highPrecision);

    return pr;

  } else {
    SylvesterMatrix<N> s(f->get<N>(), g->get<N>());

    s.print();

    return normal(s.det());
  }

  // Poly<N> p_constant = s.det();
  // Poly<N> p_qr = s.det();

  // p_qr.print();
  // printf("\n\n\n");
  // p_constant.print();
  // printf("\n");
  // exit(1);
}

template <class N>
Poly<N> SubstitutePoly1::calculate() {
  Poly2<N> fp = f->get<N>();

  N xp = x->get<N>();
  N xacc = xp;

  // degree of the variable remaining after substitution
  int d = (i == 0 ? fp.degY() : fp.degX());

  // degree of variable being substituted
  int d_other = (i == 0 ? fp.degX() : fp.degY());

  // initial coef for ppoly1 result
  vector<N> a(d + 1);
  for (int j = 0; j < d + 1; j++) a[j] = N::constant(0);

  // memoize all possible powers of the substitution variable
  vector<N> xa(d_other + 1);
  xa[0] = N::constant(1);
  for (int j = 1; j < d_other + 1; j++) {
    xa[j] = xacc;
    xacc = xacc * xp;
  }

  for (int j = 0; j < fp.size(); j++) {
    int ind = fp.m[2 * j + (1 - i)];
    int sub_ind = fp.m[2 * j + i];

    if (sub_ind == 0)
      a[ind] = a[ind] + fp.a[j];
    else
      a[ind] = a[ind] + fp.a[j] * xa[sub_ind];
  }

  return Poly<N>(a);
}

template <class N>
PV2<N> ResultantRoot::calculate() {
  /*
    if(!verify()) {
      throw noIntersectionException;
    }
  */

  N xp = x->get<N>();
  N yp = y->get<N>();

  PV2<N> p(xp, yp);

  Root2PTR r = new Root2(f, g, PV2<Parameter>(p));

  return r->get<N>();
}
/*
int LessThan::sign() {
  return (a->get<N>() - b->get<N>()).sign();
}
*/

template Poly<Parameter> SubstitutePoly1::calculate<acp::Parameter>();
template Poly<MParameter> SubstitutePoly1::calculate<acp::MParameter>();
template Poly<PParameter> SubstitutePoly1::calculate<acp::PParameter>();
template Poly<Parameter> Resultant::calculate<acp::Parameter>();
template Poly<MParameter> Resultant::calculate<acp::MParameter>();
template Poly<PParameter> Resultant::calculate<acp::PParameter>();
