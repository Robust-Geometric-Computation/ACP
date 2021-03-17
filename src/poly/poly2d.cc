#include "acp/poly/poly2d.h"

#undef VERBOSE
//#define VERBOSE

PTR<Object<PV2>> degenerateIntersection;
int yzAvoidSplit = 0;
double sAvoidSplit = 0;

bool isConsecutive(double a, double b) {
  double c = (a + b) / 2;
  return a == c || b == c;
}

bool isConsecutive(double a, double b, int yz) {
  if (isConsecutive(a, b)) return true;
  if (yzAvoidSplit == 0 || yz != yzAvoidSplit - 1) return false;
  double c = (a + b) / 2;
  return c == sAvoidSplit && isConsecutive(a, c) && isConsecutive(c, b);
}

class V2I_D : public Primitive {
  PTR<Object<PV2>> p;
  int xy;
  double d;
  DeclareSign { return p->get<N>()[xy] - d; }

 public:
  V2I_D(PTR<Object<PV2>> p, int xy, double d) : p(p), xy(xy), d(d) {}
};

vector<PTR<Object<PV2>>> boundaryRoots(PTR<Object<Poly2D>> f,
                                       const PV2<Parameter>& box, int xy,
                                       int lu) {
  PTR<Object<Scalar>> s =
      new Object<Scalar>(lu ? box[xy].ubP() : box[xy].lbP());
  PTR<Object<Poly>> f1 = new Sub1Poly2D(f, xy, s);

  int ix = (xy + 1) % 2;
  PTR<Object<Scalar>> l = new Object<Scalar>(Scalar<Parameter>(box[ix].lbP()));
  PTR<Object<Scalar>> u = new Object<Scalar>(Scalar<Parameter>(box[ix].ubP()));
  vector<PTR<Object<Scalar>>> roots1 = getRoots(f1, l, u);

  vector<PTR<Object<PV2>>> roots;
  for (int i = 0; i < roots1.size(); ++i) {
    PTR<Object<Scalar>> root1 = roots1[i];
    // int sf1r = PolyVal(f1, root1);
    // assert(sf1r == 0);
    roots.push_back(new P1to2(root1, s, xy));
  }

#ifdef VERBOSE
  if (roots.size() != 0) {
    for (int i = 0; i < roots.size(); ++i) {
      PV2<Parameter> r = roots[i]->getApprox();
      assert(f->getApprox().value(r).sign(false) == 0);
      cout << "    ";
      cout << "xy " << xy << " lu " << lu << " "
           << "root " << r.x.mid() << " " << r.y.mid() << endl;
    }
  }
#endif

  return roots;
}

vector<PTR<Object<PV2>>> boundaryRoots(PTR<Object<Poly2D>> f,
                                       const PV2<Parameter>& box) {
  vector<PTR<Object<PV2>>> all;
  for (int xy = 0; xy < 2; ++xy)
    for (int lu = 0; lu < 2; ++lu) {
      vector<PTR<Object<PV2>>> roots = boundaryRoots(f, box, xy, lu);
      all.insert(all.end(), roots.begin(), roots.end());
    }
  return all;
}

vector<PTR<Object<PV2>>> getRoots(PTR<Object<Poly2D>> f, PTR<Object<Poly2D>> g,
                                  PV2<Parameter> box) {
  vector<PTR<Object<PV2>>> roots;

  Poly2D<Parameter> fp = f->getApprox();
  Poly2D<Parameter> gp = g->getApprox();

  vector<PV2<Parameter>> boxes;
  boxes.push_back(box);
  while (boxes.size() > 0) {
    static int count = 0;
    ++count;

    PV2<Parameter> b = boxes.back();
    boxes.pop_back();
#ifdef MAIN
    bool sub = p.x.subset(b.x) && p.y.subset(b.y);
#endif
#ifdef VERBOSE
    cout << "    ";
    cout << "b " << b.x.lb() << " " << b.x.ub() << " " << b.y.lb() << " "
         << b.y.ub() << " ";

#ifdef MAIN
    if (sub) cout << "*";
#endif
    cout << endl;
#endif

    int sf = Descartes2D(f, b);
    if (sf != 0) {
#ifdef VERBOSE
      cout << "    ";
      cout << "Descartes2D(f, b) = " << sf << endl;
#endif
#ifdef MAIN
      assert(!sub);
#endif
      continue;
    }

    int sg = Descartes2D(g, b);
    if (sg != 0) {
#ifdef VERBOSE
      cout << "    ";
      cout << "Descartes2D(g, b) = " << sg << endl;
#endif
#ifdef MAIN
      assert(!sub);
#endif
      continue;
    }

    PV2<Parameter> df(fp.der(0, b), fp.der(1, b));
    PV2<Parameter> dg(gp.der(0, b), gp.der(1, b));

    Parameter cross = df.cross(dg);
    int sCross = cross.sign(false);
#ifdef VERBOSE
    cout << "    ";
    cout << "cross.sign(false) = " << sCross << endl;
#endif

    if (sCross == 0 && isConsecutive(b[0].lb(), b[0].ub(), 0) &&
        isConsecutive(b[1].lb(), b[1].ub(), 1)) {
      cout << "HACK!!!" << endl;
      cout << "b " << b.x.lb() << " " << b.x.ub() << " " << b.y.lb() << " "
           << b.y.ub() << endl;
      continue;
    }

    if (sCross != 0) {
      vector<PTR<Object<PV2>>> bRoots = boundaryRoots(f, b);
#ifdef VERBOSE
      cout << "    ";
      cout << "bRoots.size() " << bRoots.size() << endl;
#endif
      if (bRoots.size() == 0) continue;
      assert(bRoots.size() % 2 == 0);
      int nPlus = 0, nMnus = 0;
      for (int i = 0; i < bRoots.size(); ++i) {
        int s = Poly2Dval(g, bRoots[i]);
        if (s > 0)
          nPlus++;
        else if (s < 0)
          nMnus++;
        else
          assert(0);
      }

#ifdef VERBOSE
      cout << "    ";
      cout << "nPlus " << nPlus << " nMnus " << nMnus << endl;
#endif

      if (nPlus == 0 || nMnus == 0) continue;

      if (nPlus == 1 || nMnus == 1) {
        try {
          PV2<Parameter> b_in;
          PV2<Parameter> b_out = b;
          int newton_count = 0;

          do {
            b_in = b_out;
            b_out = newton<Parameter>(fp, gp, b_in);
            newton_count++;
          } while (b_out.x.subset(b_in.x) && b_out.y.subset(b_in.y));

#ifdef VERBOSE
          cout << "    ";
          cout << "Found a root!" << endl;
#endif
          if (newton_count > 1 || (isConsecutive(b[0].lb(), b[0].ub(), 0) &&
                                   isConsecutive(b[1].lb(), b[1].ub(), 1))) {
            PTR<Object<PV2>> root = new Root2D(f, g, b_out);
            roots.push_back(root);
#ifdef VERBOSE
            PV2<Parameter> rp = root->getApprox();
            cout << "    ";
            cout << rp.x.mid() << " " << rp.y.mid() << endl;
            assert(f->getApprox().value(rp).sign(false) == 0);
            assert(g->getApprox().value(rp).sign(false) == 0);
#endif
            continue;
          }
#ifdef VERBOSE
          else {
            cout << "    ";
            cout << "Newton failed" << endl;
          }
#endif
        } catch (SignException e) {
#ifdef VERBOSE
          cout << "    ";
          cout << "Newton threw an exception" << endl;
#endif
        }
      }
    }

#ifdef VERBOSE
    cout << "    ";
    cout << "splitting" << endl;
#endif
    int iMax = -1;
    for (int i = 0; i < 2; ++i)
      if (!isConsecutive(b[i].lb(), b[i].ub(), i) &&
          (iMax == -1 || (b[i].ub() - b[i].lb() > b[iMax].ub() - b[iMax].lb())))
        iMax = i;
    assert(iMax != -1);

    // A box of width 1e-6 is heuristically assumed to isolate the degenerate
    // point. This is more expensive than the approach in the paper.
    if (degenerateIntersection != 0 && b[iMax].intervalWidth() < 1e-6 &&
        V2I_D(degenerateIntersection, 0, b[0].lb()) > 0 &&
        V2I_D(degenerateIntersection, 0, b[0].ub()) < 0 &&
        V2I_D(degenerateIntersection, 1, b[1].lb()) > 0 &&
        V2I_D(degenerateIntersection, 1, b[1].ub()) < 0)
      continue;

    PV2<Parameter> b2 = b;
    double mid = b[iMax].mid();
    if (iMax == yzAvoidSplit - 1 && mid == sAvoidSplit) {
      if (!isConsecutive(b[iMax].lb(), mid))
        mid = (b[iMax].lb() + mid) / 2;
      else if (!isConsecutive(mid, b[iMax].ub()))
        mid = (b[iMax].ub() + mid) / 2;
      else
        assert(0);
    }
    b2[iMax] =
        Parameter::constant(b[iMax].lb()).interval(Parameter::constant(mid));
    assert(b2[iMax].intervalWidth() < b[iMax].intervalWidth());
    boxes.push_back(b2);

    b2[iMax] =
        Parameter::constant(mid).interval(Parameter::constant(b[iMax].ub()));
    assert(b2[iMax].intervalWidth() < b[iMax].intervalWidth());
    boxes.push_back(b2);
  }

  return roots;
}

bool containsZero(PTR<Object<Poly2D>> f, PTR<Object<Poly2D>> g,
                  PV2<Parameter>& b) {
  vector<PTR<Object<PV2>>> bRoots = boundaryRoots(f, b);
  assert(bRoots.size() % 2 == 0);
  if (bRoots.size() == 0) return false;

  int nPlus = 0, nMnus = 0;
  for (int i = 0; i < bRoots.size(); ++i) {
    int s = Poly2Dval(g, bRoots[i]);
    if (s > 0)
      nPlus++;
    else if (s < 0)
      nMnus++;
    else
      assert(0);
  }

  return nPlus % 2;
}

PV2<Parameter> subdivide(PTR<Object<Poly2D>> f, PTR<Object<Poly2D>> g,
                         PV2<Parameter>& box) {
  int iMax = 0;
  for (int i = 1; i < 2; ++i)
    if (box[i].ub() - box[i].lb() > box[iMax].ub() - box[iMax].lb()) iMax = i;

  PV2<Parameter> b = box;
  b[iMax] = Parameter::constant(box[iMax].lb())
                .interval(Parameter::constant(box[iMax].mid()));
  assert(b[iMax].intervalWidth() < box[iMax].intervalWidth());

  if (containsZero(f, g, b)) return b;

  b[iMax] = Parameter::constant(box[iMax].mid())
                .interval(Parameter::constant(box[iMax].ub()));
  assert(containsZero(f, g, b));
  assert(b[iMax].intervalWidth() < box[iMax].intervalWidth());
  return b;
}
