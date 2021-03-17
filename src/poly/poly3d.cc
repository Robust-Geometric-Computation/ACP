#include "acp/poly/poly3d.h"

vector<PTR<Object<PV3>>> boundaryRoots(PTR<Object<Poly3D>> f,
                                       PTR<Object<Poly3D>> g,
                                       const PV3<Parameter>& box, int xyz,
                                       int lu) {
  static int count;
  ++count;
  PTR<Object<Scalar>> s =
      new Object<Scalar>(lu ? box[xyz].ubP() : box[xyz].lbP());
  int ix = (xyz + 1) % 3;
  int iy = (xyz + 2) % 3;

#ifdef USE_ENCASEMENT
  PTR<Object<Poly2>> f2 = new Sub1Poly3D(f, xyz, s);
  PTR<Object<Poly2>> g2 = new Sub1Poly3D(g, xyz, s);
  PTR<Object<PV2>> l =
      new Object<PV2>(PV2<Parameter>(box[ix].lbP(), box[iy].lbP()));
  PTR<Object<PV2>> u =
      new Object<PV2>(PV2<Parameter>(box[ix].ubP(), box[iy].ubP()));
  Encasement* enc = new Encasement(l, u, rand());
  enc->print = false;
  enc->init();
  vector<PTR<Object<Poly2>>> curves = {f2, g2};
  enc->calculate(curves, false);
  vector<Root2PTR> roots2 = enc->getRoots();
  delete enc;

#ifdef CHECK_ENCASEMENT
  PTR<Object<Poly2D>> f2d = new Sub1dPoly3D(f, xyz, s);
  PTR<Object<Poly2D>> g2d = new Sub1dPoly3D(g, xyz, s);
  PV2<Parameter> box2d(box[ix], box[iy]);
  vector<PTR<Object<PV2>>> roots2d = getRoots(f2d, g2d, box2d);

  if (roots2.size() != roots2d.size()) {
    cout << "count " << count << " wrong number of roots" << endl;

    for (int i = 0; i < roots2d.size(); ++i) {
      PV2<Parameter> root2d = roots2d[i]->getApprox();
      int j = 0;
      for (; j < roots2.size(); ++j)
        if ((roots2[i]->getApprox().x - root2d.x).sign(false) == 0) break;
      if (j == roots2.size())
        cout << "missing roots2d[" << i << "] " << root2d.x.mid() << " "
             << root2d.y.mid() << endl;
    }
    assert(0);
  }
#endif

  vector<PTR<Object<PV3>>> roots;
  for (int i = 0; i < roots2.size(); ++i) {
    PTR<Object<PV2>> root2 = roots2[i];
    roots.push_back(new P2to3(root2, s, xyz));
  }
#else
  PTR<Object<Poly2D>> f2d = new Sub1dPoly3D(f, xyz, s);
  PTR<Object<Poly2D>> g2d = new Sub1dPoly3D(g, xyz, s);
  PV2<Parameter> box2d(box[ix], box[iy]);
  vector<PTR<Object<PV2>>> roots2d = getRoots(f2d, g2d, box2d);

  vector<PTR<Object<PV3>>> roots;
  for (int i = 0; i < roots2d.size(); ++i) {
    PTR<Object<PV2>> root2 = roots2d[i];
    roots.push_back(new P2to3(root2, s, xyz));
  }
#endif

#ifdef VERBOSE
  if (roots.size() != 0) {
    for (int i = 0; i < roots.size(); ++i) {
#ifndef USE_ENCASEMENT
      int sf2r = Poly2Dval(f2d, roots2d[i]);
      assert(sf2r == 0);
#endif
      int sfr = Poly3Dval(f, roots[i]);
      assert(sfr == 0);
      PV3<Parameter> r = roots[i]->getApprox();
      Parameter fr = f->getApprox().value(r);
      if (fr.sign(false) != 0) {
        cerr << fr.lb() << " " << fr.ub() << " " << fr.ub() - fr.lb() << endl;
        int sfr = Poly3Dval(f, roots[i]);
        cout << sfr << endl;
#ifndef USE_ENCASEMENT
        int sf2r = Poly2Dval(f2d, roots2d[i]);
        cout << sf2r << endl;
#endif
      }
      assert(f->getApprox().value(r).sign(false) == 0);
      assert(g->getApprox().value(r).sign(false) == 0);
#ifdef VERBOSE
      cout << "xyz " << xyz << " lu " << lu << " "
           << "root " << r.x.mid() << " " << r.y.mid() << " " << r.z.mid()
           << endl;
#endif
    }
  }
#endif

  return roots;
}

vector<PTR<Object<PV3>>> boundaryRoots(PTR<Object<Poly3D>> f,
                                       PTR<Object<Poly3D>> g,
                                       const PV3<Parameter>& box) {
  vector<PTR<Object<PV3>>> all;
  for (int xyz = 0; xyz < 3; ++xyz)
    for (int lu = 0; lu < 2; ++lu) {
      vector<PTR<Object<PV3>>> roots = boundaryRoots(f, g, box, xyz, lu);
      all.insert(all.end(), roots.begin(), roots.end());
    }
  return all;
}

vector<PTR<Object<PV3>>> getRoots(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g,
                                  PTR<Object<Poly3D>> h, PV3<Parameter> box) {
  vector<PTR<Object<PV3>>> roots;

  Poly3D<Parameter> fp = f->getApprox();
  Poly3D<Parameter> gp = g->getApprox();
  Poly3D<Parameter> hp = h->getApprox();

  vector<PV3<Parameter>> boxes;
  boxes.push_back(box);
  while (boxes.size() > 0) {
    PV3<Parameter> b = boxes.back();
    boxes.pop_back();

#ifdef VERBOSE
    // bool sub = p.x.subset(b.x) && p.y.subset(b.y) && p.z.subset(b.z);
    cout << "b " << b.x.lb() << " " << b.x.ub() << " " << b.y.lb() << " "
         << b.y.ub() << " " << b.z.lb() << " " << b.z.ub() << " ";
    // if (sub) cout << "*";
    cout << endl;
#endif

    int sf = Descartes3D(f, b);
    if (sf != 0) {
#ifdef VERBOSE
      cout << "Descartes3D(f, b) = " << sf << endl;
      // assert(!sub);
#endif
      continue;
    }

    int sg = Descartes3D(g, b);
    if (sg != 0) {
#ifdef VERBOSE
      cout << "Descartes3D(g, b) = " << sg << endl;
      // assert(!sub);
#endif
      continue;
    }

    int sh = Descartes3D(h, b);
    if (sh != 0) {
#ifdef VERBOSE
      cout << "Descartes3D(h, b) = " << sh << endl;
      // assert(!sub);
#endif
      continue;
    }

    PV3<Parameter> df(fp.der(0, b), fp.der(1, b), fp.der(2, b));
    PV3<Parameter> dg(gp.der(0, b), gp.der(1, b), gp.der(2, b));
    PV3<Parameter> dh(hp.der(0, b), hp.der(1, b), hp.der(2, b));

    Parameter triple = df.cross(dg).dot(dh);
    int sTriple = triple.sign(false);
#ifdef VERBOSE
    cout << "triple.sign(false) = " << sTriple << endl;
#endif

    if (sTriple != 0) {
      vector<PTR<Object<PV3>>> bRoots = boundaryRoots(f, g, b);
#ifdef VERBOSE
      cout << "bRoots.size() " << bRoots.size() << endl;
#endif
      if (bRoots.size() == 0) continue;
      assert(bRoots.size() % 2 == 0);
      int nPlus = 0, nMnus = 0;
      for (int i = 0; i < bRoots.size(); ++i) {
        int s = Poly3Dval(h, bRoots[i]);
        if (s > 0)
          nPlus++;
        else if (s < 0)
          nMnus++;
        else
          assert(0);
      }
#ifdef VERBOSE
      cout << "nPlus " << nPlus << " nMnus " << nMnus << endl;
#endif

      if (nPlus == 0 || nMnus == 0) continue;

      if (nPlus == 1 || nMnus == 1) {
        try {
          // Newton it down as far as possible, it must make one successful step
          // save that count for the root3d and also count successful newton
          // steps in calculate, if it's less than the current minimum that's
          // suspicious.

          // just use b!

          PV3<Parameter> b_in;
          PV3<Parameter> b_out = b;
          int newton_count = 0;

          // If the interval isn't trivial, try to shrink with Newton's
          bool trivial = (b.x.mid() == b.x.lb() || b.x.mid() == b.x.ub() ||
                          b.y.mid() == b.y.lb() || b.y.mid() == b.y.ub() ||
                          b.z.mid() == b.z.lb() || b.z.mid() == b.z.ub());

          if (!trivial) {
            // Iterate Newton's until the iteration doesn't shrink the input
            // interval.
            do {
              b_in = b_out;
              b_out = newton<Parameter>(fp, gp, hp, b_in);
              newton_count++;
            } while (b_out.x.subset(b_in.x) && b_out.y.subset(b_in.y) &&
                     b_out.z.subset(b_in.z));
          }

          // If we are actually a subset of the input interval and we've made
          // more than 1 newton step.
          if (trivial || newton_count > 1) {
            PTR<Object<PV3>> root = new Root3D(f, g, h, b_out);
#ifdef MAIN
            cout << "Found a root!" << endl;
            PV3<Parameter> rp = root->getApprox();
            cout << rp.x.mid() << " " << rp.y.mid() << " " << rp.z.mid()
                 << endl;
#endif
            roots.push_back(root);
            continue;
          }
#ifdef VERBOSE
          else
            cout << "Newton failed" << endl;
#endif
        } catch (SignException e) {
#ifdef VERBOSE
          cout << "Newton threw an exception" << endl;
#endif
        }
      }
    }

#ifdef VERBOSE
    cout << "splitting" << endl;
#endif
    int iMax = 0;
    for (int i = 1; i < 3; ++i)
      if (b[i].ub() - b[i].lb() > b[iMax].ub() - b[iMax].lb()) iMax = i;
    PV3<Parameter> b2 = b;
    b2[iMax] = Parameter::constant(b[iMax].lb())
                   .interval(Parameter::constant(b[iMax].mid()));
    boxes.push_back(b2);
    assert(b2[iMax].intervalWidth() < b[iMax].intervalWidth());

    b2[iMax] = Parameter::constant(b[iMax].mid())
                   .interval(Parameter::constant(b[iMax].ub()));
    boxes.push_back(b2);
    assert(b2[iMax].intervalWidth() < b[iMax].intervalWidth());
  }

  return roots;
}

bool containsZero(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g,
                  PTR<Object<Poly3D>> h, PV3<Parameter>& b) {
  vector<PTR<Object<PV3>>> bRoots = boundaryRoots(f, g, b);
  assert(bRoots.size() % 2 == 0);
  if (bRoots.size() == 0) return false;

  int nPlus = 0, nMnus = 0;
  for (int i = 0; i < bRoots.size(); ++i) {
    int s = Poly3Dval(h, bRoots[i]);
    if (s > 0)
      nPlus++;
    else if (s < 0)
      nMnus++;
    else
      assert(0);
  }

  // return nPlus > 0 && nMnus > 0;
  return nPlus % 2;
}

PV3<Parameter> subdivide(PTR<Object<Poly3D>> f, PTR<Object<Poly3D>> g,
                         PTR<Object<Poly3D>> h, PV3<Parameter>& box) {
  int iMax = 0;
  for (int i = 1; i < 3; ++i)
    if (box[i].ub() - box[i].lb() > box[iMax].ub() - box[iMax].lb()) iMax = i;

  PV3<Parameter> b = box;
  b[iMax] = Parameter::constant(box[iMax].lb())
                .interval(Parameter::constant(box[iMax].mid()));
  assert(b[iMax].intervalWidth() < box[iMax].intervalWidth());

  if (containsZero(f, g, h, b)) return b;

  b[iMax] = Parameter::constant(box[iMax].mid())
                .interval(Parameter::constant(box[iMax].ub()));
  assert(containsZero(f, g, h, b));
  assert(b[iMax].intervalWidth() < box[iMax].intervalWidth());
  return b;
}
