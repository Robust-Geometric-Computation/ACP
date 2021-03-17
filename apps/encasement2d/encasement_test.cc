#include "acp/dcel/dcel.h"
#include "acp/encasement2d/critical_points.h"
#include "acp/encasement2d/encasement2d.h"
#include "acp/poly/poly.h"
#include "acp/poly/poly2d.h"

#include <stdio.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

vector<PTR<Object<PV2>>> roots;
vector<PTR<Object<Poly2>>> curves;

PTR<Object<PV2>> ll;
PTR<Object<PV2>> ur;

Encasement* encasement = nullptr;

int connected_components(PTR<Object<Poly2>> f, Encasement* enc) {
  vector<PTR<Face>> faces = enc->faces;

  for (int i = 0; i < faces.size(); i++) {
    faces[i]->clear_parent();
  }

  for (int i = 0; i < faces.size(); i++) {
    PTR<Edge> e = faces[i]->e;

    do {
      if (e->intersection && e->poly == f) {
        if (e->face && e->twin->face) e->face->unionize(e->twin->face);
      }

    } while ((e = e->next) != faces[i]->e);
  }

  set<PTR<Face>> s;

  for (int i = 0; i < faces.size(); i++) {
    if (faces[i]->contains(f)) s.insert(faces[i]->find());
  }

  return s.size();
}

void verify() {
  vector<int> cc;

  for (int i = 0; i < curves.size(); i++) {
    cc.push_back(connected_components(curves[i], encasement));
  }

  Encasement* enc_ver = new Encasement(ll, ur, 27);

  enc_ver->init();

  enc_ver->calculate(curves, true);

  // vector<PTR<Object<PV2>>> roots = enc_ver->getRoots();

  for (int i = 0; i < curves.size(); i++) {
    int vcc = connected_components(curves[i], enc_ver);

    for (int j = 0; j < enc_ver->sites[i].size(); j++) {
      PV2<Parameter> pp = enc_ver->sites[i][j]->get<Parameter>();
      printf("site[%d][%d]: (%.16f, %.16f)\n", i, j, pp.x.mid(), pp.y.mid());
    }

    printf(
        "curve %d: %d connected components %d verified connected components "
        "[%s]\n",
        i, cc[i], vcc, (cc[i] == vcc ? "PASS" : "FAIL"));
  }
}

void calculate() {
  // calculate the encasement of the ellipse
  if (curves.size() > 0) {
    encasement->init();

    // printf("calculating up to %d steps\n", encasement->STEPS);
    auto t1 = std::chrono::high_resolution_clock::now();
    encasement->calculate(curves);

    auto r1 = std::chrono::high_resolution_clock::now();

    roots = encasement->getRoots();

    PV2<Parameter> apprx[roots.size()];

    for (int i = 0; i < roots.size(); i++) {
      apprx[i] = roots[i]->getApprox();
      //   printf("r[%d]: [%.16g, %.16g] [%.16g, %.16g]\n", i, p.x.lb(),
      //   p.x.ub(), p.y.lb(), p.y.ub());
      //   //printf("dist: %.16e %.16e\n", p.x.ub() - p.x.lb(), p.y.ub() -
      //   p.y.lb());
    }

    auto r2 = std::chrono::high_resolution_clock::now();
    cout << "double aprx time: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(r2 - r1)
                .count()
         << " ms\n";

    printf("intersections size: %ld\n", roots.size());

    printf("faces size: %ld\n", encasement->faces.size());

    /*
        //printf("%ld intersection points\n", roots.size());
    */

    auto t2 = std::chrono::high_resolution_clock::now();
    cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
                .count()
         << " ms\n";

    // verify();
  }
}

int main(int argc, char** argv) {
  enable();

  ll = new InputPoint(-2, -2);
  ur = new InputPoint(2, 2);

  // r[0]: [0.7726254415927478, 0.7726254415932418] [-1.220957532372445,
  // -1.220957532369735]

  // ll = new InputPoint(0.65, -1.44);
  // ur = new InputPoint(0.89, -1.00);

  if (argc < 2) {
    printf("usage: ./encasement_test <file_name>\n");
    exit(1);
  }

  int power = -27;

  if (argc > 2) {
    if (argc == 3) {
      power = atoi(argv[2]);
    } else if (argc == 6) {
      ll = new InputPoint(atof(argv[2]), atof(argv[3]));
      ur = new InputPoint(atof(argv[4]), atof(argv[5]));
    } else {
      printf("usage: ./encasement_test <file_name> <llx> <lly> <urx> <ury>\n");
      printf(
          "       ./encasement_test <file_name> <power> (perturbs by "
          "2^power)\n");
      exit(1);
    }
  }

  encasement = new Encasement(ll, ur, 27);

  ifstream inFile(argv[1]);
  string line;

  double prev_delta = Parameter::delta;
  Parameter::delta = ::pow(2.0, power);

  std::vector<PTR<Object<Poly2D>>> poly2D_curves;

  while (getline(inFile, line)) {
    istringstream buf(line);
    istream_iterator<std::string> beg(buf), end;

    int count = 0;

    vector<Parameter> coef;
    vector<int> pow;

    int max_i = 0;
    int max_j = 0;

    for (; beg != end;) {
      double c = stod((*beg));
      beg++;
      int i = stoi((*beg));
      beg++;
      int j = stoi((*beg));
      beg++;

      max_i = max(i, max_i);
      max_j = max(j, max_j);

      coef.push_back(Parameter::input(c));
      pow.push_back(i);
      pow.push_back(j);
    }

    curves.push_back(
        new Ellipse(Poly2<Parameter>(coef.size(), &coef[0], &pow[0])));

    //PTR<Object<Poly2D>> poly = new Poly2D<Parameter>(max_i, max_j);
    //for(int n = 0; n < static_cast<int>(coef.size()); ++n) {
    //  int i = pow[2*n];
    //  int j = pow[2*n + 1];
    //  poly.set(i, j, coef[n]);
    //}

  }

  Parameter::delta = prev_delta;

  try {
    calculate();
  } catch (PrecisionException e) {
    printf("Precision Exception\n");
  } catch (SignException e) {
    printf("Sign Exception\n");
  }



  //vector<PTR<Object<PV2>>> poly2d_roots = getRoots(PTR<Object<Poly2D>> f, PTR<Object<Poly2D>> g,
  //                                 PV2<Parameter> box)

  disable();
}
