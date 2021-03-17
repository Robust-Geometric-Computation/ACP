#include "acp/dcel/dcel.h"
#include "acp/encasement2d/critical_points.h"
#include "acp/encasement2d/encasement2d.h"
#include "acp/encasement2d/plot.h"
#include "acp/poly/poly.h"
#include "acp/poly/root.h"

#include <stdio.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef ENCASEMENT_GPU
#include "gpu/encasement_gpu.h"
#endif

#ifdef OSX
#include <Accelerate/Accelerate.h>
#else
extern "C" {
#include <clapack.h>
}
#endif

using namespace std;

vector<PTR<Object<PV2>>> roots;
int fsize = -1;

int subdiv = (1 << 3);

Encasement* encasement = nullptr;

int g_index = -1;
PTR<Edge> g_edge;

bool draw_intersections = false;
bool draw_crit = false;

bool dfirst = true;

bool draw_roots = false;

bool draw_triangulation = false;

vector<PTR<Triangle>> triangulation;

int zoom_index;

void calculate();
bool calculated = false;

bool draw_encasement = true;

bool draw_gpu_rects = false;

int STEPS = 0;

int root_index = -1;

vector<Poly2<Parameter>*> polys;
vector<PTR<Object<Poly2>>> curves;

PTR<Object<PV2>> click = nullptr;

typedef struct {
  double a, b, xc, yc, phi;

} Ellipse_struct;

vector<Ellipse_struct> ellipses;

void updateTexture();

void curve_trace2(PTR<Face> face);

Ellipse_struct ellipse_to_struct(Poly2<Parameter>* ellipse) {
  Ellipse_struct e;

  double a = ellipse->a[0].mid();
  double b = ellipse->a[1].mid() / 2;
  double c = ellipse->a[2].mid();
  double d = ellipse->a[3].mid() / 2;
  double f = ellipse->a[4].mid() / 2;
  double g = ellipse->a[5].mid();

  e.xc = ((c * d - b * f) / (b * b - a * c));
  e.yc = ((a * f - b * d) / (b * b - a * c));

  e.a =
      (sqrt(2 *
            (a * f * f + c * d * d + g * b * b - 2 * b * d * f - a * c * g)) /
       sqrt((b * b - a * c) * (sqrt((a - c) * (a - c) + 4 * b * b) - (a + c))));

  e.b = (sqrt(2 *
              (a * f * f + c * d * d + g * b * b - 2 * b * d * f - a * c * g)) /
         sqrt((b * b - a * c) *
              (-1 * sqrt((a - c) * (a - c) + 4 * b * b) - (a + c))));

  if (a < c) {
    e.phi = -(1 / (2 * atan2((a - c), (2 * b))));
  } else {
    e.phi = -(1 / (2 * atan2((a - c), (2 * b)))) - (0.5 * M_PI);
  }

  return e;
}

void readFile() {
  FILE* file = fopen("data.raw", "r+");

  double coef[6];

  Parameter coefficients[6];
  int powers[12];

  powers[0] = 2;
  powers[1] = 0;
  powers[2] = 1;
  powers[3] = 1;
  powers[4] = 0;
  powers[5] = 2;
  powers[6] = 1;
  powers[7] = 0;
  powers[8] = 0;
  powers[9] = 1;
  powers[10] = 0;
  powers[11] = 0;

  printf("<<<reading poly>>>\n");

  while ((fread(coef, sizeof(double), 6, file) == 6)) {
    for (int i = 0; i < 6; i++) {
      coefficients[i] = Parameter::constant(coef[i]);
      printf("%s%f", i == 0 ? "\t" : " ", coef[i]);
    }
    printf("\n");

    polys.push_back(new Poly2<Parameter>(6, coefficients, powers));

    if (fread(coef, sizeof(double), 5, file) == 5) {
      Ellipse_struct e = {coef[0], coef[1], coef[2], coef[3], coef[4]};
      ellipses.push_back(e);
    } else {
      fprintf(stderr, "Error reading ellipses\n");
    }
  }

  fclose(file);

  updateTexture();
}

void writeFile() {
  FILE* file = fopen("data.raw", "w+");

  double coef[6];

  printf("<<<writing poly>>>\n");
  for (int i = 0; i < polys.size(); i++) {
    Poly2<Parameter>* p = polys[i];
    for (int j = 0; j < 6; j++) {
      coef[j] = p->a[j].mid();
      printf("%s%f", j == 0 ? "\t" : " ", coef[j]);
    }
    printf("\n");

    if ((fwrite(coef, sizeof(double), 6, file) != 6)) {
      fprintf(stderr, "Error writing data file\n");
      exit(1);
    }
    coef[0] = ellipses[i].a;
    coef[1] = ellipses[i].b;
    coef[2] = ellipses[i].xc;
    coef[3] = ellipses[i].yc;
    coef[4] = ellipses[i].phi;

    if ((fwrite(coef, sizeof(double), 5, file) != 5)) {
      fprintf(stderr, "Error writing ellipse to data file\n");
      exit(1);
    }
  }

  fclose(file);
}

void updateEllipse();
void init_ellipse(double a, double b, double xc, double yc, double phi);

int cell_index = -1;

enum State { NONE, CENTER, A, B, PHI, NUMBER_STATES };

bool input;

State state;

extern double draw_scale;
extern double x_off;
extern double y_off;

extern int counter;

GLuint textureName;
bool useTexture = false;

int w, h;

void updateTexture() {
  return;

  printf("updating texture\n");

  GLubyte c[w][h][4];

  int val[w][h];

  int** valx = new int*[w];
  for (int i = 0; i < w; i++) valx[i] = new int[h];

  int** valy = new int*[w];
  for (int i = 0; i < w; i++) valy[i] = new int[h];

  Parameter xy[2];

  for (int i = 0; i < w; i++) {
    for (int j = 0; j < h; j++) {
      c[i][j][0] = (GLuint)0;
      c[i][j][1] = (GLuint)0;
      c[i][j][2] = (GLuint)0;
      c[i][j][3] = (GLuint)0;
    }
  }

  printf("polys size(): %ld\n", polys.size());
  printf("test_poly size(): %ld\n", Rectangle::test_poly.size());

  bool d = Rectangle::test_poly.size() > 0;

  for (int k = 0; k < (d ? Rectangle::test_poly.size() : polys.size()); k++) {
    Poly2<Parameter>* p = (d ? &Rectangle::test_poly[k] : polys[k]);

    Poly2<Parameter> px = p->derX();
    Poly2<Parameter> py = p->derY();

    for (int i = 0; i < w; i++) {
      for (int j = 0; j < h; j++) {
        double y1 = ((2.0 * ((i + 0.0) / w) - 1.0));
        double y2 = ((2.0 * ((i + 0.5) / w) - 1.0));
        double x1 = ((2.0 * ((j + 0.0) / h) - 1.0));
        double x2 = ((2.0 * ((j + 0.5) / h) - 1.0));

        int pos_count = 0;
        int neg_count = 0;

        xy[0] = Parameter::constant(x1);
        xy[1] = Parameter::constant(y1);
        int s1 = p->value(xy).sign(false);
        // xy[0] = Parameter::constant(x2); xy[1] = Parameter::constant(y1);
        // int s2 = p->value(xy).sign(false);
        // xy[0] = Parameter::constant(x1); xy[1] = Parameter::constant(y2);
        // int s3 = p->value(xy).sign(false);
        // xy[0] = Parameter::constant(x2); xy[1] = Parameter::constant(y2);
        // int s4 = p->value(xy).sign(false);

        int sx = px.value(xy).sign(false);
        int sy = py.value(xy).sign(false);

        if (abs(sx) != 1) {
          c[i][j][1] = 255;
          c[i][j][3] = 255;
        }
        if (abs(sy) != 1) {
          c[i][j][2] = 255;
          c[i][j][3] = 255;
        }

        valx[i][j] = sx;
        valy[i][j] = sy;

        // val[i][j] = s1+s2+s3+s4;
        val[i][j] = s1;

        if (abs(val[i][j]) != 1) {
          c[i][j][0] = 255;
          c[i][j][3] = 255;
        }
      }
    }

    for (int i = 0; i < w; i++) {
      for (int j = 0; j < h; j++) {
        if (val[i][j] == -1) {
          int pos_count = 0;
          for (int x = -1; x <= 1; x++) {
            for (int y = -1; y <= 1; y++) {
              if (x == 0 && y == 0) continue;
              int ix = i + x;
              int jy = j + y;
              if (ix < 0 || ix >= w || jy < 0 || jy >= h) continue;
              pos_count += (val[ix][jy] >= 0 ? 1 : 0);
            }
          }
          if (pos_count >= 2) c[i][j][0] = (GLuint)255, c[i][j][3] = 255;
        }
      }
    }
    /*
     for(int i = 0; i < w; i++) {
     for(int j = 0; j < h; j++) {
     if(valx[i][j] == -1) {
     int pos_count = 0;
     for(int x = -1; x <= 1; x++) {
     for(int y = -1; y <= 1; y++) {
     if(x == 0 && y == 0) continue;
     int ix = i+x;
     int jy = j+y;
     if(ix < 0 || ix >= w || jy < 0 || jy >= h) continue;
     pos_count += (valx[ix][jy] >= 0 ? 1 : 0);
     }
     }
     if(pos_count >= 2)
     c[i][j][1] = (GLuint)255, c[i][j][3] = 255;
     }
     }
     }


     for(int i = 0; i < w; i++) {
     for(int j = 0; j < h; j++) {
     if(valy[i][j] == -1) {
     int pos_count = 0;
     for(int x = -1; x <= 1; x++) {
     for(int y = -1; y <= 1; y++) {
     if(x == 0 && y == 0) continue;
     int ix = i+x;
     int jy = j+y;
     if(ix < 0 || ix >= w || jy < 0 || jy >= h) continue;
     pos_count += (valy[ix][jy] >= 0 ? 1 : 0);
     }
     }
     if(pos_count >= 2)
     c[i][j][2] = (GLuint)255, c[i][j][3] = 255;
     }
     }
     }
     */
  }

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE,
               c);

  printf("updated\n");
}

double r_x[6];
double r_y[6];
bool r_init = false;

void init_random_points() {
  int count = 0;

  for (double i = -0.75; i < 0.75; i += 0.5) {
    for (double j = -0.6; j < 0.6; j += 0.6) {
      // r_x[count] = i + ((rand() * 1.0)/ RAND_MAX) * 0.5;
      // r_y[count] = j + ((rand() * 1.0)/ RAND_MAX) * 0.6;
      r_x[count] = i + randomNumber(0.0, 0.5);
      r_y[count] = i + randomNumber(0.0, 0.6);
      count++;
    }
  }

  r_init = true;
}

vector<double> random_points(int count) {
  vector<double> ret;
  /*
   for(int i = 0; i < count; i++) {
   //double r = 0.2 + ((rand() * 0.8) / RAND_MAX);
   //double t = ((rand() * 2.0 * M_PI) / RAND_MAX);
   double r = randomNumber(0.2, 1);
   double t = randomNumber(0, 2.0 * M_PI);
   ret.push_back(r*cos(t));
   ret.push_back(r*sin(t));
   }
   */

  for (int i = 0; i < count; i++) {
    double x = -2 + ((rand() * 4.0) / RAND_MAX);
    double y = -2 + ((rand() * 4.0) / RAND_MAX);
    ret.push_back(x);
    ret.push_back(y);
  }

  return ret;
}

vector<int> choices;

vector<vector<double>> randomPoints;

void random_curve(int degree) {
  disable();

  // printf("degree: %d\n", degree);

  int num_points = ((degree + 1) * (degree + 2) / 2) - 1;

  // printf("num_points: %d\n", num_points);

  vector<double> p = random_points(num_points);

  randomPoints.push_back(p);

  double x[num_points];
  double y[num_points];
  for (int j = 0; j < num_points; j++) {
    x[j] = p[2 * j];
    y[j] = p[2 * j + 1];
  }

  double m[num_points * num_points];
  double b[num_points];
  int ipiv[num_points];

  for (int k = 0; k < num_points; k++) {
    int c = 0;

    for (int i = 0; i <= degree; i++) {
      for (int j = 0; j <= degree - i; j++) {
        if (i == 0 && j == 0)
          b[k] = -1;
        else {
          m[c * num_points + k] = pow(x[k], i) * pow(y[k], j);
          c++;
        }
      }
    }
  }

#ifdef OSX
  __CLPK_integer info;
#else
  int info;
#endif

#ifdef OSX
  __CLPK_integer n = num_points;
  __CLPK_integer nrhs = 1;
  __CLPK_integer lda = num_points;
  __CLPK_integer ipivm[num_points];
  __CLPK_integer ldb = num_points;
  dgesv_(&n, &nrhs, m, &lda, ipivm, b, &ldb, &info);
#else
  info = clapack_dgesv(CblasColMajor, num_points, 1, m, num_points, ipiv, b,
                       num_points);
#endif

  vector<Parameter> coef;
  vector<int> pow;

  int c = 0;

  double max_b = fabs(b[0]);

  for (int i = 1; i < num_points; i++) {
    if (fabs(b[i]) > max_b) max_b = fabs(b[i]);
  }

  for (int i = 0; i <= degree; i++) {
    for (int j = 0; j <= degree - i; j++) {
      if (i == 0 && j == 0) {
        coef.push_back(Parameter::input(1.0 / max_b));
        pow.push_back(0);
        pow.push_back(0);
        printf("%.16f 0 0 ", 1.0 / max_b);
      } else {
        coef.push_back(Parameter::input(b[c++] / max_b));
        pow.push_back(i);
        pow.push_back(j);
        printf("%.16f %d %d ", b[c - 1] / max_b, i, j);
      }
    }
  }

  printf("\n");

  polys.push_back(new Poly2<Parameter>(coef.size(), &coef[0], &pow[0]));

  // cout << "Created one curve size: " << coef.size() << endl << endl;

  enable();
}

void random_intersecting_ellipse() {
  disable();

  if (!r_init) {
    init_random_points();
  }

  if (choices.size() == 6) {
    enable();
    return;
  }

  int skip = (int)randomNumber(0, 6);
  while (std::find(choices.begin(), choices.end(), skip) != choices.end()) {
    skip = (int)(randomNumber(0, 6));
  }

  choices.push_back(skip);

  double x[5];
  double y[5];
  int i = 0;
  for (int j = 0; j < 6; j++) {
    if (j == skip) continue;
    x[i] = r_x[j];
    y[i] = r_y[j];
    i++;
  }

  double m[25];
  double b[5];
  int ipiv[5];

  // xy y^2 x y 1 = -x^2

  for (i = 0; i < 5; i++) {
    m[5 * i + 0] = x[i] * y[i];
    m[5 * i + 1] = y[i] * y[i];
    m[5 * i + 2] = x[i];
    m[5 * i + 3] = y[i];
    m[5 * i + 4] = 1;
    b[i] = -x[i] * x[i];
  }

#ifdef OSX
  __CLPK_integer info;
#else
  int info;
#endif

#ifdef OSX
  __CLPK_integer n = 5;
  __CLPK_integer nrhs = 1;
  __CLPK_integer lda = 5;
  __CLPK_integer ipivm[5];
  __CLPK_integer ldb = 5;
  dgesv_(&n, &nrhs, m, &lda, ipivm, b, &ldb, &info);
#else
  info = clapack_dgesv(CblasColMajor, 5, 1, m, 5, ipiv, b, 5);
#endif

  vector<Parameter> coef;
  vector<int> pow;

  coef.push_back(Parameter::input(1));
  pow.push_back(2);
  pow.push_back(0);
  coef.push_back(Parameter::input(b[0]));
  pow.push_back(1);
  pow.push_back(1);
  coef.push_back(Parameter::input(b[1]));
  pow.push_back(0);
  pow.push_back(2);
  coef.push_back(Parameter::input(b[2]));
  pow.push_back(1);
  pow.push_back(0);
  coef.push_back(Parameter::input(b[3]));
  pow.push_back(0);
  pow.push_back(1);
  coef.push_back(Parameter::input(b[4]));
  pow.push_back(0);
  pow.push_back(0);

  polys.push_back(new Poly2<Parameter>(coef.size(), &coef[0], &pow[0]));

  cout << "Created one curve size: " << coef.size() << endl << endl;

  enable();
}

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
  for (int i = 0; i < curves.size(); i++) {
    int cc = connected_components(curves[i], encasement);

    printf("curve %d: %d connected components\n", i, cc);
  }
}

void graphics_init_aux(int argc, char** argv) {
  /*
   vector<Parameter> coef;
   vector<int> pow;
   coef.push_back(Parameter::input(1));
   coef.push_back(Parameter::input(1));
   coef.push_back(Parameter::input(0));
   pow.push_back(1);
   pow.push_back(0);
   pow.push_back(0);
   pow.push_back(1);
   pow.push_back(0);
   pow.push_back(0);

   Poly2 f(3, &coef[0], &pow[0]);

   coef[0] = Parameter::input(-1);

   Poly2 g(3, &coef[0], &pow[0]);

   PV2<Parameter> p(Parameter::interval(-1000, 1000), Parameter::interval(-1000,
   1000));

   p = Root2::polish(p, f, g);
   */

  //  encasement = new Encasement(new InputPoint(-1, -1),
  //			      new InputPoint( 0,  0),
  //			      27);
  //
  encasement =
      new Encasement(new InputPoint(-2, -2), new InputPoint(2, 2), rand());

  w = glutGet(GLUT_WINDOW_WIDTH);
  h = glutGet(GLUT_WINDOW_HEIGHT);

  Poly2<Parameter>* ellipse;
  Ellipse_struct e;

  if (argc == 2) {
    useTexture = true;

    input = false;

    ifstream inFile(argv[1]);
    string line;

    while (getline(inFile, line)) {
      istringstream buf(line);
      istream_iterator<std::string> beg(buf), end;

      int count = 0;

      vector<Parameter> coef;
      vector<int> pow;

      for (; beg != end;) {
        double c = stod((*beg));
        beg++;
        int i = stoi((*beg));
        beg++;
        int j = stoi((*beg));
        beg++;
        // cout << c << "*(x^" << i << ")*(y^" << j << ") " << (beg != end ? "+
        // " : "");

        coef.push_back(Parameter::input(c));
        pow.push_back(i);
        pow.push_back(j);
      }
      // cout << " = 0; " << endl << endl;

      polys.push_back(new Poly2<Parameter>(coef.size(), &coef[0], &pow[0]));

      // cout << "Created one curve size: " << coef.size() << endl << endl;
    }

    // TEST
    for (int i = 0; i < polys.size(); i++) {
      for (int j = 0; j < polys[i]->size(); j++) {
        // cout << polys[i]->a[j].mid() << " " << polys[i]->m[2*j] << " " <<
        // polys[i]->m[2*j + 1] << endl;
      }
      // cout << endl;
    }
    // TEST

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glGenTextures(1, &textureName);
    glBindTexture(GL_TEXTURE_2D, textureName);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    updateTexture();

  } else if (argc == 6) {
    useTexture = true;

    input = false;

    ifstream inFile(argv[1]);
    string line;

    while (getline(inFile, line)) {
      istringstream buf(line);
      istream_iterator<std::string> beg(buf), end;

      int count = 0;

      vector<Parameter> coef;
      vector<int> pow;

      for (; beg != end;) {
        double c = stod((*beg));
        beg++;
        int i = stoi((*beg));
        beg++;
        int j = stoi((*beg));
        beg++;
        // cout << c << "*(x^" << i << ")*(y^" << j << ") " << (beg != end ? "+
        // " : "");

        coef.push_back(Parameter::input(c));
        pow.push_back(i);
        pow.push_back(j);
      }
      // cout << " = 0; " << endl << endl;

      polys.push_back(new Poly2<Parameter>(coef.size(), &coef[0], &pow[0]));

      // cout << "Created one curve size: " << coef.size() << endl << endl;
    }

    // TEST
    for (int i = 0; i < polys.size(); i++) {
      for (int j = 0; j < polys[i]->size(); j++) {
        // cout << polys[i]->a[j].mid() << " " << polys[i]->m[2*j] << " " <<
        // polys[i]->m[2*j + 1] << endl;
      }
      // cout << endl;
    }

    // TEST
    double llx = atof(argv[2]);
    double lly = atof(argv[3]);
    double urx = atof(argv[4]);
    double ury = atof(argv[5]);

    encasement =
        new Encasement(new InputPoint(atof(argv[2]), atof(argv[3])),
                       new InputPoint(atof(argv[4]), atof(argv[5])), 27);

    draw_scale = 1.0 / (fmax(fabs(llx - urx), fabs(lly - ury)));
    x_off = -1 * ((llx + urx) / 2);
    y_off = -1 * ((lly + ury) / 2);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glGenTextures(1, &textureName);
    glBindTexture(GL_TEXTURE_2D, textureName);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    updateTexture();

  } else if (argc == 13) {
    useTexture = true;

    input = false;

    Parameter* coefficients = new Parameter[6];
    int* powers = new int[12];

    for (int x = 0; x < 2; x++) {
      for (int i = 0; i < 6; i++) {
        coefficients[i] = Parameter::input(atof(argv[(i + x * 6) + 1]));
      }

      powers[0] = 2;
      powers[1] = 0;
      powers[2] = 1;
      powers[3] = 1;
      powers[4] = 0;
      powers[5] = 2;
      powers[6] = 1;
      powers[7] = 0;
      powers[8] = 0;
      powers[9] = 1;
      powers[10] = 0;
      powers[11] = 0;

      ellipse = new Poly2<Parameter>(6, coefficients, powers);

      Parameter a = ellipse->a[0];
      Parameter b = ellipse->a[1] / 2;
      Parameter c = ellipse->a[2];
      Parameter d = ellipse->a[3] / 2;
      Parameter f = ellipse->a[4] / 2;
      Parameter g = ellipse->a[5];

      e.xc = ((c * d - b * f) / (b * b - a * c)).mid();
      e.yc = ((a * f - b * d) / (b * b - a * c)).mid();

      e.a =
          ((2 * (a * f * f + c * d * d + g * b * b - 2 * b * d * f - a * c * g))
               .sqrt() /
           ((b * b - a * c) *
            (((a - c) * (a - c) + 4 * b * b).sqrt() - (a + c)))
               .sqrt())
              .mid();

      e.b =
          ((2 * (a * f * f + c * d * d + g * b * b - 2 * b * d * f - a * c * g))
               .sqrt() /
           ((b * b - a * c) *
            (-1 * ((a - c) * (a - c) + 4 * b * b).sqrt() - (a + c)))
               .sqrt())
              .mid();

      if (a < c) {
        e.phi = 1 / (2 * atan2((a - c).mid(), (2 * b).mid()));
      } else {
        e.phi = (1 / (2 * atan2((a - c).mid(), (2 * b).mid()))) + (0.5 * M_PI);
      }

      polys.push_back(ellipse);
      ellipses.push_back(e);
    }

    // TEST
    for (int i = 0; i < polys.size(); i++) {
      for (int j = 0; j < polys[i]->size(); j++) {
        cout << polys[i]->a[j].mid() << " " << polys[i]->m[2 * j] << " "
             << polys[i]->m[2 * j + 1] << endl;
      }
      cout << endl;
    }
    // TEST

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glGenTextures(1, &textureName);
    glBindTexture(GL_TEXTURE_2D, textureName);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    updateTexture();
  } else if (argc == 9) {
    useTexture = true;

    input = false;

    Parameter* coefficients = new Parameter[4];
    int* powers = new int[8];

    for (int x = 0; x < 2; x++) {
      for (int i = 0; i < 4; i++) {
        coefficients[i] = Parameter::input(atof(argv[(i + x * 4) + 1]));
      }

      powers[0] = 3;
      powers[1] = 0;
      powers[2] = 0;
      powers[3] = 2;
      powers[4] = 1;
      powers[5] = 0;
      powers[6] = 0;
      powers[7] = 0;

      ellipse = new Poly2<Parameter>(4, coefficients, powers);

      polys.push_back(ellipse);
    }

    // TEST
    for (int i = 0; i < polys.size(); i++) {
      for (int j = 0; j < polys[i]->size(); j++) {
        cout << polys[i]->a[j].mid() << " " << polys[i]->m[2 * j] << " "
             << polys[i]->m[2 * j + 1] << endl;
      }
      cout << endl;
    }
    // TEST

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glGenTextures(1, &textureName);
    glBindTexture(GL_TEXTURE_2D, textureName);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    updateTexture();

  } else {
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glGenTextures(1, &textureName);
    glBindTexture(GL_TEXTURE_2D, textureName);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    state = NONE;
    useTexture = true;
    input = false;
  }

  // encasement->STEPS = 1;
  // calculate();
}

void init_ellipse(double a, double b, double xc, double yc, double phi) {
  Parameter* coefficients = new Parameter[6];
  int* powers = new int[12];

  Poly2<Parameter>* ellipse;
  Ellipse_struct e;

  e.a = a;
  e.b = b;
  e.xc = xc;
  e.yc = yc;
  e.phi = phi;

  coefficients[0] = Parameter::input(1 / (e.a * e.a));
  coefficients[2] = Parameter::input(1 / (e.b * e.b));
  coefficients[5] = Parameter::input(1);

  powers[0] = 2;
  powers[1] = 0;
  powers[2] = 1;
  powers[3] = 1;
  powers[4] = 0;
  powers[5] = 2;
  powers[6] = 1;
  powers[7] = 0;
  powers[8] = 0;
  powers[9] = 1;
  powers[10] = 0;
  powers[11] = 0;

  // printf("created ellipse\n");

  ellipse = new Poly2<Parameter>(6, coefficients, powers);

  ellipses.push_back(e);
  polys.push_back(ellipse);
}

void random_ellipses(int n) {
  for (int i = 0; i < n; i++) {
    // double a = (rand() / ((double) RAND_MAX)) * 0.8 + 0.01;
    // double b = (rand() / ((double) RAND_MAX)) * 0.8 + 0.01;
    // double xc = (rand() / ((double) RAND_MAX)) - 0.5;
    // double yc = (rand() / ((double) RAND_MAX)) - 0.5;
    // double phi = (rand() / ((double) RAND_MAX)) * 2 * M_PI;
    double a = randomNumber(0.01, 0.801);
    double b = randomNumber(0.01, 0.801);
    double xc = randomNumber(-0.5, 0.5);
    double yc = randomNumber(-0.5, 0.5);
    double phi = randomNumber(0, 2 * M_PI);

    // printf("a: %g, b: %g, xc: %g, yc: %g, phi: %g\n", a, b, xc, yc, phi);
    init_ellipse(a, b, xc, yc, phi);
    updateEllipse();
  }
}

void menu_init() {
  glutCreateMenu(menu);
  glutAddMenuEntry("calculate", 1);
  glutAddMenuEntry("clear", 2);
  glutAddMenuEntry("intersections", 3);
  glutAddMenuEntry("iterate cells", 4);
  glutAddMenuEntry("iterate edges", 5);
  glutAddMenuEntry("random", 6);
  glutAddMenuEntry("STEP", 7);
  glutAddMenuEntry("Reload data", 8);
  glutAddMenuEntry("Report", 9);
  glutAddMenuEntry("Centroids", 10);
  glutAddMenuEntry("Update Texture", 11);
  glutAddMenuEntry("random intersecting", 12);
  glutAddMenuEntry("random curve", 13);
  glutAddMenuEntry("draw roots", 14);
  glutAddMenuEntry("draw triangulation", 15);
  glutAddMenuEntry("convert triangulation", 16);
  glutAddMenuEntry("toggle encasement", 17);
  glutAddMenuEntry("save image", 18);
#ifdef ENCASEMENT_GPU
  glutAddMenuEntry("gpu boxes", 19);
  glutAddMenuEntry("exit", 20);
#else
  glutAddMenuEntry("exit", 19);
#endif

  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void iterate_cell() {
  cell_index++;

  if (cell_index < encasement->faces.size()) {
    // printf("face: 0x%lx\n", (long)encasement->faces[cell_index]);
  }

  if (cell_index == encasement->faces.size()) cell_index = -1;
}

void iterate_edge() {
  if (encasement->faces.size() > 0) {
    if (g_index == -1) {
      g_index = 0;
      g_edge = encasement->faces[g_index]->e;
    } else {
      g_edge = g_edge->next;
      if (g_edge == encasement->faces[g_index]->e) {
        g_index++;
        if (g_index == encasement->faces.size()) {
          g_index = -1;
          return;
        }

        g_edge = encasement->faces[g_index]->e;
      }
    }
  }
}

void iterate_edge_back() {
  if (encasement->faces.size() > 0) {
    if (g_index == -1) {
      g_index = 0;
      g_edge = encasement->faces[g_index]->e;
    } else {
      g_edge = g_edge->prev;
      if (g_edge == encasement->faces[g_index]->e) {
        g_index++;
        if (g_index == encasement->faces.size()) {
          g_index = -1;
          return;
        }

        g_edge = encasement->faces[g_index]->e;
      }
    }
  }
}

void report() {
  int inter_count = 0;
  printf("\n%ld faces\n", encasement->faces.size());

  if (roots.size() == 0) {
    fsize = encasement->faces.size();
    roots = encasement->getRoots();
  }

  printf("%ld intersection points detected in total\n", roots.size());
}

void calculate() {
  // writeFile();

  // put all of the polys into ellipses and calculate
  /*
   */

  // TODO
  // make sure no re-perturbation
  // reseed properly
  for (int i = curves.size(); i < polys.size(); i++) {
    curves.push_back(new Ellipse(*polys[i]));
  }
  // calculate the encasement of the ellipse
  if (curves.size() > 0) {
    // printf("calculating up to %d steps\n", encasement->STEPS);
    auto t1 = std::chrono::high_resolution_clock::now();
    encasement->init();
    // encasement->calculate_single(curves[0]);
    encasement->calculate(curves, false);
    auto t2 = std::chrono::high_resolution_clock::now();
    calculated = true;
    // report();

    cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
                .count()
         << " ms\n";

    printf("intersections: %ld\n", encasement->getRoots().size());

    verify();
  }

  // exit(0);
}

void menu(int value) {
  int inter_count = 0;

  switch (value) {
    case 7:
      encasement->STEPS = ++STEPS;
    case 1:
      calculate();
      break;
    case 2:
      ellipses.clear();
      for (int i = 0; i < polys.size(); i++) {
        delete polys[i];
      }
      polys.clear();
      for (int i = 0; i < curves.size(); i++) {
        /*delete curves[i]*/;
      }
      curves.clear();
      encasement->init();
      encasement->STEPS = STEPS = 1;
      break;
    case 3:
      draw_intersections = !draw_intersections;
      break;
    case 4:
      iterate_cell();
      break;
    case 5:
      iterate_edge();
      break;
    case 6:
      random_ellipses(1);
      break;
    case 8:
      readFile();
      break;
    case 9:
      report();
      break;
    case 10:
      for (int i = 0; i < encasement->faces.size(); i++) {
        PV2D c = encasement->faces[i]->centroid();
        encasement->faces[i]->intersection_point = new InputPoint(c.x, c.y);
      }
      break;
    case 11:
      updateTexture();
      break;
    case 12:
      for (int i = 0; i < 6; i++) random_intersecting_ellipse();

      updateTexture();
      break;
    case 13:
      // random_curve((rand() % 3) + 2);
      // for(int i = 0; i < 4; i++)
      random_curve(((int)randomNumber(0, 3)) + 6);
      // random_curve(1);
      updateTexture();
      break;
    case 14:
      draw_roots = !draw_roots;
      break;
    case 15:
      draw_triangulation = !draw_triangulation;
      break;
    case 16:
      if (roots.size() == 0) {
        fsize = encasement->faces.size();
        roots = encasement->getRoots();
      }

      if (roots.size() == 0) {
        printf("No roots\n");
      }

      if (triangulation.size() == 0) {
        triangulation = encasement->getTriangulation(roots);
      }

      encasement->makeTriangulationEncasement(triangulation);

      break;
    case 17:
      draw_encasement = !draw_encasement;
      break;
    case 18:
      saveImage();
      break;

#ifdef ENCASEMENT_GPU
    case 19:
      draw_gpu_rects = !draw_gpu_rects;
    case 20:
#else
    case 19:
#endif
      ellipses.clear();
      for (int i = 0; i < polys.size(); i++) {
        delete polys[i];
      }
      polys.clear();

      encasement->clean();
      delete encasement;
      // printf("disabling parameter\n");
      disable();
      exit(0);
  }
  glutPostRedisplay();
}

void keyPressed(unsigned char key, int x, int y) {
  switch (key) {
    case 'l':
      if (roots.size() == 0) {
        fsize = encasement->faces.size();
        roots = encasement->getRoots();
      }
      root_index++;
      if (root_index == roots.size()) {
        root_index = -1;
        draw_scale = 1.0;
        x_off = 0.0;
        y_off = 0.0;
      } else {
        PV2<Parameter> r = roots[root_index]->getCurrentP();
        draw_scale = 1.0 / (fmax(r.x.intervalWidth(), r.y.intervalWidth()));
        x_off = -r.x.mid();
        y_off = -r.y.mid();
      }
      break;
    case '+':
      draw_scale *= 1.1;
      break;
    case '-':
      draw_scale /= 1.1;
      break;
    case '1':
      subdiv = max(2, subdiv >> 1);
      break;
    case '2':
      subdiv = min((1 << 10), subdiv << 1);
      break;
    case 'a':
      x_off += 0.1 / draw_scale;
      break;
    case 'd':
      x_off -= 0.1 / draw_scale;
      break;
    case 'w':
      y_off -= 0.1 / draw_scale;
      break;
    case 's':
      y_off += 0.1 / draw_scale;
      break;
    case 'c':
      iterate_cell();
      break;
    case 'n':
      iterate_edge();
      break;
    case 'p':
      iterate_edge_back();
      break;
    case 'i':
      if (encasement->STEPS < 1)
        encasement->STEPS = 1;
      else
        encasement->STEPS = ++STEPS;
      calculate();
      break;
    case 'o':
      encasement->STEPS = --STEPS;
      calculate();
      break;
#ifdef ENCASEMENT_GPU
    case 'g':
      draw_gpu_rects = !draw_gpu_rects;
      break;
#endif

    case 'r':
      draw_crit = !draw_crit;
      break;
    case 'q':
      dfirst = !dfirst;
      break;
  }

  glutPostRedisplay();
}

void mouse(int button, int buttonState, int x, int y) {
  if (buttonState == GLUT_DOWN) {
    if (calculated) {
      double w = glutGet(GLUT_WINDOW_WIDTH);
      double h = glutGet(GLUT_WINDOW_HEIGHT);

      double px = 2 * ((x - 0.5 * w) / w);
      double py = 2 * ((0.5 * h - y) / h);

      px = (px / draw_scale) - x_off;
      py = (py / draw_scale) - y_off;

      printf("click (%.2f %.2f)\n", px, py);

      click = point(x, y);
      printf("clicked face, finding arrangement face!\n");
    }

    glutPostRedisplay();
  }

  if (!input) return;

  if (buttonState == GLUT_DOWN) {
    PTR<Object<PV2>> p = point(x, y);

    switch (state) {
      case NONE:
        state = (State)((state + 1) % NUMBER_STATES);
        init_ellipse(0.5, 0.75, 0, 0, 0);
        return;
      case CENTER:
        ellipses.back().xc = p->getApprox(1.0).x.mid();
        ellipses.back().yc = p->getApprox(1.0).y.mid();
        break;
      case A:
        ellipses.back().a = p->getApprox(1.0).x.abs().mid();
        break;
      case B:
        ellipses.back().b = p->getApprox(1.0).y.abs().mid();
        break;
      case PHI:
        ellipses.back().phi = M_PI * (p->getApprox(1.0).y.mid());
        break;
      default:
        break;
    }

    state = (State)((state + 1) % NUMBER_STATES);
    updateEllipse();
    glutPostRedisplay();
  }
}

void mouseMoved(int x, int y) {
  double w = glutGet(GLUT_WINDOW_WIDTH);
  double h = glutGet(GLUT_WINDOW_HEIGHT);

  double px = 2 * ((x - 0.5 * w) / w);
  double py = 2 * ((0.5 * h - y) / h);

  px = (px / draw_scale) - x_off;
  py = (py / draw_scale) - y_off;
  char s[200];

  if (polys.size() >= 2) {
    double fv =
        polys[0]->value(Parameter::constant(px), Parameter::constant(py)).mid();
    double gv =
        polys[1]->value(Parameter::constant(px), Parameter::constant(py)).mid();

    sprintf(s, "encasement (%f, %f)    f = %.3f    g = %.3f", px, py, fv, gv);

  } else {
    sprintf(s, "encasement (%f %f)", px, py);
  }

  glutSetWindowTitle(s);

  if (!input) return;

  PTR<Object<PV2>> p = point(x, y);

  switch (state) {
    case NONE:
      return;
    case CENTER:
      ellipses.back().xc = p->getApprox(1.0).x.mid();
      ellipses.back().yc = p->getApprox(1.0).y.mid();
      break;
    case A:
      ellipses.back().a = p->getApprox(1.0).x.abs().mid();
      break;
    case B:
      ellipses.back().b = p->getApprox(1.0).y.abs().mid();
      break;
    case PHI:
      ellipses.back().phi = M_PI * (p->getApprox(1.0).y.mid());
      break;
    default:
      break;
  }

  updateEllipse();

  glutPostRedisplay();
}

void updateEllipse() {
  Ellipse_struct e = ellipses.back();

  double cos_phi = cos(e.phi);
  double sin_phi = sin(e.phi);

  double A =
      ((cos_phi * cos_phi) / (e.a * e.a)) + ((sin_phi * sin_phi) / (e.b * e.b));

  double B = -2 * cos_phi * sin_phi * ((1 / (e.a * e.a)) - (1 / (e.b * e.b)));

  double C =
      ((cos_phi * cos_phi) / (e.b * e.b)) + ((sin_phi * sin_phi) / (e.a * e.a));

  double D = -(2 * A * e.xc + e.yc * B);

  double E = -(2 * C * e.yc + B * e.xc);

  double F = A * e.xc * e.xc + B * e.xc * e.yc + C * e.yc * e.yc - 1;

  polys.back()->a[0] = Parameter::input(A);
  polys.back()->a[1] = Parameter::input(B);
  polys.back()->a[2] = Parameter::input(C);
  polys.back()->a[3] = Parameter::input(D);
  polys.back()->a[4] = Parameter::input(E);
  polys.back()->a[5] = Parameter::input(F);

  // printf("%.3f %.3f %.3f %.3f %.3f %.3f\n", A, B, C, D, E, F);
}

void drawEllipse(int i) {
  Ellipse_struct e = ellipses[i];

  double x;
  double y;

  vector<PTR<Object<PV2>>> points;

  double theta = 0;

  while (theta < 360) {
    x = e.xc + e.a * cos((theta / 180) * M_PI) * cos(-e.phi) -
        e.b * sin((theta / 180) * M_PI) * sin(-e.phi);
    y = e.yc + e.b * sin((theta / 180) * M_PI) * cos(-e.phi) +
        e.a * cos((theta / 180) * M_PI) * sin(-e.phi);

    points.push_back(new InputPoint(x, y));

    theta += 0.1;
  }

  draw_loop(points);

  points.clear();
}

void drawEllipse() {
  if (!useTexture) {
    for (int i = 0; i < ellipses.size(); i++) {
      drawEllipse(i);
    }
  } else {
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, textureName);

    double x1 = (-1 + x_off) * draw_scale;
    double x2 = (1 + x_off) * draw_scale;
    double y1 = (-1 + y_off) * draw_scale;
    double y2 = (1 + y_off) * draw_scale;

    // Draw a textured quad
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0);
    glVertex3f(x1, y1, 0);
    glTexCoord2f(1, 0);
    glVertex3f(x2, y1, 0);
    glTexCoord2f(1, 1);
    glVertex3f(x2, y2, 0);
    glTexCoord2f(0, 1);
    glVertex3f(x1, y2, 0);
    glEnd();

    glDisable(GL_TEXTURE_2D);
  }
}

void drawEncasement() {
  bool red = true;

  for (int i = 0; i < encasement->faces.size(); i++) {
    PTR<Face> c = encasement->faces[i];

    //    printf("cell[%d]\n", i);

    PTR<Edge> e = c->e;

    if (c->intersection_point != NULL) {
      glColor3f(1.0f, 1.0f, 0.0f);
      std::vector<PTR<Object<PV2>>> points;
      points.push_back(c->intersection_point);
      draw_points(points);
      points.clear();
    }

    PTR<Object<PV2>> o = e->tail->p;

    e = e->next;

    float f1, f2, f3;
    // f1 = (float)(rand() / (1.0 * RAND_MAX));
    // f2 = (float)(rand() / (1.0 * RAND_MAX));
    // f3 = (float)(rand() / (1.0 * RAND_MAX));

    // f1 = (float)(randomNumber(0, 1));
    // f2 = (float)(randomNumber(0, 1));
    // f3 = (float)(randomNumber(0, 1));

    f3 = 0.75f;

    // glColor4f(f1, f2, f3, 0.2f);
    glColor4f(0.0f, 0.0f, f3, 0.2f);

    do {
      /*
       if(red) {
       glColor4f(1.0f, 0.0f, 0.0f, 0.1f);
       } else {
       glColor4f(0.0f, 0.0f, 1.0f, 0.1f);
       }
       */

      draw_triangle(o, e->tail->p, e->head->p);

      e = e->next;

    } while (e->next != c->e);

    e = c->e;

    do {
      glColor3f(0.0f, 0.0f, 0.0f);

      draw_line(e->tail->p, e->head->p);

      if (e->intersection != nullptr) {
        glColor3f(0.0f, 1.0f, 0.0f);
        PTR<Object<PV2>> ip = e->intersection_point();
        std::vector<PTR<Object<PV2>>> points;
        points.push_back(ip);
        draw_points(points);
        /*delete ip*/;
        points.clear();

        //        printf("\tcurve: %ld\n", (long int) e->getPoly());
      }

      e = e->next;

    } while (e != c->e);

    red = !red;
  }
}

void drawIntersections() {
  for (int i = 0; i < encasement->faces.size(); i++) {
    if (encasement->faces[i]->intersection_point) {
      glColor3f(1.0f, 1.0f, 0.0f);
      std::vector<PTR<Object<PV2>>> points;
      points.push_back(encasement->faces[i]->intersection_point);
      draw_points(points);
      points.clear();
      // curve_trace2(encasement->faces[i]);
    } else if (encasement->curves.size() > 0) {
      // curve_trace2(encasement->faces[i]);
    }
  }
}

void drawRandomPoints() {
  for (int i = 0; i < randomPoints.size(); i++) {
    vector<double> p = randomPoints[i];
    vector<PTR<Object<PV2>>> pp;
    for (int j = 0; j < p.size() / 2; j++)
      pp.push_back(new InputPoint(p[2 * j], p[2 * j + 1]));
    // printf("pp_size: %ld\n", pp.size());
    glColor3f(1.0f, 0.0f, 1.0f);
    draw_points(pp);
    glColor3f(1.0f, 0.0f, 0.0f);
  }

  /*
   //TEST



   double d[] = {0.29844, 0.411647, 0.65085, 0.41779, 0.804897, 0.382075,
   0.678415, 0.287163, 0.35357, 0.150391, 0.0287249, 0.0136203, -0.0977568,
   -0.0812922, 0.05629, -0.117008, 0.4087, -0.110865};

   vector<double> v(d, d + sizeof(d)/sizeof(d[0]));
   vector< PTR<Object<PV2>> > pp;
   for(int j = 0; j < v.size()/2; j++) pp.push_back(new InputPoint(v[2*j],
   v[2*j+1])); glColor3f(1.0f, 0.0f, 1.0f); draw_points(pp); glColor3f(1.0f,
   0.0f, 0.0f);

   //TEST
   */
}

void debugDrawCrit(DrawItem item) {
  glutSetWindowTitle("debug face");
  float white[3] = {1.0f, 1.0f, 1.0f};
  clear_screen(white);

  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

  glLineWidth(3.0f);

  for (int i = 0; i < item.pvs.size(); i++) {
    glColor4f(item.colors[i].r, item.colors[i].g, item.colors[i].b, 1.0f);
    switch (item.types[i]) {
      case POINT:
        draw_points(item.pvs[i]);
        break;
      case LINE:
        draw_line(item.pvs[i][0], item.pvs[i][1]);
        break;
      case VECTOR:
        draw_vector(item.pvs[i][0], item.pvs[i][1]);
        break;
      case TRI:
        draw_triangle(item.pvs[i][0], item.pvs[i][1], item.pvs[i][2]);
        break;
    }
  }

  glLineWidth(linewidth);

  glutSwapBuffers();

  saveImage();
}

double debugDrawFacePoint(PTR<Face> face, PTR<Object<Poly2>> f,
                          PTR<Object<Poly2>> g, double width, double dx,
                          double dy, DrawItem item, bool draw) {
  if (!draw) return width;

  glutSetWindowTitle("debug face");
  float white[4] = {1.0f, 1.0f, 1.0f, 1.0f};
  clear_screen(white);

  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

  double ds_s = draw_scale;
  double xo_s = x_off;
  double yo_s = y_off;
  double ps_s = post_scale;

  post_scale = 1.0;
  draw_scale = 1.0 / width;
  x_off = -dx;
  y_off = -dy;

  curve_trace2(face);

  PTR<Edge> e = face->e;
  vector<PTR<Object<PV2>>> inters;
  do {
    if (e->intersection) inters.push_back(e->intersection_point());
  } while ((e = e->next) != face->e);

  glLineWidth(3.0f);

  for (int i = 0; i < item.pvs.size(); i++) {
    glColor4f(item.colors[i].r, item.colors[i].g, item.colors[i].b, 1.0f);
    switch (item.types[i]) {
      case POINT:
        draw_points(item.pvs[i]);
        break;
      case LINE:
        draw_line(item.pvs[i][0], item.pvs[i][1]);
        break;
      case VECTOR:
        draw_vector(item.pvs[i][0], item.pvs[i][1]);
        break;
      case TRI:
        draw_triangle(item.pvs[i][0], item.pvs[i][1], item.pvs[i][2]);
        break;
    }
  }

  glLineWidth(linewidth);

  glutSwapBuffers();

  draw_scale = ds_s;
  x_off = xo_s;
  y_off = yo_s;
  post_scale = ps_s;

  string s;
  cin >> s;

  if (s == "s") {
    cout << "saving image" << endl;
    saveImage();
  } else if (s == "+") {
    return -width / 2;
  } else if (s == "-") {
    return -width * 2;
  }

  return width;
}

void debugDrawFaces(vector<PTR<Face>> faces, PTR<Object<Poly2>> f,
                    PTR<Object<Poly2>> g, double zoom, double dx, double dy,
                    DrawItem item, bool draw, bool hold) {
  if (!draw) return;

  glutSetWindowTitle("debug face");
  float white[3] = {1.0f, 1.0f, 1.0f};
  float clr_alpha[4] = {0.0f, 0.0f, 0.0f, 0.0f};
  if (!hold) {
    clear_screen(white);
  }

  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

  double ds_s = draw_scale;
  double xo_s = x_off;
  double yo_s = y_off;
  double ps_s = post_scale;

  vector<PTR<Object<Scalar>>> bounds = faces[0]->bounds();

  Parameter xmin = bounds[0]->getApprox(1.0);
  Parameter ymin = bounds[1]->getApprox(1.0);
  Parameter xmax = bounds[2]->getApprox(1.0);
  Parameter ymax = bounds[3]->getApprox(1.0);

  Parameter x = xmin.interval(xmax);
  Parameter y = ymin.interval(ymax);

  post_scale = 1.0;
  draw_scale = 1.0 / (fmax(x.intervalWidth(), y.intervalWidth()) * zoom);
  // draw_scale = 1.0;
  x_off = -x.mid() + dx * x.intervalWidth();
  y_off = -y.mid() + dy * y.intervalWidth();

  for (int i = 1; i < faces.size(); i++) {
    curve_trace2(faces[i]);
  }

  PTR<Edge> e = faces[0]->e;
  vector<PTR<Object<PV2>>> inters;
  do {
    if (e->intersection) inters.push_back(e->intersection_point());
  } while ((e = e->next) != faces[0]->e);

  glLineWidth(3.0f);

  for (int i = 0; i < item.pvs.size(); i++) {
    glColor4f(item.colors[i].r, item.colors[i].g, item.colors[i].b, 1.0f);
    switch (item.types[i]) {
      case POINT:
        draw_points(item.pvs[i]);
        break;
      case LINE:
        draw_line(item.pvs[i][0], item.pvs[i][1]);
        break;
      case VECTOR:
        draw_vector(item.pvs[i][0], item.pvs[i][1]);
        break;
      case TRI:
        draw_triangle(item.pvs[i][0], item.pvs[i][1], item.pvs[i][2]);
        break;
    }
  }

  glLineWidth(linewidth);

  glutSwapBuffers();

  draw_scale = ds_s;
  x_off = xo_s;
  y_off = yo_s;
  post_scale = ps_s;

  string s;
  cin >> s;

  if (s == "s") {
    cout << "saving image" << endl;
    saveImage();
  }
}
void debugDrawFace(PTR<Face> face, PTR<Object<Poly2>> f, PTR<Object<Poly2>> g,
                   double zoom, double dx, double dy, DrawItem item, bool draw,
                   bool hold) {
  if (!draw) return;

  glutSetWindowTitle("debug face");
  float white[3] = {1.0f, 1.0f, 1.0f};
  if (!hold) clear_screen(white);

  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

  double ds_s = draw_scale;
  double xo_s = x_off;
  double yo_s = y_off;
  double ps_s = post_scale;

  vector<PTR<Object<Scalar>>> bounds = face->bounds();

  Parameter xmin = bounds[0]->getApprox(1.0);
  Parameter ymin = bounds[1]->getApprox(1.0);
  Parameter xmax = bounds[2]->getApprox(1.0);
  Parameter ymax = bounds[3]->getApprox(1.0);

  Parameter x = xmin.interval(xmax);
  Parameter y = ymin.interval(ymax);

  post_scale = 1.0;
  draw_scale = 1.0 / (fmax(x.intervalWidth(), y.intervalWidth()) * zoom);
  // draw_scale = 1.0;
  x_off = -x.mid() + dx * x.intervalWidth();
  y_off = -y.mid() + dy * y.intervalWidth();

  curve_trace2(face);

  PTR<Edge> e = face->e;
  vector<PTR<Object<PV2>>> inters;
  do {
    if (e->intersection) inters.push_back(e->intersection_point());
  } while ((e = e->next) != face->e);

  glLineWidth(3.0f);

  for (int i = 0; i < item.pvs.size(); i++) {
    glColor4f(item.colors[i].r, item.colors[i].g, item.colors[i].b, 1.0f);
    switch (item.types[i]) {
      case POINT:
        draw_points(item.pvs[i]);
        break;
      case LINE:
        draw_line(item.pvs[i][0], item.pvs[i][1]);
        break;
      case VECTOR:
        draw_vector(item.pvs[i][0], item.pvs[i][1]);
        break;
      case TRI:
        draw_triangle(item.pvs[i][0], item.pvs[i][1], item.pvs[i][2]);
        break;
    }
  }

  glLineWidth(linewidth);

  glutSwapBuffers();

  draw_scale = ds_s;
  x_off = xo_s;
  y_off = yo_s;
  post_scale = ps_s;

  string s;
  cin >> s;

  if (s == "s") {
    cout << "saving image" << endl;
    saveImage();
  }
}

// scale is first pv max(width, height)
// translation is center of first pv
void debugDraw(vector<PV2<Parameter>> pvs, vector<PTR<Object<Poly2>>> fs,
               bool draw = false) {
  if (!draw) return;

  glutSetWindowTitle("encasement");
  float white[3] = {1.0f, 1.0f, 1.0f};
  clear_screen(white);

  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

  double ds_s = draw_scale;
  double xo_s = x_off;
  double yo_s = y_off;
  double ps_s = post_scale;

  PV2<Parameter> b = pvs[0];
  post_scale = 1.0;
  draw_scale = 1.0 / fmax(b.x.intervalWidth(), b.y.intervalWidth());
  x_off = -b.x.mid();
  y_off = -b.y.mid();

  for (int i = 0; i < pvs.size(); i++) {
    PTR<Face> f = Face::make(pvs[i]);
    for (int j = 0; j < fs.size(); j++) f->intersect(fs[j]);

    curve_trace2(f);
  }

  glutSwapBuffers();

  draw_scale = ds_s;
  x_off = xo_s;
  y_off = yo_s;
  post_scale = ps_s;
}

void drawCell(int i, bool inter) {
  bool red = true;

  PTR<Face> c = encasement->faces[i];

  // printf("\ncontains %ld curves\n", c.curves.size());

  //    printf("cell[%d]\n", i);

  PTR<Edge> e = c->e;

  if (c->intersection_point != nullptr) {
    glColor3f(1.0f, 1.0f, 0.0f);
    std::vector<PTR<Object<PV2>>> points;
    points.push_back(c->intersection_point);
    draw_points(points);
    points.clear();
  }

  do {
    glColor3f(1.0f, 0.0f, 0.0f);
    if (e->poly && !useTexture) {
      drawEllipse((static_cast<Ellipse*>((Object<Poly2>*)e->poly))->id);
    } else {
      drawEllipse();
    }

    if (red) {
      glColor3f(1.0f, 0.0f, 0.0f);
    } else {
      glColor3f(0.0f, 0.0f, 1.0f);
    }

    draw_line(e->tail->p, e->head->p);

    if (inter && e->intersection != nullptr) {
      glColor3f(0.0f, 1.0f, 0.0f);
      PTR<Object<PV2>> ip = e->intersection_point();
      std::vector<PTR<Object<PV2>>> points;
      points.push_back(ip);
      draw_points(points);
      /*delete ip*/;
      points.clear();
    }

    e = e->next;

    red = !red;

  } while (e != c->e);
}

void drawEdge() {
  glColor3f(0.0f, 0.0f, 0.0f);

  // draw_line(g_edge->tail->p, g_edge->head->p);

  PTR<Object<PV2>> p = new LinePoint(g_edge->line, g_edge->t1);
  PTR<Object<PV2>> q = new LinePoint(g_edge->line, g_edge->t2);

  // printf("g_edge->line: 0x%lx\n", (long)g_edge);
  // printf("p: (%g, %g) q: (%g, %g)\n", p->get().x.mid(), p->get().y.mid(),
  // q->get().x.mid(), q->get().y.mid());
  draw_line(p, q);

  printf("g_edge: t1: %g t2: %g\n", g_edge->t1->getApprox().x.mid(),
         g_edge->t2->getApprox().x.mid());
  // printf("g_edge->tail: (%g, %g)\n", g_edge->tail->p->get().x.mid(),
  // g_edge->tail->p->get().y.mid()); printf("g_edge->head: (%g, %g)\n",
  // g_edge->head->p->get().x.mid(), g_edge->head->p->get().y.mid());

  if (g_edge->intersection != nullptr) {
    printf(" ip: %g\n", g_edge->intersection->getApprox().x.mid());
    glColor3f(0.0f, 1.0f, 0.0f);
    PTR<Object<PV2>> ip = g_edge->intersection_point();
    std::vector<PTR<Object<PV2>>> points;
    points.push_back(ip);
    draw_points(points);
    /*delete ip*/;
    points.clear();
  }

  fflush(stdout);
}

std::vector<double> bounds_list;

vector<vector<Rectangle*>> vs;
vector<std::map<int, std::set<Rectangle*>>> maps;

void drawRawIntersectionPoints() {
  // if(encasement->faces.size() != 1) return;
  if (polys.size() < 2) return;

  double xl, xu, yl, yu;

  double wi = 128;
  double hi = 128;

  if (bounds_list.size() == 0) {
    xl = encasement->faces[0]->e->tail->p->getApprox(1.0).x.lb();
    xu = encasement->faces[0]->e->tail->p->getApprox(1.0).x.ub();
    yl = encasement->faces[0]->e->tail->p->getApprox(1.0).y.lb();
    yu = encasement->faces[0]->e->tail->p->getApprox(1.0).y.ub();

    double xl2, xu2, yl2, yu2;

    printf("faces size: %ld\n", encasement->faces.size());

    PTR<Edge> e = encasement->faces[0]->e->next;
    while (e != encasement->faces[0]->e) {
      xl2 = e->tail->p->getApprox(1.0).x.lb();
      xu2 = e->tail->p->getApprox(1.0).x.ub();
      yl2 = e->tail->p->getApprox(1.0).y.lb();
      yu2 = e->tail->p->getApprox(1.0).y.ub();

      xl = min(xl, xl2);
      xu = max(xu, xu2);
      yl = min(yl, yl2);
      yu = max(yu, yu2);

      e = e->next;
    }

    xl -= Parameter::delta;
    xu += Parameter::delta;
    yl -= Parameter::delta;
    yu += Parameter::delta;

    bounds_list.push_back(xl);
    bounds_list.push_back(xu);
    bounds_list.push_back(yl);
    bounds_list.push_back(yu);
  }

  xl = bounds_list[0];
  xu = bounds_list[1];
  yl = bounds_list[2];
  yu = bounds_list[3];

  for (int i = curves.size(); i < polys.size(); i++) {
    curves.push_back(new Ellipse(*polys[i]));
  }

  PTR<Object<Poly2>> f = curves[0];
  PTR<Object<Poly2>> g = curves[1];

  double dx = (xu - xl) / wi;
  double dy = (yu - yl) / hi;

  for (int i = 0; i < wi; i++) {
    for (int j = 0; j < hi; j++) {
      double pxl = xl + (i)*dx;
      double pxu = xl + (i + 1) * dx;
      double pyl = yl + (j)*dy;
      double pyu = yl + (j + 1) * dy;

      // printf("f([%f, %f]x[%f, %f]\n", pxl, pxu, pyl, pyu);

      PTR<Object<PV2>> p = new ParameterPoint(
          new DoubleInputInterval(pxl, pxu), new DoubleInputInterval(pyl, pyu));

      // int sf = (new PolyScalar(f, p))->get().sign(false);
      // int sg = (new PolyScalar(g, p))->get().sign(false);

      int sf = quad_approx_sign(p, f);
      int sg = quad_approx_sign(p, g);

      if (sf == 0 && sg == 0) {
        glColor4f(0.8f, 0.0f, 0.0f, 0.5f);

        PTR<Object<PV2>> p = new InputPoint(PV2<Parameter>::constant(pxl, pyl));
        PTR<Object<PV2>> q = new InputPoint(PV2<Parameter>::constant(pxu, pyl));
        PTR<Object<PV2>> r = new InputPoint(PV2<Parameter>::constant(pxu, pyu));
        PTR<Object<PV2>> s = new InputPoint(PV2<Parameter>::constant(pxl, pyu));

        draw_triangle(p, q, r);
        draw_triangle(p, r, s);

        glColor4f(0.0f, 0.0f, 0.0f, 1.0f);

        draw_line(p, q);
        draw_line(q, r);
        draw_line(r, s);
        draw_line(s, p);
      }
    }
  }
}

void drawCriticalPoints() {
  // drawRawIntersectionPoints();
  // return;

  // if(encasement->faces.size() != 1) return;
  if (polys.size() < 1) return;

  double xl, xu, yl, yu;

  if (bounds_list.size() == 0) {
    xl = encasement->faces[0]->e->tail->p->getApprox(1.0).x.lb();
    xu = encasement->faces[0]->e->tail->p->getApprox(1.0).x.ub();
    yl = encasement->faces[0]->e->tail->p->getApprox(1.0).y.lb();
    yu = encasement->faces[0]->e->tail->p->getApprox(1.0).y.ub();

    double xl2, xu2, yl2, yu2;

    printf("faces size: %ld\n", encasement->faces.size());

    PTR<Edge> e = encasement->faces[0]->e->next;
    while (e != encasement->faces[0]->e) {
      xl2 = e->tail->p->getApprox(1.0).x.lb();
      xu2 = e->tail->p->getApprox(1.0).x.ub();
      yl2 = e->tail->p->getApprox(1.0).y.lb();
      yu2 = e->tail->p->getApprox(1.0).y.ub();

      xl = min(xl, xl2);
      xu = max(xu, xu2);
      yl = min(yl, yl2);
      yu = max(yu, yu2);

      e = e->next;
    }

    xl -= Parameter::delta;
    xu += Parameter::delta;
    yl -= Parameter::delta;
    yu += Parameter::delta;

    bounds_list.push_back(xl);
    bounds_list.push_back(xu);
    bounds_list.push_back(yl);
    bounds_list.push_back(yu);
  }

  xl = bounds_list[0];
  xu = bounds_list[1];
  yl = bounds_list[2];
  yu = bounds_list[3];

  printf("bounds: [%f %f] [%f %f]\n", xl, xu, yl, yu);

  printf("polys size: %ld\n", polys.size());

  static unordered_map<Object<Poly2>*, Color> color_map;

  for (int i = 0; i < polys.size(); i++) {
    // for(int i = 2; i < 3; i++) {

    vector<Rectangle*> v;
    if (vs.size() <= i) {
      v = critical_points(
          new Ellipse(*polys[i]),
          new Rectangle(PV2<Parameter>(Parameter::interval(xl, xu),
                                       Parameter::interval(yl, yu))));
      vs.push_back(v);
    } else {
      v = vs[i];
    }

    printf("v size: %ld\n", v.size());

    std::map<int, std::set<Rectangle*>> mymap;
    if (maps.size() <= i) {
      mymap = connected_components(v);
    } else {
      mymap = maps[i];
      maps.push_back(mymap);
    }

    int mcount = 0;
    float r = (float)randomNumber(0, 1);
    float g = (float)randomNumber(0, 1);
    float b = (float)randomNumber(0, 1);
    Color c;
    if (color_map.find((Object<Poly2>*)polys[i]) == color_map.end()) {
      if (i == 0)
        c = {0.0f, 1.0f, 0.0f};
      else if (i == 1)
        c = {1.0f, 0.0f, 0.0f};
      else if (i == 2)
        c = {0.0f, 0.0f, 1.0f};
      else
        c = {(float)randomNumber(0.5, 1), (float)randomNumber(0.5, 1),
             (float)randomNumber(0.5, 1)};

      color_map[(Object<Poly2>*)polys[i]] = c;
    }
    c = color_map[(Object<Poly2>*)polys[i]];
    r = (float)c.r;
    g = (float)c.g;
    b = (float)c.b;

    for (map<int, set<Rectangle*>>::iterator it = mymap.begin();
         it != mymap.end(); it++) {
      set<Rectangle*> myset = (*it).second;
      mcount++;

      // printf("%ld connected components\n", set.size());

      for (set<Rectangle*>::iterator itt = myset.begin(); itt != myset.end();
           itt++) {
        // float r = (i * 1.0) / polys.size();
        // float g = 0.34f;
        // float b = 0.23f;

        // for(int j = 0; j < v.size(); j++) {

        glColor4f(r, g, b, 0.5f);

        // Rectangle * re = v[j];

        Rectangle* re = *itt;

        PTR<Object<PV2>> p = new InputPoint(
            PV2<Parameter>::constant(re->box.x.lb(), re->box.y.lb()));
        PTR<Object<PV2>> q = new InputPoint(
            PV2<Parameter>::constant(re->box.x.ub(), re->box.y.lb()));
        PTR<Object<PV2>> r = new InputPoint(
            PV2<Parameter>::constant(re->box.x.ub(), re->box.y.ub()));
        PTR<Object<PV2>> s = new InputPoint(
            PV2<Parameter>::constant(re->box.x.lb(), re->box.y.ub()));

        draw_triangle(p, q, r);
        draw_triangle(p, r, s);

        glColor4f(0.0f, 0.0f, 0.0f, 1.0f);

        draw_line(p, q);
        draw_line(q, r);
        draw_line(r, s);
        draw_line(s, p);
      }
    }

    printf("%d connected components\n", mcount);
  }
}

#ifdef ENCASEMENT_GPU

void drawGPURects() {
  // drawRawIntersectionPoints();
  // return;

  // if(encasement->faces.size() != 1) return;
  if (polys.size() < 1) return;

  double xl, xu, yl, yu;

  if (bounds_list.size() == 0) {
    xl = encasement->faces[0]->e->tail->p->getApprox(1.0).x.lb();
    xu = encasement->faces[0]->e->tail->p->getApprox(1.0).x.ub();
    yl = encasement->faces[0]->e->tail->p->getApprox(1.0).y.lb();
    yu = encasement->faces[0]->e->tail->p->getApprox(1.0).y.ub();

    double xl2, xu2, yl2, yu2;

    printf("faces size: %ld\n", encasement->faces.size());

    PTR<Edge> e = encasement->faces[0]->e->next;
    while (e != encasement->faces[0]->e) {
      xl2 = e->tail->p->getApprox(1.0).x.lb();
      xu2 = e->tail->p->getApprox(1.0).x.ub();
      yl2 = e->tail->p->getApprox(1.0).y.lb();
      yu2 = e->tail->p->getApprox(1.0).y.ub();

      xl = min(xl, xl2);
      xu = max(xu, xu2);
      yl = min(yl, yl2);
      yu = max(yu, yu2);

      e = e->next;
    }

    xl -= Parameter::delta;
    xu += Parameter::delta;
    yl -= Parameter::delta;
    yu += Parameter::delta;

    bounds_list.push_back(xl);
    bounds_list.push_back(xu);
    bounds_list.push_back(yl);
    bounds_list.push_back(yu);
  }

  xl = bounds_list[0];
  xu = bounds_list[1];
  yl = bounds_list[2];
  yu = bounds_list[3];

  printf("bounds: [%f %f] [%f %f]\n", xl, xu, yl, yu);

  printf("polys size: %ld\n", polys.size());

  for (int i = curves.size(); i < polys.size(); i++) {
    curves.push_back(new Ellipse(*polys[i]));
  }

  vector<Rectangle*> gpu_all_rects =
      get_all_rects(curves[0], curves[1], xl, yl, xu, yu);
  vector<Rectangle*> gpu_regions =
      get_regions(curves[0], curves[1], xl, yl, xu, yu);

  float cr = 0.8f;
  float cg = 0.8f;
  float cb = 0.0f;

  for (int i = 0; i < gpu_all_rects.size(); i++) {
    Rectangle* re = gpu_all_rects[i];

    PTR<Object<PV2>> p = new InputPoint(
        PV2<Parameter>::constant(re->box.x.lb(), re->box.y.lb()));
    PTR<Object<PV2>> q = new InputPoint(
        PV2<Parameter>::constant(re->box.x.ub(), re->box.y.lb()));
    PTR<Object<PV2>> r = new InputPoint(
        PV2<Parameter>::constant(re->box.x.ub(), re->box.y.ub()));
    PTR<Object<PV2>> s = new InputPoint(
        PV2<Parameter>::constant(re->box.x.lb(), re->box.y.ub()));

    glColor4f(cr, cg, cb, 0.7f);

    draw_triangle(p, q, r);
    draw_triangle(p, r, s);

    glColor4f(0.0f, 0.0f, 0.0f, 1.0f);

    draw_line(p, q);
    draw_line(q, r);
    draw_line(r, s);
    draw_line(s, p);
  }

  for (int i = 0; i < gpu_regions.size(); i++) {
    Rectangle* re = gpu_regions[i];

    PTR<Object<PV2>> p = new InputPoint(
        PV2<Parameter>::constant(re->box.x.lb(), re->box.y.lb()));
    PTR<Object<PV2>> q = new InputPoint(
        PV2<Parameter>::constant(re->box.x.ub(), re->box.y.lb()));
    PTR<Object<PV2>> r = new InputPoint(
        PV2<Parameter>::constant(re->box.x.ub(), re->box.y.ub()));
    PTR<Object<PV2>> s = new InputPoint(
        PV2<Parameter>::constant(re->box.x.lb(), re->box.y.ub()));

    glColor4f(1.0f, 0.0f, 1.0f, 1.0f);
    glLineWidth(2.0f);

    draw_line(p, q);
    draw_line(q, r);
    draw_line(r, s);
    draw_line(s, p);
  }
}

#endif

PV2<Parameter> transform(PV2<Parameter> p) {
  return PV2<Parameter>(post_scale * ((p.x / draw_scale) - x_off),
                        post_scale * ((p.y / draw_scale) - y_off));
}

// does edge ab cross cd?
bool crosses(PV2<Parameter> a, PV2<Parameter> b, PV2<Parameter> c,
             PV2<Parameter> d) {
  return (b - a).cross(c - a).sign(false) != (b - a).cross(d - a).sign(false) &&
         (d - c).cross(a - c).sign(false) != (d - c).cross(b - c).sign(false);
}

// does the edge ab intersect the view box
bool crosses_view(PV2<Parameter> a, PV2<Parameter> b) {
  PV2<Parameter> ll = transform(PV2<Parameter>::constant(-1, -1));
  PV2<Parameter> lr = transform(PV2<Parameter>::constant(1, -1));
  PV2<Parameter> ur = transform(PV2<Parameter>::constant(1, 1));
  PV2<Parameter> ul = transform(PV2<Parameter>::constant(-1, 1));

  return crosses(a, b, ll, lr) || crosses(a, b, lr, ur) ||
         crosses(a, b, ur, ul) || crosses(a, b, ul, ll);
}

// is this point in the view box
bool inside_view(PV2<Parameter> p) {
  double x = p.x.mid();
  double y = p.y.mid();

  double lx = (-1 / draw_scale) - x_off;
  double rx = (1 / draw_scale) - x_off;
  double ly = (-1 / draw_scale) - y_off;
  double uy = (1 / draw_scale) - y_off;

  return x <= rx && x >= lx && y >= ly && y <= uy;
}

// is the view box contained entirely inside this face
bool inside_face(PTR<Face> f) {
  /*
   PV2<Parameter> _ll = transform(PV2<Parameter>::constant(-1, -1));
   PV2<Parameter> _lr = transform(PV2<Parameter>::constant( 1, -1));
   PV2<Parameter> _ur = transform(PV2<Parameter>::constant( 1,  1));
   PV2<Parameter> _ul = transform(PV2<Parameter>::constant(-1,  1));

   PTR<Object<PV2>> ll = new InputPoint(_ll);
   PTR<Object<PV2>> lr = new InputPoint(_lr);
   PTR<Object<PV2>> ur = new InputPoint(_ur);
   PTR<Object<PV2>> ul = new InputPoint(_ul);
   */

  //  PTR<Object<PV2>> ll = new TransformPoint(new
  //  InputPoint(PV2<Parameter>::input(-1, -1)) , transform ); PTR<Object<PV2>>
  //  lr = new TransformPoint(new InputPoint(PV2<Parameter>::input( 1, -1)) ,
  //  transform ); PTR<Object<PV2>> ur = new TransformPoint(new
  //  InputPoint(PV2<Parameter>::input( 1,  1)) , transform ); PTR<Object<PV2>>
  //  ul = new TransformPoint(new InputPoint(PV2<Parameter>::input(-1,  1)) ,
  //  transform );
  PTR<Object<PV2>> ll =
      new TransformPoint(new InputPoint(PV2<Parameter>::constant(-1, -1)),
                         post_scale, draw_scale, x_off, y_off);
  PTR<Object<PV2>> lr =
      new TransformPoint(new InputPoint(PV2<Parameter>::constant(1, -1)),
                         post_scale, draw_scale, x_off, y_off);
  PTR<Object<PV2>> ur =
      new TransformPoint(new InputPoint(PV2<Parameter>::constant(1, 1)),
                         post_scale, draw_scale, x_off, y_off);
  PTR<Object<PV2>> ul =
      new TransformPoint(new InputPoint(PV2<Parameter>::constant(-1, 1)),
                         post_scale, draw_scale, x_off, y_off);

  return f->inside(ll) && f->inside(lr) && f->inside(ur) && f->inside(ul);
}

// either a vertex of this face is inside of the view or an edge crosses the
// view
bool intersects_view(PTR<Face> f) {
  PTR<Edge> e = f->e;

  // some vertex of this face is in the view box, or an edge crosses
  do {
    if (inside_view(e->tail->p->getApprox(1.0)) ||
        crosses_view(e->tail->p->getApprox(1.0), e->head->p->getApprox(1.0)))
      return true;

    e = e->next;
  } while (e != f->e);

  // or the entire view box is inside this face
  if (inside_face(f)) return true;

  return false;
}

PV2<Parameter> newton(Poly2<Parameter> f, PV2<Parameter> p) {
  Poly2<Parameter> dx = f.derX();
  Poly2<Parameter> dy = f.derY();

  Parameter xy[2] = {Parameter::constant(p.x.mid()),
                     Parameter::constant(p.y.mid())};
  Parameter fp = Parameter::constant(f.value(xy).mid());

  int count = 0;

  do {
    PV2<Parameter> grad =
        PV2<Parameter>::constant(dx.value(xy).mid(), dy.value(xy).mid());

    Parameter t = -fp / (grad.dot(grad));
    p = p + t * grad;

    xy[0] = Parameter::constant(p.x.mid());
    xy[1] = Parameter::constant(p.y.mid());
    fp = Parameter::constant(f.value(xy).mid());

  } while (++count < 100 && fp.mid() > 0.00000000000001);

  return p;
}

void curve_trace_subdivide(PTR<Face> f, int points_size, int sign1 = -5,
                           int sign2 = -5, bool no_face = false) {
  PTR<Edge> e = f->e;

  points_size = max(points_size, 2);

  static unordered_map<Object<Poly2>*, Color> color_map;

  PTR<Object<Scalar>> ZERO = new InputScalar(Parameter::constant(0));

  glColor3f(1.0f, 0.0f, 0.0f);

  set<PTR<Edge>> edge_set;

  do {
    if (e->intersection) break;
  } while ((e = e->next) != f->e);

  if (e->intersection) {
    do {
      edge_set.insert(e);

      PTR<Edge> e2 = e->next;

      while (e2 != f->e && (!e2->intersection || e2->poly != e->poly ||
                            edge_set.find(e2) != edge_set.end()))
        e2 = e2->next;

      if (e2 == f->e) break;

      edge_set.insert(e2);

      vector<PTR<Object<PV2>>> points;
      points.push_back(e->intersection_point());
      points[0] = new InputPoint(points[0]->getApprox(1.0).x.mid(),
                                 points[0]->getApprox(1.0).y.mid());
      points.push_back(e2->intersection_point());
      points[1] = new InputPoint(points[1]->getApprox(1.0).x.mid(),
                                 points[1]->getApprox(1.0).y.mid());

      while (points.size() < points_size) {
        vector<PTR<Object<PV2>>> new_points;

        for (int i = 0; i < points.size() - 1; i++) {
          PTR<Object<PV2>> a = points[i];
          PTR<Object<PV2>> b = points[i + 1];

          PTR<Object<PV2>> ab = new VectorAB(a, b);
          PTR<Object<PV2>> ab_perp = new Rot90(ab);

          PTR<Object<PV2>> mid = new MidVector(a, b);
          mid = new InputPoint(mid->getApprox(1.0).x.mid(),
                               mid->getApprox(1.0).y.mid());

          PTR<Line> line = new Line(mid, ab_perp);

          PTR<Object<Poly>> line_poly = new LinePoly(line, e->poly);

          vector<PTR<Edge>> hits = f->intersection_edges(line);

          assert(hits.size() == 2);

          PTR<Object<Scalar>> h0 = new LLHit(line, hits[0]->line);
          PTR<Object<Scalar>> h1 = new LLHit(line, hits[1]->line);

          if (Poly1Sign(line_poly, h0) == Poly1Sign(line_poly, h1)) {
            vector<PTR<Object<Scalar>>> h0roots;
            vector<PTR<Object<Scalar>>> h1roots;

            if (LessThan(ZERO, h0) == -1)
              h0roots = PolySolver(line_poly).getRoots(ZERO, h0);
            else
              h0roots = PolySolver(line_poly).getRoots(h0, ZERO);

            if (LessThan(ZERO, h1) == -1)
              h1roots = PolySolver(line_poly).getRoots(ZERO, h1);
            else
              h1roots = PolySolver(line_poly).getRoots(h1, ZERO);

            mid = nullptr;
            PTR<Object<Scalar>> h0r = h0roots.size() > 0 ? h0roots[0] : nullptr;
            PTR<Object<Scalar>> h1r = h1roots.size() > 0 ? h1roots[0] : nullptr;

            if (h0roots.size() > 0 || h1roots.size() > 0) {
              PTR<Object<Scalar>> min_root =
                  (h0roots.size() <= 0 || LessThanAbs(h1r, h0r) == -1) ? h1r
                                                                       : h0r;

              mid = new LinePoint(line, min_root);
              mid = new InputPoint(mid->getApprox(1.0).x.mid(),
                                   mid->getApprox(1.0).y.mid());
            }

            if (!mid) {
              printf(
                  "ERROR didn't intersect curve in subdivisision curve "
                  "tracting!\n");
            }

          } else {
            // Experimental
            PTR<Object<Scalar>> r;
            Parameter rint;
            if (LessThan(h0, h1) == -1)
              rint = Parameter::interval(h0->getApprox(1.0).x.lb(),
                                         h1->getApprox(1.0).x.ub());
            else
              rint = Parameter::interval(h1->getApprox(1.0).x.lb(),
                                         h0->getApprox(1.0).x.ub());

            // rint = line_poly->polish(rint);

            r = new PolyRoot(line_poly, new Object<Scalar>(rint.lbP()),
                             new Object<Scalar>(rint.ubP()));
            // LessThan(r, ZERO) == 0;

            mid = new LinePoint(line, r);
            mid = new InputPoint(mid->getApprox(1.0).x.mid(),
                                 mid->getApprox(1.0).y.mid());
          }

          if (i == 0) new_points.push_back(a);

          if (!mid) {
            printf(
                "ERROR didn't intersect curve in subdivisision curve "
                "tracting!\n");
          } else {
            new_points.push_back(mid);
          }

          new_points.push_back(b);
        }

        points = new_points;
      }

      if (color_map.find((Object<Poly2>*)e->poly) == color_map.end()) {
        Color c;
        if (encasement->curves.size() > 0 && e->poly == encasement->curves[0])
          c = {0.0f, 1.0f, 0.0f};
        else if (encasement->curves.size() > 1 &&
                 e->poly == encasement->curves[1])
          c = {1.0f, 0.0f, 0.0f};
        else if (encasement->curves.size() > 2 &&
                 e->poly == encasement->curves[2])
          c = {0.0f, 0.0f, 1.0f};
        else
          c = {(float)randomNumber(0.5, 1), (float)randomNumber(0.5, 1),
               (float)randomNumber(0.5, 1)};

        color_map[(Object<Poly2>*)e->poly] = c;
      }

      Color col = color_map[(Object<Poly2>*)e->poly];

      glColor3f(col.r, col.g, col.b);

      for (int i = 0; i < points.size() - 1; i++) {
        PTR<Object<PV2>> a = points[i];
        PTR<Object<PV2>> b = points[i + 1];

        if (sign1 == 0 && sign2 == 0) {
          glColor3f(0.0f, 0.0f, 0.0f);
          glLineWidth(3.0f);
        } else if (sign1 != -5 && sign2 != -5) {
          bool swap_black = false;

          if (e->poly == f->f) {
            int sign = Sign(new PolyScalar(f->g, b));
            swap_black = (sign == sign1);
          } else {
            int sign = Sign(new PolyScalar(f->f, b));
            swap_black = (sign == sign2);
          }

          if (swap_black) {
            glColor3f(0.0f, 0.0f, 0.0f);
            glLineWidth(3.0f);
          } else {
            glColor3f(col.r, col.g, col.b);
            glLineWidth(linewidth);
          }
        }

        draw_line(a, b);
      }

      while ((e = e->next) != f->e) {
        if (e->intersection && edge_set.find(e) == edge_set.end()) break;
      }

    } while (e != f->e);
  }

  glLineWidth(linewidth);

  if (no_face) return;

  if (draw_triangulation) return;

  e = f->e;

  PTR<Object<PV2>> o = e->tail->p;

  e = e->next;

  // float f3 = (float)(rand() / (1.0 * RAND_MAX));
  // float f3 = (float)(randomNumber(0, 1));
  float f3 = 0.5f;

  glColor4f(0.0f, 0.0f, f3, 0.2f);

  // paint the face
  do {
    draw_triangle(o, e->tail->p, e->head->p);

    e = e->next;

  } while (e->next != f->e);

  if (!draw_encasement) return;

  e = f->e;

  if (!draw_roots) {
    // draw the edges, black
    do {
      glColor3f(0.0f, 0.0f, 0.0f);

      draw_line(e->tail->p, e->head->p);

      // no edge intersections for now
      if (false && e->intersection != nullptr) {
        glColor3f(0.0f, 1.0f, 0.0f);
        PTR<Object<PV2>> ip = e->intersection_point();
        std::vector<PTR<Object<PV2>>> points;
        points.push_back(ip);
        draw_points(points);
        points.clear();
      }

      e = e->next;

    } while (e != f->e);
  }
}

void curve_trace(PTR<Face> f, int first, int sign1 = -5, int sign2 = -5,
                 bool no_face = false) {
  PTR<Edge> e = f->e;

  double step_size = (0.01 / draw_scale);
  // printf("step_size: %g\ndraw_scale: %g\n\n", step_size, draw_scale);
  // double step_size = 0.01;

  static unordered_map<Object<Poly2>*, Color> color_map;

  glColor3f(1.0f, 0.0f, 0.0f);

  // trace the actual curves
  do {
    if (e->intersection != nullptr) {
      Poly2<Parameter> dx = e->poly->getApprox(1.0).derX();
      Poly2<Parameter> dy = e->poly->getApprox(1.0).derY();

      PV2<Parameter> ip = e->line->getP<Parameter>(e->intersection);
      ip = PV2<Parameter>::constant(ip.x.mid(), ip.y.mid());

      static int not_found = 0;
      if (color_map.find((Object<Poly2>*)e->poly) == color_map.end()) {
        Color c;
        if (encasement->curves.size() > 0 && e->poly == encasement->curves[0])
          c = {0.0f, 1.0f, 0.0f};
        else if (encasement->curves.size() > 1 &&
                 e->poly == encasement->curves[1])
          c = {1.0f, 0.0f, 0.0f};
        else if (encasement->curves.size() > 2 &&
                 e->poly == encasement->curves[2])
          c = {0.0f, 0.0f, 1.0f};
        else
          c = {(float)randomNumber(0.5, 1), (float)randomNumber(0.5, 1),
               (float)randomNumber(0.5, 1)};

        color_map[(Object<Poly2>*)e->poly] = c;
        not_found++;
      }

      Color col = color_map[(Object<Poly2>*)e->poly];

      glColor3f(col.r, col.g, col.b);

      int count = 0;

      bool t_select = false;

      while (count < 1000) {
        Parameter xy[] = {ip.x, ip.y};
        PV2<Parameter> v =
            PV2<Parameter>::constant(-dy.value(xy).mid(), dx.value(xy).mid())
                .unit();

        PV2<Parameter> np = ip + v * step_size;

        np = newton(e->poly->getApprox(1.0), np);

        PTR<Object<PV2>> ipp = new InputPoint(PV2<Parameter>::constant(
            ip.x.mid(), ip.y.mid()));  // this is ok because ip is constant
        PTR<Object<PV2>> npp =
            new InputPoint(PV2<Parameter>::constant(np.x.mid(), np.y.mid()));

        if (sign1 == 0 && sign2 == 0) {
          glColor3f(0.0f, 0.0f, 0.0f);
          glLineWidth(3.0f);
        } else if (sign1 != -5 && sign2 != -5) {
          bool swap_black = false;

          if (e->poly == f->f) {
            int sign = Sign(new PolyScalar(f->g, npp));
            swap_black = (sign == sign1);
          } else {
            int sign = Sign(new PolyScalar(f->f, npp));
            swap_black = (sign == sign2);
          }

          if (swap_black) {
            glColor3f(0.0f, 0.0f, 0.0f);
            glLineWidth(3.0f);
          } else {
            glColor3f(col.r, col.g, col.b);
            glLineWidth(linewidth);
          }
        }

        draw_line(ipp, npp);

        if (!f->inside(npp)) break;

        count++;

        ip = np;
      }

      count = 0;

      while (count < 1000) {
        Parameter xy[] = {ip.x, ip.y};
        PV2<Parameter> v =
            -PV2<Parameter>::constant(-dy.value(xy).mid(), dx.value(xy).mid())
                 .unit();

        PV2<Parameter> np = ip + v * step_size;

        np = newton(e->poly->getApprox(1.0), np);

        PTR<Object<PV2>> ipp = new InputPoint(PV2<Parameter>::constant(
            ip.x.mid(), ip.y.mid()));  // this is ok because ip is constant
        PTR<Object<PV2>> npp =
            new InputPoint(PV2<Parameter>::constant(np.x.mid(), np.y.mid()));

        if (sign1 == 0 && sign2 == 0) {
          glColor3f(0.0f, 0.0f, 0.0f);
          glLineWidth(3.0f);
        } else if (sign1 != -5 && sign2 != -5) {
          bool swap_black = false;

          if (e->poly == f->f) {
            int sign = Sign(new PolyScalar(f->g, npp));
            swap_black = (sign == sign1);
          } else {
            int sign = Sign(new PolyScalar(f->f, npp));
            swap_black = (sign == sign2);
          }

          if (swap_black) {
            glColor3f(0.0f, 0.0f, 0.0f);
            glLineWidth(3.0f);
          } else {
            glColor3f(col.r, col.g, col.b);
            glLineWidth(linewidth);
          }
        }

        draw_line(ipp, npp);

        if (!f->inside(npp)) break;

        count++;

        ip = np;
      }
    }

    e = e->next;

  } while (e != f->e);

  glLineWidth(linewidth);

  if (no_face) return;

  if (draw_triangulation) return;

  e = f->e;

  PTR<Object<PV2>> o = e->tail->p;

  e = e->next;

  // float f3 = (float)(rand() / (1.0 * RAND_MAX));
  // float f3 = (float)(randomNumber(0, 1));
  float f3 = 0.5f;

  if (dfirst && first == 0) {
    glColor4f(0.0f, 1.0f, 0.0f, 0.4f);
    printf("bad1 %d: %p\n", zoom_index, (Face*)f);
  } else if (!dfirst && first == 1) {
    glColor4f(1.0f, 0.0f, 0.0f, 0.4f);
    printf("bad2 %d: %p\n", zoom_index, (Face*)f);
  } else {
    // glColor4f(0.0f, 0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 0.0f, f3, 0.2f);
  }

  // paint the face
  do {
    draw_triangle(o, e->tail->p, e->head->p);

    e = e->next;

  } while (e->next != f->e);

  if (!draw_encasement) return;

  e = f->e;

  if (!draw_roots) {
    // draw the edges, black
    do {
      glColor3f(0.0f, 0.0f, 0.0f);

      draw_line(e->tail->p, e->head->p);

      // no edge intersections for now
      if (false && e->intersection != nullptr) {
        glColor3f(0.0f, 1.0f, 0.0f);
        PTR<Object<PV2>> ip = e->intersection_point();
        std::vector<PTR<Object<PV2>>> points;
        points.push_back(ip);
        draw_points(points);
        points.clear();
      }

      e = e->next;

    } while (e != f->e);
  }
}

void curve_trace2(PTR<Face> face) {
  // BAD!!!!

  // Transform point
  // constructor (point, transform)
  // calculate transforms the point

  // PTR<Object<PV2>> ll = new InputPoint(transform(PV2<Parameter>::input(-1,
  // -1))); PTR<Object<PV2>> lr = new
  // InputPoint(transform(PV2<Parameter>::input( 1, -1))); PTR<Object<PV2>> ur =
  // new InputPoint(transform(PV2<Parameter>::input( 1,  1))); PTR<Object<PV2>>
  // ul = new InputPoint(transform(PV2<Parameter>::input(-1,  1)));

  //  PTR<Object<PV2>> ll = new TransformPoint(new
  //  InputPoint(PV2<Parameter>::input(-1, -1)) , transform ); PTR<Object<PV2>>
  //  lr = new TransformPoint(new InputPoint(PV2<Parameter>::input( 1, -1)) ,
  //  transform ); PTR<Object<PV2>> ur = new TransformPoint(new
  //  InputPoint(PV2<Parameter>::input( 1,  1)) , transform ); PTR<Object<PV2>>
  //  ul = new TransformPoint(new InputPoint(PV2<Parameter>::input(-1,  1)) ,
  //  transform );

  PTR<Object<PV2>> ll =
      new TransformPoint(new InputPoint(PV2<Parameter>::constant(-1, -1)),
                         post_scale, draw_scale, x_off, y_off);
  PTR<Object<PV2>> lr =
      new TransformPoint(new InputPoint(PV2<Parameter>::constant(1, -1)),
                         post_scale, draw_scale, x_off, y_off);
  PTR<Object<PV2>> ur =
      new TransformPoint(new InputPoint(PV2<Parameter>::constant(1, 1)),
                         post_scale, draw_scale, x_off, y_off);
  PTR<Object<PV2>> ul =
      new TransformPoint(new InputPoint(PV2<Parameter>::constant(-1, 1)),
                         post_scale, draw_scale, x_off, y_off);

  PTR<Object<PV2>> bp = new MidVector(ll, lr);
  PTR<Object<PV2>> rp = new MidVector(lr, ur);
  PTR<Object<PV2>> tp = new MidVector(ul, ul);
  PTR<Object<PV2>> lp = new MidVector(ul, ll);

  PTR<Object<PV2>> bv = new VectorAB(ll, lr);
  PTR<Object<PV2>> rv = new VectorAB(lr, ur);
  PTR<Object<PV2>> tv = new VectorAB(ur, ul);
  PTR<Object<PV2>> lv = new VectorAB(ul, ll);

  PTR<Line> l1 = new Line(bp, bv);
  PTR<Line> l2 = new Line(rp, rv);
  PTR<Line> l3 = new Line(tp, tv);
  PTR<Line> l4 = new Line(lp, lv);

  // clone the face and intersect the clone with the zoom box
  PTR<Face> f_clone = face->clone();

  // each successful split creates 1 new face and 6 new edges
  // need to remove them from the global vectors right after
  // creating so they aren't involed in the encasement computation
  //
  // split always returns the face on the right, which is not the one
  // we want so we have to delete it right after
  PTR<Face> other;
  other = f_clone->split(l1);
  other = f_clone->split(l2);
  other = f_clone->split(l3);
  other = f_clone->split(l4);

  int first = (face == encasement->bad_faces.first
                   ? 0
                   : (face == encasement->bad_faces.second ? 1 : -1));

  // curve trace the cut up face
  curve_trace(f_clone, first);

  // curve_trace_subdivide(f_clone, subdiv);
}

void arrangement_faces() {
  // find the face the point is in
  PTR<Face> f = encasement->faces[0];
  int i = -1;
  while (++i < encasement->faces.size() && !encasement->faces[i]->inside(click))
    ;

  if (i >= encasement->faces.size()) return;

  PTR<Face> start = encasement->faces[i];

  queue<PTR<Face>> q;
  q.push(start);

  // find first neighboring face that contains a curve
  while (q.size() > 0) {
    start = q.front();
    q.pop();
    if (start->curves.size() > 0) break;

    PTR<Edge> e = start->e;
    do {
      q.push(e->twin->face);
    } while ((e = e->next) != start->e);
  }

  // keep track of the curve we're tracing and sign of f at p

  PTR<Face> fac = start;
  PTR<Edge> start_edge = fac->e;

  while (!start_edge->intersection) start_edge = start_edge->next;

  fac->e = start_edge;

  while (q.size() > 0) q.pop();

  set<PTR<Face>> f_set;

  q.push(fac);
  f_set.insert(fac);

  while (q.size() > 0) {
    PTR<Face> fac = q.front();
    q.pop();

    // if single curve just trace the whole thing
    if (fac->curves.size() == 0) {
      PTR<Edge> e = fac->e;

      do {
        if (e->twin->face && f_set.find(e->twin->face) == f_set.end()) {
          e->twin->face->e = e->twin;
          q.push(e->twin->face);
          f_set.insert(fac);
        }

      } while ((e = e->next) != fac->e);

    } else if (fac->curves.size() == 1) {
      curve_trace(fac, -1, 0, 0, true);
      // curve_trace_subdivide(fac, subdiv, 0, 0, true);

      PTR<Edge> e = fac->e;

      int sign_f = Sign(new PolyScalar(*(fac->curves.begin()), click));

      do {
        if (e->intersection) {
          if (e->twin->face && f_set.find(e->twin->face) == f_set.end()) {
            e->twin->face->e = e->twin;
            q.push(e->twin->face);
            f_set.insert(fac);
          }

        } else {
          if (Sign(new PolyScalar(*(fac->curves.begin()), e->tail->p)) ==
              sign_f) {
            if (e->twin->face && f_set.find(e->twin->face) == f_set.end()) {
              e->twin->face->e = e->twin;
              q.push(e->twin->face);
              f_set.insert(fac);
            }
          }
        }

      } while ((e = e->next) != fac->e);

    } else {
      int sign_f = Sign(new PolyScalar(fac->f, click));
      int sign_g = Sign(new PolyScalar(fac->g, click));

      // if multiple, find intersection of f on boundary with correct sign of g
      //             find intersection of g on boundary with correct sign of f
      //             trace from them to the intersection and then swap
      //             arrangement edge curve
      curve_trace(fac, -1, sign_g, sign_f, true);
      // curve_trace_subdivide(fac, subdiv, sign_g, sign_f, true);

      PTR<Edge> e = fac->e;

      do {
        bool f_or_g = (e->poly == fac->f);

        if (e->intersection && (Sign(new PolyScalar(f_or_g ? fac->g : fac->f,
                                                    e->intersection_point())) ==
                                (f_or_g ? sign_g : sign_f))) {
          if (e->twin->face && f_set.find(e->twin->face) == f_set.end()) {
            e->twin->face->e = e->twin;
            q.push(e->twin->face);
            f_set.insert(fac);
          }

        } else {
          if (Sign(new PolyScalar(fac->f, e->tail->p)) == sign_f &&
              Sign(new PolyScalar(fac->g, e->tail->p)) == sign_g) {
            if (e->twin->face && f_set.find(e->twin->face) == f_set.end()) {
              e->twin->face->e = e->twin;
              q.push(e->twin->face);
              f_set.insert(fac);
            }
          }
        }

      } while ((e = e->next) != fac->e);

      //       if(e->twin->face && f_set.find(e->twin->face) == f_set.end()) {
      //         e->twin->face->e = e->twin;
      //         q.push(e->twin->face);
      //       }
      //
      //       e = fac->e->next;
      //       while(!e->intersection || (e->poly == (f_or_g ? fac->f : fac->g))
      //       || (Sign(new PolyScalar(f_or_g ? fac->f : fac->g,
      //       e->intersection_point())) != sign)) {
      //         e = e->next;
      //       }
      //
      //       if(e->twin->face && f_set.find(e->twin->face) == f_set.end()) {
      //         e->twin->face->e = e->twin;
      //         q.push(e->twin->face);
      //       }
    }
  }
}

void zoom_render() {
  // box goes from (-1, -1) to (1, 1)
  // with scale and offset it's (-1 + x_off, -1 + y_off) * scale to (1 + x_off,
  // 1 + y_off) * scale
  fsize = encasement->faces.size();
  for (int i = 0; i < (fsize > -1 ? fsize : encasement->faces.size()); i++) {
    zoom_index = i;
    PTR<Face> f = encasement->faces[i];

    if (intersects_view(f)) {
      if (cell_index == -1 || i == cell_index) curve_trace2(f);
    }
  }

  if (draw_triangulation) {
    if (roots.size() == 0) {
      fsize = encasement->faces.size();
      roots = encasement->getRoots();
    }

    if (roots.size() == 0) {
      printf("Some error with finding roots\n");
      if (encasement->bad_faces.first)
        curve_trace2(encasement->bad_faces.first);
      if (encasement->bad_faces.second)
        curve_trace2(encasement->bad_faces.second);
    }

    if (triangulation.size() == 0) {
      triangulation = encasement->getTriangulation(roots);
    }

    // float f3 = (float)(rand() / (1.0 * RAND_MAX));
    // float f3 = (float)(randomNumber(0, 1));
    float f3 = 0.5f;
    glColor3f(0.0f, 0.0f, f3);

    // PTR<Object<PV2>> e_bl = new TransformPoint(encasement->bl->p, transform);
    // PTR<Object<PV2>> e_br = new TransformPoint(encasement->br->p, transform);
    // PTR<Object<PV2>> e_tr = new TransformPoint(encasement->tr->p, transform);
    // PTR<Object<PV2>> e_tl = new TransformPoint(encasement->tl->p, transform);

    // PTR<Object<PV2>> e_bl = encasement->bl->p;
    // PTR<Object<PV2>> e_br = encasement->br->p;
    // PTR<Object<PV2>> e_tr = encasement->tr->p;
    // PTR<Object<PV2>> e_tl = encasement->tl->p;

    // draw_triangle(e_bl, e_br, e_tr);
    // draw_triangle(e_bl, e_tr, e_tl);

    glColor3f(0.0f, 0.0f, 0.0f);

    queue<PTR<Triangle>> q;

    for (int i = 0; i < triangulation.size(); i++) {
      q.push(triangulation[i]);
    }

    while (q.size() > 0) {
      PTR<Triangle> t = q.front();
      q.pop();

      // if leaf, draw it. Else add children
      if (t->chld[0] == 0) {
        draw_line(t->vert[0], t->vert[1]);
        draw_line(t->vert[1], t->vert[2]);
        draw_line(t->vert[2], t->vert[0]);

      } else {
        for (int i = 0; i < 3; i++)
          if (t->chld[i]) q.push(t->chld[i]);
      }
    }

  } else if (draw_roots) {
    if (roots.size() == 0) {
      fsize = encasement->faces.size();
      roots = encasement->getRoots();
    }

    if (roots.size() == 0) {
      printf("Some error with finding roots\n");
      if (encasement->bad_faces.first)
        curve_trace2(encasement->bad_faces.first);
      if (encasement->bad_faces.second)
        curve_trace2(encasement->bad_faces.second);
    }

    glColor3f(0.0f, 1.0f, 0.0f);

    // printf("%ld roots\n", roots.size());

    for (int i = 0; i < roots.size(); i++) {
      PV2<Parameter> r = roots[i]->getCurrentP();
      try {
        // r = roots[i]->calculate();
      } catch (PrecisionException ex) {
        printf("not enough precision at double for this root");
      } catch (SignException ex) {
        printf("box too large\n");
      }

      PTR<Object<PV2>> ll =
          new InputPoint(PV2<Parameter>(r.x.lbP(), r.y.lbP()));
      PTR<Object<PV2>> lr =
          new InputPoint(PV2<Parameter>(r.x.ubP(), r.y.lbP()));
      PTR<Object<PV2>> ur =
          new InputPoint(PV2<Parameter>(r.x.ubP(), r.y.ubP()));
      PTR<Object<PV2>> ul =
          new InputPoint(PV2<Parameter>(r.x.lbP(), r.y.ubP()));

      draw_line(ll, lr);
      draw_line(lr, ur);
      draw_line(ur, ul);
      draw_line(ul, ll);
    }
  }
}

void display() {
  if (!calculated) {
    // calculate();
    // calculated = true;
    // TEST TEST TEST
    // exit(0);
  }

  glutSetWindowTitle("encasement");

  float white[4] = {1.0f, 1.0f, 1.0f, 1.0f};
  clear_screen(white);

  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

  // draw the ellipse

  // draw the encasement

  if (calculated) {
    if (draw_crit) drawCriticalPoints();

#ifdef ENCASEMENT_GPU
    if (draw_gpu_rects) drawGPURects();
#endif

    zoom_render();

    if (click) arrangement_faces();

    if (draw_intersections) drawIntersections();

    // vector< PTR<Face> > faces = encasement->faces;
    // PTR<Face> parent = nullptr;
    // for(int i = 0; i < faces.size(); i++) {
    //  if(faces[i]->curves.size() == 1 && (*faces[i]->curves.begin() ==
    //  curves[curves.size()-1])) {
    //    parent = faces[i]->find();
    //    break;
    //  }
    //}

    // if(parent) {

    //  for(int i = 0; i < faces.size(); i++) {
    //    if(faces[i]->find() == parent) {
    //      curve_trace2(faces[i]);
    //    }
    //  }

    //}

    /*
     for(int i = 0; i < roots.size(); i++) {
     PV2<Parameter> r = roots[i]->getApprox(1.0);
     printf("x: %g\n", r.x.intervalWidth());
     printf("y: %g\n\n", r.y.intervalWidth());
     }
     */

    // drawIntersections();
  } else if (encasement != NULL && cell_index == -1 && g_index == -1) {
    drawEllipse();
    // if(state == NONE)
    drawEncasement();

    if (draw_crit) drawCriticalPoints();

#ifdef ENCASEMENT_GPU
    if (draw_gpu_rects) drawGPURects();
#endif

  } else if (cell_index != -1) {
    drawCell(cell_index, true);
  } else {
    drawCell(g_index, false);
    drawEdge();
  }

  // drawRandomPoints();

  // saveImage();

  glutSwapBuffers();

  // glFlush();
}
