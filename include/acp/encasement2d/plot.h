#ifndef PLOT
#define PLOT

#include <GL/glut.h>

#include "acp/core/acp.h"
#include "acp/dcel/dcel.h"
#include "acp/encasement2d/predicates2d.h"

#include <fstream>

extern double draw_scale;
extern double x_off;
extern double y_off;
extern double post_scale;

extern int counter;

extern int linewidth;
extern float pointsize;
extern int wsize;
extern int hsize;

using namespace std;

typedef enum { POINT, LINE, VECTOR, TRI } DrawType;

typedef struct {
  float r, g, b;
} Color;

typedef struct {
  vector<Color> colors;
  vector<vector<PTR<Object<PV2>>>> pvs;
  vector<DrawType> types;
} DrawItem;

void graphics_init(int argc, char** argv, int wsize, int hsize, int linewidth,
                   float pointsize);

void graphics_init_aux(int argc, char** argv);

void menu_init();

void menu(int value);

void mouse(int button, int buttonState, int x, int y);

void mouseMoved(int x, int y);

void keyPressed(unsigned char key, int x, int y);

void reshape(int w, int h);

PTR<Object<PV2>> point(int x, int y);

void display();

void clear_screen(float* color);

void draw_points(const vector<PTR<Object<PV2>>>& pts);

void draw_lines(const vector<PTR<Object<PV2>>>& pts);

void draw_loop(const vector<PTR<Object<PV2>>>& pts);

void draw_line(PTR<Object<PV2>> t, PTR<Object<PV2>> h);

void draw_vector(PTR<Object<PV2>> t, PTR<Object<PV2>> h);

void draw_triangle(PTR<Object<PV2>> a, PTR<Object<PV2>> b, PTR<Object<PV2>> c);

void draw_box(PTR<Object<PV2>> bl, PTR<Object<PV2>> tr);

void glVertex(PTR<Object<PV2>> p);

double debugDrawFacePoint(PTR<Face> face, PTR<Object<Poly2>> f,
                          PTR<Object<Poly2>> g, double width, double dx,
                          double dy, DrawItem item, bool draw);

void debugDrawFaces(vector<PTR<Face>> face, PTR<Object<Poly2>> f,
                    PTR<Object<Poly2>> g, double zoom, double dx, double dy,
                    DrawItem item, bool draw = false, bool hold = false);

void debugDrawFace(PTR<Face> face, PTR<Object<Poly2>> f, PTR<Object<Poly2>> g,
                   double zoom, double dx, double dy, DrawItem item,
                   bool draw = false, bool hold = false);

void debugDraw(vector<PV2<Parameter>> p, vector<PTR<Object<Poly2>>> f);

void debugDrawCrit(DrawItem item);

void saveImage();
#endif
