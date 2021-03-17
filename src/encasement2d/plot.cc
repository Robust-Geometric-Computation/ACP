#include "acp/encasement2d/plot.h"
#include "acp/encasement2d/pngwriter.h"

int window1;

double draw_scale = 0.4;
double x_off = 0.0;
double y_off = 0.0;

int counter = 0;

double post_scale = 1.0;

int linewidth = 2;
float pointsize = 10.0f;
int wsize = 800;
int hsize = 800;

void saveImage() {
  int w = glutGet(GLUT_WINDOW_WIDTH);
  int h = glutGet(GLUT_WINDOW_HEIGHT);

  GLubyte* data = (GLubyte*)malloc(3 * w * h);

  glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, data);

  static long file_name = -1;
  file_name++;

  char c[50];
  memset(c, '\0', 50);
  sprintf(c, "%ld.png", file_name);

  pngwriter png(w, h, 0, c);

  for (int i = 0; i < w; i++) {
    for (int j = 0; j < h; j++) {
      int index = 3 * (i * w + j);
      double r = (unsigned char)data[index];
      double g = (unsigned char)data[index + 1];
      double b = (unsigned char)data[index + 2];
      r /= 256;
      g /= 256;
      b /= 256;
      png.plot(j, i, r, g, b);
    }
  }

  png.close();
}

void graphics_init(int argc, char** argv, int wsize, int hsize, int linewidth,
                   float pointsize) {
  glutInit(&argc, argv);
  glutInitWindowSize(wsize, hsize);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  window1 = glutCreateWindow("");
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutMotionFunc(mouseMoved);
  glutPassiveMotionFunc(mouseMoved);
  glutKeyboardFunc(keyPressed);
  glutDisplayFunc(display);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-1.0, 1.0, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glLineWidth(linewidth);
  glPointSize(pointsize);

  glEnable(GL_BLEND);

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  menu_init();
  enable();

  graphics_init_aux(argc, argv);

  glutMainLoop();
}

void reshape(int w, int h) { glViewport(0, 0, w, h); }

PTR<Object<PV2>> point(int x, int y) {
  double w = glutGet(GLUT_WINDOW_WIDTH);
  double h = glutGet(GLUT_WINDOW_HEIGHT);

  double px = 2 * ((x - 0.5 * w) / w);
  double py = 2 * ((0.5 * h - y) / h);

  px = (px / draw_scale) - x_off;
  py = (py / draw_scale) - y_off;

  return new InputPoint(px, py);
}

void clear_screen(float* color) {
  glClearColor(color[0], color[1], color[2], color[3]);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void draw_points(const vector<PTR<Object<PV2>>>& points) {
  glBegin(GL_POINTS);
  for (int i = 0; i < points.size(); ++i) glVertex(points[i]);
  glEnd();
}

void draw_lines(const vector<PTR<Object<PV2>>>& pts) {
  glBegin(GL_LINE_STRIP);
  for (int i = 0; i < pts.size(); ++i) glVertex(pts[i]);
  glEnd();
}

void draw_loop(const vector<PTR<Object<PV2>>>& pts) {
  glBegin(GL_LINE_LOOP);
  for (int i = 0; i < pts.size(); ++i) glVertex(pts[i]);
  glEnd();
}

void draw_line(PTR<Object<PV2>> t, PTR<Object<PV2>> h) {
  glBegin(GL_LINES);
  glVertex(t);
  glVertex(h);
  glEnd();
}

void draw_vector(PTR<Object<PV2>> t, PTR<Object<PV2>> h) {
  glBegin(GL_LINES);
  glVertex(t);
  glVertex(h);
  glEnd();

  PV2<Parameter> tt = t->getApprox(1.0);
  PV2<Parameter> hh = h->getApprox(1.0);

  PV2<Parameter> ht = tt - hh;
  PV2<Parameter> ht90(-ht.y, ht.x);

  double htx = ht.x.mid();
  double hty = ht.y.mid();
  double htl = sqrt(htx * htx + hty * hty);

  Parameter s = Parameter::constant(htl / 6.0);
  Parameter u = Parameter::constant(htl / 12.0);

  PV2<Parameter> l = hh + s * ht + u * ht90;
  PV2<Parameter> r = hh + s * ht - u * ht90;

  PTR<Object<PV2>> ll = new InputPoint(l.x.mid(), l.y.mid());
  PTR<Object<PV2>> rr = new InputPoint(r.x.mid(), r.y.mid());

  glBegin(GL_LINES);
  glVertex(h);
  glVertex(ll);
  glEnd();

  glBegin(GL_LINES);
  glVertex(h);
  glVertex(rr);
  glEnd();
}

void draw_triangle(PTR<Object<PV2>> a, PTR<Object<PV2>> b, PTR<Object<PV2>> c) {
  glBegin(GL_TRIANGLES);
  glVertex(a);
  glVertex(b);
  glVertex(c);
  glEnd();
}

void draw_box(PTR<Object<PV2>> bl, PTR<Object<PV2>> tr) {
  glBegin(GL_LINE_LOOP);
  glVertex(bl);
  glVertex2d(draw_scale * (tr->getApprox(1.0).x.mid() + x_off),
             draw_scale * (bl->getApprox(1.0).y.mid() + y_off));
  glVertex(tr);
  glVertex2d(draw_scale * (bl->getApprox(1.0).x.mid() + x_off),
             draw_scale * (tr->getApprox(1.0).y.mid() + y_off));
  glEnd();
}

void glVertex(PTR<Object<PV2>> p) {
  glVertex2d(draw_scale * (p->getApprox(1.0).x.mid() + x_off),
             draw_scale * (p->getApprox(1.0).y.mid() + y_off));
}

int main(int argc, char** argv) {
  srand(time(0));
  srandom(time(0));
  // srandom(time(0));
  // srand(time(0));

  // printf("enabling parameter\n");
  enable();

  // srand(1);

  graphics_init(argc, argv, wsize, hsize, linewidth, pointsize);

  return 0;
}
