#include <ctime>
#include <iomanip>
#include <sstream>

#include <gflags/gflags.h>
#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkDirectory.h>
#include <vtkDiskSource.h>
#include <vtkDoubleArray.h>
#include <vtkFeatureEdges.h>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkJPEGReader.h>
#include <vtkNamedColors.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataReader.h>
#include <vtkPropPicker.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkWindowToImageFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtksys/SystemTools.hxx>

#include "acp/encasement3d/encasement3d.h"
#include "acp/encasement3d/encasement_utils.h"
#include "acp/encasement3d/mesh3d.h"
#include "acp/linmath/pv.h"
#include "acp/poly/poly3d.h"

// Flags
DEFINE_double(extent, 0.5, "Extent of the bounding box");

DEFINE_double(mesh_epsilon, 1e-3, "Epsilon value controlling mesh fineness");

DEFINE_string(surface, "torus",
              "Specify the surfaces to mesh in a comma delimited list. "
              "[--surface={torus,cayley,clebsch,steiner,hunt,spheres}]");

DEFINE_int32(num_spheres, 1, "number of spheres");

DEFINE_string(background_image, "paraview",
              "Specify the background image to use. "
              "[--background_image={paraview,gradient}]");
DEFINE_double(angle, M_PI / 3.0,
              "angle (radians) to randomly rotate the surface");

// Globals
using namespace acp;

// Maps between ACP face and VTK actor
std::map<Encasement3D::Face*, vtkActor*> face_to_actor;
std::map<vtkActor*, Encasement3D::Face*> actor_to_face;

// Maps between ACP face and VTK mesh
std::map<Encasement3D::Face*, vtkPolyData*> face_to_mesh;
std::map<Encasement3D::Face*, vector<Mesh3D::Triangle>> face_to_acp_mesh;

double colors[6][3] = {
    {0.0, 0.8, 0.8},  // teal
    {0.8, 0.0, 0.8},  // purple
    {0.0, 0.8, 0.0},  // green
    {0.0, 0.0, 0.8},  // blue
    {0.8, 0.8, 0.0},  // yellow
    {0.3, 0.8, 0.5}   // ???
};

std::map<vtkActor*, int> actor_to_color;

// Handle mouse events
class HighlightFaceStyle : public vtkInteractorStyleTrackballCamera {
 public:
  static HighlightFaceStyle* New();
  vtkTypeMacro(HighlightFaceStyle, vtkInteractorStyleTrackballCamera);

  virtual void OnLeftButtonDown() override {
    int* clickPos = this->GetInteractor()->GetEventPosition();

    vtkSmartPointer<vtkPropPicker> picker =
        vtkSmartPointer<vtkPropPicker>::New();
    picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

    if (picker->GetActor()) {
      double rgb[3];
      picker->GetActor()->GetProperty()->GetColor(rgb);
      if (rgb[0] != selected[0] && rgb[1] != selected[1] &&
          rgb[2] != selected[2]) {
        picker->GetActor()->GetProperty()->SetColor(selected[0], selected[1],
                                                    selected[2]);
      } else {
        picker->GetActor()->GetProperty()->SetColor(
            colors[actor_to_color[picker->GetActor()]]);
      }

      Encasement3D::Face* face = actor_to_face[picker->GetActor()];
      std::cout << "face: " << face << "\n";
    }
    // Forward events
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
  }

  virtual void OnRightButtonDown() override {
    int* clickPos = this->GetInteractor()->GetEventPosition();

    vtkSmartPointer<vtkPropPicker> picker =
        vtkSmartPointer<vtkPropPicker>::New();
    picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

    if (picker->GetActor() &&
        actor_to_face.find(picker->GetActor()) != actor_to_face.end()) {
      double opacity = picker->GetActor()->GetProperty()->GetOpacity();
      picker->GetActor()->SetVisibility(0);
      if (opacity == 0.1) {
        picker->GetActor()->GetProperty()->SetOpacity(1.0);
      } else {
        picker->GetActor()->GetProperty()->SetOpacity(0.1);
      }
    } else {
      for (auto& [face, actor] : face_to_actor) {
        actor->GetProperty()->SetOpacity(1.0);
        actor->SetVisibility(1);
      }
    }
    // Forward events
    vtkInteractorStyleTrackballCamera::OnRightButtonDown();
  }

  virtual void OnKeyPress() {
    // Get the keypress.
    vtkRenderWindowInteractor* rwi = this->Interactor;
    std::string key = rwi->GetKeySym();

    // Handle an space key.
    if (key == "space") {
      std::cout << "The space key was pressed." << std::endl;

      vtkRenderer* renderer = this->GetDefaultRenderer();

      if (!renderer) {
        std::cout << "Renderer null!!!\n";
      } else {
        vtkWindow* renderWindow = renderer->GetRenderWindow();

        if (!renderWindow) {
          std::cout << "Render window null!!!\n";
        } else {
          // Screenshot
          vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
              vtkSmartPointer<vtkWindowToImageFilter>::New();
          windowToImageFilter->SetInput(renderWindow);
          windowToImageFilter->SetInputBufferTypeToRGBA();
          windowToImageFilter->ReadFrontBufferOff();
          windowToImageFilter->Update();

          vtkSmartPointer<vtkPNGWriter> writer =
              vtkSmartPointer<vtkPNGWriter>::New();
          // Timestamped screenshot filename.
          auto t = std::time(nullptr);
          auto tm = *std::localtime(&t);

          std::ostringstream oss;
          oss << "screenshot_";
          oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
          oss << ".png";
          auto filename = oss.str();
          writer->SetFileName(filename.c_str());
          writer->SetInputConnection(windowToImageFilter->GetOutputPort());
          writer->Write();
        }
      }
    }

    // Forward events
    vtkInteractorStyleTrackballCamera::OnKeyPress();
  }

 private:
  const double selected[3] = {0.9, 0.0, 0.0};
};

vtkStandardNewMacro(HighlightFaceStyle);

class UnitNormal : public Object<PV3> {
  PTR<Object<Poly3D>> p;
  PTR<Object<PV3>> v;
  DeclareCalculate(PV3) { return p->get<N>().gradient(v->get<N>()).unit(); }

 public:
  UnitNormal(PTR<Object<Poly3D>> p, PTR<Object<PV3>> v) : p(p), v(v) {}
};

class PV3withZeros : public Object<PV3> {
  class V3I : public Primitive {
    PTR<Object<PV3>> p;
    int i;
    DeclareSign { return p->get<N>()[i]; }

   public:
    V3I(PTR<Object<PV3>> p, int i) : p(p), i(i) {}
  };

  PTR<Object<PV3>> pv3;
  bool zero[3];

  DeclareCalculate(PV3) {
    PV3<N> p = pv3->get<N>();
    for (int i = 0; i < 3; ++i)
      if (zero[i]) p[i] = N(0);
    return p;
  }

 public:
  PV3withZeros(PTR<Object<PV3>> p) : pv3(p) {
    for (int i = 0; i < 3; ++i) zero[i] = (V3I(p, i) == 0);
  }
};

PTR<Object<Poly3D>> makeHunt() {
  Poly3D<Parameter> hunt(6, 6, 6);

  //   4·z6
  hunt.set(0, 0, 6, Parameter::input(4));
  // + 12·y2·z4
  hunt.set(0, 2, 4, Parameter::input(12));
  // + 12·x2·z4
  hunt.set(2, 0, 4, Parameter::input(12));
  // + 276·z4
  hunt.set(0, 0, 4, Parameter::input(276));
  // + 12·y4·z2
  hunt.set(0, 4, 2, Parameter::input(12));
  // + 24·x2·y2·z2
  hunt.set(2, 2, 2, Parameter::input(24));
  // + −528·y2·z2
  hunt.set(0, 2, 2, Parameter::input(-528));
  // + 12·x4·z2
  hunt.set(4, 0, 2, Parameter::input(12));
  // + −960·x2·z2
  hunt.set(2, 0, 2, Parameter::input(-960));
  // + 4620·z2
  hunt.set(0, 0, 2, Parameter::input(4620));
  // + 4·y6
  hunt.set(0, 6, 0, Parameter::input(4));
  // + 12·x2·y4
  hunt.set(2, 4, 0, Parameter::input(12));
  // + −129·y4
  hunt.set(0, 4, 0, Parameter::input(-129));
  // + 12·x4·y2
  hunt.set(4, 2, 0, Parameter::input(12));
  // + −150·x2·y2
  hunt.set(2, 2, 0, Parameter::input(-150));
  // + 1380·y2
  hunt.set(0, 2, 0, Parameter::input(1380));
  // +  4·x6
  hunt.set(6, 0, 0, Parameter::input(4));
  // + 87·x4
  hunt.set(4, 0, 0, Parameter::input(87));
  // + 84·x2
  hunt.set(2, 0, 0, Parameter::input(84));
  // + −4900
  hunt.set(0, 0, 0, Parameter::input(-4900));

  Poly3D<Parameter> eps(0, 0, 0);
  eps.set(0, 0, 0, Parameter::constant(1e-4));

  hunt = hunt * eps;

  for (int ix = 0; ix <= 6; ++ix)
    for (int iy = 0; iy <= 6; ++iy)
      for (int iz = 0; iz <= 6; ++iz)
        if (ix + iy + iz <= 6)
          hunt.add(
              ix, iy, iz,
              Parameter::input(0.00001 * ix + 0.00002 * iy + 0.00004 * iz));

  for (Parameter& p : hunt.a) p = Parameter::input(p.mid());

  return new Object<Poly3D>(hunt);
}

// x²y² + x²z² + y²z² + xyz
PTR<Object<Poly3D>> makeSteinerRoman() {
  Poly3D<Parameter> steiner(3, 3, 3);

  steiner.set(2, 2, 0, Parameter::input(1));
  steiner.set(2, 0, 2, Parameter::input(1));
  steiner.set(0, 2, 2, Parameter::input(1));
  steiner.set(1, 1, 1, Parameter::input(1));

  for (int ix = 0; ix <= 3; ++ix)
    for (int iy = 0; iy <= 3; ++iy)
      for (int iz = 0; iz <= 3; ++iz)
        if (ix + iy + iz <= 3)
          steiner.add(
              ix, iy, iz,
              Parameter::input(0.00001 * ix + 0.00002 * iy + 0.00004 * iz));

  for (Parameter& p : steiner.a) p = Parameter::input(p.mid());

  return new Object<Poly3D>(steiner);
}

// 3x²y + 3x²z + 3x² + 3xy² + 6xyz + 6xy + 3xz² + 6xz +
// 3x + 3y²z + 3y² + 3yz² + 6yz + 3y + 3z² + 3z
PTR<Object<Poly3D>> makeClebschDiagonalCubic() {
  Poly3D<Parameter> clebsch(3, 3, 3);

  // original:
  //   3x²y + 3x²z + 3x² + 3xy² + 6xyz + 6xy + 3xz² + 6xz +
  //   3x + 3y²z + 3y² + 3yz² + 6yz + 3y + 3z² + 3z
  //
  // rotated around x by pi/3:
  // c = cos(pi/3)    s = sin(pi/3)
  // x' = x
  // y' = cy - sz
  // z' = sy + cz
  //
  const double theta = FLAGS_angle;
  const double c = cos(theta);
  const double s = sin(theta);
  const double c2 = c * c;
  const double s2 = s * s;
  const double c3 = c2 * c;
  const double s3 = s2 * s;
  const double cs = c * s;
  const double c2s = c2 * s;
  const double cs2 = c * s2;

  // Linear terms
  // + 3x
  clebsch.set(1, 0, 0, Parameter::input(3));
  // + (3c + 3s)y
  clebsch.set(0, 1, 0, Parameter::input(3 * c + 3 * s));
  // + (3c - 3s)z
  clebsch.set(0, 0, 1, Parameter::input(3 * c - 3 * s));

  // Quadratic terms
  // + 3x²
  clebsch.set(2, 0, 0, Parameter::input(3));
  // + (3c² + 6cs + 3s²)y²
  clebsch.set(0, 2, 0, Parameter::input(3 * c2 + 6 * cs + 3 * s2));
  // + (3c² - 6cs + 3s²)z²
  clebsch.set(0, 0, 2, Parameter::input(3 * c2 - 6 * cs + 3 * s2));
  // + (6c + 6s)xy
  clebsch.set(1, 1, 0, Parameter::input(6 * c + 6 * s));
  // + (6c - 6s)xz
  clebsch.set(1, 0, 1, Parameter::input(6 * c - 6 * s));
  // + (6c² - 6s²)yz
  clebsch.set(0, 1, 1, Parameter::input(6 * c2 - 6 * s2));

  // Cubic terms
  // + (6c² - 6s²)xyz
  clebsch.set(1, 1, 1, Parameter::input(6 * c2 - 6 * s2));
  // + (3c + 3s)x²y
  clebsch.set(2, 1, 0, Parameter::input(3 * c + 3 * s));
  // + (3c - 3s)x²z
  clebsch.set(2, 0, 1, Parameter::input(3 * c - 3 * s));
  // + (3c² + 6cs + 3s²)xy²
  clebsch.set(1, 2, 0, Parameter::input(3 * c2 + 6 * cs + 3 * s2));
  // + (3c³ + 6c²s - 6cs² - 3s³)y²z
  clebsch.set(0, 2, 1,
              Parameter::input(3 * c3 + 6 * c2 * s - 6 * c * s2 - 3 * s3));
  // + (3c² - 6cs + 3s²)xz²
  clebsch.set(1, 0, 2, Parameter::input(3 * c2 - 6 * cs + 3 * s2));
  // + (3c³ - 6c²s - 6cs² + 3s³)yz²
  clebsch.set(0, 1, 2,
              Parameter::input(3 * c3 - 6 * c2 * s - 6 * c * s2 + 3 * s3));

  // + (3c²s + 3cs²)y³
  clebsch.set(0, 3, 0, Parameter::input(3 * c2 * s + 3 * c * s2));
  // + (-3c²s + 3cs²)z³
  clebsch.set(0, 0, 3, Parameter::input(-3 * c2 * s + 3 * c * s2));

  Poly3D<Parameter> eps(0, 0, 0);
  eps.set(0, 0, 0, Parameter::constant(1e-2));

  clebsch = clebsch * eps;

  for (int ix = 0; ix <= 3; ++ix)
    for (int iy = 0; iy <= 3; ++iy)
      for (int iz = 0; iz <= 3; ++iz)
        if (ix + iy + iz <= 3)
          clebsch.add(
              ix, iy, iz,
              Parameter::input(0.0001 * ix + 0.0002 * iy + 0.0004 * iz));

  for (Parameter& p : clebsch.a) p = Parameter::input(p.mid());

  return new Object<Poly3D>(clebsch);
}

// x² + y² + z² + x²z - y²z - 1
PTR<Object<Poly3D>> makeCayleyNodalCubic() {
  Poly3D<Parameter> cayley(3, 3, 3);

  // original: x² + y² + z² + x²z - y²z - 1
  //
  // rotated around x by pi/3:
  // c = cos(pi/3)    s = sin(pi/3)
  // x' = x
  // y' = cy - sz
  // z' = sy + cz
  //
  // (-sc²)y³ + (-cs²)z³ + (s)x²y + (c)x²z + (2cs² - c³)y²z + (2c²s - s³)yz² +
  // x² + y² + z² - 1

  const double theta = FLAGS_angle;
  const double c = cos(theta);
  const double s = sin(theta);

  cayley.set(0, 3, 0, Parameter::input(-s * c * c));
  cayley.set(0, 0, 3, Parameter::input(-c * s * s));
  cayley.set(2, 1, 0, Parameter::input(s));
  cayley.set(2, 0, 1, Parameter::input(c));
  cayley.set(0, 2, 1, Parameter::input(2 * c * s * s - c * c * c));
  cayley.set(0, 1, 2, Parameter::input(2 * c * c * s - s * s * s));
  cayley.set(2, 0, 0, Parameter::input(1));
  cayley.set(0, 2, 0, Parameter::input(1));
  cayley.set(0, 0, 2, Parameter::input(1));
  cayley.set(0, 0, 0, Parameter::input(-1));

  for (int ix = 0; ix <= 3; ++ix)
    for (int iy = 0; iy <= 3; ++iy)
      for (int iz = 0; iz <= 3; ++iz)
        if (ix + iy + iz <= 3)
          cayley.add(ix, iy, iz,
                     Parameter::input(0.0001 * ix + 0.0002 * iy + 0.0004 * iz));

  for (Parameter& p : cayley.a) p = Parameter::constant(p.mid());

  return new Object<Poly3D>(cayley);
}

PTR<Object<Poly3D>> makeSphere(double x, double y, double z, double r) {
  acp::disable();
  const double constant = -r * r + x * x + y * y + z * z;
  const double xlinear = -2 * x;
  const double ylinear = -2 * y;
  const double zlinear = -2 * z;
  acp::enable();

  // (x - dx)^2 + (y - dy)^2 + (z - dz)^2 - r^2 = 0
  // x^2 + y^2 + z^2 - 2 dx x - 2 dy y - 2 dz z + (-r^2 + dx^2 + dy^2 + dz^2)
  Poly3D<Parameter> sphere(2, 2, 2);
  sphere.set(0, 0, 0, Parameter::input(constant));
  sphere.set(1, 0, 0, Parameter::input(xlinear));
  sphere.set(0, 1, 0, Parameter::input(ylinear));
  sphere.set(0, 0, 1, Parameter::input(zlinear));
  sphere.set(2, 0, 0, Parameter::input(1));
  sphere.set(0, 2, 0, Parameter::input(1));
  sphere.set(0, 0, 2, Parameter::input(1));

  for (int ix = 0; ix <= 2; ++ix)
    for (int iy = 0; iy <= 2; ++iy)
      for (int iz = 0; iz <= 2; ++iz)
        if (ix + iy + iz <= 2)
          sphere.add(ix, iy, iz,
                     Parameter::input(0.0001 * ix + 0.0002 * iy + 0.0004 * iz));

  for (Parameter& p : sphere.a) p = Parameter::constant(p.mid());

  return new Object<Poly3D>(sphere);
}

PTR<Object<Poly3D>> makeTorus() {
  acp::disable();
  double rot = 0.1;
  double len = sqrt(1 + rot * rot);
  double c = 1 / len;
  double s = rot / len;
  acp::enable();

  // (𝑥2+𝑦2+𝑧2+𝑅2−𝑟2)2−4𝑅2(𝑥2+𝑦2)
  Parameter r = Parameter::input(0.2);
  Parameter R = Parameter::input(0.8);

  Poly3D<Parameter> pxyzRr(2, 2, 2);
  pxyzRr.set(0, 0, 0, R * R - r * r);
  pxyzRr.set(2, 0, 0, Parameter::input(1));
  pxyzRr.set(0, 2, 0, Parameter::input(1));
  // pxyzRr.set(0, 0, 2, Parameter::input(1));
  // r -> r (1 - z)
  // -r^2 = -r^2 + 2 r^2 z - r^2 z^2
  pxyzRr.set(0, 0, 1, r * r * 2);
  pxyzRr.set(0, 0, 2, Parameter::input(1) - r * r);

  Poly3D<Parameter> torus = pxyzRr * pxyzRr;
  Parameter m4R2 = R * R * -4;
  // torus.add(2, 0, 0, m4R2);
  // x^2 -> c^2 x^2 + 2 c s x y + s^2 y^2
  torus.add(2, 0, 0, m4R2 * c * c);
  torus.add(1, 1, 0, m4R2 * 2 * c * s);
  torus.add(0, 2, 0, m4R2 * s * s);
  torus.add(0, 0, 2, m4R2);

  int n = 20;
  for (int i = 0; i < n; ++i) {
    double z = -0.5 + i * 1.0 / n;
    PV3<Parameter> p(0, 0, Parameter::constant(z));
    Parameter v = torus.value(p);
    cout << z << " " << v.mid() << endl;
  }

  for (int ix = 0; ix <= 4; ++ix)
    for (int iy = 0; iy <= 4; ++iy)
      for (int iz = 0; iz <= 4; ++iz)
        if (ix + iy + iz <= 4)
          torus.add(ix, iy, iz,
                    Parameter::input(0.0001 * ix + 0.0002 * iy + 0.0004 * iz));

  for (Parameter& p : torus.a) p = Parameter::constant(p.mid());

  return new Object<Poly3D>(torus);
}

bool write_to_vtk_file(const std::string& filename,
                       const std::vector<vtkPolyData*>& meshes) {
  vtkSmartPointer<vtkAppendPolyData> appendFilter =
      vtkSmartPointer<vtkAppendPolyData>::New();

  for (auto& input : meshes) {
    appendFilter->AddInputData(input);
  }

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputConnection(appendFilter->GetOutputPort());
  writer->Update();

  writer->Write();
}

// Torus meshing and conversion to VTK.

vtkSmartPointer<vtkPolyData> mesh_to_poly_data_vtk(
    const vector<Mesh3D::Triangle>& mesh, PTR<Object<Poly3D>> poly) {
  // List of points and a map from point to index in `points`.
  vector<PTR<Object<PV3>>> points;
  map<PTR<Object<PV3>>, int> index;

  // Find all unique points and index them.
  for (Mesh3D::Triangle t : mesh) {
    for (int i = 0; i < 3; ++i) {
      PTR<Object<PV3>> p = t.p[i];
      if (index.find(p) == index.end()) {
        index[p] = points.size();
        points.push_back(p);
      }
    }
  }

  // Add all of the points.
  vtkSmartPointer<vtkPoints> points_vtk = vtkSmartPointer<vtkPoints>::New();

  for (Object<PV3>* p : points) {
    PV3<Parameter> pp = p->getApprox(1);
    points_vtk->InsertNextPoint(pp.x.mid(), pp.y.mid(), pp.z.mid());
  }

  // Get all of the vertex normals using the gradient.
  vtkSmartPointer<vtkDoubleArray> normals_array =
      vtkSmartPointer<vtkDoubleArray>::New();
  normals_array->SetNumberOfComponents(3);
  normals_array->SetNumberOfTuples(points.size());

  for (int i = 0; i < static_cast<int>(points.size()); ++i) {
    Object<PV3>* p = points[i];

    PTR<Object<PV3>> normal = new UnitNormal(poly, p);
    PV3<Parameter> normal_p = normal->get<Parameter>();

    double n[3] = {normal_p.x.mid(), normal_p.y.mid(), normal_p.z.mid()};
    normals_array->SetTuple(i, n);
  }

  // Add all of the triangles.
  vtkSmartPointer<vtkCellArray> triangles_vtk =
      vtkSmartPointer<vtkCellArray>::New();

  for (Mesh3D::Triangle t : mesh) {
    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
    triangle->GetPointIds()->SetId(0, index[t.p[0]]);
    triangle->GetPointIds()->SetId(1, index[t.p[1]]);
    triangle->GetPointIds()->SetId(2, index[t.p[2]]);
    triangles_vtk->InsertNextCell(triangle);
  }

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

  polydata->SetPoints(points_vtk);
  polydata->GetPointData()->SetNormals(normals_array);
  polydata->SetPolys(triangles_vtk);

  return polydata;
}

void calc_surface_encasement_and_mesh(vtkSmartPointer<vtkRenderer> renderer,
                                      const double extent,
                                      const double epsilon) {
  acp::enable();

  PV3<Parameter> box(
      Parameter::constant(-extent).interval(Parameter::constant(extent)),
      Parameter::constant(-extent).interval(Parameter::constant(extent)),
      Parameter::constant(-extent).interval(Parameter::constant(extent)));

  std::vector<PTR<Object<Poly3D>>> polys;

  if (FLAGS_surface.find("torus") != std::string::npos) {
    polys.push_back(makeTorus());
  }
  if (FLAGS_surface.find("cayley") != std::string::npos) {
    polys.push_back(makeCayleyNodalCubic());
  }
  if (FLAGS_surface.find("clebsch") != std::string::npos) {
    polys.push_back(makeClebschDiagonalCubic());
  }
  if (FLAGS_surface.find("steiner") != std::string::npos) {
    polys.push_back(makeSteinerRoman());
  }
  if (FLAGS_surface.find("hunt") != std::string::npos) {
    polys.push_back(makeHunt());
  }
  if (FLAGS_surface.find("spheres") != std::string::npos) {
    switch (FLAGS_num_spheres) {
      case 4:
        polys.push_back(makeSphere(-0.5, 0.4, -0.8, 0.9));
      case 3:
        polys.push_back(makeSphere(0.5, -0.4, 0.4, 0.9));
      case 2:
        polys.push_back(makeSphere(0.4, -0.5, 0.8, 0.9));
      case 1:
        polys.push_back(makeSphere(-0.4, 0.5, -0.4, 0.9));
    }
  }
  Timer timer;
  timer.Mark("start");
  Encasement3D* encasement = new Encasement3D(box, polys);
  timer.Mark("encasement");

  int triangle_count = 0;
  int face_count = 0;
  class Mesh3D mesh3d(encasement, epsilon);

  // Loop over all Encasement3D::F objects (except for the boundary faces).
  for (int h = 6; h < encasement->f.size(); ++h) {
    // Loop over all faces for this F.
    for (int i = 0; i < encasement->f[h]->faces.size(); ++i) {
      // Calculate the mesh.
      const vector<Mesh3D::Triangle>& mesh =
          mesh3d.getMesh(encasement->f[h]->faces[i]);
    }
  }

  // Loop over all Encasement3D::F objects (except for the boundary faces).
  for (int h = 6; h < encasement->f.size(); ++h) {
    // Loop over all faces for this F.
    for (int i = 0; i < encasement->f[h]->faces.size(); ++i) {
      // Calculate the mesh.
      const vector<Mesh3D::Triangle>& mesh =
          mesh3d.getMesh(encasement->f[h]->faces[i]);

      face_count++;
      triangle_count += mesh.size();

      face_to_acp_mesh[encasement->f[h]->faces[i]] = mesh;
    }
  }
  timer.Mark("meshing");

  std::cout << face_count << " faces\n";
  std::cout << triangle_count << " triangles\n";
  timer.Report();

  // Loop over all Encasement3D::F objects (except for the boundary faces).
  for (int h = 6; h < encasement->f.size(); ++h) {
    // Loop over all faces for this F.
    for (int i = 0; i < encasement->f[h]->faces.size(); ++i) {
      const vector<Mesh3D::Triangle>& mesh =
          face_to_acp_mesh[encasement->f[h]->faces[i]];

      // mesh3d to vtk poly data.
      vtkSmartPointer<vtkPolyData> mesh_poly_data =
          mesh_to_poly_data_vtk(mesh, encasement->bf[encasement->f[h]->jf]);

      face_to_mesh[encasement->f[h]->faces[i]] = mesh_poly_data;

      // Visualize Polys with vertex normals computed.
      // vtkSmartPointer<vtkPolyDataNormals> normalGenerator =
      //     vtkSmartPointer<vtkPolyDataNormals>::New();
      // normalGenerator->SetInputData(mesh_poly_data);
      // normalGenerator->ComputePointNormalsOn();
      // normalGenerator->ComputeCellNormalsOn();
      // normalGenerator->Update();

      vtkSmartPointer<vtkPolyDataMapper> polyMapper =
          vtkSmartPointer<vtkPolyDataMapper>::New();
      // polyMapper->SetInputConnection(normalGenerator->GetOutputPort());
      polyMapper->SetInputData(mesh_poly_data);
      vtkSmartPointer<vtkActor> polyActor = vtkSmartPointer<vtkActor>::New();
      polyActor->SetMapper(polyMapper);

      polyActor->GetProperty()->SetColor(colors[h % 6]);

      actor_to_color[polyActor] = h % 6;

      // Establish face to actor mapping.
      face_to_actor[encasement->f[h]->faces[i]] = polyActor;
      actor_to_face[polyActor] = encasement->f[h]->faces[i];

      // Visualize Edges
      vtkSmartPointer<vtkFeatureEdges> featureEdges =
          vtkSmartPointer<vtkFeatureEdges>::New();
      featureEdges->SetInputData(mesh_poly_data);
      featureEdges->BoundaryEdgesOn();
      featureEdges->FeatureEdgesOff();
      featureEdges->ManifoldEdgesOff();
      featureEdges->NonManifoldEdgesOff();
      featureEdges->Update();

      vtkSmartPointer<vtkPolyDataMapper> edgeMapper =
          vtkSmartPointer<vtkPolyDataMapper>::New();
      edgeMapper->SetInputConnection(featureEdges->GetOutputPort());
      vtkSmartPointer<vtkActor> edgeActor = vtkSmartPointer<vtkActor>::New();
      edgeActor->SetMapper(edgeMapper);
      edgeActor->GetProperty()->SetLineWidth(3);

      // Add Actors.
      renderer->AddActor(polyActor);
      renderer->AddActor(edgeActor);
    }
  }

  // For each Cell of the arrangement, collect all the faces and write them to a
  // single VTP file.
  int cell_count = 0;
  for (auto& cell : encasement->cells) {
    std::vector<vtkPolyData*> meshes;
    for (auto& shell : cell->shells) {
      for (auto& face : shell->faces) {
        meshes.push_back(face_to_mesh[face]);
      }
    }

    std::string filename = "cell_" + std::to_string(cell_count++) + ".vtp";

    write_to_vtk_file(filename, meshes);
  }

  // Alternatively, write all individual faces to VTP files.
  face_count = 0;
  for (auto& [face, mesh] : face_to_mesh) {
    std::string filename = "face_" + std::to_string(face_count++) + ".vtp";

    write_to_vtk_file(filename, {mesh});
  }

  acp::disable();
}

int main(int argc, char** argv) {
  gflags::SetUsageMessage(
      "\n  Run a single surface encasement, meshing and visualization."
      "\n  Use the options listed below to select the surface and box "
      "parameters.");

  gflags::ParseCommandLineFlags(&argc, &argv, true);

  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

  // Read the background image.
  std::string exe_path = argv[0];
  std::string bg_path = exe_path.substr(0, exe_path.find_last_of("/")) + "/" +
                        FLAGS_background_image + ".jpg";

  std::cout << "bg path: " << bg_path << "\n";

  vtkSmartPointer<vtkJPEGReader> jpegReader =
      vtkSmartPointer<vtkJPEGReader>::New();
  if (!jpegReader->CanReadFile(bg_path.c_str())) {
    std::cerr << "Error reading file " << bg_path << std::endl;
    return EXIT_FAILURE;
  }
  jpegReader->SetFileName(bg_path.c_str());
  jpegReader->Update();
  vtkImageData* imageData = jpegReader->GetOutput();

  vtkSmartPointer<vtkImageActor> imageActor =
      vtkSmartPointer<vtkImageActor>::New();
  imageActor->SetInputData(imageData);

  vtkSmartPointer<vtkRenderer> backgroundRenderer =
      vtkSmartPointer<vtkRenderer>::New();
  backgroundRenderer->AddActor(imageActor);

  // Create the render window.
  vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
  backgroundRenderer->SetLayer(0);
  backgroundRenderer->InteractiveOff();
  renderer->SetLayer(1);
  renderWindow->SetNumberOfLayers(2);
  renderWindow->AddRenderer(backgroundRenderer);
  renderWindow->AddRenderer(renderer);
  renderWindow->SetSize(2000, 1500);

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  vtkSmartPointer<HighlightFaceStyle> style =
      vtkSmartPointer<HighlightFaceStyle>::New();
  style->SetDefaultRenderer(renderer);
  renderWindowInteractor->SetInteractorStyle(style);

  // Calculate the surface encasement and mesh it, adding meshes to the
  // renderer.
  calc_surface_encasement_and_mesh(renderer, FLAGS_extent, FLAGS_mesh_epsilon);

  // Render once to figure out where the background camera will be
  renderWindow->Render();

  // Set up the background camera to fill the renderer with the image
  double origin[3];
  double spacing[3];
  int extent[6];
  imageData->GetOrigin(origin);
  imageData->GetSpacing(spacing);
  imageData->GetExtent(extent);

  vtkCamera* camera = backgroundRenderer->GetActiveCamera();
  camera->ParallelProjectionOn();

  double xc = origin[0] + 0.5 * (extent[0] + extent[1]) * spacing[0];
  double yc = origin[1] + 0.5 * (extent[2] + extent[3]) * spacing[1];
  // double xd = (extent[1] - extent[0] + 1)*spacing[0];
  double yd = (extent[3] - extent[2] + 1) * spacing[1];
  double d = camera->GetDistance();
  camera->SetParallelScale(0.5 * yd);
  camera->SetFocalPoint(xc, yc, 0.0);
  camera->SetPosition(xc, yc, d);

  // Render again to set the correct view
  renderWindow->Render();
  renderWindowInteractor->Initialize();
  renderWindowInteractor->Start();

  face_to_actor.clear();
  actor_to_face.clear();
  face_to_mesh.clear();

  gflags::ShutDownCommandLineFlags();
  return 0;
}
