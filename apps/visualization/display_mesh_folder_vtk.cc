#include <vtkActor.h>
#include <vtkDirectory.h>
#include <vtkDiskSource.h>
#include <vtkFeatureEdges.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNamedColors.h>
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
#include <vtksys/SystemTools.hxx>

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
        picker->GetActor()->GetProperty()->SetColor(white[0], white[1],
                                                    white[2]);
      }
    }
    // Forward events
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
  }

  virtual void OnRightButtonDown() override {
    int* clickPos = this->GetInteractor()->GetEventPosition();

    vtkSmartPointer<vtkPropPicker> picker =
        vtkSmartPointer<vtkPropPicker>::New();
    picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

    if (picker->GetActor()) {
      picker->GetActor()->SetVisibility(false);
    }
    // Forward events
    vtkInteractorStyleTrackballCamera::OnRightButtonDown();
  }

 private:
  const double selected[3] = {0.0, 0.6, 0.8};
  const double white[3] = {1, 1, 1};
};

vtkStandardNewMacro(HighlightFaceStyle);

// Load all vtk files in a directory.
bool load_all_vtk_files(std::string directory_name,
                        vtkSmartPointer<vtkRenderer> renderer) {
  vtkSmartPointer<vtkDirectory> directory =
      vtkSmartPointer<vtkDirectory>::New();
  int opened = directory->Open(directory_name.c_str());

  if (!opened) {
    std::cerr << "Invalid directory!\n";
    return false;
  }

  for (int i = 0; i < directory->GetNumberOfFiles(); ++i) {
    std::string file_string = directory_name;
    file_string += "/";
    file_string += directory->GetFile(i);

    std::cout << file_string << "\n";

    if (vtksys::SystemTools::GetFilenameLastExtension(file_string) == ".vtk") {
      auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
      reader->SetFileName(file_string.c_str());
      reader->Update();

      // Visualize Polys with vertex normals computed.
      vtkSmartPointer<vtkPolyDataNormals> normalGenerator =
          vtkSmartPointer<vtkPolyDataNormals>::New();
      normalGenerator->SetInputConnection(reader->GetOutputPort());
      normalGenerator->ComputePointNormalsOn();
      normalGenerator->ComputeCellNormalsOn();
      normalGenerator->Update();

      vtkSmartPointer<vtkPolyDataMapper> polyMapper =
          vtkSmartPointer<vtkPolyDataMapper>::New();
      polyMapper->SetInputConnection(normalGenerator->GetOutputPort());
      vtkSmartPointer<vtkActor> polyActor = vtkSmartPointer<vtkActor>::New();
      polyActor->SetMapper(polyMapper);

      // Visualize Edges
      vtkSmartPointer<vtkFeatureEdges> featureEdges =
          vtkSmartPointer<vtkFeatureEdges>::New();
      featureEdges->SetInputConnection(reader->GetOutputPort());
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
      edgeActor->GetProperty()->SetLineWidth(2);

      // Add Actors.
      renderer->AddActor(polyActor);
      renderer->AddActor(edgeActor);
    }
  }
}

int main(int argc, char** argv) {
  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

  vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  renderWindow->SetSize(2000, 1500);

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  vtkSmartPointer<HighlightFaceStyle> style =
      vtkSmartPointer<HighlightFaceStyle>::New();
  style->SetDefaultRenderer(renderer);
  renderWindowInteractor->SetInteractorStyle(style);

  // Load the directory of .vtk files.
  load_all_vtk_files(argv[1], renderer);

  renderWindow->Render();
  renderWindowInteractor->Initialize();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
