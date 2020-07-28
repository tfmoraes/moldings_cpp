#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>
#include <vtkAbstractPolyDataReader.h>
#include <vtkCleanPolyData.h>
#include <vtkIdList.h>
#include <vtkPLYReader.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkTriangleFilter.h>

using G_type = std::unordered_map<int, std::vector<int>>;

struct Cell {
  int id;
  int plane;
  float area;
  float center[3];
  float normals[3];
};

struct Plane {
  int id;
  float center[3];
  float normals[3];
  std::vector<Cell> cells;
};

// https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
inline bool ends_with(std::string const &value, std::string const &ending) {
  if (ending.size() > value.size())
    return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

int get_cell_edge_neighbours(vtkSmartPointer<vtkPolyData> polydata, int cell_id,
                             int v0, int v1) {
  int cell_neighbour;
  auto cells = vtkSmartPointer<vtkIdList>::New();
  polydata->GetCellEdgeNeighbors(cell_id, v0, v1, cells);
  assert(cells->GetNumberOfIds() == 1);
  for (int i = 0; i < cells->GetNumberOfIds(); i++) {
    cell_neighbour = cells->GetId(i);
    if (cell_neighbour != cell_id) {
      return cell_neighbour;
    }
  }
}

std::unique_ptr<G_type>
build_cell_connectivity_graph(vtkSmartPointer<vtkPolyData> polydata) {
  auto G = std::make_unique<G_type>();
  int number_cells = polydata->GetNumberOfCells();
  for (int cell_id = 0; cell_id < number_cells; cell_id++) {
    auto points = vtkSmartPointer<vtkIdList>::New();
    polydata->GetCellPoints(cell_id, points);
    assert(points->GetNumberOfIds() == 3);
    int p0 = points->GetId(0);
    int p1 = points->GetId(1);
    int p2 = points->GetId(2);

    int cell_neighbour_0 = get_cell_edge_neighbours(polydata, cell_id, p0, p1);
    int cell_neighbour_1 = get_cell_edge_neighbours(polydata, cell_id, p1, p2);
    int cell_neighbour_2 = get_cell_edge_neighbours(polydata, cell_id, p2, p0);

    (*G)[cell_id].push_back(cell_neighbour_0);
    (*G)[cell_neighbour_0].push_back(cell_id);

    (*G)[cell_id].push_back(cell_neighbour_0);
    (*G)[cell_neighbour_1].push_back(cell_id);

    (*G)[cell_id].push_back(cell_neighbour_2);
    (*G)[cell_neighbour_2].push_back(cell_id);
  }
  return G;
}

vtkSmartPointer<vtkPolyData> read_file(const std::string &fname) {
  vtkSmartPointer<vtkAbstractPolyDataReader> reader;
  if (ends_with(fname, ".stl")) {
    reader = vtkSmartPointer<vtkSTLReader>::New();
  } else {
    reader = vtkSmartPointer<vtkPLYReader>::New();
  }
  reader->SetFileName(fname.c_str());
  reader->Update();

  auto triangle = vtkTriangleFilter::New();
  triangle->SetInputConnection(reader->GetOutputPort());
  triangle->PassVertsOff();
  triangle->Update();

  auto clean = vtkCleanPolyData::New();
  clean->SetInputConnection(triangle->GetOutputPort());
  clean->Update();

  auto normals = vtkPolyDataNormals::New();
  normals->SetInputConnection(clean->GetOutputPort());
  normals->ComputePointNormalsOn();
  normals->SplittingOff();
  normals->ComputeCellNormalsOn();
  normals->Update();

  auto polydata = normals->GetOutput();
  polydata->BuildLinks();
  return polydata;
}

int main(int argc, char *argv[]) {
  std::string input_filename = argv[1];
  std::string output_filename = argv[2];

  auto polydata = read_file(argv[1]);
  std::unique_ptr<G_type> G = build_cell_connectivity_graph(polydata);

  std::cout << polydata << std::endl;
}
