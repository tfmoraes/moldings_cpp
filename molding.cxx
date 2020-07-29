#include <array>
#include <cassert>
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <memory>
#include <random>
#include <unordered_map>
#include <vector>
#include <vtkAbstractPolyDataReader.h>
#include <vtkCellData.h>
#include <vtkCleanPolyData.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkCellData.h>
#include <vtkPLYReader.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>
#include <vtkXMLPolyDataWriter.h>

using G_type = std::vector<std::vector<int>>;

const int WN = 10;
const double WV = 0.05;

struct Cell {
  int id;
  int plane;
  double area;
  std::array<double, 3> center;
  std::array<double, 3> normals;
};

struct Plane {
  int id;
  double area;
  std::array<double, 3> center;
  std::array<double, 3> normals;
  std::list<int> cells;
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
  return -1;
}

void build_cell_connectivity_graph(vtkSmartPointer<vtkPolyData> polydata, G_type &G) {
  int number_cells = polydata->GetNumberOfCells();
  G.resize(polydata->GetNumberOfCells());
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

    if (cell_neighbour_0 != -1) {
      G[cell_id].push_back(cell_neighbour_0);
      // G[cell_neighbour_0].push_back(cell_id);
    }

    if (cell_neighbour_1 != -1) {
      G[cell_id].push_back(cell_neighbour_1);
      // G[cell_neighbour_1].push_back(cell_id);
    }

    if (cell_neighbour_2 != -1) {
      G[cell_id].push_back(cell_neighbour_2);
      // G[cell_neighbour_2].push_back(cell_id);
    }
  }
}

std::array<double, 3> calc_cell_center(vtkSmartPointer<vtkPolyData> polydata, int cell_id) {
  auto points = vtkSmartPointer<vtkIdList>::New();
  polydata->GetCellPoints(cell_id, points);
  std::array<double, 3> center = {0, 0, 0};
  for (int i = 0; i < points->GetNumberOfIds(); i++) {
    auto vertice = polydata->GetPoint(points->GetId(i));
    center[0] += vertice[0];
    center[1] += vertice[1];
    center[2] += vertice[2];
  }

  center[0] /= 3.0;
  center[1] /= 3.0;
  center[2] /= 3.0;

  return center;
}

double calc_cell_area(vtkSmartPointer<vtkPolyData> polydata, int cell_id) {
  vtkCell *cell = polydata->GetCell(cell_id);
  vtkTriangle *triangle = dynamic_cast<vtkTriangle *>(cell);
  return triangle->ComputeArea();
}

void populate_cells(vtkSmartPointer<vtkPolyData> polydata, std::vector<Cell> &cells) {
  auto normals = polydata->GetCellData()->GetArray("Normals");
  cells.reserve(polydata->GetNumberOfCells());
  for (int cell_id = 0; cell_id < polydata->GetNumberOfCells(); cell_id++) {
    double area = calc_cell_area(polydata, cell_id);
    auto normal = normals->GetTuple(cell_id);
    auto center = calc_cell_center(polydata, cell_id);
    cells.push_back(Cell{cell_id, -1, area, center, {normal[0], normal[1], normal[2]}});
  }
}

void init_labeling(vtkSmartPointer<vtkPolyData> polydata, std::vector<Cell> &cells, std::vector<Plane> &planes, G_type G, int number_planes) {
  std::default_random_engine generator;
  std::uniform_int_distribution<size_t> distribution(0, polydata->GetNumberOfCells() - 1);
  std::deque<int> stack;
  int seed;
  int cell_id;
  for (int i = 0; i < number_planes; i++) {
    seed = distribution(generator);
    planes.push_back(Plane{});
    stack.push_back(seed);
    cells[seed].plane = i;
  }

  while (stack.size()) {
    cell_id = stack.front();
    // std::cout << cell_id << std::endl;
    stack.pop_front();
    for (auto j : G[cell_id]) {
      if (cells[j].plane == -1) {
        cells[j].plane = cells[cell_id].plane;
        stack.push_back(j);
      }
    }
  }
}

double calc_manhathan_distance(std::array<double, 3> v0, std::array<double, 3> v1) {
  return std::fabs(v0[0] - v1[0]) + std::fabs(v0[1] - v1[1]) + std::fabs(v0[2] - v1[2]);
}

std::array<double, 3> arr_diff(std::array<double, 3> v0, std::array<double, 3> v1) {
  return {(v0[0] - v1[0]), (v0[1] - v1[1]), (v0[2] - v1[2])};
}

double dot(std::array<double, 3> v0, std::array<double, 3> v1) {
  return (v0[0] * v1[0]) + (v0[1] * v1[1]) + (v0[2] * v1[2]);
}

double calc_energy(G_type &G, std::vector<Cell> &cells, std::vector<Plane> &planes, std::vector<int> plane_ids) {
  double energy = 0.0;
  double dp = 0.0;
  double vpq = 0.0;
  for (auto plane_id : plane_ids) {
    auto &plane = planes[plane_id];
    for (auto cell_id : plane.cells) {
      auto &cell = cells[cell_id];
      dp += (cell.area * (calc_manhathan_distance(plane.normals, cell.normals) + WN * dot(arr_diff(plane.center, cell.center), plane.normals)));
      for (auto j : G[cell_id]) {
        auto &cell_neighbour = cells[j];
        if (cell.plane == cell_neighbour.plane) {
          vpq += 0;
        } else {
          for (auto p : plane_ids) {
            if (cell_neighbour.plane == p) {
              vpq += (((cell.area + cell_neighbour.area) / 2.0) * (4 - calc_manhathan_distance(cell.normals, cell_neighbour.normals)));
            }
          }
        }
      }
    }
  }
  energy = dp + WV * vpq;
  return energy;
}

int update_planes(std::vector<Plane> &planes, std::vector<Cell> &cells) {
  std::vector old_planes = planes;
  int updated = 0;
  for (auto &plane : planes) {
    plane.area = 0.0;
    plane.center = {0, 0, 0};
    plane.normals = {0, 0, 0};
    plane.cells.clear();
  }

  for (auto &cell : cells) {
    planes[cell.plane].area += cell.area;

    planes[cell.plane].center[0] += (cell.center[0] * cell.area);
    planes[cell.plane].center[1] += (cell.center[1] * cell.area);
    planes[cell.plane].center[2] += (cell.center[2] * cell.area);

    planes[cell.plane].normals[0] += (cell.normals[0] * cell.area);
    planes[cell.plane].normals[1] += (cell.normals[1] * cell.area);
    planes[cell.plane].normals[2] += (cell.normals[2] * cell.area);

    planes[cell.plane].cells.push_back(cell.id);
  }

  for (auto &plane : planes) {
    auto &old_plane = old_planes[plane.id];
    plane.center[0] /= plane.area;
    plane.center[1] /= plane.area;
    plane.center[2] /= plane.area;

    plane.normals[0] /= plane.area;
    plane.normals[1] /= plane.area;
    plane.normals[2] /= plane.area;

    if (
        (plane.center[0] != old_plane.center[0]) &&
        (plane.center[1] != old_plane.center[1]) &&
        (plane.center[2] != old_plane.center[2]) &&
        (plane.normals[0] != old_plane.normals[0]) &&
        (plane.normals[1] != old_plane.normals[1]) &&
        (plane.normals[2] != old_plane.normals[2])) {
      updated = 1;
    }
  }
  return updated;
}

void swap_optimize(G_type &G, std::vector<Cell> &cells, std::vector<Plane> &planes) {
  double energy, new_energy;
  int modified = 0;
  int updated = 0;
  int step = 0;
  int sub_step = 0;
  do {
    do {
        std::cout << "step:" << step << "." << sub_step << std::endl;
        modified = 0;
      for (int p0 = 0; p0 < planes.size(); p0++) {
        for (int p1 = 0; p1 < planes.size(); p1++) {
          if (p0 == p1) {
            continue;
          }
          energy = calc_energy(G, cells, planes, {p0, p1});
          auto plane_0 = planes[p0];
          auto plane_1 = planes[p1];
          for (auto cell_id : plane_0.cells) {
            auto &cell = cells[cell_id];
            for (auto j : G[cell_id]) {
              auto &cell_neighbour = cells[j];
              if (cell_neighbour.plane == p1) {
                planes[p1].cells.remove(j);
                planes[p0].cells.push_back(j);
                new_energy = calc_energy(G, cells, planes, {p0, p1});
                if (new_energy < energy) {
                  energy = new_energy;
                  modified = 1;
                } else {
                  planes[p1].cells.push_back(j);
                  planes[p0].cells.pop_back();
                }
              }
            }
          }
        }
      }
      sub_step++;
    } while (modified);
    updated = update_planes(planes, cells);
    step++;
  } while (updated);
}

vtkSmartPointer<vtkPolyData> gen_output_polydata(vtkSmartPointer<vtkPolyData> polydata, std::vector<Cell> cells) {
  auto output_polydata = vtkSmartPointer<vtkPolyData>::New();
  output_polydata->DeepCopy(polydata);

  auto color_array1 = vtkSmartPointer<vtkIdTypeArray>::New();
  color_array1->SetName("RegionID1");
  color_array1->SetNumberOfValues(polydata->GetNumberOfCells());

  auto color_array2 = vtkSmartPointer<vtkIdTypeArray>::New();
  color_array2->SetName("RegionID2");
  color_array2->SetNumberOfValues(polydata->GetNumberOfCells());

  for (auto &cell : cells) {
    color_array1->SetValue(cell.id, cell.plane);
    color_array2->SetValue(cell.id, cell.plane);
  }
  output_polydata->GetCellData()->AddArray(color_array1);
  output_polydata->GetCellData()->AddArray(color_array2);
  return output_polydata;
}

void update_polydata(vtkSmartPointer<vtkPolyData> polydata, std::vector<Cell> cells){
  auto color_array2 = dynamic_cast<vtkIdTypeArray *>(polydata->GetCellData()->GetArray("RegionID2"));
  for (auto &cell : cells) {
    color_array2->SetValue(cell.id, cell.plane);
  }
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
  int number_planes = std::stoi(argv[3]);

  std::cout << number_planes << std::endl;

  std::vector<Cell> cells;
  std::vector<Plane> planes;
  planes.reserve(number_planes);
  G_type G;

  auto polydata = read_file(argv[1]);
  build_cell_connectivity_graph(polydata, G);
  populate_cells(polydata, cells);
  init_labeling(polydata, cells, planes, G, number_planes);
  update_planes(planes, cells);
  auto output_polydata = gen_output_polydata(polydata, cells);
  swap_optimize(G, cells, planes);

  std::cout << cells.size() << std::endl;
  std::cout << polydata << std::endl;

  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(output_filename.c_str());
  writer->SetInputData(output_polydata);
  writer->Write();
}
