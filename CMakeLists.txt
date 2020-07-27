cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(Molding)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(VTK COMPONENTS
  vtkCommonCore
  vtkCommonDataModel
  vtkIOGeometry
  vtkFiltersCore
  vtkFiltersExtraction
  vtkFiltersGeneral
  vtkFiltersSources
  vtkInteractionStyle
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingOpenGL2 QUIET)
if (NOT VTK_FOUND)
    message("Skipping Molding: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(molding MACOSX_BUNDLE molding.cxx )
  target_link_libraries(molding PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(molding MACOSX_BUNDLE molding.cxx )
  target_link_libraries(molding PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS molding
    MODULES ${VTK_LIBRARIES}
    )
endif ()
