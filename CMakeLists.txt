cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(Molding)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# https://cmake.org/cmake/help/latest/prop_tgt/CXX_STANDARD.html
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(molding MACOSX_BUNDLE molding.cxx)

find_package (Threads)
target_link_libraries (molding PUBLIC Threads::Threads)

find_package(
  VTK
  COMPONENTS CommonCore
             CommonDataModel
             IOGeometry
             IOPLY
             IOXML
             FiltersCore
             FiltersExtraction
             FiltersGeneral
             FiltersSources
             InteractionStyle
             RenderingCore
             RenderingFreeType
             RenderingOpenGL2
             RenderingGL2PSOpenGL2
  QUIET)
if(NOT VTK_FOUND)
  message("Skipping Molding: ${VTK_NOT_FOUND_MESSAGE}")
  return()
endif()
message(STATUS "VTK_VERSION: ${VTK_VERSION}")
if(VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  target_link_libraries(molding PRIVATE ${VTK_LIBRARIES})
else()
  # include all components
  target_link_libraries(molding PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(TARGETS molding MODULES ${VTK_LIBRARIES})
endif()
