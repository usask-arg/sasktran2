project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

python_add_library(
  _core
  MODULE
  main.cpp
  config.cpp
  coordinates.cpp
  derivative_mappings.cpp
  geometry.cpp
  geodetic.cpp
  grids.cpp
  surface.cpp
  atmosphere.cpp
  viewing_geometry.cpp
  engine.cpp
  mie.cpp
  output.cpp
  math.cpp
  WITH_SOABI)

target_link_libraries(_core PRIVATE pybind11::headers)
target_link_libraries(_core PRIVATE sasktran2)

install(TARGETS _core DESTINATION sasktran2)
