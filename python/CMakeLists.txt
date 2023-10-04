cmake_minimum_required(VERSION 3.4...3.18)
project(
        ${SKBUILD_PROJECT_NAME}
        VERSION ${SKBUILD_PROJECT_VERSION}
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

python_add_library(_core MODULE
        src/main.cpp
        src/config.cpp
        src/coordinates.cpp
        src/geometry.cpp
        src/grids.cpp
        src/atmosphere.cpp
        src/viewing_geometry.cpp
        src/engine.cpp
        src/output.cpp
        WITH_SOABI)

target_link_libraries(_core PRIVATE pybind11::headers)
target_link_libraries(_core PRIVATE sasktran2)
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS _core DESTINATION sasktran2)



