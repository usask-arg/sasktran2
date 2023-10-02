cmake_minimum_required(VERSION 3.4...3.18)
project(_pysasktran2)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(pybind11)
pybind11_add_module(_pysasktran2
        src/main.cpp
        src/config.cpp
        src/coordinates.cpp
        src/geometry.cpp
        src/grids.cpp
        src/atmosphere.cpp
        src/viewing_geometry.cpp
        src/engine.cpp
        src/output.cpp
        )

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(_pysasktran2
        PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})

# Check if sasktran2 is already installed
find_package(sasktran2 CONFIG)

if(sasktran2_FOUND)
    # If so we can just link it in.  This is the path taken by the conda build where the c++ library is separated
    # from the python wheel
    target_link_libraries(_pysasktran2 PRIVATE sasktran2::sasktran2)
else()
    # Alternatively we include sasktran2 directly and build it along side pysasktran2
    # This lets us just call pip install from the Python folder
    if(PROJECT_IS_TOP_LEVEL)
        set(FORCE_BOOST_STATIC OFF CACHE BOOL "Disable Boost forced static" FORCE)
        add_subdirectory(../ ${CMAKE_CURRENT_BINARY_DIR}/sasktran2)
    endif()
    target_link_libraries(_pysasktran2 PRIVATE sasktran2)
endif()

