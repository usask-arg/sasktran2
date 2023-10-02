include( "${CMAKE_CURRENT_LIST_DIR}/sasktran2Targets.cmake")

include (CMakeFindDependencyMacro)

find_package(Boost REQUIRED COMPONENTS log)
find_dependency(OpenMP)
find_dependency(Eigen3)

find_dependency(BLAS)
find_dependency(LAPACK)