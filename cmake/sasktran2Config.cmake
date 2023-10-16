include("${CMAKE_CURRENT_LIST_DIR}/sasktran2Targets.cmake")

include(CMakeFindDependencyMacro)

find_dependency(OpenMP)
find_dependency(Eigen3)

find_dependency(BLAS)
find_dependency(LAPACK)

# optional
find_package(LAPACKE CONFIG)
