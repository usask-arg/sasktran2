cmake_minimum_required(VERSION 3.20.0)
include(CMakePackageConfigHelpers) # Provides function/macro
                                   # write_basic_package_version_file
include(ExternalProject)

# Convenience
if(UNIX AND NOT APPLE)
  set(LINUX TRUE)
endif()

# Export all symbols on windows
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

if(UNIX)
  add_definitions(
    "-Wno-inconsistent-missing-override -Wno-suggest-override -Wno-comment -Wno-deprecated-declarations -Wno-psabi"
  )
  enable_language(Fortran)
endif()

execute_process(
  COMMAND git describe --tags
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_DESCRIBE_TEMP)

if(EXISTS ${GIT_DESCRIBE_TEMP})
  string(REGEX MATCH "v([0-9\\.]+)" SASKTRAN2_VERSION ${GIT_DESCRIBE_TEMP})
  string(REPLACE "v" "" SASKTRAN2_VERSION ${GIT_DESCRIBE_TEMP})
else()
  set(SASKTRAN2_VERSION 0.0.0)
endif()

project(sasktran2 VERSION ${SASKTRAN2_VERSION})

set(CMAKE_CXX_STANDARD 17)

if(APPLE)
  # This is nominally a windows definition but for some reason BOOST < 1.8.3 has
  # problems on clang without this?
  add_compile_definitions(_HAS_AUTO_PTR_ETC=0)
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add an interface target for storing build properties common to all sasktran
# modules
add_library(sasktran2BuildProperties INTERFACE)

find_package(Eigen3 CONFIG)

if(Eigen3_FOUND)
  target_link_libraries(sasktran2BuildProperties INTERFACE Eigen3::Eigen)
else()
  set(EIGEN_INSTALL_DIR "${CMAKE_BINARY_DIR}/eigen-install/")

  ExternalProject_Add(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    SOURCE_DIR "${CMAKE_BINARY_DIR}/eigen-src"
    BINARY_DIR "${CMAKE_BINARY_DIR}/eigen-build"
    INSTALL_DIR "${EIGEN_INSTALL_DIR}"
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
               -DCMAKE_BUILD_TYPE=Release)

  file(MAKE_DIRECTORY ${EIGEN_INSTALL_DIR}/include) # avoid race condition

  add_library(eigenlib INTERFACE IMPORTED GLOBAL)
  add_dependencies(eigenlib eigen)

  set_target_properties(eigenlib PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                            ${EIGEN_INSTALL_DIR}/include/eigen3)

  target_link_libraries(sasktran2BuildProperties INTERFACE eigenlib)
endif()

find_package(spdlog CONFIG)

if(spdlog_FOUND)
  target_link_libraries(sasktran2BuildProperties
                        INTERFACE spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>)
else()
  set(SPDLOG_INSTALL_DIR "${CMAKE_BINARY_DIR}/spdlog-install/")

  ExternalProject_Add(
    spdlog
    URL https://github.com/gabime/spdlog/archive/refs/tags/v1.12.0.tar.gz
    SOURCE_DIR "${CMAKE_BINARY_DIR}/spdlog-src"
    BINARY_DIR "${CMAKE_BINARY_DIR}/spdlog-build"
    INSTALL_DIR "${SPDLOG_INSTALL_DIR}"
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
               -DCMAKE_BUILD_TYPE=Release)

  file(MAKE_DIRECTORY ${SPDLOG_INSTALL_DIR}/include) # avoid race condition

  add_library(spdloglib INTERFACE IMPORTED GLOBAL)
  add_dependencies(spdloglib spdlog)

  set_target_properties(spdloglib PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                             ${SPDLOG_INSTALL_DIR}/include/)

  target_link_libraries(sasktran2BuildProperties INTERFACE spdloglib)
endif()

# add the additional compiler flags
set(CMAKE_CXX_FLAGS " ${CompilerFlags} ${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Directory to install the static library files to
set(STATIC_LIB_INSTALL_DIR lib)
set(INCLUDE_INSTALL_DIR include)
set(EXPORT_NAME sasktran2Targets)

# Try to find OpenMP
option(USE_OMP "Enables OPENMP support in SASKTRAN" ON)

if(USE_OMP)
  find_package(OpenMP REQUIRED)

  if(OpenMP_FOUND)
    target_link_libraries(sasktran2BuildProperties INTERFACE OpenMP::OpenMP_CXX)
    target_compile_definitions(sasktran2BuildProperties
                               INTERFACE SKTRAN_OPENMP_SUPPORT)
  endif()
endif()

option(DO_STREAM_TEMPLATES
       "Enables Compilation of templated number of streams in DO" OFF)

if(DO_STREAM_TEMPLATES)
  target_compile_definitions(sasktran2BuildProperties
                             INTERFACE SASKTRAN_DISCO_FULL_COMPILE)
  message(STATUS "Enabling Full DO Stream Template Compilation")
endif()

# -------- LAPACK FINDING -----------
# We support multiple blas vendors, but these can require some funny options
set(SKTRAN_BLAS_VENDOR
    "OpenBLAS"
    CACHE
      STRING
      "Sets the BLAS/LAPACK vendor. See https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors."
)
set_property(CACHE SKTRAN_BLAS_VENDOR PROPERTY STRINGS OpenBLAS Intel10_64lp
                                               Apple Generic)

if(DEFINED ENV{CONDA_BUILD_STATE})
  # we are inside conda build, force SKTRAN_BLAS_VENDOR to be Generic
  set(SKTRAN_BLAS_VENDOR "Generic")
endif()

if(SKTRAN_BLAS_VENDOR MATCHES "Intel")
  set(BLA_STATIC FALSE)
  set(BLA_SIZEOF_INTEGER 4)
elseif(SKTRAN_BLAS_VENDOR MATCHES "Apple")
  set(BLA_STATIC TRUE)
  set(BLA_SIZEOF_INTEGER 4)
  target_include_directories(
    sasktran2BuildProperties
    INTERFACE
      /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers
  )
  target_compile_definitions(sasktran2BuildProperties
                             INTERFACE SKTRAN_USE_ACCELERATE)
  target_compile_definitions(sasktran2BuildProperties
                             INTERFACE ACCELERATE_NEW_LAPACK=1)
  target_compile_definitions(sasktran2BuildProperties INTERFACE EIGEN_USE_BLAS)
else()
  # have to set this on unix even when using static openblas because of a cmake
  # fortran bug
  set(BLA_STATIC FALSE)
  set(BLA_SIZEOF_INTEGER 4)
  target_compile_definitions(sasktran2BuildProperties INTERFACE EIGEN_USE_BLAS)
endif()
set(BLA_VENDOR ${SKTRAN_BLAS_VENDOR})

set(CMAKE_FIND_LIBRARY_PREFIXES "" lib) # openblas include a "lib" prefix in
                                        # it's names
find_package(BLAS)

if(BLAS_FOUND)
  find_package(LAPACK REQUIRED)
  target_link_libraries(sasktran2BuildProperties INTERFACE BLAS::BLAS
                                                           LAPACK::LAPACK)
else()
  if(SKTRAN_BLAS_VENDOR MATCHES "OpenBLAS")
    # Try searching for openblas directly This path should only be taken with
    # the numpy precompiled openblas...
    find_package(OpenBLAS REQUIRED)
    if(NOT WIN32)
      target_link_libraries(sasktran2BuildProperties
                            INTERFACE ${OpenBLAS_LIBRARIES})
    else()
      string(REGEX REPLACE "[.]dll$" ".lib" OpenBLAS_LIBS_TEMP
                           ${OpenBLAS_LIBRARIES})
      string(REGEX REPLACE "/c/" "c:/" OpenBLAS_LIBS_TEMP2
                           ${OpenBLAS_LIBS_TEMP})
      string(REGEX REPLACE "bin" "lib" OpenBLAS_LIBS ${OpenBLAS_LIBS_TEMP2})

      target_link_libraries(sasktran2BuildProperties INTERFACE ${OpenBLAS_LIBS})
    endif()
    target_include_directories(sasktran2BuildProperties
                               INTERFACE ${OpenBLAS_INCLUDE_DIRS})
  endif()
endif()

target_compile_definitions(
  sasktran2BuildProperties
  INTERFACE
    $<$<AND:$<PLATFORM_ID:Windows>,$<STREQUAL:${BLA_VENDOR},OpenBLAS>>:__WIN64__> # for
                                                                                  # compatibility
                                                                                  # with
                                                                                  # Eigen/src/misc/blas.h
)

if(SKTRAN_BLAS_VENDOR MATCHES "Intel")
  find_path(LAPACKE_H_INCLUDE_DIR mkl_lapacke.h REQUIRED) # include directory
                                                          # for <lapacke.h>
  target_include_directories(sasktran2BuildProperties
                             INTERFACE ${LAPACKE_H_INCLUDE_DIR})
  target_compile_definitions(sasktran2BuildProperties INTERFACE SKTRAN_USE_MKL)
else()
  if(WIN32)
    # OpenBlas can give warnings about std::complex linkage on windows/visual
    # studio, but we don't use any complex lapack functions anyway so just
    # disable the warning
    target_compile_options(sasktran2BuildProperties INTERFACE /wd4190)
    target_compile_options(sasktran2BuildProperties INTERFACE /wd4005)
  endif()

  if(SKTRAN_BLAS_VENDOR MATCHES "Generic")
    find_package(LAPACKE CONFIG)
    if(LAPACKE_FOUND)
      target_link_libraries(sasktran2BuildProperties INTERFACE lapacke)
    else()
      target_link_libraries(sasktran2BuildProperties INTERFACE -llapacke)
    endif()
  endif()
endif()

if(WIN32)
  target_compile_options(sasktran2BuildProperties INTERFACE /bigobj)
  target_compile_definitions(sasktran2BuildProperties
                             INTERFACE EIGEN_STRONG_INLINE=inline)
endif()

# -------- END LAPACK FINDING -------

message(STATUS "Lapack Configuration Complete")

add_subdirectory(lib)

option(BUILD_PYTHON "True if Building Python Files" ON)
# We basically just include this for clion help
if(BUILD_PYTHON)
  project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

  set(CMAKE_CXX_STANDARD 17)

  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG REQUIRED)

  python_add_library(
    _core
    MODULE
    src/main.cpp
    src/config.cpp
    src/coordinates.cpp
    src/derivative_mappings.cpp
    src/geometry.cpp
    src/geodetic.cpp
    src/grids.cpp
    src/surface.cpp
    src/atmosphere.cpp
    src/viewing_geometry.cpp
    src/engine.cpp
    src/mie.cpp
    src/output.cpp
    src/math.cpp
    WITH_SOABI)

  target_link_libraries(_core PRIVATE pybind11::headers)
  target_link_libraries(_core PRIVATE sasktran2)

  install(TARGETS _core DESTINATION sasktran2)
else()
  install(
    TARGETS sasktran2BuildProperties
    EXPORT ${EXPORT_NAME}
    LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
    INCLUDES
    DESTINATION ${INCLUDE_INSTALL_DIR})

  write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/sasktran2ConfigVersion.cmake"
    VERSION ${SASKTRAN2_VERSION}
    COMPATIBILITY SameMajorVersion)

  install(
    EXPORT ${EXPORT_NAME}
    FILE sasktran2Targets.cmake
    NAMESPACE sasktran2::
    DESTINATION lib/cmake/sasktran2)

  install(FILES "cmake/sasktran2Config.cmake"
                "${CMAKE_BINARY_DIR}/sasktran2ConfigVersion.cmake"
          DESTINATION lib/cmake/sasktran2)

endif()
