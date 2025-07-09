cmake_minimum_required(VERSION 3.20.0)
include(CMakePackageConfigHelpers) # Provides function/macro
                                   # write_basic_package_version_file
include(ExternalProject)
include(FetchContent)

# Convenience
if(UNIX AND NOT APPLE)
  set(LINUX TRUE)
endif()

message(STATUS "CARGO TARGET: ${Rust_CARGO_TARGET}")

# Export all symbols on windows
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

if(UNIX)
  add_definitions(
    "-Wno-inconsistent-missing-override -Wno-suggest-override -Wno-comment -Wno-deprecated-declarations -Wno-psabi"
  )
  # enable_language(Fortran)
endif()

if(WIN32)
  # force MSVC Release => /O2 /DNDEBUG /Zi
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}   /O2 /DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG")
endif()

if(NOT APPLE)
  add_compile_options(-frtti) # Force RTTI on
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

project(
  sasktran2
  VERSION ${SASKTRAN2_VERSION}
  LANGUAGES C CXX)

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
  message(STATUS "Found spdlog")
  target_link_libraries(sasktran2BuildProperties
                        INTERFACE spdlog::spdlog_header_only)

  target_compile_definitions(sasktran2BuildProperties INTERFACE FMT_HEADER_ONLY)
else()

  set(SPDLOG_INSTALL_DIR "${CMAKE_BINARY_DIR}/spdlog-install/")

  ExternalProject_Add(
    spdlog
    URL https://github.com/gabime/spdlog/archive/refs/tags/v1.12.0.tar.gz
    SOURCE_DIR "${CMAKE_BINARY_DIR}/spdlog-src"
    BINARY_DIR "${CMAKE_BINARY_DIR}/spdlog-build"
    INSTALL_DIR "${SPDLOG_INSTALL_DIR}"
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
               -DCMAKE_BUILD_TYPE=Release
               -DSPDLOG_BUILD_EXAMPLES=OFF
               -DSPDLOG_BUILD_SHARED=OFF
               -DSPDLOG_BUILD_STATIC=OFF
               -DSPDLOG_BUILD_BENCH=OFF
               -DSPDLOG_BUILD_TESTS=OFF
               -DSPDLOG_FMT_EXTERNAL=OFF # <- make spdlog self-contained
               -DSPDLOG_INSTALL=ON # <- Necessary for exported CMake targets
  )

  # Ensure the include dir exists before consuming
  file(MAKE_DIRECTORY ${SPDLOG_INSTALL_DIR}/include)

  # Import the spdlog header-only target
  add_library(spdlog_header_only INTERFACE IMPORTED GLOBAL)
  add_dependencies(spdlog_header_only spdlog)

  set_target_properties(
    spdlog_header_only PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                  "${SPDLOG_INSTALL_DIR}/include")

  target_link_libraries(sasktran2BuildProperties INTERFACE spdlog_header_only)
endif()

# add the additional compiler flags
add_compile_options(${CompilerFlags} ${OpenMP_CXX_FLAGS})

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
    message(STATUS "OpenMP found, enabling support")
    target_link_libraries(sasktran2BuildProperties INTERFACE OpenMP::OpenMP_CXX)
    target_compile_definitions(sasktran2BuildProperties
                               INTERFACE SKTRAN_OPENMP_SUPPORT)
    if(NOT LINKED_LIBS)
      set(LINKED_LIBS "${OpenMP_CXX_LIBRARIES}")
    else()
      set(LINKED_LIBS "${LINKED_LIBS};${OpenMP_CXX_LIBRARIES}")
    endif()
  else()
    message(WARNING "OpenMP requested but not found, continuing without it")
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
      /Library/Developer/CommandLineTools/SDKs/MacOSX15.5.sdk/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers
  )
  target_compile_definitions(sasktran2BuildProperties
                             INTERFACE SKTRAN_USE_ACCELERATE)
  target_compile_definitions(sasktran2BuildProperties
                             INTERFACE ACCELERATE_NEW_LAPACK=1)
  # target_compile_definitions(sasktran2BuildProperties INTERFACE
  # ACCELERATE_LAPACK_ILP64=1)
  # target_compile_definitions(sasktran2BuildProperties INTERFACE
  # EIGEN_USE_BLAS)
else()
  # have to set this on unix even when using static openblas because of a cmake
  # fortran bug
  set(BLA_STATIC FALSE)
  set(BLA_SIZEOF_INTEGER 4)
  # target_compile_definitions(sasktran2BuildProperties INTERFACE
  # EIGEN_USE_BLAS)
endif()

if(${SKTRAN_BLAS_VENDOR} MATCHES "SCIPY_OPENBLAS")
  # assume that scipy openblas is installed in the python environment
  find_package(Python REQUIRED COMPONENTS Interpreter)

  execute_process(COMMAND ${Python_EXECUTABLE} -m venv
                          ${CMAKE_CURRENT_BINARY_DIR}/.openblas_venv)

  if(WIN32)
    set(OPENBLAS_PY_BIN
        ${CMAKE_CURRENT_BINARY_DIR}/.openblas_venv/Scripts/python.exe)
  else()
    set(OPENBLAS_PY_BIN ${CMAKE_CURRENT_BINARY_DIR}/.openblas_venv/bin/python)
  endif()

  execute_process(COMMAND ${OPENBLAS_PY_BIN} -m pip install scipy_openblas64)

  execute_process(
    COMMAND ${OPENBLAS_PY_BIN} -c
            "import scipy_openblas64; print(scipy_openblas64.get_include_dir())"
    OUTPUT_VARIABLE OPENBLAS_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(
    COMMAND ${OPENBLAS_PY_BIN} -c
            "import scipy_openblas64; print(scipy_openblas64.get_library())"
    OUTPUT_VARIABLE OPENBLAS_LIBRARY
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(
    COMMAND ${OPENBLAS_PY_BIN} -c
            "import scipy_openblas64; print(scipy_openblas64.get_lib_dir())"
    OUTPUT_VARIABLE OPENBLAS_LIB_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  set(OpenBLAS_DIR ${OPENBLAS_LIB_DIR}/cmake/openblas)

  message(STATUS "Found OpenBLAS: ${OPENBLAS_LIBRARY}")
  message(STATUS "Found OpenBLAS: ${OPENBLAS_INCLUDE_DIR}")
  message(STATUS "Found OpenBLAS: ${OPENBLAS_LIB_DIR}")

  set(BLA_VENDOR OpenBLAS)
  set(SKTRAN_BLAS_VENDOR OpenBLAS)

  target_compile_definitions(sasktran2BuildProperties
                             INTERFACE SKTRAN_NO_LAPACKE)
  target_compile_definitions(sasktran2BuildProperties
                             INTERFACE SKTRAN_SCIPY_BLAS)

  target_compile_definitions(sasktran2BuildProperties
                             INTERFACE BLAS_SYMBOL_PREFIX=scipy_)
  target_compile_definitions(sasktran2BuildProperties
                             INTERFACE BLAS_SYMBOL_SUFFIX=64_)
  target_compile_definitions(sasktran2BuildProperties
                             INTERFACE DHAVE_BLAS_ILP64)
  target_compile_definitions(sasktran2BuildProperties
                             INTERFACE OPENBLAS_ILP64_NAMING_SCHEME)

  target_compile_definitions(sasktran2BuildProperties INTERFACE LAPACK_ILP64)

  target_include_directories(
    sasktran2BuildProperties
    INTERFACE $<BUILD_INTERFACE:${OPENBLAS_INCLUDE_DIR}>
              $<INSTALL_INTERFACE:include>)

  set(OPENBLAS_LIB_FULL
      "${OPENBLAS_LIB_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}${OPENBLAS_LIBRARY}${CMAKE_SHARED_LIBRARY_SUFFIX}"
  )

  target_link_libraries(sasktran2BuildProperties INTERFACE ${OPENBLAS_LIB_FULL})

  if(NOT LINKED_LIBS)
    set(LINKED_LIBS "${OPENBLAS_LIB_FULL}")
  else()
    set(LINKED_LIBS "${LINKED_LIBS};${OPENBLAS_LIB_FULL}")
  endif()

  if(APPLE)
    set_target_properties(
      sasktran2BuildProperties
      PROPERTIES INSTALL_RPATH_USE_LINK_PATH FALSE
                 BUILD_WITH_INSTALL_RPATH TRUE
                 INSTALL_RPATH "" # Or set to a specific value you control
    )

    set_target_properties(sasktran2BuildProperties PROPERTIES MACOSX_RPATH OFF)
  endif()

else()

  set(BLA_VENDOR ${SKTRAN_BLAS_VENDOR})

  set(CMAKE_FIND_LIBRARY_PREFIXES "" lib) # openblas include a "lib" prefix in
                                          # it's names
  find_package(BLAS)

  if(BLAS_FOUND)
    message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
    find_package(LAPACK REQUIRED)
    target_link_libraries(sasktran2BuildProperties INTERFACE BLAS::BLAS
                                                             LAPACK::LAPACK)
    target_compile_definitions(sasktran2BuildProperties
                               INTERFACE SKTRAN_NO_LAPACKE)

    if(NOT LINKED_LIBS)
      set(LINKED_LIBS "${BLAS_LIBRARIES};${LAPACK_LIBRARIES}")
    else()
      set(LINKED_LIBS "${LINKED_LIBS};${BLAS_LIBRARIES};${LAPACK_LIBRARIES}")
    endif()

  else()
    if(SKTRAN_BLAS_VENDOR MATCHES "OpenBLAS")
      # Try searching for openblas directly This path should only be taken with
      # the numpy precompiled openblas...
      find_package(OpenBLAS REQUIRED)
      message(STATUS "Found OpenBLAS: ${OpenBLAS_LIBRARIES}")

      if(NOT WIN32)
        target_link_libraries(sasktran2BuildProperties
                              INTERFACE ${OpenBLAS_LIBRARIES})

        if(NOT LINKED_LIBS)
          set(LINKED_LIBS "${OpenBLAS_LIBRARIES}")
        else()
          set(LINKED_LIBS "${LINKED_LIBS};${OpenBLAS_LIBRARIES}")
        endif()
      else()
        string(REGEX REPLACE "[.]dll$" ".lib" OpenBLAS_LIBS_TEMP
                             ${OpenBLAS_LIBRARIES})
        string(REGEX REPLACE "/c/" "c:/" OpenBLAS_LIBS_TEMP2
                             ${OpenBLAS_LIBS_TEMP})
        string(REGEX REPLACE "bin" "lib" OpenBLAS_LIBS ${OpenBLAS_LIBS_TEMP2})

        target_link_libraries(sasktran2BuildProperties
                              INTERFACE ${OpenBLAS_LIBS})

        if(NOT LINKED_LIBS)
          set(LINKED_LIBS "${OpenBLAS_LIBRARIES}")
        else()
          set(LINKED_LIBS "${LINKED_LIBS};${OpenBLAS_LIBS}")
        endif()
      endif()
      target_include_directories(sasktran2BuildProperties
                                 INTERFACE ${OpenBLAS_INCLUDE_DIRS})

      target_compile_definitions(sasktran2BuildProperties
                                 INTERFACE SKTRAN_NO_LAPACKE)
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
    target_compile_definitions(sasktran2BuildProperties
                               INTERFACE SKTRAN_USE_MKL)
  else()
    if(WIN32)
      # OpenBlas can give warnings about std::complex linkage on windows/visual
      # studio, but we don't use any complex lapack functions anyway so just
      # disable the warning
      target_compile_options(sasktran2BuildProperties INTERFACE /wd4190)
      target_compile_options(sasktran2BuildProperties INTERFACE /wd4005)
    endif()

    if(SKTRAN_BLAS_VENDOR MATCHES "Generic")
      # find_package(LAPACKE CONFIG) if(LAPACKE_FOUND)
      # target_link_libraries(sasktran2BuildProperties INTERFACE lapacke) else()
      # target_link_libraries(sasktran2BuildProperties INTERFACE -llapacke)
      # endif()
    endif()
  endif()

  if(WIN32)
    target_compile_options(sasktran2BuildProperties INTERFACE /bigobj)
    target_compile_definitions(sasktran2BuildProperties
                               INTERFACE EIGEN_STRONG_INLINE=inline)
  endif()
  # -------- END LAPACK FINDING -------
endif()
message(STATUS "Lapack Configuration Complete")

add_subdirectory(lib)

option(INCLUDE_TRACY "Include Tracy Profiling" OFF)

if(INCLUDE_TRACY)
  target_compile_definitions(sasktran2BuildProperties INTERFACE SKTRAN_TRACY)
  option(TRACY_ENABLE "Enable Tracy Profiling" ON)
  option(TRACY_ON_DEMAND "Enable Tracy On Demand" OFF)

  FetchContent_Declare(
    tracy
    GIT_REPOSITORY https://github.com/wolfpld/tracy.git
    GIT_TAG v0.11.1
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE)

  FetchContent_MakeAvailable(tracy)

  target_link_libraries(sasktran2BuildProperties INTERFACE TracyClient)

endif()

add_subdirectory(c_api)

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

if(NOT LINKED_LIBS)
  set(LINKED_LIBS "")
endif()

file(WRITE ${CMAKE_BINARY_DIR}/libs_to_link.txt "${LINKED_LIBS}")
