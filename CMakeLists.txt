cmake_minimum_required(VERSION 3.20.0)
include(CMakePackageConfigHelpers )# Provides function/macro write_basic_package_version_file

# Convenience
if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

# Export all symbols on windows
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

if(UNIX)
    add_definitions("-Wno-inconsistent-missing-override -Wno-suggest-override -Wno-comment -Wno-deprecated-declarations -Wno-psabi")
    enable_language(Fortran)
endif()

execute_process(COMMAND git describe --tags
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_DESCRIBE_TEMP
        )

if(EXISTS ${GIT_DESCRIBE_TEMP})
    string(REGEX MATCH "v([0-9\\.]+)" SASKTRAN2_VERSION ${GIT_DESCRIBE_TEMP})
    string(REPLACE "v" "" SASKTRAN2_VERSION ${GIT_DESCRIBE_TEMP})
else()
    set(SASKTRAN2_VERSION 0.0.0)
endif()

project (sasktran2 VERSION ${SASKTRAN2_VERSION})

set(CMAKE_CXX_STANDARD 17)

if(APPLE)
    # This is nominally a windows definition but for some reason BOOST < 1.8.3 has problems on clang without this?
    add_compile_definitions(_HAS_AUTO_PTR_ETC=0)
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(FORCE_BOOST_STATIC "Force the use of boost static libraries" OFF)
# Try to find boost 1.6 or greater
set(Boost_USE_STATIC_LIBS ${FORCE_BOOST_STATIC})
find_package(Boost REQUIRED COMPONENTS log)

# Try to find OpenMP
find_package(OpenMP REQUIRED)

find_package(Eigen3 REQUIRED)

# add the additional compiler flags
set (CMAKE_CXX_FLAGS " ${CompilerFlags} ${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")


set(CMAKE_POSITION_INDEPENDENT_CODE ON)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Directory to install the static library files to
set (STATIC_LIB_INSTALL_DIR lib)
set (INCLUDE_INSTALL_DIR include)
set (EXPORT_NAME sasktran2Targets)


# Add an interface target for storing build properties common to all sasktran modules
add_library(sasktran2BuildProperties INTERFACE)

target_link_libraries(sasktran2BuildProperties INTERFACE OpenMP::OpenMP_CXX)
target_link_libraries(sasktran2BuildProperties INTERFACE Eigen3::Eigen)
target_link_libraries(sasktran2BuildProperties INTERFACE Boost::log)

if(WIN32)
    target_compile_definitions(sasktran2BuildProperties INTERFACE BOOST_USE_WINAPI_VERSION=BOOST_WINAPI_VERSION_WIN7)
endif()

# -------- LAPACK FINDING -----------
# We support multiple blas vendors, but these can require some funny options
set(SKTRAN_BLAS_VENDOR "OpenBLAS" CACHE STRING
        "Sets the BLAS/LAPACK vendor. See https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors.")
set_property(CACHE SKTRAN_BLAS_VENDOR PROPERTY STRINGS OpenBLAS Intel10_64lp Apple Generic)

if(SKTRAN_BLAS_VENDOR MATCHES "Intel")
    set(BLA_STATIC FALSE)
    set(BLA_SIZEOF_INTEGER 4)
elseif(SKTRAN_BLAS_VENDOR MATCHES "Apple")
    set(BLA_STATIC TRUE)
    set(BLA_SIZEOF_INTEGER 4)
    target_include_directories(sasktran2BuildProperties INTERFACE /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers)
    target_compile_definitions(sasktran2BuildProperties INTERFACE SKTRAN_USE_ACCELERATE)
    target_compile_definitions(sasktran2BuildProperties INTERFACE ACCELERATE_NEW_LAPACK=1)
    target_compile_definitions(sasktran2BuildProperties INTERFACE EIGEN_USE_BLAS)
else()
    # have to set this on unix even when using static openblas because of a cmake fortran bug
    set(BLA_STATIC FALSE)
    set(BLA_SIZEOF_INTEGER 4)
    target_compile_definitions(sasktran2BuildProperties INTERFACE EIGEN_USE_BLAS)
endif()
#if(NOT (SKTRAN_BLAS_VENDOR MATCHES "Generic"))
    set(BLA_VENDOR ${SKTRAN_BLAS_VENDOR})
#endif()
set(CMAKE_FIND_LIBRARY_PREFIXES "" lib)  # openblas include a "lib" prefix in it's names
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(sasktran2BuildProperties INTERFACE BLAS::BLAS LAPACK::LAPACK)
target_compile_definitions(sasktran2BuildProperties INTERFACE
        $<$<AND:$<PLATFORM_ID:Windows>,$<STREQUAL:${BLA_VENDOR},OpenBLAS>>:__WIN64__>   # for compatibility with Eigen/src/misc/blas.h
        )

if(SKTRAN_BLAS_VENDOR MATCHES "Intel")
    find_path(LAPACKE_H_INCLUDE_DIR mkl_lapacke.h REQUIRED)  # include directory for <lapacke.h>
    target_include_directories(sasktran2BuildProperties INTERFACE ${LAPACKE_H_INCLUDE_DIR})
    target_compile_definitions(sasktran2BuildProperties INTERFACE SKTRAN_USE_MKL)
else()
    if(WIN32)
        # OpenBlas can give warnings about std::complex linkage on windows/visual studio, but we don't use any complex lapack functions anyway so just disable the warning
        target_compile_options(sasktran2BuildProperties INTERFACE /wd4190)
    endif()

    if(SKTRAN_BLAS_VENDOR MATCHES "Generic")
        find_package(LAPACKE REQUIRED CONFIG)
        target_link_libraries(sasktran2BuildProperties INTERFACE lapacke)
    endif()
endif()

# -------- END LAPACK FINDING -------

message(STATUS "Lapack Configuration Complete")



add_subdirectory (src)

install (TARGETS sasktran2BuildProperties
        EXPORT ${EXPORT_NAME}
        LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
        INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR})


if(PROJECT_IS_TOP_LEVEL)
    write_basic_package_version_file( "${CMAKE_BINARY_DIR}/sasktran2ConfigVersion.cmake"
            VERSION ${SASKTRAN2_VERSION}
            COMPATIBILITY SameMajorVersion
            )

    install(EXPORT ${EXPORT_NAME}
            FILE sasktran2Targets.cmake
            NAMESPACE sasktran2::
            DESTINATION lib/cmake/sasktran2
            )

    install( FILES "cmake/sasktran2Config.cmake" "${CMAKE_BINARY_DIR}/sasktran2ConfigVersion.cmake"
            DESTINATION lib/cmake/sasktran2)

    # We basically just include this for clion help
    if(BUILD_PYTHON)
        add_subdirectory(python)
    endif()
endif()