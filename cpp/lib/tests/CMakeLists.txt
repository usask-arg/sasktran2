add_library(
  sasktran2_tests
  SHARED
  test_main.cpp
  dual/test_dual.cpp
  viewinggeometry/test_tangentaltitude.cpp
  viewinggeometry/test_groundviewing.cpp
  unitsphere/test_lebedev.cpp
  geometry/test_geometry.cpp
  grids/test_grid.cpp
  raytracing/test_spherical_shell.cpp
  raytracing/test_plane_parallel.cpp
  engine/test_engine_basic.cpp
  engine/wf/test_wf_singlescatter_scalar.cpp
  engine/wf/test_wf_vector.cpp
  engine/hr/test_hr_basic.cpp
  sktran_disco/test_util.cpp
  sktran_disco/tests_lowlevel.cpp
  sktran_disco/legacy/test_scalar.cpp
  sktran_disco/legacy/test_vector.cpp
  sktran_disco/test_crash.cpp
  twostream/test_twostream.cpp
  twostream/test_backprop.cpp
  sktran_disco/benchmark/bench_lowlevel.cpp
  engine/benchmark/bench_singlescatter.cpp
  engine/benchmark/bench_do_wf.cpp
  engine/benchmark/bench_hr.cpp
  engine/benchmark/bench_occ.cpp
  rust/test_integration.cpp
  atmosphere/surface.cpp)

# Some special logic to support both CATCH2 version 2 and 3 for now
find_package(Catch2 2 CONFIG QUIET)

target_link_libraries(
  sasktran2_tests
  PRIVATE sasktran2 $<$<PLATFORM_ID:Windows>:bcrypt>
          $<$<PLATFORM_ID:Linux,Darwin>:dl> $<$<PLATFORM_ID:Linux>:gfortran>
          $<$<PLATFORM_ID:Linux>:rt>)
if(Catch2_FOUND)
  # Found Catch2 Version 2
  # Error, no longer support catch2
  message(FATAL_ERROR "Catch2 Version 2 is no longer supported, please upgrade to Catch2 Version 3")
else()
  # Look for Catch2 Version 3
  message("Could not find Catch2 Version 2, Looking for Version 3")
  find_package(Catch2 3 CONFIG REQUIRED)
  target_link_libraries(sasktran2_tests PRIVATE Catch2::Catch2)
  target_compile_definitions(sasktran2_tests PUBLIC SKTRAN_CATCH2_VERSION3)
endif()

# Ensure all test object files are linked when consumed as a static lib
# macOS needs -force_load, Linux needs --whole-archive/--no-whole-archive
if(APPLE)
  target_link_options(
    sasktran2_tests
    INTERFACE "-Wl,-force_load,$<TARGET_FILE:sasktran2_tests>")
else()
  target_link_options(
    sasktran2_tests
    INTERFACE "-Wl,--whole-archive" "-lsasktran2_tests" "-Wl,--no-whole-archive")
endif()

if(INCLUDE_TRACY)
  target_link_libraries(sasktran2_tests INTERFACE TracyClient)
endif()

install(
  TARGETS sasktran2_tests
  EXPORT ${EXPORT_NAME}
  LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
  INCLUDES
  DESTINATION ${INCLUDE_INSTALL_DIR})
