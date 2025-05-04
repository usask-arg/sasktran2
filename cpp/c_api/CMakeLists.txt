# Set the include directories
target_include_directories(
  sasktran2BuildProperties
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include/c_api/>
            $<INSTALL_INTERFACE:include>)

# Add sub libraries
add_library(
  csasktran2
  geometry.cpp
  viewing_geometry.cpp
  atmosphere.cpp
  config.cpp
  engine.cpp
  output.cpp
  geodetic.cpp
  brdf.cpp
  deriv_mapping.cpp)

target_link_libraries(csasktran2 PUBLIC sasktran2)

set_target_properties(csasktran2 PROPERTIES OUTPUT_NAME csasktran2)

install(
  TARGETS csasktran2
  EXPORT ${EXPORT_NAME}
  LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
  INCLUDES
  DESTINATION ${INCLUDE_INSTALL_DIR})
