#include "sasktran2/viewinggeometry.h"
#include <sasktran2/test_helper.h>

#include <sasktran2.h>

TEST_CASE("Spherical Shell Raytracer - Observer Outside limb Viewing",
          "[sasktran2][raytracing]") {
    // Construct a constant spacing grid

    double dx = 1000;
    int nx = 100;
    double x0 = 0;
    double tangent_altitude = GENERATE(10000.0, 20000.0, 30405.0);
    bool include_refraction = GENERATE(false, true);

    Eigen::VectorXd grid_values(nx);
    for (int i = 0; i < nx; ++i) {
        grid_values(i) = x0 + i * dx;
    }

    sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(
        std::move(grid_values), sasktran2::grids::gridspacing::constant,
        sasktran2::grids::outofbounds::extend,
        sasktran2::grids::interpolation::linear);

    sasktran2::Coordinates coord(1, 0, 6372000);

    sasktran2::Geometry1D geo(std::move(coord), std::move(grid));

    sasktran2::raytracing::SphericalShellRayTracer raytracer(geo);

    sasktran2::viewinggeometry::TangentAltitude ray_policy(tangent_altitude, 0,
                                                           600000);

    auto viewing_ray = ray_policy.construct_ray(geo.coordinates());

    sasktran2::raytracing::TracedRay traced_ray;
    raytracer.trace_ray(viewing_ray, traced_ray, include_refraction);

    double min_alt = 1e99;
    for (const auto& layer : traced_ray.layers) {
        double entrance_altitude =
            layer.entrance.radius() - geo.coordinates().earth_radius();
        double exit_altitude =
            layer.exit.radius() - geo.coordinates().earth_radius();
        if (entrance_altitude < min_alt) {
            min_alt = entrance_altitude;
        }
    }

    REQUIRE(fabs(tangent_altitude - min_alt) < 1e-8);
}

TEST_CASE(
    "Spherical Shell Raytracer - Observer Outside limb Viewing - Refraction",
    "[sasktran2][raytracing]") {
    // Construct a constant spacing grid

    double dx = 1000;
    int nx = 66;
    double x0 = 0;
    double tangent_altitude = GENERATE(5000, 10000.0, 18000, 20000.0, 30405.0);

    Eigen::VectorXd grid_values(nx);
    for (int i = 0; i < nx; ++i) {
        grid_values(i) = x0 + i * dx;
    }

    sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(
        std::move(grid_values), sasktran2::grids::gridspacing::constant,
        sasktran2::grids::outofbounds::extend,
        sasktran2::grids::interpolation::linear);

    sasktran2::Coordinates coord(1, 0, 6372000);

    sasktran2::Geometry1D geo(std::move(coord), std::move(grid));

    geo.refractive_index() << 1.00027691, 1.00025137, 1.00022759, 1.00020559,
        1.00018526, 1.00016652, 1.00014926, 1.00013342, 1.00011888, 1.0001056,
        1.0000935, 1.00008041, 1.00006916, 1.00005949, 1.00005117, 1.00004402,
        1.00003763, 1.00003217, 1.0000275, 1.00002351, 1.0000201, 1.00001713,
        1.00001461, 1.00001246, 1.00001062, 1.00000906, 1.00000775, 1.00000664,
        1.00000568, 1.00000486, 1.00000416, 1.00000357, 1.00000306, 1.00000263,
        1.00000226, 1.00000194, 1.00000166, 1.00000143, 1.00000122, 1.00000105,
        1.0000009, 1.00000079, 1.00000069, 1.0000006, 1.00000052, 1.00000046,
        1.0000004, 1.00000035, 1.0000003, 1.00000027, 1.00000023, 1.00000021,
        1.00000018, 1.00000016, 1.00000014, 1.00000013, 1.00000011, 1.0000001,
        1.00000009, 1.00000008, 1.00000007, 1.00000006, 1.00000005, 1.00000005,
        1.00000004, 1.00000004;

    sasktran2::raytracing::SphericalShellRayTracer raytracer(geo);

    sasktran2::viewinggeometry::TangentAltitude ray_policy(tangent_altitude, 0,
                                                           600000);

    auto viewing_ray = ray_policy.construct_ray(geo.coordinates());

    sasktran2::raytracing::TracedRay traced_ray;
    raytracer.trace_ray(viewing_ray, traced_ray, true);

    double min_alt = 1e99;
    for (const auto& layer : traced_ray.layers) {
        double entrance_altitude =
            layer.entrance.radius() - geo.coordinates().earth_radius();
        double exit_altitude =
            layer.exit.radius() - geo.coordinates().earth_radius();
        if (entrance_altitude < min_alt) {
            min_alt = entrance_altitude;
        }
    }

    REQUIRE(min_alt < tangent_altitude);
}

TEST_CASE("Spherical Shell Raytracer - Observer Inside Limb Viewing",
          "[sasktran2][raytracing]") {
    // Construct a constant spacing grid

    double dx = 1000;
    int nx = 100;
    double x0 = 0;
    double tangent_altitude = GENERATE(10000.0, 20000.0, 30405.0);
    double dx_factor = GENERATE(0, 0.01, 0.2, 2, 2.2, 5.7);
    bool include_refraction = GENERATE(true, false);

    Eigen::VectorXd grid_values(nx);
    for (int i = 0; i < nx; ++i) {
        grid_values(i) = x0 + i * dx;
    }

    double entrance_altitude = grid_values(Eigen::last);
    double exit_altitude = grid_values(Eigen::last - 1);

    sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(
        std::move(grid_values), sasktran2::grids::gridspacing::constant,
        sasktran2::grids::outofbounds::extend,
        sasktran2::grids::interpolation::linear);

    sasktran2::Coordinates coord(1, 0, 6372000);

    sasktran2::Geometry1D geo(std::move(coord), std::move(grid));

    sasktran2::raytracing::SphericalShellRayTracer raytracer(geo);

    sasktran2::viewinggeometry::TangentAltitude ray_policy(
        tangent_altitude, 0, tangent_altitude + dx * dx_factor);

    auto viewing_ray = ray_policy.construct_ray(geo.coordinates());

    sasktran2::raytracing::TracedRay traced_ray;
    raytracer.trace_ray(viewing_ray, traced_ray, include_refraction);

    double min_alt = 1e99;
    for (const auto& layer : traced_ray.layers) {
        double diff = layer.exit.radius() - geo.coordinates().earth_radius() -
                      entrance_altitude;
        REQUIRE(abs(diff) < 1e-8);
        entrance_altitude =
            layer.entrance.radius() - geo.coordinates().earth_radius();
        exit_altitude = layer.exit.radius() - geo.coordinates().earth_radius();
        if (entrance_altitude < min_alt) {
            min_alt = entrance_altitude;
        }
    }

    REQUIRE(fabs(tangent_altitude - min_alt) < 1e-8);
}

TEST_CASE("Spherical Shell Raytracer - Observer Outside Ground Viewing",
          "[sasktran2][raytracing]") {
    // Construct a constant spacing grid

    double dx = 1000;
    int nx = 100;
    double x0 = 0;
    double cos_viewing_zenith = GENERATE(1, 0.5, 0.1, 0.01, 0.001);
    double dx_factor = GENERATE(0, 0.01, 0.2, 2, 2.2, 5.7);
    bool include_refraction = GENERATE(true, false);

    Eigen::VectorXd grid_values(nx);
    for (int i = 0; i < nx; ++i) {
        grid_values(i) = x0 + i * dx;
    }

    double entrance_altitude = grid_values(0);
    double exit_altitude = grid_values(0);

    sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(
        std::move(grid_values), sasktran2::grids::gridspacing::constant,
        sasktran2::grids::outofbounds::extend,
        sasktran2::grids::interpolation::linear);

    sasktran2::Coordinates coord(1, 0, 6372000);

    sasktran2::Geometry1D geo(std::move(coord), std::move(grid));

    sasktran2::raytracing::SphericalShellRayTracer raytracer(geo);

    sasktran2::viewinggeometry::GroundViewingSolar ray_policy(
        0.6, 0, cos_viewing_zenith, 200000);

    auto viewing_ray = ray_policy.construct_ray(geo.coordinates());

    sasktran2::raytracing::TracedRay traced_ray;
    raytracer.trace_ray(viewing_ray, traced_ray, include_refraction);

    double min_alt = 1e99;
    for (const auto& layer : traced_ray.layers) {
        double diff = layer.exit.radius() - geo.coordinates().earth_radius() -
                      entrance_altitude;
        REQUIRE(abs(diff) < 1e-8);
        entrance_altitude =
            layer.entrance.radius() - geo.coordinates().earth_radius();
        exit_altitude = layer.exit.radius() - geo.coordinates().earth_radius();
        if (exit_altitude < min_alt) {
            min_alt = exit_altitude;
        }
    }

    REQUIRE(fabs(min_alt) < 1e-8);
}
