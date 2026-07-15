#include <sasktran2/test_helper.h>

#include <sasktran2.h>

namespace {
    constexpr double earth_radius_m = 6'372'000.0;
    constexpr double tangent_altitude_m = 20'000.0;
    constexpr double observer_altitude_m = 200'000.0;
} // namespace

TEST_CASE("Geometry-relative tangent altitude constructs expected ray",
          "[sasktran2][viewinggeometry][tangentaltitude]") {
    const double horizontal_angle = 0.25;
    const double viewing_azimuth = 0.4;
    sasktran2::Coordinates coordinates(0.6, 1.1, earth_radius_m);
    sasktran2::viewinggeometry::TangentAltitude policy(
        tangent_altitude_m, viewing_azimuth, observer_altitude_m,
        horizontal_angle);

    const auto ray = policy.construct_ray(coordinates);
    const auto tangent_basis =
        coordinates.local_x_y_from_angles(horizontal_angle, 0.0);
    const Eigen::Vector3d expected_look =
        std::cos(viewing_azimuth) * tangent_basis.first +
        std::sin(viewing_azimuth) * tangent_basis.second;
    const Eigen::Vector3d expected_tangent =
        (earth_radius_m + tangent_altitude_m) *
        coordinates.unit_vector_from_angles(horizontal_angle, 0.0);
    const double observer_to_tangent =
        std::sqrt(std::pow(earth_radius_m + observer_altitude_m, 2) -
                  std::pow(earth_radius_m + tangent_altitude_m, 2));

    REQUIRE(ray.look_away.isApprox(expected_look, 1e-13));
    REQUIRE(ray.look_away.norm() == Catch::Approx(1.0));
    REQUIRE(ray.observer.radius() ==
            Catch::Approx(earth_radius_m + observer_altitude_m));
    REQUIRE((ray.observer.position + observer_to_tangent * ray.look_away)
                .isApprox(expected_tangent, 1e-7));
    REQUIRE(ray.relative_azimuth == Catch::Approx(viewing_azimuth));
}

TEST_CASE("Geometry-relative tangent angles follow the Geometry2D plane",
          "[sasktran2][viewinggeometry][tangentaltitude]") {
    sasktran2::Coordinates coordinates(0.2, 0.9, earth_radius_m);
    const double horizontal_angle = -0.3;

    SECTION("zero viewing azimuth follows increasing horizontal angle") {
        sasktran2::viewinggeometry::TangentAltitude policy(
            tangent_altitude_m, 0.0, observer_altitude_m, horizontal_angle);
        const auto ray = policy.construct_ray(coordinates);
        const auto tangent_basis =
            coordinates.local_x_y_from_angles(horizontal_angle, 0.0);

        REQUIRE(ray.look_away.isApprox(tangent_basis.first, 1e-13));
    }

    SECTION("pi viewing azimuth reverses the in-plane direction") {
        sasktran2::viewinggeometry::TangentAltitude policy(
            tangent_altitude_m, EIGEN_PI, observer_altitude_m,
            horizontal_angle);
        const auto ray = policy.construct_ray(coordinates);
        const auto tangent_basis =
            coordinates.local_x_y_from_angles(horizontal_angle, 0.0);

        REQUIRE(ray.look_away.isApprox(-tangent_basis.first, 1e-13));
    }

    SECTION("pi over two follows the invariant direction") {
        sasktran2::viewinggeometry::TangentAltitude policy(
            tangent_altitude_m, EIGEN_PI / 2.0, observer_altitude_m,
            horizontal_angle);
        const auto ray = policy.construct_ray(coordinates);
        const auto tangent_basis =
            coordinates.local_x_y_from_angles(horizontal_angle, 0.0);

        REQUIRE(ray.look_away.isApprox(tangent_basis.second, 1e-13));
        const auto observer_angles = coordinates.angles_from_unit_vector(
            ray.observer.position.normalized());
        REQUIRE(observer_angles.first == Catch::Approx(horizontal_angle));
    }
}

TEST_CASE("Geometry-relative tangent ray is independent of solar angles",
          "[sasktran2][viewinggeometry][tangentaltitude]") {
    sasktran2::Coordinates first_coordinates(0.8, 0.0, earth_radius_m);
    sasktran2::Coordinates second_coordinates(-0.2, 1.4, earth_radius_m);
    sasktran2::viewinggeometry::TangentAltitude policy(
        tangent_altitude_m, 0.35, observer_altitude_m, 0.2);

    const auto first = policy.construct_ray(first_coordinates);
    const auto second = policy.construct_ray(second_coordinates);

    REQUIRE(first.observer.position.isApprox(second.observer.position, 1e-12));
    REQUIRE(first.look_away.isApprox(second.look_away, 1e-12));
}

TEST_CASE("Geometry-relative tangent ray validates inputs",
          "[sasktran2][viewinggeometry][tangentaltitude]") {
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // Preserve the existing C++ policy behavior used to represent rays whose
    // geometric tangent lies below the surface.
    REQUIRE_NOTHROW(sasktran2::viewinggeometry::TangentAltitude(
        -1.0, 0.0, observer_altitude_m));
    REQUIRE_THROWS_AS(sasktran2::viewinggeometry::TangentAltitude(
                          nan, 0.0, observer_altitude_m),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(sasktran2::viewinggeometry::TangentAltitude(
                          tangent_altitude_m, 0.0, tangent_altitude_m - 1.0),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(sasktran2::viewinggeometry::TangentAltitude(
                          tangent_altitude_m, nan, observer_altitude_m),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(sasktran2::viewinggeometry::TangentAltitude(
                          tangent_altitude_m, 0.0, observer_altitude_m, nan),
                      std::invalid_argument);
}

TEST_CASE("Viewing geometry container copies geometry-relative tangent ray",
          "[sasktran2][viewinggeometry][tangentaltitude]") {
    sasktran2::viewinggeometry::TangentAltitude ray(tangent_altitude_m, 0.2,
                                                    observer_altitude_m, -0.1);
    sasktran2::viewinggeometry::ViewingGeometryContainer container;

    container.add_ray(ray);

    REQUIRE(container.observer_rays().size() == 1);
    REQUIRE(dynamic_cast<sasktran2::viewinggeometry::TangentAltitude*>(
                container.observer_rays()[0].get()) != nullptr);
}
