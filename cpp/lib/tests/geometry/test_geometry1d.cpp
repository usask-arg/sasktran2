#include <sasktran2/test_helper.h>

#include <sasktran2.h>

#include <limits>
#include <type_traits>

namespace {
    constexpr double earth_radius_m = 6'372'000.0;

    Eigen::VectorXd altitude_grid(std::initializer_list<double> values) {
        Eigen::VectorXd result(values.size());
        Eigen::Index index = 0;
        for (const double value : values) {
            result[index++] = value;
        }
        return result;
    }

    sasktran2::Geometry1D
    geometry(sasktran2::geometrytype geometry_type,
             sasktran2::grids::interpolation interpolation =
                 sasktran2::grids::interpolation::linear) {
        return sasktran2::Geometry1D(
            0.5, 0.25, earth_radius_m,
            altitude_grid({0.0, 1'000.0, 4'000.0, 10'000.0}), interpolation,
            geometry_type);
    }

    sasktran2::Location spherical_location(double altitude) {
        sasktran2::Location location;
        location.position = Eigen::Vector3d(2.0, -3.0, 4.0).normalized() *
                            (earth_radius_m + altitude);
        return location;
    }

    void require_weight(const std::pair<int, double>& actual, int index,
                        double weight) {
        REQUIRE(actual.first == index);
        REQUIRE(actual.second == Catch::Approx(weight));
    }
} // namespace

TEST_CASE("Geometry1D construction and metadata",
          "[sasktran2][geometry][geometry1d]") {
    static_assert(std::has_virtual_destructor_v<sasktran2::Geometry>);

    auto model_geometry = geometry(sasktran2::geometrytype::spherical);

    REQUIRE(model_geometry.num_atmosphere_dimensions() == 1);
    REQUIRE(model_geometry.size() == 4);
    REQUIRE(model_geometry.altitude_grid().grid().isApprox(
        altitude_grid({0.0, 1'000.0, 4'000.0, 10'000.0})));
    REQUIRE(model_geometry.refractive_index().size() == 4);
    REQUIRE(model_geometry.refractive_index().isOnes());
}

TEST_CASE("Geometry1D accepts explicit coordinates and altitude grid",
          "[sasktran2][geometry][geometry1d]") {
    sasktran2::Coordinates coordinates(Eigen::Vector3d::UnitZ(),
                                       Eigen::Vector3d::UnitX(),
                                       Eigen::Vector3d::UnitZ(), earth_radius_m,
                                       sasktran2::geometrytype::spherical);
    sasktran2::grids::AltitudeGrid grid(altitude_grid({0.0, 2'000.0, 5'000.0}),
                                        sasktran2::grids::gridspacing::variable,
                                        sasktran2::grids::outofbounds::extend,
                                        sasktran2::grids::interpolation::shell);

    sasktran2::Geometry1D model_geometry(std::move(coordinates),
                                         std::move(grid));

    REQUIRE_FALSE(model_geometry.coordinates().sun_forced_z());
    REQUIRE(model_geometry.altitude_grid().interpolation_method() ==
            sasktran2::grids::interpolation::shell);
    REQUIRE(model_geometry.altitude_grid().grid().isApprox(
        altitude_grid({0.0, 2'000.0, 5'000.0})));
}

TEST_CASE("Geometry1D rejects invalid altitude grids",
          "[sasktran2][geometry][geometry1d][validation]") {
    const auto construct = [](Eigen::VectorXd altitudes) {
        return sasktran2::Geometry1D(0.5, 0.0, earth_radius_m,
                                     std::move(altitudes),
                                     sasktran2::grids::interpolation::linear,
                                     sasktran2::geometrytype::spherical);
    };

    SECTION("empty") {
        REQUIRE_THROWS_AS(construct(Eigen::VectorXd(0)), std::runtime_error);
    }
    SECTION("single point") {
        REQUIRE_THROWS_AS(construct(altitude_grid({0.0})), std::runtime_error);
    }
    SECTION("duplicate point") {
        REQUIRE_THROWS_AS(construct(altitude_grid({0.0, 1'000.0, 1'000.0})),
                          std::runtime_error);
    }
    SECTION("decreasing") {
        REQUIRE_THROWS_AS(construct(altitude_grid({0.0, 2'000.0, 1'000.0})),
                          std::runtime_error);
    }
    SECTION("not finite") {
        REQUIRE_THROWS_AS(construct(altitude_grid(
                              {0.0, std::numeric_limits<double>::infinity()})),
                          std::runtime_error);
        REQUIRE_THROWS_AS(construct(altitude_grid(
                              {0.0, std::numeric_limits<double>::quiet_NaN()})),
                          std::runtime_error);
    }
}

TEST_CASE("Geometry1D rejects invalid coordinate inputs",
          "[sasktran2][geometry][geometry1d][validation]") {
    SECTION("solar zenith cosine outside its domain") {
        REQUIRE_THROWS_AS(
            sasktran2::Geometry1D(1.01, 0.0, earth_radius_m,
                                  altitude_grid({0.0, 1'000.0}),
                                  sasktran2::grids::interpolation::linear,
                                  sasktran2::geometrytype::spherical),
            std::runtime_error);
    }
    SECTION("non-finite earth radius") {
        REQUIRE_THROWS_AS(
            sasktran2::Geometry1D(0.5, 0.0,
                                  std::numeric_limits<double>::quiet_NaN(),
                                  altitude_grid({0.0, 1'000.0}),
                                  sasktran2::grids::interpolation::linear,
                                  sasktran2::geometrytype::spherical),
            std::runtime_error);
    }
    SECTION("unknown geometry type") {
        REQUIRE_THROWS_AS(
            sasktran2::Geometry1D(0.5, 0.0, earth_radius_m,
                                  altitude_grid({0.0, 1'000.0}),
                                  sasktran2::grids::interpolation::linear,
                                  static_cast<sasktran2::geometrytype>(99)),
            std::runtime_error);
    }
    SECTION("spherical grid reaches the coordinate origin") {
        REQUIRE_THROWS_AS(
            sasktran2::Geometry1D(0.5, 0.0, earth_radius_m,
                                  altitude_grid({-earth_radius_m, 0.0}),
                                  sasktran2::grids::interpolation::linear,
                                  sasktran2::geometrytype::spherical),
            std::runtime_error);
    }
    SECTION("non-orthogonal manual basis") {
        REQUIRE_THROWS_AS(
            sasktran2::Coordinates(Eigen::Vector3d::UnitZ(),
                                   Eigen::Vector3d::UnitZ(),
                                   Eigen::Vector3d::UnitX(), earth_radius_m,
                                   sasktran2::geometrytype::spherical),
            std::runtime_error);
    }
}

TEST_CASE("Geometry1D interpolates spherical locations",
          "[sasktran2][geometry][geometry1d][interpolation]") {
    auto model_geometry = geometry(sasktran2::geometrytype::spherical);
    std::vector<std::pair<int, double>> weights;

    model_geometry.assign_interpolation_weights(spherical_location(2'500.0),
                                                weights);

    REQUIRE(weights.size() == 2);
    require_weight(weights[0], 1, 0.5);
    require_weight(weights[1], 2, 0.5);

    SECTION("below the grid extends the first point") {
        model_geometry.assign_interpolation_weights(spherical_location(-500.0),
                                                    weights);
        REQUIRE(weights.size() == 1);
        require_weight(weights[0], 0, 1.0);
    }
    SECTION("above the grid extends the last point") {
        model_geometry.assign_interpolation_weights(
            spherical_location(15'000.0), weights);
        REQUIRE(weights.size() == 1);
        require_weight(weights[0], 3, 1.0);
    }
}

TEST_CASE("Geometry1D supports every interpolation mode",
          "[sasktran2][geometry][geometry1d][interpolation]") {
    const auto interpolation = GENERATE(sasktran2::grids::interpolation::linear,
                                        sasktran2::grids::interpolation::shell,
                                        sasktran2::grids::interpolation::lower);
    auto model_geometry =
        geometry(sasktran2::geometrytype::spherical, interpolation);
    std::vector<std::pair<int, double>> weights;

    model_geometry.assign_interpolation_weights(spherical_location(2'500.0),
                                                weights);

    if (interpolation == sasktran2::grids::interpolation::lower) {
        REQUIRE(weights.size() == 1);
        require_weight(weights[0], 1, 1.0);
    } else {
        REQUIRE(weights.size() == 2);
        require_weight(weights[0], 1, 0.5);
        require_weight(weights[1], 2, 0.5);
    }
}

TEST_CASE("Geometry1D uses vertical altitude for plane-parallel geometry",
          "[sasktran2][geometry][geometry1d][interpolation]") {
    auto model_geometry = geometry(sasktran2::geometrytype::planeparallel);
    sasktran2::Location location;
    location.position =
        Eigen::Vector3d(1'000'000.0, -2'000'000.0, earth_radius_m + 2'500.0);
    std::vector<std::pair<int, double>> weights;

    model_geometry.assign_interpolation_weights(location, weights);

    REQUIRE(weights.size() == 2);
    require_weight(weights[0], 1, 0.5);
    require_weight(weights[1], 2, 0.5);
}

TEST_CASE("Geometry1D preserves radial pseudospherical solar interpolation",
          "[sasktran2][geometry][geometry1d][interpolation][regression]") {
    auto model_geometry = geometry(sasktran2::geometrytype::pseudospherical);
    std::vector<std::pair<int, double>> weights;

    // Pseudospherical discrete ordinates uses a spherical shell ray tracer for
    // the solar Chapman-factor path, so non-exact locations are radial.
    model_geometry.assign_interpolation_weights(spherical_location(2'500.0),
                                                weights);
    REQUIRE(weights.size() == 2);
    require_weight(weights[0], 1, 0.5);
    require_weight(weights[1], 2, 0.5);

    // The LOS ray tracer is plane parallel. Its exact grid metadata remains
    // authoritative even though the Cartesian radius is not the altitude.
    sasktran2::Location los_location;
    los_location.position =
        Eigen::Vector3d(1'000'000.0, -2'000'000.0, earth_radius_m + 4'000.0);
    los_location.on_exact_altitude = true;
    los_location.lower_alt_index = 2;

    model_geometry.assign_interpolation_weights(los_location, weights);
    REQUIRE(weights.size() == 1);
    require_weight(weights[0], 2, 1.0);
}

TEST_CASE("Geometry1D validates exact-altitude cache metadata",
          "[sasktran2][geometry][geometry1d][interpolation][regression]") {
    auto model_geometry = geometry(sasktran2::geometrytype::spherical);
    std::vector<std::pair<int, double>> weights;

    SECTION("valid cache is used") {
        auto location = spherical_location(4'000.0);
        location.on_exact_altitude = true;
        location.lower_alt_index = 2;

        model_geometry.assign_interpolation_weights(location, weights);

        REQUIRE(weights.size() == 1);
        require_weight(weights[0], 2, 1.0);
    }
    SECTION("sub-metre shell reconstruction drift uses cache") {
        auto location = spherical_location(4'000.5);
        location.on_exact_altitude = true;
        location.lower_alt_index = 2;

        model_geometry.assign_interpolation_weights(location, weights);

        REQUIRE(weights.size() == 1);
        require_weight(weights[0], 2, 1.0);
    }
    SECTION("stale cache falls back to the location") {
        auto location = spherical_location(2'500.0);
        location.on_exact_altitude = true;
        location.lower_alt_index = 0;

        model_geometry.assign_interpolation_weights(location, weights);

        REQUIRE(weights.size() == 2);
        require_weight(weights[0], 1, 0.5);
        require_weight(weights[1], 2, 0.5);
    }
    SECTION("out-of-range cache falls back to the location") {
        auto location = spherical_location(2'500.0);
        location.on_exact_altitude = true;
        location.lower_alt_index = 100;

        model_geometry.assign_interpolation_weights(location, weights);

        REQUIRE(weights.size() == 2);
        require_weight(weights[0], 1, 0.5);
        require_weight(weights[1], 2, 0.5);
    }
}

TEST_CASE("Geometry1D rejects non-finite interpolation locations",
          "[sasktran2][geometry][geometry1d][interpolation][validation]") {
    auto model_geometry = geometry(sasktran2::geometrytype::spherical);
    sasktran2::Location location;
    location.position = Eigen::Vector3d(
        std::numeric_limits<double>::quiet_NaN(), 0.0, earth_radius_m);
    std::vector<std::pair<int, double>> weights;

    REQUIRE_THROWS_AS(
        model_geometry.assign_interpolation_weights(location, weights),
        std::runtime_error);
}

TEST_CASE("Geometry1D validates its refractive-index profile",
          "[sasktran2][geometry][geometry1d][validation]") {
    auto model_geometry = geometry(sasktran2::geometrytype::spherical);

    SECTION("valid profile") {
        model_geometry.refractive_index() << 1.0003, 1.0002, 1.0001, 1.0;
        REQUIRE_NOTHROW(model_geometry.validate());
    }
    SECTION("wrong size") {
        model_geometry.refractive_index().resize(3);
        REQUIRE_THROWS_AS(model_geometry.validate(), std::runtime_error);
    }
    SECTION("non-positive value") {
        model_geometry.refractive_index()[2] = 0.0;
        REQUIRE_THROWS_AS(model_geometry.validate(), std::runtime_error);
    }
    SECTION("non-finite value") {
        model_geometry.refractive_index()[2] =
            std::numeric_limits<double>::quiet_NaN();
        REQUIRE_THROWS_AS(model_geometry.validate(), std::runtime_error);
    }
}
