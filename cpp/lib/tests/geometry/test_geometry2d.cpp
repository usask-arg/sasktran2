#include <sasktran2/test_helper.h>

#include <sasktran2.h>

#include <limits>
#include <numeric>

namespace {
    constexpr double earth_radius_m = 6'372'000.0;

    Eigen::VectorXd grid(std::initializer_list<double> values) {
        Eigen::VectorXd result(values.size());
        Eigen::Index index = 0;
        for (const double value : values) {
            result[index++] = value;
        }
        return result;
    }

    sasktran2::Geometry2D
    geometry2d(sasktran2::grids::interpolation interpolation =
                   sasktran2::grids::interpolation::linear) {
        return sasktran2::Geometry2D(0.5, 0.25, earth_radius_m,
                                     grid({0.0, 10.0, 20.0}),
                                     grid({-0.5, 0.0, 0.5}), interpolation);
    }

    sasktran2::Location location(const sasktran2::Geometry2D& geometry,
                                 double altitude, double horizontal_angle,
                                 double omitted_angle = 0.0) {
        const auto& coordinates = geometry.coordinates();
        const Eigen::Vector3d in_plane =
            std::sin(horizontal_angle) * coordinates.reference_x() +
            std::cos(horizontal_angle) * coordinates.reference_z();
        const Eigen::Vector3d direction =
            std::cos(omitted_angle) * in_plane +
            std::sin(omitted_angle) * coordinates.reference_y();

        sasktran2::Location result;
        result.position = (coordinates.earth_radius() + altitude) * direction;
        return result;
    }

    void require_weight(const std::pair<int, double>& actual, int index,
                        double weight) {
        REQUIRE(actual.first == index);
        REQUIRE(actual.second == Catch::Approx(weight));
    }

    double weight_sum(const std::vector<std::pair<int, double>>& weights) {
        return std::accumulate(weights.begin(), weights.end(), 0.0,
                               [](double total, const auto& weight) {
                                   return total + weight.second;
                               });
    }
} // namespace

TEST_CASE("Geometry2D construction and metadata",
          "[sasktran2][geometry][geometry2d]") {
    auto geometry = geometry2d();

    REQUIRE(geometry.num_atmosphere_dimensions() == 2);
    REQUIRE(geometry.num_altitudes() == 3);
    REQUIRE(geometry.num_horizontal_locations() == 3);
    REQUIRE(geometry.size() == 9);
    REQUIRE(geometry.num_cells() == 4);
    REQUIRE(geometry.location_shape() == std::make_pair(3, 3));
    REQUIRE(geometry.cell_shape() == std::make_pair(2, 2));
    REQUIRE(geometry.altitude_grid().grid().isApprox(grid({0.0, 10.0, 20.0})));
    REQUIRE(geometry.horizontal_angle_grid().isApprox(grid({-0.5, 0.0, 0.5})));
    REQUIRE(geometry.altitude_grid().interpolation_method() ==
            sasktran2::grids::interpolation::linear);
}

TEST_CASE("Geometry2D accepts an explicit rotated coordinate basis",
          "[sasktran2][geometry][geometry2d]") {
    sasktran2::Coordinates coordinates(Eigen::Vector3d::UnitX(),
                                       -Eigen::Vector3d::UnitZ(),
                                       Eigen::Vector3d::UnitX(), earth_radius_m,
                                       sasktran2::geometrytype::spherical);
    sasktran2::Geometry2D geometry(std::move(coordinates), grid({0.0, 10.0}),
                                   grid({-0.5, 0.0, 0.5}),
                                   sasktran2::grids::interpolation::shell);

    REQUIRE(geometry.coordinates().reference_z().isApprox(
        Eigen::Vector3d::UnitX()));
    REQUIRE(geometry.coordinates().reference_x().isApprox(
        -Eigen::Vector3d::UnitZ()));
    REQUIRE(geometry.altitude_grid().interpolation_method() ==
            sasktran2::grids::interpolation::shell);

    const Eigen::Vector3d point = geometry.grid_location(1, 2);
    REQUIRE(point.norm() == Catch::Approx(earth_radius_m + 10.0));
    REQUIRE(
        point.normalized().isApprox(std::sin(0.5) * -Eigen::Vector3d::UnitZ() +
                                    std::cos(0.5) * Eigen::Vector3d::UnitX()));
}

TEST_CASE("Geometry2D validates its grids and coordinate type",
          "[sasktran2][geometry][geometry2d][validation]") {
    const auto construct = [](Eigen::VectorXd altitudes, Eigen::VectorXd angles,
                              sasktran2::grids::interpolation interpolation =
                                  sasktran2::grids::interpolation::linear) {
        return sasktran2::Geometry2D(0.5, 0.0, earth_radius_m,
                                     std::move(altitudes), std::move(angles),
                                     interpolation);
    };

    SECTION("altitude grid") {
        REQUIRE_THROWS_AS(construct(grid({0.0}), grid({-0.5, 0.5})),
                          std::runtime_error);
        REQUIRE_THROWS_AS(construct(grid({0.0, 0.0}), grid({-0.5, 0.5})),
                          std::runtime_error);
        REQUIRE_THROWS_AS(
            construct(grid({0.0, std::numeric_limits<double>::quiet_NaN()}),
                      grid({-0.5, 0.5})),
            std::runtime_error);
        REQUIRE_THROWS_AS(sasktran2::Geometry2D(0.5, 0.0, earth_radius_m,
                                                grid({-earth_radius_m, 0.0}),
                                                grid({-0.5, 0.5})),
                          std::runtime_error);
    }

    SECTION("horizontal grid") {
        REQUIRE_THROWS_AS(construct(grid({0.0, 10.0}), grid({0.0})),
                          std::runtime_error);
        REQUIRE_THROWS_AS(construct(grid({0.0, 10.0}), grid({0.0, 0.0})),
                          std::runtime_error);
        REQUIRE_THROWS_AS(
            construct(grid({0.0, 10.0}),
                      grid({0.0, std::numeric_limits<double>::infinity()})),
            std::runtime_error);
        REQUIRE_THROWS_AS(construct(grid({0.0, 10.0}), grid({0.0, EIGEN_PI})),
                          std::runtime_error);
    }

    SECTION("interpolation mode") {
        REQUIRE_THROWS_AS(
            construct(grid({0.0, 10.0}), grid({-0.5, 0.5}),
                      static_cast<sasktran2::grids::interpolation>(99)),
            std::runtime_error);
    }

    SECTION("only spherical coordinates") {
        sasktran2::Coordinates coordinates(
            0.5, 0.0, earth_radius_m, sasktran2::geometrytype::pseudospherical);
        REQUIRE_THROWS_AS(sasktran2::Geometry2D(std::move(coordinates),
                                                grid({0.0, 10.0}),
                                                grid({-0.5, 0.5})),
                          std::runtime_error);
    }
}

TEST_CASE("Geometry2D uses altitude-fastest location indexing",
          "[sasktran2][geometry][geometry2d][indexing]") {
    auto geometry = geometry2d();

    REQUIRE(geometry.location_index(0, 0) == 0);
    REQUIRE(geometry.location_index(2, 0) == 2);
    REQUIRE(geometry.location_index(0, 1) == 3);
    REQUIRE(geometry.location_index(2, 2) == 8);

    for (int index = 0; index < geometry.size(); ++index) {
        const auto [altitude, horizontal] = geometry.location_indices(index);
        REQUIRE(geometry.location_index(altitude, horizontal) == index);
    }

    REQUIRE_THROWS_AS(geometry.location_index(-1, 0), std::out_of_range);
    REQUIRE_THROWS_AS(geometry.location_index(3, 0), std::out_of_range);
    REQUIRE_THROWS_AS(geometry.location_index(0, 3), std::out_of_range);
    REQUIRE_THROWS_AS(geometry.location_indices(-1), std::out_of_range);
    REQUIRE_THROWS_AS(geometry.location_indices(9), std::out_of_range);
}

TEST_CASE("Geometry2D uses altitude-fastest cell indexing",
          "[sasktran2][geometry][geometry2d][indexing]") {
    auto geometry = geometry2d();

    REQUIRE(geometry.cell_index(0, 0) == 0);
    REQUIRE(geometry.cell_index(1, 0) == 1);
    REQUIRE(geometry.cell_index(0, 1) == 2);
    REQUIRE(geometry.cell_index(1, 1) == 3);

    for (int index = 0; index < geometry.num_cells(); ++index) {
        const auto [altitude, horizontal] = geometry.cell_indices(index);
        REQUIRE(geometry.cell_index(altitude, horizontal) == index);
    }

    REQUIRE_THROWS_AS(geometry.cell_index(-1, 0), std::out_of_range);
    REQUIRE_THROWS_AS(geometry.cell_index(2, 0), std::out_of_range);
    REQUIRE_THROWS_AS(geometry.cell_index(0, 2), std::out_of_range);
    REQUIRE_THROWS_AS(geometry.cell_indices(-1), std::out_of_range);
    REQUIRE_THROWS_AS(geometry.cell_indices(4), std::out_of_range);
}

TEST_CASE("Geometry2D extracts spherical altitude and horizontal angle",
          "[sasktran2][geometry][geometry2d][coordinates]") {
    auto geometry = geometry2d();
    const auto point = location(geometry, 15.0, 0.25, 0.4);

    REQUIRE(geometry.altitude_at(point) == Catch::Approx(15.0).margin(1e-8));
    REQUIRE(geometry.horizontal_angle_at(point) == Catch::Approx(0.25));

    sasktran2::Location on_omitted_axis;
    on_omitted_axis.position =
        (earth_radius_m + 15.0) * geometry.coordinates().reference_y();
    REQUIRE(geometry.horizontal_angle_at(on_omitted_axis) == 0.0);
    REQUIRE(geometry.altitude_at(on_omitted_axis) ==
            Catch::Approx(15.0).margin(1e-8));
}

TEST_CASE("Geometry2D unwraps angles around the grid center",
          "[sasktran2][geometry][geometry2d][coordinates]") {
    sasktran2::Geometry2D geometry(0.5, 0.0, earth_radius_m, grid({0.0, 10.0}),
                                   grid({3.0, 3.1, 3.2}));
    const double raw_angle = -3.1;

    REQUIRE(geometry.horizontal_angle_at(location(geometry, 5.0, raw_angle)) ==
            Catch::Approx(raw_angle + 2.0 * EIGEN_PI));
}

TEST_CASE(
    "Geometry2D cell lookup has finite altitude and extended horizontal cells",
    "[sasktran2][geometry][geometry2d][cells]") {
    auto geometry = geometry2d();

    REQUIRE(geometry.cell_indices(location(geometry, 5.0, -0.25)) ==
            std::make_optional(std::make_pair(0, 0)));
    REQUIRE(geometry.cell_indices(location(geometry, 15.0, 0.25)) ==
            std::make_optional(std::make_pair(1, 1)));

    // Exact interior nodes belong to the upper/right cell.
    REQUIRE(geometry.cell_indices(location(geometry, 10.0, 0.0)) ==
            std::make_optional(std::make_pair(1, 1)));
    // The final nodes remain in the final cell.
    REQUIRE(geometry.cell_indices(location(geometry, 20.0, 0.5)) ==
            std::make_optional(std::make_pair(1, 1)));

    REQUIRE_FALSE(
        geometry.cell_indices(location(geometry, -1.0, 0.0)).has_value());
    REQUIRE_FALSE(
        geometry.cell_indices(location(geometry, 21.0, 0.0)).has_value());

    REQUIRE(geometry.cell_indices(location(geometry, 5.0, -2.0)) ==
            std::make_optional(std::make_pair(0, 0)));
    REQUIRE(geometry.cell_indices(location(geometry, 5.0, 2.0)) ==
            std::make_optional(std::make_pair(0, 1)));
}

TEST_CASE("Geometry2D two-node horizontal grid has one global cell",
          "[sasktran2][geometry][geometry2d][cells]") {
    sasktran2::Geometry2D geometry(0.5, 0.0, earth_radius_m, grid({0.0, 10.0}),
                                   grid({-0.25, 0.25}));

    for (const double angle : {-3.0, -0.25, 0.0, 0.25, 3.0}) {
        REQUIRE(geometry.cell_indices(location(geometry, 5.0, angle)) ==
                std::make_optional(std::make_pair(0, 0)));
    }
}

TEST_CASE("Geometry2D builds bilinear interpolation weights",
          "[sasktran2][geometry][geometry2d][interpolation]") {
    auto geometry = geometry2d();
    std::vector<std::pair<int, double>> weights;

    geometry.assign_interpolation_weights(location(geometry, 5.0, -0.25),
                                          weights);

    REQUIRE(weights.size() == 4);
    require_weight(weights[0], 0, 0.25);
    require_weight(weights[1], 1, 0.25);
    require_weight(weights[2], 3, 0.25);
    require_weight(weights[3], 4, 0.25);
    REQUIRE(weight_sum(weights) == Catch::Approx(1.0));
}

TEST_CASE("Geometry2D collapses interpolation on exact axes",
          "[sasktran2][geometry][geometry2d][interpolation]") {
    auto geometry = geometry2d();
    std::vector<std::pair<int, double>> weights;

    SECTION("exact node") {
        geometry.assign_interpolation_weights(location(geometry, 10.0, 0.0),
                                              weights);
        REQUIRE(weights.size() == 1);
        require_weight(weights[0], 4, 1.0);
    }

    SECTION("exact altitude") {
        geometry.assign_interpolation_weights(location(geometry, 10.0, -0.25),
                                              weights);
        REQUIRE(weights.size() == 2);
        require_weight(weights[0], 1, 0.5);
        require_weight(weights[1], 4, 0.5);
    }

    SECTION("exact horizontal angle") {
        geometry.assign_interpolation_weights(location(geometry, 5.0, 0.0),
                                              weights);
        REQUIRE(weights.size() == 2);
        require_weight(weights[0], 3, 0.5);
        require_weight(weights[1], 4, 0.5);
    }
}

TEST_CASE("Geometry2D cell interpolation coordinates are one-sided",
          "[sasktran2][geometry2d]") {
    auto geometry = geometry2d();
    const auto corner = location(geometry, 10.0, 0.0);

    const auto lower_left =
        geometry.cell_interpolation_coordinates(corner, 0, 0);
    REQUIRE(lower_left.first == Catch::Approx(1.0));
    REQUIRE(lower_left.second == Catch::Approx(1.0));

    const auto upper_right =
        geometry.cell_interpolation_coordinates(corner, 1, 1);
    REQUIRE(upper_right.first == Catch::Approx(0.0));
    REQUIRE(upper_right.second == Catch::Approx(0.0));

    const auto left_extension = geometry.cell_interpolation_coordinates(
        location(geometry, 5.0, -2.0), 0, 0);
    REQUIRE(left_extension.second == Catch::Approx(0.0));
    const auto right_extension = geometry.cell_interpolation_coordinates(
        location(geometry, 5.0, 2.0), 0, 1);
    REQUIRE(right_extension.second == Catch::Approx(1.0));
    REQUIRE_THROWS_AS(geometry.cell_interpolation_coordinates(corner, -1, 0),
                      std::out_of_range);
    REQUIRE_THROWS_AS(geometry.cell_interpolation_coordinates(corner, 0, 2),
                      std::out_of_range);
}

TEST_CASE("Geometry2D combines each altitude interpolation mode with linear "
          "horizontal interpolation",
          "[sasktran2][geometry][geometry2d][interpolation]") {
    const auto interpolation = GENERATE(sasktran2::grids::interpolation::linear,
                                        sasktran2::grids::interpolation::shell,
                                        sasktran2::grids::interpolation::lower);
    auto geometry = geometry2d(interpolation);
    std::vector<std::pair<int, double>> weights;

    geometry.assign_interpolation_weights(location(geometry, 2.5, -0.25),
                                          weights);

    if (interpolation == sasktran2::grids::interpolation::lower) {
        REQUIRE(weights.size() == 2);
        require_weight(weights[0], 0, 0.5);
        require_weight(weights[1], 3, 0.5);
    } else if (interpolation == sasktran2::grids::interpolation::shell) {
        REQUIRE(weights.size() == 4);
        for (const auto& weight : weights) {
            REQUIRE(weight.second == Catch::Approx(0.25));
        }
    } else {
        REQUIRE(weights.size() == 4);
        require_weight(weights[0], 0, 0.375);
        require_weight(weights[1], 1, 0.125);
        require_weight(weights[2], 3, 0.375);
        require_weight(weights[3], 4, 0.125);
    }
    REQUIRE(weight_sum(weights) == Catch::Approx(1.0));
}

TEST_CASE("Geometry2D interpolation clamps both axes",
          "[sasktran2][geometry][geometry2d][interpolation]") {
    auto geometry = geometry2d();
    std::vector<std::pair<int, double>> weights;

    SECTION("below altitude") {
        geometry.assign_interpolation_weights(location(geometry, -100.0, -0.25),
                                              weights);
        REQUIRE(weights.size() == 2);
        require_weight(weights[0], 0, 0.5);
        require_weight(weights[1], 3, 0.5);
    }
    SECTION("above altitude") {
        geometry.assign_interpolation_weights(location(geometry, 100.0, -0.25),
                                              weights);
        REQUIRE(weights.size() == 2);
        require_weight(weights[0], 2, 0.5);
        require_weight(weights[1], 5, 0.5);
    }
    SECTION("below horizontal") {
        geometry.assign_interpolation_weights(location(geometry, 5.0, -2.0),
                                              weights);
        REQUIRE(weights.size() == 2);
        require_weight(weights[0], 0, 0.5);
        require_weight(weights[1], 1, 0.5);
    }
    SECTION("above horizontal") {
        geometry.assign_interpolation_weights(location(geometry, 5.0, 2.0),
                                              weights);
        REQUIRE(weights.size() == 2);
        require_weight(weights[0], 6, 0.5);
        require_weight(weights[1], 7, 0.5);
    }
    SECTION("outside both axes") {
        geometry.assign_interpolation_weights(location(geometry, -100.0, 2.0),
                                              weights);
        REQUIRE(weights.size() == 1);
        require_weight(weights[0], 6, 1.0);
    }
}

TEST_CASE("Geometry2D exactly reproduces a bilinear nodal field",
          "[sasktran2][geometry][geometry2d][interpolation]") {
    auto geometry = geometry2d();
    Eigen::VectorXd field(geometry.size());
    const auto value = [](double altitude, double angle) {
        return 1.0 + 2.0 * altitude + 3.0 * angle + 4.0 * altitude * angle;
    };

    for (int horizontal = 0; horizontal < geometry.num_horizontal_locations();
         ++horizontal) {
        for (int altitude = 0; altitude < geometry.num_altitudes();
             ++altitude) {
            field[geometry.location_index(altitude, horizontal)] =
                value(geometry.altitude_grid().grid()[altitude],
                      geometry.horizontal_angle_grid()[horizontal]);
        }
    }

    for (const auto& query :
         {std::make_pair(2.5, -0.25), std::make_pair(7.5, 0.25),
          std::make_pair(15.0, 0.1)}) {
        std::vector<std::pair<int, double>> weights;
        geometry.assign_interpolation_weights(
            location(geometry, query.first, query.second, 0.3), weights);

        double interpolated = 0.0;
        for (const auto& [index, weight] : weights) {
            interpolated += weight * field[index];
        }
        REQUIRE(interpolated ==
                Catch::Approx(value(query.first, query.second)));
    }
}

TEST_CASE("Geometry2D interpolation is polymorphic through Geometry",
          "[sasktran2][geometry][geometry2d][interpolation]") {
    auto geometry = geometry2d();
    const sasktran2::Geometry& base = geometry;
    std::vector<std::pair<int, double>> weights;

    base.assign_interpolation_weights(location(geometry, 5.0, -0.25), weights);

    REQUIRE(base.num_atmosphere_dimensions() == 2);
    REQUIRE(base.size() == 9);
    REQUIRE(weights.size() == 4);
    REQUIRE(weight_sum(weights) == Catch::Approx(1.0));
}

TEST_CASE("Geometry2D has a deterministic antipodal seam",
          "[sasktran2][geometry][geometry2d][coordinates][cells]") {
    auto geometry = geometry2d();
    const auto positive_side = location(geometry, 5.0, EIGEN_PI - 1e-6);
    const auto negative_side = location(geometry, 5.0, -EIGEN_PI + 1e-6);
    std::vector<std::pair<int, double>> weights;

    REQUIRE(geometry.cell_indices(positive_side) ==
            std::make_optional(std::make_pair(0, 1)));
    REQUIRE(geometry.cell_indices(negative_side) ==
            std::make_optional(std::make_pair(0, 0)));

    geometry.assign_interpolation_weights(positive_side, weights);
    REQUIRE(weights.size() == 2);
    require_weight(weights[0], 6, 0.5);
    require_weight(weights[1], 7, 0.5);

    geometry.assign_interpolation_weights(negative_side, weights);
    REQUIRE(weights.size() == 2);
    require_weight(weights[0], 0, 0.5);
    require_weight(weights[1], 1, 0.5);
}

TEST_CASE("Geometry2D clamps boundary roundoff",
          "[sasktran2][geometry][geometry2d][interpolation][cells]") {
    auto geometry = geometry2d();
    std::vector<std::pair<int, double>> weights;
    const auto point = location(geometry, -5e-9, -0.5 - 5e-9);

    geometry.assign_interpolation_weights(point, weights);

    REQUIRE(weights.size() == 1);
    require_weight(weights[0], 0, 1.0);
    REQUIRE(geometry.cell_indices(point) ==
            std::make_optional(std::make_pair(0, 0)));
}

TEST_CASE("Geometry2D rejects non-finite locations",
          "[sasktran2][geometry][geometry2d][validation]") {
    auto geometry = geometry2d();
    sasktran2::Location invalid;
    invalid.position = Eigen::Vector3d(std::numeric_limits<double>::quiet_NaN(),
                                       0.0, earth_radius_m);
    std::vector<std::pair<int, double>> weights;

    REQUIRE_THROWS_AS(geometry.altitude_at(invalid), std::runtime_error);
    REQUIRE_THROWS_AS(geometry.horizontal_angle_at(invalid),
                      std::runtime_error);
    REQUIRE_THROWS_AS(geometry.cell_indices(invalid), std::runtime_error);
    REQUIRE_THROWS_AS(geometry.assign_interpolation_weights(invalid, weights),
                      std::runtime_error);
}
