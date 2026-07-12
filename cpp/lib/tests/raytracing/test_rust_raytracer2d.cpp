#include <sasktran2/test_helper.h>

#include <sasktran2.h>

#ifdef SKTRAN_RUST_SUPPORT

namespace {
    constexpr double earth_radius = 10.0;
    constexpr double tolerance = 1e-9;

    Eigen::VectorXd values(std::initializer_list<double> input) {
        Eigen::VectorXd result(input.size());
        int index = 0;
        for (double value : input) {
            result[index++] = value;
        }
        return result;
    }

    sasktran2::Geometry2D geometry(
        sasktran2::grids::interpolation interpolation =
            sasktran2::grids::interpolation::linear) {
        return sasktran2::Geometry2D(
            0.6, 0.3, earth_radius, values({0.0, 10.0, 20.0}),
            values({-0.5, 0.0, 0.5}), interpolation);
    }

    Eigen::Vector3d radial_direction(double angle) {
        return {std::sin(angle), 0.0, std::cos(angle)};
    }

    sasktran2::viewinggeometry::ViewingRay
    ray(const Eigen::Vector3d& observer, const Eigen::Vector3d& look) {
        sasktran2::viewinggeometry::ViewingRay result;
        result.observer.position = observer;
        result.look_away = look.normalized();
        result.relative_azimuth = 0.0;
        return result;
    }

    void require_close(double actual, double expected,
                       double absolute_tolerance = tolerance) {
        CAPTURE(actual, expected);
        REQUIRE(std::abs(actual - expected) <= absolute_tolerance);
    }

    void require_weights_match_geometry(
        const sasktran2::Geometry2D& geometry,
        const sasktran2::Location& location) {
        std::vector<std::pair<int, double>> expected;
        geometry.assign_interpolation_weights(location, expected);
        REQUIRE(location.interpolation_weights.size() == expected.size());
        for (std::size_t index = 0; index < expected.size(); ++index) {
            REQUIRE(location.interpolation_weights[index].first ==
                    expected[index].first);
            require_close(location.interpolation_weights[index].second,
                          expected[index].second);
        }
    }
} // namespace

TEST_CASE("RustRayTracer2D traces radial layers and exposes 2D cells",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay2D traced;
    const double angle = 0.25;
    tracer.trace_ray(ray(40.0 * radial_direction(angle),
                         -radial_direction(angle)),
                     traced);

    REQUIRE(traced.ground_is_hit);
    REQUIRE(traced.is_straight);
    REQUIRE(traced.layers.size() == 2);
    REQUIRE(traced.layers[0].altitude_cell == 0);
    REQUIRE(traced.layers[0].horizontal_cell == 1);
    REQUIRE(traced.layers[1].altitude_cell == 1);
    REQUIRE(traced.layers[1].horizontal_cell == 1);
    require_close(traced.layers[0].layer_distance, 10.0);
    require_close(traced.layers[1].layer_distance, 10.0);

    for (const auto& layer : traced.layers) {
        require_weights_match_geometry(geo, layer.entrance);
        require_weights_match_geometry(geo, layer.exit);
        REQUIRE(layer.entrance.interpolation_weights.size() <= 4);
        REQUIRE(layer.exit.interpolation_weights.size() <= 4);
        REQUIRE(std::isfinite(layer.cos_sza_entrance));
        REQUIRE(std::isfinite(layer.saz_entrance));
        REQUIRE(layer.entrance.lower_alt_index >= 0);
        REQUIRE(layer.exit.lower_alt_index >= 0);
    }
}

TEST_CASE("RustRayTracer2D preserves each Geometry2D altitude interpolation mode",
          "[raytracing][rust][geometry2d]") {
    const auto interpolation = GENERATE(
        sasktran2::grids::interpolation::linear,
        sasktran2::grids::interpolation::shell,
        sasktran2::grids::interpolation::lower);
    auto geo = geometry(interpolation);
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay2D traced;
    tracer.trace_ray(ray(25.0 * radial_direction(-0.25),
                         radial_direction(-0.25)),
                     traced);

    REQUIRE(traced.layers.size() == 1);
    const auto& weights = traced.layers[0].entrance.interpolation_weights;
    REQUIRE(weights.size() ==
            (interpolation == sasktran2::grids::interpolation::lower ? 2 : 4));
    double sum = 0.0;
    for (const auto& [index, weight] : weights) {
        CAPTURE(index, weight);
        sum += weight;
    }
    require_close(sum, 1.0);
    require_weights_match_geometry(geo, traced.layers[0].entrance);
}

TEST_CASE("RustRayTracer2D splits rays at internal horizontal nodes",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay2D traced;
    tracer.trace_ray(ray(15.0 * radial_direction(-0.25),
                         Eigen::Vector3d::UnitX()),
                     traced);

    REQUIRE(traced.layers.size() == 3);
    REQUIRE(traced.layers[0].horizontal_cell == 1);
    REQUIRE(traced.layers[1].horizontal_cell == 1);
    REQUIRE(traced.layers[2].horizontal_cell == 0);
    REQUIRE(traced.layers[0].altitude_cell == 1);
    REQUIRE(traced.layers[1].altitude_cell == 0);
    REQUIRE(traced.layers[2].altitude_cell == 0);
    REQUIRE((traced.layers[1].entrance.position -
             traced.layers[2].exit.position)
                .norm() <= tolerance);
    require_close(traced.layers[1].entrance.position.x(), 0.0);
}

TEST_CASE("RustRayTracer2D horizontally extends edge cells and clamps weights",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay2D traced;

    SECTION("left edge") {
        tracer.trace_ray(ray(40.0 * radial_direction(-0.75),
                             -radial_direction(-0.75)),
                         traced);
        REQUIRE(traced.ground_is_hit);
        REQUIRE(traced.layers.size() == 2);
        for (const auto& layer : traced.layers) {
            REQUIRE(layer.horizontal_cell == 0);
            for (const auto& [index, weight] :
                 layer.entrance.interpolation_weights) {
                REQUIRE(index < geo.num_altitudes());
                REQUIRE(weight >= 0.0);
            }
        }
    }

    SECTION("right edge") {
        tracer.trace_ray(ray(40.0 * radial_direction(0.75),
                             -radial_direction(0.75)),
                         traced);
        REQUIRE(traced.ground_is_hit);
        REQUIRE(traced.layers.size() == 2);
        for (const auto& layer : traced.layers) {
            REQUIRE(layer.horizontal_cell == 1);
            for (const auto& [index, weight] :
                 layer.entrance.interpolation_weights) {
                REQUIRE(index >= 2 * geo.num_altitudes());
                REQUIRE(weight >= 0.0);
            }
        }
    }
}

TEST_CASE("RustRayTracer2D splits extended edge cells at the angular seam",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay2D traced;
    tracer.trace_ray(ray(15.0 * radial_direction(3.0),
                         -Eigen::Vector3d::UnitX()),
                     traced);

    bool found_seam = false;
    for (std::size_t index = 1; index < traced.layers.size(); ++index) {
        const auto& far = traced.layers[index - 1];
        const auto& near = traced.layers[index];
        if (far.horizontal_cell == 0 && near.horizontal_cell == 1 &&
            far.entrance.position.z() < 0.0 &&
            (far.entrance.position - near.exit.position).norm() <= tolerance) {
            require_close(far.entrance.position.x(), 0.0);
            found_seam = true;
        }
    }
    REQUIRE(found_seam);
}

TEST_CASE("RustRayTracer2D handles observers inside and outside the altitude domain",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay2D traced;

    tracer.trace_ray(ray(40.0 * radial_direction(0.25),
                         radial_direction(0.25)),
                     traced);
    REQUIRE(traced.layers.empty());
    REQUIRE_FALSE(traced.ground_is_hit);

    tracer.trace_ray(ray(25.0 * radial_direction(-0.25),
                         radial_direction(-0.25)),
                     traced);
    REQUIRE(traced.layers.size() == 1);
    REQUIRE_FALSE(traced.ground_is_hit);
    REQUIRE(traced.layers[0].type == sasktran2::raytracing::partial);
    require_close(traced.layers[0].layer_distance, 5.0);
}

TEST_CASE("RustRayTracer2D respects a rotated Geometry2D basis",
          "[raytracing][rust][geometry2d]") {
    sasktran2::Coordinates coordinates(
        Eigen::Vector3d::UnitX(), -Eigen::Vector3d::UnitZ(),
        Eigen::Vector3d::UnitY(), earth_radius);
    sasktran2::Geometry2D geo(std::move(coordinates),
                             values({0.0, 10.0, 20.0}),
                             values({-0.5, 0.0, 0.5}));
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay2D traced;
    const double angle = 0.25;
    const Eigen::Vector3d direction =
        geo.coordinates().unit_vector_from_angles(angle, 0.0);
    tracer.trace_ray(ray(40.0 * direction, -direction), traced);

    REQUIRE(traced.ground_is_hit);
    REQUIRE(traced.layers.size() == 2);
    REQUIRE(traced.layers[0].horizontal_cell == 1);
    REQUIRE(traced.layers[1].horizontal_cell == 1);
}

TEST_CASE("RustRayTracer2D supports a different altitude-only refractive profile per ray",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    const double tangent_radius = 15.0;
    const double observer_radius = 40.0;
    const auto limb_ray = ray(
        {0.0, 0.0, observer_radius},
        {tangent_radius / observer_radius, 0.0,
         -std::sqrt(1.0 - std::pow(tangent_radius / observer_radius, 2))});

    sasktran2::raytracing::TracedRay2D weak;
    sasktran2::raytracing::TracedRay2D strong;
    tracer.trace_ray(limb_ray, values({1.02, 1.01, 1.0}), weak);
    tracer.trace_ray(limb_ray, values({1.04, 1.02, 1.0}), strong);

    REQUIRE_FALSE(weak.is_straight);
    REQUIRE_FALSE(strong.is_straight);
    REQUIRE_FALSE(weak.layers.empty());
    REQUIRE_FALSE(strong.layers.empty());
    REQUIRE(std::abs(weak.tangent_radius - strong.tangent_radius) > 1e-6);
    REQUIRE(weak.layers.front().curvature_factor !=
            strong.layers.front().curvature_factor);
}

TEST_CASE("RustRayTracer2D unity refraction preserves straight geometry",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    const auto viewing_ray =
        ray({-40.0, 0.0, 20.0}, Eigen::Vector3d::UnitX());
    sasktran2::raytracing::TracedRay2D straight;
    sasktran2::raytracing::TracedRay2D unity;
    tracer.trace_ray(viewing_ray, straight);
    tracer.trace_ray(viewing_ray, values({1.0, 1.0, 1.0}), unity);

    REQUIRE(straight.is_straight);
    REQUIRE_FALSE(unity.is_straight);
    REQUIRE(straight.layers.size() == unity.layers.size());
    for (std::size_t index = 0; index < straight.layers.size(); ++index) {
        require_close(straight.layers[index].layer_distance,
                      unity.layers[index].layer_distance);
        REQUIRE(straight.layers[index].altitude_cell ==
                unity.layers[index].altitude_cell);
        REQUIRE(straight.layers[index].horizontal_cell ==
                unity.layers[index].horizontal_cell);
    }
}

TEST_CASE("RustRayTracer2D validates per-ray refractive profiles",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay2D traced;
    const auto viewing_ray =
        ray({-40.0, 0.0, 20.0}, Eigen::Vector3d::UnitX());

    REQUIRE_THROWS_AS(tracer.trace_ray(viewing_ray, values({1.0, 1.0}), traced),
                      std::invalid_argument);
    REQUIRE_THROWS(tracer.trace_ray(
        viewing_ray,
        values({1.0, std::numeric_limits<double>::quiet_NaN(), 1.0}), traced));
}

TEST_CASE("RustRayTracer2D reuses output storage without stale layers",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay2D traced;

    tracer.trace_ray(ray(40.0 * radial_direction(0.25),
                         -radial_direction(0.25)),
                     traced);
    REQUIRE(traced.layers.size() == 2);
    tracer.trace_ray(ray(40.0 * radial_direction(0.25),
                         radial_direction(0.25)),
                     traced);
    REQUIRE(traced.layers.empty());
    REQUIRE_FALSE(traced.ground_is_hit);
}

#endif
