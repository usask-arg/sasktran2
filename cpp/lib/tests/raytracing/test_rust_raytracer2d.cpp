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

    sasktran2::Geometry2D
    geometry(sasktran2::grids::interpolation interpolation =
                 sasktran2::grids::interpolation::linear) {
        return sasktran2::Geometry2D(0.6, 0.3, earth_radius,
                                     values({0.0, 10.0, 20.0}),
                                     values({-0.5, 0.0, 0.5}), interpolation);
    }

    Eigen::Vector3d radial_direction(double angle) {
        return {std::sin(angle), 0.0, std::cos(angle)};
    }

    sasktran2::viewinggeometry::ViewingRay ray(const Eigen::Vector3d& observer,
                                               const Eigen::Vector3d& look) {
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

    std::pair<int, int>
    layer_cell(const sasktran2::Geometry2D& geometry,
               const sasktran2::raytracing::TracedLayer& layer) {
        sasktran2::Location midpoint;
        midpoint.position =
            0.5 * (layer.entrance.position + layer.exit.position);
        const auto cell = geometry.cell_indices(midpoint);
        REQUIRE(cell.has_value());
        return *cell;
    }

    void require_weights_match_location(
        const sasktran2::Geometry2D& geometry,
        const sasktran2::Location& location,
        const sasktran2::raytracing::GridWeightStencilView& actual) {
        std::vector<std::pair<int, double>> expected;
        geometry.assign_interpolation_weights(location, expected);
        expected.erase(std::remove_if(expected.begin(), expected.end(),
                                      [](const auto& weight) {
                                          return weight.second == 0.0;
                                      }),
                       expected.end());
        std::vector<std::pair<int, double>> actual_nonzero;
        double weight_sum = 0.0;
        for (std::size_t index = 0; index < actual.size(); ++index) {
            weight_sum += actual[index].second;
            if (actual[index].second != 0.0) {
                actual_nonzero.push_back(actual[index]);
            }
        }
        REQUIRE(actual_nonzero.size() == expected.size());
        for (std::size_t index = 0; index < expected.size(); ++index) {
            REQUIRE(actual_nonzero[index].first == expected[index].first);
            require_close(actual_nonzero[index].second, expected[index].second);
        }
        require_close(weight_sum, 1.0);
    }
} // namespace

TEST_CASE("RustRayTracer2D traces radial layers and exposes 2D cells",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay traced;
    const double angle = 0.25;
    tracer.trace_ray(
        ray(40.0 * radial_direction(angle), -radial_direction(angle)), traced);

    REQUIRE(traced.ground_is_hit);
    REQUIRE(traced.is_straight);
    REQUIRE(traced.layers.size() == 2);
    REQUIRE(layer_cell(geo, traced.layers[0]) == std::make_pair(0, 1));
    REQUIRE(layer_cell(geo, traced.layers[1]) == std::make_pair(1, 1));
    require_close(traced.layers[0].layer_distance, 10.0);
    require_close(traced.layers[1].layer_distance, 10.0);

    for (std::size_t layer_index = 0; layer_index < traced.layers.size();
         ++layer_index) {
        const auto& layer = traced.layers[layer_index];
        require_weights_match_location(geo, layer.entrance,
                                       traced.entrance_weights(layer_index));
        require_weights_match_location(geo, layer.exit,
                                       traced.exit_weights(layer_index));
        const auto od_weights = traced.optical_depth_weights(layer_index);
        double od_weight_sum = 0.0;
        for (std::size_t index = 0; index < od_weights.size(); ++index) {
            od_weight_sum += od_weights[index].second;
        }
        require_close(od_weight_sum,
                      layer.layer_distance * layer.curvature_factor);
        REQUIRE(std::isfinite(layer.cos_sza_entrance));
        REQUIRE(std::isfinite(layer.saz_entrance));
        REQUIRE(layer.entrance.lower_alt_index >= 0);
        REQUIRE(layer.exit.lower_alt_index >= 0);
    }
}

TEST_CASE(
    "RustRayTracer2D preserves each Geometry2D altitude interpolation mode",
    "[raytracing][rust][geometry2d]") {
    const auto interpolation = GENERATE(sasktran2::grids::interpolation::linear,
                                        sasktran2::grids::interpolation::shell,
                                        sasktran2::grids::interpolation::lower);
    auto geo = geometry(interpolation);
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay traced;
    tracer.trace_ray(
        ray(25.0 * radial_direction(-0.25), radial_direction(-0.25)), traced);

    REQUIRE(traced.layers.size() == 1);
    const auto& layer = traced.layers[0];
    const auto weights = traced.entrance_weights(0);
    double weight_sum = 0.0;
    for (std::size_t index = 0; index < weights.size(); ++index) {
        weight_sum += weights[index].second;
    }
    require_close(weight_sum, 1.0);
    if (interpolation == sasktran2::grids::interpolation::lower) {
        require_close(weights[1].second, 0.0);
        require_close(weights[3].second, 0.0);
    } else if (interpolation == sasktran2::grids::interpolation::shell) {
        require_close(weights[1].second + weights[3].second, 0.5);
    }
    require_weights_match_location(geo, layer.entrance, weights);
}

TEST_CASE("RustRayTracer2D exposes conservative four-node OD weights",
          "[raytracing][rust][geometry2d]") {
    const auto interpolation = GENERATE(sasktran2::grids::interpolation::linear,
                                        sasktran2::grids::interpolation::shell,
                                        sasktran2::grids::interpolation::lower);
    const bool refracted = GENERATE(false, true);
    CAPTURE(interpolation, refracted);

    auto geo = geometry(interpolation);
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay traced;
    const auto oblique_ray =
        ray(15.0 * radial_direction(-0.25), Eigen::Vector3d::UnitX());
    if (refracted) {
        tracer.trace_ray(oblique_ray, values({1.02, 1.01, 1.0}), traced);
    } else {
        tracer.trace_ray(oblique_ray, traced);
    }

    const std::array<double, 4> field = {0.3, 1.1, 2.7, 4.2};
    for (std::size_t layer_index = 0; layer_index < traced.layers.size();
         ++layer_index) {
        const auto& layer = traced.layers[layer_index];
        const auto od_weights = traced.optical_depth_weights(layer_index);
        double weight_sum = 0.0;
        for (std::size_t index = 0; index < od_weights.size(); ++index) {
            weight_sum += od_weights[index].second;
        }
        require_close(weight_sum, layer.layer_distance * layer.curvature_factor,
                      2e-8);
        double optical_depth = 0.0;
        for (std::size_t index = 0; index < od_weights.size(); ++index) {
            REQUIRE(std::isfinite(od_weights[index].second));
            REQUIRE(od_weights[index].second >= 0.0);
            optical_depth += od_weights[index].second * field[index];
        }
        REQUIRE(std::isfinite(optical_depth));
        for (std::size_t index = 0; index < field.size(); ++index) {
            auto perturbed = field;
            constexpr double delta = 1e-6;
            perturbed[index] += delta;
            double perturbed_od = 0.0;
            for (std::size_t weight_index = 0; weight_index < od_weights.size();
                 ++weight_index) {
                perturbed_od +=
                    od_weights[weight_index].second * perturbed[weight_index];
            }
            require_close((perturbed_od - optical_depth) / delta,
                          od_weights[index].second, 2e-8);
        }
    }

    Eigen::SparseMatrix<double, Eigen::RowMajor> matrix;
    sasktran2::raytracing::construct_od_matrix(traced, geo, matrix);
    REQUIRE(matrix.rows() == traced.layers.size());
    REQUIRE(matrix.cols() == geo.size());
    for (int layer_index = 0; layer_index < traced.layers.size();
         ++layer_index) {
        const auto od_weights = traced.optical_depth_weights(layer_index);
        REQUIRE(od_weights.size() == 4);
        for (std::size_t local_index = 0; local_index < od_weights.size();
             ++local_index) {
            require_close(
                matrix.coeff(layer_index, od_weights[local_index].first),
                od_weights[local_index].second, 1e-12);
        }
    }

    Eigen::MatrixXd spectra(geo.size(), 4);
    for (int location_index = 0; location_index < spectra.rows();
         ++location_index) {
        for (int wavelength = 0; wavelength < spectra.cols(); ++wavelength) {
            spectra(location_index, wavelength) =
                1e-6 * (1.0 + 0.13 * location_index + 0.21 * wavelength);
        }
    }
    const Eigen::MatrixXd optical_depth = matrix * spectra;
    REQUIRE(optical_depth.rows() == traced.layers.size());
    REQUIRE(optical_depth.cols() == spectra.cols());
    REQUIRE(optical_depth.allFinite());
}

TEST_CASE("RustRayTracer2D splits rays at internal horizontal nodes",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay traced;
    tracer.trace_ray(
        ray(15.0 * radial_direction(-0.25), Eigen::Vector3d::UnitX()), traced);

    REQUIRE(traced.layers.size() == 3);
    REQUIRE(layer_cell(geo, traced.layers[0]) == std::make_pair(1, 1));
    REQUIRE(layer_cell(geo, traced.layers[1]) == std::make_pair(0, 1));
    REQUIRE(layer_cell(geo, traced.layers[2]) == std::make_pair(0, 0));
    REQUIRE(
        (traced.layers[1].entrance.position - traced.layers[2].exit.position)
            .norm() <= tolerance);
    require_close(traced.layers[1].entrance.position.x(), 0.0);
}

TEST_CASE("RustRayTracer2D horizontally extends edge cells and clamps weights",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay traced;

    SECTION("left edge") {
        tracer.trace_ray(
            ray(40.0 * radial_direction(-0.75), -radial_direction(-0.75)),
            traced);
        REQUIRE(traced.ground_is_hit);
        REQUIRE(traced.layers.size() == 2);
        for (std::size_t index = 0; index < traced.layers.size(); ++index) {
            REQUIRE(layer_cell(geo, traced.layers[index]).second == 0);
            const auto weights = traced.entrance_weights(index);
            require_close(weights[2].second + weights[3].second, 0.0);
        }
    }

    SECTION("right edge") {
        tracer.trace_ray(
            ray(40.0 * radial_direction(0.75), -radial_direction(0.75)),
            traced);
        REQUIRE(traced.ground_is_hit);
        REQUIRE(traced.layers.size() == 2);
        for (std::size_t index = 0; index < traced.layers.size(); ++index) {
            REQUIRE(layer_cell(geo, traced.layers[index]).second == 1);
            const auto weights = traced.entrance_weights(index);
            require_close(weights[2].second + weights[3].second, 1.0);
        }
    }
}

TEST_CASE("RustRayTracer2D splits extended edge cells at the angular seam",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay traced;
    tracer.trace_ray(
        ray(15.0 * radial_direction(3.0), -Eigen::Vector3d::UnitX()), traced);

    bool found_seam = false;
    for (std::size_t index = 1; index < traced.layers.size(); ++index) {
        const auto& before_seam = traced.layers[index - 1];
        const auto& after_seam = traced.layers[index];
        if (layer_cell(geo, before_seam).second == 0 &&
            layer_cell(geo, after_seam).second == 1 &&
            before_seam.entrance.position.z() < 0.0 &&
            (before_seam.entrance.position - after_seam.exit.position).norm() <=
                tolerance) {
            require_close(before_seam.entrance.position.x(), 0.0);
            found_seam = true;
        }
    }
    REQUIRE(found_seam);
}

TEST_CASE(
    "RustRayTracer2D handles observers inside and outside the altitude domain",
    "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay traced;

    tracer.trace_ray(ray(40.0 * radial_direction(0.25), radial_direction(0.25)),
                     traced);
    REQUIRE(traced.layers.empty());
    REQUIRE_FALSE(traced.ground_is_hit);

    tracer.trace_ray(
        ray(25.0 * radial_direction(-0.25), radial_direction(-0.25)), traced);
    REQUIRE(traced.layers.size() == 1);
    REQUIRE_FALSE(traced.ground_is_hit);
    REQUIRE(traced.layers[0].type == sasktran2::raytracing::partial);
    require_close(traced.layers[0].layer_distance, 5.0);
}

TEST_CASE("RustRayTracer2D respects a rotated Geometry2D basis",
          "[raytracing][rust][geometry2d]") {
    sasktran2::Coordinates coordinates(Eigen::Vector3d::UnitX(),
                                       -Eigen::Vector3d::UnitZ(),
                                       Eigen::Vector3d::UnitY(), earth_radius);
    sasktran2::Geometry2D geo(std::move(coordinates), values({0.0, 10.0, 20.0}),
                              values({-0.5, 0.0, 0.5}));
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay traced;
    const double angle = 0.25;
    const Eigen::Vector3d direction =
        geo.coordinates().unit_vector_from_angles(angle, 0.0);
    tracer.trace_ray(ray(40.0 * direction, -direction), traced);

    REQUIRE(traced.ground_is_hit);
    REQUIRE(traced.layers.size() == 2);
    REQUIRE(layer_cell(geo, traced.layers[0]).second == 1);
    REQUIRE(layer_cell(geo, traced.layers[1]).second == 1);
}

TEST_CASE("RustRayTracer2D supports a different altitude-only refractive "
          "profile per ray",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    const double tangent_radius = 15.0;
    const double observer_radius = 40.0;
    const auto limb_ray =
        ray({0.0, 0.0, observer_radius},
            {tangent_radius / observer_radius, 0.0,
             -std::sqrt(1.0 - std::pow(tangent_radius / observer_radius, 2))});

    sasktran2::raytracing::TracedRay weak;
    sasktran2::raytracing::TracedRay strong;
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
    const auto viewing_ray = ray({-40.0, 0.0, 20.0}, Eigen::Vector3d::UnitX());
    sasktran2::raytracing::TracedRay straight;
    sasktran2::raytracing::TracedRay unity;
    tracer.trace_ray(viewing_ray, straight);
    tracer.trace_ray(viewing_ray, values({1.0, 1.0, 1.0}), unity);

    REQUIRE(straight.is_straight);
    REQUIRE_FALSE(unity.is_straight);
    REQUIRE(straight.layers.size() == unity.layers.size());
    for (std::size_t index = 0; index < straight.layers.size(); ++index) {
        require_close(straight.layers[index].layer_distance,
                      unity.layers[index].layer_distance);
        REQUIRE(layer_cell(geo, straight.layers[index]) ==
                layer_cell(geo, unity.layers[index]));
    }
}

TEST_CASE("RustRayTracer2D validates per-ray refractive profiles",
          "[raytracing][rust][geometry2d]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::raytracing::TracedRay traced;
    const auto viewing_ray = ray({-40.0, 0.0, 20.0}, Eigen::Vector3d::UnitX());

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
    sasktran2::raytracing::TracedRay traced;

    tracer.trace_ray(
        ray(40.0 * radial_direction(0.25), -radial_direction(0.25)), traced);
    REQUIRE(traced.layers.size() == 2);
    tracer.trace_ray(ray(40.0 * radial_direction(0.25), radial_direction(0.25)),
                     traced);
    REQUIRE(traced.layers.empty());
    REQUIRE_FALSE(traced.ground_is_hit);
}

#endif
