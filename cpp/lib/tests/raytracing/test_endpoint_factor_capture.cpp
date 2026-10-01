#include "../../successive_orders/geometry.h"

#include <sasktran2/test_helper.h>

#include <sasktran2.h>

#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#ifdef SKTRAN_RUST_SUPPORT

namespace {
    using Factors = sasktran2::raytracing::LayerEndpointFactors2D;
    using TracedRay = sasktran2::raytracing::TracedRay;

    Eigen::VectorXd values(std::initializer_list<double> input) {
        Eigen::VectorXd result(input.size());
        std::size_t index = 0;
        for (double value : input) {
            result[index++] = value;
        }
        return result;
    }

    sasktran2::Geometry2D
    geometry(sasktran2::grids::interpolation interpolation =
                 sasktran2::grids::interpolation::linear) {
        return sasktran2::Geometry2D(0.6, 0.3, 10.0,
                                     values({0.0, 10.0, 20.0, 30.0}),
                                     values({-0.5, 0.0, 0.5}), interpolation);
    }

    Eigen::Vector3d radial(double angle) {
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

    std::uint64_t bits(double value) {
        std::uint64_t result;
        static_assert(sizeof(result) == sizeof(value));
        std::memcpy(&result, &value, sizeof(result));
        return result;
    }

    void require_bits(double actual, double expected) {
        CAPTURE(actual, expected);
        REQUIRE(bits(actual) == bits(expected));
    }

    void require_stencil_equal(
        const sasktran2::raytracing::GridWeightStencilView& actual,
        const sasktran2::raytracing::GridWeightStencilView& expected) {
        REQUIRE(actual.size() == expected.size());
        for (std::size_t index = 0; index < actual.size(); ++index) {
            REQUIRE(actual[index].first == expected[index].first);
            require_bits(actual[index].second, expected[index].second);
        }
    }

    void require_trace_equal(const TracedRay& actual,
                             const TracedRay& expected) {
        REQUIRE(actual.is_straight == expected.is_straight);
        REQUIRE(actual.ground_is_hit == expected.ground_is_hit);
        REQUIRE(actual.layers.size() == expected.layers.size());
        require_bits(actual.tangent_radius, expected.tangent_radius);
        for (std::size_t index = 0; index < actual.layers.size(); ++index) {
            const auto& layer = actual.layers[index];
            const auto& reference = expected.layers[index];
            REQUIRE(layer.type == reference.type);
            for (int axis = 0; axis < 3; ++axis) {
                require_bits(layer.entrance.position[axis],
                             reference.entrance.position[axis]);
                require_bits(layer.exit.position[axis],
                             reference.exit.position[axis]);
                require_bits(layer.average_look_away[axis],
                             reference.average_look_away[axis]);
            }
            require_bits(layer.layer_distance, reference.layer_distance);
            require_bits(layer.curvature_factor, reference.curvature_factor);
            require_bits(layer.od_quad_start, reference.od_quad_start);
            require_bits(layer.od_quad_end, reference.od_quad_end);
            require_bits(layer.od_quad_start_fraction,
                         reference.od_quad_start_fraction);
            require_bits(layer.od_quad_end_fraction,
                         reference.od_quad_end_fraction);
            require_bits(layer.cos_sza_entrance, reference.cos_sza_entrance);
            require_bits(layer.cos_sza_exit, reference.cos_sza_exit);
            require_bits(layer.saz_entrance, reference.saz_entrance);
            require_bits(layer.saz_exit, reference.saz_exit);
            require_stencil_equal(actual.entrance_weights(index),
                                  expected.entrance_weights(index));
            require_stencil_equal(actual.exit_weights(index),
                                  expected.exit_weights(index));
            require_stencil_equal(actual.optical_depth_weights(index),
                                  expected.optical_depth_weights(index));
        }
    }

    void require_original_factors(
        const sasktran2::Geometry2D& geo, const sasktran2::Location& location,
        const sasktran2::raytracing::GridWeightStencilView& stencil,
        const std::array<double, 2>& factors) {
        REQUIRE(stencil.size() == 4);
        const int base = stencil[0].first;
        const auto coordinates = geo.cell_interpolation_coordinates(
            location, base % geo.num_altitudes(), base / geo.num_altitudes());
        require_bits(factors[0], coordinates.first);
        require_bits(factors[1], coordinates.second);
        const double altitude_lower = 1.0 - factors[0];
        const double horizontal_lower = 1.0 - factors[1];
        const std::array<double, 4> expanded{
            horizontal_lower * altitude_lower, horizontal_lower * factors[0],
            factors[1] * altitude_lower, factors[1] * factors[0]};
        for (std::size_t index = 0; index < expanded.size(); ++index) {
            require_bits(expanded[index], stencil[index].second);
        }
    }

    void require_capture(const sasktran2::Geometry2D& geo,
                         const TracedRay& traced,
                         const std::vector<Factors>& factors) {
        REQUIRE(factors.size() == traced.layers.size());
        for (std::size_t index = 0; index < factors.size(); ++index) {
            require_original_factors(geo, traced.layers[index].entrance,
                                     traced.entrance_weights(index),
                                     factors[index].entrance);
            require_original_factors(geo, traced.layers[index].exit,
                                     traced.exit_weights(index),
                                     factors[index].exit);
        }
    }

    void require_factor_sidecar_equal(
        const std::vector<std::vector<Factors>>& actual,
        const std::vector<std::vector<Factors>>& expected) {
        REQUIRE(actual.size() == expected.size());
        for (std::size_t ray_index = 0; ray_index < actual.size();
             ++ray_index) {
            REQUIRE(actual[ray_index].size() == expected[ray_index].size());
            for (std::size_t layer = 0; layer < actual[ray_index].size();
                 ++layer) {
                for (std::size_t coordinate = 0; coordinate < 2; ++coordinate) {
                    require_bits(
                        actual[ray_index][layer].entrance[coordinate],
                        expected[ray_index][layer].entrance[coordinate]);
                    require_bits(actual[ray_index][layer].exit[coordinate],
                                 expected[ray_index][layer].exit[coordinate]);
                }
            }
        }
    }
} // namespace

TEST_CASE("2D endpoint capture preserves the ordinary trace and original "
          "factor bits",
          "[raytracing][rust][geometry2d][endpoint_factors]") {
    const auto interpolation = GENERATE(sasktran2::grids::interpolation::linear,
                                        sasktran2::grids::interpolation::shell,
                                        sasktran2::grids::interpolation::lower);
    auto geo = geometry(interpolation);
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    const std::array<sasktran2::viewinggeometry::ViewingRay, 4> inputs{
        ray(50.0 * radial(0.25), -radial(0.25)),
        ray(25.0 * radial(-0.25), radial(-0.25)),
        ray(15.0 * radial(-0.75), Eigen::Vector3d::UnitX()),
        ray(Eigen::Vector3d(-50.0, 0.0, 20.0), Eigen::Vector3d::UnitX())};
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        CAPTURE(interpolation, index);
        TracedRay ordinary;
        TracedRay captured;
        std::vector<Factors> factors;
        tracer.trace_ray(inputs[index], ordinary);
        tracer.trace_ray_with_endpoint_factors(inputs[index], captured,
                                               factors);
        REQUIRE(!ordinary.layers.empty());
        require_trace_equal(captured, ordinary);
        require_capture(geo, captured, factors);
    }
}

TEST_CASE("2D endpoint capture replaces its independent construction sidecar",
          "[raytracing][rust][geometry2d][endpoint_factors]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    const auto inward = ray(50.0 * radial(0.25), -radial(0.25));
    const auto outward = ray(50.0 * radial(0.25), radial(0.25));
    TracedRay traced;
    std::vector<Factors> factors(100, Factors{{0.123, 0.456}, {0.789, 0.234}});
    tracer.trace_ray_with_endpoint_factors(inward, traced, factors);
    REQUIRE(!factors.empty());
    require_capture(geo, traced, factors);
    const std::vector<std::vector<Factors>> saved{factors};

    tracer.trace_ray(inward, traced);
    require_factor_sidecar_equal({factors}, saved);
    tracer.trace_ray_optical_depth(inward, traced);
    require_factor_sidecar_equal({factors}, saved);

    tracer.trace_ray_with_endpoint_factors(outward, traced, factors);
    REQUIRE(traced.layers.empty());
    REQUIRE(factors.empty());
    tracer.trace_ray_with_endpoint_factors(inward, traced, factors);
    require_factor_sidecar_equal({factors}, saved);
    require_capture(geo, traced, factors);
}

TEST_CASE("Source endpoint sidecars survive LOS refresh and release with "
          "incoming traced rays",
          "[successive_orders][geometry2d][endpoint_factors]") {
    const int threads = GENERATE(1, 2);
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.num_sza = 2;
    settings.num_threads = threads;
    const sasktran2::viewinggeometry::InternalViewingGeometry empty_viewing;
    sasktran2::successive_orders::SourceGeometry1D source(tracer, geo);
    source.initialize(empty_viewing, settings);
    REQUIRE(source.incoming_endpoint_factors().empty());

    source.initialize(empty_viewing, settings, true);
    REQUIRE(source.incoming_endpoint_factors().size() ==
            source.incoming_rays().size());
    REQUIRE(!source.incoming_endpoint_factors().empty());
    for (std::size_t index = 0; index < source.incoming_rays().size();
         ++index) {
        require_capture(geo, source.incoming_rays()[index],
                        source.incoming_endpoint_factors()[index]);
    }
    const auto saved = source.incoming_endpoint_factors();
    const auto saved_offsets = source.transport_row_offsets();
    const auto saved_columns = source.transport_column_indices().to_vector();
    sasktran2::viewinggeometry::InternalViewingGeometry new_viewing;
    new_viewing.traced_rays.resize(1);
    tracer.trace_ray(ray(50.0 * radial(0.25), -radial(0.25)),
                     new_viewing.traced_rays.front());
    source.refresh_los(new_viewing);
    require_factor_sidecar_equal(source.incoming_endpoint_factors(), saved);

    source.release_incoming_traced_rays();
    REQUIRE(source.incoming_rays().empty());
    REQUIRE(source.incoming_endpoint_factors().empty());
    REQUIRE(source.incoming_endpoint_factors().capacity() == 0);
    REQUIRE(source.incoming_interpolation().size() == saved.size());
    REQUIRE(source.transport_row_offsets() == saved_offsets);
    REQUIRE(source.transport_column_indices().to_vector() == saved_columns);

    source.initialize(new_viewing, settings);
    REQUIRE(source.incoming_endpoint_factors().empty());
    source.initialize(new_viewing, settings, true);
    require_factor_sidecar_equal(source.incoming_endpoint_factors(), saved);

    // This fails after diffuse capture, during LOS interpolation construction.
    auto invalid_viewing = new_viewing;
    invalid_viewing.traced_rays.front().layers.front().entrance.position.x() =
        std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS(source.initialize(invalid_viewing, settings, true));
    REQUIRE(source.incoming_endpoint_factors().empty());
    REQUIRE(source.incoming_endpoint_factors().capacity() == 0);
    source.initialize(new_viewing, settings, true);
    require_factor_sidecar_equal(source.incoming_endpoint_factors(), saved);
}

TEST_CASE("Invalid source settings release captured endpoint factors before "
          "geometry construction",
          "[successive_orders][geometry2d][endpoint_factors]") {
    auto geo = geometry();
    sasktran2::raytracing::RustRayTracer2D tracer(geo);
    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.num_sza = 2;
    const sasktran2::viewinggeometry::InternalViewingGeometry empty_viewing;
    sasktran2::successive_orders::SourceGeometry1D source(tracer, geo);
    source.initialize(empty_viewing, settings, true);
    const auto saved_factors = source.incoming_endpoint_factors();
    REQUIRE(!saved_factors.empty());
    REQUIRE(source.incoming_endpoint_factors().capacity() > 0);

    auto invalid_settings = settings;
    invalid_settings.num_threads = 0;
    REQUIRE_THROWS_AS(source.initialize(empty_viewing, invalid_settings, true),
                      std::invalid_argument);
    REQUIRE(source.incoming_endpoint_factors().empty());
    REQUIRE(source.incoming_endpoint_factors().capacity() == 0);

    source.initialize(empty_viewing, settings, true);
    sasktran2::successive_orders::SourceGeometry1D reference(tracer, geo);
    reference.initialize(empty_viewing, settings, true);
    require_factor_sidecar_equal(source.incoming_endpoint_factors(),
                                 saved_factors);
    require_factor_sidecar_equal(source.incoming_endpoint_factors(),
                                 reference.incoming_endpoint_factors());
    REQUIRE(source.incoming_rays().size() == reference.incoming_rays().size());
    for (std::size_t ray_index = 0; ray_index < source.incoming_rays().size();
         ++ray_index) {
        require_trace_equal(source.incoming_rays()[ray_index],
                            reference.incoming_rays()[ray_index]);
        require_capture(geo, source.incoming_rays()[ray_index],
                        source.incoming_endpoint_factors()[ray_index]);
    }
    REQUIRE(source.transport_row_offsets() ==
            reference.transport_row_offsets());
    REQUIRE(source.transport_column_indices().to_vector() ==
            reference.transport_column_indices().to_vector());
}

TEST_CASE("1D source geometry ignores the optional 2D factor capture",
          "[successive_orders][geometry][endpoint_factors]") {
    sasktran2::Geometry1D geo(0.6, 0.3, 6372000.0,
                              values({0.0, 1000.0, 3000.0}));
    sasktran2::raytracing::SphericalShellRayTracer tracer(geo);
    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    sasktran2::successive_orders::SourceGeometry1D source(tracer, geo);
    const sasktran2::viewinggeometry::InternalViewingGeometry empty_viewing;
    source.initialize(empty_viewing, settings, true);
    REQUIRE(!source.incoming_rays().empty());
    REQUIRE(source.incoming_endpoint_factors().empty());
    REQUIRE(source.incoming_endpoint_factors().capacity() == 0);
}

#endif
