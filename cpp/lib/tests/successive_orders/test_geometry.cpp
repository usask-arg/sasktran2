#include "../../successive_orders/geometry.h"

#include <sasktran2/solartransmission.h>
#include <sasktran2/test_helper.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace {
    Eigen::VectorXd altitude_grid() {
        return (Eigen::Vector3d() << 0.0, 1000.0, 3000.0).finished();
    }

    sasktran2::viewinggeometry::InternalViewingGeometry
    make_los_geometry(const sasktran2::Geometry1D& geometry,
                      const sasktran2::raytracing::RayTracerBase& raytracer) {
        sasktran2::viewinggeometry::InternalViewingGeometry result;
        result.traced_rays.resize(1);
        sasktran2::viewinggeometry::ViewingRay ray;
        ray.observer.position = geometry.coordinates().reference_point(4000.0);
        ray.look_away = -ray.observer.position.normalized();
        raytracer.trace_ray(ray, result.traced_rays.front());
        return result;
    }

#ifdef SKTRAN_RUST_SUPPORT
    Eigen::VectorXd horizontal_angle_grid() {
        return (Eigen::Vector3d() << -0.4, 0.0, 0.4).finished();
    }

    sasktran2::viewinggeometry::InternalViewingGeometry
    make_los_geometry(const sasktran2::Geometry2D& geometry,
                      const sasktran2::raytracing::RustRayTracer2D& raytracer) {
        sasktran2::viewinggeometry::InternalViewingGeometry result;
        result.traced_rays.resize(1);
        sasktran2::viewinggeometry::ViewingRay ray;
        ray.observer.position =
            (geometry.coordinates().earth_radius() + 4000.0) *
            geometry.coordinates().unit_vector_from_angles(0.0, 0.0);
        ray.look_away = -ray.observer.position.normalized();
        raytracer.trace_ray(ray, result.traced_rays.front());
        return result;
    }
#endif

    class ThrowingRayTracer final
        : public sasktran2::raytracing::RayTracerBase {
      public:
        void trace_ray(const sasktran2::viewinggeometry::ViewingRay&,
                       sasktran2::raytracing::TracedRay&, bool) const override {
            throw std::runtime_error("deliberate incoming ray-tracing failure");
        }
    };

    sasktran2::viewinggeometry::InternalViewingGeometry
    make_exact_direction_los(const sasktran2::Geometry1D& geometry) {
        const std::array<Eigen::Vector3d, 4> directions{
            Eigen::Vector3d::UnitZ(), -Eigen::Vector3d::UnitZ(),
            Eigen::Vector3d::UnitX(), Eigen::Vector3d::UnitY()};
        const Eigen::Vector3d location =
            geometry.coordinates().reference_point(500.0);

        sasktran2::viewinggeometry::InternalViewingGeometry result;
        result.traced_rays.resize(directions.size());
        for (std::size_t ray_index = 0; ray_index < directions.size();
             ++ray_index) {
            auto& ray = result.traced_rays[ray_index];
            ray.observer_and_look.observer.position = location;
            ray.observer_and_look.look_away = directions[ray_index];
            ray.layers.resize(1);
            auto& layer = ray.layers.front();
            layer.entrance.position = location;
            layer.exit.position = location;
            layer.average_look_away = directions[ray_index];
            layer.cos_sza_entrance = 1.0;
            layer.cos_sza_exit = 1.0;
        }
        return result;
    }

    template <typename Weights> void require_sorted(const Weights& weights) {
        REQUIRE(std::is_sorted(weights.begin(), weights.end(),
                               [](const auto& left, const auto& right) {
                                   return left.index < right.index;
                               }));
    }

    void require_compiled_topology(
        const std::vector<sasktran2::successive_orders::RayInterpolation>& rays,
        const std::vector<int>& row_offsets,
        const std::vector<int>& column_indices, int num_source_columns) {
        REQUIRE(row_offsets.size() == rays.size() + 1);
        REQUIRE(row_offsets.front() == 0);
        REQUIRE(row_offsets.back() == static_cast<int>(column_indices.size()));
        for (std::size_t row = 0; row < rays.size(); ++row) {
            const auto& ray = rays[row];
            REQUIRE((ray.traced_ray != nullptr || ray.layers.empty() ||
                     !ray.optical_depth_weights.empty()));
            REQUIRE(ray.transport_value_offset ==
                    static_cast<std::size_t>(row_offsets[row]));
            REQUIRE(ray.transport_row_nnz ==
                    static_cast<std::uint32_t>(row_offsets[row + 1] -
                                               row_offsets[row]));
            const auto begin = column_indices.begin() + row_offsets[row];
            const auto end = column_indices.begin() + row_offsets[row + 1];
            REQUIRE(std::is_sorted(begin, end));
            REQUIRE(std::adjacent_find(begin, end) == end);
            for (std::size_t layer = 0; layer < ray.layers.size(); ++layer) {
                require_sorted(ray.atmosphere_for_layer(layer));
                const auto optical_depth = ray.optical_depth_for_layer(layer);
                for (std::size_t index = 1; index < optical_depth.size();
                     ++index) {
                    REQUIRE(optical_depth[index - 1].first <=
                            optical_depth[index].first);
                }
                const auto source = ray.source_for_layer(layer);
                for (const auto& weight : source) {
                    REQUIRE(weight.source_index >= 0);
                    REQUIRE(weight.source_index < num_source_columns);
                    REQUIRE(row_offsets[row] + weight.row_inner_index <
                            row_offsets[row + 1]);
                    REQUIRE(column_indices[row_offsets[row] +
                                           weight.row_inner_index] ==
                            weight.source_index);
                }
            }
            for (const auto& weight : ray.ground()) {
                REQUIRE(weight.source_index >= 0);
                REQUIRE(weight.source_index < num_source_columns);
                REQUIRE(
                    column_indices[row_offsets[row] + weight.row_inner_index] ==
                    weight.source_index);
            }
        }
    }

    void require_los_direction(
        const sasktran2::successive_orders::SourceGeometry1D& geometry,
        std::size_t ray_index, const Eigen::Vector3d& expected_direction) {
        const auto weights =
            geometry.los_interpolation()[ray_index].source_for_layer(0);
        REQUIRE(!weights.empty());

        double weight_sum = 0.0;
        for (const auto& weight : weights) {
            REQUIRE(std::isfinite(weight.weight));
            weight_sum += weight.weight;

            const sasktran2::successive_orders::SourcePoint* owner = nullptr;
            for (const auto& point : geometry.source_points()) {
                if (weight.source_index >= point.outgoing_offset() &&
                    weight.source_index <
                        point.outgoing_offset() + point.num_outgoing()) {
                    owner = &point;
                    break;
                }
            }
            REQUIRE(owner != nullptr);
            const int local_direction =
                weight.source_index - owner->outgoing_offset();
            const Eigen::Vector3d compiled_direction =
                owner->outgoing_sphere().get_quad_position(local_direction);
            REQUIRE(compiled_direction.dot(expected_direction) ==
                    Catch::Approx(1.0).margin(1.0e-12));
        }
        REQUIRE(weight_sum == Catch::Approx(1.0).margin(1.0e-12));
    }
} // namespace

TEST_CASE("Successive-orders 1D geometry compiles midpoint source and LOS "
          "topology",
          "[successive_orders][geometry]") {
    sasktran2::Geometry1D geometry(0.4, 0.0, 6372000.0, altitude_grid(),
                                   sasktran2::grids::interpolation::linear,
                                   sasktran2::geometrytype::spherical);
    sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    const auto los = make_los_geometry(geometry, raytracer);

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 14;
    settings.num_threads = 1;
    sasktran2::successive_orders::SourceGeometry1D source_geometry(raytracer,
                                                                   geometry);
    source_geometry.initialize(los, settings);

    REQUIRE(source_geometry.source_altitudes_m() ==
            std::vector<double>{500.0, 2000.0});
    REQUIRE(source_geometry.num_interior_points() == 2);
    REQUIRE(source_geometry.num_ground_points() == 1);
    REQUIRE(source_geometry.num_points() == 3);

    const auto& ground_point =
        source_geometry.source_point(source_geometry.num_interior_points());
    const std::array<Eigen::Vector3d, 3> ground_directions{
        Eigen::Vector3d::UnitZ(), Eigen::Vector3d::UnitX(),
        -Eigen::Vector3d::UnitZ()};
    for (const auto& direction : ground_directions) {
        std::vector<std::pair<int, double>> weights;
        int count = 0;
        ground_point.outgoing_sphere().interpolate(direction, weights, count);
        REQUIRE(count == 3);
        REQUIRE(weights.size() == 3);
        double total_weight = 0.0;
        for (const auto& [index, weight] : weights) {
            REQUIRE(index >= 0);
            REQUIRE(index < ground_point.outgoing_sphere().num_points());
            REQUIRE(std::isfinite(weight));
            REQUIRE(weight >= 0.0);
            REQUIRE(ground_point.outgoing_sphere().get_quad_position(index).dot(
                        ground_point.location().position.normalized()) > 0.0);
            total_weight += weight;
        }
        REQUIRE(total_weight == Catch::Approx(1.0).margin(1.0e-13));
    }

    REQUIRE(source_geometry.incoming_point_offsets().size() == 4);
    REQUIRE(source_geometry.outgoing_point_offsets().size() == 4);
    REQUIRE(source_geometry.incoming_rays().size() ==
            static_cast<std::size_t>(source_geometry.total_num_incoming()));
    for (std::size_t ray = 0; ray < source_geometry.incoming_rays().size();
         ++ray) {
        REQUIRE(source_geometry.incoming_interpolation()[ray].layers.size() ==
                source_geometry.incoming_rays()[ray].layers.size());
        REQUIRE(source_geometry.incoming_interpolation()[ray].ground_is_hit() ==
                source_geometry.incoming_rays()[ray].ground_is_hit);
    }
    for (const auto& point : source_geometry.source_points()) {
        require_sorted(point.atmosphere_weights());
    }

    require_compiled_topology(source_geometry.incoming_interpolation(),
                              source_geometry.transport_row_offsets(),
                              source_geometry.transport_column_indices(),
                              source_geometry.total_num_outgoing());
    require_compiled_topology(source_geometry.los_interpolation(),
                              source_geometry.los_transport_row_offsets(),
                              source_geometry.los_transport_column_indices(),
                              source_geometry.total_num_outgoing());
}

TEST_CASE("Successive-orders reduced-horizon grid supports arbitrary practical "
          "incoming counts",
          "[successive_orders][geometry][reduced_horizon]") {
    sasktran2::Geometry1D geometry(0.4, 0.0, 6372000.0, altitude_grid(),
                                   sasktran2::grids::interpolation::linear,
                                   sasktran2::geometrytype::spherical);
    sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    const auto los = make_los_geometry(geometry, raytracer);

    struct Budget {
        int points;
        int ground_rings;
        int space_rings;
    };
    constexpr std::array<Budget, 4> budgets{Budget{6, 1, 1}, Budget{37, 4, 4},
                                            Budget{110, 7, 6},
                                            Budget{257, 10, 10}};
    const double surface_radius = geometry.coordinates().earth_radius();
    const auto unique_count = [](std::vector<double> values) {
        std::sort(values.begin(), values.end());
        const auto end = std::unique(
            values.begin(), values.end(), [](double left, double right) {
                return std::abs(left - right) < 1.0e-12;
            });
        return static_cast<int>(std::distance(values.begin(), end));
    };

    for (const auto& budget : budgets) {
        DYNAMIC_SECTION(budget.points << " incoming nodes") {
            sasktran2::successive_orders::SourceGeometrySettings settings;
            settings.num_incoming = budget.points;
            settings.num_outgoing = 14;
            settings.num_threads = 1;
            settings.use_reduced_horizon_quadrature = true;
            sasktran2::successive_orders::SourceGeometry1D source_geometry(
                raytracer, geometry);
            source_geometry.initialize(los, settings);

            REQUIRE(source_geometry.num_interior_points() == 2);
            REQUIRE(&source_geometry.source_point(0).incoming_sphere() !=
                    &source_geometry.source_point(1).incoming_sphere());
            REQUIRE(&source_geometry.source_point(0).outgoing_sphere() ==
                    &source_geometry.source_point(1).outgoing_sphere());

            for (int direction = 0;
                 direction < source_geometry.source_point(0).num_outgoing();
                 ++direction) {
                const auto position = source_geometry.source_point(0)
                                          .outgoing_sphere()
                                          .get_quad_position(direction);
                REQUIRE(position.head<2>().norm() > 1.0e-8);
            }

            for (int point_index = 0;
                 point_index < source_geometry.num_interior_points();
                 ++point_index) {
                const auto& point = source_geometry.source_point(point_index);
                REQUIRE(point.num_incoming() == budget.points);
                const Eigen::Vector3d vertical =
                    point.location().position.normalized();
                const double radius = point.location().position.norm();
                const double horizon_mu = -std::sqrt(
                    1.0 - surface_radius * surface_radius / (radius * radius));
                std::vector<double> ground_mu;
                std::vector<double> space_mu;
                double weight_sum = 0.0;
                for (int direction = 0; direction < point.num_incoming();
                     ++direction) {
                    const double mu = point.incoming_sphere()
                                          .get_quad_position(direction)
                                          .dot(vertical);
                    const double weight =
                        point.incoming_sphere().quadrature_weight(direction);
                    REQUIRE(weight > 0.0);
                    weight_sum += weight;
                    (mu < horizon_mu ? ground_mu : space_mu).push_back(mu);
                    REQUIRE(
                        source_geometry
                            .incoming_interpolation()[point.incoming_offset() +
                                                      direction]
                            .ground_is_hit() == (mu < horizon_mu));
                }
                REQUIRE(unique_count(std::move(ground_mu)) ==
                        budget.ground_rings);
                REQUIRE(unique_count(std::move(space_mu)) ==
                        budget.space_rings);
                REQUIRE(weight_sum == Catch::Approx(1.0).margin(2.0e-14));
                const int exact_degree =
                    2 * std::min(budget.ground_rings, budget.space_rings) - 1;
                for (int degree = 0; degree <= exact_degree; ++degree) {
                    double actual_moment = 0.0;
                    for (int direction = 0; direction < point.num_incoming();
                         ++direction) {
                        const double mu = point.incoming_sphere()
                                              .get_quad_position(direction)
                                              .dot(vertical);
                        actual_moment +=
                            point.incoming_sphere().quadrature_weight(
                                direction) *
                            std::pow(mu, degree);
                    }
                    const double expected_moment =
                        degree % 2 == 0 ? 1.0 / (degree + 1.0) : 0.0;
                    REQUIRE(actual_moment ==
                            Catch::Approx(expected_moment).margin(2.0e-13));
                }
            }

            const auto& ground = source_geometry.source_point(
                source_geometry.num_interior_points());
            double ground_weight_sum = 0.0;
            for (int direction = 0; direction < ground.num_incoming();
                 ++direction) {
                REQUIRE(
                    ground.incoming_sphere().get_quad_position(direction).dot(
                        ground.location().position) > 0.0);
                const double weight =
                    ground.incoming_sphere().quadrature_weight(direction);
                REQUIRE(weight > 0.0);
                ground_weight_sum += weight;
            }
            REQUIRE(ground_weight_sum == Catch::Approx(0.5).margin(2.0e-14));
        }
    }

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.use_reduced_horizon_quadrature = true;
    settings.num_incoming = 5;
    REQUIRE_NOTHROW(settings.validate());
}

TEST_CASE("Successive-orders reduced-horizon grid uses the Geometry1D lower "
          "boundary",
          "[successive_orders][geometry][reduced_horizon]") {
    Eigen::VectorXd altitudes(3);
    altitudes << 5000.0, 7000.0, 10000.0;
    sasktran2::Geometry1D geometry(0.4, 0.0, 6372000.0, std::move(altitudes),
                                   sasktran2::grids::interpolation::linear,
                                   sasktran2::geometrytype::spherical);
    sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 37;
    settings.num_outgoing = 14;
    settings.num_threads = 1;
    settings.use_reduced_horizon_quadrature = true;
    sasktran2::successive_orders::SourceGeometry1D source_geometry(raytracer,
                                                                   geometry);
    sasktran2::viewinggeometry::InternalViewingGeometry viewing;
    source_geometry.initialize(viewing, settings);

    const double surface_radius = geometry.coordinates().earth_radius() +
                                  geometry.altitude_grid().grid()[0];
    REQUIRE(source_geometry.settings().use_reduced_horizon_quadrature);
    for (int ground_index = 0;
         ground_index < source_geometry.num_ground_points(); ++ground_index) {
        const auto& ground = source_geometry.source_point(
            source_geometry.num_interior_points() + ground_index);
        REQUIRE(ground.location().position.norm() ==
                Catch::Approx(surface_radius + 0.01).margin(1.0e-8));
    }
}

TEST_CASE("Successive-orders reduced-horizon incoming grid rotates with its "
          "coordinate system including at the solar poles",
          "[successive_orders][geometry][reduced_horizon]") {
    constexpr std::array<double, 3> solar_cosines{0.4, 1.0, -1.0};
    for (const double cos_sza : solar_cosines) {
        DYNAMIC_SECTION("cos_sza=" << cos_sza) {
            Eigen::VectorXd base_altitudes(3);
            base_altitudes << 0.0, 1000.0, 3000.0;
            sasktran2::Geometry1D base_geometry(
                cos_sza, 0.2, 6372000.0, std::move(base_altitudes),
                sasktran2::grids::interpolation::linear,
                sasktran2::geometrytype::spherical);

            const Eigen::Matrix3d rotation =
                (Eigen::AngleAxisd(0.7, Eigen::Vector3d::UnitY()) *
                 Eigen::AngleAxisd(-0.3, Eigen::Vector3d::UnitX()))
                    .toRotationMatrix();
            sasktran2::Coordinates rotated_coordinates(
                rotation * base_geometry.coordinates().reference_z(),
                rotation * base_geometry.coordinates().reference_x(),
                rotation * base_geometry.coordinates().sun_unit(),
                base_geometry.coordinates().earth_radius(),
                sasktran2::geometrytype::spherical);
            Eigen::VectorXd rotated_altitudes(3);
            rotated_altitudes << 0.0, 1000.0, 3000.0;
            sasktran2::grids::AltitudeGrid rotated_altitude_grid(
                std::move(rotated_altitudes),
                sasktran2::grids::gridspacing::automatic,
                sasktran2::grids::outofbounds::extend,
                sasktran2::grids::interpolation::linear);
            sasktran2::Geometry1D rotated_geometry(
                std::move(rotated_coordinates),
                std::move(rotated_altitude_grid));

            sasktran2::raytracing::SphericalShellRayTracer base_raytracer(
                base_geometry);
            sasktran2::raytracing::SphericalShellRayTracer rotated_raytracer(
                rotated_geometry);
            sasktran2::successive_orders::SourceGeometrySettings settings;
            settings.num_incoming = 37;
            settings.num_outgoing = 14;
            settings.num_threads = 1;
            settings.use_reduced_horizon_quadrature = true;
            sasktran2::successive_orders::SourceGeometry1D base_source(
                base_raytracer, base_geometry);
            sasktran2::successive_orders::SourceGeometry1D rotated_source(
                rotated_raytracer, rotated_geometry);
            sasktran2::viewinggeometry::InternalViewingGeometry viewing;
            base_source.initialize(viewing, settings);
            rotated_source.initialize(viewing, settings);

            REQUIRE(base_source.num_points() == rotated_source.num_points());
            for (int point_index = 0; point_index < base_source.num_points();
                 ++point_index) {
                const auto& base_point = base_source.source_point(point_index);
                const auto& rotated_point =
                    rotated_source.source_point(point_index);
                REQUIRE(base_point.num_incoming() ==
                        rotated_point.num_incoming());
                REQUIRE((rotation * base_point.location().position -
                         rotated_point.location().position)
                            .norm() < 1.0e-8);
                for (int direction = 0; direction < base_point.num_incoming();
                     ++direction) {
                    REQUIRE((rotation *
                                 base_point.incoming_sphere().get_quad_position(
                                     direction) -
                             rotated_point.incoming_sphere().get_quad_position(
                                 direction))
                                .norm() < 1.0e-11);
                    REQUIRE(base_point.incoming_sphere().quadrature_weight(
                                direction) ==
                            Catch::Approx(rotated_point.incoming_sphere()
                                              .quadrature_weight(direction))
                                .margin(2.0e-15));
                }
            }
        }
    }
}

TEST_CASE("Successive-orders reduced-horizon preference falls back for "
          "diffuse refraction",
          "[successive_orders][geometry][reduced_horizon]") {
    sasktran2::Geometry1D geometry(0.4, 0.0, 6372000.0, altitude_grid(),
                                   sasktran2::grids::interpolation::linear,
                                   sasktran2::geometrytype::spherical);
    sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    const auto los = make_los_geometry(geometry, raytracer);

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.num_threads = 1;
    settings.include_refraction = true;
    settings.use_reduced_horizon_quadrature = true;
    sasktran2::successive_orders::SourceGeometry1D source_geometry(raytracer,
                                                                   geometry);
    source_geometry.initialize(los, settings);

    REQUIRE_FALSE(source_geometry.settings().use_reduced_horizon_quadrature);
}

#ifdef SKTRAN_RUST_SUPPORT
TEST_CASE("Successive-orders 2D geometry uses an independent horizontal "
          "source grid",
          "[successive_orders][geometry][geometry2d]") {
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, altitude_grid(),
                                   horizontal_angle_grid(),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    const auto los = make_los_geometry(geometry, raytracer);

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.num_sza = 5;
    settings.num_threads = 1;
    sasktran2::successive_orders::SourceGeometry1D source_geometry(raytracer,
                                                                   geometry);
    source_geometry.initialize(los, settings);

    REQUIRE(source_geometry.source_altitudes_m() ==
            std::vector<double>{500.0, 2000.0});
    const std::array<double, 5> expected_horizontal_angles = {-0.4, -0.2, 0.0,
                                                              0.2, 0.4};
    REQUIRE(source_geometry.source_horizontal_angles_rad().size() ==
            expected_horizontal_angles.size());
    for (std::size_t index = 0; index < expected_horizontal_angles.size();
         ++index) {
        REQUIRE(
            source_geometry.source_horizontal_angles_rad()[index] ==
            Catch::Approx(expected_horizontal_angles[index]).margin(1.0e-13));
    }
    REQUIRE(source_geometry.num_interior_points() == 10);
    REQUIRE(source_geometry.num_ground_points() == 5);
    REQUIRE(source_geometry.num_points() == 15);

    for (int horizontal = 0; horizontal < 5; ++horizontal) {
        for (int altitude = 0; altitude < 2; ++altitude) {
            const int point_index = altitude + 2 * horizontal;
            const auto& point = source_geometry.source_point(point_index);
            REQUIRE(
                geometry.altitude_at(point.location()) ==
                Catch::Approx(source_geometry.source_altitudes_m()[altitude])
                    .margin(1.0e-8));
            REQUIRE(
                geometry.horizontal_angle_at(point.location()) ==
                Catch::Approx(
                    source_geometry.source_horizontal_angles_rad()[horizontal])
                    .margin(1.0e-13));

            std::vector<double> atmosphere_weights(geometry.size(), 0.0);
            for (const auto& weight : point.atmosphere_weights()) {
                atmosphere_weights[weight.index] += weight.weight;
            }
            REQUIRE(std::accumulate(atmosphere_weights.begin(),
                                    atmosphere_weights.end(),
                                    0.0) == Catch::Approx(1.0).margin(1.0e-13));
        }
    }

    // The source column at -0.2 radians lies halfway between atmosphere
    // columns. Combined with the midpoint altitude, it samples four native
    // atmosphere locations with equal weight.
    std::vector<double> midpoint_weights(geometry.size(), 0.0);
    for (const auto& weight :
         source_geometry.source_point(2).atmosphere_weights()) {
        midpoint_weights[weight.index] += weight.weight;
    }
    REQUIRE(midpoint_weights[geometry.location_index(0, 0)] ==
            Catch::Approx(0.25).margin(1.0e-13));
    REQUIRE(midpoint_weights[geometry.location_index(1, 0)] ==
            Catch::Approx(0.25).margin(1.0e-13));
    REQUIRE(midpoint_weights[geometry.location_index(0, 1)] ==
            Catch::Approx(0.25).margin(1.0e-13));
    REQUIRE(midpoint_weights[geometry.location_index(1, 1)] ==
            Catch::Approx(0.25).margin(1.0e-13));

    require_compiled_topology(source_geometry.incoming_interpolation(),
                              source_geometry.transport_row_offsets(),
                              source_geometry.transport_column_indices(),
                              source_geometry.total_num_outgoing());
    require_compiled_topology(source_geometry.los_interpolation(),
                              source_geometry.los_transport_row_offsets(),
                              source_geometry.los_transport_column_indices(),
                              source_geometry.total_num_outgoing());

    settings.include_refraction = true;
    sasktran2::successive_orders::SourceGeometry1D refracted(raytracer,
                                                             geometry);
    REQUIRE_THROWS_WITH(
        refracted.initialize(los, settings),
        "Geometry2D successive orders does not support diffuse-ray refraction");
}

TEST_CASE("Successive-orders reduced-horizon grid follows the Geometry2D "
          "ground boundary",
          "[successive_orders][geometry][geometry2d][reduced_horizon]") {
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, altitude_grid(),
                                   horizontal_angle_grid(),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    const auto los = make_los_geometry(geometry, raytracer);

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 37;
    settings.num_outgoing = 14;
    settings.num_sza = 3;
    settings.num_threads = 1;
    settings.use_reduced_horizon_quadrature = true;
    sasktran2::successive_orders::SourceGeometry1D source_geometry(raytracer,
                                                                   geometry);
    source_geometry.initialize(los, settings);

    const double surface_radius = geometry.coordinates().earth_radius() +
                                  geometry.altitude_grid().grid()[0];
    REQUIRE(source_geometry.num_interior_points() == 6);
    for (int point_index = 0;
         point_index < source_geometry.num_interior_points(); ++point_index) {
        const auto& point = source_geometry.source_point(point_index);
        REQUIRE(point.num_incoming() == settings.num_incoming);
        REQUIRE(point.num_outgoing() == settings.num_outgoing);
        const Eigen::Vector3d vertical = point.location().position.normalized();
        const double radius = point.location().position.norm();
        const double horizon_mu = -std::sqrt(
            1.0 - surface_radius * surface_radius / (radius * radius));
        double weight_sum = 0.0;
        for (int direction = 0; direction < point.num_incoming(); ++direction) {
            const double mu =
                point.incoming_sphere().get_quad_position(direction).dot(
                    vertical);
            weight_sum += point.incoming_sphere().quadrature_weight(direction);
            REQUIRE(source_geometry
                        .incoming_interpolation()[point.incoming_offset() +
                                                  direction]
                        .ground_is_hit() == (mu < horizon_mu));
        }
        REQUIRE(weight_sum == Catch::Approx(1.0).margin(2.0e-14));
    }

    require_compiled_topology(source_geometry.incoming_interpolation(),
                              source_geometry.transport_row_offsets(),
                              source_geometry.transport_column_indices(),
                              source_geometry.total_num_outgoing());
}

TEST_CASE("Successive-orders 2D geometry accepts explicit horizontal source "
          "angles",
          "[successive_orders][geometry][geometry2d]") {
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, altitude_grid(),
                                   horizontal_angle_grid(),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    const auto los = make_los_geometry(geometry, raytracer);

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.num_sza = 99;
    settings.horizontal_angle_grid_radians = {-0.35, -0.05, 0.12, 0.38};
    sasktran2::successive_orders::SourceGeometry1D source_geometry(raytracer,
                                                                   geometry);
    source_geometry.initialize(los, settings);

    REQUIRE(source_geometry.source_horizontal_angles_rad() ==
            settings.horizontal_angle_grid_radians);
    REQUIRE(source_geometry.num_interior_points() == 8);
    REQUIRE(source_geometry.num_ground_points() == 4);

    settings.horizontal_angle_grid_radians = {-0.3, -0.3};
    sasktran2::successive_orders::SourceGeometry1D unordered(raytracer,
                                                             geometry);
    REQUIRE_THROWS_WITH(
        unordered.initialize(los, settings),
        "Successive-orders source horizontal angles must be finite and "
        "strictly increasing");

    settings.horizontal_angle_grid_radians = {-0.5, 0.0};
    sasktran2::successive_orders::SourceGeometry1D outside(raytracer, geometry);
    REQUIRE_THROWS_WITH(
        outside.initialize(los, settings),
        "Successive-orders source horizontal angles must lie inside the "
        "Geometry2D horizontal angle range");
}

TEST_CASE("Successive-orders 2D solar table resolves source-ray endpoint OD",
          "[successive_orders][geometry][geometry2d][solartable]") {
    constexpr int num_altitudes = 16;
    constexpr int num_horizontal = 8;
    Eigen::VectorXd altitudes =
        Eigen::VectorXd::LinSpaced(num_altitudes, 0.0, 80000.0);
    Eigen::VectorXd horizontal =
        Eigen::VectorXd::LinSpaced(num_horizontal, -0.4, 0.4);
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    const auto los = make_los_geometry(geometry, raytracer);

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 14;
    settings.num_outgoing = 6;
    settings.num_sza = 3;
    settings.num_threads = 1;
    settings.altitude_grid_m.resize(7);
    for (int index = 0; index < 7; ++index) {
        settings.altitude_grid_m[index] = (index + 0.5) * 80000.0 / 7.0;
    }
    sasktran2::successive_orders::SourceGeometry1D source_geometry(raytracer,
                                                                   geometry);
    source_geometry.initialize(los, settings);
    const auto& rays = source_geometry.incoming_rays();

    sasktran2::Config config;
    sasktran2::solartransmission::SolarTransmissionTable2D table(geometry,
                                                                 raytracer);
    table.initialize_config(config);
    table.initialize_geometry(rays);
    sasktran2::solartransmission::SolarTableInterpolation interpolation;
    std::vector<bool> table_ground_hit;
    table.generate_interpolation(rays, interpolation, table_ground_hit);

    sasktran2::solartransmission::SolarTransmissionExact exact(geometry,
                                                               raytracer);
    sasktran2::solartransmission::SolarGeometryMatrix exact_matrix;
    std::vector<bool> exact_ground_hit;
    exact.generate_geometry_matrix(rays, exact_matrix, exact_ground_hit);
    REQUIRE(table_ground_hit == exact_ground_hit);

    Eigen::VectorXd extinction(geometry.size());
    for (int horizontal_index = 0; horizontal_index < num_horizontal;
         ++horizontal_index) {
        for (int altitude_index = 0; altitude_index < num_altitudes;
             ++altitude_index) {
            const double altitude =
                geometry.altitude_grid().grid()[altitude_index];
            const double angle =
                geometry.horizontal_angle_grid()[horizontal_index];
            extinction[geometry.location_index(altitude_index,
                                               horizontal_index)] =
                1.5e-5 * std::exp(-altitude / 18000.0) *
                (1.0 + 0.2 * std::sin(2.0 * EIGEN_PI * angle / 0.8));
        }
    }
    Eigen::VectorXd table_nodes(table.table_size());
    Eigen::VectorXd table_od(interpolation.rows());
    Eigen::VectorXd exact_od(exact_matrix.rows());
    table.apply(extinction, table_nodes);
    interpolation.apply(table_nodes, table_od);
    exact_matrix.multiply(extinction, exact_od);

    double maximum_absolute = 0.0;
    double maximum_relative = 0.0;
    double mean_absolute = 0.0;
    int active = 0;
    for (Eigen::Index row = 0; row < exact_od.size(); ++row) {
        if (exact_ground_hit[row]) {
            continue;
        }
        const double absolute = std::abs(table_od[row] - exact_od[row]);
        maximum_absolute = std::max(maximum_absolute, absolute);
        if (std::abs(exact_od[row]) > 1.0e-10) {
            maximum_relative =
                std::max(maximum_relative, absolute / std::abs(exact_od[row]));
        }
        mean_absolute += absolute;
        ++active;
    }
    mean_absolute /= active;
    CAPTURE(active, maximum_absolute, maximum_relative, mean_absolute);
    REQUIRE(maximum_relative < 0.06);
    REQUIRE(mean_absolute < 0.006);
}
#endif

TEST_CASE("Successive-orders default source grid preserves nonuniform midpoint "
          "interpolation",
          "[successive_orders][geometry]") {
    Eigen::VectorXd altitudes(4);
    altitudes << 0.0, 1000.0, 3000.0, 6000.0;
    sasktran2::Geometry1D geometry(0.4, 0.0, 6372000.0, std::move(altitudes),
                                   sasktran2::grids::interpolation::linear,
                                   sasktran2::geometrytype::planeparallel);
    sasktran2::raytracing::PlaneParallelRayTracer raytracer(geometry);

    sasktran2::viewinggeometry::InternalViewingGeometry los;
    los.traced_rays.resize(1);
    auto& ray = los.traced_rays.front();
    const Eigen::Vector3d location =
        geometry.coordinates().reference_point(3000.0);
    const Eigen::Vector3d direction = location.normalized();
    ray.observer_and_look.observer.position = location;
    ray.observer_and_look.look_away = direction;
    ray.layers.resize(1);
    ray.layers.front().entrance.position = location;
    ray.layers.front().exit.position = location;
    ray.layers.front().average_look_away = direction;
    ray.layers.front().cos_sza_entrance = 0.4;
    ray.layers.front().cos_sza_exit = 0.4;

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.use_reduced_horizon_quadrature = true;
    sasktran2::successive_orders::SourceGeometry1D source_geometry(raytracer,
                                                                   geometry);
    source_geometry.initialize(los, settings);

    REQUIRE_FALSE(source_geometry.settings().use_reduced_horizon_quadrature);

    REQUIRE(source_geometry.source_altitudes_m() ==
            std::vector<double>{500.0, 2000.0, 4500.0});
    std::vector<double> location_weights(3, 0.0);
    for (const auto& weight :
         source_geometry.los_interpolation().front().source_for_layer(0)) {
        int owner = -1;
        for (int point_index = 0;
             point_index < source_geometry.num_interior_points();
             ++point_index) {
            const auto& point = source_geometry.source_point(point_index);
            if (weight.source_index >= point.outgoing_offset() &&
                weight.source_index <
                    point.outgoing_offset() + point.num_outgoing()) {
                owner = point_index;
                break;
            }
        }
        REQUIRE(owner >= 0);
        location_weights[owner] += weight.weight;
    }

    REQUIRE(location_weights[0] == Catch::Approx(0.0).margin(1.0e-14));
    REQUIRE(location_weights[1] == Catch::Approx(0.6).margin(1.0e-13));
    REQUIRE(location_weights[2] == Catch::Approx(0.4).margin(1.0e-13));
}

TEST_CASE("Successive-orders 1D geometry accepts a valid explicit altitude "
          "grid and rejects invalid grids",
          "[successive_orders][geometry]") {
    sasktran2::Geometry1D geometry(0.4, 0.0, 6372000.0, altitude_grid(),
                                   sasktran2::grids::interpolation::linear,
                                   sasktran2::geometrytype::spherical);
    sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    const auto los = make_los_geometry(geometry, raytracer);

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.altitude_grid_m = {250.0, 2750.0};
    sasktran2::successive_orders::SourceGeometry1D valid(raytracer, geometry);
    valid.initialize(los, settings);
    REQUIRE(valid.source_altitudes_m() == settings.altitude_grid_m);

    settings.altitude_grid_m = {500.0, 400.0};
    sasktran2::successive_orders::SourceGeometry1D unordered(raytracer,
                                                             geometry);
    REQUIRE_THROWS_AS(unordered.initialize(los, settings),
                      std::invalid_argument);

    settings.altitude_grid_m = {-1.0, 500.0};
    sasktran2::successive_orders::SourceGeometry1D outside(raytracer, geometry);
    REQUIRE_THROWS_AS(outside.initialize(los, settings), std::invalid_argument);

    const std::string boundary_error =
        "Successive-orders source altitudes must lie strictly inside the "
        "atmosphere altitude range";
    settings.altitude_grid_m = {0.0, 500.0};
    sasktran2::successive_orders::SourceGeometry1D lower_boundary(raytracer,
                                                                  geometry);
    REQUIRE_THROWS_WITH(lower_boundary.initialize(los, settings),
                        boundary_error);

    settings.altitude_grid_m = {500.0, 3000.0};
    sasktran2::successive_orders::SourceGeometry1D upper_boundary(raytracer,
                                                                  geometry);
    REQUIRE_THROWS_WITH(upper_boundary.initialize(los, settings),
                        boundary_error);
}

TEST_CASE("Successive-orders pseudospherical geometry avoids duplicate SZA "
          "source columns",
          "[successive_orders][geometry]") {
    sasktran2::Geometry1D geometry(0.4, 0.0, 6372000.0, altitude_grid(),
                                   sasktran2::grids::interpolation::linear,
                                   sasktran2::geometrytype::pseudospherical);
    sasktran2::raytracing::PlaneParallelRayTracer raytracer(geometry);
    auto los = make_los_geometry(geometry, raytracer);
    // Even if externally supplied ray metadata spans an SZA range, a
    // pseudospherical Geometry1D source has only one physical column.
    los.traced_rays.front().layers.front().cos_sza_entrance = 0.2;
    los.traced_rays.front().layers.front().cos_sza_exit = 0.6;

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.num_sza = 4;
    sasktran2::successive_orders::SourceGeometry1D source_geometry(raytracer,
                                                                   geometry);
    source_geometry.initialize(los, settings);

    REQUIRE(source_geometry.source_cos_sza() == std::vector<double>{0.4});
    REQUIRE(source_geometry.num_interior_points() == 2);
    REQUIRE(source_geometry.num_ground_points() == 1);
}

TEST_CASE("Successive-orders 1D rejects a Geometry2D horizontal source grid",
          "[successive_orders][geometry]") {
    sasktran2::Geometry1D geometry(0.4, 0.0, 6372000.0, altitude_grid(),
                                   sasktran2::grids::interpolation::linear,
                                   sasktran2::geometrytype::spherical);
    sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    const auto los = make_los_geometry(geometry, raytracer);

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.horizontal_angle_grid_radians = {-0.1, 0.1};
    sasktran2::successive_orders::SourceGeometry1D source_geometry(raytracer,
                                                                   geometry);
    REQUIRE_THROWS_WITH(
        source_geometry.initialize(los, settings),
        "An explicit successive-orders horizontal-angle grid is supported "
        "only with Geometry2D");
}

TEST_CASE("Successive-orders geometry propagates incoming ray-tracing errors "
          "outside the parallel region",
          "[successive_orders][geometry]") {
    sasktran2::Geometry1D geometry(0.4, 0.0, 6372000.0, altitude_grid(),
                                   sasktran2::grids::interpolation::linear,
                                   sasktran2::geometrytype::spherical);
    sasktran2::raytracing::SphericalShellRayTracer valid_raytracer(geometry);
    const auto los = make_los_geometry(geometry, valid_raytracer);
    const ThrowingRayTracer throwing_raytracer;

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.num_threads = 2;
    sasktran2::successive_orders::SourceGeometry1D source_geometry(
        throwing_raytracer, geometry);

    REQUIRE_THROWS_WITH(source_geometry.initialize(los, settings),
                        "deliberate incoming ray-tracing failure");
}

TEST_CASE("Successive-orders geometry propagates LOS interpolation errors "
          "outside the parallel region",
          "[successive_orders][geometry]") {
    sasktran2::Geometry1D geometry(0.4, 0.0, 6372000.0, altitude_grid(),
                                   sasktran2::grids::interpolation::linear,
                                   sasktran2::geometrytype::spherical);
    sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    auto los = make_los_geometry(geometry, raytracer);
    los.traced_rays.front().layers.front().entrance.position.x() =
        std::numeric_limits<double>::quiet_NaN();

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.num_threads = 2;
    sasktran2::successive_orders::SourceGeometry1D source_geometry(raytracer,
                                                                   geometry);

    REQUIRE_THROWS_WITH(source_geometry.initialize(los, settings),
                        "Invalid input. Check log for more information");
}

TEST_CASE("Successive-orders spherical interpolation preserves exact axial "
          "directions",
          "[successive_orders][geometry]") {
    sasktran2::Geometry1D geometry(1.0, 0.0, 6372000.0, altitude_grid(),
                                   sasktran2::grids::interpolation::linear,
                                   sasktran2::geometrytype::spherical);
    sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    const auto los = make_exact_direction_los(geometry);

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.num_threads = 2;
    sasktran2::successive_orders::SourceGeometry1D source_geometry(raytracer,
                                                                   geometry);
    source_geometry.initialize(los, settings);

    require_los_direction(source_geometry, 0, Eigen::Vector3d::UnitZ());
    require_los_direction(source_geometry, 1, -Eigen::Vector3d::UnitZ());
    require_los_direction(source_geometry, 2, Eigen::Vector3d::UnitX());
    require_los_direction(source_geometry, 3, Eigen::Vector3d::UnitY());
}
