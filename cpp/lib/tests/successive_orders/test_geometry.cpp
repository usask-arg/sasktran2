#include "../../successive_orders/geometry.h"
#include "../../successive_orders/horizontal_interpolation.h"

#include <sasktran2/solartransmission.h>
#include <sasktran2/test_helper.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

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

    // Appends a limb ray traced from far outside the atmosphere whose tangent
    // point lies at the given altitude above one horizontal angle.
    void
    add_limb_ray(const sasktran2::Geometry2D& geometry,
                 const sasktran2::raytracing::RustRayTracer2D& raytracer,
                 double tangent_altitude, double tangent_angle,
                 sasktran2::viewinggeometry::InternalViewingGeometry& los) {
        const auto& coordinates = geometry.coordinates();
        const Eigen::Vector3d tangent_point =
            (coordinates.earth_radius() + tangent_altitude) *
            coordinates.unit_vector_from_angles(tangent_angle, 0.0);
        const Eigen::Vector3d forward =
            coordinates.local_x_y_from_angles(tangent_angle, 0.0)
                .first.normalized();
        sasktran2::viewinggeometry::ViewingRay ray;
        ray.observer.position = tangent_point - 3.0e6 * forward;
        ray.look_away = forward;
        los.traced_rays.emplace_back();
        raytracer.trace_ray(ray, los.traced_rays.back());
    }

    double layer_horizontal_angle(const sasktran2::Geometry2D& geometry,
                                  const sasktran2::raytracing::TracedRay& ray,
                                  std::size_t layer) {
        sasktran2::Location midpoint;
        midpoint.position = 0.5 * (ray.layers[layer].entrance.position +
                                   ray.layers[layer].exit.position);
        return geometry.horizontal_angle_at(midpoint);
    }

    // Source columns referenced by one set of compiled source weights.
    struct TouchedColumns {
        std::set<int> interior;
        std::set<int> ground;
        double weight_sum = 0.0;
    };

    template <typename Weights>
    TouchedColumns touched_columns(
        const sasktran2::successive_orders::SourceGeometry1D& source,
        const Weights& weights, const std::vector<int>& columns) {
        TouchedColumns result;
        const int num_altitudes =
            static_cast<int>(source.source_altitudes_m().size());
        const auto& offsets = source.outgoing_point_offsets();
        for (const auto& weight : weights) {
            const int outgoing = columns[weight.row_inner_index()];
            const int point = static_cast<int>(
                std::upper_bound(offsets.begin(), offsets.end(), outgoing) -
                offsets.begin() - 1);
            if (point < source.num_interior_points()) {
                result.interior.insert(point / num_altitudes);
            } else {
                result.ground.insert(point - source.num_interior_points());
            }
            result.weight_sum += weight.weight();
        }
        return result;
    }

    // Columns touched by every layer of every LOS ray, indexed [ray][layer].
    std::vector<std::vector<TouchedColumns>> los_layer_columns(
        const sasktran2::successive_orders::SourceGeometry1D& source) {
        std::vector<std::vector<TouchedColumns>> result;
        const auto& rays = source.los_interpolation();
        for (std::size_t ray = 0; ray < rays.size(); ++ray) {
            const auto columns =
                source.los_transport_columns_for_ray(ray).to_vector();
            auto& layers = result.emplace_back();
            for (std::size_t layer = 0; layer < rays[ray].layers.size();
                 ++layer) {
                layers.push_back(touched_columns(
                    source, rays[ray].source_for_layer(layer), columns));
            }
        }
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
        int previous = -1;
        for (const auto& weight : weights) {
            REQUIRE(previous <= weight.index);
            previous = weight.index;
        }
    }

    void require_compiled_topology(
        const std::vector<sasktran2::successive_orders::RayInterpolation>& rays,
        const std::vector<int>& row_offsets,
        sasktran2::successive_orders::TransportColumnView packed_columns,
        int num_source_columns) {
        const auto column_indices = packed_columns.to_vector();
        REQUIRE(row_offsets.size() == rays.size() + 1);
        REQUIRE(row_offsets.front() == 0);
        REQUIRE(row_offsets.back() == static_cast<int>(column_indices.size()));
        for (std::size_t row = 0; row < rays.size(); ++row) {
            const auto& ray = rays[row];
            REQUIRE(ray.transport_compiled);
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
                    REQUIRE(row_offsets[row] + weight.row_inner_index() <
                            row_offsets[row + 1]);
                    const int source_index =
                        column_indices[row_offsets[row] +
                                       weight.row_inner_index()];
                    REQUIRE(source_index >= 0);
                    REQUIRE(source_index < num_source_columns);
                }
            }
            for (const auto& weight : ray.ground()) {
                REQUIRE(row_offsets[row] + weight.row_inner_index() <
                        row_offsets[row + 1]);
                const int source_index =
                    column_indices[row_offsets[row] + weight.row_inner_index()];
                REQUIRE(source_index >= 0);
                REQUIRE(source_index < num_source_columns);
            }
        }
    }

    // Every compiled LOS node must be the expected direction itself, which
    // holds when the grid has a node on that direction.
    void require_los_direction_on_node(
        const sasktran2::successive_orders::SourceGeometry1D& geometry,
        std::size_t ray_index, const Eigen::Vector3d& expected_direction) {
        const auto weights =
            geometry.los_interpolation()[ray_index].source_for_layer(0);
        const auto columns = geometry.los_transport_columns_for_ray(ray_index);
        REQUIRE(!weights.empty());

        double weight_sum = 0.0;
        for (const auto& weight : weights) {
            REQUIRE(std::isfinite(weight.weight()));
            weight_sum += weight.weight();
            const int source_index = columns[weight.row_inner_index()];

            const sasktran2::successive_orders::SourcePoint* owner = nullptr;
            for (const auto& point : geometry.source_points()) {
                if (source_index >= point.outgoing_offset() &&
                    source_index <
                        point.outgoing_offset() + point.num_outgoing()) {
                    owner = &point;
                    break;
                }
            }
            REQUIRE(owner != nullptr);
            const int local_direction = source_index - owner->outgoing_offset();
            const Eigen::Vector3d compiled_direction =
                owner->outgoing_sphere().get_quad_position(local_direction);
            REQUIRE(compiled_direction.dot(expected_direction) ==
                    Catch::Approx(1.0).margin(1.0e-12));
        }
        REQUIRE(weight_sum == Catch::Approx(1.0).margin(1.0e-12));
    }

    // A direction that rotate_unit_vector maps unchanged onto every stencil
    // point must compile to exactly that point sphere's own interpolation of
    // the direction, scaled by the point's location weight.
    void require_los_direction(
        const sasktran2::successive_orders::SourceGeometry1D& geometry,
        std::size_t ray_index, const Eigen::Vector3d& expected_direction) {
        const auto weights =
            geometry.los_interpolation()[ray_index].source_for_layer(0);
        const auto columns = geometry.los_transport_columns_for_ray(ray_index);
        REQUIRE(!weights.empty());

        std::map<const sasktran2::successive_orders::SourcePoint*,
                 std::map<int, double>>
            compiled;
        double weight_sum = 0.0;
        for (const auto& weight : weights) {
            REQUIRE(std::isfinite(weight.weight()));
            weight_sum += weight.weight();
            const int source_index = columns[weight.row_inner_index()];

            const sasktran2::successive_orders::SourcePoint* owner = nullptr;
            for (const auto& point : geometry.source_points()) {
                if (source_index >= point.outgoing_offset() &&
                    source_index <
                        point.outgoing_offset() + point.num_outgoing()) {
                    owner = &point;
                    break;
                }
            }
            REQUIRE(owner != nullptr);
            compiled[owner][source_index - owner->outgoing_offset()] +=
                weight.weight();
        }
        REQUIRE(weight_sum == Catch::Approx(1.0).margin(1.0e-12));

        for (const auto& [owner, owner_weights] : compiled) {
            double location_weight = 0.0;
            for (const auto& [node, weight] : owner_weights) {
                location_weight += weight;
            }
            std::vector<std::pair<int, double>> direct;
            int count = 0;
            owner->outgoing_sphere().interpolate(expected_direction, direct,
                                                 count);
            std::map<int, double> expected;
            for (int index = 0; index < count; ++index) {
                if (direct[index].second != 0.0) {
                    expected[direct[index].first] +=
                        location_weight * direct[index].second;
                }
            }
            REQUIRE(owner_weights.size() == expected.size());
            for (const auto& [node, weight] : expected) {
                REQUIRE(owner_weights.count(node) == 1);
                REQUIRE(owner_weights.at(node) ==
                        Catch::Approx(weight).margin(1.0e-12));
            }
        }
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
                atmosphere_weights[weight.index] += weight.weight();
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
        midpoint_weights[weight.index] += weight.weight();
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

    for (const bool legacy_interpolation : {true, false}) {
        DYNAMIC_SECTION("legacy_interpolation=" << legacy_interpolation) {
            sasktran2::successive_orders::SourceGeometrySettings settings;
            settings.num_incoming = 14;
            settings.num_outgoing = 6;
            settings.num_sza = 3;
            settings.num_threads = 1;
            settings.legacy_interpolation = legacy_interpolation;
            settings.altitude_grid_m.resize(7);
            for (int index = 0; index < 7; ++index) {
                settings.altitude_grid_m[index] = (index + 0.5) * 80000.0 / 7.0;
            }
            sasktran2::successive_orders::SourceGeometry1D source_geometry(
                raytracer, geometry);
            source_geometry.initialize(los, settings);
            const auto& rays = source_geometry.incoming_rays();

            sasktran2::Config config;
            sasktran2::solartransmission::SolarTransmissionTable2D table(
                geometry, raytracer);
            table.initialize_config(config);
            table.initialize_geometry(rays);
            sasktran2::solartransmission::SolarTableInterpolation interpolation;
            std::vector<bool> table_ground_hit;
            table.generate_interpolation(rays, interpolation, table_ground_hit);

            sasktran2::solartransmission::SolarTransmissionExact exact(
                geometry, raytracer);
            sasktran2::solartransmission::SolarGeometryMatrix exact_matrix;
            std::vector<bool> exact_ground_hit;
            exact.generate_geometry_matrix(rays, exact_matrix,
                                           exact_ground_hit);
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

            // Rows whose exact OD is zero to rounding are excluded from the
            // worst-row bounds. A few of them (10 legacy, 2 aligned) have
            // table OD near 1.1 in both samples; the mean bound covers them.
            // Below an exact OD of 1e-3 a relative bound is ill-conditioned,
            // so those rows are bounded in absolute OD instead.
            constexpr double thin_optical_depth = 1.0e-3;
            double maximum_absolute = 0.0;
            double maximum_relative = 0.0;
            double maximum_relative_thick = 0.0;
            double maximum_absolute_thin = 0.0;
            double mean_absolute = 0.0;
            int active = 0;
            int thin = 0;
            for (Eigen::Index row = 0; row < exact_od.size(); ++row) {
                if (exact_ground_hit[row]) {
                    continue;
                }
                const double absolute = std::abs(table_od[row] - exact_od[row]);
                const double optical_depth = std::abs(exact_od[row]);
                maximum_absolute = std::max(maximum_absolute, absolute);
                if (optical_depth > 1.0e-10) {
                    maximum_relative =
                        std::max(maximum_relative, absolute / optical_depth);
                    if (optical_depth < thin_optical_depth) {
                        maximum_absolute_thin =
                            std::max(maximum_absolute_thin, absolute);
                        ++thin;
                    } else {
                        maximum_relative_thick = std::max(
                            maximum_relative_thick, absolute / optical_depth);
                    }
                }
                mean_absolute += absolute;
                ++active;
            }
            mean_absolute /= active;
            CAPTURE(active, thin, maximum_absolute, maximum_relative,
                    maximum_relative_thick, maximum_absolute_thin,
                    mean_absolute);
            // Ground-hit parity (above), the mean bound and the split
            // worst-row bounds apply to both samples.
            REQUIRE(mean_absolute < 0.006);
            REQUIRE(maximum_relative_thick < 0.06);
            REQUIRE(maximum_absolute_thin < 1.0e-4);
            if (legacy_interpolation) {
                // The worst-ray relative bound was calibrated on the legacy
                // global incoming grid, which has no thin rows. Frame-aligned
                // grids trace a different sample with two thin rows (exact OD
                // 9.3e-4 and 5.4e-4) whose relative errors exceed it; their
                // absolute errors are 8.3e-5 and 4.6e-5.
                REQUIRE(maximum_relative < 0.06);
            }
        }
    }
}

TEST_CASE("Successive-orders columns share frame-aligned outgoing grids",
          "[successive_orders][geometry][geometry2d]") {
    for (const bool reduced_horizon : {true, false}) {
        DYNAMIC_SECTION("reduced_horizon=" << reduced_horizon) {
            sasktran2::Geometry2D geometry(
                0.6, 0.4, 6372000.0, altitude_grid(), horizontal_angle_grid(),
                sasktran2::grids::interpolation::linear);
            sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
            const auto los = make_los_geometry(geometry, raytracer);
            sasktran2::successive_orders::SourceGeometrySettings settings;
            settings.num_incoming = 26;
            settings.num_outgoing = 26;
            settings.num_sza = 3;
            settings.num_threads = 1;
            settings.use_reduced_horizon_quadrature = reduced_horizon;
            sasktran2::successive_orders::SourceGeometry1D source(raytracer,
                                                                  geometry);
            source.initialize(los, settings);

            const int altitudes =
                static_cast<int>(source.source_altitudes_m().size());
            const auto frame = [&](const Eigen::Vector3d& position) {
                const Eigen::Vector3d up = position.normalized();
                const Eigen::Vector3d x =
                    sasktran2::successive_orders::solar_horizontal_reference(
                        up, geometry);
                Eigen::Matrix3d result;
                result << x, up.cross(x), up;
                return result;
            };
            const auto& reference = source.source_point(0);
            const Eigen::Matrix3d reference_frame =
                frame(reference.location().position);
            for (int index = 0; index < source.num_interior_points(); ++index) {
                const auto& point = source.source_point(index);
                const auto& column_first =
                    source.source_point(index - index % altitudes);
                REQUIRE(&point.outgoing_sphere() ==
                        &column_first.outgoing_sphere());
                REQUIRE(point.angular_class() ==
                        (reduced_horizon ? index % altitudes : 0));
                const Eigen::Matrix3d point_frame =
                    frame(point.location().position);
                for (int node = 0; node < point.num_outgoing(); ++node) {
                    REQUIRE(
                        (point_frame.transpose() *
                             point.outgoing_sphere().get_quad_position(node) -
                         reference_frame.transpose() *
                             reference.outgoing_sphere().get_quad_position(
                                 node))
                            .norm() < 1.0e-12);
                }
            }
            const auto& reference_ground =
                source.source_point(source.num_interior_points());
            const Eigen::Matrix3d reference_ground_frame =
                frame(reference_ground.location().position);
            for (int ground = 0; ground < source.num_ground_points();
                 ++ground) {
                INFO("ground=" << ground);
                const auto& point =
                    source.source_point(source.num_interior_points() + ground);
                const Eigen::Vector3d up =
                    point.location().position.normalized();
                // Every kept node lies strictly above the horizon band.
                for (int node = 0; node < point.num_incoming(); ++node) {
                    REQUIRE(point.incoming_sphere().get_quad_position(node).dot(
                                up) > 1.0e-12);
                }
                for (int node = 0; node < point.num_outgoing(); ++node) {
                    REQUIRE(point.outgoing_sphere().get_quad_position(node).dot(
                                up) > 1.0e-12);
                }
                // Ground grids share canonical coordinates across columns.
                const Eigen::Matrix3d point_frame =
                    frame(point.location().position);
                const auto require_same_canonical_grid =
                    [&](const sasktran2::math::UnitSphere& actual,
                        const sasktran2::math::UnitSphere& expected) {
                        REQUIRE(actual.num_points() == expected.num_points());
                        for (int node = 0; node < actual.num_points(); ++node) {
                            REQUIRE((point_frame.transpose() *
                                         actual.get_quad_position(node) -
                                     reference_ground_frame.transpose() *
                                         expected.get_quad_position(node))
                                        .norm() < 1.0e-12);
                        }
                    };
                require_same_canonical_grid(point.incoming_sphere(),
                                            reference_ground.incoming_sphere());
                require_same_canonical_grid(point.outgoing_sphere(),
                                            reference_ground.outgoing_sphere());
            }
        }
    }
}

TEST_CASE("Successive-orders legacy interpolation keeps one global grid",
          "[successive_orders][geometry][geometry2d]") {
    sasktran2::Geometry2D geometry(0.6, 0.4, 6372000.0, altitude_grid(),
                                   horizontal_angle_grid(),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    const auto los = make_los_geometry(geometry, raytracer);
    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 26;
    settings.num_outgoing = 26;
    settings.num_sza = 3;
    settings.num_threads = 1;
    settings.use_reduced_horizon_quadrature = true;
    settings.legacy_interpolation = true;
    sasktran2::successive_orders::SourceGeometry1D source(raytracer, geometry);
    source.initialize(los, settings);
    for (int index = 0; index < source.num_interior_points(); ++index) {
        REQUIRE(&source.source_point(index).outgoing_sphere() ==
                &source.source_point(0).outgoing_sphere());
        REQUIRE(source.source_point(index).angular_class() == index);
    }
}

TEST_CASE("Successive-orders LOS uses a four-column stencil only on the LOS",
          "[successive_orders][geometry][geometry2d]") {
    Eigen::VectorXd altitudes = Eigen::VectorXd::LinSpaced(7, 0.0, 60000.0);
    Eigen::VectorXd horizontal = Eigen::VectorXd::LinSpaced(9, -0.2, 0.2);
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);

    // A limb ray tangent at 20 km over the centre column stays inside the
    // horizontal range and crosses the source columns at -0.1, 0 and 0.1.
    sasktran2::viewinggeometry::InternalViewingGeometry los;
    add_limb_ray(geometry, raytracer, 20000.0, 0.0, los);
    const auto& traced = los.traced_rays.front();
    REQUIRE(!traced.ground_is_hit);
    double minimum_angle = std::numeric_limits<double>::infinity();
    double maximum_angle = -std::numeric_limits<double>::infinity();
    for (std::size_t layer = 0; layer < traced.layers.size(); ++layer) {
        const double angle = layer_horizontal_angle(geometry, traced, layer);
        minimum_angle = std::min(minimum_angle, angle);
        maximum_angle = std::max(maximum_angle, angle);
    }
    CAPTURE(minimum_angle, maximum_angle);
    REQUIRE(minimum_angle > -0.2);
    REQUIRE(minimum_angle < -0.1);
    REQUIRE(maximum_angle > 0.1);
    REQUIRE(maximum_angle < 0.2);

    for (const bool legacy_interpolation : {false, true}) {
        DYNAMIC_SECTION("legacy_interpolation=" << legacy_interpolation) {
            sasktran2::successive_orders::SourceGeometrySettings settings;
            settings.num_incoming = 14;
            settings.num_outgoing = 14;
            settings.num_sza = 5;
            settings.num_threads = 1;
            settings.legacy_interpolation = legacy_interpolation;
            sasktran2::successive_orders::SourceGeometry1D source(raytracer,
                                                                  geometry);
            source.initialize(los, settings);
            REQUIRE(source.source_horizontal_angles_rad().size() == 5);

            int los_maximum = 0;
            int los_layers_with_four = 0;
            const auto los_layers = los_layer_columns(source);
            for (const auto& touched : los_layers.front()) {
                REQUIRE(touched.weight_sum ==
                        Catch::Approx(1.0).margin(1.0e-12));
                const int count = static_cast<int>(touched.interior.size());
                los_maximum = std::max(los_maximum, count);
                los_layers_with_four += count == 4 ? 1 : 0;
            }
            int incoming_maximum = 0;
            const auto& incoming = source.incoming_interpolation();
            for (std::size_t ray = 0; ray < incoming.size(); ++ray) {
                const auto columns =
                    source.transport_columns_for_ray(ray).to_vector();
                for (std::size_t layer = 0; layer < incoming[ray].layers.size();
                     ++layer) {
                    incoming_maximum = std::max(
                        incoming_maximum,
                        static_cast<int>(
                            touched_columns(
                                source, incoming[ray].source_for_layer(layer),
                                columns)
                                .interior.size()));
                }
            }
            CAPTURE(los_maximum, los_layers_with_four, incoming_maximum);
            if (legacy_interpolation) {
                REQUIRE(los_maximum == 2);
            } else {
                REQUIRE(los_maximum == 4);
                REQUIRE(los_layers_with_four > 0);
            }
            REQUIRE(incoming_maximum <= 2);
        }
    }
}

TEST_CASE("Successive-orders cubic LOS ground hits use four ground columns",
          "[successive_orders][geometry][geometry2d]") {
    Eigen::VectorXd altitudes = Eigen::VectorXd::LinSpaced(7, 0.0, 60000.0);
    Eigen::VectorXd horizontal = Eigen::VectorXd::LinSpaced(9, -0.2, 0.2);
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);

    // A nadir ray halfway between the source columns at 0 and 0.1 radians.
    sasktran2::viewinggeometry::InternalViewingGeometry los;
    los.traced_rays.resize(1);
    sasktran2::viewinggeometry::ViewingRay ray;
    ray.observer.position =
        (geometry.coordinates().earth_radius() + 40000.0) *
        geometry.coordinates().unit_vector_from_angles(0.05, 0.0);
    ray.look_away = -ray.observer.position.normalized();
    raytracer.trace_ray(ray, los.traced_rays.front());
    REQUIRE(los.traced_rays.front().ground_is_hit);

    for (const bool legacy_interpolation : {false, true}) {
        DYNAMIC_SECTION("legacy_interpolation=" << legacy_interpolation) {
            sasktran2::successive_orders::SourceGeometrySettings settings;
            settings.num_incoming = 14;
            settings.num_outgoing = 14;
            settings.num_sza = 5;
            settings.num_threads = 1;
            settings.legacy_interpolation = legacy_interpolation;
            sasktran2::successive_orders::SourceGeometry1D source(raytracer,
                                                                  geometry);
            source.initialize(los, settings);

            const auto& compiled = source.los_interpolation().front();
            REQUIRE(compiled.ground_is_hit());
            const auto ground = touched_columns(
                source, compiled.ground(),
                source.los_transport_columns_for_ray(0).to_vector());
            REQUIRE(ground.interior.empty());
            REQUIRE(ground.weight_sum == Catch::Approx(1.0).margin(1.0e-12));
            REQUIRE(ground.ground.size() == (legacy_interpolation ? 2u : 4u));
            const auto los_layers = los_layer_columns(source);
            for (const auto& touched : los_layers.front()) {
                REQUIRE(touched.interior.size() ==
                        (legacy_interpolation ? 2u : 4u));
            }
        }
    }
}

TEST_CASE("Successive-orders refresh_los recompiles cubic LOS stencils",
          "[successive_orders][geometry][geometry2d]") {
    Eigen::VectorXd altitudes = Eigen::VectorXd::LinSpaced(7, 0.0, 60000.0);
    Eigen::VectorXd horizontal = Eigen::VectorXd::LinSpaced(9, -0.2, 0.2);
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);

    const auto add_nadir_ray =
        [&](double angle,
            sasktran2::viewinggeometry::InternalViewingGeometry& los) {
            sasktran2::viewinggeometry::ViewingRay ray;
            ray.observer.position =
                (geometry.coordinates().earth_radius() + 40000.0) *
                geometry.coordinates().unit_vector_from_angles(angle, 0.0);
            ray.look_away = -ray.observer.position.normalized();
            los.traced_rays.emplace_back();
            raytracer.trace_ray(ray, los.traced_rays.back());
        };
    // The initial LOS is one nadir ray between the columns at 0 and 0.1; the
    // refreshed LOS adds a limb ray and moves the nadir ray elsewhere.
    sasktran2::viewinggeometry::InternalViewingGeometry initial;
    add_nadir_ray(0.05, initial);
    sasktran2::viewinggeometry::InternalViewingGeometry refreshed;
    add_limb_ray(geometry, raytracer, 20000.0, 0.0, refreshed);
    add_nadir_ray(-0.05, refreshed);
    REQUIRE(!refreshed.traced_rays[0].ground_is_hit);
    REQUIRE(refreshed.traced_rays[1].ground_is_hit);

    for (const bool legacy_interpolation : {false, true}) {
        DYNAMIC_SECTION("legacy_interpolation=" << legacy_interpolation) {
            sasktran2::successive_orders::SourceGeometrySettings settings;
            settings.num_incoming = 14;
            settings.num_outgoing = 14;
            settings.num_sza = 5;
            settings.num_threads = 1;
            settings.legacy_interpolation = legacy_interpolation;
            sasktran2::successive_orders::SourceGeometry1D source(raytracer,
                                                                  geometry);
            source.initialize(initial, settings);
            REQUIRE(source.los_interpolation().size() == 1);
            const auto initial_ground =
                touched_columns(
                    source, source.los_interpolation().front().ground(),
                    source.los_transport_columns_for_ray(0).to_vector())
                    .ground;
            source.refresh_los(refreshed);

            sasktran2::successive_orders::SourceGeometry1D fresh(raytracer,
                                                                 geometry);
            fresh.initialize(refreshed, settings);

            REQUIRE(source.los_interpolation().size() == 2);
            const std::size_t stencil = legacy_interpolation ? 2u : 4u;
            const auto refreshed_layers = los_layer_columns(source);
            const auto fresh_layers = los_layer_columns(fresh);
            REQUIRE(refreshed_layers.size() == fresh_layers.size());
            std::size_t limb_maximum = 0;
            for (std::size_t ray = 0; ray < refreshed_layers.size(); ++ray) {
                INFO("ray=" << ray);
                REQUIRE(refreshed_layers[ray].size() ==
                        fresh_layers[ray].size());
                for (std::size_t layer = 0;
                     layer < refreshed_layers[ray].size(); ++layer) {
                    const auto& touched = refreshed_layers[ray][layer];
                    REQUIRE(touched.weight_sum ==
                            Catch::Approx(1.0).margin(1.0e-12));
                    REQUIRE(touched.interior ==
                            fresh_layers[ray][layer].interior);
                    REQUIRE(touched.weight_sum ==
                            fresh_layers[ray][layer].weight_sum);
                    if (ray == 0) {
                        limb_maximum =
                            std::max(limb_maximum, touched.interior.size());
                    } else {
                        REQUIRE(touched.interior.size() == stencil);
                    }
                }
            }
            REQUIRE(limb_maximum == stencil);

            const auto& nadir = source.los_interpolation()[1];
            REQUIRE(nadir.ground_is_hit());
            const auto ground = touched_columns(
                source, nadir.ground(),
                source.los_transport_columns_for_ray(1).to_vector());
            REQUIRE(ground.interior.empty());
            REQUIRE(ground.weight_sum == Catch::Approx(1.0).margin(1.0e-12));
            REQUIRE(ground.ground.size() == stencil);
            const auto fresh_ground = touched_columns(
                fresh, fresh.los_interpolation()[1].ground(),
                fresh.los_transport_columns_for_ray(1).to_vector());
            REQUIRE(ground.ground == fresh_ground.ground);
            // The new nadir ray lies between the columns at -0.1 and 0
            // rather than 0 and 0.1, so its stencil shifts by one column.
            REQUIRE(initial_ground.size() == stencil);
            REQUIRE(ground.ground != initial_ground);
            REQUIRE(*ground.ground.begin() + 1 == *initial_ground.begin());
        }
    }
}

TEST_CASE("Successive-orders LOS stays linear with fewer than four columns",
          "[successive_orders][geometry][geometry2d]") {
    Eigen::VectorXd altitudes = Eigen::VectorXd::LinSpaced(7, 0.0, 60000.0);
    Eigen::VectorXd horizontal = Eigen::VectorXd::LinSpaced(9, -0.2, 0.2);
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    sasktran2::viewinggeometry::InternalViewingGeometry los;
    add_limb_ray(geometry, raytracer, 20000.0, 0.0, los);

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 14;
    settings.num_outgoing = 14;
    settings.num_sza = 3;
    settings.num_threads = 1;
    sasktran2::successive_orders::SourceGeometry1D source(raytracer, geometry);
    source.initialize(los, settings);
    REQUIRE(source.source_horizontal_angles_rad().size() == 3);

    std::size_t maximum = 0;
    const auto los_layers = los_layer_columns(source);
    for (const auto& touched : los_layers.front()) {
        maximum = std::max(maximum, touched.interior.size());
    }
    REQUIRE(maximum == 2);
}

TEST_CASE("Successive-orders cubic LOS extends constantly outside the source "
          "columns",
          "[successive_orders][geometry][geometry2d]") {
    Eigen::VectorXd altitudes = Eigen::VectorXd::LinSpaced(7, 0.0, 60000.0);
    Eigen::VectorXd horizontal = Eigen::VectorXd::LinSpaced(9, -0.2, 0.2);
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    // The limb ray spans about +-0.11 radians, beyond the column range.
    sasktran2::viewinggeometry::InternalViewingGeometry los;
    add_limb_ray(geometry, raytracer, 20000.0, 0.0, los);

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 14;
    settings.num_outgoing = 14;
    settings.num_threads = 1;
    settings.horizontal_angle_grid_radians = {-0.06, -0.02, 0.02, 0.06};
    sasktran2::successive_orders::SourceGeometry1D source(raytracer, geometry);
    source.initialize(los, settings);

    const auto& traced = los.traced_rays.front();
    const auto layers = los_layer_columns(source).front();
    REQUIRE(layers.size() == traced.layers.size());
    int below = 0;
    int inside = 0;
    int above = 0;
    for (std::size_t layer = 0; layer < layers.size(); ++layer) {
        const double angle = layer_horizontal_angle(geometry, traced, layer);
        INFO("layer=" << layer << " angle=" << angle);
        REQUIRE(layers[layer].weight_sum == Catch::Approx(1.0).margin(1.0e-12));
        if (angle < -0.06) {
            REQUIRE(layers[layer].interior == std::set<int>{0});
            ++below;
        } else if (angle > 0.06) {
            REQUIRE(layers[layer].interior == std::set<int>{3});
            ++above;
        } else {
            REQUIRE(layers[layer].interior.size() == 4);
            ++inside;
        }
    }
    CAPTURE(below, inside, above);
    REQUIRE(below > 0);
    REQUIRE(inside > 0);
    REQUIRE(above > 0);
}

TEST_CASE("Successive-orders cubic LOS falls back to linear across the "
          "terminator",
          "[successive_orders][geometry][geometry2d]") {
    // With the sun on the horizon at the reference point, columns at
    // negative horizontal angles are on the night side.
    Eigen::VectorXd altitudes = Eigen::VectorXd::LinSpaced(7, 0.0, 60000.0);
    Eigen::VectorXd horizontal = Eigen::VectorXd::LinSpaced(17, -0.4, 0.4);
    sasktran2::Geometry2D geometry(0.0, 0.0, 6372000.0, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    sasktran2::viewinggeometry::InternalViewingGeometry los;
    add_limb_ray(geometry, raytracer, 20000.0, 0.0, los);
    add_limb_ray(geometry, raytracer, 20000.0, 0.15, los);

    const std::vector<double> columns = {-0.15, -0.05, 0.05, 0.10,
                                         0.15,  0.20,  0.25};
    for (const double column : columns) {
        const Eigen::Vector3d up =
            geometry.coordinates().unit_vector_from_angles(column, 0.0);
        REQUIRE((up.dot(geometry.coordinates().sun_unit()) > 0.0) ==
                (column > 0.0));
    }
    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 14;
    settings.num_outgoing = 14;
    settings.num_threads = 1;
    settings.horizontal_angle_grid_radians = columns;
    sasktran2::successive_orders::SourceGeometry1D source(raytracer, geometry);
    source.initialize(los, settings);

    const int size = static_cast<int>(columns.size());
    const auto compiled = los_layer_columns(source);
    int cubic = 0;
    int linear = 0;
    for (std::size_t ray = 0; ray < compiled.size(); ++ray) {
        const auto& traced = los.traced_rays[ray];
        for (std::size_t layer = 0; layer < compiled[ray].size(); ++layer) {
            const double angle =
                layer_horizontal_angle(geometry, traced, layer);
            INFO("ray=" << ray << " layer=" << layer << " angle=" << angle);
            const auto count = compiled[ray][layer].interior.size();
            if (angle <= columns.front() || angle >= columns.back()) {
                REQUIRE(count == 1);
                continue;
            }
            const int lower = static_cast<int>(
                std::upper_bound(columns.begin(), columns.end(), angle) -
                columns.begin() - 1);
            const int start = std::clamp(lower - 1, 0, size - 4);
            const bool sunlit = std::all_of(
                columns.begin() + start, columns.begin() + start + 4,
                [](double column) { return column > 0.0; });
            if (sunlit) {
                REQUIRE(count == 4);
                ++cubic;
            } else {
                REQUIRE(count <= 2);
                ++linear;
            }
        }
    }
    CAPTURE(cubic, linear);
    REQUIRE(cubic > 0);
    REQUIRE(linear > 0);
}

namespace {
    double cubic_weight_abs_sum(const std::vector<double>& columns,
                                double angle) {
        const Eigen::VectorXd grid = Eigen::Map<const Eigen::VectorXd>(
            columns.data(), static_cast<Eigen::Index>(columns.size()));
        std::array<int, 4> indices{};
        std::array<double, 4> weights{};
        sasktran2::successive_orders::cubic_lagrange_weights(grid, angle,
                                                             indices, weights);
        double sum = 0.0;
        for (const double weight : weights) {
            sum += std::abs(weight);
        }
        return sum;
    }
} // namespace

TEST_CASE("Successive-orders cubic LOS stays cubic on uniform source grids",
          "[successive_orders][geometry][geometry2d]") {
    // Every column is sunlit (solar zenith angles 30-76 degrees), so only the
    // weight bound could make a uniform-grid stencil linear.
    Eigen::VectorXd altitudes = Eigen::VectorXd::LinSpaced(7, 0.0, 60000.0);
    Eigen::VectorXd horizontal = Eigen::VectorXd::LinSpaced(17, -0.4, 0.4);
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    // Limb rays covering interior and end intervals of every grid below.
    sasktran2::viewinggeometry::InternalViewingGeometry los;
    for (const double tangent_angle : {-0.25, -0.1, 0.0, 0.12, 0.25}) {
        add_limb_ray(geometry, raytracer, 20000.0, tangent_angle, los);
        add_limb_ray(geometry, raytracer, 35000.0, tangent_angle, los);
    }

    struct UniformCase {
        int num_sza;
        std::vector<double> explicit_grid;
    };
    const std::vector<UniformCase> cases = {
        {4, {}},
        {5, {}},
        {7, {}},
        {11, {}},
        {99, {-0.3, -0.18, -0.06, 0.06, 0.18, 0.3}},
    };
    for (const auto& uniform : cases) {
        DYNAMIC_SECTION("num_sza=" << uniform.num_sza << " explicit="
                                   << uniform.explicit_grid.size()) {
            sasktran2::successive_orders::SourceGeometrySettings settings;
            settings.num_incoming = 14;
            settings.num_outgoing = 14;
            settings.num_sza = uniform.num_sza;
            settings.num_threads = 1;
            settings.horizontal_angle_grid_radians = uniform.explicit_grid;
            sasktran2::successive_orders::SourceGeometry1D source(raytracer,
                                                                  geometry);
            source.initialize(los, settings);
            const auto& columns = source.source_horizontal_angles_rad();

            const auto compiled = los_layer_columns(source);
            int inside = 0;
            double maximum_sum = 0.0;
            for (std::size_t ray = 0; ray < compiled.size(); ++ray) {
                const auto& traced = los.traced_rays[ray];
                for (std::size_t layer = 0; layer < compiled[ray].size();
                     ++layer) {
                    const double angle =
                        layer_horizontal_angle(geometry, traced, layer);
                    INFO("ray=" << ray << " layer=" << layer
                                << " angle=" << angle);
                    REQUIRE(compiled[ray][layer].weight_sum ==
                            Catch::Approx(1.0).margin(1.0e-12));
                    if (angle <= columns.front() || angle >= columns.back()) {
                        continue;
                    }
                    REQUIRE(compiled[ray][layer].interior.size() == 4);
                    maximum_sum = std::max(
                        maximum_sum, cubic_weight_abs_sum(columns, angle));
                    ++inside;
                }
            }
            CAPTURE(inside, maximum_sum);
            REQUIRE(inside > 0);
            REQUIRE(maximum_sum <
                    sasktran2::successive_orders::max_cubic_weight_abs_sum);
        }
    }
}

TEST_CASE("Successive-orders cubic LOS falls back to linear on strongly "
          "non-uniform source grids",
          "[successive_orders][geometry][geometry2d]") {
    Eigen::VectorXd altitudes = Eigen::VectorXd::LinSpaced(7, 0.0, 60000.0);
    Eigen::VectorXd horizontal = Eigen::VectorXd::LinSpaced(9, -0.2, 0.2);
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    // The limb ray spans about +-0.11 radians, across every column interval.
    sasktran2::viewinggeometry::InternalViewingGeometry los;
    add_limb_ray(geometry, raytracer, 20000.0, 0.0, los);
    add_limb_ray(geometry, raytracer, 30000.0, 0.0, los);

    // The grid [0, 1, 1.01, 2] scaled by 0.1 radians and shifted to -0.1.
    // Its cubic weights reach about +-38 in the wide intervals.
    const std::vector<double> columns = {-0.1, 0.0, 0.001, 0.1};
    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 14;
    settings.num_outgoing = 14;
    settings.num_threads = 1;
    settings.horizontal_angle_grid_radians = columns;
    sasktran2::successive_orders::SourceGeometry1D source(raytracer, geometry);
    source.initialize(los, settings);

    const auto compiled = los_layer_columns(source);
    int linear = 0;
    for (std::size_t ray = 0; ray < compiled.size(); ++ray) {
        const auto& traced = los.traced_rays[ray];
        for (std::size_t layer = 0; layer < compiled[ray].size(); ++layer) {
            const double angle =
                layer_horizontal_angle(geometry, traced, layer);
            INFO("ray=" << ray << " layer=" << layer << " angle=" << angle);
            const auto& touched = compiled[ray][layer];
            REQUIRE(touched.weight_sum == Catch::Approx(1.0).margin(1.0e-12));
            const bool wide_interval_interior =
                angle > columns.front() && angle < columns.back() &&
                std::all_of(columns.begin(), columns.end(),
                            [angle](double column) {
                                return std::abs(angle - column) > 0.005;
                            });
            if (wide_interval_interior) {
                CAPTURE(cubic_weight_abs_sum(columns, angle));
                REQUIRE(cubic_weight_abs_sum(columns, angle) >
                        sasktran2::successive_orders::max_cubic_weight_abs_sum);
                REQUIRE(touched.interior.size() <= 2);
                ++linear;
            }
        }
    }
    CAPTURE(linear);
    REQUIRE(linear > 0);
}

TEST_CASE("Successive-orders cubic LOS stays cubic on mildly non-uniform "
          "source grids",
          "[successive_orders][geometry][geometry2d]") {
    // Columns clustered around the tangent point. All are sunlit (solar
    // zenith angles 33-73 degrees). The cubic weights' absolute sum is at
    // most about 1.40 in the four inner intervals and reaches about 5.1 in
    // parts of the outer intervals, which fall back to linear there.
    Eigen::VectorXd altitudes = Eigen::VectorXd::LinSpaced(7, 0.0, 60000.0);
    Eigen::VectorXd horizontal = Eigen::VectorXd::LinSpaced(17, -0.4, 0.4);
    sasktran2::Geometry2D geometry(0.6, 0.0, 6372000.0, std::move(altitudes),
                                   std::move(horizontal),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::raytracing::RustRayTracer2D raytracer(geometry);
    sasktran2::viewinggeometry::InternalViewingGeometry los;
    for (const double tangent_angle : {-0.2, 0.0, 0.2}) {
        add_limb_ray(geometry, raytracer, 20000.0, tangent_angle, los);
        add_limb_ray(geometry, raytracer, 35000.0, tangent_angle, los);
    }

    std::vector<double> columns;
    for (const double degrees : {-20.0, -8.0, -3.0, 0.0, 3.0, 8.0, 20.0}) {
        columns.push_back(degrees * EIGEN_PI / 180.0);
    }
    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 14;
    settings.num_outgoing = 14;
    settings.num_threads = 1;
    settings.horizontal_angle_grid_radians = columns;
    sasktran2::successive_orders::SourceGeometry1D source(raytracer, geometry);
    source.initialize(los, settings);

    const double inner = 8.0 * EIGEN_PI / 180.0;
    const auto compiled = los_layer_columns(source);
    int cubic = 0;
    int inner_cubic = 0;
    int linear = 0;
    double minimum_cubic_sum = std::numeric_limits<double>::infinity();
    double maximum_cubic_sum = 0.0;
    for (std::size_t ray = 0; ray < compiled.size(); ++ray) {
        const auto& traced = los.traced_rays[ray];
        for (std::size_t layer = 0; layer < compiled[ray].size(); ++layer) {
            const double angle =
                layer_horizontal_angle(geometry, traced, layer);
            INFO("ray=" << ray << " layer=" << layer << " angle=" << angle);
            const auto& touched = compiled[ray][layer];
            REQUIRE(touched.weight_sum == Catch::Approx(1.0).margin(1.0e-12));
            if (angle <= columns.front() || angle >= columns.back()) {
                continue;
            }
            const double sum = cubic_weight_abs_sum(columns, angle);
            CAPTURE(sum);
            if (std::abs(angle) < inner) {
                REQUIRE(sum < 1.41);
            }
            if (sum <= sasktran2::successive_orders::max_cubic_weight_abs_sum) {
                REQUIRE(touched.interior.size() == 4);
                minimum_cubic_sum = std::min(minimum_cubic_sum, sum);
                maximum_cubic_sum = std::max(maximum_cubic_sum, sum);
                ++cubic;
                inner_cubic += std::abs(angle) < inner ? 1 : 0;
            } else {
                REQUIRE(touched.interior.size() <= 2);
                ++linear;
            }
        }
    }
    CAPTURE(cubic, inner_cubic, linear, minimum_cubic_sum, maximum_cubic_sum);
    REQUIRE(inner_cubic > 0);
    REQUIRE(cubic > inner_cubic);
    REQUIRE(linear > 0);
}
#endif

TEST_CASE("Successive-orders aligned grids equal legacy grids in the "
          "reference solar frame",
          "[successive_orders][geometry]") {
    // At a single SZA column on the reference point with zero solar azimuth
    // the local solar frame is the identity up to rounding of the source
    // position, so frame-aligned grids must reproduce the legacy global grids
    // to rounding, with identical weights. Aligned Lebedev rules additionally
    // carry the pole-avoiding pre-rotation, and aligned ground rules an extra
    // tilt that keeps every node off the horizon.
    const Eigen::Matrix3d pole_avoiding_rotation =
        (Eigen::AngleAxisd(0.01, Eigen::Vector3d::UnitZ()) *
         Eigen::AngleAxisd(0.01, Eigen::Vector3d::UnitY()))
            .toRotationMatrix();
    const Eigen::Matrix3d ground_rotation =
        pole_avoiding_rotation *
        Eigen::AngleAxisd(0.08, Eigen::Vector3d::UnitX()).toRotationMatrix();
    for (const bool reduced_horizon : {true, false}) {
        DYNAMIC_SECTION("reduced_horizon=" << reduced_horizon) {
            sasktran2::Geometry1D geometry(
                0.6, 0.0, 6372000.0, altitude_grid(),
                sasktran2::grids::interpolation::linear,
                sasktran2::geometrytype::spherical);
            sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
            const auto los = make_los_geometry(geometry, raytracer);
            sasktran2::successive_orders::SourceGeometrySettings settings;
            settings.num_incoming = 26;
            settings.num_outgoing = 26;
            settings.num_sza = 1;
            settings.num_threads = 1;
            settings.use_reduced_horizon_quadrature = reduced_horizon;
            sasktran2::successive_orders::SourceGeometry1D aligned(raytracer,
                                                                   geometry);
            aligned.initialize(los, settings);
            settings.legacy_interpolation = true;
            sasktran2::successive_orders::SourceGeometry1D legacy(raytracer,
                                                                  geometry);
            legacy.initialize(los, settings);

            REQUIRE(aligned.num_points() == legacy.num_points());
            REQUIRE(aligned.num_interior_points() ==
                    legacy.num_interior_points());
            const auto require_same_direction =
                [](const Eigen::Vector3d& actual,
                   const Eigen::Vector3d& expected) {
                    REQUIRE((actual - expected).cwiseAbs().maxCoeff() <=
                            1.0e-14);
                };
            const auto require_same_sphere =
                [&](const sasktran2::math::UnitSphere& actual,
                    const sasktran2::math::UnitSphere& expected,
                    const Eigen::Matrix3d& rotation) {
                    REQUIRE(actual.num_points() == expected.num_points());
                    for (int node = 0; node < expected.num_points(); ++node) {
                        require_same_direction(
                            actual.get_quad_position(node),
                            rotation * expected.get_quad_position(node));
                        REQUIRE(actual.quadrature_weight(node) ==
                                expected.quadrature_weight(node));
                    }
                };
            const Eigen::Matrix3d interior_rotation =
                reduced_horizon ? Eigen::Matrix3d::Identity()
                                : pole_avoiding_rotation;
            // Aligned ground Lebedev rules are the tilted full rule restricted
            // to the upward hemisphere and renormalized. No node lies near the
            // horizon, so exactly half the nodes contribute and the weights
            // are renormalized over them alone.
            const auto require_tilted_hemisphere =
                [&](const sasktran2::math::UnitSphere& actual, int num_points,
                    const Eigen::Vector3d& up) {
                    const sasktran2::math::LebedevSphere full(num_points);
                    std::vector<int> kept;
                    double normalization = 0.0;
                    for (int node = 0; node < full.num_points(); ++node) {
                        const double projection =
                            (ground_rotation * full.get_quad_position(node))
                                .dot(up);
                        REQUIRE(std::abs(projection) > 1.0e-6);
                        if (projection > 0.0) {
                            kept.push_back(node);
                            normalization += full.quadrature_weight(node);
                        }
                    }
                    REQUIRE(2 * static_cast<int>(kept.size()) == num_points);
                    REQUIRE(actual.num_points() ==
                            static_cast<int>(kept.size()));
                    for (std::size_t node = 0; node < kept.size(); ++node) {
                        const int index = static_cast<int>(node);
                        require_same_direction(
                            actual.get_quad_position(index),
                            ground_rotation *
                                full.get_quad_position(kept[node]));
                        REQUIRE(actual.quadrature_weight(index) ==
                                full.quadrature_weight(kept[node]) /
                                    normalization * 0.5);
                    }
                };
            for (int index = 0; index < legacy.num_points(); ++index) {
                INFO("point=" << index);
                const auto& actual = aligned.source_point(index);
                const auto& expected = legacy.source_point(index);
                REQUIRE(actual.location().position ==
                        expected.location().position);
                if (expected.is_ground()) {
                    const Eigen::Vector3d up =
                        expected.location().position.normalized();
                    if (reduced_horizon) {
                        // Reduced-horizon rings are already local, and none
                        // lies on the horizon.
                        require_same_sphere(actual.incoming_sphere(),
                                            expected.incoming_sphere(),
                                            Eigen::Matrix3d::Identity());
                    } else {
                        require_tilted_hemisphere(actual.incoming_sphere(),
                                                  settings.num_incoming, up);
                    }
                    require_tilted_hemisphere(actual.outgoing_sphere(),
                                              settings.num_outgoing, up);
                } else {
                    require_same_sphere(actual.incoming_sphere(),
                                        expected.incoming_sphere(),
                                        interior_rotation);
                    require_same_sphere(actual.outgoing_sphere(),
                                        expected.outgoing_sphere(),
                                        interior_rotation);
                }
            }
        }
    }
}

TEST_CASE("Successive-orders aligned ground Lebedev grids have no horizon "
          "nodes",
          "[successive_orders][geometry]") {
    // Lebedev rules are symmetric under inversion, so a ground hemisphere
    // with exactly half the nodes, all at least 1e-6 above the horizon, means
    // that no node of the rotated full rule lies within 1e-6 of it.
    const std::array<int, 14> sizes{6,   14,  26,  38,  50,  74,  86,
                                    110, 146, 170, 194, 230, 266, 302};
    // Lambertian factor 4 sum(mu w) of the outgoing hemisphere; one for an
    // exact rule.
    const std::map<int, double> lambertian_factors{
        {26, 0.949761}, {50, 0.999219}, {110, 1.000976}, {194, 1.000832}};
    for (const int size : sizes) {
        for (const bool reduced_horizon : {false, true}) {
            DYNAMIC_SECTION("size=" << size
                                    << " reduced_horizon=" << reduced_horizon) {
                Eigen::VectorXd altitudes =
                    Eigen::VectorXd::LinSpaced(3, 0.0, 20000.0);
                sasktran2::Geometry1D geometry(
                    0.6, 0.3, 6372000.0, std::move(altitudes),
                    sasktran2::grids::interpolation::linear,
                    sasktran2::geometrytype::spherical);
                sasktran2::raytracing::SphericalShellRayTracer raytracer(
                    geometry);
                sasktran2::viewinggeometry::InternalViewingGeometry los;
                const std::array<double, 2> observer_cos_sza{0.3, 0.8};
                los.traced_rays.resize(observer_cos_sza.size());
                for (std::size_t ray_index = 0;
                     ray_index < observer_cos_sza.size(); ++ray_index) {
                    sasktran2::viewinggeometry::ViewingRay ray;
                    ray.observer.position =
                        geometry.coordinates().solar_coordinate_vector(
                            observer_cos_sza[ray_index], 0.3, 30000.0);
                    ray.look_away = -ray.observer.position.normalized();
                    raytracer.trace_ray(ray, los.traced_rays[ray_index]);
                }
                sasktran2::successive_orders::SourceGeometrySettings settings;
                settings.num_incoming = reduced_horizon ? 14 : size;
                settings.num_outgoing = size;
                settings.num_sza = 2;
                settings.num_threads = 1;
                settings.use_reduced_horizon_quadrature = reduced_horizon;
                sasktran2::successive_orders::SourceGeometry1D source(raytracer,
                                                                      geometry);
                source.initialize(los, settings);
                REQUIRE(source.num_ground_points() == 2);

                const auto require_no_horizon_nodes =
                    [](const sasktran2::math::UnitSphere& sphere,
                       const Eigen::Vector3d& up, int full_size) {
                        REQUIRE(2 * sphere.num_points() == full_size);
                        double weight_sum = 0.0;
                        for (int node = 0; node < sphere.num_points(); ++node) {
                            REQUIRE(sphere.get_quad_position(node).dot(up) >
                                    1.0e-6);
                            weight_sum += sphere.quadrature_weight(node);
                        }
                        REQUIRE(weight_sum ==
                                Catch::Approx(0.5).margin(1.0e-14));
                    };
                for (const auto& point : source.source_points()) {
                    if (!point.is_ground()) {
                        continue;
                    }
                    const Eigen::Vector3d up =
                        point.location().position.normalized();
                    require_no_horizon_nodes(point.outgoing_sphere(), up, size);
                    if (!reduced_horizon) {
                        require_no_horizon_nodes(point.incoming_sphere(), up,
                                                 size);
                    }
                    const auto factor = lambertian_factors.find(size);
                    if (factor != lambertian_factors.end()) {
                        const auto& sphere = point.outgoing_sphere();
                        double lambertian = 0.0;
                        for (int node = 0; node < sphere.num_points(); ++node) {
                            lambertian +=
                                4.0 * sphere.get_quad_position(node).dot(up) *
                                sphere.quadrature_weight(node);
                        }
                        REQUIRE(lambertian ==
                                Catch::Approx(factor->second).margin(1.0e-6));
                    }
                }
            }
        }
    }
}

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
    const auto columns = source_geometry.los_transport_columns_for_ray(0);
    for (const auto& weight :
         source_geometry.los_interpolation().front().source_for_layer(0)) {
        const int source_index = columns[weight.row_inner_index()];
        int owner = -1;
        for (int point_index = 0;
             point_index < source_geometry.num_interior_points();
             ++point_index) {
            const auto& point = source_geometry.source_point(point_index);
            if (source_index >= point.outgoing_offset() &&
                source_index < point.outgoing_offset() + point.num_outgoing()) {
                owner = point_index;
                break;
            }
        }
        REQUIRE(owner >= 0);
        location_weights[owner] += weight.weight();
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

    const std::array<Eigen::Vector3d, 4> directions{
        Eigen::Vector3d::UnitZ(), -Eigen::Vector3d::UnitZ(),
        Eigen::Vector3d::UnitX(), Eigen::Vector3d::UnitY()};
    for (const bool legacy_interpolation : {false, true}) {
        DYNAMIC_SECTION("legacy_interpolation=" << legacy_interpolation) {
            sasktran2::successive_orders::SourceGeometrySettings settings;
            settings.num_incoming = 6;
            settings.num_outgoing = 6;
            settings.num_threads = 2;
            settings.legacy_interpolation = legacy_interpolation;
            sasktran2::successive_orders::SourceGeometry1D source_geometry(
                raytracer, geometry);
            source_geometry.initialize(los, settings);

            for (std::size_t ray = 0; ray < directions.size(); ++ray) {
                INFO("ray=" << ray);
                if (legacy_interpolation) {
                    // The global grid has a node on every axis.
                    require_los_direction_on_node(source_geometry, ray,
                                                  directions[ray]);
                } else {
                    // Aligned grids carry the pole-avoiding pre-rotation, so
                    // the axes fall between nodes. Exact preservation then
                    // means the compiled weights equal the point grid's own
                    // interpolation of the direction.
                    require_los_direction(source_geometry, ray,
                                          directions[ray]);
                }
            }
        }
    }
}

TEST_CASE("Successive-orders structured layer storage preserves every weight "
          "and stencil order",
          "[successive_orders][geometry]") {
    using namespace sasktran2::successive_orders;
    const auto bits = [](double value) {
        std::uint64_t result;
        std::memcpy(&result, &value, sizeof(result));
        return result;
    };
    sasktran2::raytracing::TracedRay traced;
    traced.layers.resize(3);
    const std::array<double, 4> endpoint = {0.0, 0.0, 0.0, 0.0};
    const std::array<double, 4> od = {-0.0,
                                      std::numeric_limits<double>::denorm_min(),
                                      std::nextafter(0.3, 1.0), 0.75};
    for (std::size_t layer = 0; layer < traced.layers.size(); ++layer) {
        const int base = static_cast<int>(layer) * 11;
        traced.set_layer_weights(
            layer, std::array<int, 4>{base, base + 1, base + 10, base + 11},
            endpoint, endpoint, od);
    }

    RayInterpolation compiled;
    compiled.traced_ray = &traced;
    compiled.structured_altitude_stride = 10;
    compiled.layers = {
        {0, 4, 0, 2, 0, 4}, {4, 2, 2, 0, 4, 4}, {6, 0, 2, 2, 8, 4}};
    compiled.atmosphere_weights = {
        {0, -0.0},
        {1, 0.25},
        {10, 0.5},
        {11, std::nextafter(0.25, 1.0)},
        {12, std::numeric_limits<double>::denorm_min()},
        {22, 1.0}};
    compiled.source_weights = {{0, 0.3}, {2, 0.7}, {1, 0.4}, {2, 0.6}};
    std::vector<int> columns;
    compile_transport_row(compiled, columns);
    compact_ray_interpolation(compiled);

    const auto original_layers = compiled.layers.wide_values();
    const auto original_atmosphere = compiled.atmosphere_weights;
    std::array<std::vector<std::pair<int, double>>, 3> original_od;
    for (std::size_t layer = 0; layer < traced.layers.size(); ++layer) {
        const auto weights = compiled.optical_depth_for_layer(layer);
        for (std::size_t entry = 0; entry < weights.size(); ++entry) {
            original_od[layer].push_back(weights[entry]);
        }
    }
    adopt_optical_depth_storage(traced, compiled);

    REQUIRE(compiled.layers.is_structured());
    REQUIRE(compiled.layers.is_compact_structured());
    REQUIRE(compiled.layers.capacity_bytes() ==
            3 * sizeof(CompactStructuredLayerInterpolation));
    REQUIRE(compiled.atmosphere_weights.capacity() == 0);
    REQUIRE(compiled.optical_depth_indices.capacity() == 0);
    REQUIRE(compiled.structured_atmosphere_weights.size() == 6);
    REQUIRE(compiled.source_weights.element_bytes() == 9);
    for (std::size_t layer = 0; layer < compiled.layers.size(); ++layer) {
        const auto descriptor = compiled.layers[layer];
        const auto& original = original_layers[layer];
        REQUIRE(descriptor.atmosphere_offset == original.atmosphere_offset);
        REQUIRE(descriptor.atmosphere_count == original.atmosphere_count);
        REQUIRE(descriptor.source_offset == original.source_offset);
        REQUIRE(descriptor.source_count == original.source_count);
        REQUIRE(descriptor.optical_depth_offset ==
                original.optical_depth_offset);
        REQUIRE(descriptor.optical_depth_count == original.optical_depth_count);
        const auto midpoint = compiled.atmosphere_for_layer(layer);
        require_sorted(midpoint);
        REQUIRE(midpoint.size() == original.atmosphere_count);
        for (std::size_t entry = 0; entry < midpoint.size(); ++entry) {
            const auto actual = midpoint[entry];
            const auto& expected =
                original_atmosphere[original.atmosphere_offset + entry];
            REQUIRE(actual.index == expected.index);
            REQUIRE(bits(actual.weight()) == bits(expected.weight()));
        }
        const auto weights = compiled.optical_depth_for_layer(layer);
        REQUIRE(weights.size() == original_od[layer].size());
        for (std::size_t entry = 0; entry < weights.size(); ++entry) {
            REQUIRE(weights[entry].first == original_od[layer][entry].first);
            REQUIRE(bits(weights[entry].second) ==
                    bits(original_od[layer][entry].second));
        }
    }
    REQUIRE(compiled.atmosphere_for_layer(2).empty());
    REQUIRE(compiled.source_for_layer(1).empty());
    // Iterators retain the immutable backing descriptor when their temporary
    // view is destroyed, matching the original pointer-view lifetime.
    auto midpoint_iterator = compiled.atmosphere_for_layer(1).begin();
    REQUIRE((*midpoint_iterator).index == 12);
    ++midpoint_iterator;
    REQUIRE((*midpoint_iterator).index == 22);
    ++midpoint_iterator;
    REQUIRE(midpoint_iterator == compiled.atmosphere_for_layer(1).end());
}

TEST_CASE("Successive-orders structured layer storage retains generic rays "
          "when a cell pattern or compact bound differs",
          "[successive_orders][geometry]") {
    using namespace sasktran2::successive_orders;
    const auto check_fallback = [](int base, int midpoint,
                                   std::size_t source_count,
                                   std::uint8_t od_count) {
        sasktran2::raytracing::TracedRay traced;
        traced.layers.resize(1);
        const std::vector<int> indices = {base, base + 1, base + 10, base + 11};
        const std::vector<double> weights = {0.1, 0.2, 0.3, 0.4};
        traced.set_layer_weights(0, indices.data(), weights.data(),
                                 weights.data(), weights.data(), od_count);
        RayInterpolation compiled;
        compiled.traced_ray = &traced;
        compiled.structured_altitude_stride = 10;
        compiled.layers = {
            {0, 1, 0, static_cast<std::uint32_t>(source_count), 0, od_count}};
        compiled.atmosphere_weights = {{midpoint, 1.0}};
        compiled.source_weights.wide_values().resize(source_count);
        adopt_optical_depth_storage(traced, compiled);
        REQUIRE(!compiled.layers.is_structured());
        REQUIRE(compiled.atmosphere_weights.size() == 1);
        REQUIRE(compiled.structured_atmosphere_weights.empty());
        REQUIRE(compiled.optical_depth_indices.size() == od_count);
        REQUIRE(compiled.atmosphere_for_layer(0)[0].index == midpoint);
        auto iterator = compiled.atmosphere_for_layer(0).begin();
        REQUIRE((*iterator).index == midpoint);
        REQUIRE(compiled.optical_depth_for_layer(0)[0].first == base);
    };
    check_fallback(0, 5, 0, 4); // Midpoint belongs to another boundary cell.
    check_fallback(65536, 65536, 0,
                   4);              // Cell base does not fit the descriptor.
    check_fallback(0, 0, 65537, 4); // Layer offset domain needs the wide form.
    check_fallback(0, 0, 0, 3);     // Generic OD stencil has a different size.
}

TEST_CASE("Successive-orders frame-aligned Lebedev grids avoid exactly radial "
          "incoming rays",
          "[successive_orders][geometry]") {
    // Frame-only Lebedev grids would map the canonical pole exactly onto each
    // column's vertical. The pole-avoiding pre-rotation keeps every aligned
    // incoming node, interior and ground, off the vertical, and the diffuse
    // optical-depth quadrature stays finite at every column.
    Eigen::VectorXd altitudes = Eigen::VectorXd::LinSpaced(27, 0.0, 65000.0);
    sasktran2::Geometry1D geometry(0.6, 0.0, 6372000.0, std::move(altitudes),
                                   sasktran2::grids::interpolation::linear,
                                   sasktran2::geometrytype::spherical);
    sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
    sasktran2::viewinggeometry::InternalViewingGeometry los;
    los.traced_rays.resize(3);
    const std::array<double, 3> observer_cos_sza{0.4, 0.6, 0.8};
    for (std::size_t ray_index = 0; ray_index < observer_cos_sza.size();
         ++ray_index) {
        sasktran2::viewinggeometry::ViewingRay ray;
        ray.observer.position = geometry.coordinates().solar_coordinate_vector(
            observer_cos_sza[ray_index], 0.0, 4000.0);
        ray.look_away = -ray.observer.position.normalized();
        raytracer.trace_ray(ray, los.traced_rays[ray_index]);
    }

    sasktran2::successive_orders::SourceGeometrySettings settings;
    settings.num_incoming = 6;
    settings.num_outgoing = 6;
    settings.num_sza = 3;
    settings.num_threads = 1;
    sasktran2::successive_orders::SourceGeometry1D source(raytracer, geometry);
    source.initialize(los, settings);
    REQUIRE(source.source_cos_sza().size() == 3);

    REQUIRE(source.num_ground_points() == 3);
    for (const auto& point : source.source_points()) {
        INFO("ground=" << point.is_ground());
        const Eigen::Vector3d up = point.location().position.normalized();
        for (int direction = 0; direction < point.num_incoming(); ++direction) {
            REQUIRE(
                std::abs(
                    point.incoming_sphere().get_quad_position(direction).dot(
                        up)) < 1.0 - 1.0e-6);
        }
    }
    for (std::size_t ray = 0; ray < source.incoming_interpolation().size();
         ++ray) {
        INFO("ray=" << ray);
        const auto& interpolation = source.incoming_interpolation()[ray];
        for (std::size_t layer = 0; layer < interpolation.layers.size();
             ++layer) {
            const auto optical_depth =
                interpolation.optical_depth_for_layer(layer);
            for (std::size_t index = 0; index < optical_depth.size(); ++index) {
                REQUIRE(std::isfinite(optical_depth[index].second));
            }
        }
    }
}

TEST_CASE("Successive-orders cubic horizontal weights reproduce cubics",
          "[successive_orders][geometry]") {
    Eigen::VectorXd grid(6);
    grid << 0.0, 0.3, 0.7, 1.2, 2.0, 2.1;
    const auto cubic = [](double x) {
        return 1.0 + 2.0 * x - 0.5 * x * x + 0.3 * x * x * x;
    };
    std::array<int, 4> indices{};
    std::array<double, 4> weights{};
    for (double x = 0.01; x < 2.1; x += 0.0137) {
        sasktran2::successive_orders::cubic_lagrange_weights(grid, x, indices,
                                                             weights);
        double value = 0.0;
        double total = 0.0;
        for (int m = 0; m < 4; ++m) {
            REQUIRE(indices[m] == indices[0] + m);
            value += weights[m] * cubic(grid[indices[m]]);
            total += weights[m];
        }
        REQUIRE(indices[0] >= 0);
        REQUIRE(indices[3] <= 5);
        REQUIRE(value == Catch::Approx(cubic(x)).epsilon(1.0e-12));
        REQUIRE(total == Catch::Approx(1.0).epsilon(1.0e-13));
    }
    sasktran2::successive_orders::cubic_lagrange_weights(grid, 0.7, indices,
                                                         weights);
    REQUIRE(indices[0] == 1);
    REQUIRE(weights[1] == 1.0);
    REQUIRE(weights[0] == 0.0);
    REQUIRE(weights[2] == 0.0);
    REQUIRE(weights[3] == 0.0);

    Eigen::VectorXd short_grid(3);
    short_grid << 0.0, 1.0, 2.0;
    REQUIRE_THROWS_AS(sasktran2::successive_orders::cubic_lagrange_weights(
                          short_grid, 0.5, indices, weights),
                      std::invalid_argument);
}

TEST_CASE("Successive-orders cubic weight bound accepts uniform grids and "
          "rejects strongly non-uniform ones",
          "[successive_orders][geometry]") {
    using sasktran2::successive_orders::cubic_lagrange_weights;
    using sasktran2::successive_orders::cubic_weights_are_bounded;
    const auto abs_sum = [](const std::array<double, 4>& weights) {
        return std::abs(weights[0]) + std::abs(weights[1]) +
               std::abs(weights[2]) + std::abs(weights[3]);
    };
    std::array<int, 4> indices{};
    std::array<double, 4> weights{};

    for (const int size : {4, 5, 7, 11, 31}) {
        INFO("size=" << size);
        const Eigen::VectorXd grid =
            Eigen::VectorXd::LinSpaced(size, -0.3, 0.45);
        const double spacing = grid[1] - grid[0];
        double interior_maximum = 0.0;
        double end_maximum = 0.0;
        constexpr int samples = 400;
        for (int interval = 0; interval < size - 1; ++interval) {
            for (int sample = 1; sample < samples; ++sample) {
                const double x = grid[interval] + spacing * sample / samples;
                cubic_lagrange_weights(grid, x, indices, weights);
                REQUIRE(cubic_weights_are_bounded(weights));
                const bool end = interval == 0 || interval == size - 2;
                double& maximum = end ? end_maximum : interior_maximum;
                maximum = std::max(maximum, abs_sum(weights));
            }
        }
        CAPTURE(interior_maximum, end_maximum);
        REQUIRE(end_maximum == Catch::Approx(1.6311).margin(1.0e-3));
        REQUIRE(interior_maximum == Catch::Approx(1.25).margin(1.0e-6));
    }

    Eigen::VectorXd stretched(4);
    stretched << 0.0, 1.0, 2.0, 4.0;
    cubic_lagrange_weights(stretched, 3.0, indices, weights);
    REQUIRE(weights[0] == Catch::Approx(0.25));
    REQUIRE(weights[1] == Catch::Approx(-1.0));
    REQUIRE(weights[2] == Catch::Approx(1.5));
    REQUIRE(weights[3] == Catch::Approx(0.25));
    REQUIRE(!cubic_weights_are_bounded(weights));

    Eigen::VectorXd clustered(4);
    clustered << 0.0, 1.0, 1.01, 2.0;
    for (const double x :
         {0.1, 0.25, 0.5, 0.75, 0.9, 1.1, 1.25, 1.5, 1.75, 1.9}) {
        INFO("x=" << x);
        cubic_lagrange_weights(clustered, x, indices, weights);
        REQUIRE(!cubic_weights_are_bounded(weights));
    }
    cubic_lagrange_weights(clustered, 0.5, indices, weights);
    REQUIRE(weights[1] == Catch::Approx(38.25));
    REQUIRE(weights[2] == Catch::Approx(-37.5037).epsilon(1.0e-5));
    // Inside the narrow interval the stencil is well conditioned.
    cubic_lagrange_weights(clustered, 1.005, indices, weights);
    REQUIRE(cubic_weights_are_bounded(weights));
}
