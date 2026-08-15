#include "geometry.h"

#include <sasktran2/grids.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>

#ifdef SKTRAN_OPENMP_SUPPORT
#include <omp.h>
#endif

namespace sasktran2::successive_orders {
    namespace {
        std::vector<InterpolationWeight>
        sorted_weights(const std::vector<std::pair<int, double>>& weights) {
            std::vector<InterpolationWeight> result;
            result.reserve(weights.size());
            for (const auto& [index, weight] : weights) {
                if (index < 0 || !std::isfinite(weight)) {
                    throw std::runtime_error(
                        "Invalid successive-orders point interpolation");
                }
                if (weight != 0.0) {
                    result.push_back({index, weight});
                }
            }
            std::stable_sort(result.begin(), result.end(),
                             [](const auto& left, const auto& right) {
                                 return left.index < right.index;
                             });
            std::size_t write = 0;
            for (const auto& value : result) {
                if (write != 0 && result[write - 1].index == value.index) {
                    result[write - 1].weight += value.weight;
                } else {
                    result[write++] = value;
                }
            }
            result.resize(write);
            result.erase(std::remove_if(result.begin(), result.end(),
                                        [](const auto& value) {
                                            return value.weight == 0.0;
                                        }),
                         result.end());
            return result;
        }

        int checked_add(int left, int right, const char* description) {
            if (right < 0 || left > std::numeric_limits<int>::max() - right) {
                throw std::length_error(std::string(description) +
                                        " exceeds the supported index range");
            }
            return left + right;
        }
    } // namespace

    void SourceGeometrySettings::validate() const {
        if (num_incoming <= 0 || num_outgoing <= 0) {
            throw std::invalid_argument(
                "Successive-orders angular grids must be non-empty");
        }
        if (num_sza <= 0) {
            throw std::invalid_argument(
                "Successive-orders SZA grid must be non-empty");
        }
        if (num_threads <= 0) {
            throw std::invalid_argument(
                "Successive-orders geometry requires at least one thread");
        }
    }

    SourceGeometry1D::SourceGeometry1D(
        const sasktran2::raytracing::RayTracerBase& raytracer,
        const sasktran2::Geometry1D& geometry)
        : m_raytracer(raytracer), m_geometry(geometry) {}

    SourceGeometry1D::~SourceGeometry1D() = default;

    sasktran2::grids::AltitudeGrid SourceGeometry1D::make_altitude_grid() {
        const auto& atmosphere_altitudes = m_geometry.altitude_grid().grid();
        Eigen::VectorXd altitudes;
        // Atmosphere layer midpoints inherit arbitrary spacing from the
        // atmosphere grid. Let Grid retain the constant fast path only when
        // those midpoints are actually uniform.
        const auto spacing = sasktran2::grids::gridspacing::automatic;

        if (m_settings.altitude_grid_m.empty()) {
            altitudes = 0.5 * (atmosphere_altitudes(Eigen::seq(
                                   0, Eigen::placeholders::last - 1)) +
                               atmosphere_altitudes(
                                   Eigen::seq(1, Eigen::placeholders::last)));
        } else {
            altitudes.resize(
                static_cast<Eigen::Index>(m_settings.altitude_grid_m.size()));
            const double minimum_altitude = atmosphere_altitudes[0];
            const double maximum_altitude =
                atmosphere_altitudes[atmosphere_altitudes.size() - 1];
            for (std::size_t index = 0;
                 index < m_settings.altitude_grid_m.size(); ++index) {
                const double altitude = m_settings.altitude_grid_m[index];
                if (!std::isfinite(altitude) ||
                    (index != 0 &&
                     altitude <= m_settings.altitude_grid_m[index - 1])) {
                    throw std::invalid_argument(
                        "Successive-orders source altitudes must be finite "
                        "and strictly increasing");
                }
                if (altitude <= minimum_altitude ||
                    altitude >= maximum_altitude) {
                    throw std::invalid_argument(
                        "Successive-orders source altitudes must lie strictly "
                        "inside the atmosphere altitude range");
                }
                altitudes[static_cast<Eigen::Index>(index)] = altitude;
            }
        }

        m_source_altitudes_m.assign(altitudes.data(),
                                    altitudes.data() + altitudes.size());
        return sasktran2::grids::AltitudeGrid(
            std::move(altitudes), spacing,
            sasktran2::grids::outofbounds::extend,
            sasktran2::grids::interpolation::linear);
    }

    sasktran2::grids::Grid SourceGeometry1D::make_cos_sza_grid(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        Eigen::VectorXd cos_sza;
        const auto geometry_type = m_geometry.coordinates().geometry_type();
        const auto bounds = sasktran2::raytracing::min_max_cos_sza_of_all_rays(
            internal_viewing.traced_rays);
        // Plane-parallel and pseudospherical Geometry1D source locations both
        // occupy one invariant reference column. Multiple SZA entries would
        // therefore duplicate the same physical source points, so retain one
        // point for those geometries even when num_sza is greater than one.
        const bool usable_range =
            geometry_type == sasktran2::geometrytype::spherical &&
            m_settings.num_sza > 1 && std::isfinite(bounds.first) &&
            std::isfinite(bounds.second) && bounds.first <= bounds.second &&
            bounds.second - bounds.first > 1.0e-14;

        if (usable_range) {
            cos_sza.setLinSpaced(m_settings.num_sza, bounds.first,
                                 bounds.second);
        } else {
            cos_sza.resize(1);
            cos_sza[0] = m_geometry.coordinates().cos_sza_at_reference();
        }

        m_source_cos_sza.assign(cos_sza.data(),
                                cos_sza.data() + cos_sza.size());
        return sasktran2::grids::Grid(std::move(cos_sza),
                                      sasktran2::grids::gridspacing::constant,
                                      sasktran2::grids::outofbounds::extend,
                                      sasktran2::grids::interpolation::linear);
    }

    void SourceGeometry1D::construct_source_points() {
        m_num_interior_points = m_location_interpolator->num_interior_points();
        m_num_ground_points = m_location_interpolator->num_ground_points();
        const int total_points =
            checked_add(m_num_interior_points, m_num_ground_points,
                        "Successive-orders source point count");

        m_angular_grids.clear();
        m_angular_grids.reserve(static_cast<std::size_t>(m_num_ground_points) +
                                1);
        auto volume_grid = std::make_unique<AngularGridPair>();
        volume_grid->incoming =
            std::make_unique<sasktran2::math::LebedevSphere>(
                m_settings.num_incoming);
        volume_grid->outgoing =
            std::make_unique<sasktran2::math::LebedevSphere>(
                m_settings.num_outgoing);
        m_angular_grids.push_back(std::move(volume_grid));

        for (int ground_index = 0; ground_index < m_num_ground_points;
             ++ground_index) {
            const Eigen::Vector3d location =
                m_location_interpolator->ground_location(
                    m_geometry.coordinates(), ground_index);
            auto ground_grid = std::make_unique<AngularGridPair>();
            ground_grid->incoming =
                std::make_unique<sasktran2::math::UnitSphereGround>(
                    std::make_unique<sasktran2::math::LebedevSphere>(
                        m_settings.num_incoming),
                    location);
            ground_grid->outgoing =
                std::make_unique<sasktran2::math::UnitSphereGround>(
                    std::make_unique<sasktran2::math::LebedevSphere>(
                        m_settings.num_outgoing),
                    location);
            m_angular_grids.push_back(std::move(ground_grid));
        }

        m_source_points.clear();
        m_source_points.resize(total_points);
        for (int point_index = 0; point_index < m_num_interior_points;
             ++point_index) {
            auto& point = m_source_points[point_index];
            point.m_location.position = m_location_interpolator->grid_location(
                m_geometry.coordinates(), point_index);
            point.m_incoming_sphere = m_angular_grids[0]->incoming.get();
            point.m_outgoing_sphere = m_angular_grids[0]->outgoing.get();
            point.m_is_ground = false;
        }
        for (int ground_index = 0; ground_index < m_num_ground_points;
             ++ground_index) {
            auto& point = m_source_points[m_num_interior_points + ground_index];
            point.m_location.position =
                m_location_interpolator->ground_location(
                    m_geometry.coordinates(), ground_index);
            // Avoid tracing a nominal ground point to the wrong side of the
            // lower boundary because of Cartesian roundoff.
            point.m_location.position +=
                0.01 * point.m_location.position.normalized();
            point.m_incoming_sphere =
                m_angular_grids[ground_index + 1]->incoming.get();
            point.m_outgoing_sphere =
                m_angular_grids[ground_index + 1]->outgoing.get();
            point.m_is_ground = true;
        }

        m_incoming_point_offsets.assign(
            static_cast<std::size_t>(total_points) + 1, 0);
        m_outgoing_point_offsets.assign(
            static_cast<std::size_t>(total_points) + 1, 0);
        std::vector<std::pair<int, double>> atmosphere_weights;
        for (int point_index = 0; point_index < total_points; ++point_index) {
            auto& point = m_source_points[point_index];
            point.m_incoming_offset = m_incoming_point_offsets[point_index];
            point.m_outgoing_offset = m_outgoing_point_offsets[point_index];
            m_incoming_point_offsets[point_index + 1] =
                checked_add(point.m_incoming_offset, point.num_incoming(),
                            "Successive-orders incoming direction count");
            m_outgoing_point_offsets[point_index + 1] =
                checked_add(point.m_outgoing_offset, point.num_outgoing(),
                            "Successive-orders outgoing direction count");

            m_geometry.assign_interpolation_weights(point.location(),
                                                    atmosphere_weights);
            point.m_atmosphere_weights = sorted_weights(atmosphere_weights);
        }
    }

    void SourceGeometry1D::trace_and_compile_incoming() {
        const int num_rays = total_num_incoming();
        m_incoming_viewing.traced_rays.clear();
        m_incoming_viewing.flux_observers.clear();
        m_incoming_viewing.traced_rays.resize(num_rays);
        m_incoming_interpolation.clear();
        m_incoming_interpolation.resize(num_rays);

        std::vector<sasktran2::viewinggeometry::ViewingRay> viewing_rays(
            m_settings.num_threads);
        std::vector<InterpolationScratch> interpolation_scratch(
            m_settings.num_threads);
        std::vector<std::exception_ptr> thread_exceptions(
            m_settings.num_threads);
        std::atomic<bool> failed{false};

#pragma omp parallel for num_threads(m_settings.num_threads) schedule(dynamic)
        for (int point_index = 0;
             point_index < static_cast<int>(m_source_points.size());
             ++point_index) {
#ifdef SKTRAN_OPENMP_SUPPORT
            const int thread_index = omp_get_thread_num();
#else
            const int thread_index = 0;
#endif
            if (failed.load(std::memory_order_relaxed)) {
                continue;
            }
            try {
                const auto& point = m_source_points[point_index];
                auto& viewing_ray = viewing_rays[thread_index];
                viewing_ray.observer = point.location();
                for (int direction_index = 0;
                     direction_index < point.num_incoming();
                     ++direction_index) {
                    viewing_ray.look_away =
                        point.incoming_sphere().get_quad_position(
                            direction_index);
                    const int ray_index =
                        point.incoming_offset() + direction_index;
                    auto& traced_ray =
                        m_incoming_viewing.traced_rays[ray_index];
                    m_raytracer.trace_ray(viewing_ray, traced_ray,
                                          m_settings.include_refraction);
                    compile_ray_interpolation(
                        traced_ray, m_geometry, *m_location_interpolator,
                        m_source_points, m_incoming_interpolation[ray_index],
                        interpolation_scratch[thread_index]);
                }
            } catch (...) {
                thread_exceptions[thread_index] = std::current_exception();
                failed.store(true, std::memory_order_relaxed);
            }
        }

        for (const auto& exception : thread_exceptions) {
            if (exception) {
                std::rethrow_exception(exception);
            }
        }
    }

    void SourceGeometry1D::compile_los_interpolation(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        const int num_rays =
            static_cast<int>(internal_viewing.traced_rays.size());
        m_los_interpolation.clear();
        m_los_interpolation.resize(num_rays);
        std::vector<InterpolationScratch> interpolation_scratch(
            m_settings.num_threads);
        std::vector<std::exception_ptr> thread_exceptions(
            m_settings.num_threads);
        std::atomic<bool> failed{false};

#pragma omp parallel for num_threads(m_settings.num_threads) schedule(dynamic)
        for (int ray_index = 0; ray_index < num_rays; ++ray_index) {
#ifdef SKTRAN_OPENMP_SUPPORT
            const int thread_index = omp_get_thread_num();
#else
            const int thread_index = 0;
#endif
            if (failed.load(std::memory_order_relaxed)) {
                continue;
            }
            try {
                compile_ray_interpolation(
                    internal_viewing.traced_rays[ray_index], m_geometry,
                    *m_location_interpolator, m_source_points,
                    m_los_interpolation[ray_index],
                    interpolation_scratch[thread_index]);
            } catch (...) {
                thread_exceptions[thread_index] = std::current_exception();
                failed.store(true, std::memory_order_relaxed);
            }
        }

        for (const auto& exception : thread_exceptions) {
            if (exception) {
                std::rethrow_exception(exception);
            }
        }
    }

    void SourceGeometry1D::compile_transport_topology(
        std::vector<RayInterpolation>& interpolation,
        std::vector<int>& row_offsets, std::vector<int>& column_indices) {
        row_offsets.assign(interpolation.size() + 1, 0);
        std::vector<int> columns;
        std::size_t num_columns = 0;
        for (std::size_t ray_index = 0; ray_index < interpolation.size();
             ++ray_index) {
            auto& ray = interpolation[ray_index];
            compile_transport_row(ray, columns);
            ray.transport_value_offset = num_columns;
            num_columns += columns.size();
            if (num_columns >
                static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw std::length_error(
                    "Successive-orders transport topology exceeds the "
                    "supported index range");
            }
            row_offsets[ray_index + 1] = static_cast<int>(num_columns);
        }

        column_indices.clear();
        column_indices.reserve(num_columns);
        for (auto& ray : interpolation) {
            compile_transport_row(ray, columns);
            column_indices.insert(column_indices.end(), columns.begin(),
                                  columns.end());
        }
    }

    void SourceGeometry1D::initialize(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing,
        const SourceGeometrySettings& settings) {
        settings.validate();
        m_settings = settings;

        auto altitude_grid = make_altitude_grid();
        auto cos_sza_grid = make_cos_sza_grid(internal_viewing);
        m_location_interpolator = std::make_unique<
            sasktran2::grids::AltitudeSZASourceLocationInterpolator>(
            std::move(altitude_grid), std::move(cos_sza_grid));

        construct_source_points();
        trace_and_compile_incoming();
        compile_los_interpolation(internal_viewing);
        compile_transport_topology(m_incoming_interpolation,
                                   m_transport_row_offsets,
                                   m_transport_column_indices);
        compile_transport_topology(m_los_interpolation,
                                   m_los_transport_row_offsets,
                                   m_los_transport_column_indices);
    }

} // namespace sasktran2::successive_orders
