#include "sasktran2/geometry.h"
#include "sasktran2/viewinggeometry_internal.h"
#include <algorithm>
#include <cmath>
#ifdef SKTRAN_OPENMP_SUPPORT
#include <omp.h>
#endif
#include <sasktran2/hr/diffuse_source.h>
#include <sasktran2/math/unitsphere.h>
#include <sasktran2/math/scattering.h>
#include <sasktran2/solartransmission.h>
#include <sasktran2/do_source.h>

#include <fstream>
#include <stdexcept>

#ifdef SKTRAN_RUST_SUPPORT
namespace {
    template <typename T>
    ::rust::Slice<const T> as_rust_slice(const std::vector<T>& values) {
        return {values.data(), values.size()};
    }

    template <typename T>
    ::rust::Slice<T> as_rust_mut_slice(std::vector<T>& values) {
        return {values.data(), values.size()};
    }

    ::rust::Slice<const double> as_rust_slice(const Eigen::VectorXd& values) {
        return {values.data(), static_cast<std::size_t>(values.size())};
    }

    ::rust::Slice<const std::int32_t>
    as_rust_slice(const Eigen::VectorXi& values) {
        static_assert(sizeof(int) == sizeof(std::int32_t));
        return {reinterpret_cast<const std::int32_t*>(values.data()),
                static_cast<std::size_t>(values.size())};
    }

    ::rust::Slice<double> as_rust_mut_slice(Eigen::VectorXd& values) {
        return {values.data(), static_cast<std::size_t>(values.size())};
    }
} // namespace
#endif

namespace sasktran2::hr {

    template <int NSTOKES>
    DiffuseTable<NSTOKES>::DiffuseTable(
        const sasktran2::raytracing::RayTracerBase& ray_tracer,
        const sasktran2::Geometry1D& geometry, bool use_rust_solver)
        : m_integrator(false), m_raytracer(ray_tracer), m_geometry(geometry),
          m_altitude_grid(geometry.altitude_grid()), m_geometry_1d(&geometry),
          m_geometry_2d(nullptr),
#ifdef SKTRAN_RUST_SUPPORT
          m_raytracer_2d(nullptr),
#endif
          m_use_rust_solver(use_rust_solver) {
        m_integrator.set_on_demand_optical_depth(use_rust_solver);
    }

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    DiffuseTable<NSTOKES>::DiffuseTable(
        const sasktran2::raytracing::RustRayTracer2D& ray_tracer,
        const sasktran2::Geometry2D& geometry, bool use_rust_solver)
        : m_integrator(false), m_raytracer(ray_tracer), m_geometry(geometry),
          m_altitude_grid(geometry.altitude_grid()), m_geometry_1d(nullptr),
          m_geometry_2d(&geometry), m_raytracer_2d(&ray_tracer),
          m_use_rust_solver(use_rust_solver) {
        m_integrator.set_on_demand_optical_depth(use_rust_solver);
    }
#endif

    template <int NSTOKES>
    sasktran2::grids::Grid
    DiffuseTable<NSTOKES>::generate_cos_sza_grid(double min_cos_sza,
                                                 double max_cos_sza) {
        Eigen::VectorXd cos_sza_grid_values;

        // TODO: should we have separate SZA spacings for DO/HR? Or harmonize
        // the naming?
        if (m_config->num_do_sza() > 1) {
            cos_sza_grid_values.setLinSpaced(m_config->num_do_sza(),
                                             min_cos_sza, max_cos_sza);
        } else {
            // Set to reference cos_sza
            cos_sza_grid_values.resize(1);
            cos_sza_grid_values.setConstant(
                m_geometry.coordinates().cos_sza_at_reference());
        }

        return sasktran2::grids::Grid(std::move(cos_sza_grid_values),
                                      sasktran2::grids::gridspacing::constant,
                                      sasktran2::grids::outofbounds::extend,
                                      sasktran2::grids::interpolation::linear);
    }

    template <int NSTOKES>
    sasktran2::grids::AltitudeGrid
    DiffuseTable<NSTOKES>::generate_altitude_grid() {
        const auto& configured_grid =
            m_config->successive_orders_altitude_grid_m();
        Eigen::VectorXd alt_values;
        auto spacing = sasktran2::grids::gridspacing::constant;

        if (configured_grid.empty()) {
            // Preserve the historical default: one source point at the
            // midpoint of each atmospheric layer.
            alt_values = (m_altitude_grid.grid()(
                              Eigen::seq(0, Eigen::placeholders::last - 1)) +
                          m_altitude_grid.grid()(
                              Eigen::seq(1, Eigen::placeholders::last))) /
                         2.0;
        } else {
            const double minimum_altitude = m_altitude_grid.grid()[0];
            const double maximum_altitude =
                m_altitude_grid.grid()[m_altitude_grid.grid().size() - 1];
            alt_values.resize(
                static_cast<Eigen::Index>(configured_grid.size()));
            for (std::size_t altitude_index = 0;
                 altitude_index < configured_grid.size(); ++altitude_index) {
                const double altitude = configured_grid[altitude_index];
                if (!std::isfinite(altitude) ||
                    (altitude_index > 0 &&
                     altitude <= configured_grid[altitude_index - 1])) {
                    throw std::invalid_argument(
                        "successive_orders_altitude_grid_m must contain "
                        "finite, "
                        "strictly increasing altitudes");
                }
                if (altitude < minimum_altitude ||
                    altitude > maximum_altitude) {
                    throw std::invalid_argument(
                        "successive_orders_altitude_grid_m must lie within the "
                        "atmospheric altitude grid");
                }
                alt_values[static_cast<Eigen::Index>(altitude_index)] =
                    altitude;
            }
            spacing = sasktran2::grids::gridspacing::automatic;
        }

        return sasktran2::grids::AltitudeGrid(
            std::move(alt_values), spacing,
            sasktran2::grids::outofbounds::extend,
            sasktran2::grids::interpolation::linear);
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::construct_diffuse_points() {
        // Create the sphere pairs
        // One for the atmosphere points and one for each of the ground points
        // for now

        int num_interior_spheres = 1;
        m_unit_sphere_pairs.resize(
            num_interior_spheres +
            m_location_interpolator->num_ground_points());

        // TODO: Get npoints from config, figure out if we want to use more than
        // one type of sphere
        m_unit_sphere_pairs[0] = std::make_unique<
            sasktran2::hr::IncomingOutgoingSpherePair<NSTOKES>>(
            m_config->num_do_streams(),
            std::move(std::make_unique<sasktran2::math::LebedevSphere>(
                m_config->num_hr_incoming())),
            std::move(std::make_unique<sasktran2::math::LebedevSphere>(
                m_config->num_hr_outgoing())),
            !m_use_rust_solver);

        // TODO: Same number of streams for the ground term? probably... to be
        // figured out when BRDF is implemented
        for (int i = 0; i < m_location_interpolator->num_ground_points(); ++i) {
            Eigen::Vector3d location = m_location_interpolator->ground_location(
                m_geometry.coordinates(), i);

            m_unit_sphere_pairs[num_interior_spheres + i] = std::make_unique<
                sasktran2::hr::IncomingOutgoingSpherePair<NSTOKES>>(
                m_config->num_do_streams(),
                std::make_unique<sasktran2::math::UnitSphereGround>(
                    std::move(std::make_unique<sasktran2::math::LebedevSphere>(
                        m_config->num_hr_incoming())),
                    location),
                std::make_unique<sasktran2::math::UnitSphereGround>(
                    std::move(std::make_unique<sasktran2::math::LebedevSphere>(
                        m_config->num_hr_outgoing())),
                    location),
                !m_use_rust_solver);
        }

        m_diffuse_points.resize(m_location_interpolator->num_interior_points() +
                                m_location_interpolator->num_ground_points());

        if (m_config->num_hr_full_incoming_points() > 0) {
            // We are approximating the multiple scatter source by only
            // calculating incoming quantities at a subset of the diffuse points

            // Start by setting the calculation to false at all points
            m_diffuse_point_full_calculation.resize(m_diffuse_points.size(),
                                                    false);

            // At all ground points we will do the incoming calculation
            for (int i = m_location_interpolator->num_interior_points();
                 i < m_location_interpolator->num_interior_points() +
                         m_location_interpolator->num_ground_points();
                 ++i) {
                m_diffuse_point_full_calculation[i] = true;
            }

            int num_inc_per_profile = m_config->num_hr_full_incoming_points();
            // TODO: This implicitly assumes the construction of the diffuse
            // point locations, should figure out a better way to do this
            int num_diffuse_in_profile =
                (m_location_interpolator->num_interior_points()) /
                m_config->num_do_sza();
            if (num_inc_per_profile > num_diffuse_in_profile) {
                throw std::invalid_argument(
                    "num_successive_order_points cannot exceed the number of "
                    "successive-orders source altitudes");
            }

            for (int i = 0; i < m_config->num_do_sza(); ++i) {
                // Start of the profile index
                int profile_start = i * num_diffuse_in_profile;

                for (int j = 0; j < num_inc_per_profile; ++j) {
                    // Basically want linearly spaced from 0 to end inclusive
                    int alt_index = (j * (num_diffuse_in_profile - 1)) /
                                    (num_inc_per_profile - 1);

                    m_diffuse_point_full_calculation[profile_start +
                                                     alt_index] = true;
                }
            }
        } else {
            // We are calculating incoming quantities at all points
            m_diffuse_point_full_calculation.resize(m_diffuse_points.size(),
                                                    true);
        }

        sasktran2::Location loc;

        for (int i = 0; i < m_location_interpolator->num_interior_points();
             ++i) {
            auto& point = m_diffuse_points[i];

            loc.position = m_location_interpolator->grid_location(
                m_geometry.coordinates(), i);

            point = std::make_unique<sasktran2::hr::DiffusePoint<NSTOKES>>(
                *m_unit_sphere_pairs[0], loc);
        }

        for (int i = 0; i < m_location_interpolator->num_ground_points(); ++i) {
            auto& point = m_diffuse_points[i + m_location_interpolator
                                                   ->num_interior_points()];

            loc.position = m_location_interpolator->ground_location(
                m_geometry.coordinates(), i);

            // Add 0.01m to the ground location to avoid rounding errors
            loc.position += loc.position.normalized() * 0.01;

            point = std::make_unique<sasktran2::hr::DiffusePoint<NSTOKES>>(
                *m_unit_sphere_pairs[i + num_interior_spheres], loc);
        }

        // Construct interpolators to the diffuse point locations
        m_diffuse_point_interpolation_weights.resize(m_diffuse_points.size());
        for (int i = 0; i < m_diffuse_points.size(); ++i) {
            auto& point = m_diffuse_points[i];

            m_geometry.assign_interpolation_weights(
                point->location(), m_diffuse_point_interpolation_weights[i]);
        }

        m_diffuse_incoming_index_map.resize(m_diffuse_points.size());
        m_diffuse_outgoing_index_map.resize(m_diffuse_points.size());

        int start_incoming_idx = 0;
        int start_outgoing_idx = 0;
        for (int i = 0; i < m_diffuse_points.size(); ++i) {
            m_diffuse_incoming_index_map[i] = start_incoming_idx;
            m_diffuse_outgoing_index_map[i] = start_outgoing_idx;
            if (m_diffuse_point_full_calculation[i]) {
                start_incoming_idx += m_diffuse_points[i]->num_incoming();
            }
            start_outgoing_idx += m_diffuse_points[i]->num_outgoing();
        }
        m_num_outgoing_values =
            static_cast<std::size_t>(start_outgoing_idx) * NSTOKES;

        m_internal_viewing.traced_rays.resize(start_incoming_idx);

        for (auto& storage : m_thread_storage) {
            storage.m_incoming_radiances.resize(start_incoming_idx * NSTOKES, 0,
                                                false);
            storage.m_firstorder_radiances.resize(start_incoming_idx * NSTOKES,
                                                  0, false);
            storage.m_outgoing_sources.resize(start_outgoing_idx * NSTOKES, 0,
                                              false);
            if (m_use_rust_solver && m_wavelength_block_capacity > 1) {
                storage.rust_batch_outgoing_sources.resize(
                    start_outgoing_idx * NSTOKES * m_wavelength_block_capacity);
            }

            storage.point_scattering_matrices.resize(m_diffuse_points.size());
            for (int i = 0; i < m_diffuse_points.size(); ++i) {
                if (m_diffuse_point_full_calculation[i] && !m_use_rust_solver) {
                    storage.point_scattering_matrices[i].resize(
                        m_diffuse_points[i]->num_outgoing() * NSTOKES,
                        m_diffuse_points[i]->num_incoming() * NSTOKES);
                }
            }
        }
    }

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::initialize_twostream_initial_guess_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        if (m_config->successive_orders_initialization() !=
            sasktran2::Config::SuccessiveOrdersInitialization::twostream) {
            return;
        }
        if (!m_use_rust_solver) {
            throw std::invalid_argument(
                "Two-stream successive-orders initialization requires the "
                "Rust solver");
        }
        (void)internal_viewing;

        const int num_profiles = m_location_interpolator->num_ground_points();
        const int num_interior = m_location_interpolator->num_interior_points();
        if (num_profiles < 1 || num_interior % num_profiles != 0) {
            throw std::logic_error(
                "Successive-orders source profiles have an invalid shape");
        }
        m_twostream_num_source_altitudes = num_interior / num_profiles;

        const auto& atmospheric_altitudes = m_altitude_grid.grid();
        const std::size_t num_levels =
            static_cast<std::size_t>(atmospheric_altitudes.size());
        const int num_layers = static_cast<int>(num_levels) - 1;
        const auto geometry_type = m_geometry.coordinates().geometry_type();
        const double earth_radius = m_geometry.coordinates().earth_radius();

        m_twostream_atmosphere_weights.assign(
            num_profiles,
            std::vector<std::vector<std::pair<int, double>>>(num_levels));
        m_twostream_initializers.clear();
        m_twostream_initializers.resize(m_thread_storage.size());
        if (m_twostream_initializers.empty()) {
            throw std::logic_error(
                "Two-stream initialization requires a wavelength worker");
        }
        for (auto& worker : m_twostream_initializers) {
            worker.reserve(num_profiles);
        }

        for (int profile = 0; profile < num_profiles; ++profile) {
            const int first_point = profile * m_twostream_num_source_altitudes;
            const auto& profile_location =
                m_diffuse_points[first_point]->location();
            const double cos_sza =
                m_geometry.coordinates()
                    .solar_angles_at_location(profile_location.position)
                    .first;

            Eigen::VectorXd altitude_copy = atmospheric_altitudes;
            sasktran2::Geometry1D column_geometry(
                cos_sza, 0.0, earth_radius, std::move(altitude_copy),
                sasktran2::grids::interpolation::linear, geometry_type);
            if (m_geometry_1d != nullptr &&
                m_geometry_1d->refractive_index().size() ==
                    column_geometry.refractive_index().size()) {
                column_geometry.refractive_index() =
                    m_geometry_1d->refractive_index();
            }

            std::vector<double> layer_thickness(num_layers);
            std::vector<double> chapman_factors(
                static_cast<std::size_t>(num_layers) * num_layers);
            for (int layer = 0; layer < num_layers; ++layer) {
                layer_thickness[layer] =
                    atmospheric_altitudes[num_layers - layer] -
                    atmospheric_altitudes[num_layers - layer - 1];
            }

            if (geometry_type == sasktran2::geometrytype::planeparallel) {
                const double factor = cos_sza > 0 ? 1.0 / cos_sza : -1.0;
                for (int boundary = 0; boundary < num_layers; ++boundary) {
                    for (int layer = 0; layer <= boundary; ++layer) {
                        chapman_factors[boundary * num_layers + layer] = factor;
                    }
                }
            } else {
                sasktran2::raytracing::SphericalShellRayTracer solar_tracer(
                    column_geometry);
                sasktran2::viewinggeometry::ViewingRay solar_ray;
                solar_ray.look_away = column_geometry.coordinates().sun_unit();
                sasktran2::raytracing::TracedRay traced_solar_ray;
                for (int boundary = 0; boundary < num_layers; ++boundary) {
                    const double floor_altitude =
                        atmospheric_altitudes[num_layers - boundary - 1];
                    solar_ray.observer.position =
                        column_geometry.coordinates().solar_coordinate_vector(
                            cos_sza, 0.0, floor_altitude);
                    solar_tracer.trace_ray(solar_ray, traced_solar_ray,
                                           m_config->solar_refraction());
                    if (traced_solar_ray.ground_is_hit) {
                        std::fill(chapman_factors.begin() +
                                      boundary * num_layers,
                                  chapman_factors.begin() +
                                      (boundary + 1) * num_layers,
                                  -1.0);
                        continue;
                    }
                    for (const auto& traced_layer : traced_solar_ray.layers) {
                        const double average_altitude =
                            0.5 * (traced_layer.entrance.radius() +
                                   traced_layer.exit.radius()) -
                            earth_radius;
                        const double* upper = std::upper_bound(
                            atmospheric_altitudes.data(),
                            atmospheric_altitudes.data() + num_levels,
                            average_altitude);
                        int bottom_layer =
                            upper == atmospheric_altitudes.data()
                                ? 0
                                : static_cast<int>(
                                      upper - atmospheric_altitudes.data() - 1);
                        bottom_layer =
                            std::clamp(bottom_layer, 0, num_layers - 1);
                        const int top_down_layer =
                            num_layers - 1 - bottom_layer;
                        chapman_factors[boundary * num_layers +
                                        top_down_layer] +=
                            traced_layer.layer_distance *
                            traced_layer.curvature_factor /
                            layer_thickness[top_down_layer];
                    }
                }
            }

            std::vector<std::size_t> sample_layers;
            std::vector<double> sample_fractions;
            std::vector<double> sample_cosines;
            std::vector<double> sample_azimuths;
            const int num_outgoing =
                m_diffuse_points[first_point]->num_outgoing();
            const std::size_t num_samples =
                static_cast<std::size_t>(m_twostream_num_source_altitudes) *
                num_outgoing;
            sample_layers.reserve(num_samples);
            sample_fractions.reserve(num_samples);
            sample_cosines.reserve(num_samples);
            sample_azimuths.reserve(num_samples);

            for (int altitude_index = 0;
                 altitude_index < m_twostream_num_source_altitudes;
                 ++altitude_index) {
                const auto& point =
                    *m_diffuse_points[first_point + altitude_index];
                const auto& location = point.location();
                const double altitude =
                    geometry_type == sasktran2::geometrytype::spherical
                        ? location.radius() - earth_radius
                        : location.position.z() - earth_radius;
                const double* upper = std::upper_bound(
                    atmospheric_altitudes.data(),
                    atmospheric_altitudes.data() + num_levels, altitude);
                int bottom_layer =
                    upper == atmospheric_altitudes.data()
                        ? 0
                        : static_cast<int>(upper -
                                           atmospheric_altitudes.data() - 1);
                bottom_layer = std::clamp(bottom_layer, 0, num_layers - 1);
                const double bottom = atmospheric_altitudes[bottom_layer];
                const double top = atmospheric_altitudes[bottom_layer + 1];
                const double fraction_from_top =
                    std::clamp((top - altitude) / (top - bottom), 0.0, 1.0);
                const Eigen::Vector3d local_up =
                    geometry_type == sasktran2::geometrytype::spherical
                        ? location.position.normalized()
                        : Eigen::Vector3d(0, 0, 1);

                for (int direction_index = 0;
                     direction_index < point.num_outgoing();
                     ++direction_index) {
                    const Eigen::Vector3d direction =
                        point.sphere_pair().outgoing_sphere().get_quad_position(
                            direction_index);
                    double sample_cos_sza;
                    double relative_azimuth;
                    sasktran2::raytracing::calculate_csz_saz(
                        m_geometry.coordinates().sun_unit(), location,
                        direction, sample_cos_sza, relative_azimuth,
                        geometry_type);
                    sample_layers.push_back(static_cast<std::size_t>(
                        num_layers - 1 - bottom_layer));
                    sample_fractions.push_back(fraction_from_top);
                    sample_cosines.push_back(
                        std::clamp(-direction.dot(local_up), -1.0, 1.0));
                    sample_azimuths.push_back(relative_azimuth);
                }
            }

            auto source = sasktran2::rust::twostream::new_rust_twostream_source(
                as_rust_slice(layer_thickness), as_rust_slice(chapman_factors),
                cos_sza, 0, 1);
            sasktran2::rust::twostream::set_local_source_geometry(
                *source, as_rust_slice(sample_layers),
                as_rust_slice(sample_fractions), as_rust_slice(sample_cosines),
                as_rust_slice(sample_azimuths));
            m_twostream_initializers[0].push_back(std::move(source));
            for (std::size_t worker = 1;
                 worker < m_twostream_initializers.size(); ++worker) {
                m_twostream_initializers[worker].push_back(
                    sasktran2::rust::twostream::clone_rust_twostream_source(
                        *m_twostream_initializers[0].back()));
            }

            const Eigen::Vector3d horizontal_unit =
                profile_location.position.normalized();
            for (std::size_t level = 0; level < num_levels; ++level) {
                if (m_geometry_1d != nullptr) {
                    m_twostream_atmosphere_weights[profile][level] = {
                        {static_cast<int>(level), 1.0}};
                } else {
                    sasktran2::Location location;
                    location.position =
                        (earth_radius + atmospheric_altitudes[level]) *
                        horizontal_unit;
                    m_geometry.assign_interpolation_weights(
                        location,
                        m_twostream_atmosphere_weights[profile][level]);
                }
            }
        }
    }
#endif

    template <int NSTOKES>
    std::vector<std::vector<int>> DiffuseTable<NSTOKES>::trace_incoming_rays() {
        ZoneScopedN("Trace and Compile Incoming Rays");
        int nthreads = m_config->num_threads();

        std::vector<sasktran2::viewinggeometry::ViewingRay> thread_viewing_ray;
        thread_viewing_ray.resize(nthreads);
        std::vector<std::vector<std::pair<int, double>>>
            thread_temp_location_storage(nthreads);
        std::vector<std::vector<std::pair<int, double>>>
            thread_temp_direction_storage(nthreads);
        std::vector<std::vector<std::pair<int, double>>>
            thread_temp_atmosphere_storage(nthreads);
        std::vector<std::vector<std::pair<int, std::uint16_t*>>>
            thread_sorting_helper(nthreads);
        std::vector<std::vector<int>> transport_columns(
            m_internal_viewing.traced_rays.size());
        m_diffuse_source_weights.resize(m_internal_viewing.traced_rays.size());

#pragma omp parallel for num_threads(nthreads)
        for (int i = 0; i < m_diffuse_points.size(); ++i) {
#ifdef SKTRAN_OPENMP_SUPPORT
            const int threadidx = omp_get_thread_num();
#else
            const int threadidx = 0;
#endif
            auto& viewing_ray = thread_viewing_ray[threadidx];
            if (!m_diffuse_point_full_calculation[i]) {
                continue;
            }
            viewing_ray.observer = m_diffuse_points[i]->location();
            for (int j = 0; j < m_diffuse_points[i]->num_incoming(); ++j) {
                viewing_ray.look_away =
                    m_diffuse_points[i]->incoming_direction(j);

                const int rayidx = m_diffuse_incoming_index_map[i] + j;
                auto& traced_ray = m_internal_viewing.traced_rays[rayidx];
                m_raytracer.trace_ray(viewing_ray, traced_ray,
                                      m_config->multiple_scatter_refraction());
                generate_source_interpolation_weights(
                    traced_ray, m_diffuse_source_weights[rayidx],
                    thread_temp_location_storage[threadidx],
                    thread_temp_direction_storage[threadidx],
                    thread_temp_atmosphere_storage[threadidx]);
                compile_accumulation_row(m_diffuse_source_weights[rayidx],
                                         transport_columns[rayidx],
                                         thread_sorting_helper[threadidx]);
                if (m_use_rust_solver && m_geometry_2d != nullptr &&
                    (NSTOKES == 1 ||
                     m_altitude_grid.interpolation_method() !=
                         sasktran2::grids::interpolation::lower)) {
                    m_integrator.compact_geometry_2d_ray(rayidx, traced_ray);
                }
            }
        }
        return transport_columns;
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        ZoneScopedN("Initialize HR Geometry");
        if (m_geometry_2d != nullptr) {
            Eigen::VectorXd horizontal_source_grid;
            if (m_config->num_do_sza() > 1) {
                horizontal_source_grid.setLinSpaced(
                    m_config->num_do_sza(),
                    m_geometry_2d->horizontal_angle_grid()[0],
                    m_geometry_2d->horizontal_angle_grid()
                        [m_geometry_2d->num_horizontal_locations() - 1]);
            } else {
                horizontal_source_grid.resize(1);
                horizontal_source_grid[0] =
                    0.5 * (m_geometry_2d->horizontal_angle_grid()[0] +
                           m_geometry_2d->horizontal_angle_grid()
                               [m_geometry_2d->num_horizontal_locations() - 1]);
            }
            m_location_interpolator = std::make_unique<
                sasktran2::grids::AltitudeHorizontalSourceLocationInterpolator>(
                generate_altitude_grid(), std::move(horizontal_source_grid),
                *m_geometry_2d);
        } else {
            const std::pair<double, double> min_max_cos_sza =
                sasktran2::raytracing::min_max_cos_sza_of_all_rays(
                    internal_viewing.traced_rays);
            m_location_interpolator = std::make_unique<
                sasktran2::grids::AltitudeSZASourceLocationInterpolator>(
                generate_altitude_grid(),
                generate_cos_sza_grid(min_max_cos_sza.first,
                                      min_max_cos_sza.second));
        }
        // Construct the actual diffuse points
        construct_diffuse_points();

#ifdef SKTRAN_RUST_SUPPORT
        initialize_twostream_initial_guess_geometry(internal_viewing);
#endif

        const bool use_incremental_compact_geometry =
            m_use_rust_solver && m_geometry_2d != nullptr &&
            (NSTOKES == 1 || m_altitude_grid.interpolation_method() !=
                                 sasktran2::grids::interpolation::lower);
        if (use_incremental_compact_geometry) {
            m_integrator.begin_compact_geometry_2d(
                m_internal_viewing.traced_rays, m_geometry);
        }

        // Trace all of the incoming rays
        auto transport_columns = trace_incoming_rays();

        // Set up the integrator
        if (use_incremental_compact_geometry) {
            m_integrator.finalize_compact_geometry_2d();
        } else {
            m_integrator.initialize_geometry(m_internal_viewing.traced_rays,
                                             this->m_geometry);
        }
        // And the initial sources
        // This is a little tricky, any source that is used for the incoming
        // rays needs to be initialized with the traced incoming rays

        if (use_incremental_compact_geometry) {
            std::vector<std::uint32_t> layer_counts(
                m_diffuse_source_weights.size());
            for (std::size_t rayidx = 0;
                 rayidx < m_diffuse_source_weights.size(); ++rayidx) {
                layer_counts[rayidx] = static_cast<std::uint32_t>(
                    m_diffuse_source_weights[rayidx].interior_weights.size());
            }
            auto* single_scatter =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionExact,
                    NSTOKES>*>(m_initial_sources.front());
            if (single_scatter == nullptr || m_initial_sources.size() != 1) {
                throw std::logic_error(
                    "Compact 2D successive orders requires one exact "
                    "single-scatter source");
            }
            single_scatter->initialize_geometry_compact(m_internal_viewing,
                                                        layer_counts);
        } else {
            for (auto& source : m_initial_sources) {
                source->initialize_geometry(m_internal_viewing);
            }
        }
        // But the DO Source should be initialized with the LOS rays
        if (m_config->initialize_hr_with_do()) {
            m_do_source->initialize_geometry(internal_viewing);
        }

        construct_accumulation_sparsity(transport_columns);
        std::vector<std::vector<int>>().swap(transport_columns);
#ifdef SKTRAN_RUST_SUPPORT
        release_cpp_transport_geometry();
#endif
        generate_source_interpolation_weights(internal_viewing.traced_rays,
                                              m_los_source_weights);
#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            initialize_rust_source_interpolator(internal_viewing);
        }
#endif

        if (m_config->initialize_hr_with_do()) {
            // Have to create a vector of all locations and directions
            std::vector<Eigen::Vector3d> locations, directions;
            std::vector<bool> ground_point;

            for (int i = 0; i < m_diffuse_points.size(); ++i) {
                const auto& point = m_diffuse_points[i];

                for (int j = 0; j < point->num_outgoing(); ++j) {
                    locations.push_back(point->location().position);
                    directions.push_back(point->sphere_pair()
                                             .outgoing_sphere()
                                             .get_quad_position(j));

                    if (i < m_location_interpolator->num_interior_points()) {
                        ground_point.push_back(false);
                    } else {
                        ground_point.push_back(true);
                    }
                }
            }

            m_do_source->storage().create_location_source_interpolator(
                locations, directions, ground_point,
                m_do_to_diffuse_outgoing_interpolator);
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_atmosphere = &atmosphere;

        m_integrator.initialize_atmosphere(atmosphere);

        for (auto& source : m_initial_owned_sources) {
            source->initialize_atmosphere(atmosphere);
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::initialize_atmosphere_native(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_atmosphere = &atmosphere;
        m_integrator.initialize_atmosphere(atmosphere);
        for (auto& source : m_initial_owned_sources) {
            source->initialize_atmosphere_native(atmosphere);
        }
    }

    template <int NSTOKES>
    void
    DiffuseTable<NSTOKES>::set_wavelength_block_capacity(int block_capacity) {
        if (block_capacity < 1 ||
            block_capacity > maximum_wavelength_block_size()) {
            throw std::invalid_argument(
                "Invalid successive-orders wavelength block capacity");
        }
        m_wavelength_block_capacity = block_capacity;

        if (m_use_rust_solver) {
            for (auto& storage : m_thread_storage) {
                if (block_capacity == 1) {
                    std::vector<double>().swap(
                        storage.rust_batch_outgoing_sources);
                } else {
                    storage.rust_batch_outgoing_sources.resize(
                        m_num_outgoing_values * block_capacity);
                }
            }
        }
    }

    template <int NSTOKES>
    void
    DiffuseTable<NSTOKES>::begin_forward_state_capture(int num_wavelengths) {
        if (!m_use_rust_solver || num_wavelengths < 1) {
            return;
        }
        const bool compatible_existing_state =
            m_forward_state_atmosphere == m_atmosphere &&
            m_forward_state_atmosphere_revision == m_atmosphere->revision() &&
            m_forward_states.size() ==
                static_cast<std::size_t>(num_wavelengths);
        if (!compatible_existing_state) {
            m_forward_states.resize(static_cast<std::size_t>(num_wavelengths));
            for (auto& state : m_forward_states) {
                state.valid = false;
            }
        }
        m_forward_state_atmosphere = m_atmosphere;
        m_forward_state_atmosphere_revision = m_atmosphere->revision();
        m_capture_forward_state = true;
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::end_forward_state_capture() {
        m_capture_forward_state = false;
    }

    template <int NSTOKES>
    std::size_t DiffuseTable<NSTOKES>::retained_forward_state_bytes() const {
        std::size_t bytes =
            m_forward_states.capacity() * sizeof(DiffuseTableForwardState);
        for (const auto& state : m_forward_states) {
            bytes += state.first_order_forcing.capacity() * sizeof(double);
            bytes += state.solution.capacity() * sizeof(double);
        }
        return bytes;
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::construct_accumulation_sparsity(
        const std::vector<std::vector<int>>& transport_columns) {
        ZoneScopedN("Construct Accumulation Sparsity");
        if (transport_columns.size() != m_diffuse_source_weights.size()) {
            throw std::invalid_argument(
                "Successive-orders transport row count mismatch");
        }

        std::uint64_t total_nnz_64 = 0;
        for (const auto& columns : transport_columns) {
            total_nnz_64 +=
                static_cast<std::uint64_t>(columns.size()) * NSTOKES;
        }
        if (total_nnz_64 >
            static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
            throw std::length_error(
                "Successive-orders transport storage exceeds Eigen index "
                "range");
        }
        const int total_nnz = static_cast<int>(total_nnz_64);
        m_inner_indicies.resize(total_nnz);
        m_outer_starts.resize(NSTOKES * m_diffuse_source_weights.size() + 1);
        m_outer_starts(0) = 0;

        int current_index = 0;
        for (std::size_t row = 0; row < transport_columns.size(); ++row) {
            auto& weights = m_diffuse_source_weights[row];
            const auto& columns = transport_columns[row];
            weights.accumulation_row_offset =
                static_cast<std::uint32_t>(current_index);
            weights.accumulation_row_nnz =
                static_cast<std::uint32_t>(columns.size());

            for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                for (const int column : columns) {
                    m_inner_indicies(current_index++) =
                        column * NSTOKES + stokes;
                }
                m_outer_starts(row * NSTOKES + stokes + 1) = current_index;
            }
        }

        // The legacy and lower-interpolation paths assemble into C++ storage.
        // The normal Rust path owns its wavelength-local CSR values inside the
        // per-thread solver and only retains the shared sparsity pattern here.
        const bool retain_cpp_transport_values =
            !m_use_rust_solver || m_altitude_grid.interpolation_method() ==
                                      sasktran2::grids::interpolation::lower;
        for (auto& storage : m_thread_storage) {
            if (retain_cpp_transport_values) {
                storage.accumulation_summed_values.resize(total_nnz);
            } else {
                storage.accumulation_summed_values.resize(0);
            }
        }

#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            const bool rust_owned_transport =
                m_altitude_grid.interpolation_method() !=
                sasktran2::grids::interpolation::lower;
            m_rust_transport_sparsity.clear();
            if (rust_owned_transport) {
                m_rust_transport_sparsity.push_back(
                    sasktran2::rust::successive_orders::new_transport_sparsity(
                        static_cast<std::size_t>(m_outer_starts.size() - 1),
                        m_num_outgoing_values, as_rust_slice(m_outer_starts),
                        as_rust_slice(m_inner_indicies)));
                spdlog::debug(
                    "Packed Rust successive-orders transport sparsity: {} "
                    "rows, {} values, {:.3f} MiB shared across wavelength "
                    "workers",
                    m_outer_starts.size() - 1, m_inner_indicies.size(),
                    static_cast<double>(
                        sasktran2::rust::successive_orders::
                            transport_sparsity_storage_bytes(
                                *m_rust_transport_sparsity.front())) /
                        (1024.0 * 1024.0));
            }
            std::vector<std::uint32_t> scattering_point_offsets = {0};
            std::vector<std::uint32_t> scattering_atmosphere_indices;
            std::vector<double> scattering_interpolation_weights;
            const int num_interior_points =
                m_location_interpolator->num_interior_points();
            scattering_point_offsets.reserve(
                static_cast<std::size_t>(num_interior_points) + 1);
            for (int point = 0; point < num_interior_points; ++point) {
                for (const auto& [atmosphere_index, weight] :
                     m_diffuse_point_interpolation_weights[point]) {
                    if (atmosphere_index < 0) {
                        throw std::logic_error(
                            "Negative atmosphere index in successive-orders "
                            "scattering interpolation");
                    }
                    scattering_atmosphere_indices.push_back(
                        static_cast<std::uint32_t>(atmosphere_index));
                    scattering_interpolation_weights.push_back(weight);
                }
                if (scattering_atmosphere_indices.size() >
                    std::numeric_limits<std::uint32_t>::max()) {
                    throw std::length_error(
                        "Successive-orders scattering interpolation exceeds "
                        "the packed offset range");
                }
                scattering_point_offsets.push_back(static_cast<std::uint32_t>(
                    scattering_atmosphere_indices.size()));
            }
            const std::size_t coefficient_families = NSTOKES == 1 ? 1 : 4;
            const std::size_t output_coefficients =
                static_cast<std::size_t>(m_config->num_do_streams()) *
                coefficient_families;
            m_rust_scattering_interpolators.clear();
            m_rust_scattering_interpolators.push_back(
                sasktran2::rust::successive_orders::
                    new_scattering_coefficient_interpolator(
                        as_rust_slice(scattering_point_offsets),
                        as_rust_slice(scattering_atmosphere_indices),
                        as_rust_slice(scattering_interpolation_weights),
                        static_cast<std::size_t>(m_geometry.size()),
                        output_coefficients));
            spdlog::debug(
                "Packed Rust successive-orders scattering interpolation: {} "
                "points, {} weights, {:.3f} MiB immutable",
                num_interior_points,
                sasktran2::rust::successive_orders::
                    scattering_coefficient_interpolator_num_weights(
                        *m_rust_scattering_interpolators.front()),
                static_cast<double>(
                    sasktran2::rust::successive_orders::
                        scattering_coefficient_interpolator_storage_bytes(
                            *m_rust_scattering_interpolators.front())) /
                    (1024.0 * 1024.0));
            for (int point = 0; point < num_interior_points; ++point) {
                std::vector<std::pair<int, double>>().swap(
                    m_diffuse_point_interpolation_weights[point]);
            }

            {
                auto packed_geometry = m_integrator.pack_ray_transport_geometry(
                    m_diffuse_source_weights, m_geometry);
                m_rust_ray_transports.clear();
                m_rust_ray_transports.push_back(
                    sasktran2::rust::successive_orders::new_ray_transport(
                        as_rust_slice(packed_geometry.ray_layer_offsets),
                        as_rust_slice(packed_geometry.layer_atmosphere_offsets),
                        as_rust_slice(packed_geometry.atmosphere_indices),
                        as_rust_slice(packed_geometry.optical_depth_weights),
                        as_rust_slice(packed_geometry.albedo_weights),
                        as_rust_slice(packed_geometry.entrance_weights),
                        as_rust_slice(packed_geometry.exit_weights),
                        as_rust_slice(packed_geometry.layer_distance),
                        as_rust_slice(packed_geometry.layer_start_fraction),
                        as_rust_slice(packed_geometry.layer_end_fraction),
                        as_rust_slice(packed_geometry.ray_scattering_cosine),
                        as_rust_slice(packed_geometry.ray_phase_q_factor),
                        as_rust_slice(packed_geometry.ray_phase_u_factor),
                        static_cast<std::size_t>(
                            m_config->num_singlescatter_moments()),
                        m_geometry_2d != nullptr,
                        as_rust_slice(packed_geometry.layer_source_offsets),
                        as_rust_slice(
                            packed_geometry.ray_transport_value_offsets),
                        as_rust_slice(packed_geometry.ray_transport_row_nnz),
                        static_cast<std::size_t>(NSTOKES),
                        as_rust_slice(
                            packed_geometry.source_value_inner_indices),
                        as_rust_slice(packed_geometry.source_weights),
                        as_rust_slice(packed_geometry.ray_ground_offsets),
                        as_rust_slice(
                            packed_geometry.ground_value_inner_indices),
                        as_rust_slice(packed_geometry.ground_weights)));

                m_rust_transport_ray_layer_offsets =
                    std::move(packed_geometry.ray_layer_offsets);
                const std::size_t num_layers = sasktran2::rust::
                    successive_orders::ray_transport_num_layers(
                        *m_rust_ray_transports.front());
                const std::size_t num_rays =
                    sasktran2::rust::successive_orders::ray_transport_num_rays(
                        *m_rust_ray_transports.front());
                const std::size_t geometry_bytes = sasktran2::rust::
                    successive_orders::ray_transport_storage_bytes(
                        *m_rust_ray_transports.front());
                const std::size_t atmosphere_weights = sasktran2::rust::
                    successive_orders::ray_transport_num_atmosphere_weights(
                        *m_rust_ray_transports.front());
                const std::size_t source_weights = sasktran2::rust::
                    successive_orders::ray_transport_num_source_weights(
                        *m_rust_ray_transports.front());
                const std::size_t ground_weights = sasktran2::rust::
                    successive_orders::ray_transport_num_ground_weights(
                        *m_rust_ray_transports.front());
                const std::size_t phase_geometries = sasktran2::rust::
                    successive_orders::ray_transport_num_phase_geometries(
                        *m_rust_ray_transports.front());
                const std::size_t grouped_phase_geometries =
                    sasktran2::rust::successive_orders::
                        ray_transport_num_grouped_phase_geometries(
                            *m_rust_ray_transports.front(), m_geometry.size());
                spdlog::debug(
                    "Packed {}-Stokes Rust successive-orders geometry: {} "
                    "rays, "
                    "{} phase geometries ({} grouped, {} direct), {} layers, "
                    "{} atmosphere and {}+{} transport weights, "
                    "{:.3f} MiB immutable + up to {:.3f} MiB lazy VJP "
                    "attenuation scratch across {} wavelength workers",
                    NSTOKES, num_rays, phase_geometries,
                    grouped_phase_geometries,
                    phase_geometries - grouped_phase_geometries, num_layers,
                    atmosphere_weights, source_weights, ground_weights,
                    static_cast<double>(geometry_bytes) / (1024.0 * 1024.0),
                    static_cast<double>(
                        m_thread_storage.size() *
                        ((NSTOKES == 1 ? 3 : 4) * num_layers + 3 * num_rays) *
                        sizeof(double)) /
                        (1024.0 * 1024.0),
                    m_thread_storage.size());
                if (num_rays != m_diffuse_source_weights.size() ||
                    m_rust_transport_ray_layer_offsets.size() != num_rays + 1) {
                    throw std::logic_error(
                        "Rust transport geometry size mismatch");
                }
                for (auto& storage : m_thread_storage) {
                    if (m_altitude_grid.interpolation_method() ==
                        sasktran2::grids::interpolation::lower) {
                        storage.rust_layer_optical_depth.resize(num_layers);
                        storage.rust_layer_attenuation.resize(num_layers);
                        storage.rust_layer_prefix_attenuation.resize(
                            num_layers);
                        storage.rust_ray_end_attenuation.resize(num_rays);
                    }
                    storage.rust_end_of_ray_source.resize(num_rays * NSTOKES);
                    storage.rust_end_of_ray_source_tangent.resize(num_rays *
                                                                  NSTOKES);
                    storage.rust_end_of_ray_source_gradient.resize(num_rays *
                                                                   NSTOKES);
                }
            }

            std::vector<std::size_t> scattering_output_offsets = {0};
            std::vector<std::size_t> scattering_input_offsets = {0};
            scattering_output_offsets.reserve(m_diffuse_points.size() + 1);
            scattering_input_offsets.reserve(m_diffuse_points.size() + 1);
            const std::size_t coefficient_blocks =
                static_cast<std::size_t>(num_interior_points);
            m_rust_boundary_scattering_offsets.clear();
            m_rust_boundary_scattering_offsets.reserve(m_diffuse_points.size() -
                                                       coefficient_blocks + 1);
            m_rust_boundary_scattering_offsets.push_back(0);

            for (std::size_t point_index = 0;
                 point_index < m_diffuse_points.size(); ++point_index) {
                const auto& point = m_diffuse_points[point_index];
                const std::size_t output_size =
                    static_cast<std::size_t>(point->num_outgoing() * NSTOKES);
                const std::size_t input_size =
                    static_cast<std::size_t>(point->num_incoming() * NSTOKES);
                scattering_output_offsets.push_back(
                    scattering_output_offsets.back() + output_size);
                scattering_input_offsets.push_back(
                    scattering_input_offsets.back() + input_size);
                if (point_index >= coefficient_blocks) {
                    m_rust_boundary_scattering_offsets.push_back(
                        m_rust_boundary_scattering_offsets.back() +
                        static_cast<std::size_t>(point->num_outgoing()) *
                            point->num_incoming());
                }
            }

            std::vector<double> incoming_directions;
            std::vector<double> incoming_weights;
            std::vector<double> outgoing_directions;
            const auto& spheres = m_diffuse_points.front()->sphere_pair();
            incoming_directions.reserve(spheres.incoming_sphere().num_points() *
                                        3);
            incoming_weights.reserve(spheres.incoming_sphere().num_points());
            for (int index = 0; index < spheres.incoming_sphere().num_points();
                 ++index) {
                const auto direction =
                    spheres.incoming_sphere().get_quad_position(index);
                incoming_directions.insert(incoming_directions.end(),
                                           direction.data(),
                                           direction.data() + 3);
                incoming_weights.push_back(
                    spheres.incoming_sphere().quadrature_weight(index));
            }
            outgoing_directions.reserve(spheres.outgoing_sphere().num_points() *
                                        3);
            for (int index = 0; index < spheres.outgoing_sphere().num_points();
                 ++index) {
                const auto direction =
                    spheres.outgoing_sphere().get_quad_position(index);
                outgoing_directions.insert(outgoing_directions.end(),
                                           direction.data(),
                                           direction.data() + 3);
            }

            for (auto& storage : m_thread_storage) {
                if (rust_owned_transport) {
                    storage.rust_boundary_scattering_values.resize(0);
                } else {
                    storage.rust_boundary_scattering_values.resize(
                        m_rust_boundary_scattering_offsets.back());
                }
            }

            m_los_solution_cotangents.resize(
                m_config->num_wavelength_threads());
            for (auto& wavelength_cotangents : m_los_solution_cotangents) {
                wavelength_cotangents.resize(m_config->num_source_threads());
                for (auto& cotangent : wavelength_cotangents) {
                    cotangent.setZero(
                        static_cast<Eigen::Index>(m_num_outgoing_values));
                }
            }

            m_rust_solvers.clear();
            m_rust_solvers.reserve(m_thread_storage.size());
            for (std::size_t thread = 0; thread < m_thread_storage.size();
                 ++thread) {
                if constexpr (NSTOKES == 1) {
                    if (rust_owned_transport) {
                        m_rust_solvers.push_back(
                            sasktran2::rust::successive_orders::
                                new_scalar_coefficient_successive_orders_solver_with_sparsity(
                                    *m_rust_transport_sparsity.front(),
                                    as_rust_slice(scattering_output_offsets),
                                    as_rust_slice(scattering_input_offsets),
                                    coefficient_blocks,
                                    static_cast<std::size_t>(
                                        m_config->num_do_streams()),
                                    as_rust_slice(incoming_directions),
                                    as_rust_slice(incoming_weights),
                                    as_rust_slice(outgoing_directions),
                                    static_cast<std::size_t>(
                                        m_config
                                            ->successive_orders_max_iterations()),
                                    m_config
                                        ->successive_orders_relative_tolerance(),
                                    m_config
                                        ->successive_orders_absolute_tolerance(),
                                    static_cast<std::size_t>(
                                        m_config
                                            ->successive_orders_anderson_depth()),
                                    m_config->successive_orders_damping()));
                    } else {
                        m_rust_solvers.push_back(
                            sasktran2::rust::successive_orders::
                                new_scalar_coefficient_successive_orders_solver(
                                    static_cast<std::size_t>(
                                        m_thread_storage[thread]
                                            .m_incoming_radiances.value.size()),
                                    m_num_outgoing_values,
                                    as_rust_slice(m_outer_starts),
                                    as_rust_slice(m_inner_indicies),
                                    as_rust_slice(scattering_output_offsets),
                                    as_rust_slice(scattering_input_offsets),
                                    coefficient_blocks,
                                    static_cast<std::size_t>(
                                        m_config->num_do_streams()),
                                    as_rust_slice(incoming_directions),
                                    as_rust_slice(incoming_weights),
                                    as_rust_slice(outgoing_directions),
                                    static_cast<std::size_t>(
                                        m_config
                                            ->successive_orders_max_iterations()),
                                    m_config
                                        ->successive_orders_relative_tolerance(),
                                    m_config
                                        ->successive_orders_absolute_tolerance(),
                                    static_cast<std::size_t>(
                                        m_config
                                            ->successive_orders_anderson_depth()),
                                    m_config->successive_orders_damping()));
                    }
                } else {
                    if (rust_owned_transport) {
                        m_rust_solvers.push_back(
                            sasktran2::rust::successive_orders::
                                new_vector_coefficient_successive_orders_solver_with_sparsity(
                                    *m_rust_transport_sparsity.front(),
                                    as_rust_slice(scattering_output_offsets),
                                    as_rust_slice(scattering_input_offsets),
                                    coefficient_blocks,
                                    static_cast<std::size_t>(
                                        m_config->num_do_streams()),
                                    as_rust_slice(incoming_directions),
                                    as_rust_slice(incoming_weights),
                                    as_rust_slice(outgoing_directions),
                                    static_cast<std::size_t>(
                                        m_config
                                            ->successive_orders_max_iterations()),
                                    m_config
                                        ->successive_orders_relative_tolerance(),
                                    m_config
                                        ->successive_orders_absolute_tolerance(),
                                    static_cast<std::size_t>(
                                        m_config
                                            ->successive_orders_anderson_depth()),
                                    m_config->successive_orders_damping()));
                    } else {
                        m_rust_solvers.push_back(
                            sasktran2::rust::successive_orders::
                                new_vector_coefficient_successive_orders_solver(
                                    static_cast<std::size_t>(
                                        m_thread_storage[thread]
                                            .m_incoming_radiances.value.size()),
                                    m_num_outgoing_values,
                                    as_rust_slice(m_outer_starts),
                                    as_rust_slice(m_inner_indicies),
                                    as_rust_slice(scattering_output_offsets),
                                    as_rust_slice(scattering_input_offsets),
                                    coefficient_blocks,
                                    static_cast<std::size_t>(
                                        m_config->num_do_streams()),
                                    as_rust_slice(incoming_directions),
                                    as_rust_slice(incoming_weights),
                                    as_rust_slice(outgoing_directions),
                                    static_cast<std::size_t>(
                                        m_config
                                            ->successive_orders_max_iterations()),
                                    m_config
                                        ->successive_orders_relative_tolerance(),
                                    m_config
                                        ->successive_orders_absolute_tolerance(),
                                    static_cast<std::size_t>(
                                        m_config
                                            ->successive_orders_anderson_depth()),
                                    m_config->successive_orders_damping()));
                    }
                }
            }
            if (rust_owned_transport) {
                m_outer_starts.resize(0);
                m_inner_indicies.resize(0);
                for (auto& storage : m_thread_storage) {
                    storage.m_incoming_radiances.resize(0, 0, false);
                    storage.m_firstorder_radiances.resize(0, 0, false);
                    storage.rust_boundary_scattering_values.resize(0);
                    storage.rust_boundary_scattering_value_tangent.clear();
                    storage.rust_first_order_forcing_tangent.resize(0);
                }
            }
            if (m_config->successive_orders_initialization() !=
                sasktran2::Config::SuccessiveOrdersInitialization::twostream) {
                for (auto& storage : m_thread_storage) {
                    storage.m_outgoing_sources.resize(0, 0, false);
                }
                spdlog::debug(
                    "Released {:.3f} MiB of redundant C++ "
                    "successive-orders solution storage across {} wavelength "
                    "workers",
                    static_cast<double>(m_num_outgoing_values *
                                        m_thread_storage.size() *
                                        sizeof(double)) /
                        (1024.0 * 1024.0),
                    m_thread_storage.size());
            }
        }
#endif
    }

    template <int NSTOKES>
    void
    DiffuseTable<NSTOKES>::initialize_config(const sasktran2::Config& config) {
        m_config = &config;
        m_integrator.initialize_thread_storage(m_config->num_threads(), 1);

        if (m_use_rust_solver && m_config->num_hr_full_incoming_points() > 0) {
            throw std::invalid_argument(
                "The Rust successive-orders source currently requires all "
                "diffuse points to have explicitly traced incoming rays");
        }
        if (m_use_rust_solver && m_config->initialize_hr_with_do()) {
            throw std::invalid_argument(
                "The Rust successive-orders source does not yet support a "
                "discrete-ordinates initial guess");
        }

        m_thread_storage.resize(m_config->num_wavelength_threads());
        m_deferred_jvp_transport_restore.assign(
            m_config->num_wavelength_threads(), false);
        m_last_nonconverged_vjp_warning_wavelength.assign(
            m_config->num_wavelength_threads(), -1);
        m_active_vjp_derivatives.resize(m_config->num_wavelength_threads());

        if (m_geometry_1d != nullptr) {
            auto initial_single_scatter = std::make_unique<
                sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionTable,
                    NSTOKES>>(*m_geometry_1d, m_raytracer);
            if (m_use_rust_solver) {
                initial_single_scatter->enable_table_native_products();
            }
            m_initial_owned_sources.emplace_back(
                std::move(initial_single_scatter));
#ifdef SKTRAN_RUST_SUPPORT
        } else {
            m_initial_owned_sources.emplace_back(
                std::make_unique<
                    sasktran2::solartransmission::SingleScatterSource<
                        sasktran2::solartransmission::SolarTransmissionExact,
                        NSTOKES>>(*m_geometry_2d, *m_raytracer_2d, true));
#else
        } else {
            throw std::invalid_argument(
                "Geometry2D successive orders requires Rust support");
#endif
        }

        m_initial_sources.push_back(m_initial_owned_sources[0].get());

        if (m_use_rust_solver && m_altitude_grid.interpolation_method() !=
                                     sasktran2::grids::interpolation::lower) {
            if (m_geometry_1d != nullptr) {
                auto* source = dynamic_cast<
                    sasktran2::solartransmission::SingleScatterSource<
                        sasktran2::solartransmission::SolarTransmissionTable,
                        NSTOKES>*>(m_initial_sources.front());
                if (source == nullptr) {
                    throw std::logic_error(
                        "Rust source requires table solar transmission in "
                        "Geometry1D");
                }
                source->delegate_interior_source();
            } else {
                auto* source = dynamic_cast<
                    sasktran2::solartransmission::SingleScatterSource<
                        sasktran2::solartransmission::SolarTransmissionExact,
                        NSTOKES>*>(m_initial_sources.front());
                if (source == nullptr) {
                    throw std::logic_error(
                        "Rust source requires exact solar transmission in "
                        "Geometry2D");
                }
                source->delegate_interior_source();
            }
        }

        if (m_config->initialize_hr_with_do()) {
            m_initial_owned_sources.emplace_back(
                std::make_unique<
                    sasktran2::DOSourceInterpolatedPostProcessing<NSTOKES, -1>>(
                    *m_geometry_1d, m_raytracer, false));

            m_do_source =
                static_cast<DOSourceInterpolatedPostProcessing<NSTOKES, -1>*>(
                    m_initial_owned_sources[1].get());
        }

        for (auto& source : m_initial_owned_sources) {
            source->initialize_config(config);
        }
    }

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::release_cpp_transport_geometry() {
        if (m_rust_ray_transports.empty()) {
            return;
        }
        if (m_altitude_grid.interpolation_method() ==
            sasktran2::grids::interpolation::lower) {
            return;
        }
        SInterpolator().swap(m_diffuse_source_weights);

        std::vector<std::uint32_t>().swap(m_rust_transport_ray_layer_offsets);
        m_integrator.release_interior_geometry();

        // The single-scatter source has already packed its compact ground
        // geometry. Rust owns every atmospheric transport layer after packing.
        for (auto& ray : m_internal_viewing.traced_rays) {
            sasktran2::raytracing::TracedRay retained;
            retained.observer_and_look = ray.observer_and_look;
            retained.is_straight = ray.is_straight;
            retained.ground_is_hit = ray.ground_is_hit;
            retained.tangent_radius = ray.tangent_radius;
            ray = std::move(retained);
        }
    }
#endif

    template <int NSTOKES>
    Eigen::Vector3d DiffuseTable<NSTOKES>::rotate_unit_vector(
        const Eigen::Vector3d& vector, const Eigen::Vector3d& initial_position,
        const Eigen::Vector3d& new_position) const {
        // Reconstruct the interpolation location based on relative
        // azimuth/zenith angles
        sasktran2::Location temp;
        temp.position = initial_position;
        double csz_initial, saa_initial;
        sasktran2::raytracing::calculate_csz_saz(
            m_geometry.coordinates().sun_unit(), temp, vector, csz_initial,
            saa_initial);

        return m_geometry.coordinates().look_vector_from_azimuth(
            new_position, -saa_initial, temp.cos_zenith_angle(vector));
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::generate_source_interpolation_weights(
        const sasktran2::raytracing::TracedRay& ray,
        RaySourceInterpolationWeights<NSTOKES>& ray_interpolator,
        std::vector<std::pair<int, double>>& temp_location_storage,
        std::vector<std::pair<int, double>>& temp_direction_storage,
        std::vector<std::pair<int, double>>& temp_atmosphere_storage) const {
        int num_location, num_direction;
        Eigen::Vector3d rotated_los;
        sasktran2::Location temp_location;

        ray_interpolator = {};
        ray_interpolator.interior_weights.resize(ray.layers.size());
        ray_interpolator.atmosphere_weights.reserve(ray.num_grid_weights());
        ray_interpolator.source_indices.reserve(ray.layers.size() * 12);
        ray_interpolator.source_weights.reserve(ray.layers.size() * 12);

        for (int layeridx = 0; layeridx < ray.layers.size(); ++layeridx) {
            const auto& layer = ray.layers[layeridx];
            auto& layer_interpolator =
                ray_interpolator.interior_weights[layeridx];

            temp_location.position =
                (layer.entrance.position + layer.exit.position) / 2.0;

            m_geometry.assign_interpolation_weights(temp_location,
                                                    temp_atmosphere_storage);
            const auto atmosphere_nodes = ray.entrance_weights(layeridx);
            layer_interpolator.atmosphere_offset = static_cast<std::uint32_t>(
                ray_interpolator.atmosphere_weights.size());
            layer_interpolator.atmosphere_count =
                static_cast<std::uint32_t>(atmosphere_nodes.size());
            for (std::size_t node = 0; node < atmosphere_nodes.size(); ++node) {
                const int atmosphere_index = atmosphere_nodes[node].first;
                double atmosphere_weight = 0.0;
                for (const auto& index_weight : temp_atmosphere_storage) {
                    if (index_weight.first == atmosphere_index) {
                        atmosphere_weight = index_weight.second;
                        break;
                    }
                }
                ray_interpolator.atmosphere_weights.push_back(
                    atmosphere_weight);
            }

            layer_interpolator.source_offset = static_cast<std::uint32_t>(
                ray_interpolator.source_indices.size());

            m_location_interpolator->interior_interpolation_weights(
                m_geometry.coordinates(), temp_location, temp_location_storage,
                num_location);
            for (int locidx = 0; locidx < num_location; ++locidx) {
                const auto& contributing_point =
                    m_diffuse_points[temp_location_storage[locidx].first];

                if (m_geometry.coordinates().geometry_type() ==
                    sasktran2::geometrytype::spherical) {
                    rotated_los = rotate_unit_vector(
                        layer.average_look_away, temp_location.position,
                        contributing_point->location().position);
                } else {
                    rotated_los = layer.average_look_away;
                }

                contributing_point->sphere_pair().outgoing_sphere().interpolate(
                    rotated_los, temp_direction_storage, num_direction);

                for (int diridx = 0; diridx < num_direction; ++diridx) {
                    const int direction_index =
                        temp_direction_storage[diridx].first;
                    const int location_index =
                        temp_location_storage[locidx].first;
                    if (direction_index < 0 ||
                        direction_index >= contributing_point->num_outgoing()) {
                        throw std::runtime_error(
                            "Unit-sphere interpolation returned an invalid "
                            "direction index");
                    }
#ifdef SASKTRAN_DEBUG_ASSERTS
                    if (m_diffuse_outgoing_index_map[location_index] +
                            direction_index >
                        m_diffuse_outgoing_index_map.back() +
                            m_config->num_hr_outgoing()) {
                        spdlog::error("BAD INDEX {} {}", location_index,
                                      direction_index);
                    }
#endif

                    ray_interpolator.source_indices.push_back(
                        m_diffuse_outgoing_index_map[location_index] +
                        direction_index);
                    ray_interpolator.source_weights.push_back(
                        temp_location_storage[locidx].second *
                        temp_direction_storage[diridx].second);
                }
            }
            layer_interpolator.source_count = static_cast<std::uint32_t>(
                ray_interpolator.source_indices.size() -
                layer_interpolator.source_offset);
        }

        ray_interpolator.ground_is_hit = ray.ground_is_hit;
        if (ray_interpolator.ground_is_hit && !ray.layers.empty()) {
            const auto& layer = ray.layers[0];

            temp_location.position = layer.exit.position;

            m_location_interpolator->ground_interpolation_weights(
                m_geometry.coordinates(), temp_location, temp_location_storage,
                num_location);

            for (int locidx = 0; locidx < num_location; ++locidx) {
                const auto& contributing_point =
                    m_diffuse_points[temp_location_storage[locidx].first];

                if (m_geometry.coordinates().geometry_type() ==
                    sasktran2::geometrytype::spherical) {
                    rotated_los = rotate_unit_vector(
                        -1 * layer.average_look_away, temp_location.position,
                        contributing_point->location().position);
                } else {
                    rotated_los = -1 * layer.average_look_away;
                }

                contributing_point->sphere_pair().outgoing_sphere().interpolate(
                    rotated_los, temp_direction_storage, num_direction);

                for (int diridx = 0; diridx < num_direction; ++diridx) {
                    const int direction_index =
                        temp_direction_storage[diridx].first;
                    if (direction_index < 0 ||
                        direction_index >= contributing_point->num_outgoing()) {
                        throw std::runtime_error(
                            "Ground unit-sphere interpolation returned an "
                            "invalid direction index");
                    }
                    ray_interpolator.ground_source_indices.push_back(
                        m_diffuse_outgoing_index_map
                            [temp_location_storage[locidx].first] +
                        direction_index);
                    ray_interpolator.ground_source_weights.push_back(
                        temp_location_storage[locidx].second *
                        temp_direction_storage[diridx].second);
                }
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::generate_source_interpolation_weights(
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
        SInterpolator& interpolator) const {
        ZoneScopedN("Generate Source Interpolation Weights");
        interpolator.resize(rays.size());

        const int nthreads = m_config->num_threads();
        std::vector<std::vector<std::pair<int, double>>>
            thread_temp_location_storage(nthreads);
        std::vector<std::vector<std::pair<int, double>>>
            thread_temp_direction_storage(nthreads);
        std::vector<std::vector<std::pair<int, double>>>
            thread_temp_atmosphere_storage(nthreads);

#pragma omp parallel for num_threads(nthreads)
        for (int rayidx = 0; rayidx < rays.size(); ++rayidx) {
#ifdef SKTRAN_OPENMP_SUPPORT
            const int threadidx = omp_get_thread_num();
#else
            const int threadidx = 0;
#endif
            generate_source_interpolation_weights(
                rays[rayidx], interpolator[rayidx],
                thread_temp_location_storage[threadidx],
                thread_temp_direction_storage[threadidx],
                thread_temp_atmosphere_storage[threadidx]);
        }
    }

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::initialize_rust_source_interpolator(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        ZoneScopedN("Pack Rust LOS Source Interpolator");
        if (internal_viewing.traced_rays.size() !=
            m_los_source_weights.size()) {
            throw std::logic_error(
                "LOS source interpolation ray count mismatch");
        }
        if (m_thread_storage.empty()) {
            throw std::logic_error(
                "LOS source interpolation requires thread storage");
        }

        const auto checked_u32 = [](std::size_t value,
                                    const char* description) {
            if (value > std::numeric_limits<std::uint32_t>::max()) {
                throw std::length_error(std::string(description) +
                                        " exceeds the packed offset range");
            }
            return static_cast<std::uint32_t>(value);
        };
        const auto checked_index = [&](int value, const char* description) {
            if (value < 0) {
                throw std::logic_error(std::string("Negative ") + description +
                                       " in LOS source interpolation");
            }
            return checked_u32(static_cast<std::size_t>(value), description);
        };

        std::size_t num_layers = 0;
        std::size_t num_atmosphere_weights = 0;
        std::size_t num_source_weights = 0;
        std::size_t num_ground_weights = 0;
        for (const auto& ray : m_los_source_weights) {
            num_layers += ray.interior_weights.size();
            num_atmosphere_weights += ray.atmosphere_weights.size();
            num_source_weights += ray.source_indices.size();
            num_ground_weights += ray.ground_source_indices.size();
        }

        std::vector<std::uint32_t> ray_layer_offsets;
        std::vector<std::uint32_t> layer_atmosphere_offsets;
        std::vector<std::uint32_t> atmosphere_indices;
        std::vector<double> atmosphere_weights;
        std::vector<double> optical_depth_weights;
        std::vector<std::uint32_t> layer_source_offsets;
        std::vector<std::uint32_t> source_indices;
        std::vector<double> source_weights;
        std::vector<std::uint32_t> ray_ground_offsets;
        std::vector<std::uint32_t> ground_source_indices;
        std::vector<double> ground_source_weights;

        ray_layer_offsets.reserve(m_los_source_weights.size() + 1);
        layer_atmosphere_offsets.reserve(num_layers + 1);
        atmosphere_indices.reserve(num_atmosphere_weights);
        atmosphere_weights.reserve(num_atmosphere_weights);
        optical_depth_weights.reserve(num_atmosphere_weights);
        layer_source_offsets.reserve(num_layers + 1);
        source_indices.reserve(num_source_weights);
        source_weights.reserve(num_source_weights);
        ray_ground_offsets.reserve(m_los_source_weights.size() + 1);
        ground_source_indices.reserve(num_ground_weights);
        ground_source_weights.reserve(num_ground_weights);
        ray_layer_offsets.push_back(0);
        layer_atmosphere_offsets.push_back(0);
        layer_source_offsets.push_back(0);
        ray_ground_offsets.push_back(0);

        for (std::size_t ray_index = 0; ray_index < m_los_source_weights.size();
             ++ray_index) {
            const auto& ray_interpolator = m_los_source_weights[ray_index];
            const auto& traced_ray = internal_viewing.traced_rays[ray_index];
            if (ray_interpolator.interior_weights.size() !=
                traced_ray.layers.size()) {
                throw std::logic_error(
                    "LOS source interpolation layer count mismatch");
            }

            for (std::size_t layer_index = 0;
                 layer_index < ray_interpolator.interior_weights.size();
                 ++layer_index) {
                const auto& layer =
                    ray_interpolator.interior_weights[layer_index];
                const auto entrance_weights =
                    traced_ray.entrance_weights(layer_index);
                const auto layer_optical_depth_weights =
                    traced_ray.optical_depth_weights(layer_index);
                if (entrance_weights.size() != layer.atmosphere_count) {
                    throw std::logic_error(
                        "LOS source atmosphere interpolation size mismatch");
                }
                if (layer_optical_depth_weights.size() !=
                    entrance_weights.size()) {
                    throw std::logic_error(
                        "LOS source optical-depth interpolation size "
                        "mismatch");
                }

                for (std::size_t index = 0; index < layer.atmosphere_count;
                     ++index) {
                    if (layer_optical_depth_weights[index].first !=
                        entrance_weights[index].first) {
                        throw std::logic_error(
                            "LOS source optical-depth atmosphere index "
                            "mismatch");
                    }
                    atmosphere_indices.push_back(checked_index(
                        entrance_weights[index].first, "atmosphere index"));
                    atmosphere_weights.push_back(
                        ray_interpolator
                            .atmosphere_weights[layer.atmosphere_offset +
                                                index]);
                    optical_depth_weights.push_back(
                        layer_optical_depth_weights[index].second);
                }
                layer_atmosphere_offsets.push_back(checked_u32(
                    atmosphere_indices.size(), "atmosphere weight count"));

                for (std::size_t index = 0; index < layer.source_count;
                     ++index) {
                    const std::size_t packed_index =
                        layer.source_offset + index;
                    source_indices.push_back(checked_index(
                        ray_interpolator.source_indices[packed_index],
                        "source index"));
                    source_weights.push_back(
                        ray_interpolator.source_weights[packed_index]);
                }
                layer_source_offsets.push_back(
                    checked_u32(source_indices.size(), "source weight count"));
            }
            ray_layer_offsets.push_back(checked_u32(
                layer_source_offsets.size() - 1, "LOS source layer count"));

            for (std::size_t index = 0;
                 index < ray_interpolator.ground_source_indices.size();
                 ++index) {
                ground_source_indices.push_back(
                    checked_index(ray_interpolator.ground_source_indices[index],
                                  "ground source index"));
                ground_source_weights.push_back(
                    ray_interpolator.ground_source_weights[index]);
            }
            ray_ground_offsets.push_back(checked_u32(
                ground_source_indices.size(), "ground source weight count"));
        }

        const auto solution_size = m_num_outgoing_values;
        if (solution_size == 0 || solution_size % NSTOKES != 0) {
            throw std::logic_error(
                "LOS source interpolation solution size mismatch");
        }
        m_rust_source_interpolators.clear();
        m_rust_source_interpolators.push_back(
            sasktran2::rust::successive_orders::new_source_interpolator(
                as_rust_slice(ray_layer_offsets),
                as_rust_slice(layer_atmosphere_offsets),
                as_rust_slice(atmosphere_indices),
                as_rust_slice(atmosphere_weights),
                as_rust_slice(optical_depth_weights),
                as_rust_slice(layer_source_offsets),
                as_rust_slice(source_indices), as_rust_slice(source_weights),
                as_rust_slice(ray_ground_offsets),
                as_rust_slice(ground_source_indices),
                as_rust_slice(ground_source_weights), m_geometry.size(),
                solution_size / NSTOKES, NSTOKES));

        spdlog::debug(
            "Packed Rust successive-orders LOS source interpolation: {} "
            "rays, {} layers, {:.3f} MiB immutable",
            sasktran2::rust::successive_orders::source_interpolator_num_rays(
                *m_rust_source_interpolators.front()),
            sasktran2::rust::successive_orders::source_interpolator_num_layers(
                *m_rust_source_interpolators.front()),
            static_cast<double>(sasktran2::rust::successive_orders::
                                    source_interpolator_storage_bytes(
                                        *m_rust_source_interpolators.front())) /
                (1024.0 * 1024.0));

        // The Rust source integrates the complete LOS through its whole-ray
        // callbacks. Rust therefore owns the only LOS interpolation topology
        // required by this source after packing.
        SInterpolator().swap(m_los_source_weights);
    }
#endif

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::compile_accumulation_row(
        RaySourceInterpolationWeights<NSTOKES>& weights,
        std::vector<int>& transport_columns,
        std::vector<std::pair<int, std::uint16_t*>>& sorting_helper) const {
        weights.source_accumulation_inner_indices.resize(
            weights.source_indices.size());
        weights.ground_accumulation_inner_indices.resize(
            weights.ground_source_indices.size());

        sorting_helper.clear();
        sorting_helper.reserve(weights.source_indices.size() +
                               weights.ground_source_indices.size());
        for (std::size_t index = 0; index < weights.source_indices.size();
             ++index) {
            sorting_helper.emplace_back(
                weights.source_indices[index],
                &weights.source_accumulation_inner_indices[index]);
        }
        for (std::size_t index = 0;
             index < weights.ground_source_indices.size(); ++index) {
            sorting_helper.emplace_back(
                weights.ground_source_indices[index],
                &weights.ground_accumulation_inner_indices[index]);
        }
        std::stable_sort(sorting_helper.begin(), sorting_helper.end(),
                         [](const auto& left, const auto& right) {
                             return left.first < right.first;
                         });

        transport_columns.clear();
        transport_columns.reserve(sorting_helper.size());
        int previous_column = -1;
        for (auto& [column, inner_index] : sorting_helper) {
            if (transport_columns.empty() || column != previous_column) {
                if (transport_columns.size() >
                    std::numeric_limits<std::uint16_t>::max()) {
                    throw std::length_error(
                        "Successive-orders transport row exceeds compact "
                        "inner index range");
                }
                transport_columns.push_back(column);
                previous_column = column;
            }
            *inner_index =
                static_cast<std::uint16_t>(transport_columns.size() - 1);
        }

        std::vector<int>().swap(weights.source_indices);
        std::vector<int>().swap(weights.ground_source_indices);
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::generate_scattering_matrices(int wavelidx,
                                                             int threadidx) {
        ZoneScopedN("Generate Scattering Matrices");
        auto& storage = m_thread_storage[threadidx];

#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            if (m_rust_scattering_interpolators.empty()) {
                throw std::logic_error(
                    "Rust scattering interpolator is not initialized");
            }
            const auto& coefficients = m_atmosphere->storage().leg_coeff;
            const std::size_t coefficient_stride =
                static_cast<std::size_t>(coefficients.dimension(0));
            const std::size_t num_geometry =
                static_cast<std::size_t>(coefficients.dimension(1));
            if (wavelidx < 0 || wavelidx >= coefficients.dimension(2)) {
                throw std::out_of_range(
                    "Successive-orders scattering wavelength is out of range");
            }
            const std::size_t wavelength_size =
                coefficient_stride * num_geometry;
            const double* wavelength_coefficients =
                coefficients.data() +
                static_cast<std::size_t>(wavelidx) * wavelength_size;
            sasktran2::rust::successive_orders::
                set_scattering_coefficients_from_atmosphere(
                    *m_rust_solvers[threadidx],
                    *m_rust_scattering_interpolators.front(),
                    {wavelength_coefficients, wavelength_size},
                    coefficient_stride);
            const bool rust_owned_transport =
                m_altitude_grid.interpolation_method() !=
                sasktran2::grids::interpolation::lower;
            auto boundary_values =
                rust_owned_transport
                    ? sasktran2::rust::successive_orders::
                          boundary_scattering_values_mut(
                              *m_rust_solvers[threadidx])
                    : as_rust_mut_slice(
                          storage.rust_boundary_scattering_values);
            if (boundary_values.size() !=
                m_rust_boundary_scattering_offsets.back()) {
                throw std::logic_error(
                    "Rust boundary scattering storage size mismatch");
            }
            std::fill(boundary_values.begin(), boundary_values.end(), 0.0);

#pragma omp parallel for num_threads(m_config->num_source_threads())           \
    schedule(dynamic)
            for (int ground_index = 0;
                 ground_index < m_location_interpolator->num_ground_points();
                 ++ground_index) {
                const std::size_t point_index = static_cast<std::size_t>(
                    m_location_interpolator->num_interior_points() +
                    ground_index);
                if (!m_diffuse_point_full_calculation[point_index]) {
                    continue;
                }
                const auto& point = m_diffuse_points[point_index];
                point->sphere_pair().calculate_ground_scattering_values(
                    m_atmosphere->surface(), point->location(), wavelidx,
                    boundary_values.data() +
                        m_rust_boundary_scattering_offsets[ground_index]);
            }
        } else
#endif
        {
#pragma omp parallel for num_threads(m_config->num_source_threads())           \
    schedule(dynamic)
            for (int i = 0; i < m_location_interpolator->num_interior_points();
                 ++i) {
                if (!m_diffuse_point_full_calculation[i]) {
                    continue;
                }

                const auto& point = m_diffuse_points[i];
                point->sphere_pair().calculate_scattering_matrix(
                    m_atmosphere->storage(), wavelidx,
                    m_diffuse_point_interpolation_weights[i],
                    storage.point_scattering_matrices[i].data());
            }

#pragma omp parallel for num_threads(m_config->num_source_threads())           \
    schedule(dynamic)
            for (int i = 0; i < m_location_interpolator->num_ground_points();
                 ++i) {
                if (!m_diffuse_point_full_calculation
                        [i + m_location_interpolator->num_interior_points()]) {
                    continue;
                }

                const auto& point =
                    m_diffuse_points[i + m_location_interpolator
                                             ->num_interior_points()];

                point->sphere_pair().calculate_ground_scattering_matrix(
                    m_atmosphere->surface(),
                    m_diffuse_point_interpolation_weights
                        [i + m_location_interpolator->num_interior_points()],
                    point->location(), wavelidx,
                    storage
                        .point_scattering_matrices
                            [i + m_location_interpolator->num_interior_points()]
                        .data());
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::generate_accumulation_matrix(int wavelidx,
                                                             int threadidx) {}

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::iterate_to_solution(int wavelidx,
                                                    int threadidx) {
        ZoneScopedN("Iterate to Solution");
#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            iterate_to_solution_rust(wavelidx, threadidx);
            return;
        }
#endif
        auto& storage = m_thread_storage[threadidx];

        Eigen::Map<Eigen::SparseMatrix<double, Eigen::RowMajor>>
            accumulation_matrix(
                m_thread_storage[threadidx].m_incoming_radiances.value.size(),
                m_thread_storage[threadidx].m_outgoing_sources.value.size(),
                m_inner_indicies.size(), m_outer_starts.data(),
                m_inner_indicies.data(),
                m_thread_storage[threadidx].accumulation_summed_values.data());

        Eigen::VectorXd old_outgoing_vals;

        if (m_config->initialize_hr_with_do()) {
            m_thread_storage[threadidx].m_outgoing_sources.value =
                m_do_to_diffuse_outgoing_interpolator *
                m_do_source->storage().linear_source(threadidx).value;

            if (m_config->wf_precision() ==
                sasktran2::Config::WeightingFunctionPrecision::full) {
                m_thread_storage[threadidx].m_outgoing_sources.deriv =
                    m_do_to_diffuse_outgoing_interpolator *
                    m_do_source->storage().linear_source(threadidx).deriv;
            }
        }

        for (int spher_iter = 0;
             spher_iter < m_config->num_hr_spherical_iterations();
             ++spher_iter) {
            // Apply the scattering matrices
            if (spher_iter == 0) {
                if (m_config->initialize_hr_with_do()) {
                    // Apply the accumulation matrix
                    storage.m_incoming_radiances.value =
                        accumulation_matrix * storage.m_outgoing_sources.value +
                        storage.m_firstorder_radiances.value;
                } else {
                    storage.m_incoming_radiances.value =
                        storage.m_firstorder_radiances.value;
                }
            } else {
                ZoneScopedN("Apply Accumulation Matrix");
                storage.m_incoming_radiances.value.noalias() =
                    accumulation_matrix * storage.m_outgoing_sources.value +
                    storage.m_firstorder_radiances.value;
            }

            old_outgoing_vals =
                m_thread_storage[threadidx].m_outgoing_sources.value;

#pragma omp parallel for num_threads(m_config->num_source_threads())           \
    schedule(dynamic)
            for (int i = 0; i < m_diffuse_points.size(); ++i) {
                if (!m_diffuse_point_full_calculation[i]) {
                    continue;
                }

                const auto& point = m_diffuse_points[i];

                auto incoming_seq =
                    Eigen::seq(m_diffuse_incoming_index_map[i] * NSTOKES,
                               NSTOKES * (m_diffuse_incoming_index_map[i] +
                                          point->num_incoming()) -
                                   1);
                auto outgoing_seq =
                    Eigen::seq(m_diffuse_outgoing_index_map[i] * NSTOKES,
                               NSTOKES * (m_diffuse_outgoing_index_map[i] +
                                          point->num_outgoing()) -
                                   1);

                {
                    ZoneScopedN("Apply Scattering Matrix");

                    storage.m_outgoing_sources.value(outgoing_seq).noalias() =
                        storage.point_scattering_matrices[i] *
                        storage.m_incoming_radiances.value(incoming_seq);
                }

#ifdef SASKTRAN_DEBUG_ASSERTS
                if (storage.m_outgoing_sources.value(outgoing_seq).hasNaN()) {
                    spdlog::error("NaN in outgoing point: {}", i);
                }
#endif
            }

            interpolate_sources(old_outgoing_vals, storage.m_outgoing_sources);

            if (m_config->wf_precision() ==
                    sasktran2::Config::WeightingFunctionPrecision::full &&
                m_config->initialize_hr_with_do()) {
                storage.m_outgoing_sources.deriv.array().colwise() *=
                    storage.m_outgoing_sources.value.array() /
                    old_outgoing_vals.array();
            }
        }
    }

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    void
    DiffuseTable<NSTOKES>::generate_twostream_initial_guess(int wavelidx,
                                                            int threadidx) {
        if (m_config->successive_orders_initialization() !=
            sasktran2::Config::SuccessiveOrdersInitialization::twostream) {
            return;
        }
        if (m_atmosphere == nullptr || threadidx < 0 ||
            static_cast<std::size_t>(threadidx) >=
                m_twostream_initializers.size()) {
            throw std::logic_error(
                "Two-stream successive-orders initializer is not ready");
        }

        auto& thread_storage = m_thread_storage[threadidx];
        thread_storage.m_outgoing_sources.value.setZero();
        const auto& atmosphere_storage = m_atmosphere->storage();
        const std::size_t num_levels =
            static_cast<std::size_t>(m_altitude_grid.grid().size());
        const std::size_t coefficient_families = NSTOKES == 1 ? 1 : 4;
        const std::size_t num_legendre =
            static_cast<std::size_t>(
                atmosphere_storage.leg_coeff.dimension(0)) /
            coefficient_families;
        if (num_legendre < 2) {
            throw std::invalid_argument(
                "Two-stream initialization requires at least two Legendre "
                "moments");
        }

        auto& extinction = thread_storage.twostream_extinction;
        auto& ssa = thread_storage.twostream_ssa;
        auto& legendre = thread_storage.twostream_legendre;
        auto& delta_m = thread_storage.twostream_delta_m;
        auto& emission = thread_storage.twostream_emission;
        auto& surface_albedo = thread_storage.twostream_surface_albedo;
        extinction.resize(num_levels);
        ssa.resize(num_levels);
        legendre.resize(num_legendre * num_levels);
        delta_m.resize(num_levels);
        emission.resize(num_levels);
        surface_albedo.resize(1);
        surface_albedo[0] =
            m_atmosphere->surface().brdf(wavelidx, 0, 0, 0)(0, 0) * EIGEN_PI;

        const int num_profiles = m_location_interpolator->num_ground_points();
        std::vector<double> scattering_legendre(num_legendre);
        for (int profile = 0; profile < num_profiles; ++profile) {
            for (std::size_t level = 0; level < num_levels; ++level) {
                double local_extinction = 0;
                double scattering_extinction = 0;
                double local_emission = 0;
                double scattering_delta_m = 0;
                std::fill(scattering_legendre.begin(),
                          scattering_legendre.end(), 0.0);
                for (const auto& [location_index, weight] :
                     m_twostream_atmosphere_weights[profile][level]) {
                    const double node_extinction =
                        atmosphere_storage.total_extinction(location_index,
                                                            wavelidx);
                    const double node_scattering =
                        node_extinction *
                        atmosphere_storage.ssa(location_index, wavelidx);
                    local_extinction += weight * node_extinction;
                    scattering_extinction += weight * node_scattering;
                    scattering_delta_m +=
                        weight * node_scattering *
                        atmosphere_storage.f(location_index, wavelidx);
                    local_emission +=
                        weight * atmosphere_storage.emission_source(
                                     location_index, wavelidx);
                    for (std::size_t moment = 0; moment < num_legendre;
                         ++moment) {
                        scattering_legendre[moment] +=
                            weight * node_scattering *
                            atmosphere_storage.leg_coeff(
                                moment * coefficient_families, location_index,
                                wavelidx);
                    }
                }
                extinction[level] = local_extinction;
                ssa[level] = local_extinction > 0
                                 ? scattering_extinction / local_extinction
                                 : 0;
                delta_m[level] =
                    scattering_extinction > 0
                        ? scattering_delta_m / scattering_extinction
                        : 0;
                emission[level] = local_emission;
                for (std::size_t moment = 0; moment < num_legendre; ++moment) {
                    legendre[level * num_legendre + moment] =
                        scattering_extinction > 0
                            ? scattering_legendre[moment] /
                                  scattering_extinction
                            : 0;
                }
            }

            const std::array<double, 1> surface_emission = {
                m_atmosphere->surface().emission()(wavelidx)};
            const std::array<double, 1> solar_irradiance = {
                atmosphere_storage.solar_irradiance(wavelidx)};
            auto& source = m_twostream_initializers[threadidx][profile];
            sasktran2::rust::twostream::solve(
                *source, 1, num_levels, as_rust_slice(extinction),
                as_rust_slice(ssa), as_rust_slice(legendre), num_legendre,
                as_rust_slice(delta_m), as_rust_slice(emission),
                as_rust_slice(surface_albedo),
                ::rust::Slice<const double>(surface_emission.data(), 1),
                ::rust::Slice<const double>(solar_irradiance.data(), 1), false);

            {
                const auto sampled_source =
                    sasktran2::rust::twostream::radiance(*source);
                const int first_point =
                    profile * m_twostream_num_source_altitudes;
                const int num_outgoing =
                    m_diffuse_points[first_point]->num_outgoing();
                const std::size_t num_volume_samples =
                    static_cast<std::size_t>(m_twostream_num_source_altitudes) *
                    num_outgoing;
                if (sampled_source.size() != num_volume_samples + 1) {
                    throw std::runtime_error(
                        "Rust two-stream initializer returned an invalid "
                        "source size");
                }
                std::size_t sample = 0;
                for (int altitude_index = 0;
                     altitude_index < m_twostream_num_source_altitudes;
                     ++altitude_index) {
                    const int point_index = first_point + altitude_index;
                    for (int direction_index = 0;
                         direction_index < num_outgoing; ++direction_index) {
                        const int state_index =
                            NSTOKES *
                            (m_diffuse_outgoing_index_map[point_index] +
                             direction_index);
                        thread_storage.m_outgoing_sources.value(state_index) =
                            sampled_source[sample++];
                    }
                }

                const int ground_point =
                    m_location_interpolator->num_interior_points() + profile;
                const double ground_source = sampled_source[num_volume_samples];
                for (int direction_index = 0;
                     direction_index <
                     m_diffuse_points[ground_point]->num_outgoing();
                     ++direction_index) {
                    const int state_index =
                        NSTOKES * (m_diffuse_outgoing_index_map[ground_point] +
                                   direction_index);
                    thread_storage.m_outgoing_sources.value(state_index) =
                        ground_source;
                }
            }
            sasktran2::rust::twostream::clear_output(*source);
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::pack_rust_scattering_coefficient_jvp(
        int wavelidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        if (m_rust_scattering_interpolators.empty()) {
            throw std::logic_error(
                "Rust scattering interpolator is not initialized");
        }
        const auto& derivatives = m_atmosphere->storage().d_leg_coeff;
        const std::size_t coefficient_stride =
            static_cast<std::size_t>(derivatives.dimension(0));
        const std::size_t num_geometry =
            static_cast<std::size_t>(derivatives.dimension(1));
        const std::size_t num_wavelengths =
            static_cast<std::size_t>(derivatives.dimension(2));
        const std::size_t num_groups =
            static_cast<std::size_t>(derivatives.dimension(3));
        if (wavelidx < 0 ||
            static_cast<std::size_t>(wavelidx) >= num_wavelengths) {
            throw std::out_of_range(
                "Successive-orders scattering JVP wavelength is out of range");
        }
        const std::size_t tangent_size = num_geometry * num_groups;
        const double* scattering_tangent =
            native_tangent.data() + m_atmosphere->scat_deriv_start_index();
        sasktran2::rust::successive_orders::
            set_scattering_coefficient_jvp_from_atmosphere(
                *m_rust_solvers[threadidx],
                *m_rust_scattering_interpolators.front(),
                {derivatives.data(),
                 static_cast<std::size_t>(derivatives.size())},
                coefficient_stride, num_wavelengths,
                static_cast<std::size_t>(wavelidx),
                {scattering_tangent, tangent_size});
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::pack_rust_boundary_scattering_jvp(
        int wavelidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        const bool rust_owned_transport =
            m_altitude_grid.interpolation_method() !=
            sasktran2::grids::interpolation::lower;
        auto tangent = rust_owned_transport
                           ? sasktran2::rust::successive_orders::
                                 boundary_scattering_value_tangent_mut(
                                     *m_rust_solvers[threadidx])
                           : as_rust_mut_slice(
                                 m_thread_storage[threadidx]
                                     .rust_boundary_scattering_value_tangent);
        if (tangent.size() != m_rust_boundary_scattering_offsets.back()) {
            throw std::logic_error(
                "Rust boundary scattering tangent size mismatch");
        }
        std::fill(tangent.begin(), tangent.end(), 0.0);
        const int num_surface_derivatives = m_atmosphere->surface().num_deriv();
        if (num_surface_derivatives == 0 ||
            native_tangent
                .segment(m_atmosphere->surface_deriv_start_index(),
                         num_surface_derivatives)
                .isZero(0.0)) {
            return;
        }
        std::size_t value_offset = 0;
        for (std::size_t point_index = static_cast<std::size_t>(
                 m_location_interpolator->num_interior_points());
             point_index < m_diffuse_points.size(); ++point_index) {
            const auto& point = *m_diffuse_points[point_index];
            const auto& spheres = point.sphere_pair();
            const auto vertical = point.location().position.normalized();
            for (int outgoing_index = 0; outgoing_index < point.num_outgoing();
                 ++outgoing_index) {
                const auto outgoing =
                    spheres.outgoing_sphere().get_quad_position(outgoing_index);
                const double mu_out =
                    point.location().cos_zenith_angle(outgoing);
                const auto horiz_out =
                    (outgoing - mu_out * vertical).normalized();
                for (int incoming_index = 0;
                     incoming_index < point.num_incoming(); ++incoming_index) {
                    const auto incoming =
                        spheres.incoming_sphere().get_quad_position(
                            incoming_index);
                    const double mu_in =
                        point.location().cos_zenith_angle(incoming);
                    const auto horiz_in =
                        (incoming - mu_in * vertical).normalized();
                    const double phi_diff =
                        EIGEN_PI - std::acos(horiz_in.dot(horiz_out));
                    double brdf_jvp = 0.0;
                    for (int derivative = 0;
                         derivative < m_atmosphere->surface().num_deriv();
                         ++derivative) {
                        brdf_jvp +=
                            native_tangent(
                                m_atmosphere->surface_deriv_start_index() +
                                derivative) *
                            m_atmosphere->surface().d_brdf(wavelidx, mu_in,
                                                           mu_out, phi_diff,
                                                           derivative)(0, 0);
                    }
                    tangent[value_offset++] =
                        4.0 * EIGEN_PI * brdf_jvp * mu_in *
                        spheres.incoming_sphere().quadrature_weight(
                            incoming_index);
                }
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::accumulate_rust_scattering_coefficient_vjp(
        int wavelidx, int threadidx,
        const std::vector<int>& active_scattering_groups,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        if (m_rust_scattering_interpolators.empty()) {
            throw std::logic_error(
                "Rust scattering interpolator is not initialized");
        }
        static_assert(sizeof(int) == sizeof(std::int32_t));
        const auto& derivatives = m_atmosphere->storage().d_leg_coeff;
        const std::size_t coefficient_stride =
            static_cast<std::size_t>(derivatives.dimension(0));
        const std::size_t num_geometry =
            static_cast<std::size_t>(derivatives.dimension(1));
        const std::size_t num_wavelengths =
            static_cast<std::size_t>(derivatives.dimension(2));
        const std::size_t num_groups =
            static_cast<std::size_t>(derivatives.dimension(3));
        if (wavelidx < 0 ||
            static_cast<std::size_t>(wavelidx) >= num_wavelengths) {
            throw std::out_of_range(
                "Successive-orders scattering VJP wavelength is out of range");
        }
        double* scattering_gradient =
            native_gradient.data() + m_atmosphere->scat_deriv_start_index();
        sasktran2::rust::successive_orders::
            accumulate_solver_scattering_coefficient_vjp(
                *m_rust_solvers[threadidx],
                *m_rust_scattering_interpolators.front(),
                {derivatives.data(),
                 static_cast<std::size_t>(derivatives.size())},
                coefficient_stride, num_wavelengths,
                static_cast<std::size_t>(wavelidx),
                {reinterpret_cast<const std::int32_t*>(
                     active_scattering_groups.data()),
                 active_scattering_groups.size()},
                {scattering_gradient, num_geometry * num_groups});
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::accumulate_rust_boundary_scattering_vjp(
        int wavelidx, ::rust::Slice<const double> boundary_gradient,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        std::size_t value_offset = 0;
        for (std::size_t point_index = static_cast<std::size_t>(
                 m_location_interpolator->num_interior_points());
             point_index < m_diffuse_points.size(); ++point_index) {
            const auto& point = *m_diffuse_points[point_index];
            const auto& spheres = point.sphere_pair();
            const auto vertical = point.location().position.normalized();
            for (int outgoing_index = 0; outgoing_index < point.num_outgoing();
                 ++outgoing_index) {
                const auto outgoing =
                    spheres.outgoing_sphere().get_quad_position(outgoing_index);
                const double mu_out =
                    point.location().cos_zenith_angle(outgoing);
                const auto horiz_out =
                    (outgoing - mu_out * vertical).normalized();
                for (int incoming_index = 0;
                     incoming_index < point.num_incoming(); ++incoming_index) {
                    const auto incoming =
                        spheres.incoming_sphere().get_quad_position(
                            incoming_index);
                    const double mu_in =
                        point.location().cos_zenith_angle(incoming);
                    const auto horiz_in =
                        (incoming - mu_in * vertical).normalized();
                    const double phi_diff =
                        EIGEN_PI - std::acos(horiz_in.dot(horiz_out));
                    const double scale =
                        boundary_gradient[value_offset++] * 4.0 * EIGEN_PI *
                        mu_in *
                        spheres.incoming_sphere().quadrature_weight(
                            incoming_index);
                    for (int derivative = 0;
                         derivative < m_atmosphere->surface().num_deriv();
                         ++derivative) {
                        native_gradient(
                            m_atmosphere->surface_deriv_start_index() +
                            derivative) +=
                            scale * m_atmosphere->surface().d_brdf(
                                        wavelidx, mu_in, mu_out, phi_diff,
                                        derivative)(0, 0);
                    }
                }
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::iterate_to_solution_rust(int wavelidx,
                                                         int threadidx) {
        ZoneScopedN("Rust Successive Orders Solve");
        auto& storage = m_thread_storage[threadidx];

        const std::vector<double> no_initial_guess;
        generate_twostream_initial_guess(wavelidx, threadidx);
        const bool use_twostream_initial_guess =
            m_config->successive_orders_initialization() ==
            sasktran2::Config::SuccessiveOrdersInitialization::twostream;
        const auto initial_guess =
            use_twostream_initial_guess
                ? as_rust_slice(storage.m_outgoing_sources.value)
                : as_rust_slice(no_initial_guess);
        auto& solver = m_rust_solvers[threadidx];
        const bool rust_owned_transport =
            m_altitude_grid.interpolation_method() !=
            sasktran2::grids::interpolation::lower;
        if (rust_owned_transport) {
            sasktran2::rust::successive_orders::solve_coefficients_assembled(
                *solver, initial_guess);
        } else if constexpr (NSTOKES == 1) {
            sasktran2::rust::successive_orders::solve_scalar_coefficients(
                *solver, as_rust_slice(m_outer_starts),
                as_rust_slice(m_inner_indicies),
                as_rust_slice(storage.accumulation_summed_values),
                as_rust_slice(storage.rust_boundary_scattering_values),
                as_rust_slice(storage.m_firstorder_radiances.value),
                initial_guess);
        } else {
            sasktran2::rust::successive_orders::solve_vector_coefficients(
                *solver, as_rust_slice(m_outer_starts),
                as_rust_slice(m_inner_indicies),
                as_rust_slice(storage.accumulation_summed_values),
                as_rust_slice(storage.rust_boundary_scattering_values),
                as_rust_slice(storage.m_firstorder_radiances.value),
                initial_guess);
        }
        const auto solution =
            sasktran2::rust::successive_orders::solution(*solver);
        if (solution.size() != m_num_outgoing_values) {
            throw std::runtime_error(
                "Rust successive-orders solver returned an invalid solution "
                "size");
        }

        const bool tolerance_requested =
            m_config->successive_orders_relative_tolerance() > 0 ||
            m_config->successive_orders_absolute_tolerance() > 0;
        if (tolerance_requested &&
            !sasktran2::rust::successive_orders::converged(*solver)) {
            spdlog::warn(
                "Rust successive-orders source reached {} iterations with "
                "residual {}",
                sasktran2::rust::successive_orders::iterations(*solver),
                sasktran2::rust::successive_orders::final_residual(*solver));
        } else {
            spdlog::debug(
                "Rust successive-orders source completed {} iterations; "
                "initial residual {}, final residual {}",
                sasktran2::rust::successive_orders::iterations(*solver),
                sasktran2::rust::successive_orders::initial_residual(*solver),
                sasktran2::rust::successive_orders::final_residual(*solver));
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::capture_forward_state(int wavelidx,
                                                      int threadidx) {
        if (!m_capture_forward_state || wavelidx < 0 ||
            static_cast<std::size_t>(wavelidx) >= m_forward_states.size()) {
            return;
        }
        auto& state = m_forward_states[static_cast<std::size_t>(wavelidx)];
        const auto& solver = m_rust_solvers[threadidx];
        if (!sasktran2::rust::successive_orders::converged(*solver)) {
            state.valid = false;
            return;
        }

        const auto& storage = m_thread_storage[threadidx];
        if (m_altitude_grid.interpolation_method() !=
            sasktran2::grids::interpolation::lower) {
            const auto forcing =
                sasktran2::rust::successive_orders::first_order_forcing(
                    *solver);
            state.first_order_forcing.assign(forcing.begin(), forcing.end());
        } else {
            const auto& forcing = storage.m_firstorder_radiances.value;
            state.first_order_forcing.assign(forcing.data(),
                                             forcing.data() + forcing.size());
        }
        const auto solution =
            sasktran2::rust::successive_orders::solution(*solver);
        state.solution.assign(solution.begin(), solution.end());
        state.valid = true;
    }

    template <int NSTOKES>
    void
    DiffuseTable<NSTOKES>::pack_rust_end_of_ray_source(int wavelidx,
                                                       int threadidx) const {
        auto& packed = m_thread_storage[threadidx].rust_end_of_ray_source;
        const int num_rays =
            static_cast<int>(m_internal_viewing.traced_rays.size());
        if (packed.size() != static_cast<std::size_t>(num_rays * NSTOKES)) {
            throw std::logic_error("Rust end-of-ray source size mismatch");
        }
        const int num_threads = m_config->num_source_threads();
        std::vector<sasktran2::WavelengthBlockDual<NSTOKES>> scratch(
            num_threads);
        for (auto& source_value : scratch) {
            source_value.resize(1, 0, true);
        }
        const sasktran2::WavelengthBlock<> block{wavelidx, 1};
#pragma omp parallel for num_threads(num_threads)
        for (int rayidx = 0; rayidx < num_rays; ++rayidx) {
#ifdef SKTRAN_OPENMP_SUPPORT
            const int source_threadidx =
                num_threads == 1 ? 0 : omp_get_thread_num();
#else
            const int source_threadidx = 0;
#endif
            auto& source_value = scratch[source_threadidx];
            source_value.set_zero(1);
            for (const auto* source : m_initial_sources) {
                source->end_of_ray_source(block, rayidx, threadidx,
                                          source_threadidx + threadidx,
                                          source_value);
            }
            for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                packed[rayidx * NSTOKES + stokes] =
                    source_value.value(stokes, 0);
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::pack_rust_end_of_ray_source_jvp(
        int wavelidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        auto& storage = m_thread_storage[threadidx];
        const int num_rays =
            static_cast<int>(m_internal_viewing.traced_rays.size());
        const int num_threads = m_config->num_source_threads();
        std::vector<sasktran2::RadianceJVP<NSTOKES>> scratch(num_threads);
#pragma omp parallel for num_threads(num_threads)
        for (int rayidx = 0; rayidx < num_rays; ++rayidx) {
#ifdef SKTRAN_OPENMP_SUPPORT
            const int source_threadidx =
                num_threads == 1 ? 0 : omp_get_thread_num();
#else
            const int source_threadidx = 0;
#endif
            auto& source_value = scratch[source_threadidx];
            source_value.set_zero();
            for (const auto* source : m_initial_sources) {
                source->end_of_ray_source_jvp(wavelidx, rayidx, threadidx,
                                              source_threadidx + threadidx,
                                              native_tangent, source_value);
            }
            for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                storage.rust_end_of_ray_source[rayidx * NSTOKES + stokes] =
                    source_value.value(stokes);
                storage
                    .rust_end_of_ray_source_tangent[rayidx * NSTOKES + stokes] =
                    source_value.jvp(stokes);
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::accumulate_rust_end_of_ray_source_vjp(
        int wavelidx, int threadidx,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        const auto& packed_gradient =
            m_thread_storage[threadidx].rust_end_of_ray_source_gradient;
        const int num_rays =
            static_cast<int>(m_internal_viewing.traced_rays.size());
        const int num_threads = m_config->num_source_threads();
        std::vector<Eigen::VectorXd> thread_gradients(num_threads);
        for (auto& gradient : thread_gradients) {
            gradient.setZero(native_gradient.size());
        }
#pragma omp parallel for num_threads(num_threads)
        for (int rayidx = 0; rayidx < num_rays; ++rayidx) {
#ifdef SKTRAN_OPENMP_SUPPORT
            const int source_threadidx =
                num_threads == 1 ? 0 : omp_get_thread_num();
#else
            const int source_threadidx = 0;
#endif
            Eigen::Vector<double, NSTOKES> cotangent;
            for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                cotangent(stokes) = packed_gradient[rayidx * NSTOKES + stokes];
            }
            for (const auto* source : m_initial_sources) {
                source->end_of_ray_source_vjp(
                    wavelidx, rayidx, threadidx, source_threadidx + threadidx,
                    cotangent, thread_gradients[source_threadidx]);
            }
        }
        for (const auto& gradient : thread_gradients) {
            native_gradient += gradient;
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::assemble_rust_ray_transport(
        int wavelidx, int threadidx, bool assemble_first_order) {
        if (m_rust_ray_transports.empty() || m_atmosphere == nullptr) {
            throw std::logic_error(
                "Rust transport geometry is not initialized");
        }
        auto& storage = m_thread_storage[threadidx];
        const auto& atmosphere_storage = m_atmosphere->storage();
        const std::size_t num_geometry = static_cast<std::size_t>(
            atmosphere_storage.total_extinction.rows());
        const ::rust::Slice<const double> extinction(
            atmosphere_storage.total_extinction.col(wavelidx).data(),
            num_geometry);
        const ::rust::Slice<const double> single_scatter_albedo(
            atmosphere_storage.ssa.col(wavelidx).data(), num_geometry);
        const bool rust_owned_transport =
            m_altitude_grid.interpolation_method() !=
            sasktran2::grids::interpolation::lower;
        const bool retain_layer_scratch =
            m_active_vjp_derivatives[threadidx].prepared;

        if (!assemble_first_order) {
            if (rust_owned_transport) {
                sasktran2::rust::successive_orders::
                    assemble_solver_ray_transport(
                        *m_rust_ray_transports.front(),
                        *m_rust_solvers[threadidx], extinction,
                        single_scatter_albedo, retain_layer_scratch);
                return;
            }
            const std::size_t num_layers =
                sasktran2::rust::successive_orders::ray_transport_num_layers(
                    *m_rust_ray_transports.front());
            storage.rust_layer_optical_depth.resize(num_layers);
            storage.rust_layer_attenuation.resize(num_layers);
            storage.rust_layer_prefix_attenuation.resize(num_layers);
            sasktran2::rust::successive_orders::assemble_ray_transport(
                *m_rust_ray_transports.front(), extinction,
                single_scatter_albedo,
                as_rust_mut_slice(storage.accumulation_summed_values),
                as_rust_mut_slice(storage.rust_layer_optical_depth),
                as_rust_mut_slice(storage.rust_layer_attenuation),
                as_rust_mut_slice(storage.rust_layer_prefix_attenuation),
                as_rust_mut_slice(storage.rust_ray_end_attenuation));
            return;
        }
        if (!rust_owned_transport) {
            throw std::logic_error(
                "Fused first-order transport requires Rust-owned storage");
        }

        const Eigen::VectorXd* solar_transmission = nullptr;
        if (m_geometry_1d != nullptr) {
            const auto* single_scatter =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionTable,
                    NSTOKES>*>(m_initial_sources.front());
            if (single_scatter != nullptr) {
                solar_transmission =
                    &single_scatter->solar_transmission(threadidx);
            }
        } else {
            const auto* single_scatter =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionExact,
                    NSTOKES>*>(m_initial_sources.front());
            if (single_scatter != nullptr) {
                solar_transmission =
                    &single_scatter->solar_transmission(threadidx);
            }
        }
        if (solar_transmission == nullptr || m_initial_sources.size() != 1) {
            throw std::logic_error(
                "Rust first-order forcing requires one single-scatter source");
        }

        const std::size_t coefficient_families = NSTOKES == 1 ? 1 : 4;
        const std::size_t num_phase_moments =
            static_cast<std::size_t>(atmosphere_storage.max_stored_legendre());
        if (num_phase_moments !=
            static_cast<std::size_t>(m_config->num_singlescatter_moments())) {
            throw std::logic_error("Rust phase storage size mismatch");
        }
        const ::rust::Slice<const double> legendre_coefficients(
            &atmosphere_storage.leg_coeff(0, 0, wavelidx),
            num_geometry * num_phase_moments * coefficient_families);
        static_assert(sizeof(int) == sizeof(std::int32_t));
        const ::rust::Slice<const std::int32_t> maximum_order(
            reinterpret_cast<const std::int32_t*>(
                atmosphere_storage.max_order.col(wavelidx).data()),
            num_geometry);
        pack_rust_end_of_ray_source(wavelidx, threadidx);

        sasktran2::rust::successive_orders::
            assemble_solver_ray_transport_with_first_order(
                *m_rust_ray_transports.front(), *m_rust_solvers[threadidx],
                extinction, single_scatter_albedo, legendre_coefficients,
                maximum_order, as_rust_slice(*solar_transmission),
                as_rust_slice(storage.rust_end_of_ray_source),
                retain_layer_scratch);
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::assemble_rust_ray_transport_jvp(
        int wavelidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        if (m_rust_ray_transports.empty() || m_atmosphere == nullptr) {
            throw std::logic_error(
                "Rust transport geometry is not initialized");
        }
        auto& output = m_thread_storage[threadidx];
        const auto& atmosphere = m_atmosphere->storage();
        const std::size_t num_geometry =
            static_cast<std::size_t>(atmosphere.total_extinction.rows());
        const std::size_t coefficient_families = NSTOKES == 1 ? 1 : 4;
        const std::size_t num_phase_moments =
            static_cast<std::size_t>(atmosphere.max_stored_legendre());
        const std::size_t num_phase_values =
            num_geometry * num_phase_moments * coefficient_families;
        const ::rust::Slice<const double> extinction(
            atmosphere.total_extinction.col(wavelidx).data(), num_geometry);
        const ::rust::Slice<const double> single_scatter_albedo(
            atmosphere.ssa.col(wavelidx).data(), num_geometry);
        const ::rust::Slice<const double> legendre_coefficients(
            &atmosphere.leg_coeff(0, 0, wavelidx), num_phase_values);
        static_assert(sizeof(int) == sizeof(std::int32_t));
        const ::rust::Slice<const std::int32_t> maximum_order(
            reinterpret_cast<const std::int32_t*>(
                atmosphere.max_order.col(wavelidx).data()),
            num_geometry);

        const Eigen::VectorXd* solar_transmission = nullptr;
        const Eigen::VectorXd* solar_transmission_jvp = nullptr;
        if (m_geometry_1d != nullptr) {
            const auto* single_scatter =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionTable,
                    NSTOKES>*>(m_initial_sources.front());
            if (single_scatter != nullptr) {
                solar_transmission =
                    &single_scatter->solar_transmission(threadidx);
                solar_transmission_jvp =
                    &single_scatter->solar_transmission_jvp(threadidx);
            }
        } else {
            const auto* single_scatter =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionExact,
                    NSTOKES>*>(m_initial_sources.front());
            if (single_scatter != nullptr) {
                solar_transmission =
                    &single_scatter->solar_transmission(threadidx);
                solar_transmission_jvp =
                    &single_scatter->solar_transmission_jvp(threadidx);
            }
        }
        if (solar_transmission == nullptr ||
            solar_transmission_jvp == nullptr) {
            throw std::logic_error(
                "Rust first-order JVP requires a single-scatter source");
        }

        auto& legendre_tangent = output.rust_single_scatter_legendre_tangent;
        legendre_tangent.resize(num_phase_values);
        const auto& coefficient_derivatives = atmosphere.d_leg_coeff;
        const std::size_t coefficient_stride =
            static_cast<std::size_t>(coefficient_derivatives.dimension(0));
        const std::size_t num_wavelengths =
            static_cast<std::size_t>(coefficient_derivatives.dimension(2));
        const std::size_t num_derivative_groups =
            static_cast<std::size_t>(coefficient_derivatives.dimension(3));
        const double* scattering_tangent =
            native_tangent.data() + m_atmosphere->scat_deriv_start_index();
        sasktran2::rust::successive_orders::
            calculate_atmospheric_coefficient_jvp(
                {coefficient_derivatives.data(),
                 static_cast<std::size_t>(coefficient_derivatives.size())},
                coefficient_stride, num_geometry, num_wavelengths,
                static_cast<std::size_t>(wavelidx),
                {scattering_tangent, num_geometry * num_derivative_groups},
                as_rust_mut_slice(legendre_tangent));

        const ::rust::Slice<const double> extinction_tangent(
            native_tangent.data(), num_geometry);
        const ::rust::Slice<const double> albedo_tangent(
            native_tangent.data() + m_atmosphere->ssa_deriv_start_index(),
            num_geometry);
        pack_rust_end_of_ray_source_jvp(wavelidx, threadidx, native_tangent);
        sasktran2::rust::successive_orders::
            assemble_solver_ray_transport_with_first_order_jvp(
                *m_rust_ray_transports.front(), *m_rust_solvers[threadidx],
                extinction, single_scatter_albedo, legendre_coefficients,
                maximum_order, as_rust_slice(*solar_transmission),
                extinction_tangent, albedo_tangent,
                as_rust_slice(legendre_tangent),
                as_rust_slice(*solar_transmission_jvp),
                as_rust_slice(output.rust_end_of_ray_source),
                as_rust_slice(output.rust_end_of_ray_source_tangent));
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::accumulate_rust_transport_vjp_with_first_order(
        int wavelidx, int threadidx,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        const auto& output = m_thread_storage[threadidx];
        const auto& activity = m_active_vjp_derivatives[threadidx];
        pack_rust_end_of_ray_source(wavelidx, threadidx);

        const auto& atmosphere = m_atmosphere->storage();
        const std::size_t num_geometry =
            static_cast<std::size_t>(atmosphere.total_extinction.rows());
        const std::size_t coefficient_families = NSTOKES == 1 ? 1 : 4;
        const std::size_t num_phase_moments =
            static_cast<std::size_t>(atmosphere.max_stored_legendre());
        const std::size_t num_phase_values =
            num_geometry * num_phase_moments * coefficient_families;
        const Eigen::VectorXd* solar_transmission = nullptr;
        const sasktran2::solartransmission::SingleScatterSource<
            sasktran2::solartransmission::SolarTransmissionTable, NSTOKES>*
            table_source = nullptr;
        const sasktran2::solartransmission::SingleScatterSource<
            sasktran2::solartransmission::SolarTransmissionExact, NSTOKES>*
            exact_source = nullptr;
        if (m_geometry_1d != nullptr) {
            table_source = dynamic_cast<
                const sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionTable,
                    NSTOKES>*>(m_initial_sources.front());
            if (table_source != nullptr) {
                solar_transmission =
                    &table_source->solar_transmission(threadidx);
            }
        } else {
            exact_source = dynamic_cast<
                const sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionExact,
                    NSTOKES>*>(m_initial_sources.front());
            if (exact_source != nullptr) {
                solar_transmission =
                    &exact_source->solar_transmission(threadidx);
            }
        }
        if (solar_transmission == nullptr) {
            throw std::logic_error(
                "Rust first-order VJP requires a single-scatter source");
        }

        output.rust_single_scatter_extinction_gradient.resize(num_geometry);
        output.rust_single_scatter_albedo_gradient.resize(num_geometry);
        output.rust_single_scatter_legendre_gradient.resize(num_phase_values);
        output.rust_single_scatter_solar_gradient.resize(
            solar_transmission->size());
        const ::rust::Slice<const double> extinction(
            atmosphere.total_extinction.col(wavelidx).data(), num_geometry);
        const ::rust::Slice<const double> single_scatter_albedo(
            atmosphere.ssa.col(wavelidx).data(), num_geometry);
        const ::rust::Slice<const double> legendre_coefficients(
            &atmosphere.leg_coeff(0, 0, wavelidx), num_phase_values);
        static_assert(sizeof(int) == sizeof(std::int32_t));
        const ::rust::Slice<const std::int32_t> maximum_order(
            reinterpret_cast<const std::int32_t*>(
                atmosphere.max_order.col(wavelidx).data()),
            num_geometry);
        sasktran2::rust::successive_orders::
            assemble_solver_ray_transport_with_first_order_vjp(
                *m_rust_ray_transports.front(), *m_rust_solvers[threadidx],
                extinction, single_scatter_albedo, legendre_coefficients,
                maximum_order, as_rust_slice(*solar_transmission),
                as_rust_slice(output.rust_end_of_ray_source),
                as_rust_mut_slice(
                    output.rust_single_scatter_extinction_gradient),
                as_rust_mut_slice(output.rust_single_scatter_albedo_gradient),
                as_rust_mut_slice(output.rust_single_scatter_legendre_gradient),
                as_rust_mut_slice(output.rust_single_scatter_solar_gradient),
                as_rust_mut_slice(output.rust_end_of_ray_source_gradient));

        accumulate_rust_end_of_ray_source_vjp(wavelidx, threadidx,
                                              native_gradient);

        if (activity.extinction) {
            for (std::size_t location = 0; location < num_geometry;
                 ++location) {
                native_gradient(location) +=
                    output.rust_single_scatter_extinction_gradient[location];
            }
        }
        if (activity.ssa) {
            for (std::size_t location = 0; location < num_geometry;
                 ++location) {
                native_gradient(m_atmosphere->ssa_deriv_start_index() +
                                static_cast<int>(location)) +=
                    output.rust_single_scatter_albedo_gradient[location];
            }
        }
        if (!activity.scattering_groups.empty()) {
            const auto& coefficient_derivatives = atmosphere.d_leg_coeff;
            const std::size_t coefficient_stride =
                static_cast<std::size_t>(coefficient_derivatives.dimension(0));
            const std::size_t num_wavelengths =
                static_cast<std::size_t>(coefficient_derivatives.dimension(2));
            const std::size_t num_derivative_groups =
                static_cast<std::size_t>(coefficient_derivatives.dimension(3));
            static_assert(sizeof(int) == sizeof(std::int32_t));
            sasktran2::rust::successive_orders::
                accumulate_atmospheric_coefficient_gradient(
                    as_rust_slice(output.rust_single_scatter_legendre_gradient),
                    {coefficient_derivatives.data(),
                     static_cast<std::size_t>(coefficient_derivatives.size())},
                    coefficient_stride, num_geometry, num_wavelengths,
                    static_cast<std::size_t>(wavelidx),
                    {reinterpret_cast<const std::int32_t*>(
                         activity.scattering_groups.data()),
                     activity.scattering_groups.size()},
                    {native_gradient.data() +
                         m_atmosphere->scat_deriv_start_index(),
                     num_geometry * num_derivative_groups});
        }

        const Eigen::Map<const Eigen::VectorXd> solar_gradient(
            output.rust_single_scatter_solar_gradient.data(),
            output.rust_single_scatter_solar_gradient.size());
        if (table_source != nullptr) {
            table_source->accumulate_solar_transmission_vjp(threadidx, 0,
                                                            solar_gradient);
        } else {
            exact_source->accumulate_solar_transmission_vjp(threadidx, 0,
                                                            solar_gradient);
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::restore_transport_operator(int wavelidx,
                                                           int threadidx) {
        auto& storage = m_thread_storage[threadidx];
        if (!m_rust_ray_transports.empty()) {
            assemble_rust_ray_transport(wavelidx, threadidx);
            return;
        }

        storage.accumulation_summed_values.setZero();
        const int num_threads = m_config->num_source_threads();
#pragma omp parallel for num_threads(num_threads)
        for (int rayidx = 0; rayidx < m_internal_viewing.traced_rays.size();
             ++rayidx) {
            m_integrator.emplace_accumulation_transport(
                wavelidx, rayidx, m_diffuse_source_weights,
                storage.accumulation_summed_values);
        }
    }

#endif

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::interpolate_sources(
        const Eigen::VectorXd& old_outgoing,
        sasktran2::Dual<double>& new_outgoing) {
        for (int i = 0; i < m_diffuse_points.size(); ++i) {
            if (!m_diffuse_point_full_calculation[i]) {
                // Have to interpolate sources from the above/below diffuse
                // points that do have a full calculation

                // Find the point below
                int lower_idx = i - 1;
                while (!m_diffuse_point_full_calculation[lower_idx]) {
                    --lower_idx;
                }

                int upper_idx = i + 1;
                // Find the point above
                while (!m_diffuse_point_full_calculation[upper_idx]) {
                    ++upper_idx;
                }

                double alt_above =
                    m_diffuse_points[upper_idx]->location().radius();
                double alt_below =
                    m_diffuse_points[lower_idx]->location().radius();

                double alt = m_diffuse_points[i]->location().radius();

                double w_above = (alt - alt_below) / (alt_above - alt_below);
                double w_below = 1 - w_above;

                // TODO: These adjustments work well for I, but more work is
                // needed on adjusting Q/U...
                auto above_seq = Eigen::seq(
                    NSTOKES * m_diffuse_outgoing_index_map[upper_idx],
                    NSTOKES * (m_diffuse_outgoing_index_map[upper_idx] +
                               m_diffuse_points[upper_idx]->num_outgoing()) -
                        1,
                    NSTOKES);
                auto below_seq = Eigen::seq(
                    NSTOKES * m_diffuse_outgoing_index_map[lower_idx],
                    NSTOKES * (m_diffuse_outgoing_index_map[lower_idx] +
                               m_diffuse_points[lower_idx]->num_outgoing()) -
                        1,
                    NSTOKES);

                for (int s = 0; s < NSTOKES; ++s) {
                    auto seq = Eigen::seq(
                        NSTOKES * m_diffuse_outgoing_index_map[i] + s,
                        NSTOKES * (m_diffuse_outgoing_index_map[i] +
                                   m_diffuse_points[i]->num_outgoing()) -
                            1 + s,
                        NSTOKES);
                    new_outgoing.value.array()(seq) *=
                        w_above * (new_outgoing.value.array()(above_seq) /
                                   old_outgoing.array()(above_seq)) +
                        w_below * (new_outgoing.value.array()(below_seq) /
                                   old_outgoing.array()(below_seq));
                }
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::prepare_wavelength(int wavelidx, int threadidx,
                                                   bool value_only) {
        const sasktran2::WavelengthBlock<> scalar_block{wavelidx, 1};
        for (auto& source : m_initial_owned_sources) {
            if (value_only) {
                source->calculate_value(scalar_block, threadidx);
            } else {
                source->calculate(scalar_block, threadidx);
            }
        }
        int nthreads = m_config->num_source_threads();

        if (m_config->num_hr_spherical_iterations() > 0) {
            // Calculate the first order incoming signal
            std::vector<
                sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>>
                temp_result;
            temp_result.resize(nthreads);
            bool cpp_first_order_required = true;

            m_thread_storage[threadidx].accumulation_summed_values.setZero();

#ifdef SKTRAN_RUST_SUPPORT
            bool rust_transport_assembled = false;
            bool rust_first_order_assembled = false;
            if (m_use_rust_solver && !m_rust_ray_transports.empty()) {
                rust_first_order_assembled =
                    m_altitude_grid.interpolation_method() !=
                    sasktran2::grids::interpolation::lower;
                assemble_rust_ray_transport(wavelidx, threadidx,
                                            rust_first_order_assembled);
                rust_transport_assembled = true;
                cpp_first_order_required = !rust_first_order_assembled;
            }
#endif

            if (cpp_first_order_required) {
#pragma omp parallel for num_threads(nthreads)
                for (int rayidx = 0;
                     rayidx < m_internal_viewing.traced_rays.size(); ++rayidx) {
                    int ray_threadidx;
                    if (nthreads == 1) {
                        ray_threadidx = 0;
                    } else {
#ifdef SKTRAN_OPENMP_SUPPORT
                        ray_threadidx = omp_get_thread_num();
#else
                        ray_threadidx = 0;
#endif
                    }

                    temp_result[ray_threadidx].value.setZero();
#ifdef SKTRAN_RUST_SUPPORT
                    if (rust_transport_assembled) {
                        const auto& storage = m_thread_storage[threadidx];
                        m_integrator.integrate_first_order_precomputed(
                            temp_result[ray_threadidx], m_initial_sources,
                            wavelidx, rayidx, threadidx,
                            ray_threadidx + threadidx,
                            m_rust_transport_ray_layer_offsets[rayidx],
                            storage.rust_layer_optical_depth,
                            storage.rust_layer_attenuation,
                            storage.rust_layer_prefix_attenuation,
                            storage.rust_ray_end_attenuation[rayidx]);
                    } else
#endif
                    {
                        m_integrator
                            .integrate_and_emplace_accumulation_triplets(
                                temp_result[ray_threadidx], m_initial_sources,
                                wavelidx, rayidx, threadidx,
                                ray_threadidx + threadidx,
                                m_diffuse_source_weights,
                                m_thread_storage[threadidx]
                                    .accumulation_summed_values);
                    }

                    auto first_order_segment =
                        m_thread_storage[threadidx]
                            .m_firstorder_radiances.value(
                                Eigen::seq(rayidx * NSTOKES,
                                           rayidx * NSTOKES + NSTOKES - 1));
                    first_order_segment = temp_result[ray_threadidx].value;

#ifdef SASKTRAN_DEBUG_ASSERTS
                    if (temp_result[ray_threadidx].value.hasNaN()) {
                        spdlog::error("Incoming Ray: {} has NaN", rayidx);
                    }
#endif
                }
            }

            // The Rust-owned path consumes first-order forcing directly. The
            // legacy path still stages it in the incoming-radiance buffer.
#ifdef SKTRAN_RUST_SUPPORT
            if (!rust_first_order_assembled)
#endif
            {
                m_thread_storage[threadidx].m_incoming_radiances.value =
                    m_thread_storage[threadidx].m_firstorder_radiances.value;
            }
            // Generate the scattering and accumulation matrices
            generate_scattering_matrices(wavelidx, threadidx);
            generate_accumulation_matrix(wavelidx, threadidx);
        }
    }

    template <int NSTOKES>
    void
    DiffuseTable<NSTOKES>::calculate(const sasktran2::WavelengthBlock<>& block,
                                     int threadidx) {
        m_deferred_jvp_transport_restore[threadidx] = false;
        if (block.count < 1 || block.count > maximum_wavelength_block_size() ||
            block.count > m_wavelength_block_capacity) {
            throw std::invalid_argument(
                "Invalid successive-orders wavelength block");
        }

#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            auto& storage = m_thread_storage[threadidx];
            if (m_wavelength_block_capacity == 1) {
                prepare_wavelength(block.start, threadidx);
                iterate_to_solution(block.start, threadidx);
                capture_forward_state(block.start, threadidx);
                return;
            }
            const std::size_t state_size = m_num_outgoing_values;
            for (int lane = 0; lane < block.count; ++lane) {
                prepare_wavelength(block.wavelength(lane), threadidx);
                iterate_to_solution(block.wavelength(lane), threadidx);
                capture_forward_state(block.wavelength(lane), threadidx);
                const auto solution =
                    sasktran2::rust::successive_orders::solution(
                        *m_rust_solvers[threadidx]);
                if (solution.size() != state_size) {
                    throw std::logic_error(
                        "Rust batched successive-orders solution size "
                        "mismatch");
                }
                for (std::size_t element = 0; element < state_size; ++element) {
                    storage.rust_batch_outgoing_sources
                        [element * m_wavelength_block_capacity + lane] =
                        solution[element];
                }
            }
            return;
        }
#endif

        prepare_wavelength(block.start, threadidx);
        iterate_to_solution(block.start, threadidx);
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::calculate_value(
        const sasktran2::WavelengthBlock<>& block, int threadidx) {
        m_deferred_jvp_transport_restore[threadidx] = false;
        if (block.count != 1) {
            calculate(block, threadidx);
            return;
        }
        prepare_wavelength(block.start, threadidx, true);
        iterate_to_solution(block.start, threadidx);
#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            capture_forward_state(block.start, threadidx);
        }
#endif
    }

    template <int NSTOKES>
    bool DiffuseTable<NSTOKES>::restore_forward_state(int wavelidx,
                                                      int threadidx) {
#ifdef SKTRAN_RUST_SUPPORT
        m_deferred_jvp_transport_restore[threadidx] = false;
        if (!m_use_rust_solver || m_forward_state_atmosphere != m_atmosphere ||
            m_forward_state_atmosphere_revision != m_atmosphere->revision() ||
            wavelidx < 0 ||
            static_cast<std::size_t>(wavelidx) >= m_forward_states.size() ||
            !m_forward_states[static_cast<std::size_t>(wavelidx)].valid) {
            return false;
        }

        const auto& state =
            m_forward_states[static_cast<std::size_t>(wavelidx)];
        auto& storage = m_thread_storage[threadidx];
        const bool rust_owned_transport =
            m_altitude_grid.interpolation_method() !=
            sasktran2::grids::interpolation::lower;
        const std::size_t expected_forcing_size =
            rust_owned_transport
                ? m_internal_viewing.traced_rays.size() * NSTOKES
                : static_cast<std::size_t>(
                      storage.m_firstorder_radiances.value.size());
        if (state.first_order_forcing.size() != expected_forcing_size ||
            state.solution.size() != m_num_outgoing_values) {
            return false;
        }

        const sasktran2::WavelengthBlock<> scalar_block{wavelidx, 1};
        for (auto& source : m_initial_owned_sources) {
            source->calculate(scalar_block, threadidx);
        }
        restore_transport_operator(wavelidx, threadidx);
        generate_scattering_matrices(wavelidx, threadidx);

        if (!rust_owned_transport) {
            std::copy(state.first_order_forcing.begin(),
                      state.first_order_forcing.end(),
                      storage.m_firstorder_radiances.value.data());
        }
        if (rust_owned_transport) {
            sasktran2::rust::successive_orders::
                restore_coefficient_solution_assembled(
                    *m_rust_solvers[threadidx],
                    as_rust_slice(state.first_order_forcing),
                    as_rust_slice(state.solution));
        } else {
            sasktran2::rust::successive_orders::restore_coefficient_solution(
                *m_rust_solvers[threadidx],
                as_rust_slice(storage.rust_boundary_scattering_values),
                as_rust_slice(state.solution));
        }
        return true;
#else
        (void)wavelidx;
        (void)threadidx;
        return false;
#endif
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::prepare_forward_state_for_jvp(
        int wavelidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        for (auto* source : m_initial_sources) {
            source->prepare_forward_state_for_jvp(wavelidx, threadidx,
                                                  native_tangent);
        }
    }

    template <int NSTOKES>
    bool DiffuseTable<NSTOKES>::restore_forward_state_for_jvp(
        int wavelidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
#ifdef SKTRAN_RUST_SUPPORT
        m_deferred_jvp_transport_restore[threadidx] = false;
        if (!m_use_rust_solver || m_forward_state_atmosphere != m_atmosphere ||
            m_forward_state_atmosphere_revision != m_atmosphere->revision() ||
            wavelidx < 0 ||
            static_cast<std::size_t>(wavelidx) >= m_forward_states.size() ||
            !m_forward_states[static_cast<std::size_t>(wavelidx)].valid) {
            return false;
        }

        const auto& state =
            m_forward_states[static_cast<std::size_t>(wavelidx)];
        auto& storage = m_thread_storage[threadidx];
        const bool rust_owned_transport =
            m_altitude_grid.interpolation_method() !=
            sasktran2::grids::interpolation::lower;
        const std::size_t expected_forcing_size =
            rust_owned_transport
                ? m_internal_viewing.traced_rays.size() * NSTOKES
                : static_cast<std::size_t>(
                      storage.m_firstorder_radiances.value.size());
        if (state.first_order_forcing.size() != expected_forcing_size ||
            state.solution.size() != m_num_outgoing_values) {
            return false;
        }

        const sasktran2::WavelengthBlock<> scalar_block{wavelidx, 1};
        for (auto& source : m_initial_owned_sources) {
            if (native_tangent.isZero(0.0)) {
                source->calculate_value(scalar_block, threadidx);
            } else {
                source->calculate(scalar_block, threadidx);
            }
        }
        if (!rust_owned_transport) {
            std::copy(state.first_order_forcing.begin(),
                      state.first_order_forcing.end(),
                      storage.m_firstorder_radiances.value.data());
        }
        sasktran2::rust::successive_orders::restore_solution_only(
            *m_rust_solvers[threadidx], as_rust_slice(state.solution));

        // A zero tangent needs only retained source values for LOS
        // integration. For a nonzero tangent, scattering is inexpensive to
        // restore here while transport reconstruction is fused with its JVP.
        if (!native_tangent.isZero(0.0)) {
            generate_scattering_matrices(wavelidx, threadidx);
            m_deferred_jvp_transport_restore[threadidx] = true;
        }
        return true;
#else
        (void)wavelidx;
        (void)threadidx;
        (void)native_tangent;
        return false;
#endif
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::prepare_jvp(
        int wavelidx, int wavel_threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
#ifdef SKTRAN_RUST_SUPPORT
        if (!native_products_available()) {
            throw std::logic_error(
                "Native successive-orders products require the Rust solver");
        }
        if (m_config->successive_orders_initialization() ==
                sasktran2::Config::SuccessiveOrdersInitialization::twostream &&
            m_config->successive_orders_relative_tolerance() == 0 &&
            m_config->successive_orders_absolute_tolerance() == 0 &&
            m_config->successive_orders_anderson_depth() == 0) {
            throw std::logic_error(
                "Fixed-iteration derivatives with a two-stream "
                "successive-orders initial state are not yet supported");
        }
        auto& storage = m_thread_storage[wavel_threadidx];
        const bool rust_owned_transport =
            !m_rust_ray_transports.empty() &&
            m_altitude_grid.interpolation_method() !=
                sasktran2::grids::interpolation::lower;
        if (rust_owned_transport) {
            storage.rust_transport_value_tangent.resize(0);
        } else {
            storage.rust_transport_value_tangent.resize(
                storage.accumulation_summed_values.size());
        }
        storage.rust_boundary_scattering_value_tangent.resize(
            storage.rust_boundary_scattering_values.size());
        if (rust_owned_transport) {
            storage.rust_first_order_forcing_tangent.resize(0);
        } else {
            storage.rust_first_order_forcing_tangent.resize(
                storage.m_firstorder_radiances.value.size());
        }
        storage.rust_transport_value_tangent.setZero();
        storage.rust_first_order_forcing_tangent.setZero();
        const bool restore_transport =
            m_deferred_jvp_transport_restore[wavel_threadidx];
        m_deferred_jvp_transport_restore[wavel_threadidx] = false;
        if (native_tangent.isZero(0.0)) {
            sasktran2::rust::successive_orders::clear_solution_jvp(
                *m_rust_solvers[wavel_threadidx]);
            for (auto* source : m_initial_sources) {
                source->prepare_jvp(wavelidx, wavel_threadidx, native_tangent);
            }
            return;
        }
        if (restore_transport) {
            storage.accumulation_summed_values.setZero();
        }
        pack_rust_scattering_coefficient_jvp(wavelidx, wavel_threadidx,
                                             native_tangent);
        pack_rust_boundary_scattering_jvp(wavelidx, wavel_threadidx,
                                          native_tangent);
        for (auto* source : m_initial_sources) {
            source->prepare_jvp(wavelidx, wavel_threadidx, native_tangent);
        }
        if (!sasktran2::rust::successive_orders::converged(
                *m_rust_solvers[wavel_threadidx]) &&
            m_config->successive_orders_anderson_depth() != 0) {
            spdlog::warn(
                "Calculating an approximate successive-orders JVP for "
                "wavelength index {} because the primal did not converge "
                "after {} iterations (residual {}). The result is based on "
                "the current approximate solution and may not represent the "
                "converged model.",
                wavelidx,
                sasktran2::rust::successive_orders::iterations(
                    *m_rust_solvers[wavel_threadidx]),
                sasktran2::rust::successive_orders::final_residual(
                    *m_rust_solvers[wavel_threadidx]));
        }

        if (rust_owned_transport) {
            assemble_rust_ray_transport_jvp(wavelidx, wavel_threadidx,
                                            native_tangent);
        } else {
            const int num_threads = m_config->num_source_threads();
#pragma omp parallel for num_threads(num_threads)
            for (int rayidx = 0; rayidx < m_internal_viewing.traced_rays.size();
                 ++rayidx) {
#ifdef SKTRAN_OPENMP_SUPPORT
                const int source_threadidx =
                    num_threads == 1 ? 0 : omp_get_thread_num();
#else
                const int source_threadidx = 0;
#endif
                m_integrator.integrate_and_emplace_accumulation_jvp(
                    m_initial_sources, wavelidx, rayidx, wavel_threadidx,
                    source_threadidx + wavel_threadidx,
                    m_diffuse_source_weights, native_tangent,
                    restore_transport ? &storage.accumulation_summed_values
                                      : nullptr,
                    storage.rust_transport_value_tangent,
                    storage.rust_first_order_forcing_tangent);
            }
        }
        auto& solver = m_rust_solvers[wavel_threadidx];
        if (restore_transport) {
            const auto& state =
                m_forward_states[static_cast<std::size_t>(wavelidx)];
            if (rust_owned_transport) {
                sasktran2::rust::successive_orders::
                    restore_coefficient_solution_assembled(
                        *solver, as_rust_slice(state.first_order_forcing),
                        as_rust_slice(state.solution));
            } else {
                sasktran2::rust::successive_orders::
                    restore_coefficient_solution(
                        *solver,
                        as_rust_slice(storage.rust_boundary_scattering_values),
                        as_rust_slice(state.solution));
            }
        }
        if (rust_owned_transport) {
            sasktran2::rust::successive_orders::
                linearize_coefficients_jvp_assembled(*solver);
        } else {
            sasktran2::rust::successive_orders::linearize_coefficients_jvp(
                *solver, as_rust_slice(m_outer_starts),
                as_rust_slice(m_inner_indicies),
                as_rust_slice(storage.accumulation_summed_values),
                as_rust_slice(storage.rust_transport_value_tangent),
                as_rust_slice(storage.rust_boundary_scattering_value_tangent),
                as_rust_slice(storage.m_firstorder_radiances.value),
                as_rust_slice(storage.rust_first_order_forcing_tangent));
        }
        const auto solution_jvp =
            sasktran2::rust::successive_orders::solution_jvp(*solver);
        if (solution_jvp.size() != m_num_outgoing_values) {
            throw std::runtime_error(
                "Rust successive-orders JVP returned an invalid size");
        }
#else
        (void)wavelidx;
        (void)wavel_threadidx;
        (void)native_tangent;
        throw std::logic_error("Rust support is disabled");
#endif
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::prepare_vjp(
        int wavelidx, int wavel_threadidx,
        Eigen::Ref<const Eigen::VectorXi> active_derivatives) {
        if (!native_products_available()) {
            throw std::logic_error(
                "Native successive-orders products require the Rust solver");
        }
        if (m_config->successive_orders_initialization() ==
                sasktran2::Config::SuccessiveOrdersInitialization::twostream &&
            m_config->successive_orders_relative_tolerance() == 0 &&
            m_config->successive_orders_absolute_tolerance() == 0 &&
            m_config->successive_orders_anderson_depth() == 0) {
            throw std::logic_error(
                "Fixed-iteration derivatives with a two-stream "
                "successive-orders initial state are not yet supported");
        }
        auto& activity = m_active_vjp_derivatives[wavel_threadidx];
        activity = DiffuseTableDerivativeActivity{};
        activity.prepared = true;
        const int num_geometry =
            static_cast<int>(m_atmosphere->storage().ssa.rows());
        activity.extinction = active_derivatives.head(num_geometry).any();
        activity.ssa =
            active_derivatives
                .segment(m_atmosphere->ssa_deriv_start_index(), num_geometry)
                .any();
        for (int group = 0; group < m_atmosphere->num_scattering_deriv_groups();
             ++group) {
            if (active_derivatives
                    .segment(m_atmosphere->scat_deriv_start_index() +
                                 group * num_geometry,
                             num_geometry)
                    .any()) {
                activity.scattering_groups.push_back(group);
            }
        }
        const int num_surface_derivatives = m_atmosphere->surface().num_deriv();
        activity.surface =
            num_surface_derivatives > 0 &&
            active_derivatives
                .segment(m_atmosphere->surface_deriv_start_index(),
                         num_surface_derivatives)
                .any();
        if (m_atmosphere->include_emission_derivatives()) {
            activity.emission =
                active_derivatives
                    .segment(m_atmosphere->emission_deriv_start_index(),
                             num_geometry)
                    .any();
            activity.surface_emission =
                active_derivatives(
                    m_atmosphere->surface_emission_deriv_start_index()) != 0;
        }
        for (auto* source : m_initial_sources) {
            source->prepare_vjp(wavelidx, wavel_threadidx, active_derivatives);
        }
        for (auto& cotangent : m_los_solution_cotangents[wavel_threadidx]) {
            cotangent.setZero();
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::finalize_vjp(
        int wavelidx, int wavel_threadidx,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
#ifdef SKTRAN_RUST_SUPPORT
        Eigen::VectorXd solution_cotangent = Eigen::VectorXd::Zero(
            static_cast<Eigen::Index>(m_num_outgoing_values));
        for (const auto& thread_cotangent :
             m_los_solution_cotangents[wavel_threadidx]) {
            solution_cotangent += thread_cotangent;
        }
        const auto& storage = m_thread_storage[wavel_threadidx];
        auto& solver = m_rust_solvers[wavel_threadidx];
        if (!sasktran2::rust::successive_orders::converged(*solver) &&
            m_config->successive_orders_anderson_depth() != 0 &&
            m_last_nonconverged_vjp_warning_wavelength[wavel_threadidx] !=
                wavelidx) {
            spdlog::warn(
                "Calculating an approximate successive-orders VJP for "
                "wavelength index {} because the primal did not converge "
                "after {} iterations (residual {}). The result is based on "
                "the current approximate solution and may not represent the "
                "converged model.",
                wavelidx,
                sasktran2::rust::successive_orders::iterations(*solver),
                sasktran2::rust::successive_orders::final_residual(*solver));
            m_last_nonconverged_vjp_warning_wavelength[wavel_threadidx] =
                wavelidx;
        }
        const bool rust_owned_transport =
            !m_rust_ray_transports.empty() &&
            m_altitude_grid.interpolation_method() !=
                sasktran2::grids::interpolation::lower;
        if (rust_owned_transport) {
            sasktran2::rust::successive_orders::
                linearize_coefficients_vjp_assembled(
                    *solver, as_rust_slice(solution_cotangent));
        } else {
            sasktran2::rust::successive_orders::linearize_coefficients_vjp(
                *solver, as_rust_slice(m_outer_starts),
                as_rust_slice(m_inner_indicies),
                as_rust_slice(storage.accumulation_summed_values),
                as_rust_slice(storage.m_firstorder_radiances.value),
                as_rust_slice(solution_cotangent));
        }
        const auto transport_gradient =
            sasktran2::rust::successive_orders::transport_value_gradient(
                *solver);
        const auto boundary_gradient = sasktran2::rust::successive_orders::
            boundary_scattering_value_gradient(*solver);
        const auto forcing_gradient =
            sasktran2::rust::successive_orders::first_order_forcing_gradient(
                *solver);
        const auto converged_solution =
            sasktran2::rust::successive_orders::solution(*solver);

        const auto& activity = m_active_vjp_derivatives[wavel_threadidx];
        if (activity.scattering()) {
            accumulate_rust_scattering_coefficient_vjp(
                wavelidx, wavel_threadidx, activity.scattering_groups,
                native_gradient);
        }
        if (activity.surface) {
            accumulate_rust_boundary_scattering_vjp(wavelidx, boundary_gradient,
                                                    native_gradient);
        }
        if (rust_owned_transport) {
            accumulate_rust_transport_vjp_with_first_order(
                wavelidx, wavel_threadidx, native_gradient);
        } else {
            const int num_threads = m_config->num_source_threads();
            std::vector<Eigen::VectorXd> thread_gradients(num_threads);
            for (auto& gradient : thread_gradients) {
                gradient.setZero(native_gradient.size());
            }
#pragma omp parallel for num_threads(num_threads)
            for (int rayidx = 0; rayidx < m_internal_viewing.traced_rays.size();
                 ++rayidx) {
#ifdef SKTRAN_OPENMP_SUPPORT
                const int source_threadidx =
                    num_threads == 1 ? 0 : omp_get_thread_num();
#else
                const int source_threadidx = 0;
#endif
                m_integrator.accumulate_accumulation_vjp(
                    m_initial_sources, wavelidx, rayidx, wavel_threadidx,
                    source_threadidx + wavel_threadidx,
                    m_diffuse_source_weights,
                    Eigen::Map<const Eigen::VectorXd>(
                        transport_gradient.data(), transport_gradient.size()),
                    Eigen::Map<const Eigen::VectorXd>(
                        converged_solution.data(), converged_solution.size()),
                    m_inner_indicies,
                    Eigen::Map<const Eigen::VectorXd>(forcing_gradient.data(),
                                                      forcing_gradient.size()),
                    activity.extinction, activity.ssa,
                    activity.extinction || activity.ssa ||
                        activity.scattering() || activity.emission,
                    true, thread_gradients[source_threadidx]);
            }
            for (const auto& gradient : thread_gradients) {
                native_gradient += gradient;
            }
        }
        for (const auto* source : m_initial_sources) {
            source->finalize_vjp(wavelidx, wavel_threadidx, native_gradient);
        }
        m_active_vjp_derivatives[wavel_threadidx] =
            DiffuseTableDerivativeActivity{};
#else
        (void)wavelidx;
        (void)wavel_threadidx;
        (void)native_gradient;
        throw std::logic_error("Rust support is disabled");
#endif
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::integrated_source(
        const sasktran2::WavelengthBlock<>& block, int losidx, int layeridx,
        int wavel_threadidx, int threadidx,
        const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        sasktran2::WavelengthBlockDual<NSTOKES>& block_source,
        typename SourceTermInterface<NSTOKES>::IntegrationDirection direction)
        const {
        auto& storage = m_thread_storage[wavel_threadidx];

#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            throw std::logic_error(
                "Rust successive-orders LOS integration must use the "
                "whole-ray callback");
        }
#endif

        const auto& ray_interpolator = m_los_source_weights[losidx];
        const auto& layer_interpolator =
            ray_interpolator.interior_weights[layeridx];

        for (int lane = 0; lane < block.count; ++lane) {
            const int wavelidx = block.wavelength(lane);
            sasktran2::WavelengthBlockLaneDualView<NSTOKES> source(block_source,
                                                                   lane);

            // Start by calculating ssa at the source point.
            double omega = 0;
            for (std::size_t i = layer_interpolator.atmosphere_offset;
                 i < layer_interpolator.atmosphere_offset +
                         layer_interpolator.atmosphere_count;
                 ++i) {
                const int atmosphere_index =
                    entrance_weights[i - layer_interpolator.atmosphere_offset]
                        .first;
                omega +=
                    m_atmosphere->storage().ssa(atmosphere_index, wavelidx) *
                    ray_interpolator.atmosphere_weights[i];
            }

            const double source_factor = 1 - std::exp(-shell_od.od(lane));

            for (std::size_t i = layer_interpolator.source_offset;
                 i < layer_interpolator.source_offset +
                         layer_interpolator.source_count;
                 ++i) {
                const int source_index = ray_interpolator.source_indices[i];
                const double interpolation_weight =
                    ray_interpolator.source_weights[i];

                for (int s = 0; s < NSTOKES; ++s) {
                    const std::size_t outgoing_index =
                        static_cast<std::size_t>(source_index * NSTOKES + s);
                    double outgoing_value;
                    if (m_use_rust_solver && m_wavelength_block_capacity > 1) {
                        outgoing_value =
                            storage.rust_batch_outgoing_sources
                                [outgoing_index * m_wavelength_block_capacity +
                                 lane];
                    } else {
                        outgoing_value =
                            storage.m_outgoing_sources.value(outgoing_index);
                    }
                    const double source_value =
                        outgoing_value * interpolation_weight;

                    source.value(s) += omega * source_factor * source_value;

                    if (source.derivative_size() > 0) {
                        // Now we need dJ/dthickness.
                        for (auto it = shell_od.derivative_iterator(); it;
                             ++it) {
                            source.deriv(s, it.index()) += it.value() *
                                                           (1 - source_factor) *
                                                           source_value * omega;
                        }

                        // And dJ/dssa.
                        for (std::size_t atmosphere_index =
                                 layer_interpolator.atmosphere_offset;
                             atmosphere_index <
                             layer_interpolator.atmosphere_offset +
                                 layer_interpolator.atmosphere_count;
                             ++atmosphere_index) {
                            const int native_atmosphere_index =
                                entrance_weights[atmosphere_index -
                                                 layer_interpolator
                                                     .atmosphere_offset]
                                    .first;
                            source.deriv(s,
                                         m_atmosphere->ssa_deriv_start_index() +
                                             native_atmosphere_index) +=
                                ray_interpolator
                                    .atmosphere_weights[atmosphere_index] *
                                source_factor * source_value;
                        }

                        if (this->m_config->wf_precision() ==
                                sasktran2::Config::WeightingFunctionPrecision::
                                    full &&
                            m_config->initialize_hr_with_do()) {
                            source.deriv(s, Eigen::placeholders::all) +=
                                omega * source_factor * interpolation_weight *
                                storage.m_outgoing_sources.deriv(
                                    outgoing_index, Eigen::placeholders::all);
                        }
                    }
                }
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::integrated_source_jvp(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        sasktran2::RadianceJVP<NSTOKES>& source) const {
        (void)threadidx;
        (void)layer;
        (void)exit_weights;
        const auto& storage = m_thread_storage[wavel_threadidx];
#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            throw std::logic_error(
                "Rust successive-orders LOS JVP must use the whole-ray "
                "callback");
        }
#endif
        const auto& interpolator = m_los_source_weights[losidx];
        const auto& layer_interpolator =
            interpolator.interior_weights[layeridx];
        double omega = 0.0;
        double omega_jvp = 0.0;
        for (std::size_t index = layer_interpolator.atmosphere_offset;
             index < layer_interpolator.atmosphere_offset +
                         layer_interpolator.atmosphere_count;
             ++index) {
            const int atmosphere_index =
                entrance_weights[index - layer_interpolator.atmosphere_offset]
                    .first;
            const double weight = interpolator.atmosphere_weights[index];
            omega += weight *
                     m_atmosphere->storage().ssa(atmosphere_index, wavelidx);
            omega_jvp +=
                weight * native_tangent(m_atmosphere->ssa_deriv_start_index() +
                                        atmosphere_index);
        }
        double shell_od_jvp = 0.0;
        for (auto derivative = shell_od.derivative_iterator(); derivative;
             ++derivative) {
            shell_od_jvp +=
                derivative.value() * native_tangent(derivative.index());
        }
        const double attenuation = shell_od.exp_minus_od(0);
        const double factor = 1.0 - attenuation;
        for (std::size_t index = layer_interpolator.source_offset;
             index <
             layer_interpolator.source_offset + layer_interpolator.source_count;
             ++index) {
            const int source_index = interpolator.source_indices[index];
            const double weight = interpolator.source_weights[index];
            for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                const std::size_t outgoing_index =
                    static_cast<std::size_t>(source_index * NSTOKES + stokes);
                const double outgoing =
                    storage.m_outgoing_sources.value(outgoing_index);
                const double outgoing_jvp =
                    storage.rust_solution_jvp[outgoing_index];
                source.value(stokes) += omega * factor * weight * outgoing;
                source.jvp(stokes) +=
                    weight *
                    ((omega_jvp * factor + omega * attenuation * shell_od_jvp) *
                         outgoing +
                     omega * factor * outgoing_jvp);
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::end_of_ray_source_jvp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        sasktran2::RadianceJVP<NSTOKES>& source) const {
        (void)wavelidx;
        (void)threadidx;
        (void)native_tangent;
        const auto& storage = m_thread_storage[wavel_threadidx];
#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            return;
        }
#endif
        const auto& interpolator = m_los_source_weights[losidx];
        for (std::size_t index = 0;
             index < interpolator.ground_source_indices.size(); ++index) {
            const int source_index = interpolator.ground_source_indices[index];
            const double weight = interpolator.ground_source_weights[index];
            for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                const std::size_t outgoing_index =
                    static_cast<std::size_t>(source_index * NSTOKES + stokes);
                source.value(stokes) +=
                    weight * storage.m_outgoing_sources.value(outgoing_index);
                source.jvp(stokes) +=
                    weight * storage.rust_solution_jvp[outgoing_index];
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::integrated_source_vjp(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        const Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::Vector<double, NSTOKES>> source_value,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        (void)layer;
        (void)exit_weights;
        const auto& storage = m_thread_storage[wavel_threadidx];
#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            throw std::logic_error(
                "Rust successive-orders LOS VJP must use the whole-ray "
                "callback");
        }
#endif
        const auto& interpolator = m_los_source_weights[losidx];
        const auto& layer_interpolator =
            interpolator.interior_weights[layeridx];
        double omega = 0.0;
        for (std::size_t index = layer_interpolator.atmosphere_offset;
             index < layer_interpolator.atmosphere_offset +
                         layer_interpolator.atmosphere_count;
             ++index) {
            const int atmosphere_index =
                entrance_weights[index - layer_interpolator.atmosphere_offset]
                    .first;
            omega += interpolator.atmosphere_weights[index] *
                     m_atmosphere->storage().ssa(atmosphere_index, wavelidx);
        }
        Eigen::Vector<double, NSTOKES> interpolated_outgoing =
            Eigen::Vector<double, NSTOKES>::Zero();
        const double factor = 1.0 - shell_od.exp_minus_od(0);
        for (std::size_t index = layer_interpolator.source_offset;
             index <
             layer_interpolator.source_offset + layer_interpolator.source_count;
             ++index) {
            const int source_index = interpolator.source_indices[index];
            const double weight = interpolator.source_weights[index];
            for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                const std::size_t outgoing_index =
                    static_cast<std::size_t>(source_index * NSTOKES + stokes);
                const double outgoing =
                    storage.m_outgoing_sources.value(outgoing_index);
                interpolated_outgoing(stokes) += weight * outgoing;
                m_los_solution_cotangents[wavel_threadidx]
                                         [threadidx - wavel_threadidx](
                                             outgoing_index) +=
                    weight * omega * factor * cotangent(stokes);
            }
        }
        source_value += omega * factor * interpolated_outgoing;
        const double amplitude_cotangent = cotangent.dot(interpolated_outgoing);
        for (std::size_t index = layer_interpolator.atmosphere_offset;
             index < layer_interpolator.atmosphere_offset +
                         layer_interpolator.atmosphere_count;
             ++index) {
            const int atmosphere_index =
                entrance_weights[index - layer_interpolator.atmosphere_offset]
                    .first;
            native_gradient(m_atmosphere->ssa_deriv_start_index() +
                            atmosphere_index) +=
                interpolator.atmosphere_weights[index] * factor *
                amplitude_cotangent;
        }
        const double shell_od_cotangent =
            omega * shell_od.exp_minus_od(0) * amplitude_cotangent;
        for (auto derivative = shell_od.derivative_iterator(); derivative;
             ++derivative) {
            native_gradient(derivative.index()) +=
                derivative.value() * shell_od_cotangent;
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::end_of_ray_source_vjp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        const Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        (void)wavelidx;
        (void)native_gradient;
#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            return;
        }
#endif
        const auto& interpolator = m_los_source_weights[losidx];
        for (std::size_t index = 0;
             index < interpolator.ground_source_indices.size(); ++index) {
            const int source_index = interpolator.ground_source_indices[index];
            const double weight = interpolator.ground_source_weights[index];
            for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                m_los_solution_cotangents[wavel_threadidx]
                                         [threadidx - wavel_threadidx](
                                             source_index * NSTOKES + stokes) +=
                    weight * cotangent(stokes);
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::end_of_ray_source(
        const sasktran2::WavelengthBlock<>& block, int losidx,
        int wavel_threadidx, int threadidx,
        sasktran2::WavelengthBlockDual<NSTOKES>& block_source) const {
        auto& storage = m_thread_storage[wavel_threadidx];

#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            return;
        }
#endif

        const auto& interpolator = m_los_source_weights[losidx];

        for (int lane = 0; lane < block.count; ++lane) {
            sasktran2::WavelengthBlockLaneDualView<NSTOKES> source(block_source,
                                                                   lane);
            for (std::size_t i = 0;
                 i < interpolator.ground_source_indices.size(); ++i) {
                const int source_index = interpolator.ground_source_indices[i];
                const double interpolation_weight =
                    interpolator.ground_source_weights[i];

                for (int s = 0; s < NSTOKES; ++s) {
                    const std::size_t outgoing_index =
                        static_cast<std::size_t>(source_index * NSTOKES + s);
                    double outgoing_value;
                    if (m_use_rust_solver && m_wavelength_block_capacity > 1) {
                        outgoing_value =
                            storage.rust_batch_outgoing_sources
                                [outgoing_index * m_wavelength_block_capacity +
                                 lane];
                    } else {
                        outgoing_value =
                            storage.m_outgoing_sources.value(outgoing_index);
                    }
                    source.value(s) += outgoing_value * interpolation_weight;

                    if (this->m_config->wf_precision() ==
                            sasktran2::Config::WeightingFunctionPrecision::
                                full &&
                        m_config->initialize_hr_with_do()) {
                        source.deriv(s, Eigen::placeholders::all) +=
                            interpolation_weight *
                            storage.m_outgoing_sources.deriv(
                                outgoing_index, Eigen::placeholders::all);
                    }
                }
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::start_of_ray_source(
        const sasktran2::WavelengthBlock<>& block, int losidx,
        int wavel_threadidx, int threadidx,
        sasktran2::WavelengthBlockDual<NSTOKES>& block_source) const {
        (void)threadidx;
#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            if (m_rust_source_interpolators.empty()) {
                throw std::logic_error(
                    "Rust LOS source interpolator is not initialized");
            }
            const auto& atmosphere = m_atmosphere->storage();
            const std::size_t num_atmosphere_points =
                static_cast<std::size_t>(atmosphere.total_extinction.rows());
            const std::size_t atmosphere_values =
                num_atmosphere_points * static_cast<std::size_t>(block.count);
            const bool batched = m_wavelength_block_capacity > 1;
            const auto solution =
                batched ? as_rust_slice(m_thread_storage[wavel_threadidx]
                                            .rust_batch_outgoing_sources)
                        : sasktran2::rust::successive_orders::solution(
                              *m_rust_solvers[wavel_threadidx]);
            const std::size_t solution_stride =
                batched ? static_cast<std::size_t>(m_wavelength_block_capacity)
                        : 1;
            auto output = ::rust::Slice<double>{
                block_source.value.data(),
                static_cast<std::size_t>(block_source.value.size())};
            const auto extinction = ::rust::Slice<const double>{
                atmosphere.total_extinction.col(block.start).data(),
                atmosphere_values};
            const auto single_scatter_albedo = ::rust::Slice<const double>{
                atmosphere.ssa.col(block.start).data(), atmosphere_values};

            if (block_source.derivative_size() == 0) {
                sasktran2::rust::successive_orders::accumulate_source_ray_value(
                    *m_rust_source_interpolators.front(), losidx, extinction,
                    single_scatter_albedo, block.count, solution,
                    solution_stride, output, block_source.block_capacity());
            } else {
                sasktran2::rust::successive_orders::
                    accumulate_source_ray_jacobian(
                        *m_rust_source_interpolators.front(), losidx,
                        extinction, single_scatter_albedo, block.count,
                        solution, solution_stride, output,
                        block_source.block_capacity(),
                        {block_source.deriv.data(),
                         static_cast<std::size_t>(block_source.deriv.size())},
                        block_source.derivative_size(),
                        block_source.block_capacity(),
                        m_atmosphere->ssa_deriv_start_index());
            }
            return;
        }
#else
        (void)block;
        (void)losidx;
        (void)wavel_threadidx;
        (void)block_source;
#endif
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::start_of_ray_source_jvp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        sasktran2::RadianceJVP<NSTOKES>& source) const {
        (void)threadidx;
#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            if (m_rust_source_interpolators.empty()) {
                throw std::logic_error(
                    "Rust LOS source interpolator is not initialized");
            }
            const auto& atmosphere = m_atmosphere->storage();
            const std::size_t num_atmosphere_points =
                static_cast<std::size_t>(atmosphere.total_extinction.rows());
            const std::size_t ssa_derivative_start =
                static_cast<std::size_t>(m_atmosphere->ssa_deriv_start_index());
            sasktran2::rust::successive_orders::accumulate_source_ray_jvp(
                *m_rust_source_interpolators.front(), losidx,
                {atmosphere.total_extinction.col(wavelidx).data(),
                 num_atmosphere_points},
                {atmosphere.ssa.col(wavelidx).data(), num_atmosphere_points},
                {native_tangent.data(), num_atmosphere_points},
                {native_tangent.data() + ssa_derivative_start,
                 num_atmosphere_points},
                sasktran2::rust::successive_orders::solution(
                    *m_rust_solvers[wavel_threadidx]),
                sasktran2::rust::successive_orders::solution_jvp(
                    *m_rust_solvers[wavel_threadidx]),
                {source.value.data(), NSTOKES}, {source.jvp.data(), NSTOKES});
            return;
        }
#else
        (void)wavelidx;
        (void)losidx;
        (void)wavel_threadidx;
        (void)native_tangent;
        (void)source;
#endif
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::start_of_ray_source_vjp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        const Eigen::Vector<double, NSTOKES>& value_before,
        Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        (void)value_before;
#ifdef SKTRAN_RUST_SUPPORT
        if (m_use_rust_solver) {
            if (m_rust_source_interpolators.empty()) {
                throw std::logic_error(
                    "Rust LOS source interpolator is not initialized");
            }
            const auto& atmosphere = m_atmosphere->storage();
            const std::size_t num_atmosphere_points =
                static_cast<std::size_t>(atmosphere.total_extinction.rows());
            const std::size_t ssa_derivative_start =
                static_cast<std::size_t>(m_atmosphere->ssa_deriv_start_index());
            auto& solution_cotangent =
                m_los_solution_cotangents[wavel_threadidx]
                                         [threadidx - wavel_threadidx];
            sasktran2::rust::successive_orders::accumulate_source_ray_vjp(
                *m_rust_source_interpolators.front(), losidx,
                {atmosphere.total_extinction.col(wavelidx).data(),
                 num_atmosphere_points},
                {atmosphere.ssa.col(wavelidx).data(), num_atmosphere_points},
                sasktran2::rust::successive_orders::solution(
                    *m_rust_solvers[wavel_threadidx]),
                {cotangent.data(), NSTOKES},
                as_rust_mut_slice(solution_cotangent),
                {native_gradient.data() + ssa_derivative_start,
                 num_atmosphere_points},
                {native_gradient.data(), num_atmosphere_points});
            return;
        }
#else
        (void)wavelidx;
        (void)losidx;
        (void)wavel_threadidx;
        (void)threadidx;
        (void)cotangent;
        (void)native_gradient;
#endif
    }

    template class DiffuseTable<1>;
    template class DiffuseTable<3>;

} // namespace sasktran2::hr
