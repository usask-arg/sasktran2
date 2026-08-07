#include "sasktran2/geometry.h"
#include "sasktran2/viewinggeometry_internal.h"
#include <algorithm>
#include <cmath>
#ifdef SKTRAN_OPENMP_SUPPORT
#include <omp.h>
#endif
#include "source_internal.h"
#include "solar_transmission.h"
#include <sasktran2/math/unitsphere.h>
#include <sasktran2/math/scattering.h>
#include <sasktran2/solartransmission.h>

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

namespace sasktran2::successive_orders {

    template <int NSTOKES>
    SourceImpl<NSTOKES>::SourceImpl(
        const sasktran2::raytracing::RayTracerBase& ray_tracer,
        const sasktran2::Geometry1D& geometry)
        : m_raytracer(ray_tracer), m_geometry(geometry),
          m_altitude_grid(geometry.altitude_grid()), m_geometry_1d(&geometry),
          m_geometry_2d(nullptr)
#ifdef SKTRAN_RUST_SUPPORT
          ,
          m_raytracer_2d(nullptr)
#endif
    {
    }

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    SourceImpl<NSTOKES>::SourceImpl(
        const sasktran2::raytracing::RustRayTracer2D& ray_tracer,
        const sasktran2::Geometry2D& geometry)
        : m_raytracer(ray_tracer), m_geometry(geometry),
          m_altitude_grid(geometry.altitude_grid()), m_geometry_1d(nullptr),
          m_geometry_2d(&geometry), m_raytracer_2d(&ray_tracer) {}
#endif

    template <int NSTOKES>
    sasktran2::grids::Grid
    SourceImpl<NSTOKES>::generate_cos_sza_grid(double min_cos_sza,
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
    SourceImpl<NSTOKES>::generate_altitude_grid() {
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
    void SourceImpl<NSTOKES>::construct_diffuse_points() {
        // Create the sphere pairs
        // One for the atmosphere points and one for each of the ground points
        // for now

        int num_interior_spheres = 1;
        m_unit_sphere_pairs.resize(
            num_interior_spheres +
            m_location_interpolator->num_ground_points());

        // TODO: Get npoints from config, figure out if we want to use more than
        // one type of sphere
        m_unit_sphere_pairs[0] =
            std::make_unique<IncomingOutgoingSpherePair<NSTOKES>>(
                m_config->num_do_streams(),
                std::move(std::make_unique<sasktran2::math::LebedevSphere>(
                    m_config->num_hr_incoming())),
                std::move(std::make_unique<sasktran2::math::LebedevSphere>(
                    m_config->num_hr_outgoing())),
                false);

        // TODO: Same number of streams for the ground term? probably... to be
        // figured out when BRDF is implemented
        for (int i = 0; i < m_location_interpolator->num_ground_points(); ++i) {
            Eigen::Vector3d location = m_location_interpolator->ground_location(
                m_geometry.coordinates(), i);

            m_unit_sphere_pairs[num_interior_spheres + i] =
                std::make_unique<IncomingOutgoingSpherePair<NSTOKES>>(
                    m_config->num_do_streams(),
                    std::make_unique<sasktran2::math::UnitSphereGround>(
                        std::move(
                            std::make_unique<sasktran2::math::LebedevSphere>(
                                m_config->num_hr_incoming())),
                        location),
                    std::make_unique<sasktran2::math::UnitSphereGround>(
                        std::move(
                            std::make_unique<sasktran2::math::LebedevSphere>(
                                m_config->num_hr_outgoing())),
                        location),
                    false);
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

            point = std::make_unique<DiffusePoint<NSTOKES>>(
                *m_unit_sphere_pairs[0], loc);
        }

        for (int i = 0; i < m_location_interpolator->num_ground_points(); ++i) {
            auto& point = m_diffuse_points[i + m_location_interpolator
                                                   ->num_interior_points()];

            loc.position = m_location_interpolator->ground_location(
                m_geometry.coordinates(), i);

            // Add 0.01m to the ground location to avoid rounding errors
            loc.position += loc.position.normalized() * 0.01;

            point = std::make_unique<DiffusePoint<NSTOKES>>(
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
            if (m_config->successive_orders_initialization() ==
                sasktran2::Config::SuccessiveOrdersInitialization::twostream) {
                storage.twostream_initial_guess.resize(start_outgoing_idx *
                                                       NSTOKES);
            }
            if (m_wavelength_block_capacity > 1) {
                storage.rust_batch_outgoing_sources.resize(
                    start_outgoing_idx * NSTOKES * m_wavelength_block_capacity);
            }
        }
    }

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    void SourceImpl<NSTOKES>::initialize_twostream_initial_guess_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        if (m_config->successive_orders_initialization() !=
            sasktran2::Config::SuccessiveOrdersInitialization::twostream) {
            return;
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
    std::vector<std::vector<int>> SourceImpl<NSTOKES>::trace_incoming_rays() {
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
                if (m_geometry_2d != nullptr) {
                    m_transport_geometry.compact_2d_ray(rayidx, traced_ray);
                }
            }
        }
        return transport_columns;
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        ZoneScopedN("Initialize HR Geometry");
        if (!internal_viewing.flux_observers.empty()) {
            throw std::invalid_argument(
                "The Rust successive-orders source currently supports "
                "radiances only, not flux observers");
        }
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
            m_location_interpolator =
                std::make_unique<AltitudeHorizontalSourceLocationInterpolator>(
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

        const bool use_incremental_compact_geometry = m_geometry_2d != nullptr;
        if (use_incremental_compact_geometry) {
            m_transport_geometry.begin_compact_2d(
                m_internal_viewing.traced_rays, m_geometry);
        }

        // Trace all of the incoming rays
        auto transport_columns = trace_incoming_rays();

        // Set up the integrator
        if (use_incremental_compact_geometry) {
            m_transport_geometry.finalize_compact_2d();
        } else {
            m_transport_geometry.initialize(m_internal_viewing.traced_rays);
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
        construct_accumulation_sparsity(transport_columns);
        std::vector<std::vector<int>>().swap(transport_columns);
#ifdef SKTRAN_RUST_SUPPORT
        release_cpp_transport_geometry();
#endif
        generate_source_interpolation_weights(internal_viewing.traced_rays,
                                              m_los_source_weights);
#ifdef SKTRAN_RUST_SUPPORT
        initialize_rust_source_interpolator(internal_viewing);
#endif
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_atmosphere = &atmosphere;

        for (auto& source : m_initial_owned_sources) {
            source->initialize_atmosphere(atmosphere);
        }
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::initialize_atmosphere_native(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_atmosphere = &atmosphere;
        for (auto& source : m_initial_owned_sources) {
            source->initialize_atmosphere_native(atmosphere);
        }
    }

    template <int NSTOKES>
    void
    SourceImpl<NSTOKES>::set_wavelength_block_capacity(int block_capacity) {
        if (block_capacity < 1 ||
            block_capacity > maximum_wavelength_block_size()) {
            throw std::invalid_argument(
                "Invalid successive-orders wavelength block capacity");
        }
        m_wavelength_block_capacity = block_capacity;
        int initial_source_block_capacity = 1;
#ifdef SKTRAN_RUST_SUPPORT
        if constexpr (NSTOKES == 1) {
            if (block_capacity == 4 &&
                m_config->successive_orders_initialization() ==
                    sasktran2::Config::SuccessiveOrdersInitialization::none) {
                initial_source_block_capacity = block_capacity;
            }
        }
#endif
        for (auto& source : m_initial_owned_sources) {
            source->set_wavelength_block_capacity(
                initial_source_block_capacity);
        }

        for (auto& storage : m_thread_storage) {
            if (block_capacity == 1) {
                std::vector<double>().swap(storage.rust_batch_outgoing_sources);
                std::vector<double>().swap(
                    storage.rust_batch_end_of_ray_source);
            } else {
                storage.rust_batch_outgoing_sources.resize(
                    m_num_outgoing_values * block_capacity);
                storage.rust_batch_end_of_ray_source.resize(
                    storage.rust_end_of_ray_source.size() * block_capacity);
            }
        }
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::begin_forward_state_capture(int num_wavelengths) {
        if (num_wavelengths < 1) {
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
    void SourceImpl<NSTOKES>::end_forward_state_capture() {
        m_capture_forward_state = false;
    }

    template <int NSTOKES>
    std::size_t SourceImpl<NSTOKES>::retained_forward_state_bytes() const {
        std::size_t bytes = m_forward_states.capacity() * sizeof(ForwardState);
        for (const auto& state : m_forward_states) {
            bytes += state.first_order_forcing.capacity() * sizeof(double);
            bytes += state.solution.capacity() * sizeof(double);
        }
        return bytes;
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::construct_accumulation_sparsity(
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

#ifdef SKTRAN_RUST_SUPPORT
        {
            m_rust_transport_sparsity.clear();
            m_rust_transport_sparsity.push_back(
                sasktran2::rust::successive_orders::new_transport_sparsity(
                    static_cast<std::size_t>(m_outer_starts.size() - 1),
                    m_num_outgoing_values, as_rust_slice(m_outer_starts),
                    as_rust_slice(m_inner_indicies)));
            spdlog::debug(
                "Packed Rust successive-orders transport sparsity: {} rows, "
                "{} values, {:.3f} MiB shared across wavelength workers",
                m_outer_starts.size() - 1, m_inner_indicies.size(),
                static_cast<double>(
                    sasktran2::rust::successive_orders::
                        transport_sparsity_storage_bytes(
                            *m_rust_transport_sparsity.front())) /
                    (1024.0 * 1024.0));
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
                auto packed_geometry = m_transport_geometry.pack(
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
                    packed_geometry.ray_layer_offsets.size() != num_rays + 1) {
                    throw std::logic_error(
                        "Rust transport geometry size mismatch");
                }
                for (auto& storage : m_thread_storage) {
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
                }
            }
            m_outer_starts.resize(0);
            m_inner_indicies.resize(0);
            if (m_config->successive_orders_initialization() !=
                sasktran2::Config::SuccessiveOrdersInitialization::twostream) {
                for (auto& storage : m_thread_storage) {
                    storage.twostream_initial_guess.resize(0);
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
    SourceImpl<NSTOKES>::initialize_config(const sasktran2::Config& config) {
        m_config = &config;
        if (m_config->apply_delta_scaling()) {
            throw std::invalid_argument(
                "SuccessiveOrdersRust does not support delta-M scaling");
        }
        if (m_altitude_grid.interpolation_method() ==
            sasktran2::grids::interpolation::lower) {
            throw std::invalid_argument(
                "SuccessiveOrdersRust requires linear altitude-grid "
                "interpolation; lower interpolation is only supported by the "
                "legacy SuccessiveOrders source");
        }
        if (m_config->num_hr_full_incoming_points() > 0) {
            throw std::invalid_argument(
                "The Rust successive-orders source currently requires all "
                "diffuse points to have explicitly traced incoming rays");
        }
        if (m_config->initialize_hr_with_do()) {
            throw std::invalid_argument(
                "The Rust successive-orders source does not yet support a "
                "discrete-ordinates initial guess");
        }

        m_thread_storage.resize(m_config->num_wavelength_threads());
        m_deferred_jvp_transport_restore.assign(
            m_config->num_wavelength_threads(), false);
        m_batch_jvp_zero_tangent.assign(m_config->num_wavelength_threads(),
                                        false);
        m_last_nonconverged_vjp_warning_wavelength.assign(
            m_config->num_wavelength_threads(), -1);
        m_active_vjp_derivatives.resize(m_config->num_wavelength_threads());
        m_batch_vjp_derivative_activity.resize(
            m_config->num_wavelength_threads());
        m_batch_vjp_active_derivatives.resize(
            m_config->num_wavelength_threads());

        if (m_geometry_1d != nullptr) {
            auto initial_single_scatter = std::make_unique<
                sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionTable,
                    NSTOKES>>(*m_geometry_1d, m_raytracer);
            initial_single_scatter->enable_table_native_products();
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

        if (m_geometry_1d != nullptr) {
            auto* source =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionTable,
                    NSTOKES>*>(m_initial_sources.front());
            if (source == nullptr) {
                throw std::logic_error(
                    "Rust source requires table solar transmission in "
                    "Geometry1D");
            }
            source->delegate_interior_source(
                [](const Eigen::MatrixXd& dense_path,
                   const SolarTransmission::SparseMatrix& interpolation,
                   const std::vector<bool>& ground_hit, std::size_t input_size,
                   int wavelength_workers) {
                    return SolarTransmission::create_table(
                        dense_path, interpolation, ground_hit, input_size,
                        wavelength_workers);
                });
        } else {
            auto* source =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionExact,
                    NSTOKES>*>(m_initial_sources.front());
            if (source == nullptr) {
                throw std::logic_error(
                    "Rust source requires exact solar transmission in "
                    "Geometry2D");
            }
            source->delegate_interior_source(
                [](const Eigen::MatrixXd& dense_path,
                   const SolarTransmission::SparseMatrix& sparse_path,
                   const std::vector<bool>& ground_hit, std::size_t input_size,
                   int wavelength_workers) {
                    return SolarTransmission::create_exact(
                        dense_path, sparse_path, ground_hit, input_size,
                        wavelength_workers);
                });
        }

        for (auto& source : m_initial_owned_sources) {
            source->initialize_config(config);
        }
    }

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    void SourceImpl<NSTOKES>::release_cpp_transport_geometry() {
        if (m_rust_ray_transports.empty()) {
            return;
        }
        SInterpolator().swap(m_diffuse_source_weights);

        m_transport_geometry.release();

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
    Eigen::Vector3d SourceImpl<NSTOKES>::rotate_unit_vector(
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
    void SourceImpl<NSTOKES>::generate_source_interpolation_weights(
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
    void SourceImpl<NSTOKES>::generate_source_interpolation_weights(
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
    void SourceImpl<NSTOKES>::initialize_rust_source_interpolator(
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
    void SourceImpl<NSTOKES>::compile_accumulation_row(
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
    void SourceImpl<NSTOKES>::generate_scattering_matrices(int wavelidx,
                                                           int threadidx) {
        ZoneScopedN("Generate Scattering Matrices");
        auto& storage = m_thread_storage[threadidx];

#ifdef SKTRAN_RUST_SUPPORT
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
        const std::size_t wavelength_size = coefficient_stride * num_geometry;
        const double* wavelength_coefficients =
            coefficients.data() +
            static_cast<std::size_t>(wavelidx) * wavelength_size;
        sasktran2::rust::successive_orders::
            set_scattering_coefficients_from_atmosphere(
                *m_rust_solvers[threadidx],
                *m_rust_scattering_interpolators.front(),
                {wavelength_coefficients, wavelength_size}, coefficient_stride);
        auto boundary_values =
            sasktran2::rust::successive_orders::boundary_scattering_values_mut(
                *m_rust_solvers[threadidx]);
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
                m_location_interpolator->num_interior_points() + ground_index);
            if (!m_diffuse_point_full_calculation[point_index]) {
                continue;
            }
            const auto& point = m_diffuse_points[point_index];
            point->sphere_pair().calculate_ground_scattering_values(
                m_atmosphere->surface(), point->location(), wavelidx,
                boundary_values.data() +
                    m_rust_boundary_scattering_offsets[ground_index]);
        }
#else
        (void)wavelidx;
        (void)storage;
        throw std::logic_error(
            "The Rust successive-orders source requires Rust support");
#endif
    }

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    void SourceImpl<NSTOKES>::generate_scattering_matrices_batch4(
        const sasktran2::WavelengthBlock<>& block, int threadidx) {
        if constexpr (NSTOKES != 1) {
            throw std::logic_error(
                "Four-wavelength scattering requires scalar mode");
        }
        if (block.count != 4 || m_rust_scattering_interpolators.empty()) {
            throw std::logic_error(
                "Four-wavelength Rust scattering is not initialized");
        }
        const auto& coefficients = m_atmosphere->storage().leg_coeff;
        const std::size_t coefficient_stride =
            static_cast<std::size_t>(coefficients.dimension(0));
        const std::size_t num_geometry =
            static_cast<std::size_t>(coefficients.dimension(1));
        if (block.start < 0 || block.end() > coefficients.dimension(2)) {
            throw std::out_of_range(
                "Successive-orders scattering block is out of range");
        }
        const std::size_t wavelength_size = coefficient_stride * num_geometry;
        const double* wavelength_coefficients =
            coefficients.data() +
            static_cast<std::size_t>(block.start) * wavelength_size;
        auto& solver = m_rust_solvers[threadidx];
        sasktran2::rust::successive_orders::
            set_scattering_coefficients_from_atmosphere_batch4(
                *solver, *m_rust_scattering_interpolators.front(),
                {wavelength_coefficients,
                 wavelength_size * static_cast<std::size_t>(block.count)},
                coefficient_stride);

        auto boundary_values = sasktran2::rust::successive_orders::
            batch4_boundary_scattering_values_mut(*solver);
        const std::size_t boundary_size =
            m_rust_boundary_scattering_offsets.back();
        if (boundary_values.size() !=
            boundary_size * static_cast<std::size_t>(block.count)) {
            throw std::logic_error(
                "Rust SIMD boundary scattering storage size mismatch");
        }
        auto& scratch =
            m_thread_storage[threadidx].rust_boundary_scattering_scratch;
        for (int lane = 0; lane < block.count; ++lane) {
            const int wavelength = block.wavelength(lane);
            for (int ground_index = 0;
                 ground_index < m_location_interpolator->num_ground_points();
                 ++ground_index) {
                const std::size_t point_index = static_cast<std::size_t>(
                    m_location_interpolator->num_interior_points() +
                    ground_index);
                if (!m_diffuse_point_full_calculation[point_index]) {
                    continue;
                }
                const std::size_t offset =
                    m_rust_boundary_scattering_offsets[ground_index];
                const std::size_t count =
                    m_rust_boundary_scattering_offsets[ground_index + 1] -
                    offset;
                scratch.resize(count);
                const auto& point = m_diffuse_points[point_index];
                point->sphere_pair().calculate_ground_scattering_values(
                    m_atmosphere->surface(), point->location(), wavelength,
                    scratch.data());
                for (std::size_t element = 0; element < count; ++element) {
                    boundary_values[(offset + element) * 4 + lane] =
                        scratch[element];
                }
            }
        }
    }
#endif

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::iterate_to_solution(int wavelidx, int threadidx) {
        ZoneScopedN("Iterate to Solution");
#ifdef SKTRAN_RUST_SUPPORT
        iterate_to_solution_rust(wavelidx, threadidx);
#else
        (void)wavelidx;
        (void)threadidx;
        throw std::logic_error(
            "The Rust successive-orders source requires Rust support");
#endif
    }

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    void SourceImpl<NSTOKES>::generate_twostream_initial_guess(int wavelidx,
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
        thread_storage.twostream_initial_guess.setZero();
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
                        thread_storage.twostream_initial_guess(state_index) =
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
                    thread_storage.twostream_initial_guess(state_index) =
                        ground_source;
                }
            }
            sasktran2::rust::twostream::clear_output(*source);
        }
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::pack_rust_scattering_coefficient_jvp(
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
    void SourceImpl<NSTOKES>::pack_rust_boundary_scattering_jvp(
        int wavelidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        auto tangent = sasktran2::rust::successive_orders::
            boundary_scattering_value_tangent_mut(*m_rust_solvers[threadidx]);
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
    void SourceImpl<NSTOKES>::accumulate_rust_scattering_coefficient_vjp(
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
    void SourceImpl<NSTOKES>::accumulate_rust_boundary_scattering_vjp(
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
    void SourceImpl<NSTOKES>::iterate_to_solution_rust(int wavelidx,
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
                ? as_rust_slice(storage.twostream_initial_guess)
                : as_rust_slice(no_initial_guess);
        auto& solver = m_rust_solvers[threadidx];
        sasktran2::rust::successive_orders::solve_coefficients_assembled(
            *solver, initial_guess);
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
    void SourceImpl<NSTOKES>::capture_forward_state(int wavelidx,
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

        const auto forcing =
            sasktran2::rust::successive_orders::first_order_forcing(*solver);
        state.first_order_forcing.assign(forcing.begin(), forcing.end());
        const auto solution =
            sasktran2::rust::successive_orders::solution(*solver);
        state.solution.assign(solution.begin(), solution.end());
        state.valid = true;
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::pack_rust_end_of_ray_source(int wavelidx,
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
    void SourceImpl<NSTOKES>::pack_rust_end_of_ray_source_batch4(
        const sasktran2::WavelengthBlock<>& block, int threadidx) const {
        if constexpr (NSTOKES != 1) {
            throw std::logic_error(
                "Four-wavelength end-of-ray packing requires scalar mode");
        }
        if (block.count != 4) {
            throw std::invalid_argument(
                "Rust SIMD end-of-ray packing requires four wavelengths");
        }
        auto& packed = m_thread_storage[threadidx].rust_batch_end_of_ray_source;
        const int num_rays =
            static_cast<int>(m_internal_viewing.traced_rays.size());
        if (packed.size() != static_cast<std::size_t>(num_rays * block.count)) {
            throw std::logic_error("Rust SIMD end-of-ray source size mismatch");
        }
        const int num_threads = m_config->num_source_threads();
        std::vector<sasktran2::WavelengthBlockDual<NSTOKES>> scratch(
            num_threads);
        for (auto& source_value : scratch) {
            source_value.resize(block.count, 0, true);
        }
#pragma omp parallel for num_threads(num_threads)
        for (int rayidx = 0; rayidx < num_rays; ++rayidx) {
#ifdef SKTRAN_OPENMP_SUPPORT
            const int source_threadidx =
                num_threads == 1 ? 0 : omp_get_thread_num();
#else
            const int source_threadidx = 0;
#endif
            auto& source_value = scratch[source_threadidx];
            source_value.set_zero(block.count);
            for (const auto* source : m_initial_sources) {
                source->end_of_ray_source(block, rayidx, threadidx,
                                          source_threadidx + threadidx,
                                          source_value);
            }
            for (int lane = 0; lane < block.count; ++lane) {
                packed[static_cast<std::size_t>(rayidx * block.count + lane)] =
                    source_value.value(0, lane);
            }
        }
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::pack_rust_end_of_ray_source_jvp(
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
    void SourceImpl<NSTOKES>::accumulate_rust_end_of_ray_source_vjp(
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
    void SourceImpl<NSTOKES>::assemble_rust_ray_transport(
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
        const bool retain_layer_scratch =
            m_active_vjp_derivatives[threadidx].prepared;

        if (!assemble_first_order) {
            sasktran2::rust::successive_orders::assemble_solver_ray_transport(
                *m_rust_ray_transports.front(), *m_rust_solvers[threadidx],
                extinction, single_scatter_albedo, retain_layer_scratch);
            return;
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
    void SourceImpl<NSTOKES>::assemble_rust_ray_transport_batch4(
        const sasktran2::WavelengthBlock<>& block, int threadidx) {
        if constexpr (NSTOKES != 1) {
            throw std::logic_error(
                "Four-wavelength Rust transport requires scalar mode");
        }
        if (block.count != 4 || m_rust_ray_transports.empty() ||
            m_atmosphere == nullptr) {
            throw std::logic_error(
                "Four-wavelength Rust transport is not initialized");
        }
        const auto& atmosphere_storage = m_atmosphere->storage();
        const std::size_t num_geometry = static_cast<std::size_t>(
            atmosphere_storage.total_extinction.rows());
        if (block.start < 0 ||
            block.end() > atmosphere_storage.total_extinction.cols()) {
            throw std::out_of_range(
                "Successive-orders wavelength block is out of range");
        }
        const std::size_t atmosphere_batch_size =
            num_geometry * static_cast<std::size_t>(block.count);
        const ::rust::Slice<const double> extinction(
            atmosphere_storage.total_extinction.col(block.start).data(),
            atmosphere_batch_size);
        const ::rust::Slice<const double> single_scatter_albedo(
            atmosphere_storage.ssa.col(block.start).data(),
            atmosphere_batch_size);

        const double* solar_data = nullptr;
        std::size_t solar_size = 0;
        if (m_geometry_1d != nullptr) {
            const auto* single_scatter =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionTable,
                    NSTOKES>*>(m_initial_sources.front());
            if (single_scatter != nullptr) {
                solar_data =
                    single_scatter->solar_transmission_batch_data(threadidx);
                solar_size =
                    single_scatter->solar_transmission_batch_size(threadidx);
            }
        } else {
            const auto* single_scatter =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionExact,
                    NSTOKES>*>(m_initial_sources.front());
            if (single_scatter != nullptr) {
                solar_data =
                    single_scatter->solar_transmission_batch_data(threadidx);
                solar_size =
                    single_scatter->solar_transmission_batch_size(threadidx);
            }
        }
        if (solar_data == nullptr || m_initial_sources.size() != 1) {
            throw std::logic_error(
                "Rust SIMD forcing requires one batched single-scatter "
                "source");
        }

        const std::size_t num_phase_moments =
            static_cast<std::size_t>(atmosphere_storage.max_stored_legendre());
        if (num_phase_moments !=
            static_cast<std::size_t>(m_config->num_singlescatter_moments())) {
            throw std::logic_error("Rust SIMD phase storage size mismatch");
        }
        const ::rust::Slice<const double> legendre_coefficients(
            &atmosphere_storage.leg_coeff(0, 0, block.start),
            atmosphere_batch_size * num_phase_moments);
        static_assert(sizeof(int) == sizeof(std::int32_t));
        const ::rust::Slice<const std::int32_t> maximum_order(
            reinterpret_cast<const std::int32_t*>(
                atmosphere_storage.max_order.col(block.start).data()),
            atmosphere_batch_size);
        pack_rust_end_of_ray_source_batch4(block, threadidx);
        auto& storage = m_thread_storage[threadidx];

        sasktran2::rust::successive_orders::
            assemble_solver_ray_transport_with_first_order_batch4(
                *m_rust_ray_transports.front(), *m_rust_solvers[threadidx],
                extinction, single_scatter_albedo, legendre_coefficients,
                maximum_order, {solar_data, solar_size},
                as_rust_slice(storage.rust_batch_end_of_ray_source));
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::prepare_rust_ray_transport_vjp_attenuation(
        int wavelidx, int threadidx) {
        if (m_rust_ray_transports.empty() || m_atmosphere == nullptr) {
            throw std::logic_error(
                "Rust transport geometry is not initialized");
        }
        const auto& atmosphere = m_atmosphere->storage();
        const std::size_t num_geometry =
            static_cast<std::size_t>(atmosphere.total_extinction.rows());
        sasktran2::rust::successive_orders::
            prepare_solver_ray_transport_vjp_attenuation(
                *m_rust_ray_transports.front(), *m_rust_solvers[threadidx],
                {atmosphere.total_extinction.col(wavelidx).data(),
                 num_geometry},
                {atmosphere.ssa.col(wavelidx).data(), num_geometry});
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::assemble_rust_ray_transport_jvp(
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
    void SourceImpl<NSTOKES>::stage_rust_ray_transport_jvp_batch4_lane(
        int wavelidx, int threadidx, int lane, bool zero_tangent,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        if constexpr (NSTOKES != 1) {
            throw std::logic_error(
                "Four-wavelength transport JVP requires scalar mode");
        }
        if (lane < 0 || lane >= 4 || m_rust_ray_transports.empty() ||
            m_atmosphere == nullptr) {
            throw std::invalid_argument(
                "Invalid four-wavelength transport JVP lane");
        }
        auto& output = m_thread_storage[threadidx];
        const auto& atmosphere = m_atmosphere->storage();
        const std::size_t num_geometry =
            static_cast<std::size_t>(atmosphere.total_extinction.rows());
        const std::size_t num_phase_moments =
            static_cast<std::size_t>(atmosphere.max_stored_legendre());
        const std::size_t num_phase_values = num_geometry * num_phase_moments;

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

        pack_rust_end_of_ray_source_jvp(wavelidx, threadidx, native_tangent);
        const auto& state =
            m_forward_states[static_cast<std::size_t>(wavelidx)];
        sasktran2::rust::successive_orders::
            stage_solver_ray_transport_jvp_batch4_lane(
                *m_rust_solvers[threadidx], static_cast<std::size_t>(lane),
                zero_tangent, {native_tangent.data(), num_geometry},
                {native_tangent.data() + m_atmosphere->ssa_deriv_start_index(),
                 num_geometry},
                as_rust_slice(legendre_tangent),
                as_rust_slice(output.rust_end_of_ray_source),
                as_rust_slice(output.rust_end_of_ray_source_tangent),
                as_rust_slice(state.first_order_forcing));
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::assemble_rust_ray_transport_jvp_batch4(
        const sasktran2::WavelengthBlock<>& block, int threadidx) {
        if constexpr (NSTOKES != 1) {
            throw std::logic_error(
                "Four-wavelength transport JVP requires scalar mode");
        }
        if (block.count != 4 || m_rust_ray_transports.empty() ||
            m_atmosphere == nullptr) {
            throw std::invalid_argument(
                "Invalid four-wavelength transport JVP block");
        }
        const auto& atmosphere = m_atmosphere->storage();
        const std::size_t num_geometry =
            static_cast<std::size_t>(atmosphere.total_extinction.rows());
        const std::size_t atmosphere_batch_size = num_geometry * 4;
        const std::size_t num_phase_moments =
            static_cast<std::size_t>(atmosphere.max_stored_legendre());
        const double* solar_data = nullptr;
        std::size_t solar_size = 0;
        const double* solar_tangent_data = nullptr;
        std::size_t solar_tangent_size = 0;
        if (m_geometry_1d != nullptr) {
            const auto* single_scatter =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionTable,
                    NSTOKES>*>(m_initial_sources.front());
            if (single_scatter != nullptr) {
                solar_data =
                    single_scatter->solar_transmission_batch_data(threadidx);
                solar_size =
                    single_scatter->solar_transmission_batch_size(threadidx);
                solar_tangent_data =
                    single_scatter->solar_transmission_product_batch_data(
                        threadidx);
                solar_tangent_size =
                    single_scatter->solar_transmission_product_batch_size(
                        threadidx);
            }
        } else {
            const auto* single_scatter =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionExact,
                    NSTOKES>*>(m_initial_sources.front());
            if (single_scatter != nullptr) {
                solar_data =
                    single_scatter->solar_transmission_batch_data(threadidx);
                solar_size =
                    single_scatter->solar_transmission_batch_size(threadidx);
                solar_tangent_data =
                    single_scatter->solar_transmission_product_batch_data(
                        threadidx);
                solar_tangent_size =
                    single_scatter->solar_transmission_product_batch_size(
                        threadidx);
            }
        }
        if (solar_data == nullptr || solar_tangent_data == nullptr ||
            solar_tangent_size != solar_size) {
            throw std::logic_error(
                "Rust SIMD JVP requires a batched single-scatter source");
        }
        sasktran2::rust::successive_orders::
            assemble_solver_ray_transport_with_first_order_jvp_batch4(
                *m_rust_ray_transports.front(), *m_rust_solvers[threadidx],
                {atmosphere.total_extinction.col(block.start).data(),
                 atmosphere_batch_size},
                {atmosphere.ssa.col(block.start).data(), atmosphere_batch_size},
                {&atmosphere.leg_coeff(0, 0, block.start),
                 atmosphere_batch_size * num_phase_moments},
                {reinterpret_cast<const std::int32_t*>(
                     atmosphere.max_order.col(block.start).data()),
                 atmosphere_batch_size},
                {solar_data, solar_size},
                {solar_tangent_data, solar_tangent_size});
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::assemble_rust_ray_transport_vjp_batch4(
        const sasktran2::WavelengthBlock<>& block, int threadidx) const {
        if constexpr (NSTOKES != 1) {
            throw std::logic_error(
                "Four-wavelength Rust transport VJP requires scalar mode");
        }
        if (block.count != 4 || m_rust_ray_transports.empty() ||
            m_atmosphere == nullptr) {
            throw std::invalid_argument(
                "Invalid four-wavelength transport VJP block");
        }
        const auto& atmosphere = m_atmosphere->storage();
        const std::size_t num_geometry =
            static_cast<std::size_t>(atmosphere.total_extinction.rows());
        const std::size_t atmosphere_batch_size = num_geometry * 4;
        const std::size_t num_phase_moments =
            static_cast<std::size_t>(atmosphere.max_stored_legendre());

        const double* solar_data = nullptr;
        std::size_t solar_size = 0;
        double* solar_gradient_data = nullptr;
        std::size_t solar_gradient_size = 0;
        if (m_geometry_1d != nullptr) {
            auto* single_scatter =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionTable,
                    NSTOKES>*>(m_initial_sources.front());
            if (single_scatter != nullptr) {
                solar_data =
                    single_scatter->solar_transmission_batch_data(threadidx);
                solar_size =
                    single_scatter->solar_transmission_batch_size(threadidx);
                solar_gradient_data =
                    single_scatter->solar_transmission_product_batch_data(
                        threadidx);
                solar_gradient_size =
                    single_scatter->solar_transmission_product_batch_size(
                        threadidx);
            }
        } else {
            auto* single_scatter =
                dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                    sasktran2::solartransmission::SolarTransmissionExact,
                    NSTOKES>*>(m_initial_sources.front());
            if (single_scatter != nullptr) {
                solar_data =
                    single_scatter->solar_transmission_batch_data(threadidx);
                solar_size =
                    single_scatter->solar_transmission_batch_size(threadidx);
                solar_gradient_data =
                    single_scatter->solar_transmission_product_batch_data(
                        threadidx);
                solar_gradient_size =
                    single_scatter->solar_transmission_product_batch_size(
                        threadidx);
            }
        }
        if (solar_data == nullptr || solar_gradient_data == nullptr ||
            solar_gradient_size != solar_size ||
            m_initial_sources.size() != 1) {
            throw std::logic_error(
                "Rust SIMD VJP requires one batched single-scatter source");
        }

        pack_rust_end_of_ray_source_batch4(block, threadidx);
        const auto& output = m_thread_storage[threadidx];
        sasktran2::rust::successive_orders::
            assemble_solver_ray_transport_with_first_order_vjp_batch4(
                *m_rust_ray_transports.front(), *m_rust_solvers[threadidx],
                {atmosphere.total_extinction.col(block.start).data(),
                 atmosphere_batch_size},
                {atmosphere.ssa.col(block.start).data(), atmosphere_batch_size},
                {&atmosphere.leg_coeff(0, 0, block.start),
                 atmosphere_batch_size * num_phase_moments},
                {reinterpret_cast<const std::int32_t*>(
                     atmosphere.max_order.col(block.start).data()),
                 atmosphere_batch_size},
                {solar_data, solar_size},
                as_rust_slice(output.rust_batch_end_of_ray_source),
                {solar_gradient_data, solar_gradient_size});
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::select_rust_ray_transport_vjp_batch4_lane(
        int threadidx, int lane) const {
        if constexpr (NSTOKES != 1) {
            throw std::logic_error(
                "Four-wavelength Rust transport VJP requires scalar mode");
        }
        auto& output = m_thread_storage[threadidx];
        const auto& atmosphere = m_atmosphere->storage();
        const std::size_t num_geometry =
            static_cast<std::size_t>(atmosphere.total_extinction.rows());
        const std::size_t num_phase_moments =
            static_cast<std::size_t>(atmosphere.max_stored_legendre());
        output.rust_single_scatter_extinction_gradient.resize(num_geometry);
        output.rust_single_scatter_albedo_gradient.resize(num_geometry);
        output.rust_single_scatter_legendre_gradient.resize(num_geometry *
                                                            num_phase_moments);
        sasktran2::rust::successive_orders::
            select_solver_ray_transport_vjp_batch4_lane(
                *m_rust_solvers[threadidx], static_cast<std::size_t>(lane),
                as_rust_mut_slice(
                    output.rust_single_scatter_extinction_gradient),
                as_rust_mut_slice(output.rust_single_scatter_albedo_gradient),
                as_rust_mut_slice(output.rust_single_scatter_legendre_gradient),
                as_rust_mut_slice(output.rust_end_of_ray_source_gradient));
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::accumulate_rust_transport_vjp_with_first_order(
        int wavelidx, int threadidx,
        Eigen::Ref<Eigen::VectorXd> native_gradient,
        bool transport_vjp_precomputed, bool packed_solar_vjp) const {
        const auto& output = m_thread_storage[threadidx];
        const auto& activity = m_active_vjp_derivatives[threadidx];
        if (!transport_vjp_precomputed) {
            pack_rust_end_of_ray_source(wavelidx, threadidx);
        }

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

        if (!transport_vjp_precomputed) {
            output.rust_single_scatter_extinction_gradient.resize(num_geometry);
            output.rust_single_scatter_albedo_gradient.resize(num_geometry);
            output.rust_single_scatter_legendre_gradient.resize(
                num_phase_values);
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
                    as_rust_mut_slice(
                        output.rust_single_scatter_albedo_gradient),
                    as_rust_mut_slice(
                        output.rust_single_scatter_legendre_gradient),
                    as_rust_mut_slice(
                        output.rust_single_scatter_solar_gradient),
                    as_rust_mut_slice(output.rust_end_of_ray_source_gradient));
        }

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

        if (!packed_solar_vjp) {
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
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::restore_transport_operator(int wavelidx,
                                                         int threadidx) {
        if (m_rust_ray_transports.empty()) {
            throw std::logic_error(
                "Rust successive-orders transport is not initialized");
        }
        assemble_rust_ray_transport(wavelidx, threadidx);
    }

#endif

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::prepare_wavelength(int wavelidx, int threadidx,
                                                 bool value_only) {
        const sasktran2::WavelengthBlock<> scalar_block{wavelidx, 1};
        for (auto& source : m_initial_owned_sources) {
            if (value_only) {
                source->calculate_value(scalar_block, threadidx);
            } else {
                source->calculate(scalar_block, threadidx);
            }
        }
        if (m_config->num_hr_spherical_iterations() > 0) {
#ifdef SKTRAN_RUST_SUPPORT
            if (m_rust_ray_transports.empty()) {
                throw std::logic_error(
                    "Rust successive-orders transport is not initialized");
            }
            assemble_rust_ray_transport(wavelidx, threadidx, true);
#else
            throw std::logic_error(
                "The Rust successive-orders source requires Rust support");
#endif

            // Generate the scattering coefficients and boundary terms.
            generate_scattering_matrices(wavelidx, threadidx);
        }
    }

    template <int NSTOKES>
    void
    SourceImpl<NSTOKES>::calculate(const sasktran2::WavelengthBlock<>& block,
                                   int threadidx) {
        m_deferred_jvp_transport_restore[threadidx] = false;
        if (block.count < 1 || block.count > maximum_wavelength_block_size() ||
            block.count > m_wavelength_block_capacity) {
            throw std::invalid_argument(
                "Invalid successive-orders wavelength block");
        }

#ifdef SKTRAN_RUST_SUPPORT
        {
            auto& storage = m_thread_storage[threadidx];
            if (m_wavelength_block_capacity == 1) {
                prepare_wavelength(block.start, threadidx);
                iterate_to_solution(block.start, threadidx);
                capture_forward_state(block.start, threadidx);
                return;
            }
            if constexpr (NSTOKES == 1) {
                if (block.count == 4) {
                    auto& solver = m_rust_solvers[threadidx];
                    const bool use_twostream_initial_guess =
                        m_config->successive_orders_initialization() ==
                        sasktran2::Config::SuccessiveOrdersInitialization::
                            twostream;
                    if (use_twostream_initial_guess) {
                        sasktran2::rust::successive_orders::
                            begin_coefficients_batch4(*solver);
                        for (int lane = 0; lane < block.count; ++lane) {
                            const int wavelength = block.wavelength(lane);
                            prepare_wavelength(wavelength, threadidx);
                            generate_twostream_initial_guess(wavelength,
                                                             threadidx);
                            sasktran2::rust::successive_orders::
                                stage_coefficients_batch4_lane(
                                    *solver, static_cast<std::size_t>(lane),
                                    as_rust_slice(
                                        storage.twostream_initial_guess));
                        }
                    } else {
                        try {
                            for (auto& source : m_initial_owned_sources) {
                                source->calculate(block, threadidx);
                            }
                        } catch (const std::exception& error) {
                            throw std::runtime_error(
                                std::string("SIMD single-scatter batch: ") +
                                error.what());
                        }
                        try {
                            assemble_rust_ray_transport_batch4(block,
                                                               threadidx);
                        } catch (const std::exception& error) {
                            throw std::runtime_error(
                                std::string("SIMD ray assembly: ") +
                                error.what());
                        }
                        try {
                            generate_scattering_matrices_batch4(block,
                                                                threadidx);
                        } catch (const std::exception& error) {
                            throw std::runtime_error(
                                std::string("SIMD scattering assembly: ") +
                                error.what());
                        }
                    }
                    try {
                        sasktran2::rust::successive_orders::
                            solve_coefficients_batch4(*solver);
                    } catch (const std::exception& error) {
                        throw std::runtime_error(
                            std::string("SIMD successive-orders solve: ") +
                            error.what());
                    }
                    spdlog::debug(
                        "Rust successive-orders four-wavelength SIMD "
                        "workspace: {:.3f} MiB on wavelength worker {}",
                        static_cast<double>(
                            sasktran2::rust::successive_orders::
                                batch4_workspace_bytes(*solver)) /
                            (1024.0 * 1024.0),
                        threadidx);
                    const auto packed_solution =
                        sasktran2::rust::successive_orders::batch4_solution(
                            *solver);
                    if (packed_solution.size() != m_num_outgoing_values * 4) {
                        throw std::logic_error(
                            "Rust SIMD successive-orders solution size "
                            "mismatch");
                    }
                    std::copy(packed_solution.begin(), packed_solution.end(),
                              storage.rust_batch_outgoing_sources.begin());

                    const auto packed_forcing = sasktran2::rust::
                        successive_orders::batch4_first_order_forcing(*solver);
                    const bool tolerance_requested =
                        m_config->successive_orders_relative_tolerance() > 0 ||
                        m_config->successive_orders_absolute_tolerance() > 0;
                    for (int lane = 0; lane < block.count; ++lane) {
                        const int wavelength = block.wavelength(lane);
                        const bool converged = sasktran2::rust::
                            successive_orders::batch4_converged(*solver, lane);
                        if (tolerance_requested && !converged) {
                            spdlog::warn(
                                "Rust SIMD successive-orders source reached "
                                "{} iterations with residual {} for "
                                "wavelength index {}",
                                sasktran2::rust::successive_orders::
                                    batch4_iterations(*solver, lane),
                                sasktran2::rust::successive_orders::
                                    batch4_final_residual(*solver, lane),
                                wavelength);
                        } else {
                            spdlog::debug(
                                "Rust SIMD successive-orders source "
                                "completed {} iterations; initial residual "
                                "{}, final residual {} for wavelength index "
                                "{}",
                                sasktran2::rust::successive_orders::
                                    batch4_iterations(*solver, lane),
                                sasktran2::rust::successive_orders::
                                    batch4_initial_residual(*solver, lane),
                                sasktran2::rust::successive_orders::
                                    batch4_final_residual(*solver, lane),
                                wavelength);
                        }

                        if (!m_capture_forward_state || wavelength < 0 ||
                            static_cast<std::size_t>(wavelength) >=
                                m_forward_states.size()) {
                            continue;
                        }
                        auto& state = m_forward_states[static_cast<std::size_t>(
                            wavelength)];
                        state.valid = converged;
                        if (!converged) {
                            continue;
                        }
                        state.solution.resize(m_num_outgoing_values);
                        for (std::size_t element = 0;
                             element < m_num_outgoing_values; ++element) {
                            state.solution[element] =
                                packed_solution[element * 4 + lane];
                        }
                        if (packed_forcing.size() % 4 != 0) {
                            throw std::logic_error(
                                "Rust SIMD successive-orders forcing size "
                                "mismatch");
                        }
                        state.first_order_forcing.resize(packed_forcing.size() /
                                                         4);
                        for (std::size_t element = 0;
                             element < state.first_order_forcing.size();
                             ++element) {
                            state.first_order_forcing[element] =
                                packed_forcing[element * 4 + lane];
                        }
                    }
                    return;
                }
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
    void SourceImpl<NSTOKES>::calculate_value(
        const sasktran2::WavelengthBlock<>& block, int threadidx) {
        m_deferred_jvp_transport_restore[threadidx] = false;
        if (block.count != 1) {
            calculate(block, threadidx);
            return;
        }
        prepare_wavelength(block.start, threadidx, true);
        iterate_to_solution(block.start, threadidx);
#ifdef SKTRAN_RUST_SUPPORT
        capture_forward_state(block.start, threadidx);
#endif
    }

    template <int NSTOKES>
    bool SourceImpl<NSTOKES>::restore_forward_state(int wavelidx,
                                                    int threadidx) {
#ifdef SKTRAN_RUST_SUPPORT
        m_deferred_jvp_transport_restore[threadidx] = false;
        if (m_forward_state_atmosphere != m_atmosphere ||
            m_forward_state_atmosphere_revision != m_atmosphere->revision() ||
            wavelidx < 0 ||
            static_cast<std::size_t>(wavelidx) >= m_forward_states.size() ||
            !m_forward_states[static_cast<std::size_t>(wavelidx)].valid) {
            return false;
        }

        const auto& state =
            m_forward_states[static_cast<std::size_t>(wavelidx)];
        const std::size_t expected_forcing_size =
            m_internal_viewing.traced_rays.size() * NSTOKES;
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

        sasktran2::rust::successive_orders::
            restore_coefficient_solution_assembled(
                *m_rust_solvers[threadidx],
                as_rust_slice(state.first_order_forcing),
                as_rust_slice(state.solution));
        return true;
#else
        (void)wavelidx;
        (void)threadidx;
        return false;
#endif
    }

    template <int NSTOKES>
    bool SourceImpl<NSTOKES>::restore_forward_state_block(
        const sasktran2::WavelengthBlock<>& block, int threadidx) {
#ifdef SKTRAN_RUST_SUPPORT
        m_deferred_jvp_transport_restore[threadidx] = false;
        if constexpr (NSTOKES != 1) {
            return false;
        }
        if (m_rust_ray_transports.empty() || block.count != 4 ||
            m_wavelength_block_capacity != 4 ||
            m_forward_state_atmosphere != m_atmosphere ||
            m_forward_state_atmosphere_revision != m_atmosphere->revision() ||
            block.start < 0 ||
            static_cast<std::size_t>(block.end()) > m_forward_states.size()) {
            return false;
        }
        const std::size_t expected_forcing_size =
            m_internal_viewing.traced_rays.size();
        for (int lane = 0; lane < block.count; ++lane) {
            const auto& state = m_forward_states[static_cast<std::size_t>(
                block.wavelength(lane))];
            if (!state.valid ||
                state.first_order_forcing.size() != expected_forcing_size ||
                state.solution.size() != m_num_outgoing_values) {
                return false;
            }
        }

        // The converged solution and scattering blocks are enough to begin a
        // derivative calculation. JVP reconstructs transport together with
        // its tangent lane-by-lane, while VJP defers the packed transport
        // assembly until its block hook. Avoiding it here prevents JVP from
        // assembling the primal transport twice.
        generate_scattering_matrices_batch4(block, threadidx);

        auto& packed_solution =
            m_thread_storage[threadidx].rust_batch_outgoing_sources;
        packed_solution.resize(m_num_outgoing_values * 4);
        for (int lane = 0; lane < block.count; ++lane) {
            const auto& solution = m_forward_states[static_cast<std::size_t>(
                                                        block.wavelength(lane))]
                                       .solution;
            for (std::size_t element = 0; element < m_num_outgoing_values;
                 ++element) {
                packed_solution[element * 4 + lane] = solution[element];
            }
        }
        sasktran2::rust::successive_orders::
            restore_coefficients_batch4_solution(
                *m_rust_solvers[threadidx], as_rust_slice(packed_solution));
        m_deferred_jvp_transport_restore[threadidx] = true;
        return true;
#else
        (void)block;
        (void)threadidx;
        return false;
#endif
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::prepare_forward_state_for_jvp(
        int wavelidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        for (auto* source : m_initial_sources) {
            source->prepare_forward_state_for_jvp(wavelidx, threadidx,
                                                  native_tangent);
        }
    }

    template <int NSTOKES>
    bool SourceImpl<NSTOKES>::restore_forward_state_for_jvp(
        int wavelidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
#ifdef SKTRAN_RUST_SUPPORT
        m_deferred_jvp_transport_restore[threadidx] = false;
        if (m_forward_state_atmosphere != m_atmosphere ||
            m_forward_state_atmosphere_revision != m_atmosphere->revision() ||
            wavelidx < 0 ||
            static_cast<std::size_t>(wavelidx) >= m_forward_states.size() ||
            !m_forward_states[static_cast<std::size_t>(wavelidx)].valid) {
            return false;
        }

        const auto& state =
            m_forward_states[static_cast<std::size_t>(wavelidx)];
        const std::size_t expected_forcing_size =
            m_internal_viewing.traced_rays.size() * NSTOKES;
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
    void SourceImpl<NSTOKES>::prepare_jvp(
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
        if (m_rust_ray_transports.empty()) {
            throw std::logic_error(
                "Rust successive-orders transport is not initialized");
        }
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

        assemble_rust_ray_transport_jvp(wavelidx, wavel_threadidx,
                                        native_tangent);
        auto& solver = m_rust_solvers[wavel_threadidx];
        if (restore_transport) {
            const auto& state =
                m_forward_states[static_cast<std::size_t>(wavelidx)];
            sasktran2::rust::successive_orders::
                restore_coefficient_solution_assembled(
                    *solver, as_rust_slice(state.first_order_forcing),
                    as_rust_slice(state.solution));
        }
        sasktran2::rust::successive_orders::
            linearize_coefficients_jvp_assembled(*solver);
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
    bool SourceImpl<NSTOKES>::prepare_jvp_block(
        const sasktran2::WavelengthBlock<>& block, int wavel_threadidx,
        const std::vector<Eigen::VectorXd>& native_tangents) {
#ifdef SKTRAN_RUST_SUPPORT
        if constexpr (NSTOKES != 1) {
            return false;
        }
        if (m_rust_ray_transports.empty() || block.count != 4 ||
            m_wavelength_block_capacity != 4 || native_tangents.size() != 4) {
            return false;
        }
        bool all_zero = true;
        for (int lane = 0; lane < block.count; ++lane) {
            if (native_tangents[lane].size() != m_atmosphere->num_deriv()) {
                throw std::invalid_argument(
                    "Four-wavelength JVP tangent has an invalid size");
            }
            all_zero = all_zero && native_tangents[lane].isZero(0.0);
            if (!sasktran2::rust::successive_orders::batch4_converged(
                    *m_rust_solvers[wavel_threadidx],
                    static_cast<std::size_t>(lane))) {
                return false;
            }
        }
        if (m_config->successive_orders_initialization() ==
                sasktran2::Config::SuccessiveOrdersInitialization::twostream &&
            m_config->successive_orders_relative_tolerance() == 0 &&
            m_config->successive_orders_absolute_tolerance() == 0 &&
            m_config->successive_orders_anderson_depth() == 0) {
            return false;
        }

        m_batch_jvp_zero_tangent[wavel_threadidx] = all_zero;
        m_deferred_jvp_transport_restore[wavel_threadidx] = false;
        if (all_zero) {
            return true;
        }

        // Forward-state retention does not keep the single-scatter source's
        // block workspace. Rebuild it once with the packed solar operator,
        // then select its four lanes below; this replaces four scalar primal
        // solar-path evaluations while preserving the correct block identity.
        for (auto& source : m_initial_owned_sources) {
            source->calculate(block, wavel_threadidx);
        }
        auto* table_source =
            dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                sasktran2::solartransmission::SolarTransmissionTable,
                NSTOKES>*>(m_initial_sources.front());
        auto* exact_source =
            dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                sasktran2::solartransmission::SolarTransmissionExact,
                NSTOKES>*>(m_initial_sources.front());
        if (table_source != nullptr) {
            table_source->calculate_solar_transmission_jvp_batch4(
                wavel_threadidx, native_tangents);
        } else if (exact_source != nullptr) {
            exact_source->calculate_solar_transmission_jvp_batch4(
                wavel_threadidx, native_tangents);
        } else {
            throw std::logic_error(
                "Rust SIMD JVP requires a batched single-scatter source");
        }

        auto& solver = m_rust_solvers[wavel_threadidx];
        sasktran2::rust::successive_orders::
            begin_linearize_coefficients_jvp_batch4(*solver);
        for (int lane = 0; lane < block.count; ++lane) {
            const int wavelength = block.wavelength(lane);
            const auto& tangent = native_tangents[lane];
            const sasktran2::WavelengthBlock<> scalar_block{wavelength, 1};
            for (auto& source : m_initial_owned_sources) {
                auto* block_source = source->wavelength_block_linearization();
                if (block_source == nullptr ||
                    !block_source->select_calculated_block_lane(
                        lane, wavel_threadidx)) {
                    source->prepare_forward_state_for_jvp(
                        wavelength, wavel_threadidx, tangent);
                    if (tangent.isZero(0.0)) {
                        source->calculate_value(scalar_block, wavel_threadidx);
                    } else {
                        source->calculate(scalar_block, wavel_threadidx);
                    }
                }
            }
            const bool zero_tangent = tangent.isZero(0.0);
            if (!zero_tangent) {
                pack_rust_scattering_coefficient_jvp(wavelength,
                                                     wavel_threadidx, tangent);
                pack_rust_boundary_scattering_jvp(wavelength, wavel_threadidx,
                                                  tangent);
            }
            for (auto* source : m_initial_sources) {
                if (source != m_initial_sources.front()) {
                    source->prepare_jvp(wavelength, wavel_threadidx, tangent);
                }
            }
            if (table_source != nullptr) {
                table_source->select_solar_transmission_jvp_batch4_lane(
                    lane, wavel_threadidx, tangent);
            } else {
                exact_source->select_solar_transmission_jvp_batch4_lane(
                    lane, wavel_threadidx, tangent);
            }
            stage_rust_ray_transport_jvp_batch4_lane(
                wavelength, wavel_threadidx, lane, zero_tangent, tangent);
        }
        assemble_rust_ray_transport_jvp_batch4(block, wavel_threadidx);
        sasktran2::rust::successive_orders::linearize_coefficients_jvp_batch4(
            *solver);
        return true;
#else
        (void)block;
        (void)wavel_threadidx;
        (void)native_tangents;
        return false;
#endif
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::select_jvp_block_lane(int lane,
                                                    int wavel_threadidx) {
#ifdef SKTRAN_RUST_SUPPORT
        if (lane < 0 || lane >= 4) {
            throw std::out_of_range("Invalid successive-orders JVP block lane");
        }
        if (m_batch_jvp_zero_tangent[wavel_threadidx]) {
            sasktran2::rust::successive_orders::select_coefficients_batch4_lane(
                *m_rust_solvers[wavel_threadidx],
                static_cast<std::size_t>(lane));
            sasktran2::rust::successive_orders::clear_solution_jvp(
                *m_rust_solvers[wavel_threadidx]);
        } else {
            sasktran2::rust::successive_orders::
                select_linearize_coefficients_jvp_batch4_lane(
                    *m_rust_solvers[wavel_threadidx],
                    static_cast<std::size_t>(lane));
        }
#else
        (void)lane;
        (void)wavel_threadidx;
#endif
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::prepare_vjp(
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
        activity = DerivativeActivity{};
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
    bool SourceImpl<NSTOKES>::begin_vjp_block(
        const sasktran2::WavelengthBlock<>& block, int wavel_threadidx,
        Eigen::Ref<const Eigen::VectorXi> active_derivatives) {
#ifdef SKTRAN_RUST_SUPPORT
        if constexpr (NSTOKES != 1) {
            return false;
        }
        if (m_rust_ray_transports.empty() || block.count != 4 ||
            m_wavelength_block_capacity != 4) {
            return false;
        }
        for (int lane = 0; lane < block.count; ++lane) {
            if (!sasktran2::rust::successive_orders::batch4_converged(
                    *m_rust_solvers[wavel_threadidx],
                    static_cast<std::size_t>(lane))) {
                return false;
            }
        }
        if (m_deferred_jvp_transport_restore[wavel_threadidx]) {
            for (auto& source : m_initial_owned_sources) {
                source->calculate(block, wavel_threadidx);
            }
            assemble_rust_ray_transport_batch4(block, wavel_threadidx);
            sasktran2::rust::successive_orders::
                restore_coefficients_batch4_solution(
                    *m_rust_solvers[wavel_threadidx],
                    as_rust_slice(m_thread_storage[wavel_threadidx]
                                      .rust_batch_outgoing_sources));
            m_deferred_jvp_transport_restore[wavel_threadidx] = false;
        }
        prepare_vjp(block.start, wavel_threadidx, active_derivatives);
        m_batch_vjp_derivative_activity[wavel_threadidx] =
            m_active_vjp_derivatives[wavel_threadidx];
        m_batch_vjp_active_derivatives[wavel_threadidx] = active_derivatives;
        sasktran2::rust::successive_orders::
            begin_linearize_coefficients_vjp_batch4(
                *m_rust_solvers[wavel_threadidx]);
        return true;
#else
        (void)block;
        (void)wavel_threadidx;
        (void)active_derivatives;
        return false;
#endif
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::select_vjp_block_lane(int lane,
                                                    int wavel_threadidx) {
#ifdef SKTRAN_RUST_SUPPORT
        if (lane < 0 || lane >= 4) {
            throw std::out_of_range("Invalid successive-orders VJP block lane");
        }
        for (auto& cotangent : m_los_solution_cotangents[wavel_threadidx]) {
            cotangent.setZero();
        }
        sasktran2::rust::successive_orders::select_coefficients_batch4_lane(
            *m_rust_solvers[wavel_threadidx], static_cast<std::size_t>(lane));
#else
        (void)lane;
        (void)wavel_threadidx;
#endif
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::stage_vjp_block_lane(int lane,
                                                   int wavel_threadidx) {
#ifdef SKTRAN_RUST_SUPPORT
        if (lane < 0 || lane >= 4) {
            throw std::out_of_range(
                "Invalid successive-orders staged VJP block lane");
        }
        Eigen::VectorXd solution_cotangent = Eigen::VectorXd::Zero(
            static_cast<Eigen::Index>(m_num_outgoing_values));
        for (const auto& thread_cotangent :
             m_los_solution_cotangents[wavel_threadidx]) {
            solution_cotangent += thread_cotangent;
        }
        sasktran2::rust::successive_orders::
            stage_linearize_coefficients_vjp_batch4_lane(
                *m_rust_solvers[wavel_threadidx],
                static_cast<std::size_t>(lane),
                as_rust_slice(solution_cotangent));
#else
        (void)lane;
        (void)wavel_threadidx;
#endif
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::finalize_vjp_block(
        const sasktran2::WavelengthBlock<>& block, int wavel_threadidx,
        std::vector<Eigen::VectorXd>& native_gradients) const {
#ifdef SKTRAN_RUST_SUPPORT
        if (block.count != 4 || native_gradients.size() != 4) {
            throw std::invalid_argument(
                "Invalid four-wavelength VJP finalization block");
        }
        auto& solver = m_rust_solvers[wavel_threadidx];
        sasktran2::rust::successive_orders::linearize_coefficients_vjp_batch4(
            *solver);
        assemble_rust_ray_transport_vjp_batch4(block, wavel_threadidx);
        spdlog::debug(
            "Rust successive-orders four-wavelength SIMD product workspace: "
            "{:.3f} MiB on wavelength worker {}",
            static_cast<double>(
                sasktran2::rust::successive_orders::batch4_workspace_bytes(
                    *solver)) /
                (1024.0 * 1024.0),
            wavel_threadidx);
        auto* table_source =
            dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                sasktran2::solartransmission::SolarTransmissionTable,
                NSTOKES>*>(m_initial_sources.front());
        auto* exact_source =
            dynamic_cast<sasktran2::solartransmission::SingleScatterSource<
                sasktran2::solartransmission::SolarTransmissionExact,
                NSTOKES>*>(m_initial_sources.front());
        if (table_source == nullptr && exact_source == nullptr) {
            throw std::logic_error(
                "Rust SIMD VJP requires a batched single-scatter source");
        }
        for (int lane = 0; lane < block.count; ++lane) {
            const int wavelength = block.wavelength(lane);
            m_active_vjp_derivatives[wavel_threadidx] =
                m_batch_vjp_derivative_activity[wavel_threadidx];
            const auto& active_derivatives =
                m_batch_vjp_active_derivatives[wavel_threadidx];
            // begin_vjp_block already prepared lane zero before LOS
            // cotangents were staged. Reuse that empty source-wide reverse
            // state; later lanes must reset it after the preceding finalize.
            if (lane != 0) {
                for (auto* source : m_initial_sources) {
                    source->prepare_vjp(wavelength, wavel_threadidx,
                                        active_derivatives);
                }
            }
            const sasktran2::WavelengthBlock<> scalar_block{wavelength, 1};
            for (auto& source : m_initial_owned_sources) {
                auto* block_source = source->wavelength_block_linearization();
                if (block_source == nullptr ||
                    !block_source->select_calculated_block_lane(
                        lane, wavel_threadidx)) {
                    source->calculate(scalar_block, wavel_threadidx);
                }
            }
            sasktran2::rust::successive_orders::
                select_linearize_coefficients_vjp_batch4_lane(
                    *solver, static_cast<std::size_t>(lane));
            if (table_source != nullptr) {
                table_source->select_solar_transmission_vjp_batch4_lane(
                    lane, wavel_threadidx);
            } else {
                exact_source->select_solar_transmission_vjp_batch4_lane(
                    lane, wavel_threadidx);
            }
            select_rust_ray_transport_vjp_batch4_lane(wavel_threadidx, lane);
            accumulate_rust_vjp_products(wavelength, wavel_threadidx,
                                         native_gradients[lane], false, true,
                                         true);
            if (table_source != nullptr) {
                table_source->stage_solar_transmission_vjp_batch4_lane(
                    lane, wavel_threadidx);
            } else {
                exact_source->stage_solar_transmission_vjp_batch4_lane(
                    lane, wavel_threadidx);
            }
        }
        if (table_source != nullptr) {
            table_source->finalize_solar_transmission_vjp_batch4(
                wavel_threadidx, native_gradients);
        } else {
            exact_source->finalize_solar_transmission_vjp_batch4(
                wavel_threadidx, native_gradients);
        }
        m_active_vjp_derivatives[wavel_threadidx] = DerivativeActivity{};
#else
        (void)block;
        (void)wavel_threadidx;
        (void)native_gradients;
        throw std::logic_error("Rust support is disabled");
#endif
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::finalize_vjp(
        int wavelidx, int wavel_threadidx,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
#ifdef SKTRAN_RUST_SUPPORT
        Eigen::VectorXd solution_cotangent = Eigen::VectorXd::Zero(
            static_cast<Eigen::Index>(m_num_outgoing_values));
        for (const auto& thread_cotangent :
             m_los_solution_cotangents[wavel_threadidx]) {
            solution_cotangent += thread_cotangent;
        }
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
        sasktran2::rust::successive_orders::
            linearize_coefficients_vjp_assembled(
                *solver, as_rust_slice(solution_cotangent));
        accumulate_rust_vjp_products(wavelidx, wavel_threadidx,
                                     native_gradient);
#else
        (void)wavelidx;
        (void)wavel_threadidx;
        (void)native_gradient;
        throw std::logic_error("Rust support is disabled");
#endif
    }

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    void SourceImpl<NSTOKES>::accumulate_rust_vjp_products(
        int wavelidx, int wavel_threadidx,
        Eigen::Ref<Eigen::VectorXd> native_gradient, bool reset_activity,
        bool transport_vjp_precomputed, bool packed_solar_vjp) const {
        auto& solver = m_rust_solvers[wavel_threadidx];
        const auto boundary_gradient = sasktran2::rust::successive_orders::
            boundary_scattering_value_gradient(*solver);

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
        accumulate_rust_transport_vjp_with_first_order(
            wavelidx, wavel_threadidx, native_gradient,
            transport_vjp_precomputed, packed_solar_vjp);
        for (std::size_t source_index = 0;
             source_index < m_initial_sources.size(); ++source_index) {
            if (packed_solar_vjp && source_index == 0) {
                continue;
            }
            m_initial_sources[source_index]->finalize_vjp(
                wavelidx, wavel_threadidx, native_gradient);
        }
        if (reset_activity) {
            m_active_vjp_derivatives[wavel_threadidx] = DerivativeActivity{};
        }
    }
#endif

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::integrated_source(
        const sasktran2::WavelengthBlock<>& block, int losidx, int layeridx,
        int wavel_threadidx, int threadidx,
        const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        sasktran2::WavelengthBlockDual<NSTOKES>& block_source,
        typename SourceTermInterface<NSTOKES>::IntegrationDirection direction)
        const {
        (void)block;
        (void)losidx;
        (void)layeridx;
        (void)wavel_threadidx;
        (void)threadidx;
        (void)layer;
        (void)entrance_weights;
        (void)exit_weights;
        (void)shell_od;
        (void)block_source;
        (void)direction;
        throw std::logic_error(
            "Rust successive-orders LOS integration must use the whole-ray "
            "callback");
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::integrated_source_jvp(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        sasktran2::RadianceJVP<NSTOKES>& source) const {
        (void)wavelidx;
        (void)losidx;
        (void)layeridx;
        (void)wavel_threadidx;
        (void)threadidx;
        (void)layer;
        (void)entrance_weights;
        (void)exit_weights;
        (void)shell_od;
        (void)native_tangent;
        (void)source;
        throw std::logic_error(
            "Rust successive-orders LOS JVP must use the whole-ray callback");
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::end_of_ray_source_jvp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        sasktran2::RadianceJVP<NSTOKES>& source) const {
        (void)wavelidx;
        (void)losidx;
        (void)wavel_threadidx;
        (void)threadidx;
        (void)native_tangent;
        (void)source;
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::integrated_source_vjp(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        const Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::Vector<double, NSTOKES>> source_value,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        (void)wavelidx;
        (void)losidx;
        (void)layeridx;
        (void)wavel_threadidx;
        (void)threadidx;
        (void)layer;
        (void)entrance_weights;
        (void)exit_weights;
        (void)shell_od;
        (void)cotangent;
        (void)source_value;
        (void)native_gradient;
        throw std::logic_error(
            "Rust successive-orders LOS VJP must use the whole-ray callback");
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::end_of_ray_source_vjp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        const Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        (void)wavelidx;
        (void)losidx;
        (void)wavel_threadidx;
        (void)threadidx;
        (void)cotangent;
        (void)native_gradient;
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::end_of_ray_source(
        const sasktran2::WavelengthBlock<>& block, int losidx,
        int wavel_threadidx, int threadidx,
        sasktran2::WavelengthBlockDual<NSTOKES>& block_source) const {
        (void)block;
        (void)losidx;
        (void)wavel_threadidx;
        (void)threadidx;
        (void)block_source;
    }

    template <int NSTOKES>
    void SourceImpl<NSTOKES>::start_of_ray_source(
        const sasktran2::WavelengthBlock<>& block, int losidx,
        int wavel_threadidx, int threadidx,
        sasktran2::WavelengthBlockDual<NSTOKES>& block_source) const {
        (void)threadidx;
#ifdef SKTRAN_RUST_SUPPORT
        {
            if (m_rust_source_interpolators.empty()) {
                throw std::logic_error(
                    "Rust LOS source interpolator is not initialized");
            }
            const auto& atmosphere = m_atmosphere->storage();
            const std::size_t num_atmosphere_points =
                static_cast<std::size_t>(atmosphere.total_extinction.rows());
            const std::size_t atmosphere_values =
                num_atmosphere_points * static_cast<std::size_t>(block.count);
            // A one-wavelength derivative integration follows a packed
            // primal solve after select_*_block_lane has restored that lane
            // into the scalar Rust solver.  Its local block lane is always
            // zero, so indexing the packed output here would incorrectly
            // reuse wavelength zero for lanes 1--3.
            const bool batched =
                m_wavelength_block_capacity > 1 && block.count > 1;
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
    void SourceImpl<NSTOKES>::start_of_ray_source_jvp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        sasktran2::RadianceJVP<NSTOKES>& source) const {
        (void)threadidx;
#ifdef SKTRAN_RUST_SUPPORT
        {
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
    void SourceImpl<NSTOKES>::start_of_ray_source_vjp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        const Eigen::Vector<double, NSTOKES>& value_before,
        Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        (void)value_before;
#ifdef SKTRAN_RUST_SUPPORT
        {
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

    template class SourceImpl<1>;
    template class SourceImpl<3>;

} // namespace sasktran2::successive_orders
