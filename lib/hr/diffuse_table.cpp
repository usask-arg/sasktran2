#include "sasktran2/geometry.h"
#include <algorithm>
#ifdef SKTRAN_OPENMP_SUPPORT
#include <omp.h>
#endif
#include <sasktran2/hr/diffuse_source.h>
#include <sasktran2/math/unitsphere.h>
#include <sasktran2/math/scattering.h>
#include <sasktran2/solartransmission.h>
#include <sasktran2/do_source.h>

#include <fstream>

namespace sasktran2::hr {

    template <int NSTOKES>
    DiffuseTable<NSTOKES>::DiffuseTable(
        const sasktran2::raytracing::RayTracerBase& ray_tracer,
        const sasktran2::Geometry1D& geometry)
        : m_raytracer(ray_tracer), m_geometry(geometry), m_integrator(false) {}

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
        // TODO: Decouple altitude grid here from the global geometry grid ?

        // Put diffuse points inbetween altitude levels
        // TODO: why though? sources are usually calculated at the levels not
        // inbetween
        Eigen::VectorXd alt_values =
            (m_geometry.altitude_grid().grid()(Eigen::seq(0, Eigen::last - 1)) +
             m_geometry.altitude_grid().grid()(Eigen::seq(1, Eigen::last))) /
            2.0;

        return sasktran2::grids::AltitudeGrid(
            std::move(alt_values), sasktran2::grids::gridspacing::constant,
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
                m_config->num_hr_outgoing())));

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
                    location));
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

        m_incoming_traced_rays.resize(start_incoming_idx);

        for (auto& storage : m_thread_storage) {
            storage.m_incoming_radiances.resize(start_incoming_idx * NSTOKES, 0,
                                                false);
            storage.m_firstorder_radiances.resize(start_incoming_idx * NSTOKES,
                                                  0, false);
            storage.m_outgoing_sources.resize(start_outgoing_idx * NSTOKES, 0,
                                              false);

            storage.point_scattering_matrices.resize(m_diffuse_points.size());
            for (int i = 0; i < m_diffuse_points.size(); ++i) {
                if (m_diffuse_point_full_calculation[i]) {
                    storage.point_scattering_matrices[i].resize(
                        m_diffuse_points[i]->num_outgoing() * NSTOKES,
                        m_diffuse_points[i]->num_incoming() * NSTOKES);
                }
            }
        }
    }

    template <int NSTOKES> void DiffuseTable<NSTOKES>::trace_incoming_rays() {
        int nthreads = m_config->num_threads();

        std::vector<sasktran2::viewinggeometry::ViewingRay> thread_viewing_ray;
        thread_viewing_ray.resize(nthreads);

#pragma omp parallel for num_threads(nthreads)
        for (int i = 0; i < m_diffuse_points.size(); ++i) {
#ifdef SKTRAN_OPENMP_SUPPORT
            auto& viewing_ray = thread_viewing_ray[omp_get_thread_num()];
#else
            auto& viewing_ray = thread_viewing_ray[0];
#endif
            if (!m_diffuse_point_full_calculation[i]) {
                continue;
            }
            viewing_ray.observer = m_diffuse_points[i]->location();
            for (int j = 0; j < m_diffuse_points[i]->num_incoming(); ++j) {
                viewing_ray.look_away =
                    m_diffuse_points[i]->incoming_direction(j);

                m_raytracer.trace_ray(
                    viewing_ray,
                    m_incoming_traced_rays[m_diffuse_incoming_index_map[i] + j],
                    m_config->multiple_scatter_refraction());
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& los_rays) {
        // TODO: This is the 1D case, need to add separate logic for 2d/3d

        // find the min/max SZA from the LOS rays and generate the cos_sza_grid
        std::pair<double, double> min_max_cos_sza =
            sasktran2::raytracing::min_max_cos_sza_of_all_rays(los_rays);

        // create the location interpolator
        m_location_interpolator = std::make_unique<
            sasktran2::grids::AltitudeSZASourceLocationInterpolator>(
            generate_altitude_grid(),
            generate_cos_sza_grid(min_max_cos_sza.first,
                                  min_max_cos_sza.second));

        // Construct the actual diffuse points
        construct_diffuse_points();

        // Trace all of the incoming rays
        trace_incoming_rays();

        // Set up the integrator
        m_integrator.initialize_geometry(m_incoming_traced_rays,
                                         this->m_geometry);
        // And the initial sources
        // This is a little tricky, any source that is used for the incoming
        // rays needs to be initialized with the traced incoming rays
        for (auto& source : m_initial_sources) {
            source->initialize_geometry(m_incoming_traced_rays);
        }

        // But the DO Source should be initialized with the LOS rays
        if (m_config->initialize_hr_with_do()) {
            m_do_source->initialize_geometry(los_rays);
        }

        int temp;
        generate_source_interpolation_weights(m_incoming_traced_rays,
                                              m_diffuse_source_weights,
                                              m_total_num_diffuse_weights);

        construct_accumulation_sparsity();

        generate_source_interpolation_weights(los_rays, m_los_source_weights,
                                              temp);

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
    void DiffuseTable<NSTOKES>::construct_accumulation_sparsity() {
        // These are the length of the rows
        m_outer_starts.resize(NSTOKES * m_diffuse_source_weights.size() + 1);
        m_inner_nnz.resize(NSTOKES * m_diffuse_source_weights.size());

        std::vector<Eigen::VectorXi> inner_indicies;

        inner_indicies.resize(NSTOKES * m_diffuse_source_weights.size());

        std::vector<
            std::tuple<int, std::tuple<int, double, std::array<int, NSTOKES>>*>>
            sorting_helper;

        int total_nnz = 0;
        // Construct the matrix NSTOKES rows at a time
        for (int row = 0; row < m_diffuse_source_weights.size(); ++row) {
            sorting_helper.clear();
            auto& weights = m_diffuse_source_weights[row];

            // Start by setting the index storage to be equal to the weight
            // index
            for (int i = 0; i < weights.interior_weights.size(); ++i) {
                auto& weight = weights.interior_weights[i].second;

                for (int j = 0; j < weight.size(); ++j) {
                    sorting_helper.push_back(
                        std::make_tuple(std::get<0>(weight[j]), &weight[j]));
                }
            }
            // And for the ground weights
            for (int i = 0; i < weights.ground_weights.size(); ++i) {
                auto& weight = weights.ground_weights[i];

                sorting_helper.push_back(
                    std::make_tuple(std::get<0>(weight), &weight));
            }

            // Now we have to sort on the column indices
            std::stable_sort(
                std::begin(sorting_helper), std::end(sorting_helper),
                [](const std::tuple<
                       int, std::tuple<int, double, std::array<int, NSTOKES>>*>&
                       left,
                   const std::tuple<
                       int, std::tuple<int, double, std::array<int, NSTOKES>>*>&
                       right) {
                    return std::get<0>(left) < std::get<0>(right);
                });

            if (sorting_helper.size() == 0) {
                // 0 Elements in this row, not sure what this means? Probably
                // never actually happens
                for (int s = 0; s < NSTOKES; ++s) {
                    m_inner_nnz(row * NSTOKES + s) = 0;
                }
                continue;
            }

            // Now we have to do a first pass through the sorted list to find
            // the number of unique elements (nnz)
            int nnz = 1;
            for (int i = 1; i < sorting_helper.size(); ++i) {
                if (std::get<0>(sorting_helper[i]) !=
                    std::get<0>(sorting_helper[i - 1])) {
                    ++nnz;
                }
            }

            // Go through our NSTOKES rows
            for (int s = 0; s < NSTOKES; ++s) {
                auto& inner = inner_indicies[row * NSTOKES + s];
                // First assign nnz
                m_inner_nnz(row * NSTOKES + s) = nnz;

                inner.resize(nnz);

                int inner_index = 0;
                // Now go through the sorting helper and assign our indicies
                for (int i = 0; i < sorting_helper.size(); ++i) {
                    if (i > 0) {
                        if (std::get<0>(sorting_helper[i]) !=
                            std::get<0>(sorting_helper[i - 1])) {
                            ++inner_index;
                        }
                    }
                    inner[inner_index] =
                        std::get<0>(*std::get<1>(sorting_helper[i])) * NSTOKES +
                        s;

                    std::get<2>(*std::get<1>(sorting_helper[i]))[s] =
                        inner_index + total_nnz;
                }
                total_nnz += nnz;
            }
        }

        // Now we have to copy our row inner indicies to the full matrix
        m_inner_indicies.resize(total_nnz);
        int current_index = 0;
        for (auto& index : inner_indicies) {
            for (int i = 0; i < index.size(); ++i) {
                m_inner_indicies[current_index] = index[i];
                ++current_index;
            }
        }

        // Lastly set the outerstarts vector
        m_outer_starts(0) = 0;
        for (int i = 0; i < inner_indicies.size(); ++i) {
            m_outer_starts(i + 1) =
                m_outer_starts(i) + inner_indicies[i].size();
        }

        // And we can now resize our value storage
        for (auto& storage : m_thread_storage) {
            storage.accumulation_value_storage.resize(
                m_config->num_source_threads());

            for (auto& vec : storage.accumulation_value_storage) {
                vec.resize(total_nnz);
            }
            storage.accumulation_summed_values.resize(total_nnz);
        }
    }

    template <int NSTOKES>
    void
    DiffuseTable<NSTOKES>::initialize_config(const sasktran2::Config& config) {
        m_config = &config;

        m_thread_storage.resize(m_config->num_wavelength_threads());

        m_initial_owned_sources.emplace_back(
            std::make_unique<sasktran2::solartransmission::SingleScatterSource<
                sasktran2::solartransmission::SolarTransmissionTable, NSTOKES>>(
                m_geometry, m_raytracer));

        m_initial_sources.push_back(m_initial_owned_sources[0].get());

        if (m_config->initialize_hr_with_do()) {
            m_initial_owned_sources.emplace_back(
                std::make_unique<
                    sasktran2::DOSourceInterpolatedPostProcessing<NSTOKES, -1>>(
                    m_geometry, m_raytracer, false));

            m_do_source =
                static_cast<DOSourceInterpolatedPostProcessing<NSTOKES, -1>*>(
                    m_initial_owned_sources[1].get());
        }

        for (auto& source : m_initial_owned_sources) {
            source->initialize_config(config);
        }
    }

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
            new_position, saa_initial, temp.cos_zenith_angle(vector));
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::generate_source_interpolation_weights(
        const std::vector<sasktran2::raytracing::TracedRay>& rays,
        SInterpolator& interpolator, int& total_num_weights) const {
        total_num_weights = 0;

        interpolator.resize(rays.size());

        int nthreads = m_config->num_threads();

        std::vector<std::vector<std::pair<int, double>>>
            thread_temp_location_storage;
        std::vector<std::vector<std::pair<int, double>>>
            thread_temp_direction_storage;

        thread_temp_location_storage.resize(nthreads);
        thread_temp_direction_storage.resize(nthreads);

        int num_location, num_direction;

        Eigen::Vector3d rotated_los;

        sasktran2::Location temp_location;

#pragma omp parallel for num_threads(nthreads) private(                        \
        num_location, num_direction, rotated_los, temp_location)
        for (int rayidx = 0; rayidx < rays.size(); ++rayidx) {
#ifdef SKTRAN_OPENMP_SUPPORT
            int threadidx = omp_get_thread_num();
#else
            int threadidx = 0;
#endif

            auto& temp_location_storage =
                thread_temp_location_storage[threadidx];
            auto& temp_direction_storage =
                thread_temp_direction_storage[threadidx];

            auto& ray_interpolator = interpolator[rayidx];
            const auto& ray = rays[rayidx];

            ray_interpolator.interior_weights.resize(ray.layers.size());

            for (int layeridx = 0; layeridx < ray.layers.size(); ++layeridx) {
                const auto& layer = ray.layers[layeridx];
                auto& layer_interpolator =
                    ray_interpolator.interior_weights[layeridx].second;
                auto& atmosphere_interpolator =
                    ray_interpolator.interior_weights[layeridx].first;

                temp_location.position =
                    (layer.entrance.position + layer.exit.position) / 2.0;

                m_geometry.assign_interpolation_weights(
                    temp_location, atmosphere_interpolator);

                m_location_interpolator->interior_interpolation_weights(
                    m_geometry.coordinates(), temp_location,
                    temp_location_storage, num_location);
                for (int locidx = 0; locidx < num_location; ++locidx) {
                    const auto& contributing_point =
                        m_diffuse_points[temp_location_storage[locidx].first];

                    if (m_geometry.coordinates().geometry_type() ==
                        sasktran2::geometrytype::spherical) {
                        // Multiply by -1 to get the propagation direction
                        rotated_los = rotate_unit_vector(
                            -1 * layer.average_look_away,
                            temp_location.position,
                            contributing_point->location().position);
                    } else {
                        // Don't need to rotate in plane parallel geometry
                        rotated_los = -1 * layer.average_look_away;
                    }

                    contributing_point->sphere_pair()
                        .outgoing_sphere()
                        .interpolate(rotated_los, temp_direction_storage,
                                     num_direction);

                    for (int diridx = 0; diridx < num_direction; ++diridx) {
#ifdef SASKTRAN_DEBUG_ASSERTS
                        if (m_diffuse_outgoing_index_map
                                    [temp_location_storage[locidx].first] +
                                temp_direction_storage[diridx].first >
                            m_diffuse_outgoing_index_map.back() +
                                m_config->num_hr_outgoing()) {
                            spdlog::error("BAD INDEX {} {}",
                                          temp_location_storage[locidx].first,
                                          temp_direction_storage[diridx].first);
                        }
#endif

                        layer_interpolator.emplace_back(std::make_tuple(
                            m_diffuse_outgoing_index_map
                                    [temp_location_storage[locidx].first] +
                                temp_direction_storage[diridx].first,
                            temp_location_storage[locidx].second *
                                temp_direction_storage[diridx].second,
                            std::array<int, NSTOKES>{}));
                    }
                }
                total_num_weights += num_location * num_direction;
            }

            ray_interpolator.ground_is_hit = ray.ground_is_hit;
            if (ray_interpolator.ground_is_hit) {
                const auto& layer = ray.layers[0];

                temp_location.position = layer.exit.position;

                m_location_interpolator->ground_interpolation_weights(
                    m_geometry.coordinates(), temp_location,
                    temp_location_storage, num_location);

                for (int locidx = 0; locidx < num_location; ++locidx) {
                    const auto& contributing_point =
                        m_diffuse_points[temp_location_storage[locidx].first];

                    if (m_geometry.coordinates().geometry_type() ==
                        sasktran2::geometrytype::spherical) {
                        // Multiply by -1 to get the propagation direction
                        rotated_los = rotate_unit_vector(
                            -1 * layer.average_look_away,
                            temp_location.position,
                            contributing_point->location().position);
                    } else {
                        rotated_los = -1 * layer.average_look_away;
                    }

                    contributing_point->sphere_pair()
                        .outgoing_sphere()
                        .interpolate(rotated_los, temp_direction_storage,
                                     num_direction);

                    for (int diridx = 0; diridx < num_direction; ++diridx) {
                        ray_interpolator.ground_weights.emplace_back(
                            std::make_tuple(
                                m_diffuse_outgoing_index_map
                                        [temp_location_storage[locidx].first] +
                                    temp_direction_storage[diridx].first,
                                temp_location_storage[locidx].second *
                                    temp_direction_storage[diridx].second,
                                std::array<int, NSTOKES>{}));
                    }
                }
                total_num_weights += num_location * num_direction;
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::generate_scattering_matrices(int wavelidx,
                                                             int threadidx) {
        auto& storage = m_thread_storage[threadidx];

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
        for (int i = 0; i < m_location_interpolator->num_ground_points(); ++i) {
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
                    .point_scattering_matrices[i + m_location_interpolator
                                                       ->num_interior_points()]
                    .data());
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::generate_accumulation_matrix(int wavelidx,
                                                             int threadidx) {}

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::iterate_to_solution(int wavelidx,
                                                    int threadidx) {
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

                storage.m_outgoing_sources.value(outgoing_seq).noalias() =
                    storage.point_scattering_matrices[i] *
                    storage.m_incoming_radiances.value(incoming_seq);

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
    void DiffuseTable<NSTOKES>::calculate(int wavelidx, int threadidx) {
        for (auto& source : m_initial_owned_sources) {
            source->calculate(wavelidx, threadidx);
        }

        int nthreads = m_config->num_source_threads();

        if (m_config->num_hr_spherical_iterations() > 0) {
            // Calculate the first order incoming signal
            std::vector<
                sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>>
                temp_result;
            temp_result.resize(nthreads);

            for (int i = 0; i < nthreads; ++i) {
                m_thread_storage[threadidx]
                    .accumulation_value_storage[i]
                    .setZero();
            }

#pragma omp parallel for num_threads(nthreads)
            for (int rayidx = 0; rayidx < m_incoming_traced_rays.size();
                 ++rayidx) {
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

                m_integrator.integrate_and_emplace_accumulation_triplets(
                    temp_result[ray_threadidx], m_initial_sources, wavelidx,
                    rayidx, threadidx, ray_threadidx + threadidx,
                    m_diffuse_source_weights,
                    m_thread_storage[threadidx]
                        .accumulation_value_storage[ray_threadidx]);

                m_thread_storage[threadidx].m_firstorder_radiances.value(
                    Eigen::seq(rayidx * NSTOKES,
                               rayidx * NSTOKES + NSTOKES - 1)) =
                    temp_result[ray_threadidx].value;

#ifdef SASKTRAN_DEBUG_ASSERTS
                if (temp_result.value.hasNaN()) {
                    spdlog::error("Incoming Ray: {} has NaN", rayidx);
                }
#endif
            }

            m_thread_storage[threadidx].accumulation_summed_values.setZero();
            for (int i = 0; i < nthreads; ++i) {
                m_thread_storage[threadidx].accumulation_summed_values +=
                    m_thread_storage[threadidx].accumulation_value_storage[i];
            }

            // Else, start with total first order incoming
            m_thread_storage[threadidx].m_incoming_radiances.value =
                m_thread_storage[threadidx].m_firstorder_radiances.value;

            // Generate the scattering and accumulation matrices
            generate_scattering_matrices(wavelidx, threadidx);
            generate_accumulation_matrix(wavelidx, threadidx);
        }

        // Iterate to the solution
        iterate_to_solution(wavelidx, threadidx);
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::integrated_source(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
        const sasktran2::SparseODDualView& shell_od,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {
        auto& storage = m_thread_storage[wavel_threadidx];

        auto& interpolator =
            m_los_source_weights[losidx].interior_weights[layeridx];

        // Start by calculating ssa at the source point
        double omega = 0;
        for (int i = 0; i < interpolator.first.size(); ++i) {
            auto& index_weight = interpolator.first[i];
            omega += m_atmosphere->storage().ssa(index_weight.first, wavelidx) *
                     index_weight.second;
        }

        double source_factor = (1 - exp(-1 * shell_od.od));

        for (int i = 0; i < interpolator.second.size(); ++i) {
            auto& index_weight = interpolator.second[i];

            for (int s = 0; s < NSTOKES; ++s) {
                double source_value =
                    storage.m_outgoing_sources.value(
                        (std::get<0>(index_weight) * NSTOKES + s)) *
                    std::get<1>(index_weight);

                source.value(s) += omega * source_factor * source_value;

                if (m_atmosphere->num_deriv() > 0) {
                    // Now we need dJ/dthickness
                    for (auto it = shell_od.deriv_iter; it; ++it) {
                        source.deriv(s, it.index()) += it.value() *
                                                       (1 - source_factor) *
                                                       source_value * omega;
                    }

                    // And dJ/dssa
                    for (auto& ele : interpolator.first) {
                        source.deriv(s, m_atmosphere->ssa_deriv_start_index() +
                                            ele.first) +=
                            ele.second * source_factor * source_value;
                    }

                    if (this->m_config->wf_precision() ==
                            sasktran2::Config::WeightingFunctionPrecision::
                                full &&
                        m_config->initialize_hr_with_do()) {
                        source.deriv(s, Eigen::all) +=
                            omega * source_factor * std::get<1>(index_weight) *
                            storage.m_outgoing_sources.deriv(
                                std::get<0>(index_weight) * NSTOKES + s,
                                Eigen::all);
                    }
                }
            }
        }
    }

    template <int NSTOKES>
    void DiffuseTable<NSTOKES>::end_of_ray_source(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {
        auto& interpolator = m_los_source_weights[losidx].ground_weights;
        auto& storage = m_thread_storage[wavel_threadidx];

        for (int i = 0; i < interpolator.size(); ++i) {
            auto& index_weight = interpolator[i];

            for (int s = 0; s < NSTOKES; ++s) {
                double source_value =
                    storage.m_outgoing_sources.value(
                        (std::get<0>(index_weight) * NSTOKES + s)) *
                    std::get<1>(index_weight);

                source.value(s) += source_value;

                if (this->m_config->wf_precision() ==
                        sasktran2::Config::WeightingFunctionPrecision::full &&
                    m_config->initialize_hr_with_do()) {
                    source.deriv(s, Eigen::all) +=
                        std::get<1>(index_weight) *
                        storage.m_outgoing_sources.deriv(
                            std::get<0>(index_weight) * NSTOKES + s,
                            Eigen::all);
                }
            }
        }
    }

    template class DiffuseTable<1>;
    template class DiffuseTable<3>;

} // namespace sasktran2::hr
