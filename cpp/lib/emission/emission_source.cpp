#include <sasktran2/raytracing.h>
#include "sasktran2/source_algorithms.h"
#include "sasktran2/viewinggeometry_internal.h"
#include <sasktran2/emission_source.h>
#include <sasktran2/dual.h>
#include <sasktran2/config.h>
#include <sasktran2/atmosphere/atmosphere.h>

#include <stdexcept>

namespace sasktran2::emission {
    template <int NSTOKES, Config::EmissionSource EMISSION_SOURCE_TYPE>
    void EmissionSource<NSTOKES, EMISSION_SOURCE_TYPE>::initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        // Store the rays for later
        m_los_rays = &internal_viewing.traced_rays;
        m_los_rays_2d = nullptr;
    }

    template <int NSTOKES, Config::EmissionSource EMISSION_SOURCE_TYPE>
    void EmissionSource<NSTOKES, EMISSION_SOURCE_TYPE>::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay2D>& traced_rays,
        const sasktran2::Geometry2D& geometry) {
        m_los_rays = nullptr;
        m_los_rays_2d = &traced_rays;
    }

    template <int NSTOKES, Config::EmissionSource EMISSION_SOURCE_TYPE>
    void EmissionSource<NSTOKES, EMISSION_SOURCE_TYPE>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        if (atmosphere.num_deriv() > 0 &&
            !atmosphere.include_emission_derivatives()) {
            throw std::invalid_argument(
                "EmissionSource requires emission derivative storage when "
                "atmospheric derivatives are enabled");
        }

        // Store the atmosphere for later
        m_atmosphere = &atmosphere;
    }

    template <int NSTOKES, Config::EmissionSource EMISSION_SOURCE_TYPE>
    template <typename EntranceWeights, typename ExitWeights>
    void
    EmissionSource<NSTOKES, EMISSION_SOURCE_TYPE>::integrated_source_constant(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::LayerGeometry& layer,
        const EntranceWeights& entrance_weights,
        const ExitWeights& exit_weights,
        const sasktran2::SparseODDualView& shell_od,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {
        // Integrates assuming the source is constant in the layer and
        // determined by the average of the layer boundaries

        bool calculate_derivative = source.derivative_size() > 0;

        double ssa_start = 0;
        double ssa_end = 0;
        double emission_start = 0;
        double emission_end = 0;

        // Calculate SSA and emission at the layer boundaries
        for (std::size_t index = 0; index < entrance_weights.size(); ++index) {
            const auto ele = entrance_weights[index];
            ssa_start +=
                m_atmosphere->storage().ssa(ele.first, wavelidx) * ele.second;
            emission_start +=
                m_atmosphere->storage().emission_source(ele.first, wavelidx) *
                ele.second;
        }
        for (std::size_t index = 0; index < exit_weights.size(); ++index) {
            const auto ele = exit_weights[index];
            ssa_end +=
                m_atmosphere->storage().ssa(ele.first, wavelidx) * ele.second;
            emission_end +=
                m_atmosphere->storage().emission_source(ele.first, wavelidx) *
                ele.second;
        }

        double source_factor1;
        double emission_cell;

        if constexpr (EMISSION_SOURCE_TYPE ==
                      Config::EmissionSource::standard) {
            // Average value of layer boundaries
            source_factor1 = (1 - shell_od.exp_minus_od);
            emission_cell = source_factor1 * ((1 - ssa_start) * emission_start *
                                                  layer.od_quad_start_fraction +
                                              (1 - ssa_end) * emission_end *
                                                  layer.od_quad_end_fraction);

        } else if constexpr (EMISSION_SOURCE_TYPE ==
                             Config::EmissionSource::volume_emission_rate) {
            // Volume emission rate source, so no (1 - ssa) term, and the source
            // factor is the cell path length
            source_factor1 = layer.layer_distance;
            emission_cell = source_factor1 *
                            (emission_start * layer.od_quad_start_fraction +
                             emission_end * layer.od_quad_end_fraction);

        } else {
            // Cant be here
        }

        if constexpr (NSTOKES == 1) {
            source.value.array() += emission_cell;
        } else {
            source.value(0) += emission_cell;
        }

        if (calculate_derivative) {
            if constexpr (EMISSION_SOURCE_TYPE ==
                          Config::EmissionSource::standard) {
                // For standard emission source, the derivatives are more
                // complex Now for the derivatives, start with dsource_factor
                // which is sparse
                Eigen::Ref<Eigen::Matrix<double, NSTOKES, -1>> d_ssa =
                    source.d_ssa(m_atmosphere->storage().ssa.rows());
                Eigen::Ref<Eigen::Matrix<double, NSTOKES, -1>> d_emission =
                    source.d_emission(
                        m_atmosphere->storage().ssa.rows(),
                        m_atmosphere->num_scattering_deriv_groups());

                for (auto it = shell_od.deriv_iter; it; ++it) {
                    source.deriv(0, it.index()) +=
                        it.value() * (1 - source_factor1) *
                        ((1 - ssa_start) * emission_start *
                             layer.od_quad_start_fraction +
                         (1 - ssa_end) * emission_end *
                             layer.od_quad_end_fraction);
                }

                // And the SSA/emission derivatives
                for (std::size_t index = 0; index < entrance_weights.size();
                     ++index) {
                    const auto ele = entrance_weights[index];
                    d_ssa(0, ele.first) -= ele.second * emission_start *
                                           source_factor1 *
                                           layer.od_quad_start_fraction;
                    d_emission(0, ele.first) += ele.second * (1 - ssa_start) *
                                                source_factor1 *
                                                layer.od_quad_start_fraction;
                }
                for (std::size_t index = 0; index < exit_weights.size();
                     ++index) {
                    const auto ele = exit_weights[index];
                    d_ssa(0, ele.first) -= ele.second * emission_end *
                                           source_factor1 *
                                           layer.od_quad_end_fraction;
                    d_emission(0, ele.first) += ele.second * (1 - ssa_end) *
                                                source_factor1 *
                                                layer.od_quad_end_fraction;
                }
            } else if constexpr (EMISSION_SOURCE_TYPE ==
                                 Config::EmissionSource::volume_emission_rate) {
                Eigen::Ref<Eigen::Matrix<double, NSTOKES, -1>> d_emission =
                    source.d_emission(
                        m_atmosphere->storage().ssa.rows(),
                        m_atmosphere->num_scattering_deriv_groups());
                // And the emission derivatives
                for (std::size_t index = 0; index < entrance_weights.size();
                     ++index) {
                    const auto ele = entrance_weights[index];
                    d_emission(0, ele.first) += ele.second * source_factor1 *
                                                layer.od_quad_start_fraction;
                }
                for (std::size_t index = 0; index < exit_weights.size();
                     ++index) {
                    const auto ele = exit_weights[index];
                    d_emission(0, ele.first) += ele.second * source_factor1 *
                                                layer.od_quad_end_fraction;
                }
            } else {
                // Cant be here
            }
        }
    }

    template <int NSTOKES, Config::EmissionSource EMISSION_SOURCE_TYPE>
    void EmissionSource<NSTOKES, EMISSION_SOURCE_TYPE>::integrated_source(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
        const sasktran2::SparseODDualView& shell_od,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source,
        typename SourceTermInterface<NSTOKES>::IntegrationDirection direction)
        const {
        if (layer.layer_distance < MINIMUM_SHELL_SIZE_M) {
            // Essentially an empty shell from rounding, don't have to do
            // anything
            return;
        }

        integrated_source_constant(
            wavelidx, losidx, layeridx, wavel_threadidx, threadidx, layer,
            layer.entrance_interpolation_weights,
            layer.exit_interpolation_weights, shell_od, source);
    }

    template <int NSTOKES, Config::EmissionSource EMISSION_SOURCE_TYPE>
    void EmissionSource<NSTOKES, EMISSION_SOURCE_TYPE>::integrated_source(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::StructuredLayer2D& layer,
        const sasktran2::Geometry2D& geometry,
        const sasktran2::SparseODDualView& shell_od,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source,
        typename SourceTermInterface<NSTOKES>::IntegrationDirection direction)
        const {
        if (layer.layer_distance < MINIMUM_SHELL_SIZE_M) {
            return;
        }

        const std::array<int, 4> location_indices = {
            geometry.location_index(layer.altitude_cell, layer.horizontal_cell),
            geometry.location_index(layer.altitude_cell + 1,
                                    layer.horizontal_cell),
            geometry.location_index(layer.altitude_cell,
                                    layer.horizontal_cell + 1),
            geometry.location_index(layer.altitude_cell + 1,
                                    layer.horizontal_cell + 1)};
        const auto entrance_local_weights =
            layer.entrance_interpolation.weights();
        const auto exit_local_weights = layer.exit_interpolation.weights();
        std::array<std::pair<int, double>, 4> entrance_weights;
        std::array<std::pair<int, double>, 4> exit_weights;
        for (int index = 0; index < 4; ++index) {
            entrance_weights[index] = std::make_pair(
                location_indices[index], entrance_local_weights[index]);
            exit_weights[index] = std::make_pair(location_indices[index],
                                                 exit_local_weights[index]);
        }

        integrated_source_constant(wavelidx, losidx, layeridx, wavel_threadidx,
                                   threadidx, layer, entrance_weights,
                                   exit_weights, shell_od, source);
    }

    template <int NSTOKES, Config::EmissionSource EMISSION_SOURCE_TYPE>
    void EmissionSource<NSTOKES, EMISSION_SOURCE_TYPE>::end_of_ray_source(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {
        const bool ground_is_hit =
            m_los_rays != nullptr ? m_los_rays->at(losidx).ground_is_hit
                                  : m_los_rays_2d->at(losidx).ground_is_hit;
        if (ground_is_hit) {
            double emission_surface =
                m_atmosphere->surface().emission()[wavelidx];
            if constexpr (NSTOKES == 1) {
                source.value.array() += emission_surface;
            } else {
                source.value(0) += emission_surface;
            }

            if (source.derivative_size() > 0) {
                source.deriv(
                    0, m_atmosphere->surface_emission_deriv_start_index()) += 1;
            }
        }
    }

    template class EmissionSource<1, Config::EmissionSource::standard>;
    template class EmissionSource<3, Config::EmissionSource::standard>;

    template class EmissionSource<1,
                                  Config::EmissionSource::volume_emission_rate>;
    template class EmissionSource<3,
                                  Config::EmissionSource::volume_emission_rate>;

} // namespace sasktran2::emission
