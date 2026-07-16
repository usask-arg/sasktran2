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
        m_ground_is_hit.resize(internal_viewing.traced_rays.size());
        std::transform(internal_viewing.traced_rays.begin(),
                       internal_viewing.traced_rays.end(),
                       m_ground_is_hit.begin(),
                       [](const auto& ray) { return ray.ground_is_hit; });
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
    void EmissionSource<NSTOKES, EMISSION_SOURCE_TYPE>::integrated_source(
        const sasktran2::WavelengthBlock& batch, int losidx, int layeridx,
        int wavel_threadidx, int threadidx,
        const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        sasktran2::WavelengthBlockDual<NSTOKES>& source,
        typename SourceTermInterface<NSTOKES>::IntegrationDirection direction)
        const {
        if (layer.layer_distance < MINIMUM_SHELL_SIZE_M) {
            return;
        }

        const auto attenuation = shell_od.exp_minus_od();
        for (int lane = 0; lane < batch.count; ++lane) {
            const int wavelidx = batch.wavelength(lane);
            double ssa_start = 0.0;
            double ssa_end = 0.0;
            double emission_start = 0.0;
            double emission_end = 0.0;
            for (std::size_t index = 0; index < entrance_weights.size();
                 ++index) {
                const auto weight = entrance_weights[index];
                if (weight.second == 0.0) {
                    continue;
                }
                ssa_start +=
                    m_atmosphere->storage().ssa(weight.first, wavelidx) *
                    weight.second;
                emission_start += m_atmosphere->storage().emission_source(
                                      weight.first, wavelidx) *
                                  weight.second;
            }
            for (std::size_t index = 0; index < exit_weights.size(); ++index) {
                const auto weight = exit_weights[index];
                if (weight.second == 0.0) {
                    continue;
                }
                ssa_end += m_atmosphere->storage().ssa(weight.first, wavelidx) *
                           weight.second;
                emission_end += m_atmosphere->storage().emission_source(
                                    weight.first, wavelidx) *
                                weight.second;
            }

            double source_factor;
            double emission_cell;
            if constexpr (EMISSION_SOURCE_TYPE ==
                          Config::EmissionSource::standard) {
                source_factor = 1.0 - attenuation(lane);
                emission_cell =
                    source_factor * ((1.0 - ssa_start) * emission_start *
                                         layer.od_quad_start_fraction +
                                     (1.0 - ssa_end) * emission_end *
                                         layer.od_quad_end_fraction);
            } else {
                source_factor = layer.layer_distance;
                emission_cell = source_factor *
                                (emission_start * layer.od_quad_start_fraction +
                                 emission_end * layer.od_quad_end_fraction);
            }
            source.value(0, lane) += emission_cell;

            if (source.derivative_size() == 0) {
                continue;
            }
            if constexpr (EMISSION_SOURCE_TYPE ==
                          Config::EmissionSource::standard) {
                for (auto derivative = shell_od.derivative_iterator();
                     derivative; ++derivative) {
                    source.derivative(derivative.index(), batch.count)(
                        0, lane) += derivative.value() * attenuation(lane) *
                                    ((1.0 - ssa_start) * emission_start *
                                         layer.od_quad_start_fraction +
                                     (1.0 - ssa_end) * emission_end *
                                         layer.od_quad_end_fraction);
                }
            }

            const auto add_endpoint_derivatives =
                [&](const auto& weights, double ssa, double emission,
                    double quadrature_fraction) {
                    for (std::size_t index = 0; index < weights.size();
                         ++index) {
                        const auto weight = weights[index];
                        if (weight.second == 0.0) {
                            continue;
                        }
                        if constexpr (EMISSION_SOURCE_TYPE ==
                                      Config::EmissionSource::standard) {
                            source.derivative(
                                m_atmosphere->ssa_deriv_start_index() +
                                    weight.first,
                                batch.count)(0, lane) -=
                                weight.second * emission * source_factor *
                                quadrature_fraction;
                        }
                        source.derivative(
                            m_atmosphere->emission_deriv_start_index() +
                                weight.first,
                            batch.count)(0, lane) +=
                            weight.second *
                            (EMISSION_SOURCE_TYPE ==
                                     Config::EmissionSource::standard
                                 ? 1.0 - ssa
                                 : 1.0) *
                            source_factor * quadrature_fraction;
                    }
                };
            add_endpoint_derivatives(entrance_weights, ssa_start,
                                     emission_start,
                                     layer.od_quad_start_fraction);
            add_endpoint_derivatives(exit_weights, ssa_end, emission_end,
                                     layer.od_quad_end_fraction);
        }
    }

    template <int NSTOKES, Config::EmissionSource EMISSION_SOURCE_TYPE>
    void EmissionSource<NSTOKES, EMISSION_SOURCE_TYPE>::end_of_ray_source(
        const sasktran2::WavelengthBlock& batch, int losidx,
        int wavel_threadidx, int threadidx,
        sasktran2::WavelengthBlockDual<NSTOKES>& source) const {
        if (!m_ground_is_hit.at(losidx)) {
            return;
        }
        for (int lane = 0; lane < batch.count; ++lane) {
            source.value(0, lane) +=
                m_atmosphere->surface().emission()[batch.wavelength(lane)];
            if (source.derivative_size() > 0) {
                source.derivative(
                    m_atmosphere->surface_emission_deriv_start_index(),
                    batch.count)(0, lane) += 1.0;
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
