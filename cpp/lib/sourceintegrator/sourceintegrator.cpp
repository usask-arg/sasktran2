
#include "sasktran2/source_interface.h"
#include <sasktran2/source_integrator.h>

namespace sasktran2 {
    template <int NSTOKES>
    SourceIntegrator<NSTOKES>::SourceIntegrator(bool calculate_derivatives)
        : m_calculate_derivatives(calculate_derivatives) {}

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& traced_rays,
        const Geometry& geometry) {
        // Construct the optical depth matrices.
        // This is the matrix so that matrix @ extinction = layer od, one matrix
        // for each ray Calculating this matrix beforehand makes calculating
        // derivatives easier, and removes excess computation for every
        // wavelength
        m_traced_ray_od_matrix.resize(traced_rays.size());
        for (int i = 0; i < traced_rays.size(); ++i) {
            sasktran2::raytracing::construct_od_matrix(
                traced_rays[i], geometry, m_traced_ray_od_matrix[i]);
        }

        m_shell_od.resize(traced_rays.size());
        m_exp_minus_shell_od.resize(traced_rays.size());

        m_traced_rays = &traced_rays;
        m_traced_rays_2d = nullptr;
        m_geometry_2d = nullptr;
        m_num_geometry_locations = geometry.size();
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay2D>& traced_rays,
        const Geometry2D& geometry) {
        m_traced_ray_od_matrix.resize(traced_rays.size());
        for (int i = 0; i < traced_rays.size(); ++i) {
            sasktran2::raytracing::construct_od_matrix(
                traced_rays[i], geometry, m_traced_ray_od_matrix[i]);
        }

        m_shell_od.resize(traced_rays.size());
        m_exp_minus_shell_od.resize(traced_rays.size());
        m_traced_rays = nullptr;
        m_traced_rays_2d = &traced_rays;
        m_geometry_2d = &geometry;
        m_num_geometry_locations = geometry.size();
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmo) {
        if (atmo.storage().total_extinction.rows() !=
            m_num_geometry_locations) {
            throw std::invalid_argument(
                "Atmosphere extinction size does not match ray geometry");
        }
// Multithread over LOS? or wavelength? Or just let Eigen do it?
#pragma omp parallel for
        for (int i = 0; i < m_traced_ray_od_matrix.size(); ++i) {
            m_shell_od[i].noalias() =
                m_traced_ray_od_matrix[i] * atmo.storage().total_extinction;

#ifdef SASKTRAN_DEBUG_ASSERTS
            if (!m_shell_od[i].allFinite()) {
                spdlog::error("Error calculating Layer OD for ray: ", i);
            }
#endif
        }

        m_atmosphere = &atmo;

        if (atmo.num_deriv() == 0) {
            m_calculate_derivatives = false;
        }
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::integrate(
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            radiance,
        std::vector<SourceTermInterface<NSTOKES>*> source_terms, int wavelidx,
        int rayidx, int wavel_threadidx, int threadidx) {

        bool have_to_integrate = false;
        for (const auto& source : source_terms) {
            if (source->requires_integration()) {
                have_to_integrate = true;
                if (m_traced_rays == nullptr && source->has_interior_source() &&
                    !source->supports_2d_interior_source()) {
                    throw std::invalid_argument(
                        "Interior source integration is not supported for "
                        "structured 2D rays");
                }
            }
        }

        // Add source at the end of the ray
        for (const auto& source : source_terms) {
            source->end_of_ray_source(wavelidx, rayidx, wavel_threadidx,
                                      threadidx, radiance);
        }

        if (!have_to_integrate) {
            return;
        }
        // Iterate through each OD row from the end of the ray to the observer.
        // Attenuation is geometry-independent; only interior sources require a
        // concrete layer type.
        for (int j = 0; j < m_traced_ray_od_matrix.at(rayidx).rows(); ++j) {
            sasktran2::SparseODDualView local_shell_od(
                m_shell_od.at(rayidx)(j, wavelidx),
                std::exp(-m_shell_od.at(rayidx)(j, wavelidx)),
                m_traced_ray_od_matrix.at(rayidx), j);

            attenuate_layer(radiance, local_shell_od);

            if (m_traced_rays != nullptr) {
                // Calculate all of the 1D layer sources. The capability check
                // above guarantees that 2D integration has no interior
                // sources, so it avoids geometry-specific dispatch entirely.
                const auto& layer = m_traced_rays->at(rayidx).layers.at(j);
                for (const auto& source : source_terms) {
                    source->integrated_source(
                        wavelidx, rayidx, j, wavel_threadidx, threadidx, layer,
                        local_shell_od, radiance,
                        SourceTermInterface<
                            NSTOKES>::IntegrationDirection::backward);
                }
            } else {
                const auto& layer = m_traced_rays_2d->at(rayidx).layers.at(j);
                for (const auto& source : source_terms) {
                    if (source->has_interior_source()) {
                        source->integrated_source(
                            wavelidx, rayidx, j, wavel_threadidx, threadidx,
                            layer, *m_geometry_2d, local_shell_od, radiance,
                            SourceTermInterface<
                                NSTOKES>::IntegrationDirection::backward);
                    }
                }
            }

#ifdef SASKTRAN_DEBUG_ASSERTS
            if (radiance.value.hasNaN()) {
                static bool message = false;
                if (!message) {
                    spdlog::error("One of the sources was  NaN Ray: {} layer: "
                                  "{} Layer od: {} Layer Atten Factor: {}",
                                  rayidx, j, local_shell_od.od,
                                  local_shell_od.exp_minus_od);
                    message = true;
                }
            }
#endif
        }
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::attenuate_layer(
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            radiance,
        const sasktran2::SparseODDualView& local_shell_od) const {
        // rad = rad * atten, drad = drad * atten + rad * datten. Attenuation
        // affects every derivative, while datten only affects extinction.
        if (m_calculate_derivatives) {
            for (auto it = local_shell_od.deriv_iter; it; ++it) {
                radiance.deriv(Eigen::placeholders::all, it.index()) -=
                    it.value() * radiance.value;
            }
        }

        radiance.value *= local_shell_od.exp_minus_od;
        if (m_calculate_derivatives) {
            radiance.deriv *= local_shell_od.exp_minus_od;
        }
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::integrate_optical_depth(
        Eigen::MatrixXd& optical_depth) {
        for (int i = 0; i < m_shell_od.size(); ++i) {
            optical_depth.col(i) = m_shell_od[i].colwise().sum();
        }
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::integrate_and_emplace_accumulation_triplets(
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            radiance,
        std::vector<SourceTermInterface<NSTOKES>*> source_terms, int wavelidx,
        int rayidx, int wavel_threadidx, int threadidx,
        const SInterpolator& source_interpolator,
        Eigen::VectorXd& accumulation_values) {
        ZoneScopedN("Integrate and Emplace Accumulation Triplets");
        const auto& ray = (*m_traced_rays)[rayidx];
        const auto& interpolator = source_interpolator[rayidx];

        // If we don't have to calculate derivatives then it is faster to
        // iterate over the ray backwards, i.e., from the observer to the end of
        // the atmosphere
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>
            layer_source(NSTOKES, 0);

        double current_od = 0;
        for (int j = (int)ray.layers.size() - 1; j >= 0; --j) {
            const sasktran2::raytracing::SphericalLayer& layer = ray.layers[j];

            sasktran2::SparseODDualView local_shell_od(
                m_shell_od[rayidx](j, wavelidx),
                std::exp(-m_shell_od[rayidx](j, wavelidx)),
                m_traced_ray_od_matrix[rayidx], j);
            const auto& layer_interpolator = interpolator.interior_weights[j];
            // Calculate and add the layer source to the radiance
            double atten_factor = std::exp(-current_od);

            // Calculate all of the layer sources
            layer_source.value.setZero();
            for (const auto& source : source_terms) {
                source->integrated_source(
                    wavelidx, rayidx, j, wavel_threadidx, threadidx, layer,
                    local_shell_od, layer_source,
                    SourceTermInterface<
                        NSTOKES>::IntegrationDirection::forward);
            }

            radiance.value += layer_source.value * atten_factor;

            // Assign the accumulation weights
            double omega = 0;
            for (int i = 0; i < layer_interpolator.first.size(); ++i) {
                auto& index_weight = layer_interpolator.first[i];
                omega +=
                    m_atmosphere->storage().ssa(index_weight.first, wavelidx) *
                    index_weight.second;
            }
            double source_factor =
                omega * (1 - local_shell_od.exp_minus_od) * atten_factor;

            for (const auto& ele : layer_interpolator.second) {
                for (int s = 0; s < NSTOKES; ++s) {
                    accumulation_values(std::get<2>(ele)[s]) +=
                        std::get<1>(ele) * source_factor;
                }
            }

            current_od += local_shell_od.od;
        }

        // Add source at the end of the ray
        layer_source.value.setZero();
        for (const auto& source : source_terms) {
            source->end_of_ray_source(wavelidx, rayidx, wavel_threadidx,
                                      threadidx, layer_source);
        }

        radiance.value += layer_source.value * std::exp(-1 * current_od);

        // Add ground interpolation triplets
        if (ray.ground_is_hit) {
            const auto& ground_interpolator = interpolator.ground_weights;

            for (const auto& ele : ground_interpolator) {
                for (int s = 0; s < NSTOKES; ++s) {
                    accumulation_values(std::get<2>(ele)[s]) +=
                        std::get<1>(ele) * std::exp(-1 * current_od);
                }
            }
        }
    }

    template class SourceIntegrator<1>;
    template class SourceIntegrator<3>;
} // namespace sasktran2
