
#include "sasktran2/source_interface.h"
#include <sasktran2/source_integrator.h>

namespace sasktran2 {
    template <int NSTOKES>
    SourceIntegrator<NSTOKES>::SourceIntegrator(bool calculate_derivatives)
        : m_derivatives_enabled(calculate_derivatives),
          m_calculate_derivatives(calculate_derivatives) {}

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::initialize_geometry(
        const std::vector<sasktran2::raytracing::TracedRay>& traced_rays,
        const Geometry& geometry) {
        m_use_sparse_derivative_tracking = false;
        m_attenuation_active_derivative_ranges.clear();

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

        m_traced_rays = &traced_rays;
        m_num_geometry_locations = geometry.size();
        m_num_geometry_dimensions = geometry.num_atmosphere_dimensions();
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmo) {
        m_use_sparse_derivative_tracking = false;
        m_attenuation_active_derivative_ranges.clear();

        if (atmo.storage().total_extinction.rows() !=
            m_num_geometry_locations) {
            throw std::invalid_argument(
                "Atmosphere extinction size does not match ray geometry");
        }
        if (m_wavelength_block_capacity == 1) {
            m_scalar_shell_od.resize(m_traced_ray_od_matrix.size());
            m_shell_od.clear();
        } else {
            m_shell_od.resize(m_traced_ray_od_matrix.size());
            m_scalar_shell_od.clear();
        }

// Multithread over LOS? or wavelength? Or just let Eigen do it?
#pragma omp parallel for
        for (int i = 0; i < m_traced_ray_od_matrix.size(); ++i) {
            if (m_wavelength_block_capacity == 1) {
                m_scalar_shell_od[i].noalias() =
                    m_traced_ray_od_matrix[i] * atmo.storage().total_extinction;
            } else {
                m_shell_od[i].noalias() =
                    m_traced_ray_od_matrix[i] * atmo.storage().total_extinction;
            }

#ifdef SASKTRAN_DEBUG_ASSERTS
            const bool all_finite = m_wavelength_block_capacity == 1
                                        ? m_scalar_shell_od[i].allFinite()
                                        : m_shell_od[i].allFinite();
            if (!all_finite) {
                spdlog::error("Error calculating Layer OD for ray: ", i);
            }
#endif
        }

        m_atmosphere = &atmo;

        // This object may be reused with derivative-free and derivative-enabled
        // atmospheres. Do not let a derivative-free call permanently disable
        // attenuation derivatives for later calculations.
        m_calculate_derivatives = m_derivatives_enabled && atmo.num_deriv() > 0;
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::initialize_derivative_sparsity(
        const std::vector<SourceTermInterface<NSTOKES>*>& source_terms) {
        m_use_sparse_derivative_tracking = false;
        m_attenuation_active_derivative_ranges.clear();

        if (!m_calculate_derivatives || m_traced_rays == nullptr ||
            source_terms.empty()) {
            return;
        }

        for (const auto* source : source_terms) {
            if (!source->supports_sparse_derivative_tracking()) {
                return;
            }
        }

        const auto sort_and_deduplicate = [](std::vector<int>& indices) {
            std::sort(indices.begin(), indices.end());
            indices.erase(std::unique(indices.begin(), indices.end()),
                          indices.end());
        };

        const auto assign_contiguous_ranges =
            [](const std::vector<int>& indices,
               std::vector<std::pair<int, int>>& ranges) {
                ranges.clear();
                for (const int index : indices) {
                    if (ranges.empty() ||
                        index != ranges.back().first + ranges.back().second) {
                        ranges.emplace_back(index, 1);
                    } else {
                        ++ranges.back().second;
                    }
                }
            };

        m_attenuation_active_derivative_ranges.resize(m_traced_rays->size());
        for (std::size_t ray_index = 0; ray_index < m_traced_rays->size();
             ++ray_index) {
            const auto& ray = (*m_traced_rays)[ray_index];
            auto& ray_active =
                m_attenuation_active_derivative_ranges[ray_index];
            ray_active.resize(ray.layers.size());

            std::vector<int> cumulative_active;
            for (const auto* source : source_terms) {
                source->append_end_of_ray_active_derivatives(
                    static_cast<int>(ray_index), cumulative_active);
            }
            sort_and_deduplicate(cumulative_active);

            for (std::size_t layer_index = 0; layer_index < ray.layers.size();
                 ++layer_index) {
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
                         derivative(m_traced_ray_od_matrix[ray_index],
                                    static_cast<int>(layer_index));
                     derivative; ++derivative) {
                    cumulative_active.push_back(derivative.index());
                }
                sort_and_deduplicate(cumulative_active);
                assign_contiguous_ranges(cumulative_active,
                                         ray_active[layer_index]);

                for (const auto* source : source_terms) {
                    if (source->has_interior_source()) {
                        source->append_interior_active_derivatives(
                            static_cast<int>(ray_index),
                            static_cast<int>(layer_index), cumulative_active);
                    }
                }
                sort_and_deduplicate(cumulative_active);
            }
        }

        m_use_sparse_derivative_tracking = true;
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::integrate(
        sasktran2::WavelengthBlockDual<NSTOKES>& radiance,
        const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
        const sasktran2::WavelengthBlock<>& block, int rayidx,
        int wavel_threadidx, int threadidx) {
        dispatch_wavelength_block(block, [&](const auto& fixed_block) {
            integrate_block(radiance, source_terms, fixed_block, rayidx,
                            wavel_threadidx, threadidx);
        });
    }

    template <int NSTOKES>
    template <int N>
    void SourceIntegrator<NSTOKES>::integrate_block(
        sasktran2::WavelengthBlockDual<NSTOKES>& radiance,
        const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
        const sasktran2::WavelengthBlock<N>& block, int rayidx,
        int wavel_threadidx, int threadidx) const {
        if (block.count < 1 || block.count > radiance.block_capacity()) {
            throw std::invalid_argument(
                "Wavelength block does not fit the radiance storage");
        }

        bool have_to_integrate = false;
        for (const auto* source : source_terms) {
            if (source->maximum_wavelength_block_size() < block.count) {
                throw std::invalid_argument(
                    "Source does not support the requested wavelength block "
                    "size");
            }
            if (source->requires_integration()) {
                have_to_integrate = true;
                if (source->has_interior_source() &&
                    !source->supports_geometry_dimension(
                        m_num_geometry_dimensions)) {
                    throw std::invalid_argument(
                        "Interior source integration is not supported for "
                        "the configured atmosphere dimensionality");
                }
            }
        }

        const auto dynamic_block = block.dynamic();
        for (const auto* source : source_terms) {
            source->end_of_ray_source(dynamic_block, rayidx, wavel_threadidx,
                                      threadidx, radiance);
        }

        if (!have_to_integrate) {
            return;
        }

        if constexpr (N == 1) {
            if (m_wavelength_block_capacity == 1) {
                integrate_ray(radiance, source_terms, (*m_traced_rays)[rayidx],
                              m_traced_ray_od_matrix[rayidx],
                              m_scalar_shell_od[rayidx], block, rayidx,
                              wavel_threadidx, threadidx);
                return;
            }
        }
        integrate_ray(radiance, source_terms, (*m_traced_rays)[rayidx],
                      m_traced_ray_od_matrix[rayidx], m_shell_od[rayidx], block,
                      rayidx, wavel_threadidx, threadidx);
    }

    template <int NSTOKES>
    template <int N, typename ShellODMatrix>
    void SourceIntegrator<NSTOKES>::integrate_ray(
        sasktran2::WavelengthBlockDual<NSTOKES>& radiance,
        const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
        const sasktran2::raytracing::TracedRay& ray,
        const Eigen::SparseMatrix<double, Eigen::RowMajor>& od_matrix,
        const ShellODMatrix& shell_od,
        const sasktran2::WavelengthBlock<N>& batch, int rayidx,
        int wavel_threadidx, int threadidx) const {
        if (threadidx < 0 || threadidx >= m_thread_attenuation.size()) {
            throw std::invalid_argument(
                "Source-integrator thread scratch is not initialized");
        }
        if constexpr (N == 1) {
            for (int layeridx = 0; layeridx < ray.layers.size(); ++layeridx) {
                const double* od = &shell_od(layeridx, batch.start);
                const double attenuation = std::exp(-*od);
                sasktran2::WavelengthBlockODView local_shell_od(
                    od, &attenuation, 1, od_matrix, layeridx);

                auto radiance_value = radiance.value.col(0);
                if (m_calculate_derivatives) {
                    for (auto it = local_shell_od.derivative_iterator(); it;
                         ++it) {
                        radiance.derivative(it.index(), batch).col(0) -=
                            it.value() * radiance_value;
                    }
                }

                radiance_value *= attenuation;
                if (m_calculate_derivatives) {
                    if (m_use_sparse_derivative_tracking) {
                        for (const auto& [derivative_start, derivative_count] :
                             m_attenuation_active_derivative_ranges[rayidx]
                                                                   [layeridx]) {
                            radiance.derivative_range(derivative_start,
                                                      derivative_count,
                                                      batch) *= attenuation;
                        }
                    } else {
                        radiance.deriv.col(0) *= attenuation;
                    }
                }

                const auto& layer = ray.layers[layeridx];
                const auto entrance_weights = ray.entrance_weights(layeridx);
                const auto exit_weights = ray.exit_weights(layeridx);
                for (const auto* source : source_terms) {
                    if (source->has_interior_source()) {
                        source->dispatch_integrated_source(
                            batch, rayidx, layeridx, wavel_threadidx, threadidx,
                            layer, entrance_weights, exit_weights,
                            local_shell_od, radiance,
                            SourceTermInterface<
                                NSTOKES>::IntegrationDirection::backward);
                    }
                }
            }
            return;
        }

        auto& attenuation_storage = m_thread_attenuation[threadidx];
        if (attenuation_storage.size() < batch.count) {
            attenuation_storage.resize(batch.count);
        }
        auto attenuation = wavelength_head(attenuation_storage, batch);
        for (int layeridx = 0; layeridx < ray.layers.size(); ++layeridx) {
            const auto od = [&]() {
                if constexpr (N == Eigen::Dynamic) {
                    return shell_od.block(layeridx, batch.start, 1,
                                          batch.count);
                } else {
                    return shell_od.template block<1, N>(layeridx, batch.start);
                }
            }();
            attenuation.array() = (-od.array()).exp();
            sasktran2::WavelengthBlockODView local_shell_od(
                od.data(), attenuation.data(), batch.count, od_matrix,
                layeridx);

            auto radiance_value = wavelength_left_cols(radiance.value, batch);

            if (m_calculate_derivatives) {
                for (auto it = local_shell_od.derivative_iterator(); it; ++it) {
                    radiance.derivative(it.index(), batch).array() -=
                        it.value() * radiance_value.array();
                }
            }

            radiance_value.array().rowwise() *= attenuation.array();
            if (m_calculate_derivatives) {
                if (m_use_sparse_derivative_tracking) {
                    for (const auto& [derivative_start, derivative_count] :
                         m_attenuation_active_derivative_ranges[rayidx]
                                                               [layeridx]) {
                        radiance
                            .derivative_range(derivative_start,
                                              derivative_count, batch)
                            .array()
                            .rowwise() *= attenuation.array();
                    }
                } else {
                    wavelength_left_cols(radiance.deriv, batch)
                        .array()
                        .rowwise() *= attenuation.array();
                }
            }

            const auto& layer = ray.layers[layeridx];
            const auto entrance_weights = ray.entrance_weights(layeridx);
            const auto exit_weights = ray.exit_weights(layeridx);
            for (const auto* source : source_terms) {
                if (source->has_interior_source()) {
                    source->dispatch_integrated_source(
                        batch, rayidx, layeridx, wavel_threadidx, threadidx,
                        layer, entrance_weights, exit_weights, local_shell_od,
                        radiance,
                        SourceTermInterface<
                            NSTOKES>::IntegrationDirection::backward);
                }
            }
        }
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::integrate_optical_depth(
        Eigen::MatrixXd& optical_depth) {
        const int num_rays = static_cast<int>(m_traced_ray_od_matrix.size());
        for (int i = 0; i < num_rays; ++i) {
            if (m_wavelength_block_capacity == 1) {
                optical_depth.col(i) = m_scalar_shell_od[i].colwise().sum();
            } else {
                optical_depth.col(i) = m_shell_od[i].colwise().sum();
            }
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
        sasktran2::WavelengthBlockDual<NSTOKES> layer_source;
        layer_source.resize(1, 0);

        double current_od = 0;
        for (int j = (int)ray.layers.size() - 1; j >= 0; --j) {
            const sasktran2::raytracing::SphericalLayer& layer = ray.layers[j];
            const auto entrance_weights = ray.entrance_weights(j);
            const auto exit_weights = ray.exit_weights(j);

            const double shell_od = m_wavelength_block_capacity == 1
                                        ? m_scalar_shell_od[rayidx](j, wavelidx)
                                        : m_shell_od[rayidx](j, wavelidx);
            const double shell_attenuation = std::exp(-shell_od);
            const sasktran2::WavelengthBlockODView block_shell_od(
                &shell_od, &shell_attenuation, 1,
                m_traced_ray_od_matrix[rayidx], j);
            const auto& layer_interpolator = interpolator.interior_weights[j];
            // Calculate and add the layer source to the radiance
            double atten_factor = std::exp(-current_od);

            // Calculate all of the layer sources
            layer_source.set_zero(1);
            const sasktran2::WavelengthBlock<> block{wavelidx, 1};
            for (const auto& source : source_terms) {
                source->integrated_source(
                    block, rayidx, j, wavel_threadidx, threadidx, layer,
                    entrance_weights, exit_weights, block_shell_od,
                    layer_source,
                    SourceTermInterface<
                        NSTOKES>::IntegrationDirection::forward);
            }

            radiance.value += layer_source.value.col(0) * atten_factor;

            // Assign the accumulation weights
            double omega = 0;
            for (int i = 0; i < layer_interpolator.first.size(); ++i) {
                auto& index_weight = layer_interpolator.first[i];
                omega +=
                    m_atmosphere->storage().ssa(index_weight.first, wavelidx) *
                    index_weight.second;
            }
            double source_factor =
                omega * (1 - shell_attenuation) * atten_factor;

            for (const auto& ele : layer_interpolator.second) {
                for (int s = 0; s < NSTOKES; ++s) {
                    accumulation_values(std::get<2>(ele)[s]) +=
                        std::get<1>(ele) * source_factor;
                }
            }

            current_od += shell_od;
        }

        // Add source at the end of the ray
        layer_source.set_zero(1);
        const sasktran2::WavelengthBlock<> end_block{wavelidx, 1};
        for (const auto& source : source_terms) {
            source->end_of_ray_source(end_block, rayidx, wavel_threadidx,
                                      threadidx, layer_source);
        }

        radiance.value += layer_source.value.col(0) * std::exp(-1 * current_od);

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
