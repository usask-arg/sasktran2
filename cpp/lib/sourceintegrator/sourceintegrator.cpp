
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
        m_use_compact_geometry = false;
        m_compact_rays.clear();
        m_compact_max_layers = 0;
        m_use_sparse_derivative_tracking = false;
        m_attenuation_active_derivative_ranges.clear();

        // Construct the optical depth matrices.
        // This is the matrix so that matrix @ extinction = layer od, one matrix
        // for each ray Calculating this matrix beforehand makes calculating
        // derivatives easier, and removes excess computation for every
        // wavelength
        if (m_on_demand_optical_depth) {
            m_traced_ray_od_matrix.clear();
            std::size_t max_layers = 0;
            for (const auto& ray : traced_rays) {
                max_layers = std::max(max_layers, ray.layers.size());
            }
            m_empty_od_matrix.resize(static_cast<int>(max_layers),
                                     geometry.size());
            m_empty_od_matrix.setZero();
            m_empty_od_matrix.makeCompressed();
        } else {
            m_traced_ray_od_matrix.resize(traced_rays.size());
            for (int i = 0; i < traced_rays.size(); ++i) {
                sasktran2::raytracing::construct_od_matrix(
                    traced_rays[i], geometry, m_traced_ray_od_matrix[i]);
            }
        }

        m_traced_rays = &traced_rays;
        m_num_geometry_locations = geometry.size();
        m_num_geometry_dimensions = geometry.num_atmosphere_dimensions();
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::begin_compact_geometry_2d(
        std::vector<sasktran2::raytracing::TracedRay>& traced_rays,
        const Geometry& geometry) {
        if (!m_on_demand_optical_depth ||
            geometry.num_atmosphere_dimensions() != 2) {
            throw std::logic_error(
                "Incremental compact geometry requires on-demand 2D source "
                "integration");
        }
        m_use_compact_geometry = false;
        m_compact_rays.clear();
        m_compact_rays.resize(traced_rays.size());
        m_compact_max_layers = 0;
        m_traced_ray_od_matrix.clear();
        m_empty_od_matrix.resize(0, geometry.size());
        m_traced_rays = &traced_rays;
        m_num_geometry_locations = geometry.size();
        m_num_geometry_dimensions = geometry.num_atmosphere_dimensions();
        m_use_sparse_derivative_tracking = false;
        m_attenuation_active_derivative_ranges.clear();
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::compact_geometry_2d_ray(
        int rayidx, sasktran2::raytracing::TracedRay& ray) {
        if (!ray.is_straight) {
            throw std::invalid_argument(
                "Compact 2D source integration requires straight rays");
        }
        auto& compact = m_compact_rays.at(rayidx);
        compact.layers.clear();
        compact.layers.reserve(ray.layers.size());

        std::vector<int> first_indices;
        std::vector<double> first_entrance;
        std::vector<double> first_exit;
        std::vector<double> first_optical_depth;
        if (!ray.layers.empty()) {
            const auto entrance = ray.entrance_weights(0);
            const auto exit = ray.exit_weights(0);
            const auto optical_depth = ray.optical_depth_weights(0);
            first_indices.assign(entrance.indices_data(),
                                 entrance.indices_data() + entrance.size());
            first_entrance.assign(entrance.weights_data(),
                                  entrance.weights_data() + entrance.size());
            first_exit.assign(exit.weights_data(),
                              exit.weights_data() + exit.size());
            first_optical_depth.assign(optical_depth.weights_data(),
                                       optical_depth.weights_data() +
                                           optical_depth.size());
        }

        for (std::size_t layeridx = 0; layeridx < ray.layers.size();
             ++layeridx) {
            const auto& layer = ray.layers[layeridx];
            if (layer.grid_weight_count >
                std::numeric_limits<std::uint8_t>::max()) {
                throw std::length_error(
                    "Compact 2D layer stencil exceeds 8-bit count");
            }
            compact.layers.push_back(
                {layer.layer_distance, layer.od_quad_start, layer.od_quad_end,
                 layer.od_quad_start_fraction, layer.od_quad_end_fraction,
                 layer.grid_weight_offset, layer.grid_weight_count});
        }

        ray.move_grid_weights_to(compact.weights.indices,
                                 compact.weights.entrance, compact.weights.exit,
                                 compact.weights.optical_depth);

        sasktran2::raytracing::TracedRay retained;
        retained.observer_and_look = ray.observer_and_look;
        retained.is_straight = ray.is_straight;
        retained.ground_is_hit = ray.ground_is_hit;
        retained.tangent_radius = ray.tangent_radius;
        if (!ray.layers.empty()) {
            retained.layers.push_back(ray.layers.front());
            retained.set_layer_weights(0, first_indices.data(),
                                       first_entrance.data(), first_exit.data(),
                                       first_optical_depth.data(),
                                       first_indices.size());
        }
        ray = std::move(retained);
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::finalize_compact_geometry_2d() {
        m_compact_max_layers = 0;
        for (const auto& ray : m_compact_rays) {
            m_compact_max_layers =
                std::max(m_compact_max_layers, ray.layers.size());
        }
        m_empty_od_matrix.resize(static_cast<int>(m_compact_max_layers),
                                 m_num_geometry_locations);
        m_empty_od_matrix.setZero();
        m_empty_od_matrix.makeCompressed();
        m_use_compact_geometry = true;
    }

    template <int NSTOKES>
    int SourceIntegrator<NSTOKES>::integration_num_layers(int rayidx) const {
        return m_use_compact_geometry
                   ? static_cast<int>(m_compact_rays[rayidx].layers.size())
                   : static_cast<int>((*m_traced_rays)[rayidx].layers.size());
    }

    template <int NSTOKES>
    const sasktran2::raytracing::TracedLayer&
    SourceIntegrator<NSTOKES>::integration_layer(
        int rayidx, int layeridx,
        sasktran2::raytracing::TracedLayer& scratch) const {
        if (!m_use_compact_geometry) {
            return (*m_traced_rays)[rayidx].layers[layeridx];
        }
        const auto& compact = m_compact_rays[rayidx].layers[layeridx];
        scratch = {};
        scratch.layer_distance = compact.layer_distance;
        scratch.curvature_factor = 1.0;
        scratch.od_quad_start = compact.od_quad_start;
        scratch.od_quad_end = compact.od_quad_end;
        scratch.od_quad_start_fraction = compact.od_quad_start_fraction;
        scratch.od_quad_end_fraction = compact.od_quad_end_fraction;
        const auto& retained_ray = (*m_traced_rays)[rayidx];
        if (!retained_ray.layers.empty()) {
            scratch.average_look_away =
                retained_ray.layers.front().average_look_away;
        }
        return scratch;
    }

    template <int NSTOKES>
    sasktran2::raytracing::GridWeightStencilView
    SourceIntegrator<NSTOKES>::integration_entrance_weights(
        int rayidx, int layeridx) const {
        if (!m_use_compact_geometry) {
            return (*m_traced_rays)[rayidx].entrance_weights(layeridx);
        }
        const auto& layer = m_compact_rays[rayidx].layers[layeridx];
        const auto& weights = m_compact_rays[rayidx].weights;
        return {weights.indices.data() + layer.grid_weight_offset,
                weights.entrance.data() + layer.grid_weight_offset,
                layer.grid_weight_count};
    }

    template <int NSTOKES>
    sasktran2::raytracing::GridWeightStencilView
    SourceIntegrator<NSTOKES>::integration_exit_weights(int rayidx,
                                                        int layeridx) const {
        if (!m_use_compact_geometry) {
            return (*m_traced_rays)[rayidx].exit_weights(layeridx);
        }
        const auto& layer = m_compact_rays[rayidx].layers[layeridx];
        const auto& weights = m_compact_rays[rayidx].weights;
        return {weights.indices.data() + layer.grid_weight_offset,
                weights.exit.data() + layer.grid_weight_offset,
                layer.grid_weight_count};
    }

    template <int NSTOKES>
    sasktran2::raytracing::GridWeightStencilView
    SourceIntegrator<NSTOKES>::integration_optical_depth_weights(
        int rayidx, int layeridx) const {
        if (!m_use_compact_geometry) {
            return (*m_traced_rays)[rayidx].optical_depth_weights(layeridx);
        }
        const auto& layer = m_compact_rays[rayidx].layers[layeridx];
        const auto& weights = m_compact_rays[rayidx].weights;
        return {weights.indices.data() + layer.grid_weight_offset,
                weights.optical_depth.data() + layer.grid_weight_offset,
                layer.grid_weight_count};
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
        m_atmosphere = &atmo;
        if (m_on_demand_optical_depth) {
            m_scalar_shell_od.clear();
            m_shell_od.clear();
            m_calculate_derivatives = false;
            return;
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
    void SourceIntegrator<NSTOKES>::integrate_jvp(
        sasktran2::RadianceJVP<NSTOKES>& radiance,
        const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
        int wavelength, int rayidx, int wavel_threadidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) const {
        if (m_wavelength_block_capacity != 1 || m_scalar_shell_od.empty()) {
            throw std::logic_error(
                "Native JVP requires scalar source-integrator storage");
        }
        radiance.set_zero();
        for (const auto* source : source_terms) {
            source->end_of_ray_source_jvp(wavelength, rayidx, wavel_threadidx,
                                          threadidx, native_tangent, radiance);
        }

        const auto& ray = (*m_traced_rays)[rayidx];
        const auto& od_matrix = m_traced_ray_od_matrix[rayidx];
        for (int layeridx = 0; layeridx < ray.layers.size(); ++layeridx) {
            const double od = m_scalar_shell_od[rayidx](layeridx, wavelength);
            const double attenuation = std::exp(-od);
            const sasktran2::WavelengthBlockODView shell_od(
                &od, &attenuation, 1, od_matrix, layeridx);
            double od_jvp = 0.0;
            for (auto derivative = shell_od.derivative_iterator(); derivative;
                 ++derivative) {
                od_jvp +=
                    derivative.value() * native_tangent(derivative.index());
            }
            radiance.jvp =
                attenuation * (radiance.jvp - od_jvp * radiance.value);
            radiance.value *= attenuation;

            const auto& layer = ray.layers[layeridx];
            const auto entrance_weights = ray.entrance_weights(layeridx);
            const auto exit_weights = ray.exit_weights(layeridx);
            for (const auto* source : source_terms) {
                if (source->has_interior_source()) {
                    source->integrated_source_jvp(
                        wavelength, rayidx, layeridx, wavel_threadidx,
                        threadidx, layer, entrance_weights, exit_weights,
                        shell_od, native_tangent, radiance);
                }
            }
        }
        for (const auto* source : source_terms) {
            source->start_of_ray_source_jvp(wavelength, rayidx, wavel_threadidx,
                                            threadidx, native_tangent,
                                            radiance);
        }
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::integrate_vjp(
        Eigen::Vector<double, NSTOKES>& radiance,
        const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
        int wavelength, int rayidx, int wavel_threadidx, int threadidx,
        const Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        if (m_wavelength_block_capacity != 1 || m_scalar_shell_od.empty()) {
            throw std::logic_error(
                "Native VJP requires scalar source-integrator storage");
        }
        const auto& ray = (*m_traced_rays)[rayidx];
        const auto& od_matrix = m_traced_ray_od_matrix[rayidx];
        const sasktran2::WavelengthBlock<1> block{wavelength, 1};
        sasktran2::WavelengthBlockDual<NSTOKES> state;
        state.resize(1, 0, true);
        for (const auto* source : source_terms) {
            source->end_of_ray_source(block.dynamic(), rayidx, wavel_threadidx,
                                      threadidx, state);
        }

        std::vector<Eigen::Vector<double, NSTOKES>> values_before_layer(
            ray.layers.size());
        std::vector<double> attenuation(ray.layers.size());
        for (int layeridx = 0; layeridx < ray.layers.size(); ++layeridx) {
            values_before_layer[layeridx] = state.value.col(0);
            const double od = m_scalar_shell_od[rayidx](layeridx, wavelength);
            attenuation[layeridx] = std::exp(-od);
            const sasktran2::WavelengthBlockODView shell_od(
                &od, &attenuation[layeridx], 1, od_matrix, layeridx);
            state.value.col(0) *= attenuation[layeridx];

            const auto& layer = ray.layers[layeridx];
            const auto entrance_weights = ray.entrance_weights(layeridx);
            const auto exit_weights = ray.exit_weights(layeridx);
            for (const auto* source : source_terms) {
                if (source->has_interior_source()) {
                    source->dispatch_integrated_source(
                        block, rayidx, layeridx, wavel_threadidx, threadidx,
                        layer, entrance_weights, exit_weights, shell_od, state,
                        SourceTermInterface<
                            NSTOKES>::IntegrationDirection::backward);
                }
            }
        }

        std::vector<Eigen::Vector<double, NSTOKES>> values_before_start;
        values_before_start.reserve(source_terms.size());
        for (const auto* source : source_terms) {
            values_before_start.push_back(state.value.col(0));
            source->start_of_ray_source(block.dynamic(), rayidx,
                                        wavel_threadidx, threadidx, state);
        }
        radiance = state.value.col(0);

        Eigen::Vector<double, NSTOKES> state_cotangent = cotangent;
        for (int source_index = static_cast<int>(source_terms.size()) - 1;
             source_index >= 0; --source_index) {
            source_terms[source_index]->start_of_ray_source_vjp(
                wavelength, rayidx, wavel_threadidx, threadidx,
                values_before_start[source_index], state_cotangent,
                native_gradient);
        }
        for (int layeridx = static_cast<int>(ray.layers.size()) - 1;
             layeridx >= 0; --layeridx) {
            const double od = m_scalar_shell_od[rayidx](layeridx, wavelength);
            const sasktran2::WavelengthBlockODView shell_od(
                &od, &attenuation[layeridx], 1, od_matrix, layeridx);
            const auto& layer = ray.layers[layeridx];
            const auto entrance_weights = ray.entrance_weights(layeridx);
            const auto exit_weights = ray.exit_weights(layeridx);
            for (const auto* source : source_terms) {
                if (source->has_interior_source()) {
                    source->integrated_source_vjp(
                        wavelength, rayidx, layeridx, wavel_threadidx,
                        threadidx, layer, entrance_weights, exit_weights,
                        shell_od, state_cotangent, native_gradient);
                }
            }

            const double attenuation_cotangent =
                state_cotangent.dot(values_before_layer[layeridx]);
            const double od_cotangent =
                -attenuation[layeridx] * attenuation_cotangent;
            for (auto derivative = shell_od.derivative_iterator(); derivative;
                 ++derivative) {
                native_gradient(derivative.index()) +=
                    derivative.value() * od_cotangent;
            }
            state_cotangent *= attenuation[layeridx];
        }

        for (const auto* source : source_terms) {
            source->end_of_ray_source_vjp(wavelength, rayidx, wavel_threadidx,
                                          threadidx, state_cotangent,
                                          native_gradient);
        }
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
        const int num_layers = integration_num_layers(rayidx);

        // If we don't have to calculate derivatives then it is faster to
        // iterate over the ray backwards, i.e., from the observer to the end of
        // the atmosphere
        sasktran2::WavelengthBlockDual<NSTOKES> layer_source;
        layer_source.resize(1, 0);

        double current_od = 0;
        for (int j = num_layers - 1; j >= 0; --j) {
            sasktran2::raytracing::TracedLayer layer_scratch;
            const auto& layer = integration_layer(rayidx, j, layer_scratch);
            const auto entrance_weights =
                integration_entrance_weights(rayidx, j);
            const auto exit_weights = integration_exit_weights(rayidx, j);

            double shell_od = 0.0;
            if (m_on_demand_optical_depth) {
                const auto optical_depth_weights =
                    integration_optical_depth_weights(rayidx, j);
                for (std::size_t index = 0;
                     index < optical_depth_weights.size(); ++index) {
                    const auto index_weight = optical_depth_weights[index];
                    shell_od += index_weight.second *
                                m_atmosphere->storage().total_extinction(
                                    index_weight.first, wavelidx);
                }
            } else {
                shell_od = m_wavelength_block_capacity == 1
                               ? m_scalar_shell_od[rayidx](j, wavelidx)
                               : m_shell_od[rayidx](j, wavelidx);
            }
            const double shell_attenuation = std::exp(-shell_od);
            const auto& derivative_matrix =
                m_on_demand_optical_depth ? m_empty_od_matrix
                                          : m_traced_ray_od_matrix[rayidx];
            const sasktran2::WavelengthBlockODView block_shell_od(
                &shell_od, &shell_attenuation, 1, derivative_matrix, j);
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
            const auto atmosphere_nodes =
                integration_entrance_weights(rayidx, j);
            for (std::size_t i = layer_interpolator.atmosphere_offset;
                 i < layer_interpolator.atmosphere_offset +
                         layer_interpolator.atmosphere_count;
                 ++i) {
                const auto node =
                    atmosphere_nodes[i - layer_interpolator.atmosphere_offset];
                omega += m_atmosphere->storage().ssa(node.first, wavelidx) *
                         interpolator.atmosphere_weights[i];
            }
            double source_factor =
                omega * (1 - shell_attenuation) * atten_factor;

            for (std::size_t i = layer_interpolator.source_offset;
                 i < layer_interpolator.source_offset +
                         layer_interpolator.source_count;
                 ++i) {
                for (int s = 0; s < NSTOKES; ++s) {
                    accumulation_values(
                        interpolator.source_accumulation_index(i, s)) +=
                        interpolator.source_weights[i] * source_factor;
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
            for (std::size_t i = 0;
                 i < interpolator.ground_source_weights.size(); ++i) {
                for (int s = 0; s < NSTOKES; ++s) {
                    accumulation_values(interpolator.ground_accumulation_index(
                        i, s)) += interpolator.ground_source_weights[i] *
                                  std::exp(-1 * current_od);
                }
            }
        }
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::integrate_and_emplace_accumulation_jvp(
        const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
        int wavelidx, int rayidx, int wavel_threadidx, int threadidx,
        const SInterpolator& source_interpolator,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        Eigen::Ref<Eigen::VectorXd> accumulation_value_tangent,
        Eigen::Ref<Eigen::VectorXd> first_order_forcing_tangent) const {
        const auto& ray = (*m_traced_rays)[rayidx];
        const auto& interpolator = source_interpolator[rayidx];
        double current_od = 0.0;
        double current_od_jvp = 0.0;
        const int num_layers = integration_num_layers(rayidx);

        for (int layeridx = num_layers - 1; layeridx >= 0; --layeridx) {
            const auto od_weights =
                integration_optical_depth_weights(rayidx, layeridx);
            double shell_od = 0.0;
            double shell_od_jvp = 0.0;
            for (std::size_t index = 0; index < od_weights.size(); ++index) {
                const auto weight = od_weights[index];
                shell_od +=
                    weight.second * m_atmosphere->storage().total_extinction(
                                        weight.first, wavelidx);
                shell_od_jvp += weight.second * native_tangent(weight.first);
            }
            const double shell_attenuation = std::exp(-shell_od);
            const double current_attenuation = std::exp(-current_od);
            const sasktran2::WavelengthBlockODView block_shell_od(
                &shell_od, &shell_attenuation, 1, od_weights.indices_data(),
                od_weights.weights_data(), od_weights.size());

            sasktran2::RadianceJVP<NSTOKES> layer_source;
            layer_source.set_zero();
            sasktran2::raytracing::TracedLayer layer_scratch;
            const auto& layer =
                integration_layer(rayidx, layeridx, layer_scratch);
            const auto entrance_weights =
                integration_entrance_weights(rayidx, layeridx);
            const auto exit_weights =
                integration_exit_weights(rayidx, layeridx);
            for (const auto* source : source_terms) {
                if (source->has_interior_source()) {
                    source->integrated_source_jvp(
                        wavelidx, rayidx, layeridx, wavel_threadidx, threadidx,
                        layer, entrance_weights, exit_weights, block_shell_od,
                        native_tangent, layer_source);
                }
            }

            for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                first_order_forcing_tangent(rayidx * NSTOKES + stokes) +=
                    current_attenuation *
                    (layer_source.jvp(stokes) -
                     current_od_jvp * layer_source.value(stokes));
            }

            const auto& layer_interpolator =
                interpolator.interior_weights[layeridx];
            const auto atmosphere_nodes =
                integration_entrance_weights(rayidx, layeridx);
            double omega = 0.0;
            double omega_jvp = 0.0;
            for (std::size_t index = layer_interpolator.atmosphere_offset;
                 index < layer_interpolator.atmosphere_offset +
                             layer_interpolator.atmosphere_count;
                 ++index) {
                const int atmosphere_index =
                    atmosphere_nodes[index -
                                     layer_interpolator.atmosphere_offset]
                        .first;
                const double weight = interpolator.atmosphere_weights[index];
                omega +=
                    m_atmosphere->storage().ssa(atmosphere_index, wavelidx) *
                    weight;
                omega_jvp +=
                    native_tangent(m_atmosphere->ssa_deriv_start_index() +
                                   atmosphere_index) *
                    weight;
            }
            const double source_factor_jvp =
                current_attenuation *
                (omega_jvp * (1.0 - shell_attenuation) +
                 omega * shell_attenuation * shell_od_jvp -
                 omega * (1.0 - shell_attenuation) * current_od_jvp);
            for (std::size_t index = layer_interpolator.source_offset;
                 index < layer_interpolator.source_offset +
                             layer_interpolator.source_count;
                 ++index) {
                for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                    accumulation_value_tangent(
                        interpolator.source_accumulation_index(index,
                                                               stokes)) +=
                        interpolator.source_weights[index] * source_factor_jvp;
                }
            }
            current_od += shell_od;
            current_od_jvp += shell_od_jvp;
        }

        sasktran2::RadianceJVP<NSTOKES> end_source;
        end_source.set_zero();
        for (const auto* source : source_terms) {
            source->end_of_ray_source_jvp(wavelidx, rayidx, wavel_threadidx,
                                          threadidx, native_tangent,
                                          end_source);
        }
        const double end_attenuation = std::exp(-current_od);
        for (int stokes = 0; stokes < NSTOKES; ++stokes) {
            first_order_forcing_tangent(rayidx * NSTOKES + stokes) +=
                end_attenuation * (end_source.jvp(stokes) -
                                   current_od_jvp * end_source.value(stokes));
        }

        if (ray.ground_is_hit) {
            const double ground_factor_jvp = -end_attenuation * current_od_jvp;
            for (std::size_t index = 0;
                 index < interpolator.ground_source_weights.size(); ++index) {
                for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                    accumulation_value_tangent(
                        interpolator.ground_accumulation_index(index,
                                                               stokes)) +=
                        interpolator.ground_source_weights[index] *
                        ground_factor_jvp;
                }
            }
        }
    }

    template <int NSTOKES>
    void SourceIntegrator<NSTOKES>::accumulate_accumulation_vjp(
        const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
        int wavelidx, int rayidx, int wavel_threadidx, int threadidx,
        const SInterpolator& source_interpolator,
        Eigen::Ref<const Eigen::VectorXd> accumulation_value_gradient,
        Eigen::Ref<const Eigen::VectorXd> transport_solution,
        Eigen::Ref<const Eigen::VectorXi> transport_column_indices,
        Eigen::Ref<const Eigen::VectorXd> first_order_forcing_gradient,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        const auto& ray = (*m_traced_rays)[rayidx];
        const auto& interpolator = source_interpolator[rayidx];
        const Eigen::Vector<double, NSTOKES> forcing_cotangent =
            first_order_forcing_gradient.template segment<NSTOKES>(rayidx *
                                                                   NSTOKES);
        const int num_layers = integration_num_layers(rayidx);
        std::vector<double> shell_od_cotangent(num_layers, 0.0);
        std::vector<double> current_od_cotangent(num_layers, 0.0);
        double current_od = 0.0;

        for (int traversal = 0; traversal < num_layers; ++traversal) {
            const int layeridx = num_layers - 1 - traversal;
            const auto od_weights =
                integration_optical_depth_weights(rayidx, layeridx);
            double shell_od = 0.0;
            for (std::size_t index = 0; index < od_weights.size(); ++index) {
                const auto weight = od_weights[index];
                shell_od +=
                    weight.second * m_atmosphere->storage().total_extinction(
                                        weight.first, wavelidx);
            }
            const double shell_attenuation = std::exp(-shell_od);
            const double current_attenuation = std::exp(-current_od);
            const sasktran2::WavelengthBlockODView block_shell_od(
                &shell_od, &shell_attenuation, 1, od_weights.indices_data(),
                od_weights.weights_data(), od_weights.size());

            sasktran2::WavelengthBlockDual<NSTOKES> layer_source;
            layer_source.resize(1, 0, true);
            const sasktran2::WavelengthBlock<> block{wavelidx, 1};
            sasktran2::raytracing::TracedLayer layer_scratch;
            const auto& layer =
                integration_layer(rayidx, layeridx, layer_scratch);
            const auto entrance_weights =
                integration_entrance_weights(rayidx, layeridx);
            const auto exit_weights =
                integration_exit_weights(rayidx, layeridx);
            for (const auto* source : source_terms) {
                if (source->has_interior_source()) {
                    source->integrated_source(
                        block, rayidx, layeridx, wavel_threadidx, threadidx,
                        layer, entrance_weights, exit_weights, block_shell_od,
                        layer_source,
                        SourceTermInterface<
                            NSTOKES>::IntegrationDirection::forward);
                    source->integrated_source_vjp(
                        wavelidx, rayidx, layeridx, wavel_threadidx, threadidx,
                        layer, entrance_weights, exit_weights, block_shell_od,
                        current_attenuation * forcing_cotangent,
                        native_gradient);
                }
            }
            current_od_cotangent[traversal] -=
                current_attenuation *
                forcing_cotangent.dot(layer_source.value.col(0));

            const auto& layer_interpolator =
                interpolator.interior_weights[layeridx];
            const auto atmosphere_nodes =
                integration_entrance_weights(rayidx, layeridx);
            double omega = 0.0;
            for (std::size_t index = layer_interpolator.atmosphere_offset;
                 index < layer_interpolator.atmosphere_offset +
                             layer_interpolator.atmosphere_count;
                 ++index) {
                omega +=
                    interpolator.atmosphere_weights[index] *
                    m_atmosphere->storage().ssa(
                        atmosphere_nodes[index -
                                         layer_interpolator.atmosphere_offset]
                            .first,
                        wavelidx);
            }
            double source_factor_cotangent = 0.0;
            for (std::size_t index = layer_interpolator.source_offset;
                 index < layer_interpolator.source_offset +
                             layer_interpolator.source_count;
                 ++index) {
                for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                    source_factor_cotangent +=
                        interpolator.source_weights[index] *
                        (accumulation_value_gradient.size() == 0
                             ? forcing_cotangent(stokes) *
                                   transport_solution(transport_column_indices(
                                       interpolator.source_accumulation_index(
                                           index, stokes)))
                             : accumulation_value_gradient(
                                   interpolator.source_accumulation_index(
                                       index, stokes)));
                }
            }
            const double omega_cotangent = source_factor_cotangent *
                                           (1.0 - shell_attenuation) *
                                           current_attenuation;
            for (std::size_t index = layer_interpolator.atmosphere_offset;
                 index < layer_interpolator.atmosphere_offset +
                             layer_interpolator.atmosphere_count;
                 ++index) {
                const int atmosphere_index =
                    atmosphere_nodes[index -
                                     layer_interpolator.atmosphere_offset]
                        .first;
                native_gradient(m_atmosphere->ssa_deriv_start_index() +
                                atmosphere_index) +=
                    interpolator.atmosphere_weights[index] * omega_cotangent;
            }
            shell_od_cotangent[traversal] += source_factor_cotangent * omega *
                                             shell_attenuation *
                                             current_attenuation;
            current_od_cotangent[traversal] -= source_factor_cotangent * omega *
                                               (1.0 - shell_attenuation) *
                                               current_attenuation;
            current_od += shell_od;
        }

        sasktran2::WavelengthBlockDual<NSTOKES> end_source;
        end_source.resize(1, 0, true);
        const sasktran2::WavelengthBlock<> end_block{wavelidx, 1};
        for (const auto* source : source_terms) {
            source->end_of_ray_source(end_block, rayidx, wavel_threadidx,
                                      threadidx, end_source);
            source->end_of_ray_source_vjp(
                wavelidx, rayidx, wavel_threadidx, threadidx,
                std::exp(-current_od) * forcing_cotangent, native_gradient);
        }
        double trailing_current_od_cotangent =
            -std::exp(-current_od) *
            forcing_cotangent.dot(end_source.value.col(0));

        if (ray.ground_is_hit) {
            double ground_factor_cotangent = 0.0;
            for (std::size_t index = 0;
                 index < interpolator.ground_source_weights.size(); ++index) {
                for (int stokes = 0; stokes < NSTOKES; ++stokes) {
                    ground_factor_cotangent +=
                        interpolator.ground_source_weights[index] *
                        (accumulation_value_gradient.size() == 0
                             ? forcing_cotangent(stokes) *
                                   transport_solution(transport_column_indices(
                                       interpolator.ground_accumulation_index(
                                           index, stokes)))
                             : accumulation_value_gradient(
                                   interpolator.ground_accumulation_index(
                                       index, stokes)));
                }
            }
            trailing_current_od_cotangent -=
                ground_factor_cotangent * std::exp(-current_od);
        }

        for (int traversal = num_layers - 1; traversal >= 0; --traversal) {
            shell_od_cotangent[traversal] += trailing_current_od_cotangent;
            trailing_current_od_cotangent += current_od_cotangent[traversal];
            const int layeridx = num_layers - 1 - traversal;
            const auto od_weights =
                integration_optical_depth_weights(rayidx, layeridx);
            for (std::size_t index = 0; index < od_weights.size(); ++index) {
                const auto weight = od_weights[index];
                native_gradient(weight.first) +=
                    weight.second * shell_od_cotangent[traversal];
            }
        }
    }

    template class SourceIntegrator<1>;
    template class SourceIntegrator<3>;
} // namespace sasktran2
