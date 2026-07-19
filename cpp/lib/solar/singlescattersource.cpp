#include "sasktran2/grids.h"
#include "sasktran2/raytracing.h"
#include "sasktran2/source_algorithms.h"
#include "sasktran2/source_interface.h"
#include <sasktran2/solartransmission.h>
#include <sasktran2/dual.h>
#include <sasktran2/config.h>
#include <sasktran2/atmosphere/atmosphere.h>
#include <sasktran2/math/trig.h>
#include <cmath>

namespace {
    constexpr double DENSE_GEOMETRY_THRESHOLD = 0.25;
}

namespace sasktran2::solartransmission {
    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        // Store the atmosphere for later
        m_atmosphere = &atmosphere;
        this->m_phase_handler.initialize_atmosphere(atmosphere);

        if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
            if (atmosphere.num_deriv() > 0) {
                initialize_active_derivative_indices();
            } else {
                m_active_derivative_indices.clear();
            }

            initialize_wavelength_blocks(m_wavelength_batch_capacity);
        }

        // Initialize some local memory storage
        for (int i = 0; i < m_start_source_cache.size(); ++i) {
            m_start_source_cache[i].resize(NSTOKES, atmosphere.num_deriv(),
                                           false);
            m_end_source_cache[i].resize(NSTOKES, atmosphere.num_deriv(),
                                         false);
        }
    };

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::initialize_config(
        const sasktran2::Config& config) {
        m_config = &config;

        this->m_solar_transmission.initialize_config(config);
        this->m_phase_handler.initialize_config(config);

        // Set up storage for each thread
        // m_solar_trans.resize(config.num_threads());
        m_solar_trans.resize(config.num_wavelength_threads());

        m_start_source_cache.resize(config.num_threads());
        m_end_source_cache.resize(config.num_threads());
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::calculate_single(int wavelidx,
                                                           int threadidx) {
        ZoneScopedN("Single Scatter Source Calculation");
        // Don't have to do anything here
        m_phase_handler.calculate(wavelidx, threadidx);

        // Calculate the solar transmission at each cell
        if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
            // Faster to use the dense matrix if most of the elements are
            // nonzero
            if (m_geometry_matrix.size() > 0 &&
                double(m_geometry_sparse.nonZeros()) /
                        double(m_geometry_matrix.size()) >
                    DENSE_GEOMETRY_THRESHOLD) {
                m_solar_trans[threadidx].noalias() =
                    m_geometry_matrix *
                    m_atmosphere->storage().total_extinction(
                        Eigen::placeholders::all, wavelidx);
            } else {
                m_solar_trans[threadidx].noalias() =
                    m_geometry_sparse *
                    m_atmosphere->storage().total_extinction(
                        Eigen::placeholders::all, wavelidx);
            }
        }

        if constexpr (std::is_same_v<S, SolarTransmissionTable>) {
            m_solar_trans[threadidx].noalias() =
                m_geometry_sparse * (m_solar_transmission.geometry_matrix() *
                                     m_atmosphere->storage().total_extinction(
                                         Eigen::placeholders::all, wavelidx));
        }

        m_solar_trans[threadidx] =
            exp(-m_solar_trans[threadidx].array()) *
            m_atmosphere->storage().solar_irradiance(wavelidx);
        for (int i = 0; i < m_ground_hit_flag.size(); ++i) {
            if (m_ground_hit_flag[i]) {
                m_solar_trans[threadidx][i] = 0;
            }
        }
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::initialize_wavelength_blocks(
        int block_size) {
        if constexpr (!std::is_same_v<S, SolarTransmissionExact>) {
            throw std::logic_error(
                "Solar transmission tables do not support wavelength "
                "batching");
        } else {
            m_wavelength_batch_capacity = block_size;
            m_solar_trans_batch.resize(m_config->num_wavelength_threads());
            for (auto& storage : m_solar_trans_batch) {
                storage.resize(m_geometry_sparse.rows(), block_size);
            }
            m_batch_source_cache.resize(m_config->num_threads());
            m_batch_integration_cache.resize(m_config->num_threads());
            for (auto& thread_cache : m_batch_source_cache) {
                for (auto& endpoint_cache : thread_cache) {
                    endpoint_cache.resize(block_size);
                }
            }
            for (auto& thread_cache : m_batch_integration_cache) {
                thread_cache.resize(block_size);
            }
            m_phase_handler.initialize_wavelength_blocks(block_size);
        }
    }

    template <typename S, int NSTOKES>
    template <int N>
    void SingleScatterSource<S, NSTOKES>::calculate_block(
        const sasktran2::WavelengthBlock<N>& batch, int threadidx) {
        if constexpr (!std::is_same_v<S, SolarTransmissionExact>) {
            throw std::logic_error(
                "Solar transmission tables do not support wavelength "
                "batching");
        } else {
            ZoneScopedN("Single Scatter Source Batch Calculation");
            if (batch.count > m_wavelength_batch_capacity) {
                throw std::invalid_argument(
                    "Wavelength batch exceeds single scatter storage "
                    "capacity");
            }

            m_phase_handler.calculate_block(batch, threadidx);
            auto solar_trans =
                wavelength_left_cols(m_solar_trans_batch[threadidx], batch);
            const auto extinction = wavelength_middle_cols(
                m_atmosphere->storage().total_extinction, batch);
            if (m_geometry_matrix.size() > 0 &&
                double(m_geometry_sparse.nonZeros()) /
                        double(m_geometry_matrix.size()) >
                    DENSE_GEOMETRY_THRESHOLD) {
                solar_trans.noalias() = m_geometry_matrix * extinction;
            } else {
                solar_trans.noalias() = m_geometry_sparse * extinction;
            }
            solar_trans = (-solar_trans.array()).exp().matrix();
            for (int lane = 0; lane < batch.count; ++lane) {
                solar_trans.col(lane) *=
                    m_atmosphere->storage().solar_irradiance(
                        batch.wavelength(lane));
            }
            for (int row = 0; row < m_ground_hit_flag.size(); ++row) {
                if (m_ground_hit_flag[row]) {
                    solar_trans.row(row).setZero();
                }
            }
        }
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::end_of_ray_source_single(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        sasktran2::WavelengthBlockLaneDualView<NSTOKES, 1>& source) const {
        if (m_los_ground_is_hit.at(losidx)) {
            const auto& first_layer = m_los_end_layers.at(losidx);
            // Single scatter ground source is solar_trans * cos(th) * brdf

            // Cosine of direction to the sun at the surface
            // TODO: This does not account for refraction?
            double mu_in = first_layer.exit.cos_zenith_angle(
                m_geometry.coordinates().sun_unit());
            if (mu_in <= 0.0) {
                return;
            }

            // Cosine of direction to LOS at the surface
            double mu_out = -1.0 * first_layer.exit.cos_zenith_angle(
                                       first_layer.average_look_away);

            // We already have the azimuthal difference
            double phi_diff = first_layer.saz_exit;

            Eigen::Matrix<double, NSTOKES, NSTOKES> brdf =
                m_atmosphere->surface().brdf(wavelidx, mu_in, mu_out, phi_diff);

            int exit_index = m_index_map[losidx][0];

            double solar_trans = m_solar_trans[wavel_threadidx](exit_index);

            Eigen::Vector<double, NSTOKES> source_value =
                solar_trans * brdf(Eigen::placeholders::all, 0) * mu_in;

#ifdef SASKTRAN_DEBUG_ASSERTS
            if (source_value.hasNaN()) {
                spdlog::warn(
                    "NaN detected in single scatter ground source calculation");
                source_value.setZero();
            }
#endif

            source.value.array() += source_value.array();
            if (source.deriv.size() > 0) {
                // Add on the solar transmission derivative factors
                if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
                    if (m_config->wf_precision() !=
                        sasktran2::Config::WeightingFunctionPrecision::
                            limited) {
                        // Have to apply the solar transmission derivative
                        // factors
                        for (Eigen::SparseMatrix<double,
                                                 Eigen::RowMajor>::InnerIterator
                                 it(m_geometry_sparse, exit_index);
                             it; ++it) {
                            source.deriv(Eigen::placeholders::all,
                                         it.index()) -=
                                it.value() * source_value;
                        }
                    }
                }

                for (int k = 0; k < m_atmosphere->surface().num_deriv(); ++k) {
                    // And then the surface derivative factors
                    Eigen::Matrix<double, NSTOKES, NSTOKES> brdf_deriv =
                        m_atmosphere->surface().d_brdf(wavelidx, mu_in, mu_out,
                                                       phi_diff, k);

                    source.deriv(Eigen::placeholders::all,
                                 m_atmosphere->surface_deriv_start_index() +
                                     k) +=
                        solar_trans * mu_in *
                        brdf_deriv(Eigen::placeholders::all, 0);
                }
            }
        }
    }

    template <typename S, int NSTOKES>
    template <int N>
    void SingleScatterSource<S, NSTOKES>::end_of_ray_source_block(
        const sasktran2::WavelengthBlock<N>& batch, int losidx,
        int wavel_threadidx, int threadidx,
        sasktran2::WavelengthBlockDual<NSTOKES>& source) const {
        if constexpr (!std::is_same_v<S, SolarTransmissionExact>) {
            throw std::logic_error(
                "Solar transmission tables do not support wavelength "
                "batching");
        } else {
            if (!m_los_ground_is_hit.at(losidx)) {
                return;
            }

            const auto& first_layer = m_los_end_layers.at(losidx);
            const double mu_in = first_layer.exit.cos_zenith_angle(
                m_geometry.coordinates().sun_unit());
            if (mu_in <= 0.0) {
                return;
            }
            const double mu_out = -first_layer.exit.cos_zenith_angle(
                first_layer.average_look_away);
            const double phi_diff = first_layer.saz_exit;
            const int exit_index = m_index_map[losidx][0];
            const auto solar_trans = wavelength_head(
                m_solar_trans_batch[wavel_threadidx].row(exit_index), batch);

            for (int lane = 0; lane < batch.count; ++lane) {
                const int wavelength = batch.wavelength(lane);
                const auto brdf = m_atmosphere->surface().brdf(
                    wavelength, mu_in, mu_out, phi_diff);
                const Eigen::Vector<double, NSTOKES> source_value =
                    solar_trans(lane) * brdf(Eigen::placeholders::all, 0) *
                    mu_in;
                source.value.col(lane) += source_value;

                if (source.derivative_size() == 0) {
                    continue;
                }
                if (m_config->wf_precision() !=
                    sasktran2::Config::WeightingFunctionPrecision::limited) {
                    for (Eigen::SparseMatrix<double,
                                             Eigen::RowMajor>::InnerIterator
                             derivative(m_geometry_sparse, exit_index);
                         derivative; ++derivative) {
                        source.derivative(derivative.index(), batch)
                            .col(lane) -= derivative.value() * source_value;
                    }
                }
                for (int derivative = 0;
                     derivative < m_atmosphere->surface().num_deriv();
                     ++derivative) {
                    const auto brdf_derivative = m_atmosphere->surface().d_brdf(
                        wavelength, mu_in, mu_out, phi_diff, derivative);
                    source
                        .derivative(m_atmosphere->surface_deriv_start_index() +
                                        derivative,
                                    batch)
                        .col(lane) +=
                        solar_trans(lane) * mu_in *
                        brdf_derivative(Eigen::placeholders::all, 0);
                }
            }
        }
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::append_end_of_ray_active_derivatives(
        int losidx, std::vector<int>& derivative_indices) const {
        if constexpr (!std::is_same_v<S, SolarTransmissionExact>) {
            return;
        }

        if (!m_los_ground_is_hit.at(losidx)) {
            return;
        }

        if (m_config->wf_precision() !=
            sasktran2::Config::WeightingFunctionPrecision::limited) {
            const int exit_index = m_index_map[losidx][0];
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
                     derivative(m_geometry_sparse, exit_index);
                 derivative; ++derivative) {
                derivative_indices.push_back(derivative.index());
            }
        }

        for (int derivative = 0;
             derivative < m_atmosphere->surface().num_deriv(); ++derivative) {
            derivative_indices.push_back(
                m_atmosphere->surface_deriv_start_index() + derivative);
        }
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::append_interior_active_derivatives(
        int losidx, int layeridx, std::vector<int>& derivative_indices) const {
        if constexpr (!std::is_same_v<S, SolarTransmissionExact>) {
            return;
        }

        const auto& layer_active_derivatives =
            m_active_derivative_indices[losidx][layeridx];
        const bool use_lower_interpolation =
            m_geometry_1d != nullptr &&
            m_geometry_1d->altitude_grid().interpolation_method() ==
                grids::interpolation::lower;

        const std::vector<int>* start_active;
        const std::vector<int>* end_active;
        if (use_lower_interpolation) {
            const auto& layer = (*m_traced_rays)[losidx].layers[layeridx];
            const bool use_entrance_weights = layer.r_exit > layer.r_entrance;
            start_active = &layer_active_derivatives[1][use_entrance_weights];
            end_active = &layer_active_derivatives[0][use_entrance_weights];
        } else {
            start_active = &layer_active_derivatives[1][1];
            end_active = &layer_active_derivatives[0][0];
        }

        derivative_indices.insert(derivative_indices.end(),
                                  start_active->begin(), start_active->end());
        derivative_indices.insert(derivative_indices.end(), end_active->begin(),
                                  end_active->end());
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        ZoneScopedN("Initialize Single Scatter Source Geometry");
        m_traced_rays = &internal_viewing.traced_rays;
        this->m_solar_transmission.initialize_geometry(
            internal_viewing.traced_rays);

        if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
            {
                ZoneScopedN("Single Scatter Source Exact Geometry Matrix");
                if (m_geometry_2d == nullptr) {
                    // The 1D solar geometry is usually dense.
                    this->m_solar_transmission.generate_geometry_matrix(
                        internal_viewing.traced_rays, m_geometry_matrix,
                        m_ground_hit_flag);
                    m_geometry_sparse = m_geometry_matrix.sparseView();
                } else {
#ifdef SKTRAN_RUST_SUPPORT
                    // Higher-dimensional solar paths remain sparse.
                    m_geometry_matrix.resize(0, 0);
                    m_solar_transmission.generate_geometry_matrix(
                        internal_viewing.traced_rays, m_geometry_sparse,
                        m_ground_hit_flag);
#else
                    throw std::invalid_argument(
                        "Geometry2D exact solar transmission requires Rust "
                        "support");
#endif
                }
            }
        }
        if constexpr (std::is_same_v<S, SolarTransmissionTable>) {
            this->m_solar_transmission.generate_interpolation_matrix(
                internal_viewing.traced_rays, m_geometry_sparse,
                m_ground_hit_flag);
        }

        // We need some mapping between the layers inside each ray to our
        // calculated solar transmission
        m_index_map.resize(internal_viewing.traced_rays.size());
        int c = 0;
        for (int i = 0; i < internal_viewing.traced_rays.size(); ++i) {
            m_index_map[i].resize(
                internal_viewing.traced_rays[i].layers.size());

            for (int j = 0; j < m_index_map[i].size(); ++j) {
                m_index_map[i][j] = c;
                ++c;
            }
            // Final exit layer
            ++c;
        }
        {
            ZoneScopedN("Single Scatter Source Phase Geometry");
            this->m_phase_handler.initialize_geometry(
                internal_viewing.traced_rays, m_index_map);
        }

        m_los_ground_is_hit.resize(internal_viewing.traced_rays.size());
        m_los_end_layers.resize(internal_viewing.traced_rays.size());
        for (std::size_t ray_index = 0;
             ray_index < internal_viewing.traced_rays.size(); ++ray_index) {
            const auto& ray = internal_viewing.traced_rays[ray_index];
            m_los_ground_is_hit[ray_index] = ray.ground_is_hit;
            if (!ray.layers.empty()) {
                m_los_end_layers[ray_index] = ray.layers.front();
            }
        }
    }

    template <typename S, int NSTOKES>
    void
    SingleScatterSource<S, NSTOKES>::initialize_active_derivative_indices() {
        if (m_traced_rays == nullptr || m_atmosphere == nullptr) {
            throw std::runtime_error(
                "Single scatter geometry and atmosphere must be initialized "
                "before derivative sparsity");
        }

        m_active_derivative_indices.clear();
        m_active_derivative_indices.resize(m_traced_rays->size());

        const auto build_indices =
            [&](int solar_index,
                const sasktran2::raytracing::GridWeightStencilView&
                    geometry_weights,
                std::vector<int>& result) {
                result.clear();
                result.reserve(
                    m_geometry_sparse.innerVector(solar_index).nonZeros() +
                    geometry_weights.size() *
                        (2 + m_atmosphere->num_scattering_deriv_groups()));

                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
                         it(m_geometry_sparse, solar_index);
                     it; ++it) {
                    result.push_back(it.index());
                }

                for (std::size_t index = 0; index < geometry_weights.size();
                     ++index) {
                    const auto weight = geometry_weights[index];
                    if (weight.second == 0.0) {
                        continue;
                    }

                    result.push_back(weight.first);
                    result.push_back(m_atmosphere->ssa_deriv_start_index() +
                                     weight.first);
                    for (int derivative_group = 0;
                         derivative_group <
                         m_atmosphere->num_scattering_deriv_groups();
                         ++derivative_group) {
                        result.push_back(
                            m_atmosphere->scat_deriv_start_index() +
                            derivative_group * m_geometry.size() +
                            weight.first);
                    }
                }

                std::sort(result.begin(), result.end());
                result.erase(std::unique(result.begin(), result.end()),
                             result.end());
            };

        for (std::size_t los_index = 0; los_index < m_traced_rays->size();
             ++los_index) {
            const auto& ray = (*m_traced_rays)[los_index];
            auto& los_indices = m_active_derivative_indices[los_index];
            los_indices.resize(ray.layers.size());

            for (std::size_t layer_index = 0; layer_index < ray.layers.size();
                 ++layer_index) {
                const int exit_solar_index =
                    m_index_map[los_index][layer_index];
                const int entrance_solar_index = exit_solar_index + 1;
                const auto entrance_weights = ray.entrance_weights(layer_index);
                const auto exit_weights = ray.exit_weights(layer_index);
                auto& layer_indices = los_indices[layer_index];

                build_indices(exit_solar_index, exit_weights,
                              layer_indices[0][0]);
                build_indices(exit_solar_index, entrance_weights,
                              layer_indices[0][1]);
                build_indices(entrance_solar_index, exit_weights,
                              layer_indices[1][0]);
                build_indices(entrance_solar_index, entrance_weights,
                              layer_indices[1][1]);
            }
        }
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::integrated_source_constant(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::LayerGeometry& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        sasktran2::WavelengthBlockLaneDualView<NSTOKES, 1>& source,
        typename SourceTermInterface<NSTOKES>::IntegrationDirection direction)
        const {
        ZoneScopedN("Single Scatter Source Layer Calculation");

        bool calculate_derivatives = source.derivative_size() > 0;

        // Solar transmission is available at each layer boundary. Exact
        // transmission uses analytic log-linear weights below; table-based and
        // shadow-boundary cases retain the endpoint-average fallback.
        int exit_index = m_index_map[losidx][layeridx];
        int entrance_index = m_index_map[losidx][layeridx] + 1;

        double solar_trans_exit = m_solar_trans[wavel_threadidx](exit_index);
        double solar_trans_entrance =
            m_solar_trans[wavel_threadidx](entrance_index);

        auto& start_phase = m_start_source_cache[threadidx];
        auto& end_phase = m_end_source_cache[threadidx];

        const bool use_lower_interpolation =
            m_geometry_1d != nullptr &&
            m_geometry_1d->altitude_grid().interpolation_method() ==
                grids::interpolation::lower;
        const bool use_log_linear_solar_transmission =
            std::is_same_v<S, SolarTransmissionExact> &&
            solar_trans_entrance > 0.0 && solar_trans_exit > 0.0;
        const bool use_fused_exact_derivatives =
            calculate_derivatives &&
            std::is_same_v<S, SolarTransmissionExact> &&
            !use_log_linear_solar_transmission;

        if (use_log_linear_solar_transmission) {
            const double od = shell_od.od(0);
            const auto integration_weights =
                sasktran2::sourcealgo::single_scatter_source_weights(
                    od, shell_od.exp_minus_od(0), solar_trans_entrance,
                    solar_trans_exit, layer.od_quad_start_fraction,
                    layer.od_quad_end_fraction);

            const auto* start_weights = &entrance_weights;
            const auto* end_weights = &exit_weights;
            bool start_is_entrance = true;
            bool end_is_entrance = false;
            if (use_lower_interpolation) {
                if (layer.r_exit > layer.r_entrance) {
                    end_weights = &entrance_weights;
                    end_is_entrance = true;
                } else {
                    start_weights = &exit_weights;
                    start_is_entrance = false;
                }
            }

            const Eigen::Vector<double, NSTOKES> start_value =
                accumulate_exact_scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, *start_weights, start_is_entrance, 1.0,
                    *m_atmosphere,
                    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                        m_geometry_sparse, entrance_index),
                    layer.layer_distance * integration_weights.start, 0.0,
                    source);
            const Eigen::Vector<double, NSTOKES> end_value =
                accumulate_exact_scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, *end_weights, end_is_entrance, 1.0, *m_atmosphere,
                    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                        m_geometry_sparse, exit_index),
                    layer.layer_distance * integration_weights.end, 0.0,
                    source);

            const Eigen::Vector<double, NSTOKES> source_value =
                layer.layer_distance *
                (integration_weights.start * start_value +
                 integration_weights.end * end_value);
            source.value += source_value;

            if (calculate_derivatives) {
                const Eigen::Vector<double, NSTOKES> view_od_derivative =
                    layer.layer_distance *
                    (integration_weights.d_start_d_view_od * start_value +
                     integration_weights.d_end_d_view_od * end_value);
                for (auto derivative = shell_od.derivative_iterator();
                     derivative; ++derivative) {
                    source.deriv.col(derivative.index()) +=
                        derivative.value() * view_od_derivative;
                }

                const Eigen::Vector<double, NSTOKES> solar_start_od_derivative =
                    layer.layer_distance *
                    (integration_weights.d_start_d_solar_start_od *
                         start_value +
                     integration_weights.d_end_d_solar_start_od * end_value);
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
                         derivative(m_geometry_sparse, entrance_index);
                     derivative; ++derivative) {
                    source.deriv.col(derivative.index()) +=
                        derivative.value() * solar_start_od_derivative;
                }

                const Eigen::Vector<double, NSTOKES> solar_end_od_derivative =
                    layer.layer_distance *
                    (integration_weights.d_start_d_solar_end_od * start_value +
                     integration_weights.d_end_d_solar_end_od * end_value);
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
                         derivative(m_geometry_sparse, exit_index);
                     derivative; ++derivative) {
                    source.deriv.col(derivative.index()) +=
                        derivative.value() * solar_end_od_derivative;
                }
            }

#ifdef SASKTRAN_DEBUG_ASSERTS
            if (source_value.hasNaN()) {
                spdlog::error(
                    "NaN detected in log-linear single scatter source");
            }
#endif
            return;
        }

        if (use_fused_exact_derivatives) {
            const double od = shell_od.od(0);
            double source_factor;
            double source_factor_derivative;
            if (std::abs(od) < 1e-12) {
                source_factor = 1.0;
                source_factor_derivative = -0.5;
            } else {
                source_factor = -std::expm1(-od) / od;
                source_factor_derivative =
                    1 / od - source_factor * (1 + 1 / od);
            }

            const auto* start_weights = &entrance_weights;
            const auto* end_weights = &exit_weights;
            bool start_is_entrance = true;
            bool end_is_entrance = false;
            if (use_lower_interpolation) {
                if (layer.r_exit > layer.r_entrance) {
                    end_weights = &entrance_weights;
                    end_is_entrance = true;
                } else {
                    start_weights = &exit_weights;
                    start_is_entrance = false;
                }
            }

            const Eigen::Vector<double, NSTOKES> start_value =
                accumulate_exact_scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, *start_weights, start_is_entrance,
                    solar_trans_entrance, *m_atmosphere,
                    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                        m_geometry_sparse, entrance_index),
                    source_factor * layer.od_quad_start,
                    source_factor * layer.od_quad_start, source);
            const Eigen::Vector<double, NSTOKES> end_value =
                accumulate_exact_scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, *end_weights, end_is_entrance, solar_trans_exit,
                    *m_atmosphere,
                    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                        m_geometry_sparse, exit_index),
                    source_factor * layer.od_quad_end,
                    source_factor * layer.od_quad_end, source);

            const Eigen::Vector<double, NSTOKES> source_value =
                source_factor *
                (start_value.array() * layer.od_quad_start_fraction +
                 end_value.array() * layer.od_quad_end_fraction) *
                layer.layer_distance;
            source.value += source_value;

            for (auto derivative = shell_od.derivative_iterator(); derivative;
                 ++derivative) {
                source.deriv.col(derivative.index()) +=
                    derivative.value() * source_factor_derivative *
                    (start_value * layer.od_quad_start +
                     end_value * layer.od_quad_end);
            }

#ifdef SASKTRAN_DEBUG_ASSERTS
            if (source_value.hasNaN()) {
                spdlog::error(
                    "NaN detected in fused exact single scatter source");
            }
#endif
            return;
        }

        if (use_lower_interpolation) {
            if (layer.r_exit > layer.r_entrance) {
                scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, entrance_weights, true, solar_trans_entrance,
                    *m_atmosphere,
                    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                        m_geometry_sparse, entrance_index),
                    calculate_derivatives, start_phase);

                scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, entrance_weights, true, solar_trans_exit,
                    *m_atmosphere,
                    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                        m_geometry_sparse, exit_index),
                    calculate_derivatives, end_phase);
            } else {
                scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, exit_weights, false, solar_trans_entrance,
                    *m_atmosphere,
                    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                        m_geometry_sparse, entrance_index),
                    calculate_derivatives, start_phase);

                scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, exit_weights, false, solar_trans_exit,
                    *m_atmosphere,
                    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                        m_geometry_sparse, exit_index),
                    calculate_derivatives, end_phase);
            }
        } else {
            scattering_source(
                m_phase_handler, wavel_threadidx, losidx, layeridx, wavelidx,
                entrance_weights, true, solar_trans_entrance, *m_atmosphere,
                Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                    m_geometry_sparse, entrance_index),
                calculate_derivatives, start_phase);

            scattering_source(
                m_phase_handler, wavel_threadidx, losidx, layeridx, wavelidx,
                exit_weights, false, solar_trans_exit, *m_atmosphere,
                Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                    m_geometry_sparse, exit_index),
                calculate_derivatives, end_phase);
        }

        double source_factor1;
        double d_source_factor1;
        const double od = shell_od.od(0);
        if (std::abs(od) < 1e-12) {
            source_factor1 = 1.0;
            d_source_factor1 = -0.5;
        } else {
            source_factor1 = -std::expm1(-od) / od;
            d_source_factor1 = 1 / od - source_factor1 * (1 + 1 / od);
        }
        // Note dsource_factor = d_od * (1/od - source_factor * (1 + 1/od))

        // Get the phase matrix and add on the sources
        // The source factor term will only have extinction derivatives, the
        // phase term will have local SSA/scattering derivatives and is ~dense
        // in a 1D atmosphere

        // std::cout << solar_trans_entrance << " " << solar_trans_exit << "\n";

        Eigen::Vector<double, NSTOKES> source_val =
            source_factor1 *
            (start_phase.value.array() * layer.od_quad_start_fraction +
             end_phase.value.array() * layer.od_quad_end_fraction) *
            layer.layer_distance;

#ifdef SASKTRAN_DEBUG_ASSERTS
        if (source_val.hasNaN()) {
            static bool message = false;
            if (!message) {
                spdlog::error("SS Source NaN {} {} {} {} {} {}", source_factor1,
                              layer.od_quad_start_fraction,
                              layer.od_quad_end_fraction, layer.layer_distance,
                              start_phase.value(1), end_phase.value(1));
                message = true;
            }
            if constexpr (NSTOKES == 3) {
                source_val.setConstant(0.0);
            }
        }
#endif

        source.value.array() += source_val.array();

        if (calculate_derivatives) {
            // Now for the derivatives, start with dsource_factor which is
            // sparse
            for (auto it = shell_od.derivative_iterator(); it; ++it) {
                source.deriv(Eigen::placeholders::all, it.index()).array() +=
                    it.value() * d_source_factor1 *
                    (start_phase.value.array() * layer.od_quad_start +
                     end_phase.value.array() * layer.od_quad_end);
            }
            source.deriv.array() +=
                source_factor1 * start_phase.deriv.array() *
                    layer.od_quad_start +
                source_factor1 * end_phase.deriv.array() * layer.od_quad_end;
        }

#ifdef SASKTRAN_DEBUG_ASSERTS
        if (source.value.hasNaN()) {
            static bool message = false;
            if (!message) {
                spdlog::error("SS Source NaN {} {} {} {} {} {}", source_factor1,
                              layer.od_quad_start_fraction,
                              layer.od_quad_end_fraction, layer.layer_distance,
                              start_phase.value(1), end_phase.value(1));
                message = true;
            }
        }
#endif
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::integrated_source_single(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        sasktran2::WavelengthBlockLaneDualView<NSTOKES, 1>& source,
        typename SourceTermInterface<NSTOKES>::IntegrationDirection direction)
        const {
        if (layer.layer_distance < MINIMUM_SHELL_SIZE_M) {
            // Essentially an empty shell from rounding, don't have to do
            // anything
            return;
        }

        integrated_source_constant(wavelidx, losidx, layeridx, wavel_threadidx,
                                   threadidx, layer, entrance_weights,
                                   exit_weights, shell_od, source, direction);
    }

    template <typename S, int NSTOKES>
    template <int N>
    void SingleScatterSource<S, NSTOKES>::integrated_source_block(
        const sasktran2::WavelengthBlock<N>& batch, int losidx, int layeridx,
        int wavel_threadidx, int threadidx,
        const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        sasktran2::WavelengthBlockDual<NSTOKES>& source,
        typename SourceTermInterface<NSTOKES>::IntegrationDirection direction)
        const {
        if constexpr (!std::is_same_v<S, SolarTransmissionExact>) {
            throw std::logic_error(
                "Solar transmission tables do not support wavelength "
                "batching");
        } else {
            (void)direction;
            if (layer.layer_distance < MINIMUM_SHELL_SIZE_M) {
                return;
            }

            ZoneScopedN("Single Scatter Source Batch Layer Calculation");
            const int exit_index = m_index_map[losidx][layeridx];
            const int entrance_index = exit_index + 1;
            const auto solar_trans_exit = wavelength_head(
                m_solar_trans_batch[wavel_threadidx].row(exit_index), batch);
            const auto solar_trans_entrance = wavelength_head(
                m_solar_trans_batch[wavel_threadidx].row(entrance_index),
                batch);

            auto& integration_cache = m_batch_integration_cache[threadidx];
            auto source_factor =
                wavelength_head(integration_cache.source_factor, batch);
            auto source_factor_derivative = wavelength_head(
                integration_cache.source_factor_derivative, batch);
            auto endpoint_solar_transmission_start = wavelength_head(
                integration_cache.endpoint_solar_transmission_start, batch);
            auto endpoint_solar_transmission_end = wavelength_head(
                integration_cache.endpoint_solar_transmission_end, batch);
            auto start_weight =
                wavelength_head(integration_cache.start_weight, batch);
            auto end_weight =
                wavelength_head(integration_cache.end_weight, batch);
            auto d_start_d_view_od =
                wavelength_head(integration_cache.d_start_d_view_od, batch);
            auto d_end_d_view_od =
                wavelength_head(integration_cache.d_end_d_view_od, batch);
            auto d_start_d_solar_start_od = wavelength_head(
                integration_cache.d_start_d_solar_start_od, batch);
            auto d_end_d_solar_start_od = wavelength_head(
                integration_cache.d_end_d_solar_start_od, batch);
            auto d_start_d_solar_end_od = wavelength_head(
                integration_cache.d_start_d_solar_end_od, batch);
            auto d_end_d_solar_end_od =
                wavelength_head(integration_cache.d_end_d_solar_end_od, batch);
            auto start_derivative_scale = wavelength_head(
                integration_cache.start_derivative_scale, batch);
            auto end_derivative_scale =
                wavelength_head(integration_cache.end_derivative_scale, batch);
            auto start_solar_derivative_scale = wavelength_head(
                integration_cache.start_solar_derivative_scale, batch);
            auto end_solar_derivative_scale = wavelength_head(
                integration_cache.end_solar_derivative_scale, batch);
            for (int lane = 0; lane < batch.count; ++lane) {
                const double od = shell_od.od(lane);
                if (std::abs(od) < 1e-12) {
                    source_factor(lane) = 1.0;
                    source_factor_derivative(lane) = -0.5;
                } else {
                    source_factor(lane) = -std::expm1(-od) / od;
                    source_factor_derivative(lane) =
                        1 / od - source_factor(lane) * (1 + 1 / od);
                }

                const bool use_log_linear_solar_transmission =
                    solar_trans_entrance(lane) > 0.0 &&
                    solar_trans_exit(lane) > 0.0;
                if (use_log_linear_solar_transmission) {
                    const auto weights =
                        sasktran2::sourcealgo::single_scatter_source_weights(
                            od, shell_od.exp_minus_od(lane),
                            solar_trans_entrance(lane), solar_trans_exit(lane),
                            layer.od_quad_start_fraction,
                            layer.od_quad_end_fraction);
                    endpoint_solar_transmission_start(lane) = 1.0;
                    endpoint_solar_transmission_end(lane) = 1.0;
                    start_weight(lane) = weights.start;
                    end_weight(lane) = weights.end;
                    d_start_d_view_od(lane) = weights.d_start_d_view_od;
                    d_end_d_view_od(lane) = weights.d_end_d_view_od;
                    d_start_d_solar_start_od(lane) =
                        weights.d_start_d_solar_start_od;
                    d_end_d_solar_start_od(lane) =
                        weights.d_end_d_solar_start_od;
                    d_start_d_solar_end_od(lane) =
                        weights.d_start_d_solar_end_od;
                    d_end_d_solar_end_od(lane) = weights.d_end_d_solar_end_od;
                    start_derivative_scale(lane) =
                        layer.layer_distance * weights.start;
                    end_derivative_scale(lane) =
                        layer.layer_distance * weights.end;
                    start_solar_derivative_scale(lane) = 0.0;
                    end_solar_derivative_scale(lane) = 0.0;
                } else {
                    endpoint_solar_transmission_start(lane) =
                        solar_trans_entrance(lane);
                    endpoint_solar_transmission_end(lane) =
                        solar_trans_exit(lane);
                    start_weight(lane) = 0.0;
                    end_weight(lane) = 0.0;
                    d_start_d_view_od(lane) = 0.0;
                    d_end_d_view_od(lane) = 0.0;
                    d_start_d_solar_start_od(lane) = 0.0;
                    d_end_d_solar_start_od(lane) = 0.0;
                    d_start_d_solar_end_od(lane) = 0.0;
                    d_end_d_solar_end_od(lane) = 0.0;
                    start_derivative_scale(lane) =
                        source_factor(lane) * layer.od_quad_start;
                    end_derivative_scale(lane) =
                        source_factor(lane) * layer.od_quad_end;
                    start_solar_derivative_scale(lane) =
                        start_derivative_scale(lane);
                    end_solar_derivative_scale(lane) =
                        end_derivative_scale(lane);
                }
            }

            const auto* start_weights = &entrance_weights;
            const auto* end_weights = &exit_weights;
            bool start_is_entrance = true;
            bool end_is_entrance = false;
            const bool use_lower_interpolation =
                m_geometry_1d != nullptr &&
                m_geometry_1d->altitude_grid().interpolation_method() ==
                    grids::interpolation::lower;
            if (use_lower_interpolation) {
                if (layer.r_exit > layer.r_entrance) {
                    end_weights = &entrance_weights;
                    end_is_entrance = true;
                } else {
                    start_weights = &exit_weights;
                    start_is_entrance = false;
                }
            }

            auto& start_cache = m_batch_source_cache[threadidx][0];
            accumulate_exact_scattering_source_block<NSTOKES, N>(
                m_phase_handler, wavel_threadidx, losidx, layeridx, batch,
                *start_weights, start_is_entrance,
                endpoint_solar_transmission_start, *m_atmosphere,
                Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                    m_geometry_sparse, entrance_index),
                start_derivative_scale, start_solar_derivative_scale, source,
                start_cache);
            const auto start_value =
                wavelength_left_cols(start_cache.endpoint_source, batch);

            auto& end_cache = m_batch_source_cache[threadidx][1];
            accumulate_exact_scattering_source_block<NSTOKES, N>(
                m_phase_handler, wavel_threadidx, losidx, layeridx, batch,
                *end_weights, end_is_entrance, endpoint_solar_transmission_end,
                *m_atmosphere,
                Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                    m_geometry_sparse, exit_index),
                end_derivative_scale, end_solar_derivative_scale, source,
                end_cache);
            const auto end_value =
                wavelength_left_cols(end_cache.endpoint_source, batch);

            auto integrated_value =
                wavelength_left_cols(integration_cache.integrated_value, batch);
            for (int lane = 0; lane < batch.count; ++lane) {
                if (solar_trans_entrance(lane) > 0.0 &&
                    solar_trans_exit(lane) > 0.0) {
                    integrated_value.col(lane) =
                        layer.layer_distance *
                        (start_weight(lane) * start_value.col(lane) +
                         end_weight(lane) * end_value.col(lane));
                } else {
                    integrated_value.col(lane) =
                        source_factor(lane) * layer.layer_distance *
                        (layer.od_quad_start_fraction * start_value.col(lane) +
                         layer.od_quad_end_fraction * end_value.col(lane));
                }
            }
            wavelength_left_cols(source.value, batch) += integrated_value;

            if (source.derivative_size() > 0) {
                auto endpoint_quadrature = wavelength_left_cols(
                    integration_cache.endpoint_quadrature, batch);
                for (int lane = 0; lane < batch.count; ++lane) {
                    if (solar_trans_entrance(lane) > 0.0 &&
                        solar_trans_exit(lane) > 0.0) {
                        endpoint_quadrature.col(lane) =
                            layer.layer_distance *
                            (d_start_d_view_od(lane) * start_value.col(lane) +
                             d_end_d_view_od(lane) * end_value.col(lane));
                    } else {
                        endpoint_quadrature.col(lane) =
                            source_factor_derivative(lane) *
                            (layer.od_quad_start * start_value.col(lane) +
                             layer.od_quad_end * end_value.col(lane));
                    }
                }
                for (auto derivative = shell_od.derivative_iterator();
                     derivative; ++derivative) {
                    auto target_derivative =
                        source.derivative(derivative.index(), batch);
                    for (int lane = 0; lane < batch.count; ++lane) {
                        target_derivative.col(lane) +=
                            derivative.value() * endpoint_quadrature.col(lane);
                    }
                }

                for (int lane = 0; lane < batch.count; ++lane) {
                    endpoint_quadrature.col(lane) =
                        layer.layer_distance *
                        (d_start_d_solar_start_od(lane) *
                             start_value.col(lane) +
                         d_end_d_solar_start_od(lane) * end_value.col(lane));
                }
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
                         derivative(m_geometry_sparse, entrance_index);
                     derivative; ++derivative) {
                    auto target_derivative =
                        source.derivative(derivative.index(), batch);
                    for (int lane = 0; lane < batch.count; ++lane) {
                        target_derivative.col(lane) +=
                            derivative.value() * endpoint_quadrature.col(lane);
                    }
                }

                for (int lane = 0; lane < batch.count; ++lane) {
                    endpoint_quadrature.col(lane) =
                        layer.layer_distance *
                        (d_start_d_solar_end_od(lane) * start_value.col(lane) +
                         d_end_d_solar_end_od(lane) * end_value.col(lane));
                }
                for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
                         derivative(m_geometry_sparse, exit_index);
                     derivative; ++derivative) {
                    auto target_derivative =
                        source.derivative(derivative.index(), batch);
                    for (int lane = 0; lane < batch.count; ++lane) {
                        target_derivative.col(lane) +=
                            derivative.value() * endpoint_quadrature.col(lane);
                    }
                }
            }
        }
    }

#define SASKTRAN2_INSTANTIATE_SINGLE_SCATTER_BLOCK(NSTOKES, BLOCK_SIZE)        \
    template void SingleScatterSource<SolarTransmissionExact, NSTOKES>::       \
        calculate_block<BLOCK_SIZE>(                                           \
            const sasktran2::WavelengthBlock<BLOCK_SIZE>&, int);               \
    template void SingleScatterSource<SolarTransmissionExact, NSTOKES>::       \
        end_of_ray_source_block<BLOCK_SIZE>(                                   \
            const sasktran2::WavelengthBlock<BLOCK_SIZE>&, int, int, int,      \
            sasktran2::WavelengthBlockDual<NSTOKES>&) const;                   \
    template void SingleScatterSource<SolarTransmissionExact, NSTOKES>::       \
        integrated_source_block<BLOCK_SIZE>(                                   \
            const sasktran2::WavelengthBlock<BLOCK_SIZE>&, int, int, int, int, \
            const sasktran2::raytracing::TracedLayer&,                         \
            const sasktran2::raytracing::GridWeightStencilView&,               \
            const sasktran2::raytracing::GridWeightStencilView&,               \
            const sasktran2::WavelengthBlockODView&,                           \
            sasktran2::WavelengthBlockDual<NSTOKES>&,                          \
            SourceTermInterface<NSTOKES>::IntegrationDirection) const

    SASKTRAN2_INSTANTIATE_SINGLE_SCATTER_BLOCK(1, Eigen::Dynamic);
    SASKTRAN2_INSTANTIATE_SINGLE_SCATTER_BLOCK(1, 1);
    SASKTRAN2_INSTANTIATE_SINGLE_SCATTER_BLOCK(1, 4);
    SASKTRAN2_INSTANTIATE_SINGLE_SCATTER_BLOCK(3, Eigen::Dynamic);
    SASKTRAN2_INSTANTIATE_SINGLE_SCATTER_BLOCK(3, 1);
    SASKTRAN2_INSTANTIATE_SINGLE_SCATTER_BLOCK(3, 4);

#undef SASKTRAN2_INSTANTIATE_SINGLE_SCATTER_BLOCK

    template class SingleScatterSource<SolarTransmissionExact, 1>;
    template class SingleScatterSource<SolarTransmissionExact, 3>;

    template class SingleScatterSource<SolarTransmissionTable, 1>;
    template class SingleScatterSource<SolarTransmissionTable, 3>;
} // namespace sasktran2::solartransmission
