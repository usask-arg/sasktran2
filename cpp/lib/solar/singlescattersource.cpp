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
        initialize_atmosphere_impl(atmosphere, true);
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::initialize_atmosphere_native(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        initialize_atmosphere_impl(atmosphere, false);
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::initialize_atmosphere_impl(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
        bool materialized_derivative_storage) {
        const bool reusable_volume_state =
            atmosphere.revision() != 0 && m_has_atmosphere_volume_revision &&
            m_atmosphere_instance_id == atmosphere.instance_id() &&
            m_atmosphere_volume_revision == atmosphere.volume_revision();
        if (!reusable_volume_state) {
            std::fill(m_active_wavelength_block_count.begin(),
                      m_active_wavelength_block_count.end(), 0);
        }
        // Store the atmosphere for later
        m_atmosphere = &atmosphere;
        m_atmosphere_instance_id = atmosphere.instance_id();
        m_atmosphere_volume_revision = atmosphere.volume_revision();
        m_has_atmosphere_volume_revision = atmosphere.revision() != 0;
        // Composite engines materialize only the derivative mappings selected
        // for a native product. If no volume mapping is resident, atmospheric
        // single-scatter derivatives cannot contribute to the requested JVP
        // or VJP; ground BRDF derivatives remain active independently.
        m_native_volume_linearization_active =
            materialized_derivative_storage ||
            !atmosphere.storage().derivative_mappings_const().empty();
        if (!reusable_volume_state) {
            this->m_phase_handler.initialize_atmosphere(atmosphere);
        }

        if constexpr (compact_2d_table) {
            if (materialized_derivative_storage && atmosphere.num_deriv() > 0) {
                throw std::invalid_argument(
                    "Geometry2D table single scattering supports native JVP "
                    "and VJP linearization, not a materialized Jacobian");
            }
        }

        if constexpr (exact_transmission) {
            if (materialized_derivative_storage && atmosphere.num_deriv() > 0) {
                initialize_active_derivative_indices();
            } else {
                m_active_derivative_indices.clear();
            }

            if (!reusable_volume_state) {
                initialize_wavelength_blocks(m_wavelength_batch_capacity);
            }
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

        this->m_solar_transmission->initialize_config(config);
#ifdef SKTRAN_RUST_SUPPORT
        if constexpr (exact_transmission) {
            if (m_geometry_2d != nullptr && config.solar_refraction() &&
                m_shared_solar_table_2d == nullptr) {
                m_shared_solar_table_2d =
                    std::make_shared<SolarTransmissionTable2D>(*m_geometry_2d,
                                                               *m_raytracer_2d);
            }
            if (m_shared_solar_table_2d != nullptr) {
                m_shared_solar_table_2d->initialize_config(config);
            }
        }
#endif
        this->m_phase_handler.initialize_config(config);

        // Set up storage for each thread
        // m_solar_trans.resize(config.num_threads());
        m_solar_trans.resize(config.num_wavelength_threads());
        if constexpr (compact_2d_table) {
            m_solar_table_product.resize(config.num_wavelength_threads());
            m_solar_trans_jvp.resize(config.num_wavelength_threads());
            m_solar_endpoint_cotangent.resize(config.num_threads());
        }
        m_active_wavelength_block_start.assign(config.num_wavelength_threads(),
                                               0);
        m_active_wavelength_block_count.assign(config.num_wavelength_threads(),
                                               0);

        m_start_source_cache.resize(config.num_threads());
        m_end_source_cache.resize(config.num_threads());
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::calculate_single(int wavelidx,
                                                           int threadidx) {
        ZoneScopedN("Single Scatter Source Calculation");
        if (m_active_wavelength_block_start[threadidx] == wavelidx &&
            m_active_wavelength_block_count[threadidx] == 1) {
            return;
        }
        m_active_wavelength_block_start[threadidx] = wavelidx;
        m_active_wavelength_block_count[threadidx] = 1;
        // Don't have to do anything here
        m_phase_handler.calculate(wavelidx, threadidx);

        // Calculate the solar transmission at each cell
        if constexpr (exact_transmission) {
            // Faster to use the dense matrix if most of the elements are
            // nonzero
            if (m_geometry_matrix.size() > 0 &&
                double(m_geometry_sparse.non_zeros()) /
                        double(m_geometry_matrix.size()) >
                    DENSE_GEOMETRY_THRESHOLD) {
                m_solar_trans[threadidx].noalias() =
                    m_geometry_matrix *
                    m_atmosphere->storage().total_extinction(
                        Eigen::placeholders::all, wavelidx);
            } else {
                m_solar_trans[threadidx].resize(m_geometry_sparse.rows());
                const auto extinction =
                    m_atmosphere->storage().total_extinction(
                        Eigen::placeholders::all, wavelidx);
                m_geometry_sparse.multiply(extinction,
                                           m_solar_trans[threadidx]);
            }
        }

        if constexpr (std::is_same_v<S, SolarTransmissionTable>) {
            m_solar_trans[threadidx].noalias() =
                m_geometry_sparse.standard() *
                (m_solar_transmission->geometry_matrix() *
                 m_atmosphere->storage().total_extinction(
                     Eigen::placeholders::all, wavelidx));
        }

        if constexpr (compact_2d_table) {
            auto& table_od = m_solar_table_product[threadidx];
            table_od.resize(m_solar_transmission->table_size());
            m_solar_transmission->apply(
                m_atmosphere->storage().total_extinction(
                    Eigen::placeholders::all, wavelidx),
                table_od);
            m_solar_trans[threadidx].resize(m_solar_interpolation.rows());
            m_solar_interpolation.apply(table_od, m_solar_trans[threadidx]);
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
        if constexpr (!exact_transmission) {
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
        if constexpr (!exact_transmission) {
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
            if (m_active_wavelength_block_start[threadidx] == batch.start &&
                m_active_wavelength_block_count[threadidx] == batch.count) {
                return;
            }
            m_active_wavelength_block_start[threadidx] = batch.start;
            m_active_wavelength_block_count[threadidx] = batch.count;

            m_phase_handler.calculate_block(batch, threadidx);
            auto solar_trans =
                wavelength_left_cols(m_solar_trans_batch[threadidx], batch);
            const auto extinction = wavelength_middle_cols(
                m_atmosphere->storage().total_extinction, batch);
            if (m_geometry_matrix.size() > 0 &&
                double(m_geometry_sparse.non_zeros()) /
                        double(m_geometry_matrix.size()) >
                    DENSE_GEOMETRY_THRESHOLD) {
                solar_trans.noalias() = m_geometry_matrix * extinction;
            } else {
                m_geometry_sparse.multiply(extinction, solar_trans);
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
    double SingleScatterSource<S, NSTOKES>::solar_transmission_value(
        int wavelidx, int threadidx, int solar_index) const {
        if (m_wavelength_batch_capacity == 1) {
            return m_solar_trans.at(threadidx)(solar_index);
        }
        const int lane =
            wavelidx - m_active_wavelength_block_start.at(threadidx);
        if (lane < 0 || lane >= m_active_wavelength_block_count.at(threadidx)) {
            throw std::out_of_range(
                "Single-scatter wavelength is outside the active block");
        }
        return m_solar_trans_batch.at(threadidx)(solar_index, lane);
    }

    template <typename S, int NSTOKES>
    double SingleScatterSource<S, NSTOKES>::solar_transmission_tangent(
        int wavelidx, int threadidx, int solar_index,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) const {
        if constexpr (compact_2d_table) {
            return m_solar_trans_jvp.at(threadidx)(solar_index);
        } else {
            const double solar_trans =
                solar_transmission_value(wavelidx, threadidx, solar_index);
            return -solar_trans *
                   m_geometry_sparse.row_dot(solar_index, native_tangent);
        }
    }

    template <typename S, int NSTOKES>
    bool SingleScatterSource<S, NSTOKES>::ground_scattering_geometry(
        int losidx, double& mu_in, double& mu_out, double& phi_diff) const {
        const auto& first_layer = m_los_end_layers.at(losidx);
        Eigen::Vector3d direction_to_sun = m_geometry.coordinates().sun_unit();
        if (!m_solar_propagation_directions.empty()) {
            direction_to_sun =
                -m_solar_propagation_directions[m_index_map[losidx][0]];
        }
        sasktran2::raytracing::calculate_csz_saz(
            direction_to_sun.normalized(), first_layer.exit,
            first_layer.average_look_away, mu_in, phi_diff,
            m_geometry.coordinates().geometry_type());
        mu_out =
            -first_layer.exit.cos_zenith_angle(first_layer.average_look_away);
        return mu_in > 0.0;
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::calculate_jvp(
        const sasktran2::WavelengthBlock<>& block, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        calculate(block, threadidx);
        if constexpr (compact_2d_table) {
            if (block.count != 1 ||
                native_tangent.size() != m_atmosphere->num_deriv()) {
                throw std::invalid_argument(
                    "Invalid Geometry2D table single-scatter JVP block");
            }
            if (!m_native_volume_linearization_active) {
                return;
            }
            auto& table_tangent = m_solar_table_product[threadidx];
            table_tangent.resize(m_solar_transmission->table_size());
            m_solar_transmission->apply(
                native_tangent.head(m_solar_transmission->atmosphere_size()),
                table_tangent);
            auto& endpoint_tangent = m_solar_trans_jvp[threadidx];
            endpoint_tangent.resize(m_solar_interpolation.rows());
            m_solar_interpolation.apply(table_tangent, endpoint_tangent);
            endpoint_tangent.array() *= -m_solar_trans[threadidx].array();
        }
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::calculate_vjp(
        const sasktran2::WavelengthBlock<>& block, int threadidx) {
        calculate(block, threadidx);
        if constexpr (compact_2d_table) {
            if (block.count != 1) {
                throw std::invalid_argument(
                    "Geometry2D table single-scatter VJP requires scalar "
                    "wavelength blocks");
            }
            const int first_thread = threadidx;
            const int last_thread =
                first_thread + m_config->num_source_threads();
            if (last_thread >
                static_cast<int>(m_solar_endpoint_cotangent.size())) {
                throw std::logic_error(
                    "Single-scatter VJP thread layout is invalid");
            }
            for (int thread = first_thread; thread < last_thread; ++thread) {
                m_solar_endpoint_cotangent[thread].setZero();
            }
        }
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::finalize_vjp(
        const sasktran2::WavelengthBlock<>& block, int threadidx,
        Eigen::Ref<Eigen::MatrixXd> native_gradient) const {
        if constexpr (compact_2d_table) {
            if (block.count != 1 || native_gradient.cols() != 1) {
                throw std::invalid_argument(
                    "Invalid Geometry2D table single-scatter VJP result");
            }
            if (!m_native_volume_linearization_active) {
                return;
            }
            m_solar_endpoint_cotangent_sum.setZero(
                m_solar_interpolation.rows());
            const int last_thread = threadidx + m_config->num_source_threads();
            for (int thread = threadidx; thread < last_thread; ++thread) {
                m_solar_endpoint_cotangent_sum +=
                    m_solar_endpoint_cotangent[thread];
            }
            m_solar_table_cotangent.resize(m_solar_transmission->table_size());
            m_solar_interpolation.apply_transpose(
                m_solar_endpoint_cotangent_sum, m_solar_table_cotangent);
            auto extinction_gradient = native_gradient.col(0).head(
                m_solar_transmission->atmosphere_size());
            m_solar_transmission->accumulate_transpose(
                m_solar_table_cotangent, extinction_gradient, 1.0);
        }
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::endpoint_source_jvp(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int solar_index,
        const sasktran2::raytracing::GridWeightStencilView& weights,
        bool is_entrance, Eigen::Ref<const Eigen::VectorXd> native_tangent,
        sasktran2::RadianceJVP<NSTOKES>& result) const {
        if constexpr (!native_transmission_linearization) {
            throw std::logic_error(
                "Native JVP requires exact solar transmission");
        } else {
            const auto& storage = m_atmosphere->storage();
            double extinction = 0.0;
            double ssa = 0.0;
            double extinction_jvp = 0.0;
            double ssa_jvp = 0.0;
            for (std::size_t index = 0; index < weights.size(); ++index) {
                const auto weight = weights[index];
                if (weight.second == 0.0) {
                    continue;
                }
                extinction += storage.total_extinction(weight.first, wavelidx) *
                              weight.second;
                ssa += storage.ssa(weight.first, wavelidx) * weight.second;
                if (m_native_volume_linearization_active) {
                    extinction_jvp +=
                        native_tangent(weight.first) * weight.second;
                    ssa_jvp +=
                        native_tangent(m_atmosphere->ssa_deriv_start_index() +
                                       weight.first) *
                        weight.second;
                }
            }

            const double solar_trans = solar_transmission_value(
                wavelidx, wavel_threadidx, solar_index);
            if (!m_native_volume_linearization_active) {
                const auto phase = m_phase_handler.scatter_value(
                    wavel_threadidx, losidx, layeridx, wavelidx, weights,
                    is_entrance);
                result.value =
                    extinction * ssa * solar_trans / (EIGEN_PI * 4) * phase;
                result.jvp.setZero();
                return;
            }

            const double solar_trans_jvp = solar_transmission_tangent(
                wavelidx, wavel_threadidx, solar_index, native_tangent);
            Eigen::Vector<double, NSTOKES> phase;
            Eigen::Vector<double, NSTOKES> phase_jvp;
            m_phase_handler.scatter_jvp(wavel_threadidx, losidx, layeridx,
                                        wavelidx, weights, is_entrance,
                                        native_tangent, phase, phase_jvp);

            const double scale = 1.0 / (EIGEN_PI * 4);
            const double amplitude = extinction * ssa * solar_trans * scale;
            const double amplitude_jvp =
                scale * (extinction_jvp * ssa * solar_trans +
                         extinction * ssa_jvp * solar_trans +
                         extinction * ssa * solar_trans_jvp);
            result.value = amplitude * phase;
            result.jvp = amplitude_jvp * phase + amplitude * phase_jvp;
        }
    }

    template <typename S, int NSTOKES>
    Eigen::Vector<double, NSTOKES>
    SingleScatterSource<S, NSTOKES>::endpoint_source_vjp(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, int solar_index,
        const sasktran2::raytracing::GridWeightStencilView& weights,
        bool is_entrance, const Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        if constexpr (!native_transmission_linearization) {
            throw std::logic_error(
                "Native VJP requires exact solar transmission");
        } else {
            const auto& storage = m_atmosphere->storage();
            double extinction = 0.0;
            double ssa = 0.0;
            for (std::size_t index = 0; index < weights.size(); ++index) {
                const auto weight = weights[index];
                if (weight.second == 0.0) {
                    continue;
                }
                extinction += storage.total_extinction(weight.first, wavelidx) *
                              weight.second;
                ssa += storage.ssa(weight.first, wavelidx) * weight.second;
            }

            const double solar_trans = solar_transmission_value(
                wavelidx, wavel_threadidx, solar_index);
            const Eigen::Vector<double, NSTOKES> phase =
                m_phase_handler.scatter_value(wavel_threadidx, losidx, layeridx,
                                              wavelidx, weights, is_entrance);

            const double scale = 1.0 / (EIGEN_PI * 4);
            const double amplitude = extinction * ssa * solar_trans * scale;
            const double amplitude_cotangent = cotangent.dot(phase);
            const double extinction_cotangent =
                amplitude_cotangent * ssa * solar_trans * scale;
            const double ssa_cotangent =
                amplitude_cotangent * extinction * solar_trans * scale;
            const double solar_trans_cotangent =
                amplitude_cotangent * extinction * ssa * scale;

            for (std::size_t index = 0; index < weights.size(); ++index) {
                const auto weight = weights[index];
                if (weight.second == 0.0) {
                    continue;
                }
                native_gradient(weight.first) +=
                    weight.second * extinction_cotangent;
                native_gradient(m_atmosphere->ssa_deriv_start_index() +
                                weight.first) += weight.second * ssa_cotangent;
            }
            const double solar_od_cotangent =
                -solar_trans * solar_trans_cotangent;
            if constexpr (compact_2d_table) {
                m_solar_endpoint_cotangent[threadidx](solar_index) +=
                    solar_od_cotangent;
            } else {
                m_geometry_sparse.accumulate_row(
                    solar_index, solar_od_cotangent, native_gradient);
            }
            m_phase_handler.scatter_vjp(wavel_threadidx, losidx, layeridx,
                                        wavelidx, weights, is_entrance,
                                        amplitude * cotangent, native_gradient);
            return amplitude * phase;
        }
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::end_of_ray_source_jvp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        sasktran2::RadianceJVP<NSTOKES>& source) const {
        if constexpr (!native_transmission_linearization) {
            throw std::logic_error(
                "Native JVP requires exact solar transmission");
        } else {
            if (!m_los_ground_is_hit.at(losidx)) {
                return;
            }
            double mu_in;
            double mu_out;
            double phi_diff;
            if (!ground_scattering_geometry(losidx, mu_in, mu_out, phi_diff)) {
                return;
            }
            const int solar_index = m_index_map[losidx][0];
            const double solar_trans = solar_transmission_value(
                wavelidx, wavel_threadidx, solar_index);
            double solar_trans_jvp = 0.0;
            if (m_native_volume_linearization_active &&
                m_config->wf_precision() !=
                    sasktran2::Config::WeightingFunctionPrecision::limited) {
                solar_trans_jvp = solar_transmission_tangent(
                    wavelidx, wavel_threadidx, solar_index, native_tangent);
            }
            const auto brdf = m_atmosphere->surface().brdf(
                wavelidx, mu_in, mu_out, phi_diff,
                m_los_surface_interpolation_weights.at(losidx));
            Eigen::Matrix<double, NSTOKES, NSTOKES> brdf_jvp =
                Eigen::Matrix<double, NSTOKES, NSTOKES>::Zero();
            for (int derivative = 0;
                 derivative < m_atmosphere->surface().num_deriv();
                 ++derivative) {
                brdf_jvp +=
                    native_tangent(m_atmosphere->surface_deriv_start_index() +
                                   derivative) *
                    m_atmosphere->surface().d_brdf(
                        wavelidx, mu_in, mu_out, phi_diff, derivative,
                        m_los_surface_interpolation_weights.at(losidx));
            }
            source.value +=
                solar_trans * mu_in * brdf(Eigen::placeholders::all, 0);
            source.jvp +=
                mu_in * (solar_trans_jvp * brdf(Eigen::placeholders::all, 0) +
                         solar_trans * brdf_jvp(Eigen::placeholders::all, 0));
        }
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::end_of_ray_source_vjp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        const Eigen::Vector<double, NSTOKES>&,
        Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        if constexpr (!native_transmission_linearization) {
            throw std::logic_error(
                "Native VJP requires exact solar transmission");
        } else {
            if (!m_los_ground_is_hit.at(losidx)) {
                return;
            }
            double mu_in;
            double mu_out;
            double phi_diff;
            if (!ground_scattering_geometry(losidx, mu_in, mu_out, phi_diff)) {
                return;
            }
            const int solar_index = m_index_map[losidx][0];
            const double solar_trans = solar_transmission_value(
                wavelidx, wavel_threadidx, solar_index);
            const auto brdf = m_atmosphere->surface().brdf(
                wavelidx, mu_in, mu_out, phi_diff,
                m_los_surface_interpolation_weights.at(losidx));
            const double solar_trans_cotangent =
                mu_in * cotangent.dot(brdf(Eigen::placeholders::all, 0));
            if (m_native_volume_linearization_active &&
                m_config->wf_precision() !=
                    sasktran2::Config::WeightingFunctionPrecision::limited) {
                const double solar_od_cotangent =
                    -solar_trans * solar_trans_cotangent;
                if constexpr (compact_2d_table) {
                    m_solar_endpoint_cotangent[threadidx](solar_index) +=
                        solar_od_cotangent;
                } else {
                    m_geometry_sparse.accumulate_row(
                        solar_index, solar_od_cotangent, native_gradient);
                }
            }
            for (int derivative = 0;
                 derivative < m_atmosphere->surface().num_deriv();
                 ++derivative) {
                const auto brdf_derivative = m_atmosphere->surface().d_brdf(
                    wavelidx, mu_in, mu_out, phi_diff, derivative,
                    m_los_surface_interpolation_weights.at(losidx));
                native_gradient(m_atmosphere->surface_deriv_start_index() +
                                derivative) +=
                    solar_trans * mu_in *
                    cotangent.dot(brdf_derivative(Eigen::placeholders::all, 0));
            }
        }
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::end_of_ray_source_single(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        sasktran2::WavelengthBlockLaneDualView<NSTOKES, 1>& source) const {
        if (m_los_ground_is_hit.at(losidx)) {
            // Single scatter ground source is solar_trans * cos(th) * brdf
            double mu_in;
            double mu_out;
            double phi_diff;
            if (!ground_scattering_geometry(losidx, mu_in, mu_out, phi_diff)) {
                return;
            }

            Eigen::Matrix<double, NSTOKES, NSTOKES> brdf =
                m_atmosphere->surface().brdf(
                    wavelidx, mu_in, mu_out, phi_diff,
                    m_los_surface_interpolation_weights.at(losidx));

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
                if constexpr (exact_transmission) {
                    if (m_config->wf_precision() !=
                        sasktran2::Config::WeightingFunctionPrecision::
                            limited) {
                        // Have to apply the solar transmission derivative
                        // factors
                        for (SolarGeometryMatrix::InnerIterator it(
                                 m_geometry_sparse, exit_index);
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
                        m_atmosphere->surface().d_brdf(
                            wavelidx, mu_in, mu_out, phi_diff, k,
                            m_los_surface_interpolation_weights.at(losidx));

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
        if constexpr (!exact_transmission) {
            throw std::logic_error(
                "Solar transmission tables do not support wavelength "
                "batching");
        } else {
            if (!m_los_ground_is_hit.at(losidx)) {
                return;
            }

            double mu_in;
            double mu_out;
            double phi_diff;
            if (!ground_scattering_geometry(losidx, mu_in, mu_out, phi_diff)) {
                return;
            }
            const int exit_index = m_index_map[losidx][0];
            const auto solar_trans = wavelength_head(
                m_solar_trans_batch[wavel_threadidx].row(exit_index), batch);

            for (int lane = 0; lane < batch.count; ++lane) {
                const int wavelength = batch.wavelength(lane);
                const auto brdf = m_atmosphere->surface().brdf(
                    wavelength, mu_in, mu_out, phi_diff,
                    m_los_surface_interpolation_weights.at(losidx));
                const Eigen::Vector<double, NSTOKES> source_value =
                    solar_trans(lane) * brdf(Eigen::placeholders::all, 0) *
                    mu_in;
                source.value.col(lane) += source_value;

                if (source.derivative_size() == 0) {
                    continue;
                }
                if (m_config->wf_precision() !=
                    sasktran2::Config::WeightingFunctionPrecision::limited) {
                    for (SolarGeometryMatrix::InnerIterator derivative(
                             m_geometry_sparse, exit_index);
                         derivative; ++derivative) {
                        source.derivative(derivative.index(), batch)
                            .col(lane) -= derivative.value() * source_value;
                    }
                }
                for (int derivative = 0;
                     derivative < m_atmosphere->surface().num_deriv();
                     ++derivative) {
                    const auto brdf_derivative = m_atmosphere->surface().d_brdf(
                        wavelength, mu_in, mu_out, phi_diff, derivative,
                        m_los_surface_interpolation_weights.at(losidx));
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
        if constexpr (!exact_transmission) {
            return;
        }

        if (!m_los_ground_is_hit.at(losidx)) {
            return;
        }

        if (m_config->wf_precision() !=
            sasktran2::Config::WeightingFunctionPrecision::limited) {
            const int exit_index = m_index_map[losidx][0];
            for (SolarGeometryMatrix::InnerIterator derivative(
                     m_geometry_sparse, exit_index);
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
        if constexpr (!exact_transmission) {
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
        std::fill(m_active_wavelength_block_count.begin(),
                  m_active_wavelength_block_count.end(), 0);
        m_traced_rays = &internal_viewing.traced_rays;
        this->m_solar_transmission->initialize_geometry(
            internal_viewing.traced_rays);

        m_solar_propagation_directions.clear();
        if constexpr (exact_transmission) {
            {
                ZoneScopedN("Single Scatter Source Exact Geometry Matrix");
                if (m_geometry_2d == nullptr) {
                    // The 1D solar geometry is usually dense.
                    this->m_solar_transmission->generate_geometry_matrix(
                        internal_viewing.traced_rays, m_geometry_matrix,
                        m_ground_hit_flag);
                    m_geometry_sparse.use_standard() =
                        m_geometry_matrix.sparseView();
                } else {
#ifdef SKTRAN_RUST_SUPPORT
                    // Higher-dimensional solar paths remain sparse.
                    m_geometry_matrix.resize(0, 0);
                    if (m_config->solar_refraction()) {
                        m_shared_solar_table_2d->initialize_geometry(
                            internal_viewing.traced_rays);
                        m_shared_solar_table_2d->generate_solar_geometry(
                            internal_viewing.traced_rays, m_ground_hit_flag,
                            m_solar_propagation_directions);
                        m_solar_transmission
                            ->generate_refracted_geometry_matrix(
                                internal_viewing.traced_rays,
                                m_solar_propagation_directions,
                                m_geometry_sparse, m_ground_hit_flag);
                    } else {
                        m_solar_transmission->generate_geometry_matrix(
                            internal_viewing.traced_rays, m_geometry_sparse,
                            m_ground_hit_flag);
                    }
#else
                    throw std::invalid_argument(
                        "Geometry2D exact solar transmission requires Rust "
                        "support");
#endif
                }
            }
        }
        if constexpr (std::is_same_v<S, SolarTransmissionTable>) {
            this->m_solar_transmission->generate_interpolation_matrix(
                internal_viewing.traced_rays, m_geometry_sparse.use_standard(),
                m_ground_hit_flag);
        }
        if constexpr (compact_2d_table) {
            m_solar_transmission->generate_interpolation(
                internal_viewing.traced_rays, m_solar_interpolation,
                m_ground_hit_flag,
                m_config->solar_refraction() ? &m_solar_propagation_directions
                                             : nullptr);
            auto& empty_geometry = m_geometry_sparse.use_standard();
            empty_geometry.resize(m_solar_interpolation.rows(),
                                  m_geometry.size());
            empty_geometry.setZero();
            empty_geometry.makeCompressed();
            for (auto& cotangent : m_solar_endpoint_cotangent) {
                cotangent.setZero(m_solar_interpolation.rows());
            }
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
                internal_viewing.traced_rays, m_index_map,
                m_solar_propagation_directions.empty()
                    ? nullptr
                    : &m_solar_propagation_directions);
        }

        m_los_ground_is_hit.resize(internal_viewing.traced_rays.size());
        m_los_end_layers.resize(internal_viewing.traced_rays.size());
        m_los_surface_interpolation_weights.clear();
        m_los_surface_interpolation_weights.resize(
            internal_viewing.traced_rays.size());
        for (std::size_t ray_index = 0;
             ray_index < internal_viewing.traced_rays.size(); ++ray_index) {
            const auto& ray = internal_viewing.traced_rays[ray_index];
            m_los_ground_is_hit[ray_index] = ray.ground_is_hit;
            if (!ray.layers.empty()) {
                m_los_end_layers[ray_index] = ray.layers.front();
                if (ray.ground_is_hit && m_geometry_2d != nullptr) {
                    m_geometry_2d->assign_horizontal_interpolation_weights(
                        ray.layers.front().exit,
                        m_los_surface_interpolation_weights[ray_index]);
                }
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
                    m_geometry_sparse.row_nonzeros(solar_index) +
                    geometry_weights.size() *
                        (2 + m_atmosphere->num_scattering_deriv_groups()));

                for (SolarGeometryMatrix::InnerIterator it(m_geometry_sparse,
                                                           solar_index);
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
        ZoneScopedN("Single Scatter Source Constant Calculation");

        bool calculate_derivatives = source.derivative_size() > 0;

        // Integrates assuming the source is constant in the layer and
        // determined by the average of the layer boundaries
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
        const bool use_fused_exact_derivatives =
            calculate_derivatives && exact_transmission;

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
                    SolarGeometryMatrix::InnerIterator(m_geometry_sparse,
                                                       entrance_index),
                    source_factor * layer.od_quad_start, source);
            const Eigen::Vector<double, NSTOKES> end_value =
                accumulate_exact_scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, *end_weights, end_is_entrance, solar_trans_exit,
                    *m_atmosphere,
                    SolarGeometryMatrix::InnerIterator(m_geometry_sparse,
                                                       exit_index),
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
                scattering_source(m_phase_handler, wavel_threadidx, losidx,
                                  layeridx, wavelidx, entrance_weights, true,
                                  solar_trans_entrance, *m_atmosphere,
                                  SolarGeometryMatrix::InnerIterator(
                                      m_geometry_sparse, entrance_index),
                                  calculate_derivatives, start_phase);

                scattering_source(m_phase_handler, wavel_threadidx, losidx,
                                  layeridx, wavelidx, entrance_weights, true,
                                  solar_trans_exit, *m_atmosphere,
                                  SolarGeometryMatrix::InnerIterator(
                                      m_geometry_sparse, exit_index),
                                  calculate_derivatives, end_phase);
            } else {
                scattering_source(m_phase_handler, wavel_threadidx, losidx,
                                  layeridx, wavelidx, exit_weights, false,
                                  solar_trans_entrance, *m_atmosphere,
                                  SolarGeometryMatrix::InnerIterator(
                                      m_geometry_sparse, entrance_index),
                                  calculate_derivatives, start_phase);

                scattering_source(m_phase_handler, wavel_threadidx, losidx,
                                  layeridx, wavelidx, exit_weights, false,
                                  solar_trans_exit, *m_atmosphere,
                                  SolarGeometryMatrix::InnerIterator(
                                      m_geometry_sparse, exit_index),
                                  calculate_derivatives, end_phase);
            }
        } else {
            scattering_source(m_phase_handler, wavel_threadidx, losidx,
                              layeridx, wavelidx, entrance_weights, true,
                              solar_trans_entrance, *m_atmosphere,
                              SolarGeometryMatrix::InnerIterator(
                                  m_geometry_sparse, entrance_index),
                              calculate_derivatives, start_phase);

            scattering_source(m_phase_handler, wavel_threadidx, losidx,
                              layeridx, wavelidx, exit_weights, false,
                              solar_trans_exit, *m_atmosphere,
                              SolarGeometryMatrix::InnerIterator(
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
    void SingleScatterSource<S, NSTOKES>::integrated_source_jvp(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        sasktran2::RadianceJVP<NSTOKES>& source) const {
        (void)threadidx;
        if constexpr (!native_transmission_linearization) {
            throw std::logic_error(
                "Native JVP requires exact solar transmission");
        } else {
            if (layer.layer_distance < MINIMUM_SHELL_SIZE_M) {
                return;
            }
            const int exit_index = m_index_map[losidx][layeridx];
            const int entrance_index = exit_index + 1;
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

            sasktran2::RadianceJVP<NSTOKES> start;
            sasktran2::RadianceJVP<NSTOKES> end;
            endpoint_source_jvp(wavelidx, losidx, layeridx, wavel_threadidx,
                                entrance_index, *start_weights,
                                start_is_entrance, native_tangent, start);
            endpoint_source_jvp(wavelidx, losidx, layeridx, wavel_threadidx,
                                exit_index, *end_weights, end_is_entrance,
                                native_tangent, end);

            const double od = shell_od.od(0);
            double factor;
            double factor_derivative;
            if (std::abs(od) < 1e-12) {
                factor = 1.0;
                factor_derivative = -0.5;
            } else {
                factor = -std::expm1(-od) / od;
                factor_derivative = 1 / od - factor * (1 + 1 / od);
            }
            // Match the established source-value quadrature exactly.  The
            // optical-depth endpoint coefficients are not source blending
            // weights for lower interpolation: one is the full path length
            // and the other is zero, while the source remains the average of
            // its endpoint values.
            const auto endpoint_value =
                (start.value * layer.od_quad_start_fraction +
                 end.value * layer.od_quad_end_fraction) *
                layer.layer_distance;
            source.value += factor * endpoint_value;
            if (!m_native_volume_linearization_active) {
                return;
            }

            double od_jvp = 0.0;
            for (auto derivative = shell_od.derivative_iterator(); derivative;
                 ++derivative) {
                od_jvp +=
                    derivative.value() * native_tangent(derivative.index());
            }
            const auto optical_depth_endpoint_value =
                start.value * layer.od_quad_start +
                end.value * layer.od_quad_end;
            const auto endpoint_jvp =
                start.jvp * layer.od_quad_start + end.jvp * layer.od_quad_end;
            source.jvp +=
                factor * endpoint_jvp +
                factor_derivative * od_jvp * optical_depth_endpoint_value;
        }
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::integrated_source_vjp(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        const Eigen::Vector<double, NSTOKES>&,
        Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        if constexpr (!native_transmission_linearization) {
            throw std::logic_error(
                "Native VJP requires exact solar transmission");
        } else {
            if (layer.layer_distance < MINIMUM_SHELL_SIZE_M) {
                return;
            }
            if (!m_native_volume_linearization_active) {
                return;
            }
            const int exit_index = m_index_map[losidx][layeridx];
            const int entrance_index = exit_index + 1;
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

            const double od = shell_od.od(0);
            double factor;
            double factor_derivative;
            if (std::abs(od) < 1e-12) {
                factor = 1.0;
                factor_derivative = -0.5;
            } else {
                factor = -std::expm1(-od) / od;
                factor_derivative = 1 / od - factor * (1 + 1 / od);
            }
            const auto start_value = endpoint_source_vjp(
                wavelidx, losidx, layeridx, wavel_threadidx, threadidx,
                entrance_index, *start_weights, start_is_entrance,
                factor * layer.od_quad_start * cotangent, native_gradient);
            const auto end_value = endpoint_source_vjp(
                wavelidx, losidx, layeridx, wavel_threadidx, threadidx,
                exit_index, *end_weights, end_is_entrance,
                factor * layer.od_quad_end * cotangent, native_gradient);
            const double od_cotangent =
                factor_derivative *
                cotangent.dot(start_value * layer.od_quad_start +
                              end_value * layer.od_quad_end);
            for (auto derivative = shell_od.derivative_iterator(); derivative;
                 ++derivative) {
                native_gradient(derivative.index()) +=
                    derivative.value() * od_cotangent;
            }
        }
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
        if constexpr (!exact_transmission) {
            throw std::logic_error(
                "Solar transmission tables do not support wavelength "
                "batching");
        } else {
            (void)direction;
            if (layer.layer_distance < MINIMUM_SHELL_SIZE_M) {
                return;
            }

            ZoneScopedN("Single Scatter Source Batch Constant Calculation");
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

            auto start_derivative_scale = wavelength_head(
                integration_cache.start_derivative_scale, batch);
            start_derivative_scale = source_factor * layer.od_quad_start;
            auto& start_cache = m_batch_source_cache[threadidx][0];
            accumulate_exact_scattering_source_block<NSTOKES, N>(
                m_phase_handler, wavel_threadidx, losidx, layeridx, batch,
                *start_weights, start_is_entrance, solar_trans_entrance,
                *m_atmosphere,
                SolarGeometryMatrix::InnerIterator(m_geometry_sparse,
                                                   entrance_index),
                start_derivative_scale, source, start_cache);
            const auto start_value =
                wavelength_left_cols(start_cache.endpoint_source, batch);

            auto end_derivative_scale =
                wavelength_head(integration_cache.end_derivative_scale, batch);
            end_derivative_scale = source_factor * layer.od_quad_end;
            auto& end_cache = m_batch_source_cache[threadidx][1];
            accumulate_exact_scattering_source_block<NSTOKES, N>(
                m_phase_handler, wavel_threadidx, losidx, layeridx, batch,
                *end_weights, end_is_entrance, solar_trans_exit, *m_atmosphere,
                SolarGeometryMatrix::InnerIterator(m_geometry_sparse,
                                                   exit_index),
                end_derivative_scale, source, end_cache);
            const auto end_value =
                wavelength_left_cols(end_cache.endpoint_source, batch);

            auto integrated_value =
                wavelength_left_cols(integration_cache.integrated_value, batch);
            integrated_value = start_value * layer.od_quad_start_fraction +
                               end_value * layer.od_quad_end_fraction;
            integrated_value.array().rowwise() *= source_factor.array();
            integrated_value *= layer.layer_distance;
            wavelength_left_cols(source.value, batch) += integrated_value;

            if (source.derivative_size() > 0) {
                auto endpoint_quadrature = wavelength_left_cols(
                    integration_cache.endpoint_quadrature, batch);
                endpoint_quadrature = start_value * layer.od_quad_start +
                                      end_value * layer.od_quad_end;
                for (auto derivative = shell_od.derivative_iterator();
                     derivative; ++derivative) {
                    auto target_derivative =
                        source.derivative(derivative.index(), batch);
                    target_derivative.array() +=
                        derivative.value() *
                        (endpoint_quadrature.array().rowwise() *
                         source_factor_derivative.array());
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
#ifdef SKTRAN_RUST_SUPPORT
    template class SingleScatterSource<SolarTransmissionTable2D, 1>;
    template class SingleScatterSource<SolarTransmissionTable2D, 3>;
#endif
} // namespace sasktran2::solartransmission
