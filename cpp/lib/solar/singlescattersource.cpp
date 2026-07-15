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
        m_start_active_derivative_indices.resize(config.num_threads());
        m_end_active_derivative_indices.resize(config.num_threads());

        m_solar_trans.resize(config.num_wavelength_threads());

        m_start_source_cache.resize(config.num_threads());
        m_end_source_cache.resize(config.num_threads());
    }

    template <typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::calculate(int wavelidx,
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
    void SingleScatterSource<S, NSTOKES>::end_of_ray_source(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source)
        const {
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
    void SingleScatterSource<S, NSTOKES>::initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        ZoneScopedN("Initialize Single Scatter Source Geometry");
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
    void SingleScatterSource<S, NSTOKES>::integrated_source_constant(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::LayerGeometry& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::SparseODDualView& shell_od,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source,
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
        auto& start_active_derivative_indices =
            m_start_active_derivative_indices[threadidx];
        auto& end_active_derivative_indices =
            m_end_active_derivative_indices[threadidx];

        const bool use_active_derivatives = m_geometry_2d != nullptr;
        const bool use_lower_interpolation =
            m_geometry_1d != nullptr &&
            m_geometry_1d->altitude_grid().interpolation_method() ==
                grids::interpolation::lower;

        if (use_lower_interpolation) {
            if (layer.r_exit > layer.r_entrance) {
                scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, entrance_weights, true, solar_trans_entrance,
                    *m_atmosphere,
                    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                        m_geometry_sparse, entrance_index),
                    calculate_derivatives, use_active_derivatives,
                    start_active_derivative_indices, start_phase);

                scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, entrance_weights, true, solar_trans_exit,
                    *m_atmosphere,
                    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                        m_geometry_sparse, exit_index),
                    calculate_derivatives, use_active_derivatives,
                    end_active_derivative_indices, end_phase);
            } else {
                scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, exit_weights, false, solar_trans_entrance,
                    *m_atmosphere,
                    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                        m_geometry_sparse, entrance_index),
                    calculate_derivatives, use_active_derivatives,
                    start_active_derivative_indices, start_phase);

                scattering_source(
                    m_phase_handler, wavel_threadidx, losidx, layeridx,
                    wavelidx, exit_weights, false, solar_trans_exit,
                    *m_atmosphere,
                    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                        m_geometry_sparse, exit_index),
                    calculate_derivatives, use_active_derivatives,
                    end_active_derivative_indices, end_phase);
            }
        } else {
            scattering_source(
                m_phase_handler, wavel_threadidx, losidx, layeridx, wavelidx,
                entrance_weights, true, solar_trans_entrance, *m_atmosphere,
                Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                    m_geometry_sparse, entrance_index),
                calculate_derivatives, use_active_derivatives,
                start_active_derivative_indices, start_phase);

            scattering_source(
                m_phase_handler, wavel_threadidx, losidx, layeridx, wavelidx,
                exit_weights, false, solar_trans_exit, *m_atmosphere,
                Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(
                    m_geometry_sparse, exit_index),
                calculate_derivatives, use_active_derivatives,
                end_active_derivative_indices, end_phase);
        }

        double source_factor1;
        double d_source_factor1;
        if (std::abs(shell_od.od) < 1e-12) {
            source_factor1 = 1.0;
            d_source_factor1 = -0.5;
        } else {
            source_factor1 = -std::expm1(-shell_od.od) / shell_od.od;
            d_source_factor1 =
                1 / shell_od.od - source_factor1 * (1 + 1 / shell_od.od);
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
            for (auto it = shell_od.deriv_iter; it; ++it) {
                source.deriv(Eigen::placeholders::all, it.index()).array() +=
                    it.value() * d_source_factor1 *
                    (start_phase.value.array() * layer.od_quad_start +
                     end_phase.value.array() * layer.od_quad_end);
            }
            // And add on d_phase. Structured 2D endpoints touch only a small
            // subset of the full atmosphere derivative vector, so avoid a
            // dense pass over every horizontal grid location for each layer.
            if (use_active_derivatives) {
                for (const int derivative_index :
                     start_active_derivative_indices) {
                    source.deriv.col(derivative_index) +=
                        source_factor1 * layer.od_quad_start *
                        start_phase.deriv.col(derivative_index);
                }
                for (const int derivative_index :
                     end_active_derivative_indices) {
                    source.deriv.col(derivative_index) +=
                        source_factor1 * layer.od_quad_end *
                        end_phase.deriv.col(derivative_index);
                }
            } else {
                source.deriv.array() +=
                    source_factor1 * start_phase.deriv.array() *
                        layer.od_quad_start +
                    source_factor1 * end_phase.deriv.array() *
                        layer.od_quad_end;
            }
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
    void SingleScatterSource<S, NSTOKES>::integrated_source(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::SparseODDualView& shell_od,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source,
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

    template class SingleScatterSource<SolarTransmissionExact, 1>;
    template class SingleScatterSource<SolarTransmissionExact, 3>;

    template class SingleScatterSource<SolarTransmissionTable, 1>;
    template class SingleScatterSource<SolarTransmissionTable, 3>;
} // namespace sasktran2::solartransmission
