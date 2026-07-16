#pragma once

#include "sasktran2/atmosphere/grid_storage.h"
#include "sasktran2/geometry.h"
#include <sasktran2/source_interface.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/atmosphere/atmosphere.h>
#include <sasktran2/config.h>
#include <sasktran2/dual.h>
#include <array>

namespace sasktran2::solartransmission {
    /** Generates one row of what is known as the geometry matrix.  The optical
     * depth at every necessary integration point can be written as \f$M k\f$
     * where \f$k\f$ is the extinction coefficient vector.  In 1D geometry, the
     *  matrix \f$M\f$ is densish, but it may be sparse in higher dimensions.
     * This function constructs one row of the matrix \f$M\f$, which corresponds
     * to one traced ray.
     *
     * @param row row index
     * @param traced_ray ray traced to the sun
     * @param result matrix to store the result in
     */
    inline void assign_dense_matrix_column(
        int row, const sasktran2::raytracing::TracedRay& traced_ray,
        Eigen::MatrixXd& result) {
        for (int i = 0; i < traced_ray.layers.size(); ++i) {
            const auto weights = traced_ray.optical_depth_weights(i);
            for (std::size_t j = 0; j < weights.size(); ++j) {
                const auto weight = weights[j];
                result(row, weight.first) += weight.second;
            }
        }
    }

    class SolarTransmissionBase {
      protected:
        const Geometry& m_geometry;
        const Geometry1D* m_geometry_1d = nullptr;
        const sasktran2::raytracing::RayTracerBase* m_raytracer = nullptr;
#ifdef SKTRAN_RUST_SUPPORT
        const Geometry2D* m_geometry_2d = nullptr;
        const sasktran2::raytracing::RustRayTracer2D* m_raytracer_2d = nullptr;
#endif

      public:
        SolarTransmissionBase(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : m_geometry(geometry), m_geometry_1d(&geometry),
              m_raytracer(&raytracer) {}

#ifdef SKTRAN_RUST_SUPPORT
        SolarTransmissionBase(
            const Geometry2D& geometry,
            const sasktran2::raytracing::RustRayTracer2D& raytracer)
            : m_geometry(geometry), m_geometry_2d(&geometry),
              m_raytracer_2d(&raytracer) {}
#endif
    };

    class SolarTransmissionExact : public SolarTransmissionBase {
      private:
      public:
        SolarTransmissionExact(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : SolarTransmissionBase(geometry, raytracer) {}

        virtual void initialize_config(const sasktran2::Config& config){};

        virtual void
        initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>&
                                integration_rays){};

        void generate_geometry_matrix(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            Eigen::MatrixXd& od_matrix,
            std::vector<bool>& ground_hit_flag) const;

#ifdef SKTRAN_RUST_SUPPORT
        SolarTransmissionExact(
            const Geometry2D& geometry,
            const sasktran2::raytracing::RustRayTracer2D& raytracer)
            : SolarTransmissionBase(geometry, raytracer) {}

        void generate_geometry_matrix(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            Eigen::SparseMatrix<double, Eigen::RowMajor>& od_matrix,
            std::vector<bool>& ground_hit_flag) const;
#endif
    };

    class SolarTransmissionTable : public SolarTransmissionExact {
      private:
        std::unique_ptr<sasktran2::grids::SourceLocationInterpolator>
            m_location_interpolator;
        const sasktran2::Config* m_config;
        Eigen::MatrixXd m_geometry_matrix;

        std::vector<bool> m_ground_hit_flag;

      public:
        SolarTransmissionTable(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : SolarTransmissionExact(geometry, raytracer) {}

        void initialize_config(const sasktran2::Config& config) override {
            m_config = &config;
        };

        void
        initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>&
                                integration_rays) override;

        void generate_interpolation_matrix(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            Eigen::SparseMatrix<double, Eigen::RowMajor>& interpolator,
            std::vector<bool>& ground_hit_flag) const;

        const Eigen::MatrixXd& geometry_matrix() const {
            return m_geometry_matrix;
        }
    };

    /**
     * The PhaseHandler is responsible for constructing the phase function for
     * the single scatter source term. The single scatter source term is needed
     * on a set of (wavelength, cos_angle, geometry) points.
     *
     * Construction can be done in one of two ways, depending on the
     * configuration and user input. The default is to construct the phase
     * function through the user input Legendre coefficients.
     *
     * The other option is for the user to directly input the phase function.
     * This is useful for cases where the phase function requires many terms in
     * the Legendre series.
     *
     * Storage for the phase function is handled in a semi-complicated hard to
     * understand manner. For each thread, we store the phase function (stokes,
     * internal_index) where internal_index represents a single scattering angle
     * at a single geomtry level in the atmosphere.
     *
     * The internal_index is mapped back to the geometry index through the
     * m_internal_to_geometry, and the scattering angle is determined by the
     * m_internal_to_cos_scatter.  This combination lets us calculate the phase
     * function.  But to actually use it, we need to map the entrance and exit
     * points of the ray to the internal indices.  This is done through the
     * m_geometry_entrance_to_internal and m_geometry_exit_to_internal
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> class PhaseHandler {
      private:
        using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>;

        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere;
        const sasktran2::Config* m_config;
        const sasktran2::Geometry& m_geometry;

        std::vector<std::array<double, NSTOKES>>
            m_scatter_angles; /** Full list of  scattering angles that we need.
                                 for NSTOKES =3 this is (cos_scatter, C1, C2) */

        Eigen::MatrixXd m_wigner_d00; /** Wigner D matrix for the phase function
                                         (legendre_order, scatter_angle) */
        Eigen::MatrixXd m_wigner_d02; /** Wigner D matrix for the phase function
                                         (legendre_order, scatter_angle) */

        // Internal phase functions and derivatives on the actual grid
        Eigen::Tensor<double, 3>
            m_phase; /** (stokes eq, internal_index, thread) **/
        Eigen::Tensor<double, 4>
            m_d_phase; /** (stokes eq, internal_index, deriv, thread) **/

        // Batch storage is wavelength-contiguous. The row mappings are
        // phase_component * num_internal + internal_index and
        // (deriv * num_phase_components + phase_component) * num_internal +
        // internal_index, respectively.
        std::vector<RowMajorMatrix> m_phase_batch;
        std::vector<RowMajorMatrix> m_d_phase_batch;
        int m_wavelength_batch_capacity = 0;

        std::vector<std::vector<std::vector<int>>>
            m_geometry_entrance_to_internal; /** Mapping from layer entrances to
                                                internal,
                                                [los][layer][interp_index] */
        std::vector<std::vector<std::vector<int>>>
            m_geometry_exit_to_internal;         /** Mapping from layer exits to
                                                    internal, [los][layer][interp_index]
                                                  */
        std::vector<int> m_internal_to_geometry; /** Maps the internal index
                                                      to the geometry index */
        std::vector<int>
            m_internal_to_cos_scatter; /** Determines what scattering angle to
                                          use for each internal index */

      public:
        PhaseHandler(const Geometry& geometry) : m_geometry(geometry) {}

        /**
         * Initializes the phase handler with the configuration object
         */
        void initialize_config(const sasktran2::Config& config) {
            m_config = &config;
        }

        /**
         *  Initializes the phase handler with the atmosphere object
         *
         *  @param atmosphere The atmosphere object
         */
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere);

        /**
         * Initializes the phase handler with the geometry object
         *
         * @param los_rays The traced line of sight rays
         * @param index_map The index map
         */
        void initialize_geometry(
            const std::vector<sasktran2::raytracing::TracedRay>& los_rays,
            const std::vector<std::vector<int>>& index_map);

        /**
         *   Calculates the phase function from the legendre coefficients at the
         * necessary scatter angles
         *
         *   @param threadidx The thread index
         *   @param wavelidx The wavelength index
         */
        void calculate(int wavelidx, int threadidx);

        void initialize_wavelength_blocks(int batch_size);

        void calculate_block(const sasktran2::WavelengthBlock& batch,
                             int threadidx);

        /**
         * Calculates the phase function at a given point and puts it into
         * source
         *
         * @param wavelidx The wavelength index
         * @param losidx The line of sight index
         * @param layeridx The layer index
         * @param index_weights The interpolation weights
         * @param is_entrance  True if we are at the entrance to a layer, false
         * if we are at the exit to a layer
         * @param source The source term
         */
        void scatter(int wavelidx, int losidx, int layeridx,
                     const raytracing::GridWeightStencilView& index_weights,
                     bool is_entrance,
                     sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                     NSTOKES>& source) const;

        /** Returns the interpolated phase vector and accumulates its
         * scattering-property derivatives directly into a target source. */
        Eigen::Vector<double, NSTOKES> scatter_and_accumulate_derivative(
            int threadidx, int losidx, int layeridx,
            const raytracing::GridWeightStencilView& index_weights,
            bool is_entrance, double source_amplitude, double derivative_scale,
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                target) const;

        void scatter_and_accumulate_derivative_block(
            int threadidx, int losidx, int layeridx,
            const raytracing::GridWeightStencilView& index_weights,
            bool is_entrance,
            const Eigen::Ref<const Eigen::RowVectorXd>& source_amplitude,
            const Eigen::Ref<const Eigen::RowVectorXd>& derivative_scale,
            sasktran2::WavelengthBlockDual<NSTOKES>& target,
            Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>&
                phase_result) const;

      private:
        void initialize_geometry_impl(
            const std::vector<sasktran2::raytracing::TracedRay>& los_rays,
            const std::vector<std::vector<int>>& index_map);

        void
        scatter_impl(int wavelidx, int losidx, int layeridx,
                     const raytracing::GridWeightStencilView& index_weights,
                     bool is_entrance,
                     sasktran2::Dual<double, sasktran2::dualstorage::dense,
                                     NSTOKES>& source) const;
    };

    template <int NSTOKES>
    inline void scattering_source(
        const PhaseHandler<NSTOKES>& phase_handler, int threadidx, int losidx,
        int layeridx, int wavelidx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance, double solar_trans,
        const atmosphere::Atmosphere<NSTOKES>& atmosphere,
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
            solar_trans_iter,
        bool calculate_derivatives,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            source) {
        const auto& storage = atmosphere.storage();
        source.value.setZero();

        if (calculate_derivatives) {
            source.deriv.setZero();
        }

        double ssa = 0;
        double k = 0;
        if (index_weights.size() == 2) {
            // Structured 1D layers always have two backing nodes. Keep the
            // common ray representation while retaining the direct endpoint
            // interpolation used by the previous 1D-specific stencil.
            const auto lower = index_weights[0];
            const auto upper = index_weights[1];
            if (lower.second == 0.0) {
                ssa = storage.ssa(upper.first, wavelidx) * upper.second;
                k = storage.total_extinction(upper.first, wavelidx) *
                    upper.second;
            } else if (upper.second == 0.0) {
                ssa = storage.ssa(lower.first, wavelidx) * lower.second;
                k = storage.total_extinction(lower.first, wavelidx) *
                    lower.second;
            } else {
                ssa = storage.ssa(lower.first, wavelidx) * lower.second +
                      storage.ssa(upper.first, wavelidx) * upper.second;
                k = storage.total_extinction(lower.first, wavelidx) *
                        lower.second +
                    storage.total_extinction(upper.first, wavelidx) *
                        upper.second;
            }
        } else {
            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto ele = index_weights[index];
                if (ele.second == 0.0) {
                    continue;
                }
                ssa += storage.ssa(ele.first, wavelidx) * ele.second;
                k += storage.total_extinction(ele.first, wavelidx) * ele.second;
            }
        }

        const double source_amplitude = k * ssa * solar_trans / (EIGEN_PI * 4);
        const bool use_zero_safe_derivative_path =
            calculate_derivatives && (k == 0.0 || ssa == 0.0);
        if (!use_zero_safe_derivative_path) {
            // The common path lets the phase handler work directly with the
            // scaled source, avoiding a second dense derivative pass.
            source.value(0) = source_amplitude;
            phase_handler.scatter(threadidx, losidx, layeridx, index_weights,
                                  is_entrance, source);

            if (!calculate_derivatives) {
                return;
            }
            for (auto it = solar_trans_iter; it; ++it) {
                source.deriv(Eigen::placeholders::all, it.index()) -=
                    it.value() * source.value;
            }
            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto ele = index_weights[index];
                if (ele.second == 0.0) {
                    continue;
                }
                source.deriv(Eigen::placeholders::all,
                             atmosphere.ssa_deriv_start_index() + ele.first) +=
                    ele.second * source.value / ssa;
                source.deriv(Eigen::placeholders::all, ele.first) +=
                    ele.second * source.value / k;
            }
            return;
        }

        // For derivatives at k == 0 or SSA == 0, evaluate the unit-amplitude
        // phase first so the nonzero boundary derivative remains well-defined.
        source.value(0) = 1.0;

        phase_handler.scatter(threadidx, losidx, layeridx, index_weights,
                              is_entrance, source);

        const Eigen::Vector<double, NSTOKES> phase = source.value;
        source.value *= source_amplitude;

        if (!calculate_derivatives) {
            return;
        }
        source.deriv *= source_amplitude;
        // Solar transmission derivative factors
        for (auto it = solar_trans_iter; it; ++it) {
            source.deriv(Eigen::placeholders::all, it.index()) -=
                it.value() * source.value;
        }

        // And SSA/k derivative factors
        for (std::size_t index = 0; index < index_weights.size(); ++index) {
            const auto ele = index_weights[index];
            if (ele.second == 0.0) {
                continue;
            }
            source.deriv(Eigen::placeholders::all,
                         atmosphere.ssa_deriv_start_index() + ele.first) +=
                ele.second * k * solar_trans / (EIGEN_PI * 4) * phase;

            source.deriv(Eigen::placeholders::all, ele.first) +=
                ele.second * ssa * solar_trans / (EIGEN_PI * 4) * phase;
        }
    }

    /** Evaluates one exact single-scatter endpoint and accumulates its
     * derivatives directly into the integrated ray source. This avoids
     * materializing and then copying a dense endpoint derivative buffer. */
    template <int NSTOKES>
    inline Eigen::Vector<double, NSTOKES> accumulate_exact_scattering_source(
        const PhaseHandler<NSTOKES>& phase_handler, int threadidx, int losidx,
        int layeridx, int wavelidx,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance, double solar_trans,
        const atmosphere::Atmosphere<NSTOKES>& atmosphere,
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
            solar_trans_iter,
        double derivative_scale,
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
            target) {
        const auto& storage = atmosphere.storage();
        double ssa = 0.0;
        double extinction = 0.0;
        for (std::size_t index = 0; index < index_weights.size(); ++index) {
            const auto weight = index_weights[index];
            if (weight.second == 0.0) {
                continue;
            }
            ssa += storage.ssa(weight.first, wavelidx) * weight.second;
            extinction += storage.total_extinction(weight.first, wavelidx) *
                          weight.second;
        }

        const double unscaled_amplitude = solar_trans / (EIGEN_PI * 4);
        const double source_amplitude = extinction * ssa * unscaled_amplitude;
        const Eigen::Vector<double, NSTOKES> phase =
            phase_handler.scatter_and_accumulate_derivative(
                threadidx, losidx, layeridx, index_weights, is_entrance,
                source_amplitude, derivative_scale, target);
        const Eigen::Vector<double, NSTOKES> endpoint_source =
            source_amplitude * phase;

        for (auto it = solar_trans_iter; it; ++it) {
            target.deriv.col(it.index()) -=
                derivative_scale * it.value() * endpoint_source;
        }

        for (std::size_t index = 0; index < index_weights.size(); ++index) {
            const auto weight = index_weights[index];
            if (weight.second == 0.0) {
                continue;
            }
            target.deriv.col(atmosphere.ssa_deriv_start_index() +
                             weight.first) += derivative_scale * weight.second *
                                              extinction * unscaled_amplitude *
                                              phase;
            target.deriv.col(weight.first) += derivative_scale * weight.second *
                                              ssa * unscaled_amplitude * phase;
        }

        return endpoint_source;
    }

    template <int NSTOKES> struct ExactScatteringBlockScratch {
        using BatchMatrix =
            Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>;

        Eigen::RowVectorXd ssa;
        Eigen::RowVectorXd extinction;
        Eigen::RowVectorXd unscaled_amplitude;
        Eigen::RowVectorXd source_amplitude;
        Eigen::RowVectorXd ssa_factor;
        Eigen::RowVectorXd extinction_factor;
        BatchMatrix phase;
        BatchMatrix endpoint_source;

        void resize(int capacity) {
            ssa.resize(capacity);
            extinction.resize(capacity);
            unscaled_amplitude.resize(capacity);
            source_amplitude.resize(capacity);
            ssa_factor.resize(capacity);
            extinction_factor.resize(capacity);
            phase.resize(NSTOKES, capacity);
            endpoint_source.resize(NSTOKES, capacity);
        }
    };

    template <int NSTOKES> struct ExactIntegrationBlockScratch {
        using BatchMatrix =
            Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>;

        Eigen::RowVectorXd source_factor;
        Eigen::RowVectorXd source_factor_derivative;
        Eigen::RowVectorXd start_derivative_scale;
        Eigen::RowVectorXd end_derivative_scale;
        BatchMatrix integrated_value;
        BatchMatrix endpoint_quadrature;

        void resize(int capacity) {
            source_factor.resize(capacity);
            source_factor_derivative.resize(capacity);
            start_derivative_scale.resize(capacity);
            end_derivative_scale.resize(capacity);
            integrated_value.resize(NSTOKES, capacity);
            endpoint_quadrature.resize(NSTOKES, capacity);
        }
    };

    /** Batch equivalent of accumulate_exact_scattering_source. Values for a
     * fixed atmospheric/derivative coordinate are contiguous in wavelength so
     * Eigen can vectorize the common arithmetic. */
    template <int NSTOKES>
    inline void accumulate_exact_scattering_source_block(
        const PhaseHandler<NSTOKES>& phase_handler, int threadidx, int losidx,
        int layeridx, const sasktran2::WavelengthBlock& batch,
        const raytracing::GridWeightStencilView& index_weights,
        bool is_entrance,
        const Eigen::Ref<const Eigen::RowVectorXd>& solar_trans,
        const atmosphere::Atmosphere<NSTOKES>& atmosphere,
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator
            solar_trans_iter,
        const Eigen::Ref<const Eigen::RowVectorXd>& derivative_scale,
        sasktran2::WavelengthBlockDual<NSTOKES>& target,
        ExactScatteringBlockScratch<NSTOKES>& scratch) {
        const auto& storage = atmosphere.storage();
        auto ssa = scratch.ssa.head(batch.count);
        auto extinction = scratch.extinction.head(batch.count);
        ssa.setZero();
        extinction.setZero();
        for (std::size_t index = 0; index < index_weights.size(); ++index) {
            const auto weight = index_weights[index];
            if (weight.second == 0.0) {
                continue;
            }
            for (int lane = 0; lane < batch.count; ++lane) {
                const int wavelength = batch.wavelength(lane);
                ssa(lane) +=
                    storage.ssa(weight.first, wavelength) * weight.second;
                extinction(lane) +=
                    storage.total_extinction(weight.first, wavelength) *
                    weight.second;
            }
        }

        auto unscaled_amplitude = scratch.unscaled_amplitude.head(batch.count);
        auto source_amplitude = scratch.source_amplitude.head(batch.count);
        unscaled_amplitude.array() = solar_trans.array() / (EIGEN_PI * 4);
        source_amplitude.array() =
            extinction.array() * ssa.array() * unscaled_amplitude.array();
        phase_handler.scatter_and_accumulate_derivative_block(
            threadidx, losidx, layeridx, index_weights, is_entrance,
            source_amplitude, derivative_scale, target, scratch.phase);
        auto phase = scratch.phase.leftCols(batch.count);
        auto endpoint_source = scratch.endpoint_source.leftCols(batch.count);
        endpoint_source = phase;
        endpoint_source.array().rowwise() *= source_amplitude.array();

        if (target.derivative_size() > 0) {
            for (auto derivative = solar_trans_iter; derivative; ++derivative) {
                auto target_derivative =
                    target.derivative(derivative.index(), batch.count);
                target_derivative.array() -=
                    (endpoint_source.array().rowwise() *
                     derivative_scale.array()) *
                    derivative.value();
            }

            for (std::size_t index = 0; index < index_weights.size(); ++index) {
                const auto weight = index_weights[index];
                if (weight.second == 0.0) {
                    continue;
                }
                auto ssa_derivative = target.derivative(
                    atmosphere.ssa_deriv_start_index() + weight.first,
                    batch.count);
                auto ssa_factor = scratch.ssa_factor.head(batch.count);
                ssa_factor.array() = weight.second * derivative_scale.array() *
                                     extinction.array() *
                                     unscaled_amplitude.array();
                ssa_derivative.array() +=
                    phase.array().rowwise() * ssa_factor.array();

                auto extinction_derivative =
                    target.derivative(weight.first, batch.count);
                auto extinction_factor =
                    scratch.extinction_factor.head(batch.count);
                extinction_factor.array() =
                    weight.second * derivative_scale.array() * ssa.array() *
                    unscaled_amplitude.array();
                extinction_derivative.array() +=
                    phase.array().rowwise() * extinction_factor.array();
            }
        }
    }

    template <typename S, int NSTOKES>
    class SingleScatterSource : public SourceTermInterface<NSTOKES> {
      private:
        S m_solar_transmission;
        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere;

        Eigen::MatrixXd m_geometry_matrix;
        Eigen::SparseMatrix<double, Eigen::RowMajor> m_geometry_sparse;
        std::vector<bool> m_ground_hit_flag;

        std::vector<Eigen::VectorXd> m_solar_trans;
        using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>;
        std::vector<RowMajorMatrix> m_solar_trans_batch;
        int m_wavelength_batch_capacity = 0;
        std::vector<std::vector<int>> m_index_map;

        PhaseHandler<NSTOKES> m_phase_handler;

        // [los][layer][solar endpoint][geometry endpoint], where endpoint 0 is
        // the layer exit and endpoint 1 is the layer entrance.
        using ActiveDerivativeIndices =
            std::array<std::array<std::vector<int>, 2>, 2>;
        std::vector<std::vector<ActiveDerivativeIndices>>
            m_active_derivative_indices;
        const std::vector<sasktran2::raytracing::TracedRay>* m_traced_rays =
            nullptr;

        mutable std::vector<
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>>
            m_start_source_cache;
        mutable std::vector<
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>>
            m_end_source_cache;
        mutable std::vector<std::array<ExactScatteringBlockScratch<NSTOKES>, 2>>
            m_batch_source_cache;
        mutable std::vector<ExactIntegrationBlockScratch<NSTOKES>>
            m_batch_integration_cache;

        const Geometry& m_geometry;
        const Geometry1D* m_geometry_1d = nullptr;
        const Geometry2D* m_geometry_2d = nullptr;
        const sasktran2::Config* m_config;

        std::vector<bool> m_los_ground_is_hit;
        std::vector<sasktran2::raytracing::LayerGeometry> m_los_end_layers;

        void initialize_active_derivative_indices();

        void integrated_source_constant(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::LayerGeometry& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            sasktran2::WavelengthBlockLaneDualView<NSTOKES>& source,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection
                direction) const;

      public:
        SingleScatterSource(
            const Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer)
            : m_solar_transmission(geometry, raytracer), m_geometry(geometry),
              m_geometry_1d(&geometry), m_phase_handler(geometry){};

#ifdef SKTRAN_RUST_SUPPORT
        template <typename T = S,
                  std::enable_if_t<std::is_same_v<T, SolarTransmissionExact>,
                                   int> = 0>
        SingleScatterSource(
            const Geometry2D& geometry,
            const sasktran2::raytracing::RustRayTracer2D& raytracer)
            : m_solar_transmission(geometry, raytracer), m_geometry(geometry),
              m_geometry_2d(&geometry), m_phase_handler(geometry){};
#endif

        void initialize_config(const sasktran2::Config& config) override;

        /** Here the single scatter source term initializes the internal solar
         * transmission object, usually this involves tracing the required rays
         * and setting up any internal matrices.
         *
         *  @param internal_viewing Information on the internal viewing
         * geometry, los_rays and flux observers
         */
        void initialize_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) override;

        /**
         *
         */
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override;

        /** Triggers an internal calculation of the source term.  This method is
         * called at the beginning of each 'wavelength' calculation.
         *
         * @param wavelidx Index of the wavelength being calculated
         */
        int maximum_wavelength_block_size() const override {
            return std::is_same_v<S, SolarTransmissionExact>
                       ? std::numeric_limits<int>::max()
                       : 1;
        }

      private:
        void calculate_single(int wavelidx, int threadidx);

        void initialize_wavelength_blocks(int block_size);

        void calculate_block(const sasktran2::WavelengthBlock& block,
                             int threadidx);

      public:
        void calculate(const sasktran2::WavelengthBlock& block,
                       int threadidx) override {
            if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
                calculate_block(block, threadidx);
            } else {
                calculate_single(block.start, threadidx);
            }
        }

        /** Calculates the integrated source term for a given layer.
         *
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param layeridx Raw index pointing to the layer that was previosuly
         * passed in initialize_geometry
         * @param layer The layer that we are integrating over
         * @param source The returned source term
         */
      private:
        void integrated_source_single(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            sasktran2::WavelengthBlockLaneDualView<NSTOKES>& source,
            typename SourceTermInterface<
                NSTOKES>::IntegrationDirection direction =
                SourceTermInterface<NSTOKES>::IntegrationDirection::none) const;

        void integrated_source_block(
            const sasktran2::WavelengthBlock& batch, int losidx, int layeridx,
            int wavel_threadidx, int threadidx,
            const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            sasktran2::WavelengthBlockDual<NSTOKES>& source,
            typename SourceTermInterface<
                NSTOKES>::IntegrationDirection direction =
                SourceTermInterface<NSTOKES>::IntegrationDirection::none) const;

      public:
        void integrated_source(
            const sasktran2::WavelengthBlock& block, int losidx, int layeridx,
            int wavel_threadidx, int threadidx,
            const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            sasktran2::WavelengthBlockDual<NSTOKES>& source,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection
                direction =
                    SourceTermInterface<NSTOKES>::IntegrationDirection::none)
            const override {
            if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
                integrated_source_block(block, losidx, layeridx,
                                        wavel_threadidx, threadidx, layer,
                                        entrance_weights, exit_weights,
                                        shell_od, source, direction);
            } else {
                sasktran2::WavelengthBlockLaneDualView<NSTOKES> source_lane(
                    source, 0);
                integrated_source_single(block.start, losidx, layeridx,
                                         wavel_threadidx, threadidx, layer,
                                         entrance_weights, exit_weights,
                                         shell_od, source_lane, direction);
            }
        }

        bool supports_geometry_dimension(int dimension) const override {
            return dimension == 1 ||
                   (dimension == 2 && m_geometry_2d != nullptr);
        }

        bool supports_sparse_derivative_tracking() const override {
            return std::is_same_v<S, SolarTransmissionExact>;
        }

        void append_end_of_ray_active_derivatives(
            int losidx, std::vector<int>& derivative_indices) const override;

        void append_interior_active_derivatives(
            int losidx, int layeridx,
            std::vector<int>& derivative_indices) const override;

        /** Calculates the source term at the end of the ray.  Common examples
         * of this are ground scattering, ground emission, or the solar radiance
         * if looking directly at the sun.
         *
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param surface The surface object
         * @param source The returned source term
         */
      private:
        void end_of_ray_source_single(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockLaneDualView<NSTOKES>& source) const;

        void end_of_ray_source_block(
            const sasktran2::WavelengthBlock& batch, int losidx,
            int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockDual<NSTOKES>& source) const;

      public:
        void end_of_ray_source(
            const sasktran2::WavelengthBlock& block, int losidx,
            int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockDual<NSTOKES>& source) const override {
            if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
                end_of_ray_source_block(block, losidx, wavel_threadidx,
                                        threadidx, source);
            } else {
                sasktran2::WavelengthBlockLaneDualView<NSTOKES> source_lane(
                    source, 0);
                end_of_ray_source_single(block.start, losidx, wavel_threadidx,
                                         threadidx, source_lane);
            }
        }

        /**
         * @brief Not used for the Single Scatter source.
         *
         * @param wavelidx
         * @param losidx
         * @param wavel_threadidx
         * @param threadidx
         * @param source
         */
        void start_of_ray_source(
            const sasktran2::WavelengthBlock&, int, int, int,
            sasktran2::WavelengthBlockDual<NSTOKES>&) const override {}
    };

    template <int NSTOKES>
    class OccultationSource : public SourceTermInterface<NSTOKES> {
      private:
        std::vector<bool> m_ground_is_hit;

      public:
        void initialize_config(const sasktran2::Config& config) override;

        /** Initializes any geometry information that is required for
         * calculating the source term.  This method is called after the line of
         * sight rays ar traced.
         *
         * @param internal_viewing Information on the internal viewing geometry,
         * los_rays and flux observers
         */
        void initialize_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) override;

        /**
         *
         */
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override;

        /** Triggers an internal calculation of the source term.  This method is
         * called at the beginning of each 'wavelength' calculation.
         *
         * @param wavelidx Index of the wavelength being calculated
         */
        int maximum_wavelength_block_size() const override {
            return std::numeric_limits<int>::max();
        }

        void calculate(const sasktran2::WavelengthBlock&, int) override {}

        /** Calculates the integrated source term for a given layer.
         *
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param layeridx Raw index pointing to the layer that was previosuly
         * passed in initialize_geometry
         * @param layer The layer that we are integrating over
         * @param source The returned source term
         */
        void integrated_source(
            const sasktran2::WavelengthBlock&, int, int, int, int,
            const sasktran2::raytracing::TracedLayer&,
            const sasktran2::raytracing::GridWeightStencilView&,
            const sasktran2::raytracing::GridWeightStencilView&,
            const sasktran2::WavelengthBlockODView&,
            sasktran2::WavelengthBlockDual<NSTOKES>&,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection =
                SourceTermInterface<NSTOKES>::IntegrationDirection::none)
            const override {}

        /** Calculates the source term at the end of the ray.  Common examples
         * of this are ground scattering, ground emission, or the solar radiance
         * if looking directly at the sun.
         *
         * @param wavelidx Raw index for the wavelength we are calculating
         * @param losidx Raw index pointing to the ray that was previously
         * passed in initialize_geometry
         * @param source The returned source term
         */
        void end_of_ray_source(
            const sasktran2::WavelengthBlock& block, int losidx,
            int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockDual<NSTOKES>& source) const override;

        /**
         * @brief Not used for the occultation source.
         *
         * @param wavelidx
         * @param losidx
         * @param wavel_threadidx
         * @param threadidx
         * @param source
         */
        void start_of_ray_source(
            const sasktran2::WavelengthBlock&, int, int, int,
            sasktran2::WavelengthBlockDual<NSTOKES>&) const override {}

        bool has_interior_source() const override { return false; }
    };

} // namespace sasktran2::solartransmission
