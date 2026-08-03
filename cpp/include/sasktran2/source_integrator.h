#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/source_interface.h>
#include <algorithm>

namespace sasktran2 {

    /** Geometry-only transport stencils packed for one-time transfer to the
     *  Rust successive-orders implementation. */
    struct RayTransportGeometry {
        std::vector<std::uint32_t> ray_layer_offsets;
        std::vector<std::uint32_t> layer_atmosphere_offsets;
        std::vector<std::uint32_t> atmosphere_indices;
        std::vector<double> optical_depth_weights;
        std::vector<double> albedo_weights;
        std::vector<double> entrance_weights;
        std::vector<double> exit_weights;
        std::vector<double> layer_distance;
        std::vector<double> layer_start_fraction;
        std::vector<double> layer_end_fraction;
        std::vector<double> ray_scattering_cosine;
        std::vector<double> ray_phase_q_factor;
        std::vector<double> ray_phase_u_factor;
        std::vector<std::uint32_t> ray_transport_value_offsets;
        std::vector<std::uint32_t> ray_transport_row_nnz;
        std::vector<std::uint32_t> layer_source_offsets;
        std::vector<std::uint16_t> source_value_inner_indices;
        std::vector<double> source_weights;
        std::vector<std::uint32_t> ray_ground_offsets;
        std::vector<std::uint16_t> ground_value_inner_indices;
        std::vector<double> ground_weights;
    };

    template <int NSTOKES> struct RaySourceInterpolationWeights {
        struct Layer {
            std::uint32_t atmosphere_offset = 0;
            std::uint32_t source_offset = 0;
            std::uint32_t atmosphere_count = 0;
            std::uint32_t source_count = 0;
        };

        std::vector<Layer> interior_weights;

        // Structure-of-arrays storage avoids tuple padding and replaces two
        // heap allocations per traced layer with a few allocations per ray.
        std::vector<double> atmosphere_weights;
        std::vector<int> source_indices;
        std::vector<double> source_weights;
        std::vector<std::uint16_t> source_accumulation_inner_indices;

        std::vector<int> ground_source_indices;
        std::vector<double> ground_source_weights;
        std::vector<std::uint16_t> ground_accumulation_inner_indices;
        std::uint32_t accumulation_row_offset = 0;
        std::uint32_t accumulation_row_nnz = 0;
        bool ground_is_hit = false;

        std::size_t source_accumulation_index(std::size_t source_index,
                                              int stokes) const {
            return static_cast<std::size_t>(accumulation_row_offset) +
                   static_cast<std::size_t>(stokes) * accumulation_row_nnz +
                   source_accumulation_inner_indices[source_index];
        }

        std::size_t ground_accumulation_index(std::size_t source_index,
                                              int stokes) const {
            return static_cast<std::size_t>(accumulation_row_offset) +
                   static_cast<std::size_t>(stokes) * accumulation_row_nnz +
                   ground_accumulation_inner_indices[source_index];
        }
    };

    /** Class that integrates source terms along the ray.  Note that in
     * SASKTRAN2, source terms themselves are responsible for integrating across
     * the layer, this class simply adds the source terms in each layer and
     * attenuates them by the optical depth.
     *
     *  Integration takes place in three steps.  First initialize_geometry is
     * called with the rays that will be eventually integrated to set up
     * geometry factors.
     *
     *  Next, initialize_atmosphere is called so that any optical parameters can
     * be pre-calculated, such as the OD for each layer.
     *
     *  Lastly, integrate is called on each ray individually, summing the
     * overall sources together.
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> class SourceIntegrator {
        using SInterpolator =
            std::vector<RaySourceInterpolationWeights<NSTOKES>>;

      private:
        bool m_derivatives_enabled; /**< Whether this integrator was configured
                                       to calculate derivatives */
        bool m_calculate_derivatives; /**< True if we are calculating
                                         derivatives */
        std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>>
            m_traced_ray_od_matrix; /**< Vector of matrices A such that A *
                                       atmosphere_extinction = OD for each layer
                                       in that ray */
        Eigen::SparseMatrix<double, Eigen::RowMajor>
            m_empty_od_matrix; /**< Derivative-free matrix view used when
                                  optical depth is evaluated on demand. */
        bool m_on_demand_optical_depth = false;

        using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>;
        std::vector<RowMajorMatrix>
            m_shell_od; /**< Optical depth for every ray, layer, and
                           wavelength. */
        std::vector<Eigen::MatrixXd>
            m_scalar_shell_od; /**< Column-major optical depth storage used
                                  when the configured block capacity is one. */
        int m_wavelength_block_capacity = 1;
        mutable std::vector<Eigen::RowVectorXd> m_thread_attenuation{1};
        struct AccumulationThreadScratch {
            std::vector<double> shell_od_cotangent;
            std::vector<double> current_od_cotangent;
            Eigen::Vector<double, NSTOKES> layer_source;
            sasktran2::WavelengthBlockDual<NSTOKES> primal_layer_source;
            sasktran2::WavelengthBlockDual<NSTOKES> end_source;
        };
        mutable std::vector<AccumulationThreadScratch>
            m_accumulation_thread_scratch{1};
        const std::vector<sasktran2::raytracing::TracedRay>* m_traced_rays =
            nullptr; /**< Reference to the rays we are integrating */

        struct CompactLayer {
            double layer_distance;
            double od_quad_start;
            double od_quad_end;
            double od_quad_start_fraction;
            double od_quad_end_fraction;
            std::uint32_t grid_weight_offset;
            std::uint8_t grid_weight_count;
        };
        struct CompactGridWeights {
            std::vector<int> indices;
            std::vector<double> entrance;
            std::vector<double> exit;
            std::vector<double> optical_depth;
        };
        struct CompactRay {
            std::vector<CompactLayer> layers;
            CompactGridWeights weights;
        };
        bool m_use_compact_geometry = false;
        bool m_interior_geometry_released = false;
        std::vector<CompactRay> m_compact_rays;
        std::size_t m_compact_max_layers = 0;

        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere =
            nullptr;
        int m_num_geometry_locations = 0;
        int m_num_geometry_dimensions = 1;
        bool m_use_sparse_derivative_tracking = false;
        std::vector<std::vector<std::vector<std::pair<int, int>>>>
            m_attenuation_active_derivative_ranges;

        int integration_num_layers(int rayidx) const;
        const sasktran2::raytracing::TracedLayer&
        integration_layer(int rayidx, int layeridx,
                          sasktran2::raytracing::TracedLayer& scratch) const;
        sasktran2::raytracing::GridWeightStencilView
        integration_entrance_weights(int rayidx, int layeridx) const;
        sasktran2::raytracing::GridWeightStencilView
        integration_exit_weights(int rayidx, int layeridx) const;
        sasktran2::raytracing::GridWeightStencilView
        integration_optical_depth_weights(int rayidx, int layeridx) const;

        template <int N>
        void integrate_block(
            sasktran2::WavelengthBlockDual<NSTOKES>& radiance,
            const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
            const sasktran2::WavelengthBlock<N>& block, int rayidx,
            int wavel_threadidx, int threadidx) const;

        template <int N, typename ShellODMatrix>
        void integrate_ray(
            sasktran2::WavelengthBlockDual<NSTOKES>& radiance,
            const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
            const sasktran2::raytracing::TracedRay& ray,
            const Eigen::SparseMatrix<double, Eigen::RowMajor>& od_matrix,
            const ShellODMatrix& shell_od,
            const sasktran2::WavelengthBlock<N>& batch, int rayidx,
            int wavel_threadidx, int threadidx) const;

      public:
        /**
         *
         * @param calculate_derivatives True if the source integrator should
         * calculate derivatives
         */
        SourceIntegrator(bool calculate_derivatives);

        /** Avoid persistent optical-depth matrices and wavelength caches.
         *  Intended for derivative-free callers that use
         *  integrate_and_emplace_accumulation_triplets. */
        void set_on_demand_optical_depth(bool enable) {
            m_on_demand_optical_depth = enable;
        }

        /**
         *
         * @param enable True if the source integrator should calculate
         * derivatives
         */
        void set_calculate_derivatives(bool enable) {
            m_derivatives_enabled = enable;
            m_calculate_derivatives = enable;
            m_use_sparse_derivative_tracking = false;
            m_attenuation_active_derivative_ranges.clear();
        }

        /** Initializes the geometry of the source integrator
         *
         * @param traced_rays Vector of traced rays
         * @param geometry Global geometry
         */
        void initialize_geometry(
            const std::vector<sasktran2::raytracing::TracedRay>& traced_rays,
            const Geometry& geometry);

        /** Starts incremental compact 2D geometry construction. */
        void begin_compact_geometry_2d(
            std::vector<sasktran2::raytracing::TracedRay>& traced_rays,
            const Geometry& geometry);

        /** Moves one freshly traced 2D ray into compact integration storage,
         *  retaining only its first layer for end-of-ray source handling. */
        void
        compact_geometry_2d_ray(int rayidx,
                                sasktran2::raytracing::TracedRay& traced_ray);

        /** Completes incremental compact 2D geometry construction. */
        void finalize_compact_geometry_2d();

        /** Releases layer geometry after it has been transferred to an
         * external transport implementation. End-of-ray source helpers and
         * atmosphere initialization remain available. */
        void release_interior_geometry();

        /** Initializes the atmosphere
         *
         * @param atmo
         */
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmo);

        /** Configures the active wavelength-block capacity and allocates
         * reusable thread scratch. This must precede initialize_atmosphere so
         * shell optical depth uses the matching scalar or batched layout. */
        void initialize_thread_storage(int num_threads, int block_capacity) {
            m_wavelength_block_capacity = block_capacity;
            m_thread_attenuation.resize(num_threads);
            m_accumulation_thread_scratch.resize(num_threads);
            for (auto& attenuation : m_thread_attenuation) {
                attenuation.resize(block_capacity);
            }
            for (auto& scratch : m_accumulation_thread_scratch) {
                scratch.primal_layer_source.resize(1, 0, true);
                scratch.end_source.resize(1, 0, true);
            }
        }

        /** Precomputes the cumulative derivative columns that can be nonzero
         * before each layer attenuation. This is enabled only when every
         * active source can report its derivative sparsity exactly. */
        void initialize_derivative_sparsity(
            const std::vector<SourceTermInterface<NSTOKES>*>& source_terms);

        /** Returns true when the integrator and every active source implement
         * the requested native derivative execution mode. */
        bool
        supports_linearization(sasktran2::LinearizationMode mode,
                               const std::vector<SourceTermInterface<NSTOKES>*>&
                                   source_terms) const {
            return std::all_of(
                source_terms.begin(), source_terms.end(),
                [mode](const SourceTermInterface<NSTOKES>* source) {
                    return source->supports_linearization(mode);
                });
        }

        void integrate_jvp(
            sasktran2::RadianceJVP<NSTOKES>& radiance,
            const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
            int wavelength, int rayidx, int wavel_threadidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent) const;

        void integrate_value(
            Eigen::Vector<double, NSTOKES>& radiance,
            const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
            int wavelength, int rayidx, int wavel_threadidx,
            int threadidx) const;

        void integrate_vjp(
            Eigen::Vector<double, NSTOKES>& radiance,
            const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
            int wavelength, int rayidx, int wavel_threadidx, int threadidx,
            const Eigen::Vector<double, NSTOKES>& cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const;

        /** Integrates the source terms and stores the result in radiance
         *
         * @param radiance
         * @param source_terms
         * @param wavelidx
         * @param wavel_threadidx
         * @param threadidx
         * @param rayidx
         */
        void integrate(
            sasktran2::WavelengthBlockDual<NSTOKES>& radiance,
            const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
            const sasktran2::WavelengthBlock<>& block, int rayidx,
            int wavel_threadidx, int threadidx);

        void integrate_and_emplace_accumulation_triplets(
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                radiance,
            const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
            int wavelidx, int rayidx, int wavel_threadidx, int threadidx,
            const SInterpolator& source_interpolator,
            Eigen::VectorXd& accumulation_values);

        /** Compiles traced geometry and source interpolation into flat,
         *  wavelength-independent Rust transport stencils. */
        RayTransportGeometry
        pack_ray_transport_geometry(const SInterpolator& source_interpolator,
                                    const Geometry& geometry) const;

        /** Evaluates only the exact first-order source using optical depth and
         *  attenuation already assembled by Rust. */
        void integrate_first_order_precomputed(
            sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>&
                radiance,
            const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
            int wavelidx, int rayidx, int wavel_threadidx, int threadidx,
            std::uint32_t flat_layer_offset,
            const std::vector<double>& layer_optical_depth,
            const std::vector<double>& layer_attenuation,
            const std::vector<double>& layer_prefix_attenuation,
            double ray_end_attenuation);

        /** Rebuilds only the compact transport operator for one diffuse ray.
         * This omits source evaluation and radiance bookkeeping when restoring
         * a converged forward checkpoint. */
        void emplace_accumulation_transport(
            int wavelidx, int rayidx, const SInterpolator& source_interpolator,
            Eigen::Ref<Eigen::VectorXd> accumulation_values) const;

        /** Differentiates the compact transport values and first-order
         * forcing assembled for one internal diffuse ray. */
        void integrate_and_emplace_accumulation_jvp(
            const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
            int wavelidx, int rayidx, int wavel_threadidx, int threadidx,
            const SInterpolator& source_interpolator,
            Eigen::Ref<const Eigen::VectorXd> native_tangent,
            Eigen::VectorXd* accumulation_values,
            Eigen::Ref<Eigen::VectorXd> accumulation_value_tangent,
            Eigen::Ref<Eigen::VectorXd> first_order_forcing_tangent) const;

        /** Pulls compact transport/forcing cotangents for one internal
         * diffuse ray back to native atmosphere parameters. */
        void accumulate_accumulation_vjp(
            const std::vector<SourceTermInterface<NSTOKES>*>& source_terms,
            int wavelidx, int rayidx, int wavel_threadidx, int threadidx,
            const SInterpolator& source_interpolator,
            Eigen::Ref<const Eigen::VectorXd> accumulation_value_gradient,
            Eigen::Ref<const Eigen::VectorXd> transport_solution,
            Eigen::Ref<const Eigen::VectorXi> transport_column_indices,
            Eigen::Ref<const Eigen::VectorXd> first_order_forcing_gradient,
            bool active_extinction, bool active_ssa,
            bool active_interior_source_parameters,
            bool active_transport_parameters,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const;

        /** Calculates the Optical Depth for each ray */
        void integrate_optical_depth(Eigen::MatrixXd& optical_depth);
    };
} // namespace sasktran2
