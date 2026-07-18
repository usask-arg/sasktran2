#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/source_interface.h>
#include <algorithm>

namespace sasktran2 {

    template <int NSTOKES> struct RaySourceInterpolationWeights {
        std::vector<std::pair<
            std::vector<std::pair<int, double>>,
            std::vector<std::tuple<int, double, std::array<int, NSTOKES>>>>>
            interior_weights;
        std::vector<std::tuple<int, double, std::array<int, NSTOKES>>>
            ground_weights;
        bool ground_is_hit;
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
        const std::vector<sasktran2::raytracing::TracedRay>* m_traced_rays =
            nullptr; /**< Reference to the rays we are integrating */

        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere =
            nullptr;
        int m_num_geometry_locations = 0;
        int m_num_geometry_dimensions = 1;
        bool m_use_sparse_derivative_tracking = false;
        std::vector<std::vector<std::vector<std::pair<int, int>>>>
            m_attenuation_active_derivative_ranges;

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
            for (auto& attenuation : m_thread_attenuation) {
                attenuation.resize(block_capacity);
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
            std::vector<SourceTermInterface<NSTOKES>*> source_terms,
            int wavelidx, int rayidx, int wavel_threadidx, int threadidx,
            const SInterpolator& source_interpolator,
            Eigen::VectorXd& accumulation_values);

        /** Calculates the Optical Depth for each ray */
        void integrate_optical_depth(Eigen::MatrixXd& optical_depth);
    };
} // namespace sasktran2
