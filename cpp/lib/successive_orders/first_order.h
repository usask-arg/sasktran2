#pragma once

#include "geometry.h"
#include "endpoint_stencil_storage.h"
#include "transport.h"

#include <sasktran2/solartransmission.h>
#include <sasktran2/source_integrator.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cstddef>
#include <memory>
#include <vector>

namespace sasktran2::successive_orders {

    /** Exact single-scatter illumination on the source's incoming rays.
     *
     * This is an implementation detail of the successive-orders source.  It
     * deliberately owns its source and integrator so no engine source needs to
     * be exposed or shared.  Ordinary integration is derivative-free; native
     * directional and reverse products use the exact source hooks directly.
     */
    template <int NSTOKES> class FirstOrderProvider {
      public:
        FirstOrderProvider(
            const sasktran2::Geometry1D& geometry,
            const sasktran2::raytracing::RayTracerBase& raytracer);
#ifdef SKTRAN_RUST_SUPPORT
        FirstOrderProvider(
            const sasktran2::Geometry2D& geometry,
            const sasktran2::raytracing::RustRayTracer2D& raytracer,
            std::shared_ptr<
                sasktran2::solartransmission::SolarTransmissionTable2D>
                shared_solar_table = nullptr);
#endif

        void initialize_config(const sasktran2::Config& config);
        void initialize_geometry(const SourceGeometry1D& source_geometry);
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
            bool volume_changed = true);

        /** Select the bounded spectral cache owned by this wavelength worker.
         *
         * Call before returning a cached primal or requesting native products.
         * Simultaneous workers must evaluate distinct wavelengths.
         */
        void prepare_wavelength(int wavelength, int wavelength_thread);

        int size() const { return m_num_rays * NSTOKES; }

        /** Calculates the first-order incoming radiance. */
        void calculate(int wavelength, int wavelength_thread,
                       Eigen::Ref<Eigen::VectorXd> forcing);

        bool uses_compact_scalar_kernel() const { return m_use_compact_scalar; }
        bool can_release_incoming_geometry() const {
            return m_use_compact_scalar && !m_use_lower_interpolation;
        }

        void calculate_with_transport(int wavelength, int wavelength_thread,
                                      TransportOperator& transport,
                                      Eigen::Ref<Eigen::VectorXd> forcing);

        /** Rebuilds transport and product layer caches without first-order
         * forcing. The caller must already own current, matching forcing. */
        void calculate_transport_only(int wavelength, int wavelength_thread,
                                      TransportOperator& transport);

        /** Calculates the first-order radiance and one native JVP. */
        void calculate_jvp(int wavelength, int wavelength_thread,
                           Eigen::Ref<const Eigen::VectorXd> native_tangent,
                           Eigen::Ref<Eigen::VectorXd> forcing,
                           Eigen::Ref<Eigen::VectorXd> forcing_tangent);

        void calculate_jvp_with_transport(
            int wavelength, int wavelength_thread,
            Eigen::Ref<const Eigen::VectorXd> native_tangent,
            const Eigen::VectorXd& layer_state_projection,
            const Eigen::VectorXd& ground_state_projection,
            Eigen::VectorXd& direct_transport_tangent,
            Eigen::Ref<Eigen::VectorXd> forcing,
            Eigen::Ref<Eigen::VectorXd> forcing_tangent);

        void project_transport_state(
            Eigen::Ref<const Eigen::VectorXd> transport_state,
            Eigen::VectorXd& layer_state_projection,
            Eigen::VectorXd& ground_state_projection) const;

        /** Accumulates the VJP of the first-order radiance. */
        void accumulate_vjp(int wavelength, int wavelength_thread,
                            Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
                            Eigen::Ref<Eigen::VectorXd> native_gradient);

        void accumulate_vjp_with_transport(
            int wavelength, int wavelength_thread,
            const Eigen::VectorXd& transport_state,
            Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient);

        void accumulate_vjp_with_projected_transport(
            int wavelength, int wavelength_thread,
            const Eigen::VectorXd& layer_state_projection,
            const Eigen::VectorXd& ground_state_projection,
            Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient);

        std::size_t workspace_bytes() const;

      private:
        using ExactSource = sasktran2::solartransmission::SingleScatterSource<
            sasktran2::solartransmission::SolarTransmissionExact, NSTOKES>;
        void validate_ready(int wavelength, int wavelength_thread) const;
        int ray_thread_index(int wavelength_thread) const;
        int scalar_cache_index(int wavelength) const;
        const Eigen::VectorXd& ensure_solar_transmission(int wavelength,
                                                         int wavelength_thread);
        void ensure_endpoint_medium(int wavelength);

        struct ScalarEndpoint {
            double extinction = 0.0;
            double albedo = 0.0;
            double phase = 0.0;
            double solar_transmission = 0.0;
            double source = 0.0;
        };

        struct ScalarValueTangent {
            double value = 0.0;
            double tangent = 0.0;
        };

        struct ScalarVjpScratch {
            Eigen::VectorXd optical_depth;
            Eigen::VectorXd attenuation;
            Eigen::VectorXd source_factor;
            Eigen::VectorXd prefix_attenuation;
            Eigen::VectorXd albedo;
            std::vector<ScalarEndpoint> endpoints;
            Eigen::VectorXd endpoint_cotangent;
        };

        struct ScalarLayerCache {
            Eigen::VectorXd optical_depth;
            Eigen::VectorXd source_factor;
            bool active = false;
        };

        struct ScalarEndpointMediumCache {
            Eigen::VectorXd extinction;
            Eigen::VectorXd albedo;
            bool active = false;
        };

        // These views are prepared after the worker's physical caches and
        // product scratch are sized. They live only for the current sweep.
        struct ScalarEndpointContext {
            const double* extinction = nullptr;
            const double* albedo = nullptr;
            const double* coefficients = nullptr;
            Eigen::Index coefficient_stride = 0;
            const int* phase_orders = nullptr;
            const int* maximum_orders = nullptr;
            const double* basis = nullptr;
            int endpoint_basis_stride = 0;
            const double* uniform_phase = nullptr;
            const int* endpoint_slots = nullptr;
            const double* endpoint_extinction = nullptr;
            const double* endpoint_albedo = nullptr;
            const double* solar = nullptr;
            int ssa_deriv_start = 0;
            const double* coefficient_tangent = nullptr;
            const int* tangent_orders = nullptr;
            const double* endpoint_extinction_tangent = nullptr;
            const double* endpoint_albedo_tangent = nullptr;
        };

        ScalarEndpointContext
        prepare_scalar_endpoint_context(int wavelength, int cache_slot) const;
        ScalarEndpointContext
        scalar_endpoint_context_for_ray(const ScalarEndpointContext& context,
                                        int ray) const;

        struct ScalarVolumeCache {
            // Retain only the exact layer sweep's transmission to the ground
            // for reverse products; transport and forcing live in the worker.
            Eigen::VectorXd ground_prefix;
            bool active = false;
        };

        struct ScalarPackedLayer {
            // Keep only owned scalar metadata here. Geometry stencil views are
            // reacquired from SourceGeometry1D so their backing storage and
            // bounds remain authoritative across repeated engine evaluations.
            double source_quad_start = 0.0;
            double source_quad_end = 0.0;
        };

        struct ScalarPackedRay {
            std::uint32_t layer_begin = 0;
            std::uint32_t layer_end = 0;
            std::int32_t ground_geometry = -1;
        };

        struct ScalarGroundGeometry {
            Eigen::Vector3d up;
            Eigen::Vector3d look_away;
        };

        void calculate_scalar(int wavelength, int wavelength_thread,
                              Eigen::Ref<Eigen::VectorXd> forcing,
                              TransportOperator* transport = nullptr);
        template <bool WITH_TRANSPORT, bool LOWER_INTERPOLATION,
                  bool WITH_FORCING = true>
        void calculate_scalar_impl(int wavelength, int wavelength_thread,
                                   Eigen::Ref<Eigen::VectorXd> forcing,
                                   TransportOperator* transport);
        template <bool WITH_TRANSPORT, bool WITH_FORCING = true>
        void calculate_scalar_uniform_impl(int wavelength,
                                           int wavelength_thread,
                                           Eigen::Ref<Eigen::VectorXd> forcing,
                                           TransportOperator* transport);
        void calculate_scalar_jvp(
            int wavelength, int wavelength_thread,
            Eigen::Ref<const Eigen::VectorXd> native_tangent,
            Eigen::Ref<Eigen::VectorXd> forcing,
            Eigen::Ref<Eigen::VectorXd> forcing_tangent,
            const Eigen::VectorXd* layer_state_projection = nullptr,
            const Eigen::VectorXd* ground_state_projection = nullptr,
            Eigen::VectorXd* direct_transport_tangent = nullptr);
        template <bool WITH_TRANSPORT, bool LOWER_INTERPOLATION>
        void calculate_scalar_jvp_impl(
            int wavelength, int wavelength_thread,
            Eigen::Ref<const Eigen::VectorXd> native_tangent,
            Eigen::Ref<Eigen::VectorXd> forcing,
            Eigen::Ref<Eigen::VectorXd> forcing_tangent,
            const Eigen::VectorXd* layer_state_projection,
            const Eigen::VectorXd* ground_state_projection,
            Eigen::VectorXd* direct_transport_tangent);
        void calculate_scalar_jvp_uniform_proportional(
            int wavelength, Eigen::Ref<const Eigen::VectorXd> native_tangent,
            const Eigen::VectorXd& solar_tangent,
            double extinction_direction_scale, double albedo,
            double albedo_tangent,
            const Eigen::VectorXd& layer_state_projection,
            const Eigen::VectorXd& ground_state_projection,
            Eigen::VectorXd& direct_transport_tangent,
            Eigen::Ref<Eigen::VectorXd> forcing_tangent);
        void accumulate_scalar_vjp(
            int wavelength, int wavelength_thread,
            Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient,
            const Eigen::VectorXd* transport_state,
            const Eigen::VectorXd* layer_state_projection,
            const Eigen::VectorXd* ground_state_projection);
        template <bool WITH_TRANSPORT, bool LOWER_INTERPOLATION>
        void accumulate_scalar_surface_vjp(
            int wavelength, int wavelength_thread,
            Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient,
            const Eigen::VectorXd* transport_state,
            const Eigen::VectorXd* ground_state_projection);
        template <bool WITH_TRANSPORT, bool LOWER_INTERPOLATION>
        void accumulate_scalar_vjp_dispatch(
            int wavelength, int wavelength_thread,
            Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient,
            const Eigen::VectorXd* transport_state,
            const Eigen::VectorXd* layer_state_projection,
            const Eigen::VectorXd* ground_state_projection);
        template <bool WITH_TRANSPORT, bool LOWER_INTERPOLATION,
                  bool WITH_PHASE_GRADIENT>
        void accumulate_scalar_vjp_impl(
            int wavelength, int wavelength_thread,
            Eigen::Ref<const Eigen::VectorXd> forcing_cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient,
            const Eigen::VectorXd* transport_state,
            const Eigen::VectorXd* layer_state_projection,
            const Eigen::VectorXd* ground_state_projection);

        // Callers resolve and validate their worker cache before entering the
        // ray loop; endpoint helpers reuse the prepared views in that sweep.
        template <typename Weights>
        ScalarEndpoint scalar_endpoint(const ScalarEndpointContext& context,
                                       int layer, bool entrance,
                                       int solar_index,
                                       const Weights& weights) const;
        template <bool USE_ENDPOINT_MEDIUM, typename Weights>
        ScalarValueTangent scalar_endpoint_jvp(
            const ScalarEndpointContext& context, int layer, bool entrance,
            int solar_index, const Weights& weights,
            const double* extinction_direction, const double* albedo_direction,
            const double* solar_tangent, bool uniform_albedo_direction,
            double uniform_albedo_tangent, bool phase_tangent_active) const;
        template <bool WITH_PHASE_GRADIENT, typename Weights>
        void accumulate_scalar_endpoint_vjp(
            const ScalarEndpointContext& context, int layer, bool entrance,
            int solar_index, const Weights& weights,
            const ScalarEndpoint& endpoint, double source_cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient,
            Eigen::Ref<Eigen::VectorXd> solar_gradient,
            Eigen::Ref<Eigen::VectorXd> coefficient_gradient) const;

        EIGEN_STRONG_INLINE ScalarEndpoint scalar_endpoint(
            const ScalarEndpointContext& context, int layer, bool entrance,
            int solar_index, const EndpointStencilView& weights) const {
            return weights.visit([&](const auto& values) {
                return scalar_endpoint(context, layer, entrance, solar_index,
                                       values);
            });
        }

        template <bool USE_ENDPOINT_MEDIUM>
        EIGEN_STRONG_INLINE ScalarValueTangent scalar_endpoint_jvp(
            const ScalarEndpointContext& context, int layer, bool entrance,
            int solar_index, const EndpointStencilView& weights,
            const double* extinction_direction, const double* albedo_direction,
            const double* solar_tangent, bool uniform_albedo_direction,
            double uniform_albedo_tangent, bool phase_tangent_active) const {
            return weights.visit([&](const auto& values) {
                return scalar_endpoint_jvp<USE_ENDPOINT_MEDIUM>(
                    context, layer, entrance, solar_index, values,
                    extinction_direction, albedo_direction, solar_tangent,
                    uniform_albedo_direction, uniform_albedo_tangent,
                    phase_tangent_active);
            });
        }

        template <bool WITH_PHASE_GRADIENT>
        EIGEN_STRONG_INLINE void accumulate_scalar_endpoint_vjp(
            const ScalarEndpointContext& context, int layer, bool entrance,
            int solar_index, const EndpointStencilView& weights,
            const ScalarEndpoint& endpoint, double source_cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient,
            Eigen::Ref<Eigen::VectorXd> solar_gradient,
            Eigen::Ref<Eigen::VectorXd> coefficient_gradient) const {
            weights.visit([&](const auto& values) {
                accumulate_scalar_endpoint_vjp<WITH_PHASE_GRADIENT>(
                    context, layer, entrance, solar_index, values, endpoint,
                    source_cotangent, native_gradient, solar_gradient,
                    coefficient_gradient);
            });
        }

        EndpointStencilView endpoint_weights(int solar_index) const {
            const int slot = m_endpoint_slots[solar_index];
            return m_endpoint_stencils.view(static_cast<std::size_t>(slot));
        }

        int phase_basis_slot(int ray, int solar_index) const {
            return m_endpoint_phase_basis ? solar_index : ray;
        }
        int num_phase_basis_slots() const {
            return m_endpoint_phase_basis ? m_solar_offsets.back() : m_num_rays;
        }
        bool ground_scattering_geometry(int solar_index,
                                        const ScalarGroundGeometry& ground,
                                        double& mu_in, double& mu_out,
                                        double& phi) const;
        double ground_transport_albedo(int wavelength, int ray) const;
        double ground_transport_albedo_tangent(
            int ray, Eigen::Ref<const Eigen::VectorXd> native_tangent) const;
        void accumulate_ground_transport_albedo_vjp(
            int ray, double albedo_cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const;

        const sasktran2::Geometry& m_geometry;
        const sasktran2::Geometry1D* m_geometry_1d = nullptr;
        ExactSource m_source;
        std::shared_ptr<
            sasktran2::solartransmission::SolarTransmissionTableEvaluator>
            m_solar_table;
        sasktran2::solartransmission::SolarTableInterpolation
            m_solar_interpolation;
        std::vector<bool> m_solar_ground_hit;
        std::vector<Eigen::Vector3d> m_solar_propagation_directions;
        std::vector<std::vector<std::pair<int, double>>>
            m_ground_horizontal_weights;
        sasktran2::SourceIntegrator<NSTOKES> m_integrator;
        std::vector<SourceTermInterface<NSTOKES>*> m_source_terms;
        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere =
            nullptr;
        const SourceGeometry1D* m_source_geometry = nullptr;
        int m_num_rays = 0;
        int m_num_phase_moments = 0;
        int m_num_threads = 1;
        int m_num_source_threads = 1;
        int m_num_wavelength_threads = 1;
        bool m_geometry_initialized = false;
        bool m_compact_scalar_requested = false;
        bool m_use_compact_scalar = false;
        bool m_use_lower_interpolation = false;
        bool m_solar_refraction = false;
        bool m_endpoint_phase_basis = false;

        std::vector<int> m_solar_offsets;
        std::vector<double> m_phase_basis;
        std::vector<int> m_endpoint_slots;
        EndpointStencilStorage m_endpoint_stencils;
        std::vector<ScalarPackedRay> m_scalar_packed_rays;
        std::vector<ScalarPackedLayer> m_scalar_packed_layers;
        std::vector<ScalarGroundGeometry> m_scalar_ground_geometry;
        mutable std::vector<sasktran2::WavelengthBlockDual<NSTOKES>>
            m_primal_scratch;
        mutable std::vector<Eigen::MatrixXd> m_gradient_scratch;
        mutable std::vector<
            Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>>
            m_vjp_radiance_scratch;
        mutable std::vector<
            Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>>
            m_vjp_cotangent_scratch;
        mutable std::vector<Eigen::VectorXd> m_solar_product_scratch;
        mutable std::vector<Eigen::VectorXd> m_solar_table_product_scratch;
        mutable std::vector<Eigen::VectorXd> m_phase_product_scratch;
        mutable std::vector<std::vector<int>> m_phase_order_scratch;
        mutable std::vector<Eigen::VectorXd>
            m_endpoint_extinction_tangent_scratch;
        mutable std::vector<Eigen::VectorXd> m_endpoint_albedo_tangent_scratch;
        std::vector<int> m_scalar_phase_orders;
        std::vector<unsigned char> m_uniform_phase_active;
        std::vector<std::vector<double>> m_uniform_phase_values;
        std::vector<unsigned char> m_uniform_albedo_active;
        std::vector<double> m_uniform_albedo_values;
        std::vector<Eigen::VectorXd> m_cached_solar_transmission;
        std::vector<unsigned char> m_cached_solar_active;
        std::vector<ScalarLayerCache> m_scalar_layer_cache;
        std::vector<ScalarEndpointMediumCache> m_endpoint_medium_cache;
        std::vector<ScalarVolumeCache> m_scalar_volume_cache;
        // Spectral values live in one reusable slot per wavelength worker.
        // Each active wavelength installs its own mapping before use; stale
        // mappings are deliberately left untouched by other workers.
        std::vector<int> m_scalar_cache_wavelength;
        std::vector<int> m_scalar_cache_index;
        mutable std::vector<ScalarVjpScratch> m_scalar_vjp_scratch;
    };

    extern template class FirstOrderProvider<1>;
    extern template class FirstOrderProvider<3>;

} // namespace sasktran2::successive_orders
