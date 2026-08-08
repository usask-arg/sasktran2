#pragma once

#include <sasktran2/raytracing.h>
#include <sasktran2/source_interface.h>
#include <memory>

#ifdef SKTRAN_RUST_SUPPORT
namespace sasktran2::successive_orders {
    template <int NSTOKES> class SourceImpl;

    /** Engine-facing façade for the Rust successive-orders source.
     *
     * C++ traces LOS/source rays and forwards engine lifecycle callbacks. All
     * solver state, geometry packing, wavelength batching, and linearization
     * details are hidden in the module-private implementation.
     */
    template <int NSTOKES>
    class RustSource final
        : public SourceTermInterface<NSTOKES>,
          public sasktran2::WavelengthBlockLinearizationInterface {
      public:
        RustSource(const sasktran2::raytracing::RayTracerBase& ray_tracer,
                   const sasktran2::Geometry1D& geometry);
        RustSource(const sasktran2::raytracing::RustRayTracer2D& ray_tracer,
                   const sasktran2::Geometry2D& geometry);
        ~RustSource() override;

        RustSource(const RustSource&) = delete;
        RustSource& operator=(const RustSource&) = delete;
        RustSource(RustSource&&) = delete;
        RustSource& operator=(RustSource&&) = delete;

        sasktran2::WavelengthBlockLinearizationInterface*
        wavelength_block_linearization() override {
            return this;
        }

        bool supports_geometry_dimension(int dimension) const override;
        void initialize_config(const sasktran2::Config& config) override;
        void initialize_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) override;
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override;
        void initialize_atmosphere_native(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override;

        void set_wavelength_block_capacity(int block_capacity) override;
        int maximum_wavelength_block_size() const override;
        bool supports_linearization(
            sasktran2::LinearizationMode mode) const override;
        bool requires_integration() const override;
        bool has_interior_source() const override;

        void calculate(const sasktran2::WavelengthBlock<>& block,
                       int threadidx) override;
        void calculate_value(const sasktran2::WavelengthBlock<>& block,
                             int threadidx) override;

        void begin_forward_state_capture(int num_wavelengths) override;
        void end_forward_state_capture() override;
        bool restore_forward_state(int wavelidx, int threadidx) override;
        bool
        restore_forward_state_block(const sasktran2::WavelengthBlock<>& block,
                                    int threadidx) override;
        void prepare_forward_state_for_jvp(
            int wavelidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent) override;
        bool restore_forward_state_for_jvp(
            int wavelidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent) override;
        std::size_t retained_forward_state_bytes() const override;

        void
        prepare_jvp(int wavelidx, int wavel_threadidx,
                    Eigen::Ref<const Eigen::VectorXd> native_tangent) override;
        bool prepare_jvp_block(
            const sasktran2::WavelengthBlock<>& block, int wavel_threadidx,
            const std::vector<Eigen::VectorXd>& native_tangents) override;
        void select_jvp_block_lane(int lane, int wavel_threadidx) override;

        void prepare_vjp(
            int wavelidx, int wavel_threadidx,
            Eigen::Ref<const Eigen::VectorXi> active_derivatives) override;
        bool begin_vjp_block(
            const sasktran2::WavelengthBlock<>& block, int wavel_threadidx,
            Eigen::Ref<const Eigen::VectorXi> active_derivatives) override;
        void select_vjp_block_lane(int lane, int wavel_threadidx) override;
        void stage_vjp_block_lane(int lane, int wavel_threadidx) override;
        void finalize_vjp_block(
            const sasktran2::WavelengthBlock<>& block, int wavel_threadidx,
            std::vector<Eigen::VectorXd>& native_gradients) const override;
        void finalize_vjp(
            int wavelidx, int wavel_threadidx,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const override;

        void integrated_source(
            const sasktran2::WavelengthBlock<>& block, int losidx, int layeridx,
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
            const override;
        void end_of_ray_source(
            const sasktran2::WavelengthBlock<>& block, int losidx,
            int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockDual<NSTOKES>& source) const override;
        void start_of_ray_source(
            const sasktran2::WavelengthBlock<>& block, int losidx,
            int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockDual<NSTOKES>& source) const override;

        void integrated_source_jvp(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            Eigen::Ref<const Eigen::VectorXd> native_tangent,
            sasktran2::RadianceJVP<NSTOKES>& source) const override;
        void end_of_ray_source_jvp(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent,
            sasktran2::RadianceJVP<NSTOKES>& source) const override;
        void start_of_ray_source_jvp(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent,
            sasktran2::RadianceJVP<NSTOKES>& source) const override;

        void integrated_source_vjp(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            const Eigen::Vector<double, NSTOKES>& value_before,
            Eigen::Vector<double, NSTOKES>& cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const override;
        void end_of_ray_source_vjp(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            const Eigen::Vector<double, NSTOKES>& value_before,
            Eigen::Vector<double, NSTOKES>& cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const override;
        void start_of_ray_source_vjp(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            const Eigen::Vector<double, NSTOKES>& value_before,
            Eigen::Vector<double, NSTOKES>& cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const override;

      private:
        std::unique_ptr<SourceImpl<NSTOKES>> m_impl;
    };
} // namespace sasktran2::successive_orders
#endif
