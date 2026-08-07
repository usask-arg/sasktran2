#include <sasktran2/successive_orders/source.h>

#include "source_internal.h"

namespace sasktran2::successive_orders {
    template <int NSTOKES>
    RustSource<NSTOKES>::RustSource(
        const sasktran2::raytracing::RayTracerBase& ray_tracer,
        const sasktran2::Geometry1D& geometry)
        : m_impl(std::make_unique<SourceImpl<NSTOKES>>(ray_tracer, geometry)) {}

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    RustSource<NSTOKES>::RustSource(
        const sasktran2::raytracing::RustRayTracer2D& ray_tracer,
        const sasktran2::Geometry2D& geometry)
        : m_impl(std::make_unique<SourceImpl<NSTOKES>>(ray_tracer, geometry)) {}
#endif

    template <int NSTOKES> RustSource<NSTOKES>::~RustSource() = default;

    template <int NSTOKES>
    bool RustSource<NSTOKES>::supports_geometry_dimension(int dimension) const {
        return m_impl->supports_geometry_dimension(dimension);
    }

    template <int NSTOKES>
    void
    RustSource<NSTOKES>::initialize_config(const sasktran2::Config& config) {
        m_impl->initialize_config(config);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        m_impl->initialize_geometry(internal_viewing);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_impl->initialize_atmosphere(atmosphere);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::initialize_atmosphere_native(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_impl->initialize_atmosphere_native(atmosphere);
    }

    template <int NSTOKES>
    void
    RustSource<NSTOKES>::set_wavelength_block_capacity(int block_capacity) {
        m_impl->set_wavelength_block_capacity(block_capacity);
    }

    template <int NSTOKES>
    int RustSource<NSTOKES>::maximum_wavelength_block_size() const {
        return m_impl->maximum_wavelength_block_size();
    }

    template <int NSTOKES>
    bool RustSource<NSTOKES>::supports_linearization(
        sasktran2::LinearizationMode mode) const {
        return m_impl->supports_linearization(mode);
    }

    template <int NSTOKES>
    bool RustSource<NSTOKES>::requires_integration() const {
        return m_impl->requires_integration();
    }

    template <int NSTOKES>
    bool RustSource<NSTOKES>::has_interior_source() const {
        return m_impl->has_interior_source();
    }

    template <int NSTOKES>
    void
    RustSource<NSTOKES>::calculate(const sasktran2::WavelengthBlock<>& block,
                                   int threadidx) {
        m_impl->calculate(block, threadidx);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::calculate_value(
        const sasktran2::WavelengthBlock<>& block, int threadidx) {
        m_impl->calculate_value(block, threadidx);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::begin_forward_state_capture(int num_wavelengths) {
        m_impl->begin_forward_state_capture(num_wavelengths);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::end_forward_state_capture() {
        m_impl->end_forward_state_capture();
    }

    template <int NSTOKES>
    bool RustSource<NSTOKES>::restore_forward_state(int wavelidx,
                                                    int threadidx) {
        return m_impl->restore_forward_state(wavelidx, threadidx);
    }

    template <int NSTOKES>
    bool RustSource<NSTOKES>::restore_forward_state_block(
        const sasktran2::WavelengthBlock<>& block, int threadidx) {
        return m_impl->restore_forward_state_block(block, threadidx);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::prepare_forward_state_for_jvp(
        int wavelidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        m_impl->prepare_forward_state_for_jvp(wavelidx, threadidx,
                                              native_tangent);
    }

    template <int NSTOKES>
    bool RustSource<NSTOKES>::restore_forward_state_for_jvp(
        int wavelidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        return m_impl->restore_forward_state_for_jvp(wavelidx, threadidx,
                                                     native_tangent);
    }

    template <int NSTOKES>
    std::size_t RustSource<NSTOKES>::retained_forward_state_bytes() const {
        return m_impl->retained_forward_state_bytes();
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::prepare_jvp(
        int wavelidx, int wavel_threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        m_impl->prepare_jvp(wavelidx, wavel_threadidx, native_tangent);
    }

    template <int NSTOKES>
    bool RustSource<NSTOKES>::prepare_jvp_block(
        const sasktran2::WavelengthBlock<>& block, int wavel_threadidx,
        const std::vector<Eigen::VectorXd>& native_tangents) {
        return m_impl->prepare_jvp_block(block, wavel_threadidx,
                                         native_tangents);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::select_jvp_block_lane(int lane,
                                                    int wavel_threadidx) {
        m_impl->select_jvp_block_lane(lane, wavel_threadidx);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::prepare_vjp(
        int wavelidx, int wavel_threadidx,
        Eigen::Ref<const Eigen::VectorXi> active_derivatives) {
        m_impl->prepare_vjp(wavelidx, wavel_threadidx, active_derivatives);
    }

    template <int NSTOKES>
    bool RustSource<NSTOKES>::begin_vjp_block(
        const sasktran2::WavelengthBlock<>& block, int wavel_threadidx,
        Eigen::Ref<const Eigen::VectorXi> active_derivatives) {
        return m_impl->begin_vjp_block(block, wavel_threadidx,
                                       active_derivatives);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::select_vjp_block_lane(int lane,
                                                    int wavel_threadidx) {
        m_impl->select_vjp_block_lane(lane, wavel_threadidx);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::stage_vjp_block_lane(int lane,
                                                   int wavel_threadidx) {
        m_impl->stage_vjp_block_lane(lane, wavel_threadidx);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::finalize_vjp_block(
        const sasktran2::WavelengthBlock<>& block, int wavel_threadidx,
        std::vector<Eigen::VectorXd>& native_gradients) const {
        m_impl->finalize_vjp_block(block, wavel_threadidx, native_gradients);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::finalize_vjp(
        int wavelidx, int wavel_threadidx,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        m_impl->finalize_vjp(wavelidx, wavel_threadidx, native_gradient);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::integrated_source(
        const sasktran2::WavelengthBlock<>& block, int losidx, int layeridx,
        int wavel_threadidx, int threadidx,
        const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        sasktran2::WavelengthBlockDual<NSTOKES>& source,
        typename SourceTermInterface<NSTOKES>::IntegrationDirection direction)
        const {
        m_impl->integrated_source(block, losidx, layeridx, wavel_threadidx,
                                  threadidx, layer, entrance_weights,
                                  exit_weights, shell_od, source, direction);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::end_of_ray_source(
        const sasktran2::WavelengthBlock<>& block, int losidx,
        int wavel_threadidx, int threadidx,
        sasktran2::WavelengthBlockDual<NSTOKES>& source) const {
        m_impl->end_of_ray_source(block, losidx, wavel_threadidx, threadidx,
                                  source);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::start_of_ray_source(
        const sasktran2::WavelengthBlock<>& block, int losidx,
        int wavel_threadidx, int threadidx,
        sasktran2::WavelengthBlockDual<NSTOKES>& source) const {
        m_impl->start_of_ray_source(block, losidx, wavel_threadidx, threadidx,
                                    source);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::integrated_source_jvp(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        sasktran2::RadianceJVP<NSTOKES>& source) const {
        m_impl->integrated_source_jvp(
            wavelidx, losidx, layeridx, wavel_threadidx, threadidx, layer,
            entrance_weights, exit_weights, shell_od, native_tangent, source);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::end_of_ray_source_jvp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        sasktran2::RadianceJVP<NSTOKES>& source) const {
        m_impl->end_of_ray_source_jvp(wavelidx, losidx, wavel_threadidx,
                                      threadidx, native_tangent, source);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::start_of_ray_source_jvp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent,
        sasktran2::RadianceJVP<NSTOKES>& source) const {
        m_impl->start_of_ray_source_jvp(wavelidx, losidx, wavel_threadidx,
                                        threadidx, native_tangent, source);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::integrated_source_vjp(
        int wavelidx, int losidx, int layeridx, int wavel_threadidx,
        int threadidx, const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        const Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::Vector<double, NSTOKES>> source_value,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        m_impl->integrated_source_vjp(wavelidx, losidx, layeridx,
                                      wavel_threadidx, threadidx, layer,
                                      entrance_weights, exit_weights, shell_od,
                                      cotangent, source_value, native_gradient);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::end_of_ray_source_vjp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        const Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        m_impl->end_of_ray_source_vjp(wavelidx, losidx, wavel_threadidx,
                                      threadidx, cotangent, native_gradient);
    }

    template <int NSTOKES>
    void RustSource<NSTOKES>::start_of_ray_source_vjp(
        int wavelidx, int losidx, int wavel_threadidx, int threadidx,
        const Eigen::Vector<double, NSTOKES>& value_before,
        Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::VectorXd> native_gradient) const {
        m_impl->start_of_ray_source_vjp(wavelidx, losidx, wavel_threadidx,
                                        threadidx, value_before, cotangent,
                                        native_gradient);
    }

    template class RustSource<1>;
    template class RustSource<3>;
} // namespace sasktran2::successive_orders
