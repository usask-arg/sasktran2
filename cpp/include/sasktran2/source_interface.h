#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/viewinggeometry_internal.h>
#include <sasktran2/wavelength_block.h>
#include <sasktran2/dual.h>
#include <sasktran2/config.h>
#include <sasktran2/atmosphere/atmosphere.h>

namespace sasktran2 {
    template <int NSTOKES> class SourceIntegrator;

    template <int NSTOKES> struct RadianceJVP {
        Eigen::Vector<double, NSTOKES> value =
            Eigen::Vector<double, NSTOKES>::Zero();
        Eigen::Vector<double, NSTOKES> jvp =
            Eigen::Vector<double, NSTOKES>::Zero();

        void set_zero() {
            value.setZero();
            jvp.setZero();
        }
    };

    /** Native derivative execution modes that an engine component may
     * implement. Capability reporting is kept separate from execution so
     * specialized JVP/VJP interfaces can be added without placeholder hooks. */
    enum class LinearizationMode { Jacobian = 0, JVP = 1, VJP = 2 };
    enum class LinearizationBackend {
        Unavailable = 0,
        StreamingJacobian = 1,
        Native = 2
    };
} // namespace sasktran2

/** Base interface class that provides source term functionality to the Engine
 *
 */
template <int NSTOKES> class SourceTermInterface {
  public:
    // Enum determinig the direction the source integration happens
    // Knowing the direction of integration when calculating the source can
    // enable some optimization. This is essentially a contract from the user to
    // the source term indicating the order that integrated_source will be
    // called in backward indicates we are starting at the end of the LOS, and
    // moving towards the observer forward indicates we are starting at the
    // observer and moving towards the end of the LOS none indicates that the
    // source should not assume any direction of integration
    enum class IntegrationDirection { forward, backward, none };

  private:
    template <int> friend class sasktran2::SourceIntegrator;

    template <int N>
    using FixedIntegratedSourceFunction = void (*)(
        const SourceTermInterface&, const sasktran2::WavelengthBlock<N>&,
        int losidx, int layeridx, int wavel_threadidx, int threadidx,
        const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        sasktran2::WavelengthBlockDual<NSTOKES>& source,
        IntegrationDirection direction);

    template <int N>
    static void dynamic_integrated_source_fallback(
        const SourceTermInterface& source_term,
        const sasktran2::WavelengthBlock<N>& block, int losidx, int layeridx,
        int wavel_threadidx, int threadidx,
        const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        sasktran2::WavelengthBlockDual<NSTOKES>& source,
        IntegrationDirection direction) {
        source_term.integrated_source(
            block.dynamic(), losidx, layeridx, wavel_threadidx, threadidx,
            layer, entrance_weights, exit_weights, shell_od, source, direction);
    }

    FixedIntegratedSourceFunction<1> m_integrated_source_1 =
        &dynamic_integrated_source_fallback<1>;
    FixedIntegratedSourceFunction<4> m_integrated_source_4 =
        &dynamic_integrated_source_fallback<4>;

  protected:
    template <int N>
    void set_fixed_integrated_source_dispatch(
        FixedIntegratedSourceFunction<N> function) {
        if constexpr (N == 1) {
            m_integrated_source_1 = function;
        } else if constexpr (N == 4) {
            m_integrated_source_4 = function;
        } else {
            static_assert(N == 1 || N == 4,
                          "Unsupported fixed wavelength block size");
        }
    }

  public:
    virtual ~SourceTermInterface(){};

    virtual void initialize_config(const sasktran2::Config& config){};

    /** Initializes any geometry information that is required for calculating
     * the source term.  This method is called after the line of sight rays ar
     * traced.
     *
     * @param internal_viewing Information on the internal viewing geometry,
     * los_rays and flux observers
     */
    virtual void initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing){};

    /**
     *
     */
    virtual void initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere){};

    /** Sets the wavelength block capacity negotiated by the engine for the
     * current calculation. Sources may use this to size per-block storage. */
    virtual void set_wavelength_block_capacity(int block_capacity) {}

    /** Triggers calculation for a contiguous block of wavelengths. */
    virtual void calculate(const sasktran2::WavelengthBlock<>& block,
                           int threadidx){};

    /** Maximum number of wavelengths this source can process together. */
    virtual int maximum_wavelength_block_size() const { return 1; }

    // TODO: Is Dual proper here? what about when the source term derivative is
    // sparse? Maybe it isn't that important... Should we be templated over
    // NSTOKES in the interface?

    /** Calculates the integrated source term for a given layer.
     *
     * @param block Contiguous wavelengths being integrated
     * @param losidx Raw index pointing to the ray that was previously passed
     * in initialize_geometry
     * @param layeridx Raw index pointing to the layer that was previosuly
     * passed in initialize_geometry
     * @param layer The layer that we are integrating over
     * @param source The returned source term
     */
    virtual void integrated_source(
        const sasktran2::WavelengthBlock<>&, int losidx, int layeridx,
        int wavel_threadidx, int threadidx,
        const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        sasktran2::WavelengthBlockDual<NSTOKES>& source,
        IntegrationDirection direction = IntegrationDirection::none) const = 0;

  private:
    /** Internal fixed-width hook. Sources that do not register a specialized
     * kernel fall back to the dynamic virtual block interface. */
    template <int N>
    void dispatch_integrated_source(
        const sasktran2::WavelengthBlock<N>& block, int losidx, int layeridx,
        int wavel_threadidx, int threadidx,
        const sasktran2::raytracing::TracedLayer& layer,
        const sasktran2::raytracing::GridWeightStencilView& entrance_weights,
        const sasktran2::raytracing::GridWeightStencilView& exit_weights,
        const sasktran2::WavelengthBlockODView& shell_od,
        sasktran2::WavelengthBlockDual<NSTOKES>& source,
        IntegrationDirection direction = IntegrationDirection::none) const {
        if constexpr (N == 1) {
            m_integrated_source_1(*this, block, losidx, layeridx,
                                  wavel_threadidx, threadidx, layer,
                                  entrance_weights, exit_weights, shell_od,
                                  source, direction);
        } else if constexpr (N == 4) {
            m_integrated_source_4(*this, block, losidx, layeridx,
                                  wavel_threadidx, threadidx, layer,
                                  entrance_weights, exit_weights, shell_od,
                                  source, direction);
        } else {
            integrated_source(block.dynamic(), losidx, layeridx,
                              wavel_threadidx, threadidx, layer,
                              entrance_weights, exit_weights, shell_od, source,
                              direction);
        }
    }

  public:
    /** Calculates the source term at the end of the ray.  Common examples of
     * this are ground scattering, ground emission, or the solar radiance if
     * looking directly at the sun.
     *
     * @param block Contiguous wavelengths being integrated
     * @param losidx Raw index pointing to the ray that was previously passed in
     * initialize_geometry
     * @param source The returned source term
     */
    virtual void end_of_ray_source(
        const sasktran2::WavelengthBlock<>&, int losidx, int wavel_threadidx,
        int threadidx,
        sasktran2::WavelengthBlockDual<NSTOKES>& source) const = 0;

    /** Calculates the radiance at the start of the ray, i.e., the source term
     * has done the equivalent of the integration along the ray.  This is useful
     * if the source term has it's own way of performing integration that is
     * different than the standard method used in the model.  It can also be
     * used if the source term calculates quantities that can be used to get the
     * radiance directly instead of doing integration.
     *
     *  Typically source terms will only either implement start_of_ray_source,
     * or integrated_source + end_of_ray_source and not both
     *
     * @param block Contiguous wavelengths being integrated
     * @param losidx
     * @param wavel_threadidx
     * @param threadidx
     * @param source
     */
    virtual void start_of_ray_source(
        const sasktran2::WavelengthBlock<>&, int losidx, int wavel_threadidx,
        int threadidx,
        sasktran2::WavelengthBlockDual<NSTOKES>& source) const = 0;

    virtual void
    flux(int wavelidx, int fluxidx, int wavelt_threadidx, int threadidx,
         sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>& flux,
         sasktran2::Config::FluxType flux_type) const {
        spdlog::error("Flux calculation not implemented for this source term");
    }

    virtual bool requires_integration() const { return true; }

    /** Returns true when the source contributes within atmospheric layers. */
    virtual bool has_interior_source() const { return true; }

    /** Reports native support for a derivative execution mode. */
    virtual bool
    supports_linearization(sasktran2::LinearizationMode mode) const {
        return mode == sasktran2::LinearizationMode::Jacobian;
    }

    virtual void end_of_ray_source_jvp(int, int, int, int,
                                       Eigen::Ref<const Eigen::VectorXd>,
                                       sasktran2::RadianceJVP<NSTOKES>&) const {
        throw std::logic_error("Native source JVP is not implemented");
    }

    virtual void
    integrated_source_jvp(int, int, int, int, int,
                          const sasktran2::raytracing::TracedLayer&,
                          const sasktran2::raytracing::GridWeightStencilView&,
                          const sasktran2::raytracing::GridWeightStencilView&,
                          const sasktran2::WavelengthBlockODView&,
                          Eigen::Ref<const Eigen::VectorXd>,
                          sasktran2::RadianceJVP<NSTOKES>&) const {
        throw std::logic_error(
            "Native integrated source JVP is not implemented");
    }

    /** Applies a source contribution evaluated after line-of-sight
     * integration to a value/tangent pair.  Unlike end-of-ray and interior
     * hooks, this operation may transform the radiance already accumulated by
     * earlier sources, so the source owns both the input and output state. */
    virtual void
    start_of_ray_source_jvp(int, int, int, int,
                            Eigen::Ref<const Eigen::VectorXd>,
                            sasktran2::RadianceJVP<NSTOKES>&) const {
        throw std::logic_error(
            "Native start-of-ray source JVP is not implemented");
    }

    virtual void end_of_ray_source_vjp(int, int, int, int,
                                       const Eigen::Vector<double, NSTOKES>&,
                                       Eigen::Ref<Eigen::VectorXd>) const {
        throw std::logic_error("Native source VJP is not implemented");
    }

    virtual void
    integrated_source_vjp(int, int, int, int, int,
                          const sasktran2::raytracing::TracedLayer&,
                          const sasktran2::raytracing::GridWeightStencilView&,
                          const sasktran2::raytracing::GridWeightStencilView&,
                          const sasktran2::WavelengthBlockODView&,
                          const Eigen::Vector<double, NSTOKES>&,
                          Eigen::Ref<Eigen::VectorXd>) const {
        throw std::logic_error(
            "Native integrated source VJP is not implemented");
    }

    /** Pulls a cotangent through a post-integration source contribution.
     * `value_before` is the radiance immediately before this source was
     * applied.  Implementations add parameter gradients to `native_gradient`
     * and replace `cotangent` with the cotangent of `value_before`. */
    virtual void start_of_ray_source_vjp(int, int, int, int,
                                         const Eigen::Vector<double, NSTOKES>&,
                                         Eigen::Vector<double, NSTOKES>&,
                                         Eigen::Ref<Eigen::VectorXd>) const {
        throw std::logic_error(
            "Native start-of-ray source VJP is not implemented");
    }

    /** Returns whether this source supports the requested atmosphere
     * dimensionality. Sources without an interior contribution are independent
     * of the layer grid and need not override this method. */
    virtual bool supports_geometry_dimension(int dimension) const {
        return dimension == 1 || !has_interior_source();
    }

    /** Returns true when this source can report every derivative column that
     * it may modify. The source integrator uses this information to avoid
     * scaling derivative columns that are still identically zero. */
    virtual bool supports_sparse_derivative_tracking() const { return false; }

    /** Appends derivative columns that may be modified by the source at the
     * end of a ray. */
    virtual void append_end_of_ray_active_derivatives(int,
                                                      std::vector<int>&) const {
    }

    /** Appends derivative columns that may be modified by the source within a
     * layer. */
    virtual void append_interior_active_derivatives(int, int,
                                                    std::vector<int>&) const {}
};
