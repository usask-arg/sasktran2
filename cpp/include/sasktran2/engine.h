#pragma once

#include <sasktran2/config.h>
#include <sasktran2/atmosphere/atmosphere.h>
#include <sasktran2/viewinggeometry.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/viewinggeometry_internal.h>
#include <sasktran2/solartransmission.h>
#include <sasktran2/emission_source.h>
#include <sasktran2/source_integrator.h>
#include <sasktran2/output.h>

/** Internal essentially void class that the main Sasktran2 inherits from.
 * Necessary to remove the NSTOKES template to interface with SWIG.
 *
 */
class Sasktran2Interface {
  public:
    virtual ~Sasktran2Interface() {}
};

/** The main object for SASKTRAN2 that performs the radiative transfer
 * calculation.  The model is initialized with a Config object that contains
 * user settings, a Geometry object to specify grids and coordinate information,
 *  and a ViewingGeometryContainer object that contains information on the
 * viewing geometry.
 *
 *  After construction, the model may be called using an Atmosphere object.
 *
 * @tparam NSTOKES Number of stokes parameters, either 1 or 3
 */
template <int NSTOKES> class Sasktran2 : public Sasktran2Interface {
  private:
    sasktran2::Config m_config; /**< Engine-owned configuration snapshot */
    const sasktran2::viewinggeometry::ViewingGeometryContainer&
        m_viewing_geometry; /**< Internal reference to the viewing geometry */
    const sasktran2::Geometry*
        m_geometry; /**< Internal reference to the model geometry */
    const sasktran2::Geometry1D* m_geometry_1d = nullptr;
    const sasktran2::Geometry2D* m_geometry_2d = nullptr;

    /** Optional altitude-only refractive-index profile for every structured
     * 2D line of sight. Profiles are deliberately engine state rather than
     * Geometry2D state because an orbital calculation may select a different
     * atmospheric column for every ray and refresh it between retrieval
     * iterations. */
    std::vector<Eigen::VectorXd> m_refractive_profiles_2d;

    std::unique_ptr<const sasktran2::raytracing::RayTracerBase>
        m_raytracer; /**< Ray tracer that is internally constructed */
#ifdef SKTRAN_RUST_SUPPORT
    std::unique_ptr<const sasktran2::raytracing::RustRayTracer2D>
        m_raytracer_2d;
#endif

    sasktran2::viewinggeometry::InternalViewingGeometry
        m_internal_viewing_geometry; /**< Internal viewing geometry
                                        representation */

    std::unique_ptr<sasktran2::SourceIntegrator<NSTOKES>>
        m_source_integrator; /** integrator for the source terms */

    std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>>
        m_traced_ray_od_matrix; /**< Vector of matrices A such that A *
                                   atmosphere_extinction = OD for each layer in
                                   that ray */

    // We split source terms by characterization
    std::vector<std::unique_ptr<SourceTermInterface<NSTOKES>>>
        m_source_terms; /**< All of the source terms */
    std::vector<SourceTermInterface<NSTOKES>*>
        m_los_source_terms; /**< Source terms that are integrated over the line
                               of sight */
    std::vector<SourceTermInterface<NSTOKES>*>
        m_thermal_source; /**< Thermal source terms (planck, photochemical) */

    std::vector<sasktran2::WavelengthBlockDual<NSTOKES>> m_thread_radiance;

    // Persistent native-product storage permits an external scheduler to
    // initialize a product once and dispatch independent wavelength blocks.
    // The same storage is used by the ordinary OpenMP path.
    mutable std::vector<Eigen::VectorXd> m_native_jvp_tangents;
    mutable std::vector<Eigen::MatrixXd> m_native_vjp_gradients;
    using NativeProductStokesBlock =
        Eigen::Matrix<double, NSTOKES, Eigen::Dynamic, Eigen::RowMajor>;
    mutable std::vector<NativeProductStokesBlock> m_native_vjp_radiance;
    mutable std::vector<NativeProductStokesBlock> m_native_vjp_cotangent;

    // Thread storage to avoid reallocs on the derivatives of flux
    // Can't reuse radiance storage because flux NSTOKES is always 1 for flux
    std::vector<sasktran2::Dual<double, sasktran2::dualstorage::dense, 1>>
        m_thread_flux;

    /** Internal method to calculate all terms inside the engine that are only
     * geometry dependent
     */
    void calculate_geometry(bool refresh_los_sources = false);

    /** Internal method to construct the internal ray tracer
     */
    void construct_raytracer();

    /** Internal method to construct all of the source terms from the config
     */
    void construct_source_terms();

    /** Internal method to construct the source term integrator from the config
     *
     */
    void construct_integrator();

    /** Internal method to validate that the input atmosphere is in the correct
     * format and contains all of the necessary information.
     */
    void validate_input_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) const;

    void initialize();

  public:
    /** Constructs the model
     *
     * @param config user configuration options
     * @param geometry geometry and coordinate information
     * @param viewing_rays information on the viewing geometry
     */
    Sasktran2(const sasktran2::Config& config,
              const sasktran2::Geometry1D* geometry,
              const sasktran2::viewinggeometry::ViewingGeometryContainer&
                  viewing_rays)
        : m_config(config), m_viewing_geometry(viewing_rays),
          m_geometry(geometry), m_geometry_1d(geometry) {
        initialize();
    }

    /** Constructs the structured 2D model. */
    Sasktran2(const sasktran2::Config& config,
              const sasktran2::Geometry2D* geometry,
              const sasktran2::viewinggeometry::ViewingGeometryContainer&
                  viewing_rays)
        : m_config(config), m_viewing_geometry(viewing_rays),
          m_geometry(geometry), m_geometry_2d(geometry) {
        initialize();
    }

    virtual ~Sasktran2() {}

    /** Calculates the radiance for a given Atmosphere object.  May be called
     * any number of times with different Atmosphere objects.  Results are
     * stored in output.
     *
     * @param atmosphere
     * @param output
     */
    void calculate_radiance(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
        sasktran2::Output<NSTOKES>& output, bool only_initialize = false) const;

    void
    calculate_jvp(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
                  sasktran2::OutputJVP<NSTOKES>& output) const;

    void
    initialize_jvp(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
                   sasktran2::OutputJVP<NSTOKES>& output) const;

    void calculate_jvp_wavelength_thread(sasktran2::OutputJVP<NSTOKES>& output,
                                         int wavelength, int thread_idx) const;

    void
    calculate_vjp(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
                  sasktran2::OutputVJP<NSTOKES>& output) const;

    void
    initialize_vjp(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
                   sasktran2::OutputVJP<NSTOKES>& output) const;

    void calculate_vjp_block_thread(sasktran2::OutputVJP<NSTOKES>& output,
                                    const sasktran2::WavelengthBlock<>& block,
                                    int thread_idx) const;

    int effective_wavelength_batch_size(int num_wavelengths) const;

    using RowMajorMatrix =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    /** Copy current ground-intersection horizontal stencils. Rows are LOS
     * rays and columns are Geometry2D horizontal nodes. Non-ground rays have
     * all-zero rows. */
    void assign_2d_surface_interpolation_weights(
        Eigen::Ref<RowMajorMatrix> weights) const;

    using RowMajorIntMatrix =
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    /** Report use of Geometry2D's extended horizontal edge cells by traced
     * LOS paths. Columns are `(before_first_node, after_last_node)` and rows
     * are viewing rays. Only the portions of rays inside the atmosphere are
     * considered. */
    void
    assign_2d_horizontal_edge_usage(Eigen::Ref<RowMajorIntMatrix> usage) const;

    /** Replace the per-ray refractive profiles used by a structured 2D
     * engine and rebuild all geometry-dependent integration/source data.
     * `profiles` is `(num_los, num_altitudes)`.
     */
    void set_2d_refractive_profiles(
        const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic,
                                             Eigen::Dynamic, Eigen::RowMajor>>&
            profiles);

    /** Reports whether the complete line-of-sight model supports a native
     * derivative execution mode. */
    bool supports_linearization(sasktran2::LinearizationMode mode) const {
        return m_source_integrator->supports_linearization(mode,
                                                           m_los_source_terms);
    }

    /** Selects the best executable backend for a derivative product. */
    sasktran2::LinearizationBackend
    linearization_backend(sasktran2::LinearizationMode mode) const {
        if (supports_linearization(mode)) {
            return sasktran2::LinearizationBackend::Native;
        }
        if (mode != sasktran2::LinearizationMode::Jacobian &&
            supports_linearization(sasktran2::LinearizationMode::Jacobian)) {
            return sasktran2::LinearizationBackend::StreamingJacobian;
        }
        return sasktran2::LinearizationBackend::Unavailable;
    }

    void
    calculate_radiance_block_thread(sasktran2::Output<NSTOKES>& output,
                                    const sasktran2::WavelengthBlock<>& batch,
                                    int thread_idx) const;
};
