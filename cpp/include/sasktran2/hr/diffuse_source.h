#pragma once

#include "sasktran2/viewinggeometry_internal.h"
#include <sasktran2/internal_common.h>
#include <sasktran2/source_interface.h>
#include <sasktran2/source_integrator.h>
#include <sasktran2/do_source.h>
#include <sasktran2/hr/diffuse_point.h>

#ifdef SKTRAN_RUST_SUPPORT
#include "sasktran2-core/src/successive_orders/cxx.rs.h"
#include "sasktran2-core/src/twostream/cxx.rs.h"
#endif

namespace sasktran2::hr {
    /** Thread specific storage for the diffuse table
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES> struct DiffuseTableThreadStorage {
        sasktran2::Dual<double, sasktran2::dualstorage::dense>
            m_incoming_radiances; /**< Total incoming stokes vectors, [stokes,
                                     direction, diffusegrid] */
        sasktran2::Dual<double, sasktran2::dualstorage::dense>
            m_firstorder_radiances; /**< Total first order incoming radiances
                                       (usaully single scatter + optional
                                       thermal/photochemical), [stokes,
                                       direction, diffusegrid] */
        sasktran2::Dual<double, sasktran2::dualstorage::dense>
            m_outgoing_sources; /**< Outgoing source terms, [stokes, direction,
                                   diffusegrid] */

        std::vector<Eigen::MatrixXd>
            point_scattering_matrices; /** For each point, outgoing source =
                                          scattering matrix @ incoming */

        /** CSR transport values. Each traced ray owns a disjoint row range,
         *  so source threads can fill this single vector concurrently. */
        Eigen::VectorXd accumulation_summed_values;

        std::vector<double>
            rust_boundary_scattering_values; /**< Row-major dense boundary
                                                scattering blocks. */

        Eigen::VectorXd rust_transport_value_tangent;
        std::vector<double> rust_boundary_scattering_value_tangent;
        Eigen::VectorXd rust_first_order_forcing_tangent;
        std::vector<double> rust_solution_jvp; /**< Legacy callback fallback;
                                                   empty for the Rust source. */

        /** Per-wavelength scalar attenuation scratch produced together with
         *  the transport values in Rust and consumed by the exact source. */
        std::vector<double> rust_layer_optical_depth;
        std::vector<double> rust_layer_attenuation;
        std::vector<double> rust_layer_prefix_attenuation;
        std::vector<double> rust_ray_end_attenuation;
        mutable std::vector<double> rust_end_of_ray_source;
        std::vector<double> rust_end_of_ray_source_tangent;
        mutable std::vector<double> rust_end_of_ray_source_gradient;
        std::vector<double> rust_single_scatter_legendre_tangent;
        mutable std::vector<double> rust_single_scatter_extinction_gradient;
        mutable std::vector<double> rust_single_scatter_albedo_gradient;
        mutable std::vector<double> rust_single_scatter_legendre_gradient;
        mutable std::vector<double> rust_single_scatter_solar_gradient;

        /** Wavelength-contiguous outgoing sources retained until LOS
         *  integration completes. */
        std::vector<double> rust_batch_outgoing_sources;

        /** One-profile scratch for a two-stream initial state. */
        std::vector<double> twostream_extinction;
        std::vector<double> twostream_ssa;
        std::vector<double> twostream_legendre;
        std::vector<double> twostream_delta_m;
        std::vector<double> twostream_emission;
        std::vector<double> twostream_surface_albedo;
    };

    /** Minimal state required to differentiate a converged fixed point. The
     * wavelength-dependent sparse transport values are rebuilt rather than
     * retained because they are much larger than these two vectors. */
    struct DiffuseTableForwardState {
        std::vector<double> first_order_forcing;
        std::vector<double> solution;
        bool valid = false;
    };

    struct DiffuseTableDerivativeActivity {
        bool prepared = false;
        bool extinction = false;
        bool ssa = false;
        bool surface = false;
        bool emission = false;
        bool surface_emission = false;
        std::vector<int> scattering_groups;

        bool scattering() const { return !scattering_groups.empty(); }
    };

    /** An implementation of the successive orders of scattering technique.  We
     * call this HR mostly for historic reasons.
     *
     *
     * Note: Currently this source is not perfectly linearized.
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES>
    class DiffuseTable : public SourceTermInterface<NSTOKES> {
        using SInterpolator =
            std::vector<RaySourceInterpolationWeights<NSTOKES>>;

      private:
        std::vector<DiffuseTableThreadStorage<NSTOKES>>
            m_thread_storage; /** Thread (optical) data [nthreads] */

        const sasktran2::Config*
            m_config; /**< Internal reference to the config */

        sasktran2::SourceIntegrator<NSTOKES>
            m_integrator; /**< Integrator used to get OD along the rays, as well
                             as integrate the initial sources to get the
                             incoming radiances */
        std::vector<SourceTermInterface<NSTOKES>*>
            m_initial_sources; /**< Initial sources used to get the first order
                                  incoming radiances */

        std::vector<std::unique_ptr<SourceTermInterface<NSTOKES>>>
            m_initial_owned_sources; /**< Initial sources used to get the first
                                        order incoming radiances */

        const sasktran2::raytracing::RayTracerBase&
            m_raytracer; /**< Raytracer used for the diffuse point incoming rays
                          */
        const sasktran2::Geometry& m_geometry; /** Global geometry object */
        const sasktran2::grids::AltitudeGrid&
            m_altitude_grid; /** Vertical grid shared by supported geometries */
        const sasktran2::Geometry1D* m_geometry_1d;
        const sasktran2::Geometry2D* m_geometry_2d;
#ifdef SKTRAN_RUST_SUPPORT
        const sasktran2::raytracing::RustRayTracer2D* m_raytracer_2d;
        std::vector<std::vector<
            ::rust::Box<sasktran2::rust::twostream::RustTwoStreamSource>>>
            m_twostream_initializers; /**< [wavelength worker][source profile]
                                       */
        std::vector<std::vector<std::vector<std::pair<int, double>>>>
            m_twostream_atmosphere_weights; /**< [profile][level][weight] */
        int m_twostream_num_source_altitudes = 0;
#endif

        std::unique_ptr<sasktran2::grids::SourceLocationInterpolator>
            m_location_interpolator; /** Interpolates location */

        sasktran2::viewinggeometry::InternalViewingGeometry
            m_internal_viewing; /** Traced incoming rays to all diffuse
                                points */

        std::vector<std::unique_ptr<DiffusePoint<NSTOKES>>>
            m_diffuse_points; /** Stacked vector of all interior diffuse points,
                                 including ground, interpolated using
                                 m_location_interpolator */
        std::vector<bool>
            m_diffuse_point_full_calculation; /** True if we are doing the full
                                                 incoming calculation at this
                                                 diffuse point, false if it is
                                                 interpolated from spherical
                                                 corrections */

        std::vector<
            std::unique_ptr<sasktran2::hr::IncomingOutgoingSpherePair<NSTOKES>>>
            m_unit_sphere_pairs; /** Unit sphere ownership for the diffuse
                                    points, never accessed after construction by
                                    this class */

        const sasktran2::atmosphere::Atmosphere<NSTOKES>*
            m_atmosphere; /** Reference to the atmosphere object */

        std::vector<int>
            m_diffuse_incoming_index_map; /** element i is the start index of
                                             point i in m_incoming_traced_rays
                                           */
        std::vector<int>
            m_diffuse_outgoing_index_map; /** element i is the start index of
                                             point i in outgoing directions */

        std::vector<std::vector<std::pair<int, double>>>
            m_diffuse_point_interpolation_weights; /** Interpolation mapping
                                                      from the global
                                                      coordinates to the diffuse
                                                      locations.  Mostly used to
                                                      generate the scattering
                                                      matrices. */

        SInterpolator
            m_los_source_weights; /** Interpolator mapping from line of sight
                                     points to source terms in this table */
        SInterpolator
            m_diffuse_source_weights; /** Interpolator mapping from incoming
                                         rays to source terms in this table */

        Eigen::SparseMatrix<double, Eigen::RowMajor>
            m_do_to_diffuse_outgoing_interpolator; /** Mapping from the DO
                                                      source terms to the
                                                      outgoing sphere sources */
        DOSourceInterpolatedPostProcessing<NSTOKES, -1>*
            m_do_source; /** Reference to the DO source */

        // Accumulation matrix sparsity
        Eigen::VectorXi m_inner_indicies;
        Eigen::VectorXi m_outer_starts;
        bool m_use_rust_solver;
        std::size_t m_num_outgoing_values = 0;
        int m_wavelength_block_capacity = 1;
        mutable std::vector<std::vector<Eigen::VectorXd>>
            m_los_solution_cotangents;
        mutable std::vector<int> m_last_nonconverged_vjp_warning_wavelength;
        mutable std::vector<DiffuseTableDerivativeActivity>
            m_active_vjp_derivatives;
        std::vector<DiffuseTableForwardState> m_forward_states;
        const sasktran2::atmosphere::Atmosphere<NSTOKES>*
            m_forward_state_atmosphere = nullptr;
        std::uint64_t m_forward_state_atmosphere_revision = 0;
        bool m_capture_forward_state = false;
        std::vector<bool> m_deferred_jvp_transport_restore;

#ifdef SKTRAN_RUST_SUPPORT
        mutable std::vector<::rust::Box<
            sasktran2::rust::successive_orders::RustSuccessiveOrdersSolver>>
            m_rust_solvers;
        std::vector<::rust::Box<
            sasktran2::rust::successive_orders::RustTransportSparsity>>
            m_rust_transport_sparsity;
        std::vector<
            ::rust::Box<sasktran2::rust::successive_orders::RustRayTransport>>
            m_rust_ray_transports;
        std::vector<::rust::Box<sasktran2::rust::successive_orders::
                                    RustScatteringCoefficientInterpolator>>
            m_rust_scattering_interpolators;
        std::vector<::rust::Box<
            sasktran2::rust::successive_orders::RustSourceInterpolator>>
            m_rust_source_interpolators;
        std::vector<std::uint32_t> m_rust_transport_ray_layer_offsets;
        std::vector<std::size_t> m_rust_boundary_scattering_offsets;
#endif

      private:
        sasktran2::grids::Grid generate_cos_sza_grid(double min_cos_sza,
                                                     double max_cos_sza);
        sasktran2::grids::AltitudeGrid generate_altitude_grid();

        void construct_diffuse_points();
        std::vector<std::vector<int>> trace_incoming_rays();
        void generate_scattering_matrices(int wavelidx, int threadidx);
        void generate_accumulation_matrix(int wavelidx, int threadidx);
        void prepare_wavelength(int wavelidx, int threadidx,
                                bool value_only = false);
        void iterate_to_solution(int wavelidx, int threadidx);
#ifdef SKTRAN_RUST_SUPPORT
        void initialize_twostream_initial_guess_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing);
        void generate_twostream_initial_guess(int wavelidx, int threadidx);
        void pack_rust_scattering_coefficient_jvp(
            int wavelidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent);
        void pack_rust_boundary_scattering_jvp(
            int wavelidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent);
        void accumulate_rust_scattering_coefficient_vjp(
            int wavelidx, int threadidx,
            const std::vector<int>& active_scattering_groups,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const;
        void accumulate_rust_boundary_scattering_vjp(
            int wavelidx, ::rust::Slice<const double> boundary_gradient,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const;
        void iterate_to_solution_rust(int wavelidx, int threadidx);
        void assemble_rust_ray_transport(int wavelidx, int threadidx,
                                         bool assemble_first_order = false);
        void pack_rust_end_of_ray_source(int wavelidx, int threadidx) const;
        void pack_rust_end_of_ray_source_jvp(
            int wavelidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent);
        void accumulate_rust_end_of_ray_source_vjp(
            int wavelidx, int threadidx,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const;
        void assemble_rust_ray_transport_jvp(
            int wavelidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent);
        void accumulate_rust_transport_vjp_with_first_order(
            int wavelidx, int threadidx,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const;
        void initialize_rust_source_interpolator(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing);
        void release_cpp_transport_geometry();
        void capture_forward_state(int wavelidx, int threadidx);
        void restore_transport_operator(int wavelidx, int threadidx);
#endif
        void interpolate_sources(const Eigen::VectorXd& old_outgoing,
                                 sasktran2::Dual<double>& new_outgoing);
        void generate_source_interpolation_weights(
            const std::vector<sasktran2::raytracing::TracedRay>& rays,
            SInterpolator& interpolator) const;
        void generate_source_interpolation_weights(
            const sasktran2::raytracing::TracedRay& ray,
            RaySourceInterpolationWeights<NSTOKES>& interpolator,
            std::vector<std::pair<int, double>>& temp_location_storage,
            std::vector<std::pair<int, double>>& temp_direction_storage,
            std::vector<std::pair<int, double>>& temp_atmosphere_storage) const;
        void compile_accumulation_row(
            RaySourceInterpolationWeights<NSTOKES>& weights,
            std::vector<int>& transport_columns,
            std::vector<std::pair<int, std::uint16_t*>>& sorting_helper) const;

        void construct_accumulation_sparsity(
            const std::vector<std::vector<int>>& transport_columns);

        Eigen::Vector3d
        rotate_unit_vector(const Eigen::Vector3d& vector,
                           const Eigen::Vector3d& initial_position,
                           const Eigen::Vector3d& new_position) const;

      public:
        DiffuseTable(const sasktran2::raytracing::RayTracerBase& ray_tracer,
                     const sasktran2::Geometry1D& geometry,
                     bool use_rust_solver = false);

#ifdef SKTRAN_RUST_SUPPORT
        DiffuseTable(const sasktran2::raytracing::RustRayTracer2D& ray_tracer,
                     const sasktran2::Geometry2D& geometry,
                     bool use_rust_solver = false);
#endif

        bool supports_geometry_dimension(int dimension) const override {
            return dimension == 1 ||
                   (dimension == 2 && m_geometry_2d != nullptr);
        }

        /** Initializes the config inside the source term
         *
         * @param config
         */
        virtual void initialize_config(const sasktran2::Config& config);

        /** Initializes any geometry information that is required for
         * calculating the source term.  This method is called after the line of
         * sight rays ar traced.
         *
         * @param internal_viewing Information on the internal viewing geometry,
         * los_rays and flux observers
         */
        virtual void initialize_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) override;

        /**
         *
         */
        virtual void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere);

        void initialize_atmosphere_native(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override;

        void set_wavelength_block_capacity(int block_capacity) override;

        int maximum_wavelength_block_size() const override {
#ifdef SKTRAN_RUST_SUPPORT
            // Four wavelengths match the engine's optimized fixed-width
            // integration dispatch while bounding retained source data.
            return m_use_rust_solver ? 4 : 1;
#endif
            return 1;
        }

        bool native_products_available() const {
#ifdef SKTRAN_RUST_SUPPORT
            return m_use_rust_solver && m_config != nullptr;
#else
            return false;
#endif
        }

        /** Triggers an internal calculation of the source term.  This method is
         * called at the beginning of each 'wavelength' calculation.
         *
         * @param wavelidx Index of the wavelength being calculated
         */
        void calculate(const sasktran2::WavelengthBlock<>& block,
                       int threadidx) override;

        void calculate_value(const sasktran2::WavelengthBlock<>& block,
                             int threadidx) override;

        void begin_forward_state_capture(int num_wavelengths) override;

        void end_forward_state_capture() override;

        bool restore_forward_state(int wavelidx, int threadidx) override;

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

        void prepare_vjp(
            int wavelidx, int wavel_threadidx,
            Eigen::Ref<const Eigen::VectorXi> active_derivatives) override;

        void finalize_vjp(
            int wavelidx, int wavel_threadidx,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const override;

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

        void integrated_source_vjp(
            int wavelidx, int losidx, int layeridx, int wavel_threadidx,
            int threadidx, const sasktran2::raytracing::TracedLayer& layer,
            const sasktran2::raytracing::GridWeightStencilView&
                entrance_weights,
            const sasktran2::raytracing::GridWeightStencilView& exit_weights,
            const sasktran2::WavelengthBlockODView& shell_od,
            const Eigen::Vector<double, NSTOKES>& cotangent,
            Eigen::Ref<Eigen::Vector<double, NSTOKES>> source_value,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const override;

        void end_of_ray_source_vjp(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            const Eigen::Vector<double, NSTOKES>& cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const override;

        /**
         * @brief Not used for the HR source.
         *
         * @param wavelidx
         * @param losidx
         * @param wavel_threadidx
         * @param threadidx
         * @param source
         */
        void start_of_ray_source(
            const sasktran2::WavelengthBlock<>& block, int losidx,
            int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockDual<NSTOKES>& source) const override;

        void start_of_ray_source_jvp(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent,
            sasktran2::RadianceJVP<NSTOKES>& source) const override;

        void start_of_ray_source_vjp(
            int wavelidx, int losidx, int wavel_threadidx, int threadidx,
            const Eigen::Vector<double, NSTOKES>& value_before,
            Eigen::Vector<double, NSTOKES>& cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const override;
    };
} // namespace sasktran2::hr
