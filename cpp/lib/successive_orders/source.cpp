#include "factory.h"

#include "first_order.h"
#include "fixed_point.h"
#include "geometry.h"
#include "problem.h"
#include "ray_transport.h"
#include "scattering_assembler.h"

#include <sasktran2/geometry.h>

#include <spdlog/spdlog.h>

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace sasktran2::successive_orders {
    /** Private concrete source. Engine code interacts with it only through the
     * SourceTermInterface returned by make_successive_orders_source(). */
    template <int NSTOKES>
    class SuccessiveOrdersSource final : public SourceTermInterface<NSTOKES> {
      public:
        SuccessiveOrdersSource(
            const sasktran2::raytracing::RayTracerBase& raytracer,
            const sasktran2::Geometry1D& geometry);
#ifdef SKTRAN_RUST_SUPPORT
        SuccessiveOrdersSource(
            const sasktran2::raytracing::RustRayTracer2D& raytracer,
            const sasktran2::Geometry2D& geometry,
            std::shared_ptr<
                sasktran2::solartransmission::SolarTransmissionTable2D>
                shared_solar_table = nullptr);
#endif
        ~SuccessiveOrdersSource() override;

        SuccessiveOrdersSource(const SuccessiveOrdersSource&) = delete;
        SuccessiveOrdersSource&
        operator=(const SuccessiveOrdersSource&) = delete;

        void initialize_config(const sasktran2::Config& config) override;
        void initialize_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) override;
        void refresh_los_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) override;
        void initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)
            override;

        void calculate(const sasktran2::WavelengthBlock<>& block,
                       int threadidx) override;
        void calculate_jvp(
            const sasktran2::WavelengthBlock<>& block, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent) override;
        void calculate_vjp(const sasktran2::WavelengthBlock<>& block,
                           int threadidx) override;
        void finalize_vjp(
            const sasktran2::WavelengthBlock<>& block, int threadidx,
            Eigen::Ref<Eigen::MatrixXd> native_gradient) const override;

        int maximum_wavelength_block_size() const override { return 1; }
        bool requires_integration() const override { return false; }
        bool has_interior_source() const override { return false; }
        bool supports_linearization(
            sasktran2::LinearizationMode mode) const override;
        bool supports_geometry_dimension(int dimension) const override {
            if (dimension == 1) {
                return true;
            }
#ifdef SKTRAN_RUST_SUPPORT
            return dimension == 2;
#else
            return false;
#endif
        }
        bool supports_sparse_derivative_tracking() const override {
            return true;
        }

        void integrated_source(
            const sasktran2::WavelengthBlock<>&, int, int, int, int,
            const sasktran2::raytracing::TracedLayer&,
            const sasktran2::raytracing::GridWeightStencilView&,
            const sasktran2::raytracing::GridWeightStencilView&,
            const sasktran2::WavelengthBlockODView&,
            sasktran2::WavelengthBlockDual<NSTOKES>&,
            typename SourceTermInterface<NSTOKES>::IntegrationDirection)
            const override {}

        void end_of_ray_source(
            const sasktran2::WavelengthBlock<>&, int, int, int,
            sasktran2::WavelengthBlockDual<NSTOKES>&) const override {}
        void start_of_ray_source(
            const sasktran2::WavelengthBlock<>& block, int losidx,
            int wavel_threadidx, int threadidx,
            sasktran2::WavelengthBlockDual<NSTOKES>& source) const override;

        void
        end_of_ray_source_jvp(int, int, int, int,
                              Eigen::Ref<const Eigen::VectorXd>,
                              sasktran2::RadianceJVP<NSTOKES>&) const override {
        }
        void start_of_ray_source_jvp(
            int wavelength, int losidx, int wavel_threadidx, int threadidx,
            Eigen::Ref<const Eigen::VectorXd> native_tangent,
            sasktran2::RadianceJVP<NSTOKES>& source) const override;

        void end_of_ray_source_vjp(int, int, int, int,
                                   const Eigen::Vector<double, NSTOKES>&,
                                   Eigen::Vector<double, NSTOKES>&,
                                   Eigen::Ref<Eigen::VectorXd>) const override {
        }
        void start_of_ray_source_vjp(
            int wavelength, int losidx, int wavel_threadidx, int threadidx,
            const Eigen::Vector<double, NSTOKES>& value_before,
            Eigen::Vector<double, NSTOKES>& cotangent,
            Eigen::Ref<Eigen::VectorXd> native_gradient) const override;

      private:
        class Impl;
        std::unique_ptr<Impl> m_impl;
    };

    namespace {
        // Solver tolerances may deliberately be much tighter than radiometric
        // accuracy requirements. Only warn when the unresolved fixed-point
        // update is large enough to plausibly change the diffuse radiance at
        // the approximately 0.1% level. The absolute floor handles vanishing
        // states without tying warning behavior back to the solver settings.
        void warn_if_not_converged(const FixedPointDiagnostics& diagnostics,
                                   const FixedPointSettings& settings,
                                   int wavelength, const char* calculation) {
            if (settings.convergence_enabled() &&
                diagnostics.warrants_convergence_warning()) {
                spdlog::warn(
                    "C++ successive-orders {} did not converge at wavelength "
                    "index {} after {} iterations (residual {}, requested "
                    "threshold {}, estimated relative unresolved state {}, "
                    "warning threshold {}). Returning the current result; "
                    "its accuracy may be reduced.",
                    calculation, wavelength, diagnostics.iterations,
                    diagnostics.residual_norm,
                    diagnostics.convergence_threshold,
                    diagnostics.relative_residual(),
                    FixedPointDiagnostics::material_warning_relative_tolerance);
            }
        }

        void warn_if_implicit_derivative_is_unchecked(
            const FixedPointSettings& settings, int wavelength,
            const char* calculation) {
            if (!settings.convergence_enabled()) {
                spdlog::warn(
                    "C++ successive-orders {} was requested at wavelength "
                    "index {} with convergence tolerances disabled. The "
                    "implicit derivative assumes a converged fixed point and "
                    "may differ from the derivative of the configured "
                    "fixed-count primal iteration.",
                    calculation, wavelength);
            }
        }

        void validate_single_wavelength_block(
            const sasktran2::WavelengthBlock<>& block) {
            if (block.count != 1) {
                throw std::invalid_argument(
                    "The C++ successive-orders source currently processes "
                    "one wavelength per block");
            }
        }
    } // namespace

    template <int NSTOKES> struct ScatteringAssemblerAdapter;

    template <> struct ScatteringAssemblerAdapter<1> {
        using Assembler = ScalarScatteringAssembler;

        static void
        assemble_jvp(const Assembler& assembler,
                     const sasktran2::atmosphere::Atmosphere<1>& atmosphere,
                     int wavelength,
                     Eigen::Ref<const Eigen::VectorXd> native_tangent,
                     ProblemParameterData<1>& tangent) {
            assembler.assemble_jvp(atmosphere, wavelength, native_tangent,
                                   tangent.atmospheric_coefficients,
                                   tangent.ground_values);
        }

        static void
        accumulate_vjp(const Assembler& assembler,
                       const sasktran2::atmosphere::Atmosphere<1>& atmosphere,
                       int wavelength, const ProblemParameterData<1>& gradient,
                       Eigen::Ref<Eigen::VectorXd> native_gradient) {
            assembler.accumulate_vjp(atmosphere, wavelength,
                                     gradient.atmospheric_coefficients,
                                     gradient.ground_values, native_gradient);
        }
    };

    template <> struct ScatteringAssemblerAdapter<3> {
        using Assembler = VectorScatteringAssembler;

        static void
        assemble_jvp(const Assembler& assembler,
                     const sasktran2::atmosphere::Atmosphere<3>& atmosphere,
                     int wavelength,
                     Eigen::Ref<const Eigen::VectorXd> native_tangent,
                     ProblemParameterData<3>& tangent) {
            assembler.assemble_jvp(atmosphere, wavelength, native_tangent,
                                   tangent.atmospheric_coefficients,
                                   tangent.ground_values);
        }

        static void
        accumulate_vjp(const Assembler& assembler,
                       const sasktran2::atmosphere::Atmosphere<3>& atmosphere,
                       int wavelength, const ProblemParameterData<3>& gradient,
                       Eigen::Ref<Eigen::VectorXd> native_gradient) {
            assembler.accumulate_vjp(atmosphere, wavelength,
                                     gradient.atmospheric_coefficients,
                                     gradient.ground_values, native_gradient);
        }
    };

    /** Shared scalar/vector driver.  Only scattering-parameter assembly differs
     * between the scalar and polarized coefficient-space operators; transport,
     * fixed-point solves, and native products use identical state and thread
     * ownership. */
    template <int NSTOKES> class SuccessiveOrdersImplementation {
      private:
        using Adapter = ScatteringAssemblerAdapter<NSTOKES>;
        using Assembler = typename Adapter::Assembler;
        using Atmosphere = sasktran2::atmosphere::Atmosphere<NSTOKES>;

        struct WavelengthState {
            WavelengthState(const RayTransportMap& transport_map,
                            const RayTransportMap& los_map,
                            const Assembler& assembler)
                : transport(transport_map.sparsity()),
                  los_transport(los_map.sparsity()),
                  scattering(assembler.create_operator()),
                  problem(transport, scattering) {
                workspace.resize(transport, scattering);
                forcing.resize(scattering.input_size());
                state.resize(scattering.output_size());
                state_tangent.resize(scattering.output_size());
                direct_transport_tangent.resize(scattering.input_size());
                const int los_size = los_transport.sparsity().rows() * NSTOKES;
                los_value.resize(los_size);
                los_tangent.resize(los_size);
                los_cotangent.resize(los_size);
                state_cotangent.resize(scattering.output_size());
                los_value_gradient.resize(los_transport.values().size());
                los_transport_tangent.resize(los_transport.values().size());
            }

            TransportOperator transport;
            TransportOperator los_transport;
            ScatteringOperator<NSTOKES> scattering;
            Problem<NSTOKES> problem;
            ProblemWorkspace<NSTOKES> workspace;
            ProblemParameterData<NSTOKES> tangent;
            ProblemParameterData<NSTOKES> gradient;
            RayTransportWorkspace transport_vjp_workspace;
            RayTransportWorkspace los_vjp_workspace;
            Eigen::VectorXd forcing;
            Eigen::VectorXd state;
            Eigen::VectorXd state_tangent;
            Eigen::VectorXd layer_state_projection;
            Eigen::VectorXd ground_state_projection;
            Eigen::VectorXd direct_transport_tangent;
            Eigen::VectorXd los_value;
            Eigen::VectorXd los_tangent;
            Eigen::VectorXd los_cotangent;
            Eigen::VectorXd state_cotangent;
            Eigen::VectorXd los_value_gradient;
            Eigen::VectorXd los_transport_tangent;
            Eigen::VectorXd adjoint;
            Eigen::MatrixXd jacobian;
            int active_wavelength = -1;
            int state_wavelength = -1;
            bool transport_state_projected = false;
        };

      public:
        SuccessiveOrdersImplementation(
            const sasktran2::raytracing::RayTracerBase& raytracer,
            const sasktran2::Geometry1D& geometry)
            : m_source_geometry(raytracer, geometry),
              m_first_order(geometry, raytracer) {}

#ifdef SKTRAN_RUST_SUPPORT
        SuccessiveOrdersImplementation(
            const sasktran2::raytracing::RustRayTracer2D& raytracer,
            const sasktran2::Geometry2D& geometry,
            std::shared_ptr<
                sasktran2::solartransmission::SolarTransmissionTable2D>
                shared_solar_table = nullptr)
            : m_source_geometry(raytracer, geometry),
              m_first_order(geometry, raytracer,
                            std::move(shared_solar_table)) {}
#endif

        void initialize_config(const sasktran2::Config& config) {
            invalidate_geometry();
            m_config = &config;
            m_geometry_settings.num_incoming = config.num_hr_incoming();
            m_geometry_settings.num_outgoing = config.num_hr_outgoing();
            m_geometry_settings.num_sza = config.num_do_sza();
            m_geometry_settings.num_threads = config.num_threads();
            m_geometry_settings.include_refraction =
                config.multiple_scatter_refraction();
            m_geometry_settings.altitude_grid_m =
                config.successive_orders_altitude_grid_m();
            m_geometry_settings.horizontal_angle_grid_radians =
                config.successive_orders_horizontal_angle_grid_radians();

            m_solver_settings.maximum_iterations =
                config.num_hr_spherical_iterations();
            m_solver_settings.relative_tolerance =
                config.successive_orders_relative_tolerance();
            m_solver_settings.absolute_tolerance =
                config.successive_orders_absolute_tolerance();
            m_solver_settings.anderson_depth =
                config.successive_orders_anderson_depth();
            m_solver_settings.damping = config.successive_orders_damping();
            m_solver_settings.validate();
            m_geometry_settings.validate();
            m_first_order.initialize_config(config);
        }

        void initialize_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) {
            if (m_config == nullptr) {
                throw std::logic_error(
                    "C++ successive-orders config must be initialized before "
                    "geometry");
            }
            invalidate_geometry();
            m_source_geometry.initialize(internal_viewing, m_geometry_settings);
            m_transport_map = std::make_unique<RayTransportMap>(
                m_source_geometry.incoming_interpolation(),
                m_source_geometry.total_num_outgoing(),
                m_source_geometry.transport_row_offsets(),
                m_source_geometry.transport_column_indices());
            m_los_map = std::make_unique<RayTransportMap>(
                m_source_geometry.los_interpolation(),
                m_source_geometry.total_num_outgoing(),
                m_source_geometry.los_transport_row_offsets(),
                m_source_geometry.los_transport_column_indices());
            m_scattering_assembler = std::make_unique<Assembler>(
                m_source_geometry, m_config->num_do_streams());

            const int output_size = m_los_map->num_rays() * NSTOKES;
            m_thread_los_cotangent.resize(m_config->num_threads());
            for (auto& cotangent : m_thread_los_cotangent) {
                cotangent.setZero(output_size);
            }
            m_first_order.initialize_geometry(m_source_geometry);
            if (m_first_order.can_release_incoming_geometry()) {
                m_source_geometry.release_incoming_traced_rays();
            }
            m_geometry_initialized = true;
        }

        void refresh_los_geometry(
            const sasktran2::viewinggeometry::InternalViewingGeometry&
                internal_viewing) {
            if (!m_geometry_initialized) {
                initialize_geometry(internal_viewing);
                return;
            }

            m_atmosphere = nullptr;
            m_has_atmosphere_revision = false;
            m_wavelength_state.clear();
            m_vector_state_cache.clear();
            m_source_geometry.refresh_los(internal_viewing);
            m_los_map = std::make_unique<RayTransportMap>(
                m_source_geometry.los_interpolation(),
                m_source_geometry.total_num_outgoing(),
                m_source_geometry.los_transport_row_offsets(),
                m_source_geometry.los_transport_column_indices());

            const int output_size = m_los_map->num_rays() * NSTOKES;
            m_thread_los_cotangent.resize(m_config->num_threads());
            for (auto& cotangent : m_thread_los_cotangent) {
                cotangent.setZero(output_size);
            }
        }

        void initialize_atmosphere(const Atmosphere& atmosphere) {
            if (!m_geometry_initialized) {
                throw std::logic_error(
                    "C++ successive-orders geometry must be initialized "
                    "before the atmosphere");
            }
            // Preserve directly-mutable C++ Atmosphere behavior at revision
            // zero. Python constituent builds and explicitly tracked native
            // storage opt into cache reuse by calling mark_changed().
            const bool tracked_atmosphere =
                atmosphere.revision() != 0 && m_atmosphere == &atmosphere &&
                m_has_atmosphere_revision &&
                m_atmosphere_instance_id == atmosphere.instance_id();
            if (tracked_atmosphere &&
                m_atmosphere_revision == atmosphere.revision()) {
                return;
            }
            const bool volume_changed =
                !tracked_atmosphere ||
                m_atmosphere_volume_revision != atmosphere.volume_revision();
            m_atmosphere = nullptr;
            // Scalar workspaces are compact enough to retain one complete
            // state per wavelength. Vector workspaces stay thread-local, but
            // their converged solution is cached separately below so native
            // products do not have to repeat the primal solve.
            if constexpr (NSTOKES == 1) {
                if (static_cast<int>(m_wavelength_state.size()) !=
                    atmosphere.num_wavel()) {
                    create_wavelength_states(atmosphere.num_wavel());
                }
            } else {
                if (static_cast<int>(m_wavelength_state.size()) !=
                    m_config->num_wavelength_threads()) {
                    create_wavelength_states(
                        m_config->num_wavelength_threads());
                }
                if (static_cast<int>(m_vector_state_cache.size()) !=
                    atmosphere.num_wavel()) {
                    m_vector_state_cache.clear();
                    m_vector_state_cache.resize(atmosphere.num_wavel());
                }
            }
            for (auto& state : m_wavelength_state) {
                state->active_wavelength = -1;
            }
            m_first_order.initialize_atmosphere(atmosphere, volume_changed);
            m_atmosphere = &atmosphere;
            m_atmosphere_instance_id = atmosphere.instance_id();
            m_atmosphere_revision = atmosphere.revision();
            m_atmosphere_volume_revision = atmosphere.volume_revision();
            m_has_atmosphere_revision = true;
        }

        void calculate(const sasktran2::WavelengthBlock<>& block,
                       int threadidx) {
            validate_single_wavelength_block(block);
            auto& work = prepare_primal(block.start, threadidx);
            if (m_atmosphere->num_deriv() > 0) {
                warn_if_implicit_derivative_is_unchecked(
                    m_solver_settings, block.start, "full Jacobian");
                calculate_jacobian(block.start, threadidx, work);
            }
        }

        void calculate_jvp(const sasktran2::WavelengthBlock<>& block,
                           int threadidx,
                           Eigen::Ref<const Eigen::VectorXd> native_tangent) {
            validate_single_wavelength_block(block);
            validate_calculation(block.start, threadidx);
            if (native_tangent.size() != m_atmosphere->num_deriv()) {
                throw std::invalid_argument(
                    "Invalid C++ successive-orders native tangent size");
            }
            auto& work = wavelength_state(block.start, threadidx);
            const bool compact_scalar =
                m_first_order.uses_compact_scalar_kernel();
            work.tangent.resize(work.transport, work.scattering,
                                !compact_scalar);
            work.tangent.set_zero();
            prepare_primal(block.start, threadidx);
            if (native_tangent.isZero(0.0)) {
                work.state_tangent.setZero();
                work.los_tangent.setZero();
            } else {
                warn_if_implicit_derivative_is_unchecked(m_solver_settings,
                                                         block.start, "JVP");
                if (compact_scalar) {
                    if (!work.transport_state_projected) {
                        m_first_order.project_transport_state(
                            work.state, work.layer_state_projection,
                            work.ground_state_projection);
                        work.transport_state_projected = true;
                    }
                    m_first_order.calculate_jvp_with_transport(
                        block.start, threadidx, native_tangent,
                        work.layer_state_projection,
                        work.ground_state_projection,
                        work.direct_transport_tangent, work.forcing,
                        work.tangent.forcing);
                } else {
                    m_first_order.calculate_jvp(block.start, threadidx,
                                                native_tangent, work.forcing,
                                                work.tangent.forcing);
                }
                if (!compact_scalar) {
                    m_transport_map->assemble_jvp(
                        *m_atmosphere, block.start, native_tangent,
                        work.tangent.transport_values);
                }
                Adapter::assemble_jvp(*m_scattering_assembler, *m_atmosphere,
                                      block.start, native_tangent,
                                      work.tangent);
                const auto jvp_diagnostics = work.problem.solve_jvp(
                    work.forcing, work.state, work.tangent, work.state_tangent,
                    m_solver_settings, work.workspace,
                    compact_scalar ? &work.direct_transport_tangent : nullptr);
                warn_if_not_converged(jvp_diagnostics, m_solver_settings,
                                      block.start, "JVP solve");
                m_los_map->assemble_jvp(*m_atmosphere, block.start,
                                        native_tangent,
                                        work.los_transport_tangent);
                work.los_transport.template apply_jvp_stokes<NSTOKES>(
                    work.state, work.state_tangent, work.los_transport_tangent,
                    work.los_tangent);
            }
            work.jacobian.resize(0, 0);
            work.active_wavelength = block.start;
        }

        void calculate_vjp(const sasktran2::WavelengthBlock<>& block,
                           int threadidx) {
            validate_single_wavelength_block(block);
            prepare_primal(block.start, threadidx);
            wavelength_state(block.start, threadidx).jacobian.resize(0, 0);
            reset_vjp_cotangents(threadidx);
        }

        void finalize_vjp(const sasktran2::WavelengthBlock<>& block,
                          int threadidx,
                          Eigen::Ref<Eigen::MatrixXd> native_gradient) {
            validate_single_wavelength_block(block);
            validate_active(block.start, threadidx);
            if (native_gradient.rows() != m_atmosphere->num_deriv() ||
                native_gradient.cols() != 1) {
                throw std::invalid_argument(
                    "Invalid C++ successive-orders native VJP storage");
            }
            auto& work = wavelength_state(block.start, threadidx);
            work.los_cotangent.setZero();
            const auto [first_thread, last_thread] =
                source_thread_range(threadidx);
            for (int thread = first_thread; thread < last_thread; ++thread) {
                work.los_cotangent += m_thread_los_cotangent[thread];
            }
            auto gradient = native_gradient.col(0);
            warn_if_implicit_derivative_is_unchecked(m_solver_settings,
                                                     block.start, "VJP");
            reverse(block.start, threadidx, work, work.los_cotangent, gradient);
        }

        void start_of_ray_source(
            const sasktran2::WavelengthBlock<>& block, int losidx,
            int wavel_threadidx,
            sasktran2::WavelengthBlockDual<NSTOKES>& source) const {
            validate_single_wavelength_block(block);
            validate_active(block.start, wavel_threadidx);
            validate_los(losidx);
            const auto& work = wavelength_state(block.start, wavel_threadidx);
            const int output_offset = losidx * NSTOKES;
            source.value.col(0) +=
                work.los_value.template segment<NSTOKES>(output_offset);
            if (work.jacobian.rows() != 0) {
                if (source.derivative_size() != work.jacobian.cols()) {
                    throw std::invalid_argument(
                        "C++ successive-orders Jacobian storage mismatch");
                }
                for (int derivative = 0; derivative < source.derivative_size();
                     ++derivative) {
                    source.derivative(derivative, 1).col(0) +=
                        work.jacobian.col(derivative)
                            .template segment<NSTOKES>(output_offset);
                }
            }
        }

        void
        start_of_ray_source_jvp(int wavelength, int losidx, int wavel_threadidx,
                                sasktran2::RadianceJVP<NSTOKES>& source) const {
            validate_active(wavelength, wavel_threadidx);
            validate_los(losidx);
            const auto& work = wavelength_state(wavelength, wavel_threadidx);
            const int output_offset = losidx * NSTOKES;
            source.value +=
                work.los_value.template segment<NSTOKES>(output_offset);
            source.jvp +=
                work.los_tangent.template segment<NSTOKES>(output_offset);
        }

        void start_of_ray_source_vjp(
            int wavelength, int losidx, int wavel_threadidx, int threadidx,
            Eigen::Vector<double, NSTOKES>& cotangent) const {
            validate_active(wavelength, wavel_threadidx);
            validate_los(losidx);
            if (threadidx < 0 ||
                threadidx >= static_cast<int>(m_thread_los_cotangent.size())) {
                throw std::out_of_range(
                    "C++ successive-orders VJP thread is out of range");
            }
            m_thread_los_cotangent[threadidx].template segment<NSTOKES>(
                losidx * NSTOKES) += cotangent;
        }

        bool supports_linearization(sasktran2::LinearizationMode) const {
            return true;
        }

      private:
        void create_wavelength_states(int count) {
            if (count < 1 || m_transport_map == nullptr ||
                m_los_map == nullptr || m_scattering_assembler == nullptr) {
                throw std::logic_error(
                    "Cannot allocate C++ successive-orders wavelength state");
            }
            m_wavelength_state.clear();
            m_vector_state_cache.clear();
            m_wavelength_state.reserve(count);
            for (int index = 0; index < count; ++index) {
                m_wavelength_state.emplace_back(
                    std::make_unique<WavelengthState>(
                        *m_transport_map, *m_los_map, *m_scattering_assembler));
            }
        }

        int wavelength_state_index(int wavelength, int threadidx) const {
            if constexpr (NSTOKES == 1) {
                return wavelength;
            } else {
                return threadidx;
            }
        }

        WavelengthState& wavelength_state(int wavelength, int threadidx) {
            const int index = wavelength_state_index(wavelength, threadidx);
            if (index < 0 ||
                index >= static_cast<int>(m_wavelength_state.size())) {
                throw std::logic_error(
                    "C++ successive-orders wavelength cache is invalid");
            }
            return *m_wavelength_state[index];
        }

        const WavelengthState& wavelength_state(int wavelength,
                                                int threadidx) const {
            const int index = wavelength_state_index(wavelength, threadidx);
            if (index < 0 ||
                index >= static_cast<int>(m_wavelength_state.size())) {
                throw std::logic_error(
                    "C++ successive-orders wavelength cache is invalid");
            }
            return *m_wavelength_state[index];
        }

        void invalidate_geometry() {
            m_geometry_initialized = false;
            m_atmosphere = nullptr;
            m_has_atmosphere_revision = false;
            m_wavelength_state.clear();
            m_thread_los_cotangent.clear();
            m_scattering_assembler.reset();
            m_los_map.reset();
            m_transport_map.reset();
        }

        std::pair<int, int> source_thread_range(int wavelength_thread) const {
            const int first_thread = wavelength_thread;
            const int last_thread =
                first_thread + m_config->num_source_threads();
            if (first_thread < 0 ||
                last_thread > static_cast<int>(m_thread_los_cotangent.size())) {
                throw std::logic_error(
                    "C++ successive-orders thread layout is invalid");
            }
            return {first_thread, last_thread};
        }

        void validate_calculation(int wavelength, int threadidx) const {
            if (m_atmosphere == nullptr || threadidx < 0 ||
                threadidx >= m_config->num_wavelength_threads() ||
                wavelength < 0 || wavelength >= m_atmosphere->num_wavel()) {
                throw std::invalid_argument(
                    "Invalid C++ successive-orders calculation");
            }
        }

        void validate_active(int wavelength, int threadidx) const {
            validate_calculation(wavelength, threadidx);
            if (wavelength_state(wavelength, threadidx).active_wavelength !=
                wavelength) {
                throw std::logic_error(
                    "C++ successive-orders wavelength state is not active");
            }
        }

        void validate_los(int losidx) const {
            if (losidx < 0 || losidx >= m_los_map->num_rays()) {
                throw std::out_of_range(
                    "C++ successive-orders LOS index is out of range");
            }
        }

        void assemble_values(int wavelength, WavelengthState& work) {
            if (!m_first_order.uses_compact_scalar_kernel()) {
                m_transport_map->assemble_values(*m_atmosphere, wavelength,
                                                 work.transport);
            }
            m_los_map->assemble_values(*m_atmosphere, wavelength,
                                       work.los_transport);
            m_scattering_assembler->assemble_values(*m_atmosphere, wavelength,
                                                    work.scattering);
        }

        WavelengthState& prepare_primal(int wavelength, int threadidx) {
            validate_calculation(wavelength, threadidx);
            auto& work = wavelength_state(wavelength, threadidx);
            if (work.active_wavelength == wavelength) {
                return work;
            }
            work.active_wavelength = -1;
            work.transport_state_projected = false;
            assemble_values(wavelength, work);
            if (m_first_order.uses_compact_scalar_kernel()) {
                m_first_order.calculate_with_transport(
                    wavelength, threadidx, work.transport, work.forcing);
            } else {
                m_first_order.calculate(wavelength, threadidx, work.forcing);
            }
            // Retrieval loops usually perturb the atmosphere only slightly,
            // so a previous converged state is a useful initial guess. Fixed
            // iteration mode must instead start from zero: otherwise its
            // result would depend on earlier engine calculations. Preserve
            // the historical zero-state semantics when no iterations are
            // requested as well.
            if constexpr (NSTOKES == 1) {
                if (work.state_wavelength != wavelength ||
                    m_solver_settings.maximum_iterations == 0 ||
                    !m_solver_settings.convergence_enabled()) {
                    work.state.setZero();
                }
            } else {
                const auto& cached_state = m_vector_state_cache[wavelength];
                if (m_solver_settings.maximum_iterations > 0 &&
                    m_solver_settings.convergence_enabled() &&
                    cached_state.size() == work.problem.state_size()) {
                    work.state = cached_state;
                } else {
                    work.state.setZero();
                }
            }
            const auto diagnostics = work.problem.solve(
                work.forcing, work.state, m_solver_settings, work.workspace);
            if constexpr (NSTOKES == 3) {
                if (m_solver_settings.maximum_iterations > 0 &&
                    m_solver_settings.convergence_enabled()) {
                    m_vector_state_cache[wavelength] = work.state;
                }
            }
            work.state_wavelength =
                m_solver_settings.maximum_iterations == 0 ? -1 : wavelength;
            warn_if_not_converged(diagnostics, m_solver_settings, wavelength,
                                  "primal solve");
            work.los_transport.template apply_stokes<NSTOKES>(work.state,
                                                              work.los_value);
            work.los_tangent.setZero();
            work.jacobian.resize(0, 0);
            work.active_wavelength = wavelength;
            return work;
        }

        void reverse(int wavelength, int threadidx, WavelengthState& work,
                     Eigen::Ref<const Eigen::VectorXd> los_cotangent,
                     Eigen::Ref<Eigen::VectorXd> native_gradient,
                     bool emit_convergence_warning = true) {
            work.los_transport.template apply_vjp_stokes<NSTOKES>(
                work.state, los_cotangent, work.state_cotangent,
                work.los_value_gradient);
            m_los_map->accumulate_vjp(*m_atmosphere, wavelength,
                                      work.los_value_gradient, native_gradient,
                                      work.los_vjp_workspace);
            const auto diagnostics = work.problem.solve_vjp(
                work.forcing, work.state, work.state_cotangent, work.gradient,
                m_solver_settings, work.workspace, work.adjoint,
                !m_first_order.uses_compact_scalar_kernel());
            if (emit_convergence_warning) {
                warn_if_not_converged(diagnostics, m_solver_settings,
                                      wavelength, "VJP solve");
            }
            if (!m_first_order.uses_compact_scalar_kernel()) {
                m_transport_map->accumulate_vjp(
                    *m_atmosphere, wavelength, work.gradient.transport_values,
                    native_gradient, work.transport_vjp_workspace);
            }
            Adapter::accumulate_vjp(*m_scattering_assembler, *m_atmosphere,
                                    wavelength, work.gradient, native_gradient);
            if (m_first_order.uses_compact_scalar_kernel()) {
                if (work.transport_state_projected) {
                    m_first_order.accumulate_vjp_with_projected_transport(
                        wavelength, threadidx, work.layer_state_projection,
                        work.ground_state_projection, work.gradient.forcing,
                        native_gradient);
                } else {
                    m_first_order.accumulate_vjp_with_transport(
                        wavelength, threadidx, work.state,
                        work.gradient.forcing, native_gradient);
                }
            } else {
                m_first_order.accumulate_vjp(wavelength, threadidx,
                                             work.gradient.forcing,
                                             native_gradient);
            }
        }

        void calculate_jacobian(int wavelength, int threadidx,
                                WavelengthState& work) {
            const int outputs = m_los_map->num_rays() * NSTOKES;
            const int derivatives = m_atmosphere->num_deriv();
            work.jacobian.resize(outputs, derivatives);
            work.los_cotangent.setZero();
            Eigen::VectorXd native_gradient(derivatives);
            if (m_first_order.uses_compact_scalar_kernel() &&
                !work.transport_state_projected) {
                m_first_order.project_transport_state(
                    work.state, work.layer_state_projection,
                    work.ground_state_projection);
                work.transport_state_projected = true;
            }
            for (int output = 0; output < outputs; ++output) {
                work.los_cotangent(output) = 1.0;
                native_gradient.setZero();
                reverse(wavelength, threadidx, work, work.los_cotangent,
                        native_gradient, output == 0);
                work.jacobian.row(output) = native_gradient.transpose();
                work.los_cotangent(output) = 0.0;
            }
        }

        void reset_vjp_cotangents(int wavelength_thread) {
            const auto [first_thread, last_thread] =
                source_thread_range(wavelength_thread);
            for (int thread = first_thread; thread < last_thread; ++thread) {
                m_thread_los_cotangent[thread].setZero();
            }
        }

        const sasktran2::Config* m_config = nullptr;
        const Atmosphere* m_atmosphere = nullptr;
        std::uint64_t m_atmosphere_instance_id = 0;
        std::uint64_t m_atmosphere_revision = 0;
        std::uint64_t m_atmosphere_volume_revision = 0;
        bool m_has_atmosphere_revision = false;
        SourceGeometrySettings m_geometry_settings;
        FixedPointSettings m_solver_settings;
        SourceGeometry1D m_source_geometry;
        FirstOrderProvider<NSTOKES> m_first_order;
        std::unique_ptr<RayTransportMap> m_transport_map;
        std::unique_ptr<RayTransportMap> m_los_map;
        std::unique_ptr<Assembler> m_scattering_assembler;
        std::vector<std::unique_ptr<WavelengthState>> m_wavelength_state;
        std::vector<Eigen::VectorXd> m_vector_state_cache;
        mutable std::vector<Eigen::VectorXd> m_thread_los_cotangent;
        bool m_geometry_initialized = false;
    };

    template <>
    class SuccessiveOrdersSource<1>::Impl
        : public SuccessiveOrdersImplementation<1> {
      public:
        using SuccessiveOrdersImplementation<1>::SuccessiveOrdersImplementation;
    };

    template <>
    class SuccessiveOrdersSource<3>::Impl
        : public SuccessiveOrdersImplementation<3> {
      public:
        using SuccessiveOrdersImplementation<3>::SuccessiveOrdersImplementation;
    };

    template <int NSTOKES>
    SuccessiveOrdersSource<NSTOKES>::SuccessiveOrdersSource(
        const sasktran2::raytracing::RayTracerBase& raytracer,
        const sasktran2::Geometry1D& geometry)
        : m_impl(std::make_unique<Impl>(raytracer, geometry)) {}

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    SuccessiveOrdersSource<NSTOKES>::SuccessiveOrdersSource(
        const sasktran2::raytracing::RustRayTracer2D& raytracer,
        const sasktran2::Geometry2D& geometry,
        std::shared_ptr<sasktran2::solartransmission::SolarTransmissionTable2D>
            shared_solar_table)
        : m_impl(std::make_unique<Impl>(raytracer, geometry,
                                        std::move(shared_solar_table))) {}
#endif

    template <int NSTOKES>
    SuccessiveOrdersSource<NSTOKES>::~SuccessiveOrdersSource() = default;

    template <int NSTOKES>
    void SuccessiveOrdersSource<NSTOKES>::initialize_config(
        const sasktran2::Config& config) {
        m_impl->initialize_config(config);
    }

    template <int NSTOKES>
    void SuccessiveOrdersSource<NSTOKES>::initialize_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        m_impl->initialize_geometry(internal_viewing);
    }

    template <int NSTOKES>
    void SuccessiveOrdersSource<NSTOKES>::refresh_los_geometry(
        const sasktran2::viewinggeometry::InternalViewingGeometry&
            internal_viewing) {
        m_impl->refresh_los_geometry(internal_viewing);
    }

    template <int NSTOKES>
    void SuccessiveOrdersSource<NSTOKES>::initialize_atmosphere(
        const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) {
        m_impl->initialize_atmosphere(atmosphere);
    }

    template <int NSTOKES>
    void SuccessiveOrdersSource<NSTOKES>::calculate(
        const sasktran2::WavelengthBlock<>& block, int threadidx) {
        m_impl->calculate(block, threadidx);
    }

    template <int NSTOKES>
    void SuccessiveOrdersSource<NSTOKES>::calculate_jvp(
        const sasktran2::WavelengthBlock<>& block, int threadidx,
        Eigen::Ref<const Eigen::VectorXd> native_tangent) {
        m_impl->calculate_jvp(block, threadidx, native_tangent);
    }

    template <int NSTOKES>
    void SuccessiveOrdersSource<NSTOKES>::calculate_vjp(
        const sasktran2::WavelengthBlock<>& block, int threadidx) {
        m_impl->calculate_vjp(block, threadidx);
    }

    template <int NSTOKES>
    void SuccessiveOrdersSource<NSTOKES>::finalize_vjp(
        const sasktran2::WavelengthBlock<>& block, int threadidx,
        Eigen::Ref<Eigen::MatrixXd> native_gradient) const {
        m_impl->finalize_vjp(block, threadidx, native_gradient);
    }

    template <int NSTOKES>
    bool SuccessiveOrdersSource<NSTOKES>::supports_linearization(
        sasktran2::LinearizationMode mode) const {
        return m_impl->supports_linearization(mode);
    }

    template <int NSTOKES>
    void SuccessiveOrdersSource<NSTOKES>::start_of_ray_source(
        const sasktran2::WavelengthBlock<>& block, int losidx,
        int wavel_threadidx, int,
        sasktran2::WavelengthBlockDual<NSTOKES>& source) const {
        m_impl->start_of_ray_source(block, losidx, wavel_threadidx, source);
    }

    template <int NSTOKES>
    void SuccessiveOrdersSource<NSTOKES>::start_of_ray_source_jvp(
        int wavelength, int losidx, int wavel_threadidx, int,
        Eigen::Ref<const Eigen::VectorXd>,
        sasktran2::RadianceJVP<NSTOKES>& source) const {
        m_impl->start_of_ray_source_jvp(wavelength, losidx, wavel_threadidx,
                                        source);
    }

    template <int NSTOKES>
    void SuccessiveOrdersSource<NSTOKES>::start_of_ray_source_vjp(
        int wavelength, int losidx, int wavel_threadidx, int threadidx,
        const Eigen::Vector<double, NSTOKES>&,
        Eigen::Vector<double, NSTOKES>& cotangent,
        Eigen::Ref<Eigen::VectorXd>) const {
        m_impl->start_of_ray_source_vjp(wavelength, losidx, wavel_threadidx,
                                        threadidx, cotangent);
    }

    template class SuccessiveOrdersSource<1>;
    template class SuccessiveOrdersSource<3>;

    template <int NSTOKES>
    std::unique_ptr<SourceTermInterface<NSTOKES>> make_successive_orders_source(
        const sasktran2::raytracing::RayTracerBase& raytracer,
        const sasktran2::Geometry1D& geometry) {
        return std::make_unique<SuccessiveOrdersSource<NSTOKES>>(raytracer,
                                                                 geometry);
    }

    template std::unique_ptr<SourceTermInterface<1>>
    make_successive_orders_source<1>(
        const sasktran2::raytracing::RayTracerBase&,
        const sasktran2::Geometry1D&);
    template std::unique_ptr<SourceTermInterface<3>>
    make_successive_orders_source<3>(
        const sasktran2::raytracing::RayTracerBase&,
        const sasktran2::Geometry1D&);

#ifdef SKTRAN_RUST_SUPPORT
    template <int NSTOKES>
    std::unique_ptr<SourceTermInterface<NSTOKES>> make_successive_orders_source(
        const sasktran2::raytracing::RustRayTracer2D& raytracer,
        const sasktran2::Geometry2D& geometry,
        std::shared_ptr<sasktran2::solartransmission::SolarTransmissionTable2D>
            shared_solar_table) {
        return std::make_unique<SuccessiveOrdersSource<NSTOKES>>(
            raytracer, geometry, std::move(shared_solar_table));
    }

    template std::unique_ptr<SourceTermInterface<1>>
    make_successive_orders_source<1>(
        const sasktran2::raytracing::RustRayTracer2D&,
        const sasktran2::Geometry2D&,
        std::shared_ptr<
            sasktran2::solartransmission::SolarTransmissionTable2D>);
    template std::unique_ptr<SourceTermInterface<3>>
    make_successive_orders_source<3>(
        const sasktran2::raytracing::RustRayTracer2D&,
        const sasktran2::Geometry2D&,
        std::shared_ptr<
            sasktran2::solartransmission::SolarTransmissionTable2D>);
#endif

} // namespace sasktran2::successive_orders
