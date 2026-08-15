#pragma once

#include "fixed_point.h"
#include "scattering.h"
#include "transport.h"

#include <Eigen/Core>

#include <cstddef>

namespace sasktran2::successive_orders {

    template <int NSTOKES> struct ProblemParameterData;

    template <> struct ProblemParameterData<1> {
        Eigen::VectorXd forcing;
        Eigen::VectorXd transport_values;
        Eigen::MatrixXd atmospheric_coefficients;
        Eigen::VectorXd ground_values;

        void resize(const TransportOperator& transport,
                    const ScatteringOperator<1>& scattering,
                    bool include_transport = true);
        void set_zero();
        std::size_t storage_bytes() const;
    };

    template <> struct ProblemParameterData<3> {
        Eigen::VectorXd forcing;
        Eigen::VectorXd transport_values;
        Eigen::MatrixXd atmospheric_coefficients;
        Eigen::VectorXd ground_values;

        void resize(const TransportOperator& transport,
                    const ScatteringOperator<3>& scattering,
                    bool include_transport = true);
        void set_zero();
        std::size_t storage_bytes() const;
    };

    template <int NSTOKES> struct ProblemWorkspace {
        Eigen::VectorXd incoming;
        Eigen::VectorXd auxiliary_incoming;
        Eigen::VectorXd auxiliary_state;
        Eigen::VectorXd direct_state;
        ScatteringWorkspace<NSTOKES> scattering;
        FixedPointWorkspace fixed_point;

        void resize(const TransportOperator&,
                    const ScatteringOperator<NSTOKES>& scattering_operator) {
            incoming.resize(scattering_operator.input_size());
            auxiliary_incoming.resize(scattering_operator.input_size());
            auxiliary_state.resize(scattering_operator.output_size());
            direct_state.resize(scattering_operator.output_size());
            scattering = scattering_operator.make_workspace();
        }

        std::size_t storage_bytes() const {
            return static_cast<std::size_t>(
                       incoming.size() + auxiliary_incoming.size() +
                       auxiliary_state.size() + direct_state.size()) *
                       sizeof(double) +
                   scattering.storage_bytes();
        }
    };

    /** Matrix-free successive-orders system
     *
     *     x = S (b + T x).
     *
     * The transport and scattering objects own only the current atmosphere's
     * values. This class owns no large matrices and uses the same operations
     * for the primal, implicit forward product, and implicit adjoint product.
     */
    template <int NSTOKES> class Problem {
      public:
        Problem(TransportOperator& transport,
                ScatteringOperator<NSTOKES>& scattering)
            : m_transport(&transport), m_scattering(&scattering) {
            if (transport.sparsity().rows() * NSTOKES !=
                    scattering.input_size() ||
                transport.sparsity().columns() * NSTOKES !=
                    scattering.output_size()) {
                throw std::invalid_argument(
                    "successive-orders transport/scattering dimensions do not "
                    "match");
            }
        }

        int state_size() const { return m_scattering->output_size(); }
        int incoming_size() const { return m_scattering->input_size(); }

        void apply(Eigen::Ref<const Eigen::VectorXd> state,
                   Eigen::Ref<const Eigen::VectorXd> forcing,
                   Eigen::Ref<Eigen::VectorXd> result,
                   ProblemWorkspace<NSTOKES>& workspace) const {
            validate_state_and_forcing(state, forcing, result);
            m_transport->template apply_stokes<NSTOKES>(state,
                                                        workspace.incoming);
            workspace.incoming += forcing;
            m_scattering->apply(workspace.incoming, result,
                                workspace.scattering);
        }

        void apply_linear(Eigen::Ref<const Eigen::VectorXd> state,
                          Eigen::Ref<Eigen::VectorXd> result,
                          ProblemWorkspace<NSTOKES>& workspace) const {
            if (state.size() != state_size() || result.size() != state_size()) {
                throw std::invalid_argument(
                    "invalid successive-orders linear operator dimensions");
            }
            m_transport->template apply_stokes<NSTOKES>(state,
                                                        workspace.incoming);
            m_scattering->apply(workspace.incoming, result,
                                workspace.scattering);
        }

        void
        apply_linear_transpose(Eigen::Ref<const Eigen::VectorXd> state,
                               Eigen::Ref<Eigen::VectorXd> result,
                               ProblemWorkspace<NSTOKES>& workspace) const {
            if (state.size() != state_size() || result.size() != state_size()) {
                throw std::invalid_argument(
                    "invalid successive-orders transpose operator dimensions");
            }
            m_scattering->apply_transpose(state, workspace.incoming,
                                          workspace.scattering);
            m_transport->template apply_transpose_stokes<NSTOKES>(
                workspace.incoming, result);
        }

        FixedPointDiagnostics
        solve(Eigen::Ref<const Eigen::VectorXd> forcing, Eigen::VectorXd& state,
              const FixedPointSettings& settings,
              ProblemWorkspace<NSTOKES>& workspace) const {
            prepare_workspace(workspace);
            if (state.size() != state_size()) {
                state.setZero(state_size());
            }
            return FixedPointSolver::solve(
                state,
                [&](const Eigen::VectorXd& input, Eigen::VectorXd& output) {
                    apply(input, forcing, output, workspace);
                },
                settings, workspace.fixed_point);
        }

        FixedPointDiagnostics solve_jvp(
            Eigen::Ref<const Eigen::VectorXd> forcing,
            Eigen::Ref<const Eigen::VectorXd> state,
            const ProblemParameterData<NSTOKES>& tangent,
            Eigen::VectorXd& state_tangent, const FixedPointSettings& settings,
            ProblemWorkspace<NSTOKES>& workspace,
            const Eigen::VectorXd* direct_transport_tangent = nullptr) const {
            prepare_workspace(workspace);
            validate_parameter_data(tangent);
            if (state.size() != state_size() ||
                forcing.size() != incoming_size()) {
                throw std::invalid_argument(
                    "invalid successive-orders JVP primal dimensions");
            }

            // Direct parameter effect: S_p(b + T x) + S(db + dT x).
            m_transport->template apply_stokes<NSTOKES>(state,
                                                        workspace.incoming);
            workspace.incoming += forcing;
            if (direct_transport_tangent != nullptr) {
                if (direct_transport_tangent->size() != incoming_size()) {
                    throw std::invalid_argument(
                        "invalid direct successive-orders transport JVP "
                        "size");
                }
                workspace.auxiliary_incoming = *direct_transport_tangent;
            } else {
                if (tangent.transport_values.size() !=
                    m_transport->values().size()) {
                    throw std::invalid_argument(
                        "missing successive-orders transport JVP values");
                }
                m_transport->template apply_value_jvp_stokes<NSTOKES>(
                    state, tangent.transport_values,
                    workspace.auxiliary_incoming);
            }
            workspace.auxiliary_incoming += tangent.forcing;
            scattering_jvp(workspace.incoming, workspace.auxiliary_incoming,
                           tangent, workspace.direct_state, workspace);

            if (workspace.direct_state.isZero(0.0)) {
                state_tangent.setZero(state_size());
                FixedPointDiagnostics diagnostics;
                diagnostics.termination = FixedPointTermination::tolerance;
                diagnostics.residual_norm = 0.0;
                return diagnostics;
            }
            // The first fixed-point image of a zero tangent is exactly the
            // direct parameter response. Start there and avoid spending an
            // iteration rediscovering a value already available.
            state_tangent = workspace.direct_state;
            return FixedPointSolver::solve(
                state_tangent,
                [&](const Eigen::VectorXd& input, Eigen::VectorXd& output) {
                    apply_linear(input, output, workspace);
                    output += workspace.direct_state;
                },
                settings, workspace.fixed_point);
        }

        FixedPointDiagnostics
        solve_vjp(Eigen::Ref<const Eigen::VectorXd> forcing,
                  Eigen::Ref<const Eigen::VectorXd> state,
                  Eigen::Ref<const Eigen::VectorXd> state_cotangent,
                  ProblemParameterData<NSTOKES>& gradient,
                  const FixedPointSettings& settings,
                  ProblemWorkspace<NSTOKES>& workspace,
                  Eigen::VectorXd& adjoint,
                  bool materialize_transport_gradient = true) const {
            prepare_workspace(workspace);
            gradient.resize(*m_transport, *m_scattering,
                            materialize_transport_gradient);
            gradient.set_zero();
            if (state.size() != state_size() ||
                state_cotangent.size() != state_size() ||
                forcing.size() != incoming_size()) {
                throw std::invalid_argument(
                    "invalid successive-orders VJP primal dimensions");
            }

            // As in the forward product, the adjoint right-hand side is the
            // first image of a zero initial iterate.
            adjoint = state_cotangent;
            const auto diagnostics = FixedPointSolver::solve(
                adjoint,
                [&](const Eigen::VectorXd& input, Eigen::VectorXd& output) {
                    apply_linear_transpose(input, output, workspace);
                    output += state_cotangent;
                },
                settings, workspace.fixed_point);

            m_transport->template apply_stokes<NSTOKES>(state,
                                                        workspace.incoming);
            workspace.incoming += forcing;
            scattering_vjp(workspace.incoming, adjoint,
                           workspace.auxiliary_incoming, gradient, workspace);
            if (materialize_transport_gradient) {
                m_transport->template apply_vjp_stokes<NSTOKES>(
                    state, workspace.auxiliary_incoming,
                    workspace.auxiliary_state, gradient.transport_values);
            }
            gradient.forcing = workspace.auxiliary_incoming;
            return diagnostics;
        }

      private:
        void prepare_workspace(ProblemWorkspace<NSTOKES>& workspace) const {
            if (workspace.incoming.size() != incoming_size() ||
                workspace.direct_state.size() != state_size()) {
                workspace.resize(*m_transport, *m_scattering);
            }
        }

        void validate_state_and_forcing(
            Eigen::Ref<const Eigen::VectorXd> state,
            Eigen::Ref<const Eigen::VectorXd> forcing,
            Eigen::Ref<const Eigen::VectorXd> result) const {
            if (state.size() != state_size() ||
                forcing.size() != incoming_size() ||
                result.size() != state_size()) {
                throw std::invalid_argument(
                    "invalid successive-orders problem dimensions");
            }
        }

        void validate_parameter_data(
            const ProblemParameterData<NSTOKES>& data) const;
        void scattering_jvp(Eigen::Ref<const Eigen::VectorXd> incoming,
                            Eigen::Ref<const Eigen::VectorXd> incoming_tangent,
                            const ProblemParameterData<NSTOKES>& tangent,
                            Eigen::Ref<Eigen::VectorXd> state_tangent,
                            ProblemWorkspace<NSTOKES>& workspace) const;
        void scattering_vjp(Eigen::Ref<const Eigen::VectorXd> incoming,
                            Eigen::Ref<const Eigen::VectorXd> state_cotangent,
                            Eigen::Ref<Eigen::VectorXd> incoming_cotangent,
                            ProblemParameterData<NSTOKES>& gradient,
                            ProblemWorkspace<NSTOKES>& workspace) const;

        TransportOperator* m_transport;
        ScatteringOperator<NSTOKES>* m_scattering;
    };

    template <>
    void Problem<1>::validate_parameter_data(
        const ProblemParameterData<1>& data) const;
    template <>
    void Problem<3>::validate_parameter_data(
        const ProblemParameterData<3>& data) const;
    template <>
    void Problem<1>::scattering_jvp(
        Eigen::Ref<const Eigen::VectorXd> incoming,
        Eigen::Ref<const Eigen::VectorXd> incoming_tangent,
        const ProblemParameterData<1>& tangent,
        Eigen::Ref<Eigen::VectorXd> state_tangent,
        ProblemWorkspace<1>& workspace) const;
    template <>
    void Problem<3>::scattering_jvp(
        Eigen::Ref<const Eigen::VectorXd> incoming,
        Eigen::Ref<const Eigen::VectorXd> incoming_tangent,
        const ProblemParameterData<3>& tangent,
        Eigen::Ref<Eigen::VectorXd> state_tangent,
        ProblemWorkspace<3>& workspace) const;
    template <>
    void Problem<1>::scattering_vjp(
        Eigen::Ref<const Eigen::VectorXd> incoming,
        Eigen::Ref<const Eigen::VectorXd> state_cotangent,
        Eigen::Ref<Eigen::VectorXd> incoming_cotangent,
        ProblemParameterData<1>& gradient,
        ProblemWorkspace<1>& workspace) const;
    template <>
    void Problem<3>::scattering_vjp(
        Eigen::Ref<const Eigen::VectorXd> incoming,
        Eigen::Ref<const Eigen::VectorXd> state_cotangent,
        Eigen::Ref<Eigen::VectorXd> incoming_cotangent,
        ProblemParameterData<3>& gradient,
        ProblemWorkspace<3>& workspace) const;

    extern template class Problem<1>;
    extern template class Problem<3>;

} // namespace sasktran2::successive_orders
