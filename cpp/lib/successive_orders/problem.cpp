#include "problem.h"

namespace sasktran2::successive_orders {
    namespace {
        std::size_t vector_bytes(const Eigen::VectorXd& vector) {
            return static_cast<std::size_t>(vector.size()) * sizeof(double);
        }
        std::size_t matrix_bytes(const Eigen::MatrixXd& matrix) {
            return static_cast<std::size_t>(matrix.size()) * sizeof(double);
        }
    } // namespace

    void
    ProblemParameterData<1>::resize(const TransportOperator& transport,
                                    const ScatteringOperator<1>& scattering,
                                    bool include_transport) {
        forcing.resize(scattering.input_size());
        transport_values.resize(include_transport ? transport.values().size()
                                                  : 0);
        atmospheric_coefficients.resize(
            scattering.atmospheric_coefficients().rows(),
            scattering.atmospheric_coefficients().cols());
        ground_values.resize(scattering.ground_values().size());
    }

    void ProblemParameterData<1>::set_zero() {
        forcing.setZero();
        transport_values.setZero();
        atmospheric_coefficients.setZero();
        ground_values.setZero();
    }

    std::size_t ProblemParameterData<1>::storage_bytes() const {
        return vector_bytes(forcing) + vector_bytes(transport_values) +
               matrix_bytes(atmospheric_coefficients) +
               vector_bytes(ground_values);
    }

    void
    ProblemParameterData<3>::resize(const TransportOperator& transport,
                                    const ScatteringOperator<3>& scattering,
                                    bool include_transport) {
        forcing.resize(scattering.input_size());
        transport_values.resize(include_transport ? transport.values().size()
                                                  : 0);
        atmospheric_coefficients.resize(
            scattering.atmospheric_coefficients().rows(),
            scattering.atmospheric_coefficients().cols());
        ground_values.resize(scattering.ground_values().size());
    }

    void ProblemParameterData<3>::set_zero() {
        forcing.setZero();
        transport_values.setZero();
        atmospheric_coefficients.setZero();
        ground_values.setZero();
    }

    std::size_t ProblemParameterData<3>::storage_bytes() const {
        return vector_bytes(forcing) + vector_bytes(transport_values) +
               matrix_bytes(atmospheric_coefficients) +
               vector_bytes(ground_values);
    }

    template <>
    void Problem<1>::validate_parameter_data(
        const ProblemParameterData<1>& data) const {
        if (data.forcing.size() != incoming_size() ||
            (data.transport_values.size() != 0 &&
             data.transport_values.size() != m_transport->values().size()) ||
            data.atmospheric_coefficients.rows() !=
                m_scattering->atmospheric_coefficients().rows() ||
            data.atmospheric_coefficients.cols() !=
                m_scattering->atmospheric_coefficients().cols() ||
            data.ground_values.size() != m_scattering->ground_values().size()) {
            throw std::invalid_argument("invalid scalar successive-orders "
                                        "parameter product dimensions");
        }
    }

    template <>
    void Problem<3>::validate_parameter_data(
        const ProblemParameterData<3>& data) const {
        if (data.forcing.size() != incoming_size() ||
            (data.transport_values.size() != 0 &&
             data.transport_values.size() != m_transport->values().size()) ||
            data.atmospheric_coefficients.rows() !=
                m_scattering->atmospheric_coefficients().rows() ||
            data.atmospheric_coefficients.cols() !=
                m_scattering->atmospheric_coefficients().cols() ||
            data.ground_values.size() != m_scattering->ground_values().size()) {
            throw std::invalid_argument("invalid vector successive-orders "
                                        "parameter product dimensions");
        }
    }

    template <>
    void Problem<1>::scattering_jvp(
        Eigen::Ref<const Eigen::VectorXd> incoming,
        Eigen::Ref<const Eigen::VectorXd> incoming_tangent,
        const ProblemParameterData<1>& tangent,
        Eigen::Ref<Eigen::VectorXd> state_tangent,
        ProblemWorkspace<1>& workspace) const {
        m_scattering->apply_jvp(
            incoming, incoming_tangent, tangent.atmospheric_coefficients,
            tangent.ground_values, state_tangent, workspace.scattering);
    }

    template <>
    void Problem<3>::scattering_jvp(
        Eigen::Ref<const Eigen::VectorXd> incoming,
        Eigen::Ref<const Eigen::VectorXd> incoming_tangent,
        const ProblemParameterData<3>& tangent,
        Eigen::Ref<Eigen::VectorXd> state_tangent,
        ProblemWorkspace<3>& workspace) const {
        m_scattering->apply_jvp(
            incoming, incoming_tangent, tangent.atmospheric_coefficients,
            tangent.ground_values, state_tangent, workspace.scattering);
    }

    template <>
    void Problem<1>::scattering_vjp(
        Eigen::Ref<const Eigen::VectorXd> incoming,
        Eigen::Ref<const Eigen::VectorXd> state_cotangent,
        Eigen::Ref<Eigen::VectorXd> incoming_cotangent,
        ProblemParameterData<1>& gradient,
        ProblemWorkspace<1>& workspace) const {
        m_scattering->apply_vjp(incoming, state_cotangent, incoming_cotangent,
                                gradient.atmospheric_coefficients,
                                gradient.ground_values, workspace.scattering);
    }

    template <>
    void Problem<3>::scattering_vjp(
        Eigen::Ref<const Eigen::VectorXd> incoming,
        Eigen::Ref<const Eigen::VectorXd> state_cotangent,
        Eigen::Ref<Eigen::VectorXd> incoming_cotangent,
        ProblemParameterData<3>& gradient,
        ProblemWorkspace<3>& workspace) const {
        m_scattering->apply_vjp(incoming, state_cotangent, incoming_cotangent,
                                gradient.atmospheric_coefficients,
                                gradient.ground_values, workspace.scattering);
    }

    template class Problem<1>;
    template class Problem<3>;

} // namespace sasktran2::successive_orders
