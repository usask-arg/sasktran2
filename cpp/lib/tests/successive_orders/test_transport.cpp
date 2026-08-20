#include "../../successive_orders/transport.h"

#include <sasktran2/test_helper.h>

namespace {
    sasktran2::successive_orders::TransportSparsity test_sparsity() {
        return {4, {0, 2, 5, 7}, {0, 2, 0, 1, 3, 1, 2}};
    }
} // namespace

TEST_CASE("Successive-orders CSR transport applies forward and transpose",
          "[successive_orders][transport]") {
    const auto sparsity = test_sparsity();
    sasktran2::successive_orders::TransportOperator transport(sparsity);
    transport.values() << 0.2, 0.4, -0.1, 0.3, 0.8, 0.5, -0.25;

    const Eigen::VectorXd state =
        (Eigen::VectorXd(4) << 0.6, -0.2, 0.9, 0.3).finished();
    const Eigen::VectorXd cotangent =
        (Eigen::VectorXd(3) << -0.4, 0.7, 0.1).finished();
    Eigen::VectorXd incoming(3);
    Eigen::VectorXd state_cotangent(4);
    transport.apply(state, incoming);
    transport.apply_transpose(cotangent, state_cotangent);

    REQUIRE(incoming.dot(cotangent) ==
            Catch::Approx(state.dot(state_cotangent)).epsilon(1.0e-14));
}

TEST_CASE("Successive-orders CSR transport JVP and VJP are adjoint",
          "[successive_orders][transport]") {
    const auto sparsity = test_sparsity();
    sasktran2::successive_orders::TransportOperator transport(sparsity);
    transport.values() << 0.2, 0.4, -0.1, 0.3, 0.8, 0.5, -0.25;

    const Eigen::VectorXd state = Eigen::VectorXd::Random(4);
    const Eigen::VectorXd state_tangent = Eigen::VectorXd::Random(4);
    const Eigen::VectorXd value_tangent = Eigen::VectorXd::Random(7);
    const Eigen::VectorXd incoming_cotangent = Eigen::VectorXd::Random(3);
    Eigen::VectorXd incoming_tangent(3);
    Eigen::VectorXd state_cotangent(4);
    Eigen::VectorXd value_gradient(7);
    transport.apply_jvp(state, state_tangent, value_tangent, incoming_tangent);
    transport.apply_vjp(state, incoming_cotangent, state_cotangent,
                        value_gradient);

    const double forward = incoming_tangent.dot(incoming_cotangent);
    const double reverse =
        state_tangent.dot(state_cotangent) + value_tangent.dot(value_gradient);
    REQUIRE(forward == Catch::Approx(reverse).epsilon(1.0e-13));
}

TEST_CASE("Successive-orders CSR transport shares geometry across Stokes",
          "[successive_orders][transport]") {
    const auto sparsity = test_sparsity();
    sasktran2::successive_orders::TransportOperator transport(sparsity);
    transport.values() << 0.2, 0.4, -0.1, 0.3, 0.8, 0.5, -0.25;

    constexpr int nstokes = 3;
    const Eigen::VectorXd state = Eigen::VectorXd::Random(4 * nstokes);
    const Eigen::VectorXd state_tangent = Eigen::VectorXd::Random(4 * nstokes);
    const Eigen::VectorXd value_tangent = Eigen::VectorXd::Random(7);
    const Eigen::VectorXd incoming_cotangent =
        Eigen::VectorXd::Random(3 * nstokes);
    Eigen::VectorXd incoming(3 * nstokes);
    Eigen::VectorXd state_cotangent(4 * nstokes);
    Eigen::VectorXd incoming_tangent(3 * nstokes);
    Eigen::VectorXd value_gradient(7);

    transport.apply_stokes<nstokes>(state, incoming);
    transport.apply_transpose_stokes<nstokes>(incoming_cotangent,
                                              state_cotangent);
    REQUIRE(incoming.dot(incoming_cotangent) ==
            Catch::Approx(state.dot(state_cotangent)).epsilon(1.0e-13));

    transport.apply_jvp_stokes<nstokes>(state, state_tangent, value_tangent,
                                        incoming_tangent);
    transport.apply_vjp_stokes<nstokes>(state, incoming_cotangent,
                                        state_cotangent, value_gradient);
    REQUIRE(incoming_tangent.dot(incoming_cotangent) ==
            Catch::Approx(state_tangent.dot(state_cotangent) +
                          value_tangent.dot(value_gradient))
                .epsilon(1.0e-13));
}
