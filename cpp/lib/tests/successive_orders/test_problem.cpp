#include "../../successive_orders/problem.h"

#include <sasktran2/math/unitsphere.h>
#include <sasktran2/test_helper.h>

#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace {
    using namespace sasktran2::successive_orders;

    TransportSparsity dense_sparsity(int size) {
        std::vector<int> row_offsets(static_cast<std::size_t>(size) + 1);
        std::vector<int> column_indices;
        column_indices.reserve(static_cast<std::size_t>(size) * size);
        for (int row = 0; row < size; ++row) {
            row_offsets[row] = static_cast<int>(column_indices.size());
            for (int column = 0; column < size; ++column) {
                column_indices.push_back(column);
            }
        }
        row_offsets[size] = static_cast<int>(column_indices.size());
        return {size, std::move(row_offsets), std::move(column_indices)};
    }

    ScalarAngularBasis scalar_basis() {
        sasktran2::math::LebedevSphere incoming(6);
        sasktran2::math::LebedevSphere outgoing(6);
        return {incoming, outgoing, 3};
    }

    FixedPointSettings tight_settings() {
        FixedPointSettings settings;
        settings.maximum_iterations = 150;
        settings.relative_tolerance = 1.0e-13;
        settings.absolute_tolerance = 1.0e-14;
        settings.anderson_depth = 3;
        return settings;
    }

    template <int NSTOKES>
    Eigen::MatrixXd build_linear_matrix(const Problem<NSTOKES>& problem,
                                        ProblemWorkspace<NSTOKES>& workspace) {
        Eigen::MatrixXd result(problem.state_size(), problem.state_size());
        Eigen::VectorXd input = Eigen::VectorXd::Zero(problem.state_size());
        Eigen::VectorXd output(problem.state_size());
        for (int column = 0; column < problem.state_size(); ++column) {
            input(column) = 1.0;
            problem.apply_linear(input, output, workspace);
            result.col(column) = output;
            input(column) = 0.0;
        }
        return result;
    }

    struct ScalarProblemFixture {
        ScalarProblemFixture()
            : sparsity(dense_sparsity(8)), transport(sparsity),
              scattering(ScatteringBlockLayout(1, 1, 6, 6, 2, 2, 1),
                         scalar_basis()),
              problem(transport, scattering),
              forcing(Eigen::VectorXd::LinSpaced(8, 0.04, 0.11)) {
            for (int row = 0; row < 8; ++row) {
                for (int column = 0; column < 8; ++column) {
                    transport.values()(row * 8 + column) =
                        0.002 * (1 + (row + 2 * column) % 5);
                }
            }
            scattering.atmospheric_coefficients() << 0.22, 0.035, -0.012;
            Eigen::Matrix2d ground;
            ground << 0.12, 0.025, -0.015, 0.09;
            scattering.set_ground_block(0, ground);
            workspace.resize(transport, scattering);
        }

        ProblemParameterData<1> tangent() const {
            ProblemParameterData<1> result;
            result.resize(transport, scattering);
            result.forcing = Eigen::VectorXd::LinSpaced(8, -0.025, 0.018);
            for (int index = 0; index < result.transport_values.size();
                 ++index) {
                result.transport_values(index) = 0.0007 * ((index % 7) - 3);
            }
            result.atmospheric_coefficients << 0.018, -0.011, 0.006;
            result.ground_values << 0.012, -0.008, 0.004, 0.009;
            return result;
        }

        TransportSparsity sparsity;
        TransportOperator transport;
        ScatteringOperator<1> scattering;
        Problem<1> problem;
        ProblemWorkspace<1> workspace;
        Eigen::VectorXd forcing;
    };

    double parameter_inner_product(const ProblemParameterData<1>& left,
                                   const ProblemParameterData<1>& right) {
        return left.forcing.dot(right.forcing) +
               left.transport_values.dot(right.transport_values) +
               (left.atmospheric_coefficients.array() *
                right.atmospheric_coefficients.array())
                   .sum() +
               left.ground_values.dot(right.ground_values);
    }
} // namespace

TEST_CASE("Successive-orders scalar problem matches a dense fixed-point solve",
          "[successive_orders][problem]") {
    ScalarProblemFixture fixture;
    const Eigen::MatrixXd linear =
        build_linear_matrix(fixture.problem, fixture.workspace);
    Eigen::VectorXd direct(fixture.problem.state_size());
    Eigen::VectorXd zero = Eigen::VectorXd::Zero(fixture.problem.state_size());
    fixture.problem.apply(zero, fixture.forcing, direct, fixture.workspace);
    const Eigen::VectorXd expected =
        (Eigen::MatrixXd::Identity(fixture.problem.state_size(),
                                   fixture.problem.state_size()) -
         linear)
            .partialPivLu()
            .solve(direct);

    Eigen::VectorXd state = zero;
    const auto diagnostics = fixture.problem.solve(
        fixture.forcing, state, tight_settings(), fixture.workspace);
    REQUIRE(diagnostics.converged());
    REQUIRE(state.isApprox(expected, 2.0e-11));
}

TEST_CASE("Successive-orders scalar implicit JVP matches finite differences",
          "[successive_orders][problem][jvp]") {
    ScalarProblemFixture fixture;
    Eigen::VectorXd state = Eigen::VectorXd::Zero(fixture.problem.state_size());
    REQUIRE(
        fixture.problem
            .solve(fixture.forcing, state, tight_settings(), fixture.workspace)
            .converged());
    const ProblemParameterData<1> tangent = fixture.tangent();
    Eigen::VectorXd state_tangent;
    REQUIRE(fixture.problem
                .solve_jvp(fixture.forcing, state, tangent, state_tangent,
                           tight_settings(), fixture.workspace)
                .converged());

    const Eigen::VectorXd transport_values = fixture.transport.values();
    const Eigen::MatrixXd coefficients =
        fixture.scattering.atmospheric_coefficients();
    const Eigen::VectorXd ground_values = fixture.scattering.ground_values();
    constexpr double step = 1.0e-6;

    fixture.transport.values() =
        transport_values + step * tangent.transport_values;
    fixture.scattering.atmospheric_coefficients() =
        coefficients + step * tangent.atmospheric_coefficients;
    fixture.scattering.ground_values() =
        ground_values + step * tangent.ground_values;
    Eigen::VectorXd plus = state;
    const Eigen::VectorXd plus_forcing =
        fixture.forcing + step * tangent.forcing;
    REQUIRE(fixture.problem
                .solve(plus_forcing, plus, tight_settings(), fixture.workspace)
                .converged());

    fixture.transport.values() =
        transport_values - step * tangent.transport_values;
    fixture.scattering.atmospheric_coefficients() =
        coefficients - step * tangent.atmospheric_coefficients;
    fixture.scattering.ground_values() =
        ground_values - step * tangent.ground_values;
    Eigen::VectorXd minus = state;
    const Eigen::VectorXd minus_forcing =
        fixture.forcing - step * tangent.forcing;
    REQUIRE(
        fixture.problem
            .solve(minus_forcing, minus, tight_settings(), fixture.workspace)
            .converged());

    fixture.transport.values() = transport_values;
    fixture.scattering.atmospheric_coefficients() = coefficients;
    fixture.scattering.ground_values() = ground_values;
    const Eigen::VectorXd finite_difference = (plus - minus) / (2.0 * step);
    REQUIRE((state_tangent - finite_difference).norm() <=
            2.0e-7 * std::max(1.0, finite_difference.norm()));
}

TEST_CASE("Successive-orders scalar implicit VJP is adjoint to the JVP",
          "[successive_orders][problem][jvp][vjp]") {
    ScalarProblemFixture fixture;
    Eigen::VectorXd state = Eigen::VectorXd::Zero(fixture.problem.state_size());
    REQUIRE(
        fixture.problem
            .solve(fixture.forcing, state, tight_settings(), fixture.workspace)
            .converged());
    const ProblemParameterData<1> tangent = fixture.tangent();
    Eigen::VectorXd state_tangent;
    REQUIRE(fixture.problem
                .solve_jvp(fixture.forcing, state, tangent, state_tangent,
                           tight_settings(), fixture.workspace)
                .converged());

    const Eigen::VectorXd state_cotangent =
        Eigen::VectorXd::LinSpaced(fixture.problem.state_size(), -0.4, 0.65);
    ProblemParameterData<1> gradient;
    Eigen::VectorXd adjoint;
    REQUIRE(fixture.problem
                .solve_vjp(fixture.forcing, state, state_cotangent, gradient,
                           tight_settings(), fixture.workspace, adjoint)
                .converged());

    const double forward = state_tangent.dot(state_cotangent);
    const double reverse = parameter_inner_product(tangent, gradient);
    REQUIRE(forward == Catch::Approx(reverse).epsilon(3.0e-10));
}

TEST_CASE(
    "Successive-orders coefficient vector problem supports primal products",
    "[successive_orders][problem][vector]") {
    sasktran2::math::LebedevSphere sphere(6);
    auto basis = std::make_shared<const VectorAngularBasis>(sphere, sphere, 3);
    TransportSparsity sparsity = dense_sparsity(sphere.num_points());
    TransportOperator transport(sparsity);
    for (int index = 0; index < transport.values().size(); ++index) {
        transport.values()(index) = 0.002 * (1 + index % 5);
    }
    ScatteringOperator<3> scattering(
        ScatteringBlockLayout(1, 0, sphere.num_points(), sphere.num_points(), 1,
                              1, 3),
        std::move(basis));
    Eigen::MatrixXd coefficients = Eigen::MatrixXd::Zero(1, 12);
    coefficients(0, 0) = 0.18;
    coefficients(0, 4) = 0.025;
    coefficients(0, 5) = -0.012;
    coefficients(0, 6) = 0.009;
    coefficients(0, 7) = 0.006;
    scattering.set_atmospheric_coefficients(coefficients);
    Problem<3> problem(transport, scattering);
    ProblemWorkspace<3> workspace;
    workspace.resize(transport, scattering);
    const Eigen::VectorXd forcing =
        Eigen::VectorXd::LinSpaced(problem.incoming_size(), 0.02, 0.09);

    const Eigen::MatrixXd linear = build_linear_matrix(problem, workspace);
    Eigen::VectorXd direct(problem.state_size());
    Eigen::VectorXd zero = Eigen::VectorXd::Zero(problem.state_size());
    problem.apply(zero, forcing, direct, workspace);
    const Eigen::VectorXd expected =
        (Eigen::MatrixXd::Identity(problem.state_size(), problem.state_size()) -
         linear)
            .partialPivLu()
            .solve(direct);
    Eigen::VectorXd state = zero;
    REQUIRE(
        problem.solve(forcing, state, tight_settings(), workspace).converged());
    REQUIRE(state.isApprox(expected, 2.0e-12));

    ProblemParameterData<3> tangent;
    tangent.resize(transport, scattering);
    tangent.set_zero();
    tangent.forcing =
        Eigen::VectorXd::LinSpaced(problem.incoming_size(), -0.01, 0.015);
    Eigen::VectorXd state_tangent;
    REQUIRE(problem
                .solve_jvp(forcing, state, tangent, state_tangent,
                           tight_settings(), workspace)
                .converged());
    const Eigen::VectorXd state_cotangent =
        Eigen::VectorXd::LinSpaced(problem.state_size(), 0.1, 0.6);
    ProblemParameterData<3> gradient;
    Eigen::VectorXd adjoint;
    REQUIRE(problem
                .solve_vjp(forcing, state, state_cotangent, gradient,
                           tight_settings(), workspace, adjoint)
                .converged());
    REQUIRE(
        state_tangent.dot(state_cotangent) ==
        Catch::Approx(tangent.forcing.dot(gradient.forcing)).epsilon(2.0e-11));
}
