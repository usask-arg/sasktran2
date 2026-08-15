#include "../../successive_orders/scattering.h"

#include <sasktran2/math/unitsphere.h>
#include <sasktran2/test_helper.h>

#include <stdexcept>
#include <memory>
#include <utility>
#include <vector>

namespace {
    sasktran2::successive_orders::ScatteringOperator<1> scalar_operator() {
        constexpr int atmospheric_points = 2;
        constexpr int ground_points = 1;
        constexpr int num_coefficients = 5;
        sasktran2::math::LebedevSphere incoming(14);
        sasktran2::math::LebedevSphere outgoing(26);
        sasktran2::successive_orders::ScalarAngularBasis basis(
            incoming, outgoing, num_coefficients);
        sasktran2::successive_orders::ScatteringBlockLayout layout(
            atmospheric_points, ground_points, incoming.num_points(),
            outgoing.num_points(), 3, 2, 1);
        sasktran2::successive_orders::ScatteringOperator<1> scattering(
            std::move(layout), std::move(basis));

        Eigen::MatrixXd coefficients(atmospheric_points, num_coefficients);
        coefficients << 1.0, 0.35, -0.2, 0.08, 0.01, 0.9, -0.15, 0.12, 0.03,
            -0.005;
        scattering.set_atmospheric_coefficients(coefficients);
        Eigen::MatrixXd ground(2, 3);
        ground << 0.2, -0.1, 0.5, 0.4, 0.3, -0.25;
        scattering.set_ground_block(0, ground);
        return scattering;
    }
} // namespace

TEST_CASE("Successive-orders scattering layout is point-major with explicit "
          "ground dimensions",
          "[successive_orders][scattering]") {
    const sasktran2::successive_orders::ScatteringBlockLayout layout(2, 1, 4, 5,
                                                                     2, 3, 1);
    const std::vector<int> expected_input_offsets{0, 4, 8, 10};
    const std::vector<int> expected_output_offsets{0, 5, 10, 13};
    REQUIRE(layout.input_offsets() == expected_input_offsets);
    REQUIRE(layout.output_offsets() == expected_output_offsets);
    REQUIRE(layout.atmospheric_blocks() == 2);
    REQUIRE(layout.ground_blocks() == 1);
    REQUIRE(layout.input_directions(2) == 2);
    REQUIRE(layout.output_directions(2) == 3);
}

TEST_CASE("Scalar successive-orders scattering combines coefficient and "
          "dense ground blocks",
          "[successive_orders][scattering]") {
    auto scattering = scalar_operator();
    auto workspace = scattering.make_workspace();

    Eigen::VectorXd incoming(scattering.input_size());
    for (int index = 0; index < incoming.size(); ++index) {
        incoming(index) = 0.15 + 0.02 * index;
    }
    Eigen::VectorXd actual(scattering.output_size());
    scattering.apply(incoming, actual, workspace);

    Eigen::MatrixXd atmospheric_incoming(2, 14);
    atmospheric_incoming.row(0) = incoming.segment(0, 14).transpose();
    atmospheric_incoming.row(1) = incoming.segment(14, 14).transpose();
    Eigen::MatrixXd atmospheric_outgoing(2, 26);
    Eigen::MatrixXd moments;
    scattering.angular_basis().apply(atmospheric_incoming,
                                     scattering.atmospheric_coefficients(),
                                     atmospheric_outgoing, moments);
    Eigen::VectorXd expected(scattering.output_size());
    expected.segment(0, 26) = atmospheric_outgoing.row(0).transpose();
    expected.segment(26, 26) = atmospheric_outgoing.row(1).transpose();
    expected.tail(2) = scattering.ground_block(0) * incoming.tail(3);

    REQUIRE(actual.isApprox(expected, 2.0e-13));
    const auto memory = scattering.memory_usage(&workspace);
    REQUIRE(memory.atmospheric_value_bytes == 2 * 5 * sizeof(double));
    REQUIRE(memory.boundary_value_bytes == 2 * 3 * sizeof(double));
    REQUIRE(memory.angular_basis_bytes > 0);
    REQUIRE(memory.workspace_bytes > 0);
    REQUIRE(memory.total_bytes() ==
            memory.operator_bytes() + memory.workspace_bytes);
}

TEST_CASE("Scalar successive-orders scattering transpose is adjoint",
          "[successive_orders][scattering]") {
    auto scattering = scalar_operator();
    auto workspace = scattering.make_workspace();
    const Eigen::VectorXd incoming =
        Eigen::VectorXd::Random(scattering.input_size());
    const Eigen::VectorXd outgoing_cotangent =
        Eigen::VectorXd::Random(scattering.output_size());
    Eigen::VectorXd outgoing(scattering.output_size());
    Eigen::VectorXd incoming_cotangent(scattering.input_size());

    scattering.apply(incoming, outgoing, workspace);
    scattering.apply_transpose(outgoing_cotangent, incoming_cotangent,
                               workspace);
    REQUIRE(outgoing.dot(outgoing_cotangent) ==
            Catch::Approx(incoming.dot(incoming_cotangent)).epsilon(2.0e-12));
}

TEST_CASE("Scalar successive-orders scattering JVP and VJP include phase and "
          "ground values",
          "[successive_orders][scattering]") {
    auto scattering = scalar_operator();
    auto workspace = scattering.make_workspace();
    const Eigen::VectorXd incoming =
        Eigen::VectorXd::Random(scattering.input_size());
    const Eigen::VectorXd incoming_tangent =
        Eigen::VectorXd::Random(scattering.input_size());
    const Eigen::MatrixXd coefficient_tangent =
        Eigen::MatrixXd::Random(scattering.atmospheric_coefficients().rows(),
                                scattering.atmospheric_coefficients().cols());
    const Eigen::VectorXd ground_tangent =
        Eigen::VectorXd::Random(scattering.ground_value_size());
    const Eigen::VectorXd outgoing_cotangent =
        Eigen::VectorXd::Random(scattering.output_size());

    Eigen::VectorXd outgoing_tangent(scattering.output_size());
    scattering.apply_jvp(incoming, incoming_tangent, coefficient_tangent,
                         ground_tangent, outgoing_tangent, workspace);

    Eigen::VectorXd incoming_cotangent(scattering.input_size());
    Eigen::MatrixXd coefficient_gradient(
        scattering.atmospheric_coefficients().rows(),
        scattering.atmospheric_coefficients().cols());
    Eigen::VectorXd ground_gradient(scattering.ground_value_size());
    scattering.apply_vjp(incoming, outgoing_cotangent, incoming_cotangent,
                         coefficient_gradient, ground_gradient, workspace);

    const double forward = outgoing_tangent.dot(outgoing_cotangent);
    const double reverse =
        incoming_tangent.dot(incoming_cotangent) +
        (coefficient_tangent.array() * coefficient_gradient.array()).sum() +
        ground_tangent.dot(ground_gradient);
    REQUIRE(forward == Catch::Approx(reverse).epsilon(3.0e-12));
}

TEST_CASE("Scalar successive-orders scattering skips inactive parameter "
          "directions exactly",
          "[successive_orders][scattering][jvp]") {
    auto scattering = scalar_operator();
    auto workspace = scattering.make_workspace();
    const Eigen::VectorXd incoming =
        Eigen::VectorXd::Random(scattering.input_size());
    const Eigen::VectorXd incoming_tangent =
        Eigen::VectorXd::Random(scattering.input_size());
    const Eigen::MatrixXd coefficient_tangent =
        Eigen::MatrixXd::Zero(scattering.atmospheric_coefficients().rows(),
                              scattering.atmospheric_coefficients().cols());
    const Eigen::VectorXd ground_tangent =
        Eigen::VectorXd::Zero(scattering.ground_value_size());
    Eigen::VectorXd actual(scattering.output_size());
    Eigen::VectorXd expected(scattering.output_size());

    scattering.apply_jvp(incoming, incoming_tangent, coefficient_tangent,
                         ground_tangent, actual, workspace);
    scattering.apply(incoming_tangent, expected, workspace);
    REQUIRE(actual.isApprox(expected, 2.0e-13));
}

TEST_CASE("Scalar successive-orders scattering detects trailing zero "
          "coefficients",
          "[successive_orders][scattering]") {
    auto scattering = scalar_operator();
    auto coefficients = scattering.atmospheric_coefficients();
    coefficients.rightCols(2).setZero();
    scattering.set_atmospheric_coefficients(coefficients);
    REQUIRE(scattering.active_coefficients() == 3);

    coefficients(0, 4) = 0.01;
    scattering.set_atmospheric_coefficients(coefficients);
    REQUIRE(scattering.active_coefficients() == 5);
}

TEST_CASE("Coefficient vector successive-orders scattering products are dual",
          "[successive_orders][scattering]") {
    sasktran2::math::LebedevSphere incoming_sphere(6);
    sasktran2::math::LebedevSphere outgoing_sphere(14);
    auto basis = std::make_shared<
        const sasktran2::successive_orders::VectorAngularBasis>(
        incoming_sphere, outgoing_sphere, 3);
    sasktran2::successive_orders::ScatteringBlockLayout layout(
        1, 1, incoming_sphere.num_points(), outgoing_sphere.num_points(), 1, 2,
        3);
    sasktran2::successive_orders::ScatteringOperator<3> scattering(
        std::move(layout), std::move(basis));
    scattering.set_atmospheric_coefficients(Eigen::MatrixXd::Random(1, 12));
    scattering.set_ground_block(0, Eigen::MatrixXd::Random(6, 3));
    auto workspace = scattering.make_workspace();

    const Eigen::VectorXd incoming =
        Eigen::VectorXd::Random(scattering.input_size());
    const Eigen::VectorXd incoming_tangent =
        Eigen::VectorXd::Random(scattering.input_size());
    const Eigen::MatrixXd atmospheric_tangent =
        Eigen::MatrixXd::Random(scattering.atmospheric_coefficients().rows(),
                                scattering.atmospheric_coefficients().cols());
    const Eigen::VectorXd ground_tangent =
        Eigen::VectorXd::Random(scattering.ground_value_size());
    const Eigen::VectorXd outgoing_cotangent =
        Eigen::VectorXd::Random(scattering.output_size());

    Eigen::VectorXd outgoing(scattering.output_size());
    Eigen::VectorXd transposed(scattering.input_size());
    scattering.apply(incoming, outgoing, workspace);
    scattering.apply_transpose(outgoing_cotangent, transposed, workspace);
    REQUIRE(outgoing.dot(outgoing_cotangent) ==
            Catch::Approx(incoming.dot(transposed)).epsilon(2.0e-13));

    Eigen::VectorXd outgoing_tangent(scattering.output_size());
    scattering.apply_jvp(incoming, incoming_tangent, atmospheric_tangent,
                         ground_tangent, outgoing_tangent, workspace);
    Eigen::VectorXd incoming_gradient(scattering.input_size());
    Eigen::MatrixXd atmospheric_gradient(
        scattering.atmospheric_coefficients().rows(),
        scattering.atmospheric_coefficients().cols());
    Eigen::VectorXd ground_gradient(scattering.ground_value_size());
    scattering.apply_vjp(incoming, outgoing_cotangent, incoming_gradient,
                         atmospheric_gradient, ground_gradient, workspace);
    REQUIRE(
        outgoing_tangent.dot(outgoing_cotangent) ==
        Catch::Approx(
            incoming_tangent.dot(incoming_gradient) +
            (atmospheric_tangent.array() * atmospheric_gradient.array()).sum() +
            ground_tangent.dot(ground_gradient))
            .epsilon(3.0e-13));
}

TEST_CASE("Successive-orders scattering rejects mismatched dimensions",
          "[successive_orders][scattering]") {
    REQUIRE_THROWS_AS(sasktran2::successive_orders::ScatteringBlockLayout(
                          2, 1, std::vector<int>{0, 2}, std::vector<int>{0, 2}),
                      std::invalid_argument);

    auto scattering = scalar_operator();
    auto workspace = scattering.make_workspace();
    Eigen::VectorXd bad_input(scattering.input_size() - 1);
    Eigen::VectorXd output(scattering.output_size());
    REQUIRE_THROWS_AS(scattering.apply(bad_input, output, workspace),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(scattering.ground_block(1), std::out_of_range);
}
