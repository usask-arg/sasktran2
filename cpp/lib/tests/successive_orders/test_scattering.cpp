#include "../../successive_orders/scattering.h"

#include <sasktran2/math/unitsphere.h>
#include <sasktran2/test_helper.h>

#include <array>
#include <cstring>
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

    sasktran2::successive_orders::ScatteringOperator<1>
    scalar_point_basis_operator(bool shared_synthesis = false) {
        constexpr int atmospheric_points = 2;
        constexpr int ground_points = 1;
        constexpr int num_coefficients = 5;
        sasktran2::math::LebedevSphere incoming(14);
        sasktran2::math::LebedevSphere outgoing(26);
        std::vector<std::shared_ptr<
            const sasktran2::successive_orders::ScalarAngularBasis>>
            bases;
        for (int point = 0; point < atmospheric_points; ++point) {
            bases.push_back(
                std::make_shared<
                    const sasktran2::successive_orders::ScalarAngularBasis>(
                    incoming, outgoing, num_coefficients));
        }
        sasktran2::successive_orders::ScatteringBlockLayout layout(
            atmospheric_points, ground_points, incoming.num_points(),
            outgoing.num_points(), 3, 2, 1);
        sasktran2::successive_orders::ScatteringOperator<1> scattering(
            std::move(layout), std::move(bases), shared_synthesis);

        Eigen::MatrixXd coefficients(atmospheric_points, num_coefficients);
        coefficients << 1.0, 0.35, -0.2, 0.08, 0.01, 0.9, -0.15, 0.12, 0.03,
            -0.005;
        scattering.set_atmospheric_coefficients(coefficients);
        Eigen::MatrixXd ground(2, 3);
        ground << 0.2, -0.1, 0.5, 0.4, 0.3, -0.25;
        scattering.set_ground_block(0, ground);
        return scattering;
    }

    sasktran2::successive_orders::ScatteringOperator<3>
    vector_point_basis_operator() {
        constexpr int atmospheric_points = 2;
        constexpr int ground_points = 1;
        constexpr int num_coefficients = 3;
        sasktran2::math::LebedevSphere incoming(6);
        sasktran2::math::LebedevSphere outgoing(14);
        std::vector<std::shared_ptr<
            const sasktran2::successive_orders::VectorAngularBasis>>
            bases;
        for (int point = 0; point < atmospheric_points; ++point) {
            bases.push_back(
                std::make_shared<
                    const sasktran2::successive_orders::VectorAngularBasis>(
                    incoming, outgoing, num_coefficients));
        }
        sasktran2::successive_orders::ScatteringBlockLayout layout(
            atmospheric_points, ground_points, incoming.num_points(),
            outgoing.num_points(), 1, 2, 3);
        sasktran2::successive_orders::ScatteringOperator<3> scattering(
            std::move(layout), std::move(bases));
        scattering.set_atmospheric_coefficients(
            Eigen::MatrixXd::Random(atmospheric_points, 4 * num_coefficients));
        scattering.set_ground_block(0, Eigen::MatrixXd::Random(6, 3));
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

TEST_CASE("Scalar successive-orders point angular bases preserve forward and "
          "reverse products",
          "[successive_orders][scattering]") {
    auto scattering = scalar_point_basis_operator();
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

    Eigen::VectorXd outgoing(scattering.output_size());
    Eigen::VectorXd incoming_transpose(scattering.input_size());
    scattering.apply(incoming, outgoing, workspace);
    scattering.apply_transpose(outgoing_cotangent, incoming_transpose,
                               workspace);
    REQUIRE(outgoing.dot(outgoing_cotangent) ==
            Catch::Approx(incoming.dot(incoming_transpose)).epsilon(2.0e-12));

    Eigen::VectorXd outgoing_tangent(scattering.output_size());
    scattering.apply_jvp(incoming, incoming_tangent, coefficient_tangent,
                         ground_tangent, outgoing_tangent, workspace);
    Eigen::VectorXd incoming_gradient(scattering.input_size());
    Eigen::MatrixXd coefficient_gradient(
        scattering.atmospheric_coefficients().rows(),
        scattering.atmospheric_coefficients().cols());
    Eigen::VectorXd ground_gradient(scattering.ground_value_size());
    scattering.apply_vjp(incoming, outgoing_cotangent, incoming_gradient,
                         coefficient_gradient, ground_gradient, workspace);
    REQUIRE(
        outgoing_tangent.dot(outgoing_cotangent) ==
        Catch::Approx(
            incoming_tangent.dot(incoming_gradient) +
            (coefficient_tangent.array() * coefficient_gradient.array()).sum() +
            ground_tangent.dot(ground_gradient))
            .epsilon(3.0e-12));
}

TEST_CASE("Scalar scattering scratch follows reused primal and derivative "
          "batches without changing products",
          "[successive_orders][scattering][workspace]") {
    const auto require_identical = [](const auto& actual,
                                      const auto& expected) {
        REQUIRE(actual.rows() == expected.rows());
        REQUIRE(actual.cols() == expected.cols());
        REQUIRE(std::memcmp(actual.data(), expected.data(),
                            actual.size() * sizeof(double)) == 0);
    };
    for (const bool shared_synthesis : {false, true}) {
        CAPTURE(shared_synthesis);
        auto scattering = scalar_point_basis_operator(shared_synthesis);
        Eigen::MatrixXd coefficients = scattering.atmospheric_coefficients();
        coefficients.rightCols(2).setZero();
        scattering.set_atmospheric_coefficients(coefficients);
        REQUIRE(scattering.active_coefficients() == 3);
        auto reused = scattering.make_workspace();
        const auto& basis = scattering.angular_basis();
        const std::size_t packed_primal_bytes =
            static_cast<std::size_t>(scattering.layout().atmospheric_blocks()) *
            (basis.input_size() + basis.output_size()) * sizeof(double);
        REQUIRE(reused.storage_bytes() == packed_primal_bytes);

        const Eigen::VectorXd incoming =
            Eigen::VectorXd::LinSpaced(scattering.input_size(), -0.3, 0.7);
        const Eigen::VectorXd incoming_tangent =
            Eigen::VectorXd::LinSpaced(scattering.input_size(), 0.2, -0.4);
        const Eigen::VectorXd outgoing_cotangent =
            Eigen::VectorXd::LinSpaced(scattering.output_size(), -0.8, 0.4);
        Eigen::MatrixXd coefficient_tangent = Eigen::MatrixXd::Zero(
            scattering.layout().atmospheric_blocks(), basis.num_coefficients());
        // A parameter product can activate modes beyond the primal's order.
        coefficient_tangent.col(4).setConstant(0.2);
        const Eigen::VectorXd ground_tangent =
            Eigen::VectorXd::Zero(scattering.ground_value_size());
        Eigen::VectorXd actual(scattering.output_size());
        Eigen::VectorXd expected(scattering.output_size());
        {
            auto fresh = scattering.make_workspace();
            scattering.apply(incoming, actual, reused);
            scattering.apply(incoming, expected, fresh);
            require_identical(actual, expected);
        }
        {
            auto fresh = scattering.make_workspace();
            Eigen::VectorXd actual_transpose(scattering.input_size());
            Eigen::VectorXd expected_transpose(scattering.input_size());
            scattering.apply_transpose(outgoing_cotangent, actual_transpose,
                                       reused);
            scattering.apply_transpose(outgoing_cotangent, expected_transpose,
                                       fresh);
            require_identical(actual_transpose, expected_transpose);
        }
        {
            auto fresh = scattering.make_workspace();
            scattering.apply_jvp(incoming, incoming_tangent,
                                 coefficient_tangent, ground_tangent, actual,
                                 reused);
            scattering.apply_jvp(incoming, incoming_tangent,
                                 coefficient_tangent, ground_tangent, expected,
                                 fresh);
            require_identical(actual, expected);
        }
        {
            auto fresh = scattering.make_workspace();
            Eigen::VectorXd actual_input(scattering.input_size());
            Eigen::VectorXd expected_input(scattering.input_size());
            Eigen::MatrixXd actual_coefficients(coefficient_tangent.rows(),
                                                coefficient_tangent.cols());
            Eigen::MatrixXd expected_coefficients(coefficient_tangent.rows(),
                                                  coefficient_tangent.cols());
            Eigen::VectorXd actual_ground(ground_tangent.size());
            Eigen::VectorXd expected_ground(ground_tangent.size());
            scattering.apply_vjp(incoming, outgoing_cotangent, actual_input,
                                 actual_coefficients, actual_ground, reused);
            scattering.apply_vjp(incoming, outgoing_cotangent, expected_input,
                                 expected_coefficients, expected_ground, fresh);
            require_identical(actual_input, expected_input);
            require_identical(actual_coefficients, expected_coefficients);
            require_identical(actual_ground, expected_ground);
        }
        {
            auto fresh = scattering.make_workspace();
            scattering.apply(incoming, actual, reused);
            scattering.apply(incoming, expected, fresh);
            require_identical(actual, expected);
            coefficient_tangent.setZero();
            scattering.apply_jvp(incoming, incoming_tangent,
                                 coefficient_tangent, ground_tangent, actual,
                                 reused);
            scattering.apply(incoming_tangent, expected, fresh);
            require_identical(actual, expected);
        }
    }
}

TEST_CASE("Scalar point input views match packed strided scattering bitwise",
          "[successive_orders][scattering][workspace]") {
    constexpr int points = 5;
    const std::array<std::array<int, 3>, 2> dimensions{std::array{110, 110, 16},
                                                       std::array{14, 26, 5}};
    for (const auto& [incoming_directions, outgoing_directions,
                      coefficients_count] : dimensions) {
        CAPTURE(incoming_directions, outgoing_directions, coefficients_count);
        sasktran2::math::LebedevSphere incoming_sphere(incoming_directions);
        sasktran2::math::LebedevSphere outgoing_sphere(outgoing_directions);
        std::vector<std::shared_ptr<
            const sasktran2::successive_orders::ScalarAngularBasis>>
            bases;
        bases.push_back(std::make_shared<
                        const sasktran2::successive_orders::ScalarAngularBasis>(
            incoming_sphere, outgoing_sphere, coefficients_count));
        for (int point = 1; point < points; ++point) {
            bases.push_back(
                std::make_shared<
                    const sasktran2::successive_orders::ScalarAngularBasis>(
                    incoming_sphere, *bases.front()));
        }
        for (const bool shared_synthesis : {false, true}) {
            CAPTURE(shared_synthesis);
            sasktran2::successive_orders::ScatteringBlockLayout layout(
                points, 1, incoming_directions, outgoing_directions, 3, 2, 1);
            sasktran2::successive_orders::ScatteringOperator<1> scattering(
                std::move(layout), bases, shared_synthesis);
            Eigen::MatrixXd ground(2, 3);
            ground << 0.2, -0.1, 0.5, 0.4, 0.3, -0.25;
            scattering.set_ground_block(0, ground);
            const Eigen::VectorXd incoming =
                Eigen::VectorXd::LinSpaced(scattering.input_size(), -0.9, 0.7);
            auto workspace = scattering.make_workspace();
            REQUIRE(scattering.memory_usage().angular_basis_bytes ==
                    points * bases.front()->analysis_storage_bytes() +
                        bases.front()->synthesis_storage_bytes());
            for (const int active :
                 {1, 3, coefficients_count / 2, coefficients_count}) {
                CAPTURE(active);
                Eigen::MatrixXd coefficients =
                    Eigen::MatrixXd::Zero(points, coefficients_count);
                for (int point = 0; point < points; ++point) {
                    for (int degree = 0; degree < active; ++degree) {
                        coefficients(point, degree) =
                            0.17 * (point + 1) - 0.03 * degree;
                    }
                }
                scattering.set_atmospheric_coefficients(coefficients);
                REQUIRE(scattering.active_coefficients() == active);

                // Independently retain the original column-major packing and
                // strided per-point analysis, including the batched synthesis.
                Eigen::MatrixXd packed_input(points, incoming_directions);
                Eigen::MatrixXd packed_output(points, outgoing_directions);
                for (int point = 0; point < points; ++point) {
                    packed_input.row(point) =
                        incoming
                            .segment(point * incoming_directions,
                                     incoming_directions)
                            .transpose();
                }
                Eigen::MatrixXd moments(points, active * active);
                Eigen::MatrixXd auxiliary_moments;
                for (int point = 0; point < points; ++point) {
                    if (shared_synthesis) {
                        bases[point]->analyze_active(
                            packed_input.middleRows(point, 1),
                            coefficients.middleRows(point, 1), active,
                            auxiliary_moments);
                        moments.row(point) = auxiliary_moments.row(0);
                    } else {
                        bases[point]->apply_active(
                            packed_input.middleRows(point, 1),
                            coefficients.middleRows(point, 1), active,
                            packed_output.middleRows(point, 1),
                            auxiliary_moments);
                    }
                }
                if (shared_synthesis) {
                    bases.front()->synthesize_active(moments, active,
                                                     packed_output);
                }
                Eigen::VectorXd expected(scattering.output_size());
                for (int point = 0; point < points; ++point) {
                    expected.segment(point * outgoing_directions,
                                     outgoing_directions) =
                        packed_output.row(point).transpose();
                }
                expected.tail(2).noalias() =
                    scattering.ground_block(0) * incoming.tail(3);
                Eigen::VectorXd actual(scattering.output_size());
                scattering.apply(incoming, actual, workspace);
                REQUIRE(std::memcmp(actual.data(), expected.data(),
                                    actual.size() * sizeof(double)) == 0);

                const Eigen::VectorXd cotangent = Eigen::VectorXd::LinSpaced(
                    scattering.output_size(), -0.45, 0.81);
                Eigen::MatrixXd packed_cotangent(points, outgoing_directions);
                for (int point = 0; point < points; ++point) {
                    packed_cotangent.row(point) =
                        cotangent
                            .segment(point * outgoing_directions,
                                     outgoing_directions)
                            .transpose();
                }
                Eigen::MatrixXd packed_transpose(points, incoming_directions);
                for (int point = 0; point < points; ++point) {
                    bases[point]->apply_transpose_active(
                        packed_cotangent.middleRows(point, 1),
                        coefficients.middleRows(point, 1), active,
                        packed_transpose.middleRows(point, 1),
                        auxiliary_moments);
                }
                Eigen::VectorXd expected_transpose(scattering.input_size());
                for (int point = 0; point < points; ++point) {
                    expected_transpose.segment(point * incoming_directions,
                                               incoming_directions) =
                        packed_transpose.row(point).transpose();
                }
                expected_transpose.tail(3).noalias() =
                    scattering.ground_block(0).transpose() * cotangent.tail(2);
                Eigen::VectorXd actual_transpose(scattering.input_size());
                scattering.apply_transpose(cotangent, actual_transpose,
                                           workspace);
                REQUIRE(std::memcmp(
                            actual_transpose.data(), expected_transpose.data(),
                            actual_transpose.size() * sizeof(double)) == 0);

                const Eigen::VectorXd incoming_tangent =
                    Eigen::VectorXd::LinSpaced(scattering.input_size(), -0.31,
                                               0.29);
                Eigen::MatrixXd packed_tangent(points, incoming_directions);
                for (int point = 0; point < points; ++point) {
                    packed_tangent.row(point) =
                        incoming_tangent
                            .segment(point * incoming_directions,
                                     incoming_directions)
                            .transpose();
                }
                Eigen::VectorXd ground_tangent = Eigen::VectorXd::LinSpaced(
                    scattering.ground_value_size(), -0.21, 0.33);
                Eigen::Map<const sasktran2::successive_orders::
                               ScatteringOperator<1>::RowMajorMatrix>
                    ground_tangent_block(ground_tangent.data(), 2, 3);
                for (const bool higher_order_tangent : {false, true}) {
                    CAPTURE(higher_order_tangent);
                    Eigen::MatrixXd coefficient_tangent =
                        Eigen::MatrixXd::Zero(points, coefficients_count);
                    const int tangent_active =
                        higher_order_tangent ? coefficients_count : active;
                    for (int point = 0; point < points; ++point) {
                        for (int degree = 0; degree < tangent_active;
                             ++degree) {
                            coefficient_tangent(point, degree) =
                                0.11 * (point + 1) + 0.013 * degree;
                        }
                    }
                    Eigen::MatrixXd tangent_moments;
                    for (int point = 0; point < points; ++point) {
                        bases[point]->apply_jvp_active(
                            packed_input.middleRows(point, 1),
                            packed_tangent.middleRows(point, 1),
                            coefficients.middleRows(point, 1),
                            coefficient_tangent.middleRows(point, 1),
                            tangent_active, packed_output.middleRows(point, 1),
                            auxiliary_moments, tangent_moments);
                    }
                    Eigen::VectorXd expected_jvp(scattering.output_size());
                    for (int point = 0; point < points; ++point) {
                        expected_jvp.segment(point * outgoing_directions,
                                             outgoing_directions) =
                            packed_output.row(point).transpose();
                    }
                    expected_jvp.tail(2).noalias() =
                        scattering.ground_block(0) * incoming_tangent.tail(3);
                    expected_jvp.tail(2).noalias() +=
                        ground_tangent_block * incoming.tail(3);
                    Eigen::VectorXd actual_jvp(scattering.output_size());
                    scattering.apply_jvp(incoming, incoming_tangent,
                                         coefficient_tangent, ground_tangent,
                                         actual_jvp, workspace);
                    REQUIRE(std::memcmp(actual_jvp.data(), expected_jvp.data(),
                                        actual_jvp.size() * sizeof(double)) ==
                            0);
                }

                Eigen::MatrixXd expected_coefficient_gradient(
                    points, coefficients_count);
                Eigen::MatrixXd cotangent_moments;
                for (int point = 0; point < points; ++point) {
                    bases[point]->apply_vjp(
                        packed_input.middleRows(point, 1),
                        coefficients.middleRows(point, 1),
                        packed_cotangent.middleRows(point, 1),
                        packed_transpose.middleRows(point, 1),
                        expected_coefficient_gradient.middleRows(point, 1),
                        auxiliary_moments, cotangent_moments);
                }
                Eigen::VectorXd expected_input_gradient(
                    scattering.input_size());
                for (int point = 0; point < points; ++point) {
                    expected_input_gradient.segment(point * incoming_directions,
                                                    incoming_directions) =
                        packed_transpose.row(point).transpose();
                }
                expected_input_gradient.tail(3).noalias() =
                    scattering.ground_block(0).transpose() * cotangent.tail(2);
                Eigen::VectorXd expected_ground_gradient(
                    scattering.ground_value_size());
                Eigen::Map<sasktran2::successive_orders::ScatteringOperator<
                    1>::RowMajorMatrix>
                    expected_ground_block(expected_ground_gradient.data(), 2,
                                          3);
                expected_ground_block.noalias() =
                    cotangent.tail(2) * incoming.tail(3).transpose();
                Eigen::VectorXd actual_input_gradient(scattering.input_size());
                Eigen::MatrixXd actual_coefficient_gradient(points,
                                                            coefficients_count);
                Eigen::VectorXd actual_ground_gradient(
                    scattering.ground_value_size());
                scattering.apply_vjp(incoming, cotangent, actual_input_gradient,
                                     actual_coefficient_gradient,
                                     actual_ground_gradient, workspace);
                REQUIRE(std::memcmp(actual_input_gradient.data(),
                                    expected_input_gradient.data(),
                                    actual_input_gradient.size() *
                                        sizeof(double)) == 0);
                REQUIRE(std::memcmp(actual_coefficient_gradient.data(),
                                    expected_coefficient_gradient.data(),
                                    actual_coefficient_gradient.size() *
                                        sizeof(double)) == 0);
                REQUIRE(std::memcmp(actual_ground_gradient.data(),
                                    expected_ground_gradient.data(),
                                    actual_ground_gradient.size() *
                                        sizeof(double)) == 0);
            }
        }
    }
}

TEST_CASE("Scalar input-only VJP preserves full configured transpose bitwise",
          "[successive_orders][scattering][vjp][workspace]") {
    using namespace sasktran2::successive_orders;
    constexpr int points = 3;
    const std::array<std::array<int, 3>, 2> dimensions{std::array{110, 110, 16},
                                                       std::array{14, 26, 5}};
    for (const auto& [input_directions, output_directions, degrees] :
         dimensions) {
        CAPTURE(input_directions, output_directions, degrees);
        sasktran2::math::LebedevSphere incoming_sphere(input_directions);
        sasktran2::math::LebedevSphere outgoing_sphere(output_directions);
        auto basis = std::make_shared<const ScalarAngularBasis>(
            incoming_sphere, outgoing_sphere, degrees);
        for (const int basis_layout : {0, 1, 2}) {
            CAPTURE(basis_layout);
            const ScatteringBlockLayout layout(points, 1, input_directions,
                                               output_directions, 3, 2, 1);
            auto make_scattering = [&]() {
                if (basis_layout == 0) {
                    return ScatteringOperator<1>(layout, basis);
                }
                std::vector<std::shared_ptr<const ScalarAngularBasis>> bases;
                bases.push_back(basis);
                for (int point = 1; point < points; ++point) {
                    bases.push_back(std::make_shared<const ScalarAngularBasis>(
                        incoming_sphere, *basis));
                }
                return ScatteringOperator<1>(layout, std::move(bases),
                                             basis_layout == 2);
            };
            auto scattering = make_scattering();
            Eigen::MatrixXd ground(2, 3);
            ground << 0.2, -0.1, 0.5, 0.4, 0.3, -0.25;
            scattering.set_ground_block(0, ground);
            const Eigen::VectorXd incoming =
                Eigen::VectorXd::LinSpaced(scattering.input_size(), -0.9, 0.7);
            const Eigen::VectorXd cotangent = Eigen::VectorXd::LinSpaced(
                scattering.output_size(), -0.45, 0.81);
            auto full_workspace = scattering.make_workspace();
            auto input_workspace = scattering.make_workspace();
            for (const int active : {1, 3, degrees / 2, degrees}) {
                CAPTURE(active);
                Eigen::MatrixXd coefficients =
                    Eigen::MatrixXd::Zero(points, degrees);
                for (int point = 0; point < points; ++point) {
                    for (int degree = 0; degree < degrees; ++degree) {
                        coefficients(point, degree) =
                            degree < active
                                ? 0.17 * (point + 1) - 0.03 * degree
                                : ((point + degree) % 2 == 0 ? 0.0 : -0.0);
                    }
                }
                scattering.set_atmospheric_coefficients(coefficients);
                REQUIRE(scattering.active_coefficients() == active);
                Eigen::VectorXd expected(scattering.input_size());
                Eigen::MatrixXd coefficient_gradient(points, degrees);
                Eigen::VectorXd ground_gradient(scattering.ground_value_size());
                scattering.apply_vjp(incoming, cotangent, expected,
                                     coefficient_gradient, ground_gradient,
                                     full_workspace);
                // Exercise fresh and reused scratch, including an earlier
                // active-width fixed-point transpose on the same workspace.
                for (const bool reuse_scratch : {false, true}) {
                    CAPTURE(reuse_scratch);
                    if (!reuse_scratch) {
                        input_workspace = scattering.make_workspace();
                    }
                    Eigen::VectorXd actual(scattering.input_size());
                    if (reuse_scratch) {
                        scattering.apply_transpose(cotangent, actual,
                                                   input_workspace);
                    }
                    scattering.apply_input_vjp(cotangent, actual,
                                               input_workspace);
                    REQUIRE(std::memcmp(actual.data(), expected.data(),
                                        actual.size() * sizeof(double)) == 0);
                }
            }
            Eigen::VectorXd wrong_size(scattering.input_size() - 1);
            REQUIRE_THROWS_AS(scattering.apply_input_vjp(cotangent, wrong_size,
                                                         input_workspace),
                              std::invalid_argument);
        }
    }
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

TEST_CASE("Vector successive-orders point angular bases preserve forward and "
          "reverse products",
          "[successive_orders][scattering][vector]") {
    auto scattering = vector_point_basis_operator();
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

    Eigen::VectorXd outgoing(scattering.output_size());
    Eigen::VectorXd incoming_transpose(scattering.input_size());
    scattering.apply(incoming, outgoing, workspace);
    scattering.apply_transpose(outgoing_cotangent, incoming_transpose,
                               workspace);
    REQUIRE(outgoing.dot(outgoing_cotangent) ==
            Catch::Approx(incoming.dot(incoming_transpose)).epsilon(2.0e-12));

    Eigen::VectorXd outgoing_tangent(scattering.output_size());
    scattering.apply_jvp(incoming, incoming_tangent, coefficient_tangent,
                         ground_tangent, outgoing_tangent, workspace);
    Eigen::VectorXd incoming_gradient(scattering.input_size());
    Eigen::MatrixXd coefficient_gradient(
        scattering.atmospheric_coefficients().rows(),
        scattering.atmospheric_coefficients().cols());
    Eigen::VectorXd ground_gradient(scattering.ground_value_size());
    scattering.apply_vjp(incoming, outgoing_cotangent, incoming_gradient,
                         coefficient_gradient, ground_gradient, workspace);
    REQUIRE(
        outgoing_tangent.dot(outgoing_cotangent) ==
        Catch::Approx(
            incoming_tangent.dot(incoming_gradient) +
            (coefficient_tangent.array() * coefficient_gradient.array()).sum() +
            ground_tangent.dot(ground_gradient))
            .epsilon(3.0e-12));
    REQUIRE(scattering.memory_usage().angular_basis_bytes > 0);
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
