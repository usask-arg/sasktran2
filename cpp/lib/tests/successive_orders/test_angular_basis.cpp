#include "../../successive_orders/angular_basis.h"

#include <sasktran2/math/unitsphere.h>
#include <sasktran2/math/wigner.h>
#include <sasktran2/test_helper.h>

#include <algorithm>
#include <cmath>

TEST_CASE("Scalar successive-orders coefficient scattering matches Legendre "
          "matrix",
          "[successive_orders][scattering]") {
    constexpr int num_coefficients = 6;
    sasktran2::math::LebedevSphere incoming_sphere(14);
    sasktran2::math::LebedevSphere outgoing_sphere(26);
    sasktran2::successive_orders::ScalarAngularBasis basis(
        incoming_sphere, outgoing_sphere, num_coefficients);

    constexpr int num_blocks = 3;
    Eigen::MatrixXd incoming(num_blocks, incoming_sphere.num_points());
    Eigen::MatrixXd coefficients(num_blocks, num_coefficients);
    for (int point = 0; point < num_blocks; ++point) {
        for (int direction = 0; direction < incoming.cols(); ++direction) {
            incoming(point, direction) = 0.2 + 0.01 * direction + 0.03 * point;
        }
        for (int degree = 0; degree < num_coefficients; ++degree) {
            coefficients(point, degree) =
                std::exp(-0.35 * degree) * (1.0 + 0.05 * point);
        }
    }

    Eigen::MatrixXd actual(num_blocks, outgoing_sphere.num_points());
    Eigen::MatrixXd moments;
    basis.apply(incoming, coefficients, actual, moments);

    Eigen::MatrixXd expected =
        Eigen::MatrixXd::Zero(actual.rows(), actual.cols());
    sasktran2::math::WignerDCalculator legendre(0, 0);
    for (int point = 0; point < num_blocks; ++point) {
        for (int outgoing = 0; outgoing < outgoing_sphere.num_points();
             ++outgoing) {
            for (int incoming_index = 0;
                 incoming_index < incoming_sphere.num_points();
                 ++incoming_index) {
                const double cosine = std::clamp(
                    incoming_sphere.get_quad_position(incoming_index)
                        .dot(outgoing_sphere.get_quad_position(outgoing)),
                    -1.0, 1.0);
                double phase = 0.0;
                for (int degree = 0; degree < num_coefficients; ++degree) {
                    phase += coefficients(point, degree) *
                             legendre.d(std::acos(cosine), degree);
                }
                expected(point, outgoing) +=
                    incoming_sphere.quadrature_weight(incoming_index) * phase *
                    incoming(point, incoming_index);
            }
        }
    }

    INFO("maximum absolute error = "
         << (actual - expected).cwiseAbs().maxCoeff());
    REQUIRE(actual.isApprox(expected, 2.0e-12));
}

TEST_CASE("Scalar successive-orders scattering JVP and VJP are adjoint",
          "[successive_orders][scattering]") {
    constexpr int num_coefficients = 5;
    sasktran2::math::LebedevSphere incoming_sphere(14);
    sasktran2::math::LebedevSphere outgoing_sphere(26);
    sasktran2::successive_orders::ScalarAngularBasis basis(
        incoming_sphere, outgoing_sphere, num_coefficients);

    constexpr int num_blocks = 4;
    Eigen::MatrixXd incoming =
        Eigen::MatrixXd::Random(num_blocks, incoming_sphere.num_points());
    Eigen::MatrixXd incoming_tangent =
        Eigen::MatrixXd::Random(num_blocks, incoming_sphere.num_points());
    Eigen::MatrixXd coefficients =
        Eigen::MatrixXd::Random(num_blocks, num_coefficients);
    Eigen::MatrixXd coefficient_tangent =
        Eigen::MatrixXd::Random(num_blocks, num_coefficients);
    Eigen::MatrixXd output_tangent(num_blocks, outgoing_sphere.num_points());
    Eigen::MatrixXd moments;
    Eigen::MatrixXd tangent_moments;
    basis.apply_jvp(incoming, incoming_tangent, coefficients,
                    coefficient_tangent, output_tangent, moments,
                    tangent_moments);

    Eigen::MatrixXd output_cotangent =
        Eigen::MatrixXd::Random(num_blocks, outgoing_sphere.num_points());
    Eigen::MatrixXd incoming_cotangent(num_blocks,
                                       incoming_sphere.num_points());
    Eigen::MatrixXd coefficient_gradient(num_blocks, num_coefficients);
    Eigen::MatrixXd analyzed;
    Eigen::MatrixXd moment_cotangent;
    basis.apply_vjp(incoming, coefficients, output_cotangent,
                    incoming_cotangent, coefficient_gradient, analyzed,
                    moment_cotangent);

    const double forward =
        (output_tangent.array() * output_cotangent.array()).sum();
    const double reverse =
        (incoming_tangent.array() * incoming_cotangent.array()).sum() +
        (coefficient_tangent.array() * coefficient_gradient.array()).sum();
    REQUIRE(forward == Catch::Approx(reverse).epsilon(2.0e-12));
}
