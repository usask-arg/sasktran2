#include "../../successive_orders/angular_basis.h"
#include "../../successive_orders/geometry.h"

#include <sasktran2/math/unitsphere.h>
#include <sasktran2/math/wigner.h>
#include <sasktran2/test_helper.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
    // Keep the original scalar basis implementation independent from the
    // production vector recurrence, including repeated calculator.d calls.
    void original_scalar_values(const sasktran2::math::UnitSphere& sphere,
                                int num_coefficients, Eigen::MatrixXd& values,
                                std::vector<int>& mode_degrees) {
        const int num_modes = num_coefficients * num_coefficients;
        values.setZero(sphere.num_points(), num_modes);
        mode_degrees.assign(num_modes, 0);
        Eigen::VectorXd degrees(num_coefficients);
        for (int direction_index = 0; direction_index < sphere.num_points();
             ++direction_index) {
            const Eigen::Vector3d direction =
                sphere.get_quad_position(direction_index).normalized();
            const double theta =
                std::acos(std::clamp(direction.z(), -1.0, 1.0));
            const double phi = std::atan2(direction.y(), direction.x());
            for (int order = 0; order < num_coefficients; ++order) {
                sasktran2::math::WignerDCalculator calculator(order, 0);
                for (int degree = 0; degree < num_coefficients; ++degree) {
                    degrees(degree) = calculator.d(theta, degree);
                }
                for (int degree = order; degree < num_coefficients; ++degree) {
                    const int mode_start = degree * degree;
                    mode_degrees[mode_start] = degree;
                    if (order == 0) {
                        values(direction_index, mode_start) = degrees(degree);
                    } else {
                        const double scale = std::sqrt(2.0) * degrees(degree);
                        const int cosine_mode = mode_start + 2 * order - 1;
                        const int sine_mode = mode_start + 2 * order;
                        mode_degrees[cosine_mode] = degree;
                        mode_degrees[sine_mode] = degree;
                        values(direction_index, cosine_mode) =
                            scale * std::cos(order * phi);
                        values(direction_index, sine_mode) =
                            scale * std::sin(order * phi);
                    }
                }
            }
        }
    }

    struct OriginalScalarBasis {
        Eigen::MatrixXd incoming_values;
        Eigen::MatrixXd analysis;
        Eigen::MatrixXd synthesis;
        std::vector<int> mode_degrees;

        OriginalScalarBasis(const sasktran2::math::UnitSphere& incoming,
                            const sasktran2::math::UnitSphere& outgoing,
                            int num_coefficients) {
            original_scalar_values(incoming, num_coefficients, incoming_values,
                                   mode_degrees);
            std::vector<int> outgoing_degrees;
            original_scalar_values(outgoing, num_coefficients, synthesis,
                                   outgoing_degrees);
            REQUIRE(mode_degrees == outgoing_degrees);
            analysis = incoming_values.transpose();
            for (int direction = 0; direction < incoming.num_points();
                 ++direction) {
                analysis.col(direction) *=
                    incoming.quadrature_weight(direction);
            }
        }

        void multiply(Eigen::MatrixXd& moments,
                      const Eigen::MatrixXd& coefficients) const {
            for (int mode = 0; mode < moments.cols(); ++mode) {
                moments.col(mode).array() *=
                    coefficients.col(mode_degrees[mode]).array();
            }
        }
    };

    void require_same_bits(const Eigen::MatrixXd& actual,
                           const Eigen::MatrixXd& expected) {
        REQUIRE(actual.rows() == expected.rows());
        REQUIRE(actual.cols() == expected.cols());
        REQUIRE(std::memcmp(actual.data(), expected.data(),
                            actual.size() * sizeof(double)) == 0);
    }

    const Eigen::MatrixXd& scalar_synthesis_generation(
        const sasktran2::successive_orders::ScalarAngularBasis& basis) {
        // The internal header defines this identity as m_synthesis.get().
        // Inspect its exact coefficients before products can hide zero signs.
        return *static_cast<const Eigen::MatrixXd*>(
            basis.synthesis_storage_id());
    }

    void require_original_scalar_products(
        const sasktran2::math::UnitSphere& incoming_sphere,
        const sasktran2::math::UnitSphere& outgoing_sphere,
        int num_coefficients) {
        const OriginalScalarBasis original(incoming_sphere, outgoing_sphere,
                                           num_coefficients);
        sasktran2::successive_orders::ScalarAngularBasis basis(
            incoming_sphere, outgoing_sphere, num_coefficients);
        require_same_bits(scalar_synthesis_generation(basis),
                          original.synthesis);
        sasktran2::successive_orders::ScalarAngularBasis incoming_generation(
            incoming_sphere, incoming_sphere, num_coefficients);
        require_same_bits(scalar_synthesis_generation(incoming_generation),
                          original.incoming_values);
        sasktran2::successive_orders::ScalarAngularBasis shared(incoming_sphere,
                                                                basis);

        for (const int points : {1, 3}) {
            CAPTURE(points);
            Eigen::MatrixXd incoming(points, incoming_sphere.num_points());
            Eigen::MatrixXd incoming_tangent(incoming.rows(), incoming.cols());
            Eigen::MatrixXd coefficients(points, num_coefficients);
            Eigen::MatrixXd coefficient_tangent(points, num_coefficients);
            Eigen::MatrixXd cotangent(points, outgoing_sphere.num_points());
            for (Eigen::Index index = 0; index < incoming.size(); ++index) {
                incoming.data()[index] = std::sin(0.019 * (index + 1));
                incoming_tangent.data()[index] = std::cos(0.031 * (index + 2));
            }
            for (Eigen::Index index = 0; index < coefficients.size(); ++index) {
                coefficients.data()[index] = 0.2 - 0.017 * index;
                coefficient_tangent.data()[index] = -0.1 + 0.011 * index;
            }
            for (Eigen::Index index = 0; index < cotangent.size(); ++index) {
                cotangent.data()[index] = std::cos(0.023 * (index + 3));
            }
            const int modes = num_coefficients * num_coefficients;
            Eigen::MatrixXd expected_moments(points, modes);
            expected_moments.noalias() =
                incoming * original.analysis.topRows(modes).transpose();
            const Eigen::MatrixXd analyzed_jvp = expected_moments;
            original.multiply(expected_moments, coefficients);
            Eigen::MatrixXd expected_output(points,
                                            outgoing_sphere.num_points());
            expected_output.noalias() =
                expected_moments *
                original.synthesis.leftCols(modes).transpose();
            Eigen::MatrixXd expected_tangent_moments(points, modes);
            expected_tangent_moments.noalias() =
                incoming_tangent * original.analysis.topRows(modes).transpose();
            for (int mode = 0; mode < expected_moments.cols(); ++mode) {
                const int degree = original.mode_degrees[mode];
                expected_tangent_moments.col(mode).array() =
                    expected_tangent_moments.col(mode).array() *
                        coefficients.col(degree).array() +
                    analyzed_jvp.col(mode).array() *
                        coefficient_tangent.col(degree).array();
            }
            Eigen::MatrixXd expected_output_tangent(
                points, outgoing_sphere.num_points());
            expected_output_tangent.noalias() =
                expected_tangent_moments *
                original.synthesis.leftCols(modes).transpose();
            Eigen::MatrixXd analyzed(points, modes);
            analyzed.noalias() = incoming * original.analysis.transpose();
            Eigen::MatrixXd expected_moment_cotangent(points, modes);
            expected_moment_cotangent.noalias() =
                cotangent * original.synthesis;
            Eigen::MatrixXd expected_coefficient_gradient =
                Eigen::MatrixXd::Zero(points, num_coefficients);
            for (int mode = 0; mode < expected_moments.cols(); ++mode) {
                const int degree = original.mode_degrees[mode];
                expected_coefficient_gradient.col(degree).array() +=
                    analyzed.col(mode).array() *
                    expected_moment_cotangent.col(mode).array();
                expected_moment_cotangent.col(mode).array() *=
                    coefficients.col(degree).array();
            }
            Eigen::MatrixXd expected_input_gradient(points, incoming.cols());
            expected_input_gradient.noalias() =
                expected_moment_cotangent * original.analysis;
            for (const auto* actual_basis : {&basis, &shared}) {
                Eigen::MatrixXd actual_output(points,
                                              outgoing_sphere.num_points());
                Eigen::MatrixXd actual_moments, actual_tangent_moments;
                actual_basis->apply(incoming, coefficients, actual_output,
                                    actual_moments);
                require_same_bits(actual_moments, expected_moments);
                require_same_bits(actual_output, expected_output);
                actual_basis->apply_jvp(incoming, incoming_tangent,
                                        coefficients, coefficient_tangent,
                                        actual_output, actual_moments,
                                        actual_tangent_moments);
                require_same_bits(actual_tangent_moments,
                                  expected_tangent_moments);
                require_same_bits(actual_output, expected_output_tangent);
                Eigen::MatrixXd actual_input_gradient(points, incoming.cols());
                Eigen::MatrixXd actual_coefficient_gradient(points,
                                                            num_coefficients);
                actual_basis->apply_vjp(incoming, coefficients, cotangent,
                                        actual_input_gradient,
                                        actual_coefficient_gradient,
                                        actual_moments, actual_tangent_moments);
                require_same_bits(actual_moments, analyzed);
                require_same_bits(actual_tangent_moments,
                                  expected_moment_cotangent);
                require_same_bits(actual_input_gradient,
                                  expected_input_gradient);
                require_same_bits(actual_coefficient_gradient,
                                  expected_coefficient_gradient);
            }
        }
    }

    class SignedZeroSphere final : public sasktran2::math::UnitSphere {
      public:
        int num_points() const override { return 6; }
        Eigen::Vector3d get_quad_position(int index) const override {
            const std::array<Eigen::Vector3d, 6> positions{
                Eigen::Vector3d{1.0, 0.0, 0.0},
                Eigen::Vector3d{1.0, -0.0, 0.0},
                Eigen::Vector3d{-1.0, 0.0, 0.0},
                Eigen::Vector3d{-1.0, -0.0, 0.0},
                Eigen::Vector3d{0.0, -0.0, 1.0},
                Eigen::Vector3d{0.0, 0.0, -1.0}};
            return positions.at(index);
        }
        double quadrature_weight(int) const override { return 1.0 / 6.0; }
        void interpolate(const Eigen::Vector3d&,
                         std::vector<std::pair<int, double>>&,
                         int&) const override {
            throw std::logic_error(
                "signed-zero basis fixture has no interpolation");
        }
    };
} // namespace

TEST_CASE("Scalar vector Wigner recurrence preserves original basis and "
          "products bitwise",
          "[successive_orders][scattering][wigner]") {
    for (const int directions : {14, 26, 110}) {
        CAPTURE(directions);
        sasktran2::math::LebedevSphere full(directions);
        for (const int coefficients : {1, 3, 8, 16}) {
            CAPTURE(coefficients);
            require_original_scalar_products(full, full, coefficients);
        }
        Eigen::VectorXd altitudes(3);
        altitudes << 0.0, 1000.0, 3000.0;
        sasktran2::Geometry1D geometry(0.4, 0.0, 6372000.0,
                                       std::move(altitudes),
                                       sasktran2::grids::interpolation::linear,
                                       sasktran2::geometrytype::spherical);
        sasktran2::raytracing::SphericalShellRayTracer raytracer(geometry);
        sasktran2::successive_orders::SourceGeometrySettings settings;
        settings.num_incoming = directions;
        settings.num_outgoing = directions;
        settings.use_reduced_horizon_quadrature = true;
        sasktran2::successive_orders::SourceGeometry1D source_geometry(
            raytracer, geometry);
        source_geometry.initialize({}, settings);
        REQUIRE(source_geometry.settings().use_reduced_horizon_quadrature);
        for (const auto& point : source_geometry.source_points()) {
            if (point.is_ground()) {
                continue;
            }
            for (const int coefficients : {1, 3, 8, 16}) {
                CAPTURE(coefficients);
                require_original_scalar_products(point.incoming_sphere(),
                                                 point.outgoing_sphere(),
                                                 coefficients);
            }
        }
    }
}

TEST_CASE("Scalar vector Wigner recurrence preserves sine mode signed zeros",
          "[successive_orders][scattering][wigner]") {
    const SignedZeroSphere sphere;
    Eigen::MatrixXd original_values;
    std::vector<int> degrees;
    original_scalar_values(sphere, 16, original_values, degrees);
    REQUIRE(original_values(0, 3) == 0.0);
    REQUIRE(original_values(1, 3) == 0.0);
    REQUIRE(std::signbit(original_values(0, 3)) !=
            std::signbit(original_values(1, 3)));
    for (const int coefficients : {1, 3, 8, 16}) {
        CAPTURE(coefficients);
        require_original_scalar_products(sphere, sphere, coefficients);
    }
}

TEST_CASE("Scalar outgoing angular transforms share an immutable generation",
          "[successive_orders][scattering]") {
    constexpr int coefficients_count = 6;
    constexpr int points = 3;
    sasktran2::math::LebedevSphere original_incoming(14);
    sasktran2::math::LebedevSphere incoming_sphere(26);
    sasktran2::math::LebedevSphere outgoing_sphere(50);
    auto original =
        std::make_unique<sasktran2::successive_orders::ScalarAngularBasis>(
            original_incoming, outgoing_sphere, coefficients_count);
    sasktran2::successive_orders::ScalarAngularBasis shared(incoming_sphere,
                                                            *original);
    sasktran2::successive_orders::ScalarAngularBasis independent(
        incoming_sphere, outgoing_sphere, coefficients_count);
    REQUIRE(shared.synthesis_storage_id() == original->synthesis_storage_id());
    REQUIRE(shared.synthesis_storage_id() !=
            independent.synthesis_storage_id());
    REQUIRE(shared.storage_bytes() ==
            shared.analysis_storage_bytes() + shared.synthesis_storage_bytes());
    original.reset();

    Eigen::MatrixXd incoming(points, incoming_sphere.num_points());
    Eigen::MatrixXd tangent(incoming.rows(), incoming.cols());
    Eigen::MatrixXd coefficients(points, coefficients_count);
    Eigen::MatrixXd coefficient_tangent(points, coefficients_count);
    Eigen::MatrixXd cotangent(points, outgoing_sphere.num_points());
    for (Eigen::Index i = 0; i < incoming.size(); ++i) {
        incoming.data()[i] = std::sin(0.013 * (i + 1));
        tangent.data()[i] = std::cos(0.017 * (i + 2));
    }
    for (Eigen::Index i = 0; i < coefficients.size(); ++i) {
        coefficients.data()[i] = 0.2 + 0.003 * i;
        coefficient_tangent.data()[i] = -0.02 + 0.001 * i;
    }
    for (Eigen::Index i = 0; i < cotangent.size(); ++i) {
        cotangent.data()[i] = std::cos(0.023 * (i + 3));
    }
    const auto same_bits = [](const Eigen::MatrixXd& left,
                              const Eigen::MatrixXd& right) {
        return left.rows() == right.rows() && left.cols() == right.cols() &&
               std::memcmp(left.data(), right.data(),
                           left.size() * sizeof(double)) == 0;
    };
    Eigen::MatrixXd actual(points, outgoing_sphere.num_points());
    Eigen::MatrixXd expected(actual.rows(), actual.cols());
    Eigen::MatrixXd actual_moments, expected_moments;
    shared.apply(incoming, coefficients, actual, actual_moments);
    independent.apply(incoming, coefficients, expected, expected_moments);
    REQUIRE(same_bits(actual, expected));

    Eigen::MatrixXd actual_tangent, expected_tangent;
    shared.apply_jvp(incoming, tangent, coefficients, coefficient_tangent,
                     actual, actual_moments, actual_tangent);
    independent.apply_jvp(incoming, tangent, coefficients, coefficient_tangent,
                          expected, expected_moments, expected_tangent);
    REQUIRE(same_bits(actual, expected));

    Eigen::MatrixXd actual_input(points, incoming_sphere.num_points());
    Eigen::MatrixXd expected_input(actual_input.rows(), actual_input.cols());
    Eigen::MatrixXd actual_coefficients(points, coefficients_count);
    Eigen::MatrixXd expected_coefficients(points, coefficients_count);
    shared.apply_vjp(incoming, coefficients, cotangent, actual_input,
                     actual_coefficients, actual_moments, actual_tangent);
    independent.apply_vjp(incoming, coefficients, cotangent, expected_input,
                          expected_coefficients, expected_moments,
                          expected_tangent);
    REQUIRE(same_bits(actual_input, expected_input));
    REQUIRE(same_bits(actual_coefficients, expected_coefficients));
}

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
