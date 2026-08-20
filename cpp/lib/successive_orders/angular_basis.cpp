#include "angular_basis.h"

#include <sasktran2/math/scattering.h>
#include <sasktran2/math/unitsphere.h>
#include <sasktran2/math/wigner.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace sasktran2::successive_orders {
    namespace {
        void fill_basis(const sasktran2::math::UnitSphere& sphere,
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
                    for (int degree = order; degree < num_coefficients;
                         ++degree) {
                        const int mode_start = degree * degree;
                        mode_degrees[mode_start] = degree;
                        if (order == 0) {
                            values(direction_index, mode_start) =
                                degrees(degree);
                        } else {
                            const double scale =
                                std::sqrt(2.0) * degrees(degree);
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

        void fill_wigner_degrees(sasktran2::math::WignerDCalculator& calculator,
                                 double theta, int num_coefficients,
                                 std::vector<double>& values) {
#ifdef SKTRAN_RUST_SUPPORT
            // The Rust-backed wrapper takes an output count, matching the
            // vector length used by PhaseHandler.
            calculator.vec_d_emplace(theta, num_coefficients, values.data());
#else
            // The legacy header-only fallback exposes an inclusive maximum
            // degree instead. Avoid relying on that differing convention.
            for (int degree = 0; degree < num_coefficients; ++degree) {
                values[degree] = calculator.d(theta, degree);
            }
#endif
        }
    } // namespace

    ScalarAngularBasis::ScalarAngularBasis(
        const sasktran2::math::UnitSphere& incoming,
        const sasktran2::math::UnitSphere& outgoing, int num_coefficients)
        : m_num_coefficients(num_coefficients) {
        if (num_coefficients < 1 || incoming.num_points() < 1 ||
            outgoing.num_points() < 1) {
            throw std::invalid_argument(
                "invalid scalar successive-orders angular basis dimensions");
        }

        Eigen::MatrixXd incoming_basis;
        fill_basis(incoming, num_coefficients, incoming_basis, m_mode_degrees);
        std::vector<int> outgoing_degrees;
        fill_basis(outgoing, num_coefficients, m_synthesis, outgoing_degrees);
        if (outgoing_degrees != m_mode_degrees) {
            throw std::logic_error(
                "inconsistent successive-orders angular basis ordering");
        }

        m_analysis = incoming_basis.transpose();
        for (int direction = 0; direction < incoming.num_points();
             ++direction) {
            m_analysis.col(direction) *= incoming.quadrature_weight(direction);
        }
    }

    std::size_t ScalarAngularBasis::storage_bytes() const {
        return static_cast<std::size_t>(m_analysis.size() +
                                        m_synthesis.size()) *
                   sizeof(double) +
               m_mode_degrees.capacity() * sizeof(int);
    }

    void ScalarAngularBasis::validate_blocks(
        Eigen::Ref<const Eigen::MatrixXd> incoming,
        Eigen::Ref<const Eigen::MatrixXd> coefficients,
        Eigen::Ref<const Eigen::MatrixXd> outgoing) const {
        if (incoming.rows() != coefficients.rows() ||
            outgoing.rows() != coefficients.rows() ||
            incoming.cols() != input_size() ||
            outgoing.cols() != output_size() ||
            coefficients.cols() != num_coefficients()) {
            throw std::invalid_argument(
                "invalid scalar successive-orders scattering block dimensions");
        }
    }

    void ScalarAngularBasis::multiply_coefficients(
        Eigen::MatrixXd& moments,
        Eigen::Ref<const Eigen::MatrixXd> coefficients) const {
        for (int mode = 0; mode < num_modes(); ++mode) {
            moments.col(mode).array() *=
                coefficients.col(m_mode_degrees[mode]).array();
        }
    }

    void
    ScalarAngularBasis::apply(Eigen::Ref<const Eigen::MatrixXd> incoming,
                              Eigen::Ref<const Eigen::MatrixXd> coefficients,
                              Eigen::Ref<Eigen::MatrixXd> outgoing,
                              Eigen::MatrixXd& moment_workspace) const {
        apply_active(incoming, coefficients, num_coefficients(), outgoing,
                     moment_workspace);
    }

    void ScalarAngularBasis::apply_active(
        Eigen::Ref<const Eigen::MatrixXd> incoming,
        Eigen::Ref<const Eigen::MatrixXd> coefficients, int active_coefficients,
        Eigen::Ref<Eigen::MatrixXd> outgoing,
        Eigen::MatrixXd& moment_workspace) const {
        validate_blocks(incoming, coefficients, outgoing);
        if (active_coefficients < 1 ||
            active_coefficients > num_coefficients()) {
            throw std::invalid_argument(
                "invalid active scalar scattering coefficient count");
        }
        const int active_modes = active_coefficients * active_coefficients;
        moment_workspace.resize(incoming.rows(), active_modes);
        moment_workspace.noalias() =
            incoming * m_analysis.topRows(active_modes).transpose();
        for (int mode = 0; mode < active_modes; ++mode) {
            moment_workspace.col(mode).array() *=
                coefficients.col(m_mode_degrees[mode]).array();
        }
        outgoing.noalias() =
            moment_workspace * m_synthesis.leftCols(active_modes).transpose();
    }

    void ScalarAngularBasis::apply_transpose(
        Eigen::Ref<const Eigen::MatrixXd> outgoing,
        Eigen::Ref<const Eigen::MatrixXd> coefficients,
        Eigen::Ref<Eigen::MatrixXd> incoming,
        Eigen::MatrixXd& moment_workspace) const {
        apply_transpose_active(outgoing, coefficients, num_coefficients(),
                               incoming, moment_workspace);
    }

    void ScalarAngularBasis::apply_transpose_active(
        Eigen::Ref<const Eigen::MatrixXd> outgoing,
        Eigen::Ref<const Eigen::MatrixXd> coefficients, int active_coefficients,
        Eigen::Ref<Eigen::MatrixXd> incoming,
        Eigen::MatrixXd& moment_workspace) const {
        validate_blocks(incoming, coefficients, outgoing);
        if (active_coefficients < 1 ||
            active_coefficients > num_coefficients()) {
            throw std::invalid_argument(
                "invalid active scalar scattering coefficient count");
        }
        const int active_modes = active_coefficients * active_coefficients;
        moment_workspace.resize(outgoing.rows(), active_modes);
        moment_workspace.noalias() =
            outgoing * m_synthesis.leftCols(active_modes);
        for (int mode = 0; mode < active_modes; ++mode) {
            moment_workspace.col(mode).array() *=
                coefficients.col(m_mode_degrees[mode]).array();
        }
        incoming.noalias() =
            moment_workspace * m_analysis.topRows(active_modes);
    }

    void ScalarAngularBasis::apply_jvp(
        Eigen::Ref<const Eigen::MatrixXd> incoming,
        Eigen::Ref<const Eigen::MatrixXd> incoming_tangent,
        Eigen::Ref<const Eigen::MatrixXd> coefficients,
        Eigen::Ref<const Eigen::MatrixXd> coefficient_tangent,
        Eigen::Ref<Eigen::MatrixXd> outgoing_tangent,
        Eigen::MatrixXd& moment_workspace,
        Eigen::MatrixXd& tangent_workspace) const {
        apply_jvp_active(incoming, incoming_tangent, coefficients,
                         coefficient_tangent, num_coefficients(),
                         outgoing_tangent, moment_workspace, tangent_workspace);
    }

    void ScalarAngularBasis::apply_jvp_active(
        Eigen::Ref<const Eigen::MatrixXd> incoming,
        Eigen::Ref<const Eigen::MatrixXd> incoming_tangent,
        Eigen::Ref<const Eigen::MatrixXd> coefficients,
        Eigen::Ref<const Eigen::MatrixXd> coefficient_tangent,
        int active_coefficients, Eigen::Ref<Eigen::MatrixXd> outgoing_tangent,
        Eigen::MatrixXd& moment_workspace,
        Eigen::MatrixXd& tangent_workspace) const {
        validate_blocks(incoming, coefficients, outgoing_tangent);
        if (incoming_tangent.rows() != incoming.rows() ||
            incoming_tangent.cols() != incoming.cols() ||
            coefficient_tangent.rows() != coefficients.rows() ||
            coefficient_tangent.cols() != coefficients.cols()) {
            throw std::invalid_argument(
                "invalid scalar successive-orders scattering JVP dimensions");
        }
        if (active_coefficients < 1 ||
            active_coefficients > num_coefficients()) {
            throw std::invalid_argument(
                "invalid active scalar scattering JVP coefficient count");
        }

        const int active_modes = active_coefficients * active_coefficients;
        moment_workspace.resize(incoming.rows(), active_modes);
        tangent_workspace.resize(incoming.rows(), active_modes);
        moment_workspace.noalias() =
            incoming * m_analysis.topRows(active_modes).transpose();
        tangent_workspace.noalias() =
            incoming_tangent * m_analysis.topRows(active_modes).transpose();
        for (int mode = 0; mode < active_modes; ++mode) {
            const int degree = m_mode_degrees[mode];
            tangent_workspace.col(mode).array() =
                tangent_workspace.col(mode).array() *
                    coefficients.col(degree).array() +
                moment_workspace.col(mode).array() *
                    coefficient_tangent.col(degree).array();
        }
        outgoing_tangent.noalias() =
            tangent_workspace * m_synthesis.leftCols(active_modes).transpose();
    }

    void ScalarAngularBasis::apply_vjp(
        Eigen::Ref<const Eigen::MatrixXd> incoming,
        Eigen::Ref<const Eigen::MatrixXd> coefficients,
        Eigen::Ref<const Eigen::MatrixXd> outgoing_cotangent,
        Eigen::Ref<Eigen::MatrixXd> incoming_cotangent,
        Eigen::Ref<Eigen::MatrixXd> coefficient_gradient,
        Eigen::MatrixXd& analyzed_workspace,
        Eigen::MatrixXd& moment_cotangent_workspace) const {
        validate_blocks(incoming, coefficients, outgoing_cotangent);
        if (incoming_cotangent.rows() != incoming.rows() ||
            incoming_cotangent.cols() != incoming.cols() ||
            coefficient_gradient.rows() != coefficients.rows() ||
            coefficient_gradient.cols() != coefficients.cols()) {
            throw std::invalid_argument(
                "invalid scalar successive-orders scattering VJP dimensions");
        }

        analyzed_workspace.resize(incoming.rows(), num_modes());
        moment_cotangent_workspace.resize(incoming.rows(), num_modes());
        analyzed_workspace.noalias() = incoming * m_analysis.transpose();
        moment_cotangent_workspace.noalias() = outgoing_cotangent * m_synthesis;

        coefficient_gradient.setZero();
        for (int mode = 0; mode < num_modes(); ++mode) {
            const int degree = m_mode_degrees[mode];
            coefficient_gradient.col(degree).array() +=
                analyzed_workspace.col(mode).array() *
                moment_cotangent_workspace.col(mode).array();
            moment_cotangent_workspace.col(mode).array() *=
                coefficients.col(degree).array();
        }
        incoming_cotangent.noalias() = moment_cotangent_workspace * m_analysis;
    }

    std::size_t VectorAngularWorkspace::storage_bytes() const {
        std::size_t values = static_cast<std::size_t>(
            stacked_input.size() + stacked_moments.size() +
            real_products.size() + imaginary_products.size());
        for (const auto& matrix : stokes_input) {
            values += static_cast<std::size_t>(matrix.size());
        }
        for (const auto& matrix : stokes_output) {
            values += static_cast<std::size_t>(matrix.size());
        }
        for (const auto& matrix : moments) {
            values += static_cast<std::size_t>(matrix.size());
        }
        return values * sizeof(double);
    }

    VectorAngularBasis::VectorAngularBasis(
        const sasktran2::math::UnitSphere& incoming,
        const sasktran2::math::UnitSphere& outgoing, int num_coefficients)
        : m_num_coefficients(num_coefficients),
          m_mode_degrees(static_cast<std::size_t>(num_coefficients) *
                             num_coefficients,
                         0) {
        if (num_coefficients < 1 || incoming.num_points() < 1 ||
            outgoing.num_points() < 1) {
            throw std::invalid_argument(
                "invalid vector successive-orders angular basis dimensions");
        }
        const int modes = num_modes();
        for (int channel = 0; channel < 3; ++channel) {
            m_analysis_real[channel].setZero(modes, incoming.num_points());
            m_analysis_imaginary[channel].setZero(modes, incoming.num_points());
            m_synthesis_real[channel].setZero(outgoing.num_points(), modes);
            m_synthesis_imaginary[channel].setZero(outgoing.num_points(),
                                                   modes);
        }

        const std::array<int, 3> spins{0, 2, -2};
        std::vector<Eigen::Vector3d> incoming_directions(incoming.num_points());
        std::vector<Eigen::Vector3d> outgoing_directions(outgoing.num_points());
        for (int index = 0; index < incoming.num_points(); ++index) {
            incoming_directions[index] =
                incoming.get_quad_position(index).normalized();
        }
        for (int index = 0; index < outgoing.num_points(); ++index) {
            outgoing_directions[index] =
                outgoing.get_quad_position(index).normalized();
        }

        // Reuse one recurrence calculator for an entire angular order and
        // evaluate all degrees in one call. With the Rust-backed calculator
        // this also removes tens of thousands of small FFI calls at setup.
        for (int order = 1 - num_coefficients; order < num_coefficients;
             ++order) {
            std::array<sasktran2::math::WignerDCalculator, 3> calculators{
                sasktran2::math::WignerDCalculator(order, spins[0]),
                sasktran2::math::WignerDCalculator(order, spins[1]),
                sasktran2::math::WignerDCalculator(order, spins[2])};
            std::array<std::vector<double>, 3> degree_values;
            for (auto& values : degree_values) {
                values.resize(num_coefficients);
            }
            for (int direction_index = 0;
                 direction_index < incoming.num_points(); ++direction_index) {
                const auto& direction = incoming_directions[direction_index];
                const double theta =
                    std::acos(std::clamp(direction.z(), -1.0, 1.0));
                const double phi = std::atan2(direction.y(), direction.x());
                const double cosine = std::cos(order * phi);
                const double sine = std::sin(order * phi);
                const double weight =
                    incoming.quadrature_weight(direction_index);
                for (int channel = 0; channel < 3; ++channel) {
                    fill_wigner_degrees(calculators[channel], theta,
                                        num_coefficients,
                                        degree_values[channel]);
                }
                for (int degree = std::abs(order); degree < num_coefficients;
                     ++degree) {
                    const int mode = degree * degree + order + degree;
                    m_mode_degrees[mode] = degree;
                    for (int channel = 0; channel < 3; ++channel) {
                        const double value = degree_values[channel][degree];
                        m_analysis_real[channel](mode, direction_index) =
                            weight * value * cosine;
                        m_analysis_imaginary[channel](mode, direction_index) =
                            -weight * value * sine;
                    }
                }
            }
            for (int direction_index = 0;
                 direction_index < outgoing.num_points(); ++direction_index) {
                const auto& direction = outgoing_directions[direction_index];
                const double theta =
                    std::acos(std::clamp(direction.z(), -1.0, 1.0));
                const double phi = std::atan2(direction.y(), direction.x());
                const double cosine = std::cos(order * phi);
                const double sine = std::sin(order * phi);
                for (int channel = 0; channel < 3; ++channel) {
                    fill_wigner_degrees(calculators[channel], theta,
                                        num_coefficients,
                                        degree_values[channel]);
                }
                for (int degree = std::abs(order); degree < num_coefficients;
                     ++degree) {
                    const int mode = degree * degree + order + degree;
                    for (int channel = 0; channel < 3; ++channel) {
                        const double value = degree_values[channel][degree];
                        m_synthesis_real[channel](direction_index, mode) =
                            value * cosine;
                        m_synthesis_imaginary[channel](direction_index, mode) =
                            value * sine;
                    }
                }
            }
        }

        sasktran2::math::WignerDCalculator d00(0, 0);
        sasktran2::math::WignerDCalculator d22(2, 2);
        sasktran2::math::WignerDCalculator d02(0, 2);
        sasktran2::math::WignerDCalculator d2m2(2, -2);
        for (int output_index = 0; output_index < outgoing.num_points();
             ++output_index) {
            const auto& outgoing_direction = outgoing_directions[output_index];
            for (int input_index = 0; input_index < incoming.num_points();
                 ++input_index) {
                const auto& incoming_direction =
                    incoming_directions[input_index];
                const double direction_cosine = std::clamp(
                    incoming_direction.dot(outgoing_direction), -1.0, 1.0);
                const bool special_frame =
                    std::sqrt(std::max(0.0, 1.0 - direction_cosine *
                                                      direction_cosine)) <
                        1.0e-6 ||
                    incoming_direction.head<2>().norm() < 1.0e-8 ||
                    outgoing_direction.head<2>().norm() < 1.0e-8;
                if (!special_frame) {
                    continue;
                }

                FrameCorrection correction;
                correction.input_index = input_index;
                correction.output_index = output_index;
                correction.values.assign(
                    static_cast<std::size_t>(num_coefficients) * 4 * 9, 0.0);
                double theta = 0.0;
                double c1 = 0.0;
                double c2 = 0.0;
                double s1 = 0.0;
                double s2 = 0.0;
                int negation = 1;
                sasktran2::math::stokes_scattering_factors(
                    -incoming_direction, -outgoing_direction, theta, c1, c2, s1,
                    s2, negation);
                const double weight = incoming.quadrature_weight(input_index);
                bool nonzero = false;
                for (int degree = 0; degree < num_coefficients; ++degree) {
                    const double value00 = d00.d(theta, degree);
                    const double value02 = d02.d(theta, degree);
                    const double value22 = d22.d(theta, degree);
                    const double value2m2 = d2m2.d(theta, degree);
                    const double wigner_add = value22 + value2m2;
                    const double wigner_minus = value22 - value2m2;
                    const int first_mode = degree * degree;
                    const int last_mode = (degree + 1) * (degree + 1);
                    for (int family = 0; family < 4; ++family) {
                        std::array<double, 4> greek{};
                        greek[family] = 1.0;
                        Eigen::Matrix3d expected = Eigen::Matrix3d::Zero();
                        expected(0, 0) = greek[0] * value00;
                        expected(1, 0) = -greek[3] * c2 * value02;
                        expected(2, 0) = -greek[3] * s2 * value02;
                        expected(0, 1) = -greek[3] * c1 * value02;
                        expected(0, 2) = greek[3] * s1 * value02;
                        expected(1, 1) =
                            0.5 * greek[1] *
                                (c1 * c2 * wigner_add -
                                 s1 * s2 * wigner_minus) +
                            0.5 * greek[2] *
                                (c1 * c2 * wigner_minus - s1 * s2 * wigner_add);
                        expected(2, 1) =
                            0.5 * greek[1] *
                                (c1 * s2 * wigner_add +
                                 s1 * c2 * wigner_minus) +
                            0.5 * greek[2] *
                                (c1 * s2 * wigner_minus + s1 * c2 * wigner_add);
                        expected(1, 2) =
                            -0.5 * greek[1] *
                                (s1 * c2 * wigner_add +
                                 c1 * s2 * wigner_minus) -
                            0.5 * greek[2] *
                                (s1 * c2 * wigner_minus + c1 * s2 * wigner_add);
                        expected(2, 2) = 0.5 * greek[1] *
                                             (-s1 * s2 * wigner_add +
                                              c1 * c2 * wigner_minus) +
                                         0.5 * greek[2] *
                                             (-s1 * s2 * wigner_minus +
                                              c1 * c2 * wigner_add);

                        for (int input_stokes = 0; input_stokes < 3;
                             ++input_stokes) {
                            std::array<std::complex<double>, 3> result{};
                            for (int mode = first_mode; mode < last_mode;
                                 ++mode) {
                                const std::complex<double> analysis0(
                                    m_analysis_real[0](mode, input_index),
                                    m_analysis_imaginary[0](mode, input_index));
                                const std::complex<double> analysis1(
                                    m_analysis_real[1](mode, input_index),
                                    m_analysis_imaginary[1](mode, input_index));
                                const std::complex<double> analysis2(
                                    m_analysis_real[2](mode, input_index),
                                    m_analysis_imaginary[2](mode, input_index));
                                std::complex<double> intensity{};
                                std::complex<double> plus{};
                                std::complex<double> minus{};
                                if (input_stokes == 0) {
                                    intensity = analysis0;
                                } else if (input_stokes == 1) {
                                    plus = analysis1;
                                    minus = analysis2;
                                } else {
                                    plus = std::complex<double>(0.0, 1.0) *
                                           analysis1;
                                    minus = std::complex<double>(0.0, -1.0) *
                                            analysis2;
                                }
                                const std::array<std::complex<double>, 3> mixed{
                                    greek[0] * intensity -
                                        0.5 * greek[3] * (plus + minus),
                                    -greek[3] * intensity +
                                        0.5 * ((greek[1] + greek[2]) * plus +
                                               (greek[1] - greek[2]) * minus),
                                    -greek[3] * intensity +
                                        0.5 * ((greek[1] - greek[2]) * plus +
                                               (greek[1] + greek[2]) * minus)};
                                for (int channel = 0; channel < 3; ++channel) {
                                    const std::complex<double> synthesis(
                                        m_synthesis_real[channel](output_index,
                                                                  mode),
                                        m_synthesis_imaginary[channel](
                                            output_index, mode));
                                    result[channel] +=
                                        synthesis * mixed[channel];
                                }
                            }
                            const std::array<double, 3> factored{
                                result[0].real(),
                                0.5 * (result[1].real() + result[2].real()),
                                0.5 * (result[1].imag() - result[2].imag())};
                            for (int output_stokes = 0; output_stokes < 3;
                                 ++output_stokes) {
                                const int index = ((degree * 4 + family) * 3 +
                                                   output_stokes) *
                                                      3 +
                                                  input_stokes;
                                correction.values[index] =
                                    weight *
                                        expected(output_stokes, input_stokes) -
                                    factored[output_stokes];
                                nonzero = nonzero ||
                                          std::abs(correction.values[index]) >
                                              1.0e-15;
                            }
                        }
                    }
                }
                if (nonzero) {
                    m_frame_corrections.push_back(std::move(correction));
                }
            }
        }
    }

    std::size_t VectorAngularBasis::storage_bytes() const {
        std::size_t values = 0;
        for (int channel = 0; channel < 3; ++channel) {
            values +=
                static_cast<std::size_t>(m_analysis_real[channel].size() +
                                         m_analysis_imaginary[channel].size() +
                                         m_synthesis_real[channel].size() +
                                         m_synthesis_imaginary[channel].size());
        }
        std::size_t correction_bytes =
            m_frame_corrections.capacity() * sizeof(FrameCorrection);
        for (const auto& correction : m_frame_corrections) {
            correction_bytes += correction.values.capacity() * sizeof(double);
        }
        return values * sizeof(double) +
               m_mode_degrees.capacity() * sizeof(int) + correction_bytes;
    }

    void VectorAngularBasis::validate_blocks(
        Eigen::Ref<const Eigen::MatrixXd> incoming,
        Eigen::Ref<const Eigen::MatrixXd> coefficients,
        Eigen::Ref<const Eigen::MatrixXd> outgoing) const {
        if (incoming.rows() != coefficients.rows() ||
            outgoing.rows() != coefficients.rows() ||
            incoming.cols() != input_size() ||
            outgoing.cols() != output_size() ||
            coefficients.cols() != 4 * num_coefficients()) {
            throw std::invalid_argument(
                "invalid vector successive-orders scattering block dimensions");
        }
    }

    void VectorAngularBasis::analyze(Eigen::Ref<const Eigen::MatrixXd> incoming,
                                     int active_modes,
                                     VectorAngularWorkspace& workspace) const {
        const int blocks = static_cast<int>(incoming.rows());
        for (int component = 0; component < 3; ++component) {
            workspace.stokes_input[component].resize(blocks,
                                                     input_directions());
        }
        for (int block = 0; block < blocks; ++block) {
            for (int direction = 0; direction < input_directions();
                 ++direction) {
                for (int component = 0; component < 3; ++component) {
                    workspace.stokes_input[component](block, direction) =
                        incoming(block, 3 * direction + component);
                }
            }
        }

        for (auto& moment : workspace.moments) {
            moment.resize(blocks, active_modes);
        }
        workspace.moments[0].noalias() =
            workspace.stokes_input[0] *
            m_analysis_real[0].topRows(active_modes).transpose();
        workspace.moments[1].noalias() =
            workspace.stokes_input[0] *
            m_analysis_imaginary[0].topRows(active_modes).transpose();

        workspace.stacked_input.resize(2 * blocks, input_directions());
        workspace.stacked_input.topRows(blocks) = workspace.stokes_input[1];
        workspace.stacked_input.bottomRows(blocks) = workspace.stokes_input[2];
        for (int channel = 1; channel < 3; ++channel) {
            workspace.real_products.noalias() =
                workspace.stacked_input *
                m_analysis_real[channel].topRows(active_modes).transpose();
            workspace.imaginary_products.noalias() =
                workspace.stacked_input *
                m_analysis_imaginary[channel].topRows(active_modes).transpose();
            const auto q_real = workspace.real_products.topRows(blocks);
            const auto u_real = workspace.real_products.bottomRows(blocks);
            const auto q_imaginary =
                workspace.imaginary_products.topRows(blocks);
            const auto u_imaginary =
                workspace.imaginary_products.bottomRows(blocks);
            if (channel == 1) {
                workspace.moments[2] = q_real - u_imaginary;
                workspace.moments[3] = q_imaginary + u_real;
            } else {
                workspace.moments[4] = q_real + u_imaginary;
                workspace.moments[5] = q_imaginary - u_real;
            }
        }
    }

    void VectorAngularBasis::apply_active(
        Eigen::Ref<const Eigen::MatrixXd> incoming,
        Eigen::Ref<const Eigen::MatrixXd> coefficients, int active_coefficients,
        Eigen::Ref<Eigen::MatrixXd> outgoing,
        VectorAngularWorkspace& workspace) const {
        validate_blocks(incoming, coefficients, outgoing);
        if (active_coefficients < 1 ||
            active_coefficients > num_coefficients()) {
            throw std::invalid_argument(
                "invalid active vector scattering coefficient count");
        }
        const int blocks = static_cast<int>(incoming.rows());
        const int active_modes = active_coefficients * active_coefficients;
        analyze(incoming, active_modes, workspace);

        for (int block = 0; block < blocks; ++block) {
            for (int mode = 0; mode < active_modes; ++mode) {
                const int degree = m_mode_degrees[mode];
                const double a1 = coefficients(block, 4 * degree);
                const double a2 = coefficients(block, 4 * degree + 1);
                const double a3 = coefficients(block, 4 * degree + 2);
                const double b1 = coefficients(block, 4 * degree + 3);
                const double intensity_real = workspace.moments[0](block, mode);
                const double intensity_imag = workspace.moments[1](block, mode);
                const double plus_real = workspace.moments[2](block, mode);
                const double plus_imag = workspace.moments[3](block, mode);
                const double minus_real = workspace.moments[4](block, mode);
                const double minus_imag = workspace.moments[5](block, mode);
                workspace.moments[0](block, mode) =
                    a1 * intensity_real - 0.5 * b1 * (plus_real + minus_real);
                workspace.moments[1](block, mode) =
                    a1 * intensity_imag - 0.5 * b1 * (plus_imag + minus_imag);
                workspace.moments[2](block, mode) =
                    -b1 * intensity_real +
                    0.5 * ((a2 + a3) * plus_real + (a2 - a3) * minus_real);
                workspace.moments[3](block, mode) =
                    -b1 * intensity_imag +
                    0.5 * ((a2 + a3) * plus_imag + (a2 - a3) * minus_imag);
                workspace.moments[4](block, mode) =
                    -b1 * intensity_real +
                    0.5 * ((a2 - a3) * plus_real + (a2 + a3) * minus_real);
                workspace.moments[5](block, mode) =
                    -b1 * intensity_imag +
                    0.5 * ((a2 - a3) * plus_imag + (a2 + a3) * minus_imag);
            }
        }

        outgoing.setZero();
        workspace.stacked_moments.resize(2 * blocks, active_modes);
        for (int channel = 0; channel < 3; ++channel) {
            workspace.stacked_moments.topRows(blocks) =
                workspace.moments[2 * channel];
            workspace.stacked_moments.bottomRows(blocks) =
                workspace.moments[2 * channel + 1];
            workspace.real_products.noalias() =
                workspace.stacked_moments *
                m_synthesis_real[channel].leftCols(active_modes).transpose();
            workspace.imaginary_products.noalias() =
                workspace.stacked_moments * m_synthesis_imaginary[channel]
                                                .leftCols(active_modes)
                                                .transpose();
            for (int block = 0; block < blocks; ++block) {
                for (int direction = 0; direction < output_directions();
                     ++direction) {
                    const double real =
                        workspace.real_products(block, direction) -
                        workspace.imaginary_products(block + blocks, direction);
                    const double imaginary =
                        workspace.imaginary_products(block, direction) +
                        workspace.real_products(block + blocks, direction);
                    if (channel == 0) {
                        outgoing(block, 3 * direction) = real;
                    } else if (channel == 1) {
                        outgoing(block, 3 * direction + 1) = 0.5 * real;
                        outgoing(block, 3 * direction + 2) = 0.5 * imaginary;
                    } else {
                        outgoing(block, 3 * direction + 1) += 0.5 * real;
                        outgoing(block, 3 * direction + 2) -= 0.5 * imaginary;
                    }
                }
            }
        }
        apply_frame_corrections(incoming, coefficients, active_coefficients,
                                outgoing);
    }

    void VectorAngularBasis::analyze_synthesis_cotangent(
        Eigen::Ref<const Eigen::MatrixXd> outgoing, int active_modes,
        VectorAngularWorkspace& workspace) const {
        const int blocks = static_cast<int>(outgoing.rows());
        for (int component = 0; component < 3; ++component) {
            workspace.stokes_output[component].resize(blocks,
                                                      output_directions());
        }
        for (int block = 0; block < blocks; ++block) {
            for (int direction = 0; direction < output_directions();
                 ++direction) {
                for (int component = 0; component < 3; ++component) {
                    workspace.stokes_output[component](block, direction) =
                        outgoing(block, 3 * direction + component);
                }
            }
        }
        for (auto& moment : workspace.moments) {
            moment.resize(blocks, active_modes);
        }
        workspace.moments[0].noalias() =
            workspace.stokes_output[0] *
            m_synthesis_real[0].leftCols(active_modes);
        workspace.moments[1].noalias() =
            -workspace.stokes_output[0] *
            m_synthesis_imaginary[0].leftCols(active_modes);

        workspace.stacked_input.resize(2 * blocks, output_directions());
        workspace.stacked_input.topRows(blocks) = workspace.stokes_output[1];
        workspace.stacked_input.bottomRows(blocks) = workspace.stokes_output[2];
        for (int channel = 1; channel < 3; ++channel) {
            workspace.real_products.noalias() =
                workspace.stacked_input *
                m_synthesis_real[channel].leftCols(active_modes);
            workspace.imaginary_products.noalias() =
                workspace.stacked_input *
                m_synthesis_imaginary[channel].leftCols(active_modes);
            const auto q_real = workspace.real_products.topRows(blocks);
            const auto u_real = workspace.real_products.bottomRows(blocks);
            const auto q_imaginary =
                workspace.imaginary_products.topRows(blocks);
            const auto u_imaginary =
                workspace.imaginary_products.bottomRows(blocks);
            if (channel == 1) {
                workspace.moments[2] = 0.5 * (q_real + u_imaginary);
                workspace.moments[3] = 0.5 * (u_real - q_imaginary);
            } else {
                workspace.moments[4] = 0.5 * (q_real - u_imaginary);
                workspace.moments[5] = -0.5 * (u_real + q_imaginary);
            }
        }
    }

    void VectorAngularBasis::apply_transpose_active(
        Eigen::Ref<const Eigen::MatrixXd> outgoing,
        Eigen::Ref<const Eigen::MatrixXd> coefficients, int active_coefficients,
        Eigen::Ref<Eigen::MatrixXd> incoming,
        VectorAngularWorkspace& workspace) const {
        validate_blocks(incoming, coefficients, outgoing);
        if (active_coefficients < 1 ||
            active_coefficients > num_coefficients()) {
            throw std::invalid_argument(
                "invalid active vector scattering transpose coefficient count");
        }
        const int blocks = static_cast<int>(outgoing.rows());
        const int active_modes = active_coefficients * active_coefficients;
        analyze_synthesis_cotangent(outgoing, active_modes, workspace);

        for (int block = 0; block < blocks; ++block) {
            for (int mode = 0; mode < active_modes; ++mode) {
                const int degree = m_mode_degrees[mode];
                const double a1 = coefficients(block, 4 * degree);
                const double a2 = coefficients(block, 4 * degree + 1);
                const double a3 = coefficients(block, 4 * degree + 2);
                const double b1 = coefficients(block, 4 * degree + 3);
                const double intensity_real = workspace.moments[0](block, mode);
                const double intensity_imag = workspace.moments[1](block, mode);
                const double plus_real = workspace.moments[2](block, mode);
                const double plus_imag = workspace.moments[3](block, mode);
                const double minus_real = workspace.moments[4](block, mode);
                const double minus_imag = workspace.moments[5](block, mode);
                workspace.moments[0](block, mode) =
                    a1 * intensity_real - b1 * (plus_real + minus_real);
                workspace.moments[1](block, mode) =
                    a1 * intensity_imag - b1 * (plus_imag + minus_imag);
                workspace.moments[2](block, mode) =
                    -0.5 * b1 * intensity_real +
                    0.5 * ((a2 + a3) * plus_real + (a2 - a3) * minus_real);
                workspace.moments[3](block, mode) =
                    -0.5 * b1 * intensity_imag +
                    0.5 * ((a2 + a3) * plus_imag + (a2 - a3) * minus_imag);
                workspace.moments[4](block, mode) =
                    -0.5 * b1 * intensity_real +
                    0.5 * ((a2 - a3) * plus_real + (a2 + a3) * minus_real);
                workspace.moments[5](block, mode) =
                    -0.5 * b1 * intensity_imag +
                    0.5 * ((a2 - a3) * plus_imag + (a2 + a3) * minus_imag);
            }
        }

        incoming.setZero();
        workspace.stacked_moments.resize(2 * blocks, active_modes);
        for (int channel = 0; channel < 3; ++channel) {
            workspace.stacked_moments.topRows(blocks) =
                workspace.moments[2 * channel];
            workspace.stacked_moments.bottomRows(blocks) =
                workspace.moments[2 * channel + 1];
            workspace.real_products.noalias() =
                workspace.stacked_moments *
                m_analysis_real[channel].topRows(active_modes);
            workspace.imaginary_products.noalias() =
                workspace.stacked_moments *
                m_analysis_imaginary[channel].topRows(active_modes);
            for (int block = 0; block < blocks; ++block) {
                for (int direction = 0; direction < input_directions();
                     ++direction) {
                    const double real_real =
                        workspace.real_products(block, direction);
                    const double imaginary_real =
                        workspace.real_products(block + blocks, direction);
                    const double real_imaginary =
                        workspace.imaginary_products(block, direction);
                    const double imaginary_imaginary =
                        workspace.imaginary_products(block + blocks, direction);
                    if (channel == 0) {
                        incoming(block, 3 * direction) =
                            real_real + imaginary_imaginary;
                    } else if (channel == 1) {
                        incoming(block, 3 * direction + 1) +=
                            real_real + imaginary_imaginary;
                        incoming(block, 3 * direction + 2) +=
                            imaginary_real - real_imaginary;
                    } else {
                        incoming(block, 3 * direction + 1) +=
                            real_real + imaginary_imaginary;
                        incoming(block, 3 * direction + 2) +=
                            real_imaginary - imaginary_real;
                    }
                }
            }
        }
        apply_frame_corrections_transpose(outgoing, coefficients,
                                          active_coefficients, incoming);
    }

    void VectorAngularBasis::accumulate_coefficient_vjp(
        Eigen::Ref<const Eigen::MatrixXd> incoming,
        Eigen::Ref<const Eigen::MatrixXd> outgoing_cotangent,
        Eigen::Ref<Eigen::MatrixXd> coefficient_gradient,
        VectorAngularWorkspace& workspace) const {
        if (incoming.rows() != outgoing_cotangent.rows() ||
            incoming.cols() != input_size() ||
            outgoing_cotangent.cols() != output_size() ||
            coefficient_gradient.rows() != incoming.rows() ||
            coefficient_gradient.cols() != 4 * num_coefficients()) {
            throw std::invalid_argument(
                "invalid vector scattering coefficient VJP dimensions");
        }
        const int blocks = static_cast<int>(incoming.rows());
        const int modes = num_modes();
        analyze(incoming, modes, workspace);
        std::array<Eigen::MatrixXd, 6> analyzed = workspace.moments;
        analyze_synthesis_cotangent(outgoing_cotangent, modes, workspace);
        coefficient_gradient.setZero();
        const auto real_inner = [](double left_real, double left_imag,
                                   double right_real, double right_imag) {
            return left_real * right_real + left_imag * right_imag;
        };
        for (int block = 0; block < blocks; ++block) {
            for (int mode = 0; mode < modes; ++mode) {
                const int degree = m_mode_degrees[mode];
                const double ir = analyzed[0](block, mode);
                const double ii = analyzed[1](block, mode);
                const double pr = analyzed[2](block, mode);
                const double pi = analyzed[3](block, mode);
                const double mr = analyzed[4](block, mode);
                const double mi = analyzed[5](block, mode);
                const double iar = workspace.moments[0](block, mode);
                const double iai = workspace.moments[1](block, mode);
                const double par = workspace.moments[2](block, mode);
                const double pai = workspace.moments[3](block, mode);
                const double mar = workspace.moments[4](block, mode);
                const double mai = workspace.moments[5](block, mode);
                coefficient_gradient(block, 4 * degree) +=
                    real_inner(iar, iai, ir, ii);
                coefficient_gradient(block, 4 * degree + 1) +=
                    0.5 * (real_inner(par, pai, pr + mr, pi + mi) +
                           real_inner(mar, mai, pr + mr, pi + mi));
                coefficient_gradient(block, 4 * degree + 2) +=
                    0.5 * (real_inner(par, pai, pr - mr, pi - mi) +
                           real_inner(mar, mai, mr - pr, mi - pi));
                coefficient_gradient(block, 4 * degree + 3) +=
                    real_inner(iar, iai, -0.5 * (pr + mr), -0.5 * (pi + mi)) +
                    real_inner(par, pai, -ir, -ii) +
                    real_inner(mar, mai, -ir, -ii);
            }
        }
        accumulate_frame_correction_vjp(incoming, outgoing_cotangent,
                                        coefficient_gradient);
    }

    void VectorAngularBasis::apply_frame_corrections(
        Eigen::Ref<const Eigen::MatrixXd> incoming,
        Eigen::Ref<const Eigen::MatrixXd> coefficients, int active_coefficients,
        Eigen::Ref<Eigen::MatrixXd> outgoing) const {
        for (int degree = 0; degree < active_coefficients; ++degree) {
            for (int family = 0; family < 4; ++family) {
                const int coefficient_index = 4 * degree + family;
                if (coefficients.col(coefficient_index).isZero(0.0)) {
                    continue;
                }
                const int matrix_start = coefficient_index * 9;
                const auto coefficient_values =
                    coefficients.col(coefficient_index).array();
                for (const auto& correction : m_frame_corrections) {
                    const int input_start = 3 * correction.input_index;
                    const int output_start = 3 * correction.output_index;
                    for (int row = 0; row < 3; ++row) {
                        for (int column = 0; column < 3; ++column) {
                            const double value =
                                correction
                                    .values[matrix_start + 3 * row + column];
                            if (value != 0.0) {
                                outgoing.col(output_start + row).array() +=
                                    value * coefficient_values *
                                    incoming.col(input_start + column).array();
                            }
                        }
                    }
                }
            }
        }
    }

    void VectorAngularBasis::apply_frame_corrections_transpose(
        Eigen::Ref<const Eigen::MatrixXd> outgoing,
        Eigen::Ref<const Eigen::MatrixXd> coefficients, int active_coefficients,
        Eigen::Ref<Eigen::MatrixXd> incoming) const {
        for (int degree = 0; degree < active_coefficients; ++degree) {
            for (int family = 0; family < 4; ++family) {
                const int coefficient_index = 4 * degree + family;
                if (coefficients.col(coefficient_index).isZero(0.0)) {
                    continue;
                }
                const int matrix_start = coefficient_index * 9;
                const auto coefficient_values =
                    coefficients.col(coefficient_index).array();
                for (const auto& correction : m_frame_corrections) {
                    const int input_start = 3 * correction.input_index;
                    const int output_start = 3 * correction.output_index;
                    for (int row = 0; row < 3; ++row) {
                        for (int column = 0; column < 3; ++column) {
                            const double value =
                                correction
                                    .values[matrix_start + 3 * row + column];
                            if (value != 0.0) {
                                incoming.col(input_start + column).array() +=
                                    value * coefficient_values *
                                    outgoing.col(output_start + row).array();
                            }
                        }
                    }
                }
            }
        }
    }

    void VectorAngularBasis::accumulate_frame_correction_vjp(
        Eigen::Ref<const Eigen::MatrixXd> incoming,
        Eigen::Ref<const Eigen::MatrixXd> outgoing_cotangent,
        Eigen::Ref<Eigen::MatrixXd> coefficient_gradient) const {
        for (int degree = 0; degree < num_coefficients(); ++degree) {
            for (int family = 0; family < 4; ++family) {
                const int coefficient_index = 4 * degree + family;
                const int matrix_start = coefficient_index * 9;
                for (const auto& correction : m_frame_corrections) {
                    const int input_start = 3 * correction.input_index;
                    const int output_start = 3 * correction.output_index;
                    for (int row = 0; row < 3; ++row) {
                        for (int column = 0; column < 3; ++column) {
                            const double value =
                                correction
                                    .values[matrix_start + 3 * row + column];
                            if (value != 0.0) {
                                coefficient_gradient.col(coefficient_index)
                                    .array() +=
                                    value *
                                    incoming.col(input_start + column).array() *
                                    outgoing_cotangent.col(output_start + row)
                                        .array();
                            }
                        }
                    }
                }
            }
        }
    }

} // namespace sasktran2::successive_orders
