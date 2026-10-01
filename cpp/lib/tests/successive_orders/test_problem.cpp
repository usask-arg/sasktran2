#include "../../successive_orders/problem.h"

#include <sasktran2.h>
#include <sasktran2/math/unitsphere.h>
#include <sasktran2/test_helper.h>

#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <cstring>
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

    ScalarAngularBasis scalar_basis(int coefficients = 3) {
        sasktran2::math::LebedevSphere incoming(6);
        sasktran2::math::LebedevSphere outgoing(6);
        return {incoming, outgoing, coefficients};
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
        ScalarProblemFixture(int coefficients = 3)
            : sparsity(dense_sparsity(8)), transport(sparsity),
              scattering(ScatteringBlockLayout(1, 1, 6, 6, 2, 2, 1),
                         scalar_basis(coefficients)),
              problem(transport, scattering),
              forcing(Eigen::VectorXd::LinSpaced(8, 0.04, 0.11)) {
            for (int row = 0; row < 8; ++row) {
                for (int column = 0; column < 8; ++column) {
                    transport.values()(row * 8 + column) =
                        0.002 * (1 + (row + 2 * column) % 5);
                }
            }
            auto& phase = scattering.atmospheric_coefficients();
            phase.setZero();
            phase.leftCols(3) << 0.22, 0.035, -0.012;
            scattering.set_active_coefficients(3);
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

TEST_CASE("Successive-orders JVP preserves the primal active phase modes",
          "[successive_orders][problem][jvp][vjp][cache]") {
    ScalarProblemFixture fixture(4);
    auto settings = tight_settings();
    settings.maximum_iterations = 2;
    settings.relative_tolerance = 0.0;
    settings.absolute_tolerance = 0.0;
    const auto require_same = [](const auto& actual, const auto& expected) {
        REQUIRE(actual.rows() == expected.rows());
        REQUIRE(actual.cols() == expected.cols());
        REQUIRE(actual.allFinite());
        if (actual.size() != 0) {
            REQUIRE(std::memcmp(actual.data(), expected.data(),
                                actual.size() * sizeof(double)) == 0);
        }
    };
    const auto& scattering = std::as_const(fixture.scattering);
    REQUIRE(scattering.active_coefficients() == 3);
    REQUIRE(scattering.atmospheric_coefficients()(0, 3) == 0.0);
    const Eigen::MatrixXd coefficients = scattering.atmospheric_coefficients();
    const Eigen::MatrixXd linear =
        build_linear_matrix(fixture.problem, fixture.workspace);
    const Eigen::VectorXd cotangent =
        Eigen::VectorXd::LinSpaced(fixture.problem.state_size(), -0.4, 0.65);
    Eigen::VectorXd transpose_before(fixture.problem.state_size());
    fixture.problem.apply_linear_transpose(cotangent, transpose_before,
                                           fixture.workspace);
    Eigen::VectorXd state = Eigen::VectorXd::Zero(fixture.problem.state_size());
    REQUIRE(fixture.problem
                .solve(fixture.forcing, state, settings, fixture.workspace)
                .iterations == 2);
    const Eigen::VectorXd primal = state;
    ProblemParameterData<1> gradient_before;
    Eigen::VectorXd adjoint_before;
    REQUIRE(fixture.problem
                .solve_vjp(fixture.forcing, state, cotangent, gradient_before,
                           settings, fixture.workspace, adjoint_before)
                .iterations == 2);

    ProblemParameterData<1> tangent;
    tangent.resize(fixture.transport, scattering);
    tangent.set_zero();
    // Degree three belongs to the derivative but is absent from the primal.
    // Product validation must not activate that physical phase coefficient.
    tangent.atmospheric_coefficients(0, 3) = 0.01;
    Eigen::VectorXd state_tangent;
    REQUIRE(fixture.problem
                .solve_jvp(fixture.forcing, state, tangent, state_tangent,
                           settings, fixture.workspace)
                .iterations == 2);
    REQUIRE(scattering.active_coefficients() == 3);
    REQUIRE_FALSE(state_tangent.isZero(0.0));
    require_same(state, primal);
    require_same(scattering.atmospheric_coefficients(), coefficients);
    require_same(build_linear_matrix(fixture.problem, fixture.workspace),
                 linear);
    Eigen::VectorXd transpose_after(fixture.problem.state_size());
    fixture.problem.apply_linear_transpose(cotangent, transpose_after,
                                           fixture.workspace);
    require_same(transpose_after, transpose_before);

    ProblemParameterData<1> gradient_after;
    Eigen::VectorXd adjoint_after;
    REQUIRE(fixture.problem
                .solve_vjp(fixture.forcing, state, cotangent, gradient_after,
                           settings, fixture.workspace, adjoint_after)
                .iterations == 2);
    require_same(adjoint_after, adjoint_before);
    require_same(gradient_after.forcing, gradient_before.forcing);
    require_same(gradient_after.transport_values,
                 gradient_before.transport_values);
    require_same(gradient_after.atmospheric_coefficients,
                 gradient_before.atmospheric_coefficients);
    require_same(gradient_after.ground_values, gradient_before.ground_values);
}

TEST_CASE(
    "Successive-orders scalar VJP can omit discarded scattering parameters",
    "[successive_orders][problem][vjp]") {
    const auto require_same = [](const auto& actual, const auto& expected) {
        REQUIRE(actual.rows() == expected.rows());
        REQUIRE(actual.cols() == expected.cols());
        if (actual.size() != 0) {
            REQUIRE(std::memcmp(actual.data(), expected.data(),
                                actual.size() * sizeof(double)) == 0);
        }
    };
    for (const int active : {1, 3}) {
        CAPTURE(active);
        ScalarProblemFixture fixture;
        Eigen::MatrixXd coefficients = Eigen::MatrixXd::Zero(1, 3);
        coefficients(0, 0) = 0.22;
        if (active == 3) {
            coefficients(0, 1) = 0.035;
            coefficients(0, 2) = -0.012;
        }
        fixture.scattering.set_atmospheric_coefficients(coefficients);
        Eigen::VectorXd state =
            Eigen::VectorXd::Zero(fixture.problem.state_size());
        REQUIRE(fixture.problem
                    .solve(fixture.forcing, state, tight_settings(),
                           fixture.workspace)
                    .converged());
        const Eigen::VectorXd cotangent = Eigen::VectorXd::LinSpaced(
            fixture.problem.state_size(), -0.4, 0.65);
        for (const int iterations : {0, 2, 150}) {
            CAPTURE(iterations);
            auto settings = tight_settings();
            settings.maximum_iterations = iterations;
            if (iterations != 150) {
                settings.relative_tolerance = 0.0;
                settings.absolute_tolerance = 0.0;
            }
            for (const bool materialize_transport : {false, true}) {
                CAPTURE(materialize_transport);
                ProblemWorkspace<1> full_workspace, input_workspace;
                full_workspace.resize(fixture.transport, fixture.scattering);
                input_workspace.resize(fixture.transport, fixture.scattering);
                ProblemParameterData<1> full_gradient, input_gradient;
                Eigen::VectorXd full_adjoint, input_adjoint;
                const auto full_diagnostics = fixture.problem.solve_vjp(
                    fixture.forcing, state, cotangent, full_gradient, settings,
                    full_workspace, full_adjoint, materialize_transport);
                // Prefill the output to catch stale parameter gradients when
                // a generic caller alternates full and input-only products.
                input_gradient = full_gradient;
                const auto input_diagnostics = fixture.problem.solve_vjp(
                    fixture.forcing, state, cotangent, input_gradient, settings,
                    input_workspace, input_adjoint, materialize_transport,
                    true);
                REQUIRE(input_diagnostics.termination ==
                        full_diagnostics.termination);
                REQUIRE(input_diagnostics.iterations ==
                        full_diagnostics.iterations);
                REQUIRE(input_diagnostics.residual_norm ==
                        full_diagnostics.residual_norm);
                require_same(input_adjoint, full_adjoint);
                require_same(input_gradient.forcing, full_gradient.forcing);
                require_same(input_gradient.transport_values,
                             full_gradient.transport_values);
                REQUIRE(input_gradient.atmospheric_coefficients.rows() == 1);
                REQUIRE(input_gradient.atmospheric_coefficients.cols() == 3);
                REQUIRE(input_gradient.atmospheric_coefficients.isZero(0.0));
                REQUIRE(input_gradient.ground_values.size() ==
                        fixture.scattering.ground_value_size());
                REQUIRE(input_gradient.ground_values.isZero(0.0));
                REQUIRE_FALSE(
                    full_gradient.atmospheric_coefficients.isZero(0.0));
                REQUIRE_FALSE(full_gradient.ground_values.isZero(0.0));
                // The default remains the complete generic parameter VJP.
                fixture.problem.solve_vjp(
                    fixture.forcing, state, cotangent, input_gradient, settings,
                    input_workspace, input_adjoint, materialize_transport);
                require_same(input_gradient.forcing, full_gradient.forcing);
                require_same(input_gradient.transport_values,
                             full_gradient.transport_values);
                require_same(input_gradient.atmospheric_coefficients,
                             full_gradient.atmospheric_coefficients);
                require_same(input_gradient.ground_values,
                             full_gradient.ground_values);
            }
        }
    }
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

    ProblemParameterData<3> unchanged_gradient;
    Eigen::VectorXd unchanged_adjoint;
    REQUIRE(problem
                .solve_vjp(forcing, state, state_cotangent, unchanged_gradient,
                           tight_settings(), workspace, unchanged_adjoint, true,
                           true)
                .converged());
    REQUIRE(std::memcmp(unchanged_gradient.forcing.data(),
                        gradient.forcing.data(),
                        gradient.forcing.size() * sizeof(double)) == 0);
    REQUIRE(std::memcmp(unchanged_gradient.atmospheric_coefficients.data(),
                        gradient.atmospheric_coefficients.data(),
                        gradient.atmospheric_coefficients.size() *
                            sizeof(double)) == 0);
}

#ifdef SKTRAN_RUST_SUPPORT
TEST_CASE("Successive-orders scalar input-only VJP preserves mapped native "
          "gradients after atmosphere updates",
          "[successive_orders][problem][engine][vjp]") {
    constexpr int altitudes_count = 9;
    constexpr int horizontal_count = 2;
    constexpr int locations_count = altitudes_count * horizontal_count;
    constexpr int wavelengths = 2;
    constexpr int rays = 3;
    sasktran2::Config config;
    config.set_num_threads(1);
    config.set_single_scatter_source(
        sasktran2::Config::SingleScatterSource::exact);
    config.set_multiple_scatter_source(
        sasktran2::Config::MultipleScatterSource::successive_orders);
    config.set_num_hr_incoming(14);
    config.set_num_hr_outgoing(14);
    config.set_num_hr_spherical_iterations(30);
    config.set_successive_orders_relative_tolerance(1.0e-8);
    config.set_successive_orders_absolute_tolerance(1.0e-12);
    config.set_num_do_streams(8);
    config.set_num_do_sza(horizontal_count);
    config.set_num_singlescatter_moments(8);
    config.set_apply_delta_scaling(false);
    Eigen::VectorXd altitudes =
        Eigen::VectorXd::LinSpaced(altitudes_count, 0.0, 40000.0);
    Eigen::VectorXd horizontal_angles(2);
    horizontal_angles << -0.4, 0.4;
    sasktran2::Geometry2D geometry(0.55, 0.15, 6372000.0, std::move(altitudes),
                                   std::move(horizontal_angles),
                                   sasktran2::grids::interpolation::linear);
    sasktran2::viewinggeometry::ViewingGeometryContainer viewing;
    viewing.observer_rays().emplace_back(
        std::make_unique<sasktran2::viewinggeometry::GroundViewingSolar>(
            0.55, 0.35, 0.72, 100000.0));
    viewing.observer_rays().emplace_back(
        std::make_unique<sasktran2::viewinggeometry::TangentAltitudeSolar>(
            10000.0, -0.4, 100000.0, 0.55));
    viewing.observer_rays().emplace_back(
        std::make_unique<sasktran2::viewinggeometry::TangentAltitudeSolar>(
            22500.0, 0.6, 100000.0, 0.55));
    const auto require_same = [](const auto& actual, const auto& expected) {
        REQUIRE(actual.size() == expected.size());
        REQUIRE(actual.allFinite());
        REQUIRE(std::memcmp(actual.data(), expected.data(),
                            actual.size() * sizeof(double)) == 0);
    };
    for (const bool spatial_surface : {false, true}) {
        CAPTURE(spatial_surface);
        sasktran2::atmosphere::Atmosphere<1> no_phase(wavelengths, geometry,
                                                      config, true);
        sasktran2::atmosphere::Atmosphere<1> zero_phase(wavelengths, geometry,
                                                        config, true);
        for (auto* atmosphere : {&no_phase, &zero_phase}) {
            atmosphere->storage().resize_derivatives(
                atmosphere == &no_phase ? 0 : 1);
            for (int wavelength = 0; wavelength < wavelengths; ++wavelength) {
                for (int horizontal = 0; horizontal < horizontal_count;
                     ++horizontal) {
                    for (int altitude = 0; altitude < altitudes_count;
                         ++altitude) {
                        const int location =
                            geometry.location_index(altitude, horizontal);
                        atmosphere->storage().total_extinction(location,
                                                               wavelength) =
                            (1.4e-5 * std::exp(-altitude / 1.7) + 2.0e-9) *
                            (0.9 + 0.2 * wavelength) *
                            (1.0 + 0.07 * horizontal);
                        atmosphere->storage().ssa(location, wavelength) =
                            0.90 - 0.008 * wavelength - 0.002 * altitude -
                            0.005 * horizontal;
                    }
                }
            }
            atmosphere->storage().leg_coeff.chip(0, 0).setConstant(1.0);
            atmosphere->storage().leg_coeff.chip(1, 0).setConstant(0.08);
            atmosphere->storage().leg_coeff.chip(2, 0).setConstant(0.5);
            atmosphere->surface().brdf_args().row(0).setConstant(0.12);
            if (spatial_surface) {
                atmosphere->surface().set_spatial_lambertian_albedo(
                    Eigen::MatrixXd::Constant(2, wavelengths, 0.12));
            }
            auto& mapping =
                atmosphere->storage().get_derivative_mapping("retrieval_state");
            mapping.allocate_extinction_derivatives();
            mapping.allocate_ssa_derivatives();
            mapping.native_mapping().d_extinction->setConstant(1.0e-5);
            mapping.native_mapping().d_ssa->setConstant(1.0e-2);
            if (atmosphere == &zero_phase) {
                mapping.allocate_legendre_derivatives();
                mapping.native_mapping().d_legendre->setZero();
                mapping.native_mapping().scat_factor->setOnes();
                atmosphere->storage().finalize_scattering_derivatives(0);
            }
            atmosphere->mark_changed();
        }
        REQUIRE(no_phase.num_scattering_deriv_groups() == 0);
        REQUIRE(zero_phase.num_scattering_deriv_groups() == 1);
        REQUIRE(no_phase.surface().has_spatial_lambertian_albedo() ==
                spatial_surface);
        const Eigen::MatrixXd initial_extinction =
            no_phase.storage().total_extinction;
        const Eigen::MatrixXd initial_ssa = no_phase.storage().ssa;
        Sasktran2<1> input_engine(config, &geometry, viewing);
        Sasktran2<1> full_engine(config, &geometry, viewing);
        const auto calculate = [&](Sasktran2<1>& engine,
                                   const auto& atmosphere) {
            Eigen::VectorXd radiance =
                Eigen::VectorXd::Zero(wavelengths * rays);
            const Eigen::VectorXd cotangent =
                Eigen::VectorXd::LinSpaced(wavelengths * rays, 0.4, 1.1);
            Eigen::VectorXd gradient = Eigen::VectorXd::Zero(locations_count);
            Eigen::Map<Eigen::VectorXd> radiance_map(radiance.data(),
                                                     radiance.size());
            Eigen::Map<const Eigen::VectorXd> cotangent_map(cotangent.data(),
                                                            cotangent.size());
            Eigen::Map<Eigen::VectorXd> gradient_map(gradient.data(),
                                                     gradient.size());
            sasktran2::OutputVJP<1> output(radiance_map, cotangent_map);
            output.set_derivative_gradient_memory("retrieval_state",
                                                  gradient_map);
            engine.calculate_vjp(atmosphere, output);
            output.finalize();
            return std::make_pair(std::move(radiance), std::move(gradient));
        };
        for (int evaluation = 0; evaluation < 4; ++evaluation) {
            CAPTURE(evaluation);
            for (auto* atmosphere : {&no_phase, &zero_phase}) {
                if (evaluation == 1) {
                    atmosphere->surface().brdf_args().row(0).setConstant(0.23);
                    if (spatial_surface) {
                        atmosphere->surface().set_spatial_lambertian_albedo(
                            Eigen::MatrixXd::Constant(2, wavelengths, 0.23));
                    }
                    atmosphere->mark_surface_changed();
                } else if (evaluation == 2) {
                    atmosphere->storage().total_extinction *= 1.25;
                    atmosphere->storage().ssa.array() -= 0.04;
                    atmosphere->mark_changed();
                } else if (evaluation == 3) {
                    atmosphere->storage().total_extinction = initial_extinction;
                    atmosphere->storage().ssa = initial_ssa;
                    atmosphere->surface().brdf_args().row(0).setConstant(0.12);
                    if (spatial_surface) {
                        atmosphere->surface().set_spatial_lambertian_albedo(
                            Eigen::MatrixXd::Constant(2, wavelengths, 0.12));
                    }
                    atmosphere->mark_changed();
                }
            }
            const auto input_only = calculate(input_engine, no_phase);
            const auto full = calculate(full_engine, zero_phase);
            require_same(input_only.first, full.first);
            require_same(input_only.second, full.second);
        }
    }
}
#endif
