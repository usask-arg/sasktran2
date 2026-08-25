#include "../../successive_orders/fixed_point.h"

#include <sasktran2/test_helper.h>

#include <Eigen/LU>

#include <utility>
#include <vector>

namespace {
    sasktran2::successive_orders::FixedPointSettings test_settings() {
        sasktran2::successive_orders::FixedPointSettings settings;
        settings.maximum_iterations = 100;
        settings.relative_tolerance = 1.0e-12;
        settings.absolute_tolerance = 1.0e-14;
        return settings;
    }
} // namespace

TEST_CASE("Successive-orders fixed point converges to affine solution",
          "[successive_orders][solver]") {
    Eigen::Matrix2d transport;
    transport << 0.7, 0.05, 0.1, 0.4;
    const Eigen::Vector2d forcing(1.0, 0.25);
    const Eigen::Vector2d expected =
        (Eigen::Matrix2d::Identity() - transport).inverse() * forcing;

    Eigen::VectorXd state = Eigen::Vector2d::Zero();
    sasktran2::successive_orders::FixedPointWorkspace workspace;
    const auto diagnostics =
        sasktran2::successive_orders::FixedPointSolver::solve(
            state,
            [&](const Eigen::VectorXd& input, Eigen::VectorXd& output) {
                output.noalias() = forcing + transport * input;
            },
            test_settings(), workspace);

    REQUIRE(diagnostics.converged());
    REQUIRE(state.isApprox(expected, 1.0e-10));
}

TEST_CASE("Successive-orders Anderson history is reusable",
          "[successive_orders][solver]") {
    sasktran2::successive_orders::FixedPointWorkspace workspace;
    auto settings = test_settings();
    settings.anderson_depth = 3;

    for (double contraction : {0.85, 0.65}) {
        Eigen::VectorXd state = Eigen::VectorXd::Zero(8);
        Eigen::VectorXd forcing = Eigen::VectorXd::LinSpaced(8, 0.1, 0.8);
        const auto diagnostics =
            sasktran2::successive_orders::FixedPointSolver::solve(
                state,
                [&](const Eigen::VectorXd& input, Eigen::VectorXd& output) {
                    output = forcing + contraction * input;
                },
                settings, workspace);
        REQUIRE(diagnostics.converged());
        REQUIRE(state.isApprox(forcing / (1.0 - contraction), 1.0e-9));
    }
}

TEST_CASE("Successive-orders zero tolerances preserve fixed iteration mode",
          "[successive_orders][solver]") {
    auto settings = test_settings();
    settings.maximum_iterations = 4;
    settings.relative_tolerance = 0.0;
    settings.absolute_tolerance = 0.0;
    settings.anderson_depth = 0;

    Eigen::VectorXd state = Eigen::VectorXd::Zero(1);
    sasktran2::successive_orders::FixedPointWorkspace workspace;
    const auto diagnostics =
        sasktran2::successive_orders::FixedPointSolver::solve(
            state,
            [](const Eigen::VectorXd& input, Eigen::VectorXd& output) {
                output(0) = 1.0 + 0.5 * input(0);
            },
            settings, workspace);

    REQUIRE_FALSE(diagnostics.converged());
    REQUIRE(diagnostics.iterations == 4);
    REQUIRE(state(0) == Catch::Approx(1.875));
}

TEST_CASE("Successive-orders diagnostics separate material residuals from "
          "solver tolerance",
          "[successive_orders][solver]") {
    auto settings = test_settings();
    settings.maximum_iterations = 1;
    settings.relative_tolerance = 1.0e-12;
    settings.absolute_tolerance = 0.0;
    settings.anderson_depth = 0;

    sasktran2::successive_orders::FixedPointWorkspace workspace;
    for (const auto& test_case : std::vector<std::pair<double, bool>>{
             {5.0e-4, false}, {2.0e-3, true}}) {
        const double increment = test_case.first;
        const bool materially_unconverged = test_case.second;
        Eigen::VectorXd state = Eigen::VectorXd::Ones(1);
        const auto diagnostics =
            sasktran2::successive_orders::FixedPointSolver::solve(
                state,
                [increment](const Eigen::VectorXd&, Eigen::VectorXd& output) {
                    output(0) = 1.0 + increment;
                },
                settings, workspace);

        REQUIRE_FALSE(diagnostics.converged());
        REQUIRE(diagnostics.state_scale == Catch::Approx(1.0 + increment));
        REQUIRE(diagnostics.warrants_convergence_warning() ==
                materially_unconverged);
    }
}

TEST_CASE("Successive-orders damped Anderson uses undamped state history",
          "[successive_orders][solver]") {
    auto settings = test_settings();
    settings.maximum_iterations = 2;
    settings.relative_tolerance = 0.0;
    settings.absolute_tolerance = 0.0;
    settings.damping = 0.5;
    settings.anderson_depth = 1;

    Eigen::VectorXd state = Eigen::VectorXd::Zero(1);
    sasktran2::successive_orders::FixedPointWorkspace workspace;
    const auto diagnostics =
        sasktran2::successive_orders::FixedPointSolver::solve(
            state,
            [](const Eigen::VectorXd& input, Eigen::VectorXd& output) {
                output(0) = 1.0 + 0.5 * input(0);
            },
            settings, workspace);

    // x_1 = 0.5, delta_x = 0.5, delta_f = -0.25, and gamma = -3.
    // Type-II damping gives
    // x_2 = x_1 + beta*f_1 - gamma*(delta_x + beta*delta_f) = 2.
    // Scaling delta_x by beta as well would instead give 1.25.
    REQUIRE_FALSE(diagnostics.converged());
    REQUIRE(diagnostics.iterations == 2);
    REQUIRE(state(0) == Catch::Approx(2.0).margin(1.0e-12));
}

TEST_CASE("Successive-orders implicit transpose uses the same solver",
          "[successive_orders][solver][vjp]") {
    Eigen::Matrix3d transport;
    transport << 0.35, 0.04, 0.01, 0.12, 0.45, 0.02, 0.03, 0.08, 0.25;
    const Eigen::Vector3d forcing(0.2, 0.5, -0.1);
    const Eigen::Vector3d expected =
        (Eigen::Matrix3d::Identity() - transport.transpose()).inverse() *
        forcing;

    Eigen::VectorXd adjoint = Eigen::Vector3d::Zero();
    sasktran2::successive_orders::FixedPointWorkspace workspace;
    const auto diagnostics =
        sasktran2::successive_orders::FixedPointSolver::solve(
            adjoint,
            [&](const Eigen::VectorXd& input, Eigen::VectorXd& output) {
                output.noalias() = forcing + transport.transpose() * input;
            },
            test_settings(), workspace);

    REQUIRE(diagnostics.converged());
    REQUIRE(adjoint.isApprox(expected, 1.0e-10));
}
