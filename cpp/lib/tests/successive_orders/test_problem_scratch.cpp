#include "../../successive_orders/problem_scratch.h"

#include <sasktran2.h>
#include <sasktran2/math/unitsphere.h>
#include <sasktran2/test_helper.h>

#include <array>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <future>
#include <limits>
#include <mutex>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
    using namespace sasktran2::successive_orders;

    TransportSparsity pooled_dense_sparsity(int size) {
        std::vector<int> offsets(static_cast<std::size_t>(size) + 1);
        std::vector<int> columns;
        for (int row = 0; row < size; ++row) {
            offsets[row] = static_cast<int>(columns.size());
            for (int column = 0; column < size; ++column) {
                columns.push_back(column);
            }
        }
        offsets[size] = static_cast<int>(columns.size());
        return {size, std::move(offsets), std::move(columns)};
    }

    ScalarAngularBasis pooled_basis(int moments) {
        sasktran2::math::LebedevSphere incoming(6), outgoing(6);
        return {incoming, outgoing, moments};
    }

    struct PooledProblemFixture {
        PooledProblemFixture(int atmospheric_points, int moments)
            : sparsity(pooled_dense_sparsity(6 * atmospheric_points + 2)),
              transport(sparsity),
              scattering(
                  ScatteringBlockLayout(atmospheric_points, 1, 6, 6, 2, 2, 1),
                  pooled_basis(moments)),
              problem(transport, scattering),
              forcing(Eigen::VectorXd::LinSpaced(problem.incoming_size(), 0.04,
                                                 0.11)) {
            for (int index = 0; index < transport.values().size(); ++index) {
                transport.values()(index) = 0.001 * (1 + index % 5);
            }
            auto& coefficients = scattering.atmospheric_coefficients();
            for (int point = 0; point < coefficients.rows(); ++point) {
                for (int degree = 0; degree < coefficients.cols(); ++degree) {
                    coefficients(point, degree) = degree == 0
                                                      ? 0.22 + 0.01 * point
                                                      : 0.021 / (degree + 1);
                }
            }
            Eigen::Matrix2d ground;
            ground << 0.12, 0.025, -0.015, 0.09;
            scattering.set_ground_block(0, ground);
        }

        ProblemParameterData<1> tangent() const {
            ProblemParameterData<1> result;
            result.resize(transport, scattering);
            result.forcing = Eigen::VectorXd::LinSpaced(problem.incoming_size(),
                                                        -0.025, 0.018);
            for (int index = 0; index < result.transport_values.size();
                 ++index) {
                result.transport_values(index) = 0.0007 * ((index % 7) - 3);
            }
            for (int point = 0; point < result.atmospheric_coefficients.rows();
                 ++point) {
                for (int degree = 0;
                     degree < result.atmospheric_coefficients.cols();
                     ++degree) {
                    result.atmospheric_coefficients(point, degree) =
                        0.002 * (point - degree + 1);
                }
            }
            result.ground_values = Eigen::VectorXd::LinSpaced(
                result.ground_values.size(), -0.008, 0.012);
            return result;
        }

        TransportSparsity sparsity;
        TransportOperator transport;
        ScatteringOperator<1> scattering;
        Problem<1> problem;
        Eigen::VectorXd forcing;
    };

    template <typename Actual, typename Expected>
    void pooled_require_bits(const Actual& actual, const Expected& expected) {
        REQUIRE(actual.rows() == expected.rows());
        REQUIRE(actual.cols() == expected.cols());
        if (actual.size() != 0) {
            REQUIRE(actual.allFinite());
            REQUIRE(expected.allFinite());
            REQUIRE(std::memcmp(actual.data(), expected.data(),
                                actual.size() * sizeof(double)) == 0);
        }
    }

    void pooled_require_diagnostics(const FixedPointDiagnostics& actual,
                                    const FixedPointDiagnostics& expected) {
        REQUIRE(actual.termination == expected.termination);
        REQUIRE(actual.iterations == expected.iterations);
        REQUIRE(actual.residual_norm == expected.residual_norm);
        REQUIRE(actual.state_scale == expected.state_scale);
        REQUIRE(actual.convergence_threshold == expected.convergence_threshold);
    }

    FixedPointSettings pooled_settings(int maximum_iterations, int depth) {
        FixedPointSettings result;
        result.maximum_iterations = maximum_iterations;
        result.anderson_depth = depth;
        result.relative_tolerance = maximum_iterations == 150 ? 1.0e-13 : 0.0;
        result.absolute_tolerance = maximum_iterations == 150 ? 1.0e-14 : 0.0;
        return result;
    }
} // namespace

TEST_CASE("Scalar problem leases retain owning workspace between calls",
          "[successive_orders][problem_scratch]") {
    std::uint64_t id;
    const double* incoming;
    const double* mapped;
    std::size_t bytes;
    {
        ScalarProblemWorkspaceLease lease;
        REQUIRE(lease.shared());
        id = lease.allocation_id();
        auto& workspace = lease.workspace();
        workspace.incoming.setConstant(19, 0.125);
        workspace.auxiliary_incoming.setConstant(19, -0.375);
        workspace.direct_state.setConstant(23, 0.625);
        workspace.fixed_point.resize(23, 4);
        workspace.fixed_point.mapped().setConstant(0.875);
        incoming = workspace.incoming.data();
        mapped = workspace.fixed_point.mapped().data();
        bytes = lease.storage_bytes();
        REQUIRE(bytes == lease.problem_bytes() + lease.fixed_point_bytes());
    }
    ScalarProblemWorkspaceLease lease;
    REQUIRE(lease.shared());
    REQUIRE(lease.allocation_id() == id);
    REQUIRE(lease.workspace().incoming.data() == incoming);
    REQUIRE(lease.workspace().fixed_point.mapped().data() == mapped);
    REQUIRE(lease.workspace().incoming.isConstant(0.125));
    REQUIRE(lease.workspace().auxiliary_incoming.isConstant(-0.375));
    REQUIRE(lease.workspace().direct_state.isConstant(0.625));
    REQUIRE(lease.storage_bytes() == bytes);
}

TEST_CASE("Fixed-point workspace replacement preserves storage on failure",
          "[successive_orders][problem_scratch]") {
    const auto overflow_size = std::numeric_limits<Eigen::Index>::max();
    REQUIRE(static_cast<std::size_t>(overflow_size) >
            std::numeric_limits<std::size_t>::max() / sizeof(double));
    FixedPointWorkspace workspace;
    workspace.resize(7, 3);
    workspace.mapped().setConstant(0.25);
    workspace.residual().setConstant(-0.375);
    const auto* mapped = workspace.mapped().data();
    const auto* residual = workspace.residual().data();
    const auto bytes = workspace.storage_bytes();
    REQUIRE_THROWS_AS(workspace.resize(overflow_size, 0), std::bad_alloc);
    REQUIRE(workspace.storage_bytes() == bytes);
    REQUIRE(workspace.mapped().data() == mapped);
    REQUIRE(workspace.residual().data() == residual);
    REQUIRE(workspace.mapped().isConstant(0.25));
    REQUIRE(workspace.residual().isConstant(-0.375));
    workspace.resize(7, 3);
    REQUIRE(workspace.mapped().data() == mapped);
    REQUIRE(workspace.residual().data() == residual);

    // Repair one companion while the previously guarded mapped shape agrees.
    workspace.residual().resize(2);
    workspace.resize(7, 3);
    REQUIRE(workspace.mapped().size() == 7);
    REQUIRE(workspace.residual().size() == 7);
    REQUIRE(workspace.storage_bytes() == bytes);

    Eigen::VectorXd state = Eigen::VectorXd::Constant(7, 0.125);
    Eigen::VectorXd independent_state = state;
    FixedPointWorkspace independent;
    const auto map = [](const Eigen::VectorXd& input, Eigen::VectorXd& output) {
        output = 0.3 * input.array() + 0.2;
    };
    const auto actual =
        FixedPointSolver::solve(state, map, pooled_settings(150, 3), workspace);
    const auto expected = FixedPointSolver::solve(
        independent_state, map, pooled_settings(150, 3), independent);
    pooled_require_bits(state, independent_state);
    pooled_require_diagnostics(actual, expected);
}

TEST_CASE("Problem preparation repairs all companion vector dimensions",
          "[successive_orders][problem_scratch]") {
    PooledProblemFixture fixture(2, 3);
    ScalarProblemWorkspaceLease lease;
    auto& workspace = lease.workspace();
    workspace.resize(fixture.transport, fixture.scattering);
    workspace.auxiliary_incoming.resize(0);
    workspace.auxiliary_state.resize(0);
    const auto settings = pooled_settings(150, 3);
    Eigen::VectorXd state = Eigen::VectorXd::Zero(fixture.problem.state_size());
    Eigen::VectorXd independent_state = state;
    ProblemWorkspace<1> independent;
    const auto actual =
        fixture.problem.solve(fixture.forcing, state, settings, workspace);
    const auto expected = fixture.problem.solve(
        fixture.forcing, independent_state, settings, independent);
    REQUIRE(workspace.incoming.size() == fixture.problem.incoming_size());
    REQUIRE(workspace.auxiliary_incoming.size() ==
            fixture.problem.incoming_size());
    REQUIRE(workspace.auxiliary_state.size() == fixture.problem.state_size());
    REQUIRE(workspace.direct_state.size() == fixture.problem.state_size());
    pooled_require_bits(state, independent_state);
    pooled_require_diagnostics(actual, expected);
}

TEST_CASE("Scalar problem workspace reuse preserves complete native products",
          "[successive_orders][problem_scratch][linearization]") {
    // Equal input/output shapes with different phase orders exercise
    // scattering scratch changes that Problem::prepare_workspace cannot see.
    const std::array<std::pair<int, int>, 5> shapes{
        std::pair<int, int>{1, 3}, {2, 3}, {1, 2}, {1, 3}, {2, 2}};
    std::uint64_t shared_id = 0;
    for (const auto [points, moments] : shapes) {
        CAPTURE(points, moments);
        PooledProblemFixture fixture(points, moments);
        ProblemWorkspace<1> independent;
        const auto tangent = fixture.tangent();
        const auto original_transport = fixture.transport.values().eval();
        const auto original_coefficients =
            fixture.scattering.atmospheric_coefficients().eval();
        const auto original_ground = fixture.scattering.ground_values().eval();
        const auto original_forcing = fixture.forcing.eval();
        for (const int iterations : {150, 0, 2, 150}) {
            CAPTURE(iterations);
            for (const int depth : {3, 0, 5}) {
                CAPTURE(depth);
                const auto settings = pooled_settings(iterations, depth);
                Eigen::VectorXd independent_state = Eigen::VectorXd::LinSpaced(
                    fixture.problem.state_size(), -0.03, 0.02);
                Eigen::VectorXd pooled_state = independent_state;
                for (int update = 0; update < 3; ++update) {
                    CAPTURE(update);
                    const double scale = update == 1 ? 0.02 : 0.0;
                    fixture.transport.values() =
                        original_transport + scale * tangent.transport_values;
                    fixture.scattering.atmospheric_coefficients() =
                        original_coefficients +
                        scale * tangent.atmospheric_coefficients;
                    fixture.scattering.ground_values() =
                        original_ground + scale * tangent.ground_values;
                    fixture.forcing =
                        original_forcing + scale * tangent.forcing;
                    FixedPointDiagnostics pooled_diagnostics;
                    {
                        ScalarProblemWorkspaceLease lease;
                        REQUIRE(lease.shared());
                        if (shared_id == 0)
                            shared_id = lease.allocation_id();
                        REQUIRE(lease.allocation_id() == shared_id);
                        pooled_diagnostics =
                            fixture.problem.solve(fixture.forcing, pooled_state,
                                                  settings, lease.workspace());
                    }
                    const auto independent_diagnostics = fixture.problem.solve(
                        fixture.forcing, independent_state, settings,
                        independent);
                    pooled_require_diagnostics(pooled_diagnostics,
                                               independent_diagnostics);
                    pooled_require_bits(pooled_state, independent_state);
                    Eigen::VectorXd pooled_jvp, independent_jvp;
                    {
                        ScalarProblemWorkspaceLease lease;
                        pooled_diagnostics = fixture.problem.solve_jvp(
                            fixture.forcing, pooled_state, tangent, pooled_jvp,
                            settings, lease.workspace());
                    }
                    const auto independent_jvp_diagnostics =
                        fixture.problem.solve_jvp(
                            fixture.forcing, independent_state, tangent,
                            independent_jvp, settings, independent);
                    pooled_require_diagnostics(pooled_diagnostics,
                                               independent_jvp_diagnostics);
                    pooled_require_bits(pooled_jvp, independent_jvp);
                    const Eigen::VectorXd cotangent =
                        Eigen::VectorXd::LinSpaced(fixture.problem.state_size(),
                                                   -0.4, 0.65);
                    ProblemParameterData<1> pooled_gradient,
                        independent_gradient;
                    Eigen::VectorXd pooled_adjoint, independent_adjoint;
                    {
                        ScalarProblemWorkspaceLease lease;
                        pooled_diagnostics = fixture.problem.solve_vjp(
                            fixture.forcing, pooled_state, cotangent,
                            pooled_gradient, settings, lease.workspace(),
                            pooled_adjoint);
                    }
                    const auto independent_vjp_diagnostics =
                        fixture.problem.solve_vjp(
                            fixture.forcing, independent_state, cotangent,
                            independent_gradient, settings, independent,
                            independent_adjoint);
                    pooled_require_diagnostics(pooled_diagnostics,
                                               independent_vjp_diagnostics);
                    pooled_require_bits(pooled_adjoint, independent_adjoint);
                    pooled_require_bits(pooled_gradient.forcing,
                                        independent_gradient.forcing);
                    pooled_require_bits(pooled_gradient.transport_values,
                                        independent_gradient.transport_values);
                    pooled_require_bits(
                        pooled_gradient.atmospheric_coefficients,
                        independent_gradient.atmospheric_coefficients);
                    pooled_require_bits(pooled_gradient.ground_values,
                                        independent_gradient.ground_values);
                }
            }
        }
    }
}

TEST_CASE("A zero direct JVP cannot leak pooled history into the next solve",
          "[successive_orders][problem_scratch][linearization]") {
    PooledProblemFixture fixture(1, 3);
    Eigen::VectorXd state =
        Eigen::VectorXd::LinSpaced(fixture.problem.state_size(), -0.03, 0.02);
    Eigen::VectorXd independent_state = state;
    ProblemWorkspace<1> independent;
    std::uint64_t id;
    Eigen::VectorXd previous_residual;
    FixedPointDiagnostics pooled_diagnostics;
    {
        ScalarProblemWorkspaceLease lease;
        id = lease.allocation_id();
        // Leave nonempty Anderson history from a deliberately incomplete
        // solve, rather than priming the pool with an already fixed state.
        pooled_diagnostics = fixture.problem.solve(
            fixture.forcing, state, pooled_settings(2, 3), lease.workspace());
        previous_residual = lease.workspace().fixed_point.residual();
    }
    const auto independent_diagnostics = fixture.problem.solve(
        fixture.forcing, independent_state, pooled_settings(2, 3), independent);
    pooled_require_diagnostics(pooled_diagnostics, independent_diagnostics);
    pooled_require_bits(state, independent_state);
    auto zero_tangent = fixture.tangent();
    zero_tangent.set_zero();
    Eigen::VectorXd pooled_jvp =
        Eigen::VectorXd::Constant(fixture.problem.state_size(), 19.0);
    Eigen::VectorXd independent_jvp = pooled_jvp;
    {
        ScalarProblemWorkspaceLease lease;
        REQUIRE(lease.allocation_id() == id);
        pooled_diagnostics = fixture.problem.solve_jvp(
            fixture.forcing, state, zero_tangent, pooled_jvp,
            pooled_settings(150, 3), lease.workspace());
        // The shortcut does not enter FixedPointSolver. Its old residual is
        // still present, and the following solve must reset history itself.
        pooled_require_bits(lease.workspace().fixed_point.residual(),
                            previous_residual);
    }
    const auto independent_jvp_diagnostics = fixture.problem.solve_jvp(
        fixture.forcing, independent_state, zero_tangent, independent_jvp,
        pooled_settings(150, 3), independent);
    pooled_require_diagnostics(pooled_diagnostics, independent_jvp_diagnostics);
    REQUIRE(pooled_diagnostics.converged());
    REQUIRE(pooled_diagnostics.iterations == 0);
    REQUIRE(pooled_jvp.isZero(0.0));
    pooled_require_bits(pooled_jvp, independent_jvp);

    fixture.transport.values() *= 1.07;
    fixture.scattering.atmospheric_coefficients() *= 0.91;
    fixture.forcing *= 1.03;
    ProblemWorkspace<1> fresh_independent;
    // The reference has no prior solver history. Both retain the identical
    // physical warm seed while solving the changed problem to tolerance.
    independent_state = state;
    {
        ScalarProblemWorkspaceLease lease;
        REQUIRE(lease.allocation_id() == id);
        pooled_diagnostics = fixture.problem.solve(
            fixture.forcing, state, pooled_settings(150, 3), lease.workspace());
    }
    const auto fresh_diagnostics =
        fixture.problem.solve(fixture.forcing, independent_state,
                              pooled_settings(150, 3), fresh_independent);
    REQUIRE(pooled_diagnostics.converged());
    pooled_require_diagnostics(pooled_diagnostics, fresh_diagnostics);
    pooled_require_bits(state, independent_state);
}

TEST_CASE("Nested scalar problem solves cannot overwrite outer scratch",
          "[successive_orders][problem_scratch][reentrant]") {
    PooledProblemFixture inner(2, 2);
    ScalarProblemWorkspaceLease outer;
    auto& workspace = outer.workspace();
    workspace.incoming.setConstant(7, 0.125);
    workspace.direct_state.setConstant(11, -0.375);
    const double* incoming = workspace.incoming.data();
    const auto original_id = outer.allocation_id();
    bool invoked = false;
    Eigen::VectorXd state = Eigen::VectorXd::Constant(11, 0.4);
    FixedPointSolver::solve(
        state,
        [&](const Eigen::VectorXd& input, Eigen::VectorXd& output) {
            ScalarProblemWorkspaceLease nested;
            REQUIRE_FALSE(nested.shared());
            REQUIRE(nested.allocation_id() != original_id);
            Eigen::VectorXd nested_state =
                Eigen::VectorXd::Zero(inner.problem.state_size());
            inner.problem.solve(inner.forcing, nested_state,
                                pooled_settings(2, 5), nested.workspace());
            REQUIRE(workspace.incoming.data() == incoming);
            REQUIRE(workspace.incoming.isConstant(0.125));
            REQUIRE(workspace.direct_state.isConstant(-0.375));
            output = 0.25 * input;
            invoked = true;
        },
        pooled_settings(2, 3), workspace.fixed_point);
    REQUIRE(invoked);
    REQUIRE(workspace.incoming.data() == incoming);
    REQUIRE(workspace.incoming.isConstant(0.125));
}

TEST_CASE("Scalar problem workspace leases release after exceptions",
          "[successive_orders][problem_scratch]") {
    std::uint64_t id = 0;
    try {
        ScalarProblemWorkspaceLease lease;
        id = lease.allocation_id();
        Eigen::VectorXd state = Eigen::VectorXd::Ones(13);
        FixedPointSolver::solve(
            state,
            [](const Eigen::VectorXd&, Eigen::VectorXd&) {
                throw std::runtime_error("Solve interrupted");
            },
            pooled_settings(2, 4), lease.workspace().fixed_point);
        FAIL("The interrupted solve must throw");
    } catch (const std::runtime_error&) {
    }
    ScalarProblemWorkspaceLease lease;
    REQUIRE(lease.shared());
    REQUIRE(lease.allocation_id() == id);
    PooledProblemFixture fixture(1, 3);
    Eigen::VectorXd state = Eigen::VectorXd::Zero(fixture.problem.state_size());
    REQUIRE(fixture.problem
                .solve(fixture.forcing, state, pooled_settings(150, 3),
                       lease.workspace())
                .converged());
}

TEST_CASE("Concurrent OS threads own distinct scalar problem workspaces",
          "[successive_orders][problem_scratch]") {
    std::mutex mutex;
    std::condition_variable ready;
    int arrivals = 0;
    const auto solve = [&](int points) {
        PooledProblemFixture fixture(points, 3);
        ScalarProblemWorkspaceLease lease;
        Eigen::VectorXd state =
            Eigen::VectorXd::Zero(fixture.problem.state_size());
        const auto diagnostics = fixture.problem.solve(
            fixture.forcing, state, pooled_settings(150, 3), lease.workspace());
        const Eigen::VectorXd saved = lease.workspace().incoming;
        {
            std::unique_lock<std::mutex> lock(mutex);
            ++arrivals;
            ready.notify_all();
            ready.wait(lock, [&]() { return arrivals == 2; });
        }
        return std::make_pair(
            lease.allocation_id(),
            lease.shared() && diagnostics.converged() &&
                std::memcmp(saved.data(), lease.workspace().incoming.data(),
                            saved.size() * sizeof(double)) == 0);
    };
    auto first = std::async(std::launch::async, solve, 1);
    auto second = std::async(std::launch::async, solve, 2);
    const auto first_result = first.get();
    const auto second_result = second.get();
    REQUIRE(first_result.first != second_result.first);
    REQUIRE(first_result.second);
    REQUIRE(second_result.second);
}
