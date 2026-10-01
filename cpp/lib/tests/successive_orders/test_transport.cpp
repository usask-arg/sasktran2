#include "../../successive_orders/transport.h"

#include <sasktran2/test_helper.h>
#include <cstring>

namespace {
    sasktran2::successive_orders::TransportSparsity test_sparsity() {
        return {4, {0, 2, 5, 7}, {0, 2, 0, 1, 3, 1, 2}};
    }

    void require_same_bits(const Eigen::VectorXd& left,
                           const Eigen::VectorXd& right) {
        REQUIRE(left.size() == right.size());
        REQUIRE(std::memcmp(left.data(), right.data(),
                            left.size() * sizeof(double)) == 0);
    }

    template <int NSTOKES> void compare_column_width_products() {
        using sasktran2::successive_orders::TransportOperator;
        using sasktran2::successive_orders::TransportSparsity;
        constexpr int columns = 65536;
        const std::vector<int> offsets{0, 3, 5, 6};
        const std::vector<int> indices{0, 7, 65535, 7, 65535, 65535};
        const TransportSparsity compact(columns, offsets, indices);
        const TransportSparsity wide(columns + 1, offsets, indices);
        REQUIRE(compact.compact_column_indices());
        REQUIRE_FALSE(wide.compact_column_indices());
        REQUIRE(wide.storage_bytes() - compact.storage_bytes() ==
                indices.size() * sizeof(std::uint16_t));
        TransportOperator compact_operator(compact);
        TransportOperator wide_operator(wide);
        compact_operator.values() << 0.23, -0.41, 0.19, 0.0, 0.36, -0.0;
        wide_operator.values() = compact_operator.values();
        const Eigen::VectorXd state =
            Eigen::VectorXd::LinSpaced(columns * NSTOKES, -0.2, 0.7);
        const Eigen::VectorXd tangent =
            Eigen::VectorXd::LinSpaced(columns * NSTOKES, 0.3, -0.5);
        Eigen::VectorXd wide_state((columns + 1) * NSTOKES);
        Eigen::VectorXd wide_tangent((columns + 1) * NSTOKES);
        wide_state.head(state.size()) = state;
        wide_tangent.head(tangent.size()) = tangent;
        wide_state.tail(NSTOKES).setConstant(19.0);
        wide_tangent.tail(NSTOKES).setConstant(-13.0);
        const Eigen::VectorXd cotangent =
            Eigen::VectorXd::LinSpaced(3 * NSTOKES, -0.6, 0.4);
        const Eigen::VectorXd value_tangent =
            Eigen::VectorXd::LinSpaced(indices.size(), 0.05, -0.02);
        Eigen::VectorXd compact_result(3 * NSTOKES);
        Eigen::VectorXd wide_result(3 * NSTOKES);
        Eigen::VectorXd compact_reverse(columns * NSTOKES);
        Eigen::VectorXd wide_reverse((columns + 1) * NSTOKES);
        Eigen::VectorXd compact_gradient(indices.size());
        Eigen::VectorXd wide_gradient(indices.size());

        compact_operator.apply_stokes<NSTOKES>(state, compact_result);
        wide_operator.apply_stokes<NSTOKES>(wide_state, wide_result);
        require_same_bits(compact_result, wide_result);
        compact_operator.apply_transpose_stokes<NSTOKES>(cotangent,
                                                         compact_reverse);
        wide_operator.apply_transpose_stokes<NSTOKES>(cotangent, wide_reverse);
        require_same_bits(compact_reverse,
                          wide_reverse.head(compact_reverse.size()));
        REQUIRE(wide_reverse.tail(NSTOKES).isZero(0.0));
        compact_operator.apply_jvp_stokes<NSTOKES>(
            state, tangent, value_tangent, compact_result);
        wide_operator.apply_jvp_stokes<NSTOKES>(wide_state, wide_tangent,
                                                value_tangent, wide_result);
        require_same_bits(compact_result, wide_result);
        compact_operator.apply_value_jvp_stokes<NSTOKES>(state, value_tangent,
                                                         compact_result);
        wide_operator.apply_value_jvp_stokes<NSTOKES>(wide_state, value_tangent,
                                                      wide_result);
        require_same_bits(compact_result, wide_result);
        compact_operator.apply_vjp_stokes<NSTOKES>(
            state, cotangent, compact_reverse, compact_gradient);
        wide_operator.apply_vjp_stokes<NSTOKES>(wide_state, cotangent,
                                                wide_reverse, wide_gradient);
        require_same_bits(compact_reverse,
                          wide_reverse.head(compact_reverse.size()));
        require_same_bits(compact_gradient, wide_gradient);
    }
} // namespace

TEST_CASE("Successive-orders CSR compact indices validate the 16-bit boundary",
          "[successive_orders][transport][compact_indices]") {
    using sasktran2::successive_orders::TransportSparsity;
    const TransportSparsity compact(65536, {0, 2}, {0, 65535});
    REQUIRE(compact.column_indices().element_bytes() == 2);
    REQUIRE(compact.column_indices()[1] == 65535);
    const TransportSparsity wide(65537, {0, 2}, {0, 65536});
    REQUIRE(wide.column_indices().element_bytes() == sizeof(int));
    REQUIRE(wide.column_indices()[1] == 65536);
    REQUIRE_THROWS_AS((TransportSparsity(65536, {0, 1}, {65536})),
                      std::invalid_argument);
    REQUIRE_THROWS_AS((TransportSparsity(65536, {0, 1}, {-1})),
                      std::invalid_argument);
    REQUIRE_THROWS_AS((TransportSparsity(65537, {0, 1}, {65537})),
                      std::invalid_argument);
    REQUIRE_THROWS_AS((TransportSparsity(65536, {0, 2}, {7, 7})),
                      std::invalid_argument);
    REQUIRE_THROWS_AS((TransportSparsity(65536, {0, 2}, {7, 0})),
                      std::invalid_argument);
    REQUIRE_THROWS_AS((TransportSparsity(65536, {0, -1}, {})),
                      std::invalid_argument);
    REQUIRE_THROWS_AS((TransportSparsity(65536, {0, 3}, {0, 7})),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(compact.column_indices().subview(1, 2),
                      std::out_of_range);
    REQUIRE(TransportSparsity().column_indices().to_vector().empty());
}

TEST_CASE("Successive-orders compact CSR products preserve accumulation bits",
          "[successive_orders][transport][compact_indices]") {
    compare_column_width_products<1>();
    compare_column_width_products<3>();
}

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
