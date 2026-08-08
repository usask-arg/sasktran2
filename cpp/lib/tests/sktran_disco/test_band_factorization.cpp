#include <sasktran2/test_helper.h>

#include <sktran_disco/sktran_do.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {
    struct BandSystem {
        lapack_int size;
        lapack_int bandwidth;
        lapack_int right_hand_sides;
        lapack_int leading_dimension;
        std::vector<double> matrix;
        std::vector<double> right_hand_side;
        std::vector<double> expected;
    };

    BandSystem make_band_system(lapack_int size, lapack_int bandwidth,
                                lapack_int right_hand_sides,
                                bool require_pivoting = false) {
        BandSystem result{
            size,
            bandwidth,
            right_hand_sides,
            3 * bandwidth + 1,
            std::vector<double>((3 * bandwidth + 1) * size, -71.25),
            std::vector<double>(size * right_hand_sides, 0.0),
            std::vector<double>(size * right_hand_sides, 0.0)};
        const lapack_int diagonal_row = 2 * bandwidth;
        for (lapack_int rhs = 0; rhs < right_hand_sides; ++rhs) {
            for (lapack_int row = 0; row < size; ++row) {
                result.expected[rhs * size + row] =
                    std::sin(0.13 * static_cast<double>((rhs + 1) * (row + 1)));
            }
        }
        for (lapack_int column = 0; column < size; ++column) {
            const lapack_int first_row =
                std::max<lapack_int>(0, column - bandwidth);
            const lapack_int last_row =
                std::min<lapack_int>(size - 1, column + bandwidth);
            for (lapack_int row = first_row; row <= last_row; ++row) {
                double value =
                    0.015 *
                    std::sin(static_cast<double>((row + 3) * (column + 5)));
                if (row == column) {
                    value += 2.5 + 0.001 * static_cast<double>(row);
                }
                if (require_pivoting && column % 5 == 0 && column + 1 < size) {
                    if (row == column) {
                        value = 0.25;
                    } else if (row == column + 1) {
                        value = 3.0;
                    }
                }
                const lapack_int band_row = diagonal_row + row - column;
                result.matrix[column * result.leading_dimension + band_row] =
                    value;
                for (lapack_int rhs = 0; rhs < right_hand_sides; ++rhs) {
                    result.right_hand_side[rhs * size + row] +=
                        value * result.expected[rhs * size + column];
                }
            }
        }
        return result;
    }

    int solve_with_configured_lapack(BandSystem& system,
                                     std::vector<lapack_int>& pivots) {
#if defined(SKTRAN_USE_ACCELERATE) || defined(SKTRAN_NO_LAPACKE)
        lapack_int info;
        dgbsv_(&system.size, &system.bandwidth, &system.bandwidth,
               &system.right_hand_sides, system.matrix.data(),
               &system.leading_dimension, pivots.data(),
               system.right_hand_side.data(), &system.size, &info);
        return info;
#else
        return LAPACKE_dgbsv(LAPACK_COL_MAJOR, system.size, system.bandwidth,
                             system.bandwidth, system.right_hand_sides,
                             system.matrix.data(), system.leading_dimension,
                             pivots.data(), system.right_hand_side.data(),
                             system.size);
#endif
    }

    void benchmark_band_solve(lapack_int size, lapack_int bandwidth) {
        const auto source = make_band_system(size, bandwidth, 1, true);

        BENCHMARK_ADVANCED("configured LAPACK")
        (Catch::Benchmark::Chronometer meter) {
            std::vector<BandSystem> systems(meter.runs(), source);
            std::vector<std::vector<lapack_int>> pivots(
                meter.runs(), std::vector<lapack_int>(size));
            meter.measure([&](int run) {
                return solve_with_configured_lapack(systems[run], pivots[run]);
            });
        };

        BENCHMARK_ADVANCED("unblocked C++")
        (Catch::Benchmark::Chronometer meter) {
            std::vector<BandSystem> systems(meter.runs(), source);
            std::vector<std::vector<lapack_int>> pivots(
                meter.runs(), std::vector<lapack_int>(size));
            meter.measure([&](int run) {
                auto& system = systems[run];
                return sasktran_disco::la::dgbsv_unblocked(
                    system.size, system.bandwidth, system.bandwidth,
                    system.right_hand_sides, system.matrix.data(),
                    system.leading_dimension, pivots[run].data(),
                    system.right_hand_side.data(), system.size);
            });
        };
    }
} // namespace

TEST_CASE("Unblocked band solve reproduces known solutions",
          "[sktran_do][band_factorization]") {
    for (const lapack_int bandwidth :
         {lapack_int{5}, lapack_int{11}, lapack_int{17}}) {
        auto system = make_band_system(6 * bandwidth + 7, bandwidth, 3);
        std::vector<lapack_int> pivots(system.size);
        const int info = sasktran_disco::la::dgbsv_unblocked(
            system.size, system.bandwidth, system.bandwidth,
            system.right_hand_sides, system.matrix.data(),
            system.leading_dimension, pivots.data(),
            system.right_hand_side.data(), system.size);
        REQUIRE(info == 0);
        for (std::size_t index = 0; index < system.expected.size(); ++index) {
            REQUIRE(std::abs(system.right_hand_side[index] -
                             system.expected[index]) < 2.0e-13);
        }
    }
}

TEST_CASE("Unblocked band factorization matches configured LAPACK",
          "[sktran_do][band_factorization]") {
    for (const lapack_int bandwidth :
         {lapack_int{5}, lapack_int{8}, lapack_int{11}, lapack_int{17},
          lapack_int{23}, lapack_int{26}, lapack_int{29}}) {
        auto custom = make_band_system(10 * bandwidth + 3, bandwidth, 3, true);
        auto configured_lapack = custom;
        std::vector<lapack_int> custom_pivots(custom.size);
        std::vector<lapack_int> configured_lapack_pivots(custom.size);

        const int custom_info = sasktran_disco::la::dgbsv_unblocked(
            custom.size, custom.bandwidth, custom.bandwidth,
            custom.right_hand_sides, custom.matrix.data(),
            custom.leading_dimension, custom_pivots.data(),
            custom.right_hand_side.data(), custom.size);
        const int configured_lapack_info = solve_with_configured_lapack(
            configured_lapack, configured_lapack_pivots);

        REQUIRE(custom_info == configured_lapack_info);
#if defined(SKTRAN_USE_ACCELERATE) || defined(SKTRAN_USE_MKL)
        for (std::size_t index = 0; index < custom.right_hand_side.size();
             ++index) {
            REQUIRE(custom.right_hand_side[index] ==
                    Catch::Approx(configured_lapack.right_hand_side[index])
                        .epsilon(2.0e-13)
                        .margin(2.0e-13));
        }
#else
        REQUIRE(custom_pivots == configured_lapack_pivots);
        REQUIRE(custom.matrix == configured_lapack.matrix);
        REQUIRE(custom.right_hand_side == configured_lapack.right_hand_side);
#endif
    }
}

TEST_CASE("Unblocked band factorization reports singular pivots",
          "[sktran_do][band_factorization]") {
    constexpr lapack_int size = 12;
    constexpr lapack_int bandwidth = 3;
    constexpr lapack_int leading_dimension = 3 * bandwidth + 1;
    std::vector<double> matrix(leading_dimension * size, 0.0);
    std::vector<lapack_int> pivots(size);
    REQUIRE(sasktran_disco::la::dgbtf2_unblocked(
                size, bandwidth, bandwidth, matrix.data(), leading_dimension,
                pivots.data()) == 1);
}

TEST_CASE("Unblocked band factorization dispatches for narrow DO systems",
          "[sktran_do][band_factorization]") {
    REQUIRE_FALSE(sasktran_disco::la::use_unblocked_band_factorization(2, 2));
#if defined(SKTRAN_USE_ACCELERATE) || defined(SKTRAN_USE_MKL)
    REQUIRE_FALSE(sasktran_disco::la::use_unblocked_band_factorization(5, 5));
    REQUIRE_FALSE(sasktran_disco::la::use_unblocked_band_factorization(29, 29));
#else
    REQUIRE(sasktran_disco::la::use_unblocked_band_factorization(5, 5));
    REQUIRE(sasktran_disco::la::use_unblocked_band_factorization(29, 29));
#endif
    REQUIRE_FALSE(sasktran_disco::la::use_unblocked_band_factorization(32, 32));
    REQUIRE_FALSE(sasktran_disco::la::use_unblocked_band_factorization(11, 12));
}

TEST_CASE("Band solve performance", "[.benchmark][band_solve]") {
    SECTION("scalar 4 streams, 20 layers") { benchmark_band_solve(4 * 20, 5); }
    SECTION("scalar 4 streams, 100 layers") {
        benchmark_band_solve(4 * 100, 5);
    }
    SECTION("scalar 8 streams, 20 layers") { benchmark_band_solve(8 * 20, 11); }
    SECTION("scalar 8 streams, 100 layers") {
        benchmark_band_solve(8 * 100, 11);
    }
    SECTION("scalar 16 streams, 20 layers") {
        benchmark_band_solve(16 * 20, 23);
    }
    SECTION("scalar 16 streams, 100 layers") {
        benchmark_band_solve(16 * 100, 23);
    }
    SECTION("IQU 2 streams, 20 layers") { benchmark_band_solve(6 * 20, 8); }
    SECTION("IQU 2 streams, 100 layers") { benchmark_band_solve(6 * 100, 8); }
    SECTION("IQU 4 streams, 20 layers") { benchmark_band_solve(12 * 20, 17); }
    SECTION("IQU 4 streams, 100 layers") { benchmark_band_solve(12 * 100, 17); }
    SECTION("IQU 8 streams, 20 layers") { benchmark_band_solve(24 * 20, 35); }
    SECTION("IQU 8 streams, 100 layers") { benchmark_band_solve(24 * 100, 35); }
    SECTION("bandwidth 26, 100 layers") { benchmark_band_solve(18 * 100, 26); }
    SECTION("bandwidth 29, 100 layers") { benchmark_band_solve(20 * 100, 29); }
    SECTION("bandwidth 32, 100 layers") { benchmark_band_solve(22 * 100, 32); }
}
