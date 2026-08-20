#include <sasktran2/test_helper.h>
#include <sasktran2/runtime_backend_tuner.h>

#include <sktran_disco/sktran_do.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
    bool set_environment_variable(
        const std::string& name,
        const std::optional<std::string>& value = std::nullopt) {
#if defined(_WIN32)
        return ::_putenv_s(name.c_str(), value ? value->c_str() : "") == 0;
#else
        return value ? ::setenv(name.c_str(), value->c_str(), 1) == 0
                     : ::unsetenv(name.c_str()) == 0;
#endif
    }

    class ScopedEnvironmentVariable {
      public:
        ScopedEnvironmentVariable(std::string name,
                                  std::optional<std::string> value)
            : m_name(std::move(name)) {
            if (const char* original = std::getenv(m_name.c_str());
                original != nullptr) {
                m_original = original;
            }
            if (!set_environment_variable(m_name, value)) {
                throw std::runtime_error("Failed to set environment variable");
            }
        }

        ~ScopedEnvironmentVariable() {
            (void)set_environment_variable(m_name, m_original);
        }

        ScopedEnvironmentVariable(const ScopedEnvironmentVariable&) = delete;
        ScopedEnvironmentVariable&
        operator=(const ScopedEnvironmentVariable&) = delete;

      private:
        std::string m_name;
        std::optional<std::string> m_original;
    };

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
        return sasktran_disco::la::dgbsv_configured(
            system.size, system.bandwidth, system.bandwidth,
            system.right_hand_sides, system.matrix.data(),
            system.leading_dimension, pivots.data(),
            system.right_hand_side.data(), system.size);
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

    void benchmark_runtime_selection(int num_stokes, int num_streams,
                                     int num_layers) {
        BENCHMARK("construction-time backend selection") {
            sasktran2::Config config;
            config.set_num_stokes(num_stokes);
            config.set_num_do_streams(num_streams);
            config.set_multiple_scatter_source(
                sasktran2::Config::MultipleScatterSource::discrete_ordinates);
            sasktran2::detail::RuntimeBackendTuner::resolve(config, num_layers);
            return sasktran2::detail::RuntimeBackendTuner::
                resolved_banded_lu_backend(config);
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
#if !defined(SKTRAN_USE_ACCELERATE) && !defined(SKTRAN_USE_MKL)
        REQUIRE(custom_pivots == configured_lapack_pivots);
#endif
        for (std::size_t index = 0; index < custom.right_hand_side.size();
             ++index) {
            REQUIRE(custom.right_hand_side[index] ==
                    Catch::Approx(configured_lapack.right_hand_side[index])
                        .epsilon(2.0e-13)
                        .margin(2.0e-13));
        }
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

TEST_CASE("Runtime band factorization policy prefers LAPACK for close results",
          "[sktran_do][band_factorization]") {
    using Backend = sasktran2::detail::BandedLUBackend;
    using Tuner = sasktran2::detail::RuntimeBackendTuner;

    REQUIRE(Tuner::select_banded_lu_backend(100.0, 80.0) == Backend::unblocked);
    REQUIRE(Tuner::select_banded_lu_backend(100.0, 90.0) == Backend::lapack);
    REQUIRE(Tuner::select_banded_lu_backend(100.0, 95.0) == Backend::lapack);
    REQUIRE(Tuner::select_banded_lu_backend(100.0, 105.0) == Backend::lapack);
    REQUIRE(Tuner::select_banded_lu_backend(0.0, 1.0) == Backend::lapack);
}

TEST_CASE("Runtime band factorization tuner covers every DISCO source path",
          "[sktran_do][band_factorization]") {
    using Backend = sasktran2::detail::BandedLUBackend;
    using Source = sasktran2::Config::MultipleScatterSource;
    using Tuner = sasktran2::detail::RuntimeBackendTuner;

    const ScopedEnvironmentVariable enable_unblocked(
        "SASKTRAN2_DISABLE_DO_UNBLOCKED_BAND_LU", std::nullopt);
    const ScopedEnvironmentVariable force_unblocked(
        "SASKTRAN2_DO_BANDED_LU_BACKEND", std::string("unblocked"));

    SECTION("direct scalar two-stream source honors the explicit override") {
        sasktran2::Config config;
        config.set_num_stokes(1);
        config.set_num_do_streams(2);
        config.set_multiple_scatter_source(Source::discrete_ordinates);

        Tuner::resolve(config, 20);

        REQUIRE(Tuner::resolved_banded_lu_backend(config) ==
                Backend::unblocked);
    }

    SECTION("HR initialized with discrete ordinates is tuned") {
        sasktran2::Config config;
        config.set_num_do_streams(8);
        config.set_multiple_scatter_source(Source::successive_orders_legacy);
        config.set_initialize_hr_with_do(true);

        Tuner::resolve(config, 20);

        REQUIRE(Tuner::resolved_banded_lu_backend(config) ==
                Backend::unblocked);
    }

    SECTION("HR without discrete-ordinates initialization is not tuned") {
        sasktran2::Config config;
        config.set_num_do_streams(8);
        config.set_multiple_scatter_source(Source::successive_orders_legacy);
        config.set_initialize_hr_with_do(false);

        Tuner::resolve(config, 20);

        REQUIRE(Tuner::resolved_banded_lu_backend(config) == Backend::lapack);
    }
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

TEST_CASE("Runtime backend selection overhead",
          "[.benchmark][band_backend_selection]") {
    SECTION("scalar 4 streams, 20 layers") {
        benchmark_runtime_selection(1, 4, 20);
    }
    SECTION("scalar 8 streams, 100 layers") {
        benchmark_runtime_selection(1, 8, 100);
    }
    SECTION("scalar 16 streams, 100 layers") {
        benchmark_runtime_selection(1, 16, 100);
    }
    SECTION("IQU 4 streams, 100 layers") {
        benchmark_runtime_selection(3, 4, 100);
    }
}
