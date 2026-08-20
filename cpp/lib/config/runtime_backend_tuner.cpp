#include "sasktran2/runtime_backend_tuner.h"

#include "sktran_disco/sktran_do.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
    using Backend = sasktran2::detail::BandedLUBackend;
    using Clock = std::chrono::steady_clock;

    struct BandSystem {
        lapack_int size;
        lapack_int bandwidth;
        lapack_int leading_dimension;
        std::vector<double> matrix;
        std::vector<double> right_hand_side;
        std::vector<lapack_int> pivots;

        BandSystem(lapack_int system_size, lapack_int system_bandwidth)
            : size(system_size), bandwidth(system_bandwidth),
              leading_dimension(3 * system_bandwidth + 1),
              right_hand_side(system_size), pivots(system_size) {
            const auto matrix_size =
                static_cast<std::size_t>(leading_dimension) *
                static_cast<std::size_t>(size);
            if (size <= 0 || bandwidth < 0 ||
                matrix_size / static_cast<std::size_t>(size) !=
                    static_cast<std::size_t>(leading_dimension)) {
                throw std::length_error("Invalid band benchmark dimensions");
            }
            matrix.resize(matrix_size);
        }

        void reset() {
            std::fill(matrix.begin(), matrix.end(), 0.0);
            std::uint64_t random_state = 0x243f6a8885a308d3ULL;
            const lapack_int diagonal_row = 2 * bandwidth;
            for (lapack_int column = 0; column < size; ++column) {
                const lapack_int first_row =
                    std::max<lapack_int>(0, column - bandwidth);
                const lapack_int last_row =
                    std::min<lapack_int>(size - 1, column + bandwidth);
                for (lapack_int row = first_row; row <= last_row; ++row) {
                    random_state = random_state * 6364136223846793005ULL + 1;
                    const double random_value =
                        static_cast<double>((random_state >> 40) & 0xffffU) /
                            65535.0 -
                        0.5;
                    double value = 0.03 * random_value;
                    if (row == column) {
                        value +=
                            2.5 + 0.001 * static_cast<double>(column % 101);
                    }
                    // Exercise the row-interchange path periodically, as the
                    // physical BVP is not guaranteed to be diagonal dominant.
                    if (column % 17 == 0 && column + 1 < size) {
                        if (row == column) {
                            value = 0.25;
                        } else if (row == column + 1) {
                            value = 3.0;
                        }
                    }
                    const lapack_int band_row = diagonal_row + row - column;
                    matrix[column * leading_dimension + band_row] = value;
                }
            }
            for (lapack_int row = 0; row < size; ++row) {
                right_hand_side[row] =
                    std::sin(0.13 * static_cast<double>(row + 1));
            }
        }
    };

    int solve(BandSystem& system, Backend backend) {
        if (backend == Backend::unblocked) {
            return sasktran_disco::la::dgbsv_unblocked(
                system.size, system.bandwidth, system.bandwidth, 1,
                system.matrix.data(), system.leading_dimension,
                system.pivots.data(), system.right_hand_side.data(),
                system.size);
        }
        return sasktran_disco::la::dgbsv_configured(
            system.size, system.bandwidth, system.bandwidth, 1,
            system.matrix.data(), system.leading_dimension,
            system.pivots.data(), system.right_hand_side.data(), system.size);
    }

    int timed_solve(BandSystem& system, Backend backend,
                    double& elapsed_nanoseconds) {
        system.reset();
        const auto start = Clock::now();
        const int info = solve(system, backend);
        const auto stop = Clock::now();
        elapsed_nanoseconds =
            std::chrono::duration<double, std::nano>(stop - start).count();
        return info;
    }

    double median(std::vector<double> values) {
        const auto middle = values.begin() + values.size() / 2;
        std::nth_element(values.begin(), middle, values.end());
        return *middle;
    }

    struct BenchmarkResult {
        double lapack_nanoseconds = 0.0;
        double unblocked_nanoseconds = 0.0;
        bool valid = false;
    };

    BenchmarkResult benchmark_band_solvers(lapack_int size,
                                           lapack_int bandwidth) {
        BandSystem system(size, bandwidth);

        // Warm both implementations and use the same synthetic system as a
        // cheap correctness check before trusting the timing result.
        system.reset();
        const int lapack_warm_info = solve(system, Backend::lapack);
        const auto lapack_solution = system.right_hand_side;
        system.reset();
        const int unblocked_warm_info = solve(system, Backend::unblocked);
        if (lapack_warm_info != 0 || unblocked_warm_info != lapack_warm_info) {
            return {};
        }
        for (std::size_t index = 0; index < lapack_solution.size(); ++index) {
            const double scale =
                std::max(1.0, std::abs(lapack_solution[index]));
            if (std::abs(system.right_hand_side[index] -
                         lapack_solution[index]) > 1.0e-11 * scale) {
                return {};
            }
        }

        std::vector<double> lapack_times;
        std::vector<double> unblocked_times;
        lapack_times.reserve(7);
        unblocked_times.reserve(7);

        double lapack_time;
        double unblocked_time;
        if (timed_solve(system, Backend::lapack, lapack_time) != 0 ||
            timed_solve(system, Backend::unblocked, unblocked_time) != 0) {
            return {};
        }
        lapack_times.push_back(lapack_time);
        unblocked_times.push_back(unblocked_time);

        const double slower_time = std::max(lapack_time, unblocked_time);
        const int sample_count = slower_time < 2.0e5   ? 7
                                 : slower_time < 2.0e6 ? 5
                                 : slower_time < 1.0e7 ? 3
                                                       : 1;
        for (int sample = 1; sample < sample_count; ++sample) {
            int lapack_info;
            int unblocked_info;
            if (sample % 2 == 0) {
                lapack_info = timed_solve(system, Backend::lapack, lapack_time);
                unblocked_info =
                    timed_solve(system, Backend::unblocked, unblocked_time);
            } else {
                unblocked_info =
                    timed_solve(system, Backend::unblocked, unblocked_time);
                lapack_info = timed_solve(system, Backend::lapack, lapack_time);
            }
            if (lapack_info != 0 || unblocked_info != 0) {
                return {};
            }
            lapack_times.push_back(lapack_time);
            unblocked_times.push_back(unblocked_time);
        }

        return {median(lapack_times), median(unblocked_times), true};
    }

    bool environment_requests(const char* value, const char* requested) {
        return value != nullptr && std::strcmp(value, requested) == 0;
    }
} // namespace

namespace sasktran2::detail {
    BandedLUBackend RuntimeBackendTuner::select_banded_lu_backend(
        double lapack_nanoseconds, double unblocked_nanoseconds) {
        if (!std::isfinite(lapack_nanoseconds) ||
            !std::isfinite(unblocked_nanoseconds) ||
            lapack_nanoseconds <= 0.0 || unblocked_nanoseconds <= 0.0) {
            return Backend::lapack;
        }
        return unblocked_nanoseconds < 0.90 * lapack_nanoseconds
                   ? Backend::unblocked
                   : Backend::lapack;
    }

    BandedLUBackend
    RuntimeBackendTuner::resolved_banded_lu_backend(const Config& config) {
        return config.m_runtime_backends.do_banded_lu;
    }

    void RuntimeBackendTuner::resolve(Config& config, int num_layers) {
        config.m_runtime_backends.do_banded_lu = Backend::lapack;
        const auto multiple_scatter_source = config.multiple_scatter_source();
        const bool uses_discrete_ordinates =
            multiple_scatter_source ==
                Config::MultipleScatterSource::discrete_ordinates ||
            (multiple_scatter_source ==
                 Config::MultipleScatterSource::successive_orders_legacy &&
             config.initialize_hr_with_do());
        if (!uses_discrete_ordinates || num_layers <= 0) {
            return;
        }

        // Preserve the original opt-out while also providing a deterministic
        // developer override for reproducing either path.
        if (std::getenv("SASKTRAN2_DISABLE_DO_UNBLOCKED_BAND_LU") != nullptr) {
            return;
        }
        const char* override_backend =
            std::getenv("SASKTRAN2_DO_BANDED_LU_BACKEND");
        if (environment_requests(override_backend, "lapack")) {
            return;
        }
        if (environment_requests(override_backend, "unblocked")) {
            config.m_runtime_backends.do_banded_lu = Backend::unblocked;
            return;
        }

        const int num_stokes = config.num_stokes();
        const int num_streams = config.num_do_streams();
        // Automatic tuning is unnecessary for the common scalar two-stream
        // path, which uses its dedicated pentadiagonal solver. The explicit
        // override above remains available to exercise generic paths.
        if (num_stokes == 1 && num_streams == 2) {
            return;
        }

        const auto size =
            static_cast<long long>(num_streams) * num_stokes * num_layers;
        const auto bandwidth = 3LL * num_streams * num_stokes / 2LL - 1LL;
        if (size <= 0 || bandwidth <= 2 ||
            size > std::numeric_limits<lapack_int>::max() ||
            bandwidth > std::numeric_limits<lapack_int>::max()) {
            return;
        }

        try {
            const auto result =
                benchmark_band_solvers(static_cast<lapack_int>(size),
                                       static_cast<lapack_int>(bandwidth));
            if (!result.valid) {
                return;
            }
            config.m_runtime_backends.do_banded_lu = select_banded_lu_backend(
                result.lapack_nanoseconds, result.unblocked_nanoseconds);
            spdlog::debug(
                "DISCO band LU autotune (N={}, bandwidth={}): LAPACK {:.1f} "
                "us, unblocked {:.1f} us; selected {}",
                size, bandwidth, result.lapack_nanoseconds / 1000.0,
                result.unblocked_nanoseconds / 1000.0,
                config.m_runtime_backends.do_banded_lu == Backend::unblocked
                    ? "unblocked"
                    : "LAPACK");
        } catch (const std::exception& error) {
            spdlog::debug("DISCO band LU autotune failed; using LAPACK: {}",
                          error.what());
        }
    }
} // namespace sasktran2::detail
