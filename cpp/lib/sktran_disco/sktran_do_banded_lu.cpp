#include "sktran_disco/sktran_do.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>

namespace sasktran_disco::la {
    int dgbtf2_unblocked(lapack_int n, lapack_int lower_bandwidth,
                         lapack_int upper_bandwidth, double* matrix,
                         lapack_int leading_dimension, lapack_int* pivots) {
        if (n < 0) {
            return -1;
        }
        if (lower_bandwidth < 0) {
            return -2;
        }
        if (upper_bandwidth < 0) {
            return -3;
        }
        if (leading_dimension < 2 * lower_bandwidth + upper_bandwidth + 1) {
            return -5;
        }
        if (n == 0) {
            return 0;
        }

        const lapack_int diagonal_row = lower_bandwidth + upper_bandwidth;
        for (lapack_int column = upper_bandwidth + 1;
             column < std::min(diagonal_row, n); ++column) {
            const lapack_int first_row = diagonal_row - column;
            std::fill(matrix + column * leading_dimension + first_row,
                      matrix + column * leading_dimension + lower_bandwidth,
                      0.0);
        }

        lapack_int info = 0;
        lapack_int last_updated_column = 0;
        lapack_int fill_column = diagonal_row;
        for (lapack_int column = 0; column < n; ++column) {
            if (fill_column < n) {
                std::fill_n(matrix + fill_column * leading_dimension,
                            lower_bandwidth, 0.0);
                ++fill_column;
            }

            const lapack_int lower_count =
                std::min(lower_bandwidth, n - 1 - column);
            const lapack_int column_start =
                column * leading_dimension + diagonal_row;
            lapack_int pivot_offset = 0;
            double pivot_absolute = std::abs(matrix[column_start]);
            for (lapack_int offset = 1; offset <= lower_count; ++offset) {
                const double candidate =
                    std::abs(matrix[column_start + offset]);
                if (candidate > pivot_absolute) {
                    pivot_absolute = candidate;
                    pivot_offset = offset;
                }
            }
            pivots[column] = column + pivot_offset + 1;

            const double pivot = matrix[column_start + pivot_offset];
            if (pivot == 0.0) {
                if (info == 0) {
                    info = column + 1;
                }
                continue;
            }

            last_updated_column = std::max(
                last_updated_column,
                std::min(column + upper_bandwidth + pivot_offset, n - 1));
            if (pivot_offset != 0) {
                for (lapack_int offset = 0;
                     offset <= last_updated_column - column; ++offset) {
                    const lapack_int pivot_index =
                        (column + offset) * leading_dimension + diagonal_row +
                        pivot_offset - offset;
                    const lapack_int diagonal_index =
                        (column + offset) * leading_dimension + diagonal_row -
                        offset;
                    std::swap(matrix[pivot_index], matrix[diagonal_index]);
                }
            }

            if (lower_count == 0) {
                continue;
            }
            const double diagonal = matrix[column_start];
            if (std::abs(diagonal) >= std::numeric_limits<double>::min()) {
                const double inverse_diagonal = 1.0 / diagonal;
                for (lapack_int offset = 1; offset <= lower_count; ++offset) {
                    matrix[column_start + offset] *= inverse_diagonal;
                }
            } else {
                for (lapack_int offset = 1; offset <= lower_count; ++offset) {
                    matrix[column_start + offset] /= diagonal;
                }
            }

            for (lapack_int column_offset = 1;
                 column_offset <= last_updated_column - column;
                 ++column_offset) {
                const lapack_int target_start =
                    (column + column_offset) * leading_dimension +
                    diagonal_row - column_offset;
                const double upper = matrix[target_start];
                if (upper == 0.0) {
                    continue;
                }
                const double scaled_upper = -upper;
                for (lapack_int offset = 1; offset <= lower_count; ++offset) {
                    matrix[target_start + offset] +=
                        matrix[column_start + offset] * scaled_upper;
                }
            }
        }
        return info;
    }

    int dgbsv_unblocked(lapack_int n, lapack_int lower_bandwidth,
                        lapack_int upper_bandwidth, lapack_int right_hand_sides,
                        double* matrix, lapack_int leading_dimension,
                        lapack_int* pivots, double* right_hand_side,
                        lapack_int right_hand_side_leading_dimension) {
        if (n < 0) {
            return -1;
        }
        if (lower_bandwidth < 0) {
            return -2;
        }
        if (upper_bandwidth < 0) {
            return -3;
        }
        if (right_hand_sides < 0) {
            return -4;
        }
        if (leading_dimension < 2 * lower_bandwidth + upper_bandwidth + 1) {
            return -6;
        }
        if (right_hand_side_leading_dimension < std::max<lapack_int>(1, n)) {
            return -9;
        }
        lapack_int info = dgbtf2_unblocked(n, lower_bandwidth, upper_bandwidth,
                                           matrix, leading_dimension, pivots);
        if (info != 0 || n == 0 || right_hand_sides == 0) {
            return info;
        }

#if defined(SKTRAN_USE_ACCELERATE) || defined(SKTRAN_NO_LAPACKE)
        char transpose = 'N';
        dgbtrs_(&transpose, &n, &lower_bandwidth, &upper_bandwidth,
                &right_hand_sides, matrix, &leading_dimension, pivots,
                right_hand_side, &right_hand_side_leading_dimension, &info);
#else
        info = LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N', n, lower_bandwidth,
                              upper_bandwidth, right_hand_sides, matrix,
                              leading_dimension, pivots, right_hand_side,
                              right_hand_side_leading_dimension);
#endif
        return info;
    }

    bool use_unblocked_band_factorization(lapack_int lower_bandwidth,
                                          lapack_int upper_bandwidth) {
#if defined(SKTRAN_USE_ACCELERATE) || defined(SKTRAN_USE_MKL)
        // Keep vendor-tuned implementations until this crossover has been
        // measured for them independently.
        return false;
#else
        static const bool disabled =
            std::getenv("SASKTRAN2_DISABLE_DO_UNBLOCKED_BAND_LU") != nullptr;
        // Production DO bandwidths advance in steps of three. Benchmarks show
        // the custom unblocked loop wins through 29, while the configured
        // LAPACK's blocked path catches up at the next bandwidth (32).
        return !disabled && lower_bandwidth == upper_bandwidth &&
               lower_bandwidth > 2 && lower_bandwidth < 32;
#endif
    }
} // namespace sasktran_disco::la
