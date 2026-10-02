#pragma once

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace sasktran2::successive_orders {
    /** Horizontal interpolation of the stored source between columns. */
    enum class HorizontalInterpolation { linear, cubic };

    /** Largest sum of absolute cubic weights accepted for one stencil.
     *
     * Lagrange weights grow without bound on strongly non-uniform grids, for
     * example 0.25, -1, 1.5, 0.25 on [0, 1, 2, 4] at 3 and about +-38 on
     * [0, 1, 1.01, 2] at 0.5, which amplifies any error in the stored
     * source. Stencils whose absolute weights sum to more than this use
     * linear weights instead. On a uniform grid the sum is at most 1.25 in
     * interior intervals and about 1.63 in the end intervals, so uniform
     * grids always stay cubic.
     */
    inline constexpr double max_cubic_weight_abs_sum = 2.0;

    /** Whether four-point weights satisfy max_cubic_weight_abs_sum. */
    inline bool
    cubic_weights_are_bounded(const std::array<double, 4>& weights) {
        double sum = 0.0;
        for (const double weight : weights) {
            sum += std::abs(weight);
        }
        return sum <= max_cubic_weight_abs_sum;
    }

    /** Four-point Lagrange weights on a strictly increasing grid.
     *
     * Requires at least four nodes (std::invalid_argument otherwise) and
     * grid[0] < x < grid[n-1]. The stencil is centred on the interval
     * containing x and shifted inward at the ends. Weights reproduce cubic
     * polynomials and are exactly one/zero at nodes.
     */
    inline void cubic_lagrange_weights(const Eigen::VectorXd& grid, double x,
                                       std::array<int, 4>& indices,
                                       std::array<double, 4>& weights) {
        const int size = static_cast<int>(grid.size());
        if (size < 4) {
            throw std::invalid_argument(
                "Cubic Lagrange weights require at least four grid nodes");
        }
        const int lower = static_cast<int>(
            std::upper_bound(grid.data(), grid.data() + size, x) - grid.data() -
            1);
        const int start = std::clamp(lower - 1, 0, size - 4);
        for (int m = 0; m < 4; ++m) {
            double weight = 1.0;
            for (int n = 0; n < 4; ++n) {
                if (n != m) {
                    weight *= (x - grid[start + n]) /
                              (grid[start + m] - grid[start + n]);
                }
            }
            indices[m] = start + m;
            weights[m] = weight;
        }
    }
} // namespace sasktran2::successive_orders
