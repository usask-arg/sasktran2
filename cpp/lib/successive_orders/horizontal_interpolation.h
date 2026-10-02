#pragma once

#include <Eigen/Core>

#include <algorithm>
#include <array>

namespace sasktran2::successive_orders {
    /** Horizontal interpolation of the stored source between columns. */
    enum class HorizontalInterpolation { linear, cubic };

    /** Four-point Lagrange weights on a strictly increasing grid.
     *
     * Requires at least four nodes and grid[0] < x < grid[n-1]. The stencil
     * is centred on the interval containing x and shifted inward at the ends.
     * Weights reproduce cubic polynomials and are exactly one/zero at nodes.
     */
    inline void cubic_lagrange_weights(const Eigen::VectorXd& grid, double x,
                                       std::array<int, 4>& indices,
                                       std::array<double, 4>& weights) {
        const int size = static_cast<int>(grid.size());
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
