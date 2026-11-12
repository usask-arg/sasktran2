#pragma once

#ifdef SKTRAN_CATCH2_VERSION3
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_all.hpp>

#else
#include <catch2/catch.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <limits>

namespace sasktran2::testing {
    /** One unit in the last place at the magnitude of x. */
    inline double ulp(double x) {
        x = std::abs(x);
        if (!std::isfinite(x)) {
            return std::numeric_limits<double>::infinity();
        }
        return std::nextafter(x, std::numeric_limits<double>::infinity()) - x;
    }

    /** Whether actual and expected differ by at most max_ulps units in the
     * last place of the largest of |scale|, |actual| and |expected|.
     *
     * Two code paths that perform the same operations in the same order give
     * identical bits only if the compiler also rounds them identically, which
     * GCC 15 builds on CI runners do not always do. For an accumulated sum,
     * pass the sum of the absolute values of its terms as scale, so that
     * cancellation in the result does not shrink the tolerance below the
     * rounding of the terms.
     */
    inline bool within_ulps(double actual, double expected, double max_ulps,
                            double scale = 0.0) {
        if (actual == expected) {
            return true;
        }
        if (!std::isfinite(actual) || !std::isfinite(expected)) {
            return false;
        }
        const double magnitude =
            std::max({std::abs(scale), std::abs(actual), std::abs(expected)});
        return std::abs(actual - expected) <= max_ulps * ulp(magnitude);
    }

    /** The difference between actual and expected in units in the last place
     * of the larger of the two, for failure messages. */
    inline double ulp_difference(double actual, double expected) {
        return std::abs(actual - expected) /
               ulp(std::max(std::abs(actual), std::abs(expected)));
    }
} // namespace sasktran2::testing
