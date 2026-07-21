#include <sasktran2/test_helper.h>
#include <sasktran2/source_algorithms.h>

#include <algorithm>
#include <cmath>

#ifdef SKTRAN_CATCH2_VERSION3

TEST_CASE("Linear source weights match analytic quadrature",
          "[sourceintegrator][sourcealgo]") {
    for (const double od : {1e-12, 1e-8, 1e-4, 0.1, 1.0, 10.0}) {
        const double attenuation = std::exp(-od);
        const auto weights =
            sasktran2::sourcealgo::linear_source_weights(od, attenuation);
        double expected_start;
        double expected_end;
        if (od < 1e-2) {
            const double od2 = od * od;
            const double od3 = od2 * od;
            const double od4 = od3 * od;
            const double od5 = od4 * od;
            const double od6 = od5 * od;
            expected_start = od / 2.0 - od2 / 6.0 + od3 / 24.0 - od4 / 120.0 +
                             od5 / 720.0 - od6 / 5040.0;
            expected_end = od / 2.0 - od2 / 3.0 + od3 / 8.0 - od4 / 30.0 +
                           od5 / 144.0 - od6 / 840.0;
        } else {
            const double phi = -std::expm1(-od) / od;
            expected_start = 1.0 - phi;
            expected_end = phi - attenuation;
        }

        REQUIRE(weights.start == Catch::Approx(expected_start).epsilon(1e-11));
        REQUIRE(weights.end == Catch::Approx(expected_end).epsilon(1e-11));
        REQUIRE(weights.start + weights.end ==
                Catch::Approx(-std::expm1(-od)).epsilon(1e-12));

        const double step = std::max(1e-7 * od, 1e-10);
        const auto above = sasktran2::sourcealgo::linear_source_weights(
            od + step, std::exp(-(od + step)));
        const auto below = sasktran2::sourcealgo::linear_source_weights(
            od - step, std::exp(-(od - step)));
        REQUIRE(weights.d_start_d_od ==
                Catch::Approx((above.start - below.start) / (2.0 * step))
                    .epsilon(2e-6));
        REQUIRE(weights.d_end_d_od ==
                Catch::Approx((above.end - below.end) / (2.0 * step))
                    .epsilon(2e-6));
    }
}

TEST_CASE("Linear source weights retain curved-layer quadrature",
          "[sourceintegrator][sourcealgo]") {
    const double od = 1e-10;
    const auto weights = sasktran2::sourcealgo::linear_source_weights(
        od, std::exp(-od), 0.25, 0.75);

    REQUIRE(weights.start / od == Catch::Approx(0.25).margin(1e-10));
    REQUIRE(weights.end / od == Catch::Approx(0.75).margin(1e-10));
}

#endif
