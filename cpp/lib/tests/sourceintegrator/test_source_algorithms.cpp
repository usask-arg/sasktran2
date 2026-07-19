#include <sasktran2/test_helper.h>
#include <sasktran2/source_algorithms.h>

#include <algorithm>
#include <array>
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

TEST_CASE("Single-scatter weights include log-linear solar transmission",
          "[sourceintegrator][sourcealgo]") {
    constexpr double start_fraction = 0.35;
    constexpr double end_fraction = 0.65;

    for (const auto& optical_depths :
         {std::array{0.4, 0.2, 0.7}, std::array{0.4, 0.7, 0.3},
          std::array{1e-8, 0.4, 0.4}}) {
        const double view_od = optical_depths[0];
        const double solar_start_od = optical_depths[1];
        const double solar_end_od = optical_depths[2];
        const double solar_start = std::exp(-solar_start_od);
        const double solar_end = std::exp(-solar_end_od);
        const auto evaluate = [](double current_view_od,
                                 double current_solar_start_od,
                                 double current_solar_end_od) {
            return sasktran2::sourcealgo::single_scatter_source_weights(
                current_view_od, std::exp(-current_view_od),
                std::exp(-current_solar_start_od),
                std::exp(-current_solar_end_od), start_fraction, end_fraction);
        };
        const auto weights =
            sasktran2::sourcealgo::single_scatter_source_weights(
                view_od, std::exp(-view_od), solar_start, solar_end,
                start_fraction, end_fraction);

        const double combined_od = view_od + solar_end_od - solar_start_od;
        const double expected_constant_source =
            std::abs(combined_od) < 1e-10
                ? solar_start
                : solar_start * -std::expm1(-combined_od) / combined_od;
        REQUIRE(weights.start + weights.end ==
                Catch::Approx(expected_constant_source).epsilon(2e-12));

        constexpr double step = 1e-6;
        const auto view_above =
            evaluate(view_od + step, solar_start_od, solar_end_od);
        const auto view_below =
            evaluate(view_od - step, solar_start_od, solar_end_od);
        const auto start_above =
            evaluate(view_od, solar_start_od + step, solar_end_od);
        const auto start_below =
            evaluate(view_od, solar_start_od - step, solar_end_od);
        const auto end_above =
            evaluate(view_od, solar_start_od, solar_end_od + step);
        const auto end_below =
            evaluate(view_od, solar_start_od, solar_end_od - step);

        REQUIRE(
            weights.d_start_d_view_od ==
            Catch::Approx((view_above.start - view_below.start) / (2.0 * step))
                .epsilon(2e-7));
        REQUIRE(weights.d_end_d_view_od ==
                Catch::Approx((view_above.end - view_below.end) / (2.0 * step))
                    .epsilon(2e-7));
        REQUIRE(weights.d_start_d_solar_start_od ==
                Catch::Approx((start_above.start - start_below.start) /
                              (2.0 * step))
                    .epsilon(2e-7));
        REQUIRE(
            weights.d_end_d_solar_start_od ==
            Catch::Approx((start_above.end - start_below.end) / (2.0 * step))
                .epsilon(2e-7));
        REQUIRE(
            weights.d_start_d_solar_end_od ==
            Catch::Approx((end_above.start - end_below.start) / (2.0 * step))
                .epsilon(2e-7));
        REQUIRE(weights.d_end_d_solar_end_od ==
                Catch::Approx((end_above.end - end_below.end) / (2.0 * step))
                    .epsilon(2e-7));
    }
}

TEST_CASE("Single-scatter weights remain finite for a large negative rate",
          "[sourceintegrator][sourcealgo]") {
    const double view_od = 1.0;
    const auto weights = sasktran2::sourcealgo::single_scatter_source_weights(
        view_od, std::exp(-view_od), 1e-300, 1.0);

    REQUIRE(std::isfinite(weights.start));
    REQUIRE(std::isfinite(weights.end));
    REQUIRE(std::isfinite(weights.d_start_d_view_od));
    REQUIRE(std::isfinite(weights.d_end_d_view_od));
}

#endif
