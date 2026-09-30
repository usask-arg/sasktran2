#include <sasktran2/solartransmission.h>
#include <sasktran2/test_helper.h>

#include <cstring>

namespace {
    using Rows = std::vector<std::vector<std::pair<int, double>>>;
    using sasktran2::solartransmission::SolarTableInterpolation;

    void require_solar_bits(const Eigen::VectorXd& actual,
                            const Eigen::VectorXd& reference) {
        REQUIRE(actual.size() == reference.size());
        if (actual.size() != 0) {
            REQUIRE(std::memcmp(actual.data(), reference.data(),
                                actual.size() * sizeof(double)) == 0);
        }
    }

    void check_solar_products(SolarTableInterpolation& interpolation,
                              const Rows& rows, int columns) {
        Eigen::VectorXd table(columns);
        Eigen::VectorXd cotangent(rows.size());
        if (columns != 0) {
            table.setLinSpaced(-0.35, 0.47);
        }
        if (!rows.empty()) {
            cotangent.setLinSpaced(-0.17, 0.29);
        }
        Eigen::VectorXd expected_forward(rows.size());
        Eigen::VectorXd expected_reverse = Eigen::VectorXd::Zero(columns);
        for (std::size_t row = 0; row < rows.size(); ++row) {
            double result = 0.0;
            for (const auto& [column, weight] : rows[row]) {
                if (weight == 0.0) {
                    continue;
                }
                result += weight * table(column);
                expected_reverse(column) += weight * cotangent(row);
            }
            expected_forward(row) = result;
        }
        Eigen::VectorXd forward(rows.size());
        Eigen::VectorXd reverse(columns);
        interpolation.apply(table, forward);
        interpolation.apply_transpose(cotangent, reverse);
        require_solar_bits(forward, expected_forward);
        require_solar_bits(reverse, expected_reverse);
        REQUIRE(forward.dot(cotangent) ==
                Catch::Approx(table.dot(reverse)).margin(2.0e-12));
        // Repeated products do not change storage, cursors, or row ordering.
        forward.setConstant(31.0);
        reverse.setConstant(-37.0);
        interpolation.apply(table, forward);
        interpolation.apply_transpose(cotangent, reverse);
        require_solar_bits(forward, expected_forward);
        require_solar_bits(reverse, expected_reverse);
    }

    void initialize_solar(SolarTableInterpolation& interpolation,
                          const Rows& rows, int columns) {
        interpolation.initialize(rows.size(), columns, 1024);
        for (const auto& row : rows) {
            interpolation.append_row(row);
        }
        interpolation.finalize();
        interpolation.finalize();
    }
} // namespace

TEST_CASE("Solar interpolation row-relative indices preserve complete products",
          "[solar_interpolation][compact_indices]") {
    SolarTableInterpolation interpolation;
    // The relative span includes 65535 exactly, while the global columns are
    // wider. Unsorted duplicate columns retain their original update order.
    const Rows rows{{},
                    {{65535, 0.3},
                     {131070, -0.17},
                     {70000, 0.29},
                     {65535, 0.09},
                     {90000, 0.0},
                     {65536, -0.0}},
                    {{131071, 0.11},
                     {131060, -0.03},
                     {131065, 0.09},
                     {131071, 0.4},
                     {131065, 0.001},
                     {131060, -0.2}}};
    initialize_solar(interpolation, rows, 131072);
    REQUIRE(interpolation.relative_column_indices());
    REQUIRE(interpolation.compact_row_counts());
    REQUIRE(interpolation.non_zeros() == 10);
    check_solar_products(interpolation, rows, 131072);
}

TEST_CASE(
    "Solar interpolation falls back at row-count and index-span boundaries",
    "[solar_interpolation][compact_indices]") {
    SolarTableInterpolation interpolation;
    for (const int count : {255, 256}) {
        for (const bool wide_span : {false, true}) {
            Rows rows(2);
            for (int entry = 0; entry < count; ++entry) {
                const int column = wide_span && entry % 2 != 0 ? 65536 : 0;
                rows[1].emplace_back(column, entry % 2 != 0 ? -0.25 : 0.125);
            }
            initialize_solar(interpolation, rows, 65537);
            REQUIRE(interpolation.compact_row_counts() == (count == 255));
            REQUIRE(interpolation.relative_column_indices() == !wide_span);
            REQUIRE(interpolation.non_zeros() == count);
            check_solar_products(interpolation, rows, 65537);
        }
    }
    // Changing dimensions and formats releases stale bases and row counts.
    const Rows empty_rows(2);
    initialize_solar(interpolation, empty_rows, 0);
    REQUIRE_FALSE(interpolation.relative_column_indices());
    REQUIRE(interpolation.compact_row_counts());
    check_solar_products(interpolation, empty_rows, 0);
    interpolation.clear();
    REQUIRE(interpolation.rows() == 0);
    REQUIRE(interpolation.cols() == 0);
    REQUIRE(interpolation.storage_bytes() == 0);
    REQUIRE_FALSE(interpolation.relative_column_indices());
    REQUIRE_FALSE(interpolation.compact_row_counts());
    initialize_solar(interpolation, {}, 0);
    check_solar_products(interpolation, {}, 0);
}

TEST_CASE("Solar interpolation validates original columns before narrowing",
          "[solar_interpolation][compact_indices]") {
    SolarTableInterpolation interpolation;
    interpolation.initialize(1, 65536, 4);
    REQUIRE_THROWS_AS(interpolation.append_row({{65536, 0.1}}),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(interpolation.append_row({{-1, 0.1}}),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(interpolation.finalize(), std::logic_error);
    const Rows sparse_rows{{}, {{65535, 0.2}}};
    initialize_solar(interpolation, sparse_rows, 65536);
    // One narrow index would save less than two row bases cost.
    REQUIRE_FALSE(interpolation.relative_column_indices());
    check_solar_products(interpolation, sparse_rows, 65536);
}
