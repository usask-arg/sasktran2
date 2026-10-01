#include "../../successive_orders/endpoint_stencil_storage.h"

#include <sasktran2/test_helper.h>

#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

namespace {
    using namespace sasktran2::successive_orders;

    std::uint64_t bits(double value) {
        std::uint64_t result;
        std::memcpy(&result, &value, sizeof(result));
        return result;
    }

    double from_bits(std::uint64_t value) {
        double result;
        std::memcpy(&result, &value, sizeof(result));
        return result;
    }

    std::pair<int, double> pair(const InterpolationWeight& value) {
        return {value.index, value.weight()};
    }
    std::pair<int, double> pair(const std::pair<int, double>& value) {
        return value;
    }

    void
    require_original_stencils(const EndpointStencilStorage& storage,
                              const std::vector<int>& offsets,
                              const std::vector<InterpolationWeight>& weights) {
        REQUIRE(storage.size() == offsets.size() - 1);
        REQUIRE(storage.weight_count() == weights.size());
        storage.visit([&](const auto& stencils) {
            for (std::size_t slot = 0; slot < storage.size(); ++slot) {
                const auto values = stencils[slot];
                REQUIRE(values.size() == offsets[slot + 1] - offsets[slot]);
                for (std::size_t index = 0; index < values.size(); ++index) {
                    const auto value = pair(values[index]);
                    const auto& original = weights[offsets[slot] + index];
                    REQUIRE(value.first == original.index);
                    REQUIRE(bits(value.second) == bits(original.weight()));
                }
            }
        });
        for (std::size_t slot = 0; slot < storage.size(); ++slot) {
            storage.view(slot).visit([&](const auto& values) {
                for (std::size_t index = 0; index < values.size(); ++index) {
                    const auto value = pair(values[index]);
                    const auto& original = weights[offsets[slot] + index];
                    REQUIRE(value.first == original.index);
                    REQUIRE(bits(value.second) == bits(original.weight()));
                }
            });
        }
    }

    void
    require_original_products(const EndpointStencilStorage& storage,
                              const std::vector<int>& offsets,
                              const std::vector<InterpolationWeight>& weights) {
        storage.visit([&](const auto& stencils) {
            for (std::size_t slot = 0; slot < storage.size(); ++slot) {
                double expected = 0.0;
                double expected_tangent = 0.0;
                for (int index = offsets[slot]; index < offsets[slot + 1];
                     ++index) {
                    const auto& value = weights[index];
                    expected += value.weight() * (0.013 * value.index + 0.07);
                    expected_tangent +=
                        value.weight() * (-0.003 * value.index + 0.01);
                }
                double actual = 0.0;
                double actual_tangent = 0.0;
                const auto values = stencils[slot];
                for (std::size_t index = 0; index < values.size(); ++index) {
                    const auto value = pair(values[index]);
                    actual += value.second * (0.013 * value.first + 0.07);
                    actual_tangent +=
                        value.second * (-0.003 * value.first + 0.01);
                }
                REQUIRE(bits(actual) == bits(expected));
                REQUIRE(bits(actual_tangent) == bits(expected_tangent));
            }
        });
    }

    std::vector<InterpolationWeight>
    corners(int base, int stride, const std::array<double, 4>& weights) {
        return {{base, weights[0]},
                {base + 1, weights[1]},
                {base + stride, weights[2]},
                {base + stride + 1, weights[3]}};
    }
} // namespace

TEST_CASE("Structured endpoint stencils preserve all indices and double bits",
          "[successive_orders][endpoint_storage]") {
    for (const int base : {12, 65535, 65536, 70000}) {
        const int stride = 10;
        CAPTURE(base);
        const std::vector<int> original_offsets{0, 4, 8};
        auto original = corners(base, stride, {0.125, -0.0, -0.375, 0.5});
        const auto special =
            corners(base == 65535 ? base - 20 : base + 20, stride,
                    {from_bits(1), from_bits(0x8000000000000001ULL),
                     from_bits(0x7ff8000000000123ULL), +0.0});
        original.insert(original.end(), special.begin(), special.end());
        auto offsets = original_offsets;
        auto weights = original;
        EndpointStencilStorage storage;
        storage.assign(std::move(offsets), std::move(weights), stride, 140000);
        REQUIRE(storage.is_structured());
        REQUIRE(storage.base_bits() == (base < 65536 ? 16 : 32));
        REQUIRE(storage.storage_bytes() == 2 * (base < 65536 ? 34 : 36));
        require_original_stencils(storage, original_offsets, original);

        original.resize(4);
        offsets = {0, 4};
        weights = original;
        storage.assign(std::move(offsets), std::move(weights), stride, 140000);
        require_original_products(storage, {0, 4}, original);
    }
}

TEST_CASE("Endpoint storage keeps generic stencil order and reinitializes",
          "[successive_orders][endpoint_storage]") {
    EndpointStencilStorage storage;
    const std::vector<std::vector<InterpolationWeight>> cases{
        {{3, 0.125}},
        {{3, 0.125}, {2, 0.875}},
        {{3, 0.125}, {2, 0.25}, {1, 0.625}},
        {{0, 0.125}, {10, 0.25}, {1, 0.375}, {11, 0.25}},
        {{9, 0.125}, {10, 0.25}, {19, 0.375}, {20, 0.25}},
    };
    for (const auto& original : cases) {
        CAPTURE(original.size());
        auto offsets = std::vector<int>{0, static_cast<int>(original.size())};
        const auto original_offsets = offsets;
        auto weights = original;
        storage.assign(std::move(offsets), std::move(weights), 10, 100);
        REQUIRE_FALSE(storage.is_structured());
        require_original_stencils(storage, original_offsets, original);
        require_original_products(storage, original_offsets, original);
    }
    // A single generic stencil keeps the whole engine in the original form.
    auto mixed = corners(0, 10, {0.125, 0.25, 0.375, 0.25});
    mixed.push_back({2, 1.0});
    auto original_mixed = mixed;
    auto offsets = std::vector<int>{0, 4, 5};
    storage.assign(std::move(offsets), std::move(mixed), 10, 100);
    REQUIRE_FALSE(storage.is_structured());
    require_original_stencils(storage, {0, 4, 5}, original_mixed);

    auto structured = corners(0, 10, {0.125, 0.25, 0.375, 0.25});
    offsets = {0, 4};
    storage.assign(std::move(offsets), std::move(structured), 10, 100);
    REQUIRE(storage.is_structured());
    storage.clear();
    REQUIRE(storage.empty());
    REQUIRE(storage.storage_bytes() == 0);
    REQUIRE_THROWS_AS(storage.view(0), std::out_of_range);

    offsets = {0, 4};
    std::vector<InterpolationWeight> invalid_cell{
        {std::numeric_limits<int>::max(), 0.25},
        {1, 0.25},
        {2, 0.25},
        {3, 0.25}};
    storage.assign(std::move(offsets), std::move(invalid_cell), 1, 100);
    REQUIRE_FALSE(storage.is_structured());
}

TEST_CASE("Endpoint storage rejects inconsistent offsets",
          "[successive_orders][endpoint_storage]") {
    EndpointStencilStorage storage;
    REQUIRE_THROWS_AS(storage.assign({0, 5}, {{0, 1.0}}, 10, 100),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(storage.assign({0, 1, 0}, {}, 10, 100),
                      std::invalid_argument);
}
