#include "../../successive_orders/interpolation.h"

#include <sasktran2/test_helper.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

namespace {
    using namespace sasktran2::successive_orders;

    std::uint64_t layer_weight_bits(double value) {
        std::uint64_t result;
        std::memcpy(&result, &value, sizeof(result));
        return result;
    }

    void require_structured_descriptor(
        const StructuredLayerInterpolation& actual,
        const StructuredLayerInterpolation& expected) {
        REQUIRE(actual.atmosphere_offset == expected.atmosphere_offset);
        REQUIRE(actual.source_offset == expected.source_offset);
        REQUIRE(actual.cell_base == expected.cell_base);
        REQUIRE(actual.atmosphere_mask == expected.atmosphere_mask);
        REQUIRE(actual.reserved == expected.reserved);
    }

    void require_layer_descriptor(const LayerInterpolation& actual,
                                  const LayerInterpolation& expected) {
        REQUIRE(actual.atmosphere_offset == expected.atmosphere_offset);
        REQUIRE(actual.atmosphere_count == expected.atmosphere_count);
        REQUIRE(actual.source_offset == expected.source_offset);
        REQUIRE(actual.source_count == expected.source_count);
        REQUIRE(actual.optical_depth_offset == expected.optical_depth_offset);
        REQUIRE(actual.optical_depth_count == expected.optical_depth_count);
    }

    struct LayerProducts {
        double state = 0.0;
        double tangent = 0.0;
        std::array<double, 40> extinction_gradient{};
        std::array<double, 40> albedo_gradient{};
    };

    // These products use the native layer traversal and interpolation update
    // order. The storage representation changes only how indices are read.
    LayerProducts layer_products(const RayInterpolation& ray) {
        LayerProducts result;
        for (int layer = static_cast<int>(ray.layers.size()) - 1; layer >= 0;
             --layer) {
            const auto midpoint = ray.atmosphere_for_layer(layer);
            double albedo = 0.0;
            double albedo_tangent = 0.0;
            for (const auto& weight : midpoint) {
                albedo += weight.weight() * (0.17 + 0.001 * weight.index);
                albedo_tangent +=
                    weight.weight() * (-0.007 + 0.0002 * weight.index);
            }
            const auto od = ray.optical_depth_for_layer(layer);
            double optical_depth = 0.0;
            double optical_depth_tangent = 0.0;
            for (std::size_t index = 0; index < od.size(); ++index) {
                const auto [location, weight] = od[index];
                optical_depth += weight * (0.013 + 0.0003 * location);
                optical_depth_tangent += weight * (-0.003 + 0.0001 * location);
            }
            const double factor = albedo * optical_depth;
            const double factor_tangent =
                albedo_tangent * optical_depth + albedo * optical_depth_tangent;
            ray.source_for_layer(layer).visit([&](const auto& sources) {
                for (const auto& source : sources) {
                    result.state += source.weight() * factor;
                    result.tangent += source.weight() * factor_tangent;
                }
            });
        }
        for (std::size_t layer = 0; layer < ray.layers.size(); ++layer) {
            const double cotangent = 0.033 - 0.005 * layer;
            const auto midpoint = ray.atmosphere_for_layer(layer);
            for (const auto& weight : midpoint) {
                result.albedo_gradient[weight.index] +=
                    weight.weight() * cotangent;
            }
            const auto od = ray.optical_depth_for_layer(layer);
            for (std::size_t index = 0; index < od.size(); ++index) {
                const auto [location, weight] = od[index];
                result.extinction_gradient[location] += weight * cotangent;
            }
        }
        return result;
    }

    void require_layer_product_bits(const LayerProducts& actual,
                                    const LayerProducts& expected) {
        REQUIRE(layer_weight_bits(actual.state) ==
                layer_weight_bits(expected.state));
        REQUIRE(layer_weight_bits(actual.tangent) ==
                layer_weight_bits(expected.tangent));
        for (std::size_t index = 0; index < actual.extinction_gradient.size();
             ++index) {
            REQUIRE(layer_weight_bits(actual.extinction_gradient[index]) ==
                    layer_weight_bits(expected.extinction_gradient[index]));
            REQUIRE(layer_weight_bits(actual.albedo_gradient[index]) ==
                    layer_weight_bits(expected.albedo_gradient[index]));
        }
    }
} // namespace

TEST_CASE("Successive-orders six-byte structured descriptors preserve all "
          "corner masks and full source and cell indices",
          "[successive_orders][layer_storage]") {
    for (std::uint8_t mask = 0; mask < 16; ++mask) {
        const StructuredLayerInterpolation original{4095, 65535, 65535, mask,
                                                    0};
        LayerInterpolationStorage storage;
        storage.assign_structured({original}, 65536);
        REQUIRE(storage.is_structured());
        REQUIRE(storage.is_compact_structured());
        REQUIRE(storage.element_bytes() == 6);
        REQUIRE(storage.capacity_bytes() == 6);
        require_structured_descriptor(storage.structured_layer(0), original);
        unsigned count = 0;
        for (unsigned bit = 0; bit < 4; ++bit) {
            count += (mask >> bit) & 1;
        }
        require_layer_descriptor(storage[0], {4095, count, 65535, 1, 0, 4});
    }
    const StructuredLayerInterpolation zeros{};
    require_structured_descriptor(
        CompactStructuredLayerInterpolation(zeros).expanded(), zeros);
}

TEST_CASE("Successive-orders structured descriptor overflow keeps the whole "
          "ray in its eight-byte representation",
          "[successive_orders][layer_storage]") {
    const auto check_fallback = [](StructuredLayerInterpolation boundary) {
        const std::vector<StructuredLayerInterpolation> original{
            {0, 0, 65535, 15, 0}, boundary};
        LayerInterpolationStorage storage;
        storage.assign_structured(
            std::vector<StructuredLayerInterpolation>(original), 65536);
        REQUIRE(storage.is_structured());
        REQUIRE_FALSE(storage.is_compact_structured());
        REQUIRE(storage.element_bytes() == 8);
        REQUIRE(storage.capacity_bytes() == original.size() * 8);
        for (std::size_t index = 0; index < original.size(); ++index) {
            require_structured_descriptor(storage.structured_layer(index),
                                          original[index]);
        }
        unsigned count = 0;
        for (unsigned bit = 0; bit < 4; ++bit) {
            count += (boundary.atmosphere_mask >> bit) & 1;
        }
        require_layer_descriptor(storage[1],
                                 {boundary.atmosphere_offset, count,
                                  boundary.source_offset,
                                  65536U - boundary.source_offset, 4, 4});
        REQUIRE_THROWS_AS(CompactStructuredLayerInterpolation(boundary),
                          std::out_of_range);
    };
    check_fallback({4096, 65535, 65535, 15, 0});
    check_fallback({65535, 65535, 65535, 0, 0});
    check_fallback({4095, 65535, 65535, 15, 1});
}

TEST_CASE("Successive-orders invalid structured corner masks preserve the "
          "previous decodable storage",
          "[successive_orders][layer_storage]") {
    LayerInterpolationStorage storage;
    const StructuredLayerInterpolation original{4096, 65535, 65535, 15, 0};
    storage.assign_structured({original}, 65536);
    for (const std::uint8_t mask : {16, 255}) {
        REQUIRE_THROWS_AS(
            storage.assign_structured(
                {{0, 0, 0, 0, 0}, {4095, 65535, 65535, mask, 0}}, 65535),
            std::out_of_range);
        REQUIRE(storage.size() == 1);
        REQUIRE_FALSE(storage.is_compact_structured());
        require_structured_descriptor(storage.structured_layer(0), original);
        require_layer_descriptor(storage[0], {4096, 4, 65535, 1, 0, 4});
    }
}

TEST_CASE("Successive-orders structured descriptors retain immutable values "
          "through copying, moves and representation changes",
          "[successive_orders][layer_storage]") {
    LayerInterpolationStorage storage;
    const std::vector<StructuredLayerInterpolation> original{
        {0, 0, 0, 3, 0}, {2, 0, 1, 0, 0}, {2, 65535, 65535, 15, 0}};
    storage.assign_structured(
        std::vector<StructuredLayerInterpolation>(original), 65536);
    const auto held = storage.structured_layer(2);
    auto copied = storage;
    auto moved = std::move(copied);
    REQUIRE(copied.empty());
    for (std::size_t index = 0; index < original.size(); ++index) {
        require_structured_descriptor(moved.structured_layer(index),
                                      original[index]);
    }
    storage.assign_structured({{4096, 65535, 65535, 15, 7}}, 65536);
    REQUIRE_FALSE(storage.is_compact_structured());
    require_structured_descriptor(held, original.back());
    storage = {{1, 2, 3, 4, 5, 6}};
    REQUIRE_FALSE(storage.is_structured());
    REQUIRE(storage.element_bytes() == 24);
    require_layer_descriptor(storage[0], {1, 2, 3, 4, 5, 6});
    storage.resize(2);
    REQUIRE(storage.size() == 2);
    storage.assign_structured({{0, 0, 0, 0, 0}}, 0);
    REQUIRE(storage.is_compact_structured());
    require_layer_descriptor(storage[0], {0, 0, 0, 0, 0, 4});
    storage.assign_structured({}, 0);
    REQUIRE(storage.empty());
    REQUIRE(storage.is_structured());
    REQUIRE(storage.is_compact_structured());
    storage.shrink_to_fit();
    REQUIRE(storage.capacity_bytes() == 0);
}

TEST_CASE("Successive-orders compact descriptor views preserve midpoint and "
          "optical-depth bits and complete native interpolation products",
          "[successive_orders][layer_storage]") {
    RayInterpolation compact;
    compact.structured_altitude_stride = 10;
    compact.layers.assign_structured(
        {{0, 0, 0, 5, 0}, {2, 2, 10, 15, 0}, {6, 2, 20, 0, 0}}, 4);
    compact.structured_atmosphere_weights = {
        -0.0,
        std::nextafter(0.25, 1.0),
        std::numeric_limits<double>::denorm_min(),
        0.3,
        0.4,
        0.3};
    compact.optical_depth_weights = {0.1, -0.0, 0.2, 0.7, 0.3, 0.1,
                                     0.4, 0.2,  0.1, 0.2, 0.3, 0.4};
    compact.source_weights = {{0, 0.2}, {1, 0.8}, {1, 0.3}, {2, 0.7}};

    RayInterpolation generic;
    generic.layers = {
        {0, 2, 0, 2, 0, 4}, {2, 4, 2, 0, 4, 4}, {6, 0, 2, 2, 8, 4}};
    generic.atmosphere_weights = {
        {0, compact.structured_atmosphere_weights[0]},
        {10, compact.structured_atmosphere_weights[1]},
        {10, compact.structured_atmosphere_weights[2]},
        {11, compact.structured_atmosphere_weights[3]},
        {20, compact.structured_atmosphere_weights[4]},
        {21, compact.structured_atmosphere_weights[5]}};
    generic.optical_depth_indices = {0,  1,  10, 11, 10, 11,
                                     20, 21, 20, 21, 30, 31};
    generic.optical_depth_weights = compact.optical_depth_weights;
    generic.source_weights = compact.source_weights;

    REQUIRE(compact.layers.is_compact_structured());
    REQUIRE_FALSE(generic.layers.is_structured());
    for (std::size_t layer = 0; layer < generic.layers.size(); ++layer) {
        require_layer_descriptor(compact.layers[layer], generic.layers[layer]);
        const auto actual = compact.atmosphere_for_layer(layer);
        const auto expected = generic.atmosphere_for_layer(layer);
        REQUIRE(actual.size() == expected.size());
        for (std::size_t index = 0; index < actual.size(); ++index) {
            REQUIRE(actual[index].index == expected[index].index);
            REQUIRE(layer_weight_bits(actual[index].weight()) ==
                    layer_weight_bits(expected[index].weight()));
        }
        const auto actual_od = compact.optical_depth_for_layer(layer);
        const auto expected_od = generic.optical_depth_for_layer(layer);
        for (std::size_t index = 0; index < actual_od.size(); ++index) {
            REQUIRE(actual_od[index].first == expected_od[index].first);
            REQUIRE(layer_weight_bits(actual_od[index].second) ==
                    layer_weight_bits(expected_od[index].second));
        }
    }
    auto iterator = compact.atmosphere_for_layer(0).begin();
    REQUIRE((*iterator).index == 0);
    REQUIRE(layer_weight_bits((*iterator).weight()) == layer_weight_bits(-0.0));
    ++iterator;
    REQUIRE((*iterator).index == 10);
    require_layer_product_bits(layer_products(compact),
                               layer_products(generic));
}
