#include "../../successive_orders/source_weight_storage.h"

#include <sasktran2/test_helper.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
    using namespace sasktran2::successive_orders;

    // A minimal wide record lets the fixture exercise arbitrary IEEE payloads.
    struct CodecWideWeight {
        CodecWideWeight(int index, double weight)
            : m_index(static_cast<std::uint32_t>(index)) {
            std::memcpy(m_weight.data(), &weight, sizeof(weight));
        }
        std::uint32_t row_inner_index() const { return m_index; }
        double weight() const {
            double result;
            std::memcpy(&result, m_weight.data(), sizeof(result));
            return result;
        }

      private:
        std::uint32_t m_index;
        std::array<std::byte, sizeof(double)> m_weight;
    };
    static_assert(sizeof(CodecWideWeight) == 12);
    using CodecStorage = SourceWeightStorage<CodecWideWeight>;

    std::uint64_t codec_bits(double value) {
        std::uint64_t result;
        std::memcpy(&result, &value, sizeof(result));
        return result;
    }

    double codec_from_bits(std::uint64_t value) {
        double result;
        std::memcpy(&result, &value, sizeof(result));
        return result;
    }

    std::uint64_t positive_payload(unsigned exponent, std::uint64_t mantissa) {
        return (static_cast<std::uint64_t>(exponent) << 52) |
               (mantissa & 0x000fffffffffffffULL);
    }

    std::vector<CodecWideWeight> boundary_payloads(unsigned base) {
        std::vector<CodecWideWeight> result;
        result.reserve(544);
        // Both ends of a 15-exponent window force its base. The mantissas are
        // different full-bit payloads, not rounded decimal literals.
        for (int index = 0; index < 512; ++index) {
            const unsigned exponent = base + (index % 2 == 0 ? 0 : 14);
            const auto mantissa =
                index % 2 == 0
                    ? static_cast<std::uint64_t>(index + 1)
                    : 0x000fffffffffffffULL - static_cast<std::uint64_t>(index);
            result.emplace_back(
                index % 2 == 0 ? 0 : 255,
                codec_from_bits(positive_payload(exponent, mantissa)));
        }
        const std::array<std::uint64_t, 17> special{
            0x0000000000000000ULL,
            0x8000000000000000ULL,
            0x0000000000000001ULL,
            0x000fffffffffffffULL,
            0x0010000000000000ULL,
            0x8000000000000001ULL,
            0xbfd5555555555555ULL,
            0x7fefffffffffffffULL,
            0x7ff0000000000000ULL,
            0xfff0000000000000ULL,
            0x7ff0000000000001ULL,
            0x7ff8123456789abcULL,
            0xfff8123456789abcULL,
            positive_payload(base, 0x00076543210fedcbULL),
            positive_payload(base + 14, 0x000123456789abcdULL),
            positive_payload(base + 15, 0x0009876543210abcULL),
            positive_payload(base == 0 ? 31 : base - 1, 0x000abcdef0123456ULL)};
        for (std::size_t index = 0; index < special.size(); ++index) {
            result.emplace_back(static_cast<int>(index % 2 == 0 ? 255 : 0),
                                codec_from_bits(special[index]));
            REQUIRE(codec_bits(result.back().weight()) == special[index]);
        }
        return result;
    }

    CodecStorage storage_from(const std::vector<CodecWideWeight>& values) {
        CodecStorage result;
        result.wide_values() = values;
        return result;
    }

    void
    require_original_records(const CodecStorage& storage,
                             const std::vector<CodecWideWeight>& original) {
        REQUIRE(storage.size() == original.size());
        storage.visit([&](const auto& values) {
            REQUIRE(values.size() == original.size());
            std::size_t index = 0;
            for (const auto& value : values) {
                REQUIRE(value.row_inner_index() ==
                        original[index].row_inner_index());
                REQUIRE(codec_bits(value.weight()) ==
                        codec_bits(original[index].weight()));
                ++index;
            }
            REQUIRE(index == original.size());
        });
        std::size_t index = 0;
        for (const auto value : storage) {
            REQUIRE(value.row_inner_index() ==
                    original[index].row_inner_index());
            REQUIRE(codec_bits(value.weight()) ==
                    codec_bits(original[index].weight()));
            ++index;
        }
        REQUIRE(index == original.size());
        if (!original.empty()) {
            const auto last = storage.view(original.size() - 1, 1);
            REQUIRE(codec_bits(last[0].weight()) ==
                    codec_bits(original.back().weight()));
            REQUIRE(last[0].row_inner_index() ==
                    original.back().row_inner_index());
        }
        REQUIRE(storage.view(storage.size(), 0).empty());
        REQUIRE(storage.view(storage.size(), 0).begin() ==
                storage.view(storage.size(), 0).end());
        REQUIRE_THROWS_AS(storage.view(storage.size() + 1, 0),
                          std::out_of_range);
        REQUIRE_THROWS_AS(storage.view(0, storage.size() + 1),
                          std::out_of_range);
        REQUIRE_THROWS_AS(storage[storage.size()], std::out_of_range);
        REQUIRE_THROWS_AS(*storage.end(), std::out_of_range);
    }

    struct CodecProducts {
        std::array<double, 256> transport{};
        std::array<double, 256> transport_jvp{};
        std::array<double, 256> state_vjp{};
        double projected_state = 0.0;
        double projected_state_jvp = 0.0;
        double factor_vjp = 0.0;
    };

    template <typename Weights>
    void accumulate_products(const Weights& sources, int layer,
                             CodecProducts& result) {
        const double factor = 0.03125 * (layer + 1) + 0.007;
        const double factor_tangent = -0.004 * layer + 0.002;
        const double forcing_gradient = 0.37 - 0.021 * layer;
        for (const auto& source : sources) {
            const auto slot = source.row_inner_index();
            const double state = 0.014 * slot - 1.19;
            const double state_tangent = -0.003 * slot + 0.023;
            const double value_gradient = 0.019 * slot - 1.7;
            // Keep the source gather and transport/linearization update
            // expressions and their iteration order identical in both paths.
            result.transport[slot] += source.weight() * factor;
            result.transport_jvp[slot] += source.weight() * factor_tangent;
            result.projected_state += source.weight() * state;
            result.projected_state_jvp += source.weight() * state_tangent;
            result.factor_vjp +=
                source.weight() * forcing_gradient * value_gradient;
            result.state_vjp[slot] +=
                source.weight() * factor * forcing_gradient;
        }
    }

    void require_product_bits(const CodecProducts& actual,
                              const CodecProducts& expected) {
        for (std::size_t slot = 0; slot < actual.transport.size(); ++slot) {
            REQUIRE(codec_bits(actual.transport[slot]) ==
                    codec_bits(expected.transport[slot]));
            REQUIRE(codec_bits(actual.transport_jvp[slot]) ==
                    codec_bits(expected.transport_jvp[slot]));
            REQUIRE(codec_bits(actual.state_vjp[slot]) ==
                    codec_bits(expected.state_vjp[slot]));
        }
        REQUIRE(codec_bits(actual.projected_state) ==
                codec_bits(expected.projected_state));
        REQUIRE(codec_bits(actual.projected_state_jvp) ==
                codec_bits(expected.projected_state_jvp));
        REQUIRE(codec_bits(actual.factor_vjp) ==
                codec_bits(expected.factor_vjp));
    }
} // namespace

TEST_CASE("Seven-byte source weight codec preserves exact payloads and slots",
          "[successive_orders][source_weight_codec]") {
    for (const unsigned base : {0U, 1010U, 2032U}) {
        CAPTURE(base);
        const auto original = boundary_payloads(base);
        auto storage = storage_from(original);
        storage.narrow(256);
        REQUIRE(storage.is_encoded());
        REQUIRE(storage.is_compact());
        REQUIRE(storage.view().is_encoded());
        REQUIRE(storage.view().index_bytes() == 1);
        REQUIRE(storage.encoded_base_exponent() == base);
        REQUIRE(storage.encoded_record_count() == original.size());
        REQUIRE(storage.element_bytes() == 8);
        REQUIRE(storage.capacity_bytes() == storage.encoded_storage_bytes());
        REQUIRE(storage.encoded_storage_bytes() +
                    CodecStorage::encoded_header_growth_bytes() <
                original.size() * 9);
        REQUIRE(storage.encoded_escape_count() > 0);
        std::size_t expected_escapes = 0;
        for (const auto& value : original) {
            const auto representation = codec_bits(value.weight());
            const auto exponent = (representation >> 52) & 2047;
            const bool encoded = (representation >> 63) == 0 &&
                                 representation != 0 && exponent != 2047 &&
                                 exponent >= base && exponent < base + 15;
            expected_escapes += !encoded;
        }
        REQUIRE(storage.encoded_escape_count() == expected_escapes);
        require_original_records(storage, original);
        const auto bytes = storage.capacity_bytes();
        const auto escapes = storage.encoded_escape_count();
        storage.narrow(256);
        storage.shrink_to_fit();
        REQUIRE(storage.is_encoded());
        REQUIRE(storage.capacity_bytes() <= bytes);
        REQUIRE(storage.encoded_escape_count() == escapes);
        require_original_records(storage, original);

        // A view is a value containing its decode context, not a reference
        // to a temporary owner view created by storage.view().
        const auto saved_view = storage.view(3, 19);
        saved_view.visit([&](const auto& values) {
            for (std::size_t index = 0; index < values.size(); ++index) {
                REQUIRE(codec_bits(values[index].weight()) ==
                        codec_bits(original[index + 3].weight()));
            }
        });
        const auto saved_proxy =
            storage.view(original.size() - 1, 1)
                .visit([&](const auto& values)
                           -> std::optional<EncodedSourceWeight> {
                    if constexpr (std::is_same_v<decltype(values[0]),
                                                 EncodedSourceWeight>) {
                        return values[0];
                    } else {
                        return std::nullopt;
                    }
                });
        REQUIRE(saved_proxy.has_value());
        REQUIRE(codec_bits(saved_proxy->weight()) ==
                codec_bits(original.back().weight()));
        REQUIRE(saved_proxy->row_inner_index() ==
                original.back().row_inner_index());
        auto copied = storage;
        auto moved = std::move(copied);
        REQUIRE(copied.empty());
        REQUIRE(copied.encoded_record_count() == 0);
        REQUIRE(copied.encoded_escape_count() == 0);
        REQUIRE(copied.view().begin() == copied.view().end());
        require_original_records(moved, original);
        copied = {{0, 0.375}};
        copied.narrow(256);
        REQUIRE_FALSE(copied.is_encoded());
        REQUIRE(codec_bits(copied[0].weight()) == codec_bits(0.375));
    }
}

TEST_CASE("Source codec profitability fallback and reset preserve records",
          "[successive_orders][source_weight_codec]") {
    const std::array<std::vector<CodecWideWeight>, 5> cases{
        std::vector<CodecWideWeight>{},
        std::vector<CodecWideWeight>(512, CodecWideWeight(255, +0.0)),
        std::vector<CodecWideWeight>(512, CodecWideWeight(0, -0.375)),
        std::vector<CodecWideWeight>(8, CodecWideWeight(255, 0.125)),
        std::vector<CodecWideWeight>{{0, 0.125}, {255, 0.875}}};
    for (const auto& original : cases) {
        CAPTURE(original.size());
        auto storage = storage_from(original);
        storage.narrow(256);
        REQUIRE_FALSE(storage.is_encoded());
        REQUIRE(storage.element_bytes() == 9);
        REQUIRE(storage.view().index_bytes() == 1);
        require_original_records(storage, original);
    }
    EncodedSourceWeightStorage empty_encoded;
    const SourceWeightView<CodecWideWeight> empty_encoded_view(empty_encoded, 0,
                                                               0);
    REQUIRE(empty_encoded_view.empty());
    REQUIRE(empty_encoded_view.begin() == empty_encoded_view.end());
    empty_encoded_view.visit([&](const auto& values) {
        REQUIRE(values.empty());
        REQUIRE(values.begin() == values.end());
    });
    EncodedSourceWeightStorage invalid_encoded;
    invalid_encoded.record_count = 1;
    REQUIRE_THROWS_AS(SourceWeightView<CodecWideWeight>(invalid_encoded, 0, 0),
                      std::out_of_range);
    invalid_encoded.words.resize(1);
    invalid_encoded.base_exponent = 2033;
    REQUIRE_THROWS_AS(SourceWeightView<CodecWideWeight>(invalid_encoded, 0, 1),
                      std::out_of_range);
    const std::vector<CodecWideWeight> no_escape_original(
        512, CodecWideWeight(255, 0.12345678901234567));
    auto no_escape = storage_from(no_escape_original);
    no_escape.narrow(256);
    REQUIRE(no_escape.is_encoded());
    REQUIRE(no_escape.encoded_escape_count() == 0);
    require_original_records(no_escape, no_escape_original);
    const auto original = boundary_payloads(1010);
    for (const unsigned row_size : {257U, 65536U, 65537U}) {
        CAPTURE(row_size);
        auto storage = storage_from(original);
        storage.narrow(row_size);
        REQUIRE_FALSE(storage.is_encoded());
        REQUIRE(storage.element_bytes() == (row_size <= 65536 ? 10 : 12));
        require_original_records(storage, original);
    }
    auto storage = storage_from(original);
    storage.narrow(256);
    REQUIRE(storage.is_encoded());
    storage = {{255, -0.0}, {0, 0.25}};
    REQUIRE_FALSE(storage.is_encoded());
    REQUIRE(storage.wide_values().size() == 2);
    REQUIRE(codec_bits(storage[0].weight()) == 0x8000000000000000ULL);
    storage.narrow(256);
    REQUIRE_FALSE(storage.is_encoded());
    storage = {};
    REQUIRE_FALSE(storage.is_encoded());
    REQUIRE(storage.empty());
    REQUIRE(storage.view().begin() == storage.view().end());
    storage.wide_values() = original;
    storage.narrow(256);
    REQUIRE(storage.is_encoded());
    require_original_records(storage, original);

    auto invalid = storage_from(original);
    invalid.wide_values().emplace_back(256, 0.75);
    REQUIRE_THROWS_AS(invalid.narrow(256), std::out_of_range);
    REQUIRE_FALSE(invalid.is_encoded());
    REQUIRE_FALSE(invalid.is_compact());
    REQUIRE(invalid.wide_values().size() == original.size() + 1);
    REQUIRE(invalid.wide_values().back().row_inner_index() == 256);
}

TEST_CASE("Encoded source visitors preserve product and reduction order bits",
          "[successive_orders][source_weight_codec]"
          "[linearization]") {
    std::vector<CodecWideWeight> original;
    for (int index = 0; index < 768; ++index) {
        const unsigned exponent = 1010 + static_cast<unsigned>(index % 15);
        const auto mantissa =
            0x0000123456789abcULL +
            static_cast<std::uint64_t>(index) * 0x0000000789abcdefULL;
        original.emplace_back(
            (index * 37) % 256,
            codec_from_bits(positive_payload(exponent, mantissa)));
    }
    // Sparse escapes also pass through arithmetic and shared-slot reductions.
    original.emplace_back(255, -0.3123456789012345);
    original.emplace_back(0, -0.0);
    original.emplace_back(255, 0.000000000000001);
    auto storage = storage_from(original);
    storage.narrow(256);
    REQUIRE(storage.is_encoded());
    REQUIRE(storage.encoded_escape_count() >= 3);
    CodecProducts expected;
    CodecProducts actual;
    std::size_t begin = 0;
    for (int layer = 0; begin < original.size(); ++layer) {
        const std::size_t count =
            std::min<std::size_t>(97, original.size() - begin);
        const auto wide = SourceWeightArrayView<CodecWideWeight>(
            original.data() + begin, count);
        accumulate_products(wide, layer, expected);
        storage.view(begin, count).visit([&](const auto& values) {
            accumulate_products(values, layer, actual);
        });
        begin += count;
    }
    require_product_bits(actual, expected);
}
