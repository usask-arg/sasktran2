#include "../../successive_orders/endpoint_stencil_storage.h"

#include <sasktran2/test_helper.h>

#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
    using namespace sasktran2::successive_orders;

    class EndpointTestRounding {
      public:
        EndpointTestRounding() : m_original(std::fegetround()) {
            REQUIRE(m_original != -1);
            REQUIRE(std::fesetround(FE_TONEAREST) == 0);
        }
        ~EndpointTestRounding() { std::fesetround(m_original); }

      private:
        int m_original;
    };

    std::uint64_t endpoint_factor_bits(double value) {
        std::uint64_t result;
        std::memcpy(&result, &value, sizeof(result));
        return result;
    }
    double endpoint_factor_from_bits(std::uint64_t value) {
        double result;
        std::memcpy(&result, &value, sizeof(result));
        return result;
    }

    // An independent copy of the tracer's original six statements produces
    // the saved coefficients; the fixture does not use the candidate decoder.
    std::array<double, 4>
    original_endpoint_products(const std::array<double, 2>& factors) {
        const double altitude_upper = factors[0];
        const double horizontal_upper = factors[1];
        const double altitude_lower = 1.0 - altitude_upper;
        const double horizontal_lower = 1.0 - horizontal_upper;
        std::array<double, 4> result;
        result[0] = horizontal_lower * altitude_lower;
        result[1] = horizontal_lower * altitude_upper;
        result[2] = horizontal_upper * altitude_lower;
        result[3] = horizontal_upper * altitude_upper;
        return result;
    }

    struct EndpointOriginals {
        std::vector<int> offsets{0};
        std::vector<InterpolationWeight> weights;
        std::vector<std::array<double, 2>> factors;
    };

    EndpointOriginals
    endpoint_originals(const std::vector<std::array<double, 2>>& factors,
                       int stride, const std::vector<int>& bases) {
        REQUIRE(factors.size() == bases.size());
        EndpointOriginals result;
        result.factors = factors;
        for (std::size_t slot = 0; slot < factors.size(); ++slot) {
            const auto weights = original_endpoint_products(factors[slot]);
            for (int corner = 0; corner < 4; ++corner) {
                result.weights.emplace_back(bases[slot] + (corner & 1) +
                                                (corner >> 1) * stride,
                                            weights[corner]);
            }
            result.offsets.push_back(static_cast<int>(result.weights.size()));
        }
        return result;
    }

    EndpointStencilStorage
    assign_endpoint_originals(const EndpointOriginals& originals, int stride,
                              int geometry_size, bool include_factors = true) {
        EndpointStencilStorage result;
        auto offsets = originals.offsets;
        auto weights = originals.weights;
        auto factors = include_factors ? originals.factors
                                       : std::vector<std::array<double, 2>>{};
        result.assign(std::move(offsets), std::move(weights), stride,
                      geometry_size, std::move(factors));
        return result;
    }

    template <typename Entry>
    std::pair<int, double> endpoint_factor_pair(const Entry& entry) {
        if constexpr (std::is_same_v<Entry, InterpolationWeight>) {
            return {entry.index, entry.weight()};
        } else {
            return entry;
        }
    }

    template <typename View>
    void require_endpoint_view(const View& view, const EndpointOriginals& saved,
                               std::size_t slot) {
        const auto begin = saved.offsets[slot];
        const auto end = saved.offsets[slot + 1];
        REQUIRE(view.size() == static_cast<std::size_t>(end - begin));
        for (std::size_t index = 0; index < view.size(); ++index) {
            const auto [location, weight] = endpoint_factor_pair(view[index]);
            REQUIRE(location == saved.weights[begin + index].index);
            REQUIRE(
                endpoint_factor_bits(weight) ==
                endpoint_factor_bits(saved.weights[begin + index].weight()));
        }
        REQUIRE_THROWS_AS(view[view.size()], std::out_of_range);
    }

    void require_endpoint_originals(const EndpointStencilStorage& storage,
                                    const EndpointOriginals& originals) {
        REQUIRE(storage.size() == originals.offsets.size() - 1);
        REQUIRE(storage.weight_count() == originals.weights.size());
        storage.visit([&](const auto& stencils) {
            for (std::size_t slot = 0; slot < stencils.size(); ++slot) {
                require_endpoint_view(stencils[slot], originals, slot);
            }
        });
        for (std::size_t slot = 0; slot < storage.size(); ++slot) {
            storage.view(slot).visit([&](const auto& weights) {
                require_endpoint_view(weights, originals, slot);
            });
        }
        REQUIRE_THROWS_AS(storage.view(storage.size()), std::out_of_range);
    }

    struct EndpointProducts {
        std::vector<double> extinction;
        std::vector<double> albedo;
        std::vector<double> tangent;
        std::array<double, 16> extinction_gradient{};
        std::array<double, 16> albedo_gradient{};
    };

    EndpointProducts endpoint_products(const EndpointStencilStorage& storage) {
        EndpointProducts result;
        storage.visit([&](const auto& stencils) {
            for (std::size_t slot = 0; slot < stencils.size(); ++slot) {
                const auto weights = stencils[slot];
                double extinction = 0.0;
                double albedo = 0.0;
                double tangent = 0.0;
                const double cotangent = 0.039 + 0.0027 * slot;
                for (std::size_t index = 0; index < weights.size(); ++index) {
                    const auto [location, weight] =
                        endpoint_factor_pair(weights[index]);
                    extinction +=
                        weight * std::ldexp(1.019 + 0.173 * location,
                                            location % 2 == 0 ? 400 : -400);
                    albedo += weight * (0.117 - 0.13 * location);
                    tangent += weight * (-0.0037 + 0.0013 * location);
                    if (weight != 0.0) {
                        result.extinction_gradient[location] +=
                            weight * cotangent;
                        result.albedo_gradient[location] +=
                            weight * cotangent * (0.71 + 0.03 * location);
                    }
                }
                result.extinction.push_back(extinction);
                result.albedo.push_back(albedo);
                result.tangent.push_back(tangent);
            }
        });
        return result;
    }

    void require_endpoint_product_bits(const EndpointProducts& actual,
                                       const EndpointProducts& expected) {
        REQUIRE(actual.extinction.size() == expected.extinction.size());
        for (std::size_t index = 0; index < actual.extinction.size(); ++index) {
            REQUIRE(endpoint_factor_bits(actual.extinction[index]) ==
                    endpoint_factor_bits(expected.extinction[index]));
            REQUIRE(endpoint_factor_bits(actual.albedo[index]) ==
                    endpoint_factor_bits(expected.albedo[index]));
            REQUIRE(endpoint_factor_bits(actual.tangent[index]) ==
                    endpoint_factor_bits(expected.tangent[index]));
        }
        for (std::size_t index = 0; index < actual.extinction_gradient.size();
             ++index) {
            REQUIRE(endpoint_factor_bits(actual.extinction_gradient[index]) ==
                    endpoint_factor_bits(expected.extinction_gradient[index]));
            REQUIRE(endpoint_factor_bits(actual.albedo_gradient[index]) ==
                    endpoint_factor_bits(expected.albedo_gradient[index]));
        }
    }
} // namespace

TEST_CASE("Factored endpoint stencils preserve original IEEE coefficients at "
          "fraction boundaries and all typed consumer products",
          "[successive_orders][endpoint_factors]") {
    EndpointTestRounding rounding;
    const std::array<double, 8> boundaries{
        0.0,
        -0.0,
        1.0,
        std::nextafter(1.0, 0.0),
        std::numeric_limits<double>::denorm_min(),
        std::numeric_limits<double>::min(),
        0.5,
        std::nextafter(0.5, 0.0)};
    std::vector<std::array<double, 2>> factors;
    std::vector<int> bases;
    for (double altitude : boundaries) {
        for (double horizontal : boundaries) {
            factors.push_back({altitude, horizontal});
            bases.push_back(2 * static_cast<int>(bases.size() % 4));
        }
    }
    const auto originals = endpoint_originals(factors, 2, bases);
    auto factored = assign_endpoint_originals(originals, 2, 12);
    const auto independent = assign_endpoint_originals(originals, 2, 12, false);
    REQUIRE(factored.is_factored());
    REQUIRE(factored.is_structured());
    REQUIRE(factored.base_bits() == 16);
    REQUIRE(factored.storage_bytes() == factors.size() * 18);
    REQUIRE(independent.storage_bytes() == factors.size() * 34);
    factored.prepare_for_current_rounding();
    REQUIRE(factored.is_factored());
    require_endpoint_originals(factored, originals);
    require_endpoint_product_bits(endpoint_products(factored),
                                  endpoint_products(independent));
}

TEST_CASE("Factored endpoint stencils select full sixteen or thirty-two bit "
          "bases without truncating corner indices",
          "[successive_orders][endpoint_factors]") {
    EndpointTestRounding rounding;
    for (int base : {65535, 65536}) {
        const auto originals =
            endpoint_originals({{0.13159, 0.271828}}, 3, {base});
        const auto storage = assign_endpoint_originals(originals, 3, 65544);
        REQUIRE(storage.is_factored());
        REQUIRE(storage.base_bits() == (base == 65535 ? 16 : 32));
        REQUIRE(storage.storage_bytes() == (base == 65535 ? 18 : 20));
        require_endpoint_originals(storage, originals);
    }
}

TEST_CASE("Factored endpoint adoption rejects a whole provider when capture "
          "provenance or any original coefficient differs",
          "[successive_orders][endpoint_factors]") {
    EndpointTestRounding rounding;
    const auto originals =
        endpoint_originals({{0.37, 0.41}, {0.13, 0.27}}, 2, {0, 4});
    const auto check_fallback = [&](EndpointOriginals altered) {
        const auto storage = assign_endpoint_originals(altered, 2, 12);
        REQUIRE_FALSE(storage.is_factored());
        REQUIRE(storage.is_structured());
        REQUIRE(storage.base_bits() == 16);
        require_endpoint_originals(storage, altered);
    };
    auto altered = originals;
    const double old = altered.weights[5].weight();
    altered.weights[5] = {altered.weights[5].index, std::nextafter(old, 1.0)};
    check_fallback(altered);
    altered = originals;
    altered.factors.clear();
    check_fallback(altered);
    altered = originals;
    altered.factors.pop_back();
    check_fallback(altered);
    altered = originals;
    std::swap(altered.factors[0], altered.factors[1]);
    check_fallback(altered);
    for (double invalid : {-0.01, std::nextafter(1.0, 2.0),
                           std::numeric_limits<double>::infinity(),
                           endpoint_factor_from_bits(0x7ff8123456789abcULL)}) {
        altered = originals;
        altered.factors[1][0] = invalid;
        check_fallback(altered);
    }
    altered = originals;
    altered.weights[0] = {0, endpoint_factor_from_bits(0x7ff8123456789abcULL)};
    altered.weights[1] = {1, -std::numeric_limits<double>::infinity()};
    altered.weights[2] = {2, -0.0};
    check_fallback(altered);
}

TEST_CASE("Endpoint factors preserve generic and uncaptured native inputs",
          "[successive_orders][endpoint_factors]") {
    EndpointTestRounding rounding;
    auto originals = endpoint_originals({{0.37, 0.41}}, 2, {0});
    originals.weights.erase(originals.weights.begin() + 1);
    originals.offsets.back() = 3;
    const auto generic = assign_endpoint_originals(originals, 2, 12);
    REQUIRE_FALSE(generic.is_structured());
    REQUIRE_FALSE(generic.is_factored());
    require_endpoint_originals(generic, originals);

    originals = endpoint_originals({{0.37, 0.41}}, 2, {0});
    std::swap(originals.weights[1], originals.weights[2]);
    const auto reordered = assign_endpoint_originals(originals, 2, 12);
    REQUIRE_FALSE(reordered.is_structured());
    require_endpoint_originals(reordered, originals);

    originals = endpoint_originals({{0.37, 0.41}}, 2, {0});
    const auto uncaptured = assign_endpoint_originals(originals, 2, 12, false);
    REQUIRE(uncaptured.is_structured());
    REQUIRE_FALSE(uncaptured.is_factored());
    require_endpoint_originals(uncaptured, originals);
    EndpointStencilStorage empty;
    empty.assign({}, {}, 2, 12, {});
    REQUIRE(empty.empty());
    REQUIRE_FALSE(empty.is_factored());
    REQUIRE_THROWS_AS(empty.assign({}, {{0, 1.0}}, 2, 12, {}),
                      std::invalid_argument);
}

TEST_CASE("Factored endpoint views own their expanded coefficients through "
          "temporary destruction, moves and storage reset",
          "[successive_orders][endpoint_factors]") {
    EndpointTestRounding rounding;
    const auto originals =
        endpoint_originals({{0.37, 0.41}, {0.13, 0.27}}, 2, {0, 4});
    auto storage = assign_endpoint_originals(originals, 2, 12);
    auto copied = storage;
    auto moved = std::move(copied);
    require_endpoint_originals(moved, originals);
    const auto held = storage.view(1);
    auto held_copy = held;
    auto held_move = std::move(held_copy);
    storage.clear();
    moved.clear();
    held.visit([&](const auto& weights) {
        require_endpoint_view(weights, originals, 1);
    });
    held_move.visit([&](const auto& weights) {
        require_endpoint_view(weights, originals, 1);
    });
    REQUIRE_THROWS_AS(held.visit([](const auto&) {
        throw std::runtime_error("Visitor fixture");
    }),
                      std::runtime_error);
    storage = assign_endpoint_originals(originals, 2, 12, false);
    REQUIRE_FALSE(storage.is_factored());
    storage = assign_endpoint_originals(originals, 2, 12);
    REQUIRE(storage.is_factored());
    require_endpoint_originals(storage, originals);
}

TEST_CASE("Changed rounding permanently materializes original endpoint "
          "coefficient bits and restores the caller's mode before products",
          "[successive_orders][endpoint_factors]") {
    EndpointTestRounding rounding;
    for (int base : {0, 65536}) {
        for (int changed_mode : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
            REQUIRE(std::fesetround(FE_TONEAREST) == 0);
            const auto originals = endpoint_originals(
                {{0.137, 0.271}, {0.731, 0.197}}, 3, {base, base + 3});
            auto storage =
                assign_endpoint_originals(originals, 3, base == 0 ? 18 : 65550);
            const auto independent = assign_endpoint_originals(
                originals, 3, base == 0 ? 18 : 65550, false);
            REQUIRE(storage.is_factored());
            const int base_bits = storage.base_bits();
            REQUIRE(std::fesetround(changed_mode) == 0);
            storage.prepare_for_current_rounding();
            REQUIRE(std::fegetround() == changed_mode);
            REQUIRE_FALSE(storage.is_factored());
            REQUIRE(storage.is_structured());
            REQUIRE(storage.base_bits() == base_bits);
            REQUIRE(storage.storage_bytes() ==
                    originals.factors.size() * (base_bits == 16 ? 34 : 36));
            require_endpoint_originals(storage, originals);
            if (base == 0) {
                require_endpoint_product_bits(endpoint_products(storage),
                                              endpoint_products(independent));
            }
            REQUIRE(std::fesetround(FE_TONEAREST) == 0);
            storage.prepare_for_current_rounding();
            REQUIRE_FALSE(storage.is_factored());
            require_endpoint_originals(storage, originals);
        }
    }
}

TEST_CASE("Typed endpoint callbacks preserve every representation and reject "
          "invalid slots before invoking the consumer",
          "[successive_orders][endpoint_factors][visit_view]") {
    EndpointTestRounding rounding;
    // Generic, structured16/32, and factored16/32 use the same two stencils.
    const int representation = GENERATE(0, 1, 2, 3, 4);
    const bool wide = representation == 2 || representation == 4;
    const bool factored = representation >= 3;
    const int base = wide ? 65536 : 0;
    auto originals =
        endpoint_originals({{0.137, 0.271}, {-0.0, 1.0}}, 3, {base, base + 3});
    if (representation == 0) {
        originals.weights.erase(originals.weights.begin() + 5);
        originals.weights.erase(originals.weights.begin() + 1);
        originals.offsets = {0, 3, 6};
    }
    auto storage =
        assign_endpoint_originals(originals, 3, wide ? 65550 : 18, factored);
    REQUIRE(storage.is_factored() == factored);
    REQUIRE(storage.base_bits() == (representation == 0 ? 0 : wide ? 32 : 16));
    std::size_t calls = 0;
    for (std::size_t slot = 0; slot < storage.size(); ++slot) {
        const auto returned =
            storage.visit_view(slot, [&](const auto& weights) {
                ++calls;
                using View = std::decay_t<decltype(weights)>;
                REQUIRE((std::is_same_v<View, FactoredEndpointStencilView>) ==
                        factored);
                REQUIRE(
                    (std::is_same_v<View,
                                    InterpolationView<InterpolationWeight>>) ==
                    (representation == 0));
                require_endpoint_view(weights, originals, slot);
                const auto legacy = storage.view(slot);
                legacy.visit([&](const auto& expected) {
                    REQUIRE(weights.size() == expected.size());
                    for (std::size_t index = 0; index < weights.size();
                         ++index) {
                        const auto actual_entry =
                            endpoint_factor_pair(weights[index]);
                        const auto expected_entry =
                            endpoint_factor_pair(expected[index]);
                        REQUIRE(actual_entry.first == expected_entry.first);
                        REQUIRE(endpoint_factor_bits(actual_entry.second) ==
                                endpoint_factor_bits(expected_entry.second));
                    }
                });
                return weights.size();
            });
        REQUIRE(returned ==
                static_cast<std::size_t>(originals.offsets[slot + 1] -
                                         originals.offsets[slot]));
    }
    REQUIRE(calls == storage.size());
    const auto rejected = [&](const auto&) {
        ++calls;
        return 0;
    };
    REQUIRE_THROWS_AS(storage.visit_view(storage.size(), rejected),
                      std::out_of_range);
    REQUIRE_THROWS_AS(
        storage.visit_view(std::numeric_limits<std::size_t>::max(), rejected),
        std::out_of_range);
    REQUIRE(calls == storage.size());
    REQUIRE_THROWS_AS(storage.visit_view(0,
                                         [](const auto&) {
                                             throw std::runtime_error(
                                                 "Typed visitor fixture");
                                         }),
                      std::runtime_error);
    storage.clear();
    REQUIRE_THROWS_AS(storage.visit_view(0, rejected), std::out_of_range);
    REQUIRE(calls == 2);
    storage =
        assign_endpoint_originals(originals, 3, wide ? 65550 : 18, factored);
    storage.visit_view(1, [&](const auto& weights) {
        require_endpoint_view(weights, originals, 1);
    });
}

TEST_CASE("Factored typed callback results own coefficients through clearing "
          "and reassigning their storage",
          "[successive_orders][endpoint_factors][visit_view]") {
    EndpointTestRounding rounding;
    const int base = GENERATE(0, 65536);
    const auto originals = endpoint_originals({{0.137, 0.271}, {0.731, 0.197}},
                                              3, {base, base + 3});
    auto storage =
        assign_endpoint_originals(originals, 3, base == 0 ? 18 : 65550);
    REQUIRE(storage.is_factored());
    const auto returned = storage.visit_view(1, [&](const auto& weights) {
        REQUIRE((std::is_same_v<std::decay_t<decltype(weights)>,
                                FactoredEndpointStencilView>));
        const auto owned = EndpointStencilView(weights);
        storage.clear();
        // The concrete callback view remains valid after its owner is cleared.
        require_endpoint_view(weights, originals, 1);
        return owned;
    });
    REQUIRE(storage.empty());
    storage =
        assign_endpoint_originals(originals, 3, base == 0 ? 18 : 65550, false);
    REQUIRE_FALSE(storage.is_factored());
    auto copied = returned;
    auto moved = std::move(copied);
    returned.visit([&](const auto& weights) {
        require_endpoint_view(weights, originals, 1);
    });
    moved.visit([&](const auto& weights) {
        require_endpoint_view(weights, originals, 1);
    });
}
