#include "../../successive_orders/ray_transport.h"

#include <sasktran2/test_helper.h>

#include <cmath>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
    constexpr int num_locations = 3;
    constexpr int num_wavelengths = 2;
    constexpr int num_source_columns = 4;

    struct RayTransportFixture {
        RayTransportFixture()
            : atmosphere(sasktran2::atmosphere::AtmosphereGridStorageFull<1>(
                             num_wavelengths, num_locations, 1),
                         sasktran2::atmosphere::Surface<1>(num_wavelengths),
                         true) {
            atmosphere.storage().total_extinction.col(0) << 0.3, 0.5, 0.2;
            atmosphere.storage().total_extinction.col(1) << 0.1, 0.4, 0.7;
            atmosphere.storage().ssa.col(0) << 0.8, 0.6, 0.4;
            atmosphere.storage().ssa.col(1) << 0.35, 0.55, 0.75;

            interpolation.resize(2);
            interpolation[0].layers = {{0, 2, 0, 2, 0, 2}, {2, 2, 2, 2, 2, 2}};
            interpolation[0].atmosphere_weights = {
                {0, 0.25}, {1, 0.75}, {1, 0.4}, {2, 0.6}};
            interpolation[0].source_weights = {
                {0, 0.25}, {2, 0.75}, {1, 0.6}, {2, 0.4}};
            interpolation[0].optical_depth_indices = {0, 1, 1, 2};
            interpolation[0].optical_depth_weights = {0.7, 0.2, 0.4, 1.1};
            interpolation[0].ground_weights = {{0, 0.2}, {3, 0.8}};
            interpolation[0].ground_hit = true;
            interpolation[0].transport_value_offset = 0;
            interpolation[0].transport_row_nnz = 4;

            interpolation[1].layers = {{0, 2, 0, 2, 0, 2}};
            interpolation[1].atmosphere_weights = {{0, 0.5}, {2, 0.5}};
            interpolation[1].source_weights = {{1, 0.7}, {3, 0.3}};
            interpolation[1].optical_depth_indices = {0, 2};
            interpolation[1].optical_depth_weights = {0.3, 0.5};
            interpolation[1].transport_value_offset = 4;
            interpolation[1].transport_row_nnz = 2;

            std::vector<int> columns;
            for (auto& ray : interpolation) {
                sasktran2::successive_orders::compile_transport_row(ray,
                                                                    columns);
            }
        }

        RayTransportFixture(const RayTransportFixture&) = delete;
        RayTransportFixture& operator=(const RayTransportFixture&) = delete;

        sasktran2::successive_orders::RayTransportMap make_map() const {
            return {interpolation, num_source_columns, row_offsets,
                    column_indices};
        }

        void perturb(int wavelength, const Eigen::VectorXd& native_tangent,
                     double scale) {
            atmosphere.storage().total_extinction.col(wavelength) +=
                scale * native_tangent.head(num_locations);
            atmosphere.storage().ssa.col(wavelength) +=
                scale * native_tangent.segment(
                            atmosphere.ssa_deriv_start_index(), num_locations);
        }

        std::vector<sasktran2::successive_orders::RayInterpolation>
            interpolation;
        const std::vector<int> row_offsets{0, 4, 6};
        const std::vector<int> column_indices{0, 1, 2, 3, 1, 3};
        sasktran2::atmosphere::Atmosphere<1> atmosphere;
    };

    Eigen::VectorXd expected_values_at_wavelength_one() {
        const double layer_zero_od = 0.7 * 0.1 + 0.2 * 0.4;
        const double layer_one_od = 0.4 * 0.4 + 1.1 * 0.7;
        const double layer_zero_ssa = 0.25 * 0.35 + 0.75 * 0.55;
        const double layer_one_ssa = 0.4 * 0.55 + 0.6 * 0.75;

        const double layer_one_factor =
            layer_one_ssa * (1.0 - std::exp(-layer_one_od));
        const double layer_zero_factor = std::exp(-layer_one_od) *
                                         layer_zero_ssa *
                                         (1.0 - std::exp(-layer_zero_od));
        const double ground_factor = std::exp(-(layer_one_od + layer_zero_od));

        const double second_ray_od = 0.3 * 0.1 + 0.5 * 0.7;
        const double second_ray_ssa = 0.5 * 0.35 + 0.5 * 0.75;
        const double second_ray_factor =
            second_ray_ssa * (1.0 - std::exp(-second_ray_od));

        Eigen::VectorXd expected(6);
        expected << 0.25 * layer_zero_factor + 0.2 * ground_factor,
            0.6 * layer_one_factor,
            0.75 * layer_zero_factor + 0.4 * layer_one_factor,
            0.8 * ground_factor, 0.7 * second_ray_factor,
            0.3 * second_ray_factor;
        return expected;
    }
} // namespace

TEST_CASE(
    "Source slot storage preserves exact payloads at every width boundary",
    "[successive_orders][ray_transport][storage]") {
    using namespace sasktran2::successive_orders;
    const std::array<std::uint64_t, 8> payloads{
        0x0000000000000000ULL, 0x8000000000000000ULL, 0x0000000000000001ULL,
        0x8000000000000001ULL, 0x3fd5555555555555ULL, 0x7fefffffffffffffULL,
        0x7ff0000000000000ULL, 0x7ff8123456789abcULL};
    const auto double_bits = [](double value) {
        std::uint64_t bits;
        std::memcpy(&bits, &value, sizeof(bits));
        return bits;
    };
    for (const std::uint32_t row_size : {256U, 257U, 65536U, 65537U}) {
        CAPTURE(row_size);
        SourceInterpolationStorage storage;
        for (std::size_t index = 0; index < payloads.size(); ++index) {
            double weight;
            std::memcpy(&weight, &payloads[index], sizeof(weight));
            storage.wide_values().emplace_back(
                index + 1 == payloads.size() ? row_size - 1 : index, weight);
        }
        const auto expected = storage.wide_values();
        const auto* wide_pointer = storage.wide_values().data();
        storage.narrow(row_size);
        const std::size_t bytes = row_size <= 256     ? 9
                                  : row_size <= 65536 ? 10
                                                      : 12;
        const unsigned char index_bytes = row_size <= 256     ? 1
                                          : row_size <= 65536 ? 2
                                                              : 4;
        REQUIRE(storage.element_bytes() == bytes);
        REQUIRE(storage.is_compact() == (row_size <= 65536));
        REQUIRE(storage.capacity_bytes() == storage.capacity() * bytes);
        REQUIRE(storage.view().index_bytes() == index_bytes);
        REQUIRE(storage.size() == expected.size());
        storage.visit([&](const auto& values) {
            REQUIRE(sizeof(values[0]) == bytes);
            for (std::size_t index = 0; index < values.size(); ++index) {
                REQUIRE(values[index].row_inner_index() ==
                        expected[index].row_inner_index());
                REQUIRE(double_bits(values[index].weight()) == payloads[index]);
            }
        });
        std::size_t index = 0;
        for (const auto value : storage) {
            REQUIRE(value.row_inner_index() ==
                    expected[index].row_inner_index());
            REQUIRE(double_bits(value.weight()) == payloads[index]);
            ++index;
        }
        REQUIRE(index == storage.size());
        const auto slice = storage.view(1, 2);
        REQUIRE(slice.size() == 2);
        REQUIRE(double_bits(slice[0].weight()) == payloads[1]);
        REQUIRE(double_bits(slice[1].weight()) == payloads[2]);
        REQUIRE(storage.view(storage.size(), 0).empty());
        REQUIRE(storage.view(storage.size(), 0).begin() ==
                storage.view(storage.size(), 0).end());
        REQUIRE_THROWS_AS(storage.view(storage.size() + 1, 0),
                          std::out_of_range);
        REQUIRE_THROWS_AS(storage.view(1, storage.size()), std::out_of_range);
        REQUIRE_THROWS_AS(storage[storage.size()], std::out_of_range);
        REQUIRE_THROWS_AS(*storage.end(), std::out_of_range);
        storage.narrow(row_size);
        REQUIRE(storage.element_bytes() == bytes);
        if (storage.is_compact()) {
            REQUIRE_THROWS_AS(storage.wide_values(), std::logic_error);
            REQUIRE_THROWS_AS(storage.reserve(20), std::logic_error);
        } else {
            REQUIRE(storage.wide_values().data() == wide_pointer);
        }
    }

    REQUIRE_THROWS_AS(CompactSourceWeight<std::uint8_t>(256, 0.3),
                      std::out_of_range);
    REQUIRE_THROWS_AS(CompactSourceWeight<std::uint16_t>(65536, 0.3),
                      std::out_of_range);
    for (const std::uint32_t row_size : {256U, 65536U}) {
        SourceInterpolationStorage invalid{{0, -0.0},
                                           {static_cast<int>(row_size), 0.3}};
        REQUIRE_THROWS_AS(invalid.narrow(row_size), std::out_of_range);
        REQUIRE_FALSE(invalid.is_compact());
        REQUIRE(invalid.wide_values().size() == 2);
        REQUIRE(invalid[1].row_inner_index() == row_size);
        REQUIRE(double_bits(invalid[0].weight()) == payloads[1]);
    }
    SourceInterpolationStorage negative{{-1, 0.4}};
    REQUIRE_THROWS_AS(negative.narrow(256), std::out_of_range);
    REQUIRE(negative.wide_values()[0].source_index() == -1);
    SourceInterpolationWeight invalid_slot(0, 0.5);
    REQUIRE_THROWS_AS(
        invalid_slot.set_row_inner_index(
            static_cast<std::uint32_t>(std::numeric_limits<int>::max()) + 1U),
        std::length_error);
    SourceInterpolationStorage empty;
    empty.narrow(0);
    REQUIRE(empty.empty());
    REQUIRE(empty.view().begin() == empty.view().end());
    REQUIRE_THROWS_AS(empty[0], std::out_of_range);
}

TEST_CASE("Source slot widths preserve complete ray transport products bitwise",
          "[successive_orders][ray_transport][storage][linearization]") {
    using namespace sasktran2::successive_orders;
    const auto require_same = [](const Eigen::VectorXd& actual,
                                 const Eigen::VectorXd& expected) {
        REQUIRE(actual.size() == expected.size());
        REQUIRE(std::memcmp(actual.data(), expected.data(),
                            actual.size() * sizeof(double)) == 0);
    };
    for (const int row_size : {256, 257, 65536, 65537}) {
        CAPTURE(row_size);
        RayTransportFixture fixture;
        auto wide_rays = fixture.interpolation;
        wide_rays[0].transport_row_nnz = row_size;
        wide_rays[1].transport_value_offset = row_size;
        wide_rays[0].source_weights.wide_values()[1] = {row_size - 1, -0.0};
        wide_rays[0].source_weights.wide_values()[3] = {row_size - 1,
                                                        0.4123456789012345};
        wide_rays[0].ground_weights.wide_values()[1] = {row_size - 1,
                                                        0.8123456789012345};
        std::vector<int> columns(row_size + 2);
        std::iota(columns.begin(), columns.begin() + row_size, 0);
        columns[row_size] = 1;
        columns[row_size + 1] = 3;
        const TransportSparsity sparsity(row_size, {0, row_size, row_size + 2},
                                         columns);
        auto narrow_rays = wide_rays;
        for (auto& ray : narrow_rays) {
            compact_ray_interpolation(ray);
        }
        REQUIRE(narrow_rays[0].source_weights.element_bytes() ==
                (row_size <= 256     ? 9
                 : row_size <= 65536 ? 10
                                     : 12));
        REQUIRE(narrow_rays[0].ground_weights.element_bytes() ==
                narrow_rays[0].source_weights.element_bytes());
        const RayTransportMap wide_map(wide_rays, sparsity);
        const RayTransportMap narrow_map(narrow_rays, sparsity);
        TransportOperator wide_transport(wide_map.sparsity()),
            narrow_transport(narrow_map.sparsity());
        const Eigen::VectorXd tangent = Eigen::VectorXd::LinSpaced(
            fixture.atmosphere.num_deriv(), -0.07, 0.09);
        const Eigen::VectorXd value_gradient =
            Eigen::VectorXd::LinSpaced(sparsity.nonzeros(), -0.21, 0.31);
        for (const int wavelength : {0, 1}) {
            CAPTURE(wavelength);
            wide_map.assemble_values(fixture.atmosphere, wavelength,
                                     wide_transport);
            narrow_map.assemble_values(fixture.atmosphere, wavelength,
                                       narrow_transport);
            require_same(narrow_transport.values(), wide_transport.values());
            Eigen::VectorXd wide_jvp(sparsity.nonzeros()),
                narrow_jvp(sparsity.nonzeros());
            wide_map.assemble_jvp(fixture.atmosphere, wavelength, tangent,
                                  wide_jvp);
            narrow_map.assemble_jvp(fixture.atmosphere, wavelength, tangent,
                                    narrow_jvp);
            require_same(narrow_jvp, wide_jvp);
            Eigen::VectorXd wide_gradient =
                Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
            Eigen::VectorXd narrow_gradient = wide_gradient;
            RayTransportWorkspace wide_workspace, narrow_workspace;
            wide_map.accumulate_vjp(fixture.atmosphere, wavelength,
                                    value_gradient, wide_gradient,
                                    wide_workspace);
            narrow_map.accumulate_vjp(fixture.atmosphere, wavelength,
                                      value_gradient, narrow_gradient,
                                      narrow_workspace);
            require_same(narrow_gradient, wide_gradient);
            const Eigen::VectorXd state =
                Eigen::VectorXd::LinSpaced(row_size, -0.4, 0.6);
            Eigen::VectorXd wide_value(2), narrow_value(2);
            wide_transport.apply(state, wide_value);
            narrow_transport.apply(state, narrow_value);
            require_same(narrow_value, wide_value);
        }
        // Width conversion trusts finalized metadata; the owning map validates
        // actual row membership for every typed representation.
        auto invalid_rays = narrow_rays;
        invalid_rays[0].source_weights = {{row_size, 0.5}};
        if (row_size <= 256) {
            invalid_rays[0].source_weights.narrow(257);
        } else if (row_size < 65536) {
            invalid_rays[0].source_weights.narrow(65536);
        }
        REQUIRE_THROWS_AS(RayTransportMap(invalid_rays, sparsity),
                          std::invalid_argument);
    }
    RayTransportFixture invalid_byte;
    invalid_byte.interpolation[0].source_weights = {{4, 0.5}};
    invalid_byte.interpolation[0].source_weights.narrow(4);
    REQUIRE(invalid_byte.interpolation[0].source_weights.element_bytes() == 9);
    REQUIRE_THROWS_AS(invalid_byte.make_map(), std::invalid_argument);
}

TEST_CASE("Successive-orders packed ray transport assembles layer and ground "
          "values",
          "[successive_orders][ray_transport]") {
    RayTransportFixture fixture;
    const auto map = fixture.make_map();
    sasktran2::successive_orders::TransportOperator transport(map.sparsity());

    map.assemble_values(fixture.atmosphere, 1, transport);

    REQUIRE(map.num_rays() == 2);
    REQUIRE(map.maximum_layers() == 2);
    REQUIRE(map.sparsity().row_offsets() == fixture.row_offsets);
    REQUIRE(map.sparsity().column_indices() == fixture.column_indices);
    REQUIRE(map.sparsity().row_offsets().data() != fixture.row_offsets.data());
    REQUIRE(map.sparsity().column_indices().data() !=
            fixture.column_indices.data());
    REQUIRE(transport.values().isApprox(expected_values_at_wavelength_one(),
                                        2.0e-14));
}

TEST_CASE("Successive-orders ray transport retains its shared CSR generation",
          "[successive_orders][ray_transport][ownership]") {
    using sasktran2::successive_orders::RayTransportMap;
    using sasktran2::successive_orders::TransportOperator;
    using sasktran2::successive_orders::TransportSparsity;

    RayTransportFixture fixture;
    std::unique_ptr<RayTransportMap> map;
    const int* shared_offsets = nullptr;
    const int* shared_columns = nullptr;
    {
        TransportSparsity generation(num_source_columns, fixture.row_offsets,
                                     fixture.column_indices);
        const auto copied_generation = generation;
        shared_offsets = generation.row_offsets().data();
        shared_columns = generation.column_indices().data();
        REQUIRE(copied_generation.row_offsets().data() == shared_offsets);
        REQUIRE(copied_generation.column_indices().data() == shared_columns);

        map = std::make_unique<RayTransportMap>(fixture.interpolation,
                                                generation);
        REQUIRE(map->sparsity().row_offsets().data() == shared_offsets);
        REQUIRE(map->sparsity().column_indices().data() == shared_columns);

        generation = TransportSparsity(num_source_columns, {0, 0, 0}, {});
        REQUIRE(generation.nonzeros() == 0);
        REQUIRE(copied_generation.column_indices() == fixture.column_indices);
        REQUIRE(map->sparsity().column_indices() == fixture.column_indices);
    }

    // Both source handles are gone; the map keeps the previous generation.
    REQUIRE(map->sparsity().row_offsets().data() == shared_offsets);
    REQUIRE(map->sparsity().column_indices().data() == shared_columns);
    REQUIRE(map->sparsity().row_offsets() == fixture.row_offsets);
    REQUIRE(map->sparsity().column_indices() == fixture.column_indices);
    TransportOperator transport(map->sparsity());
    map->assemble_values(fixture.atmosphere, 1, transport);
    REQUIRE(transport.values().isApprox(expected_values_at_wavelength_one(),
                                        2.0e-14));
}

TEST_CASE("Successive-orders ray transport applies spatial albedo at the "
          "ground intersection",
          "[successive_orders][ray_transport][ground][linearization]") {
    RayTransportFixture fixture;
    constexpr int wavelength = 1;
    Eigen::MatrixXd albedo(3, num_wavelengths);
    albedo << 0.1, 0.2, 0.3, 0.4, 0.7, 0.8;
    fixture.atmosphere.surface().set_spatial_lambertian_albedo(albedo);
    fixture.interpolation[0].ground_horizontal_weights = {{0, 0.25}, {2, 0.75}};
    const auto map = fixture.make_map();

    sasktran2::successive_orders::TransportOperator transport(map.sparsity());
    map.assemble_values(fixture.atmosphere, wavelength, transport);
    const double layer_zero_od = 0.7 * 0.1 + 0.2 * 0.4;
    const double layer_one_od = 0.4 * 0.4 + 1.1 * 0.7;
    const double transmission = std::exp(-(layer_zero_od + layer_one_od));
    const double intersection_albedo = 0.25 * 0.2 + 0.75 * 0.8;
    REQUIRE(transport.values()(3) ==
            Catch::Approx(0.8 * transmission * intersection_albedo)
                .epsilon(2.0e-14));

    Eigen::VectorXd tangent =
        Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
    tangent(fixture.atmosphere.surface_deriv_start_index()) = 0.4;
    tangent(fixture.atmosphere.surface_deriv_start_index() + 2) = -0.2;
    Eigen::VectorXd analytic(map.sparsity().nonzeros());
    map.assemble_jvp(fixture.atmosphere, wavelength, tangent, analytic);

    constexpr double step = 1.0e-6;
    sasktran2::successive_orders::TransportOperator above(map.sparsity());
    sasktran2::successive_orders::TransportOperator below(map.sparsity());
    Eigen::MatrixXd direction = Eigen::MatrixXd::Zero(3, num_wavelengths);
    direction(0, wavelength) = 0.4;
    direction(2, wavelength) = -0.2;
    fixture.atmosphere.surface().set_spatial_lambertian_albedo(
        albedo + step * direction);
    map.assemble_values(fixture.atmosphere, wavelength, above);
    fixture.atmosphere.surface().set_spatial_lambertian_albedo(
        albedo - step * direction);
    map.assemble_values(fixture.atmosphere, wavelength, below);
    fixture.atmosphere.surface().set_spatial_lambertian_albedo(albedo);
    const Eigen::VectorXd finite_difference =
        (above.values() - below.values()) / (2.0 * step);
    REQUIRE(analytic.isApprox(finite_difference, 2.0e-9));

    const Eigen::VectorXd value_gradient =
        (Eigen::VectorXd(6) << 0.2, -0.35, 0.5, 0.1, -0.4, 0.3).finished();
    Eigen::VectorXd native_gradient =
        Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
    sasktran2::successive_orders::RayTransportWorkspace workspace;
    map.accumulate_vjp(fixture.atmosphere, wavelength, value_gradient,
                       native_gradient, workspace);
    REQUIRE(analytic.dot(value_gradient) ==
            Catch::Approx(tangent.dot(native_gradient)).epsilon(2.0e-13));
}

TEST_CASE("Successive-orders packed ray transport JVP matches finite "
          "differences",
          "[successive_orders][ray_transport][linearization]") {
    RayTransportFixture fixture;
    const auto map = fixture.make_map();
    const int wavelength = 1;

    Eigen::VectorXd tangent =
        Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
    tangent.head(num_locations) << 0.08, -0.03, 0.06;
    tangent.segment(fixture.atmosphere.ssa_deriv_start_index(), num_locations)
        << -0.04,
        0.07, 0.02;

    Eigen::VectorXd analytic(map.sparsity().nonzeros());
    map.assemble_jvp(fixture.atmosphere, wavelength, tangent, analytic);

    constexpr double step = 1.0e-6;
    sasktran2::successive_orders::TransportOperator above(map.sparsity());
    sasktran2::successive_orders::TransportOperator below(map.sparsity());
    fixture.perturb(wavelength, tangent, step);
    map.assemble_values(fixture.atmosphere, wavelength, above);
    fixture.perturb(wavelength, tangent, -2.0 * step);
    map.assemble_values(fixture.atmosphere, wavelength, below);
    fixture.perturb(wavelength, tangent, step);

    const Eigen::VectorXd finite_difference =
        (above.values() - below.values()) / (2.0 * step);
    REQUIRE(analytic.isApprox(finite_difference, 2.0e-9));
}

TEST_CASE("Successive-orders packed ray transport VJP is adjoint and matches "
          "finite differences",
          "[successive_orders][ray_transport][linearization]") {
    RayTransportFixture fixture;
    const auto map = fixture.make_map();
    const int wavelength = 1;
    const Eigen::VectorXd value_gradient =
        (Eigen::VectorXd(6) << 0.2, -0.35, 0.5, 0.1, -0.4, 0.3).finished();
    Eigen::VectorXd tangent =
        Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
    tangent.head(num_locations) << 0.08, -0.03, 0.06;
    tangent.segment(fixture.atmosphere.ssa_deriv_start_index(), num_locations)
        << -0.04,
        0.07, 0.02;

    Eigen::VectorXd value_tangent(map.sparsity().nonzeros());
    map.assemble_jvp(fixture.atmosphere, wavelength, tangent, value_tangent);
    Eigen::VectorXd native_gradient =
        Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
    sasktran2::successive_orders::RayTransportWorkspace workspace;
    map.accumulate_vjp(fixture.atmosphere, wavelength, value_gradient,
                       native_gradient, workspace);

    REQUIRE(value_tangent.dot(value_gradient) ==
            Catch::Approx(tangent.dot(native_gradient)).epsilon(2.0e-13));
    REQUIRE(workspace.storage_bytes() ==
            static_cast<std::size_t>(5 * map.maximum_layers()) *
                sizeof(double));

    constexpr double step = 1.0e-6;
    for (int derivative = 0; derivative < 2 * num_locations; ++derivative) {
        Eigen::VectorXd coordinate =
            Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
        coordinate(derivative) = 1.0;
        sasktran2::successive_orders::TransportOperator above(map.sparsity());
        sasktran2::successive_orders::TransportOperator below(map.sparsity());
        fixture.perturb(wavelength, coordinate, step);
        map.assemble_values(fixture.atmosphere, wavelength, above);
        fixture.perturb(wavelength, coordinate, -2.0 * step);
        map.assemble_values(fixture.atmosphere, wavelength, below);
        fixture.perturb(wavelength, coordinate, step);
        const double finite_difference =
            value_gradient.dot(above.values() - below.values()) / (2.0 * step);
        REQUIRE(native_gradient(derivative) ==
                Catch::Approx(finite_difference).margin(2.0e-9));
    }

    REQUIRE(
        native_gradient.tail(fixture.atmosphere.num_deriv() - 2 * num_locations)
            .isZero());
}

TEST_CASE("Successive-orders packed ray transport preserves very thin layer "
          "sources and derivatives",
          "[successive_orders][ray_transport][linearization]") {
    RayTransportFixture fixture;
    const auto map = fixture.make_map();
    constexpr int wavelength = 0;
    constexpr double extinction = 1.0e-18;
    fixture.atmosphere.storage()
        .total_extinction.col(wavelength)
        .setConstant(extinction);

    sasktran2::successive_orders::TransportOperator transport(map.sparsity());
    map.assemble_values(fixture.atmosphere, wavelength, transport);

    const double optical_depth = 0.8 * extinction;
    const double source_fraction = -std::expm1(-optical_depth);
    const double albedo = 0.6;
    REQUIRE(transport.values()(4) > 0.0);
    REQUIRE(transport.values()(4) ==
            Catch::Approx(0.7 * albedo * source_fraction).margin(1.0e-30));
    REQUIRE(transport.values()(5) ==
            Catch::Approx(0.3 * albedo * source_fraction).margin(1.0e-30));

    Eigen::VectorXd tangent =
        Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
    tangent.segment(fixture.atmosphere.ssa_deriv_start_index(), num_locations)
        .setOnes();
    Eigen::VectorXd value_tangent(map.sparsity().nonzeros());
    map.assemble_jvp(fixture.atmosphere, wavelength, tangent, value_tangent);
    REQUIRE(value_tangent(4) > 0.0);
    REQUIRE(value_tangent(4) ==
            Catch::Approx(0.7 * source_fraction).margin(1.0e-30));
    REQUIRE(value_tangent(5) ==
            Catch::Approx(0.3 * source_fraction).margin(1.0e-30));

    Eigen::VectorXd value_gradient =
        Eigen::VectorXd::Zero(map.sparsity().nonzeros());
    value_gradient(4) = 1.0;
    Eigen::VectorXd native_gradient =
        Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
    sasktran2::successive_orders::RayTransportWorkspace workspace;
    map.accumulate_vjp(fixture.atmosphere, wavelength, value_gradient,
                       native_gradient, workspace);

    REQUIRE(native_gradient(fixture.atmosphere.ssa_deriv_start_index()) > 0.0);
    REQUIRE(value_tangent.dot(value_gradient) ==
            Catch::Approx(tangent.dot(native_gradient)).margin(1.0e-30));
}
