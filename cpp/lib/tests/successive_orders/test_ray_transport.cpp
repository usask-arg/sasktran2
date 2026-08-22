#include "../../successive_orders/ray_transport.h"

#include <sasktran2/test_helper.h>

#include <cmath>
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
                {0, 0.25, 0}, {2, 0.75, 2}, {1, 0.6, 1}, {2, 0.4, 2}};
            interpolation[0].optical_depth_indices = {0, 1, 1, 2};
            interpolation[0].optical_depth_weights = {0.7, 0.2, 0.4, 1.1};
            interpolation[0].ground_weights = {{0, 0.2, 0}, {3, 0.8, 3}};
            interpolation[0].ground_hit = true;
            interpolation[0].transport_value_offset = 0;
            interpolation[0].transport_row_nnz = 4;

            interpolation[1].layers = {{0, 2, 0, 2, 0, 2}};
            interpolation[1].atmosphere_weights = {{0, 0.5}, {2, 0.5}};
            interpolation[1].source_weights = {{1, 0.7, 0}, {3, 0.3, 1}};
            interpolation[1].optical_depth_indices = {0, 2};
            interpolation[1].optical_depth_weights = {0.3, 0.5};
            interpolation[1].transport_value_offset = 4;
            interpolation[1].transport_row_nnz = 2;
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
