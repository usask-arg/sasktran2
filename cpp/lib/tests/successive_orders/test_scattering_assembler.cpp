#include "../../successive_orders/scattering_assembler.h"

#include <sasktran2/hr/diffuse_point.h>
#include <sasktran2/test_helper.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
    using namespace sasktran2::successive_orders;

    Eigen::MatrixXd
    materialize_atmospheric_block(const ScatteringOperator<3>& scattering,
                                  int block) {
        const int rows = scattering.layout().output_block_size(block);
        const int columns = scattering.layout().input_block_size(block);
        Eigen::MatrixXd result(rows, columns);
        Eigen::VectorXd incoming =
            Eigen::VectorXd::Zero(scattering.input_size());
        Eigen::VectorXd outgoing(scattering.output_size());
        auto workspace = scattering.make_workspace();
        for (int column = 0; column < columns; ++column) {
            incoming(scattering.input_offsets()[block] + column) = 1.0;
            scattering.apply(incoming, outgoing, workspace);
            result.col(column) =
                outgoing.segment(scattering.output_offsets()[block], rows);
            incoming(scattering.input_offsets()[block] + column) = 0.0;
        }
        return result;
    }

    template <int NSTOKES>
    class FiniteAzimuthBRDF final
        : public sasktran2::atmosphere::brdf::BRDF<NSTOKES> {
      public:
        Eigen::Matrix<double, NSTOKES, NSTOKES>
        brdf(double, double, double phi_difference,
             Eigen::Ref<const Eigen::VectorXd> args) const override {
            if (!std::isfinite(phi_difference)) {
                throw std::runtime_error("BRDF received a non-finite azimuth");
            }
            Eigen::Matrix<double, NSTOKES, NSTOKES> result =
                Eigen::Matrix<double, NSTOKES, NSTOKES>::Zero();
            result(0, 0) =
                args(0) * (1.0 + 0.25 * std::cos(phi_difference)) / EIGEN_PI;
            return result;
        }

        int num_deriv() const override { return 0; }
        int num_args() const override { return 1; }
    };

    sasktran2::Geometry1D make_geometry() {
        sasktran2::Coordinates coordinates(
            0.6, 0.2, 6372000.0, sasktran2::geometrytype::planeparallel);
        Eigen::VectorXd altitudes(3);
        altitudes << 0.0, 1000.0, 2000.0;
        sasktran2::grids::AltitudeGrid altitude_grid(
            std::move(altitudes), sasktran2::grids::gridspacing::constant,
            sasktran2::grids::outofbounds::extend,
            sasktran2::grids::interpolation::linear);
        return {std::move(coordinates), std::move(altitude_grid)};
    }

    struct AssemblyFixture {
        AssemblyFixture()
            : geometry(make_geometry()), raytracer(geometry),
              source_geometry(raytracer, geometry),
              atmosphere(sasktran2::atmosphere::AtmosphereGridStorageFull<1>(
                             num_wavelengths, geometry.size(), num_legendre),
                         sasktran2::atmosphere::Surface<1>(num_wavelengths),
                         true) {
            SourceGeometrySettings settings;
            settings.num_incoming = 6;
            settings.num_outgoing = 6;
            settings.num_sza = 1;
            settings.num_threads = 1;
            settings.altitude_grid_m = {500.0, 1500.0};
            sasktran2::viewinggeometry::InternalViewingGeometry viewing;
            source_geometry.initialize(viewing, settings);

            atmosphere.storage().resize_derivatives(num_scattering_groups);
            for (int wavelength = 0; wavelength < num_wavelengths;
                 ++wavelength) {
                for (int location = 0; location < geometry.size(); ++location) {
                    for (int degree = 0; degree < num_legendre; ++degree) {
                        atmosphere.storage().leg_coeff(degree, location,
                                                       wavelength) =
                            0.12 + 0.025 * degree + 0.04 * location +
                            0.015 * wavelength;
                        for (int group = 0; group < num_scattering_groups;
                             ++group) {
                            atmosphere.storage().d_leg_coeff(
                                degree, location, wavelength, group) =
                                0.01 * (1 + degree + 2 * location + 3 * group +
                                        wavelength);
                        }
                    }
                }
            }
            atmosphere.surface().brdf_args().row(0) << 0.24, 0.37;
        }

        Eigen::VectorXd native_tangent() const {
            Eigen::VectorXd tangent =
                Eigen::VectorXd::Zero(atmosphere.num_deriv());
            tangent.head(geometry.size()) << 0.2, -0.1, 0.05;
            tangent.segment(atmosphere.ssa_deriv_start_index(), geometry.size())
                << -0.03,
                0.04, 0.07;
            for (int group = 0; group < num_scattering_groups; ++group) {
                tangent.segment(atmosphere.scat_deriv_start_index() +
                                    group * geometry.size(),
                                geometry.size())
                    << 0.13 - 0.02 * group,
                    -0.08 + 0.03 * group, 0.05 + 0.01 * group;
            }
            tangent(atmosphere.surface_deriv_start_index()) = -0.17;
            return tangent;
        }

        void perturb(int wavelength, const Eigen::VectorXd& tangent,
                     double scale) {
            for (int location = 0; location < geometry.size(); ++location) {
                for (int degree = 0; degree < num_legendre; ++degree) {
                    double direction = 0.0;
                    for (int group = 0; group < num_scattering_groups;
                         ++group) {
                        direction +=
                            tangent(atmosphere.scat_deriv_start_index() +
                                    group * geometry.size() + location) *
                            atmosphere.storage().d_leg_coeff(degree, location,
                                                             wavelength, group);
                    }
                    atmosphere.storage().leg_coeff(
                        degree, location, wavelength) += scale * direction;
                }
            }
            atmosphere.surface().brdf_args()(0, wavelength) +=
                scale * tangent(atmosphere.surface_deriv_start_index());
        }

        void enable_delta_m_scaling(int order) {
            atmosphere.storage().total_extinction.setConstant(1.0e-4);
            atmosphere.storage().ssa.setConstant(0.9);
            for (int wavelength = 0; wavelength < num_wavelengths;
                 ++wavelength) {
                for (int location = 0; location < geometry.size(); ++location) {
                    atmosphere.storage().leg_coeff(0, location, wavelength) =
                        1.0;
                    for (int group = 0; group < num_scattering_groups;
                         ++group) {
                        atmosphere.storage().d_leg_coeff(
                            0, location, wavelength, group) = 0.0;
                    }
                }
            }
            atmosphere.apply_delta_m_scaling(order);
        }

        void perturb_delta_scaled(int wavelength,
                                  const Eigen::VectorXd& tangent,
                                  double scale) {
            for (int location = 0; location < geometry.size(); ++location) {
                double f_direction = 0.0;
                for (int group = 0; group < num_scattering_groups; ++group) {
                    const double direction =
                        tangent(atmosphere.scat_deriv_start_index() +
                                group * geometry.size() + location);
                    f_direction += direction * atmosphere.storage().d_f(
                                                   location, wavelength, group);
                    for (int degree = 0; degree < num_legendre; ++degree) {
                        atmosphere.storage().leg_coeff(degree, location,
                                                       wavelength) +=
                            scale * direction *
                            atmosphere.storage().d_leg_coeff(degree, location,
                                                             wavelength, group);
                    }
                }
                atmosphere.storage().f(location, wavelength) +=
                    scale * f_direction;
            }
            atmosphere.surface().brdf_args()(0, wavelength) +=
                scale * tangent(atmosphere.surface_deriv_start_index());
        }

        static constexpr int num_wavelengths = 2;
        static constexpr int num_legendre = 4;
        static constexpr int num_scattering_groups = 2;

        sasktran2::Geometry1D geometry;
        sasktran2::raytracing::PlaneParallelRayTracer raytracer;
        SourceGeometry1D source_geometry;
        sasktran2::atmosphere::Atmosphere<1> atmosphere;
    };

    class CopiedUnitSphere final : public sasktran2::math::UnitSphere {
      public:
        explicit CopiedUnitSphere(const sasktran2::math::UnitSphere& source) {
            m_directions.reserve(source.num_points());
            m_weights.reserve(source.num_points());
            for (int index = 0; index < source.num_points(); ++index) {
                m_directions.push_back(source.get_quad_position(index));
                m_weights.push_back(source.quadrature_weight(index));
            }
        }

        int num_points() const override {
            return static_cast<int>(m_weights.size());
        }
        Eigen::Vector3d get_quad_position(int index) const override {
            return m_directions.at(index);
        }
        double quadrature_weight(int index) const override {
            return m_weights.at(index);
        }
        void interpolate(const Eigen::Vector3d&,
                         std::vector<std::pair<int, double>>&,
                         int&) const override {
            throw std::logic_error(
                "CopiedUnitSphere interpolation is not used by this test");
        }

      private:
        std::vector<Eigen::Vector3d> m_directions;
        std::vector<double> m_weights;
    };

    struct VectorAssemblyFixture {
        VectorAssemblyFixture()
            : geometry(make_geometry()), raytracer(geometry),
              source_geometry(raytracer, geometry),
              atmosphere(sasktran2::atmosphere::AtmosphereGridStorageFull<3>(
                             num_wavelengths, geometry.size(), num_legendre),
                         sasktran2::atmosphere::Surface<3>(num_wavelengths),
                         true) {
            SourceGeometrySettings settings;
            settings.num_incoming = 6;
            settings.num_outgoing = 6;
            settings.num_sza = 1;
            settings.num_threads = 1;
            settings.altitude_grid_m = {500.0, 1500.0};
            sasktran2::viewinggeometry::InternalViewingGeometry viewing;
            source_geometry.initialize(viewing, settings);

            atmosphere.storage().resize_derivatives(num_scattering_groups);
            for (int wavelength = 0; wavelength < num_wavelengths;
                 ++wavelength) {
                for (int location = 0; location < geometry.size(); ++location) {
                    for (int degree = 0; degree < num_legendre; ++degree) {
                        for (int greek = 0; greek < 4; ++greek) {
                            const int coefficient = 4 * degree + greek;
                            atmosphere.storage().leg_coeff(
                                coefficient, location, wavelength) =
                                0.09 + 0.017 * degree + 0.011 * greek +
                                0.023 * location + 0.007 * wavelength;
                            for (int group = 0; group < num_scattering_groups;
                                 ++group) {
                                atmosphere.storage().d_leg_coeff(
                                    coefficient, location, wavelength, group) =
                                    0.004 * (1 + coefficient + 2 * location +
                                             3 * group + wavelength);
                            }
                        }
                    }
                }
            }
            atmosphere.surface().brdf_args().row(0) << 0.19, 0.31;
        }

        Eigen::VectorXd native_tangent() const {
            Eigen::VectorXd tangent =
                Eigen::VectorXd::Zero(atmosphere.num_deriv());
            tangent.head(geometry.size()) << 0.1, -0.03, 0.02;
            tangent.segment(atmosphere.ssa_deriv_start_index(), geometry.size())
                << -0.04,
                0.06, 0.01;
            for (int group = 0; group < num_scattering_groups; ++group) {
                tangent.segment(atmosphere.scat_deriv_start_index() +
                                    group * geometry.size(),
                                geometry.size())
                    << 0.12 - 0.03 * group,
                    -0.07 + 0.02 * group, 0.045 + 0.01 * group;
            }
            tangent(atmosphere.surface_deriv_start_index()) = -0.14;
            return tangent;
        }

        void perturb(int wavelength, const Eigen::VectorXd& tangent,
                     double scale) {
            for (int location = 0; location < geometry.size(); ++location) {
                for (int coefficient = 0; coefficient < 4 * num_legendre;
                     ++coefficient) {
                    double direction = 0.0;
                    for (int group = 0; group < num_scattering_groups;
                         ++group) {
                        direction +=
                            tangent(atmosphere.scat_deriv_start_index() +
                                    group * geometry.size() + location) *
                            atmosphere.storage().d_leg_coeff(
                                coefficient, location, wavelength, group);
                    }
                    atmosphere.storage().leg_coeff(
                        coefficient, location, wavelength) += scale * direction;
                }
            }
            atmosphere.surface().brdf_args()(0, wavelength) +=
                scale * tangent(atmosphere.surface_deriv_start_index());
        }

        void enable_delta_m_scaling(int order) {
            atmosphere.storage().total_extinction.setConstant(1.0e-4);
            atmosphere.storage().ssa.setConstant(0.9);
            for (int wavelength = 0; wavelength < num_wavelengths;
                 ++wavelength) {
                for (int location = 0; location < geometry.size(); ++location) {
                    for (int greek = 0; greek < 3; ++greek) {
                        atmosphere.storage().leg_coeff(greek, location,
                                                       wavelength) = 1.0;
                        for (int group = 0; group < num_scattering_groups;
                             ++group) {
                            atmosphere.storage().d_leg_coeff(
                                greek, location, wavelength, group) = 0.0;
                        }
                    }
                }
            }
            atmosphere.apply_delta_m_scaling(order);
        }

        void perturb_delta_scaled(int wavelength,
                                  const Eigen::VectorXd& tangent,
                                  double scale) {
            for (int location = 0; location < geometry.size(); ++location) {
                double f_direction = 0.0;
                for (int group = 0; group < num_scattering_groups; ++group) {
                    const double direction =
                        tangent(atmosphere.scat_deriv_start_index() +
                                group * geometry.size() + location);
                    f_direction += direction * atmosphere.storage().d_f(
                                                   location, wavelength, group);
                    for (int coefficient = 0; coefficient < 4 * num_legendre;
                         ++coefficient) {
                        atmosphere.storage().leg_coeff(coefficient, location,
                                                       wavelength) +=
                            scale * direction *
                            atmosphere.storage().d_leg_coeff(
                                coefficient, location, wavelength, group);
                    }
                }
                atmosphere.storage().f(location, wavelength) +=
                    scale * f_direction;
            }
            atmosphere.surface().brdf_args()(0, wavelength) +=
                scale * tangent(atmosphere.surface_deriv_start_index());
        }

        Eigen::MatrixXd legacy_block(int point_index, int wavelength,
                                     int num_coefficients) const {
            const auto& point = source_geometry.source_point(point_index);
            sasktran2::hr::IncomingOutgoingSpherePair<3> legacy(
                num_coefficients,
                std::make_unique<CopiedUnitSphere>(point.incoming_sphere()),
                std::make_unique<CopiedUnitSphere>(point.outgoing_sphere()));
            std::vector<std::pair<int, double>> weights;
            weights.reserve(point.atmosphere_weights().size());
            for (const auto& weight : point.atmosphere_weights()) {
                weights.emplace_back(weight.index, weight.weight);
            }
            Eigen::MatrixXd result(3 * point.num_outgoing(),
                                   3 * point.num_incoming());
            if (point.is_ground()) {
                legacy.calculate_ground_scattering_matrix(
                    atmosphere.surface(), weights, point.location(), wavelength,
                    result.data());
            } else {
                legacy.calculate_scattering_matrix(
                    atmosphere.storage(), wavelength, weights, result.data());
            }
            return result;
        }

        static constexpr int num_wavelengths = 2;
        static constexpr int num_legendre = 4;
        static constexpr int num_scattering_groups = 2;

        sasktran2::Geometry1D geometry;
        sasktran2::raytracing::PlaneParallelRayTracer raytracer;
        SourceGeometry1D source_geometry;
        sasktran2::atmosphere::Atmosphere<3> atmosphere;
    };
} // namespace

TEST_CASE("Scalar scattering assembler maps atmosphere points and legacy "
          "ground normalization",
          "[successive_orders][scattering_assembler]") {
    AssemblyFixture fixture;
    constexpr int wavelength = 1;
    ScalarScatteringAssembler assembler(fixture.source_geometry, 3);
    auto scattering = assembler.create_operator();
    assembler.assemble_values(fixture.atmosphere, wavelength, scattering);

    REQUIRE(scattering.input_offsets() ==
            fixture.source_geometry.incoming_point_offsets());
    REQUIRE(scattering.output_offsets() ==
            fixture.source_geometry.outgoing_point_offsets());
    for (int point = 0; point < fixture.source_geometry.num_interior_points();
         ++point) {
        for (int degree = 0; degree < assembler.num_coefficients(); ++degree) {
            double expected = 0.0;
            for (const auto& weight :
                 fixture.source_geometry.source_point(point)
                     .atmosphere_weights()) {
                expected +=
                    weight.weight * fixture.atmosphere.storage().leg_coeff(
                                        degree, weight.index, wavelength);
            }
            REQUIRE(scattering.atmospheric_coefficients()(point, degree) ==
                    Catch::Approx(expected).epsilon(2.0e-14));
        }
    }

    REQUIRE(fixture.source_geometry.num_ground_points() == 1);
    const int ground_point = fixture.source_geometry.num_interior_points();
    const auto& point = fixture.source_geometry.source_point(ground_point);
    const auto ground = scattering.ground_block(0);
    const double albedo =
        fixture.atmosphere.surface().brdf_args()(0, wavelength);
    for (int input = 0; input < point.num_incoming(); ++input) {
        const double mu_in = point.location().cos_zenith_angle(
            point.incoming_sphere().get_quad_position(input));
        const double expected =
            4.0 * albedo * mu_in *
            point.incoming_sphere().quadrature_weight(input);
        for (int output = 0; output < point.num_outgoing(); ++output) {
            REQUIRE(ground(output, input) ==
                    Catch::Approx(expected).epsilon(2.0e-14));
        }
    }
    REQUIRE(assembler.storage_bytes() >
            scattering.memory_usage().angular_basis_bytes);
}

TEST_CASE("Spatial Lambertian ground source points use unit albedo",
          "[successive_orders][scattering_assembler][ground][linearization]") {
    SECTION("scalar") {
        AssemblyFixture fixture;
        Eigen::MatrixXd albedo(2, AssemblyFixture::num_wavelengths);
        albedo << 0.15, 0.25, 0.65, 0.75;
        fixture.atmosphere.surface().set_spatial_lambertian_albedo(albedo);
        ScalarScatteringAssembler assembler(fixture.source_geometry, 3);
        auto first = assembler.create_operator();
        assembler.assemble_values(fixture.atmosphere, 1, first);

        fixture.atmosphere.surface().set_spatial_lambertian_albedo(
            Eigen::MatrixXd::Constant(2, AssemblyFixture::num_wavelengths,
                                      0.95));
        auto second = assembler.create_operator();
        assembler.assemble_values(fixture.atmosphere, 1, second);
        REQUIRE(first.ground_values().isApprox(second.ground_values(), 0.0));

        Eigen::VectorXd tangent =
            Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
        tangent.segment(fixture.atmosphere.surface_deriv_start_index(), 2)
            << 0.3,
            -0.2;
        Eigen::MatrixXd coefficient_tangent;
        Eigen::VectorXd ground_tangent;
        assembler.assemble_jvp(fixture.atmosphere, 1, tangent,
                               coefficient_tangent, ground_tangent);
        REQUIRE(ground_tangent.isZero());

        Eigen::VectorXd native_gradient =
            Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
        assembler.accumulate_vjp(
            fixture.atmosphere, 1,
            Eigen::MatrixXd::Zero(fixture.source_geometry.num_interior_points(),
                                  assembler.num_coefficients()),
            Eigen::VectorXd::Ones(assembler.ground_value_size()),
            native_gradient);
        REQUIRE(native_gradient
                    .segment(fixture.atmosphere.surface_deriv_start_index(), 2)
                    .isZero());
    }

    SECTION("vector") {
        VectorAssemblyFixture fixture;
        Eigen::MatrixXd albedo(2, VectorAssemblyFixture::num_wavelengths);
        albedo << 0.15, 0.25, 0.65, 0.75;
        fixture.atmosphere.surface().set_spatial_lambertian_albedo(albedo);
        VectorScatteringAssembler assembler(fixture.source_geometry, 3);
        auto first = assembler.create_operator();
        assembler.assemble_values(fixture.atmosphere, 1, first);

        fixture.atmosphere.surface().set_spatial_lambertian_albedo(
            Eigen::MatrixXd::Constant(2, VectorAssemblyFixture::num_wavelengths,
                                      0.95));
        auto second = assembler.create_operator();
        assembler.assemble_values(fixture.atmosphere, 1, second);
        REQUIRE(first.ground_values().isApprox(second.ground_values(), 0.0));
    }
}

TEST_CASE("Ground scattering supplies finite azimuths for vertical directions",
          "[successive_orders][scattering_assembler][ground]") {
    SECTION("scalar") {
        AssemblyFixture fixture;
        const int ground_point = fixture.source_geometry.num_interior_points();
        const auto& point = fixture.source_geometry.source_point(ground_point);
        const Eigen::Vector3d normal = point.location().position.normalized();
        bool has_vertical_direction = false;
        for (int input = 0; input < point.num_incoming(); ++input) {
            const Eigen::Vector3d direction =
                point.incoming_sphere().get_quad_position(input);
            has_vertical_direction =
                has_vertical_direction ||
                (direction - direction.dot(normal) * normal).norm() < 1.0e-14;
        }
        REQUIRE(has_vertical_direction);

        fixture.atmosphere.surface().set_brdf_object(
            std::make_shared<FiniteAzimuthBRDF<1>>());
        fixture.atmosphere.surface().brdf_args().setConstant(0.2);
        ScalarScatteringAssembler assembler(fixture.source_geometry, 3);
        auto scattering = assembler.create_operator();
        REQUIRE_NOTHROW(
            assembler.assemble_values(fixture.atmosphere, 0, scattering));
        REQUIRE(scattering.ground_values().allFinite());
    }

    SECTION("vector") {
        VectorAssemblyFixture fixture;
        fixture.atmosphere.surface().set_brdf_object(
            std::make_shared<FiniteAzimuthBRDF<3>>());
        fixture.atmosphere.surface().brdf_args().setConstant(0.2);
        VectorScatteringAssembler assembler(fixture.source_geometry, 3);
        auto scattering = assembler.create_operator();
        REQUIRE_NOTHROW(
            assembler.assemble_values(fixture.atmosphere, 0, scattering));
        REQUIRE(scattering.ground_values().allFinite());
    }
}

TEST_CASE("Scalar scattering assembler native JVP matches finite differences",
          "[successive_orders][scattering_assembler][jvp]") {
    AssemblyFixture fixture;
    constexpr int wavelength = 1;
    ScalarScatteringAssembler assembler(fixture.source_geometry, 3);
    const Eigen::VectorXd tangent = fixture.native_tangent();
    Eigen::MatrixXd coefficient_tangent;
    Eigen::VectorXd ground_tangent;
    assembler.assemble_jvp(fixture.atmosphere, wavelength, tangent,
                           coefficient_tangent, ground_tangent);

    constexpr double step = 1.0e-5;
    fixture.perturb(wavelength, tangent, step);
    auto above = assembler.create_operator();
    assembler.assemble_values(fixture.atmosphere, wavelength, above);
    fixture.perturb(wavelength, tangent, -2.0 * step);
    auto below = assembler.create_operator();
    assembler.assemble_values(fixture.atmosphere, wavelength, below);
    fixture.perturb(wavelength, tangent, step);

    const Eigen::MatrixXd coefficient_difference =
        (above.atmospheric_coefficients() - below.atmospheric_coefficients()) /
        (2.0 * step);
    const Eigen::VectorXd ground_difference =
        (above.ground_values() - below.ground_values()) / (2.0 * step);
    REQUIRE(
        (coefficient_tangent - coefficient_difference).cwiseAbs().maxCoeff() <
        2.0e-10);
    REQUIRE(ground_tangent.isApprox(ground_difference, 2.0e-10));
}

TEST_CASE("Scalar scattering assembler native VJP is adjoint",
          "[successive_orders][scattering_assembler][jvp][vjp]") {
    AssemblyFixture fixture;
    constexpr int wavelength = 0;
    ScalarScatteringAssembler assembler(fixture.source_geometry, 3);
    const Eigen::VectorXd tangent = fixture.native_tangent();
    Eigen::MatrixXd coefficient_tangent;
    Eigen::VectorXd ground_tangent;
    assembler.assemble_jvp(fixture.atmosphere, wavelength, tangent,
                           coefficient_tangent, ground_tangent);

    const Eigen::MatrixXd coefficient_gradient =
        Eigen::MatrixXd::Random(fixture.source_geometry.num_interior_points(),
                                assembler.num_coefficients());
    const Eigen::VectorXd ground_gradient =
        Eigen::VectorXd::Random(assembler.ground_value_size());
    Eigen::VectorXd native_gradient =
        Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
    assembler.accumulate_vjp(fixture.atmosphere, wavelength,
                             coefficient_gradient, ground_gradient,
                             native_gradient);

    const double parameter_forward =
        (coefficient_tangent.array() * coefficient_gradient.array()).sum() +
        ground_tangent.dot(ground_gradient);
    REQUIRE(parameter_forward ==
            Catch::Approx(tangent.dot(native_gradient)).epsilon(3.0e-13));
    REQUIRE(native_gradient.head(2 * fixture.geometry.size()).isZero());
}

TEST_CASE("Scalar scattering assembler completes delta-M scaling",
          "[successive_orders][scattering_assembler][delta_m][jvp][vjp]") {
    AssemblyFixture fixture;
    constexpr int wavelength = 1;
    constexpr int delta_m_order = 3;
    fixture.enable_delta_m_scaling(delta_m_order);

    ScalarScatteringAssembler assembler(fixture.source_geometry, delta_m_order);
    auto scattering = assembler.create_operator();
    assembler.assemble_values(fixture.atmosphere, wavelength, scattering);

    const auto& storage = fixture.atmosphere.storage();
    for (int point = 0; point < fixture.source_geometry.num_interior_points();
         ++point) {
        for (int degree = 0; degree < assembler.num_coefficients(); ++degree) {
            double expected = 0.0;
            for (const auto& weight :
                 fixture.source_geometry.source_point(point)
                     .atmosphere_weights()) {
                const double f = storage.f(weight.index, wavelength);
                expected +=
                    weight.weight *
                    (storage.leg_coeff(degree, weight.index, wavelength) -
                     (2.0 * degree + 1.0) * f / (1.0 - f));
            }
            REQUIRE(scattering.atmospheric_coefficients()(point, degree) ==
                    Catch::Approx(expected).epsilon(2.0e-14));
        }
        REQUIRE(scattering.atmospheric_coefficients()(point, 0) ==
                Catch::Approx(1.0).epsilon(2.0e-14));
    }

    const Eigen::VectorXd tangent = fixture.native_tangent();
    Eigen::MatrixXd coefficient_tangent;
    Eigen::VectorXd ground_tangent;
    assembler.assemble_jvp(fixture.atmosphere, wavelength, tangent,
                           coefficient_tangent, ground_tangent);
    REQUIRE(coefficient_tangent.col(0).isZero(2.0e-14));

    constexpr double step = 1.0e-5;
    fixture.perturb_delta_scaled(wavelength, tangent, step);
    auto above = assembler.create_operator();
    assembler.assemble_values(fixture.atmosphere, wavelength, above);
    fixture.perturb_delta_scaled(wavelength, tangent, -2.0 * step);
    auto below = assembler.create_operator();
    assembler.assemble_values(fixture.atmosphere, wavelength, below);
    fixture.perturb_delta_scaled(wavelength, tangent, step);

    const Eigen::MatrixXd coefficient_difference =
        (above.atmospheric_coefficients() - below.atmospheric_coefficients()) /
        (2.0 * step);
    const Eigen::VectorXd ground_difference =
        (above.ground_values() - below.ground_values()) / (2.0 * step);
    REQUIRE(
        (coefficient_tangent - coefficient_difference).cwiseAbs().maxCoeff() <
        2.0e-10);
    REQUIRE(ground_tangent.isApprox(ground_difference, 2.0e-10));

    const Eigen::MatrixXd coefficient_gradient =
        Eigen::MatrixXd::Random(fixture.source_geometry.num_interior_points(),
                                assembler.num_coefficients());
    const Eigen::VectorXd ground_gradient =
        Eigen::VectorXd::Random(assembler.ground_value_size());
    Eigen::VectorXd native_gradient =
        Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
    assembler.accumulate_vjp(fixture.atmosphere, wavelength,
                             coefficient_gradient, ground_gradient,
                             native_gradient);
    const double parameter_forward =
        (coefficient_tangent.array() * coefficient_gradient.array()).sum() +
        ground_tangent.dot(ground_gradient);
    REQUIRE(parameter_forward ==
            Catch::Approx(tangent.dot(native_gradient)).epsilon(3.0e-13));
}

TEST_CASE("Vector scattering assembler matches legacy atmospheric and ground "
          "blocks",
          "[successive_orders][scattering_assembler][vector]") {
    VectorAssemblyFixture fixture;
    constexpr int wavelength = 1;
    constexpr int num_coefficients = 3;
    VectorScatteringAssembler assembler(fixture.source_geometry,
                                        num_coefficients);
    auto scattering = assembler.create_operator();
    assembler.assemble_values(fixture.atmosphere, wavelength, scattering);

    std::vector<int> expected_input_offsets =
        fixture.source_geometry.incoming_point_offsets();
    std::vector<int> expected_output_offsets =
        fixture.source_geometry.outgoing_point_offsets();
    for (auto& offset : expected_input_offsets) {
        offset *= 3;
    }
    for (auto& offset : expected_output_offsets) {
        offset *= 3;
    }
    REQUIRE(scattering.input_offsets() == expected_input_offsets);
    REQUIRE(scattering.output_offsets() == expected_output_offsets);

    for (int point = 0; point < fixture.source_geometry.num_interior_points();
         ++point) {
        const Eigen::MatrixXd expected =
            fixture.legacy_block(point, wavelength, num_coefficients);
        const Eigen::MatrixXd actual =
            materialize_atmospheric_block(scattering, point);
        INFO("point=" << point << " max error="
                      << (actual - expected).cwiseAbs().maxCoeff());
        REQUIRE(actual.isApprox(expected, 5.0e-13));
    }
    for (int ground = 0; ground < fixture.source_geometry.num_ground_points();
         ++ground) {
        const int point =
            fixture.source_geometry.num_interior_points() + ground;
        const Eigen::MatrixXd expected =
            fixture.legacy_block(point, wavelength, num_coefficients);
        REQUIRE(scattering.ground_block(ground).isApprox(expected, 5.0e-13));
    }

    const auto& atmospheric_point = fixture.source_geometry.source_point(0);
    const std::size_t legacy_dense_basis_bytes =
        static_cast<std::size_t>(4 * num_coefficients) * 3 *
        atmospheric_point.num_outgoing() * 3 *
        atmospheric_point.num_incoming() * sizeof(double);
    REQUIRE(assembler.storage_bytes() < legacy_dense_basis_bytes);
    REQUIRE(scattering.memory_usage().atmospheric_value_bytes ==
            static_cast<std::size_t>(assembler.coefficient_value_size()) *
                sizeof(double));
}

TEST_CASE("Vector scattering assembler native JVP matches finite differences",
          "[successive_orders][scattering_assembler][vector][jvp]") {
    VectorAssemblyFixture fixture;
    constexpr int wavelength = 1;
    VectorScatteringAssembler assembler(fixture.source_geometry, 3);
    const Eigen::VectorXd tangent = fixture.native_tangent();
    Eigen::MatrixXd atmospheric_tangent;
    Eigen::VectorXd ground_tangent;
    assembler.assemble_jvp(fixture.atmosphere, wavelength, tangent,
                           atmospheric_tangent, ground_tangent);

    constexpr double step = 1.0e-5;
    fixture.perturb(wavelength, tangent, step);
    auto above = assembler.create_operator();
    assembler.assemble_values(fixture.atmosphere, wavelength, above);
    fixture.perturb(wavelength, tangent, -2.0 * step);
    auto below = assembler.create_operator();
    assembler.assemble_values(fixture.atmosphere, wavelength, below);
    fixture.perturb(wavelength, tangent, step);

    const Eigen::MatrixXd atmospheric_difference =
        (above.atmospheric_coefficients() - below.atmospheric_coefficients()) /
        (2.0 * step);
    const Eigen::VectorXd ground_difference =
        (above.ground_values() - below.ground_values()) / (2.0 * step);
    REQUIRE(
        (atmospheric_tangent - atmospheric_difference).cwiseAbs().maxCoeff() <
        3.0e-10);
    REQUIRE((ground_tangent - ground_difference).cwiseAbs().maxCoeff() <
            3.0e-10);
}

TEST_CASE("Vector scattering assembler native VJP is adjoint",
          "[successive_orders][scattering_assembler][vector][jvp][vjp]") {
    VectorAssemblyFixture fixture;
    constexpr int wavelength = 0;
    VectorScatteringAssembler assembler(fixture.source_geometry, 3);
    const Eigen::VectorXd tangent = fixture.native_tangent();
    Eigen::MatrixXd atmospheric_tangent;
    Eigen::VectorXd ground_tangent;
    assembler.assemble_jvp(fixture.atmosphere, wavelength, tangent,
                           atmospheric_tangent, ground_tangent);

    const Eigen::MatrixXd atmospheric_gradient =
        Eigen::MatrixXd::Random(fixture.source_geometry.num_interior_points(),
                                4 * assembler.num_coefficients());
    const Eigen::VectorXd ground_gradient =
        Eigen::VectorXd::Random(assembler.ground_value_size());
    Eigen::VectorXd native_gradient =
        Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
    assembler.accumulate_vjp(fixture.atmosphere, wavelength,
                             atmospheric_gradient, ground_gradient,
                             native_gradient);

    const double parameter_forward =
        (atmospheric_tangent.array() * atmospheric_gradient.array()).sum() +
        ground_tangent.dot(ground_gradient);
    REQUIRE(parameter_forward ==
            Catch::Approx(tangent.dot(native_gradient)).epsilon(5.0e-12));
    REQUIRE(native_gradient.head(2 * fixture.geometry.size()).isZero());
}

TEST_CASE("Vector scattering assembler completes delta-M scaling",
          "[successive_orders][scattering_assembler][vector][delta_m][jvp]"
          "[vjp]") {
    VectorAssemblyFixture fixture;
    constexpr int wavelength = 1;
    constexpr int delta_m_order = 3;
    fixture.enable_delta_m_scaling(delta_m_order);

    VectorScatteringAssembler assembler(fixture.source_geometry, delta_m_order);
    auto scattering = assembler.create_operator();
    assembler.assemble_values(fixture.atmosphere, wavelength, scattering);

    const auto& storage = fixture.atmosphere.storage();
    for (int point = 0; point < fixture.source_geometry.num_interior_points();
         ++point) {
        for (int coefficient = 0;
             coefficient < scattering.atmospheric_coefficients().cols();
             ++coefficient) {
            const int degree = coefficient / 4;
            const bool subtract_delta_peak = coefficient % 4 != 3;
            double expected = 0.0;
            for (const auto& weight :
                 fixture.source_geometry.source_point(point)
                     .atmosphere_weights()) {
                const double f = storage.f(weight.index, wavelength);
                const double correction =
                    subtract_delta_peak ? (2.0 * degree + 1.0) * f / (1.0 - f)
                                        : 0.0;
                expected +=
                    weight.weight *
                    (storage.leg_coeff(coefficient, weight.index, wavelength) -
                     correction);
            }
            REQUIRE(scattering.atmospheric_coefficients()(point, coefficient) ==
                    Catch::Approx(expected).epsilon(2.0e-14));
        }
        for (int greek = 0; greek < 3; ++greek) {
            REQUIRE(scattering.atmospheric_coefficients()(point, greek) ==
                    Catch::Approx(1.0).epsilon(2.0e-14));
        }
    }

    const Eigen::VectorXd tangent = fixture.native_tangent();
    Eigen::MatrixXd coefficient_tangent;
    Eigen::VectorXd ground_tangent;
    assembler.assemble_jvp(fixture.atmosphere, wavelength, tangent,
                           coefficient_tangent, ground_tangent);
    for (int greek = 0; greek < 3; ++greek) {
        REQUIRE(coefficient_tangent.col(greek).isZero(2.0e-14));
    }

    constexpr double step = 1.0e-5;
    fixture.perturb_delta_scaled(wavelength, tangent, step);
    auto above = assembler.create_operator();
    assembler.assemble_values(fixture.atmosphere, wavelength, above);
    fixture.perturb_delta_scaled(wavelength, tangent, -2.0 * step);
    auto below = assembler.create_operator();
    assembler.assemble_values(fixture.atmosphere, wavelength, below);
    fixture.perturb_delta_scaled(wavelength, tangent, step);

    const Eigen::MatrixXd coefficient_difference =
        (above.atmospheric_coefficients() - below.atmospheric_coefficients()) /
        (2.0 * step);
    const Eigen::VectorXd ground_difference =
        (above.ground_values() - below.ground_values()) / (2.0 * step);
    REQUIRE(
        (coefficient_tangent - coefficient_difference).cwiseAbs().maxCoeff() <
        3.0e-10);
    REQUIRE((ground_tangent - ground_difference).cwiseAbs().maxCoeff() <
            3.0e-10);

    const Eigen::MatrixXd coefficient_gradient =
        Eigen::MatrixXd::Random(fixture.source_geometry.num_interior_points(),
                                4 * assembler.num_coefficients());
    const Eigen::VectorXd ground_gradient =
        Eigen::VectorXd::Random(assembler.ground_value_size());
    Eigen::VectorXd native_gradient =
        Eigen::VectorXd::Zero(fixture.atmosphere.num_deriv());
    assembler.accumulate_vjp(fixture.atmosphere, wavelength,
                             coefficient_gradient, ground_gradient,
                             native_gradient);
    const double parameter_forward =
        (coefficient_tangent.array() * coefficient_gradient.array()).sum() +
        ground_tangent.dot(ground_gradient);
    REQUIRE(parameter_forward ==
            Catch::Approx(tangent.dot(native_gradient)).epsilon(5.0e-12));
}
