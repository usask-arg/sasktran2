#include <sasktran2.h>
#include <sasktran2/test_helper.h>

TEST_CASE("Batched native C++ outputs match scalar assignment",
          "[sasktran2][output]") {
    constexpr int nwavel = 5;
    constexpr int nalt = 21;

    sasktran2::Coordinates coordinates(0.6, 0.2, 6372000,
                                       sasktran2::geometrytype::spherical);
    Eigen::VectorXd altitude_values(nalt);
    for (int altitude = 0; altitude < nalt; ++altitude) {
        altitude_values(altitude) = altitude * 5000.0;
    }
    sasktran2::grids::AltitudeGrid altitude_grid(
        std::move(altitude_values), sasktran2::grids::gridspacing::constant,
        sasktran2::grids::outofbounds::extend,
        sasktran2::grids::interpolation::linear);
    sasktran2::Geometry1D geometry(std::move(coordinates),
                                  std::move(altitude_grid));

    sasktran2::viewinggeometry::ViewingGeometryContainer viewing_geometry;
    viewing_geometry.observer_rays().emplace_back(
        std::make_unique<sasktran2::viewinggeometry::GroundViewingSolar>(
            0.6, -0.4, 0.5, 200000.0));
    viewing_geometry.observer_rays().emplace_back(
        std::make_unique<sasktran2::viewinggeometry::TangentAltitudeSolar>(
            17500.0, 0.3, 200000.0, 0.6));

    sasktran2::Config scalar_config;
    sasktran2::Config batch_config = scalar_config;
    batch_config.set_wavelength_batch_size(4);

    sasktran2::atmosphere::Atmosphere<1> atmosphere(
        nwavel, geometry, scalar_config, true);
    for (int wavelength = 0; wavelength < nwavel; ++wavelength) {
        for (int altitude = 0; altitude < nalt; ++altitude) {
            atmosphere.storage().total_extinction(altitude, wavelength) =
                2.0e-5 * std::exp(-altitude / 4.0) *
                    (1.0 + 0.1 * wavelength) +
                1.0e-10;
            atmosphere.storage().ssa(altitude, wavelength) =
                0.92 - 0.01 * wavelength;
        }
    }
    atmosphere.storage().leg_coeff.chip(0, 0).setConstant(1.0);
    atmosphere.storage().leg_coeff.chip(2, 0).setConstant(0.5);
    atmosphere.surface().brdf_args().row(0).setConstant(0.2);

    auto& mapping =
        atmosphere.storage().get_derivative_mapping("constituent");
    mapping.allocate_extinction_derivatives();
    mapping.allocate_ssa_derivatives();
    mapping.native_mapping().d_extinction->setConstant(0.7);
    mapping.native_mapping().d_ssa->setConstant(0.3);

    Sasktran2<1> scalar_engine(scalar_config, &geometry, viewing_geometry);
    Sasktran2<1> batch_engine(batch_config, &geometry, viewing_geometry);

    sasktran2::OutputIdealDense<1> scalar_dense;
    sasktran2::OutputIdealDense<1> batch_dense;
    scalar_engine.calculate_radiance(atmosphere, scalar_dense);
    batch_engine.calculate_radiance(atmosphere, batch_dense);
    REQUIRE(batch_dense.radiance().value.isApprox(
        scalar_dense.radiance().value, 1e-12));
    REQUIRE(batch_dense.radiance().deriv.isApprox(
        scalar_dense.radiance().deriv, 1e-12));

    sasktran2::OutputDerivMapped<1> scalar_mapped;
    sasktran2::OutputDerivMapped<1> batch_mapped;
    scalar_engine.calculate_radiance(atmosphere, scalar_mapped);
    batch_engine.calculate_radiance(atmosphere, batch_mapped);
    REQUIRE(batch_mapped.radiance().value.isApprox(
        scalar_mapped.radiance().value, 1e-12));
    REQUIRE(batch_mapped.derivatives().at("constituent")
                .isApprox(scalar_mapped.derivatives().at("constituent"),
                          1e-12));
}
