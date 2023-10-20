import numpy as np
import sasktran2 as sk


def test_scattering_db_construction():
    mie = sk.optical.database.OpticalDatabaseGenericScatterer(
        sk.appconfig.database_root().joinpath("cross_sections/mie/sulfate_test.nc")
    )

    config = sk.Config()

    altitude_grid = np.arange(0, 65001, 1000)

    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )

    wavel = np.arange(350, 500, 0.1)

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    _ = mie.atmosphere_quantities(
        atmo=atmosphere,
        lognormal_median_radius=np.ones_like(atmosphere.model_geometry.altitudes())
        * 105,
    ).extinction[10, :]
