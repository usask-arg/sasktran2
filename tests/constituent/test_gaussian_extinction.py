from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def _test_scenarios():
    config = sk.Config()

    altitude_grid = np.arange(0, 65001, 1000.0)

    geometry = sk.Geometry1D(
        0.6,
        0,
        6327000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )

    viewing_geo = sk.ViewingGeometry()

    for tan_alt in np.arange(10000, 60000, 2000):
        viewing_geo.add_ray(sk.TangentAltitudeSolar(tan_alt, 0, 600000, 0.6))

    wavel = np.array([310, 330, 550, 600])

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere.surface.albedo[:] = 0.3

    scen = []

    scen.append(
        {
            "config": config,
            "geometry": geometry,
            "viewing_geo": viewing_geo,
            "atmosphere": atmosphere,
        }
    )

    return scen


def test_extinction_cloud_optical_depth_wf_native_grid():

    scens = _test_scenarios()

    cloud_height_m = 12500
    cloud_width_fwhm = 50000
    vertical_optical_depth = 0.1
    vertical_optical_depth_wavel_nm = 550

    for scen in scens:
        altitude_grid = scen["atmosphere"].model_geometry.altitudes()
        atmosphere = scen["atmosphere"]

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), altitude_grid, np.ones_like(altitude_grid) * 1e-6
        )

        mie = sk.optical.database.OpticalDatabaseGenericScatterer(
            sk.database.StandardDatabase().path("cross_sections/mie/sulfate_test.nc")
        )

        atmosphere["cloud"] = sk.constituent.GaussianHeightExtinction(
            optical_property=mie,
            altitudes_m=altitude_grid,
            cloud_height_m=cloud_height_m,
            cloud_width_fwhm_m=cloud_width_fwhm,
            vertical_optical_depth=vertical_optical_depth,
            vertical_optical_depth_wavel_nm=vertical_optical_depth_wavel_nm,
            lognormal_median_radius=np.ones_like(altitude_grid) * 105.0,
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf_scalar(
            atmosphere["cloud"].vertical_optical_depth,
            0.00001,
            engine,
            atmosphere,
            "wf_cloud_vertical_optical_depth",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_cloud_vertical_optical_depth"],
            radiance["wf_cloud_vertical_optical_depth_numeric"],
            wf_dim="vertical_optical_depth",
            decimal=4,
        )


def test_extinction_cloud_height_wf_native_grid():
    scens = _test_scenarios()

    cloud_height_m = 12500
    cloud_width_fwhm = 50000
    vertical_optical_depth = 0.1
    vertical_optical_depth_wavel_nm = 550

    for scen in scens:
        altitude_grid = scen["atmosphere"].model_geometry.altitudes()
        atmosphere = scen["atmosphere"]

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), altitude_grid, np.ones_like(altitude_grid) * 1e-6
        )

        mie = sk.optical.database.OpticalDatabaseGenericScatterer(
            sk.database.StandardDatabase().path("cross_sections/mie/sulfate_test.nc")
        )

        atmosphere["cloud"] = sk.constituent.GaussianHeightExtinction(
            optical_property=mie,
            altitudes_m=altitude_grid,
            cloud_height_m=cloud_height_m,
            cloud_width_fwhm_m=cloud_width_fwhm,
            vertical_optical_depth=vertical_optical_depth,
            vertical_optical_depth_wavel_nm=vertical_optical_depth_wavel_nm,
            lognormal_median_radius=np.ones_like(altitude_grid) * 105.0,
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf_scalar(
            atmosphere["cloud"].cloud_height_m,
            0.0001,
            engine,
            atmosphere,
            "wf_cloud_height_m",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_cloud_height_m"],
            radiance["wf_cloud_height_m_numeric"],
            wf_dim="cloud_height_m",
            decimal=4,
        )


def test_extinction_cloud_width_wf_native_grid():
    scens = _test_scenarios()

    cloud_height_m = 12500
    cloud_width_fwhm = 50000
    vertical_optical_depth = 0.1
    vertical_optical_depth_wavel_nm = 550

    for scen in scens:
        altitude_grid = scen["atmosphere"].model_geometry.altitudes()
        atmosphere = scen["atmosphere"]

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), altitude_grid, np.ones_like(altitude_grid) * 1e-6
        )

        mie = sk.optical.database.OpticalDatabaseGenericScatterer(
            sk.database.StandardDatabase().path("cross_sections/mie/sulfate_test.nc")
        )

        atmosphere["cloud"] = sk.constituent.GaussianHeightExtinction(
            optical_property=mie,
            altitudes_m=altitude_grid,
            cloud_height_m=cloud_height_m,
            cloud_width_fwhm_m=cloud_width_fwhm,
            vertical_optical_depth=vertical_optical_depth,
            vertical_optical_depth_wavel_nm=vertical_optical_depth_wavel_nm,
            lognormal_median_radius=np.ones_like(altitude_grid) * 105.0,
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf_scalar(
            atmosphere["cloud"].cloud_width_fwhm_m,
            0.00001,
            engine,
            atmosphere,
            "wf_cloud_width_fwhm_m",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_cloud_width_fwhm_m"],
            radiance["wf_cloud_width_fwhm_m_numeric"],
            wf_dim="cloud_width_fwhm_m",
            decimal=4,
        )


def test_extinction_cloud_radius_wf_native_grid():
    scens = _test_scenarios()

    cloud_height_m = 12500
    cloud_width_fwhm = 50000
    vertical_optical_depth = 0.1
    vertical_optical_depth_wavel_nm = 550

    for scen in scens:
        altitude_grid = scen["atmosphere"].model_geometry.altitudes()
        atmosphere = scen["atmosphere"]

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), altitude_grid, np.ones_like(altitude_grid) * 1e-6
        )

        mie = sk.optical.database.OpticalDatabaseGenericScatterer(
            sk.database.StandardDatabase().path("cross_sections/mie/sulfate_test.nc")
        )

        atmosphere["cloud"] = sk.constituent.GaussianHeightExtinction(
            optical_property=mie,
            altitudes_m=altitude_grid,
            cloud_height_m=cloud_height_m,
            cloud_width_fwhm_m=cloud_width_fwhm,
            vertical_optical_depth=vertical_optical_depth,
            vertical_optical_depth_wavel_nm=vertical_optical_depth_wavel_nm,
            lognormal_median_radius=np.ones_like(altitude_grid) * 105.0,
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere["cloud"].lognormal_median_radius,
            0.00001,
            engine,
            atmosphere,
            "wf_cloud_lognormal_median_radius",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_cloud_lognormal_median_radius"],
            radiance["wf_cloud_lognormal_median_radius_numeric"],
            wf_dim="cloud_altitude",
            decimal=4,
        )


def test_extinction_cloud_optical_depth_wf_interpolated_grid():

    scens = _test_scenarios()
    new_altitudes = np.array([0, 10000, 30000, 70000])

    cloud_height_m = 12500
    cloud_width_fwhm = 50000
    vertical_optical_depth = 0.1
    vertical_optical_depth_wavel_nm = 550

    for scen in scens:
        altitude_grid = scen["atmosphere"].model_geometry.altitudes()
        atmosphere = scen["atmosphere"]

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), altitude_grid, np.ones_like(altitude_grid) * 1e-6
        )

        mie = sk.optical.database.OpticalDatabaseGenericScatterer(
            sk.database.StandardDatabase().path("cross_sections/mie/sulfate_test.nc")
        )

        atmosphere["cloud"] = sk.constituent.GaussianHeightExtinction(
            optical_property=mie,
            altitudes_m=new_altitudes,
            cloud_height_m=cloud_height_m,
            cloud_width_fwhm_m=cloud_width_fwhm,
            vertical_optical_depth=vertical_optical_depth,
            vertical_optical_depth_wavel_nm=vertical_optical_depth_wavel_nm,
            lognormal_median_radius=np.ones_like(new_altitudes) * 105.0,
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf_scalar(
            atmosphere["cloud"].vertical_optical_depth,
            0.00001,
            engine,
            atmosphere,
            "wf_cloud_vertical_optical_depth",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_cloud_vertical_optical_depth"],
            radiance["wf_cloud_vertical_optical_depth_numeric"],
            wf_dim="vertical_optical_depth",
            decimal=4,
        )


def test_extinction_cloud_height_wf_interpolated_grid():
    scens = _test_scenarios()
    new_altitudes = np.array([0, 10000, 30000, 70000])

    cloud_height_m = 12500
    cloud_width_fwhm = 50000
    vertical_optical_depth = 0.1
    vertical_optical_depth_wavel_nm = 550

    for scen in scens:
        altitude_grid = scen["atmosphere"].model_geometry.altitudes()
        atmosphere = scen["atmosphere"]

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), altitude_grid, np.ones_like(altitude_grid) * 1e-6
        )

        mie = sk.optical.database.OpticalDatabaseGenericScatterer(
            sk.database.StandardDatabase().path("cross_sections/mie/sulfate_test.nc")
        )

        atmosphere["cloud"] = sk.constituent.GaussianHeightExtinction(
            optical_property=mie,
            altitudes_m=new_altitudes,
            cloud_height_m=cloud_height_m,
            cloud_width_fwhm_m=cloud_width_fwhm,
            vertical_optical_depth=vertical_optical_depth,
            vertical_optical_depth_wavel_nm=vertical_optical_depth_wavel_nm,
            lognormal_median_radius=np.ones_like(new_altitudes) * 105.0,
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf_scalar(
            atmosphere["cloud"].cloud_height_m,
            0.001,
            engine,
            atmosphere,
            "wf_cloud_height_m",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_cloud_height_m"],
            radiance["wf_cloud_height_m_numeric"],
            wf_dim="cloud_height_m",
            decimal=4,
        )


def test_extinction_cloud_width_wf_interpolated_grid():
    scens = _test_scenarios()
    new_altitudes = np.array([0, 10000, 30000, 70000])

    cloud_height_m = 12500
    cloud_width_fwhm = 50000
    vertical_optical_depth = 0.1
    vertical_optical_depth_wavel_nm = 550

    for scen in scens:
        altitude_grid = scen["atmosphere"].model_geometry.altitudes()
        atmosphere = scen["atmosphere"]

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), altitude_grid, np.ones_like(altitude_grid) * 1e-6
        )

        mie = sk.optical.database.OpticalDatabaseGenericScatterer(
            sk.database.StandardDatabase().path("cross_sections/mie/sulfate_test.nc")
        )

        atmosphere["cloud"] = sk.constituent.GaussianHeightExtinction(
            optical_property=mie,
            altitudes_m=new_altitudes,
            cloud_height_m=cloud_height_m,
            cloud_width_fwhm_m=cloud_width_fwhm,
            vertical_optical_depth=vertical_optical_depth,
            vertical_optical_depth_wavel_nm=vertical_optical_depth_wavel_nm,
            lognormal_median_radius=np.ones_like(new_altitudes) * 105.0,
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf_scalar(
            atmosphere["cloud"].cloud_width_fwhm_m,
            0.00001,
            engine,
            atmosphere,
            "wf_cloud_width_fwhm_m",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_cloud_width_fwhm_m"],
            radiance["wf_cloud_width_fwhm_m_numeric"],
            wf_dim="cloud_width_fwhm_m",
            decimal=4,
        )


@pytest.mark.skip()
def test_extinction_cloud_radius_wf_interpolated_grid():
    scens = _test_scenarios()
    new_altitudes = np.array([0, 10000, 30000, 70000])

    cloud_height_m = 12500
    cloud_width_fwhm = 50000
    vertical_optical_depth = 0.1
    vertical_optical_depth_wavel_nm = 550

    for scen in scens:
        altitude_grid = scen["atmosphere"].model_geometry.altitudes()
        atmosphere = scen["atmosphere"]

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(), altitude_grid, np.ones_like(altitude_grid) * 1e-6
        )

        mie = sk.optical.database.OpticalDatabaseGenericScatterer(
            sk.database.StandardDatabase().path("cross_sections/mie/sulfate_test.nc")
        )

        atmosphere["cloud"] = sk.constituent.GaussianHeightExtinction(
            optical_property=mie,
            altitudes_m=new_altitudes,
            cloud_height_m=cloud_height_m,
            cloud_width_fwhm_m=cloud_width_fwhm,
            vertical_optical_depth=vertical_optical_depth,
            vertical_optical_depth_wavel_nm=vertical_optical_depth_wavel_nm,
            lognormal_median_radius=np.ones_like(new_altitudes) * 105.0,
        )

        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = sk.test_util.wf.numeric_wf(
            atmosphere["cloud"].lognormal_median_radius,
            0.00001,
            engine,
            atmosphere,
            "wf_cloud_lognormal_median_radius",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_cloud_lognormal_median_radius"],
            radiance["wf_cloud_lognormal_median_radius_numeric"],
            wf_dim="cloud_altitude",
            decimal=2,
        )
