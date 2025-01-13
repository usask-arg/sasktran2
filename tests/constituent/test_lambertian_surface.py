from __future__ import annotations

import numpy as np
import sasktran2 as sk


def _test_scenarios():
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.num_streams = 4
    config.delta_m_scaling = False

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

    wavel = np.array([310, 330, 350, 600])

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

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


def test_lambertian_surface_constant():
    """
    Checks that we can create a lambertian surface that is constant in wavelength
    """
    scens = _test_scenarios()

    for scen in scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        scen["atmosphere"]["surface"] = sk.constituent.LambertianSurface(0.3)

        _ = engine.calculate_radiance(scen["atmosphere"])

        np.testing.assert_allclose(scen["atmosphere"].surface.albedo, 0.3)


def test_lambertian_surface_interp_wavelength():
    """
    Checks that we can create a lambertian surface that is interpolated in wavelength
    """
    scens = _test_scenarios()

    for scen in scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        scen["atmosphere"]["surface"] = sk.constituent.LambertianSurface(
            albedo=[0.1, 0.3], wavelengths_nm=[280, 600]
        )

        _ = engine.calculate_radiance(scen["atmosphere"])

        np.testing.assert_allclose(
            scen["atmosphere"].surface.albedo.flatten(),
            [0.11875, 0.13125, 0.14375, 0.3],
        )


def test_lambertian_surface_interp_wavenumber():
    """
    Checks that we can create a lambertian surface that is interpolated in wavenumber
    """
    scens = _test_scenarios()

    for scen in scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        scen["atmosphere"]["surface"] = sk.constituent.LambertianSurface(
            albedo=[0.1, 0.3], wavenumbers_cminv=[10000, 30000]
        )

        _ = engine.calculate_radiance(scen["atmosphere"])

        np.testing.assert_allclose(
            scen["atmosphere"].surface.albedo.flatten(), [0, 0, 0.28571429, 0.16666667]
        )


def test_lambertian_surface_native():
    """
    Checks that we can create a lambertian surface that is on the same grid as the atmospheric wavelength
    """
    scens = _test_scenarios()

    for scen in scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        scen["atmosphere"]["surface"] = sk.constituent.LambertianSurface(
            albedo=[0.1, 0.8, 0.9, 0.3]
        )

        _ = engine.calculate_radiance(scen["atmosphere"])

        np.testing.assert_allclose(
            scen["atmosphere"].surface.albedo.flatten(), [0.1, 0.8, 0.9, 0.3]
        )


def test_lambertian_surface_wf_native():
    """
    Checks that the lambertian surface weighting functions are correct when specified on the native atmospheric wavelength grid
    """
    scens = _test_scenarios()

    for scen in scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        scen["atmosphere"]["surface"] = sk.constituent.LambertianSurface(
            albedo=[0.1, 0.8, 0.9, 0.3]
        )

        radiance = sk.test_util.wf.numeric_wf(
            scen["atmosphere"]["surface"].albedo,
            0.001,
            engine,
            scen["atmosphere"],
            "wf_surface_albedo",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_surface_albedo"],
            radiance["wf_surface_albedo_numeric"],
            "surface_wavelength",
            decimal=6,
        )


def test_lambertian_surface_wf_scalar():
    """
    Checks that the lambertian surface weighting functions are correct for a scalar albedo
    """
    scens = _test_scenarios()

    for scen in scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        scen["atmosphere"]["surface"] = sk.constituent.LambertianSurface(albedo=0.3)

        radiance = sk.test_util.wf.numeric_wf(
            scen["atmosphere"]["surface"].albedo,
            0.001,
            engine,
            scen["atmosphere"],
            "wf_surface_albedo",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_surface_albedo"],
            radiance["wf_surface_albedo_numeric"],
            "surface_wavelength",
            decimal=6,
        )


def test_lambertian_surface_wf_interpolated():
    """
    Checks that the lambertian surface weighting functions are correct when we are interpolating over wavelength
    """
    scens = _test_scenarios()

    for scen in scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        scen["atmosphere"]["surface"] = sk.constituent.LambertianSurface(
            albedo=[0.1, 0.3], wavelengths_nm=[280, 600]
        )

        radiance = sk.test_util.wf.numeric_wf(
            scen["atmosphere"]["surface"].albedo,
            0.001,
            engine,
            scen["atmosphere"],
            "wf_surface_albedo",
        )

        sk.test_util.wf.validate_wf(
            radiance["wf_surface_albedo"],
            radiance["wf_surface_albedo_numeric"],
            "surface_wavelength",
            decimal=6,
        )
