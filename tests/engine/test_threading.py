from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
from sasktran2.util import build_info

testdata = [
    (sk.ThreadingLib.Rayon, sk.ThreadingModel.Wavelength),
    (sk.ThreadingLib.OpenMP, sk.ThreadingModel.Wavelength),
    (sk.ThreadingLib.OpenMP, sk.ThreadingModel.Source),
]


@pytest.mark.parametrize(("threading_lib", "threading_model"), testdata)
def test_threading(threading_lib, threading_model):
    if (
        not build_info.openmp_support_enabled()
        and threading_lib == sk.ThreadingLib.OpenMP
    ):
        pytest.skip("OpenMP support is not enabled in this build of sasktran2.")

    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.num_streams = 2
    config.threading_lib = threading_lib
    config.threading_model = threading_model

    model_geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0,
        earth_radius_m=6372000,
        altitude_grid_m=np.arange(0, 65001, 1000.0),
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )

    viewing_geo = sk.ViewingGeometry()

    for alt in [10000, 20000, 30000, 40000]:
        ray = sk.TangentAltitudeSolar(
            tangent_altitude_m=alt,
            relative_azimuth=0,
            observer_altitude_m=200000,
            cos_sza=0.6,
        )
        viewing_geo.add_ray(ray)

    wavel = np.arange(280.0, 800.0, 10)
    atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(),
        model_geometry.altitudes(),
        np.ones_like(model_geometry.altitudes()) * 1e-6,
    )

    engine = sk.Engine(config, model_geometry, viewing_geo)
    rad_ss = engine.calculate_radiance(atmosphere)

    config.num_threads = 2
    engine = sk.Engine(config, model_geometry, viewing_geo)
    rad_threaded = engine.calculate_radiance(atmosphere)
    np.testing.assert_allclose(
        rad_ss["radiance"].to_numpy(),
        rad_threaded["radiance"].to_numpy(),
        rtol=1e-8,
    )


@pytest.mark.parametrize(
    "threading_model",
    [sk.ThreadingModel.Wavelength, sk.ThreadingModel.Source],
)
def test_native_linearization_threading(threading_model):
    if not build_info.openmp_support_enabled():
        pytest.skip("OpenMP support is not enabled in this build of sasktran2.")

    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.threading_lib = sk.ThreadingLib.OpenMP
    config.threading_model = threading_model
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6_372_000.0,
        np.arange(0.0, 30_001.0, 5_000.0),
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    viewing = sk.ViewingGeometry()
    for cos_viewing_zenith in [0.5, 0.65, 0.8, 0.95]:
        viewing.add_ray(sk.GroundViewingSolar(0.6, 0.2, cos_viewing_zenith, 100_000.0))
    atmosphere = sk.Atmosphere(
        geometry, config, wavelengths_nm=np.array([300.0, 350.0, 400.0])
    )
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["surface"] = sk.constituent.LambertianSurface(0.3)

    config.num_threads = 1
    serial = sk.Engine(config, geometry, viewing).linearize(atmosphere)
    tangent = serial.tangent_template[["pressure_pa"]]
    tangent["pressure_pa"].data[:] = np.linspace(0.2, 0.8, atmosphere.num_locations)
    serial_value = serial.value.copy(deep=True)
    serial_jvp = serial.jvp(tangent).copy(deep=True)
    cotangent = serial.value * 0 + 0.7
    serial_vjp = serial.vjp(cotangent)["pressure_pa"].copy(deep=True)

    config.num_threads = 2
    threaded = sk.Engine(config, geometry, viewing).linearize(atmosphere)
    np.testing.assert_allclose(threaded.value, serial_value, rtol=1e-12)
    np.testing.assert_allclose(threaded.jvp(tangent), serial_jvp, rtol=1e-12)
    np.testing.assert_allclose(
        threaded.vjp(cotangent)["pressure_pa"], serial_vjp, rtol=1e-12
    )
