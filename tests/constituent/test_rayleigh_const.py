from __future__ import annotations

import numpy as np
import sasktran2 as sk
from sasktran2._core_rust import Rayleigh as RustRayleigh
from sasktran2.constituent.rayleigh import Rayleigh as PyRayleigh


def test_bates_identical():
    """
    Verifies that the Rust Rayleigh Constituent gives identical radiances and
    weighting functions compared to the Python Rayleigh constituent when using
    the Bates parameterization.
    """
    for stokes in [1, 3]:
        config = sk.Config()
        config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
        config.num_streams = 4
        config.num_stokes = stokes

        model_geometry = sk.Geometry1D(
            cos_sza=0.6,
            solar_azimuth=0,
            earth_radius_m=6372000,
            altitude_grid_m=np.arange(0, 65001, 1000),
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

        atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
            sk.optical.O3DBM(),
            model_geometry.altitudes(),
            np.ones_like(model_geometry.altitudes()) * 1e-6,
        )

        engine = sk.Engine(config, model_geometry, viewing_geo)

        atmosphere["rayleigh"] = PyRayleigh()
        rad_py = engine.calculate_radiance(atmosphere)

        atmosphere["rayleigh"] = RustRayleigh()
        rad_rust = engine.calculate_radiance(atmosphere)

        # Verify the radiances and weighting functions are identical
        for var in rad_py:
            np.testing.assert_array_almost_equal(
                (rad_py[var] / rad_rust[var]).to_numpy(), 1.0, decimal=8
            )
