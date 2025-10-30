from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_manual_constituent():
    """
    Verifies that when using the manual constituent, the resulting atmospheric internal
    quantities are identical to what is input
    """
    for stokes in [1, 3]:
        config = sk.Config()
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

        atmosphere = sk.Atmosphere(model_geometry, config, numwavel=30)

        ext = np.zeros_like(atmosphere.storage.total_extinction) + 0.3
        ssa = np.zeros_like(atmosphere.storage.ssa) + 0.4
        leg = np.zeros_like(atmosphere.storage.leg_coeff) + 0.1

        const = sk.constituent.Manual(extinction=ext, ssa=ssa, legendre_moments=leg)

        atmosphere["const"] = const

        # Trigger the internal calculation
        atmosphere.internal_object()

        np.testing.assert_array_almost_equal(
            atmosphere.storage.total_extinction, const.extinction
        )
        np.testing.assert_array_almost_equal(atmosphere.storage.ssa, const.ssa)
        np.testing.assert_array_almost_equal(
            atmosphere.storage.leg_coeff, const.leg_coeff
        )

        np.testing.assert_array_almost_equal(atmosphere.storage.total_extinction, ext)
        np.testing.assert_array_almost_equal(atmosphere.storage.ssa, ssa)
        np.testing.assert_array_almost_equal(atmosphere.storage.leg_coeff, leg)
