from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_disort_thermal_only_example():
    """
    This is an example generated from DISORT test case 7a (thermal emission only) with a few
    modifications
    """
    od = 100.0
    ssa = 0.95
    g = 0.75

    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    config.emission_source = sk.EmissionSource.DiscreteOrdinates
    config.num_streams = 16
    config.num_singlescatter_moments = 17

    model_geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0,
        earth_radius_m=6372000,
        altitude_grid_m=np.array([0, 1000.0]),
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.PlaneParallel,
    )

    viewing_geo = sk.ViewingGeometry()
    viewing_geo.add_ray(
        sk.GroundViewingSolar(
            0.6, 0.0, 1.0, 200000.0
        )
    )

    atmosphere = sk.Atmosphere(model_geometry, config, numwavel=1, calculate_derivatives=False)

    atmosphere.storage.total_extinction[:] = od / 1000.0
    atmosphere.storage.ssa[:] = ssa
    atmosphere.storage.solar_irradiance[:] = 0.0
    atmosphere.storage.emission_source[:] = 1.09657540E-05

    for l in range(17):
        coeff = g**l * (2 * l + 1)
        atmosphere.leg_coeff.a1[l][:] = coeff

    atmosphere._atmosphere.apply_delta_m_scaling(config.num_streams)

    engine = sk.Engine(config, model_geometry, viewing_geo)

    rad = engine.calculate_radiance(atmosphere)

    np.testing.assert_allclose(rad["radiance"].values[0, 0, 0], 7.93075833E-06, rtol=1e-6)
