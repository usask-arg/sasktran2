from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import sasktran2 as sk


def test_limb_from_tangent_parameters():
    pytest.importorskip("astropy")
    solar = sk.solar.SolarGeometryHandlerAstropy()

    geo = sk.viewinggeo.LimbVertical.from_tangent_parameters(
        solar,
        np.arange(0, 65000, 1000.0),
        20,
        30,
        pd.Timestamp("2024-11-12 20:00:00"),
        200000,
        0,
    )

    _ = geo.model_geometry(np.arange(0, 65000, 1000.0))


def test_chain_with_viewing_container():
    pytest.importorskip("astropy")
    solar = sk.solar.SolarGeometryHandlerAstropy()

    geo = sk.viewinggeo.LimbVertical.from_tangent_parameters(
        solar,
        np.arange(0, 65000, 1000.0),
        20,
        30,
        pd.Timestamp("2024-11-12 20:00:00"),
        200000,
        0,
    )

    model_geo = geo.model_geometry(np.arange(0, 65000, 1000.0))

    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.num_streams = 2

    wavel = np.arange(280.0, 800.0, 10)
    atmosphere = sk.Atmosphere(model_geo, config, wavelengths_nm=wavel)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    atmosphere["rayleigh"] = sk.constituent.Rayleigh()

    atmosphere["ozone"] = sk.constituent.VMRAltitudeAbsorber(
        sk.optical.O3DBM(),
        model_geo.altitudes(),
        np.ones_like(model_geo.altitudes()) * 1e-6,
    )

    engine = sk.Engine(config, model_geo, geo)

    radiance = engine.calculate_radiance(atmosphere)

    assert "tangent_altitude" in radiance
