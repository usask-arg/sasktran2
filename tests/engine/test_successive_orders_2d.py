from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr

EARTH_RADIUS_M = 6_372_000.0
ALTITUDES_M = np.array([0.0, 10_000.0, 30_000.0])
HORIZONTAL_ANGLES_RAD = np.array([-0.6, 0.0, 0.6])


def successive_orders_config(
    *,
    num_stokes: int = 1,
    single_scatter_source: sk.SingleScatterSource = sk.SingleScatterSource.Exact,
    multiple_scatter_source: sk.MultipleScatterSource = (
        sk.MultipleScatterSource.SuccessiveOrders
    ),
) -> sk.Config:
    config = sk.Config()
    config.num_threads = 1
    config.num_stokes = num_stokes
    config.num_streams = 4
    config.num_singlescatter_moments = 4
    config.single_scatter_source = single_scatter_source
    config.multiple_scatter_source = multiple_scatter_source
    config.occultation_source = sk.OccultationSource.NoSource
    config.emission_source = sk.EmissionSource.NoSource
    config.num_successive_orders_incoming = 6
    config.num_successive_orders_outgoing = 6
    config.num_successive_orders_iterations = 2
    config.successive_orders_relative_tolerance = 0.0
    config.successive_orders_absolute_tolerance = 0.0
    config.delta_m_scaling = False
    return config


def geometry2d() -> sk.Geometry2D:
    return sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=EARTH_RADIUS_M,
        altitude_grid_m=ALTITUDES_M,
        horizontal_angle_grid_radians=HORIZONTAL_ANGLES_RAD,
    )


def viewing_geometry() -> sk.ViewingGeometry:
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.GroundViewingSolar(
            cos_sza=0.6,
            relative_azimuth=0.2,
            cos_viewing_zenith=0.7,
            observer_altitude_m=100_000.0,
        )
    )
    return viewing


def atmosphere(
    geometry: sk.Geometry2D,
    config: sk.Config,
    *,
    horizontal_slope: float = 0.0,
    calculate_derivatives: bool = False,
) -> sk.Atmosphere:
    result = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
        calculate_derivatives=calculate_derivatives,
    )
    horizontal, altitude = np.meshgrid(
        HORIZONTAL_ANGLES_RAD, ALTITUDES_M, indexing="ij"
    )
    result.storage.total_extinction[:, 0] = (
        1.5e-5 * np.exp(-altitude / 15_000.0) * (1.0 + horizontal_slope * horizontal)
    ).ravel()
    result.storage.ssa[:, 0] = (0.88 - 0.08 * altitude / ALTITUDES_M[-1]).ravel()
    result.leg_coeff.a1[0] = 1.0
    result.leg_coeff.a1[2] = 0.3
    if config.num_stokes == 3:
        result.leg_coeff.a2[2] = 3.0
        result.leg_coeff.b1[2] = -np.sqrt(6.0) / 2.0
    result.surface.albedo[:] = 0.1
    return result


def calculate(
    geometry: sk.Geometry2D,
    *,
    num_stokes: int = 1,
    horizontal_slope: float = 0.0,
    single_source: sk.SingleScatterSource = sk.SingleScatterSource.Exact,
    source: sk.MultipleScatterSource = sk.MultipleScatterSource.SuccessiveOrders,
) -> np.ndarray:
    config = successive_orders_config(
        num_stokes=num_stokes,
        single_scatter_source=single_source,
        multiple_scatter_source=source,
    )
    result = sk.Engine(config, geometry, viewing_geometry()).calculate_radiance(
        atmosphere(geometry, config, horizontal_slope=horizontal_slope)
    )
    return result.radiance.values


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_2d_successive_orders_is_finite_and_adds_multiple_scatter(num_stokes: int):
    geometry = geometry2d()
    multiple = calculate(
        geometry,
        num_stokes=num_stokes,
        single_source=sk.SingleScatterSource.NoSource,
    )

    assert multiple.shape == (1, 1, num_stokes)
    assert np.all(np.isfinite(multiple))
    assert multiple[0, 0, 0] > 0.0
    if num_stokes == 1:
        combined = calculate(geometry, num_stokes=num_stokes)
        assert np.all(np.isfinite(combined))
        assert combined[0, 0, 0] > multiple[0, 0, 0]


def test_2d_successive_orders_uses_horizontal_atmospheric_structure():
    geometry = geometry2d()
    config = successive_orders_config(
        single_scatter_source=sk.SingleScatterSource.NoSource
    )
    engine = sk.Engine(config, geometry, viewing_geometry())
    uniform_multiple = engine.calculate_radiance(
        atmosphere(geometry, config)
    ).radiance.values
    varying_multiple = engine.calculate_radiance(
        atmosphere(geometry, config, horizontal_slope=0.8)
    ).radiance.values

    assert not np.isclose(
        varying_multiple.item(), uniform_multiple.item(), rtol=1.0e-4, atol=0.0
    )


def test_2d_successive_orders_native_products_are_adjoint():
    geometry = geometry2d()
    config = successive_orders_config(
        single_scatter_source=sk.SingleScatterSource.NoSource
    )
    config.num_successive_orders_iterations = 60
    config.successive_orders_relative_tolerance = 1.0e-11
    config.successive_orders_absolute_tolerance = 1.0e-13
    config.successive_orders_anderson_depth = 3
    linearization = sk.Engine(config, geometry, viewing_geometry()).linearize(
        atmosphere(geometry, config, horizontal_slope=0.4, calculate_derivatives=True)
    )

    assert linearization.backends == {
        "jvp": sk.LinearizationBackend.Native,
        "vjp": sk.LinearizationBackend.Native,
    }
    tangent = linearization.tangent_template[["extinction", "ssa"]]
    tangent["extinction"].data[:] = np.linspace(
        -2.0e-7, 3.0e-7, tangent["extinction"].size
    ).reshape(tangent["extinction"].shape)
    tangent["ssa"].data[:] = np.linspace(-0.015, 0.02, tangent["ssa"].size).reshape(
        tangent["ssa"].shape
    )
    jvp = linearization.jvp(tangent)

    cotangent = xr.ones_like(linearization.value)
    gradient = linearization.vjp(cotangent, parameters=("extinction", "ssa"))
    lhs = float((jvp * cotangent).sum())
    rhs = float(
        (tangent["extinction"] * gradient["extinction"]).sum()
        + (tangent["ssa"] * gradient["ssa"]).sum()
    )

    assert np.all(np.isfinite(jvp))
    assert np.all(np.isfinite(gradient["extinction"]))
    assert np.all(np.isfinite(gradient["ssa"]))
    np.testing.assert_allclose(lhs, rhs, rtol=3.0e-8, atol=3.0e-11)


def test_2d_successive_orders_rejects_diffuse_refraction():
    config = successive_orders_config()
    config.multiple_scatter_refraction = True

    with pytest.raises(NotImplementedError, match="diffuse-ray refraction"):
        sk.Engine(config, geometry2d(), viewing_geometry())


def test_2d_rejects_legacy_successive_orders_source():
    config = successive_orders_config(
        multiple_scatter_source=sk.MultipleScatterSource.SuccessiveOrdersLegacy
    )

    with pytest.raises(NotImplementedError, match="successive-orders"):
        sk.Engine(config, geometry2d(), viewing_geometry())
