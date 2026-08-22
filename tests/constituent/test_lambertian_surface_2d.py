from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def geometry2d() -> sk.Geometry2D:
    return sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=np.array([0.0, 10_000.0, 30_000.0]),
        horizontal_angle_grid_radians=np.array([-0.3, 0.3]),
    )


def orbital_geometry() -> sk.OrbitalPlaneGeometry:
    angles = np.array([-0.3, 0.3])
    track = 6_372_000.0 * np.column_stack(
        (np.sin(angles), np.zeros_like(angles), np.cos(angles))
    )
    return sk.OrbitalPlaneGeometry(
        6_372_000.0,
        np.array([0.0, 10_000.0, 30_000.0]),
        track,
    )


def atmosphere(geometry: sk.Geometry2D) -> sk.Atmosphere:
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    return sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([400.0, 500.0, 600.0]),
        legendre_derivative=False,
    )


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"albedo": np.array([[-0.1], [0.2]])}, "between zero and one"),
        ({"albedo": np.array([[np.nan], [0.2]])}, "finite"),
        (
            {
                "albedo": np.ones((2, 2)),
                "wavelengths_nm": np.array([500.0, 400.0]),
            },
            "strictly increasing",
        ),
        (
            {
                "albedo": np.ones((2, 2)),
                "wavelengths_nm": np.array([0.0, 500.0]),
            },
            "positive",
        ),
        (
            {
                "albedo": np.ones((2, 3)),
                "wavelengths_nm": np.array([400.0, 500.0]),
            },
            "spectral-grid length",
        ),
        ({"albedo": 0.2, "out_of_bounds_mode": "invalid"}, "out_of_bounds_mode"),
    ],
)
def test_lambertian_surface_2d_rejects_invalid_constructor_inputs(kwargs, match):
    with pytest.raises(ValueError, match=match):
        sk.constituent.LambertianSurface2D(**kwargs)


def test_lambertian_surface_2d_requires_orbital_geometry_not_generic_geometry2d():
    atmo = atmosphere(geometry2d())
    atmo["surface"] = sk.constituent.LambertianSurface2D(np.array([0.2, 0.3]))
    with pytest.raises(TypeError, match="OrbitalPlaneGeometry"):
        atmo.internal_object()


def test_lambertian_surface_2d_interpolates_spatial_spectra():
    geometry = orbital_geometry()
    atmo = atmosphere(geometry)
    surface = sk.constituent.LambertianSurface2D(
        np.array([[0.1, 0.3], [0.5, 0.7]]),
        wavelengths_nm=np.array([400.0, 600.0]),
    )
    atmo["surface"] = surface

    atmo.internal_object()

    np.testing.assert_allclose(
        atmo._orbital_lambertian_surface,
        np.array([[0.1, 0.2, 0.3], [0.5, 0.6, 0.7]]),
    )


def test_lambertian_surface_2d_reports_missing_atmosphere_spectral_coordinates():
    geometry = orbital_geometry()
    config = sk.Config()
    atmo = sk.Atmosphere(geometry, config, numwavel=2, legendre_derivative=False)
    atmo["surface"] = sk.constituent.LambertianSurface2D(
        np.array([[0.1, 0.2], [0.3, 0.4]]),
        wavelengths_nm=np.array([400.0, 500.0]),
    )

    with pytest.raises(ValueError, match="Atmosphere constructed with"):
        atmo.internal_object()


def test_lambertian_surface_2d_uses_explicit_spectral_point_layout_without_grid():
    geometry = orbital_geometry()
    config = sk.Config()
    atmo = sk.Atmosphere(geometry, config, numwavel=2, legendre_derivative=False)
    atmo["surface"] = sk.constituent.LambertianSurface2D(
        np.array([[0.1, 0.2], [0.3, 0.4]])
    )

    atmo.internal_object()
    layout = atmo.surface_derivative_output_layout("wf_surface_albedo")

    assert layout is not None
    assert layout[0] == ("orbital_position", "surface_spectral_point")
    np.testing.assert_array_equal(layout[2]["surface_spectral_point"], [0, 1])


def test_lambertian_surface_2d_labels_wavenumber_input_as_wavenumber():
    geometry = orbital_geometry()
    config = sk.Config()
    atmo = sk.Atmosphere(
        geometry,
        config,
        wavenumber_cminv=np.array([10_000.0, 20_000.0]),
        legendre_derivative=False,
    )
    atmo["surface"] = sk.constituent.LambertianSurface2D(
        np.array([[0.1, 0.2], [0.3, 0.4]]),
        wavenumbers_cminv=np.array([10_000.0, 20_000.0]),
    )

    atmo.internal_object()
    layout = atmo.surface_derivative_output_layout("wf_surface_albedo")

    assert layout is not None
    assert layout[0] == ("orbital_position", "surface_wavenumber")
    np.testing.assert_array_equal(layout[2]["surface_wavenumber"], [10_000.0, 20_000.0])


def test_lambertian_surface_2d_rechecks_mutable_albedo_and_clears_when_replaced():
    geometry = orbital_geometry()
    atmo = atmosphere(geometry)
    surface = sk.constituent.LambertianSurface2D(np.array([0.2, 0.3]))
    atmo["surface"] = surface
    atmo.internal_object()
    assert atmo._orbital_lambertian_surface is not None

    surface.albedo[0] = 2.0
    with pytest.raises(ValueError, match="between zero and one"):
        atmo.internal_object()

    surface.albedo[0] = 0.2
    atmo["surface"] = sk.constituent.LambertianSurface(0.4)
    atmo.internal_object()
    assert atmo._orbital_lambertian_surface is None
