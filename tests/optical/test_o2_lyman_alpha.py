from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def _test_atmosphere(wavelengths_nm: np.ndarray) -> sk.Atmosphere:
    config = sk.Config()
    altitude_grid_m = np.array([0.0, 50_000.0])
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6_372_000.0,
        altitude_grid_m,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.PlaneParallel,
    )

    return sk.Atmosphere(geometry, config, wavelengths_nm=wavelengths_nm)


def test_o2_lyman_alpha_cross_section_database_shape_and_values():
    optical = sk.optical.O2LymanAlpha()
    center = optical.WAVELENGTH_NM
    half_width = 5.0e-4
    wavelengths_nm = np.array(
        [
            center - half_width,
            center,
            center + half_width,
            center + 2.0 * half_width,
        ]
    )
    atmosphere = _test_atmosphere(wavelengths_nm)

    quantities = optical.atmosphere_quantities(atmosphere)

    assert quantities.extinction.shape == (2, 4)
    np.testing.assert_allclose(quantities.extinction[:, 0], 0.0, atol=0.0)
    np.testing.assert_allclose(
        quantities.extinction[:, 1],
        optical.EFFECTIVE_CROSS_SECTION_M2,
        rtol=1e-12,
        atol=0.0,
    )
    np.testing.assert_allclose(quantities.extinction[:, 2:], 0.0, atol=0.0)


def test_o2_lyman_alpha_effective_cross_section_matches_toa_rate():
    optical = sk.optical.O2LymanAlpha()

    rate = optical.EFFECTIVE_CROSS_SECTION_M2 * optical.TOA_FLUX_PHOTONS_M2_S

    np.testing.assert_allclose(rate, optical.TOA_RATE_S, rtol=1e-15, atol=0.0)
    np.testing.assert_allclose(
        optical.EFFECTIVE_CROSS_SECTION_M2 * 1.0e4,
        1.0625e-20,
        rtol=1e-15,
        atol=0.0,
    )


def test_o2_lyman_alpha_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="cross_section_m2 must be non-negative"):
        sk.optical.O2LymanAlpha(cross_section_m2=-1.0)

    with pytest.raises(ValueError, match="half_width_nm must be positive"):
        sk.optical.O2LymanAlpha(half_width_nm=0.0)
