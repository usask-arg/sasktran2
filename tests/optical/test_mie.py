from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
from numpy.testing import assert_allclose


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_variable_radius_mie_matches_fixed_radius(num_stokes):
    """Varying sizes must retain spectral axes and scattering cross sections."""
    z = np.array([0.0, 10000.0, 20000.0])
    radii = np.array([80.0, 150.0, 110.0])
    config = sk.Config()
    config.num_stokes = num_stokes
    config.num_singlescatter_moments = 32
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6372000.0,
        z,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmo = sk.Atmosphere(
        geometry, config, wavelengths_nm=np.array([525.0, 756.0, 1021.0])
    )
    refractive = sk.mie.refractive.RefractiveIndex(
        lambda wavelength: np.full_like(wavelength, 1.45 - 0.01j, dtype=complex),
        "test_absorbing_sulfate",
    )
    optical = sk.optical.Mie(
        sk.mie.LogNormalDistribution().freeze(mode_width=1.6), refractive
    )
    variable = optical.atmosphere_quantities(atmo, median_radius=radii)
    assert variable.leg_coeff.shape == atmo.storage.leg_coeff.shape
    assert (variable.ssa > 0).all()
    assert (variable.ssa < variable.extinction).all()
    xs = optical.cross_sections(atmo.wavelengths_nm, z, median_radius=radii)
    assert_allclose(variable.extinction, xs.extinction, rtol=1e-9)
    assert_allclose(variable.ssa / variable.extinction, xs.ssa, rtol=1e-9)
    for i, radius in enumerate(radii):
        fixed = sk.optical.Mie(
            sk.mie.LogNormalDistribution().freeze(median_radius=radius, mode_width=1.6),
            refractive,
        ).atmosphere_quantities(atmo)
        # Integrating a batch uses a common size quadrature, whereas a frozen
        # distribution uses its own grid; agreement is limited by quadrature.
        assert_allclose(variable.extinction[i], fixed.extinction[i], rtol=1e-4)
        assert_allclose(variable.ssa[i], fixed.ssa[i], rtol=1e-4)
        assert_allclose(
            variable.leg_coeff[:, i], fixed.leg_coeff[:, i], rtol=1e-3, atol=5e-5
        )
