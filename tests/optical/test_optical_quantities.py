from __future__ import annotations

import numpy as np
import sasktran2 as sk
from sasktran2._core_rust import PyScatteringDatabaseDim2


def test_native_derivative_legendre_alias_uses_python_axis_order():
    wavenumbers = np.array([10000.0, 11000.0, 12000.0])
    asymmetry = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
    database = PyScatteringDatabaseDim2.from_asymmetry_parameter(
        np.ones_like(asymmetry),
        np.ones_like(asymmetry),
        asymmetry,
        4,
        wavenumbers,
        np.array([0.0, 1.0]),
        ["radius"],
    )
    config = sk.Config()
    config.num_singlescatter_moments = 4
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6372000.0,
        np.array([0.0, 1000.0]),
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(geometry, config, wavenumber_cminv=wavenumbers)
    atmosphere.temperature_k = np.full(2, 250.0)
    atmosphere.pressure_pa = np.full(2, 50000.0)
    derivative = database.optical_derivatives(atmosphere, radius=np.array([0.3, 0.7]))[
        "radius"
    ]
    coefficients = derivative.d_leg_coeff
    assert coefficients.shape == (4, 2, 3)
    expected = np.array(
        [
            (2 * order + 1) * (asymmetry[1] ** order - asymmetry[0] ** order)
            for order in range(4)
        ]
    )
    np.testing.assert_allclose(
        coefficients, np.broadcast_to(expected[:, None, :], (4, 2, 3))
    )
    assert np.shares_memory(coefficients, derivative.leg_coeff)
    # The NumPy view must retain its owner after the derivative wrapper is gone.
    del derivative
    np.testing.assert_allclose(coefficients[:, 0, :], expected)
