from __future__ import annotations

import numpy as np
import sasktran2 as sk


def test_air_vacuum_conversions():
    """
    Verifies that if we convert to vacuum and back, we get the same thing
    """
    wavel = np.arange(280, 10000, 0.1)

    conv_back = sk.optical.air_wavelength_to_vacuum_wavelength(
        sk.optical.vacuum_wavelength_to_air_wavelength(wavel)
    )

    np.testing.assert_allclose(wavel, conv_back)
