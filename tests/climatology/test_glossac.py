from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


@pytest.mark.skip()
def test_glossac_loading():
    """
    Verifies that the glossac data can be loaded in correctly
    """
    _ = sk.climatology.glossac.load_glossac_raw_data()


@pytest.mark.skip()
def test_stratospheric_background():
    """
    Checks the stratospheric background pulling for a variety of conditions
    """
    _ = sk.climatology.glossac.stratospheric_background(
        1, 0.0, np.arange(0, 65000, 1000.0), 525.0
    )
