from __future__ import annotations

import numpy as np
from sasktran2.optical.rayleigh import rayleigh_cross_section_bates


def test_rayleigh_bates():
    _ = rayleigh_cross_section_bates(np.array([0.350]))
