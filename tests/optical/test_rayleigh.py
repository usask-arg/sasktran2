import numpy as np
from sasktran2.optical.rayleigh import rayleigh_cross_section_bates


def test_rayleigh_bates():
    wavelengths_um = np.arange(0.1, 5, 0.0001)

    xs = rayleigh_cross_section_bates(np.array([0.350]))
