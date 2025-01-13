from __future__ import annotations

import numpy as np
import sasktran2 as sk
from scipy.special import eval_legendre


def test_wigner_construction():
    _ = sk.util.WignerD(0, 0)


def test_wigner_vectorization():
    wd = sk.util.WignerD(0, 0)

    theta = np.arange(0, np.pi, 0.01)

    _ = wd.d(theta, 10)


def test_wigner_against_legendre():
    wig = sk.util.WignerD(0, 0)

    theta = np.arange(0, np.pi, 0.01)

    for legidx in range(20):
        true = eval_legendre(legidx, np.cos(theta))

        ours = wig.d(theta, legidx)

        np.testing.assert_array_almost_equal(true, ours)
