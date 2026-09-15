from __future__ import annotations

import numpy as np
import pytest
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


@pytest.mark.parametrize(("m", "n"), [(0, 0), (2, 2), (2, -2), (0, 2), (1, 0), (-1, 2)])
def test_wigner_all_orders_and_angles(m, n):
    wig = sk.util.WignerD(m, n)
    # Include endpoints, nearly forward/backward scattering and a strided grid.
    theta = np.array([0.0, 1e-7, 0.4, 1.2, 2.0, np.pi - 1e-7, np.pi])
    storage = np.repeat(theta, 2)
    result = wig.d_all(storage[::2], 513)
    assert result.shape == (513, len(theta))
    for order in [0, 1, 2, 3, 16, 64, 256, 512]:
        np.testing.assert_allclose(
            result[order], wig.d(theta, order), atol=3e-11, rtol=3e-11
        )
        np.testing.assert_array_equal(wig.d(storage[::2], order), wig.d(theta, order))
        if m == n == 0:
            np.testing.assert_allclose(
                result[order], eval_legendre(order, np.cos(theta)), atol=3e-11
            )


def test_wigner_order_selection_and_empty_grids():
    wig = sk.util.WignerD(2, -2)
    orders = np.array([8, 2, 0, 5, 8, -1], dtype=np.int32)
    expected = np.array([wig.d(np.array([0.7]), int(order))[0] for order in orders])
    np.testing.assert_allclose(wig.d_vec(0.7, orders), expected, atol=1e-14)
    assert wig.d_vec(0.7, np.array([], dtype=np.int32)).shape == (0,)
    assert wig.d_all(np.array([]), 5).shape == (5, 0)
    assert wig.d_all(np.array([0.7]), 0).shape == (0, 1)
