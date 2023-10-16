import numpy as np
import sasktran2 as sk


def test_mie_construction():
    mie = sk.mie.LinearizedMie()

    size_param = np.array([3.0, 4.0, 5.0])
    cos_angles = np.linspace(-1, 1, 100)
    refractive_index = 1.5 + 0.0j

    output = mie.calculate(size_param, refractive_index, cos_angles, True)

    np.testing.assert_allclose(output.size_parameter, size_param)
    np.testing.assert_allclose(output.cos_angles, cos_angles)
