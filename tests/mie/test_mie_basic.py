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


def test_qext_q_Sca_s1_s2():
    mie = sk.mie.LinearizedMie()
    refractive_index = 0.75
    size_param = np.array(
        [0.101, 0.2]
    )  # not checking values for 0.2, just there to show that we can calculate with multiple size param and angles
    angles = np.arange(0, 181, 30)
    cos_angles = np.cos(angles * np.pi / 180)
    output = mie.calculate(size_param, refractive_index, cos_angles, True)

    np.testing.assert_allclose(output.values.Qext[0], 0.000008, atol=1e-06, rtol=0)
    np.testing.assert_allclose(output.values.Qsca[0], 0.000008, atol=1e-06, rtol=0)

    np.testing.assert_allclose(
        output.values.S1[0][0].real, 2.048754e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][0].imag, -1.756419e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][0].real, 2.048754e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][0].imag, -1.756419e-04, atol=1e-06, rtol=0
    )

    np.testing.assert_allclose(
        output.values.S1[0][1].real, 2.048754e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][1].imag, -1.755965e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][1].real, 1.774273e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][1].imag, -1.520629e-04, atol=1e-06, rtol=0
    )

    np.testing.assert_allclose(
        output.values.S1[0][2].real, 2.048753e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][2].imag, -1.754726e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][2].real, 1.024377e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][2].imag, -8.771198e-05, atol=1e-06, rtol=0
    )

    np.testing.assert_allclose(
        output.values.S1[0][3].real, 2.048751e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][3].imag, -1.753033e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][3].real, 1.845057e-15, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][3].imag, -3.247270e-08, atol=1e-06, rtol=0
    )

    np.testing.assert_allclose(
        output.values.S1[0][4].real, 2.048750e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][4].imag, -1.751341e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][4].real, -1.024375e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][4].imag, 8.759147e-05, atol=1e-06, rtol=0
    )

    np.testing.assert_allclose(
        output.values.S1[0][5].real, 2.048749e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][5].imag, -1.750103e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][5].real, -1.774269e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][5].imag, 1.515715e-04, atol=1e-06, rtol=0
    )

    np.testing.assert_allclose(
        output.values.S1[0][4].real, 2.048749e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S1[0][5].imag, -1.749650e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][6].real, -2.048749e-08, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        output.values.S2[0][6].imag, 1.749650e-04, atol=1e-06, rtol=0
    )
