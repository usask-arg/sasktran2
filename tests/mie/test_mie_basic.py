import matplotlib.pyplot as plt
import numpy as np
import sasktran2 as sk
import scipy.interpolate as interp
import xarray as xr

# from sasktran.legendre.wigner import WignerD
from scipy import interpolate
from scipy.special import roots_legendre


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


def test_mie_integration():
    # mie = sk.mie.LinearizedMie()
    # refractive_index = 1 - 0.1j
    # size_param = np.array(
    #   [1]
    # )  # not checking values for 0.2, just there to show that we can calculate with multiple size param and angles
    # angles = np.linspace(0, 180, 3)
    # cos_angles = np.cos(angles * np.pi / 180)
    # output = mie.calculate(size_param, refractive_index, cos_angles, True)
    from scipy.stats import rv_continuous

    class gaussian_gen(rv_continuous):
        "Gaussian distribution"

        def _pdf(self, x):
            mu = 500  # 500
            sig = 100  # 100
            return np.exp(-(((x - mu) / sig) ** 2) / 2.0) / np.sqrt(2.0 * np.pi) / sig

    prob_dist = gaussian_gen(name="gaussian")

    def refrac_index_fn():
        return 1.0 - 0.1j

    wavelengths = np.array([400, 500])
    my_mie = sk.mie.LinearizedMie()

    out = sk.mie.integrate_mie(
        my_mie,
        prob_dist,
        refrac_index_fn,
        wavelengths,
        num_angles=3,
        num_quad=20,
        compute_coeffs=True,
        num_coeffs=2,
    )
    np.testing.assert_allclose(out["p11"][0][0].data, 8.11553105e01, atol=1e-06, rtol=0)
    np.testing.assert_allclose(
        out["p11"][0][1].data, 1.16641980e-02, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        out["p11"][0][2].data, 5.64571171e-03, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(out["p12"][0][0].data, 0, atol=1e-06, rtol=0)
    np.testing.assert_allclose(out["p12"][0][1].data, 0.01117573, atol=1e-06, rtol=0)
    np.testing.assert_allclose(out["p12"][0][2].data, 0, atol=1e-06, rtol=0)
    np.testing.assert_allclose(out["p33"][0][0].data, 8.11553105e01, atol=1e-06, rtol=0)
    np.testing.assert_allclose(
        out["p33"][0][1].data, -6.33399061e-04, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(
        out["p33"][0][2].data, -5.64571171e-03, atol=1e-06, rtol=0
    )
    np.testing.assert_allclose(out["p34"][0][0].data, 0, atol=1e-06, rtol=0)
    np.testing.assert_allclose(out["p34"][0][1].data, -0.00099857, atol=1e-06, rtol=0)
    np.testing.assert_allclose(out["p34"][0][2].data, 0, atol=1e-06, rtol=0)

    np.testing.assert_allclose(
        out["xs_total"][0].data, 1077477.99015462, atol=1e-03, rtol=0
    )
    np.testing.assert_allclose(
        out["xs_scattering"][0].data, 348484.73075379, atol=1e-03, rtol=0
    )
    np.testing.assert_allclose(
        out["xs_absorption"][0].data, 728993.25940083, atol=1e-03, rtol=0
    )


def test_external_calc():
    data = xr.open_dataset(
        "\\\\utls\\SasktranFiles\\BaumIceCrystals\\GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc"
    )
    from sasktran2.legendre import compute_greek_coefficients

    num_coeffs = 10000
    wavelengths = data["wavelengths"].data

    lm_a1_master = np.zeros((len(wavelengths), num_coeffs, len(data["nDeff"].data)))
    lm_a2_master = np.zeros((len(wavelengths), num_coeffs, len(data["nDeff"].data)))
    lm_a3_master = np.zeros((len(wavelengths), num_coeffs, len(data["nDeff"].data)))
    lm_a4_master = np.zeros((len(wavelengths), num_coeffs, len(data["nDeff"].data)))
    lm_b1_master = np.zeros((len(wavelengths), num_coeffs, len(data["nDeff"].data)))
    lm_b2_master = np.zeros((len(wavelengths), num_coeffs, len(data["nDeff"].data)))

    # put total cross section in units of m ^2
    xs_total = data["extinction_efficiency"].data * data["total_area"].data / 1e6 / 1e6
    xs_scat = xs_total * data["single_scattering_albedo"].data

    angles = data["phase_angles"].data
    # angles_interp = np.linspace(0,180,10000)
    diff = angles[1:] - angles[:-1]
    sub_divisions = 16
    angles_interp = angles
    for i in range(sub_divisions - 1):
        angles_interp = np.append(
            angles_interp, angles[:-1] + (i + 1) * diff / sub_divisions
        )

    angles_interp = np.sort(angles_interp)
    for diameter_idx in data["nDeff"].data:
        p11 = data["p11_phase_function"].data[:, diameter_idx, :].transpose()
        p11_s = interp.CubicSpline(angles, p11, axis=1)
        p11_interp = p11_s(angles_interp)
        # rest of 'phase functions' are normalized by p11, multiply by p11 to return to full phase function
        p12 = data["p21_phase_function"].data[:, diameter_idx, :].transpose() * p11
        p12_s = interp.CubicSpline(angles, p12, axis=1)
        p12_interp = p12_s(angles_interp)

        p22 = data["p22_phase_function"].data[:, diameter_idx, :].transpose() * p11
        p22_s = interp.CubicSpline(angles, p22, axis=1)
        p22_interp = p22_s(angles_interp)

        p33 = data["p33_phase_function"].data[:, diameter_idx, :].transpose() * p11
        p33_s = interp.CubicSpline(angles, p33, axis=1)
        p33_interp = p33_s(angles_interp)

        p34 = -1 * data["p43_phase_function"].data[:, diameter_idx, :].transpose() * p11
        p34_s = interp.CubicSpline(angles, p34, axis=1)
        p34_interp = p34_s(angles_interp)

        p44 = data["p44_phase_function"].data[:, diameter_idx, :].transpose() * p11
        p44_s = interp.CubicSpline(angles, p44, axis=1)
        p44_interp = p44_s(angles_interp)

        lm_a1, lm_a2, lm_a3, lm_a4, lm_b1, lm_b2 = compute_greek_coefficients(
            p11=p11_interp,
            p12=p12_interp,
            p22=p22_interp,
            p33=p33_interp,
            p34=p34_interp,
            p44=p44_interp,
            angle_grid=angles_interp,
            num_coeff=num_coeffs,
        )
        lm_a1_master[:, :, diameter_idx] = lm_a1
        lm_a2_master[:, :, diameter_idx] = lm_a2
        lm_a3_master[:, :, diameter_idx] = lm_a3
        lm_a4_master[:, :, diameter_idx] = lm_a4
        lm_b1_master[:, :, diameter_idx] = lm_b1
        lm_b2_master[:, :, diameter_idx] = lm_b2

    # wavelengths in nm
    coeffs_output = xr.Dataset(
        {
            "lm_a1": (["wavelength", "num_terms", "nDeff"], lm_a1_master),
            "lm_a2": (["wavelength", "num_terms", "nDeff"], lm_a2_master),
            "lm_a3": (["wavelength", "num_terms", "nDeff"], lm_a3_master),
            "lm_a4": (["wavelength", "num_terms", "nDeff"], lm_a4_master),
            "lm_b1": (["wavelength", "num_terms", "nDeff"], lm_b1_master),
            "lm_b2": (["wavelength", "num_terms", "nDeff"], lm_b2_master),
            "xs_total": (["wavelength", "nDeff"], xs_total),
            "xs_scat": (["wavelength", "nDeff"], xs_scat),
        },
        coords={
            "wavelength": wavelengths * 1e3,
            "num_terms": np.arange(num_coeffs),
            "nDeff": data["nDeff"].data,
        },
    )

    coeffs_output.to_netcdf(
        "\\\\utls\\SasktranFiles\\BaumIceCrystals\\GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix_processed.nc"
    )


def test_dan_calc():
    ds = xr.open_dataset(
        "\\\\utls\\SasktranFiles\\BaumIceCrystals\\GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc"
    )

    phase = (
        ds.isel(nDeff=0, nWaveLen=50)["p11_phase_function"]
        .to_numpy()
        .astype(np.float64)[::-1]
    )
    cos_angles = np.cos(np.deg2rad(ds["phase_angles"].to_numpy().astype(np.float64)))[
        ::-1
    ]

    c = 0.995
    n = 4000

    nodes, weights = roots_legendre(n)
    # splitting into peak and non-peak

    nodes_left = (c - (-1)) / 2 * nodes + (c + (-1)) / 2
    weights_left = (c - (-1)) / 2 * weights

    nodes_right = (1 - c) / 2 * nodes + (1 + c) / 2
    weights_right = (1 - c) / 2 * weights

    all_weights = np.concatenate([weights_left, weights_right])
    cos_angle_grid = np.concatenate([nodes_left, nodes_right])

    # wig = WignerD(0, 0)
    wig = sk.util.WignerD(0, 0)

    l_vals = np.zeros((n, 2 * n))
    for i in range(n):
        l_vals[i] = wig.d(np.arccos(cos_angle_grid), i)
        l_vals[i] *= all_weights / (2.0 / (2.0 * i + 1))

    phase_interp = interpolate.PchipInterpolator(cos_angles, phase)(cos_angle_grid)

    a1 = l_vals @ phase_interp

    plt.plot(a1)
    plt.yscale("log")
