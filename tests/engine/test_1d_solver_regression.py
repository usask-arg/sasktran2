from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def _setup_1d(source: str, num_stokes: int, calculate_derivatives: bool):
    config = sk.Config()
    config.num_threads = 1
    config.num_stokes = num_stokes
    config.num_streams = 8
    config.num_singlescatter_moments = 16
    config.num_sza = 2
    config.output_los_optical_depth = True

    if source == "discrete_ordinates":
        config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
        config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
        num_wavelengths = 3
    elif source == "successive_orders":
        config.single_scatter_source = sk.SingleScatterSource.Exact
        config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
        config.num_successive_orders_iterations = 3
        config.num_successive_orders_incoming = 26
        config.num_successive_orders_outgoing = 26
        config.init_successive_orders_with_discrete_ordinates = False
        num_wavelengths = 2
    else:
        config.single_scatter_source = sk.SingleScatterSource.Exact
        config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
        num_wavelengths = 1

    altitude_grid_m = np.linspace(0.0, 60_000.0, 25)
    cos_sza = 0.42
    geometry = sk.Geometry1D(
        cos_sza=cos_sza,
        solar_azimuth=0.35,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=altitude_grid_m,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )

    viewing_geometry = sk.ViewingGeometry()
    viewing_geometry.add_ray(sk.GroundViewingSolar(cos_sza, -0.7, 0.32, 200_000.0))
    viewing_geometry.add_ray(sk.GroundViewingSolar(cos_sza, 0.4, 0.78, 200_000.0))
    viewing_geometry.add_ray(
        sk.TangentAltitudeSolar(12_345.0, -0.35, 200_000.0, cos_sza)
    )
    viewing_geometry.add_ray(
        sk.TangentAltitudeSolar(27_123.0, 0.65, 200_000.0, cos_sza)
    )

    atmosphere = sk.Atmosphere(
        geometry,
        config,
        calculate_derivatives=calculate_derivatives,
        numwavel=num_wavelengths,
    )
    altitude_factor = np.exp(-altitude_grid_m / 7_500.0)[:, np.newaxis]
    spectral_factor = np.linspace(0.72, 1.35, num_wavelengths)[np.newaxis, :]
    atmosphere.storage.total_extinction[:] = (
        2.4e-5 * altitude_factor + 1.0e-9
    ) * spectral_factor
    atmosphere.storage.ssa[:] = (
        0.91
        + 0.025 * np.exp(-altitude_grid_m / 18_000.0)[:, np.newaxis]
        - 0.01 * np.linspace(0.0, 1.0, num_wavelengths)[np.newaxis, :]
    )
    atmosphere.leg_coeff.a1[0, :, :] = 1.0
    atmosphere.leg_coeff.a1[1, :, :] = 0.08
    atmosphere.leg_coeff.a1[2, :, :] = 0.5
    if num_stokes == 3:
        atmosphere.leg_coeff.a2[2, :, :] = 3.0
        atmosphere.leg_coeff.b1[2, :, :] = -np.sqrt(6.0) / 2.0
    atmosphere.surface.albedo[:] = np.linspace(0.08, 0.31, num_wavelengths)

    return sk.Engine(config, geometry, viewing_geometry), atmosphere


@pytest.mark.parametrize(
    ("source", "num_stokes", "calculate_derivatives", "expected_radiance"),
    [
        (
            "discrete_ordinates",
            1,
            True,
            np.array(
                [
                    [
                        0.007325034080894165,
                        0.003911222331552876,
                        0.015579706096361565,
                        0.0038272940915240537,
                    ],
                    [
                        0.015246498141181428,
                        0.009498372748059328,
                        0.024575183549304623,
                        0.007446741685194283,
                    ],
                    [
                        0.023780265249057513,
                        0.016356191357612092,
                        0.031606480325730144,
                        0.01139900426172521,
                    ],
                ]
            )[:, :, np.newaxis],
        ),
        (
            "successive_orders",
            3,
            False,
            np.array(
                [
                    [
                        [
                            0.03485848278767426,
                            -0.0012438171986974716,
                            -0.01313439668613118,
                        ],
                        [
                            0.018305818732869187,
                            0.005400715559354511,
                            0.005256718675979713,
                        ],
                        [
                            0.1009590532632837,
                            -0.002179747695663314,
                            -0.013811326250175521,
                        ],
                        [
                            0.023692987924552825,
                            -0.0024894020696831325,
                            0.006213385479745916,
                        ],
                    ],
                    [
                        [
                            0.06499456949466546,
                            -0.0024770876498211134,
                            -0.019698355214306774,
                        ],
                        [
                            0.04381448739788437,
                            0.008772698896414867,
                            0.008513519953347171,
                        ],
                        [
                            0.12467626601396885,
                            -0.0031932263361659366,
                            -0.016039333469915443,
                        ],
                        [
                            0.04355957469824446,
                            -0.004374975287664437,
                            0.010797214016246497,
                        ],
                    ],
                ]
            ),
        ),
    ],
)
def test_1d_solver_radiance_regression(
    source, num_stokes, calculate_derivatives, expected_radiance
):
    engine, atmosphere = _setup_1d(source, num_stokes, calculate_derivatives)
    result = engine.calculate_radiance(atmosphere)

    # Iterative solver results vary at the sub-ppm level across compiler and
    # BLAS combinations; this remains tight enough to catch material drift.
    np.testing.assert_allclose(
        result.radiance.values, expected_radiance, rtol=5e-7, atol=2e-13
    )

    expected_optical_depth = np.array(
        [
            [
                0.4046672641890289,
                0.16756248699288395,
                1.8429222725097874,
                0.2573659860040063,
            ],
            [
                0.581709192271729,
                0.2408710750522707,
                2.64920076673282,
                0.3699636048807591,
            ],
            [
                0.7587511203544293,
                0.31417966311165746,
                3.455479260955851,
                0.4825612237575118,
            ],
        ]
    )
    if source == "successive_orders":
        expected_optical_depth = expected_optical_depth[[0, 2]]
    np.testing.assert_allclose(
        result.los_optical_depth.values,
        expected_optical_depth,
        rtol=2e-13,
        atol=1e-13,
    )


def test_1d_polarized_b1_phase_coefficient_wf():
    """The polarized b1 phase derivative must not be truncated to zero."""
    engine, atmosphere = _setup_1d("single_scatter", 3, True)
    result = engine.calculate_radiance(atmosphere)

    altitude_index = 5
    coefficient_order = 2
    step = 1e-6
    original = atmosphere.leg_coeff.b1[coefficient_order, altitude_index, 0]

    atmosphere.leg_coeff.b1[coefficient_order, altitude_index, 0] = original + step
    radiance_plus = engine.calculate_radiance(atmosphere).radiance.values.copy()
    atmosphere.leg_coeff.b1[coefficient_order, altitude_index, 0] = original - step
    radiance_minus = engine.calculate_radiance(atmosphere).radiance.values.copy()
    atmosphere.leg_coeff.b1[coefficient_order, altitude_index, 0] = original

    numeric_wf = (radiance_plus - radiance_minus) / (2 * step)
    analytic_wf = result.wf_leg_coeff_11.values[altitude_index, 0]
    np.testing.assert_allclose(analytic_wf, numeric_wf[0], rtol=1e-7, atol=1e-10)
