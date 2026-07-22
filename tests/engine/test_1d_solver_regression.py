from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def _setup_1d(
    source: str,
    num_stokes: int,
    calculate_derivatives: bool,
    *,
    num_wavelengths: int | None = None,
    wavelength_batch_size: int = 1,
    num_threads: int = 1,
    threading_lib: sk.ThreadingLib | None = None,
    threading_model: sk.ThreadingModel | None = None,
    emission_source: sk.EmissionSource = sk.EmissionSource.NoSource,
    occultation_source: sk.OccultationSource = sk.OccultationSource.NoSource,
    solar_transmission: sk.SingleScatterSolarTransmission = (
        sk.SingleScatterSolarTransmission.Exact
    ),
):
    config = sk.Config()
    config.num_threads = num_threads
    config.num_stokes = num_stokes
    config.num_streams = 8
    config.num_singlescatter_moments = 16
    config.num_sza = 2
    config.output_los_optical_depth = True
    config.wavelength_batch_size = wavelength_batch_size
    config.emission_source = emission_source
    config.occultation_source = occultation_source
    if threading_lib is not None:
        config.threading_lib = threading_lib
    if threading_model is not None:
        config.threading_model = threading_model

    if source == "discrete_ordinates":
        config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
        config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
        default_num_wavelengths = 3
    elif source == "successive_orders":
        config.single_scatter_source = sk.SingleScatterSource.Exact
        config.single_scatter_source_quadrature = True
        config.single_scatter_solar_transmission = solar_transmission
        config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
        config.num_successive_orders_iterations = 3
        config.num_successive_orders_incoming = 26
        config.num_successive_orders_outgoing = 26
        config.init_successive_orders_with_discrete_ordinates = False
        default_num_wavelengths = 2
    else:
        config.single_scatter_source = sk.SingleScatterSource.Exact
        config.single_scatter_source_quadrature = True
        config.single_scatter_solar_transmission = solar_transmission
        config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
        default_num_wavelengths = 1

    if num_wavelengths is None:
        num_wavelengths = default_num_wavelengths

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
    if emission_source != sk.EmissionSource.NoSource:
        atmosphere.storage.emission_source[:] = (
            3.0e-4 * np.exp(-altitude_grid_m / 12_000.0)[:, np.newaxis]
        ) * np.linspace(0.7, 1.4, num_wavelengths)[np.newaxis, :]
        atmosphere.surface.emission[:] = np.linspace(1.0e-4, 7.0e-4, num_wavelengths)

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
                            0.03493718059465559,
                            -0.00124129381179148,
                            -0.01317730962094016,
                        ],
                        [
                            0.01833445832490121,
                            0.00541996960345787,
                            0.00527420441413781,
                        ],
                        [
                            0.10101664567724669,
                            -0.00217722682907605,
                            -0.01381960846066964,
                        ],
                        [
                            0.02369440387264463,
                            -0.00248951297454448,
                            0.00621379400634047,
                        ],
                    ],
                    [
                        [
                            0.06518818242260963,
                            -0.00247087959465485,
                            -0.01980392993454027,
                        ],
                        [
                            0.04389054236089,
                            0.00882382970897506,
                            0.00855995487348336,
                        ],
                        [
                            0.12481456275373466,
                            -0.00318717297589828,
                            -0.01605922155292107,
                        ],
                        [
                            0.04356467734399301,
                            -0.0043753749550128,
                            0.01079868622169172,
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

    # Solver and traced-path results vary at the sub-ppm level across compiler
    # and BLAS combinations; this remains tight enough to catch material drift.
    cross_platform_rtol = 5e-7
    np.testing.assert_allclose(
        result.radiance.values,
        expected_radiance,
        rtol=cross_platform_rtol,
        atol=2e-13,
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
        rtol=cross_platform_rtol,
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


def test_exact_single_scatter_radiance_is_independent_of_derivative_mode():
    """Enabling weighting functions must not change the dispatched radiance."""
    value_engine, value_atmosphere = _setup_1d("single_scatter", 3, False)
    wf_engine, wf_atmosphere = _setup_1d("single_scatter", 3, True)

    value = value_engine.calculate_radiance(value_atmosphere).radiance.values
    with_wf = wf_engine.calculate_radiance(wf_atmosphere).radiance.values
    np.testing.assert_allclose(with_wf, value, rtol=2e-13, atol=2e-14)


@pytest.mark.parametrize(
    "solar_transmission",
    [
        sk.SingleScatterSolarTransmission.Exact,
        sk.SingleScatterSolarTransmission.RayTable,
    ],
)
def test_gauss8_single_scatter_extinction_wf_matches_finite_difference(
    solar_transmission: sk.SingleScatterSolarTransmission,
):
    """Fixed Gauss-8 solar, viewing-OD, and local-source WFs are analytic."""
    engine, atmosphere = _setup_1d(
        "single_scatter", 1, True, solar_transmission=solar_transmission
    )
    result = engine.calculate_radiance(atmosphere)

    for altitude_index in (0, 5, 12, 20, 24):
        original = atmosphere.storage.total_extinction[altitude_index, 0]
        step = max(1e-10, original * 2e-5)
        atmosphere.storage.total_extinction[altitude_index, 0] = original + step
        above = engine.calculate_radiance(atmosphere).radiance.values.copy()
        atmosphere.storage.total_extinction[altitude_index, 0] = original - step
        below = engine.calculate_radiance(atmosphere).radiance.values.copy()
        atmosphere.storage.total_extinction[altitude_index, 0] = original

        numeric = (above - below) / (2 * step)
        analytic = result.wf_extinction.values[altitude_index, 0]
        np.testing.assert_allclose(analytic, numeric[0], rtol=2e-7, atol=2e-9)


@pytest.mark.parametrize(
    ("num_stokes", "calculate_derivatives"),
    [(1, False), (1, True), (3, False), (3, True)],
)
def test_exact_single_scatter_wavelength_batch_parity(
    num_stokes: int, calculate_derivatives: bool
):
    """Batch execution, including a partial tail, matches the scalar path."""
    scalar_engine, scalar_atmosphere = _setup_1d(
        "single_scatter",
        num_stokes,
        calculate_derivatives,
        num_wavelengths=7,
    )
    batch_engine, batch_atmosphere = _setup_1d(
        "single_scatter",
        num_stokes,
        calculate_derivatives,
        num_wavelengths=7,
        wavelength_batch_size=4,
    )

    scalar = scalar_engine.calculate_radiance(scalar_atmosphere)
    batched = batch_engine.calculate_radiance(batch_atmosphere)

    assert scalar.data_vars.keys() == batched.data_vars.keys()
    # GEMV and GEMM solar-transmission paths differ by a few final bits. The
    # log-transmission weights expose that difference most strongly in WFs.
    parity_rtol = 1e-11 if calculate_derivatives else 2e-13
    for variable in scalar.data_vars:
        np.testing.assert_allclose(
            batched[variable].values,
            scalar[variable].values,
            rtol=parity_rtol,
            atol=2e-14,
        )


@pytest.mark.parametrize(
    "emission_source",
    [sk.EmissionSource.Standard, sk.EmissionSource.VolumeEmissionRate],
)
def test_mixed_sources_wavelength_batch_parity(emission_source):
    options = {
        "num_wavelengths": 7,
        "emission_source": emission_source,
        "occultation_source": sk.OccultationSource.Standard,
    }
    scalar_engine, scalar_atmosphere = _setup_1d("single_scatter", 3, True, **options)
    batch_engine, batch_atmosphere = _setup_1d(
        "single_scatter", 3, True, wavelength_batch_size=4, **options
    )

    scalar = scalar_engine.calculate_radiance(scalar_atmosphere)
    batched = batch_engine.calculate_radiance(batch_atmosphere)
    for variable in scalar.data_vars:
        np.testing.assert_allclose(
            batched[variable].values,
            scalar[variable].values,
            rtol=1e-11,
            atol=2e-14,
        )


@pytest.mark.parametrize("wavelength_batch_size", [1, 4])
def test_surface_albedo_wf_with_atmospheric_emission(wavelength_batch_size):
    """Atmospheric-emission derivatives must not displace surface derivatives."""
    engine, atmosphere = _setup_1d(
        "single_scatter",
        1,
        True,
        num_wavelengths=5,
        wavelength_batch_size=wavelength_batch_size,
        emission_source=sk.EmissionSource.Standard,
    )
    result = engine.calculate_radiance(atmosphere)

    step = 1e-5
    atmosphere.surface.albedo[:] += step
    radiance_above = engine.calculate_radiance(atmosphere).radiance.values.copy()
    atmosphere.surface.albedo[:] -= 2 * step
    radiance_below = engine.calculate_radiance(atmosphere).radiance.values.copy()
    atmosphere.surface.albedo[:] += step

    numeric_wf = (radiance_above - radiance_below) / (2 * step)
    np.testing.assert_allclose(
        result.wf_albedo.values,
        numeric_wf,
        rtol=1e-7,
        atol=1e-10,
    )


def test_rayon_wavelength_batch_parity():
    scalar_engine, scalar_atmosphere = _setup_1d(
        "single_scatter", 1, True, num_wavelengths=9
    )
    batch_engine, batch_atmosphere = _setup_1d(
        "single_scatter",
        1,
        True,
        num_wavelengths=9,
        wavelength_batch_size=4,
        num_threads=2,
        threading_lib=sk.ThreadingLib.Rayon,
    )

    scalar = scalar_engine.calculate_radiance(scalar_atmosphere)
    batched = batch_engine.calculate_radiance(batch_atmosphere)
    for variable in scalar.data_vars:
        np.testing.assert_allclose(
            batched[variable].values,
            scalar[variable].values,
            rtol=1e-11,
            atol=2e-13,
        )


def test_unsupported_sources_fall_back_to_scalar_wavelengths():
    scalar_engine, scalar_atmosphere = _setup_1d(
        "discrete_ordinates", 1, False, num_wavelengths=3
    )
    configured_engine, configured_atmosphere = _setup_1d(
        "discrete_ordinates",
        1,
        False,
        num_wavelengths=3,
        wavelength_batch_size=8,
    )

    scalar = scalar_engine.calculate_radiance(scalar_atmosphere)
    configured = configured_engine.calculate_radiance(configured_atmosphere)
    for variable in scalar.data_vars:
        np.testing.assert_array_equal(
            configured[variable].values, scalar[variable].values
        )


def test_mixed_batch_and_nonbatch_sources_fall_back_to_scalar_wavelengths():
    """A non-batched source caps exact single scatter to the scalar kernel."""
    scalar_engine, scalar_atmosphere = _setup_1d(
        "successive_orders", 1, True, num_wavelengths=3
    )
    configured_engine, configured_atmosphere = _setup_1d(
        "successive_orders",
        1,
        True,
        num_wavelengths=3,
        wavelength_batch_size=8,
    )

    scalar = scalar_engine.calculate_radiance(scalar_atmosphere)
    configured = configured_engine.calculate_radiance(configured_atmosphere)
    for variable in scalar.data_vars:
        np.testing.assert_array_equal(
            configured[variable].values, scalar[variable].values
        )


def test_rayon_source_threading_uses_cpp_source_threads():
    """Rayon only partitions wavelengths; source threading stays in C++."""
    scalar_engine, scalar_atmosphere = _setup_1d(
        "single_scatter", 1, True, num_wavelengths=5
    )
    threaded_engine, threaded_atmosphere = _setup_1d(
        "single_scatter",
        1,
        True,
        num_wavelengths=5,
        wavelength_batch_size=4,
        num_threads=2,
        threading_lib=sk.ThreadingLib.Rayon,
        threading_model=sk.ThreadingModel.Source,
    )

    scalar = scalar_engine.calculate_radiance(scalar_atmosphere)
    threaded = threaded_engine.calculate_radiance(threaded_atmosphere)
    for variable in scalar.data_vars:
        np.testing.assert_allclose(
            threaded[variable].values,
            scalar[variable].values,
            rtol=1e-11,
            atol=2e-13,
        )
