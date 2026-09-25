"""Independent physics checks for resonance optical properties, without downloads."""

from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2.optical.resonance import (
    LineResonance,
    resonance_phase_function,
    resonance_polarizability,
)
from scipy.constants import atomic_mass, c, e, epsilon_0, h, k, m_e
from scipy.integrate import trapezoid


def _database(**overrides):
    """A synthetic stable lower state and one electric-dipole transition."""
    fields = {
        "wavelength_nm": [500.0],
        "oscillator_strength": [0.4],
        "einstein_a_s": [1.0e6],
        "upper_total_a_s": [1.0e6],
        "lower_energy_cminv": [0.0],
        "lower_j": [0.0],
        "upper_j": [1.0],
    }
    fields.update(overrides)
    count = max(np.size(value) for value in fields.values())
    data = {
        name: ("line", np.broadcast_to(value, (count,)).astype(float).copy())
        for name, value in fields.items()
    }
    data["energy_cminv"] = ("state", [0.0])
    data["statistical_weight"] = ("state", [1.0])
    return xr.Dataset(
        data,
        attrs={
            "mass_amu": 23.0,
            "species": "synthetic",
            "wavelength_medium": "vacuum",
        },
    )


def _frequency_grid(wavelength_nm=500.0, temperature_k=200.0):
    center_hz = c / (wavelength_nm * 1e-9)
    doppler_std = center_hz / c * np.sqrt(k * temperature_k / (23.0 * atomic_mass))
    frequencies = center_hz + np.linspace(-100, 100, 20001) * doppler_std
    return frequencies, c / frequencies * 1e9


def _evaluate(database, wavelengths_nm, temperatures=(200.0,), **kwargs):
    temperatures = np.atleast_1d(temperatures)
    return LineResonance(db=database, line_wing_cutoff_nm=None).cross_sections(
        np.asarray(wavelengths_nm),
        np.arange(len(temperatures), dtype=float) * 1000,
        temperature_k=temperatures,
        **kwargs,
    )


@pytest.mark.parametrize("temperature", [150.0, 600.0])
def test_frequency_integral_obeys_oscillator_strength_sum_rule(temperature):
    frequencies, wavelengths = _frequency_grid(temperature_k=temperature)
    result = _evaluate(_database(), wavelengths, [temperature])
    # The sum rule is in Hz, not angular frequency or wavelength.
    expected_area = e**2 / (4 * epsilon_0 * m_e * c) * 0.4
    actual_area = trapezoid(result.extinction[0], frequencies)
    np.testing.assert_allclose(actual_area, expected_area, rtol=3e-6)
    np.testing.assert_allclose(result.ssa, 1.0)


def test_doppler_width_increases_as_sqrt_temperature():
    database = _database(einstein_a_s=[1.0], upper_total_a_s=[1.0])
    center_hz = c / 500e-9
    cold_sigma = center_hz / c * np.sqrt(k * 150 / (23 * atomic_mass))
    wavelength = c / (center_hz + np.array([0, 1, 2]) * cold_sigma) * 1e9
    result = _evaluate(database, wavelength, [150, 600])
    np.testing.assert_allclose(
        result.extinction[0, 0] / result.extinction[1, 0], 2.0, rtol=1e-8
    )
    np.testing.assert_allclose(
        result.extinction[0, 1] / result.extinction[0, 0], np.exp(-0.5), rtol=1e-8
    )
    np.testing.assert_allclose(
        result.extinction[1, 2] / result.extinction[1, 0], np.exp(-0.5), rtol=1e-8
    )


@pytest.mark.parametrize(
    ("lower_j", "upper_j", "expected"),
    [(0, 1, 1), (0.5, 0.5, 0), (0.5, 1.5, 0.5), (1, 1, 0.25), (1, 0, 0)],
)
def test_known_resonance_polarizabilities(lower_j, upper_j, expected):
    np.testing.assert_allclose(resonance_polarizability(lower_j, upper_j), expected)


def test_molecular_large_j_limits():
    lower = np.array([1e5, 1e5, 1e5])
    upper = lower + np.array([-1, 0, 1])
    np.testing.assert_allclose(
        resonance_polarizability(lower, upper), [0.1, 0.4, 0.1], rtol=4e-5
    )


def test_phase_normalization_and_known_angular_values():
    mu, weights = np.polynomial.legendre.leggauss(8)
    for polarizability in [0.0, 0.25, 0.5, 1.0]:
        phase = resonance_phase_function(mu, polarizability)
        np.testing.assert_allclose(weights @ phase, 2.0, atol=1e-14)
        np.testing.assert_allclose(
            resonance_phase_function(np.array([-1.0, 0.0, 1.0]), polarizability),
            [1 + polarizability / 2, 1 - polarizability / 4, 1 + polarizability / 2],
        )


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_branching_controls_scattering_and_overlap_polarization(num_stokes):
    database = _database(
        wavelength_nm=[500.0, 500.0],
        lower_j=[0.5, 0.5],
        upper_j=[0.5, 1.5],
        einstein_a_s=[2e7, 8e7],
        upper_total_a_s=[1e8, 1e8],
    )
    database["statistical_weight"][:] = 2.0
    result = _evaluate(
        database, [499.999, 500, 500.001], num_stokes=num_stokes, num_moments=5
    )
    np.testing.assert_allclose(result.ssa, 0.5)
    # Equal absorption strengths, but D2 has four times the elastic return rate.
    np.testing.assert_allclose(result.polarizability, 0.4)
    legendre = sk.polarization.LegendreStorageView(result.leg_coeff, num_stokes)
    np.testing.assert_allclose(legendre.a1[0], 1.0)
    np.testing.assert_allclose(legendre.a1[2], 0.2)
    np.testing.assert_allclose(legendre.a1[[1, 3, 4]], 0.0)
    if num_stokes == 3:
        np.testing.assert_allclose(legendre.a2[2], 1.2)
        np.testing.assert_allclose(legendre.b1[2], np.sqrt(6) * 0.2)
        np.testing.assert_allclose(legendre.a3, 0.0)


def test_lte_population_includes_nonabsorbing_states():
    database = _database()
    database = database.drop_dims("state")
    database["energy_cminv"] = ("state", [0, 100])
    database["statistical_weight"] = ("state", [1, 3])
    temperatures = np.array([150.0, 600.0])
    populated = _evaluate(database, [500.0], temperatures).extinction[:, 0]
    ground_only = _evaluate(_database(), [500.0], temperatures).extinction[:, 0]
    expected = 1 / (1 + 3 * np.exp(-100 * 100 * h * c / (k * temperatures)))
    np.testing.assert_allclose(populated / ground_only, expected)


def _atmosphere(num_stokes, wavelengths=(500.0,), *, calculate_derivatives=False):
    config = sk.Config()
    config.num_stokes = num_stokes
    config.num_streams = 4
    config.num_singlescatter_moments = 4
    config.delta_m_scaling = False
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    altitudes = np.arange(0, 60001, 10000, dtype=float)
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6372000.0,
        altitudes,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.asarray(wavelengths),
        calculate_derivatives=calculate_derivatives,
    )
    atmosphere.temperature_k = np.full_like(altitudes, 200.0)
    atmosphere.pressure_pa = np.full_like(altitudes, 1e-10)
    return atmosphere, config, geometry, altitudes


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_number_density_adapter_and_engine_match_classical_rayleigh(num_stokes):
    atmosphere, config, geometry, altitudes = _atmosphere(num_stokes)
    optical = LineResonance(db=_database())
    quantities = optical.atmosphere_quantities(atmosphere)
    # Give each atmosphere the same thin, homogeneous scattering coefficient.
    desired_extinction = 1e-7
    density = desired_extinction / quantities.extinction[:, 0]
    atmosphere["metal"] = sk.constituent.NumberDensityScatterer(
        optical, altitudes, density
    )
    rayleigh_atmosphere, _, _, _ = _atmosphere(num_stokes)
    air_density = rayleigh_atmosphere.pressure_pa / (
        k * rayleigh_atmosphere.temperature_k
    )
    rayleigh_atmosphere["rayleigh"] = sk.constituent.Rayleigh(
        method="manual",
        wavelengths_nm=np.array([499.0, 501.0]),
        xs=np.full(2, desired_extinction / air_density[0]),
        king_factor=np.ones(2),
    )
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.TangentAltitudeSolar(20000, 0, 200000, 0.6))
    engine = sk.Engine(config, geometry, viewing)
    actual = engine.calculate_radiance(atmosphere).radiance.to_numpy()
    expected = engine.calculate_radiance(rayleigh_atmosphere).radiance.to_numpy()
    assert np.all(np.isfinite(actual))
    assert np.all(actual[..., 0] > 0)
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-15)


@pytest.mark.parametrize("temperature", [0.0, -100.0, np.nan, np.inf])
def test_invalid_temperature_rejected(temperature):
    with pytest.raises(ValueError, match="temperature"):
        _evaluate(_database(), [500.0], [temperature])


@pytest.mark.parametrize("wavelength", [0.0, -500.0, np.nan, np.inf])
def test_invalid_wavelength_rejected(wavelength):
    with pytest.raises(ValueError, match="wavelength"):
        _evaluate(_database(), [wavelength])


@pytest.mark.parametrize(
    ("lower_j", "upper_j"), [(-1, 0), (0, 0), (0.2, 1.2), (0, 2), (np.nan, 1)]
)
def test_invalid_dipole_transition_rejected(lower_j, upper_j):
    with pytest.raises(ValueError, match=r"angular|dipole|quantum|transition"):
        resonance_polarizability(lower_j, upper_j)


@pytest.mark.parametrize("partition_table", [False, True])
def test_temperature_derivatives_include_profiles_populations_and_phase(
    partition_table,
):
    database = _database(
        wavelength_nm=[500.0, 500.0006],
        lower_energy_cminv=[0.0, 100.0],
        lower_j=[0.0, 1.0],
        upper_j=[1.0, 1.0],
        oscillator_strength=[0.3, 0.7],
        einstein_a_s=[1e7, 2e7],
        upper_total_a_s=[2e7, 8e7],
    ).drop_dims("state")
    if partition_table:
        database["partition_temperature_k"] = (
            "partition_temperature",
            [100.0, 300.0, 700.0],
        )
        database["partition_function"] = ("partition_temperature", [1.5, 3.0, 4.0])
    else:
        database["energy_cminv"] = ("state", [0.0, 100.0])
        database["statistical_weight"] = ("state", [1.0, 3.0])
    optical = LineResonance(db=database)
    # Includes a far-wing point exercising the stable asymptotic derivative.
    atmosphere, _, _, _ = _atmosphere(3, [500.0004, 499.95, 500.0])
    temperature = np.linspace(150.0, 550.0, len(atmosphere.temperature_k))
    atmosphere.temperature_k = temperature
    derivative = optical.optical_derivatives(atmosphere)["temperature_k"]
    step = 0.002
    atmosphere.temperature_k = temperature + step
    plus = optical.atmosphere_quantities(atmosphere)
    atmosphere.temperature_k = temperature - step
    minus = optical.atmosphere_quantities(atmosphere)
    for name in ("extinction", "ssa", "leg_coeff"):
        finite_difference = (getattr(plus, name) - getattr(minus, name)) / (2 * step)
        analytical = getattr(derivative, f"d_{name}")
        tolerance = np.max(np.abs(finite_difference)) * 1e-8
        np.testing.assert_allclose(
            analytical, finite_difference, rtol=1e-5, atol=tolerance
        )
    atmosphere.temperature_k = temperature
    standalone = optical.cross_section_derivatives(
        atmosphere.wavelengths_nm,
        atmosphere.model_geometry.altitudes(),
        temperature_k=temperature,
    )["temperature_k"]
    np.testing.assert_allclose(
        standalone.reshape(derivative.d_extinction.shape), derivative.d_extinction
    )


def test_wing_cutoff_removes_opacity_without_renormalizing_profile():
    wavelengths = np.array([500.2, 500.0, 499.8])
    optical = LineResonance(db=_database())
    truncated = optical.cross_sections(wavelengths, [0], temperature_k=200)
    complete = _evaluate(_database(), wavelengths)
    assert np.all(complete.extinction > 0)
    assert np.all(truncated.extinction[0, [0, 2]] == 0)
    np.testing.assert_allclose(truncated.extinction[0, 1], complete.extinction[0, 1])
    assert np.all(np.isfinite(truncated.ssa))
    assert np.all(np.isfinite(truncated.leg_coeff))


def test_upper_lifetime_restricts_elastic_return():
    database = _database(einstein_a_s=[2e7], upper_total_a_s=[1e8])
    result = _evaluate(database, [500.0])
    np.testing.assert_allclose(result.ssa, 0.2)
    invalid = _database(einstein_a_s=[1e8], upper_total_a_s=[2e7])
    with pytest.raises(ValueError, match="upper_total_a_s"):
        LineResonance(db=invalid)


def test_netcdf_file_preserves_spectroscopy(tmp_path):
    filename = tmp_path / "synthetic.nc"
    dataset = _database()
    dataset.to_netcdf(filename)
    optical = LineResonance(db_filepath=filename)
    xr.testing.assert_identical(optical.database, dataset)
    actual = optical.cross_sections([500.0], [0], temperature_k=200)
    expected = _evaluate(dataset, [500.0])
    np.testing.assert_allclose(actual.extinction, expected.extinction)


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_partial_return_atmosphere_and_radiance_weighting_functions(num_stokes):
    database = _database(
        wavelength_nm=[500.0, 500.0006],
        lower_energy_cminv=[0.0, 100.0],
        lower_j=[0.0, 1.0],
        upper_j=[1.0, 1.0],
        oscillator_strength=[0.3, 0.7],
        einstein_a_s=[1e7, 2e7],
        upper_total_a_s=[2e7, 8e7],
    ).drop_dims("state")
    database["energy_cminv"] = ("state", [0.0, 100.0])
    database["statistical_weight"] = ("state", [1.0, 3.0])
    atmosphere, config, geometry, altitudes = _atmosphere(
        num_stokes, [500.0, 500.0004], calculate_derivatives=True
    )
    optical = LineResonance(db=database)
    temperatures = np.linspace(180.0, 250.0, len(altitudes))
    density = np.linspace(8e7, 1e8, len(altitudes))
    atmosphere.temperature_k = temperatures.copy()
    constituent = sk.constituent.NumberDensityScatterer(
        optical, altitudes, density.copy()
    )
    atmosphere["metal"] = constituent
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.TangentAltitudeSolar(20000, 0.4, 200000, 0.6))
    engine = sk.Engine(config, geometry, viewing)
    baseline = engine.calculate_radiance(atmosphere)
    standalone = optical.cross_sections(
        atmosphere.wavelengths_nm, altitudes, temperature_k=temperatures
    )
    # The native adapter must supply scattering cross sections, not an albedo
    # that would silently be clipped to one by atmosphere normalization.
    np.testing.assert_allclose(atmosphere.storage.ssa, standalone.ssa)
    np.testing.assert_allclose(
        atmosphere.storage.total_extinction, density[:, None] * standalone.extinction
    )
    for variable, original, step in [
        ("number_density", density, 1e4),
        ("temperature_k", temperatures, 0.002),
    ]:
        for index in [2, 4]:
            perturbed = original.copy()
            perturbed[index] += step
            if variable == "number_density":
                constituent.number_density = perturbed
            else:
                atmosphere.temperature_k = perturbed
            plus = engine.calculate_radiance(atmosphere).radiance
            perturbed = original.copy()
            perturbed[index] -= step
            if variable == "number_density":
                constituent.number_density = perturbed
            else:
                atmosphere.temperature_k = perturbed
            minus = engine.calculate_radiance(atmosphere).radiance
            numerical = ((plus - minus) / (2 * step)).to_numpy()
            analytical = (
                baseline.wf_metal_number_density.isel(metal_altitude=index)
                if variable == "number_density"
                else baseline.wf_temperature_k.isel(altitude=index)
            ).to_numpy()
            np.testing.assert_allclose(
                analytical,
                numerical,
                rtol=2e-5,
                atol=np.max(np.abs(numerical)) * 1e-7,
            )
            if variable == "number_density":
                constituent.number_density = original.copy()
            else:
                atmosphere.temperature_k = original.copy()


@pytest.mark.parametrize("rayleigh_background", [False, True])
def test_off_line_weighting_functions_vanish(rayleigh_background):
    atmosphere, config, geometry, altitudes = _atmosphere(
        3, [500.0, 501.0], calculate_derivatives=True
    )
    atmosphere["metal"] = sk.constituent.NumberDensityScatterer(
        LineResonance(db=_database(einstein_a_s=[2e7], upper_total_a_s=[1e8])),
        altitudes,
        np.full(len(altitudes), 1e8),
        temperature_k=np.full(len(altitudes), 200.0),
    )
    if rayleigh_background:
        atmosphere["rayleigh"] = sk.constituent.Rayleigh(
            method="manual",
            wavelengths_nm=np.array([499.0, 502.0]),
            xs=np.full(2, 1e-7),
            king_factor=np.ones(2),
        )
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.TangentAltitudeSolar(20000, 0.4, 200000, 0.6))
    result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)
    for name in ("wf_metal_temperature_k", "wf_metal_number_density"):
        assert np.all(np.isfinite(result[name]))
        np.testing.assert_allclose(result[name].isel(wavelength=1), 0.0, atol=1e-30)


def test_evaluation_rejects_wavelengths_outside_prepared_spectroscopy():
    database = _database()
    database.attrs.update(wavelength_min_nm=490.0, wavelength_max_nm=510.0)
    optical = LineResonance(db=database)
    for outside in [489.99, 510.01]:
        with pytest.raises(ValueError, match=r"[Ww]avelength"):
            optical.cross_sections([500.0, outside], [0.0], temperature_k=200)
    result = optical.cross_sections([490.0, 500.0, 510.0], [0.0], temperature_k=200)
    assert np.all(np.isfinite(result.extinction))


@pytest.mark.parametrize("num_stokes", [1, 3])
@pytest.mark.parametrize("temperature_source", ["constituent", "state"])
def test_temperature_weighting_function_on_distinct_constituent_grid(
    num_stokes, temperature_source
):
    atmosphere, config, geometry, native_altitudes = _atmosphere(
        num_stokes, [500.0, 500.0004], calculate_derivatives=True
    )
    altitudes = np.array([0.0, 20000.0, 40000.0, 60000.0])
    temperature = (
        np.array([180.0, 210.0, 230.0, 270.0])
        if temperature_source == "constituent"
        else np.linspace(180.0, 270.0, len(native_altitudes))
    )
    if temperature_source == "state":
        atmosphere.temperature_k = temperature.copy()
        # The implicit state Jacobian combines all constituents' responses.
        atmosphere["background"] = sk.constituent.Rayleigh(
            method="manual",
            wavelengths_nm=np.array([499.0, 501.0]),
            xs=np.full(2, 2e-18),
            king_factor=np.ones(2),
        )
    kwargs = (
        {"temperature_k": temperature.copy()}
        if temperature_source == "constituent"
        else {}
    )
    constituent = sk.constituent.NumberDensityScatterer(
        LineResonance(db=_database()),
        altitudes,
        np.array([1e7, 1e8, 4e7, 3e8]),
        **kwargs,
    )
    atmosphere["metal"] = constituent
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.TangentAltitudeSolar(20000, 0.4, 200000, 0.6))
    engine = sk.Engine(config, geometry, viewing)
    baseline = engine.calculate_radiance(atmosphere)
    for index in [1, 2]:
        step = np.zeros_like(temperature)
        step[index] = 0.002
        if temperature_source == "constituent":
            constituent.temperature_k = temperature + step
        else:
            atmosphere.temperature_k = temperature + step
        plus = engine.calculate_radiance(atmosphere).radiance
        if temperature_source == "constituent":
            constituent.temperature_k = temperature - step
        else:
            atmosphere.temperature_k = temperature - step
        minus = engine.calculate_radiance(atmosphere).radiance
        numerical = ((plus - minus) / 0.004).to_numpy()
        analytical = (
            baseline.wf_metal_temperature_k.isel(metal_altitude=index)
            if temperature_source == "constituent"
            else baseline.wf_temperature_k.isel(altitude=index)
        ).to_numpy()
        np.testing.assert_allclose(analytical, numerical, rtol=2e-6, atol=1e-15)


@pytest.mark.parametrize(
    "temperature_source", ["native", "altitude", "scalar", "state", "state_profile"]
)
def test_native_2d_temperature_weighting_functions(temperature_source):
    _, config, _, _ = _atmosphere(3)
    altitudes = np.array([0.0, 10000.0, 30000.0])
    geometry = sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6372000.0,
        altitude_grid_m=altitudes,
        horizontal_angle_grid_radians=np.array([-0.3, 0.3]),
    )
    atmosphere = sk.Atmosphere(
        geometry, config, wavelengths_nm=np.array([500.0, 500.0003])
    )
    density = np.array([[1e7, 2e7, 0.0], [3e7, 0.0, 5e7]])
    temperature = {
        "native": np.array([[180.0, 210.0, 230.0], [260.0, 190.0, 205.0]]),
        "altitude": np.array([180.0, 210.0, 230.0]),
        "scalar": np.array(210.0),
        "state": np.array([[180.0, 210.0, 230.0], [260.0, 190.0, 205.0]]),
        "state_profile": np.array([180.0, 210.0, 230.0]),
    }[temperature_source]
    state_source = temperature_source in ("state", "state_profile")
    atmosphere.temperature_k = (
        temperature.copy()
        if state_source
        else np.broadcast_to(temperature, density.shape).copy()
    )
    # Density Jacobians at absent-metal nodes require a scattering background.
    atmosphere.pressure_pa = np.ones_like(altitudes)
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    kwargs = {} if state_source else {"temperature_k": temperature}
    constituent = sk.constituent.NumberDensityScatterer2D(
        LineResonance(db=_database()), density, **kwargs
    )
    atmosphere["metal"] = constituent
    viewing = sk.ViewingGeometry()
    viewing.add_ray(
        sk.TangentAltitude(
            tangent_altitude_m=5000,
            observer_altitude_m=100000,
            horizontal_angle_radians=-0.1,
            viewing_azimuth_radians=np.pi / 2,
        )
    )
    engine = sk.Engine(config, geometry, viewing)
    baseline = engine.calculate_radiance(atmosphere)
    derivative = (
        baseline.wf_temperature_k if state_source else baseline.wf_metal_temperature_k
    )
    assert np.all(np.isfinite(derivative))
    index = {
        "native": (0, 1),
        "altitude": (1,),
        "scalar": (),
        "state": (0, 1),
        "state_profile": (1,),
    }[temperature_source]
    step = np.zeros_like(temperature)
    step[index] = 0.002
    if state_source:
        atmosphere.temperature_k = temperature + step
    else:
        constituent.temperature_k = temperature + step
    plus = engine.calculate_radiance(atmosphere).radiance
    if state_source:
        atmosphere.temperature_k = temperature - step
    else:
        constituent.temperature_k = temperature - step
    minus = engine.calculate_radiance(atmosphere).radiance
    spatial_dims = [
        dim for dim in derivative.dims if dim not in ("wavelength", "los", "stokes")
    ]
    assert len(spatial_dims) == len(index)
    analytical = derivative.isel(dict(zip(spatial_dims, index, strict=True)))
    np.testing.assert_allclose(
        analytical, (plus - minus) / 0.004, rtol=2e-6, atol=1e-15
    )
    if temperature_source == "native":
        # Changing temperature cannot produce opacity where this species is absent.
        np.testing.assert_allclose(
            derivative.isel(horizontal_angle=0, altitude=2), 0.0, atol=1e-30
        )


@pytest.mark.parametrize("num_stokes", [1, 3])
@pytest.mark.parametrize("zero_profile", [False, True])
@pytest.mark.parametrize("absorbing_background", [False, True])
def test_density_jacobian_rejects_zero_total_scattering(
    num_stokes, zero_profile, absorbing_background
):
    class AbsorbingBackground(sk.constituent.base.Constituent):
        def add_to_atmosphere(self, atmo):
            atmo.storage.total_extinction[:] += 1e-7

        def register_derivative(self, atmo, name):
            pass

    atmosphere, _, _, altitudes = _atmosphere(num_stokes, calculate_derivatives=True)
    density = np.zeros_like(altitudes) if zero_profile else np.full_like(altitudes, 1e8)
    density[3] = 0.0
    atmosphere["metal"] = sk.constituent.NumberDensityScatterer(
        LineResonance(db=_database()), altitudes, density
    )
    if absorbing_background:
        atmosphere["absorption"] = AbsorbingBackground()
    with pytest.raises(
        ValueError, match="wf_metal_number_density: total scattering is zero"
    ):
        atmosphere.internal_object()


def test_native_2d_density_jacobian_rejects_zero_total_scattering():
    _, config, _, _ = _atmosphere(3)
    geometry = sk.Geometry2D(
        cos_sza=0.6,
        solar_azimuth=0.0,
        earth_radius_m=6372000.0,
        altitude_grid_m=np.array([0.0, 10000.0, 30000.0]),
        horizontal_angle_grid_radians=np.array([-0.3, 0.3]),
    )
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=np.array([500.0]))
    atmosphere.temperature_k = np.full((2, 3), 200.0)
    density = np.full((2, 3), 1e8)
    density[1, 1] = 0.0
    atmosphere["metal"] = sk.constituent.NumberDensityScatterer2D(
        LineResonance(db=_database()), density
    )
    with pytest.raises(
        ValueError, match="wf_metal_number_density: total scattering is zero"
    ):
        atmosphere.internal_object()


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_zero_density_with_background_has_correct_one_sided_jacobian(num_stokes):
    atmosphere, config, geometry, altitudes = _atmosphere(
        num_stokes, calculate_derivatives=True
    )
    atmosphere.pressure_pa = np.ones_like(altitudes)
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    density = np.zeros_like(altitudes)
    constituent = sk.constituent.NumberDensityScatterer(
        LineResonance(db=_database()), altitudes, density
    )
    atmosphere["metal"] = constituent
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.TangentAltitudeSolar(20000, 0.4, 200000, 0.6))
    engine = sk.Engine(config, geometry, viewing)
    baseline = engine.calculate_radiance(atmosphere)
    step = density.copy()
    step[3] = 1e3
    constituent.number_density = step
    perturbed = engine.calculate_radiance(atmosphere)
    numerical = (perturbed.radiance - baseline.radiance) / step[3]
    analytical = baseline.wf_metal_number_density.isel(metal_altitude=3)
    assert np.max(np.abs(numerical)) > 0
    np.testing.assert_allclose(analytical, numerical, rtol=1e-4, atol=1e-17)


def test_zero_density_forward_only_and_unsupported_profile_locations_are_valid():
    atmosphere, _, _, altitudes = _atmosphere(1)
    atmosphere["metal"] = sk.constituent.NumberDensityScatterer(
        LineResonance(db=_database()), altitudes, np.zeros_like(altitudes)
    )
    atmosphere.internal_object()
    np.testing.assert_array_equal(atmosphere.storage.total_extinction, 0.0)

    atmosphere, _, _, altitudes = _atmosphere(1, calculate_derivatives=True)
    atmosphere["metal"] = sk.constituent.NumberDensityScatterer(
        LineResonance(db=_database()),
        altitudes[2:5],
        np.full(3, 1e8),
        out_of_bounds_mode="zero",
    )
    atmosphere.internal_object()
    assert np.all(np.isfinite(atmosphere.storage.total_extinction))


def test_local_database_loads_separate_emission_templates(tmp_path):
    folder = tmp_path / "spectroscopy" / "metals" / "emission"
    folder.mkdir(parents=True)
    dataset = xr.Dataset(
        {"relative_photon_spectrum": ("wavelength_nm", [0.0, 1.0, 0.0])},
        coords={"wavelength_nm": [500.0, 600.0, 700.0]},
    )
    dataset.to_netcdf(folder / "FeO.nc")
    database = sk.database.MetalSpectroscopyDatabase(db_root=tmp_path)
    assert database.available_species(kind="emission") == ["FeO"]
    assert database.available_species(kind="atomic") == []
    xr.testing.assert_identical(database.load_ds("FeO", kind="emission"), dataset)


@pytest.mark.parametrize("kind", ["atomic", "molecular", "emission"])
def test_default_database_downloads_missing_files_and_reuses_cache(
    tmp_path, monkeypatch, kind
):
    monkeypatch.setattr(sk.appconfig, "database_root", lambda: tmp_path)
    requests = []

    def download(_self, key):
        requests.append(key)
        target = tmp_path / key
        target.parent.mkdir(parents=True, exist_ok=True)
        _database().to_netcdf(target)
        return target

    monkeypatch.setattr(sk.database.StandardDatabase, "path", download)
    database = sk.database.MetalSpectroscopyDatabase()
    expected_key = f"spectroscopy/metals/{kind}/synthetic.nc"
    assert database.path("synthetic", kind=kind) == tmp_path / expected_key
    xr.testing.assert_identical(database.load_ds("synthetic", kind=kind), _database())
    assert requests == [expected_key]


def test_explicit_database_root_is_offline(tmp_path, monkeypatch):
    def unexpected_download(*args, **kwargs):
        pytest.fail("Explicit local databases must not download files")

    monkeypatch.setattr(sk.database.StandardDatabase, "path", unexpected_download)
    database = sk.database.MetalSpectroscopyDatabase(db_root=tmp_path)
    with pytest.raises(FileNotFoundError, match="not installed"):
        database.path("Na_I")


@pytest.mark.parametrize("auxiliary_override", [False, True])
def test_disabled_state_temperature_derivative_preserves_auxiliary_parameters(
    auxiliary_override,
):
    _, config, geometry, altitudes = _atmosphere(1)
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=np.array([500.0]),
        temperature_derivative=False,
    )
    temperature = np.full(len(altitudes), 200.0)
    atmosphere.temperature_k = temperature
    kwargs = {"temperature_k": temperature} if auxiliary_override else {}
    atmosphere["metal"] = sk.constituent.NumberDensityScatterer(
        LineResonance(db=_database()),
        altitudes,
        np.full(len(altitudes), 1e8),
        **kwargs,
    )
    atmosphere.internal_object()
    mappings = atmosphere.storage.derivative_mapping_names()
    assert ("wf_metal_temperature_k" in mappings) == auxiliary_override


@pytest.mark.parametrize("distinct_grid", [False, True])
@pytest.mark.parametrize("mixed_scattering", [False, True])
def test_extinction_normalized_resonance_retains_temperature_conversion_derivative(
    distinct_grid, mixed_scattering
):
    atmosphere, config, geometry, altitudes = _atmosphere(
        3, [500.0, 500.0004], calculate_derivatives=True
    )
    if distinct_grid:
        altitudes = altitudes[::2]
    database = _database()
    if mixed_scattering:
        database = _database(
            wavelength_nm=[500.0, 500.0003],
            lower_j=[0.5, 0.5],
            upper_j=[0.5, 1.5],
            einstein_a_s=[1e6, 5e5],
            upper_total_a_s=[1e6, 1e6],
        )
        database["statistical_weight"][:] = 2.0
        atmosphere.pressure_pa = np.ones_like(geometry.altitudes())
        atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    temperature = np.linspace(180.0, 270.0, len(altitudes))
    constituent = sk.constituent.ExtinctionScatterer(
        LineResonance(db=database),
        altitudes,
        np.linspace(1e-8, 2e-8, len(altitudes)),
        500.0,
        temperature_k=temperature.copy(),
    )
    atmosphere["metal"] = constituent
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.TangentAltitudeSolar(20000, 0.4, 200000, 0.6))
    engine = sk.Engine(config, geometry, viewing)
    baseline = engine.calculate_radiance(atmosphere)
    linearization = engine.linearize(atmosphere)
    tangent = xr.Dataset(
        {"metal_temperature_k": ("metal_altitude", np.linspace(-1, 1, len(altitudes)))}
    )
    expected_jvp = (baseline.wf_metal_temperature_k * tangent.metal_temperature_k).sum(
        "metal_altitude"
    )
    xr.testing.assert_allclose(linearization.jvp(tangent), expected_jvp)
    cotangent = xr.ones_like(baseline.radiance)
    expected_vjp = (baseline.wf_metal_temperature_k * cotangent).sum(
        baseline.radiance.dims
    )
    xr.testing.assert_allclose(
        linearization.vjp(cotangent).metal_temperature_k, expected_vjp
    )
    for index in [1, 2]:
        step = np.zeros_like(temperature)
        step[index] = 0.002
        constituent.temperature_k = temperature + step
        plus = engine.calculate_radiance(atmosphere).radiance
        constituent.temperature_k = temperature - step
        minus = engine.calculate_radiance(atmosphere).radiance
        derivative = baseline.wf_metal_temperature_k.isel(metal_altitude=index)
        np.testing.assert_allclose(
            derivative, (plus - minus) / 0.004, rtol=2e-6, atol=1e-14
        )
        # Cross-section changes are canceled by the reference-extinction
        # normalization at 500 nm, including its temperature-dependent density.
        if not distinct_grid and not mixed_scattering:
            np.testing.assert_allclose(derivative.isel(wavelength=0), 0.0, atol=1e-14)
