from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2._core_rust import PyMieIntegrator, WignerD
from scipy.special import roots_legendre
from tools.databases.baum_ice_crystals import (
    COEFFICIENT_NAMES,
    MODEL_FILES,
    ConvergenceError,
    _WignerBasisSequence,
    adaptive_transform_batch,
    build_transform_operator,
    calculate_greek_coefficients,
    create_runtime_subset,
    generate_database,
    reconstruct_absolute_phase_matrix,
    reconstruct_phase_matrix,
)


def _known_coefficients(num_moments: int = 20) -> np.ndarray:
    coefficients = np.zeros((6, num_moments), dtype=np.float64)
    coefficients[0, :3] = [1.0, 0.3, 0.05]
    coefficients[1, 2] = 0.2
    coefficients[2, 2] = 0.1
    coefficients[3, :3] = [0.8, 0.0, 0.03]
    coefficients[4, 2] = 0.05
    coefficients[5, 2] = -0.02
    return coefficients


def test_reconstruct_absolute_baum_phase_matrix_signs():
    p11 = np.array([2.0, 4.0])
    phase = reconstruct_absolute_phase_matrix(
        p11,
        np.array([0.1, -0.2]),
        np.array([0.8, 0.7]),
        np.array([0.6, 0.5]),
        np.array([0.3, -0.4]),
        np.array([0.9, 0.85]),
    )

    np.testing.assert_allclose(phase[0], p11)
    np.testing.assert_allclose(phase[1], p11 * np.array([0.1, -0.2]))
    np.testing.assert_allclose(phase[4], -p11 * np.array([0.3, -0.4]))


@pytest.mark.parametrize(("m", "n"), [(0, 0), (2, 2), (2, -2), (0, 2)])
def test_wigner_basis_sequence_matches_independent_binding(m: int, n: int):
    theta = np.linspace(0.0, np.pi, 31)
    sequence = _WignerBasisSequence(m, n, theta)
    binding = WignerD(m, n)

    actual = np.concatenate((sequence.block(0, 7), sequence.block(7, 16384)))
    for order in (0, 2, 39, 511, 2047, 8191, 16383):
        expected = np.asarray(binding.d(theta, order))
        np.testing.assert_allclose(actual[order], expected, rtol=2e-13, atol=2e-14)


def test_adaptive_transform_reconstructs_498_angle_phase_matrix():
    angles = np.linspace(0.0, 180.0, 498)
    coefficients = _known_coefficients()
    phase = reconstruct_phase_matrix(
        coefficients, np.cos(np.deg2rad(angles)), block_size=4
    )
    operator = build_transform_operator(
        angles,
        max_moments=20,
        block_size=4,
        quadrature_oversampling=2.0,
        minimum_quadrature_order=4,
    )

    result = adaptive_transform_batch(
        phase,
        operator,
        block_size=4,
        relative_rms_tolerance=1e-4,
        max_normalized_tolerance=1e-3,
    )

    assert result.required_moments[0] == 8
    assert np.all(result.converged)
    np.testing.assert_allclose(
        result.coefficients[0, :, :3], coefficients[:, :3], rtol=2e-5, atol=2e-6
    )
    assert np.all(result.rms_error <= 1e-4)
    assert np.all(result.max_normalized_error <= 1e-3)


@pytest.mark.parametrize(
    ("size_parameter", "refractive_index"),
    [(2.3, 1.33 - 1e-5j), (9.0, 1.5 - 0.02j)],
)
def test_baum_coefficient_signs_match_rust_mie(
    size_parameter: float, refractive_index: complex
):
    num_angles = 96
    num_moments = 32
    cos_angles, angle_weights = roots_legendre(num_angles)
    integrator = PyMieIntegrator(cos_angles, num_moments, 1)

    xs_total = np.zeros(1)
    xs_scattering = np.zeros(1)
    mie_phase = [np.zeros((1, num_angles)) for _ in range(4)]
    mie_coefficients = [np.zeros((1, num_moments)) for _ in range(6)]
    integrator.integrate(
        500.0,
        refractive_index,
        np.array([size_parameter]),
        np.ones((1, 1)),
        np.ones(1),
        angle_weights,
        xs_total,
        xs_scattering,
        *mie_phase,
        *mie_coefficients,
    )
    p11, p12, p33, p34 = mie_phase

    # Encode the Mie sphere in the Baum Pij/P11 convention, including P43=-P34.
    baum_phase = reconstruct_absolute_phase_matrix(
        p11,
        p12 / p11,
        np.ones_like(p11),
        p33 / p11,
        -p34 / p11,
        p33 / p11,
    )
    np.testing.assert_allclose(baum_phase[:, 4], p34, rtol=0, atol=1e-14)

    transformed = calculate_greek_coefficients(
        baum_phase,
        cos_angles,
        angle_weights,
        num_moments,
        block_size=8,
        normalize=True,
    )
    mie_coefficients = np.stack(mie_coefficients, axis=1)
    mie_coefficients /= mie_coefficients[:, 0:1, 0:1]

    # Independently evaluate the defining integrals for an order where d02 is
    # nonzero. This prevents the converter and Mie accumulator from passing with
    # mutually consistent but reversed off-diagonal signs.
    order = 2
    d02 = np.asarray(WignerD(0, 2).d(np.arccos(cos_angles), order))
    coefficient_weight = (2 * order + 1) / 2
    a1_normalization = np.sum(angle_weights * p11) / 2
    expected_b1 = (
        coefficient_weight * np.sum(angle_weights * d02 * p12) / a1_normalization
    )
    expected_b2 = (
        -coefficient_weight * np.sum(angle_weights * d02 * p34) / a1_normalization
    )
    np.testing.assert_allclose(transformed[0, 4, order], expected_b1, atol=1e-13)
    np.testing.assert_allclose(transformed[0, 5, order], expected_b2, atol=1e-13)

    # B1 and B2 are both appreciably nonzero, so this comparison detects either
    # off-diagonal sign being reversed.
    assert np.max(np.abs(mie_coefficients[:, 4])) > 1e-4
    assert np.max(np.abs(mie_coefficients[:, 5])) > 1e-4
    np.testing.assert_allclose(transformed, mie_coefficients, rtol=1e-10, atol=1e-12)


def _write_synthetic_baum_sources(directory: Path) -> None:
    angles = np.linspace(0.0, 180.0, 498)
    cos_angles = np.cos(np.deg2rad(angles))
    first = _known_coefficients(20)
    second = first.copy()
    second[0, 6] = 0.01
    coefficient_grid = np.stack((first, second))
    phase = reconstruct_phase_matrix(coefficient_grid, cos_angles, block_size=4)
    p11 = phase[:, 0]

    source_phase = {
        "p11_phase_function": p11,
        "p21_phase_function": phase[:, 1] / p11,
        "p22_phase_function": phase[:, 2] / p11,
        "p33_phase_function": phase[:, 3] / p11,
        "p43_phase_function": -phase[:, 4] / p11,
        "p44_phase_function": phase[:, 5] / p11,
    }
    variables = {
        name: (
            ("nScatAng", "nDeff", "nWaveLen"),
            np.transpose(values, (1, 0))[:, np.newaxis, :].astype(np.float32),
        )
        for name, values in source_phase.items()
    }
    variables.update(
        {
            "wavelengths": (("nWaveLen",), np.array([0.5, 0.7], dtype=np.float32)),
            "effective_diameter": (
                ("nDeff",),
                np.array([10.0], dtype=np.float32),
            ),
            "phase_angles": (("nScatAng",), angles.astype(np.float32)),
            "asymmetry_parameter": (
                ("nDeff", "nWaveLen"),
                np.array([[0.1, 0.1]], dtype=np.float32),
            ),
            "single_scattering_albedo": (
                ("nDeff", "nWaveLen"),
                np.array([[0.8, 0.9]], dtype=np.float32),
            ),
            "extinction_coefficient_over_iwc": (
                ("nDeff", "nWaveLen"),
                np.array([[2.0, 3.0]], dtype=np.float32),
            ),
        }
    )
    source = xr.Dataset(variables)
    for filename in MODEL_FILES.values():
        source.to_netcdf(directory / filename)


def test_generate_baum_database_schema_and_zero_padding(tmp_path):
    _write_synthetic_baum_sources(tmp_path)
    output = tmp_path / "baum.nc"
    generate_database(
        tmp_path,
        output,
        moment_block_size=4,
        safety_cap=16,
        wavelength_batch_size=2,
    )

    with xr.open_dataset(output) as database:
        assert tuple(database.particle_model.values) == tuple(MODEL_FILES)
        assert database.sizes["legendre"] == 12
        assert database["lm_a1"].dtype == np.float32
        assert database["lm_a1"].encoding["chunksizes"] == (1, 1, 2, 64)
        assert set(COEFFICIENT_NAMES).issubset(database.data_vars)
        np.testing.assert_allclose(
            database.xs_scattering, database.xs_total * [[0.8, 0.9]]
        )
        assert np.all(database.required_moments.isel(wavelength_nm=0) == 8)
        assert np.all(database.required_moments.isel(wavelength_nm=1) == 12)
        np.testing.assert_array_equal(
            database.lm_a1.isel(wavelength_nm=0, legendre=slice(8, None)), 0.0
        )
        assert float(database.phase_relative_rms_error.max()) <= 1e-4
        assert float(database.phase_max_normalized_error.max()) <= 1e-3
        assert "P34=-P43" in database.attrs["sasktran2_phase_matrix_convention"]
        assert database.attrs["database_variant"] == "full"
        assert database.attrs["stored_moments"] == 12
        assert database.attrs["full_available_moments"] == 12
        assert database.attrs["phase_reconstruction_all_converged"] == 1
        assert np.all(database.phase_reconstruction_converged == 1)
        settings = database.attrs["generation_settings"]
        assert '"max_normalized_tolerance": 0.001' in settings


def test_create_baum_runtime_subset_is_exact_prefix(tmp_path):
    _write_synthetic_baum_sources(tmp_path)
    full_path = tmp_path / "baum_full.nc"
    subset_path = tmp_path / "baum.nc"
    generate_database(
        tmp_path,
        full_path,
        moment_block_size=4,
        safety_cap=16,
        wavelength_batch_size=2,
    )

    create_runtime_subset(full_path, subset_path, max_moments=8)

    with xr.open_dataset(full_path) as full, xr.open_dataset(subset_path) as subset:
        assert subset.sizes["legendre"] == 8
        assert subset.attrs["database_variant"] == "default"
        assert subset.attrs["stored_moments"] == 8
        assert subset.attrs["full_available_moments"] == 12
        assert subset.attrs["standard_database_key"].endswith(
            "baum_ice_crystals_v3_6.nc"
        )
        assert subset.attrs["full_database_key"].endswith(
            "baum_ice_crystals_v3_6_full.nc"
        )
        np.testing.assert_array_equal(subset.required_moments, full.required_moments)
        for name in COEFFICIENT_NAMES:
            np.testing.assert_array_equal(
                subset[name], full[name].isel(legendre=slice(0, 8))
            )

    with pytest.raises(FileExistsError):
        create_runtime_subset(full_path, subset_path, max_moments=8)
    with pytest.raises(ValueError, match="only contains 12"):
        create_runtime_subset(full_path, tmp_path / "too_many.nc", max_moments=13)


def test_adaptive_transform_reports_safety_cap_failure():
    angles = np.linspace(0.0, 180.0, 498)
    coefficients = _known_coefficients()
    coefficients[0, 10] = 0.02
    phase = reconstruct_phase_matrix(coefficients, np.cos(np.deg2rad(angles)))
    operator = build_transform_operator(angles, max_moments=8, block_size=4)

    with pytest.raises(ConvergenceError, match="did not converge by 8 moments"):
        adaptive_transform_batch(phase, operator, block_size=4)

    capped = adaptive_transform_batch(
        phase, operator, block_size=4, allow_unconverged_at_cap=True
    )
    assert capped.required_moments[0] == 8
    assert not capped.converged[0]
    assert capped.coefficients.shape[-1] == 8


@pytest.fixture()
def baum_runtime_database(tmp_path) -> Path:
    models = np.asarray(tuple(MODEL_FILES), dtype=object)
    diameters = np.array([10.0, 20.0])
    wavelengths = np.array([500.0, 700.0])
    moments = np.arange(300)
    shape = (len(models), len(diameters), len(wavelengths), len(moments))
    xs_total = np.empty(shape[:-1])
    for model_index in range(len(models)):
        xs_total[model_index] = (
            1e-3
            + model_index * 1e-4
            + diameters[:, np.newaxis] * 1e-6
            + wavelengths[np.newaxis, :] * 1e-9
        )
    coefficients = {name: np.zeros(shape) for name in COEFFICIENT_NAMES}
    coefficients["lm_a1"][..., 0] = 1.0
    coefficients["lm_a1"][..., 1] = 0.3
    coefficients["lm_a2"][..., 2] = 0.2
    coefficients["lm_a3"][..., 2] = 0.1
    coefficients["lm_a4"][..., 0] = 0.8
    coefficients["lm_b1"][..., 2] = 0.05
    coefficients["lm_b2"][..., 2] = -0.02
    dimensions = (
        "particle_model",
        "effective_diameter_um",
        "wavelength_nm",
        "legendre",
    )
    database = xr.Dataset(
        {
            "xs_total": (dimensions[:-1], xs_total),
            "xs_scattering": (dimensions[:-1], 0.8 * xs_total),
            **{name: (dimensions, value) for name, value in coefficients.items()},
        },
        coords={
            "particle_model": models,
            "effective_diameter_um": diameters,
            "wavelength_nm": wavelengths,
            "legendre": moments,
        },
    )
    path = tmp_path / "baum_runtime.nc"
    database.to_netcdf(path)
    return path


def _atmosphere(num_moments: int, num_stokes: int = 1) -> sk.Atmosphere:
    config = sk.Config()
    config.num_singlescatter_moments = num_moments
    config.num_stokes = num_stokes
    altitudes = np.array([0.0, 1000.0, 2000.0])
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6372000.0,
        altitudes,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    return sk.Atmosphere(geometry, config, wavelengths_nm=np.array([500.0, 700.0]))


@pytest.mark.parametrize("particle_model", tuple(MODEL_FILES))
def test_baum_ice_crystal_model_selection_and_rust_interpolation(
    baum_runtime_database: Path, particle_model: str
):
    optical_property = sk.optical.BaumIceCrystal(
        particle_model=particle_model,
        db_filepath=baum_runtime_database,
    )
    assert optical_property.particle_model == particle_model
    assert optical_property.max_moments == 256
    assert optical_property.available_moments == 300

    diameters = np.array([12.0, 15.0, 18.0])
    quantities = optical_property.cross_sections(
        np.array([500.0]),
        np.arange(3.0),
        effective_diameter_um=diameters,
    )
    model_index = tuple(MODEL_FILES).index(particle_model)
    expected = 1e-3 + model_index * 1e-4 + diameters * 1e-6 + 500.0e-9
    np.testing.assert_allclose(quantities.extinction[:, 0], expected, rtol=1e-12)
    np.testing.assert_allclose(quantities.ssa[:, 0], 0.8 * expected, rtol=1e-12)

    derivatives = optical_property.cross_section_derivatives(
        np.array([500.0]),
        np.arange(3.0),
        effective_diameter_um=diameters,
    )
    np.testing.assert_allclose(derivatives["effective_diameter_um"], 1e-6)


def test_baum_ice_crystal_moment_limits_and_diameter_validation(
    baum_runtime_database: Path,
):
    with pytest.raises(ValueError, match="positive"):
        sk.optical.BaumIceCrystal(max_moments=0, db_filepath=baum_runtime_database)
    with pytest.raises(TypeError, match="integer or None"):
        sk.optical.BaumIceCrystal(max_moments=12.5, db_filepath=baum_runtime_database)
    with pytest.raises(ValueError, match="only contains 300"):
        sk.optical.BaumIceCrystal(max_moments=301, db_filepath=baum_runtime_database)

    full = sk.optical.BaumIceCrystal(
        max_moments=None, db_filepath=baum_runtime_database
    )
    assert full.max_moments == full.available_moments == 300

    limited = sk.optical.BaumIceCrystal(
        max_moments=4, db_filepath=baum_runtime_database
    )
    with pytest.raises(ValueError, match=r"within \[10.0, 20.0\]"):
        limited.cross_sections(
            np.array([500.0]),
            np.array([0.0]),
            effective_diameter_um=np.array([9.0]),
        )
    with pytest.raises(ValueError, match="Atmosphere requests 5"):
        limited.atmosphere_quantities(
            _atmosphere(5), effective_diameter_um=np.full(3, 15.0)
        )


def test_baum_standard_database_selects_artifact_from_moment_limit(
    baum_runtime_database: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    default_database = tmp_path / "baum_default.nc"
    create_runtime_subset(baum_runtime_database, default_database, max_moments=256)
    requested_keys = []

    def standard_database_path(_database, key):
        requested_keys.append(key)
        if key.endswith("_full.nc"):
            return baum_runtime_database
        return default_database

    monkeypatch.setattr(
        "sasktran2.optical.baum.StandardDatabase.path", standard_database_path
    )

    default = sk.optical.BaumIceCrystal()
    assert default.max_moments == default.available_moments == 256
    assert requested_keys[-1].endswith("baum_ice_crystals_v3_6.nc")

    limited = sk.optical.BaumIceCrystal(max_moments=128)
    assert limited.max_moments == 128
    assert limited.available_moments == 256
    assert requested_keys[-1].endswith("baum_ice_crystals_v3_6.nc")

    with pytest.warns(UserWarning, match="full Baum ice-crystal database"):
        expanded = sk.optical.BaumIceCrystal(max_moments=257)
    assert expanded.max_moments == 257
    assert expanded.available_moments == 300
    assert requested_keys[-1].endswith("baum_ice_crystals_v3_6_full.nc")

    with pytest.warns(UserWarning, match="full Baum ice-crystal database"):
        complete = sk.optical.BaumIceCrystal(max_moments=None)
    assert complete.max_moments == complete.available_moments == 300
    assert requested_keys[-1].endswith("baum_ice_crystals_v3_6_full.nc")


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_baum_ice_crystal_extinction_scatterer(
    baum_runtime_database: Path, num_stokes: int
):
    optical_property = sk.optical.BaumIceCrystal(
        max_moments=16, db_filepath=baum_runtime_database
    )
    atmosphere = _atmosphere(16, num_stokes)
    altitudes = atmosphere.model_geometry.altitudes()
    constituent = sk.constituent.ExtinctionScatterer(
        optical_property,
        altitudes,
        np.full_like(altitudes, 1e-6),
        500.0,
        effective_diameter_um=np.array([10.0, 15.0, 20.0]),
    )
    constituent.add_to_atmosphere(atmosphere)

    diameters = np.array([10.0, 15.0, 20.0])
    xs_500 = 1e-3 + diameters * 1e-6 + 500.0e-9
    xs_700 = 1e-3 + diameters * 1e-6 + 700.0e-9
    expected_extinction = np.column_stack((np.full(3, 1e-6), 1e-6 * xs_700 / xs_500))
    np.testing.assert_allclose(
        atmosphere.storage.total_extinction, expected_extinction, rtol=1e-12
    )
    assert np.all(np.isfinite(atmosphere.storage.leg_coeff))


def test_baum_ice_crystal_engine_smoke(baum_runtime_database: Path):
    config = sk.Config()
    config.num_streams = 4
    config.num_singlescatter_moments = 16
    altitudes = np.arange(0.0, 30001.0, 5000.0)
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6372000.0,
        altitudes,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    viewing_geometry = sk.ViewingGeometry()
    viewing_geometry.add_ray(sk.TangentAltitudeSolar(10000.0, 0.0, 100000.0, 0.6))
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=np.array([500.0]))
    atmosphere["ice"] = sk.constituent.ExtinctionScatterer(
        sk.optical.BaumIceCrystal(max_moments=16, db_filepath=baum_runtime_database),
        altitudes,
        np.full_like(altitudes, 1e-7),
        500.0,
        effective_diameter_um=np.full_like(altitudes, 15.0),
    )

    radiance = sk.Engine(config, geometry, viewing_geometry).calculate_radiance(
        atmosphere
    )
    assert np.all(np.isfinite(radiance.radiance))
