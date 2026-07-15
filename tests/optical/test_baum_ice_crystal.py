from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2._core_rust import PyMieIntegrator, WignerD
from scipy.special import roots_legendre

PARTICLE_MODELS = (
    "general_habit_mixture",
    "aggregate_solid_columns",
    "solid_columns",
)
COEFFICIENT_NAMES = ("lm_a1", "lm_a2", "lm_a3", "lm_a4", "lm_b1", "lm_b2")


def _reconstruct_baum_phase_matrix(
    p11: np.ndarray,
    p21_over_p11: np.ndarray,
    p22_over_p11: np.ndarray,
    p33_over_p11: np.ndarray,
    p43_over_p11: np.ndarray,
    p44_over_p11: np.ndarray,
) -> np.ndarray:
    return np.stack(
        (
            p11,
            p11 * p21_over_p11,
            p11 * p22_over_p11,
            p11 * p33_over_p11,
            -p11 * p43_over_p11,
            p11 * p44_over_p11,
        ),
        axis=1,
    )


def _direct_greek_transform(
    phase: np.ndarray,
    cos_angles: np.ndarray,
    angle_weights: np.ndarray,
    num_moments: int,
) -> np.ndarray:
    """Independent coefficient oracle used only for the Rust Mie convention test."""
    theta = np.arccos(cos_angles)
    wigners = {
        "00": WignerD(0, 0),
        "22": WignerD(2, 2),
        "2m2": WignerD(2, -2),
        "02": WignerD(0, 2),
    }
    coefficients = np.zeros((phase.shape[0], 6, num_moments))
    for order in range(num_moments):
        factor = (2 * order + 1) / 2
        basis = {
            name: factor * angle_weights * np.asarray(wigner.d(theta, order))
            for name, wigner in wigners.items()
        }
        coefficients[:, 0, order] = phase[:, 0] @ basis["00"]
        coefficients[:, 3, order] = phase[:, 5] @ basis["00"]
        coefficients[:, 4, order] = phase[:, 1] @ basis["02"]
        coefficients[:, 5, order] = -(phase[:, 4] @ basis["02"])
        plus = (phase[:, 2] + phase[:, 3]) @ basis["22"]
        minus = (phase[:, 2] - phase[:, 3]) @ basis["2m2"]
        coefficients[:, 1, order] = (plus + minus) / 2
        coefficients[:, 2, order] = (plus - minus) / 2

    coefficients /= coefficients[:, :1, :1]
    return coefficients


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
    baum_p43 = -p34 / p11
    baum_phase = _reconstruct_baum_phase_matrix(
        p11,
        p12 / p11,
        np.ones_like(p11),
        p33 / p11,
        baum_p43,
        p33 / p11,
    )
    np.testing.assert_allclose(p11 * baum_p43, -p34, rtol=0, atol=1e-14)
    np.testing.assert_allclose(baum_phase[:, 4], p34, rtol=0, atol=1e-14)

    transformed = _direct_greek_transform(
        baum_phase,
        cos_angles,
        angle_weights,
        num_moments,
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


@pytest.fixture()
def baum_runtime_database(tmp_path) -> Path:
    models = np.asarray(PARTICLE_MODELS, dtype=object)
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


def _write_moment_subset(
    source_path: Path, output_path: Path, max_moments: int
) -> None:
    with xr.open_dataset(source_path) as source:
        source.isel(legendre=slice(0, max_moments)).load().to_netcdf(output_path)


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


@pytest.mark.parametrize("particle_model", PARTICLE_MODELS)
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
    model_index = PARTICLE_MODELS.index(particle_model)
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
    _write_moment_subset(baum_runtime_database, default_database, max_moments=256)
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
