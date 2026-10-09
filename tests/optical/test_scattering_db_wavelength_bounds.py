from __future__ import annotations

import warnings

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2.optical.database import (
    OpticalDatabaseGenericScatterer,
    OpticalDatabaseGenericScattererRust,
)

DB_WAVELENGTHS_NM = np.array([1300.0, 1350.0, 1400.0])
DB_XS_M2 = np.array([3.0e-13, 2.0e-13, 1.0e-13])
DB_SSA = np.array([0.9, 0.8, 0.7])
ALTITUDES_M = np.arange(0.0, 65001.0, 5000.0)
RANGE_MESSAGE = "wavelength range of 1300 to 1400 nm"
WARNING_MESSAGE = "so the scatterer contributes nothing there"


def _scattering_dataset() -> xr.Dataset:
    # Henyey-Greenstein moments, with more terms than the atmosphere requests
    order = np.arange(32)
    a1 = np.tile((2 * order + 1) * 0.3**order, (len(DB_WAVELENGTHS_NM), 1))
    zeros = np.zeros_like(a1)

    dims = ("wavelength_nm", "legendre")
    return xr.Dataset(
        {
            "xs_total": (("wavelength_nm",), DB_XS_M2),
            "xs_scattering": (("wavelength_nm",), DB_XS_M2 * DB_SSA),
            "lm_a1": (dims, a1),
            "lm_a2": (dims, zeros),
            "lm_a3": (dims, zeros),
            "lm_a4": (dims, zeros),
            "lm_b1": (dims, zeros),
            "lm_b2": (dims, zeros),
        },
        coords={"wavelength_nm": DB_WAVELENGTHS_NM},
    )


def _rust_generic(tmp_path, **kwargs):  # noqa: ARG001
    return OpticalDatabaseGenericScattererRust(db=_scattering_dataset(), **kwargs)


def _xarray_generic(tmp_path, **kwargs):
    path = tmp_path / "scatterer.nc"
    _scattering_dataset().to_netcdf(path)
    return OpticalDatabaseGenericScatterer(path, **kwargs)


def _henyey_greenstein(tmp_path, **kwargs):  # noqa: ARG001
    return sk.optical.HenyeyGreenstein.from_parameters(
        wavelength_nm=DB_WAVELENGTHS_NM,
        xs_total=DB_XS_M2,
        ssa=DB_SSA,
        g=np.full_like(DB_XS_M2, 0.3),
        max_num_moments=16,
        **kwargs,
    )


def _mie_database(tmp_path, **kwargs):
    return sk.database.MieDatabase(
        sk.mie.distribution.LogNormalDistribution().freeze(
            mode_width=1.6, median_radius=80
        ),
        sk.mie.refractive.H2SO4(),
        DB_WAVELENGTHS_NM,
        db_root=tmp_path,
        **kwargs,
    )


SCATTERING_DATABASES = pytest.mark.parametrize(
    "make_db",
    [_rust_generic, _xarray_generic, _henyey_greenstein, _mie_database],
    ids=["rust_generic", "xarray_generic", "henyey_greenstein", "mie_database"],
)


def _atmosphere(wavelengths_nm: np.ndarray) -> sk.Atmosphere:
    config = sk.Config()
    config.num_streams = 4
    # The range checks must not depend on the engine's input validation
    config.input_validation_mode = sk.InputValidationMode.Disabled

    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6372000.0,
        ALTITUDES_M,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavelengths_nm)
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    return atmosphere


def _aerosol_extinction(db, wavelengths_nm: np.ndarray) -> np.ndarray:
    """Storage extinction of an atmosphere containing only an aerosol defined at 1350 nm"""
    atmosphere = _atmosphere(wavelengths_nm)
    atmosphere["aerosol"] = sk.constituent.ExtinctionScatterer(
        db, ALTITUDES_M, np.full_like(ALTITUDES_M, 1e-5), 1350.0
    )
    atmosphere.internal_object()
    return np.array(atmosphere.storage.total_extinction)


@SCATTERING_DATABASES
def test_wavelength_range_and_default_mode(make_db, tmp_path):
    db = make_db(tmp_path)

    assert db.wavelength_range_nm == (1300.0, 1400.0)
    assert db.wavelength_out_of_bounds_mode == "warn"


def test_wavelength_range_is_none_without_tabulated_grid():
    mie = sk.optical.Mie(
        sk.mie.distribution.LogNormalDistribution().freeze(
            mode_width=1.6, median_radius=80
        ),
        sk.mie.refractive.H2SO4(),
    )

    assert mie.wavelength_range_nm is None


@SCATTERING_DATABASES
def test_invalid_wavelength_out_of_bounds_mode(make_db, tmp_path):
    with pytest.raises(ValueError, match="wavelength_out_of_bounds_mode"):
        make_db(tmp_path, wavelength_out_of_bounds_mode="clip")

    db = make_db(tmp_path)
    with pytest.raises(ValueError, match="wavelength_out_of_bounds_mode"):
        db.wavelength_out_of_bounds_mode = "clip"
    assert db.wavelength_out_of_bounds_mode == "warn"


@SCATTERING_DATABASES
@pytest.mark.parametrize("mode", ["raise", "warn", "zero", "extend"])
def test_extinction_scatterer_rejects_out_of_range_extinction_wavelength(
    make_db, mode, tmp_path
):
    # The reference wavelength defines the whole profile, so it must be tabulated
    # even when out of range model wavelengths are accepted
    db = make_db(tmp_path, wavelength_out_of_bounds_mode=mode)

    with pytest.raises(ValueError, match="extinction_wavelength_nm=745 nm") as err:
        sk.constituent.ExtinctionScatterer(
            db, ALTITUDES_M, np.full_like(ALTITUDES_M, 1e-5), 745.0
        )
    assert RANGE_MESSAGE in str(err.value)


@SCATTERING_DATABASES
def test_extinction_scatterer_accepts_database_endpoints(make_db, tmp_path):
    db = make_db(tmp_path, wavelength_out_of_bounds_mode="raise")

    for wavelength_nm in DB_WAVELENGTHS_NM[[0, -1]]:
        constituent = sk.constituent.ExtinctionScatterer(
            db, ALTITUDES_M, np.full_like(ALTITUDES_M, 1e-5), wavelength_nm
        )
        factors = constituent._extinction_to_numden_factors
        assert np.all(np.isfinite(factors))
        assert np.all(factors > 0)


@SCATTERING_DATABASES
def test_out_of_range_atmosphere_wavelengths_warn_by_default(make_db, tmp_path):
    db = make_db(tmp_path)

    with pytest.warns(UserWarning, match=RANGE_MESSAGE) as record:
        extinction = _aerosol_extinction(db, np.array([470.0, 745.0, 1020.0, 1350.0]))

    messages = [str(w.message) for w in record if RANGE_MESSAGE in str(w.message)]
    assert "3 of 4 requested wavelengths (spanning 470 to 1020 nm)" in messages[0]
    assert WARNING_MESSAGE in messages[0]
    # Out of range wavelengths are zero filled, as with "zero"
    np.testing.assert_array_equal(extinction[:, :3], 0.0)
    np.testing.assert_allclose(extinction[:, 3], 1e-5, rtol=1e-6)


@SCATTERING_DATABASES
def test_out_of_range_warning_is_reported_once_per_calculation(make_db, tmp_path):
    # add_to_atmosphere and register_derivative both evaluate the optical property, but
    # users with the default warning filter should see one warning for the same wavelengths
    atmosphere = _atmosphere(np.array([1350.0, 1500.0]))
    assert atmosphere.calculate_derivatives
    atmosphere["aerosol"] = sk.constituent.ExtinctionScatterer(
        make_db(tmp_path), ALTITUDES_M, np.full_like(ALTITUDES_M, 1e-5), 1350.0
    )

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("default")
        atmosphere.internal_object()

    assert len([w for w in record if RANGE_MESSAGE in str(w.message)]) == 1


@SCATTERING_DATABASES
def test_out_of_range_atmosphere_wavelengths_raise(make_db, tmp_path):
    db = make_db(tmp_path, wavelength_out_of_bounds_mode="raise")

    with pytest.raises(ValueError, match=RANGE_MESSAGE) as err:
        _aerosol_extinction(db, np.array([470.0, 745.0, 1020.0, 1350.0]))
    assert "3 of 4 requested wavelengths (spanning 470 to 1020 nm)" in str(err.value)


@SCATTERING_DATABASES
@pytest.mark.parametrize("mode", ["raise", "warn"])
def test_number_density_scatterer_out_of_range_wavelengths(make_db, mode, tmp_path):
    atmosphere = _atmosphere(np.array([1350.0, 1500.0]))
    atmosphere["aerosol"] = sk.constituent.NumberDensityScatterer(
        make_db(tmp_path, wavelength_out_of_bounds_mode=mode),
        ALTITUDES_M,
        np.full_like(ALTITUDES_M, 1e6),
    )

    if mode == "raise":
        with pytest.raises(ValueError, match=RANGE_MESSAGE):
            atmosphere.internal_object()
    else:
        with pytest.warns(UserWarning, match=RANGE_MESSAGE):
            atmosphere.internal_object()


@SCATTERING_DATABASES
def test_in_range_atmosphere_wavelengths_unaffected(make_db, tmp_path):
    wavelengths_nm = np.array([1300.0, 1325.0, 1400.0])

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        default = _aerosol_extinction(make_db(tmp_path), wavelengths_nm)
        raised = _aerosol_extinction(
            make_db(tmp_path, wavelength_out_of_bounds_mode="raise"), wavelengths_nm
        )
    zero = _aerosol_extinction(
        make_db(tmp_path, wavelength_out_of_bounds_mode="zero"), wavelengths_nm
    )

    assert np.all(np.isfinite(default))
    assert np.all(default > 0)
    np.testing.assert_array_equal(default, zero)
    np.testing.assert_array_equal(default, raised)


@SCATTERING_DATABASES
def test_zero_mode_omits_out_of_range_wavelengths_silently(make_db, tmp_path):
    db = make_db(tmp_path, wavelength_out_of_bounds_mode="zero")

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        extinction = _aerosol_extinction(db, np.array([1200.0, 1350.0, 1500.0]))

    assert not [w for w in record if RANGE_MESSAGE in str(w.message)]

    np.testing.assert_array_equal(extinction[:, [0, 2]], 0.0)
    np.testing.assert_allclose(extinction[:, 1], 1e-5, rtol=1e-6)


@SCATTERING_DATABASES
def test_extend_mode_holds_nearest_database_wavelength(make_db, tmp_path):
    db = make_db(tmp_path, wavelength_out_of_bounds_mode="extend")

    extended = _aerosol_extinction(db, np.array([1200.0, 1500.0]))
    edges = _aerosol_extinction(db, np.array([1300.0, 1400.0]))

    np.testing.assert_allclose(extended, edges, rtol=1e-12)
    assert np.all(extended > 0)


@SCATTERING_DATABASES
def test_mode_can_be_changed_after_construction(make_db, tmp_path):
    db = make_db(tmp_path)
    wavelengths_nm = np.array([1350.0, 1500.0])

    db.wavelength_out_of_bounds_mode = "raise"
    with pytest.raises(ValueError, match=RANGE_MESSAGE):
        _aerosol_extinction(db, wavelengths_nm)

    db.wavelength_out_of_bounds_mode = "zero"
    assert db.wavelength_out_of_bounds_mode == "zero"
    np.testing.assert_array_equal(_aerosol_extinction(db, wavelengths_nm)[:, 1], 0.0)


@SCATTERING_DATABASES
def test_cross_sections_out_of_range(make_db, tmp_path):
    db = make_db(tmp_path)

    with pytest.warns(UserWarning, match=RANGE_MESSAGE):
        extinction = db.cross_sections(np.array([745.0]), ALTITUDES_M).extinction
    np.testing.assert_array_equal(extinction, 0.0)

    db.wavelength_out_of_bounds_mode = "raise"
    with pytest.raises(ValueError, match=RANGE_MESSAGE):
        db.cross_sections(np.array([745.0]), ALTITUDES_M)


def test_extinction_scatterer_rejects_zero_cross_section():
    # Zero cross section inside the tabulated range has no valid number density either
    db = sk.optical.HenyeyGreenstein.from_parameters(
        wavelength_nm=DB_WAVELENGTHS_NM,
        xs_total=np.array([1.0e-13, 0.0, 1.0e-13]),
        ssa=DB_SSA,
        g=np.full_like(DB_XS_M2, 0.3),
        max_num_moments=16,
    )

    with pytest.raises(ValueError, match="must be finite and positive") as err:
        sk.constituent.ExtinctionScatterer(
            db, ALTITUDES_M, np.full_like(ALTITUDES_M, 1e-5), 1350.0
        )
    message = str(err.value)
    assert "extinction_wavelength_nm=1350 nm" in message
    assert f"{len(ALTITUDES_M)} of {len(ALTITUDES_M)} values" in message
    assert "1300 to 1400 nm" in message


def test_gaussian_height_extinction_rejects_out_of_range_wavelength(tmp_path):
    atmosphere = _atmosphere(np.array([1350.0]))
    atmosphere["aerosol"] = sk.constituent.GaussianHeightExtinction(
        _rust_generic(tmp_path),
        height_m=20000.0,
        width_fwhm_m=5000.0,
        vertical_optical_depth=0.01,
        vertical_optical_depth_wavel_nm=745.0,
        altitudes_m=ALTITUDES_M,
    )

    with pytest.raises(ValueError, match="vertical_optical_depth_wavel_nm=745 nm"):
        atmosphere.internal_object()
