from __future__ import annotations

import hashlib
import io
from unittest.mock import Mock

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from sasktran2.climatology import ers
from sasktran2.constants import K_BOLTZMANN
from sasktran2.database import ers as ers_database
from sasktran2.optical.base import OpticalProperty, OpticalQuantities

SCENARIO = {"month": 7, "latitude_degrees": 45, "local_time_hours": 9.5}


@pytest.fixture()
def raw_dataset():
    """Distinct values on every axis catch swapped or ignored selectors."""
    coords = {
        "month": [1, 4, 7, 10],
        "hour": [9.5, 21.5],
        "lev": [0, 10, 20],
        "lat": [-80, -45, 0, 45, 80],
        "solmin_solmax": [0, 1],
        "volcanism": [0, 1],
    }
    data = xr.Dataset(coords=coords, attrs={"DOI": "10.5281/zenodo.10022129"})
    data.lev.attrs["units"] = "km"
    families = {
        "o3": ("month", "hour", "lev", "lat"),
        "h2o": ("month", "hour", "lev", "lat"),
        "o": ("month", "hour", "lev", "lat", "solmin_solmax"),
        "no": ("month", "hour", "lev", "lat", "solmin_solmax"),
        "cfc11": ("month", "lev", "lat"),
        "so2": ("month", "volcanism", "lev", "lat"),
        "h2so4m_c": ("month", "volcanism", "lev", "lat"),
    }
    for name, dims in families.items():
        shape = tuple(len(coords[d]) for d in dims)
        values = (np.arange(np.prod(shape)).reshape(shape) + 1) * 1e-9
        units = "ug/m3" if name == "h2so4m_c" else "mol mol-1"
        data[f"{name}_mean"] = xr.DataArray(values, dims=dims, attrs={"units": units})
        data[f"{name}_std"] = xr.DataArray(
            values / 10, dims=dims, attrs={"units": units}
        )
    dims = ("month", "hour", "lev", "lat", "solmin_solmax")
    shape = tuple(len(coords[d]) for d in dims)
    temperature = np.arange(np.prod(shape)).reshape(shape) / 10 + 220
    pressure = np.broadcast_to(
        np.array([1e5, 1e4, 1e3])[None, None, :, None, None], shape
    )
    data["temperature_mean"] = xr.DataArray(
        temperature, dims=dims, attrs={"units": "K"}
    )
    data["temperature_std"] = xr.DataArray(
        np.ones(shape), dims=dims, attrs={"units": "K"}
    )
    data["pressure_mean"] = xr.DataArray(
        pressure.copy(), dims=dims, attrs={"units": "Pa"}
    )
    data["surface_pressure_mean"] = data.pressure_mean.isel(lev=0, drop=True) * 0.8
    data.surface_pressure_mean.attrs["units"] = "Pa"
    data["airmolmass"] = xr.full_like(data.temperature_mean, 28.9)
    data.airmolmass.attrs["units"] = "gr/mol"
    data["lat_bnds"] = (
        ("bnds", "lat"),
        [[-90, -55, -20, 35, 70], [-70, -35, 20, 55, 90]],
    )
    return data


@pytest.fixture()
def local_file(raw_dataset, tmp_path):
    path = tmp_path / "ers.nc"
    raw_dataset.to_netcdf(path)
    return path


class _Absorber(OpticalProperty):
    def atmosphere_quantities(self, atmo, **kwargs):
        cross_sections = np.full((atmo.num_locations, atmo.num_wavel), 1e-22)
        return OpticalQuantities(
            extinction=cross_sections, ssa=np.zeros_like(cross_sections)
        )


def _atmosphere(altitudes=(0, 5000, 10000, 15000, 20000)):
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.single_scatter_source = sk.SingleScatterSource.Exact
    geometry = sk.Geometry1D(
        0.6,
        0,
        6372000,
        np.array(altitudes, dtype=float),
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    return sk.Atmosphere(geometry, config, wavelengths_nm=np.array([500.0])), config


def test_local_source_preserves_raw_values_and_closes_file(local_file):
    result = ers.load_dataset(path=local_file)
    local_file.unlink()
    assert result.lev.attrs["units"] == "km"
    assert result.hour.dtype.kind in "fi"
    assert result.o3_mean.shape == (4, 2, 3, 5)
    assert result.attrs["source_mode"] == "local"
    assert result.attrs["source_md5"]
    assert "2030" in str(result.attrs)


def test_exact_profile_all_dimension_families(local_file, raw_dataset):
    result = ers.profile(
        **SCENARIO,
        solar_activity="maximum",
        volcanic_activity="enhanced",
        path=local_file,
    )
    expected = raw_dataset.sel(month=7, hour=9.5, lat=45, solmin_solmax=1, volcanism=1)
    for name in (
        "o3_mean",
        "o3_std",
        "cfc11_mean",
        "no_mean",
        "so2_mean",
        "h2so4m_c_mean",
    ):
        np.testing.assert_array_equal(result[name], expected[name])
        assert result[name].dims == ("altitude_m",)
        assert result[name].attrs["units"] == raw_dataset[name].attrs["units"]
    np.testing.assert_array_equal(result.altitude_m, [0, 10000, 20000])
    np.testing.assert_array_equal(result.temperature_k, expected.temperature_mean)
    np.testing.assert_array_equal(result.pressure_pa, expected.pressure_mean)
    np.testing.assert_array_equal(result.lat_bnds, [35, 55])
    assert result.surface_pressure_pa == 80000
    assert result.attrs["solar_activity"] == "maximum"
    assert result.attrs["selected_month"] == 7
    assert "requires_verification" in result.airmolmass.attrs["quality_flags"]


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("month", 2),
        ("month", np.nan),
        ("latitude_degrees", 90),
        ("local_time_hours", 12),
        ("solar_activity", "average"),
        ("volcanic_activity", "low"),
    ],
)
def test_no_silent_nearest_selection(local_file, key, value):
    with pytest.raises(ValueError, match="ERS"):
        ers.profile(**{**SCENARIO, key: value}, path=local_file)


def test_species_alias_and_absent_dimensions(local_file):
    result = ers.profile(**SCENARIO, species="F11", path=local_file)
    other = ers.profile(
        **{**SCENARIO, "local_time_hours": 21.5},
        species="CFCl3",
        solar_activity="maximum",
        volcanic_activity="enhanced",
        path=local_file,
    )
    np.testing.assert_array_equal(result.cfc11_mean, other.cfc11_mean)
    assert "o3_mean" not in result
    assert "h2so4m_c_mean" not in result
    vmr = ers.constituent("F11", _Absorber(), **SCENARIO, path=local_file)
    np.testing.assert_array_equal(vmr.vmr, result.cfc11_mean)
    with pytest.raises(ValueError, match="Unsupported ERS gas"):
        ers.constituent("H2SO4", _Absorber(), **SCENARIO, path=local_file)
    with pytest.raises(ValueError, match="does not contain"):
        ers.profile(**SCENARIO, species="HDO", path=local_file)


def test_raw_quality_flags_and_explicit_clipping(raw_dataset, local_file):
    raw_dataset.h2o_mean.loc[{"lev": 20}] = -1e-8
    raw_dataset.o_mean[:] = 0
    raw_dataset.to_netcdf(local_file)
    result = ers.profile(**SCENARIO, path=local_file)
    assert result.h2o_mean[-1] == -1e-8
    assert result.h2o_mean.attrs["quality_flags"] == "negative_values"
    assert result.o_mean.attrs["quality_flags"] == "all_zero_profile"
    with pytest.raises(ValueError, match="negative VMRs"):
        ers.constituent("H2O", _Absorber(), **SCENARIO, path=local_file)
    with pytest.warns(UserWarning, match="Clipped 1 negative ERS h2o"):
        vmr = ers.constituent(
            "H2O", _Absorber(), **SCENARIO, negative_vmr="clip", path=local_file
        )
    assert vmr.vmr[-1] == 0
    assert ers.load_dataset(path=local_file).h2o_mean.min() < 0
    with pytest.raises(ValueError, match="all-zero atomic-oxygen"):
        ers.constituent("O", _Absorber(), **SCENARIO, path=local_file)


def test_add_to_atmosphere_interpolation_and_single_load(local_file, monkeypatch):
    loader = Mock(wraps=ers.load_dataset)
    monkeypatch.setattr(ers, "load_dataset", loader)
    atmosphere, _ = _atmosphere()
    ers.add_to_atmosphere(
        atmosphere, {"O3": _Absorber(), "F11": _Absorber()}, **SCENARIO, path=local_file
    )
    assert loader.call_count == 1
    np.testing.assert_allclose(
        atmosphere.pressure_pa, 10 ** np.array([5, 4.5, 4, 3.5, 3])
    )
    native = ers.profile(**SCENARIO, path=local_file)
    target = atmosphere.model_geometry.altitudes()
    np.testing.assert_allclose(
        atmosphere.temperature_k,
        np.interp(target, native.altitude_m, native.temperature_k),
    )
    atmosphere.internal_object()
    vmr = np.interp(target, native.altitude_m, native.o3_mean + native.cfc11_mean)
    expected = (
        vmr * 1e-22 * atmosphere.pressure_pa / (K_BOLTZMANN * atmosphere.temperature_k)
    )
    np.testing.assert_allclose(atmosphere.storage.total_extinction[:, 0], expected)


def test_only_relevant_altitudes_are_validated(raw_dataset, local_file):
    raw_dataset.h2o_mean.loc[{"lev": 20}] = -1e-8
    raw_dataset.to_netcdf(local_file)
    atmosphere, _ = _atmosphere((0, 5000, 10000))
    ers.add_to_atmosphere(atmosphere, {"H2O": _Absorber()}, **SCENARIO, path=local_file)
    np.testing.assert_array_equal(atmosphere["H2O"].altitudes_m, [0, 10000])
    atmosphere, _ = _atmosphere((0, 5000, 11000))
    # The negative upper bracketing level is still needed at 11 km.
    with pytest.raises(ValueError, match="negative VMRs"):
        ers.add_to_atmosphere(
            atmosphere, {"H2O": _Absorber()}, **SCENARIO, path=local_file
        )


@pytest.mark.parametrize(
    ("variable", "value", "message"),
    [
        ("h2o_mean", -1, "negative VMRs"),
        ("h2o_mean", np.nan, "nonfinite"),
        ("h2o_mean", 2, "greater than one"),
        ("pressure_mean", 0, "finite and positive"),
        ("pressure_mean", 1e5, "decrease with altitude"),
        ("temperature_mean", -1, "finite and positive"),
    ],
)
def test_validation_precedes_atmosphere_mutation(
    raw_dataset, local_file, variable, value, message
):
    raw_dataset[variable][:] = value
    raw_dataset.to_netcdf(local_file)
    atmosphere, _ = _atmosphere()
    atmosphere.temperature_k = np.full(5, 270.0)
    atmosphere.pressure_pa = np.full(5, 90000.0)
    original = ers.constituent("O3", _Absorber(), **SCENARIO, path=local_file)
    atmosphere["O3"] = original
    with pytest.raises(ValueError, match=message):
        ers.add_to_atmosphere(
            atmosphere,
            {"O3": _Absorber(), "H2O": _Absorber()},
            **SCENARIO,
            path=local_file,
        )
    np.testing.assert_array_equal(atmosphere.temperature_k, np.full(5, 270.0))
    np.testing.assert_array_equal(atmosphere.pressure_pa, np.full(5, 90000.0))
    assert atmosphere["O3"] is original
    assert atmosphere["H2O"] is None


def test_bounds_state_only_and_existing_humidity(local_file):
    atmosphere, _ = _atmosphere((-1000, 10000, 21000))
    with pytest.raises(ValueError, match="outside ERS range"):
        ers.add_to_atmosphere(atmosphere, {}, **SCENARIO, path=local_file)
    ers.add_to_atmosphere(
        atmosphere, {}, **SCENARIO, out_of_bounds_mode="extend", path=local_file
    )
    np.testing.assert_allclose(atmosphere.pressure_pa, [1e5, 1e4, 1e3])
    atmosphere.specific_humidity = np.full(3, 0.01)
    with pytest.raises(ValueError, match="specific_humidity"):
        ers.add_to_atmosphere(
            atmosphere, {"O3": _Absorber()}, **SCENARIO, path=local_file
        )
    atmosphere.specific_humidity = np.zeros(3)
    old_temperature = atmosphere.temperature_k.copy()
    old_pressure = atmosphere.pressure_pa.copy()
    ers.add_to_atmosphere(
        atmosphere,
        {"O3": _Absorber()},
        **SCENARIO,
        path=local_file,
        out_of_bounds_mode="extend",
        set_pressure_temperature=False,
    )
    np.testing.assert_array_equal(atmosphere.temperature_k, old_temperature)
    np.testing.assert_array_equal(atmosphere.pressure_pa, old_pressure)


@pytest.mark.parametrize(
    ("variable", "units"), [("lev", "m"), ("o3_mean", "ppmv"), ("pressure_mean", "hPa")]
)
def test_units_are_checked(raw_dataset, local_file, variable, units):
    raw_dataset[variable].attrs["units"] = units
    raw_dataset.to_netcdf(local_file)
    with pytest.raises(ValueError, match="must have units"):
        ers.profile(**SCENARIO, path=local_file)


def test_invalid_source_grid(raw_dataset, local_file):
    raw_dataset = raw_dataset.assign_coords(lev=[0, 20, 10])
    raw_dataset.lev.attrs["units"] = "km"
    raw_dataset.to_netcdf(local_file)
    with pytest.raises(ValueError, match="strictly increasing"):
        ers.profile(**SCENARIO, path=local_file)


def test_cache_download_reuse_and_corruption(local_file, tmp_path, monkeypatch):
    payload = local_file.read_bytes()
    digest = hashlib.md5(payload, usedforsecurity=False).hexdigest()
    monkeypatch.setitem(ers_database._FILES, "v07", ("CAIRT_ERS_v07.nc", digest))
    download = Mock(side_effect=lambda *args, **kwargs: io.BytesIO(payload))
    monkeypatch.setattr(ers_database.urllib.request, "urlopen", download)
    db = ers_database.ERSDatabase(db_root=tmp_path)
    path = db.path()
    assert path == tmp_path / "climatology/ers/v07/CAIRT_ERS_v07.nc"
    assert db.path() == path
    assert download.call_count == 1
    assert download.call_args.args[0].endswith("/CAIRT_ERS_v07.nc/content")
    ds = ers.load_dataset(db_root=tmp_path)
    assert ds.attrs["source_mode"] == "verified_cache"
    assert ds.attrs["source_md5"] == digest
    xr.testing.assert_equal(
        db.load_ds(), xr.load_dataset(local_file, decode_timedelta=False)
    )
    path.write_bytes(b"broken")
    with pytest.raises(OSError, match="checksum mismatch"):
        db.path()
    assert download.call_count == 1
    db.clear()
    assert not path.exists()


@pytest.mark.parametrize("failure", ["checksum", "interrupted"])
def test_failed_download_does_not_leave_a_cache(tmp_path, monkeypatch, failure):
    class InterruptedResponse(io.BytesIO):
        def read(self, *args, **kwargs):
            super().read(*args, **kwargs)
            msg = "interrupted download"
            raise OSError(msg)

    response = (
        io.BytesIO(b"wrong checksum")
        if failure == "checksum"
        else InterruptedResponse(b"partial")
    )
    monkeypatch.setattr(
        ers_database.urllib.request, "urlopen", Mock(return_value=response)
    )
    db = ers_database.ERSDatabase(db_root=tmp_path)
    with pytest.raises(OSError, match=r"checksum mismatch|interrupted download"):
        db.path()
    assert list((tmp_path / "climatology/ers/v07").iterdir()) == []


def test_invalid_version_and_source_options(local_file, tmp_path):
    with pytest.raises(ValueError, match="Unsupported ERS version"):
        ers.load_dataset("latest", path=local_file)
    with pytest.raises(ValueError, match="Unsupported ERS version"):
        ers_database.ERSDatabase("latest", db_root=tmp_path)
    with pytest.raises(ValueError, match="either path or db_root"):
        ers.load_dataset(path=local_file, db_root=tmp_path)


def test_radiance_and_vmr_derivative(local_file):
    atmosphere, config = _atmosphere()
    ers.add_to_atmosphere(atmosphere, {"O3": _Absorber()}, **SCENARIO, path=local_file)
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.GroundViewingSolar(0.6, 0, 0.8, 200000))
    engine = sk.Engine(config, atmosphere.model_geometry, viewing)
    result = engine.calculate_radiance(atmosphere)
    assert np.isfinite(result.radiance).all()
    assert (result.radiance > 0).all()
    original = atmosphere["O3"].vmr.copy()
    step = 1e-10
    atmosphere["O3"].vmr[1] += step
    plus = engine.calculate_radiance(atmosphere).radiance.to_numpy().copy()
    atmosphere["O3"].vmr[:] = original
    atmosphere["O3"].vmr[1] -= step
    minus = engine.calculate_radiance(atmosphere).radiance.to_numpy().copy()
    atmosphere["O3"].vmr[:] = original
    np.testing.assert_allclose(
        result.wf_O3_vmr.isel(O3_altitude=1),
        (plus - minus) / (2 * step),
        rtol=1e-5,
        atol=1e-8,
    )
