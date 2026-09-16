from __future__ import annotations

import hashlib

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from numpy.testing import assert_allclose, assert_array_equal
from sasktran2.climatology import stratospheric_aerosol as aerosol
from sasktran2.database import stratospheric_aerosol as database
from sasktran2.optical.base import OpticalProperty, OpticalQuantities
from scipy.integrate import trapezoid


@pytest.fixture()
def catalogue():
    z = np.arange(17000.0, 32001.0, 500.0)
    waves = np.array([384.0, 448.0, 520.0, 601.0, 676.0, 756.0, 869.0, 1021.0, 1543.0])
    names = [
        f"{band}_{tier}"
        for band in ("sh_midlat", "tropical", "nh_midlat")
        for tier in ("low", "typical", "elevated", "extreme")
    ]
    extinction = (
        1e-6
        * np.exp(-(((z - 22000) / 4000) ** 2))
        * (1 + 0.15 * np.cos((z - z[0]) / 500 * np.pi))
    )
    raw = (
        extinction[None, :, None]
        * np.arange(1.0, 13.0)[:, None, None]
        * (756 / waves)[None, None, :]
    )
    radius = np.broadcast_to(80 + (z - z[0]) / 500, (12, len(z))).copy()
    return xr.Dataset(
        {
            "raw_extinction_per_m": (
                ("scenario", "altitude_m", "wavelength_nm"),
                raw,
                {"units": "m-1"},
            ),
            "raw_extinction_uncertainty_per_m": (
                ("scenario", "altitude_m", "wavelength_nm"),
                raw / 10,
                {"units": "m-1"},
            ),
            "raw_median_radius_nm": (
                ("scenario", "altitude_m"),
                radius,
                {"units": "nm"},
            ),
            "raw_median_radius_uncertainty_nm": (
                ("scenario", "altitude_m"),
                radius / 10,
                {"units": "nm"},
            ),
            "observed_valid": (
                ("scenario", "altitude_m"),
                np.ones((12, len(z)), dtype=np.int8),
            ),
            "source_aerosol_flag": (
                ("scenario", "altitude_m", "wavelength_nm"),
                np.full_like(raw, 2),
            ),
            "reference_upper_scale_height_m": (
                "scenario",
                np.repeat([3000.0, 4000.0, 3000.0], 4),
                {"units": "m"},
            ),
            "observed_bottom_m": ("scenario", np.full(12, z[0]), {"units": "m"}),
            "observed_top_m": ("scenario", np.full(12, z[-1]), {"units": "m"}),
            "event_id": ("scenario", [f"event_{n}" for n in range(12)]),
        },
        coords={
            "scenario": names,
            "altitude_m": ("altitude_m", z, {"units": "m"}),
            "wavelength_nm": ("wavelength_nm", waves, {"units": "nm"}),
        },
        attrs={
            "catalogue_version": "v1",
            "schema_version": 1,
            "reference_wavelength_nm": 756.0,
            "mode_width": 1.6,
        },
    )


@pytest.fixture()
def source(catalogue, tmp_path):
    path = tmp_path / "custom.nc"
    catalogue.to_netcdf(path)
    return path


def test_raw_access_and_exact_selection(source, catalogue):
    data = aerosol.load_dataset(path=source)
    assert_array_equal(data.raw_extinction_per_m, catalogue.raw_extinction_per_m)
    assert (
        data.attrs["source_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    )
    assert data.attrs["source_mode"] == "local"
    listing = aerosol.scenarios(path=source)
    assert listing.sizes == {"scenario": 12}
    assert_array_equal(listing.event_id, catalogue.event_id)
    with pytest.raises(ValueError, match="exact aerosol scenario"):
        aerosol.profile("tropical", path=source)
    with pytest.raises(ValueError, match="either path or db_root"):
        aerosol.load_dataset(path=source, db_root=source.parent)
    with pytest.raises(ValueError, match="Unsupported"):
        aerosol.load_dataset("v2", path=source)


def test_smoothing_preserves_native_aod_and_raw_pairs(source, catalogue):
    p = aerosol.profile("tropical_extreme", path=source)
    raw = catalogue.sel(scenario="tropical_extreme")
    z = raw.altitude_m.values
    extinction = raw.raw_extinction_per_m.sel(wavelength_nm=756).values
    assert_allclose(p.core_aod, trapezoid(extinction, z), rtol=1e-14)
    assert_allclose(
        trapezoid(p.smoothed_observed_extinction_per_m, z), p.core_aod, rtol=1e-14
    )
    assert_array_equal(p.raw_extinction_per_m, raw.raw_extinction_per_m)
    assert_array_equal(
        p.raw_median_radius_uncertainty_nm, raw.raw_median_radius_uncertainty_nm
    )
    assert_array_equal(p.median_radius_nm.sel(altitude_m=z), raw.raw_median_radius_nm)
    assert not np.allclose(
        p.smoothed_observed_extinction_per_m, extinction, rtol=1e-4, atol=0
    )
    unchanged = aerosol.profile(
        "tropical_extreme", path=source, smoothing_fwhm_m=0, altitudes_m=z
    )
    assert_array_equal(unchanged.extinction_per_m, extinction)


def test_extensions_continuity_scale_heights_and_integrals(source):
    z = np.unique(np.r_[np.arange(-1000.0, 100001.0, 10), 17000 - 1e-4, 32000 + 1e-4])
    p = aerosol.profile(
        "tropical_extreme", path=source, altitudes_m=z, ground_altitude_m=1000
    )
    e = p.extinction_per_m
    assert_array_equal(e.where(p.altitude_m <= 1000, drop=True), 0)
    assert (e.where(p.altitude_m > 1000, drop=True) > 0).all()
    assert_allclose(e.sel(altitude_m=17000 - 1e-4), e.sel(altitude_m=17000), rtol=1e-7)
    assert_allclose(e.sel(altitude_m=32000 + 1e-4), e.sel(altitude_m=32000), rtol=1e-7)
    assert_allclose(e.sel(altitude_m=36000) / e.sel(altitude_m=32000), np.exp(-1))
    assert_allclose(
        trapezoid(e.sel(altitude_m=slice(1000, 17000)), z[(z >= 1000) & (z <= 17000)]),
        p.lower_extension_aod,
        rtol=3e-6,
    )
    assert_allclose(
        trapezoid(e.sel(altitude_m=slice(32000, None)), z[z >= 32000]),
        p.upper_extension_aod_to_infinity,
        rtol=3e-6,
    )
    assert_allclose(
        p.total_aod_to_infinity,
        p.core_aod + p.lower_extension_aod + p.upper_extension_aod_to_infinity,
    )
    assert_array_equal(
        p.median_radius_nm.sel(altitude_m=[-1000, 0, 17000]), [80, 80, 80]
    )
    assert_array_equal(p.median_radius_nm.sel(altitude_m=[32000, 100000]), [110, 110])
    assert_array_equal(
        p.region.sel(altitude_m=[0, 17000, 32000, 100000]), [-1, 0, 0, 1]
    )


def test_shared_tail_and_override_independent_of_smoothing(source):
    for tier in ("low", "typical", "elevated", "extreme"):
        for width in (0, 1500, 3000):
            p = aerosol.profile(f"tropical_{tier}", path=source, smoothing_fwhm_m=width)
            assert p.attrs["upper_scale_height_m"] == 4000
            assert_allclose(
                p.extinction_per_m.sel(altitude_m=40000)
                / p.extinction_per_m.sel(altitude_m=36000),
                np.exp(-1),
            )
    a = aerosol.profile("tropical_extreme", path=source)
    b = aerosol.profile("tropical_extreme", path=source, upper_scale_height_m=3200)
    assert b.attrs["upper_scale_height_source"] == "override"
    assert_array_equal(
        a.smoothed_observed_extinction_per_m, b.smoothed_observed_extinction_per_m
    )
    assert_allclose(
        b.upper_extension_aod_to_infinity / a.upper_extension_aod_to_infinity, 0.8
    )
    assert_allclose(
        b.extinction_per_m.sel(altitude_m=40000)
        / b.extinction_per_m.sel(altitude_m=32000),
        np.exp(-8000 / 3200),
    )


def test_zero_extensions_and_extreme_lower_scale_heights(source):
    p = aerosol.profile(
        "sh_midlat_low", path=source, lower_extension="zero", upper_extension="zero"
    )
    assert_array_equal(p.extinction_per_m.where(p.region != 0, drop=True), 0)
    assert p.lower_extension_aod == p.upper_extension_aod_to_infinity == 0
    assert (p.extinction_per_m.where(p.region == 0, drop=True) > 0).all()
    for scale in (0.01, 1e12):
        p = aerosol.profile("sh_midlat_low", path=source, lower_scale_height_m=scale)
        assert np.isfinite(p.extinction_per_m).all()
        assert float(p.lower_extension_aod) >= 0
    assert_allclose(
        p.lower_extension_aod,
        p.extinction_per_m.sel(altitude_m=17000) * 17000 / 2,
        rtol=1e-8,
    )


@pytest.mark.parametrize(
    "option",
    [
        {"smoothing_fwhm_m": -1},
        {"smoothing_fwhm_m": np.inf},
        {"upper_scale_height_m": "auto"},
        {"upper_scale_height_m": 0},
        {"lower_scale_height_m": -2000},
        {"lower_scale_height_m": np.nan},
        {"ground_altitude_m": 17000},
        {"ground_altitude_m": np.nan},
        {"lower_extension": "extend"},
        {"upper_extension": "reference"},
        {"altitudes_m": [1000]},
        {"altitudes_m": [1000, 0]},
        {"altitudes_m": [0, 0]},
        {"altitudes_m": [0, np.nan]},
    ],
)
def test_invalid_options(source, option):
    with pytest.raises(ValueError, match=r"must|Invalid"):
        aerosol.profile("sh_midlat_low", path=source, **option)


@pytest.mark.parametrize(
    ("change", "match"),
    [
        ("gap", "contiguous"),
        ("extinction", "positive"),
        ("radius", "radii"),
        ("flags", "cloud/invalid"),
        ("units", "units"),
        ("bounds", "bounds"),
        ("grid", "regular native"),
        ("schema", "schema"),
    ],
)
def test_invalid_source_not_silently_repaired(catalogue, tmp_path, change, match):
    if change == "gap":
        catalogue.observed_valid[0, 10] = 0
    elif change == "extinction":
        catalogue.raw_extinction_per_m[0, 10, :] = np.nan
    elif change == "radius":
        catalogue.raw_median_radius_nm[0, 10] = 590
    elif change == "flags":
        catalogue.source_aerosol_flag[0, 10, :] = 4
    elif change == "units":
        catalogue.raw_extinction_per_m.attrs["units"] = "km-1"
    elif change == "bounds":
        catalogue.observed_top_m[0] = 30000
    elif change == "grid":
        heights = catalogue.altitude_m.values.copy()
        heights[10] += 10
        catalogue = catalogue.assign_coords(
            altitude_m=("altitude_m", heights, {"units": "m"})
        )
    else:
        catalogue.attrs["schema_version"] = 2
    path = tmp_path / "invalid.nc"
    catalogue.to_netcdf(path)
    with pytest.raises(ValueError, match=match):
        aerosol.profile("sh_midlat_low", path=path)


class _RadiusDependentOptics(OpticalProperty):
    def cross_sections(self, wavelengths_nm, altitudes_m, **kwargs):
        xs = np.broadcast_to(
            np.asarray(kwargs["median_radius"])[:, None] * 1e-14,
            (len(altitudes_m), len(wavelengths_nm)),
        ).copy()
        return OpticalQuantities(extinction=xs, ssa=np.ones_like(xs))

    def atmosphere_quantities(self, atmo, **kwargs):
        q = self.cross_sections(
            atmo.wavelengths_nm, atmo.model_geometry.altitudes(), **kwargs
        )
        q.ssa *= q.extinction
        q.leg_coeff = np.zeros_like(atmo.storage.leg_coeff)
        q.leg_coeff[0] = 1
        return q


def test_constituent_normalization_in_atmosphere(source):
    z = np.arange(0.0, 50001.0, 1000.0)
    p = aerosol.profile("tropical_typical", path=source, altitudes_m=z)
    optical = _RadiusDependentOptics()
    c = aerosol.constituent("tropical_typical", optical, path=source, altitudes_m=z)
    assert isinstance(c, sk.constituent.ExtinctionScatterer)
    assert_array_equal(c.extinction_per_m, p.extinction_per_m)
    config = sk.Config()
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6372000.0,
        z,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=np.array([756.0]))
    c.add_to_atmosphere(atmosphere)
    assert_allclose(
        atmosphere.storage.total_extinction[:, 0],
        p.extinction_per_m,
        rtol=1e-12,
        atol=1e-25,
    )
    assert_allclose(atmosphere.storage.ssa[:, 0], p.extinction_per_m)


def test_default_mie_constituent(source, monkeypatch):
    # Exercise real Mie integration without downloading the refractive index.
    monkeypatch.setattr(
        sk.mie.refractive,
        "H2SO4",
        lambda: sk.mie.refractive.RefractiveIndex(
            lambda wavelength: np.full_like(wavelength, 1.45, dtype=complex),
            "test_sulfate",
        ),
    )
    z = np.array([0.0, 17000.0, 22000.0, 32000.0, 50000.0])
    c = aerosol.constituent("tropical_typical", path=source, altitudes_m=z)
    xs = c._optical_property.cross_sections(
        np.array([525.0, 756.0, 1021.0]), altitudes_m=z, median_radius=c.median_radius
    )
    assert np.isfinite(xs.extinction).all()
    assert (xs.extinction > 0).all()
    assert_allclose(
        c.number_density * xs.extinction[:, 1], c.extinction_per_m, rtol=1e-9
    )
    config = sk.Config()
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6372000.0,
        z,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(
        geometry, config, wavelengths_nm=np.array([525.0, 756.0, 1021.0])
    )
    atmosphere["aerosol"] = c
    atmosphere.internal_object()
    assert_allclose(
        atmosphere.storage.total_extinction[:, 1], c.extinction_per_m, rtol=1e-9
    )
    # Nonabsorbing sulfate must retain unit albedo after atmosphere assembly.
    assert_allclose(atmosphere.storage.ssa[c.extinction_per_m > 0], 1.0, atol=1e-12)


def test_packaged_catalogue_and_cache(tmp_path):
    db = database.StratosphericAerosolDatabase(db_root=tmp_path)
    path = db.path()
    data = aerosol.load_dataset(db_root=tmp_path)
    assert data.sizes["scenario"] == 12
    assert float(data.raw_extinction_per_m.wavelength_nm.sel(wavelength_nm=756)) == 756
    for name in data.scenario.values:
        p = aerosol.profile(str(name), db_root=tmp_path)
        assert np.isfinite(p.extinction_per_m).all()
        assert (p.extinction_per_m >= 0).all()
        assert float(p.extinction_per_m.sel(altitude_m=0)) == 0
        assert float(p.observed_bottom_m) < 18000 or float(p.observed_top_m) > 30000
        assert p.attrs["upper_scale_height_m"] > 0
    for band in ("sh_midlat", "tropical", "nh_midlat"):
        group = data.sel(
            scenario=[
                f"{band}_{tier}" for tier in ("low", "typical", "elevated", "extreme")
            ]
        )
        assert len(np.unique(group.reference_upper_scale_height_m)) == 1
        assert (np.diff(group.selection_aod_18_30) > 0).all()
    path.write_bytes(b"damaged")
    with pytest.raises(OSError, match="checksum mismatch"):
        db.path()
    db.clear()
    assert not path.exists()
    assert db.path().exists()
    with pytest.raises(ValueError, match="contains only"):
        db.path("../other.nc")
    with pytest.raises(ValueError, match="Unsupported"):
        database.StratosphericAerosolDatabase("v2", db_root=tmp_path)
