from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk
import xarray as xr
from scipy.constants import c, h
from scipy.integrate import trapezoid


def _atmosphere(wavelengths=None, wavenumbers=None, num_stokes=1):
    config = sk.Config()
    config.num_stokes = num_stokes
    config.emission_source = sk.EmissionSource.VolumeEmissionRate
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    geometry = sk.Geometry1D(
        cos_sza=-0.6,
        solar_azimuth=0,
        earth_radius_m=6372000,
        altitude_grid_m=np.linspace(0, 100000, 11),
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    atmo = sk.Atmosphere(
        geometry,
        config,
        wavelengths_nm=None if wavelengths is None else np.array(wavelengths),
        wavenumber_cminv=None if wavenumbers is None else np.array(wavenumbers),
    )
    atmo.temperature_k = np.full(11, 200.0)
    atmo.pressure_pa = np.ones(11)
    return atmo, config, geometry


def _source(**kwargs):
    return sk.constituent.SpectralVolumeEmissionRate(
        np.array([0, 50000, 100000]),
        np.array([0.0, 8.0, 2.0]),
        np.array([500.0, 550.0, 600.0]),
        np.array([0.0, 1.0, 0.0]),
        **kwargs,
    )


@pytest.mark.parametrize("wavenumber_space", [False, True])
def test_template_conserves_band_photons_in_both_spectral_coordinates(wavenumber_space):
    if wavenumber_space:
        coordinate = np.linspace(1e7 / 600, 1e7 / 500, 10001)
        atmosphere, _, _ = _atmosphere(wavenumbers=coordinate)
        assert atmosphere.spectral_coordinate == "wavenumber_cminv"
    else:
        coordinate = np.linspace(500, 600, 1001)
        atmosphere, _, _ = _atmosphere(wavelengths=coordinate)
        assert atmosphere.spectral_coordinate == "wavelength_nm"
    _source().add_to_atmosphere(atmosphere)
    actual = (
        trapezoid(atmosphere.storage.emission_source, coordinate, axis=1) * 4 * np.pi
    )
    np.testing.assert_allclose(
        actual,
        np.interp(atmosphere._native_altitudes(), [0, 50000, 100000], [0, 8, 2]),
        rtol=1e-7,
    )


def test_fitting_subset_does_not_renormalize_ver_and_outside_support_is_zero():
    broad, _, _ = _atmosphere(wavelengths=[450.0, 500.0, 525.0, 550.0, 600.0, 650.0])
    subset, _, _ = _atmosphere(wavelengths=[525.0, 550.0])
    constituent = _source()
    constituent.add_to_atmosphere(broad)
    constituent.add_to_atmosphere(subset)
    np.testing.assert_allclose(
        subset.storage.emission_source, broad.storage.emission_source[:, 2:4]
    )
    np.testing.assert_array_equal(broad.storage.emission_source[:, [0, -1]], 0)
    np.testing.assert_allclose(subset.storage.emission_source[5, 1] * 4 * np.pi, 8 / 50)


def test_energy_source_converts_each_photon_with_its_wavelength():
    photons, _, _ = _atmosphere(wavelengths=[525.0, 575.0])
    energy, _, _ = _atmosphere(wavelengths=[525.0, 575.0])
    _source().add_to_atmosphere(photons)
    _source(emission_units="energy").add_to_atmosphere(energy)
    np.testing.assert_allclose(
        energy.storage.emission_source,
        photons.storage.emission_source * h * c / (photons.wavelengths_nm * 1e-9),
        atol=0,
    )


@pytest.mark.parametrize("num_stokes", [1, 3])
def test_free_ver_radiance_is_linear_and_weighting_function_matches_finite_difference(
    num_stokes,
):
    atmosphere, config, geometry = _atmosphere(
        wavelengths=[525.0, 550.0, 575.0], num_stokes=num_stokes
    )
    constituent = _source()
    atmosphere["band"] = constituent
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.TangentAltitudeSolar(30000, 0, 200000, -0.6))
    engine = sk.Engine(config, geometry, viewing)
    baseline = engine.calculate_radiance(atmosphere)
    original = constituent.photon_ver
    for index in range(3):
        perturbed = original.copy()
        perturbed[index] += 0.01
        constituent.photon_ver = perturbed
        plus = engine.calculate_radiance(atmosphere).radiance.to_numpy()
        numerical = (plus - baseline.radiance.to_numpy()) / 0.01
        analytical = baseline.wf_band_photon_ver.isel(band_altitude=index).to_numpy()
        np.testing.assert_allclose(analytical, numerical, rtol=1e-9, atol=1e-8)
    constituent.photon_ver = original * 2
    doubled = engine.calculate_radiance(atmosphere).radiance.to_numpy()
    np.testing.assert_allclose(doubled, 2 * baseline.radiance.to_numpy())
    assert np.all(doubled[..., 0] > 0)
    if num_stokes == 3:
        np.testing.assert_array_equal(doubled[..., 1:], 0)
    constituent.photon_ver = np.zeros_like(original)
    zero = engine.calculate_radiance(atmosphere)
    np.testing.assert_array_equal(zero.radiance, 0)
    np.testing.assert_allclose(zero.wf_band_photon_ver, baseline.wf_band_photon_ver)


@pytest.mark.parametrize("species", ["FeO", "NiO"])
def test_local_template_loading_preserves_provisional_metadata(tmp_path, species):
    folder = tmp_path / "spectroscopy" / "metals" / "emission"
    folder.mkdir(parents=True)
    dataset = xr.Dataset(
        {"photon_spectrum": ("wavelength_nm", [0, 1, 0])},
        coords={"wavelength_nm": [500, 550, 600]},
        attrs={"wavelength_medium": "unspecified", "attribution": "provisional"},
    )
    dataset.to_netcdf(folder / f"{species}.nc")
    database = sk.database.MetalSpectroscopyDatabase(db_root=tmp_path)
    factory = getattr(sk.constituent, f"{species}VolumeEmissionRate")
    constituent = factory([0, 100000], [1, 1], db=database)
    assert constituent.metadata == dataset.attrs
    np.testing.assert_allclose(
        trapezoid(constituent.photon_spectrum, [500, 550, 600]), 1
    )


@pytest.mark.parametrize(
    ("wavelengths", "shape"),
    [
        ([500, 500], [1, 1]),
        ([500, 600], [0, 0]),
        ([500, 600], [-1, 1]),
        ([500, 600], [1, np.nan]),
        ([500], [1]),
        ([500, 600], [1]),
    ],
)
def test_invalid_templates_rejected(wavelengths, shape):
    with pytest.raises(ValueError, match=r"wavelengths_nm|photon_spectrum"):
        sk.constituent.SpectralVolumeEmissionRate([0, 1], [1, 1], wavelengths, shape)


def test_integrated_spectral_grid_rejected_instead_of_misnormalizing():
    atmosphere, _, _ = _atmosphere(wavelengths=[500.0, 550.0, 600.0])
    atmosphere._spectral_integration_mode = (
        sk.SpectralGridMode.AtmosphereIntegratedLineShape
    )
    with pytest.raises(ValueError, match="monochromatic"):
        _source().add_to_atmosphere(atmosphere)
