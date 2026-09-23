from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
import sasktran2 as sk
import sasktran2.atmosphere as atmosphere_module
import xarray as xr
from sasktran2._core_rust import LineDatabaseType, PyLineAbsorber
from sasktran2.optical.hitran import LineAbsorber


@pytest.fixture()
def line_db(tmp_path):
    """Small local O2 database with distinct rotational energies in each band."""
    records = []
    for band, center, energy, j in [
        ("0-0", 13100.0, 20.0, 1),
        ("0-0", 13100.5, 400.0, 5),
        ("1-1", 13100.1, 1600.0, 1),
        ("1-1", 13100.6, 2100.0, 5),
        ("1-0", 14500.0, 30.0, 1),
        ("1-0", 14500.4, 380.0, 5),
    ]:
        upper, lower = band.split("-")
        record = (
            f"{7:2d}{1:1d}{center:12.6f}{1e-27:10.3e}{0.1:10.3e}"
            f"{0.06:5.3f}{0.10:5.3f}{energy:10.4f}{0.7:4.2f}{0.003:8.6f}"
            f"{'b ' + upper:>15}{'X ' + lower:>15}{j!s:>15}{'':15}"
            + " " * 19
            + f"{float(2 * j + 1):7.1f}{1.0:7.1f}"
        )
        assert len(record) == 160
        records.append(record)
    (tmp_path / "O2.data").write_text("\n".join(records) + "\n")
    return SimpleNamespace(path=lambda _species: tmp_path)


def scenario(
    *,
    temperature_derivative=True,
    calculate_derivatives=True,
    coordinate="wavenumber",
    spectral_mode=sk.SpectralGridMode.Monochromatic,
    single_scatter=False,
    grid=None,
):
    config = sk.Config()
    config.emission_source = sk.EmissionSource.VolumeEmissionRate
    config.single_scatter_source = (
        sk.SingleScatterSource.Exact
        if single_scatter
        else sk.SingleScatterSource.NoSource
    )
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.spectral_grid_mode = spectral_mode
    altitudes = np.array([0.0, 20_000.0, 40_000.0, 60_000.0])
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6_372_000,
        altitudes,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    if grid is None:
        grid = np.linspace(13099.8, 13100.8, 501)
        if coordinate == "wavelength":
            grid = np.sort(1e7 / grid)
    spectral = {
        "wavenumber_cminv" if coordinate == "wavenumber" else "wavelengths_nm": grid
    }
    atmo = sk.Atmosphere(
        geometry,
        config,
        **spectral,
        calculate_derivatives=calculate_derivatives,
        temperature_derivative=temperature_derivative,
        pressure_derivative=False,
        specific_humidity_derivative=False,
    )
    atmo.temperature_k = np.array([280.0, 220.0, 250.0, 210.0])
    atmo.pressure_pa = np.array([1e4, 2e3, 200.0, 20.0])
    viewing = sk.ViewingGeometry()
    for tangent in [10_000.0, 30_000.0]:
        viewing.add_ray(sk.TangentAltitudeSolar(tangent, 0.0, 200_000.0, 0.6))
    return atmo, sk.Engine(config, geometry, viewing)


def constituent(db, band="0-0", model="einstein_a_branching", **kwargs):
    # Deliberately use a different grid from atmospheric T, and include a zero VER.
    return sk.constituent.O2BandEmissionRate(
        np.array([0.0, 30_000.0, 60_000.0]),
        np.array([0.0, 3.0, 1.0]),
        band=band,
        db=db,
        line_weight_model=model,
        **kwargs,
    )


def add_self_absorption(atmo, db):
    absorber = LineAbsorber.__new__(LineAbsorber)
    absorber._internal = PyLineAbsorber(
        LineDatabaseType.HITRAN,
        "O2",
        db.path("O2").as_posix(),
        py_tips=lambda _m, _i, t: t**1.5,
        py_molmass=lambda _m, _i: 31.9988,
    )
    atmo["o2"] = sk.constituent.VMRAltitudeAbsorber(
        absorber, atmo.model_geometry.altitudes(), np.full(4, 0.21)
    )


@pytest.mark.parametrize("model", ["einstein_a_branching", "hitran_line_strength"])
@pytest.mark.parametrize("band", ["0-0", "1-1", "1-0"])
def test_band_source_conservation_and_temperature_derivative(line_db, model, band):
    grid = np.linspace(14499.8, 14500.8, 501) if band == "1-0" else None
    atmo, _ = scenario(grid=grid)
    emission = constituent(line_db, band, model)
    emission.add_to_atmosphere(atmo)
    emission.register_derivative(atmo, "band")
    ver_mapping = atmo.storage.get_derivative_mapping("wf_band_photon_ver")
    t_mapping = atmo.storage.get_derivative_mapping("wf_band_temperature_k")
    expected_ver = np.interp(
        atmo.model_geometry.altitudes(), emission.altitudes_m, emission.photon_ver
    )
    np.testing.assert_allclose(
        np.trapezoid(atmo.storage.emission_source, atmo.wavenumbers_cminv, axis=1)
        * 4
        * np.pi,
        expected_ver,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        np.trapezoid(ver_mapping.d_emission, atmo.wavenumbers_cminv, axis=1)
        * 4
        * np.pi,
        1.0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        np.trapezoid(t_mapping.d_emission, atmo.wavenumbers_cminv, axis=1),
        0.0,
        atol=1e-12,
    )
    analytic = t_mapping.d_emission.copy()
    step = 0.001
    for sign in [1, -1]:
        atmo.temperature_k += sign * step
        atmo.storage.emission_source[:] = 0.0
        emission.add_to_atmosphere(atmo)
        if sign == 1:
            above = atmo.storage.emission_source.copy()
        else:
            below = atmo.storage.emission_source.copy()
        atmo.temperature_k -= sign * step
    np.testing.assert_allclose(
        analytic, (above - below) / (2 * step), rtol=2e-6, atol=1e-12
    )
    # Independent normalized weights must follow the temperature, not constructor state.
    weights = emission.line_weights([180.0, 280.0])
    np.testing.assert_allclose(weights.sum(axis=1), 1.0)
    assert not np.allclose(weights[0], weights[1])


@pytest.mark.parametrize("coordinate", ["wavenumber", "wavelength"])
@pytest.mark.parametrize("self_absorption", [False, True])
def test_radiance_ver_and_temperature_jacobians(line_db, coordinate, self_absorption):
    atmo, engine = scenario(coordinate=coordinate)
    bands = {
        name: constituent(line_db, band)
        for name, band in [("b00", "0-0"), ("b11", "1-1")]
    }
    for name, emission in bands.items():
        atmo[name] = emission
    if self_absorption:
        add_self_absorption(atmo, line_db)
    result = engine.calculate_radiance(atmo)
    assert np.max(result.radiance.to_numpy()) > 0.0
    for name, emission in bands.items():
        for level in range(len(emission.photon_ver)):
            original = emission.photon_ver[level]
            # Radiance is exactly linear in VER at fixed optical properties;
            # a unit step avoids subtractive cancellation under the other band.
            step = 1.0
            emission.photon_ver[level] = original + step
            above = engine.calculate_radiance(atmo).radiance
            emission.photon_ver[level] = original
            baseline = engine.calculate_radiance(atmo).radiance
            numeric = (above - baseline) / step
            analytic = result[f"wf_{name}_photon_ver"].isel({f"{name}_altitude": level})
            np.testing.assert_allclose(analytic, numeric, rtol=1e-6, atol=1e-6)
    for level in range(len(atmo.temperature_k)):
        original = atmo.temperature_k[level]
        step = 0.001
        atmo.temperature_k[level] = original + step
        above = engine.calculate_radiance(atmo).radiance
        atmo.temperature_k[level] = original - step
        below = engine.calculate_radiance(atmo).radiance
        atmo.temperature_k[level] = original
        np.testing.assert_allclose(
            result.wf_temperature_k.isel(altitude=level),
            (above - below) / (2 * step),
            rtol=2e-5,
            atol=1e-5,
        )


@pytest.mark.parametrize("model", ["einstein_a_branching", "hitran_line_strength"])
def test_joint_ver_temperature_linearization_products(line_db, model):
    atmo, engine = scenario(single_scatter=True)
    bands = {
        "b00": constituent(line_db, model=model),
        "b11": sk.constituent.O2BandEmissionRate(
            [0.0, 60_000.0],
            [0.5, 1.5],
            band="1-1",
            line_weight_model=model,
            db=line_db,
        ),
    }
    for name, emission in bands.items():
        atmo[name] = emission
    add_self_absorption(atmo, line_db)
    atmo["rayleigh"] = sk.constituent.Rayleigh()

    # Exercise independent VER grids and the shared atmospheric T parameter,
    # including absorption and solar scattering in the same linearization.
    directions = {
        "b00_photon_ver": [0.0, 0.3, -0.1],
        "b11_photon_ver": [0.2, -0.3],
        "temperature_k": [-0.8, 0.4, 0.3, -0.5],
    }
    linearization = engine.linearize(atmo)
    tangent = linearization.tangent_template[list(directions)]
    for name, values in directions.items():
        tangent[name].data[:] = values
    analytic = linearization.jvp(tangent)
    cotangent = xr.DataArray(
        np.linspace(-0.7, 1.0, analytic.size).reshape(analytic.shape),
        dims=analytic.dims,
        coords=analytic.coords,
    )
    gradient = linearization.vjp(cotangent, parameters=tuple(directions))
    np.testing.assert_allclose(
        float((analytic * cotangent).sum()),
        sum(float((gradient[name] * tangent[name]).sum()) for name in directions),
        rtol=1e-11,
    )

    temperature = atmo.temperature_k.copy()
    profiles = {name: emission.photon_ver.copy() for name, emission in bands.items()}
    step = 0.001
    values = []
    for sign in [1, -1]:
        atmo.temperature_k = (
            temperature + sign * step * tangent.temperature_k.to_numpy()
        )
        for name, emission in bands.items():
            emission.photon_ver = (
                profiles[name] + sign * step * tangent[f"{name}_photon_ver"].to_numpy()
            )
        values.append(engine.calculate_radiance(atmo).radiance)
    numeric = (values[0] - values[1]) / (2 * step)
    np.testing.assert_allclose(analytic, numeric, rtol=2e-5, atol=1e-5)


def population(db):
    # This temperature intentionally differs from atmospheric T.
    state = xr.Dataset(
        {
            "temperature": ("altitude", [300.0, 310.0, 320.0]),
            "O2(b)": ("altitude", [1.0, 10.0, 5.0]),
            "O2(b, v=1)": ("altitude", [3.0, 4.0, 2.0]),
        },
        coords={"altitude": [0.0, 30_000.0, 60_000.0]},
    )
    return sk.constituent.PopulationEmissionRate(state, db=db)


@pytest.mark.parametrize(
    ("name", "index"),
    [
        ("photon_ver", None),
        ("altitudes_m", None),
        ("wavelengths_nm", None),
        ("weights", None),
        ("line_list_photon_ver", 0),
        ("line_list_photon_ver", 1),
        ("line_list_wavelengths_nm", 0),
        ("line_list_wavelengths_nm", 1),
        ("line_list_weights", 0),
        ("line_list_weights", 1),
    ],
)
def test_population_inspection_arrays_reject_mutation(line_db, name, index):
    pop = population(line_db)
    inspection = getattr(pop, name)
    if index is not None:
        inspection = inspection(index)
    expected = inspection.copy()

    # These construction-time arrays no longer control the emitted source.
    # Legacy retrieval updates must fail instead of silently editing a snapshot.
    with pytest.raises(ValueError, match="read-only"):
        inspection[...] = 0.0
    with pytest.raises(ValueError, match="WRITEABLE"):
        inspection.setflags(write=True)

    # A retained inspection view must keep its native owner alive.
    del pop
    np.testing.assert_array_equal(inspection, expected)


def test_population_conversion_uses_same_band_path(line_db):
    pop = population(line_db)
    bands = pop.to_band_emissions()
    assert set(bands) == {"0-0", "1-1", "1-0"}
    np.testing.assert_allclose(bands["0-0"].photon_ver, np.array([1, 10, 5]) * 0.0758)
    np.testing.assert_allclose(bands["1-1"].photon_ver, np.array([3, 4, 2]) * 0.07)
    atmo, engine = scenario()
    atmo["population"] = pop
    combined = engine.calculate_radiance(atmo)
    direct_atmo, direct_engine = scenario()
    for name, band in bands.items():
        direct_atmo[name] = band
    direct = direct_engine.calculate_radiance(direct_atmo)
    np.testing.assert_allclose(combined.radiance, direct.radiance, rtol=1e-14)
    np.testing.assert_allclose(
        combined.wf_temperature_k, direct.wf_temperature_k, rtol=1e-12, atol=1e-12
    )
    original_ver = bands["0-0"].photon_ver.copy()
    bands["0-0"].photon_ver = [0, 0, 0]
    np.testing.assert_array_equal(
        engine.calculate_radiance(atmo).radiance, combined.radiance
    )
    assert not np.allclose(
        direct_engine.calculate_radiance(direct_atmo).radiance, direct.radiance
    )
    bands["0-0"].photon_ver[:] = original_ver
    np.testing.assert_array_equal(
        direct_engine.calculate_radiance(direct_atmo).radiance, direct.radiance
    )


@pytest.mark.parametrize(
    ("calculate_derivatives", "temperature_derivative"), [(False, True), (True, False)]
)
def test_derivative_flags_preserve_values(
    line_db, calculate_derivatives, temperature_derivative
):
    atmo, engine = scenario(
        calculate_derivatives=calculate_derivatives,
        temperature_derivative=temperature_derivative,
    )
    atmo["band"] = constituent(line_db)
    result = engine.calculate_radiance(atmo)
    assert "wf_temperature_k" not in result
    assert ("wf_band_photon_ver" in result) == calculate_derivatives
    reference_atmo, reference_engine = scenario()
    reference_atmo["band"] = constituent(line_db)
    np.testing.assert_array_equal(
        result.radiance, reference_engine.calculate_radiance(reference_atmo).radiance
    )


def test_band_spectra_and_derivatives_reduce_on_fine_grid(line_db, monkeypatch):
    fine_atmo, _ = scenario()
    fine_grid = sk.basis.Grid.from_triangles(fine_atmo.wavenumbers_cminv)
    coarse_grid = sk.basis.Grid.from_triangles(np.linspace(13099.8, 13100.8, 21))
    native_atmosphere = atmosphere_module.PyAtmosphere

    def with_fine_grid(*args):
        return native_atmosphere(*args[:7], fine_grid._internal_object(), *args[8:])

    with monkeypatch.context() as patch:
        patch.setattr(atmosphere_module, "PyAtmosphere", with_fine_grid)
        coarse_atmo, _ = scenario(
            grid=np.linspace(13099.8, 13100.8, 21),
            spectral_mode=sk.SpectralGridMode.AtmosphereIntegratedLineShape,
        )
    emission = constituent(line_db)
    for atmo in [fine_atmo, coarse_atmo]:
        emission.add_to_atmosphere(atmo)
        emission.register_derivative(atmo, "band")
    mapping = fine_grid.mapping_to(coarse_grid, normalize=False)
    np.testing.assert_allclose(
        coarse_atmo.storage.emission_source,
        fine_atmo.storage.emission_source @ mapping.T,
        rtol=2e-12,
        atol=1e-14,
    )
    for name in ["wf_band_photon_ver", "wf_band_temperature_k"]:
        fine = fine_atmo.storage.get_derivative_mapping(name).d_emission
        coarse = coarse_atmo.storage.get_derivative_mapping(name).d_emission
        np.testing.assert_allclose(coarse, fine @ mapping.T, rtol=2e-12, atol=1e-14)


@pytest.mark.parametrize("bad_temperature", [0.0, -1.0, np.nan])
def test_invalid_temperatures_raise(line_db, bad_temperature):
    atmo, _ = scenario()
    atmo.temperature_k[1] = bad_temperature
    with pytest.raises(ValueError, match="positive and finite"):
        constituent(line_db).add_to_atmosphere(atmo)


def test_ver_update_and_band_validation(line_db):
    emission = constituent(line_db)
    profile_view = emission.photon_ver
    emission.photon_ver = [1.0, 2.0, 3.0]
    np.testing.assert_array_equal(emission.photon_ver, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(profile_view, [1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="length"):
        emission.photon_ver = [1.0]
    with pytest.raises(ValueError, match="non-negative"):
        emission.photon_ver = [-1.0, 2.0, 3.0]
    with pytest.raises(ValueError, match="Unsupported O2 band"):
        constituent(line_db, band="2-0")


@pytest.mark.parametrize(
    ("mode", "expected"),
    [("zero", [0.0, 1.5, 2.5, 0.0]), ("extend", [1.0, 1.5, 2.5, 3.0])],
)
def test_ver_altitude_bounds(line_db, mode, expected):
    emission = sk.constituent.O2BandEmissionRate(
        [10_000, 50_000], [1, 3], db=line_db, out_of_bounds_mode=mode
    )
    atmo, _ = scenario()
    emission.add_to_atmosphere(atmo)
    emission.register_derivative(atmo, "band")
    np.testing.assert_allclose(
        np.trapezoid(atmo.storage.emission_source, atmo.wavenumbers_cminv, axis=1)
        * 4
        * np.pi,
        expected,
        atol=1e-12,
    )
    # The same interpolation must control both the value and the VER Jacobian.
    mapping = atmo.storage.get_derivative_mapping("wf_band_photon_ver")
    np.testing.assert_allclose(mapping.interpolator @ emission.photon_ver, expected)


@pytest.mark.parametrize("altitudes", [[], [1.0, 0.0], [0.0, 0.0], [0.0, np.nan]])
def test_invalid_altitude_grid(line_db, altitudes):
    with pytest.raises(ValueError, match="altitude"):
        sk.constituent.O2BandEmissionRate(
            altitudes, np.ones(len(altitudes)), db=line_db
        )
