from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def _case(*, thermal=False, spherical=False, solar=True, derivatives=True):
    config = sk.Config()
    config.num_streams = 2
    config.num_threads = 1
    config.wavelength_batch_size = 8
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.TwoStream
    config.emission_source = (
        sk.EmissionSource.TwoStream if thermal else sk.EmissionSource.VolumeEmissionRate
    )
    altitudes = np.arange(0.0, 40_001.0, 5_000.0)
    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6_371_000.0,
        altitudes,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical if spherical else sk.GeometryType.PlaneParallel,
    )
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.GroundViewingSolar(0.6, 0.3, 0.7, 200_000.0))
    viewing.add_ray(sk.GroundViewingSolar(0.6, -0.4, 0.35, 200_000.0))
    if spherical:
        viewing.add_ray(sk.TangentAltitudeSolar(12_000.0, 0.4, 200_000.0, 0.6))
    atmosphere = sk.Atmosphere(
        geometry, config, numwavel=9, calculate_derivatives=derivatives
    )
    atmosphere.storage.total_extinction[:] = 2.0e-5
    atmosphere.storage.ssa[:] = 0.8
    atmosphere.leg_coeff.a1[0] = 1.0
    atmosphere.leg_coeff.a1[1] = 0.6
    atmosphere.storage.emission_source[:] = 1.0e-5 * np.array(
        [0.0, 0.5, 1.0, 1.5, 0.0, 1.0, 0.5, 2.0, 1.0]
    )
    if thermal:
        atmosphere.storage.emission_source[:] /= atmosphere.storage.total_extinction * (
            1.0 - atmosphere.storage.ssa
        )
    atmosphere.surface.albedo[:] = 0.25
    atmosphere.surface.emission[:] = 0.1
    atmosphere.storage.solar_irradiance[:] = 1.1 if solar else 0.0
    return config, geometry, viewing, atmosphere


@pytest.mark.parametrize("spherical", [False, True])
def test_volume_emission_matches_equivalent_thermal_source(spherical):
    results = []
    for thermal in [False, True]:
        config, geometry, viewing, atmosphere = _case(
            thermal=thermal, spherical=spherical
        )
        results.append(
            sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)
        )
    np.testing.assert_allclose(
        results[0].radiance, results[1].radiance, rtol=1.0e-9, atol=1.0e-12
    )
    np.testing.assert_allclose(
        results[0].wf_emission.sum("altitude").isel(wavelength=[1, 2, 3, 5, 6, 7, 8])
        * 4.0e-6,
        results[1].wf_emission.sum("altitude").isel(wavelength=[1, 2, 3, 5, 6, 7, 8]),
        rtol=1.0e-8,
        atol=1.0e-12,
    )


@pytest.mark.parametrize("extinction", [0.0, 1.0e-12, 2.0e-5])
def test_volume_emission_absorbing_slab(extinction):
    config, geometry, viewing, atmosphere = _case(solar=False)
    atmosphere.storage.total_extinction[:] = extinction
    atmosphere.storage.ssa[:] = 0.0
    atmosphere.surface.albedo[:] = 0.0
    atmosphere.surface.emission[:] = 0.0
    result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)
    path_length = 40_000.0 / np.array([0.7, 0.35])
    optical_depth = extinction * path_length
    integral = path_length.copy()
    if extinction > 0.0:
        integral *= -np.expm1(-optical_depth) / optical_depth
    expected = atmosphere.storage.emission_source[0, :, None] * integral
    np.testing.assert_allclose(result.radiance[..., 0], expected, rtol=1.0e-10)
    np.testing.assert_allclose(
        result.wf_emission.sum("altitude")[..., 0],
        np.broadcast_to(integral, expected.shape),
        rtol=1.0e-10,
    )
    assert all(np.isfinite(value).all() for value in result.data_vars.values())


@pytest.mark.parametrize("vacuum_layers", [False, True])
def test_volume_emission_conservative_flux(vacuum_layers):
    config, geometry, _, atmosphere = _case(solar=False)
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.GroundViewingSolar(0.6, 0.0, 0.5, 200_000.0))
    atmosphere.storage.ssa[:] = 1.0
    if vacuum_layers:
        atmosphere.storage.total_extinction[4:] = 0.0
    atmosphere.surface.albedo[:] = 1.0
    atmosphere.surface.emission[:] = 0.0
    result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)
    # With a reflecting lower boundary and no absorption, all 4*pi*VER
    # escapes at the top. The two-stream upward flux is 2*pi*mu*I(mu).
    expected = 4.0 * 40_000.0 * atmosphere.storage.emission_source[0]
    np.testing.assert_allclose(result.radiance[:, 0, 0], expected, rtol=2.0e-8)
    assert all(np.isfinite(value).all() for value in result.data_vars.values())


def test_volume_emission_matches_discrete_ordinates():
    config, geometry, viewing, atmosphere = _case(solar=False)
    twostream = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)
    config, geometry, viewing, atmosphere = _case(thermal=True, solar=False)
    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.emission_source = sk.EmissionSource.DiscreteOrdinates
    discrete_ordinates = sk.Engine(config, geometry, viewing).calculate_radiance(
        atmosphere
    )
    np.testing.assert_allclose(
        twostream.radiance, discrete_ordinates.radiance, rtol=1.0e-9, atol=1.0e-12
    )


def test_volume_emission_vacuum_limb_matches_direct_source():
    results = []
    for scattering in [True, False]:
        config, geometry, viewing, atmosphere = _case(spherical=True, solar=False)
        if not scattering:
            config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
        atmosphere.storage.total_extinction[:] = 0.0
        atmosphere.storage.ssa[:] = 0.0
        atmosphere.surface.albedo[:] = 0.0
        atmosphere.surface.emission[:] = 0.0
        atmosphere.storage.emission_source[:] *= np.arange(9)[:, None]
        results.append(
            sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)
        )
    for name in ["radiance", "wf_emission"]:
        np.testing.assert_allclose(results[0][name], results[1][name], rtol=1.0e-12)


@pytest.mark.parametrize(
    ("spherical", "delta_m"), [(False, False), (True, False), (False, True)]
)
def test_volume_emission_jacobians(spherical, delta_m):
    config, geometry, viewing, atmosphere = _case(spherical=spherical)
    if delta_m:
        config.delta_m_scaling = True
        atmosphere.leg_coeff.a1[2] = 0.25
        atmosphere.leg_coeff.a1[3] = 0.1
    if spherical:
        config.los_refraction = True
        geometry.refractive_index = 1.0 + 3.0e-4 * np.exp(
            -geometry.altitudes() / 7_000.0
        )
    # A localized source includes zero-to-positive transitions and entire
    # zero-emission wavelengths. VER derivatives must remain linear there.
    atmosphere.storage.emission_source[:3] = 0.0
    atmosphere.storage.emission_source[6:] = 0.0
    atmosphere.storage.total_extinction[:] *= np.linspace(1.0, 0.2, 9)[:, None]
    mapping = atmosphere.surface.get_derivative_mapping("wf_surface_emission")
    mapping.d_emission[:] = 1.0
    engine = sk.Engine(config, geometry, viewing)

    def calculate():
        # The raw interface delta-scales these arrays in place. Restore the
        # unscaled inputs so each finite difference applies scaling once.
        fields = [
            atmosphere.storage.total_extinction,
            atmosphere.storage.ssa,
            atmosphere.storage.leg_coeff,
        ]
        originals = [field.copy() for field in fields]
        output = engine.calculate_radiance(atmosphere)
        for field, original in zip(fields, originals, strict=True):
            field[:] = original
        return output

    result = calculate()
    for name, values, step in [
        ("wf_emission", atmosphere.storage.emission_source, 1.0e-8),
        ("wf_extinction", atmosphere.storage.total_extinction, 1.0e-9),
        ("wf_ssa", atmosphere.storage.ssa, 1.0e-5),
        ("wf_leg_coeff_1", atmosphere.leg_coeff.a1[1], 1.0e-5),
        ("wf_leg_coeff_2", atmosphere.leg_coeff.a1[2], 1.0e-5),
    ]:
        for level in [0, 4, 8]:
            original = values[level].copy()
            values[level] = original + step
            plus = calculate().radiance
            # Emission is nonnegative; use a forward difference, which is
            # exact for this linear source, at its zero values.
            if name == "wf_emission":
                minus = result.radiance
                denominator = step
            else:
                values[level] = original - step
                minus = calculate().radiance
                denominator = 2.0 * step
            values[level] = original
            np.testing.assert_allclose(
                result[name].isel(altitude=level),
                (plus - minus) / denominator,
                rtol=2.0e-5,
                atol=2.0e-8,
                err_msg=f"{name}, level {level}",
            )
    for name, values in [
        ("wf_albedo", atmosphere.surface.albedo),
        ("wf_surface_emission", atmosphere.surface.emission),
    ]:
        original = values.copy()
        values[:] = original + 1.0e-5
        plus = calculate().radiance
        values[:] = original - 1.0e-5
        minus = calculate().radiance
        values[:] = original
        np.testing.assert_allclose(
            result[name], (plus - minus) / 2.0e-5, rtol=1.0e-6, atol=1.0e-10
        )


@pytest.mark.parametrize("spherical", [False, True])
@pytest.mark.parametrize("single_scatter", [False, True])
def test_volume_emission_superposition_and_batching(spherical, single_scatter):
    config, geometry, viewing, atmosphere = _case(spherical=spherical)
    if single_scatter:
        config.single_scatter_source = sk.SingleScatterSource.Exact
    atmosphere.storage.emission_source[::2] = 0.0
    engine = sk.Engine(config, geometry, viewing)
    combined = engine.calculate_radiance(atmosphere)
    atmosphere.storage.solar_irradiance[:] = 0.0
    emitted = engine.calculate_radiance(atmosphere)
    atmosphere.storage.emission_source[:] = 0.0
    atmosphere.surface.emission[:] = 0.0
    atmosphere.storage.solar_irradiance[:] = 1.1
    solar = engine.calculate_radiance(atmosphere)
    for name in ["radiance", "wf_extinction", "wf_ssa", "wf_leg_coeff_1", "wf_albedo"]:
        np.testing.assert_allclose(
            combined[name], emitted[name] + solar[name], rtol=1.0e-10, atol=1.0e-12
        )
    for batch_size, threads, derivatives in [(1, 1, True), (4, 3, True), (8, 1, False)]:
        config, geometry, viewing, atmosphere = _case(
            spherical=spherical, derivatives=derivatives
        )
        config.wavelength_batch_size = batch_size
        config.num_threads = threads
        if single_scatter:
            config.single_scatter_source = sk.SingleScatterSource.Exact
        atmosphere.storage.emission_source[::2] = 0.0
        result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)
        for name in result.data_vars:
            np.testing.assert_allclose(
                result[name], combined[name], rtol=1.0e-11, atol=1.0e-12
            )


@pytest.mark.parametrize("spectral_template", [False, True])
def test_volume_emission_constituent_with_rayleigh(spectral_template):
    config, geometry, viewing, _ = _case(spherical=True)
    atmosphere = sk.Atmosphere(
        geometry, config, wavelengths_nm=np.linspace(550.0, 558.0, 9)
    )
    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    altitude = np.array([10_000.0, 20_000.0, 30_000.0])
    ver = np.array([0.0, 1.0, 0.0])
    if spectral_template:
        constituent = sk.constituent.SpectralVolumeEmissionRate(
            altitude, ver, np.array([552.0, 554.0, 556.0]), np.array([0.0, 1.0, 0.0])
        )
        property_name = "photon_ver"
    else:
        constituent = sk.constituent.MonochromaticVolumeEmissionRate(
            altitude, ver, 554.0
        )
        property_name = "ver"
    atmosphere["emitter"] = constituent
    engine = sk.Engine(config, geometry, viewing)
    result = engine.calculate_radiance(atmosphere)
    for level in range(ver.size):
        perturbed = ver.copy()
        perturbed[level] += 1.0e-3
        setattr(constituent, property_name, perturbed)
        plus = engine.calculate_radiance(atmosphere)
        np.testing.assert_allclose(
            result[f"wf_emitter_{property_name}"].isel(emitter_altitude=level),
            (plus.radiance - result.radiance) / 1.0e-3,
            rtol=1.0e-8,
            atol=1.0e-7,
        )
    # Solar illumination remains present outside the emitting line/band.
    assert (result.radiance.isel(wavelength=[0, -1]) > 0.0).all()
