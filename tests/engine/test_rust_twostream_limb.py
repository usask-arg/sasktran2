from __future__ import annotations

import numpy as np
import pytest
import sasktran2 as sk


def _config(*, solar: bool, thermal: bool, refracted: bool = False) -> sk.Config:
    config = sk.Config()
    config.num_threads = 1
    config.num_streams = 2
    config.num_sza = 3
    config.do_backprop = True
    config.los_refraction = refracted
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = (
        sk.MultipleScatterSource.TwoStream
        if solar
        else sk.MultipleScatterSource.NoSource
    )
    config.emission_source = (
        sk.EmissionSource.TwoStream if thermal else sk.EmissionSource.NoSource
    )
    config.two_stream_backend = sk.TwoStreamBackend.Rust
    return config


def _case(config: sk.Config):
    altitude_grid = np.arange(0.0, 60_001.0, 5_000.0)
    geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.2,
        earth_radius_m=6_371_000.0,
        altitude_grid_m=altitude_grid,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    if config.los_refraction:
        geometry.refractive_index = 1.0 + 3.0e-4 * np.exp(-altitude_grid / 7_000.0)

    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.TangentAltitudeSolar(15_000.0, 0.4, 200_000.0, 0.6))
    viewing.add_ray(sk.TangentAltitudeSolar(35_000.0, -0.3, 200_000.0, 0.6))
    viewing.add_ray(sk.GroundViewingSolar(0.6, 0.3, 0.65, 200_000.0))

    atmosphere = sk.Atmosphere(
        geometry,
        config,
        calculate_derivatives=True,
        numwavel=9,
    )
    spectral = 0.8 + 0.05 * np.arange(9)[np.newaxis, :]
    atmosphere.storage.total_extinction[:] = (
        2.0e-5 * np.exp(-altitude_grid[:, np.newaxis] / 8_000.0) + 1.0e-8
    ) * spectral
    atmosphere.storage.ssa[:] = 0.85
    atmosphere.leg_coeff.a1[0, :, :] = 1.0
    atmosphere.leg_coeff.a1[1, :, :] = 0.2
    atmosphere.surface.albedo[:] = 0.2
    atmosphere.storage.solar_irradiance[:] = 1.1
    if config.emission_source != sk.EmissionSource.NoSource:
        atmosphere.storage.emission_source[:] = (
            2.0 - 0.03 * np.arange(altitude_grid.size)[:, np.newaxis]
        )
        atmosphere.surface.emission[:] = 1.5
    return geometry, viewing, atmosphere


def test_rust_spherical_solar_is_close_to_two_stream_interpolated_do():
    rust_config = _config(solar=True, thermal=False)
    geometry, viewing, atmosphere = _case(rust_config)
    rust = sk.Engine(rust_config, geometry, viewing).calculate_radiance(atmosphere)

    do_config = _config(solar=True, thermal=False)
    do_config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    do_geometry, do_viewing, do_atmosphere = _case(do_config)
    discrete_ordinates = sk.Engine(
        do_config, do_geometry, do_viewing
    ).calculate_radiance(do_atmosphere)

    np.testing.assert_allclose(
        rust["radiance"], discrete_ordinates["radiance"], rtol=2.0e-2, atol=1.0e-12
    )


def test_rust_twostream_delta_m_matches_two_stream_discrete_ordinates():
    def calculate(multiple_scatter_source, backend):
        config = sk.Config()
        config.num_threads = 1
        config.num_streams = 2
        config.num_singlescatter_moments = 4
        config.do_backprop = True
        config.delta_m_scaling = True
        config.single_scatter_source = sk.SingleScatterSource.NoSource
        config.multiple_scatter_source = multiple_scatter_source
        config.two_stream_backend = backend
        config.wavelength_batch_size = 4

        altitude_grid = np.arange(0.0, 40_001.0, 5_000.0)
        geometry = sk.Geometry1D(
            cos_sza=0.6,
            solar_azimuth=0.2,
            earth_radius_m=6_371_000.0,
            altitude_grid_m=altitude_grid,
            interpolation_method=sk.InterpolationMethod.LinearInterpolation,
            geometry_type=sk.GeometryType.PlaneParallel,
        )
        viewing = sk.ViewingGeometry()
        viewing.add_ray(sk.GroundViewingSolar(0.6, 0.3, 0.7, 200_000.0))
        viewing.add_ray(sk.GroundViewingSolar(0.6, -0.4, 0.35, 200_000.0))

        num_wavelengths = 9
        atmosphere = sk.Atmosphere(
            geometry,
            config,
            calculate_derivatives=True,
            numwavel=num_wavelengths,
        )
        spectral = 0.8 + 0.04 * np.arange(num_wavelengths)[np.newaxis, :]
        atmosphere.storage.total_extinction[:] = (
            2.0e-5 * np.exp(-altitude_grid[:, np.newaxis] / 8_000.0) + 1.0e-8
        ) * spectral
        atmosphere.storage.ssa[:] = 0.87
        g = 0.62 + 0.01 * np.arange(num_wavelengths)[np.newaxis, :] / num_wavelengths
        atmosphere.leg_coeff.a1[0, :, :] = 1.0
        atmosphere.leg_coeff.a1[1, :, :] = 3.0 * g
        atmosphere.leg_coeff.a1[2, :, :] = 5.0 * g**2
        atmosphere.leg_coeff.a1[3, :, :] = 7.0 * g**3
        atmosphere.surface.albedo[:] = 0.2
        atmosphere.storage.solar_irradiance[:] = 1.1

        return sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)

    rust = calculate(
        sk.MultipleScatterSource.TwoStream,
        sk.TwoStreamBackend.Rust,
    )
    discrete_ordinates = calculate(
        sk.MultipleScatterSource.DiscreteOrdinates,
        sk.TwoStreamBackend.Cpp,
    )

    assert rust.data_vars.keys() == discrete_ordinates.data_vars.keys()
    for name in rust.data_vars:
        np.testing.assert_allclose(
            rust[name],
            discrete_ordinates[name],
            rtol=2.0e-8,
            atol=2.0e-12,
            err_msg=f"delta-M mismatch for {name}",
        )


def test_rust_spherical_thermal_and_combined_refracted_sources_are_finite():
    separate = []
    for solar, thermal in [(True, False), (False, True)]:
        config = _config(solar=solar, thermal=thermal, refracted=True)
        geometry, viewing, atmosphere = _case(config)
        result = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)
        assert all(np.isfinite(value).all() for value in result.data_vars.values())
        separate.append(result["radiance"])

    config = _config(solar=True, thermal=True, refracted=True)
    geometry, viewing, atmosphere = _case(config)
    combined = sk.Engine(config, geometry, viewing).calculate_radiance(atmosphere)
    assert all(np.isfinite(value).all() for value in combined.data_vars.values())
    np.testing.assert_allclose(
        combined["radiance"], separate[0] + separate[1], rtol=1e-12
    )


@pytest.mark.parametrize("refracted", [False, True])
def test_rust_spherical_multithread_matches_serial(refracted: bool):
    serial_config = _config(solar=True, thermal=True, refracted=refracted)
    geometry, viewing, atmosphere = _case(serial_config)
    serial = sk.Engine(serial_config, geometry, viewing).calculate_radiance(atmosphere)

    parallel_config = _config(solar=True, thermal=True, refracted=refracted)
    parallel_config.num_threads = 4
    parallel_geometry, parallel_viewing, parallel_atmosphere = _case(parallel_config)
    parallel = sk.Engine(
        parallel_config, parallel_geometry, parallel_viewing
    ).calculate_radiance(parallel_atmosphere)

    assert parallel.data_vars.keys() == serial.data_vars.keys()
    for name in serial.data_vars:
        np.testing.assert_allclose(
            parallel[name],
            serial[name],
            rtol=1.0e-12,
            atol=1.0e-13,
            err_msg=f"threading changed {name}",
        )
