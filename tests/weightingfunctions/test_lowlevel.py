from __future__ import annotations

import numpy as np
import sasktran2 as sk
from sasktran2.test_util.wf import validate_wf


def _ground_scenarios() -> list:
    scens = []

    altitude_grid = np.arange(0, 65001, 5000.0)
    viewing_geo = sk.ViewingGeometry()
    viewing_geo.add_ray(sk.GroundViewingSolar(0.6, 0, 0.8, 200000))

    # Single scatter spherical
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact

    geometry = sk.Geometry1D(
        0.6,
        0,
        6372000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )

    scens.append(
        {
            "config": config,
            "geometry": geometry,
            "viewing_geo": viewing_geo,
            "atmosphere": sk.test_util.scenarios.default_pure_scattering_atmosphere(
                config, geometry, 0.8, albedo=0.5
            ),
        }
    )

    # Multiple scatter spherical
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.num_sza = 2

    geometry = sk.Geometry1D(
        0.6,
        0,
        6372000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )

    scens.append(
        {
            "config": config,
            "geometry": geometry,
            "viewing_geo": viewing_geo,
            "atmosphere": sk.test_util.scenarios.default_pure_scattering_atmosphere(
                config, geometry, 0.8, albedo=0.5
            ),
        }
    )

    # Two stream Source
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.TwoStream

    config.num_streams = 2

    geometry = sk.Geometry1D(
        0.6,
        0,
        6372000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.PlaneParallel,
    )

    scens.append(
        {
            "config": config,
            "geometry": geometry,
            "viewing_geo": viewing_geo,
            "atmosphere": sk.test_util.scenarios.default_pure_scattering_atmosphere(
                config, geometry, 0.8, albedo=0.5
            ),
        }
    )

    # Two stream Source with pseudo-spherical
    config = sk.Config()
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.TwoStream

    config.num_streams = 2

    geometry = sk.Geometry1D(
        0.01,
        0,
        6372000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.PseudoSpherical,
    )

    scens.append(
        {
            "config": config,
            "geometry": geometry,
            "viewing_geo": viewing_geo,
            "atmosphere": sk.test_util.scenarios.default_pure_scattering_atmosphere(
                config, geometry, 0.8, albedo=0.5
            ),
        }
    )

    return scens


def _raw_scenarios() -> list:
    """
    Defines the set of scenarios that we later on test dI/dk, dI/dssa, dI/dleg_coeff, dI/dalbedo on

    Right now we test for

    Configurations

    - Single scatter only
    - DO Source

    Viewing geometries

    - Limb viewing
    - Nadir viewing

    Model Geometry

    - 1D spherical

    Atmosphere

    - ficticious atmosphere with only Rayleigh scattering

    """
    configs = []

    configs.append(sk.Config())
    configs[-1].multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    configs[-1].num_streams = 4
    configs[-1].single_scatter_source = sk.SingleScatterSource.NoSource
    configs[-1].num_stokes = 3

    configs.append(sk.Config())
    configs[-1].multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    configs[-1].num_streams = 4

    configs.append(sk.Config())

    for config in configs:
        config.num_singlescatter_moments = 4

    geometrys = []
    altitude_grid = np.arange(0, 65001, 5000.0)

    geometrys.append(
        sk.Geometry1D(
            0.6,
            0,
            6327000,
            altitude_grid,
            sk.InterpolationMethod.LinearInterpolation,
            sk.GeometryType.Spherical,
        )
    )

    viewing_geos = []
    viewing_geos.append(sk.ViewingGeometry())

    for tan_alt in np.arange(10000, 60000, 2000):
        viewing_geos[-1].add_ray(sk.TangentAltitudeSolar(tan_alt, 0, 600000, 0.6))

    viewing_geos.append(sk.ViewingGeometry())
    viewing_geos[-1].add_ray(sk.GroundViewingSolar(0.6, 0, 0.8, 200000))

    scen = []

    for config in configs:
        for geometry in geometrys:
            for viewing_geo in viewing_geos:
                scen.append(
                    {
                        "config": config,
                        "geometry": geometry,
                        "viewing_geo": viewing_geo,
                        "atmosphere": sk.test_util.scenarios.default_pure_scattering_atmosphere(
                            config, geometry, 0.8
                        ),
                    }
                )

    # Add an additional scenario that is for plane parallel geometry and the DO source
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates

    geometry = sk.Geometry1D(
        0.6,
        0,
        6372000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.PseudoSpherical,
    )

    scen.append(
        {
            "config": config,
            "geometry": geometry,
            "viewing_geo": viewing_geos[-1],
            "atmosphere": sk.test_util.scenarios.default_pure_scattering_atmosphere(
                config, geometry, 0.8, albedo=0.5
            ),
        }
    )

    # Add an additional scenario that is for plane parallel geometry and the DO source, using backprop
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    config.do_backprop = True
    config.num_streams = 4
    config.num_singlescatter_moments = 4

    geometry = sk.Geometry1D(
        0.6,
        0,
        6372000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.PseudoSpherical,
    )

    scen.append(
        {
            "config": config,
            "geometry": geometry,
            "viewing_geo": viewing_geos[-1],
            "atmosphere": sk.test_util.scenarios.default_pure_scattering_atmosphere(
                config, geometry, 0.8, albedo=0.2
            ),
        }
    )

    # Add an additional scenario that is for plane parallel geometry and the DO source, using backprop
    # And with nstokes=3
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    config.do_backprop = True
    config.num_stokes = 3
    config.num_streams = 4
    config.num_singlescatter_moments = 4

    geometry = sk.Geometry1D(
        0.6,
        0,
        6372000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.PseudoSpherical,
    )

    scen.append(
        {
            "config": config,
            "geometry": geometry,
            "viewing_geo": viewing_geos[-1],
            "atmosphere": sk.test_util.scenarios.default_pure_scattering_atmosphere(
                config, geometry, 0.8, albedo=0.2
            ),
        }
    )

    # two stream backprop
    config = sk.Config()
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    config.do_backprop = True
    config.num_stokes = 1
    config.num_streams = 2
    config.num_singlescatter_moments = 4

    geometry = sk.Geometry1D(
        0.6,
        0,
        6372000,
        altitude_grid,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.PseudoSpherical,
    )

    scen.append(
        {
            "config": config,
            "geometry": geometry,
            "viewing_geo": viewing_geos[-1],
            "atmosphere": sk.test_util.scenarios.default_pure_scattering_atmosphere(
                config, geometry, 0.8, albedo=0.2
            ),
        }
    )

    return scen + _ground_scenarios()


def test_wf_extinction():
    """
    Verifies that the extinction derivative is accurate to 5 decimal places
    """
    D_FRACTION = 1e-4
    test_scens = _raw_scenarios()

    for scen in test_scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = engine.calculate_radiance(scen["atmosphere"])

        numeric_wf = np.zeros_like(radiance["wf_extinction"].values)

        for i in range(len(radiance["altitude"])):
            dk = D_FRACTION * scen["atmosphere"].storage.total_extinction[i]

            scen["atmosphere"].storage.total_extinction[i] += dk
            radiance_above = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].storage.total_extinction[i] -= 2 * dk
            radiance_below = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].storage.total_extinction[i] += dk

            numeric_wf[i, :, :, :] = (
                radiance_above["radiance"].to_numpy()
                - radiance_below["radiance"].to_numpy()
            ) / (2 * dk)

        radiance["wf_extinction_numeric"] = (
            ["altitude", "wavelength", "los", "stokes"],
            numeric_wf,
        )
        validate_wf(
            radiance["wf_extinction"], radiance["wf_extinction_numeric"], decimal=5
        )


def test_wf_ssa():
    """
    Verifies that the single scatter albedo derivative is accurate to 6 decimal places
    """
    D_FRACTION = 1e-5
    test_scens = _raw_scenarios()

    for scen in test_scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = engine.calculate_radiance(scen["atmosphere"])

        numeric_wf = np.zeros_like(radiance["wf_ssa"].values)

        for i in range(len(radiance["altitude"])):
            d_ssa = D_FRACTION * scen["atmosphere"].storage.ssa[i]

            scen["atmosphere"].storage.ssa[i] += d_ssa
            radiance_above = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].storage.ssa[i] -= 2 * d_ssa
            radiance_below = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].storage.ssa[i] += d_ssa

            numeric_wf[i, :, :, :] = (
                radiance_above["radiance"].to_numpy()
                - radiance_below["radiance"].to_numpy()
            ) / (2 * d_ssa)

        radiance["wf_ssa_numeric"] = (
            ["altitude", "wavelength", "los", "stokes"],
            numeric_wf,
        )
        validate_wf(radiance["wf_ssa"], radiance["wf_ssa_numeric"], decimal=6)


def test_wf_legendre():
    """
    Verifies that the legendre coefficient derivatives are accurate to some decimal places
    """
    D_L = 1e-4
    test_scens = _raw_scenarios()

    for scen in test_scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = engine.calculate_radiance(scen["atmosphere"])

        for i in range(scen["atmosphere"].storage.leg_coeff.shape[0]):
            if i == 0:
                continue

            numeric_wf = np.zeros_like(radiance[f"wf_leg_coeff_{i}"].values)

            for j in range(len(radiance["altitude"])):
                scen["atmosphere"].storage.leg_coeff[i, j] += D_L
                radiance_above = engine.calculate_radiance(scen["atmosphere"])

                scen["atmosphere"].storage.leg_coeff[i, j] -= 2 * D_L
                radiance_below = engine.calculate_radiance(scen["atmosphere"])

                scen["atmosphere"].storage.leg_coeff[i, j] += D_L

                numeric_wf[j, :, :, :] = (
                    radiance_above["radiance"].to_numpy()
                    - radiance_below["radiance"].to_numpy()
                ) / (2 * D_L)

            radiance[f"wf_leg_coeff_{i}_numeric"] = (
                ["altitude", "wavelength", "los", "stokes"],
                numeric_wf,
            )
            validate_wf(
                radiance[f"wf_leg_coeff_{i}"],
                radiance[f"wf_leg_coeff_{i}_numeric"],
                decimal=6,
            )


def test_wf_albedo():
    D_ALBEDO = 1e-5
    test_scens = _raw_scenarios()

    for scen in test_scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = engine.calculate_radiance(scen["atmosphere"])

        numeric_wf = np.zeros_like(radiance["wf_albedo"].values)

        if scen["atmosphere"].surface.albedo[0] > D_ALBEDO:
            # Do central difference
            scen["atmosphere"].surface.albedo[:] += D_ALBEDO
            radiance_above = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].surface.albedo[:] -= 2 * D_ALBEDO
            radiance_below = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].surface.albedo[:] += D_ALBEDO

            numeric_wf[:, :, :] = (
                radiance_above["radiance"].to_numpy()
                - radiance_below["radiance"].to_numpy()
            ) / (2 * D_ALBEDO)

            val_decimal = 6
        else:
            # Have to do forward difference
            scen["atmosphere"].surface.albedo[:] += D_ALBEDO
            radiance_above = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].surface.albedo[:] -= D_ALBEDO

            numeric_wf[:, :, :] = (
                radiance_above["radiance"].to_numpy() - radiance["radiance"].to_numpy()
            ) / (D_ALBEDO)

            val_decimal = 4

        radiance["wf_albedo_numeric"] = (
            ["wavelength", "los", "stokes"],
            numeric_wf,
        )
        validate_wf(
            radiance["wf_albedo"],
            radiance["wf_albedo_numeric"],
            decimal=val_decimal,
            wf_dim="wavelength",
        )
