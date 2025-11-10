from __future__ import annotations

import numpy as np
import sasktran2 as sk
from sasktran2.test_util.wf import validate_wf


def _raw_scenarios() -> list:
    """ """
    configs = []

    configs.append(sk.Config())
    configs[-1].multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    configs[-1].num_streams = 4
    configs[-1].single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    configs[-1].num_stokes = 1

    geometrys = []
    altitude_grid = np.arange(0, 65001, 5000.0)

    geometrys.append(
        sk.Geometry1D(
            0.6,
            0,
            6327000,
            altitude_grid,
            sk.InterpolationMethod.LinearInterpolation,
            sk.GeometryType.PlaneParallel,
        )
    )

    viewing_geos = []
    viewing_geos.append(sk.ViewingGeometry())

    for alt in np.arange(10000, 60000, 2000):
        viewing_geos[-1].add_flux_observer(sk.FluxObserverSolar(0.6, alt))

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

    return scen


def test_wf_extinction():
    """
    Verifies that the extinction derivative is accurate to 5 decimal places
    """
    D_FRACTION = 1e-4
    test_scens = _raw_scenarios()

    flux_names = ["upwelling_flux", "downwelling_flux"]

    for scen in test_scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = engine.calculate_radiance(scen["atmosphere"])

        numeric_wf = {}
        for flux_name in flux_names:
            numeric_wf[flux_name] = np.zeros_like(
                radiance[f"wf_extinction_{flux_name}"].values
            )

        for i in range(len(radiance["altitude"])):
            dk = D_FRACTION * scen["atmosphere"].storage.total_extinction[i]

            scen["atmosphere"].storage.total_extinction[i] += dk
            radiance_above = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].storage.total_extinction[i] -= 2 * dk
            radiance_below = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].storage.total_extinction[i] += dk

            for flux_name in flux_names:
                numeric_wf[flux_name][i, :, :] = (
                    radiance_above[flux_name].to_numpy()
                    - radiance_below[flux_name].to_numpy()
                ) / (2 * dk)

        for flux_name in flux_names:
            radiance[f"wf_extinction_{flux_name}_numeric"] = (
                ["altitude", "wavelength", "flux_location"],
                numeric_wf[flux_name],
            )

            validate_wf(
                radiance[f"wf_extinction_{flux_name}"],
                radiance[f"wf_extinction_{flux_name}_numeric"],
                decimal=5,
            )


def test_wf_ssa():
    """
    Verifies that the ssa derivative is accurate to 5 decimal places
    """
    D_FRACTION = 1e-4
    test_scens = _raw_scenarios()

    flux_names = ["upwelling_flux", "downwelling_flux"]

    for scen in test_scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = engine.calculate_radiance(scen["atmosphere"])

        numeric_wf = {}
        for flux_name in flux_names:
            numeric_wf[flux_name] = np.zeros_like(
                radiance[f"wf_ssa_{flux_name}"].values
            )

        for i in range(len(radiance["altitude"])):
            dk = D_FRACTION * scen["atmosphere"].storage.ssa[i]

            scen["atmosphere"].storage.ssa[i] += dk
            radiance_above = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].storage.ssa[i] -= 2 * dk
            radiance_below = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].storage.ssa[i] += dk

            for flux_name in flux_names:
                numeric_wf[flux_name][i, :, :] = (
                    radiance_above[flux_name].to_numpy()
                    - radiance_below[flux_name].to_numpy()
                ) / (2 * dk)

        for flux_name in flux_names:
            radiance[f"wf_ssa_{flux_name}_numeric"] = (
                ["altitude", "wavelength", "flux_location"],
                numeric_wf[flux_name],
            )

            validate_wf(
                radiance[f"wf_ssa_{flux_name}"],
                radiance[f"wf_ssa_{flux_name}_numeric"],
                decimal=5,
            )


def test_wf_legendre():
    """
    Verifies that the legendre coefficient derivatives are accurate to some decimal places
    """
    D_L = 1e-4
    test_scens = _raw_scenarios()

    flux_names = ["upwelling_flux", "downwelling_flux"]

    for scen in test_scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = engine.calculate_radiance(scen["atmosphere"])

        for i in range(scen["atmosphere"].storage.leg_coeff.shape[0]):
            if i == 0:
                continue

            numeric_wf = {}
            for flux_name in flux_names:
                numeric_wf[flux_name] = np.zeros_like(
                    radiance[f"wf_leg_coeff_{i}_{flux_name}"].values
                )

            for j in range(len(radiance["altitude"])):
                scen["atmosphere"].storage.leg_coeff[i, j] += D_L
                radiance_above = engine.calculate_radiance(scen["atmosphere"])

                scen["atmosphere"].storage.leg_coeff[i, j] -= 2 * D_L
                radiance_below = engine.calculate_radiance(scen["atmosphere"])

                scen["atmosphere"].storage.leg_coeff[i, j] += D_L

                for flux_name in flux_names:
                    numeric_wf[flux_name][j, :, :] = (
                        radiance_above[flux_name].to_numpy()
                        - radiance_below[flux_name].to_numpy()
                    ) / (2 * D_L)

            for flux_name in flux_names:
                radiance[f"wf_leg_coeff_{i}_{flux_name}_numeric"] = (
                    ["altitude", "wavelength", "flux_location"],
                    numeric_wf[flux_name],
                )
                validate_wf(
                    radiance[f"wf_leg_coeff_{i}_{flux_name}"],
                    radiance[f"wf_leg_coeff_{i}_{flux_name}_numeric"],
                    decimal=5,
                )


def test_wf_albedo():
    D_ALBEDO = 1e-5
    test_scens = _raw_scenarios()

    flux_names = ["upwelling_flux", "downwelling_flux"]

    for scen in test_scens:
        engine = sk.Engine(scen["config"], scen["geometry"], scen["viewing_geo"])

        radiance = engine.calculate_radiance(scen["atmosphere"])

        numeric_wf = {}
        for flux_name in flux_names:
            numeric_wf[flux_name] = np.zeros_like(
                radiance[f"wf_albedo_{flux_name}"].values
            )

        if scen["atmosphere"].surface.albedo[0] > D_ALBEDO:
            # Do central difference
            scen["atmosphere"].surface.albedo[:] += D_ALBEDO
            radiance_above = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].surface.albedo[:] -= 2 * D_ALBEDO
            radiance_below = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].surface.albedo[:] += D_ALBEDO

            for flux_name in flux_names:
                numeric_wf[flux_name][:, :] = (
                    radiance_above[flux_name].to_numpy()
                    - radiance_below[flux_name].to_numpy()
                ) / (2 * D_ALBEDO)

            val_decimal = 6
        else:
            # Have to do forward difference
            scen["atmosphere"].surface.albedo[:] += D_ALBEDO
            radiance_above = engine.calculate_radiance(scen["atmosphere"])

            scen["atmosphere"].surface.albedo[:] -= D_ALBEDO

            for flux_name in flux_names:
                numeric_wf[flux_name][:, :] = (
                    radiance_above[flux_name].to_numpy()
                    - radiance[flux_name].to_numpy()
                ) / (D_ALBEDO)

            val_decimal = 4

        for flux_name in flux_names:
            radiance[f"wf_albedo_{flux_name}_numeric"] = (
                ["wavelength", "flux_location"],
                numeric_wf[flux_name],
            )
            validate_wf(
                radiance[f"wf_albedo_{flux_name}"],
                radiance[f"wf_albedo_{flux_name}_numeric"],
                decimal=val_decimal,
                wf_dim="wavelength",
            )
