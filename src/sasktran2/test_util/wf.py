from __future__ import annotations

from matplotlib.pylab import f
import numpy as np
import xarray as xr

import sasktran2 as sk


def validate_wf(analytic, numerical, wf_dim="altitude", decimal=6):
    max_by_alt = np.abs(analytic).max(dim=wf_dim)

    max_by_alt.to_numpy()[max_by_alt.to_numpy() == 0] = 1e99

    rel_diff = (analytic - numerical) / max_by_alt

    nonzero_analytic = np.abs(analytic.to_numpy()) > 1e-99
    nonzero_numerical = np.abs(numerical.to_numpy()) > 1e-99

    np.testing.assert_array_almost_equal(
        rel_diff.to_numpy()[nonzero_analytic & nonzero_numerical], 0, decimal=decimal
    )


def numeric_wf(
    input_var: np.array,
    fractional_change: float,
    engine: sk.Engine,
    atmosphere: sk.Atmosphere,
    analytic_wf_name: str,
    calc_vars=None,
) -> xr.Dataset:
    if calc_vars is None:
        calc_vars = ["radiance"]

    base_radiance = engine.calculate_radiance(atmosphere)

    central_diff_wf = {v: np.zeros_like(base_radiance[analytic_wf_name + ("" if v == "radiance" else f"_{v}")].to_numpy()) for v in calc_vars}

    for i in range(len(input_var)):
        dx = input_var[i] * fractional_change

        if dx == 0:
            dx = np.nanmean(input_var) * fractional_change

        input_var[i] += dx
        radiance_above = engine.calculate_radiance(atmosphere)

        if input_var[i] >= dx:
            # central diff
            input_var[i] -= 2 * dx
            radiance_below = engine.calculate_radiance(atmosphere)
            input_var[i] += dx

            for v in calc_vars:
                central_diff_wf[v][i] = (
                    radiance_above[v] - radiance_below[v]
                ) / (2 * dx)
        else:
            # forward diff
            for v in calc_vars:
                central_diff_wf[v][i] = (
                    radiance_above[v] - base_radiance[v]
                ) / (dx)

            input_var[i] -= dx
    for v in calc_vars:
        append_name = "" if v == "radiance" else f"_{v}"
        analytic_wf_name_full = analytic_wf_name + append_name
        base_radiance[analytic_wf_name_full + "_numeric"] = (
            base_radiance[analytic_wf_name_full].dims,
            central_diff_wf[v],
        )

    return base_radiance


def numeric_wf_scalar(
    input_var: np.ndarray,
    fractional_change: float,
    engine: sk.Engine,
    atmosphere: sk.Atmosphere,
    analytic_wf_name: str,
) -> xr.Dataset:
    base_radiance = engine.calculate_radiance(atmosphere)

    dx = input_var * fractional_change

    if dx == 0:
        dx = np.nanmean(input_var) * fractional_change

    input_var += dx
    radiance_above = engine.calculate_radiance(atmosphere)

    if input_var >= dx:
        # central diff
        input_var -= 2 * dx
        radiance_below = engine.calculate_radiance(atmosphere)
        input_var += dx

        central_diff_wf = (radiance_above["radiance"] - radiance_below["radiance"]) / (
            2 * dx
        )
    else:
        # forward diff
        central_diff_wf = (radiance_above["radiance"] - base_radiance["radiance"]) / (
            dx
        )

        input_var -= dx

    base_radiance[analytic_wf_name + "_numeric"] = central_diff_wf

    return base_radiance
