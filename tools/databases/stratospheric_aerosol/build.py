"""Build the fixed v1 SAGE-derived catalogue; never run at package import time.

Usage: python tools/databases/stratospheric_aerosol/build.py --psd-root PATH
       --sage-root PATH --output src/sasktran2/_data/stratospheric_aerosol
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import sasktran2 as sk
import xarray as xr
from scipy.integrate import trapezoid
from scipy.ndimage import gaussian_filter1d
from scipy.stats import theilslopes

BANDS = {"sh_midlat": (-55, -35), "tropical": (-20, 20), "nh_midlat": (35, 55)}
TIERS = {
    "low": (0.10, 0.03),
    "typical": (0.50, 0.03),
    "elevated": (0.90, 0.03),
    "extreme": (0.99, 0.005),
}
CHANNELS = [756, 869, 1021, 1543]
BAD_EVENTS = ["2024030913SS", "2024030915SS", "2024030917SS"]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def level_valid(data):
    e = data.aerosol_extinction.sel(nominal_aerosol_wavelength=756).values
    de = data.aerosol_extinction_uncertainty.sel(nominal_aerosol_wavelength=756).values
    r = data.aerosol_median_radius.values
    dr = data.aerosol_median_1sigma.values
    with np.errstate(invalid="ignore", divide="ignore"):
        return (
            np.isfinite(e)
            & (e > 0)
            & np.isfinite(de)
            & (de >= 0)
            & (de / e < 0.5)
            & np.isfinite(r)
            & (r > 10.01)
            & (r < 589.99)
            & np.isfinite(dr)
            & (dr >= 0)
            & (dr / r < 0.5)
        )


def load_sources(psd_root, sage_root):
    parts, sources = [], []
    # Pin the initial catalogue's archive period, independently of future months.
    files = sorted(psd_root.glob("AEROSOL-PSD-L2-SAGE-III-ISS-SASK-V1.0.0-*-fv0001.nc"))
    files = [p for p in files if "201706" <= p.name.split("-")[-2] <= "202607"]
    if len(files) != 110:
        raise ValueError(f"Expected 110 months June 2017–July 2026, found {len(files)}")
    for n, path in enumerate(files):
        month = path.name.split("-")[-2]
        matches = sorted(sage_root.glob(f"g3b_smnc_6.*_{month}.nc"))
        if not matches:
            raise FileNotFoundError(f"No original solar SAGE file for {month}")
        original = matches[-1]
        with xr.open_dataset(path) as psd, xr.open_dataset(original) as sage:
            if float(psd.assumed_mode_width) != 1.6:
                raise ValueError("Expected fixed width 1.6")
            for key, unit in [("aerosol_extinction", "km-1"), ("altitude", "km")]:
                if sage[key].attrs.get("unit") != unit:
                    raise ValueError(f"Unexpected source unit: {key}")
            ids = sage.event_id.values.tolist()
            if len(set(ids)) != len(ids):
                raise ValueError("Duplicate original event IDs")
            lookup = {eid: i for i, eid in enumerate(ids)}
            source = sage.isel(datetime=[lookup[eid] for eid in psd.event_id.values])
            names = [
                "aerosol_median_radius",
                "aerosol_median_1sigma",
                "aerosol_extinction",
                "aerosol_extinction_uncertainty",
                "aerosol_wavelength",
            ]
            part = (
                psd[names]
                .sel(altitude=slice(10, 40), nominal_aerosol_wavelength=CHANNELS)
                .load()
            )
            source = source.sel(
                altitude=part.altitude, nominal_aerosol_wavelength=CHANNELS
            )
            if not np.array_equal(
                part.aerosol_extinction.values,
                source.aerosol_extinction.values,
                equal_nan=True,
            ):
                raise ValueError(f"Source extinction mismatch: {path.name}")
            part["source_flags"] = (
                part.aerosol_extinction.dims,
                source.derived_aerosol_flag.values,
            )
            part["source_index"] = ("datetime", np.full(part.sizes["datetime"], n))
            parts.append(part)
            sources.append(
                {
                    "psd_file": path.name,
                    "psd_sha256": sha256(path),
                    "sage_file": original.name,
                    "sage_sha256": sha256(original),
                    "sage_version": sage.attrs["product_version"],
                }
            )
        if (n + 1) % 10 == 0:
            print(f"Verified {n + 1}/{len(files)} source months", flush=True)
    data = xr.concat(
        parts, dim="datetime", data_vars="minimal", coords="minimal", compat="equals"
    )
    if len(set(data.event_id.values)) != data.sizes["datetime"]:
        raise ValueError("Duplicate PSD event IDs")
    return data, sources


def calibration(data, indices):
    z = data.altitude.values
    e = data.aerosol_extinction.sel(nominal_aerosol_wavelength=756).values[indices]
    flags = data.source_flags.values[indices]
    months = data.datetime.values[indices].astype("datetime64[M]")
    core = (z >= 18) & (z <= 30)
    # Raw fits avoid the artificial smoothing boundary at 30 km. The other
    # fits quantify sensitivity on the same population and background flags.
    upper = (z >= 26) & (z <= 30)
    background = (flags[:, upper] == 2).all(axis=(1, 2))
    e, months = e[background], months[background]
    variants = {
        "raw_27_30": (np.log(e[:, core]), z[core], 27, 30),
        "raw_26_29": (np.log(e[:, core]), z[core], 26, 29),
        "raw_26_30": (np.log(e[:, core]), z[core], 26, 30),
        "smoothed_27_30": (
            gaussian_filter1d(
                np.log(e[:, core]), 1.5 / 2.354820045 / 0.5, axis=1, mode="reflect"
            ),
            z[core],
            27,
            30,
        ),
    }
    result = {}
    for name, (values, heights, lo, hi) in variants.items():
        window = (heights >= lo) & (heights <= hi)
        rates = -np.array(
            [theilslopes(v[window], heights[window]).slope for v in values]
        )
        finite = np.isfinite(rates)
        monthly = []
        # Aggregate before testing the sign, to avoid selecting only individual
        # negative slopes out of noisy upper-level observations.
        for month in np.unique(months):
            group = rates[(months == month) & finite]
            if len(group) >= 5:
                monthly.append(
                    {
                        "month": str(month),
                        "count": len(group),
                        "decay_rate_per_km": float(np.median(group)),
                    }
                )
        if len(monthly) < 24:
            raise ValueError("Insufficient calibration months")
        monthly_rates = np.array([m["decay_rate_per_km"] for m in monthly])
        rate = float(np.median(monthly_rates))
        if rate <= 0:
            raise ValueError("Reference calibration must decay")
        result[name] = {
            "scale_height_m": 1000 / rate,
            "monthly_rates": monthly,
            "n_profiles": len(rates),
            "n_non_decreasing_profiles": int(np.sum(rates <= 0)),
            "monthly_decay_rate_quartiles_per_km": np.quantile(
                monthly_rates, [0.25, 0.5, 0.75]
            ).tolist(),
        }
    return result


def build(psd_root, sage_root, output):
    data, sources = load_sources(psd_root, sage_root)
    z = data.altitude.values
    core = (z >= 18) & (z <= 30)
    ec = data.aerosol_extinction.sel(nominal_aerosol_wavelength=756).values[:, core]
    rc = data.aerosol_median_radius.values[:, core]
    valid = (
        level_valid(data)[:, core].all(axis=1)
        & np.isin(data.source_flags.values[:, core], [2, 3]).all(axis=(1, 2))
        & np.isfinite(data.aerosol_tropopause_height.values)
        & (data.aerosol_tropopause_height.values + 1 <= 18)
        & ~np.isin(data.event_id.values, BAD_EVENTS)
    )
    tau = trapezoid(ec, z[core], axis=1)
    # Check spectral closure with the same sulfate optical assumptions as PSD retrieval.
    mie = sk.database.MieDatabase(
        sk.mie.LogNormalDistribution().freeze(mode_width=1.6),
        sk.mie.refractive.H2SO4(),
        np.arange(300.0, 1600.0, 10.0),
        median_radius=np.arange(10.0, 600.0, 10.0),
    )
    with xr.open_dataset(mie.path()) as table:
        xs = (
            table.xs_total.interp(wavelength_nm=CHANNELS)
            .transpose("median_radius", "wavelength_nm")
            .values
        )
        radii = table.median_radius.values
    pred = np.stack(
        [
            np.interp(rc, radii, xs[:, j]) / np.interp(rc, radii, xs[:, 0])
            for j in range(1, 4)
        ],
        axis=-1,
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        ratios = data.aerosol_extinction.values[:, core, 1:] / ec[:, :, None]
    # Only evaluate complete, accepted profiles; invalid data are not repaired.
    residual = np.full(len(valid), np.nan)
    residual[valid] = np.median(
        np.abs(ratios[valid] - pred[valid]) / pred[valid], axis=(1, 2)
    )
    valid &= residual < 0.3
    chosen, bands = [], {}
    lat = data.latitude.values
    for band, (south, north) in BANDS.items():
        pool = np.flatnonzero(valid & (lat >= south) & (lat <= north))
        lo, hi = np.quantile(tau[pool], [0.1, 0.6])
        regular = pool[(tau[pool] >= lo) & (tau[pool] <= hi)]
        cal = calibration(data, regular)
        # Round to 100 m: precision is limited by variability and methodology.
        height = float(np.round(cal["raw_27_30"]["scale_height_m"] / 100) * 100)
        bands[band] = {
            "population": len(pool),
            "regular_loading_aod_cutoffs": [float(lo), float(hi)],
            "reference_upper_scale_height_m": height,
            "calibration": cal,
        }
        print(
            band,
            "accepted",
            len(pool),
            "reference H",
            height,
            "sensitivity",
            {k: round(v["scale_height_m"]) for k, v in cal.items()},
            flush=True,
        )
        for tier, (quantile, half_width) in TIERS.items():
            lo, hi = np.quantile(
                tau[pool], [quantile - half_width, quantile + half_width]
            )
            group = pool[(tau[pool] >= lo) & (tau[pool] <= hi)]
            le, radius = np.log(ec[group]), rc[group]
            score = np.mean(((le - np.median(le, axis=0)) / 0.5) ** 2, axis=1)
            score += np.mean(((radius - np.median(radius, axis=0)) / 50) ** 2, axis=1)
            index = int(group[np.argmin(score)])
            chosen.append(
                {
                    "scenario": f"{band}_{tier}",
                    "index": index,
                    "latitude_band": band,
                    "loading": tier,
                    "selection_aod_18_30": float(tau[index]),
                    "loading_percentile": float(100 * np.mean(tau[pool] <= tau[index])),
                    "candidate_count": len(group),
                    "spectral_fractional_residual": float(residual[index]),
                    "reference_upper_scale_height_m": height,
                }
            )
    profiles = []
    for case in chosen:
        index = case.pop("index")
        source = sources[int(data.source_index[index])]
        eid = str(data.event_id[index].item())
        with (
            xr.open_dataset(psd_root / source["psd_file"]) as psd,
            xr.open_dataset(sage_root / source["sage_file"]) as sage,
        ):
            p = psd.isel(
                datetime=int(np.flatnonzero(psd.event_id.values == eid)[0])
            ).load()
            s = sage.isel(
                datetime=int(np.flatnonzero(sage.event_id.values == eid)[0])
            ).sel(altitude=p.altitude)
            flag = s.derived_aerosol_flag.sel(
                nominal_aerosol_wavelength=CHANNELS
            ).values
            good = (
                level_valid(p)
                & np.isin(flag, [2, 3]).all(axis=1)
                & (p.altitude.values >= float(p.aerosol_tropopause_height) + 1)
            )
            # Extend the complete common core only through contiguous valid data.
            a = int(np.flatnonzero(p.altitude.values == 18)[0])
            b = int(np.flatnonzero(p.altitude.values == 30)[0])
            while a > 0 and good[a - 1]:
                a -= 1
            while b + 1 < len(good) and good[b + 1]:
                b += 1
            keep = slice(a, b + 1)
            p = p.isel(altitude=keep)
            if not np.array_equal(
                p.aerosol_extinction.values,
                s.aerosol_extinction.isel(altitude=keep).values,
                equal_nan=True,
            ):
                raise ValueError("Selected full-spectrum source mismatch")
            alt = p.altitude.values.astype(float) * 1000
            dims = ("altitude_m", "wavelength_nm")
            profile = xr.Dataset(
                {
                    "raw_extinction_per_m": (
                        dims,
                        p.aerosol_extinction.values.astype(float) / 1000,
                        {"units": "m-1"},
                    ),
                    "raw_extinction_uncertainty_per_m": (
                        dims,
                        p.aerosol_extinction_uncertainty.values.astype(float) / 1000,
                        {"units": "m-1"},
                    ),
                    "raw_median_radius_nm": (
                        "altitude_m",
                        p.aerosol_median_radius.values,
                        {"units": "nm"},
                    ),
                    "raw_median_radius_uncertainty_nm": (
                        "altitude_m",
                        p.aerosol_median_1sigma.values,
                        {"units": "nm"},
                    ),
                    "source_aerosol_flag": (
                        dims,
                        s.derived_aerosol_flag.isel(altitude=keep).values,
                    ),
                    "actual_wavelength_nm": (
                        "wavelength_nm",
                        p.aerosol_wavelength.values,
                        {"units": "nm"},
                    ),
                    "observed_valid": ("altitude_m", np.ones(len(alt), dtype=np.int8)),
                },
                coords={
                    "altitude_m": ("altitude_m", alt, {"units": "m"}),
                    "wavelength_nm": (
                        "wavelength_nm",
                        p.nominal_aerosol_wavelength.values,
                        {"units": "nm"},
                    ),
                },
            )
            case.update(
                event_id=eid,
                datetime=str(p.datetime.values),
                latitude_degrees=float(p.latitude),
                longitude_degrees=float(p.longitude),
                tropopause_altitude_m=float(p.aerosol_tropopause_height) * 1000,
                observed_bottom_m=float(alt[0]),
                observed_top_m=float(alt[-1]),
                **source,
            )
            for key, value in case.items():
                if key != "scenario":
                    profile[key] = xr.DataArray(value)
            profile = profile.expand_dims(scenario=[case["scenario"]])
            profiles.append(profile)
            print(case["scenario"], eid, "valid bounds", alt[[0, -1]], flush=True)
    catalogue = xr.concat(profiles, dim="scenario", join="outer")
    catalogue["observed_valid"] = catalogue.observed_valid.fillna(0).astype(np.int8)
    catalogue.attrs.update(
        catalogue_version="v1",
        schema_version=1,
        reference_wavelength_nm=756.0,
        mode_width=1.6,
        source="USask SAGE III-ISS particle size v1.0.0 fv0001; NASA SAGE III-ISS solar v6",
        archive_period="2017-06/2026-07",
        selection_interval_m=[18000.0, 30000.0],
        calibration_method="Median raw log-extinction Theil-Sen decay rate per month (>=5 profiles), median across months; H rounded to 100 m",
        source_radius_units="nm, verified against USask retrieval implementation; missing source unit attributes",
        limitations="Fixed-width sulfate reference cases, not a climatology for arbitrary dates or smoke optical properties; formal errors omit systematic/model uncertainty",
    )
    for key in [
        "observed_bottom_m",
        "observed_top_m",
        "tropopause_altitude_m",
        "reference_upper_scale_height_m",
    ]:
        catalogue[key].attrs["units"] = "m"
    output.mkdir(parents=True, exist_ok=True)
    destination = output / "stratospheric_aerosol_v1.nc"
    catalogue.to_netcdf(
        destination,
        engine="netcdf4",
        encoding={
            k: {"zlib": True, "complevel": 6}
            for k, v in catalogue.data_vars.items()
            if v.dtype.kind not in "OU"
        },
    )
    report = {
        "catalogue_sha256": sha256(destination),
        "spectral_check_optical_table_file": mie.path().name,
        "spectral_check_optical_table_sha256": sha256(mie.path()),
        "archive_profiles": data.sizes["datetime"],
        "screened_profiles": int(valid.sum()),
        "sources": sources,
        "bands": bands,
        "scenarios": chosen,
        "calibration_background_flags": "flag 2 in all four near-IR channels throughout 26–30 km",
        "selection_screen": "All 18–30 km levels: positive extinction; finite radius 10.01–589.99 nm; nonnegative formal errors below 50%; flags 2/3 in four near-IR channels; 1 km above aerosol tropopause; spectral residual <30%",
        "selection_weighting": "Observation-weighted loading quantiles; actual paired profile nearest candidate median log extinction and radius",
        "references": [
            "https://doi.org/10.1029/JD092iD03p03051",
            "https://doi.org/10.1364/AO.8.000893",
            "https://doi.org/10.5194/amtd-5-5993-2012",
        ],
    }
    (output / "build_report_v1.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        "Wrote",
        destination,
        destination.stat().st_size,
        "bytes",
        report["catalogue_sha256"],
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--psd-root", type=Path, required=True)
    parser.add_argument("--sage-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.psd_root, args.sage_root, args.output)
