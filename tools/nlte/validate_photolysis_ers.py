"""Compare sasktran2.photolysis rates with the TUV rates of the CAIRT ERS archive.

The ERS ``npar`` files hold photolysis rates computed with KIT's TUV 5.4
(MIPAS preprocessor) for each scenario. This script rebuilds a scenario's
atmosphere from its ``pt`` and ``vmr`` files, runs
:class:`sasktran2.photolysis.ActinicFlux`, and plots both sets of rates.

The archive does not record solar zenith angles or surface albedo. Day
scenarios are at 9.5 h local solar time; the zenith angle and Earth-Sun
distance are computed for the 15th of the month, which matters only where the
atmosphere is optically thick. The albedo (default 0.2) mainly affects the
visible Chappuis band and so the O3 -> O(3P) rate.

    python tools/nlte/validate_photolysis_ers.py /Volumes/T9/data/cairt_ers_kopra/ERS_kopra_ascii april+00 out.png [albedo]
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import sasktran2 as sk

sys.path.insert(0, str(Path(__file__).parent))
import kopra_prf  # noqa: E402

K_BOLTZMANN = 1.380649e-23
MONTH_DAY_OF_YEAR = {"january": 15, "april": 105, "july": 196, "october": 288}


def cos_solar_zenith(
    day_of_year: int, latitude_deg: float, local_time_h: float
) -> float:
    declination = np.radians(23.44) * np.sin(2 * np.pi * (284 + day_of_year) / 365)
    hour_angle = np.radians(15.0 * (local_time_h - 12.0))
    lat = np.radians(latitude_deg)
    return float(
        np.sin(lat) * np.sin(declination)
        + np.cos(lat) * np.cos(declination) * np.cos(hour_angle)
    )


def earth_sun_distance_au(day_of_year: int) -> float:
    return float(1.0 - 0.01672 * np.cos(np.radians(0.9856 * (day_of_year - 4))))


def ers_atmosphere(root: Path, scenario: str, top_km: float = 150.0) -> xr.Dataset:
    pt = kopra_prf.read_pt(root / f"{scenario}_pt.prf")
    vmr = kopra_prf.read_vmr(root / f"{scenario}_vmr.prf")
    keep = pt.altitude_km <= top_km
    pt, vmr = pt.sel(altitude_km=keep), vmr.sel(altitude_km=keep)
    pressure = pt["pressure_hpa"].to_numpy() * 100.0
    temperature = pt["temperature_k"].to_numpy()
    total = pressure / (K_BOLTZMANN * temperature)
    data = {
        "temperature_k": ("altitude", temperature),
        "pressure_pa": ("altitude", pressure),
    }
    for species in ("O2", "O3", "N2", "NO2"):
        data[species] = (
            "altitude",
            vmr["vmr_ppmv"].sel(species=species).to_numpy() * 1e-6 * total,
        )
    return xr.Dataset(data, coords={"altitude": pt.altitude_km.to_numpy() * 1e3})


def main(root: str, scenario: str, output: str, albedo: str = "0.2") -> None:
    root = Path(root)
    month, lat = scenario[:-3], float(scenario[-3:])
    day = MONTH_DAY_OF_YEAR[month]
    cos_sza = cos_solar_zenith(day, lat, 9.5)

    atmosphere = ers_atmosphere(root, scenario)
    flux = sk.photolysis.ActinicFlux(atmosphere["altitude"].to_numpy()).calculate(
        atmosphere,
        cos_sza=cos_sza,
        albedo=float(albedo),
        earth_sun_distance_au=earth_sun_distance_au(day),
    )
    total_o2 = [
        sk.photolysis.Photolysis(
            "J_O2_CONTINUUM", "O2", wavelength_range_nm=(122.0, 245.0)
        ),
        sk.photolysis.LinePhotolysis(
            "J_O2_LYA_TOTAL",
            sk.photolysis.LYMAN_ALPHA_WAVELENGTH_NM,
            sk.photolysis.presets.O2_LYMAN_ALPHA_CROSS_SECTION_M2,
            sk.photolysis.presets.LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S,
        ),
    ]
    ours = sk.photolysis.photolysis_rates(
        flux, [*sk.photolysis.presets.oxygen_photolysis(), *total_o2]
    )
    ours["J_O2"] = ours["J_O2_CONTINUUM"] + ours["J_O2_LYA_TOTAL"]
    ers = kopra_prf.read_npar(root / f"{scenario}_npar.prf").interp(
        altitude_km=atmosphere["altitude"].to_numpy() / 1e3
    )

    z = atmosphere["altitude"].to_numpy() / 1e3
    pairs = [
        ("O3 -> O(1D)", ours["J_O3_O1D"], ers["J_O3_1"]),
        ("O3 -> O(3P)", ours["J_O3_O3P"], ers["J_O3_2"]),
        ("O2 photolysis, all channels", ours["J_O2"], ers["J_O2"]),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(13, 9), sharey=True)
    for column, (label, mine, theirs) in enumerate(pairs):
        top, bottom = axes[0, column], axes[1, column]
        top.semilogx(theirs, z, label="ERS (TUV 5.4)")
        top.semilogx(mine, z, "--", label="sasktran2.photolysis")
        top.set_title(label)
        top.set_xlabel("J [s$^{-1}$]")
        top.legend()
        bottom.plot(mine.to_numpy() / theirs.to_numpy(), z)
        bottom.axvline(1.0, color="k", lw=0.5)
        bottom.set_xlim(0.0, 2.0)
        bottom.set_xlabel("sasktran2 / ERS")
    axes[0, 0].set_ylabel("Altitude [km]")
    axes[1, 0].set_ylabel("Altitude [km]")
    fig.suptitle(f"{scenario} day: cos SZA = {cos_sza:.3f}, albedo = {albedo}")
    fig.tight_layout()
    fig.savefig(output, dpi=120)

    print(f"{scenario}: cos SZA {cos_sza:.3f}, albedo {albedo}")
    print("altitude  O(1D) ours/ERS  O(3P) ours/ERS  O2 ours/ERS")
    for km in (20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120):
        i = int(np.argmin(abs(z - km)))
        print(
            f"{km:6d} km  {float(ours['J_O3_O1D'][i] / ers['J_O3_1'][i]):13.3f}"
            f"  {float(ours['J_O3_O3P'][i] / ers['J_O3_2'][i]):14.3f}"
            f"  {float(ours['J_O2'][i] / ers['J_O2'][i]):11.3f}"
        )


if __name__ == "__main__":
    main(*sys.argv[1:5])
