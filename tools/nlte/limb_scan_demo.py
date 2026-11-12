"""Daytime limb scan, 280-800 nm, with photochemical emission.

An OSIRIS-like scan through a CAIRT ERS day atmosphere: Rayleigh-scattered
sunlight with O3, O2 and NO2 absorption, and O2(b) (A, B and gamma bands),
OH A-X fluorescence (with OH self-absorption), O(1S) (557.7 and 297.2 nm)
and O(1D) (630.0 and 636.4 nm) emission from
sasktran2.nlte.add_photochemical_species. No instrument convolution.

    python tools/nlte/limb_scan_demo.py /Volumes/T9/data/cairt_ers_kopra/ERS_kopra_ascii april+00 out_prefix
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import sasktran2 as sk

sys.path.insert(0, str(Path(__file__).parent))
import kopra_prf  # noqa: E402
import validate_photolysis_ers as ers  # noqa: E402

K_BOLTZMANN = 1.380649e-23
SPECIES = ["O2(b)", "OH(A)", "O(1S)", "O(1D)"]
TANGENTS_KM = np.arange(10.0, 100.1, 5.0)


def main(root: str, scenario: str, out_prefix: str) -> None:
    root = Path(root)
    month, lat = scenario[:-3], float(scenario[-3:])
    day = ers.MONTH_DAY_OF_YEAR[month]
    cos_sza = ers.cos_solar_zenith(day, lat, 9.5)
    source = ers.ers_atmosphere(root, scenario, top_km=150.0)
    vmr = kopra_prf.read_vmr(root / f"{scenario}_vmr.prf")
    z_src = source["altitude"].to_numpy()
    air = source["pressure_pa"].to_numpy() / (
        K_BOLTZMANN * source["temperature_k"].to_numpy()
    )

    def density(species):
        return (
            vmr["vmr_ppmv"].sel(species=species).interp(altitude_km=z_src / 1e3)
            * 1e-6
            * air
        ).to_numpy()

    background = xr.Dataset(
        {name: ("altitude", density(name)) for name in ("O", "OH", "CO2")},
        coords={"altitude": z_src},
    )

    altitude = np.arange(0.0, 150_001.0, 1000.0)
    wavelength = sk.nlte.emission_wavelength_grid(SPECIES)
    config = sk.Config()
    config.emission_source = sk.EmissionSource.VolumeEmissionRate
    config.num_threads = 8
    geometry = sk.Geometry1D(
        cos_sza,
        0.0,
        6_372_000.0,
        altitude,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    viewing = sk.ViewingGeometry()
    for tangent in TANGENTS_KM:
        viewing.add_ray(sk.TangentAltitudeSolar(tangent * 1e3, 0.0, 600e3, cos_sza))

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavelength)
    atmosphere.temperature_k = np.interp(
        altitude, z_src, source["temperature_k"].to_numpy()
    )
    atmosphere.pressure_pa = np.exp(
        np.interp(altitude, z_src, np.log(source["pressure_pa"].to_numpy()))
    )
    optics = {
        "O3": sk.optical.O3DBM(),
        "O2": sk.optical.HITRANAbsorber("O2"),
        "NO2": sk.optical.NO2Vandaele(),
    }
    for name, optical in optics.items():
        atmosphere[name] = sk.constituent.VMRAltitudeAbsorber(
            optical, z_src, source[name].to_numpy() / air
        )
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["solar"] = sk.constituent.SolarIrradiance(photon_units=True)
    engine = sk.Engine(config, geometry, viewing)

    start = time.time()
    clear = (
        engine.calculate_radiance(atmosphere)["radiance"]
        .to_numpy()
        .reshape(wavelength.size, -1)
    )
    t_clear = time.time() - start
    start = time.time()
    solution = sk.nlte.add_photochemical_species(
        atmosphere,
        SPECIES,
        cos_sza=cos_sza,
        background=background,
        albedo=0.2,
        earth_sun_distance_au=ers.earth_sun_distance_au(day),
    )
    t_chem = time.time() - start
    start = time.time()
    total = (
        engine.calculate_radiance(atmosphere)["radiance"]
        .to_numpy()
        .reshape(wavelength.size, -1)
    )
    t_total = time.time() - start
    print(
        f"{wavelength.size} wavelengths x {TANGENTS_KM.size} tangents: "
        f"scattering {t_clear:.0f} s, photochemistry {t_chem:.0f} s, with emission {t_total:.0f} s"
    )

    result = xr.Dataset(
        {
            "radiance": (("wavelength", "tangent_altitude_km"), total),
            "radiance_no_emission": (("wavelength", "tangent_altitude_km"), clear),
        },
        coords={"wavelength": wavelength, "tangent_altitude_km": TANGENTS_KM},
        attrs={
            "units": "photons m^-2 s^-1 nm^-1 sr^-1",
            "scenario": scenario,
            "cos_sza": cos_sza,
        },
    )
    result.to_netcdf(f"{out_prefix}.nc")
    solution.drop_vars(
        [v for v in solution.data_vars if solution[v].dtype == object]
    ).to_netcdf(f"{out_prefix}_photochemistry.nc")

    # Bin to 1 nm for display only (not an instrument model).
    edges = np.arange(280.0, 800.01, 1.0)
    centres = 0.5 * (edges[1:] + edges[:-1])

    def binned(values):
        widths = np.gradient(wavelength)
        index = np.digitize(wavelength, edges) - 1
        ok = (index >= 0) & (index < centres.size)
        out = np.zeros((centres.size, values.shape[1]))
        for j in range(values.shape[1]):
            out[:, j] = np.bincount(
                index[ok], weights=(values[ok, j] * widths[ok]), minlength=centres.size
            )
        return out

    total_b, clear_b = binned(total), binned(clear)
    fig, axes = plt.subplots(2, 1, figsize=(12, 9))
    for km in (30, 50, 70, 90):
        j = int(np.argmin(abs(TANGENTS_KM - km)))
        axes[0].semilogy(centres, total_b[:, j], label=f"{km} km")
        axes[1].plot(centres, total_b[:, j] / clear_b[:, j], label=f"{km} km")
    axes[0].set_ylabel("Radiance, 1 nm bins [photons m-2 s-1 sr-1]")
    axes[0].legend()
    axes[1].set_ylabel("With emission / without")
    axes[1].set_xlabel("Wavelength [nm]")
    axes[1].set_yscale("log")
    axes[1].legend()
    fig.suptitle(
        f"{scenario}: cos SZA {cos_sza:.2f}; emission from O2(b), OH(A), O(1S), O(1D)"
    )
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=120)


if __name__ == "__main__":
    main(*sys.argv[1:4])
