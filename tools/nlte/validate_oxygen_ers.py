"""Compare daytime oxygen populations from sasktran2 with GRANADA (CAIRT ERS).

Runs :func:`sasktran2.nlte.add_photochemical_species` on an ERS day
scenario, with the scenario's temperature, pressure, O2, O3, N2, CO2 and
atomic oxygen, and compares O(1D), O2(a), O2(b) and O2(b, v=1) with the
GRANADA results in the archive (``npar`` and ``ratio`` files).

GRANADA gives O2 states as r = n / n_LTE, with n_LTE normalised over the
modelled X(v=0-35), a(v=0-5) and b(v=0-2) states using a degeneracy and
energy per state from a level file that is not in the archive. Energies here
come from the spectroscopic constants of Huber and Herzberg (1979); the
degeneracy is the electronic one (X 3, a 2, b 1) unless ``--no-degeneracy``
is given, which uses 1 for every state.

    python tools/nlte/validate_oxygen_ers.py /Volumes/T9/data/cairt_ers_kopra/ERS_kopra_ascii \\
        april+00 out.png [--no-degeneracy]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import sasktran2 as sk

sys.path.insert(0, str(Path(__file__).parent))
import kopra_prf  # noqa: E402
import validate_photolysis_ers as ers  # noqa: E402

K_BOLTZMANN = 1.380649e-23
C2_CM_K = 1.4387769
#: Te, we, wexe [cm^-1] and electronic degeneracy (Huber and Herzberg 1979).
O2_CONSTANTS = {
    "X": (0.0, 1580.19, 11.98, 3),
    "a": (7918.1, 1483.50, 12.90, 2),
    "b": (13195.1, 1432.77, 14.00, 1),
}
STATES = {"O2(b)": "b0", "O2(b, v=1)": "b1", "O2(a)": "a0"}


def granada_o2(root: Path, scenario: str, altitude_m, temperature_k, n_o2, degeneracy):
    """GRANADA O2 state densities [m^-3] by 'b0', 'a0', ... on ``altitude_m``."""
    ratio = kopra_prf.read_ratio(root / f"{scenario}_ratio.prf")
    o2 = ratio.where(ratio.species_id == 71, drop=True).interp(
        altitude_km=altitude_m / 1e3
    )
    labels = [str(label).split()[:2] for label in o2["label"].to_numpy()]

    def term(state, v):
        te, we, wexe, _ = O2_CONSTANTS[state]
        return te + we * (v + 0.5) - wexe * (v + 0.5) ** 2

    energy = np.array([term(s, int(v)) - term("X", 0) for s, v in labels])
    g = np.array([O2_CONSTANTS[s][3] if degeneracy else 1.0 for s, _ in labels])
    weights = g[:, np.newaxis] * np.exp(
        -C2_CM_K * energy[:, np.newaxis] / temperature_k[np.newaxis, :]
    )
    density = o2["ratio"].to_numpy() * weights / weights.sum(axis=0) * n_o2
    return {f"{s}{v}": density[i] for i, (s, v) in enumerate(labels)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("root", type=Path)
    parser.add_argument("scenario")
    parser.add_argument("output")
    parser.add_argument("--albedo", type=float, default=0.2)
    parser.add_argument("--no-degeneracy", action="store_true")
    args = parser.parse_args()

    month, lat = args.scenario[:-3], float(args.scenario[-3:])
    day = ers.MONTH_DAY_OF_YEAR[month]
    cos_sza = ers.cos_solar_zenith(day, lat, 9.5)
    source = ers.ers_atmosphere(args.root, args.scenario, top_km=150.0)
    vmr = kopra_prf.read_vmr(args.root / f"{args.scenario}_vmr.prf")
    altitude = source["altitude"].to_numpy()
    temperature = source["temperature_k"].to_numpy()
    air = source["pressure_pa"].to_numpy() / (K_BOLTZMANN * temperature)

    def density(species):
        return (
            vmr["vmr_ppmv"].sel(species=species).interp(altitude_km=altitude / 1e3)
            * 1e-6
            * air
        ).to_numpy()

    background = xr.Dataset(
        {name: ("altitude", density(name)) for name in ("O", "O2", "O3", "N2", "CO2")},
        coords={"altitude": altitude},
    )
    config = sk.Config()
    config.emission_source = sk.EmissionSource.VolumeEmissionRate
    geometry = sk.Geometry1D(
        cos_sza,
        0.0,
        6_372_000.0,
        altitude,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=np.array([762.0]))
    atmosphere.temperature_k = temperature
    atmosphere.pressure_pa = source["pressure_pa"].to_numpy()
    solution = sk.nlte.add_photochemical_species(
        atmosphere,
        ["O2(b)"],
        cos_sza=cos_sza,
        background=background,
        albedo=args.albedo,
        earth_sun_distance_au=ers.earth_sun_distance_au(day),
    )

    granada = granada_o2(
        args.root,
        args.scenario,
        altitude,
        temperature,
        density("O2"),
        degeneracy=not args.no_degeneracy,
    )
    npar = kopra_prf.read_npar(args.root / f"{args.scenario}_npar.prf").interp(
        altitude_km=altitude / 1e3
    )
    reference = {
        "O(1D)": npar["O1D_D"].to_numpy() * 1e6,
        **{state: granada[key] for state, key in STATES.items()},
    }

    z_km = altitude / 1e3
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    print(f"{args.scenario}: cos SZA {cos_sza:.3f}; sasktran2 / GRANADA")
    print("altitude  " + "  ".join(f"{name:>11s}" for name in reference))
    for km in (40, 50, 60, 70, 80, 90, 100, 110):
        i = int(np.argmin(abs(z_km - km)))
        cells = [
            float(solution["density"].sel(state=name)[i]) / values[i]
            for name, values in reference.items()
        ]
        print(f"{km:5d} km  " + "  ".join(f"{c:11.3f}" for c in cells))
    for name, values in reference.items():
        ours = solution["density"].sel(state=name).to_numpy()
        (line,) = axes[0].semilogx(ours * 1e-6, z_km, label=name)
        axes[0].semilogx(values * 1e-6, z_km, "--", color=line.get_color())
        axes[1].plot(ours / values, z_km, color=line.get_color(), label=name)
    axes[0].set_xlabel("Density [cm^-3]; solid sasktran2, dashed GRANADA")
    axes[0].set_ylabel("Altitude [km]")
    axes[0].set_ylim(30, 120)
    axes[0].legend()
    axes[1].axvline(1.0, color="k", lw=0.5)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("sasktran2 / GRANADA")
    fig.suptitle(
        f"{args.scenario}: cos SZA {cos_sza:.3f}, mechanism {solution.attrs['mechanism']}"
    )
    fig.tight_layout()
    fig.savefig(args.output, dpi=120)


if __name__ == "__main__":
    main()
