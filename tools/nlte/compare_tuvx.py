"""Compare sasktran2.photolysis with TUV-x on identical inputs.

Builds a CAIRT ERS scenario's atmosphere, runs TUV-x through
``tuvx_reference.py`` in a separate environment that has ``musica``, and
runs sasktran2 with the same profiles, solar zenith angle, albedo and solar
spectrum (sasktran2's, integrated onto the TUV-x wavelength bins). TUV-x's
aerosols are off; it is also run with its own solar spectrum to show that
effect separately.

It also runs :class:`sasktran2.photolysis.TUVActinicFlux`, the TUV mode:
with TUV-x's data and two streams against TUV-x, which isolates the
radiative-transfer solvers, and with sasktran2's data against the full
resolution calculation, which isolates the TUV binning and O2 band
parameterisations.

    python tools/nlte/compare_tuvx.py /Volumes/T9/data/cairt_ers_kopra/ERS_kopra_ascii \\
        april+00 path/to/musica-venv/bin/python out_prefix
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import sasktran2 as sk

sys.path.insert(0, str(Path(__file__).parent))
import validate_photolysis_ers as ers  # noqa: E402

K_BOLTZMANN = 1.380649e-23
ALBEDO = 0.2


def tuvx_wavelength_edges(musica_python: str) -> np.ndarray:
    # Read the edges from a live calculator: the grids returned by
    # v54.wavelength_grid() on their own do not hold valid data.
    code = (
        "import musica.tuvx.v54 as v; import numpy as np; import sys; "
        "t = v.get_tuvx_calculator(); g = t.get_grid_map(); "
        "np.savetxt(sys.stdout, np.array(g['wavelength', 'nm'].edges, copy=True))"
    )
    output = subprocess.run(
        [musica_python, "-c", code], check=True, capture_output=True, text=True
    ).stdout
    return np.array(output.split(), dtype=float)


def binned(wavelength, values, edges) -> np.ndarray:
    """Integrate a fine spectrum over each bin; ``values`` has wavelength first."""
    out = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        fine = np.linspace(lo, hi, 400)
        sampled = np.array(
            [np.interp(fine, wavelength, v) for v in np.atleast_2d(values.T)]
        )
        out.append(np.trapezoid(sampled, fine, axis=-1))
    return np.array(out)


def run_tuvx(musica_python, atmosphere, cos_sza, solar_per_bin, edges, workdir):
    total = atmosphere["pressure_pa"] / (K_BOLTZMANN * atmosphere["temperature_k"])
    inputs = xr.Dataset(
        {
            "temperature_k": atmosphere["temperature_k"],
            "air": total,
            "O2": atmosphere["O2"],
            "O3": atmosphere["O3"],
        },
        attrs={
            "sza_deg": float(np.degrees(np.arccos(cos_sza))),
            "albedo": ALBEDO,
            "earth_sun_distance_au": 1.0,
        },
    )
    if solar_per_bin is not None:
        inputs["solar_flux_per_bin"] = ("wavelength_bin", solar_per_bin)
        inputs = inputs.assign_coords(wavelength_edge=("wavelength_edge", edges))
    inputs_path = Path(workdir) / "tuvx_inputs.nc"
    outputs_path = Path(workdir) / f"tuvx_outputs_{solar_per_bin is None}.nc"
    inputs.to_netcdf(inputs_path)
    script = Path(__file__).parent / "tuvx_reference.py"
    subprocess.run(
        [musica_python, str(script), str(inputs_path), str(outputs_path)], check=True
    )
    return xr.load_dataset(outputs_path)


def main(root: str, scenario: str, musica_python: str, out_prefix: str) -> None:
    root = Path(root)
    month, lat = scenario[:-3], float(scenario[-3:])
    cos_sza = ers.cos_solar_zenith(ers.MONTH_DAY_OF_YEAR[month], lat, 9.5)
    atmosphere = ers.ers_atmosphere(root, scenario, top_km=150.0)

    flux = sk.photolysis.ActinicFlux(atmosphere["altitude"].to_numpy()).calculate(
        atmosphere, cos_sza=cos_sza, albedo=ALBEDO
    )
    total_o2 = [
        sk.photolysis.Photolysis("J_O2_CONT", "O2", wavelength_range_nm=(122.0, 245.0)),
        sk.photolysis.LymanAlphaPhotolysis(
            "J_O2_LYA", sk.photolysis.presets.LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S
        ),
    ]
    ours = sk.photolysis.photolysis_rates(
        flux, [*sk.photolysis.presets.oxygen_photolysis()[:2], *total_o2]
    )
    ours["J_O2"] = ours["J_O2_CONT"] + ours["J_O2_LYA"]

    # TUV mode: TUV-x data and two streams (against TUV-x with its own sun),
    # and sasktran2 data (against the full-resolution calculation).
    altitude = atmosphere["altitude"].to_numpy()
    tuv_mode = {
        "tuv-x": sk.photolysis.photolysis_rates(
            sk.photolysis.TUVActinicFlux(
                altitude, data="tuv-x", num_streams=2
            ).calculate(atmosphere, cos_sza=cos_sza, albedo=ALBEDO),
            sk.photolysis.presets.tuvx_v54_photolysis(),
        ),
        "sasktran2": sk.photolysis.photolysis_rates(
            sk.photolysis.TUVActinicFlux(altitude).calculate(
                atmosphere, cos_sza=cos_sza, albedo=ALBEDO
            ),
            [
                *sk.photolysis.presets.oxygen_photolysis(excitation=False)[:2],
                sk.photolysis.Photolysis("J_O2", "O2"),
            ],
        ),
    }

    edges = tuvx_wavelength_edges(musica_python)
    wavelength = flux["wavelength"].to_numpy()
    solar_per_bin = binned(wavelength, flux["solar_flux"].to_numpy(), edges)[:, 0]

    with tempfile.TemporaryDirectory() as workdir:
        tuvx_same_sun = run_tuvx(
            musica_python, atmosphere, cos_sza, solar_per_bin, edges, workdir
        )
        tuvx_own_sun = run_tuvx(
            musica_python, atmosphere, cos_sza, None, edges, workdir
        )

    z_tuvx = tuvx_same_sun["altitude_km"].to_numpy()
    z_ours = atmosphere["altitude"].to_numpy() / 1e3

    def ours_at(name, rates=ours):
        return np.interp(z_tuvx, z_ours, rates[name].to_numpy())

    ratios = {
        name: {
            "sasktran2 / TUV-x, same sun": ours_at(name) / tuvx_same_sun[name],
            "TUV mode, TUV-x data, 2 streams / TUV-x": (
                ours_at(name, tuv_mode["tuv-x"]) / tuvx_own_sun[name]
            ),
            "TUV mode / sasktran2, sasktran2 data": (
                ours_at(name, tuv_mode["sasktran2"]) / ours_at(name)
            ),
        }
        for name in ("J_O3_O1D", "J_O3_O3P", "J_O2")
    }

    names = {
        "J_O3_O1D": "O3 -> O(1D)",
        "J_O3_O3P": "O3 -> O(3P)",
        "J_O2": "O2, all channels",
    }
    fig, axes = plt.subplots(1, 4, figsize=(17, 5))
    for ax in axes[1:3]:
        ax.sharey(axes[0])
    for ax, (name, label) in zip(axes[:3], names.items(), strict=True):
        for (case, ratio), style in zip(
            ratios[name].items(), ("-", "--", ":"), strict=True
        ):
            ax.plot(ratio, z_tuvx, style, label=case)
        ax.axvline(1.0, color="k", lw=0.5)
        ax.set_xlim(0.5, 1.5)
        ax.set_title(label)
        ax.legend(fontsize="small")
    axes[0].set_ylabel("Altitude [km]")

    # Actinic flux by bin at a few altitudes, same solar spectrum.
    ours_binned = (
        binned(wavelength, flux["actinic_flux"].to_numpy(), edges)
        / np.diff(edges)[:, np.newaxis]
    )
    tuvx_total = tuvx_same_sun["actinic_flux"].sum("component")
    centres = 0.5 * (edges[:-1] + edges[1:])
    for km in (20, 30, 40, 60, 90):
        mine = ours_binned[:, int(np.argmin(abs(z_ours - km)))]
        theirs = tuvx_total.sel(altitude_km=km, method="nearest").to_numpy()
        valid = theirs > 1e-6 * theirs.max()
        axes[3].plot(centres[valid], mine[valid] / theirs[valid], label=f"{km} km")
    axes[3].axhline(1.0, color="k", lw=0.5)
    axes[3].set_ylim(0.5, 1.5)
    axes[3].set_xlabel("Wavelength [nm]")
    axes[3].set_title("Actinic flux, sasktran2 / TUV-x")
    axes[3].legend()
    fig.suptitle(
        f"{scenario}: cos SZA {cos_sza:.3f}, albedo {ALBEDO}, no aerosol, TUV-x v5.4"
    )
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=120)

    print(f"{scenario}: cos SZA {cos_sza:.3f}, albedo {ALBEDO}")
    print(
        "Ratios: (a) sasktran2 / TUV-x, same sun; (b) sasktran2 / TUV-x, TUV-x sun;\n"
        "(c) TUV mode with TUV-x data and 2 streams / TUV-x;\n"
        "(d) TUV mode / full resolution, sasktran2 data"
    )
    print(
        "altitude   O(1D) a / b / c / d             O(3P) a / b / c / d             O2 a / b / c / d"
    )
    for km in (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120):
        i = int(np.argmin(abs(z_tuvx - km)))
        row = []
        for name in names:
            values = [
                ours_at(name)[i] / float(tuvx_same_sun[name][i]),
                ours_at(name)[i] / float(tuvx_own_sun[name][i]),
                *(float(r[i]) for r in list(ratios[name].values())[1:]),
            ]
            row.append(" ".join(f"{v:6.3f}" for v in values))
        print(f"{km:6d} km   " + "    ".join(row))
    print("binned actinic flux sasktran2/TUV-x (same sun), wavelength bands:")
    for lo, hi in (
        (175, 200),
        (200, 242),
        (250, 290),
        (295, 320),
        (320, 400),
        (400, 700),
    ):
        band = (centres >= lo) & (centres < hi)
        row = []
        for km in (20, 30, 40, 60, 90):
            mine = ours_binned[band, int(np.argmin(abs(z_ours - km)))].sum()
            theirs = (
                tuvx_total.sel(altitude_km=km, method="nearest").to_numpy()[band].sum()
            )
            row.append(f"{km}km {mine / theirs:5.3f}")
        print(f"  {lo}-{hi} nm: " + "  ".join(row))


if __name__ == "__main__":
    main(*sys.argv[1:5])
