"""Break the sasktran2.photolysis / TUV-x differences down by cause.

Uses the same inputs as ``compare_tuvx.py`` but only O2, O3 and Rayleigh
scattering, which is what TUV-x v5.4 includes:

1. Solar spectra: TUV-x's built-in extraterrestrial flux against TSIS-1 HSRS,
   per TUV-x bin.
2. Photolysis rates by wavelength region. TUV-x rates are linear in the input
   solar flux, so masking the flux to one region at a time gives that
   region's contribution; sasktran2 rates are integrated over the same
   ranges.
3. Direct-beam transmission per bin: TUV-x's direct component divided by its
   input flux, against sasktran2's exact flux-weighted transmission computed
   from its own cross sections. Differences here are cross sections.
4. Diffuse actinic flux (total minus direct) per bin: differences here are
   the radiative-transfer solver and the Rayleigh cross section.

    python tools/nlte/diagnose_tuvx.py /Volumes/T9/data/cairt_ers_kopra/ERS_kopra_ascii \\
        april+00 path/to/musica-venv/bin/python out_prefix
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import sasktran2 as sk
from sasktran2.optical.rayleigh import rayleigh_cross_section_bates

sys.path.insert(0, str(Path(__file__).parent))
import compare_tuvx as cmp  # noqa: E402
import validate_photolysis_ers as ers  # noqa: E402

K_BOLTZMANN = 1.380649e-23
LYMAN_ALPHA_BIN = (121.4, 121.9)
#: Regions by TUV-x bin edges [nm]; O2 and O3 rates are split over all of them.
REGIONS = {
    "Lyman-alpha": LYMAN_ALPHA_BIN,
    "SR continuum": (120.0, 175.4),
    "SR bands": (175.4, 204.1),
    "Herzberg": (204.1, 243.902),
    "Hartley": (243.902, 302.5),
    "Huggins": (302.5, 342.5),
    "Chappuis": (342.5, 735.0),
}
ALTITUDES_KM = (20, 30, 40, 50, 60, 70, 80, 90, 100, 110)
#: Longward of this, O2 absorption is excitation, which TUV-x does not count
#: as photolysis.
O2_DISSOCIATION_NM = 243.902


def tuvx_solar_per_bin(musica_python: str) -> np.ndarray:
    code = (
        "import musica.tuvx.v54 as v; import numpy as np; import sys; "
        "t = v.get_tuvx_calculator(); np.savetxt(sys.stdout, np.asarray("
        "t.get_profile_map()['extraterrestrial flux', 'photon cm-2 s-1'].midpoint_values))"
    )
    out = cmp.subprocess.run(
        [musica_python, "-c", code], check=True, capture_output=True, text=True
    ).stdout
    return np.array(out.split(), dtype=float) * 1.0e4  # photons m^-2 s^-1 per bin


def region_mask(edges: np.ndarray, lo: float, hi: float, exclude=None) -> np.ndarray:
    lower, upper = edges[:-1], edges[1:]
    mask = (lower >= lo - 1e-9) & (upper <= hi + 1e-9)
    if exclude is not None:
        mask &= ~((lower >= exclude[0] - 1e-9) & (upper <= exclude[1] + 1e-9))
    return mask


def sasktran2_region_rates(flux: xr.Dataset) -> dict[str, xr.Dataset]:
    """O3 channel and total O2 rates of each region."""
    out = {}
    for name, (lo, hi) in REGIONS.items():
        reactions = []
        if name == "Lyman-alpha":
            reactions.append(
                sk.photolysis.LymanAlphaPhotolysis(
                    "J_O2", sk.photolysis.presets.LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S
                )
            )
        elif hi <= O2_DISSOCIATION_NM:
            # The SR continuum region excludes the Lyman-alpha bin, which the
            # line reaction above covers.
            lo_eff = LYMAN_ALPHA_BIN[1] if name == "SR continuum" else lo
            reactions.append(
                sk.photolysis.Photolysis("J_O2", "O2", wavelength_range_nm=(lo_eff, hi))
            )
        reactions += [
            sk.photolysis.Photolysis(
                "J_O3_O1D", "O3", sk.photolysis.o3_o1d_matsumi2002, (lo, hi)
            ),
            sk.photolysis.Photolysis(
                "J_O3_O3P", "O3", sk.photolysis.o3_o3p_matsumi2002, (lo, hi)
            ),
        ]
        rates = xr.Dataset(
            {"J_O2": ("altitude", np.zeros(flux["altitude"].size))},
            coords={"altitude": flux["altitude"]},
        )
        for reaction in reactions:
            try:
                rates[reaction.name] = sk.photolysis.photolysis_rates(flux, [reaction])[
                    reaction.name
                ]
            except ValueError:
                rates[reaction.name] = ("altitude", np.zeros(flux["altitude"].size))
        out[name] = rates
    return out


def direct_transmission(flux, atmosphere, cos_sza) -> np.ndarray:
    """Exact direct-beam transmission (wavelength, altitude) from sasktran2's cross sections."""
    z = flux["altitude"].to_numpy()
    wavelength = flux["wavelength"].to_numpy()
    total = atmosphere["pressure_pa"] / (K_BOLTZMANN * atmosphere["temperature_k"])
    rayleigh, _ = rayleigh_cross_section_bates(wavelength / 1.0e3)
    extinction = rayleigh[:, np.newaxis] * total.to_numpy()[np.newaxis, :]
    for species in flux["species"].to_numpy():
        extinction = extinction + (
            flux["cross_section"].sel(species=species).to_numpy()
            * atmosphere[species].to_numpy()[np.newaxis, :]
        )
    # Optical depth above each level, linear extinction between levels.
    layer = 0.5 * (extinction[:, 1:] + extinction[:, :-1]) * np.diff(z)[np.newaxis, :]
    above = np.concatenate(
        [np.cumsum(layer[:, ::-1], axis=1)[:, ::-1], np.zeros((wavelength.size, 1))],
        axis=1,
    )
    transmission = np.exp(-above / cos_sza)
    # The Lyman-alpha line is not resolved on this grid; use the Chabrillat
    # and Kockarts line reduction factor (Rayleigh is negligible there).
    i_lya = int(np.argmin(abs(wavelength - sk.photolysis.LYMAN_ALPHA_WAVELENGTH_NM)))
    transmission[i_lya] = sk.photolysis.lyman_alpha_reduction_factor(
        flux["slant_column"].sel(species="O2").to_numpy()
    )
    return transmission


def main(root: str, scenario: str, musica_python: str, out_prefix: str) -> None:
    root = Path(root)
    month, lat = scenario[:-3], float(scenario[-3:])
    cos_sza = ers.cos_solar_zenith(ers.MONTH_DAY_OF_YEAR[month], lat, 9.5)
    atmosphere = ers.ers_atmosphere(root, scenario, top_km=120.0)[
        ["temperature_k", "pressure_pa", "O2", "O3"]
    ]
    z_km = atmosphere["altitude"].to_numpy() / 1e3

    flux = sk.photolysis.ActinicFlux(atmosphere["altitude"].to_numpy()).calculate(
        atmosphere, cos_sza=cos_sza, albedo=cmp.ALBEDO
    )
    wavelength = flux["wavelength"].to_numpy()
    edges = cmp.tuvx_wavelength_edges(musica_python)
    centres = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    hsrs_per_bin = cmp.binned(wavelength, flux["solar_flux"].to_numpy(), edges)[:, 0]
    # The Lyman-alpha line is a single grid point; give its bin the line flux.
    lya_bin = region_mask(edges, *LYMAN_ALPHA_BIN)
    hsrs_per_bin[lya_bin] = sk.photolysis.presets.LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S

    # 1. Solar spectra
    tuvx_sun = tuvx_solar_per_bin(musica_python)

    # 2. Region contributions
    ours_regions = sasktran2_region_rates(flux)
    tuvx_regions = {}
    with tempfile.TemporaryDirectory() as workdir:
        for name, (lo, hi) in REGIONS.items():
            exclude = LYMAN_ALPHA_BIN if name == "SR continuum" else None
            masked = np.where(region_mask(edges, lo, hi, exclude), hsrs_per_bin, 0.0)
            tuvx_regions[name] = cmp.run_tuvx(
                musica_python, atmosphere, cos_sza, masked, edges, workdir
            )
        tuvx_full = cmp.run_tuvx(
            musica_python, atmosphere, cos_sza, hsrs_per_bin, edges, workdir
        )

    z_tuvx = tuvx_full["altitude_km"].to_numpy()

    def at(values, km):
        return float(np.interp(km, z_km, np.asarray(values)))

    def tuvx_at(values, km):
        return float(np.interp(km, z_tuvx, np.asarray(values)))

    print(
        f"{scenario}: cos SZA {cos_sza:.3f}, albedo {cmp.ALBEDO}, O2 + O3 + Rayleigh only"
    )

    print("\n1. Solar spectrum, TUV-x / HSRS integrated over each region:")
    for name, (lo, hi) in REGIONS.items():
        exclude = LYMAN_ALPHA_BIN if name == "SR continuum" else None
        m = region_mask(edges, lo, hi, exclude)
        print(f"   {name:14s} {tuvx_sun[m].sum() / hsrs_per_bin[m].sum():6.3f}")

    print("\n2. Rate contributions by region at 40 / 70 / 90 km: fraction of the")
    print("   sasktran2 total, and sasktran2 / TUV-x for that region")
    for rate in ("J_O2", "J_O3_O1D", "J_O3_O3P"):
        print(f"   {rate}")
        totals = {
            km: sum(at(ours_regions[n][rate], km) for n in REGIONS)
            for km in (40, 70, 90)
        }
        for name in REGIONS:
            cells = []
            for km in (40, 70, 90):
                mine = at(ours_regions[name][rate], km)
                theirs = tuvx_at(tuvx_regions[name][rate], km)
                ratio = mine / theirs if theirs > 0 else np.nan
                cells.append(f"{mine / totals[km]:6.1%} {ratio:6.3f}")
            print(f"     {name:14s} " + "   ".join(cells))

    # 3. Direct-beam transmission per bin
    t_direct = direct_transmission(flux, atmosphere, cos_sza)
    sun = flux["solar_flux"].to_numpy()
    ours_direct = cmp.binned(wavelength, sun[:, np.newaxis] * t_direct, edges)
    ours_direct_t = ours_direct / cmp.binned(wavelength, sun, edges)[:, :1]
    ours_direct_t[lya_bin] = t_direct[
        int(np.argmin(abs(wavelength - sk.photolysis.LYMAN_ALPHA_WAVELENGTH_NM)))
    ]
    tuvx_direct_t = (
        tuvx_full["actinic_flux"].sel(component="direct").to_numpy()
        / (hsrs_per_bin / widths)[:, np.newaxis]
    )

    print("\n3. Direct-beam transmission, sasktran2 / TUV-x, by region and altitude")
    for name, (lo, hi) in REGIONS.items():
        exclude = LYMAN_ALPHA_BIN if name == "SR continuum" else None
        m = region_mask(edges, lo, hi, exclude)
        cells = []
        for km in (40, 60, 70, 80, 90, 100):
            i_ours = int(np.argmin(abs(z_km - km)))
            i_tuvx = int(np.argmin(abs(z_tuvx - km)))
            weights = hsrs_per_bin[m]
            mine = (weights * ours_direct_t[m, i_ours]).sum()
            theirs = (weights * tuvx_direct_t[m, i_tuvx]).sum()
            cells.append(f"{km}km {mine / theirs if theirs > 1e-30 else np.nan:7.3f}")
        print(f"   {name:14s} " + " ".join(cells))

    # 4. Diffuse actinic flux per bin
    total_binned = cmp.binned(wavelength, flux["actinic_flux"].to_numpy(), edges)
    ours_diffuse = total_binned - ours_direct
    tuvx_diffuse = (
        tuvx_full["actinic_flux"]
        .sel(component=["upwelling", "downwelling"])
        .sum("component")
        * widths[:, np.newaxis]
    ).to_numpy()
    print("\n4. Diffuse actinic flux, sasktran2 / TUV-x, by region and altitude")
    for name, (lo, hi) in list(REGIONS.items())[3:]:
        m = region_mask(edges, lo, hi)
        cells = []
        for km in (20, 40, 60, 90, 110):
            i_ours = int(np.argmin(abs(z_km - km)))
            i_tuvx = int(np.argmin(abs(z_tuvx - km)))
            cells.append(
                f"{km}km {ours_diffuse[m, i_ours].sum() / tuvx_diffuse[m, i_tuvx].sum():6.3f}"
            )
        print(f"   {name:14s} " + " ".join(cells))

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    axes[0, 0].plot(centres, tuvx_sun / hsrs_per_bin)
    axes[0, 0].axhline(1.0, color="k", lw=0.5)
    axes[0, 0].set_ylim(0.6, 1.6)
    axes[0, 0].set_title("Solar flux per bin, TUV-x / TSIS-1 HSRS")
    axes[0, 0].set_xlabel("Wavelength [nm]")

    for name in ("Lyman-alpha", "SR continuum", "SR bands", "Herzberg"):
        mine = np.interp(z_tuvx, z_km, ours_regions[name]["J_O2"].to_numpy())
        theirs = tuvx_regions[name]["J_O2"].to_numpy()
        ok = theirs > 1e-3 * np.nanmax(theirs)
        axes[0, 1].plot(mine[ok] / theirs[ok], z_tuvx[ok], label=name)
    axes[0, 1].axvline(1.0, color="k", lw=0.5)
    axes[0, 1].set_xlim(0.0, 3.0)
    axes[0, 1].set_title("O2 photolysis by region, sasktran2 / TUV-x")
    axes[0, 1].set_ylabel("Altitude [km]")
    axes[0, 1].legend()

    for km in (60, 70, 80, 90, 100):
        i_ours = int(np.argmin(abs(z_km - km)))
        i_tuvx = int(np.argmin(abs(z_tuvx - km)))
        m = (centres < 245.0) & (tuvx_direct_t[:, i_tuvx] > 1e-6)
        axes[1, 0].plot(
            centres[m],
            ours_direct_t[m, i_ours] / tuvx_direct_t[m, i_tuvx],
            ".-",
            label=f"{km} km",
        )
    axes[1, 0].axhline(1.0, color="k", lw=0.5)
    axes[1, 0].set_ylim(0.0, 3.0)
    axes[1, 0].set_title("Direct transmission, sasktran2 / TUV-x (< 245 nm)")
    axes[1, 0].set_xlabel("Wavelength [nm]")
    axes[1, 0].legend()

    for km in (20, 40, 60, 90):
        i_ours = int(np.argmin(abs(z_km - km)))
        i_tuvx = int(np.argmin(abs(z_tuvx - km)))
        m = (centres > 250.0) & (
            tuvx_diffuse[:, i_tuvx] > 1e-6 * tuvx_diffuse[:, i_tuvx].max()
        )
        axes[1, 1].plot(
            centres[m],
            ours_diffuse[m, i_ours] / tuvx_diffuse[m, i_tuvx],
            label=f"{km} km",
        )
    axes[1, 1].axhline(1.0, color="k", lw=0.5)
    axes[1, 1].set_ylim(0.6, 1.4)
    axes[1, 1].set_title("Diffuse actinic flux, sasktran2 / TUV-x (> 250 nm)")
    axes[1, 1].set_xlabel("Wavelength [nm]")
    axes[1, 1].legend()

    fig.suptitle(f"{scenario}: sources of the sasktran2 / TUV-x differences")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=120)


if __name__ == "__main__":
    main(*sys.argv[1:5])
