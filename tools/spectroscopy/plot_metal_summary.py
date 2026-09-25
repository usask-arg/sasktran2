# ruff: noqa: T201, RUF001
"""Plot the installed metal spectroscopy and free-VER templates.

Run with the project Python environment. No downloads or database writes occur.
The broad-band cross-section views use integrated line areas and an illustrative
Gaussian display kernel, not an instrument response or a radiance calculation.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import html
import json
import warnings
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sasktran2 as sk
import xarray as xr
from matplotlib.colors import Normalize
from sasktran2.optical.resonance import LineResonance, resonance_phase_function
from scipy.constants import atomic_mass, c, k
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.ndimage import gaussian_filter1d

INK = "#172b40"
MUTED = "#566779"
BLUE = "#17689c"
TEAL = "#078177"
ORANGE = "#ce672b"
PURPLE = "#84529b"
RED = "#b44151"
COLORS = [BLUE, TEAL, ORANGE, PURPLE, RED]
MOLECULES = ["AlO", "MgO", "CaO", "TiO"]
PRIORITY = [
    "Na_I",
    "K_I",
    "Li_I",
    "Mg_I",
    "Mg_II",
    "Ca_I",
    "Ca_II",
    "Fe_I",
    "Ni_I",
    "Al_I",
    "Cr_I",
    "Mn_I",
    "Ti_I",
    "Co_I",
    "Cu_I",
    "Zn_I",
    "Si_I",
    "Fe_II",
    "V_I",
    "Rb_I",
    "Sr_I",
    "Ba_I",
]
DX = 0.02
GRID = np.arange(260.0, 830.0 + DX / 2, DX)
WINDOW = (GRID >= 274 - DX / 4) & (GRID <= 810 + DX / 4)
WAVELENGTH = GRID[WINDOW]


def label(species):
    return species.replace("_", " ")


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "text.color": INK,
            "axes.labelcolor": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#c0cbd5",
            "axes.titleweight": "bold",
            "axes.labelsize": 11,
            "axes.titlesize": 13,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "legend.frameon": False,
            "legend.fontsize": 10,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "axes.formatter.use_mathtext": True,
        }
    )


def decorate(ax, *, spectral=True):
    ax.grid(True, color="#e6edf2", lw=0.65)
    ax.set_axisbelow(True)
    if spectral:
        ax.set_xlim(274, 810)
        ax.set_xticks([280, 400, 500, 600, 700, 800])
        ax.set_xlabel("Vacuum wavelength [nm]")


def heading(fig, title, subtitle):
    fig.text(0.065, 0.971, title, fontsize=23, fontweight="bold", va="top")
    fig.text(0.065, 0.936, subtitle, fontsize=11.5, color=MUTED, va="top")


def footer(fig, text):
    fig.text(
        0.065, 0.019, text, fontsize=9.5, color=MUTED, va="bottom", linespacing=1.55
    )


def broaden(wavelength, area, fwhm):
    """Deposit each line's m² nm area, then apply a unit-area display kernel."""
    position = (wavelength - GRID[0]) / DX
    left = np.floor(position).astype(int)
    fraction = position - left
    if np.any(left < 0) or np.any(left + 1 >= len(GRID)):
        msg = "Plot grid does not enclose the source lines"
        raise ValueError(msg)
    values = np.bincount(left, weights=area * (1 - fraction), minlength=len(GRID))
    values += np.bincount(left + 1, weights=area * fraction, minlength=len(GRID))
    result = gaussian_filter1d(
        values / DX,
        fwhm / np.sqrt(8 * np.log(2)) / DX,
        mode="constant",
        truncate=6,
    )
    np.testing.assert_allclose(result.sum() * DX, area.sum(), rtol=2e-12, atol=0)
    return result[WINDOW]


def summarize(path, temperature, fwhm):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", UserWarning)
        optical = LineResonance(db_filepath=path)
    # Read source coordinates separately to avoid copying the full database
    # through the public database property (millions of molecular transitions).
    with xr.open_dataset(path) as ds:
        wavelength = ds.wavelength_nm.values
        strength = optical.line_strengths([temperature])[0]
        area = strength * wavelength**2 / (c * 1e9)
        branch = ds.einstein_a_s.values / ds.upper_total_a_s.values
        w2 = sk.optical.resonance.resonance_polarizability(
            ds.lower_j.values, ds.upper_j.values
        )
        ext = broaden(wavelength, area, fwhm)
        scattering = broaden(wavelength, area * branch, fwhm)
        polarized = broaden(wavelength, area * branch * w2, fwhm)
        in_band = (wavelength >= 274) & (wavelength <= 810)
        strongest = np.argmax(np.where(in_band, area, -np.inf))
        minimum = float(ds.attrs.get("wavelength_min_nm", 274))
        maximum = float(ds.attrs.get("wavelength_max_nm", 810))
        valid = (minimum <= WAVELENGTH) & (maximum >= WAVELENGTH)
        metadata = {
            "species": path.stem,
            "kind": path.parent.name,
            "lines_in_file": ds.sizes["line"],
            "isotopologue": ds.attrs.get("isotopologue", "atomic isotope mixture"),
            "source_wavelength_min_nm": minimum,
            "source_wavelength_max_nm": maximum,
            "strongest_lte_line_nm": float(wavelength[strongest]),
            "strongest_line_w2": float(w2[strongest]),
            "strongest_line_elastic_fraction": float(branch[strongest]),
            "display_peak_extinction_m2": float(np.max(ext[valid])),
            "band_line_area_m2_nm": float(area[in_band].sum()),
            "incomplete_decay_warning": any(
                "decay sums" in str(w.message) for w in caught
            ),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    omega = np.divide(scattering, ext, out=np.zeros_like(ext), where=ext > 0)
    effective_w2 = np.divide(
        polarized, scattering, out=np.zeros_like(ext), where=scattering > 0
    )
    for values in (ext, scattering, omega, effective_w2):
        values[~valid] = np.nan
    np.testing.assert_array_less(omega[valid], 1 + 1e-12)
    np.testing.assert_array_less(effective_w2[valid], 1 + 1e-12)
    del optical, strength, area, branch, w2
    gc.collect()
    return {
        "ext": ext,
        "scattering": scattering,
        "omega": omega,
        "w2": effective_w2,
        "metadata": metadata,
    }


def emission_data(root):
    result = {}
    for species in ["FeO", "NiO"]:
        path = root / "emission" / f"{species}.nc"
        with xr.open_dataset(path) as ds:
            wavelength = ds.wavelength_nm.values
            spectrum = ds.photon_spectrum.values
            np.testing.assert_allclose(trapezoid(spectrum, wavelength), 1, rtol=1e-10)
            result[species] = {
                "wavelength": wavelength,
                "spectrum": spectrum,
                "metadata": dict(ds.attrs),
            }
    return result


def save(fig, output, stem, title, caption, gallery):
    for extension in ("png", "svg"):
        fig.savefig(output / f"{stem}.{extension}", dpi=175)
    plt.close(fig)
    gallery.append({"stem": stem, "title": title, "caption": caption})
    print(f"Saved {stem}", flush=True)


def phase_curves(ax, *, polarization=False):
    angle = np.linspace(0, 180, 721)
    mu = np.cos(np.deg2rad(angle))
    cases = [
        (0, "D1: W₂ = 0"),
        (0.1, "High-J P/R: W₂ ≈ 0.1"),
        (0.4, "High-J Q: W₂ ≈ 0.4"),
        (0.5, "D2: W₂ = 0.5"),
        (1, "J = 0 → 1: W₂ = 1"),
    ]
    for (w2, name), color in zip(cases, COLORS, strict=True):
        values = resonance_phase_function(mu, w2)
        if polarization:
            values = 3 * w2 * (1 - mu**2) / (4 - w2 + 3 * w2 * mu**2)
        ax.plot(angle, values, color=color, lw=2, label=name)
    ax.set_xlim(0, 180)
    ax.set_xticks(np.arange(0, 181, 30))
    ax.set_xlabel("Scattering angle [degrees]")
    ax.set_ylabel("Single-scatter |Q| / I" if polarization else "Phase function P₁₁")
    ax.set_ylim((0, 1.05) if polarization else (0.7, 1.55))
    decorate(ax, spectral=False)


def emission_curves(ax, emission):
    for species, color in [("FeO", ORANGE), ("NiO", BLUE)]:
        item = emission[species]
        ax.plot(
            item["wavelength"],
            item["spectrum"],
            color=color,
            lw=1.8,
            label=(
                "FeO · numerical model"
                if species == "FeO"
                else "NiO · provisional, 5 nm trace"
            ),
        )
        if species == "NiO":
            ax.scatter(item["wavelength"], item["spectrum"], color=color, s=9, zorder=3)
    ax.set_xlim(420, 730)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Published wavelength [nm; medium unspecified]")
    ax.set_ylabel("Normalized photon spectrum [nm⁻¹]")
    ax.legend(loc="upper right")
    decorate(ax, spectral=False)


def overview(profiles, emission, output, temperature, fwhm, gallery):
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.subplots_adjust(
        left=0.09, right=0.965, bottom=0.115, top=0.867, hspace=0.42, wspace=0.25
    )
    heading(
        fig,
        "Metal spectroscopy · the overall picture",
        f"Local database overview  |  LTE cross sections at {temperature:g} K  |  Ordinary mesospheric candidates",
    )
    ax = axes[0, 0]
    species = ["Mg_II", "Mg_I", "Ni_I", "Fe_I", "Ca_II", "Ca_I", "Na_I", "Li_I", "K_I"]
    for row, key in enumerate(species):
        data = profiles[key]["ext"]
        curve = data / np.nanmax(data)
        ax.fill_between(
            WAVELENGTH,
            row + 0.8,
            row + 0.8 - 0.77 * curve,
            color=COLORS[row % len(COLORS)],
            alpha=0.9,
        )
    ax.set_yticks(np.arange(len(species)) + 0.25, [label(s) for s in species])
    ax.set_ylim(-0.15, len(species))
    ax.invert_yaxis()
    decorate(ax)
    ax.grid(axis="y", visible=False)
    ax.set_title("A  Selected atomic signatures", loc="left", pad=12)
    ax.text(
        0,
        1.01,
        f"Each row has its own peak normalization · {fwhm:g} nm display blur",
        transform=ax.transAxes,
        fontsize=9,
        color=MUTED,
    )
    ax = axes[0, 1]
    for species, color in zip(MOLECULES, COLORS, strict=False):
        y = profiles[species]["ext"]
        ax.plot(
            WAVELENGTH,
            np.where(y > 1e-28, y, np.nan),
            color=color,
            lw=1.5,
            label=species,
        )
    ax.set_yscale("log")
    ax.set_ylim(1e-27, 3e-18)
    ax.set_ylabel("Extinction cross section [m² / molecule]")
    ax.legend(ncol=2, loc="upper right")
    ax.set_title("B  Molecular absorption bands", loc="left", pad=12)
    decorate(ax)
    ax = axes[1, 0]
    phase_curves(ax)
    ax.legend(ncol=1, fontsize=8.5, loc="upper center")
    ax.set_title("C  Resonance phase-function family", loc="left", pad=12)
    ax = axes[1, 1]
    emission_curves(ax, emission)
    ax.set_title("D  Free-VER emission templates", loc="left", pad=12)
    footer(
        fig,
        f"A–B: {fwhm:g} nm FWHM Gaussian display broadening of integrated line areas; not a radiance or detectability calculation.\nC: isolated E1 lines; no hyperfine or magnetic effects. D: unit photon integral over each stored band; no abundance or chemical yield assumed.",
    )
    save(
        fig,
        output,
        "00_summary",
        "The overall picture",
        "Atomic line locations, molecular absorption, resonance angular shapes, and independently fitted chemical emission.",
        gallery,
    )


def inventory(profiles, output, temperature, fwhm, gallery):
    atomic = [key for key in PRIORITY if key in profiles]
    atomic += sorted(
        key for key in profiles if key not in atomic and key not in MOLECULES
    )
    order = atomic + MOLECULES
    relative = np.array(
        [profiles[key]["ext"] / np.nanmax(profiles[key]["ext"]) for key in order]
    )
    shown = np.ma.masked_where(
        ~np.isfinite(relative) | (relative < 1e-4),
        np.log10(np.maximum(relative, 1e-300)),
    )
    fig, ax = plt.subplots(figsize=(16, 17))
    fig.subplots_adjust(left=0.125, right=0.90, bottom=0.13, top=0.905)
    heading(
        fig,
        "Wavelength inventory · every installed species",
        f"{len(atomic)} atomic/ionic spectra + 4 molecular line lists  |  {temperature:g} K  |  {fwhm:g} nm display blur",
    )
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#f1f4f7")
    result = ax.imshow(
        shown,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        norm=Normalize(-4, 0),
        extent=[274, 810, len(order) - 0.5, -0.5],
        rasterized=True,
    )
    ax.set_yticks(np.arange(len(order)), [label(s) for s in order], fontsize=10)
    ax.set_xlabel("Vacuum wavelength [nm]", labelpad=10)
    ax.set_xticks([274, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 810])
    ax.tick_params(axis="y", length=0)
    for row in np.arange(len(order)) - 0.5:
        ax.axhline(row, color="white", lw=0.6)
    ax.axhline(len(atomic) - 0.5, color=INK, lw=1.5)
    for row, species in enumerate(order):
        ax.text(
            1.015,
            row,
            f"{profiles[species]['metadata']['lines_in_file']:,}",
            transform=ax.get_yaxis_transform(),
            fontsize=8.5,
            color=MUTED,
            va="center",
        )
    ax.text(
        1.015, 1.015, "Lines in file", transform=ax.transAxes, fontsize=9, color=MUTED
    )
    colorbar = fig.colorbar(
        result, ax=ax, orientation="horizontal", pad=0.065, fraction=0.025, aspect=55
    )
    colorbar.set_ticks(
        [-4, -3, -2, -1, 0], labels=["10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "1"]
    )
    colorbar.set_label(
        "Cross section relative to each species' own maximum; rows do not compare absolute strength"
    )
    footer(
        fig,
        "All cached species are shown, including auxiliary atomic-data candidates with no established ambient detectability.\nPale cells: <10⁻⁴ of that row's maximum or unavailable coverage. AlO coverage starts at 285.7 nm. Emission templates appear separately.",
    )
    save(
        fig,
        output,
        "01_all_species",
        "Complete wavelength inventory",
        "All 51 cross-section datasets. Each row is normalized independently; faint species are not made more detectable by this normalization.",
        gallery,
    )


def atomic_panels(profiles, output, temperature, fwhm, gallery):
    fig, axes = plt.subplots(6, 4, figsize=(18, 16), sharex=True, sharey=True)
    fig.subplots_adjust(
        left=0.075, right=0.97, bottom=0.105, top=0.89, hspace=0.49, wspace=0.18
    )
    heading(
        fig,
        "Atomic cross sections · a common absolute scale",
        f"22 selected candidates  |  {temperature:g} K  |  {fwhm:g} nm FWHM Gaussian display broadening",
    )
    for ax, species in zip(axes.flat, PRIORITY, strict=False):
        item = profiles[species]
        for key, color, linestyle in [("ext", BLUE, "-"), ("scattering", ORANGE, "--")]:
            ax.plot(
                WAVELENGTH,
                np.where(item[key] > 1e-27, item[key], np.nan),
                color=color,
                ls=linestyle,
                lw=1.2,
            )
        ax.set_yscale("log")
        ax.set_ylim(1e-26, 3e-17)
        ax.set_yticks([1e-26, 1e-23, 1e-20, 1e-17])
        ax.set_title(label(species), loc="left", fontsize=12, pad=8)
        meta = item["metadata"]
        if meta["display_peak_extinction_m2"] < 1e-26:
            ax.text(
                0.5,
                0.47,
                f"Peak {meta['display_peak_extinction_m2']:.1e} m²\n(below the shared axis range)",
                transform=ax.transAxes,
                ha="center",
                fontsize=9,
                color=MUTED,
            )
        ax.text(
            1,
            1.045,
            f"strongest {meta['strongest_lte_line_nm']:.2f} nm",
            transform=ax.transAxes,
            ha="right",
            fontsize=8.5,
            color=MUTED,
        )
        decorate(ax)
        ax.set_xticks([300, 500, 700])
        ax.set_xlabel("")
        ax.tick_params(labelsize=9)
    for ax in axes.flat[len(PRIORITY) :]:
        ax.axis("off")
    key = axes.flat[-2]
    key.plot([], [], color=BLUE, label="Extinction")
    key.plot([], [], color=ORANGE, ls="--", label="Same-line elastic return")
    key.legend(loc="upper left", fontsize=11)
    key.text(
        0.02,
        0.3,
        "Other radiative branches\nleave the incident wavelength.\nNo solar spectrum or density applied.",
        transform=key.transAxes,
        fontsize=10,
        color=MUTED,
        va="top",
    )
    fig.text(
        0.025,
        0.5,
        "Cross section [m² / atom or ion]",
        rotation=90,
        va="center",
        fontsize=13,
    )
    fig.text(0.50, 0.069, "Vacuum wavelength [nm]", ha="center", fontsize=12)
    footer(
        fig,
        "Solid and dashed curves overlap for nearly closed transitions. Available atomic decay sums can overestimate elastic return.\nCommon vertical scale; portions below 10⁻²⁶ m² are not shown. Actual resolved line peaks are much higher than these display-broadened values.",
    )
    save(
        fig,
        output,
        "02_atomic_cross_sections",
        "Absolute atomic cross sections",
        "Extinction and elastic return on one shared scale, so weak and strong per-atom signals can be compared without assuming abundances.",
        gallery,
    )


def molecular_panels(profiles, output, temperature, fwhm, gallery):
    fig, axes = plt.subplots(4, 2, figsize=(16, 13), sharex=True)
    fig.subplots_adjust(
        left=0.10, right=0.965, bottom=0.12, top=0.875, hspace=0.50, wspace=0.22
    )
    heading(
        fig,
        "Molecular absorption, branching and angular response",
        f"ExoMol line lists  |  {temperature:g} K LTE  |  {fwhm:g} nm display blur  |  One named isotopologue per molecule",
    )
    for row, species in enumerate(MOLECULES):
        item = profiles[species]
        ax = axes[row, 0]
        ax.plot(
            WAVELENGTH,
            np.where(item["ext"] > 1e-28, item["ext"], np.nan),
            color=BLUE,
            lw=1.5,
            label="Extinction",
        )
        ax.plot(
            WAVELENGTH,
            np.where(item["scattering"] > 1e-28, item["scattering"], np.nan),
            color=ORANGE,
            lw=1.4,
            ls="--",
            label="Elastic return",
        )
        ax.set_yscale("log")
        ax.set_ylim(1e-27, 3e-18)
        ax.set_ylabel("Cross section [m²]")
        ax.set_title(
            f"{species}  ·  {item['metadata']['isotopologue']}", loc="left", fontsize=12
        )
        decorate(ax)
        if row == 0:
            ax.legend(loc="upper right", ncol=2, fontsize=9)
        ax = axes[row, 1]
        bright = item["ext"] > np.nanmax(item["ext"]) * 1e-3
        ax.plot(
            WAVELENGTH,
            np.where(bright, item["omega"], np.nan),
            color=ORANGE,
            label="Elastic-return fraction",
            lw=1.3,
        )
        ax.plot(
            WAVELENGTH,
            np.where(bright, item["w2"], np.nan),
            color=TEAL,
            label="Effective W₂",
            lw=1.3,
        )
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_title(
            f"{species}  ·  ratios of smoothed line contributions",
            loc="left",
            fontsize=12,
        )
        decorate(ax)
        if row == 0:
            ax.legend(loc="upper right", fontsize=9)
        if row < 3:
            axes[row, 0].set_xlabel("")
            axes[row, 1].set_xlabel("")
    footer(
        fig,
        "Ratios shown only where extinction exceeds 10⁻³ of that molecule's maximum. W₂ is weighted by same-line elastic scattering.\nThese LTE absorption/elastic-return models do not include fluorescence between different wavelengths or chemical emission. AlO has no coverage below 285.7 nm.",
    )
    save(
        fig,
        output,
        "03_molecular_cross_sections",
        "Molecular bands and elastic branching",
        "Cross sections plus spectral changes in the same-line return fraction and scattering-weighted polarizability.",
        gallery,
    )


def phases(output, gallery):
    fig, axes = plt.subplots(1, 2, figsize=(16, 7.7))
    fig.subplots_adjust(left=0.085, right=0.965, bottom=0.27, top=0.84, wspace=0.22)
    heading(
        fig,
        "What the resonance phase functions look like",
        "The angular shape is set by each transition's W₂; it is independent of its absolute line strength",
    )
    for ax, polarization in zip(axes, [False, True], strict=True):
        phase_curves(ax, polarization=polarization)
    axes[0].set_title("Intensity redistribution", loc="left")
    axes[1].set_title("Polarization of initially unpolarized light", loc="left")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.52, 0.13),
        ncol=3,
        fontsize=10.5,
    )
    footer(
        fig,
        "P₁₁ has unit spherical mean (integral over solid angle = 4π). The right panel is single scattering in the scattering-plane basis.\nIsolated electric-dipole approximation: no hyperfine/isotope splitting, magnetic fields, lower-state alignment or line interference.\nChemical VER sources are assumed isotropic and unpolarized at emission; their amplitude does not require a resonance phase function.",
    )
    save(
        fig,
        output,
        "04_phase_functions",
        "Resonance phase and polarization",
        "The whole isolated-line phase family, including D1, D2, classical Rayleigh and high-J molecular limits.",
        gallery,
    )


def resolved(root, output, gallery):
    targets = [
        ("Na_I", 589.1583),
        ("K_I", 770.1084),
        ("Li_I", 670.961),
        ("Mg_I", 285.2964),
        ("Ca_II", 393.4777),
        ("Fe_I", 372.0993),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    fig.subplots_adjust(
        left=0.075, right=0.97, bottom=0.17, top=0.86, hspace=0.43, wspace=0.28
    )
    heading(
        fig,
        "Resolved line profiles · temperature and intrinsic width",
        "Voigt profiles from the runtime optical properties  |  No display or instrument smoothing",
    )
    for ax, (species, target) in zip(axes.flat, targets, strict=True):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            optical = LineResonance(
                root / "atomic" / f"{species}.nc", line_wing_cutoff_nm=None
            )
        dataset = optical.database
        line = int(np.argmin(np.abs(dataset.wavelength_nm.values - target)))
        center = float(dataset.wavelength_nm.values[line])
        sigma_nm = (
            center * np.sqrt(k * 300 / (dataset.attrs["mass_amu"] * atomic_mass)) / c
        )
        offset = np.linspace(-6 * sigma_nm, 6 * sigma_nm, 1801)
        result = optical.cross_sections(
            center + offset,
            np.array([0, 1, 2]),
            temperature_k=np.array([150, 200, 300]),
        )
        for index, (temperature, color) in enumerate(
            zip([150, 200, 300], [BLUE, TEAL, ORANGE], strict=True)
        ):
            ax.plot(
                offset * 1000,
                result.extinction[index],
                color=color,
                lw=1.9,
                label=f"{temperature} K",
            )
        branch = float(
            dataset.einstein_a_s.values[line] / dataset.upper_total_a_s.values[line]
        )
        if branch < 0.995:
            ax.plot(
                offset * 1000,
                result.extinction[1] * result.ssa[1],
                color=INK,
                ls="--",
                lw=1.3,
                label="Elastic, 200 K",
            )
        w2 = float(
            sk.optical.resonance.resonance_polarizability(
                dataset.lower_j.values[line], dataset.upper_j.values[line]
            )
        )
        ax.set_title(f"{label(species)}  ·  {center:.5f} nm", loc="left", fontsize=12)
        ax.text(
            0.97,
            0.93,
            f"W₂ = {w2:.3g}\nElastic fraction ≈ {branch:.3f}",
            transform=ax.transAxes,
            fontsize=9,
            ha="right",
            va="top",
            color=MUTED,
        )
        ax.set_xlabel("Wavelength offset from line center [pm]")
        ax.set_ylabel("Extinction cross section [m²]")
        ax.set_ylim(bottom=0)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        decorate(ax, spectral=False)
    handles, labels = axes[1, 1].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="lower center", bbox_to_anchor=(0.52, 0.073), ncol=4
    )
    footer(
        fig,
        "Each panel has its own vertical scale. Widths combine thermal Doppler and available natural decay rates; no pressure broadening.\nVacuum line centers; the atomic model resolves fine-structure transitions, not hyperfine/isotope components. Incomplete decay sums limit branching accuracy.",
    )
    save(
        fig,
        output,
        "05_resolved_line_profiles",
        "Intrinsic line profiles at 150, 200 and 300 K",
        "Resolved cross sections show the very narrow structure hidden by the broad-band overview, and how it changes with temperature.",
        gallery,
    )


def emissions(emission, output, gallery):
    fig, axes = plt.subplots(
        2, 1, figsize=(15, 10), sharex=True, gridspec_kw={"height_ratios": [1.3, 1]}
    )
    fig.subplots_adjust(left=0.10, right=0.965, bottom=0.16, top=0.86, hspace=0.27)
    heading(
        fig,
        "Chemical emission · the actual templates used for free VER",
        "Both integrate to one photon over their own stored wavelength interval; the fitted VER sets the amplitude",
    )
    emission_curves(axes[0], emission)
    axes[0].set_xlabel("")
    for species, color in [("FeO", ORANGE), ("NiO", BLUE)]:
        item = emission[species]
        cumulative = cumulative_trapezoid(
            item["spectrum"], item["wavelength"], initial=0
        )
        axes[1].plot(item["wavelength"], cumulative, color=color, lw=2, label=species)
        axes[1].axvline(item["wavelength"][0], color=color, ls=":", lw=1, alpha=0.7)
        axes[1].axvline(item["wavelength"][-1], color=color, ls=":", lw=1, alpha=0.7)
    axes[1].set_ylim(0, 1.04)
    axes[1].set_ylabel("Cumulative fraction of band photons")
    axes[1].set_xlabel("Published wavelength [nm; medium unspecified]")
    axes[1].legend(loc="upper left")
    decorate(axes[1], spectral=False)
    footer(
        fig,
        "FeO: published numerical model, 560–719.9 nm, 0.1 nm sampling (sampling is not accuracy). NiO: approximate 5 nm digitization, 430–670 nm.\nNiO attribution is provisional. Neither spectrum has a specified air/vacuum convention. No additional smoothing is applied here.\nVER refers to the stored band; a fit to a smaller wavelength interval retains the same normalization. No chemistry, abundance or full-band yield is inferred.",
    )
    save(
        fig,
        output,
        "06_emission_templates",
        "FeO and provisional NiO VER templates",
        "Published photon shapes and their cumulative band fractions. Source uncertainties and different wavelength supports are explicit.",
        gallery,
    )


def write_gallery(output, gallery, temperature, fwhm):
    cards = []
    for item in gallery:
        stem = item["stem"]
        cards.append(
            f"<article><h2>{html.escape(item['title'])}</h2><p>{html.escape(item['caption'])}</p>"
            f'<a href="{stem}.png"><img src="{stem}.png" alt="{html.escape(item["title"])}"></a>'
            f'<p class="links"><a href="{stem}.png">PNG</a><a href="{stem}.svg">Vector SVG</a></p></article>'
        )
    document = f"""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Metal spectroscopy summary plots</title><style>
body{{font:16px/1.65 system-ui,sans-serif;color:#172b40;background:#edf2f6;margin:0}}main{{max-width:1250px;margin:auto;padding:50px 24px}}
h1{{font-size:38px;line-height:1.15;margin-bottom:16px}}h2{{font-size:23px;margin-bottom:4px}}p{{max-width:950px;color:#566779}}
article{{background:white;border:1px solid #dce5ed;border-radius:12px;padding:24px;margin:28px 0}}img{{width:100%;height:auto}}a{{color:#17689c}}
.links{{display:flex;gap:24px;margin:8px 0 0}}.badge{{font-size:13px;letter-spacing:.08em;text-transform:uppercase;color:#078177}}
</style><main><div class="badge">SASKTRAN2 · codex/metal-spectroscopy</div><h1>Metal spectroscopy, plotted</h1>
<p>All 47 atomic/ionic datasets, four molecular line lists, and two emission templates from the local database.
Cross-section overviews use {temperature:g} K LTE and an illustrative {fwhm:g} nm Gaussian blur. Intrinsic line profiles and emission templates have no added blur.
These plots compare spectroscopy; predicted OSIRIS brightness requires species profiles, solar illumination, transfer and instrument response.</p>
<p><a href="plot_metadata.json">Calculation metadata</a> · <a href="spectral_summaries.npz">Numerical spectra (NumPy)</a> · <a href="README.md">Methods and sources</a></p>
{"".join(cards)}</main></html>"""
    (output / "index.html").write_text(document)
    (output / "README.md").write_text(
        f"""# Metal spectroscopy plots

Reproduce: `python tools/spectroscopy/plot_metal_summary.py --output {output}`.
The script reads the normal local SASKTRAN2 database without downloads or changes.

- Cross sections: {temperature:g} K LTE. All profiles are per atom/ion or the named molecular isotopologue, without abundance scaling.
- Display broadening: unit-area Gaussian with {fwhm:g} nm FWHM. Frequency-integrated line strengths are converted to wavelength areas using S_lambda = S_nu lambda_nm² / (c × 10⁹), deposited conservatively on a 0.02 nm grid, then blurred. This narrow-line approximation is only for the broad overview, not radiative transfer or an instrument model.
- Elastic return uses A_ul / Gamma_u. W₂ is weighted by elastic scattering, after the same display broadening. Incomplete atomic decay sums can overestimate return fractions.
- Resolved profiles call the runtime Voigt implementation at 150, 200, 300 K with full wings. Hyperfine/isotope structure is unresolved.
- FeO/NiO curves reproduce the stored photon templates, with no new smoothing. Each has unit area over its own support. Source wavelength medium is unspecified; NiO is a provisional, approximate 5 nm digitization.
- The inventory's colors are normalized independently by species, and are not a detectability ranking. AlO wavelengths below 285.7 nm are unsupported.
- Numerical arrays are in `spectral_summaries.npz`; each species has `_ext`, `_scattering`, `_omega`, `_w2` arrays on `wavelength_nm`. Unsupported coverage is NaN. Emission arrays use `_emission_wavelength_nm` and `_photon_spectrum`. Individual source hashes and peak statistics are in `plot_metadata.json`.

Sources: [NIST ASD](https://physics.nist.gov/asd), [ExoMol](https://exomol.com/data/molecules/),
[FeO numerical supplement](https://acp.copernicus.org/articles/17/4177/2017/),
[NiO model](https://acp.copernicus.org/articles/11/9595/2011/),
[NiO interpretation caveat](https://acp.copernicus.org/articles/24/1143/2024/).
For the full source and license record, see `tools/spectroscopy/README.md`, `MOLECULAR.md`, `EMISSION.md` and the local database provenance.
"""
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path("artifacts/metal_spectroscopy_summary")
    )
    parser.add_argument(
        "--database-root",
        type=Path,
        default=sk.appconfig.database_root() / "spectroscopy/metals",
    )
    parser.add_argument("--temperature", type=float, default=200)
    parser.add_argument("--display-fwhm-nm", type=float, default=1)
    args = parser.parse_args()
    if not 0 < args.display_fwhm_nm <= 3:
        msg = "Display FWHM must be positive and at most 3 nm for this plot grid"
        raise ValueError(msg)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    style()
    profiles = {}
    for kind in ("atomic", "molecular"):
        for path in sorted((args.database_root / kind).glob("*.nc")):
            profiles[path.stem] = summarize(
                path, args.temperature, args.display_fwhm_nm
            )
        print(f"Loaded {kind} line summaries", flush=True)
    emission = emission_data(args.database_root)
    arrays = {"wavelength_nm": WAVELENGTH}
    arrays.update(
        {
            f"{species}_{key}": item[key]
            for species, item in profiles.items()
            for key in ("ext", "scattering", "omega", "w2")
        }
    )
    arrays.update(
        {
            f"{species}_emission_wavelength_nm": item["wavelength"]
            for species, item in emission.items()
        }
    )
    arrays.update(
        {
            f"{species}_photon_spectrum": item["spectrum"]
            for species, item in emission.items()
        }
    )
    np.savez_compressed(output / "spectral_summaries.npz", **arrays)
    metadata = {
        "temperature_k": args.temperature,
        "display_gaussian_fwhm_nm": args.display_fwhm_nm,
        "database_root": str(args.database_root),
        "species": [item["metadata"] for item in profiles.values()],
        "emission": {key: item["metadata"] for key, item in emission.items()},
    }
    (output / "plot_metadata.json").write_text(
        json.dumps(metadata, indent=2, default=lambda value: value.item()) + "\n"
    )
    gallery = []
    overview(
        profiles, emission, output, args.temperature, args.display_fwhm_nm, gallery
    )
    inventory(profiles, output, args.temperature, args.display_fwhm_nm, gallery)
    atomic_panels(profiles, output, args.temperature, args.display_fwhm_nm, gallery)
    molecular_panels(profiles, output, args.temperature, args.display_fwhm_nm, gallery)
    phases(output, gallery)
    resolved(args.database_root, output, gallery)
    emissions(emission, output, gallery)
    write_gallery(output, gallery, args.temperature, args.display_fwhm_nm)
    print(output / "index.html", flush=True)


if __name__ == "__main__":
    main()
