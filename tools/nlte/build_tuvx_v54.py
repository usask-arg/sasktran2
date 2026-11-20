"""Build the TUV-x v5.4 tables used by ``sasktran2.photolysis.TUVActinicFlux``.

Needs the ``musica`` package (TUV-x through MUSICA), which is not a sasktran2
dependency; run this script in its own environment::

    python -m venv musica-venv && musica-venv/bin/pip install musica==0.17.1 xarray netcdf4
    musica-venv/bin/python tools/nlte/build_tuvx_v54.py --output tuvx_v54.nc

The product belongs at ``photolysis/tuvx_v54.nc`` in the sasktran2 standard
database. On the 156 bins of the TUV-x v5.4 wavelength grid, it holds:

- the extraterrestrial photon flux per bin;
- the O3 cross section and O(1D) and O(3P) quantum yields from 180 to 300 K
  in 1 K steps, which covers every temperature knot of the TUV-x O3 data
  (TUV-x interpolates linearly between them and clamps outside 203-295 K);
- the O2 cross section outside the Lyman-alpha and Schumann-Runge band bins;
- the Rayleigh cross section of air;
- the Koppers and Murtagh (1996) Schumann-Runge band coefficients.

Cross sections and yields are taken from TUV-x itself, from the diagnostic
files it writes when ``enable diagnostics`` is set, so they carry TUV-x's own
regridding and temperature interpolation. The band constants are read from
the TUV-x source at the version MUSICA wraps.

TUV-x and its data are Apache-2.0 licensed (Copyright (C) 2020 National
Center for Atmospheric Research); see https://github.com/NCAR/tuv-x.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import musica
import musica.backend
import musica.tuvx.v54 as v54
import numpy as np
import xarray as xr
from musica.tuvx.grid_map import GridMap
from musica.tuvx.profile import Profile
from musica.tuvx.profile_map import ProfileMap
from musica.tuvx.radiator_map import RadiatorMap
from musica.tuvx.tuvx import TUVX

#: The TUV-x commit wrapped by MUSICA 0.17.1 (TUV-x 0.17.0).
TUVX_COMMIT = "ebec1cb47084fd102e338fb896c6d2b80fbcd292"
LA_SR_BANDS_URL = (
    f"https://raw.githubusercontent.com/NCAR/tuv-x/{TUVX_COMMIT}/src/la_sr_bands.F90"
)
TEMPERATURES_K = np.arange(180.0, 300.1, 1.0)
O3_REACTIONS = {"o1d": "O3+hv->O2+O(1D)", "o3p": "O3+hv->O2+O(3P)"}
O2_REACTION = "O2+hv->O+O"
CM2_TO_M2 = 1.0e-4


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fortran_array(source: str, name: str) -> np.ndarray:
    """Values of the first ``name(...) = (/ ... /)`` parameter array in Fortran source."""
    match = re.search(
        rf"\b{name}\s*\([^)]*\)\s*=\s*&?\s*\(/(.*?)/\)", source, flags=re.S | re.I
    )
    if match is None:
        msg = f"{name} not found in la_sr_bands.F90"
        raise ValueError(msg)
    body = re.sub(r"!.*", "", match.group(1)).replace("&", " ")
    return np.array(
        [float(v.replace("_dk", "").replace("D", "E")) for v in body.split(",")]
    )


def fortran_scalar(source: str, name: str) -> float:
    match = re.search(rf"\b{name}\s*=\s*([0-9.eEdD+-]+)(_dk)?", source)
    if match is None:
        msg = f"{name} not found in la_sr_bands.F90"
        raise ValueError(msg)
    return float(match.group(1).replace("D", "E"))


def read_band_constants(source_dir: Path) -> dict:
    path = source_dir / "la_sr_bands.F90"
    if not path.exists():
        urllib.request.urlretrieve(LA_SR_BANDS_URL, path)
    source = path.read_text()
    return {
        "source_sha256": sha256(path),
        "lyman_alpha_edges_nm": fortran_array(source, "wlla"),
        "srb_edges_nm": fortran_array(source, "wlsrb"),
        "default_cross_section_cm2": fortran_array(source, "xslod"),
        "num_coefficients": int(fortran_scalar(source, "nPoly")),
        "log_column_limits": (
            fortran_scalar(source, "kLowerLimit"),
            fortran_scalar(source, "kUpperLimit"),
        ),
        "reference_temperature_k": fortran_scalar(source, "T0"),
    }


def read_chebyshev_coefficients(path: Path, num_coefficients: int) -> tuple:
    """The A and B coefficient blocks of O2_parameters.txt, each (coefficient, band)."""
    lines = path.read_text().splitlines()
    blocks = []
    for start in (2, 4 + num_coefficients):
        rows = lines[start : start + num_coefficients]
        blocks.append(
            np.array([[float(v) for v in row.split(",") if v.strip()] for row in rows])
        )
    return tuple(blocks)


def read_diagnostic(path: Path, num_bins: int, num_levels: int) -> np.ndarray:
    """A TUV-x (level, wavelength) double array from an unformatted record."""
    raw = path.read_bytes()
    size = int(np.frombuffer(raw[:4], dtype=np.int32)[0])
    values = np.frombuffer(raw[4 : 4 + size], dtype=np.float64)
    # Column-major (level, wavelength) is row-major (wavelength, level).
    return values.reshape(num_bins, num_levels).T


def run_diagnostics(workdir: Path) -> tuple:
    """Run TUV-x v5.4 once with level temperatures TEMPERATURES_K."""
    config_path = Path(v54.config_file_path())
    config = json.loads(config_path.read_text())
    config["photolysis"]["enable diagnostics"] = True
    (workdir / "output").mkdir()
    os.symlink(config_path.parent / "data", workdir / "data")
    diagnostic_config = workdir / "tuv_5_4_diagnostics.json"
    diagnostic_config.write_text(json.dumps(config))

    grids = GridMap()
    grids["height", "km"] = v54.height_grid()
    grids["wavelength", "nm"] = v54.wavelength_grid()
    heights = grids["height", "km"]
    wavelengths = grids["wavelength", "nm"]
    if len(heights.edges) != TEMPERATURES_K.size:
        msg = "the v5.4 height grid no longer has one level per temperature"
        raise ValueError(msg)
    profiles = ProfileMap()
    for name in ("air", "O2", "O3"):
        profiles[name, "molecule cm-3"] = v54.profile(name, heights)
    profiles["temperature", "K"] = Profile(
        name="temperature",
        units="K",
        grid=heights,
        edge_values=TEMPERATURES_K,
        midpoint_values=0.5 * (TEMPERATURES_K[1:] + TEMPERATURES_K[:-1]),
    )
    profiles["surface albedo", "none"] = v54.profile("surface albedo", wavelengths)
    solar = v54.profile("extraterrestrial flux", wavelengths)
    profiles["extraterrestrial flux", "photon cm-2 s-1"] = solar
    radiators = RadiatorMap()
    radiators["aerosol"] = v54.radiator("aerosol", heights, wavelengths)

    tuvx = TUVX(
        grid_map=grids,
        profile_map=profiles,
        radiator_map=radiators,
        config_path=str(diagnostic_config),
    )
    cwd = Path.cwd()
    os.chdir(workdir)
    try:
        tuvx.run(sza=np.radians(30.0), earth_sun_distance=1.0)
    finally:
        os.chdir(cwd)

    # MUSICA arrays are views into memory owned by the map objects: keep the
    # maps alive and copy.
    grid_map, profile_map = tuvx.get_grid_map(), tuvx.get_profile_map()
    edges = np.array(grid_map["wavelength", "nm"].edges, dtype=float, copy=True)
    solar_per_bin = np.array(
        profile_map["extraterrestrial flux", "photon cm-2 s-1"].midpoint_values,
        dtype=float,
        copy=True,
    )
    num_bins = edges.size - 1

    def diagnostic(reaction, kind):
        return read_diagnostic(
            workdir / "output" / f"{reaction}.{kind}.new", num_bins, TEMPERATURES_K.size
        )

    return edges, solar_per_bin, diagnostic, config_path.parent / "data"


def main(output: Path, source_dir: Path) -> None:
    source_dir.mkdir(parents=True, exist_ok=True)
    band = read_band_constants(source_dir)

    with tempfile.TemporaryDirectory() as tmp:
        edges, solar_per_bin, diagnostic, data_dir = run_diagnostics(Path(tmp))
        o3_cross_section = diagnostic(O3_REACTIONS["o1d"], "xsect")
        if not np.allclose(diagnostic(O3_REACTIONS["o3p"], "xsect"), o3_cross_section):
            msg = "the O3 channels have different cross sections"
            raise ValueError(msg)
        yields = {k: diagnostic(r, "qyld") for k, r in O3_REACTIONS.items()}
        o2_cross_section = diagnostic(O2_REACTION, "xsect")

    a_coefficients, b_coefficients = read_chebyshev_coefficients(
        data_dir / "cross_sections" / "O2_parameters.txt", band["num_coefficients"]
    )

    centres = 0.5 * (edges[:-1] + edges[1:])
    lyman_alpha = (edges[:-1] >= band["lyman_alpha_edges_nm"][0] - 1e-6) & (
        edges[1:] <= band["lyman_alpha_edges_nm"][-1] + 1e-6
    )
    srb = (edges[:-1] >= band["srb_edges_nm"][0] - 1e-6) & (
        edges[1:] <= band["srb_edges_nm"][-1] + 1e-6
    )
    if lyman_alpha.sum() != 1 or srb.sum() != band["srb_edges_nm"].size - 1:
        msg = "the band edges of la_sr_bands.F90 do not match the v5.4 grid"
        raise ValueError(msg)
    parameterised = lyman_alpha | srb
    # Elsewhere the O2 cross section is the same at every level.
    if not np.allclose(
        o2_cross_section[:, ~parameterised], o2_cross_section[0, ~parameterised]
    ):
        msg = "the O2 cross section depends on level outside the parameterised bins"
        raise ValueError(msg)
    o2 = np.where(parameterised, np.nan, o2_cross_section[0])

    # Rayleigh cross section of air at bin centres (TUV-x rayliegh.F90,
    # after Nicolet 1984), lambda in micrometres.
    wavelength_um = 1.0e-3 * centres
    power = np.where(
        wavelength_um <= 0.55,
        3.6772 + 0.389 * wavelength_um + 0.09426 / wavelength_um,
        4.04,
    )
    rayleigh = 4.02e-28 / wavelength_um**power

    sources = {
        "la_sr_bands.F90": band["source_sha256"],
        **{
            name: sha256(data_dir / "cross_sections" / name)
            for name in (
                "O2_1.nc",
                "O2_parameters.txt",
                "O3_1.nc",
                "O3_2.nc",
                "O3_3.nc",
                "O3_4.nc",
            )
        },
        "tuv_5_4.json": sha256(Path(v54.config_file_path())),
    }
    ds = xr.Dataset(
        {
            "solar_flux": (("wavelength",), solar_per_bin / CM2_TO_M2),
            "o2_cross_section": (("wavelength",), o2 * CM2_TO_M2),
            "o3_cross_section": (
                ("temperature_k", "wavelength"),
                o3_cross_section * CM2_TO_M2,
            ),
            "o3_o1d_quantum_yield": (("temperature_k", "wavelength"), yields["o1d"]),
            "o3_o3p_quantum_yield": (("temperature_k", "wavelength"), yields["o3p"]),
            "rayleigh_cross_section": (("wavelength",), rayleigh * CM2_TO_M2),
            "lyman_alpha_bin": (("wavelength",), lyman_alpha),
            "schumann_runge_bin": (("wavelength",), srb),
            "srb_chebyshev_a": (("chebyshev", "srb_band"), a_coefficients),
            "srb_chebyshev_b": (("chebyshev", "srb_band"), b_coefficients),
            "srb_default_cross_section": (
                ("srb_band",),
                band["default_cross_section_cm2"] * CM2_TO_M2,
            ),
        },
        coords={
            "wavelength": centres,
            "wavelength_edge": ("wavelength_edge", edges),
            "temperature_k": TEMPERATURES_K,
        },
        attrs={
            "title": "TUV-x v5.4 photolysis tables for sasktran2.photolysis",
            "source": (
                "TUV-x {tuvx} through MUSICA {musica}, v5.4 configuration "
                "(tuv_5_4.json); tools/nlte/build_tuvx_v54.py"
            ).format(
                tuvx=musica.backend.get_backend()._tuvx._get_tuvx_version(),
                musica=musica.__version__,
            ),
            "tuvx_commit": TUVX_COMMIT,
            "source_sha256": json.dumps(sources),
            "license": "Apache-2.0",
            "attribution": (
                "TUV-x, Copyright (C) 2020 National Center for Atmospheric "
                "Research, Apache License 2.0, https://github.com/NCAR/tuv-x"
            ),
            "references": (
                "Koppers, G. A. A., and D. P. Murtagh (1996), Model studies of the "
                "influence of O2 photodissociation parameterizations in the "
                "Schumann-Runge bands on ozone related photolysis in the upper "
                "atmosphere, Ann. Geophys., 14, 68-79. "
                "Nicolet, M. (1984), On the molecular scattering in the terrestrial "
                "atmosphere: an empirical formula for its calculation in the "
                "homosphere, Planet. Space Sci., 32, 1467-1468."
            ),
            "srb_log_column_limits": np.array(band["log_column_limits"]),
            "srb_reference_temperature_k": band["reference_temperature_k"],
            "history": f"{datetime.now(timezone.utc).isoformat()} built",
        },
    )
    ds["solar_flux"].attrs.update(
        units="photons m^-2 s^-1",
        long_name="extraterrestrial photon flux per bin at 1 AU",
    )
    for name in ("o2_cross_section", "o3_cross_section", "rayleigh_cross_section"):
        ds[name].attrs["units"] = "m^2"
    ds["o2_cross_section"].attrs[
        "comment"
    ] = "NaN in the Lyman-alpha and Schumann-Runge band bins, which are parameterised"
    ds["srb_default_cross_section"].attrs["units"] = "m^2"
    ds["srb_chebyshev_a"].attrs["comment"] = (
        "sigma = exp(A (T - T0) + B) [cm^2] with A and B Chebyshev series in ln N "
        "[ln cm^-2] over srb_log_column_limits"
    )
    ds["wavelength"].attrs["units"] = "nm"
    ds["wavelength_edge"].attrs["units"] = "nm"
    ds.to_netcdf(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=Path("tuvx_v54.nc"))
    parser.add_argument("--source-dir", type=Path, default=Path("tuvx_sources"))
    args = parser.parse_args()
    main(args.output, args.source_dir)
