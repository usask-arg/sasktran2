"""Build band-limited FeO and NiO photon-emission templates for free VER fits.

FeO uses a published numerical model spectrum; NiO uses an explicitly approximate
digitization of a published 5 nm-averaged mesospheric model. Neither is an
absorption cross section, an excitation model, or a conversion to metal density.

Run with the project environment, e.g.::

    python tools/spectroscopy/build_metal_emission_templates.py --output <directory>

The directory should be ``database/spectroscopy/metals/emission``. Source files,
unmodified numerical tables, digitization coordinates, and a NiO overlay are
retained below it. Downloads are pinned by SHA-256. The source wavelength media
are not specified in the publications, so the coordinates are explicitly marked
``unspecified``; no air-to-vacuum conversion is silently assumed.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import logging
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

FE_SOURCE = (
    "https://acp.copernicus.org/articles/17/4177/2017/acp-17-4177-2017-supplement.zip"
)
NI_SOURCE = "https://acp.copernicus.org/articles/11/9595/2011/acp-11-9595-2011.pdf"
FE_HASH = "02f41b4d7ddbc03ba777f8e3c84fffb67cf0656fab8c5872cceb61c3599d495f"
NI_HASH = "d35d0aaf465bdc8d02fe1c3f073b8f8e51b06fac807c8fc0c29533a42bfd34c0"

# Reviewed manual trace, Evans et al. 2011, Fig. 4 (PDF p. 4, second image).
# Original image is 512 x 384 pixels, top-left origin. The SHORT-DASH curve is
# NiO, the dot-dot-dash curve is FeO, long dash is NO2, solid is the total.
# Samples are every 5 nm from 430 through 670 nm. Gaps within dashed strokes
# were interpolated by eye along the identified curve, not other components.
# These are pixel measurements, not claimed original numerical model values.
NI_PIXEL_Y = (
    263,
    261,
    263,
    259,
    249,
    236,
    238,
    231,
    217,
    204,
    211,
    197,
    169,
    206,
    204,
    155,
    185,
    175,
    161,
    162,
    184,
    184,
    172,
    207,
    184,
    164,
    190,
    197,
    199,
    197,
    211,
    216,
    200,
    182,
    186,
    201,
    199,
    168,
    199,
    179,
    201,
    189,
    214,
    190,
    204,
    212,
    207,
    226,
    223,
)


def download(url: str, destination: Path, expected_sha256: str) -> Path:
    """Cache a source, checking that it is the inspected published version."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        temporary = destination.with_suffix(destination.suffix + ".part")
        with urllib.request.urlopen(url, timeout=120) as response:
            temporary.write_bytes(response.read())
        actual = hashlib.sha256(temporary.read_bytes()).hexdigest()
        if actual != expected_sha256:
            temporary.unlink()
            msg = f"Source changed at {url}: inspect it before updating its checksum"
            raise ValueError(msg)
        temporary.replace(destination)
    actual = hashlib.sha256(destination.read_bytes()).hexdigest()
    if actual != expected_sha256:
        msg = f"Cached source checksum mismatch: {destination}"
        raise ValueError(msg)
    return destination


def write_template(
    output: Path,
    species: str,
    wavelength: np.ndarray,
    spectrum: np.ndarray,
    attributes: dict,
) -> Path:
    """Normalize only over the supplied band, using trapezoidal integration."""
    if (
        not np.all(np.isfinite(spectrum))
        or not np.all(np.isfinite(wavelength))
        or np.any(spectrum < 0)
        or np.any(np.diff(wavelength) <= 0)
    ):
        msg = "Invalid wavelength or photon-spectrum values"
        raise ValueError(msg)
    normalization = np.trapezoid(spectrum, wavelength)
    if normalization <= 0:
        msg = "Photon spectrum has zero integral"
        raise ValueError(msg)
    dataset = xr.Dataset(
        {
            "photon_spectrum": (
                "wavelength_nm",
                spectrum / normalization,
                {"units": "nm-1", "long_name": "Normalized relative photon spectrum"},
            ),
        },
        coords={"wavelength_nm": ("wavelength_nm", wavelength, {"units": "nm"})},
        attrs={
            "species": species,
            "schema_version": 1,
            "quantity": "relative photon emission per wavelength",
            "wavelength_medium": "unspecified",
            "normalization": "trapezoidal integral over wavelength_nm equals one",
            "ver_definition": "photons per volume per second integrated over the stored band",
            "wavelength_min_nm": float(wavelength[0]),
            "wavelength_max_nm": float(wavelength[-1]),
            "outside_band": "unsupported; the band-limited template is zero outside its band",
            "emission_angular_model": "isotropic and unpolarized assumed by the VER constituent",
            "chemistry": "none; the VER amplitude is an independent retrieval parameter",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "license": "CC BY 3.0",
            "license_url": "https://creativecommons.org/licenses/by/3.0/",
            **attributes,
        },
    )
    path = output / f"{species}.nc"
    temporary = path.with_suffix(".nc.part")
    dataset.to_netcdf(temporary)
    temporary.replace(path)
    path.with_suffix(".provenance.json").write_text(
        json.dumps(dataset.attrs, indent=2) + "\n"
    )
    return path


def build_feo(output: Path) -> Path:
    """Read the actual Gattinger model column archived with the 2017 paper."""
    raw = output / "raw" / "unterguggenberger2017"
    archive = download(FE_SOURCE, raw / "acp-17-4177-2017-supplement.zip", FE_HASH)
    with zipfile.ZipFile(archive) as source:
        table = source.read("plot_data/fig3.dat")
        (raw / "fig3.dat").write_bytes(table)
        (raw / "source_readme.txt").write_bytes(
            source.read("plot_data/readme_Unterguggenberger_etal_ACP_2017.txt")
        )
    values = np.loadtxt(io.BytesIO(table))
    # The final row (720 nm) contains zero in ALL eleven spectral columns,
    # unlike the nonzero preceding spectrum. Treat it as an unsupported endpoint
    # rather than introducing an artificial sharp fall to zero at 720 nm.
    if not np.all(values[-1, 1:] == 0) or not np.isclose(values[-1, 0], 0.7200):
        msg = "Unexpected FeO endpoint; inspect the source table"
        raise ValueError(msg)
    return write_template(
        output,
        "FeO",
        values[:-1, 0] * 1000,
        values[:-1, 11],
        {
            "source": FE_SOURCE,
            "source_sha256": FE_HASH,
            "source_member": "plot_data/fig3.dat; columns 1 and 12",
            "source_description": "Gattinger et al. (2011) theoretical FeO orange-band spectrum",
            "references": (
                "https://doi.org/10.1139/P11-003; "
                "https://doi.org/10.5194/acp-17-4177-2017; "
                "https://doi.org/10.1029/2010GL045310"
            ),
            "attribution": "Gattinger et al. (2011), as archived by Unterguggenberger et al. (2017)",
            "source_sampling_nm": 0.1,
            "source_resolution": "Not specified for the archived column; 0.1 nm is sampling, not accuracy",
            "source_wavelength_accuracy_nm_approx": 1 / 3,
            "template_quality": "published numerical preliminary spectral model",
            "modifications": "micrometres converted to nm; all-zero final row omitted; band normalization",
            "limitations": (
                "Only the published 560-719.9 nm window is present. The paper reports "
                "significant model/observation shape differences near the main peak. "
                "Upper-state populations are fixed by the source model, not LTE. "
                "The source does not specify air/vacuum wavelength convention. "
                "Do not interpret this amplitude as the full-band photon yield or Fe abundance."
            ),
        },
    )


def build_nio(output: Path) -> Path:
    """Reproduce and retain the inspected approximate NiO digitization."""
    from PIL import ImageDraw  # noqa: PLC0415
    from pypdf import PdfReader  # noqa: PLC0415

    raw = output / "raw" / "evans2011"
    pdf = download(NI_SOURCE, raw / "acp-11-9595-2011.pdf", NI_HASH)
    figure = PdfReader(pdf).pages[3].images[1].image.convert("RGB")
    if figure.size != (512, 384):
        msg = "Source figure dimensions changed"
        raise ValueError(msg)
    figure.save(raw / "figure4.png")
    wavelength = np.arange(430.0, 675.0, 5.0)
    pixels_y = np.asarray(NI_PIXEL_Y, dtype=float)
    pixels_x = 82 + (wavelength - 350) * (480 - 82) / (650 - 350)
    # The axis zero is y=263. Scale cancels under band normalization, but retain
    # the published photon units for transparent digitization reproducibility.
    spectrum = (263 - pixels_y) * (6e7 / (263 - 34))
    np.savetxt(
        raw / "figure4_nio_digitization.csv",
        np.column_stack([wavelength, pixels_x, pixels_y, spectrum]),
        delimiter=",",
        header="wavelength_nm,pixel_x,pixel_y,source_photon_brightness",
        comments="",
    )
    overlay = figure.resize((1536, 1152))
    draw = ImageDraw.Draw(overlay)
    coordinates = [
        (float(x * 3), float(y * 3)) for x, y in zip(pixels_x, pixels_y, strict=True)
    ]
    for x, y in coordinates:
        draw.ellipse((x - 2, y - 2, x + 2, y + 2), fill="red")
    overlay.save(raw / "figure4_nio_digitization_overlay.png")
    return write_template(
        output,
        "NiO",
        wavelength,
        spectrum,
        {
            "source": NI_SOURCE,
            "source_sha256": NI_HASH,
            "source_description": "Figure 4 short-dash NiO mesospheric model, PDF page 4, second image",
            "references": (
                "https://doi.org/10.5194/acp-11-9595-2011; "
                "https://doi.org/10.1139/p11-068; "
                "https://doi.org/10.5194/acp-24-1143-2024"
            ),
            "attribution": "Evans, Gattinger, Broadfoot and Llewellyn (2011), Figure 4",
            "template_quality": "approximate raster digitization of a preliminary model",
            "source_averaging_nm": 5.0,
            "source_sampling_nm": 5.0,
            "digitization_wavelength_pixel_nm": 300 / 398,
            "digitization_intensity_reading_pixels_approx": 2.0,
            "modifications": "manual trace at 5 nm knots; pixel-axis calibration; band normalization",
            "limitations": (
                "Approximate 512x384 raster-figure digitization; cannot support line-resolved "
                "or 1 nm spectroscopy. Source model is already averaged in 5 nm intervals. "
                "The model assumes equal excited vibrational populations for mesopause pressure. "
                "Wavelength medium is unspecified. Only 430-670 nm is supported, not the "
                "total molecular emission. Source interpretation is provisional: Noll et al. "
                "2024 discuss uncertainty and lack of a strong NiO signature in later data. "
                "Blend with NO2 and FeO must be fitted; this template is not a measured abundance."
            ),
        },
    )


def main() -> None:
    """Build the selected species in a supplied local database directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--species", nargs="+", choices=("FeO", "NiO"), default=("FeO", "NiO")
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    args.output.mkdir(parents=True, exist_ok=True)
    for species in args.species:
        path = {"FeO": build_feo, "NiO": build_nio}[species](args.output)
        logging.info("Wrote %s", path)


if __name__ == "__main__":
    main()
