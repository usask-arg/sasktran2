"""Regression fixtures for the photochemistry models.

The fixtures pin the behaviour of the photochemistry code as it was before it
moved into the ``sasktran2-nlte`` crate. They check that refactors do not
change results; they are not a statement of physical truth. Regenerate them
deliberately, after a reviewed behaviour change, with

    python tests/photchem/test_photchem_goldens.py --regenerate
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import xarray as xr
from sasktran2.photchem import LYMAN_ALPHA_WAVELENGTH_NM, Yankovsky

GOLDEN_FILE = Path(__file__).parent / "goldens" / "yankovsky.npz"

# The fixtures were recorded with LAPACK dgesv; the crate uses a pure-Rust LU,
# so agreement is to round-off amplified by the conditioning of the rate
# matrices rather than bit-for-bit.
RTOL = 1.0e-9
ATOL_FRACTION = 1.0e-12


def _photon_flux_toa(wavelength_nm: np.ndarray) -> np.ndarray:
    """5778 K blackbody photon flux at 1 AU, photons m^-2 s^-1 nm^-1."""
    h = 6.62607015e-34
    c = 2.99792458e8
    k = 1.380649e-23
    wavelength_m = wavelength_nm * 1.0e-9
    radiance = (
        2.0
        * c
        / wavelength_m**4
        / np.expm1(h * c / (wavelength_m * k * 5778.0))
        * 1.0e-9
    )
    solid_angle = 6.794e-5
    return radiance * solid_angle


def _inputs() -> xr.Dataset:
    """Deterministic stand-in for the output of ``sasktran2.photchem.actinic_flux``.

    The spectra and profiles are smooth analytic shapes with realistic
    magnitudes so that every code path of the solver is exercised with a rate
    matrix of representative conditioning.
    """
    altitude = np.arange(50.0e3, 130.0e3 + 1.0, 5.0e3)
    z_km = altitude / 1.0e3

    wavelength = np.unique(
        np.concatenate(
            [np.arange(120.0, 1280.0 + 0.5, 1.0), [LYMAN_ALPHA_WAVELENGTH_NM]]
        )
    )

    temperature = 200.0 + 40.0 * np.exp(-(((z_km - 50.0) / 12.0) ** 2))
    temperature += 60.0 * np.clip((z_km - 95.0) / 35.0, 0.0, None) ** 2

    total = 2.5e25 * np.exp(-z_km / 7.0)
    n2_density = 0.78 * total
    o2_density = 0.21 * total
    co2_density = 4.0e-4 * total
    o3_density = 1.0e18 * np.exp(-(((z_km - 40.0) / 9.0) ** 2)) + 2.0e14 * np.exp(
        -(((z_km - 92.0) / 6.0) ** 2)
    )
    o_density = 4.0e17 * np.exp(-(((z_km - 96.0) / 9.0) ** 2)) + 1.0e14

    wl = wavelength[:, np.newaxis]
    o2_xs = (
        1.0e-21 * np.exp(-(((wl - 145.0) / 18.0) ** 2))
        + 6.0e-28 * np.exp(-(((wl - 220.0) / 25.0) ** 2))
        + 2.0e-27 * np.exp(-(((wl - 762.0) / 3.0) ** 2))
        + 2.0e-28 * np.exp(-(((wl - 689.0) / 3.0) ** 2))
        + 1.0e-29 * np.exp(-(((wl - 1270.0) / 4.0) ** 2))
    ) * (1.0 + 2.0e-4 * (temperature[np.newaxis, :] - 200.0))
    o3_xs = (
        1.1e-21 * np.exp(-(((wl - 255.0) / 20.0) ** 2))
        + 4.5e-25 * np.exp(-(((wl - 600.0) / 60.0) ** 2))
    ) * (1.0 + 1.0e-4 * (temperature[np.newaxis, :] - 200.0))

    # Vertical columns above each altitude for a simple attenuation model.
    scale_height_m = 7.0e3
    o2_column = o2_density * scale_height_m
    o3_column = o3_density * 4.0e3
    airmass = 1.0 / 0.3
    tau = (
        o2_xs * o2_column[np.newaxis, :] + o3_xs * o3_column[np.newaxis, :]
    ) * airmass
    actinic_flux = _photon_flux_toa(wavelength)[:, np.newaxis] * np.exp(-tau)

    lyman_alpha_flux = 3.2e15 * np.exp(-1.0e-24 * o2_column * airmass)

    return xr.Dataset(
        {
            "actinic_flux": (("wavelength", "altitude"), actinic_flux),
            "lyman_alpha_actinic_flux": (("altitude",), lyman_alpha_flux),
            "o2_xs": (("wavelength", "altitude"), o2_xs),
            "o3_xs": (("wavelength", "altitude"), o3_xs),
            "temperature": (("altitude",), temperature),
            "o2_density": (("altitude",), o2_density),
            "o3_density": (("altitude",), o3_density),
            "n2_density": (("altitude",), n2_density),
            "co2_density": (("altitude",), co2_density),
            "o_density": (("altitude",), o_density),
        },
        coords={"altitude": altitude, "wavelength": wavelength},
    )


def _outputs() -> dict[str, np.ndarray]:
    inputs = _inputs()
    model = Yankovsky()

    state = model.solve(inputs)
    emissions = model.emissions(state=state)
    green_line = model.oxygen_green_line_mcdade(inputs)

    outputs = {}
    for prefix, ds in (
        ("state", state),
        ("emissions", emissions),
        ("green_line", green_line),
    ):
        for name, values in ds.data_vars.items():
            outputs[f"{prefix}::{name}"] = values.to_numpy()
    return outputs


def test_yankovsky_matches_golden_outputs():
    expected = np.load(GOLDEN_FILE)
    actual = _outputs()

    assert sorted(actual) == sorted(expected.files)
    for key in expected.files:
        reference = expected[key]
        atol = ATOL_FRACTION * np.max(np.abs(reference), initial=0.0)
        np.testing.assert_allclose(
            actual[key], reference, rtol=RTOL, atol=atol, err_msg=key
        )


if __name__ == "__main__":
    if "--regenerate" not in sys.argv:
        sys.exit("Pass --regenerate to overwrite the golden fixtures.")
    GOLDEN_FILE.parent.mkdir(parents=True, exist_ok=True)
    np.savez(GOLDEN_FILE, **_outputs())
    print(f"Wrote {GOLDEN_FILE}")
