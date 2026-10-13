"""The O2 UV cross-section builder (offline) and the O2UV optical property."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
import sasktran2 as sk


def _builder():
    path = Path(__file__).resolve().parents[2] / "tools/spectroscopy/build_o2_uv.py"
    spec = importlib.util.spec_from_file_location("build_o2_uv", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write_minschwaner(path: Path, scale: float) -> None:
    rows = [
        f" {nu:.1f}  {1e-6 * scale:.3E} {2e-4 * scale:.3E} {3e-2 * scale:.3E}  0.1  200."
        for nu in (49000.5, 51000.0, 57000.0)
    ]
    header = [
        " 6",
        " wavenum, a0, a1, a2, maxerr, temp(maxerr)",
        " fit to xsec vs ((T-100)/10)^2",
    ]
    path.write_text("\n".join(header + rows) + "\n")


def _sources(tmp_path: Path) -> tuple[Path, Path]:
    builder = _builder()
    for name, scale in zip(builder.MINSCHWANER_FILES, (1.0, 2.0, 3.0), strict=True):
        _write_minschwaner(tmp_path / name, scale)
    (tmp_path / builder.YOSHINO_FILE).write_text(
        "205\t7.35e-24\r\n220\t4.46e-24\r\n240\t1.01e-24\r\n"
    )
    cfa = tmp_path / "O2SCHRUNG"
    cfa.write_text(
        "   O2-SR 49000.0000 61000.0000      4     90.0 101325.00  0.000E+00  0  X\n"
        " 1.000E-20 2.000E-20 3.000E-20 4.000E-20\n"
        "   O2-SR 49000.0000 61000.0000      4    300.0 101325.00  0.000E+00  0  X\n"
        " 2.000E-20 4.000E-20 6.000E-20 8.000E-20\n"
    )
    return tmp_path, cfa


def test_minschwaner_polynomial_and_temperature_ranges():
    builder = _builder()
    coeffs = [
        (190.0, np.array([[1.0, 2.0, 3.0]])),
        (np.inf, np.array([[0.0, 0.0, 5.0]])),
    ]
    # T = 120 K: X = 4, sigma = 1*16 + 2*4 + 3 = 27 (1e-20 cm^2).
    np.testing.assert_allclose(
        builder.minschwaner_cross_section_m2(coeffs, 120.0), 27.0e-24
    )
    np.testing.assert_allclose(
        builder.minschwaner_cross_section_m2(coeffs, 250.0), 5.0e-24
    )


def test_cfa_reader_splits_temperature_blocks(tmp_path):
    builder = _builder()
    _, cfa = _sources(tmp_path)
    blocks = builder.read_cfa_xsc(cfa)
    assert sorted(blocks) == [90.0, 300.0]
    nu, xs = blocks[300.0]
    np.testing.assert_allclose(nu, [49000.0, 53000.0, 57000.0, 61000.0])
    np.testing.assert_allclose(xs, [2e-24, 4e-24, 6e-24, 8e-24])


def test_herzberg_extension():
    builder = _builder()
    table_nm = np.array([205.0, 220.0, 240.0])
    table_m2 = np.array([7.35e-28, 4.46e-28, 1.01e-28])
    wavelengths = np.array([190.0, 200.0, 205.0, 241.2, 243.0])
    xs = builder.herzberg_continuum_m2(1e7 / wavelengths, table_nm, table_m2)
    np.testing.assert_allclose(xs, [0.0, 7.35e-28, 7.35e-28, 0.5 * 1.01e-28, 0.0])


def test_build_combines_the_three_regions(tmp_path):
    builder = _builder()
    source_dir, cfa = _sources(tmp_path)
    ds = builder.build(source_dir, cfa)

    assert ds["xs"].dims == ("temperature_k", "wavenumber_cminv")
    assert np.all(np.diff(ds["wavenumber_cminv"].to_numpy()) > 0)
    xs = ds["xs"].sel(temperature_k=195.0)
    # Bands: mid coefficients (scale 2) at X = 90.25, plus the continuum under them.
    x = ((195.0 - 100.0) / 10.0) ** 2
    band = (2e-6 * x**2 + 4e-4 * x + 6e-2) * 1e-24
    np.testing.assert_allclose(xs.sel(wavenumber_cminv=51000.0), band + 7.35e-28)
    # Continuum above the bands: CfA interpolated between 90 K and 300 K.
    np.testing.assert_allclose(
        xs.sel(wavenumber_cminv=61000.0), 8e-24 * (0.5 + 0.5 * (195.0 - 90.0) / 210.0)
    )
    assert "source_sha256" in ds.attrs


def test_o2uv_optical_property():
    try:
        o2 = sk.optical.O2UV()
    except OSError:
        pytest.skip("O2 UV table not in the local sasktran2 database")

    altitudes = np.array([0.0, 50.0e3])
    geometry = sk.Geometry1D(
        0.5,
        0.0,
        6371000.0,
        altitudes,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.PseudoSpherical,
    )
    atmosphere = sk.Atmosphere(
        geometry, sk.Config(), np.array([145.0, 180.0, 220.0, 250.0])
    )
    atmosphere.temperature_k = np.array([288.0, 200.0])
    atmosphere.pressure_pa = np.array([101325.0, 80.0])

    xs_cm2 = np.asarray(o2.atmosphere_quantities(atmosphere).extinction) * 1e4
    # Schumann-Runge continuum peak, temperature-dependent bands, Herzberg
    # continuum, and nothing beyond the dissociation threshold.
    assert 1.0e-17 < xs_cm2[0, 0] < 2.0e-17
    assert xs_cm2[0, 1] > xs_cm2[1, 1] > 0.0
    np.testing.assert_allclose(xs_cm2[:, 2], 4.46e-24, rtol=1e-6)
    assert np.all(xs_cm2[:, 3] == 0.0)
