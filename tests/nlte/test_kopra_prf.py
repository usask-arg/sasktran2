"""Offline checks for the KOPRA .prf readers in tools/nlte."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


def _reader():
    path = Path(__file__).resolve().parents[2] / "tools/nlte/kopra_prf.py"
    spec = importlib.util.spec_from_file_location("kopra_prf", path)
    module = importlib.util.module_from_spec(spec)
    # dataclasses resolve their module through sys.modules.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


RATIO = """$
 3

$
  0.000000E+00  1.000000E+00  2.000000E+00

$
 T

$
2

$        X      0      0.0000
  71
   1
  1.000000E+00  1.000000E+00  1.000000E+00

$        X     35      0.0000
  71
   2
  1.000000E+00  2.487055+100  6.273800-101
"""

MIXER_PT = """Mixer generated PT profiles
number of levels
$
 3
altitude [km]
$
  0.0  1.0  2.0
pressure [hPa]
$  in_pt.dat
  1.000E+03  8.900E+02  7.900E+02
temperature [K]
$  in_pt.dat
  2.880E+02  2.810E+02  2.750E+02
"""

VMR = """$ Number of vertical profile points
          3

$ Altitudes (km)
  0.00000e+00  1.00000e+00  2.00000e+00

$ Number of profiles
          2

$ (ppmv) H2O (HITRAN)    source: ERS v7
       1
  1.0e+04  2.0e+03  3.0e+02

$ (ppmv) O(1D)           source: ERS v7
     151
  0.0e+00  1.0e-12  2.0e-12
"""

NPAR = """Mixer generated NLTE  profile collection
number of levels
$
 3
altitude [km]
$
  0.0  1.0  2.0
number of nlte profiles given below
$
  1
nlte parameter number
profile
$    in_npar.dat
    3 J_O3
  6.3E-04  6.5E-04  6.6E-04
"""


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_ratio_reader_parses_fortran_three_digit_exponents(tmp_path):
    ds = _reader().read_ratio(_write(tmp_path, "ratio.prf", RATIO))

    np.testing.assert_array_equal(ds.altitude_km, [0.0, 1.0, 2.0])
    np.testing.assert_array_equal(ds.species_id, [71, 71])
    np.testing.assert_array_equal(ds.state_index, [1, 2])
    assert ds.label.values[1] == "X     35      0.0000"
    np.testing.assert_allclose(ds.ratio[1], [1.0, 2.487055e100, 6.2738e-101])


def test_pt_reader_handles_mixer_labels_before_blocks(tmp_path):
    ds = _reader().read_pt(_write(tmp_path, "pt.prf", MIXER_PT))

    np.testing.assert_allclose(ds.pressure_hpa, [1000.0, 890.0, 790.0])
    np.testing.assert_allclose(ds.temperature_k, [288.0, 281.0, 275.0])


def test_vmr_reader_names_species_and_sources(tmp_path):
    ds = _reader().read_vmr(_write(tmp_path, "vmr.prf", VMR))

    assert list(ds.species.values) == ["H2O", "O(1D)"]
    np.testing.assert_array_equal(ds.species_id, [1, 151])
    assert list(ds.source.values) == ["ERS v7", "ERS v7"]
    np.testing.assert_allclose(ds.vmr_ppmv.sel(species="O(1D)"), [0.0, 1e-12, 2e-12])


def test_npar_reader_keys_profiles_by_name(tmp_path):
    ds = _reader().read_npar(_write(tmp_path, "npar.prf", NPAR))

    np.testing.assert_allclose(ds["J_O3"], [6.3e-4, 6.5e-4, 6.6e-4])
    assert ds["J_O3"].attrs["parameter_number"] == 3
