from __future__ import annotations

import numpy as np
import pytest
from sasktran2.optical import AERLineAbsorber

species_list = [
    "H2O",
    "CO2",
    "O3",
    "N2O",
    "CO",
    "CH4",
    "O2",
    "NO",
    "SO2",
    "NO2",
    "NH3",
    "HNO3",
    "OH",
    "HF",
    "HCL",
    "HBR",
    "HI",
    "CLO",
    "OCS",
    "H2CO",
    "HOCL",
    "N2",
    "HCN",
    "CH3CL",
    "H2O2",
    "C2H2",
    "C2H6",
    "PH3",
    "COF2",
    "SF6",
    "H2S",
    "HCOOH",
    "HO2",
    "O",
    "CLONO2",
    "NO+",
    "HOBr",
    "C2H4",
    "CH3OH",
    "CH3Br",
    "CH3CN",
    "CF4",
    "C4H2",
    "HC3N",
    "H2",
    "CS",
    "SO3",
]


@pytest.mark.parametrize("species", species_list)
def test_aer_line_loading(species):
    pytest.importorskip("zenodo_get")

    if species == "H2":
        pytest.skip("H2 data is misformed and not fully supported")

    _ = AERLineAbsorber(species)


@pytest.mark.parametrize("species", species_list)
def test_aer_line_cross_section(species):
    pytest.importorskip("zenodo_get")
    pytest.importorskip("hapi")

    if species == "O":
        pytest.skip("Atomic oxygen partition sum data missing?")
    if species == "H2":
        pytest.skip("H2 data is misformed and not fully supported")

    absorber = AERLineAbsorber(species)
    wavenumbers = np.arange(1000, 2000, 10.0)
    _ = absorber.cross_sections(
        wavenumbers, [0.0], pressure_pa=[101325], temperature_k=[273.15]
    )
