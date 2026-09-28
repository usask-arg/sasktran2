from __future__ import annotations

from sasktran2.mie.refractive import Water, _from_aria_file


def test_aria_refractive_index():

    _from_aria_file("H2O_Mcgarragh_2018", download=True, extract=True)

    water = Water(source="aria")
    k = water.refractive_index("550.0")
    assert k == 1.333 - 1.96e-9j
