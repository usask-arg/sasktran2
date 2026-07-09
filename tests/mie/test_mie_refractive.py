from __future__ import annotations

from sasktran2.mie.refractive import Water


def test_aria_refractive_index():
    water = Water(source="aria")
    k = water.refractive_index("550.0")
    assert k == 1.333 - 1.96e-9j
