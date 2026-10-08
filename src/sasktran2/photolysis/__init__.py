"""Actinic flux and photolysis rates, in the spirit of TUV.

:class:`ActinicFlux` runs SASKTRAN2's discrete-ordinates engine with flux
observers at every altitude; :func:`photolysis_rates` integrates cross
sections and quantum yields against the result. The atmosphere is an
``xr.Dataset`` on an ``altitude`` coordinate [m] with ``temperature_k``,
``pressure_pa`` and number densities [m^-3] named by species id, the same
convention :func:`sasktran2.nlte.solve` uses, so one dataset drives both::

    flux = sk.photolysis.ActinicFlux(altitudes_m).calculate(atmosphere, cos_sza=0.5)
    rates = sk.photolysis.presets.oxygen_yankovsky_rates(flux)
    mechanism = sk.nlte.Mechanism.bundled("oxygen_yankovsky")
    solution = sk.nlte.solve(mechanism, atmosphere, rates)
"""

from __future__ import annotations

from . import presets
from .flux import (
    DEFAULT_OPTICAL_PROPERTIES,
    LYMAN_ALPHA_WAVELENGTH_NM,
    O2_LINE_WINDOWS_NM,
    O2_SCHUMANN_RUNGE_BANDS_NM,
    ActinicFlux,
    airglow_wavelength_grid,
    default_optical_properties,
    slant_columns,
)
from .quantum_yields import o3_o1d_matsumi2002, o3_o3p_matsumi2002
from .rates import (
    LinePhotolysis,
    LymanAlphaPhotolysis,
    Photolysis,
    lyman_alpha_o2_rate_per_photon_m2,
    lyman_alpha_reduction_factor,
    photolysis_rates,
)

__all__ = [
    "DEFAULT_OPTICAL_PROPERTIES",
    "LYMAN_ALPHA_WAVELENGTH_NM",
    "O2_LINE_WINDOWS_NM",
    "O2_SCHUMANN_RUNGE_BANDS_NM",
    "ActinicFlux",
    "LinePhotolysis",
    "LymanAlphaPhotolysis",
    "Photolysis",
    "airglow_wavelength_grid",
    "default_optical_properties",
    "lyman_alpha_o2_rate_per_photon_m2",
    "lyman_alpha_reduction_factor",
    "o3_o1d_matsumi2002",
    "o3_o3p_matsumi2002",
    "photolysis_rates",
    "presets",
    "slant_columns",
]
