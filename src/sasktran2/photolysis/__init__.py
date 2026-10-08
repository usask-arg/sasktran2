"""Actinic flux and photolysis rates, in the spirit of TUV.

:class:`ActinicFlux` runs SASKTRAN2's discrete-ordinates engine with flux
observers at every altitude; :func:`photolysis_rates` integrates cross
sections and quantum yields against the result. The atmosphere is an
``xr.Dataset`` on an ``altitude`` coordinate [m] with ``temperature_k``,
``pressure_pa`` and number densities [m^-3] named by species id, the same
convention :func:`sasktran2.nlte.solve` uses, so one dataset drives both::

    flux = sk.photolysis.ActinicFlux(altitudes_m).calculate(atmosphere, cos_sza=0.5)
    rates = sk.photolysis.photolysis_rates(flux, sk.photolysis.presets.oxygen_photolysis())
    mechanism = sk.nlte.Mechanism.bundled("oxygen")
    solution = sk.nlte.solve(mechanism, atmosphere, rates)

:func:`sasktran2.nlte.add_photochemical_species` does all of this for a
radiative-transfer atmosphere and adds the emission.
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
from .tuv import (
    TUV_OPTICAL_PROPERTIES,
    TUVActinicFlux,
    TUVXQuantumYield,
    chabrillat_kockarts_cross_section,
    hsrs_photons_per_bin,
    koppers_murtagh_cross_section,
    tuvx_v54_tables,
    tuvx_v54_wavelength_edges,
)

__all__ = [
    "DEFAULT_OPTICAL_PROPERTIES",
    "LYMAN_ALPHA_WAVELENGTH_NM",
    "O2_LINE_WINDOWS_NM",
    "O2_SCHUMANN_RUNGE_BANDS_NM",
    "TUV_OPTICAL_PROPERTIES",
    "ActinicFlux",
    "LinePhotolysis",
    "LymanAlphaPhotolysis",
    "Photolysis",
    "TUVActinicFlux",
    "TUVXQuantumYield",
    "airglow_wavelength_grid",
    "chabrillat_kockarts_cross_section",
    "default_optical_properties",
    "hsrs_photons_per_bin",
    "koppers_murtagh_cross_section",
    "lyman_alpha_o2_rate_per_photon_m2",
    "lyman_alpha_reduction_factor",
    "o3_o1d_matsumi2002",
    "o3_o3p_matsumi2002",
    "photolysis_rates",
    "presets",
    "slant_columns",
    "tuvx_v54_tables",
    "tuvx_v54_wavelength_edges",
]
