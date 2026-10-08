"""Photolysis and photoexcitation rates for the bundled ``sasktran2-nlte`` mechanisms."""

from __future__ import annotations

import xarray as xr

from .flux import O2_LINE_WINDOWS_NM
from .quantum_yields import o3_o1d_matsumi2002, o3_o3p_matsumi2002
from .rates import LymanAlphaPhotolysis, Photolysis, photolysis_rates

#: Line-integrated solar Lyman-alpha photon flux at 1 AU [photons m^-2 s^-1].
LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S = 3.2e15
#: O(1D) yield of O2 photolysis at Lyman-alpha.
O2_LYMAN_ALPHA_O1D_YIELD = 0.53

_MATSUMI = "O3DBM cross sections; O(1D) yield of Matsumi et al. (2002)"
_O2 = "AER O2 line cross sections, resolved at 0.001 nm"


def oxygen_photolysis() -> list[Photolysis | LymanAlphaPhotolysis]:
    """Physical channel rates for O3 and O2 photolysis and O2 photoexcitation.

    ===============  ============================================================
    ``J_O3_O1D``     O3 -> O2 + O(1D)
    ``J_O3_O3P``     O3 -> O2 + O(3P)
    ``J_O2_SRC``     O2 -> O(3P) + O(1D) in the Schumann-Runge continuum (130-175 nm)
    ``J_O2_LYA``     O2 -> O(3P) + O(1D) at Lyman-alpha
    ``J_O2_EXC_B0``  O2(X) -> O2(b, v=0), A band
    ``J_O2_EXC_B1``  O2(X) -> O2(b, v=1), B band
    ``J_O2_EXC_B2``  O2(X) -> O2(b, v=2), gamma band
    ``J_O2_EXC_A0``  O2(X) -> O2(a, v=0), 1.27 um band
    ===============  ============================================================

    ``J_O2_LYA`` uses the Chabrillat and Kockarts (1997) slant-column
    parameterisation. The excitation rates integrate all O2 absorption in each
    band's window of :data:`O2_LINE_WINDOWS_NM`.
    """
    excitation = {
        "J_O2_EXC_B0": "b-X(0,0) A band",
        "J_O2_EXC_B1": "b-X(1,0) B band",
        "J_O2_EXC_B2": "b-X(2,0) gamma band",
        "J_O2_EXC_A0": "a-X(0,0) 1.27 um band",
    }
    return [
        Photolysis("J_O3_O1D", "O3", o3_o1d_matsumi2002, reference=_MATSUMI),
        Photolysis("J_O3_O3P", "O3", o3_o3p_matsumi2002, reference=_MATSUMI),
        Photolysis(
            "J_O2_SRC",
            "O2",
            wavelength_range_nm=(130.0, 175.0),
            reference="O2SchumannRunge continuum cross sections, unit O(1D) yield",
        ),
        LymanAlphaPhotolysis(
            "J_O2_LYA",
            toa_line_flux_photons_m2_s=LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S,
            quantum_yield=O2_LYMAN_ALPHA_O1D_YIELD,
        ),
        *(
            Photolysis(
                name,
                "O2",
                wavelength_range_nm=O2_LINE_WINDOWS_NM[window],
                max_grid_spacing_nm=0.001,
                reference=_O2,
            )
            for name, window in excitation.items()
        ),
    ]


#: Fractions of the O(1D) channel producing O2(a, v=0..5), and the O2(X, v)
#: levels sharing the O(3P) channel equally, as in the legacy Yankovsky
#: model. They are product distributions, kept here only to build that
#: mechanism's per-channel rate inputs.
YANKOVSKY_O2A_FRACTIONS = tuple(
    q / 0.90 for q in (0.441, 0.135, 0.135, 0.072, 0.072, 0.045)
)
YANKOVSKY_O2X_LEVELS = range(1, 36)


def oxygen_yankovsky_rates(flux: xr.Dataset) -> xr.Dataset:
    """The rate inputs of the bundled ``oxygen_yankovsky`` mechanism [s^-1].

    ``flux`` is the output of :meth:`ActinicFlux.calculate`. The O3 rates
    split the physical O(1D) and O(3P) channels with
    :data:`YANKOVSKY_O2A_FRACTIONS` and :data:`YANKOVSKY_O2X_LEVELS`.
    """
    channels = photolysis_rates(flux, oxygen_photolysis())
    rates = channels.drop_vars(["J_O3_O1D", "J_O3_O3P"])
    for v, fraction in enumerate(YANKOVSKY_O2A_FRACTIONS):
        rates[f"J_O3_A{v}"] = channels["J_O3_O1D"] * fraction
    for v in YANKOVSKY_O2X_LEVELS:
        rates[f"J_O3_X{v}"] = channels["J_O3_O3P"] / len(YANKOVSKY_O2X_LEVELS)
    for name in rates.data_vars:
        rates[name].attrs["units"] = "s^-1"
    return rates
