"""Photolysis and photoexcitation rates for the bundled ``sasktran2-nlte`` mechanisms."""

from __future__ import annotations

from .flux import O2_LINE_WINDOWS_NM
from .quantum_yields import (
    O3O1DYield,
    o2_o1s_yield,
    o3_o1d_matsumi2002,
    o3_o3p_matsumi2002,
)
from .rates import LymanAlphaPhotolysis, Photolysis
from .tuv import TUVXQuantumYield

#: Line-integrated solar Lyman-alpha photon flux at 1 AU [photons m^-2 s^-1].
LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S = 3.2e15
#: O(1D) yield of O2 photolysis at Lyman-alpha.
O2_LYMAN_ALPHA_O1D_YIELD = 0.53

_MATSUMI = "O3DBM cross sections; O(1D) yield of Matsumi et al. (2002)"
_O2 = "AER O2 line cross sections, resolved at 0.001 nm"
_MATSUMI_YV2020 = (
    "O3DBM cross sections; O(1D) yield of Matsumi et al. (2002); O2(a, v) split "
    "of Yankovsky and Vorobeva (2020)"
)


def oxygen_photolysis(
    excitation: bool = True,
) -> list[Photolysis | LymanAlphaPhotolysis]:
    """Physical channel rates for O3 and O2 photolysis and O2 photoexcitation.

    =================  ==========================================================
    ``J_O3_O1D``       O3 -> O2 + O(1D)
    ``J_O3_O3P``       O3 -> O2 + O(3P)
    ``J_O3_O1D_A{v}``  O3 -> O2(a, v=0-5) + O(1D), wavelength-dependent split
    ``J_O3_O1D_X``     O3 -> O2(X) + O(1D), the spin-forbidden channel beyond 310 nm
    ``J_O2_SRC``       O2 -> O(3P) + O(1D) in the Schumann-Runge continuum (130-175 nm)
    ``J_O2_LYA``       O2 -> O(3P) + O(1D) at Lyman-alpha
    ``J_O2_EXC_B0``    O2(X) -> O2(b, v=0), A band
    ``J_O2_EXC_B1``    O2(X) -> O2(b, v=1), B band
    ``J_O2_EXC_B2``    O2(X) -> O2(b, v=2), gamma band
    ``J_O2_EXC_A0``    O2(X) -> O2(a, v=0), 1.27 um band
    =================  ==========================================================

    ``J_O2_LYA`` uses the Chabrillat and Kockarts (1997) slant-column
    parameterisation. The excitation rates integrate all O2 absorption in each
    band's window of :data:`O2_LINE_WINDOWS_NM`, and need a line-resolving
    grid; ``excitation=False`` leaves them out, e.g. for
    :class:`TUVActinicFlux`.
    """
    bands = {
        "J_O2_EXC_B0": "b-X(0,0) A band",
        "J_O2_EXC_B1": "b-X(1,0) B band",
        "J_O2_EXC_B2": "b-X(2,0) gamma band",
        "J_O2_EXC_A0": "a-X(0,0) 1.27 um band",
    }
    photolysis = [
        Photolysis("J_O3_O1D", "O3", o3_o1d_matsumi2002, reference=_MATSUMI),
        Photolysis("J_O3_O3P", "O3", o3_o3p_matsumi2002, reference=_MATSUMI),
        *(
            Photolysis(f"J_O3_O1D_A{v}", "O3", O3O1DYield(v), reference=_MATSUMI_YV2020)
            for v in range(6)
        ),
        Photolysis("J_O3_O1D_X", "O3", O3O1DYield("X"), reference=_MATSUMI_YV2020),
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
    ]
    if not excitation:
        return photolysis
    return [
        *photolysis,
        *(
            Photolysis(
                name,
                "O2",
                wavelength_range_nm=O2_LINE_WINDOWS_NM[window],
                max_grid_spacing_nm=0.001,
                reference=_O2,
            )
            for name, window in bands.items()
        ),
    ]


def tuvx_v54_photolysis() -> list[Photolysis]:
    """The O2 and O3 photolysis reactions of TUV-x v5.4, for :class:`TUVActinicFlux`.

    ==============  ========================================================
    ``J_O2``        O2 -> O + O, all bins, unit yield (TUV-x ``O2+hv->O+O``)
    ``J_O3_O1D``    O3 -> O2 + O(1D), TUV-x quantum yields
    ``J_O3_O3P``    O3 -> O2 + O(3P), TUV-x quantum yields
    ==============  ========================================================

    The quantum yields are TUV-x's own, tabulated by bin, so these rates are
    only valid on the TUV-x v5.4 bins.
    """
    reference = "TUV-x v5.4"
    return [
        Photolysis("J_O2", "O2", reference=reference),
        Photolysis("J_O3_O1D", "O3", TUVXQuantumYield("o3_o1d"), reference=reference),
        Photolysis("J_O3_O3P", "O3", TUVXQuantumYield("o3_o3p"), reference=reference),
    ]


def green_line_photolysis() -> list[Photolysis]:
    """``J_O2_O1S``: O2 photodissociation into O(1S), 81-121 nm.

    Needs an actinic flux down to 81 nm (the default grid of
    :func:`airglow_wavelength_grid`) with O2, N2 and O absorbing.
    """
    return [
        Photolysis(
            "J_O2_O1S",
            "O2",
            o2_o1s_yield,
            wavelength_range_nm=(80.0, 121.0),
            reference="Leiden O2 cross sections; GLOW O(1S) yields",
        )
    ]
