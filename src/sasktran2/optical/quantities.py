from __future__ import annotations

from sasktran2._core_rust import PyOpticalQuantities


class OpticalQuantities:
    _py_oq: PyOpticalQuantities

    def __init__(self, py_oq: PyOpticalQuantities):
        self._py_oq = py_oq

    @property
    def ssa(self):
        """Single scattering albedo"""
        return self._py_oq.ssa

    @property
    def extinction(self):
        """Extinction"""
        return self._py_oq.extinction

    @property
    def leg_coeff(self):
        """Legendre coefficients"""
        return self._py_oq.leg_coeff.transpose([2, 0, 1])

    # Deprecated properties? Not sure why these ever existed
    @property
    def d_ssa(self):
        return self.ssa

    @property
    def d_extinction(self):
        return self.extinction

    @property
    def d_leg_coeff(self):
        return self.leg_coeff
