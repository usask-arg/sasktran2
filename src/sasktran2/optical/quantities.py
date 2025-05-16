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