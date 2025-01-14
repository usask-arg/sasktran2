from __future__ import annotations

import numpy as np


class LegendreStorageView:
    def __init__(self, leg_coeff_stacked: np.ndarray, nstokes: int):
        """
        Convenience class to allow easier interfacing with the stacked Legendre coefficient format used
        internally in SASKTRAN2 for performance reasons.

        Internally SASKTRAN2 contains an array of Legendre coefficient (ldim, wavelength, geometric location),
        where ldim is a stacked Legendre dimension.  In the scalar case this is simply the Legendre coefficients
        :math:`\beta_l`, but in the `nstokes=3` case this dimension is :math:`a_{1,1}, a_{2, 1}, a_{3, 1}, b_{1, 1}, a_{1, 2},...`
        where the second coefficient is the expansion index.

        This class provides accessors to set/view :math:`a_1, a_2, a_3, b_1`.


        Parameters
        ----------
        leg_coeff_stacked : np.ndarray
            The raw stacked legendre coefficient data
        nstokes : int
            The number of stokes parameters being included.  Currently only 1 and 3 are supported.

        Raises
        ------
        ValueError
            If nstokes is not 1 or 3
        """
        if nstokes == 1:
            num_coeff = 1
        elif nstokes == 3:
            num_coeff = 4
        else:
            msg = "Only nstokes=1 and nstokes=3 are supported by LegendreStorageView"
            raise ValueError(msg)

        self._raw = leg_coeff_stacked.view()

        self._a1 = self._raw[::num_coeff].view()

        if nstokes == 3:
            self._a2 = self._raw[1::num_coeff].view()
            self._a3 = self._raw[2::num_coeff].view()
            self._b1 = self._raw[3::num_coeff].view()
        else:
            self._a2 = None
            self._a3 = None
            self._b1 = None

    @property
    def a1(self):
        return self._a1

    @property
    def a2(self):
        return self._a2

    @property
    def a3(self):
        return self._a3

    @property
    def b1(self):
        return self._b1
