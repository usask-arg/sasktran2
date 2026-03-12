from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2.util.interpolation import linear_interpolating_matrix

from .base import Constituent


class MTCKDContinuum(Constituent):
    def __init__(
        self, h2o_name: str = "H2O", co2_name: str = "CO2", o3_name: str = "O3"
    ):
        """
        The MT-CKD continuum absorption model. Requires the saskran2_ext package.

        Parameters
        ----------
        h2o_name : str, optional
            The name of the H2O constituent in the atmosphere., by default "H2O"
        co2_name : str, optional
            The name of the CO2 constituent in the atmosphere., by default "CO2"
        o3_name : str, optional
            The name of the O3 constituent in the atmosphere., by default "O3"
        """
        self._h2o_name = h2o_name
        self._co2_name = co2_name
        self._o3_name = o3_name
        self._mtckd_wavenumbers = np.arange(
            -10, 19910, 10
        )  # grid hardcoded in fortran module

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        """
        Parameters
        ----------
        atmo : sk.Atmosphere


        :meta private:
        """
        try:
            from sasktran2_ext import mt_ckd
        except ImportError:
            msg = "sasktran2_ext is required to use the MT-CKD continuum constituent"
            raise ImportError(msg)  # noqa: B904

        if atmo.wavelengths_nm is None:
            msg = "It is required to give the Atmosphere object wavelengths to use the continuum constituent"
            raise ValueError(msg)

        if atmo.pressure_pa is None:
            msg = "It is required to set the pressure_pa property in the Atmosphere object to use the continuum constituent"
            raise ValueError(msg)

        if atmo.temperature_k is None:
            msg = "It is required to set the temperature_k property in the Atmosphere object to use the continuum constituent"
            raise ValueError(msg)

        if atmo[self._h2o_name] is None:
            msg = f"It is required to add an {self._h2o_name} constituent to the Atmosphere object to use the continuum constituent"
            raise ValueError(msg)

        if atmo[self._co2_name] is None:
            msg = f"It is required to add an {self._co2_name} constituent to the Atmosphere object to use the continuum constituent"
            raise ValueError(msg)

        if atmo[self._o3_name] is None:
            msg = f"It is required to add an {self._o3_name} constituent to the Atmosphere object to use the continuum constituent"
            raise ValueError(msg)

        # interpolate vmrs to altitude grid
        alts = atmo.model_geometry.altitudes()

        h2o_vmr = np.interp(
            alts,
            atmo[self._h2o_name].altitudes_m,
            atmo[self._h2o_name].vmr,
            left=atmo[self._h2o_name].vmr[0],
            right=atmo[self._h2o_name].vmr[-1],
        )
        co2_vmr = np.interp(
            alts,
            atmo[self._co2_name].altitudes_m,
            atmo[self._co2_name].vmr,
            left=atmo[self._co2_name].vmr[0],
            right=atmo[self._co2_name].vmr[-1],
        )
        o3_vmr = np.interp(
            alts,
            atmo[self._o3_name].altitudes_m,
            atmo[self._o3_name].vmr,
            left=atmo[self._o3_name].vmr[0],
            right=atmo[self._o3_name].vmr[-1],
        )

        # mt_ckd returns optical depth. set path length to 1.0 m (100.0 cm) so return value is equivalent to absorption in m^-1
        continuum_absorption = mt_ckd(
            atmo.pressure_pa,
            atmo.temperature_k,
            h2o_vmr,
            co2_vmr,
            o3_vmr,
            100.0,  # path length in cm
        )

        # remove unused portion of array
        continuum_absorption = continuum_absorption[:, 0 : len(self._mtckd_wavenumbers)]

        # interpolate continuum to atmosphere grid
        # f = interp1d(
        #     self._mtckd_wavenumbers,
        #     continuum_absorption,
        #     axis=1,
        #     bounds_error=False,
        #     fill_value=0,
        # )
        # atmo.storage.total_extinction[:] += f(atmo.wavenumbers_cminv)
        interp_matrix = linear_interpolating_matrix(
            self._mtckd_wavenumbers,
            atmo.wavenumbers_cminv,
            "zero",
        )
        atmo.storage.total_extinction[:] += (
            np.nan_to_num(continuum_absorption) @ interp_matrix.T
        )

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        try:
            from sasktran2_ext import mt_ckd
        except ImportError:
            msg = "sasktran2_ext is required to use the MT-CKD continuum constituent"
            raise ImportError(msg)  # noqa: B904

        interp_matrix = linear_interpolating_matrix(
            self._mtckd_wavenumbers,
            atmo.wavenumbers_cminv,
            "zero",
        )

        p_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_pressure_pa")
        t_mapping = atmo.storage.get_derivative_mapping(f"wf_{name}_temperature_k")
        h2o_mapping = atmo.storage.get_derivative_mapping(
            f"wf_{name}_{self._h2o_name}_vmr"
        )
        co2_mapping = atmo.storage.get_derivative_mapping(
            f"wf_{name}_{self._co2_name}_vmr"
        )
        o3_mapping = atmo.storage.get_derivative_mapping(
            f"wf_{name}_{self._o3_name}_vmr"
        )

        pressure_pa = np.copy(atmo.pressure_pa)

        temperature_k = np.copy(atmo.temperature_k)

        alts = atmo.model_geometry.altitudes()

        h2o_vmr = np.interp(
            alts,
            atmo[self._h2o_name].altitudes_m,
            atmo[self._h2o_name].vmr,
            left=atmo[self._h2o_name].vmr[0],
            right=atmo[self._h2o_name].vmr[-1],
        )
        co2_vmr = np.interp(
            alts,
            atmo[self._co2_name].altitudes_m,
            atmo[self._co2_name].vmr,
            left=atmo[self._co2_name].vmr[0],
            right=atmo[self._co2_name].vmr[-1],
        )
        o3_vmr = np.interp(
            alts,
            atmo[self._o3_name].altitudes_m,
            atmo[self._o3_name].vmr,
            left=atmo[self._o3_name].vmr[0],
            right=atmo[self._o3_name].vmr[-1],
        )

        base_cont = mt_ckd(
            pressure_pa,
            temperature_k,
            h2o_vmr,
            co2_vmr,
            o3_vmr,
            100.0,
        )[:, 0 : len(self._mtckd_wavenumbers)]

        fractional_change = 1e-4
        for input_var, d_mapping in zip(
            [pressure_pa, temperature_k, h2o_vmr, co2_vmr, o3_vmr],
            [p_mapping, t_mapping, h2o_mapping, co2_mapping, o3_mapping],
            strict=True,
        ):
            dx = input_var * fractional_change

            input_var += dx
            cont_above = mt_ckd(
                pressure_pa,
                temperature_k,
                h2o_vmr,
                co2_vmr,
                o3_vmr,
                100.0,
            )[:, 0 : len(self._mtckd_wavenumbers)]

            # if input_var[i] >= dx:
            if True:
                # central diff
                input_var -= 2 * dx
                cont_below = mt_ckd(
                    pressure_pa,
                    temperature_k,
                    h2o_vmr,
                    co2_vmr,
                    o3_vmr,
                    100.0,
                )[:, 0 : len(self._mtckd_wavenumbers)]
                input_var += dx

                central_diff_wf = (cont_above - cont_below) / (2 * dx[:, np.newaxis])
            else:
                # forward diff
                central_diff_wf = (cont_above - base_cont) / (dx[:, np.newaxis])

                input_var -= dx

            # f = interp1d(
            #     self._mtckd_wavenumbers,
            #     central_diff_wf,
            #     axis=1,
            #     bounds_error=False,
            #     fill_value=0,
            # )
            # d_mapping.d_extinction[:] += f(atmo.wavenumbers_cminv)
            d_mapping.d_extinction[:] += (
                np.nan_to_num(central_diff_wf) @ interp_matrix.T
            )

            # not sure if correct
            d_ssa = (
                np.nan_to_num(central_diff_wf)
                @ interp_matrix.T
                * (-atmo.storage.ssa)
                / atmo.storage.total_extinction
            )
            d_mapping.d_ssa[:] += d_ssa

            d_mapping.interp_dim = "altitude"
            d_mapping.interpolator = np.eye(len(alts))
