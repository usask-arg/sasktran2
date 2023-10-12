import abc
from typing import Union

import numpy as np
import xarray as xr

import sasktran2 as sk


class Output(abc.ABC):
    """
    An abstract class which defines how output from the radiative transfer Engine is handled.

    The Engine first calls :py:meth:`~internal_output` to get a `pybind11` object which can be passed
    to the internal radiative transfer model. After the calculation is completed, :py:meth:`~post_process`
    is called to get the exact output container that is passed back to the user.

    """
    @abc.abstractmethod
    def internal_output(self):
        pass

    @abc.abstractmethod
    def post_process(self, atmo: sk.Atmosphere, geometry: sk.Geometry1D, viewing_geo: sk.ViewingGeometry):
        pass


class OutputIdeal(Output):
    def __init__(self, nstokes: int):
        """
        The default output container used by SASKTRAN2.  Here radiances (and derivatives) are returned
        back densely for every input line of sight and wavelength without modification.

        Parameters
        ----------
        nstokes : int
            Number of stokes parameters to store for the output quantities

        Raises
        ------
        ValueError
            If nstokes is not 1 or 3
        """
        self._nstokes = nstokes

        if nstokes == 1:
            self._output = sk.OutputIdealStokes_1()
        elif nstokes == 3:
            self._output = sk.OutputIdealStokes_3()
        else:
            msg = 'nstokes must be 1 or 3'
            raise ValueError(msg)

    def internal_output(self) -> Union[sk.OutputIdealStokes_1, sk.OutputIdealStokes_3]:
        """
        The internal output object that can be passed to the internal engine

        Returns
        -------
        Union[sk.OutputIdealStokes_1, sk.OutputIdealStokes_3]
        """
        return self._output

    def post_process(self, atmo: sk.Atmosphere, geometry: sk.Geometry1D, viewing_geo: sk.ViewingGeometry):
        """
        Converts the raw output values to a more usable format.  Also performs the mapping of raw model
        derivatives to derivatives of the constituent input quantities.

        Parameters
        ----------
        atmo : sk.Atmosphere
            Atmosphere object after calculation has been performed
        geometry : sk.Geometry1D
            Geometry object
        viewing_geo : sk.ViewingGeometry
            Viewing geometry

        Returns
        -------
        xr.Dataset
        """
        # TODO: A lot of this reshaping and organizing should be done inside C++

        # Reshape and organize output
        result = xr.Dataset()
        if atmo.wavelengths_nm is not None:
            result.coords['wavelength'] = atmo.wavelengths_nm

        result.coords['stokes'] = ['I', 'Q', 'U', 'V'][:atmo.nstokes]
        result.coords['altitude'] = geometry.altitudes()

        radiance = self._output.radiance.reshape(atmo.nstokes, -1, len(viewing_geo.observer_rays))
        result['radiance'] = (['stokes', 'wavelength', 'los'], radiance)

        natmo_grid = atmo.storage.total_extinction.shape[0]
        nwavel = radiance.shape[1]

        d_radiance_raw = self._output.d_radiance.reshape(atmo.nstokes, nwavel, len(viewing_geo.observer_rays), -1)
        d_radiance_k = d_radiance_raw[:, :, :, :natmo_grid]
        d_radiance_ssa = d_radiance_raw[:, :, :,  natmo_grid:2*natmo_grid]

        dk_scaled_by_dk = (1 - atmo.storage.f * atmo.unscaled_ssa)
        dssa_scaled_by_dssa = (1 - atmo.storage.f) / dk_scaled_by_dk + atmo.storage.f * atmo.storage.ssa / dk_scaled_by_dk

        for constituent_name, deriv_mapping in atmo.deriv_mappings.items():
            if deriv_mapping is None:
                continue

            for deriv_name, mapping in deriv_mapping.items():
                if mapping.summable:
                    name_to_place_result = 'wf_' + deriv_name
                else:
                    name_to_place_result = 'wf_' + constituent_name + '_' + deriv_name

                # Start by calculating the derivative with respect to the quantity on the native grid
                np_deriv = (mapping.native_grid_mapping.d_extinction * dk_scaled_by_dk - atmo.storage.f*mapping.native_grid_mapping.d_ssa*atmo.unscaled_extinction).T[np.newaxis, :, np.newaxis, :] * d_radiance_k +\
                (mapping.native_grid_mapping.d_ssa*dssa_scaled_by_dssa).T[np.newaxis, :, np.newaxis, :]*d_radiance_ssa

                if mapping.log_radiance_space:
                    np_deriv /= radiance[:, :, :, np.newaxis]

                if name_to_place_result in result:
                    result[name_to_place_result].values += np_deriv
                else:
                    result[name_to_place_result] = (['stokes', 'wavelength', 'los', 'altitude'], np_deriv)

        return result
