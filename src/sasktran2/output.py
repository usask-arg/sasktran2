import abc
import sasktran2 as sk
import numpy as np
import xarray as xr


class Output(abc.ABC):
    @abc.abstractmethod
    def internal_output(self):
        pass

    @abc.abstractmethod
    def post_process(self, atmo: sk.Atmosphere, geometry: sk.Geometry1D, viewing_geo: sk.ViewingGeometry):
        pass


class OutputIdeal(Output):
    def __init__(self, nstokes: int):
        self._nstokes = nstokes

        if nstokes == 1:
            self._output = sk.OutputIdealStokes_1()
        elif nstokes == 3:
            self._output = sk.OutputIdealStokes_3()
        else:
            msg = 'nstokes must be 1 or 3'
            raise ValueError(msg)

    def internal_output(self):
        return self._output

    def post_process(self, atmo: sk.Atmosphere, geometry: sk.Geometry1D, viewing_geo: sk.ViewingGeometry):
        # TODO: A lot of this reshaping and organizing should be done inside C++

        # Reshape and organize output
        result = xr.Dataset()
        if atmo.wavelengths_nm is not None:
            result.coords['wavelength'] = atmo.wavelengths_nm

        result.coords['stokes'] = ['I', 'Q', 'U', 'V'][:atmo.nstokes]
        result.coords['altitude'] = geometry.altitudes()

        radiance = self._output.radiance.reshape(atmo.nstokes, len(viewing_geo.observer_rays), -1)
        result['radiance'] = (['stokes', 'los', 'wavelength'], radiance)

        natmo_grid = atmo.storage.total_extinction.shape[0]
        nwavel = radiance.shape[2]

        d_radiance_raw = self._output.d_radiance.reshape(atmo.nstokes, len(viewing_geo.observer_rays), nwavel, -1).transpose((0, 1, 3, 2))
        d_radiance_k = d_radiance_raw[:, :, :natmo_grid, :]
        d_radiance_ssa = d_radiance_raw[:, :, natmo_grid:2*natmo_grid, :]

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
                np_deriv = (mapping.native_grid_mapping.d_extinction * dk_scaled_by_dk - atmo.storage.f*mapping.native_grid_mapping.d_ssa*atmo.unscaled_extinction)[np.newaxis, np.newaxis, :, :] * d_radiance_k +\
                (mapping.native_grid_mapping.d_ssa*dssa_scaled_by_dssa)[np.newaxis, np.newaxis, :, :]*d_radiance_ssa

                if name_to_place_result in result:
                    result[name_to_place_result].values += np_deriv
                else:
                    result[name_to_place_result] = (['stokes', 'los', 'altitude', 'wavelength'], np_deriv)

        return result
