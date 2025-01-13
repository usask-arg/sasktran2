from __future__ import annotations

import abc

import numpy as np
import xarray as xr
from scipy import sparse

import sasktran2 as sk


def map_derivative(mapping, np_deriv: np.ndarray, dims: list[str]) -> xr.DataArray:
    if mapping.interpolator is None:
        return xr.DataArray(np_deriv, dims=dims)
    return (
        xr.DataArray(
            (
                np_deriv.reshape(-1, np_deriv.shape[-1])
                @ sparse.csc_matrix(mapping.interpolator)
            ).reshape(np.concatenate((np_deriv.shape[:-1], [-1]))),
            dims=dims,
        )
        .transpose(*[dims[-1], *dims[:-1]])
        .rename({dims[-1]: mapping.interp_dim})
    )


def map_surface_derivative(
    mapping, np_deriv: np.ndarray, dims: list[str]
) -> xr.DataArray:
    if mapping.interpolator is None:
        return xr.DataArray(np_deriv, dims=dims)
    return xr.DataArray(
        np.einsum(
            "ijk, il->lijk",
            np_deriv,
            mapping.interpolator,
            optimize=True,
        ),
        dims=[mapping.interp_dim, *dims],
    )


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
    def post_process(
        self,
        atmo: sk.Atmosphere,
        geometry: sk.Geometry1D,
        viewing_geo: sk.ViewingGeometry,
    ):
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
            msg = "nstokes must be 1 or 3"
            raise ValueError(msg)

    def internal_output(self) -> sk.OutputIdealStokes_1 | sk.OutputIdealStokes_3:
        """
        The internal output object that can be passed to the internal engine

        Returns
        -------
        Union[sk.OutputIdealStokes_1, sk.OutputIdealStokes_3]
        """
        return self._output

    def post_process(
        self,
        atmo: sk.Atmosphere,
        geometry: sk.Geometry1D,
        viewing_geo: sk.ViewingGeometry,
    ):
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
        # TODO: A lot of this reshaping and organizing should be done inside C++, not here
        # Maybe some of this should go into the derivative mapping classes too...
        # Have to be a little careful since this can be a slow part of the code

        # Reshape and organize output
        result = xr.Dataset()
        if atmo.wavelengths_nm is not None:
            result.coords["wavelength"] = atmo.wavelengths_nm

        result.coords["stokes"] = ["I", "Q", "U", "V"][: atmo.nstokes]
        result.coords["altitude"] = geometry.altitudes()

        radiance = self._output.radiance.reshape(
            -1, len(viewing_geo.observer_rays), atmo.nstokes
        )
        result["radiance"] = (["wavelength", "los", "stokes"], radiance)

        natmo_grid = atmo.storage.total_extinction.shape[0]
        nwavel = radiance.shape[0]

        d_radiance_raw = self._output.d_radiance.reshape(
            nwavel, len(viewing_geo.observer_rays), atmo.nstokes, -1
        )
        d_radiance_k = d_radiance_raw[:, :, :, :natmo_grid]
        d_radiance_ssa = d_radiance_raw[:, :, :, natmo_grid : 2 * natmo_grid]

        d_radiance_albedo = d_radiance_raw[:, :, :, -atmo.surface.brdf.num_deriv :]

        dk_scaled_by_dk = 1 - atmo.storage.f * atmo.unscaled_ssa
        dssa_scaled_by_dssa = (
            1 - atmo.storage.f
        ) / dk_scaled_by_dk + atmo.storage.f * atmo.storage.ssa / dk_scaled_by_dk

        for deriv_name, mapping in atmo.storage.derivative_mappings.items():
            if mapping.assign_name != "":
                name_to_place_result = mapping.assign_name
            else:
                name_to_place_result = deriv_name

            np_deriv = (
                mapping.d_extinction * dk_scaled_by_dk
                - atmo.storage.f * mapping.d_ssa * atmo.unscaled_extinction
            ).T[:, np.newaxis, np.newaxis, :] * d_radiance_k + (
                mapping.d_ssa * dssa_scaled_by_dssa
            ).T[
                :, np.newaxis, np.newaxis, :
            ] * d_radiance_ssa

            if mapping.is_scattering_derivative:
                d_radiance_scat = d_radiance_raw[
                    :,
                    :,
                    :,
                    (2 + mapping.scat_deriv_index)
                    * natmo_grid : (3 + mapping.scat_deriv_index)
                    * natmo_grid,
                ]

                np_deriv += (
                    d_radiance_scat
                    * mapping.scat_factor.T[:, np.newaxis, np.newaxis, :]
                )

                if atmo.applied_delta_m_order:
                    # Change in f for the scatterer
                    d_f = (
                        atmo.storage.d_f[:, :, mapping.scat_deriv_index]
                        * mapping.scat_factor
                    )

                    # Change in k* due to a change in f
                    np_deriv -= (d_f * atmo.unscaled_ssa * atmo.unscaled_extinction).T[
                        :, np.newaxis, np.newaxis, :
                    ] * d_radiance_k

                    # Change in ssa* due to a change in f
                    np_deriv += (
                        d_f
                        * atmo.unscaled_ssa
                        / (1 - atmo.unscaled_ssa * atmo.storage.f)
                        * (atmo.storage.ssa - 1)
                    ).T[:, np.newaxis, np.newaxis, :] * d_radiance_ssa

            if name_to_place_result in result:
                result[name_to_place_result] += map_derivative(
                    mapping, np_deriv, ["wavelength", "los", "stokes", "altitude"]
                )
            else:
                result[name_to_place_result] = map_derivative(
                    mapping, np_deriv, ["wavelength", "los", "stokes", "altitude"]
                )

        for deriv_name, mapping in atmo.surface.derivative_mappings.items():
            name_to_place_result = deriv_name

            np_deriv = d_radiance_albedo * mapping.d_brdf[:, np.newaxis, np.newaxis, :]

            mapped_derivative = map_surface_derivative(
                mapping, np_deriv.sum(axis=-1), ["wavelength", "los", "stokes"]
            )
            result[name_to_place_result] = mapped_derivative

        return result


class OutputDerivMapped(Output):
    def __init__(self, nstokes: int):
        """
        The default output container used by SASKTRAN2.  Here radiances  are returned
        back densely for every input line of sight without modification, and derivatives are mapped
        based on the atmospheric derivative mappings.

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
            self._output = sk.OutputDerivMappedStokes_1()
        elif nstokes == 3:
            self._output = sk.OutputDerivMappedStokes_3()
        else:
            msg = "nstokes must be 1 or 3"
            raise ValueError(msg)

    def internal_output(
        self,
    ) -> sk.OutputDerivMappedStokes_1 | sk.OutputDerivMappedStokes_3:
        """
        The internal output object that can be passed to the internal engine

        Returns
        -------
        Union[sk.OutputDerivMappedStokes_1, sk.OutputDerivMappedStokes_3]
        """
        return self._output

    def post_process(
        self,
        atmo: sk.Atmosphere,
        geometry: sk.Geometry1D,
        viewing_geo: sk.ViewingGeometry,
    ):
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
        # Reshape and organize output
        result = xr.Dataset()
        if atmo.wavelengths_nm is not None:
            result.coords["wavelength"] = atmo.wavelengths_nm

        result.coords["stokes"] = ["I", "Q", "U", "V"][: atmo.nstokes]
        result.coords["altitude"] = geometry.altitudes()

        radiance = self._output.radiance.reshape(
            -1, len(viewing_geo.observer_rays), atmo.nstokes
        )
        result["radiance"] = (["wavelength", "los", "stokes"], radiance)

        for deriv_name, mapping in atmo.storage.derivative_mappings.items():
            if mapping.assign_name != "":
                name_to_place_result = mapping.assign_name
            else:
                name_to_place_result = deriv_name

            # Make sure to copy before reshaping because the internal array is
            # Fortran ordered
            np_deriv = (
                self._output.deriv_map[deriv_name]
                .copy()
                .reshape(
                    -1, len(viewing_geo.observer_rays), atmo.nstokes, mapping.num_output
                )
                .transpose([3, 0, 1, 2])
            )

            if name_to_place_result in result:
                result[name_to_place_result] += xr.DataArray(
                    np_deriv.copy(),
                    dims=[
                        mapping.interp_dim if mapping.interp_dim != "" else "altitude",
                        "wavelength",
                        "los",
                        "stokes",
                    ],
                )
            else:
                result[name_to_place_result] = xr.DataArray(
                    np_deriv.copy(),
                    dims=[
                        mapping.interp_dim if mapping.interp_dim != "" else "altitude",
                        "wavelength",
                        "los",
                        "stokes",
                    ],
                )
        # And the surface derivatives
        for deriv_name, mapping in atmo.surface.derivative_mappings.items():
            name_to_place_result = deriv_name

            # Make sure to copy before reshaping because the internal array is
            # Fortran ordered
            np_deriv = (
                self._output.surface_deriv_map[deriv_name]
                .copy()
                .reshape(-1, len(viewing_geo.observer_rays), atmo.nstokes)
            )

            mapped_derivative = map_surface_derivative(
                mapping, np_deriv, ["wavelength", "los", "stokes"]
            )
            if mapping.interp_dim == "dummy":
                mapped_derivative = mapped_derivative.isel(**{mapping.interp_dim: 0})
            result[name_to_place_result] = mapped_derivative

        return result
