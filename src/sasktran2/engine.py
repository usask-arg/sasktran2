from __future__ import annotations

import numpy as np
import xarray as xr

import sasktran2 as sk
from sasktran2._core_rust import PyEngine
from sasktran2.viewinggeo.base import ViewingGeometryContainer


def map_surface_derivative(
    mapping, np_deriv: np.ndarray, dims: list[str]
) -> xr.DataArray:
    if mapping.interpolator is None or len(mapping.interpolator) == 0:
        return xr.DataArray(np_deriv, dims=dims)
    return xr.DataArray(
        np.einsum(
            "ij..., il->lij...",
            np_deriv,
            mapping.interpolator,
            optimize=True,
        ),
        dims=[mapping.interp_dim, *dims],
    )


class Engine:
    _engine: PyEngine

    def __init__(
        self,
        config: sk.Config,
        geometry: sk.Geometry1D,
        viewing_geometry: sk.ViewingGeometry,
    ):
        """
        An Engine is the main class that handles the radiative transfer calculation.  The calculation takes
        place in two components.

        First, upon construction of the Engine, the majority of the geometry information is computed and
        cached.

        The main calculation takes place when calling :py:meth:`~calculate_radiance` with an
        :py:class:`sasktran2.Atmosphere` object where the actual radiative transfer calculation
        is performed.

        Parameters
        ----------
        config : sk.Config
            Configuration object
        model_geometry : sk.Geometry1D
            Geometry for the model
        viewing_geo : sk.ViewingGeometry
            Viewing geometry
        """
        self._engine = PyEngine(
            config._config, geometry._geometry, viewing_geometry._viewing_geometry
        )
        self._config = config
        self._geometry = geometry
        self._viewing_geometry = viewing_geometry

    def calculate_radiance(self, atmosphere: sk.Atmosphere) -> xr.Dataset:
        """
        Performs the radiative transfer calculation for a given atmosphere

        Parameters
        ----------
        atmosphere : sk.Atmosphere
            The atmosphere object containing the atmospheric profile and constituents

        Returns
        -------
        xr.Dataset
            An xarray dataset containing the radiance and derivatives
        """
        output = self._engine.calculate_radiance(atmosphere.internal_object())

        out_ds = xr.Dataset()

        out_ds["radiance"] = xr.DataArray(
            output.radiance,
            dims=["wavelength", "los", "stokes"],
        )

        flux_types = ["upwelling", "downwelling"]
        if len(self._viewing_geometry.flux_observers) > 0:
            # TODO: Grab this from the config

            for i, flux_type in enumerate(flux_types):
                out_ds[f"{flux_type}_flux"] = xr.DataArray(
                    output.flux[i],
                    dims=["wavelength", "flux_location"],
                )

        if atmosphere.wavelengths_nm is not None:
            out_ds.coords["wavelength"] = atmosphere.wavelengths_nm

        out_ds.coords["stokes"] = ["I", "Q", "U", "V"][: len(out_ds.stokes)]

        for k, v in output.d_radiance.items():
            mapping = atmosphere.storage.get_derivative_mapping(k)

            name = k if mapping.assign_name == "" else mapping.assign_name

            if name in out_ds:
                out_ds[name] += v
            else:
                out_ds[name] = xr.DataArray(
                    v,
                    dims=[mapping.interp_dim, "wavelength", "los", "stokes"],
                )
        for k, v in output.d_radiance_surf.items():
            mapping = atmosphere.surface.get_derivative_mapping(k)

            mapped_derivative = map_surface_derivative(
                mapping, v, ["wavelength", "los", "stokes"]
            )
            if mapping.interp_dim == "dummy":
                mapped_derivative = mapped_derivative.isel(**{mapping.interp_dim: 0})
            out_ds[k] = mapped_derivative

        for k, v in output.d_flux.items():
            mapping = atmosphere.storage.get_derivative_mapping(k)

            base_name = k if mapping.assign_name == "" else mapping.assign_name

            for i, flux_type in enumerate(flux_types):
                name = f"{base_name}_{flux_type}_flux"

                if name in out_ds:
                    out_ds[name] += v[:, i]
                else:
                    out_ds[name] = xr.DataArray(
                        v[:, i],
                        dims=[mapping.interp_dim, "wavelength", "flux_location"],
                    )

        for k, v in output.d_flux_surf.items():
            mapping = atmosphere.surface.get_derivative_mapping(k)

            base_name = k

            for i, flux_type in enumerate(flux_types):
                name = f"{base_name}_{flux_type}_flux"

                mapped_derivative = map_surface_derivative(
                    mapping, v[i], ["wavelength", "flux_location"]
                )
                if mapping.interp_dim == "dummy":
                    mapped_derivative = mapped_derivative.isel(
                        **{mapping.interp_dim: 0}
                    )
                out_ds[name] = mapped_derivative

        if isinstance(self._viewing_geometry, ViewingGeometryContainer):
            out_ds = self._viewing_geometry.add_geometry_to_radiance(out_ds)

        if self._config.output_los_optical_depth:
            los_od = output.los_optical_depth
            out_ds["los_optical_depth"] = xr.DataArray(
                los_od,
                dims=["wavelength", "los"],
            )

        return out_ds
