from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import xarray as xr

import sasktran2 as sk
from sasktran2._core_rust import PyEngine
from sasktran2.linearization import (
    Linearization,
    LinearizationBackend,
    _ParameterSpec,
    _semantic_parameter_name,
)
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


def atmosphere_derivative_dataarray(
    atmosphere: sk.Atmosphere,
    mapping_name: str,
    mapping,
    derivative: np.ndarray,
    trailing_dims: list[str],
) -> xr.DataArray:
    output_shape = atmosphere.derivative_output_shape(mapping_name)
    if output_shape is None:
        return xr.DataArray(
            derivative,
            dims=[mapping.interp_dim, *trailing_dims],
        )

    if len(output_shape) != 2:
        msg = f"Unsupported structured atmosphere derivative shape: {output_shape}"
        raise ValueError(msg)
    return xr.DataArray(
        derivative.reshape((*output_shape, *derivative.shape[1:])),
        dims=["horizontal_angle", "altitude", *trailing_dims],
    )


@dataclass(frozen=True)
class _LinearizationRegistry:
    specs: dict[str, _ParameterSpec]
    tangent_template: xr.Dataset
    volume_names: dict[str, tuple[str, ...]]
    surface_names: dict[str, tuple[str, ...]]
    volume_sizes: dict[str, int]
    surface_sizes: dict[str, int]
    output_names: dict[str, str]
    log_parameters: frozenset[str]


class Engine:
    _engine: PyEngine

    def __init__(
        self,
        config: sk.Config,
        geometry: sk.Geometry1D | sk.Geometry2D,
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
        model_geometry : sk.Geometry1D | sk.Geometry2D
            Geometry for the model
        viewing_geo : sk.ViewingGeometry
            Viewing geometry
        """
        if isinstance(geometry, sk.Geometry2D):
            if (
                config.single_scatter_source
                not in (
                    sk.SingleScatterSource.NoSource,
                    sk.SingleScatterSource.Exact,
                )
                or config.multiple_scatter_source != sk.MultipleScatterSource.NoSource
                or config.emission_source
                not in (
                    sk.EmissionSource.NoSource,
                    sk.EmissionSource.Standard,
                    sk.EmissionSource.VolumeEmissionRate,
                )
            ):
                msg = (
                    "Geometry2D Engine currently supports exact single scattering, "
                    "occultation, standard emission, and volume emission rate "
                    "sources with multiple scattering disabled"
                )
                raise NotImplementedError(msg)
            if viewing_geometry.flux_observers:
                msg = "Geometry2D Engine does not yet support flux observers"
                raise NotImplementedError(msg)
            if config.los_refraction:
                msg = (
                    "Geometry2D Engine does not yet accept per-ray refractive-index "
                    "profiles"
                )
                raise NotImplementedError(msg)

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
        result, _ = self._calculate_radiance(atmosphere)
        return result

    def linearize(self, atmosphere: sk.Atmosphere) -> Linearization:
        """Construct the local linear radiance model for an atmosphere.

        JVP and VJP operations use the most specialized backend supported by
        all active line-of-sight sources. They do not allocate the complete
        structured radiance Jacobian. The Jacobian is calculated and cached
        only if requested.

        Parameters
        ----------
        atmosphere : sk.Atmosphere
            Atmosphere defining the linearization point. It must have been
            constructed with calculate_derivatives=True.

        Returns
        -------
        Linearization
            The radiance and local derivative operations.
        """
        self._validate_atmosphere_geometry(atmosphere)
        if not atmosphere.calculate_derivatives:
            msg = (
                "Engine.linearize requires an atmosphere constructed with "
                "calculate_derivatives=True"
            )
            raise ValueError(msg)
        if not self._engine._supports_linearization(0):
            msg = "The configured engine does not support full-Jacobian linearization"
            raise NotImplementedError(msg)

        native_atmosphere = atmosphere.internal_object()
        revision = atmosphere.revision
        initial_output = self._engine._calculate_jvp(native_atmosphere, {}, {})
        value = self._radiance_dataarray(initial_output.radiance, atmosphere)
        registry = self._linearization_registry(atmosphere)
        backend_names = {
            1: LinearizationBackend.StreamingJacobian,
            2: LinearizationBackend.Native,
        }
        backends = {
            mode: backend_names[self._engine._linearization_backend(mode_index)]
            for mode, mode_index in (("jvp", 1), ("vjp", 2))
        }

        def load_jacobian() -> xr.Dataset:
            result, _ = self._calculate_radiance(
                atmosphere, internal_atmosphere=native_atmosphere
            )
            jacobian: dict[str, xr.DataArray] = {}
            for parameter, output_name in registry.output_names.items():
                block = result[output_name]
                if not registry.specs[parameter].parameter_dims:
                    scalar_dims = tuple(
                        dim for dim in block.dims if dim not in value.dims
                    )
                    for dim in scalar_dims:
                        if block.sizes[dim] != 1:
                            msg = (
                                f"Scalar parameter {parameter!r} has a non-scalar "
                                f"Jacobian dimension {dim!r}"
                            )
                            raise ValueError(msg)
                        block = block.isel({dim: 0}, drop=True)
                if parameter in registry.log_parameters:
                    block = block * result["radiance"]
                jacobian[parameter] = block
            return xr.Dataset(jacobian)

        def evaluate_jvp(tangent: xr.Dataset) -> xr.DataArray:
            volume_tangents: dict[str, np.ndarray] = {}
            surface_tangents: dict[str, np.ndarray] = {}
            for parameter in tangent.data_vars:
                values = np.ascontiguousarray(
                    tangent[parameter].values, dtype=np.float64
                ).reshape(-1)
                for name in registry.volume_names.get(parameter, ()):
                    volume_tangents[name] = values
                for name in registry.surface_names.get(parameter, ()):
                    surface_tangents[name] = values
            output = self._engine._calculate_jvp(
                native_atmosphere, volume_tangents, surface_tangents
            )
            return self._radiance_dataarray(output.jvp, atmosphere)

        def evaluate_vjp(
            cotangent: xr.DataArray, parameters: tuple[str, ...]
        ) -> xr.Dataset:
            if not parameters:
                return xr.Dataset()
            volume_sizes = {
                name: registry.volume_sizes[name]
                for parameter in parameters
                for name in registry.volume_names.get(parameter, ())
            }
            surface_sizes = {
                name: registry.surface_sizes[name]
                for parameter in parameters
                for name in registry.surface_names.get(parameter, ())
            }
            output = self._engine._calculate_vjp(
                native_atmosphere,
                np.ascontiguousarray(cotangent.values, dtype=np.float64),
                volume_sizes,
                surface_sizes,
            )
            gradients = {
                parameter: xr.zeros_like(registry.tangent_template[parameter])
                for parameter in parameters
            }
            for parameter in parameters:
                for name in registry.volume_names.get(parameter, ()):
                    gradients[parameter].data += np.asarray(
                        output.derivative_gradients[name]
                    ).reshape(gradients[parameter].shape)
                for name in registry.surface_names.get(parameter, ()):
                    gradients[parameter].data += np.asarray(
                        output.surface_gradients[name]
                    ).reshape(gradients[parameter].shape)
            return xr.Dataset(gradients)

        return Linearization._from_engine(
            value,
            registry.specs,
            registry.tangent_template,
            atmosphere,
            revision,
            load_jacobian,
            evaluate_jvp,
            evaluate_vjp,
            backends,
            (self, native_atmosphere),
        )

    def _radiance_dataarray(
        self, radiance: np.ndarray, atmosphere: sk.Atmosphere
    ) -> xr.DataArray:
        dataset = xr.Dataset(
            {"radiance": xr.DataArray(radiance, dims=["wavelength", "los", "stokes"])}
        )
        if atmosphere.wavelengths_nm is not None:
            dataset.coords["wavelength"] = atmosphere.wavelengths_nm
        dataset.coords["stokes"] = ["I", "Q", "U", "V"][: len(dataset.stokes)]
        if isinstance(self._viewing_geometry, ViewingGeometryContainer):
            dataset = self._viewing_geometry.add_geometry_to_radiance(dataset)
        return dataset["radiance"]

    def _linearization_registry(
        self, atmosphere: sk.Atmosphere
    ) -> _LinearizationRegistry:
        specs: dict[str, _ParameterSpec] = {}
        templates: dict[str, xr.DataArray] = {}
        volume_names: dict[str, list[str]] = {}
        surface_names: dict[str, list[str]] = {}
        volume_sizes: dict[str, int] = {}
        surface_sizes: dict[str, int] = {}
        output_names: dict[str, str] = {}
        log_parameters: set[str] = set()

        def template_for(dims: tuple[str, ...], shape: tuple[int, ...]) -> xr.DataArray:
            coords: dict[str, np.ndarray] = {}
            for dim, size in zip(dims, shape, strict=True):
                if dim == "wavelength" and atmosphere.wavelengths_nm is not None:
                    coords[dim] = atmosphere.wavelengths_nm
                elif dim == "altitude" and isinstance(
                    atmosphere.model_geometry, sk.Geometry2D
                ):
                    coordinate = atmosphere.model_geometry.altitudes()
                    if len(coordinate) == size:
                        coords[dim] = coordinate
                elif dim == "horizontal_angle" and isinstance(
                    atmosphere.model_geometry, sk.Geometry2D
                ):
                    coordinate = atmosphere.model_geometry.horizontal_angles()
                    if len(coordinate) == size:
                        coords[dim] = coordinate
            return xr.DataArray(np.zeros(shape), dims=dims, coords=coords)

        def register(
            internal_name: str,
            output_name: str,
            spec: _ParameterSpec,
            template: xr.DataArray,
            *,
            surface: bool,
            log_radiance_space: bool = False,
        ) -> None:
            parameter = _semantic_parameter_name(output_name)
            if parameter in specs:
                if specs[parameter] != spec or output_names[parameter] != output_name:
                    msg = (
                        "Multiple derivative mappings for semantic parameter "
                        f"{parameter!r} have incompatible layouts"
                    )
                    raise ValueError(msg)
                if (parameter in log_parameters) != log_radiance_space:
                    msg = (
                        f"Derivative mappings assigned to {parameter!r} mix "
                        "radiance and log-radiance output spaces"
                    )
                    raise ValueError(msg)
            else:
                specs[parameter] = spec
                templates[parameter] = template
                output_names[parameter] = output_name
            if log_radiance_space:
                log_parameters.add(parameter)
            target = surface_names if surface else volume_names
            target.setdefault(parameter, []).append(internal_name)

        for internal_name in atmosphere.storage.derivative_mapping_names():
            mapping = atmosphere.storage.get_derivative_mapping(internal_name)
            interpolator = np.asarray(mapping.interpolator)
            size = (
                int(interpolator.shape[1])
                if interpolator.size
                else atmosphere.num_locations
            )
            structured_shape = atmosphere.derivative_output_shape(internal_name)
            if structured_shape is not None:
                dims = ("horizontal_angle", "altitude")
                shape = tuple(structured_shape)
            else:
                dims = (mapping.interp_dim,)
                shape = (size,)
            output_name = (
                internal_name if mapping.assign_name == "" else mapping.assign_name
            )
            register(
                internal_name,
                output_name,
                _ParameterSpec(dims),
                template_for(dims, shape),
                surface=False,
                log_radiance_space=mapping.log_radiance_space,
            )
            volume_sizes[internal_name] = size

        for internal_name in atmosphere.surface.derivative_mapping_names():
            mapping = atmosphere.surface.get_derivative_mapping(internal_name)
            interpolator = np.asarray(mapping.interpolator)
            if interpolator.size:
                size = int(interpolator.shape[1])
                if mapping.interp_dim == "dummy" or size == 1:
                    dims: tuple[str, ...] = ()
                    shape: tuple[int, ...] = ()
                else:
                    dims = (mapping.interp_dim,)
                    shape = (size,)
                spec = _ParameterSpec(dims)
            else:
                size = atmosphere.num_wavel
                dims = ("wavelength",)
                shape = (size,)
                spec = _ParameterSpec(dims, ("wavelength",))
            register(
                internal_name,
                internal_name,
                spec,
                template_for(dims, shape),
                surface=True,
            )
            surface_sizes[internal_name] = size

        return _LinearizationRegistry(
            specs,
            xr.Dataset(templates),
            {name: tuple(names) for name, names in volume_names.items()},
            {name: tuple(names) for name, names in surface_names.items()},
            volume_sizes,
            surface_sizes,
            output_names,
            frozenset(log_parameters),
        )

    def _validate_atmosphere_geometry(self, atmosphere: sk.Atmosphere) -> None:
        if isinstance(self._geometry, sk.Geometry2D) != isinstance(
            atmosphere.model_geometry, sk.Geometry2D
        ):
            msg = (
                "Engine and atmosphere geometry dimensions do not match: "
                f"{type(self._geometry).__name__} != "
                f"{type(atmosphere.model_geometry).__name__}"
            )
            raise ValueError(msg)
        if isinstance(self._geometry, sk.Geometry2D) and (
            atmosphere.model_geometry is not self._geometry
        ):
            msg = (
                "A Geometry2D atmosphere must use the same Geometry2D object "
                "that was supplied to the Engine"
            )
            raise ValueError(msg)

    def _calculate_radiance(
        self,
        atmosphere: sk.Atmosphere,
        internal_atmosphere=None,
    ) -> tuple[xr.Dataset, dict[str, _ParameterSpec]]:
        self._validate_atmosphere_geometry(atmosphere)

        if internal_atmosphere is None:
            internal_atmosphere = atmosphere.internal_object()
        output = self._engine.calculate_radiance(internal_atmosphere)

        out_ds = xr.Dataset()
        radiance_derivative_specs: dict[str, _ParameterSpec] = {}

        out_ds["radiance"] = xr.DataArray(
            output.radiance,
            dims=["wavelength", "los", "stokes"],
        )

        flux_map = {
            0: "upwelling",
            1: "downwelling",
            2: "actinic",
            3: "divergence",
        }
        flux_types = [flux_map[int(ft)] for ft in self._config.flux_types]

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

            mapped_derivative = atmosphere_derivative_dataarray(
                atmosphere,
                k,
                mapping,
                v,
                ["wavelength", "los", "stokes"],
            )

            if name in out_ds:
                out_ds[name] += mapped_derivative
            else:
                out_ds[name] = mapped_derivative
            spec = _ParameterSpec(tuple(mapped_derivative.dims[:-3]))
            previous_spec = radiance_derivative_specs.setdefault(name, spec)
            if previous_spec != spec:
                msg = (
                    f"Radiance derivative mappings assigned to {name!r} have "
                    "incompatible parameter dimensions"
                )
                raise ValueError(msg)
        for k, v in output.d_radiance_surf.items():
            mapping = atmosphere.surface.get_derivative_mapping(k)

            mapped_derivative = map_surface_derivative(
                mapping, v, ["wavelength", "los", "stokes"]
            )
            if mapping.interp_dim == "dummy":
                mapped_derivative = mapped_derivative.isel(**{mapping.interp_dim: 0})
            out_ds[k] = mapped_derivative

            if mapping.interpolator is None or len(mapping.interpolator) == 0:
                spec = _ParameterSpec(("wavelength",), ("wavelength",))
            else:
                spec = _ParameterSpec(tuple(mapped_derivative.dims[:-3]))
            radiance_derivative_specs[k] = spec

        for k, v in output.d_flux.items():
            mapping = atmosphere.storage.get_derivative_mapping(k)

            base_name = k if mapping.assign_name == "" else mapping.assign_name

            for i, flux_type in enumerate(flux_types):
                name = f"{base_name}_{flux_type}_flux"

                mapped_derivative = atmosphere_derivative_dataarray(
                    atmosphere,
                    k,
                    mapping,
                    v[:, i],
                    ["wavelength", "flux_location"],
                )

                if name in out_ds:
                    out_ds[name] += mapped_derivative
                else:
                    out_ds[name] = mapped_derivative

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

        if isinstance(atmosphere.model_geometry, sk.Geometry2D):
            if "horizontal_angle" in out_ds.dims:
                out_ds.coords["horizontal_angle"] = (
                    atmosphere.model_geometry.horizontal_angles()
                )
            if "altitude" in out_ds.dims:
                out_ds.coords["altitude"] = atmosphere.model_geometry.altitudes()

        if self._config.output_los_optical_depth:
            los_od = output.los_optical_depth
            out_ds["los_optical_depth"] = xr.DataArray(
                los_od,
                dims=["wavelength", "los"],
            )

        return out_ds, radiance_derivative_specs
