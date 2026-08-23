from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass

import numpy as np
import xarray as xr

import sasktran2 as sk
from sasktran2._core_rust import PyEngine
from sasktran2.linearization import (
    Linearization,
    LinearizationBackend,
    _ParameterSpec,
    _selected_parameters,
    _semantic_parameter_name,
)
from sasktran2.viewinggeo.base import ViewingGeometryContainer

_ORBITAL_MAPPED_SURFACE_PREFIX = "__orbital_mapped_surface__"


def map_surface_derivative(
    mapping,
    np_deriv: np.ndarray,
    dims: list[str],
    structured_layout: (
        tuple[tuple[str, ...], tuple[int, ...], dict[str, np.ndarray]] | None
    ) = None,
) -> xr.DataArray:
    if mapping.interpolator is None or len(mapping.interpolator) == 0:
        return xr.DataArray(np_deriv, dims=dims)
    mapped = np.einsum(
        "ij..., il->lij...",
        np_deriv,
        mapping.interpolator,
        optimize=True,
    )
    if structured_layout is None:
        return xr.DataArray(mapped, dims=[mapping.interp_dim, *dims])
    parameter_dims, parameter_shape, coords = structured_layout
    return xr.DataArray(
        mapped.reshape((*parameter_shape, *mapped.shape[1:])),
        dims=[*parameter_dims, *dims],
        coords=coords,
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
    horizontal_dimension = getattr(
        atmosphere.model_geometry, "horizontal_dimension", "horizontal_angle"
    )
    return xr.DataArray(
        derivative.reshape((*output_shape, *derivative.shape[1:])),
        dims=[horizontal_dimension, "altitude", *trailing_dims],
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
                    sk.SingleScatterSource.Table,
                )
                or config.multiple_scatter_source
                not in (
                    sk.MultipleScatterSource.NoSource,
                    sk.MultipleScatterSource.SuccessiveOrders,
                )
                or config.emission_source
                not in (
                    sk.EmissionSource.NoSource,
                    sk.EmissionSource.Standard,
                    sk.EmissionSource.VolumeEmissionRate,
                )
            ):
                msg = (
                    "Geometry2D Engine currently supports exact and table single "
                    "scattering, successive-orders multiple scattering, "
                    "occultation, standard emission, and volume emission rate "
                    "sources"
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
            if (
                config.multiple_scatter_source
                == sk.MultipleScatterSource.SuccessiveOrders
                and config.multiple_scatter_refraction
            ):
                msg = (
                    "Geometry2D successive orders does not support diffuse-ray "
                    "refraction"
                )
                raise NotImplementedError(msg)
            horizontal_source_angles = (
                config.successive_orders_horizontal_angle_grid_radians
            )
            if (
                config.multiple_scatter_source
                == sk.MultipleScatterSource.SuccessiveOrders
                and horizontal_source_angles is not None
            ):
                horizontal_grid = geometry.horizontal_angles()
                outside = (horizontal_source_angles < horizontal_grid[0]) | (
                    horizontal_source_angles > horizontal_grid[-1]
                )
                if np.any(outside):
                    msg = (
                        "successive_orders_horizontal_angle_grid_radians must lie "
                        "inside the Geometry2D horizontal angle range "
                        f"[{horizontal_grid[0]}, {horizontal_grid[-1]}]"
                    )
                    raise ValueError(msg)
        elif (
            isinstance(geometry, sk.Geometry1D)
            and config.multiple_scatter_source
            == sk.MultipleScatterSource.SuccessiveOrders
            and config.successive_orders_horizontal_angle_grid_radians is not None
        ):
            msg = (
                "successive_orders_horizontal_angle_grid_radians is supported "
                "only with Geometry2D"
            )
            raise ValueError(msg)

        self._engine = PyEngine(
            config._config, geometry._geometry, viewing_geometry._viewing_geometry
        )
        self._config = config
        self._geometry = geometry
        self._viewing_geometry = viewing_geometry

    def calculate_radiance(
        self,
        atmosphere: sk.Atmosphere,
        *,
        derivatives: bool | None = None,
    ) -> xr.Dataset:
        """
        Performs the radiative transfer calculation for a given atmosphere

        Parameters
        ----------
        atmosphere : sk.Atmosphere
            The atmosphere object containing the atmospheric profile and constituents
        derivatives : bool, optional
            Whether to calculate and return the registered radiance derivatives.
            The default follows ``atmosphere.calculate_derivatives``. Set this
            to ``False`` for a radiance-only calculation without derivative
            output allocations while retaining a derivative-enabled atmosphere
            for later JVP or VJP calls.

        Returns
        -------
        xr.Dataset
            An xarray dataset containing the radiance and derivatives
        """
        if derivatives is None:
            derivatives = atmosphere.calculate_derivatives
        elif not isinstance(derivatives, bool | np.bool_):
            msg = "derivatives must be a boolean or None"
            raise TypeError(msg)
        result, _ = self._calculate_radiance(
            atmosphere, include_derivatives=bool(derivatives)
        )
        return result

    def linearize(
        self,
        atmosphere: sk.Atmosphere,
        *,
        prepare_parameters: Iterable[str] | None = None,
    ) -> Linearization:
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
        prepare_parameters : iterable of str, optional
            Parameters expected in the first derivative product. Supplying
            this hint lets composite engines prepare and retain only those
            derivative workspaces while calculating the radiance value, so a
            following JVP or VJP can reuse source terms. It does not restrict
            the parameters exposed by the returned linearization. By default,
            derivative workspaces remain lazy.

        Returns
        -------
        Linearization
            The radiance and local derivative operations.
        """
        return self._linearize(atmosphere, prepare_parameters=prepare_parameters)

    def _linearize(
        self,
        atmosphere: sk.Atmosphere,
        *,
        internal_atmosphere=None,
        validate_session: Callable[[], None] | None = None,
        prepare_parameters: Iterable[str] | None = None,
    ) -> Linearization:
        """Construct a linearization from an optional materialized atmosphere."""
        self._validate_atmosphere_geometry(atmosphere)
        if not atmosphere.calculate_derivatives:
            msg = (
                "Engine.linearize requires an atmosphere constructed with "
                "calculate_derivatives=True"
            )
            raise ValueError(msg)
        derivative_backends = {
            mode: self._engine._linearization_backend(mode_index)
            for mode, mode_index in (("jvp", 1), ("vjp", 2))
        }
        if any(backend == 0 for backend in derivative_backends.values()):
            msg = "The configured engine does not support JVP/VJP linearization"
            raise NotImplementedError(msg)
        jacobian_supported = self._engine._supports_linearization(0)

        native_atmosphere = (
            atmosphere.internal_object()
            if internal_atmosphere is None
            else internal_atmosphere
        )
        registry = self._linearization_registry(atmosphere)
        prepared_parameters = (
            _selected_parameters(
                prepare_parameters,
                registry.specs,
                operation="linearization preparation",
            )
            if prepare_parameters is not None
            else ()
        )
        prepared_volume_mappings = list(
            dict.fromkeys(
                name
                for parameter in prepared_parameters
                for name in registry.volume_names.get(parameter, ())
            )
        )
        prepared_surface_mappings = list(
            dict.fromkeys(
                name
                for parameter in prepared_parameters
                for name in registry.surface_names.get(parameter, ())
            )
        )
        revision = atmosphere.revision
        prepared_volume_tangents = {
            name: np.zeros(registry.volume_sizes[name], dtype=np.float64)
            for name in prepared_volume_mappings
        }
        prepared_surface_tangents = {
            name: np.zeros(registry.surface_sizes[name], dtype=np.float64)
            for name in prepared_surface_mappings
        }
        initial_output = self._engine._calculate_jvp(
            native_atmosphere,
            prepared_volume_tangents,
            prepared_surface_tangents,
        )
        value = self._radiance_dataarray(initial_output.radiance, atmosphere)
        backend_names = {
            1: LinearizationBackend.StreamingJacobian,
            2: LinearizationBackend.Native,
        }
        backends = {
            mode: backend_names[backend]
            for mode, backend in derivative_backends.items()
        }

        def load_jacobian() -> xr.Dataset:
            if validate_session is not None:
                validate_session()
            if not jacobian_supported:
                msg = (
                    "The configured engine supports derivative products but "
                    "cannot materialize the full Jacobian"
                )
                raise NotImplementedError(msg)
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
            if validate_session is not None:
                validate_session()
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
            if validate_session is not None:
                validate_session()
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
                elif dim == "orbital_position" and isinstance(
                    atmosphere.model_geometry, sk.Geometry2D
                ):
                    coordinate = getattr(
                        atmosphere.model_geometry,
                        "orbital_positions",
                        np.arange(size),
                    )
                    if len(coordinate) == size:
                        coords[dim] = coordinate
                        orbital_coordinates = getattr(
                            atmosphere.model_geometry,
                            "_xarray_orbital_coordinates",
                            None,
                        )
                        if callable(orbital_coordinates):
                            coords.update(orbital_coordinates())
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
                dims = (
                    getattr(
                        atmosphere.model_geometry,
                        "horizontal_dimension",
                        "horizontal_angle",
                    ),
                    "altitude",
                )
                shape = tuple(structured_shape)
            elif mapping.interp_dim == "dummy" and size == 1:
                dims = ()
                shape = ()
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
            structured_layout = atmosphere.surface_derivative_output_layout(
                internal_name
            )
            if structured_layout is not None:
                dims, shape, coords = structured_layout
                size = int(np.prod(shape))
                spec = _ParameterSpec(dims)
            elif interpolator.size:
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
                (
                    xr.DataArray(np.zeros(shape), dims=dims, coords=coords)
                    if structured_layout is not None
                    else template_for(dims, shape)
                ),
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
        *,
        include_derivatives: bool = True,
    ) -> tuple[xr.Dataset, dict[str, _ParameterSpec]]:
        self._validate_atmosphere_geometry(atmosphere)

        if internal_atmosphere is None:
            internal_atmosphere = atmosphere.internal_object()
        output = (
            self._engine.calculate_radiance(internal_atmosphere)
            if include_derivatives
            else self._engine._calculate_radiance_only(internal_atmosphere)
        )

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
            if k.startswith(_ORBITAL_MAPPED_SURFACE_PREFIX):
                name = k.removeprefix(_ORBITAL_MAPPED_SURFACE_PREFIX)
                mapping = atmosphere.surface.get_derivative_mapping(name)
                structured_layout = atmosphere.surface_derivative_output_layout(name)
                if structured_layout is None:
                    mapped_derivative = xr.DataArray(
                        v,
                        dims=[mapping.interp_dim, "wavelength", "los", "stokes"],
                    )
                    if mapping.interp_dim == "dummy":
                        mapped_derivative = mapped_derivative.isel(dummy=0)
                    spec = _ParameterSpec(tuple(mapped_derivative.dims[:-3]))
                else:
                    parameter_dims, parameter_shape, coords = structured_layout
                    mapped_derivative = xr.DataArray(
                        v.reshape((*parameter_shape, *v.shape[1:])),
                        dims=[
                            *parameter_dims,
                            "wavelength",
                            "los",
                            "stokes",
                        ],
                        coords=coords,
                    )
                    spec = _ParameterSpec(parameter_dims)
                out_ds[name] = mapped_derivative
                radiance_derivative_specs[name] = spec
                continue
            mapping = atmosphere.storage.get_derivative_mapping(k)

            name = k if mapping.assign_name == "" else mapping.assign_name

            mapped_derivative = atmosphere_derivative_dataarray(
                atmosphere,
                k,
                mapping,
                v,
                ["wavelength", "los", "stokes"],
            )
            if mapping.interp_dim == "dummy" and "dummy" in mapped_derivative.dims:
                mapped_derivative = mapped_derivative.isel(dummy=0, drop=True)

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
                mapping,
                v,
                ["wavelength", "los", "stokes"],
                atmosphere.surface_derivative_output_layout(k),
            )
            if mapping.interp_dim == "dummy":
                mapped_derivative = mapped_derivative.isel(**{mapping.interp_dim: 0})
            out_ds[k] = mapped_derivative

            structured_layout = atmosphere.surface_derivative_output_layout(k)
            if structured_layout is not None:
                spec = _ParameterSpec(structured_layout[0])
            elif mapping.interpolator is None or len(mapping.interpolator) == 0:
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
                    mapping,
                    v[i],
                    ["wavelength", "flux_location"],
                    atmosphere.surface_derivative_output_layout(k),
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
            if "orbital_position" in out_ds.dims:
                out_ds.coords["orbital_position"] = getattr(
                    atmosphere.model_geometry,
                    "orbital_positions",
                    np.arange(out_ds.sizes["orbital_position"]),
                )
                orbital_coordinates = getattr(
                    atmosphere.model_geometry,
                    "_xarray_orbital_coordinates",
                    None,
                )
                if callable(orbital_coordinates):
                    out_ds = out_ds.assign_coords(orbital_coordinates())

        if self._config.output_los_optical_depth:
            los_od = output.los_optical_depth
            out_ds["los_optical_depth"] = xr.DataArray(
                los_od,
                dims=["wavelength", "los"],
            )

        return out_ds, radiance_derivative_specs
