from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
import xarray as xr


@dataclass(frozen=True)
class _ParameterSpec:
    """Internal description of a parameter axis in a Jacobian block."""

    parameter_dims: tuple[str, ...]
    diagonal_dims: tuple[str, ...] = ()


def _semantic_parameter_name(weighting_function_name: str) -> str:
    if weighting_function_name.startswith("wf_"):
        return weighting_function_name[3:]
    return weighting_function_name


def _dimension_coordinate(array: xr.DataArray, dim: str) -> xr.DataArray | None:
    if dim not in array.coords:
        return None

    coordinate = array.coords[dim]
    if coordinate.dims != (dim,):
        return None
    return coordinate


def _validated_array(
    array: xr.DataArray,
    reference: xr.DataArray,
    expected_dims: tuple[str, ...],
    description: str,
) -> xr.DataArray:
    if not isinstance(array, xr.DataArray):
        msg = f"{description} must be an xarray.DataArray"
        raise TypeError(msg)

    if len(array.dims) != len(expected_dims) or set(array.dims) != set(expected_dims):
        msg = f"{description} must have dimensions {expected_dims}; got {array.dims}"
        raise ValueError(msg)

    array = array.transpose(*expected_dims)
    for dim in expected_dims:
        if array.sizes[dim] != reference.sizes[dim]:
            msg = (
                f"{description} dimension {dim!r} must have size "
                f"{reference.sizes[dim]}; got {array.sizes[dim]}"
            )
            raise ValueError(msg)

        expected_coordinate = _dimension_coordinate(reference, dim)
        if expected_coordinate is None:
            continue

        actual_coordinate = _dimension_coordinate(array, dim)
        if actual_coordinate is None:
            msg = f"{description} is missing the coordinate for dimension {dim!r}"
            raise ValueError(msg)
        if not actual_coordinate.equals(expected_coordinate):
            msg = f"{description} coordinate for dimension {dim!r} does not match"
            raise ValueError(msg)

    if not np.issubdtype(array.dtype, np.number) or np.issubdtype(
        array.dtype, np.complexfloating
    ):
        msg = f"{description} must contain real numeric values; got {array.dtype}"
        raise ValueError(msg)

    return array


def _parameter_template(block: xr.DataArray, spec: _ParameterSpec) -> xr.DataArray:
    """Construct a zero tangent carrying only a parameter's input grid."""
    coordinates = {}
    for dim in spec.parameter_dims:
        coordinate = _dimension_coordinate(block, dim)
        if coordinate is not None:
            coordinates[dim] = coordinate

    shape = tuple(block.sizes[dim] for dim in spec.parameter_dims)
    return xr.DataArray(
        np.zeros(shape, dtype=block.dtype),
        dims=spec.parameter_dims,
        coords=coordinates,
    )


class _FullJacobianBackend:
    """Linearization backend that contracts a cached structured Jacobian."""

    def __init__(
        self,
        value: xr.DataArray,
        jacobian: xr.Dataset,
        specs: Mapping[str, _ParameterSpec],
        tangent_template: xr.Dataset,
    ) -> None:
        self.value = value
        self.jacobian = jacobian
        self.specs = MappingProxyType(dict(specs))
        self.tangent_template = tangent_template

    def jvp(self, tangent: xr.Dataset) -> xr.DataArray:
        if not isinstance(tangent, xr.Dataset):
            msg = "The JVP tangent must be an xarray.Dataset"
            raise TypeError(msg)

        unknown = set(tangent.data_vars) - set(self.specs)
        if unknown:
            names = ", ".join(sorted(unknown))
            msg = f"Unknown tangent parameter(s): {names}"
            raise ValueError(msg)

        result = xr.zeros_like(self.value)
        for name in tangent.data_vars:
            spec = self.specs[name]
            block = self.jacobian[name]
            direction = _validated_array(
                tangent[name],
                self.tangent_template[name],
                spec.parameter_dims,
                f"Tangent parameter {name!r}",
            )

            contribution = block * direction
            contraction_dims = tuple(
                dim for dim in spec.parameter_dims if dim not in spec.diagonal_dims
            )
            if contraction_dims:
                contribution = contribution.sum(dim=contraction_dims)

            contribution = contribution.transpose(*self.value.dims)
            result = result + contribution

        return result

    def vjp(self, cotangent: xr.DataArray) -> xr.Dataset:
        cotangent = _validated_array(
            cotangent,
            self.value,
            tuple(self.value.dims),
            "The VJP radiance cotangent",
        )

        gradients: dict[str, xr.DataArray] = {}
        output_dims = tuple(self.value.dims)
        for name, spec in self.specs.items():
            weighted = self.jacobian[name] * cotangent
            contraction_dims = tuple(
                dim for dim in output_dims if dim not in spec.diagonal_dims
            )
            if contraction_dims:
                weighted = weighted.sum(dim=contraction_dims)
            gradients[name] = weighted.transpose(*spec.parameter_dims)

        return xr.Dataset(gradients)


class StaleLinearizationError(RuntimeError):
    """Raised when an operation requires an atmosphere that has changed."""


class _EngineBackend:
    """Streaming engine backend with lazy full-Jacobian materialization."""

    def __init__(
        self,
        value: xr.DataArray,
        specs: Mapping[str, _ParameterSpec],
        tangent_template: xr.Dataset,
        atmosphere,
        revision: int,
        jacobian_loader: Callable[[], xr.Dataset],
        jvp_evaluator: Callable[[xr.Dataset], xr.DataArray],
        vjp_evaluator: Callable[[xr.DataArray], xr.Dataset],
        owner,
    ) -> None:
        self.value = value
        self.specs = MappingProxyType(dict(specs))
        self.tangent_template = tangent_template
        self._atmosphere = atmosphere
        self._revision = revision
        self._jacobian_loader = jacobian_loader
        self._jvp_evaluator = jvp_evaluator
        self._vjp_evaluator = vjp_evaluator
        self._owner = owner
        self._jacobian: xr.Dataset | None = None

    def _require_current(self) -> None:
        if self._atmosphere.revision != self._revision:
            msg = (
                "The atmosphere has changed since this linearization was "
                "created. Call engine.linearize(atmosphere) again."
            )
            raise StaleLinearizationError(msg)

    @property
    def jacobian(self) -> xr.Dataset:
        if self._jacobian is None:
            self._require_current()
            self._jacobian = self._jacobian_loader()
        return self._jacobian

    def jvp(self, tangent: xr.Dataset) -> xr.DataArray:
        if not isinstance(tangent, xr.Dataset):
            msg = "The JVP tangent must be an xarray.Dataset"
            raise TypeError(msg)

        unknown = set(tangent.data_vars) - set(self.specs)
        if unknown:
            names = ", ".join(sorted(unknown))
            msg = f"Unknown tangent parameter(s): {names}"
            raise ValueError(msg)

        normalized: dict[str, xr.DataArray] = {}
        for name in tangent.data_vars:
            normalized[name] = _validated_array(
                tangent[name],
                self.tangent_template[name],
                self.specs[name].parameter_dims,
                f"Tangent parameter {name!r}",
            )

        self._require_current()
        return self._jvp_evaluator(xr.Dataset(normalized))

    def vjp(self, cotangent: xr.DataArray) -> xr.Dataset:
        cotangent = _validated_array(
            cotangent,
            self.value,
            tuple(self.value.dims),
            "The VJP radiance cotangent",
        )
        self._require_current()
        return self._vjp_evaluator(cotangent)


class Linearization:
    """Radiance and its local linear model at a built atmosphere state.

    Engine-created instances stream Jacobian-vector and vector-Jacobian
    contractions through the radiative-transfer output path. The complete
    structured Jacobian is materialized only when requested.
    """

    def __init__(self, backend: _FullJacobianBackend | _EngineBackend) -> None:
        self._backend = backend

    @classmethod
    def _from_engine(
        cls,
        value: xr.DataArray,
        specs: Mapping[str, _ParameterSpec],
        tangent_template: xr.Dataset,
        atmosphere,
        revision: int,
        jacobian_loader: Callable[[], xr.Dataset],
        jvp_evaluator: Callable[[xr.Dataset], xr.DataArray],
        vjp_evaluator: Callable[[xr.DataArray], xr.Dataset],
        owner,
    ) -> Linearization:
        if not specs:
            msg = (
                "No radiance derivative mappings are available. Construct the "
                "atmosphere with calculate_derivatives=True and register at "
                "least one differentiable parameter."
            )
            raise ValueError(msg)
        return cls(
            _EngineBackend(
                value,
                specs,
                tangent_template,
                atmosphere,
                revision,
                jacobian_loader,
                jvp_evaluator,
                vjp_evaluator,
                owner,
            )
        )

    @classmethod
    def _from_result(
        cls,
        result: xr.Dataset,
        weighting_function_specs: Mapping[str, _ParameterSpec],
    ) -> Linearization:
        jacobian_variables: dict[str, xr.DataArray] = {}
        parameter_specs: dict[str, _ParameterSpec] = {}
        tangent_variables: dict[str, xr.DataArray] = {}

        for weighting_function_name, spec in weighting_function_specs.items():
            parameter_name = _semantic_parameter_name(weighting_function_name)
            if parameter_name in jacobian_variables:
                msg = (
                    "Multiple weighting functions map to the semantic parameter "
                    f"name {parameter_name!r}"
                )
                raise ValueError(msg)
            jacobian_variables[parameter_name] = result[weighting_function_name]
            parameter_specs[parameter_name] = spec
            tangent_variables[parameter_name] = _parameter_template(
                result[weighting_function_name], spec
            )

        if not jacobian_variables:
            msg = (
                "No radiance derivative mappings are available. Construct the "
                "atmosphere with calculate_derivatives=True and register at "
                "least one differentiable parameter."
            )
            raise ValueError(msg)

        backend = _FullJacobianBackend(
            result["radiance"],
            xr.Dataset(jacobian_variables),
            parameter_specs,
            xr.Dataset(tangent_variables),
        )
        return cls(backend)

    @property
    def value(self) -> xr.DataArray:
        """Radiance at the linearization point."""
        return self._backend.value

    @property
    def jacobian(self) -> xr.Dataset:
        """Structured radiance Jacobian, keyed by semantic parameter name."""
        return self._backend.jacobian

    @property
    def parameters(self) -> tuple[str, ...]:
        """Parameters accepted by jvp and returned by vjp."""
        return tuple(self._backend.specs)

    @property
    def parameter_dims(self) -> Mapping[str, tuple[str, ...]]:
        """Expected tangent dimensions for every parameter."""
        return MappingProxyType(
            {name: spec.parameter_dims for name, spec in self._backend.specs.items()}
        )

    @property
    def tangent_template(self) -> xr.Dataset:
        """Zero tangent containing every parameter's input grid.

        The returned dataset describes parameter shapes, dimensions, and
        dimension coordinates without requiring access to :attr:`jacobian`.
        It is an independent copy that callers may fill or subset before
        passing it to :meth:`jvp`.
        """
        return self._backend.tangent_template.copy(deep=True)

    def jvp(self, tangent: xr.Dataset) -> xr.DataArray:
        """Apply the radiance Jacobian to one labeled tangent direction."""
        return self._backend.jvp(tangent)

    def vjp(self, cotangent: xr.DataArray) -> xr.Dataset:
        """Apply the transpose radiance Jacobian to one radiance cotangent."""
        return self._backend.vjp(cotangent)
