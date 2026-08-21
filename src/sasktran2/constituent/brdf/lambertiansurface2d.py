from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2.atmosphere import Atmosphere
from sasktran2.util.interpolation import linear_interpolating_matrix

from ..base import Constituent
from . import PyLambertian


class LambertianSurface2D(Constituent):
    """Lambertian albedo on an orbital-position surface grid.

    ``albedo`` may be a scalar, a gray ``(orbital_position,)`` field, or a
    ``(orbital_position, spectral_point)`` field. The last form is native on
    the atmosphere wavelength grid unless wavelengths or wavenumbers are
    supplied for interpolation. Each local orbital engine retains the gathered
    horizontal field and linearly interpolates it at a surface intersection.
    """

    def __init__(
        self,
        albedo: np.ndarray | float,
        wavelengths_nm: np.ndarray | None = None,
        wavenumbers_cminv: np.ndarray | None = None,
        out_of_bounds_mode: str = "zero",
    ) -> None:
        super().__init__()
        if wavelengths_nm is not None and wavenumbers_cminv is not None:
            msg = "Specify wavelengths_nm or wavenumbers_cminv, not both"
            raise ValueError(msg)
        albedo = np.asarray(albedo, dtype=np.float64)
        if albedo.ndim > 2:
            msg = (
                "albedo must be scalar, (orbital_position,), or "
                "(orbital_position, spectral_point)"
            )
            raise ValueError(msg)
        if albedo.size == 0 or np.any(~np.isfinite(albedo)):
            msg = "albedo must contain finite values"
            raise ValueError(msg)
        if np.any((albedo < 0) | (albedo > 1)):
            msg = "Lambertian albedo must be between zero and one"
            raise ValueError(msg)
        if albedo.ndim != 2 and (
            wavelengths_nm is not None or wavenumbers_cminv is not None
        ):
            msg = "A supplied spectral grid requires 2D spatial-spectral albedo"
            raise ValueError(msg)
        if out_of_bounds_mode not in ("zero", "extend"):
            msg = "out_of_bounds_mode must be 'zero' or 'extend'"
            raise ValueError(msg)

        self._albedo = albedo.copy()
        self._wavelengths_nm = self._validate_spectral_grid(
            wavelengths_nm, "wavelengths_nm"
        )
        self._wavenumbers_cminv = self._validate_spectral_grid(
            wavenumbers_cminv, "wavenumbers_cminv"
        )
        supplied_grid = (
            self._wavelengths_nm
            if self._wavelengths_nm is not None
            else self._wavenumbers_cminv
        )
        if supplied_grid is not None and supplied_grid.size != albedo.shape[1]:
            msg = (
                "The supplied spectral-grid length must match the albedo "
                f"spectral dimension: {supplied_grid.size} != {albedo.shape[1]}"
            )
            raise ValueError(msg)
        self._out_of_bounds_mode = out_of_bounds_mode
        self._native_field: np.ndarray | None = None
        self._spectral_interpolator: np.ndarray | None = None

    @staticmethod
    def _validate_spectral_grid(
        value: np.ndarray | None, name: str
    ) -> np.ndarray | None:
        if value is None:
            return None
        value = np.asarray(value, dtype=np.float64)
        if (
            value.ndim != 1
            or value.size == 0
            or np.any(~np.isfinite(value))
            or np.any(value <= 0)
            or np.any(np.diff(value) <= 0)
        ):
            msg = (
                f"{name} must be a finite, positive, non-empty, strictly "
                "increasing array"
            )
            raise ValueError(msg)
        return value.copy()

    @property
    def volume_spatial_mode(self) -> str:
        return "native_2d"

    @property
    def albedo(self) -> np.ndarray:
        return self._albedo

    @albedo.setter
    def albedo(self, value: np.ndarray | float) -> None:
        value = np.asarray(value, dtype=np.float64)
        if value.shape != self._albedo.shape:
            msg = f"albedo must retain shape {self._albedo.shape}; got {value.shape}"
            raise ValueError(msg)
        if np.any(~np.isfinite(value)) or np.any((value < 0) | (value > 1)):
            msg = "Lambertian albedo must be finite and between zero and one"
            raise ValueError(msg)
        self._albedo = value.copy()

    def _validate_atmosphere(self, atmo: Atmosphere) -> int:
        if np.any(~np.isfinite(self._albedo)) or np.any(
            (self._albedo < 0) | (self._albedo > 1)
        ):
            msg = "Lambertian albedo must be finite and between zero and one"
            raise ValueError(msg)
        geometry = atmo.model_geometry
        if not isinstance(geometry, sk.OrbitalPlaneGeometry):
            msg = "LambertianSurface2D requires an OrbitalPlaneGeometry atmosphere"
            raise TypeError(msg)
        num_orbital_positions = atmo.volume_shape[0]
        if self._albedo.ndim >= 1 and self._albedo.shape[0] != num_orbital_positions:
            msg = (
                "albedo orbital dimension does not match the atmosphere: "
                f"{self._albedo.shape[0]} != {num_orbital_positions}"
            )
            raise ValueError(msg)
        return num_orbital_positions

    def _spectral_matrix(self, atmo: Atmosphere) -> np.ndarray:
        if self._albedo.ndim < 2:
            return np.ones((atmo.num_wavel, 1))
        num_spectral = self._albedo.shape[1]
        if self._wavelengths_nm is not None:
            if atmo.wavelengths_nm is None:
                msg = (
                    "LambertianSurface2D wavelength interpolation requires an "
                    "Atmosphere constructed with wavelengths_nm or wavenumber_cminv"
                )
                raise ValueError(msg)
            if self._wavelengths_nm.size != num_spectral:
                msg = "wavelengths_nm length must match albedo spectral_point"
                raise ValueError(msg)
            return linear_interpolating_matrix(
                self._wavelengths_nm,
                atmo.wavelengths_nm,
                self._out_of_bounds_mode,
            )
        if self._wavenumbers_cminv is not None:
            if atmo.wavenumbers_cminv is None:
                msg = (
                    "LambertianSurface2D wavenumber interpolation requires an "
                    "Atmosphere constructed with wavelengths_nm or wavenumber_cminv"
                )
                raise ValueError(msg)
            if self._wavenumbers_cminv.size != num_spectral:
                msg = "wavenumbers_cminv length must match albedo spectral_point"
                raise ValueError(msg)
            return linear_interpolating_matrix(
                self._wavenumbers_cminv,
                atmo.wavenumbers_cminv,
                self._out_of_bounds_mode,
            )
        if num_spectral != atmo.num_wavel:
            msg = (
                "The albedo spectral dimension must match the atmosphere wavelength "
                "grid when no interpolation grid is supplied"
            )
            raise ValueError(msg)
        return np.eye(atmo.num_wavel)

    def add_to_atmosphere(self, atmo: Atmosphere) -> None:
        num_orbital_positions = self._validate_atmosphere(atmo)
        atmo.surface.brdf = PyLambertian(atmo.nstokes)
        spectral = self._spectral_matrix(atmo)
        if self._albedo.ndim == 0:
            field = np.full(
                (num_orbital_positions, atmo.num_wavel), self._albedo.item()
            )
        elif self._albedo.ndim == 1:
            field = np.broadcast_to(
                self._albedo[:, np.newaxis],
                (num_orbital_positions, atmo.num_wavel),
            ).copy()
        else:
            field = self._albedo @ spectral.T
        self._native_field = np.ascontiguousarray(field, dtype=np.float64)
        self._spectral_interpolator = spectral
        # The composite engine gathers this field into each persistent group.
        # The homogeneous value keeps the ordinary native surface storage valid.
        atmo.surface.brdf_args[0, :] = self._native_field.mean(axis=0)
        atmo._orbital_lambertian_surface = self._native_field
        atmo._orbital_lambertian_spectral_interpolator = np.ascontiguousarray(
            spectral, dtype=np.float64
        )
        atmo._orbital_lambertian_spatial_parameters = self._albedo.ndim != 0

    def register_derivative(self, atmo: Atmosphere, name: str) -> None:
        num_orbital_positions = self._validate_atmosphere(atmo)
        if self._native_field is None or self._spectral_interpolator is None:
            self.add_to_atmosphere(atmo)

        mapping_name = f"wf_{name}_albedo"
        atmo._orbital_lambertian_derivative_name = mapping_name
        mapping = atmo.surface.get_derivative_mapping(mapping_name)
        mapping.d_brdf[:] = 1.0
        if self._albedo.ndim == 0:
            mapping.interpolator = np.ones((atmo.num_wavel, 1))
            mapping.interp_dim = "dummy"
            return

        if self._albedo.ndim == 1:
            mapping.interpolator = np.ones((atmo.num_wavel, num_orbital_positions))
            dims = ("orbital_position",)
            shape = (num_orbital_positions,)
            coords = atmo.model_geometry._xarray_orbital_coordinates()
        else:
            num_spectral = self._albedo.shape[1]
            interpolator = np.zeros(
                (atmo.num_wavel, num_orbital_positions * num_spectral)
            )
            for orbital_index in range(num_orbital_positions):
                start = orbital_index * num_spectral
                interpolator[:, start : start + num_spectral] = (
                    self._spectral_interpolator
                )
            mapping.interpolator = interpolator
            if self._wavenumbers_cminv is not None:
                spectral_dim = f"{name}_wavenumber"
            elif self._wavelengths_nm is not None or atmo.wavelengths_nm is not None:
                spectral_dim = f"{name}_wavelength"
            else:
                spectral_dim = f"{name}_spectral_point"
            dims = ("orbital_position", spectral_dim)
            shape = (num_orbital_positions, num_spectral)
            if self._wavelengths_nm is not None:
                spectral_coord = self._wavelengths_nm
            elif self._wavenumbers_cminv is not None:
                spectral_coord = self._wavenumbers_cminv
            elif atmo.wavelengths_nm is not None:
                spectral_coord = atmo.wavelengths_nm
            else:
                spectral_coord = np.arange(num_spectral)
            coords = {
                **atmo.model_geometry._xarray_orbital_coordinates(),
                spectral_dim: spectral_coord,
            }
        mapping.interp_dim = f"{name}_orbital_surface"
        atmo._surface_derivative_output_layouts[mapping_name] = (
            dims,
            shape,
            coords,
        )
