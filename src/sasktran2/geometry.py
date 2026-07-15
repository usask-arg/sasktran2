from __future__ import annotations

import numpy as np

from sasktran2._core_rust import (
    GeometryType,
    InterpolationMethod,
    PyGeometry1D,
    PyGeometry2D,
)


class Geometry1D:
    _geometry: PyGeometry1D

    def __init__(
        self,
        cos_sza: float,
        solar_azimuth: float,
        earth_radius_m: float,
        altitude_grid_m: np.ndarray,
        interpolation_method: InterpolationMethod = InterpolationMethod.LinearInterpolation,
        geometry_type: GeometryType = GeometryType.Spherical,
    ):
        """
        A geometry object that depends only on altitude.

        Parameters
        ----------
        cos_sza : float
            Cosine of the solar zenith angle at the "reference point".  This should be thought of as the point in the atmosphere
            where the radiative transfer calculation is most accurate.
        solar_azimuth : float
            Solar azimuth angle at the "reference point".  Most of the time you want to set this to 0.0 and only use the
            azimuth angle in the line of sight declaration.
        earth_radius_m : float
            Radius of the Earth in meters.  This is used for spherical ray-tracing.
        altitude_grid_m : np.ndarray
            Altitude grid in [m] above the surface of the Earth.
        interpolation_method : InterpolationMethod, optional
            The interpolation method to use inbetween the grid points, by default InterpolationMethod.LinearInterpolation
        geometry_type : GeometryType, optional
            The type of geometry to use in the calculation, i.e. Spherical, PlaneParallel, by default GeometryType.Spherical
        """
        self._geometry = PyGeometry1D(
            cos_sza,
            solar_azimuth,
            earth_radius_m,
            np.atleast_1d(altitude_grid_m).astype(np.float64),
            interpolation_method,
            geometry_type,
        )

    def altitudes(self) -> np.ndarray:
        return self._geometry.altitudes()

    @property
    def refractive_index(self) -> np.ndarray:
        return self._geometry.refractive_index

    @refractive_index.setter
    def refractive_index(self, value: np.ndarray) -> None:
        self._geometry.refractive_index = value


class Geometry2D:
    """A structured spherical geometry varying in altitude and one angle.

    Atmospheric locations are flattened with altitude varying fastest. The
    :attr:`shape` is ordered as ``(horizontal, altitude)``.
    """

    _geometry: PyGeometry2D

    def __init__(
        self,
        cos_sza: float,
        solar_azimuth: float,
        earth_radius_m: float,
        altitude_grid_m: np.ndarray,
        horizontal_angle_grid_radians: np.ndarray,
        interpolation_method: InterpolationMethod = InterpolationMethod.LinearInterpolation,
    ):
        """Construct a structured spherical 2D geometry.

        Parameters
        ----------
        cos_sza
            Cosine of the solar zenith angle at the reference point.
        solar_azimuth
            Solar azimuth at the reference point in radians.
        earth_radius_m
            Spherical Earth radius in metres.
        altitude_grid_m
            Strictly increasing altitude locations in metres.
        horizontal_angle_grid_radians
            Strictly increasing horizontal angles measured from reference-z
            toward reference-x. Their total span must be less than pi.
        interpolation_method
            Altitude interpolation method. Horizontal interpolation is linear.
        """
        altitude_grid_m = np.ascontiguousarray(
            np.atleast_1d(altitude_grid_m), dtype=np.float64
        )
        horizontal_angle_grid_radians = np.ascontiguousarray(
            np.atleast_1d(horizontal_angle_grid_radians), dtype=np.float64
        )
        self._geometry = PyGeometry2D(
            cos_sza,
            solar_azimuth,
            earth_radius_m,
            altitude_grid_m,
            horizontal_angle_grid_radians,
            interpolation_method,
        )

    def altitudes(self) -> np.ndarray:
        """Return a copy of the altitude grid in metres."""
        return self._geometry.altitudes()

    def horizontal_angles(self) -> np.ndarray:
        """Return a copy of the horizontal angle grid in radians."""
        return self._geometry.horizontal_angles()

    @property
    def shape(self) -> tuple[int, int]:
        """Location shape ordered as ``(horizontal, altitude)``."""
        return self._geometry.location_shape()

    def location_index(self, altitude_index: int, horizontal_index: int) -> int:
        """Return the altitude-fastest flattened location index."""
        return self._geometry.location_index(altitude_index, horizontal_index)
