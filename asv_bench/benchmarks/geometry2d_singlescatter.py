from __future__ import annotations

import numpy as np
import sasktran2 as sk

from . import parameterized


@parameterized(
    ["geometry_kind", "calculate_derivatives"],
    (["1d", "2d"], [False, True]),
)
class GeometrySingleScatterComparison:
    """Track the calculation-only cost of equivalent 1D and 2D scenes."""

    def setup(self, geometry_kind, calculate_derivatives):
        config = sk.Config()
        config.num_threads = 1
        config.single_scatter_source = sk.SingleScatterSource.Exact
        config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
        config.occultation_source = sk.OccultationSource.NoSource
        config.emission_source = sk.EmissionSource.NoSource
        config.num_singlescatter_moments = 6

        altitudes_m = np.linspace(0.0, 100_000.0, 31)
        if geometry_kind == "1d":
            geometry = sk.Geometry1D(
                cos_sza=0.6,
                solar_azimuth=0.0,
                earth_radius_m=6_372_000.0,
                altitude_grid_m=altitudes_m,
                interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                geometry_type=sk.GeometryType.Spherical,
            )
        else:
            geometry = sk.Geometry2D(
                cos_sza=0.6,
                solar_azimuth=0.0,
                earth_radius_m=6_372_000.0,
                altitude_grid_m=altitudes_m,
                horizontal_angle_grid_radians=np.linspace(-0.5, 0.5, 41),
            )

        viewing = sk.ViewingGeometry()
        viewing.add_ray(
            sk.TangentAltitudeSolar(
                tangent_altitude_m=20_000.0,
                relative_azimuth=0.0,
                observer_altitude_m=200_000.0,
                cos_sza=0.6,
            )
        )

        atmosphere = sk.Atmosphere(
            geometry,
            config,
            wavelengths_nm=np.array([600.0]),
            calculate_derivatives=calculate_derivatives,
            pressure_derivative=False,
            temperature_derivative=False,
            specific_humidity_derivative=False,
            legendre_derivative=False,
        )
        extinction = 1.0e-5 * np.exp(-altitudes_m / 30_000.0)
        if geometry_kind == "2d":
            extinction = np.tile(extinction, geometry.shape[0])
        atmosphere.storage.total_extinction[:, 0] = extinction
        atmosphere.storage.ssa[:, 0] = 0.7
        atmosphere.leg_coeff.a1[0, :, 0] = 1.0
        atmosphere.leg_coeff.a1[2, :, 0] = 0.3

        self._atmosphere = atmosphere
        self._engine = sk.Engine(config, geometry, viewing)

    def time_calculate_radiance(self, _geometry_kind, _calculate_derivatives):
        self._engine.calculate_radiance(self._atmosphere)
