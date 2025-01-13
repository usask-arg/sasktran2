from __future__ import annotations

import numpy as np
import sasktran2 as sk

from . import parameterized

nwav = 1000


@parameterized(["nlyr", "calc_deriv"], ([2, 20, 100], [False, True]))
class TwoStreamNadirPlaneParallel:
    def setup(self, nlyr, calc_deriv):
        cos_sza = 0.5

        config = sk.Config()
        config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
        config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates

        config.num_streams = 2
        config.num_singlescatter_moments = 4
        config.num_stokes = 1

        model_geometry = sk.Geometry1D(
            cos_sza=cos_sza,
            solar_azimuth=0,
            earth_radius_m=6372000,
            altitude_grid_m=np.linspace(0, 100000, nlyr + 1),
            interpolation_method=sk.InterpolationMethod.LinearInterpolation,
            geometry_type=sk.GeometryType.PlaneParallel,
        )

        viewing_geo = sk.ViewingGeometry()

        viewing_geo.add_ray(sk.GroundViewingSolar(cos_sza, np.deg2rad(30), 0.02, 2.0))
        viewing_geo.add_ray(sk.GroundViewingSolar(cos_sza, np.deg2rad(60), 0.92, 2.0))

        self._atmosphere = sk.Atmosphere(
            model_geometry, config, calculate_derivatives=calc_deriv, numwavel=nwav
        )

        self._atmosphere.storage.total_extinction[:] = 0.5 / 1.0
        self._atmosphere.storage.ssa[:] = 1

        self._atmosphere.leg_coeff.a1[0, :, :] = 1
        self._atmosphere.leg_coeff.a1[2, :, :] = 0.5

        self._atmosphere.surface.albedo[:] = 0.0

        self._engine = sk.Engine(config, model_geometry, viewing_geo)

    def time_two_stream_nadir(self, nlyr, calc_deriv):  # noqa: ARG002
        self._engine.calculate_radiance(self._atmosphere)
