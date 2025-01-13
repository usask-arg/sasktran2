from __future__ import annotations

import numpy as np
import sasktran2 as sk

from . import _skip_slow, parameterized

nwav = 1000
nlyr = 100
nlos = 50


@parameterized(["nmoments", "calc_deriv"], ([4, 16], [False, True]))
class LimbSingleScatter:
    def setup(self, nmoments, calc_deriv):
        if nmoments == 16 and calc_deriv:
            _skip_slow()

        cos_sza = 0.5

        config = sk.Config()
        config.multiple_scatter_source = sk.MultipleScatterSource.NoSource

        config.num_singlescatter_moments = nmoments
        config.num_stokes = 1

        model_geometry = sk.Geometry1D(
            cos_sza=cos_sza,
            solar_azimuth=0,
            earth_radius_m=6372000,
            altitude_grid_m=np.linspace(0, 100000, nlyr + 1),
            interpolation_method=sk.InterpolationMethod.LinearInterpolation,
            geometry_type=sk.GeometryType.Spherical,
        )

        viewing_geo = sk.ViewingGeometry()

        for i in range(nlos):
            viewing_geo.add_ray(
                sk.TangentAltitudeSolar(i * (100000 / (nlos + 1)), 0, 200000, cos_sza)
            )

        self._atmosphere = sk.Atmosphere(
            model_geometry, config, calculate_derivatives=calc_deriv, numwavel=nwav
        )

        self._atmosphere.storage.total_extinction[:] = 0.5 / 1.0
        self._atmosphere.storage.ssa[:] = 1

        self._atmosphere.leg_coeff.a1[0, :, :] = 1
        self._atmosphere.leg_coeff.a1[2, :, :] = 0.5

        self._atmosphere.surface.albedo[:] = 0.0

        self._engine = sk.Engine(config, model_geometry, viewing_geo)

    def time_limb_single_scatter(self, nmoments, calc_deriv):  # noqa: ARG002
        self._engine.calculate_radiance(self._atmosphere)
