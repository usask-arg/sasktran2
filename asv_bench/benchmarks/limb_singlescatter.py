from __future__ import annotations

import numpy as np
import sasktran2 as sk
import xarray as xr

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


@parameterized(["batch_size", "calc_deriv"], ([1, 4, 8, 16], [False, True]))
class LimbSingleScatterWavelengthBatch:
    """Compare opt-in wavelength batch widths on the limb workload."""

    def setup(self, batch_size, calc_deriv):
        cos_sza = 0.5

        config = sk.Config()
        config.num_threads = 1
        config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
        config.num_singlescatter_moments = 4
        config.num_stokes = 1
        config.wavelength_batch_size = batch_size

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
        self._atmosphere.storage.total_extinction[:] = 0.5
        self._atmosphere.storage.ssa[:] = 1
        self._atmosphere.leg_coeff.a1[0, :, :] = 1
        self._atmosphere.leg_coeff.a1[2, :, :] = 0.5
        self._atmosphere.surface.albedo[:] = 0.0
        self._engine = sk.Engine(config, model_geometry, viewing_geo)

    def time_limb_single_scatter_batch(self, batch_size, calc_deriv):  # noqa: ARG002
        self._engine.calculate_radiance(self._atmosphere)


@parameterized(["batch_size"], ([1, 32],))
class LimbSingleScatterLinearizationProducts:
    """Track exact-single-scatter product and full-Jacobian execution."""

    def setup(self, batch_size):
        config = sk.Config()
        config.num_threads = 1
        config.single_scatter_source = sk.SingleScatterSource.Exact
        config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
        config.num_singlescatter_moments = 4
        config.wavelength_batch_size = batch_size

        geometry = sk.Geometry1D(
            cos_sza=0.5,
            solar_azimuth=0.0,
            earth_radius_m=6_372_000.0,
            altitude_grid_m=np.linspace(0.0, 80_000.0, 41),
            interpolation_method=sk.InterpolationMethod.LinearInterpolation,
            geometry_type=sk.GeometryType.Spherical,
        )
        viewing = sk.ViewingGeometry()
        for altitude_m in np.linspace(0.0, 40_000.0, 20):
            viewing.add_ray(sk.TangentAltitudeSolar(altitude_m, 0.0, 200_000.0, 0.5))

        def atmosphere():
            result = sk.Atmosphere(
                geometry,
                config,
                calculate_derivatives=True,
                numwavel=256,
                pressure_derivative=False,
                temperature_derivative=False,
                specific_humidity_derivative=False,
            )
            result.storage.total_extinction[:] = 1.0e-5
            result.storage.ssa[:] = 0.8
            result.leg_coeff.a1[0] = 1.0
            result.leg_coeff.a1[2] = 0.3
            return result

        self._product_atmosphere = atmosphere()
        self._product_engine = sk.Engine(config, geometry, viewing)
        self._linearization = self._product_engine.linearize(self._product_atmosphere)
        self._tangent = self._linearization.tangent_template[["extinction"]]
        self._tangent["extinction"].data[:] = 1.0
        self._cotangent = xr.ones_like(self._linearization.value)

        # Full-Jacobian timing mutates an atmosphere revision while rebuilding
        # constituent mappings. Keep it isolated from the reusable product
        # linearization so ASV method order cannot make that object stale.
        self._jacobian_atmosphere = atmosphere()
        self._jacobian_engine = sk.Engine(config, geometry, viewing)

    def time_jvp(self, _batch_size):
        self._linearization.jvp(self._tangent)

    def time_vjp(self, _batch_size):
        self._linearization.vjp(self._cotangent)

    def time_materialized_jacobian(self, _batch_size):
        self._jacobian_engine.calculate_radiance(self._jacobian_atmosphere)
