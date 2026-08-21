from __future__ import annotations

import numpy as np
import sasktran2 as sk
import xarray as xr

from . import _skip_slow

EARTH_RADIUS_M = 6_372_000.0
TOP_ALTITUDE_M = 80_000.0

NUM_ATMOSPHERE_HORIZONTAL = 40
NUM_ATMOSPHERE_VERTICAL = 80
NUM_SOURCE_HORIZONTAL = 5
NUM_SOURCE_VERTICAL = 25
NUM_LOS_HORIZONTAL = 40
NUM_LOS_VERTICAL = 60
NUM_ANGULAR = 110
NUM_THREADS = 1

ATMOSPHERE_HORIZONTAL_ANGLES_RAD = np.linspace(-0.4, 0.4, NUM_ATMOSPHERE_HORIZONTAL)
ATMOSPHERE_ALTITUDES_M = np.linspace(0.0, TOP_ALTITUDE_M, NUM_ATMOSPHERE_VERTICAL)
SOURCE_ALTITUDES_M = (np.arange(NUM_SOURCE_VERTICAL, dtype=np.float64) + 0.5) * (
    TOP_ALTITUDE_M / NUM_SOURCE_VERTICAL
)
LOS_HORIZONTAL_ANGLES_RAD = np.linspace(-0.2, 0.2, NUM_LOS_HORIZONTAL)
LOS_TANGENT_ALTITUDES_M = (np.arange(NUM_LOS_VERTICAL, dtype=np.float64) + 0.5) * (
    TOP_ALTITUDE_M / NUM_LOS_VERTICAL
)


class SuccessiveOrders2DRepeatedVJP:
    """Large 2D source-grid workload with repeated atmosphere updates."""

    timeout = 1200
    number = 1
    repeat = 3
    warmup_time = 0

    def setup(self):
        _skip_slow()

        config = sk.Config()
        config.num_threads = NUM_THREADS
        config.log_level = sk.LogLevel.Error
        config.num_stokes = 1
        config.num_streams = 4
        config.num_singlescatter_moments = 4
        config.single_scatter_source = sk.SingleScatterSource.NoSource
        config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
        config.occultation_source = sk.OccultationSource.NoSource
        config.emission_source = sk.EmissionSource.NoSource
        config.num_sza = NUM_SOURCE_HORIZONTAL
        config.successive_orders_altitude_grid_m = SOURCE_ALTITUDES_M
        config.num_successive_orders_incoming = NUM_ANGULAR
        config.num_successive_orders_outgoing = NUM_ANGULAR
        config.delta_m_scaling = False

        geometry = sk.Geometry2D(
            cos_sza=0.6,
            solar_azimuth=0.0,
            earth_radius_m=EARTH_RADIUS_M,
            altitude_grid_m=ATMOSPHERE_ALTITUDES_M,
            horizontal_angle_grid_radians=ATMOSPHERE_HORIZONTAL_ANGLES_RAD,
        )
        viewing = sk.ViewingGeometry()
        for horizontal_angle in LOS_HORIZONTAL_ANGLES_RAD:
            for tangent_altitude_m in LOS_TANGENT_ALTITUDES_M:
                viewing.add_ray(
                    sk.TangentAltitude(
                        tangent_altitude_m=tangent_altitude_m,
                        observer_altitude_m=200_000.0,
                        horizontal_angle_radians=horizontal_angle,
                        viewing_azimuth_radians=0.0,
                    )
                )

        atmosphere = sk.Atmosphere(
            geometry,
            config,
            wavelengths_nm=np.array([600.0]),
            calculate_derivatives=True,
            pressure_derivative=False,
            temperature_derivative=False,
            specific_humidity_derivative=False,
            legendre_derivative=False,
        )
        horizontal, altitude = np.meshgrid(
            ATMOSPHERE_HORIZONTAL_ANGLES_RAD,
            ATMOSPHERE_ALTITUDES_M,
            indexing="ij",
        )
        self._base_extinction = (
            1.5e-5
            * np.exp(-altitude / 18_000.0)
            * (1.0 + 0.2 * np.sin(2.0 * np.pi * horizontal / 0.8))
        ).reshape(-1, 1)
        self._base_ssa = (
            0.92
            - 0.12 * altitude / TOP_ALTITUDE_M
            + 0.02 * np.cos(2.0 * np.pi * horizontal / 0.8)
        ).reshape(-1, 1)
        atmosphere.storage.total_extinction[:] = self._base_extinction
        atmosphere.storage.ssa[:] = self._base_ssa
        atmosphere.leg_coeff.a1[0] = 1.0
        atmosphere.leg_coeff.a1[2] = 0.3
        atmosphere.surface.albedo[:] = 0.1
        atmosphere.mark_changed()

        self._atmosphere = atmosphere
        self._engine = sk.Engine(config, geometry, viewing)
        self._linearization = self._engine.linearize(atmosphere)
        self._cotangent = xr.ones_like(self._linearization.value)
        self._cotangent.data[:] = np.linspace(0.5, 1.5, self._cotangent.size).reshape(
            self._cotangent.shape
        )
        self._update_index = 0

    def _change_atmosphere(self):
        self._update_index += 1
        sign = 1.0 if self._update_index % 2 else -1.0
        self._atmosphere.storage.total_extinction[:] = self._base_extinction * (
            1.0 + sign * 0.01
        )
        self._atmosphere.storage.ssa[:] = np.clip(
            self._base_ssa + sign * 0.002, 0.0, 1.0
        )
        self._atmosphere.mark_changed()

    def time_vjp_repeated_same_atmosphere(self):
        self._linearization.vjp(self._cotangent, parameters=("extinction", "ssa"))

    def time_linearize_changed_atmosphere(self):
        self._change_atmosphere()
        self._linearization = self._engine.linearize(self._atmosphere)

    def time_linearize_and_vjp_changed_atmosphere(self):
        self._change_atmosphere()
        self._linearization = self._engine.linearize(self._atmosphere)
        self._linearization.vjp(self._cotangent, parameters=("extinction", "ssa"))

    def peakmem_retained_linearization(self):
        pass

    def peakmem_vjp_changed_atmosphere(self):
        self._change_atmosphere()
        self._linearization = self._engine.linearize(self._atmosphere)
        self._linearization.vjp(self._cotangent, parameters=("extinction", "ssa"))
