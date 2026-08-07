from __future__ import annotations

from typing import ClassVar

import numpy as np
import sasktran2 as sk
import xarray as xr


def _scenario(
    num_stokes: int,
    *,
    anderson_depth: int = 0,
    max_iterations: int = 50,
    relative_tolerance: float = 1.0e-8,
):
    config = sk.Config()
    config.num_threads = 1
    config.num_stokes = num_stokes
    config.num_streams = 8
    config.num_singlescatter_moments = 16
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrdersRust
    config.num_successive_orders_incoming = 350
    config.num_successive_orders_outgoing = 350
    config.successive_orders_max_iterations = max_iterations
    config.successive_orders_relative_tolerance = relative_tolerance
    config.successive_orders_absolute_tolerance = 1.0e-12
    config.successive_orders_anderson_depth = anderson_depth

    geometry = sk.Geometry1D(
        0.6,
        0.0,
        6_372_000.0,
        np.arange(0.0, 30_001.0, 5_000.0),
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    viewing = sk.ViewingGeometry()
    viewing.add_ray(sk.GroundViewingSolar(0.6, 0.2, 0.8, 100_000.0))
    viewing.add_ray(sk.TangentAltitudeSolar(12_500.0, -0.3, 100_000.0, 0.6))
    atmosphere = sk.test_util.scenarios.default_pure_scattering_atmosphere(
        config, geometry, 0.82, albedo=0.3
    )
    return sk.Engine(config, geometry, viewing), atmosphere


class SuccessiveOrdersProducts350:
    params: ClassVar = ([1, 3], [0, 3])
    param_names: ClassVar = ["num_stokes", "anderson_depth"]
    number = 1
    repeat = 5
    timeout = 300

    def setup(self, num_stokes, anderson_depth):
        engine, atmosphere = _scenario(num_stokes, anderson_depth=anderson_depth)
        self.engine = engine
        self.atmosphere = atmosphere
        self.linearization = engine.linearize(atmosphere)

        phase_name = "leg_coeff_2" if num_stokes == 1 else "leg_coeff_8"
        names = ("extinction", "ssa", phase_name, "albedo")
        self.tangent = self.linearization.tangent_template[list(names)]
        self.tangent["extinction"].data[:] = atmosphere.storage.total_extinction[:, 0]
        self.tangent["ssa"].data[:] = 0.1
        self.tangent[phase_name].data[:] = 0.1
        self.tangent["albedo"].data[...] = 0.1
        self.cotangent = xr.ones_like(self.linearization.value)
        self.parameters = names

    def time_jvp(self, num_stokes, anderson_depth):  # noqa: ARG002
        self.linearization.jvp(self.tangent)

    def time_vjp(self, num_stokes, anderson_depth):  # noqa: ARG002
        self.linearization.vjp(self.cotangent, parameters=self.parameters)

    def time_full_jacobian(self, num_stokes, anderson_depth):  # noqa: ARG002
        # Construct a fresh local model because the Jacobian property is
        # intentionally cached after its first materialization.
        _ = self.engine.linearize(self.atmosphere).jacobian


class SuccessiveOrdersProducts350PeakMemory:
    """Full-Jacobian memory without an unrelated retained linearization."""

    params: ClassVar = ([1, 3], [0, 3])
    param_names: ClassVar = ["num_stokes", "anderson_depth"]
    timeout = 300

    def setup(self, num_stokes, anderson_depth):
        self.engine, self.atmosphere = _scenario(
            num_stokes, anderson_depth=anderson_depth
        )

    def peakmem_full_jacobian(self, num_stokes, anderson_depth):  # noqa: ARG002
        return self.engine.linearize(self.atmosphere).jacobian


def _production_atmosphere(
    config: sk.Config,
    geometry: sk.Geometry1D,
    num_wavelengths: int,
    *,
    calculate_derivatives: bool,
):
    altitude_grid_m = geometry.altitudes()
    atmosphere = sk.Atmosphere(
        geometry,
        config,
        calculate_derivatives=calculate_derivatives,
        numwavel=num_wavelengths,
    )
    spectral = np.linspace(0.85, 1.15, num_wavelengths)[np.newaxis, :]
    atmosphere.storage.total_extinction[:] = (
        2.5e-5 * np.exp(-altitude_grid_m[:, np.newaxis] / 7_500.0) + 1.0e-9
    ) * spectral
    atmosphere.storage.ssa[:] = 0.93
    atmosphere.leg_coeff.a1[0, :, :] = 1.0
    atmosphere.leg_coeff.a1[1, :, :] = 0.08
    atmosphere.leg_coeff.a1[2, :, :] = 0.5
    atmosphere.surface.albedo[:] = 0.2
    return atmosphere


def _production_scene():
    """Representative retrieval-loop geometry and configuration."""
    config = sk.Config()
    config.num_threads = 2
    config.num_stokes = 1
    config.num_streams = 16
    config.num_singlescatter_moments = 16
    config.single_scatter_source = sk.SingleScatterSource.Exact
    config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrdersRust
    config.num_successive_orders_incoming = 110
    config.num_successive_orders_outgoing = 110
    config.num_successive_orders_source_profiles = 1
    config.successive_orders_altitude_grid_m = np.array(
        [
            500.0,
            5_000.0,
            10_000.0,
            12_500.0,
            15_000.0,
            17_500.0,
            20_000.0,
            22_500.0,
            25_000.0,
            27_500.0,
            30_000.0,
            35_000.0,
            45_000.0,
            55_000.0,
        ]
    )
    config.successive_orders_max_iterations = 50
    config.successive_orders_relative_tolerance = 1.0e-6
    config.successive_orders_absolute_tolerance = 1.0e-12
    config.successive_orders_anderson_depth = 3
    config.los_refraction = True
    config.threading_model = sk.ThreadingModel.Wavelength
    config.threading_lib = sk.ThreadingLib.OpenMP

    altitude_grid_m = np.linspace(0.0, 70_000.0, 80)
    geometry = sk.Geometry1D(
        cos_sza=0.6,
        solar_azimuth=0.2,
        earth_radius_m=6_372_000.0,
        altitude_grid_m=altitude_grid_m,
        interpolation_method=sk.InterpolationMethod.LinearInterpolation,
        geometry_type=sk.GeometryType.Spherical,
    )
    geometry.refractive_index = 1.0 + 3.0e-4 * np.exp(-altitude_grid_m / 7_000.0)

    viewing = sk.ViewingGeometry()
    for index, tangent_altitude_m in enumerate(np.linspace(500.0, 55_000.0, 80)):
        viewing.add_ray(
            sk.TangentAltitudeSolar(
                tangent_altitude_m,
                -0.4 + 0.8 * index / 79,
                100_000.0,
                0.6,
            )
        )

    return config, geometry, viewing


def _production_inputs(num_wavelengths: int, *, calculate_derivatives: bool):
    """Representative retrieval-loop workload from the source design audit."""
    config, geometry, viewing = _production_scene()
    atmosphere = _production_atmosphere(
        config,
        geometry,
        num_wavelengths,
        calculate_derivatives=calculate_derivatives,
    )
    return config, geometry, viewing, atmosphere


class SuccessiveOrdersProductionSetup:
    """Geometry/setup cost for the production-shaped scalar workload."""

    params: ClassVar = [4, 10]
    param_names: ClassVar = ["num_wavelengths"]
    number = 1
    repeat = 3
    timeout = 600

    def setup(self, num_wavelengths):  # noqa: ARG002
        config, geometry, viewing = _production_scene()
        self.config = config
        self.geometry = geometry
        self.viewing = viewing

    def time_engine_setup(self, num_wavelengths):  # noqa: ARG002
        sk.Engine(self.config, self.geometry, self.viewing)

    def peakmem_engine_setup(self, num_wavelengths):  # noqa: ARG002
        return sk.Engine(self.config, self.geometry, self.viewing)


class SuccessiveOrdersProductionScalar:
    """80 LOS, 80 atmosphere levels, 110x110 angles, and 4-10 wavelengths."""

    params: ClassVar = [4, 10]
    param_names: ClassVar = ["num_wavelengths"]
    number = 1
    repeat = 3
    timeout = 600

    def setup(self, num_wavelengths):
        config, geometry, viewing, atmosphere = _production_inputs(
            num_wavelengths, calculate_derivatives=True
        )
        self.atmosphere = atmosphere
        self.primal_atmosphere = _production_atmosphere(
            config,
            geometry,
            num_wavelengths,
            calculate_derivatives=False,
        )
        self.engine = sk.Engine(config, geometry, viewing)
        self.primal_extinction = self.primal_atmosphere.storage.total_extinction.copy()
        self.linearization = self.engine.linearize(atmosphere)
        self.tangent = self.linearization.tangent_template[["extinction", "ssa"]]
        self.tangent["extinction"].data[:] = atmosphere.storage.total_extinction[:, 0]
        self.tangent["ssa"].data[:] = 0.1
        self.cotangent = xr.ones_like(self.linearization.value)
        self._state_scale = 1.0

    def time_primal(self, num_wavelengths):  # noqa: ARG002
        self.engine.calculate_radiance(self.primal_atmosphere)

    def time_primal_changed_state(self, num_wavelengths):  # noqa: ARG002
        self._state_scale = 1.0001 if self._state_scale == 1.0 else 1.0
        self.primal_atmosphere.storage.total_extinction[:] = (
            self.primal_extinction * self._state_scale
        )
        self.primal_atmosphere.mark_changed()
        self.engine.calculate_radiance(self.primal_atmosphere)

    def time_jvp(self, num_wavelengths):  # noqa: ARG002
        self.linearization.jvp(self.tangent)

    def time_vjp(self, num_wavelengths):  # noqa: ARG002
        self.linearization.vjp(
            self.cotangent,
            parameters=("extinction", "ssa"),
        )


class SuccessiveOrdersProductionPrimalPeakMemory:
    """Primal memory without derivative state retained by the fixture."""

    params: ClassVar = [4, 10]
    param_names: ClassVar = ["num_wavelengths"]
    timeout = 600

    def setup(self, num_wavelengths):
        config, geometry, viewing, self.atmosphere = _production_inputs(
            num_wavelengths, calculate_derivatives=False
        )
        self.engine = sk.Engine(config, geometry, viewing)

    def peakmem_primal(self, num_wavelengths):  # noqa: ARG002
        return self.engine.calculate_radiance(self.atmosphere)


class SuccessiveOrdersProductionVjpPeakMemory:
    """VJP memory without a second primal atmosphere retained by the fixture."""

    params: ClassVar = [4, 10]
    param_names: ClassVar = ["num_wavelengths"]
    timeout = 600

    def setup(self, num_wavelengths):
        config, geometry, viewing, atmosphere = _production_inputs(
            num_wavelengths, calculate_derivatives=True
        )
        self.engine = sk.Engine(config, geometry, viewing)
        self.linearization = self.engine.linearize(atmosphere)
        self.cotangent = xr.ones_like(self.linearization.value)

    def peakmem_vjp(self, num_wavelengths):  # noqa: ARG002
        return self.linearization.vjp(
            self.cotangent,
            parameters=("extinction", "ssa"),
        )
