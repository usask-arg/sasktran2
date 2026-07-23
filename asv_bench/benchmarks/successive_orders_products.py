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
