"""Compare line absorption and emission radiance with temperature derivatives.

Run with a release build on each revision:
    pixi run python tools/benchmarks/line_temperature_derivatives.py --output /tmp/lines.json

Uses a local synthetic O2 line list and a power-law partition function, with no
downloads or line mixing. Timings include allocations and atmospheric assembly.
Temperature-off radiance retains VMR derivatives; all-off disables all Jacobians.
The combined optical case falls back to separate calls on older revisions.
Output includes timing samples and an NPZ of values/Jacobians for comparison.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import numpy as np
import sasktran2 as sk
from sasktran2._core_rust import LineDatabaseType, PyLineAbsorber
from sasktran2.constituent.base import Constituent
from sasktran2.optical.hitran import LineAbsorber


class FixedEmission(Constituent):
    def add_to_atmosphere(self, atmo):
        profile = np.exp(-(((atmo.model_geometry.altitudes() - 60_000) / 20_000) ** 2))
        atmo.storage.emission_source[:] += profile[:, None]

    def register_derivative(self, atmo, name):
        pass


def measure(cases, blocks, repeats):
    for operation in cases.values():
        for _ in range(3):
            operation()
    samples = {name: [] for name in cases}
    names = list(cases)
    rng = np.random.default_rng(8729)
    for _ in range(blocks):
        for index in rng.permutation(len(names)):
            name = names[index]
            start = perf_counter()
            for _ in range(repeats):
                cases[name]()
            samples[name].append((perf_counter() - start) * 1000 / repeats)
    return {
        name: {"median_ms": float(np.median(values)), "samples_ms": values}
        for name, values in samples.items()
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--lines", type=int, default=100)
    parser.add_argument("--levels", type=int, default=40)
    parser.add_argument("--spectral-points", type=int, default=4001)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--blocks", type=int, default=11)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if (
        min(args.lines, args.spectral_points, args.threads, args.blocks, args.repeats)
        < 1
    ):
        parser.error(
            "line, spectral point, thread, block and repeat counts must be positive"
        )
    if args.levels < 2:
        parser.error("at least two altitude levels are required")

    with TemporaryDirectory() as directory:
        records = [
            f"{7:2d}{1:1d}{center:12.6f}{1e-24:10.3e}{0.1:10.3e}"
            f"{0.06:5.3f}{0.10:5.3f}{i * 20:10.4f}{0.7:4.2f}{0.003:8.6f}"
            + " " * 79
            + f"{1.0:7.1f}{1.0:7.1f}"
            for i, center in enumerate(np.linspace(13100, 13108, args.lines))
        ]
        (Path(directory) / "O2.data").write_text("\n".join(records) + "\n")
        absorber = LineAbsorber.__new__(LineAbsorber)
        absorber._internal = PyLineAbsorber(
            LineDatabaseType.HITRAN,
            "O2",
            directory,
            py_tips=lambda _mol, _iso, t: t**1.5,
            py_molmass=lambda _mol, _iso: 31.9988,
        )

    config = sk.Config()
    config.num_threads = args.threads
    config.single_scatter_source = sk.SingleScatterSource.NoSource
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    config.emission_source = sk.EmissionSource.VolumeEmissionRate
    altitudes = np.linspace(0, 120_000, args.levels)
    wavenumbers = np.linspace(13100, 13108, args.spectral_points)
    vmr = np.full(args.levels, 0.21)
    geometry = sk.Geometry1D(
        0.6,
        0,
        6_372_000,
        altitudes,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.Spherical,
    )
    viewing = sk.ViewingGeometry()
    for tangent in [10_000, 40_000, 70_000]:
        viewing.add_ray(sk.TangentAltitudeSolar(tangent, 0, 200_000, 0.6))
    engine = sk.Engine(config, geometry, viewing)

    def atmosphere(derivatives, temperature_derivative):
        atmo = sk.Atmosphere(
            geometry,
            config,
            wavenumber_cminv=wavenumbers,
            calculate_derivatives=derivatives,
            temperature_derivative=temperature_derivative,
            pressure_derivative=False,
            specific_humidity_derivative=False,
            legendre_derivative=False,
        )
        atmo.temperature_k = np.linspace(290, 190, args.levels)
        atmo.pressure_pa = np.geomspace(101325, 0.1, args.levels)
        atmo["o2"] = sk.constituent.VMRAltitudeAbsorber(absorber, altitudes, vmr)
        atmo["emission"] = FixedEmission()
        return atmo

    all_off = atmosphere(False, False)
    temperature_off = atmosphere(True, False)
    temperature_on = atmosphere(True, True)

    def forward():
        return absorber.atmosphere_quantities(temperature_on, vmr=vmr)

    def derivative():
        return absorber.optical_derivatives(temperature_on, vmr=vmr)

    def combined():
        if hasattr(absorber._internal, "atmosphere_quantities_and_derivatives"):
            return absorber.atmosphere_quantities_and_derivatives(
                temperature_on, vmr=vmr
            )
        return forward(), derivative()

    quantities, derivatives = combined()
    values = {
        "cross_section": forward().cross_section,
        "combined_cross_section": quantities.cross_section,
        "d_cross_section_dT": derivatives["temperature_k"].cross_section,
    }
    np.testing.assert_allclose(
        values["combined_cross_section"],
        values["cross_section"],
        rtol=1e-10,
        atol=1e-35,
    )
    reference = engine.calculate_radiance(all_off)
    result = engine.calculate_radiance(temperature_on)
    np.testing.assert_array_equal(reference.radiance, result.radiance)
    values.update({name: value.values for name, value in result.data_vars.items()})

    cases = {
        "optical_forward": forward,
        "optical_derivative": derivative,
        "optical_separate": lambda: (forward(), derivative()),
        "optical_combined": combined,
        "radiance_all_off": lambda: engine.calculate_radiance(all_off),
        "radiance_temperature_off": lambda: engine.calculate_radiance(temperature_off),
        "radiance_temperature_on": lambda: engine.calculate_radiance(temperature_on),
    }
    metadata = vars(args).copy()
    metadata.pop("output")
    metadata.update({"rays": 3, "partition_function": "T**1.5", "line_mixing": False})
    timings = measure(cases, args.blocks, args.repeats)
    for name, timing in timings.items():
        print(f"{name:30s} {timing['median_ms']:8.3f} ms", flush=True)
    if args.output:
        args.output.write_text(
            json.dumps({"metadata": metadata, "timings": timings}, indent=2)
        )
        np.savez(args.output.with_suffix(".npz"), **values)


if __name__ == "__main__":
    main()
