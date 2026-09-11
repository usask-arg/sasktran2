"""Benchmark Greek transforms, Mie setup, and serial/threaded integration.

Run on each revision with the same release build and machine:
    pixi run python tools/benchmarks/mie_greek.py --output /tmp/mie-benchmark.json

Select a larger transform without running particle-size integration:
    pixi run python tools/benchmarks/mie_greek.py --orders 10000 --skip-integration

The optional output saves timings in JSON and numerical results alongside it in
an NPZ file, so performance comparisons can also check numerical equivalence.
Integration timings exclude the separately measured angular-basis setup.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from sasktran2._core_rust import PyMieIntegrator
from sasktran2.legendre import compute_greek_coefficients
from scipy.special import roots_legendre


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--orders", type=int, nargs="+", default=[64, 256, 512])
    parser.add_argument("--skip-integration", action="store_true")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if any(order < 1 for order in args.orders):
        parser.error("--orders must be positive")
    timings = {}
    values = {}

    def measure(name, calculate):
        result = calculate()  # Warm up before timing.
        elapsed = []
        for _ in range(args.repeats):
            # In particular, do not retain a previous integrator's large basis
            # while constructing the next one, or time its destruction.
            del result
            start = perf_counter()
            result = calculate()
            elapsed.append(perf_counter() - start)
        timings[name] = min(elapsed)
        print(f"{name:28s} {1000 * timings[name]:10.3f} ms", flush=True)
        return result

    angles = np.linspace(0, 180, 1801)
    x = np.cos(np.deg2rad(angles))
    phase = [
        (1 + i / 10 + x**2 + np.exp(10 * (x - 1)))[None, :]
        * np.linspace(1, 2, 8)[:, None]
        for i in range(6)
    ]
    for order in args.orders:
        result = measure(
            f"greek_{order}",
            lambda: compute_greek_coefficients(*phase, angles, order),
        )
        values[f"greek_{order}"] = np.stack(result)
        cos_angles, _ = roots_legendre(2 * order)
        measure(f"init_{order}", lambda: PyMieIntegrator(cos_angles, order, 1))

    for distributions in [] if args.skip_integration else [1, 8]:
        order = 128
        cos_angles, weights = roots_legendre(2 * order)
        size = np.linspace(0.01, 100, 512)
        pdf = np.array([np.exp(-size / (5 + i)) for i in range(distributions)])
        for threads in [1, 4]:
            integrator = PyMieIntegrator(cos_angles, order, threads)

            def integrate():
                outputs = (
                    [np.zeros(distributions) for _ in range(2)]
                    + [np.zeros((distributions, len(cos_angles))) for _ in range(4)]
                    + [np.zeros((distributions, order)) for _ in range(6)]
                )
                integrator.integrate(
                    550.0,
                    1.5 - 0.01j,
                    size,
                    pdf,
                    np.ones_like(size),
                    weights,
                    *outputs,
                )
                return outputs

            result = measure(f"mie_{distributions}dist_{threads}threads", integrate)
            for i, value in enumerate(result):
                values[f"mie_{distributions}_{threads}_{i}"] = value

    if args.output:
        args.output.write_text(json.dumps(timings, indent=2) + "\n")
        np.savez(args.output.with_suffix(".npz"), **values)


if __name__ == "__main__":
    main()
