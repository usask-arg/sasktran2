"""Generate SASKTRAN2's MT_CKD 4.3 server-side data artifacts.

This script contains no AER coefficient data or wrapper source. It requires a
separately obtained MT_CKD checkout and an externally built probe library,
both governed by their respective license terms.
"""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import math
import os
import struct
import subprocess
from collections.abc import Iterable, Iterator
from itertools import zip_longest
from pathlib import Path

MT_CKD_COMMIT = "1dad6e29363a2dc75d0eb642b7321b5db51079cc"
LBLRTM_COMMIT = "f66bb596e1686e5ff162a6fadc014e91d3ce01d4"
COEFFICIENT_SHA256 = "06b5600ecb0f5a3417c46d226555bc516f5a55ae93b19b6e7b5431e50f8f0459"
REFERENCE_SHA256 = "0fdc39229ae020f2979511acfc81600ddfaa0b9bc0bbb03bb86a5ce3d01a9760"

SPECTRAL_POINTS = 1991
FIELD_COUNT = 22
WAVENUMBERS = tuple(10.0 * index for index in range(SPECTRAL_POINTS))
RADCN2 = 1.438_775_2
ALOSMT = 2.686_777_5e19
XLOSMT = 2.686_75e19


def zip_equal(*iterables: Iterable[float]) -> Iterator[tuple[float, ...]]:
    """Python 3.9-compatible equivalent of zip(..., strict=True)."""
    sentinel = object()
    for values in zip_longest(*iterables, fillvalue=sentinel):
        if sentinel in values:
            message = "iterables have different lengths"
            raise ValueError(message)
        yield values


def command_output(command: list[str], cwd: Path | None = None) -> str:
    return subprocess.run(
        command,
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def verify_checkout(checkout: Path) -> None:
    expected_paths = (
        checkout / "build" / "make_cntnm",
        checkout / "data" / "absco-ref_wv-mt-ckd.nc",
        checkout / "LBLRTM" / "src" / "contnm.f90",
    )
    missing = [str(path) for path in expected_paths if not path.exists()]
    if missing:
        message = f"MT_CKD checkout is incomplete; missing: {', '.join(missing)}"
        raise RuntimeError(message)

    revisions = (
        (checkout, MT_CKD_COMMIT, "MT_CKD"),
        (checkout / "LBLRTM", LBLRTM_COMMIT, "LBLRTM submodule"),
    )
    for repository, expected, label in revisions:
        actual = command_output(["git", "rev-parse", "HEAD"], repository)
        if actual != expected:
            message = f"expected {label} revision {expected}, found {actual}"
            raise RuntimeError(message)

    checkout_source_diff = subprocess.run(
        ["git", "diff", "--quiet", "--", "src", "data"],
        cwd=checkout,
        check=False,
    )
    if checkout_source_diff.returncode != 0:
        message = "the pinned MT_CKD source or data directory has local modifications"
        raise RuntimeError(message)
    source_diff = subprocess.run(
        ["git", "diff", "--quiet", "--", "src"],
        cwd=checkout / "LBLRTM",
        check=False,
    )
    if source_diff.returncode != 0:
        message = "the pinned LBLRTM source directory has local modifications"
        raise RuntimeError(message)


class Probe:
    def __init__(self, library: Path, data_directory: Path) -> None:
        library = library.resolve()
        os.chdir(data_directory)
        self._library = ctypes.CDLL(str(library))
        self._library.set_inputs.argtypes = [ctypes.c_double] * 6
        self._library.set_flags.argtypes = [ctypes.c_double] * 7
        self._library.run_mtckd.argtypes = []
        self._library.get_absrb.argtypes = [ctypes.c_int]
        self._library.get_absrb.restype = ctypes.c_double
        self._library.get_n2_data.argtypes = [ctypes.c_int, ctypes.c_int]
        self._library.get_n2_data.restype = ctypes.c_double

    def set_inputs(self, state: tuple[float, ...]) -> None:
        self._library.set_inputs(*state)

    def set_flags(self, flags: list[float]) -> None:
        self._library.set_flags(*flags)

    def run(self) -> None:
        self._library.run_mtckd()

    def absorption(self, index: int) -> float:
        return self._library.get_absrb(index)

    def n2_data(self, field: int, index: int) -> float:
        return self._library.get_n2_data(field, index)


def radiation(wavenumber: float, temperature: float) -> float:
    ratio = wavenumber / (temperature / RADCN2)
    if ratio <= 0.01:
        return 0.5 * ratio * wavenumber
    if ratio <= 10.0:
        exponential = math.exp(-ratio)
        return wavenumber * (1.0 - exponential) / (1.0 + exponential)
    return wavenumber


def atmospheric_state(
    pressure_mb: float,
    temperature: float,
    h2o_vmr: float,
    co2_vmr: float,
    o3_vmr: float,
    path_length_cm: float,
) -> dict[str, float]:
    initial_total = (
        ALOSMT * (pressure_mb / 1013.0) * (273.0 / temperature) * path_length_cm
    )
    dry = initial_total * (1.0 - h2o_vmr)
    h2o = initial_total if abs(h2o_vmr - 1.0) < 1.0e-5 else h2o_vmr * dry
    co2 = co2_vmr * dry
    o3 = o3_vmr * dry
    o2 = 0.21 * dry
    broadening = (0.78 + 0.009) * dry
    total = broadening + h2o + co2 + o3 + o2
    x_h2o = h2o / total
    x_o2 = o2 / total
    x_n2 = 1.0 - x_h2o - x_o2
    return {
        "h2o": h2o,
        "co2": co2,
        "o3": o3,
        "o2": o2,
        "n2": x_n2 * total,
        "total": total,
        "xh": x_h2o,
        "xo2": x_o2,
        "xn2": x_n2,
        "rho": (pressure_mb / 1013.0) * (296.0 / temperature),
        "amagat": (pressure_mb / 1013.0) * (273.0 / temperature),
    }


def isolated_component(
    probe: Probe,
    flag_index: int,
    pressure_mb: float = 1013.0,
    temperature: float = 296.0,
    h2o_vmr: float = 0.01,
    co2_vmr: float = 400.0e-6,
    o3_vmr: float = 2.0e-6,
    path_length_cm: float = 100.0,
) -> tuple[list[float], dict[str, float]]:
    state = (
        pressure_mb,
        temperature,
        h2o_vmr,
        co2_vmr,
        o3_vmr,
        path_length_cm,
    )
    flags = [0.0] * 7
    flags[flag_index] = 1.0
    probe.set_inputs(state)
    probe.set_flags(flags)
    probe.run()
    values = []
    for index, wavenumber in enumerate(WAVENUMBERS, 1):
        radiation_value = radiation(wavenumber, temperature)
        value = probe.absorption(index)
        values.append(0.0 if radiation_value == 0.0 else value / radiation_value)
    return values, atmospheric_state(*state)


def divide(values: list[float], factor: float) -> list[float]:
    return [value / factor for value in values]


def exponent(base: list[float], other: list[float], denominator: float) -> list[float]:
    return [
        math.log(second / first) / denominator if first > 0.0 and second > 0.0 else 0.0
        for first, second in zip_equal(base, other)
    ]


def solve_two(
    first_values: list[float],
    first_a: float,
    first_b: float,
    second_values: list[float],
    second_a: float,
    second_b: float,
) -> tuple[list[float], list[float]]:
    determinant = first_a * second_b - second_a * first_b
    a = [
        (first * second_b - second * first_b) / determinant
        for first, second in zip_equal(first_values, second_values)
    ]
    b = [
        (first_a * second - second_a * first) / determinant
        for first, second in zip_equal(first_values, second_values)
    ]
    return a, b


def coefficient_fields(probe: Probe) -> list[list[float]]:
    fields: list[list[float]] = []

    self_296, state = isolated_component(probe, 0, temperature=296.0)
    self_base = divide(self_296, state["h2o"] * state["xh"] * state["rho"])
    self_250, state = isolated_component(probe, 0, temperature=250.0)
    self_250 = divide(self_250, state["h2o"] * state["xh"] * state["rho"])
    self_exponent = exponent(self_base, self_250, math.log(296.0 / 250.0))
    foreign, state = isolated_component(probe, 1)
    foreign_base = divide(
        foreign,
        state["h2o"] * (1.0 - state["xh"]) * state["rho"],
    )
    fields.extend((self_base, self_exponent, foreign_base))

    co2_246, state = isolated_component(probe, 2, temperature=246.0)
    co2_base = divide(co2_246, state["co2"] * state["rho"] * 1.0e-20)
    co2_296, state = isolated_component(probe, 2, temperature=296.0)
    co2_296 = divide(co2_296, state["co2"] * state["rho"] * 1.0e-20)
    co2_exponent = exponent(co2_base, co2_296, math.log(296.0 / 246.0))
    fields.extend((co2_base, co2_exponent))

    o3_values = []
    for temperature in (273.15, 293.15, 253.15):
        values, state = isolated_component(probe, 3, temperature=temperature)
        o3_values.append(divide(values, state["o3"] * 1.0e-20))
    o3_constant, o3_plus, o3_minus = o3_values
    o3_linear = [(plus - minus) / 40.0 for plus, minus in zip_equal(o3_plus, o3_minus)]
    o3_quadratic = [
        (plus + minus - 2.0 * constant) / 800.0
        for plus, minus, constant in zip_equal(
            o3_plus,
            o3_minus,
            o3_constant,
        )
    ]
    fields.extend((o3_constant, o3_linear, o3_quadratic))

    o2_296, state = isolated_component(probe, 4, temperature=296.0)
    o2_fundamental = divide(o2_296, state["o2"] * 1.0e-20 * state["amagat"])
    o2_250, state = isolated_component(probe, 4, temperature=250.0)
    o2_250 = divide(o2_250, state["o2"] * 1.0e-20 * state["amagat"])
    o2_energy = exponent(
        o2_fundamental,
        o2_250,
        (1.0 / 296.0) - (1.0 / 250.0),
    )
    o2_fundamental = [
        value if wavenumber < 3000.0 else 0.0
        for wavenumber, value in zip_equal(WAVENUMBERS, o2_fundamental)
    ]
    o2_energy = [
        value if wavenumber < 3000.0 else 0.0
        for wavenumber, value in zip_equal(WAVENUMBERS, o2_energy)
    ]
    fields.extend((o2_fundamental, o2_energy))

    o2_all, state = isolated_component(probe, 4, temperature=296.0)
    infrared1_factor = (
        (state["o2"] / XLOSMT)
        * state["amagat"]
        * (state["xo2"] / 0.446 + 0.3 * state["xn2"] / 0.446 + state["xh"])
    )
    infrared1 = divide(o2_all, infrared1_factor)
    infrared1 = [
        value if 7000.0 <= wavenumber < 9000.0 else 0.0
        for wavenumber, value in zip_equal(WAVENUMBERS, infrared1)
    ]
    infrared2_factor = (
        (state["o2"] / state["total"])
        * (1.0 / 0.209)
        * (state["o2"] * 1.0e-20)
        * state["rho"]
    )
    infrared2 = divide(o2_all, infrared2_factor)
    infrared2 = [
        value if 9000.0 <= wavenumber < 12000.0 else 0.0
        for wavenumber, value in zip_equal(WAVENUMBERS, infrared2)
    ]
    infrared3_factor = (state["o2"] / XLOSMT) * state["amagat"]
    infrared3 = divide(o2_all, infrared3_factor)
    infrared3 = [
        value if 12000.0 <= wavenumber < 14000.0 else 0.0
        for wavenumber, value in zip_equal(WAVENUMBERS, infrared3)
    ]
    visible_factor = (
        (state["o2"] / state["total"]) * (state["o2"] * 1.0e-20) * state["amagat"]
    )
    visible = divide(o2_all, visible_factor)
    visible = [
        value if wavenumber >= 14000.0 else 0.0
        for wavenumber, value in zip_equal(WAVENUMBERS, visible)
    ]
    fields.extend((infrared1, infrared2, infrared3, visible))

    for temperature in (296.0, 220.0):
        equations = []
        for h2o_vmr in (0.01, 0.2):
            values, state = isolated_component(
                probe,
                5,
                temperature=temperature,
                h2o_vmr=h2o_vmr,
            )
            tau = (state["n2"] / XLOSMT) * state["amagat"]
            equations.append((divide(values, tau), 1.0 - state["xo2"], state["xo2"]))
        rotation, rotation_o2 = solve_two(*equations[0], *equations[1])
        rotation_scale = [
            1.0 + (o2_value / value) * 0.21 / 0.79 if value != 0.0 else 1.0
            for value, o2_value in zip_equal(rotation, rotation_o2)
        ]
        rotation = [
            value if wavenumber < 1000.0 else 0.0
            for wavenumber, value in zip_equal(WAVENUMBERS, rotation)
        ]
        rotation_scale = [
            value if wavenumber < 1000.0 else 1.0
            for wavenumber, value in zip_equal(WAVENUMBERS, rotation_scale)
        ]
        fields.extend((rotation, rotation_scale))

    for field in range(3):
        raw = [probe.n2_data(field, index) for index in range(1, 229)]
        fields.append(raw + [0.0] * (SPECTRAL_POINTS - len(raw)))

    n2_all, state = isolated_component(probe, 5, temperature=296.0)
    overtone_factor = (
        (state["n2"] / XLOSMT)
        * state["amagat"]
        * (state["xn2"] + state["xo2"] + state["xh"])
    )
    overtone = divide(n2_all, overtone_factor)
    overtone = [
        value if 4000.0 <= wavenumber < 6000.0 else 0.0
        for wavenumber, value in zip_equal(WAVENUMBERS, overtone)
    ]
    fields.append(overtone)

    if len(fields) != FIELD_COUNT:
        message = f"expected {FIELD_COUNT} fields, generated {len(fields)}"
        raise RuntimeError(message)
    return fields


def write_coefficients(path: Path, fields: list[list[float]]) -> None:
    payload = bytearray(b"MTCKD43\0")
    payload.extend(struct.pack("<II", SPECTRAL_POINTS, len(fields)))
    for field in fields:
        if len(field) != SPECTRAL_POINTS:
            message = "coefficient field has the wrong spectral length"
            raise RuntimeError(message)
        payload.extend(struct.pack(f"<{SPECTRAL_POINTS}d", *field))
    path.write_bytes(payload)


def write_reference(path: Path, probe: Probe) -> None:
    scenarios = (
        (1013.25, 288.15, 0.01, 400.0e-6, 2.0e-6, 100.0),
        (500.0, 250.0, 0.002, 420.0e-6, 1.0e-6, 100.0),
        (80.0, 220.0, 5.0e-5, 380.0e-6, 5.0e-6, 37.5),
        (900.0, 310.0, 0.2, 1.0e-3, 1.0e-5, 250.0),
    )
    payload = bytearray()
    for scenario in scenarios:
        probe.set_inputs(scenario)
        probe.set_flags([1.0] * 7)
        probe.run()
        values = [probe.absorption(index) for index in range(1, SPECTRAL_POINTS + 1)]
        values[0] = 0.0
        payload.extend(struct.pack(f"<{SPECTRAL_POINTS}d", *values))
    path.write_bytes(payload)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require_hash(path: Path, expected: str) -> None:
    actual = sha256(path)
    print(f"{actual}  {path}")  # noqa: T201
    if actual != expected:
        message = f"{path} does not reproduce the expected SHA-256 {expected}"
        raise RuntimeError(message)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mt-ckd-checkout",
        type=Path,
        required=True,
        help="separately obtained AER-RC/MT_CKD 4.3 checkout",
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        required=True,
        help="directory in which to write the two server artifacts",
    )
    parser.add_argument(
        "--probe-library",
        type=Path,
        required=True,
        help="externally built, ABI-compatible MT_CKD probe library",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    checkout = args.mt_ckd_checkout.resolve()
    output_directory = args.output_directory.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    verify_checkout(checkout)

    probe = Probe(args.probe_library, checkout / "data")
    coefficient_path = output_directory / "mt_ckd_4_3.bin"
    reference_path = output_directory / "mt_ckd_4_3_reference.bin"
    write_coefficients(coefficient_path, coefficient_fields(probe))
    write_reference(reference_path, probe)
    require_hash(coefficient_path, COEFFICIENT_SHA256)
    require_hash(reference_path, REFERENCE_SHA256)


if __name__ == "__main__":
    main()
