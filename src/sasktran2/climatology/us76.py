from __future__ import annotations

import numpy as np

import sasktran2 as sk
from sasktran2.units import celsius_to_kelvin

_ALTS = np.array(
    [
        -1000,
        0,
        1000,
        2000,
        3000,
        4000,
        5000,
        6000,
        7000,
        8000,
        9000,
        10000,
        15000,
        20000,
        25000,
        30000,
        40000,
        50000,
        60000,
        70000,
        80000,
    ]
)

_TEMPERATURE_C = np.array(
    [
        21.50,
        15.00,
        8.50,
        2.00,
        -4.49,
        -10.98,
        -17.47,
        -23.96,
        -30.45,
        -36.94,
        -43.42,
        -49.90,
        -56.50,
        -56.50,
        -51.60,
        -46.64,
        -22.80,
        -2.5,
        -26.13,
        -53.57,
        -74.51,
    ]
)

_PRESSURE = np.array(
    [
        11.39,
        10.13,
        8.988,
        7.950,
        7.012,
        6.166,
        5.405,
        4.722,
        4.111,
        3.565,
        3.080,
        2.650,
        1.211,
        0.5529,
        0.2549,
        0.1197,
        0.0287,
        0.007978,
        0.002196,
        0.00052,
        0.00011,
    ]
)


def _log_pressure_with_top_extrapolation(altitudes_m: np.ndarray) -> np.ndarray:
    log_pressure = np.log(_PRESSURE * 1e4)
    log_pressure_interp = np.interp(
        altitudes_m,
        _ALTS,
        log_pressure,
        left=log_pressure[0],
    )

    top_mask = altitudes_m > _ALTS[-1]
    if np.any(top_mask):
        top_slope = (log_pressure[-1] - log_pressure[-2]) / (_ALTS[-1] - _ALTS[-2])
        log_pressure_interp[top_mask] = log_pressure[-1] + top_slope * (
            altitudes_m[top_mask] - _ALTS[-1]
        )

    return log_pressure_interp


def add_us76_standard_atmosphere(atmo: sk.Atmosphere):
    geo = atmo.model_geometry

    atmo.pressure_pa = np.exp(_log_pressure_with_top_extrapolation(geo.altitudes()))
    atmo.temperature_k = np.interp(
        geo.altitudes(),
        _ALTS,
        celsius_to_kelvin(_TEMPERATURE_C),
        left=celsius_to_kelvin(_TEMPERATURE_C[0]),
        right=celsius_to_kelvin(_TEMPERATURE_C[-1]),
    )
