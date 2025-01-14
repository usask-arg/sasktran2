from __future__ import annotations

import numpy as np


def celsius_to_kelvin(temperature_c: np.array):
    return temperature_c + 273.15


def kelvin_to_celsius(temperature_k: np.array):
    return temperature_k - 273.15


def wavlength_nm_to_wavenumber_cminv(wavelengths_nm: np.array):
    return 1e7 / wavelengths_nm


def wavenumber_cminv_to_wavlength_nm(wavenumber_cminv: np.array):
    return 1e7 / wavenumber_cminv
