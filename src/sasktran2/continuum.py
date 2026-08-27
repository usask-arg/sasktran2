from __future__ import annotations

import hashlib
from functools import cache
from pathlib import Path

import numpy as np

import sasktran2 as sk
from sasktran2._core_rust import _mt_ckd, _mt_ckd_linearized
from sasktran2.constituent.base import Constituent
from sasktran2.util.interpolation import linear_interpolating_matrix

__all__ = [
    "MT_CKD_AER_LICENSE",
    "MT_CKD_WAVENUMBERS_CM_INV",
    "MTCKDContinuum",
    "mt_ckd",
    "mt_ckd_linearized",
]

_MT_CKD_DATA_KEY = "continuum/mt_ckd_4_3.bin"
_MT_CKD_LICENSE_KEY = "continuum/MT_CKD_LICENSE.txt"
_MT_CKD_DATA_SHA256 = "06b5600ecb0f5a3417c46d226555bc516f5a55ae93b19b6e7b5431e50f8f0459"
_MT_CKD_LICENSE_SHA256 = (
    "3648e3e8a231da514f923409b9b5e41be5e1746e92692af3ce3ac89bc6803b69"
)

MT_CKD_WAVENUMBERS_CM_INV = np.arange(0.0, 19_901.0, 10.0)
"""The native 1,991-point MT_CKD wavenumber grid in cm^-1."""
MT_CKD_WAVENUMBERS_CM_INV.setflags(write=False)

MT_CKD_AER_LICENSE = """Copyright ©, Atmospheric and Environmental Research, Inc., 2022

All rights reserved. This source code was developed as part of the LBLRTM
software and is designed for scientific and research purposes. Atmospheric
and Environmental Research Inc. (AER) grants USER the right to download,
install, use and copy this software and data for scientific and research
purposes only. This software and data may be redistributed as long as this
copyright notice is reproduced on any copy made and appropriate acknowledgment
is given to AER. This software and data or any modified version of this
software and data may not be incorporated into proprietary software or data or
commercial software or data offered for sale without the express written
consent of AER. This software and data are provided as is without any expressed
or implied warranties.

Address questions to: aer_contnm@aer.com

General reference: Mlawer et al. (2012),
https://doi.org/10.1098/rsta.2011.0295"""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _verified_standard_database_file(key: str, expected_sha256: str) -> Path:
    path = sk.database.StandardDatabase().path(key)
    if path is None or not path.exists():
        msg = f"Failed to download the required MT_CKD file {key!r}"
        raise OSError(msg)

    actual_sha256 = _sha256(path)
    if actual_sha256 != expected_sha256:
        path.unlink()
        msg = (
            f"The downloaded MT_CKD file {key!r} failed checksum validation; "
            f"expected {expected_sha256}, found {actual_sha256}. The invalid "
            "cached file has been removed."
        )
        raise OSError(msg)
    return path


@cache
def _mt_ckd_data_file() -> Path:
    local_data_file = sk.appconfig.database_root().joinpath(_MT_CKD_DATA_KEY)
    if not local_data_file.exists():
        print(  # noqa: T201
            "MT_CKD 4.3 coefficient data is licensed separately from "
            "SASKTRAN2 under the following AER terms:\n\n"
            f"{MT_CKD_AER_LICENSE}\n"
        )

    _verified_standard_database_file(
        _MT_CKD_LICENSE_KEY,
        _MT_CKD_LICENSE_SHA256,
    )
    return _verified_standard_database_file(
        _MT_CKD_DATA_KEY,
        _MT_CKD_DATA_SHA256,
    )


def mt_ckd(
    pressure_pa: np.ndarray,
    temperature_k: np.ndarray,
    vmr_h2o: np.ndarray,
    vmr_co2: np.ndarray,
    vmr_o3: np.ndarray,
    *,
    include_rayleigh: bool = True,
) -> np.ndarray:
    """Calculate the MT_CKD 4.3 extinction coefficient in m^-1.

    The returned columns correspond to
    :data:`MT_CKD_WAVENUMBERS_CM_INV`. The calculation uses a fixed one-metre
    reference column internally, so no geometric path length is an input.

    By default the result includes the historical Rayleigh-scattering term
    returned by the standalone MT_CKD wrapper, preserving low-level
    compatibility with :mod:`sasktran2_ext`. Set ``include_rayleigh=False``
    when only continuum absorption is required. :class:`MTCKDContinuum` always
    excludes Rayleigh so it can be combined safely with SASKTRAN2's Rayleigh
    constituent.

    The separately licensed AER coefficient data is downloaded from the
    SASKTRAN2 database on first use. Its scientific/research-only license is
    printed before that download, stored beside the cached data, and reproduced
    in :data:`MT_CKD_AER_LICENSE`. It is not covered by SASKTRAN2's MIT license.
    """
    return _mt_ckd(
        pressure_pa,
        temperature_k,
        vmr_h2o,
        vmr_co2,
        vmr_o3,
        str(_mt_ckd_data_file()),
        include_rayleigh,
    )


def mt_ckd_linearized(
    pressure_pa: np.ndarray,
    temperature_k: np.ndarray,
    vmr_h2o: np.ndarray,
    vmr_co2: np.ndarray,
    vmr_o3: np.ndarray,
    *,
    include_rayleigh: bool = True,
) -> tuple[np.ndarray, ...]:
    """Calculate MT_CKD extinction and its five analytic Jacobians.

    The return values are extinction coefficient in m^-1 followed by
    derivatives with respect to pressure, temperature, H2O VMR, CO2 VMR, and
    O3 VMR. Their columns correspond to
    :data:`MT_CKD_WAVENUMBERS_CM_INV`. On first use, the separately licensed
    AER coefficient data is downloaded as described in :func:`mt_ckd`; it is
    not covered by SASKTRAN2's MIT license.

    ``include_rayleigh`` has the same compatibility behavior as in
    :func:`mt_ckd`.
    """
    return _mt_ckd_linearized(
        pressure_pa,
        temperature_k,
        vmr_h2o,
        vmr_co2,
        vmr_o3,
        str(_mt_ckd_data_file()),
        include_rayleigh,
    )


class MTCKDContinuum(Constituent):
    """MT_CKD 4.3 continuum absorption with analytic linearization.

    The constructor and atmosphere requirements match
    ``sasktran2_ext.continuum.MTCKDContinuum``. The numerical-difference
    options are retained for source compatibility but are no longer used;
    pressure, temperature, H2O, CO2, and O3 derivatives are evaluated
    analytically by the Rust implementation.

    This constituent adds continuum absorption only. MT_CKD's historical
    Rayleigh term is excluded so scattering can be supplied by SASKTRAN2's
    dedicated Rayleigh constituent without double counting.

    MT_CKD's AER coefficient data is not covered by SASKTRAN2's MIT license.
    On first use it is downloaded from the SASKTRAN2 database after its
    scientific/research-only terms are printed, and the AER notice is stored
    beside the cached data. See :data:`MT_CKD_AER_LICENSE` for the full terms.

    Parameters
    ----------
    h2o_name : str, optional
        Name of the atmosphere's H2O constituent.
    co2_name : str, optional
        Name of the atmosphere's CO2 constituent.
    o3_name : str, optional
        Name of the atmosphere's O3 constituent.
    numeric_wf_fractional_change : float, optional
        Retained for compatibility with ``sasktran2_ext``.
    numeric_wf_central_difference : bool, optional
        Retained for compatibility with ``sasktran2_ext``.
    """

    def __init__(
        self,
        h2o_name: str = "H2O",
        co2_name: str = "CO2",
        o3_name: str = "O3",
        numeric_wf_fractional_change: float = 1e-2,
        numeric_wf_central_difference: bool = True,
    ) -> None:
        self._h2o_name = h2o_name
        self._co2_name = co2_name
        self._o3_name = o3_name
        self._mtckd_wavenumbers = MT_CKD_WAVENUMBERS_CM_INV
        self._fractional_change = numeric_wf_fractional_change
        self._central_difference = numeric_wf_central_difference
        self._linearized_cache: tuple[np.ndarray, ...] | None = None
        self._wavenumber_interpolator: np.ndarray | None = None
        self._vmr_interpolators: tuple[np.ndarray, ...] | None = None

    def _inputs(
        self, atmo: sk.Atmosphere
    ) -> tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...]]:
        if atmo.wavelengths_nm is None:
            msg = "It is required to give the Atmosphere object wavelengths to use the continuum constituent"
            raise ValueError(msg)
        if atmo.pressure_pa is None:
            msg = "It is required to set the pressure_pa property in the Atmosphere object to use the continuum constituent"
            raise ValueError(msg)
        if atmo.temperature_k is None:
            msg = "It is required to set the temperature_k property in the Atmosphere object to use the continuum constituent"
            raise ValueError(msg)

        constituents = []
        for species_name in (self._h2o_name, self._co2_name, self._o3_name):
            if atmo[species_name] is None:
                msg = f"It is required to add an {species_name} constituent to the Atmosphere object to use the continuum constituent"
                raise ValueError(msg)
            constituents.append(atmo[species_name])

        altitude_grid = atmo.model_geometry.altitudes()
        vmrs = []
        interpolators = []
        for constituent in constituents:
            interpolator = linear_interpolating_matrix(
                constituent.altitudes_m,
                altitude_grid,
                "zero",
            )
            interpolators.append(interpolator)
            vmrs.append(interpolator @ constituent.vmr)
        return tuple(vmrs), tuple(interpolators)

    def add_to_atmosphere(self, atmo: sk.Atmosphere) -> None:
        (h2o_vmr, co2_vmr, o3_vmr), self._vmr_interpolators = self._inputs(atmo)
        self._wavenumber_interpolator = linear_interpolating_matrix(
            self._mtckd_wavenumbers,
            atmo.wavenumbers_cminv,
            "zero",
        )
        inputs = (
            np.ascontiguousarray(atmo.pressure_pa, dtype=np.float64),
            np.ascontiguousarray(atmo.temperature_k, dtype=np.float64),
            np.ascontiguousarray(h2o_vmr, dtype=np.float64),
            np.ascontiguousarray(co2_vmr, dtype=np.float64),
            np.ascontiguousarray(o3_vmr, dtype=np.float64),
        )
        if atmo.calculate_derivatives:
            self._linearized_cache = mt_ckd_linearized(*inputs, include_rayleigh=False)
            extinction = self._linearized_cache[0]
        else:
            self._linearized_cache = None
            extinction = mt_ckd(*inputs, include_rayleigh=False)
        atmo.storage.total_extinction[:] += (
            np.nan_to_num(extinction) @ self._wavenumber_interpolator.T
        )

    def register_derivative(self, atmo: sk.Atmosphere, name: str) -> None:
        if (
            self._linearized_cache is None
            or self._wavenumber_interpolator is None
            or self._vmr_interpolators is None
        ):
            msg = "MTCKDContinuum.add_to_atmosphere must be called before register_derivative"
            raise RuntimeError(msg)

        altitude_count = len(atmo.model_geometry.altitudes())
        native_grid_interpolator = np.eye(altitude_count)
        derivative_specs = (
            (
                "pressure_pa",
                "altitude",
                native_grid_interpolator,
                "wf_pressure_pa",
            ),
            (
                "temperature_k",
                "altitude",
                native_grid_interpolator,
                "wf_temperature_k",
            ),
            (
                f"{self._h2o_name}_vmr",
                f"{self._h2o_name}_altitude",
                self._vmr_interpolators[0],
                f"wf_{self._h2o_name}_vmr",
            ),
            (
                f"{self._co2_name}_vmr",
                f"{self._co2_name}_altitude",
                self._vmr_interpolators[1],
                f"wf_{self._co2_name}_vmr",
            ),
            (
                f"{self._o3_name}_vmr",
                f"{self._o3_name}_altitude",
                self._vmr_interpolators[2],
                f"wf_{self._o3_name}_vmr",
            ),
        )
        for (derivative_name, interp_dim, interpolator, assign_name), derivative in zip(
            derivative_specs,
            self._linearized_cache[1:],
            strict=True,
        ):
            if (
                derivative_name == "pressure_pa"
                and not atmo.calculate_pressure_derivative
            ):
                continue
            if (
                derivative_name == "temperature_k"
                and not atmo.calculate_temperature_derivative
            ):
                continue

            extinction_derivative = (
                np.nan_to_num(derivative) @ self._wavenumber_interpolator.T
            )
            mapping = atmo.storage.get_derivative_mapping(
                f"wf_{name}_{derivative_name}"
            )
            mapping.d_extinction[:] += extinction_derivative
            ssa_derivative = np.zeros_like(extinction_derivative)
            np.divide(
                -extinction_derivative * atmo.storage.ssa,
                atmo.storage.total_extinction,
                out=ssa_derivative,
                where=atmo.storage.total_extinction != 0.0,
            )
            mapping.d_ssa[:] += ssa_derivative
            mapping.interp_dim = interp_dim
            mapping.interpolator = interpolator
            mapping.assign_name = assign_name
