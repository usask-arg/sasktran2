from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.constants as const
import xarray as xr
from scipy.integrate import trapezoid

import sasktran2 as sk
from sasktran2.database import MetalSpectroscopyDatabase
from sasktran2.util.interpolation import linear_interpolating_matrix

from .base import Constituent


class SpectralVolumeEmissionRate(Constituent):
    """An isotropic, unpolarized source with a fixed photon spectral shape.

    ``photon_ver`` is the total number of photons emitted per cubic metre per
    second into all directions, integrated over the *entire supplied template*.
    The template is linearly interpolated and normalized once over that support;
    restricting the atmosphere's wavelength range does not renormalize the VER.
    Neither a chemical excitation model nor an emitter number density is needed.

    Parameters
    ----------
    altitudes_m : np.ndarray
        Strictly increasing profile altitudes in metres.
    photon_ver : np.ndarray
        Nonnegative band-integrated photon VER at each profile altitude.
    wavelengths_nm : np.ndarray
        Strictly increasing vacuum wavelengths of the template in nm.
    photon_spectrum : np.ndarray
        Nonnegative relative photon spectral density per nm. At least two
        samples and a positive integral are required. Zero outside its support.
    out_of_bounds_mode : str
        Altitude interpolation: ``"zero"`` (default) or ``"extend"``.
    emission_units : str
        ``"photons"`` yields photon radiance. ``"energy"`` converts each photon
        using hc/lambda, for radiance compatible with the default solar source.
        Wavelength and wavenumber atmosphere grids use densities per nm and
        per cm^-1, respectively. Only monochromatic spectral grids are supported.

    Notes
    -----
    The derivative is ``wf_<name>_photon_ver``. A fixed template does not describe
    changes in rotational/vibrational excitation, spectral redistribution or
    instrumental convolution. Resolve its structure with the atmosphere grid,
    then apply the instrument line shape to the calculated radiance.
    """

    def __init__(
        self,
        altitudes_m: np.ndarray,
        photon_ver: np.ndarray,
        wavelengths_nm: np.ndarray,
        photon_spectrum: np.ndarray,
        *,
        out_of_bounds_mode: str = "zero",
        emission_units: str = "photons",
    ):
        self._altitudes_m = self._grid(altitudes_m, "altitudes_m", minimum=1)
        self._wavelengths_nm = self._grid(wavelengths_nm, "wavelengths_nm")
        if np.any(self._wavelengths_nm <= 0):
            msg = "Template wavelengths must be positive"
            raise ValueError(msg)
        self.photon_ver = photon_ver
        spectrum = self._nonnegative(photon_spectrum, "photon_spectrum")
        if spectrum.shape != self._wavelengths_nm.shape:
            msg = "photon_spectrum must match wavelengths_nm"
            raise ValueError(msg)
        integral = trapezoid(spectrum, self._wavelengths_nm)
        if not np.isfinite(integral) or integral <= 0:
            msg = "photon_spectrum must have a positive finite integral"
            raise ValueError(msg)
        self._photon_spectrum = spectrum / integral
        if out_of_bounds_mode not in ("zero", "extend"):
            msg = "out_of_bounds_mode must be 'zero' or 'extend'"
            raise ValueError(msg)
        if emission_units not in ("photons", "energy"):
            msg = "emission_units must be 'photons' or 'energy'"
            raise ValueError(msg)
        self._out_of_bounds_mode = out_of_bounds_mode
        self._emission_units = emission_units

    @staticmethod
    def _grid(values, name, minimum=2):
        array = np.array(values, dtype=float, copy=True)
        if (
            array.ndim != 1
            or len(array) < minimum
            or not np.all(np.isfinite(array))
            or np.any(np.diff(array) <= 0)
        ):
            msg = f"{name} must be a finite, strictly increasing 1D grid"
            raise ValueError(msg)
        return array

    @staticmethod
    def _nonnegative(values, name):
        array = np.array(values, dtype=float, copy=True)
        if array.ndim != 1 or not np.all(np.isfinite(array)) or np.any(array < 0):
            msg = f"{name} must be a finite, nonnegative 1D array"
            raise ValueError(msg)
        return array

    @property
    def photon_ver(self) -> np.ndarray:
        return self._photon_ver.copy()

    @photon_ver.setter
    def photon_ver(self, values):
        values = self._nonnegative(values, "photon_ver")
        if values.shape != self._altitudes_m.shape:
            msg = "photon_ver must match altitudes_m"
            raise ValueError(msg)
        self._photon_ver = values

    @property
    def altitudes_m(self) -> np.ndarray:
        return self._altitudes_m.copy()

    @property
    def wavelengths_nm(self) -> np.ndarray:
        return self._wavelengths_nm.copy()

    @property
    def photon_spectrum(self) -> np.ndarray:
        """Normalized photon spectral density per nm on the template grid."""
        return self._photon_spectrum.copy()

    def _source(self, atmo):
        if atmo._config.emission_source != sk.EmissionSource.VolumeEmissionRate:
            msg = "Set config.emission_source to VolumeEmissionRate"
            raise ValueError(msg)
        if atmo.spectral_integration_mode != sk.SpectralGridMode.Monochromatic:
            msg = "SpectralVolumeEmissionRate requires a monochromatic grid"
            raise ValueError(msg)
        if atmo.spectral_coordinate is None:
            msg = "Supply atmosphere wavelengths or wavenumbers"
            raise ValueError(msg)
        wavelength = atmo.wavelengths_nm
        if not np.all(np.isfinite(wavelength)) or np.any(wavelength <= 0):
            msg = "Atmosphere wavelengths must be positive and finite"
            raise ValueError(msg)
        shape = np.interp(
            wavelength, self._wavelengths_nm, self._photon_spectrum, left=0, right=0
        )
        if self._emission_units == "energy":
            shape *= const.h * const.c / (wavelength * 1e-9)
        if atmo.spectral_coordinate == "wavenumber_cminv":
            shape *= wavelength**2 / 1e7
        interpolator = linear_interpolating_matrix(
            self._altitudes_m, atmo._native_altitudes(), self._out_of_bounds_mode
        )
        return interpolator, shape / (4 * np.pi)

    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        interpolator, shape = self._source(atmo)
        atmo.storage.emission_source[:] += (interpolator @ self._photon_ver)[
            :, np.newaxis
        ] * shape

    def register_derivative(self, atmo: sk.Atmosphere, name: str):
        interpolator, shape = self._source(atmo)
        full_name = f"wf_{name}_photon_ver"
        mapping = atmo.storage.get_derivative_mapping(full_name)
        # Native output assembly expects both optical derivative arrays even
        # for a constituent that contributes only to the emission source.
        mapping.d_extinction[:] = 0
        mapping.d_ssa[:] = 0
        mapping.d_emission[:] = shape[np.newaxis, :]
        mapping.interpolator = interpolator
        mapping.interp_dim = f"{name}_altitude"
        mapping.assign_name = full_name


class MetalVolumeEmissionRate(SpectralVolumeEmissionRate):
    """Free VER using a cached metal emission template.

    Templates can represent a provisional spectral identification, and may be
    much coarser than an instrument's resolution. Inspect ``metadata`` before
    interpreting a fitted VER. The source wavelength convention and uncertainty
    are retained there; an unspecified convention is not converted to vacuum.
    See :class:`SpectralVolumeEmissionRate` for units and normalization.
    """

    def __init__(
        self,
        species: str,
        altitudes_m: np.ndarray,
        photon_ver: np.ndarray,
        *,
        db: MetalSpectroscopyDatabase | None = None,
        db_filepath: str | Path | None = None,
        out_of_bounds_mode: str = "zero",
        emission_units: str = "photons",
    ):
        if db is not None and db_filepath is not None:
            msg = "Supply either db or db_filepath, not both"
            raise ValueError(msg)
        path = (
            db_filepath
            if db_filepath is not None
            else (db if db is not None else MetalSpectroscopyDatabase()).path(
                species, kind="emission"
            )
        )
        with xr.open_dataset(path) as source:
            dataset = source.load()
        medium = dataset.attrs.get("wavelength_medium")
        if medium not in ("vacuum", "unspecified"):
            msg = "Emission templates require vacuum or unspecified wavelengths"
            raise ValueError(msg)
        self._metadata = dict(dataset.attrs)
        super().__init__(
            altitudes_m,
            photon_ver,
            dataset["wavelength_nm"].values,
            dataset["photon_spectrum"].values,
            out_of_bounds_mode=out_of_bounds_mode,
            emission_units=emission_units,
        )

    @property
    def metadata(self) -> dict:
        """Source provenance, support, wavelength convention and limitations."""
        return self._metadata.copy()


class FeOVolumeEmissionRate(MetalVolumeEmissionRate):
    """Free VER for the FeO orange-band photon template."""

    def __init__(self, altitudes_m, photon_ver, **kwargs):
        super().__init__("FeO", altitudes_m, photon_ver, **kwargs)


class NiOVolumeEmissionRate(MetalVolumeEmissionRate):
    """Free VER for a provisional NiO spectral template; attribution is uncertain."""

    def __init__(self, altitudes_m, photon_ver, **kwargs):
        super().__init__("NiO", altitudes_m, photon_ver, **kwargs)
