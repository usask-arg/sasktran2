"""Actinic flux from the SASKTRAN2 discrete-ordinates engine."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import xarray as xr

import sasktran2 as sk
from sasktran2.constants import K_BOLTZMANN

LYMAN_ALPHA_WAVELENGTH_NM = 121.567

#: Windows [nm] resolved at the line resolution of
#: :func:`airglow_wavelength_grid`, where O2 line structure matters for both
#: attenuation and photoexcitation.
O2_LINE_WINDOWS_NM = {
    "b-X(0,0) A band": (752.0, 776.0),
    "b-X(1,0) B band": (675.0, 705.0),
    "b-X(2,0) gamma band": (626.0, 634.0),
    "a-X(0,0) 1.27 um band": (1260.0, 1280.0),
}

#: O2 Schumann-Runge bands [nm], resolved at the band resolution of
#: :func:`airglow_wavelength_grid` to match the 0.5 cm^-1 cross sections of
#: :class:`sasktran2.optical.O2UV`.
O2_SCHUMANN_RUNGE_BANDS_NM = (175.4, 204.1)

#: Vacuum-ultraviolet range [nm], resolved at the VUV resolution of
#: :func:`airglow_wavelength_grid` to sample the O2 window structure and the
#: N2 bands.
VUV_RANGE_NM = (80.0, 121.9)

#: Solar spectrum: WHI 2008 quiet Sun below 116 nm, TSIS-1 HSRS above.
SOLAR_SOURCE = "solar_irradiance_whi2008_hsrs_composite"


EARTH_RADIUS_M = 6371000.0


def slant_columns(
    altitude_m,
    densities_m3,
    cos_sza: float,
    earth_radius_m: float = EARTH_RADIUS_M,
    num_points: int = 2000,
) -> np.ndarray:
    """Columns [m^-2] along the straight path to the sun from each altitude.

    ``densities_m3`` is (species, altitude) on ``altitude_m``; densities vary
    exponentially between levels and vanish above the top level. For
    ``cos_sza < 0`` the path passes through its tangent point; paths that
    meet the Earth's surface have an infinite column. Refraction is ignored.
    """
    altitude_m = np.asarray(altitude_m, dtype=float)
    log_density = np.log(np.maximum(np.atleast_2d(densities_m3), 1e-300))
    r_top = earth_radius_m + altitude_m[-1]
    u = np.linspace(0.0, 1.0, num_points)

    def column(s, r0):
        r = np.sqrt(r0**2 + s**2 + 2.0 * r0 * s * cos_sza)
        h = np.clip(r - earth_radius_m, altitude_m[0], altitude_m[-1])
        n = np.exp(np.array([np.interp(h, altitude_m, ld) for ld in log_density]))
        return np.trapezoid(n, s, axis=1)

    columns = np.empty((log_density.shape[0], altitude_m.size))
    for i, z in enumerate(altitude_m):
        r0 = earth_radius_m + z
        b = r0 * cos_sza
        s_top = -b + np.sqrt(b * b + r_top**2 - r0**2)
        if cos_sza >= 0.0:
            # Densest at the start; cluster samples there.
            columns[:, i] = column(s_top * u**2, r0)
        elif r0 * np.sqrt(1.0 - cos_sza**2) < earth_radius_m:
            columns[:, i] = np.inf
        else:
            # Down to the tangent point, then up; cluster at the tangent.
            s_tangent = -b
            columns[:, i] = column(s_tangent * (1.0 - (1.0 - u) ** 2), r0) + column(
                s_tangent + (s_top - s_tangent) * u**2, r0
            )
    return columns


def _closed_arange(start: float, stop: float, step: float) -> np.ndarray:
    return np.arange(start, stop + step / 2.0, step)


def airglow_wavelength_grid(
    range_nm: tuple[float, float] = (80.0, 1280.0),
    resolution_nm: float = 0.1,
    line_resolution_nm: float = 0.001,
    band_resolution_nm: float = 0.002,
    vuv_resolution_nm: float = 0.01,
) -> np.ndarray:
    """Wavelength grid [nm] for photolysis and O2 photoexcitation.

    ``resolution_nm`` spacing over ``range_nm``, ``line_resolution_nm``
    inside :data:`O2_LINE_WINDOWS_NM`, ``band_resolution_nm`` over
    :data:`O2_SCHUMANN_RUNGE_BANDS_NM`, ``vuv_resolution_nm`` over
    :data:`VUV_RANGE_NM`, plus Lyman-alpha exactly.
    """
    parts = [
        _closed_arange(*range_nm, resolution_nm),
        np.array([LYMAN_ALPHA_WAVELENGTH_NM]),
    ]
    lo, hi = O2_SCHUMANN_RUNGE_BANDS_NM
    if lo >= range_nm[0] and hi <= range_nm[1]:
        parts.append(_closed_arange(lo, hi, band_resolution_nm))
    lo, hi = max(VUV_RANGE_NM[0], range_nm[0]), min(VUV_RANGE_NM[1], range_nm[1])
    if lo < hi:
        parts.append(_closed_arange(lo, hi, vuv_resolution_nm))
    parts.extend(
        _closed_arange(lo, hi, line_resolution_nm)
        for lo, hi in O2_LINE_WINDOWS_NM.values()
        if lo >= range_nm[0] and hi <= range_nm[1]
    )
    return np.unique(np.round(np.concatenate(parts), decimals=6))


#: Factories for the cross sections of each absorber, keyed by species id.
#: They are built only for species present in an atmosphere.
DEFAULT_OPTICAL_PROPERTIES = {
    "O3": sk.optical.O3DBM,
    "O2": lambda: sk.optical.HITRANAbsorber("O2")
    + sk.optical.O2UV()
    + sk.optical.O2LymanAlpha()
    + sk.optical.VUVAbsorber("O2"),
    "N2": lambda: sk.optical.VUVAbsorber("N2"),
    "O(3P)": lambda: sk.optical.VUVAbsorber("O"),
    "NO2": sk.optical.NO2Vandaele,
}


def default_optical_properties() -> dict:
    """The default cross sections of every absorber, keyed by species id."""
    return {species: make() for species, make in DEFAULT_OPTICAL_PROPERTIES.items()}


class ActinicFlux:
    """Spherically averaged actinic flux on an altitude grid.

    Uses SASKTRAN2's discrete-ordinates source in pseudo-spherical geometry
    with flux observers at every altitude: Rayleigh scattering, the
    absorbers present in the atmosphere, and a Lambertian surface.

    Parameters
    ----------
    altitudes_m
        Model and output altitude grid [m].
    wavelengths_nm
        Wavelength grid [nm]; :func:`airglow_wavelength_grid` by default.
    num_streams
        Discrete-ordinates streams.
    optical_properties
        Cross sections by species id, replacing entries of
        :data:`DEFAULT_OPTICAL_PROPERTIES`.
    num_threads
        Engine threads.
    """

    def __init__(
        self,
        altitudes_m,
        wavelengths_nm=None,
        num_streams: int = 4,
        optical_properties: Mapping | None = None,
        num_threads: int = 8,
    ):
        self.altitudes_m = np.asarray(altitudes_m, dtype=float)
        self.wavelengths_nm = (
            airglow_wavelength_grid()
            if wavelengths_nm is None
            else np.asarray(wavelengths_nm, dtype=float)
        )
        self.num_streams = num_streams
        self.optical_properties = dict(optical_properties or {})
        self.num_threads = num_threads

    def _optical_property(self, species: str):
        if species not in self.optical_properties:
            self.optical_properties[species] = DEFAULT_OPTICAL_PROPERTIES[species]()
        return self.optical_properties[species]

    def calculate(
        self,
        atmosphere: xr.Dataset,
        cos_sza: float,
        albedo: float = 0.0,
        earth_sun_distance_au: float = 1.0,
    ) -> xr.Dataset:
        """Actinic flux for one solar zenith angle.

        Parameters
        ----------
        atmosphere
            On an ``altitude`` coordinate [m]: ``temperature_k``,
            ``pressure_pa``, and number densities [m^-3] named by species id.
            Species with an entry in ``optical_properties`` absorb; others
            are ignored.
        cos_sza
            Cosine of the solar zenith angle.
        albedo
            Lambertian surface albedo.
        earth_sun_distance_au
            Scales the solar spectrum by its inverse square.

        Returns
        -------
        xr.Dataset
            ``actinic_flux`` (wavelength, altitude) and the top-of-atmosphere
            ``solar_flux`` (wavelength) [photons m^-2 s^-1 nm^-1];
            ``cross_section`` (species, wavelength, altitude) [m^2] of the
            absorbers used; their ``slant_column`` (species, altitude)
            [m^-2] along the solar path; ``temperature_k`` (altitude).
        """
        config = sk.Config()
        config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
        config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
        config.flux_types = [sk.FluxType.Actinic]
        config.num_streams = self.num_streams
        config.num_forced_azimuth = 1
        config.num_threads = self.num_threads

        geometry = sk.Geometry1D(
            cos_sza,
            0.0,
            6371000.0,
            self.altitudes_m,
            sk.InterpolationMethod.LinearInterpolation,
            sk.GeometryType.PseudoSpherical,
        )
        viewing = sk.ViewingGeometry()
        for altitude in self.altitudes_m:
            viewing.add_flux_observer(sk.FluxObserverSolar(cos_sza, altitude))

        atmo = sk.Atmosphere(
            geometry, config, self.wavelengths_nm, calculate_derivatives=False
        )
        source_altitude = atmosphere["altitude"].to_numpy()
        temperature = atmosphere["temperature_k"].to_numpy()
        pressure = atmosphere["pressure_pa"].to_numpy()
        atmo.temperature_k = np.interp(self.altitudes_m, source_altitude, temperature)
        atmo.pressure_pa = np.exp(
            np.interp(self.altitudes_m, source_altitude, np.log(pressure))
        )

        total_density = pressure / (K_BOLTZMANN * temperature)
        known = [
            *DEFAULT_OPTICAL_PROPERTIES,
            *(
                s
                for s in self.optical_properties
                if s not in DEFAULT_OPTICAL_PROPERTIES
            ),
        ]
        absorbers = [s for s in known if s in atmosphere]
        for species in absorbers:
            vmr = atmosphere[species].to_numpy() / total_density
            atmo[species] = sk.constituent.VMRAltitudeAbsorber(
                self._optical_property(species), source_altitude, vmr
            )
        atmo["rayleigh"] = sk.constituent.Rayleigh()
        atmo["surface"] = sk.constituent.LambertianSurface(albedo)
        atmo["solar"] = sk.constituent.SolarIrradiance(
            photon_units=True, mode="average", source=SOLAR_SOURCE
        )

        engine = sk.Engine(config, geometry, viewing)
        radiance = engine.calculate_radiance(atmo)

        scale = 1.0 / earth_sun_distance_au**2
        cross_section = (
            np.stack(
                [
                    self._optical_property(s).atmosphere_quantities(atmo).extinction.T
                    for s in absorbers
                ]
            )
            if absorbers
            else np.zeros((0, self.wavelengths_nm.size, self.altitudes_m.size))
        )
        densities = (
            np.zeros((0, self.altitudes_m.size))
            if not absorbers
            else np.array(
                [
                    np.exp(
                        np.interp(
                            self.altitudes_m,
                            source_altitude,
                            np.log(np.maximum(atmosphere[s].to_numpy(), 1e-300)),
                        )
                    )
                    for s in absorbers
                ]
            )
        )
        ds = xr.Dataset(
            {
                "actinic_flux": (
                    ("wavelength", "altitude"),
                    radiance["actinic_flux"].to_numpy() * scale,
                ),
                "solar_flux": (
                    ("wavelength",),
                    np.array(atmo.storage.solar_irradiance) * scale,
                ),
                "cross_section": (("species", "wavelength", "altitude"), cross_section),
                "slant_column": (
                    ("species", "altitude"),
                    (
                        slant_columns(self.altitudes_m, densities, cos_sza)
                        if absorbers
                        else np.zeros((0, self.altitudes_m.size))
                    ),
                ),
                "temperature_k": (("altitude",), np.array(atmo.temperature_k)),
            },
            coords={
                "altitude": self.altitudes_m,
                "wavelength": self.wavelengths_nm,
                "species": absorbers,
            },
            attrs={
                "cos_sza": cos_sza,
                "albedo": albedo,
                "earth_sun_distance_au": earth_sun_distance_au,
                "num_streams": self.num_streams,
            },
        )
        ds["actinic_flux"].attrs["units"] = "photons m^-2 s^-1 nm^-1"
        ds["solar_flux"].attrs["units"] = "photons m^-2 s^-1 nm^-1"
        ds["cross_section"].attrs["units"] = "m^2"
        ds["slant_column"].attrs["units"] = "m^-2"
        return ds
