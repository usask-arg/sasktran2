from __future__ import annotations

import numpy as np
import xarray as xr

import sasktran2 as sk

LYMAN_ALPHA_WAVELENGTH_NM = 121.567
LYMAN_ALPHA_TOA_RATE_S = 3.40e-9
LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S = 3.2e15

ACTINIC_FLUX_BASE_WAVELENGTH_RANGE_NM = (120.0, 1280.0)
ACTINIC_FLUX_BASE_RESOLUTION_NM = 0.1
ACTINIC_FLUX_O2_LINE_RESOLUTION_NM = 0.001
ACTINIC_FLUX_O2_LINE_BANDS_NM = (
    (679.0, 699.0),  # O2 B-band, centered near 689 nm.
    (752.0, 776.0),  # O2 A-band, centered near 762 nm.
    (1260.0, 1280.0),  # O2 singlet-delta band, centered near 1.27 um.
)


def _closed_arange(start: float, stop: float, step: float) -> np.ndarray:
    return np.arange(start, stop + step / 2.0, step)


def _actinic_flux_wavelength_grid() -> np.ndarray:
    grid_parts = [
        _closed_arange(
            ACTINIC_FLUX_BASE_WAVELENGTH_RANGE_NM[0],
            ACTINIC_FLUX_BASE_WAVELENGTH_RANGE_NM[1],
            ACTINIC_FLUX_BASE_RESOLUTION_NM,
        ),
        np.array([LYMAN_ALPHA_WAVELENGTH_NM]),
    ]

    grid_parts.extend(
        _closed_arange(start, stop, ACTINIC_FLUX_O2_LINE_RESOLUTION_NM)
        for start, stop in ACTINIC_FLUX_O2_LINE_BANDS_NM
    )

    return np.unique(np.round(np.concatenate(grid_parts), decimals=6))


def actinic_flux():
    altitudes_m = np.arange(0.0, 130.0e3, 1000.0)
    cos_sza = 0.5

    species_list = ["o3", "o2", "n2"]

    config = sk.Config()

    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.flux_types = [
        sk.FluxType.Actinic,
        sk.FluxType.Downwelling,
        sk.FluxType.Upwelling,
    ]
    config.num_streams = 4
    config.num_forced_azimuth = 1
    config.num_threads = 8

    geometry = sk.Geometry1D(
        cos_sza,
        0.0,
        6371000.0,
        altitudes_m,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.PseudoSpherical,
    )

    wavel_nm = _actinic_flux_wavelength_grid()
    lyman_alpha_index = np.where(wavel_nm == LYMAN_ALPHA_WAVELENGTH_NM)[0][0]

    viewing = sk.ViewingGeometry()

    for alt in altitudes_m:
        viewing.add_flux_observer(sk.FluxObserverSolar(cos_sza, alt))

    atmosphere = sk.Atmosphere(geometry, config, wavel_nm, calculate_derivatives=False)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    optical = {
        "o3": sk.optical.O3DBM(),
        "o2": sk.optical.AERLineAbsorber("O2")
        + sk.optical.O2SchumannRunge()
        + sk.optical.O2LymanAlpha(),
        "n2": sk.optical.AERLineAbsorber("N2"),
    }

    atmosphere["o3"] = sk.climatology.mipas.constituent("o3", optical["o3"])
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["o2"] = sk.climatology.mipas.constituent("o2", optical["o2"])
    atmosphere["n2"] = sk.climatology.mipas.constituent("n2", optical["n2"])
    atmosphere["solar"] = sk.constituent.SolarIrradiance(
        photon_units=True, mode="average"
    )

    engine = sk.Engine(config, geometry, viewing)

    flux = engine.calculate_radiance(atmosphere)

    actinic_flux = flux["actinic_flux"].to_numpy()  # (wavelength, altitude)

    output = xr.Dataset(coords={"altitude": altitudes_m, "wavelength": wavel_nm})
    output["actinic_flux"] = (("wavelength", "altitude"), actinic_flux)

    toa_spectral_flux = atmosphere.storage.solar_irradiance[lyman_alpha_index]
    output["lyman_alpha_actinic_flux"] = (
        ["altitude"],
        actinic_flux[lyman_alpha_index, :]
        / toa_spectral_flux
        * LYMAN_ALPHA_TOA_FLUX_PHOTONS_M2_S,
    )

    output["temperature"] = (["altitude"], atmosphere.temperature_k)

    for species in species_list:
        optical_property = optical[species]

        oq = optical_property.atmosphere_quantities(atmosphere)

        output[species + "_xs"] = (("wavelength", "altitude"), oq.extinction.T)

        numden = (
            np.interp(
                altitudes_m, atmosphere[species].altitudes_m, atmosphere[species].vmr
            )
            * atmosphere.state_equation.dry_air_numberdensity["N"]
        )

        output[species + "_density"] = (["altitude"], numden)

    output["o_density"] = sk.climatology.atomic_oxygen.number_density(
        output["altitude"].to_numpy(), np.datetime64("2010-01-01"), -50.0
    )
    output["co2_density"] = (
        400e-6 * atmosphere.state_equation.dry_air_numberdensity["N"]
    )

    return output


if __name__ == "__main__":
    actinic_flux()


from .models import Yankovsky  # noqa: E402

__all__ = ["Yankovsky", "actinic_flux"]
