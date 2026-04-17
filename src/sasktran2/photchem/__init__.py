from re import A

import sasktran2 as sk
import numpy as np
import xarray as xr

from sasktran2.viewinggeo.wrappers import ViewingGeometry


def actinic_flux():
    altitudes_m = np.arange(0.0, 100.0e3, 1000.0)
    cos_sza = 0.5
    bands = [
        {"name": "o3",
         "species": "o3",
         "wavelegth_range_nm": (280, 800),},
        {"name": "o2_762",
         "species": "o2",
         "wavelegth_range_nm": (755, 775),},
        {"name": "o2_689",
         "species": "o2",
         "wavelegth_range_nm": (685, 700),},
    ]

    species_list = ["o3", "o2"]

    config = sk.Config()

    config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
    config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
    config.flux_types = [sk.FluxType.Actinic, sk.FluxType.Upwelling, sk.FluxType.Downwelling]
    config.num_streams = 4
    config.num_forced_azimuth = 1

    geometry = sk.Geometry1D(
        cos_sza,
        0.0,
        6371000.0,
        altitudes_m,
        sk.InterpolationMethod.LinearInterpolation,
        sk.GeometryType.PlaneParallel
    )

    wavel_nm = np.arange(280, 800, 0.01)

    viewing = sk.ViewingGeometry()

    for alt in altitudes_m:
        viewing.add_flux_observer(
            sk.FluxObserverSolar(
                cos_sza, alt
            )
        )

    atmosphere = sk.Atmosphere(geometry, config, wavel_nm, calculate_derivatives=False)

    sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

    optical = {
        "o3": sk.optical.O3DBM(),
        "h2o": sk.optical.AERLineAbsorber("H2O"),
        "o2": sk.optical.AERLineAbsorber("O2")
    }

    atmosphere["o3"] = sk.climatology.mipas.constituent("o3", optical["o3"])
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["h2o"] = sk.climatology.mipas.constituent("h2o", optical["h2o"])
    atmosphere["o2"] = sk.climatology.mipas.constituent("o2", optical["o2"])
    atmosphere["solar"] = sk.constituent.SolarIrradiance(photon_units=True, mode="average")

    engine = sk.Engine(config, geometry, viewing)

    flux = engine.calculate_radiance(atmosphere)

    actinic_flux = flux["actinic_flux"].to_numpy() # (wavelength, altitude)

    output = xr.Dataset(coords={"altitude": altitudes_m, "wavelength": wavel_nm})
    output["actinic_flux"] = (("wavelength", "altitude"), actinic_flux)


    for species in species_list:
        optical_property = optical[species]

        oq = optical_property.atmosphere_quantities(atmosphere)

        output[species + "_xs"] = (("wavelength", "altitude"), oq.extinction.T)

    return output

if __name__ == "__main__":
    actinic_flux()
