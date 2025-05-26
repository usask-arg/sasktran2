---
file_format: mystnb
---

(_source_occultation)=
# Occultation Sources
SASKTRAN2 can be used in pure occultation mode, i.e., to calculate transmittances instead of radiances.
To do this, we disable all other sources and enable the occultation source.

```{code-cell}
import sasktran2 as sk

config = sk.Config()

config.single_scatter_source = sk.SingleScatterSource.NoSource
config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
config.emission_source = sk.EmissionSource.NoSource

config.occultation_source = sk.OccultationSource.Standard
```


```{code-cell}
import numpy as np

model_geometry = sk.Geometry1D(cos_sza=0.6,
                                solar_azimuth=0,
                                earth_radius_m=6372000,
                                altitude_grid_m=np.arange(0, 65001, 1000),
                                interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                geometry_type=sk.GeometryType.Spherical)

viewing_geo = sk.ViewingGeometry()

for alt in [10000, 20000, 30000, 40000]:
    ray = sk.TangentAltitudeSolar(tangent_altitude_m=alt,
                                    relative_azimuth=0,
                                    observer_altitude_m=200000,
                                    cos_sza=0.6)
    viewing_geo.add_ray(ray)

wavel = np.arange(280.0, 800.0, 0.01)
atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

atmosphere['rayleigh'] = sk.constituent.Rayleigh()

atmosphere['ozone'] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())
atmosphere['no2'] = sk.climatology.mipas.constituent("NO2", sk.optical.NO2Vandaele())

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

engine = sk.Engine(config, model_geometry, viewing_geo)

output = engine.calculate_radiance(atmosphere)


output['radiance'].isel(los=0).plot()
```

Note that in occultation mode the model produces transmittances, i.e., values between 0 and 1 where 0 is completely attenuated and 1 is no attenuation.

## Relevant Configuration Options

```{eval-rst}
.. autosummary::

  sasktran2.Config.occultation_source

```
