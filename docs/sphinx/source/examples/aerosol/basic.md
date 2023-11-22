---
file_format: mystnb
---

(_example_aerosol_basic)=
# Basic Aerosol Inclusion
```{code-cell}
import sasktran2 as sk
import numpy as np
import matplotlib.pyplot as plt

config = sk.Config()

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

wavel = np.arange(280.0, 800.0, 0.1)
atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

atmosphere['rayleigh'] = sk.constituent.Rayleigh()
atmosphere['ozone'] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())
atmosphere['no2'] = sk.climatology.mipas.constituent("NO2", sk.optical.NO2Vandaele())

engine = sk.Engine(config, model_geometry, viewing_geo)
output_no_aerosol = engine.calculate_radiance(atmosphere)
```
Now we add in aerosol

```{code-cell}
atmosphere["strat_aerosol"] = sk.test_util.scenarios.test_aerosol_constituent(
    model_geometry.altitudes()
)

output_with_aerosol = engine.calculate_radiance(atmosphere)
```

```{code-cell}
plt.plot(
    wavel, output_no_aerosol['radiance'].isel(los=0, stokes=0),
    wavel, output_with_aerosol['radiance'].isel(los=0, stokes=0)
)
```
