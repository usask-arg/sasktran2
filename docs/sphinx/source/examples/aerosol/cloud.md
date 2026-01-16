---
file_format: mystnb
---

(_example_aerosol_cloud)=
# Gaussian Cloud

Clouds can be included in the atmosphere with `sk.constituent.GaussianHeightExtinction` which
creates a gaussian shaped extinction profile. The extinction profile is defined by a height, width, and
total optical depth that are input by the user. An optical property must also be provided to the
constituent constructor, along with any inputs needed by the optical property.

```{code-cell}
import sasktran2 as sk
import numpy as np
import matplotlib.pyplot as plt

config = sk.Config()

altitude_grid_m = np.arange(0, 65001, 1000)

model_geometry = sk.Geometry1D(cos_sza=0.6,
                                solar_azimuth=0,
                                earth_radius_m=6372000,
                                altitude_grid_m=altitude_grid_m,
                                interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                geometry_type=sk.GeometryType.Spherical)

viewing_geo = sk.ViewingGeometry()

for alt in np.arange(0, 60000, 1000):
    ray = sk.TangentAltitudeSolar(tangent_altitude_m=alt,
                                    relative_azimuth=0,
                                    observer_altitude_m=200000,
                                    cos_sza=0.6)
    viewing_geo.add_ray(ray)

wavel = np.array([310, 330, 550, 600])
atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

atmosphere['rayleigh'] = sk.constituent.Rayleigh()
atmosphere['ozone'] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())
atmosphere['no2'] = sk.climatology.mipas.constituent("NO2", sk.optical.NO2Vandaele())
```

For this example we specify the cloud optical properties using a Mie
scattering optical database. This optical property requires that we input
particle size through the `lognormal_median_radius` argument.

```{code-cell}
mie = sk.optical.database.OpticalDatabaseGenericScatterer(
    sk.database.StandardDatabase().path("cross_sections/mie/sulfate_test.nc")
)
atmosphere['cloud'] = sk.constituent.GaussianHeightExtinction(
    optical_property=mie,
    height_m=15500,
    width_fwhm_m=4000,
    vertical_optical_depth=0.1,
    vertical_optical_depth_wavel_nm=550,
    altitudes_m=altitude_grid_m,
    lognormal_median_radius=np.ones_like(altitude_grid_m) * 105.,
)

engine = sk.Engine(config, model_geometry, viewing_geo)
output = engine.calculate_radiance(atmosphere)
```

```{code-cell}
_ = output['radiance'].isel(stokes=0).plot.line(y='los')
```

Weighting functions are calculated for height, width, and optical depth as well as for any
additional inputs to the optical property. For example, the weighting function for the
cloud width can be plotted as follows:

```{code-cell}
_ = output['wf_cloud_width_fwhm_m'].isel(stokes=0, cloud_width_fwhm_m=0).plot.line(y='los')
```
