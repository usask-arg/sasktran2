---
file_format: mystnb
---

(_example_amf)=
# Calculating Air Mass Factors
SASKTRAN2 contains built in capability to calculate box air mass factors.

```{note}
Air mass factors in SASKTRAN2 are estimated using the analytic weighting function capability implemented inside the model.
This means they are exact as long as the sources you are using are linearized, which is currently all sources except for
`sk.MultipleScatterSource.SuccessiveOrders`.  If you use the successive orders source the AMFs are only approximate.
```

We can set up a standard nadir viewing
calculation

```{code-cell}
import sasktran2 as sk
import numpy as np
import matplotlib.pyplot as plt

config = sk.Config()
config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates

model_geometry = sk.Geometry1D(cos_sza=0.6,
                               solar_azimuth=0,
                               earth_radius_m=6372000,
                               altitude_grid_m=np.arange(0, 65001, 1000),
                               interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                               geometry_type=sk.GeometryType.Spherical)

viewing_geo = sk.ViewingGeometry()

viewing_geo.add_ray(
    sk.GroundViewingSolar(
        cos_sza=0.6,
        relative_azimuth=0,
        cos_viewing_zenith=0.8,
        observer_altitude_m=200000,
    )
)

wavel = np.array([290, 310, 330, 350, 400, 600])
atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

atmosphere.surface.albedo[:] = 0.3

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

atmosphere["rayleigh"] = sk.constituent.Rayleigh()
atmosphere['ozone'] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())
atmosphere['no2'] = sk.climatology.mipas.constituent("NO2", sk.optical.NO2Vandaele())

engine = sk.Engine(config, model_geometry, viewing_geo)
```

At this point we have a nadir viewing scene setup with Rayleigh scattering, and Ozone + NO2 absorption.
We have added nothing specific to calculating AMFs.  To enable the calculation of AMFs, we add

```{code-cell}
atmosphere["air_mass_factor"] = sk.constituent.AirMassFactor()
```

The {py:class}`sasktran2.constituent.AirMassFactor` is a dummy constituent that does not actually modify the atmospheric state,
but it provides the necessary configuration to enable the calculation of air mass factors.  If we then do the calculation and
look at the results,

```{code-cell}
radiance = engine.calculate_radiance(atmosphere)

print(radiance)
```
We see that there is a `air_mass_factor` field that contains the air mass factors, we can look at a few of them

```{code-cell}
radiance["air_mass_factor"].isel(stokes=0, los=0).sel(wavelength=310).plot(y="altitude")
radiance["air_mass_factor"].isel(stokes=0, los=0).sel(wavelength=350).plot(y="altitude")
radiance["air_mass_factor"].isel(stokes=0, los=0).sel(wavelength=600).plot(y="altitude")

plt.legend(["310 nm", "350 nm", "600 nm"])
plt.ylabel("Altitude [m]")
plt.xlabel("Air Mass Factor")
plt.title("")
```
