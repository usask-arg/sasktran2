---
file_format: mystnb
---

(_users_flux)=
# Calculating Fluxes
Fluxes can be internally calcluated in `sasktran2` by specifying a location to calculate
the fluxes at with the {py:meth}`sasktran2.ViewingGeometry.add_flux_observer` method.

```{note}
Currently flux output is only supported when using the `DiscreteOrdinates` source, in `PlaneParallel` mode, i.e.

```python
config = sk.Config()
config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates

model_geometry = sk.Geometry1D(...,
                               geometry_type=sk.GeometryType.PlaneParallel)


```

Flux calculations are setup the exact same way as standard radiance calculations, except we tell the
the model to output fluxes instead of radiances. For example,


```{code-cell}
import sasktran2 as sk
import numpy as np
import matplotlib.pyplot as plt

config = sk.Config()
config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates

config.num_streams = 2
config.delta_m_scaling = True
config.num_forced_azimuth = 1

model_geometry = sk.Geometry1D(cos_sza=0.6,
                               solar_azimuth=0,
                               earth_radius_m=6372000,
                               altitude_grid_m=np.arange(0, 65001, 1000),
                               interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                               geometry_type=sk.GeometryType.PlaneParallel)

viewing_geo = sk.ViewingGeometry()

# Add a flux observer at every layer boundary
for alt in model_geometry.altitudes():
    viewing_geo.add_flux_observer(
        sk.FluxObserverSolar(0.6, alt)
    )

wavel = np.arange(280, 2000, 1.0)
atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

atmosphere["surface"] = sk.constituent.LambertianSurface(0.3)

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

atmosphere["rayleigh"] = sk.constituent.Rayleigh()
atmosphere['ozone'] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())
atmosphere['no2'] = sk.climatology.mipas.constituent("NO2", sk.optical.NO2Vandaele())

atmosphere['solar'] = sk.constituent.SolarIrradiance(mode="average")

engine = sk.Engine(config, model_geometry, viewing_geo)

rad = engine.calculate_radiance(atmosphere)

# TOA Flux
rad["upwelling_flux"].isel(flux_location=-1).plot()
```

By default both upwelling and downwelling fluxes are calculated for every "flux observer" location.

```{code-cell}
print(rad)
```

In addition, derivatives of the flux with respect to the input atmospheric parameters are calculated.
E.g. we can look at how changes in ozone influence the TOA outgoing spectral flux,

```{code-cell}
rad["wf_ozone_vmr_upwelling_flux"].isel(flux_location=-1).plot(y="ozone_altitude", cmap="Blues_r")
```

Or how the surface albedo influences the net flux at TOA,

```{code-cell}
(rad["wf_surface_albedo_upwelling_flux"] - rad["wf_surface_albedo_downwelling_flux"]).isel(flux_location=-1).plot(x="wavelength")
```


```{note}
By default, the flux calculation includes the downwelling flux from the direct solar beam. To disable this
you can set

```python
config.single_scatter_source = sk.SingleScatterSource.NoSource

```
