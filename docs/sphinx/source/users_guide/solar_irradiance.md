---
file_format: mystnb
---

(_users_solar_irradiance)=
# Solar Irradiance/Absolute Units

By default SASKTRAN2 will calculate the radiance in units of `/steradian`,
which is essentially the same as assuming the input solar irradiance at
the top of the atmosphere is 1.

```{code-cell}
import sasktran2 as sk
import numpy as np
import matplotlib.pyplot as plt


config = sk.Config()
config.multiple_scatter_source = sk.MultipleScatterSource.TwoStream
config.num_streams = 2

model_geometry = sk.Geometry1D(cos_sza=0.6,
                               solar_azimuth=0,
                               earth_radius_m=6372000,
                               altitude_grid_m=np.arange(0, 65001, 2000),
                               interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                               geometry_type=sk.GeometryType.PseudoSpherical)

viewing_geo = sk.ViewingGeometry()

viewing_geo.add_ray(
    sk.GroundViewingSolar(
        cos_sza=0.6,
        relative_azimuth=0,
        cos_viewing_zenith=0.8,
        observer_altitude_m=200000,
    )
)

wavel = np.arange(200, 1500, 0.01)
atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel, calculate_derivatives=False)

atmosphere["surface"] = sk.constituent.LambertianSurface(0.3)

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

atmosphere["rayleigh"] = sk.constituent.Rayleigh()
atmosphere['ozone'] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())
atmosphere['no2'] = sk.climatology.mipas.constituent("NO2", sk.optical.NO2Vandaele())

engine = sk.Engine(config, model_geometry, viewing_geo)

rad = engine.calculate_radiance(atmosphere)

rad.isel(stokes=0, los=0)["radiance"].plot()
plt.ylabel("Sun Normalized Radiance /ster")
plt.xlabel("Wavelength [nm]")
```

And we see we have a very smooth radiance structure where the ozone and NO2 absorption structures are visible.
Typically we refer to this quantity as "Sun Normalized Radiance" because it is equal to the true radiance divided
by the top of the atmosphere solar irradiance.  To compute absolute radiances, we have to multiply by the solar irradiance.
This could be done outside of the calculation by the user, but SASKTRAN2 also contains the capability to include
the solar irradiance internally as a constitent, {py:class}`sasktran2.constituent.SolarIrradiance`.  Let's add it
to the atmosphere and repeat the calculation,


```{code-cell}
atmosphere["solar_irradiance"] = sk.constituent.SolarIrradiance()

rad = engine.calculate_radiance(atmosphere)

rad.isel(stokes=0, los=0)["radiance"].plot()
plt.ylabel("Radiance W/m^2/ster")
plt.xlabel("Wavelength [nm]")
```

and now the units are in absolute radiances (W/m^2/nm/ster) and we see the standard solar features in the UV/VIS/NIR.

## Default Solar Model
The default solar irradiance model is that of

Coddington, O. M., Richard, E. C., Harber, D., Pilewskie, P., Woods, T. N., Snow, M., et al. (2023). Version 2 of the TSIS-1 Hybrid Solar Reference Spectrum and Extension to the Full Spectrum. Earth and Space Science, 10, e2022EA002637. https://doi.org/10.1029/2022EA002637

which is valid from 200 nm to 200 microns and integrates exactly to the expected solar flux of 1362.8 W/m^2, and is given at an extremely high spectral resolution of 0.001 nm.

## Resolution Options
The default option is to sample the solar irradiance at the native high resolution, which may not be desirable for some applications.
Several other options are provided to average or integrate the solar spectrum across the wavelength grid.  Let's first set up a low resolution
calculation to demonstrate

```{code-cell}
wavel = np.arange(200, 1500, 0.2)
atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel, calculate_derivatives=False)

atmosphere["surface"] = sk.constituent.LambertianSurface(0.3)

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

atmosphere["rayleigh"] = sk.constituent.Rayleigh()
atmosphere['ozone'] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())
atmosphere['no2'] = sk.climatology.mipas.constituent("NO2", sk.optical.NO2Vandaele())
atmosphere["solar_irradiance"] = sk.constituent.SolarIrradiance()

rad_lowres = engine.calculate_radiance(atmosphere)

rad.isel(stokes=0, los=0)["radiance"].plot()
rad_lowres.isel(stokes=0, los=0)["radiance"].plot()
plt.ylabel("Radiance W/m^2/nm/ster")
plt.xlabel("Wavelength [nm]")
plt.xlim(400, 402)
plt.legend(["High res", "Low res"])
```

and we can see that the two curves agree at the low-resolution sample points, but the low-resolution calculation
may not be a good representation of the average solar flux.

### mode="average"
The solar irradiance constituent has the option to average the solar irradiance at the resolution of the
calculation.

```{code-cell}
atmosphere["solar_irradiance"] = sk.constituent.SolarIrradiance(mode="average")

rad_average = engine.calculate_radiance(atmosphere)

rad.isel(stokes=0, los=0)["radiance"].plot()
rad_lowres.isel(stokes=0, los=0)["radiance"].plot()
rad_average.isel(stokes=0, los=0)["radiance"].plot()
plt.ylabel("Radiance W/m^2/nm/ster")
plt.xlabel("Wavelength [nm]")
plt.xlim(400, 402)
plt.legend(["High res", "Low res", "Average"])
```

which provides a much more reasonable representation of the solar spectrum at the resolution of the calculation.
By default the solar spectrum is degraded to match the calculation resolution, i.e., so that integrals over the wavelength
dimension can be done using the calculation grid.
Optionally we can also manually specify the resolution the solar spectrum is degraded to,

```{code-cell}
atmosphere["solar_irradiance"] = sk.constituent.SolarIrradiance(mode="average", resolution=0.1)

rad_res = engine.calculate_radiance(atmosphere)

rad.isel(stokes=0, los=0)["radiance"].plot()
rad_lowres.isel(stokes=0, los=0)["radiance"].plot()
rad_average.isel(stokes=0, los=0)["radiance"].plot()
rad_res.isel(stokes=0, los=0)["radiance"].plot()
plt.ylabel("Radiance W/m^2/ster")
plt.xlabel("Wavelength [nm]")
plt.xlim(400, 402)
plt.legend(["High res", "Low res", "Average", "0.1 nm"])
```

Manually specifying the resolution is most useful when you are doing calculations not over
a wavelength grid.

### mode="integrate"
The `mode="integrate"` option allows the model to integrate over the spectral band rather than take the average
solar irradiance.  Here the output will be in units of radiance (W/m^2/ster) rather than spectral radiance (W/m^2/nm/ster)

```{code-cell}
atmosphere["solar_irradiance"] = sk.constituent.SolarIrradiance(mode="integrate")

rad_int = engine.calculate_radiance(atmosphere)

rad_int.isel(stokes=0, los=0)["radiance"].plot()
plt.ylabel("radiance W/m^2/ster")
plt.xlabel("Wavelength [nm]")
```

This is most useful when you want to integrate over the spectral dimension, because then this can be done just by summing all
of the spectral grid points. For example, we could integrate our originally high resolution spectral radiance,

```{code-cell}
np.trapezoid(rad.isel(stokes=0, los=0)["radiance"].to_numpy(), x=rad["wavelength"].to_numpy())
```

and compare it to our low resolution sum

```{code-cell}
float(rad_int.isel(stokes=0, los=0)["radiance"].sum().to_numpy())
```

And see that they are quite close.
