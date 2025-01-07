---
file_format: mystnb
---

(_users_emissions)=
# Emissions

To include any type of emissions as part of the radiation source, they must first be enabled in the configuration.

```{code-cell}
import sasktran2 as sk

config = sk.Config()
config.emission_source = sk.EmissionSource.Standard
```

By default, {py:attr}`sasktran2.Config.emission_source` is set to `sk.EmissionSource.NoSource`, where emissions are not added to
the source and any emissions present in the model atmosphere will be ignored.

## Thermal Emissions

After enabling emissions in the configuration, we must also specify the type of emissions to include in the
atmosphere. Currently, thermal emissions from the surface and atmosphere are supported. First we will
enable emission sources in the configuration and disable the single scatter source.

```{code-cell}
import matplotlib.pyplot as plt
import numpy as np
import sasktran2 as sk

config = sk.Config()
config.single_scatter_source = sk.SingleScatterSource.NoSource
config.emission_source = sk.EmissionSource.Standard
```

Now, set up a line of sight that ends at the ground.

```{code-cell}
model_geometry = sk.Geometry1D(cos_sza=0.6,
                                solar_azimuth=0,
                                earth_radius_m=6372000,
                                altitude_grid_m=np.arange(0, 100001, 1000),
                                interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                geometry_type=sk.GeometryType.Spherical)

viewing_geo = sk.ViewingGeometry()

ray = sk.GroundViewingSolar(0.6, 0, 0.8, 200000)
viewing_geo.add_ray(ray)
```

Next, we will add an absorbing/emitting species to the atmosphere. For this example we will use
methane, with the cross sections calculated using the HITRAN database. Note that the cross section
calculation will take several minutes the first time it runs but the result is saved in a
database cache in order to speed up future calculations. The `hitran-api` python package must
be installed to download the database files.

```{code-cell}
:tags: ["remove-cell"]
wavelengths = np.arange(7370, 7380, 0.01)
hitran_db = sk.optical.database.OpticalDatabaseGenericAbsorber(sk.database.StandardDatabase().path("hitran/CH4/sasktran2/6f006bcb051f81fc57d1bd09315589bfe77b4348.nc"))
```

```{code}
wavelengths = np.arange(7370, 7380, 0.01)

hitran_db = sk.database.HITRANDatabase(
    molecule="CH4",
    start_wavenumber=1355,
    end_wavenumber=1357,
    wavenumber_resolution=0.01,
    reduction_factor=1,
    backend="sasktran2",
    profile="voigt"
)
```

Now we will create the atmopshere object and enable both line of sight and surface emissions.
For the ground emission, we need to specify the surface temperature and emissivity.
Wavelength-dependent emissivities can be specified by passing an array of values to the
{py:class}`sasktran2.constituent.SurfaceThermalEmission` constructor. The array must be the same length as the
wavelength array given to the {py:class}`sasktran2.Atmosphere` object.

```{code-cell}
atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavelengths, calculate_derivatives=False)

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)
atmosphere['ch4'] = sk.climatology.mipas.constituent('CH4', hitran_db)
atmosphere['emission'] = sk.constituent.ThermalEmission()
atmosphere['surface_emission'] = atmosphere['surface_emission'] = sk.constituent.SurfaceThermalEmission(temperature_k=300, emissivity=0.9)
```

Lastly, perform the calculation and plot the result.

```{code-cell}
engine = sk.Engine(config, model_geometry, viewing_geo)
rad = engine.calculate_radiance(atmosphere)

rad.isel(stokes=0)["radiance"].plot()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Radiance [W/m^2/nm/ster]')

```

## Combining Thermal Emissions and Solar Irradiance

If we wish to perform calculations that include both thermal emissions and scattered sunlight, we must take care
to ensure that the units of the thermal and solar sources are in agreement. By default, the solar
source terms are normalized by the irradiance at the top of the atmosphere, but thermal emissions
are input to sasktran with absolute radiance units (W/m^2/nm/ster). To combine these sources we must
add the solar irradiance constituent to the atmosphere which changes the units to match the thermal emissions.

```{code-cell}
atmosphere["solar_irradiance"] = sk.constituent.SolarIrradiance()
```

See [Solar Irradiance](_users_solar_irradiance) documentation for more details.
