---
file_format: mystnb
---

(_users_mie)=
# Mie Scattering Databases

```{code-cell}
import sasktran2 as sk
import numpy as np
```

SASKTRAN2 includes the ability to generate scattering properties for Mie scatterers.  Generally a Mie
scattering atmospheric constituent is described by a particle size distribution, and a refractive index.
Since Mie scattering calculations can be slow, SASKTRAN2 does not contain perform Mie calculations
alongside the radiative transfer calculation, instead databases are generated using internal tools that
can then be re-used.

Let's start by creating a Mie scattering database for sulfate aerosols following a log-normal particle
distribution parameterized by a median radius and mode width.

```{code-cell}
mie_db = sk.database.MieDatabase(
    sk.mie.distribution.LogNormalDistribution(),
    sk.mie.refractive.H2SO4(),
    wavelengths_nm=np.arange(270, 1000, 50.0),
    median_radius=[100, 150, 200],
    mode_width=[1.5, 1.7],
)
```

The main class which provides the database functionality is the {py:class}`sasktran2.database.MieDatabase` object.
When it is created for the first time the local database will be generated. Any subsequent instantiations of the
object will re-use the cached database.  In creating the database we provided a {py:class}`sasktran2.mie.distribution.LogNormalDistribution`
object which defines the particle size distribution for the scatterer, as well as a {py:class}`sasktran2.mie.refractive.H2SO4`
which in this case is the refractive index function for sulfates.  We also specified the wavelengths to create the database
at, as well as the particle size parameters that should be used.

We can look at the created database
```{code-cell}
mie_db.load_ds()
```
Which contains the scattering and absorption cross sections, as well as the Legendre expansion moments of the phase function
that are necessary for the radiative transfer calculation.

## Including Mie Scatterers in the Radiative Transfer Calculation
The Mie database can be included and used as an optical property with the standard {py:class}`sasktran2.constituent.NumberDensityScatterer` and
{py:class}`sasktran2.constituent.ExtinctionScatterer` constituents.

First we set up the example calculation without aerosol

```{code-cell}
import matplotlib.pyplot as plt

config = sk.Config()

altitudes_m = np.arange(0, 65001, 1000)

model_geometry = sk.Geometry1D(cos_sza=0.6,
                               solar_azimuth=0,
                               earth_radius_m=6372000,
                               altitude_grid_m=altitudes_m,
                               interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                               geometry_type=sk.GeometryType.Spherical)

viewing_geo = sk.ViewingGeometry()

for alt in [10000, 20000, 30000, 40000]:
    ray = sk.TangentAltitudeSolar(tangent_altitude_m=alt,
                                    relative_azimuth=45*np.pi/180,
                                    observer_altitude_m=200000,
                                    cos_sza=0.6)
    viewing_geo.add_ray(ray)

wavel = np.arange(280.0, 800.0, 1)

atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

atmosphere['rayleigh'] = sk.constituent.Rayleigh()
atmosphere['ozone'] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())

engine = sk.Engine(config, model_geometry, viewing_geo)

radiance_no_aerosol = engine.calculate_radiance(atmosphere)
```

Now we can add in aerosol, let's set the extinction to be a constant 1e-7 per metre until 30 km,
and then zero above that.

```{code-cell}
aero_ext = np.zeros(len(altitudes_m))
aero_ext[0:30] = 1e-7

atmosphere['aerosol'] = sk.constituent.ExtinctionScatterer(
    mie_db,
    altitudes_m=altitudes_m,
    extinction_per_m=aero_ext,
    extinction_wavelength_nm=745,
    median_radius=np.ones_like(aero_ext)*120,
    mode_width=np.ones_like(aero_ext)*1.6
)

radiance_with_aerosol = engine.calculate_radiance(atmosphere)
```

And we can plot the results

```{code-cell}
plt.plot(
    wavel, radiance_no_aerosol['radiance'].isel(los=0, stokes=0),
    wavel, radiance_with_aerosol['radiance'].isel(los=0, stokes=0)
)
```

If we look at the dataset with aerosol included,

```{code-cell}
radiance_with_aerosol
```

We can see that the model has also calculated the derivatives with respect to the particle size
parameters, `median_radius` and `mode_width`.
