---
file_format: mystnb
---

(_users_mie)=
# Mie Scattering Databases

```{code-cell}
import sasktran2 as sk
import numpy as np
import matplotlib.pyplot as plt
```

SASKTRAN2 includes the ability to generate scattering properties for Mie scatterers.  Generally a Mie
scattering atmospheric constituent is described by a particle size distribution, and a refractive index.
Since Mie scattering calculations can be slow, SASKTRAN2 does not contain perform Mie calculations
alongside the radiative transfer calculation, instead databases are generated using internal tools that
can then be re-used.

Mie scatterers are parameterized by two quantities, the refractive index, and a particle size distribution.
SASKTRAN2 contains several built in refractive index functions (see [Refractive Index Containers](../api/mie.md#refractive-index-containers)), let's create one for sulfates,
{py:class}`sasktran2.mie.refractive.H2SO4`

```{code-cell}
refrac = sk.mie.refractive.H2SO4()

wavelengths = np.arange(280, 2000, 0.1)
plt.plot(wavelengths, np.real(refrac.refractive_index(wavelengths)))
plt.xlabel("Wavelength [nm]")
plt.ylabel("Real Part of Refractive Index")
```
it is also possible to define your own refractive index entirely using the {py:class}`sasktran2.mie.refractive.RefractiveIndex`
object.

The second quantity that we need is a particle size distribution (see [Particle Size Distributions](../api/mie.md#particle-size-distributions)).  Here we will create a log-normal particle size distribution

```{code-cell}
distribution = sk.mie.distribution.LogNormalDistribution()
```

Every particle size distribution is parameterized by several parameters, we can look at them explicity,

```{code-cell}
distribution.args()
```

And we can also inspect the distribution for specific values of the parameters (median_radius, mode_width)

```{code-cell}
radii = np.arange(1, 500, 0.1)
plt.plot(radii, distribution.distribution(median_radius=160, mode_width=1.6).pdf(radii))
plt.xlabel("Radius [nm]")
plt.ylabel("PDF")
```

Now that we have a particle size distribution and a refractive index, we can create our Mie scattering database.

```{code-cell}
mie_db = sk.database.MieDatabase(
    distribution,
    refrac,
    wavelengths_nm=np.arange(270, 1000, 50.0),
    median_radius=[100, 150, 200],
    mode_width=[1.5, 1.7],
)
```

The database is a function of wavelength, and any arguments of the particle size distribution, `median_radius` and `mode_width` in our case.
When it is created for the first time the local database will be generated. Any subsequent instantiations of the
object will re-use the cached database.


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
and then zero above that.  When we specify our Mie scatterer, we have to make sure to also specify
any arguments of the particle size distribution (`median_radius` and `mode_width`) as a function of `altitudes_m`.

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


## Freezing Distribution Parameters
The above method of creating Mie scattering databases is useful if you want to run calculations over a wide range of the
particle size distribution parameters (`median_radius` and `mode_width`) in our example.  But there are many applications
where some or all of these values are fixed and known.  In this case we can "freeze" the parameters of the distribution,

```{code-cell}
distribution = sk.mie.distribution.LogNormalDistribution().freeze(median_radius=160, mode_width=1.6)
```

Now when we look at the `args` property of the distribution,

```{code-cell}
distribution.args()
```

We see that it is empty, the distribution is not a function of any parameters because they are all frozen.

We can construct the Mie scattering database and inspect it,

```{code-cell}
mie_db = sk.database.MieDatabase(
    distribution,
    refrac,
    wavelengths_nm=np.arange(270, 1000, 50.0),
)

mie_db.load_ds()
```

We can see that the Mie properties are no longer a function of the particle size distribution parameters since
they are assumed to be known. In addition when constructing the database we did not have to specify them.

Similarly, to add aerosol in to our SASKTRAN2 atmosphere we no longer have to specify the frozen particle
size parameters,

```{code-cell}
atmosphere['aerosol'] = sk.constituent.ExtinctionScatterer(
    mie_db,
    altitudes_m=altitudes_m,
    extinction_per_m=aero_ext,
    extinction_wavelength_nm=745,
)
```

In this example we froze all of the distribution parameters, but we could have frozen only one of them

```{code-cell}
sk.mie.distribution.LogNormalDistribution().freeze(mode_width=1.6).args()
```

Using this distribution we would have to specify `median_radius` when constructing the database and when
adding scatterers to the atmsophere.
