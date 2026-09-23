---
file_format: mystnb
---

(_users_photochemical_emission)=
# Photochemical Emission
SASKTRAN2 has basic support for some photochemical based emission sources.  This is enabled by setting

```{code-cell}
import sasktran2 as sk

config = sk.Config()
config.emission_source = sk.EmissionSource.VolumeEmissionRate
```

note that it is not currently possible to combine thermal emissions with other photo-chemical based emissions.

## O2 band VER and temperature profiles

`O2BandEmissionRate` accepts the total photon VER of an individual vibrational
band, independently of an excited-state population model. This is useful when
retrieving VER and temperature profiles directly:

```python
# altitude_ver_m, ver_00, and ver_11 are one-dimensional arrays.
atmosphere["o2_00"] = sk.constituent.O2BandEmissionRate(
    altitude_ver_m, ver_00, band="0-0"
)
atmosphere["o2_11"] = sk.constituent.O2BandEmissionRate(
    altitude_ver_m, ver_11, band="1-1"
)
atmosphere.temperature_k = temperature_on_model_grid

# Update these during the retrieval; rebuilding the constituent is unnecessary.
atmosphere["o2_00"].photon_ver = updated_ver_00
result = engine.calculate_radiance(atmosphere)

dI_dVER00 = result["wf_o2_00_photon_ver"]
dI_dVER11 = result["wf_o2_11_photon_ver"]
dI_dT = result["wf_temperature_k"]
```

VER is in photons m^-3 s^-1, integrated over all directions and the selected
band's lines. The constituent supplies the isotropic factor of `1 / (4 pi)`.
Supported bands are A-band `0-0` and `1-1`, and B-band `1-0`. The VER altitude
grid may differ from the atmospheric grid. VER is linearly interpolated first;
rotational line weights and Doppler widths are then evaluated at each model
location's current atmospheric temperature. Outside the VER grid, the default
is zero emission; `out_of_bounds_mode="extend"` uses the nearest endpoint.

Temperature changes the normalized relative line intensities and the Doppler
widths, holding each band VER fixed. The rotational distribution follows the
atmospheric temperature, without imposing a thermal ratio between vibrational
bands. The photon emission integrated over a fully resolved band is therefore
independent of temperature. The weighting functions use the VER input grid
(`o2_00_altitude` and `o2_11_altitude` above), and the atmospheric temperature
grid (`altitude`). Both `einstein_a_branching` (default) and
`hitran_line_strength` line-weight models support these derivatives.

The shared atmospheric temperature derivative holds pressure and gas VMRs fixed
as well as band VER. It includes the resulting gas-density change through the
ideal gas law. A retrieval that updates pressure with temperature must also
apply the chain rule using `wf_pressure_pa`. No photochemical population
derivative is implied by an independent VER parameter.

Add an O2 absorption constituent separately to model self-absorption. With
temperature derivatives enabled, its absorption contribution and the emission
contributions are summed into the same `wf_temperature_k`. Pressure broadening
is included in that absorption calculation; emitted lines currently use
Doppler profiles. Setting `temperature_derivative=False` on the atmosphere
retains VER derivatives and skips emission temperature-derivative evaluation.
Setting `calculate_derivatives=False` skips all derivative registration.

### Spectral sampling and radiance units

Resolve the emission and absorption lines on the model wavelength grid before
convolving the calculated radiance and both Jacobians with an instrument's
spectral response. A 1 nm instrument resolution does not permit a 1 nm model
grid: self-absorption acts on the narrow lines before the instrument averages
them. Validate the grid over the full temperature range allowed by the
retrieval. `AtmosphereIntegratedLineShape` averages atmospheric properties;
it does not replace convolution of the radiance when self-absorption is present.

With `wavelengths_nm`, these sources produce photon radiance in photons
m^-2 s^-1 sr^-1 nm^-1. For a calculation that also includes scattered sunlight,
use `SolarIrradiance(photon_units=True)` so the two contributions have the same
units. Apply any conversion to energy radiance or solar normalization
consistently to the modeled radiance, Jacobians, and observations.

The derivative dimensions describe the input grid but do not automatically
carry its physical altitude coordinates. Assign them before selecting by height:

```python
result = result.assign_coords(
    altitude=atmosphere.model_geometry.altitudes(),
    o2_00_altitude=atmosphere["o2_00"].altitudes_m,
    o2_11_altitude=atmosphere["o2_11"].altitudes_m,
)
```

### Initializing from populations

The population interface converts each upper-state population to band VER using
the corresponding band Einstein-A coefficient, then uses the same band emission
calculation. To initialize an independent VER retrieval from those populations:

```python
population_emission = sk.constituent.PopulationEmissionRate(population_dataset)
bands = population_emission.to_band_emissions()
atmosphere["o2_00"] = bands["0-0"]
atmosphere["o2_11"] = bands["1-1"]
```

These are independent copies; add them instead of the population constituent to
avoid counting emission twice. They can be updated without changing the input
population dataset. The population constituent itself also contributes a
temperature derivative, at fixed supplied populations and fixed band Einstein-A
coefficients; it does not differentiate a photochemical model.

The population interface's inspection arrays (`photon_ver`, `altitudes_m`,
`wavelengths_nm`, `weights`, and the `line_list_*` methods) are read-only views
of the combined A/B-band spectra calculated from the input dataset's temperature.
Attempts to modify them raise `ValueError`; for retrieval updates, modify the
`photon_ver` of a band returned by `to_band_emissions()` instead. Actual source
calculations use the current atmospheric temperature. Use a band's
`line_weights(temperature_k)` method to inspect its normalized line weights at
another temperature.

## Monochromatic Sources
Many photochemical sources in the atmosphere are essentially monochromatic, and can be included by using the
{py:class}`sasktran2.constituent.MonochromaticVolumeEmissionRate` constituent.

```{code-cell}
import sasktran2 as sk
import numpy as np
import matplotlib.pyplot as plt

config = sk.Config()
config.emission_source = sk.EmissionSource.VolumeEmissionRate

model_geometry = sk.Geometry1D(cos_sza=-0.6,
                                solar_azimuth=0,
                                earth_radius_m=6372000,
                                altitude_grid_m=np.arange(0, 120001, 1000),
                                interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                geometry_type=sk.GeometryType.Spherical)

viewing_geo = sk.ViewingGeometry()

for alt in [95000]:
    ray = sk.TangentAltitudeSolar(tangent_altitude_m=alt,
                                    relative_azimuth=0,
                                    observer_altitude_m=200000,
                                    cos_sza=-0.6)
    viewing_geo.add_ray(ray)

wavel = np.arange(556.0, 560.0, 0.01)
atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

# Oxygen green line VER profile
altitude = np.array([
    140, 135, 130, 125, 120, 115, 110,
    105, 100,  95,  90,  85
]).astype(np.float64)[::-1] * 1000

# VER in ph cm^-3 s^-1
VER = np.array([
     400,   600,   800,  1100,  1500,  2000,  2600,
    3200,  3800,  4300,  3600,   900
])[::-1] / 100 # convert to ph cm^-2 m^-1

atmosphere["rayleigh"] = sk.constituent.Rayleigh()
atmosphere["ver"] = sk.constituent.MonochromaticVolumeEmissionRate(altitude, VER, 557.7)

engine = sk.Engine(config, model_geometry, viewing_geo)
output = engine.calculate_radiance(atmosphere)

output["radiance"].isel(los=0).plot()
```

Note that since the calculation must be performed on a finite resolution spectral grid, SASKTRAN internally
"widens" the monochromatic line based on the resolution of the calculation.  This is done so that integrals
over the line produce the correct integrated radiance profile.
