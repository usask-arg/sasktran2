---
file_format: mystnb
---

(_example_sodium)=
# High-resolution sodium resonance scattering

This example calculates limb-scattered sunlight across the sodium D2 and D1
lines. We first calculate a Rayleigh background, add a sodium layer, and compare
the resulting spectra at three tangent altitudes. The calculation includes
polarization and the reference solar spectrum.

The atmosphere is an illustrative 70–120 km shell with a 200 K temperature and
a Gaussian sodium layer. We use single scattering to keep the example short;
light scattered by the lower atmosphere and multiple scattering within the layer
are outside this example.

## Resolve the sodium doublet

Use the vacuum line centers from the sodium database. Each window extends
0.02 nm either side of its line, with 0.00005 nm spacing. The two windows are
calculated together, without filling the gap between the D lines.

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt
import sasktran2 as sk

sodium = sk.optical.Sodium()
lines = sodium.database
centers_nm = np.sort(lines.wavelength_nm.values[
    (lines.wavelength_nm.values > 588.0) & (lines.wavelength_nm.values < 591.0)
])
line_names = ["D2", "D1"]
offset_nm = np.linspace(-0.02, 0.02, 801)
windows_nm = centers_nm[:, None] + offset_nm
wavelength_nm = windows_nm.ravel()
print(f"D2: {centers_nm[0]:.7f} nm; D1: {centers_nm[1]:.7f} nm (vacuum)")
```

The decay-data warning concerns other sodium lines in the database. The D-line
return branches used in this example have been audited as complete.

## Set up the atmosphere and viewing rays

The pressure is 1 Pa at 80 km and falls with a 6 km scale height. The sodium
number density peaks at 3 × 10⁹ m⁻³ at 92 km, with a 5 km standard deviation.
These simple profiles make the example reproducible without an external
climatology.

```{code-cell}
altitude_m = np.arange(70_000.0, 120_001.0, 1_000.0)
number_density = 3e9 * np.exp(-0.5 * ((altitude_m - 92_000.0) / 5_000.0) ** 2)
tangent_altitudes_m = np.array([80_000.0, 90_000.0, 100_000.0])

config = sk.Config()
config.num_stokes = 3
config.num_singlescatter_moments = 4
config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
config.delta_m_scaling = False
geometry = sk.Geometry1D(
    cos_sza=0.6,
    solar_azimuth=0.0,
    earth_radius_m=6_372_000.0,
    altitude_grid_m=altitude_m,
    interpolation_method=sk.InterpolationMethod.LinearInterpolation,
    geometry_type=sk.GeometryType.Spherical,
)
viewing = sk.ViewingGeometry()
for tangent in tangent_altitudes_m:
    viewing.add_ray(sk.TangentAltitudeSolar(
        tangent_altitude_m=tangent,
        relative_azimuth=0.4,
        observer_altitude_m=600_000.0,
        cos_sza=0.6,
    ))

atmosphere = sk.Atmosphere(
    geometry, config, wavelengths_nm=wavelength_nm, calculate_derivatives=False
)
atmosphere.temperature_k = np.full_like(altitude_m, 200.0)
atmosphere.pressure_pa = np.exp(-(altitude_m - 80_000.0) / 6_000.0)
atmosphere["rayleigh"] = sk.constituent.Rayleigh()
atmosphere["sun"] = sk.constituent.SolarIrradiance(mode="sample")

fig, ax = plt.subplots(figsize=(4.5, 4), layout="constrained")
ax.plot(number_density / 1e9, altitude_m / 1000, label="Sodium layer")
for tangent in tangent_altitudes_m:
    ax.axhline(tangent / 1000, color="0.6", linestyle=":", linewidth=1)
ax.set(xlabel="Na number density [10⁹ m⁻³]", ylabel="Altitude [km]")
ax.legend()
ax.grid(alpha=0.25)
plt.show()
```

The dotted lines mark the tangent altitudes. A limb ray samples the layer above
its tangent point as well as the density at that point.

## Calculate the background and add sodium

Both calculations use identical atmospheric and solar inputs. Adding the
constituent changes the scattering and attenuation along the rays.

```{code-cell}
engine = sk.Engine(config, geometry, viewing)
background = engine.calculate_radiance(atmosphere)

atmosphere["sodium"] = sk.constituent.NumberDensityScatterer(
    sodium, altitude_m, number_density
)
with_sodium = engine.calculate_radiance(atmosphere)
```

Plot each line in its own window. Solid curves include sodium; dashed curves
show the Rayleigh background. The vertical scale is spectral radiance because
we included `SolarIrradiance` in energy units.

```{code-cell}
fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), sharey=True, layout="constrained")
colors = ["tab:blue", "tab:orange", "tab:green"]
for line_index, (ax, name, center) in enumerate(zip(axes, line_names, centers_nm)):
    window = slice(line_index * len(offset_nm), (line_index + 1) * len(offset_nm))
    for los, (tangent, color) in enumerate(zip(tangent_altitudes_m, colors)):
        ax.plot(
            offset_nm * 1000,
            with_sodium.radiance.isel(wavelength=window, los=los).sel(stokes="I"),
            color=color, label=f"{tangent / 1000:.0f} km",
        )
        ax.plot(
            offset_nm * 1000,
            background.radiance.isel(wavelength=window, los=los).sel(stokes="I"),
            color=color, linestyle="--", linewidth=1,
        )
    ax.set(title=f"Na {name}: {center:.4f} nm", xlabel="Offset from line center [pm]")
    ax.set_yscale("log")
    ax.grid(alpha=0.25)
axes[0].set_ylabel("Spectral radiance [W m⁻² nm⁻¹ sr⁻¹]")
axes[1].legend(title="Tangent altitude")
plt.show()
```

The sodium features are narrow even on this wavelength scale. Their relative
brightness depends on the incident solar absorption spectrum, the different
line strengths and phase functions, and attenuation within the layer.

## Inspect the illumination and polarization

The stored solar irradiance is the reference spectrum sampled at the calculation
wavelengths. A fine metal-line grid does not add resolution to that solar
reference; its sampling remains a limit on the illumination model.

```{code-cell}
fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), layout="constrained")
for line_index, (name, color) in enumerate(zip(line_names, colors)):
    window = slice(line_index * len(offset_nm), (line_index + 1) * len(offset_nm))
    axes[0].plot(offset_nm * 1000, atmosphere.storage.solar_irradiance[window],
                 color=color, label=name)
    stokes = with_sodium.radiance.isel(wavelength=window, los=1)
    linear_polarization = np.hypot(stokes.sel(stokes="Q"), stokes.sel(stokes="U"))
    axes[1].plot(offset_nm * 1000, 100 * linear_polarization / stokes.sel(stokes="I"),
                 color=color, label=name)
axes[0].set(xlabel="Offset from line center [pm]", ylabel="Solar irradiance [W m⁻² nm⁻¹]")
axes[1].set(xlabel="Offset from line center [pm]", ylabel="Linear polarization [%]",
            title="90 km tangent altitude")
for ax in axes:
    ax.legend()
    ax.grid(alpha=0.25)
plt.show()
```

In this fine-structure approximation, D1 resonance scattering is unpolarized;
D2 can produce linear polarization. The total signal also contains polarized
Rayleigh scattering. Hyperfine structure and frequency redistribution are not
included, so these curves illustrate the implemented model rather than a full
high-resolution sodium fluorescence calculation.

An instrument response should be applied to the calculated Stokes radiances,
before forming a polarization ratio. See the
[metal resonance guide](../../users_guide/metal_resonance.md) for cross-section
and phase-function plots, and [Supported metals](../../api/supported_metals.md)
for other available species.
