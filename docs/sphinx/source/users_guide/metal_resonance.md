---
file_format: mystnb
---

(_users_metal_resonance)=
# Metal resonance scattering

To add a metal layer, choose an optical property such as
{py:class}`sasktran2.optical.Sodium` and give it a number-density profile using
{py:class}`sasktran2.constituent.NumberDensityScatterer`. The optical property
supplies the line cross sections and scattering phase function; the constituent
supplies the amount of metal at each altitude.

This page uses sodium to show what those inputs look like. All wavelengths are
**vacuum wavelengths in nm**. The familiar sodium D lines near 589.0 and 589.6 nm
in air are near 589.1583 and 589.7558 nm in vacuum.

## Choose an optical property

Creating the sodium property downloads its spectroscopy on first use and caches
it in the standard SASKTRAN2 database. Here we read the two D-line centers from
that spectroscopy rather than rounding them to a plotting grid.

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt
import sasktran2 as sk

sodium = sk.optical.Sodium()
lines = sodium.database
d_lines = lines.where(
    (lines.wavelength_nm > 588.0) & (lines.wavelength_nm < 591.0), drop=True
).sortby("wavelength_nm")
centers_nm = d_lines.wavelength_nm.to_numpy()
line_names = ["D2", "D1"]

for name, center in zip(line_names, centers_nm):
    print(f"{name}: {center:.7f} nm (vacuum)")
```

The decay-data warning from `Sodium()` concerns other lines in its database;
the two D lines used here have audited complete return branches.

Other species work the same way: for example, `sk.optical.Potassium()` or
`sk.optical.MagnesiumIon()`. The [supported metals table](../api/supported_metals.md)
lists the available names. To use your own file, pass `db_filepath` to the
optical-property constructor; `db` accepts an in-memory `xarray.Dataset`.

## Plot the line cross sections

Metal lines are much narrower than a typical broadband wavelength grid. We use
a spacing of 0.00002 nm here and compare three temperatures. Increasing the
temperature broadens each line and lowers its peak.

```{code-cell}
offset_nm = np.linspace(-0.004, 0.004, 401)
fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True, layout="constrained")

for ax, center, name in zip(axes, centers_nm, line_names):
    for temperature in [150.0, 200.0, 300.0]:
        properties = sodium.cross_sections(
            center + offset_nm,
            altitudes_m=np.array([90_000.0]),
            temperature_k=np.array([temperature]),
        )
        scattering = properties.extinction[0] * properties.ssa[0]
        ax.plot(offset_nm * 1000, scattering, label=f"{temperature:.0f} K")
    ax.set(title=f"Na {name}: {center:.4f} nm", xlabel="Offset from line center [pm]")
    ax.grid(alpha=0.25)
axes[0].set_ylabel("Scattering cross section [m² atom⁻¹]")
axes[1].legend()
plt.show()
```

`cross_sections()` returns arrays shaped `(altitude, wavelength)`.
`extinction` is the total cross section in m², and `ssa` is the fraction that
returns by scattering at the same wavelength. Their product is the scattering
cross section. For the sodium D lines this fraction is essentially one; other
species can emit some of the absorbed light at different wavelengths.

## Plot the scattering direction dependence

The two D lines also have different phase functions. D1 is isotropic in this
fine-structure model, while D2 scatters more strongly in the forward and backward
directions. The phase function below has a spherical average of one.

```{code-cell}
from sasktran2.optical.resonance import resonance_phase_function

at_center = sodium.cross_sections(
    centers_nm,
    altitudes_m=np.array([90_000.0]),
    temperature_k=np.array([200.0]),
)
angle_deg = np.linspace(0, 180, 361)
fig, ax = plt.subplots(figsize=(6, 3.4), layout="constrained")
for name, w2 in zip(line_names, at_center.polarizability[0]):
    phase = resonance_phase_function(np.cos(np.deg2rad(angle_deg)), w2)
    ax.plot(angle_deg, phase, label=name)
ax.set(xlabel="Scattering angle [degrees]", ylabel="Phase function", ylim=(0.8, 1.3))
ax.legend()
ax.grid(alpha=0.25)
plt.show()
```

SASKTRAN2 uses the corresponding phase matrix automatically in a radiance
calculation. Set `config.num_stokes = 3` when polarization is needed.

## Add a sodium layer to an atmosphere

Here is a complete atmosphere setup for a simple mesospheric layer. The Gaussian
profile peaks at 92 km. The temperature and pressure are illustrative profiles;
replace them with your atmospheric state for a particular observation.

```{code-cell}
altitude_m = np.arange(70_000.0, 120_001.0, 1_000.0)
config = sk.Config()
config.num_singlescatter_moments = 4
config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
geometry = sk.Geometry1D(
    cos_sza=0.6,
    solar_azimuth=0.0,
    earth_radius_m=6_372_000.0,
    altitude_grid_m=altitude_m,
    interpolation_method=sk.InterpolationMethod.LinearInterpolation,
    geometry_type=sk.GeometryType.Spherical,
)
atmosphere = sk.Atmosphere(
    geometry,
    config,
    wavelengths_nm=centers_nm[0] + offset_nm,
    calculate_derivatives=False,
)
atmosphere.temperature_k = np.full_like(altitude_m, 200.0)
atmosphere.pressure_pa = np.exp(-(altitude_m - 80_000.0) / 6_000.0)
atmosphere["rayleigh"] = sk.constituent.Rayleigh()

number_density = 3e9 * np.exp(-0.5 * ((altitude_m - 92_000.0) / 5_000.0) ** 2)
atmosphere["sodium"] = sk.constituent.NumberDensityScatterer(
    sodium, altitude_m, number_density
)

fig, ax = plt.subplots(figsize=(4.5, 4), layout="constrained")
ax.plot(number_density / 1e9, altitude_m / 1000)
ax.set(xlabel="Na number density [10⁹ m⁻³]", ylabel="Altitude [km]")
ax.grid(alpha=0.25)
plt.show()
```

The number density refers to neutral sodium atoms. The optical property uses
`atmosphere.temperature_k` to calculate thermal line broadening and lower-state
populations. It does not calculate how much sodium is ionized or chemically
bound in molecules. For a molecular property, supply the density of its named
isotopologue.

To finish this setup with viewing rays, a solar spectrum, and radiance plots,
follow the [high-resolution sodium example](../examples/metals/sodium.md).

## Choosing a useful calculation

Resolve the metal lines before applying an instrument response to the calculated
radiance. A coarse wavelength grid can miss the lines entirely. The default
line-wing cutoff is ±0.1 nm; `line_wing_cutoff_nm=None` keeps the complete Voigt
wings.

This model describes scattering back to the absorbing lower state. It does not
include fluorescence between different lines, frequency redistribution,
hyperfine/isotope splitting, or magnetic effects. Molecular band spectra can
therefore have emission that this scattering calculation does not represent.

For derivatives, include the physical background scattering, as above: a fitted
density cannot start at exactly zero *total* scattering. Temperature supplied
through the atmosphere contributes to `wf_temperature_k`; an explicit
`temperature_k` profile on the constituent gives `wf_<name>_temperature_k`.

See the [model and data-format reference](../api/metal_resonance.md) for the
equations and detailed assumptions, and
[Supported metals and spectroscopy sources](../api/supported_metals.md) for
the complete available inventory.
