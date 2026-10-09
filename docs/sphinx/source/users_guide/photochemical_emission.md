---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
  language: python
mystnb:
  output_stderr: remove
---

(_users_airglow)=
# Photochemical Emission (Airglow)

In the daytime mesosphere, sunlight produces excited states of oxygen and OH that
emit in the UV and visible: the O2 A, B and gamma bands, the O(1S) green line, the
O(1D) red line, and OH A-X resonance fluorescence near 308 nm. These states are not
in local thermodynamic equilibrium, so their populations come from photochemistry
rather than from the temperature.

SASKTRAN2 handles this in two layers:

- {py:mod}`sasktran2.nlte` solves kinetic mechanisms for the steady-state
  populations of excited states and the photon volume emission rate (VER) of each
  transition.
- {py:func}`sasktran2.nlte.add_photochemical_species` does everything needed for a
  radiance calculation in one call: it computes the actinic flux and photolysis
  rates (see {ref}`_users_photolysis`), solves the kinetics, and adds the emission
  to an {py:class}`sasktran2.Atmosphere` as volume emission rate constituents.

## A background atmosphere

We use the same daytime equatorial profile from the CAIRT Extended Reference
Scenarios as in {ref}`_users_photolysis`.

```{code-cell} ipython3
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import sasktran2 as sk

K_BOLTZMANN = 1.380649e-23
COS_SZA = 0.6

ers = sk.climatology.ers.profile(
    month=4, latitude_degrees=0, local_time_hours=9.5,
    species=["O3", "O2", "O", "NO2", "CO2"],
)
ers = ers.sel(altitude_m=ers.altitude_m <= 150e3)
z = ers["altitude_m"].to_numpy()
temperature = ers["temperature_k"].to_numpy()
pressure = ers["pressure_pa"].to_numpy()
air = pressure / (K_BOLTZMANN * temperature)

chemistry = xr.Dataset(
    {
        "temperature_k": ("altitude", temperature),
        "pressure_pa": ("altitude", pressure),
        "O3": ("altitude", ers["o3_mean"].to_numpy() * air),
        "O2": ("altitude", ers["o2_mean"].to_numpy() * air),
        "N2": ("altitude", 0.7808 * air),
        "CO2": ("altitude", ers["co2_mean"].to_numpy() * air),
        "O(3P)": ("altitude", ers["o_mean"].to_numpy() * air),
    },
    coords={"altitude": z},
)
```

## Kinetic mechanisms

A {py:class}`sasktran2.nlte.Mechanism` lists the excited states, the background
species whose densities are inputs, the processes that couple them, and the
per-molecule rates (photolysis and photoexcitation) that must be supplied. Two are
bundled:

- `oxygen`: O(1D), O2(b, v=0-2) and O2(a, v=0-5) from O3 and O2 photolysis and
  resonant absorption in the O2 atmospheric bands, with JPL 19-5 rates and the
  Yankovsky et al. (2019) vibrational product distributions.
- `oxygen_green`: O(1S) and the Barth mechanism precursor, from O2
  photodissociation and O recombination (McDade et al. 1986).

```{code-cell} ipython3
mechanism = sk.nlte.Mechanism.bundled("oxygen")
print(mechanism)
print("states:", mechanism.states)
print("rate inputs:", mechanism.rate_inputs)
```

The rates come from an actinic flux calculation;
{py:func}`sasktran2.photolysis.presets.oxygen_photolysis` gives exactly the rates this
mechanism needs. {py:func}`sasktran2.nlte.solve` then gives the steady-state
densities, every process rate, and the photon VER of each radiative transition.

```{code-cell} ipython3
altitudes = np.arange(0.0, 150_001.0, 5_000.0)
flux = sk.photolysis.ActinicFlux(altitudes).calculate(chemistry, cos_sza=COS_SZA, albedo=0.3)
rates = sk.photolysis.photolysis_rates(
    flux, sk.photolysis.presets.oxygen_photolysis()
)
rates["P_O1D_ION"] = xr.zeros_like(rates["J_O2_SRC"])  # no ionospheric source here

solution = sk.nlte.solve(mechanism, chemistry.interp(altitude=altitudes), rates)

fig, ax = plt.subplots(figsize=(6, 5))
for state in ("O(1D)", "O2(b)", "O2(b, v=1)", "O2(a)"):
    ax.semilogx(solution["density"].sel(state=state), altitudes / 1e3, label=state)
ax.set_xlabel("Number density [m$^{-3}$]")
ax.set_ylabel("Altitude [km]")
ax.set_xlim(1e6, 1e17)
ax.legend();
```

{py:func}`sasktran2.nlte.budget` splits a state's production and loss by process.
Here are the processes that make up at least 10% of the production or loss of
O2(b, v=0) somewhere in the profile:

```{code-cell} ipython3
o2b = sk.nlte.budget(mechanism, solution, "O2(b)")
production = o2b.where(o2b > 0).sum("process")
fig, ax = plt.subplots(figsize=(7, 5))
for process in o2b["process"].values:
    rate = o2b.sel(process=process)
    if float((abs(rate) / production).max()) > 0.1:
        ax.semilogx(abs(rate), altitudes / 1e3, "-" if rate.max() > 0 else "--", label=process)
ax.set_xlabel("Production (solid), loss (dashed) [m$^{-3}$ s$^{-1}$]")
ax.set_ylabel("Altitude [km]")
ax.set_xlim(1e4, None)
ax.legend(fontsize=8);
```

Through the mesosphere O2(b, v=0) comes in comparable parts from resonant
absorption of sunlight in the A band (`o2_excitation_b0`), from O(1D) + O2, and from
collisional relaxation of O2(b, v=1), which is itself made by O(1D) + O2 and by
B-band absorption. Quenching by N2 is the main loss below about 90 km, and A-band
emission above.

## Emission in a limb calculation

{py:func}`sasktran2.nlte.add_photochemical_species` takes a set-up
{py:class}`sasktran2.Atmosphere` and a list of emitting species:

```{code-cell} ipython3
for name, emitter in sk.nlte.PHOTOCHEMICAL_SPECIES.items():
    print(f"{name:6s} {emitter.description}")
```

The background densities for the kinetics come from the `background` dataset and
from the atmosphere's {py:class}`~sasktran2.constituent.VMRAltitudeAbsorber`
constituents; atomic oxygen must be in `background`. Two settings matter:

- the emission only reaches the radiance with
  `config.emission_source = sk.EmissionSource.VolumeEmissionRate`;
- the emission is in photons, so the solar spectrum must be too:
  `sk.constituent.SolarIrradiance(photon_units=True)`.

Line emission is much narrower than the usual model wavelength spacing.
{py:func}`sasktran2.nlte.emission_wavelength_grid` builds a grid that resolves the
O2 bands and the Doppler profile of each line, and is coarse elsewhere.

Here is a helper that sets up a limb calculation over a wavelength window:

```{code-cell} ipython3
model_altitudes = np.arange(0.0, 150_001.0, 1_000.0)
tangents_km = np.arange(30.0, 110.1, 10.0)


def limb_setup(wavelengths):
    config = sk.Config()
    config.emission_source = sk.EmissionSource.VolumeEmissionRate
    config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
    geometry = sk.Geometry1D(
        COS_SZA, 0.0, 6_372_000.0, model_altitudes,
        sk.InterpolationMethod.LinearInterpolation, sk.GeometryType.Spherical,
    )
    viewing = sk.ViewingGeometry()
    for tangent in tangents_km:
        viewing.add_ray(sk.TangentAltitudeSolar(tangent * 1e3, 0.0, 600e3, COS_SZA))

    atmosphere = sk.Atmosphere(geometry, config, wavelengths_nm=wavelengths)
    atmosphere.temperature_k = np.interp(model_altitudes, z, temperature)
    atmosphere.pressure_pa = np.exp(np.interp(model_altitudes, z, np.log(pressure)))
    vmr = lambda key: np.interp(model_altitudes, z, ers[f"{key}_mean"].to_numpy())
    atmosphere["O3"] = sk.constituent.VMRAltitudeAbsorber(sk.optical.O3DBM(), model_altitudes, vmr("o3"))
    atmosphere["O2"] = sk.constituent.VMRAltitudeAbsorber(sk.optical.HITRANAbsorber("O2"), model_altitudes, vmr("o2"))
    atmosphere["rayleigh"] = sk.constituent.Rayleigh()
    atmosphere["solar"] = sk.constituent.SolarIrradiance(photon_units=True)
    return sk.Engine(config, geometry, viewing), atmosphere


background = chemistry[["O(3P)", "CO2"]]
# The actinic flux on a coarser grid than the model, for speed.
actinic_flux = sk.photolysis.ActinicFlux(altitudes)
```

Multiple scattering is switched off to keep this page fast; it matters little for
the emission above the stratosphere, but adds to the scattered background.

### The O2 A band

```{code-cell} ipython3
wavelengths = sk.nlte.emission_wavelength_grid(["O2(b)"], range_nm=(758.0, 772.0))
engine, atmosphere = limb_setup(wavelengths)
without = engine.calculate_radiance(atmosphere)

sk.nlte.add_photochemical_species(
    atmosphere, ["O2(b)"], cos_sza=COS_SZA, background=background, actinic_flux=actinic_flux
)
with_emission = engine.calculate_radiance(atmosphere)

fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for ax, km in zip(axes, (50, 90)):
    ax.plot(wavelengths, without["radiance"].sel(los=list(tangents_km).index(km)), label="no emission")
    ax.plot(wavelengths, with_emission["radiance"].sel(los=list(tangents_km).index(km)), lw=0.8, label="with O2(b)")
    ax.set_title(f"{km} km tangent")
    ax.set_xlabel("Wavelength [nm]")
axes[0].set_ylabel("Radiance [photons m$^{-2}$ s$^{-1}$ nm$^{-1}$ sr$^{-1}$]")
axes[0].set_yscale("log")
axes[0].legend();
```

Without emission the band appears in absorption at 50 km. With it, the
brightest lines are an order of magnitude above the scattered sunlight at 50 km,
and two at 90 km.

### The green and red lines

`"O(1S)"` adds the 557.7 nm green line and its 297.2 nm partner, and `"O(1D)"` the
630.0 and 636.4 nm lines. The `oxygen_green` mechanism needs O2(a) from `oxygen`,
which is solved first. Above about 105 km (O(1S)) and 150 km (O(1D)) ionospheric
sources dominate; supply them with `ionospheric_o1s_production` and
`ionospheric_o1d_production` if they matter.

```{code-cell} ipython3
wavelengths = sk.nlte.emission_wavelength_grid(["O(1S)"], range_nm=(557.0, 558.5))
engine, atmosphere = limb_setup(wavelengths)
green = sk.nlte.add_photochemical_species(
    atmosphere, ["O(1S)"], cos_sza=COS_SZA, background=background, actinic_flux=actinic_flux
)
radiance = engine.calculate_radiance(atmosphere)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].semilogx(green["photon_ver"].sel(transition="o1s_green_line") * 1e-6, altitudes / 1e3)
axes[0].set_xlim(1e-1, None)
axes[0].set_xlabel("557.7 nm VER [photons cm$^{-3}$ s$^{-1}$]")
axes[0].set_ylabel("Altitude [km]")
for i in (5, 7, 8):
    axes[1].plot(wavelengths, radiance["radiance"].sel(los=i), label=f"{tangents_km[i]:.0f} km")
axes[1].set_xlim(557.84, 557.94)
axes[1].set_xlabel("Wavelength [nm]")
axes[1].legend()
fig.tight_layout();
```

### OH fluorescence and self-absorption

`"OH(A)"` adds solar resonance fluorescence of OH A-X, computed line by line from
HITRAN, from an OH density in `background`. By default it also adds OH as a line
absorber, so that OH absorbs both the scattered sunlight and its own emission along
the line of sight. The ERS profiles have no OH, so we use a simple daytime-like
profile peaking near 75 km.

```{code-cell} ipython3
oh = 1.0e13 * np.exp(-0.5 * ((z - 75e3) / 12e3) ** 2) + 1.0e10 * np.exp(-z / 20e3)
background_oh = chemistry[["O(3P)", "CO2"]].assign(OH=("altitude", oh))

wavelengths = sk.nlte.emission_wavelength_grid(["OH(A)"], range_nm=(306.0, 311.0))
engine, atmosphere = limb_setup(wavelengths)
without = engine.calculate_radiance(atmosphere)["radiance"]
sk.nlte.add_photochemical_species(
    atmosphere, ["OH(A)"], cos_sza=COS_SZA, background=background_oh, actinic_flux=actinic_flux
)
with_oh = engine.calculate_radiance(atmosphere)["radiance"]

i = list(tangents_km).index(70)
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(wavelengths, without.sel(los=i), label="no OH")
ax.plot(wavelengths, with_oh.sel(los=i), lw=0.6, label="OH fluorescence and absorption")
ax.set_xlim(306.5, 310.0)
ax.set_xlabel("Wavelength [nm]")
ax.set_ylabel("Radiance [photons m$^{-2}$ s$^{-1}$ nm$^{-1}$ sr$^{-1}$]")
ax.set_title("70 km tangent")
ax.legend();
```

Pass `oh_self_absorption=False` to add only the emission.

## Putting it together

All of the species can be added at once, and
{py:func}`~sasktran2.nlte.emission_wavelength_grid` with the default range gives a
280-800 nm grid for them, of about 76 000 wavelengths. The script
`tools/nlte/limb_scan_demo.py` in the repository runs a full OSIRIS-like scan this
way, in about two minutes:

```python
species = ["O2(b)", "OH(A)", "O(1S)", "O(1D)"]
wavelengths = sk.nlte.emission_wavelength_grid(species)
...
solution = sk.nlte.add_photochemical_species(
    atmosphere, species, cos_sza=cos_sza, background=background
)
```

The returned dataset holds the solution of every mechanism involved, with the rates
that drove it; `sk.nlte.budget` works on it directly.
