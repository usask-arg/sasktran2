---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
  language: python
mystnb:
  output_stderr: remove
---

(_users_photolysis)=
# Photolysis and Actinic Flux

{py:mod}`sasktran2.photolysis` computes the actinic flux, the spherically
integrated radiance that drives photochemistry, with the SASKTRAN2
discrete-ordinates engine, and integrates it against cross sections and quantum
yields to give photolysis and photoexcitation rates, in the spirit of the TUV model.

The default calculation covers 80-1280 nm. It resolves the O2 band structure that
matters for the mesosphere: the O2 A, B, gamma and 1.27 µm bands at 0.001 nm, the
Schumann-Runge bands at 0.002 nm, and the vacuum ultraviolet below 122 nm at
0.01 nm. Lyman-alpha O2 photolysis follows the Chabrillat and Kockarts (1997)
slant-column parameterisation.

## The atmosphere

Photolysis works on an `xarray.Dataset` on an `altitude` coordinate in metres, with
`temperature_k`, `pressure_pa` and number densities [m^-3] named by species id. The
same dataset convention is used by {py:func}`sasktran2.nlte.solve`. Here we take a
daytime equatorial profile from the CAIRT Extended Reference Scenarios.

```{code-cell} ipython3
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import sasktran2 as sk

ers = sk.climatology.ers.profile(
    month=4, latitude_degrees=0, local_time_hours=9.5,
    species=["O3", "O2", "O", "NO2", "CO2"],
)
z = ers["altitude_m"].to_numpy()
air = ers["pressure_pa"].to_numpy() / (1.380649e-23 * ers["temperature_k"].to_numpy())
atmosphere = xr.Dataset(
    {
        "temperature_k": ("altitude", ers["temperature_k"].to_numpy()),
        "pressure_pa": ("altitude", ers["pressure_pa"].to_numpy()),
        **{
            name: ("altitude", ers[f"{key}_mean"].to_numpy() * air)
            for name, key in [("O3", "o3"), ("O2", "o2"), ("NO2", "no2"), ("CO2", "co2"), ("O(3P)", "o")]
        },
        "N2": ("altitude", 0.7808 * air),
    },
    coords={"altitude": z},
)
atmosphere = atmosphere.sel(altitude=atmosphere.altitude <= 150e3)
```

## Actinic flux

{py:class}`sasktran2.photolysis.ActinicFlux` sets up the engine with flux observers
at every altitude. The absorbers are the species of the dataset that have default
cross sections (O3, O2, N2, NO2 and atomic O); Rayleigh scattering and a Lambertian
surface are always included. The O2 line absorption comes from HITRAN, downloaded on
first use, which needs the optional `hitran-api` package (`pip install sasktran2[hapi]`).

```{code-cell} ipython3
altitudes = np.arange(0.0, 150_001.0, 5_000.0)
calculator = sk.photolysis.ActinicFlux(altitudes)
flux = calculator.calculate(atmosphere, cos_sza=0.6, albedo=0.3)
flux
```

The ratio of the actinic flux to the top-of-atmosphere solar flux shows where
sunlight is absorbed, and how much the scattered light adds:

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(9, 4))
ratio = flux["actinic_flux"] / flux["solar_flux"]
for km in (20, 50, 80, 110):
    ratio.sel(altitude=km * 1e3).sel(wavelength=slice(100, 700)).plot(ax=ax, label=f"{km} km")
ax.set_yscale("log")
ax.set_ylim(1e-6, 3)
ax.set_xlabel("Wavelength [nm]")
ax.set_ylabel("Actinic flux / solar flux")
ax.set_title("")
ax.legend();
```

## Photolysis rates

{py:func}`sasktran2.photolysis.photolysis_rates` integrates flux × cross section ×
quantum yield for a list of reactions. {py:func}`sasktran2.photolysis.presets.oxygen_photolysis`
gives the rates that drive the bundled `oxygen` mechanism: O3 photolysis split into
the O(1D) and O(3P) channels (with the Matsumi et al. 2002 O(1D) yield), O2
photolysis in the Schumann-Runge continuum and at Lyman-alpha, and resonant
excitation of O2(b) and O2(a) in the atmospheric bands.

```{code-cell} ipython3
rates = sk.photolysis.photolysis_rates(flux, sk.photolysis.presets.oxygen_photolysis())

fig, ax = plt.subplots(figsize=(6, 5))
for name in ("J_O3_O1D", "J_O3_O3P", "J_O2_SRC", "J_O2_LYA", "J_O2_EXC_B0"):
    ax.semilogx(rates[name], rates["altitude"] / 1e3, label=name)
ax.set_xlabel("Rate [s$^{-1}$]")
ax.set_ylabel("Altitude [km]")
ax.set_xlim(1e-12, 1e-1)
ax.legend();
```

Any absorber in the flux dataset can be used for a reaction. A quantum yield is a
constant or a function of wavelength [nm] and temperature [K]; here NO2 photolysis
with a unit yield below its 398 nm threshold:

```{code-cell} ipython3
no2 = sk.photolysis.Photolysis(
    "J_NO2", "NO2", lambda wavelength, temperature: 1.0 * (wavelength < 398.0)
)
j_no2 = sk.photolysis.photolysis_rates(flux, [no2])["J_NO2"]
print(j_no2.sel(altitude=[0, 30e3, 60e3]).values)
```

## TUV mode

{py:class}`sasktran2.photolysis.TUVActinicFlux` runs the same engine at the
resolution of TUV-x v5.4: its 156 wavelength bins, homogeneous layers, and its O2
parameterisations (Chabrillat and Kockarts at Lyman-alpha, Koppers and Murtagh in
the Schumann-Runge bands). With `data="tuv-x"` it also uses TUV-x's own cross
sections, solar spectrum and quantum yields. It takes a fraction of a second, and
is useful for comparisons with TUV or where the full resolution is not needed.

```{code-cell} ipython3
tuv_flux = sk.photolysis.TUVActinicFlux(altitudes, data="tuv-x").calculate(
    atmosphere, cos_sza=0.6, albedo=0.3
)
tuv_rates = sk.photolysis.photolysis_rates(tuv_flux, sk.photolysis.presets.tuvx_v54_photolysis())

fig, ax = plt.subplots(figsize=(6, 5))
ax.semilogx(rates["J_O3_O1D"], altitudes / 1e3, label="full resolution")
ax.semilogx(tuv_rates["J_O3_O1D"], altitudes / 1e3, "--", label="TUV mode, TUV-x data")
ax.set_xlabel("O3 -> O(1D) [s$^{-1}$]")
ax.set_ylabel("Altitude [km]")
ax.legend();
```

In this example the TUV mode with TUV-x data is 9-13% higher above 50 km, and up to
25% higher between 30 and 45 km. Above the ozone layer most of the difference comes
from the solar spectrum: TUV-x's extraterrestrial flux is 10-14% above TSIS-1 HSRS in
the Hartley band. With `data="sasktran2"` the TUV mode uses SASKTRAN2's cross
sections and solar spectrum averaged over the TUV bins.
