(_users_metal_emission)=
# Retrieving metal emission rates

A chemical nightglow feature can be retrieved as a free volume emission rate
(VER), without knowing the reaction rate, excited-state population or metal
abundance. {py:class}`sasktran2.constituent.SpectralVolumeEmissionRate` combines a
fixed photon spectrum with a retrieved altitude profile. Each component supplies
an isotropic, unpolarized emission source, so it does not need an absorption cross
section or a resonance phase function. Atmospheric scattering and absorption can
still act on those photons through the other constituents.

## Local FeO and provisional NiO templates

{py:class}`sasktran2.constituent.FeOVolumeEmissionRate` loads the numerical FeO
model archived in the [Unterguggenberger et al. (2017) supplement](https://acp.copernicus.org/articles/17/4177/2017/).
It covers **560–719.9 nm**, sampled every 0.1 nm. That sampling is not a claim of
0.1 nm model accuracy: the source notes approximate line positions and differences
between model and observation. The upper-state populations are fixed by this
template; changing the VER does not change its spectrum.

{py:class}`sasktran2.constituent.NiOVolumeEmissionRate` loads an **approximate
digitization**, every 5 nm over **430–670 nm**, of the isolated NiO component in
[Evans et al. (2011), Fig. 4](https://acp.copernicus.org/articles/11/9595/2011/).
The source was already averaged over 5 nm. This is a coarse exploratory fitting
component, with the original figure and trace overlay retained for inspection.
It cannot recover the missing fine spectral structure. The attribution to NiO
is **provisional**: [Noll et al. (2024)](https://acp.copernicus.org/articles/24/1143/2024/)
find the NiO spectral and atmospheric interpretation problematic. A fitted
amplitude alone does not identify NiO.

Both publications leave the air/vacuum convention of these particular templates
unspecified. The files preserve the published coordinates and explicitly record
that uncertainty; no undocumented air-to-vacuum conversion is applied. Inspect
`constituent.metadata` for provenance, resolution and limitations. This differs
from the resonance line databases, whose coordinates are vacuum wavelengths.

```python
import numpy as np
import sasktran2 as sk

config = sk.Config()
config.emission_source = sk.EmissionSource.VolumeEmissionRate
# Use this config when constructing the atmosphere and engine.

altitude_m = np.arange(70_000.0, 110_001.0, 1_000.0)
photon_ver = 1e6 * np.exp(-0.5 * ((altitude_m - 90_000.0) / 4_000.0) ** 2)
feo = sk.constituent.FeOVolumeEmissionRate(
    altitude_m, photon_ver, emission_units="energy"
)
atmosphere["feo"] = feo
# The radiance result includes wf_feo_photon_ver on feo_altitude.
```

`photon_ver` is in **photons m⁻³ s⁻¹ integrated over the stored band and all 4π
steradians**. It is not the photon yield over unobserved wavelengths, nor the
number density of FeO or NiO. The constituent divides by 4π internally.
`emission_units="energy"` converts photons with hc/λ for consistency with the
default solar irradiance units. The default `"photons"` gives photon radiance;
any other source added in that case must use matching photon units.

Normalization is over the entire template, independently of the wavelengths
used in a fit. Fitting only part of the band therefore retains the same VER
definition. The template is zero outside its supplied support; that convention
does not claim that the physical species has no emission elsewhere. On a
wavelength grid the source density is per nm; on a wavenumber grid the appropriate
Jacobian gives a density per cm⁻¹. Only monochromatic spectral grids are supported.

The free profile can be updated through `feo.photon_ver`. Its analytic weighting
function includes interpolation from that profile grid to the atmosphere grid,
including a zero starting profile. The template itself remains fixed. Fit
spectral-shape uncertainty separately if it matters to the intended retrieval.

## Other chemical emitters and custom spectra

For another published, measured or provisional photon template, use
`SpectralVolumeEmissionRate(altitudes_m, photon_ver, wavelengths_nm, photon_spectrum)`.
The relative spectrum must be nonnegative, per nm, and on an increasing vacuum
wavelength grid. It is normalized internally. Use separate constituents for
independently retrieved components. There is no requirement to supply chemistry
or oscillator strengths for these source terms.

Atomic nightglow such as Na D emission can instead use independent
{py:class}`sasktran2.constituent.MonochromaticVolumeEmissionRate` constituents at
the vacuum line centers. With `line_shape="doppler"`, supply the emitting atom's
mass through `emitter_molecular_weight_g_per_mol`; do not assume a fixed doublet
ratio without evidence. Its generic `ver` parameter can be supplied in photon
units, provided all sources use consistent units.

The [candidate survey](metal_candidates.md) records other plausible ordinary-layer
emitters and distinguishes a missing spectral shape from missing chemistry.
The local templates are under `spectroscopy/metals/emission/`, prepared with
`tools/spectroscopy/build_metal_emission_templates.py`. Source data, licenses,
checksums and modifications accompany the hosted files. Default template loaders
fetch missing files through `StandardDatabase` and retain them in the local cache.
