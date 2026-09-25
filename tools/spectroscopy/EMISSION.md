# Metal emission templates and ordinary-layer candidates

Free VER fitting needs a spectral shape and an amplitude profile. It does not
require a chemical production rate, absolute absorption cross section, or a
conversion from emitted photons to metal abundance. The supplied components are
isotropic and unpolarized. They describe photons integrated over the stored
spectral interval, not unobserved parts of a molecular band.

## Reproduce the local files

Run `build_metal_emission_templates.py --output <database_root>/spectroscopy/metals/emission`
with NumPy, xarray, a NetCDF backend, Pillow and pypdf available. Downloads are
pinned by SHA-256. The builder retains sources and per-file provenance, plus an
image overlay and pixel measurements for the approximate NiO trace.

| File | Source and retained interval | Limitations |
| --- | --- | --- |
| `FeO.nc` | [Unterguggenberger et al. (2017) supplement](https://acp.copernicus.org/articles/17/4177/2017/acp-17-4177-2017-supplement.zip), `plot_data/fig3.dat`, column 12: Gattinger model, 560–719.9 nm, 0.1 nm sampling | Numerical preliminary model with fixed excited populations. Sampling is not spectral accuracy; source reports approximate line positions and model/observation differences. All-zero 720 nm endpoint omitted. |
| `NiO.nc` | [Evans et al. (2011), Fig. 4](https://acp.copernicus.org/articles/11/9595/2011/), isolated short-dash NiO component, 430–670 nm, sampled every 5 nm | Manual raster digitization of a preliminary model already averaged over 5 nm. Approximate reading uncertainty two pixels, wavelength pixel 0.754 nm. Attribution is provisional. |

Both are relative **photon** spectra normalized to unit integral per nm. FeO's
archived column is the normalized model compared to the paper's photon spectrum;
the NiO figure explicitly has photon-brightness units. The source wavelength
medium is unspecified, so the published coordinates are retained and that
uncertainty is explicit. Neither file silently claims vacuum wavelengths.
Both sources carry CC BY 3.0. Retain attribution, modifications and source hashes.

NiO's coarse template cannot recover structure at OSIRIS's finer spectral
resolution. Do not apply its original 5 nm averaging again. Later
[Noll et al. (2024)](https://acp.copernicus.org/articles/24/1143/2024/) observations
do not support a strong blue NiO contribution and raise uncertainties about the
older attribution. A nonzero fitted coefficient is not a species identification.

## Chemical channels reconsidered

| Species/channel | Evidence and remaining spectral requirement |
| --- | --- |
| Na I D lines | Ordinary chemical nightglow; fit independent line VERs. The D2/D1 ratio varies ([Slanger et al. 2005](https://doi.org/10.1029/2005JD006078)). NaO is a chemical intermediate, not evidence of an observed NaO molecular band ([Chapman-mechanism study](https://doi.org/10.1038/356414a0)). |
| K I D lines | Ordinary chemical nightglow; fit independent line VERs ([Noll et al. 2019](https://doi.org/10.1029/2018JD030044)). KO is an intermediate. |
| Li I / LiO | Li I remains an atomic candidate. LiO intermediacy alone does not establish an ambient LiO emission band; no usable molecular photon template found. |
| OFeOH, FeO₂ | Noll2024 discusses energetically possible ambient iron chemiluminescence. OFeOH has suitable low electronic states but no available emission spectrum. The gap is spectral shape. |
| CaOH, O₂CaOH | Real reservoir candidates. [Gómez Martín and Plane (2017)](https://doi.org/10.1021/acsearthspacechem.7b00072) measure laboratory CaO chemiluminescence and CaOH laser fluorescence; this does not identify ambient CaOH nightglow. [Calcium reservoir model](https://doi.org/10.5194/acp-18-14799-2018). |
| MgOH | Ambient reservoir and laboratory UV spectroscopy candidate; no ordinary-layer photon emission template found. [Primary spectroscopy](https://ntrs.nasa.gov/api/citations/20150020854/downloads/20150020854.pdf). |
| CoO | Laboratory ozone chemiluminescence exists ([Burgard et al. 2006](https://doi.org/10.1366/000370206775382730)); no ambient detection established here. Laboratory pressure/source-specific populations cannot define a unique mesospheric shape. |
| NaOH, KOH, AlOH and larger reservoirs | Keep as possible spectroscopy targets; no ordinary-layer visible photon templates established in this search. Artificial releases are outside this survey's scope. |
| AlO, MgO, CaO, TiO | Supplied absorption/elastic-return line lists are useful. Chemical emission requires assumed or fitted excited-state populations; absorption strengths are not automatically emission weights. |

The existing `MonochromaticVolumeEmissionRate` supports free atomic line VERs
with an explicitly chosen emitter mass. The new `SpectralVolumeEmissionRate`
accepts any supplied photon template; new sources can be added independently of
the resonance optical properties.

## Empirical continuum data

Noll2024 provides a numerical continuum decomposition at
[Zenodo 8335837](https://doi.org/10.5281/zenodo.8335837), CC BY 4.0. Its visual
component includes FeO and potentially other emitters. It is useful as an empirical
continuum template but cannot be labeled a pure FeO or OFeOH spectrum. The
[PALACE model](https://gmd.copernicus.org/articles/18/4353/2025/) provides additional
continuum data. These mixed components have not been installed under species names.
