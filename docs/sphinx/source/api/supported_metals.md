(_api_supported_metals)=
# Supported metals and spectroscopy sources

The standard database supplies the following atomic, molecular, and emission
datasets. Named optical properties are in `sasktran2.optical`; any listed atomic
key can also be passed to {py:class}`sasktran2.optical.AtomicResonance`, for example
`AtomicResonance("Ag_I")`. I denotes a neutral atom and II a singly ionized atom.

## Atomic resonance properties

All 47 atomic/ionic datasets use [NIST ASD 5.12](https://physics.nist.gov/asd),
including its transition-specific source references. The prepared files cover
selected lines in **274–810 nm vacuum** and temperatures up to **1000 K**.
These are cold-atmosphere subsets, not complete high-temperature spectra.

| Element | Available species keys | Named optical properties |
| --- | --- | --- |
| Silver | `Ag_I` | — |
| Aluminium | `Al_I` | {py:class}`Aluminium <sasktran2.optical.Aluminium>` |
| Barium | `Ba_I`, `Ba_II` | {py:class}`Barium <sasktran2.optical.Barium>` (I) |
| Beryllium | `Be_I`, `Be_II` | — |
| Bismuth | `Bi_I` | — |
| Calcium | `Ca_I`, `Ca_II` | {py:class}`Calcium <sasktran2.optical.Calcium>`, {py:class}`CalciumIon <sasktran2.optical.CalciumIon>` |
| Cadmium | `Cd_I` | — |
| Cobalt | `Co_I` | {py:class}`Cobalt <sasktran2.optical.Cobalt>` |
| Chromium | `Cr_I` | {py:class}`Chromium <sasktran2.optical.Chromium>` |
| Caesium | `Cs_I` | — |
| Copper | `Cu_I` | {py:class}`Copper <sasktran2.optical.Copper>` |
| Iron | `Fe_I`, `Fe_II` | {py:class}`Iron <sasktran2.optical.Iron>`, {py:class}`IronIon <sasktran2.optical.IronIon>` |
| Gallium | `Ga_I` | — |
| Germanium | `Ge_I` | — |
| Hafnium | `Hf_I`, `Hf_II` | — |
| Indium | `In_I` | — |
| Potassium | `K_I` | {py:class}`Potassium <sasktran2.optical.Potassium>` |
| Lithium | `Li_I` | {py:class}`Lithium <sasktran2.optical.Lithium>` |
| Magnesium | `Mg_I`, `Mg_II` | {py:class}`Magnesium <sasktran2.optical.Magnesium>`, {py:class}`MagnesiumIon <sasktran2.optical.MagnesiumIon>` |
| Manganese | `Mn_I` | {py:class}`Manganese <sasktran2.optical.Manganese>` |
| Molybdenum | `Mo_I` | — |
| Sodium | `Na_I` | {py:class}`Sodium <sasktran2.optical.Sodium>` |
| Nickel | `Ni_I` | {py:class}`Nickel <sasktran2.optical.Nickel>` |
| Lead | `Pb_I` | — |
| Rubidium | `Rb_I` | {py:class}`Rubidium <sasktran2.optical.Rubidium>` |
| Scandium | `Sc_I`, `Sc_II` | — |
| Silicon | `Si_I` | {py:class}`Silicon <sasktran2.optical.Silicon>` |
| Tin | `Sn_I` | — |
| Strontium | `Sr_I`, `Sr_II` | {py:class}`Strontium <sasktran2.optical.Strontium>` (I) |
| Tantalum | `Ta_I` | — |
| Titanium | `Ti_I`, `Ti_II` | {py:class}`Titanium <sasktran2.optical.Titanium>` (I) |
| Vanadium | `V_I`, `V_II` | — |
| Tungsten | `W_I`, `W_II` | — |
| Yttrium | `Y_I`, `Y_II` | — |
| Zinc | `Zn_I` | {py:class}`Zinc <sasktran2.optical.Zinc>` |

## Molecular resonance properties

Molecular keys can also be passed to {py:class}`sasktran2.optical.MolecularResonance`.
Each molecular dataset represents one isotopologue and supports temperatures up
to **500 K**. Listed intervals are rounded inward from the file limits, in
vacuum nm. Coverage includes the source's selected electronic band systems;
it does not imply complete opacity throughout the interval.

| Species | Optical property | Vacuum interval [nm] | Source and version |
| --- | --- | --- | --- |
| ²⁷Al¹⁶O | {py:class}`AluminiumOxide <sasktran2.optical.AluminiumOxide>` / `MolecularResonance("AlO")` | 285.715–820 | [ExoMol ATP](https://exomol.com/data/molecules/AlO/27Al-16O/ATP/), 20210622; [line-list paper](https://doi.org/10.1093/mnras/stv507), [2021 update](https://doi.org/10.1093/mnras/stab2525) |
| ²⁴Mg¹⁶O | {py:class}`MagnesiumOxide <sasktran2.optical.MagnesiumOxide>` / `MolecularResonance("MgO")` | 270.271–820 | [ExoMol LiTY](https://exomol.com/data/molecules/MgO/24Mg-16O/LiTY/), 20241211; [line-list paper](https://doi.org/10.1093/mnras/stz912), [2024 update](https://doi.org/10.1093/rasti/rzae037) |
| ⁴⁰Ca¹⁶O | {py:class}`CalciumOxide <sasktran2.optical.CalciumOxide>` / `MolecularResonance("CaO")` | 400.001–820 | [ExoMol VBATHY](https://exomol.com/data/molecules/CaO/40Ca-16O/VBATHY/), 20230220; [line-list paper](https://doi.org/10.1093/mnras/stv2858) |
| ⁴⁸Ti¹⁶O | {py:class}`TitaniumOxide <sasktran2.optical.TitaniumOxide>` / `MolecularResonance("TiO")` | 333.334–820 | [ExoMol Toto](https://exomol.com/data/molecules/TiO/48Ti-16O/Toto/), 20240509; [line-list paper](https://doi.org/10.1093/mnras/stz1818), [2024 update](https://doi.org/10.1093/rasti/rzae037) |

## Emission templates

These classes are in `sasktran2.constituent`. They provide fixed photon spectra
with independently specified volume emission rates, rather than resonance
cross sections. The original wavelength coordinates are retained; their source
wavelength medium is unspecified.

| Species | Constituent | Interval [nm] | Source and qualification |
| --- | --- | --- | --- |
| FeO | {py:class}`FeOVolumeEmissionRate <sasktran2.constituent.FeOVolumeEmissionRate>` | 560–719.9 | Gattinger numerical model archived in the [Unterguggenberger et al. (2017) supplement](https://acp.copernicus.org/articles/17/4177/2017/acp-17-4177-2017-supplement.zip), `fig3.dat`, column 12; fixed excited-state populations. |
| NiO | {py:class}`NiOVolumeEmissionRate <sasktran2.constituent.NiOVolumeEmissionRate>` | 430–670 | Approximate 5 nm digitization of [Evans et al. (2011), Fig. 4](https://acp.copernicus.org/articles/11/9595/2011/); provisional model and attribution, questioned by [Noll et al. (2024)](https://acp.copernicus.org/articles/24/1143/2024/). |

The [model reference](metal_resonance.md) describes the resonance approximation
and data format. Dataset metadata retain the source URLs, versions, selection
limits, and attribution; the [database manifest](https://arg.usask.ca/sasktranfiles/sasktran2_db/v_latest/spectroscopy/metals/manifest.json)
lists the distributed files and checksums.
