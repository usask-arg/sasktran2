(_users_metal_candidates)=
# Ordinary mesospheric metal candidates for OSIRIS

This survey concerns ambient mesospheric layers sustained by the normal cosmic-dust
input. Transient meteor trails, spacecraft plumes, and artificial releases are
outside its detection rationale. The nominal OSIRIS spectrograph interval is
[274–810 nm](https://research-groups.usask.ca/osiris/index.php). A line inside that
interval does not establish detectability: abundance, chemistry, solar irradiance
at the resolved pumping line, throughput, backgrounds, and measurement noise also
matter. No instrument noise calculation has been performed here.

## Atomic and ionic shortlist

The wavelengths below are **vacuum nm**, rounded from the locally retained
[NIST Atomic Spectra Database](https://physics.nist.gov/asd) records. They are
examples, not the complete line lists. Each file contains oscillator strengths,
angular momenta, available decay rates, lower-state energies, and source references.
Use the [resonance guide](metal_resonance.md) for cross sections and phase functions.

| Species / optical property | Example vacuum wavelengths (nm) | Assessment for ordinary layers |
| --- | --- | --- |
| Na I / `Sodium` | 589.1583, 589.7558 | Established OSIRIS resonance retrieval. |
| K I / `Potassium` | 766.7009, 770.1084 | Established OSIRIS retrieval using D1; O₂ complicates D2. |
| Li I / `Lithium` | 670.961, 670.976 | Observed ambient layer; much weaker abundance than Na. |
| Mg I / `Magnesium` | 285.2964 | Satellite metal-layer target; assess OSIRIS UV sensitivity. |
| Mg II / `MagnesiumIon` | 279.6352, 280.3531 | Satellite ion-layer target close to OSIRIS's UV edge. |
| Ca I / `Calcium` | 422.7918 | Observed ambient layer; available decay sum is incomplete. |
| Ca II / `CalciumIon` | 393.4777, 396.9591 | Observed ion layer; return fractions are about 0.93. |
| Fe I / `Iron` | 302.1519, 372.0993 | Established ambient metal layer with many pumping lines. |
| Ni I / `Nickel` | 337.0531, 341.5744, 352.5544 | Observed ambient layer; low-lying excited levels contribute. |
| Al I / `Aluminium` | 308.3046, 309.3606, 394.5122 | Atomic Al rapidly oxidizes; AlO is a more relevant reservoir. |
| Cr I / `Chromium` | 357.9708, 359.4511, 360.6357 | Trace meteoric screening candidate; no OSIRIS detection established here. |
| Mn I / `Manganese` | 279.5641, 279.9094, 280.1907 | Trace candidate with strong lines near the UV edge. |
| Ti I, Co I, Cu I | 363.6498; 352.7852; 324.8473, 327.4896 | Trace candidates requiring abundance and radiance assessment. |
| Zn I / `Zinc` | 307.6791 | Weak intercombination line; auxiliary natural trace candidate, without established OSIRIS detectability. |
| Rb I, Sr I, Ba I | 780.2415, 794.9789; 460.8624; 553.7018 | Speculative trace candidates, not demonstrated OSIRIS targets. |
| Fe II / `IronIon`, Si I / `Silicon` | 276.0150, 277.6159; 300.7615, 302.0884 | Retained cold-layer lines are weak; stronger ground-term resonances are outside the band. |

The direct OSIRIS evidence is from
[Gumbel et al. (2007), sodium](https://doi.org/10.1029/2006GL028687) and
[Dawkins et al. (2014), potassium](https://doi.org/10.1002/2014GL060801).
The latter also documents the O₂ interference and instrument/solar-background
constraints on potassium. Mg and Mg⁺ limb retrievals were demonstrated with
[SCIAMACHY](https://acp.copernicus.org/articles/8/1963/2008/); that is evidence for
the atmospheric target, not an OSIRIS detection claim.

For less familiar observed layers, see
[Gerding et al. (2019), nickel](https://doi.org/10.1029/2018GL080701) and
[Gerding et al. (2025), lithium](https://doi.org/10.1029/2025GL118710).
The ambient Li observations/model comparison are relevant here; enhanced plume
events are outside the present scope. The broader chemical and observational
context is reviewed by [Plane et al. (2015)](https://doi.org/10.1021/cr500501m).

## Metal-bearing molecules

| Species | Available property/data | Evidence and remaining limitation |
| --- | --- | --- |
| AlO | `AluminiumOxide`, ExoMol ATP | Natural-layer chemistry and a lidar upper limit support investigating it. The ordinary layer was not detected in that lidar study. |
| MgO | `MagnesiumOxide`, ExoMol LiTY | Spectroscopic coverage is available; ambient OSIRIS detectability is unestablished. |
| CaO | `CalciumOxide`, ExoMol VBATHY | Spectroscopic coverage is available; ambient OSIRIS detectability is unestablished. |
| TiO | `TitaniumOxide`, ExoMol TOTO | Auxiliary molecular spectroscopy; no ambient-layer detection claim. |
| FeO | `FeOVolumeEmissionRate`, numerical photon template | OSIRIS orange-band nightglow; free band VER over 560–719.9 nm. Absolute absorption and chemistry are unnecessary for this fit. |
| NiO | `NiOVolumeEmissionRate`, approximate photon template | Provisional identification; 5 nm figure digitization over 430–670 nm supports exploratory VER fitting. Later work questions a substantial NiO contribution. |
| OFeOH, FeO₂ | Candidate inventory | Energetically plausible chemical emission; no usable species-specific photon spectrum found. A spectral-template gap remains. |
| CoO, CaOH, MgOH | Candidate inventory | Laboratory optical signals or ambient reservoir chemistry motivate a search, but no validated ordinary-layer photon template was found. |
| VO, other hydroxides and larger metal-bearing species | Candidate inventory | No validated ordinary-layer emission template or OSIRIS resonance property supplied. |

The ambient aluminium analysis is
[Plane et al. (2021)](https://doi.org/10.1029/2020JA028792).
The strong AlO blue-band feature near 484 nm overlaps the order-sorter region:
the [OSIRIS NiO analysis](https://acp.copernicus.org/articles/11/9595/2011/)
excluded 480–530 nm. Other AlO bands and the actual instrument calibration must
therefore be considered. The ATP list also has a short-wavelength completeness
limit near 285.7 nm, so nominal OSIRIS coverage cannot be assumed for that dataset.

The molecular optical properties compute LTE absorption and same-transition elastic
return with line-specific angular-momentum phase functions. They **do not predict
total molecular fluorescence**, which is fed by absorption in other lines, or
chemical nightglow. Independently retrieved emission sources are supported for
[FeO's observed orange bands](https://doi.org/10.1029/2010GL045310) and
[the proposed NiO component](https://acp.copernicus.org/articles/11/9595/2011/).
See the [VER guide](metal_emission.md) for the templates, normalization and
uncertainties. A fitted VER requires neither an excitation rate nor an abundance
conversion. [Noll et al. (2024)](https://acp.copernicus.org/articles/24/1143/2024/)
support FeO's 595 nm feature, question a strong NiO contribution, and identify
OFeOH as another possible chemical emitter without an available spectrum.

## Atomic chemical nightglow

Na D and K D lines also have ordinary-layer chemical nightglow contributions.
Use independent `MonochromaticVolumeEmissionRate` constituents to retrieve their
line VERs, with the wavelengths and atom masses from the spectroscopy files.
The emission amplitude need not follow a chemical model or the resonance
cross section. In particular, do not fix the Na D ratio to 2:1 by default;
[Slanger et al. (2005)](https://doi.org/10.1029/2005JD006078) document its
variability. For potassium see
[Noll et al. (2019)](https://doi.org/10.1029/2018JD030044).
NaO and KO participate in the chemical pathways but this evidence does not
establish separate detected molecular emission bands.

## Installed inventory and data gaps

The initial atomic screen queried 37 elements in neutral and singly ionized form:
74 spectra, yielding 47 usable files and 1432 selected E1 transitions. This broad
screen is larger than the recommended atmospheric shortlist. The remaining
elements are auxiliary spectroscopy, with no assertion that an ordinary layer is
detectable. List the actual local files using:

```python
from sasktran2.database import MetalSpectroscopyDatabase

database = MetalSpectroscopyDatabase()
atomic_keys = database.available_species()
molecular_keys = database.available_species(kind="molecular")
emission_keys = database.available_species(kind="emission")
```

`AtomicResonance("Cr_I")` and `MolecularResonance("AlO")` select independent
properties. Additional charge states use the same interface; number density
always refers to that particular species, not to the total elemental inventory.

Missing oscillator strengths, unresolved classifications, and out-of-band
resonances are recorded in the local atomic `catalog.json`. Inconsistent NIST
Einstein-A/oscillator-strength pairs for Zn I and Ti II were reconciled against
primary transition-probability tables. Oscillator strengths are now derived from
the verified A values and statistical weights; the original values and explicit
reconciliation ledger are retained. One coarsely rounded W I lifetime branch
remains quarantined. Missing database records
are not physical upper limits on abundance or signal.

The repository's `tools/spectroscopy/README.md`, `MOLECULAR.md`, `EMISSION.md`,
and `screening_snapshot.json` preserve the survey, preparation commands, sources,
licenses, and exclusions. Runtime data are in
`<database_root>/spectroscopy/metals/{atomic,molecular,emission}/`. Default loaders
fetch missing files through `StandardDatabase`; an explicit database root remains
local-only. Source notices, attribution, preparation records, and SHA-256 checksums
accompany the hosted files. The molecular data retain CC BY-SA 4.0 and the emission
templates CC BY 3.0; no open-data license is assigned here to the NIST-derived data.
