# Molecular spectroscopy for ordinary mesospheric metal layers

Scope: ambient natural layers from the continuous meteoric input. Meteor trails,
bolides, spacecraft reentry and artificial releases are not evidence of ambient
OSIRIS detectability here. Presence of a molecular line list establishes only
spectroscopic feasibility, not abundance or a successful detection.

## Priority assessment

| Molecule | Ambient evidence | Appropriate model / present status |
| --- | --- | --- |
| AlO | WACCM-Al predicts an ambient layer peaking near 89 km; a lidar attempt gave an upper limit rather than a detection. | Strongest additional molecular solar-resonance candidate. ATP line/level data prepared. |
| FeO | Orange-band nightglow detected directly by OSIRIS. | Numerical Gattinger model archived by Unterguggenberger2017 supplied for free band VER (560–719.9 nm); no absolute absorption or excitation model is needed. |
| NiO | Proposed blue/red chemiluminescent contribution to OSIRIS and GLO-1 ambient nightglow; later work questions its significance. | Approximate 5 nm digitization of the isolated Evans2011 model supplied for exploratory free VER (430–670 nm), with provisional attribution. |
| MgO, CaO | Intermediates in ambient metal chemistry; no ambient OSIRIS molecule detection established in this research. | ExoMol LiTY/VBATHY spectroscopy prepared as exploratory candidates. Detection feasibility requires a chemistry/column and radiance calculation. |
| SiO | Modeled ambient silicon reservoir, with strong electronic bands mainly ultraviolet. | New SiOUVenIR list includes electronic transitions; original EBJT list is infrared only. A useful OSIRIS fluorescence system is not established here. Data not prepared. |
| TiO, VO | Visible spectra exist but no convincing ordinary mesospheric layer evidence established here. | Not justified as priority ambient-layer targets. TiO already downloaded/prepared as optional exploratory data; no large VO HyVO download performed. |
| OFeOH, FeO2 | Proposed excited products in ambient iron chemistry (Noll2024). | No usable species-specific emission spectral shape found. Keep as spectral-template gaps; missing absolute absorption is not a blocker for free VER. |
| AlOH, CaOH, MgOH, FeOH, NaOH, KOH, NaO, KO, LiO, metal carbonates/dioxides | Some are atmospheric reservoirs/intermediates. | No validated ordinary-layer visible emission template established. NaO and KO chemistry leads to atomic Na/K nightglow, which can already use independently retrieved atomic line VERs. |

The wavelength window is nominal, not a sensitivity curve. Evans et al. (2011)
excluded 480–530 nm around the OSIRIS order sorter in their nightglow analysis.
AlO's main 484.23 nm air bandhead lies in this interval. Other AlO bands and the
actual calibrated instrument response must be considered before predicting a
detection.

## Primary atmospheric references

- Plane et al. (2021), *Meteor-Ablated Aluminum in the Mesosphere-Lower Thermosphere*,
  https://doi.org/10.1029/2020JA028792 and author repository
  https://eprints.whiterose.ac.uk/id/eprint/170006/ . Modeled ambient AlO peaks
  near 89 km, night/day density about 10/60 cm^-3, and lidar upper limit 60 cm^-3.
  This is a candidate prediction, not an ambient AlO detection.
- Evans et al. (2010), *Discovery of the FeO orange bands in the terrestrial night
  airglow spectrum obtained with OSIRIS on the Odin spacecraft*,
  https://doi.org/10.1029/2010GL045310 . Observed emission is chemiluminescence.
- Evans et al. (2011), *The observation of chemiluminescent NiO* emissions in the
  laboratory and in the night airglow*, https://doi.org/10.5194/acp-11-9595-2011 .
  Emission extends longward of 440 nm; their molecular model is explicitly
  preliminary and uses nonthermal vibrational populations. It is not an absolute
  absorption cross-section table.
- Unterguggenberger et al. (2017), https://doi.org/10.5194/acp-17-4177-2017 .
  The numerical supplement `plot_data/fig3.dat`, column 12, contains the Gattinger
  FeO model. The template builder retains this original table and provenance.
- Noll et al. (2024), https://doi.org/10.5194/acp-24-1143-2024 . Robust FeO main
  feature, substantial uncertainty in NiO attribution, and proposed OFeOH emission.
  Free VER fits remain useful without adopting a chemical excitation model.
- Fjodorow et al. (2021), *Determination of gas-phase absorption cross-sections of
  FeO in a shock tube using intracavity absorption spectroscopy near 611 nm*,
  https://doi.org/10.1016/j.proci.2020.06.251 . Covers only 16316–16353 cm^-1;
  measured cross sections were obtained around 2200 K and 1.3 bar. The abstract
  reports individual oscillator strengths, but no complete cold-layer table was
  accessible. Do not apply those hot peak cross sections directly at 200 K.

## Data construction

`tools/spectroscopy/build_exomol_metals.py` downloads current source states,
transition lists, partition functions and definitions from ExoMol. The output
under `spectroscopy/metals/molecular/` includes single-isotopologue spectra for
AlO (27Al16O ATP), MgO (24Mg16O LiTY), CaO (40Ca16O VBATHY), and TiO (48Ti16O Toto).
The default spectral selection is 270–820 nm vacuum, limited further by the
source's maximum wavenumber. The operational lower-state-energy cutoff is
5000 cm^-1 and the default maximum temperature is 500 K. Full spectral subsets
are retained in `full/`, and original downloads in `raw/exomol/`.

All radiative upper/lower decay sums use the full original transition file,
before wavelength or population filtering. Source level weights include nuclear
spin; LTE populations must use `lower_statistical_weight`, not merely `2J+1`,
with the provided partition sum. No isotope-abundance factor is applied: number
density refers to the named isotopologue. Full state arrays remain available for
consistent partition sums. An omitted LTE population fraction is recorded; it is
not an opacity error bound.

The source is Einstein-A spectroscopy, not a measured phase database. An isolated
E1 angular momentum model may estimate the elastic return phase if it states its
assumptions. Only A_ul / sum_l A_ul returns to the original lower state; other
branches change wavelength. Molecular fluorescence, non-LTE pumping, overlapping
line interference, hyperfine structure, collisions and magnetic effects are not
made valid by assigning a Rayleigh phase or unit albedo.

The published models contain selected low electronic states. In particular, ATP
contains AlO X, A and B states; higher electronic systems are omitted. Coverage
inside a file's wavelength limits therefore does not establish complete molecular
UV-visible opacity. The source references define the supported band systems.

## Spectroscopic references and data license

- AlO ATP: https://doi.org/10.1093/mnras/stv507 ; empirical update
  https://doi.org/10.1093/mnras/stab2525
- MgO LiTY: https://doi.org/10.1093/mnras/stz912
- CaO VBATHY: https://doi.org/10.1093/mnras/stv2858
- TiO Toto: https://doi.org/10.1093/mnras/stz1818
- MgO/TiO updated state energies: https://doi.org/10.1093/rasti/rzae037
- VO current recommended dataset is hyperfine-resolved HyVO:
  https://doi.org/10.1093/mnras/stae542 ; the older VOMYT should not be silently
  presented as the current best dataset.
- ExoMol data license: https://exomol.com/data/licence/ specifies CC BY-SA 4.0.
  Preserve source attribution, modifications, source hashes and license when
  later hosting the derived data. Data licensing does not relicense project code.
