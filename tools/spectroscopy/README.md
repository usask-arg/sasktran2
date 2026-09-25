# Atomic metal spectroscopy acquisition

This catalog supports a search for **ordinary natural mesospheric metal layers**
with Odin/OSIRIS. An atomic transition in the instrument band is a necessary
condition, not evidence that the species is detectable. No transient meteor,
spacecraft reentry, or artificial-release abundance is used to justify a candidate.

## Evidence and priorities

| Species | Evidence relevant to ordinary mesospheric layers |
| --- | --- |
| Na I, K I | Atomic layers retrieved from OSIRIS dayglow. |
| Li I, Ca I, Ca II, Fe I, Ni I | Natural layers observed with resonance lidar or other atmospheric measurements; OSIRIS sensitivity still needs calculation. |
| Mg I, Mg II | Natural layers retrieved from SCIAMACHY at approximately 285 and 280 nm; OSIRIS ultraviolet sensitivity is an additional constraint. |
| Al, Si, Cr, Mn, Ti, Co, Cu, Zn, V | Meteoritic-composition screening candidates. Atmospheric partitioning, low abundance, and usable solar-pumped lines need assessment. |
| Rb, Sr, Ba | Trace meteoric candidates; ordinary-layer OSIRIS detectability is unestablished. |
| Other downloaded elements | Auxiliary atomic-data screen only; their presence in the cache is not an ambient-detectability claim. |

The established OSIRIS results are [Na (Gumbel et al., 2007)](https://doi.org/10.1029/2006GL028687)
and [K (Dawkins et al., 2014)](https://doi.org/10.1002/2014GL060801).
[Scharringhausen et al. (2008)](https://acp.copernicus.org/articles/8/1963/2008/)
demonstrate Mg/Mg+ limb retrievals with SCIAMACHY.
[Collins et al. (2015)](https://doi.org/10.1002/2014GL062716) report Ni lidar
observations and summarize earlier metal-layer measurements.
[Plane et al. (2015)](https://doi.org/10.1021/cr500501m) review the chemistry,
abundances, and observational history of natural mesospheric metals.

FeO has an OSIRIS **chemiluminescence** detection
([Evans et al. 2010](https://doi.org/10.1029/2010GL045310)); the NiO component proposed
by [Evans et al. (2011)](https://acp.copernicus.org/articles/11/9595/2011/)
remains provisional in light of later observations. These emission analyses do
not establish atomic Fe/Ni resonance signals or molecular resonance cross sections.
See [the molecular source assessment](MOLECULAR.md) for AlO and other ambient-layer
candidates and ExoMol data preparation, and [the emission assessment](EMISSION.md)
for the supplied FeO/NiO free-VER templates and remaining spectral-shape gaps.

## Summary plots

```sh
pixi run python tools/spectroscopy/plot_metal_summary.py
```

This reads the installed data and writes a seven-figure gallery under
`artifacts/metal_spectroscopy_summary/`, with PNG and vector SVG versions,
numerical NumPy arrays, methods and source hashes. The inventory includes all
47 atomic/ionic files and all four molecular line lists. Separate figures show
absolute cross sections, elastic branching, phase functions, resolved line
profiles and the FeO/NiO photon templates.

The default cross-section overview uses 200 K and a 1 nm FWHM Gaussian display
kernel. `--temperature` and `--display-fwhm-nm` change these choices. This is
illustrative spectral broadening of line areas, not an instrument response or
radiance simulation. Resolved profiles use the runtime Voigt calculation at
150/200/300 K without this blur. Emission templates receive no added smoothing.
The inventory normalizes each species independently; absolute plots use m² per
particle with no abundance scaling. `--output` selects another plot directory.

## Reproduction and stored data

```sh
pixi run python tools/spectroscopy/build_atomic_metals.py \
  --output "<database_root>/spectroscopy/metals/atomic"
```

`--species Na_I Mg_II` limits acquisition; `--refresh` replaces snapshots.
Existing valid snapshots are reused. Downloads and derived NetCDF files are
written atomically. The script uses the public NIST forms with curl and requires
the project's numpy, xarray, and netCDF4 environment.

Each species has an independent NetCDF file. `catalog.json` records every screened
spectrum, including unsuccessful screens, rejected records, request URLs, hashes,
and dates. Original tables remain in `raw/nist`. Principal metadata and assumptions
are embedded in each NetCDF file; NIST reference identifiers remain attached to
individual transitions. `screening_snapshot.json` in this directory is the checked
support inventory from the initial acquisition, not a replacement for the local data.

The source is [NIST ASD version 5.12](https://www.nist.gov/pml/atomic-spectra-database),
Kramida, Ralchenko, Reader, and NIST ASD Team (2024),
[DOI 10.18434/T4W30F](https://doi.org/10.18434/T4W30F), accessed 2026-09-24.
ASD's official page carries a U.S. Department of Commerce copyright notice; these
files are not labeled public domain or assigned an open-data license here.
The [NIST SRD copyright statement](https://www.nist.gov/open/license) applies
separately from the project code license. Source notices, citations, and the
per-line bibliographic identifiers are retained; this preparation does not grant
an additional license to the underlying NIST data.

The query requests vacuum wavelengths and all listed decay channels at all
wavelengths, with upper energy at most 50000 cm−1. Optical records are then selected
for 274–810 nm, lower energy at most 5000 cm−1, electric-dipole type, positive A and
oscillator strength, and resolved electronic angular momenta. Stable NIST level IDs
identify branches and partition states. See [ASD field definitions](https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html).
Masses represent terrestrial isotope mixtures using [CIAAW values](https://ciaaw.org/atomic-weights.htm).

No usable record in this screen can mean out-of-band resonance, missing
classification, missing transition strength, or a rejected inconsistent record.
It does not prove absence of emission at other excitation conditions. In particular,
strong Fe II and Si I ground-term lines lie below the selected band, and their
remaining in-band cold-layer lines are weak. Be I's strong resonance is also outside
the band. Zr and Nb have no classified transitions returned by this bounded query.

## Physics and limitations

For an isolated electric-dipole transition and an unpolarized lower state, the
electronic polarizability is

\[
W_2=3(2J_u+1)
\begin{Bmatrix}1&1&2\\J_u&J_u&J_l\end{Bmatrix}^{2}.
\]

This gives 1 for J=0→1, 0 for an alkali D1 line, and 1/2 for D2. With the usual
phase-function normalization, P11(μ)=1+(W2/2)P2(μ). The general redistribution
theory and limits of the isolated-line treatment are given by
[Belluzzi and Trujillo Bueno (2014), equations 1–3](https://arxiv.org/abs/1403.1701).
Hyperfine and isotope structure, magnetic fields, lower-state alignment, collisions,
and interference between electronic levels are not represented by this electronic
W2. These omissions matter for precision resonance polarimetry, particularly alkalis.

A normalized frequency profile has integrated cross section π r_e c f. Thermal
Gaussian standard deviation is (ν0/c)√(kT/m); natural Lorentz half-width is
(Γu+Γl)/(4π). Multiply by the lower-level population fraction. An LTE partition
formed only from classified levels below 5000 cm−1 is a cold-atmosphere
approximation, not a general high-temperature model.

`upper_total_a_s` and `lower_total_a_s` are sums of distinct **available** radiative
decays, including forbidden transitions. They are lower bounds when rates are
missing. `upper_decay_data_complete` is set only for principal E1-closed doublets
and the Mg I singlet resonance; negligible forbidden channels are not certified.
For other lines, Aul divided by the available sum is an upper-bound estimate of
the same-line return probability, not proof of conservative scattering. Ca I is
an important example: a tabulated ratio of one does not establish a closed upper
level. Raman decay into a different lower level produces a different wavelength.

NIST's f values are derived quantities, not independent measurements. An audit
found conversion discrepancies for Ti II and Zn I. Corroboration against
[NIST's zinc persistent lines](https://physics.nist.gov/PhysRefData/Handbook/Tables/zinctable3_a.htm),
[NIST's titanium table](https://physics.nist.gov/PhysRefData/Handbook/Tables/titaniumtable4.htm)
and [Wood et al. (2013), Table 4](https://cdsarc.cds.unistra.fr/ftp/J/ApJS/208/27/table4.dat)
supports the A values at their stated scale. `atomic_strength_reconciliations.json`
records the nine reconciled decay records; these restore Zn I 307.679063 nm and
Ti II 335.036506/344.529311 nm to the operational set. Oscillator strengths now
derive consistently from A, wavelength and resolved degeneracies. Raw source f
values, discrepancies and per-line reconciliation flags remain in the files.
Corroborating source snapshots and hashes live under `atomic/raw/reconciliation/`.
One unresolved, coarsely rounded W I decay record remains quarantined, affecting
a lifetime sum rather than an in-band absorption line. The final atomic inventory
contains 47 species files and 1432 lines from 74 screened spectra.

For chemical emission, `build_metal_emission_templates.py --output
"<database_root>/spectroscopy/metals/emission"` prepares a numerical FeO template
and an explicitly approximate NiO template, with original source files and trace
provenance. These support free band VERs without an abundance/excitation model;
see `EMISSION.md` and the source metadata for limitations.

Line centers must be resolved in the radiative-transfer calculation before
instrument convolution. Using cross sections already smoothed to OSIRIS's spectral
resolution can lose saturation and solar Fraunhofer-line effects. Instrument band
coverage alone cannot substitute for a radiance and uncertainty calculation.

## Standard database distribution

The prepared runtime files use the standard database prefix
`spectroscopy/metals/{atomic,molecular,emission}/`. The default
`MetalSpectroscopyDatabase` downloads missing files through `StandardDatabase`.
An explicit `db_root` selects local-only access. Availability lists describe the
local cache, not the full remote inventory.

On the SasktranFiles share, the distribution is under
`sasktran2_db/v_latest/spectroscopy/metals/`. Alongside the 53 runtime NetCDF files
are source records, attribution and license notices, and a `manifest.json` with
SHA-256 hashes. The large ExoMol `raw/` and `full/` acquisition intermediates are
not needed at runtime and are excluded; their URLs and hashes remain in the
molecular provenance records. Generated summary figures stay in the ignored
`artifacts/metal_spectroscopy_summary/` directory and can be regenerated with
`plot_metal_summary.py`.
