# sasktran2-nlte: Non-LTE Emission Plan

This note plans a non-LTE excited-state population capability for SASKTRAN2, aimed at forward simulation of airglow and dayglow in the UV, visible and near-IR. The review of the GRANADA IR non-LTE code that motivated this work is in `nlte_granada_review.md`. That review also outlines a possible later IR extension.

## Decisions

| Topic | Decision |
|---|---|
| Scope | Forward simulation only. No retrievals or Jacobians in this plan. |
| Priority emissions | Oxygen first (O2 and atomic O), then OH and NO. H2O enters as an input (HOx source, quencher, absorber), not as an emitter. IR vibrational non-LTE is deferred. |
| Isolation | The kinetics are a black box in a new crate, `rust/sasktran2-nlte`, with a small, well-defined API. It has no dependency on `sasktran2-rs`, `sasktran2-core` or `sasktran2-sys`. |
| Photolysis | A TUV-like actinic flux and photolysis interface, `sasktran2.photolysis`, built on the SASKTRAN2 engine (discrete ordinates). It supplies photolysis and excitation rates to the crate. |
| Existing code | The current `photchem` module migrates into the new crate, with compatibility shims during the transition. |
| Background chemistry | Prescribed from climatologies (CAIRT ERS, MSIS, WACCM-derived). The crate solves only excited-state populations. |
| Provenance | Clean-room implementation. Algorithms come from the literature. Every rate constant, yield and A-value carries a citation in data files. No translation of GRANADA/KOPRA (LGPL-2.1) or TUV 5.4 Fortran. |

## Architecture

```text
profiles: T, p, background number densities (O, O2, N2, O3, H, OH, HO2, NO, N, H2O, ...)
   |                                   |
   v                                   |
sasktran2.photolysis  (Python on the SASKTRAN2 engine)
   actinic flux F(lambda, z)  ->  photolysis rates J_k(z), line excitation rates g_l(z)
   |                                   |
   |  rates by name                    |
   v                                   v
sasktran2-nlte  (Rust crate, black box)
   mechanism file (TOML, cited) + column inputs
   ->  excited-state densities, process budgets, transition and line emission rates
   |
   |  xarray Dataset (transition / line photon VER)
   v
emission glue  (existing sasktran2 constituents)
   ->  emission_source on the spectral grid
   v
sasktran2 Engine: limb and nadir radiances, self-absorption, resonance scattering
```

- **Contracts.** At the Python level, the boundaries are xarray Dataset schemas. At the Rust level they are plain structs. Nothing else crosses a boundary.
- **Crate dependencies.** `sasktran2-nlte` is a leaf crate (ndarray, nalgebra, serde, toml, thiserror). `sasktran2-rs` depends on it only for emission glue, and `sasktran2-py-ext` exposes bindings. The crate never imports SASKTRAN2 types.

## Piece 1: `sasktran2.photolysis`

This is a TUV-equivalent actinic flux and photolysis module built on the SASKTRAN2 engine. It grows out of today's `sasktran2.photchem.actinic_flux()`.

```python
calc = sk.photolysis.ActinicFlux(altitudes_m, grid="airglow", streams=4)
flux = calc.calculate(state, sun=sk.photolysis.Sun(cos_sza=0.3), albedo=0.3)
# actinic_flux(wavelength, altitude), direct and diffuse parts, photons m-2 s-1 nm-1

rates = sk.photolysis.PhotolysisRates(["O3 -> O(1D) + O2(a)", "O2 -> O(3P) + O(1D)"])
j = rates.calculate(flux, state)                 # j_<id>(altitude) in s-1
g = sk.photolysis.LineExcitation(lines).calculate(state, sun)   # g_<id>(altitude) in s-1
```

### Spectral grids

- Presets: a base grid as today (120–1280 nm at 0.1 nm), plus high-resolution windows where line structure matters:
  - O2 A, B and 1.27 µm bands (already present);
  - OH A–X, 306–320 nm;
  - NO γ, 200–250 nm.

### Radiative transfer

- Discrete ordinates in pseudo-spherical geometry with flux observers, as used today.
- Optionally the two-stream solver, once it can output fluxes.

### Cross sections and quantum yields

- A reaction library with temperature (and, where needed, pressure) dependence.
- Data sources:
  - JPL Evaluation 19 (JPL Publication 19-5);
  - the TUV-x data set (Apache-2.0);
  - SASKTRAN2's existing O3 DBM, O2 Schumann–Runge and Lyman-α optics, so photolysis and the forward model stay consistent.
- First reactions:
  - O2 (Schumann–Runge continuum, Lyman-α, Herzberg continuum);
  - O3 channels (Hartley, Huggins, Chappuis);
  - H2O (Lyman-α, 175–190 nm);
  - NO2;
  - NO predissociation in the δ bands.

### Line excitation ("g-factors")

- For each line, the line cross section is integrated against the attenuated solar flux.
- This uses a hybrid scheme:
  - the direct beam is computed line by line on each line's own profile, which is cheap with Beer–Lambert along the solar path;
  - the diffuse field comes from coarse-resolution DO.

### Sun

- Solar zenith angle from time and location.
- Earth–Sun distance.
- TSIS-1 HSRS v2 (already used by `SolarIrradiance`).
- Line-integrated Lyman-α flux.

### Replacing the TOA scaling

The current photo-reactions scale photolysis rates to fixed top-of-atmosphere values (`toa_rate_constant`). These are replaced by computed J values. The scaling remains available as an option for parity tests.

### Validation

- TUV-x, run via the MUSICA Python package or a local build.
- TUV 5.4, used only as a reference executable, from the KIT repository on the T9 drive.
- ERS `npar` profiles, which are TUV photolysis rates used as GRANADA inputs.

### Engine changes

Each would be a separate small PR, made only if needed:
- flux output from the two-stream solver;
- lifting the forced wavelength block size of 1 when flux observers are present;
- confirming DO pseudo-spherical direct-beam attenuation at twilight (SZA > 90°), or adding a spherical direct beam from the SASKTRAN2 ray tracer.

## Piece 2: the `sasktran2-nlte` crate

### Responsibilities

Given a mechanism and column inputs, the crate computes:
- excited-state number densities;
- production and loss budgets per process;
- photon volume emission rates per radiative transition, and optionally per line.

The crate does no radiative transfer, does not read atmospheric databases, and does not use SASKTRAN2 types.

### Mechanism files

Mechanisms are TOML files bundled in `rust/sasktran2-nlte/mechanisms/` and versioned with the crate.
- Every rate, yield and A-value has a reference key.
- Units are explicit and converted to SI on load.

The numbers below are illustrative only:

```toml
[meta]
name = "oxygen"
version = "0.1.0"

[references.jpl19]
citation = "Burkholder et al. (2019), JPL Publication 19-5"

[[state]]
id = "O2(b,v=0)"
species = "O2"
energy_cm = 13120.0
degeneracy = 1

[[reaction]]
id = "o2b_quench_n2"
reactants = ["O2(b,v=0)", "N2"]
products = [["O2(a,v=0)", 1.0], ["N2", 1.0]]
rate = { law = "arrhenius", a = 1.0e-15, ea_over_r = 0.0, units = "cm3 s-1" }
reference = "jpl19"

[[photo_reaction]]
id = "o3_hartley_o1d"
reactant = "O3"
products = [["O(1D)", 1.0], ["O2(a,v=0)", 1.0]]
rate_input = "j_o3_o1d"          # supplied by sasktran2.photolysis

[[transition]]
id = "o2_a_band"
upper = "O2(b,v=0)"
lower = "O2(X,v=0)"
einstein_a = { value = 0.08, units = "s-1" }
reference = "hitran"
```

**Rate laws:**
- constant;
- Arrhenius with a power-law temperature term;
- JPL termolecular fall-off;
- JPL chemical activation;
- Troe;
- tabulated k(T);
- Landau–Teller;
- reverse rates from detailed balance, using state energies and degeneracies.

Parameterisations that cannot be expressed as data, such as the McDade/Barth green-line form, are named plugins implemented in Rust. Each plugin is documented and cited.

**Nascent distributions** are product yield lists, for example H + O3 → OH(v=6..9).

**Validation on load:**
- unknown species or states;
- unit errors;
- reactions that do not balance;
- a `required_inputs()` listing that names every background density and rate the mechanism needs.

### Column inputs

Each column provides:
- altitude, T and p (giving [M]);
- background densities by species id;
- photolysis and excitation rates by name.

Missing inputs produce an error that names them.

### Solver

- **Steady state, per altitude.** The unknowns are the excited states.
  - With no excited–excited products the system is linear and needs one small dense LU solve.
  - Bilinear terms are handled with Newton iteration, starting from the linear solution.
  - The LU is pure Rust (nalgebra), replacing the current LAPACK FFI call.
- **Transient integration** for slow states such as O2(a¹Δ) (lifetime about 72 minutes) at twilight. It uses an implicit integrator (TR-BDF2 or Rosenbrock), with photolysis rates supplied as a time series.
- **Diagnostics.** Residual norms and positivity checks are returned with the solution.
- **Parallelism.** Altitudes and columns run in parallel with rayon. There is no global state.

### Outputs

Rust:

```rust
pub struct Solution {
    state_density,           // [state, altitude]
    transition_photon_ver,   // [transition, altitude]
    budget,                  // [process, altitude]
    residual,                // [altitude]
}
```

Python returns the same content as an xarray Dataset with matching names, coordinates and units attributes.

### Line-resolved emission

The crate accepts generic line records:

```rust
LineRecord { upper_state, lower_state, wavenumber_cm, einstein_a, g_upper, e_upper_cm, labels }
```

It also takes a rotational model for each line:
- thermal at T;
- two-temperature (the OH high-N tail);
- explicit per-line upper-state populations (for example, driven by fluorescence excitation rates).

It returns the line photon VER per line and altitude. The line loaders (HITRAN through `OpticalLineDB`, ExoMol, PGOPHER/MoLLIST) live outside the crate and produce `LineRecord`s.

### Python API

```python
mech = sk.nlte.Mechanism.bundled("oxygen")      # or Mechanism.from_file(path)
mech.required_inputs()
sol = sk.nlte.solve(mech, state, rates)         # xr.Dataset
atmosphere["o2_a"] = sk.nlte.emission_constituent(sol, "o2_a_band", lines=o2_lines)
```

## Piece 3: emission glue

**Existing constituents to reuse:**
- `PopulationEmissionRate`;
- `LineListVolumeEmissionRate`, generalised to take the molecular mass per species (it currently hard-codes O2);
- `SpectralVolumeEmissionRate` templates;
- `MonochromaticVolumeEmissionRate`;
- the metal resonance optics.

**New pieces:**
- `LineRecord` adapters for HITRAN, ExoMol NO and the OH A–X line list;
- band templates for the O2 Herzberg/Chamberlain systems and the NO2 continuum, where no line list exists.

**Self-absorption and scattering:**
- Self-absorption along the line of sight is handled by the forward model's line absorbers.
- Resonance scattering of O2 A-band, OH A–X and NO γ emission can reuse the resonance optics from PR #302. Whether it is needed will be checked per band.

## Migrating `photchem`

| Current | Destination |
|---|---|
| `rust/sasktran2-rs/src/photchem/types.rs` (`Molecule`, reaction parsers) | `sasktran2-nlte` (state ids, reaction-string parsing) |
| `photchem/models.rs` (`PhotochemicalModel`, steady-state solve) | `sasktran2-nlte` solver (nalgebra LU instead of `bindings::lapack::dgesv`) |
| `photchem/models.rs` (`Yankovsky`) | `mechanisms/oxygen_yankovsky.toml` |
| `photchem/models.rs` (`calculate_photolysis_rate`, `wavelength_bin_widths`) | `sasktran2-nlte` numerics, called from `sasktran2.photolysis` |
| `photchem/emission.rs` (transitions, McDade green line, populations → VER, branching weights) | `sasktran2-nlte` |
| `photchem/emission.rs` (HITRAN O2 band builders using `OpticalLineDB`) | stay in `sasktran2-rs` as `OpticalLine → LineRecord` adapters |
| `emission/o2.rs`, `constituent/types/population_emission_rate.rs` | stay in `sasktran2-rs` and import from the crate |
| `sasktran2-py-ext/src/photchem/yankovsky.rs` (`PyYankovsky`) | generic `PyMechanism` / `PySolution`, with `PyYankovsky` as a temporary shim |
| Python `sasktran2.photchem` (`actinic_flux`, `Yankovsky`) | `sasktran2.photolysis.ActinicFlux` and `sasktran2.nlte`; `photchem` re-exports with a `DeprecationWarning` |

The migration happens in four steps:
1. Capture golden outputs from the current code: `Yankovsky.solve`, `emissions`, the McDade green line, and `PopulationEmissionRate` radiances.
2. Move the code mechanically with no behaviour change. Goldens must match to round-off.
3. Re-express Yankovsky as TOML and check parity against the goldens. The current Rust code has no literature citations and uses TOA-scaled photolysis rates, so recovering the provenance of each constant is part of this step.
4. Deprecate the old entry points.

## Emission targets and data sources

The kinetics sources are a starting survey. Exact citations are fixed per entry in the mechanism files.

| Emission | Excitation | Kinetics sources | Spectroscopy | Validation |
|---|---|---|---|---|
| O2(b) A and B bands (762, 689 nm) | Day: resonant solar absorption, O(¹D)+O2 energy transfer, O3/O2 photolysis chain. Night: O+O+M (Barth two-step) | Yankovsky & Manuilova (2006) and updates, JPL 19-5, McDade et al. (1986) | HITRAN O2 (in use) | ERS GRANADA O2(b), OSIRIS A-band (Sheese et al. 2010), existing tests |
| O2(a) 1.27 µm | O3 Hartley photolysis, b→a quenching, resonant absorption. Transient at twilight | As above | HITRAN O2 | ERS GRANADA O2(a), OSIRIS IRI, SABER 1.27 µm |
| O(¹S) 557.7 nm | Night: Barth (McDade et al. 1986, existing). Day: O2 photodissociation, O2⁺ recombination, N2(A) transfer | McDade et al. (1986), Slanger & Copeland (2003) review, JPL 19-5 | NIST ASD | Existing tests |
| O(¹D) 630/636.4 nm | Thermospheric: O2⁺+e recombination, O2 photolysis, N(²D)+O2. Needs ionospheric inputs | Literature survey | NIST ASD | Lower priority |
| O2 Herzberg I/II and Chamberlain (UV–blue nightglow) | O+O+M yields to A, A′ and c states; quenching by O, O2, N2 | McDade et al. (1986), Stegman & Murtagh (1991), Slanger & Copeland (2003) | No HITRAN coverage; band templates first | Literature spectra |
| OH Meinel (visible/near-IR) | H+O3 → OH(v≤9) nascent; HO2+O minor; quenching by O2 (multi-quantum), N2, O | Adler-Golden (1997), Sharma et al. (2015), Kalogerakis et al. (2016), Panka et al. (2017), JPL 19-5 | HITRAN OH (Brooke et al. 2016); thermal plus non-thermal rotational | ERS GRANADA OH(v=0–10), SABER 1.6/2.0 µm, OSIRIS |
| OH A–X (308 nm) dayglow | Solar resonance fluorescence; quenching at lower altitudes | Literature g-factors | Yousefi et al. (2018) A–X line list (MoLLIST), or LIFBASE | MAHRSI and SHIMMER literature |
| NO γ/δ dayglow (UV) | Solar resonance fluorescence; self-absorption in the (0,0) and (1,0) bands | Stevens (1995) and related g-factor work | ExoMol XABC (Qu et al. 2021) | Literature g-factors, SNOE |
| NO δ/β nightglow | N(⁴S)+O → NO(C, B) recombination | Rate sourcing to be done | ExoMol XABC | Literature |
| NO2 continuum (visible/near-IR) | NO+O(+M) → NO2* chemiluminescence | Laboratory rate and spectrum; OSIRIS (Gattinger et al. 2009) | Template spectrum | OSIRIS |
| H2O (input only) | HOx source (J(H2O), O(¹D)+H2O), quencher, absorber | JPL 19-5 | HITRAN | — |

**Background profiles:**
- CAIRT ERS: VMRs including O, O(¹D), H, OH, HO2, N(⁴S), N(²D), O3 and H2O;
- MSIS for the thermosphere;
- an empirical NO model.

### CAIRT ERS as a reference

The ERS archive (Zenodo 8256025, staged at `/Volumes/T9/data/cairt_ers_kopra`) contains complete GRANADA runs. Coverage:
- 0–200 km at 1 km;
- four months × five latitudes;
- day and night;
- solar-activity and volcanic variants;
- kinetic-perturbation runs.

| Content | Use in this plan |
|---|---|
| O2 X(v=0–35), a(v=0–5), b(v=0–2) population ratios | The same state set as the current Yankovsky model. A direct reference for the oxygen mechanism (phases 2 and 5). |
| OH X(v=0–10) population ratios, per spin component | Reference for the OH Meinel mechanism (phase 6). |
| Photolysis rates by product channel: O3 (→ O2(a, v=0–5) and → O2(X, v=0–35)), O2 (including the O(¹D) channel), NO2; plus O(¹D) density | Reference for `sasktran2.photolysis` (phase 3). These rates can also drive the crate directly, so kinetics can be validated independently of the photolysis module. |
| Kinetic perturbation runs (O2-O-vt, O+O2+M, O3-Oeq, NO2+hv, ...) | Sensitivity checks against GRANADA's response. |

GRANADA's $r$ is relative to LTE populations normalised over only the modelled states. Convert it before comparing (see `nlte_granada_review.md`).

## Data management

**Development data** lives on the external drive, under `/Volumes/T9/data/<dataset>/`. Each dataset has a README with its source URL or DOI, licence, checksum and download date.

**Shipped data:**
- Mechanism TOMLs are bundled in the crate.
- Tests use tiny extracts.
- Large tables use SASKTRAN2's existing database download mechanism with checksums.

| Dataset | Licence and use |
|---|---|
| CAIRT ERS (Zenodo 8256025, staged at `/Volumes/T9/data/cairt_ers_kopra`, MD5 verified) | CC-BY-4.0. Validation and background profiles. |
| TUV-x cross sections and quantum yields | Apache-2.0. Attribution and NOTICE required if redistributed. |
| JPL 19-5 | NASA publication. Values transcribed with citation. |
| HITRAN | Existing HAPI download path. |
| ExoMol, MoLLIST, LIFBASE | Terms to confirm before redistribution. |
| GRANADA/KOPRA, TUV 5.4 Fortran | Reference only. No code reuse. |

## Status

**Phases 0 and 1 are complete.**

- **Crate.** `rust/sasktran2-nlte` holds the former `photchem` reaction types, the Yankovsky model, photolysis-rate integration and the emission code. It depends only on anyhow, ndarray and nalgebra.
- **Linear solver.** The LAPACK `dgesv` call was replaced by a nalgebra LU solve.
- **What stays in `sasktran2-rs`.** `sasktran2_rs::photchem` re-exports the crate. Its `emission` module keeps the HITRAN adapters, now free functions: `oxygen_a_band_from_hitran`, `oxygen_b_band_from_hitran` and `emission_band_from_hitran_lines`.
- **No API change.** The Python API and `sasktran2-py-ext` are unchanged.
- **Regression fixtures.**
  - `tests/photchem/test_photchem_goldens.py` and `tests/photchem/goldens/yankovsky.npz` capture the pre-move solver, emission and green-line outputs on a deterministic synthetic actinic-flux input.
  - After the move, the largest relative deviation is 3e-12, from the LU change. The test tolerance is 1e-9.
- **ERS reader.** `tools/nlte/kopra_prf.py` reads ERS `.prf` files: p/T, VMR, `npar`, ratio, and the Mixer variants. It handles Fortran three-digit exponents. Its tests are in `tests/nlte/test_kopra_prf.py`.

**Phase 2: the steady-state mechanism engine is in place.**

- **Format.** Mechanisms are TOML files, documented in `nlte_mechanism_format.md`. The Rust `toml` crate adds five small pure-Rust packages.
- **Validation on load:**
  - species declarations;
  - citation keys;
  - unit/order consistency;
  - yield sums;
  - unknown keys;
  - `for_v` families with `{v±k}` templating.
- **Solver.** One linear solve per level, returning state densities, per-process rates (photon VER for transitions) and production/loss budgets.
- **Python.** `sasktran2.nlte` provides `Mechanism`, `solve` and `budget`.
- **Yankovsky.** It is now the bundled `oxygen_yankovsky` mechanism, translated entry by entry. Its rate inputs follow the ERS names (`J_O3_A{v}`, `J_O3_X{v}`), so ERS photolysis rates can be fed in directly.
- **Parity.**
  - A one-off comparison against the hand-coded solver used 200 random levels spanning 9 decades in density and rates. The old solutions satisfy the new equations to 5e-16 relative residual.
  - The regression fixtures still pass, with a worst difference of 3e-12.
- **Removed code.** The hand-coded reaction list, `PhotochemicalModel`, `ChemicalReaction` and `MoleculeMap` are gone.
- **What remains of the old Yankovsky code.** It now only supplies the legacy top-of-atmosphere-scaled photolysis rates (`Yankovsky.photolysis_rates`) until phase 3. `Yankovsky.solve` returns the same per-state Dataset as before, and `Yankovsky.solve_full` returns the full `nlte.solve` result.
- **Deferred:**
  - Time-dependent solves (O2(a) at twilight).
  - Processes with two state reactants. Loading rejects them.
  - `photchem` deprecation warnings, which wait until `sasktran2.photolysis` replaces the legacy rates in phase 3.

**Phase 3, first step: `sasktran2.photolysis`.**

- **`ActinicFlux`.** Discrete ordinates in pseudo-spherical geometry, with flux observers at every altitude. Inputs: Rayleigh scattering, the absorbers present in the atmosphere, a Lambertian surface, and the Earth–Sun distance. Output: actinic flux, top-of-atmosphere flux and the cross sections used.
- **Default grid.** The default `airglow_wavelength_grid` (120–1280 nm at 0.1 nm, 0.001 nm in the O2 A, B, γ and 1.27 µm windows, plus Lyman-α) takes about 7 s for 130 levels.
- **Rate definitions.**
  - `Photolysis`: the integral of flux × σ × φ(λ, T), with an optional check on grid resolution.
  - `LinePhotolysis`: one unresolved solar line, with a constant effective cross section.
  - `LymanAlphaPhotolysis`: O2 at Lyman-α, from the slant O2 column (Chabrillat & Kockarts 1997).
  - `photolysis_rates` returns an xarray Dataset of named rates.
- **Shared atmosphere.** One atmosphere Dataset (species ids, SI units) drives both `photolysis` and `nlte.solve`.
- **Quantum yield.** The O(¹D) yield from O3 follows Matsumi et al. (2002), with the coefficients checked against the TUV 5.4 implementation.
- **Presets.** `presets.oxygen_photolysis` gives physical channel rates. `presets.oxygen_yankovsky_rates` gives the 47 rate inputs of `oxygen_yankovsky`, using the legacy O2(a, v) and O2(X, v) product splits.
- **Validation against ERS** (TUV 5.4; `tools/nlte/validate_photolysis_ers.py`, april+00, albedo 0.2):

  | Region | O3 → O(¹D) | O3 → O(³P) | O2 photolysis, total |
  |---|---|---|---|
  | Above 60 km | +4% | within 1% | too high at 70–100 km (up to 1.57×), 7% low at 120 km |
  | Below 60 km | +22–26% at 30–40 km | within 7% | 2.5–4× too low (20–50 km) |

  - The stratospheric O(¹D) difference cannot be removed by changing the solar zenith angle, so it is not just a geometry assumption.
  - Total O2 photolysis is too low below 60 km because the Herzberg continuum is missing.
  - It is too high at 70–100 km because the 0.1 nm Schumann–Runge bands under-attenuate.
- **Computed versus legacy rates.** These change the oxygen populations by a few percent at 50–70 km. Up to ~2× less O(¹D) at 90–100 km, mainly because the legacy O(¹D) rate also counted Schumann–Runge band absorption.
- **The `Yankovsky` shim still uses the legacy rates.**

**Phase 3, O2 cross sections.**

- **New optical property.** `sk.optical.O2UV` covers 130–242.4 nm, from 130 to 500 K:

  | Region | Source |
  |---|---|
  | Below 175.44 nm | CfA Schumann–Runge continuum |
  | 175.44–204.08 nm | Schumann–Runge bands, resolved at 0.5 cm⁻¹ from the Minschwaner et al. (1992) polynomials |
  | 205–242.4 nm | Yoshino et al. (1988) Herzberg continuum, also extended under the bands to 194 nm |

  It replaces `O2SchumannRunge` in the photolysis defaults. Resolving the bands directly suits sasktran2's discrete-ordinates calculation better than an effective-cross-section parameterisation such as Koppers & Murtagh, which applies only to the direct beam.
- **Build and hosting.** `tools/spectroscopy/build_o2_uv.py` builds the table from the original sources, recording their SHA-256 hashes. The table belongs at `cross_sections/o2/o2_uv.nc` in the sasktran2 standard database; until it is uploaded, `O2UV` raises `OSError` and its test skips.
- **Grid.** `airglow_wavelength_grid` adds 0.002 nm sampling over the bands, giving 107k points in total.
- **Effect on total O2 photolysis against ERS:**
  - before (Herzberg continuum missing, bands averaged at 0.1 nm): 0.24–1.57 over 20–120 km;
  - now: 1.11–1.17 at 40–80 km, 1.25–1.33 at 90–100 km and 0.98–1.07 at 110–120 km;
  - still 1.5–2.1 at 20–30 km, where pressure-induced Herzberg absorption is missing and O3 → O(¹D) already disagrees.

**Phase 3, TUV-x comparison and a discrete-ordinates flux fix.**

- **Setup.** `tools/nlte/compare_tuvx.py` runs TUV-x (MUSICA v0.17.1, v5.4 configuration) through `tools/nlte/tuvx_reference.py` in a separate environment. Both models get the same ERS atmosphere, SZA, albedo and solar spectrum, with TUV-x aerosols off.
- **The bug it found.** The discrete-ordinates flux observers added the direct beam at the *ceiling* of the observer's layer, so actinic and downwelling fluxes were reported one model level too high. In a pure absorber this was a 9–20% error at τ = 0.5. It caused the 22–26% stratospheric O(¹D) excess seen earlier.
- **The fix** (`do_source_planeparallel.cpp`) attenuates the beam from the ceiling to the observer and carries the derivatives. The existing flux derivative tests pass, and a new exact test (`test_direct_beam_flux_is_attenuated_to_the_observer`) holds to 1e-10 at grid levels and between them. The fix also corrects `photchem.actinic_flux`.
- **After the fix, against TUV-x with the same solar spectrum:**

  | Rate | Agreement |
  |---|---|
  | O3 → O(¹D) | 0.99–1.02 over 20–120 km |
  | O3 → O(³P) | 1.01–1.04 |
  | Total O2 | 1.00–1.05 at 30–60 km and 100–110 km; 1.11–1.28 at 70–90 km |

  - The 70–90 km O2 excess is not in the flux; binned 175–200 nm flux agrees to 2%. The diagnosis below traces most of it to the comparison harness.
  - Actinic flux agrees to within 2–4% above 40 km at all wavelengths. Below 30 km in the Hartley and Schumann–Runge regions, where the flux is negligible, the models differ.
- **Solar spectrum.** TUV-x's own extraterrestrial flux gives 12.5% more O(¹D) production above 50 km than TSIS-1 HSRS.
- **Against ERS after the fix:** O3 rates within 4–7% over 20–120 km, and total O2 within 5–12% over 30–120 km.

**Phase 3, TUV-x diagnosis.** `tools/nlte/diagnose_tuvx.py` splits the comparison by wavelength region (TUV-x rates are linear in the solar flux, so it is masked one region at a time), and separates the direct beam from the diffuse flux. Atmosphere: O2, O3 and Rayleigh only.

- **Harness bug.** `tuvx_reference.py` wrote the ERS profiles into an existing TUV-x calculator through MUSICA's zero-copy arrays. That skips TUV-x's profile update, so the Lyman-α and Schumann–Runge band parameterisations kept the O2 column of the built-in v5.4 profile, which was 15–20% larger at 70–90 km. The radiators did see the new profiles. The script now builds every profile before constructing the calculator, and the inferred TUV-x O2 slant column matches ours to 0.3% at 60–85 km.
- **Lyman-α.** `LymanAlphaPhotolysis` evaluates Chabrillat & Kockarts (1997) on a straight-line spherical slant O2 column (`slant_column` in the `ActinicFlux` output), as TUV does. It agrees with TUV-x to 0.3% at 70–90 km; the constant effective cross section it replaces was 20% high at 70 km. `LinePhotolysis` now also applies the Earth–Sun distance.
- **By region, sasktran2 / TUV-x with the same sun:**

  | Region | Result | Cause |
  |---|---|---|
  | Hartley, Herzberg, Chappuis (O3) | within 1% | — |
  | Huggins O3 → O(³P) | −3% | O3 cross sections or yield |
  | Herzberg O2 | +1% | — |
  | Schumann–Runge bands O2 | 1.00 at 40 km, 1.04 at 70 km, 1.09 at 90 km | Koppers & Murtagh against resolved bands. Direct-beam band transmission agrees to 1–2%, and the result is grid-converged and only ±4% for ±20 K. |
  | Schumann–Runge continuum O2 | 0.99 at 100 km, 2.7 at 90 km | Where the continuum is optically thick, CfA against Brasseur & Solomon cross sections (7–24% apart) give large ratios. It is 15% of total O2 photolysis at 90 km and negligible below. |
  | 121.9–130 nm O2 | missing in sasktran2 | `O2UV` starts at 130 nm. 10% of TUV-x's continuum rate at 90 km, 1–2% above 100 km. |

- **Diffuse flux.** Ours is 15–20% higher in the Hartley band, 21% lower near 310 nm and 3–4% higher in the visible. Direct-beam transmission agrees to 2% down to the surface, so this is the solver: TUV-x uses 2-stream delta-Eddington. Against 16 streams at 60 km, 2-stream DO is 25% low at 290 nm but 8–11% high at 310–315 nm, the same change of sign; 4-stream DO is within 6% at 300–400 nm. The effect on photolysis rates is under 1% above 30 km. Rayleigh cross sections (Bates against Nicolet) agree to 0.6%.
- **Top boundary.** TUV-x adds the exo column to its 119–120 km layer, so its top-edge rates are unattenuated; above 110 km it is not a useful reference.
- **Solar spectrum.** TUV-x's spectrum is 10–14% above TSIS-1 HSRS in the Schumann–Runge continuum, Herzberg and Hartley regions, and within 1.3% elsewhere.
- **Result.** With the same sun, O3 → O(¹D) agrees within 1% and O3 → O(³P) within 1–3% from 30 to 120 km. Total O2 agrees within 0.4–4% from 40 to 100 km, except +16% at 90 km. Against ERS, total O2 is now within 1–6% from 40 to 120 km (12% at 90 km), down from 5–12%.

**Phase 3, TUV mode.** `TUVActinicFlux` runs at TUV's resolution with TUV's parameterisations, on the same discrete-ordinates engine:

- **Spectral and vertical treatment.** The 156 TUV-x v5.4 bins (120–735 nm), with homogeneous layers between levels, as in TUV. Layer columns assume exponential variation; the geometry uses `LowerInterpolation`, and everything goes in through a `Manual` constituent. Rates are sums over bins; `photolysis_rates` switches to them when the flux has a `wavelength_edge` coordinate.
- **O2 parameterisations.**
  - Lyman-α: Chabrillat & Kockarts.
  - The 17 SR band bins: Koppers & Murtagh.
  - Both use effective cross sections from the slant O2 column, turned into layer optical depths as in TUV-x `la_sr_bands.F90`, guards included. One deviation: where TUV's band layer optical depth turns slightly negative near the top of the fit range (about −1e-6), it is clipped to zero.
- **Data**, selected with `data=`:
  - `"sasktran2"` (default): `O3DBM`, `O2UV` and `NO2Vandaele`, averaged over each bin at 0.05 nm sampling; Bates Rayleigh; the HSRS spectrum integrated over each bin.
  - `"tuv-x"`: TUV-x's own O2, O3 and Rayleigh cross sections, extraterrestrial flux and O3 quantum yields (`TUVXQuantumYield`, `presets.tuvx_v54_photolysis`).
- **The TUV-x table.** `tools/nlte/build_tuvx_v54.py` (run with `musica`) builds `photolysis/tuvx_v54.nc`, which is now in the standard database (0.5 MB, Apache-2.0 with attribution). It doesn't regrid anything itself: it runs TUV-x once with diagnostics on and level temperatures of 180–300 K, then reads TUV-x's binned cross sections and yields. That covers every O3 temperature knot, so the table reproduces TUV-x exactly.
- **Cost.** A 151-level calculation takes 0.1–0.5 s, against several seconds on the line-resolved grid.
- **Against TUV-x** (`compare_tuvx.py`, april+00):

  | Comparison | O3 → O(¹D) | O3 → O(³P) | Total O2 |
  |---|---|---|---|
  | TUV mode, TUV-x data, 2 streams / TUV-x | within 0.2%, 30–110 km; +5% at 10 km | within 0.4%, 30–110 km | within 0.2%, 30–100 km; +1.5% at 110 km |
  | TUV mode / line-resolved, both sasktran2 data | within 0.6% | 1–2.5% low | within 1.5% at 30–70 km; 0.97, 0.91, 0.97 at 80, 90, 100 km |

  - The first row isolates the solvers. The remaining differences are where diffuse light matters (the troposphere) and at TUV-x's 120 km lid.
  - The second row is the cost of TUV's resolution and parameterisations. For O3 → O(³P), the cause is entirely the grid ending at 735 nm: 1.2–2.8% of the line-resolved rate comes from longer wavelengths. For O2, it is Koppers & Murtagh against the resolved bands, plus bin-averaged continuum cross sections where the continuum is optically thick.
- **Bug found along the way.** MUSICA arrays are views into memory owned by their map objects. The comparison scripts now keep the maps alive and copy, which explains the earlier "garbage" `wavelength_grid()` edges.

**Phase 4, first step: `sk.nlte.add_photochemical_species`.**

- **Interface.** `add_photochemical_species(atmosphere, ["O2(b)"], cos_sza=..., background=...)` takes an `sk.Atmosphere` that already has its state and absorbers. It then:
  - runs `ActinicFlux` and the oxygen rate presets;
  - solves the bundled mechanism;
  - adds a `PopulationEmissionRate` constituent (`"O2(b) emission"`: the A band with its 1-1 hot band, and the B band);
  - returns the solution and rates for inspection.

  Inputs come from `background`, then from the atmosphere's `VMRAltitudeAbsorber` constituents, then constant N2 and CO2 mixing ratios. Atomic oxygen must be in `background`.
- **Engine setting.** Emission needs `config.emission_source = sk.EmissionSource.VolumeEmissionRate`; the function warns otherwise.
- **Not yet available.** `O2(a)` (no a-X band emission yet) and `O(1S)` (the mechanism has no source) raise `NotImplementedError`. The γ band (2-0) is not in `PopulationEmissionRate` yet.
- **Checks.**
  - Optically thin limb radiance equals ∫VER ds/4π to 0.01%.
  - A daytime ERS limb case (april+00, 755–775 nm at 0.001 nm) takes 10 s for the photochemistry. The A band is 7× the in-band Rayleigh radiance at a 60 km tangent and 290× at 90 km.
- **Against GRANADA** (`tools/nlte/validate_oxygen_ers.py`, april+00). The O2 conversion assumes electronic degeneracies in GRANADA's LTE weights; the level file that would confirm this is not in the archive.

  | Altitude | O(¹D) | O2(a) | O2(b) |
  |---|---|---|---|
  | 40–110 km | 0.83–1.06 | 0.60–0.93 | 3.4–4.7× too high at 40–70 km; 1.0–1.5 at 90–110 km |

  O(¹D) agreeing confirms the photolysis source. The O2(b) excess comes from a missing reaction: the legacy mechanism has no O2(b, v=0) + N2 quenching.

**Phase 4, the `oxygen` mechanism.**

- **Mechanisms.** `oxygen_yankovsky` is frozen as the exact legacy translation. The new bundled `oxygen` mechanism is audited process by process against:
  - the newest Yankovsky database (Yankovsky et al. 2019 supplement, Yankovsky and Vorobeva 2020);
  - JPL 19-5 (unchanged in 25-1);
  - HITRAN.

  The audit is in `nlte_oxygen_audit.md`. `add_photochemical_species` uses `oxygen`.
- **Mechanism changes.**
  - It takes the physical rates of `presets.oxygen_photolysis`; product yields live in the mechanism. `photo_reaction` entries can now branch into `channels`.
  - O2(X, v) and O(1S) are dropped.
  - O2(b, v=0) + N2 is added.
  - b-X Einstein coefficients come from HITRAN (A band 8.75e-2 s⁻¹, B band 7.34e-3 s⁻¹).
- **Emission.** Band emission now comes from the mechanism's transition VERs through `O2BandEmissionRate`: A band 0-0 and 1-1, and B band 1-0. The Einstein coefficients therefore live in one place.
- **B-band fix.** `PopulationEmissionRate` used 7.0e-2 s⁻¹ (the 1-1 value) for the B band; it now uses 7.34e-3.
- **Against GRANADA** (three day scenarios):
  - O(¹D): 0.95–1.0 at 40–60 km.
  - O2(a): 0.81–0.88 at 40–60 km and 0.69–0.71 at 70–80 km.
  - O2(b): 0.51–0.58 at 40–80 km.

  The O2(b) ratio stays constant while A-band pumping goes from 9% to 80% of its production, which points to GRANADA's loss term (N2 quenching) or its population convention. GRANADA's rate files are not available to check.

**Phase 3, remaining:**

- O2 cross sections for 121.9–130 nm (between Lyman-α and `O2UV`).
- Pressure-induced Herzberg absorption.
- SZA > 90° (twilight).
- Per-transition O2 excitation from HITRAN-filtered lines, with a hybrid high-resolution direct beam.
- NO2, H2O and NO δ-band reactions.
- Switching `Yankovsky` to computed rates, and deprecating `photchem`.

## Phases

Each phase is one or two reviewable PRs.

**0. Staging and goldens**
- Capture `photchem` golden outputs as test fixtures.
- Stage the remaining datasets on T9 with READMEs. ERS is already staged.
- Write a small reader for ERS `.prf` files to support validation.

**1. Crate skeleton and mechanical migration.** No behaviour change.
- Acceptance:
  - `tests/photchem` passes unchanged;
  - goldens match to round-off;
  - `cargo clippy --all-targets --all-features` is clean;
  - the crate's `Cargo.toml` has no SASKTRAN2 dependencies.

**2. Generic mechanism engine**
- TOML schema, loader and validation; rate-law library; steady-state and transient solvers; budgets; the `sasktran2.nlte` Python API; Yankovsky as TOML; `photchem` shims.
- Acceptance:
  - parity with the goldens;
  - per-rate-law unit tests against hand-computed values;
  - closed-form tests (a two-state system, conservation).

**3. `sasktran2.photolysis`**
- Actinic flux API, reaction library, line excitation rates, twilight check.
- Acceptance:
  - in an optically thin TOA case, J equals the integral of σφF;
  - agreement with TUV-x and ERS `npar` within documented tolerances;
  - Yankovsky driven by computed J values.

**4. Emission glue**
- `LineRecord` adapters, rotational models, templates, per-species mass in line-list VER.
- Acceptance:
  - optically thin limb radiance matches the VER line integral divided by 4π;
  - existing A-band tests pass.

**5. Oxygen complete**
- O2(b) and O2(a) day and night, O(¹S) day and night, Herzberg/Chamberlain bands.
- O(¹D) only if ionospheric inputs are available.

**6. OH**
- Meinel bands and A–X dayglow.

**7. NO**
- γ/δ dayglow, N+O nightglow, NO2 continuum.

**8. Later**
- IR vibrational non-LTE (see `nlte_granada_review.md`).

## Open questions and risks

- **Twilight.** DO pseudo-spherical support for SZA > 90° needs checking. A spherical direct beam may be required.
- **Cost.** High-resolution windows are expensive; the O2 bands at 0.001 nm already need about 75k points. The hybrid direct-beam scheme is the mitigation.
- **Radiative coupling within bands is neglected.** This means earthshine and re-absorbed airglow in the A-band. The error should be quantified once the basic chain works.
- **Rotational non-LTE.** The OH high-N tail and NO γ fluorescence structure need line-level hooks from the start.
- **Data gaps.**
  - NO δ/β nightglow rates.
  - O2 Herzberg/Chamberlain spectroscopy.
  - Ionospheric inputs for O(¹D).
- **Twilight consistency.** Prescribed background chemistry (O3, H, O) must be consistent across day and night at twilight.
- **Transient input format.** O2(a) transient runs need a time history of J. The input format is still to be defined.
