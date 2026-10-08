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
