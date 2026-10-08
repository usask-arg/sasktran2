# GRANADA Review (IR Non-LTE Reference)

This note reviews the GRANADA non-LTE population code (Funke et al. 2012, JQSRT 113, 1771) as shipped in the KIT `kopra_granada_rcp` repository. It records what would matter for a later IR vibrational non-LTE extension of SASKTRAN2. The active plan, focused on UV, visible and near-IR emission, is in `nlte_plan.md`.

GRANADA references below are relative to the root of that repository, e.g. `Kopra/source_10.0.0/modules/nltmodel_m.f90:3228`. SASKTRAN2 references are relative to this repository.

## Summary

- **Size and layout.** GRANADA is about 44k lines of Fortran 90. The `Granada/` directory holds only a 71-line driver. The physics lives in the KOPRA library:
  - `nltm*` modules: statistical equilibrium and the process library;
  - `nmrt*` modules: radiative transfer for the excitation rates;
  - `ifnlte_m`: the interface back to KOPRA's limb model.
- **Algorithm.**
  - The prognostic variable is the population ratio $r = n/n_\mathrm{LTE}$ per vibrational level and altitude.
  - The statistical equilibrium equations (SEE) are solved per altitude.
  - Lambda iteration handles weak or overlapping bands; a Curtis-matrix solve handles optically thick bands.
  - Spectra are computed by rescaling precomputed LTE line opacities with smooth NLTE factors.
- **Weaknesses of the implementation:**
  - about 85% of the 8.6k-line process library is near-duplicate parameterisations;
  - a hand-configured iteration schedule;
  - species coupled by lagged block Gauss–Seidel;
  - local finite-difference Jacobians with the radiation field frozen;
  - module-level mutable state that prevents running cases in parallel;
  - several bugs (listed below).
- **Licence.** GRANADA and KOPRA are LGPL-2.1, while SASKTRAN2 is MIT. Any reimplementation must be clean-room. Reading GRANADA output files, such as the CAIRT ERS ratio files, is unaffected.
- **Relevance to the UV/visible plan.** The process library is mostly IR vibrational kinetics. A few parts are useful literature pointers:
  - the O2 electronic chain (O3 photolysis → O2(a,b), O(¹D) → O2(b));
  - OH(v) relaxation;
  - the nascent OH(v) distribution from H + O3;
  - NO(v) from N + O2.

## What GRANADA does

### Unknowns and level model

- **Level identity.** A level is `(HITRAN molecule*10 + isotope, vibrational class index)`. Its quantum numbers are parsed from HITRAN global quanta according to a per-molecule "vibrational class" (`nltmvibclass_m.f90`).
  - Class 5 is CO2 (`v1 v2 l2 v3 r`). Class 6 covers H2O and O3. Class 3 covers NO, OH and ClO with spin and electronic labels. Class 2 is O2 with electronic labels. Class 10 covers CH4 and HNO3.
- **Level data** come from a "spectroscopy file" (`nltminput_m.f90:2547-2895`): energy, degeneracy, and per-transition Einstein A.
  - Band Einstein coefficients are optionally recomputed from HITRAN line sums as temperature-dependent (altitude-dependent) values (`nmrtabco_m.f90:3388-3500`).
- **Unknown.** The fractional population $x_k = n_k/N$ of each level. The stored quantity is $r_k = x_k / x_k^\mathrm{LTE}$.
  - $x^\mathrm{LTE}$ is normalised over only the levels in the state vector, not the full vibrational partition function (`nltmodel_m.f90:1153, 1241`).
- **Level groups.** Levels can be grouped and assumed in mutual LTE: summed $g$, a user-supplied effective energy, and g-weighted A (`nltminput_m.f90:1284-1373`).
- **Rotational NLTE** (NO, OH, CO) is optional, applies above `alt_rot`, and adds one unknown per (v, spin, J).

### Statistical equilibrium

Per altitude, GRANADA assembles a small dense system $M x = b$ (`see_calc`, `nltmodel_m.f90:416`):

| Contribution | Source | Form |
|---|---|---|
| Radiative | `rates_rd`, `:1940` | Emission $A + B\,\bar J$, absorption $(g_u/g_l) B \bar J$. Band B is corrected for NLTE stimulated emission. |
| Vibrational–translational and vibrational–vibrational | `rates_vt`, `:2304`, through `proclib_vt` | Rate coefficient × collider density. The reverse rate comes from detailed balance at kinetic $T$. VV terms are multiplied by the partner's population, which makes them bilinear across species. |
| Rotational–translational | `rates_rt`, `:2249` | Energy-corrected power-gap law or strong-collision model. Rotational NLTE only. |
| Chemistry and photolysis | `rates_ch`, `:2659` | Production goes to the right-hand side with a nascent vibrational distribution. Loss goes to the diagonal. |

The system is solved in one of two ways:

- **Gauss–Seidel**, ordered by descending energy, followed by renormalisation.
- **LAPACK `dgesv`**, with the row for the lowest state replaced by number conservation (`:608-665`). That row replacement silently drops an excited level's balance if the ground state is not in the sub-calculation.

### Radiative excitation

- **Exchanged quantity.** Per band, `rflx` $= 4\pi c\,\bar J$, where
  $$\bar J_b(z) = \int \phi_b(\nu, z)\, J(\nu, z)\, d\nu, \qquad \phi_b = \frac{k_b(\nu)/\nu}{\int k_b/\nu\, d\nu}.$$
  - $k_b$ is the NLTE absorption coefficient of the target band.
  - The field is split into upwelling, downwelling and direct-solar parts (`nmrtradtra_m.f90:812-937`).
- **Geometry and grids.**
  - Plane-parallel diffuse field with $N$ midpoint-rule zenith angles, or a single diffusivity stream (`init_ang`, `:3963`).
  - KOPRA's irregular line-by-line grid.
  - Scattering is ignored, apart from optional Lambertian reflection of the direct sun.
- **Lambda iteration** (types 10–14): line-by-line $\bar J$ with the total opacity of all species in the microwindow, so overlap is exact. Populations from the previous iteration are used.
- **Curtis matrix** (types 20–32): built analytically from transmission differences, not by perturbation (`radtrans_cm`, `:1115`).
  - $\bar J_k = \sum_j C_{kj} R_j + v_{1,k} F_\uparrow + v_{2,k} F_\downarrow + v_{3,k} F_\odot$, with $R_j = n_u g_l/(n_l g_u)$.
  - Only the target band's own opacity is included.
  - The SEE for one upper state is then solved jointly across altitude (`curtis_matrix`, `nltmodel_m.f90:3228`, paper eqs. 24–25).
  - C is rebuilt once per main iteration.
- **"Comb" variant** (types 30–32): lines are binned by $\log S$ into pseudo-lines with Curtis–Godson-averaged widths.
- **Layer source.** A heuristic blend, $S_\mathrm{eff} = w S_\mathrm{exit} + (1-w) S_\mathrm{CG}$ with $w = 1 - \mathcal{T}^{0.3}$ (`:872`).
- **Solar pumping.** Straight-line spherical paths without refraction. The slant optical depth is built as airmass factors times vertical layer optical depths (`:940-1080`). The solar spectrum is a brightness-temperature parameterisation with optional Fraunhofer lines.

### Iteration control

- The user configures it all in the setup file:
  - "calcs": independent runs, executed in sequence;
  - SEE and RTE "subcalcs";
  - ordered "groups" of subcalcs, each sub-iterated to convergence;
  - an outer main iteration.
- Species are always solved one at a time. VV partners use populations from the previous pass.
- The only joint multi-species step is an optional elimination that makes one VV partner state implicit inside a Curtis-matrix solve (`nltmodel_m.f90:482-545, 3269-3473`).
- Non-convergence only produces a warning.

### Jacobians

- `deri_calc` (`nltmodel_m.f90:4498`) differentiates the converged local SEE implicitly: $dx = M^{-1}(\partial b - \partial M\,x)$.
- $\partial M$ and $\partial b$ come from:
  - ±0.5 K central differences for temperature;
  - ±1% for the target VMR, which only affects chemical dilution;
  - linearity in the scale factor for process scales;
  - ±0.5% central differences for process parameters.
- During the temperature perturbation the radiation field is replaced by its local Curtis-diagonal part.
- **The Jacobians are therefore diagonal in altitude.** They include no $\partial \bar J/\partial T$ from other altitudes, no cross-species response, and no temperature dependence of the band coefficients or of C.

### What KOPRA consumes

- **Interface.** `ifnlte_m` passes $r$ per KOPRA NLTE state per altitude. The arrays are named `Tvib`, but they hold $r$. It also passes the Jacobians `dr_dT`, `dr_dvmr` and `dr_dx`.
- **Optional extras:**
  - a temperature parameterisation for horizontal gradients, $r(T) = r_\mathrm{th} + (r_0 - r_\mathrm{th}) e^{c_2 E_\mathrm{eff}(1/T - 1/T_0)}$;
  - five additional GRANADA runs at SZA, SZA ± 2.5° and SZA ± 10°, interpolated by local SZA along the line of sight.
- **Per line in the limb radiative transfer** (`radtra_m.f90:2330-2364, 5883-5925`), with $x = c_2\nu/T$:
  $$k = k^\mathrm{LTE}\,\alpha,\quad \alpha = \frac{r_l - r_u e^{-x}}{1 - e^{-x}},\qquad j = k^\mathrm{LTE}\, r_u\, B_\nu(T),\qquad S = \frac{j}{k} = \frac{2hc^2\nu^3}{(r_l/r_u)e^{x} - 1}.$$
  - Negative opacity from population inversion is handled explicitly.
  - Rotational LTE at the kinetic temperature is assumed within each vibrational level.

### Process library and inputs

- **Process library.** `nltmproclib_m.f90` (8642 lines) dispatches about 150 vibrational-translational / vibrational-vibrational parameterisations and about 35 chemical ones by name string, each with a positional `par(:)` array.
  - The distinct physics is roughly:
    - 15 temperature-dependence forms: constant, Arrhenius, power law, Landau–Teller $a + b e^{c/T^{1/3}}$, `exppow`, tabulated, and so on;
    - about 10 selection-rule idioms on $\Delta v_i$, polyad and electronic state;
    - 5 quantum-number scaling rules;
    - 4 nascent-distribution families;
    - about 8 bespoke schemes: CO2 Fermi mixing and branching, O3 Landau–Teller/SSH, the CH4 cascade, OH Varandas, the O/O3 photochemical-equilibrium side effect, and Mars tables.
  - A table-driven evaluator would be roughly 1k lines plus data files.
- **Inputs.** A KOPRA-style `$x.y` main input, a setup file per calc, a process file per species, a spectroscopy file, profile files (p, T, VMR, partner $r$, parameter profiles such as $J(\mathrm{O_3})$ and $[\mathrm{O(^1D)}]$), and a solar spectrum.
- **Missing data.** The repository contains no setup, process or spectroscopy files and no reference outputs.
  - The CAIRT ERS archive (Zenodo 8256025) provides GRANADA outputs and inputs: ratio files, p/T, VMRs including OH, HO2, N(⁴S), N(²D), H and O(¹D), and TUV photolysis rates in `npar` files.
  - It does not include the process or setup files.

## Assessment

### Keep (good design)

- **$r$ as the prognostic variable**, with LTE line shapes computed once and NLTE entering only as per-line multipliers ($\alpha$, $r_u$).
- **Detailed balance** applied automatically for every reverse rate, including VV partner populations.
- **A Curtis-matrix (linearised) solve for optically thick bands.** Lambda iteration stalls for CO2 15 µm.
- **Level grouping**, and band Einstein coefficients that depend on temperature, computed from line sums.
- **A per-band, profile-weighted $\bar J$** as the radiative coupling quantity.

### Improve (do not port)

| GRANADA | Replacement |
|---|---|
| Hand-configured calc/subcalc/group schedule; sequential species with lagged VV coupling | One joint Newton solve over all coupled levels, species and altitudes |
| Curtis matrix uses only the target band's opacity; overlap handled only by Lambda iteration | Cross-band Curtis matrices computed with the total opacity |
| Midpoint-rule angles plus the empirical $\mathcal{T}^{0.3}$ source blend | Exact angular integration (exponential integrals) with linear-in-$\tau$ layer sources, or Gauss–Legendre in $\mu$ |
| Local finite-difference Jacobians with $\bar J$ frozen | Implicit-function-theorem Jacobians from the converged Newton system, non-local and exact |
| Name-dispatched, hard-coded rate routines | Table-driven rate forms plus a few named plugins; analytic $dk/dT$ |
| Five-SZA ratio trick and the $(r_\mathrm{th}, E_\mathrm{eff})$ temperature parameterisation for limb inhomogeneity | Independent 1D column solves feeding a 2D population field into Geometry2D |
| Module-level mutable state | Pure functions, parallel over columns and SZA |

## Issues found in GRANADA

I confirmed these in the source:

1. **`dα/dT` has the wrong sign** in GRANADA's `alphan` (`nmrtradtra_m.f90:3199, 3202, 3205`) compared with KOPRA's correct copy (`radtra_m.f90:5774-5780`). The correct form is $d\alpha/dT = e^{-x}(x/T)(\alpha - r_u)/(1-e^{-x})$. The bug is latent: GRANADA never requests `dx_dT`.
2. **`radv2r` copies the upwelling component into all three field components**, and has no "band not found" skip, so rotational lines of other bands receive the last checked band's field (`nltmodel_m.f90:4205-4226`).
3. **Out-of-bounds write in `deri_calc`.** The `if (iprof/=0)` guard covers only the `poplte` assignment, so the `Tvib(iprof)` write runs with index 0 (`nltmodel_m.f90:4619-4624`, repeated at `4670-4675` and `4719-4723`).
4. **Derivative skip test uses a signed comparison.** `any(aux2_mat > 1e-60)` skips the $\partial M\,x$ term when $\partial M$ has only non-positive entries, as with pure-loss processes (`nltmodel_m.f90:4826`; `4953`).

The subagent reviews reported these but I have not re-verified them:

- The VV rate cache is keyed without the partner index (`nltmproclib_m.f90:428-429, 805-819`).
- Uninitialised variables in `vt_expB_ho_v2v4` and `vt_pow_v3_hb_deltaE`, and `par(6)` reused for two meanings in `vt_pow_v3tov2_branch_t`.
- The chemistry temperature derivative reuses cached rates computed at the unperturbed temperature.
- VV elimination assumes the partner is row 1.
- `MakeSolspec` indexes the regular grid instead of the irregular one (`nmrtsolar_m.f90:287`).

These are worth reporting upstream to KIT/IAA.

### Comparison gotchas

When comparing against GRANADA output:

- **$r$ normalisation.** GRANADA's $r$ is relative to a truncated partition function. SASKTRAN2 line strengths use the full TIPS $Q(T)$. Convert with $r_\mathrm{full} = r_\mathrm{GRANADA}\,Q_\mathrm{vib}/\sum_{k \in \mathrm{state\ vector}} g_k e^{-c_2E_k/T}$.
- **Units.** GRANADA's `rflx` is $4\pi c\bar J$, with $B = A/(8\pi h c^3\nu^3)$.
- **Naming.** KOPRA's `Tvib` arrays hold $r$, not temperatures.

## Future IR extension outline

If IR vibrational non-LTE is added later, it belongs in the `sasktran2-nlte` crate. It would reuse the state, process and mechanism model from `nlte_plan.md` and add five things.

### 1. NLTE line absorber in `sasktran2-rs`

- Per-line $\alpha$ is applied where line strengths are scaled (`rust/sasktran2-rs/src/optical/line/db.rs:126`, `adjusted_parameters`).
- The absorber accumulates the excess emission $\Delta j = k^\mathrm{LTE} B_\nu(T)(r_u - \alpha)$, which is identically zero when $r \equiv 1$.
- This needs an engine emission mode that combines the LTE source function (`standard`) with an additive emission coefficient (`volume_emission_rate`).

### 2. Radiative-rate kernel in the crate

- No-scattering mean-intensity operator, profile-weighted per band.
- Cross-band Curtis matrices computed with the total opacity.
- Solar airmass factors from SASKTRAN2's Geometry1D ray tracing.
- Validated against DO actinic flux with thermal emission.

### 3. Joint Newton solver

- Unknowns: populations of every level and species at every altitude.
- Curtis-linearised radiative terms.
- Lower-state opacity updated in an outer iteration.

### 4. Jacobians

- From the converged system via the implicit function theorem.
- Delivered as one derivative mapping per level, with a dense interpolator $\partial r_\ell(z)/\partial T(z')$. Interpolators are dense `Eigen::MatrixXd` (`cpp/include/sasktran2/derivative_mapping.h:87`), so no engine change is needed.

### 5. Line-to-level mapping

- HITRAN global-quanta parsers per class.
- HITRAN loader records already keep quanta, $A$ and $g$ (`rust/sasktran2-rs/src/optical/line/db.rs:11-31`).
