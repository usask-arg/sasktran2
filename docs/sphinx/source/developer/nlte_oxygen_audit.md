# Oxygen mechanism audit

This page records the choices behind the bundled `oxygen` mechanism (`rust/sasktran2-nlte/mechanisms/oxygen.toml`). That mechanism replaces `oxygen_yankovsky` for new work. `oxygen_yankovsky` stays frozen as the exact translation of the legacy photochem model, for its golden tests and the `Yankovsky` shim.

## Sources

| Key | Source | Use |
|---|---|---|
| ym2019 | Yankovsky, Vorobeva and Manuilova (2019), Adv. Space Res. 64, 1948, supplement Tables S1 and S2, "status on March 2019" | Newest complete rate and yield table of the Yankovsky model. It covers O(1D), O2(b, v=0-2) and O2(a, v=0-5). Read in full; the main text is paywalled and was not read. |
| yv2020 | Yankovsky and Vorobeva (2020), Atmosphere 11, 116 | Newest description of the model ("version 2019"): O2(a, v) yields from O3 photolysis, measured Einstein coefficients, and radiative lifetimes. Read in full. |
| ym2006 | Yankovsky and Manuilova (2006), Ann. Geophys. 24, 2823 | The source of `oxygen_yankovsky`: every value there matches its Tables 1-5, except the errors listed below. |
| jpl19 | JPL 19-5 (Burkholder et al. 2019) | Rate evaluations. The entries used are unchanged in JPL 25-1 (2025). |
| hitran | HITRAN2020 16O2 lines | b-X band Einstein A values, computed from the line list. |
| nist | NIST ASD | O(1D) radiative rates. |

No newer daytime kinetics from the Yankovsky group were found for 2021-2026.

**Rules:**
- JPL for rates it evaluates.
- ym2019 for products and yields, and for processes JPL does not cover.
- HITRAN for the b-X band Einstein coefficients, so the mechanism's band VERs match the line lists used for absorption and emission.

## Scope changes from `oxygen_yankovsky`

- **Physical photolysis inputs.** The rate inputs are `presets.oxygen_photolysis` (`J_O3_O1D`, `J_O2_SRC`, `J_O2_LYA`, `J_O2_EXC_*`), with product distributions as yields in the mechanism. The legacy file instead takes 47 rates pre-split by product level.
- **No O2(X, v).** Vibrationally excited ground-state products go to O2. Nothing feeds back from them to O2(a) or O2(b). Their legacy rates had the errors below, and the newest model replaces them with a multi-quantum cascade that is not tabulated.
- **No O(1S).** It had no source in the legacy file, the Yankovsky model has none, and daytime sources (O2 photolysis below 133 nm, photoelectrons, O2+ recombination) have not been reviewed.
- **Not included:** H2O quenching of O2(b) (about 2% of the loss at 6 ppmv), and reverse energy transfer (negligible for these states).

## Errors found in `oxygen_yankovsky`

These are kept in that file for parity:

| Entry | Problem | Effect |
|---|---|---|
| none | O2(b, v=0) + N2 quenching is missing | The main O2(b) loss below 90 km is absent, so daytime O2(b) is 3-5× too high at 40-80 km |
| `o3_o_v{v}` | 5.6e-11 instead of ym2006's 5.6e-12, and every v channel runs at the full rate (yields sum to 30) | O2(X, v) only |
| `o2x{v}_o2_vv` | 2.6e-13 applied to v=3-35 on top of the v=4-20 family | O2(X, v) only |
| `o3_o2x_v{v}` (preset) | O3 → O(3P) shared equally over O2(X, v=1-35) | O2(X, v) only |
| `o2b1_o` | Products O2(b, v=0); ym2019 has O2 + O(3P) | Above 90 km only |

There is also an emission-side error, now fixed: `PopulationEmissionRate` used A = 7.0e-2 s⁻¹ for the B band (b-X 1-0), which is the 1-1 value. HITRAN gives 7.34e-3, close to the measured 7.2e-3. B-band emission was therefore 10× too bright.

## Process table

Rates in cm³ s⁻¹ (T in K) unless noted; ✓ marks the value chosen.

### Photolysis and excitation

| Process | Legacy | ym2019 / yv2020 | JPL 19-5 | `oxygen` |
|---|---|---|---|---|
| O3 → O(1D) + O2(a, v) | 6 inputs; yields 0.49, 0.15, 0.15, 0.08, 0.08, 0.05 (ym2006, 254 nm) | λ-dependent (yv2020 Eqs. 3-4); at 254 nm 0.362, 0.276, 0.115, 0.075, 0.076, 0.095 | Φ(O1D) only | ✓ yv2020, wavelength-dependent: one rate per level (`J_O3_O1D_A{v}`), plus `J_O3_O1D_X` for the spin-forbidden O(1D) + O2(X) channel beyond 310 nm. Eq. 2 of yv2020 does not reproduce its Table 1, so x(λ) is interpolated between the tabulated thresholds. O2(a, v≥1) relax to v=0 within milliseconds, so only the vibrational distribution depends on this. |
| O2 SRC → O(3P) + O(1D) | yield 1 | yield 1 | Φ = 1 at 139-175 nm | ✓ |
| O2 Lyman-α → O(1D) | Φ = 0.53 in the rate | 0.48-0.58 | 0.44 ± 0.05 at line centre | ✓ Φ = 0.53 in `J_O2_LYA` |
| O2 + hν → O2(b, v=0-2), O2(a) | line-by-line | top-of-atmosphere 5.35e-9, 2.94e-10, 7.94e-12, 1.54e-10 s⁻¹ (ym2006) | — | ✓ line-by-line; ours at the top of the atmosphere: 6.2e-9, 3.7e-10, 1.2e-11, 1.2e-10 |

### O(1D)

| Process | Legacy | ym2019 | JPL 19-5 | `oxygen` |
|---|---|---|---|---|
| + O2 | 3.2e-11 exp(67/T); b(v=1) 0.40, b(v=0) 0.55, a 0.05 | 3.3e-11 exp(55/T); b(v=1) 0.8, b(v=0) 0.2 (Pejakovic et al. 2014) | 3.3e-11 exp(55/T); total b 0.8 ± 0.2 | ✓ JPL rate and yields: b 0.8 split 4:1 between v=1 and v=0 (0.64, 0.16), a(v=0) 0.2. |
| + N2 | 2.0e-11 exp(107/T) | 2.15e-11 exp(110/T) | 2.15e-11 exp(110/T) | ✓ JPL (+9% vs legacy) |
| + O3 | 2.4e-10 → 2 O2 | 2.4e-10, 50% → O2 + 2O | 2.4e-10, 50/50 | ✓ JPL |
| + CO2 | — | 7.5e-11 exp(115/T) | 7.5e-11 exp(115/T) | ✓ JPL (new; negligible) |
| + O(3P) | 4.0e-12 | 2.2e-11 (Kalogerakis et al. 2009) | — | ✓ ym2019 (matters above 95 km) |
| radiative | 9.0e-3 s⁻¹ | — | — | ✓ NIST 630.0 nm 5.63e-3 and 636.4 nm 1.82e-3 s⁻¹ |

### O2(b, v)

| Process | Legacy | ym2019 / yv2020 | JPL 19-5 | `oxygen` |
|---|---|---|---|---|
| b0 radiative | A band 7.58e-2 s⁻¹ | measured 0-0 8.93e-2, 0-1 4.67e-3, b-a 1.2e-3 | — | ✓ HITRAN 0-0 8.75e-2; yv2020 0-1 and Noxon |
| **b0 + N2** | **missing** | 2.2e-15 (Yankovsky and Manuilova 2018); products 0.5 a(v=2), 0.5 X(v=9) | 1.8e-15 exp(45/T) | ✓ JPL rate, ym2019 products |
| b0 + O2 | 3.9e-17; a(v=0-3) yields | same (Klingshirn and Maier 1985) | 3.9e-17 | ✓ |
| b0 + CO2 | 4.2e-13 → a0 | 4.4e-13; a ≥ 0.9 | 4.2e-13 | ✓ JPL, yield 1 |
| b0 + O3 | 2.2e-11; 0.3 → a0 | 3.5e-11 exp(-135/T); 0.7 → 2O2 + O | 3.5e-11 exp(-135/T) | ✓ JPL rate, ym2019 channels |
| b0 + O(3P) | 8e-14; 0.75 → a0 | same (Hadj-Ziane et al. 1992) | 8e-14 (factor 5) | ✓ |
| b1 radiative | 7.0e-2 s⁻¹ (to X(1)) | 1-0 7.20e-3, 1-1 7.01e-2 | — | ✓ HITRAN 1-0 7.34e-3, 1-1 6.99e-2 |
| b1 + O2 → b0 | 4.2e-11 exp(-312/T) | same (Hwang et al. 1999) | — | ✓ |
| b1 + N2 → b0 | 5.0e-13 | < 7e-13 | — | ✓ 7e-13 (limit used as value) |
| b1 + O3 | 3.0e-10 | < 3e-10 | — | ✓ 3e-10 |
| b1 + O(3P) | 4.5e-12 → b0 | 4.5e-12 → O2 | — | ✓ ym2019 |
| b1 + CO2 → b0 | — | 9e-13 (220 K) | — | ✓ |
| b2 radiative | 5.4e-2 s⁻¹ (to X(2)) | lifetime 11.75 s | — | ✓ HITRAN 2-0 (γ) 2.56e-4, 2-1 9.7e-3; remainder 7.5e-2 |
| b2 + O2 → b0 | 1.2e-11 exp(-596/T) | 2.3e-11 exp(-691/T) (Hwang et al. 1999) | — | ✓ ym2019 |
| b2 + O(3P) | 1.1e-11 → b1 | 1.07e-11 → b1 0.5, b0 0.5 | — | ✓ ym2019 |
| b2 + O3 | 2.9e-10 | same | — | ✓ |
| b2 + N2 | 2.0e-14 → b1 | 8.0e-15 (110 K) → b0 | — | ✓ ym2019 |
| b2 + CO2 → b1 | — | 3.0e-12 exp(-158/T) | — | ✓ |

### O2(a, v)

| Process | Legacy | ym2019 / yv2020 | JPL 19-5 | `oxygen` |
|---|---|---|---|---|
| a0 radiative (1.27 µm) | 2.58e-4 s⁻¹ | measured 2.26e-4 | — | ✓ 2.26e-4 |
| a0 + O2 | 3.6e-18 exp(-220/T) | same | same | ✓ |
| a0 + O3 | 5.2e-11 exp(-2840/T) | same | same | ✓ |
| a0 + O(3P) | 6.5e-17 | 1e-16 (Yankovsky et al. 2016; factor 3) | < 2e-16 | ✓ ym2019 |
| a0 + N2, + CO2 | 1e-20, — | ≤ 1.4e-19, ≤ 2e-20 | < 1e-20, < 2e-20 | omitted (upper limits, negligible) |
| a(v) + O2 → a0 | 5.6e-11 (v=1), 3.6e-11 (v=2-5) | same (Pejakovic et al. 2011) | — | ✓ |
| a1 + O3 | 4.7e-12 | same (Klais et al. 1980) | — | ✓ |
| a(v≥1) + O(3P) | 1e-14 (guess) | < 4e-13 (v=1), no data above | — | omitted (no data, negligible against O2) |

## Against GRANADA

GRANADA (Funke et al. 2012) produced the CAIRT ERS populations. Its source on T9 contains only rate-law shapes. The O2 rate constants, yields, level energies and degeneracies live in input files that are not in the code, the ERS archive or the git history. What could be established:
- O(1D) is an input to GRANADA, from ERS v7; it is not solved for.
- O(1S) is not modelled.
- Its O3 → O2(a, v) yields differ from the ones above. Above 50 km, a0-a5 = 0.555, 0.263, 0.104, 0.046, 0.018, 0.015.
- Its exchange processes are reversible.
- Its LTE weights very likely include the electronic degeneracy (X 3, a 2, b 1). With degeneracy 1 for all states, the implied A-band pumping would be 5-6× too high. The comparisons use this convention.

`tools/nlte/validate_oxygen_ers.py`, with the `oxygen` mechanism, for three day scenarios (april+00, july+45, october-45):

| Altitude | O(1D) | O2(a) | O2(b) |
|---|---|---|---|
| 40-60 km | 0.95-1.00 | 0.81-0.88 | 0.51-0.55 |
| 70-80 km | 0.90-1.07 | 0.69-0.71 | 0.51-0.58 |
| 90-100 km | 0.83-1.19 | 0.58-0.83 | 0.68-0.85 |

The O2(b) ratio is nearly constant from 40 to 80 km. Over that range, the share of O2(b, v=0) production from A-band pumping rises from 9% to 80%, while N2 quenching stays 79-89% of the loss. A factor of about 1.9 common to all of it points to the loss, or to the population convention, not to production.

GRANADA's O2(b) + N2 rate would explain it if it were about half the JPL and ym2019 values (which agree to 5%). Its value cannot be read without the input files.

### GRANADA's published scheme

Funke et al. (2012), Table 4, lists GRANADA's O2 processes, taken mainly from ym2006. The O2(a)/O2(b) values that matter at 40-80 km agree with ours:

| Process | Funke et al. (2012) | `oxygen` |
|---|---|---|
| b0 + N2 | 1.05e-15 → a(v=2) + N2(1), plus 1.05e-15 → X(v=9); total 2.1e-15 | 1.8e-15 exp(45/T) (2.2e-15 at 200 K) |
| b0 + O2 | total 3.85e-17, a(v=0-3) | 3.9e-17 |
| b0 + CO2 | 4.2e-13 → a0 | same |
| b0 + O3 | 6.6e-12 → a0 | 3.5e-11 exp(-135/T) (1.8e-11 at 200 K; at most 6% of the loss) |
| b0 + O | 2.0e-14 → X, 6.0e-14 → a0 | same |
| b1 + O | 4.5e-12 → b0 | → O2 (only above 90 km) |
| b2 + O2 | 1.2e-11 exp(-596/T) | 2.3e-11 exp(-691/T) |
| a0 + O | 6.5e-17 | 1e-16 |

The O(1D) + O2 yields are not tabulated; the code takes them from input files.

The paper also states that radiative transfer in the a-X and b-X bands is computed line by line, with a modified Curtis matrix. That includes absorption of upwelling and emitted band radiation, which we do not include, but radiation is only 1-12% of the O2(b) loss at 40-80 km.

With the loss terms matching, the factor of about 1.9 does not come from N2 quenching.

**The population convention is not the cause.** Fig. 9 of the paper plots daytime number densities directly:
- O2(b, v=0) is about 1.5e6 cm⁻³, nearly flat from 30 to 90 km;
- O2(a, v=0) is 4-5e9 cm⁻³ at 40-60 km.

These agree with the ERS ratios converted with electronic degeneracies, so that convention is right.

**The difference is the A-band pumping.** Where N2 quenching dominates the loss, n_b ≈ J [O2] / (k_N2 [N2]). For about 1.5e6 cm⁻³ that needs J ≈ 1.2e-8 s⁻¹ at every altitude from 30 to 90 km. That is:
- about 2× our top-of-atmosphere rate of 6.2e-9 s⁻¹, which lies within published values (5.35e-9 in ym2006);
- about 10× ours at 40 km, where the strong A-band lines are optically thick to the direct beam.

Our 0.52-0.58× ratio is constant because pumping in GRANADA dominates O2(b) production at all these altitudes. Ours is O(1D)-dominated at 40 km and pumping-dominated at 80 km. We keep our line-by-line pumping; GRANADA's is not reproduced here and is not explained by the paper.

### GRANADA source

An audit of the source (Kopra/source_10.0.0/modules) found no physics error in the O2 path:
- **Solar term and absorption rate.** The initial solar term and the absorption rate are correctly normalised. With the top-of-atmosphere beam they give 6.2e-9 s⁻¹ for the A band, the same as ours.
- **Rate laws.** These are `p1` or `p1·exp(p2/T)`, with correct detailed balance for the reverse rates.

**Where attenuation happens.** The solar beam is attenuated line by line (Voigt lines, slant paths) only inside a radiative-transfer sub-calculation. It stays at the unattenuated top-of-atmosphere value at every altitude in three cases:
- when the b-X band is not part of such a calculation;
- for calculation types 11, 13, 21, 31 and 41, which never compute the solar geometry;
- below a sub-calculation's lower altitude limit, where it is frozen at its value just above that limit.

**Not checkable.** The ERS setup, process and spectroscopy files are not available. So it is open which route the ERS runs took, and whether a configuration value (degeneracies, band A, line reduction factor) supplies the remaining factor of about 2. GRANADA's code also does not stop the O(1D) + O2 yields from summing above 1.

The flat, roughly 1.2e-8 s⁻¹ pumping inferred above fits the unattenuated routes. The ERS O2(b) is therefore not a reference for daytime O2(b) below about 80 km.

### Update: JPL yield and wavelength-dependent O2(a, v)

With the JPL O2(b) yield of 0.8 and the wavelength-dependent O2(a, v) split, april+00 gives the following against GRANADA (sasktran2 / GRANADA):

| Altitude | O(1D) | O2(a) | O2(b) |
|---|---|---|---|
| 40-60 km | 0.97-0.99 | 0.83-0.89 | 0.44-0.47 |
| 70-80 km | 0.94-0.99 | 0.70 | 0.50-0.55 |

O2(b) drops by about 18% at 40 km, where the O(1D) route dominates, and is unchanged at 80 km, where pumping dominates.

## Open items

- GRANADA's A-band pumping (about 1.2e-8 s⁻¹, nearly unattenuated to 30 km), about 2-10× ours.
- γ-band (2-0) and a-X (1.27 µm) emission constituents.
- O(1S) daytime sources.
