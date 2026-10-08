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
| O3 → O(1D) + O2(a, v) | 6 inputs; yields 0.49, 0.15, 0.15, 0.08, 0.08, 0.05 (ym2006, 254 nm) | λ-dependent (yv2020 Eqs. 3-4); at 254 nm 0.362, 0.276, 0.115, 0.075, 0.076, 0.095 | Φ(O1D) only | ✓ yv2020 at 254 nm, on `J_O3_O1D`. Eq. 2 of yv2020 does not reproduce its Table 1, so the tabulated thresholds are used. O2(a, v≥1) relax to v=0 within milliseconds, so only the vibrational distribution depends on this. |
| O2 SRC → O(3P) + O(1D) | yield 1 | yield 1 | Φ = 1 at 139-175 nm | ✓ |
| O2 Lyman-α → O(1D) | Φ = 0.53 in the rate | 0.48-0.58 | 0.44 ± 0.05 at line centre | ✓ Φ = 0.53 in `J_O2_LYA` |
| O2 + hν → O2(b, v=0-2), O2(a) | line-by-line | top-of-atmosphere 5.35e-9, 2.94e-10, 7.94e-12, 1.54e-10 s⁻¹ (ym2006) | — | ✓ line-by-line; ours at the top of the atmosphere: 6.2e-9, 3.7e-10, 1.2e-11, 1.2e-10 |

### O(1D)

| Process | Legacy | ym2019 | JPL 19-5 | `oxygen` |
|---|---|---|---|---|
| + O2 | 3.2e-11 exp(67/T); b(v=1) 0.40, b(v=0) 0.55, a 0.05 | 3.3e-11 exp(55/T); b(v=1) 0.8, b(v=0) 0.2 (Pejakovic et al. 2014) | 3.3e-11 exp(55/T); total b 0.8 ± 0.2 | ✓ JPL rate, ym2019 yields. Open: JPL's total b yield would lower the O(1D) route by 20%. |
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

## Open items

- O(1D) + O2 total O2(b) yield: 1.0 here (ym2019) against 0.8 ± 0.2 (JPL).
- The GRANADA O2(b) factor of about 1.9; needs its input files or Funke et al. (2012).
- Wavelength-dependent O2(a, v) yields from O3 photolysis (yv2020, or GRANADA's).
- γ-band (2-0) and a-X (1.27 µm) emission constituents.
- O(1S) daytime sources.
