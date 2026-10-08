# sasktran2-nlte Mechanism Format

A mechanism is a TOML file. It lists the excited states to solve for, the background species whose densities are inputs, and the processes that couple them. The `sasktran2-nlte` crate validates the file on load and solves for steady-state populations. Bundled mechanisms live in `rust/sasktran2-nlte/mechanisms/`.

```python
import sasktran2 as sk

mechanism = sk.nlte.Mechanism.bundled("oxygen_yankovsky")   # or Mechanism.from_file(path)
solution = sk.nlte.solve(mechanism, atmosphere, rates)
sk.nlte.budget(mechanism, solution, "O2(b)")                # production and loss by process
```

## Model

Every process has an event rate:

$$R = c \prod_{\text{reactants}} n_i.$$

- **The coefficient $c$** is one of:
  - a temperature-dependent rate coefficient (reactions);
  - a supplied per-molecule rate (photolysis);
  - an Einstein A coefficient (radiative transitions).
- **Reactants.** Each reactant occurrence is consumed once per event.
- **Products.** A channel with yield $y$ produces its products at $yR$.
- **Steady state.** For each state, production equals loss.
- **Linearity.** A process may have at most one state reactant, so the system is linear in the states and each level needs one linear solve. Loading rejects anything else.

Units are SI throughout:

| Quantity | Units |
|---|---|
| Number densities | m⁻³ |
| Supplied rates and Einstein A | s⁻¹ |
| Rate coefficients | converted on load from the units given in the file |

## File layout

The numbers in this and later examples are illustrative, not vetted data.

```toml
[mechanism]
name = "example"
version = "0.1.0"
description = "..."
background = ["O2", "O3", "N2", "O(3P)"]   # densities supplied as inputs [m^-3]

[references]
jpl19 = "Burkholder et al. (2019), JPL Publication 19-5"
nist = "NIST Atomic Spectra Database"

[[state]]
id = "O(1D)"

[[state]]
id = "O2(a)"

[[reaction]]
id = "o1d_n2"
reactants = ["O(1D)", "N2"]
products = ["O(3P)", "N2"]
rate = { law = "arrhenius", a = 2.15e-11, ea_over_r = -110.0, units = "cm3 s-1" }
reference = "jpl19"

[[photo_reaction]]
id = "o3_hartley"
reactant = "O3"
products = ["O(1D)", "O2(a)"]
rate_input = "J_O3_O1D"            # supplied per altitude [s^-1], quantum yield included
reference = "jpl19"

[[transition]]
id = "red_line"
upper = "O(1D)"
lower = "O(3P)"
einstein_a_s = 5.6e-3
wavelength_nm = 630.0              # optional
reference = "nist"
```

Rules that apply to every entry:

- Every process needs a `reference` key that exists in `[references]`.
- Process ids must be unique.
- Unknown keys are an error. This catches typos such as `wavelenght_nm`.

## Species ids

- **Ids are free-form strings.** Each must be declared, either as a `[[state]]` or in `mechanism.background`.
- **`M` is special.** When `M` is a background species and no density is supplied, it is derived from `pressure_pa` and temperature.
- **Ground vibrational level.** A `v=0` qualifier is dropped on load, so `O2(b, v=0)` and `O2(b)` name the same species. Write either one.

## Reactions

**Products.** Give exactly one of:
- `products`, a single channel with yield 1; or
- `channels`, a list of branches sharing one total rate:

```toml
[[reaction]]
id = "o1d_o2"
reactants = ["O(1D)", "O2"]
rate = { law = "arrhenius", a = 3.2e-11, ea_over_r = -67.0, units = "cm3 s-1" }
channels = [
  { yield = 0.40, products = ["O2(b, v=1)", "O(3P)"] },
  { yield = 0.55, products = ["O2(b)", "O(3P)"] },
  { yield = 0.05, products = ["O2(a)", "O(3P)"] },
]
reference = "..."
```

**Yields.**
- Each yield is between 0 and 1, and the yields of a reaction sum to at most 1.
- If they sum to less than 1, the remaining events go to untracked products. The loss still uses the full rate.
- A product listed twice in a channel is produced twice per event.

**Rate laws.** Two laws are available:

| `law` | Parameters | $k(T)$ |
|---|---|---|
| `constant` | `value` | `value` |
| `arrhenius` | `a`, `n` (default 0), `t0` (default 300 K), `ea_over_r` (default 0) | $a\,(T/t_0)^n \exp(-E_a/RT)$ |

**Units.** `units` is required. It must match the number of reactants:

| Reactants | Accepted units |
|---|---|
| 1 | `s-1` |
| 2 | `cm3 s-1` or `m3 s-1` |
| 3 | `cm6 s-1` or `m6 s-1` |

## Families

Many processes repeat over vibrational levels. Any `[[state]]`, `[[reaction]]`, `[[photo_reaction]]` or `[[transition]]` can take `for_v = [lo, hi]`, which expands the entry once for each `v` in the inclusive range.

- **Placeholders.** `{v}`, `{v-1}`, `{v+2}` and so on can appear in any string field. Using one without `for_v` is an error, as is a level below zero.
- **`exp_v`.** In a family, a rate may also take `exp_v`, which multiplies it by $\exp(\texttt{exp\_v}\cdot v)$.
- **Ids.** Give each family a templated id so the expanded ids stay unique.

```toml
[[state]]
id = "O2(X, v={v})"
for_v = [1, 35]

[[reaction]]
id = "o2x{v}_o2_vv"
for_v = [4, 20]
reactants = ["O2(X, v={v})", "O2"]
products = ["O2(X, v={v-1})", "O2(X, v=1)"]
rate = { law = "constant", value = 1.3e-12, exp_v = -0.31, units = "cm3 s-1" }
reference = "..."
```

## Solver inputs and outputs

### Inputs to `sasktran2.nlte.solve(mechanism, atmosphere, rates)`

| Argument | Contents |
|---|---|
| `atmosphere` | `temperature_k` on a single dimension (usually `altitude`), and one variable per background species named by its id, in m⁻³. Optionally `pressure_pa`, used only to derive `M`. |
| `rates` | One profile per name in `mechanism.rate_inputs`, in s⁻¹. |

Missing inputs raise `ValueError`, listing every missing name.

### Output Dataset

| Variable | Dimensions | Units |
|---|---|---|
| `density` | state, altitude | m⁻³ |
| `production`, `loss` | state, altitude | m⁻³ s⁻¹ |
| `process_rate` | process, altitude | m⁻³ s⁻¹ (radiative processes: photons) |
| `photon_ver` | transition, altitude | photons m⁻³ s⁻¹ |
| `relative_residual` | altitude | max \|production − loss\| / max production |

The process coordinates carry each process's `process_kind` and `process_reference`. The transition coordinates carry `transition_upper`, `transition_lower` and `transition_wavelength_nm`.

## Not yet supported

- **Nonlinear systems.** Processes with two state reactants, such as energy pooling, are rejected on load.
- **Time dependence.** There are no transient solutions; O2(a¹Δ) at twilight needs one. This is planned.
