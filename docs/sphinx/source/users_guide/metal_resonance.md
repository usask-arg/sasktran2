(_users_metal_resonance)=
# Metal resonance scattering

Metal atoms, ions, and diatomic molecules absorb sunlight in narrow spectral lines.
`LineResonance` calculates temperature-dependent line cross sections and the
polarized phase matrix for elastic return to the absorbing lower state. Separate
species classes share this implementation and use species-specific spectroscopy.
The supplied spectroscopy is selected for ordinary mesospheric metal layers;
it is not a hot-plume, meteor-trail, or artificial-release model.
The [candidate survey](metal_candidates.md) separates established atmospheric
targets, speculative candidates, and missing spectroscopy.
Chemical nightglow can be fitted independently using the
[free VER emission constituents](metal_emission.md), without an excitation model
or an absolute absorption cross section.

## Add a metal layer

Use a resonance optical property with
{py:class}`sasktran2.constituent.NumberDensityScatterer`. Number densities are in
m⁻³ and refer to the total population of the named atomic ionization stage or
molecular species. Lower-state populations follow a Boltzmann distribution at the
atmospheric temperature; the optical property does not calculate ionization or
chemical equilibrium.

When temperature comes from `atmosphere.temperature_k`, the metal contribution
is included in the shared `wf_temperature_k` on the atmosphere's state grid.
Passing an explicit `temperature_k` profile to `NumberDensityScatterer` instead
produces `wf_<name>_temperature_k` on that constituent's altitude grid. Both
include the temperature dependence of the line shape and level populations.

```python
import numpy as np
import sasktran2 as sk

sodium = sk.optical.Sodium()
altitude_m = np.arange(70_000.0, 115_001.0, 1_000.0)
number_density = 3e9 * np.exp(-0.5 * ((altitude_m - 92_000.0) / 5_000.0) ** 2)

# `atmosphere` already has geometry, vacuum wavelengths, and temperature_k.
atmosphere["sodium"] = sk.constituent.NumberDensityScatterer(
    sodium, altitude_m, number_density
)
```

A local NetCDF file can be supplied explicitly as `Sodium(db_filepath=path)`.
For another species use `AtomicResonance(species, db_filepath=path)` or
`MolecularResonance(species, db_filepath=path)`. An in-memory `xarray.Dataset` can
be passed with `db=dataset`, including to the common `LineResonance` base.
Missing default files are downloaded through `StandardDatabase` and cached under
`spectroscopy/metals/`. For local-only access, use
`database = sk.database.MetalSpectroscopyDatabase(db_root=...)` and pass
`db_filepath=database.path("Na_I")` to `Sodium`. The optical property's `db`
argument accepts an in-memory dataset, not a database manager.
Molecular datasets represent one named isotopologue; isotope-abundance scaling
is not applied automatically.

## Evaluate optical properties directly

```python
wavelength_nm = np.linspace(589.156, 589.160, 2001)
properties = sodium.cross_sections(
    wavelength_nm,
    np.array([90_000.0]),
    temperature_k=np.array([200.0]),
    num_stokes=3,
    num_moments=3,
)

extinction_m2 = properties.extinction
elastic_scattering_m2 = properties.extinction * properties.ssa
w2 = properties.polarizability
```

Cross sections have shape `(altitude, wavelength)` and units m² per atom or
molecule. The phase coefficients use the standard SASKTRAN2 stacked convention.
At least three Legendre orders are required to represent the electric-dipole
phase function. Scalar and three-Stokes calculations are supported.
The internal `atmosphere_quantities` interface follows the engine convention:
its field named `ssa` holds the scattering cross section, whereas the direct
`cross_sections` interface shown above returns a dimensionless albedo.

All wavelengths are **vacuum wavelengths in nm**. Familiar tabulated air
wavelengths such as the sodium lines near 588.995 and 589.592 nm must be converted
before using the optical property.

## Line profile and branching

For line $l\rightarrow u$, the frequency-integrated absorption cross section is

$$
S_{lu}=\frac{e^2}{4\epsilon_0 m_e c}f_{lu}.
$$

The monochromatic cross section is the lower-state population fraction times
$S_{lu}$ and a normalized Voigt profile. Doppler broadening depends on temperature
and species mass; the natural linewidth uses the total radiative upper-state
and, when supplied, lower-state decay rates. Broadening by collisions, bulk
velocity shifts, saturation, and
magnetic-field effects are not included.

By default, each profile is truncated at 0.1 nm on either side of its line center
without renormalization. Set `line_wing_cutoff_nm=None` to retain complete Voigt
wings. This cutoff controls evaluation cost; it is not an instrument response.

Only the branch returning to the original lower state contributes to elastic
scattering:

$$
\sigma_{\rm elastic}=\sum_{lu}\sigma_{lu}\frac{A_{ul}}{\Gamma_u},
\qquad
\omega=\frac{\sigma_{\rm elastic}}{\sum_{lu}\sigma_{lu}}.
$$

Here $\Gamma_u$ includes radiative decays to all available lower states, including
lines outside the evaluated wavelength interval. A return branching fraction
below one represents removal from the current wavelength. It does not imply that
the missing photons are destroyed: they can appear as fluorescence at another
wavelength. Incomplete atomic decay data are identified in the source data and
produce a warning because missing decay channels can overestimate the elastic
return fraction.

For an isolated line with unpolarized lower magnetic substates, the phase matrix
for Stokes I, Q, U is

$$
\mathbf P=(1-W_2)\mathbf E_{11}+W_2\mathbf P_{\rm Rayleigh}.
$$

The angular-momentum quantum numbers determine $W_2$. A $J_l=0\rightarrow J_u=1$
transition gives the classical Rayleigh matrix. The fine-structure sodium D1 and
D2 transitions have $W_2=0$ and $1/2$, respectively. The scalar phase function is

$$
P_{11}(\mu)=1+\frac{W_2}{4}(3\mu^2-1),\qquad\mu=\cos\theta,
$$

normalized so that its integral over solid angle is $4\pi$. Overlapping lines are
combined using their elastic scattering cross sections as weights. These
relations follow the resonance-scattering treatment in
[Langowski et al. (2016), Eq. 4 and Table 1](https://amt.copernicus.org/articles/9/295/2016/).

The same angular-momentum treatment applies to isolated diatomic rotational
lines. It gives different polarizabilities for P, Q, and R branches, rather than
assigning one Rayleigh fraction to an entire molecular band. See
[Berdyugina, Stenflo and Gandorfer (2002), Eqs. 15–17](https://doi.org/10.1051/0004-6361:20020587).

## Limits for retrievals

The optical property represents the **elastic return component**. It does not
transfer photons between wavelengths or include emission pumped through other
lines. This distinction is especially important for AlO and other molecules,
where a common upper state can decay to many rotational and vibrational levels.
The modeled elastic component alone is not the total molecular fluorescence
spectrum. Fluorescence and Raman branches require a source calculation that
connects incident and emitted wavelengths.
The molecular data describe selected published band systems, not guaranteed
complete UV–visible opacity. For example, AlO ATP represents the X, A and B
electronic states; additional higher electronic systems are not supplied.

The atomic phase matrices use isolated fine-structure transitions. Hyperfine
structure, isotope shifts, interference between upper levels, and ground-state
polarization are not resolved. These effects can change narrow line profiles and
polarization even if an instrument cannot resolve the individual components.
For the general treatment see
[Stenflo (1997)](https://ethz.ch/content/dam/ethz/special-interest/phys/particle-physics/cosmologygroup-dam/People/StenfloPDFs/AA324_344_1997.pdf).

Thermal motion also redistributes photon frequency during scattering in the
observer frame. Independent monochromatic calculations approximate this process;
they do not solve a partial-frequency-redistribution problem.

Resolve the intrinsic metal lines and the incident solar Fraunhofer spectrum
before applying the instrument spectral response. Evaluating a narrow line only
at a coarse instrument wavelength grid can miss nearly all of its integrated
signal. Convolving the optical cross section first generally does not reproduce
line transfer and self-absorption. The retained spectroscopy identifies candidate
signals; wavelength coverage alone does not establish OSIRIS detectability.

For density retrievals, the engine's scattering normalization has a singular
derivative at locations with exactly zero total scattering. The calculation now
raises an explicit error if a fitted density can add scattering there. Include
the physical ambient Rayleigh background or start from a positive metal density
at fitted locations. Forward calculations with `calculate_derivatives=False`,
zero-cross-section wavelengths, and locations outside a fitted profile remain
supported. This limitation does not apply to the independently retrieved VER
source: its derivatives remain defined at zero emission.

The temperature Jacobians include both the native-grid optical response and,
for extinction-normalized profiles, the source-grid conversion to number
density. Single-scattering finite-difference checks cover different profile
grids, overlapping lines, and mixed scattering. The existing discrete-ordinates
layer interpolation approximates phase-function Jacobians; multiple-scattering
retrievals still require grid-convergence and finite-difference checks for the
chosen atmosphere.

## Local spectroscopy format

Each NetCDF file describes a single species or molecular isotopologue. Required
line arrays use the dimension `line`:

| Variable | Meaning |
| --- | --- |
| `wavelength_nm` | Positive vacuum line-center wavelength, nm |
| `oscillator_strength` | Positive absorption oscillator strength $f_{lu}$ |
| `lower_energy_cminv` | Lower-state energy above the ground state, cm⁻¹ |
| `lower_j`, `upper_j` | Total electronic or rotational angular momentum |
| `einstein_a_s` | Spontaneous return rate $A_{ul}$, s⁻¹ |
| `upper_total_a_s` | Total spontaneous decay rate from the upper state, s⁻¹ |

Optional line arrays are `lower_total_a_s`, default zero, and
`lower_statistical_weight`, default $2J_l+1$. ExoMol nuclear-spin statistical
weights must use the same convention in both the line weight and partition
function. `upper_decay_data_complete` records whether all upper-state radiative
decay channels are represented.

The dataset must also provide either one-dimensional `energy_cminv` and
`statistical_weight` state arrays, or a table of `partition_temperature_k` and
`partition_function`. State arrays include nonabsorbing states needed to
normalize the population. Tabulated partition functions are interpolated in
logarithmic temperature and partition function and cannot be extrapolated.

Required metadata are a positive `mass_amu` attribute and vacuum wavelength
convention (`wavelength_medium="vacuum"`). Keep species, isotopologue, source
URLs, source version, citations, selection limits, and redistribution permissions
in dataset metadata. Optional `temperature_min_k` and `temperature_max_k`
attributes constrain evaluation to the dataset's intended temperature range.
Likewise, `wavelength_min_nm` and `wavelength_max_nm` constrain evaluation to
the supplied spectral interval; values outside it raise an error rather than
assuming missing spectroscopy has zero opacity.

The prepared mesospheric subsets retain lines with lower-state energies up to
5000 cm⁻¹. Their upper temperature limits are 1000 K for atoms and ions and
500 K for molecules. Atomic partition sums contain classified states up to
5000 cm⁻¹; they are cold-atmosphere approximations. Molecular files retain the
full source partition-state data. Upper-state decay sums are formed before the
absorption-line selection, so out-of-band decay channels still affect return
probabilities.
