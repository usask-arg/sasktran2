(_api_photochemistry)=
# Photolysis and Photochemistry

For executable examples, see {ref}`_users_photolysis` and {ref}`_users_airglow`.
The mechanism file format is described in the developer documentation,
`docs/sphinx/source/developer/nlte_mechanism_format.md`.

## Actinic flux and photolysis rates

```{eval-rst}
.. autosummary::
    :toctree: generated/

    sasktran2.photolysis.ActinicFlux
    sasktran2.photolysis.TUVActinicFlux
    sasktran2.photolysis.Photolysis
    sasktran2.photolysis.LinePhotolysis
    sasktran2.photolysis.LymanAlphaPhotolysis
    sasktran2.photolysis.TUVXQuantumYield
    sasktran2.photolysis.photolysis_rates
    sasktran2.photolysis.default_optical_properties
    sasktran2.photolysis.airglow_wavelength_grid
    sasktran2.photolysis.slant_columns
```

## Rate presets and quantum yields

```{eval-rst}
.. autosummary::
    :toctree: generated/

    sasktran2.photolysis.presets.oxygen_photolysis
    sasktran2.photolysis.presets.green_line_photolysis
    sasktran2.photolysis.presets.tuvx_v54_photolysis
    sasktran2.photolysis.quantum_yields.O3O1DYield
    sasktran2.photolysis.quantum_yields.o2_o1s_yield
    sasktran2.photolysis.quantum_yields.o3_o1d_matsumi2002
    sasktran2.photolysis.quantum_yields.o3_o3p_matsumi2002
```

## Parameterisations

```{eval-rst}
.. autosummary::
    :toctree: generated/

    sasktran2.photolysis.chabrillat_kockarts_cross_section
    sasktran2.photolysis.koppers_murtagh_cross_section
    sasktran2.photolysis.lyman_alpha_reduction_factor
```

## Excited-state kinetics

```{eval-rst}
.. autosummary::
    :toctree: generated/

    sasktran2.nlte.Mechanism
    sasktran2.nlte.solve
    sasktran2.nlte.budget
    sasktran2.nlte.add_photochemical_species
    sasktran2.nlte.emission_wavelength_grid
```
