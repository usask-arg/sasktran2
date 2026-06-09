# Photochemistry Emission Notes

## Implemented state

- Photochemical model output can now be converted to line-list volume emission through `sasktran2.constituent.PopulationEmissionRate`.
- `Yankovsky.solve(...)` remains the population-producing photochemical model; O2 band emission construction is now delegated to the standalone Rust-backed population-emission constituent.
- The oxygen green line uses the `O(1S) -> O(1D)` transition at `557.7 nm` with the scaffolded Einstein A coefficient in Rust.
- `Yankovsky.oxygen_green_line_mcdade()` computes a direct green-line VER using the McDade/Murtagh short-form Barth parameterization from `T`, `[O]`, `[O2]`, and `[N2]`.
- `Yankovsky.oxygen_green_line_mcdade_constituent()` converts that parameterized VER into a monochromatic volume-emission constituent.
- The oxygen A-band uses HITRAN O2 line metadata for resolved `b 0 -> X 0` and `b 1 -> X 1` transitions between `759 nm` and `776 nm`.
- The oxygen B-band uses HITRAN O2 line metadata for resolved `b 1 -> X 0` transitions between `680 nm` and `700 nm`.
- The Rust line-list volume emission constituent maps population-specific photon VER profiles onto Doppler-broadened emission lines and into the SASKTRAN2 spectral grid.
- O2 line weights can be evaluated with normalized Einstein-A branching or a HITRAN line-strength LTE fallback.
- `Yankovsky.emission_constituents()` and `Yankovsky.add_emissions_to_atmosphere()` remain thin compatibility wrappers; O2 band emission itself is no longer calculated in Python.

## Data sources

- Einstein A values, transition wavenumbers, lower-state energies, quantum labels, and statistical weights come from the local HITRAN O2 line database through `HITRANLineDatabase`.
- The line-list loader preserves those metadata fields in Rust `OpticalLine` records so emission code does not need to parse raw HITRAN text.
- Current A/B-band filtering uses available HITRAN global/local quantum labels and wavelength bounds to select `O2(b, v=0) -> O2(X, v=0)`, `O2(b, v=1) -> O2(X, v=1)`, and `O2(b, v=1) -> O2(X, v=0)` transitions.
- Branching ratios are normalized within each resolved upper-state identifier. HITRAN O2 often lacks enough upper local quantum detail, so the current implementation falls back to a resolved-level surrogate when needed.
- The McDade green-line parameterization follows the short-form equation reproduced by Murtagh et al. (1990) and summarized by Lednyts'kyy and von Savigny (2020): `C0 = 0`, `C1 = 211`, `C2 = 15`, `A558 = 1.18 s^-1`, `A1S = 1.35 s^-1`, `k1 = 4.7e-33 * (300 / T)^2 cm^6 s^-1`, and `3k5 = 4e-12 * exp(-865 / T) cm^3 s^-1`.
- Public inputs to the McDade implementation use `m^-3`; the rate expression is evaluated internally in `cm^-3` and converted back to `photons m^-3 s^-1`.

## Current limitations

- The O2 band source still assumes vibrational populations supplied by the photochemical model with rotational LTE within each vibrational state; full NLTE rotational populations are not solved yet.
- The `O2(b, v=1) -> O2(X, v=0)` B-band branch currently uses a scaffolded branch A-value and needs a sourced value before absolute B-band intensity is final.
- LTE line weights are a placeholder for testing spectral redistribution and radiance plumbing. They should be replaced or complemented by rotationally resolved NLTE state populations later.
- The green-line `O(1S)` state is still scaffolded in the full photochemical state solve; the McDade pathway currently provides a direct empirical VER rather than a coupled O(1S) source/loss state.
- Altitude-dependent line-weight derivatives are currently supported only when the source altitude grid matches the native model altitude grid.

## Next steps

- Compare McDade green-line limb profiles against a reference atmosphere or published nightglow profile.
- Add physically sourced `O(1S)` production and quenching pathways for the green line when the full coupled chemistry is expanded.
- Replace the O2 band LTE rotational placeholder with an NLTE population provider that produces upper-state populations on the emission altitude grid.
- Tighten HITRAN state identification once the NLTE state model has an explicit rotational-vibrational indexing convention.
- Add validated reference calculations for green-line and A-band limb radiances once the science inputs are fixed.
