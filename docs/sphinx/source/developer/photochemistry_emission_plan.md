# Photochemistry Emission Notes

## Implemented state

- Photochemical model output can now be converted to photon volume emission rates through `Yankovsky.emissions()`.
- The oxygen green line uses the `O(1S) -> O(1D)` transition at `557.7 nm` with the scaffolded Einstein A coefficient in Rust.
- `Yankovsky.oxygen_green_line_mcdade()` computes a direct green-line VER using the McDade/Murtagh short-form Barth parameterization from `T`, `[O]`, `[O2]`, and `[N2]`.
- `Yankovsky.oxygen_green_line_mcdade_constituent()` converts that parameterized VER into a monochromatic volume-emission constituent.
- The oxygen A-band uses HITRAN O2 line metadata for resolved `b 0 -> X 0` transitions between `759 nm` and `772 nm`.
- A new line-list volume emission constituent maps a total photon VER profile onto many emission lines and into the SASKTRAN2 spectral grid.
- A-band line weights can be fixed from normalized Einstein A values or evaluated as altitude-dependent LTE placeholder weights using upper-state energies, statistical weights, and per-upper-state branching ratios.
- `Yankovsky.emission_constituents()` returns the green-line and A-band constituents without requiring callers to wire the source objects manually.
- `Yankovsky.add_emissions_to_atmosphere()` attaches those constituents directly to an `Atmosphere`; when A-band temperature weights are requested implicitly, atmosphere temperature is interpolated onto the emission altitude grid.

## Data sources

- Einstein A values, transition wavenumbers, lower-state energies, quantum labels, and statistical weights come from the local HITRAN O2 line database through `HITRANLineDatabase`.
- The line-list loader preserves those metadata fields in Rust `OpticalLine` records so emission code does not need to parse raw HITRAN text.
- Current A-band filtering uses available HITRAN global/local quantum labels and wavelength bounds to select the `O2(b, v=0) -> O2(X, v=0)` band.
- Branching ratios are normalized within each resolved upper-state identifier. HITRAN O2 often lacks enough upper local quantum detail, so the current implementation falls back to a resolved-level surrogate when needed.
- The McDade green-line parameterization follows the short-form equation reproduced by Murtagh et al. (1990) and summarized by Lednyts'kyy and von Savigny (2020): `C0 = 0`, `C1 = 211`, `C2 = 15`, `A558 = 1.18 s^-1`, `A1S = 1.35 s^-1`, `k1 = 4.7e-33 * (300 / T)^2 cm^6 s^-1`, and `3k5 = 4e-12 * exp(-865 / T) cm^3 s^-1`.
- Public inputs to the McDade implementation use `m^-3`; the rate expression is evaluated internally in `cm^-3` and converted back to `photons m^-3 s^-1`.

## Current limitations

- The A-band source still assumes a total `O2(b)` population profile supplied by the photochemical model; full NLTE rotational-vibrational populations are not solved yet.
- LTE line weights are a placeholder for testing spectral redistribution and radiance plumbing. They should be replaced or complemented by NLTE state populations later.
- The green-line `O(1S)` state is still scaffolded in the full photochemical state solve; the McDade pathway currently provides a direct empirical VER rather than a coupled O(1S) source/loss state.
- Altitude-dependent line-weight derivatives are currently supported only when the source altitude grid matches the native model altitude grid.

## Next steps

- Compare McDade green-line limb profiles against a reference atmosphere or published nightglow profile.
- Add physically sourced `O(1S)` production and quenching pathways for the green line when the full coupled chemistry is expanded.
- Replace the A-band LTE placeholder with an NLTE population provider that produces upper-state populations on the emission altitude grid.
- Tighten HITRAN state identification once the NLTE state model has an explicit rotational-vibrational indexing convention.
- Add validated reference calculations for green-line and A-band limb radiances once the science inputs are fixed.
