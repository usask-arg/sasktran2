# Lyman-alpha Photolysis Notes

## Implemented state

- `PhotoReaction` supports explicit monochromatic line metadata through `with_line_center_nm`.
- The Yankovsky O2 Lyman-alpha reaction is enabled at `121.567 nm` with an effective cross section.
- Continuum photolysis still uses the existing wavelength-bin integration rule.
- Line photolysis uses line-integrated actinic flux and computes `J = q * Phi_line * sigma_eff`.
- `src/sasktran2/photchem/__init__.py` includes `121.567 nm` exactly in the generated actinic-flux wavelength grid.
- `src/sasktran2/photchem/__init__.py` also exports `lyman_alpha_actinic_flux`, a line-integrated flux channel scaled from the modeled actinic factor at 121.567 nm.
- The photolysis-rate calculation is implemented in Rust core code and called from the PyO3 wrapper.

## Scientific inputs to confirm

- Adopt the Lyman-alpha center wavelength, nominally 121.567 nm in vacuum.
- Confirm the O2 effective cross section and quantum yield at Lyman-alpha. The current implementation uses a quantum yield of `1.0`.
- Confirm whether attenuation should use the SASKTRAN2 actinic flux at 121.567 nm, a Beer-Lambert parameterization from TOA, or a hybrid fallback when high-resolution actinic flux is unavailable.
- Confirm the adopted line-integrated TOA flux. The current implementation uses `3.2e15 photons m^-2 s^-1`, which gives `sigma_eff = 1.0625e-24 m^2` from the existing `3.40e-9 s^-1` TOA rate.

## Test coverage

- `PhotoReaction` line-center metadata.
- Yankovsky Lyman-alpha reaction configuration.
- Existing central-difference wavelength-bin widths.
- Continuum photolysis integration over a selected wavelength range.
- Negative actinic-flux and cross-section clamping.
- Line photolysis interpolation between grid points.
- Exact-grid-point line photolysis.
- Effective-cross-section line photolysis.
- Lyman-alpha TOA rate and effective-cross-section consistency.
- Error handling for missing line wavelength coverage.
- Error handling for mismatched spectral array shapes.

## Remaining design question

The line treatment avoids dependence on the `0.001 nm` grid spacing. The remaining science decision is whether `3.2e15 photons m^-2 s^-1` is the right nominal TOA line-integrated flux for the intended reference conditions, or whether this should be supplied as a time-dependent solar input.
