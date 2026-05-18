# Lyman-alpha Photolysis Notes

## Implemented state

- `PhotoReaction` supports explicit monochromatic line metadata through `with_line_center_nm`.
- The Yankovsky O2 Lyman-alpha reaction is enabled at `121.567 nm`.
- Continuum photolysis still uses the existing wavelength-bin integration rule.
- Line photolysis interpolates actinic flux and cross section to the line center, then computes `J = q * F_line * sigma_line`.
- `src/sasktran2/photchem/__init__.py` includes `121.567 nm` exactly in the generated actinic-flux wavelength grid.
- The photolysis-rate calculation is implemented in Rust core code and called from the PyO3 wrapper.

## Scientific inputs to confirm

- Adopt the Lyman-alpha center wavelength, nominally 121.567 nm in vacuum.
- Confirm the O2 cross section and quantum yield at Lyman-alpha. The current implementation uses a quantum yield of `1.0`.
- Confirm whether attenuation should use the SASKTRAN2 actinic flux at 121.567 nm, a Beer-Lambert parameterization from TOA, or a hybrid fallback when high-resolution actinic flux is unavailable.
- Confirm whether `actinic_flux` at the line center should be interpreted as line-integrated photon flux or whether an explicit line-flux input should be added.

## Test coverage

- `PhotoReaction` line-center metadata.
- Yankovsky Lyman-alpha reaction configuration.
- Existing central-difference wavelength-bin widths.
- Continuum photolysis integration over a selected wavelength range.
- Negative actinic-flux and cross-section clamping.
- Line photolysis interpolation between grid points.
- Exact-grid-point line photolysis.
- Error handling for missing line wavelength coverage.
- Error handling for mismatched spectral array shapes.

## Remaining design question

The line treatment avoids dependence on the `0.001 nm` grid spacing, but the flux units need a final decision. If SASKTRAN2 `actinic_flux` remains a spectral quantity per nm, Lyman-alpha should probably receive a separate line-integrated flux input rather than using the interpolated spectral value directly.
