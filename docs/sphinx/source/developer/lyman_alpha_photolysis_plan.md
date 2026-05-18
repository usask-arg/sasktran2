# Lyman-alpha Photolysis Plan

## Current state

- `src/sasktran2/photchem/__init__.py` builds an actinic-flux dataset on a dense wavelength grid from 120 nm to 800 nm and exports species cross sections as `{species}_xs`.
- `rust/sasktran2-py-ext/src/photchem/yankovsky.rs` computes photolysis rates by integrating `quantum_yield * actinic_flux * cross_section * delta_wavelength` over each reaction wavelength range.
- `rust/sasktran2-rs/src/photchem/types.rs` already parses banded photon reactions such as `hv(lyman-alpha)`.
- `rust/sasktran2-rs/src/photchem/models.rs` already contains an O2 Lyman-alpha reaction, but it is effectively disabled with `with_quantum_yield(0.0)` and has no wavelength range or monochromatic handling.

## Scientific inputs to confirm

- Adopt the Lyman-alpha center wavelength, nominally 121.567 nm in vacuum.
- Confirm whether the model should use a monochromatic line flux, a finite line profile, or the nearest available spectral grid point.
- Confirm the O2 cross section and quantum yield at Lyman-alpha. The existing reaction has a TOA rate constant of `3.40e-9 s^-1`, but the current integration path does not use TOA rate constants.
- Confirm whether attenuation should use the SASKTRAN2 actinic flux at 121.567 nm, a Beer-Lambert parameterization from TOA, or a hybrid fallback when high-resolution actinic flux is unavailable.

## Implementation outline

1. Add first-class band metadata for monochromatic photolysis:
   - Extend `PhotoReaction` with an optional line center, or reuse `wavelength_range_nm` with a very narrow band only if the numerical behavior is acceptable.
   - Prefer an explicit `with_line_center_nm(121.567)` API so line reactions are not coupled to the wavelength grid spacing.

2. Update `Yankovsky::new()`:
   - Enable the O2 Lyman-alpha branch with the confirmed quantum yield.
   - Attach the Lyman-alpha line center or a validated finite integration range.
   - Keep the existing SRC range separate from the Lyman-alpha reaction to avoid double counting around 121.6 nm.

3. Update photolysis-rate calculation in `PyYankovsky::solve`:
   - For continuum reactions, keep the current wavelength-bin integration.
   - For monochromatic reactions, interpolate `actinic_flux` and `{species}_xs` to the line center for each altitude, then compute `J = q * F_line * sigma_line`.
   - If the chosen line representation is a finite profile, integrate the normalized profile against flux and cross section instead.
   - Validate input wavelength coverage and return a clear Python error if Lyman-alpha is requested but the dataset does not cover the line.

4. Update the actinic-flux helper:
   - Ensure the generated wavelength grid includes 121.567 nm exactly, or provide a small targeted grid around it for interpolation.
   - Verify `O2SchumannRunge` or the active O2 optical database provides valid cross sections at Lyman-alpha; add or switch databases if not.

5. Add tests:
   - Rust unit test for `PhotoReaction` line metadata and required-rate naming.
   - Rust/Python-extension test with a synthetic dataset where the expected Lyman-alpha `J` value is analytically known.
   - Python smoke test confirming the actinic-flux helper includes Lyman-alpha coverage and produces finite O2 cross sections there.
   - Regression test that continuum SRC integration is unchanged.

## Open design question

The main decision is whether Lyman-alpha should be treated as a true spectral line independent of grid spacing, or as a narrow wavelength interval on the existing actinic-flux grid. The line treatment is less fragile and avoids dependence on the current `0.001 nm` grid spacing, but it requires explicit interpolation and clearer units for line photon flux.
