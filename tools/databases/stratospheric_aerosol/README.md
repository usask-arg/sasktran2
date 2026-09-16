# Building the stratospheric aerosol catalogue

The runtime API uses the small packaged `stratospheric_aerosol_v1.nc`. This
offline builder is needed only to reproduce or revise that catalogue.

```bash
python tools/databases/stratospheric_aerosol/build.py \
  --psd-root /path/to/SAGEIII_ISS/aerosol/particle_size/v1.0.0 \
  --sage-root /path/to/SAGE_III_ISS/monthly_v6-0 \
  --output src/sasktran2/_data/stratospheric_aerosol
```

Use the repository's Python environment. The spectral consistency check uses
the same H2SO4 refractive index, width 1.6, and radius/wavelength Mie table as
the USask retrieval. SASKTRAN2 may generate that table or download its standard
refractive-index input on first use. The builder reads 110 PSD monthly files
from June 2017 through July 2026 and their matching original SAGE solar files.
It matches events by ID and verifies that copied extinction agrees exactly.

`build_report_v1.json` records source filenames, SHA-256 digests, selection
populations, event identities, calibration monthly rates and sensitivity
results. Reproduction requires those exact source revisions. Source updates
require scientific review and a new derived catalogue version; they must not
silently replace an existing runtime catalogue. The builder chooses the latest
available original revision per month and fails on extinction mismatches.

Selection uses 756 nm AOD over a common 18–30 km interval. For every level in
that interval, require positive finite extinction; finite median radius strictly
inside 10.01–589.99 nm; nonnegative formal radius/extinction uncertainties below
50%; source aerosol flags 2 or 3 at 756/869/1021/1543 nm; and at least 1 km above
the source aerosol tropopause. Three known anomalous March 9, 2024 events are
excluded explicitly. The median absolute fractional mismatch of the three
near-IR extinction ratios with the retrieved sulfate size model must be below
30%. This is a consistency screen, not independent validation of the retrieval.

Within each latitude band, low/typical/elevated use ±3 percentile points around
P10/P50/P90; extreme uses P98.5–P99.5. Choose the actual paired profile nearest
the candidate group's median log extinction and median radius, with distance
scales 0.5 in log extinction and 50 nm in radius. Retain the selected profile's
contiguous valid interval around the common core, extending through reliable
source measurements on either side. No internal gaps are bridged. All nine
measured spectra and formal errors within that interval remain in the catalogue.

Upper-tail calibration uses the screened P10–P60 population and requires flag
2 throughout 26–30 km in all four near-IR channels. The default uses raw
log-extinction Theil–Sen slopes over 27–30 km, monthly median rates with at
least five profiles, then an unweighted median of monthly rates. Individual
non-decreasing slopes are retained when aggregating, avoiding a sign-selection
bias in noisy data. The aggregate must decrease and have at least 24 contributing
months. Convert rate to scale height and round to 100 m. The report also records
raw fits over 26–29 and 26–30 km and the former 1.5 km-smoothed 27–30 km fit.

Radius units are verified from the USask retrieval implementation because the
PSD files omit them. Extinction is converted from km^-1 to m^-1; altitude from
km to m; radius remains in nm. Do not use the copied source descriptions to infer
units or independent information content: the 1 km radius retrieval is stored
on a 0.5 km grid, and formal uncertainties omit model/systematic uncertainty.

After a reviewed rebuild, update the pinned digest in
`src/sasktran2/database/stratospheric_aerosol.py`, run the climatology and
scatterer tests, execute the documentation example, and check representative
spectra/radiances including the 3.2 km and 4 km upper-tail overrides. The example
documents literature comparisons and the limitations of modeled extensions.
