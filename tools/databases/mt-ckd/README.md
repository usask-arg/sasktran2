# MT_CKD 4.3 data generation

This directory records the extraction and validation recipe used to generate
SASKTRAN2's separately distributed MT_CKD 4.3 coefficient and reference files.
It contains no AER coefficient data, wrapper source, or generated artifacts.

The input MT_CKD checkout, its data, and the generated files are governed by
AER's license rather than SASKTRAN2's MIT license. Obtain MT_CKD directly from
AER-RC and review its terms before running this tool.

## Pinned inputs

- AER-RC/MT_CKD tag `4.3`, commit
  `1dad6e29363a2dc75d0eb642b7321b5db51079cc`
- Its LBLRTM submodule at
  `f66bb596e1686e5ff162a6fadc014e91d3ce01d4`
- An externally built shared library exposing the probe ABI described below.
  Its source is deliberately not redistributed in this MIT-licensed
  repository.

The published hashes were reproduced with GNU Fortran 14.2.0, NetCDF 4.9.2,
and NetCDF-Fortran 4.6.1 on macOS arm64. Other compatible toolchains may be
numerically identical without being byte-identical; the script treats a hash
mismatch as a failure so that such differences are explicit.

Obtain the pinned source with its submodules and separately prepare a compatible
probe library under the applicable upstream terms. Then run:

```shell
python tools/databases/mt-ckd/generate_mt_ckd_4_3.py \
    --mt-ckd-checkout /path/to/MT_CKD \
    --probe-library /path/to/libmtckd_probe.dylib \
    --output-directory /Volumes/SasktranFiles/sasktran2_db/v_latest/continuum
```

The script verifies the source revisions and refuses modified MT_CKD source,
MT_CKD data, or LBLRTM source. Through the external probe, it enables each
continuum component separately and reconstructs the 22 fields in the exact
order consumed by `rust/sasktran2-rs/src/continuum/mt_ckd.rs`. It then runs all
components, including the standalone model's historical Rayleigh term, for the
four validation atmospheres.

The probe library must export the C symbols `set_inputs`, `set_flags`,
`run_mtckd`, `get_absrb`, and `get_n2_data` with the signatures configured in
`generate_mt_ckd_4_3.py`. Supplying and building that adapter remains outside
this repository's licensing scope. The final hashes are the authoritative
compatibility check for both the adapter and compiler toolchain.

Successful output must have these hashes:

- `mt_ckd_4_3.bin`:
  `06b5600ecb0f5a3417c46d226555bc516f5a55ae93b19b6e7b5431e50f8f0459`
- `mt_ckd_4_3_reference.bin`:
  `0fdc39229ae020f2979511acfc81600ddfaa0b9bc0bbb03bb86a5ce3d01a9760`

The coefficient file begins with `MTCKD43\0`, followed by little-endian `u32`
values for the 1,991-point spectral length and 22-field count, then the fields
as field-major little-endian `f64` values. The reference file is four
consecutive 1,991-element little-endian `f64` spectra.
