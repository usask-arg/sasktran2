# Resident OMPS memory validation

`omps_memory_probe.py` runs the saved OMPS orbital-plane ozone objective and
gradient with all time groups resident. `omps_memory_compare.py` reads the saved
arrays without running radiative transfer. These tools require the external
`data-analysis` and `arg-processing` repositories and OMPS inputs; they are not
a self-contained benchmark shipped with SASKTRAN2.

## Prerequisites

- Use the tomography processor's pinned virtual-environment Python. Keep its
  executable path rather than resolving its symlink to the base interpreter.
- The saved guard summary supplies the exact reference command, input paths and
  Python executable. The default summary is
  `data-analysis/outputs/omps_paper_5113/investigation/memory_logs/source_memory_ozone_ms5_in110_out110.summary.json`.
- The analysis repository supplies the physical-footprint guard and shared
  calculation lock. Run only one calculation at a time. The probe uses one RT
  thread and a 27 GiB ceiling (28.991 decimal GB), including compressed memory;
  RSS alone is not the acceptance measurement.
- Build and extract a candidate wheel separately, preserving the pinned
  installation. Match the reference compiler and numerical dependencies when
  checking numerical equivalence. Record the build's source hashes, compiler,
  dependencies, wheel hash and native binary hash.
- Candidate extraction directories contain `sasktran2/`. Supplying both
  `--wheel` and `--build-provenance` checks the imported native binary and every
  extracted package file against the supplied evidence before launching RT.
  A build record must include `build_completed_utc`, `binary_sha256` and
  `wheel_sha256`; these fields supplement the recorded build configuration.

## Run and compare

Run from the SASKTRAN2 checkout with `OMPS_PYTHON` set to the pinned executable.
Prepare the candidate package, wheel and build record at the example paths
below. The output directory must be fresh and beneath `build/omps-memory`.

```sh
"$OMPS_PYTHON" tools/benchmarks/omps_memory_probe.py \
  --candidate-package build/omps-memory/candidate-package \
  --wheel build/omps-memory/candidate.whl \
  --build-provenance build/omps-memory/build-provenance.json \
  --output build/omps-memory/candidate-ms11 --label candidate-ms11 \
  --columns 11 --check-updates --profile-geometry
```

For a different repository layout, supply both `--analysis-root` and
`--reference-summary`; changing the analysis root does not relocate the default
summary. Omitting the candidate, wheel and build options imports the pinned
reference package. Its installed binary has no commit attestation, so compare
its recorded runtime hash rather than assuming a source commit.

The probe keeps 110 incoming and outgoing directions and all time groups
resident. The default saved scene has 11 wavelengths, 25 vertical source points
and 26 time groups. `--columns` selects the horizontal source count.
`--check-updates` adds perturbed and
restored objective/gradient evaluations and initial/restored native JVP/VJP
products. `--profile-geometry` enables native memory diagnostic records.

Compare against a completed reference probe at the same resolution:

```sh
"$OMPS_PYTHON" tools/benchmarks/omps_memory_compare.py \
  --reference build/omps-memory/reference-ms11/reference-ms11 \
  --trial build/omps-memory/candidate-ms11/candidate-ms11 \
  --reference-summary build/omps-memory/reference-ms11/guard.summary.json \
  --trial-summary build/omps-memory/candidate-ms11/guard.summary.json \
  --require-updates --output build/omps-memory/comparison.json
```

`--require-updates` requires both complete histories. `--initial-only` instead
compares the saved initial evaluation and cannot certify whole-probe completion.
Different angular or spatial resolutions are rejected even when arrays are
close. Reference equivalence is reported separately from each probe's internal
adjoint and restored-state checks; an inherited failing check remains visible.
The default comparison tolerances are `rtol=1e-10` and `atol=1e-12`. The exit
status uses these tolerances; the JSON also reports bitwise equality for every
array.

## Validated checkpoint

The native implementation at `60de1b4e` was validated on OMPS orbit 50430
(21 July 2021), with 159 images, 10,176 rays and 15,390 ozone parameters.
Both resident probes completed under the guard at the resolution above.

| Horizontal source columns | Peak physical footprint, decimal GB | Full validation probe, seconds |
| --- | ---: | ---: |
| 5 | 9.747 | 308.06 |
| 11 | 18.245 | 638.18 |

All 17 corresponding saved arrays were bitwise identical to the preceding
validated checkpoint at each matching column count, including complete
radiances, ozone gradients, atmosphere updates and native products. Against
the pinned reference at five columns, radiances and ozone gradients were
bitwise identical; two arbitrary VJP arrays differed by at most `7.11e-15`.

A fresh unchanged 11-column control used 18.427 GB and 636.29 seconds. The
candidate's runtime was 0.297% longer in this pair; individual probes do not
establish a statistical runtime bound or causal speedup. Existing strict
adjoint and warm-start restoration residuals were unchanged. These fixed-state
and update probes do not establish full iterative-retrieval capacity.

Validation included 276 targeted Python tests, 131 C++ cases with 101,539
assertions, the focused Rust derivative-storage test, Clippy and pre-commit.
Detailed arrays, guard logs and source/build attestation records remain in the
local ignored `build/omps-memory` archive; they are not bundled with the PR.
