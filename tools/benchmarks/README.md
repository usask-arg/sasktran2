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

# Successive-orders source interpolation

`so_interpolation_benchmark.py` compares the default successive-orders source
interpolation with the legacy behaviour
(`Config.successive_orders_legacy_interpolation = True`). The default uses
frame-aligned angular grids and cubic Geometry2D line-of-sight interpolation.
It is self-contained, scalar (`num_stokes = 1`) and uses two synthetic scenes:

- a standard Geometry2D limb scan: tangents 10–50 km, observer at 600 km,
  350/525/750 nm, Rayleigh + ozone + aerosol, albedo 0.3, a 0.5° atmosphere grid
  over ±20°, 25 source altitudes. It is run for eight dayside solar geometries;
- an orbital-plane track with vertical limb scans, two images per local engine,
  and a solar zenith angle of 35–75° along the track.

For each scheme and horizontal source-column count it reports four quantities:

- the maximum relative multiple-scatter error against that scheme's own
  dense-column reference at the same angular resolution: 81 columns for the
  limb scan and 31 for the orbital-plane scene;
- the wall time of engine and atmosphere construction plus radiance, with and
  without the full Jacobian (orbital-plane timings also include the geometry
  construction). Each timing is a single run, not an average;
- the peak RSS of a single-threaded run in a fresh process;
- the diffuse and line-of-sight source-weight memory from
  `SASKTRAN2_PROFILE_MEMORY`. `parse_memory` tells the two apart by order: it
  assumes their interpolation records alternate, diffuse first, as they do when
  each geometry is constructed in turn.

```sh
python tools/benchmarks/so_interpolation_benchmark.py --output build/so-interp
```

`--directions`, `--columns`, `--cases`, `--reference-columns`,
`--skip-jacobian`, `--skip-memory` and `--skip-orbital` control the run, which
takes about 30 minutes at the default 110 directions on an Apple M4 Pro with
8 threads. Results are written to `results.json` and `results.md`.

## Results at 110 directions

These were measured at commit `fce16966` on an Apple M4 Pro, with 8 threads for
the timings and 1 thread for the memory probes.

Standard Geometry2D limb scan, maximum relative multiple-scatter error:

| Case | Legacy 11 columns | Default 5 columns | Default 7 columns | Default 11 columns |
| --- | ---: | ---: | ---: | ---: |
| SZA 30°, sun ahead | 7.4e-4 | 4.1e-4 | 2.4e-4 | 1.1e-4 |
| SZA 60°, sun ahead | 1.7e-3 | 1.1e-3 | 6.7e-4 | 3.4e-4 |
| SZA 60°, sun behind | 1.5e-3 | 1.1e-3 | 6.7e-4 | 3.3e-4 |
| SZA 60°, sun out of plane | 9.0e-4 | 3.7e-4 | 2.3e-4 | 1.1e-4 |
| SZA 70°, oblique sun | 1.8e-3 | 1.4e-3 | 8.1e-4 | 4.1e-4 |
| SZA 75°, sun ahead | 5.6e-3 | 3.2e-2 | 2.2e-3 | 1.1e-3 |
| SZA 75°, sun behind | 3.7e-3 | 9.3e-3 | 2.2e-3 | 1.0e-3 |
| SZA 80°, sun out of plane | 7.7e-4 | 4.5e-4 | 2.7e-4 | 1.3e-4 |

Seven default columns beat eleven legacy columns in every case. Five do so
except at SZA 75° with the sun in the plane, where the domain edges reach the
terminator. Cost scales with the column count and is the same for both
schemes:

| Columns | Radiance s | Radiance + Jacobian s | Peak RSS MB |
| ---: | ---: | ---: | ---: |
| 5 | 2.3 | 3.0–3.2 | 1015–1022 |
| 7 | 3.2–3.4 | 4.2–4.4 | 1277–1340 |
| 11 | 4.9–5.2 | 6.4–6.9 | 1863–1994 |

Orbital-plane scene (16 images × 9 tangent altitudes):

| Scheme | Columns | Max MS error | Seconds | Peak RSS MB |
| --- | ---: | ---: | ---: | ---: |
| legacy | 11 | 1.2e-2 | 8.5 | 6348 |
| default | 5 | 1.9e-3 | 3.8 | 3029 |
| default | 7 | 9.0e-4 | 5.0 | 3994 |
| default | 11 | 4.0e-4 | 8.3 | 5911 |

Diffuse-ray source weights dominate memory and are unchanged by the default
interpolation. Cubic line-of-sight weights roughly double the line-of-sight
weight memory, from 1.9 MB to 4.5 MB for the 144 orbital-plane lines of sight.

The errors above measure horizontal (column) discretization. The absolute
multiple-scatter error at 110 directions is dominated by angular
discretization: against a 194-direction, 81-column reference it is about
3–4% for both schemes. For SZA 70° with an oblique sun it is 3.7–4.0e-2 for
legacy and 3.4e-2 for default. For SZA 60° with the sun ahead it is
3.0–3.3e-2 and 3.2e-2. The two schemes' dense limits converge with angular
resolution: they differ by 1.3e-3 at 194 directions, against 1.3e-2 at 110 for
the oblique case.
