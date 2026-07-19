# Source-integration accuracy and performance

The source-integration harness compares the same continuous atmosphere sampled on
successively finer altitude grids. It covers standard and volume emission, exact
single scatter, solar and thermal discrete ordinates, successive orders, and the
existing C++ solar and thermal two-stream sources.

## Recommended runs

Use one thread for stable algorithmic comparisons. The full command is deliberately
explicit because the dense references and derivative checks are computationally
expensive.

```bash
pixi run python docs/performance_book/source/internal_validation/source_integration_harness.py \
  --scenarios all \
  --derivatives \
  --finite-difference \
  --timing-samples 15 \
  --num-threads 1 \
  --output build/source-integration/baseline
```

For quick iteration, select one source and omit derivatives:

```bash
pixi run python docs/performance_book/source/internal_validation/source_integration_harness.py \
  --scenarios single_scatter_exact \
  --spacings 2000,1000,500,250 \
  --output build/source-integration/quick
```

For spherical single scatter, use the grid-linear profile to isolate layer-source
integration error from error caused by sampling a curved atmospheric profile on a
coarse grid:

```bash
pixi run python docs/performance_book/source/internal_validation/source_integration_harness.py \
  --scenarios single_scatter_exact \
  --profile-shape grid-linear \
  --spacings 2000,1000,500,250,125 \
  --output build/source-integration/single-scatter-integration-only
```

Use the default `--profile-shape continuous` for end-to-end altitude-grid
convergence. Baseline and candidate reports must use the same profile shape.

To enforce regression gates against a saved report:

```bash
pixi run python docs/performance_book/source/internal_validation/source_integration_harness.py \
  --scenarios all \
  --derivatives \
  --timing-samples 15 \
  --baseline build/source-integration/baseline/report.json \
  --output build/source-integration/candidate \
  --fail-on-gate
```

The output directory contains a complete JSON report and flat CSV files for radiance
accuracy, per-viewing-geometry accuracy, directional derivatives, finite-difference
checks, and scenario summaries. `viewing_accuracy.csv` keeps limb tangent rays and
ground-looking rays separate so an aggregate improvement cannot hide a regression in
either geometry.

## Metrics

For candidate radiance $I$ and the confirmed dense result $I_r$, pointwise error is

$$
e_i = \frac{|I_i-I_{r,i}|}{\max(|I_{r,i}|, 10^{-12}\max|I_r|)}.
$$

Each grid reports maximum absolute error, maximum normalized error, normalized RMS,
and normalized 95th percentile. It also reports the observed refinement order and the
minimum layer count that reaches the requested tolerance. The same statistics are
applied to native weighting functions after contraction with uniform, linear, and
Gaussian altitude perturbations. Optional central differences test those contracted
derivatives independently; bounded profiles use a one-sided difference when needed.

The default coarse spacings are 2000, 1000, 500, 250, and 125 m. Emission and exact
single scatter use 25 m references confirmed against 12.5 m. The more expensive
diffuse solvers use 100 m references confirmed against 50 m. The confirmation error
is reported separately so a result cannot pass by comparison to an unconverged
"truth" calculation.

Timing uses three warmups by default, then reports per-calculation median, median
absolute deviation, minimum, and layer-wavelength-LOS throughput. Keep the same
machine, build type, thread count, scenarios, and timing sample count when comparing
reports.

## Actionable gates

The defaults are intended as merge guards, not claims about final scientific
requirements:

- Dense-reference maximum normalized error no greater than $10^{-4}$.
- No candidate grid's maximum error more than 10% worse than its saved baseline.
- Median calculation time no more than 2% slower than its saved baseline.
- A targeted accuracy change should reduce median grid error by at least 2x and the
  worst targeted case by at least 5x; inspect the CSV rather than averaging unrelated
  source families together.
- Analytic limiting cases should be added as ordinary unit tests with radiance error
  below $10^{-12}$ and directional-derivative error below $10^{-9}$ where the closed
  form is numerically well conditioned.

If the reference-convergence gate fails, halve both dense spacings before judging the
candidate algorithm. If timing noise exceeds the proposed speed change, increase the
sample count rather than relaxing the regression limit.
