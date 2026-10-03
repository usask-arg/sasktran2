# Frame-aligned source grids and cubic LOS interpolation for successive orders

Status: implemented on branch `claude/successive-orders-optimization` (2026-10-02).

## Problem

Geometry2D successive orders stores the multiple-scatter source on a grid of
altitudes and horizontal columns and interpolates it along every diffuse and
observer ray. Memory and runtime are roughly linear in the number of columns
(OMPS orbital plane, peak physical footprint at the validated checkpoint
`60de1b4e` in `tools/benchmarks/README.md`: 5 columns 9.747 GB, 11 columns
18.245 GB; PR #304 quoted 18.854 GB for 11 columns), so the goal is to reach today's 11-column accuracy with far fewer columns for SZA
below about 85 degrees. Terminator conditions are out of scope.

## Findings that drive the design

Prototype measurements, taken before the shipped implementation, on a limb
scan (tangents 10-50 km, 350/525/750 nm, Rayleigh, ozone, aerosol, albedo 0.3,
25 source altitudes, 50 directions) against 161-column references:

1. Error is dominated by interpolation between columns, not by the solution
   at the columns.
2. The outgoing Lebedev grid is fixed in global coordinates. Each column sees
   it rotated differently relative to its local vertical and solar azimuth, so
   the same local direction is interpolated from different nodes and weights at
   neighbouring columns. This produces a column-count-independent error floor
   (2-7e-3 of the multiple-scatter radiance) that more columns or higher-order
   horizontal interpolation cannot remove.
3. Rotating each column's outgoing grid into its local frame removes the floor
   without changing angular accuracy (at 41 columns the error against a
   194-direction truth is 3.2e-2 for both frames at 110 directions).
4. With aligned grids, 4-point cubic horizontal interpolation of the LOS source
   is 3-20x more accurate than linear at the same column count. Seven columns
   with aligned grids and cubic LOS interpolation are at least as accurate as
   today's 11 linear columns in every dayside case tested. In the prototype
   seven aligned/cubic columns peaked at 1350 MB RSS against 1884 MB for
   eleven legacy columns at 110 directions; the saving comes from the smaller
   column count, not from cubic weights being free. The shipped cubic LOS
   weights do add memory; see section 4.
5. Cubic interpolation of diffuse incoming rays is a further accuracy gain but
   doubles transport weights; it is not part of this change.

## Scope

In scope:

- Frame-aligned outgoing (and, in Lebedev mode, incoming) angular grids for
  spherical Geometry1D and Geometry2D successive orders, scalar and vector.
- Angular-basis sharing that is correct for per-column grids and exploits the
  rotation invariance of the scalar operator.
- Cubic horizontal interpolation of the observer-LOS source for Geometry2D,
  including orbital-plane group engines.
- One `Config` switch restoring the legacy behaviour.

Out of scope: cubic diffuse-ray interpolation, automatic column placement,
terminator-specific interpolation, angular-interpolation improvements.

## Design

### 1. Frame-aligned angular grids (`geometry.cpp`)

The local frame of a source point at position `p` is the rotation

```
local_solar_frame(p) = [x, up x x, up],  up = p / |p|,  x = solar_horizontal_reference(up)
```

`solar_horizontal_reference` already defines the frame used by
`rotate_unit_vector` and by `ReducedHorizonSphere`, so aligned grids agree with
the existing LOS direction mapping: a direction mapped to a column keeps its
local zenith and solar azimuth, and therefore has the same canonical
coordinates at every column.

Every aligned interior Lebedev rule is rotated by
`local_solar_frame(p) * pole_avoiding_rotation()`. The fixed pole-avoiding
rotation keeps nodes off the canonical poles and mirror ties, so no node lies on
the local vertical. This is not "the legacy sphere rotated": legacy grids applied
`pole_avoiding_rotation()` to the reduced-horizon outgoing rule only, and the
legacy pre-rotation for plain Lebedev grids was the identity. Aligned ground
Lebedev rules are rotated by `local_solar_frame(p) * ground_pre_rotation()`,
where `ground_pre_rotation() = pole_avoiding_rotation() * Rx(0.08 rad)`; see
the ground paragraph below. A
`RotatedLebedevSphere` wraps a canonical `LebedevSphere` and a rotation; its
positions are `R * q` and `interpolate` evaluates the canonical sphere at
`R^T d`.

All interior points in one column share `up` and the solar reference, hence one
rotation and one shared sphere object. Ground points use the frame at their own
column location for both the incoming and outgoing spheres.

Ground hemispheres, legacy and aligned, follow main's rule (#306): nodes within
`1e-12 * |location|` of the horizon are excluded but keep half their weight in
the normalization. The ground point's local vertical is the canonical z axis,
so `pole_avoiding_rotation()` alone would leave the canonical y-axis Lebedev
pair exactly on the horizon, and that rule would then change plain-Lebedev
ground reflection by about 10% at 26 nodes compared with plain exclusion. The
extra 0.08 rad tilt about x keeps every node of every Lebedev rule up to 302
points at least 1e-2 from the horizon (0.01 rad would leave the (1, 1, 0)
nodes about 3.5e-7 from it), so aligned ground hemispheres contain exactly half
the nodes and none is ever within the tolerance. The Lambertian factor
`4 sum(mu w)` of the tilted outgoing hemisphere, 1 for an exact rule, is
0.950, 0.9992, 1.0010 and 1.0008 at 26, 50, 110 and 194 nodes. Legacy grids
keep the identity orientation and therefore main's half-weighted horizon
nodes.

- Reduced-horizon mode (default): outgoing spheres become per-column aligned
  spheres; incoming `ReducedHorizonSphere`s are unchanged (already local).
- Lebedev mode: incoming and outgoing spheres are both per-column aligned
  spheres with the same rotation.
- Ground points: aligned plain-Lebedev outgoing spheres, and aligned
  plain-Lebedev incoming spheres in Lebedev mode.
- Plane-parallel and pseudo-spherical geometries keep the current global grids;
  all their points share one local frame, so alignment cannot change
  interpolation consistency.

One frame per column is valid only if interior points are altitude-fastest
within columns that share one direction, and equal altitude indices share one
radius. A premise check at construction throws `std::logic_error` otherwise.
`SourceGeometry1D` always builds its source points in that order, so the check
guards future changes and cannot be reached through public interfaces.

`SourcePoint` gains an `angular_class` index: interior points with equal index
have identical grids up to a rigid rotation of the same canonical nodes. It
equals the altitude index in reduced-horizon mode and is zero in Lebedev mode.
Ground points have `angular_class == -1` and are not shared.

Bitwise equality with legacy in an identity frame is not claimed: the aligned
grids always include the pole-avoiding rotation (and the ground tilt), so
aligned 1D plain-Lebedev single-column results differ from legacy at about
1e-3. The legacy switch is the
only bitwise-reproducing path.

### 2. Angular-basis sharing (`scattering_assembler.cpp`, `scattering.cpp`)

Scalar: the operator matrix depends only on the angles between incoming and
outgoing nodes, which a common rigid rotation preserves. One
`ScalarAngularBasis` per angular class is therefore exact for every column of
that class. All classes share one synthesis transform. This reduces analysis
storage from one matrix per point to one per altitude.

Vector: Q/U are defined relative to global spin-harmonic frames and transport
copies Stokes vectors without rotation, so bases cannot be shared across
columns. Each column gets its own synthesis and each point its own analysis.
The forward apply batches synthesis per group of points that share an outgoing
basis instead of assuming one global synthesis. This removes the hard-coded
`point_bases_share_synthesis = true`, the cause of the prototype's JVP/VJP
mismatch. The non-reduced vector path also moves to per-column bases.

### 3. Cubic LOS interpolation (`AltitudeAngleSourceLocationInterpolator`)

- 4-point Lagrange weights in horizontal angle on the (possibly non-uniform)
  column grid, stencil `[i-1, i+2]` clamped to `[0, n-4]` for the interval
  containing the angle (`horizontal_interpolation.h`).
- Constant extension outside the column range; linear when fewer than four
  columns exist.
- Applies to observer LOS layer sources and to LOS ground-hit weights, which use
  the same four-column stencil. Diffuse incoming rays and ground forcing remain
  bilinear.
- The LOS stencil is provided by a second, LOS-only location interpolator
  (`m_los_location_interpolator`, null when unused) rather than a mode argument
  on the weight functions. `compile_los_interpolation` uses it when present and
  the diffuse-ray interpolator otherwise. No global state; safe under OpenMP.
- Sunlit-stencil guard: the stencil is cubic only if the solar zenith angle is
  below 90 degrees at every one of its four columns (geometry-only,
  `up . sun_unit > 0`); otherwise that location uses linear weights. This
  avoids the negative radiances seen near the terminator at coarse column
  spacing, where Lagrange weights mix a dim night-side column with bright
  sunlit ones. The guard reduces rather than guarantees freedom from such
  artefacts: a dim but sunlit column (cos SZA near zero) next to bright ones
  can still give a negative interpolated source.
- Weight bound: the stencil is also cubic only if its absolute weights sum to
  at most `max_cubic_weight_abs_sum` = 2 (`horizontal_interpolation.h`);
  otherwise that location uses linear weights. Lagrange weights grow without
  bound on strongly non-uniform explicit grids (0.25, -1, 1.5, 0.25 on
  `[0, 1, 2, 4]` at 3; about +-38 on `[0, 1, 1.01, 2]` at 0.5), which the
  sunlit guard does not catch. On a uniform grid the sum is at most 1.25 in
  interior intervals and about 1.63 in the end intervals, so uniform grids
  always stay cubic. On the tangent-clustered grid `[-20, -8, -3, 0, 3, 8, 20]`
  degrees the sum stays below 1.40 in the inner four intervals and reaches
  about 5.1 in parts of the outer two, which fall back to linear there.
- Geometry1D (cos SZA columns) is unchanged.

Weights are geometry-only, so primal, JVP, VJP and full-Jacobian paths need no
new derivative terms. Negative weights are already supported by the source
weight encoding (escape-encoded, which is where the extra memory comes from).

### 4. Configuration and compatibility

`Config.successive_orders_legacy_interpolation` (bool, default `False`) restores
global grids, the legacy basis layout and linear LOS interpolation. It is
plumbed through C++ `Config`, the C API, `sasktran2-sys`, `sasktran2-rs`,
`sasktran2-py-ext`, `sasktran2.Config`, and the orbital-plane structural config
signature.

Default results change at the angular-discretization level. In the
benchmark, the default and legacy dense-column limits differ by up to about
1e-2 in multiple-scatter radiance at 110 directions (1.3e-2 in the oblique
case), shrinking to about 1e-3 at 194 directions. This applies to:

- Geometry2D with several columns;
- spherical Geometry1D, including single-column plain-Lebedev cases (about
  1e-3 against legacy) and any solar azimuth. Aligned 1D results are invariant
  to the solar-azimuth convention to about 1e-7 (a 7e-8 jump that
  single-scatter-only runs share), against 8e-4 to 7e-3 (scalar) and 3e-2 to
  4e-2 (vector) with legacy grids in the invariance test, so Geometry1D runs
  with a non-zero solar azimuth can change by more than the
  angular-discretization level;
- multi-SZA Geometry1D.

Measured accuracy: a dayside convergence test gives a maximum relative
multiple-scatter error of 5.8e-5 for aligned + cubic at 7 columns, against
5.1e-4 for legacy at 11 columns. Typical dayside limb cases at 50 directions
give 2e-4 to 2.5e-3 (aligned + cubic, 7 columns) against 2.5e-3 to 7.5e-3
(legacy, 11 columns). In the scalar benchmark at 110 directions
(`tools/benchmarks/README.md`; dayside cases with SZA 30-80 degrees), seven
default columns beat eleven legacy columns by factors of 1.7-3.9 in maximum
multiple-scatter error, at about 30-35% less time and memory. The benefit
depends on every LOS stencil column being sunlit, so the domain width and the
terminator position matter.

Costs: frame alignment has no memory cost for scalar calculations (bases are
shared per angular class) or with reduced-horizon quadrature (bases were
already per point). Plain-Lebedev vector calculations previously shared one
angular basis and now build one per column. Cubic LOS doubles the number of LOS
source weights, and bytes grow by more than that because negative weights are
escape-encoded. The per-LOS cost depends on the scene: in the orbital-plane
benchmark the LOS source weights grow from 1.9 MB to 4.5 MB for 144 LOS, and
in a reviewer's 2000-LOS Geometry2D case from 23.5 MB to 50.5 MB, about 12 KB
to 25 KB per LOS (about +135 MB for 10k such LOS). Diffuse-ray weights, which
dominate orbital-plane memory, are unchanged.

### Related fixes

- The C++ spherical ray tracer clamps `1 - cos^2` at zero for exactly radial
  rays (`7638de57`).
- A typo in an unreachable branch of `add_od_quadrature` (the `t1 < t0`
  near-radial OD quadrature) is fixed (`494d2420`).
- Merged from main: ground hemispheres exclude roundoff horizon nodes with half
  weight (#306), and the 2D tracer no longer truncates upward rays from ground
  endpoints (#307). The latter caused the vector Geometry2D solar-azimuth jump
  and the vector-vs-scalar I bias listed as known limitations before the
  merge.

## Testing

- C++: aligned sphere positions and weights; equal canonical coordinates of a
  rotated LOS direction at different columns; cubic weights reproduce cubic
  polynomials exactly, sum to one, fall back correctly at edges and with fewer
  than four columns; scalar basis sharing reproduces the per-point operator;
  vector grouped synthesis equals per-point apply; aligned ground Lebedev
  rules up to 302 points have no node within 1e-6 of the horizon;
  `refresh_los` recompiles cubic LOS stencils and ground-hit stencils.
- Python: existing successive-orders, 2D, orbital-plane and linearization tests
  pass (expected-value updates only where results legitimately change);
  JVP/VJP adjoint and finite-difference checks for scalar and vector with both
  quadrature modes; horizontal convergence regression showing seven
  aligned/cubic columns within the error of eleven legacy columns for a dayside
  case; #306 regressions in both interpolation modes.
- Legacy reproduces main: the legacy switch is bitwise identical to main except
  where main produced NaN for exactly radial rays (ray-tracer clamp in a shared
  header). This was verified manually on 224 arrays against main `942d5494`
  (1D, 2D and orbital plane; scalar and vector; both quadratures; Lambertian
  and MODIS surfaces; radiance and all weighting functions).
  `test_legacy_interpolation_reproduces_main` pins a small set of main's
  radiances (1D scalar/vector with reduced horizon on and off, 2D with seven
  columns) at `rtol = 1e-10`.

## Benchmarks (deliverable)

Timing, peak memory and accuracy versus legacy for 5, 7 and 11 columns at 110
directions: a standard Geometry2D limb scan across several dayside solar
geometries, plus an orbital-plane run. A reproducible script will live under
`tools/benchmarks/`.

## Risks

- Vector per-column synthesis can be slower than today's single shared
  synthesis. Vector timing was not benchmarked; the benchmark is scalar.
- Regression baselines that pin multiple-scatter values will change and must be
  regenerated deliberately, never loosened.

## Known limitations / follow-ups

Left open by #307 and pre-existing on main:

- Successive-orders forcing Stokes-frame mismatch: `PhaseHandler` uses
  local-vertical reference frames while `VectorAngularBasis` uses global-z
  frames.
- Geometry2D single-scatter error when an LOS ground hit lands on a grid
  corner.
- The MODIS/SnowKokhanovsky BRDF evaluates an unclamped `sqrt(1 - mu*mu)`.

Other:

- Per-column vector synthesis sharing could save memory.
- `geometry.cpp` could be split.
