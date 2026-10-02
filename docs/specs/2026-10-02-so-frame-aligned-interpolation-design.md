# Frame-aligned source grids and cubic LOS interpolation for successive orders

Status: implemented on branch `claude/successive-orders-optimization` (2026-10-02).

## Problem

Geometry2D successive orders stores the multiple-scatter source on a grid of
altitudes and horizontal columns and interpolates it along every diffuse and
observer ray. Memory and runtime are roughly linear in the number of columns
(PR #304: OMPS orbital plane, 5 columns 9.75 GB, 11 columns 18.85 GB), so the
goal is to reach today's 11-column accuracy with far fewer columns for SZA
below about 85 degrees. Terminator conditions are out of scope.

## Findings that drive the design

Measured on a limb scan (tangents 10-50 km, 350/525/750 nm, Rayleigh, ozone,
aerosol, albedo 0.3, 25 source altitudes, 50 directions) against 161-column
references:

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
   today's 11 linear columns in every dayside case tested, with no measurable
   memory cost in the prototype (1350 MB vs 1884 MB peak RSS at 110
   directions). The shipped cubic LOS weights do add memory; see section 4.
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

Every aligned Lebedev rule, interior and ground, is rotated by
`local_solar_frame(p) * pole_avoiding_rotation()`. The fixed pole-avoiding
rotation keeps nodes off the canonical poles and mirror ties, so no node lies on
the local vertical. This is not "the legacy sphere rotated": legacy grids applied
`pole_avoiding_rotation()` to the reduced-horizon outgoing rule only, and the
legacy pre-rotation for plain Lebedev grids was the identity. A
`RotatedLebedevSphere` wraps a canonical `LebedevSphere` and a rotation; its
positions are `R * q` and `interpolate` evaluates the canonical sphere at
`R^T d`.

All interior points in one column share `up` and the solar reference, hence one
rotation and one shared sphere object. Ground points use the frame at their own
column location for both the incoming and outgoing spheres. Aligned ground
hemispheres reject nodes within `1e-12 * |location|` of the horizon; the legacy
ground grids keep a zero tolerance.

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

`SourcePoint` gains an `angular_class` index: interior points with equal index
have identical grids up to a rigid rotation of the same canonical nodes. It
equals the altitude index in reduced-horizon mode and is zero in Lebedev mode.
Ground points have `angular_class == -1` and are not shared.

Bitwise equality with legacy in an identity frame is not claimed: the aligned
grids always include the pole-avoiding rotation, so aligned 1D plain-Lebedev
single-column results differ from legacy at about 1e-3. The legacy switch is the
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

Default results change at the angular-discretization level, up to about 1e-2 at
26 directions and much less at 110, for:

- Geometry2D with several columns;
- spherical Geometry1D, including single-column plain-Lebedev cases (about
  1e-3 against legacy) and any solar azimuth. Aligned 1D results are invariant
  to the solar-azimuth convention to about 1e-8, against up to 4e-2 with legacy
  grids;
- multi-SZA Geometry1D.

Measured accuracy: a dayside convergence test gives a maximum relative
multiple-scatter error of 5.8e-5 for aligned + cubic at 7 columns, against
5.1e-4 for legacy at 11 columns. Typical dayside limb cases at 50 directions
give 2e-4 to 2.5e-3 (aligned + cubic, 7 columns) against 2.5e-3 to 7.5e-3
(legacy, 11 columns).

Costs: frame alignment has no memory cost. Cubic LOS doubles the LOS source
weights, and bytes grow about 2.15x because negative weights are
escape-encoded: roughly 12 KB to 25 KB per LOS, about +135 MB for 10k LOS.
Diffuse-ray weights, which dominate orbital-plane memory, are unchanged.

### Related fixes

- The C++ spherical ray tracer clamps `1 - cos^2` at zero for exactly radial
  rays (`7638de57`).
- A typo in an unreachable branch of `add_od_quadrature` (the `t1 < t0`
  near-radial OD quadrature) is fixed (`494d2420`).

## Testing

- C++: aligned sphere positions and weights; equal canonical coordinates of a
  rotated LOS direction at different columns; cubic weights reproduce cubic
  polynomials exactly, sum to one, fall back correctly at edges and with fewer
  than four columns; scalar basis sharing reproduces the per-point operator;
  vector grouped synthesis equals per-point apply.
- Python: existing successive-orders, 2D, orbital-plane and linearization tests
  pass (expected-value updates only where results legitimately change);
  JVP/VJP adjoint and finite-difference checks for scalar and vector with both
  quadrature modes; legacy switch reproduces `main` bitwise; horizontal
  convergence regression showing seven aligned/cubic columns within the error
  of eleven legacy columns for a dayside case.

## Benchmarks (deliverable)

Timing, peak memory and accuracy versus legacy for 5, 7 and 11 columns at 110
directions: a standard Geometry2D limb scan across several dayside solar
geometries, plus an orbital-plane run. A reproducible script will live under
`tools/benchmarks/`.

## Risks

- Vector per-column synthesis can be slower than today's single shared
  synthesis; measured in the benchmarks.
- Regression baselines that pin multiple-scatter values will change and must be
  regenerated deliberately, never loosened.

## Known limitations / follow-ups

- A pre-existing vector Geometry2D jump of about 1e-3 under 1e-9 rad
  solar-azimuth perturbations, also present on `main`, is being investigated
  separately.
- Vector I differs from scalar I by 0.15-1.4% in 2D successive orders even with
  decoupled polarization; pre-existing.
- The MODIS/SnowKokhanovsky BRDF evaluates an unclamped `sqrt(1 - mu*mu)`; being
  investigated separately.
- Per-column vector synthesis sharing could save memory.
- `geometry.cpp` could be split.
