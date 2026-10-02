# Frame-aligned source grids and cubic LOS interpolation for successive orders

Status: approved direction (2026-10-02). Branch: `claude/successive-orders-optimization`.

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
   memory cost (1350 MB vs 1884 MB peak RSS at 110 directions).
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
R(p) = [x, up x x, up],  up = p / |p|,  x = solar_horizontal_reference(up)
```

`solar_horizontal_reference` already defines the frame used by
`rotate_unit_vector` and by `ReducedHorizonSphere`, so aligned grids agree with
the existing LOS direction mapping: a direction mapped to a column keeps its
local zenith and solar azimuth, and therefore has the same canonical
coordinates at every column.

A `FrameAlignedSphere` wraps a canonical `LebedevSphere` and a rotation; its
positions are `R * q` and `interpolate` evaluates the canonical sphere at
`R^T d`. It replaces `PoleAvoidingLebedevSphere`, whose fixed rotation becomes
the canonical pre-rotation `R(p) * R_pole` (pole avoidance is retained for the
reduced-horizon and vector paths exactly as today).

All interior points in one column share `up` and the solar reference, hence one
rotation and one shared sphere object. Ground points use their own column's
rotation for both the wrapped incoming and outgoing spheres.

- Reduced-horizon mode (default): outgoing spheres become per-column aligned
  spheres; incoming `ReducedHorizonSphere`s are unchanged (already local).
- Lebedev mode: incoming and outgoing spheres are both per-column aligned
  spheres with the same rotation.
- Plane-parallel and pseudo-spherical geometries keep the current global grids;
  all their points share one local frame, so alignment cannot change
  interpolation consistency.

`SourcePoint` gains an `angular_class` index: points with equal index have
identical grids up to a rigid rotation of the same canonical nodes. It equals
the altitude index in reduced-horizon mode and is zero in Lebedev mode; ground
points get their own classes.

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
  containing the angle.
- Constant extension outside the column range; linear when fewer than four
  columns exist.
- Applies to observer LOS layer sources and LOS ground-hit weights. Diffuse
  incoming rays and ground forcing remain bilinear.
- The mode is an explicit argument of the location interpolator's weight
  functions and of `compile_ray_interpolation`, chosen by the caller for LOS
  versus incoming compilation. No global state; safe under OpenMP.
- Geometry1D (cos SZA columns) is unchanged.

Weights are geometry-only, so primal, JVP, VJP and full-Jacobian paths need no
new derivative terms. Negative weights are already supported by the source
weight encoding.

### 4. Configuration and compatibility

`Config.successive_orders_legacy_interpolation` (bool, default `False`) restores
global grids, the legacy basis layout and linear LOS interpolation. It is
plumbed through C++ `Config`, the C API, `sasktran2-sys`, `sasktran2-rs`,
`sasktran2-py-ext`, `sasktran2.Config`, and the orbital-plane structural config
signature.

Default results change at the angular-discretization level:

- Geometry2D with several columns, about 1e-3 at 50 directions;
- spherical Geometry1D with a solar azimuth other than zero;
- multi-SZA Geometry1D.

Spherical Geometry1D with one column and zero solar azimuth is bitwise unchanged
because its frame is the identity.

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
