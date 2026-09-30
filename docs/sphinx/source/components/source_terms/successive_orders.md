---
file_format: mystnb
---

(_source_succesive_orders)=
# The Successive Orders of Scattering Source
The successive orders of scattering source is an implementation of the multiple scattering source using the successive orders
of scattering method. It can be enabled with

```{code-cell}
import sasktran2 as sk

config = sk.Config()

config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
```

## Advantages

 - Fully accounts for sphericity of the atmosphere
 - Supports full weighting-function calculations and memory-efficient native JVP/VJP products

## Disadvantages

 - Can use a large amount of RAM for large, polarized calculations
 - May require sub-layering in optically thick areas of the atmosphere


## Iteration Options

The successive orders method starts from the exact single-scatter illumination and repeatedly
applies its transport and scattering operators. The maximum iteration count is
configured with
{py:attr}`sasktran2.Config.num_successive_orders_iterations`. By default it
stops early when the configured absolute or relative residual tolerance is
met, and Anderson acceleration is enabled.

Set both tolerances to zero to perform exactly the configured number of
iterations.

The explicit {py:attr}`sasktran2.Config.successive_orders_altitude_grid_m`
option decouples the source grid from the atmosphere altitude grid. When it is
not set, the source uses the midpoints of the atmosphere layers.

## Memory Use

Scalar calculations keep transport values and solver workspace per wavelength
worker. The compact scalar first-order caches also use one reusable slot per
worker, so their storage scales with active wavelength concurrency. Geometry
and transport topology are shared across wavelengths. Derivative-only buffers
are allocated when needed. Completed geometry arrays release unused capacity,
and compact interpolation weight records omit alignment padding while
preserving all 64 bits of each double.

Each wavelength retains its forward diffuse solution and compact scalar direct
forcing. Native JVP/VJP products reuse these values, and tolerance-controlled
solves can use the diffuse solution as a warm start
after atmosphere updates. Fixed-iteration solves continue to start from zero
after updates.

Reusing a worker for another wavelength recomputes the wavelength-dependent
transport and first-order derivative quantities. Forcing is reused only while
the atmosphere revision and geometry remain current. This trades repeated
assembly for lower
resident memory while preserving the angular and spatial grids and each
wavelength's physical values. Surface-only updates also recompute volume
transport and first-order forcing instead of retaining duplicate volume
buffers.

The optional
{py:attr}`sasktran2.Config.successive_orders_transport_cache_wavelengths`
setting pins transport values for the first N wavelength indices on each scalar
compact worker. Its default, zero, keeps the current active-only policy. Cached
and active vectors exchange ownership, so retaining all N wavelengths adds
N-1 transport vectors per worker. Partial caches can additionally retain an
active wavelength outside the pinned set. Counts are clamped to the available
wavelengths; increasing them trades resident memory for less transport assembly.
Cache entries are invalidated on atmosphere and geometry updates. Calculation
order, warm starts and gradient accumulation order are unchanged. Polarized and
noncompact paths retain their existing storage policy.

Finalized source interpolation uses byte or 16-bit CSR slots when the row fits,
with a 32-bit fallback for larger rows. Verified structured 2D ray cells reconstruct
their original corner indices from a compact descriptor, retaining every original
double coefficient and its order. Other ray stencils keep explicit indices.

Orbital native products preserve absent phase derivative mappings when copying
requested parameters into local atmospheres. Changes in phase-mapping presence
rebuild the local derivative storage. For scalar calculations with no native
phase derivatives and a spatial Lambertian surface, VJP assembly omits unused
scattering parameter gradients while retaining the full configured angular
width for the forcing cotangent and the same adjoint iterations.

Geometry and ray transport maps share immutable CSR generations rather than
copying their index arrays. A geometry refresh creates a new generation;
existing handles keep the old generation alive until their users release it.
Point-specific scalar incoming transforms also share the outgoing transform
when their outgoing sphere is the same object. Scattering reads contiguous
point inputs directly, and angular routines size temporary moment arrays to
the batch they actually use.

## Structured 2D Geometry

The successive-orders source supports horizontally varying atmospheres defined
with {py:class}`sasktran2.Geometry2D`. Its source grid is independent of the
atmosphere grid:

- {py:attr}`sasktran2.Config.num_sza` selects the number of evenly spaced
  horizontal source columns spanning the atmosphere's horizontal-angle grid
  when no explicit horizontal source grid is supplied.
- {py:attr}`sasktran2.Config.successive_orders_horizontal_angle_grid_radians`
  selects the exact local Geometry2D angles of the horizontal source columns.
  In an orbital-plane engine, zero is the center of each group's fitted local
  plane. The supplied points must lie within every local group's horizontal
  grid.
- {py:attr}`sasktran2.Config.successive_orders_altitude_grid_m` selects the
  vertical source locations. When it is not set, atmosphere-layer midpoints are
  used.

Optional
{py:attr}`sasktran2.Config.successive_orders_incoming_directions_by_altitude`
and
{py:attr}`sasktran2.Config.successive_orders_outgoing_directions_by_altitude`
arrays select a direction count for each resolved source altitude. Empty arrays
or `None` retain the uniform counts. Every horizontal column uses the same
altitude profile, and ground points keep the uniform settings and existing
hemisphere rules.
Horizon-fitted incoming rules require at least six directions. Explicit outgoing
counts must match supported quadrature rules. Profile lengths are checked when
the source altitude grid is resolved.

An explicit uniform profile follows the existing uniform numerical path. A
nonuniform profile changes angular discretization and requires a convergence
study against a uniform reference. Keep the profile fixed throughout retrieval
so native JVP/VJP products differentiate the same discrete model after atmosphere
updates. Shared outgoing transforms are retained for identical quadrature rules.

The incoming direct beam is obtained from a solar-characteristic table
parameterized by altitude, solar zenith angle, and off-plane azimuth. Setting
{py:attr}`sasktran2.Config.solar_refraction` bends these solar characteristics
using the refractive-index profile stored on the geometry. When exact or table
single scattering is also enabled, the sources share this table.

Geometry2D successive orders provides native JVP and VJP calculations, which
avoid allocating the complete structured radiance Jacobian. Refraction of the
diffuse multiple-scatter rays is not yet supported, so
{py:attr}`sasktran2.Config.multiple_scatter_refraction` must remain disabled.
Line-of-sight refraction and flux observers are also not supported with
Geometry2D.

## Relevant Configuration Options

```{eval-rst}
.. autosummary::

  sasktran2.Config.multiple_scatter_source
  sasktran2.Config.num_sza
  sasktran2.Config.num_successive_orders_iterations
  sasktran2.Config.num_successive_orders_incoming
  sasktran2.Config.num_successive_orders_outgoing
  sasktran2.Config.successive_orders_relative_tolerance
  sasktran2.Config.successive_orders_absolute_tolerance
  sasktran2.Config.successive_orders_anderson_depth
  sasktran2.Config.successive_orders_damping
  sasktran2.Config.successive_orders_altitude_grid_m
  sasktran2.Config.successive_orders_transport_cache_wavelengths
  sasktran2.Config.successive_orders_incoming_directions_by_altitude
  sasktran2.Config.successive_orders_outgoing_directions_by_altitude
  sasktran2.Config.successive_orders_horizontal_angle_grid_radians
  sasktran2.Config.num_stokes

```
