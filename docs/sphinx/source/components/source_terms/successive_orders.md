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

Scalar calculations keep wavelength-dependent transport values and compact
first-order caches in reusable wavelength-worker slots. Geometry and transport
topology are shared across wavelengths. Derivative-only buffers are allocated
when needed. These changes preserve the configured angular and spatial grids.

Each wavelength retains its forward diffuse solution and compact scalar direct
forcing. Native JVP/VJP products reuse these values, and tolerance-controlled
solves can use the diffuse solution as a warm start
after atmosphere updates. Fixed-iteration solves continue to start from zero
after updates.

Reusing a worker for another wavelength recomputes the wavelength-dependent
transport and first-order derivative quantities. Forcing is reused only while
the atmosphere revision and geometry remain current. This trades repeated
assembly for lower resident memory while preserving each wavelength's physical
values. Surface-only updates also recompute volume
transport and first-order forcing instead of retaining duplicate volume
buffers.

Finalized source interpolation uses byte or 16-bit CSR slots when the row fits,
with a 32-bit fallback for larger rows. Compact weight records reconstruct every
original double bit pattern; weights that cannot use the compact representation
retain all eight bytes. Verified structured 2D ray cells reconstruct their original
corner indices from a compact descriptor. Other ray stencils keep explicit
indices. Neither representation combines or reorders contributions.

Orbital native products preserve absent phase derivative mappings when copying
requested parameters into local atmospheres. Changes in phase-mapping presence
rebuild the local derivative storage. For scalar calculations with no native
phase derivatives and a spatial Lambertian surface, VJP assembly omits unused
scattering parameter gradients while retaining the full configured angular
width for the forcing cotangent and the same adjoint iterations.

Uniform-phase values are stored only for wavelengths whose phase functions are
spatially uniform. Native phase derivative scratch is allocated only when phase
mappings are present; ozone-only native products omit those unused phase sums.
Changes in atmosphere volume and mapping presence refresh these buffers before
reuse.

For eligible scalar 2D calculations with one thread, finalized endpoint stencils
store two original interpolation coordinates and one verified cell base.
Factoring is adopted only when expanding those coordinates reproduces every
original weight bit for the entire provider. Ineligible providers retain their
four original weights or explicit stencils. If the floating-point rounding mode
changes, the provider materializes the original weights before reuse.

Shared transport CSR columns use 16 bits when the complete source grid fits.
Solar interpolation uses
one-byte row counts when each row has at most 255 entries, and row-relative
16-bit column indices when every span fits and the representation saves memory.
Wider rows and indices retain their original forms. Each product selects its
storage format before its accumulation loop, preserving the original entry
order. Finalization also releases unused construction capacity.

When many solar interpolation rows have the same ordered column indices, the
immutable topology stores each pattern once and gives each row a compact pattern
ID. The double-valued weights remain independent and in their original order.
Pattern IDs use 16 bits when possible, with wider IDs or the existing row
representations as fallbacks. Pattern storage is selected only when it reduces
retained memory; it does not combine or reorder floating-point contributions.

Scalar first-order products, diffuse solves and line-of-sight derivative
products lease temporary workspace from calling-thread arenas. Resident local
engines can reuse this scratch because no product retains a view after its call
returns. First-order arrays grow to the largest requested size and products use
only their active ranges. Line-of-sight arrays follow the active transport shape;
solver history is reset for each solve. Source threads receive separate solar
cotangent ranges, and nested calls
use private scratch while the arena is busy. Physical transmission and medium
caches, derivative lifetimes, and each wavelength's warm-start state remain
owned by their engines.

Geometry and ray transport maps share immutable CSR generations rather than
copying their index arrays. A geometry refresh creates a new generation;
existing handles keep the old generation alive until their users release it.
Point-specific scalar incoming transforms also share the outgoing transform
when their outgoing sphere is the same object. Scattering reads contiguous
point inputs directly, and angular routines size temporary moment arrays to
the batch they actually use. Finalization releases unused construction capacity
and the temporary endpoint-coordinate capture used to verify factoring.

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
  sasktran2.Config.successive_orders_horizontal_angle_grid_radians
  sasktran2.Config.num_stokes

```
