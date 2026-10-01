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
