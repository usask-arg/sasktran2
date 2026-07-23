---
file_format: mystnb
---

(_source_succesive_orders)=
# The Successive Orders of Scattering Source

SASKTRAN2 has two implementations of successive orders of scattering. The
legacy C++ source is selected with `SuccessiveOrders`, while the newer
fixed-point solver implemented in Rust is selected with
`SuccessiveOrdersRust`.

```{code-cell}
import sasktran2 as sk

config = sk.Config()

config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrdersRust
```

The Rust source uses C++-supplied traced geometry and supports both
{py:class}`sasktran2.Geometry1D` and {py:class}`sasktran2.Geometry2D`.

## Advantages

- Accounts for spherical transport rather than using a plane-parallel
  approximation.
- Supports scalar and polarized calculations.
- The Rust source supports native Jacobian-vector and vector-Jacobian
  products.

## Disadvantages

- Runtime and memory grow with the number of source locations and angular
  quadrature points.
- May require a finer source grid or atmospheric sub-layering in optically
  thick regions.
- The Rust source does not currently form a complete dense Jacobian.

## Source Location Grid

By default, one volume-source altitude is placed at the midpoint of every
atmospheric layer. The source grid can instead be specified independently:

```{code-cell}
import numpy as np

config.successive_orders_altitude_grid_m = np.array(
    [500.0, 2_500.0, 7_500.0, 15_000.0, 30_000.0, 50_000.0]
)
```

The values are in metres and must be finite, strictly increasing, and inside
the atmospheric altitude range. Irregular spacing and a single source
altitude are supported. Setting the option to `None` restores the layer
midpoints. Ground boundary sources are represented separately.

Using fewer source altitudes reduces the number of traced diffuse rays and
the transport state, often substantially reducing setup time, iteration time,
and memory. It also coarsens the vertical representation of the multiple
scatter source, so the grid should be refined until the required radiance
accuracy is reached.

In `Geometry2D`, {py:attr}`sasktran2.Config.num_sza` sets the number of
horizontal source profiles. It is independent of the number of atmospheric
horizontal grid points.

{py:attr}`sasktran2.Config.num_successive_order_points` is a different,
legacy-only option. It subsamples the locations where incoming radiances are
fully traced without changing the source altitude grid. The Rust source
requires its default value of `-1`.

## Initialization Options

The legacy method uses a configured number of iterations, set with
{py:attr}`sasktran2.Config.num_successive_orders_iterations`. It can optionally
be initialized from a discrete-ordinates solution with
{py:attr}`sasktran2.Config.init_successive_orders_with_discrete_ordinates`.

The Rust method stops when its absolute or relative residual tolerance is
reached, or at `successive_orders_max_iterations`. Setting both tolerances to
zero selects a fixed iteration count. Anderson acceleration is controlled by
`successive_orders_anderson_depth`; zero selects damped Picard iteration.

## Linearization

The Rust source provides native JVP and VJP calculations. In converged mode,
the products are calculated from the converged fixed-point system, so Anderson
acceleration can still be used for the primal solution. Fixed-iteration mode
differentiates the configured iterations. A complete dense Jacobian is not
currently available for this source.

The legacy source does not provide the same exact native-product
linearization.

## Relevant Configuration Options

```{eval-rst}
.. autosummary::

  sasktran2.Config.multiple_scatter_source
  sasktran2.Config.num_sza
  sasktran2.Config.successive_orders_altitude_grid_m
  sasktran2.Config.num_successive_orders_incoming
  sasktran2.Config.num_successive_orders_outgoing
  sasktran2.Config.successive_orders_max_iterations
  sasktran2.Config.successive_orders_relative_tolerance
  sasktran2.Config.successive_orders_absolute_tolerance
  sasktran2.Config.successive_orders_anderson_depth
  sasktran2.Config.successive_orders_damping
  sasktran2.Config.num_successive_orders_iterations
  sasktran2.Config.num_successive_order_points
  sasktran2.Config.init_successive_orders_with_discrete_ordinates
  sasktran2.Config.num_stokes

```
