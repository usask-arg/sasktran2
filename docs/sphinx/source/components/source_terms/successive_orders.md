---
file_format: mystnb
---

(_source_succesive_orders)=
# The Successive Orders of Scattering Source
SASKTRAN2 provides two implementations of the successive orders of scattering
method. The newer, self-contained C++ source can be enabled with

```{code-cell}
import sasktran2 as sk

config = sk.Config()

config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrdersCpp
```

The existing implementation remains available as
`sk.MultipleScatterSource.SuccessiveOrders`.

## Advantages

 - Multiple scattering that fully accounts for spherical source geometry
 - Full-Jacobian calculations and native JVP/VJP products
 - An explicit source-altitude grid can reduce cost independently of the atmosphere grid

## Disadvantages

 - Currently supports one-dimensional atmospheres only
 - Can use a large amount of RAM for large polarized calculations
 - May require sub-layering in optically thick areas of the atmosphere


## Iteration and convergence

The C++ source starts from the exact single-scatter illumination and repeatedly
applies its transport and scattering operators. The maximum iteration count is
configured with
{py:attr}`sasktran2.Config.num_successive_orders_iterations`. By default it
stops early when the configured absolute or relative residual tolerance is
met, and Anderson acceleration is enabled.

Set both tolerances to zero to perform exactly the configured number of
iterations. In this mode each new atmospheric state starts from zero so the
result is independent of engine history. Native JVPs and VJPs differentiate
the converged fixed-point system; when tolerances are disabled they need not
equal the derivative of the finite iteration sequence.

The explicit {py:attr}`sasktran2.Config.successive_orders_altitude_grid_m`
option decouples the source grid from the atmosphere altitude grid. When it is
not set, the source uses the midpoints of the atmosphere layers.

## Legacy initialization option

The discrete-ordinates initialization controlled by
{py:attr}`sasktran2.Config.init_successive_orders_with_discrete_ordinates`
applies only to the legacy `SuccessiveOrders` source. It is not used by
`SuccessiveOrdersCpp`.

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
  sasktran2.Config.init_successive_orders_with_discrete_ordinates
  sasktran2.Config.num_stokes

```
