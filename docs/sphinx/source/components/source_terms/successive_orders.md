---
file_format: mystnb
---

(_source_succesive_orders)=
# The Successive Orders of Scattering Source
The successive orders of scattering source is an implementation of the multiple scattering source using the successive orders
of scattering method.  It can be enabled with

```{code-cell}
import sasktran2 as sk

config = sk.Config()

config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
```

## Advantages

 - The only available multiple scattering source that fully accounts for sphericity of the atmosphere

## Disadvantages

 - Weighting function calculations are only approximate
 - Can use a large amount of RAM for large, polarized, calculations
 - May require sub-layering in optically thick areas of the atmosphere


## Initialization Options
The successive orders method is an iterative method which repeatedly refines the solution.  When the method is initialized
with the single scatter source function, each iteration physically represents including an additional "order of scatter" into
the solution.  The number of these iterations can be configured with {py:attr}`sasktran2.Config.num_successive_orders_iterations`.

Alternatively, the successive orders method can be initialized with the source function calculated with a discrete ordinates solution.
This has two advantages: the first is that fewer iterations are required to reach convergence, typically only 1--2 instead of sometimes 50;
the second advantage is that the discrete ordinates method can provide a reasonable approximation for the weighting functions, whereas
full weighting function calculations in a successive orders model are extremely computationally expensive.
This option can be enabled with the {py:attr}`sasktran2.Config.init_successive_orders_with_discrete_ordinates` option.  If it is set then
all of the standard discrete ordinates options apply to the source used to initialize the solution.

## Notes

 - The succesive orders of scattering source is not fully linearized, which means that weighting functions returned back will only be approximate.

## Relevant Configuration Options

```{eval-rst}
.. autosummary::

  sasktran2.Config.multiple_scatter_source
  sasktran2.Config.num_sza
  sasktran2.Config.num_successive_orders_iterations
  sasktran2.Config.init_successive_orders_with_discrete_ordinates
  sasktran2.Config.num_stokes

```
