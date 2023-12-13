---
file_format: mystnb
---

(_source_single_scatter)=
# Single Scatter Sources
Single scattering in the atmosphere is defined as the direct solar beam attenuated to a point in the atmosphere,
and then scattered into a specific direction, e.g., back into an instrument's line of sight.  We also consider
direct reflection of the Earth's surface when viewing in the nadir direction to be considered single scattering
for convenience.

Generally the recommended way to enable single scattering is with

```{code-cell}
import sasktran2 as sk

config = sk.Config()

config.single_scatter_source = sk.SingleScatterSource.Exact
```

here `Exact` refers to how the solar attenuation is calculated.  In the `Exact` mode, a ray is traced
towards the sun everytime the single scatter source is required at a new point, i.e., the solar attenuation
is calculated "exactly".

The other available single scatter source,

```{code-cell}
config.single_scatter_source = sk.SingleScatterSource.Table
```

differs in that the solar attenuation is pre-computed on an appropriate grid, and then interpolated whenever
it is requested.  This can reduce accuracy in some cases, but also offer computational efficiency improvements.
Generally the advantages of using `Table` over `Exact` are minimal, and we only recommended experimenting
with the table option if you have a situation where many lines of sight are requested (>100) and the majority
of the calculation time is suspected to be inside the single scatter source.

Single scattering can be explicity disabled with

```{code-cell}
config.single_scatter_source = sk.SingleScatterSource.NoSource
```

## Additonal Notes

 - All available single scatter sources are perfectly linearized, capable of producing weighting functions to machine precision
 - All single scatter sources support polarized calculations (`nstokes=3`)

## Relevant Configuration Options
```{eval-rst}
.. autosummary::

  sasktran2.Config.single_scatter_source
  sasktran2.Config.num_stokes

```
