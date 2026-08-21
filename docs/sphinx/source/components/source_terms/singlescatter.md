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
towards the sun every time the single scatter source is required at a new point, i.e., the solar attenuation
is calculated "exactly". With solar refraction in a 2D atmosphere, a shared characteristic table supplies
the local bent direction and the optical depth along that characteristic is retraced at every point.

The other available single scatter source,

```{code-cell}
config.single_scatter_source = sk.SingleScatterSource.Table
```

differs in that the solar attenuation is pre-computed on an appropriate grid, and then interpolated whenever
it is requested. For a 2D atmosphere this is a three-dimensional table in altitude, solar zenith angle, and
off-plane azimuth. It stores incremental characteristic stencils rather than a cumulative optical-depth row
for every requested point. This can reduce accuracy in some cases, but substantially reduces setup time when
many lines of sight are requested or when the table is shared with successive orders.

Single scattering can be explicity disabled with

```{code-cell}
config.single_scatter_source = sk.SingleScatterSource.NoSource
```

## Additonal Notes

 - The exact source has native JVP and VJP implementations and can also materialize a complete Jacobian.
 - The 2D table source has native JVP and VJP implementations. To retain its memory advantage it does not materialize a complete Jacobian; use {py:meth}`sasktran2.Engine.linearize` and its `jvp`/`vjp` methods.
 - The 1D table source continues to use the materialized-Jacobian backend.
 - All single scatter sources support polarized calculations (`nstokes=3`)

## Relevant Configuration Options
```{eval-rst}
.. autosummary::

  sasktran2.Config.single_scatter_source
  sasktran2.Config.num_stokes

```
