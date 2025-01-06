---
file_format: mystnb
---

(_source_emission)=
# Emission Sources
Emission sources are enabled with

```{code-cell}
import sasktran2 as sk

config = sk.Config()

config.emission_source = sk.EmissionSource.Standard
```

Note that for this option to have any effect, emissions must be included in the model atmosphere, for example,
by using the {py:class}`sasktran2.constituent.ThermalEmission` and
{py:class}`sasktran2.constituent.SurfaceThermalEmission` constituents.

## Notes

Currently, emissions can not be included in the multiple scattering source.
They are only considered along the line of sight and from the surface when
the line of sight intersects the ground.

## Relevant Configuration Options

```{eval-rst}
.. autosummary::

  sasktran2.Config.emission_source

```
