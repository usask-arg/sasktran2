---
file_format: mystnb
---

(_sources)=
# Source Terms
In SASKTRAN2, source terms are entirely determined through configuration options, using the
{py:class}`sasktran2.Config` object

```{code-cell}
import sasktran2 as sk

config = sk.Config()
```


## Single Scatter Source Terms
Single scatter source terms are enabled by default, however you can also explicitly enable them with

```{code-cell}
config.single_scatter_source = sk.SingleScatterSource.Exact
```

To disable the single scatter source you can set

```{code-cell}
config.single_scatter_source = sk.SingleScatterSource.NoSource
```

## Multiple Scatter Source Terms
SASKTRAN2 contains two options to compute the multiple scatter source, the discrete ordinates technique,

```{code-cell}
config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
```

and then successive orders of scattering method

```{code-cell}
config.multiple_scatter_source = sk.MultipleScatterSource.SuccessiveOrders
```

The default is to disable multiple scattering, which can be explicity set through

```{code-cell}
config.multiple_scatter_source = sk.MultipleScatterSource.NoSource
```

## Occultation Source Terms
SASKTRAN2 can include an "occultation" source, where a solar source term is added at the end of the line of sight.
By setting

```{code-cell}
config.occultation_source = sk.OccultationSource.Standard
```

A source of identically 1 is added to the end of each line of sight.  The default is to include no
occultation sources, which can be explicity set through

```{code-cell}
config.occultation_source = sk.OccultationSource.NoSource
```

## Emission Source Terms
SASKTRAN2 can include emission sources along the sight and from the surface for lines of sight that intersect
the ground, by setting

```{code-cell}
config.emission_source = sk.EmissionSource.Standard
```

The default is to include no emission sources, which can be explicitly set through

```{code-cell}
config.emission_source = sk.EmissionSource.NoSource
```

Note that for this option to have any effect, emissions must be included in the model atmosphere, for example,
by using the {py:class}`sasktran2.constituent.ThermalEmission` and
{py:class}`sasktran2.constituent.SurfaceThermalEmission` constituents.

## Detailed Descriptions

```{eval-rst}
.. toctree::
   :maxdepth: 2

   source_terms/singlescatter
   source_terms/discrete_ordinates
   source_terms/occultation
   source_terms/successive_orders
   source_terms/emission
```
