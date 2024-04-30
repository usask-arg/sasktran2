---
sd_hide_title: true
---

# SASKTRAN2

::::{grid}
:reverse:
:gutter: 2 1 1 1
:margin: 4 4 1 1

:::{grid-item}
:columns: 4

```{image} ./_static/sasktran-dark.svg
:width: 250px
:class: sd-m-auto transparent-image only-light
:name: landing-page-logo-light
```

```{image} ./_static/sasktran-light.svg
:width: 250px
:class: sd-m-auto transparent-image only-dark
:name: landing-page-logo-dark
```
:::

:::{grid-item}
:columns: 8
:class: sd-fs-3

The SASKTRAN Radiative Transfer Framework (Version 2)

:::

::::

SASKTRAN2 is an atmospheric radiative transfer model developed at the University of Saskatchewan. It is a full rewrite of the original [SASKTRAN Model](https://github.com/usask-arg/sasktran)
providing

 - A simpler, more powerful, user interface
 - Full capability to calculate weighting functions, i.e. derivatives of output radiative quantities with respect to input atmospheric properties
 - Greatly improved efficiency for hyperspectral calculations

SASKTRAN2 is still in active development, and things are subject to change quickly, however it is already usable for many applications.
After [Installing](_installation) the model, we recommend starting with the [Quick Start Guide](_quickstart) to get a feel for the core concepts of the interface.

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item-card} {octicon}`paper-airplane;1.5em;sd-mr-1` Installation
:link: _installation
:link-type: ref

SASKTRAN2 is available both as a `conda` and `pip` package on most platforms.

+++
[Learn more »](installation)
:::

:::{grid-item-card} {octicon}`plug;1.5em;sd-mr-1` Quick Start Guide
:link: _quickstart
:link-type: ref

The quick start guide will help you set up your first radiative transfer calculation.

+++
[Learn more »](quickstart)
:::

:::{grid-item-card} {octicon}`tools;1.5em;sd-mr-1` User's Guide
:link: _weighting_functions
:link-type: ref

The user's guide demonstrates SASKTRAN2's features through example, including powerful features
such as calculating derivatives of output values with respect to input values.

+++
[Learn more »](users_guide)
:::

:::{grid-item-card} {octicon}`code;1.5em;sd-mr-1` API Reference
:link: _api_reference
:link-type: ref

A full reference to SASKTRAN2's API.  This section assumes you are already
familiar with the core concepts of the model.

+++
[Learn more »](api_reference)
:::


::::

## License
SASKTRAN2 is made available under the MIT license (see [License](https://github.com/usask-arg/sasktran2/blob/main/license.md))

We request that users of the model contact the authors before publishing results using SASKTRAN2, and that
the following publications are acknowledged:

Zawada, D. J., Dueck, S. R., Rieger, L. A., Bourassa, A. E., Lloyd, N. D., and Degenstein, D. A.: High-resolution and Monte Carlo additions to the SASKTRAN radiative transfer model, Atmos. Meas. Tech., 8, 2609-2623, https://doi.org/10.5194/amt-8-2609-2015, 2015.

Bourassa, A. E., Degenstein, D. A., and Llewellyn, E. J.: SASKTRAN: A Spherical Geometry Radiative Transfer Code for Efficient Estimation of Limb Scattered Sunlight, J Quant Spectrosc Radiat Trans, Volume 109, Issue 1, 52-73, https://doi.org/10.1016/j.jqsrt.2007.07.007, 2008.



```{toctree}
:maxdepth: 4
:hidden:
:caption: For Users

installation
quickstart
users_guide
api_reference
changelog
```

```{toctree}
:maxdepth: 4
:hidden:
:caption: Extending SASKTRAN2

extending/constituent.md
```

```{toctree}
:maxdepth: 4
:hidden:
:caption: Developer Documentation

developer/compiling.md
developer/contributing.md
developer/codespaces.md
```
