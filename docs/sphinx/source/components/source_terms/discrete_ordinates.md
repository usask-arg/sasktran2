---
file_format: mystnb
---

(_source_discrete_ordinates)=
# The Discrete Ordinates Source
The discrete ordinates method is an implementation of multiple scattering that is enabled with

```{code-cell}
import sasktran2 as sk

config = sk.Config()

config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
```

## Advantages

 - "exact" solution in a plane-parallel medium of stacked homogenous layers
 - Correct solutions in layers of large optical depth without sub-layering
 - Computationally efficient for most cases

## Disadvantages

 - May not be suitable for some limb-viewing applications because of errors in the plane-parallel assumption
 - Only suitable for 1-D atmospheres
 - May not be suitable for twilight viewing conditions (solar zenith angle > 85)

## Corrections in a Spherical Atmosphere
As noted, the discrete ordinates technique is a solution to the radiative transfer equation in a plane-parallel medium,
however several corrections are implemented to improve the solution in spherical atmospheres.

The first is the pesudo-spherical approximation, which accounts for sphericity when calculating the solar attenuation
factor.  This correction has very little computational overhead and is thus always enabled when appropriate for the
given model geometry.

The second correction is more involved and attempts to correct for line of sight sphericity effects.  In a spherical
atmosphere, the solar zenith angle changes along any given line of sight, but the discrete ordinates solution is only
valid for a single solar zenith angle.  To account for this multiple discrete ordinates calculations are performed at
different solar zenith angles.  The resulting sources are then interpolated in cosine of solar zenith angle.
The number of such calculations is controlled with the {py:attr}`sasktran2.Config.num_sza` attribute.  The required value
of this attribute to obtain an accurate solution is very application dependent, we recommend increasing it until the solution
stops changing.  A rough rule of thumb is that for all ground viewing applications, 1-2, is usually sufficient.  For limb viewing
applications usually somewhere in the range 1-7 is good enough.


## Notes

 - The discrete ordinates source is perfectly linearized, i.e., weighting functions are calculated to machine precision

## Relevant Configuration Options

```{eval-rst}
.. autosummary::

  sasktran2.Config.multiple_scatter_source
  sasktran2.Config.num_streams
  sasktran2.Config.num_sza
  sasktran2.Config.num_stokes
  sasktran2.Config.do_backprop

```
