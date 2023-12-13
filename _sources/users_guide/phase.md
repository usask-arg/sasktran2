---
file_format: mystnb
---

(_users_phase)=
# Phase Function Specification
Generally in SASKTRAN2 the phase function is specified in terms of it's Legendre coefficient expansion.
The number of coefficients in the expansion determines the accuracy of the calculation, however a different
number of expansion terms are used in the single scatter and multiple scatter calculations.  The number
of expansion terms for each is controlled by {py:attr}`sasktran2.Config.num_singlescatter_moments` and
{py:attr}`sasktran2.Config.num_streams` respectively.

## Delta-m Scaling
SASKTRAN2 has the capability to perform Delta-m scaling, where strong forward peaks of the phase functions are
handled in an approximate fashion.  Generally the approximation performs very well for situations where you
are not directly observing the forward scatter peak of a strong forward scatterer such as a cloud.

Delta-m scaling is disabled by default in SASKTRAN2,
however it can be explicitly enabled with

```{code}
import sasktran2 as sk

config = sk.Config()

config.delta_m_scaling = True # False
```

As with most things, the performance of delta-m scaling depends heavily on the specific application.
The best way to determine if it is useful for your specific calculation is to compare your results
to a calculation with a large number of streams.
