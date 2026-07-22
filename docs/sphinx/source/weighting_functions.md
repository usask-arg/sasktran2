---
file_format: mystnb
---

(_weighting_functions)=
# Weighting Functions

One of the key features of SASKTRAN2 is that it supports algorithmic calculation of weighting functions,
i.e. derivatives of the output parameters with respect to input parameters.  In essence, we are calculating

```{math}
    w(x) = \frac{\partial I}{\partial x}
```

where {math}`x` is an input parameter to the model, such as ozone VMR at the 10 km level, or the Lambertian
surface reflectance at 350 nm.

## How to Calculate Weighting Functions
For the most part, weighting functions are calculated automatically, completely opaquely to the user.
For example, if we set up a simple calculation,

```{code-cell} ipython3
import sasktran2 as sk
import numpy as np
import xarray as xr

config = sk.Config()

model_geometry = sk.Geometry1D(cos_sza=0.6,
                                solar_azimuth=0,
                                earth_radius_m=6372000,
                                altitude_grid_m=np.arange(0, 65001, 1000),
                                interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                geometry_type=sk.GeometryType.Spherical)

viewing_geo = sk.ViewingGeometry()

for alt in [10000, 20000, 30000, 40000]:
    ray = sk.TangentAltitudeSolar(tangent_altitude_m=alt,
                                    relative_azimuth=0,
                                    observer_altitude_m=200000,
                                    cos_sza=0.6)
    viewing_geo.add_ray(ray)

wavel = np.arange(280.0, 800.0, 0.1)

atmosphere = sk.Atmosphere(model_geometry, config, wavelengths_nm=wavel)

sk.climatology.us76.add_us76_standard_atmosphere(atmosphere)

atmosphere['rayleigh'] = sk.constituent.Rayleigh()

atmosphere['ozone'] = sk.climatology.mipas.constituent("O3", sk.optical.O3DBM())

atmosphere['no2'] = sk.climatology.mipas.constituent("NO2", sk.optical.NO2Vandaele())

engine = sk.Engine(config, model_geometry, viewing_geo)
```

and then inspect the output

```{code-cell} ipython3
radiance = engine.calculate_radiance(atmosphere)
print(radiance)
```

We see several variables in the output: [`wf_pressure_pa`, `wf_temperature_k`, `wf_ozone_vmr`, `wf_no2_vmr`].
These are the derivatives of the `radiance` field with respect to the input variables.
All of these variables have the same dimensions as the radiance, plus an additional dimension.  The additional
dimension is the input dimension, e.g. we specified ozone VMR on an altitude grid thus the input dimension is
this altitude grid.

We can take a look at one of the weighting functions,

```{code-cell} ipython3
radiance["wf_ozone_vmr"].isel(los=1, stokes=0).sel(wavelength=330, method="nearest").plot(y="ozone_altitude")
```

## Jacobian-vector and vector-Jacobian products

The same derivatives are available as a local linear model of the radiance:

    linearization = engine.linearize(atmosphere)

    # A labeled perturbation. Parameters not included here are held fixed.
    tangent = linearization.tangent_template[["ozone_vmr"]]
    tangent["ozone_vmr"].data[:] = ozone_perturbation
    radiance_perturbation = linearization.jvp(tangent)

    # Propagate a scalar objective's radiance sensitivity back to every input.
    parameter_gradients = linearization.vjp(radiance_sensitivity)

    # Retrievals can request only their active parameters. Unrequested
    # derivative mappings are then excluded from the VJP calculation.
    retrieval_gradients = linearization.vjp(
        radiance_sensitivity,
        parameters=("ozone_vmr",),
    )

The `linearization.value` property is the radiance at the linearization point.
JVP and VJP operations contract the native derivative rows as they are produced,
so they do not allocate the complete mapped Jacobian.

The read-only `linearization.backends` mapping reports whether each product
uses a `LinearizationBackend.Native`, `StreamingJacobian`, or
`MaterializedJacobian` implementation. A streaming-Jacobian product avoids
allocating the complete mapped Jacobian, but repeats the full derivative
propagation for every call. Iterative solvers may therefore prefer a native
product backend, or explicitly materialize the Jacobian once when only a
streaming backend is available.

Accessing `linearization.jacobian` materializes and caches one labeled block per
semantic parameter (for example, `ozone_vmr` rather than `wf_ozone_vmr`). The
`tangent_template` property describes the shape, dimensions, and coordinates of
every parameter without materializing the Jacobian. Flux and diagnostic outputs
are not part of this interface.

A linearization keeps a reference to the atmosphere; it does not copy it. The
atmosphere's revision is recorded when the linearization is created. Rebuilding
a constituent-backed atmosphere advances this revision automatically. When
modifying a borrowed raw-storage NumPy view, call `atmosphere.mark_changed()`
after the modification. Once the revision changes, new JVPs, VJPs, and a
first-time Jacobian materialization raise `StaleLinearizationError`; already
cached `value` and `jacobian` objects remain readable. Create a new
linearization to evaluate derivatives at the changed state.

Trust-region solvers may need to evaluate a trial atmosphere while retaining
the linearization at the currently accepted state. Use two distinct
`Atmosphere` instances for these states. Linearisations backed by distinct
atmospheres can coexist and be evaluated sequentially with the same `Engine`;
rebuilding one atmosphere does not invalidate a linearization backed by the
other. Concurrent calls on a shared engine are not supported.

Weighting functions configured in log-radiance space are converted back to
radiance derivatives by this interface, so `jvp`, `vjp`, and `jacobian` always
describe the derivative of `linearization.value`. The legacy
`calculate_radiance()` output is unchanged.

## The Details

Weighting functions in SASKTRAN2 are calculated analytically, evaluating the derivative of the radiative
transfer equation alongside the radiative transfer equation itself in code.  This is similar to the
"forward mode" of auto-differentiation schemes, however all derivatives in SASKTRAN2 are hand-crafted
rather than automatically generated for maximum computational efficiency.  Generally you can expect a
calculation including full derivatives to be somewhere in the range of 2-5x slower than the calculation
without derivatives.

The accuracy of the weighting functions depends on which source terms are included in the calculation.
With the exception of the `sk.MultipleScatterSource.SuccessiveOrders` source, all other sources are
"perfectly linearized", which means that the weighting function calculation is accurate to relative machine
precision.  Typically we estimate this to be around 6 decimal places for most applications.
