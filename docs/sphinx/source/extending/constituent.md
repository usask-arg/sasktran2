
(_extend_constituent)=
# Adding Constituents
Put simply, a constituent in SASKTRAN2 is simply an object that modifies the atmospheric state in some way.
This means that it is somewhat straightforward for users to implement their own constituents in the case that
one of SASKTRAN2's built in ones does not fit their needs.

In this guide, we will implement a new constituent that acts as a spectrally flat absorber, and has a constant VMR
in altitude.

## Initializing the Constituent
To start, we create an object that inherits from `sasktran2.constituent.base.Constituent`, the abstract base class
that defines the constituent interface.

```python
import sasktran2 as sk


class SpectrallyFlatAbsorber(sk.constituent.base.Constituent):
    def __init__(self, vmr: float):
        self._vmr = vmr

```

## Modifying the Atmospheric State
The main function of the constituent is to modify the atmosphere, this is done through the abstract method
`add_to_atmosphere` which is called before construction.  For our constituent we have to modify the
total atmospheric extinction,

```python
    def add_to_atmosphere(self, atmo: sk.Atmosphere):
        number_density = sk.optical.pressure_temperature_to_numberdensity(
            atmo.pressure_pa, atmo.temperature_k
        )

        atmo.storage.total_extinction += (
            * (number_density * self._vmr)[:, np.newaxis]
        )
```

Here we have done a couple of things, first we have taken pressure and temperature from the atmospheric state in order
to get the background number density.

## Enabling Weighting Function Calculations
The last piece of the constituent is telling the model how to propagate internal derivatives (in this case with respect to extinction,
and single scatter albedo) to derivatives of our input quantities (the constant vmr value).  This is done through the
method `register_derivative`.

Propagating the derivative is tricky to get correct, and so we recommend taking some care and really understand what is going on.
SASKTRAN2 internally calculates derivatives with respect to the fundamental input parameters, in this case the input parameters that are
affected by our constituent are the total extinction, {math}`k`, and the single scatter albedo, {math}`\omega`.  You may wonder why
the single scatter albedo is affected since our `add_to_atmosphere` function did not modify it directly, but since single scatter albedo
at location index {math}`i` is defined as,

```{math}
\omega_i = \frac{\sum_j \omega_{i,j} k_{i,j}}{\sum_j k_{i,j}}
```
where the sum {math}`j` is over all constituents, there is an implicit dependence.

SASKTRAN2 calculates the quantities,
```{math}
\frac{\partial I}{\partial k_i}, \frac{\partial I}{\partial \omega_i},
```
and our desired derivative can be written as,
```{math}
\frac{\partial I}{\partial x} = \sum_i \frac{\partial k_i}{\partial x} \frac{\partial I}{\partial k_i} + \sum_i \frac{\partial \omega_i}{\partial x} \frac{\partial I}{\partial \omega_i},
```
which means we need the quantities {math}`\partial k_i / \partial x` and {math}`\partial \omega_i / \partial x`.  These are the quantities our constituent must supply
in order to propagate the derivative correctly.

The extinction derivative is straightforward, since {math}`k_i = \sum_j k_{i, j}`,

```{math}
\frac{\partial k_i}{\partial x} = n_i,
```
where {math}`n_i` is our background number density.  The single scattering albedo derivative is a little more involved, but from a simple chain rule

```{math}
\frac{\partial \omega_i}{\partial x} = \frac{-n_i \omega_i}{k_i}
```
where we have used the fact that the single scatter albedo for our absorber is 0.

Now we can put this all together,

```python
   def register_derivative(self, atmo: sk.Atmosphere, name: str):
        number_density = sk.optical.pressure_temperature_to_numberdensity(
            atmo.pressure_pa, atmo.temperature_k
        )
        interp_matrix = np.ones_like(atmo.model_geometry.altitudes())
        derivs = {}

        derivs["vmr"] = sk.atmosphere.InterpolatedDerivativeMapping(
           sk.atmosphere.NativeGridDerivative(
                d_extinction=number_density[:, np.newaxis],
                d_ssa=(-atmo.storage.ssa)
                / atmo.storage.total_extinction
                * number_density[:, np.newaxis],
            ),
            interpolating_matrix=interp_matrix,
            interp_dim="altitude",
            result_dim=f"{name}_altitude",
        )

        return derivs
```

Note that here we have done two things.  The first is the construction of a `sk.atmosphere.NativeGridDerivative` which applies our
derivative factors, but does not sum them together.  The `sk.atmosphere.InterpolatedDerivativeMapping` takes in our summation
vector (in this case, all ones), and tells the model to perform the sum.
