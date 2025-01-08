
(_users_guide)=
# Users Guide


[Configuring the Model and Setting up Your Calculation](required_components.md)
: Information on how to set up the model for your specific calculation

[Solar Irradiance](users_guide/solar_irradiance.md)
: Absolute vs relative radiances, and how to include the solar irradiance

[Atmospheric and Surface Emissions](users_guide/emissions.md)
: Enabling emission source and how to include thermal emissions

[Calculating Weighting Functions](weighting_functions.md)
: How to calculate derivatives (weighting functions) of output quantities with respect to input quantities

[Polarized Calculations](users_guide/polarization.md)
: Enabling calculation of the full Stokes Vector

[Setting the Surface BRDF](users_guide/brdf.md)
: How to use the built in surface reflectance models

[Including Mie Scatterers](users_guide/mie.md)
: Adding Mie scattering particles to your atmospheric state

[Ray Tracing and Refraction](users_guide/refraction.md)
: Including refractive effects in the ray-tracing

[Built in Geodesy Tools](users_guide/geodetic.md)
: Methods to work with datums such as WGS84

[Input Validation](users_guide/input_validation.md)
: Information on how user input values are automatically validated

[More Examples](examples/examples.rst)
: Miscallaneous examples of model usage


```{toctree}
:maxdepth: 4
:hidden:

required_components
weighting_functions
users_guide/solar_irradiance.md
users_guide/emissions.md
users_guide/mie.md
users_guide/refraction.md
users_guide/polarization.md
users_guide/brdf.md
users_guide/performance.md
users_guide/phase.md
users_guide/geodetic.md
users_guide/input_validation.md
examples/examples
```
