---
file_format: mystnb
---

(_users_input_validation)=
# Input validation
```{code-cell} ipython3
import sasktran2 as sk
```

By default, input parameters to SASKTRAN2 are validated at runtime, which can help to catch common mistakes.
For example, let's try to setup a calculation where the geometry object does not contain an ascending altitude grid,

```{code-cell} ipython3
import sasktran2 as sk
import numpy as np

config = sk.Config()

model_geometry = sk.Geometry1D(cos_sza=0.6,
                                solar_azimuth=0,
                                earth_radius_m=6372000,
                                altitude_grid_m=[0, 1000, 500, 2000],
                                interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                geometry_type=sk.GeometryType.Spherical)

```

Input validation can cause a slight runtime performance overhead.
In cases where you are confident that the input parameters are valid, the expensive components of input
validation can be disabled by setting

```{code-cell} ipython3
config.input_validation_mode = sk.InputValidationMode.Disabled
```
