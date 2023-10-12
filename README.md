# SASKTRAN2
The SASKTRAN radiative transfer framework is a radiative transfer tool developed at the University of Saskatchewan. Originally designed for use with the OSIRIS instrument (https://research-groups.usask.ca/osiris/) it has since evolved to be applicable to a large variety of applications. SASKTRAN is a full framework and not just a radiative transfer model, as such it contains databases or interfaces to standard climatologies and species optical properties.

SASKTRAN2 is a full re-implementation of the original SASKTRAN framework with large computational efficiency
improvements, full linearizations of atmospheric input properties, and an improved Python interface.

## Installation
The preferred method to install SASKTRAN2 is through the pre-compiled Conda package

```
conda install -c usask-arg -c conda-forge sasktran2
```
these packages are made available for Python versions 3.8, 3.9, 3.10, 3.11, 3.12 on Windows/Linux/Mac platforms.
For Mac, both x86_64 and Arm packages are available.

Wheels are also built for the same platforms and can be installed through,
```
pip install sasktran2
```

SASKTRAN2 can also be built directly from source,
```
pip install .
```

When building from source it is required that a Blas/LAPACK implementation is findable by CMake.

## Usage
Documentation can be found at https://usask-arg.github.io/sasktran2/

## License
SASKTRAN2 is made available under the MIT license subject to the Commons Clause condition (see license.md). Effectively this is a MIT license restricted for academic and educational use, for commercial use please contact the package authors. Commerical level support may also be available for specific applications.

## Acknowledgement
We request that users of the model contact the authors before publishing results using SASKTRAN, and that the following publications are acknowledged:

Zawada, D. J., Dueck, S. R., Rieger, L. A., Bourassa, A. E., Lloyd, N. D., and Degenstein, D. A.: High-resolution and Monte Carlo additions to the SASKTRAN radiative transfer model, Atmos. Meas. Tech., 8, 2609-2623, https://doi.org/10.5194/amt-8-2609-2015, 2015.

Bourassa, A. E., Degenstein, D. A., and Llewellyn, E. J.: SASKTRAN: A Spherical Geometry Radiative Transfer Code for Efficient Estimation of Limb Scattered Sunlight, J Quant Spectrosc Radiat Trans, Volume 109, Issue 1, 52-73, https://doi.org/10.1016/j.jqsrt.2007.07.007, 2008.