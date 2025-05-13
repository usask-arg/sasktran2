# SASKTRAN
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/sasktran2/badges/version.svg)](https://anaconda.org/conda-forge/sasktran2)
[![Available on pypi](https://img.shields.io/pypi/v/sasktran2.svg)](https://pypi.python.org/pypi/sasktran2/)
[![Documentation Status](https://readthedocs.org/projects/sasktran2/badge/?version=latest)](https://sasktran2.readthedocs.io/en/latest/?badge=latest)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/usask-arg/sasktran2/main.svg)](https://results.pre-commit.ci/latest/github/usask-arg/sasktran2/main)



The SASKTRAN radiative transfer framework is a radiative transfer tool developed at the University of Saskatchewan. Originally designed for use with the OSIRIS instrument (https://research-groups.usask.ca/osiris/) it has since evolved to be applicable to a large variety of applications. SASKTRAN is a full framework and not just a radiative transfer model, as such it contains databases or interfaces to standard climatologies and species optical properties.

SASKTRAN2 is a full re-implementation of the original SASKTRAN framework with large computational efficiency
improvements, full linearizations of atmospheric input properties, and an improved Python interface.

## Installation
The preferred method to install SASKTRAN2 is through the pre-compiled Conda package

```
conda install -c conda-forge sasktran2
```
these packages are made available for Python versions 3.10, 3.11, 3.12, 3.13 on Windows/Linux/Mac platforms.
For Mac, both x86_64 and Arm packages are available.
For Linux, arm/ppc are also supported.

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
Documentation can be found at https://sasktran2.readthedocs.io/

## License
SASKTRAN2 is made available under the MIT license.

## Acknowledgement
We request that users of the model contact the authors before publishing results using SASKTRAN, and that the following publications are acknowledged:

Zawada, D. J., Dueck, S. R., Rieger, L. A., Bourassa, A. E., Lloyd, N. D., and Degenstein, D. A.: High-resolution and Monte Carlo additions to the SASKTRAN radiative transfer model, Atmos. Meas. Tech., 8, 2609-2623, https://doi.org/10.5194/amt-8-2609-2015, 2015.

Bourassa, A. E., Degenstein, D. A., and Llewellyn, E. J.: SASKTRAN: A Spherical Geometry Radiative Transfer Code for Efficient Estimation of Limb Scattered Sunlight, J Quant Spectrosc Radiat Trans, Volume 109, Issue 1, 52-73, https://doi.org/10.1016/j.jqsrt.2007.07.007, 2008.
