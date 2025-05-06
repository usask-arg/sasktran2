
(_installation)=
# Installation

The recommended way to install SASKTRAN2 is through `pip/uv`

    pip install sasktran2

Packages are also provided on `conda-forge`

    conda install conda-forge::sasktran2



```{note}
Python 3.10 support for conda-forge SASKTRAN2 was dropped in 2025.01.0
```


# Supported Platforms
|   | macOS Intel | macOS Apple Silicon | Windows 64bit | manylinux x86_64 | manylinux aarch64 |
|---------------|----|-----|-----|-----|-----|
| Pypi CPython 3.10  | ✅ | ✅  | ✅  | ✅  | ✅  |
| Pypi CPython 3.11  | ✅ | ✅  | ✅  | ✅  | ✅  |
| Pypi CPython 3.12  | ✅ | ✅  | ✅  | ✅  | ✅  |
| Pypi CPython 3.13  | ✅ | ✅  | ✅  | ✅  | ✅  |
| conda-forge Py 3.11  | ✅ | ✅  | ✅  | ✅  | ✅  |
| conda-forge Py 3.12  | ✅ | ✅  | ✅  | ✅  | ✅  |
| conda-forge Py 3.13  | ✅ | ✅  | ✅  | ✅  | ✅  |

## Mac OMP Errors
If you encounter the error `OMP: Error #15: Initializing libomp.dylib, but found libomp.dylib already initialized.` it is most likely
that you are trying to use the `sasktran2` wheel package (through `pip`) inside a conda environment.  The best solution
is to instead install the conda-forge package.

Alternatively, if that is not possible, or you are getting the error outside of a conda environment, you can try adding

```{code}
import os

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
```

To the top of your script.
