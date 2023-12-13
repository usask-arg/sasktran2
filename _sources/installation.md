
(_installation)=
# Installation

The recommended way to install SASKTRAN2 is through `conda`

    conda install -c usask-arg -c conda-forge sasktran2

Wheels are also made available through `pip`

    pip install sasktran2


# Nightly Builds

The latest nightly version of SASKTRAN2 is made available through::

    conda install -c usask-arg-nightly -c conda-forge sasktran2

# Supported Platforms
|   | macOS Intel | macOS Apple Silicon | Windows 64bit | manylinux x86_64 |
|---------------|----|-----|-----|-----|
| Pip CPython 3.10  | ✅¹ | ✅¹  | ✅  | ✅  |
| Pip CPython 3.11  | ✅¹ | ✅¹  | ✅  | ✅  |
| Pip CPython 3.12  | ❌ | ❌  | ❌  | ❌  |
| conda-forge Py 3.10  | ✅ | ✅  | ✅  | ✅  |
| conda-forge Py 3.11  | ✅ | ✅  | ✅  | ✅  |
| conda-forge Py 3.12  | ✅ | ✅  | ✅  | ✅  |

<sup>¹ Open-MP Support is Disabled</sup><br>
