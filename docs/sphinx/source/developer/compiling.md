
(_dev_compiling)=
# Compiling from Source
In order to develop SASKTRAN2 locally you will need

- cmake
- C++ compiler
- A BLAS/Lapack Library
- Optionally Eigen3, spdlog, and Catch2

If you know how to install all of these things on your platform of choice, great, do that and it should work. If you don't know
where to begin, read on.

```{toctree}
:maxdepth: 4
:hidden:

windows.md
ide/vscode.md
```

## Local Installation
```{note}
The following assumes you have cloned the repository https://github.com/usask-arg/sasktran2 and have a terminal open
from inside the folder.
```
The easiest way to set up a local development environment for SASKTRAN2 is to use the provided `conda` environment
files

On windows,
```
conda create -f conda/dev_env_windows.yml
```

And on any other platform,
```
conda create -f conda/dev_env.yml
```

This will create a `conda` environment that you can activate with
```
conda activate sasktran2-dev-env
```

If you then run
```
pip install -e .
```
this will compile the code and install the package locally.

# Local Development
Local development for SASKTRAN2 can be done on a wide varity of tools and platforms, and providing
instructions on how to set up your local development environment is a work in progress.  The following
pages may be useful depending on your exact environment

[Visual Studio Code Configuration](ide/vscode.md)
: Instructions on how to set up visual studio code for local development

[Windows Specific Gotchas](windows.md)
: Oddities we have noticed while trying to develop the model on Windows
