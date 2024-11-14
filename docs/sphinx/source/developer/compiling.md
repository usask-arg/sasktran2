
(_dev_compiling)=
# Compiling from Source
SASKTRAN2 is primarily `c++` code, and thus local development requires the capability to compile the model.  Below are various
different ways you can install the dependencies to do this, for most users using `pixi` is the easiest option.

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

### Using pixi
The easiest way to setup a `sasktran2` development environment is through the package manager [pixi](https://pixi.sh/latest/).
Assuming you have `pixi` installed, you can create a development environment on any platform through

```
pixi install
```

Then the code can be compiled through

```
pixi run build
```

You can then verify your environment is correct by running the provided tests

```
pixi run test
```

### Using the conda development environnment
Alternatively you can create a `conda` environment directly from provided development environment files

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
