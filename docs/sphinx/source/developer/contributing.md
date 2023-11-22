
(_dev_contributing)=
# Contributing

## Setting up the repository

To contribute to SASKTRAN2, start by creating a fork of the main repository.  Go to https://github.com/usask-arg/sasktran2 and press the fork button on the top right.  This will create a copy of the repository in your local Github namespace.

Your new fork will open and then you can clone the forked repository.  The command will look something like

```
git clone git@github.com:YOURUSERNAME/sasktran2.git
```

## Developing a feature

To create a new feature, you should work on a feature branch.  Ideally each branch is isolated to a single feature.  You can make a new branch with

```
git checkout -b shiny-new-feature
```

which will then create a new branch with the name `shiny-new-feature`.  The branch only exists on your local fork.
Don't worry about how many commits you make to the branch or how messy they are, when it is time to merge the branch into
the main repository all of the commits will be squashed into a single one and you will have the opportunity to write a new
description.

## Creating a pull request
Once you have made a commit to your branch and pushed the changes to your fork, you are ready to make a pull request.  This can
be done from the `pull requests` page on the main repository.  If your pull request is still a work in progress, prepend the title with `WIP:`

In addition to the feature, we require that all pull requests contain tests for the feature (if applicable), as well as updates to the documentation (if applicable)


## Code Style
For Python code we follow the black code style (https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html), and for C++ code we do not directly follow a published style.  C++ code is auto-formatted through `clang-format` so typically you do not have to worry about it.

For both we recommend using the provided pre-commit hook to automatically both format and check your code.  These can be installed with

```
pip install pre-commit
```

You can run the formatters manually by typing

```
pre-commit run -a
```

You will see output that looks like

```
check for added large files..............................................Passed
check for case conflicts.................................................Passed
check for merge conflicts................................................Passed
check for broken symlinks............................(no files to check)Skipped
check yaml...............................................................Passed
debug statements (python)................................................Passed
fix end of files.........................................................Passed
mixed line ending........................................................Passed
fix requirements.txt.................................(no files to check)Skipped
trim trailing whitespace.................................................Passed
black....................................................................Passed
ruff.....................................................................Passed
Tabs remover.............................................................Passed
cmake-format.............................................................Passed
clang-format.............................................................Passed
```

If you want, you can then run

```
pre-commit install
```

which will set up the pre commit hooks to run automatically every-time you commit to the repository.

We highly recommend using the pre commit hooks. On every pull request these checks are automatically ran, and the code will not be merged in
if any fail.  Using the pre commit hooks locally saves everyone time.
