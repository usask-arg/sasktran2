
(_dev_vscode)=
# Visual Studio Code Configuration
First run

`python -c "import os; print(os.environ['CONDA_PREFIX'])"`

which will print out a directory.  Make note of this directory.

## Settings
-----------
Open the folder in vscode. Do ctrl+shift+P to open the settings search menu, and select "Preferences: Open Workspace Settings (JSON). Add the following to the file:

```
    "cmake.configureArgs": [
        "-DCMAKE_PREFIX_PATH=CONDA_ENV_ROOT",
        "-DSKTRAN_BLAS_VENDOR=OpenBLAS",
        "-DCOPY_TESTING_DLLS=ON",
        "-DBUILD_TESTS=ON"
    ]
```
where `CONDA_ENV_ROOT` is the path found from above.  NOTE on windows it is necessary to append `/Library` to the path.


## Useful Extensions

### C/C++ Extension Pack
VS Code will likely install this by default. It includes almost mandatory extensions to develop C/C++ code

### C++ TestMate
Adds support for the Catch2 tests used within the project, so you can run and debug individual tests.

### Python and Pylance
Almost mandatory extensions for Python development, VS Code likely will ask to automatically install them when you open a Python file.

### Formatting Tools
Several extensions can help with automatic formatting:

 - `autoDocstring - Python Docstring Generator`: Helps you create Python docstrings faster
 - `Black Formatter`: Automatically formats Python code with the style the project uses
 - `Clang-Format`: Automatically formats c++ code with the style the project uses
