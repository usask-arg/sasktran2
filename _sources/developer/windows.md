
(_dev_windows)=
# Local Windows Development
Here we describe several specific things we have noticed when trying to develop SASKTRAN2 on windows.

## Debug Builds on Windows
To build the code in debug mode on windows you also need a debug version of catch2.  You can run

```
conda activate sasktran2-dev-env
git clone git@github.com:catchorg/Catch2.git

cd catch2

cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=CONDA_ENV_ROOT
cmake --build build --config Debug --target Install
cmake --build build --config Release --target Install
```
where `CONDA_ENV_ROOT` is the path from the command

```
python -c "import os; print(os.environ['CONDA_PREFIX'] + '/Library')"`
```
