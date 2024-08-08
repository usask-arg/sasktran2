# Copied from numpy, allows for binary installation of openblas
set -xe

PROJECT_DIR="${1:-$PWD}"


if [[ $RUNNER_OS == "Windows" ]]; then
    # delvewheel is the equivalent of delocate/auditwheel for windows.
    python -m pip install delvewheel wheel

    # Need pkg-config to use scipy openblas
    choco install pkgconfiglite
fi


if [[ $RUNNER_OS == "Linux" ]]; then
    # Install same version of gfortran as the openblas-libs builds
    if [[ $PLATFORM == "macosx-arm64" ]]; then
        PLAT="arm64"
    fi
    source $PROJECT_DIR/ci/cibw/gfortran_utils.sh
    install_gfortran
    pip install "delocate==0.10.4"

    cd ~ && wget -nv https://github.com/xianyi/OpenBLAS/releases/download/v0.3.18/OpenBLAS-0.3.18.tar.gz && tar xf OpenBLAS-0.3.18.tar.gz
    mkdir ~/OpenBLAS-0.3.18/build && cd ~/OpenBLAS-0.3.18/build && cmake ../ -DBUILD_SHARED_LIBS:bool=OFF -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true -DLAPACKE:BOOL=on && cmake --build . && cmake --install .
fi

if [[ $RUNNER_OS == "macOS" ]]; then
    # Install same version of gfortran as the openblas-libs builds
    if [[ $PLATFORM == "macosx-arm64" ]]; then
        PLAT="arm64"
    fi
    source $PROJECT_DIR/ci/cibw/gfortran_utils.sh
    install_gfortran
    pip install "delocate==0.10.4"

    cd ~ && wget -nv https://github.com/xianyi/OpenBLAS/releases/download/v0.3.18/OpenBLAS-0.3.18.tar.gz && tar xf OpenBLAS-0.3.18.tar.gz
    mkdir ~/OpenBLAS-0.3.18/build && cd ~/OpenBLAS-0.3.18/build && cmake ../ -DBUILD_SHARED_LIBS:bool=OFF -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true -DLAPACKE:BOOL=on && cmake --build . && cmake --install .
fi
