# Copied from numpy, allows for binary installation of openblas

set -xe

PROJECT_DIR="$1"
PLATFORM=$(PYTHONPATH=tools python -c "import openblas_support; print(openblas_support.get_plat())")

# Install Openblas
if [[ $RUNNER_OS == "Linux" || $RUNNER_OS == "macOS" ]] ; then
    basedir=$(python tools/openblas_support.py)
    if [[ $RUNNER_OS == "macOS" && $PLATFORM == "macosx-arm64" ]]; then
        # /usr/local/lib doesn't exist on cirrus-ci runners
        sudo mkdir -p /usr/local/lib /usr/local/include /usr/local/lib/cmake/openblas
        sudo mkdir -p /opt/arm64-builds/lib /opt/arm64-builds/include
        sudo chown -R $USER /opt/arm64-builds
        cp -r $basedir/lib/* /opt/arm64-builds/lib
        cp $basedir/include/* /opt/arm64-builds/include
        sudo cp -r $basedir/lib/* /usr/local/lib
        sudo cp $basedir/include/* /usr/local/include
    else
        cp -r $basedir/lib/* /usr/local/lib
        cp $basedir/include/* /usr/local/include
    fi
elif [[ $RUNNER_OS == "Windows" ]]; then
    # delvewheel is the equivalent of delocate/auditwheel for windows.
    python -m pip install delvewheel
    python -m pip install wheel

    # make the DLL available for tools/wheels/repair_windows.sh. If you change
    # this location you need to alter that script.
    mkdir -p /c/opt/openblas/openblas_dll

    mkdir -p /c/opt/32/lib/pkgconfig
    mkdir -p /c/opt/64/lib/pkgconfig
    target=$(python -c "import tools.openblas_support as obs; plat=obs.get_plat(); ilp64=obs.get_ilp64(); target=f'openblas_{plat}.zip'; obs.download_openblas(target, plat, ilp64);print(target)")
    if [[ $PLATFORM == 'win-32' ]]; then
        # 32-bit openBLAS
        # Download 32 bit openBLAS and put it into c/opt/32/lib
        unzip -o -d /c/opt/ $target
        cp /c/opt/32/bin/*.dll /c/opt/openblas/openblas_dll
    else
        # 64-bit openBLAS
        unzip -o -d /c/opt/ $target
        if [[ -f /c/opt/64/lib/pkgconfig/openblas64.pc ]]; then
            # As of v0.3.23, the 64-bit interface has a openblas64.pc file,
            # but this is wrong. It should be openblas.pc
            cp /c/opt/64/lib/pkgconfig/openblas{64,}.pc
        fi
        cp /c/opt/64/bin/*.dll /c/opt/openblas/openblas_dll
    fi
fi

if [[ $RUNNER_OS == "macOS" ]]; then
    # Install same version of gfortran as the openblas-libs builds
    if [[ $PLATFORM == "macosx-arm64" ]]; then
        PLAT="arm64"
    fi
    source $PROJECT_DIR/ci/cibw/gfortran_utils.sh
    install_gfortran
    pip install "delocate==0.10.4"

    # Also install libomp
    brew install libomp
    if [[ $PLATFORM == "macosx-arm64" ]]; then
        sudo cp -r /opt/homebrew/opt/libomp/* /opt/homebrew/
    else
        sudo cp -r /usr/local/opt/libomp/* /usr/local/
    fi
fi
