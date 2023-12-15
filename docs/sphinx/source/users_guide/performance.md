---
file_format: mystnb
---

(_users_performance)=
# Performance Tips

```{code-cell}
import sasktran2 as sk

config = sk.Config()
```

## Multi-Threading Performance
By default SASKTRAN2 runs on a single thread. To change this behaviour you can set

```{code}
config.num_threads = 8
```
which will let the model run on 8 threads.  For maximum performance
it is reccomended to set {py:attr}`config.num_threads` to the number of `physical`
cores on your machine, which may be less than the number of logical cores.

SASKTRAN2 is most efficiently multi-threaded over the wavelength dimension, which
works very well when the number of wavelengths is high (a factor of at least 2 larger than
the number of threads).  When the number of wavelengths is small, it may be better
to multi-thread over other parts of the code instead.  This behaviour can be changed by
setting

```{code}
config.threading_model = sk.ThreadingModel.Source
```
which will attempt to multi-thread the calculation over the source function calculation.
This threading is usually significantly less efficient than the wavelength dimension
multi-threading, but has the advantage that it works well for a small number of wavelengths, and also
may use significantly less RAM.

```{note}
The `sk.ThreadingModel.Source` works well for the `sk.MultipleScatterSource.SuccessiveOrders` source,
but may not offer much improvement for the `sk.MultipleScatterSource.DiscreteOrdinates` source.
```

## Installation Method
For maximum performance, we recommend installing SASKTRAN2 through the `conda` packages rather than
the `pip` wheels.  The conda packages are able to use accelerated LAPACK/BLAS libraries for your system
rather than the bundled OpenBLAS library in the wheels.  You can check what BLAS vendor you have installed
by running `conda list | grep libblas` from inside your conda environment.  You should see something
like `20_osxarm64_accelerate`.  Which BLAS library offers the best performance is very dependent on what
platform you are on as well as the specific calculation you are performing.  If you want the absolute best
performance we recommend trying a few different one.
